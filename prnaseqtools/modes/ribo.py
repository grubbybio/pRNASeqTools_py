"""
Ribo-seq analysis mode — RIBO Taper pipeline.

Workflow:
  1. Bowtie2 decontamination (rRNA, tRNA, snRNA, snoRNA)
  2. RNA-seq STAR 2-pass mapping + StringTie transcriptome assembly
  3. gffcompare → R filter novel transcripts → updated GTF with gene_biotype
  4. RSEM quantification → R filter expressed isoforms (TPM > 0)
  5. STAR re-mapping (RNA-seq + Ribo-seq) with expressed annotation
  6. RIBO Taper annotation + ORF detection
  7. RSEM CDS quantification
  8. P-site analysis & visualization (metagene plots)
  9. Translation Efficiency (TE) = Ribo TPM / RNA TPM per gene

Requires: Bowtie2, STAR, samtools, bedtools, StringTie, gffcompare, RSEM,
          RIBO Taper scripts, R (dplyr)
"""

import os
import sys
import re
import math
import glob as globmod
import shutil
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input, _parse_to_dict
from prnaseqtools.functions import download_sra, unzip_file, _tee, run_cmd


# ── Helper: detect last completed step from log ────────────────────────────

def _detect_last_step(search_dir=None, exclude_log=None):
    """Detect the last completed step from log files.

    Args:
        search_dir: Directory to scan (defaults to os.getcwd()).
        exclude_log: Log filename to skip (the current run's log).
                     Pass log_{start_time}.txt to avoid reading the
                     current (almost empty) log as the resume source.

    Returns:
        tuple: (last_step, log_path) where last_step is int (0 if none found),
               log_path is the file path or None
    """
    search_dir = search_dir or os.getcwd()
    log_files = globmod.glob(os.path.join(search_dir, "log_*.txt"))
    if not log_files:
        return (0, None)

    # Exclude the current run's log so we find the PREVIOUS run
    candidates = [f for f in log_files if os.path.basename(f) != exclude_log] if exclude_log else log_files
    if not candidates:
        return (0, None)

    # Find the most recently modified log file (among previous runs)
    candidates.sort(key=os.path.getmtime, reverse=True)
    latest_log = candidates[0]

    last_step = 0
    with open(latest_log) as f:
        for line in f:
            match = re.search(r'STEP\s+(\d+)\s+COMPLETE', line)
            if match:
                step = int(match.group(1))
                if step > last_step:
                    last_step = step

    return (last_step, latest_log)


# R packages required by the TE statistical testing module (Step 11)
_TE_R_PACKAGES = ['emmeans', 'car', 'multcomp']


def _check_r_packages(tee):
    """Check and auto-install the R packages required for TE statistics.

    Called once at pipeline start. Uses scripts/checkPackages.R for
    installation when packages are missing.
    """
    rscript = shutil.which("Rscript")
    if not rscript:
        tee.write("  WARNING: Rscript not found - TE statistics may fail\n")
        return

    check_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        'scripts', 'checkPackages.R')

    # Build an inline R expression that reports missing packages
    pkg_str = ','.join(_TE_R_PACKAGES)
    probe = (
        f'missing <- c(); '
        f'for (p in strsplit("{pkg_str}", ",")[[1]]) '
        f'if (!requireNamespace(p, quietly=TRUE)) missing <- c(missing, p); '
        f'cat(if (length(missing)) paste(missing, collapse=",") else "")'
    )
    try:
        r = subprocess.run([rscript, '--vanilla', '-e', probe],
                           capture_output=True, text=True, timeout=120)
        missing = r.stdout.strip()
        if not missing:
            tee.write(f"  R packages OK: {pkg_str}\n")
            return
        tee.write(f"  Missing R packages: {missing} - attempting installation...\n")
        if os.path.exists(check_script):
            run_cmd(
                f'Rscript --vanilla "{check_script}" '
                f'--packages "{missing}"')
        else:
            run_cmd(
                f"Rscript --vanilla -e "
                f"'if(!requireNamespace(\"BiocManager\",quietly=TRUE))"
                f"install.packages(\"BiocManager\");"
                f"install.packages(c({','.join(repr(p) for p in _TE_R_PACKAGES)}))'")
    except Exception as exc:
        tee.write(f"  WARNING: R package check failed: {exc}\n")


# ── Main entry point ─────────────────────────────────────────────────────

def run(opts):
    """Main entry point for Ribo-seq analysis (RIBO Taper pipeline)."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))

    # ── RIBO Taper-specific options ──────────────────────────────────
    contam = opts.get('contam')                # contamination fasta
    if not contam:
        contam = os.path.join(prefix, "reference", f"{genome}_contam4.fa")
        if not os.path.exists(contam):
            sys.exit(
                f"Contamination file not found: {contam}\n"
                f"  Please provide --contam or place {genome}_contam4.fa "
                f"in the reference directory."
            )
    ribo_len = opts.get('ribo_len', '24,25,26,27,28')
    cutoffs = opts.get('cutoffs', '8,9,10,11,12')
    tpm_threshold = opts.get('tpm_threshold', 0)
    ribotaper_path = opts.get('ribotaper')
    ribotaper_env = opts.get('ribotaper_env')

    def _ribo_bash(cmd):
        """Wrap a bash command with conda run if ribotaper_env is set."""
        if ribotaper_env:
            return f"conda run -n {ribotaper_env} bash {cmd}"
        return f"bash {cmd}"

    # ── Parse RNA-seq inputs ─────────────────────────────────────────
    # rna_groups records [(group_name, n_replicates), ...] in input order
    # (control first, then treatments) so TE calculation can pair samples
    # by experimental role instead of raw list position.
    rna_groups = []
    rna_ctrl_dict = _parse_to_dict(opts.get('rna_control', ''))
    rna_tags, rna_files, rna_pars = parse_input(rna_ctrl_dict)
    if len(rna_pars) >= 2:
        rna_groups.append((str(rna_pars[0]), int(rna_pars[1])))

    rna_trt = opts.get('rna_treatment')
    if rna_trt:
        if isinstance(rna_trt, str):
            rna_trt = [rna_trt]
        for trt_str in rna_trt:
            trt_dict = _parse_to_dict(trt_str)
            rt_tags, rt_files, rt_pars = parse_input(trt_dict)
            if len(rt_pars) >= 2:
                rna_groups.append((str(rt_pars[0]), int(rt_pars[1])))
            rna_tags.extend(rt_tags)
            rna_files.extend(rt_files)

    all_rna_tags = list(rna_tags)

    # ── Parse Ribo-seq inputs ────────────────────────────────────────
    ribo_groups = []
    ribo_ctrl_dict = _parse_to_dict(opts.get('ribo_control', ''))
    ribo_tags, ribo_files, ribo_pars = parse_input(ribo_ctrl_dict)
    if len(ribo_pars) >= 2:
        ribo_groups.append((str(ribo_pars[0]), int(ribo_pars[1])))

    ribo_trt = opts.get('ribo_treatment')
    if ribo_trt:
        if isinstance(ribo_trt, str):
            ribo_trt = [ribo_trt]
        for rbt_str in ribo_trt:
            rbt_dict = _parse_to_dict(rbt_str)
            rbt_tags, rbt_files, rbt_pars = parse_input(rbt_dict)
            if len(rbt_pars) >= 2:
                ribo_groups.append((str(rbt_pars[0]), int(rbt_pars[1])))
            ribo_tags.extend(rbt_tags)
            ribo_files.extend(rbt_files)

    all_ribo_tags = list(ribo_tags)

    # ── Reference paths ──────────────────────────────────────────────
    ref_dir = os.path.join(prefix, "reference")
    gtf_path = os.path.join(ref_dir, f"{genome}_genes.gtf")
    gff_path = os.path.join(ref_dir, f"{genome}_genes.gff")
    fasta_path = os.path.join(ref_dir, f"{genome}_chr_all.fasta")

    # Prefer GTF; if missing, convert from GFF
    if os.path.exists(gtf_path):
        pass  # GTF already exists
    elif os.path.exists(gff_path):
        tee.write(f"  GTF not found, converting GFF → GTF with gffread...\n")
        run_cmd(f"gffread -T -o {gtf_path} {gff_path}")
    else:
        sys.exit(
            f"Reference annotation not found: neither {gtf_path} nor {gff_path} exists.\n"
            f"  Please provide {genome}_genes.gtf or {genome}_genes.gff in {ref_dir}")

    if not os.path.exists(fasta_path):
        sys.exit(f"Reference FASTA not found: {fasta_path}")

    # Use GTF as the annotation reference throughout
    anno_path = gtf_path

    # ── Utility: find bedtools ───────────────────────────────────────
    bedtools_bin = shutil.which("bedtools") or "bedtools"

    # ── Shared directory variables ───────────────────────────────────
    expressed_gtf = f"{genome}_expressed.gtf"
    updated_gtf = f"{genome}_updated.gtf"
    rsem_dir = "RSEM_results"
    ribo_map_dir = "STAR_Ribo_map"
    star_rna_new_map = "STAR_RNA_map_new"
    contam_index = "Contam"
    ribo_anno_dir = "RiboTaper_annotation"
    # Defined here (not inside step 3) so step 4 can use them when resuming
    # with last_step >= 3 or --restart-step >= 4.
    assembled_gtf_dir = "assembledGTF"
    gffcomp_prefix = os.path.join(assembled_gtf_dir, genome)

    # ── Auto-detect last completed step from log ─────────────────────
    _start_time = opts.get('_start_time', 0)
    current_log = f"log_{_start_time}.txt"
    last_step, resume_log = _detect_last_step(exclude_log=current_log)

    # ── Manual restart-step override ──────────────────────────────
    _restart = opts.get('restart_step')
    if _restart is not None:
        if not 1 <= _restart <= 12:
            sys.exit(f"--restart-step must be 1-12, got {_restart}")
        last_step = _restart - 1  # step N → set last_step to N-1 so step N re-runs
        tee.write(f"\n  Manual restart from step {_restart} "
                  f"(--restart-step={_restart})\n\n")
    elif last_step > 0:
        tee.write(f"\n  Resuming from step {last_step + 1} "
                  f"(last completed: step {last_step})\n")
        tee.write(f"  Log file: {resume_log}\n\n")
    else:
        tee.write("  No previous log found, starting from step 1.\n\n")

        # ── Verify intermediate files if resuming from step 6+ ──
    if last_step >= 6:
        if not os.path.exists(expressed_gtf):
            tee.write(f"  WARNING: {expressed_gtf} not found - steps 4-5 may need to be re-run\n")
        if not os.path.exists(updated_gtf):
            tee.write(f"  WARNING: {updated_gtf} not found - steps 4-5 may need to be re-run\n")

    # ── Check & install required R packages (for TE statistics) ──────
    _check_r_packages(tee)

    # ==================================================================
    # STEP 1 — Build Bowtie2 contamination index
    # ==================================================================
    if last_step < 1:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 1: Building Bowtie2 contamination index\n")
        tee.write("=" * 60 + "\n")

        # Cleanup previous incomplete step 1
        for pat in [f"{contam_index}*.ebwt", f"{contam_index}*.bt2"]:
            for fname in globmod.glob(pat):
                os.unlink(fname)

        run_cmd(
            f"bowtie2-build {contam} {contam_index}", quiet=True)
        tee.write("  Contamination index built.\n")
        tee.write("  STEP 1 COMPLETE\n")

    # ==================================================================
    # STEP 2 — Preprocess Ribo-seq reads (decontamination)
    # ==================================================================
    if last_step < 2:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 2: Preprocessing Ribo-seq reads\n")
        tee.write("=" * 60 + "\n")

        # Cleanup previous incomplete step 2 outputs
        for tag in all_ribo_tags:
            for pat in [f"{tag}.noContam*", f"{tag}.mapped_and_unmapped*",
                         f"{tag}.bothEndsUnmapped*", f"{tag}.unmapped*"]:
                for fname in globmod.glob(pat):
                    if os.path.exists(fname):
                        os.unlink(fname)

        for i in range(len(ribo_tags)):
            tag = ribo_tags[i]
            fpath = ribo_files[i]

            tee.write(f"\n  --- Ribo-seq sample: {tag} ---\n")

            if ',' not in fpath:
                # Download / unzip (SRA or single file)
                sra_results = download_sra(fpath, thread)
                is_paired = len(sra_results) > 1

                if is_paired:
                    unzip_file(sra_results[0], f"{tag}_r1")
                    unzip_file(sra_results[1], f"{tag}_r2")
                    r1_fq = f"{tag}_r1.fastq"
                    r2_fq = f"{tag}_r2.fastq"
                else:
                    unzip_file(sra_results[0], tag)
                    r1_fq = f"{tag}.fastq"
            else:
                # Explicit paired-end (file1,file2)
                is_paired = True
                f1, f2 = fpath.split(',')
                unzip_file(f1, f"{tag}_r1")
                unzip_file(f2, f"{tag}_r2")
                r1_fq = f"{tag}_r1.fastq"
                r2_fq = f"{tag}_r2.fastq"

            # Adaptor trimming
            if adaptor:
                if is_paired:
                    tee.write(f"  Trimming {tag} (paired-end)...\n")
                    run_cmd(
                        f"cutadapt -j {thread} -m 18 --trim-n "
                        f"-a {adaptor} -A {adaptor} "
                        f"-o {tag}_trimmed_r1.fastq -p {tag}_trimmed_r2.fastq "
                        f"{r1_fq} {r2_fq}")
                    os.rename(f"{tag}_trimmed_r1.fastq", r1_fq)
                    os.rename(f"{tag}_trimmed_r2.fastq", r2_fq)
                else:
                    tee.write(f"  Trimming {tag}...\n")
                    run_cmd(
                        f"cutadapt -j {thread} -m 18 --discard-untrimmed --trim-n "
                        f"-a {adaptor} -o {tag}_trimmed.fastq {r1_fq}")
                    os.rename(f"{tag}_trimmed.fastq", r1_fq)

            # Bowtie2 decontamination
            tee.write(f"  Removing contamination...\n")
            if is_paired:
                run_cmd(
                    f"bowtie2 -L 20 -p {thread} -x {contam_index} "
                    f"-1 {r1_fq} -2 {r2_fq} "
                    f"-S {tag}.mapped_and_unmapped.sam")
            else:
                run_cmd(
                    f"bowtie2 -L 20 -p {thread} -x {contam_index} "
                    f"-U {r1_fq} "
                    f"-S {tag}.mapped_and_unmapped.sam")

            # Convert SAM to BAM
            run_cmd(
                f"samtools view -bS -o {tag}.mapped_and_unmapped.bam "
                f"{tag}.mapped_and_unmapped.sam")

            # Extract unmapped reads
            if is_paired:
                # -f 12: both ends unmapped; -F 256: discard secondary alignments
                unmapped_flag = "-f 12 -F 256"
                run_cmd(
                    f"samtools view -b {unmapped_flag} "
                    f"-o {tag}.bothEndsUnmapped.bam "
                    f"{tag}.mapped_and_unmapped.bam")
                # Split into R1/R2 fastq
                run_cmd(
                    f"samtools sort -n -o {tag}.bothEndsUnmapped_sorted.bam "
                    f"{tag}.bothEndsUnmapped.bam")
                run_cmd(
                    f"bedtools bamtofastq -i {tag}.bothEndsUnmapped_sorted.bam "
                    f"-fq {tag}.r1.fastq -fq2 {tag}.r2.fastq")
                # Reuse r1_fq/r2_fq pointing to clean files
                os.rename(f"{tag}.r1.fastq", r1_fq)
                os.rename(f"{tag}.r2.fastq", r2_fq)
            else:
                # -f 4: unmapped; -F 256: discard secondary
                run_cmd(
                    f"samtools view -b -f 4 -F 256 "
                    f"-o {tag}.unmapped.bam "
                    f"{tag}.mapped_and_unmapped.bam")
                run_cmd(
                    f"samtools sort -n -o {tag}.unmapped_sorted.bam "
                    f"{tag}.unmapped.bam")
                run_cmd(
                    f"bedtools bamtofastq -i {tag}.unmapped_sorted.bam "
                    f"-fq {tag}.clean.fastq")
                os.rename(f"{tag}.clean.fastq", r1_fq)

            # Gzip clean fastq
            if is_paired:
                run_cmd(f"gzip -f {r1_fq}")
                run_cmd(f"gzip -f {r2_fq}")
                # Store as noContam paired files for later use
                os.rename(f"{r1_fq}.gz", f"{tag}.noContam_r1.fastq.gz")
                os.rename(f"{r2_fq}.gz", f"{tag}.noContam_r2.fastq.gz")
            else:
                run_cmd(f"gzip -f {r1_fq}")
                os.rename(f"{r1_fq}.gz", f"{tag}.noContam.fastq.gz")

            # Cleanup intermediate files
            for pat in [f"{tag}.mapped_and_unmapped.*",
                         f"{tag}.bothEndsUnmapped*",
                         f"{tag}.unmapped*",
                         f"{tag}.clean.fastq",
                         f"{tag}.r1.fastq",
                         f"{tag}.r2.fastq"]:
                for fname in globmod.glob(pat):
                    if os.path.exists(fname):
                        os.unlink(fname)

            if not is_paired and os.path.exists(f"{tag}.fastq"):
                os.unlink(f"{tag}.fastq")
            if is_paired and os.path.exists(f"{tag}_r1.fastq"):
                os.unlink(f"{tag}_r1.fastq")
            if is_paired and os.path.exists(f"{tag}_r2.fastq"):
                os.unlink(f"{tag}_r2.fastq")

        tee.write("\n  Ribo-seq decontamination complete.\n")
        tee.write("  STEP 2 COMPLETE\n")

    # ==================================================================
    # STEP 3 — RNA-seq STAR mapping + StringTie transcriptome assembly
    # ==================================================================
    if last_step < 3:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 3: RNA-seq transcriptome assembly (STAR + StringTie)\n")
        tee.write("=" * 60 + "\n")

        # Cleanup previous incomplete step 3
        for d in ["STAR_RNA_index", "STAR_SJ", "STAR_RNA_2nd", "assembledGTF"]:
            if os.path.exists(d):
                import shutil as _shutil
                _shutil.rmtree(d, ignore_errors=True)
        for tag in rna_tags:
            for d in [f"star_{tag}_1st"]:
                if os.path.exists(d):
                    import shutil as _shutil
                    _shutil.rmtree(d, ignore_errors=True)

        # Build STAR genome index for RNA-seq
        star_rna_index = "STAR_RNA_index"
        if os.path.exists(star_rna_index):
            run_cmd(f"rm -rf {star_rna_index}")
        os.makedirs(star_rna_index, exist_ok=True)

        run_cmd(
            f"STAR --runThreadN {thread} --runMode genomeGenerate "
            f"--genomeDir {star_rna_index} --genomeFastaFiles {fasta_path} "
            f"--sjdbGTFfile {anno_path} --sjdbOverhang 99 "
            f"--limitGenomeGenerateRAM 64000000000")

        # STAR 1st pass for each RNA-seq sample
        star_sj_dir = "STAR_SJ"
        os.makedirs(star_sj_dir, exist_ok=True)

        for tag, fpath in zip(rna_tags, rna_files):
            tee.write(f"\n  --- RNA-seq STAR 1st pass: {tag} ---\n")

            if ',' not in fpath:
                sra_r = download_sra(fpath, thread)
                rna_paired = len(sra_r) > 1

                if rna_paired:
                    unzip_file(sra_r[0], f"{tag}_r1")
                    unzip_file(sra_r[1], f"{tag}_r2")
                    run_cmd(f"gzip -f {tag}_r1.fastq")
                    run_cmd(f"gzip -f {tag}_r2.fastq")
                    read_files = f"{tag}_r1.fastq.gz {tag}_r2.fastq.gz"
                else:
                    unzip_file(sra_r[0], tag)
                    run_cmd(f"gzip -f {tag}.fastq")
                    read_files = f"{tag}.fastq.gz"
            else:
                # Explicit paired-end (file1,file2)
                rna_paired = True
                f1, f2 = fpath.split(',')
                unzip_file(f1, f"{tag}_r1")
                unzip_file(f2, f"{tag}_r2")
                run_cmd(f"gzip -f {tag}_r1.fastq")
                run_cmd(f"gzip -f {tag}_r2.fastq")
                read_files = f"{tag}_r1.fastq.gz {tag}_r2.fastq.gz"

            star_out_1st = f"star_{tag}_1st"
            os.makedirs(star_out_1st, exist_ok=True)

            run_cmd(
                f"STAR --runThreadN {thread} --genomeDir {star_rna_index} "
                f"--readFilesCommand zcat --readFilesIn {read_files} "
                f"--alignIntronMax 5000 --alignIntronMin 15 "
                f"--outFilterMismatchNmax 2 --outFilterMultimapNmax 20 "
                f"--outFilterType BySJout --alignSJoverhangMin 8 "
                f"--alignSJDBoverhangMin 2 "
                f"--outSAMtype BAM SortedByCoordinate "
                f"--quantMode TranscriptomeSAM --outSAMmultNmax 1 "
                f"--outMultimapperOrder Random "
                f"--outFileNamePrefix {star_out_1st}/star_{tag}_ ")

            # Collect splice junctions
            sj_src = f"{star_out_1st}/star_{tag}_SJ.out.tab"
            if os.path.exists(sj_src):
                os.rename(sj_src, os.path.join(star_sj_dir, f"{tag}_SJ.out.tab"))

            # Remove transcriptome BAM (not needed for assembly)
            tx_bam = f"{star_out_1st}/star_{tag}_Aligned.toTranscriptome.out.bam"
            if os.path.exists(tx_bam):
                os.unlink(tx_bam)

        # STAR 2nd pass with all splice junctions
        sj_files = ' '.join(sorted([
            os.path.join(star_sj_dir, f)
            for f in os.listdir(star_sj_dir) if f.endswith('.tab')
        ]))

        star_out_2nd = "STAR_RNA_2nd"
        os.makedirs(star_out_2nd, exist_ok=True)

        for tag in rna_tags:
            tee.write(f"\n  --- RNA-seq STAR 2nd pass: {tag} ---\n")

            r1_gz = f"{tag}_r1.fastq.gz"
            r2_gz = f"{tag}_r2.fastq.gz"
            if os.path.exists(r2_gz):
                read_files = f"{r1_gz} {r2_gz}"
            else:
                read_files = f"{tag}.fastq.gz"

            run_cmd(
                f"STAR --runThreadN {thread} --genomeDir {star_rna_index} "
                f"--sjdbFileChrStartEnd {sj_files} "
                f"--readFilesCommand zcat --readFilesIn {read_files} "
                f"--alignIntronMax 5000 --alignIntronMin 15 "
                f"--outFilterMismatchNmax 2 --outFilterMultimapNmax 20 "
                f"--outFilterType BySJout --alignSJoverhangMin 8 "
                f"--alignSJDBoverhangMin 2 "
                f"--outSAMtype BAM SortedByCoordinate "
                f"--quantMode TranscriptomeSAM --outSAMmultNmax 1 "
                f"--outMultimapperOrder Random "
                f"--outFileNamePrefix {star_out_2nd}/star_{tag}_ ")

            bam_file = f"{star_out_2nd}/star_{tag}_Aligned.sortedByCoord.out.bam"
            run_cmd(f"samtools index {bam_file}")

            tx_bam = f"{star_out_2nd}/star_{tag}_Aligned.toTranscriptome.out.bam"
            if os.path.exists(tx_bam):
                os.unlink(tx_bam)

        # ── StringTie per sample ─────────────────────────────────────────
        tee.write("\n  --- StringTie assembly ---\n")
        os.makedirs(assembled_gtf_dir, exist_ok=True)

        merge_list_paths = []
        for tag in rna_tags:
            tee.write(f"  StringTie: {tag}\n")
            bam_file = f"{star_out_2nd}/star_{tag}_Aligned.sortedByCoord.out.bam"
            gtf_out = os.path.join(assembled_gtf_dir, f"{tag}.gtf")

            run_cmd(
                f"stringtie --rf -p {thread} -G {anno_path} "
                f"-o {gtf_out} -l {tag} {bam_file}")
            merge_list_paths.append(gtf_out)

        # Merge transcripts
        tee.write("\n  Merging transcripts...\n")
        merge_list_file = os.path.join(assembled_gtf_dir, "mergeList.txt")
        with open(merge_list_file, 'w') as fh:
            for p in merge_list_paths:
                fh.write(p + '\n')

        merged_gtf = os.path.join(assembled_gtf_dir, f"{genome}_merged.gtf")
        run_cmd(
            f"stringtie --merge -p {thread} -T 0.05 -G {anno_path} "
            f"-o {merged_gtf} {merge_list_file}")

        # gffcompare
        tee.write("  Running gffcompare...\n")
        run_cmd(
            f"gffcompare -V -r {anno_path} -o {gffcomp_prefix} {merged_gtf}")
        tee.write("  STEP 3 COMPLETE\n")

    # ==================================================================
    # STEP 4 — R: Filter novel transcripts + add gene_biotype
    # ==================================================================
    if last_step < 4:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 4: Filtering novel transcripts + gene_biotype annotation\n")
        tee.write("=" * 60 + "\n")

        # Cleanup previous incomplete step 4
        if os.path.exists(updated_gtf):
            os.unlink(updated_gtf)
        # Also clean gffcompare outputs that might be stale
        for pat in [os.path.join("assembledGTF", f"{genome}.*")]:
            for fname in globmod.glob(pat):
                if os.path.exists(fname):
                    os.unlink(fname)

        annotated_gtf = f"{gffcomp_prefix}.annotated.gtf"
        if not os.path.exists(annotated_gtf):
            # gffcompare might use a different naming convention
            candidates = globmod.glob(os.path.join(assembled_gtf_dir, "*.annotated.gtf"))
            if candidates:
                annotated_gtf = candidates[0]
            else:
                sys.exit("gffcompare annotated GTF not found")

        updated_gtf = f"{genome}_updated.gtf"

        run_cmd(
            f"Rscript --vanilla {prefix}/scripts/ribotaper_filter_gtf.R "
            f"{annotated_gtf} {anno_path} {updated_gtf}")
        tee.write(f"  Updated GTF written to: {updated_gtf}\n")
        tee.write("  STEP 4 COMPLETE\n")

    # ==================================================================
    # STEP 5 — RSEM quantification + filter expressed isoforms
    # ==================================================================
    if last_step < 5:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 5: RSEM quantification + expressed isoform filtering\n")
        tee.write("=" * 60 + "\n")

        # Cleanup previous incomplete step 5
        for d in ["RSEM_index", "RSEM_results"]:
            if os.path.exists(d):
                import shutil as _shutil
                _shutil.rmtree(d, ignore_errors=True)
        if os.path.exists(expressed_gtf):
            os.unlink(expressed_gtf)

        # Build RSEM index
        rsem_index = "RSEM_index"
        os.makedirs(rsem_index, exist_ok=True)
        star_bin = shutil.which("STAR") or "STAR"

        run_cmd(
            f"rsem-prepare-reference --gtf {updated_gtf} "
            f"--star --star-path {os.path.dirname(star_bin)} "
            f"--star-sjdboverhang 99 -p {thread} "
            f"{fasta_path} {rsem_index}/RNA", quiet=True)

        # Run RSEM for each RNA-seq sample
        rsem_dir = "RSEM_results"
        os.makedirs(rsem_dir, exist_ok=True)

        for tag in rna_tags:
            tee.write(f"\n  RSEM: {tag}\n")

            r1_gz = f"{tag}_r1.fastq.gz"
            r2_gz = f"{tag}_r2.fastq.gz"
            if os.path.exists(r2_gz):
                paired_flag = "--paired-end"
                read_files = f"{r1_gz} {r2_gz}"
            else:
                paired_flag = ""
                read_files = f"{tag}.fastq.gz"

            # Use RSEM's built-in STAR alignment
            run_cmd(
                f"rsem-calculate-expression {paired_flag} "
                f"--star --star-path {os.path.dirname(star_bin)} "
                f"--star-gzipped-read-file "
                f"-p {thread} --time --strandedness reverse "
                f"--no-bam-output "
                f"{read_files} {rsem_index}/RNA {rsem_dir}/{tag}", quiet=True)

        # Filter expressed isoforms with R
        tee.write(f"\n  Filtering expressed isoforms (mean TPM > {tpm_threshold})...\n")
        expressed_gtf = f"{genome}_expressed.gtf"

        rsem_files_arg = ' '.join([
            f"{rsem_dir}/{tag}.isoforms.results" for tag in rna_tags
        ])
        run_cmd(
            f"Rscript --vanilla {prefix}/scripts/ribotaper_filter_rsem.R "
            f"{updated_gtf} {expressed_gtf} {tpm_threshold} {rsem_files_arg}")
        tee.write(f"  Expressed GTF written to: {expressed_gtf}\n")
        tee.write("  STEP 5 COMPLETE\n")

    # ==================================================================
    # STEP 6 — STAR re-mapping with expressed annotation
    # ==================================================================
    if last_step < 6:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 6: STAR re-mapping with expressed annotation\n")
        tee.write("=" * 60 + "\n")

        # Cleanup previous incomplete step 6
        for d in ["STAR_Ribo_index", "STAR_Ribo_map",
                   "STAR_RNA_index_new", "STAR_RNA_map_new"]:
            if os.path.exists(d):
                import shutil as _shutil
                _shutil.rmtree(d, ignore_errors=True)

        # 6a — Build Ribo-seq STAR index (sjdbOverhang based on read length)
        star_ribo_idx = "STAR_Ribo_index"
        os.makedirs(star_ribo_idx, exist_ok=True)
        ribo_overhang = min(int(ribo_len.split(',')[0]) - 1, 99)

        run_cmd(
            f"STAR --runThreadN {thread} --runMode genomeGenerate "
            f"--genomeDir {star_ribo_idx} --genomeFastaFiles {fasta_path} "
            f"--sjdbGTFfile {expressed_gtf} --sjdbOverhang {ribo_overhang} "
            f"--limitGenomeGenerateRAM 64000000000")

        # Map Ribo-seq reads
        os.makedirs(ribo_map_dir, exist_ok=True)

        for tag in all_ribo_tags:
            tee.write(f"\n  STAR Ribo-seq: {tag}\n")
            fq_gz = f"{tag}.noContam.fastq.gz"
            if not os.path.exists(fq_gz):
                # Try paired naming
                r1 = f"{tag}.noContam_r1.fastq.gz"
                if os.path.exists(r1):
                    fq_gz = r1
                else:
                    tee.write(f"  Warning: Ribo-seq fastq not found for {tag}, skipping\n")
                    continue

            run_cmd(
                f"STAR --runThreadN {thread} --genomeDir {star_ribo_idx} "
                f"--alignEndsType EndToEnd --readFilesCommand zcat "
                f"--readFilesIn {fq_gz} "
                f"--alignIntronMax 5000 --alignIntronMin 15 "
                f"--outFilterMismatchNmax 2 --outFilterMultimapNmax 20 "
                f"--outFilterType BySJout --alignSJoverhangMin 4 "
                f"--alignSJDBoverhangMin 1 "
                f"--outSAMtype BAM SortedByCoordinate "
                f"--quantMode TranscriptomeSAM --outSAMmultNmax 1 "
                f"--outMultimapperOrder Random "
                f"--outFileNamePrefix {ribo_map_dir}/star_{tag}_ ")

    # ── STEP 6b: RNA-seq STAR mapping ──
    if last_step < 6:
        # 6b — Build updated RNA-seq STAR index
        star_rna_new_idx = "STAR_RNA_index_new"
        os.makedirs(star_rna_new_idx, exist_ok=True)

        run_cmd(
            f"STAR --runThreadN {thread} --runMode genomeGenerate "
            f"--genomeDir {star_rna_new_idx} --genomeFastaFiles {fasta_path} "
            f"--sjdbGTFfile {expressed_gtf} --sjdbOverhang 99 "
            f"--limitGenomeGenerateRAM 64000000000")

        # Map RNA-seq with new annotation
        star_rna_new_map = "STAR_RNA_map_new"
        os.makedirs(star_rna_new_map, exist_ok=True)

        for tag in all_rna_tags:
            tee.write(f"\n  STAR RNA-seq (updated): {tag}\n")
            r1_gz = f"{tag}_r1.fastq.gz"
            r2_gz = f"{tag}_r2.fastq.gz"
            if os.path.exists(r2_gz):
                read_files = f"{r1_gz} {r2_gz}"
            else:
                read_files = f"{tag}.fastq.gz"

            run_cmd(
                f"STAR --runThreadN {thread} --genomeDir {star_rna_new_idx} "
                f"--readFilesCommand zcat --readFilesIn {read_files} "
                f"--alignIntronMax 5000 --alignIntronMin 15 "
                f"--outFilterMismatchNmax 2 --outFilterMultimapNmax 20 "
                f"--outFilterType BySJout --alignSJoverhangMin 2 "
                f"--alignSJDBoverhangMin 1 "
                f"--outSAMtype BAM SortedByCoordinate "
                f"--quantMode TranscriptomeSAM --outSAMmultNmax 1 "
                f"--outMultimapperOrder Random "
                f"--outFileNamePrefix {star_rna_new_map}/star_{tag}_ ")
        tee.write("  STEP 6 COMPLETE\n")

    # ==================================================================
    # STEP 7 — RIBO Taper annotation
    # ==================================================================
    if last_step < 7:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 7: Creating RIBO Taper annotation files\n")
        tee.write("=" * 60 + "\n")

        # Cleanup previous incomplete step 7
        if os.path.exists(ribo_anno_dir):
            import shutil as _shutil
            _shutil.rmtree(ribo_anno_dir, ignore_errors=True)

        os.makedirs(ribo_anno_dir, exist_ok=True)

        create_anno_script = os.path.join(
            ribotaper_path, "create_annotations_files.bash"
        )
        if not os.path.exists(create_anno_script):
            sys.exit(f"RIBO Taper script not found: {create_anno_script}")

        from prnaseqtools.reference import check_and_build_indices
        check_and_build_indices(prefix, genome=genome, mode='ribo', tee=tee)

        # RIBO Taper create_annotations_files.bash usage:
        #   ./create_annotation_files.bash <gtf_file> <genome_fasta_file> \
        #       <use_ccdsid?> <use_appris?> <dest_folder>
        # use_ccdsid: use CCDS IDs when mapping (true/false)
        # use_appris: use APPRIS principal transcripts (true/false)
        run_cmd(
            _ribo_bash(f"{create_anno_script} {expressed_gtf} {fasta_path} "
            f"false false {ribo_anno_dir}"))
        tee.write("  RIBO Taper annotation created.\n")
        tee.write("  STEP 7 COMPLETE\n")

    # ==================================================================
    # STEP 8 — Merge BAMs + Run RIBO Taper
    # ==================================================================
    if last_step < 8:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 8: Merging BAMs + Running RIBO Taper\n")
        tee.write("=" * 60 + "\n")

        # Cleanup previous incomplete step 8
        for f in ["RNA_merged.bam", "RNA_merged.bam.bai",
                   "Ribo_merged.bam", "Ribo_merged.bam.bai"]:
            if os.path.exists(f):
                os.unlink(f)
        # Clean per-sample R1-only BAMs
        for tag in all_ribo_tags:
            r1 = f"{ribo_map_dir}/star_{tag}_Aligned.sortedByCoord.out_r1.bam"
            for f in [r1, r1 + '.bai']:
                if os.path.exists(f):
                    os.unlink(f)
        # Clean RIBO Taper output files (P_sites, ORFs, etc.)
        for pat in ["P_sites*", "ORF*", "start_stops*",
                     "all_tracks*", "reads*", "offsets*",
                     "Ribo_metaplots*"]:
            for fname in globmod.glob(pat):
                if os.path.exists(fname):
                    try:
                        if os.path.isdir(fname):
                            import shutil as _shutil
                            _shutil.rmtree(fname, ignore_errors=True)
                        else:
                            os.unlink(fname)
                    except OSError:
                        pass

        # Merge RNA-seq BAMs
        tee.write("\n  Merging RNA-seq BAMs...\n")
        rna_bam_list = [
            f"{star_rna_new_map}/star_{tag}_Aligned.sortedByCoord.out.bam"
            for tag in all_rna_tags
        ]
        rna_bam_list = [b for b in rna_bam_list if os.path.exists(b)]

        rna_merged = "RNA_merged.bam"
        _merge_or_copy_bams(rna_bam_list, rna_merged, tee=tee)

        # Merge Ribo-seq BAMs (R1 only for RIBO Taper)
        tee.write("  Merging Ribo-seq BAMs...\n")
        ribo_bam_list = [
            f"{ribo_map_dir}/star_{tag}_Aligned.sortedByCoord.out.bam"
            for tag in all_ribo_tags
        ]
        ribo_bam_list = [b for b in ribo_bam_list if os.path.exists(b)]

        ribo_r1_bams = []
        for bam in ribo_bam_list:
            r1_bam = bam.replace('.bam', '_r1.bam')
            run_cmd(
                f"samtools view -h -b -F 128 {bam} | "
                f"samtools sort -o {r1_bam}")
            run_cmd(f"samtools index {r1_bam}")
            ribo_r1_bams.append(r1_bam)

        ribo_merged = "Ribo_merged.bam"
        _merge_or_copy_bams(ribo_r1_bams, ribo_merged, tee=tee)

        # ── Generate metaplots for parameter selection ────────────────
        tee.write("\n  Generating metaplots for read-length/cutoff selection...\n")
        metaplot_script = os.path.join(
            ribotaper_path, "create_metaplots.bash"
        )
        start_stops_bed = os.path.join(ribo_anno_dir, "start_stops_FAR.bed")
        if not os.path.exists(metaplot_script):
            tee.write(f"  Warning: create_metaplots.bash not found at "
                      f"{metaplot_script}, skipping metaplots\n")
        elif not os.path.exists(start_stops_bed):
            tee.write(f"  Warning: start_stops_FAR.bed not found at "
                      f"{start_stops_bed}, skipping metaplots\n")
        else:
            metaplot_name = "Ribo_metaplots"
            run_cmd(
                _ribo_bash(
                    f"{metaplot_script} {ribo_merged} "
                    f"{start_stops_bed} {metaplot_name}"))
            tee.write(f"  Metaplots generated: {metaplot_name}*.pdf\n")

        # ── Interactive parameter confirmation ────────────────────────
        tee.write("\n  ── RIBO Taper parameter setup ──\n")
        tee.write(f"  Current ribo-len: {ribo_len}\n")
        tee.write(f"  Current cutoffs:  {cutoffs}\n")
        try:
            _new_len = input(
                f"  Press Enter to use default ribo-len='{ribo_len}', "
                f"or enter new comma-separated lengths: ").strip()
            if _new_len:
                if all(re.fullmatch(r"\d+", x.strip()) for x in _new_len.split(',')):
                    ribo_len = ",".join(x.strip() for x in _new_len.split(','))
                    tee.write(f"  Updated ribo-len: {ribo_len}\n")
                else:
                    tee.write(f"  Warning: invalid ribo-len '{_new_len}' "
                              f"(comma-separated positive integers expected), "
                              f"keeping {ribo_len}\n")
            _new_cut = input(
                f"  Press Enter to use default cutoffs='{cutoffs}', "
                f"or enter new comma-separated cutoffs: ").strip()
            if _new_cut:
                if all(re.fullmatch(r"\d+", x.strip()) for x in _new_cut.split(',')):
                    cutoffs = ",".join(x.strip() for x in _new_cut.split(','))
                    tee.write(f"  Updated cutoffs: {cutoffs}\n")
                else:
                    tee.write(f"  Warning: invalid cutoffs '{_new_cut}' "
                              f"(comma-separated positive integers expected), "
                              f"keeping {cutoffs}\n")

            # Validate matching counts
            _rl_n = len([x for x in ribo_len.split(',') if x.strip()])
            _co_n = len([x for x in cutoffs.split(',') if x.strip()])
            if _rl_n != _co_n:
                sys.exit(
                    f"ribo-len ({_rl_n} items) and cutoffs ({_co_n} items) "
                    f"must have the same number of values!\n"
                    f"  ribo-len:  {ribo_len}\n"
                    f"  cutoffs:   {cutoffs}"
                )
        except (EOFError, KeyboardInterrupt):
            tee.write(f"  Using defaults: ribo-len={ribo_len}, cutoffs={cutoffs}\n")

        # Run RIBO Taper
        tee.write("\n  Running RIBO Taper ORF detection...\n")
        ribotaper_script = os.path.join(ribotaper_path, "Ribotaper.sh")
        if not os.path.exists(ribotaper_script):
            sys.exit(f"RIBO Taper script not found: {ribotaper_script}")

        run_cmd(
            _ribo_bash(f"{ribotaper_script} {ribo_merged} {rna_merged} "
            f"{ribo_anno_dir} {ribo_len} {cutoffs} {thread}"))

        tee.write("\n  RIBO Taper analysis complete.\n")
        tee.write("  STEP 8 COMPLETE\n")

    # ── Plot Ribo-seq read length distributions ──
    # Generate when steps 1-8 are run (last_step < 8).
    # Skip when running step 9 only (last_step = 8, mapping already done).
    if last_step < 8:
        _len_pdf = "Ribo_length_distributions.pdf"
        if os.path.exists(_len_pdf):
            os.unlink(_len_pdf)
        _ribo_bam_files = [
            os.path.join(ribo_map_dir, f"star_{tag}_Aligned.sortedByCoord.out.bam")
            for tag in all_ribo_tags
        ]
        _ribo_bam_files = [b for b in _ribo_bam_files if os.path.exists(b)]
        if _ribo_bam_files:
            tee.write("\n  Plotting Ribo-seq read length distributions...\n")
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                from matplotlib.backends.backend_pdf import PdfPages
                _plot_ribo_lengths(ribo_map_dir, all_ribo_tags, tee)
            except ImportError:
                tee.write("  Warning: matplotlib not available, skipping length distribution plots\n")
        else:
            tee.write(f"\n  No Ribo-seq BAMs found in {ribo_map_dir}, skipping length distribution plots.\n")

    # ==================================================================
    # STEP 9 — P-site analysis & visualization
    # ==================================================================
    if last_step < 9:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 9: P-site analysis & visualization\n")
        tee.write("=" * 60 + "\n")

        # Cleanup previous incomplete step 9 outputs
        for pat in ["P_site_analysis.pdf", "P_site_frame_distribution.pdf",
                     "P_site_metagene.pdf"]:
            if os.path.exists(pat):
                os.unlink(pat)
        # Clean per-sample P-site split files
        for pat in ["P_sites_*.txt", "psite_*.txt"]:
            for fname in globmod.glob(pat):
                if os.path.exists(fname):
                    os.unlink(fname)
        # Extensionless per-sample outputs from previous runs. Exact names
        # only: P_sites_all* must survive — it is this step's input file.
        for tag in all_ribo_tags:
            _ps_out = f"P_sites_{tag}"
            if os.path.exists(_ps_out):
                os.unlink(_ps_out)

        _plot_psite_analysis(ribo_map_dir, all_ribo_tags, ribo_anno_dir, tee)
        tee.write("  STEP 9 COMPLETE\n")

    # ==================================================================
    # STEP 10 — Translation Efficiency (TE) calculation (per sample)
    # ==================================================================
    if last_step < 10:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 10: Translation Efficiency (TE) calculation\n")
        tee.write("=" * 60 + "\n")

        te_dir = "TE_results"
        os.makedirs(te_dir, exist_ok=True)

        # Build RSEM index if not already done (reuse from step 5)
        # Different RSEM versions lay out the STAR index differently: some
        # put it directly under RSEM_index/, others in RSEM_index/RNA/.
        rsem_index = "RSEM_index"
        _rsem_idx_markers = [
            os.path.join(rsem_index, "sjdbList.out.tab"),
            os.path.join(rsem_index, "RNA", "sjdbList.out.tab"),
        ]
        if not any(os.path.exists(p) for p in _rsem_idx_markers):
            tee.write("  Building RSEM index for TE calculation...\n")
            star_bin = shutil.which("STAR") or "STAR"
            run_cmd(
                f"rsem-prepare-reference --gtf {expressed_gtf} "
                f"--star --star-path {os.path.dirname(star_bin)} "
                f"--star-sjdboverhang 99 -p {thread} "
                f"{fasta_path} {rsem_index}/RNA", quiet=True)

        # Helper: detect paired-end status
        def _is_paired_bam(bam_path):
            try:
                r = subprocess.run(
                    f"samtools view -F 260 {bam_path} | head -1 | cut -f2",
                    shell=True, capture_output=True, text=True
                )
                flag = int(r.stdout.strip()) if r.stdout.strip() else 0
                return bool(flag & 1)
            except Exception:
                return False

        # ── Per-sample TE calculation ──
        # Pair Ribo-seq and RNA-seq samples by experimental role
        # (control↔control, treatment i↔treatment i), replicates positionally
        # within each group.
        te_pairs = _pair_te_samples(rna_groups, ribo_groups,
                                    all_rna_tags, all_ribo_tags, tee)
        if not te_pairs:
            tee.write("  Warning: No valid Ribo-seq/RNA-seq sample pairs, skipping TE\n")
        else:
            tee.write(f"  Found {len(te_pairs)} sample pair(s)\n")

            for i, (ribo_tag, rna_tag) in enumerate(te_pairs):
                pair_name = f"{ribo_tag}__{rna_tag}"
                pair_dir = os.path.join(te_dir, pair_name)
                os.makedirs(pair_dir, exist_ok=True)

                tee.write(f"\n  ── Pair {i+1}: {pair_name} ──\n")

                # Find Ribo-seq BAM
                ribo_bam_candidates = [
                    os.path.join(ribo_map_dir, f"star_{ribo_tag}_Aligned.toTranscriptome.out.bam"),
                    os.path.join(ribo_map_dir, f"{ribo_tag}_Aligned.toTranscriptome.out.bam"),
                ]
                ribo_bam = next((p for p in ribo_bam_candidates if os.path.exists(p)), None)

                # Find RNA-seq BAM
                rna_bam_candidates = [
                    os.path.join(star_rna_new_map, f"star_{rna_tag}_Aligned.toTranscriptome.out.bam"),
                    os.path.join(star_rna_new_map, f"{rna_tag}_Aligned.toTranscriptome.out.bam"),
                ]
                rna_bam = next((p for p in rna_bam_candidates if os.path.exists(p)), None)

                if not ribo_bam:
                    tee.write(f"    Warning: Ribo-seq BAM not found for {ribo_tag}, skipping\n")
                    continue
                if not rna_bam:
                    tee.write(f"    Warning: RNA-seq BAM not found for {rna_tag}, skipping\n")
                    continue

                tee.write(f"    Ribo-seq BAM: {ribo_bam}\n")
                tee.write(f"    RNA-seq BAM:  {rna_bam}\n")

                # Run RSEM on Ribo-seq BAM
                tee.write(f"    Running RSEM on {ribo_tag} (Ribo-seq)...\n")
                rsem_ribo_out = os.path.join(pair_dir, "ribo")
                ribo_pe = "--paired-end" if _is_paired_bam(ribo_bam) else ""
                run_cmd(
                    f"rsem-calculate-expression --alignments "
                    f"{ribo_pe} -p {thread} --time "
                    f"{os.path.abspath(ribo_bam)} "
                    f"{rsem_index}/RNA {rsem_ribo_out}", quiet=True)

                # Run RSEM on RNA-seq BAM
                tee.write(f"    Running RSEM on {rna_tag} (RNA-seq)...\n")
                rsem_rna_out = os.path.join(pair_dir, "rna")
                rna_pe = "--paired-end" if _is_paired_bam(rna_bam) else ""
                run_cmd(
                    f"rsem-calculate-expression --alignments "
                    f"{rna_pe} -p {thread} --time "
                    f"{os.path.abspath(rna_bam)} "
                    f"{rsem_index}/RNA {rsem_rna_out}", quiet=True)

                # Compute TE (saves te_table.tsv in pair_dir)
                tee.write(f"    Computing TE for {pair_name}...\n")
                try:
                    _compute_te(
                        os.path.join(pair_dir, "ribo.genes.results"),
                        os.path.join(pair_dir, "rna.genes.results"),
                        pair_dir, tee)
                except Exception as _te_err:
                    tee.write(f"    ERROR: TE computation failed for {pair_name}: "
                              f"{_te_err}\n")
                    tee.write("    Continuing with remaining pairs...\n")

        tee.write(f"\n  TE calculation complete. Results saved to: {te_dir}/\n")
        tee.write("  STEP 10 COMPLETE\n")

    # ==================================================================
    # STEP 11 — TE plotting (per-sample + combined summary)
    # ==================================================================
    if last_step < 11:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 11: TE plotting\n")
        tee.write("=" * 60 + "\n")

        te_dir = "TE_results"
        if not os.path.isdir(te_dir):
            tee.write("  No TE_results directory found, skipping TE plotting\n")
        else:
            # Load all TE tables from disk
            all_te_tables = {}
            for pair_name in sorted(os.listdir(te_dir)):
                pair_dir = os.path.join(te_dir, pair_name)
                if not os.path.isdir(pair_dir):
                    continue
                te_file = os.path.join(pair_dir, "te_table.tsv")
                if not os.path.exists(te_file):
                    tee.write(f"  Warning: No te_table.tsv for {pair_name}, skipping\n")
                    continue
                # Parse TSV back into list of dicts
                table = []
                with open(te_file) as f:
                    header = f.readline().strip().split('\t')
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        fields = line.split('\t')
                        row = {}
                        for h, val in zip(header, fields):
                            try:
                                row[h] = float(val) if val != 'NA' else float('nan')
                            except ValueError:
                                row[h] = val
                        table.append(row)

                if table:
                    all_te_tables[pair_name] = table
                    tee.write(f"  Loaded {pair_name}: {len(table)} genes\n")

            if not all_te_tables:
                tee.write("  No TE tables found, skipping plotting\n")
            else:
                # Plot per-sample TE
                for pair_name, te_table in all_te_tables.items():
                    pair_dir = os.path.join(te_dir, pair_name)
                    tee.write(f"  Plotting {pair_name}...\n")
                    _plot_te(te_table, pair_dir, tee)

                # Generate combined summary
                tee.write("  Generating combined TE summary...\n")
                _generate_te_summary(all_te_tables, te_dir, tee, rna_groups=rna_groups)

        tee.write("  STEP 11 COMPLETE\n")

    # ==================================================================
    # STEP 12 — Final output summary
    # ==================================================================
    if last_step < 12:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 12: Pipeline complete\n")
        tee.write("=" * 60 + "\n")

        tee.write(f"\n  Output files:\n")
        tee.write(f"    Updated GTF:          {updated_gtf}\n")
        tee.write(f"    Expressed GTF:        {expressed_gtf}\n")
        tee.write(f"    RSEM results:         {rsem_dir}/\n")
        tee.write(f"    RIBO Taper results:   (current directory)\n")
        tee.write(f"    RIBO Taper annotation:{ribo_anno_dir}/\n")
        tee.write(f"    TE results:           TE_results/\n")
        tee.write("  STEP 12 COMPLETE\n")

    # Cleanup contamination index
    for pat in [f"{contam_index}*.ebwt", f"{contam_index}*.bt2"]:
        for fname in globmod.glob(pat):
            if os.path.exists(fname):
                os.unlink(fname)

    tee.write("\nDone.\n")


# ── Helpers ──────────────────────────────────────────────────────────────


def _merge_or_copy_bams(bam_list, out_path, tee=None):
    """Merge BAMs into out_path; copy directly when there is only one input.

    `samtools merge` requires >= 2 input files on most versions, and an
    empty list fails outright, so handle both cases explicitly.
    """
    if not bam_list:
        sys.exit(f"  Error: no input BAMs to merge into {out_path}")
    if len(bam_list) == 1:
        if tee:
            tee.write(f"    Only one BAM ({bam_list[0]}), copying instead of merging\n")
        shutil.copy(bam_list[0], out_path)
    else:
        run_cmd(
            f"samtools merge -f {out_path} {' '.join(bam_list)}")
        sorted_tmp = out_path.replace('.bam', '_merged_sorted.bam')
        run_cmd(
            f"samtools sort -o {sorted_tmp} {out_path}")
        os.rename(sorted_tmp, out_path)
    run_cmd(f"samtools index {out_path}")


def _plot_ribo_lengths(ribo_map_dir, tags, tee):
    """Plot read length distribution for each Ribo-seq BAM into one PDF."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    bam_files = [
        os.path.join(ribo_map_dir, f"star_{tag}_Aligned.sortedByCoord.out.bam")
        for tag in tags
    ]
    bam_files = [(t, b) for t, b in zip(tags, bam_files) if os.path.exists(b)]

    if not bam_files:
        tee.write("  No Ribo-seq BAM files found for plotting.\n")
        return

    tee.write(f"  Found {len(bam_files)} BAM file(s), generating length distribution PDF...\n")

    pdf_path = "Ribo_length_distributions.pdf"
    with PdfPages(pdf_path) as pdf:
        # Page 1: combined histogram (all samples overlayed)
        fig, ax = plt.subplots(figsize=(10, 6))
        _all_lens = []
        for tag, bam in bam_files:
            lengths = _extract_read_lengths(bam)
            if lengths:
                _all_lens.extend(lengths)
                ax.hist(lengths, bins=range(min(lengths), max(lengths) + 2),
                        alpha=0.5, label=tag, density=True, align='left')
        ax.set_xlabel("Read Length (bp)")
        ax.set_ylabel("Density")
        ax.set_title("Ribo-seq Read Length Distribution (all samples)")
        ax.legend(fontsize=8)
        ax.set_xlim(0, 100)
        pdf.savefig(fig)
        plt.close(fig)

        # Per-sample pages
        for tag, bam in bam_files:
            lengths = _extract_read_lengths(bam)
            if not lengths:
                continue
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.hist(lengths, bins=range(min(lengths), max(lengths) + 2),
                    color='steelblue', edgecolor='white', density=True, align='left')
            ax.set_xlabel("Read Length (bp)")
            ax.set_ylabel("Density")
            ax.set_title(f"{tag} — Read Length Distribution")
            ax.set_xlim(0, 100)
            ax.axvline(x=28, color='red', linestyle='--', label='28 nt (typical Ribo-seq)')
            ax.legend()
            stats_text = (
                f"N reads: {len(lengths):,}\n"
                f"Mean: {sum(lengths)/len(lengths):.1f} bp\n"
                f"Median: {sorted(lengths)[len(lengths)//2]} bp"
            )
            ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
                    fontsize=9, verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            pdf.savefig(fig)
            plt.close(fig)

    tee.write(f"  Length distribution plot saved to: {pdf_path}\n")


def _extract_read_lengths(bam_path):
    """Extract read lengths from a BAM file using samtools.
    
    Uses SAM text output (not binary BAM) and filters to primary
    mapped reads (-F 260 = exclude unmapped + secondary alignments).
    """
    try:
        result = subprocess.run(
            f"samtools view -F 260 {bam_path} | awk '{{print length($10)}}'",
            shell=True, capture_output=True, text=True
        )
        lengths = [int(x) for x in result.stdout.strip().split('\n') if x]
        return lengths
    except Exception:
        return []


# ── Translation Efficiency (TE) helpers ─────────────────────────────────

def _pair_te_samples(rna_groups, ribo_groups, rna_tags, ribo_tags, tee):
    """Pair Ribo-seq/RNA-seq tags for TE calculation, aligned by role.

    Groups come from parse_input pars: [(group_name, n_replicates), ...]
    in input order (control first, then treatments).

    Strategy:
    - Same number of groups on both sides → pair k-th group with k-th group
      and pair replicates positionally within each group. Group names that
      differ between RNA and Ribo only produce a warning (role alignment is
      assumed from input order).
    - Different group structure → refuse to guess: report both structures
      and return no pairs.

    Returns:
        list of (ribo_tag, rna_tag) pairs.
    """
    if sum(n for _, n in rna_groups) != len(rna_tags) or \
            sum(n for _, n in ribo_groups) != len(ribo_tags):
        tee.write("  ERROR: group sizes do not add up to sample tag counts, "
                  "skipping TE\n")
        return []

    if len(rna_groups) != len(ribo_groups):
        tee.write("  ERROR: RNA-seq and Ribo-seq sample groups differ - "
                  "cannot safely pair samples for TE\n")
        tee.write(f"    RNA groups:  {rna_groups}\n")
        tee.write(f"    Ribo groups: {ribo_groups}\n")
        tee.write("    Please give control and treatments matching roles, e.g. "
                  "--rc ctrl --rt trt1 --bc ctrl --bt trt1\n")
        return []

    for (rna_name, _), (ribo_name, _) in zip(rna_groups, ribo_groups):
        if str(rna_name).lower() != str(ribo_name).lower():
            tee.write(f"  Warning: group names differ between RNA ('{rna_name}') "
                      f"and Ribo ('{ribo_name}') - pairing by input order\n")

    pairs = []
    rna_offset = 0
    ribo_offset = 0
    for g_idx, ((_, rna_n), (_, ribo_n)) in enumerate(zip(rna_groups, ribo_groups)):
        n_pair = min(rna_n, ribo_n)
        if rna_n != ribo_n:
            tee.write(f"  Warning: group {g_idx + 1} ({rna_groups[g_idx][0]}): "
                      f"{ribo_n} ribo vs {rna_n} rna samples - pairing first "
                      f"{n_pair} replicate(s)\n")
        for j in range(n_pair):
            pairs.append((ribo_tags[ribo_offset + j], rna_tags[rna_offset + j]))
        rna_offset += rna_n
        ribo_offset += ribo_n

    return pairs


def _compute_te(ribo_results_file, rna_results_file, output_dir, tee):
    """Compute Translation Efficiency (TE) = Ribo TPM / RNA TPM per gene.
    
    Reads RSEM gene-level results files for Ribo-seq and RNA-seq,
    matches genes, and computes TE ratios.
    
    Args:
        ribo_results_file: Path to RSEM ribo.genes.results
        rna_results_file: Path to RSEM rna.genes.results
        output_dir: Directory to write TE results
        tee: Tee writer for logging
    
    Returns:
        list of dicts with gene_id, ribo_tpm, rna_tpm, te, log2te or None
    """
    if not os.path.exists(ribo_results_file):
        tee.write(f"  Warning: {ribo_results_file} not found\n")
        return None
    if not os.path.exists(rna_results_file):
        tee.write(f"  Warning: {rna_results_file} not found\n")
        return None

    # Parse RSEM gene results
    # Format: gene_id, transcript_id(s), length, effective_length,
    #         expected_count, TPM, FPKM, IsoPct, ReprLeft, ...
    def parse_rsem_genes(filepath):
        """Parse RSEM .genes.results file into dict: gene_id -> {tpm, fpkm, count}"""
        genes = {}
        with open(filepath) as f:
            header = f.readline().strip().split('\t')
            # Find column indices
            col_idx = {h: i for i, h in enumerate(header)}
            for line in f:
                line = line.strip()
                if not line:
                    continue
                fields = line.split('\t')
                gene_id = fields[col_idx['gene_id']]
                tpm = float(fields[col_idx['TPM']]) if 'TPM' in col_idx else 0.0
                fpkm = float(fields[col_idx['FPKM']]) if 'FPKM' in col_idx else 0.0
                count = float(fields[col_idx['expected_count']]) if 'expected_count' in col_idx else 0.0
                genes[gene_id] = {'tpm': tpm, 'fpkm': fpkm, 'count': count}
        return genes

    tee.write(f"  Parsing RSEM results...\n")
    ribo_genes = parse_rsem_genes(ribo_results_file)
    rna_genes = parse_rsem_genes(rna_results_file)

    tee.write(f"    Ribo-seq genes: {len(ribo_genes):,}\n")
    tee.write(f"    RNA-seq genes: {len(rna_genes):,}\n")

    # Match genes and compute TE
    common_genes = set(ribo_genes.keys()) & set(rna_genes.keys())
    tee.write(f"    Common genes: {len(common_genes):,}\n")

    if not common_genes:
        tee.write("  Warning: No common genes between Ribo-seq and RNA-seq RSEM results\n")
        return None

    te_table = []
    for gene_id in sorted(common_genes):
        r_info = ribo_genes[gene_id]
        t_info = rna_genes[gene_id]
        
        ribo_tpm = r_info['tpm']
        rna_tpm = t_info['tpm']
        ribo_fpkm = r_info['fpkm']
        rna_fpkm = t_info['fpkm']
        ribo_count = r_info['count']
        rna_count = t_info['count']

        # Compute TE: use TPM ratio
        if rna_tpm > 0 and ribo_tpm > 0:
            te = ribo_tpm / rna_tpm
            log2te = math.log2(te)
        elif rna_tpm > 0:
            te = 0.0
            log2te = float('-inf')
        else:
            te = float('nan')
            log2te = float('nan')

        te_table.append({
            'gene_id': gene_id,
            'ribo_tpm': ribo_tpm,
            'rna_tpm': rna_tpm,
            'ribo_fpkm': ribo_fpkm,
            'rna_fpkm': rna_fpkm,
            'ribo_count': ribo_count,
            'rna_count': rna_count,
            'te': te,
            'log2te': log2te,
        })

    # Write TE table
    te_output = os.path.join(output_dir, "te_table.tsv")
    with open(te_output, 'w') as f:
        header = ['gene_id', 'ribo_tpm', 'rna_tpm', 'ribo_fpkm', 'rna_fpkm',
                  'ribo_count', 'rna_count', 'te', 'log2te']
        f.write('\t'.join(header) + '\n')
        for row in te_table:
            f.write('\t'.join([
                str(row['gene_id']),
                f"{row['ribo_tpm']:.4f}",
                f"{row['rna_tpm']:.4f}",
                f"{row['ribo_fpkm']:.4f}",
                f"{row['rna_fpkm']:.4f}",
                f"{row['ribo_count']:.1f}",
                f"{row['rna_count']:.1f}",
                f"{row['te']:.4f}" if not math.isnan(row['te']) else 'NA',
                f"{row['log2te']:.4f}" if not math.isnan(row['log2te']) and row['log2te'] != float('-inf') else 'NA',
            ]) + '\n')

    tee.write(f"  TE table written: {te_output} ({len(te_table):,} genes)\n")

    # Summary statistics
    log2tes = [r['log2te'] for r in te_table
               if not math.isnan(r['log2te']) and r['log2te'] != float('-inf')]
    if log2tes:
        import numpy as np
        arr = np.array(log2tes)
        tee.write(f"  TE summary (log2 scale):\n")
        tee.write(f"    Mean:   {np.mean(arr):.3f}\n")
        tee.write(f"    Median: {np.median(arr):.3f}\n")
        tee.write(f"    Std:    {np.std(arr):.3f}\n")
        tee.write(f"    Range:  [{np.min(arr):.3f}, {np.max(arr):.3f}]\n")
        # Percentiles
        p5, p25, p75, p95 = np.percentile(arr, [5, 25, 75, 95])
        tee.write(f"    P5: {p5:.3f}, P25: {p25:.3f}, P75: {p75:.3f}, P95: {p95:.3f}\n")

    return te_table


def _plot_te(te_table, output_dir, tee):
    """Plot Translation Efficiency (TE) distribution.
    
    Creates a PDF with:
    - Page 1: Histogram of log2(TE) with statistics
    - Page 2: Scatter plot of Ribo TPM vs RNA TPM (log-log) with TE contours
    
    Args:
        te_table: List of dicts from _compute_te
        output_dir: Directory for output PDF
        tee: Tee writer for logging
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        import numpy as np

        te_pdf = os.path.join(output_dir, "TE_distribution.pdf")
        log2tes = np.array([r['log2te'] for r in te_table
                           if not math.isnan(r['log2te']) and r['log2te'] != float('-inf')])

        if len(log2tes) == 0:
            tee.write("  Warning: No valid TE values to plot\n")
            return

        with PdfPages(te_pdf) as pdf:
            # Page 1: log2(TE) histogram
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.hist(log2tes, bins=80, color='steelblue', edgecolor='white',
                    density=True, alpha=0.85)
            
            # Add median line
            median_val = np.median(log2tes)
            ax.axvline(x=median_val, color='red', linestyle='--', linewidth=1.5,
                      label=f'Median: {median_val:.2f}')
            # Add mean line
            mean_val = np.mean(log2tes)
            ax.axvline(x=mean_val, color='green', linestyle=':', linewidth=1.5,
                      label=f'Mean: {mean_val:.2f}')
            # Add zero line (TE = 1, no translational change)
            ax.axvline(x=0, color='black', linestyle='-', linewidth=1, alpha=0.5,
                      label='TE = 1 (no change)')
            
            ax.set_xlabel('log₂(TE) = log₂(Ribo TPM / RNA TPM)')
            ax.set_ylabel('Density')
            ax.set_title(f'Translation Efficiency Distribution\n({len(log2tes):,} genes)')
            ax.legend(fontsize=8)
            
            # Add stats text box
            stats_text = (
                f"N genes: {len(log2tes):,}\n"
                f"Mean: {mean_val:.3f}\n"
                f"Median: {median_val:.3f}\n"
                f"Std: {np.std(log2tes):.3f}\n"
                f"Range: [{np.min(log2tes):.2f}, {np.max(log2tes):.2f}]"
            )
            ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
                    fontsize=9, verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            pdf.savefig(fig)
            plt.close(fig)

            # Page 2: Ribo TPM vs RNA TPM scatter
            valid_for_scatter = [r for r in te_table
                                if r['ribo_tpm'] > 0 and r['rna_tpm'] > 0]
            if valid_for_scatter:
                fig, ax = plt.subplots(figsize=(8, 8))
                
                ribo_tpms = np.array([r['ribo_tpm'] for r in valid_for_scatter])
                rna_tpms = np.array([r['rna_tpm'] for r in valid_for_scatter])
                log2tes_scatter = np.array([r['log2te'] for r in valid_for_scatter])
                
                # Scatter plot with color = TE
                sc = ax.scatter(np.log10(rna_tpms + 0.01), np.log10(ribo_tpms + 0.01),
                               c=log2tes_scatter, cmap='RdYlGn', alpha=0.4, s=5,
                               vmin=-3, vmax=3)
                plt.colorbar(sc, ax=ax, label='log₂(TE)')
                
                # Diagonal line (TE = 1)
                lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                       max(ax.get_xlim()[1], ax.get_ylim()[1])]
                ax.plot(lims, lims, 'k--', linewidth=0.8, alpha=0.5, label='TE = 1')
                
                ax.set_xlabel('log₁₀(RNA TPM + 0.01)')
                ax.set_ylabel('log₁₀(Ribo TPM + 0.01)')
                ax.set_title('Ribo TPM vs RNA TPM\n(colored by log₂ TE)')
                ax.legend(fontsize=8, loc='upper left')
                pdf.savefig(fig)
                plt.close(fig)

        tee.write(f"  TE plot saved: {te_pdf}\n")

    except ImportError as e:
        tee.write(f"  Warning: Missing required package for TE plotting ({e})\n")
    except Exception as e:
        tee.write(f"  Warning: TE plotting failed: {e}\n")


def _te_group_stats(wide_tsv, output_dir, tee, alpha=0.05,
                    rna_groups=None):
    """Statistical testing of TE across conditions via external R script.

    Design assumption: multiple conditions, each with biological replicates.
    Per-gene log2(TE) values are paired across conditions (same gene = unit).

    Runs scripts/TE_stats.R which:
      1. Fits linear model log2te ~ condition + gene_id (paired design)
      2. Uses emmeans for estimated marginal means + Tukey-adjusted contrasts
      3. Uses agricolae::HSD.test for compact letter display
      4. Outputs TE_emmeans.tsv, TE_contrasts.tsv, TE_letters.tsv,
         TE_summary_stats.tsv in output_dir

    Args:
        wide_tsv: Path to wide-format TE table (te_combined.tsv) with
                  columns: gene_id, {pair1}_log2te, {pair2}_log2te, ...
        output_dir: Directory for output files (R script writes here too)
        tee: Tee writer for logging
        alpha: Significance threshold (default 0.05)
        rna_groups: Optional list of [(group_name, n_replicates), ...] from
                    input parameters. If None, uses each pair as one group.

    Returns:
        (anova_f, anova_p, letters) — letters maps pair_name/condition ->
        letter string; empty dict when statistics could not be run.
    """
    letters = {}

    # Locate the R script relative to this file
    prefix = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    r_script = os.path.join(prefix, 'scripts', 'TE_stats.R')
    if not os.path.exists(r_script):
        tee.write(f"  ERROR: {r_script} not found - skipping R stats\n")
        return (None, None, letters)

    if not os.path.exists(wide_tsv):
        tee.write(f"  ERROR: {wide_tsv} not found - skipping R stats\n")
        return (None, None, letters)

    tee.write(f"  Using wide TE table for R stats: {wide_tsv}\n")

    # Group args: derive condition structure from rna_groups if provided,
    # otherwise treat each sample pair as its own group with 1 replicate.
    pair_names = []
    with open(wide_tsv) as f:
        header = f.readline().strip().split('\t')
        pair_names = [h[:-len('_log2te')] for h in header
                      if h.endswith('_log2te')]

    group_args = []
    if rna_groups and len(rna_groups) > 0:
        # Count how many pairs map to each group based on n_replicates
        col_idx = 0
        for gname, n_rep in rna_groups:
            take = min(n_rep, len(pair_names) - col_idx)
            if take <= 0:
                break
            # Use the first replicate's condition name as the group label
            label = str(gname).split('__')[0]  # strip ribo part if present
            group_args.append(f"{label}={take}")
            col_idx += take
        # Handle any leftover pairs
        while col_idx < len(pair_names):
            left = pair_names[col_idx]
            label = left.split('__')[0]
            remaining_same_label = sum(
                1 for p in pair_names[col_idx:] if p.split('__')[0] == label)
            group_args.append(f"{label}={remaining_same_label}")
            col_idx += remaining_same_label
    else:
        for p in pair_names:
            group_args.append(f"{p}=1")

    if len(group_args) < 2:
        tee.write("  Only 1 condition available, skipping statistical test\n")
        return (None, None, letters)

    cmd = (
        f'Rscript --vanilla "{r_script}" '
        f'"{wide_tsv}" "{output_dir}" "{alpha}" {" ".join(group_args)}'
    )
    run_cmd(cmd)

    # Parse back the letters result to annotate the boxplot
    letters_file = os.path.join(output_dir, "TE_letters.tsv")
    anova_p = None
    anova_f = None

    if os.path.exists(letters_file):
        try:
            with open(letters_file) as f:
                header = f.readline().strip().split('\t')
                cond_idx = header.index('condition') if 'condition' in header else 0
                letter_idx = header.index('letters') if 'letters' in header else 1
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) <= max(cond_idx, letter_idx):
                        continue
                    letters[fields[cond_idx]] = fields[letter_idx]
            tee.write(f"  Loaded significance letters: "
                      f"{', '.join(f'{k}={v}' for k, v in sorted(letters.items()))}\n")
        except Exception as e:
            tee.write(f"  Warning: Failed to parse {letters_file}: {e}\n")

    # Try parsing ANOVA p-value from log if R wrote a summary file
    stats_log = os.path.join(output_dir, "TE_anova.tsv")
    if os.path.exists(stats_log):
        try:
            with open(stats_log) as f:
                next(f)
                row = f.readline().strip().split('\t')
                if len(row) >= 3:
                    anova_f = float(row[1])
                    anova_p = float(row[2])
        except Exception:
            pass

    return (anova_f, anova_p, letters)


def _generate_te_summary(all_te_tables, output_dir, tee, rna_groups=None):
    """Generate combined TE summary across all sample pairs.
    
    Creates:
    1. Combined TE table (all pairs merged)
    2. Summary PDF with boxplot of log2(TE) per sample pair, annotated
       with significance letters (one-way ANOVA + Tukey HSD, shared letter
       = not significantly different at alpha = 0.05)
    3. Heatmap of TE statistics per sample pair
    
    Args:
        all_te_tables: dict of pair_name -> te_table (list of dicts)
        output_dir: Directory for output files
        tee: Tee writer for logging
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        import numpy as np

        if not all_te_tables:
            return

        pair_names = list(all_te_tables.keys())
        
        # Build wide-format table: one row per gene, columns per sample pair
        # Collect all gene_ids across all pairs
        all_genes = set()
        for table in all_te_tables.values():
            for row in table:
                all_genes.add(row['gene_id'])
        
        # Build lookup: gene_id -> pair_name -> data
        gene_pair_data = {}
        for pair_name, table in all_te_tables.items():
            for row in table:
                gid = row['gene_id']
                if gid not in gene_pair_data:
                    gene_pair_data[gid] = {}
                gene_pair_data[gid][pair_name] = row
        
        # Build header: gene_id + per-pair columns
        value_cols = ['ribo_tpm', 'rna_tpm', 'te', 'log2te']
        header = ['gene_id']
        for pn in pair_names:
            for vc in value_cols:
                header.append(f"{pn}_{vc}")
        
        # Write wide-format table
        combined_file = os.path.join(output_dir, "te_combined.tsv")
        n_genes = len(all_genes)
        with open(combined_file, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for gid in sorted(all_genes):
                parts = [gid]
                for pn in pair_names:
                    if gid in gene_pair_data and pn in gene_pair_data[gid]:
                        row = gene_pair_data[gid][pn]
                        for vc in value_cols:
                            val = row[vc]
                            if isinstance(val, float) and (math.isnan(val) or val == float('-inf')):
                                parts.append('NA')
                            else:
                                parts.append(f"{val:.4f}")
                    else:
                        parts.extend(['NA'] * len(value_cols))
                f.write('\t'.join(parts) + '\n')
        tee.write(f"  Combined TE table (wide format): {combined_file} ({n_genes:,} genes x {len(pair_names)} pairs)\n")

        # Summary PDF
        summary_pdf = os.path.join(output_dir, "TE_summary.pdf")
        with PdfPages(summary_pdf) as pdf:
            # Collect valid log2(TE) data per pair
            valid_data = {}
            for pair_name in pair_names:
                table = all_te_tables[pair_name]
                log2tes = [r['log2te'] for r in table
                          if not math.isnan(r['log2te']) and r['log2te'] != float('-inf')]
                if log2tes:
                    valid_data[pair_name] = np.array(log2tes)

            if valid_data:
                labels = list(valid_data.keys())
                n_pairs = len(labels)

                # Statistical test via external R script (emmeans + agricolae)
                _anova_f, anova_p, sig_letters = _te_group_stats(
                    combined_file, output_dir, tee, alpha=0.05,
                    rna_groups=rna_groups)

                # Page 1: Boxplot of log2(TE) per sample pair
                fig, ax = plt.subplots(figsize=(max(8, n_pairs * 1.5), 6))
                data_list = [valid_data[k] for k in labels]
                
                bp = ax.boxplot(data_list, labels=[l[:30] for l in labels],
                              patch_artist=True, widths=0.5)
                
                # Color boxes
                colors = plt.cm.Set2(np.linspace(0, 1, n_pairs))
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
                
                # Significance letters above each box (shared letter =
                # not significantly different by Tukey HSD)
                _ymin = min(float(v.min()) for v in valid_data.values())
                _ymax = max(float(v.max()) for v in valid_data.values())
                _pad = (_ymax - _ymin) * 0.08 if _ymax > _ymin else 1.0
                ax.set_ylim(_ymin - 0.2 * _pad, _ymax + 1.6 * _pad)
                for xi, lab in enumerate(labels):
                    ls = sig_letters.get(lab, '')
                    if ls:
                        ax.text(xi + 1, _ymax + 0.4 * _pad, ls,
                                ha='center', va='bottom', fontsize=12,
                                fontweight='bold', color='black')
                
                ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
                ax.set_ylabel('log₂(TE)')
                _anova_txt = (f' (ANOVA p = {anova_p:.2e})'
                              if anova_p is not None else '')
                _letter_txt = ('\n(letters: shared = not significant, '
                               'Tukey HSD α=0.05)') if sig_letters else ''
                ax.set_title(f'Translation Efficiency per Sample Pair{_anova_txt}'
                             f'{_letter_txt}')
                plt.xticks(rotation=45, ha='right', fontsize=7)
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

                # Page 2: Heatmap of mean log2(TE) per pair
                medians = [np.median(valid_data[k]) for k in valid_data]
                means = [np.mean(valid_data[k]) for k in valid_data]
                
                fig, ax = plt.subplots(figsize=(max(6, n_pairs * 1.2), 3))
                heat_data = np.array([means])  # 1 row: mean log2(TE) per pair
                im = ax.imshow(heat_data, aspect='auto', cmap='RdYlGn',
                             vmin=-2, vmax=2)
                ax.set_xticks(range(n_pairs))
                ax.set_xticklabels([l[:30] for l in labels], rotation=45, ha='right', fontsize=7)
                ax.set_yticks([0])
                ax.set_yticklabels(['Mean log₂(TE)'])
                ax.set_title('Mean Translation Efficiency per Sample Pair')
                plt.colorbar(im, ax=ax, label='log₂(TE)')
                
                # Add text annotations
                for j, (m, md) in enumerate(zip(means, medians)):
                    ax.text(j, 0, f'{m:.2f}\n(med:{md:.2f})',
                           ha='center', va='center', fontsize=7,
                           color='black' if abs(m) < 1.5 else 'white')
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

        tee.write(f"  TE summary plot saved: {summary_pdf}\n")

    except ImportError as e:
        tee.write(f"  Warning: Missing required package for TE summary ({e})\n")
    except Exception as e:
        tee.write(f"  Warning: TE summary generation failed: {e}\n")


def _read_psite_tracks(psite_file):
    """Read P_sites_all (BED14) file produced by RIBO Taper.
    
    Enforces one P-site per read (by read_name). When a read has multiple
    P-site entries, only the first is kept.
    
    BED14 format:
    1. chrom             2. chromStart (P-site)    3. chromEnd
    4. name (read_id)      5. score                     6. strand
    7. chrom (full read)   8. chromStart (full read)    9. RGB color
    10. block count        11. block sizes             12. block starts
    13. block ends         14. read length
    
    Frame is NOT a direct column — it must be computed by comparing
    the P-site position to the nearest start codon.
    """
    psites = []
    seen_reads = set()  # Track unique read names — each read gets exactly one P-site
    skipped_dupes = 0
    with open(psite_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            fields = line.split()
            if len(fields) < 7:
                continue
            # Extract key fields (P-site position: columns 1-3)
            read_name = fields[3]      # column 4: name/read_id
            
            # Enforce one P-site per read — skip duplicate read names
            if read_name in seen_reads:
                skipped_dupes += 1
                continue
            seen_reads.add(read_name)
            
            chrom = fields[0]          # column 1: chrom (P-site)
            pos = int(fields[1])       # column 2: chromStart (P-site position)
            chrom_end = int(fields[2]) if len(fields) > 2 and fields[2].lstrip('-').isdigit() else pos + 1
            strand = fields[5]         # column 6: strand
            
            # Full read coordinates (columns 7-8)
            full_read_chrom = fields[6] if len(fields) > 6 else chrom
            full_read_pos = int(fields[7]) if len(fields) > 7 and fields[7].lstrip('-').isdigit() else pos
            
            # RGB color (column 9)
            rgb = fields[8] if len(fields) > 8 else '0,0,0'
            
            # Block info (columns 10-13)
            block_count = int(fields[9]) if len(fields) > 9 and fields[9].lstrip('-').isdigit() else 0
            block_sizes = fields[10] if len(fields) > 10 else '0'
            block_starts = fields[11] if len(fields) > 11 else '0'
            block_ends = fields[12] if len(fields) > 12 else '0'
            
            # Read length (column 14)
            length = int(fields[13]) if len(fields) > 13 and fields[13].lstrip('-').isdigit() else 0
            
            # frame is not a direct column; set to 0, computed later if needed
            frame = 0
            
            entry = {
                'read_name': read_name,
                'chrom': chrom,
                'pos': pos,
                'chrom_end': chrom_end,
                'frame': frame,
                'length': length,
                'strand': strand,
                'full_read_chrom': full_read_chrom,
                'full_read_pos': full_read_pos,
                'rgb': rgb,
                'block_count': block_count,
                'block_sizes': block_sizes,
                'block_starts': block_starts,
                'block_ends': block_ends,
                'fields': fields,
                'line': line,
            }
            psites.append(entry)
    
    if skipped_dupes > 0:
        import sys
        print(f"  Note: Skipped {skipped_dupes} duplicate P-sites (same read had multiple entries)",
              file=sys.stderr)
    
    return psites, skipped_dupes


def _normalize_strand(strand):
    """Normalize strand to +/-."""
    s = strand.strip()
    if s in ('+', '-'):
        return s
    if s == '.':
        return '+'  # default
    return '+'


def _assign_frames(psites, transcr_exons_file, cds_coords_file, tmp_dir, tee=None):
    """Assign reading frames to P-sites using transcript coordinates.
    
    Uses cds_coords_transcripts (transcript-level CDS start/stop) and
    transcr_exons_ccds.bed (genomic exon positions) to compute frames.
    
    Steps:
    1. Write P-sites to temp BED file
    2. Filter transcr_exons_ccds.bed for .1 transcripts → temp BED
    3. Run bedtools intersect to assign transcript IDs to P-sites
    4. Build transcript exon lookup + read CDS coords
    5. For each P-site, convert genomic → transcript pos, compute frame
    
    Returns (valid_psites, debug_info).
    """
    # Step 1: Write P-sites to temp BED (name = index into psites list)
    psite_bed = os.path.join(tmp_dir, "psites_for_intersect.bed")
    with open(psite_bed, 'w') as f:
        for i, psite in enumerate(psites):
            chrom = psite['chrom']
            pos = psite['pos']  # BED0-based chromStart
            chrom_end = psite.get('chrom_end', pos + 1)
            strand = _normalize_strand(psite['strand'])
            f.write(f"{chrom}\t{pos}\t{chrom_end}\tpsite_{i}\t0\t{strand}\n")
    
    # Step 2: Use all transcript exons from transcr_exons_ccds.bed
    gene_bed = os.path.join(tmp_dir, "gene_exons_all.bed")
    with open(transcr_exons_file) as fin, open(gene_bed, 'w') as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            fields = line.split('\t')
            if len(fields) < 6:
                continue
            fout.write(line + '\n')
    
    # Step 3: Run bedtools intersect to assign transcript IDs
    def _af_warn(msg):
        if tee is not None:
            tee.write(msg)
        else:
            print(msg, file=sys.stderr)

    if shutil.which("bedtools") is None:
        _af_warn("  Warning: bedtools not found in PATH - frame assignment "
                 "will fail (all frames = -1)\n")

    psite_tx = {}  # psite_idx -> transcript_id
    try:
        intersect_cmd = f"bedtools intersect -a {psite_bed} -b {gene_bed} -s -wb"
        result = subprocess.run(intersect_cmd, shell=True, capture_output=True, text=True)
        n_primary = 0
        if result.returncode != 0:
            _af_warn(f"  Warning: bedtools intersect (-s) failed: "
                     f"{(result.stderr or '').strip()[:300]}\n")
        elif result.stdout:
            for line in result.stdout.strip().split('\n'):
                if not line.strip():
                    continue
                fields = line.split('\t')
                # -wb: 6 (psite A) + 6 (gene_bed B) = 12 columns
                # B columns: chrom, start, end, transcript_id (.1), gene_id, strand
                if len(fields) >= 12:
                    psite_idx = int(fields[3].replace('psite_', ''))
                    tx_id = fields[9]   # column 4 of B = transcript_id
                    if tx_id:
                        psite_tx[psite_idx] = tx_id
                        n_primary += 1
        
        # Retry without strand if no matches
        if n_primary == 0:
            if result.returncode == 0:
                _af_warn("  Warning: no strand-matched overlaps; retrying "
                         "bedtools intersect without -s\n")
            intersect_cmd = f"bedtools intersect -a {psite_bed} -b {gene_bed} -wb"
            result2 = subprocess.run(intersect_cmd, shell=True, capture_output=True, text=True)
            if result2.returncode != 0:
                _af_warn(f"  Warning: bedtools intersect fallback failed: "
                         f"{(result2.stderr or '').strip()[:300]}\n")
            elif result2.stdout:
                for line in result2.stdout.strip().split('\n'):
                    if not line.strip():
                        continue
                    fields = line.split('\t')
                    if len(fields) >= 12:
                        psite_idx = int(fields[3].replace('psite_', ''))
                        tx_id = fields[9]
                        if tx_id:
                            psite_tx[psite_idx] = tx_id
            else:
                _af_warn("  Warning: no P-sites overlapped transcript exons "
                         "(with or without strand matching)\n")
    except Exception as exc:
        _af_warn(f"  Warning: frame assignment via bedtools raised an error: "
                 f"{exc}\n")
    
    # Step 4: Build transcript exon lookup from filtered .1 BED
    tx_exons = {}  # transcript_id -> (starts[], ends[], strand)
    with open(gene_bed) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            gf = line.split('\t')
            if len(gf) < 6:
                continue
            tx_id = gf[3]
            tstart = int(gf[1])
            tend = int(gf[2])
            strand = gf[5]
            if tx_id not in tx_exons:
                tx_exons[tx_id] = [[], [], strand]
            tx_exons[tx_id][0].append(tstart)
            tx_exons[tx_id][1].append(tend)
    
    # Sort exon intervals per transcript
    for tx_id in tx_exons:
        tx_exons[tx_id][0].sort()
        tx_exons[tx_id][1].sort()
    
    # Step 5: Read CDS coordinates (transcript-level) for ALL transcripts
    # Format: transcript_id  start_codon  stop_codon
    cds_coords = {}  # transcript_id -> (cds_start, cds_stop)
    if cds_coords_file and os.path.exists(cds_coords_file):
        with open(cds_coords_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 3:
                    continue
                tx_id = fields[0]
                try:
                    cds_start = int(fields[1])
                    cds_stop = int(fields[2])
                    cds_coords[tx_id] = (cds_start, cds_stop)
                except ValueError:
                    continue
    
    # Step 6: Compute frames
    valid_psites = []
    for i, psite in enumerate(psites):
        chrom = psite['chrom']
        pos = psite['pos']  # BED0-based position
        strand = psite['strand']
        
        tx_id = psite_tx.get(i)
        if tx_id is None:
            psite['frame'] = -1
            continue
        
        # Get CDS start for this transcript
        if tx_id not in cds_coords:
            psite['frame'] = -1
            continue
        
        cds_start_tx, cds_stop = cds_coords[tx_id]
        
        # Get transcript exon structure
        if tx_id not in tx_exons:
            psite['frame'] = -1
            continue
        
        tx_starts, tx_ends, tx_strand = tx_exons[tx_id]
        
        # Convert genomic position → transcript position
        psite_tx_pos = _genomic_to_transcript(pos, tx_starts, tx_ends, tx_strand)
        
        # Compute frame (cds_coords is 1-based, psite_tx_pos is 0-based)
        frame = (psite_tx_pos - (cds_start_tx - 1)) % 3
        psite['frame'] = frame
        psite['psite_tx_pos'] = psite_tx_pos
        psite['tx_id'] = tx_id
        psite['cds_start'] = cds_start_tx
        psite['cds_stop'] = cds_stop
        valid_psites.append(psite)
    
    # Clean up temp BED files
    try:
        os.unlink(psite_bed)
        os.unlink(gene_bed)
    except Exception:
        pass
    
    return valid_psites, len(psite_tx)


def _genomic_to_transcript(genomic_pos, starts, ends, strand):
    """Convert a genomic position to transcript position given exon structure.
    
    For + strand: transcript starts at the first exon's start.
    For - strand: transcript starts at the last exon's end (reversed).
    """
    if strand == '+':
        cumulative = 0
        for s, e in zip(starts, ends):
            if genomic_pos < e:
                offset = max(0, genomic_pos - s)
                return cumulative + offset
            cumulative += (e - s)
        return cumulative  # beyond last exon
    else:
        # Minus strand: transcript is reversed
        # Iterate exons from right to left
        cumulative = 0
        for s, e in reversed(list(zip(starts, ends))):
            if genomic_pos >= s:
                offset = max(0, e - genomic_pos)
                return cumulative + offset
            cumulative += (e - s)
        return cumulative  # beyond first (leftmost) exon


def _read_transcript_lookups(transcr_exons_file, cds_coords_file):
    """Build per-transcript exon structure and CDS coordinate lookups.
    
    Uses ALL transcripts (no .1 filter).
    
    Returns:
        tx_exons: {transcript_id: (starts[], ends[], strand)}
        cds_coords: {transcript_id: (cds_start, cds_stop)}
    """
    tx_exons = {}  # transcript_id -> (starts[], ends[], strand)
    if transcr_exons_file and os.path.exists(transcr_exons_file):
        with open(transcr_exons_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                    continue
                fields = line.split('\t')
                if len(fields) < 6:
                    continue
                tx_id = fields[3]
                tstart = int(fields[1])
                tend = int(fields[2])
                strand = fields[5]
                if tx_id not in tx_exons:
                    tx_exons[tx_id] = [[], [], strand]
                tx_exons[tx_id][0].append(tstart)
                tx_exons[tx_id][1].append(tend)
        # Sort exon intervals per transcript
        for tx_id in tx_exons:
            tx_exons[tx_id][0].sort()
            tx_exons[tx_id][1].sort()
    
    cds_coords = {}  # transcript_id -> (cds_start, cds_stop)
    if cds_coords_file and os.path.exists(cds_coords_file):
        with open(cds_coords_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 3:
                    continue
                tx_id = fields[0]
                try:
                    cds_start = int(fields[1])
                    cds_stop = int(fields[2])
                    cds_coords[tx_id] = (cds_start, cds_stop)
                except ValueError:
                    continue
    
    return tx_exons, cds_coords


def _compute_frame_for_psite(psite_genomic_pos, transcript_id, tx_exons, cds_coords):
    """Compute reading frame for a single P-site given its transcript_id.
    
    Uses the transcript's exon structure and CDS coordinates.
    Returns frame (0, 1, or 2) or -1 if computation fails.
    """
    if transcript_id not in tx_exons or transcript_id not in cds_coords:
        return -1
    
    tx_starts, tx_ends, tx_strand = tx_exons[transcript_id]
    cds_start_tx, _ = cds_coords[transcript_id]
    
    # Convert genomic position → transcript position
    psite_tx_pos = _genomic_to_transcript(psite_genomic_pos, tx_starts, tx_ends, tx_strand)
    
    # Compute frame (cds_coords is 1-based, psite_tx_pos is 0-based)
    frame = (psite_tx_pos - (cds_start_tx - 1)) % 3
    return frame


def _split_psites_by_sample(psites, tags, bam_dir, tee):
    """Split P-sites by sample based on read name mapping.
    
    Uses the sample BAM files to identify which P-sites belong to which sample.
    Reads all mapped read names from each BAM, then matches P-sites by read name.
    Returns dict: tag -> list of psite dicts
    
    Reports unmatched P-sites (read_name not found in any sample BAM).
    """
    # First, extract read names from ALL samples
    tee.write("    Extracting read names from sample BAMs...\n")
    sample_read_sets = {}
    all_bam_reads = set()  # Union of all reads across all BAMs
    
    for tag in tags:
        # Check both STAR and Ribo BAM paths
        bam_candidates = [
            os.path.join(bam_dir, f"star_{tag}_Aligned.sortedByCoord.out.bam"),
            os.path.join(bam_dir, f"{tag}_Ribo.bam"),
            os.path.join(bam_dir, f"{tag}.bam"),
        ]
        bam_path = next((p for p in bam_candidates if os.path.exists(p)), None)
        
        if not bam_path:
            tee.write(f"    BAM not found for {tag}, skipping\n")
            sample_read_sets[tag] = set()
            continue
        
        tee.write(f"    Reading {tag} ({bam_path})...\n")
        try:
            result = subprocess.run(
                f"samtools view -F 260 {bam_path} | cut -f1",
                shell=True, capture_output=True, text=True
            )
            read_names = set(result.stdout.strip().split('\n')) - {'read_name'}
            sample_read_sets[tag] = read_names
            all_bam_reads.update(read_names)
            tee.write(f"    {tag}: {len(read_names):,} mapped reads\n")
        except Exception:
            tee.write(f"    Warning: Failed to read {bam_path}\n")
            sample_read_sets[tag] = set()
    
    # Detect read names present in more than one sample BAM. These make
    # per-sample assignment ambiguous (a P-site can only be given to the
    # first matching sample) and cause false dedup in _read_psite_tracks.
    from collections import Counter
    _name_owner_count = Counter()
    for rn_set in sample_read_sets.values():
        _name_owner_count.update(rn_set)
    cross_sample_names = {rn for rn, cnt in _name_owner_count.items() if cnt > 1}
    if cross_sample_names:
        n_cross = len(cross_sample_names)
        total_unique = len(_name_owner_count)
        tee.write(f"    Warning: {n_cross:,} read names appear in more than one "
                  f"sample BAM ({100 * n_cross / max(1, total_unique):.1f}% of "
                  f"unique BAM read names)\n")
        examples = sorted(cross_sample_names)[:3]
        if examples:
            tee.write(f"      e.g. {examples}\n")
        tee.write("      P-sites with these names are assigned to the first "
                  "matching sample only - use unique read names (e.g. SRA "
                  "accessions) for reliable per-sample splitting\n")

    # Build set of all P-site read names
    all_psite_reads = {p['read_name'] for p in psites}
    
    # Report matching stats
    n_unmatched_reads = len(all_psite_reads - all_bam_reads)
    n_shared_reads = len(all_psite_reads & all_bam_reads)
    tee.write(f"    Read name matching: {len(all_psite_reads):,} unique P-site reads, "
              f"{len(all_bam_reads):,} unique BAM reads\n")
    tee.write(f"    Matched: {n_shared_reads:,}, Unmatched (not in any BAM): {n_unmatched_reads:,}\n")
    if n_unmatched_reads > 0 and len(all_psite_reads) > 0:
        pct = 100 * n_unmatched_reads / len(all_psite_reads)
        tee.write(f"    Warning: {pct:.1f}% of P-site read names not found in any sample BAM\n")
        # Show a few examples of unmatched read names for debugging
        unmatched_examples = list(all_psite_reads - all_bam_reads)[:5]
        if unmatched_examples:
            tee.write(f"    Example unmatched read names: {unmatched_examples[:3]}\n")
    
    # Assign P-sites to samples
    sample_psites = {tag: [] for tag in tags}
    seen = {tag: set() for tag in tags}  # track unique keys per sample
    n_unassigned = 0
    
    tee.write(f"    Assigning {len(psites):,} P-sites to samples...\n")
    
    for psite in psites:
        rn = psite['read_name']
        # Deduplication key: chrom + pos + read_name + length + strand
        key = (psite['chrom'], psite['pos'], rn, psite['length'], psite['strand'])
        assigned = False
        for tag in tags:
            read_set = sample_read_sets.get(tag, set())
            if rn in read_set:
                if key not in seen[tag]:
                    seen[tag].add(key)
                    sample_psites[tag].append(psite)
                assigned = True
                break  # assign to first matching sample only
        if not assigned:
            n_unassigned += 1
    
    for tag in tags:
        n = len(sample_psites[tag])
        if n > 0:
            tee.write(f"    {tag}: assigned {n:,} P-sites\n")
    
    tee.write(f"    Total assigned: {sum(len(v) for v in sample_psites.values()):,}, "
              f"unassigned (read not in any BAM): {n_unassigned:,}\n")
    
    return sample_psites


# ── Metagene via bedtools closest (RiboTaper-style) ────────────────


def _run_metagene_bedtools(psite_file, start_stop_bed, tmp_dir, tee):
    """Run bedtools closest to map P-sites to nearest start/stop codons.
    
    Follows the RiboTaper approach (metag.R + create_metaplots.bash):
    - Input A: P_sites_all (BED14 format)
    - Input B: start_stop_FAR.bed (BED6 format with start/stop codon positions)
      Expected columns: chrom, start(0-based), end, name(start_codon/stop_codon),
                        score(gene_id), strand
    - Uses bedtools closest -s (strand-specific) -t "first" (first nearest)
    - Output: 14 + 6 = 20 columns (BED14 from A + BED6 from B)
    - Computes distance per RiboTaper formula:
        + strand: distance = psite_start - codon_start
        - strand: distance = codon_end - psite_end
    - Groups by read length for per-length metagene plotting
    
    Returns dict:
      {
        'by_length': {length: {'start': [distances], 'stop': [distances]}},  # legacy
        'raw_entries': [{'read_id': str, 'distance': int, 'codon_class': str, 'read_length': int}, ...],
        'total_psites': int,
        'total_matched': int,
      }
    """
    tee.write("  Running bedtools closest for metagene analysis...\n")
    
    # Step 1: Prepare P_sites_all as BED file (use directly)
    # The P_sites_all file is already in BED14 format, which bedtools can read.
    # We need to ensure it's sorted for bedtools closest.
    # Use version sort (-V) for chromosome names (chr1, chr2, ..., chr10)
    psite_sorted = os.path.join(tmp_dir, "psite_sorted.bed")
    sort_cmd = f"sort -k1,1V -k2,2n {psite_file} > {psite_sorted}"
    try:
        subprocess.run(sort_cmd, shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError:
        tee.write("  Warning: Failed to sort P_sites, using unsorted file\n")
        psite_sorted = psite_file
    
    # Step 2: Sort start_stop_FAR.bed
    ss_sorted = os.path.join(tmp_dir, "start_stop_sorted.bed")
    sort_cmd2 = f"sort -k1,1V -k2,2n {start_stop_bed} > {ss_sorted}"
    try:
        subprocess.run(sort_cmd2, shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError:
        tee.write("  Warning: Failed to sort start_stop file\n")
        ss_sorted = start_stop_bed
    
    # Step 3: Validate start_stop_FAR.bed format
    tee.write(f"  Validating start_stop file format...\n")
    sample_lines = []
    with open(ss_sorted) as f:
        for i, line in enumerate(f):
            if i >= 5:
                break
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('track'):
                sample_lines.append(line)
    
    if sample_lines:
        fields = sample_lines[0].split('\t')
        tee.write(f"  start_stop columns: {len(fields)}, sample: {fields[:6]}\n")
        if len(fields) < 4:
            tee.write("  Warning: start_stop_FAR.bed has < 4 columns, might not be valid BED format\n")
    
    # Step 4: Run bedtools closest
    closest_out = os.path.join(tmp_dir, "metagene_closest.bed")
    closest_cmd = (
        f"bedtools closest -a {psite_sorted} -b {ss_sorted} -s -t \"first\" > {closest_out}"
    )
    tee.write(f"  Command: {closest_cmd}\n")
    try:
        subprocess.run(closest_cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        tee.write(f"  Warning: bedtools closest failed: {e.stderr}\n")
        return None
    
    if not os.path.exists(closest_out) or os.path.getsize(closest_out) == 0:
        tee.write("  Warning: bedtools closest produced no output\n")
        return None
    
    tee.write(f"  bedtools closest output: {os.path.getsize(closest_out)} bytes\n")
    
    # Step 5: Parse the bedtools closest output
    # Format: BED14 (from A) + BED6 (from B) = 20 columns
    # Columns 1-14: from P_sites_all
    #   1: chrom, 2: start (P-site), 3: end, 4: read_id, 5: score, 6: strand
    #   7: chrom_full, 8: start_full, 9: RGB, 10: block_count
    #   11: block_sizes, 12: block_starts, 13: block_ends, 14: read_length
    # Columns 15-20: from start_stop_FAR.bed
    #   15: chrom, 16: start (codon), 17: end (codon)
    #   18: name (start_codon/stop_codon), 19: score (gene_id), 20: strand
    
    by_length = {}  # length -> {'start': [distances], 'stop': [distances]}
    raw_entries = []  # list of dicts with read_id, distance, codon_class, read_length
    total_psites = 0
    total_matched = 0
    
    with open(closest_out) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 20:
                continue
            
            total_psites += 1
            
            # Extract P-site info (cols 1-14)
            psite_start = int(fields[1]) if fields[1].lstrip('-').isdigit() else 0
            psite_end = int(fields[2]) if fields[2].lstrip('-').isdigit() else psite_start + 1
            read_id = fields[3]  # read name / ID
            psite_strand = fields[5]
            read_length = int(fields[13]) if len(fields) > 13 and fields[13].lstrip('-').isdigit() else 0
            
            # Extract codon info (cols 15-20 from start_stop_FAR.bed)
            # start_stop_FAR.bed format: chrom, start, end, name, score=transcript_id, strand
            codon_start = int(fields[15]) if fields[15].lstrip('-').isdigit() else 0
            codon_end = int(fields[16]) if fields[16].lstrip('-').isdigit() else codon_start + 3
            codon_type = fields[17]  # "start_codon" or "stop_codon" (BED name)
            transcript_id = fields[18] if len(fields) > 18 else ''  # BED score = transcript_id
            
            # Classify codon type (handle various naming conventions)
            codon_type_lower = codon_type.lower()
            if 'start' in codon_type_lower:
                codon_class = 'start'
            elif 'stop' in codon_type_lower or 'end' in codon_type_lower:
                codon_class = 'stop'
            else:
                continue  # Skip if codon type is not recognized
            
            # Compute distance (RiboTaper formula)
            # For + strand: distance = psite_start - codon_start
            # For - strand: distance = codon_end - psite_end
            if psite_strand == '+':
                distance = psite_start - codon_start
            else:
                distance = codon_end - psite_end
            
            total_matched += 1
            
            # Group by read length (legacy, keep for stats)
            if read_length not in by_length:
                by_length[read_length] = {'start': [], 'stop': []}
            by_length[read_length][codon_class].append(distance)
            
            # Store raw entry for per-sample/per-frame plotting
            raw_entries.append({
                'read_id': read_id,
                'distance': distance,
                'codon_class': codon_class,
                'read_length': read_length,
                'transcript_id': transcript_id,
                'psite_start': psite_start,
                'psite_strand': psite_strand,
            })
    
    tee.write(f"  Total P-sites in closest output: {total_psites}\n")
    tee.write(f"  Matched to start/stop codons: {total_matched}\n")
    tee.write(f"  Read lengths found: {sorted(by_length.keys())}\n")
    
    return {
        'by_length': by_length,
        'raw_entries': raw_entries,
        'total_psites': total_psites,
        'total_matched': total_matched,
    }


def _plot_psite_analysis(ribo_map_dir, tags, ribo_anno_dir, tee):
    """Run P-site analysis using RIBO Taper pre-computed outputs.
    
    Splits P-sites by sample, computes frame distribution and metagene plots.
    Outputs: Ribo_psite_analysis.pdf
    """
    tee.write("\n  Analyzing P-site distribution...\n")
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        import numpy as np
    except ImportError:
        tee.write("  Warning: matplotlib/numpy not available, skipping P-site analysis\n")
        return
    
    # Find RIBO Taper P-site output (check exact names first, then glob)
    psite_file = None
    candidates = [
        "P_sites_all",
        "P_sites_all.txt",
        "P_sites_all_tracks_exonsccds",
        "P_sites_all_tracks_exonsccds.txt",
        "P_sites_all_tracks",
        "P_sites",
    ]
    for c in candidates:
        if os.path.exists(c):
            psite_file = c
            break
    
    # Fallback: try glob for P_sites*
    if not psite_file:
        import glob as _glob
        _glob_candidates = sorted(_glob.glob("P_sites*"))
        # Prefer P_sites_all over other matches
        for g in _glob_candidates:
            if "P_sites_all" in g and not g.endswith(".txt"):
                psite_file = g
                break
        if not psite_file and _glob_candidates:
            psite_file = _glob_candidates[0]
    
    if not psite_file:
        tee.write(f"  Warning: RIBO Taper P-site file not found\n")
        return
    
    tee.write(f"  Reading P-site file: {psite_file}\n")
    psites, n_skipped_dupes = _read_psite_tracks(psite_file)
    if not psites:
        tee.write("  Warning: No P-sites found\n")
        return
    
    tee.write(f"  Found {len(psites)} P-sites total")
    if n_skipped_dupes > 0:
        tee.write(f" (skipped {n_skipped_dupes} duplicate P-sites from reads with multiple entries)")
    tee.write("\n")
    
    # Find start_stop_FAR.bed for bedtools closest approach
    start_stop_bed = None
    ss_candidates = [
        os.path.join(ribo_anno_dir, "start_stop_FAR.bed"),
        os.path.join(ribo_anno_dir, "start_stops_FAR.bed"),
        os.path.join(ribo_anno_dir, "start_stop.bed"),
        os.path.join(ribo_anno_dir, "start_stops.bed"),
        "start_stop_FAR.bed",
        "start_stops_FAR.bed",
    ]
    for ss in ss_candidates:
        if os.path.exists(ss):
            start_stop_bed = ss
            break
    
    # Also find transcr_exons_ccds.bed and cds_coords_transcripts for frame computation
    transcr_exons_file = None
    te_candidates = [
        os.path.join(ribo_anno_dir, "transcr_exons_ccds.bed"),
        os.path.join(ribo_anno_dir, "transcr_exons.bed"),
        "transcr_exons_ccds.bed",
        "transcr_exons.bed",
    ]
    for te in te_candidates:
        if os.path.exists(te):
            transcr_exons_file = te
            break
    
    cds_coords_file = None
    cds_candidates = [
        os.path.join(ribo_anno_dir, "cds_coords_transcripts"),
        os.path.join(ribo_anno_dir, "cds_coords_transcripts.txt"),
        "cds_coords_transcripts",
        "cds_coords_transcripts.txt",
    ]
    for cc in cds_candidates:
        if os.path.exists(cc):
            cds_coords_file = cc
            break
    
    # Compute frames if possible (for frame distribution plot)
    frame_valid_psites = []
    if transcr_exons_file and cds_coords_file:
        tee.write(f"  Loading transcript exon structures from: {transcr_exons_file}\n")
        tee.write(f"  Loading CDS coordinates from: {cds_coords_file}\n")
        tee.write("  Assigning frames (transcript-level coordinates, all transcripts)...\n")
        frame_valid_psites, n_matched = _assign_frames(
            psites, transcr_exons_file, cds_coords_file, ribo_map_dir, tee=tee)
        tee.write(f"  bedtools matched: {n_matched} P-sites\n")
        n_excluded = len(psites) - len(frame_valid_psites)
        frames_list = [p['frame'] for p in frame_valid_psites]
        tee.write(f"  Excluded {n_excluded} P-sites outside transcript exons\n")
        tee.write(f"  Valid P-sites (frame-assigned): {len(frame_valid_psites)}\n")
        if frames_list:
            fc = np.bincount(frames_list, minlength=3)[:3]
            tee.write(f"  Frame distribution: 0={fc[0]}, 1={fc[1]}, 2={fc[2]}\n")
    else:
        if not transcr_exons_file:
            tee.write("  Warning: transcr_exons_ccds.bed not found, skipping frame computation\n")
        if not cds_coords_file:
            tee.write("  Warning: cds_coords_transcripts not found, skipping frame computation\n")
        tee.write("  Frame assignment skipped, using all P-sites\n")
    
    # Split P-sites by sample (uses ALL psites, not just frame-valid ones)
    # _assign_frames modifies psite dicts in-place, so frame info is preserved
    tee.write("  Splitting P-sites by sample...\n")
    sample_psites = _split_psites_by_sample(
        psites, tags, ribo_map_dir, tee
    )
    
    # ── Run bedtools closest for metagene plotting (RiboTaper-style) ──
    metagene_data = None
    if start_stop_bed:
        tee.write(f"\n  Running bedtools closest with {start_stop_bed}...\n")
        metagene_data = _run_metagene_bedtools(
            psite_file, start_stop_bed, ribo_map_dir, tee)
        if metagene_data:
            tee.write(f"  Metagene data ready: {len(metagene_data['by_length'])} read lengths\n")
    else:
        tee.write("  Warning: start_stop_FAR.bed not found, skipping metagene plots\n")
    
    # Generate plots
    pdf_path = "Ribo_psite_analysis.pdf"
    tee.write(f"  Generating plots -> {pdf_path}\n")
    
    with PdfPages(pdf_path) as pdf:
        # ── Page 1: Frame distribution summary (all samples combined) ──
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Combined frame bar chart
        all_frames = []
        for tag, sp in sample_psites.items():
            for p in sp:
                if 'frame' in p and p['frame'] >= 0:
                    all_frames.append(p['frame'])
        
        if all_frames:
            frame_counts = np.bincount(all_frames, minlength=3)[:3]
            colors = ['#2ecc71', '#e74c3c', '#3498db']
            axes[0].bar(['Frame 0', 'Frame +1', 'Frame +2'],
                       frame_counts, color=colors, edgecolor='white')
            axes[0].set_ylabel('Count')
            axes[0].set_title('P-site Frame Distribution (All Samples)')
            total = frame_counts.sum()
            if total > 0:
                axes[0].legend([f'Frame {i}: {c} ({100*c/total:.1f}%)'
                              for i, c in enumerate(frame_counts)],
                             fontsize=8, loc='upper right')
            else:
                axes[0].legend([f'Frame {i}: {c}'
                              for i, c in enumerate(frame_counts)],
                             fontsize=8, loc='upper right')
        
        # Per-sample frame heatmap
        valid_samples = [t for t in tags if sample_psites.get(t)]
        if valid_samples:
            heatmap_data = np.zeros((3, len(valid_samples)))
            for j, tag in enumerate(valid_samples):
                frames = [p['frame'] for p in sample_psites[tag] if 'frame' in p and p['frame'] >= 0]
                if frames:
                    counts = np.bincount(frames, minlength=3)[:3]
                    total = counts.sum()
                    if total > 0:
                        heatmap_data[:, j] = counts / total
            
            im = axes[1].imshow(heatmap_data, aspect='auto', cmap='YlOrRd')
            axes[1].set_xticks(range(len(valid_samples)))
            axes[1].set_xticklabels(valid_samples, rotation=45, ha='right', fontsize=7)
            axes[1].set_yticks([0, 1, 2])
            axes[1].set_yticklabels(['Frame 0', 'Frame +1', 'Frame +2'])
            axes[1].set_title('Frame Proportion per Sample')
            plt.colorbar(im, ax=axes[1], label='Proportion')
        else:
            axes[1].text(0.5, 0.5, 'No P-sites assigned', ha='center', va='center')
            axes[1].set_title('Frame Proportion per Sample')
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        
        # ── Per-sample metagene plots (bedtools closest, RiboTaper-style) ──
        if metagene_data and metagene_data.get('raw_entries'):
            tee.write("\n  Preparing per-sample metagene data...\n")
            raw_entries = metagene_data['raw_entries']
            
            # Build read_id → sample_tag mapping from sample_psites
            # (for sample assignment only; frame is computed from transcript_id)
            read_to_sample = {}  # read_id -> tag
            for tag in tags:
                for psite_entry in sample_psites.get(tag, []):
                    rn = psite_entry.get('read_name', '')
                    if rn and rn not in read_to_sample:
                        read_to_sample[rn] = tag
            
            tee.write(f"  Read→sample mapping: {len(read_to_sample)} entries\n")
            
            # Build transcript lookups (ALL transcripts, no .1 filter)
            tx_exons, cds_coords = _read_transcript_lookups(transcr_exons_file, cds_coords_file)
            tee.write(f"  Transcript exon lookups: {len(tx_exons)} transcripts\n")
            tee.write(f"  CDS coordinate lookups: {len(cds_coords)} transcripts\n")
            
            # Compute frame for each raw_entry using its matched transcript_id
            # (from start_stop_FAR.bed match, not limited to .1)
            n_frame0 = n_frame1 = n_frame2 = n_unassigned = 0
            for entry in raw_entries:
                tid = entry.get('transcript_id', '')
                pstart = entry.get('psite_start', 0)
                if tid and tid in tx_exons and tid in cds_coords:
                    frame = _compute_frame_for_psite(pstart, tid, tx_exons, cds_coords)
                    entry['frame'] = frame
                    if frame == 0:
                        n_frame0 += 1
                    elif frame == 1:
                        n_frame1 += 1
                    elif frame == 2:
                        n_frame2 += 1
                    else:
                        n_unassigned += 1
                else:
                    entry['frame'] = -1
                    n_unassigned += 1
                # Assign sample from read_to_sample
                entry['sample'] = read_to_sample.get(entry.get('read_id', ''), None)
            
            tee.write(f"  Frame assignment: 0={n_frame0}, 1={n_frame1}, 2={n_frame2}, unassigned={n_unassigned}\n")
            
            # Define sample order: all tags that have data + "All" combined
            sample_order = ['All']
            for tag in tags:
                if sample_psites.get(tag):
                    sample_order.append(tag)
            
            # Initialize per-sample data structure
            # For each sample/codon: count[frame][distance]
            # distance range: -30..30 (61 positions)
            DIST_MIN, DIST_MAX = -30, 30
            N_DIST = DIST_MAX - DIST_MIN + 1  # 61
            
            # Structure: data[sample][codon] = [[counts for dist -30..30] for frame 0..2]
            # Plus 'unassigned' list for entries without frame/sample
            metagene_by_sample = {}  # sample -> {'start': [3 lists of 61], 'stop': [3 lists of 61]}
            for s in sample_order:
                metagene_by_sample[s] = {
                    'start': [np.zeros(N_DIST) for _ in range(3)],
                    'stop': [np.zeros(N_DIST) for _ in range(3)],
                    'unassigned_start': [],
                    'unassigned_stop': [],
                }
            
            # Process raw entries (frame already computed per-entry)
            n_assigned = 0
            n_unassigned_count = 0
            for entry in raw_entries:
                dist = entry['distance']
                codon = entry['codon_class']  # 'start' or 'stop'
                frame = entry.get('frame', -1)
                sample = entry.get('sample', None)
                
                # Filter to -30..30 range
                if dist < DIST_MIN or dist > DIST_MAX:
                    continue
                
                dist_idx = dist - DIST_MIN  # 0..60
                
                if 0 <= frame <= 2:
                    # Frame assigned: assign to "All" and specific sample
                    targets = ['All']
                    if sample and sample in metagene_by_sample:
                        targets.append(sample)
                    for s in targets:
                        if s in metagene_by_sample:
                            frame_counts = metagene_by_sample[s][codon][frame]
                            frame_counts[dist_idx] += 1
                    n_assigned += 1
                else:
                    # Frame unassigned
                    key = 'unassigned_start' if codon == 'start' else 'unassigned_stop'
                    metagene_by_sample['All'][key].append(dist)
                    if sample and sample in metagene_by_sample:
                        metagene_by_sample[sample][key].append(dist)
                    n_unassigned_count += 1
            
            tee.write(f"  Assigned (frame+sample): {n_assigned}, unassigned: {n_unassigned}\n")
            
            dist_positions = np.arange(DIST_MIN, DIST_MAX + 1)  # -30..30
            
            def _draw_sample_page(ax_row, sample_name, sample_data):
                """Draw a 2x2 metagene page for one sample.
                ax_row[0,0] = Start Codon (line plot, per-frame)
                ax_row[0,1] = Stop Codon (line plot, per-frame)
                ax_row[1,0] = Overlay (Start vs Stop total line plot)
                ax_row[1,1] = Stats
                """
                # ── Start Codon (line plot only) ──
                ax_start = ax_row[0, 0]
                _draw_codon_line_panel(ax_start, dist_positions, sample_data, 'start',
                                       f'Start Codon — {sample_name}')
                
                # ── Stop Codon (line plot only) ──
                ax_stop = ax_row[0, 1]
                _draw_codon_line_panel(ax_stop, dist_positions, sample_data, 'stop',
                                       f'Stop Codon — {sample_name}')
                
                # ── Overlay (total counts line plot) ──
                ax_overlay = ax_row[1, 0]
                start_total = np.zeros(N_DIST)
                stop_total = np.zeros(N_DIST)
                for f in range(3):
                    start_total += sample_data['start'][f]
                    stop_total += sample_data['stop'][f]
                # Add unassigned
                for d in sample_data.get('unassigned_start', []):
                    if DIST_MIN <= d <= DIST_MAX:
                        start_total[d - DIST_MIN] += 1
                for d in sample_data.get('unassigned_stop', []):
                    if DIST_MIN <= d <= DIST_MAX:
                        stop_total[d - DIST_MIN] += 1
                
                ax_overlay.plot(dist_positions, start_total, color='#3498db',
                               linewidth=2, marker='o', markersize=3, label='Start', alpha=0.85)
                ax_overlay.plot(dist_positions, stop_total, color='#e74c3c',
                               linewidth=2, marker='s', markersize=3, label='Stop', alpha=0.85)
                ax_overlay.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
                ax_overlay.set_xlabel('Distance from Codon (nt)')
                ax_overlay.set_ylabel('Count')
                ax_overlay.set_title(f'Overlay — {sample_name}')
                ax_overlay.set_xlim(DIST_MIN, DIST_MAX)
                ax_overlay.legend(fontsize=8, loc='upper right')
                
                # ── Stats ──
                ax_stats = ax_row[1, 1]
                ax_stats.axis('off')
                n_start_total = int(start_total.sum())
                n_stop_total = int(stop_total.sum())
                lines = [
                    f'Sample: {sample_name}',
                    f'',
                    f'Start codon (in range): {n_start_total}',
                    f'  Frame 0: {int(sample_data["start"][0].sum())}',
                    f'  Frame 1: {int(sample_data["start"][1].sum())}',
                    f'  Frame 2: {int(sample_data["start"][2].sum())}',
                    f'',
                    f'Stop codon (in range): {n_stop_total}',
                    f'  Frame 0: {int(sample_data["stop"][0].sum())}',
                    f'  Frame 1: {int(sample_data["stop"][1].sum())}',
                    f'  Frame 2: {int(sample_data["stop"][2].sum())}',
                ]
                ax_stats.text(0.1, 0.95, '\n'.join(lines), fontsize=10,
                             verticalalignment='top', fontfamily='monospace',
                             transform=ax_stats.transAxes)
                ax_stats.set_title(f'Stats — {sample_name}')
            
            def _draw_codon_line_panel(ax, dist_pos, sample_data, codon_type, title):
                """Draw line plot (total only) for one codon type.
                No frame-specific lines — total counts only.
                """
                frame_data_list = sample_data[codon_type]  # list of 3 arrays
                
                # Total line (sum of all frames + unassigned)
                total_counts = sum(frame_data_list)
                unassigned_key = 'unassigned_' + codon_type
                for d in sample_data.get(unassigned_key, []):
                    if DIST_MIN <= d <= DIST_MAX:
                        total_counts[d - DIST_MIN] += 1
                if total_counts.sum() > 0:
                    ax.plot(dist_pos, total_counts, color='#2c3e50',
                           linewidth=1.8, marker='.', markersize=3,
                           alpha=0.85, label='Total')
                
                ax.axvline(x=0, color='red', linestyle='--', linewidth=1, alpha=0.6)
                ax.set_xlabel('Distance from Codon (nt)')
                ax.set_ylabel('Count')
                n_total = int(total_counts.sum())
                ax.set_title(f'{title} (n={n_total})')
                ax.set_xlim(DIST_MIN, DIST_MAX)
                ax.set_xticks(range(DIST_MIN, DIST_MAX + 1, 5))
                ax.legend(fontsize=7, loc='upper right')
            
            # ── Generate all pages ──
            pages_generated = 0
            for sample_name in sample_order:
                sample_data = metagene_by_sample[sample_name]
                # Check if there's any data for this sample
                has_data = (
                    sample_data['start'][0].sum() + sample_data['start'][1].sum() + sample_data['start'][2].sum() +
                    sample_data['stop'][0].sum() + sample_data['stop'][1].sum() + sample_data['stop'][2].sum() +
                    len(sample_data.get('unassigned_start', [])) +
                    len(sample_data.get('unassigned_stop', []))
                ) > 0
                
                if not has_data:
                    tee.write(f"  Skipping {sample_name}: no metagene data\n")
                    continue
                
                fig, axes_2x2 = plt.subplots(2, 2, figsize=(14, 10))
                fig.suptitle(
                    f'Metagene Profile — {sample_name}',
                    fontsize=14, fontweight='bold')
                
                _draw_sample_page(axes_2x2, sample_name, sample_data)
                
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                pages_generated += 1
            
            tee.write(f"  Generated {pages_generated} metagene pages (All + per-sample)\n")
        else:
            tee.write("  No metagene data available, skipping metagene plots\n")
    
    # Save per-sample P-site files (same format as P_sites_all + frame column)
    tee.write("\n  Saving per-sample P-site files...\n")
    for tag, assignments in sample_psites.items():
        if not assignments:
            continue
        out_file = f"P_sites_{tag}"
        with open(out_file, 'w') as f:
            for p in assignments:
                line = p.get('line', '')
                frame = p.get('frame', -1)
                # Append frame as last column (-1, 0, 1, or 2)
                f.write(f"{line}\t{frame}\n")
        tee.write(f"    {out_file}: {len(assignments)} P-sites\n")
    
    # Clean up temp files from bedtools closest
    for tmp_f in ['psite_sorted.bed', 'start_stop_sorted.bed', 'metagene_closest.bed']:
        tmp_path = os.path.join(ribo_map_dir, tmp_f)
        if os.path.exists(tmp_path):
            try:
                os.unlink(tmp_path)
            except Exception:
                pass
    
    tee.write(f"\n  P-site analysis complete. Plots saved to: {pdf_path}\n")