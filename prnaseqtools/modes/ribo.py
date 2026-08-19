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

Requires: Bowtie2, STAR, samtools, bedtools, StringTie, gffcompare, RSEM,
          RIBO Taper scripts, R (dplyr)
"""

import os
import sys
import re
import glob as globmod
import shutil
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input, _parse_to_dict
from prnaseqtools.functions import download_sra, unzip_file, _tee, run_cmd


# ── Helper: detect last completed step from log ────────────────────────────

def _detect_last_step():
    """Detect the last completed step from log files in the current directory.

    Returns:
        tuple: (last_step, log_path) where last_step is int (0 if none found),
               log_path is the file path or None
    """
    log_files = globmod.glob(os.path.join(os.getcwd(), "log_*.txt"))
    if not log_files:
        return (0, None)

    # Find the most recently modified log file
    log_files.sort(key=os.path.getmtime, reverse=True)
    latest_log = log_files[0]

    last_step = 0
    with open(latest_log) as f:
        for line in f:
            match = re.search(r'STEP\s+(\d+):', line)
            if match:
                step = int(match.group(1))
                if step > last_step:
                    last_step = step

    return (last_step, latest_log)


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
    ribo_len = opts.get('ribo_len', '24,25,26,27,28,29,30,31')
    cutoffs = opts.get('cutoffs', '8,9,10,11,12,13,14,15')
    tpm_threshold = opts.get('tpm_threshold', 0)
    ribotaper_path = opts.get('ribotaper')
    ribotaper_env = opts.get('ribotaper_env')

    def _ribo_bash(cmd):
        """Wrap a bash command with conda run if ribotaper_env is set."""
        if ribotaper_env:
            return f"conda run -n {ribotaper_env} bash {cmd}"
        return f"bash {cmd}"

    # ── Parse RNA-seq inputs ─────────────────────────────────────────
    rna_ctrl_dict = _parse_to_dict(opts.get('rna_control', ''))
    rna_tags, rna_files, _ = parse_input(rna_ctrl_dict)

    rna_trt = opts.get('rna_treatment')
    if rna_trt:
        if isinstance(rna_trt, str):
            rna_trt = [rna_trt]
        for trt_str in rna_trt:
            trt_dict = _parse_to_dict(trt_str)
            rt_tags, rt_files, _ = parse_input(trt_dict)
            rna_tags.extend(rt_tags)
            rna_files.extend(rt_files)

    all_rna_tags = list(rna_tags)

    # ── Parse Ribo-seq inputs ────────────────────────────────────────
    ribo_ctrl_dict = _parse_to_dict(opts.get('ribo_control', ''))
    ribo_tags, ribo_files, _ = parse_input(ribo_ctrl_dict)

    ribo_trt = opts.get('ribo_treatment')
    if ribo_trt:
        if isinstance(ribo_trt, str):
            ribo_trt = [ribo_trt]
        for rbt_str in ribo_trt:
            rbt_dict = _parse_to_dict(rbt_str)
            rbt_tags, rbt_files, _ = parse_input(rbt_dict)
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

    # ── Auto-detect last completed step from log ─────────────────────
    last_step, resume_log = _detect_last_step()
    if last_step > 0:
        tee.write(f"\n  Resuming from step {last_step + 1} (last completed: step {last_step})\n")
        tee.write(f"  Log file: {resume_log}\n\n")
    else:
        tee.write("  No previous log found, starting from step 1.\n\n")

    # ── Verify intermediate files if resuming from step 6+ ──
    if last_step >= 6:
        if not os.path.exists(expressed_gtf):
            tee.write(f"  WARNING: {expressed_gtf} not found - steps 4-5 may need to be re-run\n")
        if not os.path.exists(updated_gtf):
            tee.write(f"  WARNING: {updated_gtf} not found - steps 4-5 may need to be re-run\n")

    # ==================================================================
    # STEP 1 — Build Bowtie2 contamination index
    # ==================================================================
    if last_step < 1:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 1: Building Bowtie2 contamination index\n")
        tee.write("=" * 60 + "\n")

        run_cmd(
            f"bowtie2-build {contam} {contam_index}")
        tee.write("  Contamination index built.\n")

    # ==================================================================
    # STEP 2 — Preprocess Ribo-seq reads (decontamination)
    # ==================================================================
    if last_step < 2:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 2: Preprocessing Ribo-seq reads\n")
        tee.write("=" * 60 + "\n")

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

    # ==================================================================
    # STEP 3 — RNA-seq STAR mapping + StringTie transcriptome assembly
    # ==================================================================
    if last_step < 3:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 3: RNA-seq transcriptome assembly (STAR + StringTie)\n")
        tee.write("=" * 60 + "\n")

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
        assembled_gtf_dir = "assembledGTF"
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
        gffcomp_prefix = os.path.join(assembled_gtf_dir, genome)
        run_cmd(
            f"gffcompare -V -r {anno_path} -o {gffcomp_prefix} {merged_gtf}")

    # ==================================================================
    # STEP 4 — R: Filter novel transcripts + add gene_biotype
    # ==================================================================
    if last_step < 4:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 4: Filtering novel transcripts + gene_biotype annotation\n")
        tee.write("=" * 60 + "\n")

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

    # ==================================================================
    # STEP 5 — RSEM quantification + filter expressed isoforms
    # ==================================================================
    if last_step < 5:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 5: RSEM quantification + expressed isoform filtering\n")
        tee.write("=" * 60 + "\n")

        # Build RSEM index
        rsem_index = "RSEM_index"
        os.makedirs(rsem_index, exist_ok=True)
        star_bin = shutil.which("STAR") or "STAR"

        run_cmd(
            f"rsem-prepare-reference --gtf {updated_gtf} "
            f"--star --star-path {os.path.dirname(star_bin)} "
            f"--star-sjdboverhang 99 -p {thread} "
            f"{fasta_path} {rsem_index}/RNA")

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
                f"{read_files} {rsem_index}/RNA {rsem_dir}/{tag}")

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

    # ==================================================================
    # STEP 6 — STAR re-mapping with expressed annotation
    # ==================================================================
    if last_step < 6:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 6: STAR re-mapping with expressed annotation\n")
        tee.write("=" * 60 + "\n")

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
        ribo_map_dir = "STAR_Ribo_map"
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

    # ── Plot Ribo-seq read length distributions (always runs when BAMs exist) ──
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

    # ==================================================================
    # STEP 7 — RIBO Taper annotation
    # ==================================================================
    if last_step < 7:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 7: Creating RIBO Taper annotation files\n")
        tee.write("=" * 60 + "\n")

        ribo_anno_dir = "RiboTaper_annotation"
        os.makedirs(ribo_anno_dir, exist_ok=True)

        create_anno_script = os.path.join(
            ribotaper_path, "create_annotations_files.bash"
        )
        if not os.path.exists(create_anno_script):
            sys.exit(f"RIBO Taper script not found: {create_anno_script}")

        from prnaseqtools.reference import check_and_build_indices
        check_and_build_indices(prefix, genome=genome, tee=tee)

        run_cmd(
            _ribo_bash(f"{create_anno_script} {expressed_gtf} {fasta_path} "
            f"false false {ribo_anno_dir} {bedtools_bin} "
            f"{ribotaper_path}/"))
        tee.write("  RIBO Taper annotation created.\n")

    # ==================================================================
    # STEP 8 — Merge BAMs + Run RIBO Taper
    # ==================================================================
    if last_step < 8:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 8: Merging BAMs + Running RIBO Taper\n")
        tee.write("=" * 60 + "\n")

        # Merge RNA-seq BAMs
        tee.write("\n  Merging RNA-seq BAMs...\n")
        rna_bam_list = [
            f"{star_rna_new_map}/star_{tag}_Aligned.sortedByCoord.out.bam"
            for tag in all_rna_tags
        ]
        rna_bam_list = [b for b in rna_bam_list if os.path.exists(b)]
        rna_bams_arg = ' '.join(rna_bam_list)

        rna_merged = "RNA_merged.bam"
        run_cmd(
            f"samtools merge -f {rna_merged} {rna_bams_arg}")
        run_cmd(
            f"samtools sort -o RNA_merged_sorted.bam {rna_merged}")
        os.rename("RNA_merged_sorted.bam", rna_merged)
        run_cmd(f"samtools index {rna_merged}")

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
        ribo_bams_arg = ' '.join(ribo_r1_bams)
        run_cmd(
            f"samtools merge -f {ribo_merged} {ribo_bams_arg}")
        run_cmd(
            f"samtools sort -o Ribo_merged_sorted.bam {ribo_merged}")
        os.rename("Ribo_merged_sorted.bam", ribo_merged)
        run_cmd(f"samtools index {ribo_merged}")

        # Run RIBO Taper
        tee.write("\n  Running RIBO Taper ORF detection...\n")
        ribotaper_script = os.path.join(ribotaper_path, "Ribotaper.sh")
        if not os.path.exists(ribotaper_script):
            sys.exit(f"RIBO Taper script not found: {ribotaper_script}")

        run_cmd(
            _ribo_bash(f"{ribotaper_script} {ribo_merged} {rna_merged} "
            f"{ribo_anno_dir} {ribo_len} {cutoffs} "
            f"{ribotaper_path}/ {bedtools_bin} {thread}"))

        tee.write("\n  RIBO Taper analysis complete.\n")

    # ==================================================================
    # STEP 9 — P-site analysis & visualization
    # ==================================================================
    if last_step < 9:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 9: P-site analysis & visualization\n")
        tee.write("=" * 60 + "\n")

        _plot_psite_analysis(ribo_map_dir, all_ribo_tags, ribo_anno_dir, ribo_files, tee, thread)

    # ==================================================================
    # STEP 10 — Final output summary
    # ==================================================================
    if last_step < 10:
        tee.write("\n" + "=" * 60 + "\n")
        tee.write("STEP 10: Pipeline complete\n")
        tee.write("=" * 60 + "\n")

        tee.write(f"\n  Output files:\n")
        tee.write(f"    Updated GTF:          {updated_gtf}\n")
        tee.write(f"    Expressed GTF:        {expressed_gtf}\n")
        tee.write(f"    RSEM results:         {rsem_dir}/\n")
        tee.write(f"    RIBO Taper results:   (current directory)\n")
        tee.write(f"    RIBO Taper annotation:{ribo_anno_dir}/\n")

    # Cleanup contamination index
    for pat in [f"{contam_index}*.ebwt", f"{contam_index}*.bt2"]:
        for fname in globmod.glob(pat):
            if os.path.exists(fname):
                os.unlink(fname)

    tee.write("\nDone.\n")


# ── Helpers ──────────────────────────────────────────────────────────────


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
        for tag, bam in bam_files:
            lengths = _extract_read_lengths(bam)
            if lengths:
                ax.hist(lengths, bins=50, alpha=0.5, label=tag, density=True)
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
            ax.hist(lengths, bins=50, color='steelblue', edgecolor='white', density=True)
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
    """Extract read lengths from a BAM file using samtools."""
    import subprocess
    try:
        result = subprocess.run(
            f"samtools view -b {bam_path} | awk '{{print length($10)}}'",
            shell=True, capture_output=True, text=True
        )
        lengths = [int(x) for x in result.stdout.strip().split('\n') if x]
        return lengths
    except Exception:
        return []


# ── P-site analysis helpers (using RIBO Taper output) ────────────────


def _read_psite_tracks(psite_file):
    """Read RIBO Taper P_sites_all_tracks_exonsccds file.
    
    Returns list of dicts with keys: read_name, chrom, pos, frame, length, etc.
    Auto-detects delimiter (tab or space).
    """
    psites = []
    with open(psite_file) as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 4:
                continue
            entry = {
                'read_name': fields[0],
                'chrom': fields[1],
                'pos': int(fields[2]),
            }
            if len(fields) >= 5:
                entry['frame'] = int(fields[3]) if fields[3].lstrip('-').isdigit() else 0
            if len(fields) >= 6:
                entry['length'] = int(fields[4]) if fields[4].lstrip('-').isdigit() else 0
            if len(fields) >= 7:
                entry['strand'] = fields[5]
            psites.append(entry)
    return psites


def _read_start_stops_far(bed_file):
    """Read start_stops_FAR.bed from RIBO Taper reference.

    Extracts only the first transcript (transcript .1) per gene.
    Returns dict: chrom -> {gene: (start, end)} for starts and stops.
    """
    starts = {}
    stops = {}
    with open(bed_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            transcript_info = fields[3]

            if chrom not in starts:
                starts[chrom] = {}
                stops[chrom] = {}

            # Only keep first transcript (.1): e.g. AT1G01010.1_start / AT1G01010.1_stop
            parts = transcript_info.rsplit('_', 1)
            tx_name = parts[0]
            if not tx_name.endswith('.1'):
                continue

            gene = tx_name.split('.')[0]

            if transcript_info.endswith('_start'):
                starts[chrom][gene] = (start, end)
            elif transcript_info.endswith('_stop'):
                stops[chrom][gene] = (start, end)

    return starts, stops


def _get_sample_read_names(fastq_path, max_reads=100000):
    """Get read names from a FASTQ/FASTQ.gz file to identify sample origin.
    
    Returns set of read name prefixes (first part of @readname before space).
    """
    import gzip
    
    read_names = set()
    try:
        if fastq_path.endswith('.gz'):
            with gzip.open(fastq_path, 'rt') as f:
                for i, line in enumerate(f):
                    if i % 4 == 0:
                        name = line.strip()[1:].split()[0]
                        read_names.add(name)
                    if len(read_names) >= max_reads:
                        break
        else:
            with open(fastq_path) as f:
                for i, line in enumerate(f):
                    if i % 4 == 0:
                        name = line.strip()[1:].split()[0]
                        read_names.add(name)
                    if len(read_names) >= max_reads:
                        break
    except Exception:
        pass
    return read_names


def _get_sample_bam_read_names(bam_path, max_reads=100000):
    """Get read names from a BAM file using samtools view.
    
    Returns set of read names.
    """
    import subprocess
    try:
        result = subprocess.run(
            f"samtools view {bam_path} | cut -f1 | head -n {max_reads}",
            shell=True, capture_output=True, text=True
        )
        return set(result.stdout.strip().split('\n')) - {'read_name'}
    except Exception:
        return set()


def _split_psites_by_sample(psites, tags, fastq_files, bam_dir, tee):
    """Split P-sites by sample based on read name mapping.
    
    Uses the sample BAM files to identify which P-sites belong to which sample.
    Returns dict: tag -> list of psite dicts
    """
    sample_psites = {tag: [] for tag in tags}
    
    for tag in tags:
        bam_path = os.path.join(bam_dir, f"star_{tag}_Aligned.sortedByCoord.out.bam")
        if not os.path.exists(bam_path):
            tee.write(f"    BAM not found for {tag}, trying FASTQ...\n")
            continue
        
        tee.write(f"    Getting read names for {tag}...\n")
        read_names = _get_sample_bam_read_names(bam_path)
        
        if not read_names:
            tee.write(f"    Warning: No reads found for {tag}\n")
            continue
        
        tee.write(f"    Found {len(read_names)} reads for {tag}, assigning P-sites...\n")
        
        # Build a lookup: check if psite read_name is in sample's read names
        # Use prefix matching since P-site file may use truncated names
        read_prefixes = set()
        for rn in read_names:
            # Try to extract the unique part (before any suffix like /1, /2)
            parts = rn.rsplit('/', 1)
            read_prefixes.add(parts[0])
            read_prefixes.add(rn)  # also keep full name
        
        for psite in psites:
            rn = psite['read_name']
            if rn in read_prefixes or rn in read_names:
                sample_psites[tag].append(psite)
    
    return sample_psites


def _plot_psite_analysis(ribo_map_dir, tags, ribo_anno_dir, fastq_files, tee, threads):
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
    
    # Find RIBO Taper P-site output
    psite_file = None
    candidates = [
        "P_sites_all_tracks_exonsccds",
        "P_sites_all_tracks_exonsccds.txt",
        "P_sites_all_tracks",
        "P_sites",
    ]
    for c in candidates:
        if os.path.exists(c):
            psite_file = c
            break
    
    if not psite_file:
        tee.write(f"  Warning: RIBO Taper P-site file not found\n")
        return
    
    tee.write(f"  Reading P-site file: {psite_file}\n")
    psites = _read_psite_tracks(psite_file)
    if not psites:
        tee.write("  Warning: No P-sites found\n")
        return
    
    tee.write(f"  Found {len(psites)} P-sites total\n")
    
    # Find start/stop codon file
    starts_stops_file = None
    ss_candidates = [
        os.path.join(ribo_anno_dir, "start_stops_FAR.bed"),
        os.path.join(ribo_anno_dir, "start_stops.bed"),
        "start_stops_FAR.bed",
        "start_stops.bed",
    ]
    for s in ss_candidates:
        if os.path.exists(s):
            starts_stops_file = s
            break
    
    start_codons = {}
    stop_codons = {}
    if starts_stops_file:
        tee.write(f"  Reading start/stop codons: {starts_stops_file}\n")
        start_codons, stop_codons = _read_start_stops_far(starts_stops_file)
        tee.write(f"  Found start codons on {len(start_codons)} chromosomes\n")
        tee.write(f"  Found stop codons on {len(stop_codons)} chromosomes\n")
    else:
        tee.write("  Warning: start_stops_FAR.bed not found, metagene plots will be empty\n")
    
    # Split P-sites by sample
    tee.write("  Splitting P-sites by sample...\n")
    sample_psites = _split_psites_by_sample(
        psites, tags, fastq_files, ribo_map_dir, tee
    )
    
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
                if 'frame' in p:
                    all_frames.append(p['frame'])
        
        if all_frames:
            frame_counts = np.bincount(all_frames, minlength=3)[:3]
            colors = ['#2ecc71', '#e74c3c', '#3498db']
            bars = axes[0].bar(['Frame 0', 'Frame +1', 'Frame +2'],
                              frame_counts, color=colors, edgecolor='white')
            axes[0].set_ylabel('Count')
            axes[0].set_title('P-site Frame Distribution (All Samples)')
            total = frame_counts.sum()
            axes[0].legend([f'Frame {i}: {c} ({100*c/total:.1f}%)'
                          for i, c in enumerate(frame_counts)],
                         fontsize=8, loc='upper right')
        
        # Per-sample frame heatmap
        valid_samples = [t for t in tags if sample_psites.get(t)]
        if valid_samples:
            heatmap_data = np.zeros((3, len(valid_samples)))
            for j, tag in enumerate(valid_samples):
                frames = [p['frame'] for p in sample_psites[tag] if 'frame' in p]
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
        
        # ── Per-sample pages: frame bar + metagene ──
        for tag in tags:
            assignments = sample_psites.get(tag, [])
            if not assignments:
                continue
            
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
            
            # Frame distribution
            frames = [p['frame'] for p in assignments if 'frame' in p]
            if frames:
                frame_counts = np.bincount(frames, minlength=3)[:3]
                colors = ['#2ecc71', '#e74c3c', '#3498db']
                total = frame_counts.sum()
                bars = ax1.bar(['F0', 'F+1', 'F+2'],
                              frame_counts, color=colors, edgecolor='white')
                ax1.set_title(f'{tag} — Frame Distribution')
                ax1.set_ylabel('Count')
                for i, (bar, count) in enumerate(zip(bars, frame_counts)):
                    pct = 100 * count / total if total > 0 else 0
                    ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                            f'{pct:.1f}%', ha='center', va='bottom', fontsize=9)
            else:
                ax1.text(0.5, 0.5, 'No frame data', ha='center', va='center', transform=ax1.transAxes)
                ax1.set_title(f'{tag} — Frame Distribution')
            
            # Metagene: near start codon
            if start_codons:
                start_positions = []
                for psite in assignments:
                    chrom = psite.get('chrom', '')
                    pos = psite.get('pos', 0)
                    if chrom in start_codons:
                        for gene, (s, e) in start_codons[chrom].items():
                            mid = (s + e) // 2
                            offset = pos - mid
                            if -200 <= offset <= 300:
                                start_positions.append(offset)
                
                if start_positions:
                    bins = np.arange(-200, 301, 5)
                    ax2.hist(start_positions, bins=bins, density=True,
                            color='steelblue', edgecolor='white', alpha=0.7)
                    ax2.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
                    ax2.text(0.5, -0.15, 'Start Codon', transform=ax2.get_xaxis_transform(),
                            ha='center', fontsize=8, color='red')
            ax2.set_title(f'{tag} — P-sites near Start Codon')
            ax2.set_xlabel('Offset from Start Codon (bp)')
            ax2.set_ylabel('Density')
            
            # Metagene: near stop codon
            if stop_codons:
                stop_positions = []
                for psite in assignments:
                    chrom = psite.get('chrom', '')
                    pos = psite.get('pos', 0)
                    if chrom in stop_codons:
                        for gene, (s, e) in stop_codons[chrom].items():
                            mid = (s + e) // 2
                            offset = pos - mid
                            if -300 <= offset <= 100:
                                stop_positions.append(offset)
                
                if stop_positions:
                    bins = np.arange(-300, 101, 5)
                    ax3.hist(stop_positions, bins=bins, density=True,
                            color='darkorange', edgecolor='white', alpha=0.7)
                    ax3.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
                    ax3.text(0.5, -0.15, 'Stop Codon', transform=ax3.get_xaxis_transform(),
                            ha='center', fontsize=8, color='red')
            ax3.set_title(f'{tag} — P-sites near Stop Codon')
            ax3.set_xlabel('Offset from Stop Codon (bp)')
            ax3.set_ylabel('Density')
            
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    
    # Save per-sample P-site files
    tee.write("\n  Saving per-sample P-site files...\n")
    for tag, assignments in sample_psites.items():
        if not assignments:
            continue
        out_file = f"P_sites_{tag}.txt"
        with open(out_file, 'w') as f:
            f.write("read_name\tchrom\tpos\tframe\tlength\n")
            for p in assignments:
                f.write(f"{p.get('read_name','')}\t{p.get('chrom','')}\t{p.get('pos',0)}\t{p.get('frame',0)}\t{p.get('length',0)}\n")
        tee.write(f"    {out_file}: {len(assignments)} P-sites\n")
    
    tee.write(f"\n  P-site analysis complete. Plots saved to: {pdf_path}\n")