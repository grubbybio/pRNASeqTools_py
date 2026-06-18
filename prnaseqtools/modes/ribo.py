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
import glob as globmod
import shutil
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input
from prnaseqtools.functions import download_sra, unzip_file, _tee, run_cmd


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
    ribo_len = opts.get('ribo_len', '24,25,26,27,28')
    cutoffs = opts.get('cutoffs', '8,9,10,11,12')
    tpm_threshold = opts.get('tpm_threshold', 0)
    ribotaper_path = opts.get('ribotaper')
    if not ribotaper_path:
        ribotaper_path = _find_ribotaper(prefix)

    # ── Parse RNA-seq inputs ─────────────────────────────────────────
    rna_ctrl_dict = _parse_to_dict(opts.get('rna_control', ''))
    rna_tags, rna_files, _ = parse_input(rna_ctrl_dict)

    rna_trt = opts.get('rna_treatment')
    if rna_trt:
        trt_dict = _parse_to_dict(rna_trt)
        rt_tags, rt_files, _ = parse_input(trt_dict)
        rna_tags.extend(rt_tags)
        rna_files.extend(rt_files)

    all_rna_tags = list(rna_tags)

    # ── Parse Ribo-seq inputs ────────────────────────────────────────
    ribo_ctrl_dict = _parse_to_dict(opts.get('ribo_control', ''))
    ribo_tags, ribo_files, _ = parse_input(ribo_ctrl_dict)

    ribo_trt = opts.get('ribo_treatment')
    if ribo_trt:
        rbt_dict = _parse_to_dict(ribo_trt)
        rbt_tags, rbt_files, _ = parse_input(rbt_dict)
        ribo_tags.extend(rbt_tags)
        ribo_files.extend(rbt_files)

    all_ribo_tags = list(ribo_tags)

    # ── Reference paths ──────────────────────────────────────────────
    gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
    fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")

    for p in (gff_path, fasta_path):
        if not os.path.exists(p):
            sys.exit(f"Reference file not found: {p}")

    # ── Utility: find bedtools ───────────────────────────────────────
    bedtools_bin = shutil.which("bedtools") or "bedtools"

    # ==================================================================
    # STEP 1 — Build Bowtie2 contamination index
    # ==================================================================
    tee.write("\n" + "=" * 60 + "\n")
    tee.write("STEP 1: Building Bowtie2 contamination index\n")
    tee.write("=" * 60 + "\n")

    contam_index = "Contam"
    run_cmd(
        f"bowtie2-build {contam} {contam_index}")
    tee.write("  Contamination index built.\n")

    # ==================================================================
    # STEP 2 — Preprocess Ribo-seq reads (decontamination)
    # ==================================================================
    tee.write("\n" + "=" * 60 + "\n")
    tee.write("STEP 2: Preprocessing Ribo-seq reads\n")
    tee.write("=" * 60 + "\n")

    for i in range(len(ribo_tags)):
        tag = ribo_tags[i]
        fpath = ribo_files[i]

        tee.write(f"\n  --- Ribo-seq sample: {tag} ---\n")

        # Download / unzip
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
    tee.write("\n" + "=" * 60 + "\n")
    tee.write("STEP 3: RNA-seq transcriptome assembly (STAR + StringTie)\n")
    tee.write("=" * 60 + "\n")

    # Build STAR genome index for RNA-seq
    star_rna_index = "STAR_RNA_index"
    if os.path.exists(star_rna_index):
        run_cmd(f"rm -rf {star_rna_index}", shell=True)
    os.makedirs(star_rna_index, exist_ok=True)

    run_cmd(
        f"STAR --runThreadN {thread} --runMode genomeGenerate "
        f"--genomeDir {star_rna_index} --genomeFastaFiles {fasta_path} "
        f"--sjdbGTFfile {gff_path} --sjdbOverhang 99 "
        f"--limitGenomeGenerateRAM 64000000000")

    # STAR 1st pass for each RNA-seq sample
    star_sj_dir = "STAR_SJ"
    os.makedirs(star_sj_dir, exist_ok=True)

    for tag, fpath in zip(rna_tags, rna_files):
        tee.write(f"\n  --- RNA-seq STAR 1st pass: {tag} ---\n")

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
            f"stringtie --rf -p {thread} -G {gff_path} "
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
        f"stringtie --merge -p {thread} -T 0.05 -G {gff_path} "
        f"-o {merged_gtf} {merge_list_file}")

    # gffcompare
    tee.write("  Running gffcompare...\n")
    gffcomp_prefix = os.path.join(assembled_gtf_dir, genome)
    run_cmd(
        f"gffcompare -V -r {gff_path} -o {gffcomp_prefix} {merged_gtf}")

    # ==================================================================
    # STEP 4 — R: Filter novel transcripts + add gene_biotype
    # ==================================================================
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
        f"{annotated_gtf} {gff_path} {updated_gtf}")
    tee.write(f"  Updated GTF written to: {updated_gtf}\n")

    # ==================================================================
    # STEP 5 — RSEM quantification + filter expressed isoforms
    # ==================================================================
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
    tee.write("\n" + "=" * 60 + "\n")
    tee.write("STEP 7: Creating RIBO Taper annotation files\n")
    tee.write("=" * 60 + "\n")

    ribo_anno_dir = "RiboTaper_annotation"
    os.makedirs(ribo_anno_dir, exist_ok=True)

    create_anno_script = os.path.join(
        ribotaper_path, "scripts", "create_annotations_files.bash"
    )
    if not os.path.exists(create_anno_script):
        sys.exit(f"RIBO Taper script not found: {create_anno_script}")

    # Ensure genome fasta is indexed (fai)
    fai_file = f"{fasta_path}.fai"
    if not os.path.exists(fai_file):
        run_cmd(f"samtools faidx {fasta_path}")

    run_cmd(
        f"bash {create_anno_script} {expressed_gtf} {fasta_path} "
        f"false false {ribo_anno_dir} {bedtools_bin} "
        f"{ribotaper_path}/scripts/")
    tee.write("  RIBO Taper annotation created.\n")

    # ==================================================================
    # STEP 8 — Merge BAMs + Run RIBO Taper
    # ==================================================================
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

    # Merge Ribo-seq BAMs
    tee.write("  Merging Ribo-seq BAMs...\n")
    ribo_bam_list = [
        f"{ribo_map_dir}/star_{tag}_Aligned.sortedByCoord.out.bam"
        for tag in all_ribo_tags
    ]
    ribo_bam_list = [b for b in ribo_bam_list if os.path.exists(b)]
    ribo_bams_arg = ' '.join(ribo_bam_list)

    ribo_merged = "Ribo_merged.bam"
    run_cmd(
        f"samtools merge -f {ribo_merged} {ribo_bams_arg}")
    run_cmd(
        f"samtools sort -o Ribo_merged_sorted.bam {ribo_merged}")
    os.rename("Ribo_merged_sorted.bam", ribo_merged)
    run_cmd(f"samtools index {ribo_merged}")

    # Run RIBO Taper
    tee.write("\n  Running RIBO Taper ORF detection...\n")
    ribotaper_script = os.path.join(ribotaper_path, "scripts", "Ribotaper.sh")
    if not os.path.exists(ribotaper_script):
        sys.exit(f"RIBO Taper script not found: {ribotaper_script}")

    run_cmd(
        f"bash {ribotaper_script} {ribo_merged} {rna_merged} "
        f"{ribo_anno_dir} {ribo_len} {cutoffs} "
        f"{ribotaper_path}/scripts/ {bedtools_bin} {thread}")

    tee.write("\n  RIBO Taper analysis complete.\n")

    # ==================================================================
    # STEP 9 — Final output summary
    # ==================================================================
    tee.write("\n" + "=" * 60 + "\n")
    tee.write("STEP 9: Pipeline complete\n")
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

def _parse_to_dict(arg_str):
    """Parse 'name=value' string to dict."""
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _find_ribotaper(prefix):
    """Locate RIBO Taper installation directory."""
    candidates = [
        os.path.join(prefix, "RiboTaper_v1.3"),
        os.path.expanduser("~/RiboTaper_v1.3"),
        os.path.expanduser("~/Software/RiboTaper_v1.3"),
        os.path.expanduser("~/software/RiboTaper_v1.3"),
    ]
    for cand in candidates:
        script = os.path.join(cand, "scripts", "Ribotaper.sh")
        if os.path.exists(script):
            return cand

    sys.exit(
        "Cannot find RIBO Taper installation.\n"
        "Please specify --ribotaper PATH or install RIBO Taper from:\n"
        "  https://github.com/hsinyenwu/RiboTaper"
    )
