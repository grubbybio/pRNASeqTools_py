"""
mRNA-seq analysis mode.
Mirrors the Perl mode::mrna module.
HISAT2 alignment → featureCounts → DESeq2 DE analysis.
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import (parse_input, _parse_to_dict,
                                        _resolve_path)
from prnaseqtools.functions import (
    download_sra, unzip_file, _tee, run_cmd, gzip_fastq, preserve_raw_fastq,
    try_use_shared_hisat_index, save_hisat_index_to_reference,
    cleanup_alignment_outputs, hisat2_build_index,
    hisat2_align_single, hisat2_align_paired,
)


def run(opts):
    """Main entry point for mRNA-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    # prefix must point to the project root so that prefix/reference/{genome}_Genome
    # resolves to ~/software/pRNASeqTools_py/reference/ (not prnaseqtools/reference/)
    prefix = opts.get('prefix', '/Users/cjyou/software/pRNASeqTools_py')
    foldchange = opts.get('foldchange', 2.0)
    pvalue = opts.get('pvalue', 0.01)
    fdr = opts.get('fdr', 1.0)
    run_mode = opts.get('run_mode', 'whole')
    control = opts.get('control', '')
    treatment = opts.get('treatment')
    norm = opts.get('deseq2_norm', 'DESeq2')
    mask = opts.get('mask')
    total = opts.get('total', False)
    seq_strategy = opts.get('seq_strategy')
    genome_size = opts.get('genome_size', 10)
    skip_fastp = opts.get('skip_fastp', False)

    mapping = run_mode in ('whole', 'mapping-only')
    do_count = run_mode != 'count-table'
    do_de = run_mode != 'mapping-only'

    # Parse inputs
    control_dict = _parse_to_dict(control)
    tags, files, pars = parse_input(control_dict)

    if treatment:
        for t in (treatment if isinstance(treatment, list) else [treatment]):
            treatment_dict = _parse_to_dict(t)
            t_tags, t_files, t_pars = parse_input(treatment_dict)
            tags.extend(t_tags)
            files.extend(t_files)
            pars.extend(t_pars)

    par_str = ' '.join(pars)

    if mapping:
        # Mask setup
        if mask:
            mask_path = _resolve_path(mask)
            os.symlink(mask_path, "mask.fa")
            run_cmd("bowtie-build -q mask.fa mask")

        tee.write("\nBuilding HISAT2 genome index ...\n")

        # Define reference paths first — they are also needed by the gffread call
        # below even when try_use_shared_hisat_index returns True.
        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        hisat_idx = "HisatIdx"  # HISAT2 索引前缀（不是目录）

        if not try_use_shared_hisat_index(genome, prefix, hisat_idx, tee):
            tee.write("\nBuilding HISAT2 index from FASTA + known splice sites ...\n")
            # 用 gffread 先把 GFF → GTF（extract_splice_sites.py 偏好 GTF）
            run_cmd(
                f"gffread -T -o _hisat_genes.gtf -g {fasta_path} {gff_path}")
            hisat2_build_index(fasta_path, "_hisat_genes.gtf", hisat_idx, thread, tee)
            save_hisat_index_to_reference(genome, prefix, hisat_idx, tee)
            os.unlink("_hisat_genes.gtf")

        if total:
            run_cmd(
                f"gffread -T -o {genome}_genes.gtf -g {fasta_path} {gff_path}")
        else:
            run_cmd(
                f"gffread -T -C -o {genome}_genes.gtf -g {fasta_path} {gff_path}")

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nProcessing {tag}...\n")

            if ',' not in fpath:
                # Single-end or paired from SRA
                sra_results, from_sra = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    if seq_strategy is None:
                        seq_strategy = 'single'
                    unzip_file(sra_results[0], tag)
                    preserve_raw_fastq(from_sra, tag, paired=False)

                    # fastp trimming (skippable via --skip-fastp)
                    if skip_fastp:
                        tee.write("  [SKIP] fastp (--skip-fastp)\n")
                    else:
                        tee.write("  Trimming (fastp)...\n")
                        fastp_adaptor = f"--adapter_sequence {adaptor}" if adaptor else ""
                        run_cmd(
                            f"fastp -w {thread} --length_required 20 "
                            f"-i {tag}.fastq -o {tag}_trimmed.fastq "
                            f"-h {tag}.html -j {tag}.json {fastp_adaptor}")
                        os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")
                        # clean up fastp reports
                        for rep in (f"{tag}.html", f"{tag}.json"):
                            if os.path.exists(rep):
                                os.unlink(rep)

                    if mask:
                        run_cmd(
                            f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                            f"{tag}.fastq {tag}.mask.out")
                        os.rename("tmp.fastq", f"{tag}.fastq")
                        if os.path.exists(f"{tag}.mask.out"):
                            os.unlink(f"{tag}.mask.out")

                    hisat2_align_single(tag, f"{tag}.fastq", hisat_idx, thread,
                        extra_args="-k 50 --mp 3,3 --score-min L,0,-0.3",
                        nov_sj_file=f"{tag}_SJ.out.tab")
                    if from_sra:
                        gzip_fastq(f"{tag}.raw.fastq")
                else:
                    # Paired from SRA
                    if seq_strategy is None:
                        seq_strategy = 'paired'
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    preserve_raw_fastq(from_sra, tag, paired=True)

                    # fastp trimming (skippable via --skip-fastp)
                    if skip_fastp:
                        tee.write("  [SKIP] fastp (--skip-fastp)\n")
                    else:
                        tee.write("  Trimming (fastp)...\n")
                        fastp_adaptor = f"--adapter_sequence {adaptor} --adapter_sequence_r2 {adaptor}" if adaptor else ""
                        run_cmd(
                            f"fastp -w {thread} --length_required 20 "
                            f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                            f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq "
                            f"-h {tag}.html -j {tag}.json {fastp_adaptor}")
                        os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                        os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                        for rep in (f"{tag}.html", f"{tag}.json"):
                            if os.path.exists(rep):
                                os.unlink(rep)

                    if mask:
                        run_cmd(
                            f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                            f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq {tag}.mask.out")
                        if os.path.exists("tmp_1.fastq"):
                            os.rename("tmp_1.fastq", f"{tag}_R1.fastq")
                        if os.path.exists("tmp_2.fastq"):
                            os.rename("tmp_2.fastq", f"{tag}_R2.fastq")
                        if os.path.exists(f"{tag}.mask.out"):
                            os.unlink(f"{tag}.mask.out")

                    hisat2_align_paired(tag, f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                        hisat_idx, thread,
                        extra_args="-k 50 --mp 3,3 --score-min L,0,-0.3",
                        nov_sj_file=f"{tag}_SJ.out.tab")
                    if from_sra:
                        for _f in (f"{tag}_R1.raw.fastq", f"{tag}_R2.raw.fastq"):
                            gzip_fastq(_f)
            else:
                # Explicit paired-end (local files; no SRA preservation)
                f1, f2 = fpath.split(',')
                if seq_strategy is None:
                    seq_strategy = 'paired'
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")

                # fastp trimming (skippable via --skip-fastp)
                if skip_fastp:
                    tee.write("  [SKIP] fastp (--skip-fastp)\n")
                else:
                    tee.write("  Trimming (fastp)...\n")
                    fastp_adaptor = f"--adapter_sequence {adaptor} --adapter_sequence_r2 {adaptor}" if adaptor else ""
                    run_cmd(
                        f"fastp -w {thread} --length_required 20 "
                        f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                        f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq "
                        f"-h {tag}.html -j {tag}.json {fastp_adaptor}")
                    os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                    os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                    for rep in (f"{tag}.html", f"{tag}.json"):
                        if os.path.exists(rep):
                            os.unlink(rep)

                if mask:
                    run_cmd(
                        f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                        f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq {tag}.mask.out")
                    if os.path.exists("tmp_1.fastq"):
                        os.rename("tmp_1.fastq", f"{tag}_R1.fastq")
                    if os.path.exists("tmp_2.fastq"):
                        os.rename("tmp_2.fastq", f"{tag}_R2.fastq")
                    if os.path.exists(f"{tag}.mask.out"):
                        os.unlink(f"{tag}.mask.out")

                tee.write("\nMapping...\n")
                hisat2_align_paired(tag, f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                    hisat_idx, thread,
                    extra_args="-k 50 --mp 3,3 --score-min L,0,-0.3",
                    nov_sj_file=f"{tag}_SJ.out.tab")

            # HISAT2 直接输出 {tag}.bam，不再需要重命名 STAR 的 _Aligned.sortedByCoord.out.bam
            if not os.path.exists(f"{tag}.bam"):
                _aligned_bam = f"{tag}_Aligned.sortedByCoord.out.bam"  # 兼容旧残留
                if os.path.exists(_aligned_bam):
                    os.rename(_aligned_bam, f"{tag}.bam")

            # Show HISAT2 mapping stats
            if os.path.exists(f"{tag}_Log.txt"):
                with open(f"{tag}_Log.txt") as lf:
                    tee.write(lf.read())

            _count_bam(tag, thread, seq_strategy, prefix, genome)

        # Cleanup
        for fname in (f"{genome}_genes.gtf",):
            if os.path.exists(fname):
                os.unlink(fname)
        # Per-tag HISAT2/STAR outputs, fastp intermediates, and fastp reports
        cleanup_alignment_outputs(extra_globs=(f"{genome}_chr_all*", "total.count*", "_hisat_genes.gtf"))
        if mask:
            for mf in globmod.glob("mask*"):
                os.unlink(mf)
        # NOTE: HISAT2 index 文件 (HisatIdx.*.ht2) 不在这里删——可能是 symlink 到 reference/{genome}_HisatIdx
        # STAR 旧版是目录 "Genome"，HISAT2 是前缀文件组，都在 cleanup_alignment_outputs 里不碰

        if do_de:
            tee.write(f"\nFinding DEG...\nFold Change\t{foldchange}\tP Value\t{pvalue}\tFDR\t{fdr}\n")
            run_cmd(
                f"Rscript --vanilla {prefix}/scripts/DEG.R {norm} {pvalue} {fdr} "
                f"{foldchange} {prefix} {genome} {par_str}")
            # ── GO enrichment for each comparison ────────────────────────────
            _run_go_enrichment(prefix, genome, tee)

    else:
        # No mapping: symlink existing files or process BAMs
        for pre in tags:
            if do_count:
                os.symlink(f"../{pre}.bam", f"{pre}.bam")
            else:
                os.symlink(f"../{pre}.txt", f"{pre}.txt")

        if do_count:
            if total:
                gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
                fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
                run_cmd(
                    f"gffread -T -o {genome}_genes.gtf -g {fasta_path} {gff_path}")
            else:
                gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
                fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
                run_cmd(
                    f"gffread -T -C -o {genome}_genes.gtf -g {fasta_path} {gff_path}")
            for pre in tags:
                _count_bam(pre, thread, seq_strategy, prefix, genome)
            if os.path.exists(f"{genome}_genes.gtf"):
                os.unlink(f"{genome}_genes.gtf")

        tee.write(f"\nFinding DEG...\nFold Change\t{foldchange}\tP Value\t{pvalue}\tFDR\t{fdr}\n")
        run_cmd(
            f"Rscript --vanilla {prefix}/scripts/DEG.R {norm} {pvalue} {fdr} "
            f"{foldchange} {prefix} {genome} {par_str}")
        # ── GO enrichment for each comparison ────────────────────────────
        _run_go_enrichment(prefix, genome, tee)

        # Clean up symlinks
        ext = "bam" if do_count else "txt"
        for pre in tags:
            fname = f"{pre}.{ext}"
            if os.path.exists(fname):
                os.unlink(fname)


def _count_bam(tag, thread, seq_strategy, prefix, genome):
    """Count reads in BAM using featureCounts."""
    run_cmd(f"samtools index {tag}.bam")
    run_cmd(
        f"bamCoverage -b {tag}.bam --skipNAs -bs 5 -p {thread} "
        f"--minMappingQuality 10 --normalizeUsing CPM -o {tag}.bw")

    tee = _tee()
    tee.write("\nStart counting...\n")

    fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
    gtf_file = f"{genome}_genes.gtf"

    if seq_strategy == 'single':
        run_cmd(
            f"featureCounts -T {thread} -O -G {fasta_path} -s 0 "
            f"-a {gtf_file} -o total.count {tag}.bam")
    elif seq_strategy == 'paired':
        run_cmd(
            f"featureCounts -T {thread} -p --countReadPairs -BCO -G {fasta_path} "
            f"-s 0 -a {gtf_file} -o total.count {tag}.bam")
    else:
        # Default to single
        run_cmd(
            f"featureCounts -T {thread} -O -G {fasta_path} -s 0 "
            f"-a {gtf_file} -o total.count {tag}.bam")

    count_data = {}
    count_sum = 0
    with open("total.count") as fh:
        fh.readline()  # skip header
        fh.readline()  # skip header
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) >= 7:
                count_data[cols[0]] = {
                    'exon': int(cols[6]),
                    'length': int(cols[5])
                }
                count_sum += int(cols[6])

    tee.write(f"\nRead Count: {count_sum}\n")

    with open(f"{tag}.txt", 'w') as fh:
        fh.write("Gene\tCount\tLength\n")
        for name in sorted(count_data.keys()):
            fh.write(f"{name}\t{count_data[name]['exon']}\t{count_data[name]['length']}\n")


def _run_go_enrichment(prefix, genome, tee):
    """Run GO enrichment for every comparison found in the working dir.

    Finds all ``*.total.upregulated.csv`` files produced by DEG.R,
    derives the comparison name, and invokes scripts/GO_enrichment.R.
    """
    import glob as _glob
    up_files = _glob.glob("*.total.upregulated.csv")
    if not up_files:
        tee.write("\nNo DEG comparison files found — skipping GO enrichment.\n")
        return
    tee.write("\n" + "=" * 60 + "\n")
    tee.write("GO Enrichment Analysis\n")
    tee.write("=" * 60 + "\n")
    for uf in sorted(up_files):
        # "{comparison}.total.upregulated.csv" → "{comparison}"
        comparison = uf.replace(".total.upregulated.csv", "")
        tee.write(f"\n  Processing: {comparison}\n")
        run_cmd(
            f"Rscript --vanilla {prefix}/scripts/GO_enrichment.R "
            f"{genome} {prefix} {comparison}")
    tee.write("\nGO enrichment done.\n")
