"""
CLIP-seq analysis mode.
HISAT2 alignment → CLIPper peak finding → differential peak analysis.
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path
from collections import defaultdict

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input, _parse_to_dict
from prnaseqtools.functions import (download_sra, unzip_file, _tee, run_cmd,
                                     gzip_fastq, preserve_raw_fastq,
                                     try_use_shared_hisat_index,
                                     save_hisat_index_to_reference,
                                     hisat2_build_index,
                                     hisat2_align_single, hisat2_align_paired,
                                     cleanup_alignment_outputs)


def run(opts):
    """Main entry point for CLIP-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)
    foldchange = opts.get('foldchange', 2.0)
    pvalue = opts.get('pvalue', 0.05)
    skip_fastp = opts.get('skip_fastp', False)

    # CLIP-seq HISAT2 参数：严格唯一比对，小内含子 (CLIP 是 RBP 结合位点)
    _HISAT_EXTRA = "-k 1 --mp 3,3 --np 1 --score-min L,0,-0.3 --max-intronlen 5000"

    control_dict = _parse_to_dict(opts.get('control', ''))
    tags, files, pars = parse_input(control_dict)

    if opts.get('treatment'):
        for t in (opts['treatment'] if isinstance(opts['treatment'], list) else [opts['treatment']]):
            treatment_dict = _parse_to_dict(t)
            t_tags, t_files, t_pars = parse_input(treatment_dict)
            tags.extend(t_tags)
            files.extend(t_files)
            pars.extend(t_pars)

    if not nomapping:
        if not adaptor and not skip_fastp:
            sys.exit("Please specify the 3' adaptor (or use --skip-fastp).")

        tee.write("\nBuilding HISAT2 genome index ...\n")
        hisat_idx = "HisatIdx"
        if not try_use_shared_hisat_index(genome, prefix, hisat_idx, tee):
            # HISAT2 index 是前缀文件组 (HisatIdx.1.ht2 ...)，不是目录
            gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
            fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
            hisat2_build_index(fasta_path, gff_path, hisat_idx, thread, tee)
            save_hisat_index_to_reference(genome, prefix, hisat_idx, tee)

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nMapping {tag}...\n")

            if ',' not in fpath:
                sra_results, from_sra = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    unzip_file(sra_results[0], tag)
                    preserve_raw_fastq(from_sra, tag, paired=False)
                    if skip_fastp:
                        tee.write("  [SKIP] fastp (--skip-fastp)\n")
                        fq_in = f"{tag}.fastq"
                    else:
                        run_cmd(
                            f"fastp -w {thread} --length_required 18 "
                            f"-i {tag}.fastq -o {tag}_trimmed.fastq "
                            f"--adapter_sequence {adaptor}")
                        fq_in = f"{tag}_trimmed.fastq"
                    hisat2_align_single(
                        tag, fq_in, hisat_idx, thread,
                        extra_args=_HISAT_EXTRA,
                        nov_sj_file=f"{tag}_SJ.out.tab")
                    # 提取 unmapped reads 到 fastq（兼容旧版 cleanup）
                    if os.path.exists(f"{tag}.bam"):
                        run_cmd(
                            f"samtools view -b -f 4 {tag}.bam | "
                            f"samtools fastq - > {tag}.unmapped.fastq 2>/dev/null || true")
                    if not skip_fastp and os.path.exists(f"{tag}_trimmed.fastq"):
                        os.unlink(f"{tag}_trimmed.fastq")
                    if from_sra:
                        gzip_fastq(f"{tag}.raw.fastq")
                else:
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    preserve_raw_fastq(from_sra, tag, paired=True)
                    if skip_fastp:
                        tee.write("  [SKIP] fastp (--skip-fastp)\n")
                        r1_in, r2_in = f"{tag}_R1.fastq", f"{tag}_R2.fastq"
                    else:
                        run_cmd(
                            f"fastp -w {thread} --length_required 18 "
                            f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                            f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq "
                            f"--adapter_sequence {adaptor} "
                            f"--adapter_sequence_r2 AAAAAAGAAAAAA")
                        r1_in, r2_in = f"{tag}_R1_trimmed.fastq", f"{tag}_R2_trimmed.fastq"
                    hisat2_align_paired(
                        tag, r1_in, r2_in,
                        hisat_idx, thread,
                        extra_args=_HISAT_EXTRA,
                        nov_sj_file=f"{tag}_SJ.out.tab")
                    if os.path.exists(f"{tag}.bam"):
                        run_cmd(
                            f"samtools view -b -f 4 {tag}.bam | "
                            f"samtools fastq - > {tag}.unmapped.fastq 2>/dev/null || true")
                    if not skip_fastp:
                        for _f in (f"{tag}_R1_trimmed.fastq", f"{tag}_R2_trimmed.fastq"):
                            if os.path.exists(_f):
                                os.unlink(_f)
                    if from_sra:
                        for _f in (f"{tag}_R1.raw.fastq", f"{tag}_R2.raw.fastq"):
                            gzip_fastq(_f)
            else:
                f1, f2 = fpath.split(',')
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")
                if skip_fastp:
                    tee.write("  [SKIP] fastp (--skip-fastp)\n")
                    r1_in, r2_in = f"{tag}_R1.fastq", f"{tag}_R2.fastq"
                else:
                    run_cmd(
                        f"fastp -w {thread} --length_required 18 "
                        f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                        f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq "
                        f"--adapter_sequence {adaptor} "
                        f"--adapter_sequence_r2 AAAAAAGAAAAAA")
                    r1_in, r2_in = f"{tag}_R1_trimmed.fastq", f"{tag}_R2_trimmed.fastq"
                hisat2_align_paired(
                    tag, r1_in, r2_in,
                    hisat_idx, thread,
                    extra_args=_HISAT_EXTRA,
                    nov_sj_file=f"{tag}_SJ.out.tab")
                if os.path.exists(f"{tag}.bam"):
                    run_cmd(
                        f"samtools view -b -f 4 {tag}.bam | "
                        f"samtools fastq - > {tag}.unmapped.fastq 2>/dev/null || true")
                if not skip_fastp:
                    for _f in (f"{tag}_R1_trimmed.fastq", f"{tag}_R2_trimmed.fastq"):
                        if os.path.exists(_f):
                            os.unlink(_f)

            # HISAT2 直接输出 {tag}.bam，不需要 rename
            # HISAT2 日志
            if os.path.exists(f"{tag}_Log.txt"):
                with open(f"{tag}_Log.txt") as lf:
                    tee.write(lf.read())

            run_cmd(
                f"bamCoverage -b {tag}.bam --skipNAs -bs 5 -p {thread} "
                f"--minMappingQuality 10 "
                f"--normalizeUsing CPM -o {tag}.bw")

            tee.write("\nAlignment completed!\nFinding peaks...\n")
            run_cmd(
                f"clipper -b {tag}.bam -s {genome} --FDR=0.01 --minreads=2 "
                f"--processors={thread} --threshold-method=binomial --min_width=20 "
                f"-o {tag}.fitted_clusters.bed -v")

        cleanup_alignment_outputs(extra_globs=(f"{genome}_chr_all*",
                                                f"{genome}*.bt2",
                                                "unmapped*.fastq"))
        # 清理 HISAT2 index 临时文件（如果不是 symlink）
        import glob as _glob
        for f in _glob.glob(f"{hisat_idx}.*.ht2") + _glob.glob(f"{hisat_idx}.htl"):
            if not os.path.islink(f):
                try:
                    os.unlink(f)
                except OSError:
                    pass

        if not mappingonly and len(pars) > 1:
            _differential_analysis(pars, prefix, pvalue, foldchange, tee)
    else:
        for pre in tags:
            os.symlink(f"../{pre}.nf", f"{pre}.nf")
            os.symlink(f"../{pre}.bam", f"{pre}.bam")
            os.symlink(f"../{pre}.fitted_clusters.bed", f"{pre}.fitted_clusters.bed")

        if len(pars) > 2:
            _differential_analysis(pars, prefix, pvalue, foldchange, tee)

        for fname in globmod.glob("*.bed"):
            os.unlink(fname)
        for fname in globmod.glob("*.nf"):
            os.unlink(fname)
        for fname in globmod.glob("*.bam"):
            if os.path.exists(fname):
                os.unlink(fname)
        for fname in globmod.glob("*.bam.bai"):
            if os.path.exists(fname):
                os.unlink(fname)


def _differential_analysis(pars, prefix, pvalue, foldchange, tee):
    """Overlap peaks and find differential peaks."""
    tee.write("\nOverlapping peaks...\n")

    p_list = list(pars)
    all_tags = []
    command = ""

    while len(p_list) >= 2:
        group_name = p_list.pop(0)
        n_reps = int(p_list.pop(0))
        for s in range(1, n_reps + 1):
            tag = f"{group_name}_{s}"
            all_tags.append(tag)
            command += f"{tag}.fitted_clusters.bed "

    par_str = ' '.join(pars)

    # Merge peaks
    run_cmd(
        f"cat {command} | sort -k1,1 -k2,2n | bedtools merge "
        f"-c 4,6 -o collapse,distinct -s -i - > tmp.bed")

    with open("tmp.bed") as fh_peak, open("ref.bed", 'w') as fh_ref:
        for line in fh_peak:
            cols = line.strip().split('\t')
            if len(cols) < 6:
                continue

            peak_names = cols[3].split(',') if ',' in cols[3] else [cols[3]]
            pgenes = {}
            total_sum = 0

            for peak in peak_names:
                sample_parts = peak.split('_')
                if sample_parts:
                    pgenes[sample_parts[0]] = ''
                    total_sum += int(sample_parts[2]) if len(sample_parts) > 2 else 0

            output_gene = '_'.join(sorted(pgenes.keys()))
            if total_sum >= 10:
                fh_ref.write(
                    f"{cols[0]}\t{cols[1]}\t{cols[2]}\t{output_gene}\t"
                    f"{cols[3]}_{total_sum}\t{cols[4]}\n"
                )

    # Count reads per sample
    tee.write("Counting reads in each sample...\n")
    for tag in all_tags:
        run_cmd(
            f"bedtools intersect -c -s -a ref.bed -b {tag}.bam | "
            f"awk '{{print $1\"_\"$2\"_\"$3\"_\"$6\"_\"$4\"\\t\"$7}}' > {tag}.txt")

    tee.write("Finding differential peaks...\n")
    run_cmd(
        f"Rscript --vanilla {prefix}/scripts/CLIP.R {pvalue} {foldchange} {par_str}")

    for fname in ("tmp.bed", "ref.bed"):
        if os.path.exists(fname):
            os.unlink(fname)
