"""
TS-CLIP-seq analysis mode.
Similar to CLIP-seq but with 3' PolyC removal before mapping.
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input, _parse_to_dict
from prnaseqtools.functions import (download_sra, unzip_file, rmvc, _tee, run_cmd,
                                     try_use_shared_hisat_index,
                                     save_hisat_index_to_reference,
                                     hisat2_build_index,
                                     hisat2_align_single, hisat2_align_paired,
                                     cleanup_alignment_outputs)


def run(opts):
    """Main entry point for TS-CLIP-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    skip_fastp = opts.get('skip_fastp', False)
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)
    foldchange = opts.get('foldchange', 2.0)
    pvalue = opts.get('pvalue', 0.05)

    # TS-CLIP-seq HISAT2 参数
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
        if not adaptor:
            sys.exit("Please specify the 3' adaptor!")

        tee.write("\nBuilding HISAT2 genome index ...\n")
        hisat_idx = "HisatIdx"
        if not try_use_shared_hisat_index(genome, prefix, hisat_idx, tee):
            gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
            fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
            hisat2_build_index(fasta_path, gff_path, hisat_idx, thread, tee)
            save_hisat_index_to_reference(genome, prefix, hisat_idx, tee)

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nMapping {tag}...\n")

            if ',' not in fpath:
                sra_results, _from_sra = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    unzip_file(sra_results[0], tag)
                    if skip_fastp:
                        tee.write("  [SKIP] fastp (--skip-fastp)\n")
                    else:
                        run_cmd(
                            f"fastp -w {thread} --length_required 20 "
                            f"-i {tag}.fastq -o {tag}_trimmed.fastq "
                            f"--adapter_sequence {adaptor}")
                    rmvc(tag)  # rmvc 输出到 {tag}.fastq
                    hisat2_align_single(
                        tag, f"{tag}.fastq", hisat_idx, thread,
                        extra_args=_HISAT_EXTRA,
                        nov_sj_file=f"{tag}_SJ.out.tab")
                    for fname in (f"{tag}.fastq", f"{tag}_trimmed.fastq"):
                        if os.path.exists(fname):
                            os.unlink(fname)
                else:
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    if skip_fastp:
                        tee.write("  [SKIP] fastp (--skip-fastp)\n")
                    else:
                        run_cmd(
                            f"fastp -w {thread} --length_required 20 "
                            f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                            f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq "
                            f"--adapter_sequence {adaptor} "
                            f"--adapter_sequence_r2 AAAAAAGAAAAAA")
                    rmvc(f"{tag}_R1", f"{tag}_R2")  # 输出到 {tag}_R1.fastq / {tag}_R2.fastq
                    hisat2_align_paired(
                        tag, f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                        hisat_idx, thread,
                        extra_args=_HISAT_EXTRA,
                        nov_sj_file=f"{tag}_SJ.out.tab")
                    for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                                  f"{tag}_R1_trimmed.fastq", f"{tag}_R2_trimmed.fastq"):
                        if os.path.exists(fname):
                            os.unlink(fname)
            else:
                f1, f2 = fpath.split(',')
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")
                if skip_fastp:
                    tee.write("  [SKIP] fastp (--skip-fastp)\n")
                else:
                    run_cmd(
                        f"fastp -w {thread} --length_required 20 "
                        f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                        f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq "
                        f"--adapter_sequence {adaptor} "
                        f"--adapter_sequence_r2 AAAAAAGAAAAAA")
                rmvc(f"{tag}_R1", f"{tag}_R2")
                hisat2_align_paired(
                    tag, f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                    hisat_idx, thread,
                    extra_args=_HISAT_EXTRA,
                    nov_sj_file=f"{tag}_SJ.out.tab")
                for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                              f"{tag}_R1_trimmed.fastq", f"{tag}_R2_trimmed.fastq"):
                    if os.path.exists(fname):
                        os.unlink(fname)

            # HISAT2 直接输出 {tag}.bam
            if os.path.exists(f"{tag}_Log.txt"):
                with open(f"{tag}_Log.txt") as lf:
                    tee.write(lf.read())

            run_cmd(
                f"bamCoverage -b {tag}.bam --skipNAs -bs 5 -p {thread} "
                f"--minMappingQuality 10 "
                f"--normalizeUsing CPM -o {tag}.bw")

            tee.write("\nFinding peaks...\n")
            run_cmd(
                f"clipper -b {tag}.bam -s {genome} --FDR=0.01 --minreads=2 "
                f"--processors={thread} --threshold-method=binomial --min_width=20 "
                f"-o {tag}.fitted_clusters.bed -v")

        cleanup_alignment_outputs(extra_globs=(f"{genome}_chr_all*",
                                                f"{genome}*.bt2"))
        # 清理 HISAT2 index 临时文件
        import glob as _glob
        for f in _glob.glob(f"{hisat_idx}.*.ht2") + _glob.glob(f"{hisat_idx}.htl"):
            if not os.path.islink(f):
                try:
                    os.unlink(f)
                except OSError:
                    pass

        if not mappingonly and len(pars) > 1:
            _ts_ana(pars, prefix, pvalue, foldchange, tee)
    else:
        for pre in tags:
            os.symlink(f"../{pre}.nf", f"{pre}.nf")
            os.symlink(f"../{pre}.bam", f"{pre}.bam")
            os.symlink(f"../{pre}.fitted_clusters.bed", f"{pre}.fitted_clusters.bed")

        if len(pars) > 2:
            _ts_ana(pars, prefix, pvalue, foldchange, tee)

        for fname in globmod.glob("*.bed"):
            os.unlink(fname)
        for fname in globmod.glob("*.nf") + globmod.glob("*.bam*"):
            os.unlink(fname)


def _ts_ana(pars, prefix, pvalue, foldchange, tee):
    """TS-specific peak overlapping (merge distinct on col 4 instead of col 4,6)."""
    tee.write("\nStart to overlap peaks...\n")

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

    run_cmd(
        f"cat {command} | sort -k1,1 -k2,2n | bedtools merge "
        f"-c 4 -o collapse -s -i - > tmp.bed")

    with open("tmp.bed") as fh_peak, open("ref.bed", 'w') as fh_ref:
        for line in fh_peak:
            cols = line.strip().split('\t')
            if len(cols) < 5:
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
                    f"{cols[0]}\t{cols[1]}\t{cols[2]}\t{output_gene}\t{total_sum}\t{cols[4]}\n"
                )

    tee.write("\nCounting reads in each sample...\n")
    for tag in all_tags:
        run_cmd(
            f"bedtools intersect -c -s -a ref.bed -b {tag}.bam | "
            f"awk '{{print $1\"_\"$2\"_\"$3\"_\"$6\"_\"$4\"\\t\"$7}}' > {tag}.txt")

    tee.write("\nFinding differential peaks...\n")
    run_cmd(
        f"Rscript --vanilla {prefix}/scripts/CLIP.R {pvalue} {foldchange} "
        + ' '.join(pars))

    for fname in ("tmp.bed", "ref.bed"):
        if os.path.exists(fname):
            os.unlink(fname)
