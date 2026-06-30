"""
ChIP-seq analysis mode.
bowtie2 alignment → Genrich / MACS3 peak calling → Peak QC (FRiP, avg length, etc.)
Group comparison (bdgdiff) is implemented in the tf module.
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input
from prnaseqtools.functions import download_sra, unzip_file, _tee, run_cmd


def run(opts):
    """Main entry point for ChIP-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    peak_caller = opts.get('peak_caller', 'genrich')
    genome_size = opts.get('genome_size')
    if peak_caller == 'macs3' and not genome_size:
        sys.exit("--genome-size is required for MACS3 (e.g. 1.35e8 for ath)")
    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    auc = opts.get('auc', 20)
    qvalue = opts.get('qvalue', 1.0)
    pvalue = opts.get('pvalue', 0.01)
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)
    seq_strategy = opts.get('seq_strategy', 'paired')

    tags = []
    files = []

    control_tags = []
    input_opt = opts.get('control')
    if isinstance(input_opt, list):
        input_opt = input_opt[0] if input_opt else ''
    genrich_input = ""
    if input_opt:
        input_dict = _parse_to_dict(input_opt)
        i_tags, i_files, i_pars = parse_input(input_dict)
        control_tags.extend(i_tags)
        tags.extend(i_tags)
        files.extend(i_files)
        genrich_input = "-c " + ','.join(f"{t}.sorted.name.bam" for t in i_tags)

    # Parse IP (treatment) — 支持多组: -p group1=file1 -p group2=file2
    treatment_groups = []
    ip_opt = opts.get('treatment', '')
    if isinstance(ip_opt, str):
        ip_opt = [ip_opt]
    if ip_opt:
        for treat_str in ip_opt:
            t_dict = _parse_to_dict(treat_str)
            i_tags, i_files, i_pars = parse_input(t_dict)
            treatment_groups.append({
                'ip_tags': list(i_tags),
                'ip_files': list(i_files),
                'ip_tag': i_tags[0] if i_tags else 'ip',
            })
            tags.extend(i_tags)
            files.extend(i_files)

    if not treatment_groups:
        sys.exit("At least one --treatment (-p) is required.")

    if not nomapping:
        tee.write("\nBuilding index...\n")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        os.symlink(fasta_path, f"{genome}_chr_all.fasta")
        run_cmd(
            f"bowtie2-build --threads {thread} -q {genome}_chr_all.fasta {genome}_chr_all")

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nMapping {tag}...\n")

            if ',' not in fpath:
                sra_results = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    seq_strategy = 'single'
                    unzip_file(sra_results[0], tag)
                    if adaptor:
                        tee.write("\nTrimming...\n")
                        run_cmd(
                            f"cutadapt -j {thread} -m 50 --trim-n -a {adaptor} "
                            f"-o {tag}_trimmed.fastq {tag}.fastq")
                        os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")
                    run_cmd(
                        f"bowtie2 -p {thread} -x {genome}_chr_all -U {tag}.fastq "
                        f"-S {tag}.sam")
                    if os.path.exists(f"{tag}.fastq"):
                        os.unlink(f"{tag}.fastq")
                else:
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    if adaptor:
                        run_cmd(
                            f"cutadapt -j {thread} -m 50 --trim-n -a {adaptor} -A {adaptor} "
                            f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                            f"{tag}_R1.fastq {tag}_R2.fastq")
                        os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                        os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                    run_cmd(
                        f"bowtie2 -p {thread} -x {genome}_chr_all "
                        f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq -S {tag}.sam")
                    for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq"):
                        if os.path.exists(fname):
                            os.unlink(fname)
            else:
                f1, f2 = fpath.split(',')
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")
                if adaptor:
                    run_cmd(
                        f"cutadapt -j {thread} -m 50 --trim-n -a {adaptor} -A {adaptor} "
                        f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                        f"{tag}_R1.fastq {tag}_R2.fastq")
                    os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                    os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                run_cmd(
                    f"bowtie2 -p {thread} -x {genome}_chr_all "
                    f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq -S {tag}.sam")
                for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq"):
                    if os.path.exists(fname):
                        os.unlink(fname)

            # Process BAM
            run_cmd(
                f"samtools view -Sb -q 10 --threads {thread} {tag}.sam > {tag}.bam")
            os.unlink(f"{tag}.sam")
            run_cmd(f"samtools sort -n -o {tag}.sorted.name.bam {tag}.bam")
            run_cmd(f"samtools sort -o {tag}.sorted.bam {tag}.bam")
            run_cmd(f"samtools index {tag}.sorted.bam")
            os.unlink(f"{tag}.bam")
            run_cmd(
                f"bamCoverage -b {tag}.sorted.bam --skipNAs -bs 5 -p {thread} "
                f"--ignoreDuplicates --minMappingQuality 10 --normalizeUsing CPM -o {tag}.bw")
            tee.write("\nMapping completed!\n")

        # Cleanup
        for fname in globmod.glob(f"{genome}_chr_all*"):
            os.unlink(fname)
        if os.path.exists("igv.log"):
            os.unlink("igv.log")

        if not mappingonly:
            for g in treatment_groups:
                if peak_caller == 'macs3':
                    _run_macs3(g['ip_tags'], control_tags,
                               g['ip_tag'], seq_strategy, genome_size,
                               qvalue, pvalue, tee)
                else:
                    genrich_ip = "-t " + ','.join(
                        f"{t}.sorted.name.bam" for t in g['ip_tags'])
                    _run_genrich(genrich_ip, genrich_input, g['ip_tag'],
                                 seq_strategy, qvalue, pvalue, auc, tee)

    else:
        for pre in tags:
            os.symlink(f"../{pre}.sorted.name.bam", f"{pre}.sorted.name.bam")

        for g in treatment_groups:
            genrich_ip = "-t " + ','.join(
                f"{t}.sorted.name.bam" for t in g['ip_tags'])
            _run_genrich(genrich_ip, genrich_input, g['ip_tag'],
                         seq_strategy, qvalue, pvalue, auc, tee)

        for pre in tags:
            if os.path.exists(f"{pre}.sorted.name.bam"):
                os.unlink(f"{pre}.sorted.name.bam")


def _parse_to_dict(arg_str):
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _peak_qc(ip_tag, narrow_file, ip_bam_list, tee, caller='MACS3'):
    """Run peak QC: count, length, FRiP, and save report."""
    if not os.path.exists(narrow_file):
        tee.write(f"\nWarning: {narrow_file} not found, skipping QC.\n")
        return

    tee.write(f"\n{'='*60}\n")
    tee.write(f"Peak Quality Control — {ip_tag} ({caller})\n")
    tee.write(f"{'='*60}\n")

    # Parse peaks
    peak_count = 0
    total_length = 0
    try:
        with open(narrow_file) as f:
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= 3:
                    start, end = int(cols[1]), int(cols[2])
                    total_length += end - start
                    peak_count += 1
    except Exception as e:
        tee.write(f"Error reading peak file: {e}\n")
        return

    if peak_count == 0:
        tee.write("No peaks found.\n")
        return

    avg_length = total_length / peak_count

    # FRiP (每个生物学重复单独计算)
    tee.write("\nCalculating FRiP (Fraction of Reads in Peaks)...\n")
    frip_results = []
    total_reads_all = 0
    reads_in_peaks_all = 0
    for bam in ip_bam_list:
        if not os.path.exists(bam):
            continue
        label = os.path.basename(bam).replace('.sorted.bam', '').replace('.sorted.dedup.bam', '')
        res = subprocess.run(
            ['samtools', 'view', '-c', bam],
            capture_output=True, text=True)
        try:
            n_total = int(res.stdout.strip())
        except ValueError:
            continue
        res2 = subprocess.run(
            ['samtools', 'view', '-c', bam, '-L', narrow_file],
            capture_output=True, text=True)
        try:
            n_peaks = int(res2.stdout.strip())
        except ValueError:
            continue
        frip_rep = n_peaks / n_total if n_total > 0 else 0
        frip_results.append((label, n_total, n_peaks, frip_rep))
        total_reads_all += n_total
        reads_in_peaks_all += n_peaks

    frip_pooled = reads_in_peaks_all / total_reads_all if total_reads_all > 0 else 0

    # Output
    tee.write(f"\n{'─'*40}\n")
    tee.write(f"  Peak QC Summary\n")
    tee.write(f"{'─'*40}\n")
    tee.write(f"  Total peaks           : {peak_count}\n")
    tee.write(f"  Total peak length (bp): {total_length:,}\n")
    tee.write(f"  Average peak length   : {avg_length:.1f} bp\n")
    for label, n_total, n_peaks, frip_rep in frip_results:
        tee.write(f"\n  ── {label} ──\n")
        tee.write(f"    Total reads  : {n_total:,}\n")
        tee.write(f"    Reads in peaks: {n_peaks:,}\n")
        tee.write(f"    FRiP         : {frip_rep:.4f} ({frip_rep*100:.2f}%)\n")
    tee.write(f"\n  ── Pooled ──\n")
    tee.write(f"    FRiP         : {frip_pooled:.4f} ({frip_pooled*100:.2f}%)\n")
    tee.write(f"{'─'*40}\n")

    # Save
    qc_file = f"{ip_tag}_peak_qc.txt"
    with open(qc_file, 'w') as f:
        f.write(f"Sample\t{ip_tag}\n")
        f.write(f"Peak_caller\t{caller}\n")
        f.write(f"Total_peaks\t{peak_count}\n")
        f.write(f"Total_peak_length_bp\t{total_length}\n")
        f.write(f"Average_peak_length_bp\t{avg_length:.1f}\n")
        for label, n_total, n_peaks, frip_rep in frip_results:
            f.write(f"Replicate\t{label}\n")
            f.write(f"{label}_total_reads\t{n_total}\n")
            f.write(f"{label}_reads_in_peaks\t{n_peaks}\n")
            f.write(f"{label}_FRiP\t{frip_rep:.4f}\n")
        f.write(f"Pooled_FRiP\t{frip_pooled:.4f}\n")
    tee.write(f"QC report saved to: {qc_file}\n")


def _run_macs3(ip_tags, control_tags,
               ip_tag, seq_strategy, genome_size, qvalue, pvalue, tee):
    """Run MACS3 peak calling and generate peak QC report."""
    fmt = "BAMPE" if seq_strategy == 'paired' else "BAM"

    ip_bam_list = [f"{t}.sorted.bam" for t in ip_tags]
    ctrl_bam_list = [f"{t}.sorted.bam" for t in control_tags] if control_tags else []
    ip_bam = ' '.join(ip_bam_list)
    ctrl_bam = ' '.join(ctrl_bam_list) if ctrl_bam_list else ""

    # ── callpeak ─────────────────────────────────────────────────────────
    tee.write(f"\nMACS3 callpeak — {ip_tag}\n")
    cmd = f"macs3 callpeak -t {ip_bam} -f {fmt} -g {genome_size} -n {ip_tag}"
    if ctrl_bam:
        cmd += f" -c {ctrl_bam}"
    if qvalue < 1:
        tee.write(f"  Q-value threshold: {qvalue}\n")
        cmd += f" -q {qvalue}"
    else:
        tee.write(f"  P-value threshold: {pvalue}\n")
        cmd += f" -p {pvalue}"
    cmd += " --bdg"
    run_cmd(cmd)

    # ── Peak QC ──────────────────────────────────────────────────────────
    _peak_qc(ip_tag, f"{ip_tag}_peaks.narrowPeak", ip_bam_list, tee, caller='MACS3')

    # Cleanup bedGraph files
    for pat in ["*_treat_pileup.bdg", "*_control_lambda.bdg"]:
        for fname in globmod.glob(pat):
            os.unlink(fname)

    tee.write(f"\nMACS3 peak calling and QC completed.\n")


def _run_genrich(genrich_ip, genrich_input, ip_tag, seq_strategy, qvalue, pvalue, auc, tee):
    """Run Genrich peak calling + QC."""
    command = f"Genrich -r -v {genrich_ip} {genrich_input} -o {ip_tag}.narrowPeak.txt"

    if seq_strategy == 'single':
        tee.write("\nSingle-end mode is On...")
        command += " -y"

    if qvalue < 1:
        tee.write(f"\nFinding peaks...\nAUC\t{auc}\tQ Value\t{qvalue}\n")
        command += f" -q {qvalue}"
    else:
        tee.write(f"\nFinding peaks...\nAUC\t{auc}\tP Value\t{pvalue}\n")
        command += f" -p {pvalue}"

    run_cmd(command)

    # ── Peak QC ──────────────────────────────────────────────────────────
    # Reconstruct sorted.bam names from the genrich -t argument
    ip_bam_list = []
    for part in genrich_ip.replace('-t ', '').split(','):
        base = part.strip().replace('.sorted.name.bam', '')
        if base:
            ip_bam_list.append(f"{base}.sorted.bam")
    _peak_qc(ip_tag, f"{ip_tag}.narrowPeak.txt", ip_bam_list, tee, caller='Genrich')