"""
ATAC-seq analysis mode.
bowtie2 alignment → Genrich peak calling (ATAC-specific parameters).
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input
from prnaseqtools.functions import download_sra, unzip_file, _tee


def run(opts):
    """Main entry point for ATAC-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    auc = opts.get('auc', 20)
    qvalue = opts.get('qvalue', 1.0)
    pvalue = opts.get('pvalue', 0.01)
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)

    tags = []
    files = []

    input_opt = opts.get('control')
    genrich_input = ""
    if input_opt:
        input_dict = _parse_to_dict(input_opt)
        i_tags, i_files, _ = parse_input(input_dict)
        tags.extend(i_tags)
        files.extend(i_files)
        genrich_input = "-c " + ','.join(f"{t}.sorted.name.bam" for t in i_tags)

    ip_opt = opts.get('treatment', '')
    ip_dict = _parse_to_dict(ip_opt)
    i_tags, i_files, _ = parse_input(ip_dict)
    tags.extend(i_tags)
    files.extend(i_files)
    genrich_ip = "-t " + ','.join(f"{t}.sorted.name.bam" for t in i_tags)
    ip_tag = i_tags[0] if i_tags else 'ip'

    if not nomapping:
        tee.write("\nBuilding index...\n")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        os.symlink(fasta_path, f"{genome}_chr_all.fasta")
        subprocess.run(
            f"bowtie2-build -q {genome}_chr_all.fasta {genome}_chr_all",
            shell=True, check=True
        )

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nMapping {tag}...\n")

            if ',' not in fpath:
                sra_results = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    unzip_file(sra_results[0], tag)
                    if adaptor:
                        tee.write("\nTrimming...\n")
                        subprocess.run(
                            f"cutadapt -j {thread} -m 50 --trim-n -a {adaptor} "
                            f"-o {tag}_trimmed.fastq {tag}.fastq 2>&1",
                            shell=True, check=True
                        )
                        os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")
                    subprocess.run(
                        f"bowtie2 -p {thread} -x {genome}_chr_all -U {tag}.fastq "
                        f"-S {tag}.sam 2>&1",
                        shell=True, check=True
                    )
                    if os.path.exists(f"{tag}.fastq"):
                        os.unlink(f"{tag}.fastq")
                else:
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    if adaptor:
                        subprocess.run(
                            f"cutadapt -j {thread} -m 50 --trim-n -a {adaptor} -A {adaptor} "
                            f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                            f"{tag}_R1.fastq {tag}_R2.fastq 2>&1",
                            shell=True, check=True
                        )
                        os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                        os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                    subprocess.run(
                        f"bowtie2 -p {thread} -x {genome}_chr_all "
                        f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq -S {tag}.sam 2>&1",
                        shell=True, check=True
                    )
                    for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq"):
                        if os.path.exists(fname):
                            os.unlink(fname)
            else:
                f1, f2 = fpath.split(',')
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")
                if adaptor:
                    subprocess.run(
                        f"cutadapt -j {thread} -m 50 --trim-n -a {adaptor} -A {adaptor} "
                        f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                        f"{tag}_R1.fastq {tag}_R2.fastq 2>&1",
                        shell=True, check=True
                    )
                    os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                    os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                subprocess.run(
                    f"bowtie2 -p {thread} -x {genome}_chr_all "
                    f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq -S {tag}.sam 2>&1",
                    shell=True, check=True
                )
                for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq"):
                    if os.path.exists(fname):
                        os.unlink(fname)

            subprocess.run(
                f"samtools view -Sb -q 10 --threads {thread} {tag}.sam > {tag}.bam",
                shell=True, check=True
            )
            os.unlink(f"{tag}.sam")
            subprocess.run(f"samtools sort -n -o {tag}.sorted.name.bam {tag}.bam", shell=True, check=True)
            subprocess.run(f"samtools sort -o {tag}.sorted.bam {tag}.bam", shell=True, check=True)
            subprocess.run(f"samtools index {tag}.sorted.bam", shell=True, check=True)
            os.unlink(f"{tag}.bam")
            subprocess.run(
                f"bamCoverage -b {tag}.sorted.bam --skipNAs -bs 5 -p {thread} "
                f"--ignoreDuplicates --minMappingQuality 10 --normalizeUsing CPM -o {tag}.bw",
                shell=True, check=True
            )
            tee.write("\nMapping completed!\n")

        for fname in globmod.glob(f"{genome}_chr_all*"):
            os.unlink(fname)
        if os.path.exists("igv.log"):
            os.unlink("igv.log")

        if not mappingonly and input_opt:
            _run_genrich(genrich_ip, genrich_input, ip_tag, qvalue, pvalue, auc, tee)

    else:
        for pre in tags:
            os.symlink(f"../{pre}.sorted.name.bam", f"{pre}.sorted.name.bam")

        _run_genrich(genrich_ip, genrich_input, ip_tag, qvalue, pvalue, auc, tee)

        for pre in tags:
            if os.path.exists(f"{pre}.sorted.name.bam"):
                os.unlink(f"{pre}.sorted.name.bam")


def _parse_to_dict(arg_str):
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _run_genrich(genrich_ip, genrich_input, ip_tag, qvalue, pvalue, auc, tee):
    """Run Genrich for ATAC peak calling."""
    command = (
        f"Genrich -j -y -r -v -a 20 -e chrC,chrM {genrich_ip} {genrich_input} "
        f"-o {ip_tag}.ATAC.narrowPeak.txt"
    )

    if qvalue < 1:
        tee.write(f"\nFinding peaks...\nAUC\t{auc}\tQ Value\t{qvalue}\n")
        command += f" -q {qvalue}"
    else:
        tee.write(f"\nFinding peaks...\nAUC\t{auc}\tP Value\t{pvalue}\n")
        command += f" -p {pvalue}"

    subprocess.run(command + " 2>&1", shell=True, check=True)
