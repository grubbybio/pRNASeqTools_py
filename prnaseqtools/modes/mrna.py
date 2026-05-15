"""
mRNA-seq analysis mode.
Mirrors the Perl mode::mrna module.
STAR alignment → featureCounts → DESeq2 DE analysis.
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
    """Main entry point for mRNA-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    foldchange = opts.get('foldchange', 2.0)
    pvalue = opts.get('pvalue', 0.01)
    fdr = opts.get('fdr', 1.0)
    run_mode = opts.get('run_mode', 1)
    control = opts.get('control')
    treatment = opts.get('treatment')
    norm = opts.get('deseq2_norm', 'DESeq2')
    mask = opts.get('mask')
    total = opts.get('total', False)
    seq_strategy = opts.get('seq_strategy')
    genome_size = opts.get('genome_size', 10)

    mapping = run_mode not in (3,)
    do_count = run_mode != 4
    do_de = run_mode != 2

    # Parse inputs
    control_dict = _parse_to_dict(control)
    tags, files, pars = parse_input(control_dict)

    if treatment:
        treatment_dict = _parse_to_dict(treatment)
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
            subprocess.run("bowtie-build -q mask.fa mask", shell=True, check=True)

        tee.write("\nBuilding STAR genome index ...\n")

        if os.path.exists("Genome"):
            subprocess.run("rm -rf Genome", shell=True)
        os.makedirs("Genome", exist_ok=True)

        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")

        subprocess.run(
            f"STAR --runThreadN {thread} --genomeDir Genome --runMode genomeGenerate "
            f"--genomeSAindexNbases {genome_size} --genomeFastaFiles {fasta_path} "
            f"--sjdbGTFfile {gff_path} --sjdbGTFtagExonParentTranscript Parent "
            f"--sjdbGTFtagExonParentGene ID --limitGenomeGenerateRAM 64000000000",
            shell=True, check=True
        )

        if total:
            subprocess.run(
                f"gffread -T -o {genome}_genes.gtf -g {fasta_path} {gff_path}",
                shell=True, check=True
            )
        else:
            subprocess.run(
                f"gffread -T -C -o {genome}_genes.gtf -g {fasta_path} {gff_path}",
                shell=True, check=True
            )

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nProcessing {tag}...\n")

            if ',' not in fpath:
                # Single-end or paired from SRA
                sra_results = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    seq_strategy = 'single'
                    unzip_file(sra_results[0], tag)
                    if adaptor:
                        subprocess.run(
                            f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} "
                            f"-o {tag}_trimmed.fastq {tag}.fastq 2>&1",
                            shell=True, check=True
                        )
                        os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")
                    if mask:
                        subprocess.run(
                            f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                            f"{tag}.fastq {tag}.mask.out 2>&1",
                            shell=True, check=True
                        )
                        os.rename("tmp.fastq", f"{tag}.fastq")
                        if os.path.exists(f"{tag}.mask.out"):
                            os.unlink(f"{tag}.mask.out")

                    subprocess.run(
                        f"STAR --runMode alignReads --genomeDir Genome --alignIntronMax 5000 "
                        f"--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate "
                        f"--limitBAMsortRAM 10000000000 --outSAMmultNmax 1 "
                        f"--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 "
                        f"--runThreadN {thread} --readFilesIn {tag}.fastq 2>&1",
                        shell=True, check=True
                    )
                    if os.path.exists(f"{tag}.fastq"):
                        os.unlink(f"{tag}.fastq")
                else:
                    # Paired from SRA
                    seq_strategy = 'paired'
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    if adaptor:
                        subprocess.run(
                            f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} -A {adaptor} "
                            f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                            f"{tag}_R1.fastq {tag}_R2.fastq 2>&1",
                            shell=True, check=True
                        )
                        os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                        os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                    if mask:
                        subprocess.run(
                            f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                            f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq {tag}.mask.out 2>&1",
                            shell=True, check=True
                        )
                        if os.path.exists("tmp_1.fastq"):
                            os.rename("tmp_1.fastq", f"{tag}_R1.fastq")
                        if os.path.exists("tmp_2.fastq"):
                            os.rename("tmp_2.fastq", f"{tag}_R2.fastq")
                        if os.path.exists(f"{tag}.mask.out"):
                            os.unlink(f"{tag}.mask.out")

                    subprocess.run(
                        f"STAR --runMode alignReads --genomeDir Genome --alignIntronMax 5000 "
                        f"--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate "
                        f"--limitBAMsortRAM 10000000000 --outSAMmultNmax 1 "
                        f"--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 "
                        f"--runThreadN {thread} --readFilesIn {tag}_R1.fastq {tag}_R2.fastq 2>&1",
                        shell=True, check=True
                    )
                    for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq"):
                        if os.path.exists(fname):
                            os.unlink(fname)
            else:
                # Explicit paired-end
                f1, f2 = fpath.split(',')
                seq_strategy = 'paired'
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")
                if adaptor:
                    subprocess.run(
                        f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} -A {adaptor} "
                        f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                        f"{tag}_R1.fastq {tag}_R2.fastq 2>&1",
                        shell=True, check=True
                    )
                    os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                    os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                if mask:
                    subprocess.run(
                        f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                        f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq {tag}.mask.out 2>&1",
                        shell=True, check=True
                    )
                    if os.path.exists("tmp_1.fastq"):
                        os.rename("tmp_1.fastq", f"{tag}_R1.fastq")
                    if os.path.exists("tmp_2.fastq"):
                        os.rename("tmp_2.fastq", f"{tag}_R2.fastq")
                    if os.path.exists(f"{tag}.mask.out"):
                        os.unlink(f"{tag}.mask.out")

                tee.write("\nMapping...\n")
                subprocess.run(
                    f"STAR --runMode alignReads --genomeDir Genome --alignIntronMax 5000 "
                    f"--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 "
                    f"--outSAMmultNmax 1 --outFilterMultimapNmax 50 "
                    f"--outFilterMismatchNoverLmax 0.1 --runThreadN {thread} "
                    f"--readFilesIn {tag}_R1.fastq {tag}_R2.fastq 2>&1",
                    shell=True, check=True
                )
                for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq"):
                    if os.path.exists(fname):
                        os.unlink(fname)

            os.rename("Aligned.sortedByCoord.out.bam", f"{tag}.bam")

            # Show mapping stats
            if os.path.exists("Log.final.out"):
                with open("Log.final.out") as lf:
                    tee.write(lf.read())

            _count_bam(tag, thread, seq_strategy, prefix, genome)

        # Cleanup
        for fname in ("Log.out", "Log.progress.out", "Log.final.out", "SJ.out.tab",
                      f"{genome}_genes.gtf"):
            if os.path.exists(fname):
                os.unlink(fname)
        for fname in globmod.glob("total.count*"):
            os.unlink(fname)
        if mask:
            for mf in globmod.glob("mask*"):
                os.unlink(mf)
        if os.path.exists("Genome"):
            subprocess.run("rm -rf Genome", shell=True)

        if do_de:
            tee.write(f"\nFinding DEG...\nFold Change\t{foldchange}\tP Value\t{pvalue}\tFDR\t{fdr}\n")
            subprocess.run(
                f"Rscript --vanilla {prefix}/scripts/DEG.R {norm} {pvalue} {fdr} "
                f"{foldchange} {prefix} {genome} {par_str}",
                shell=True, check=True
            )

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
                subprocess.run(
                    f"gffread -T -o {genome}_genes.gtf -g {fasta_path} {gff_path}",
                    shell=True, check=True
                )
            else:
                gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
                fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
                subprocess.run(
                    f"gffread -T -C -o {genome}_genes.gtf -g {fasta_path} {gff_path}",
                    shell=True, check=True
                )
            for pre in tags:
                _count_bam(pre, thread, seq_strategy, prefix, genome)
            if os.path.exists(f"{genome}_genes.gtf"):
                os.unlink(f"{genome}_genes.gtf")

        tee.write(f"\nFinding DEG...\nFold Change\t{foldchange}\tP Value\t{pvalue}\tFDR\t{fdr}\n")
        subprocess.run(
            f"Rscript --vanilla {prefix}/scripts/DEG.R {norm} {pvalue} {fdr} "
            f"{foldchange} {prefix} {genome} {par_str}",
            shell=True, check=True
        )

        # Clean up symlinks
        if do_count:
            for fname in globmod.glob("*_?.bam"):
                os.unlink(fname)
        else:
            for fname in globmod.glob("*_?.txt"):
                os.unlink(fname)


def _parse_to_dict(arg_str):
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _resolve_path(filepath):
    if filepath.startswith('~/'):
        return os.path.expanduser(filepath)
    elif not filepath.startswith('/'):
        return os.path.abspath(os.path.join('..', filepath))
    return filepath


def _count_bam(tag, thread, seq_strategy, prefix, genome):
    """Count reads in BAM using featureCounts."""
    subprocess.run(f"samtools index {tag}.bam", shell=True, check=True)
    subprocess.run(
        f"bamCoverage -b {tag}.bam --skipNAs -bs 5 -p {thread} "
        f"--minMappingQuality 10 --normalizeUsing CPM -o {tag}.bw",
        shell=True, check=True
    )

    tee = _tee()
    tee.write("\nStart counting...\n")

    fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
    gtf_file = f"{genome}_genes.gtf"

    if seq_strategy == 'single':
        subprocess.run(
            f"featureCounts -T {thread} -O -G {fasta_path} -s 0 "
            f"-a {gtf_file} -o total.count {tag}.bam 2>&1",
            shell=True, check=True
        )
    elif seq_strategy == 'paired':
        subprocess.run(
            f"featureCounts -T {thread} -p --countReadPairs -BCO -G {fasta_path} "
            f"-s 0 -a {gtf_file} -o total.count {tag}.bam 2>&1",
            shell=True, check=True
        )
    else:
        # Default to single
        subprocess.run(
            f"featureCounts -T {thread} -O -G {fasta_path} -s 0 "
            f"-a {gtf_file} -o total.count {tag}.bam 2>&1",
            shell=True, check=True
        )

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
