"""
TS-CLIP-seq analysis mode.
Similar to CLIP-seq but with 3' PolyC removal before mapping.
Reuses CLIP mode's _differential_analysis().
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input
from prnaseqtools.functions import download_sra, unzip_file, rmvc, _tee
from prnaseqtools.modes.clip import _differential_analysis


def run(opts):
    """Main entry point for TS-CLIP-seq analysis."""
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

    control_dict = _parse_to_dict(opts.get('control', ''))
    tags, files, pars = parse_input(control_dict)

    if opts.get('treatment'):
        treatment_dict = _parse_to_dict(opts.get('treatment', ''))
        t_tags, t_files, t_pars = parse_input(treatment_dict)
        tags.extend(t_tags)
        files.extend(t_files)
        pars.extend(t_pars)

    if not nomapping:
        if not adaptor:
            sys.exit("Please specify the 3' adaptor!")

        tee.write("\nBuilding STAR genome index ...\n")
        if os.path.exists("Genome"):
            subprocess.run("rm -rf Genome", shell=True)
        os.makedirs("Genome", exist_ok=True)

        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")

        subprocess.run(
            f"STAR --runThreadN {thread} --genomeDir Genome --runMode genomeGenerate "
            f"--genomeFastaFiles {fasta_path} --sjdbGTFfile {gff_path} "
            f"--sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene ID",
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
                    subprocess.run(
                        f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} "
                        f"-o {tag}_trimmed.fastq {tag}.fastq 2>&1",
                        shell=True, check=True
                    )
                    rmvc(tag)

                    subprocess.run(
                        f"STAR --genomeDir Genome --alignIntronMax 5000 "
                        f"--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 "
                        f"--outReadsUnmapped Fastx --outSAMmultNmax 1 "
                        f"--outFilterMismatchNoverLmax 0.1 --runThreadN {thread} "
                        f"--readFilesIn {tag}.fastq 2>&1",
                        shell=True, check=True
                    )
                    if os.path.exists("Unmapped.out.mate1"):
                        os.rename("Unmapped.out.mate1", f"{tag}.unmapped.fastq")
                    for fname in (f"{tag}.fastq", f"{tag}_trimmed.fastq"):
                        if os.path.exists(fname):
                            os.unlink(fname)
                else:
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    subprocess.run(
                        f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} -A AAAAAAGAAAAAA "
                        f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                        f"{tag}_R1.fastq {tag}_R2.fastq 2>&1",
                        shell=True, check=True
                    )
                    rmvc(f"{tag}_R1", f"{tag}_R2")

                    subprocess.run(
                        f"STAR --genomeDir Genome --alignIntronMax 5000 "
                        f"--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 "
                        f"--outReadsUnmapped Fastx --outSAMmultNmax 1 "
                        f"--outFilterMismatchNoverLmax 0.1 --runThreadN {thread} "
                        f"--readFilesIn {tag}_R1.fastq {tag}_R2.fastq 2>&1",
                        shell=True, check=True
                    )
                    if os.path.exists("Unmapped.out.mate1"):
                        os.rename("Unmapped.out.mate1", f"{tag}.unmapped_R1.fastq")
                    if os.path.exists("Unmapped.out.mate2"):
                        os.rename("Unmapped.out.mate2", f"{tag}.unmapped_R2.fastq")
                    for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                                  f"{tag}_R1_trimmed.fastq", f"{tag}_R2_trimmed.fastq"):
                        if os.path.exists(fname):
                            os.unlink(fname)
            else:
                f1, f2 = fpath.split(',')
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")
                subprocess.run(
                    f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} -A AAAAAAGAAAAAA "
                    f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                    f"{tag}_R1.fastq {tag}_R2.fastq 2>&1",
                    shell=True, check=True
                )
                rmvc(f"{tag}_R1", f"{tag}_R2")

                subprocess.run(
                    f"STAR --genomeDir Genome --alignIntronMax 5000 "
                    f"--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 "
                    f"--outReadsUnmapped Fastx --outSAMmultNmax 1 "
                    f"--outFilterMismatchNoverLmax 0.1 --runThreadN {thread} "
                    f"--readFilesIn {tag}_R1.fastq {tag}_R2.fastq 2>&1",
                    shell=True, check=True
                )
                if os.path.exists("Unmapped.out.mate1"):
                    os.rename("Unmapped.out.mate1", f"{tag}.unmapped_R1.fastq")
                if os.path.exists("Unmapped.out.mate2"):
                    os.rename("Unmapped.out.mate2", f"{tag}.unmapped_R2.fastq")
                for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                              f"{tag}_R1_trimmed.fastq", f"{tag}_R2_trimmed.fastq"):
                    if os.path.exists(fname):
                        os.unlink(fname)

            os.rename("Aligned.sortedByCoord.out.bam", f"{tag}.bam")
            if os.path.exists("Log.final.out"):
                with open("Log.final.out") as lf:
                    tee.write(lf.read())

            subprocess.run(f"samtools index {tag}.bam", shell=True, check=True)
            subprocess.run(
                f"bamCoverage -b {tag}.bam --skipNAs -bs 5 -p {thread} "
                f"--minMappingQuality 10 --ignoreDuplicates --normalizeUsing CPM -o {tag}.bw",
                shell=True, check=True
            )

            tee.write("\nFinding peaks...\n")
            subprocess.run(
                f"clipper -b {tag}.bam -s {genome} --FDR=0.01 --minreads=2 "
                f"--processors={thread} --threshold-method=binomial --min_width=20 "
                f"-o {tag}.fitted_clusters.bed -v 2>&1",
                shell=True, check=True
            )

        for fname in ("Log.out", "Log.progress.out", "Log.final.out", "SJ.out.tab"):
            if os.path.exists(fname):
                os.unlink(fname)
        if os.path.exists("Genome"):
            subprocess.run("rm -rf Genome", shell=True)

        if not mappingonly and len(pars) > 1:
            _ts_ana(pars, prefix, pvalue, foldchange, tee)
    else:
        for pre in tags:
            os.symlink(f"../{pre}.nf", f"{pre}.nf")
            os.symlink(f"../{pre}.bam", f"{pre}.bam")
            os.symlink(f"../{pre}.fitted_clusters.bed", f"{pre}.fitted_clusters.bed")

        _ts_ana(pars, prefix, pvalue, foldchange, tee)

        for fname in globmod.glob("*.bed"):
            os.unlink(fname)
        for fname in globmod.glob("*.nf") + globmod.glob("*.bam*"):
            os.unlink(fname)


def _parse_to_dict(arg_str):
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


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

    subprocess.run(
        f"cat {command} | sort -k1,1 -k2,2n | bedtools merge "
        f"-c 4 -o collapse -s -i - > tmp.bed",
        shell=True, check=True
    )

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
                    f"{cols[0]}\t{cols[1]}\t{cols[2]}\t{output_gene}\t{total_sum}\t{cols[3]}\n"
                )

    tee.write("\nCounting reads in each sample...\n")
    for tag in all_tags:
        subprocess.run(
            f"bedtools intersect -c -s -a ref.bed -b {tag}.bam | "
            f"awk '{{print $1\"_\"$2\"_\"$3\"_\"$6\"_\"$4\"\\t\"$7}}' > {tag}.txt",
            shell=True, check=True
        )

    tee.write("\nFinding differential peaks...\n")
    subprocess.run(
        f"Rscript --vanilla {prefix}/scripts/CLIP.R {pvalue} {foldchange} "
        + ' '.join(p_list),
        shell=True, check=True
    )

    for fname in ("tmp.bed", "ref.bed"):
        if os.path.exists(fname):
            os.unlink(fname)
