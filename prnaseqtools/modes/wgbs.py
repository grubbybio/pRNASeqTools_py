"""
WGBS-seq (whole genome bisulfite sequencing) analysis mode.
Bismark alignment → methylation extraction → DMRcaller.
"""

import os
import sys
import glob as globmod
import subprocess
import math
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input
from prnaseqtools.functions import download_sra, unzip_file, _tee


def run(opts):
    """Main entry point for WGBS analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)
    binsize = opts.get('binsize', 100)
    min_c = opts.get('min_c', 4)

    control_dict = _parse_to_dict(opts.get('control', ''))
    tags, files, pars = parse_input(control_dict)

    if opts.get('treatment'):
        treatment_dict = _parse_to_dict(opts.get('treatment', ''))
        t_tags, t_files, t_pars = parse_input(treatment_dict)
        tags.extend(t_tags)
        files.extend(t_files)
        pars.extend(t_pars)

    par_str = ' '.join(pars)

    if not nomapping:
        tee.write("\nBuilding indices...\n")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        os.symlink(fasta_path, f"{genome}.fasta")
        subprocess.run("bismark_genome_preparation . 2>&1", shell=True, check=True)

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nMapping {tag}...\n")

            if ',' not in fpath:
                sra_results = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    unzip_file(sra_results[0], tag)
                    if adaptor:
                        tee.write("\nStart trimming...\r")
                        subprocess.run(
                            f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} "
                            f"-o {tag}_trimmed.fastq {tag}.fastq 2>&1",
                            shell=True, check=True
                        )
                        os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")

                    subprocess.run(
                        f"bismark -p {thread} -N 1 . {tag}.fastq 2>&1",
                        shell=True, check=True
                    )
                    subprocess.run(
                        f"deduplicate_bismark -s --bam {tag}_bismark_bt2.bam 2>&1",
                        shell=True, check=True
                    )
                    os.rename(f"{tag}_bismark_bt2.deduplicated.bam", f"{tag}.bam")
                    subprocess.run(
                        f"bismark_methylation_extractor --parallel {thread} -s --bedGraph "
                        f"--cutoff 4 --cytosine_report --CX --genome_folder . {tag}.bam 2>&1",
                        shell=True, check=True
                    )
                    for fname in (f"{tag}.fastq", f"{tag}_bismark_bt2.bam"):
                        if os.path.exists(fname):
                            os.unlink(fname)
                else:
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

                    subprocess.run(
                        f"bismark -p {thread} -N 1 . -1 {tag}_R1.fastq -2 {tag}_R2.fastq 2>&1",
                        shell=True, check=True
                    )
                    subprocess.run(
                        f"deduplicate_bismark -p --bam {tag}_R1_bismark_bt2_pe.bam 2>&1",
                        shell=True, check=True
                    )
                    os.rename(f"{tag}_R1_bismark_bt2_pe.deduplicated.bam", f"{tag}.bam")
                    subprocess.run(
                        f"bismark_methylation_extractor --parallel {thread} -p --bedGraph "
                        f"--cutoff 4 --cytosine_report --CX --genome_folder . {tag}.bam 2>&1",
                        shell=True, check=True
                    )
                    for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq", f"{tag}_R1_bismark_bt2_pe.bam"):
                        if os.path.exists(fname):
                            os.unlink(fname)
            else:
                f1, f2 = fpath.split(',')
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

                subprocess.run(
                    f"bismark -p {thread} -N 1 . -1 {tag}_R1.fastq -2 {tag}_R2.fastq 2>&1",
                    shell=True, check=True
                )
                subprocess.run(
                    f"deduplicate_bismark -p --bam {tag}_R1_bismark_bt2_pe.bam 2>&1",
                    shell=True, check=True
                )
                os.rename(f"{tag}_R1_bismark_bt2_pe.deduplicated.bam", f"{tag}.bam")
                subprocess.run(
                    f"bismark_methylation_extractor --parallel {thread} -p --bedGraph "
                    f"--cutoff 4 --cytosine_report --CX --genome_folder . {tag}.bam 2>&1",
                    shell=True, check=True
                )
                for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq", f"{tag}_R1_bismark_bt2_pe.bam"):
                    if os.path.exists(fname):
                        os.unlink(fname)

                # Generate bedgraph
                subprocess.run(
                    f"awk '{{OFS=\"\\t\";if($4+$5>0){{"
                    f"if($6==\"CG\"){{print $1,$2,$2+1,$4/($4+$5) > \"{tag}.CG.bedgraph\"}}; "
                    f"if($6==\"CHG\"){{print $1,$2,$2+1,$4/($4+$5) > \"{tag}.CHG.bedgraph\"}}; "
                    f"if($6==\"CHH\"){{print $1,$2,$2+1,$4/($4+$5) > \"{tag}.CHH.bedgraph\"}}}}}}' "
                    f"{tag}.CX_report.txt",
                    shell=True, check=True
                )
                subprocess.run(
                    f"sort -k1,1 -k2,2n {tag}.CX_report.txt > tmp",
                    shell=True, check=True
                )
                subprocess.run(f"bgzip -c tmp > {tag}.CX_report.txt.gz", shell=True, check=True)
                subprocess.run(f"tabix -C -p vcf {tag}.CX_report.txt.gz", shell=True, check=True)

            tee.write("\nAlignment finished...\n")
            _bin_methylation(tag, binsize, min_c, tee)

            for fname in globmod.glob(f"C*_O?_*{tag}*.txt"):
                os.unlink(fname)
            if os.path.exists("tmp"):
                os.unlink("tmp")

        if not mappingonly and len(pars) > 1:
            tee.write("\nPerforming DMRcaller...\n")
            subprocess.run(
                f"Rscript --vanilla {prefix}/scripts/DMRcaller.R {thread} {par_str}",
                shell=True, check=True
            )

        subprocess.run("rm -rf Bisulfite_Genome", shell=True)
        for fname in ([f"{genome}.fasta"] + globmod.glob("*.ebwt")):
            if os.path.exists(fname):
                os.unlink(fname)

    else:
        for pre in tags:
            os.symlink(f"../{pre}.CX_report.txt.gz", f"{pre}.CX_report.txt.gz")

        tee.write("\nPerforming DMRcaller...\n")
        subprocess.run(
            f"Rscript --vanilla {prefix}/scripts/DMRcaller.R {thread} {par_str}",
            shell=True, check=True
        )
        for fname in globmod.glob("*.CX_report.txt.gz"):
            os.unlink(fname)


def _parse_to_dict(arg_str):
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _bin_methylation(tag, binsize, min_c, tee):
    """Bin methylation data into windows."""
    contexts = ['CG', 'CHG', 'CHH']
    ccount_total = {c: 0 for c in contexts}
    ctcount_total = {c: 0 for c in contexts}
    cavg = {c: 0.0 for c in contexts}
    num = {c: 0 for c in contexts}
    covtotal = {c: 0 for c in contexts}
    sum_vals = {c: 0 for c in contexts}
    totalbin = {c: 0 for c in contexts}
    passedbin = {c: 0 for c in contexts}

    fhs = {}
    for cont, ctx in enumerate(contexts):
        fhs[ctx] = open(f"{tag}.bin.{binsize}.{ctx}.txt", 'w')
        fhs[ctx].write("bin\tccount\tctcount\tno.cytosin\tno.qualified.coverage\n")

    chr_name = None
    flag = binsize

    report_file = f"{tag}.CX_report.txt"
    if not os.path.exists(report_file):
        for fh in fhs.values():
            fh.close()
        return

    with open(report_file) as fh:
        for line_num, line in enumerate(fh):
            cols = line.strip().split('\t')
            if len(cols) < 6:
                continue

            if line_num == 0:
                chr_name = cols[0]

            curr_chr = cols[0]
            if curr_chr != chr_name:
                # Output bin
                for cont, ctx in enumerate(contexts):
                    avg_meth = 0.0
                    if num[ctx] > 0:
                        avg_meth = cavg[ctx] / num[ctx]
                        total_reads = ccount_total[ctx] + ctcount_total[ctx]
                        ccount_total[ctx] = math.floor(total_reads * avg_meth + 0.5)
                        ctcount_total[ctx] = math.floor(total_reads * (1 - avg_meth) + 0.5)

                    bin_idx = flag // binsize
                    fhs[ctx].write(f"{chr_name}_{bin_idx}\t{ccount_total[ctx]}\t{ctcount_total[ctx]}\t{num[ctx]}\t{covtotal[ctx]}\n")

                flag = binsize
                chr_name = curr_chr
                tee.write(f"Working on {chr_name}...\n")
                continue

            pos = int(cols[1])
            ctx_type = cols[5]
            ccount = int(cols[3])
            ctcount = int(cols[4])
            sum_vals[ctx_type] += 1

            while pos >= flag:
                for cont, ctx in enumerate(contexts):
                    avg_meth = 0.0
                    if num[ctx] > 0:
                        avg_meth = cavg[ctx] / num[ctx]
                        total_reads = ccount_total[ctx] + ctcount_total[ctx]
                        ccount_total[ctx] = math.floor(total_reads * avg_meth + 0.5)
                        ctcount_total[ctx] = math.floor(total_reads * (1 - avg_meth) + 0.5)

                    bin_idx = flag // binsize
                    fhs[ctx].write(f"{curr_chr}_{bin_idx}\t{ccount_total[ctx]}\t{ctcount_total[ctx]}\t{num[ctx]}\t{covtotal[ctx]}\n")
                    totalbin[ctx] += 1
                    if covtotal[ctx] >= 4:
                        passedbin[ctx] += 1

                    ccount_total[ctx] = 0
                    ctcount_total[ctx] = 0
                    cavg[ctx] = 0.0
                    num[ctx] = 0
                    covtotal[ctx] = 0
                flag += binsize

            if ccount + ctcount < min_c:
                continue

            ccount_total[ctx_type] += ccount
            ctcount_total[ctx_type] += ctcount
            num[ctx_type] += 1
            cur_avg = ccount / (ccount + ctcount) if (ccount + ctcount) > 0 else 0
            cavg[ctx_type] += cur_avg
            if ctcount + ccount >= 4:
                covtotal[ctx_type] += 1

    # Last bin
    for cont, ctx in enumerate(contexts):
        avg_meth = 0.0
        if num[ctx] > 0:
            avg_meth = cavg[ctx] / num[ctx]
            total_reads = ccount_total[ctx] + ctcount_total[ctx]
            ccount_total[ctx] = math.floor(total_reads * avg_meth + 0.5)
            ctcount_total[ctx] = math.floor(total_reads * (1 - avg_meth) + 0.5)

        bin_idx = flag // binsize
        fhs[ctx].write(f"{chr_name}_{bin_idx}\t{ccount_total[ctx]}\t{ctcount_total[ctx]}\t{num[ctx]}\t{covtotal[ctx]}\n")

    for ctx in contexts:
        tee.write(f"{tag}\t{ctx}\t{totalbin[ctx]}\t{passedbin[ctx]}\n")
        fhs[ctx].close()
