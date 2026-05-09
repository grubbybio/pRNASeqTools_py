"""
PhasiRNA analysis mode.
ShortStack alignment → phasing score calculation.
"""

import os
import sys
import glob as globmod
import subprocess
import math
from pathlib import Path
from collections import defaultdict

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input
from prnaseqtools.functions import download_sra, unzip_file, _tee
from prnaseqtools import reference as ref


def run(opts):
    """Main entry point for phasiRNA analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    mmap = opts.get('mmap', 'u')
    nomapping = opts.get('no_mapping', False)
    binsize = opts.get('binsize', 100)
    norm = opts.get('norm', 'rRNA,total')
    norms = norm.split(',')
    period = opts.get('period', 21)
    phasing_score_cutoff = opts.get('phasingscore', 50)

    # Parse inputs
    control_dict = _parse_to_dict(opts.get('control', ''))
    tags, files, pars = parse_input(control_dict)

    if opts.get('treatment'):
        treatment_dict = _parse_to_dict(opts.get('treatment', ''))
        t_tags, t_files, t_pars = parse_input(treatment_dict)
        tags.extend(t_tags)
        files.extend(t_files)
        pars.extend(t_pars)

    if not nomapping:
        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nMapping {tag}...\n")

            sra_results = download_sra(fpath, thread)
            unzip_file(sra_results[0], tag)

            if adaptor:
                tee.write(f"\nTrimming {tag}...\n")
                subprocess.run(
                    f"cutadapt -j {thread} -m 18 -M 42 --discard-untrimmed --trim-n "
                    f"-a {adaptor} -o {tag}_trimmed.fastq {tag}.fastq 2>&1",
                    shell=True, check=True
                )
                os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")

            # rRNA filtering
            lsu_rRNA = os.path.join(prefix, "reference", "lsu_rrna")
            subprocess.run(
                f"bowtie -v 2 -a -p {thread} -t {lsu_rRNA} "
                f"{tag}.fastq {tag}.rRNA.out 2>&1",
                shell=True, check=True
            )
            subprocess.run(
                f"awk -F \"\\t\" 'BEGIN{{x=0}}{{if(/{genome}/){{x++}}}}END{{print \"rRNA\\t\"x}}' "
                f"{tag}.rRNA.out > {tag}.nf",
                shell=True, check=True
            )

            # ShortStack alignment
            fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
            subprocess.run(
                f"ShortStack --outdir ShortStack_{tag} --align_only --bowtie_m 1000 "
                f"--ranmax 50 --mmap {mmap} --mismatches 0 --bowtie_cores {thread} "
                f"--nohp --readfile {tag}.fastq --genomefile {fasta_path} 2>&1",
                shell=True, check=True
            )

            for fname in globmod.glob(f"{tag}*.fastq"):
                os.unlink(fname)
            if os.path.exists(f"{tag}.rRNA.out"):
                os.unlink(f"{tag}.rRNA.out")

            tee.write("\nAlignment Completed!\n")

            # Process SAM/BAM
            subprocess.run(f"samtools view -h ShortStack_{tag}/{tag}.bam > {tag}", shell=True, check=True)
            subprocess.run(
                f"awk '{{if($0~/^@/) print > (FILENAME\".unmapped.sam\"); "
                f"if($10!=\"*\" && $3!=\"*\") print > (FILENAME\".sam\"); "
                f"if($10!=\"*\" && $3==\"*\") print > (FILENAME\".unmapped.sam\")}}' {tag}",
                shell=True, check=True
            )
            subprocess.run(f"samtools view -Sb {tag}.unmapped.sam > {tag}.unmapped.bam", shell=True, check=True)
            subprocess.run(f"samtools view -Sb {tag}.sam > {tag}.bam", shell=True, check=True)
            subprocess.run(
                f"awk '{{n++}}END{{print \"total\\t\"n}}' {tag}.sam >> {tag}.nf",
                shell=True, check=True
            )

            for fname in (tag, f"{tag}.unmapped.sam", f"{tag}.sam"):
                if os.path.exists(fname):
                    os.unlink(fname)
            if os.path.exists(f"ShortStack_{tag}"):
                subprocess.run(f"rm -rf ShortStack_{tag}", shell=True)

        _phasi_analysis(pars, period, norms, prefix, genome, binsize, phasing_score_cutoff, tee)
    else:
        for pre in tags:
            os.symlink(f"../{pre}.nf", f"{pre}.nf")
            os.symlink(f"../{pre}.bam", f"{pre}.bam")

        _phasi_analysis(pars, period, norms, prefix, genome, binsize, phasing_score_cutoff, tee)

        for fname in globmod.glob("*.bam"):
            os.unlink(fname)
        for fname in globmod.glob("*.nf"):
            os.unlink(fname)


def _parse_to_dict(arg_str):
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _phasi_analysis(pars, period, norms, prefix, genome, binsize, phasing_score_cutoff, tee):
    """Calculate phasing scores for each bin."""
    p = list(pars)
    while len(p) >= 2:
        sample = p.pop(0)
        num = int(p.pop(0))

        tee.write("\nphasiRNA analysis...\n")

        count = 0
        data1 = defaultdict(lambda: defaultdict(int))  # forward strand
        data2 = defaultdict(lambda: defaultdict(int))  # reverse strand
        data0 = defaultdict(lambda: defaultdict(int))  # all reads
        normc = defaultdict(int)

        # Merge BAMs
        command = f"samtools merge {sample}.bam "
        for j in range(1, num + 1):
            with open(f"{sample}_{j}.nf") as fh:
                for line in fh:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        normc[parts[0]] += int(parts[1])
            command += f"{sample}_{j}.bam "

        subprocess.run(command, shell=True, check=True)
        subprocess.run(f"samtools view -h {sample}.bam > {sample}.sam", shell=True, check=True)

        tee.write("Reading SAM ...\n")

        with open(f"{sample}.sam") as sam_fh, open(f"{sample}.phasi.sam", 'w') as pha_fh:
            for line in sam_fh:
                line = line.strip()
                if not line or line.startswith('@'):
                    pha_fh.write(line + '\n')
                    continue

                cols = line.split('\t')
                if len(cols) < 6:
                    continue

                read_id, flag, chr_name, pos, mapq, cigar = cols[0], int(cols[1]), cols[2], int(cols[3]), cols[4], cols[5]

                data0[chr_name][pos] += 1

                if cigar == f"{period}M":
                    pha_fh.write(line + '\n')
                    if flag == 0:
                        data1[chr_name][pos][0] += 1
                        count += 1
                    elif flag == 16:
                        data2[chr_name][pos][16] += 1
                        count += 1

                if count > 0 and count % 100000 == 0:
                    tee.write(f"Read counts: {count}\r")

        subprocess.run(
            f"samtools view -Sb {sample}.phasi.sam > {sample}.phasi.bam",
            shell=True, check=True
        )
        for fname in (f"{sample}.phasi.sam", f"{sample}.sam", f"{sample}.bam"):
            if os.path.exists(fname):
                os.unlink(fname)

        tee.write(f"Read counts: {count}\nCombining strands ...\n")

        # Combine strands with 2-nt overhang
        for chr_name in data2:
            for pos in data2[chr_name]:
                data1[chr_name][pos + 2][0] += data2[chr_name][pos][16]

        for mnorm in norms:
            data3 = defaultdict(lambda: defaultdict(int))
            rc = normc.get(mnorm, 1)
            if rc == 0:
                rc = 1

            for chr_name in data1:
                for pos in data1[chr_name]:
                    data3[chr_name][pos][0] = data1[chr_name][pos][0] * 1000000.0 / rc

            tee.write("Calculating phasing score ...\n")

            phase_res = defaultdict(lambda: defaultdict(float))
            hash2 = defaultdict(lambda: defaultdict(float))

            for chr_name in sorted(data3.keys()):
                tee.write(f"{chr_name} ...\n")
                for pos in sorted(data3[chr_name].keys()):
                    n = 0
                    phased = 0.0
                    unphased = 0.0
                    xx = 0

                    for cycle in range(10):
                        new_pos = pos + period * cycle
                        xx += data0[chr_name].get(new_pos, 0)

                        if data3[chr_name].get(new_pos, {}).get(0, 0) > 0:
                            n += 1
                            phased += data3[chr_name][new_pos][0]

                        for unphased_nt in range(1, period):
                            new_pos += 1
                            xx += data0[chr_name].get(new_pos, 0)
                            if data3[chr_name].get(new_pos, {}).get(0, 0) > 0:
                                n += 1
                                unphased += data3[chr_name][new_pos][0]

                    ratio = (phased + unphased) / ((xx + 1) * 1000000.0 / rc)
                    phasing_score = 0

                    if n >= 3 and phased >= 10 and ratio >= 0.25:
                        phasing_score = (n - 2) * math.log(1 + 10 * phased / (unphased + 1))

                    phase_res[chr_name][pos] = phasing_score

            tee.write("Outputting result ...\n")

            with open(f"{sample}.{mnorm}.phasiRNA.txt", 'w') as fh_out, \
                 open(f"{sample}.{mnorm}.phasiRNA.bedgraph", 'w') as fh_bg, \
                 open(f"{sample}.{mnorm}.plot.phasiRNA.bedgraph", 'w') as fh_plot:

                fh_out.write("#BIN,Phase_score\n")
                for chr_name in sorted(phase_res.keys()):
                    for pos in sorted(phase_res[chr_name].keys()):
                        score = phase_res[chr_name][pos]
                        if score > 0:
                            fh_plot.write(f"{chr_name}\t{pos}\t{pos}\t{score}\n")
                            if score >= phasing_score_cutoff and data3[chr_name].get(pos, {}).get(0, 0) >= 1:
                                fh_bg.write(f"{chr_name}\t{pos}\t{pos}\t{score}\n")
                                bin_idx = pos // binsize
                                if bin_idx not in hash2[chr_name] or score > hash2[chr_name][bin_idx]:
                                    hash2[chr_name][bin_idx] = score

                for chr_name in sorted(hash2.keys()):
                    for bin_idx in sorted(hash2[chr_name].keys()):
                        fh_out.write(f"{chr_name}_{bin_idx},{hash2[chr_name][bin_idx]}\n")

            # Annotate results
            ann = ref.build_annotation(prefix, genome, binsize)
            _annotate_result(f"{sample}.{mnorm}.phasiRNA.txt", ann)

        tee.write("\nphasiRNA analysis complete!\n\n")


def _annotate_result(result_file, ann):
    """Add genome annotation to phasiRNA result file."""
    tmp_file = "tmp3"
    with open(result_file) as fh_in, open(tmp_file, 'w') as fh_out:
        for line in fh_in:
            line = line.strip()
            if not line or line.startswith('#'):
                fh_out.write(line + '\n')
                continue
            cols = line.split(',')
            if not cols:
                continue
            key = cols[0].strip('"')
            if key in ann:
                fh_out.write(f"{line},{ann[key]}\n")
            else:
                fh_out.write(f"{line},NA\n")
    os.rename(tmp_file, result_file)
