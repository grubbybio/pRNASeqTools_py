"""
RiboMeth-seq analysis mode.
HISAT2 alignment to transcriptome → RNAmodR analysis.
"""

import os
import re
import sys
import glob as globmod
import subprocess
from pathlib import Path
from collections import defaultdict

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import (parse_input, _parse_to_dict,
                                        _resolve_path)
from prnaseqtools.functions import (download_sra, unzip_file, _tee, run_cmd,
                                     gzip_fastq, preserve_raw_fastq,
                                     cleanup_alignment_outputs,
                                     hisat2_build_index, hisat2_align_single,
                                     hisat2_align_paired)
from prnaseqtools import reference as ref


def run(opts):
    """Main entry point for RiboMeth-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    nomapping = opts.get('no_mapping', False)
    reference = opts.get('ref_file', 'genome')
    read_length = opts.get('readlength', 50)
    coverage = opts.get('coverage', 1000)
    adaptor2 = opts.get('adaptor2', '1')
    skip_fastp = opts.get('skip_fastp', False)

    control_dict = _parse_to_dict(opts.get('control', ''))
    tags, files, pars = parse_input(control_dict)

    if opts.get('treatment'):
        for t in (opts['treatment'] if isinstance(opts['treatment'], list) else [opts['treatment']]):
            treatment_dict = _parse_to_dict(t)
            t_tags, t_files, t_pars = parse_input(treatment_dict)
            tags.extend(t_tags)
            files.extend(t_files)
            pars.extend(t_pars)

    par_str = ' '.join(pars)

    if not nomapping:
        tee.write("\nBuilding HISAT2 reference index ...\n")

        gene_seq = defaultdict(lambda: {'transcriptID': '1', 'seq': ''})
        gene_info = {}

        if reference == 'genome':
            gene_info = ref.get_gene_info(prefix, genome)
            with open("transcripts.fa") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith('>'):
                        m = re.match(r'(\w+)\.(\d+)', line[1:])
                        if m:
                            gene_id = m.group(1)
                            gene_seq[gene_id]['transcriptID'] = m.group(2)
                        else:
                            gene_id = line[1:]
                    else:
                        gene_seq[gene_id]['seq'] += line
            os.unlink("transcripts.fa")
        else:
            ref_path = _resolve_path(reference)
            src = f"{ref_path}.fa" if os.path.exists(f"{ref_path}.fa") else ref_path
            with open(src) as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith('>'):
                        gene_id = line[1:].split()[0]
                    else:
                        gene_seq[gene_id]['seq'] += line
            for gene_id in gene_seq:
                gene_info[gene_id] = len(gene_seq[gene_id]['seq'])

        # Build reference files
        with open("reference.fa", 'w') as ref_fh, open("reference.gff", 'w') as gff_fh:
            for gene_id in sorted(gene_seq.keys()):
                ref_fh.write(f">{gene_id}\n{gene_seq[gene_id]['seq']}\n")
                length = gene_info.get(gene_id, len(gene_seq[gene_id]['seq']))
                tid = gene_seq[gene_id]['transcriptID']
                gff_fh.write(f"{gene_id}\treference\tgene\t1\t{length}\t.\t+\t.\tID={gene_id}\n")
                gff_fh.write(f"{gene_id}\treference\tmRNA\t1\t{length}\t.\t+\t.\t"
                           f"ID={gene_id}.{tid};Parent={gene_id}\n")
                gff_fh.write(f"{gene_id}\treference\texon\t1\t{length}\t.\t+\t.\t"
                           f"ID={gene_id}:exon:1;Parent={gene_id}.{tid}\n")

        # HISAT2 index — task-specific (reference.fa/gff 是过滤后的基因集，
        # 不兼容 genome-wide HISAT2 index；每次 run 都重建)
        hisat_idx = "HisatIdx"
        hisat2_build_index("reference.fa", "reference.gff", hisat_idx, thread, tee=tee)

        hisat_extra = "-k 50 --mp 3,3 --score-min L,0,-0.3"

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nMapping {tag}...\n")

            if ',' not in fpath:
                sra_results, from_sra = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    unzip_file(sra_results[0], tag)
                    preserve_raw_fastq(from_sra, tag, paired=False)
                    if adaptor and not skip_fastp:
                        tee.write("Trimming...\n")
                        run_cmd(
                            f"fastp -w {thread} --length_required 15 "
                            f"-i {tag}.fastq -o {tag}_trimmed.fastq --adapter_sequence {adaptor}")
                        os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")
                    hisat2_align_single(tag, f"{tag}.fastq", hisat_idx, thread,
                                        extra_args=hisat_extra)
                    if from_sra:
                        gzip_fastq(f"{tag}.raw.fastq")
                else:
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    preserve_raw_fastq(from_sra, tag, paired=True)
                    if adaptor and not skip_fastp:
                        tee.write("Trimming...\n")
                        run_cmd(
                            f"fastp -w {thread} --length_required 15 "
                            f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                            f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq --adapter_sequence {adaptor} --adapter_sequence_r2 {adaptor2}")
                        os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                        os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                    hisat2_align_paired(tag, f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                                        hisat_idx, thread,
                                        extra_args=hisat_extra)
                    if from_sra:
                        for _f in (f"{tag}_R1.raw.fastq", f"{tag}_R2.raw.fastq"):
                            gzip_fastq(_f)
            else:
                f1, f2 = fpath.split(',')
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")
                if adaptor and not skip_fastp:
                    run_cmd(
                        f"fastp -w {thread} --length_required 15 "
                        f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                        f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq --adapter_sequence {adaptor} --adapter_sequence_r2 {adaptor2}")
                    os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                    os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                hisat2_align_paired(tag, f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                                    hisat_idx, thread,
                                    extra_args=hisat_extra)

            # HISAT2 直接输出 {tag}.bam — 转成 SAM 供后续覆盖度计算
            run_cmd(f"samtools view -h {tag}.bam > {tag}.sam")
            if os.path.exists(f"{tag}_Log.txt"):
                with open(f"{tag}_Log.txt") as lf:
                    tee.write(lf.read())

            # Count reads per gene
            output_gene = defaultdict(int)
            with open(f"{tag}.sam") as fh:
                for line in fh:
                    if line.startswith('@'):
                        continue
                    cols = line.strip().split('\t')
                    if len(cols) >= 3 and cols[2] != '*':
                        output_gene[cols[2]] += 1

            # Filter genes by coverage
            with open(f"{tag}.fa", 'w') as ref_out, open(f"{tag}.gff", 'w') as gff_out:
                for gene_id in sorted(output_gene.keys()):
                    min_reads = coverage * len(gene_seq.get(gene_id, {}).get('seq', '')) / 2
                    if output_gene[gene_id] >= min_reads:
                        length = gene_info.get(gene_id, len(gene_seq.get(gene_id, {}).get('seq', '')))
                        tid = gene_seq.get(gene_id, {}).get('transcriptID', '1')
                        ref_out.write(f">{gene_id}\n{gene_seq.get(gene_id, {}).get('seq', '')}\n")
                        gff_out.write(f"{gene_id}\treference\tgene\t1\t{length}\t.\t+\t.\tID={gene_id}\n")
                        gff_out.write(f"{gene_id}\treference\tmRNA\t1\t{length}\t.\t+\t.\t"
                                     f"ID={gene_id}.{tid};Parent={gene_id}\n")
                        gff_out.write(f"{gene_id}\treference\texon\t1\t{length}\t.\t+\t.\t"
                                     f"ID={gene_id}:exon:1;Parent={gene_id}.{tid}\n")

            # Calculate read ends
            ends = defaultdict(lambda: defaultdict(int))
            with open(f"{tag}.sam") as fh_in, open(f"{tag}.filtered.sam", 'w') as fh_out:
                for line in fh_in:
                    if line.startswith('@'):
                        fh_out.write(line)
                        continue
                    cols = line.strip().split('\t')
                    if len(cols) < 6:
                        continue
                    ref_name = cols[2]
                    if output_gene.get(ref_name, 0) < coverage * len(gene_seq.get(ref_name, {}).get('seq', '')) / 2:
                        continue

                    cigar = cols[5]
                    m = re.match(r'^(\d+)M$', cigar)
                    if not m:
                        continue
                    map_len = int(m.group(1))
                    flag = cols[1]
                    pos = int(cols[3])

                    if int(flag) & 16 == 0:
                        ends[ref_name][pos] += 1
                        if map_len < read_length:
                            ends[ref_name][pos + map_len] += 1
                    else:
                        ends[ref_name][pos + map_len] += 1
                        if map_len < read_length:
                            ends[ref_name][max(pos - 1, 1)] += 1

                    fh_out.write(line)

            # Write wig file
            with open(f"{tag}.ends.wig", 'w') as fh:
                for chr_name in sorted(ends.keys()):
                    fh.write(f"variableStep chrom={chr_name}\n")
                    for site in sorted(ends[chr_name].keys()):
                        fh.write(f"{site}\t{ends[chr_name][site]}\n")

            run_cmd(f"samtools view -Sb {tag}.filtered.sam > {tag}.filtered.bam")
            run_cmd(f"samtools index {tag}.filtered.bam")

            for fname in (f"{tag}.bam", f"{tag}.bam.bai", f"{tag}.sam", f"{tag}.filtered.sam"):
                if os.path.exists(fname):
                    os.unlink(fname)

            tee.write("\nCalculating RiboMeth Scores...\n")
            run_cmd(
                f"Rscript --vanilla {prefix}/scripts/RNAmodR.R {tag} {coverage}")

        # Cleanup
        for fname in ("reference.fa", "reference.gff"):
            if os.path.exists(fname):
                os.unlink(fname)
        cleanup_alignment_outputs(extra_globs=(f"{genome}_chr_all*", f"{genome}*.bt2", "reference.fa*", "reference.gff*",
                                                 "HisatIdx*.ht2", "HisatIdx.htl", "_HisatIdx_splicesites.txt"))
