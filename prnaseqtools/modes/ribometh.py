"""
RiboMeth-seq analysis mode.
STAR alignment to transcriptome → RNAmodR analysis.
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path
from collections import defaultdict

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input
from prnaseqtools.functions import download_sra, unzip_file, _tee
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
        tee.write("\nBuilding STAR reference index ...\n")

        gene_seq = defaultdict(lambda: {'transcriptID': '1', 'seq': ''})
        gene_info = {}

        if reference == 'genome':
            gene_info = ref.get_gene_info(prefix, genome)
            with open("transcripts.fa") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith('>'):
                        m = __import__('re').match(r'(\w+)\.(\d+)', line[1:])
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
            os.symlink(f"{ref_path}.fa" if os.path.exists(f"{ref_path}.fa") else ref_path, "reference.fa")
            with open("reference.fa") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith('>'):
                        gene_id = line[1:]
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

        # STAR index
        if os.path.exists("Genome"):
            subprocess.run("rm -rf Genome", shell=True)
        os.makedirs("Genome", exist_ok=True)

        subprocess.run(
            f"STAR --runThreadN {thread} --genomeDir Genome --runMode genomeGenerate "
            f"--genomeFastaFiles reference.fa --sjdbGTFfile reference.gff "
            f"--sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene ID "
            f"--limitGenomeGenerateRAM 64000000000 --genomeSAindexNbases 5",
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
                        tee.write("Trimming...\n")
                        subprocess.run(
                            f"cutadapt -j {thread} -m 15 --trim-n -a {adaptor} "
                            f"-o {tag}_trimmed.fastq {tag}.fastq 2>&1",
                            shell=True, check=True
                        )
                        os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")
                    subprocess.run(
                        f"STAR --genomeDir Genome --seedSearchStartLmax 15 "
                        f"--outSAMtype SAM SortedByCoordinate --limitBAMsortRAM 10000000000 "
                        f"--outSAMmultNmax 1 --outFilterMultimapNmax 50 "
                        f"--outFilterMismatchNoverLmax 0.1 --runThreadN {thread} "
                        f"--readFilesIn {tag}.fastq 2>&1",
                        shell=True, check=True
                    )
                    if os.path.exists(f"{tag}.fastq"):
                        os.unlink(f"{tag}.fastq")
                else:
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    if adaptor:
                        tee.write("Trimming...\n")
                        subprocess.run(
                            f"cutadapt -j {thread} -m 15 --trim-n -a {adaptor} -A {adaptor2} "
                            f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                            f"{tag}_R1.fastq {tag}_R2.fastq 2>&1",
                            shell=True, check=True
                        )
                        os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                        os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                    subprocess.run(
                        f"STAR --genomeDir Genome --seedSearchStartLmax 15 "
                        f"--outSAMtype SAM SortedByCoordinate --limitBAMsortRAM 10000000000 "
                        f"--outSAMmultNmax 1 --outFilterMultimapNmax 50 "
                        f"--outFilterMismatchNoverLmax 0.1 --runThreadN {thread} "
                        f"--readFilesIn {tag}_R1.fastq {tag}_R2.fastq 2>&1",
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
                        f"cutadapt -j {thread} -m 15 --trim-n -a {adaptor} -A {adaptor2} "
                        f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                        f"{tag}_R1.fastq {tag}_R2.fastq 2>&1",
                        shell=True, check=True
                    )
                    os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                    os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                subprocess.run(
                    f"STAR --genomeDir Genome --seedSearchStartLmax 15 "
                    f"--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 "
                    f"--outSAMmultNmax 1 --outFilterMultimapNmax 50 "
                    f"--outFilterMismatchNoverLmax 0.1 --runThreadN {thread} "
                    f"--readFilesIn {tag}_R1.fastq {tag}_R2.fastq 2>&1",
                    shell=True, check=True
                )
                for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq"):
                    if os.path.exists(fname):
                        os.unlink(fname)

            subprocess.run(f"samtools view -h Aligned.sortedByCoord.out.bam > {tag}.sam", shell=True, check=True)
            if os.path.exists("Log.final.out"):
                with open("Log.final.out") as lf:
                    tee.write(lf.read())

            # Count reads per gene
            output_gene = defaultdict(int)
            with open(f"{tag}.sam") as fh:
                for line in fh:
                    if line.startswith('@'):
                        continue
                    cols = line.strip().split('\t')
                    if len(cols) >= 3:
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
                    m = __import__('re').match(r'^(\d+)M$', cigar)
                    if not m:
                        continue
                    map_len = int(m.group(1))
                    flag = cols[1]
                    pos = int(cols[3])

                    if not flag or flag == '0':
                        ends[ref_name][pos] += 1
                        if map_len < read_length:
                            ends[ref_name][pos + map_len] += 1
                    else:
                        ends[ref_name][pos + map_len] += 1
                        if map_len < read_length:
                            ends[ref_name][pos - 1] += 1

                    fh_out.write(line)

            # Write wig file
            with open(f"{tag}.ends.wig", 'w') as fh:
                for chr_name in sorted(ends.keys()):
                    fh.write(f"variableStep chrom={chr_name}\n")
                    for site in sorted(ends[chr_name].keys()):
                        fh.write(f"{site}\t{ends[chr_name][site]}\n")

            subprocess.run(f"samtools view -Sb {tag}.filtered.sam > {tag}.filtered.bam", shell=True, check=True)
            subprocess.run(f"samtools index {tag}.filtered.bam", shell=True, check=True)

            for fname in ("Aligned.sortedByCoord.out.bam", f"{tag}.sam", f"{tag}.filtered.sam"):
                if os.path.exists(fname):
                    os.unlink(fname)

            tee.write("\nCalculating RiboMeth Scores...\n")
            subprocess.run(
                f"Rscript --vanilla {prefix}/scripts/RNAmodR.R {tag} {coverage}",
                shell=True, check=True
            )

        # Cleanup
        for fname in ("reference.fa", "reference.gff"):
            if os.path.exists(fname):
                os.unlink(fname)
        for fname in globmod.glob("*.bt2") + globmod.glob("Log.*") + ["SJ.out.tab"]:
            if os.path.exists(fname):
                os.unlink(fname)
        if os.path.exists("Genome"):
            subprocess.run("rm -rf Genome", shell=True)


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
