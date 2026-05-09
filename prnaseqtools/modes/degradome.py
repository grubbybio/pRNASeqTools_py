"""
Degradome-seq (PARE/GMUCT) analysis mode.
STAR alignment to transcriptome + genome → sPARTA peak finding → CRI calculation.
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
from prnaseqtools.functions import download_sra, unzip_file, revcomp, _tee
from prnaseqtools import reference as ref


def run(opts):
    """Main entry point for degradome-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)
    targets = opts.get('targets', 'all')
    sirna = opts.get('si_rnas', 'none')

    # Parse control
    control_dict = _parse_to_dict(opts.get('control', ''))
    tags, files, _ = parse_input(control_dict)

    pars_list = []
    treatment = opts.get('treatment')
    if treatment:
        treatment_dict = _parse_to_dict(treatment)
        t_tags, t_files, t_pars = parse_input(treatment_dict)
        tags.extend(t_tags)
        files.extend(t_files)
        pars_list = t_pars

    par_str = ' '.join(tags)

    if not nomapping:
        ref.get_gene_info(prefix, genome)
        os.rename("transcripts.fa", f"{genome}.fa")

        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        subprocess.run(f"gffread -T {gff_path} -o {genome}.gtf", shell=True, check=True)

        # Build transcriptome index
        if os.path.exists("Genome"):
            subprocess.run("rm -rf Genome", shell=True)
        os.makedirs("Genome", exist_ok=True)
        subprocess.run(
            f"STAR --runThreadN {thread} --genomeDir Genome --runMode genomeGenerate "
            f"--genomeFastaFiles {genome}.fa --limitGenomeGenerateRAM 64000000000",
            shell=True, check=True
        )

        # Build genome index
        if os.path.exists("Genome2"):
            subprocess.run("rm -rf Genome2", shell=True)
        os.makedirs("Genome2", exist_ok=True)
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        subprocess.run(
            f"STAR --runThreadN {thread} --genomeDir Genome2 --runMode genomeGenerate "
            f"--genomeFastaFiles {fasta_path} --sjdbGTFfile {gff_path} "
            f"--sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene ID "
            f"--limitGenomeGenerateRAM 64000000000",
            shell=True, check=True
        )

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nWorking on {tag}...\n")

            sra_results = download_sra(fpath, thread)
            unzip_file(sra_results[0], tag)

            if adaptor:
                tee.write("\nTrimming...\n")
                subprocess.run(
                    f"cutadapt -j {thread} -m 18 --trim-n -a {adaptor} "
                    f"-o {tag}_trimmed.fastq {tag}.fastq 2>&1",
                    shell=True, check=True
                )
                os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")

            tee.write("\nStart mapping...\n")

            # Map to transcriptome
            subprocess.run(
                f"STAR --genomeDir Genome --outSAMtype BAM SortedByCoordinate "
                f"--limitBAMsortRAM 10000000000 --outSAMmultNmax 1 "
                f"--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 "
                f"--limitOutSJcollapsed 10000000 --limitIObufferSize 280000000 "
                f"--runThreadN {thread} --readFilesIn {tag}.fastq 2>&1",
                shell=True, check=True
            )
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

            # Map to genome
            subprocess.run(
                f"STAR --genomeDir Genome2 --outSAMtype BAM SortedByCoordinate "
                f"--limitBAMsortRAM 10000000000 --outSAMmultNmax 1 "
                f"--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 "
                f"--limitOutSJcollapsed 10000000 --limitIObufferSize 280000000 "
                f"--runThreadN {thread} --readFilesIn {tag}.fastq 2>&1",
                shell=True, check=True
            )
            os.rename("Aligned.sortedByCoord.out.bam", f"{tag}.genomic.bam")
            if os.path.exists("Log.final.out"):
                with open("Log.final.out") as lf:
                    tee.write(lf.read())

            subprocess.run(f"samtools index {tag}.genomic.bam", shell=True, check=True)
            subprocess.run(
                f"bamCoverage -b {tag}.genomic.bam --skipNAs -bs 5 -p {thread} "
                f"--minMappingQuality 10 --ignoreDuplicates --normalizeUsing CPM -o {tag}.genomic.bw",
                shell=True, check=True
            )

            tee.write("\nAlignment Completed!\n")

            # Create library file
            fq_demux = defaultdict(int)
            with open(f"{tag}.fastq") as fh:
                for line_num, line in enumerate(fh):
                    if line_num % 4 == 1:
                        fq_demux[line.strip()] += 1

            with open(f"{tag}_lib.txt", 'w') as fh:
                for seq, count in fq_demux.items():
                    fh.write(f"{seq}\t{count}\n")

            if os.path.exists(f"{tag}.fastq"):
                os.unlink(f"{tag}.fastq")

        if not mappingonly:
            tee.write("Finding peaks...\n")
            _find_peaks(thread, prefix, genome, sirna, tags)

            tee.write("Calculating CRIs...\n")
            subprocess.run(
                f"Rscript --vanilla {prefix}/scripts/ribo.R {genome} {par_str}",
                shell=True, check=True
            )
            _calculate_cri(tags, targets)
            subprocess.run(
                f"Rscript --vanilla {prefix}/scripts/CRI.R {par_str}",
                shell=True, check=True
            )

        # Cleanup
        for fname in globmod.glob("Log.*"):
            os.unlink(fname)
        for fname in ("SJ.out.tab", f"{genome}.fa", f"{genome}.gtf", f"{genome}.fa.fai"):
            if os.path.exists(fname):
                os.unlink(fname)
        for d in ("Genome", "Genome2"):
            if os.path.exists(d):
                subprocess.run(f"rm -rf {d}", shell=True)

    else:
        # No mapping mode
        for pre in tags:
            os.symlink(f"../{pre}.bam", f"{pre}.bam")
            os.symlink(f"../{pre}_lib.txt", f"{pre}_lib.txt")

        tee.write("Finding peaks...\n")
        ref.get_gene_info(prefix, genome)
        os.rename("transcripts.fa", f"{genome}.fa")
        _find_peaks(thread, prefix, genome, sirna, tags)

        tee.write("Calculating CRIs...\n")
        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        subprocess.run(f"gffread -T {gff_path} -o {genome}.gtf", shell=True, check=True)
        subprocess.run(
            f"Rscript --vanilla {prefix}/scripts/ribo.R {genome} {par_str}",
            shell=True, check=True
        )
        _calculate_cri(tags, targets)
        subprocess.run(
            f"Rscript --vanilla {prefix}/scripts/CRI.R {par_str}",
            shell=True, check=True
        )

        for fname in globmod.glob("*_lib.txt"):
            os.unlink(fname)
        for fname in (f"{genome}.gtf", f"{genome}.fa"):
            if os.path.exists(fname):
                os.unlink(fname)
        for fname in globmod.glob("*.bam"):
            os.unlink(fname)


def _parse_to_dict(arg_str):
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _calculate_cri(tags, targets):
    """Calculate cleavage ratio index for each transcript."""
    for tag in tags:
        frame = defaultdict(lambda: defaultdict(int))
        csv_file = f"{tag}.csv"
        if not os.path.exists(csv_file):
            continue

        with open(csv_file) as fh:
            fh.readline()  # header
            for line in fh:
                cols = line.strip().split(',')
                if len(cols) >= 10 and cols[9] == 'cds':
                    tmp2 = int(cols[1]) - int(cols[5])
                    if tmp2 > 0:
                        frame[cols[0]][tmp2 % 3] += 1
                        frame[cols[0]]['total'] += 1

        with open(f"{tag}_CRI.txt", 'w') as fh:
            if targets == 'all':
                for transcript in sorted(frame.keys()):
                    if frame[transcript]['total'] >= 100:
                        f0 = frame[transcript].get(0, 0)
                        f1 = frame[transcript].get(1, 0)
                        f2 = frame[transcript].get(2, 0)
                        log_val = math.log((f1 + f2 + 1) / (2 * f0 + 1), 2)
                        fh.write(f"{transcript}\t{f0}\t{f1}\t{f2}\t{log_val}\n")
            else:
                with open(targets) as tar_fh:
                    for target in tar_fh:
                        target = target.strip()
                        if target in frame and frame[target]['total'] >= 100:
                            f0 = frame[target].get(0, 0)
                            f1 = frame[target].get(1, 0)
                            f2 = frame[target].get(2, 0)
                            log_val = math.log((f1 + f2 + 1) / (2 * f0 + 1), 2)
                            fh.write(f"{target}\t{f0}\t{f1}\t{f2}\t{log_val}\n")


def _find_peaks(thread, prefix, genome, sirna, tags):
    """Run sPARTA for peak finding."""
    os.makedirs("sparta", exist_ok=True)
    os.chdir("sparta")

    os.symlink(f"../{genome}.fa", f"{genome}.fa")

    # Create miRNA FASTA
    fas = ref.read_fasta(prefix, genome)
    mirna_data = defaultdict(lambda: {'name': ''})
    mir_gff = os.path.join(prefix, "reference", f"{genome}_miRNA_miRNA_star.gff")

    with open(mir_gff) as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            miseq = fas.get(cols[0], "")[int(cols[3]) - 1:int(cols[4])]
            if cols[6] == '-':
                miseq = revcomp(miseq)
            mirna_data[miseq]['name'] += cols[8] + ';'

    with open(f"{genome}_miRNA.fa", 'w') as fh:
        for miseq in sorted(mirna_data.keys()):
            mirna_data[miseq]['name'] = mirna_data[miseq]['name'].rstrip(';')
            fh.write(f">{mirna_data[miseq]['name']}\n{miseq}\n")

    if sirna != 'none':
        subprocess.run(f"cat {sirna} >> {genome}_miRNA.fa", shell=True, check=True)

    # Symlink sPARTA
    sparta_path = os.path.join(prefix, "sPARTA.py")
    os.symlink(sparta_path, "sPARTA.py")

    # Symlink library files
    tags_lib = []
    for tag in tags:
        lib_file = f"{tag}_lib.txt"
        os.symlink(f"../{lib_file}", lib_file)
        tags_lib.append(lib_file)

    libs_str = ' '.join(tags_lib)
    subprocess.run(
        f"python3 sPARTA.py -accel {thread} -featureFile {genome}.fa "
        f"-genomeFeature 0 -miRNAFile {genome}_miRNA.fa -libs {libs_str} "
        f"-minTagLen 18 -tarPred -tarScore --tag2FASTA --map2DD --validate",
        shell=True, check=True
    )

    # Move output files
    if os.path.exists("output"):
        for fname in os.listdir("output"):
            src = os.path.join("output", fname)
            dst = os.path.join("..", fname)
            if os.path.exists(dst):
                os.unlink(dst)
            os.rename(src, dst)
        os.rmdir("output")

    os.chdir("..")
    subprocess.run("rm -rf sparta", shell=True)
