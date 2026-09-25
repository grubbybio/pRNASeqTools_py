"""
Degradome-seq (PARE/GMUCT) analysis mode.
HISAT2 alignment to transcriptome + genome → sPARTA peak finding → CRI calculation.
"""

import os
import sys
import glob as globmod
import subprocess
import math
from pathlib import Path
from collections import defaultdict

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input, _parse_to_dict
from prnaseqtools.functions import (download_sra, unzip_file, revcomp, _tee, run_cmd,
                                     gzip_fastq, preserve_raw_fastq,
                                     cleanup_alignment_outputs,
                                     hisat2_build_index, hisat2_align_single,
                                     try_use_shared_hisat_index, save_hisat_index_to_reference)
from prnaseqtools import reference as ref


def run(opts):
    """Main entry point for degradome-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    # prefix must point to the project root so that prefix/reference/{genome}_Genome
    # resolves to ~/software/pRNASeqTools_py/reference/ (not prnaseqtools/reference/)
    prefix = opts.get('prefix', '/Users/cjyou/software/pRNASeqTools_py')
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)
    targets = opts.get('targets', 'all')
    sirna = opts.get('si_rnas', 'none')

    # Parse control
    control_dict = _parse_to_dict(opts.get('control', ''))
    tags, files, pars = parse_input(control_dict)

    treatment = opts.get('treatment')
    if treatment:
        for t in (treatment if isinstance(treatment, list) else [treatment]):
            treatment_dict = _parse_to_dict(t)
            t_tags, t_files, t_pars = parse_input(treatment_dict)
            tags.extend(t_tags)
            files.extend(t_files)
            pars.extend(t_pars)

    tag_str = ' '.join(tags)

    if not nomapping:
        ref.get_gene_info(prefix, genome)
        os.rename("transcripts.fa", f"{genome}.fa")

        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        run_cmd(f"gffread -T {gff_path} -o {genome}.gtf")
        ref.primary_transcript(prefix, genome)

        # Build transcriptome HISAT2 index (special — uses {genome}.fa, not the standard
        # {genome}_chr_all.fasta; task-local 不共享，因为跟 genome-wide 不同)
        tx_idx = "TxIdx"
        tee.write("\nBuilding transcriptome HISAT2 index ...\n")
        hisat2_build_index(f"{genome}.fa", f"{genome}.gtf", tx_idx, thread, tee=tee)

        # Build/use genome HISAT2 index (共享 via reference/{genome}_HisatIdx)
        gn_idx = "GnIdx"
        tee.write("Building/using genome HISAT2 index ...\n")
        if not try_use_shared_hisat_index(genome, prefix, gn_idx, tee):
            fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
            hisat2_build_index(fasta_path, gff_path, gn_idx, thread, tee=tee)
            save_hisat_index_to_reference(genome, prefix, gn_idx, tee)

        hisat_extra = "--max-hits 50 --mp 3,3 --score-min L,0,-0.3"

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nWorking on {tag}...\n")

            sra_results, from_sra = download_sra(fpath, thread)
            unzip_file(sra_results[0], tag)
            preserve_raw_fastq(from_sra, tag, paired=False)

            # fastp trimming — always run (auto-detects adapter if --adaptor not provided)
            tee.write("\nTrimming (fastp)...\n")
            fastp_adaptor = f"--adapter_sequence {adaptor}" if adaptor else ""
            run_cmd(
                f"fastp -w {thread} --length_required 18 "
                f"-i {tag}.fastq -o {tag}_trimmed.fastq {fastp_adaptor}")
            os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")
            for rep in (f"{tag}.html", f"{tag}.json"):
                if os.path.exists(rep):
                    os.unlink(rep)

            tee.write("\nStart mapping...\n")

            # Map to transcriptome
            hisat2_align_single(tag, f"{tag}.fastq", tx_idx, thread,
                                extra_args=hisat_extra)
            if os.path.exists(f"{tag}_Log.txt"):
                with open(f"{tag}_Log.txt") as lf:
                    tee.write(lf.read())

            run_cmd(f"samtools index {tag}.bam")
            run_cmd(
                f"bamCoverage -b {tag}.bam --skipNAs -bs 5 -p {thread} "
                f"--minMappingQuality 10 --ignoreDuplicates --normalizeUsing CPM -o {tag}.bw")

            # Map to genome — 用 tag={tag}.genomic 让输出自动变成 {tag}.genomic.bam
            hisat2_align_single(f"{tag}.genomic", f"{tag}.fastq", gn_idx, thread,
                                extra_args=hisat_extra)
            if os.path.exists(f"{tag}.genomic_Log.txt"):
                with open(f"{tag}.genomic_Log.txt") as lf:
                    tee.write(lf.read())

            run_cmd(f"samtools index {tag}.genomic.bam")
            run_cmd(
                f"bamCoverage -b {tag}.genomic.bam --skipNAs -bs 5 -p {thread} "
                f"--minMappingQuality 10 --ignoreDuplicates --normalizeUsing CPM -o {tag}.genomic.bw")

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

            if from_sra:
                gzip_fastq(f"{tag}.raw.fastq")

        if not mappingonly:
            tee.write("Finding peaks...\n")
            _find_peaks(thread, prefix, genome, sirna, tags, tee)

            tee.write("Calculating CRIs...\n")
            os.makedirs("genomic_bam", exist_ok=True)
            for fname in globmod.glob("*.genomic.bam*"):
                os.rename(fname, os.path.join("genomic_bam", fname))
            run_cmd(
                f"Rscript --vanilla {prefix}/scripts/ribo.R {genome} {tag_str}")
            _calculate_cri(tags, targets)
            run_cmd(
                f"Rscript --vanilla {prefix}/scripts/CRI.R {tag_str}")

        # Cleanup
        # 注：GnIdx 可能是 symlink → reference/{genome}_HisatIdx（共享 index），cleanup_alignment_outputs
        # 不会碰 *.ht2 / *.htl 文件，这里额外清 TxIdx（task-local，每次 run 都重建）
        import glob as _glob
        for pat in ("TxIdx*.ht2", "TxIdx.htl", "_TxIdx_splicesites.txt"):
            for f in _glob.glob(pat):
                os.unlink(f)
        cleanup_alignment_outputs(extra_globs=(f"{genome}.fa", f"{genome}.gtf", f"{genome}.fa.fai", f"{genome}*.fai"))
        for fname in (f"{genome}.fa", f"{genome}.gtf", f"{genome}.PrimaryTranscript.gtf", f"{genome}.fa.fai"):
            if os.path.exists(fname):
                os.unlink(fname)

    else:
        # No mapping mode
        for pre in tags:
            os.symlink(f"../{pre}.bam", f"{pre}.bam")
            os.symlink(f"../{pre}_lib.txt", f"{pre}_lib.txt")

        tee.write("Finding peaks...\n")
        ref.get_gene_info(prefix, genome)
        os.rename("transcripts.fa", f"{genome}.fa")
        _find_peaks(thread, prefix, genome, sirna, tags, tee)

        tee.write("Calculating CRIs...\n")
        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        run_cmd(f"gffread -T {gff_path} -o {genome}.gtf")
        ref.primary_transcript(prefix, genome)
        run_cmd(
            f"Rscript --vanilla {prefix}/scripts/ribo.R {genome} {tag_str}")
        _calculate_cri(tags, targets)
        run_cmd(
            f"Rscript --vanilla {prefix}/scripts/CRI.R {tag_str}")

        for fname in globmod.glob("*_lib.txt"):
            os.unlink(fname)
        for fname in (f"{genome}.gtf", f"{genome}.PrimaryTranscript.gtf", f"{genome}.fa"):
            if os.path.exists(fname):
                os.unlink(fname)
        for fname in globmod.glob("*.bam"):
            os.unlink(fname)


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


def _find_peaks(thread, prefix, genome, sirna, tags, tee):
    """Run sPARTA for peak finding."""
    sparta_path = os.path.join(prefix, "sPARTA.py")
    if not os.path.exists(sparta_path):
        tee.write(f"\nERROR: sPARTA.py not found at {sparta_path} — please install sPARTA before running degradome analysis.\n")
        sys.exit(1)

    os.makedirs("sparta", exist_ok=True)
    os.chdir("sparta")

    try:
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
            run_cmd(f"cat {sirna} >> {genome}_miRNA.fa")

        # Symlink sPARTA
        os.symlink(sparta_path, "sPARTA.py")

        # Symlink library files
        tags_lib = []
        for tag in tags:
            lib_file = f"{tag}_lib.txt"
            os.symlink(f"../{lib_file}", lib_file)
            tags_lib.append(lib_file)

        libs_str = ' '.join(tags_lib)
        run_cmd(
            f"python3 sPARTA.py -accel {thread} -featureFile {genome}.fa "
            f"-genomeFeature 0 -miRNAFile {genome}_miRNA.fa -libs {libs_str} "
            f"-minTagLen 18 -tarPred -tarScore --tag2FASTA --map2DD --validate")

        # Move output files
        if os.path.exists("output"):
            for fname in os.listdir("output"):
                src = os.path.join("output", fname)
                dst = os.path.join("..", fname)
                if os.path.exists(dst):
                    os.unlink(dst)
                os.rename(src, dst)
            os.rmdir("output")
    finally:
        os.chdir("..")
        run_cmd("rm -rf sparta")
