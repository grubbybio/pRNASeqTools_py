"""
ATAC-seq analysis mode.
bowtie2 alignment → Genrich / MACS3 peak calling (+ MACS3 bdgdiff).
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input, _parse_to_dict
from prnaseqtools.functions import download_sra, unzip_file, _tee, run_cmd, gzip_fastq, preserve_raw_fastq


def run(opts):
    """Main entry point for ATAC-seq analysis."""
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
    cutoff = opts.get('cutoff', 3.0)
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)
    tss_distance = opts.get('tss_distance', 3000)
    skip_fastp   = opts.get('skip_fastp', False)
    skip_insert_size = opts.get('skip_insert_size', False)
    skip_peak_metaplot = opts.get('skip_peak_metaplot', False)
    metaplot_upstream = opts.get('metaplot_upstream', 2000)
    metaplot_downstream = opts.get('metaplot_downstream', 2000)

    tags = []
    files = []

    input_opt = opts.get('control')
    if isinstance(input_opt, list):
        input_opt = input_opt[0] if input_opt else ''
    genrich_input = ""
    if input_opt:
        input_dict = _parse_to_dict(input_opt)
        i_tags, i_files, _ = parse_input(input_dict)
        tags.extend(i_tags)
        files.extend(i_files)
        genrich_input = "-c " + ','.join(f"{t}.sorted.name.bam" for t in i_tags)

    # Parse IP (treatment) — pool all --treatment entries as one IP group
    ip_opt = opts.get('treatment', '')
    if isinstance(ip_opt, str):
        ip_opt = [ip_opt] if ip_opt else []
    ip_tags = []
    for ip_str in ip_opt:
        ip_dict = _parse_to_dict(ip_str)
        i_tags, i_files, _ = parse_input(ip_dict)
        ip_tags.extend(i_tags)
        tags.extend(i_tags)
        files.extend(i_files)
    genrich_ip = "-t " + ','.join(f"{t}.sorted.name.bam" for t in ip_tags) if ip_tags else ""
    ip_tag = ip_tags[0] if ip_tags else 'ip'

    if not ip_tags:
        sys.exit("At least one --treatment (-p) is required.")

    # Parse group2 for MACS3 bdgdiff
    ip2_tags = []
    if peak_caller == 'macs3':
        treatment2_opt = opts.get('treatment2')
        if treatment2_opt:
            t2_dict = _parse_to_dict(treatment2_opt)
            t2_tags, t2_files, _ = parse_input(t2_dict)
            ip2_tags.extend(t2_tags)
            tags.extend(t2_tags)
            files.extend(t2_files)

    if not nomapping:
        tee.write("\nBuilding index...\n")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        if not os.path.exists(fasta_path):
            sys.exit(f"Reference fasta not found: {fasta_path}")
        if not os.path.lexists(f"{genome}_chr_all.fasta"):
            os.symlink(fasta_path, f"{genome}_chr_all.fasta")
        run_cmd(
            f"bowtie2-build --threads {thread} -q {genome}_chr_all.fasta {genome}_chr_all")

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nMapping {tag}...\n")

            if ',' not in fpath:
                sra_results, from_sra = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    unzip_file(sra_results[0], tag)
                    preserve_raw_fastq(from_sra, tag, paired=False)
                    # fastp trimming (skippable via --skip-fastp)
                    if skip_fastp:
                        tee.write("  [SKIP] fastp (--skip-fastp)\n")
                    else:
                        tee.write("  Trimming (fastp)...\n")
                        fastp_adaptor = f"--adapter_sequence {adaptor}" if adaptor else ""
                        run_cmd(
                            f"fastp -w {thread} --length_required 50 "
                            f"-i {tag}.fastq -o {tag}_trimmed.fastq "
                            f"-h {tag}.html -j {tag}.json {fastp_adaptor}")
                        os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")
                        for rep in (f"{tag}.html", f"{tag}.json"):
                            if os.path.exists(rep):
                                os.unlink(rep)
                    run_cmd(
                        f"bowtie2 -p {thread} -x {genome}_chr_all -U {tag}.fastq "
                        f"-S {tag}.sam")
                    if from_sra:
                        gzip_fastq(f"{tag}.raw.fastq")
                else:
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    preserve_raw_fastq(from_sra, tag, paired=True)
                    # fastp trimming (skippable via --skip-fastp)
                    if skip_fastp:
                        tee.write("  [SKIP] fastp (--skip-fastp)\n")
                    else:
                        fastp_adaptor = f"--adapter_sequence {adaptor} --adapter_sequence_r2 {adaptor}" if adaptor else ""
                        run_cmd(
                            f"fastp -w {thread} --length_required 50 "
                            f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                            f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq "
                        f"-h {tag}.html -j {tag}.json {fastp_adaptor}")
                        os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                        os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                        for rep in (f"{tag}.html", f"{tag}.json"):
                            if os.path.exists(rep):
                                os.unlink(rep)
                    run_cmd(
                        f"bowtie2 -p {thread} -x {genome}_chr_all "
                        f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq -S {tag}.sam")
                    if from_sra:
                        for _f in (f"{tag}_R1.raw.fastq", f"{tag}_R2.raw.fastq"):
                            gzip_fastq(_f)
            else:
                f1, f2 = fpath.split(',')
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")
                # fastp trimming (skippable via --skip-fastp)
                if skip_fastp:
                    tee.write("  [SKIP] fastp (--skip-fastp)\n")
                else:
                    fastp_adaptor = f"--adapter_sequence {adaptor} --adapter_sequence_r2 {adaptor}" if adaptor else ""
                    run_cmd(
                        f"fastp -w {thread} --length_required 50 "
                        f"-i {tag}_R1.fastq -I {tag}_R2.fastq "
                        f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq "
                        f"-h {tag}.html -j {tag}.json {fastp_adaptor}")
                    os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                    os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                    for rep in (f"{tag}.html", f"{tag}.json"):
                        if os.path.exists(rep):
                            os.unlink(rep)
                run_cmd(
                    f"bowtie2 -p {thread} -x {genome}_chr_all "
                    f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq -S {tag}.sam")

            run_cmd(
                f"samtools view -Sb -q 10 --threads {thread} {tag}.sam > {tag}.bam")
            os.unlink(f"{tag}.sam")
            run_cmd(f"samtools sort -n -o {tag}.sorted.name.bam {tag}.bam")
            run_cmd(f"samtools sort -o {tag}.sorted.bam {tag}.bam")
            run_cmd(f"samtools index {tag}.sorted.bam")
            os.unlink(f"{tag}.bam")
            run_cmd(
                f"bamCoverage -b {tag}.sorted.bam -bs 5 -p {thread} "
                f"--Offset 1 --smoothLength 150 "
                f"--normalizeUsing CPM -o {tag}.bw")
            tee.write("\nMapping completed!\n")

        for fname in globmod.glob(f"{genome}_chr_all*"):
            os.unlink(fname)
        if os.path.exists("igv.log"):
            os.unlink("igv.log")

        if not mappingonly:
            # ── Insert size distribution ───────────────────────────────────
            if not skip_insert_size:
                _run_insert_size_plot(ip_tags, tee, thread, out_prefix=ip_tag)

            if peak_caller == 'macs3':
                _run_macs3(ip_tags, ip2_tags, ip_tag, genome_size,
                           qvalue, pvalue, tee, cutoff=cutoff,
                           genome=genome, prefix=prefix,
                           tss_distance=tss_distance, thread=thread,
                           skip_peak_metaplot=skip_peak_metaplot,
                           metaplot_upstream=metaplot_upstream,
                           metaplot_downstream=metaplot_downstream)
            else:
                _run_genrich(genrich_ip, genrich_input, ip_tag,
                             qvalue, pvalue, auc, tee,
                             genome=genome, prefix=prefix,
                             tss_distance=tss_distance, thread=thread,
                             skip_peak_metaplot=skip_peak_metaplot,
                             metaplot_upstream=metaplot_upstream,
                             metaplot_downstream=metaplot_downstream)

    else:
        for pre in tags:
            if not os.path.lexists(f"{pre}.sorted.name.bam"):
                os.symlink(f"../{pre}.sorted.name.bam", f"{pre}.sorted.name.bam")
            if not os.path.lexists(f"{pre}.sorted.bam"):
                os.symlink(f"../{pre}.sorted.bam", f"{pre}.sorted.bam")
            if (os.path.exists(f"../{pre}.sorted.bam.bai")
                    and not os.path.lexists(f"{pre}.sorted.bam.bai")):
                os.symlink(f"../{pre}.sorted.bam.bai", f"{pre}.sorted.bam.bai")
            # bigWig tracks are required for peak meta-plot/heatmap
            if (os.path.exists(f"../{pre}.bw")
                    and not os.path.lexists(f"{pre}.bw")):
                os.symlink(f"../{pre}.bw", f"{pre}.bw")

        # ── Insert size distribution ───────────────────────────────────────
        if not skip_insert_size:
            _run_insert_size_plot(ip_tags, tee, thread, out_prefix=ip_tag)

        if peak_caller == 'macs3':
            _run_macs3(ip_tags, ip2_tags, ip_tag, genome_size,
                       qvalue, pvalue, tee, cutoff=cutoff,
                       genome=genome, prefix=prefix,
                       tss_distance=tss_distance, thread=thread,
                       skip_peak_metaplot=skip_peak_metaplot,
                       metaplot_upstream=metaplot_upstream,
                       metaplot_downstream=metaplot_downstream)
        else:
            _run_genrich(genrich_ip, genrich_input, ip_tag, qvalue, pvalue, auc, tee,
                         genome=genome, prefix=prefix, tss_distance=tss_distance,
                         thread=thread,
                         skip_peak_metaplot=skip_peak_metaplot,
                         metaplot_upstream=metaplot_upstream,
                         metaplot_downstream=metaplot_downstream)

        for pre in tags:
            for bam in (f"{pre}.sorted.name.bam", f"{pre}.sorted.bam",
                        f"{pre}.sorted.bam.bai", f"{pre}.bw"):
                if os.path.lexists(bam):
                    os.unlink(bam)


def _peak_stats(peak_file, label, tee):
    """Print basic peak count and length stats for a peak file."""
    if not os.path.exists(peak_file):
        return
    count = 0
    total_bp = 0
    with open(peak_file) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 3:
                s, e = int(cols[1]), int(cols[2])
                total_bp += e - s
                count += 1
    if count > 0:
        tee.write(f"  {label}: {count:,} peaks, avg {total_bp/count:.0f} bp\n")


def _run_insert_size_plot(tags, tee, thread, out_prefix="atac"):
    """Insert/fragment size distribution using deepTools bamPEFragmentSize.

    For paired-end ATAC-seq this reports true fragment (insert) sizes and
    highlights the nucleosome phasing pattern (mono-/di-/tri-nucleosome).
    Single-end libraries produce a read-length distribution instead and may
    be skipped by the tool; in that case a warning is logged and execution
    continues.
    """
    bam_list = [f"{t}.sorted.bam" for t in tags if os.path.exists(f"{t}.sorted.bam")]
    if not bam_list:
        tee.write("  [SKIP] No sorted BAM files found; skipping insert size plot.\n")
        return

    labels = ' '.join(t for t in tags if os.path.exists(f"{t}.sorted.bam"))
    out_png = f"{out_prefix}_insert_size.png"
    out_tsv = f"{out_prefix}_insert_size.tsv"

    cmd = (
        f"bamPEFragmentSize -b {' '.join(bam_list)} "
        f"--samplesLabel {labels} "
        f"--histogram {out_png} "
        f"--outRawFragmentLengths {out_tsv} "
        f"--numberOfProcessors {thread} "
        f"--binSize 5 "
        f"--plotTitle 'ATAC-seq insert size distribution'"
    )
    tee.write("\nGenerating insert size distribution plot (bamPEFragmentSize)...\n")
    try:
        run_cmd(cmd, check=True)
        tee.write(f"  Insert size plot : {out_png}\n")
        tee.write(f"  Raw lengths      : {out_tsv}\n")
    except subprocess.CalledProcessError:
        tee.write("  [WARN] bamPEFragmentSize failed (single-end data or tool error); "
                  "skipping insert size plot.\n")


def _run_peak_metaplot(tags, peak_file, peak_label, tee, thread,
                       upstream=2000, downstream=2000):
    """Heatmap (and default profile) of bigWig signal over peak centres.

    Uses computeMatrix reference-point (peak centre) followed by
    plotHeatmap with default settings. Every tag in ``tags`` that has a
    corresponding ``{tag}.bw`` file contributes one track, so differential
    groups can be compared on the same peak set.
    """
    if not os.path.exists(peak_file):
        tee.write(f"  [SKIP] Peak file {peak_file} not found; skipping meta-plot.\n")
        return

    bw_list = [f"{t}.bw" for t in tags if os.path.exists(f"{t}.bw")]
    if not bw_list:
        tee.write("  [SKIP] No bigWig files found; skipping meta-plot & heatmap.\n")
        return

    labels = ' '.join(t for t in tags if os.path.exists(f"{t}.bw"))
    matrix_file = f"{peak_label}_peak_matrix.gz"
    heatmap_png = f"{peak_label}_peak_heatmap.png"

    tee.write(f"\nComputing signal matrix over {peak_label} peaks "
              f"({upstream} bp upstream / {downstream} bp downstream)...\n")
    cmd_matrix = (
        f"computeMatrix reference-point "
        f"--referencePoint center "
        f"-b {upstream} -a {downstream} "
        f"-R {peak_file} "
        f"-S {' '.join(bw_list)} "
        f"--samplesLabel {labels} "
        f"-o {matrix_file} "
        f"-p {thread} "
        f"--skipZeros"
    )
    try:
        run_cmd(cmd_matrix, check=True)
    except subprocess.CalledProcessError:
        tee.write("  [WARN] computeMatrix failed; skipping meta-plot & heatmap.\n")
        return

    tee.write("  Plotting heatmap (plotHeatmap)...\n")
    run_cmd(
        f"plotHeatmap -m {matrix_file} "
        f"-out {heatmap_png} "
        f"--plotTitle 'Signal over {peak_label} peaks' "
        f"--yAxisLabel 'Normalized signal (CPM)' "
        f"--colorMap Blues",
        check=False,
    )

    tee.write(f"  Heatmap : {heatmap_png}\n")
    tee.write(f"  Matrix  : {matrix_file}\n")


def _run_macs3(ip_tags, ip2_tags, ip_tag, genome_size, qvalue, pvalue, tee,
               cutoff=3.0, genome=None, prefix=None, tss_distance=3000,
               thread=4, skip_peak_metaplot=False,
               metaplot_upstream=2000, metaplot_downstream=2000):
    """Run MACS3 peak calling for ATAC-seq (+ optional bdgdiff)."""
    ip_bam = ' '.join(f"{t}.sorted.bam" for t in ip_tags)

    # ATAC-specific parameters: no model, shift, extsize
    tee.write(f"\nMACS3 callpeak — Group 1: {ip_tag}\n")
    cmd1 = (
        f"macs3 callpeak -t {ip_bam} "
        f"-f BAM -g {genome_size} -n {ip_tag} "
        f"--nomodel --extsize 200 --shift -100"
    )
    if qvalue < 1:
        tee.write(f"  Q-value threshold: {qvalue}\n")
        cmd1 += f" -q {qvalue}"
    else:
        tee.write(f"  P-value threshold: {pvalue}\n")
        cmd1 += f" -p {pvalue}"
    cmd1 += " --bdg"
    run_cmd(cmd1)

    # Group 2 (if provided)
    if ip2_tags:
        ip2_bam = ' '.join(f"{t}.sorted.bam" for t in ip2_tags)
        ip2_tag = ip2_tags[0]

        tee.write(f"\nMACS3 callpeak — Group 2: {ip2_tag}\n")
        cmd2 = (
            f"macs3 callpeak -t {ip2_bam} "
            f"-f BAM -g {genome_size} -n {ip2_tag} "
            f"--nomodel --extsize 200 --shift -100"
        )
        if qvalue < 1:
            cmd2 += f" -q {qvalue}"
        else:
            cmd2 += f" -p {pvalue}"
        cmd2 += " --bdg"
        run_cmd(cmd2)

        # bdgdiff for ATAC (treatment-only, no controls)
        tee.write(f"\nMACS3 bdgdiff — {ip_tag} vs {ip2_tag}\n")
        t1_bdg = f"{ip_tag}_treat_pileup.bdg"
        c1_bdg = f"{ip_tag}_control_lambda.bdg"
        t2_bdg = f"{ip2_tag}_treat_pileup.bdg"
        c2_bdg = f"{ip2_tag}_control_lambda.bdg"

        diff_prefix = f"diff_{ip_tag}_vs_{ip2_tag}"
        diff_cmd = (
            f"macs3 bdgdiff "
            f"--t1 {t1_bdg} --c1 {c1_bdg} "
            f"--t2 {t2_bdg} --c2 {c2_bdg} "
            f"--o-prefix {diff_prefix}"
        )
        tee.write(f"  Log-likelihood ratio cutoff: {cutoff}\n")
        diff_cmd += f" -C {cutoff}"
        run_cmd(diff_cmd)

        tee.write(f"\nDifferential peaks output:\n")
        tee.write(f"  {diff_prefix}_cond1.bed  (enriched in group1)\n")
        tee.write(f"  {diff_prefix}_cond2.bed  (enriched in group2)\n")
        tee.write(f"  {diff_prefix}_common.bed (common peaks)\n")

        # ── Peak stats ────────────────────────────────────────────────────
        _peak_stats(f"{ip_tag}_peaks.narrowPeak", ip_tag, tee)
        _peak_stats(f"{ip2_tag}_peaks.narrowPeak", ip2_tag, tee)
        _peak_stats(f"{diff_prefix}_cond1.bed", f"diff_{ip_tag}_enriched", tee)
        _peak_stats(f"{diff_prefix}_cond2.bed", f"diff_{ip2_tag}_enriched", tee)
        _peak_stats(f"{diff_prefix}_common.bed", "diff_common", tee)

        # ── ChIPseeker annotation ─────────────────────────────────────────
        diff_peaks = [
            (f"{ip_tag}_peaks.narrowPeak",   f"{ip_tag}_peaks"),
            (f"{ip2_tag}_peaks.narrowPeak",  f"{ip2_tag}_peaks"),
            (f"{diff_prefix}_cond1.bed",     f"diff_{ip_tag}_enriched"),
            (f"{diff_prefix}_cond2.bed",     f"diff_{ip2_tag}_enriched"),
            (f"{diff_prefix}_common.bed",    "diff_common"),
        ]
        for pf, pn in diff_peaks:
            if os.path.exists(pf):
                _run_chipseeker(genome, prefix, pf, pn, tee, tss_distance)

        # ── Peak meta-plot & heatmap ──────────────────────────────────────
        if not skip_peak_metaplot:
            all_tags = ip_tags + ip2_tags
            # Each group's own peaks → that group's bigWigs
            _run_peak_metaplot(
                ip_tags, f"{ip_tag}_peaks.narrowPeak", f"{ip_tag}_peaks", tee,
                thread=thread, upstream=metaplot_upstream,
                downstream=metaplot_downstream)
            _run_peak_metaplot(
                ip2_tags, f"{ip2_tag}_peaks.narrowPeak", f"{ip2_tag}_peaks", tee,
                thread=thread, upstream=metaplot_upstream,
                downstream=metaplot_downstream)
            # Differential / common peaks → all tracks for cross-group comparison
            _run_peak_metaplot(
                all_tags, f"{diff_prefix}_cond1.bed", f"diff_{ip_tag}_enriched", tee,
                thread=thread, upstream=metaplot_upstream,
                downstream=metaplot_downstream)
            _run_peak_metaplot(
                all_tags, f"{diff_prefix}_cond2.bed", f"diff_{ip2_tag}_enriched", tee,
                thread=thread, upstream=metaplot_upstream,
                downstream=metaplot_downstream)
            _run_peak_metaplot(
                all_tags, f"{diff_prefix}_common.bed", "diff_common", tee,
                thread=thread, upstream=metaplot_upstream,
                downstream=metaplot_downstream)

        # Cleanup bedGraph files
        for pat in ["*_treat_pileup.bdg", "*_control_lambda.bdg"]:
            for fname in globmod.glob(pat):
                os.unlink(fname)
    else:
        tee.write(f"\nMACS3 peak calling complete.\n")
        tee.write(f"Output: {ip_tag}_peaks.narrowPeak\n")
        # ── Peak stats ────────────────────────────────────────────────────
        _peak_stats(f"{ip_tag}_peaks.narrowPeak", ip_tag, tee)
        # ── ChIPseeker annotation ─────────────────────────────────────────
        _run_chipseeker(genome, prefix, f"{ip_tag}_peaks.narrowPeak",
                        f"{ip_tag}_peaks", tee, tss_distance)
        # ── Peak meta-plot & heatmap ──────────────────────────────────────
        if not skip_peak_metaplot:
            _run_peak_metaplot(
                ip_tags, f"{ip_tag}_peaks.narrowPeak", f"{ip_tag}_peaks", tee,
                thread=thread, upstream=metaplot_upstream,
                downstream=metaplot_downstream)
        # Cleanup bedGraph files
        for pat in ["*_treat_pileup.bdg", "*_control_lambda.bdg"]:
            for fname in globmod.glob(pat):
                os.unlink(fname)


def _run_genrich(genrich_ip, genrich_input, ip_tag, qvalue, pvalue, auc, tee,
                 genome=None, prefix=None, tss_distance=3000, thread=4,
                 skip_peak_metaplot=False, metaplot_upstream=2000,
                 metaplot_downstream=2000):
    """Run Genrich for ATAC peak calling."""
    command = (
        f"Genrich -j -y -r -v -a {auc} -e chrC,chrM {genrich_ip} {genrich_input} "
        f"-o {ip_tag}.ATAC.narrowPeak.txt"
    )

    if qvalue < 1:
        tee.write(f"\nFinding peaks...\nAUC\t{auc}\tQ Value\t{qvalue}\n")
        command += f" -q {qvalue}"
    else:
        tee.write(f"\nFinding peaks...\nAUC\t{auc}\tP Value\t{pvalue}\n")
        command += f" -p {pvalue}"

    run_cmd(command)

    # ── Peak stats ────────────────────────────────────────────────────────
    _peak_stats(f"{ip_tag}.ATAC.narrowPeak.txt", ip_tag, tee)

    # ── ChIPseeker annotation ─────────────────────────────────────────────
    _run_chipseeker(genome, prefix, f"{ip_tag}.ATAC.narrowPeak.txt",
                    f"{ip_tag}_peaks", tee, tss_distance)

    # ── Peak meta-plot & heatmap ──────────────────────────────────────────
    if not skip_peak_metaplot:
        # Reconstruct treatment tags from the genrich -t argument
        ip_tags = [
            p.strip().replace('.sorted.name.bam', '')
            for p in genrich_ip.replace('-t ', '').split(',')
            if p.strip()
        ]
        _run_peak_metaplot(
            ip_tags, f"{ip_tag}.ATAC.narrowPeak.txt", f"{ip_tag}_peaks", tee,
            thread=thread, upstream=metaplot_upstream,
            downstream=metaplot_downstream,
        )


def _run_chipseeker(genome, prefix, peak_file, peak_name, tee, tss_distance=3000):
    """Call chipseeker.R for peak annotation and GO enrichment on a single peak set."""
    gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
    if not os.path.exists(gff_path):
        tee.write(f"  GFF not found: {gff_path}, skipping ChIPseeker annotation.\n")
        return
    chipseeker_r = os.path.join(prefix, "scripts", "chipseeker.R")
    if not os.path.exists(chipseeker_r):
        tee.write(f"  Script not found: {chipseeker_r}, skipping.\n")
        return
    run_cmd(
        f"Rscript --vanilla {chipseeker_r} "
        f"{genome} {prefix} {peak_file} {peak_name} {tss_distance}"
    )