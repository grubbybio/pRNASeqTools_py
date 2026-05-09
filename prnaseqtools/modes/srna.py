"""
Small RNA-seq analysis mode.
Mirrors the Perl mode::srna module.
ShortStack alignment → BAM processing → counting → DESeq2 DE analysis.
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path
from collections import defaultdict

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input
from prnaseqtools.functions import download_sra, unzip_file, revcomp, _tee
from prnaseqtools import reference as ref


def run(opts):
    """Main entry point for sRNA-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    mmap = opts.get('mmap', 'u')
    control = opts.get('control')
    treatment = opts.get('treatment')
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)
    foldchange = opts.get('foldchange', 1.5)
    pvalue = opts.get('pvalue', 0.01)
    binsize = opts.get('binsize', 100)
    norm = opts.get('norm', 'rRNA,total')
    norms = norm.split(',')
    mask = opts.get('mask')
    run_mode = opts.get('run_mode', 'bulk')
    pattern = opts.get('pattern', 'NNNNNNNNCA')
    promoter_length = opts.get('promoter', 1000)
    spikein = opts.get('spike_in')

    # Parse control samples
    control_dict = _parse_to_dict(control)
    tags, files, pars = parse_input(control_dict)

    # Parse treatment samples if provided
    if treatment:
        treatment_dict = _parse_to_dict(treatment)
        t_tags, t_files, t_pars = parse_input(treatment_dict)
        tags.extend(t_tags)
        files.extend(t_files)
        pars.extend(t_pars)

    par_str = ' '.join(pars)

    if not nomapping:
        ref.split_gff(prefix, genome, promoter_length)

        # Mask setup
        if mask:
            mask_path = _resolve_path(mask)
            os.symlink(mask_path, "mask.fa")
            subprocess.run("bowtie-build -q mask.fa mask", shell=True, check=True)

        # Spike-in setup
        if spikein:
            spikein_path = _resolve_path(spikein)
            os.symlink(spikein_path, "spikein.fa")
            subprocess.run("bowtie-build -q spikein.fa spikein", shell=True, check=True)

        for i in range(len(tags)):
            tag = tags[i]
            fpath = files[i]

            tee.write(f"\nMapping {tag}...\n")

            # SRA download if needed
            sra_results = download_sra(fpath, thread)
            unzip_file(sra_results[0], tag)

            # UMI extraction for single-cell mode
            if run_mode == 'sc':
                tee.write(f"\nTrimming {tag}...\n")
                subprocess.run(
                    f"umi_tools extract -p {pattern} -I {tag}.fastq -S {tag}.fq",
                    shell=True, check=True
                )
                if adaptor:
                    subprocess.run(
                        f"cutadapt -j {thread} -m 18 -M 42 --discard-untrimmed --trim-n "
                        f"-a {adaptor} -o {tag}_trimmed.fastq {tag}.fq 2>&1",
                        shell=True, check=True
                    )

                # Deduplication
                _umi_dedup(tag)
                os.unlink(f"{tag}_trimmed.fastq")
                if os.path.exists(f"{tag}.fq"):
                    os.unlink(f"{tag}.fq")
            elif run_mode == 'bulk':
                if adaptor:
                    tee.write(f"\nTrimming {tag}...\n")
                    subprocess.run(
                        f"cutadapt -j {thread} -m 18 -M 42 --discard-untrimmed --trim-n "
                        f"-a {adaptor} -o {tag}_trimmed.fastq {tag}.fastq 2>&1",
                        shell=True, check=True
                    )
                    os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")

            # Mask filtering
            if mask:
                subprocess.run(
                    f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                    f"{tag}.fastq {tag}.mask.out 2>&1",
                    shell=True, check=True
                )
                os.rename("tmp.fastq", f"{tag}.fastq")
                if os.path.exists(f"{tag}.mask.out"):
                    os.unlink(f"{tag}.mask.out")

            # Spike-in counting
            if spikein:
                subprocess.run(
                    f"bowtie -v 0 -a -p {thread} -t spikein {tag}.fastq "
                    f"{tag}.spikein.out 2>&1",
                    shell=True, check=True
                )
                subprocess.run(
                    f"awk -F \"\\t\" 'length($5)==13 && $2==\"+\"{{print $3}}' "
                    f"{tag}.spikein.out | sort | uniq -c | awk '{{print $2\"\\t\"$1}}' > {tag}.spikein",
                    shell=True, check=True
                )

            tee.write("\nStart mapping...\n")

            # rRNA filtering
            lsu_rRNA = os.path.join(prefix, "reference", "lsu_rrna")
            subprocess.run(
                f"bowtie -v 2 -a -p {thread} -t {lsu_rRNA} "
                f"{tag}.fastq {tag}.rRNA.out 2>&1",
                shell=True, check=True
            )
            subprocess.run(
                f"awk -F \"\\t\" 'BEGIN{{x=0;y=0;z=0}}"
                f"{{if($2==\"+\"){{if(/{genome}_LSU/){{x++}};"
                f"if(/{genome}_SSU/){{y++}};if(/{genome}_U6/){{z++}}}}}}"
                f"END{{print \"rRNA\\t\"x\"\\nSSU\\t\"y\"\\nU6\\t\"z}}' "
                f"{tag}.rRNA.out > {tag}.nf",
                shell=True, check=True
            )

            # ShortStack alignment
            subprocess.run(
                f"ShortStack --outdir ShortStack_{tag} --align_only --mmap {mmap} "
                f"--threads {thread} --nohp --readfile {tag}.fastq "
                f"--genomefile {prefix}/reference/{genome}_chr_all.fasta 2>&1",
                shell=True, check=True
            )

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

            # Cleanup
            for fname in (tag, f"{tag}.unmapped.sam", f"{tag}.sam", f"{tag}.fastq", f"{tag}.rRNA.out"):
                if os.path.exists(fname):
                    os.unlink(fname)
            if os.path.exists(f"ShortStack_{tag}"):
                subprocess.run(f"rm -rf ShortStack_{tag}", shell=True, check=True)

            tee.write("\nConverting BAM to BED...\n")

            subprocess.run(f"bamToBed -bed12 -i {tag}.bam > {tag}.bed", shell=True, check=True)

            tee.write("\nGenerating individual files...\n")

            subprocess.run(
                f"awk -F \"\\t\" '{{a=substr(FILENAME,1,length(FILENAME)-3);"
                f"if($11>=18 && $11 <= 26) {{print $0 > (a$11\".bed\")}}}}' {tag}.bed",
                shell=True, check=True
            )
            subprocess.run(
                f"awk '!a[$4]++' {tag}.bed | awk '{{print $11}}' | sort | uniq -c | "
                f"awk '{{OFS=\"\\t\"; print $2, $1}}' > {tag}.len_dist.txt",
                shell=True, check=True
            )
            subprocess.run(
                f"awk '{{if($1<=28){{n+=$2}}}}END{{print \"total\\t\"n}}' "
                f"{tag}.len_dist.txt >> {tag}.nf",
                shell=True, check=True
            )

            tee.write("\nLength distribution summary done!\nCounting start...\n")

            # Read normalization factors
            normhash = {}
            with open(f"{tag}.nf") as nf_fh:
                for line in nf_fh:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        normhash[parts[0]] = int(parts[1])

            for mnorm in norms:
                if mnorm in normhash:
                    _count(mnorm, normhash[mnorm], prefix, genome, binsize, tag, tee)

            tee.write("Counting Completed!\n")

            # Clean up bed files
            for bed_file in globmod.glob(f"{tag}*.bed"):
                os.unlink(bed_file)

        # Clean up mask files
        if mask:
            for mf in globmod.glob("mask*"):
                os.unlink(mf)

        # Statistical analysis
        if not mappingonly and len(pars) > 1:
            for mnorm in norms:
                _stat_analysis(mnorm, prefix, genome, foldchange, pvalue, binsize,
                              promoter_length, par_str, tee)
    else:
        # No-mapping mode: symlink existing files
        for pre in tags:
            os.symlink(f"../{pre}.nf", f"{pre}.nf")
            original_dir = ".."
            for fname in os.listdir(original_dir):
                if fname.startswith(pre) and fname.endswith('count') and 'norm' not in fname:
                    os.symlink(f"../{fname}", fname)

        for mnorm in norms:
            _stat_analysis(mnorm, prefix, genome, foldchange, pvalue, binsize,
                          promoter_length, par_str, tee)

        for fname in globmod.glob("*.count"):
            os.unlink(fname)
        for fname in globmod.glob("*.nf"):
            os.unlink(fname)

    # Clean up GFF files
    for fname in globmod.glob("*.gff"):
        if fname not in ("gene.gff", "te.gff", "promoter.gff"):
            os.unlink(fname)


def _parse_to_dict(arg_str):
    """Parse 'name=value' string to dict."""
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _resolve_path(filepath):
    """Resolve file path (handle ~, relative paths)."""
    if filepath.startswith('~/'):
        return os.path.expanduser(filepath)
    elif not filepath.startswith('/'):
        return os.path.abspath(os.path.join('..', filepath))
    return filepath


def _umi_dedup(tag):
    """UMI deduplication for single-cell mode."""
    dedup = defaultdict(lambda: defaultdict(dict))
    trimmed_file = f"{tag}_trimmed.fastq"

    with open(trimmed_file) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            n = i % 4
            if n == 3:  # quality
                dedup[out_seq_names[2]][out_seq]['quality'] = line
            elif n == 1:  # sequence
                out_seq = line
                dedup[out_seq_names[2]][out_seq]['id'] = out_seq_names[0]
                dedup[out_seq_names[2]][out_seq]['count'] = \
                    dedup[out_seq_names[2]][out_seq].get('count', 0) + 1
            elif n == 0:  # header
                seq_names = line.split()
                out_seq_names = seq_names[0].split('_')

    with open(f"{tag}.fastq", 'w') as fh_out:
        for umi in dedup:
            for output_seq in dedup[umi]:
                data = dedup[umi][output_seq]
                fh_out.write(
                    f"{data['id']}_{umi}_{data['count']} {seq_names[1]}\n"
                    f"{output_seq}\n+\n{data['quality']}\n"
                )


def _count(mnorm, rc, prefix, genome, binsize, tag, tee):
    """Count reads in bins, genes, TEs, promoters, and miRNAs."""
    import math

    fas = ref.read_fasta(prefix, genome)

    # miRNA counting
    mir_gff = os.path.join(prefix, "reference", f"{genome}_miRNA_miRNA_star.gff")
    subprocess.run(
        f"bedtools intersect -a {mir_gff} -b {tag}.bam -wa -f 1 -r -c | "
        f"awk '{{print $9\"\\t\"$10}}' > {tag}.miRNA.tmp",
        shell=True, check=True
    )

    count_data = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    )

    with open(f"{tag}.miRNA.tmp") as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) >= 2:
                count_data['mir_raw'][cols[0]] = int(cols[1])

    if os.path.exists(f"{tag}.miRNA.tmp"):
        os.unlink(f"{tag}.miRNA.tmp")

    # miRNA annotation counting
    mirna_data = defaultdict(lambda: {'name': '', 'count': 0})
    with open(mir_gff) as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            miseq = fas.get(cols[0], "")[int(cols[3]) - 1:int(cols[4])]
            if cols[6] == '-':
                miseq = revcomp(miseq)
            mirna_data[miseq]['name'] += cols[8] + ';'
            mirna_data[miseq]['count'] += count_data['mir_raw'].get(cols[8], 0)

    with open(f"{tag}.miRNA.annotated.count", 'w') as fh:
        for miseq in sorted(mirna_data.keys()):
            mirna_data[miseq]['name'] = mirna_data[miseq]['name'].rstrip(';')
            fh.write(f"{mirna_data[miseq]['name']}\t{mirna_data[miseq]['count']}\n")

    # Process per-length BED files
    bed_files = [f for f in os.listdir('.') if f.startswith(tag) and f.endswith('.bed') and any(c.isdigit() for c in f.replace(tag, '').replace('.bed', ''))]

    for sbed in bed_files:
        nrc = 1000000.0 / rc

        # Strand-specific and combined bedgraph generation
        _make_bedgraph(sbed, prefix, genome, nrc, mnorm)

        # miRNA counting per length
        subprocess.run(
            f"bedtools intersect -a {mir_gff} -b {sbed} -wa -f 0.95 -c | "
            f"awk -v x={rc} '{{print $9\"\\t\"$10 * 1000000 / x}}' > {sbed}.tmp",
            shell=True, check=True
        )
        with open(f"{sbed}.tmp") as fh:
            for line in fh:
                cols = line.strip().split('\t')
                if len(cols) >= 2:
                    count_data['mir'][cols[0]][sbed]['n'] = float(cols[1])

        subprocess.run(
            f"bedtools intersect -a {mir_gff} -b {sbed} -wa -f 0.95 -c | "
            f"awk '{{print $9\"\\t\"$10}}' > {sbed}.tmp",
            shell=True, check=True
        )
        with open(f"{sbed}.tmp") as fh:
            for line in fh:
                cols = line.strip().split('\t')
                if len(cols) >= 2:
                    count_data['mir'][cols[0]][sbed]['r'] = int(cols[1])

        if os.path.exists(f"{sbed}.tmp"):
            os.unlink(f"{sbed}.tmp")

        # Bin counting
        with open(sbed) as fh:
            for line in fh:
                cols = line.strip().split('\t')
                if len(cols) < 12:
                    continue
                bin_idx = (int(cols[1]) + 1) // binsize
                count_data['r100'][cols[0]][bin_idx][sbed]['r'] += 1

        # Fill missing bins
        chrom_lengths = ref.read_chromosome_lengths(prefix, genome, binsize)
        for chr_name in sorted(chrom_lengths.keys()):
            for bi in range(chrom_lengths[chr_name] + 1):
                if chr_name not in count_data['r100'] or bi not in count_data['r100'][chr_name] or sbed not in count_data['r100'][chr_name][bi]:
                    count_data['r100'][chr_name][bi][sbed]['r'] = 0
                count_data['r100'][chr_name][bi][sbed]['n'] = \
                    count_data['r100'][chr_name][bi][sbed].get('r', 0) * 1000000.0 / rc

        # Gene counting
        subprocess.run(
            f"bedtools intersect -a gene.gff -b {sbed} -wa -c > {sbed}.gene.tmp",
            shell=True, check=True
        )
        with open(f"{sbed}.gene.tmp") as fh:
            for line in fh:
                cols = line.strip().split('\t')
                if len(cols) >= 10:
                    count_data['gene'][cols[8]][sbed]['r'] = int(cols[9])
                    count_data['gene'][cols[8]][sbed]['n'] = int(cols[9]) * 1000000.0 / rc
        if os.path.exists(f"{sbed}.gene.tmp"):
            os.unlink(f"{sbed}.gene.tmp")

        # TE counting
        if os.path.exists("te.gff"):
            subprocess.run(
                f"bedtools intersect -a te.gff -b {sbed} -wa -c > {sbed}.te.tmp",
                shell=True, check=True
            )
            with open(f"{sbed}.te.tmp") as fh:
                for line in fh:
                    cols = line.strip().split('\t')
                    if len(cols) >= 10:
                        count_data['te'][cols[8]][sbed]['r'] = int(cols[9])
                        count_data['te'][cols[8]][sbed]['n'] = int(cols[9]) * 1000000.0 / rc
            if os.path.exists(f"{sbed}.te.tmp"):
                os.unlink(f"{sbed}.te.tmp")

        # Promoter counting
        if os.path.exists("promoter.gff"):
            subprocess.run(
                f"bedtools intersect -a promoter.gff -b {sbed} -wa -c > {sbed}.promoter.tmp",
                shell=True, check=True
            )
            with open(f"{sbed}.promoter.tmp") as fh:
                for line in fh:
                    cols = line.strip().split('\t')
                    if len(cols) >= 10:
                        count_data['promoter'][cols[8]][sbed]['r'] = int(cols[9])
                        count_data['promoter'][cols[8]][sbed]['n'] = int(cols[9]) * 1000000.0 / rc
            if os.path.exists(f"{sbed}.promoter.tmp"):
                os.unlink(f"{sbed}.promoter.tmp")

    # Write count files
    _write_count_files(count_data, tag, mnorm, bed_files, prefix, genome, fas)


def _make_bedgraph(sbed, prefix, genome, nrc, mnorm):
    """Generate bedgraph and bigwig files."""
    fai = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta.fai")

    for strand_suffix, strand_flag, sign in [('p', '+', ''), ('n', '-', '-')]:
        bg_file = sbed.replace('.bed', f'{strand_suffix}.bedgraph')
        subprocess.run(
            f"bedtools genomecov -split -strand {strand_flag} -bg -i {sbed} -g {fai} > {bg_file}",
            shell=True, check=True
        )

        bg_norm = bg_file.replace('.bedgraph', f'.{mnorm}.bedgraph')
        bw_norm = bg_file.replace('.bedgraph', f'.{mnorm}.bw')
        subprocess.run(
            f"bedtools genomecov -split -strand {strand_flag} -scale {sign}{nrc} "
            f"-bg -i {sbed} -g {fai} > {bg_norm}",
            shell=True, check=True
        )
        subprocess.run(f"bedGraphToBigWig {bg_norm} {fai} {bw_norm}", shell=True, check=True)
        os.unlink(bg_file)
        os.unlink(bg_norm)

    # Combined (non-strand-specific)
    bg_file = sbed.replace('.bed', '.bedgraph')
    subprocess.run(
        f"bedtools genomecov -split -bg -i {sbed} -g {fai} > {bg_file}",
        shell=True, check=True
    )
    bg_norm = bg_file.replace('.bedgraph', f'.{mnorm}.bedgraph')
    bw_norm = bg_file.replace('.bedgraph', f'.{mnorm}.bw')
    subprocess.run(
        f"bedtools genomecov -split -scale {nrc} -bg -i {sbed} -g {fai} > {bg_norm}",
        shell=True, check=True
    )
    subprocess.run(f"bedGraphToBigWig {bg_norm} {fai} {bw_norm}", shell=True, check=True)
    os.unlink(bg_file)
    os.unlink(bg_norm)


def _write_count_files(count_data, tag, mnorm, bed_files, prefix, genome, fas):
    """Write count output files."""
    # Bin counts
    with open(f"{tag}.count", 'w') as fh1, open(f"{tag}.{mnorm}.norm.count", 'w') as fh2:
        for chr_name in sorted(count_data.get('r100', {}).keys()):
            for bi in sorted(count_data['r100'][chr_name].keys()):
                fh1.write(f"{chr_name}_{bi}")
                fh2.write(f"{chr_name}_{bi}")
                for sbed in sorted(bed_files):
                    r_val = count_data['r100'][chr_name][bi].get(sbed, {}).get('r', 0)
                    n_val = count_data['r100'][chr_name][bi].get(sbed, {}).get('n', 0)
                    fh1.write(f"\t{r_val}")
                    fh2.write(f"\t{n_val}")
                fh1.write('\n')
                fh2.write('\n')

    # Gene counts
    with open(f"{tag}.gene.count", 'w') as fh1, open(f"{tag}.gene.{mnorm}.norm.count", 'w') as fh2:
        for gene_name in sorted(count_data.get('gene', {}).keys()):
            fh1.write(gene_name)
            fh2.write(gene_name)
            for sbed in sorted(bed_files):
                r_val = count_data['gene'][gene_name].get(sbed, {}).get('r', 0)
                n_val = count_data['gene'][gene_name].get(sbed, {}).get('n', 0)
                fh1.write(f"\t{r_val}")
                fh2.write(f"\t{n_val}")
            fh1.write('\n')
            fh2.write('\n')

    # TE counts
    with open(f"{tag}.TE.count", 'w') as fh1, open(f"{tag}.TE.{mnorm}.norm.count", 'w') as fh2:
        for gene_name in sorted(count_data.get('te', {}).keys()):
            fh1.write(gene_name)
            fh2.write(gene_name)
            for sbed in sorted(bed_files):
                r_val = count_data['te'][gene_name].get(sbed, {}).get('r', 0)
                n_val = count_data['te'][gene_name].get(sbed, {}).get('n', 0)
                fh1.write(f"\t{r_val}")
                fh2.write(f"\t{n_val}")
            fh1.write('\n')
            fh2.write('\n')

    # Promoter counts
    with open(f"{tag}.promoter.count", 'w') as fh1, open(f"{tag}.promoter.{mnorm}.norm.count", 'w') as fh2:
        for gene_name in sorted(count_data.get('promoter', {}).keys()):
            fh1.write(gene_name)
            fh2.write(gene_name)
            for sbed in sorted(bed_files):
                r_val = count_data['promoter'][gene_name].get(sbed, {}).get('r', 0)
                n_val = count_data['promoter'][gene_name].get(sbed, {}).get('n', 0)
                fh1.write(f"\t{r_val}")
                fh2.write(f"\t{n_val}")
            fh1.write('\n')
            fh2.write('\n')

    # miRNA counts
    mir_gff = os.path.join(prefix, "reference", f"{genome}_miRNA_miRNA_star.gff")
    mirna_data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    with open(mir_gff) as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            miseq = fas.get(cols[0], "")[int(cols[3]) - 1:int(cols[4])]
            if cols[6] == '-':
                miseq = revcomp(miseq)
            mirna_data[miseq]['name'] = cols[8]

            for sbed in bed_files:
                if cols[8] in count_data.get('mir', {}):
                    mirna_data[miseq][sbed]['r'] += count_data['mir'][cols[8]].get(sbed, {}).get('r', 0)
                    mirna_data[miseq][sbed]['n'] += count_data['mir'][cols[8]].get(sbed, {}).get('n', 0)

    with open(f"{tag}.miRNA.count", 'w') as fh1, open(f"{tag}.miRNA.{mnorm}.norm.count", 'w') as fh2:
        for miseq in sorted(mirna_data.keys()):
            mir_name = ""
            for mir_id_data in mirna_data[miseq].values():
                if isinstance(mir_id_data, str):
                    mir_name += mir_id_data + ';'
            mir_name = mir_name.rstrip(';')
            fh1.write(mir_name)
            fh2.write(mir_name)
            for sbed in sorted(bed_files):
                r_val = mirna_data[miseq].get(sbed, {}).get('r', 0)
                n_val = mirna_data[miseq].get(sbed, {}).get('n', 0)
                fh1.write(f"\t{r_val}")
                fh2.write(f"\t{n_val}")
            fh1.write('\n')
            fh2.write('\n')


def _stat_analysis(mnorm, prefix, genome, foldchange, pvalue, binsize,
                   promoter_length, par_str, tee):
    """Run statistical analyses via R scripts."""
    tee.write(f"\nDSR analysis...\nNormalization {mnorm}\tFold Change {foldchange}\tP Value {pvalue}\n")

    subprocess.run(
        f"Rscript --vanilla {prefix}/scripts/DSR.R {mnorm} {pvalue} {foldchange} {par_str}",
        shell=True, check=True
    )

    # Generate bigwig from hyper/hypo CSV
    csv_files = [f for f in os.listdir('.') if f.endswith('.csv') and mnorm in f and 'bin' in f]
    hyper_hypo = [f for f in csv_files if 'hyper' in f or 'hypo' in f]

    fai = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta.fai")
    for hcsv in hyper_hypo:
        bg_file = hcsv.replace('.csv', '.bedgraph')
        with open(hcsv) as fh_in, open(bg_file, 'w') as fh_out:
            for line in fh_in:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                cols = line.split(',')
                if not cols:
                    continue
                m = __import__('re').match(r'(\w+)_(\d+)', cols[0].strip('"'))
                if m:
                    chr_name = m.group(1)
                    start = int(m.group(2)) * binsize
                    end = start + binsize - 1
                    fh_out.write(f"{chr_name}\t{start}\t{end}\t{cols[2] if len(cols) > 2 else '0'}\n")

        bw_file = hcsv.replace('.csv', '.bw')
        subprocess.run(f"bedGraphToBigWig {bg_file} {fai} {bw_file}", shell=True, check=True)
        os.unlink(bg_file)

    # Annotate CSV files
    ann = ref.build_annotation(prefix, genome, binsize, promoter_length)
    for csv_file in csv_files:
        _annotate_csv(csv_file, ann)

    # miRNA DE, gene DE, TE DE, promoter DE
    tee.write(f"\nDE miRNA analysis...\nNormalization {mnorm}\tFold Change {foldchange}\tP Value {pvalue}\n")
    subprocess.run(
        f"Rscript --vanilla {prefix}/scripts/DEM.R {mnorm} {pvalue} {foldchange} {par_str}",
        shell=True, check=True
    )

    tee.write(f"\nDS gene analysis...\nNormalization {mnorm}\tFold Change {foldchange}\tP Value {pvalue}\n")
    subprocess.run(
        f"Rscript --vanilla {prefix}/scripts/DSG.R {mnorm} {pvalue} {foldchange} {par_str}",
        shell=True, check=True
    )

    tee.write(f"\nDS TE analysis...\nNormalization {mnorm}\tFold Change {foldchange}\tP Value {pvalue}\n")
    subprocess.run(
        f"Rscript --vanilla {prefix}/scripts/DST.R {mnorm} {pvalue} {foldchange} {par_str}",
        shell=True, check=True
    )

    tee.write(f"\nDS Promoter analysis...\nNormalization {mnorm}\tFold Change {foldchange}\tP Value {pvalue}\n")
    subprocess.run(
        f"Rscript --vanilla {prefix}/scripts/DSP.R {mnorm} {pvalue} {foldchange} {par_str}",
        shell=True, check=True
    )


def _annotate_csv(csv_file, ann):
    """Add genome annotation to CSV file."""
    tmp_file = "tmp4"
    with open(csv_file) as fh_in, open(tmp_file, 'w') as fh_out:
        for line in fh_in:
            line = line.strip()
            if not line:
                continue
            cols = line.split(',')
            if not cols:
                continue
            key = cols[0].strip('"')
            if key in ann:
                fh_out.write(f"{line},{ann[key]}\n")
            else:
                fh_out.write(f"{line},Intergenic\n")
    os.rename(tmp_file, csv_file)
