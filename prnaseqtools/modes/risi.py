"""
risiRNA analysis mode.
ShortStack alignment to rDNA genome → counting → DESeq2 DE analysis.
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
    """Main entry point for risiRNA analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    mmap = opts.get('mmap', 'u')
    nomapping = opts.get('no_mapping', False)
    mappingonly = opts.get('mapping_only', False)
    foldchange = opts.get('foldchange', 1.5)
    pvalue = opts.get('pvalue', 0.01)
    binsize = opts.get('binsize', 10)
    mnorm = opts.get('norm', 'total')

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

            tee.write("\nStart mapping...\n")
            rDNA_fasta = os.path.join(prefix, "reference", f"{genome}_rDNA_chr_all.fasta")
            subprocess.run(
                f"ShortStack --outdir ShortStack_{tag} --align_only --bowtie_m 1000 "
                f"--ranmax 50 --mmap {mmap} --mismatches 1 --bowtie_cores {thread} "
                f"--nohp --readfile {tag}.fastq --genomefile {rDNA_fasta} 2>&1",
                shell=True, check=True
            )

            tee.write("\nAlignment Completed!\n")

            subprocess.run(f"samtools view -h ShortStack_{tag}/{tag}.bam > {tag}", shell=True, check=True)
            subprocess.run(
                f"awk '{{if($0~/^@/) print > (FILENAME\".unmapped.sam\"); "
                f"if($10!=\"*\" && $3!=\"*\") print > (FILENAME\".sam\"); "
                f"if($10!=\"*\" && $3==\"*\") print > (FILENAME\".unmapped.sam\")}}' {tag}",
                shell=True, check=True
            )
            subprocess.run(f"samtools view -Sb {tag}.unmapped.sam > {tag}.unmapped.bam", shell=True, check=True)
            subprocess.run(f"samtools view -Sb {tag}.sam > {tag}.bam", shell=True, check=True)

            for fname in (tag, f"{tag}.unmapped.sam", f"{tag}.sam", f"{tag}.fastq"):
                if os.path.exists(fname):
                    os.unlink(fname)
            if os.path.exists(f"ShortStack_{tag}"):
                subprocess.run(f"rm -rf ShortStack_{tag}", shell=True)

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
                f"awk '{{n+=$2}}END{{print \"total\\t\"n}}' {tag}.len_dist.txt >> {tag}.nf",
                shell=True, check=True
            )

            tee.write("\nLength distribution summary done!\nCounting start...\n")

            normhash = {}
            with open(f"{tag}.nf") as fh:
                for line in fh:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        normhash[parts[0]] = int(parts[1])

            if mnorm in normhash:
                _count_risi(mnorm, normhash[mnorm], prefix, genome, binsize, tag, tee)

            tee.write("Counting Completed!\n")
            for bed_file in globmod.glob(f"{tag}*.bed"):
                os.unlink(bed_file)

        if not mappingonly and len(pars) > 1:
            _stat_analysis(mnorm, prefix, genome, foldchange, pvalue, binsize, par_str, tee)
    else:
        for pre in tags:
            os.symlink(f"../{pre}.nf", f"{pre}.nf")
            for fname in os.listdir(".."):
                if fname.startswith(pre) and fname.endswith('count') and 'norm' not in fname:
                    os.symlink(f"../{fname}", fname)

        _stat_analysis(mnorm, prefix, genome, foldchange, pvalue, binsize, par_str, tee)

        for fname in globmod.glob("*.count"):
            os.unlink(fname)
        for fname in globmod.glob("*.nf"):
            os.unlink(fname)

    for fname in globmod.glob("*.gff"):
        os.unlink(fname)


def _parse_to_dict(arg_str):
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _count_risi(mnorm, rc, prefix, genome, binsize, tag, tee):
    """Count reads in rDNA bins and features."""
    count_data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    lengths = ref.read_chromosome_lengths(prefix, f"{genome}_rDNA", binsize)

    bed_files = [f for f in os.listdir('.') if f.startswith(tag) and any(c.isdigit() for c in f.replace(tag, '').replace('.bed', ''))]

    for sbed in bed_files:
        nrc = 1000000.0 / rc
        fai = os.path.join(prefix, "reference", f"{genome}_rDNA_chr_all.fasta.fai")

        _make_bedgraph(sbed, prefix, genome, nrc, mnorm, fai)

        # Bin counting
        leng = int(''.join(c for c in sbed.replace(tag, '').replace('.bed', '') if c.isdigit()) or 0)
        with open(sbed) as fh:
            for line in fh:
                cols = line.strip().split('\t')
                if len(cols) < 6:
                    continue
                bin_idx = (int(cols[1]) + 1) // binsize
                count_data['r10'][cols[0]][bin_idx][cols[5]][leng] += 1

        # Feature counting
        rDNA_gff = os.path.join(prefix, "reference", f"{genome}_rDNA.gff")
        subprocess.run(
            f"bedtools intersect -a {sbed} -b {rDNA_gff} -wb > {sbed}.tmp",
            shell=True, check=True
        )

        tmp_hash = defaultdict(lambda: {'strand': '', 'feature': set(), 'length': 0})
        with open(f"{sbed}.tmp") as fh:
            for line in fh:
                cols = line.strip().split('\t')
                if len(cols) >= 21:
                    tmp_hash[cols[3]]['strand'] = cols[5]
                    tmp_hash[cols[3]]['feature'].add(cols[20])
                    tmp_hash[cols[3]]['length'] = cols[10]

        for read_id, data in tmp_hash.items():
            feature_str = ','.join(sorted(data['feature']))
            count_data['feature'][feature_str][data['strand']][int(data['length']) if data['length'] else 0] += 1

        if os.path.exists(f"{sbed}.tmp"):
            os.unlink(f"{sbed}.tmp")

    # Write count files
    _write_risi_counts(count_data, tag, mnorm, rc, bed_files, lengths)


def _make_bedgraph(sbed, prefix, genome, nrc, mnorm, fai):
    """Generate rDNA bedgraph."""
    for suffix, strand, sign in [('p', '+', ''), ('n', '-', '-')]:
        bg = sbed.replace('.bed', f'{suffix}.bedgraph')
        subprocess.run(f"bedtools genomecov -split -strand {strand} -bg -i {sbed} -g {fai} > {bg}", shell=True, check=True)
        bg_norm = bg.replace('.bedgraph', f'.{mnorm}.bedgraph')
        bw_norm = bg.replace('.bedgraph', f'.{mnorm}.bw')
        subprocess.run(f"bedtools genomecov -split -strand {strand} -scale {sign}{nrc} -bg -i {sbed} -g {fai} > {bg_norm}", shell=True, check=True)
        subprocess.run(f"bedGraphToBigWig {bg_norm} {fai} {bw_norm}", shell=True, check=True)
        os.unlink(bg)
        os.unlink(bg_norm)

    bg = sbed.replace('.bed', '.bedgraph')
    subprocess.run(f"bedtools genomecov -split -bg -i {sbed} -g {fai} > {bg}", shell=True, check=True)
    bg_norm = bg.replace('.bedgraph', f'.{mnorm}.bedgraph')
    bw_norm = bg.replace('.bedgraph', f'.{mnorm}.bw')
    subprocess.run(f"bedtools genomecov -split -scale {nrc} -bg -i {sbed} -g {fai} > {bg_norm}", shell=True, check=True)
    subprocess.run(f"bedGraphToBigWig {bg_norm} {fai} {bw_norm}", shell=True, check=True)
    os.unlink(bg)
    os.unlink(bg_norm)


def _write_risi_counts(count_data, tag, mnorm, rc, bed_files, lengths):
    """Write risiRNA count files."""
    with open(f"{tag}.count", 'w') as fh1, open(f"{tag}.{mnorm}.norm.count", 'w') as fh2:
        for chr_name in sorted(count_data.get('r10', {}).keys()):
            max_bin = lengths.get(chr_name, 0)
            for bi in range(max_bin + 1):
                for strand in ('+', '-'):
                    fh1.write(f"{chr_name}_{bi}_{strand}")
                    fh2.write(f"{chr_name}_{bi}_{strand}")
                    for leng in range(18, 27):
                        r_val = count_data['r10'].get(chr_name, {}).get(bi, {}).get(strand, {}).get(leng, 0)
                        fh1.write(f"\t{r_val}")
                        fh2.write(f"\t{r_val * 1000000.0 / rc}")
                    fh1.write('\n')
                    fh2.write('\n')

    with open(f"{tag}.feature.count", 'w') as fh1, open(f"{tag}.feature.{mnorm}.norm.count", 'w') as fh2:
        for feature in sorted(count_data.get('feature', {}).keys()):
            for strand in ('+', '-'):
                fh1.write(f"{feature}_{strand}")
                fh2.write(f"{feature}_{strand}")
                for leng in range(18, 27):
                    r_val = count_data['feature'].get(feature, {}).get(strand, {}).get(leng, 0)
                    fh1.write(f"\t{r_val}")
                    fh2.write(f"\t{r_val * 1000000.0 / rc}")
                fh1.write('\n')
                fh2.write('\n')


def _stat_analysis(mnorm, prefix, genome, foldchange, pvalue, binsize, par_str, tee):
    """Run DSR and DSF analyses."""
    tee.write(f"\nDSR analysis...\nNormalization {mnorm}\tFold Change {foldchange}\tP Value {pvalue}\n")
    subprocess.run(
        f"Rscript --vanilla {prefix}/scripts/DSR.R {mnorm} {pvalue} {foldchange} {par_str}",
        shell=True, check=True
    )

    fai = os.path.join(prefix, "reference", f"{genome}_rDNA_chr_all.fasta.fai")
    csv_files = [f for f in os.listdir('.') if f.endswith('.csv') and mnorm in f and 'bin' in f]
    for hcsv in [f for f in csv_files if 'hyper' in f or 'hypo' in f]:
        bg_file = hcsv.replace('.csv', '.bedgraph')
        with open(hcsv) as fh_in, open(bg_file, 'w') as fh_out:
            for line in fh_in:
                line = line.strip()
                if not line:
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

    tee.write(f"\nDS feature analysis...\nNormalization {mnorm}\tFold Change {foldchange}\tP Value {pvalue}\n")
    subprocess.run(
        f"Rscript --vanilla {prefix}/scripts/DSF.R {mnorm} {pvalue} {foldchange} {par_str}",
        shell=True, check=True
    )
