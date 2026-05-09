"""
Truncation and tailing analysis for miRNAs.
Perl mode::tt equivalent.
bowtie alignment with iterative mismatches → read categorization → bubble plots.
"""

import os
import sys
import subprocess
from pathlib import Path
from collections import defaultdict

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import parse_input
from prnaseqtools.functions import download_sra, unzip_file, revcomp, _tee
from prnaseqtools import reference as ref


def run(opts):
    """Main entry point for truncation/tailing analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    mmap = opts.get('mmap', 'u')

    control_dict = _parse_to_dict(opts.get('control', ''))
    tags, files, pars = parse_input(control_dict)

    if opts.get('treatment'):
        treatment_dict = _parse_to_dict(opts.get('treatment', ''))
        t_tags, t_files, t_pars = parse_input(treatment_dict)
        tags.extend(t_tags)
        files.extend(t_files)
        pars.extend(t_pars)

    fas = ref.read_fasta(prefix, genome)

    for i in range(len(tags)):
        tag = tags[i]
        fpath = files[i]

        tee.write(f"\nMapping {tag}...\n")

        sra_results = download_sra(fpath, thread)
        unzip_file(sra_results[0], tag)

        if adaptor:
            tee.write("\nStart trimming...\n")
            subprocess.run(
                f"cutadapt -j {thread} -m 14 -M 42 --discard-untrimmed --trim-n "
                f"-a {adaptor} -o {tag}_trimmed.fastq {tag}.fastq 2>&1",
                shell=True, check=True
            )
            os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")

        # Iterative bowtie mapping with increasing mismatches
        ref_idx = os.path.join(prefix, "reference", f"{genome}_chr_all")
        subprocess.run(
            f"bowtie -v 0 -p {thread} -t --un {tag}.unmapped_0.fastq "
            f"--al {tag}.mapped.fastq {ref_idx} {tag}.fastq {tag}.out 2>&1",
            shell=True, check=True
        )

        for p in range(1, 9):
            j = p - 1
            with open(f"{tag}.unmapped_{j}.fastq") as fq_in, open("tmp.fastq", 'w') as fq_out:
                for n, line in enumerate(fq_in):
                    line = line.strip()
                    if n % 4 == 0:
                        parts = line.split()
                        row0 = parts[0].split('_')
                        seq_id = row0[0]
                        ss = row0[1] if len(row0) > 1 else ''
                        seq = None
                    elif n % 4 == 1:
                        seq = line[:-1]
                        ss += line[-1].lower()
                    elif n % 4 == 3:
                        qua = line[:-1]
                        if len(seq) >= 14:
                            fq_out.write(f"{seq_id}_{ss}\n{seq}\n+\n{qua}\n")

            os.unlink(f"{tag}.unmapped_{j}.fastq")

            subprocess.run(
                f"bowtie -v 0 -p {thread} -t --un {tag}.unmapped_{p}.fastq "
                f"--al {tag}.mapped_{p}.fastq {ref_idx} tmp.fastq {tag}.out 2>&1",
                shell=True, check=True
            )

        # Concatenate mapped reads
        subprocess.run(f"cat {tag}.mapped*.fastq > {tag}_edited.fastq", shell=True, check=True)
        for fname in ["tmp.fastq"]:
            if os.path.exists(fname):
                os.unlink(fname)

        # ShortStack alignment
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        subprocess.run(
            f"ShortStack --outdir {tag}tmp --align_only --bowtie_m 1000 "
            f"--ranmax 50 --mmap {mmap} --mismatches 0 --bowtie_cores {thread} "
            f"--nohp --readfile {tag}_edited.fastq --genomefile {fasta_path} 2>&1",
            shell=True, check=True
        )

        subprocess.run(f"samtools view -h {tag}tmp/{tag}_edited.bam > {tag}", shell=True, check=True)
        subprocess.run(
            f"awk '{{if($10!=\"*\" && $3!=\"*\") print > (FILENAME\".edited.sam\")}}' {tag}",
            shell=True, check=True
        )
        subprocess.run(f"samtools view -Sb {tag}.edited.sam > {tag}.edited.bam", shell=True, check=True)

        subprocess.run(f"rm -rf {tag}tmp", shell=True)
        for fname in (f"{tag}_edited.fastq", f"{tag}.fastq", tag, f"{tag}.edited.sam"):
            if os.path.exists(fname):
                os.unlink(fname)
        for fname in [f for f in os.listdir('.') if f.startswith(f"{tag}.mapped")]:
            os.unlink(fname)

        tee.write("\nAlignment Completed!\n")

        # miRNA intersection
        mir_gff = os.path.join(prefix, "reference", f"{genome}_miRNA_miRNA_star.gff")
        subprocess.run(
            f"bedtools intersect -wo -s -a {mir_gff} -b {tag}.edited.bam > {tag}.out",
            shell=True, check=True
        )

        mir_data = ref.read_mirna_gff(prefix, genome)

        with open(f"{tag}.out") as fh:
            for line in fh:
                cols = line.strip().split('\t')
                if len(cols) < 13:
                    continue
                parts = cols[12].split('_')
                tail = revcomp(parts[1]) if len(parts) > 1 and parts[1] else ''
                key = f"{cols[10]}_{cols[11]}_{tail}"
                if 'read' not in mir_data.get(cols[8], {}):
                    mir_data.setdefault(cols[8], {})['read'] = defaultdict(int)
                mir_data[cols[8]]['read'][key] += 1

        # Output
        flank = 25
        out_data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

        with open(f"{tag}.seq.out", 'w') as fh_out:
            for mm in sorted(mir_data.keys()):
                if mm not in mir_data or 'chromosome' not in mir_data.get(mm, {}):
                    continue
                mdata = mir_data[mm]
                lim1 = mdata['start'] - flank
                lim2 = mdata['end'] + flank
                length = mdata['end'] - mdata['start'] + 1
                refseq = fas.get(mdata['chromosome'], '')[lim1:lim2 + 1]
                miseq = fas.get(mdata['chromosome'], '')[mdata['start'] - 1:mdata['end']]

                if mdata['strand'] == '+':
                    fh_out.write(f"{refseq}\t{mm}\n{'.' * (flank - 1)}{miseq}{'.' * (flank + 1)}\tREF\n")
                else:
                    refseq = revcomp(refseq)
                    miseq = revcomp(miseq)
                    fh_out.write(f"{refseq}\t{mm}\n{'.' * (flank + 1)}{miseq}{'.' * (flank - 1)}\tREF\n")

                for read_key, rcount in mdata.get('read', {}).items():
                    features = read_key.split('_')
                    ftail = features[2] if len(features) > 2 else ''
                    fstart = int(features[0])
                    fend = int(features[1])

                    if mdata['strand'] == '+':
                        seq = fas.get(mdata['chromosome'], '')[fstart:fend]
                        flank1 = fstart - lim1
                        flank2 = lim2 - fend - len(ftail)

                        if mdata['start'] == fstart + 1:
                            if length - fend + fstart - 1 < 0:
                                out_data[mm][fend - fstart - length + len(ftail)][0] += rcount
                            else:
                                out_data[mm][len(ftail)][length - fend + fstart] += rcount
                    else:
                        seq = revcomp(fas.get(mdata['chromosome'], '')[fstart:fend])
                        flank1 = lim2 - fend + 1
                        flank2 = fstart - lim1 - len(ftail) - 1

                        if mdata['end'] == fend:
                            if length - fend + fstart - 1 < 0:
                                out_data[mm][fend - fstart - length + len(ftail)][0] += rcount
                            else:
                                out_data[mm][len(ftail)][length - fend + fstart] += rcount

                    fh_out.write(f"{'.' * flank1}{seq}{ftail}{'.' * (flank2 + 1)}\t{rcount}\n")

        # Write categorised output
        with open(f"{tag}.out", 'w') as fh:
            for mm in sorted(out_data.keys()):
                for tailing in range(0, 9):
                    for truncation in range(0, 9):
                        count = out_data[mm].get(tailing, {}).get(truncation, 0)
                        fh.write(f"{mm}\t{tailing}\t{truncation}\t{count}\n")

        # Also process MIR GFF
        _process_mir(prefix, genome, tag, fas)

    par_str = ' '.join(map(str, pars))
    subprocess.run(
        f"Rscript --vanilla {prefix}/scripts/bubble_plot.R {par_str}",
        shell=True, check=True
    )


def _parse_to_dict(arg_str):
    parts = arg_str.split('=')
    if len(parts) == 2:
        return {parts[0]: parts[1]}
    return {}


def _process_mir(prefix, genome, tag, fas):
    """Process MIR (precursor) GFF for truncation/tailing analysis."""
    mir_gff = os.path.join(prefix, "reference", f"{genome}_MIR.gff")
    mirna_gff = os.path.join(prefix, "reference", f"{genome}_miRNA_miRNA_star.gff")

    if not os.path.exists(mir_gff):
        return

    mir_data = ref.read_mirna_gff(prefix, genome)
    miR_data = {}

    with open(mir_gff) as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            miR_data[cols[8]] = {
                'chromosome': cols[0],
                'strand': cols[6],
                'start': int(cols[3]),
                'end': int(cols[4]),
            }

    subprocess.run(
        f"bedtools intersect -wo -s -a {mir_gff} -b {tag}.edited.bam > {tag}.out2",
        shell=True, check=True
    )

    with open(f"{tag}.out2") as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) < 13:
                continue
            parts = cols[12].split('_')
            tail = revcomp(parts[1]) if len(parts) > 1 and parts[1] else ''
            key = f"{cols[10]}_{cols[11]}_{tail}"
            if 'read' not in miR_data.get(cols[8], {}):
                miR_data.setdefault(cols[8], {})['read'] = defaultdict(int)
            miR_data[cols[8]]['read'][key] += 1

    flank = 25
    with open(f"{tag}.seq.out2", 'w') as fh_out:
        for mm in sorted(miR_data.keys()):
            if mm not in miR_data:
                continue
            mdata = miR_data[mm]
            lim1 = mdata['start'] - flank
            lim2 = mdata['end'] + flank
            refseq = fas.get(mdata['chromosome'], '')[lim1:lim2 + 1]

            # Find miRNA positions within MIR
            locs = []
            for mip in sorted(mir_data.keys()):
                if mm.upper() in mip.upper():
                    locs.append(mir_data[mip]['start'])
                    locs.append(mir_data[mip]['end'])
            locs = sorted(set(locs))

            miseq1 = fas.get(mdata['chromosome'], '')[locs[0] - 1:locs[1]]
            miseq2 = fas.get(mdata['chromosome'], '')[locs[2] - 1:locs[3]] if len(locs) > 2 else ''

            if mdata['strand'] == '+':
                if len(locs) > 2:
                    dots = f"{'.' * (locs[0] - lim1 - 1)}{miseq1}{'.' * (locs[2] - locs[1] - 1)}{miseq2}{'.' * (lim2 - locs[3] + 1)}"
                else:
                    dots = f"{'.' * (locs[0] - lim1 - 1)}{miseq1}{'.' * (lim2 - locs[1] + 1)}"
            else:
                refseq = revcomp(refseq)
                miseq1 = revcomp(miseq1)
                miseq2 = revcomp(miseq2) if miseq2 else ''
                if len(locs) > 2:
                    dots = f"{'.' * (lim2 - locs[3] + 1)}{miseq2}{'.' * (locs[2] - locs[1] - 1)}{miseq1}{'.' * (locs[0] - lim1 - 1)}"
                else:
                    dots = f"{'.' * (lim2 - locs[1] + 1)}{miseq1}{'.' * (locs[0] - lim1 - 1)}"

            fh_out.write(f"{refseq}\t{mm}\n{dots}\tREF\n")

            for read_key, rcount in mdata.get('read', {}).items():
                features = read_key.split('_')
                ftail = features[2] if len(features) > 2 else ''
                fstart = int(features[0])
                fend = int(features[1])

                if mdata['strand'] == '+':
                    seq = fas.get(mdata['chromosome'], '')[fstart:fend]
                    flank1 = fstart - lim1
                    flank2 = lim2 - fend - len(ftail)
                else:
                    seq = revcomp(fas.get(mdata['chromosome'], '')[fstart:fend])
                    flank1 = lim2 - fend + 1
                    flank2 = fstart - lim1 - len(ftail) - 1

                fh_out.write(f"{'.' * flank1}{seq}{ftail}{'.' * (flank2 + 1)}\t{rcount}\n")

    if os.path.exists(f"{tag}.out2"):
        os.unlink(f"{tag}.out2")
