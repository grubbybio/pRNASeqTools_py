"""
Reference genome handling for pRNASeqTools.
Mirrors the Perl Ref.pm module.
Provides: GFF parsing, FASTA reading, exon extraction, annotation, etc.
"""

import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from prnaseqtools.functions import revcomp


def _tee():
    try:
        from prnaseqtools.cli import TEE
        if TEE:
            return TEE
    except (ImportError, AttributeError):
        pass
    return sys.stderr


def _prefix():
    """Get the package root directory."""
    return str(Path(__file__).resolve().parent.parent)


def read_gff(prefix, genome):
    """
    Parse gene GFF file.
    Returns dict: {chr: {bin_index: {gene_id: {start, end, exon}}} }
    """
    gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
    index = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    current_id = ""

    with open(gff_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue

            feature = cols[2]
            # Skip UTR, c_transcript, region, protein, CDS, non-gene RNAs
            if 'UTR' in feature or 'c_transcript' in feature or 'region' in feature:
                continue
            if feature == 'protein' or feature == 'CDS' or ('RNA' in feature and 'iRNA' not in feature and 'miRNA' not in feature):
                continue

            if 'gene' in feature:
                m = __import__('re').search(r'^ID=([^;]+);', cols[8])
                if m:
                    current_id = m.group(1)
                    ind = int(cols[3]) // 100000
                    chr_name = cols[0]
                    for offset in (-1, 0, 1):
                        idx_key = ind + offset
                        index[chr_name][idx_key][current_id]['start'] = cols[3]
                        index[chr_name][idx_key][current_id]['end'] = cols[4]

            elif f"Parent={current_id}" in cols[8]:
                ind = int(cols[3]) // 100000
                chr_name = cols[0]
                for offset in (-1, 0, 1):
                    idx_key = ind + offset
                    if 'exon' not in index[chr_name][idx_key][current_id]:
                        index[chr_name][idx_key][current_id]['exon'] = ""
                    index[chr_name][idx_key][current_id]['exon'] += f"{cols[3]}\t{cols[4]};"

    # Convert defaultdict to regular dict for return
    return {k: {kk: dict(vv) for kk, vv in v.items()} for k, v in index.items()}


def read_fasta(prefix, genome):
    """
    Read genome FASTA file.
    Returns dict: {chromosome_name: sequence}
    """
    fas_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
    fas = {}
    current_name = None

    with open(fas_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                current_name = line[1:].split()[0]
                fas[current_name] = ""
            else:
                fas[current_name] += line

    return fas


def read_exons(prefix, genome):
    """
    Extract exons using gffread, return dict: {gene: {transcript: sequence}}.
    """
    exons_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
    gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")

    subprocess.run(
        f"gffread -O -w exons.fa -g {exons_path} {gff_path}",
        shell=True, check=True
    )

    exon_data = defaultdict(lambda: defaultdict(str))
    gene = None
    tran = 0

    with open("exons.fa") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                m = __import__('re').match(r'(\w+)\.(\d+)', header)
                if m:
                    gene = m.group(1)
                    tran = m.group(2)
                else:
                    gene = header
                    tran = 0
            else:
                exon_data[gene][tran] += line

    os.unlink("exons.fa")
    return dict(exon_data)


def get_gene_info(prefix, genome):
    """
    Get gene length info (longest transcript per gene).
    Writes 'transcripts.fa'.
    Returns dict: {gene: length}
    """
    exon_data = read_exons(prefix, genome)
    gene_info = {}

    with open("transcripts.fa", 'w') as tra_fh:
        for gene in sorted(exon_data.keys()):
            longest = 0
            longest_seq = ""
            longest_tran = None

            for tran, seq in exon_data[gene].items():
                if len(seq) > longest:
                    longest = len(seq)
                    longest_seq = seq
                    longest_tran = tran

            gene_info[gene] = longest
            if longest_tran and longest_tran != 0:
                tra_fh.write(f">{gene}.{longest_tran}\n{longest_seq}\n")
            else:
                tra_fh.write(f">{gene}\n{longest_seq}\n")

    return gene_info


def read_mirna_gff(prefix, genome):
    """
    Parse miRNA GFF file.
    Returns dict: {mirna_id: {chromosome, strand, start, end}}
    """
    gff_path = os.path.join(prefix, "reference", f"{genome}_miRNA_miRNA_star.gff")
    mir_data = {}

    with open(gff_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            mir_data[cols[8]] = {
                'chromosome': cols[0],
                'strand': cols[6],
                'start': int(cols[3]),
                'end': int(cols[4]),
            }

    return mir_data


def split_gff(prefix, genome, promoter_length=1000):
    """
    Split GFF into gene.gff, te.gff, promoter.gff.
    """
    gene_gff = os.path.join(prefix, "reference", f"{genome}_genes.gff")
    te_gff = os.path.join(prefix, "reference", f"{genome}_transposons.gff")

    # Process genes
    with open(gene_gff) as fh_in, \
         open("gene.gff", 'w') as fh_gene, \
         open("promoter.gff", 'w') as fh_prom:

        for line in fh_in:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')

            if 'gene' in cols[2]:
                m = __import__('re').search(r'ID=([A-Za-z0-9_.]+);', cols[8])
                if not m:
                    continue
                name = m.group(1)

                if 'Note=transposable_element_gene;' in cols[8]:
                    new_id = f"{name}_TEG"
                else:
                    new_id = name

                cols[8] = new_id
                fh_gene.write('\t'.join(cols) + '\n')

                # Promoter
                cols[8] = f"{new_id}_promoter"
                if cols[6] == '+':
                    cols[4] = str(int(cols[3]) - 1)
                    new_start = max(1, int(cols[3]) - promoter_length)
                    cols[3] = str(new_start)
                else:
                    cols[3] = str(int(cols[4]) + 1)
                    cols[4] = str(int(cols[4]) + promoter_length)

                fh_prom.write('\t'.join(cols) + '\n')

    # Process TEs
    if os.path.exists(te_gff):
        with open(te_gff) as fh_in, open("te.gff", 'w') as fh_te:
            for line in fh_in:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')

                if 'transposable_element' in cols[2]:
                    m = __import__('re').search(r'ID=([A-Za-z0-9_.]+);', cols[8])
                    if m:
                        cols[8] = m.group(1)
                        fh_te.write('\t'.join(cols) + '\n')


def build_annotation(prefix, genome, binsize=100, promoter_length=1000):
    """
    Build genome annotation by bins.
    Returns dict: {bin_id: annotation_string}
    """
    ann_path = os.path.join(prefix, "reference", f"{genome}.{binsize}.annotation")

    if os.path.exists(ann_path):
        ann = {}
        with open(ann_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')
                ann[cols[0]] = '\t'.join(cols[1:])
        return ann

    # Build from scratch
    split_gff(prefix, genome, promoter_length)

    ann = defaultdict(str)

    # Gene annotations
    with open("gene.gff") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            start_bin = int(cols[3]) // binsize
            end_bin = int(cols[4]) // binsize
            for i in range(start_bin, end_bin + 1):
                bin_id = f"{cols[0]}_{i}"
                ann[bin_id] += f"GENE:{cols[8]};"

    # TE annotations
    if os.path.exists("te.gff"):
        with open("te.gff") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')
                start_bin = int(cols[3]) // binsize
                end_bin = int(cols[4]) // binsize
                for i in range(start_bin, end_bin + 1):
                    bin_id = f"{cols[0]}_{i}"
                    ann[bin_id] += f"TE:{cols[8]};"

    # miRNA annotations
    mir_gff = os.path.join(prefix, "reference", f"{genome}_miRNA_miRNA_star.gff")
    if os.path.exists(mir_gff):
        with open(mir_gff) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')
                start_bin = int(cols[3]) // binsize
                end_bin = int(cols[4]) // binsize
                for i in range(start_bin, end_bin + 1):
                    bin_id = f"{cols[0]}_{i}"
                    ann[bin_id] += f"{cols[8]};"

    # Promoter annotations
    if os.path.exists("promoter.gff"):
        with open("promoter.gff") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')
                start_bin = int(cols[3]) // binsize
                end_bin = int(cols[4]) // binsize
                for i in range(start_bin, end_bin + 1):
                    bin_id = f"{cols[0]}_{i}"
                    ann[bin_id] += f"PROMOTER:{cols[8]};"

    # Clean up temp files
    for f in ("gene.gff", "te.gff", "promoter.gff"):
        if os.path.exists(f):
            os.unlink(f)

    # Strip trailing semicolons and write cache
    result = {}
    with open(ann_path, 'w') as fh_out:
        for bin_id in sorted(ann.keys()):
            ann_str = ann[bin_id].rstrip(';')
            result[bin_id] = ann_str
            fh_out.write(f"{bin_id}\t{ann_str}\n")

    return result


def read_chromosome_lengths(prefix, genome, binsize=100):
    """
    Read chromosome lengths from .fai file.
    Returns dict: {chromosome: num_bins}
    """
    fai_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta.fai")
    lengths = {}

    with open(fai_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            lengths[cols[0]] = int(cols[1]) // binsize

    return lengths


def read_gene_annotation(prefix, genome):
    """
    Read gene functional annotation (GO, Mapman).
    Returns dict: {gene_id: annotation_string}
    """
    gann = {}

    bin_path = os.path.join(prefix, "reference", f"{genome}.BIN")
    fun_path = os.path.join(prefix, "reference", f"{genome}.functional.annotation")

    if not os.path.exists(bin_path) or not os.path.exists(fun_path):
        return gann

    # Read BIN file
    with open(bin_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                gann[parts[0]] = parts[1] + ';'

    for gene in gann:
        gann[gene] = gann[gene].rstrip(';')

    # Read functional annotation
    with open(fun_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            for j in range(1, 5):
                if j < len(cols) and cols[j] == "":
                    cols[j] = "NA"

            m = __import__('re').match(r'(\w+)\.1$', cols[0])
            if m:
                gene_id = m.group(1)
                ann_str = f',"{cols[1]}","{cols[2]}","{cols[3]}","{cols[4]}"'
                if gene_id in gann:
                    gann[gene_id] += ann_str
                else:
                    gann[gene_id] = f"NA{ann_str}"

    return gann


def primary_transcript(prefix, genome):
    """
    Identify primary transcript (longest CDS) using getPrimaryTranscript.py.
    """
    script = os.path.join(prefix, "scripts", "getPrimaryTranscript.py")
    gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
    output_file = f"{genome}.PrimaryTranscript.txt"

    subprocess.run(
        f"python3 {script} {gff_path} > {output_file}",
        shell=True, check=True
    )

    # Parse primary transcript list
    primary = set()
    with open(output_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) >= 2:
                primary.add(cols[1])

    # Filter GTF to primary transcripts
    gtf_file = f"{genome}.gtf"
    out_file = f"{genome}.PrimaryTranscript.gtf"

    if not os.path.exists(gtf_file):
        return

    with open(gtf_file) as fh_in, open(out_file, 'w') as fh_out:
        for line in fh_in:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue

            m = __import__('re').search(r'transcript_id "(\w+\.\d+)"', cols[8])
            if m and m.group(1) in primary:
                fh_out.write(line + '\n')

    os.unlink(output_file)
