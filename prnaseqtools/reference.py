"""
Reference genome handling for pRNASeqTools.
Mirrors the Perl Ref.pm module.
Provides: GFF parsing, FASTA reading, exon extraction, annotation, etc.
"""

import os
import re
import subprocess
import tempfile
from collections import defaultdict


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
                m = re.search(r'^ID=([^;]+);', cols[8])
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
    if not os.path.exists(fas_path):
        raise FileNotFoundError(
            f"Reference FASTA not found: {fas_path}\n"
            f"Please place the genome file in reference/{genome}_chr_all.fasta"
        )
    fas = {}
    current_name = None
    chunks = []

    with open(fas_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    fas[current_name] = "".join(chunks)
                    chunks = []
                current_name = line[1:].split()[0]
            else:
                chunks.append(line)
        if current_name is not None:
            fas[current_name] = "".join(chunks)

    return fas


def read_exons(prefix, genome):
    """
    Extract exons using gffread, return dict: {gene: {transcript: sequence}}.
    """
    exons_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
    gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")

    with tempfile.TemporaryDirectory() as tmpdir:
        exons_fa = os.path.join(tmpdir, "exons.fa")
        subprocess.run(
            f"gffread -O -w {exons_fa} -g {exons_path} {gff_path}",
            shell=True, check=True, timeout=600
        )

        exon_data = defaultdict(lambda: defaultdict(str))
        gene = None
        tran = 0

        with open(exons_fa) as fh:
            chunks = []
            for line in fh:
                line = line.strip()
                if line.startswith('>'):
                    if gene is not None:
                        exon_data[gene][tran] = "".join(chunks)
                        chunks = []
                    header = line[1:]
                    m = re.match(r'(\w+)\.(\d+)', header)
                    if m:
                        gene = m.group(1)
                        tran = m.group(2)
                    else:
                        gene = header
                        tran = 0
                else:
                    chunks.append(line)
            if gene is not None:
                exon_data[gene][tran] = "".join(chunks)

    return dict(exon_data)


def get_gene_info(prefix, genome):
    """
    Get gene length info (longest transcript per gene).
    Writes 'transcripts.fa' in CWD for downstream pipeline use.
    Returns dict: {gene: length}
    """
    exon_data = read_exons(prefix, genome)
    gene_info = {}

    try:
        fh = open("transcripts.fa", 'w')
    except IOError:
        raise IOError(
            "Cannot write transcripts.fa to current directory. "
            "Ensure you have write permission in the working directory."
        )
    with fh as tra_fh:
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
    if not os.path.exists(gff_path):
        raise FileNotFoundError(
            f"miRNA GFF not found: {gff_path}\n"
            f"Please place the miRNA annotation in reference/{genome}_miRNA_miRNA_star.gff"
        )
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


def split_gff(prefix, genome, promoter_length=1000, workdir=None):
    """
    Split GFF into gene.gff, te.gff, promoter.gff.

    Args:
        workdir: output directory (None = current working directory)
    """
    gene_gff = os.path.join(prefix, "reference", f"{genome}_genes.gff")
    te_gff = os.path.join(prefix, "reference", f"{genome}_transposons.gff")

    def _path(name):
        return os.path.join(workdir, name) if workdir else name

    # Process genes
    with open(gene_gff) as fh_in, \
         open(_path("gene.gff"), 'w') as fh_gene, \
         open(_path("promoter.gff"), 'w') as fh_prom:

        for line in fh_in:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')

            if 'gene' in cols[2]:
                m = re.search(r'ID=([A-Za-z0-9_.]+);', cols[8])
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
        with open(te_gff) as fh_in, open(_path("te.gff"), 'w') as fh_te:
            for line in fh_in:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')

                if 'transposable_element' in cols[2]:
                    m = re.search(r'ID=([A-Za-z0-9_.]+);', cols[8])
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

    # Build from scratch (temp files managed by TemporaryDirectory)
    with tempfile.TemporaryDirectory() as tmpdir:
        split_gff(prefix, genome, promoter_length, workdir=tmpdir)

        ann = defaultdict(str)

        # Gene annotations
        with open(os.path.join(tmpdir, "gene.gff")) as fh:
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
        te_gff = os.path.join(tmpdir, "te.gff")
        if os.path.exists(te_gff):
            with open(te_gff) as fh:
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
        prom_gff = os.path.join(tmpdir, "promoter.gff")
        if os.path.exists(prom_gff):
            with open(prom_gff) as fh:
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

        # Write cache (outside tmpdir so it persists)
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
    if not os.path.exists(fai_path):
        # Try auto-building .fai index
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(
                f"Reference FASTA not found: {fasta_path}"
            )
        subprocess.run(
            f"samtools faidx {fasta_path}",
            shell=True, check=True, timeout=300
        )
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

            m = re.match(r'(\w+)\.1$', cols[0])
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

    with tempfile.TemporaryDirectory() as tmpdir:
        output_file = os.path.join(tmpdir, f"{genome}.PrimaryTranscript.txt")
        subprocess.run(
            f"python3 {script} {gff_path} > {output_file}",
            shell=True, check=True, timeout=600
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

            m = re.search(r'transcript_id "(\w+\.\d+)"', cols[8])
            if m and m.group(1) in primary:
                fh_out.write(line + '\n')


# ═══════════════════════════════════════════════════════════════════════════
# Mapping index management — auto-detect and build missing indices
# ═══════════════════════════════════════════════════════════════════════════

INDEX_BUILDERS = {
    'fai': {
        'name': 'Samtools FAI index',
        'build_cmd': 'samtools faidx {fasta}',
        'check_file': '{fasta}.fai',
        'location': 'reference',
        'modes': None,  # all modes
    },
    'lsu_rrna': {
        'name': 'LSU rRNA Bowtie1 index',
        'build_cmd': 'bowtie-build -q {prefix}/reference/{genome}_lsu_rrna.fasta '
                     '{prefix}/reference/lsu_rrna',
        'check_file': '{prefix}/reference/lsu_rrna.1.ebwt',
        'location': 'reference',
        'modes': ['srna', 'phasi'],
    },
    'bowtie1_genome': {
        'name': 'Bowtie1 genome index',
        'build_cmd': 'bowtie-build -q {prefix}/reference/{genome}_chr_all.fasta '
                     '{prefix}/reference/{genome}_chr_all',
        'check_file': '{prefix}/reference/{genome}_chr_all.1.ebwt',
        'location': 'reference',
        'modes': ['tt'],
    },
}


def _detect_genomes(prefix):
    """Discover available genomes in the reference directory.
    Returns list of genome names (e.g. ['ath', 'osa'])."""
    ref_dir = os.path.join(prefix, "reference")
    genomes = set()
    if not os.path.isdir(ref_dir):
        return []
    pat = re.compile(r'^(.+)_chr_all\.fasta$')
    for fname in os.listdir(ref_dir):
        m = pat.match(fname)
        if m:
            genomes.add(m.group(1))
    return sorted(genomes)


def check_and_build_indices(prefix, genome=None, mode=None, tee=None):
    """
    Check and auto-build missing mapping indices.

    Args:
        prefix: project root directory
        genome: genome name (e.g. 'ath'); None = auto-detect all genomes in reference/
        mode:   analysis mode (e.g. 'srna', 'chip'); None = check all indices
        tee:    output handle (defaults to stderr)

    Returns:
        dict: {genome: {index_name: bool}} or {index_name: bool} if single genome
    """
    if tee is None:
        import sys
        tee = sys.stderr

    if genome:
        genomes = [genome]
    else:
        genomes = _detect_genomes(prefix)

    if not genomes:
        tee.write("  No reference genomes found in reference/ (expected *_chr_all.fasta)\n")
        return {}

    tee.write(f"\nChecking mapping indices for: {', '.join(genomes)}\n")

    all_results = {}
    for g in genomes:
        ref_dir = os.path.join(prefix, "reference")
        fasta = os.path.join(ref_dir, f"{g}_chr_all.fasta")

        if not os.path.exists(fasta):
            tee.write(f"  {g}: genome FASTA not found, skipping index build\n")
            continue

        results = {}
        for key, cfg in INDEX_BUILDERS.items():
            mode_list = cfg.get('modes')
            if mode and mode_list and mode not in mode_list:
                continue

            check_file = cfg['check_file'].format(prefix=prefix, genome=g,
                                                   fasta=fasta)
            if os.path.exists(check_file):
                tee.write(f"  [{g}] {cfg['name']}: found\n")
                results[key] = True
                continue

            tee.write(f"  [{g}] {cfg['name']}: not found, building...\n")
            try:
                cmd = cfg['build_cmd'].format(prefix=prefix, genome=g,
                                               fasta=fasta)
                subprocess.run(cmd, shell=True, check=True, timeout=3600)
                tee.write(f"  [{g}] {cfg['name']}: built successfully\n")
                results[key] = True
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                tee.write(f"  [{g}] {cfg['name']}: build failed ({e})\n")
                results[key] = False

        all_results[g] = results

    # Flatten if single genome
    if genome and genome in all_results:
        return all_results[genome]
    return all_results
