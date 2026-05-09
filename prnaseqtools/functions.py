"""
Utility functions for pRNASeqTools.
Mirrors the Perl Function.pm module.
"""

import os
import re
import math
import subprocess
import sys


# ── Helper to get the TEE output handle ──────────────────────────────────
def _tee():
    """Return the TeeHandler for writing to log + stderr."""
    try:
        from prnaseqtools.cli import TEE
        if TEE:
            return TEE
    except (ImportError, AttributeError):
        pass
    return sys.stderr


# ── SRA download ─────────────────────────────────────────────────────────
def download_sra(srr, threads=4):
    """
    Download SRA file if input is an SRR accession.
    Returns file path(s): single-end returns (file,), paired-end returns (r1, r2).
    """
    m = re.search(r'([SED]RR\d+)', srr)
    if not m or os.path.exists(srr):
        return (srr,)

    srr_id = m.group(1)
    tee = _tee()
    tee.write("Downloading...\n")

    subprocess.run(
        f"fasterq-dump -p --threads {threads} --split-3 {srr_id}",
        shell=True, check=True
    )

    if os.path.exists(f"{srr_id}_1.fastq"):
        if os.path.exists(f"{srr_id}.fastq"):
            os.unlink(f"{srr_id}.fastq")
        return (f"{srr_id}_1.fastq", f"{srr_id}_2.fastq")
    else:
        return (f"{srr_id}.fastq",)


# ── Statistics ───────────────────────────────────────────────────────────
def sum_values(values):
    """Sum of a list of numbers."""
    return sum(values)


def average(values):
    """Arithmetic mean."""
    if not values:
        return 0
    return sum(values) / len(values)


def stdev(values):
    """Sample standard deviation."""
    if len(values) < 2:
        return 0
    avg = average(values)
    sqtotal = sum((avg - v) ** 2 for v in values)
    return math.sqrt(sqtotal / (len(values) - 1))


# ── Reverse complement ───────────────────────────────────────────────────
COMPLEMENT_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(seq):
    """Reverse complement of a DNA sequence."""
    return seq[::-1].translate(COMPLEMENT_TABLE)


# ── File decompression / unzip ───────────────────────────────────────────
def unzip_file(filepath, tag):
    """
    Decompress or copy sequence file to <tag>.fastq.
    Supports .bz2, .gz, .gtz, .fastq, .fq.
    """
    tee = _tee()

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Please provide the correct seq file: {filepath}")

    target = f"{tag}.fastq"

    if filepath.endswith('.bz2'):
        tee.write("Decompressing...\n")
        subprocess.run(f"bzip2 -dc {filepath} > {target}", shell=True, check=True)
    elif filepath.endswith('.gz'):
        tee.write("Decompressing...\n")
        subprocess.run(f"gzip -dc {filepath} > {target}", shell=True, check=True)
    elif filepath.endswith('.gtz'):
        tee.write("Decompressing...\n")
        subprocess.run(f"gtz -dc {filepath} > {target}", shell=True, check=True)
    elif filepath.endswith(('.fastq', '.fq')):
        if filepath != target:
            tee.write("Renaming...\n")
            subprocess.run(f"cp {filepath} {target}", shell=True, check=True)
            if filepath.startswith('SRR') or filepath.startswith('ERR') or filepath.startswith('DRR'):
                os.unlink(filepath)
        else:
            tee.write("Backing up...\n")
            subprocess.run(f"cp {filepath} {target}.bak", shell=True, check=True)
    else:
        raise ValueError("Please provide the seq file in correct formats!")

    tee.write("Completed!\n")


# ── Remove 3' PolyC and reverse complement ───────────────────────────────
def rmvc(tag_r1, tag_r2=None):
    """
    Remove 3' PolyC from trimmed reads and reverse complement.
    Used by TS-CLIP-seq mode.
    """
    tee = _tee()
    tee.write("Start to remove the 3' PolyC...\n")

    trimmed_file = f"{tag_r1}_trimmed.fastq"
    if not os.path.exists(trimmed_file):
        tee.write(f"Warning: {trimmed_file} not found, skipping PolyC removal\n")
        return

    # Process R1
    rmvc_dict = {}
    with open(trimmed_file) as fh:
        seqname = None
        for i, line in enumerate(fh):
            line = line.strip()
            n = i % 4
            if n == 3:  # quality
                if seqname and seqname in rmvc_dict:
                    clen = rmvc_dict[seqname].get('clen', 0)
                    rmvc_dict[seqname]['outq'] = line[-clen:][::-1] if clen else line
            elif n == 1:  # sequence
                m = re.match(r'^([A-Z]+[ATG])(C+)$', line)
                if m:
                    rmvc_dict[seqname] = {
                        'outseq': revcomp(m.group(1)),
                        'clen': len(m.group(2))
                    }
            elif n == 0:  # header
                parts = line.split()
                seqname = parts[0]

    # Process R2 if paired
    if tag_r2:
        rmvc_dict2 = {}
        trimmed_file2 = f"{tag_r2}_trimmed.fastq"
        with open(trimmed_file2) as fh:
            seqname = None
            for i, line in enumerate(fh):
                line = line.strip()
                n = i % 4
                if n == 3:
                    if seqname and seqname in rmvc_dict2:
                        clen = rmvc_dict2[seqname].get('clen', 0)
                        rmvc_dict2[seqname]['outq'] = line[clen:][::-1] if clen else line
                elif n == 1:
                    m = re.match(r'^(G+)([ATC][A-Z]+)$', line)
                    if m:
                        rmvc_dict2[seqname] = {
                            'outseq': revcomp(m.group(2)),
                            'clen': len(m.group(1))
                        }
                elif n == 0:
                    parts = line.split()
                    seqname = parts[0]

    # Write output
    with open(f"{tag_r1}.fastq", 'w') as out1:
        out2 = open(f"{tag_r2}.fastq", 'w') if tag_r2 else None
        try:
            for seqname, data in rmvc_dict.items():
                if 'outq' not in data:
                    continue
                if tag_r2:
                    if seqname in rmvc_dict2 and 'outq' in rmvc_dict2[seqname]:
                        out1.write(f"{seqname}\n{data['outseq']}\n+\n{data['outq']}\n")
                        out2.write(f"{seqname}\n{rmvc_dict2[seqname]['outseq']}\n+\n{rmvc_dict2[seqname]['outq']}\n")
                else:
                    out1.write(f"{seqname}\n{data['outseq']}\n+\n{data['outq']}\n")
        finally:
            if out2:
                out2.close()

    tee.write("The 3' PolyC is removed!\n")


# ── Fastq to fasta ───────────────────────────────────────────────────────
def fastq2fasta(input_file, output_file):
    """Convert FASTQ to FASTA format."""
    with open(input_file) as fin, open(output_file, 'w') as fout:
        for i, line in enumerate(fin):
            n = i % 4
            if n == 0:
                fout.write(f">{line[1:]}")
            elif n == 1:
                fout.write(line)
