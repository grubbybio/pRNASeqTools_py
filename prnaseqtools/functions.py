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


def run_cmd(cmd, tee=None, check=True, **kwargs):
    """Run a shell command, tee-ing stdout+stderr to log file and terminal.

    All subprocess output streams through the TeeHandler so alignment tool
    logs appear in both the pipeline log file and stderr in real-time.
    """
    if tee is None:
        tee = _tee()
    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, **kwargs
    )
    for line in proc.stdout:
        tee.write(line)
    ret = proc.wait()
    if check and ret != 0:
        raise subprocess.CalledProcessError(ret, cmd)
    return ret


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

    run_cmd(
        f"fasterq-dump -p --threads {threads} --split-3 {srr_id}"
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
        run_cmd(f"bzip2 -dc {filepath} > {target}")
    elif filepath.endswith('.gz'):
        tee.write("Decompressing...\n")
        run_cmd(f"gzip -dc {filepath} > {target}")
    elif filepath.endswith('.gtz'):
        tee.write("Decompressing...\n")
        run_cmd(f"gtz -dc {filepath} > {target}")
    elif filepath.endswith(('.fastq', '.fq')):
        if filepath != target:
            tee.write("Renaming...\n")
            run_cmd(f"cp {filepath} {target}")
            if filepath.startswith('SRR') or filepath.startswith('ERR') or filepath.startswith('DRR'):
                os.unlink(filepath)
        else:
            tee.write("Backing up...\n")
            run_cmd(f"cp {filepath} {target}.bak")
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


# ── ShortStack condensed BAM support ─────────────────────────────────────
def bam_is_condensed(bam_path):
    """Check if BAM uses ShortStack >= 4.1 condensed format (XW:i:n tags)."""
    result = subprocess.run(
        f"samtools view {bam_path} | head -5",
        shell=True, capture_output=True, text=True
    )
    return 'XW:i:' in result.stdout


def expand_bed_by_xw(bed_path, bam_path):
    """Expand a BED file by duplicating lines according to XW:i:n tags in BAM.
    Streams SAM → expands each BED line by its own XW value.
    Overwrites the original BED file."""
    run_cmd(
        f"samtools view {bam_path} | "
        f"awk '{{xw=1; for(i=12;i<=NF;i++) "
        f"if($i ~ /^XW:i:/) {{xw=substr($i,6)+0; break}}; "
        f"print xw}}' > {bed_path}.xw && "
        f"paste {bed_path} {bed_path}.xw | "
        f"awk -F '\\t' '{{xw=$(NF); "
        f"for(i=1;i<=xw;i++) {{for(j=1;j<NF;j++) "
        f"printf \"%s%s\", $j, (j<NF-1?\"\\t\":\"\"); print \"\"}}}}' "
        f"> {bed_path}.tmp && "
        f"mv {bed_path}.tmp {bed_path} && "
        f"rm {bed_path}.xw"
    )


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
