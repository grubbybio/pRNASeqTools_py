"""
Utility functions for pRNASeqTools.
Mirrors the Perl Function.pm module.
"""

import os
import re
import math
import shlex
import subprocess
import sys
import shutil


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


def run_cmd(cmd, tee=None, check=True, quiet=False, **kwargs):
    """Run a shell command, tee-ing stdout+stderr to log file and terminal.

    All subprocess output streams through the TeeHandler so alignment tool
    logs appear in both the pipeline log file and stderr in real-time.
    On failure the captured output is included in the exception so R / tool
    error messages are visible in the pipeline log.

    When quiet=True, command output is captured but NOT streamed to the
    terminal or log; it is only written to the tee if the command fails.
    """
    if tee is None:
        tee = _tee()
    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, **kwargs
    )
    out_lines = []
    if quiet:
        # Capture silently; only surface output on failure
        for line in proc.stdout:
            out_lines.append(line)
        remainder = proc.stdout.read()
        if remainder:
            out_lines.append(remainder)
        ret = proc.wait()
        if check and ret != 0:
            tail = "".join(out_lines[-40:]) if out_lines else "(no output)"
            tee.write(f"  Command failed: {cmd}\n{tail}\n")
            raise subprocess.CalledProcessError(
                ret, cmd, output=tail, stderr=None
            )
        return ret
    for line in proc.stdout:
        tee.write(line)
        out_lines.append(line)
    # Drain any trailing output buffered without a newline
    remainder = proc.stdout.read()
    if remainder:
        tee.write(remainder)
        out_lines.append(remainder)
    ret = proc.wait()
    if check and ret != 0:
        tail = "".join(out_lines[-40:]) if out_lines else "(no output)"
        raise subprocess.CalledProcessError(
            ret, cmd, output=tail, stderr=None
        )
    return ret


# ── SRA download ─────────────────────────────────────────────────────────
def download_sra(srr, threads=4):
    """
    Download SRA file if input is an SRR accession.
    Two-step: prefetch (parallel FTP download) → fasterq-dump (local conversion).
    Returns (file_paths_tuple, from_sra_bool).
      - from_sra=False: input was a local file or already-existing path
      - from_sra=True : input was an SRR accession and was downloaded via prefetch
    """
    import shutil
    m = re.search(r'([SED]RR\d+)', srr)
    if not m or os.path.exists(srr):
        return (srr,), False

    srr_id = m.group(1)
    tee = _tee()
    tee.write("Downloading...\n")

    # Step 1: prefetch — parallel FTP download of .sra file (much faster than fasterq-dump direct)
    # prefetch creates: .sra_cache/<SRRxxxx>/<SRRxxxx>.sra  (accession subdir)
    cache_dir = ".sra_cache"
    os.makedirs(cache_dir, exist_ok=True)
    run_cmd(f"prefetch -O {cache_dir} {shlex.quote(srr_id)}")

    # Step 2: locate .sra file — prefetch nests it under accession subdir
    sra_path = None
    for root, dirs, files in os.walk(cache_dir):
        for f in files:
            if f.endswith('.sra'):
                sra_path = os.path.join(root, f)
                break
        if sra_path:
            break
    if not sra_path:
        raise FileNotFoundError(f"prefetch failed: no .sra file found under {cache_dir}")

    # Step 3: fasterq-dump — local conversion from .sra to FASTQ
    # Use absolute path to avoid VDB relative-path resolution issues
    sra_abs = os.path.abspath(sra_path)
    run_cmd(f"fasterq-dump -p --threads {threads} --split-files {shlex.quote(sra_abs)}")

    # Step 4: clean up intermediate .sra and cache dir
    shutil.rmtree(cache_dir, ignore_errors=True)

    # Step 5: detect actual output files (fasterq-dump may skip 0-length reads)
    has_r1 = os.path.exists(f"{srr_id}_1.fastq")
    has_r2 = os.path.exists(f"{srr_id}_2.fastq")

    if has_r1 and has_r2:
        # True paired-end
        if os.path.exists(f"{srr_id}.fastq"):
            os.unlink(f"{srr_id}.fastq")
        return (f"{srr_id}_1.fastq", f"{srr_id}_2.fastq"), True
    elif has_r1 and not has_r2:
        # PE data but R2 was all 0-length / empty — R1 only
        tee.write("Warning: R2 reads are empty (0-length), treating as single-end\n")
        return (f"{srr_id}_1.fastq",), True
    else:
        # Single-end (or fallback to interleaved)
        if os.path.exists(f"{srr_id}.fastq"):
            return (f"{srr_id}.fastq",), True
        raise FileNotFoundError(
            f"fasterq-dump produced no FASTQ output for {srr_id}. "
            f"Expected {srr_id}.fastq or {srr_id}_1.fastq"
        )


def cleanup_alignment_outputs(extra_globs=()):
    """Remove per-tag aligner outputs (HISAT2 or STAR), fastp intermediates,
    and pipeline reports.

    Standard cleanup covers:
      - HISAT2: *.sam, *_Log.txt, *_SJ.out.tab
      - STAR (legacy): *_Log.*, *_Aligned.*, *_Unmapped.*
      - fastp trimmed outputs (*_trimmed.fastq / *.fq)
      - fastp reports (*.html / *.json)
      - Unzipped trimmed fastq ({tag}.fastq / {tag}_R1.fastq / {tag}_R2.fastq).
        Raw SRA unzipped files (*.raw.fastq) are skipped (see preserve_raw_fastq).
        Files matching *.raw.fastq / *.unmapped*.fastq are preserved.
    Pass extra_globs to add mode-specific patterns.
    Idempotent: missing files are silently skipped.
    """
    import glob as globmod
    patterns = [
        # HISAT2 / STAR common
        "*.sam", "*_Log.txt", "*_SJ.out.tab",
        # STAR legacy
        "*_Log.out", "*_Log.progress.out", "*_Log.final.out",
        "*_Aligned.sortedByCoord.out.bam", "*_Aligned.out.sam",
        "*_Unmapped.out.mate1", "*_Unmapped.out.mate2",
        "*_ReadsUnmapped*.fastq", "*_ReadsUnmapped*.fq",
        # fastp intermediates & reports
        "*_trimmed.fastq", "*_trimmed.fq",
        "*.fastq", "*.fq",
        "*.html", "*.json",
    ]
    patterns.extend(extra_globs)
    for pattern in patterns:
        for fname in globmod.glob(pattern):
            if not os.path.exists(fname):
                continue
            # Defensive: never delete *.raw.fastq / *.raw.fq / *.unmapped*.fastq
            # (preservation targets)
            if fname.endswith((".raw.fastq", ".raw.fq")) or ".unmapped" in fname:
                continue
            os.unlink(fname)


def preserve_raw_fastq(from_sra, tag, paired=False):
    """
    After unzip_file() produced {tag}.fastq (and optionally {tag}_R1.fastq / {tag}_R2.fastq),
    if the source was SRA/prefetch, move the raw unzipped file aside as {tag}.raw.fastq
    (and {tag}_R1.raw.fastq / {tag}_R2.raw.fastq), then re-create {tag}.fastq as a symlink
    pointing to the raw file. Downstream tools (fastp, STAR, bowtie2) follow the symlink
    transparently; later os.rename() onto {tag}.fastq replaces the symlink with the
    trimmed file, leaving the raw file intact for later gzip preservation.

    No-op when from_sra is False (local files are not preserved).
    """
    if not from_sra:
        return
    targets = [f"{tag}_R1.fastq", f"{tag}_R2.fastq"] if paired else [f"{tag}.fastq"]
    for t in targets:
        if not os.path.exists(t) or os.path.islink(t):
            continue
        raw = t.replace(".fastq", ".raw.fastq")
        os.rename(t, raw)
        os.symlink(raw, t)


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
        run_cmd(f"bzip2 -dc {shlex.quote(filepath)} > {shlex.quote(target)}")
    elif filepath.endswith('.gz'):
        tee.write("Decompressing...\n")
        run_cmd(f"gzip -dc {shlex.quote(filepath)} > {shlex.quote(target)}")
    elif filepath.endswith('.gtz'):
        tee.write("Decompressing...\n")
        run_cmd(f"gtz -dc {shlex.quote(filepath)} > {shlex.quote(target)}")
    elif filepath.endswith(('.fastq', '.fq')):
        if filepath != target:
            tee.write("Renaming...\n")
            os.rename(filepath, target)
        else:
            tee.write("Backing up...\n")
            run_cmd(f"cp {shlex.quote(filepath)} {shlex.quote(target)}.bak")
    else:
        raise ValueError("Please provide the seq file in correct formats!")

    tee.write("Completed!\n")


# ── Remove 3' PolyC and reverse complement ───────────────────────────────


# ── Save raw fastq as gzip instead of deleting ────────────────────────────────
def gzip_fastq(filepath):
    """
    Compress a FASTQ/FASTQ.gz file in-place using gzip.
    Safe: no-op if the file doesn't exist.
    """
    if not os.path.exists(filepath):
        return
    tee = _tee()
    tee.write(f"Saving compressed: {filepath}.gz\n")
    run_cmd(f"gzip -f {shlex.quote(filepath)}")


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
                    rmvc_dict[seqname]['outq'] = line[:-clen][::-1] if clen else line
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


# ── HISAT2 shared index helpers ──────────────────────────────────────────
# HISAT2 索引是一组 HisatIdx.1.ht2 ~ HisatIdx.N.ht2 的文件组，
# 用 '.1.ht2' 作为完整性检查标记。

def _hisat_index_prefix(hisat_idx):
    """Normalize HISAT2 index prefix — strip trailing markers if present."""
    for suffix in ('.1.ht2', '.ht2', '.htl'):
        if hisat_idx.endswith(suffix):
            return hisat_idx[:-len(suffix)]
    return hisat_idx


def _is_hisat_index_complete(hisat_idx):
    """Check if HISAT2 index prefix has all files (.1.ht2 or .htl)."""
    p = _hisat_index_prefix(hisat_idx)
    return os.path.exists(p + '.1.ht2') or os.path.exists(p + '.htl')


def _symlink_hisat_index_files(ref_idx, local_idx):
    """Symlink all HISAT2 index files from ref_idx prefix → local_idx prefix."""
    import glob as _glob
    ref_prefix = _hisat_index_prefix(ref_idx)
    local_prefix = _hisat_index_prefix(local_idx)
    for pattern in (ref_prefix + '.*.ht2', ref_prefix + '.htl'):
        for rf in _glob.glob(pattern):
            lf = rf.replace(ref_prefix, local_prefix)
            if not os.path.exists(lf):
                os.symlink(rf, lf)


def try_use_shared_hisat_index(genome, prefix, local_idx='HisatIdx', tee=None):
    """
    Check if a pre-built HISAT2 index exists at reference/{genome}_HisatIdx.
    If yes → symlink local_idx → reference files, return True.
    If no  → return False (caller should run hisat2-build normally).
    """
    if tee is None:
        tee = sys.stderr

    # Case 1: local index files already present
    if _is_hisat_index_complete(local_idx):
        tee.write(f"  复用本地 HISAT2 index ({local_idx})\n")
        return True

    # Case 2: shared index exists in reference/ → symlink
    ref_idx = os.path.join(prefix, 'reference', f'{genome}_HisatIdx')
    if _is_hisat_index_complete(ref_idx):
        tee.write(f"  发现共享 HISAT2 index: {ref_idx}\n")
        _symlink_hisat_index_files(ref_idx, local_idx)
        tee.write(f"  已创建 symlink 到共享 index\n")
        return True

    return False


def save_hisat_index_to_reference(genome, prefix, local_idx='HisatIdx', tee=None):
    """After hisat2-build, move local_idx files → reference/{genome}_HisatIdx
    and create symlinks back. Future runs reuse the shared index."""
    if tee is None:
        tee = sys.stderr

    if not _is_hisat_index_complete(local_idx):
        tee.write(f"  save_hisat_index: {local_idx} 不完整，跳过\n")
        return

    ref_idx = os.path.join(prefix, 'reference', f'{genome}_HisatIdx')
    import glob as _glob
    local_prefix = _hisat_index_prefix(local_idx)
    local_files = _glob.glob(local_prefix + '.*.ht2') + _glob.glob(local_prefix + '.htl')

    if _is_hisat_index_complete(ref_idx):
        tee.write(f"  reference/{genome}_HisatIdx 已存在 → 保留本地副本\n")
        return

    os.makedirs(os.path.dirname(ref_idx), exist_ok=True)
    try:
        for lf in local_files:
            rf = lf.replace(local_prefix, ref_idx)
            shutil.move(lf, rf)
            os.symlink(rf, lf)
        tee.write(f"  已保存共享 HISAT2 index: {ref_idx}\n")
    except Exception as e:
        tee.write(f"  保存到 reference 失败（{e}）→ 保留本地 index\n")


# ── HISAT2 alignment helpers ──────────────────────────────────────────────
def hisat2_build_index(fasta_path, gtf_path, hisat_idx, thread, tee=None):
    """Build HISAT2 index with known splice sites from GTF.

    Args:
        fasta_path: genome FASTA path
        gtf_path: annotation GTF (or GFF3, converted by extract_splice_sites.py)
        hisat_idx: output index prefix
        thread: thread count
    """
    if tee is None:
        tee = sys.stderr

    splicesites = f"_{hisat_idx}_splicesites.txt"
    # HISAT2 自带 extract_splice_sites.py，或 bioconda 的 HISAT2 包也带
    # 也可从 GFF 手动提取 GT/AG 边界——但用自带脚本最稳
    try:
        run_cmd(
            f"extract_splice_sites.py {gtf_path} > {splicesites}",
            tee=tee)
    except Exception:
        # fallback: 不带 splice sites 也能跑，只是比对率可能略低
        tee.write("  Warning: extract_splice_sites.py 不可用，跳过已知 splice sites\n")
        splicesites = None

    hisat2_build_cmd = f"hisat2-build -p {thread}"
    if splicesites:
        hisat2_build_cmd += f" --ss {splicesites}"
    hisat2_build_cmd += f" {fasta_path} {hisat_idx}"
    run_cmd(hisat2_build_cmd, tee=tee)

    # 清理临时 splice sites 文件
    if splicesites and os.path.exists(splicesites):
        os.unlink(splicesites)


def hisat2_align_single(tag, fastq, hisat_idx, thread, extra_args='',
                        nov_sj_file=None, known_sj_file=None, out_dir=None):
    """HISAT2 single-end alignment → sorted BAM {tag}.bam.

    Args:
        out_dir:  if given, BAM/log files are written here (creates subdir).
    """
    out_prefix = tag if out_dir is None else os.path.join(out_dir, tag)
    nov_sj_opt = f"--novel-splicesite-outfile {out_prefix}_SJ.out.tab" if nov_sj_file else ""
    known_sj_opt = f"--known-splicesite-infile {known_sj_file}" if known_sj_file else ""
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    run_cmd(
        f"hisat2 -p {thread} -x {hisat_idx} -U {fastq} "
        f"-S {out_prefix}.sam {nov_sj_opt} {known_sj_opt} {extra_args} "
        f"2>{out_prefix}_Log.txt")
    run_cmd(f"samtools sort -@ {thread} -o {out_prefix}.bam {out_prefix}.sam")
    run_cmd(f"samtools index {out_prefix}.bam")
    os.unlink(f"{out_prefix}.sam")


def hisat2_align_paired(tag, fq1, fq2, hisat_idx, thread, extra_args='',
                        nov_sj_file=None, known_sj_file=None, out_dir=None):
    """HISAT2 paired-end alignment → sorted BAM {tag}.bam."""
    out_prefix = tag if out_dir is None else os.path.join(out_dir, tag)
    nov_sj_opt = f"--novel-splicesite-outfile {out_prefix}_SJ.out.tab" if nov_sj_file else ""
    known_sj_opt = f"--known-splicesite-infile {known_sj_file}" if known_sj_file else ""
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    run_cmd(
        f"hisat2 -p {thread} -x {hisat_idx} -1 {fq1} -2 {fq2} "
        f"-S {out_prefix}.sam {nov_sj_opt} {known_sj_opt} {extra_args} "
        f"2>{out_prefix}_Log.txt")
    run_cmd(f"samtools sort -@ {thread} -o {out_prefix}.bam {out_prefix}.sam")
    run_cmd(f"samtools index {out_prefix}.bam")
    os.unlink(f"{out_prefix}.sam")
