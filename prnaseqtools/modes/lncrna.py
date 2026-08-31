"""
lncRNA-seq analysis mode.

Pipeline:
  1. STAR mapping (same as mrna mode)
  2. StringTie transcript assembly → per-sample GTF
  3. Merge all sample GTFs with reference (gffcompare) → full transcriptome
  4. Classify novel transcripts:
       - overlap reference → coding or known ncRNA
       - no overlap (intergenic / anti-sense / intronic) → novel
       - FEELnc (bioconda) / PLEK2 (deep-learning) to distinguish coding vs lncRNA
  5. Quantify featureCounts on full merged GTF
  6. DESeq2 differential expression, split by transcript class
  7. Output summary tables + visualization

Inputs: FASTQ / BAM (existing STAR-aligned BAM supported).
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.input_parser import (parse_input, _parse_to_dict,
                                        _resolve_path)
from prnaseqtools.functions import (_tee, run_cmd, download_sra, unzip_file)


def run(opts):
    """Main entry point for lncRNA-seq analysis."""
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    foldchange = opts.get('foldchange', 2.0)
    pvalue = opts.get('pvalue', 0.01)
    fdr = opts.get('fdr', 1.0)
    run_mode = opts.get('run_mode', 'whole')
    control = opts.get('control', '')
    treatment = opts.get('treatment')
    mask = opts.get('mask')
    seq_strategy = opts.get('seq_strategy')
    genome_size = opts.get('genome_size', 10)
    classifier = opts.get('classifier', 'feelnc')   # feelnc | plek2 | fallback
    feelnc_dir = opts.get('feelnc_dir', None)
    plek2_dir = opts.get('plek2_dir',
                         os.path.expanduser('~/PLEK2'))
    min_fpkm = opts.get('min_fpkm', 0.5)            # min FPKM per transcript
    min_tx_len = opts.get('min_tx_len', 200)         # min transcript length

    mapping = run_mode in ('whole', 'mapping-only')
    do_assemble = run_mode in ('whole', 'mapping-only', 'assemble-only')
    do_quantify = run_mode != 'mapping-only'
    do_de = run_mode in ('whole', 'de-only', 'count-table')

    # ── Reference paths ─────────────────────────────────────────────────
    ref_dir = os.path.join(prefix, "reference")
    ref_gff = os.path.join(ref_dir, f"{genome}_genes.gff")
    ref_fasta = os.path.join(ref_dir, f"{genome}_chr_all.fasta")
    ref_gtf = f"{genome}_genes.gtf"

    # ── Parse inputs ────────────────────────────────────────────────────
    control_dict = _parse_to_dict(control)
    tags, files, pars = parse_input(control_dict)

    if treatment:
        for t in (treatment if isinstance(treatment, list) else [treatment]):
            treatment_dict = _parse_to_dict(t)
            t_tags, t_files, t_pars = parse_input(treatment_dict)
            tags.extend(t_tags)
            files.extend(t_files)
            pars.extend(t_pars)

    par_str = ' '.join(pars)

    # ── Stage 1: Mapping ────────────────────────────────────────────────
    # Reuse mrna-style STAR mapping for each sample.
    # After this stage we expect {tag}.bam for each sample.
    if mapping:
        _do_mapping(tags, files, pars, seq_strategy, adaptor, mask,
                    thread, genome_size, ref_fasta, ref_gff, prefix, tee)

    # ── Stage 2: StringTie transcript assembly ──────────────────────────
    if do_assemble:
        _do_assembly(tags, thread, ref_gtf, prefix, genome,
                     min_fpkm, tee)

        # Merge all per-sample GTFs + reference
        merged_gtf = _do_merge(tags, ref_gtf, thread, tee)

        # gffcompare vs reference
        class_gtf = _do_gffcompare(merged_gtf, ref_gtf, tee)

        # Classify novel transcripts (coding vs lncRNA)
        full_gtf = _do_classification(class_gtf, ref_gtf, ref_fasta,
                                       classifier, feelnc_dir, plek2_dir,
                                       min_tx_len, thread, tee)

    # ── Stage 3: Quantification ─────────────────────────────────────────
    if do_quantify:
        if do_assemble:
            full_gtf = "transcripts_final.gtf"
        else:
            # User provides pre-built merged GTF via --gtf
            full_gtf = opts.get('gtf')
            if not full_gtf:
                # Fallback: reference only
                full_gtf = ref_gtf
                tee.write(
                    f"Warning: --gtf not provided and no assembly done; "
                    f"quantifying against reference GTF only ({full_gtf})\n"
                )
            full_gtf = _resolve_path(full_gtf)

        # Ensure ref GTF exists
        if not os.path.exists(ref_gtf):
            run_cmd(f"gffread -T -C -o {ref_gtf} -g {ref_fasta} {ref_gff}")

        for tag in tags:
            if mapping or do_assemble:
                # BAM from mapping stage
                bam = f"{tag}.bam"
            else:
                # User-supplied BAM
                bam = os.path.basename(tag) + ".bam"  # assume symlinked
                if not os.path.exists(bam):
                    os.symlink(f"../{tag}.bam", bam)

            _quantify_one(tag, bam, thread, full_gtf, ref_fasta,
                          seq_strategy, tee)

    # ── Stage 4: DE analysis (DESeq2, split by tx class) ────────────────
    if do_de:
        if not os.path.exists(ref_gtf):
            run_cmd(f"gffread -T -C -o {ref_gtf} -g {ref_fasta} {ref_gff}")

        full_gtf = "transcripts_final.gtf"
        if not os.path.exists(full_gtf):
            # Fall back to ref-only GTF
            tee.write(
                f"Warning: {full_gtf} not found; using {ref_gtf} for DE\n"
            )
            full_gtf = ref_gtf

        tee.write(
            f"\nDifferential expression (FC≥{foldchange}, "
            f"P<{pvalue}, FDR<{fdr})\n"
        )

        script_path = os.path.join(prefix, "scripts", "lncrna.R")
        cmd = (
            f"Rscript --vanilla {script_path} "
            f"{foldchange} {pvalue} {fdr} {full_gtf} {par_str}"
        )
        tee.write(f"Running: {cmd}\n")
        run_cmd(cmd)


# ═══════════════════════════════════════════════════════════════════════
# Stage helpers
# ═══════════════════════════════════════════════════════════════════════

def _do_mapping(tags, files, pars, seq_strategy, adaptor, mask,
                thread, genome_size, ref_fasta, ref_gff, prefix, tee):
    """Run STAR mapping for every sample. Mirrors mrna._run mapping block."""

    tee.write("\nBuilding STAR genome index ...\n")
    if os.path.exists("Genome"):
        run_cmd("rm -rf Genome")
    os.makedirs("Genome", exist_ok=True)

    run_cmd(
        f"STAR --runThreadN {thread} --genomeDir Genome --runMode genomeGenerate "
        f"--genomeSAindexNbases {genome_size} --genomeFastaFiles {ref_fasta} "
        f"--sjdbGTFfile {ref_gff} --sjdbGTFtagExonParentTranscript Parent "
        f"--sjdbGTFtagExonParentGene ID --limitGenomeGenerateRAM 64000000000")

    if mask:
        mask_path = _resolve_path(mask)
        os.symlink(mask_path, "mask.fa")
        run_cmd("bowtie-build -q mask.fa mask")

    for i, tag in enumerate(tags):
        fpath = files[i]
        tee.write(f"\nMapping {tag}...\n")

        # ── Resolve fastqs ──
        if ',' not in fpath:
            # Single: could be SRA or local fastq
            sra = download_sra(fpath, thread)
            if len(sra) == 1:
                seq_strategy = 'single'
                unzip_file(sra[0], tag)
                if adaptor:
                    run_cmd(
                        f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} "
                        f"-o {tag}_trimmed.fastq {tag}.fastq")
                    os.rename(f"{tag}_trimmed.fastq", f"{tag}.fastq")
                if mask:
                    run_cmd(
                        f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                        f"{tag}.fastq {tag}.mask.out")
                    os.rename("tmp.fastq", f"{tag}.fastq")

                run_cmd(
                    f"STAR --runMode alignReads --genomeDir Genome "
                    f"--alignIntronMax 5000 --outReadsUnmapped Fastx "
                    f"--outSAMtype BAM SortedByCoordinate "
                    f"--limitBAMsortRAM 10000000000 --outSAMmultNmax 1 "
                    f"--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 "
                    f"--runThreadN {thread} --readFilesIn {tag}.fastq")
                if os.path.exists(f"{tag}.fastq"):
                    os.unlink(f"{tag}.fastq")
            else:
                # Paired SRA
                seq_strategy = 'paired'
                unzip_file(sra[0], f"{tag}_R1")
                unzip_file(sra[1], f"{tag}_R2")
                if adaptor:
                    run_cmd(
                        f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} "
                        f"-A {adaptor} -o {tag}_R1_trimmed.fastq "
                        f"-p {tag}_R2_trimmed.fastq "
                        f"{tag}_R1.fastq {tag}_R2.fastq")
                    os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                    os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
                if mask:
                    run_cmd(
                        f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                        f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq {tag}.mask.out")
                    if os.path.exists("tmp_1.fastq"):
                        os.rename("tmp_1.fastq", f"{tag}_R1.fastq")
                    if os.path.exists("tmp_2.fastq"):
                        os.rename("tmp_2.fastq", f"{tag}_R2.fastq")

                run_cmd(
                    f"STAR --runMode alignReads --genomeDir Genome "
                    f"--alignIntronMax 5000 --outReadsUnmapped Fastx "
                    f"--outSAMtype BAM SortedByCoordinate "
                    f"--limitBAMsortRAM 10000000000 --outSAMmultNmax 1 "
                    f"--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 "
                    f"--runThreadN {thread} --readFilesIn "
                    f"{tag}_R1.fastq {tag}_R2.fastq")
                for fn in (f"{tag}_R1.fastq", f"{tag}_R2.fastq"):
                    if os.path.exists(fn):
                        os.unlink(fn)
        else:
            # Explicit local paired-end
            f1, f2 = fpath.split(',')
            seq_strategy = 'paired'
            unzip_file(f1, f"{tag}_R1")
            unzip_file(f2, f"{tag}_R2")
            if adaptor:
                run_cmd(
                    f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} "
                    f"-A {adaptor} -o {tag}_R1_trimmed.fastq "
                    f"-p {tag}_R2_trimmed.fastq "
                    f"{tag}_R1.fastq {tag}_R2.fastq")
                os.rename(f"{tag}_R1_trimmed.fastq", f"{tag}_R1.fastq")
                os.rename(f"{tag}_R2_trimmed.fastq", f"{tag}_R2.fastq")
            if mask:
                run_cmd(
                    f"bowtie -v 0 -a --un tmp.fastq -p {thread} -t mask "
                    f"-1 {tag}_R1.fastq -2 {tag}_R2.fastq {tag}.mask.out")
                if os.path.exists("tmp_1.fastq"):
                    os.rename("tmp_1.fastq", f"{tag}_R1.fastq")
                if os.path.exists("tmp_2.fastq"):
                    os.rename("tmp_2.fastq", f"{tag}_R2.fastq")

            run_cmd(
                f"STAR --runMode alignReads --genomeDir Genome "
                f"--alignIntronMax 5000 --outSAMtype BAM SortedByCoordinate "
                f"--limitBAMsortRAM 10000000000 --outSAMmultNmax 1 "
                f"--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.1 "
                f"--runThreadN {thread} --readFilesIn "
                f"{tag}_R1.fastq {tag}_R2.fastq")
            for fn in (f"{tag}_R1.fastq", f"{tag}_R2.fastq"):
                if os.path.exists(fn):
                    os.unlink(fn)

        os.rename("Aligned.sortedByCoord.out.bam", f"{tag}.bam")

        if os.path.exists("Log.final.out"):
            with open("Log.final.out") as lf:
                tee.write(lf.read())

    # Cleanup
    for fname in ("Log.out", "Log.progress.out", "Log.final.out", "SJ.out.tab"):
        if os.path.exists(fname):
            os.unlink(fname)
    if os.path.exists("Genome"):
        run_cmd("rm -rf Genome")


def _do_assembly(tags, thread, ref_gtf, prefix, genome, min_fpkm, tee):
    """Run StringTie on each BAM; output per-sample GTF."""

    tee.write("\n━━━━ Stage 2 — StringTie transcript assembly ━━━━\n")

    if not os.path.exists(ref_gtf):
        # Convert ref GFF → GTF
        ref_gff_orig = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        ref_fasta = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        run_cmd(
            f"gffread -T -C -o {ref_gtf} -g {ref_fasta} {ref_gff_orig}")

    # StringTie uses annotation as a guide; novel transcripts still emitted
    for tag in tags:
        bam = f"{tag}.bam"
        if not os.path.exists(bam):
            tee.write(f"  ✗ {bam} missing — skip StringTie for {tag}\n")
            continue

        tee.write(f"  StringTie {tag} ...\n")
        run_cmd(
            f"stringtie {bam} -G {ref_gtf} -o {tag}_assembled.gtf "
            f"-p {thread} -f {min_fpkm} -l {tag} "
            f"--merge -A {tag}_gene_abund.tab "
            f"-C {tag}_cov.gtf")

    # Clean auxiliary StringTie outputs we don't need right now
    for fn in globmod.glob("*_cov.gtf") + globmod.glob("*_gene_abund.tab"):
        os.unlink(fn)


def _do_merge(tags, ref_gtf, thread, tee):
    """Merge all per-sample StringTie GTFs + reference GTF → merged.gtf."""

    tee.write("\n━━━━ Stage 2b — Merge transcriptomes (StringTie --merge) ━━━━\n")

    sample_gtfs = [f"{t}_assembled.gtf" for t in tags
                   if os.path.exists(f"{t}_assembled.gtf")]

    if not sample_gtfs:
        tee.write("  No assembled GTFs — abort merge\n")
        sys.exit(1)

    with open("mergelist.txt", 'w') as fh:
        for g in sample_gtfs:
            fh.write(g + "\n")

    merged_out = "merged.gtf"
    run_cmd(
        f"stringtie --merge -G {ref_gtf} -o {merged_out} "
        f"-p {thread} mergelist.txt")

    tee.write(f"  Merged GTF → {merged_out} "
              f"({_count_gtf_transcripts(merged_out)} tx)\n")
    return merged_out


def _do_gffcompare(merged_gtf, ref_gtf, tee):
    """
    Run gffcompare merged.gtf vs reference.

    gffcompare produces 'merged.gtf.class_code' file + 'gffcmp.annotated.gtf'
    class codes we care about:
      '=' complete match to reference isoform
      'c' contained in reference
      'j' multi-exon w/ at least one junction shared with ref (potential isoform)
      'e' single exon transcript overlapping ref exon
      'i' fully contained in reference intron (intronic)
      'o' generic exonic overlap with ref
      'p' "polymerase run-on" — same strand extending known tx end
      'r' repeat
      'u' intergenic (no overlap with any ref feature)
      'x' anti-sense overlap with ref
      's' anti-sense total intronic
      '.' unknown / no overlap
    """

    tee.write("\n━━━━ Stage 2c — gffcompare classify vs reference ━━━━\n")

    out_prefix = "gffcmp"
    run_cmd(
        f"gffcompare -r {ref_gtf} -o {out_prefix} {merged_gtf}")

    class_gtf = f"{out_prefix}.annotated.gtf"
    class_tbl = f"{out_prefix}.tmap"          # transcript map w/ class codes
    class_code = f"{out_prefix}.class_code.tsv"  # per-transcript class

    if not os.path.exists(class_gtf):
        tee.write(f"  ✗ gffcompare output missing ({class_gtf})\n")
        sys.exit(1)

    tee.write(f"  gffcompare annotated GTF → {class_gtf}\n")
    tee.write(f"  Transcripts: {_count_gtf_transcripts(class_gtf)}\n")

    _summarize_class_codes(class_code, tee)
    return class_gtf


def _summarize_class_codes(class_code_file, tee):
    """Print class-code summary."""
    if not os.path.exists(class_code_file):
        return
    codes = {}
    with open(class_code_file) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                c = parts[0]
                codes[c] = codes.get(c, 0) + 1

    tee.write("  Class-code summary:\n")
    for c, n in sorted(codes.items(), key=lambda x: -x[1]):
        label = _class_code_label(c)
        tee.write(f"    {c:>3s}  {n:>6d}  {label}\n")


def _class_code_label(code):
    m = {
        '=': 'match ref isoform (known)',
        'c': 'contained in ref (known)',
        'j': 'multi-exon w/ shared junction',
        'e': 'single-exon overlap ref',
        'i': 'intronic (in ref intron)',
        'o': 'exonic overlap ref',
        'p': 'polymerase run-on',
        'r': 'repeat',
        'u': 'intergenic (novel)',
        'x': 'anti-sense overlap ref',
        's': 'anti-sense intronic',
        '.': 'unknown / no overlap',
    }
    return m.get(code, 'other')


def _do_classification(class_gtf, ref_gtf, ref_fasta, classifier,
                       feelnc_dir, plek2_dir, min_tx_len, thread, tee):
    """
    Split novel transcripts into coding vs lncRNA.

    Strategy:
      1. Extract NOVEL tx from gffcompare output  (class in u, x, i, s, ., p, o, e, j)
      2. gffread → novel fasta + novel.gtf
      3. If FEELnc / PLEK2 available → run classifier
         Else fall back to: ORF length + overlap ref protein coding genes as proxy
      4. Merge known coding, known ncRNA, novel coding, novel lncRNA
         → transcripts_final.gtf  (with extra 'class' attribute:
            protein_coding | known_ncRNA | novel_coding | novel_lncRNA)
    """

    tee.write("\n━━━━ Stage 2d — Coding/lncRNA classification ━━━━\n")

    # ── Parse class codes from gffcompare output ──────────────────────
    class_code_file = "gffcmp.class_code.tsv"
    novel_class_set = {'u', 'x', 'i', 's', '.', 'p', 'o', 'e', 'j'}

    tx_info = {}  # tid -> {'code': ..., 'ref_gene': ...}
    if os.path.exists(class_code_file):
        with open(class_code_file) as fh:
            next(fh, None)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    tid, code = parts[1], parts[0]
                    ref_gene = parts[2] if len(parts) > 2 else ''
                    tx_info[tid] = {'code': code, 'ref_gene': ref_gene}

    # ── Extract novel transcripts into novel.gtf ──────────────────────
    novel_tids = set()
    with open(class_gtf) as fh, open("novel.gtf", 'w') as out:
        for line in fh:
            if line.startswith('#'):
                out.write(line)
                continue
            attrs = _parse_gtf_attrs(line)
            tid = attrs.get('transcript_id')
            if tid and tid in tx_info:
                code = tx_info[tid]['code']
                if code in novel_class_set:
                    # Length filter
                    length = _tx_length_from_gtf_line(line)
                    if length >= min_tx_len:
                        novel_tids.add(tid)
                        out.write(line)

    tee.write(f"  Novel transcripts retained (≥ {min_tx_len} bp): "
              f"{len(novel_tids)}\n")

    # ── Extract novel transcript sequences ────────────────────────────
    novel_fasta = "novel_transcripts.fa"
    run_cmd(
        f"gffread novel.gtf -g {ref_fasta} -w {novel_fasta} "
        f"-M --extract-local-tags")

    # ── Classify novel fasta → coding potential ───────────────────────
    # Dict: tid -> 'novel_coding' | 'novel_lncRNA'
    novel_class = _run_classifier(
        novel_fasta, classifier, feelnc_dir, plek2_dir, thread, tee
    )

    # ── Annotate merged.gtf with class attribute ──────────────────────
    # First, tag known transcripts from reference
    ref_tx_ids = _read_ref_tx_ids(ref_gtf)
    known_coding, known_ncrna = ref_tx_ids

    # Build class lookup
    tx_class = {}
    for tid in ref_tx_ids[0]:
        tx_class[tid] = 'protein_coding'
    for tid in ref_tx_ids[1]:
        tx_class[tid] = 'known_ncRNA'

    for tid, cls in novel_class.items():
        tx_class[tid] = cls

    # Write transcripts_final.gtf with class attribute
    _write_final_gtf(class_gtf, tx_class, tee)

    return "transcripts_final.gtf"


def _run_classifier(novel_fasta, classifier, feelnc_dir, plek2_dir,
                    thread, tee):
    """Run FEELnc / PLEK2; return dict tid -> 'novel_coding|novel_lncRNA'."""

    tee.write(f"  Classifying novel transcripts with: {classifier}\n")

    out = {}

    if classifier == 'feelnc':
        # ── FEELnc (bioconda-installed or user-supplied dir) ──────────
        # bioconda 安装后 FEELnc_cmd_classifier.pl 自带默认资源路径，
        # 只有当用户手动装旧版 FEELnc 时才需要传 feelnc_dir。
        probe = subprocess.run(
            "which FEELnc_cmd_classifier.pl 2>&1 || echo NOT_FOUND",
            shell=True, capture_output=True, text=True)

        if probe.returncode == 0 and probe.stdout.strip() != 'NOT_FOUND':
            # bioconda 版 — 不带 -m/-k 参数（自动用内置资源）
            cmd = f"FEELnc_cmd_classifier.pl -i {novel_fasta} -o feelnc_out"
            if feelnc_dir:
                cmd += (f" -m {feelnc_dir}/lnc_RNA.dat "
                        f"-k {feelnc_dir}/codon_usage.table")
            cmd += f" -p {thread}"

            run_cmd(cmd, check=False)

            feelnc_res = os.path.join("feelnc_out",
                                      "novel_transcripts.fa.classifier")
            if os.path.exists(feelnc_res):
                out = _parse_feelnc(feelnc_res, tee)
            else:
                tee.write("    FEELnc output missing → fallback heuristic\n")
        else:
            tee.write("    FEELnc not found in PATH → fallback heuristic\n")

    elif classifier == 'plek2':
        # ── PLEK2 (deep-learning, plant model) ────────────────────────
        # _ensure_plek2_installed() auto-clones + downloads h5 model
        plek2_script = _ensure_plek2_installed(plek2_dir, tee)

        if plek2_script:
            # PLEK2.py writes a 'results' file of 0/1 labels per line,
            # order = input fasta order.  We need to read fasta order too.
            tee.write(f"    Running PLEK2.py -m pl (plant model)...\n")
            with open(novel_fasta) as fh:
                fasta_order = [
                    line[1:].split()[0]
                    for line in fh if line.startswith('>')
                ]

            plek2_workdir = os.path.dirname(plek2_script) or '.'
            run_cmd(
                f"python3 {plek2_script} -i {novel_fasta} -m pl",
                check=False, cwd=plek2_workdir)

            results_file = os.path.join(plek2_workdir, "results")
            if os.path.exists(results_file):
                with open(results_file) as fh:
                    labels = [l.strip() for l in fh if l.strip()]
                # Pair up: (tid, label) where 0=lncRNA, 1=coding
                for tid, lab in zip(fasta_order, labels):
                    out[tid] = ('novel_lncRNA' if lab == '0'
                                else 'novel_coding')
            else:
                tee.write("    PLEK2 results file missing → fallback\n")
        else:
            tee.write("    PLEK2 install failed → fallback heuristic\n")

    # ── Fallback heuristic: ORF coverage + length ─────────────────────
    if not out:
        tee.write("    Using ORF-based fallback classification\n")
        out = _orf_heuristic(novel_fasta, tee)

    # Summarize
    n_ld = sum(1 for v in out.values() if v == 'novel_lncRNA')
    n_cd = sum(1 for v in out.values() if v == 'novel_coding')
    tee.write(f"  Novel → lncRNA: {n_ld} | coding: {n_cd}\n")

    return out


def _parse_feelnc(feelnc_res, tee):
    """Parse FEELnc_cmd_classifier.pl output file."""
    result = {}
    with open(feelnc_res) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                tid = parts[0].split()[0].lstrip('>')
                cat = parts[1].strip().lower()
                result[tid] = (
                    'novel_lncRNA' if 'non' in cat
                    else 'novel_coding'
                )
    return result


def _ensure_plek2_installed(plek2_dir, tee):
    """
    Ensure PLEK2 is fully installed at plek2_dir.

    Steps:
      1. Git clone https://github.com/emanlee/plek2 if not present
      2. pip install keras==2.4.3 tensorflow==2.4.1 bio pandas regex numpy==1.19.2
      3. Download PLEK2_model_v3.tar.gz from SourceForge
      4. bunzip2 both .h5 model files (Arabidopsis + vertebrate)

    Returns path to PLEK2.py if successful, None otherwise.
    """
    import shutil

    plek2_py = os.path.join(plek2_dir, "PLEK2.py")
    model_pl = os.path.join(plek2_dir, "Coding_Net_kmer6_orf_Arabidopsis.h5")

    # Already fully installed?
    if os.path.exists(plek2_py) and os.path.exists(model_pl):
        tee.write(f"    PLEK2 already installed at {plek2_dir}\n")
        return plek2_py

    tee.write(f"    Setting up PLEK2 at {plek2_dir} ...\n")

    # Step 1: clone
    if not os.path.exists(plek2_dir):
        os.makedirs(plek2_dir, exist_ok=True)
        ret = subprocess.run(
            ["git", "clone", "--depth", "1",
             "https://github.com/emanlee/plek2.git", plek2_dir],
            capture_output=True, text=True)
        if ret.returncode != 0:
            tee.write(f"    ✗ git clone failed: {ret.stderr[:200]}\n")
            return None

    # Step 2: pip deps  (do best-effort; keras 2.4.3 needs tf 2.4.1 + numpy 1.19.2)
    deps = [
        'numpy==1.19.2', 'pandas', 'regex', 'bio',
        'keras==2.4.3', 'tensorflow==2.4.1'
    ]
    for dep in deps:
        try:
            subprocess.run(
                [sys.executable, '-m', 'pip', 'install', '-q', dep],
                timeout=180, check=False)
        except (subprocess.TimeoutExpired, Exception):
            pass

    # Step 3: download PLEK2_model_v3.tar.gz from SourceForge
    # URL: https://sourceforge.net/projects/plek2/files/PLEK2_model_v3.tar.gz/download
    model_tar = os.path.join(plek2_dir, "PLEK2_model_v3.tar.gz")
    if not os.path.exists(os.path.join(plek2_dir, "Coding_Net_kmer6_orf.h5")):
        if not os.path.exists(model_tar):
            tee.write("    Downloading PLEK2 model from SourceForge...\n")
            ret = subprocess.run(
                ["curl", "-sL", "-o", model_tar,
                 "https://sourceforge.net/projects/plek2/files/"
                 "PLEK2_model_v3.tar.gz/download"],
                capture_output=True, text=True, timeout=300)
            if ret.returncode != 0 or not os.path.exists(model_tar):
                tee.write("    ✗ Model download failed (network?) — "
                          "try manually downloading to "
                          f"{plek2_dir}/PLEK2_model_v3.tar.gz\n")
                return None

        # Extract
        subprocess.run(
            f"tar xzf {model_tar} -C {plek2_dir}",
            shell=True, capture_output=True)

    # Step 4: decompress .bz2 models
    for bz in globmod.glob(os.path.join(plek2_dir, "*.bz2")):
        subprocess.run(f"bunzip2 -f {bz}", shell=True, capture_output=True)

    # Verify model files
    h5_count = len(globmod.glob(os.path.join(plek2_dir, "*.h5")))
    if h5_count == 0:
        tee.write("    ✗ No .h5 model files after extraction\n")
        return None

    tee.write(f"    PLEK2 ready: {plek2_py} + {h5_count} model(s)\n")
    return plek2_py


def _orf_heuristic(fasta_path, tee):
    """
    Simple ORF-based classifier:
      - ORF length ≥ 100 aa (300 nt) AND ORF covers ≥ 50% of transcript → coding
      - Otherwise → lncRNA
    Good enough as a fallback for plant work when no proper classifier installed.
    """
    out = {}
    seqs = _read_fasta(fasta_path)

    for tid, seq in seqs.items():
        orf_len, cov = _find_best_orf(seq)
        if orf_len >= 300 and cov >= 0.50:
            out[tid] = 'novel_coding'
        else:
            out[tid] = 'novel_lncRNA'

    return out


def _find_best_orf(seq):
    """Return (best_orf_nt_length, coverage_of_full_transcript)."""
    seq = seq.upper()
    stop_codons = {'TAA', 'TAG', 'TGA'}
    best = 0

    for frame in range(3):
        # Plus strand
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon == 'ATG':
                # Walk until stop
                j = i + 3
                while j + 3 <= len(seq):
                    if seq[j:j+3] in stop_codons:
                        break
                    j += 3
                orf_len = j - i
                if orf_len > best:
                    best = orf_len

        # Reverse complement frame
        rc = _revcomp(seq)
        for i in range(frame, len(rc) - 2, 3):
            codon = rc[i:i+3]
            if codon == 'ATG':
                j = i + 3
                while j + 3 <= len(rc):
                    if rc[j:j+3] in stop_codons:
                        break
                    j += 3
                orf_len = j - i
                if orf_len > best:
                    best = orf_len

    cov = best / len(seq) if len(seq) > 0 else 0
    return best, cov


def _revcomp(s):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
            'N': 'N', 'U': 'A'}
    return ''.join(comp.get(c, 'N') for c in reversed(s))


def _read_fasta(path):
    fasta = {}
    name, seqs = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith('>'):
                if name:
                    fasta[name] = ''.join(seqs)
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
        if name:
            fasta[name] = ''.join(seqs)
    return fasta


def _read_ref_tx_ids(ref_gtf):
    """Return (set of protein_coding tx ids, set of known ncRNA tx ids)."""
    coding, ncrna = set(), set()
    if not os.path.exists(ref_gtf):
        return coding, ncrna

    with open(ref_gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            if parts[2] not in ('transcript',):
                continue
            attrs = _parse_gtf_attrs(line)
            tid = attrs.get('transcript_id')
            biotype = (
                attrs.get('transcript_biotype') or
                attrs.get('gene_biotype') or
                attrs.get('gene_type') or
                attrs.get('transcript_type') or
                ''
            ).lower()
            if not tid:
                continue
            if 'coding' in biotype or biotype == '':
                # Default unknown → coding (plants often lack biotype attrs)
                coding.add(tid)
            elif any(k in biotype for k in ('lnc', 'ncrna', 'long_non')):
                ncrna.add(tid)
            else:
                coding.add(tid)

    return coding, ncrna


def _write_final_gtf(class_gtf, tx_class, tee):
    """Add 'class=protein_coding|known_ncRNA|novel_coding|novel_lncRNA'
    to every transcript line and write transcripts_final.gtf.

    Also emit class-level summary tables.
    """

    summary = {'protein_coding': 0, 'known_ncRNA': 0,
               'novel_coding': 0, 'novel_lncRNA': 0, 'unassigned': 0}

    with open(class_gtf) as fh, open("transcripts_final.gtf", 'w') as out, \
         open("transcript_classification.tsv", 'w') as cls:
        cls.write("transcript_id\tclass\tgene_id\n")
        for line in fh:
            if line.startswith('#'):
                out.write(line)
                continue

            attrs = _parse_gtf_attrs(line)
            tid = attrs.get('transcript_id', '')
            gid = attrs.get('gene_id', '')
            cls_name = tx_class.get(tid, 'unassigned')

            summary[cls_name] = summary.get(cls_name, 0) + 1

            # Append class attribute
            attrs['class'] = cls_name
            out.write(_write_gtf_line(line, attrs))

            if tid:
                cls.write(f"{tid}\t{cls_name}\t{gid}\n")

    tee.write("  transcripts_final.gtf:\n")
    for k, v in summary.items():
        tee.write(f"    {k:>14s}: {v}\n")


def _write_gtf_line(line, new_attrs):
    """Rewrite GTF line with updated attributes."""
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 9:
        return line
    attrs_str = '; '.join(f'{k} "{v}"' for k, v in new_attrs.items() if v) + ';'
    parts[8] = attrs_str
    return '\t'.join(parts) + '\n'


def _parse_gtf_attrs(line):
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 9:
        return {}
    attrs = {}
    for kv in parts[8].rstrip(';').split(';'):
        kv = kv.strip()
        if ' ' in kv:
            k, _, v = kv.partition(' ')
            attrs[k] = v.strip('"')
    return attrs


def _tx_length_from_gtf_line(line):
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 9:
        return 0
    # For transcript lines, compute from exon ranges.
    if parts[2] == 'exon':
        return int(parts[4]) - int(parts[3]) + 1
    return 0


def _count_gtf_transcripts(gtf_path):
    n = 0
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3 and parts[2] == 'transcript':
                n += 1
    return n


def _quantify_one(tag, bam, thread, gtf_file, ref_fasta, seq_strategy, tee):
    """Run featureCounts on BAM using the full merged GTF."""

    tee.write(f"\n  Quantifying {tag} ...\n")
    run_cmd(f"samtools index {bam}")

    if seq_strategy == 'paired':
        run_cmd(
            f"featureCounts -T {thread} -p --countReadPairs -BCO "
            f"-G {ref_fasta} -s 0 -t transcript "
            f"-a {gtf_file} -o {tag}_counts.txt {bam}")
    else:
        run_cmd(
            f"featureCounts -T {thread} -O -G {ref_fasta} -s 0 "
            f"-t transcript -a {gtf_file} -o {tag}_counts.txt {bam}")

    # Clean featureCounts output → two columns: transcript_id count
    with open(f"{tag}_counts.txt") as fh, open(f"{tag}.txt", 'w') as out:
        fh.readline(); fh.readline()  # skip header
        out.write("Transcript\tCount\n")
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 7:
                out.write(f"{cols[0]}\t{cols[6]}\n")
