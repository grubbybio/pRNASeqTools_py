"""
Single-cell RNA-seq analysis mode.
STARsolo alignment/quantification → Cell count matrix → Seurat analysis → Pseudotime.

Input formats supported:
  - FASTQ (raw, SRA, single-end, paired-end)
  - BAM (auto-detects cell-tagged vs. bulk, quantifies accordingly)
  - Count matrix (.rds, .tsv, .csv, .txt, .mtx/.mtx.gz)
"""

import os
import sys
import re
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.input_parser import (parse_input, _parse_to_dict,
                                        _resolve_path)
from prnaseqtools.functions import download_sra, unzip_file, _tee, run_cmd


def run(opts):
    """Main entry point for single-cell RNA-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    run_mode = opts.get('run_mode', 'whole')
    control = opts.get('control', '')
    treatment = opts.get('treatment')
    seq_strategy = opts.get('seq_strategy')
    genome_size = opts.get('genome_size', 10)
    min_cells = opts.get('min_cells', 3)
    min_features = opts.get('min_features', 200)
    max_features = opts.get('max_features', 5000)
    pct_mt = opts.get('pct_mt', 20)
    n_pcs = opts.get('n_pcs', 30)
    n_clusters = opts.get('n_clusters', 0)
    resolution = opts.get('resolution', 0.5)
    marker_min_pct = opts.get('marker_min_pct', 0.25)
    marker_logfc = opts.get('marker_logfc', 0.25)
    pseudotime = opts.get('pseudotime', False)
    doublet_rate = opts.get('doublet_rate', 0)
    integration = opts.get('integration', 'seurat')  # seurat | harmony | none

    mapping = run_mode in ('whole', 'mapping-only')
    do_count = run_mode != 'count-table'
    do_analysis = run_mode != 'mapping-only'

    # Parse inputs — track group membership explicitly
    # control = group 0, each -p = group 1, 2, 3, ...
    control_dict = _parse_to_dict(control)
    tags, files, pars = parse_input(control_dict)
    groups = [0] * len(tags)  # control samples → group 0

    group_id = 1
    if treatment:
        for t in (treatment if isinstance(treatment, list) else [treatment]):
            treatment_dict = _parse_to_dict(t)
            t_tags, t_files, t_pars = parse_input(treatment_dict)
            tags.extend(t_tags)
            files.extend(t_files)
            pars.extend(t_pars)
            groups.extend([group_id] * len(t_tags))
            group_id += 1

    par_str = ' '.join(pars)

    # Resolve all file paths
    resolved_files = []
    for fp in files:
        if fp.startswith('/'):
            resolved_files.append(fp)
        else:
            resolved_files.append(_resolve_path(fp))

    # ── Detect input types for each sample ────────────────────────────────
    input_types = {}  # tag -> 'fastq' | 'bam' | 'count'

    for i in range(len(tags)):
        tag = tags[i]
        fpath = resolved_files[i]
        input_types[tag] = _detect_input_type(tag, fpath, tee)

    # Report group assignments
    tee.write("\nSample group assignments:\n")
    for i, tag in enumerate(tags):
        if groups[i] == 0:
            grp_label = f'control (group 0)'
        else:
            grp_label = f'treatment_{groups[i]} (group {groups[i]})'
        tee.write(f"  {tag} → {grp_label}\n")

    # If all are count tables, skip mapping entirely
    all_counts = all(t == 'count' for t in input_types.values())
    any_mapping = any(t in ('fastq', 'bam') for t in input_types.values())

    if mapping and any_mapping:
        # ── Build STAR index (shared for FASTQ and cell-tagged BAM) ──
        tee.write("\nBuilding STARsolo genome index ...\n")

        if os.path.exists("Genome"):
            run_cmd("rm -rf Genome")
        os.makedirs("Genome", exist_ok=True)

        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")

        run_cmd(
            f"STAR --runThreadN {thread} --genomeDir Genome --runMode genomeGenerate "
            f"--genomeSAindexNbases {genome_size} --genomeFastaFiles {fasta_path} "
            f"--sjdbGTFfile {gff_path} --sjdbGTFtagExonParentTranscript Parent "
            f"--sjdbGTFtagExonParentGene ID --limitGenomeGenerateRAM 64000000000")

        gtf_file = f"{genome}_genes.gtf"
        run_cmd(f"gffread -T -o {gtf_file} -g {fasta_path} {gff_path}")
        os.makedirs("Solo_out", exist_ok=True)

        for i in range(len(tags)):
            tag = tags[i]
            fpath = resolved_files[i]
            itype = input_types[tag]

            tee.write(f"\n{'─' * 40}\nProcessing {tag} (input: {itype})...\n")

            if itype == 'count':
                # Skip mapping for count inputs
                tee.write("  Skipping mapping (count matrix input).\n")
                continue

            if itype == 'bam':
                bam_kind = _detect_bam_kind(fpath, tee)

                if bam_kind == 'bulk':
                    sys.exit(
                        f"Error: '{fpath}' is a bulk BAM (no cell barcodes detected).\n"
                        f"  The 'sc' module requires single-cell data with cell barcodes.\n"
                        f"  For bulk RNA-seq BAM analysis, please use the 'mrna' module instead:\n"
                        f"    pRNASeqTools mrna --mode_mrna bam -c {tag}={fpath} ...")

                # cell-tagged BAM: symlink + quantify
                bam_src = fpath
                if not os.path.isabs(bam_src):
                    bam_src = os.path.abspath(bam_src)
                os.symlink(bam_src, f"{tag}.bam")
                run_cmd(f"samtools index {tag}.bam")

                tee.write("  Cell-tagged BAM detected → quantifying with STARsolo.\n")
                _quantify_celltagged_bam(tag, thread, gtf_file, prefix, genome)
                continue

            # ── FASTQ input ──
            sra_results = download_sra(fpath, thread)
            if len(sra_results) == 1:
                seq_strategy = 'single'
                unzip_file(sra_results[0], tag)
                input_r1 = f"{tag}.fastq"
                input_r2 = None
            else:
                seq_strategy = 'paired'
                unzip_file(sra_results[0], f"{tag}_R1")
                unzip_file(sra_results[1], f"{tag}_R2")
                input_r1 = f"{tag}_R1.fastq"
                input_r2 = f"{tag}_R2.fastq"

            # Adapter trimming
            if adaptor:
                tee.write(f"  Trimming adapters ({tag})...\n")
                if input_r2:
                    run_cmd(
                        f"cutadapt -j {thread} -m 20 --trim-n "
                        f"-a {adaptor} -A {adaptor} "
                        f"-o {tag}_R1_trimmed.fastq -p {tag}_R2_trimmed.fastq "
                        f"{input_r1} {input_r2}")
                    os.rename(f"{tag}_R1_trimmed.fastq", input_r1)
                    os.rename(f"{tag}_R2_trimmed.fastq", input_r2)
                else:
                    run_cmd(
                        f"cutadapt -j {thread} -m 20 --trim-n -a {adaptor} "
                        f"-o {tag}_trimmed.fastq {input_r1}")
                    os.rename(f"{tag}_trimmed.fastq", input_r1)

            # STARsolo alignment and quantification
            tee.write(f"  Running STARsolo alignment ({tag})...\n")
            solo_outdir = f"Solo_out/{tag}"

            star_cmd = (
                f"STAR --runMode alignReads --genomeDir Genome "
                f"--soloType SmartSeq "
                f"--outSAMtype BAM SortedByCoordinate "
                f"--limitBAMsortRAM 10000000000 "
                f"--outSAMmultNmax 1 "
                f"--outFilterMultimapNmax 50 "
                f"--outFilterMismatchNoverLmax 0.1 "
                f"--runThreadN {thread} "
                f"--soloOutDir {solo_outdir} "
                f"--soloCBwhitelist None "
                f"--soloCellFilter EmptyDroplets "
                f"--quantMode GeneCounts "
                f"--soloFeatures Gene GeneFull"
            )

            if input_r2:
                star_cmd += f" --readFilesIn {input_r1} {input_r2}"
            else:
                star_cmd += f" --readFilesIn {input_r1}"

            run_cmd(star_cmd)

            # Show mapping stats
            if os.path.exists("Log.final.out"):
                with open("Log.final.out") as lf:
                    tee.write(lf.read())

            # Clean up raw FASTQ files
            for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                         f"{tag}.fastq"):
                if os.path.exists(fname):
                    os.unlink(fname)

        # Clean up STAR files
        for fname in ("Log.out", "Log.progress.out", "Log.final.out", "SJ.out.tab"):
            if os.path.exists(fname):
                os.unlink(fname)
        if os.path.exists("Genome"):
            run_cmd("rm -rf Genome")

        # ── Convert STARsolo output and BAM counts to .rds ───────────────
        tee.write("\nPreparing count matrices for Seurat...\n")
        gtf_file = f"{genome}_genes.gtf"
        _prepare_all_counts(tags, input_types, solo_dir="Solo_out",
                            gtf_file=gtf_file, prefix=prefix, genome=genome,
                            thread=thread)

        if do_analysis:
            _run_seurat_analysis(tags, groups, prefix, genome, thread,
                                 min_cells=min_cells, min_features=min_features,
                                 max_features=max_features, pct_mt=pct_mt,
                                 n_pcs=n_pcs, n_clusters=n_clusters,
                                 resolution=resolution,
                                 marker_min_pct=marker_min_pct,
                                 marker_logfc=marker_logfc,
                                 pseudotime=pseudotime,
                                 doublet_rate=doublet_rate,
                                 integration=integration)

    elif mapping and not any_mapping:
        # Mapping requested but all inputs are count tables
        tee.write("\nAll inputs are count tables — skipping mapping step.\n")
        gtf_file = f"{genome}_genes.gtf"
        if not os.path.exists(gtf_file):
            gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
            fasta_path = os.path.join(prefix, "reference",
                                      f"{genome}_chr_all.fasta")
            run_cmd(f"gffread -T -o {gtf_file} -g {fasta_path} {gff_path}")

        # Convert plain text count matrices to .rds
        for i in range(len(tags)):
            tag = tags[i]
            fpath = resolved_files[i]
            _import_count_matrix(tag, fpath, tee)

        if do_analysis:
            _run_seurat_analysis(tags, groups, prefix, genome, thread,
                                 min_cells=min_cells, min_features=min_features,
                                 max_features=max_features, pct_mt=pct_mt,
                                 n_pcs=n_pcs, n_clusters=n_clusters,
                                 resolution=resolution,
                                 marker_min_pct=marker_min_pct,
                                 marker_logfc=marker_logfc,
                                 pseudotime=pseudotime,
                                 doublet_rate=doublet_rate,
                                 integration=integration)

    else:
        # ── No mapping mode (count-table or mapping-only with existing data) ─
        tee.write("\nNo-mapping mode — using existing files.\n")
        gtf_file = f"{genome}_genes.gtf"
        if not os.path.exists(gtf_file):
            gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
            fasta_path = os.path.join(prefix, "reference",
                                      f"{genome}_chr_all.fasta")
            run_cmd(f"gffread -T -o {gtf_file} -g {fasta_path} {gff_path}")

        for i in range(len(tags)):
            tag = tags[i]
            fpath = resolved_files[i]
            itype = input_types[tag]

            tee.write(f"\n  Processing {tag} (input: {itype})...\n")

            if itype == 'count':
                # Import plain text or .rds count matrix
                _import_count_matrix(tag, fpath, tee)

            elif itype == 'bam':
                bam_kind = _detect_bam_kind(fpath, tee)

                if bam_kind == 'bulk':
                    sys.exit(
                        f"Error: '{fpath}' is a bulk BAM (no cell barcodes detected).\n"
                        f"  The 'sc' module requires single-cell data with cell barcodes.\n"
                        f"  For bulk RNA-seq BAM analysis, please use the 'mrna' module instead:\n"
                        f"    pRNASeqTools mrna --mode_mrna bam -c {tag}={fpath} ...")

                bam_src = fpath if os.path.isabs(fpath) else f"../{fpath}"
                os.symlink(bam_src, f"{tag}.bam")
                run_cmd(f"samtools index {tag}.bam")
                tee.write("  Cell-tagged BAM detected → quantifying with STARsolo.\n")
                _quantify_celltagged_bam(tag, thread, gtf_file, prefix, genome)
                # Convert the resulting count to .rds
                _finalize_count_matrices([tag], gtf_file, tee)

            elif itype == 'fastq':
                sys.exit(
                    f"FASTQ file '{fpath}' found but --no-mapping mode "
                    f"cannot process raw reads.\n"
                    f"  Please run without --no-mapping for FASTQ inputs, "
                    f"or provide a pre-aligned BAM or count matrix.")

        if do_count:
            # Convert any remaining count outputs
            _finalize_count_matrices(tags, gtf_file, tee)

        if do_analysis:
            _run_seurat_analysis(tags, groups, prefix, genome, thread,
                                 min_cells=min_cells, min_features=min_features,
                                 max_features=max_features, pct_mt=pct_mt,
                                 n_pcs=n_pcs, n_clusters=n_clusters,
                                 resolution=resolution,
                                 marker_min_pct=marker_min_pct,
                                 marker_logfc=marker_logfc,
                                 pseudotime=pseudotime,
                                 doublet_rate=doublet_rate,
                                 integration=integration)

        # Clean up symlinks
        for fname in globmod.glob("*.bam"):
            if os.path.islink(fname):
                os.unlink(fname)
        for fname in globmod.glob("*.tsv") + globmod.glob("*.csv"):
            if os.path.islink(fname):
                os.unlink(fname)


# ═══════════════════════════════════════════════════════════════════════════
# Input type detection
# ═══════════════════════════════════════════════════════════════════════════

def _detect_input_type(tag, fpath, tee):
    """Detect input type: 'fastq', 'bam', or 'count'."""
    lower = fpath.lower()

    # Count matrix formats
    count_exts = ('.rds', '.tsv', '.csv', '.txt', '.mtx', '.mtx.gz',
                  '.mtx.bgz')
    for ext in count_exts:
        if lower.endswith(ext):
            if not os.path.exists(fpath):
                tee.write(f"  Warning: Count file not found: {fpath}\n")
            return 'count'

    # BAM format
    if lower.endswith('.bam'):
        return 'bam'

    # FASTQ / SRA
    fastq_exts = ('.fastq', '.fq', '.fastq.gz', '.fq.gz',
                  '.fastq.bz2', '.fq.bz2', '.gtz', '.sra')
    for ext in fastq_exts:
        if lower.endswith(ext):
            return 'fastq'

    # SRA accession
    if re.match(r'^[SED]RR\d+$', os.path.basename(fpath)):
        return 'fastq'

    # Default: treat as FASTQ/SRA
    tee.write(f"  Unknown input format for {tag} ({fpath}), treating as FASTQ.\n")
    return 'fastq'


def _detect_bam_kind(bam_path, tee):
    """
    Detect whether a BAM file already has cell barcodes (CB tag).

    Returns:
        'cell_tagged' - BAM has CB/BZ tags (already cell-level, e.g. STARsolo output)
        'bulk'        - Regular BAM without cell barcodes
    """
    tee.write(f"  Analyzing BAM: {os.path.basename(bam_path)}\n")

    # Check for CB tag (cell barcode) in first 1000 reads
    cmd = (
        f"samtools view {bam_path} | head -1000 | "
        f"grep -o 'CB:Z:[^\\t]*' | head -1"
    )
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True,
                                text=True, timeout=60)
        if result.stdout.strip():
            tee.write("  → Cell-tagged BAM detected (CB tag found)\n")
            return 'cell_tagged'
    except subprocess.TimeoutExpired:
        tee.write("  → Timeout scanning BAM, assuming bulk BAM\n")
        return 'bulk'

    # Also check for BZ tag (STARsolo cell ID)
    cmd2 = (
        f"samtools view {bam_path} | head -1000 | "
        f"grep -o 'BZ:Z:[^\\t]*' | head -1"
    )
    try:
        result2 = subprocess.run(cmd2, shell=True, capture_output=True,
                                 text=True, timeout=60)
        if result2.stdout.strip():
            tee.write("  → Cell-tagged BAM detected (BZ tag found)\n")
            return 'cell_tagged'
    except subprocess.TimeoutExpired:
        pass

    # Check if BAM has proper read groups with cell information
    header_cmd = f"samtools view -H {bam_path}"
    try:
        hdr = subprocess.run(header_cmd, shell=True, capture_output=True,
                              text=True, timeout=30)
        # Look for cell barcode in RG field
        if 'CB:' in hdr.stdout or 'cell' in hdr.stdout.lower():
            tee.write("  → Cell-tagged BAM detected (cell info in header)\n")
            return 'cell_tagged'
    except subprocess.TimeoutExpired:
        pass

    # Check for @RG or @CB in header
    try:
        id_cmd = (
            f"samtools view -H {bam_path} | grep -c 'CB:'"
        )
        result3 = subprocess.run(id_cmd, shell=True, capture_output=True,
                                 text=True, timeout=30)
        if result3.stdout.strip() and int(result3.stdout.strip()) > 0:
            tee.write("  → Cell-tagged BAM detected (CB in @RG)\n")
            return 'cell_tagged'
    except (subprocess.TimeoutExpired, ValueError):
        pass

    # If we get here, it's a bulk BAM
    tee.write("  → Bulk BAM (no cell barcodes detected)\n")
    return 'bulk'


# ═══════════════════════════════════════════════════════════════════════════
# BAM quantification
# ═══════════════════════════════════════════════════════════════════════════

def _quantify_celltagged_bam(tag, thread, gtf_file, prefix, genome):
    """Quantify a cell-tagged BAM using STARsolo quant mode."""
    tee = _tee()
    bam_file = f"{tag}.bam"

    # STARsolo quant mode: needs genome index. Build if missing.
    if not os.path.exists("Genome"):
        tee.write("  Building STAR genome index for quantification...\n")
        os.makedirs("Genome", exist_ok=True)
        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")
        run_cmd(
            f"STAR --runThreadN {thread} --genomeDir Genome "
            f"--runMode genomeGenerate "
            f"--genomeSAindexNbases 10 --genomeFastaFiles {fasta_path} "
            f"--sjdbGTFfile {gtf_file} --limitGenomeGenerateRAM 64000000000"
        )

    tee.write(f"  Running STARsolo quant for {tag}...\n")
    solo_outdir = f"Solo_out/{tag}"
    os.makedirs(solo_outdir, exist_ok=True)

    # Use STARsolo in soloQuant mode (quantify an existing BAM)
    star_cmd = (
        f"STAR --runMode soloQuant "
        f"--soloInputSAMfile {bam_file} "
        f"--soloOutDir {solo_outdir} "
        f"--soloCBwhitelist None "
        f"--soloCellFilter EmptyDroplets "
        f"--quantMode GeneCounts "
        f"--soloFeatures Gene GeneFull "
        f"--runThreadN {thread} "
        f"--outSAMtype BAM SortedByCoordinate "
        f"--limitBAMsortRAM 10000000000"
    )

    run_cmd(star_cmd)

    # Show stats
    log_file = os.path.join(solo_outdir, "Log.final.out")
    if os.path.exists(log_file):
        with open(log_file) as lf:
            tee.write(lf.read())


# ═══════════════════════════════════════════════════════════════════════════
# Count matrix import (plain text → .rds)
# ═══════════════════════════════════════════════════════════════════════════

def _import_count_matrix(tag, fpath, tee):
    """
    Import a count matrix from various formats and convert to .rds.

    Supported formats:
      - .rds          : R native format (already ready)
      - .tsv/.txt     : Tab-separated (genes × cells, first column = gene IDs)
      - .csv          : Comma-separated
      - .mtx/.mtx.gz  : MatrixMarket sparse format (accompanies genes.tsv, barcodes.tsv)
    """
    lower = fpath.lower()
    out_rds = f"{tag}_counts.rds"

    if lower.endswith('.rds'):
        # Already .rds, just symlink/copy
        if os.path.isabs(fpath):
            src = fpath
        else:
            src = os.path.abspath(fpath)
        if os.path.exists(out_rds):
            os.unlink(out_rds)
        os.symlink(src, out_rds)
        tee.write(f"  Using existing .rds file: {src}\n")
        return

    if lower.endswith('.mtx') or lower.endswith('.mtx.gz'):
        # MatrixMarket format
        mtx_dir = os.path.dirname(fpath)
        base = os.path.splitext(os.path.basename(fpath))[0]
        # Look for companion files
        barcodes_file = os.path.join(mtx_dir, "barcodes.tsv.gz")
        features_file = os.path.join(mtx_dir, "features.tsv.gz")

        if not os.path.exists(barcodes_file):
            # Try uncompressed
            barcodes_file = os.path.join(mtx_dir, "barcodes.tsv")
            features_file = os.path.join(mtx_dir, "features.tsv")

        if not os.path.exists(barcodes_file):
            # If no barcodes file, use generic names
            barcodes_file = ""
            features_file = ""

        tee.write(f"  Converting MatrixMarket {fpath} to Seurat format...\n")

        r_cmd = (
            f'Rscript --vanilla -e \''
            f'mtx <- as.matrix(Matrix::readMM("{fpath}")); '
        )
        if features_file and os.path.exists(features_file):
            r_cmd += f'features <- read.delim("{features_file}", header=FALSE); '
            r_cmd += f'rownames(mtx) <- if(ncol(features) >= 2) features$V2 else features$V1; '
        else:
            r_cmd += f'rownames(mtx) <- paste0("Gene", seq_len(nrow(mtx))); '

        if barcodes_file and os.path.exists(barcodes_file):
            r_cmd += f'barcodes <- read.delim("{barcodes_file}", header=FALSE); '
            r_cmd += f'colnames(mtx) <- paste0("{tag}_", barcodes$V1); '
        else:
            r_cmd += f'colnames(mtx) <- paste0("{tag}_Cell", seq_len(ncol(mtx))); '

        r_cmd += f'saveRDS(mtx, "{out_rds}"); '
        r_cmd += f'cat("Converted MatrixMarket to {out_rds} with", nrow(mtx), '
        r_cmd += f'"genes and", ncol(mtx), "cells\\n")'
        r_cmd += f'\''

        # Try to load Matrix package, fall back to simple reader
        run_cmd(r_cmd)
        return

    # Tab-separated (.tsv, .txt) or comma-separated (.csv)
    sep = ',' if lower.endswith('.csv') else '\\t'
    header = 'TRUE'  # Assume header row with cell IDs
    row_names_col = '1'  # First column is gene names

    tee.write(f"  Converting plain text {fpath} to Seurat format...\n")

    r_cmd = (
        f'Rscript --vanilla -e \''
        f'data <- read.table("{fpath}", header={header}, sep="{sep}", '
        f'row.names={row_names_col}, stringsAsFactors=FALSE, check.names=FALSE); '
        f'counts <- as.matrix(data); '
        f'# Ensure numeric '
        f'storage.mode(counts) <- "integer"; '
        f'saveRDS(counts, "{out_rds}"); '
        f'cat("Converted", basename("{fpath}"), "with", nrow(counts), '
        f'"genes and", ncol(counts), "cells to {out_rds}\\n")'
        f'\''
    )
    run_cmd(r_cmd)


# ═══════════════════════════════════════════════════════════════════════════
# Count matrix preparation
# ═══════════════════════════════════════════════════════════════════════════

def _prepare_all_counts(tags, input_types, solo_dir="Solo_out", gtf_file=None,
                         prefix=None, genome=None, thread=None):
    """Prepare all count matrices (from STARsolo, BAM counts, or plain text) to .rds."""
    tee = _tee()

    for tag in tags:
        itype = input_types.get(tag, 'fastq')

        if itype == 'count':
            # Already handled by _import_count_matrix in no-mapping path
            continue

        # Check for STARsolo output (from alignment or soloQuant)
        matrix_file = os.path.join(solo_dir, tag, "Gene", "raw", "matrix.mtx.gz")

        if os.path.exists(matrix_file):
            # STARsolo output → convert
            barcodes_file = os.path.join(solo_dir, tag, "Gene", "raw",
                                          "barcodes.tsv.gz")
            features_file = os.path.join(solo_dir, tag, "Gene", "raw",
                                          "features.tsv.gz")

            tee.write(f"  Converting STARsolo output for {tag}...\n")

            r_cmd = (
                f'Rscript --vanilla -e \''
                f'mtx <- as.matrix(Matrix::readMM("{matrix_file}")); '
            )
            if os.path.exists(features_file):
                r_cmd += f'features <- read.delim("{features_file}", header=FALSE); '
                r_cmd += f'rownames(mtx) <- features$V2; '
            else:
                r_cmd += f'rownames(mtx) <- paste0("Gene", seq_len(nrow(mtx))); '

            if os.path.exists(barcodes_file):
                r_cmd += f'barcodes <- read.delim("{barcodes_file}", header=FALSE); '
                r_cmd += f'colnames(mtx) <- paste0("{tag}_", barcodes$V1); '
            else:
                r_cmd += f'colnames(mtx) <- paste0("{tag}_Cell", seq_len(ncol(mtx))); '

            r_cmd += f'saveRDS(mtx, "{tag}_counts.rds"); '
            r_cmd += f'cat("Saved {tag}_counts.rds with", nrow(mtx), '
            r_cmd += f'"genes and", ncol(mtx), "cells\\n")'
            r_cmd += f'\''
            run_cmd(r_cmd)

    # Generate combined matrix if multiple samples
    rds_files = [f"{t}_counts.rds" for t in tags if os.path.exists(f"{t}_counts.rds")]
    if len(rds_files) > 1:
        tee.write("\nGenerating combined count matrix...\n")
        tag_list = ",".join(t for t in tags if os.path.exists(f"{t}_counts.rds"))
        r_cmd = (
            f'Rscript --vanilla -e \''
            f'tags <- strsplit("{tag_list}", ",")[[1]]; '
            f'combined <- NULL; '
            f'for (t in tags) {{'
            f'  mtx <- readRDS(paste0(t, "_counts.rds")); '
            f'  if (is.null(combined)) {{ combined <- mtx; }} '
            f'  else {{ combined <- cbind(combined, mtx); }}'
            f'}}; '
            f'saveRDS(combined, "combined_counts.rds"); '
            f'cat("Saved combined_counts.rds with", ncol(combined), "cells\\n")'
            f'\''
        )
        run_cmd(r_cmd)


def _finalize_count_matrices(tags, gtf_file, tee):
    """Safety check: warn about any tags still missing .rds files."""
    for tag in tags:
        if not os.path.exists(f"{tag}_counts.rds"):
            tee.write(f"  Warning: No count matrix (.rds) for {tag}.\n")


# ═══════════════════════════════════════════════════════════════════════════
# Seurat analysis dispatch
# ═══════════════════════════════════════════════════════════════════════════

def _run_seurat_analysis(tags, groups, prefix, genome, thread,
                         min_cells=3, min_features=200, max_features=5000,
                         pct_mt=20, n_pcs=30, n_clusters=0,
                         resolution=0.5, marker_min_pct=0.25,
                         marker_logfc=0.25, pseudotime=False,
                         doublet_rate=0, integration='seurat'):
    """Run Seurat analysis pipeline via R script."""
    tee = _tee()

    # Filter to tags that have .rds files
    valid_tags = [t for t in tags if os.path.exists(f"{t}_counts.rds")]
    if not valid_tags:
        tee.write("\nError: No count matrices (.rds) found for Seurat analysis.\n")
        tee.write("  Ensure count matrices were generated from mapping or imported.\n")
        sys.exit(1)

    tee.write(f"\nRunning Seurat analysis on {len(valid_tags)} sample(s)...\n")
    tee.write(f"  Integration method: {integration}\n")

    # Build sample info (group assignment)
    # control=0, each -p = group 1, 2, ...
    sample_info = []
    for i, tag in enumerate(tags):
        grp = groups[i] if i < len(groups) else i
        sample_info.append(f"{tag}:{grp}")
    sample_str = ",".join(sample_info)

    # Build tag list for R (only valid tags)
    tag_str = ",".join(valid_tags)

    args_str = (
        f"{prefix} {genome} {thread} "
        f"{min_cells} {min_features} {max_features} {pct_mt} "
        f"{n_pcs} {n_clusters} {resolution} "
        f"{marker_min_pct} {marker_logfc} "
        f"{int(pseudotime)} {doublet_rate} {integration} "
        f"{tag_str} {sample_str}"
    )

    script_path = os.path.join(prefix, "scripts", "sc_analysis.R")
    if not os.path.exists(script_path):
        tee.write(
            f"Error: Seurat analysis script not found: {script_path}\n"
            f"Ensure scripts/sc_analysis.R exists.\n")
        sys.exit(1)

    run_cmd(f"Rscript --vanilla {script_path} {args_str}")

    tee.write("\nSeurat analysis completed!\n")
    tee.write("Output files:\n")
    tee.write("  - seurat_clusters.rds: Seurat object with clustering\n")
    tee.write("  - cluster_markers.csv: Marker genes for each cluster\n")
    tee.write("  - umap_coordinates.csv: UMAP coordinates\n")
    if pseudotime:
        tee.write("  - pseudotime_results.rds: Monocle3 pseudotime results\n")
