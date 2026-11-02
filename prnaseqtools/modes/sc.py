"""
Single-cell RNA-seq analysis mode.
HISAT2 alignment + featureCounts per-cell quantification → Cell count matrix → Seurat analysis → Pseudotime.

Input formats supported:
  - FASTQ (raw, SRA, single-end, paired-end)
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
from prnaseqtools.functions import (download_sra, unzip_file, _tee, run_cmd,
                                     try_use_shared_hisat_index, save_hisat_index_to_reference,
                                     cleanup_alignment_outputs,
                                     hisat2_build_index, hisat2_align_single,
                                     hisat2_align_paired)


def run(opts):
    """Main entry point for single-cell RNA-seq analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    adaptor = opts.get('adaptor')
    skip_fastp = opts.get('skip_fastp', False)
    # prefix must point to the project root so that prefix/reference/{genome}_HisatIdx
    # resolves to ~/software/pRNASeqTools_py/reference/ (not prnaseqtools/reference/)
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
        # ── Build/use HISAT2 genome index (shared for FASTQ and cell-tagged BAM) ──
        tee.write("\nBuilding/using HISAT2 genome index ...\n")

        # Define reference paths first — they are also needed by the gffread call
        gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")

        hisat_idx = "HisatIdx"
        if not try_use_shared_hisat_index(genome, prefix, hisat_idx, tee):
            hisat2_build_index(fasta_path, gff_path, hisat_idx, thread, tee=tee)
            save_hisat_index_to_reference(genome, prefix, hisat_idx, tee)

        gtf_file = f"{genome}_genes.gtf"
        run_cmd(f"gffread -T -o {gtf_file} -g {fasta_path} {gff_path}")

        # HISAT2 output directory for per-cell BAMs (FASTQ input only)
        hisat_bam_dir = "Hisat_out"
        os.makedirs(hisat_bam_dir, exist_ok=True)
        cell_bams = {}  # tag → bam path

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
                sys.exit(
                    f"Error: BAM input '{fpath}' is no longer supported by the 'sc' module.\n"
                    f"  STARsolo (STAR CB_UMI_Simple) has been removed from this project.\n"
                    f"  Please re-run alignment in the 'sc' module with FASTQ inputs,\n"
                    f"  or provide a pre-computed count matrix (.rds / .mtx / .tsv).")

            # ── FASTQ input ──
            if ',' in fpath:
                f1, f2 = fpath.split(',')
                if seq_strategy is None:
                    seq_strategy = 'paired'
                unzip_file(f1, f"{tag}_R1")
                unzip_file(f2, f"{tag}_R2")
                input_r1 = f"{tag}_R1.fastq"
                input_r2 = f"{tag}_R2.fastq"
            else:
                sra_results, _from_sra = download_sra(fpath, thread)
                if len(sra_results) == 1:
                    if seq_strategy is None:
                        seq_strategy = 'single'
                    unzip_file(sra_results[0], tag)
                    input_r1 = f"{tag}.fastq"
                    input_r2 = None
                else:
                    if seq_strategy is None:
                        seq_strategy = 'paired'
                    unzip_file(sra_results[0], f"{tag}_R1")
                    unzip_file(sra_results[1], f"{tag}_R2")
                    input_r1 = f"{tag}_R1.fastq"
                    input_r2 = f"{tag}_R2.fastq"

            # Adapter trimming
            if adaptor and not skip_fastp:
                tee.write(f"  Trimming adapters ({tag})...\n")
                if input_r2:
                    run_cmd(
                        f"fastp -w {thread} --length_required 20 "
                        f"-i {input_r1} -I {input_r2} "
                        f"-o {tag}_R1_trimmed.fastq -O {tag}_R2_trimmed.fastq --adapter_sequence {adaptor} --adapter_sequence_r2 {adaptor}")
                    os.rename(f"{tag}_R1_trimmed.fastq", input_r1)
                    os.rename(f"{tag}_R2_trimmed.fastq", input_r2)
                else:
                    run_cmd(
                        f"fastp -w {thread} --length_required 20 "
                        f"-i {input_r1} -o {tag}_trimmed.fastq --adapter_sequence {adaptor}")
                    os.rename(f"{tag}_trimmed.fastq", input_r1)

            # HISAT2 alignment
            tee.write(f"  Running HISAT2 alignment ({tag})...\n")

            if input_r2:
                hisat2_align_paired(tag, input_r1, input_r2,
                                    hisat_idx, thread,
                                    out_dir=hisat_bam_dir,
                                    extra_args="-k 20 --mp 3,3 --score-min L,0,-0.4")
            else:
                hisat2_align_single(tag, input_r1, hisat_idx, thread,
                                    out_dir=hisat_bam_dir,
                                    extra_args="-k 20 --mp 3,3 --score-min L,0,-0.4")

            # Show mapping stats
            log_f = os.path.join(hisat_bam_dir, f"{tag}_Log.txt")
            if os.path.exists(log_f):
                with open(log_f) as lf:
                    tee.write(lf.read())

            cell_bams[tag] = os.path.join(hisat_bam_dir, f"{tag}.bam")

            # Clean up raw FASTQ files
            for fname in (f"{tag}_R1.fastq", f"{tag}_R2.fastq",
                         f"{tag}.fastq"):
                if os.path.exists(fname):
                    os.unlink(fname)

        # Clean up HISAT2 temp outputs
        cleanup_alignment_outputs(extra_globs=(f"{genome}_chr_all*", f"{genome}*.bt2"))
        # 不要删 HisatIdx（可能有共享 symlink）
        for fname in globmod.glob("*.bam"):
            if os.path.islink(fname):
                os.unlink(fname)

        # ── featureCounts 定量（FASTQ cell 统一跑）──────────────────
        tee.write("\nQuantifying with featureCounts...\n")
        gtf_file = f"{genome}_genes.gtf"
        bam_paths = [b for b in cell_bams.values() if os.path.exists(b)]
        bam_labels = [t for t, b in cell_bams.items() if os.path.exists(b)]

        if bam_paths:
            tee.write(f"  featureCounts on {len(bam_paths)} cell BAMs...\n")
            fc_out = "featureCounts.txt"
            run_cmd(
                f"featureCounts -a {gtf_file} -o {fc_out} "
                f"-T {thread} -s 0 -p "
                f"{' '.join(bam_paths)}"
            )

            # featureCounts 输出整理成 TSV (gene × cells, header = tag 名)
            tsv_out = "hisat_counts.tsv"
            with open(fc_out) as fh:
                lines = fh.readlines()
            # 跳过 # comment 行（第一行），从第二行开始读
            data_lines = lines[1:]  # 第二行是 header，第三行起是数据
            header = "Geneid\t" + "\t".join(bam_labels) + "\n"
            out = [header]
            for dl in data_lines[1:]:  # data_lines[0] 是原始 header，跳过
                parts = dl.strip().split('\t')
                if len(parts) < 8:
                    continue
                out.append(parts[0] + "\t" + "\t".join(parts[7:]))
            with open(tsv_out, 'w') as of:
                of.writelines(out)

            tee.write(f"  Saved {tsv_out} ({len(out)-1} genes × {len(bam_labels)} cells)\n")

            # 每个 cell 单独转 .rds —— 生成只含一列的 TSV，喂给 _import_count_matrix
            for idx, tag in enumerate(bam_labels):
                onecol_tsv = f"{tag}_onecol.tsv"
                with open(onecol_tsv, 'w') as oc:
                    oc.write("Geneid\t" + tag + "\n")
                    for dl in data_lines[1:]:
                        parts = dl.strip().split('\t')
                        if len(parts) < 8 + idx:
                            continue
                        oc.write(parts[0] + "\t" + parts[7 + idx] + "\n")
                _import_count_matrix(tag, onecol_tsv, tee)
                os.unlink(onecol_tsv)

        # ── 准备 count matrices ──────────────────────────────────────
        # FASTQ + count matrix 输入 → .rds, Seurat 下游合并
        _prepare_all_counts(tags, input_types,
                            files=dict(zip(tags, resolved_files)),
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
                sys.exit(
                    f"Error: BAM input '{fpath}' is no longer supported by the 'sc' module.\n"
                    f"  STARsolo (STAR CB_UMI_Simple) has been removed from this project.\n"
                    f"  Please provide a pre-computed count matrix (.rds / .mtx / .tsv)\n"
                    f"  when using --no-mapping mode.")

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

def _prepare_all_counts(tags, input_types, files=None, solo_dir=None,
                         gtf_file=None, prefix=None, genome=None, thread=None):
    """Prepare all count matrices (from FASTQ alignment or pre-computed files) to .rds."""
    tee = _tee()

    for tag in tags:
        itype = input_types.get(tag, 'fastq')

        if itype == 'count':
            _import_count_matrix(tag, files[tag], tee)
            continue


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
