"""
CiPS uORF analysis mode.
Identifies translated upstream ORFs (uORFs) in 5' UTRs using Ribo-seq P-site data.

Workflow:
  1. Load expressed GTF + reference FASTA → extract 5' UTRs
  2. Find all potential uORFs (ATG start, in-frame) using ORFik
  3. Compute frame-specific Ribo-seq P-site counts per uORF
  4. Filter translated uORFs (minimum uORFs: 1aa, longer uORFs: ≥2aa)
  5. Deduplicate overlapping uORFs, annotate with gene descriptions

Requires R packages: ORFik, GenomicFeatures, GenomicRanges, Biostrings,
          Rsamtools, dplyr, future.apply, openxlsx
Inputs:
  - Expressed GTF (from RIBO Taper pipeline, ribo mode output)
  - Reference genome FASTA (indexed)
  - P-site count file (from RiboTaper: count chr start strand)
"""

import os
import sys
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.functions import _tee


def run(opts):
    """Main entry point for CiPS uORF analysis."""
    opts = validate_options(opts)
    tee = _tee()

    thread = opts.get('thread', 4)
    genome = opts.get('genome', 'ath')
    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))

    # ── Input paths ───────────────────────────────────────────────────
    gtf_path = opts.get('gtf')
    fasta_path = opts.get('fasta')
    psite_file = opts.get('psite')
    output_prefix = opts.get('output', 'cips')

    if not gtf_path or not os.path.exists(gtf_path):
        # Try default: expressed GTF from ribo pipeline
        default_gtf = f"{genome}_expressed.gtf"
        if os.path.exists(default_gtf):
            gtf_path = default_gtf
        else:
            sys.exit(
                "GTF file not specified. Use --gtf or run in a ribo output directory.\n"
                "The ribo pipeline produces <genome>_expressed.gtf"
            )

    if not fasta_path:
        fasta_path = os.path.join(prefix, "reference", f"{genome}_chr_all.fasta")

    if not psite_file or not os.path.exists(psite_file):
        sys.exit(
            "P-site file not specified. Use --psite.\n"
            "This file is produced by RiboTaper and has columns: count chr start strand"
        )

    # ── Filtering thresholds ──────────────────────────────────────────
    min_inframe_counts = opts.get('min_inframe_counts', 10)
    min_inframe_perc = opts.get('min_inframe_perc', 50)
    min_psite_perc = opts.get('min_psite_perc', 30)
    gene_desc = opts.get('gene_desc', '')

    # ── Validate inputs ───────────────────────────────────────────────
    for label, fpath in [("GTF", gtf_path), ("FASTA", fasta_path),
                          ("P-site file", psite_file)]:
        if not os.path.exists(fpath):
            sys.exit(f"{label} not found: {fpath}")

    fai_file = f"{fasta_path}.fai"
    if not os.path.exists(fai_file):
        tee.write("Indexing FASTA...\n")
        subprocess.run(f"samtools faidx {fasta_path}", shell=True, check=True)

    # ── Build Rscript command ─────────────────────────────────────────
    r_script = os.path.join(prefix, "scripts", "cips_uORF.R")
    if not os.path.exists(r_script):
        sys.exit(f"R script not found: {r_script}")

    cmd = (
        f"Rscript --vanilla {r_script} "
        f"{gtf_path} {fasta_path} {psite_file} {output_prefix} "
        f"{min_inframe_counts} {min_inframe_perc} {min_psite_perc}"
    )
    if gene_desc and os.path.exists(gene_desc):
        cmd += f" {gene_desc}"

    tee.write("=== CiPS uORF Analysis ===\n")
    tee.write(f"GTF:           {gtf_path}\n")
    tee.write(f"FASTA:         {fasta_path}\n")
    tee.write(f"P-site file:   {psite_file}\n")
    tee.write(f"Output prefix: {output_prefix}\n")
    tee.write(f"Min in-frame counts: {min_inframe_counts}\n")
    tee.write(f"Min in-frame %:      {min_inframe_perc}\n")
    tee.write(f"Min P-site %:        {min_psite_perc}\n")
    if gene_desc:
        tee.write(f"Gene desc:     {gene_desc}\n")
    tee.write("\n")

    # ── Run R script ──────────────────────────────────────────────────
    tee.write("Running CiPS uORF analysis (R)...\n")
    subprocess.run(cmd, shell=True, check=True)

    tee.write("\nCiPS analysis complete.\n")
    tee.write(f"Output files:\n")
    tee.write(f"  {output_prefix}_minimum_uORF.xlsx\n")
    tee.write(f"  {output_prefix}_uORF.xlsx\n")
    tee.write(f"  {output_prefix}_uORF_dedup.xlsx\n")
    tee.write(f"  {output_prefix}_uORF_overlaps.xlsx\n")
