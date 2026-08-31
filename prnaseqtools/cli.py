#!/usr/bin/env python3
"""
pRNASeqTools - Main CLI entry point using argparse (stdlib, zero-dependency).
Mirrors the MooseX::App::Command structure of the original Perl version.
"""

import os
import sys
import time
import argparse
from pathlib import Path

from prnaseqtools import __version__
from prnaseqtools.precheck import check_dependencies
from prnaseqtools.logging_setup import setup_logging, finalize_logging

# ── Global state ─────────────────────────────────────────────────────────
PREFIX = str(Path(__file__).resolve().parent.parent)
START_TIME = int(time.time())
TEE = None


# ── Common argument groups ───────────────────────────────────────────────
def add_common_args(parser, control_required=True):
    """Add arguments shared across all subcommands."""
    parser.add_argument('--outdir', '-o', default='./out', help='Output directory.')
    parser.add_argument('--genome', '-g', default='ath',
                        help='Genome: ath, osa, b73, gma, smo, bra, w22')
    parser.add_argument('--thread', '-t', default=4, type=int, help='Threads.')
    parser.add_argument('--adaptor', '-a', default=None, help="3' adaptor sequence.")
    if control_required:
        parser.add_argument('--control', '-c', required=True,
                            help='Control: name=file1+file2...')
    else:
        parser.add_argument('--control', '-c', default=None,
                            help='Control: name=file1+file2...')
    parser.add_argument('--treatment', '-p', default=None, action='append',
                        help='Treatment: name=file1+file2... (repeat for multiple)')
    parser.add_argument('--mask', default=None, help='Mask sequences (fasta).')
    return parser


def add_mapping_args(parser):
    """Add nomapping/mappingonly arguments."""
    parser.add_argument('--no-mapping', action='store_true', help='Only statistics.')
    parser.add_argument('--mapping-only', action='store_true', help='Only mapping.')
    return parser


# ── Main argument parser ─────────────────────────────────────────────────
def build_parser():
    """Build the full argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        prog='pRNASeqTools',
        description='Integrated High-throughput Sequencing Data Analysis for Plant (Python3)',
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('--auto-install', action='store_true', default=True,
                        help='Auto-install missing dependencies (default)')
    parser.add_argument('--no-auto-install', action='store_false', dest='auto_install',
                        help='Disable automatic dependency installation')

    sub = parser.add_subparsers(dest='mode', title='Analysis modes')

    # srna
    p = sub.add_parser('srna', help='Small RNA-seq analysis')
    add_common_args(p)
    add_mapping_args(p)
    p.add_argument('--mmap', default='u', help='Multi-map: u, n, f, r')
    p.add_argument('--foldchange', default=1.5, type=float, help='FC threshold')
    p.add_argument('--pvalue', default=0.01, type=float, help='P-value threshold')
    p.add_argument('--norm', default='rRNA,total', help='Normalization')
    p.add_argument('--binsize', default=100, type=int, help='Window size')
    p.add_argument('--mode_srna', default='bulk', dest='run_mode', help='bulk or sc')
    p.add_argument('--pattern', default='NNNNNNNNCA', help='UMI pattern')
    p.add_argument('--promoter', default=1000, type=int, help='Promoter length')
    p.add_argument('--spike-in', default=None, help='Spike-in sequences')

    # mrna
    p = sub.add_parser('mrna', help='mRNA-seq analysis')
    add_common_args(p)
    p.add_argument('--total', action='store_true', help='Total RNA-seq (incl. ncRNA)')
    p.add_argument('--mode_mrna', default='whole', dest='run_mode',
                   choices=['whole', 'mapping-only', 'bam', 'count-table'],
                   help='whole=mapping+count+DE, mapping-only=mapping+count, '
                        'bam=from BAM, count-table=from count table')
    p.add_argument('--foldchange', default=2.0, type=float, help='FC threshold')
    p.add_argument('--pvalue', default=0.01, type=float, help='P-value threshold')
    p.add_argument('--fdr', default=1.0, type=float, help='FDR threshold')
    p.add_argument('--deseq2norm', default='DESeq2', dest='deseq2_norm',
                   help='DESeq2 or RPM')
    p.add_argument('--seqstrategy', default=None, dest='seq_strategy',
                   help='single or paired')
    p.add_argument('--genomesize', default=10, type=int, dest='genome_size',
                   help='genomeSAindexNbases')

    # lncrna (lncRNA-seq: de novo transcript assembly + coding/lncRNA classification)
    p = sub.add_parser('lncrna', help='lncRNA-seq analysis '
                                       '(de novo assembly + classification + DE)')
    add_common_args(p)
    p.add_argument('--mode_lncrna', default='whole', dest='run_mode',
                   choices=['whole', 'mapping-only', 'assemble-only',
                            'de-only', 'count-table'],
                   help='whole=mapping+assemble+classify+quantify+DE, '
                        'mapping-only=STAR+BAM only, '
                        'assemble-only=STAR→StringTie→classify only, '
                        'de-only=skip mapping/assembly, '
                        'count-table=use existing GTF + BAM/count tables')
    p.add_argument('--foldchange', default=2.0, type=float, help='FC threshold')
    p.add_argument('--pvalue', default=0.01, type=float, help='P-value threshold')
    p.add_argument('--fdr', default=1.0, type=float, help='FDR threshold')
    p.add_argument('--seqstrategy', default=None, dest='seq_strategy',
                   help='single or paired')
    p.add_argument('--genomesize', default=10, type=int, dest='genome_size',
                   help='genomeSAindexNbases')
    p.add_argument('--classifier', default='feelnc',
                   choices=['feelnc', 'plek2', 'fallback'],
                   help='Coding/lncRNA classifier. '
                        'feelnc=FEELnc (bioconda, recommended for plants), '
                        'plek2=PLEK2 (deep-learning, plant model), '
                        'fallback=ORF heuristic (no external deps)')
    p.add_argument('--feelnc-dir', dest='feelnc_dir', default=None,
                   help='Path to FEELnc resource dir '
                        '(only needed if FEELnc was NOT installed via bioconda; '
                        'conda-installed FEELnc auto-finds its own resources)')
    p.add_argument('--plek2-dir', dest='plek2_dir', default=None,
                   help='Path to PLEK2 install dir '
                        '(default: ~/PLEK2)')
    p.add_argument('--min-fpkm', dest='min_fpkm', default=0.5, type=float,
                   help='Min FPKM for StringTie transcript filtering (0.5)')
    p.add_argument('--min-tx-len', dest='min_tx_len', default=200, type=int,
                   help='Min transcript length (bp) for lncRNA (200)')
    p.add_argument('--gtf', default=None,
                   help='Custom GTF (use in assemble-only / de-only '
                        'to skip merging)')

    # sc (single-cell RNA-seq)
    p = sub.add_parser('sc', help='Single-cell RNA-seq analysis')
    add_common_args(p)
    p.add_argument('--mode_sc', default='whole', dest='run_mode',
                   choices=['whole', 'mapping-only', 'count-table'],
                   help='whole=mapping+count+analysis, '
                        'mapping-only=mapping+count, '
                        'count-table=from count table')
    p.add_argument('--seqstrategy', default=None, dest='seq_strategy',
                   help='single or paired')
    p.add_argument('--genomesize', default=10, type=int, dest='genome_size',
                   help='genomeSAindexNbases')
    p.add_argument('--mincells', default=3, type=int, dest='min_cells',
                   help='Min cells per gene (Seurat QC)')
    p.add_argument('--minfeatures', default=200, type=int, dest='min_features',
                   help='Min features per cell (Seurat QC)')
    p.add_argument('--maxfeatures', default=5000, type=int, dest='max_features',
                   help='Max features per cell (Seurat QC)')
    p.add_argument('--pctmt', default=20, type=float, dest='pct_mt',
                   help='Max percent mitochondrial content')
    p.add_argument('--npcs', default=30, type=int, dest='n_pcs',
                   help='Number of PCs for dimensionality reduction')
    p.add_argument('--nclusters', default=0, type=int, dest='n_clusters',
                   help='Number of clusters (0 = auto)')
    p.add_argument('--resolution', default=0.5, type=float,
                   help='Clustering resolution (0.4-1.2)')
    p.add_argument('--markerminpct', default=0.25, type=float,
                   dest='marker_min_pct',
                   help='Min percent cells expressing marker')
    p.add_argument('--markerlogfc', default=0.25, type=float,
                   dest='marker_logfc',
                   help='Min log2 fold change for markers')
    p.add_argument('--pseudotime', action='store_true',
                   help='Run pseudotime analysis (Monocle3)')
    p.add_argument('--doublet', default=0, type=float, dest='doublet_rate',
                   help='Expected doublet rate for DoubletFinder (0 = skip)')
    p.add_argument('--integration', default='seurat',
                   choices=['seurat', 'harmony', 'none'],
                   help='Batch integration method: seurat=Seurat anchors, '
                        'harmony=Harmony, none=no batch correction')

    # degradome
    p = sub.add_parser('degradome', help='Degradome-seq analysis')
    add_common_args(p)
    add_mapping_args(p)
    p.add_argument('--targets', default='all', help='Transcript list for CRI')
    p.add_argument('--sirnas', default='none', dest='si_rnas',
                   help='Extra siRNAs for targets')

    # phasi
    p = sub.add_parser('phasi', help='phasiRNA analysis')
    add_common_args(p)
    p.add_argument('--no-mapping', action='store_true', help='Only statistics')
    p.add_argument('--mmap', default='u', help='Multi-map method')
    p.add_argument('--norm', default='rRNA,total', help='Normalization')
    p.add_argument('--period', default=21, type=int, help='Period size (19-26)')
    p.add_argument('--binsize', default=100, type=int, help='Window size')
    p.add_argument('--phasingscore', default=50, type=float, help='Score cutoff')

    # tt
    p = sub.add_parser('tt', help='Truncation/tailing analysis')
    add_common_args(p)
    p.add_argument('--mmap', default='u', help='Multi-map method')

    # ribo (RIBO Taper pipeline)
    p = sub.add_parser('ribo', help='Ribo-seq analysis (RIBO Taper pipeline)')
    add_common_args(p, control_required=False)
    p.add_argument('--rna-control', '-rc', required=True,
                   help='RNA-seq control: name=file1+file2...')
    p.add_argument('--rna-treatment', '-rt', default=None, action='append',
                   help='RNA-seq treatment: name=file1+file2... (repeatable)')
    p.add_argument('--ribo-control', '-bc', required=True,
                   help='Ribo-seq control: name=file1+file2...')
    p.add_argument('--ribo-treatment', '-bt', default=None, action='append',
                   help='Ribo-seq treatment: name=file1+file2... (repeatable)')
    p.add_argument('--contam', default=None,
                   help='Contamination fasta file(s) for Bowtie2 index '
                        '(default: reference/{genome}_contam4.fa)')
    p.add_argument('--ribo-len', default='24,25,26,27,28',
                   help='Ribo-seq read lengths (comma-separated); '
                        'count must match cutoffs')
    p.add_argument('--cutoffs', default='8,9,10,11,12',
                   help='RIBO Taper cutoffs (comma-separated); '
                        'count must match ribo-len')
    p.add_argument('--ribotaper', default=os.path.expanduser('~/software/ribotaper/bin'),
                   help='Path to RIBO Taper installation directory')
    p.add_argument('--ribotaper-env', default='ribotaper',
                   help='Conda environment name for RIBO Taper')
    p.add_argument('--tpm-threshold', default=0, type=float,
                   help='Mean TPM threshold for isoform filtering (default: 0)')
    p.add_argument('--restart-step', type=int, default=None,
                   help='Force restart from this step (1-12), '
                        'overrides auto-detection')

    # cips (CiPS uORF analysis)
    p = sub.add_parser('cips', help='CiPS uORF analysis (translated uORF detection)')
    add_common_args(p, control_required=False)
    p.add_argument('--gtf', default=None,
                   help='Expressed GTF (default: <genome>_expressed.gtf)')
    p.add_argument('--fasta', default=None,
                   help='Reference genome FASTA (default: reference/<genome>_chr_all.fasta)')
    p.add_argument('--psite', required=True,
                   help='P-site count file (count chr start strand)')
    p.add_argument('--output', default='cips',
                   help='Output prefix (default: cips)')
    p.add_argument('--min-inframe-counts', default=10, type=int,
                   help='Min in-frame Ribo-seq counts (default: 10)')
    p.add_argument('--min-inframe-perc', default=50, type=float,
                   help='Min in-frame count percentage (default: 50)')
    p.add_argument('--min-psite-perc', default=30, type=float,
                   help='Min in-frame P-site %% for longer uORFs (default: 30)')
    p.add_argument('--gene-desc', default=None,
                   help='Gene description Excel for annotation')

    # chip
    p = sub.add_parser('chip', help='ChIP-seq analysis')
    add_common_args(p, control_required=False)
    add_mapping_args(p)
    p.add_argument('--peak-caller', default='genrich', choices=['genrich', 'macs3'],
                   help='Peak caller: genrich or macs3 (default: genrich)')
    p.add_argument('--genome-size', default=None,
                   help='Effective genome size for MACS3 (e.g. 1.35e8 for ath)')
    p.add_argument('--auc', default=20, type=float, help='AUC threshold (Genrich only)')
    p.add_argument('--qvalue', default=1.0, type=float, help='Q-value threshold')
    p.add_argument('--pvalue', default=0.01, type=float, help='P-value threshold')
    p.add_argument('--seqstrategy', default='paired', dest='seq_strategy',
                   help='single or paired')
    p.add_argument('--tss-distance', default=3000, type=int,
                   help='TSS distance for ChIPseeker annotation (default: 3000)')


    # atac
    p = sub.add_parser('atac', help='ATAC-seq analysis')
    add_common_args(p, control_required=False)
    add_mapping_args(p)
    p.add_argument('--peak-caller', default='genrich', choices=['genrich', 'macs3'],
                   help='Peak caller: genrich or macs3 (default: genrich)')
    p.add_argument('--genome-size', default=None,
                   help='Effective genome size for MACS3 (e.g. 1.35e8 for ath)')
    p.add_argument('--auc', default=20, type=float, help='AUC threshold (Genrich only)')
    p.add_argument('--qvalue', default=1.0, type=float, help='Q-value threshold')
    p.add_argument('--pvalue', default=0.01, type=float, help='P-value threshold')
    p.add_argument('--treatment2', default=None,
                   help='Second ATAC group for MACS3 bdgdiff: name=file1+file2...')
    p.add_argument('--tss-distance', default=3000, type=int,
                   help='TSS distance for ChIPseeker annotation (default: 3000)')

    # wgbs
    p = sub.add_parser('wgbs', help='WGBS-seq analysis')
    add_common_args(p)
    add_mapping_args(p)
    p.add_argument('--binsize', default=100, type=int, help='Window size')
    p.add_argument('--minc', default=4, type=int, dest='min_c',
                   help='Min reads per cytosine')

    # clip
    p = sub.add_parser('clip', help='CLIP-seq analysis')
    add_common_args(p)
    add_mapping_args(p)
    p.add_argument('--foldchange', default=2.0, type=float, help='FC threshold')
    p.add_argument('--pvalue', default=0.05, type=float, help='P-value threshold')

    # ts
    p = sub.add_parser('ts', help='TS-CLIP-seq analysis')
    add_common_args(p)
    add_mapping_args(p)
    p.add_argument('--foldchange', default=2.0, type=float, help='FC threshold')
    p.add_argument('--pvalue', default=0.05, type=float, help='P-value threshold')

    # ribometh
    p = sub.add_parser('ribometh', help='RiboMeth-seq analysis')
    add_common_args(p)
    add_mapping_args(p)
    p.add_argument('--reference', default='genome', dest='ref_file',
                   help='Reference transcriptome fasta')
    p.add_argument('--readlength', default=50, type=int, help='Raw read length')
    p.add_argument('--coverage', default=1000, type=int, help='Min coverage')
    p.add_argument('--adaptor2', default='1', help='Adaptor for read 2')

    # risi
    p = sub.add_parser('risi', help='risiRNA analysis')
    add_common_args(p)
    p.add_argument('--mmap', default='u', help='Multi-map method')
    p.add_argument('--no-mapping', action='store_true', help='Only statistics')
    p.add_argument('--foldchange', default=1.5, type=float, help='FC threshold')
    p.add_argument('--pvalue', default=0.01, type=float, help='P-value threshold')
    p.add_argument('--mapping-only', action='store_true', help='Only mapping')
    p.add_argument('--norm', default='total', help='Normalization')
    p.add_argument('--binsize', default=10, type=int, help='Window size')

    # tf
    p = sub.add_parser('tf', help='Two-factor DE analysis')
    p.add_argument('--outdir', '-o', default='./out', help='Output directory.')
    p.add_argument('--genome', '-g', default='ath',
                   help='Genome: ath, osa, b73, gma, smo, bra, w22')
    p.add_argument('--thread', '-t', default=4, type=int, help='Threads.')
    p.add_argument('--adaptor', '-a', default=None, help="3' adaptor sequence.")
    p.add_argument('--control', '-c', default=None,
                   help='groupName=label1,N1,label2,N2  '
                        '(e.g. Col=0,2,7,2 for mrna/srna; '
                        'WT=input,3,IP,3 for chip)')
    p.add_argument('--treatment', '-p', default=None,
                   help='groupName=label1,N1,label2,N2  '
                        '(e.g. Mut=0,2,7,2 for mrna/srna; '
                        'KO=input,3,IP,3 for chip)')
    p.add_argument('--foldchange', default=1.5, type=float, help='FC threshold')
    p.add_argument('--pvalue', default=0.05, type=float, help='P-value threshold')
    p.add_argument('--mode_tf', default='mrna', dest='run_mode',
                   help='mrna, srna, or chip')
    p.add_argument('--norm', default='rRNA,total', help='Normalization')
    p.add_argument('--binsize', default=100, type=int, help='Window size')
    p.add_argument('--deseq2norm', default='DESeq2', dest='deseq2_norm',
                   help='DESeq2 or RPM')
    # ChIP bdgdiff options
    p.add_argument('--genome-size', default=None,
                   help='Effective genome size for MACS3 (e.g. 1.35e8 for ath)')
    p.add_argument('--cutoff', default=3, type=float,
                   help='log2 fold-change cutoff for bdgdiff (default: 3)')
    p.add_argument('--qvalue', default=1.0, type=float, help='Q-value threshold')
    p.add_argument('--seqstrategy', default='paired', dest='seq_strategy',
                   help='single or paired')
    p.add_argument('--tss-distance', default=3000, type=int,
                   help='TSS distance for ChIPseeker annotation (default: 3000)')

    return parser


# ── Dispatch table ───────────────────────────────────────────────────────
MODE_RUNNERS = {
    'srna':       'prnaseqtools.modes.srna',
    'mrna':       'prnaseqtools.modes.mrna',
    'lncrna':     'prnaseqtools.modes.lncrna',
    'sc':         'prnaseqtools.modes.sc',
    'degradome':  'prnaseqtools.modes.degradome',
    'phasi':      'prnaseqtools.modes.phasi',
    'tt':         'prnaseqtools.modes.tt',
    'ribo':       'prnaseqtools.modes.ribo',
    'cips':       'prnaseqtools.modes.cips',
    'chip':       'prnaseqtools.modes.chip',
    'atac':       'prnaseqtools.modes.atac',
    'wgbs':       'prnaseqtools.modes.wgbs',
    'clip':       'prnaseqtools.modes.clip',
    'ts':         'prnaseqtools.modes.ts',
    'ribometh':   'prnaseqtools.modes.ribometh',
    'risi':       'prnaseqtools.modes.risi',
    'tf':         'prnaseqtools.modes.tf',
}


def main():
    """Main entry point."""
    global TEE

    # Log files are ephemeral by default — only keep them when a real
    # analysis pipeline runs to completion.  Empty invocations, --help,
    # --version, and argument errors all clean up after themselves.
    _cleanup_log = True

    # Setup logging
    TEE = setup_logging(PREFIX, START_TIME)

    try:
        # Parse all args in one pass (two-pass breaks subparser dispatch)
        parser = build_parser()
        args = parser.parse_args()

        if not args.mode:
            parser.print_help()
            sys.exit(0)

        # ── Real analysis confirmed — keep the log ──────────────────
        _cleanup_log = False

        auto_install = getattr(args, 'auto_install', True)
        check_dependencies(auto_install=auto_install, mode=args.mode)

        # Convert namespace to dict
        opts = vars(args)
        opts['prefix'] = PREFIX
        opts['_start_time'] = START_TIME

        # Normalize option names (Click -> argparse differences)
        _normalize_opts(opts)

        # Dispatch to mode
        module_name = MODE_RUNNERS.get(args.mode)
        if not module_name:
            TEE.write(f"Unknown mode: {args.mode}\n")
            sys.exit(1)

        import importlib
        mod = importlib.import_module(module_name)
        mod.run(opts)

    except SystemExit:
        raise
    except Exception as e:
        TEE.write(f"Error: {e}\n")
        import traceback
        traceback.print_exc(file=TEE if TEE else sys.stderr)
        sys.exit(1)
    finally:
        if TEE:
            finalize_logging(TEE, START_TIME, cleanup=_cleanup_log)


def _normalize_opts(opts):
    """Normalize argparse option names to match internal naming convention."""
    # Map argparse dest names to expected internal names
    renames = {
        'no_mapping': 'no_mapping',
        'mapping_only': 'mapping_only',
    }
    # The rest should already match since we used dest= appropriately


if __name__ == '__main__':
    main()