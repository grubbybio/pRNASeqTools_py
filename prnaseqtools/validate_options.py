"""
Option validation for pRNASeqTools.
Mirrors the Perl validate_options.pm module.
"""

import os
import sys
import shutil
from pathlib import Path


def validate_options(opts):
    """Validate and normalize command-line options. Raises SystemExit on error."""

    outdir = opts.get('outdir', './out')
    if os.path.exists(outdir):
        sys.exit('Output directory exists! Please specify another output directory!')
    os.makedirs(outdir, exist_ok=True)

    # Move log file to output directory
    log_file = f"log_{opts.get('_start_time', 0)}.txt"
    if os.path.exists(log_file):
        shutil.move(log_file, os.path.join(outdir, log_file))

    # Change to output directory
    os.chdir(outdir)

    # Adaptor shortcuts
    adaptor_map = {
        '1': 'AGATCGGAAGAGC',
        '2': 'TGGAATTCTCGGG',
        '3': 'CTGTCTCTTATAC',
    }
    if opts.get('adaptor') and opts['adaptor'] in adaptor_map:
        opts['adaptor'] = adaptor_map[opts['adaptor']]

    if opts.get('adaptor2') and opts['adaptor2'] in adaptor_map:
        opts['adaptor2'] = adaptor_map[opts['adaptor2']]
    elif opts.get('adaptor2') == '3':
        opts['adaptor2'] = 'CTGTCTCTTATA'

    # Thread validation
    thread = opts.get('thread')
    if thread is not None:
        thread_str = str(thread)
        if not thread_str.isdigit() or int(thread) < 1:
            sys.exit('Please use appropriate threads!')

    # P-value validation
    pvalue = opts.get('pvalue')
    if pvalue is not None and pvalue > 1:
        sys.exit('Please use an appropriate P value!')

    # FDR validation
    fdr = opts.get('fdr')
    if fdr is not None and fdr > 1:
        sys.exit('Please use an appropriate FDR value!')

    # Fold change validation
    foldchange = opts.get('foldchange')
    if foldchange is not None and foldchange < 1.5:
        sys.exit('Please use an appropriate fold change!')

    # Multi-map method validation
    mmap = opts.get('mmap')
    if mmap is not None and mmap not in 'ufrn':
        sys.exit('Please use a supported strategy for mapping!')

    # Length validation
    length = opts.get('length')
    if length is not None and (length < 18 or length > 42):
        sys.exit('Please specify an appropriate length of preferred small RNAs!')

    # nomapping + mappingonly conflict
    if opts.get('no_mapping') and opts.get('mapping_only'):
        sys.exit('Parameter conflict: nomapping and mappingonly!')

    # DESeq2Norm validation
    deseq2_norm = opts.get('deseq2_norm')
    if deseq2_norm is not None and deseq2_norm not in ('DESeq2', 'RPM'):
        sys.exit('Method not supported!')

    # Sequencing strategy validation
    seq_strategy = opts.get('seq_strategy')
    if seq_strategy is not None and seq_strategy not in ('single', 'paired'):
        sys.exit('Please specify an appropriate sequencing strategy!')

    # Mask validation
    mask = opts.get('mask')
    if mask is not None and not mask.endswith(('.fasta', '.fa')):
        sys.exit('Please specify a fasta file for mask!')

    # Spike-in validation
    spike_in = opts.get('spike_in')
    if spike_in is not None and not spike_in.endswith(('.fasta', '.fa')):
        sys.exit('Please specify a fasta file for spike-in!')

    # Style validation
    style = opts.get('style')
    if style is not None and style not in ('histone', 'factor', 'tss'):
        sys.exit('Please select the correct style!')

    # Targets validation
    targets = opts.get('targets')
    if targets is not None and targets != 'all' and not os.path.exists(targets):
        sys.exit('Cannot find the target list!')

    # BAM mode requires seq strategy
    run_mode = opts.get('run_mode')
    if run_mode == 3 and not opts.get('seq_strategy'):
        sys.exit('Please provide the library type when input files are in the bam format!')

    # siRNA validation
    si_rnas = opts.get('si_rnas')
    if si_rnas is not None and si_rnas != 'none' and not os.path.exists(si_rnas):
        sys.exit('Cannot find the siRNA list in fasta format!')

    # Store start_time for use by downstream code
    from prnaseqtools.cli import START_TIME
    opts['_start_time'] = START_TIME

    return opts
