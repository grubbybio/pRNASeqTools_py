"""
Option validation for pRNASeqTools.
Mirrors the Perl validate_options.pm module.
"""

import os
import sys
import shutil

from prnaseqtools.input_parser import ADAPTOR_ALIASES


def validate_options(opts):
    """Validate and normalize command-line options. Raises SystemExit on error."""

    outdir = opts.get('outdir', './out')
    mode = opts.get('mode', '')
    if os.path.exists(outdir):
        if mode == 'ribo':
            # Ribo-seq auto-resume: output dir is allowed, but it must
            # contain log files from a previous run to be a valid resume
            # target.  If no previous logs exist this is likely an
            # unrelated directory → refuse to continue.
            _start_time = opts.get('_start_time', 0)
            _current_log = f"log_{_start_time}.txt"
            _prev_logs = [
                f for f in os.listdir(outdir)
                if f.startswith('log_') and f.endswith('.txt')
                and f != _current_log
            ]
            if not _prev_logs:
                sys.exit(
                    f"Output directory '{outdir}' exists but contains "
                    f"no previous log files.\n"
                    f"  This does not appear to be a valid Ribo-seq "
                    f"resume target.\n"
                    f"  Please specify a different --outdir or run the "
                    f"pipeline from scratch."
                )
            # Valid resume target — warn user before continuing
            import time
            print(
                f"Output directory '{outdir}' already exists.\n"
                f"  Found {len(_prev_logs)} previous log file(s) → "
                f"resuming from the last completed step.\n"
                f"  Press Ctrl+C now to abort, or wait 5 seconds "
                f"to continue...",
                flush=True
            )
            try:
                time.sleep(5)
            except KeyboardInterrupt:
                sys.exit("Aborted by user.")
        else:
            sys.exit('Output directory exists! Please specify another output directory!')
    os.makedirs(outdir, exist_ok=True)

    # Move log file to output directory
    log_file = f"log_{opts.get('_start_time', 0)}.txt"
    if os.path.exists(log_file):
        shutil.move(log_file, os.path.join(outdir, log_file))

    # Change to output directory
    os.chdir(outdir)

    # Adaptor shortcuts (ADAPTOR_ALIASES shared from input_parser)
    adaptor = opts.get('adaptor')
    if adaptor and adaptor.lower() in ADAPTOR_ALIASES:
        opts['adaptor'] = ADAPTOR_ALIASES[adaptor.lower()]

    adaptor2 = opts.get('adaptor2')
    if adaptor2 and adaptor2.lower() in ADAPTOR_ALIASES:
        opts['adaptor2'] = ADAPTOR_ALIASES[adaptor2.lower()]

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

    # nomapping + mappingonly conflict (ribo mode uses auto-resume instead)
    if mode != 'ribo' and opts.get('no_mapping') and opts.get('mapping_only'):
        sys.exit('Parameter conflict: nomapping and mappingonly!')

    # Ribo-seq: validate ribo-len and cutoffs have the same number of items
    if mode == 'ribo':
        _rl = opts.get('ribo_len', '')
        _co = opts.get('cutoffs', '')
        _rl_n = len([x for x in _rl.split(',') if x.strip()])
        _co_n = len([x for x in _co.split(',') if x.strip()])
        if _rl_n != _co_n:
            sys.exit(
                f"ribo-len ({_rl_n} items) and cutoffs ({_co_n} items) "
                f"must have the same number of comma-separated values!\n"
                f"  ribo-len:  {_rl}\n"
                f"  cutoffs:   {_co}"
            )

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

    # Targets validation
    targets = opts.get('targets')
    if targets is not None and targets != 'all' and not os.path.exists(targets):
        sys.exit('Cannot find the target list!')

    # BAM mode requires seq strategy
    run_mode = opts.get('run_mode')
    if run_mode == 'bam' and not opts.get('seq_strategy'):
        sys.exit('Please provide the library type when input files are in the bam format!')

    # siRNA validation
    si_rnas = opts.get('si_rnas')
    if si_rnas is not None and si_rnas != 'none' and not os.path.exists(si_rnas):
        sys.exit('Cannot find the siRNA list in fasta format!')

    # Single-cell specific validation
    if mode == 'sc':
        # Validate Seurat QC parameters
        min_features = opts.get('min_features')
        max_features = opts.get('max_features')
        if min_features is not None and max_features is not None:
            if min_features >= max_features:
                sys.exit('minfeatures must be less than maxfeatures!')

        pct_mt = opts.get('pct_mt')
        if pct_mt is not None and (pct_mt < 0 or pct_mt > 100):
            sys.exit('pctmt must be between 0 and 100!')

        n_pcs = opts.get('n_pcs')
        if n_pcs is not None and n_pcs < 1:
            sys.exit('npcs must be at least 1!')

        resolution = opts.get('resolution')
        if resolution is not None and (resolution < 0.1 or resolution > 2.0):
            sys.exit('resolution should be between 0.1 and 2.0!')

        run_mode_sc = opts.get('run_mode')
        if run_mode_sc not in ('whole', 'mapping-only', 'count-table'):
            sys.exit('mode_sc must be one of: whole, mapping-only, count-table!')

        # count-table mode requires pre-existing count files
        if run_mode_sc == 'count-table':
            control_val = opts.get('control', '')
            if not control_val:
                sys.exit('count-table mode requires --control parameter!')

    return opts
