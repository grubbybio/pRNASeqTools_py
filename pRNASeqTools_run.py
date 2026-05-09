#!/usr/bin/env python3
"""
pRNASeqTools - Integrated High-throughput Sequencing Data Analysis for Plant
Python3 rewrite. Entry point script.

Usage: pRNASeqTools_run <mode> [OPTIONS]

Modes: srna, mrna, degradome, phasi, tt, ribo, chip, atac, wgbs, clip, ts, ribometh, risi, tf
"""

import sys
import os

_package_dir = os.path.dirname(os.path.abspath(__file__))
if _package_dir not in sys.path:
    sys.path.insert(0, _package_dir)

from prnaseqtools.cli import main

if __name__ == '__main__':
    main()
