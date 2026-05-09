"""
Dependency checker for pRNASeqTools — with auto-install support.
Checks that all required external tools are installed, and can
automatically install missing ones via mamba/conda/pip/R/git.
"""

import subprocess
import sys
import os
import re
from pathlib import Path


def _get_version(cmd, pattern=None):
    """Run a command and extract version using regex pattern."""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True,
                                text=True, timeout=30)
        output = result.stdout + result.stderr
        if pattern:
            m = re.search(pattern, output)
            if m:
                return m.group(1) if m.lastindex else m.group(0)
        return output.strip()
    except Exception:
        return None


# ═══════════════════════════════════════════════════════════════════════════
# 1. Check definitions — maps binary name → version check + install info
# ═══════════════════════════════════════════════════════════════════════════

_CHECK_DEFS = {
    'cutadapt': {
        'cmd': 'cutadapt --version',
        'pattern': r'^(\d+\.\d+)',
        'required': True,
        'msg': 'Adapter trimming tool',
    },
    'samtools': {
        'cmd': 'samtools --version',
        'pattern': r'samtools (\d+\.\d+)',
        'required': True,
        'msg': 'SAM/BAM processing',
    },
    'bowtie': {
        'cmd': 'bowtie --version',
        'pattern': r'bowtie.+ (\d+\.\d+\.\d+)',
        'required': True,
        'msg': 'Short-read aligner (sRNA)',
    },
    'bowtie2': {
        'cmd': 'bowtie2 --version',
        'pattern': r'version (\d+\.\d+\.\d+)',
        'required': True,
        'msg': 'Short-read aligner (ChIP/ATAC)',
    },
    'featureCounts': {
        'cmd': 'featureCounts -v 2>&1',
        'pattern': r'v(\d+\.\d+\.\d+)',
        'required': True,
        'msg': 'Read counting',
    },
    'ShortStack': {
        'cmd': 'ShortStack --version',
        'pattern': r'ShortStack.+[34]',
        'required': True,
        'msg': 'Small RNA alignment',
        'mode_only': ['srna', 'phasi', 'tt', 'risi'],
    },
    'bedtools': {
        'cmd': 'bedtools --version',
        'pattern': r'bedtools v(\d+\.\d+\.\d+)',
        'required': True,
        'msg': 'Genome interval arithmetic',
    },
    'R': {
        'cmd': 'R --version',
        'pattern': r'R version (\d+\.\d+\.\d+)',
        'required': True,
        'msg': 'R statistical computing',
    },
    'STAR': {
        'cmd': 'STAR --version',
        'pattern': r'(\d+\.\d+\.\d+)',
        'required': True,
        'msg': 'RNA-seq aligner',
    },
    'fasterq-dump': {
        'cmd': 'fasterq-dump -h',
        'pattern': r'fasterq-dump.+(\d+\.\d+\.\d+)',
        'required': False,  # only needed if SRR input
        'msg': 'SRA data download',
    },
    'gffread': {
        'cmd': 'gffread --version 2>&1',
        'pattern': r'(\d+\.\d+\.\d+)',
        'required': True,
        'msg': 'GFF/GTF conversion',
    },
    'deeptools': {
        'cmd': 'deeptools --version 2>&1',
        'pattern': r'(\d+\.\d+\.\d+)',
        'required': True,
        'msg': 'Coverage track generation',
    },
    'Genrich': {
        'cmd': 'Genrich --version 2>&1',
        'pattern': r'version (\d+\.\d+\.\d+)',
        'required': False,
        'msg': 'ChIP/ATAC peak calling',
        'mode_only': ['chip', 'atac'],
    },
    'umi_tools': {
        'cmd': 'umi_tools -v',
        'pattern': r'(\d+\.\d+\.\d+)',
        'required': False,
        'msg': 'UMI extraction (single-cell)',
        'mode_only': ['srna'],
    },
    'bedGraphToBigWig': {
        'cmd': 'bedGraphToBigWig 2>&1',
        'pattern': r'bedGraphToBigWig v (\d+)',
        'required': True,
        'msg': 'UCSC bigWig conversion',
    },
    'Bismark': {
        'cmd': 'bismark --version 2>&1',
        'pattern': r'Bismark Version: v(\d+\.\d+\.\d+)',
        'required': False,
        'msg': 'Bisulfite alignment (WGBS)',
        'mode_only': ['wgbs'],
    },
    'clipper': {
        'cmd': 'clipper -h 2>&1',
        'pattern': r'^Usage',
        'required': False,
        'msg': 'CLIP-seq peak caller',
        'mode_only': ['clip', 'ts'],
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# 2. Main check function
# ═══════════════════════════════════════════════════════════════════════════

def check_dependencies(auto_install=True, mode=None, interactive=True):
    """
    Check all external dependencies.  Optionally auto-install missing ones.

    Args:
        auto_install: if True, attempt automatic installation of missing tools
        mode:        current analysis mode (e.g. 'srna'); skips mode-only tools
        interactive: if True, prompt before installing (ignored if not auto_install)

    Raises SystemExit if a required tool is missing and cannot be installed.
    """
    tee = sys.stderr
    try:
        from prnaseqtools.cli import TEE
        if TEE:
            tee = TEE
    except (ImportError, AttributeError):
        pass

    prefix = str(Path(__file__).resolve().parent.parent)

    tee.write("\nChecking dependent software...\n")

    # ── Phase 1: check all tools, collect missing ─────────────────────
    missing = []
    for name, cfg in _CHECK_DEFS.items():
        # Skip mode-only tools when mode doesn't need them
        mode_only = cfg.get('mode_only')
        if mode_only and mode not in mode_only:
            continue

        version = _get_version(cfg['cmd'], cfg.get('pattern'))
        if version:
            tee.write(f"  {name} version {version}\n")
        else:
            if cfg.get('required', True):
                tee.write(f"  {name}: MISSING — {cfg['msg']}\n")
                missing.append((name, cfg))
            else:
                tee.write(f"  {name}: not found (optional)\n")

    # ── Phase 2: auto-install missing tools ───────────────────────────
    if missing and auto_install:
        from prnaseqtools.auto_install import DEPENDENCY_REGISTRY, install_missing

        # Enrich missing list with registry install info
        enriched = []
        for name, cfg in missing:
            info = DEPENDENCY_REGISTRY.get(name, {})
            if info:
                enriched.append((name, info))
            else:
                tee.write(f"  {name}: no auto-install recipe available\n")
                enriched.append((name, {'install_msg': cfg['msg']}))

        installed_count, failed_list = install_missing(enriched, tee, interactive)

        # Re-check tools that were installed
        if installed_count > 0:
            tee.write("\nRe-checking after installation...\n")
            still_missing = []
            for name, cfg in missing:
                if name in failed_list:
                    still_missing.append((name, cfg))
                    continue
                version = _get_version(cfg['cmd'], cfg.get('pattern'))
                if version:
                    tee.write(f"  {name} version {version}  ✓\n")
                else:
                    still_missing.append((name, cfg))
            missing = still_missing

    # ── Phase 3: check R packages ─────────────────────────────────────
    rscript_path = os.path.join(prefix, "scripts", "checkPackages.R")
    if os.path.exists(rscript_path):
        # The bundled R script handles its own installation
        try:
            subprocess.run(
                f"Rscript --vanilla {rscript_path}",
                shell=True, timeout=1800
            )
        except subprocess.TimeoutExpired:
            tee.write("Warning: R package check timed out\n")
        except Exception:
            tee.write("Warning: R package check failed\n")
    else:
        tee.write("Warning: scripts/checkPackages.R not found\n")

    # ── Phase 4: fail on remaining required missing ───────────────────
    if missing:
        required_missing = [name for name, cfg in missing if cfg.get('required', True)]
        if required_missing:
            tee.write(
                f"\nERROR: Required dependencies still missing: "
                f"{', '.join(required_missing)}\n"
            )
            tee.write("Install them manually or run with --auto-install\n")
            sys.exit(1)

    tee.write("Precheck completed!\n\n")
