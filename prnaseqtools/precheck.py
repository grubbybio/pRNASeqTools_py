"""
Dependency checker for pRNASeqTools — with auto-install support.
Checks that all required external tools are installed, and can
automatically install missing ones via mamba/conda/pip/R/git.
"""

import subprocess
import sys
import os
import re
import json
import hashlib
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

_CACHE_DIR = Path(__file__).resolve().parent.parent / '.cache'
_CACHE_FILE = _CACHE_DIR / 'precheck_cache.json'


def _get_cache_key(mode):
    """Generate a cache key based on environment and mode."""
    env_info = []
    
    env_info.append(os.environ.get('PATH', ''))
    env_info.append(os.environ.get('CONDA_PREFIX', ''))
    env_info.append(os.environ.get('VIRTUAL_ENV', ''))
    env_info.append(sys.executable)
    
    key = hashlib.md5(('|'.join(env_info) + f'|{mode}').encode()).hexdigest()
    return key


def _load_cache():
    """Load cached dependency check results."""
    if not _CACHE_FILE.exists():
        return {}
    try:
        with open(_CACHE_FILE, 'r') as f:
            return json.load(f)
    except Exception:
        return {}


def _save_cache(cache):
    """Save dependency check results to cache."""
    _CACHE_DIR.mkdir(exist_ok=True)
    with open(_CACHE_FILE, 'w') as f:
        json.dump(cache, f)


def _get_conda_packages():
    """Get list of installed packages in current conda environment."""
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    if not conda_prefix:
        return ''
    
    try:
        result = subprocess.run(
            ['conda', 'list', '--prefix', conda_prefix, '--json'],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode == 0:
            pkgs = json.loads(result.stdout)
            pkg_list = [f"{p['name']}={p['version']}" for p in pkgs]
            pkg_list.sort()
            return '\n'.join(pkg_list)
    except Exception:
        pass
    
    return ''


def _get_env_fingerprint():
    """Generate a fingerprint of the current environment."""
    parts = []
    parts.append(os.environ.get('PATH', ''))
    parts.append(os.environ.get('CONDA_PREFIX', ''))
    parts.append(os.environ.get('VIRTUAL_ENV', ''))
    parts.append(sys.executable)
    parts.append(_get_conda_packages())
    
    return hashlib.md5('|'.join(parts).encode()).hexdigest()


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
        'pattern': r'ShortStack (\d+\.\d+\.\d+)',
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
    'macs3': {
        'cmd': 'macs3 --version',
        'pattern': r'macs3 (\d+\.\d+\.\d+)',
        'required': False,
        'msg': 'Peak caller (ChIP/ATAC/TF)',
        'mode_only': ['chip', 'atac', 'tf'],
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
    'bgzip': {
        'cmd': 'bgzip --version 2>&1',
        'pattern': r'(\d+\.\d+)',
        'required': False,
        'msg': 'BGZip compression (WGBS)',
        'mode_only': ['wgbs'],
    },
    'tabix': {
        'cmd': 'tabix --version 2>&1',
        'pattern': r'(\d+\.\d+)',
        'required': False,
        'msg': 'Tabix indexing (WGBS)',
        'mode_only': ['wgbs'],
    },
    'clipper': {
        'cmd': 'clipper -h 2>&1',
        'pattern': r'^Usage',
        'required': False,
        'msg': 'CLIP-seq peak caller',
        'mode_only': ['clip', 'ts'],
    },
    'stringtie': {
        'cmd': 'stringtie --version',
        'pattern': r'version (\d+\.\d+\.\d+)',
        'required': False,
        'msg': 'Transcriptome assembly (Ribo-seq)',
        'mode_only': ['ribo'],
    },
    'gffcompare': {
        'cmd': 'gffcompare --version',
        'pattern': r'v(\d+\.\d+\.\d+)',
        'required': False,
        'msg': 'GTF comparison (Ribo-seq)',
        'mode_only': ['ribo'],
    },
    'rsem-prepare-reference': {
        'cmd': 'rsem-prepare-reference -h 2>&1',
        'pattern': r'v(\d+\.\d+\.\d+)',
        'required': False,
        'msg': 'RSEM reference preparation (Ribo-seq)',
        'mode_only': ['ribo'],
    },
    'rsem-calculate-expression': {
        'cmd': 'rsem-calculate-expression -h 2>&1',
        'pattern': r'v(\d+\.\d+\.\d+)',
        'required': False,
        'msg': 'RSEM expression quantification (Ribo-seq)',
        'mode_only': ['ribo'],
    },
    'sPARTA.py': {
        'cmd': 'python3 -c "import numpy, scipy; print(\'sPARTA deps OK\')" 2>&1',
        'pattern': r'sPARTA deps (OK)',
        'required': False,
        'msg': 'sPARTA degradome-seq peak caller (needs numpy, scipy)',
        'mode_only': ['degradome'],
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

    cache_key = _get_cache_key(mode)
    cache = _load_cache()
    
    if cache_key in cache:
        cached_env = cache[cache_key].get('env_fingerprint')
        current_env = _get_env_fingerprint()
        if cached_env == current_env:
            tee.write("\nSkipping dependency check (cache valid)...\n")
            tee.write("Precheck completed!\n\n")
            return

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
    from prnaseqtools.auto_install import install_r_packages
    try:
        install_r_packages(prefix, mode=mode, tee=tee)
    except Exception:
        tee.write("Warning: R package check failed\n")

    # ── Phase 3.5: ensure mapping indices exist ───────────────────────
    from prnaseqtools.reference import check_and_build_indices
    try:
        check_and_build_indices(prefix, mode=mode, tee=tee)
    except Exception:
        tee.write("Warning: Index check failed\n")

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

    cache[cache_key] = {
        'env_fingerprint': _get_env_fingerprint(),
        'mode': mode,
    }
    _save_cache(cache)

    tee.write("Precheck completed!\n\n")
