"""
Automatic dependency installer for pRNASeqTools.
Detects missing external tools and installs them via the best available
package manager (mamba > conda > pip > manual).

Installation sources per tool:
  bioconda   — cutadapt, samtools, bowtie, bowtie2, star, subread,
               bedtools, sra-tools, gffread, deeptools, genrich,
               umi_tools, ucsc-bedgraphtobigwig, bismark, shortstack,
               stringtie, gffcompare, rsem, htslib, numpy, scipy
  conda-forge — r-base, r-essentials, bioconductor-*
  pip         — clipper (CLIPper), numpy, scipy (fallback)
  R           — batch install via checkPackages.R (BiocManager/devtools)
  git+make    — CLIPper (last resort)
  manual      — RiboTaper (GitHub: hsinyenwu/RiboTaper)
"""

import subprocess
import sys
import os
import shutil


# ═══════════════════════════════════════════════════════════════════════════
# 1. Package manager detection
# ═══════════════════════════════════════════════════════════════════════════

def _which(cmd):
    """Return path to executable or None."""
    return shutil.which(cmd)


def _detect_pm():
    """Return the best available package manager: 'mamba', 'conda', or None."""
    if _which('mamba'):
        return 'mamba'
    if _which('conda'):
        return 'conda'
    return None


def _pm_install(pm, packages, channel=None):
    """Install packages via mamba or conda. Returns True on success."""
    if isinstance(packages, str):
        packages = [packages]
    cmd = [pm, 'install', '-y', '-q']
    if channel:
        cmd += ['-c', channel]
    cmd += packages
    try:
        subprocess.run(cmd, check=True, timeout=600)
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False


# ═══════════════════════════════════════════════════════════════════════════
# 2. Dependency registry
# ═══════════════════════════════════════════════════════════════════════════
#
# Each entry keys on the binary name.  Fields:
#   pkg            — conda-forge / bioconda package name (may differ from binary)
#   channel        — conda channel
#   pip            — PyPI package name (for pip fallback)
#   github         — GitHub URL for manual clone+install
#   mode_only      — list of modes that need this tool; None means always needed
#   install_msg    — shown to user before installing

DEPENDENCY_REGISTRY = {
    # ── Core tools (always required) ───────────────────────────────────
    'cutadapt': {
        'pkg': 'cutadapt', 'channel': 'bioconda',
        'install_msg': 'Adapter trimming tool',
    },
    'samtools': {
        'pkg': 'samtools', 'channel': 'bioconda',
        'install_msg': 'SAM/BAM processing',
    },
    'bowtie': {
        'pkg': 'bowtie', 'channel': 'bioconda',
        'install_msg': 'Short-read aligner (sRNA, DNA)',
    },
    'bowtie2': {
        'pkg': 'bowtie2', 'channel': 'bioconda',
        'install_msg': 'Short-read aligner (ChIP, ATAC)',
    },
    'featureCounts': {
        'pkg': 'subread', 'channel': 'bioconda',
        'install_msg': 'Read counting (ships inside subread)',
    },
    'bedtools': {
        'pkg': 'bedtools', 'channel': 'bioconda',
        'install_msg': 'Genome interval arithmetic',
    },
    'R': {
        'pkg': 'r-base', 'channel': 'conda-forge',
        'install_msg': 'R statistical computing',
    },
    'STAR': {
        'pkg': 'star', 'channel': 'bioconda',
        'install_msg': 'RNA-seq aligner',
    },
    'gffread': {
        'pkg': 'gffread', 'channel': 'bioconda',
        'install_msg': 'GFF/GTF conversion',
    },
    'deeptools': {
        'pkg': 'deeptools', 'channel': 'bioconda',
        'install_msg': 'Coverage track tools (bamCoverage)',
    },
    'bedGraphToBigWig': {
        'pkg': 'ucsc-bedgraphtobigwig', 'channel': 'bioconda',
        'install_msg': 'UCSC bigWig conversion',
    },
    'bgzip': {
        'pkg': 'htslib', 'channel': 'bioconda',
        'install_msg': 'Tabix/BGZip compression (WGBS)',
        'mode_only': ['wgbs'],
    },
    'tabix': {
        'pkg': 'htslib', 'channel': 'bioconda',
        'install_msg': 'Tabix indexing (WGBS)',
        'mode_only': ['wgbs'],
    },

    # ── Mode-specific tools ────────────────────────────────────────────
    'ShortStack': {
        'pkg': 'shortstack', 'channel': 'bioconda',
        'install_msg': 'Small RNA alignment & annotation',
        'mode_only': ['srna', 'phasi', 'tt', 'risi'],
    },
    'fasterq-dump': {
        'pkg': 'sra-tools', 'channel': 'bioconda',
        'install_msg': 'SRA data download',
    },
    'Genrich': {
        'pkg': 'genrich', 'channel': 'bioconda',
        'install_msg': 'ChIP/ATAC peak caller',
        'mode_only': ['chip', 'atac'],
    },
    'umi_tools': {
        'pkg': 'umi_tools', 'channel': 'bioconda',
        'install_msg': 'UMI extraction (single-cell sRNA)',
        'mode_only': ['srna'],
    },
    'Bismark': {
        'pkg': 'bismark', 'channel': 'bioconda',
        'install_msg': 'Bisulfite alignment (WGBS)',
        'mode_only': ['wgbs'],
    },
    'clipper': {
        'pkg': None,  # not in conda
        'pip': 'clipper',
        'github': 'https://github.com/YeoLab/clipper.git',
        'install_msg': 'CLIP-seq peak caller',
        'mode_only': ['clip', 'ts'],
    },

    # ── Peak callers ────────────────────────────────────────────────
    'macs3': {
        'pkg': 'macs3', 'channel': 'bioconda',
        'install_msg': 'MACS3 peak caller (ChIP/ATAC/TF)',
        'mode_only': ['chip', 'atac', 'tf'],
    },

    # ── RIBO Taper / Ribo-seq tools ──────────────────────────────────
    'stringtie': {
        'pkg': 'stringtie', 'channel': 'bioconda',
        'install_msg': 'Transcriptome assembly',
        'mode_only': ['ribo'],
    },
    'gffcompare': {
        'pkg': 'gffcompare', 'channel': 'bioconda',
        'install_msg': 'GTF comparison & class-code annotation',
        'mode_only': ['ribo'],
    },
    'rsem-prepare-reference': {
        'pkg': 'rsem', 'channel': 'bioconda',
        'install_msg': 'RNA-seq expression quantification (RSEM)',
        'mode_only': ['ribo'],
    },

    # ── Python libs (sPARTA.py) ────────────────────────────────────────
    'numpy': {
        'pkg': 'numpy', 'channel': 'conda-forge',
        'pip': 'numpy',
        'install_msg': 'Numerical Python (sPARTA)',
    },
    'scipy': {
        'pkg': 'scipy', 'channel': 'conda-forge',
        'pip': 'scipy',
        'install_msg': 'Scientific Python (sPARTA)',
    },

}


# ═══════════════════════════════════════════════════════════════════════════
# 2b. R package registry — single source of truth for R package installs
# ═══════════════════════════════════════════════════════════════════════════
#
# Each entry keys on the R package name.  Fields:
#   pkg            — R package name (required)
#   source         — 'bioconductor', 'github', or 'cran'
#   github_repo    — GitHub repo path (source=github)
#   github_ref     — branch/tag (optional, e.g. 'devel')
#   extra          — additional packages to install alongside
#   mode_only      — list of modes that need this package; None = all modes
#   install_msg    — shown to user

R_PACKAGE_REGISTRY = {
    'dplyr': {
        'pkg': 'dplyr',
        'source': 'cran',
        'install_msg': 'Data manipulation (R/CRAN)',
    },
    'DESeq2': {
        'pkg': 'DESeq2',
        'source': 'bioconductor',
        'install_msg': 'Differential expression (R/Bioconductor)',
    },
    'DMRcaller': {
        'pkg': 'DMRcaller',
        'source': 'bioconductor',
        'install_msg': 'Differential methylation (R/Bioconductor)',
    },
    'pheatmap': {
        'pkg': 'pheatmap',
        'source': 'bioconductor',
        'install_msg': 'Heatmap plotting (R)',
    },
    'RNAmodR.RiboMethSeq': {
        'pkg': 'RNAmodR.RiboMethSeq',
        'source': 'bioconductor',
        'install_msg': 'RiboMeth-seq analysis (R/Bioconductor)',
    },
    'riboWaltz': {
        'pkg': 'riboWaltz',
        'source': 'github',
        'github_repo': 'LabTranslationalArchitectomics/riboWaltz',
        'install_msg': 'Degradome CRI analysis (R/GitHub)',
        'mode_only': ['degradome'],
    },
    'ORFik': {
        'pkg': 'ORFik',
        'source': 'bioconductor',
        'install_msg': 'uORF detection & analysis (R/Bioconductor)',
        'mode_only': ['cips'],
    },
    'NMF': {
        'pkg': 'NMF',
        'source': 'github',
        'github_repo': 'renozao/NMF',
        'github_ref': 'devel',
        'install_msg': 'Non-negative matrix factorization (R)',
    },
    'Seurat': {
        'pkg': 'Seurat',
        'source': 'github',
        'github_repo': 'satijalab/seurat',
        'extra': ['uwot'],
        'install_msg': 'Single-cell analysis (R/GitHub)',
    },
}


def get_r_packages_for_mode(mode=None):
    """
    Return list of (pkg_name, info) tuples for the given mode.
    If mode is None, returns ALL registered R packages.
    """
    pkgs = []
    for name, info in R_PACKAGE_REGISTRY.items():
        mode_only = info.get('mode_only')
        if mode_only and mode not in mode_only:
            continue
        pkgs.append((name, info))
    return pkgs


# ═══════════════════════════════════════════════════════════════════════════
# 3. Installer dispatch
# ═══════════════════════════════════════════════════════════════════════════

def install_tool(name, info, tee=None):
    """
    Install a single missing tool.
    Tries in order: mamba/conda → pip → GitHub clone.
    Returns True if installation succeeded.
    """
    if tee is None:
        tee = sys.stderr

    pm = _detect_pm()
    msg = info.get('install_msg', name)
    tee.write(f"\n  Auto-installing {name} ({msg})...\n")

    # 1. Try conda/mamba
    if pm and info.get('pkg'):
        channel = info.get('channel')
        tee.write(f"  Trying {pm} install {info['pkg']}...\n")
        if _pm_install(pm, info['pkg'], channel):
            tee.write(f"  {name} installed via {pm}\n")
            return True

    # 2. Try pip
    if info.get('pip'):
        tee.write(f"  Trying pip install {info['pip']}...\n")
        try:
            subprocess.run(
                [sys.executable, '-m', 'pip', 'install', '-q', info['pip']],
                check=True, timeout=300
            )
            tee.write(f"  {name} installed via pip\n")
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            pass

    # 3. Try GitHub clone + manual install (CLIPper)
    github_url = info.get('github')
    if github_url and info.get('pkg') is None:
        tee.write(f"  Cloning {github_url}...\n")
        repo_name = github_url.rstrip('/').split('/')[-1].replace('.git', '')
        clone_dir = os.path.join('/tmp', repo_name)
        if os.path.exists(clone_dir):
            shutil.rmtree(clone_dir, ignore_errors=True)
        try:
            subprocess.run(
                ['git', 'clone', '--depth', '1', github_url, clone_dir],
                check=True, timeout=120
            )
            subprocess.run(
                [sys.executable, 'setup.py', 'install'],
                cwd=clone_dir, check=True, timeout=120
            )
            tee.write(f"  {name} installed from GitHub\n")
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            pass
        finally:
            shutil.rmtree(clone_dir, ignore_errors=True)

    tee.write(f"  Failed to auto-install {name}\n")
    return False


def install_missing(missing_list, tee=None, interactive=True):
    """
    Install a list of missing tool names.
    
    Args:
        missing_list: list of (name, info) tuples from DEPENDENCY_REGISTRY
        tee: output handle
        interactive: if True, prompt user before each install
    
    Returns:
        (installed_count, failed_list)
    """
    if tee is None:
        tee = sys.stderr

    if not missing_list:
        return 0, []

    if interactive:
        names = [name for name, _ in missing_list]
        tee.write(f"\nMissing dependencies: {', '.join(names)}\n")
        tee.write("Attempt automatic installation? [Y/n] ")
        tee.flush()
        try:
            answer = input().strip().lower()
            if answer and answer not in ('y', 'yes', ''):
                tee.write("Skipping auto-install. Please install manually.\n")
                return 0, [name for name, _ in missing_list]
        except EOFError:
            pass  # non-interactive, proceed

    installed = 0
    failed = []

    for name, info in missing_list:
        if install_tool(name, info, tee):
            installed += 1
        else:
            failed.append(name)

    if installed:
        tee.write(f"\nAuto-installed {installed} package(s)\n")
    if failed:
        tee.write(f"Could not auto-install: {', '.join(failed)}\n")

    return installed, failed


# ═══════════════════════════════════════════════════════════════════════════
# 4. R package batch install (with retry)
# ═══════════════════════════════════════════════════════════════════════════

def install_r_packages(prefix, mode=None, tee=None):
    """
    Run the bundled checkPackages.R to install missing R packages.
    Filters by analysis mode for efficient, targeted installation.

    Args:
        prefix: project root (contains scripts/checkPackages.R)
        mode:   current analysis mode (e.g. 'srna', 'cips'); None = install all
        tee:    output handle

    Returns:
        True if all packages installed successfully, False otherwise.
    """
    if tee is None:
        tee = sys.stderr

    rscript_path = os.path.join(prefix, "scripts", "checkPackages.R")
    if not os.path.exists(rscript_path):
        tee.write("Warning: scripts/checkPackages.R not found\n")
        return False

    pkgs = get_r_packages_for_mode(mode)
    if not pkgs:
        tee.write("\nNo R packages required for this mode.\n")
        return True

    pkg_names = ','.join(name for name, _ in pkgs)
    tee.write("\nInstalling R packages (this may take several minutes)...\n")
    try:
        subprocess.run(
            f"Rscript --vanilla {rscript_path} --packages {pkg_names}",
            shell=True, check=True, timeout=1800
        )
        tee.write("R packages installed successfully\n")
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        tee.write("Some R packages could not be installed automatically\n")
        return False