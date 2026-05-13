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
  R           — DESeq2, DMRcaller, riboWaltz (degradome), NMF, Seurat,
               pheatmap, RNAmodR.RiboMethSeq
  git+make    — CLIPper (last resort)
  manual      — RiboTaper (GitHub: hsinyenwu/RiboTaper)
"""

import subprocess
import sys
import os
import shutil
import re
from pathlib import Path


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
    except subprocess.CalledProcessError:
        return False


# ═══════════════════════════════════════════════════════════════════════════
# 2. Dependency registry
# ═══════════════════════════════════════════════════════════════════════════
#
# Each entry keys on the binary name.  Fields:
#   pkg         — conda-forge / bioconda package name (may differ from binary)
#   channel     — conda channel
#   pip         — PyPI package name (for pip fallback)
#   r_pkg       — R package name (for BiocManager::install)
#   github      — GitHub URL for manual clone+install
#   mode_only   — list of modes that need this tool; None means always needed
#   install_msg — shown to user before installing

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
        'install_msg': 'MACS3 peak caller (ChIP/ATAC)',
        'mode_only': ['chip', 'atac'],
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

    # ── R packages (checked separately) ───────────────────────────────
    'R::DESeq2': {
        'r_pkg': 'DESeq2',
        'install_msg': 'Differential expression (R/Bioconductor)',
    },
    'R::DMRcaller': {
        'r_pkg': 'DMRcaller',
        'install_msg': 'Differential methylation (R/Bioconductor)',
    },
    'R::pheatmap': {
        'r_pkg': 'pheatmap',
        'install_msg': 'Heatmap plotting (R)',
    },
    'R::RNAmodR.RiboMethSeq': {
        'r_pkg': 'RNAmodR.RiboMethSeq',
        'install_msg': 'RiboMeth-seq analysis (R/Bioconductor)',
    },
    'R::riboWaltz': {
        'r_pkg': 'riboWaltz',
        'github': 'https://github.com/LabTranslationalArchitectomics/riboWaltz',
        'install_msg': 'Degradome CRI analysis (R/GitHub)',
        'mode_only': ['degradome'],
    },
    # ── CiPS uORF analysis R packages ────────────────────────────────
    'R::ORFik': {
        'r_pkg': 'ORFik',
        'install_msg': 'uORF detection & analysis (R/Bioconductor)',
        'mode_only': ['cips'],
    },
    'R::NMF': {
        'r_pkg': 'NMF',
        'github': 'https://github.com/renozao/NMF',
        'install_msg': 'Non-negative matrix factorization (R)',
    },
    'R::Seurat': {
        'r_pkg': 'Seurat',
        'github': 'https://github.com/satijalab/seurat',
        'install_msg': 'Single-cell analysis (R/GitHub)',
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# 3. Installer dispatch
# ═══════════════════════════════════════════════════════════════════════════

def install_tool(name, info, tee=None):
    """
    Install a single missing tool.
    Tries in order: mamba/conda → pip → R → GitHub clone.
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
        except subprocess.CalledProcessError:
            pass

    # 3. Try R / Bioconductor
    if info.get('r_pkg'):
        r_pkg = info['r_pkg']
        tee.write(f"  Trying R BiocManager::install('{r_pkg}')...\n")
        try:
            r_cmd = (
                f'Rscript -e \'if(!require("{r_pkg}", quietly=TRUE))'
                f'{{BiocManager::install("{r_pkg}", update=FALSE, ask=FALSE)}}\''
            )
            # Handle GitHub-hosted R packages
            if info.get('github') and 'riboWaltz' in r_pkg:
                r_cmd = (
                    f'Rscript -e \'if(!require("{r_pkg}", quietly=TRUE))'
                    f'{{devtools::install_github("LabTranslationalArchitectomics/riboWaltz", '
                    f'dependencies=TRUE)}}\''
                )
            elif info.get('github') and 'NMF' in r_pkg:
                r_cmd = (
                    f'Rscript -e \'if(!require("{r_pkg}", quietly=TRUE))'
                    f'{{devtools::install_github("renozao/NMF@devel", '
                    f'dependencies=TRUE)}}\''
                )
            elif info.get('github') and 'Seurat' in r_pkg:
                r_cmd = (
                    f'Rscript -e \'if(!require("{r_pkg}", quietly=TRUE))'
                    f'{{remotes::install_github("satijalab/seurat")}}\''
                )

            subprocess.run(r_cmd, shell=True, check=True, timeout=600)
            tee.write(f"  {name} installed via R\n")
            return True
        except subprocess.CalledProcessError:
            pass

    # 4. Try GitHub clone + manual install (CLIPper)
    github_url = info.get('github')
    if github_url and info.get('pkg') is None and not info.get('r_pkg'):
        tee.write(f"  Cloning {github_url}...\n")
        repo_name = github_url.rstrip('/').split('/')[-1].replace('.git', '')
        try:
            subprocess.run(
                ['git', 'clone', '--depth', '1', github_url, f'/tmp/{repo_name}'],
                check=True, timeout=120
            )
            subprocess.run(
                [sys.executable, 'setup.py', 'install'],
                cwd=f'/tmp/{repo_name}', check=True, timeout=120
            )
            subprocess.run(['rm', '-rf', f'/tmp/{repo_name}'], check=False)
            tee.write(f"  {name} installed from GitHub\n")
            return True
        except subprocess.CalledProcessError:
            pass

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

def install_r_packages(prefix, tee=None):
    """
    Run the bundled checkPackages.R to install all missing R packages.
    This is the most reliable way — it uses BiocManager + devtools.
    """
    if tee is None:
        tee = sys.stderr

    rscript_path = os.path.join(prefix, "scripts", "checkPackages.R")
    if not os.path.exists(rscript_path):
        tee.write("Warning: scripts/checkPackages.R not found\n")
        return False

    tee.write("\nInstalling R packages (this may take several minutes)...\n")
    try:
        subprocess.run(
            f"Rscript --vanilla {rscript_path}",
            shell=True, check=True, timeout=1800
        )
        tee.write("R packages installed successfully\n")
        return True
    except subprocess.CalledProcessError:
        tee.write("Some R packages could not be installed automatically\n")
        return False