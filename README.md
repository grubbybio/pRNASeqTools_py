# pRNASeqTools

**Integrated High-throughput Sequencing Data Analysis for Plant** — Python3 rewrite of the original Perl pipeline.

Author: Dr. Chenjiang You  
Python rewrite: DeepSeek  

---

## Overview

pRNASeqTools is a comprehensive NGS data analysis pipeline for plant genomics. It provides **16 analysis modes** covering small RNA, mRNA, epigenomics, ribosome profiling, and more — all accessible through a single command-line interface.

### Analysis Modes

| Mode | Description | Key Tools |
|------|-------------|-----------|
| `srna` | Small RNA-seq | ShortStack, bowtie |
| `mrna` | mRNA-seq (DE analysis) | STAR, featureCounts, DESeq2 |
| `degradome` | Degradome-seq (PARE/GMUCT) | STAR, sPARTA, riboWaltz |
| `phasi` | phasiRNA analysis | ShortStack |
| `tt` | Truncation/tailing | ShortStack |
| `ribo` | **Ribo-seq (RIBO Taper)** | Bowtie2, STAR, StringTie, RSEM, RiboTaper |
| `cips` | **CiPS uORF analysis** | ORFik, GenomicFeatures (R) |
| `chip` | ChIP-seq | bowtie2, Genrich, deepTools |
| `atac` | ATAC-seq | bowtie2, Genrich, deepTools |
| `wgbs` | Whole-genome bisulfite | Bismark, DMRcaller |
| `clip` | CLIP-seq | STAR, CLIPper |
| `ts` | TS-CLIP-seq | STAR, CLIPper |
| `ribometh` | RiboMeth-seq | STAR, RNAmodR |
| `risi` | risiRNA analysis | ShortStack |
| `tf` | Two-factor DE | DESeq2 |

---

## Quick Start

### 1. Install dependencies

```bash
# Create conda environment
mamba env create -f environment.yaml
mamba activate prnaseqtools

# Install R packages from Bioconductor/GitHub
Rscript scripts/checkPackages.R
```

### 2. Prepare reference files

Place reference genome files in `reference/`:

```
reference/
├── ath_genes.gff          # Gene annotation (GFF3)
├── ath_chr_all.fasta      # Genome sequence
├── ath_chr_all.fasta.fai  # FASTA index
├── ath_miRNA_miRNA_star.gff  # miRNA annotation
└── ...
```

Supported genomes: `ath` (Arabidopsis), `osa` (Rice), `b73` (Maize), `gma` (Soybean), `smo`, `bra`, `w22`

### 3. Run analysis

```bash
# Small RNA-seq
python pRNASeqTools_run.py srna -c "WT=data/WT.fq" -p "mut=data/mut.fq"

# mRNA-seq differential expression
python pRNASeqTools_run.py mrna -c "WT=SRR111111" -t "mut=SRR222222"

# Ribo-seq (RIBO Taper pipeline)
python pRNASeqTools_run.py ribo \
  --rna-control "WT=SRR111111" \
  --ribo-control "Ribo=SRR333333" \
  --contam "rRNA.fasta,tRNA.fasta,snRNA.fasta" \
  --ribotaper ~/RiboTaper_v1.3

# ChIP-seq
python pRNASeqTools_run.py chip -c "IP=data/ip.bam" -t "Input=data/input.bam"
```

---

## RIBO Taper Pipeline (`ribo` mode)

The `ribo` mode implements a full RIBO Taper workflow for translated ORF detection:

```
Step 1 — Bowtie2 contamination removal (rRNA, tRNA, snRNA, snoRNA)
Step 2 — Preprocess Ribo-seq reads
Step 3 — RNA-seq STAR 2-pass + StringTie transcriptome assembly
Step 4 — gffcompare → novel transcript filtering + gene_biotype
Step 5 — RSEM quantification → expressed isoform filtering (TPM)
Step 6 — STAR re-mapping with expressed annotation
Step 7 — RIBO Taper annotation files
Step 8 — Merge BAMs + Ribotaper.sh ORF detection
Step 9 — Output summary
```

**Additional dependency:** [RiboTaper](https://github.com/hsinyenwu/RiboTaper) (manual install)

```bash
git clone https://github.com/hsinyenwu/RiboTaper.git ~/RiboTaper_v1.3
```

### RiboTaper CLI options

| Option | Default | Description |
|--------|---------|-------------|
| `--rna-control` | *(required)* | RNA-seq control: `name=file1+file2...` |
| `--ribo-control` | *(required)* | Ribo-seq control: `name=file1+file2...` |
| `--contam` | *(required)* | Contamination fasta for Bowtie2 index |
| `--ribo-len` | `24,25,26,27,28` | Ribo-seq read lengths |
| `--cutoffs` | `8,9,10,11,12` | RIBO Taper cutoffs |
| `--tpm-threshold` | `0` | Mean TPM threshold for isoform filtering |
| `--ribotaper` | *(auto-detect)* | Path to RIBO Taper installation |

---

## CiPS uORF Analysis (`cips` mode)

Downstream analysis after RIBO Taper: identifies **translated upstream ORFs (uORFs)** in 5' UTRs using Ribo-seq P-site periodicity.

```
Input: GTF + FASTA + P-site counts
  → Extract 5' UTRs
  → Find all uORFs (ORFik)
  → Compute frame-specific P-site counts (parallel)
  → Filter translated uORFs (minimum: 1aa, longer: ≥2aa)
  → Deduplicate overlapping uORFs
Output: 4 Excel files
```

### CiPS CLI options

| Option | Default | Description |
|--------|---------|-------------|
| `--gtf` | `<genome>_expressed.gtf` | Expressed GTF from ribo pipeline |
| `--fasta` | `reference/<genome>_chr_all.fasta` | Reference genome |
| `--psite` | *(required)* | P-site count file |
| `--min-inframe-counts` | `10` | Min in-frame Ribo-seq counts |
| `--min-inframe-perc` | `50` | Min in-frame percentage |
| `--min-psite-perc` | `30` | Min P-site % (longer uORFs) |
| `--gene-desc` | *(none)* | Gene description Excel |

---

## Common Options (all modes)

| Option | Description |
|--------|-------------|
| `--outdir`, `-o` | Output directory (default: `./out`) |
| `--genome`, `-g` | Genome: `ath`, `osa`, `b73`, `gma`, `smo`, `bra`, `w22` |
| `--thread`, `-t` | Number of threads (default: 4) |
| `--adaptor`, `-a` | 3' adaptor sequence |
| `--control`, `-c` | Control samples: `name=file1+file2...` |
| `--treatment`, `-p` | Treatment samples |

---

## Dependencies

### Core tools (via conda)

**Aligners:** STAR, bowtie, bowtie2, Bismark, ShortStack  
**Processing:** cutadapt, samtools, bedtools, gffread, deepTools  
**Counting:** featureCounts (subread), RSEM, StringTie  
**Peak calling:** Genrich  
**SRA:** sra-tools  
**R:** r-base, DESeq2, DMRcaller, pheatmap, dplyr, devtools

### Mode-specific R packages

| Package | Mode | Source |
|---------|------|--------|
| riboWaltz | degradome | GitHub |
| RNAmodR.RiboMethSeq | ribometh | Bioconductor |
| NMF | — | GitHub |
| Seurat | srna (sc) | GitHub |
| ORFik | cips | Bioconductor |

### Manual install

- **CLIPper** — `clip` / `ts` modes: `pip install clipper`
- **RiboTaper** — `ribo` mode: `git clone https://github.com/hsinyenwu/RiboTaper.git`
- **Reference files** — contact the author

---

## Input Format

Samples are specified as `name=source` pairs. Sources can be:

- **SRA accessions:** `WT=SRR123456` (auto-download via fasterq-dump)
- **Local files:** `WT=data/sample.fq` or `WT=data/sample.fq.gz`
- **Multiple replicates:** `WT=rep1.fq+rep2.fq+rep3.fq`
- **Paired-end:** handled automatically for SRA accessions

---

## Project Structure

```
pRNASeqTools_py/
├── pRNASeqTools_run.py        # Entry point
├── environment.yaml           # Conda environment
├── prnaseqtools/
│   ├── cli.py                 # CLI definition (argparse)
│   ├── auto_install.py        # Dependency auto-installer
│   ├── validate_options.py    # Input validation
│   ├── input_parser.py        # Sample specification parser
│   ├── functions.py           # Utilities (download, unzip, revcomp)
│   ├── reference.py           # Genome reference handling
│   ├── precheck.py            # Dependency checker
│   ├── logging_setup.py       # Logging
│   └── modes/                 # Analysis mode modules
│       ├── srna.py
│       ├── mrna.py
│       ├── degradome.py
│       ├── phasi.py
│       ├── tt.py
│       ├── ribo.py            # RIBO Taper pipeline
│       ├── cips.py            # CiPS uORF analysis
│       ├── chip.py
│       ├── atac.py
│       ├── wgbs.py
│       ├── clip.py
│       ├── ts.py
│       ├── ribometh.py
│       ├── risi.py
│       └── tf.py
├── scripts/                   # R analysis scripts
│   ├── checkPackages.R
│   ├── ribo.R
│   ├── ribotaper_filter_gtf.R
│   ├── ribotaper_filter_rsem.R
│   ├── cips_uORF.R
│   └── ...
└── reference/                 # Genome files (not versioned)
```

---

## License

Contact the author for licensing information.
