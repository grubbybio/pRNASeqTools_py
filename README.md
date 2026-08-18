# pRNASeqTools

**Integrated High-throughput Sequencing Data Analysis for Plant** — Python3 rewrite of the original Perl pipeline.

Author: Dr. Chenjiang You  
Version: 1.0.0

---

## Overview

pRNASeqTools is a comprehensive NGS data analysis pipeline for plant genomics. It provides **15 analysis modes** covering small RNA, mRNA, epigenomics, ribosome profiling, and more — all accessible through a single command-line interface.

### Analysis Modes

| Mode | Description | Key Tools |
|------|-------------|-----------|
| `srna` | Small RNA-seq — bulk & single-cell | ShortStack, bowtie, cutadapt, umi_tools |
| `mrna` | mRNA-seq — DE analysis | STAR, featureCounts, DESeq2 |
| `degradome` | Degradome-seq (PARE/GMUCT) — miRNA target cleavage | STAR (dual-alignment), sPARTA, riboWaltz |
| `phasi` | phasiRNA analysis — phased secondary siRNA detection | ShortStack, bowtie |
| `tt` | miRNA truncation/tailing analysis | bowtie (iterative), ShortStack |
| `ribo` | Ribo-seq — translated ORF detection (RIBO Taper pipeline) | Bowtie2, STAR, StringTie, RSEM, RiboTaper |
| `cips` | CiPS uORF analysis — translated upstream ORF detection | ORFik, GenomicFeatures (R) |
| `chip` | ChIP-seq — peak calling & differential peaks | bowtie2, Genrich / MACS3, deepTools |
| `atac` | ATAC-seq — open chromatin analysis | bowtie2, Genrich / MACS3, deepTools |
| `wgbs` | Whole-genome bisulfite — differential methylation | Bismark, DMRcaller |
| `clip` | CLIP-seq — protein-RNA interaction | STAR, CLIPper |
| `ts` | TS-CLIP-seq — target-specific CLIP | STAR, CLIPper |
| `ribometh` | RiboMeth-seq — 2'-O-methylation analysis | STAR, RNAmodR.RiboMethSeq |
| `risi` | risiRNA analysis | ShortStack, bowtie |
| `tf` | Two-factor DE — multi-condition differential analysis | DESeq2 |

---

## Quick Start

### 1. Install dependencies

```bash
# Create conda environment
conda env create -f environment.yaml
conda activate prnaseqtools

# Install R packages from Bioconductor/GitHub
Rscript scripts/checkPackages.R
```

The pipeline also supports **auto-install** of missing dependencies (enabled by default). Add `--no-auto-install` to disable.

### 2. Prepare reference files

Place reference genome files in `reference/`:

```
reference/
├── ath_genes.gff              # Gene annotation (GFF3)
├── ath_chr_all.fasta          # Genome sequence
├── ath_chr_all.fasta.fai      # FASTA index
├── ath_miRNA_miRNA_star.gff   # miRNA annotation
└── ...
```

Supported genomes: `ath` (Arabidopsis), `osa` (Rice), `b73` (Maize), `gma` (Soybean), `smo`, `bra`, `w22`

### 3. Run analysis

```bash
# Small RNA-seq (bulk mode)
python pRNASeqTools_run.py srna -c "WT=data/WT.fq" -p "mut=data/mut.fq"

# Small RNA-seq (single-cell mode)
python pRNASeqTools_run.py srna --mode_srna sc -c "WT=data/WT.fq"

# mRNA-seq differential expression
python pRNASeqTools_run.py mrna -c "WT=SRR111111" -t "mut=SRR222222"

# mRNA-seq — start from BAM files
python pRNASeqTools_run.py mrna -c "WT=data/WT.bam" -t "mut=data/mut.bam" \
  --mode_mrna bam --seqstrategy paired

# Degradome-seq — miRNA target cleavage analysis
python pRNASeqTools_run.py degradome -c "WT=data/degradome.fq" \
  --adaptor TGGAATTCTCGGG

# Ribo-seq (RIBO Taper pipeline)
python pRNASeqTools_run.py ribo \
  --rna-control "WT=SRR111111" \
  --ribo-control "Ribo=SRR333333" \
  --contam "rRNA.fasta,tRNA.fasta,snRNA.fasta"

# Ribo-seq — paired-end data
python pRNASeqTools_run.py ribo \
  --rna-control "WT=WT_r1.fq,WT_r2.fq" \
  --ribo-control "Ribo=Ribo_r1.fq,Ribo_r2.fq" \
  --contam "rRNA.fasta,tRNA.fasta"

# Ribo-seq — paired-end with multiple replicates
python pRNASeqTools_run.py ribo \
  --rna-control "WT=WT_r1_1.fq,WT_r2_1.fq+WT_r1_2.fq,WT_r2_2.fq" \
  --rna-treatment "mut=mut_r1_1.fq,mut_r2_1.fq+mut_r1_2.fq,mut_r2_2.fq" \
  --ribo-control "Ribo=Ribo_r1.fq,Ribo_r2.fq" \
  --ribo-treatment "Rmut=Rmut_r1.fq,Rmut_r2.fq" \
  --contam "rRNA.fasta,tRNA.fasta"

# Ribo-seq — multiple treatment groups
python pRNASeqTools_run.py ribo \
  --rna-control "WT=SRR111111" \
  --rna-treatment "mutA=SRR222222" --rna-treatment "mutB=SRR333333" \
  --ribo-control "Ribo=SRR444444" \
  --ribo-treatment "RmutA=SRR555555" --ribo-treatment "RmutB=SRR666666" \
  --contam "rRNA.fasta,tRNA.fasta"

# ChIP-seq (Genrich, default)
python pRNASeqTools_run.py chip --treatment "IP=data/ip.bam" --control "Input=data/input.bam"

# ChIP-seq (MACS3 + bdgdiff — two-group differential peaks)
python pRNASeqTools_run.py tf -c "WT=input,3,IP,3" -p "KO=input,3,IP,3" \
  --mode_tf chip --genome-size 1.35e8

# ATAC-seq (MACS3)
python pRNASeqTools_run.py atac --peak-caller macs3 --genome-size 1.35e8 \
  --treatment "ATAC=data/atac.bam"

# WGBS — differential methylation
python pRNASeqTools_run.py wgbs -c "WT=data/wt.fq" -p "mut=data/mut.fq"

# Two-factor DE (ChIP mode)
python pRNASeqTools_run.py tf -c "WT=input,3,IP,3" -p "KO=input,3,IP,3" \
  --mode_tf chip --genome-size 1.35e8
```

---

## Mode Details

### `srna` — Small RNA-seq

Full pipeline for small RNA-seq, supporting both bulk and single-cell (UMI-based) modes.

**Pipeline steps:**
1. SRA download (if needed) → cutadapt adapter trimming (retain 18–42 nt)
2. UMI extraction & deduplication (sc mode via `umi_tools`)
3. Optional: mask filtering, spike-in quantification (bowtie)
4. rRNA/SSU/U6 filtering via bowtie → normalization factors
5. ShortStack alignment to genome
6. BAM → BED conversion, length stratification (18–26 nt), length distribution
7. Multi-dimensional counting: bin (sliding window), gene, TE, promoter, miRNA
8. RPM normalization + bedGraph / bigWig generation
9. DESeq2 differential expression analysis (5 modes: DSR, DEM, DSG, DST, DSP)

| Option | Default | Description |
|--------|---------|-------------|
| `--mode_srna` | `bulk` | Run mode: `bulk` or `sc` |
| `--pattern` | `NNNNNNNNCA` | UMI pattern (sc mode) |
| `--mmap` | `u` | ShortStack multimap strategy |
| `--norm` | `rRNA,total` | Normalization methods (comma-separated) |
| `--binsize` | `100` | Window size for bin counts |
| `--promoter` | `1000` | Promoter region length |
| `--foldchange` | `1.5` | Fold-change cutoff |
| `--pvalue` | `0.01` | P-value cutoff |
| `--mask` | — | Mask FASTA for filtering |
| `--spike-in` | — | Spike-in FASTA for quantification |
| `--no-mapping` | — | Skip alignment, statistics only |
| `--mapping-only` | — | Alignment only, skip statistics |

**Adaptor aliases:** `truseq` / `illumina` / `srna` → `TGGAATTCTCGGG`, `neb` → `AGATCGGAAGAGC`, `nextera` → `CTGTCTCTTATAC`

---

### `mrna` — mRNA-seq

STAR-based mRNA-seq with featureCounts quantification and DESeq2 differential expression.

**Pipeline steps:**
1. STAR genome index + gffread GFF→GTF conversion
2. Per-sample: SRA download → cutadapt trimming → optional mask filtering → STAR alignment
3. samtools index + bamCoverage (CPM bigWig)
4. featureCounts gene-level quantification
5. DESeq2 differential expression (DEG.R)

| Option | Default | Description |
|--------|---------|-------------|
| `--mode_mrna` | `whole` | `whole`=full pipeline, `mapping-only`=alignment+count, `bam`=BAM→DE, `count-table`=count-table→DE |
| `--seqstrategy` | — | `single` or `paired` |
| `--total` | — | Total RNA mode (include ncRNA in GTF) |
| `--genomesize` | `10` | STAR genomeSAindexNbases |
| `--deseq2norm` | `DESeq2` | Normalization: `DESeq2` or `RPM` |
| `--foldchange` | `2.0` | Fold-change cutoff |
| `--fdr` | `1.0` | FDR cutoff |
| `--mask` | — | Mask FASTA for filtering |

---

### `degradome` — Degradome-seq (PARE/GMUCT)

Dual-alignment strategy (transcriptome + genome) with sPARTA peak-calling and CRI-based miRNA target cleavage analysis.

**Pipeline steps:**
1. Build transcriptome STAR index + genome STAR index with splice junctions
2. Per-sample: cutadapt trimming → transcriptome STAR → genome STAR alignment
3. Deduplicate reads, create library file (read counts)
4. sPARTA: build miRNA FASTA → target prediction, scoring, validation
5. CRI calculation based on CDS frame distribution (riboWaltz)

| Option | Default | Description |
|--------|---------|-------------|
| `--targets` | `all` | Transcript list for CRI analysis |
| `--sirnas` | `none` | Additional siRNA FASTA for targets |
| `--no-mapping` | — | Skip alignment |
| `--mapping-only` | — | Alignment only |

---

### `phasi` — phasiRNA Analysis

Identifies phased secondary siRNA (phasiRNA) loci with phasing score calculation.

**Pipeline steps:**
1. Per-sample: cutadapt → bowtie rRNA filtering → ShortStack (1000 multimaps, 0 mismatches)
2. Merge BAMs by group → extract exact-match reads
3. Merge plus/minus strands (minus shifted +2 nt)
4. Sliding window (10 periods) phasing score calculation
5. Output bedGraph + annotated results

| Option | Default | Description |
|--------|---------|-------------|
| `--period` | `21` | Phasing period size (19–26) |
| `--phasingscore` | `50` | Phasing score cutoff |
| `--mmap` | `u` | ShortStack multimap strategy |
| `--norm` | `rRNA,total` | Normalization |
| `--binsize` | `100` | Window size |
| `--no-mapping` | — | Skip alignment |

---

### `tt` — miRNA Truncation/Tailing

Iterative bowtie alignment (0–8 mismatches) to characterize miRNA truncation and tailing patterns.

**Pipeline steps:**
1. 0-mismatch bowtie → extract matched reads
2. Iterate: 1–8 mismatches bowtie → classify truncation/tailing variants
3. ShortStack alignment for remaining reads
4. bedtools intersect with miRNA loci → bubble plot visualization

| Option | Default | Description |
|--------|---------|-------------|
| `--mmap` | `u` | ShortStack multimap strategy |

---

### `ribo` — Ribo-seq (RIBO Taper Pipeline)

Full RIBO Taper workflow for translated ORF detection from Ribo-seq data.

**Pipeline steps:**
1. Bowtie2 contamination removal (rRNA, tRNA, snRNA, snoRNA)
2. Preprocess Ribo-seq reads
3. RNA-seq STAR 2-pass + StringTie transcriptome assembly
4. gffcompare → novel transcript filtering + gene_biotype annotation
5. RSEM quantification → expressed isoform filtering (TPM threshold)
6. STAR re-mapping with expressed annotation
7. RIBO Taper annotation files
8. Merge BAMs + Ribotaper.sh ORF detection
9. Output summary

| Option | Default | Description |
|--------|---------|-------------|
| `--rna-control` | *(required)* | RNA-seq control: `name=file1+file2...` |
| `--rna-treatment` | — | RNA-seq treatment: `name=file1+file2...` (repeatable) |
| `--ribo-control` | *(required)* | Ribo-seq control: `name=file1+file2...` |
| `--ribo-treatment` | — | Ribo-seq treatment: `name=file1+file2...` (repeatable) |
| `--contam` | *(required)* | Contamination FASTA for Bowtie2 index |
| `--ribo-len` | `24,25,26,27,28` | Ribo-seq read lengths |
| `--cutoffs` | `8,9,10,11,12` | RIBO Taper cutoffs |
| `--tpm-threshold` | `0` | Mean TPM threshold |
| `--ribotaper` | `~/software/ribotaper/bin` | Path to RIBO Taper installation |
| `--ribotaper-env` | `ribotaper` | Conda environment for RIBO Taper |

**Additional dependency:** [RiboTaper](https://github.com/hsinyenwu/RiboTaper) (manual install)

```bash
git clone https://github.com/hsinyenwu/RiboTaper.git ~/RiboTaper_v1.3
```

---

### `cips` — CiPS uORF Analysis

Downstream analysis after RIBO Taper: detects **translated upstream ORFs (uORFs)** in 5' UTRs using Ribo-seq P-site periodicity.

**Pipeline steps:**
1. Extract 5' UTRs from expressed GTF
2. Find all uORFs with ORFik
3. Compute frame-specific P-site counts (parallel)
4. Filter translated uORFs (min 1 aa; longer uORFs: ≥2 aa)
5. Deduplicate overlapping uORFs
6. Output 4 Excel files

| Option | Default | Description |
|--------|---------|-------------|
| `--gtf` | `<genome>_expressed.gtf` | Expressed GTF from ribo pipeline |
| `--fasta` | `reference/<genome>_chr_all.fasta` | Reference genome FASTA |
| `--psite` | *(required)* | P-site count file (`count chr start strand`) |
| `--min-inframe-counts` | `10` | Min in-frame Ribo-seq counts |
| `--min-inframe-perc` | `50` | Min in-frame percentage |
| `--min-psite-perc` | `30` | Min P-site % (longer uORFs) |
| `--gene-desc` | — | Gene description Excel for annotation |

---

### `chip` / `atac` — ChIP-seq / ATAC-seq

Both modes support **Genrich** (default) and **MACS3** peak calling, plus MACS3 `bdgdiff` for two-group differential peak analysis.

**ChIP-seq with MACS3 + bdgdiff:**
```
BAM → macs3 callpeak (per group, with Input) → macs3 bdgdiff → diff peaks
```

**ATAC-seq with MACS3 + bdgdiff:**
```
BAM → macs3 callpeak (per group, ATAC-specific params) → macs3 bdgdiff → diff peaks
```

| Option | Default | Description |
|--------|---------|-------------|
| `--peak-caller` | `genrich` | `genrich` or `macs3` |
| `--genome-size` | — | Effective genome size for MACS3 (e.g. `1.35e8`) |
| `--auc` | `20` | AUC threshold (Genrich) |
| `--qvalue` | `1.0` | Q-value cutoff |
| `--pvalue` | `0.01` | P-value cutoff |
| `--tss-distance` | `3000` | TSS distance for ChIPseeker annotation |
| `--no-mapping` | — | Skip alignment |
| `--mapping-only` | — | Alignment only |

**MACS3 bdgdiff output:**
- `{tag}_peaks.narrowPeak` — per-group peaks
- `diff_{g1}_vs_{g2}_cond1.bed` — group 1-specific
- `diff_{g1}_vs_{g2}_cond2.bed` — group 2-specific
- `diff_{g1}_vs_{g2}_common.bed` — shared peaks

---

### `wgbs` — Whole-Genome Bisulfite Sequencing

Bismark-based alignment and DMRcaller differential methylation analysis.

**Pipeline steps:**
1. Bismark genome preparation (if needed)
2. Per-sample: cutadapt trimming → Bismark alignment → deduplication → methylation extraction
3. Merge CpG reports → DMRcaller differential methylation

| Option | Default | Description |
|--------|---------|-------------|
| `--binsize` | `100` | Window size (bp) |
| `--minc` | `4` | Min reads per cytosine |
| `--no-mapping` | — | Skip alignment |
| `--mapping-only` | — | Alignment only |

---

### `clip` / `ts` — CLIP-seq / TS-CLIP-seq

STAR alignment with CLIPper peak calling for protein-RNA interaction sites.

| Option | Default | Description |
|--------|---------|-------------|
| `--foldchange` | `2.0` | Fold-change cutoff |
| `--pvalue` | `0.05` | P-value cutoff |
| `--no-mapping` | — | Skip alignment |
| `--mapping-only` | — | Alignment only |

---

### `ribometh` — RiboMeth-seq

STAR alignment to reference transcripts with RNAmodR.RiboMethSeq for 2'-O-methylation analysis.

**Pipeline steps:**
1. Build STAR index from reference transcriptome
2. Per-sample: cutadapt trimming → STAR alignment
3. Parse CIGAR, compute coverage and 5'/3' end distributions
4. RNAmodR.RiboMethSeq analysis

| Option | Default | Description |
|--------|---------|-------------|
| `--reference` | `genome` | Reference transcriptome FASTA |
| `--readlength` | `50` | Raw read length |
| `--coverage` | `1000` | Minimum coverage threshold |
| `--adaptor2` | `1` | Adaptor for read 2 |

---

### `risi` — risiRNA Analysis

ShortStack-based small RNA analysis optimized for risiRNA detection.

| Option | Default | Description |
|--------|---------|-------------|
| `--mmap` | `u` | ShortStack multimap strategy |
| `--norm` | `total` | Normalization method |
| `--binsize` | `10` | Window size |
| `--foldchange` | `1.5` | Fold-change cutoff |
| `--pvalue` | `0.01` | P-value cutoff |
| `--no-mapping` | — | Skip alignment |
| `--mapping-only` | — | Alignment only |

---

### `tf` — Two-Factor DE Analysis

Multi-condition differential expression/accessibility analysis with DESeq2. Supports three sub-modes:

- `--mode_tf mrna` — gene expression (mRNA-seq count table)
- `--mode_tf srna` — small RNA expression
- `--mode_tf chip` — ChIP-seq differential peaks (MACS3 bdgdiff)

**Sample format:** `groupName=label1,N1,label2,N2`  
Example: `WT=input,3,IP,3` (3 input replicates, 3 IP replicates)

| Option | Default | Description |
|--------|---------|-------------|
| `--mode_tf` | `mrna` | `mrna`, `srna`, or `chip` |
| `--foldchange` | `1.5` | Fold-change cutoff |
| `--pvalue` | `0.05` | P-value cutoff |
| `--norm` | `rRNA,total` | Normalization (srna mode) |
| `--binsize` | `100` | Window size (srna mode) |
| `--deseq2_norm` | `DESeq2` | Normalization (mrna mode) |
| `--genome-size` | — | Genome size for MACS3 (chip mode) |
| `--cutoff` | `3` | log2FC cutoff for bdgdiff (chip mode) |
| `--seq_strategy` | `paired` | Sequencing strategy (chip mode) |
| `--tss-distance` | `3000` | TSS distance for ChIPseeker |

---

## Common Options (All Modes)

| Option | Description |
|--------|-------------|
| `--outdir`, `-o` | Output directory (default: `./out`) |
| `--genome`, `-g` | Genome: `ath`, `osa`, `b73`, `gma`, `smo`, `bra`, `w22` |
| `--thread`, `-t` | Number of threads (default: `4`) |
| `--adaptor`, `-a` | 3' adaptor sequence (supports aliases: `truseq`, `illumina`, `srna`, `neb`, `nextera`) |
| `--control`, `-c` | Control samples: `name=file1+file2...` |
| `--treatment`, `-p` | Treatment samples (repeatable) |
| `--auto-install` | Auto-install missing dependencies (default) |
| `--no-auto-install` | Disable automatic dependency installation |

---

## Input Format

Samples are specified as `name=source` pairs. Sources can be:

- **SRA accessions:** `WT=SRR123456` (auto-download via `fasterq-dump`)
- **Local files:** `WT=data/sample.fq` or `WT=data/sample.fq.gz`
- **Multiple replicates:** `WT=rep1.fq+rep2.fq+rep3.fq` (concatenated)
- **Paired-end:** automatically detected for SRA accessions

---

## Project Structure

```
pRNASeqTools_py/
├── pRNASeqTools_run.py          # Entry point
├── environment.yaml             # Conda environment definition
├── README.md
├── LICENSE
├── prnaseqtools/
│   ├── __init__.py              # Package metadata (version 1.0.0)
│   ├── cli.py                   # CLI definition (argparse) + dispatch table
│   ├── auto_install.py          # Dependency auto-installer (conda/pip/git)
│   ├── validate_options.py      # Input validation + adaptor aliases
│   ├── input_parser.py          # Sample specification parser
│   ├── functions.py             # Utilities (run_cmd, download, unzip, revcomp)
│   ├── reference.py             # Genome reference handling (GFF parsing, indexing)
│   ├── precheck.py              # External dependency checker
│   ├── logging_setup.py         # Logging setup (tee to file + stderr)
│   └── modes/                   # Analysis mode implementations
│       ├── srna.py              # Small RNA-seq
│       ├── mrna.py              # mRNA-seq
│       ├── degradome.py         # Degradome-seq
│       ├── phasi.py             # phasiRNA analysis
│       ├── tt.py                # Truncation/tailing
│       ├── ribo.py              # RIBO Taper pipeline
│       ├── cips.py              # CiPS uORF analysis
│       ├── chip.py              # ChIP-seq
│       ├── atac.py              # ATAC-seq
│       ├── wgbs.py              # WGBS-seq
│       ├── clip.py              # CLIP-seq
│       ├── ts.py                # TS-CLIP-seq
│       ├── ribometh.py          # RiboMeth-seq
│       ├── risi.py              # risiRNA analysis
│       └── tf.py                # Two-factor DE
├── scripts/                     # R analysis scripts (24 R + 1 Python)
│   ├── checkPackages.R          # R package installer (BiocManager/GitHub/CRAN)
│   ├── DEG.R                    # mRNA differential expression
│   ├── DSR.R / DEM.R / DSG.R    # sRNA DE: repeat/miRNA/gene-level
│   ├── DST.R / DSP.R            # sRNA DE: TE/promoter-level
│   ├── CRI.R                    # Degradome cleavage ratio index
│   ├── ribo.R                   # Ribo-seq frame analysis
│   ├── ribotaper_filter_gtf.R   # GTF filtering for RIBO Taper
│   ├── ribotaper_filter_rsem.R  # RSEM isoform filtering
│   ├── cips_uORF.R              # CiPS uORF detection
│   ├── CLIP.R                   # CLIP-seq peak analysis
│   ├── RNAmodR.R                # RiboMeth-seq analysis
│   ├── DMRcaller.R              # WGBS differential methylation
│   ├── bubble_plot.R            # miRNA truncation/tailing visualization
│   ├── chipseeker.R             # ChIP/ATAC peak annotation
│   ├── DSF.R                    # sRNA DE: fold-change analysis
│   ├── tf_gene.R                # Two-factor DE (gene)
│   ├── tf_mirna.R               # Two-factor DE (miRNA)
│   ├── tf_mrna.R                # Two-factor DE (mRNA)
│   ├── tf_promoter.R            # Two-factor DE (promoter)
│   ├── tf_srna.R                # Two-factor DE (sRNA)
│   ├── tf_te.R                  # Two-factor DE (TE)
│   └── getPrimaryTranscript.py  # Primary transcript extraction
└── reference/                   # Genome files (not versioned, contact author)
```

---

## Dependencies

### Core tools (via conda)

| Category | Tools |
|----------|-------|
| **Runtimes** | Python ≥3.9, R ≥4.0, Perl |
| **Aligners** | STAR ≥2.7, bowtie, bowtie2, Bismark, ShortStack ≥3.0 |
| **Processing** | cutadapt, samtools ≥1.0, htslib, bedtools, gffread, deepTools |
| **Counting** | featureCounts (subread), RSEM, StringTie |
| **Peak calling** | Genrich, MACS3 |
| **Utilities** | sra-tools, gffcompare, umi_tools, ucsc-bedgraphtobigwig |
| **Python** | numpy, scipy (stdlib-only pipeline; numpy/scipy for sPARTA) |
| **R (conda)** | DESeq2, DMRcaller, RNAmodR.RiboMethSeq, pheatmap, dplyr, devtools |

### R packages (via checkPackages.R)

| Package | Source | Mode |
|---------|--------|------|
| riboWaltz | GitHub: `LabTranslationalArchitectomics/riboWaltz` | degradome |
| NMF | GitHub: `renozao/NMF` (devel) | — |
| Seurat | GitHub: `satijalab/seurat` | srna (sc) |
| ORFik | Bioconductor | cips |

### Manual install

- **CLIPper** — `clip` / `ts` modes: `git clone https://github.com/YeoLab/clipper.git && python setup.py install`
- **RiboTaper** — `ribo` mode: `git clone https://github.com/hsinyenwu/RiboTaper.git ~/RiboTaper_v1.3`
- **Reference files** — contact the author

---

## License

Contact the author for licensing information.
