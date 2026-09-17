# pRNASeqTools

**Integrated High-throughput Sequencing Data Analysis for Plant** — Python3 rewrite of the original Perl pipeline.

Author: Dr. Chenjiang You  
Version: 1.0.0

---

## Overview

pRNASeqTools is a comprehensive NGS data analysis pipeline for plant genomics. It provides **18 analysis modes** covering small RNA, mRNA, lncRNA, single-cell, epigenomics, ribosome profiling, and more — all accessible through a single command-line interface.

### Analysis Modes

| Mode | Description | Key Tools |
|------|-------------|-----------|
| `srna` | Small RNA-seq — bulk & single-cell | ShortStack, bowtie, fastp, umi_tools |
| `mrna` | mRNA-seq — DE analysis | STAR, featureCounts, DESeq2 |
| `lncrna` | lncRNA-seq — transcriptome assembly & lncRNA classification | STAR, StringTie, FEELnc / PLEK2 / ORF-fallback |
| `sc` | Single-cell RNA-seq — Seurat clustering + marker genes | STARsolo, Seurat, Harmony |
| `degradome` | Degradome-seq (PARE/GMUCT) — miRNA target cleavage | STAR (dual-alignment), sPARTA, riboWaltz |
| `phasi` | phasiRNA analysis — phased secondary siRNA detection | ShortStack, bowtie |
| `tt` | miRNA truncation/tailing analysis | bowtie (iterative), ShortStack |
| `ribo` | Ribo-seq — translated ORF detection + Translation Efficiency | Bowtie2, STAR, StringTie, RSEM, RiboTaper, DESeq2 |
| `cips` | CiPS uORF analysis — translated upstream ORF detection | ORFik, GenomicFeatures (R) |
| `chip` | ChIP-seq — peak calling | bowtie2, Genrich / MACS3 |
| `atac` | ATAC-seq — open chromatin analysis | bowtie2, Genrich / MACS3, deepTools |
| `tf` | Two-factor DE — multi-condition differential analysis (mRNA/sRNA/ChIP) | DESeq2, DiffBind, bdgdiff, ChIPseeker |
| `wgbs` | Whole-genome bisulfite — differential methylation | Bismark, DMRcaller |
| `clip` | CLIP-seq — protein-RNA interaction | STAR, CLIPper |
| `ts` | TS-CLIP-seq — target-specific CLIP | STAR, CLIPper |
| `ribometh` | RiboMeth-seq — 2'-O-methylation analysis | STAR, RNAmodR.RiboMethSeq |
| `risi` | risiRNA analysis | ShortStack, bowtie |

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

The pipeline also supports **auto-install** of missing dependencies (enabled by default). Add `--no-auto-install` to disable. Auto-install covers core tools (conda), R Bioconductor packages (DESeq2, DiffBind, ChIPseeker, clusterProfiler, enrichplot, DMRcaller, etc.), and genome annotation libraries.

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
  --ribo-control "Ribo=SRR333333"

# Ribo-seq — custom contamination
python pRNASeqTools_run.py ribo \
  --rna-control "WT=SRR111111" \
  --ribo-control "Ribo=SRR333333" \
  --contam "rRNA.fasta,tRNA.fasta,snRNA.fasta"

# Ribo-seq — auto-resume from last completed step
python pRNASeqTools_run.py ribo \
  --rna-control "WT=SRR111111" \
  --ribo-control "Ribo=SRR333333" \
  --outdir ./my_ribo_run

# Ribo-seq — force restart from step 8 (re-run RIBO Taper with new params)
python pRNASeqTools_run.py ribo \
  --rna-control "WT=SRR111111" \
  --ribo-control "Ribo=SRR333333" \
  --restart-step 8

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

# ChIP-seq differential peaks — DiffBind + dual-factor + mito normalization (RECOMMENDED)
python pRNASeqTools_run.py tf -c "WT=input,2,IP,2" -p "KO=input,2,IP,2" \
  --mode_tf chip --genome ath --chip-norm mito

# ChIP-seq differential peaks — MACS3 bdgdiff (no replicates)
python pRNASeqTools_run.py tf -c "WT=input,2,IP,2" -p "KO=input,2,IP,2" \
  --mode_tf chip --chip-method bdgdiff --genome-size 1.35e8

# ATAC-seq (MACS3)
python pRNASeqTools_run.py atac --peak-caller macs3 --genome-size 1.35e8 \
  --treatment "ATAC=data/atac.bam"
```

---

## Mode Details

### `srna` — Small RNA-seq

Full pipeline for small RNA-seq, supporting both bulk and single-cell (UMI-based) modes.

**Pipeline steps:**
1. SRA download (if needed) → fastp adapter trimming (retain 18–42 nt)
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
2. Per-sample: SRA download → fastp trimming → optional mask filtering → STAR alignment
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

**Output files:**

| File | Description |
|------|-------------|
| `{tag}.bam` | STAR 排序 BAM |
| `{tag}.bw` | CPM 归一化 bigWig (bamCoverage, bs=5, MAPQ≥10) |
| `{tag}.txt` | Per-sample count table (`Gene / Count / Length`)，DESeq2 输入 |
| `DEG_overview.pdf` | 3-panel 合页：top 1000 基因 heatmap + PCA plot + 样本距离热力图 |
| `{treat}vs{ctrl}.total.csv` | 完整 DESeq2 结果：`baseMean / log2FC / lfcSE / stat / pvalue / padj` + **每个样本的 TPM** |
| `{treat}vs{ctrl}.total.upregulated.csv` | padj < FDR 且 log2FC ≥ log₂(`foldchange`) 的上调基因 + TPM 列 |
| `{treat}vs{ctrl}.total.downregulated.csv` | padj < FDR 且 log2FC ≤ -log₂(`foldchange`) 的下调基因 + TPM 列 |
| `{treat}vs{ctrl}.total.bin.txt` | 基因功能 bin Fisher's exact test（仅当 `reference/{genome}.BIN` 存在时生成） |

> **TPM 计算**：从 `{tag}.txt` 的 `Length` 列提取基因长度，`TPM = (counts / len_kb) / colSum(counts / len_kb) × 10⁶`。

---

### `lncrna` — lncRNA-seq

De novo transcriptome assembly and lncRNA classification pipeline. Supports three classifier backends and an ORF-based heuristic fallback.

**Pipeline steps (mode-dependent):**
1. STAR genome index + fastp trimming (if `mode_lncrna` includes mapping)
2. Per-sample STAR alignment → BAM files
3. StringTie transcriptome assembly (per-sample + merge)
4. gffcompare → novel transcript filtering (class codes `u`, `i`, `x`, `o`)
5. lncRNA classification:
   - **FEELnc** (default): trained on user genome, `train() → classify() → filter()`
   - **PLEK2**: 2-kmer SVM, runs via `python3 PLEK2.py -i <fasta> -m pl`
   - **ORF-fallback**: ≥300 nt ORF + ≥50% transcript coverage = coding
6. featureCounts quantification → DESeq2 differential expression

**Modes (`--mode_lncrna`):**
| Mode | Description |
|------|-------------|
| `whole` | Full pipeline: mapping → assembly → classification → DE (default) |
| `mapping-only` | STAR alignment only |
| `assemble-only` | Assembly + classification, no mapping |
| `de-only` | Classification + DE from pre-assembled GTF |
| `count-table` | DE analysis from pre-computed count table |

| Option | Default | Description |
|--------|---------|-------------|
| `--mode_lncrna` | `whole` | Pipeline mode (see above) |
| `--seqstrategy` | — | `single` or `paired` |
| `--classifier` | `feelnc` | `feelnc`, `plek2`, or `fallback` (ORF heuristic) |
| `--feelnc-dir` | — | Custom FEELnc install path (if not in conda) |
| `--plek2-dir` | `~/PLEK2` | Custom PLEK2 install path |
| `--min-fpkm` | `0.5` | Min FPKM for expressed transcripts |
| `--min-tx-len` | `200` | Min transcript length (bp) |
| `--gtf` | — | Pre-assembled GTF (for `assemble-only` / `de-only`) |
| `--foldchange` | `2.0` | Fold-change cutoff |
| `--pvalue` | `0.01` | P-value cutoff |
| `--fdr` | `1.0` | FDR cutoff |
| `--genomesize` | `10` | STAR genomeSAindexNbases |

**Dependencies:**
- FEELnc via conda: `conda install -c bioconda feelnc`
- PLEK2 via manual clone: `git clone https://github.com/emanlee/plek2 ~/PLEK2` (requires `keras==2.4.3`, `tensorflow==2.4.1`, `numpy==1.19.2`)
- gffcompare, StringTie, STAR, featureCounts

---

### `sc` — Single-Cell RNA-seq

STARsolo-based mapping followed by Seurat integration, clustering, and marker gene identification. Supports multiple integration backends.

**Pipeline steps:**
1. STAR genome index + fastp trimming (if input is FASTQ)
2. STARsolo mapping + quantification (`--quantMode GeneCounts`, cell barcode CB/BZ tags)
3. Seurat object creation (or direct import from `.rds` / `.mtx` / `.csv` / STARsolo dir)
4. Per-sample QC: filter low-cells/low-features/high-mt cells
5. Per-sample normalization + variable feature detection
6. **Integration** (optional):
   - **Seurat anchors** (`integration=seurat`, default): `FindIntegrationAnchors` → `IntegrateData`
   - **Harmony** (`integration=harmony`, faster for large datasets): `HarmonyMatrix`
   - **none**: skip integration (per-sample or single-sample analysis)
7. PCA → UMAP → clustering (Louvain/Leiden)
8. Marker gene identification (`FindAllMarkers` / `FindClusterMarkers`)
9. Optional: Doublet detection, pseudotime analysis

**Input formats supported:**
- FASTQ → STARsolo
- Cell-tagged BAM (CB:Z:/BZ:Z: tags) → STARsolo soloQuant
- Count matrices: `.rds` (R native), `.mtx` / `.mtx.gz` (MatrixMarket), `.tsv`/`.txt`/`.csv`
- STARsolo output directory

| Option | Default | Description |
|--------|---------|-------------|
| `--mode_sc` | `whole` | `whole`, `mapping-only`, or `count-table` |
| `--seqstrategy` | — | `single` or `paired` |
| `--mincells` | `3` | Min cells a gene must be expressed in |
| `--minfeatures` | `200` | Min features per cell |
| `--maxfeatures` | `5000` | Max features per cell |
| `--pctmt` | `20` | Max % mitochondrial content |
| `--npcs` | `30` | Number of PCs for dimensionality reduction |
| `--nclusters` | `0` | Target number of clusters (0 = auto-optimize) |
| `--resolution` | `0.5` | Clustering resolution |
| `--markerminpct` | `0.25` | Min % of cells expressing marker |
| `--markerlogfc` | `0.25` | Min log₂FC for markers |
| `--pseudotime` | — | Enable pseudotime analysis |
| `--doublet` | `0` | Doublet detection (0=off, or expected % e.g. 7.5) |
| `--integration` | `seurat` | `seurat`, `harmony`, or `none` |

**Dependencies:**
- R: `Seurat` (GitHub: `satijalab/seurat`), `harmony` (Bioconductor)
- Python: `STAR`, `fastp`

---

### `degradome` — Degradome-seq (PARE/GMUCT)

Dual-alignment strategy (transcriptome + genome) with sPARTA peak-calling and CRI-based miRNA target cleavage analysis.

**Pipeline steps:**
1. Build transcriptome STAR index + genome STAR index with splice junctions
2. Per-sample: fastp trimming → transcriptome STAR → genome STAR alignment
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
1. Per-sample: fastp → bowtie rRNA filtering → ShortStack (1000 multimaps, 0 mismatches)
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

Full RIBO Taper workflow for translated ORF detection from Ribo-seq data, followed by **Translation Efficiency (TE)** analysis.

**Pipeline steps (12 steps total):**
1. Bowtie2 contamination removal (rRNA, tRNA, snRNA, snoRNA)
2. Preprocess Ribo-seq reads
3. RNA-seq STAR 2-pass + StringTie transcriptome assembly
4. gffcompare → novel transcript filtering + gene_biotype annotation
5. RSEM quantification → expressed isoform filtering (TPM threshold)
6. STAR re-mapping with expressed annotation (Ribo-seq + RNA-seq)
7. RIBO Taper annotation files (create_annotations_files.bash)
8. Merge BAMs → metaplots → interactive parameter confirmation → Ribotaper.sh ORF detection
9. P-site analysis & visualization (frame distribution + bedtools-closest metagene plots)
10. **Translation Efficiency (TE)** calculation per sample pair (RSEM-based ribo/RNA ratio) + **rPS index** (start codon P-sites / other-CDS P-sites)
11. TE statistical testing (DESeq2) + visualization (TE_stats.R)
12. Final output summary

**Auto-resume:** The pipeline automatically detects the last completed step from log files and resumes from where it left off. Use `--restart-step N` to force restart from a specific step (1–12).

**Key features:**
- Length distribution plots (`Ribo_length_distributions.pdf`) generated after mapping (steps 1–8); skipped when running step 9+ only
- Metaplots generated before RIBO Taper for parameter selection
- Interactive confirmation of `ribo-len` and `cutoffs` after metaplots
- Frame computation uses transcript_id from `start_stop_FAR.bed` match — not limited to `.1` transcripts
- Metagene plots use bedtools closest (`P_sites_all` vs `start_stop_FAR.bed`) with RiboTaper's distance formula
- Per-sample metagene plots with frame-colored stacked histograms (green→yellow→red gradient) + line plot overlays
- Combined "All reads" page plus per-sample pages (2×2 layout: start codon, stop codon, overlay, stats)
- GTF used preferentially; GFF auto-converted if needed via `gffread -T`
- Ribo-seq BAM filtered to R1-only; RNA-seq BAM retains paired-end
- TE calculation pairs samples by position: `all_ribo_tags[i] ↔ all_rna_tags[i]`
- rPS index computed per Ribo-seq sample (`rPS_results/rps_table_{tag}.tsv`)

**TE statistical methods (`--te-method`):**

| Method | Description |
|--------|-------------|
| `separate` | Two independent DESeq2: Ribo-seq DE vs RNA-seq DE, TE change = `delta_logFC` difference |
| `joint` | Single dual-factor DESeq2: `~ condition + data_type + condition:data_type`, interaction term is TE change |
| `both` | Run both methods (default, combined output) |
| `none` | Skip DESeq2, TE plotting + rPS only |

| Option | Default | Description |
|--------|---------|-------------|
| `--rna-control` | *(required)* | RNA-seq control: `name=file1+file2...` |
| `--rna-treatment` | — | RNA-seq treatment: `name=file1+file2...` (repeatable) |
| `--ribo-control` | *(required)* | Ribo-seq control: `name=file1+file2...` |
| `--ribo-treatment` | — | Ribo-seq treatment: `name=file1+file2...` (repeatable) |
| `--contam` | `reference/{genome}_contam4.fa` | Contamination FASTA for Bowtie2 index |
| `--ribo-len` | `24,25,26,27,28` | Ribo-seq read lengths (must match cutoffs count) |
| `--cutoffs` | `8,9,10,11,12` | RIBO Taper cutoffs (must match ribo-len count) |
| `--tpm-threshold` | `0` | Mean TPM threshold for expressed isoform filtering |
| `--ribotaper` | `~/software/ribotaper/bin` | Path to RIBO Taper installation |
| `--ribotaper-env` | `ribotaper` | Conda environment for RIBO Taper |
| `--restart-step` | — | Force restart from step N (1–12), overrides auto-detection |
| `--te-method` | `both` | TE change detection: `separate`, `joint`, `both`, `none` |

**TE output structure:**
```
outdir/
├── TE_results/
│   ├── RiboTag__RNATag/          # Per sample pair (ribo ↔ rna by position)
│   │   ├── te_table.tsv          # RSEM TE (ribo RPM / rna RPM) per gene
│   │   └── ...
│   ├── te_combined.tsv           # All pairs merged, per-sample TE values
│   ├── te_combined_log2.tsv      # log2-transformed combined table
│   └── TE_summary.pdf            # TE boxplot + gene-level summary
├── rPS_results/
│   ├── rps_table_{tag}.tsv       # Per Ribo-seq sample rPS index
│   └── rps_summary.tsv           # All samples merged
├── ribotaper_results/            # RIBO Taper ORF detection output
└── metagene_plots.pdf            # P-site metagene visualization
```

**Additional dependency:** [RiboTaper](https://github.com/hsinyenwu/RiboTaper) (manual install)

```bash
git clone https://github.com/hsinyenwu/RiboTaper.git ~/RiboTaper_v1.3
```

**TE DESeq2 requires these R packages** (installed via `checkPackages.R`):
- `emmeans`, `car`, `agricolae`, `multcomp`, `ggpubr` — Tukey HSD and post-hoc tests
- `ComplexHeatmap` — combined TE heatmaps
- `DESeq2` — differential testing (separate + joint models)

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

Both modes support **Genrich** (default) and **MACS3** peak calling. ChIP-seq differential peak analysis via `tf --mode_tf chip` supports two methods: **DiffBind** (default, DESeq2-based with replicates) and **MACS3 bdgdiff** (no replicates).

#### Single-sample peak calling (`chip` mode)

```
BAM → Genrich/MACS3 → narrowPeak → BED → ChIPseeker annotation → GO enrichment
```

#### Differential peak calling (`tf --mode_tf chip`, DiffBind — default)

> **DiffBind 3.x compatibility note:** DiffBind 3.x's `dba.count()` merges Input reads into each IP sample, producing a count matrix with only IP columns. To support the full dual-factor model (`~ Condition + Factor + Condition:Factor`, which requires separate Input columns), this pipeline uses **manual counting via `bedtools multicov`** from all dedup BAMs (IP + Input), bypassing DiffBind's internal count matrix. DiffBind is still used for its excellent **peak merging logic** (`minOverlap=2` consensus peaks).

**Pipeline steps:**
1. BAM resolution (auto-find `.sorted.bam` / `.sorted.dedup.bam`)
2. **Picard MarkDuplicates** — auto-generate dedup BAMs (with `AddOrReplaceReadGroups`)
3. Per-sample MACS3 `callpeak` → consensus peaks via DiffBind
4. **Spike-in normalization** (mito/chloro/rDNA reads from `idxstats` or `samtools view -c region`)
   - Scale factor = `ip_spike_in / max(ip_spike_in)` — **only IP samples**, Input not scaled
5. **bedtools multicov** — count reads in each consensus peak across ALL BAMs (IP + Input)
6. **DESeq2** — dual-factor model `~ Condition + Factor + Condition:Factor` (or affinity `~ Condition` for IP-only)
7. Export consensus / UP / DOWN peak BED files
8. **ChIPseeker** annotation + **clusterProfiler** GO enrichment (三类 peaks 分别分析)
9. **bigWig generation** (bamCoverage with `--scaleFactor = ip_scale` for IP, 1.0 for Input)

**Normalization options (spike-in based, recommended over DESeq2 default MoR):**

| Method | Statistic | Description |
|--------|-----------|-------------|
| `deseq2` | DESeq2 median-of-ratios | Default, data-driven |
| `total` | Total mapped reads | Per-sample library size |
| `mito` | Mitochondrial reads (chrM/MT) | Spike-in: Input ratio ≈ stable across conditions |
| `chloro` | Chloroplast reads (chrC/chrCP) | Spike-in: Input ratio ≈ stable across conditions |
| `rdna` | rDNA region reads (chr2:0–10500, chr3:14193500–14204500 for Arabidopsis) | Precise genomic intervals |

**DESeq2 sizeFactor direction (critical):**
- `sizeFactor = ip_scale` (smaller mito → smaller sizeFactor → `normalized = raw / sizeFactor` = upregulated in DESeq2)
- `bamCoverage --scaleFactor = 1/ip_scale` (complementary, IGV BW matches DESeq2 direction)

**Dual-factor interaction term (the core output):**
```
Design: ~ Condition + Factor + Condition:Factor
Coefficients: β₀ (baseline Input) + β₁ (Input condition diff) + β₂ (IP vs Input enrichment in ref) + β₃ (IP enrichment diff between conditions)
β₃ > 0, padj significant → treatment-specific binding
β₃ < 0, padj significant → control-specific binding
```

**DiffBind output (`diffbind_results/`):**
- `samplesheet.csv` — DiffBind sample sheet
- `norm_factors.tsv` — IP scaling factors (for non-deseq2 norm)
- `peak_counts.tsv` — bedtools multicov counts (chr, start, end, ... + per-sample columns)
- `sample_metadata.tsv` — condition/factor per sample
- `consensus_peaks.bed` — DiffBind merged peaks
- `DiffBind_dual_factor_interaction.tsv` — DESeq2 results with `category` (UP/DOWN/NS)
- `DiffBind_dual_factor_volcano.pdf` — volcano plot
- `DiffBind_consensus_peaks.bed` / `DiffBind_UP_peaks.bed` / `DiffBind_DOWN_peaks.bed` — BED files for annotation
- `DiffBind_*_annotation.txt` / `*_annotation_pie.pdf` / `*_go_enrichment.txt` / `*_go_dotplot.pdf` — ChIPseeker outputs
- `bw/` — normalized bigWig files (IP scaled, Input unscaled)

| Option | Default | Description |
|--------|---------|-------------|
| `--peak-caller` | `genrich` | `genrich` or `macs3` (single-sample mode) |
| `--chip-method` | `diffbind` | `diffbind` (DESeq2 with replicates) or `bdgdiff` (MACS3, no replicates) |
| `--chip-analysis` | `dual_factor` | `dual_factor` (IP+Input joint DESeq2, interaction term) or `affinity` (IP-only) |
| `--chip-norm` | `deseq2` | `deseq2`, `total`, `mito`, `chloro`, `rdna` |
| `--genome-size` | — | Effective genome size for MACS3 (`bdgdiff` mode) |
| `--auc` | `20` | AUC threshold (Genrich) |
| `--qvalue` | `1.0` | Q-value cutoff |
| `--pvalue` | `0.05` | P-value cutoff for significance |
| `--foldchange` | `1.5` | log₂FC cutoff for significance |
| `--tss-distance` | `3000` | TSS distance for ChIPseeker annotation |
| `--no-mapping` | — | Skip alignment |
| `--mapping-only` | — | Alignment only |

#### MACS3 bdgdiff output (legacy mode, `--chip-method bdgdiff`):
- `{tag}_peaks.narrowPeak` — per-group peaks
- `diff_{g1}_vs_{g2}_cond1.bed` — group 1-specific
- `diff_{g1}_vs_{g2}_cond2.bed` — group 2-specific
- `diff_{g1}_vs_{g2}_common.bed` — shared peaks

---

### `wgbs` — Whole-Genome Bisulfite Sequencing

Bismark-based alignment and DMRcaller differential methylation analysis.

**Pipeline steps:**
1. Bismark genome preparation (if needed)
2. Per-sample: fastp trimming → Bismark alignment → deduplication → methylation extraction
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
2. Per-sample: fastp trimming → STAR alignment
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
- `--mode_tf chip` — ChIP-seq differential peaks (**DiffBind** default, MACS3 bdgdiff legacy)

**Sample format:** `groupName=label1,N1,label2,N2`  
Example: `WT=input,2,IP,2` (2 input replicates, 2 IP replicates)

| Option | Default | Description |
|--------|---------|-------------|
| `--mode_tf` | `mrna` | `mrna`, `srna`, or `chip` |
| `--foldchange` | `1.5` | Fold-change cutoff |
| `--pvalue` | `0.05` | P-value cutoff |
| `--norm` | `rRNA,total` | Normalization (srna mode) |
| `--binsize` | `100` | Window size (srna mode) |
| `--deseq2_norm` | `DESeq2` | Normalization (mrna mode) |
| `--chip-method` | `diffbind` | `diffbind` (DESeq2, with replicates) or `bdgdiff` (MACS3, no replicates) |
| `--chip-analysis` | `dual_factor` | `dual_factor` (IP+Input interaction term) or `affinity` (IP-only) |
| `--chip-norm` | `deseq2` | `deseq2`, `total`, `mito`, `chloro`, `rdna` (spike-in based) |
| `--genome-size` | — | Genome size for MACS3 (`bdgdiff` mode) |
| `--cutoff` | `3` | log₂FC cutoff for bdgdiff |
| `--seq_strategy` | `paired` | Sequencing strategy |
| `--tss-distance` | `3000` | TSS distance for ChIPseeker annotation |

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

- **SRA accessions:** `WT=SRR123456` (auto-download via `prefetch` + local `fasterq-dump` conversion)
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
│       ├── lncrna.py            # lncRNA-seq (new)
│       ├── sc.py                # Single-cell RNA-seq (new)
│       ├── degradome.py         # Degradome-seq
│       ├── phasi.py             # phasiRNA analysis
│       ├── tt.py                # Truncation/tailing
│       ├── ribo.py              # RIBO Taper + TE pipeline
│       ├── cips.py              # CiPS uORF analysis
│       ├── chip.py              # ChIP-seq peak calling
│       ├── atac.py              # ATAC-seq
│       ├── wgbs.py              # WGBS-seq
│       ├── clip.py              # CLIP-seq
│       ├── ts.py                # TS-CLIP-seq
│       ├── ribometh.py          # RiboMeth-seq
│       ├── risi.py              # risiRNA analysis
│       └── tf.py                # Two-factor DE (mRNA/sRNA/ChIP)
├── scripts/                     # R analysis scripts (30 R + 1 Python)
│   ├── checkPackages.R          # R package installer (BiocManager/GitHub/CRAN)
│   ├── DEG.R                    # mRNA differential expression
│   ├── DSR.R / DEM.R / DSG.R    # sRNA DE: repeat/miRNA/gene-level
│   ├── DST.R / DSP.R / DSF.R    # sRNA DE: TE/promoter/fold-change
│   ├── CRI.R                    # Degradome cleavage ratio index
│   ├── ribo.R                   # Ribo-seq frame analysis
│   ├── ribotaper_filter_gtf.R   # GTF filtering for RIBO Taper
│   ├── ribotaper_filter_rsem.R  # RSEM isoform filtering
│   ├── TE_DA.R                  # TE differential abundance (separate DESeq2)
│   ├── TE_DS.R                  # TE differential synthesis (joint DESeq2)
│   ├── TE_stats.R               # TE summary statistics + visualization
│   ├── cips_uORF.R              # CiPS uORF detection
│   ├── lncrna.R                 # lncRNA classification (FEELnc/PLEK2/orf)
│   ├── sc_analysis.R            # Single-cell Seurat integration + clustering
│   ├── CLIP.R                   # CLIP-seq peak analysis
│   ├── RNAmodR.R                # RiboMeth-seq analysis
│   ├── DMRcaller.R              # WGBS differential methylation
│   ├── bubble_plot.R            # miRNA truncation/tailing visualization
│   ├── chip_diffbind.R          # ChIP-seq DiffBind + DESeq2 dual-factor model
│   ├── chipseeker.R             # ChIP/ATAC peak annotation + GO enrichment
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
| **Aligners** | STAR ≥2.7 (STARsolo), bowtie, bowtie2, Bismark, ShortStack ≥3.0 |
| **Processing** | fastp, samtools ≥1.0, htslib, bedtools, gffread, deepTools, picard-slim (openjdk ≥17) |
| **Counting** | featureCounts (subread), RSEM, StringTie |
| **Peak calling** | Genrich, MACS3 |
| **ChIP-seq** | DiffBind ≥3.0, ChIPseeker, clusterProfiler, enrichplot |
| **lncRNA** | FEELnc (bioconda), gffcompare |
| **Utilities** | sra-tools, umi_tools, ucsc-bedgraphtobigwig |
| **Python** | numpy, scipy (stdlib-only pipeline; numpy/scipy for sPARTA) |
| **R (conda)** | DESeq2, DMRcaller, RNAmodR.RiboMethSeq, pheatmap, dplyr, devtools, emmeans, car, agricolae, multcomp, ggpubr, ComplexHeatmap |

### R packages (via checkPackages.R)

| Package | Source | Mode |
|---------|--------|------|
| riboWaltz | GitHub: `LabTranslationalArchitectomics/riboWaltz` | degradome |
| NMF | GitHub: `renozao/NMF` (devel) | — |
| Seurat | GitHub: `satijalab/seurat` | `sc`, `srna` (sc mode) |
| ORFik | Bioconductor | cips |
| harmony | Bioconductor | `sc` |
| DiffBind | Bioconductor ≥3.0 | `chip`, `tf` (chip mode) |
| ChIPseeker | Bioconductor | `chip`, `atac`, `tf` |
| clusterProfiler | Bioconductor | `chip`, `atac`, `tf` |
| RNAmodR.RiboMethSeq | Bioconductor | `ribometh` |
| DMRcaller | Bioconductor | `wgbs` |

### Manual install

- **CLIPper** — `clip` / `ts` modes: `git clone https://github.com/YeoLab/clipper.git && python setup.py install`
- **RiboTaper** — `ribo` mode: `git clone https://github.com/hsinyenwu/RiboTaper.git ~/RiboTaper_v1.3`
- **PLEK2** — `lncrna` mode (optional): `git clone https://github.com/emanlee/plek2 ~/PLEK2`
- **Reference files** — contact the author

---

## License

Contact the author for licensing information.
