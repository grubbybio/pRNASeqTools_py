#!/usr/bin/env Rscript
#
# ChIPseeker Peak Annotation & clusterProfiler GO Enrichment
#
# Processes a SINGLE peak file — designed to be called once per peak set,
# making it reusable by both tf (5 peak sets) and chip (1 peak set) modules.
#
# Usage:
#   Rscript chipseeker.R <genome> <prefix> <peak_file> <peak_name> [tss_distance]
#
# Example:
#   Rscript chipseeker.R ath /path/to/pkg WT_peaks.narrowPeak WT_peaks
#   Rscript chipseeker.R ath /path/to/pkg WT_peaks.narrowPeak WT_peaks 2000
#
# Arguments:
#   genome       - genome identifier (e.g. ath → uses org.At.tair.db)
#   prefix       - project prefix path (contains reference/ directory)
#   peak_file    - peak file path (BED/narrowPeak format)
#   peak_name    - label used for output file naming
#   tss_distance - distance upstream of TSS for annotation (default: 3000)
#
# Output files (written to CWD):
#   - {peak_name}_annotation.txt       — ChIPseeker annotation table
#   - {peak_name}_go_enrichment.txt    — GO enrichment results (if OrgDb available)
#   - {peak_name}_annotation_pie.pdf   — Annotation pie chart
#   - {peak_name}_go_dotplot.pdf       — GO enrichment dotplot

suppressPackageStartupMessages({
    library(ChIPseeker)
    library(GenomicFeatures)
    # 新包：Bioc 3.20+ 的 makeTxDb*() 函数已整体迁移到 txdbmaker
    if (requireNamespace("txdbmaker", quietly = TRUE)) library(txdbmaker)
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
})

# ── Parse command-line arguments ──────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript chipseeker.R <genome> <prefix> <peak_file> <peak_name> [tss_distance]")
}

genome       <- args[1]
prefix       <- args[2]
peak_file    <- args[3]
peak_name    <- args[4]
tss_distance <- ifelse(length(args) >= 5, as.numeric(args[5]), 3000)

# ── Check peak file exists ────────────────────────────────────────────────
if (!file.exists(peak_file)) {
    stop("Peak file not found: ", peak_file)
}

# ── Resolve GFF path ──────────────────────────────────────────────────────
gff_path <- file.path(prefix, "reference", paste0(genome, "_genes.gff"))
if (!file.exists(gff_path)) {
    stop("GFF file not found: ", gff_path)
}

# ── Build TxDb from GFF ──────────────────────────────────────────────────
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("ChIPseeker Peak Annotation —", peak_name, "\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("  Genome   :", genome, "\n")
cat("  GFF      :", gff_path, "\n")
cat("  Peak file:", peak_file, "\n")
cat("  Peak name:", peak_name, "\n")
cat("  TSS region: [-", tss_distance, ", 0]\n\n", sep = "")

cat("Building TxDb from GFF...\n")
# ── TxDb 构建（Bioc 版本兼容 + 缓存 RDS 避免重复构建） ──
# tf.py 对 5 个 peak set 逐个调用 chipseeker.R，每次独立起 R 进程 ——
# 用缓存（基于 gff mtime）让后续调用直接 readRDS，秒级返回。
#
# Bioc 3.20+ (GenomicFeatures ≥1.57): makeTxDb*() 已迁移到 txdbmaker 新包
# Bioc 3.18–3.19: makeTxDbFromGRanges 在 GenomicFeatures 导出
# Bioc ≤3.17:     makeTxDbFromGFF 在 GenomicFeatures 导出
buildTxDb <- function(gff_path) {
  # 1. 首选：官方新包 txdbmaker (Bioc 3.20+)
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    return(txdbmaker::makeTxDbFromGFF(gff_path))
  }
  # 2. 次选：GenomicFeatures 直接导出 makeTxDbFromGFF (旧版)
  if (exists("makeTxDbFromGFF", where = "package:GenomicFeatures", mode = "function")) {
    return(GenomicFeatures::makeTxDbFromGFF(gff_path))
  }
  # 3. 最后：中间版本可能只有 makeTxDbFromGRanges
  if (exists("makeTxDbFromGRanges", where = "package:GenomicFeatures", mode = "function")) {
    return(GenomicFeatures::makeTxDbFromGRanges(rtracklayer::import(gff_path)))
  }
  if ("makeTxDbFromGRanges" %in% ls(getNamespace("GenomicFeatures"), all.names = TRUE)) {
    return(GenomicFeatures:::makeTxDbFromGRanges(rtracklayer::import(gff_path)))
  }
  stop("无法构建 TxDb：请安装 txdbmaker (Bioc 3.20+) 或更新 GenomicFeatures")
}

# ── 缓存策略：缓存 GRanges 中间产物，每次重建 TxDb ──
# 为什么不缓存 TxDb 对象？
#   saveRDS(txdb) → readRDS() 后 TxDb 底层 Rcpp 指针失效，报 bad_weak_ptr。
# 缓存 GRanges 是稳定的 base 对象，从 GRanges → TxDb 只需 1-2 秒
# （远快于重新解析 GFF 的 10-20 秒）。
cache_file <- paste0(gff_path, ".gff_gr_cache.rds")
gff_mtime  <- file.info(gff_path)$mtime

makeTxDbFromGffPath <- function(gff_path) {
  # 统一入口：优先 txdbmaker::makeTxDbFromGRanges，次选旧版 API
  gr <- rtracklayer::import(gff_path)
  # rtracklayer::import 结果必须持有到 TxDb 构建完成（避免 GC 回收底层指针）
  on.exit(rm(gr), add = TRUE)
  gr_hold <- gr  # 显式保持引用
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdb <- txdbmaker::makeTxDbFromGRanges(gr)
  } else if (exists("makeTxDbFromGRanges", where = "package:GenomicFeatures", mode = "function")) {
    txdb <- GenomicFeatures::makeTxDbFromGRanges(gr)
  } else if (exists("makeTxDbFromGFF", where = "package:GenomicFeatures", mode = "function")) {
    txdb <- GenomicFeatures::makeTxDbFromGFF(gff_path)
  } else {
    txdb <- GenomicFeatures:::makeTxDbFromGRanges(gr)
  }
  txdb
}

if (file.exists(cache_file) && !is.na(gff_mtime)) {
  cache_info <- tryCatch(readRDS(cache_file), error = function(e) NULL)
  if (!is.null(cache_info) && identical(cache_info$mtime, gff_mtime)) {
    cat("  GRanges loaded from cache, building TxDb...\n")
    gr_cached <- cache_info$gr
    txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) {
              txdbmaker::makeTxDbFromGRanges(gr_cached)
            } else {
              GenomicFeatures::makeTxDbFromGRanges(gr_cached)
            }
    cat("  TxDb ready (from cached GRanges):", basename(cache_file), "\n\n")
  } else {
    cat("  GRanges cache stale/broken, rebuilding TxDb from GFF...\n")
    txdb <- makeTxDbFromGffPath(gff_path)
    gr_new <- rtracklayer::import(gff_path)
    saveRDS(list(mtime = gff_mtime, gr = gr_new), cache_file)
    cat("  GRanges cached to:", basename(cache_file), "\n\n")
  }
} else {
  cat("  GRanges cache miss, building TxDb from GFF...\n")
  txdb <- makeTxDbFromGffPath(gff_path)
  if (!is.na(gff_mtime)) {
    gr_new <- rtracklayer::import(gff_path)
    saveRDS(list(mtime = gff_mtime, gr = gr_new), cache_file)
    cat("  GRanges cached to:", basename(cache_file), "\n\n")
  } else {
    cat("  TxDb created (no cache, gff mtime unknown)\n\n")
  }
}

# ── Annotate peaks ────────────────────────────────────────────────────────
cat("Annotating peaks...\n")

# 检查文件是否为空
if (!file.exists(peak_file) || file.info(peak_file)$size == 0) {
    cat("  SKIP: peak file missing or empty (0 peaks).\n")
    quit(status = 0, save = "no")
}

# 读 peak 文件 —— 双重策略：ChIPseeker::readPeakFile → 失败则手动构造 GRanges
# readPeakFile 有时会因列数/NA 问题报 "start/end cannot contain NAs"
peaks <- tryCatch(
    readPeakFile(peak_file, header = FALSE),
    error = function(e) {
        cat("  [WARN] readPeakFile failed:", e$message, "\n")
        cat("         Falling back to manual GRanges construction...\n")
        # 手动读，过滤 NA chr/start/end 行
        dt <- data.table::fread(peak_file, header = FALSE,
                                colClasses = list(character = 1L))
        # narrowPeak/BED 都是前三列 chr/start/end
        ok <- !is.na(dt$V1) & !is.na(dt$V2) & !is.na(dt$V3) &
              nzchar(as.character(dt$V1))
        if (sum(ok) == 0) {
            cat("  SKIP: no valid peaks after NA filtering.\n")
            quit(status = 0, save = "no")
        }
        if (!all(ok)) {
            cat("  Filtered out", sum(!ok), "rows with NA chr/start/end\n")
        }
        gr <- GenomicRanges::GRanges(
            seqnames = dt$V1[ok],
            ranges   = IRanges::IRanges(start = as.integer(dt$V2[ok]) + 1L,
                                         end   = as.integer(dt$V3[ok]))
        )
        cat("  Manually loaded", length(gr), "peaks as GRanges\n")
        gr
    }
)

# 兼容 data.frame 和 GRanges (ChIPseeker 版本差异)
if (is.data.frame(peaks)) {
    n_peaks <- nrow(peaks)
} else if (methods::is(peaks, "GRanges")) {
    n_peaks <- length(peaks)
} else {
    n_peaks <- tryCatch(nrow(peaks), error = function(e) NA)
}
if (is.na(n_peaks) || n_peaks == 0) {
    cat("  SKIP: no valid peaks (class=", class(peaks)[1],
        ", n_peaks=", n_peaks, ").\n", sep = "")
    quit(status = 0, save = "no")
}
cat("  Peaks to annotate:", n_peaks, "\n")

# 只有 data.frame 需要手动赋 colnames (GRanges 已有标准列名)
if (is.data.frame(peaks)) {
    peak_colnames <- if (grepl("\\.narrowPeak$", peak_file)) {
        c("chr", "start", "end", "name", "score", "strand",
          "signalValue", "pValue", "qValue", "peak")
    } else {
        c("chr", "start", "end")
    }
    # 只赋前面几列 (如果 peaks 有更多列就保留原名)
    ncn <- min(length(peak_colnames), ncol(peaks))
    colnames(peaks)[1:ncn] <- peak_colnames[1:ncn]
}

annotation <- annotatePeak(
    peaks,
    TxDb          = txdb,
    tssRegion     = c(-tss_distance, 0),
    addFlankGeneInfo = TRUE,
    flankDistance = tss_distance,
    sameStrand    = FALSE,
    ignoreOverlap = FALSE
)

# ── Write annotation table ───────────────────────────────────────────────
anno_df <- as.data.frame(annotation)
out_txt <- paste0(peak_name, "_annotation.txt")
write.table(anno_df, file = out_txt, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
cat("  ->", out_txt, "\n")

# ── Annotation pie chart ──────────────────────────────────────────────────
cat("\nGenerating annotation pie chart...\n")
out_pie <- paste0(peak_name, "_annotation_pie.pdf")
pdf(out_pie, width = 10, height = 8)
tryCatch(
    plotAnnoPie(annotation),
    error = function(e) cat("  Warning: could not generate pie chart:", e$message, "\n")
)
dev.off()
cat("  ->", out_pie, "\n")

# ── Extract gene IDs ──────────────────────────────────────────────────────
gids <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
cat("\n  Genes annotated:", length(gids), "\n")

if (length(gids) == 0) {
    cat("  No genes found in annotation, skipping GO enrichment.\n")
    q()
}

# ── GO Enrichment Analysis ────────────────────────────────────────────────
# Static OrgDb mapping (genomes whose GFF IDs match official Bioconductor OrgDb)
org_db_map <- list(
    ath = "org.At.tair.db"
)

# Species scientific names for AnnotationHub fallback (Ensembl Plants species)
species_name_map <- list(
    osa      = "Oryza sativa",
    b73      = "Zea mays",
    zma      = "Zea mays",
    w22      = "Zea mays",
    gma      = "Glycine max",
    zh13     = "Glycine max",
    hx3      = "Glycine max",
    w82njau  = "Glycine max",
    bra      = "Brassica rapa",
    ahy      = "Arachis hypogaea",
    cla      = "Citrullus lanatus",
    cme      = "Cucumis melo",
    csa      = "Cucumis sativus",
    mpo      = "Marchantia polymorpha",
    sly      = "Solanum lycopersicum",
    tae      = "Triticum aestivum"
)

# Attempt to resolve OrgDb
org_obj <- NULL
org_key <- "ENTREZID"

# Step 1: Static mapping (official pre-installed OrgDb packages)
static_pkg <- org_db_map[[genome]]
if (!is.null(static_pkg)) {
    cat("\nLoading organism database:", static_pkg, "...\n")
    if (suppressPackageStartupMessages(require(static_pkg, character.only = TRUE))) {
        org_obj <- get(static_pkg)
        if (genome == "ath") org_key <- "TAIR"
    } else {
        cat("  Warning: package", static_pkg, "not installed.\n")
    }
}

# Step 2: AnnotationHub fallback (for species in Ensembl Plants)
if (is.null(org_obj)) {
    species <- species_name_map[[genome]]
    if (!is.null(species)) {
        cat("\nSearching AnnotationHub for '", species, "'...\n", sep = "")
        tryCatch({
            suppressPackageStartupMessages(library(AnnotationHub))
            ah <- AnnotationHub()
            ah_query <- query(ah, c("OrgDb", species))
            if (length(ah_query) > 0) {
                org_obj <- ah_query[[length(ah_query)]]
                cat("  -> Retrieved OrgDb via AnnotationHub for", species, "\n")
            } else {
                cat("  -> No OrgDb found for", species, "in AnnotationHub.\n")
            }
        }, error = function(e) {
            cat("  -> AnnotationHub query failed:", e$message, "\n")
        })
    } else {
        cat("\n  No species name mapping for genome '", genome, "'.\n", sep = "")
    }
}

# Step 3: Skip if no OrgDb available
if (is.null(org_obj)) {
    cat("\nNo OrgDb available for genome '", genome, "', skipping GO enrichment.\n", sep = "")
}

# ── Run enrichGO ─────────────────────────────────────────────────────────
if (!is.null(org_obj)) {
    cat("  Running enrichGO for:", peak_name, "(", length(gids), "genes )...\n")
    cat("  Using keyType:", org_key, "\n")
    ego <- tryCatch(
        enrichGO(
            gene          = gids,
            OrgDb         = org_obj,
            keyType       = org_key,
            ont           = "ALL",
            pvalueCutoff  = 0.05,
            qvalueCutoff  = 0.1,
            readable      = TRUE
        ),
        error = function(e) {
            cat("    enrichGO failed:", e$message, "\n")
            return(NULL)
        }
    )
    if (!is.null(ego) && nrow(ego@result) > 0) {
        out_go <- paste0(peak_name, "_go_enrichment.txt")
        write.table(as.data.frame(ego), file = out_go,
                    sep = "\t", quote = FALSE, row.names = FALSE)
        cat("    ->", out_go, "\n")

        # GO dotplot
        out_dot <- paste0(peak_name, "_go_dotplot.pdf")
        cat("\nGenerating GO dotplot...\n")
        pdf(out_dot, width = 12, height = 8)
        tryCatch(
            print(dotplot(ego, showCategory = 10) +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                  ggtitle(paste("GO Enrichment —", peak_name))),
            error = function(e) cat("  Warning: could not generate dotplot:", e$message, "\n")
        )
        dev.off()
        cat("  ->", out_dot, "\n")
    } else {
        cat("    No significant GO terms found.\n")
    }
}

cat("\nChIPseeker annotation and GO enrichment completed for", peak_name, ".\n")
