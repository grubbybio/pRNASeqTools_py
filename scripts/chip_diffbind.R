#!/usr/bin/env Rscript
# ==============================================================================
# chip_diffbind.R — DiffBind differential binding analysis (manual count mode)
#
# DiffBind 3.x 的 dba.count() 会把 Input reads 合并到对应 IP sample,
# dual_factor 模型 (~ Condition + Factor + Condition:Factor) 无法估计
# IP vs Input 交互项. 因此我们:
#
#   1. 让 DiffBind 只生成 consensus peaks (dba + dba.overlap, 不 count)
#   2. tf.py 里用 bedtools multicov 从所有 BAM (IP + Input) 手动统计
#      每个 peak 的 reads → peak_counts.tsv
#   3. 本脚本读 peak_counts.tsv + samplesheet.csv 构造 DESeq2 模型
#
# Args:
#   [1] samplesheet.csv   含所有 IP 样本 (SampleID,Condition,Factor,Replicate,
#                         bamReads,bamControl,Peaks)
#   [2] out_dir
#   [3] genome
#   [4] analysis          affinity | dual_factor
#   [5] norm              deseq2 | total | mito | chloro
#   [6] pvalue
#   [7] foldchange
#   [8] norm_factors.tsv  (optional, for non-deseq2)
#   [9] peak_counts.tsv   (bedtools multicov 输出, 第 7+ 列是各样本 counts,
#                         列名 = 样本标签, 按 IP1, Input1, IP2, Input2... 顺序)
#
# Prerequisites: DiffBind, DESeq2
# ==============================================================================

suppressPackageStartupMessages({
  library(DiffBind)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 7)

ss_file     <- args[1]
out_dir     <- args[2]
genome      <- args[3]
analysis    <- args[4]
norm_mode   <- args[5]
pval        <- as.numeric(args[6])
fc          <- as.numeric(args[7])
nf_file     <- if (length(args) >= 8 && nzchar(args[8])) args[8] else NULL
counts_file <- if (length(args) >= 9 && nzchar(args[9])) args[9] else NULL
meta_file   <- if (length(args) >= 10 && nzchar(args[10])) args[10] else NULL

if (is.null(counts_file) || !file.exists(counts_file)) {
  counts_file <- file.path(out_dir, "peak_counts.tsv")
}
if (!file.exists(counts_file)) {
  stop("peak_counts.tsv not found: ", counts_file)
}
if (is.null(meta_file) || !file.exists(meta_file)) {
  meta_file <- file.path(out_dir, "sample_metadata.tsv")
}
if (!file.exists(meta_file)) {
  stop("sample_metadata.tsv not found: ", meta_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("═══ DiffBind Analysis (manual counts) ═══\n")
cat("  Sample sheet :", ss_file, "\n")
cat("  Counts file  :", counts_file, "\n")
cat("  Output dir   :", out_dir, "\n")
cat("  Analysis     :", analysis, "\n")
cat("  Norm mode    :", norm_mode, "\n")
cat("  P-value      :", pval, "\n")
cat("  Fold-change  :", fc, "\n")
if (!is.null(nf_file)) cat("  Norm factors :", nf_file, "\n")

# ── 1. Read samplesheet (IP samples only; Input 标签在 counts.tsv 头部) ─────
samples <- read.csv(ss_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("\n  Samplesheet (IP samples):", nrow(samples), "\n")
cat("  Conditions:", paste(unique(samples$Condition), collapse = ", "), "\n")

# ── 2. Read peak_counts.tsv ───────────────────────────────────────────────────
# 结构: chr, start, end, name, score, strand, {IP1, Input1, IP2, Input2, ...}
counts_df <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
cat("  Peak counts:", nrow(counts_df), "peaks x", ncol(counts_df) - 6,
    "samples (cols 7+)\n")

# 分离 peak 坐标 和 count matrix
peak_coords <- counts_df[, 1:6]
counts_mat <- as.matrix(counts_df[, 7:ncol(counts_df), drop = FALSE])
storage.mode(counts_mat) <- "integer"   # DESeq2 要求
sample_labels <- colnames(counts_mat)
cat("  Sample labels:", paste(sample_labels, collapse = ", "), "\n")

# ── 3. 读 sample_metadata.tsv (由 tf.py Step 3.8e 生成, 不用 regex 猜) ───────
meta_full <- read.delim(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("  Full sample metadata (IP + Input):\n")
print(meta_full)

# 按 sample_labels 顺序提取 condition/factor
rownames(meta_full) <- meta_full$sample
condition_vec <- meta_full[sample_labels, "condition"]
factor_vec    <- meta_full[sample_labels, "factor"]
stopifnot(!any(is.na(condition_vec)), !any(is.na(factor_vec)))

cat("\n  Sample meta table:\n")
print(table(condition_vec, factor_vec))

# ── 4. Normalization ──────────────────────────────────────────────────────────
cat("\n  ── Normalization:", norm_mode, "──\n")

custom_size_factors <- NULL  # 全样本的 size factor 向量 (IP + Input)
if (norm_mode == "deseq2") {
  cat("  Using DESeq2 default size factors (median-of-ratios)\n")
} else if (norm_mode %in% c("total", "mito", "chloro")) {
  if (is.null(nf_file) || !file.exists(nf_file)) {
    cat("  WARNING: norm_factors.tsv not found, falling back to DESeq2 default\n")
  } else {
    nf <- read.delim(nf_file, stringsAsFactors = FALSE)
    cat("  Norm factors loaded:", nrow(nf), "IP samples\n")
    # nf$ip_scale 是按 IP sample 顺序排的
    # 需要按 sample_labels 的顺序扩展: IP → nf$ip_scale, Input → 1.0
    if (!("ip_scale" %in% colnames(nf)) || !("sample" %in% colnames(nf))) {
      cat("  WARNING: norm_factors.tsv needs 'sample' and 'ip_scale' columns\n")
    } else {
      ip_scale_map <- setNames(nf$ip_scale, nf$sample)
      custom_size_factors <- vapply(sample_labels, function(lab) {
        if (lab %in% names(ip_scale_map)) {
          # 文献 scale = ip_org / max(ip_org)
          # DESeq2 sizeFactor 也是这个 (除数语义: normalized = raw / sizeFactor)
          # 生物学含义:
          #   mito 少的样本(WT_1, scale=0.330) → sizeFactor 小 → normalized 变大 ↑
          #   mito 多的样本(ago1_27_2, scale=1.000) → sizeFactor 大 → normalized 不变
          #   → WT normalized > ago1_27 normalized ✓ (符合生物学预期)
          ip_scale_map[[lab]]
        } else {
          1.0  # Input samples
        }
      }, numeric(1))
      cat("  DESeq2 sizeFactors (= ip_scale, 除数语义):\n")
      cat("  IP:   ", paste(round(custom_size_factors[factor_vec=="IP"], 3), collapse = ", "), "\n")
      cat("  Input:", paste(round(custom_size_factors[factor_vec=="Input"], 3), collapse = ", "), "\n")
      cat("  文献方法 (乘数): ip_scale = max/ip_org\n")
      cat("  WT_1 example: mito=1732, scale=0.330 → DESeq2: raw/0.330=raw×3.03 ↑ (WT>ago1_27)\n")
    }
  }
}

# ── 5. Run DESeq2 ─────────────────────────────────────────────────────────────
cat("\n  ── DESeq2 design ──\n")
cat("  Analysis:", analysis, "\n")

if (analysis == "affinity") {
  # IP-only: ~ Condition
  ip_idx <- which(factor_vec == "IP")
  cat("  Using IP-only:", length(ip_idx), "samples\n")

  cnt <- counts_mat[, ip_idx, drop = FALSE]
  cond <- factor(condition_vec[ip_idx])
  cond <- relevel(cond, ref = as.character(cond[1]))
  keep <- rowSums(cnt >= 10) >= 2
  cat("  Peaks after filter:", sum(keep), "of", nrow(cnt), "\n")
  cnt <- cnt[keep, , drop = FALSE]
  storage.mode(cnt) <- "integer"

  dds <- DESeqDataSetFromMatrix(cnt,
                                 colData = data.frame(Condition = cond,
                                                      row.names = colnames(cnt)),
                                 design = ~ Condition)
  if (!is.null(custom_size_factors)) {
    sf_ip <- custom_size_factors[ip_idx]
    sizeFactors(dds) <- sf_ip
    cat("  Using custom size factors (IP-scaled)\n")
  }
  dds <- DESeq(dds)
  cond_levels <- levels(cond)
  res <- results(dds, contrast = c("Condition", cond_levels[2], cond_levels[1]))

} else if (analysis == "dual_factor") {
  # All samples (IP + Input): ~ Condition + Factor + Condition:Factor
  cat("  Using ALL samples:", ncol(counts_mat),
      "(", sum(factor_vec == "IP"), "IP,",
      sum(factor_vec == "Input"), "Input)\n")
  cat("  Design: ~ Condition + Factor + Condition:Factor\n")

  cnt <- counts_mat
  keep <- rowSums(cnt >= 10) >= 2
  cat("  Peaks after filter:", sum(keep), "of", nrow(cnt), "\n")
  cnt <- cnt[keep, , drop = FALSE]
  storage.mode(cnt) <- "integer"

  cond_f <- factor(make.names(condition_vec),
                   levels = make.names(unique(condition_vec)))
  # Factor: Input 在前, IP 在后 (reference = Input)
  fact_f <- factor(make.names(factor_vec),
                   levels = c(make.names("Input"), make.names("IP")))

  dds <- DESeqDataSetFromMatrix(cnt,
                                 colData = data.frame(Condition = cond_f,
                                                      Factor = fact_f,
                                                      row.names = colnames(cnt)),
                                 design = ~ Condition + Factor + Condition:Factor)
  if (!is.null(custom_size_factors)) {
    sizeFactors(dds) <- custom_size_factors
    cat("  Using custom size factors (IP-scaled, Input=1.0)\n")
  }
  dds <- DESeq(dds)

  coefs <- resultsNames(dds)
  cat("  Design coefficients:", paste(coefs, collapse = ", "), "\n")
  # DESeq2 交互项分隔符可能是 ":" 或 "."
  interact_coef <- grep("[:.]Factor", coefs, value = TRUE)
  if (length(interact_coef) == 0) {
    cat("  [WARN] No interaction term, using main effect Factor (IP vs Input)\n")
    interact_coef <- grep("^Factor", coefs, value = TRUE)[1]
  }
  cat("  Testing:", interact_coef[1], "\n")
  res <- results(dds, name = interact_coef[1])
}

# ── 6. Annotate + save ────────────────────────────────────────────────────────
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)
# peak coordinates
res_df <- cbind(peak_coords[keep, ], res_df)
res_df$category <- ifelse(res_df$padj < pval & !is.na(res_df$padj),
                           ifelse(res_df$log2FoldChange > 0, "UP", "DOWN"),
                           "NS")
if (analysis == "dual_factor") {
  res_df$comparison <- paste("Interaction:", interact_coef[1])
  out_file <- file.path(out_dir, "DiffBind_dual_factor_interaction.tsv")
  volc_pdf <- file.path(out_dir, "DiffBind_dual_factor_volcano.pdf")
} else {
  res_df$comparison <- paste(cond_levels[2], "vs", cond_levels[1])
  out_file <- file.path(out_dir, paste0("DiffBind_affinity_",
                                         cond_levels[2], "_vs_",
                                         cond_levels[1], ".tsv"))
  volc_pdf <- file.path(out_dir, "DiffBind_affinity_volcano.pdf")
  map_pdf  <- file.path(out_dir, "DiffBind_affinity_MAPlot.pdf")
}
write.table(res_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n  Results:", out_file, "\n")

# ── 7. Volcano / MA ────────────────────────────────────────────────────────────
pdf(volc_pdf, width = 7, height = 6)
plot(res_df$log2FoldChange, -log10(res_df$padj + 1e-300),
     pch = 16, cex = 0.6, col = "grey60",
     xlab = expression(log[2]~FC), ylab = expression(-log[10]~padj),
     main = if (analysis == "dual_factor") "Dual-Factor Interaction"
            else paste0("Affinity: ", cond_levels[2], " vs ", cond_levels[1]))
sig_up <- res_df$category == "UP"; sig_dn <- res_df$category == "DOWN"
points(res_df$log2FoldChange[sig_up], -log10(res_df$padj[sig_up] + 1e-300),
       pch = 16, cex = 0.7, col = "red")
points(res_df$log2FoldChange[sig_dn], -log10(res_df$padj[sig_dn] + 1e-300),
       pch = 16, cex = 0.7, col = "blue")
abline(v = c(-log2(fc), log2(fc)), h = -log10(pval), lty = 2, col = "grey30")
legend("topright",
       legend = c(paste0("Up (", sum(sig_up), ")"),
                  paste0("Down (", sum(sig_dn), ")"),
                  paste0("NS (", sum(!sig_up & !sig_dn), ")")),
       col = c("red", "blue", "grey60"), pch = 16)
dev.off()
cat("  Volcano:", basename(volc_pdf), "\n")

if (analysis == "affinity") {
  pdf(map_pdf, width = 7, height = 6)
  plotMA(res, ylim = c(-4, 4), main = paste0("MA: ", cond_levels[2], " vs ", cond_levels[1]))
  dev.off()
  cat("  MA plot:", basename(map_pdf), "\n")
}

# ── 8. Export BED files for ChIPseeker annotation ──────────────────────────────
# 三类: consensus (全部 peaks), UP (显著上调), DOWN (显著下调)
# BED format: chr, start, end  (ChIPseeker readPeakFile 可直接读)
bed_consensus <- file.path(out_dir, "DiffBind_consensus_peaks.bed")
bed_up        <- file.path(out_dir, "DiffBind_UP_peaks.bed")
bed_down      <- file.path(out_dir, "DiffBind_DOWN_peaks.bed")

write.table(res_df[, c("chr", "start", "end")],
            bed_consensus, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
write.table(res_df[res_df$category == "UP", c("chr", "start", "end")],
            bed_up, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
write.table(res_df[res_df$category == "DOWN", c("chr", "start", "end")],
            bed_down, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
cat("\n  BED for annotation:\n")
cat("    consensus :", nrow(res_df), "peaks →", basename(bed_consensus), "\n")
cat("    UP        :", sum(res_df$category == "UP"), "peaks →", basename(bed_up), "\n")
cat("    DOWN      :", sum(res_df$category == "DOWN"), "peaks →", basename(bed_down), "\n")
if (sum(res_df$category == "UP") == 0 && sum(res_df$category == "DOWN") == 0) {
  cat("  NOTE: No significant peaks (UP/DOWN); only consensus will be annotated.\n")
}

# ── 9. Summary ────────────────────────────────────────────────────────────────
cat("\n═══ Summary ═══\n")
cat("  Analysis   :", analysis, "\n")
cat("  Norm       :", norm_mode, "\n")
cat("  Design     :", if (analysis == "affinity") "~ Condition (IP-only)"
                    else "~ Condition + Factor + Condition:Factor", "\n")
cat("  Peaks kept :", nrow(cnt), "\n")
cat("  Significant:", sum(res_df$category != "NS"),
    "(UP=", sum(sig_up), "DOWN=", sum(sig_dn), ")\n")
cat("  Output dir :", out_dir, "\n")
cat("Done.\n")