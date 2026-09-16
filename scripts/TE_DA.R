#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# TE_DA.R — Separate DESeq2 (方法 B), 多 treatment vs 同一 control
#
# 原理:
#   control = 用户在 -rc/-bc 里写的分组名 (第一个 group)
#   对每个 treatment 组:
#     1. 取 control + 该 treatment 的 ribo/rna count 子集
#     2. 分别跑 DESeq2 ribo DE 和 rna DE
#     3. delta_logFC = logFC_ribo - logFC_rna
#     4. 输出: TE_{ctrl}vs_{trt}_{all|sig|summary|volcano}.{tsv|pdf}
#
# 输入:
#   args[1]: ribo_count_matrix.tsv
#   args[2]: rna_count_matrix.tsv
#   args[3]: colData.tsv
#   args[4]: output_dir
#   args[5]: alpha
#   args[6]: control_group  (关键! 用户给的 control 组名, 空则自动识别)
# ─────────────────────────────────────────────────────────────────────────────

options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)

ribo_file    <- args[1]
rna_file     <- args[2]
coldata_file <- args[3]
output_dir   <- args[4]
alpha        <- as.numeric(if (length(args) >= 5) args[5] else "0.05")
control_name <- if (length(args) >= 6 && nzchar(args[6])) args[6] else NULL
te_dir_arg <- if (length(args) >= 7 && nzchar(args[7])) args[7] else NULL
te_combined_file <- if (!is.null(te_dir_arg)) file.path(te_dir_arg, "te_combined.tsv") else NULL

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))

message("═══════════════════════════════════════════")
message("  TE_DA: Separate DESeq2 (多 treatment × 1 control)")
message("═══════════════════════════════════════════")

# ── 读取 ─────────────────────────────────────────────────────────────────────
# 读 count matrix — 防御性 aggregate 重复 gene_id
.rd1 <- read.delim(ribo_file, check.names = FALSE)
.rd2 <- read.delim(rna_file,  check.names = FALSE)
.rd3 <- read.delim(coldata_file, check.names = FALSE,
                   stringsAsFactors = TRUE)

rownames(.rd1) <- .rd1[, 1]; .rd1 <- .rd1[, -1, drop = FALSE]
rownames(.rd2) <- .rd2[, 1]; .rd2 <- .rd2[, -1, drop = FALSE]
rownames(.rd3) <- .rd3[, 1]; .rd3 <- .rd3[, -1, drop = FALSE]

# 防御: 重复 gene_id 会让 DESeq2 / merge 报错 — aggregate 合并
.fix_dup <- function(mat, label) {
  dups <- duplicated(rownames(mat))
  if (any(dups)) {
    message("  WARNING: ", label, " has ", sum(dups), " duplicate gene_ids — summing counts")
    agg <- aggregate(. ~ rownames(mat), data = as.data.frame(mat), FUN = sum)
    rownames(agg) <- agg[, 1]; agg <- agg[, -1, drop = FALSE]
    mat <- as.matrix(agg)
  }
  mat
}
ribo_mat <- .fix_dup(.rd1, "ribo_counts.tsv")
rna_mat  <- .fix_dup(.rd2, "rna_counts.tsv")
coldata  <- as.data.frame(.rd3)

common_genes <- intersect(rownames(ribo_mat), rownames(rna_mat))
ribo_mat <- ribo_mat[common_genes, , drop = FALSE]
rna_mat  <- rna_mat[common_genes, , drop = FALSE]

conditions <- levels(coldata$condition)
message("所有 condition: ", paste(conditions, collapse = ", "))

# ── 确定 control ─────────────────────────────────────────────────────────────
if (is.null(control_name)) {
  # 没给 → 尝试匹配常见名
  patterns <- c("ctrl", "control", "wt", "wild", "normal", "untreated",
                "mock", "vehicle", "dmso")
  for (pat in patterns) {
    hit <- grep(pat, conditions, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) { control_name <- hit[1]; break }
  }
  if (is.null(control_name)) control_name <- conditions[1]
}
message("Control 组: ", control_name)

treatments <- setdiff(conditions, control_name)
message("Treatment 组: ", paste(treatments, collapse = ", "))

if (length(treatments) == 0) {
  stop("只有一个 condition, 无法比较. 需要至少 2 组.")
}

# ── 辅助函数: 跑一次 pairwise TE 比较 ────────────────────────────────────────
run_pairwise <- function(trt_name, ribo_mat, rna_mat, coldata,
                         control_name, alpha, out_dir) {
  pair_tag <- paste0(gsub("[^A-Za-z0-9_]", "_", control_name),
                     "_vs_",
                     gsub("[^A-Za-z0-9_]", "_", trt_name))
  message("\n─────────────────────────────────────────────")
  message("比较: ", control_name, " vs ", trt_name, "  (tag=", pair_tag, ")")

  # 只保留这两个 condition 的样本
  keep_samples <- rownames(coldata)[coldata$condition %in% c(control_name, trt_name)]
  sub_coldata  <- coldata[keep_samples, , drop = FALSE]

  # 拆分 ribo/rna
  ribo_samples <- rownames(sub_coldata)[sub_coldata$data_type == "ribo"]
  rna_samples  <- rownames(sub_coldata)[sub_coldata$data_type == "rna"]

  sub_ribo_mat <- ribo_mat[, ribo_samples, drop = FALSE]
  sub_rna_mat  <- rna_mat[,  rna_samples,  drop = FALSE]

  sub_ribo_cd  <- sub_coldata[ribo_samples, , drop = FALSE]
  sub_rna_cd   <- sub_coldata[rna_samples,  , drop = FALSE]

  # relevel: control 为 ref
  sub_ribo_cd$condition <- relevel(sub_ribo_cd$condition, ref = control_name)
  sub_rna_cd$condition  <- relevel(sub_rna_cd$condition,  ref = control_name)

  message("  Ribo samples (", ncol(sub_ribo_mat), "): ", paste(ribo_samples, collapse = ", "))
  message("  RNA  samples (", ncol(sub_rna_mat),  "): ", paste(rna_samples,  collapse = ", "))

  if (ncol(sub_ribo_mat) < 2 || ncol(sub_rna_mat) < 2) {
    message("  SKIP: 每组至少需要 2 个样本 (每个 data_type)")
    return(NULL)
  }

  # ── DESeq2 Ribo-seq ──
  tryCatch({
    dds_ribo <- DESeqDataSetFromMatrix(
      countData = round(sub_ribo_mat),
      colData   = sub_ribo_cd,
      design    = ~ condition
    )
    keep_r <- rowSums(counts(dds_ribo) >= 10) >= ncol(dds_ribo) / 2
    dds_ribo <- dds_ribo[keep_r, ]
    dds_ribo <- DESeq(dds_ribo, quiet = TRUE)
    res_ribo <- results(dds_ribo)
  }, error = function(e) {
    message("  DESeq2 Ribo error: ", e$message)
    return(NULL)
  })

  # ── DESeq2 RNA-seq ──
  tryCatch({
    dds_rna <- DESeqDataSetFromMatrix(
      countData = round(sub_rna_mat),
      colData   = sub_rna_cd,
      design    = ~ condition
    )
    keep_r <- rowSums(counts(dds_rna) >= 10) >= ncol(dds_rna) / 2
    dds_rna <- dds_rna[keep_r, ]
    dds_rna <- DESeq(dds_rna, quiet = TRUE)
    res_rna <- results(dds_rna)
  }, error = function(e) {
    message("  DESeq2 RNA error: ", e$message)
    return(NULL)
  })

  if (is.null(res_ribo) || is.null(res_rna)) return(NULL)

  # ── 合并 ──
  df_ribo <- as.data.frame(res_ribo)
  df_rna  <- as.data.frame(res_rna)
  colnames(df_ribo) <- paste0("ribo_", colnames(df_ribo))
  colnames(df_rna)  <- paste0("rna_",  colnames(df_rna))

  merged <- merge(df_ribo, df_rna, by = "row.names", all = TRUE)
  colnames(merged)[1] <- "gene_id"
  merged$delta_logFC <- merged$ribo_log2FoldChange - merged$rna_log2FoldChange
  merged$comparison  <- pair_tag

  # ── 分类 ──
  merged$category <- "NONE"
  sig_ribo <- !is.na(merged$ribo_padj) & merged$ribo_padj < alpha
  sig_rna  <- !is.na(merged$rna_padj)  & merged$rna_padj  < alpha
  merged$category[ sig_ribo & !sig_rna & merged$delta_logFC > 0]  <- "TRANSLATION_UP"
  merged$category[ sig_ribo & !sig_rna & merged$delta_logFC < 0]  <- "TRANSLATION_DOWN"
  merged$category[!sig_ribo &  sig_rna]                            <- "TRANSCRIPTION_ONLY"
  merged$category[ sig_ribo &  sig_rna &
                   sign(merged$ribo_log2FoldChange) == sign(merged$rna_log2FoldChange)] <- "COORDINATED"
  merged$category[ sig_ribo &  sig_rna &
                   sign(merged$ribo_log2FoldChange) != sign(merged$rna_log2FoldChange)] <- "RIBO_ONLY"

  merged <- merged[order(merged$delta_logFC), ]

  # ── 附加每个基因在每个样本 pair 中的 TE 值 ──
  # te_combined.tsv 由 STEP 11 的 _generate_te_summary 自动生成
  # 列名格式: {ribo_tag}__{rna_tag}_te, {ribo_tag}__{rna_tag}_log2te, ...
  if (!is.null(te_combined_file) && file.exists(te_combined_file)) {
    tryCatch({
      tc <- read.delim(te_combined_file, check.names = FALSE, stringsAsFactors = FALSE)
      # 只保留 gene_id + _te / _log2te 后缀的列 (排除 _ribo_tpm, _rna_tpm 等)
      all_cols <- colnames(tc)
      keep_cols <- c("gene_id", grep("_(te|log2te)$", all_cols, value = TRUE))
      # 安全: 确保 keep_cols 都存在
      keep_cols <- keep_cols[keep_cols %in% all_cols]
      if (length(keep_cols) > 1) {
        tc_sub <- tc[, keep_cols, drop = FALSE]
        merged <- merge(merged, tc_sub, by = "gene_id", all.x = TRUE)
        message("  Merged ", length(keep_cols) - 1, " TE/log2TE columns from te_combined.tsv")
      }
    }, error = function(e) message("  WARNING: TE merge skipped: ", e))
  }

  # ── 输出 ──
  write.table(merged, file.path(out_dir, paste0("TE_", pair_tag, "_all.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  sig_df <- merged[merged$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN"), ]
  write.table(sig_df, file.path(out_dir, paste0("TE_", pair_tag, "_sig.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  summary_df <- as.data.frame(table(merged$category), stringsAsFactors = FALSE)
  colnames(summary_df) <- c("category", "n_genes")
  write.table(summary_df, file.path(out_dir, paste0("TE_", pair_tag, "_summary.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  message("  总计 ", nrow(merged), " genes | ",
          sum(sig_df$category == "TRANSLATION_UP"), " TE up | ",
          sum(sig_df$category == "TRANSLATION_DOWN"), " TE down")

  # ── 火山图 ──
  tryCatch({
    combined_p <- -log10(pmax(merged$ribo_padj, merged$rna_padj, na.rm = TRUE))
    pdf(file.path(out_dir, paste0("TE_", pair_tag, "_volcano.pdf")), 7, 6)
    plot(merged$delta_logFC, combined_p, pch = 16, cex = 0.7, col = "grey70",
         xlab = expression(Delta*log[2]("FC"[ribo] - "FC"[rna])),
         ylab = expression(-log[10]("padj")),
         main = paste0("TE: ", control_name, " vs ", trt_name))
    sig_idx <- merged$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN")
    points(merged$delta_logFC[sig_idx], combined_p[sig_idx],
           pch = 16, cex = 1,
           col = ifelse(merged$delta_logFC[sig_idx] > 0, "red", "blue"))
    abline(v = 0, lty = 2, col = "grey30")
    abline(h = -log10(alpha), lty = 3, col = "grey50")
    legend("topright",
           legend = c(paste0("TE up (", sum(merged$category == "TRANSLATION_UP"), ")"),
                      paste0("TE down (", sum(merged$category == "TRANSLATION_DOWN"), ")")),
           col = c("red", "blue"), pch = 16)
    dev.off()
  }, error = function(e) message("  volcano plot skipped: ", e$message))

  return(merged)
}

# ── 循环: control vs 每个 treatment ──────────────────────────────────────────
all_results <- list()
for (trt in treatments) {
  res <- run_pairwise(trt, ribo_mat, rna_mat, coldata,
                      control_name, alpha, output_dir)
  if (!is.null(res)) all_results[[trt]] <- res
}

# ── 合并 summary ──
if (length(all_results) > 1) {
  combined <- do.call(rbind, all_results)
  write.table(combined, file.path(output_dir, "TE_separate_all_combined.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  message("\n═══ Combined all pairwise results ═══ ", nrow(combined), " genes")
}

message("\n═══ TE_DA complete (", length(all_results), " comparisons) ═══ ")
