#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# TE_DS.R — Joint DESeq2 (方法 C), 多 treatment × 1 control
#
# 每次只取 control + 一个 treatment 的样本做联合双因素 DESeq2:
#   design = ~ condition + data_type + condition:data_type
#   交互项 = delta_logFC (TE 变化)
#
# 输入:
#   args[1]: ribo_count_matrix.tsv
#   args[2]: rna_count_matrix.tsv
#   args[3]: colData.tsv
#   args[4]: output_dir
#   args[5]: alpha
#   args[6]: control_group
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
message("  TE_DS: Joint DESeq2 (多 treatment × 1 control)")
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
if (length(treatments) == 0) stop("只有一个 condition, 无法比较")

# ── 辅助函数 ─────────────────────────────────────────────────────────────────
run_pairwise_joint <- function(trt_name, ribo_mat, rna_mat, coldata,
                               control_name, alpha, out_dir) {
  pair_tag <- paste0(gsub("[^A-Za-z0-9_]", "_", control_name),
                     "_vs_",
                     gsub("[^A-Za-z0-9_]", "_", trt_name))
  message("\n─────────────────────────────────────────────")
  message("比较: ", control_name, " vs ", trt_name)

  # 只取 control + trt
  keep_samples <- rownames(coldata)[coldata$condition %in% c(control_name, trt_name)]
  sub_coldata  <- coldata[keep_samples, , drop = FALSE]

  ribo_samples <- rownames(sub_coldata)[sub_coldata$data_type == "ribo"]
  rna_samples  <- rownames(sub_coldata)[sub_coldata$data_type == "rna"]

  if (length(ribo_samples) < 2 || length(rna_samples) < 2) {
    message("  SKIP: 每组 (ribo/rna) 至少需要 2 个样本")
    return(NULL)
  }

  message("  Ribo (", length(ribo_samples), "): ", paste(ribo_samples, collapse = ", "))
  message("  RNA  (", length(rna_samples),  "): ", paste(rna_samples,  collapse = ", "))

  # 合并 count matrix
  sub_ribo <- ribo_mat[, ribo_samples, drop = FALSE]
  sub_rna  <- rna_mat[,  rna_samples,  drop = FALSE]
  combined_counts <- cbind(round(sub_ribo), round(sub_rna))
  # 按 sub_coldata 顺序
  combined_counts <- combined_counts[, rownames(sub_coldata), drop = FALSE]

  # relevel
  sub_coldata$condition <- relevel(sub_coldata$condition, ref = control_name)
  sub_coldata$data_type <- relevel(sub_coldata$data_type, ref = "rna")

  # ── Joint DESeq2 ──
  tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = combined_counts,
      colData   = sub_coldata,
      design    = ~ condition + data_type + condition:data_type
    )
    keep <- rowSums(counts(dds) >= 10) >= ncol(dds) / 2
    dds <- dds[keep, ]
    dds <- DESeq(dds, quiet = TRUE)

    rn <- resultsNames(dds)
    message("  resultsNames: ", paste(rn, collapse = ", "))

    # 交互项
    int_name <- rn[grepl(":", rn) | grepl("condition.*data_type|data_type.*condition", rn)]
    if (length(int_name) == 0) int_name <- rn[length(rn)]

    res_int   <- results(dds, name = int_name)
    int_label <- int_name

    # condition 主效应
    cond_name <- rn[grepl("^condition", rn) & !grepl(":", rn)][1]
    res_cond  <- if (!is.na(cond_name)) results(dds, name = cond_name) else NULL

  }, error = function(e) {
    message("  DESeq2 error: ", e$message)
    return(NULL)
  })

  if (is.null(res_int)) return(NULL)

  # ── 合并 ──
  df_int <- as.data.frame(res_int)
  colnames(df_int) <- paste0("TE_int_", colnames(df_int))
  df_int$delta_logFC <- df_int$TE_int_log2FoldChange
  df_int$comparison  <- pair_tag

  if (!is.null(res_cond)) {
    df_cond <- as.data.frame(res_cond)
    colnames(df_cond) <- paste0("rna_", colnames(df_cond))
    df_int <- merge(df_int, df_cond, by = "row.names", all.x = TRUE)
    colnames(df_int)[1] <- "gene_id"
  } else {
    df_int$gene_id <- rownames(df_int)
  }

  # ── 分类 ──
  df_int$category <- "NONE"
  sig_int <- !is.na(df_int$TE_int_padj) & df_int$TE_int_padj < alpha
  sig_rna <- !is.na(df_int$rna_padj)    & df_int$rna_padj    < alpha
  df_int$category[ sig_int & df_int$delta_logFC > 0]  <- "TRANSLATION_UP"
  df_int$category[ sig_int & df_int$delta_logFC < 0]  <- "TRANSLATION_DOWN"
  df_int$category[!sig_int & sig_rna]                  <- "TRANSCRIPTION_ONLY"
  df_int <- df_int[order(df_int$delta_logFC), ]

  # ── 附加每个基因在每个样本 pair 中的 TE 值 ──
  if (!is.null(te_combined_file) && file.exists(te_combined_file)) {
    tryCatch({
      tc <- read.delim(te_combined_file, check.names = FALSE, stringsAsFactors = FALSE)
      all_cols <- colnames(tc)
      keep_cols <- c("gene_id", grep("_(te|log2te)$", all_cols, value = TRUE))
      keep_cols <- keep_cols[keep_cols %in% all_cols]
      if (length(keep_cols) > 1) {
        tc_sub <- tc[, keep_cols, drop = FALSE]
        df_int <- merge(df_int, tc_sub, by = "gene_id", all.x = TRUE)
        message("  Merged ", length(keep_cols) - 1, " TE/log2TE columns from te_combined.tsv")
      }
    }, error = function(e) message("  WARNING: TE merge skipped: ", e))
  }

  # ── 输出 ──
  write.table(df_int, file.path(out_dir, paste0("TE_", pair_tag, "_all.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  sig_df <- df_int[df_int$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN"), ]
  write.table(sig_df, file.path(out_dir, paste0("TE_", pair_tag, "_sig.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  summary_df <- as.data.frame(table(df_int$category), stringsAsFactors = FALSE)
  colnames(summary_df) <- c("category", "n_genes")
  write.table(summary_df, file.path(out_dir, paste0("TE_", pair_tag, "_summary.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  message("  总计 ", nrow(df_int), " genes | ",
          sum(sig_df$category == "TRANSLATION_UP"), " TE up | ",
          sum(sig_df$category == "TRANSLATION_DOWN"), " TE down")

  # ── 火山图 ──
  tryCatch({
    pdf(file.path(out_dir, paste0("TE_", pair_tag, "_volcano.pdf")), 7, 6)
    plot(df_int$delta_logFC, -log10(df_int$TE_int_padj),
         pch = 16, cex = 0.7, col = "grey70",
         xlab = expression(Delta*log[2]("FC"[TE])),
         ylab = expression(-log[10]("padj"[interaction])),
         main = paste0("TE (joint): ", control_name, " vs ", trt_name))
    sig_idx <- df_int$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN")
    points(df_int$delta_logFC[sig_idx], -log10(df_int$TE_int_padj[sig_idx]),
           pch = 16, cex = 1,
           col = ifelse(df_int$delta_logFC[sig_idx] > 0, "red", "blue"))
    abline(v = 0, lty = 2, col = "grey30")
    abline(h = -log10(alpha), lty = 3, col = "grey50")
    legend("topright",
           legend = c(paste0("TE up (", sum(df_int$category == "TRANSLATION_UP"), ")"),
                      paste0("TE down (", sum(df_int$category == "TRANSLATION_DOWN"), ")")),
           col = c("red", "blue"), pch = 16)
    dev.off()
  }, error = function(e) message("  volcano plot skipped"))

  # ── MA plot ──
  tryCatch({
    mean_expr <- rowMeans(counts(dds, normalized = TRUE))
    pdf(file.path(out_dir, paste0("TE_", pair_tag, "_MAPlot.pdf")), 7, 6)
    plot(log2(mean_expr + 1), df_int$delta_logFC,
         pch = 16, cex = 0.7, col = "grey70",
         xlab = expression(log[2]("mean normalized count + 1")),
         ylab = expression(Delta*log[2]("FC"[TE])),
         main = paste0("MA (joint): ", control_name, " vs ", trt_name))
    sig_idx <- df_int$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN")
    points(log2(mean_expr + 1)[sig_idx], df_int$delta_logFC[sig_idx],
           pch = 16, cex = 1,
           col = ifelse(df_int$delta_logFC[sig_idx] > 0, "red", "blue"))
    abline(h = 0, lty = 2, col = "grey30")
    dev.off()
  }, error = function(e) message("  MA plot skipped"))

  return(df_int)
}

# ── 循环 ──
all_results <- list()
for (trt in treatments) {
  res <- run_pairwise_joint(trt, ribo_mat, rna_mat, coldata,
                            control_name, alpha, output_dir)
  if (!is.null(res)) all_results[[trt]] <- res
}

if (length(all_results) > 1) {
  combined <- do.call(rbind, all_results)
  write.table(combined, file.path(output_dir, "TE_joint_all_combined.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

message("\n═══ TE_DS complete (", length(all_results), " comparisons) ═══ ")
