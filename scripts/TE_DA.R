#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# TE_DA.R — Translational Efficiency changes via separate DESeq2 (方法 B)
#
# 原理:
#   1. 分别对 Ribo-seq count matrix 和 RNA-seq count matrix 做 DESeq2
#   2. 同一基因, treatment vs control:
#      logFC_ribo = log2(treatment_ribo / control_ribo)
#      logFC_rna  = log2(treatment_rna  / control_rna)
#      delta_logFC = logFC_ribo - logFC_rna   ← TE 变化的估计
#   3. 显著性: DESeq2 ribo DE 的 padj + DESeq2 rna DE 的 padj
#      如果两者都显著且方向一致 → TE 变化
#      (更严格: DESeq2 交互项检验, 见 TE_DS.R)
#
# 输入格式:
#   args[1]: ribo_count_matrix.tsv   (rows=genes, cols=samples, header)
#   args[2]: rna_count_matrix.tsv    (同上)
#   args[3]: colData.tsv             (sample  condition  data_type)
#   args[4]: output_dir
#   args[5]: alpha (default 0.05)
#   args[6]: lfc_threshold (default 0)
#
# colData.tsv 示例:
#   sample     condition   data_type
#   ribo_ctrl1 ctrl        ribo
#   ribo_ctrl2 ctrl        ribo
#   ribo_trt1  treatment   ribo
#   rna_ctrl1  ctrl        rna
#   rna_ctrl2  ctrl        rna
#   rna_trt1   treatment   rna
#
# 输出:
#   TE_separate_all.tsv       — 所有基因的合并结果
#   TE_separate_sig.tsv       — 显著变化的 TE 基因 (padj<alpha)
#   TE_separate_summary.tsv   — 分类汇总
#   TE_separate_volcano.pdf   — 火山图
# ─────────────────────────────────────────────────────────────────────────────

options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript TE_DA.R <ribo_count.tsv> <rna_count.tsv> ",
       "<colData.tsv> <output_dir> [alpha] [lfc_threshold]")
}

ribo_file    <- args[1]
rna_file     <- args[2]
coldata_file <- args[3]
output_dir   <- args[4]
alpha        <- as.numeric(if (length(args) >= 5) args[5] else "0.05")
lfc_thresh   <- as.numeric(if (length(args) >= 6) args[6] else "0")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))

message("═══════════════════════════════════════════")
message("  TE_DA: Separate DESeq2 (方法 B)")
message("═══════════════════════════════════════════")

# ── 读取数据 ─────────────────────────────────────────────────────────────────
ribo_mat <- as.matrix(read.delim(ribo_file, row.names = 1, check.names = FALSE))
rna_mat  <- as.matrix(read.delim(rna_file,  row.names = 1, check.names = FALSE))
coldata  <- read.delim(coldata_file, row.names = 1, check.names = FALSE,
                        stringsAsFactors = TRUE)

message("Ribo-seq genes: ", nrow(ribo_mat), "  samples: ", ncol(ribo_mat))
message("RNA-seq genes:  ", nrow(rna_mat),  "  samples: ", ncol(rna_mat))
message("colData samples:", nrow(coldata))

# ── 取交集基因 ───────────────────────────────────────────────────────────────
common_genes <- intersect(rownames(ribo_mat), rownames(rna_mat))
message("Common genes: ", length(common_genes))
ribo_mat <- ribo_mat[common_genes, , drop = FALSE]
rna_mat  <- rna_mat[common_genes, , drop = FALSE]

# ── 拆分 colData ─────────────────────────────────────────────────────────────
ribo_samples <- rownames(coldata)[coldata$data_type == "ribo"]
rna_samples  <- rownames(coldata)[coldata$data_type == "rna"]

# 确保 count matrix 列顺序与 colData 一致
ribo_mat <- ribo_mat[, ribo_samples, drop = FALSE]
rna_mat  <- rna_mat[, rna_samples,  drop = FALSE]

ribo_coldata <- coldata[ribo_samples, , drop = FALSE]
rna_coldata  <- coldata[rna_samples,  , drop = FALSE]

# condition 设 ref (第一个 levels 的 control 组)
# 让 R 自动用字母顺序, 但如果有 ctrl/control 就优先设 ref
conditions <- levels(coldata$condition)
ref_guess <- grep("ctrl|control", conditions, ignore.case = TRUE, value = TRUE)
if (length(ref_guess) > 0) {
  ribo_coldata$condition <- relevel(ribo_coldata$condition, ref = ref_guess[1])
  rna_coldata$condition  <- relevel(rna_coldata$condition,  ref = ref_guess[1])
  message("Reference condition: ", ref_guess[1])
}

# ── DESeq2 on Ribo-seq ───────────────────────────────────────────────────────
message("\n── DESeq2 on Ribo-seq ──")
dds_ribo <- DESeqDataSetFromMatrix(
  countData = round(ribo_mat),
  colData   = ribo_coldata,
  design    = ~ condition
)
keep_ribo <- rowSums(counts(dds_ribo) >= 10) >= ncol(dds_ribo) / 2
dds_ribo  <- dds_ribo[keep_ribo, ]
message("  Filtered genes (ribo): ", nrow(dds_ribo))
dds_ribo <- DESeq(dds_ribo, quiet = TRUE)
res_ribo <- results(dds_ribo)
message("  Ribo DE genes (padj<0.05): ",
        sum(!is.na(res_ribo$padj) & res_ribo$padj < 0.05))

# ── DESeq2 on RNA-seq ────────────────────────────────────────────────────────
message("\n── DESeq2 on RNA-seq ──")
dds_rna <- DESeqDataSetFromMatrix(
  countData = round(rna_mat),
  colData   = rna_coldata,
  design    = ~ condition
)
keep_rna <- rowSums(counts(dds_rna) >= 10) >= ncol(dds_rna) / 2
dds_rna  <- dds_rna[keep_rna, ]
message("  Filtered genes (rna): ", nrow(dds_rna))
dds_rna <- DESeq(dds_rna, quiet = TRUE)
res_rna <- results(dds_rna)
message("  RNA DE genes (padj<0.05): ",
        sum(!is.na(res_rna$padj) & res_rna$padj < 0.05))

# ── 合并结果 ─────────────────────────────────────────────────────────────────
message("\n── Merging results ──")
df_ribo <- as.data.frame(res_ribo)
df_rna  <- as.data.frame(res_rna)
colnames(df_ribo) <- paste0("ribo_", colnames(df_ribo))
colnames(df_rna)  <- paste0("rna_",  colnames(df_rna))

merged <- merge(df_ribo, df_rna, by = "row.names", all = TRUE)
colnames(merged)[1] <- "gene_id"

# TE 变化 = ribo lfc - rna lfc
merged$delta_logFC <- merged$ribo_log2FoldChange - merged$rna_log2FoldChange

# ── 分类 ─────────────────────────────────────────────────────────────────────
merged$category <- "NONE"
# 分类规则:
#   TRANSLATION_UP/DOWN:  ribo DE 显著 + rna 不显著 → TE 变化
#   TRANSCRIPTION_ONLY:    rna DE 显著 + ribo 不显著 → 只是转录变化
#   RIBO_ONLY:             ribo DE 显著 + rna DE 显著但方向相反 → 相互抵消/补偿
#   COORDINATED:           ribo DE 显著 + rna DE 显著 + 方向一致 → 协调变化
sig_ribo <- !is.na(merged$ribo_padj) & merged$ribo_padj < alpha
sig_rna  <- !is.na(merged$rna_padj)  & merged$rna_padj  < alpha

merged$category[sig_ribo & !sig_rna & merged$delta_logFC >  lfc_thresh] <- "TRANSLATION_UP"
merged$category[sig_ribo & !sig_rna & merged$delta_logFC < -lfc_thresh] <- "TRANSLATION_DOWN"
merged$category[!sig_ribo &  sig_rna]                                   <- "TRANSCRIPTION_ONLY"
merged$category[ sig_ribo &  sig_rna &
                 sign(merged$ribo_log2FoldChange) == sign(merged$rna_log2FoldChange)] <- "COORDINATED"
merged$category[ sig_ribo &  sig_rna &
                 sign(merged$ribo_log2FoldChange) != sign(merged$rna_log2FoldChange)] <- "RIBO_ONLY"

# ── 输出 ─────────────────────────────────────────────────────────────────────
merged <- merged[order(merged$delta_logFC), ]

all_out <- file.path(output_dir, "TE_separate_all.tsv")
write.table(merged, all_out, sep = "\t", row.names = FALSE, quote = FALSE)
message("All results written: ", all_out, " (", nrow(merged), " genes)")

sig_out <- file.path(output_dir, "TE_separate_sig.tsv")
sig_df  <- merged[merged$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN"), ]
write.table(sig_df, sig_out, sep = "\t", row.names = FALSE, quote = FALSE)
message("Significant TE changes: ", nrow(sig_df),
        " (", sum(sig_df$category == "TRANSLATION_UP"), " up, ",
        sum(sig_df$category == "TRANSLATION_DOWN"), " down)")

summary_df <- as.data.frame(table(merged$category), stringsAsFactors = FALSE)
colnames(summary_df) <- c("category", "n_genes")
summary_df$category <- factor(summary_df$category,
  levels = c("TRANSLATION_UP", "TRANSLATION_DOWN", "TRANSCRIPTION_ONLY",
             "COORDINATED", "RIBO_ONLY", "NONE"))
summary_df <- summary_df[order(summary_df$category), ]
write.table(summary_df,
            file.path(output_dir, "TE_separate_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Summary:")
print(summary_df)

# ── 火山图 (delta_logFC vs -log10 combined padj) ────────────────────────────
# combined_padj = Fisher's method (ribo + rna)
tryCatch({
  combined_p <- -log10(pmax(merged$ribo_padj, merged$rna_padj, na.rm = TRUE))
  pdf(file.path(output_dir, "TE_separate_volcano.pdf"), 7, 6)
  plot(merged$delta_logFC, -log10(combined_p),
       pch = 16, cex = 0.7, col = "grey70",
       xlab = expression(Delta*log[2]("FC"[ribo] - "FC"[rna])),
       ylab = expression(-log[10]("padj")),
       main = "TE change (Separate DESeq2)")
  # 标记显著 TE 变化
  sig_idx <- merged$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN")
  points(merged$delta_logFC[sig_idx], -log10(combined_p[sig_idx]),
         pch = 16, cex = 1,
         col = ifelse(merged$delta_logFC[sig_idx] > 0, "red", "blue"))
  abline(v = 0, lty = 2, col = "grey30")
  abline(h = -log10(alpha), lty = 3, col = "grey50")
  legend("topright",
         legend = c(paste0("TE up (", sum(merged$category == "TRANSLATION_UP"), ")"),
                    paste0("TE down (", sum(merged$category == "TRANSLATION_DOWN"), ")")),
         col = c("red", "blue"), pch = 16)
  dev.off()
  message("Volcano plot: TE_separate_volcano.pdf")
}, error = function(e) {
  message("Warning: volcano plot failed - ", e$message)
})

message("\n═══ TE_DA complete ═══")
