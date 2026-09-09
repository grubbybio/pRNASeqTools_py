#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# TE_DS.R — Translational Efficiency changes via joint dual-factor DESeq2 (方法 C)
#
# 原理:
#   把 Ribo-seq 和 RNA-seq 的 count matrix 合并成一个大 matrix, 每个样本加两个因子:
#     condition (ctrl vs treatment) × data_type (ribo vs rna)
#
#   model:  log(count) ~ condition + data_type + condition:data_type
#
#   交互项 condition:data_type 的系数就是 TE 变化的直接估计:
#     beta_interaction = (logFC_ribo - logFC_rna)  =  delta_logFC
#
#   优势:
#     ✅ 一个模型同时检验转录变化、翻译变化和交互项
#     ✅ 共享离散度估计 (ribo/rna 基因间离散度相似)
#     ✅ 对每个基因一个 delta_logFC + 一个 padj (统计效力更高)
#     ✅ 不需要先分别做两次 DE 再手工合并
#
#   要求:
#     ⚠️ 每个实验条件都必须同时有 ribo-seq 和 rna-seq (paired design)
#     ⚠️ 如果某个条件只有 ribo 或只有 rna → 方法 C 不适用, 应降级到方法 B
#
# 输入:
#   args[1]: ribo_count_matrix.tsv
#   args[2]: rna_count_matrix.tsv
#   args[3]: colData.tsv   (sample  condition  data_type)
#   args[4]: output_dir
#   args[5]: alpha (default 0.05)
#   args[6]: lfc_threshold (default 0)
#
# 输出:
#   TE_joint_all.tsv        — 所有基因
#   TE_joint_sig.tsv        — 交互项显著的 TE 变化基因
#   TE_joint_summary.tsv    — 分类汇总
#   TE_joint_volcano.pdf    — 交互项火山图
# ─────────────────────────────────────────────────────────────────────────────

options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript TE_DS.R <ribo_count.tsv> <rna_count.tsv> ",
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
message("  TE_DS: Joint dual-factor DESeq2 (方法 C)")
message("═══════════════════════════════════════════")

# ── 读取数据 ─────────────────────────────────────────────────────────────────
ribo_mat <- as.matrix(read.delim(ribo_file, row.names = 1, check.names = FALSE))
rna_mat  <- as.matrix(read.delim(rna_file,  row.names = 1, check.names = FALSE))
coldata  <- read.delim(coldata_file, row.names = 1, check.names = FALSE,
                        stringsAsFactors = TRUE)

message("Ribo-seq genes: ", nrow(ribo_mat), "  samples: ", ncol(ribo_mat))
message("RNA-seq genes:  ", nrow(rna_mat),  "  samples: ", ncol(rna_mat))
message("colData samples:", nrow(coldata))

# ── 配对检查 ─────────────────────────────────────────────────────────────────
# 方法 C 要求每个 condition 同时有 ribo 和 rna
message("\n── Paired design check ──")
conditions     <- levels(coldata$condition)
data_types     <- levels(coldata$data_type)
ribo_per_cond  <- table(coldata$condition[coldata$data_type == "ribo"])
rna_per_cond   <- table(coldata$condition[coldata$data_type == "rna"])
paired_ok <- TRUE
for (c in conditions) {
  nr <- as.integer(ribo_per_cond[c])
  nn <- as.integer(rna_per_cond[c])
  message(sprintf("  %-15s: ribo=%d  rna=%d %s",
                  c, nr, nn, ifelse(nr == nn, "✓", "⚠ unpaired")))
  if (nr != nn) paired_ok <- FALSE
}

if (!paired_ok) {
  message("\n⚠ 警告: 方法 C 要求每个 condition 有相同数目的 ribo/rna 重复.")
  message("  可以继续运行, 但交互项统计解释需要小心.")
}

# ── 合并 count matrix ────────────────────────────────────────────────────────
common_genes <- intersect(rownames(ribo_mat), rownames(rna_mat))
message("\nCommon genes: ", length(common_genes))

ribo_mat <- ribo_mat[common_genes, , drop = FALSE]
rna_mat  <- rna_mat[common_genes, , drop = FALSE]

# 合并成一个大 matrix, 列顺序 = colData row order
all_samples <- rownames(coldata)
counts_combined <- cbind(round(ribo_mat), round(rna_mat))
# 按 colData 顺序重新排列列
counts_combined <- counts_combined[, all_samples, drop = FALSE]
message("Combined count matrix: ", nrow(counts_combined), " genes × ",
        ncol(counts_combined), " samples")

# ── 构造 DESeq2 对象 ─────────────────────────────────────────────────────────
# condition 设 ref
ref_guess <- grep("ctrl|control", conditions, ignore.case = TRUE, value = TRUE)
if (length(ref_guess) > 0) {
  coldata$condition   <- relevel(coldata$condition, ref = ref_guess[1])
  coldata$data_type   <- relevel(coldata$data_type, ref = "rna")
  message("Reference: condition=", ref_guess[1], ", data_type=rna")
}

message("\n── DESeq2: ~ condition + data_type + condition:data_type ──")
dds <- DESeqDataSetFromMatrix(
  countData = counts_combined,
  colData   = coldata,
  design    = ~ condition + data_type + condition:data_type
)

# 过滤: 至少一半样本 (ribo+rna 合并) 有 >=10 count
keep <- rowSums(counts(dds) >= 10) >= ncol(dds) / 2
dds  <- dds[keep, ]
message("  Filtered genes: ", nrow(dds), "/", nrow(counts_combined))

dds <- DESeq(dds, quiet = TRUE)

# ── 提取结果 ─────────────────────────────────────────────────────────────────
# resultsNames: "Intercept", "condition_treatment_vs_ctrl",
#               "data_type_ribo_vs_rna",
#               "conditiontreatment.data_typeribo" (交互项!)
rn <- resultsNames(dds)
message("\nResults names: ", paste(rn, collapse = ", "))

# 找交互项
interaction_name <- rn[grepl(":", rn) | grepl("condition.*data_type|data_type.*condition", rn)]
if (length(interaction_name) == 0) {
  message("⚠ 交互项名不明确, 用最后一个 results()")
  interaction_name <- rn[length(rn)]
}
message("Interaction term: ", interaction_name)

# 1. 交互项 = delta_logFC (TE 变化)
res_interaction <- results(dds, name = interaction_name)

# 2. condition 主效应 = RNA 转录变化
cond_name <- rn[grepl("^condition", rn) & !grepl(":", rn)][1]
res_condition <- results(dds, name = cond_name)

# 3. data_type 主效应 = ribo vs rna 总体差异
dt_name <- rn[grepl("^data_type", rn) & !grepl(":", rn)][1]
res_datatype <- results(dds, name = dt_name)

# ── 合并 ─────────────────────────────────────────────────────────────────────
df_int <- as.data.frame(res_interaction)
df_cond <- as.data.frame(res_condition)
df_dt <- as.data.frame(res_datatype)

colnames(df_int)  <- paste0("TE_int_", colnames(df_int))
colnames(df_cond) <- paste0("rna_",  colnames(df_cond))
colnames(df_dt)   <- paste0("ribo_", colnames(df_dt))

merged <- merge(df_int, df_cond, by = "row.names", all = TRUE)
merged <- merge(merged, df_dt,    by.x = "Row.names", by.y = "row.names", all = TRUE)
colnames(merged)[1] <- "gene_id"

# delta_logFC = 交互项系数 (已经是 TE 变化的估计!)
merged$delta_logFC <- merged$TE_int_log2FoldChange

# ── 分类 ─────────────────────────────────────────────────────────────────────
merged$category <- "NONE"
sig_int  <- !is.na(merged$TE_int_padj)  & merged$TE_int_padj  < alpha
sig_rna  <- !is.na(merged$rna_padj)     & merged$rna_padj     < alpha
sig_ribo <- !is.na(merged$ribo_padj)    & merged$ribo_padj    < alpha

merged$category[ sig_int &  merged$delta_logFC >  lfc_thresh] <- "TRANSLATION_UP"
merged$category[ sig_int &  merged$delta_logFC < -lfc_thresh] <- "TRANSLATION_DOWN"
merged$category[!sig_int &  sig_rna]                                   <- "TRANSCRIPTION_ONLY"
merged$category[ sig_int &  sig_rna & sign(merged$rna_log2FoldChange) == sign(merged$delta_logFC)] <- "COORDINATED"

# ── 输出 ─────────────────────────────────────────────────────────────────────
merged <- merged[order(merged$delta_logFC), ]

all_out <- file.path(output_dir, "TE_joint_all.tsv")
write.table(merged, all_out, sep = "\t", row.names = FALSE, quote = FALSE)
message("\nAll results written: ", all_out, " (", nrow(merged), " genes)")

sig_out <- file.path(output_dir, "TE_joint_sig.tsv")
sig_df  <- merged[merged$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN"), ]
write.table(sig_df, sig_out, sep = "\t", row.names = FALSE, quote = FALSE)
message("Significant TE changes: ", nrow(sig_df),
        " (", sum(sig_df$category == "TRANSLATION_UP"), " up, ",
        sum(sig_df$category == "TRANSLATION_DOWN"), " down)")

summary_df <- as.data.frame(table(merged$category), stringsAsFactors = FALSE)
colnames(summary_df) <- c("category", "n_genes")
summary_df$category <- factor(summary_df$category,
  levels = c("TRANSLATION_UP", "TRANSLATION_DOWN", "TRANSCRIPTION_ONLY",
             "COORDINATED", "NONE"))
summary_df <- summary_df[order(summary_df$category), ]
write.table(summary_df,
            file.path(output_dir, "TE_joint_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Summary:")
print(summary_df)

# ── 火山图 ──────────────────────────────────────────────────────────────────
tryCatch({
  pdf(file.path(output_dir, "TE_joint_volcano.pdf"), 7, 6)
  plot(merged$delta_logFC, -log10(merged$TE_int_padj),
       pch = 16, cex = 0.7, col = "grey70",
       xlab = expression(Delta*log[2]("FC"[TE])),
       ylab = expression(-log[10]("padj"[interaction])),
       main = "TE change (Joint DESeq2, interaction term)")
  sig_idx <- merged$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN")
  points(merged$delta_logFC[sig_idx], -log10(merged$TE_int_padj[sig_idx]),
         pch = 16, cex = 1,
         col = ifelse(merged$delta_logFC[sig_idx] > 0, "red", "blue"))
  abline(v = 0, lty = 2, col = "grey30")
  abline(h = -log10(alpha), lty = 3, col = "grey50")
  legend("topright",
         legend = c(paste0("TE up (", sum(merged$category == "TRANSLATION_UP"), ")"),
                    paste0("TE down (", sum(merged$category == "TRANSLATION_DOWN"), ")")),
         col = c("red", "blue"), pch = 16)
  dev.off()
  message("Volcano plot: TE_joint_volcano.pdf")
}, error = function(e) {
  message("Warning: volcano plot failed - ", e$message)
})

# ── MA plot (mean expr vs delta_logFC) ────────────────────────────────────────
tryCatch({
  mean_expr <- rowMeans(counts(dds, normalized = TRUE))
  pdf(file.path(output_dir, "TE_joint_MAPlot.pdf"), 7, 6)
  plot(log2(mean_expr + 1), merged$delta_logFC,
       pch = 16, cex = 0.7, col = "grey70",
       xlab = expression(log[2]("mean normalized count + 1")),
       ylab = expression(Delta*log[2]("FC"[TE])),
       main = "MA plot (Joint DESeq2)")
  sig_idx <- merged$category %in% c("TRANSLATION_UP", "TRANSLATION_DOWN")
  points(log2(mean_expr + 1)[sig_idx], merged$delta_logFC[sig_idx],
         pch = 16, cex = 1,
         col = ifelse(merged$delta_logFC[sig_idx] > 0, "red", "blue"))
  abline(h = 0, lty = 2, col = "grey30")
  dev.off()
  message("MA plot: TE_joint_MAPlot.pdf")
}, error = function(e) {
  message("Warning: MA plot failed - ", e$message)
})

message("\n═══ TE_DS complete ═══")
