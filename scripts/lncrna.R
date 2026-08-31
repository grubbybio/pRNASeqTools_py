# ── lncRNA-seq DE analysis (DESeq2, split by transcript class) ───
# Usage:
#   Rscript --vanilla lncrna.R <foldchange> <pvalue> <fdr>
#     <full_merged_gtf> <tag1>,<tag2>,...
#
# Expected per-sample files:  {tag}.txt  (Transcript\tCount)
# GTF must have 'class' attribute: protein_coding | known_ncRNA |
#                                  novel_coding | novel_lncRNA
#
# Outputs:
#   total.xxx.csv / .pdf        — combined (all classes)
#   protein_coding.xxx.csv      — coding transcripts only
#   novel_lncRNA.xxx.csv        — novel lncRNA only
#   known_ncRNA.xxx.csv         — known ncRNA only
#   summary_table.csv           — class-level counts per comparison
#   class_barplot.pdf           — DE transcript counts by class

options(warn = -1)
suppressMessages(library(DESeq2))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

args = commandArgs(trailingOnly = TRUE)
fc.thresh   <- as.numeric(args[1])
p.thresh    <- as.numeric(args[2])
fdr.thresh  <- as.numeric(args[3])
merged_gtf  <- args[4]
par_str     <- args[-c(1:4)]

# ── Parse {tag},{replicates},{tag},{replicates},... ─────────────────────
# Mirror mrna/DEG.R convention:  args[odd]=genotype, args[even]=replicates
genotype   <- par_str[c(TRUE, FALSE)]
replicates <- as.integer(par_str[c(FALSE, TRUE)])
genotype   <- genotype[!is.na(genotype) & genotype != ""]
replicates <- replicates[!is.na(replicates)]

if (length(genotype) < 2) {
  stop("Need at least 2 groups (control + treatment)")
}

message("╔═══════════════════════════════════════════════════════╗")
message("║  lncRNA-seq Differential Expression Analysis             ║")
message("╚═══════════════════════════════════════════════════════╝")
message(paste("  Groups:   ", paste(genotype, replicates, collapse = ", ")))
message(paste("  FC ≥", fc.thresh, "P <", p.thresh, "FDR <", fdr.thresh))

# ── 1. Load per-sample count tables ──────────────────────────────────────
all_tags <- character()
for (i in seq_along(genotype)) {
  for (k in seq_len(replicates[i])) {
    all_tags <- c(all_tags, paste0(genotype[i], "_", k))
  }
}

count_list <- list()
for (tag in all_tags) {
  fn <- paste0(tag, ".txt")
  if (!file.exists(fn)) {
    message(paste("  ✗ Missing count file:", fn))
    next
  }
  df <- read.table(fn, header = TRUE, row.names = 1, as.is = TRUE,
                   check.names = FALSE, comment.char = "")
  colnames(df) <- tag
  count_list[[tag]] <- df
}

if (length(count_list) == 0) stop("✗ No count files found!")

# Merge into matrix
mat <- do.call(cbind, lapply(count_list, function(x) x[, 1, drop = FALSE]))
mat <- mat[rowSums(mat) >= 1, , drop = FALSE]
message(paste("  Transcripts (raw):", nrow(mat)))

# ── 2. Parse merged GTF → transcript class lookup ───────────────────────
class_map <- character()     # names = transcript_id, values = class
if (file.exists(merged_gtf)) {
  con <- file(merged_gtf, open = "r")
  while ((line <- readLines(con, n = 1, warn = FALSE)) != "") {
    if (grepl("^#", line)) next
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 9 || parts[3] != "transcript") next

    attrs <- parts[9]
    tid   <- sub('.*transcript_id "([^"]+)".*', "\\1", attrs)
    cls   <- sub('.*class "([^"]+)".*', "\\1", attrs)

    if (nchar(cls) == 0) cls <- "unclassified"
    class_map[tid] <- cls
  }
  close(con)
} else {
  message(paste("  ⚠ GTF not found:", merged_gtf,
                "→ all transcripts marked 'unclassified'"))
}

# Fallback for missing tid
class_map <- class_map[rownames(mat)]
class_map[is.na(class_map)] <- "unclassified"

# Group into 4 major classes
class_bucket <- dplyr::recode(
  class_map,
  protein_coding = "protein_coding",
  novel_coding   = "protein_coding",
  known_ncRNA    = "known_ncRNA",
  novel_lncRNA   = "novel_lncRNA",
  .default       = "other"
)
names(class_bucket) <- names(class_map)

message("  Transcript class distribution (merged GTF):")
print(table(class_bucket))

# ── 3. DESeq2 on full matrix ───────────────────────────────────────────
grp_vec <- rep(genotype, replicates)
coldata <- data.frame(row.names = all_tags,
                      genotype = factor(grp_vec, levels = genotype))

dds <- DESeqDataSetFromMatrix(countData = round(mat),
                               colData   = coldata,
                               design    = ~ genotype)
keep <- rowMeans(counts(dds, normalized = TRUE)) >= 1
dds <- dds[keep, ]
dds <- DESeq(dds)

# ── 4. Per-class subset + DE + visualization ────────────────────────────
summary_list <- list()

DEG <- function(dds_obj, from_class, class_label, summary_list) {
  tids <- rownames(dds_obj)[class_bucket[rownames(dds_obj)] == class_label]

  if (length(tids) == 0) {
    message(paste0("  [", class_label, "] 0 transcripts — skip"))
    summary_list[[class_label]] <- data.frame(
      comparison = character(),
      total = integer(), up = integer(), down = integer())
    return(summary_list)
  }

  dds_sub <- dds_obj[tids, , drop = FALSE]
  rld_sub <- rlog(dds_sub, blind = TRUE)

  for (j in 2:length(genotype)) {
    g1 <- genotype[1]; g2 <- genotype[j]
    contrast <- c("genotype", g2, g1)
    res <- results(dds_sub, contrast = contrast,
                   pAdjustMethod = "BH")

    fc_pass   <- abs(log2(res$padj))  # placeholder for log2FC col access
    up_mask   <- res$padj < fdr.thresh & !is.na(res$padj) &
                 res$log2FoldChange >  log(fc.thresh, 2)
    down_mask <- res$padj < fdr.thresh & !is.na(res$padj) &
                 res$log2FoldChange < -log(fc.thresh, 2)
    n_up   <- sum(up_mask,   na.rm = TRUE)
    n_down <- sum(down_mask, na.rm = TRUE)
    n_tot  <- nrow(res)

    out_fname <- paste0(class_label, "_", g2, "vs", g1, ".csv")
    out <- as.data.frame(res)
    out$transcript_id <- rownames(out)
    out$class <- class_bucket[rownames(res)]
    write.csv(out, out_fname, row.names = FALSE, quote = FALSE)

    summary_list[[class_label]] <- rbind(
      summary_list[[class_label]],
      data.frame(comparison = paste0(g2, "_vs_", g1),
                 total = n_tot, up = n_up, down = n_down,
                 row.names = NULL)
    )

    # Volcano plot
    vol_df <- data.frame(
      tid = rownames(res),
      log2FC = res$log2FoldChange,
      padj  = -log10(res$padj)
    )
    pdf(paste0(class_label, "_", g2, "vs", g1, "_volcano.pdf"),
        width = 6, height = 5)
    plot(vol_df$log2FC, vol_df$padj, pch = 16, cex = 0.4,
         col = "grey50",
         xlab = "log2 FC", ylab = "-log10 FDR",
         main = paste0(class_label, " — ", g2, " vs ", g1))
    points(vol_df$log2FC[up_mask],   vol_df$padj[up_mask],
           pch = 16, cex = 0.5, col = "red")
    points(vol_df$log2FC[down_mask], vol_df$padj[down_mask],
           pch = 16, cex = 0.5, col = "blue")
    abline(h = -log10(fdr.thresh), v = c(-log(fc.thresh, 2),
                                          log(fc.thresh, 2)),
           col = "grey30", lty = 2)
    dev.off()

    message(paste0("  [", class_label, "] ", g2, " vs ", g1,
                   ": total=", n_tot, "  up=", n_up, "  down=", n_down))
  }

  return(summary_list)
}

# Run DE on 4 transcript classes
summary_list <- DEG(dds, "protein_coding", "protein_coding", summary_list)
summary_list <- DEG(dds, "known_ncRNA",    "known_ncRNA",    summary_list)
summary_list <- DEG(dds, "novel_lncRNA",   "novel_lncRNA",   summary_list)
summary_list <- DEG(dds, "other",          "other",          summary_list)

# ── 5. Combined "all transcripts" comparison (for legacy) ──────────────
for (j in 2:length(genotype)) {
  g1 <- genotype[1]; g2 <- genotype[j]
  contrast <- c("genotype", g2, g1)
  res <- results(dds, contrast = contrast, pAdjustMethod = "BH")
  up_mask   <- res$padj < fdr.thresh & !is.na(res$padj) &
               res$log2FoldChange >  log(fc.thresh, 2)
  down_mask <- res$padj < fdr.thresh & !is.na(res$padj) &
               res$log2FoldChange < -log(fc.thresh, 2)
  out <- as.data.frame(res)
  out$transcript_id <- rownames(out)
  out$class <- class_bucket[rownames(res)]
  write.csv(out, paste0("total_", g2, "vs", g1, ".csv"),
            row.names = FALSE, quote = FALSE)
}

# ── 6. Summary table ────────────────────────────────────────────────────
summary_table <- do.call(rbind, lapply(names(summary_list), function(cls) {
  df <- summary_list[[cls]]
  df$class <- cls
  df[, c("class", "comparison", "total", "up", "down")]
}))
write.csv(summary_table, "summary_table.csv", row.names = FALSE, quote = FALSE)
message("\n  Summary table written → summary_table.csv")

# ── 7. Class barplot ───────────────────────────────────────────────────
if (nrow(summary_table) > 0) {
  pdf("class_barplot.pdf", width = 8, height = 5)
  dat <- summary_table %>%
    filter(!is.na(up)) %>%
    mutate(up_down = up + down)

  # Reshape for grouped barplot
  ptab <- xtabs(cbind(up, down) ~ class + comparison, data = summary_table)
  barplot(t(ptab), beside = TRUE,
          col = c("red3", "blue3"),
          legend.text = TRUE,
          xlab = "Transcript class",
          ylab = "Count",
          main = "DE transcripts by class")
  dev.off()
  message("  Barplot written → class_barplot.pdf")
}

# ── 8. PCA + heatmap (protein_coding subset) ───────────────────────────
rld <- rlog(dds, blind = TRUE)

pdf("total_PCA.pdf", width = 6, height = 5)
plotPCA(rld, intgroup = "genotype", ntop = 1000)
dev.off()

# Heatmap of top variable coding transcripts
coding_rows <- rownames(mat)[class_bucket[rownames(mat)] == "protein_coding"]
if (length(coding_rows) > 0) {
  sub_m <- mat[coding_rows, , drop = FALSE]
  rv <- rowVars <- apply(log2(sub_m + 1), 1, var)
  top_idx <- order(rv, decreasing = TRUE)[seq_len(min(1000, nrow(sub_m)))]
  sub_m <- sub_m[top_idx, , drop = FALSE]

  ann_col <- data.frame(genotype = coldata$genotype,
                        row.names = rownames(coldata))
  pdf("top1000_coding_heatmap.pdf", width = 8, height = 10)
  pheatmap(log2(sub_m + 1), annotation_col = ann_col,
           show_rownames = FALSE,
           color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
           main = "Top 1000 variable protein-coding transcripts")
  dev.off()
  message("  Heatmap written → top1000_coding_heatmap.pdf")
}

# ── Summary ─────────────────────────────────────────────────────────────
message("\n╔═══════════════════════════════════════════════════════╗")
message("║  lncRNA-seq DE completed                                 ║")
message("╚═══════════════════════════════════════════════════════╝")
message(paste0("  Input transcripts: ", nrow(mat)))
message(paste0("  Classes: ", paste(unique(class_bucket), collapse = ", ")))
message(paste0("  Groups:  ", paste(genotype, collapse = ", ")))
message("  Key outputs:")
message("    summary_table.csv        — class-level DE summary")
message("    total_*.csv              — all transcripts combined DE")
message("    protein_coding_*.csv     — protein-coding DE")
message("    novel_lncRNA_*.csv       — novel lncRNA DE")
message("    known_ncRNA_*.csv        — known ncRNA DE")
message("    class_barplot.pdf        — DE counts by class")
message("    *_volcano.pdf            — per-class volcano plots")
message("\n✓ Done.")
