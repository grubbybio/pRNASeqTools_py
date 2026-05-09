#!/usr/bin/env Rscript
#
# ribotaper_filter_rsem.R — Filter GTF to expressed isoforms by TPM threshold
# based on RSEM quantification results.
#
# Usage:
#   Rscript ribotaper_filter_rsem.R <input_gtf> <output_gtf> <tpm_threshold> <isoform1.results> [...]
#
# Arguments:
#   input_gtf     — updated GTF with gene_biotype (from ribotaper_filter_gtf.R)
#   output_gtf    — output GTF with only expressed isoforms + gene rows
#   tpm_threshold — minimum mean TPM to keep an isoform (e.g. 0 or 0.5)
#   isoform*.results — RSEM isoform-level quantification files
#
# Filtering:
#   - Calculate mean TPM across all samples
#   - Keep isoforms with mean TPM > threshold (user-specified)
#   - Add gene-level rows for each gene
#   - Annotate biotype from input GTF
#
# Output GTF format compatible with RIBO Taper.

suppressMessages(library(dplyr))
options(warn = -1)

# ── Parse arguments ──────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript ribotaper_filter_rsem.R <input_gtf> <output_gtf> <tpm_threshold> <isoform.results> ...")
}

input_gtf     <- args[1]
output_gtf    <- args[2]
tpm_threshold <- as.numeric(args[3])
rsem_files    <- args[-(1:3)]

cat("Input GTF:", input_gtf, "\n")
cat("TPM threshold:", tpm_threshold, "\n")
cat("RSEM files:", length(rsem_files), "\n")
for (f in rsem_files) cat("  ", f, "\n")

# ── Load RSEM isoform results ────────────────────────────────────────────
tpm_list <- list()
for (i in seq_along(rsem_files)) {
  f <- rsem_files[i]
  if (!file.exists(f)) {
    cat("  WARNING: file not found, skipping:", f, "\n")
    next
  }
  d <- read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_name <- paste0("S", i)
  tpm_list[[sample_name]] <- data.frame(
    transcript_id = d$transcript_id,
    TPM = as.numeric(d$TPM),
    stringsAsFactors = FALSE
  )
}

if (length(tpm_list) == 0) {
  stop("No valid RSEM files found")
}

# ── Calculate mean TPM ───────────────────────────────────────────────────
tpm_df <- Reduce(function(x, y) merge(x, y, by = "transcript_id", all = TRUE), tpm_list)
tpm_cols <- grep("^TPM", names(tpm_df), value = TRUE)

for (col in tpm_cols) {
  tpm_df[[col]][is.na(tpm_df[[col]])] <- 0
}

tpm_df$TPM_mean <- rowMeans(tpm_df[, tpm_cols, drop = FALSE], na.rm = TRUE)

cat(sprintf("  Total isoforms: %d\n", nrow(tpm_df)))
cat(sprintf("  Expressed (mean TPM > %s): %d\n", tpm_threshold, sum(tpm_df$TPM_mean > tpm_threshold)))

# ── Filter expressed isoforms ────────────────────────────────────────────
expressed_ids <- tpm_df$transcript_id[tpm_df$TPM_mean > tpm_threshold]
cat(sprintf("  Unique expressed transcript IDs: %d\n", length(expressed_ids)))

# ── Read input GTF ───────────────────────────────────────────────────────
gtf <- read.delim(input_gtf,
                  header = FALSE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")

# Filter to features we care about: mRNA, exon, CDS
gtf <- gtf %>% filter(V3 %in% c("mRNA", "exon", "CDS"))

# Extract transcript_id
gtf$tx_id <- sapply(gtf$V9, function(x) {
  tx <- grep("transcript_id", unlist(strsplit(x, ";")), value = TRUE)
  if (length(tx) == 0) return(NA_character_)
  tx <- gsub("\"", "", gsub("transcript_id ", "", tx[1]))
  gsub(" ", "", tx, fixed = TRUE)
})

# Add gene_biotype to CDS rows (needed by RSEM/STAR)
cds_idx <- which(gtf$V3 == "CDS")
if (length(cds_idx) > 0) {
  gtf$V9[cds_idx] <- paste0(gtf$V9[cds_idx], ' gene_biotype "protein_coding";')
}

# Extract biotype from attribute column
gtf$biotype <- sapply(gtf$V9, function(x) {
  bt <- grep("gene_biotype", unlist(strsplit(x, ";")), value = TRUE)
  if (length(bt) == 0) return("protein_coding")
  gsub("\"", "", gsub(" gene_biotype ", "", bt[1]))
})

cat(sprintf("  GTF entries (mRNA/exon/CDS): %d\n", nrow(gtf)))

# ── Filter to expressed isoforms ─────────────────────────────────────────
gtf_expressed <- gtf %>% filter(tx_id %in% expressed_ids)
cat(sprintf("  Expressed GTF entries: %d\n", nrow(gtf_expressed)))

# ── Add gene-level rows ──────────────────────────────────────────────────
# Extract gene_id (everything before the last dot in transcript_id)
gtf_expressed$gene_id <- sub("^(.*)[.].*", "\\1", gtf_expressed$tx_id)

# Split by gene and add gene row
gtf_by_gene <- split(gtf_expressed, gtf_expressed$gene_id)

add_gene_row <- function(x) {
  mrna_rows <- x[x$V3 == "mRNA", , drop = FALSE]
  if (nrow(mrna_rows) == 0) return(x)

  gene_row <- mrna_rows[1, ]
  gene_row$V3 <- "gene"
  gene_row$V4 <- min(mrna_rows$V4)
  gene_row$V5 <- max(mrna_rows$V5)
  gene_row$V9 <- paste0('gene_id "', x$gene_id[1],
                        '"; gene_biotype "', x$biotype[1], '";')

  rbind(gene_row, x)
}

gtf_with_genes_list <- lapply(gtf_by_gene, add_gene_row)
gtf_final <- bind_rows(gtf_with_genes_list)

# Fix spacing before transcript_biotype
gtf_final$V9 <- gsub("transcript_biotype", " transcript_biotype", gtf_final$V9)

cat(sprintf("  Final GTF entries (with gene rows): %d\n", nrow(gtf_final)))

# ── Write output ─────────────────────────────────────────────────────────
write.table(gtf_final[, 1:9], output_gtf,
            quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")
cat("  Wrote:", output_gtf, "\n")
