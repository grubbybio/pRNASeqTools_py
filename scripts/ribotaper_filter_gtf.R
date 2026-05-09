#!/usr/bin/env Rscript
#
# ribotaper_filter_gtf.R — Filter novel transcripts from gffcompare output
# and add gene_biotype annotation for RIBO Taper compatibility.
#
# Usage:
#   Rscript ribotaper_filter_gtf.R <annotated_gtf> <original_gtf> <output_gtf>
#
# Arguments:
#   annotated_gtf  — gffcompare .annotated.gtf (merged + reference)
#   original_gtf   — original reference GTF/GFF
#   output_gtf     — output GTF with novel transcripts + gene_biotype
#
# Novel transcript class codes selected:
#   i — fully contained within reference intron
#   x — exonic overlap on opposite strand
#   y — contains a reference transcript
#   o — generic exonic overlap
#   u — intergenic
#   s — intronic overlap on opposite strand
#
# Output GTF format compatible with RIBO Taper, RSEM, and STAR.

suppressMessages(library(dplyr))
options(warn = -1)

# ── Parse arguments ──────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript ribotaper_filter_gtf.R <annotated_gtf> <original_gtf> <output_gtf>")
}

annotated_gtf <- args[1]
original_gtf  <- args[2]
output_gtf    <- args[3]

cat("Reading annotated GTF:", annotated_gtf, "\n")
cat("Reading original GTF:", original_gtf, "\n")

# ── Read annotated GTF from gffcompare ───────────────────────────────────
gtf <- read.delim(annotated_gtf,
                  header = FALSE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")

# Split by transcript
is_transcript <- gtf$V3 == "transcript"
gtf_list <- split(gtf, cumsum(is_transcript))

# Extract class_code for each transcript
gtf_class <- sapply(gtf_list, function(x) {
  cc <- grep("class_code", unlist(strsplit(x[1, 9], ";")), value = TRUE)
  if (length(cc) == 0) return(NA_character_)
  gsub("\"", "", gsub(" class_code ", "", cc))
})

# Select novel class codes: i, x, y, o, u, s
novel_codes <- c("i", "x", "y", "o", "u", "s")
novel_idx <- which(gtf_class %in% novel_codes)
cat(sprintf("  Found %d novel transcripts (class codes: %s)\n",
            length(novel_idx), paste(novel_codes, collapse = ",")))

if (length(novel_idx) == 0) {
  cat("  No novel transcripts found. Copying original GTF as output.\n")
  file.copy(original_gtf, output_gtf, overwrite = TRUE)
  quit(save = "no", status = 0)
}

gtf_novel_list <- gtf_list[novel_idx]

# Order by original index
novel_names <- as.numeric(names(gtf_novel_list))
gtf_novel_list <- gtf_novel_list[order(novel_names)]

# Convert list to data.frame
gtf_novel <- bind_rows(gtf_novel_list)

# ── Extract transcript_id, gene_id, class_code ──────────────────────────
gtf_novel$transcript_id <- sapply(seq_len(nrow(gtf_novel)), function(i) {
  gsub("\"", "", gsub("transcript_id ", "",
       unlist(strsplit(gtf_novel$V9[i], "; "))[1]))
})

gtf_novel$gene_id <- sapply(seq_len(nrow(gtf_novel)), function(i) {
  parts <- unlist(strsplit(gtf_novel$V9[i], "; "))
  gi <- grep("gene_id ", parts, value = TRUE)
  if (length(gi) > 0) gsub("\"", "", gsub("gene_id ", "", gi[1]))
  else NA_character_
})

gtf_novel$class_code <- sapply(seq_len(nrow(gtf_novel)), function(i) {
  cc <- grep("class_code", unlist(strsplit(gtf_novel$V9[i], ";")), value = TRUE)
  if (length(cc) == 0) return(NA_character_)
  gsub("\"", "", gsub(" class_code ", "", cc))
})

# Fill NA class_codes from previous row (propagate)
for (i in seq_len(nrow(gtf_novel))) {
  if (is.na(gtf_novel$class_code[i]) && i > 1) {
    gtf_novel$class_code[i] <- gtf_novel$class_code[i - 1]
  }
}

# ── Add gene rows and gene_biotype for novel genes ──────────────────────
novel_list <- split(gtf_novel, gtf_novel$gene_id)

gene_rows_fn <- function(x) {
  # Add transcript_biotype and gene_biotype to all rows
  x$V9 <- paste0(x$V9, ' transcript_biotype "ncRNA"; gene_biotype "ncRNA";')

  # Create gene row from first mRNA/transcript
  mrna_rows <- x[x$V3 %in% c("transcript", "mRNA"), , drop = FALSE]
  if (nrow(mrna_rows) == 0) return(x[, 1:9])

  gene_row <- mrna_rows[1, ]
  gene_row$V3 <- "gene"
  gene_row$V4 <- min(mrna_rows$V4)
  gene_row$V5 <- max(mrna_rows$V5)
  gene_row$V9 <- paste0('gene_id "', x$gene_id[1],
                        '"; gene_biotype "ncRNA"; class_code ',
                        x$class_code[1], ';')

  rbind(gene_row, x)[, 1:9]
}

gtf_novel_with_genes <- lapply(novel_list, gene_rows_fn)
gtf_novel_final <- do.call("rbind", gtf_novel_with_genes)

# Change "transcript" to "mRNA" for RIBO Taper compatibility
gtf_novel_final$V3[gtf_novel_final$V3 == "transcript"] <- "mRNA"

# ── Read original GTF ────────────────────────────────────────────────────
original <- read.delim(original_gtf,
                       header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")

# ── Combine original + novel ─────────────────────────────────────────────
combined <- rbind(original, gtf_novel_final)
cat(sprintf("  Combined GTF: %d lines (original: %d, novel: %d)\n",
            nrow(combined), nrow(original), nrow(gtf_novel_final)))

# ── Write output ─────────────────────────────────────────────────────────
write.table(combined[, 1:9], output_gtf,
            quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")
cat("  Wrote:", output_gtf, "\n")
