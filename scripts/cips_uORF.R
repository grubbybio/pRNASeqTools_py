#!/usr/bin/env Rscript
#
# cips_uORF.R — CiPS uORF analysis pipeline
# Identify translated upstream ORFs (uORFs) in 5' UTRs using Ribo-seq P-sites.
#
# Usage:
#   Rscript cips_uORF.R <gtf> <fasta> <psite_file> <output_prefix>
#       [min_inframe_counts] [min_inframe_perc] [min_psite_perc] [gene_desc_xlsx]
#
# Arguments:
#   gtf               — Expressed GTF (post-RiboTaper, with gene_biotype)
#   fasta             — Reference genome FASTA (indexed, with .fai)
#   psite_file        — P-site count file (columns: count, chr, start, strand)
#   output_prefix     — Output file prefix
#   min_inframe_counts — Minimum in-frame Ribo-seq counts (default: 10)
#   min_inframe_perc   — Minimum in-frame count percentage (default: 50)
#   min_psite_perc     — Minimum in-frame P-site percentage for longer uORFs (default: 30)
#   gene_desc_xlsx     — Optional gene description Excel for annotation
#
# Outputs:
#   <output_prefix>_minimum_uORF.xlsx  — Minimum uORFs (1 aa)
#   <output_prefix>_uORF.xlsx          — All translated uORFs (≥2 aa)
#   <output_prefix>_uORF_dedup.xlsx    — Deduplicated (overlapping removed)
#   <output_prefix>_uORF_overlaps.xlsx — Overlapping uORF pairs

suppressMessages({
  library(dplyr)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
  library(ORFik)
  library(future.apply)
  library(openxlsx)
})
options(warn = -1)

# ── Parse arguments ──────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript cips_uORF.R <gtf> <fasta> <psite_file> <output_prefix> \\
       [min_inframe_counts] [min_inframe_perc] [min_psite_perc] [gene_desc_xlsx]")
}

gtf_path           <- args[1]
fasta_path         <- args[2]
psite_file         <- args[3]
output_prefix      <- args[4]
min_inframe_counts <- ifelse(length(args) >= 5, as.numeric(args[5]), 10)
min_inframe_perc   <- ifelse(length(args) >= 6, as.numeric(args[6]), 50)
min_psite_perc     <- ifelse(length(args) >= 7, as.numeric(args[7]), 30)
gene_desc_xlsx     <- ifelse(length(args) >= 8, args[8], NA_character_)

cat("=== CiPS uORF Analysis ===\n")
cat("GTF:           ", gtf_path, "\n")
cat("FASTA:         ", fasta_path, "\n")
cat("P-site file:   ", psite_file, "\n")
cat("Output prefix: ", output_prefix, "\n")
cat("Min in-frame counts:", min_inframe_counts, "\n")
cat("Min in-frame %:     ", min_inframe_perc, "\n")
cat("Min P-site %:       ", min_psite_perc, "\n")
cat("Gene desc xlsx:", ifelse(is.na(gene_desc_xlsx), "(none)", gene_desc_xlsx), "\n")

# ── Load reference data ──────────────────────────────────────────────────
cat("\n--- Loading reference data ---\n")
FA <- FaFile(fasta_path)
txdb <- makeTxDbFromGFF(file = gtf_path, format = "gtf",
                         dataSource = "reference", organism = "unknown")
exonByTx <- exonsBy(txdb, by = "tx", use.names = TRUE)
cat(sprintf("  Transcripts: %d\n", length(exonByTx)))

# Extract 5' UTRs
fiveUTR <- fiveUTRsByTranscript(txdb, use.names = TRUE)
fiveUTR_seqs <- extractTranscriptSeqs(FA, fiveUTR)

# Find all ORFs in 5' UTRs
fiveUTR_ORFs <- findMapORFs(fiveUTR, fiveUTR_seqs,
                             startCodon = "ATG", longestORF = TRUE,
                             groupByTx = FALSE)
uORF_length <- sum(width(fiveUTR_ORFs))
fiveUTR_ORFs_seq <- extractTranscriptSeqs(FA, fiveUTR_ORFs)

cat(sprintf("  5' UTRs:  %d\n", length(fiveUTR)))
cat(sprintf("  uORFs found: %d\n", length(fiveUTR_ORFs)))

# ── Load P-site data ─────────────────────────────────────────────────────
cat("\n--- Loading P-site data ---\n")
Ribo1 <- read.delim(file = psite_file, header = FALSE,
                     stringsAsFactors = FALSE, sep = "\t")
colnames(Ribo1) <- c("count", "chr", "start", "strand")
Ribo1$end <- Ribo1$start
RiboR <- makeGRangesFromDataFrame(Ribo1, keep.extra.columns = TRUE)
cat(sprintf("  P-site positions: %d\n", nrow(Ribo1)))

# ── Core function: frame-specific P-site counts per uORF ─────────────────
XAA <- function(x, y) {
  # x: index into fiveUTR_ORFs_xaa
  # y: number of amino acids
  STRAND <- as.character(strand(fiveUTR_ORFs_xaa[[x]])[1])
  if (STRAND == "+") {
    z <- fiveUTR_ORFs_xaa[[x]]
  } else if (STRAND == "-") {
    z <- rev(fiveUTR_ORFs_xaa[[x]])
  }
  aa <- unlist(Map(`:`, start(z), end(z)))
  gr <- GRanges(seqnames = as.character(seqnames(z))[1],
                ranges = IRanges(start = aa, width = 1),
                strand = STRAND)

  if (STRAND == "+") {
    gr$frame <- rep(1:3, time = length(aa) / 3)
    gr1 <- gr[1:(3 * (y - 1))]
    gr2 <- gr[(3 * (y - 1) + 1):(3 * y)]
    gr3 <- gr[(3 * y + 1):(3 * (y + 1))]
  } else {
    gr$frame <- rep(3:1, time = length(aa) / 3)
    gr1 <- gr[(3 * (y + 1)):7]
    gr2 <- gr[6:4]
    gr3 <- gr[3:1]
  }

  df0 <- data.frame(count = c(0, 0, 0), frame = c(1, 2, 3))

  # Helper: overlap RiboR with grX and compute frame counts
  compute_region <- function(gr_target) {
    ranges <- subsetByOverlaps(RiboR, gr_target)
    hits <- findOverlaps(RiboR, gr_target)
    xframe <- unlist(IntegerList(split(gr_target$frame[subjectHits(hits)],
                                       queryHits(hits))))
    mcols(ranges) <- DataFrame(mcols(ranges), xframe)
    df <- data.frame(count = ranges$count, frame = ranges$xframe)
    df <- rbind(df, df0)
    df_t <- df %>%
      group_by(frame) %>%
      summarise(total_counts = sum(count), .groups = "drop") %>%
      as.data.frame()
    freq_t <- table(df$frame) %>% as.data.frame()
    df_t$frame_freq <- freq_t$Freq - 1
    df_t
  }

  df1t <- compute_region(gr1)
  df2t <- compute_region(gr2)
  df3t <- compute_region(gr3)

  DF <- cbind(df1t, df2t[, 2:3], df3t[, 2:3])
  colnames(DF) <- c("Frame", "Counts1", "FrameFreq1",
                    "CountsLast", "FrameFreqLast",
                    "CountStop", "FrameFreqStop")
  DF
}

# ── Helper functions for counting ────────────────────────────────────────
InFrameCounts <- function(x, A) {
  N <- A[x]
  sum(N[[1]][1, 2], N[[1]][1, 4], N[[1]][2, 4])
}

TotalCounts <- function(x, A) {
  N <- A[x]
  sum(N[[1]][, c(2, 4)])
}

InFrameSites <- function(x, A) {
  N <- A[x]
  sum(N[[1]][1, 3], N[[1]][1, 5], N[[1]][2, 5])
}

TotalSites <- function(x, A) {
  N <- A[x]
  sum(N[[1]][, c(3, 5)])
}

# ── Process each AA length ───────────────────────────────────────────────
cat("\n--- Processing uORFs by length ---\n")

Length <- as.data.frame(table(uORF_length), stringsAsFactors = FALSE)
Length$aa <- (as.numeric(Length$uORF_length) / 3) - 1
AA_vec <- Length$aa
cat(sprintf("  Unique AA lengths: %d (range: %d - %d)\n",
            length(AA_vec), min(AA_vec), max(AA_vec)))

dir.create("cips_Rdata", showWarnings = FALSE)

plan(multisession)

# Process each AA length
for (i in AA_vec) {
  gc()
  cat(sprintf("  AA = %d ... ", i))
  fiveUTR_ORFs_xaa <- fiveUTR_ORFs[uORF_length == (i + 1) * 3]
  n_orf <- length(fiveUTR_ORFs_xaa)
  cat(sprintf("%d ORFs\n", n_orf))

  if (n_orf == 0) next

  B <- future_lapply(seq_len(n_orf), function(x) XAA(x, y = i))
  names(B) <- names(fiveUTR_ORFs_xaa)
  save(B, file = paste0("cips_Rdata/AA_", i, ".RData"))
}

# ── Filter translated uORFs ──────────────────────────────────────────────
cat("\n--- Filtering translated uORFs ---\n")

InFrameCountSites <- function(z) {
  fiveUTR_ORFs_xaa <- fiveUTR_ORFs[uORF_length == (z + 1) * 3]
  load(paste0("cips_Rdata/AA_", z, ".RData"))

  if (length(B) == 0) return(data.frame())

  gene_id <- substr(names(B), 1, 9)
  InFrameC_AA <- unlist(lapply(seq_len(length(B)), function(x) InFrameCounts(x, A = B)))
  TotalC_AA   <- unlist(lapply(seq_len(length(B)), function(x) TotalCounts(x, A = B)))
  InFrameS_AA <- unlist(lapply(seq_len(length(B)), function(x) InFrameSites(x, A = B)))
  TotalS_AA   <- unlist(lapply(seq_len(length(B)), function(x) TotalSites(x, A = B)))
  STRAND <- sapply(seq_len(length(B)), function(x) {
    as.character(strand(fiveUTR_ORFs_xaa[names(B)[x]])[[1]][1])
  })

  INFO <- data.frame(
    ORF_id       = names(B),
    gene_id      = gene_id,
    InFrameC     = InFrameC_AA,
    TotalC       = TotalC_AA,
    InFrameC_perc = round(InFrameC_AA / TotalC_AA * 100, 1),
    InFrameS     = InFrameS_AA,
    TotalS       = TotalS_AA,
    InFrameS_perc = round(InFrameS_AA / TotalS_AA * 100, 1),
    Strand       = STRAND,
    AA           = z,
    stringsAsFactors = FALSE
  )
  INFO[is.na(INFO)] <- 0

  if (z == 1) {
    # Minimum uORFs: only need InFrameC and InFrameC_perc
    INFO2 <- INFO %>% filter(InFrameC >= min_inframe_counts,
                              InFrameC_perc > min_inframe_perc)
  } else {
    # Longer uORFs: also check InFrameS_perc
    INFO2 <- INFO %>% filter(InFrameC >= min_inframe_counts,
                              InFrameC_perc > min_inframe_perc,
                              InFrameS_perc > min_psite_perc)
  }
  INFO2
}

results_list <- list()
for (i in AA_vec) {
  res <- tryCatch(InFrameCountSites(z = i), error = function(e) data.frame())
  if (nrow(res) > 0) {
    results_list[[as.character(i)]] <- res
    cat(sprintf("  AA=%d: %d translated uORFs\n", i, nrow(res)))
  }
}

all_results <- do.call("rbind", results_list)
cat(sprintf("\n  Total translated uORFs: %d\n", nrow(all_results)))

if (nrow(all_results) == 0) {
  cat("No translated uORFs found with current thresholds. Exiting.\n")
  quit(save = "no", status = 0)
}

# ── Separate minimum uORFs (1aa) and longer uORFs ────────────────────────
minimum_uORF <- all_results %>% filter(AA == 1)
tuORF_list   <- all_results %>% filter(AA >= 2)

cat(sprintf("  Minimum uORFs (1aa): %d\n", nrow(minimum_uORF)))
cat(sprintf("  Longer uORFs (>=2aa): %d\n", nrow(tuORF_list)))

# ── Add coordinates for minimum uORFs ────────────────────────────────────
if (nrow(minimum_uORF) > 0) {
  cat("\n--- Computing minimum uORF coordinates ---\n")

  # Minimum uORFs computing (same code for AA==1 case)
  fiveUTR_ORFs_1aa <- fiveUTR_ORFs[uORF_length == (1 + 1) * 3]

  minimum_uORF$tx_id <- sub("_.*", "", minimum_uORF$ORF_id)
  minimum_uORF$chr <- substr(minimum_uORF$ORF_id, 3, 3)

  orf_ids <- minimum_uORF$ORF_id
  tx_ids <- gsub("\\_.*", "", orf_ids)

  # Transcript coordinates
  uORF_Tx_START <- function(x) {
    TxCord <- mapToTranscripts(fiveUTR_ORFs_1aa[[orf_ids[x]]],
                                exonByTx[tx_ids[x]], ignore.strand = FALSE)
    start(TxCord[1])
  }
  uORF_Tx_END <- function(x) {
    TxCord <- mapToTranscripts(fiveUTR_ORFs_1aa[[orf_ids[x]]],
                                exonByTx[tx_ids[x]], ignore.strand = FALSE)
    end(TxCord[length(TxCord)])
  }

  tx_start <- sapply(seq_len(length(orf_ids)), uORF_Tx_START)
  tx_end   <- sapply(seq_len(length(orf_ids)), uORF_Tx_END)
  minimum_uORF$uORF_Tx_range <- paste0(tx_start, "-", tx_end)

  # Genomic coordinates
  uORF_Ge_START <- function(x) {
    s <- as.character(strand(fiveUTR_ORFs[orf_ids[x]])[[1]])[1]
    ifelse(s == "+",
           start(unlist(fiveUTR_ORFs[orf_ids[x]]))[1],
           end(unlist(fiveUTR_ORFs[orf_ids[x]]))[1])
  }
  uORF_Ge_END <- function(x) {
    s <- as.character(strand(fiveUTR_ORFs[orf_ids[x]])[[1]])[1]
    ifelse(s == "+",
           max(end(unlist(fiveUTR_ORFs[orf_ids[x]]))),
           min(start(unlist(fiveUTR_ORFs[orf_ids[x]]))))
  }

  ge_start <- sapply(seq_len(length(orf_ids)), uORF_Ge_START)
  ge_end   <- sapply(seq_len(length(orf_ids)), uORF_Ge_END)
  minimum_uORF$uORF_Ge_range <- paste0(ge_start, "-", ge_end)

  # Peptide sequence
  minimum_uORF$peptide_seq <- "M*"
}

# ── Add coordinates for longer uORFs ─────────────────────────────────────
if (nrow(tuORF_list) > 0) {
  cat("\n--- Computing longer uORF coordinates and sequences ---\n")
  tuORF_df <- data.frame()

  for (aa_val in unique(tuORF_list$AA)) {
    fiveUTR_ORFs_xaa <- fiveUTR_ORFs[uORF_length == (aa_val + 1) * 3]
    sub_df <- tuORF_list %>% filter(AA == aa_val)
    sub_df$tx_id <- sub("_.*", "", sub_df$ORF_id)
    sub_df$chr <- substr(sub_df$ORF_id, 3, 3)

    orf_ids <- sub_df$ORF_id
    tx_ids <- gsub("\\_.*", "", orf_ids)

    uORF_Tx_START <- function(x) {
      TxCord <- mapToTranscripts(fiveUTR_ORFs_xaa[[orf_ids[x]]],
                                  exonByTx[tx_ids[x]], ignore.strand = FALSE)
      start(TxCord[1])
    }
    uORF_Tx_END <- function(x) {
      TxCord <- mapToTranscripts(fiveUTR_ORFs_xaa[[orf_ids[x]]],
                                  exonByTx[tx_ids[x]], ignore.strand = FALSE)
      end(TxCord[length(TxCord)])
    }

    tx_start <- sapply(seq_len(length(orf_ids)), uORF_Tx_START)
    tx_end   <- sapply(seq_len(length(orf_ids)), uORF_Tx_END)
    sub_df$uORF_Tx_range <- paste0(tx_start, "-", tx_end)

    uORF_Ge_START <- function(x) {
      s <- as.character(strand(fiveUTR_ORFs[orf_ids[x]])[[1]])[1]
      ifelse(s == "+",
             start(unlist(fiveUTR_ORFs[orf_ids[x]]))[1],
             end(unlist(fiveUTR_ORFs[orf_ids[x]]))[1])
    }
    uORF_Ge_END <- function(x) {
      s <- as.character(strand(fiveUTR_ORFs[orf_ids[x]])[[1]])[1]
      ifelse(s == "+",
             max(end(unlist(fiveUTR_ORFs[orf_ids[x]]))),
             min(start(unlist(fiveUTR_ORFs[orf_ids[x]]))))
    }

    ge_start <- sapply(seq_len(length(orf_ids)), uORF_Ge_START)
    ge_end   <- sapply(seq_len(length(orf_ids)), uORF_Ge_END)
    sub_df$uORF_Ge_range <- paste0(ge_start, "-", ge_end)

    tuORF_df <- rbind(tuORF_df, sub_df)
  }

  # Add peptide sequences
  fiveUTR_ORFs_peptide_df <- translate(fiveUTR_ORFs_seq) %>% as.data.frame()
  colnames(fiveUTR_ORFs_peptide_df) <- c("peptide_seq")
  fiveUTR_ORFs_peptide_df$ORF_id <- rownames(fiveUTR_ORFs_peptide_df)
  tuORF_df <- left_join(tuORF_df, fiveUTR_ORFs_peptide_df, by = "ORF_id")
}

# ── Deduplication: merge overlapping uORFs ───────────────────────────────
if (nrow(tuORF_list) > 0) {
  cat("\n--- Deduplicating overlapping uORFs ---\n")

  # Group by genomic range, merge ORF IDs
  tuORF_dedup <- tuORF_df %>%
    group_by(uORF_Ge_range, peptide_seq) %>%
    mutate(ORF_id_combine = toString(ORF_id)) %>%
    as.data.frame()

  # Also create a version keeping duplicates
  tuORF_with_dup <- tuORF_df %>%
    group_by(uORF_Ge_range) %>%
    mutate(ORF_id_combine = toString(ORF_id)) %>%
    as.data.frame()

  cat(sprintf("  Before dedup: %d uORFs\n", nrow(tuORF_df)))
  cat(sprintf("  After dedup:  %d unique genomic ranges\n",
              length(unique(tuORF_dedup$uORF_Ge_range))))

  # ── Find overlapping uORFs within same gene ───────────────────────────
  cat("\n--- Finding uORF overlaps ---\n")
  dup_genes <- unique(tuORF_with_dup$gene_id[duplicated(tuORF_with_dup$gene_id)])
  cat(sprintf("  Genes with multiple uORFs: %d\n", length(dup_genes)))

  if (length(dup_genes) > 0) {
    overlaps_list <- list()
    uORF_Overlaps_sum_ORF_id <- list()
    idx <- 0

    for (g in dup_genes) {
      g_df <- tuORF_with_dup %>% filter(gene_id == g)
      orf_ids <- g_df$ORF_id
      if (length(orf_ids) < 2) next

      co <- countOverlaps(fiveUTR_ORFs[orf_ids]) - 1
      sco <- sum(co)
      if (sco > 0) {
        idx <- idx + 1
        overlaps_list[[idx]] <- data.frame(
          gene_id = g,
          overlap_count = sco,
          stringsAsFactors = FALSE
        )
        uORF_Overlaps_sum_ORF_id[[idx]] <- names(co)[which(co > 0)]
      }
    }

    if (length(uORF_Overlaps_sum_ORF_id) > 0) {
      # Remove NULLs
      uORF_Overlaps_sum_ORF_id <- Filter(Negate(is.null), uORF_Overlaps_sum_ORF_id)

      # Build overlap dataframe
      uorf_len <- sapply(seq_along(uORF_Overlaps_sum_ORF_id),
                         function(x) length(uORF_Overlaps_sum_ORF_id[[x]]))

      overlap_dfs <- list()
      for (n in unique(uorf_len)) {
        idx_n <- which(uorf_len == n)
        if (length(idx_n) == 0) next
        mat <- do.call(rbind, uORF_Overlaps_sum_ORF_id[idx_n])
        df <- as.data.frame(mat, stringsAsFactors = FALSE)
        df$gene_id <- sub("\\..*", "", df$V1)
        df$tx_id   <- sub("_.*", "", df$V1)
        overlap_cols <- paste0("V", seq_len(n))
        df$uORFs_overlapped <- apply(df[, overlap_cols, drop = FALSE], 1,
                                      paste, collapse = ", ")
        overlap_dfs[[as.character(n)]] <- df
      }

      if (length(overlap_dfs) > 0) {
        overlap_final <- do.call("rbind", overlap_dfs)
        overlap_col_names <- grep("^V\\d+$", names(overlap_final), value = TRUE)
        keep_cols <- setdiff(names(overlap_final), overlap_col_names)
        keep_cols <- c("gene_id", "tx_id", "uORFs_overlapped")
        overlap_final <- overlap_final[, keep_cols, drop = FALSE]
      } else {
        overlap_final <- data.frame()
      }
    } else {
      overlap_final <- data.frame()
    }
  } else {
    overlap_final <- data.frame()
  }
}

# ── Annotate with gene descriptions (if provided) ────────────────────────
annotate_genes <- function(df) {
  if (is.na(gene_desc_xlsx) || !file.exists(gene_desc_xlsx) || nrow(df) == 0) {
    return(df)
  }
  cat("  Adding gene descriptions...\n")
  Gene_desc <- read.xlsx(gene_desc_xlsx, sheet = 1)
  df_out <- left_join(df, Gene_desc, by = "gene_id")
  if ("gene_model_type" %in% names(df_out)) {
    df_out$gene_model_type <- NULL
  }
  rownames(df_out) <- NULL
  df_out
}

# ── Write output ─────────────────────────────────────────────────────────
cat("\n--- Writing output ---\n")

# Minimum uORFs
if (nrow(minimum_uORF) > 0) {
  min_out <- annotate_genes(minimum_uORF)
  min_out_file <- paste0(output_prefix, "_minimum_uORF.xlsx")
  write.xlsx(min_out, file = min_out_file)
  cat(sprintf("  Minimum uORFs: %s (%d rows)\n", min_out_file, nrow(min_out)))
}

# All uORFs (deduplicated)
if (nrow(tuORF_list) > 0) {
  dedup_out <- annotate_genes(tuORF_dedup)
  dedup_file <- paste0(output_prefix, "_uORF_dedup.xlsx")
  write.xlsx(dedup_out, file = dedup_file)
  cat(sprintf("  Deduplicated uORFs: %s (%d rows)\n", dedup_file, nrow(dedup_out)))

  # All uORFs (with duplicates)
  dup_out <- annotate_genes(tuORF_with_dup)
  all_file <- paste0(output_prefix, "_uORF.xlsx")
  write.xlsx(dup_out, file = all_file)
  cat(sprintf("  All uORFs: %s (%d rows)\n", all_file, nrow(dup_out)))
}

# Overlaps
if (exists("overlap_final") && nrow(overlap_final) > 0) {
  overlap_file <- paste0(output_prefix, "_uORF_overlaps.xlsx")
  write.xlsx(overlap_final, file = overlap_file)
  cat(sprintf("  uORF overlaps: %s (%d rows)\n", overlap_file, nrow(overlap_final)))
}

# Cleanup temp Rdata
unlink("cips_Rdata", recursive = TRUE)

cat("\n=== CiPS analysis complete ===\n")
