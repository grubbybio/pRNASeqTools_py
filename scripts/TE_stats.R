#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# TE_stats.R — Statistical testing of Translation Efficiency (TE)
#
# Design: multiple conditions (control first), each with biological replicates.
# Per-gene log2(TE) values are paired across conditions (same gene = same unit).
#
# Method:
#   1. Build long-format table: gene_id | condition | replicate | sample_pair
#   2. Filter genes with enough valid data across conditions
#   3. Linear mixed-effects style paired analysis via lm with gene as blocking:
#        log2(TE) ~ condition + gene   (fixed effects, paired by gene)
#      followed by emmeans for estimated marginal means + pairwise contrasts
#      and agricolae::HSD.test for compact letter display.
#   4. Also report per-condition median/mean log2(TE) and per-pair Wilcoxon.
#
# Usage:
#   Rscript TE_stats.R <te_wide_tsv> <output_dir> <alpha> [group1=rep1 group2=rep2 ...]
#
#   te_wide_tsv: wide-format TE table produced by the pipeline
#                columns: gene_id, {pair1}_log2te, {pair2}_log2te, ...
#   group specs: e.g. "ctrl__ctrl=2 trt__trt=2" — condition_name=n_replicates,
#                where each replicate maps to a consecutive set of pair columns.
#
# Outputs in output_dir:
#   TE_emmeans.tsv       - estimated marginal means per condition
#   TE_contrasts.tsv     - pairwise contrasts with p-values
#   TE_letters.tsv       - compact letter display (shared letter = not sig)
#   TE_summary_stats.tsv - per-condition summary stats
# ─────────────────────────────────────────────────────────────────────────────

options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript TE_stats.R <te_wide_tsv> <output_dir> <alpha> [group=n ...]")
}

input_file <- args[1]
output_dir <- args[2]
alpha <- as.numeric(args[3])
group_args <- args[-c(1:3)]  # remaining: "condition=n_replicates"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load libraries ───────────────────────────────────────────────────────────
suppressMessages({
  library(emmeans)
  library(multcomp)
})

message("Reading input: ", input_file)
dat <- read.delim(input_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!"gene_id" %in% colnames(dat)) {
  stop("Input file must have a 'gene_id' column")
}
gene_ids <- dat$gene_id
dat$gene_id <- NULL  # keep only value columns

# Identify log2te columns (e.g., "pairname_log2te")
value_cols <- colnames(dat)
log2te_cols <- grep("_log2te$", value_cols, value = TRUE)

if (length(log2te_cols) < 2) {
  message("ERROR: need at least 2 sample pairs to compare, found ",
          length(log2te_cols))
  quit(save = "no", status = 0)
}

# ── Parse group structure ────────────────────────────────────────────────────
# Each group spec is "condition_name=n_replicates". Column order in log2te_cols
# must match the order groups were given; replicates within a group map onto
# consecutive columns.
group_names <- character()
group_reps <- integer()
for (g in group_args) {
  parts <- strsplit(g, "=", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    stop("Invalid group spec '", g, "' - expected 'condition=n'")
  }
  group_names <- c(group_names, parts[1])
  group_reps <- c(group_reps, as.integer(parts[2]))
}

total_pairs <- sum(group_reps)
if (total_pairs != length(log2te_cols)) {
  message("WARNING: number of log2te columns (", length(log2te_cols),
          ") does not match sum of replicates (", total_pairs, ").")
  n_use <- min(total_pairs, length(log2te_cols))
} else {
  n_use <- total_pairs
}

# Map column index -> (condition, replicate_within_condition)
col_group <- character(n_use)
col_rep <- integer(n_use)
idx <- 1
for (g_idx in seq_along(group_names)) {
  for (r in seq_len(group_reps[g_idx])) {
    if (idx > n_use) break
    col_group[idx] <- group_names[g_idx]
    col_rep[idx] <- r
    idx <- idx + 1
  }
}
names(col_group) <- names(col_rep) <- log2te_cols[seq_len(n_use)]

# ── Build long-format dataframe ──────────────────────────────────────────────
message("Building long-format table...")
long_list <- lapply(seq_len(n_use), function(k) {
  cname <- log2te_cols[k]
  vals <- suppressWarnings(as.numeric(dat[[cname]]))
  data.frame(
    gene_id = gene_ids,
    condition = rep(col_group[k], length(gene_ids)),
    replicate = rep(col_rep[k], length(gene_ids)),
    log2te = vals,
    stringsAsFactors = FALSE
  )
})
long_dat <- do.call(rbind, long_list)
long_dat <- long_dat[!is.na(long_dat$log2te) & is.finite(long_dat$log2te), ]

# Require genes present in at least 2 conditions to be usable
gene_cond_counts <- tapply(long_dat$condition, long_dat$gene_id,
                           function(x) length(unique(x)))
usable_genes <- names(gene_cond_counts)[gene_cond_counts >= 2]
long_dat <- long_dat[long_dat$gene_id %in% usable_genes, ]
message("Genes used after filtering: ", length(usable_genes))

n_conditions <- length(unique(long_dat$condition))
if (n_conditions < 2) {
  message("ERROR: fewer than 2 conditions with data; no statistical test possible.")
  quit(save = "no", status = 0)
}

long_dat$condition <- factor(long_dat$condition, levels = unique(group_names))

# ── Summary statistics per condition ────────────────────────────────────────
summary_by_cond <- aggregate(log2te ~ condition, data = long_dat,
                             FUN = function(x)
                               c(median = round(median(x), 4),
                                 mean = round(mean(x), 4),
                                 sd = round(sd(x), 4),
                                 n_genes = length(x)))
write.table(summary_by_cond,
            file.path(output_dir, "TE_summary_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Paired analysis: log2(TE) ~ condition + gene ────────────────────────────
message("Fitting linear model: log2te ~ condition + gene (paired design)...")
model_fit <- tryCatch({
  fit <- lm(log2te ~ condition + gene_id, data = long_dat)
  fit
}, error = function(e) {
  message("LM failed: ", e$message,
          " - falling back to condition-only model.")
  return(lm(log2te ~ condition, data = long_dat))
})

# ANOVA table using car::Anova (Type II) for robustness to imbalance
anova_tbl <- tryCatch({
  car::Anova(model_fit, type = "II")
}, error = function(e) {
  anova(model_fit)
})
message("\nANOVA:")
print(anova_tbl)

cond_row_idx <- which(grepl("^condition", rownames(anova_tbl)))
p_anova <- if (length(cond_row_idx)) anova_tbl[[ncol(anova_tbl)]][cond_row_idx[1]] else NA
f_anova <- NA
if (length(cond_row_idx)) {
  # Type II table columns: Sum Sq, Df, F value, Pr(>F); base R: Df, Sum Sq, Mean Sq, F, Pr(>F)
  f_col_idx <- grep("F", colnames(anova_tbl))[1]
  if (!is.na(f_col_idx)) f_anova <- anova_tbl[cond_row_idx[1], f_col_idx]
}

# ── Estimated marginal means + Tukey-adjusted contrasts via emmeans ─────────
message("\nComputing estimated marginal means & Tukey contrasts (emmeans)...")
emm_result <- tryCatch({
  emm <- emmeans::emmeans(model_fit, ~ condition)
  contrast_res <- pairs(emm, adjust = "tukey")
  list(emm = emm, contrast = contrast_res)
}, error = function(e) {
  message("emmeans on full model failed: ", e$message,
          " - trying condition-only model...")
  fit_simple <- lm(log2te ~ condition, data = long_dat)
  emm <- emmeans::emmeans(fit_simple, ~ condition)
  contrast_res <- pairs(emm, adjust = "tukey")
  list(emm = emm, contrast = contrast_res)
})

emm_df <- as.data.frame(emm_result$emm)
contr_df <- as.data.frame(emm_result$contrast)
write.table(emm_df[, intersect(c("condition", "emmean", "SE", "lower.CL", "upper.CL"), colnames(emm_df))],
            file.path(output_dir, "TE_emmeans.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(contr_df, file.path(output_dir, "TE_contrasts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("\nEstimated marginal means:")
print(emm_df)
message("\nPairwise contrasts (Tukey-adjusted):")
print(contr_df)

# ── Compact Letter Display via multcomp::cld on emmeans ─────────────────────
message("\nComputing compact letter display (multcomp::cld)...")
letters_result <- tryCatch({
  cld_res <- multcomp::cld(emm_result$emm, alpha = alpha,
                           Letters = letters, sorted = TRUE)
  cld_df <- as.data.frame(cld_res)
  # .group column holds the letters (with spaces)
  data.frame(condition = cld_df$condition,
             letters = gsub(" ", "", cld_df$.group),
             stringsAsFactors = FALSE)
}, error = function(e) {
  message("multcomp::cld failed: ", e$message,
          " - deriving letters from emmeans contrasts directly.")
  # Fallback: greedy letter assignment from significant contrasts
  groups <- levels(long_dat$condition)
  sig_mat <- matrix(FALSE, length(groups), length(groups),
                    dimnames = list(groups, groups))
  for (i in seq_len(nrow(contr_df))) {
    parts <- trimws(strsplit(as.character(contr_df$contrast[i]), "-",
                             fixed = TRUE)[[1]])
    pval_i <- contr_df$p.value[i]
    if (length(parts) == 2 && is.finite(pval_i) && pval_i < alpha &&
        parts[1] %in% groups && parts[2] %in% groups) {
      sig_mat[parts[1], parts[2]] <- TRUE
      sig_mat[parts[2], parts[1]] <- TRUE
    }
  }
  # Assign letters greedily in descending emmean order
  ordered_groups <- emm_df$condition[order(-emm_df$emmean)]
  assigned <- setNames(vector("list", length(ordered_groups)), ordered_groups)
  n_letters <- 0
  for (g in ordered_groups) {
    for (li in seq_len(max(n_letters, 1))) {
      members <- names(which(vapply(assigned, function(x) li %in% x, logical(1))))
      if (!any(vapply(members, function(m) sig_mat[g, m], logical(1)))) {
        assigned[[g]] <- c(assigned[[g]], li)
        break
      }
    }
    if (length(assigned[[g]]) == 0) {
      n_letters <- n_letters + 1
      assigned[[g]] <- n_letters
    }
  }
  data.frame(condition = names(assigned),
             letters = vapply(assigned, function(x)
               paste(letters[sort(unique(unlist(x)))], collapse = ""),
             character(1)),
             stringsAsFactors = FALSE)
})
write.table(letters_result, file.path(output_dir, "TE_letters.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("\nCompact letter display:")
print(letters_result)

# ── Report overall pipeline statistics ───────────────────────────────────────
message("\n===== TE Statistical Test Summary =====")
message("Conditions tested: ", n_conditions)
message("Genes analyzed:   ", length(usable_genes))
if (is.finite(p_anova)) {
  message(sprintf("ANOVA (condition effect): F = %.3f, p = %.3e",
                  f_anova, p_anova))
}
message("Letter display written to: ",
        file.path(output_dir, "TE_letters.tsv"))
message("Full results directory: ", output_dir)
message("Done.")
