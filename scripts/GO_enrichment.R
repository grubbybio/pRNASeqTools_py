#!/usr/bin/env Rscript
#
# GO Enrichment Analysis for DEGs (mRNA-seq / DESeq2 output)
#
# For a given comparison (e.g. "AlvsWT"), this script reads the
# upregulated / downregulated DEG CSV files produced by DEG.R,
# runs clusterProfiler::enrichGO, and exports:
#   - {comparison}.upregulated.go_enrichment.txt
#   - {comparison}.upregulated.go_dotplot.pdf
#   - {comparison}.downregulated.go_enrichment.txt
#   - {comparison}.downregulated.go_dotplot.pdf
#
# Usage:
#   Rscript GO_enrichment.R <genome> <prefix> <comparison_name>
#
# Example:
#   Rscript GO_enrichment.R ath /path/to/pkg AlvsWT

options(warn = -1)

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
})

# ── Parse arguments ─────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript GO_enrichment.R <genome> <prefix> <comparison_name>")
}

genome        <- args[1]
prefix        <- args[2]
comparison    <- args[3]   # e.g. "AlvsWT"

cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("GO Enrichment Analysis — ", comparison, "\n", sep = "")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("  Genome     :", genome, "\n")
cat("  Prefix     :", prefix, "\n")
cat("  Comparison :", comparison, "\n\n")

# ── Resolve OrgDb (mirror chipseeker.R logic) ───────────────────────────
org_db_map <- list(
    ath = "org.At.tair.db"
)
species_name_map <- list(
    osa     = "Oryza sativa",
    b73     = "Zea mays",
    zma     = "Zea mays",
    w22     = "Zea mays",
    gma     = "Glycine max",
    zh13    = "Glycine max",
    hx3     = "Glycine max",
    w82njau = "Glycine max",
    bra     = "Brassica rapa",
    ahy     = "Arachis hypogaea",
    cla     = "Citrullus lanatus",
    cme     = "Cucumis melo",
    csa     = "Cucumis sativus",
    mpo     = "Marchantia polymorpha",
    sly     = "Solanum lycopersicum",
    tae     = "Triticum aestivum"
)

org_obj <- NULL
org_key <- "ENTREZID"

static_pkg <- org_db_map[[genome]]
if (!is.null(static_pkg)) {
    cat("Loading organism database:", static_pkg, "...\n")
    if (suppressPackageStartupMessages(require(static_pkg, character.only = TRUE))) {
        org_obj <- get(static_pkg)
        if (genome == "ath") org_key <- "TAIR"
        cat("  -> OK (keyType =", org_key, ")\n")
    } else {
        cat("  Warning: package", static_pkg, "not installed.\n")
    }
}

if (is.null(org_obj)) {
    species <- species_name_map[[genome]]
    if (!is.null(species)) {
        cat("Searching AnnotationHub for '", species, "'...\n", sep = "")
        tryCatch({
            suppressPackageStartupMessages(library(AnnotationHub))
            ah <- AnnotationHub()
            ah_query <- query(ah, c("OrgDb", species))
            if (length(ah_query) > 0) {
                org_obj <- ah_query[[length(ah_query)]]
                cat("  -> Retrieved OrgDb via AnnotationHub for", species, "\n")
            } else {
                cat("  -> No OrgDb found for", species, "in AnnotationHub.\n")
            }
        }, error = function(e) {
            cat("  -> AnnotationHub query failed:", e$message, "\n")
        })
    } else {
        cat("  No species name mapping for genome '", genome,
            "' — skipping GO enrichment.\n", sep = "")
    }
}

if (is.null(org_obj)) {
    cat("\nWARNING: Could not resolve OrgDb for genome '", genome,
        "'. GO enrichment skipped.\n", sep = "")
    quit(status = 0, save = "no")
}

# ── Load DEG CSV files ──────────────────────────────────────────────────
load_deg <- function(fname) {
    if (!file.exists(fname)) return(NULL)
    df <- read.csv(fname, row.names = 1, check.names = FALSE)
    rownames(df)
}

genes_up   <- load_deg(paste0(comparison, ".total.upregulated.csv"))
genes_down <- load_deg(paste0(comparison, ".total.downregulated.csv"))

cat("  Up DEG     :", length(genes_up),   "genes\n")
cat("  Down DEG   :", length(genes_down), "genes\n\n")

# ── Helper: run enrichGO + save results + dotplot ──────────────────────
run_go <- function(gene_ids, label, org_obj, org_key) {
    if (is.null(gene_ids) || length(gene_ids) == 0) {
        cat("  [", label, "] No genes — skipped.\n", sep = "")
        return(invisible(NULL))
    }

    cat("  Running enrichGO for [", label, "] (", length(gene_ids),
        " genes)...\n", sep = "")

    ego <- tryCatch(
        enrichGO(
            gene         = gene_ids,
            OrgDb        = org_obj,
            keyType      = org_key,
            ont          = "ALL",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            readable     = TRUE
        ),
        error = function(e) {
            cat("    enrichGO failed:", e$message, "\n")
            return(NULL)
        }
    )

    if (is.null(ego) || nrow(ego@result) == 0) {
        cat("    No significant GO terms found.\n")
        return(invisible(NULL))
    }

    cat("    Found", nrow(ego@result), "significant GO terms.\n")

    # ── Save enrichment table ──
    out_txt <- paste0(comparison, ".", label, ".go_enrichment.txt")
    write.table(as.data.frame(ego), file = out_txt,
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("    ->", out_txt, "\n")

    # ── Dotplot (BP / CC / MF facets, top 10 each) ──
    out_dot <- paste0(comparison, ".", label, ".go_dotplot.pdf")
    pdf(out_dot, width = 14, height = 7)
    tryCatch(
        print(dotplot(ego, showCategory = 10, split = "ONTOLOGY") +
                  facet_grid(ONTOLOGY ~ ., scales = "free_y") +
                  ggtitle(paste0("GO Enrichment — ", comparison, " (", label, ")")) +
                  theme(axis.text.y = element_text(size = 8))),
        error = function(e) cat("    Warning: dotplot failed:", e$message, "\n")
    )
    dev.off()
    cat("    ->", out_dot, "\n")

    invisible(ego)
}

# ── Run GO enrichment ─────────────────────────────────────────────────────
for (item in list(
    list(genes_up,   "upregulated"),
    list(genes_down, "downregulated")
)) {
    run_go(item[[1]], item[[2]], org_obj, org_key)
    cat("\n")
}

cat("GO enrichment completed for", comparison, ".\n")
