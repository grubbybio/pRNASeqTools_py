#!/usr/bin/env Rscript
#
# ChIPseeker Peak Annotation & clusterProfiler GO Enrichment
#
# Processes a SINGLE peak file — designed to be called once per peak set,
# making it reusable by both tf (5 peak sets) and chip (1 peak set) modules.
#
# Usage:
#   Rscript chipseeker.R <genome> <prefix> <peak_file> <peak_name> [tss_distance]
#
# Example:
#   Rscript chipseeker.R ath /path/to/pkg WT_peaks.narrowPeak WT_peaks
#   Rscript chipseeker.R ath /path/to/pkg WT_peaks.narrowPeak WT_peaks 2000
#
# Arguments:
#   genome       - genome identifier (e.g. ath → uses org.At.tair.db)
#   prefix       - project prefix path (contains reference/ directory)
#   peak_file    - peak file path (BED/narrowPeak format)
#   peak_name    - label used for output file naming
#   tss_distance - distance upstream of TSS for annotation (default: 3000)
#
# Output files (written to CWD):
#   - {peak_name}_annotation.txt       — ChIPseeker annotation table
#   - {peak_name}_go_enrichment.txt    — GO enrichment results (if OrgDb available)
#   - {peak_name}_annotation_pie.pdf   — Annotation pie chart
#   - {peak_name}_go_dotplot.pdf       — GO enrichment dotplot

suppressPackageStartupMessages({
    library(ChIPseeker)
    library(GenomicFeatures)
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
})

# ── Parse command-line arguments ──────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript chipseeker.R <genome> <prefix> <peak_file> <peak_name> [tss_distance]")
}

genome       <- args[1]
prefix       <- args[2]
peak_file    <- args[3]
peak_name    <- args[4]
tss_distance <- ifelse(length(args) >= 5, as.numeric(args[5]), 3000)

# ── Check peak file exists ────────────────────────────────────────────────
if (!file.exists(peak_file)) {
    stop("Peak file not found: ", peak_file)
}

# ── Resolve GFF path ──────────────────────────────────────────────────────
gff_path <- file.path(prefix, "reference", paste0(genome, "_genes.gff"))
if (!file.exists(gff_path)) {
    stop("GFF file not found: ", gff_path)
}

# ── Build TxDb from GFF ──────────────────────────────────────────────────
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("ChIPseeker Peak Annotation —", peak_name, "\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("  Genome   :", genome, "\n")
cat("  GFF      :", gff_path, "\n")
cat("  Peak file:", peak_file, "\n")
cat("  Peak name:", peak_name, "\n")
cat("  TSS region: [-", tss_distance, ", 0]\n\n", sep = "")

cat("Building TxDb from GFF...\n")
txdb <- makeTxDbFromGFF(gff_path)
cat("  TxDb created successfully.\n\n")

# ── Annotate peaks ────────────────────────────────────────────────────────
cat("Annotating peaks...\n")
peaks <- readPeakFile(peak_file, header = FALSE)

# Check if this is a narrowPeak file (has score column)
peak_colnames <- if (grepl("\\.narrowPeak$", peak_file)) {
    c("chr", "start", "end", "name", "score", "strand",
      "signalValue", "pValue", "qValue", "peak")
} else {
    c("chr", "start", "end")
}
colnames(peaks)[1:length(peak_colnames)] <- peak_colnames[1:ncol(peaks)]

annotation <- annotatePeak(
    peaks,
    TxDb          = txdb,
    tssRegion     = c(-tss_distance, 0),
    addFlankGeneInfo = TRUE,
    flankDistance = tss_distance,
    sameStrand    = FALSE,
    ignoreOverlap = FALSE
)

# ── Write annotation table ───────────────────────────────────────────────
anno_df <- as.data.frame(annotation)
out_txt <- paste0(peak_name, "_annotation.txt")
write.table(anno_df, file = out_txt, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
cat("  ->", out_txt, "\n")

# ── Annotation pie chart ──────────────────────────────────────────────────
cat("\nGenerating annotation pie chart...\n")
out_pie <- paste0(peak_name, "_annotation_pie.pdf")
pdf(out_pie, width = 10, height = 8)
tryCatch(
    plotAnnoPie(annotation),
    error = function(e) cat("  Warning: could not generate pie chart:", e$message, "\n")
)
dev.off()
cat("  ->", out_pie, "\n")

# ── Extract gene IDs ──────────────────────────────────────────────────────
gids <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
cat("\n  Genes annotated:", length(gids), "\n")

if (length(gids) == 0) {
    cat("  No genes found in annotation, skipping GO enrichment.\n")
    q()
}

# ── GO Enrichment Analysis ────────────────────────────────────────────────
# Static OrgDb mapping (genomes whose GFF IDs match official Bioconductor OrgDb)
org_db_map <- list(
    ath = "org.At.tair.db"
)

# Species scientific names for AnnotationHub fallback (Ensembl Plants species)
species_name_map <- list(
    osa      = "Oryza sativa",
    b73      = "Zea mays",
    zma      = "Zea mays",
    w22      = "Zea mays",
    gma      = "Glycine max",
    zh13     = "Glycine max",
    hx3      = "Glycine max",
    w82njau  = "Glycine max",
    bra      = "Brassica rapa",
    ahy      = "Arachis hypogaea",
    cla      = "Citrullus lanatus",
    cme      = "Cucumis melo",
    csa      = "Cucumis sativus",
    mpo      = "Marchantia polymorpha",
    sly      = "Solanum lycopersicum",
    tae      = "Triticum aestivum"
)

# Attempt to resolve OrgDb
org_obj <- NULL
org_key <- "ENTREZID"

# Step 1: Static mapping (official pre-installed OrgDb packages)
static_pkg <- org_db_map[[genome]]
if (!is.null(static_pkg)) {
    cat("\nLoading organism database:", static_pkg, "...\n")
    if (suppressPackageStartupMessages(require(static_pkg, character.only = TRUE))) {
        org_obj <- get(static_pkg)
        if (genome == "ath") org_key <- "TAIR"
    } else {
        cat("  Warning: package", static_pkg, "not installed.\n")
    }
}

# Step 2: AnnotationHub fallback (for species in Ensembl Plants)
if (is.null(org_obj)) {
    species <- species_name_map[[genome]]
    if (!is.null(species)) {
        cat("\nSearching AnnotationHub for '", species, "'...\n", sep = "")
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
        cat("\n  No species name mapping for genome '", genome, "'.\n", sep = "")
    }
}

# Step 3: Skip if no OrgDb available
if (is.null(org_obj)) {
    cat("\nNo OrgDb available for genome '", genome, "', skipping GO enrichment.\n", sep = "")
}

# ── Run enrichGO ─────────────────────────────────────────────────────────
if (!is.null(org_obj)) {
    cat("  Running enrichGO for:", peak_name, "(", length(gids), "genes )...\n")
    cat("  Using keyType:", org_key, "\n")
    ego <- tryCatch(
        enrichGO(
            gene          = gids,
            OrgDb         = org_obj,
            keyType       = org_key,
            ont           = "ALL",
            pvalueCutoff  = 0.05,
            qvalueCutoff  = 0.1,
            readable      = TRUE
        ),
        error = function(e) {
            cat("    enrichGO failed:", e$message, "\n")
            return(NULL)
        }
    )
    if (!is.null(ego) && nrow(ego@result) > 0) {
        out_go <- paste0(peak_name, "_go_enrichment.txt")
        write.table(as.data.frame(ego), file = out_go,
                    sep = "\t", quote = FALSE, row.names = FALSE)
        cat("    ->", out_go, "\n")

        # GO dotplot
        out_dot <- paste0(peak_name, "_go_dotplot.pdf")
        cat("\nGenerating GO dotplot...\n")
        pdf(out_dot, width = 12, height = 8)
        tryCatch(
            print(dotplot(ego, showCategory = 10) +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                  ggtitle(paste("GO Enrichment —", peak_name))),
            error = function(e) cat("  Warning: could not generate dotplot:", e$message, "\n")
        )
        dev.off()
        cat("  ->", out_dot, "\n")
    } else {
        cat("    No significant GO terms found.\n")
    }
}

cat("\nChIPseeker annotation and GO enrichment completed for", peak_name, ".\n")
