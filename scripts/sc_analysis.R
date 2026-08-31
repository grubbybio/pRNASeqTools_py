# ── Single-cell RNA-seq analysis with Seurat (multi-group integration) ───
# Usage:
#   Rscript sc_analysis.R <prefix> <genome> <thread> <min_cells> <min_features>
#     <max_features> <pct_mt> <n_pcs> <n_clusters> <resolution>
#     <marker_min_pct> <marker_logfc> <pseudotime> <doublet_rate>
#     <integration> <tag_list> <sample_info>
#
#   integration: seurat | harmony | none
#   tag_list:    comma-separated tags that have count matrices
#   sample_info: tag:group,tag:group,...  (control=0, each -p = 1,2,...)
#
# Pipeline:
#   Phase 1 — Per-group QC & normalization (each group independent)
#   Phase 2 — Batch integration  (Seurat anchors OR Harmony OR none)
#   Phase 3 — Unified clustering, markers, DE, pseudotime

options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)
prefix           <- args[1]
genome           <- args[2]
thread           <- as.integer(args[3])
min_cells        <- as.integer(args[4])
min_features     <- as.integer(args[5])
max_features     <- as.integer(args[6])
pct_mt           <- as.numeric(args[7])
n_pcs            <- as.integer(args[8])
n_clusters       <- as.integer(args[9])
resolution       <- as.numeric(args[10])
marker_min_pct   <- as.numeric(args[11])
marker_logfc     <- as.numeric(args[12])
run_pseudotime   <- as.integer(args[13])
doublet_rate     <- as.numeric(args[14])
integration      <- args[15]          # seurat | harmony | none
tag_list_str     <- args[16]
sample_info_str  <- args[17]

# ── Parse tag / group info ────────────────────────────────────────────────
tags         <- strsplit(tag_list_str, ",")[[1]]
sample_parts <- strsplit(sample_info_str, ",")[[1]]
info_tags    <- sapply(sample_parts, function(x) strsplit(x, ":")[[1]][1])
groups_num   <- sapply(sample_parts, function(x) as.integer(strsplit(x, ":")[[1]][2]))

# Map each tag to its group id
tag_to_group <- setNames(groups_num, info_tags)

message("╔═══════════════════════════════════════════════════════╗")
message("║  Single-cell RNA-seq Analysis (Multi-Group Pipeline)   ║")
message("╚═══════════════════════════════════════════════════════╝")
message(paste0("  Samples (", length(tags), "): ", paste(tags, collapse = ", ")))
message(paste0("  Groups:  ", paste(groups_num, collapse = ", ")))
message(paste0("  Integration method: ", integration))
message("")

# ── Load packages ────────────────────────────────────────────────────────
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))

# Optional: Harmony
use_harmony <- FALSE
if (integration == "harmony") {
  if (requireNamespace("harmony", quietly = TRUE)) {
    suppressMessages(library(harmony))
    use_harmony <- TRUE
  } else {
    message("  ⚠ Harmony not installed — falling back to Seurat integration")
    integration <- "seurat"
  }
}

# ── Phase 1: Load, QC, per-group normalization ────────────────────────────
message("\n━━━━━━ Phase 1 — Load & per-group QC ━━━━━━\n")

seurat_list <- list()

for (tag in tags) {
  counts_file <- paste0(tag, "_counts.rds")
  group_id    <- tag_to_group[[tag]]
  group_label <- if (group_id == 0) "control" else paste0("treatment_", group_id)

  if (!file.exists(counts_file)) {
    message(paste("  ✗ Skipping", tag, "(counts file not found)"))
    next
  }

  message(paste("  Loading", tag, "→", group_label))
  counts <- readRDS(counts_file)

  seu <- CreateSeuratObject(
    counts  = counts,
    project = tag,
    min.cells    = min_cells,
    min.features = min_features
  )

  seu$group  <- group_id
  seu$group_label <- group_label
  seu$sample <- tag

  # ── QC: mitochondrial / chloroplast genes ──
  mt_genes <- grep("^MT-", rownames(seu), value = TRUE)
  if (length(mt_genes) == 0) {
    mt_genes <- grep("^ATMg", rownames(seu), value = TRUE)     # Arabidopsis mt
    if (length(mt_genes) == 0) {
      mt_genes <- grep("^mito", rownames(seu), value = TRUE, ignore.case = TRUE)
    }
  }

  if (length(mt_genes) > 0) {
    seu$pct_mt <- PercentageFeatureSet(seu, features = mt_genes)
  } else {
    # Chloroplast fallback for plants
    cp_genes <- grep("^ATC", rownames(seu), value = TRUE)
    if (length(cp_genes) > 0) {
      seu$pct_mt <- PercentageFeatureSet(seu, features = cp_genes)
    } else {
      seu$pct_mt <- 0
    }
  }

  # ── Apply per-cell filtering ──
  n_before <- ncol(seu)
  seu <- subset(seu,
                subset = nFeature_RNA > min_features &
                         nFeature_RNA < max_features &
                         pct_mt < pct_mt)
  n_after <- ncol(seu)
  message(paste("    QC:", n_before, "→", n_after, "cells"))

  # ── Per-group normalization ──
  seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                        scale.factor = 10000)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

  # ── Optional doublet removal ──
  if (doublet_rate > 0 && requireNamespace("DoubletFinder", quietly = TRUE)) {
    seu <- ScaleData(seu, verbose = FALSE) %>%
           RunPCA(npcs = n_pcs, verbose = FALSE) %>%
           FindNeighbors(dims = 1:n_pcs) %>%
           FindClusters(resolution = 0.5) %>%
           RunUMAP(dims = 1:n_pcs)

    sweep.res.list <- DoubletFinder::paramSweep_v3(seu, PCs = 1:n_pcs)
    sweep.pK       <- DoubletFinder::paramSweepStatFile(sweep.res.list)
    optimal.pK     <- sweep.pK[which.max(sweep.pK$metric.RX), "pK"]

    seu <- DoubletFinder::doubletFinder_v3(
      seu, PCs = 1:n_pcs, pK = as.numeric(optimal.pK),
      nExp = round(doublet_rate * ncol(seu))
    )
    n_bef <- ncol(seu)
    seu  <- subset(seu, subset = DF.classifications_ == "Singlet")
    message(paste("    Doublet:", n_bef, "→", ncol(seu), "cells"))
  }

  seurat_list[[tag]] <- seu
}

if (length(seurat_list) == 0) stop("✗ No valid samples after loading!")
message(paste("\n  Loaded", length(seurat_list), "sample(s) for integration"))

# ── QC visualization (before integration) ─────────────────────────────────
combined_pre <- merge(seurat_list[[1]], y = seurat_list[-1])
pdf("qc_violin_pre_integration.pdf", width = 14, height = 6)
print(VlnPlot(combined_pre,
              features = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
              group.by = "group_label", ncol = 3))
dev.off()
rm(combined_pre); gc()


# ══════════════════════════════════════════════════════════════════════════
# Phase 2 — Batch integration
# ══════════════════════════════════════════════════════════════════════════

# Determine unique "batches" (each tag = 1 batch; or by group if user wants)
# We use sample (tag) as batch unit for finest-grained correction
batch_labels <- unique(sapply(seurat_list, function(s) unique(s$sample)))
n_batches    <- length(batch_labels)

message(paste("\n━━━━━━ Phase 2 — Batch integration (", integration, ") ━━━━━━", sep = ""))
message(paste("  Batches:", paste(batch_labels, collapse = ", ")))

integrated <- NULL

if (n_batches == 1 || integration == "none") {
  # ── No integration needed (single batch or disabled) ──
  message("  → Skipping integration (single batch or 'none' requested)")
  integrated <- merge(seurat_list[[1]], y = seurat_list[-1])
  integrated <- ScaleData(integrated, features = rownames(integrated),
                          verbose = FALSE)

} else if (integration == "seurat") {
  # ── Seurat anchor-based integration ──
  message("  → Running Seurat FindIntegrationAnchors + IntegrateData")

  # Step 1: Select features consistently variable across batches
  features <- SelectIntegrationFeatures(object.list = seurat_list,
                                        nfeatures = 3000)

  # Step 2: Find anchors (reciprocal PCA for speed on large datasets)
  anchors <- FindIntegrationAnchors(
    object.list = seurat_list,
    dims        = 1:n_pcs,
    features    = features,
    reduction   = "rpca"
  )

  # Step 3: Integrate
  integrated <- IntegrateData(anchorset = anchors, dims = 1:n_pcs,
                              features.to.integrate = rownames(seurat_list[[1]]))

  message(paste("  Integrated", ncol(integrated), "cells across",
                n_batches, "batches"))

} else if (integration == "harmony" && use_harmony) {
  # ── Harmony integration ──
  message("  → Running Harmony (via harmony::RunHarmony)")

  combined <- merge(seurat_list[[1]], y = seurat_list[-1])
  combined <- ScaleData(combined, features = rownames(combined),
                        verbose = FALSE) %>%
              RunPCA(npcs = n_pcs, verbose = FALSE)

  integrated <- RunHarmony(
    object       = combined,
    group.by.vars = "sample",    # batch = sample tag
    reduction    = "pca",
    dims.use     = 1:n_pcs,
    assay.use    = "RNA",
    plot.convergence = FALSE
  )

  # Rename harmony reduction so downstream code is consistent
  DefaultAssay(integrated) <- "RNA"
  message(paste("  Integrated", ncol(integrated), "cells via Harmony"))
}

# ── Save intermediate per-group objects for diagnostic plots ────────────
saveRDS(seurat_list, "per_group_objects.rds")
message("  Saved per-group objects → per_group_objects.rds")


# ══════════════════════════════════════════════════════════════════════════
# Phase 3 — Unified analysis on integrated data
# ══════════════════════════════════════════════════════════════════════════
message("\n━━━━━━ Phase 3 — Unified clustering & downstream ━━━━━━\n")

# Use the integrated assay if present, else RNA
if ("integrated" %in% Assays(integrated)) {
  DefaultAssay(integrated) <- "integrated"
  message("  Using 'integrated' assay for PCA/clustering")
} else if ("harmony" %in% names(integrated@reductions)) {
  DefaultAssay(integrated) <- "RNA"
  message("  Using Harmony-corrected PCA reduction")
} else {
  DefaultAssay(integrated) <- "RNA"
  message("  Using PCA from scaled RNA assay")
}

# ── Dimensionality reduction ──
integrated <- ScaleData(integrated, verbose = FALSE)

# Which reduction to use as input for clustering / UMAP?
if ("harmony" %in% names(integrated@reductions)) {
  # Harmony: use the harmony embeddings directly
  use_reduction <- "harmony"
  message(paste("  Using Harmony reduction (", ncol(integrated@reductions$harmony),
                "dims) as input"))

  integrated <- RunUMAP(integrated, reduction = "harmony",
                         dims = 1:ncol(integrated@reductions$harmony),
                         verbose = FALSE)
  integrated <- FindNeighbors(integrated, reduction = "harmony",
                               dims = 1:ncol(integrated@reductions$harmony))
} else {
  # Seurat integrated or no integration: run PCA + UMAP
  integrated <- RunPCA(integrated, npcs = n_pcs, verbose = FALSE)

  pdf("elbow_plot_integrated.pdf", width = 8, height = 6)
  print(ElbowPlot(integrated, ndims = n_pcs))
  dev.off()

  use_reduction <- "pca"
  integrated <- RunUMAP(integrated, dims = 1:n_pcs, verbose = FALSE)
  integrated <- FindNeighbors(integrated, dims = 1:n_pcs)
}

# ── Clustering ──
message(paste("\n  Clustering (resolution =", resolution, ")..."))
integrated <- FindClusters(integrated, resolution = resolution,
                            group.singletons = TRUE)

n_clusters_found <- length(unique(integrated$seurat_clusters))
message(paste("  Found", n_clusters_found, "clusters"))

# ── UMAP plots ──
pdf("umap_clusters_integrated.pdf", width = 10, height = 8)
print(DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE) + ggtitle("Integrated Clusters"))
dev.off()

pdf("umap_groups_integrated.pdf", width = 10, height = 8)
print(DimPlot(integrated, reduction = "umap", group.by = "group_label",
              label = TRUE, repel = TRUE) + ggtitle("Group Labels on Integrated UMAP"))
dev.off()

pdf("umap_samples_integrated.pdf", width = 10, height = 8)
print(DimPlot(integrated, reduction = "umap", group.by = "sample",
              label = FALSE) + ggtitle("Sample Origin (Batch) on Integrated UMAP"))
dev.off()

# Save UMAP coordinates
umap_coords <- Embeddings(integrated, "umap")
write.csv(umap_coords, "umap_coordinates_integrated.csv", quote = TRUE)

# ── Cluster markers (join layers first if integrated assay) ──
message("\n  Finding cluster marker genes...")

if ("integrated" %in% Assays(integrated)) {
  DefaultAssay(integrated) <- "integrated"
}

all_markers <- FindAllMarkers(
  integrated,
  only.pos           = TRUE,
  min.pct            = marker_min_pct,
  logfc.threshold    = marker_logfc
)

write.csv(all_markers, "cluster_markers.csv", row.names = FALSE)
message(paste("  Found", nrow(all_markers), "marker-cluster pairs"))

# ── Top markers heatmap ──
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)

if (nrow(top_markers) > 0) {
  # Use RNA assay for expression heatmap
  DefaultAssay(integrated) <- "RNA"
  pdf("top_markers_heatmap.pdf", width = 12, height = 10)
  print(DoHeatmap(integrated, features = top_markers$gene) + NoLegend())
  dev.off()

  pdf("feature_plots_top_markers.pdf", width = 12, height = 8)
  print(FeaturePlot(integrated,
                    features = head(top_markers$gene, 10),
                    reduction = "umap", ncol = 2))
  dev.off()
}

# ── Cluster × Group composition ─────────────────────────────────────────
message("\n  Cluster composition by group:")

cluster_summary <- integrated@meta.data %>%
  group_by(seurat_clusters, group_label) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = group_label, values_from = cell_count,
                     values_fill = list(cell_count = 0))

write.csv(cluster_summary, "cluster_group_composition.csv", row.names = FALSE)
print(cluster_summary)

# Also write per-cluster summary
cluster_summary_total <- integrated@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    total_cells       = n(),
    dominant_group    = names(which.max(table(group_label))),
    dominant_group_pct = round(
      max(table(group_label)) / n() * 100, 1),
    .groups = "drop"
  )
write.csv(cluster_summary_total, "cluster_summary.csv", row.names = FALSE)

# ── Cross-group DE per cluster ──────────────────────────────────────────
# For each cluster, do a statistical test of group differences
message("\n  Testing cross-group DE (cluster × group)...")

# Use RNA assay for DE tests
DefaultAssay(integrated) <- "RNA"

de_results <- data.frame()

for (cl in unique(integrated$seurat_clusters)) {
  cl_cells <- rownames(integrated@meta.data[
    integrated$seurat_clusters == cl, ])
  if (length(cl_cells) < 10) next

  # Subset to this cluster and test group differences
  cl_subset <- subset(integrated, cells = cl_cells)

  # If we have exactly 2 groups, use wilcox; otherwise use ANOVA
  group_ids <- unique(cl_subset$group)
  if (length(group_ids) >= 2) {
    try({
      # Find all markers with group as identity
      de <- FindAllMarkers(
        cl_subset,
        group.by         = "group",
        only.pos         = FALSE,
        min.pct          = 0.05,
        logfc.threshold  = 0,
        test.use         = "wilcox"
      )
      de$cluster <- cl
      de_results <- rbind(de_results, de)
    }, silent = TRUE)
  }
}

if (nrow(de_results) > 0) {
  write.csv(de_results, "cross_group_DE.csv", row.names = FALSE)
  message(paste("  Found", nrow(de_results), "cross-group DE entries"))
} else {
  message("  No cross-group DE (too few groups/cells per cluster)")
}


# ══════════════════════════════════════════════════════════════════════════
# Phase 3b — Pseudotime analysis (optional)
# ══════════════════════════════════════════════════════════════════════════

if (run_pseudotime == 1) {
  message("\n━━━━━━ Phase 3b — Pseudotime (Monocle3) ━━━━━━\n")

  if (requireNamespace("monocle3", quietly = TRUE)) {
    counts_data <- GetAssayData(integrated, assay = "RNA", slot = "counts")
    meta_data   <- integrated@meta.data

    cds <- monocle3::new_cell_data_set(
      expression_data = counts_data,
      cell_metadata   = meta_data
    )

    cds <- monocle3::preprocess_cds(cds, num_dim = n_pcs)
    cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
    cds <- monocle3::learn_graph(cds)
    cds <- monocle3::order_cells(cds)

    pseudotime <- monocle3::pseudotime(cds)
    integrated$pseudotime <- pseudotime

    saveRDS(cds, "monocle3_results.rds")

    # Pseudotime on UMAP
    pt_plot <- data.frame(
      UMAP_1      = umap_coords[, 1],
      UMAP_2      = umap_coords[, 2],
      pseudotime  = pseudotime,
      cluster     = integrated$seurat_clusters,
      group_label = integrated$group_label
    )

    pdf("pseudotime_umap.pdf", width = 10, height = 8)
    print(
      ggplot(pt_plot, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
        geom_point(size = 0.5) +
        scale_color_viridis_c() +
        theme_minimal() +
        labs(title = "Pseudotime on Integrated UMAP")
    )
    dev.off()

    # Pseudotime distribution by group
    pdf("pseudotime_by_group.pdf", width = 10, height = 6)
    print(
      ggplot(pt_plot, aes(x = pseudotime, fill = group_label)) +
        geom_density(alpha = 0.5) +
        theme_minimal() +
        labs(title = "Pseudotime distribution by group")
    )
    dev.off()

    write.csv(
      data.frame(cell = colnames(integrated),
                 pseudotime = pseudotime,
                 cluster    = integrated$seurat_clusters,
                 group_label = integrated$group_label,
                 sample     = integrated$sample),
      "pseudotime_values.csv", row.names = FALSE
    )
    message("  Pseudotime analysis done")
  } else {
    message("  ⚠ monocle3 not installed — skipping pseudotime")
  }
}


# ══════════════════════════════════════════════════════════════════════════
# Save final objects
# ══════════════════════════════════════════════════════════════════════════

message("\n━━━━━━ Saving final results ━━━━━━\n")

saveRDS(integrated, "seurat_integrated.rds")

cluster_assignments <- data.frame(
  cell        = colnames(integrated),
  cluster     = integrated$seurat_clusters,
  group       = integrated$group,
  group_label = integrated$group_label,
  sample      = integrated$sample
)
write.csv(cluster_assignments, "cluster_assignments.csv", row.names = FALSE)

expr_matrix <- GetAssayData(integrated, assay = "RNA", slot = "data")
saveRDS(expr_matrix, "normalized_expression_integrated.rds")


# ── Summary ────────────────────────────────────────────────────────────────
message("\n╔═══════════════════════════════════════════════════════╗")
message("║  Single-cell analysis completed                       ║")
message("╚═══════════════════════════════════════════════════════╝")
message(paste0("  Total cells after QC + integration: ", ncol(integrated)))
message(paste0("  Total genes:                        ", nrow(integrated)))
message(paste0("  Clusters found:                     ", n_clusters_found))
message(paste0("  Groups:                             ", length(unique(integrated$group_label))))
message(paste0("  Integration method:                 ", integration))
message("")
message("  Output files:")
message("    ── Core ──")
message("    • seurat_integrated.rds        : Integrated Seurat object")
message("    • cluster_markers.csv          : Marker genes per cluster")
message("    • cluster_assignments.csv      : Cell-to-cluster/group mapping")
message("    • cluster_summary.csv          : Per-cluster stats")
message("    • cluster_group_composition.csv: Cluster × Group contingency table")
message("    ── Cross-group DE ──")
message("    • cross_group_DE.csv           : DE genes per cluster across groups")
message("    ── QC & integration diagnostics ──")
message("    • qc_violin_pre_integration.pdf")
message("    • per_group_objects.rds       : Per-group Seurat objects (pre-integration)")
message("    ── Visualization ──")
message("    • umap_clusters_integrated.pdf")
message("    • umap_groups_integrated.pdf")
message("    • umap_samples_integrated.pdf")
message("    • elbow_plot_integrated.pdf")
message("    • top_markers_heatmap.pdf")
message("    • feature_plots_top_markers.pdf")

if (run_pseudotime == 1 && requireNamespace("monocle3", quietly = TRUE)) {
  message("    ── Pseudotime ──")
  message("    • monocle3_results.rds")
  message("    • pseudotime_umap.pdf")
  message("    • pseudotime_by_group.pdf")
  message("    • pseudotime_values.csv")
}

message("\n✓ Done.")
