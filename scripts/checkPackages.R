# ── Parse command-line arguments ────────────────────────────────────────────
# Usage:
#   Rscript checkPackages.R                        # install ALL R packages
#   Rscript checkPackages.R --packages pkg1,pkg2   # install only listed packages
args <- commandArgs(trailingOnly = TRUE)

pkg_list <- NULL
if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--packages" && i < length(args)) {
      pkg_list <- strsplit(args[i + 1], ",")[[1]]
      pkg_list <- trimws(pkg_list)
      break
    }
  }
}

if (is.null(pkg_list)) {
  message("\nChecking R packages (all)...")
} else {
  message(paste("\nChecking R packages:", paste(pkg_list, collapse = ", ")))
}

# ── Bootstrap prerequisites ────────────────────────────────────────────────
options(repos = structure(c(CRAN = "https://cloud.r-project.org")))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# ── Install dispatch ───────────────────────────────────────────────────────
# Each function receives a package name and returns invisibly.
# Returns FALSE (invisibly) if the package was already installed.

install_bioconductor <- function(pkg) {
  if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
    message(paste("Installing", pkg, "(Bioconductor)..."))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    message(paste(pkg, "already installed"))
  }
}

install_github <- function(pkg, repo, ref = NULL, extra = NULL) {
  if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
    message(paste("Installing", pkg, "(GitHub:", repo, ")..."))
    devtools::install_github(repo, ref = ref, dependencies = TRUE)
  } else {
    # Version check for packages that require a specific branch
    if (!is.null(ref)) {
      current_ref <- packageDescription(pkg)$GithubRef
      if (is.null(current_ref) || current_ref != ref) {
        message(paste("Updating", pkg, "to", ref, "branch..."))
        devtools::install_github(repo, ref = ref, dependencies = TRUE)
      } else {
        message(paste(pkg, "already installed (", ref, "branch)"))
      }
    } else {
      message(paste(pkg, "already installed"))
    }
  }
  # Install extra companion packages
  if (!is.null(extra)) {
    for (ep in extra) {
      if (!suppressMessages(require(ep, character.only = TRUE, quietly = TRUE))) {
        message(paste("Installing companion package:", ep, "..."))
        # uwot is on CRAN
        install.packages(ep)
      }
    }
  }
}

install_cran <- function(pkg) {
  if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
    message(paste("Installing", pkg, "(CRAN)..."))
    install.packages(pkg)
  } else {
    message(paste(pkg, "already installed"))
  }
}

# ── Single-package dispatcher ──────────────────────────────────────────────
install_one <- function(pkg) {
  switch(pkg,
    "dplyr"                   = install_cran("dplyr"),
    "DESeq2"                  = install_bioconductor("DESeq2"),
    "DMRcaller"               = install_bioconductor("DMRcaller"),
    "ORFik"                   = install_bioconductor("ORFik"),
    "pheatmap"                = install_bioconductor("pheatmap"),
    "ComplexHeatmap"          = install_bioconductor("ComplexHeatmap"),
    "RNAmodR.RiboMethSeq"     = install_bioconductor("RNAmodR.RiboMethSeq"),
    "emmeans"                 = install_cran("emmeans"),
    "car"                     = install_cran("car"),
    "agricolae"               = install_cran("agricolae"),
    "multcomp"                = install_cran("multcomp"),
    "ggpubr"                  = install_cran("ggpubr"),
    "NMF"                     = install_github("NMF", "renozao/NMF", ref = "devel"),
    "riboWaltz"               = install_github("riboWaltz", "LabTranslationalArchitectomics/riboWaltz"),
    "Seurat"                  = install_github("Seurat", "satijalab/seurat", extra = c("uwot")),
    {
      message(paste("Unknown package:", pkg, "- skipping"))
    }
  )
}

# ── Main ───────────────────────────────────────────────────────────────────
if (is.null(pkg_list)) {
  # Install all known packages
  all_pkgs <- c("dplyr", "DESeq2", "DMRcaller", "ORFik", "pheatmap",
                 "ComplexHeatmap", "RNAmodR.RiboMethSeq",
                 "emmeans", "car", "agricolae", "multcomp", "ggpubr",
                 "NMF", "riboWaltz", "Seurat")
  pkg_list <- all_pkgs
}

for (pkg in pkg_list) {
  tryCatch(
    install_one(pkg),
    error = function(e) message(paste("ERROR installing", pkg, ":", e$message))
  )
}

message("R package check completed!")
