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

# R 包安装策略:
#   - conda (environment.yaml): DiffBind/ChIPseeker/pkgmaker/registry/rngtools/gridBase/dendextend/viridis/locfit/ps
#   - BiocManager::install: Bioconductor 包 (DESeq2/ORFik/DMRcaller/pheatmap/ComplexHeatmap/RNAmodR.RiboMethSeq/EBSeq)
#   - install.packages: CRAN 包 (dplyr/emmeans/car/AlgDesign/agricolae/multcomp/ggpubr/harmony)
#   - devtools::install_github: GitHub 私有包 (monocle3/NMF/riboWaltz/Seurat/DoubletFinder)
# pak 保留在 environment.yaml 里, 但不主动调用 (pak 会重装 conda 已装的依赖如 locfit/ps)
if (!requireNamespace("pak", quietly = TRUE)) {
  message("pak 未检测到，正在尝试自动安装...")
  .pak_ok <- tryCatch({
    install.packages("pak")
    requireNamespace("pak", quietly = TRUE)
  }, error = function(e) { FALSE })

  if (!isTRUE(.pak_ok)) {
    message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    message("❌ pak 自动安装失败，请手动安装后重跑：")
    message("")
    message("  # 在 R 里执行：")
    message("  install.packages(\"pak\")")
    message("")
    message("  # arm64 conda R 没有匹配的官方二进制，pak 只能源码编译安装，")
    message("  # 需要先确保有 C 编译器，二选一：")
    message("  #   xcode-select --install")
    message("  #   conda install -c conda-forge c-compiler")
    message("")
    message("  # 或用我们的预编译二进制 (官方推荐)：")
    message("  install.packages(\"pak\", repos = sprintf(")
    message("    \"https://r-lib.github.io/p/pak/stable/%s/%s/%s\",")
    message("    .Platform$pkgType, R.Version()$os, R.Version()$arch))")
    message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    quit(status = 1)
  }
  message("✓ pak 安装成功\n")
}

# ── Install dispatch ───────────────────────────────────────────────────────
# 安装函数: BiocManager (Bioconductor), install.packages (CRAN), devtools (GitHub)

install_bioconductor <- function(pkg) {
  if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
    message(paste("Installing", pkg, "(Bioconductor)..."))
    # BiocManager::install 不做重依赖求解, 只装目标包
    # 避免 pak 重装 conda 已装的 locfit/ps 等依赖
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    message(paste(pkg, "already installed"))
  }
}

install_github <- function(pkg, repo, ref = NULL, extra = NULL, remotes = FALSE) {
  if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
    message(paste("Installing", pkg, "(GitHub:", repo, ")..."))
    devtools::install_github(repo, ref = ref, dependencies = TRUE, remotes = remotes)
  } else {
    # 版本检查：如果指定了 GitHub ref，验证当前安装的 ref 是否匹配
    if (!is.null(ref)) {
      current_ref <- tryCatch(packageDescription(pkg)$GithubRef,
                               error = function(e) NULL)
      if (is.null(current_ref) || current_ref != ref) {
        message(paste("Updating", pkg, "to", ref, "branch..."))
        devtools::install_github(repo, ref = ref, dependencies = TRUE, remotes = remotes)
      } else {
        message(paste(pkg, "already installed (", ref, "branch)"))
      }
    } else {
      message(paste(pkg, "already installed"))
    }
  }
  # 安装辅助 companion 包 (如 Seurat 的 uwot — 也是 CRAN)
  if (!is.null(extra)) {
    for (ep in extra) {
      if (!suppressMessages(require(ep, character.only = TRUE, quietly = TRUE))) {
        message(paste("Installing companion package:", ep, "..."))
        install.packages(ep, repos = structure(c(CRAN = "https://cloud.r-project.org")))
      }
    }
  }
}

install_cran <- function(pkg) {
  if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
    message(paste("Installing", pkg, "(CRAN)..."))
    install.packages(pkg, repos = structure(c(CRAN = "https://cloud.r-project.org")))
  } else {
    message(paste(pkg, "already installed"))
  }
}


# ── monocle3 专用安装 ──
# monocle3 在 GitHub Remotes 里声明了两个额外的私有依赖:
#   - github::bnprks/BPCells/r
#   - github::cole-trapnell-lab/speedglm
# devtools 默认不拉 Remotes, 所以手动先装它们
install_monocle3 <- function() {
  # 1. BPCells (monocle3 Remotes 里的, 必须先装)
  if (!suppressMessages(requireNamespace("BPCells", quietly = TRUE))) {
    message("Installing BPCells (monocle3 Remotes)...")
    devtools::install_github("bnprks/BPCells/r", dependencies = TRUE)
  }
  # 2. speedglm fork (monocle3 Remotes 里的)
  if (!suppressMessages(requireNamespace("speedglm", quietly = TRUE))) {
    message("Installing speedglm (monocle3 Remotes)...")
    devtools::install_github("cole-trapnell-lab/speedglm", dependencies = TRUE)
  }
  # 3. monocle3 本身 (加 remotes=TRUE 再次确认)
  devtools::install_github("cole-trapnell-lab/monocle3")
}

# ── Single-package dispatcher ──────────────────────────────────────────────
install_one <- function(pkg) {
  switch(pkg,
    "dplyr"                   = install_cran("dplyr"),
    "DESeq2"                  = install_bioconductor("DESeq2"),
    "DMRcaller"               = install_bioconductor("DMRcaller"),
    "ORFik"                   = install_bioconductor("ORFik"),
    "monocle3"                = install_monocle3(),
    "pheatmap"                = install_bioconductor("pheatmap"),
    "ComplexHeatmap"          = install_bioconductor("ComplexHeatmap"),
    "RNAmodR.RiboMethSeq"     = install_bioconductor("RNAmodR.RiboMethSeq"),
    "emmeans"                 = install_cran("emmeans"),
    "car"                     = install_cran("car"),
    # AlgDesign 已移出 conda（conda-forge 无 osx-arm64 构建），
    # 是 agricolae 的硬依赖，须从 CRAN 源码装（需要 C 编译器）
    "AlgDesign"               = install_cran("AlgDesign"),
    "agricolae"               = install_cran("agricolae"),
    "multcomp"                = install_cran("multcomp"),
    "ggpubr"                  = install_cran("ggpubr"),
    "NMF"                     = install_github("NMF", "renozao/NMF", ref = "devel"),
    "riboWaltz"               = install_github("riboWaltz", "LabTranslationalArchitectomics/riboWaltz"),
    "Seurat"                  = install_github("Seurat", "satijalab/seurat", extra = c("uwot")),
    # EBSeq 已移出 conda（bioconda 最高只编到 r-base 4.3，与全环境 R 4.5 冲突；
    # rsem 用 conda install --no-deps 安装绕开了它）。
    # 仅当需要运行 rsem-run-ebseq 时才装：
    #   Rscript checkPackages.R --packages EBSeq
    "EBSeq"                   = install_bioconductor("EBSeq"),
    "DoubletFinder"           = install_github("DoubletFinder", "chris-mcginnis-ucsf/DoubletFinder"),
    "harmony"                 = install_cran("harmony"),
    {
      message(paste("Unknown package:", pkg, "- skipping"))
    }
  )
}


# ── Main ───────────────────────────────────────────────────────────────────
# conda 负责安装的 R 包 — 只做 require() 加载检查，不在这里安装
# （安装来源：environment.yaml 的 bioconductor-* / r-* conda 包）
CONDA_R_PACKAGES <- c(
  # Bioconductor (bioconda, osx-arm64 有构建)
  "DiffBind", "ChIPseeker", "ORFik",
  # Pure R (conda-forge) — NMF / ComplexHeatmap 依赖链
  "pkgmaker", "registry", "rngtools", "gridBase", "dendextend", "viridis"
)
# osx-arm64 注意事项:
#   - monocle3:    bioconda 只有 r36_0 老版本, 必须 GitHub 源码装 (cole-trapnell-lab/monocle3)
#                   Remotes: bnprks/BPCells/r + cole-trapnell-lab/speedglm
#   - ORFik:    bioconda 有 arm64 构建, 已加入 CONDA_R_PACKAGES
#   - DMRcaller/agricolae/AlgDesign: bioconda 无 arm64 构建, 由 BiocManager/install.packages 源码装

# 1. 先检查 conda 包是否可加载
conda_failed <- character(0)
for (pkg in CONDA_R_PACKAGES) {
  ok <- tryCatch(suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE)),
                 error = function(e) FALSE)
  if (isTRUE(ok)) {
    message(paste(pkg, "✓ (conda)"))
  } else {
    message(paste(pkg, "✗ — 请确保 environment.yaml 包含对应 conda 包后："))
    conda_failed <- c(conda_failed, pkg)
  }
}

if (length(conda_failed) > 0) {
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  message(paste("❌ conda 包加载失败：", paste(conda_failed, collapse = ", ")))
  message("请执行以下任一命令修复后重跑 checkPackages.R：")
  message("")
  message("  # 更新/重建 conda 环境：")
  message("  conda env update -f environment.yaml --prune")
  message("")
  message("  # 或手动补装缺失包：")
  # R named vector 里的横杠名必须用反引号/双引号, 否则会被解析为减法
  conda_pkgs <- c(
    "pkgmaker"      = "r-pkgmaker",
    "registry"      = "r-registry",
    "rngtools"      = "r-rngtools",
    "gridBase"     = "r-gridbase",
    "dendextend"    = "r-dendextend",
    "viridis"       = "r-viridis",
    "DiffBind"      = "bioconductor-diffbind",
    "ChIPseeker"    = "bioconductor-chipseeker",
    "ORFik"         = "bioconductor-orfik"
  )
  missing_args <- paste(conda_pkgs[conda_failed], collapse = " ")
  message(paste0("  conda install -c conda-forge -c bioconda ", missing_args))
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  quit(status = 1)
}

# 2. R 原生安装负责的包 — 自动安装 + 加载
if (is.null(pkg_list)) {
  # Install all known packages (DiffBind/ChIPseeker/pure-R-deps 由 conda 负责)
  # 注意：EBSeq 不在默认清单中（仅 rsem-run-ebseq 需要，见 install_one 注释），
  #       需要时单独执行 --packages EBSeq
  all_pkgs <- c("dplyr", "DESeq2", "ORFik", "monocle3", "DMRcaller", "pheatmap",
                 "ComplexHeatmap", "RNAmodR.RiboMethSeq",
                 "emmeans", "car", "AlgDesign", "agricolae", "multcomp", "ggpubr",
                 "NMF", "riboWaltz", "Seurat", "EBSeq", "DoubletFinder", "harmony")
  pkg_list <- all_pkgs
}

failed <- character(0)
for (pkg in pkg_list) {
  ok <- tryCatch(
    {
      install_one(pkg)
      suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))
    },
    error = function(e) {
      message(paste("ERROR installing", pkg, ":", e$message))
      # NMF 没有 CRAN 版本且必须跑 devel 分支，在 GHActions / 防火墙
      # 限制 / 无编译环境的机器上经常失败，这里给一条明确的手动指引。
      if (pkg == "NMF") {
        message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        message("⚠️  NMF 安装失败 — 请手动安装 devel 分支：")
        message("  1. 前往 https://github.com/renozao/NMF/tree/devel 下载源码")
        message("  2. 解压后执行：R CMD INSTALL NMF-devel")
        message("  3. 或直接在 R 里：pak::pak('renozao/NMF@devel')")
        message("                 / devtools::install_github('renozao/NMF', ref='devel')")
        message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
      }
      FALSE
    }
  )
  if (!isTRUE(ok)) {
    failed <- c(failed, pkg)
  }
}

if (length(failed) > 0) {
  message(paste("FAILED packages:", paste(failed, collapse = ", ")))
  quit(status = 1)
}

message("R package check completed!")
