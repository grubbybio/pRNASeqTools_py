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
# 镜像自动选择（国内 bioconductor.org / cloud.r-project.org 连接不稳定时用）
# 可以通过环境变量 R_MIRROR 强制指定: TUNA | USTC | ALIBABA | OFFICIAL
set_mirrors <- function() {
  # ── BiocManager 标准 repos 名称（源码里写死的）─────────────────────
  #   BioCsoft / BioCann / BioCexp / BioCworkflows / BioCbooks / CRAN
  #
  # TUNA 镜像只有 release 路径 (packages/release/bioc), 没有版本号路径
  # BiocManager 的 BioC_mirror option 会自动拼 packages/<version>/...
  #   → TUNA 没有版本号路径, 所以不能设 BioC_mirror!
  #   → 只手动设完整的 repos, 用 release 路径
  #   → 关掉 check_repositories 警告 (因为我们覆盖了标准 repos)

  .mirror_sets <- list(
    TUNA = list(
      CRAN_root = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
      BioC_root = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/release"
    ),
    USTC = list(
      CRAN_root = "https://mirrors.ustc.edu.cn/CRAN/",
      BioC_root = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/release"  # 中科大没同步 BioC
    ),
    ALIBABA = list(
      CRAN_root = "https://mirrors.aliyun.com/CRAN/",
      BioC_root = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/release"  # 阿里云没同步 BioC
    ),
    OFFICIAL = list(
      CRAN_root = "https://cloud.r-project.org",
      BioC_root = "https://bioconductor.org/packages/release"
    )
  )

  # ── 探测函数 ────────────────────────────────────────────────────────
  .ping <- function(url) {
    .ok <- FALSE
    tryCatch({
      .tf <- tempfile()
      download.file(url, .tf, method = "libcurl", quiet = TRUE, timeout = 5)
      .ok <- file.exists(.tf) && file.info(.tf)$size > 1000
      unlink(.tf)
    }, error = function(e) NULL)
    .ok
  }

  .apply_mirror <- function(.m, .label) {
    # 关键: 不设 BioC_mirror, 防止 BiocManager 用版本号拼 URL
    #       直接手动设完整的标准 repos
    options(repos = structure(c(
      CRAN        = .m$CRAN_root,
      BioCsoft    = paste0(.m$BioC_root, "/bioc"),
      BioCann     = paste0(.m$BioC_root, "/data/annotation"),
      BioCexp     = paste0(.m$BioC_root, "/data/experiment"),
      BioCworkflows = paste0(.m$BioC_root, "/workflows")
    )))
    # 关掉 "replaces Bioconductor standard repositories" 警告
    options(BiocManager.check_repositories = FALSE)
    message(sprintf("🌐 镜像 = %s → BioC=%s", .label, .m$BioC_root))
  }

  # ── 环境变量强制覆盖 ────────────────────────────────────────────────
  .force <- Sys.getenv("R_MIRROR", "")
  if (.force %in% names(.mirror_sets)) {
    .apply_mirror(.mirror_sets[[.force]], .force)
    return(invisible())
  }

  # ── 自动探测 ────────────────────────────────────────────────────────
  for (.name in c("TUNA", "USTC", "ALIBABA")) {
    .m <- .mirror_sets[[.name]]
    .ok <- all(c(
      .ping(paste0(.m$CRAN_root, "src/contrib/PACKAGES.gz")),
      .ping(paste0(.m$BioC_root, "/bioc/src/contrib/PACKAGES.gz")),
      .ping(paste0(.m$BioC_root, "/data/annotation/src/contrib/PACKAGES.gz"))
    ))
    if (.ok) {
      .apply_mirror(.m, .name)
      return(invisible())
    }
  }

  .apply_mirror(.mirror_sets$OFFICIAL, "OFFICIAL")
}
set_mirrors()

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# R 包安装策略:
#   - conda (environment.yaml): pkgmaker/registry/rngtools/gridBase/dendextend/viridis/locfit/ps (DiffBind: bioconda R 4.5 构建, 无 post-link)
#   - BiocManager::install: Bioconductor 包 (ORFik/DMRcaller/EBSeq/ChIPseeker/clusterProfiler)
#     rhdf5 已移 environment.yaml (Seurat/monocle3 依赖, conda 预编译 R 4.5 构建无 hdf5 系统库需求)
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
    message("  # 若源码编译失败（需要 C/C++ 编译器），确保系统有编译工具链：")
    message("  #   macOS:   xcode-select --install")
    message("  #           或 conda install -c conda-forge c-compiler cxx-compiler")
    message("  #   Linux:   conda install -c conda-forge c-compiler cxx-compiler")
    message("  #           或 apt install build-essential / yum groupinstall 'Development Tools'")
    message("")
    message("  # 或用预编译二进制 (官方推荐)：")
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
    # ⚠️ 不能用 BiocManager::install():
    #    1) 它不接受 repos 参数 (直接 stop)
    #    2) 内部 BiocManager::repositories() 会 c(bioc_versioned, our_repos)[!duplicated]
    #       把我们的 TUNA release 路径覆盖成它内部拼的 packages/<version>/bioc
    #       TUNA 没有版本号路径, 只有 packages/release/bioc
    #    直接用 install.packages(repos=getOption("repos")) 就够了
    #    install.packages 会按包类型自动从 BioCsoft / BioCann / CRAN 对应 repos 找
    suppressMessages(
      suppressWarnings(
        install.packages(pkg, repos = getOption("repos"),
                         update = FALSE, ask = FALSE, Ncpus = 1,
                         quiet = TRUE)
      )
    )
    # ── Mirror miss fallback ────────────────────────────────────────────
    # 国内镜像 (TUNA 等) 的 Bioconductor 同步可能滞后，缺核心包
    # 比如 Bioc 3.22 时期 TUNA 缺 BiocVersion / BiocBaseUtils
    # 导致 AnnotationHub 等装不上。当首次尝试失败时，
    # 自动用 BiocManager::repositories() (官方版本化 URL) 重试一次
    if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
      message(paste("  Mirror fallback: retrying", pkg, "from official Bioconductor..."))
      tryCatch({
        bm_repos <- BiocManager::repositories()
        suppressMessages(
          suppressWarnings(
            install.packages(pkg, repos = bm_repos,
                             update = FALSE, ask = FALSE, Ncpus = 1,
                             quiet = TRUE)
          )
        )
      }, error = function(e) {
        message(paste("  Fallback failed:", e$message))
      })
    }
  } else {
    message(paste(pkg, "already installed"))
  }
}

install_github <- function(pkg, repo, ref = NULL, extra = NULL, remotes = FALSE) {
  if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
    message(paste("Installing", pkg, "(GitHub:", repo, ")..."))
    suppressMessages(
      suppressWarnings(
        devtools::install_github(repo, ref = ref, dependencies = TRUE,
                                 remotes = remotes, quiet = TRUE)
      )
    )
  } else {
    # 版本检查：如果指定了 GitHub ref，验证当前安装的 ref 是否匹配
    if (!is.null(ref)) {
      current_ref <- tryCatch(packageDescription(pkg)$GithubRef,
                               error = function(e) NULL)
      if (is.null(current_ref) || current_ref != ref) {
        message(paste("Updating", pkg, "to", ref, "branch..."))
        suppressMessages(
          suppressWarnings(
            devtools::install_github(repo, ref = ref, dependencies = TRUE,
                                     remotes = remotes, quiet = TRUE)
          )
        )
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
        suppressMessages(suppressWarnings(install.packages(ep, quiet = TRUE)))
      }
    }
  }
}

install_cran <- function(pkg) {
  if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
    message(paste("Installing", pkg, "(CRAN)..."))
    suppressMessages(suppressWarnings(install.packages(pkg, quiet = TRUE)))
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
  # 1. BPCells — 已在 environment.yaml (r-bpcells / bioconductor-beachmat 预编译)
  # 2. 设置 C++17 标准 — assortHEAD (tatami 底层) 用了 std::exception_ptr (C++11),
  #    某些编译器默认 CXX_STD 没设, 导致 beachmat 等源码装时报:
  #      error: no member named 'exception_ptr' in namespace 'std'
  #    C++17 覆盖 C++11, 一劳永逸
  # CXX_STD = C++17 是 R 4.5+ 的标准入口, 同时覆盖 PKG_CXXFLAGS / PKG_CXX11FLAGS / PKG_CXX17FLAGS 双保险
  old_env_vars <- c("CXX_STD", "CXX17FLAGS", "PKG_CXXFLAGS", "PKG_CXX11FLAGS", "PKG_CXX17FLAGS")
  old_env <- sapply(old_env_vars, function(v) Sys.getenv(v, unset = NA))
  Sys.setenv(
    CXX_STD      = "C++17",
    CXX17FLAGS   = "-std=c++17",
    PKG_CXXFLAGS = "-std=c++17",
    PKG_CXX11FLAGS = "-std=c++17",
    PKG_CXX17FLAGS = "-std=c++17"
  )

  # 3. speedglm fork (monocle3 Remotes 里的) — 纯 R 包, 不受 C++ 影响
  if (!suppressMessages(requireNamespace("speedglm", quietly = TRUE))) {
    message("Installing speedglm (monocle3 Remotes)...")
    suppressMessages(suppressWarnings(
      devtools::install_github("cole-trapnell-lab/speedglm", dependencies = FALSE)
    ))
  }
  # 4. monocle3 本身
  install_bioconductor("batchelor")
  suppressMessages(suppressWarnings(
    devtools::install_github("cole-trapnell-lab/monocle3", dependencies = FALSE)
  ))

}


# ── Single-package dispatcher ──────────────────────────────────────────────
install_one <- function(pkg) {
  switch(pkg,
    "dplyr"                   = install_cran("dplyr"),
    "DESeq2"                  = install_bioconductor("DESeq2"),
    "ChIPseeker"              = install_bioconductor("ChIPseeker"),
    # Bioc 3.20+: GenomicFeatures 1.57 把 makeTxDb*() 迁出，独立包 txdbmaker
    "txdbmaker"               = install_bioconductor("txdbmaker"),
    "clusterProfiler"         = install_bioconductor("clusterProfiler"),
    "enrichplot"              = install_bioconductor("enrichplot"),
    "AnnotationHub"           = install_bioconductor("AnnotationHub"),
    "org.At.tair.db"          = install_bioconductor("org.At.tair.db"),
    "DMRcaller"               = install_bioconductor("DMRcaller"),
    "ORFik"                   = install_bioconductor("ORFik"),
    "monocle3"                = install_monocle3(),
    "pheatmap"                = install_bioconductor("pheatmap"),
    "ComplexHeatmap"          = install_bioconductor("ComplexHeatmap"),
    "RNAmodR.RiboMethSeq"     = install_bioconductor("RNAmodR.RiboMethSeq"),
    "emmeans"                 = install_cran("emmeans"),
    "car"                     = install_cran("car"),
    # AlgDesign 不在 conda 预装列表，是 agricolae 的硬依赖，
    # 须从 CRAN 源码装（需要 C 编译器）
    "AlgDesign"               = install_cran("AlgDesign"),
    "agricolae"               = install_cran("agricolae"),
    "multcomp"                = install_cran("multcomp"),
    "ggpubr"                  = install_cran("ggpubr"),
    "NMF"                     = install_github("NMF", "renozao/NMF", ref = "devel"),
    "riboWaltz"               = install_github("riboWaltz", "LabTranslationalArchitectomics/riboWaltz"),
    "Seurat"                  = message("Seurat 由 conda 安装 (r-seurat, R 4.5 预编译)"),
    # EBSeq 不在默认清单（仅 rsem-run-ebseq 需要），BiocManager 自动匹配当前 R 的 Bioconductor release 安装。
    # 若 rsem-run-ebseq 报 EBSeq missing，单独装：
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
  # Bioconductor (bioconda, R 4.5 构建无 post-link)
  "DiffBind",
  # Bioc 3.20+: GenomicFeatures 1.57 把 makeTxDb*() 迁出到 txdbmaker
  # bioconda 有预编译 bioconductor-txdbmaker，优先用 conda
  "txdbmaker",
  # Pure R (conda-forge) — NMF / ComplexHeatmap 依赖链
  "pkgmaker", "registry", "rngtools", "gridBase", "dendextend", "viridis"
)
# 跨平台注意事项:
#   - monocle3:    GitHub 源码装 (cole-trapnell-lab/monocle3)
#                   Remotes: bnprks/BPCells/r + cole-trapnell-lab/speedglm
#                   环境变量 CXX_STD=C++17 已在 install_monocle3() 里自动设置
#   - ORFik:       bioconda 没有构建, 由 BiocManager 安装
#   - DMRcaller/agricolae/AlgDesign: 部分平台无预编译, 由 BiocManager/install.packages 源码装

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
    "DiffBind"      = "bioconductor-diffbind",
    # Bioc 3.20+ makeTxDb*() 新宿主
    "txdbmaker"     = "bioconductor-txdbmaker",
    "Seurat"        = "r-seurat",
    "units"         = "r-units",
    "BPCells"       = "r-bpcells",
    "sf"                    = "r-sf",
    "beachmat"              = "bioconductor-beachmat",
    "assorthead"            = "bioconductor-assorthead",
    "SingleCellExperiment"  = "bioconductor-singlecellexperiment",
    "SummarizedExperiment"  = "bioconductor-summarizedexperiment",
    "ggrastr"              = "r-ggrastr",
    "spdep"                = "r-spdep",
    "pkgmaker"      = "r-pkgmaker",
    "registry"      = "r-registry",
    "rngtools"      = "r-rngtools",
    "gridBase"     = "r-gridbase",
    "dendextend"    = "r-dendextend",
    "viridis"       = "r-viridis",
  )
  missing_args <- paste(conda_pkgs[conda_failed], collapse = " ")
  message(paste0("  conda install -c conda-forge -c bioconda ", missing_args))
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  quit(status = 1)
}

# 2. R 原生安装负责的包 — 自动安装 + 加载
if (is.null(pkg_list)) {
  # Install all known packages (DiffBind 已回归 conda; ChIPseeker/clusterProfiler 走 BiocManager 镜像)
  # 注意：EBSeq 不在默认清单中（仅 rsem-run-ebseq 需要，见 install_one 注释），
  #       需要时单独执行 --packages EBSeq
  all_pkgs <- c("dplyr", "DESeq2", "ChIPseeker", "clusterProfiler",
                 "enrichplot", "AnnotationHub", "org.At.tair.db",
                 "ORFik", "monocle3", "DMRcaller", "pheatmap",
                 "ComplexHeatmap", "RNAmodR.RiboMethSeq",
                 "emmeans", "car", "AlgDesign", "agricolae", "multcomp", "ggpubr",
                 "NMF", "riboWaltz", "Seurat", "EBSeq", "DoubletFinder", "harmony",
                 # Bioc 3.20+ 新包: GenomicFeatures 1.57 把 makeTxDb*() 迁出
                 "txdbmaker")
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
