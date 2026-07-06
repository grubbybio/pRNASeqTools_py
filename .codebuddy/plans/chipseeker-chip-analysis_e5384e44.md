---
name: chipseeker-chip-analysis
overview: 修复 tf.py 中 chip 模式的 bug，参考桌面 chipseeker.R 添加 R 脚本对 5 类 peak 进行 ChIPseeker 注释和 clusterProfiler GO 富集分析
todos:
  - id: fix-tf-py-bugs
    content: 修复 tf.py 中 chip 模块的 bugs：简化 _parse_control_fragments、去除冗余 qvalue 重定义、优化 cutoff 格式处理
    status: completed
  - id: create-chipseeker-r
    content: 创建 scripts/tf_chipseeker.R 脚本，实现 ChIPseeker 注释 + clusterProfiler GO 富集分析，支持 5 类 peak 文件输入
    status: completed
  - id: integrate-r-call
    content: 在 tf.py chip 模块 QC 之后、清理之前，集成调用 tf_chipseeker.R 的代码，包含 GFF 文件存在性检查
    status: completed
    dependencies:
      - fix-tf-py-bugs
      - create-chipseeker-r
---

## 需求概述

修改 `prnaseqtools/modes/tf.py` 中的 chip 分析模块，在分析流程最后对 5 类 peak 结果使用 ChIPseeker 和 clusterProfiler 进行注释和功能富集分析。同时修复脚本中存在的 bug。

## 核心功能

1. **Bug 修复**:

- 简化 `_parse_control_fragments` 函数：外层 `if` 已确认行内容，内层循环冗余，直接取 `parts[-1]` 即可
- 去除 chip 块内冗余的 `qvalue` 变量重定义（L191）
- 确保 `diff_cutoff` 值格式与 MACS3 bdgdiff 输出文件名一致
- 在调用 R 脚本前检查 GFF 参考文件是否存在

2. **ChIPseeker + clusterProfiler 分析**:

- 新建 `scripts/chipseeker.R` 脚本，接收基因组名、prefix、组名、cutoff 等参数
- 读取 chip 模式产出的 5 类 peak 文件（2 个 callpeak narrowPeak + 3 个 bdgdiff BED）
- 使用 `makeTxDbFromGFF` 从参考 GFF 构建 TxDb 对象
- 使用 `annotatePeak` 对每类 peak 进行基因组注释（TSS 区域 -3000~0，侧翼距离 3000）
- 提取注释结果的 geneId，运行 `enrichGO` 进行 GO 富集分析
- 保存注释表格、注释统计图、GO 富集点图

3. **集成到 tf.py**:

- 在 QC 分析之后、清理之前，调用 `Rscript --vanilla` 执行 `chipseeker.R`
- 沿用现有 R 脚本调用模式（如 `tf_mrna.R`、`tf_srna.R`）

## 技术方案

### 技术栈

- **Python 3**: 主流程编排（已有）
- **R**: ChIPseeker, clusterProfiler, GenomicFeatures, enrichplot, ggplot2（新增 R 脚本）
- **参考文件**: `{prefix}/reference/{genome}_genes.gff`（已有）

### 实现策略

**Bug 修复方案**:

1. **`_parse_control_fragments` 简化**: 将内层循环+回退的复杂结构替换为直接取 `parts[-1]` 后 `int()` 转换的简洁逻辑。MACS3 输出格式稳定，外层 `'fragments after filtering in control' in line` 已精确定位到目标行。

2. **冗余变量清理**: 删除 chip 块内 L191 的 `qvalue` 重定义，改用外层统一变量。外层已有 L29 的 `qvalue` 定义，内层重复获取相同值无意义。

3. **diff_cutoff 格式一致性**: 保留 `{diff_cutoff:.1f}` 格式化方式（argparse `type=float` 保证值为 float），与 MACS3 bdgdiff 输出文件名 `_c{value:.1f}_` 格式一致。

4. **GFF 存在性检查**: 在调用 R 脚本前，通过 `os.path.exists()` 确认 `{prefix}/reference/{genome}_genes.gff` 存在。

**R 脚本设计**:

- 参数接收: `genome`、`prefix`、`group1_name`、`group2_name`、`cutoff`
- 5 个 peak 文件路径由参数在 R 脚本内拼接
- 基因组注释数据库映射: `ath` -> `org.At.tair.db`（keyType="TAIR"），其他基因组若缺少对应 OrgDb 则跳过 GO 富集

**集成方案**:

在 tf.py chip 块中 QC 代码之后、清理代码之前插入 R 脚本调用：

```python
tee.write(f"\n{'='*60}\n")
tee.write("ChIPseeker Peak Annotation & GO Enrichment\n")
tee.write(f"{'='*60}\n")
gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
if not os.path.exists(gff_path):
    tee.write(f"  GFF not found: {gff_path}, skipping annotation.\n")
else:
    run_cmd(
        f"Rscript --vanilla {prefix}/scripts/chipseeker.R "
        f"{genome} {prefix} {group1_name} {group2_name} {diff_cutoff}"
    )
```

### 架构设计

```
┌─────────────────────────────────────────────────────────────────┐
│                      tf.py (chip mode)                          │
│                                                                  │
│  MACS3 callpeak (G1) → MACS3 callpeak (G2)                      │
│       ↓                         ↓                                │
│  Parse fragment counts from output                               │
│       ↓                                                          │
│  MACS3 bdgdiff (differential peak calling)                       │
│       ↓                                                          │
│  Peak QC (5 peak sets)                                           │
│       ↓                                                          │
│  ── [新增] ChIPseeker + clusterProfiler 注释分析 ──              │
│       ↓                                                          │
│  Cleanup temp files                                              │
└─────────────────────────────────────────────────────────────────┘
         ↕ Rscript --vanilla
┌─────────────────────────────────────────────────────────────────┐
│                scripts/chipseeker.R                            │
│                                                                  │
│  1. makeTxDbFromGFF(gff_path) → TxDb                            │
│  2. readPeakFile() for 5 peak files                             │
│  3. annotatePeak() on each peak set                             │
│  4. extract geneId → enrichGO()                                 │
│  5. Save annotation tables + plots                               │
└─────────────────────────────────────────────────────────────────┘
```

### 实现细节

#### 目录结构

```
project-root/
├── prnaseqtools/
│   └── modes/
│       └── tf.py                    # [MODIFY] 在 chip 块 QC 后添加 R 脚本调用，修复 bugs
└── scripts/
    └── chipseeker.R              # [NEW] ChIPseeker + clusterProfiler 分析脚本，供 tf / chip 模块共用
```

#### 关键数据结构

R 脚本参数设计：

```
# Usage: Rscript chipseeker.R <genome> <prefix> <group1_name> <group2_name> <cutoff>
#
# Peak 文件自动拼接:
#   1. {g1}_peaks.narrowPeak
#   2. {g2}_peaks.narrowPeak
#   3. diff_{g1}_vs_{g2}_c{cutoff}_cond1.bed
#   4. diff_{g1}_vs_{g2}_c{cutoff}_cond2.bed
#   5. diff_{g1}_vs_{g2}_c{cutoff}_common.bed
#
# 输出文件:
#   - {name}_annotation.txt          # 注释结果表
#   - {name}_go_enrichment.txt       # GO 富集结果
#   - annotation_barplot.pdf         # 注释类型分布柱状图
#   - go_dotplot.pdf                 # GO 富集点图
```

OrgDb 映射逻辑：

- `ath` -> `org.At.tair.db` (keyType = "TAIR")
- 其他基因组: 尝试加载 `org.{genome}.db`，失败则跳过 GO 富集，仅输出注释结果

### 性能考量

- R 脚本在每个 peak 集上独立运行 `annotatePeak`，5 个 peak 集共享同一 TxDb 对象
- `makeTxDbFromGFF` 在 GFF 较大时较慢，但只在脚本启动时构建一次
- GO 富集对每个 peak 集的基因列表独立运行，基因列表非空时执行