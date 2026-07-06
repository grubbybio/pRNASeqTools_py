---
name: fix-wgbs-py-bugs
overview: 修复 wgbs.py 中的多个漏洞：os.symlink 重复运行崩溃、rename 无检查、bedgraph 生成在所有分支缺失、_bin_methylation 初始化 bug、死代码和清理模式优化等。
todos:
  - id: fix-symlink-and-fasta-check
    content: 修复 symlink 重复运行崩溃（line 50、159）和 fasta_path 存在性检查（line 50），同时添加 import shutil
    status: completed
  - id: fix-rename-checks
    content: 修复 cutadapt 和 deduplicate 后 rename 无存在性检查（line 68、74、89-90、96、112-113、119），添加 warning 日志
    status: completed
  - id: promote-bedgraph-bgzip
    content: 将 bedgraph/bgzip/tabix 代码块从 comma 分支提升到 for 循环尾部共用（line 127-137），确保所有分支生成 CX_report.txt.gz
    status: completed
  - id: fix-bin-methylation-bugs
    content: 修复 _bin_methylation 中 chr_name 初始化逻辑（line 207-208）和删除 sum_vals 死代码（line 183、233）
    status: completed
  - id: fix-glob-and-rmrf
    content: 优化 glob 清理模式（line 142）和替换 rm -rf 为 shutil.rmtree（line 152）
    status: completed
---

## 用户需求

理解 prnaseqtools/modes/wgbs.py 脚本逻辑，并修复其中存在的 bug。

## 产品概述

该脚本是 pRNASeqTools 的 WGBS 全基因组亚硫酸氢盐测序分析模块，流程为：Bismark 基因组索引 → 样本比对（单端/双端SRA或本地文件） → deduplicate → methylation extraction → bin 甲基化统计 → DMRcaller 差异甲基化分析。支持 mapping_only 和 no_mapping 两种跳过模式。

## 核心修复目标

- **崩溃修复**：symlink 重复运行 FileExistsError、rename 无存在性检查
- **功能补全**：bedgraph/bgzip/tabix 提升到所有分支共用，确保 no_mapping 和 DMRcaller 正常工作
- **逻辑修正**：_bin_methylation 首行染色体名初始化 bug、死代码清理
- **健壮性增强**：fasta_path 存在性检查、glob 清理模式精确化

## 技术栈

- Python 3，复用现有 `run_cmd`、`_tee`、`unzip_file`、`download_sra` 等工具函数
- 新增 `import shutil` 用于替代 `rm -rf` shell 命令
- 不引入外部依赖

## 实现方案

### 修复策略

遵循最小改动原则，在保持现有代码结构和 Bismark 流程的基础上逐点修复。bedgraph/bgzip/tabix 代码块从 comma 分支内部提升到循环体尾部统一执行。

### 修复清单

| # | 严重程度 | 行号 | 问题 | 修复方式 |
| --- | --- | --- | --- | --- |
| 1 | P0 | 50 | `os.symlink()` 重复运行抛 `FileExistsError` | 先 `os.unlink` 再 symlink |
| 2 | P0 | 50 | `fasta_path` 无存在性检查 | 添加 `os.path.exists` 检查 |
| 3 | P0 | 68/89-90/112-113 | cutadapt 后 rename 无检查 | 添加 `os.path.exists` + warning |
| 4 | P0 | 74/96/119 | deduplicate 后 rename 无检查 | 添加 `os.path.exists` + warning |
| 5 | P0 | 127-137 | bedgraph/bgzip/tabix 仅在 comma 分支 | 提升到 for 循环尾部，三分支共用 |
| 6 | P1 | 207-208 | `if line_num == 0` 首行被跳过时 chr_name 不初始化 | 改为 `if chr_name is None` |
| 7 | P1 | 183/233 | `sum_vals` 递增后从未读取 | 删除定义和累加语句 |
| 8 | P1 | 142 | `C*_O?_*{tag}*.txt` 模式过宽 | 优化为 `*_OT_*{tag}*.txt *_OB_*{tag}*.txt` |
| 9 | P2 | 152 | `rm -rf Bisulfite_Genome` shell=True | 替换为 `shutil.rmtree` |


### 关键设计决策

**bedgraph/bgzip/tabix 提升方案**：将第127-137行整体移到 `else` 块外（与第139行 `tee.write` 同级），使单端SRA、双端SRA、本地双端三个分支都能生成 `CX_report.txt.gz`。DMRcaller.R 第23行读取 `{genotype}_{k}.CX_report.txt.gz`，no_mapping 模式第159行也依赖该文件，因此所有分支必须生成。

**`_bin_methylation` chr_name 初始化修复**：原代码 `if line_num == 0: chr_name = cols[0]` 中 `line_num` 包含被跳过的行（`len(cols) < 6`），若文件首行恰好被跳过则 `chr_name` 保持 `None`，后续 `if curr_chr != chr_name` 对所有行都成立（因为 `None != 任意值`），导致每个新染色体都误判为切换、bin 被过早写出。修复为 `if chr_name is None: chr_name = cols[0]`，只初始化一次。

**no_mapping 模式 symlink**：第159行 `os.symlink(f"../{pre}.CX_report.txt.gz", ...)` 同样存在重复运行崩溃风险，一并修复。

### 目录结构

```
prnaseqtools/
└── modes/
    └── wgbs.py    # [MODIFY] 唯一修改文件，涉及约20处代码变更
```