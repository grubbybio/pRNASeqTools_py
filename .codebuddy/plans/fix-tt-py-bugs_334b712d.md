---
name: fix-tt-py-bugs
overview: 修复 tt.py 中的多个漏洞和边界条件问题，包括列表越界、负索引序列提取错误、染色体缺失处理、命令注入风险等。
todos:
  - id: fix-negative-indices
    content: 修复负索引边界问题：第147行和256行 lim1 加上 max(0, ...)，第150行和258行增加序列长度上界保护
    status: completed
  - id: fix-locs-empty-crash
    content: 修复 _process_mir 中 locs 空列表崩溃：在第266行后添加 len(locs)
    status: completed
---

## 用户需求

理解 `prnaseqtools/modes/tt.py` 脚本逻辑，并修复其中存在的漏洞和潜在bug。

## 产品概览

该脚本是 pRNASeqTools 的 miRNA 截断/加尾（truncation/tailing）分析模块，流程为：bowtie 迭代比对（0-8错配）→ 合并 mapped reads → ShortStack 比对 → BAM 处理 → bedtools 与 miRNA/MIR GFF 取交集 → reads 分类统计 → R 语言 bubble plot 可视化。

## 核心修复目标

- **崩溃防护**：修复 `_process_mir` 中 `locs` 空列表导致的 IndexError
- **数据正确性**：修复负索引导致的错误序列提取
- **安全性**：缓解 `shell=True` 命令注入风险
- **健壮性**：添加坐标越界检查、文件名冲突修复、染色体缺失校验、无映射文件容错

## 技术栈

- Python 3，无外部依赖变更
- 复用现有 `run_cmd`、`revcomp`、`_tee` 等工具函数
- 使用标准库 `shlex.quote()` 缓解命令注入

## 实现方案

### 修复策略

本次修复遵循最小改动原则，在保持现有代码结构和 Perl 迁移风格的前提下逐点修复。不使用 list-based `subprocess` 参数重构（避免大规模变更），而是对 shell 命令中的用户可控变量使用 `shlex.quote()` 安全转义。

### 关键修复点（按严重程度排序）

**P0 — 崩溃类**：

1. `_process_mir` 函数（第262-269行）：`locs` 为空的保护。在 `sorted(set(locs))` 后添加 `if len(locs) < 2: continue`，同时在后续对 `locs[0]`、`locs[1]` 的访问前确保列表长度足够。

**P0 — 数据正确性类**：

2. 负索引边界（第147行、第256行）：`lim1 = mdata['start'] - flank` 改为 `lim1 = max(0, mdata['start'] - flank)`。当 miRNA/MIR 的 start 坐标小于 flank(25) 时，负索引会导致 Python 切片从序列末尾计数，提取到完全错误的片段。同样需要保护 `lim2` 不超过染色体长度。

**P1 — 安全性类**：

3. 命令注入（第110行及多处）：对 `tag`、`adaptor` 等用户可控变量使用 `shlex.quote()` 转义后拼入命令字符串。最危险的是 `rm -rf {tag}tmp`（第110行），若 tag 含 `; rm -rf / #` 等字符将导致任意命令执行。

**P1 — 数据完整性类**：

4. 坐标越界（第167/177/294/298行）：对 `fstart`、`fend` 添加 `min()`/`max()` 边界裁剪，确保不超出 `fas.get(chromosome, '')` 返回的序列长度。
5. 染色体缺失校验（第150/258行及多处）：在 `fas.get()` 返回空字符串时打印 warning 并跳过该记录，避免静默数据丢失。
6. 空映射容错（第93行）：在 `cat {tag}.mapped*.fastq` 前用 `glob` 检查是否存在匹配文件，无文件时跳过或创建空文件以避免 shell 报错。

**P2 — 健壮性类**：

7. 输出文件命名冲突（第190行 vs 第121行）：bedtools 第121行输出到 `{tag}.out`，第126行读取，第190行又覆盖写入 `{tag}.out`。将分类输出改为 `{tag}.categorized.out` 或类似名称，避免与 bedtools 原始输出混淆。

### 架构设计

无需架构变更。所有修改集中在 `tt.py` 单个文件的两个函数：`run()` 和 `_process_mir()`。每个修复点独立、可验证，不引入新的依赖或数据流变化。

### 目录结构

```
prnaseqtools/
└── modes/
    └── tt.py    # [MODIFY] 唯一修改文件，涉及约15处代码变更
```

### 修改清单

| 行号 | 问题 | 修复方式 |
| --- | --- | --- |
| 110 | `rm -rf {tag}tmp` 命令注入 | `shlex.quote(tag)` |
| 147 | `lim1 = mdata['start'] - flank` 负索引 | `max(0, mdata['start'] - flank)` |
| 150 | `fas.get(...)[lim1:lim2+1]` 负索引 | 同上 + 序列长度边界 |
| 167 | `fas.get(...)[fstart:fend]` 越界 | 添加 `min(fend, seq_len)` 边界 |
| 177 | 同上，负链 | 同上 |
| 256 | `lim1 = mdata['start'] - flank` 负索引 | `max(0, mdata['start'] - flank)` |
| 258 | `fas.get(...)[lim1:lim2+1]` 负索引 | 同上 + 序列长度边界 |
| 262-269 | `locs` 空列表 IndexError | `if len(locs) < 2: continue` |
| 268 | `locs[0]`/`locs[1]` 访问 | 增加 len 检查 |
| 269 | `locs[2]`/`locs[3]` 访问 | 已有 `len(locs) > 2` 保护，补充注释 |
| 93 | `cat *.fastq` 无文件 | 先 glob 检查再执行 |
| 121/190 | 输出覆盖输入 | 第190行改为 `{tag}.categorized.out` |
| 多处 | 染色体缺失 | `fas.get()` 后检查非空 |
| 多处 | `tag`/`adaptor` 命令注入 | `shlex.quote()` 包裹 |


### 实现注意事项

- **性能**：`shlex.quote()` 和 `max(0, ...)` 为 O(1) 操作，无性能影响
- **日志**：染色体缺失和坐标越界使用 `tee.write()` 输出 warning，不中断流程
- **向后兼容**：输出文件名变更为 `{tag}.categorized.out`，需确认下游 R 脚本 `bubble_plot.R` 的输入是否依赖固定文件名
- **blast radius**：仅修改 `tt.py`，不影响其他模式