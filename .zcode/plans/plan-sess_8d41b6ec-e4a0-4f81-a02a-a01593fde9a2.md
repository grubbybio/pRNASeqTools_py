## 目标

修复 prnaseqtools/modes/ribo.py 审查出的可确认 bug（含 cli.py 帮助文本失配）。原「RSEM 索引检查路径错误」一项经用户确认为非 bug（本机 RSEM 布局下 `RSEM_index/sjdbList.out.tab` 正确），已移除。

## 修复清单

### 1. 修复 resume ≥ step 4 的 NameError（高危）
将 `assembled_gtf_dir = "assembledGTF"` 与 `gffcomp_prefix = os.path.join(assembled_gtf_dir, genome)` 提升到主流程顶部共享目录变量区（现 ribo.py:166-173 一带），删除 step 3 内的重复定义（492、520 行）。

### 2. TE 样本分组配对（高危）
修改 step 10（974-984 行附近）：利用 `parse_input` 返回的第三元素 `pars`（组名+数量交替）解析出 control/treatment 组结构；按「control↔control、treatment i↔treatment i」配对并要求对应组 tag 数一致；结构不一致时 tee 显式 WARNING 并跳过无法确定的对，不再按全局下标盲目配对。

### 3. P-site 跨样本 read_name 冲突告警（高危-健壮性）
`_split_psites_by_sample` 中统计同时命中 ≥2 个样本 read set 的 read_name 数量，tee 报告数量/占比与示例；冲突比例高时醒目提示结果可能受影响。去重与归属核心逻辑不动（唯一依据仍是 read_name）。

### 4.（可选加固）RSEM 索引存在性检查双路径兼容
保留现有 `RSEM_index/sjdbList.out.tab` 检查为主，追加检测 `RSEM_index/RNA/sjdbList.out.tab` 作为备选（不同 RSEM 版本会建一级同名子目录）；任一存在即视为索引可用，避免版本差异导致无谓重建。默认行为不变。

### 5. samtools merge 单样本/空列表保护（中危）
新增小工具函数 `_merge_or_copy_bams(bam_list, out_prefix)`：
- len == 0 → sys.exit 明确报错；
- len == 1 → shutil.copy + samtools index；
- len ≥ 2 → 现有 merge + sort + index 流程。
替换 step 8 中 RNA 与 Ribo 两处调用。

### 6. step 9 清理与输出文件名对齐（中危）
step 9 开始清理时对每个 `all_ribo_tags` 精确删除旧的无扩展名 `f"P_sites_{tag}"`；不得删除作为 step 9 输入的 `P_sites_all*` 文件。

### 7. `_compute_te` 异常保护（中危）
step 10 循环内 `_compute_te(...)` 调用处包 try/except Exception，失败时 tee.write 原因并继续下一对，避免丢 STEP 10 COMPLETE 标记导致整步重跑。

### 8. `_assign_frames` 静默失败改为可见（中危）
移除裸 `except: pass`（1756-1757 行）：bedtools 不存在/返回码非 0/解析异常均 tee.write 明确告警（注明后续 frame 将全为 -1）。

### 9. cli.py 帮助文本同步
`--restart-step` help "(1-10)" → "(1-12)"（cli.py:147）。

### 10. 交互输入校验（低危）
step 8 交互 prompt 接受新值后校验各项均为正整数（正则），非法则提示并沿用原默认值。

### 不改行为、仅记录说明
- `create_metaplots.bash` 参数顺序 / start_stops_FAR.bed 拼写不一致：尝试在本机 ribotaper 目录核对安装脚本；能核实则统一字符串常量，核实不了保持现状（step 9 已双拼写兼容）。
- `--discard-untrimmed` 单端有/双端无的不对称；`--rf` + `--strandedness reverse` 写死：写入审查报告不改代码。

## 验证
1. `python -m py_compile prnaseqtools/modes/ribo.py prnaseqtools/cli.py`；
2. 控制流走查关键 resume 路径（last_step=3→step4、last_step=11→step10 等），确认所有跨块变量已定义；
3. 输出 git diff 汇总供复核。