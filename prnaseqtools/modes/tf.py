"""
Two-factor DE analysis mode.
Runs inside srna/mrna output folders.
"""

import os
import sys
import shutil
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools.functions import _tee, run_cmd
from prnaseqtools import reference as ref


def run(opts):
    """Main entry point for two-factor DE analysis."""
    opts = validate_options(opts)
    tee = _tee()

    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    foldchange = opts.get('foldchange', 1.5)
    pvalue = opts.get('pvalue', 0.05)
    control = opts.get('control', '')
    treatment = opts.get('treatment', '')
    if isinstance(treatment, list):
        treatment = treatment[0]
    run_mode = opts.get('run_mode', 'mrna')
    norm = opts.get('norm', 'rRNA,total')
    norms = norm.split(',')
    binsize = opts.get('binsize', 100)
    genome = opts.get('genome', 'ath')
    deseq2_norm = opts.get('deseq2_norm', 'DESeq2')
    qvalue = opts.get('qvalue', 1.0)
    chip_method = opts.get('chip_method', 'diffbind')
    chip_analysis = opts.get('chip_analysis', 'dual_factor')
    chip_norm = opts.get('chip_norm', 'deseq2')
    qvalue_label = f"Q value = {qvalue}" if qvalue < 1 else f"P value = {pvalue}"

    tee.write(f"Two-factor comparison between {control} and {treatment}\n"
              f"Fold change = {foldchange} {qvalue_label}\n")

    # ── 通用解析函数: 解析 "groupName=label1,N1,label2,N2" 格式 ──────────
    # 对所有模式 (mrna/srna/chip) 使用同一格式。
    # 返回 (group_name, [(label1, rep_count1), (label2, rep_count2), ...],
    #         flat_list_for_R_scripts)
    def _parse_spec(spec):
        parts = spec.split('=')
        if len(parts) != 2:
            sys.exit(
                "Format: groupName=label1,N1,label2,N2  "
                f"(got: {spec})")
        group_name = parts[0]
        fields = parts[1].split(',')
        if len(fields) < 2 or len(fields) % 2 != 0:
            sys.exit(
                "Format: groupName=label1,N1,label2,N2  "
                f"(got: {spec})")
        pairs = []
        flat = [group_name]
        for i in range(0, len(fields), 2):
            label = fields[i]
            try:
                n = int(fields[i + 1])
            except ValueError:
                sys.exit(
                    f"Replicate count must be integer (got: {fields[i + 1]})")
            if n < 1:
                sys.exit("Replicate count must be >= 1")
            pairs.append((label, n))
            flat.append(label)
            flat.append(str(n))
        return group_name, pairs, flat

    # ── 解析两组实验条件 ────────────────────────────────────────────────
    g1_name, g1_pairs, g1_flat = _parse_spec(control)
    g2_name, g2_pairs, g2_flat = _parse_spec(treatment)

    if len(g1_pairs) != len(g2_pairs):
        sys.exit("Please provide paired data! Both groups must have "
                 "the same number of condition labels.")

    par_str = ' '.join(g1_flat + g2_flat)

    # ── 生成所有样本标签 ─────────────────────────────────────────────────
    tags = []
    for name, pairs in [(g1_name, g1_pairs), (g2_name, g2_pairs)]:
        for label, n in pairs:
            for rep in range(1, n + 1):
                tags.append(f"{name}_{label}_{rep}")

    if run_mode == 'mrna':
        for pre in tags:
            os.symlink(f"../{pre}.txt", f"{pre}.txt")

        run_cmd(
            f"Rscript --vanilla {prefix}/scripts/tf_mrna.R "
            f"{deseq2_norm} {pvalue} {foldchange} {par_str}")

        for fname in globmod.glob("*_?.txt"):
            os.unlink(fname)

    elif run_mode == 'srna':
        for pre in tags:
            os.symlink(f"../{pre}.nf", f"{pre}.nf")
            for fname in os.listdir(".."):
                if fname.startswith(pre) and fname.endswith('count') and 'norm' not in fname:
                    os.symlink(f"../{fname}", fname)

        for mnorm in norms:
            run_cmd(
                f"Rscript --vanilla {prefix}/scripts/tf_srna.R "
                f"{mnorm} {pvalue} {foldchange} {par_str}")

            # Generate bedgraph from results
            csv_files = [f for f in os.listdir('.') if f.endswith('.csv') and mnorm in f and 'bin' in f]
            for hcsv in [f for f in csv_files if 'hyper' in f or 'hypo' in f]:
                bg_file = hcsv.replace('.csv', '.bedgraph')
                with open(hcsv) as fh_in, open(bg_file, 'w') as fh_out:
                    for line in fh_in:
                        line = line.strip()
                        if not line:
                            continue
                        cols = line.split(',')
                        if not cols:
                            continue
                        m = __import__('re').match(r'(\w+)_(\d+)', cols[0].strip('"'))
                        if m:
                            chr_name = m.group(1)
                            start = int(m.group(2)) * binsize
                            end = start + binsize - 1
                            fh_out.write(f"{chr_name}\t{start}\t{end}\t{cols[2] if len(cols) > 2 else '0'}\n")

            # Annotate
            ann = ref.build_annotation(prefix, genome, binsize)
            for csv_file in csv_files:
                tmp_file = "tmp4"
                with open(csv_file) as fh_in, open(tmp_file, 'w') as fh_out:
                    for line in fh_in:
                        line = line.strip()
                        if not line:
                            continue
                        cols = line.split(',')
                        if not cols:
                            continue
                        key = cols[0].strip('"')
                        if key in ann:
                            fh_out.write(f"{line},{ann[key]}\n")
                        else:
                            fh_out.write(f"{line},Intergenic\n")
                os.rename(tmp_file, csv_file)

            run_cmd(
                f"Rscript --vanilla {prefix}/scripts/tf_mirna.R "
                f"{mnorm} {pvalue} {foldchange} {par_str}")

        for fname in globmod.glob("*.count"):
            os.unlink(fname)
        for fname in globmod.glob("*.nf"):
            os.unlink(fname)

    elif run_mode == 'chip':
        # ── Chip 模式辅助: 从 pairs 提取 input/IP 标签列表 ──────────────
        def _tags_from_pairs(group_name, pairs):
            """按约定 pairs[0]=input, pairs[1]=IP, 生成标签列表。"""
            if len(pairs) < 2:
                sys.exit(f"{group_name}: need at least 2 condition labels "
                         "(input_label,N,IP_label,M)")
            in_label, in_n = pairs[0]
            ip_label, ip_n = pairs[1]
            input_tags = [f"{group_name}_{in_label}_{r}"
                          for r in range(1, in_n + 1)]
            ip_tags = [f"{group_name}_{ip_label}_{r}"
                       for r in range(1, ip_n + 1)]
            return input_tags, ip_tags
        # ── ChIP-seq differential peak calling ───────────────────────────
        # 参数格式与 mrna/srna 统一:
        #   --control "WT=input,3,IP,3"
        #   --treatment "KO=input,2,IP,2"
        # BAM 文件: {group}_{label}_{rep}.sorted.bam
        tee.write(f"ChIP-seq differential peak calling (method={chip_method})\n")

        # Helper: 统计 BAM 里的 total reads 或 mito reads
        def _count_bam_reads(bam_file, mito_mode=False, organelle_mode=False):
            """统计 BAM 中的总 reads / 线粒体 reads / 叶绿体 reads。

            统一用 samtools idxstats（不需要 .bai index 文件），一次调用同时拿到:
              - 所有染色体名（= header 里有但无 reads 的不会出现）
              - 每个染色体的 reads 数
            格式: chr\\tlength\\tmapped_reads\\tunmapped_reads

            参数:
              mito_mode (bool, deprecated): 兼容旧调用，传 True 时等价于
                                            organelle_mode='mito'
              organelle_mode (False|'mito'|'chloro'): 要统计的细胞器 reads
                                                    False 表示总 reads
            """
            import subprocess as _sp
            # 兼容旧接口
            if mito_mode and not organelle_mode:
                organelle_mode = 'mito'
            res = _sp.run(['samtools', 'idxstats', bam_file],
                          capture_output=True, text=True)
            if res.returncode != 0:
                # idxstats 失败（BAM 损坏）→ 回退 samtools view -c
                res2 = _sp.run(['samtools', 'view', '-c', bam_file],
                               capture_output=True, text=True)
                try:
                    return int(res2.stdout.strip())
                except ValueError:
                    return 0
            # 解析 idxstats: chr\tlength\tmapped\tunmapped
            chr_counts = {}
            for line in res.stdout.splitlines():
                parts = line.split('\t')
                if len(parts) >= 3:
                    chr_counts[parts[0]] = int(parts[2])

            if not chr_counts:
                return 0
            total = sum(chr_counts.values())

            if not organelle_mode:
                return total

            # 细胞器识别：从 idxstats 的 chr 名里找（不需要额外 samtools 调用！）
            # 植物常用命名:
            #   线粒体: chrM, MT, chrMT, mitochondria
            #   叶绿体: chrC, chrCP, chloroplast, Pt, cp
            organelle_patterns = {
                'mito':   ('chrm', 'mt', 'chrmt', 'mitochondria',
                           'mito', 'mitochondr'),
                'chloro': ('chrc', 'chrcp', 'chloroplast', 'pt', 'cp',
                           'plastid', 'pltd'),
            }
            organelle_fallback = {
                'mito':   ('chrM', 'MT', 'chrMT', 'mitochondria', 'ChrM', 'mito'),
                'chloro': ('chrC', 'chrCP', 'chloroplast', 'ChrC', 'Pt', 'cp'),
            }
            # rDNA normalization: 拟南芥 chr2:0-10500 + chr3:14193500-14204500
            # 精确区间统计 (不用 idxstats, 用 samtools view -c 支持区间)
            rdna_regions = ['chr2:0-10500', 'chr3:14193500-14204500']
            if organelle_mode == 'rdna':
                import subprocess as _sp2
                rdna_total = 0
                for region in rdna_regions:
                    res = _sp2.run(
                        ['samtools', 'view', '-c', bam_file, region],
                        capture_output=True, text=True)
                    try:
                        rdna_total += int(res.stdout.strip())
                    except ValueError:
                        pass
                return rdna_total
            patterns = organelle_patterns.get(organelle_mode, ())
            org_chrs = []
            for c in chr_counts:
                cl = c.lower()
                if cl in patterns:
                    org_chrs.append(c)
                elif organelle_mode == 'mito' and (
                        'mito' in cl or cl.endswith('_mt') or cl.endswith('-mt')):
                    org_chrs.append(c)
                elif organelle_mode == 'chloro' and (
                        'chloro' in cl or 'plastid' in cl):
                    org_chrs.append(c)
            if not org_chrs:
                # 回退候选列表 (mito/chloro)
                for cand in organelle_fallback.get(organelle_mode, ()):
                    if cand in chr_counts:
                        org_chrs.append(cand)
            return sum(chr_counts.get(c, 0) for c in org_chrs)

        # ── 公共：参数 + 解析 --control / --treatment → tags → BAM ────────
        seq_strategy = opts.get('seq_strategy', 'paired')
        genome_size = opts.get('genome_size')
        if not genome_size:
            sys.exit("--genome-size is required for ChIP analysis "
                     "(e.g. 1.35e8 for ath)")
        fmt = "BAMPE" if seq_strategy == 'paired' else "BAM"
        tss_distance = opts.get('tss_distance', 3000)

        # 用统一的 _parse_spec 解析，按约定: pair[0]=Input标签, pair[1]=IP标签
        control_opt = opts.get('control', '')
        if isinstance(control_opt, list):
            control_opt = control_opt[0] if control_opt else ''
        group1_name, g1_pairs, _ = _parse_spec(control_opt)

        treatment_opt = opts.get('treatment', '')
        if isinstance(treatment_opt, list):
            treatment_opt = treatment_opt[0] if treatment_opt else ''
        group2_name, g2_pairs, _ = _parse_spec(treatment_opt)

        # 从 pairs 中提取 input 和 IP 标签列表
        g1_input_tags, g1_ip_tags = _tags_from_pairs(group1_name, g1_pairs)
        g2_input_tags, g2_ip_tags = _tags_from_pairs(group2_name, g2_pairs)

        tee.write(f"  Group 1 ({group1_name}):\n")
        for t in g1_input_tags:
            tee.write(f"    Input: {t}.sorted.bam\n")
        for t in g1_ip_tags:
            tee.write(f"    IP:    {t}.sorted.bam\n")
        tee.write(f"  Group 2 ({group2_name}):\n")
        for t in g2_input_tags:
            tee.write(f"    Input: {t}.sorted.bam\n")
        for t in g2_ip_tags:
            tee.write(f"    IP:    {t}.sorted.bam\n")

        # ── 检查 BAM 文件是否存在 ─────────────────────────────────────────
        all_tags = []
        for t in g1_input_tags:
            all_tags.append((t, f"{group1_name} Input"))
        for t in g1_ip_tags:
            all_tags.append((t, f"{group1_name} IP"))
        for t in g2_input_tags:
            all_tags.append((t, f"{group2_name} Input"))
        for t in g2_ip_tags:
            all_tags.append((t, f"{group2_name} IP"))

        # ── 查找 BAM 文件（支持 .sorted.bam 或 .sorted.dedup.bam）───────────
        def _resolve_bam(tag, desc):
            """返回 dedup 后的 BAM 文件路径。

            优先级:
              1. {tag}.sorted.dedup.bam   — 已 dedup 的首选（cwd 或 ../）
              2. {tag}.sorted.bam        — 用 sorted.bam 生成 .dedup.bam
              3. 报错退出

            Picard 3.x 要求 read group (RG) 存在,
            因此生成 dedup 前先调 AddOrReplaceReadGroups。

            注意:
              1. Picard 3.5.0 在 macOS 挂载卷 (/Volumes) 上
                 会因 writability 检查失败而崩溃。绕开方法:
                 Picard 输出全部写到 /tmp, 然后 shutil.copy
                 到 sorted.bam 源文件所在目录 (即 ../)。
              2. 新生成的 dedup BAM 放在原 BAM 同级目录
                 (如 ../{tag}.sorted.dedup.bam), 不在 cwd 里。
                 cwd 里只放 symlink 供下游工具访问。
              3. 函数返回的是 **真实路径** (用于 samtools idxstats
                 等需要真实文件位置的命令)。
            """
            # ── 1. 优先用现成的 dedup (cwd 或 ../) ──
            for prefix in ('.', '..'):
                bam = f"{prefix}/{tag}.sorted.dedup.bam"
                if os.path.exists(bam):
                    real = os.path.realpath(bam)
                    # cwd 只有 symlink → 把真实路径 return 出去
                    if prefix == '.':
                        return real
                    return real

            # ── 2. 找 sorted.bam ──
            sorted_bam = None
            sorted_bam_dir = '.'
            for prefix in ('.', '..'):
                cand = f"{prefix}/{tag}.sorted.bam"
                if os.path.exists(cand):
                    sorted_bam = cand
                    sorted_bam_dir = prefix
                    break
            if sorted_bam is None:
                sys.exit(f"BAM file not found: {tag}.sorted.bam ({desc})")

            # Picard 严格拒绝 '..' 相对路径 → 先 symlink 到当前目录
            if sorted_bam.startswith('..'):
                local_sorted = f"{tag}.sorted.bam"
                if not os.path.exists(local_sorted):
                    os.symlink(sorted_bam, local_sorted)
                sorted_bam = local_sorted

            # ── 3. Picard: /tmp 临时目录跑 ──
            picard_tmp = "/tmp/picard_dedup"
            os.makedirs(picard_tmp, exist_ok=True)
            rg_bam_tmp = f"{picard_tmp}/{tag}.sorted.rg.bam"
            dedup_bam_tmp = f"{picard_tmp}/{tag}.sorted.dedup.bam"
            metrics_tmp = f"{picard_tmp}/{tag}.sorted.dup_metrics.txt"

            # 3a. Add read group (Picard 3.x 必需要)
            if not os.path.exists(rg_bam_tmp):
                tee.write(f"    Picard AddOrReplaceReadGroups: "
                          f"{os.path.basename(sorted_bam)} → tmp rg\n")
                run_cmd(
                    f"picard AddOrReplaceReadGroups "
                    f"I={sorted_bam} "
                    f"O={rg_bam_tmp} "
                    f"RGID={tag} "
                    f"RGLB={tag} "
                    f"RGPL=ILLUMINA "
                    f"RGPU={tag} "
                    f"RGSM={tag} "
                    f"CREATE_INDEX=true "
                    f"VALIDATION_STRINGENCY=LENIENT"
                )

            # 3b. 去重
            tee.write(f"    Picard MarkDuplicates: tmp rg → tmp dedup\n")
            run_cmd(
                f"picard MarkDuplicates "
                f"I={rg_bam_tmp} "
                f"O={dedup_bam_tmp} "
                f"M={metrics_tmp} "
                f"REMOVE_DUPLICATES=true "
                f"CREATE_INDEX=true "
                f"VALIDATION_STRINGENCY=LENIENT"
            )

            # ── 4. 把 dedup BAM 放到 sorted.bam 同级目录 ──
            #    用户要求: dedup 放在 "原文件夹" (即 sorted.bam 所在目录,
            #    通常是 ..), 不在 cwd (output dir).
            #    兼容: 如果 ../ 不可写 (如只读挂载卷), fallback 留在 cwd.
            target_dir = sorted_bam_dir  # '.' or '..'
            target_dedup = f"{target_dir}/{tag}.sorted.dedup.bam"
            target_metrics = f"{target_dir}/{tag}.sorted.dup_metrics.txt"

            tee.write(f"    Placing dedup in: {target_dir}/ "
                      f"(sorted.bam 同级目录)\n")

            def _safe_copy(src, dst):
                """Copy src → dst, 处理跨设备 / 权限失败. 失败返回 False."""
                try:
                    if os.path.exists(dst):
                        os.remove(dst)
                    shutil.copy(src, dst)
                    # 验证 size 一致
                    if (os.path.exists(dst) and
                            os.path.getsize(dst) == os.path.getsize(src)):
                        return True
                    # size 不一致 → 当作失败
                    try:
                        os.remove(dst)
                    except OSError:
                        pass
                    return False
                except (OSError, IOError):
                    return False

            def _final_path(cwd_name, target_dir_name):
                """决定 dedup BAM 最终路径: 优先 target_dir, fallback cwd."""
                target = f"{target_dir_name}/{cwd_name}"
                if _safe_copy(cwd_name, target):
                    return os.path.realpath(target), target
                # fallback: 留在 cwd
                return os.path.abspath(cwd_name), cwd_name

            # Step 4a: 先把 Picard 的 /tmp 输出搬到 cwd (可写)
            cwd_files = {
                f"{tag}.sorted.dedup.bam": dedup_bam_tmp,
                f"{tag}.sorted.dedup.bam.bai": dedup_bam_tmp + ".bai",
                f"{tag}.sorted.dup_metrics.txt": metrics_tmp,
            }
            for cwd_name, src in cwd_files.items():
                if not os.path.exists(src):
                    continue
                if os.path.exists(cwd_name):
                    try:
                        os.remove(cwd_name)
                    except OSError:
                        pass
                shutil.copy(src, cwd_name)
                try:
                    os.remove(src)
                except OSError:
                    pass

            # Step 4b: 把 cwd 里的文件搬到 sorted.bam 同级目录
            final_bam, final_bam_loc = _final_path(
                f"{tag}.sorted.dedup.bam", target_dir)
            final_bai, final_bai_loc = _final_path(
                f"{tag}.sorted.dedup.bam.bai", target_dir)
            final_metrics, final_metrics_loc = _final_path(
                f"{tag}.sorted.dup_metrics.txt", target_dir)
            tee.write(f"      BAM:      {final_bam}\n")
            tee.write(f"      BAI:      {final_bai}\n")
            tee.write(f"      Metrics:  {final_metrics}\n")

            # Step 4c: 清理 cwd 里的中间副本 (如果已搬到 target_dir)
            for cwd_name in cwd_files:
                if not os.path.exists(cwd_name):
                    continue
                # 如果 target_dir 也有同名文件, 且 cwd 这份是 __最后__ 唯一文件,
                # 留作 fallback; 否则可删.
                target = f"{target_dir}/{cwd_name}"
                if (os.path.exists(target) and
                        os.path.getsize(target) ==
                        os.path.getsize(cwd_name)):
                    try:
                        os.remove(cwd_name)
                    except OSError:
                        pass

            # ── 5. cwd 里建 symlink 方便下游相对路径访问 ──
            cwd_link = f"{tag}.sorted.dedup.bam"
            if not os.path.exists(cwd_link):
                os.symlink(final_bam, cwd_link)

            # ── 6. 清理临时 sorted.bam symlink ──
            tmp_sorted = f"{tag}.sorted.bam"
            if (os.path.islink(tmp_sorted) and
                    os.path.realpath(tmp_sorted) != sorted_bam):
                try:
                    os.remove(tmp_sorted)
                except OSError:
                    pass

            return final_bam

        tag_to_bam = {}
        for tag, desc in all_tags:
            tag_to_bam[tag] = _resolve_bam(tag, desc)

        # ── Helper: 构建 MACS3 -t / -c 参数 ────────────────────────────────
        def _bam_args(tags):
            return ' '.join(tag_to_bam[t] for t in tags)

        if chip_method == "bdgdiff":


            # ── Helper: 运行命令并同时捕获输出行 ─────────────────────────────
            def _run_and_capture(cmd):
                """运行 shell 命令，输出到 tee，返回所有输出行。"""
                proc = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, text=True)
                lines = []
                for line in proc.stdout:
                    tee.write(line)
                    lines.append(line.rstrip('\n'))
                ret = proc.wait()
                if ret != 0:
                    sys.exit(f"Command failed (exit={ret}): {cmd}")
                return lines

            # ── Step 1: MACS3 callpeak on Group 1 ────────────────────────────
            tee.write(f"\nMACS3 callpeak — Group 1: {group1_name}\n")
            cmd1 = (
                f"macs3 callpeak "
                f"-t {_bam_args(g1_ip_tags)} "
                f"-c {_bam_args(g1_input_tags)} "
                f"-f {fmt} -g {genome_size} -n {group1_name}"
            )
            if qvalue < 1:
                tee.write(f"  Q-value threshold: {qvalue}\n")
                cmd1 += f" -q {qvalue}"
            else:
                tee.write(f"  P-value threshold: {pvalue}\n")
                cmd1 += f" -p {pvalue}"
            cmd1 += " --bdg"
            out1 = _run_and_capture(cmd1)

            # ── Step 2: MACS3 callpeak on Group 2 ────────────────────────────
            tee.write(f"\nMACS3 callpeak — Group 2: {group2_name}\n")
            cmd2 = (
                f"macs3 callpeak "
                f"-t {_bam_args(g2_ip_tags)} "
                f"-c {_bam_args(g2_input_tags)} "
                f"-f {fmt} -g {genome_size} -n {group2_name}"
            )
            if qvalue < 1:
                cmd2 += f" -q {qvalue}"
            else:
                cmd2 += f" -p {pvalue}"
            cmd2 += " --bdg"
            out2 = _run_and_capture(cmd2)

            # ── Step 3: 从 MACS3 输出读取 d1/d2 ─────────────────────────────
            def _parse_control_fragments(lines):
                """从 MACS3 日志中提取 'fragments after filtering in control' 后的数字。"""
                for line in lines:
                    if 'fragments after filtering in control' in line:
                        parts = line.strip().split()
                        try:
                            return int(parts[-1].replace(',', ''))
                        except ValueError:
                            continue
                return None

            fragments1 = _parse_control_fragments(out1)
            fragments2 = _parse_control_fragments(out2)

            if fragments1 is None or fragments2 is None:
                sys.exit("Could not determine fragment counts from MACS3 output. "
                         "Make sure MACS3 callpeak ran successfully.")

            d1 = fragments1 / 1_000_000
            d2 = fragments2 / 1_000_000
            tee.write(f"\n  {group1_name} control fragments: {fragments1:,} "
                      f"({d1:.2f}M)\n")
            tee.write(f"  {group2_name} control fragments: {fragments2:,} "
                      f"({d2:.2f}M)\n")

            # ── Step 4: MACS3 bdgdiff ────────────────────────────────────────
            tee.write(f"\nMACS3 bdgdiff — {group1_name} vs {group2_name}\n")
            t1_bdg = f"{group1_name}_treat_pileup.bdg"
            c1_bdg = f"{group1_name}_control_lambda.bdg"
            t2_bdg = f"{group2_name}_treat_pileup.bdg"
            c2_bdg = f"{group2_name}_control_lambda.bdg"

            # 验证 bedGraph 文件存在
            for bdg_file in [t1_bdg, c1_bdg, t2_bdg, c2_bdg]:
                if not os.path.exists(bdg_file):
                    sys.exit(f"bedGraph file not found: {bdg_file}")

            diff_prefix = f"diff_{group1_name}_vs_{group2_name}"
            # MACS3 bdgdiff 的 -C 默认值为 0，输出文件名为 {prefix}_c{cutoff}_{cond}.bed
            diff_cutoff = opts.get('cutoff', 3)
            diff_cmd = (
                f"macs3 bdgdiff "
                f"--t1 {t1_bdg} --c1 {c1_bdg} "
                f"--t2 {t2_bdg} --c2 {c2_bdg} "
                f"--d1 {d1:.2f} --d2 {d2:.2f} "
                f"-C {diff_cutoff:.1f} "
                f"--o-prefix {diff_prefix}"
            )
            run_cmd(diff_cmd)

            # ── Peak QC ──────────────────────────────────────────────────────
            def _qc_peak_set(label, bed_file, ip_bam_tags, caller):
                """对一组 peak BED 文件进行质控并保存报告。"""
                if not os.path.exists(bed_file):
                    tee.write(f"  {bed_file} not found, skipping.\n")
                    return
                peak_count = 0
                total_bp = 0
                with open(bed_file) as f:
                    for line in f:
                        cols = line.strip().split('\t')
                        if len(cols) >= 3:
                            s, e = int(cols[1]), int(cols[2])
                            total_bp += e - s
                            peak_count += 1
                if peak_count == 0:
                    tee.write(f"  {bed_file}: no peaks.\n")
                    return
                avg = total_bp / peak_count

                # FRiP（每个生物学重复单独计算）
                frip_results = []
                total_reads_all = 0
                reads_in_all = 0
                for t in ip_bam_tags:
                    bam = tag_to_bam.get(t)
                    if not bam or not os.path.exists(bam):
                        continue
                    label_rep = t
                    res = subprocess.run(
                        ['samtools', 'view', '-c', bam],
                        capture_output=True, text=True)
                    try:
                        n_total = int(res.stdout.strip())
                    except ValueError:
                        continue
                    res2 = subprocess.run(
                        ['samtools', 'view', '-c', bam, '-L', bed_file],
                        capture_output=True, text=True)
                    try:
                        n_peaks = int(res2.stdout.strip())
                    except ValueError:
                        continue
                    frip_rep = n_peaks / n_total if n_total > 0 else 0
                    frip_results.append((label_rep, n_total, n_peaks, frip_rep))
                    total_reads_all += n_total
                    reads_in_all += n_peaks
                frip_pooled = reads_in_all / total_reads_all if total_reads_all > 0 else 0

                tee.write(f"\n  ── {label} ({caller}) ──\n")
                tee.write(f"    Peaks              : {peak_count}\n")
                tee.write(f"    Total length (bp)  : {total_bp:,}\n")
                tee.write(f"    Avg length         : {avg:.1f} bp\n")
                for lr, nt, np, fr in frip_results:
                    tee.write(f"\n    ── {lr} ──\n")
                    tee.write(f"      Total reads  : {nt:,}\n")
                    tee.write(f"      Reads in peaks: {np:,}\n")
                    tee.write(f"      FRiP         : {fr:.4f} ({fr*100:.2f}%)\n")
                tee.write(f"\n    ── Pooled ──\n")
                tee.write(f"      FRiP         : {frip_pooled:.4f} ({frip_pooled*100:.2f}%)\n")

                qc_file = f"{label}_peak_qc.txt"
                with open(qc_file, 'w') as f:
                    f.write(f"Sample\t{label}\n")
                    f.write(f"Peak_caller\t{caller}\n")
                    f.write(f"Total_peaks\t{peak_count}\n")
                    f.write(f"Total_peak_length_bp\t{total_bp}\n")
                    f.write(f"Average_peak_length_bp\t{avg:.1f}\n")
                    for lr, nt, np, fr in frip_results:
                        f.write(f"Replicate\t{lr}\n")
                        f.write(f"{lr}_total_reads\t{nt}\n")
                        f.write(f"{lr}_reads_in_peaks\t{np}\n")
                        f.write(f"{lr}_FRiP\t{fr:.4f}\n")
                    f.write(f"Pooled_FRiP\t{frip_pooled:.4f}\n")
                tee.write(f"    QC saved: {qc_file}\n")

            tee.write(f"\n{'='*60}\n")
            tee.write(f"Peak Quality Control\n")
            tee.write(f"{'='*60}\n")

            # 各组 callpeak 结果质控
            for tag, name, ip_tags in [
                    (f"{group1_name}_peaks.narrowPeak", group1_name, g1_ip_tags),
                    (f"{group2_name}_peaks.narrowPeak", group2_name, g2_ip_tags)]:
                if os.path.exists(tag):
                    _qc_peak_set(f"{name}_peaks", tag, ip_tags, 'MACS3')

            # 差异 peak 质控（实际文件含 cutoff 值: {prefix}_c{cutoff}_cond1.bed）
            diff_suffix = f"c{diff_cutoff:.1f}"
            for suffix, label, ip_tags in [
                    ('cond1', f'diff_{group1_name}_enriched', g1_ip_tags),
                    ('cond2', f'diff_{group2_name}_enriched', g2_ip_tags),
                    ('common', 'diff_common', g1_ip_tags + g2_ip_tags)]:
                bed = f"{diff_prefix}_{diff_suffix}_{suffix}.bed"
                _qc_peak_set(label, bed, ip_tags, 'MACS3_bdgdiff')

            # ── ChIPseeker Peak Annotation & GO Enrichment (逐 peak 文件调用) ──
            tee.write(f"\n{'='*60}\n")
            tee.write("ChIPseeker Peak Annotation & GO Enrichment\n")
            tee.write(f"{'='*60}\n")
            gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
            if not os.path.exists(gff_path):
                tee.write(f"  GFF not found: {gff_path}, skipping annotation.\n")
            else:
                chipseeker_r = os.path.join(prefix, "scripts", "chipseeker.R")
                peak_calls = [
                    (f"{group1_name}_peaks.narrowPeak", f"{group1_name}_peaks"),
                    (f"{group2_name}_peaks.narrowPeak", f"{group2_name}_peaks"),
                    (f"{diff_prefix}_{diff_suffix}_cond1.bed", f"diff_{group1_name}_enriched"),
                    (f"{diff_prefix}_{diff_suffix}_cond2.bed", f"diff_{group2_name}_enriched"),
                    (f"{diff_prefix}_{diff_suffix}_common.bed", "diff_common"),
                ]
                for peak_file, peak_name in peak_calls:
                    if os.path.exists(peak_file):
                        run_cmd(
                            f"Rscript --vanilla {chipseeker_r} "
                            f"{genome} {prefix} {peak_file} {peak_name} {tss_distance}"
                        )

            # ── Cleanup bedGraph 文件和替身 BAM ─────────────────────────────
            for pat in ["*_treat_pileup.bdg", "*_control_lambda.bdg"]:
                for fname in globmod.glob(pat):
                    os.unlink(fname)
            for tag, _ in all_tags:
                bam = f"{tag}.sorted.bam"
                if os.path.islink(bam):
                    os.unlink(bam)
                bam_dedup = f"{tag}.sorted.dedup.bam"
                if os.path.islink(bam_dedup):
                    os.unlink(bam_dedup)

            tee.write(f"\nDifferential peaks output:\n")
            tee.write(f"  {diff_prefix}_{diff_suffix}_cond1.bed  (enriched in {group1_name})\n")
            tee.write(f"  {diff_prefix}_{diff_suffix}_cond2.bed  (enriched in {group2_name})\n")
            tee.write(f"  {diff_prefix}_{diff_suffix}_common.bed (common peaks)\n")
            tee.write(f"ChIP-seq bdgdiff analysis completed!\n")

        elif chip_method == "diffbind":
            # DiffBind differential binding analysis
            if len(g1_ip_tags) < 2 or len(g2_ip_tags) < 2:
                sys.exit("DiffBind requires >= 2 biological replicates per group "
                         "for IP samples. Use --chip-method bdgdiff instead.")
            if len(g1_input_tags) < 2 or len(g2_input_tags) < 2:
                sys.exit("DiffBind requires >= 2 biological replicates per group "
                         "for Input samples.")

            tee.write("\n  Using DiffBind (DESeq2-based differential binding)\n")

            # ── MACS3 callpeak 每个样本单独 call ──
            tee.write("\n  Step 1: MACS3 callpeak (per-sample)\n")

            diffbind_dir = "diffbind_results"
            os.makedirs(diffbind_dir, exist_ok=True)

            def _callpeak_one(tag, ip_bam, input_bam):
                peak = f"{tag}_peaks.narrowPeak"
                if os.path.exists(peak):
                    tee.write(f"    {tag}: peaks exist, skipping\n")
                    return peak
                cmd = (
                    f"macs3 callpeak -t {ip_bam} -c {input_bam} "
                    f"-f {fmt} -g {genome_size} -n {tag}"
                )
                if qvalue < 1:
                    cmd += f" -q {qvalue}"
                else:
                    cmd += f" -p {pvalue}"
                run_cmd(cmd)
                return peak

            # 按位置配对 Input→IP
            def _pairs(input_tags, ip_tags):
                n = min(len(input_tags), len(ip_tags))
                return [(ip_tags[i], tag_to_bam[ip_tags[i]], tag_to_bam[input_tags[i]])
                         for i in range(n)]

            # all_pairs 格式: (input_tag, ip_tag, ip_bam, input_bam)
            all_pairs = (
                [(in_t, ipt, tag_to_bam[ipt], tag_to_bam[in_t])
                 for in_t, ipt in zip(g1_input_tags, g1_ip_tags)]
                + [(in_t, ipt, tag_to_bam[ipt], tag_to_bam[in_t])
                    for in_t, ipt in zip(g2_input_tags, g2_ip_tags)]
            )
            for _in_t, ipt, ip_b, inp_b in all_pairs:
                _callpeak_one(ipt, ip_b, inp_b)

            # ── 构建 sample sheet ──
            tee.write("\n  Step 2: Building DiffBind sample sheet\n")
            ss_file = os.path.join(diffbind_dir, "samplesheet.csv")
            rows = ["SampleID,Condition,Factor,Replicate,bamReads,bamControl,Peaks"]
            for _in_t, ipt, ip_b, inp_b in all_pairs:
                # tag = WT_IP_1 → group=WT, rep=1
                parts = ipt.rsplit('_', 1)
                rep = parts[1]
                remainder = parts[0]
                if remainder.endswith('_IP'):
                    group = remainder[:-3]
                elif remainder.endswith('_Input'):
                    group = remainder[:-6]
                else:
                    group = remainder
                peak_path = os.path.abspath(f"{ipt}_peaks.narrowPeak")
                rows.append(f"{ipt},{group},IP,{rep},{os.path.abspath(ip_b)},{os.path.abspath(inp_b)},{peak_path}")
            with open(ss_file, 'w') as f:
                f.write('\n'.join(rows) + '\n')
            tee.write(f"    Sample sheet: {ss_file}\n")

            # ── normalization factors ──
            # 文献 (mito/chloro 归一化) 的做法:
            #   1. 计算每个 IP 的 organelle reads (mt 或 chloro)
            #   2. scaling_factor = IP_organelle / max(IP_organelle across all IPs)
            #   3. normalized_IP_count = total_mapped_IP × scaling_factor
            #   4. Input 完全不修正 (ratios 接近一致, 用于 DESeq2 中)
            #   5. DESeq2 size factors 改用 IP 列做 median-of-ratios
            #   norm=total 时: 不算 scaling, 直接用总 mapped reads
            nf_file = None
            organelle_map = {'mito': 'mito', 'chloro': 'chloro',
                                'rdna': 'rdna',
                                'total': False}
            organelle_mode = organelle_map.get(chip_norm, False)

            if organelle_mode is not False or chip_norm == 'total':
                tee.write(f"\n  Step 3: Norm factors ({chip_norm})\n")
                nf_file = os.path.join(diffbind_dir, "norm_factors.tsv")
                # TSV 列:
                #   sample   - IP 样本名
                #   ip_total - IP 总 mapped reads
                #   ip_org   - IP 的 organelle reads (mito/chloro)
                #   ip_scale - IP 缩放因子 (mito/chloro 时)
                #   ip_norm  - 归一化后 IP count (total 时 = ip_total;
                #             mito/chloro 时 = ip_total × ip_scale)
                #   input_total - Input 总 mapped reads (不修正)
                with open(nf_file, 'w') as nf:
                    nf.write("sample\tip_total\tip_org\t"
                             "ip_scale\tip_norm\tinput_total\n")

                    # Pass 1: 收集 IP organelle reads, 找 max
                    ip_records = []
                    ip_orgs = []
                    for _in_t, ipt, ip_b, inp_b in all_pairs:
                        ip_total = _count_bam_reads(ip_b, organelle_mode=False)
                        if organelle_mode:
                            ip_org = _count_bam_reads(ip_b,
                                                      organelle_mode=organelle_mode)
                        else:
                            ip_org = ip_total  # total mode: scaling=1
                        in_total = _count_bam_reads(inp_b, organelle_mode=False)
                        ip_records.append((ipt, ip_total, ip_org, in_total,
                                           ip_b, inp_b))
                        ip_orgs.append(ip_org if organelle_mode else 1)

                    # scaling factor = IP_organelle / max(IP_organelle across all IPs)
                    if organelle_mode and any(o > 0 for o in ip_orgs):
                        max_ip_org = max(ip_orgs)
                    else:
                        max_ip_org = 1.0  # 全 0 或 total mode → scaling=1

                    # Pass 2: 写文件
                    for (ipt, ip_total, ip_org, in_total, ip_b, inp_b), \
                            denom in zip(ip_records, ip_orgs):
                        scale = (denom / max_ip_org) if max_ip_org > 0 else 1.0
                        ip_norm = int(round(ip_total * scale))
                        label = chip_norm
                        tee.write(f"    {ipt}: {label} "
                                  f"ip_total={ip_total:,} "
                                  f"ip_org={ip_org:,} "
                                  f"scale={scale:.3f} "
                                  f"ip_norm={ip_norm:,} "
                                  f"input_total={in_total:,} (raw, no scaling)\n")
                        nf.write(f"{ipt}\t{ip_total}\t{ip_org}\t"
                                 f"{scale:.6f}\t{ip_norm}\t{in_total}\n")
                tee.write(f"    Norm factors: {nf_file}\n")

            # ── 生成 BW 文件 (Step 3.5) ──
            # 用 dedup BAM → BW.
            #   IP samples  -- 用 norm_factors.tsv 里的 ip_scale 缩放
            #                  (mito/chloro 时 < 1; total/deseq2 时 = 1.0)
            #   Input samples -- 不缩放 (ratios across inputs 接近一致, 文献方法)
            # 仅当有 norm_factors.tsv 时才生成 BW (即 chip_norm != deseq2)
            if nf_file and os.path.exists(nf_file):
                tee.write("\n  Step 3.5: Generating BigWig files "
                          "from dedup BAM\n")
                bw_dir = os.path.join(diffbind_dir, "bw")
                os.makedirs(bw_dir, exist_ok=True)

                # 读 ip_scale (按 sample 索引)
                scale_map = {}
                with open(nf_file) as f:
                    next(f)  # skip header
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            scale_map[parts[0]] = float(parts[3])

                # IP BW: bw_scale = 1 / ip_scale (DESeq2 sizeFactor 的倒数)
                #   DESeq2 normalized = raw / sizeFactor = raw / ip_scale
                #   bamCoverage norm  = coverage × bw_scale = coverage × (1/ip_scale)
                #   → 两边数值完全一致!
                #
                #   生物学含义:
                #     WT_1 mito少 → ip_scale=0.330 → bw_scale=3.03 → BW 信号 ×3.03 ↑
                #     ago1_27 mito多 → ip_scale=1.000 → bw_scale=1.000 → 不变
                #     → IGV 里看到 WT peak 比 ago1_27 高 ✓ (符合生物学预期)
                for _in_t, ipt, ip_b, inp_b in all_pairs:
                    scale = scale_map.get(ipt, 1.0)   # ip_scale = ip_org/max
                    bw_scale = 1.0 / scale if scale > 0 else 1.0  # DESeq2 sizeFactor 的倒数
                    bw_out = os.path.join(bw_dir, f"{ipt}.bw")
                    ip_real = os.path.realpath(ip_b)
                    ip_bai = ip_real + ".bai"
                    # 显式 index (兼容 deeptools/pysam)
                    if not os.path.exists(ip_bai) or \
                            os.path.getmtime(ip_bai) < \
                            os.path.getmtime(ip_real):
                        run_cmd(f"samtools index {ip_real}")
                    tee.write(f"    IP {ipt}: nf_scale={scale:.3f} "
                              f"bw_scale={bw_scale:.3f} → {bw_out}\n")
                    run_cmd(
                        f"bamCoverage "
                        f"--bam {ip_real} "
                        f"--outFileName {bw_out} "
                        f"--outFileFormat bigwig "
                        f"--scaleFactor {bw_scale:.6f} "
                        f"--normalizeUsing None "
                        f"--binSize 10 "
                        f"--extendReads 200 "
                        f"--numberOfProcessors 4"
                    )

                # Input BW: 不缩放 (scaleFactor=1.0)
                for in_t, ipt, ip_b, inp_b in all_pairs:
                    bw_out = os.path.join(bw_dir, f"{in_t}.bw")
                    in_real = os.path.realpath(inp_b)
                    in_bai = in_real + ".bai"
                    # 显式 index
                    if not os.path.exists(in_bai) or \
                            os.path.getmtime(in_bai) < \
                            os.path.getmtime(in_real):
                        run_cmd(f"samtools index {in_real}")
                    tee.write(f"    Input {in_t}: no scaling → {bw_out}\n")
                    run_cmd(
                        f"bamCoverage "
                        f"--bam {in_real} "
                        f"--outFileName {bw_out} "
                        f"--outFileFormat bigwig "
                        f"--scaleFactor 1.0 "
                        f"--normalizeUsing None "
                        f"--binSize 10 "
                        f"--extendReads 200 "
                        f"--numberOfProcessors 4"
                    )

                tee.write(f"    BigWig files: {bw_dir}\n")

            # ── Step 3.8: 手动从所有 BAM 统计 peak counts ──
            # DiffBind 3.x 的 dba.count() 把 Input reads 合并到 IP,
            # dual_factor 模型无法估计 IP vs Input 交互项.
            # 所以用 bedtools multicov 从所有 BAM (IP + Input) 直接数
            # DiffBind consensus peaks 的 reads.
            # 先让 DiffBind 只生成 consensus peaks (不 count)
            tee.write("\n  Step 3.8: bedtools multicov from all BAMs "
                      "(manual count for dual_factor)\n")

            # 3.8a) 让 DiffBind 用 dba + dba.blacklist + dba.count 生成 consensus peaks
            #        (dba.count 会按 minOverlap=2 合并 peaks,
            #        但我们不使用它的 counts —— 只用 $merged peak 位置)
            # 用一个小 R 段生成 consensus peak BED:
            peaks_bed = os.path.join(diffbind_dir, "consensus_peaks.bed")
            if not os.path.exists(peaks_bed):
                # 写临时 R script 生成 consensus peaks
                tmp_r = os.path.join(diffbind_dir, "_make_peaks.R")
                with open(tmp_r, 'w') as tr:
                    tr.write(f'''
suppressPackageStartupMessages(library(DiffBind))
args <- commandArgs(trailingOnly=TRUE)
ss_file <- args[1]; out_bed <- args[2]
samples <- read.csv(ss_file, stringsAsFactors=FALSE, check.names=FALSE)
peakcat <- dba(sampleSheet=samples, peakCaller="macs", peakFormat="narrow",
    config=data.frame(AnalysisMethod=DBA_DESEQ2, th=0.05,
                      DataType=DBA_DATA_GRANGES, RunParallel=TRUE,
                      minQCth=15, fragmentSize=125, bCorPlot=FALSE,
                      reportInit="DBA", bUsePval=FALSE, design=TRUE,
                      doBlacklist=FALSE, doGreylist=FALSE))
# 用 dba.count 生成 consensus (minOverlap=2, 但丢弃它的 counts)
peakcat <- dba.count(peakcat, minOverlap=2, score=DBA_SCORE_NORMALIZED)
# peakcat$merged = matrix [CHR_idx, START, END]
# 需要 chrmap 把 idx 转回 chr 名字
chrmap <- peakcat$chrmap
merged <- peakcat$merged
bed <- data.frame(
    chr = chrmap[merged[, 1]],
    start = merged[, 2] - 1,
    end   = merged[, 3],
    name  = paste0("peak_", seq_len(nrow(merged))),
    score = ".",
    strand = "."
)
write.table(bed, out_bed, sep="\\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
cat("Consensus peaks:", nrow(bed), "\\n")
''')
                run_cmd(f"Rscript --vanilla {tmp_r} {ss_file} {peaks_bed}")
            else:
                tee.write(f"    Consensus peaks exist, skipping: {peaks_bed}\n")

            # 3.8b) 收集所有 BAM (IP + Input)
            all_bams = []       # [(sample_label, bam_path), ...]
            sample_labels = []  # 保持顺序
            for in_t, ipt, ip_b, inp_b in all_pairs:
                all_bams.append((ipt, os.path.realpath(ip_b)))
                all_bams.append((in_t, os.path.realpath(inp_b)))
                # 保证 .bai 存在
                for _, br in [(ipt, os.path.realpath(ip_b)),
                              (in_t, os.path.realpath(inp_b))]:
                    bai = br + ".bai"
                    if not os.path.exists(bai) or \
                            os.path.getmtime(bai) < os.path.getmtime(br):
                        run_cmd(f"samtools index {br}")

            # 3.8c) bedtools multicov
            # 输入: BED + 所有 indexed BAM
            # 输出 TSV: chrom,start,end,name,score,strand,count1,count2,...
            counts_tsv = os.path.join(diffbind_dir, "peak_counts.tsv")
            bam_list_str = ' '.join(f'"{b}"' for _, b in all_bams)
            tee.write(f"    Counting peaks in {len(all_bams)} BAMs...\n")
            run_cmd(
                f"bedtools multicov -bams {bam_list_str} -bed {peaks_bed} "
                f"-q 10 > {counts_tsv}"
            )
            # 检查 bedtools multicov 是否可用 (有些版本用 -bedFiles)
            import subprocess as _sp
            # 如果上面命令失败, 尝试 -bedFiles 替代 -bed
            if not os.path.exists(counts_tsv) or os.path.getsize(counts_tsv) == 0:
                tee.write("    [FALLBACK] trying bedtools multicov -bedFiles ...\n")
                run_cmd(
                    f"bedtools multicov -bams {bam_list_str} -bedFiles {peaks_bed} "
                    f"-q 10 > {counts_tsv}"
                )
            # 如果还是不行, 输出错误诊断
            if not os.path.exists(counts_tsv) or os.path.getsize(counts_tsv) == 0:
                run_cmd("bedtools multicov --help 2>&1 | head -20")
                sys.exit("bedtools multicov failed; check bedtools version "
                         "(需要 >= 2.26)")

            # 3.8d) 把 TSV 头加回去 (bedtools 不写 header)
            header_line = '\t'.join(
                ["chr", "start", "end", "name", "score", "strand"] +
                [s for s, _ in all_bams]
            )
            with open(counts_tsv) as f:
                body = f.read()
            with open(counts_tsv, 'w') as f:
                f.write(header_line + '\n' + body)

            tee.write(f"    Peak counts: {counts_tsv} "
                      f"({sum(1 for _ in open(counts_tsv))-1} peaks x "
                      f"{len(all_bams)} samples)\n")

            # 3.8e) 写 sample_metadata.tsv (直接给 R 读, 不用 regex 猜)
            meta_tsv = os.path.join(diffbind_dir, "sample_metadata.tsv")
            with open(meta_tsv, 'w') as mf:
                mf.write("sample\tcondition\tfactor\n")
                for in_t, ipt, _ip_b, _inp_b in all_pairs:
                    # 从 IP 标签提取 group (和 samplesheet 一致)
                    parts = ipt.rsplit('_', 1)
                    remainder = parts[0]
                    if remainder.endswith('_IP'):
                        group = remainder[:-3]
                    elif remainder.endswith('_Input'):
                        group = remainder[:-6]
                    else:
                        group = remainder
                    # IP 样本
                    mf.write(f"{ipt}\t{group}\tIP\n")
                    # Input 样本 (沿用 IP 的 condition)
                    mf.write(f"{in_t}\t{group}\tInput\n")
            tee.write(f"    Sample metadata: {meta_tsv}\n")

            # 3.8f) 清理临时文件
            tmp_r_clean = os.path.join(diffbind_dir, "_make_peaks.R")
            if os.path.exists(tmp_r_clean):
                os.remove(tmp_r_clean)

            # ── 跑 DiffBind Rscript ──
            tee.write("\n  Step 4: Running DiffBind\n")
            tee.write(f"    analysis={chip_analysis}, norm={chip_norm}\n")
            dbr = os.path.join(prefix, "scripts", "chip_diffbind.R")
            if not os.path.exists(dbr):
                sys.exit(f"DiffBind R script not found: {dbr}")

            cmd_parts = [
                "Rscript", "--vanilla", dbr,
                ss_file, diffbind_dir, genome,
                chip_analysis, chip_norm,
                str(pvalue), str(foldchange),
                nf_file or "",
                os.path.join(diffbind_dir, "peak_counts.tsv"),
                os.path.join(diffbind_dir, "sample_metadata.tsv"),
            ]
            run_cmd(' '.join(cmd_parts))

            # ── ChIPseeker Peak Annotation & GO Enrichment (三类 peaks) ──
            tee.write(f"\n{'='*60}\n")
            tee.write("Step 5: ChIPseeker Peak Annotation & GO Enrichment\n")
            tee.write(f"{'='*60}\n")
            gff_path = os.path.join(prefix, "reference", f"{genome}_genes.gff")
            chipseeker_r = os.path.join(prefix, "scripts", "chipseeker.R")
            peak_sets = [
                # (BED file, label for output naming)
                ("DiffBind_consensus_peaks.bed", "DiffBind_consensus"),
                ("DiffBind_UP_peaks.bed",        "DiffBind_UP"),
                ("DiffBind_DOWN_peaks.bed",      "DiffBind_DOWN"),
            ]
            anno_ok = 0
            for bed_name, peak_label in peak_sets:
                bed_path = os.path.join(diffbind_dir, bed_name)
                if not os.path.exists(bed_path):
                    tee.write(f"  [SKIP] {bed_name} not found.\n")
                    continue
                # 如果 BED 是空的 (该类 peak 不存在)
                with open(bed_path) as _bf:
                    _n = sum(1 for _ in _bf)
                if _n == 0:
                    tee.write(f"  [SKIP] {peak_label}: 0 peaks.\n")
                    continue
                if not os.path.exists(gff_path):
                    tee.write(f"  [SKIP] GFF not found: {gff_path}\n")
                    break
                tee.write(f"\n  Annotating {peak_label} ({_n} peaks)...\n")
                try:
                    run_cmd(
                        f"Rscript --vanilla {chipseeker_r} "
                        f"{genome} {prefix} {bed_path} {peak_label} {tss_distance}"
                    )
                    anno_ok += 1
                    # 把 chipseeker 输出移到 diffbind_dir (默认写到 CWD)
                    for suf in ["_annotation.txt", "_annotation_pie.pdf",
                                "_go_enrichment.txt", "_go_dotplot.pdf"]:
                        _src = f"{peak_label}{suf}"
                        if os.path.exists(_src):
                            os.rename(_src, os.path.join(diffbind_dir, _src))
                except subprocess.CalledProcessError:
                    tee.write(f"  [WARN] Annotation failed for {peak_label}.\n")

            tee.write(f"\n  Annotated peak sets: {anno_ok}/{len(peak_sets)}\n")
            tee.write(f"  (Results in {diffbind_dir}/)\n")

            tee.write(f"\n{'='*60}\n")
            tee.write(f"DiffBind results: {diffbind_dir}/\n")
            tee.write(f"{'='*60}\n")

