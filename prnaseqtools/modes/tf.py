"""
Two-factor DE analysis mode.
Runs inside srna/mrna output folders.
"""

import os
import sys
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.validate_options import validate_options
from prnaseqtools import reference as ref


def run(opts):
    """Main entry point for two-factor DE analysis."""
    opts = validate_options(opts)
    tee = sys.stderr

    try:
        from prnaseqtools.functions import _tee, run_cmd
        tee = _tee()
    except Exception:
        pass

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
        # ── ChIP-seq differential peak calling via MACS3 bdgdiff ──────────
        # 参数格式与 mrna/srna 统一:
        #   --control "WT=input,3,IP,3"
        #   --treatment "KO=input,2,IP,2"
        # BAM 文件: {group}_{label}_{rep}.sorted.bam
        # 例如: WT_input_1.sorted.bam, WT_IP_1.sorted.bam
        tee.write("ChIP-seq differential peak calling mode (MACS3 bdgdiff)\n")

        genome_size = opts.get('genome_size')
        if not genome_size:
            sys.exit("--genome-size is required for ChIP bdgdiff (e.g. 1.35e8 for ath)")

        qvalue = opts.get('qvalue', 1.0)
        pvalue = opts.get('pvalue', 0.05)
        seq_strategy = opts.get('seq_strategy', 'paired')
        fmt = "BAMPE" if seq_strategy == 'paired' else "BAM"

        # ── 解析 --control / --treatment ─────────────────────────────────
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

        # ── 检查 BAM 文件是否存在 ────────────────────────────────────────
        all_tags = []
        for t in g1_input_tags:
            all_tags.append((t, f"{group1_name} Input"))
        for t in g1_ip_tags:
            all_tags.append((t, f"{group1_name} IP"))
        for t in g2_input_tags:
            all_tags.append((t, f"{group2_name} Input"))
        for t in g2_ip_tags:
            all_tags.append((t, f"{group2_name} IP"))

        # ── 查找 BAM 文件（支持 .sorted.bam 或 .sorted.dedup.bam）─────────
        def _resolve_bam(tag, desc):
            """返回存在的 BAM 文件名，优先 .sorted.bam，其次 .sorted.dedup.bam。"""
            for suffix in ('.sorted.bam', '.sorted.dedup.bam'):
                bam = f"{tag}{suffix}"
                if os.path.exists(bam):
                    return bam
                if os.path.exists(os.path.join('..', bam)):
                    os.symlink(os.path.join('..', bam), bam)
                    return bam
            sys.exit(f"BAM file not found: {tag}.sorted.bam or "
                     f"{tag}.sorted.dedup.bam ({desc})")

        tag_to_bam = {}
        for tag, desc in all_tags:
            tag_to_bam[tag] = _resolve_bam(tag, desc)

        # ── Helper: 构建 MACS3 -t / -c 参数 ──────────────────────────────
        def _bam_args(tags):
            return ' '.join(tag_to_bam[t] for t in tags)

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
                    for i, p in enumerate(parts):
                        if 'fragments' in p and 'after' in parts[i+1:i+2]:
                            num_str = parts[-1].replace(',', '')
                            try:
                                return int(num_str)
                            except ValueError:
                                continue
                    # Fallback: take last field
                    try:
                        return int(parts[-1].replace(',', ''))
                    except ValueError:
                        pass
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
            f"-C {diff_cutoff} "
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
