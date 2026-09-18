"""
Alternative Splicing (AS) analysis mode.

Pipeline (基于 crescent 的 DAS + DTU, 但把 rMATS 换成了 MAJIQ):

  所有 majiq/voila 命令通过 `conda run -n majiq_academic` 跨 env 调用.

  1. [可选] STAR mapping → sorted BAM (如果输入 FASTQ)
  2. MAJIQ build
     - GFF3 annotation + sorted BAMs → splicegraph.sql + {sample}.majiq
     - 如果用户只有 GTF → gffread -T GFF3 转换
     - 生成 MAJIQ config.ini ([info] + [experiments] INI 格式)
  3. MAJIQ quantifier:
     - deltapsi : ctrl vs treatment 两组间差异 dPSI
     - heterogen: 多样本各自 PSI 独立比较 (单样品含重复场景)
  4. VOILA (MAJIQ 自带):
     - voila tsv     → 差异剪接事件 TSV (带 PSI / dPSI / padj)
     - voila modulize → splicing modules 分类
     - voila view     → 交互式浏览器可视化

Inputs: FASTQ 或 现有 sorted BAM
control = group 0, treatment = group 1, 2, ...

与 crescent 原流程的差异:
  - rMATS  → MAJIQ (支持 de novo LSV 检测, 更好处理植物等非模式生物)
  - 增加 heterogen 模式 (每个样本独立 PSI, 组间统计检验)
  - VOILA tsv / modulize 替代手动 R 可视化

命令行入口: prnaseqtools as [--classifier {deltapsi,heterogen,both}]
"""

import os
import sys
# MAJIQ 独立 conda env 名 (跨 env 调用)
MAJIQ_ENV_NAME = 'majiq_academic'


def _majiq_cmd(subcommand):
    """返回 conda run -n <env> <subcommand>, 用于跨 env 调 majiq/voila."""
    return f"conda run -n {MAJIQ_ENV_NAME} {subcommand}"
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.input_parser import (parse_input, _parse_to_dict,
                                        _resolve_path)
from prnaseqtools.functions import (_tee, run_cmd)


# ═══════════════════════════════════════════════════════════════════════════
# 主入口
# ═══════════════════════════════════════════════════════════════════════════

def run(opts):
    """主入口: 可变剪接 (MAJIQ + VOILA) 分析."""
    tee = _tee()

    thread      = opts.get('thread', 4)
    genome      = opts.get('genome', 'ath')
    prefix      = opts.get('prefix',
                           str(Path(__file__).resolve().parent.parent))
    run_mode    = opts.get('run_mode', 'whole')
    seq_strategy = opts.get('seq_strategy')
    genome_size = opts.get('genome_size', 10)

    classifier  = opts.get('classifier', 'deltapsi')  # deltapsi|heterogen|both
    strandness  = opts.get('strandness', 'none')      # none|forward|reverse
    min_exp     = opts.get('min_experiments', 1)     # build 时最低 junction 支持
    psi_thresh  = opts.get('psi_threshold', 0.05)    # dPSI 显著性阈值
    view        = opts.get('view', False)             # 是否启动 voila view
    majiq_dir   = opts.get('majiq_dir', None)         # MAJIQ 安装目录 (None=PATH)

    control     = opts.get('control', '')
    treatment   = opts.get('treatment')

    mapping     = run_mode in ('whole', 'mapping-only')
    do_build    = run_mode in ('whole', 'mapping-only', 'build-only')
    do_quantify = run_mode in ('whole', 'quantify-only')
    do_voila    = run_mode in ('whole', 'quantify-only', 'voila-only')

    # ── Reference paths ─────────────────────────────────────────────────
    ref_dir  = os.path.join(prefix, "reference")
    ref_gff  = os.path.join(ref_dir, f"{genome}_genes.gff")
    ref_fasta = os.path.join(ref_dir, f"{genome}_chr_all.fasta")
    ref_gtf  = os.path.join(ref_dir, f"{genome}_genes.gtf")

    # MAJIQ 需要 GFF3 格式 annotation
    majiq_gff3 = os.path.join(ref_dir, f"{genome}_genes_majiq.gff3")
    _ensure_gff3(ref_gff, ref_gtf, majiq_gff3, tee)

    # ── 输出目录 ──────────────────────────────────────────────────────
    work_dir     = os.path.join(prefix, f"results_as_{genome}")
    bam_dir      = os.path.join(work_dir, "bams")
    build_dir    = os.path.join(work_dir, "build")
    quant_dir    = os.path.join(work_dir, "quantify")
    voila_dir    = os.path.join(work_dir, "voila")
    for d in (work_dir, bam_dir, build_dir, quant_dir, voila_dir):
        os.makedirs(d, exist_ok=True)

    # ── 解析样本 ──────────────────────────────────────────────────────
    tee.write("\n═══ AS / MAJIQ 分析启动 ═══\n")
    control_dict = _parse_to_dict(control)
    tags, files, pars = parse_input(control_dict)
    group_map = {t: 0 for t in tags}   # control → group 0

    if treatment:
        for gi, t_item in enumerate(
            (treatment if isinstance(treatment, list) else [treatment]), start=1
        ):
            td = _parse_to_dict(t_item)
            t_tags, t_files, t_pars = parse_input(td)
            tags.extend(t_tags); files.extend(t_files); pars.extend(t_pars)
            for t in t_tags:
                group_map[t] = gi      # treatment → group 1, 2, ...

    tee.write(f"  样本数: {len(tags)}\n")
    for t, gi in group_map.items():
        tee.write(f"    {t} → group {gi}\n")

    # ── Stage 1: Mapping (STAR → sorted BAM) ─────────────────────────
    existing_bams = _find_sorted_bams(bam_dir, tags)
    if mapping and len(existing_bams) < len(tags):
        _do_mapping(tags, files, pars, ref_fasta, genome_size,
                    bam_dir, thread, seq_strategy, tee,
                    majiq_dir=majiq_dir, reuse_existing=existing_bams)
        existing_bams = _find_sorted_bams(bam_dir, tags)

    if run_mode == 'mapping-only':
        tee.write("\n[mapping-only] BAM ready → 停止\n")
        return

    if len(existing_bams) < len(tags):
        tee.write(f"\nERROR: 只有 {len(existing_bams)}/{len(tags)} 个样本有 sorted BAM\n")
        sys.exit(1)

    # ── Stage 2: MAJIQ build ─────────────────────────────────────────
    if do_build:
        _do_majiq_build(majiq_gff3, bam_dir, build_dir,
                        strandness, min_exp, tee, majiq_dir=majiq_dir)

    if run_mode == 'build-only':
        tee.write("\n[build-only] splicegraph + .majiq ready → 停止\n")
        return

    # ── Stage 3: MAJIQ quantify ──────────────────────────────────────

    # ── Stage 3: MAJIQ psi-coverage + psi + quantify ─────────────
    if do_quantify:
        splicegraph = os.path.join(build_dir, "splicegraph.sql")
        _do_psi_coverage(splicegraph, build_dir, quant_dir,
                         bam_dir, thread, tee)
        _do_psi(quant_dir, thread, tee)
        _do_majiq_quantify(group_map, quant_dir, splicegraph,
                           classifier, psi_thresh, tee, thread=thread,
                           majiq_dir=majiq_dir)

    if run_mode == 'quantify-only':
        tee.write("\n[quantify-only] .voila 文件 ready → 停止\n")
        return

    # ── Stage 4: VOILA ─────────────────────────────────────────────
    if do_voila:
        splicegraph = os.path.join(build_dir, "splicegraph.sql")
        if not os.path.exists(splicegraph):
            splicegraph = _find_splicegraph(build_dir)
        _do_voila(splicegraph, quant_dir, voila_dir,
                  view, tee, majiq_dir=majiq_dir)


# ═══════════════════════════════════════════════════════════════════════════
# MAJIQ helper: 确保 GFF3 annotation 存在
# ═══════════════════════════════════════════════════════════════════════════

def _ensure_gff3(ref_gff, ref_gtf, majiq_gff3, tee):
    """
    MAJIQ 只接受 GFF3 格式。
    如果用户提供的 ref_gff 已经是 .gff3 → 直接用；
    如果是 .gff 或 .gtf → 用 gffread 转 GFF3。
    """
    if os.path.exists(majiq_gff3):
        tee.write(f"  MAJIQ annotation: {majiq_gff3} (cached)\n")
        return majiq_gff3

    # 检查现成的 GFF3（不同扩展名）
    for ext in ('gff3', 'GFF3', 'gff.GFF3'):
        if os.path.exists(ref_gff + ext):
            return ref_gff + ext

    tee.write("  转换 annotation 为 GFF3 格式 (MAJIQ 要求)...\n")
    src = ref_gff if os.path.exists(ref_gff) else ref_gtf
    if not os.path.exists(src):
        tee.write(f"  ERROR: 找不到 annotation ({ref_gff} / {ref_gtf})\n")
        sys.exit(1)

    # 用 gffread 转 GFF3
    cmd = f"gffread {src} -o {majiq_gff3} -T GFF3"
    tee.write(f"    {cmd}\n")
    run_cmd(cmd, tee)
    tee.write(f"    ✓ {majiq_gff3}\n")
    return majiq_gff3


# ═══════════════════════════════════════════════════════════════════════════
# Stage 1: STAR mapping
# ═══════════════════════════════════════════════════════════════════════════

def _do_mapping(tags, files, pars, ref_fasta, genome_size,
                bam_dir, thread, seq_strategy, tee, majiq_dir=None,
                reuse_existing=None):
    """STAR → sorted BAM."""
    reuse_existing = reuse_existing or set()
    # 生成 STAR index (如果不存在)
    star_idx = os.path.join(os.path.dirname(ref_fasta), "STAR_idx_" +
                            os.path.basename(ref_fasta).split('_chr')[0])
    if not os.path.exists(os.path.join(star_idx, "SAindex")):
        tee.write("  构建 STAR index...\n")
        cmd = (f"STAR --runMode genomeGenerate "
               f"--genomeDir {star_idx} "
               f"--genomeFastaFiles {ref_fasta} "
               f"--runThreadN {thread} "
               f"--genomeSAindexNbases {genome_size} --sjdbOverhang 99 --limitGenomeGenerateRAM 64000000000")
        run_cmd(cmd, tee)

    for i, tag in enumerate(tags):
        bam_out = os.path.join(bam_dir, f"{tag}.sorted.bam")
        if bam_out in reuse_existing:
            tee.write(f"    [{tag}] 已有 BAM → 跳过\n")
            continue

        fq1, fq2 = files[i] if isinstance(files[i], (tuple, list)) \
            else (files[i], None)
        fq_arg = f"{fq1} {fq2}" if fq2 else fq1

        out_prefix = os.path.join(bam_dir, f"{tag}.")
        tee.write(f"    [{tag}] STAR mapping...\n")
        cmd = (f"STAR --genomeDir {star_idx} --seedSearchStartLmax 25 "
               f"--readFilesIn {fq_arg} "
               f"--readFilesCommand zcat "
               f"--runThreadN {thread} "
               f"--outSAMtype BAM SortedByCoordinate "
               f"--outSAMstrandField intronMotif "
               f"--outBAMsortingThreadN {thread} "
               f"--outFileNamePrefix {out_prefix}")
        if pars[i]:
            cmd += f" {pars[i]}"
        run_cmd(cmd, tee)

        # 重命名 STAR 默认输出的 sorted BAM
        default_bam = os.path.join(bam_dir, f"{tag}.Aligned.sortedByCoord.out.bam")
        final_bam   = os.path.join(bam_dir, f"{tag}.sorted.bam")
        if os.path.exists(default_bam) and not os.path.exists(final_bam):
            os.rename(default_bam, final_bam)
        # index
        run_cmd(f"samtools index {final_bam}", tee)


def _find_sorted_bams(bam_dir, tags):
    """返回已存在的 sorted BAM 集合."""
    found = set()
    for tag in tags:
        bam = os.path.join(bam_dir, f"{tag}.sorted.bam")
        if os.path.exists(bam):
            found.add(bam)
    return found


# ═══════════════════════════════════════════════════════════════════════════
# Stage 2: MAJIQ build
# ═══════════════════════════════════════════════════════════════════════════

def _do_majiq_build(gff3, bam_dir, build_dir, strandness, min_exp, tee,
                    majiq_dir=None):
    """
    MAJIQ build: GFF3 + sorted BAMs → splicegraph.sql + {sample}.majiq.

    先生成 config.ini (INI 格式, MAJIQ 要求):
      [info]
      bamdirs = <bam_dir>
      sjdirs  = <空或之前的 build_dir>
      genome  = pRNASeqTools_AS
      strandness = none

      [experiments]
      sample_A = sample_A
      sample_B = sample_B
      ...

    所有样本独立 build group (min-exp=1 → 每个样本自身足够支持).
    如果用户传了 >1 个重复的同一 treatment group,
    可以考虑把重复合到一个 group 里让 splicegraph 更稳健.
    """
    tee.write("\n═══ MAJIQ build ═══\n")
    os.makedirs(build_dir, exist_ok=True)

    config_ini = os.path.join(build_dir, "majiq_build_config.ini")

    # 扫描 bam_dir 里所有 .sorted.bam
    bam_files = sorted(globmod.glob(os.path.join(bam_dir, "*.sorted.bam")))
    if not bam_files:
        tee.write(f"  ERROR: {bam_dir} 里没有 .sorted.bam\n")
        sys.exit(1)

    # 提取样本 prefix (去掉 .sorted.bam)
    samples = [os.path.basename(b)[:-len(".sorted.bam")] for b in bam_files]

    # 生成 config.ini
    cfg = configparser.ConfigParser()
    cfg['info'] = {
        'bamdirs': bam_dir,
        'sjdirs':  build_dir,
        'genome':  'pRNASeqTools_as',
        'strandness': strandness,
    }
    # 每个样本独立 experiment (prefix = sample 名)
    exp_section = {}
    for s in samples:
        exp_section[s] = s
    cfg['experiments'] = exp_section

    with open(config_ini, 'w') as f:
        cfg.write(f)
    tee.write(f"  config.ini: {config_ini}\n")
    tee.write(f"  样本: {', '.join(samples)}\n")

    # 构建 splicegraph.sql + .majiq per sample
    cmd = (f"conda run -n majiq_academic majiq build {gff3} "
           f"-o {build_dir} "
           f"-c {config_ini} "
           f"--min-experiments {min_exp}")
    tee.write(f"  {cmd}\n")
    run_cmd(cmd, tee)

    # 验证输出
    splicegraph = os.path.join(build_dir, "splicegraph.sql")
    if not os.path.exists(splicegraph):
        tee.write(f"  ERROR: 没找到 {splicegraph}\n")
        sys.exit(1)
    majiq_out = globmod.glob(os.path.join(build_dir, "*.majiq"))
    tee.write(f"  ✓ splicegraph.sql + {len(majiq_out)} 个 .majiq 文件\n")




def _glob_sj(build_dir):
    """返回 build 输出的 .sj 文件字典 (新版 majiq build 输出 SJ)."""
    out = {}
    for f in sorted(globmod.glob(os.path.join(build_dir, "*.sj"))):
        sample = os.path.basename(f)[:-len(".sj")]
        out[sample] = f
    return out


def _find_splicegraph(build_dir):
    for f in globmod.glob(os.path.join(build_dir, "splicegraph*.sql")):
        return f
    return None


def _do_psi_coverage(splicegraph, build_dir, quant_dir, bam_dir, thread, tee):
    """
    majiq psi-coverage: 合并所有 .sj → 一个 psi-coverage 文件.
    新版 build 输出 .sj (旧版 .majiq).
    """
    sj_files = sorted(globmod.glob(os.path.join(build_dir, "*.sj")))
    if not sj_files:
        tee.write("  WARNING: 没有 .sj, 降级用 .majiq\n")
        sj_files = sorted(globmod.glob(os.path.join(build_dir, "*.majiq")))
    if not sj_files:
        tee.write(f"  ERROR: {build_dir} 里没有 .sj 或 .majiq\n")
        sys.exit(1)

    os.makedirs(quant_dir, exist_ok=True)
    psicov = os.path.join(quant_dir, "all.psicov")
    tee.write(f"  psi-coverage: {len(sj_files)} 个 .sj → {os.path.basename(psicov)}\n")
    cmd = (f"conda run -n majiq_academic majiq psi-coverage {splicegraph} {psicov} "
           f"{' '.join(sj_files)} "
           f"--nthreads {thread} --overwrite")
    tee.write(f"    {cmd}\n")
    run_cmd(cmd, tee)
    return psicov


def _do_psi(quant_dir, thread, tee):
    """majiq psi: PSI 量化 (默认每样本独立)."""
    psicov = os.path.join(quant_dir, "all.psicov")
    if not os.path.exists(psicov):
        tee.write(f"  ERROR: {psicov} 不存在\n")
        sys.exit(1)
    tee.write("  psi quantification...\n")
    cmd = f"conda run -n majiq_academic majiq psi {psicov} --nthreads {thread} --overwrite"
    tee.write(f"    {cmd}\n")
    run_cmd(cmd, tee)


def _do_majiq_quantify(group_map, quant_dir, splicegraph,
                       classifier, psi_thresh, tee, thread=4,
                       majiq_dir=None):
    """
    MAJIQ quantifier (新版): psi → deltapsi/heterogen.
    deltapsi/heterogen 参数: -psi1/-psi2 + -splicegraph + --output-voila.
    """
    tee.write("\n═══ MAJIQ quantify ═══\n")

    psi_files = sorted(globmod.glob(os.path.join(quant_dir, "*.psi")))
    if not psi_files:
        psi_files = sorted(globmod.glob(os.path.join(quant_dir, "*.zarr")))
    if not psi_files:
        tee.write(f"  ERROR: {quant_dir} 里没找到 psi 文件\n")
        sys.exit(1)

    # 按 group_id 分组
    by_group = {}
    for tag, gid in group_map.items():
        matches = [f for f in psi_files if tag in os.path.basename(f)]
        by_group.setdefault(gid, []).extend(matches)

    ctrl_psi = by_group.get(0, [])
    if not ctrl_psi:
        tee.write("  ERROR: group 0 (control) 无 psi 文件\n")
        sys.exit(1)

    ctrl_label = 'CTRL'
    for trt_gid in sorted(g for g in by_group if g != 0):
        trt_psi = by_group[trt_gid]
        if not trt_psi:
            continue
        trt_label = f"TREAT{trt_gid}"
        combo_dir = os.path.join(quant_dir, f"{ctrl_label}_vs_{trt_label}")
        os.makedirs(combo_dir, exist_ok=True)

        tee.write(f"\n  {ctrl_label} vs {trt_label} "
                  f"(ctrl={len(ctrl_psi)} psi, trt={len(trt_psi)} psi)\n")

        if classifier in ('deltapsi', 'both'):
            _run_deltapsi(ctrl_psi, trt_psi, ctrl_label, trt_label,
                          splicegraph, combo_dir, tee)

        if classifier in ('heterogen', 'both'):
            _run_heterogen(ctrl_psi, trt_psi, ctrl_label, trt_label,
                           splicegraph, combo_dir, tee)


def _run_deltapsi(grp1_psi, grp2_psi, n1, n2, splicegraph, out_dir, tee):
    """新版 majiq deltapsi: -psi1/-psi2 + -splicegraph."""
    tee.write(f"    [deltapsi] {n1} vs {n2}...\n")
    cmd = (f"conda run -n majiq_academic majiq deltapsi -psi1 {' '.join(grp1_psi)} "
           f"-psi2 {' '.join(grp2_psi)} "
           f"-n {n1} {n2} "
           f"--splicegraph {splicegraph} "
           f"--output-voila {out_dir}/deltapsi.voila "
           f"--output-tsv {out_dir}/deltapsi.tsv "
           f"--min-experiments 1")
    tee.write(f"      {cmd}\n")
    run_cmd(cmd, tee)
    tee.write(f"      ✓ deltapsi done\n")


def _run_heterogen(grp1_psi, grp2_psi, n1, n2, splicegraph, out_dir, tee):
    """新版 majiq heterogen: -psi1/-psi2 + -splicegraph."""
    tee.write(f"    [heterogen] {n1} vs {n2}...\n")
    cmd = (f"conda run -n majiq_academic majiq heterogen -psi1 {' '.join(grp1_psi)} "
           f"-psi2 {' '.join(grp2_psi)} "
           f"-n {n1} {n2} "
           f"--splicegraph {splicegraph} "
           f"--output-voila {out_dir}/heterogen.voila "
           f"--output-tsv {out_dir}/heterogen.tsv "
           f"--min-experiments 1")
    tee.write(f"      {cmd}\n")
    run_cmd(cmd, tee)
    tee.write(f"      ✓ heterogen done\n")

"""
Alternative Splicing (AS) analysis mode.

Pipeline (基于 crescent 的 DAS + DTU, 但把 rMATS 换成了 MAJIQ):

  所有 majiq/voila 命令通过 `conda run -n majiq_academic` 跨 env 调用.

  1. [可选] STAR mapping → sorted BAM (如果输入 FASTQ)
  2. MAJIQ build
     - GFF3 annotation + sorted BAMs → splicegraph.sql + {sample}.majiq
     - 如果用户只有 GTF → gffread -T GFF3 转换
     - 生成 MAJIQ config.ini ([info] + [experiments] INI 格式)
  3. MAJIQ quantifier:
     - deltapsi : ctrl vs treatment 两组间差异 dPSI
     - heterogen: 多样本各自 PSI 独立比较 (单样品含重复场景)
  4. VOILA (MAJIQ 自带):
     - voila tsv     → 差异剪接事件 TSV (带 PSI / dPSI / padj)
     - voila modulize → splicing modules 分类
     - voila view     → 交互式浏览器可视化

Inputs: FASTQ 或 现有 sorted BAM
control = group 0, treatment = group 1, 2, ...

与 crescent 原流程的差异:
  - rMATS  → MAJIQ (支持 de novo LSV 检测, 更好处理植物等非模式生物)
  - 增加 heterogen 模式 (每个样本独立 PSI, 组间统计检验)
  - VOILA tsv / modulize 替代手动 R 可视化

命令行入口: prnaseqtools as [--classifier {deltapsi,heterogen,both}]
"""

import os
import sys
# MAJIQ 独立 conda env 名 (跨 env 调用)
MAJIQ_ENV_NAME = 'majiq_academic'


def _majiq_cmd(subcommand):
    """返回 conda run -n <env> <subcommand>, 用于跨 env 调 majiq/voila."""
    return f"conda run -n {MAJIQ_ENV_NAME} {subcommand}"
import glob as globmod
import subprocess
from pathlib import Path

from prnaseqtools.input_parser import (parse_input, _parse_to_dict,
                                        _resolve_path)
from prnaseqtools.functions import (_tee, run_cmd)

def _do_voila(splicegraph, quant_dir, voila_dir, view, tee, majiq_dir=None):
    """
    voila tsv / modulize / view.

    对 quant_dir 下每个 ctrl_vs_XXX 子目录:
      voila tsv splicegraph.sql *.voila -f results.tsv --show-all
      voila modulize splicegraph.sql *.voila -d modules_dir --show-all
    如果用户指定 --view: voila view splicegraph.sql *.voila
    """
    tee.write("\n═══ VOILA ═══\n")
    if not splicegraph:
        tee.write("  ERROR: 找不到 splicegraph.sql\n")
        sys.exit(1)
    tee.write(f"  splicegraph: {splicegraph}\n")

    # 收集所有 quantifier 输出的 .voila 文件
    all_voila = sorted(globmod.glob(
        os.path.join(quant_dir, "**", "*.voila"), recursive=True))
    if not all_voila:
        tee.write("  WARNING: 没找到 .voila 文件\n")
        return

    tee.write(f"  .voila files ({len(all_voila)}): "
              + ", ".join(os.path.basename(v) for v in all_voila[:5])
              + ("..." if len(all_voila) > 5 else ""))

    # ── voila tsv (每个比较单独一个 TSV) ─────────────────────────
    for subdir in sorted(globmod.glob(os.path.join(quant_dir, "*_vs_*"))):
        voila_in = sorted(globmod.glob(os.path.join(subdir, "*.voila")))
        if not voila_in:
            continue
        label = os.path.basename(subdir)
        out_tsv = os.path.join(voila_dir, f"{label}_tsv.tsv")
        tee.write(f"    [tsv] {label}...\n")
        cmd = (f"conda run -n majiq_academic voila tsv {splicegraph} {' '.join(voila_in)} "
               f"-f {out_tsv} --show-all")
        run_cmd(cmd, tee)

        # ── voila modulize ───────────────────────────────────────
        mod_dir = os.path.join(voila_dir, f"{label}_modules")
        os.makedirs(mod_dir, exist_ok=True)
        tee.write(f"    [modulize] {label}...\n")
        cmd = (f"conda run -n majiq_academic voila modulize {splicegraph} {' '.join(voila_in)} "
               f"-d {mod_dir} --show-all")
        run_cmd(cmd, tee)

    # ── voila view (交互式浏览器) ───────────────────────────────
    if view:
        tee.write("\n  启动 VOILA view (Ctrl+C 停止)...\n")
        tee.write("  访问本地端口查看交互式剪接图谱.\n")
        cmd = (f"conda run -n majiq_academic voila view {splicegraph} {' '.join(all_voila)}")
        tee.write(f"    {cmd}\n")
        # 前台运行, 不返回 (让用户 Ctrl+C)
        try:
            subprocess.run(cmd, shell=True)
        except KeyboardInterrupt:
            tee.write("\n  VOILA view stopped.\n")


# ═══════════════════════════════════════════════════════════════════════════
# 直接运行入口 (测试 / CLI 分发用)
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # 允许 python -m prnaseqtools.modes.as --help 等简单调用
    print("AS / MAJIQ mode — 通过 'prnaseqtools as' 调用")
