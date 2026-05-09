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
        from prnaseqtools.functions import _tee
        tee = _tee()
    except Exception:
        pass

    prefix = opts.get('prefix', str(Path(__file__).resolve().parent.parent))
    foldchange = opts.get('foldchange', 1.5)
    pvalue = opts.get('pvalue', 0.05)
    control = opts.get('control', '')
    treatment = opts.get('treatment', '')
    run_mode = opts.get('run_mode', 'mrna')
    norm = opts.get('norm', 'rRNA,total')
    norms = norm.split(',')
    binsize = opts.get('binsize', 100)
    genome = opts.get('genome', 'ath')
    deseq2_norm = opts.get('deseq2_norm', 'DESeq2')

    tee.write(f"Two-factor comparison between {control} and {treatment}\n"
              f"Fold change = {foldchange} P value = {pvalue}\n")

    # Parse control/treatment as genotype,time1,N1,time2,N2 format
    control_parts = control.split('=')
    treatment_parts = treatment.split('=')
    control_genotype = control_parts[0]
    treatment_genotype = treatment_parts[0]

    controls = control_parts[1].split(',') if len(control_parts) > 1 else []
    treatments = treatment_parts[1].split(',') if len(treatment_parts) > 1 else []

    if len(controls) != len(treatments):
        sys.exit("Please provide paired data!")

    par_list = [control_genotype] + controls + [treatment_genotype] + treatments
    par_str = ' '.join(par_list)

    tags = []
    for c in [iter(controls)] * 2:
        try:
            time_label = next(iter(controls))
            n_rep = int(next(iter(controls)))
            for n in range(1, n_rep + 1):
                tags.append(f"{control_genotype}_{time_label}_{n}")
        except StopIteration:
            pass

    for c in [iter(treatments)] * 2:
        try:
            time_label = next(iter(treatments))
            n_rep = int(next(iter(treatments)))
            for n in range(1, n_rep + 1):
                tags.append(f"{treatment_genotype}_{time_label}_{n}")
        except StopIteration:
            pass

    # Actually parse correctly (remove above and redo)
    tags = []
    c_list = list(controls)
    while len(c_list) >= 2:
        time_label = c_list.pop(0)
        n_rep = int(c_list.pop(0))
        for n in range(1, n_rep + 1):
            tags.append(f"{control_genotype}_{time_label}_{n}")

    t_list = list(treatments)
    while len(t_list) >= 2:
        time_label = t_list.pop(0)
        n_rep = int(t_list.pop(0))
        for n in range(1, n_rep + 1):
            tags.append(f"{treatment_genotype}_{time_label}_{n}")

    if run_mode == 'mrna':
        for pre in tags:
            os.symlink(f"../{pre}.txt", f"{pre}.txt")

        subprocess.run(
            f"Rscript --vanilla {prefix}/scripts/tf_mrna.R "
            f"{deseq2_norm} {pvalue} {foldchange} {par_str}",
            shell=True, check=True
        )

        for fname in globmod.glob("*_?.txt"):
            os.unlink(fname)

    elif run_mode == 'srna':
        for pre in tags:
            os.symlink(f"../{pre}.nf", f"{pre}.nf")
            for fname in os.listdir(".."):
                if fname.startswith(pre) and fname.endswith('count') and 'norm' not in fname:
                    os.symlink(f"../{fname}", fname)

        for mnorm in norms:
            subprocess.run(
                f"Rscript --vanilla {prefix}/scripts/tf_srna.R "
                f"{mnorm} {pvalue} {foldchange} {par_str}",
                shell=True, check=True
            )

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

            subprocess.run(
                f"Rscript --vanilla {prefix}/scripts/tf_mirna.R "
                f"{mnorm} {pvalue} {foldchange} {par_str}",
                shell=True, check=True
            )

        for fname in globmod.glob("*.count"):
            os.unlink(fname)
        for fname in globmod.glob("*.nf"):
            os.unlink(fname)
