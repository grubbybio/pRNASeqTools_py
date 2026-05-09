"""
Input file parser for pRNASeqTools.
Parses control/treatment sample specifications into tags, file paths, and replicate counts.
Mirrors the Perl input.pm module.
"""

import os
import re
from pathlib import Path


def parse_input(sample_dict):
    """
    Parse a sample dictionary {name: 'file1+file2+file3'} or {name: '3'}.

    Returns:
        tags: list of sample tags like ['control_1', 'control_2', ...]
        files: list of file paths
        pars: list of [name, num_replicates, ...] for R script use
    """
    tags = []
    files = []
    pars = []

    for input_tag, input_value in sample_dict.items():
        if re.match(r'^\d', input_tag):
            raise ValueError("Please do not name the group with a numeric initial")

        pars.append(input_tag)

        input_files = input_value.split('+')

        if input_files[0] == '':
            raise ValueError("Please provide enough biological replicates")

        if not re.match(r'^\d+$', input_files[0]):
            # File paths provided
            for i in range(1, len(input_files) + 1):
                fpath = input_files[i - 1]

                if not fpath.startswith('/'):
                    if ',' not in fpath:
                        if fpath.startswith('~/'):
                            fpath = os.path.expanduser(fpath)
                        else:
                            fpath = os.path.abspath(os.path.join('..', fpath))
                    else:
                        f1, f2 = fpath.split(',')
                        if f1.startswith('~/'):
                            f1 = os.path.expanduser(f1)
                        else:
                            f1 = os.path.abspath(os.path.join('..', f1))
                        if f2.startswith('~/'):
                            f2 = os.path.expanduser(f2)
                        else:
                            f2 = os.path.abspath(os.path.join('..', f2))
                        fpath = f"{f1},{f2}"

                tags.append(f"{input_tag}_{i}")
                files.append(fpath)

            pars.append(str(len(input_files)))
        else:
            # Just replicate count
            pars.append(input_files[0])
            for i in range(1, int(input_files[0]) + 1):
                tags.append(f"{input_tag}_{i}")

    return tags, files, pars
