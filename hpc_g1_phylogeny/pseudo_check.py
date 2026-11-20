#!/usr/bin/env python3
"""pseudo_check.py -- usable BED and stats of one pseudo-genome. Needs numpy.

Usable bases are the non-N runs of the pseudo-genome written by `bcftools consensus` (01_prep).
Inputs: --pseudo (FASTA with .fai) and --dump (bcftools query: CHROM POS REF ALT GT SRC at the base-1 sites).
Outputs: --out-usable (BED) and --out-stats (sample, usable_bp, n_base1_sites, n_gt_missing).
"""
import argparse

import numpy as np

from g1lib import Faidx, log, write_tsv

N = ord('N')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--pseudo', required=True)
    ap.add_argument('--dump', required=True)
    ap.add_argument('--out-usable', required=True)
    ap.add_argument('--out-stats', required=True)
    a = ap.parse_args()

    n_sites = n_missing = 0
    with open(a.dump) as fh:
        for line in fh:
            n_sites += 1
            if line.rstrip('\n').split('\t')[4] in ('.', './.', '.|.'):
                n_missing += 1

    pse = Faidx(a.pseudo)
    usable_rows = []
    usable_bp = 0
    for contig in pse.index:
        p = np.frombuffer(pse.fetch(contig), dtype=np.uint8)
        obs_use = p != N
        d = np.diff(np.concatenate(([0], obs_use.view(np.int8), [0])))
        starts, ends = np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]
        usable_bp += int((ends - starts).sum())
        usable_rows.extend((contig, int(s), int(e)) for s, e in zip(starts, ends))
    usable_rows.sort(key=lambda x: (x[0], x[1]))
    with open(a.out_usable, 'w') as fh:
        for c, s, e in usable_rows:
            fh.write(f'{c}\t{s}\t{e}\n')
    write_tsv(a.out_stats, [{'sample': a.sample, 'usable_bp': usable_bp, 'n_base1_sites': n_sites,
                             'n_gt_missing': n_missing}], ['sample', 'usable_bp', 'n_base1_sites', 'n_gt_missing'])
    log(f'    {a.sample}: usable {usable_bp} bp')


if __name__ == '__main__':
    main()
