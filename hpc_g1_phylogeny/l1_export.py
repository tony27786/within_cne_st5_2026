#!/usr/bin/env python3
"""l1_export.py -- L1 alignment and -fconst counts from pseudo-genomes. Needs numpy.

Candidate columns: every position of the analyzable region (minus --exclude-bed); a column is kept when at most floor(max_n_frac x n_taxa) taxa are N. Class i (all A/C/G/T identical) is counted in -fconst; class ii (one determined base, some N) and class iii (two or more bases) are written to the alignment.
Outputs: <prefix>.aln.fasta, <prefix>.fconst.txt, <prefix>.columns.tsv (chrom, pos, class, n_N, bases, in_R) and <prefix>.summary.tsv.
"""
import argparse

import numpy as np

from g1lib import Faidx, Intervals, die, log, read_genome, write_tsv

N, A, C, G, T = (ord(x) for x in 'NACGT')


def member(ivset, contig, pos0):
    """Vectorized membership of 0-based positions in an Intervals set."""
    iv = ivset.intervals(contig)
    if not iv:
        return np.zeros(len(pos0), dtype=bool)
    st = np.array([s for s, _ in iv]); en = np.array([e for _, e in iv])
    k = np.searchsorted(st, pos0, side='right') - 1
    return (k >= 0) & (pos0 < en[np.clip(k, 0, None)])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--taxa', required=True, help='TSV taxon<TAB>pseudo-genome FASTA (indexed), alignment order')
    ap.add_argument('--region', required=True)
    ap.add_argument('--genome', required=True, help='genome.txt (.fai order)')
    ap.add_argument('--max-n-frac', type=float, required=True)
    ap.add_argument('--exclude-bed', default='none')
    ap.add_argument('--r-bed', default='none')
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()

    taxa = []
    with open(a.taxa) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) == 2 and f[0]:
                taxa.append((f[0], f[1]))
    n = len(taxa)
    max_n = int(a.max_n_frac * n + 1e-9)
    region, excl, rset = Intervals.from_bed(a.region), Intervals.from_bed(a.exclude_bed), Intervals.from_bed(a.r_bed)
    fas = [Faidx(p) for _, p in taxa]

    fconst = {A: 0, C: 0, G: 0, T: 0}
    blocks, col_rows = [], []
    cnt = {'n_taxa': n, 'max_N_per_column': max_n, 'region_positions': 0, 'excluded_by_exclude_bed': 0, 'excluded_missing': 0,
           'class_i': 0, 'class_ii': 0, 'class_iii': 0, 'R_positions_in_region': 0, 'R_excluded_missing': 0,
           'R_in_class_i': 0, 'R_in_class_ii': 0, 'R_in_class_iii': 0}
    for contig, length in read_genome(a.genome):
        ivs = region.intervals(contig)
        if not ivs:
            continue
        cols = np.concatenate([np.arange(s, e) for s, e in ivs])
        cnt['region_positions'] += len(cols)
        drop = member(excl, contig, cols)
        cnt['excluded_by_exclude_bed'] += int(drop.sum())
        cols = cols[~drop]
        if len(cols) == 0:
            continue
        m = np.empty((n, len(cols)), dtype=np.uint8)
        for i, fa in enumerate(fas):
            seq = np.frombuffer(fa.fetch(contig), dtype=np.uint8)
            if len(seq) != length:
                die(f'{taxa[i][0]}: {contig} has {len(seq)} bases, genome says {length}')
            m[i] = seq[cols]
        other = ~np.isin(m, np.frombuffer(b'ACGTN', dtype=np.uint8))
        if other.any():
            die(f'{contig}: characters other than A/C/G/T/N in the pseudo-genomes')
        n_n = (m == N).sum(axis=0)
        present = np.stack([(m == b).any(axis=0) for b in (A, C, G, T)])
        n_b = present.sum(axis=0)
        inc = n_n <= max_n
        c1 = inc & (n_n == 0) & (n_b == 1)
        c2 = inc & (n_n > 0) & (n_b == 1)
        c3 = inc & (n_b >= 2)
        if (inc & (n_b == 0)).any():
            die(f'{contig}: an included column has no determined base (n_taxa {n}, max N {max_n})')
        cnt['excluded_missing'] += int((~inc).sum())
        cnt['class_i'] += int(c1.sum()); cnt['class_ii'] += int(c2.sum()); cnt['class_iii'] += int(c3.sum())
        for k, b in enumerate((A, C, G, T)):
            fconst[b] += int((c1 & present[k]).sum())
        in_r = member(rset, contig, cols)
        cnt['R_positions_in_region'] += int(in_r.sum())
        cnt['R_excluded_missing'] += int((in_r & ~inc).sum())
        cnt['R_in_class_i'] += int((in_r & c1).sum()); cnt['R_in_class_ii'] += int((in_r & c2).sum()); cnt['R_in_class_iii'] += int((in_r & c3).sum())
        sel = c2 | c3
        if sel.any():
            blocks.append(m[:, sel])
            letters = np.array(list('ACGT'))
            for j in np.nonzero(sel)[0]:
                col_rows.append({'chrom': contig, 'pos': int(cols[j]) + 1, 'class': 'iii' if c3[j] else 'ii', 'n_N': int(n_n[j]),
                                 'bases': ''.join(letters[present[:, j]]), 'in_R': 'yes' if in_r[j] else 'no'})
        log(f'    {contig}: {len(cols)} columns -> i {int(c1.sum())}, ii {int(c2.sum())}, iii {int(c3.sum())}, excluded (N > {max_n}) {int((~inc).sum())}')

    aln = np.concatenate(blocks, axis=1) if blocks else np.empty((n, 0), dtype=np.uint8)
    with open(f'{a.out_prefix}.aln.fasta', 'wb') as fh:
        for i, (name, _) in enumerate(taxa):
            fh.write(b'>' + name.encode() + b'\n' + aln[i].tobytes() + b'\n')
    fc = f'{fconst[A]},{fconst[C]},{fconst[G]},{fconst[T]}'
    with open(f'{a.out_prefix}.fconst.txt', 'w') as fh:
        fh.write(fc + '\n')
    write_tsv(f'{a.out_prefix}.columns.tsv', col_rows, ['chrom', 'pos', 'class', 'n_N', 'bases', 'in_R'])
    cnt['alignment_columns'] = aln.shape[1]
    cnt['fconst'] = fc
    cnt['R_in_L1_any_class'] = cnt['R_in_class_i'] + cnt['R_in_class_ii'] + cnt['R_in_class_iii']
    write_tsv(f'{a.out_prefix}.summary.tsv', [{'metric': k, 'value': v} for k, v in cnt.items()], ['metric', 'value'])
    log(f'    alignment {n} x {aln.shape[1]} (ii {cnt["class_ii"]}, iii {cnt["class_iii"]}); -fconst {fc}; R in L1: {cnt["R_in_L1_any_class"]}')


if __name__ == '__main__':
    main()
