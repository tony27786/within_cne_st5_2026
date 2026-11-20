#!/usr/bin/env python3
"""pca_kin.py -- PCA and kinship matrix K from the L1 class-iii columns. Needs numpy, scipy.

Each class-iii column gives one 0/1 column per non-major base (N imputed with the column mean), standardized to W; PCA is the SVD of W and K = W W^T / m (haploid VanRaden), written in the square TSV format pyseer reads with --similarities.
Also reports the Spearman rho of PC1-PC3 with each taxon's usable rate (usable bp / analyzable bp) and the off-diagonal distribution of K.
"""
import argparse

import numpy as np
from scipy import stats

from g1lib import log, read_fasta_alignment, read_tsv, write_tsv

LET = np.frombuffer(b'ACGT', dtype=np.uint8)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--aln', required=True)
    ap.add_argument('--columns', required=True)
    ap.add_argument('--usable-stats', required=True)
    ap.add_argument('--qc-gate', required=True)
    ap.add_argument('--analyzable-bp', type=int, required=True)
    ap.add_argument('--n-pcs', type=int, default=10)
    ap.add_argument('--rho-max', type=float, default=0.20)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()

    names, seqs = read_fasta_alignment(a.aln)
    cls = [r['class'] for r in read_tsv(a.columns)]
    m = np.stack([np.frombuffer(seqs[n], dtype=np.uint8) for n in names])
    if m.shape[1] != len(cls):
        raise SystemExit(f'[ERROR] alignment has {m.shape[1]} columns, columns table {len(cls)}')
    m = m[:, np.array([c == 'iii' for c in cls], dtype=bool)]
    cols = []
    for j in range(m.shape[1]):
        col = m[:, j]
        det = col != ord('N')
        counts = np.array([(col == b).sum() for b in LET])
        major = int(np.argmax(counts))          # argmax takes the first maximum: A<C<G<T on ties
        for k in range(4):
            if k == major or counts[k] == 0:
                continue
            x = np.where(det, (col == LET[k]).astype(float), np.nan)
            p = np.nanmean(x)
            x = np.where(det, x, p)
            cols.append((x - p) / np.sqrt(p * (1 - p)))
    w = np.column_stack(cols)
    n_tax, n_var = w.shape
    u, s, _ = np.linalg.svd(w, full_matrices=False)
    scores = u * s
    var = s ** 2 / (s ** 2).sum()
    k_mat = w @ w.T / n_var

    npc = min(a.n_pcs, len(s))
    write_tsv(f'{a.out_dir}/pca_scores.tsv',
              [{'sample': n, **{f'PC{i + 1}': f'{scores[t, i]:.6f}' for i in range(npc)}} for t, n in enumerate(names)],
              ['sample'] + [f'PC{i + 1}' for i in range(npc)])
    write_tsv(f'{a.out_dir}/pca_variance.tsv', [{'PC': f'PC{i + 1}', 'explained_variance': f'{var[i]:.6f}'} for i in range(npc)],
              ['PC', 'explained_variance'])
    with open(f'{a.out_dir}/K.tsv', 'w') as fh:
        fh.write('sample\t' + '\t'.join(names) + '\n')
        for t, n in enumerate(names):
            fh.write(n + '\t' + '\t'.join(f'{v:.8f}' for v in k_mat[t]) + '\n')

    usable = {r['sample']: int(r['usable_bp']) / a.analyzable_bp for r in read_tsv(a.usable_stats)}
    b2 = {r['sample']: float(r['callable_frac_union_mask']) for r in read_tsv(a.qc_gate)}
    chk, warn = [], False
    for i in range(min(3, npc)):
        for label, d in (('final_usable_rate', usable), ('base2_callable_rate', b2)):
            rho, p = stats.spearmanr(scores[:, i], [d[n] for n in names])
            ok = abs(rho) < a.rho_max
            if label == 'final_usable_rate' and not ok:
                warn = True
            chk.append({'PC': f'PC{i + 1}', 'rate': label, 'spearman_rho': f'{rho:.4f}', 'p': f'{p:.3g}',
                        'abs_rho_lt_0.20': 'yes' if ok else 'NO'})
    write_tsv(f'{a.out_dir}/pc_vs_usable_rate.tsv', chk, ['PC', 'rate', 'spearman_rho', 'p', 'abs_rho_lt_0.20'])

    off = k_mat[np.triu_indices(n_tax, 1)]
    qs = (0, 1, 5, 25, 50, 75, 95, 99, 100)
    k_rows = [{'statistic': f'p{q}', 'value': f'{np.percentile(off, q):.6f}'} for q in qs]
    k_rows += [{'statistic': 'mean', 'value': f'{off.mean():.6f}'}, {'statistic': 'sd', 'value': f'{off.std():.6f}'},
               {'statistic': 'diag_mean', 'value': f'{np.diag(k_mat).mean():.6f}'},
               {'statistic': 'n_taxa', 'value': n_tax}, {'statistic': 'n_class_iii_columns', 'value': m.shape[1]},
               {'statistic': 'n_coded_variants', 'value': n_var}]
    write_tsv(f'{a.out_dir}/K_offdiag_summary.tsv', k_rows, ['statistic', 'value'])
    hist, edges = np.histogram(off, bins=50)
    write_tsv(f'{a.out_dir}/K_offdiag_hist.tsv', [{'bin_lo': f'{edges[i]:.6f}', 'bin_hi': f'{edges[i + 1]:.6f}', 'n_pairs': int(hist[i])}
                                                 for i in range(len(hist))], ['bin_lo', 'bin_hi', 'n_pairs'])
    log(f'    PCA/K: {n_tax} taxa x {n_var} coded variants ({m.shape[1]} class-iii columns); PC1-3 explain '
        f'{", ".join(f"{v:.3f}" for v in var[:3])}; K off-diagonal median {np.median(off):.4f}')
    if warn:
        log('    [WARN] |Spearman rho| >= 0.20 between a PC (1-3) and the final usable rate: handle missingness before interpreting the PCs')


if __name__ == '__main__':
    main()
