#!/usr/bin/env python3
"""sens_compare.py -- descriptive comparison of the L1 outputs with and without R. Needs numpy, dendropy. No thresholds.

Trees: normalized Robinson-Foulds distance and supported splits (SH-aLRT / UFBoot minima) in both / only one tree. PCA: |r| of PC1-PC3. K: r of the off-diagonal elements.
"""
import argparse

import dendropy
import numpy as np
from dendropy.calculate import treecompare

from g1lib import read_tsv, write_tsv, log


def splits(tree, names, alrt_min, ufb_min):
    leaves = frozenset(names)
    below, out = {}, set()
    for nd in tree.postorder_node_iter():
        below[nd] = frozenset([nd.taxon.label]) if nd.is_leaf() else frozenset().union(*(below[c] for c in nd.child_nodes()))
        if nd.is_leaf() or nd.parent_node is None or not nd.label:
            continue
        f = nd.label.split('/')
        try:
            alrt, ufb = float(f[-2]), float(f[-1])
        except (ValueError, IndexError):
            continue
        if alrt >= alrt_min and ufb >= ufb_min:
            s = below[nd]
            out.add(frozenset(leaves - s) if names[0] in s else s)
    return out


def read_k(path):
    with open(path) as fh:
        names = fh.readline().rstrip('\n').split('\t')[1:]
        m = np.array([[float(x) for x in line.rstrip('\n').split('\t')[1:]] for line in fh])
    return names, m


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for k in ('tree', 'pca', 'kin'):
        ap.add_argument(f'--main-{k}', required=True)
        ap.add_argument(f'--sens-{k}', required=True)
    ap.add_argument('--alrt-min', type=float, required=True)
    ap.add_argument('--ufboot-min', type=float, required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(path=a.main_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True, rooting='force-unrooted')
    t2 = dendropy.Tree.get(path=a.sens_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True, rooting='force-unrooted')
    names = sorted(t.label for t in tns)
    s1, s2 = splits(t1, names, a.alrt_min, a.ufboot_min), splits(t2, names, a.alrt_min, a.ufboot_min)
    t1.encode_bipartitions(); t2.encode_bipartitions()
    rf = treecompare.symmetric_difference(t1, t2)
    rf_norm = rf / (2 * (len(names) - 3))

    p1 = {r['sample']: r for r in read_tsv(a.main_pca)}
    p2 = {r['sample']: r for r in read_tsv(a.sens_pca)}
    common = [s for s in p1 if s in p2]
    rows = [{'metric': 'RF_symmetric_difference', 'value': rf}, {'metric': 'RF_normalized', 'value': round(rf_norm, 6)},
            {'metric': 'supported_splits_main', 'value': len(s1)}, {'metric': 'supported_splits_R_excluded', 'value': len(s2)},
            {'metric': 'supported_splits_in_both', 'value': len(s1 & s2)}, {'metric': 'supported_only_main', 'value': len(s1 - s2)},
            {'metric': 'supported_only_R_excluded', 'value': len(s2 - s1)}]
    for i in (1, 2, 3):
        x = np.array([float(p1[s][f'PC{i}']) for s in common]); y = np.array([float(p2[s][f'PC{i}']) for s in common])
        rows.append({'metric': f'PC{i}_abs_pearson_r', 'value': round(abs(float(np.corrcoef(x, y)[0, 1])), 6)})
    n1, k1 = read_k(a.main_kin); n2, k2 = read_k(a.sens_kin)
    order = [n2.index(n) for n in n1]
    k2 = k2[np.ix_(order, order)]
    iu = np.triu_indices(len(n1), 1)
    rows.append({'metric': 'K_offdiag_pearson_r', 'value': round(float(np.corrcoef(k1[iu], k2[iu])[0, 1]), 6)})
    write_tsv(a.out, rows, ['metric', 'value'])
    for r in rows:
        log(f"    {r['metric']:<32} {r['value']}")


if __name__ == '__main__':
    main()
