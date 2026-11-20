#!/usr/bin/env python3
"""asr_branches.py -- per-branch counts from the IQ-TREE ancestral reconstruction and tree statistics. Needs numpy, scipy, dendropy.

strip   copy of a tree without internal node labels (the ancestral run gets a clean tree)
count   inputs: .state (internal nodes), the --ancestral tree file, the L1 alignment and its columns table, the ML tree with support labels. Node state = argmax of p_A..p_T (a tie = undetermined; tip N = undetermined); per branch n_sub (both ends determined and different), n_unres and n_sub_lowpost; zero-SNP branches, polytomies, topology parameters and supported clades.
Outputs: branches.tsv, topology_parameters.tsv, supported_clades.tsv, ml_splits.tsv
"""
import argparse
import collections
import statistics

import dendropy
import numpy as np
from scipy import stats

from g1lib import die, log, read_fasta_alignment, read_tsv, write_tsv

BASE_IDX = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3}


def load_tree(path):
    return dendropy.Tree.get(path=path, schema='newick', preserve_underscores=True)


def cmd_strip(a):
    t = load_tree(a.tree)
    for nd in t.preorder_node_iter():
        if not nd.is_leaf():
            nd.label = None
    t.write(path=a.out, schema='newick', suppress_rooting=True, unquoted_underscores=True)


def leaf_sets(tree):
    """{node: frozenset(leaf labels below it)} (postorder)."""
    below = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            below[nd] = frozenset([nd.taxon.label])
        else:
            below[nd] = frozenset().union(*(below[c] for c in nd.child_nodes()))
    return below


def split_key(s, all_leaves, ref):
    """Canonical side of a bipartition: the side without the reference taxon."""
    return frozenset(all_leaves - s) if ref in s else s


def parse_support(label):
    if not label:
        return None, None
    f = label.split('/')
    try:
        if len(f) >= 2:
            return float(f[-2]), float(f[-1])
        return None, float(f[0])
    except ValueError:
        return None, None


def cmd_count(a):
    names, seqs = read_fasta_alignment(a.aln)
    cls = np.array([r['class'] for r in read_tsv(a.columns)])
    ncols = len(cls)
    tips = {}
    for n in names:
        arr = np.frombuffer(seqs[n], dtype=np.uint8)
        if len(arr) != ncols:
            die(f'{n}: alignment length {len(arr)} != columns table {ncols}')
        st = np.full(ncols, -1, dtype=np.int8)
        for code, k in BASE_IDX.items():
            st[arr == code] = k
        tips[n] = st

    tree = load_tree(a.asr_tree)
    internal = [nd for nd in tree.preorder_node_iter() if not nd.is_leaf()]
    lab = {nd: nd.label for nd in internal}
    if any(not l for l in lab.values()):
        die('the --ancestral tree file has unlabeled internal nodes')
    idx = {l: i for i, l in enumerate(lab.values())}
    istate = np.full((len(internal), ncols), -2, dtype=np.int8)     # -2 = not seen yet
    imaxp = np.zeros((len(internal), ncols), dtype=np.float32)
    n_rows = 0
    with open(a.state) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('Node\t') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            i, j = idx.get(f[0]), int(f[1]) - 1
            if i is None:
                die(f'.state node {f[0]} not in the tree')
            if not 0 <= j < ncols:
                die(f'.state site {j + 1} outside the alignment ({ncols} columns)')
            p = [float(x) for x in f[3:7]]
            mx = max(p)
            istate[i, j] = -1 if p.count(mx) > 1 else p.index(mx)
            imaxp[i, j] = mx
            n_rows += 1
    if (istate == -2).any():
        die('.state does not cover every internal node x column')
    log(f'    .state: {n_rows} rows, {len(internal)} internal nodes x {ncols} columns; tied maxima: {int((istate == -1).sum())}')

    def node_state(nd):
        return tips[nd.taxon.label] if nd.is_leaf() else istate[idx[lab[nd]]]

    def node_maxp(nd):
        return np.ones(ncols, dtype=np.float32) if nd.is_leaf() else imaxp[idx[lab[nd]]]

    # ML supports mapped by bipartition
    all_leaves = frozenset(names)
    ref = names[0]
    ml = load_tree(a.ml_tree)
    ml_below = leaf_sets(ml)
    support = {}
    ml_rows = []
    for nd in ml.preorder_node_iter():
        if nd.parent_node is None or nd.is_leaf():
            continue
        alrt, ufb = parse_support(nd.label)
        key = split_key(ml_below[nd], all_leaves, ref)
        support[key] = (alrt, ufb)
        ml_rows.append({'n_side': len(key), 'sh_alrt': alrt, 'ufboot': ufb, 'side': ','.join(sorted(key))})
    write_tsv(f'{a.out_dir}/ml_splits.tsv', ml_rows, ['n_side', 'sh_alrt', 'ufboot', 'side'])

    below = leaf_sets(tree)
    edges = []
    for nd in tree.preorder_node_iter():
        par = nd.parent_node
        if par is None:
            continue
        sp, sc = node_state(par), node_state(nd)
        det = (sp >= 0) & (sc >= 0)
        sub = det & (sp != sc)
        low = sub & (np.minimum(node_maxp(par), node_maxp(nd)) < a.lowpost)
        key = split_key(below[nd], all_leaves, ref)
        alrt, ufb = (None, None) if nd.is_leaf() else support.get(key, (None, None))
        if not nd.is_leaf() and key not in support:
            die('an internal branch of the --ancestral tree has no matching split in the ML tree (topologies differ)')
        edges.append({'parent': lab[par], 'child': nd.taxon.label if nd.is_leaf() else lab[nd],
                      'terminal': 'yes' if nd.is_leaf() else 'no', 'n_tips_below': len(below[nd]),
                      'n_sub': int(sub.sum()), 'n_unres': int((~det).sum()), 'n_sub_lowpost': int(low.sum()),
                      'sh_alrt': '' if alrt is None else alrt, 'ufboot': '' if ufb is None else ufb,
                      'tips_below': ','.join(sorted(below[nd])) if not nd.is_leaf() else '',
                      '_par': par, '_nd': nd, '_key': key})
    write_tsv(f'{a.out_dir}/branches.tsv', edges,
              ['parent', 'child', 'terminal', 'n_tips_below', 'n_sub', 'n_unres', 'n_sub_lowpost', 'sh_alrt', 'ufboot', 'tips_below'])

    iii = cls == 'iii'

    # zero-SNP branches, collapse, polytomies
    inner = [e for e in edges if e['terminal'] == 'no']
    term = [e for e in edges if e['terminal'] == 'yes']
    zero_det = [e for e in inner if e['n_sub'] == 0 and e['n_unres'] == 0]
    zero_unres = [e for e in inner if e['n_sub'] == 0 and e['n_unres'] > 0]

    def polytomies(collapse):
        parent = {}   # union-find over node objects; a node absent from the dict is its own root
        def find(x):
            while x in parent:
                x = parent[x]
            return x
        cset = {id(e['_nd']) for e in collapse}
        for e in collapse:
            ra, rb = find(e['_par']), find(e['_nd'])
            if ra is not rb:
                parent[rb] = ra
        deg = collections.Counter()
        for e in edges:
            if id(e['_nd']) in cset:
                continue
            for x in (e['_par'], e['_nd']):
                if not x.is_leaf():
                    deg[find(x)] += 1
        kids = [d - 1 for d in deg.values() if d >= 4]
        return len(kids), kids

    n_poly, kids = polytomies(zero_det)
    n_poly_s, kids_s = polytomies(zero_det + zero_unres)
    ins, tes = [e['n_sub'] for e in inner], [e['n_sub'] for e in term]
    ks = stats.ks_2samp(ins, tes)
    gate = [r for r in read_tsv(a.qc_gate) if r['role'] == 'representative']
    d_ref = [float(r['d_raw_to_ref']) for r in gate]
    med = statistics.median(d_ref)
    mad = statistics.median(abs(x - med) for x in d_ref)
    params = [
        ('variable_sites_for_tree (class iii)', int(iii.sum())),
        ('alignment_columns (class ii + iii)', ncols),
        ('internal_branches', len(inner)), ('terminal_branches', len(term)),
        ('internal_branch_n_sub_median', statistics.median(ins)), ('terminal_branch_n_sub_median', statistics.median(tes)),
        ('ks_internal_vs_terminal_statistic', round(float(ks.statistic), 6)), ('ks_internal_vs_terminal_p', f'{ks.pvalue:.3g}'),
        ('internal_zero_determined (collapsed)', len(zero_det)), ('internal_zero_with_unresolved', len(zero_unres)),
        ('collapsed_fraction_of_all_branches', round(len(zero_det) / len(edges), 6)),
        ('polytomies', n_poly), ('polytomy_children_median', statistics.median(kids) if kids else ''),
        ('polytomy_children_max', max(kids) if kids else ''),
        ('polytomies_if_unresolved_zero_also_collapsed', n_poly_s),
        ('terminal_branch_n_sub_total', sum(tes)),
        ('substitutions_with_maxpost_lt_0.95', sum(e['n_sub_lowpost'] for e in edges)),
        ('distance_to_H99_median (G0 d_raw_to_ref, 160 representatives)', med),
        ('distance_to_H99_MAD', mad), ('distance_to_H99_SD', round(statistics.stdev(d_ref), 2)),
    ]
    write_tsv(f'{a.out_dir}/topology_parameters.tsv', [{'parameter': k, 'value': v} for k, v in params], ['parameter', 'value'])

    sup = []
    for e in inner:
        if e['sh_alrt'] != '' and e['ufboot'] != '' and e['sh_alrt'] >= a.alrt_min and e['ufboot'] >= a.ufboot_min:
            sup.append({'n_tips_below': e['n_tips_below'], 'sh_alrt': e['sh_alrt'], 'ufboot': e['ufboot'], 'n_sub': e['n_sub'],
                        'n_unres': e['n_unres'], 'tips_below': e['tips_below']})
    write_tsv(f'{a.out_dir}/supported_clades.tsv', sup, ['n_tips_below', 'sh_alrt', 'ufboot', 'n_sub', 'n_unres', 'tips_below'])
    log(f'    branches: {len(inner)} internal (median n_sub {statistics.median(ins)}), {len(term)} terminal (median {statistics.median(tes)}); '
        f'collapsed {len(zero_det)}, zero-with-unresolved {len(zero_unres)}; polytomies {n_poly} ({n_poly_s} if both collapsed); '
        f'supported clades {len(sup)}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('strip')
    s.add_argument('--tree', required=True)
    s.add_argument('--out', required=True)
    c = sub.add_parser('count')
    c.add_argument('--state', required=True)
    c.add_argument('--asr-tree', required=True)
    c.add_argument('--ml-tree', required=True)
    c.add_argument('--aln', required=True)
    c.add_argument('--columns', required=True)
    c.add_argument('--qc-gate', required=True)
    c.add_argument('--lowpost', type=float, required=True)
    c.add_argument('--alrt-min', type=float, required=True)
    c.add_argument('--ufboot-min', type=float, required=True)
    c.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    {'strip': cmd_strip, 'count': cmd_count}[a.cmd](a)


if __name__ == '__main__':
    main()
