#!/usr/bin/env python3
"""l0_distances.py -- L0 pairwise distances. Needs numpy.

extract   alignment of the given positions (VCF order) cut out of each sample's pseudo-genome -> FASTA (input of snp-dists)
assemble  pairwise table of one analysis: D_raw (snp-dists), N_shared = |U_i ∩ U_j| (usable bases, through the complements M_i = A' - U_i and one bedtools multiinter), N_median, D_norm = D_raw / N_shared x N_median and the gate N_shared >= gate_frac x N_median
"""
import argparse
import collections
import itertools
import statistics

import numpy as np

from g1lib import Faidx, die, log, read_names, read_samples, read_tsv, write_tsv, pair_key



# extract

def cmd_extract(a):
    positions = []
    with open(a.positions) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 2 and f[0]:
                positions.append((f[0], int(f[1])))
    by_contig = collections.OrderedDict()
    for k, (c, p) in enumerate(positions):
        by_contig.setdefault(c, []).append((k, p - 1))
    names = read_names(a.samples)
    with open(a.out, 'wb') as out:
        for s in names:
            fa = Faidx(f'{a.fasta_dir}/{s}.fa')
            col = bytearray(len(positions))
            for c, items in by_contig.items():
                seq = fa.fetch(c)
                for k, p0 in items:
                    col[k] = seq[p0]
            out.write(b'>' + s.encode() + b'\n' + bytes(col) + b'\n')
    log(f'    extract: {len(names)} sequences x {len(positions)} positions -> {a.out}')


# helpers

def parse_kv(items):
    out = {}
    for it in items:
        k, v = it.split('=', 1)
        out[k] = v
    return out


def read_molten(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) != 3 or f[0] == f[1]:
                continue
            d[pair_key(f[0], f[1])] = int(f[2])
    return d


def bed_bp(path):
    s = 0
    with open(path) as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) >= 3:
                s += int(f[2]) - int(f[1])
    return s


def shared_matrix(names, multi_path, aprime_bp):
    """N_shared matrix from the multiinter output of the M_i (complements of the usable bases inside A')."""
    idx = {n: i for i, n in enumerate(names)}
    m = np.zeros((len(names), len(names)), dtype=np.int64)
    with open(multi_path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 5 or f[0] == 'chrom':
                continue
            ids = [idx[x] for x in f[4].split(',')]
            ln = int(f[2]) - int(f[1])
            ix = np.array(ids)
            m[np.ix_(ix, ix)] += ln
    diag = np.diag(m).copy()
    return aprime_bp - diag[:, None] - diag[None, :] + m, diag


# assemble

def cmd_assemble(a):
    samples = read_samples(a.sample_table)
    meta = {r['sample']: r for r in read_tsv(a.metadata)}
    role = {r['sample']: r['role'] for r in samples}
    callset = {r['sample']: r['callset'] for r in samples}
    order = {r['sample']: i for i, r in enumerate(samples)}
    sets, dists, multi, aprime, stats = (parse_kv(x) for x in (a.set_samples, a.dists, a.multi, a.aprime, a.stats))

    shared, d_raw = {}, {}
    for cs in sets:
        names = read_names(sets[cs])
        ap_bp = bed_bp(aprime[cs])
        mat, diag = shared_matrix(names, multi[cs], ap_bp)
        usable = {r['sample']: int(r['usable_bp']) for r in read_tsv(stats[cs])}
        bad = [n for i, n in enumerate(names) if ap_bp - diag[i] != usable[n]]
        if bad:
            die(f'{cs}: |A\'| - |M_i| differs from the usable bp of the pseudo-genome for {len(bad)} samples (e.g. {bad[:3]})')
        for i, j in itertools.combinations(range(len(names)), 2):
            shared[pair_key(names[i], names[j])] = int(mat[i, j])
        d_raw.update(read_molten(dists[cs]))
        log(f'    {cs}: |A\'| = {ap_bp} bp, {len(names)} samples, usable bp min {min(usable.values())} max {max(usable.values())}')

    # Pair list
    reps = [r['sample'] for r in samples if r['role'] == 'representative']
    pairs = [(x, y, 'rep_rep') for x, y in itertools.combinations(reps, 2)]
    by_patient = collections.defaultdict(list)
    for r in samples:
        if r['role'] in ('representative', 'within_patient_extra'):
            by_patient[meta[r['sample']]['patient_id']].append(r['sample'])
    multi_pat = {p: v for p, v in by_patient.items() if len(v) > 1}
    wp = [(x, y, 'within_patient') for v in multi_pat.values() for x, y in itertools.combinations(sorted(v, key=order.get), 2)]
    pairs += wp
    for r in samples:
        if r['role'] == 'culture_replicate_of_Cn169':
            pairs.append(('Cn169', r['sample'], 'culture_replicate'))
    h99 = [r['sample'] for r in samples if r['role'] == 'technical_replicate']
    pairs += [(x, y, 'h99_replicate') for x, y in itertools.combinations(h99, 2)]

    rep_shared = [shared[pair_key(x, y)] for x, y, t in pairs if t == 'rep_rep']
    n_median = statistics.median(rep_shared)
    gate_min = a.gate_frac * n_median
    rows = []
    for x, y, t in pairs:
        k = pair_key(x, y)
        if k not in d_raw or k not in shared:
            die(f'pair {x}-{y} missing from the distance or shared tables (callsets {callset[x]}/{callset[y]})')
        ns, dr = shared[k], d_raw[k]
        rows.append({'sample_a': x, 'sample_b': y, 'callset': callset[x], 'pair_type': t,
                     'patient_id': meta[x]['patient_id'] if t == 'within_patient' else '',
                     'D_raw': dr, 'N_shared': ns, 'D_norm': f'{dr / ns * n_median:.4f}',
                     'gate': 'PASS' if ns >= gate_min else 'FAIL'})
    cols = ['sample_a', 'sample_b', 'callset', 'pair_type', 'patient_id', 'D_raw', 'N_shared', 'D_norm', 'gate']
    write_tsv(f'{a.out_dir}/pairwise_L0.tsv', rows, cols)
    n_fail = collections.Counter(r['pair_type'] for r in rows if r['gate'] == 'FAIL')
    with open(f'{a.out_dir}/n_median.tsv', 'w') as fh:
        fh.write(f'metric\tvalue\nN_median\t{n_median}\ngate_min_N_shared\t{gate_min:.1f}\nrep_pairs\t{len(rep_shared)}\n')
        for t in ('rep_rep', 'within_patient', 'culture_replicate', 'h99_replicate'):
            fh.write(f'gate_fail_{t}\t{n_fail.get(t, 0)}\n')
    log(f'    N_median = {n_median}; gate failures by pair type: {dict(n_fail) or "none"}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    e = sub.add_parser('extract')
    e.add_argument('--fasta-dir', required=True)
    e.add_argument('--samples', required=True)
    e.add_argument('--positions', required=True)
    e.add_argument('--out', required=True)
    s = sub.add_parser('assemble')
    s.add_argument('--sample-table', required=True)
    s.add_argument('--metadata', required=True)
    for k in ('set-samples', 'dists', 'multi', 'aprime', 'stats'):
        s.add_argument(f'--{k}', nargs='+', required=True, help='callset=path')
    s.add_argument('--gate-frac', type=float, required=True)
    s.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    {'extract': cmd_extract, 'assemble': cmd_assemble}[a.cmd](a)


if __name__ == '__main__':
    main()
