#!/usr/bin/env python3
"""rulers_address.py -- rulers A / B / C, threshold anchors and SNP address of one analysis. Needs numpy, networkx.

Input: pairwise_L0.tsv of the analysis (l0_distances.py assemble); only gate-PASS pairs enter summaries, anchors and SNP-address edges.
Ruler A: H99 replicate pairs; ruler B: within-patient pairs (sampling dates, same / different admission and specimen, per-isolate intermediate-fraction flags); ruler C: representative pairs.
Anchors: t_tech (95th percentile of ruler A) and t_event (max D_raw of same-day different-site pairs, rounded up to a multiple of 5). SNP address: single-linkage clusters of the representatives at each ladder threshold, with diameter and bridge strains, on D_raw and D_norm.
"""
import argparse
import collections
import datetime
import math

import networkx as nx
import numpy as np

from g1lib import die, log, read_samples, read_tsv, write_tsv, pair_key

QS = (0, 5, 25, 50, 75, 95, 100)


def qdesc(vals, prefix=''):
    if not vals:
        return {f'{prefix}n': 0}
    v = np.asarray(vals, dtype=float)
    d = {f'{prefix}n': len(vals)}
    for q in QS:
        d[f'{prefix}p{q}'] = round(float(np.percentile(v, q, method='linear')), 4)
    return d


def interval_class(days):
    if days == 0:
        return '0'
    if days <= 7:
        return '1-7'
    if days <= 30:
        return '8-30'
    return '>30'


def snp_address(reps, pairs, metric, ladder, out_dir, suffix):
    """pairs: {(a,b): (value, gate_pass)} for representative pairs."""
    order = {s: i for i, s in enumerate(reps)}
    level_rows, cluster_rows = [], []
    labels = {s: [] for s in reps}
    for t in ladder:
        g = nx.Graph()
        g.add_nodes_from(reps)
        for (x, y), (v, ok) in pairs.items():
            if ok and v <= t:
                g.add_edge(x, y)
        comps = sorted((sorted(c, key=order.get) for c in nx.connected_components(g)), key=lambda c: order[c[0]])
        sizes = collections.Counter(len(c) for c in comps)
        n_chain = n_bridge = 0
        for cid, members in enumerate(comps, start=1):
            for s in members:
                labels[s].append(str(cid))
            if len(members) < 2:
                continue
            inside = [pairs[pair_key(x, y)] for i, x in enumerate(members) for y in members[i + 1:]]
            diam = max((v for v, ok in inside if ok), default=0)
            n_fail_inside = sum(1 for _, ok in inside if not ok)
            bridges = sorted(nx.articulation_points(g.subgraph(members)), key=order.get)
            chain = diam > 2 * t
            n_chain += chain
            n_bridge += len(bridges)
            cluster_rows.append({'threshold': t, 'cluster_id': cid, 'size': len(members), 'diameter': diam,
                                 'diameter_gt_2t': 'yes' if chain else 'no', 'gate_fail_pairs_inside': n_fail_inside,
                                 'bridge_strains': ','.join(bridges), 'members': ','.join(members)})
        level_rows.append({'threshold': t, 'n_clusters': len(comps), 'n_singletons': sizes.get(1, 0),
                           'n_multi': len(comps) - sizes.get(1, 0), 'max_size': max(sizes),
                           'size_distribution': ';'.join(f'{k}:{sizes[k]}' for k in sorted(sizes)),
                           'n_diameter_gt_2t': n_chain, 'n_bridge_strains': n_bridge})
    addr_rows = []
    for s in reps:   # address coarse -> fine, e.g. t100.t50.t25.t10.t5
        ids = labels[s]
        addr_rows.append({'sample': s, **{f't{t}': ids[i] for i, t in enumerate(ladder)},
                          'address': '.'.join(ids[i] for i in reversed(range(len(ladder))))})
    write_tsv(f'{out_dir}/snp_address_levels{suffix}.tsv', level_rows, list(level_rows[0].keys()))
    write_tsv(f'{out_dir}/snp_address_clusters{suffix}.tsv', cluster_rows,
              ['threshold', 'cluster_id', 'size', 'diameter', 'diameter_gt_2t', 'gate_fail_pairs_inside', 'bridge_strains', 'members'])
    write_tsv(f'{out_dir}/snp_address{suffix}.tsv', addr_rows, ['sample'] + [f't{t}' for t in ladder] + ['address'])
    return level_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--pairs', required=True)
    ap.add_argument('--sample-table', required=True)
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--qc-gate', required=True)
    ap.add_argument('--usable-stats', required=True, help='stats of the main pseudo-genome set (sample, usable_bp)')
    ap.add_argument('--ladder', required=True)
    ap.add_argument('--t-tech-percentile', type=float, required=True)
    ap.add_argument('--t-event-round', type=int, required=True)
    ap.add_argument('--low-usable-frac', type=float, required=True)
    ap.add_argument('--label', required=True)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    ladder = [int(x) for x in a.ladder.split()]

    samples = read_samples(a.sample_table)
    order = {r['sample']: i for i, r in enumerate(samples)}
    reps = [r['sample'] for r in samples if r['role'] == 'representative']
    meta = {r['sample']: r for r in read_tsv(a.metadata)}
    gate = {r['sample']: r for r in read_tsv(a.qc_gate)}
    usable = {r['sample']: int(r['usable_bp']) for r in read_tsv(a.usable_stats)}
    low_thr = a.low_usable_frac * float(np.median([usable[s] for s in reps]))
    rows = read_tsv(a.pairs)
    for r in rows:
        r['D_raw'] = int(r['D_raw']); r['N_shared'] = int(r['N_shared']); r['D_norm'] = float(r['D_norm']); r['ok'] = r['gate'] == 'PASS'
    by_type = collections.defaultdict(list)
    for r in rows:
        by_type[r['pair_type']].append(r)

    # Ruler A + Cn169 / Cn169b
    a_rows = [{'sample_a': r['sample_a'], 'sample_b': r['sample_b'], 'D_raw': r['D_raw'], 'N_shared': r['N_shared'],
               'D_norm': r['D_norm'], 'gate': r['gate']} for r in by_type['h99_replicate'] + by_type['culture_replicate']]
    write_tsv(f'{a.out_dir}/ruler_A_pairs.tsv', a_rows, ['sample_a', 'sample_b', 'D_raw', 'N_shared', 'D_norm', 'gate'])
    h_ok = [r['D_raw'] for r in by_type['h99_replicate'] if r['ok']]
    if len(by_type['h99_replicate']) != 10:
        die(f'expected 10 H99 pairs, got {len(by_type["h99_replicate"])}')
    t_tech = float(np.percentile(np.asarray(h_ok, dtype=float), a.t_tech_percentile, method='linear')) if h_ok else float('nan')

    # Ruler B
    b_rows = []
    for r in by_type['within_patient']:
        x, y = sorted((r['sample_a'], r['sample_b']), key=lambda s: (meta[s]['sampling_date'], order[s]))
        dx = datetime.date.fromisoformat(meta[x]['sampling_date'])
        dy = datetime.date.fromisoformat(meta[y]['sampling_date'])
        days = (dy - dx).days
        b_rows.append({
            'patient_id': meta[x]['patient_id'], 'isolate_1': x, 'isolate_2': y,
            'sampling_date_1': meta[x]['sampling_date'], 'sampling_date_2': meta[y]['sampling_date'], 'delta_days': days,
            'interval_class': interval_class(days),
            'admission_1': meta[x]['admission_no'], 'admission_2': meta[y]['admission_no'],
            'admission_axis': 'same' if meta[x]['admission_no'] == meta[y]['admission_no'] else 'different',
            'specimen_1': meta[x]['specimen'], 'specimen_2': meta[y]['specimen'],
            'site_axis': 'same' if meta[x]['specimen'] == meta[y]['specimen'] else 'different',
            'D_raw': r['D_raw'], 'N_shared': r['N_shared'], 'D_norm': r['D_norm'], 'gate': r['gate'],
            'iaf_frac_1': gate[x]['intermediate_af_frac'], 'iaf_den_1': gate[x]['n_sites_ad10_measured'],
            'iaf_frac_2': gate[y]['intermediate_af_frac'], 'iaf_den_2': gate[y]['n_sites_ad10_measured'],
            'low_usable_1': 'yes' if usable[x] < low_thr else 'no', 'low_usable_2': 'yes' if usable[y] < low_thr else 'no',
            '_ok': r['ok']})
    b_rows.sort(key=lambda r: (r['patient_id'], r['sampling_date_1'], r['isolate_1'], r['isolate_2']))
    b_cols = [k for k in b_rows[0] if not k.startswith('_')]
    write_tsv(f'{a.out_dir}/ruler_B_pairs.tsv', b_rows, b_cols)
    cat_rows = []
    for axis in ('admission_axis', 'site_axis', 'interval_class'):
        groups = collections.defaultdict(list)
        for r in b_rows:
            if r['_ok']:
                groups[r[axis]].append(r)
        for val in sorted(groups):
            g = groups[val]
            cat_rows.append({'axis': axis, 'category': val, 'n_patients': len({r['patient_id'] for r in g}),
                             **qdesc([r['D_raw'] for r in g], 'D_raw_'), **{k: v for k, v in qdesc([r['D_norm'] for r in g], 'D_norm_').items() if k != 'D_norm_n'}})
    cat_rows.append({'axis': 'all', 'category': 'all', 'n_patients': len({r['patient_id'] for r in b_rows if r['_ok']}),
                     **qdesc([r['D_raw'] for r in b_rows if r['_ok']], 'D_raw_'),
                     **{k: v for k, v in qdesc([r['D_norm'] for r in b_rows if r['_ok']], 'D_norm_').items() if k != 'D_norm_n'}})
    write_tsv(f'{a.out_dir}/ruler_B_by_category.tsv', cat_rows, list(cat_rows[-1].keys()))
    pat = collections.defaultdict(list)
    for r in b_rows:
        if r['_ok']:
            pat[r['patient_id']].append(r)
    p_rows = []
    for p in sorted(pat):
        g = pat[p]
        strains = {r['isolate_1'] for r in g} | {r['isolate_2'] for r in g}
        p_rows.append({'patient_id': p, 'n_strains': len(strains), 'n_pairs': len(g),
                       'admission_axis': ','.join(sorted({r['admission_axis'] for r in g})),
                       'D_raw_median': float(np.median([r['D_raw'] for r in g])), 'D_raw_max': max(r['D_raw'] for r in g),
                       'D_norm_median': round(float(np.median([r['D_norm'] for r in g])), 4), 'D_norm_max': max(r['D_norm'] for r in g)})
    write_tsv(f'{a.out_dir}/ruler_B_by_patient.tsv', p_rows, list(p_rows[0].keys()))
    pat_summary = [{'statistic': 'per-patient median D_raw', **qdesc([r['D_raw_median'] for r in p_rows])},
                   {'statistic': 'per-patient max D_raw', **qdesc([r['D_raw_max'] for r in p_rows])}]
    write_tsv(f'{a.out_dir}/ruler_B_patient_summary.tsv', pat_summary, list(pat_summary[0].keys()))

    # Ruler C
    c_ok = [r for r in by_type['rep_rep'] if r['ok']]
    c_sum = [{'metric': 'D_raw', **qdesc([r['D_raw'] for r in c_ok])}, {'metric': 'D_norm', **qdesc([r['D_norm'] for r in c_ok])},
             {'metric': 'gate_fail_pairs', 'n': len(by_type['rep_rep']) - len(c_ok)}]
    write_tsv(f'{a.out_dir}/ruler_C_summary.tsv', c_sum, list(c_sum[0].keys()))
    hist = collections.Counter(r['D_raw'] for r in c_ok)
    write_tsv(f'{a.out_dir}/ruler_C_hist.tsv', [{'D_raw': k, 'n_pairs': hist[k]} for k in sorted(hist)], ['D_raw', 'n_pairs'])

    # Anchors
    ev = [r['D_raw'] for r in b_rows if r['_ok'] and r['delta_days'] == 0 and r['site_axis'] == 'different']
    t_event = a.t_event_round * math.ceil(max(ev) / a.t_event_round) if ev else None
    b_vals = [r['D_raw'] for r in b_rows if r['_ok']]
    c_vals = [r['D_raw'] for r in c_ok]
    lo, hi = max(min(b_vals), min(c_vals)), min(max(b_vals), max(c_vals))
    anchors = [
        {'anchor': 't_tech', 'value': round(t_tech, 4), 'note': f'p{a.t_tech_percentile:g} (linear) of {len(h_ok)} H99 pair D_raw'},
        {'anchor': 't_event', 'value': '' if t_event is None else t_event,
         'note': f'max D_raw {max(ev) if ev else "NA"} over {len(ev)} same-day different-site pairs, rounded up to {a.t_event_round}'},
        {'anchor': 'B_C_overlap', 'value': f'{lo}-{hi}' if lo <= hi else 'none', 'note': 'overlap of the ruler B and ruler C D_raw ranges'},
        {'anchor': 'C_pairs_le_maxB', 'value': round(sum(v <= max(b_vals) for v in c_vals) / len(c_vals), 6), 'note': 'fraction of ruler C pairs with D_raw <= max ruler B'},
        {'anchor': 'B_pairs_ge_minC', 'value': round(sum(v >= min(c_vals) for v in b_vals) / len(b_vals), 6), 'note': 'fraction of ruler B pairs with D_raw >= min ruler C'},
    ]
    for name, vals in (('ruler_B', b_vals), ('ruler_C', c_vals)):
        for k, v in qdesc(vals).items():
            anchors.append({'anchor': f'{name}_{k}', 'value': v, 'note': 'D_raw, gate-PASS pairs'})
    write_tsv(f'{a.out_dir}/anchors.tsv', anchors, ['anchor', 'value', 'note'])

    # SNP address (D_raw main, D_norm sensitivity)
    rp = {pair_key(r['sample_a'], r['sample_b']): r for r in by_type['rep_rep']}
    lv = snp_address(reps, {k: (r['D_raw'], r['ok']) for k, r in rp.items()}, 'D_raw', ladder, a.out_dir, '')
    snp_address(reps, {k: (r['D_norm'], r['ok']) for k, r in rp.items()}, 'D_norm', ladder, a.out_dir, '_dnorm')

    log(f'    [{a.label}] ruler A: H99 D_raw {sorted(h_ok)} -> t_tech {t_tech:.2f}; Cn169-Cn169b '
        f'{[r["D_raw"] for r in by_type["culture_replicate"]]}')
    log(f'    [{a.label}] ruler B: {len(b_rows)} pairs ({len(b_vals)} gate PASS), D_raw median {np.median(b_vals):.1f}, max {max(b_vals)}; t_event {t_event}')
    log(f'    [{a.label}] ruler C: {len(c_vals)} pairs, D_raw median {np.median(c_vals):.1f}, range {min(c_vals)}-{max(c_vals)}')
    for r in lv:
        log(f'    [{a.label}] SNP address t={r["threshold"]}: {r["n_clusters"]} clusters ({r["n_multi"]} with >= 2 strains, max size {r["max_size"]})')


if __name__ == '__main__':
    main()
