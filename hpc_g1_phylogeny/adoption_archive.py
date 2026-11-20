#!/usr/bin/env python3
"""adoption_archive.py -- cluster adoption table. Standard library only.

Candidates: every distinct member set that forms a cluster of >= 2 representatives at some SNP-address threshold (D_raw, main analysis). Per candidate it reports whether the set is one side of a supported split of the L1 tree and a cluster at >= 2 adjacent ladder thresholds (crit3), plus diameter, chain flag, bridge strains, sampling years, HIV status and usable rate.
The adoption decision (column user_decision) is left empty. Outputs: adoption_candidates.tsv and adoption_reference.tsv.
"""
import argparse
import collections
import datetime
import statistics

from g1lib import read_samples, read_tsv, write_tsv, log


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--clusters', required=True, help='snp_address_clusters.tsv (main, D_raw)')
    ap.add_argument('--supported', required=True, help='supported_clades.tsv (04_tree)')
    ap.add_argument('--sample-table', required=True)
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--usable-stats', required=True)
    ap.add_argument('--analyzable-bp', type=int, required=True)
    ap.add_argument('--ladder', required=True)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    ladder = [int(x) for x in a.ladder.split()]

    reps = [r['sample'] for r in read_samples(a.sample_table) if r['role'] == 'representative']
    allreps = frozenset(reps)
    meta = {r['sample']: r for r in read_tsv(a.metadata)}
    usable = {r['sample']: int(r['usable_bp']) / a.analyzable_bp for r in read_tsv(a.usable_stats)}
    sup = set()
    for r in read_tsv(a.supported):
        s = frozenset(r['tips_below'].split(','))
        sup.add(s)
        sup.add(allreps - s)

    by_set = collections.defaultdict(dict)
    for r in read_tsv(a.clusters):
        by_set[frozenset(r['members'].split(','))][int(r['threshold'])] = r
    rows = []
    for members, lv in sorted(by_set.items(), key=lambda kv: (-len(kv[0]), sorted(kv[0]))):
        ts = sorted(lv)
        adjacent = any(ladder.index(t2) - ladder.index(t1) == 1 for t1, t2 in zip(ts, ts[1:]))
        in_tree = members in sup
        years = sorted(datetime.date.fromisoformat(meta[s]['sampling_date']).year for s in members)
        hiv = collections.Counter(meta[s]['hiv_status'] for s in members)
        rows.append({
            'n_members': len(members), 'thresholds': ','.join(map(str, ts)),
            'crit3_supported_split': 'yes' if in_tree else 'no', 'crit3_adjacent_thresholds': 'yes' if adjacent else 'no',
            'crit3': 'yes' if in_tree and adjacent else 'no',
            'diameter_by_threshold': ';'.join(f"{t}:{lv[t]['diameter']}" for t in ts),
            'diameter_gt_2t_any': 'yes' if any(lv[t]['diameter_gt_2t'] == 'yes' for t in ts) else 'no',
            'bridge_strains_by_threshold': ';'.join(f"{t}:{lv[t]['bridge_strains'] or '-'}" for t in ts),
            'year_min': years[0], 'year_median': statistics.median(years), 'year_max': years[-1],
            'hiv_positive': hiv.get('positive', 0), 'hiv_negative': hiv.get('negative', 0),
            'usable_rate_median': round(statistics.median(usable[s] for s in members), 6),
            'members': ','.join(sorted(members, key=reps.index)), 'user_decision': ''})
    cols = ['n_members', 'thresholds', 'crit3_supported_split', 'crit3_adjacent_thresholds', 'crit3', 'diameter_by_threshold',
            'diameter_gt_2t_any', 'bridge_strains_by_threshold', 'year_min', 'year_median', 'year_max', 'hiv_positive', 'hiv_negative',
            'usable_rate_median', 'members', 'user_decision']
    write_tsv(f'{a.out_dir}/adoption_candidates.tsv', rows, cols)
    ref = [{'reference': '160 representatives', 'year_min': min(datetime.date.fromisoformat(meta[s]['sampling_date']).year for s in reps),
            'year_max': max(datetime.date.fromisoformat(meta[s]['sampling_date']).year for s in reps),
            'hiv_positive': sum(meta[s]['hiv_status'] == 'positive' for s in reps), 'hiv_negative': sum(meta[s]['hiv_status'] == 'negative' for s in reps),
            'usable_rate_median': round(statistics.median(usable[s] for s in reps), 6)}]
    write_tsv(f'{a.out_dir}/adoption_reference.tsv', ref, list(ref[0].keys()))

    n3 = sum(r['crit3'] == 'yes' for r in rows)
    log(f'    adoption candidates: {len(rows)} member sets; criterion 3 met by {n3}')


if __name__ == '__main__':
    main()
