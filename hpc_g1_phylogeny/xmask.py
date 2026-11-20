#!/usr/bin/env python3
"""xmask.py -- site-exclusion region X_c of one callset, masked for every sample of the callset.

X_c = union of (a) every SNP record of the site-filtered VCF whose FILTER is not PASS and (b) every non-SNP record of the callset, each spanning [POS, POS + len(REF) - 1] (a symbolic allele with INFO/END spans [POS, END]).
Non-SNP records come from --nonsnp (a VCF holding only non-SNP records) or from --raw minus --snp. Output: BED (0-based, merged, lexicographic contig order) and a summary TSV.
"""
import argparse
import re

from g1lib import die, log, merge, open_text, write_tsv

END_RE = re.compile(r'(?:^|;)END=(\d+)')


def records(path):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t', 8)
            yield f


def span(f, stats):
    pos, ref, alts = int(f[1]), f[3], f[4].split(',')
    end = pos + len(ref) - 1
    if any(a.startswith('<') for a in alts):
        m = END_RE.search(f[7])
        stats['symbolic'] += 1
        if m:
            end = max(end, int(m.group(1)))
            stats['symbolic_with_end'] += 1
    return f[0], pos - 1, end   # 0-based half-open


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--name', required=True)
    ap.add_argument('--sitefilt', required=True)
    ap.add_argument('--nonsnp')
    ap.add_argument('--raw')
    ap.add_argument('--snp')
    ap.add_argument('--out-bed', required=True)
    ap.add_argument('--out-summary', required=True)
    a = ap.parse_args()
    if bool(a.nonsnp) == bool(a.raw and a.snp):
        die('give either --nonsnp or both --raw and --snp')

    stats = {'symbolic': 0, 'symbolic_with_end': 0}
    ivs = {}
    n_fail = n_pass = 0
    for f in records(a.sitefilt):
        if f[6] == 'PASS':
            n_pass += 1
            continue
        n_fail += 1
        c, s, e = span(f, stats)
        ivs.setdefault(c, []).append((s, e))

    n_nonsnp = 0
    if a.nonsnp:
        src = records(a.nonsnp)
        snp_keys = None
    else:
        snp_keys = {(f[0], f[1], f[3], f[4]) for f in records(a.snp)}
        src = records(a.raw)
    n_star_only = 0
    for f in src:
        if snp_keys is not None and (f[0], f[1], f[3], f[4]) in snp_keys:
            continue
        n_nonsnp += 1
        if all(x == '*' for x in f[4].split(',')):
            n_star_only += 1
        c, s, e = span(f, stats)
        ivs.setdefault(c, []).append((s, e))

    out = []
    for c in sorted(ivs):
        for s, e in merge(ivs[c]):
            out.append((c, s, e))
    with open(a.out_bed, 'w') as fh:
        for c, s, e in out:
            fh.write(f'{c}\t{s}\t{e}\n')
    bp = sum(e - s for _, s, e in out)
    summary = [
        {'metric': 'callset', 'value': a.name},
        {'metric': 'sitefilt_pass_records', 'value': n_pass},
        {'metric': 'sitefilt_failed_records', 'value': n_fail},
        {'metric': 'nonsnp_records', 'value': n_nonsnp},
        {'metric': 'nonsnp_star_only_records', 'value': n_star_only},
        {'metric': 'symbolic_allele_records', 'value': stats['symbolic']},
        {'metric': 'symbolic_with_END', 'value': stats['symbolic_with_end']},
        {'metric': 'x_intervals', 'value': len(out)},
        {'metric': 'x_bp', 'value': bp},
    ]
    write_tsv(a.out_summary, summary, ['metric', 'value'])
    log(f'    X_{a.name}: {n_fail} failed SNP records + {n_nonsnp} non-SNP records -> {len(out)} intervals, {bp} bp'
        f' (symbolic alleles: {stats["symbolic"]}, with END: {stats["symbolic_with_end"]})')


if __name__ == '__main__':
    main()
