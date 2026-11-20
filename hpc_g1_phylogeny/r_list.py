#!/usr/bin/env python3
"""r_list.py -- R list: positions where >= --min-strains strains carry an intermediate-allele-fraction genotype.

Input --new-geno: per-genotype listing (CHROM POS SAMPLE) made by 01_prep with bcftools from the main base 1 (FORMAT/ADV = M, sum(AD) >= 10, max(AD)/sum(AD) in [0.2, 0.8]).
Outputs (--out-dir): R.tsv, R.bed (one base per line, lexicographic order) and iaf_positions_by_strain_count.tsv.
"""
import argparse
import collections

from g1lib import die, log, write_tsv


def read_geno(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) == 3 and f[0]:
                rows.append((f[0], int(f[1]), f[2]))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--new-geno', required=True)
    ap.add_argument('--min-strains', type=int, required=True)
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()

    new = read_geno(a.new_geno)
    if len(set(new)) != len(new):
        die('duplicate (CHROM, POS, SAMPLE) rows in a genotype listing')

    new_site = collections.Counter((c, p) for c, p, _ in new)
    r_sites = {k for k, n in new_site.items() if n >= a.min_strains}

    # Distribution of the number of strains per position (descriptive).
    bins = [(1, 1), (2, 9), (10, 95), (96, 192)]
    dist = []
    for lo, hi in bins:
        sites = [k for k, n in new_site.items() if lo <= n <= hi]
        dist.append({'n_strains': f'{lo}-{hi}' if lo != hi else str(lo), 'positions': len(sites),
                     'genotypes': sum(new_site[k] for k in sites)})
    write_tsv(f'{a.out_dir}/iaf_positions_by_strain_count.tsv', dist, ['n_strains', 'positions', 'genotypes'])

    rows = sorted(((c, p, new_site[(c, p)]) for c, p in r_sites), key=lambda x: (x[0], x[1]))
    with open(f'{a.out_dir}/R.tsv', 'w') as fh:
        fh.write('chrom\tpos\tn_strains\n')
        for c, p, n in rows:
            fh.write(f'{c}\t{p}\t{n}\n')
    with open(f'{a.out_dir}/R.bed', 'w') as fh:   # one base per line, lexicographic order (bedtools)
        for c, p, _ in rows:
            fh.write(f'{c}\t{p - 1}\t{p}\n')
    log(f'    R = {len(rows)} positions with >= {a.min_strains} strains (of {len(new_site)} positions with >= 1)')


if __name__ == '__main__':
    main()
