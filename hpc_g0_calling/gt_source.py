#!/usr/bin/env python3
"""gt_source.py -- where did one sample's genotype at each base-1 site come from? Standard library only.

For every site of a callset (VCF order) it reads the sample's own gVCF and reports the source (V variant record, B reference block, D spanned by an upstream deletion, - not covered), the covering record (start, end, MIN_DP, DP, GQ, REF, ALT, AD) and whether the base is callable in base 2.
Output: gzipped TSV with one line per site in the order of --sites, read by gt_mask.py in lockstep with the VCF.
"""
import argparse
import bisect
import gzip
import sys
from collections import deque


def open_text(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def load_callable(path):
    starts, ends = {}, {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            c, a, b = line.split("\t")[:3]
            starts.setdefault(c, []).append(int(a))
            ends.setdefault(c, []).append(int(b))
    for c in starts:
        order = sorted(range(len(starts[c])), key=lambda i: starts[c][i])
        starts[c] = [starts[c][i] for i in order]
        ends[c] = [ends[c][i] for i in order]
    return starts, ends


def is_callable(starts, ends, chrom, pos):
    s = starts.get(chrom)
    if not s:
        return 0
    i = bisect.bisect_right(s, pos - 1) - 1          # BED is 0-based half-open, pos is 1-based
    return 1 if i >= 0 and pos - 1 < ends[chrom][i] else 0


def gvcf_records(path, rank):
    """Yield (rank, start, end, is_block, min_dp, dp, gq, ref, alt, ad) for each gVCF record."""
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] not in rank:
                continue
            start = int(f[1])
            end = start + len(f[3]) - 1
            for kv in f[7].split(";"):
                if kv.startswith("END="):
                    end = int(kv[4:])
                    break
            is_block = f[4] == "<NON_REF>"
            fmt = f[8].split(":")
            vals = f[9].split(":")
            get = lambda k: vals[fmt.index(k)] if k in fmt and fmt.index(k) < len(vals) else "."
            yield rank[f[0]], start, end, is_block, get("MIN_DP"), get("DP"), get("GQ"), f[3], f[4], get("AD")


def main():
    ap = argparse.ArgumentParser(description="Genotype source (variant record / reference block) per site for one sample")
    ap.add_argument("--gvcf", required=True)
    ap.add_argument("--sites", required=True, help="chrom<TAB>pos, VCF order (.fai order)")
    ap.add_argument("--callable", required=True, help="base 2 BED(.gz) of the same sample")
    ap.add_argument("--genome", required=True, help="genome.txt (contig<TAB>length, .fai order)")
    ap.add_argument("--out", required=True, help="output .tsv.gz")
    a = ap.parse_args()

    rank = {}
    with open(a.genome) as fh:
        for i, line in enumerate(fh):
            if line.strip():
                rank[line.split("\t")[0]] = i
    sites = []
    with open(a.sites) as fh:
        for line in fh:
            if line.strip():
                c, p = line.split("\t")[:2]
                sites.append((rank[c], int(p), c))
    for k in range(1, len(sites)):
        if sites[k][:2] < sites[k - 1][:2]:
            sys.exit(f"gt_source.py: sites are not in .fai order at {sites[k][2]}:{sites[k][1]}")
    cs, ce = load_callable(a.callable)

    recs = gvcf_records(a.gvcf, rank)
    active = deque()          # records that may still cover upcoming sites
    pending = next(recs, None)
    counts = {"V": 0, "B": 0, "D": 0, "-": 0}
    with gzip.open(a.out, "wt") as out:
        out.write("chrom\tpos\tsrc\trec_start\trec_end\tmin_dp\trec_dp\trec_gq\tcallable\trec_ref\trec_alt\trec_ad\n")
        for r, p, c in sites:
            while pending is not None and (pending[0], pending[1]) <= (r, p):
                active.append(pending)
                pending = next(recs, None)
            active = deque(x for x in active if (x[0], x[2]) >= (r, p))
            starting_var = [x for x in active if x[0] == r and x[1] == p and not x[3]]
            blocks = [x for x in active if x[0] == r and x[3] and x[1] <= p <= x[2]]
            spanning = [x for x in active if x[0] == r and not x[3] and x[1] < p <= x[2]]
            if starting_var:
                src, rec = "V", starting_var[0]
            elif blocks:
                src, rec = "B", blocks[0]
            elif spanning:
                src, rec = "D", spanning[0]
            else:
                src, rec = "-", None
            counts[src] += 1
            if rec is None:
                fields, tail = [".", ".", ".", ".", "."], [".", ".", "."]
            else:
                fields, tail = [str(rec[1]), str(rec[2]), rec[4], rec[5], rec[6]], [rec[7], rec[8], rec[9]]
            out.write("\t".join([c, str(p), src] + fields + [str(is_callable(cs, ce, c, p))] + tail) + "\n")
    print(f"gt_source.py: {len(sites)} sites  V={counts['V']} B={counts['B']} D={counts['D']} uncovered={counts['-']}", file=sys.stderr)


if __name__ == "__main__":
    main()
