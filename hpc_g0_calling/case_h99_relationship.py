#!/usr/bin/env python3
"""case_h99_relationship.py -- X1 side check of the case-report isolates (Cn76, Cn84) and the H99 controls. Standard library only.

Input: the X1 VCF (7 samples joint-genotyped and filtered like base 1, in a callset of their own). Output (--out-dir): pairwise_distance.tsv (SNP differences inside the region callable in both samples), pairwise_diff_sites.tsv, site_patterns.tsv and h99_allele_carriers.tsv (exact-allele carrier patterns), sequencing_metadata.tsv (run / flowcell / lane from the FASTQ headers) and summary.txt.
Observations only: no inference on within-host evolution, reinfection, mixed infection or contamination.
"""
import argparse
import bisect
import gzip
import itertools
import os
import re
import sys

H99_A = "CnH99"
H99_B = ("CnH9904", "CnH9905", "CnH9906", "CnH9907")


def open_text(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def read_bed(path):
    iv = {}
    with open_text(path) as fh:
        for line in fh:
            if line.strip():
                c, a, b = line.split("\t")[:3]
                iv.setdefault(c, []).append((int(a), int(b)))
    for c in iv:
        iv[c].sort()
        merged = []
        for a, b in iv[c]:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        iv[c] = merged
    return iv


def intersect(x, y):
    out = {}
    for c in set(x) & set(y):
        i = j = 0
        a, b = x[c], y[c]
        res = []
        while i < len(a) and j < len(b):
            lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
            if lo < hi:
                res.append((lo, hi))
            if a[i][1] < b[j][1]:
                i += 1
            else:
                j += 1
        out[c] = res
    return out


def total_bp(iv):
    return sum(b - a for v in iv.values() for a, b in v)


def nearest_distance(iv, c, pos):
    v = iv.get(c, [])
    if not v:
        return ""
    starts = [a for a, _ in v]
    i = bisect.bisect_right(starts, pos - 1)
    best = None
    for k in (i - 1, i):
        if 0 <= k < len(v):
            a, b = v[k]
            d = 0 if a <= pos - 1 < b else min(abs(pos - 1 - a), abs(pos - 1 - (b - 1)))
            best = d if best is None else min(best, d)
    return str(best)


def parse_header(h):
    """Instrument, run, flowcell, lane, index from a FASTQ header (Illumina or MGI style)."""
    h = h.strip().lstrip("@")
    first, _, rest = h.partition(" ")
    p = first.split(":")
    if len(p) >= 7:
        idx = rest.split(":")[-1] if rest.count(":") >= 3 else ""
        return "illumina", p[0], p[1], p[2], p[3], idx
    m = re.match(r"^([A-Z]*\d+)L(\d+)C\d+R\d+", first)
    if m:
        return "mgi", "", "", m.group(1), m.group(2), ""
    return "unknown", "", "", "", "", ""


def in_iv(starts, ends, pos):
    """1-based pos inside sorted, merged 0-based half-open intervals."""
    i = bisect.bisect_right(starts, pos - 1) - 1
    return i >= 0 and pos - 1 < ends[i]


def trim(ref, alt):
    """Trim the common suffix so the same SNP written with a longer REF compares equal."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return ref, alt


def main():
    ap = argparse.ArgumentParser(description="X1: Cn76 / Cn84 / H99 relationship on one set of reliable sites")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--callable-dir", required=True)
    ap.add_argument("--analyzable", required=True)
    ap.add_argument("--mask", required=True)
    ap.add_argument("--main-raw", required=True, help="05 main raw joint VCF (unfiltered)")
    ap.add_argument("--samples-tsv", required=True)
    ap.add_argument("--raw-root", required=True)
    ap.add_argument("--st5-case", default="Cn84")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    # read X1 VCF
    names = []
    sites = []   # (chrom, pos, ref, alt, qual, info, g); g[sample] = (allele or None, AD, SRC, ADV)
    with open_text(a.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                names = f[9:]
                continue
            alle = [f[3]] + f[4].split(",")
            fmt = f[8].split(":")
            idx = {k: fmt.index(k) for k in ("AD", "SRC", "ADV") if k in fmt}
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            g = {}
            for n, v in zip(names, f[9:]):
                q = v.split(":")
                gt = q[0]
                get = lambda k: q[idx[k]] if k in idx and idx[k] < len(q) else "."
                g[n] = (None if gt in (".", "./.", ".|.") else alle[int(gt)], get("AD"), get("SRC"), get("ADV"))
            sites.append((f[0], int(f[1]), f[3], f[4], f[5], info, g))
    h99 = [n for n in names if n == H99_A or n in H99_B]
    print(f"case_h99_relationship.py: {len(sites)} sites, samples {names}", file=sys.stderr)

    # pairwise distances inside the region callable in both samples
    analyzable = read_bed(a.analyzable)
    cal = {n: intersect(read_bed(os.path.join(a.callable_dir, f"{n}.cna3.bed.gz")), analyzable) for n in names}
    analyzable_bp = total_bp(analyzable)
    pair_rows, diff_rows = [], []
    for x, y in itertools.combinations(names, 2):
        shared = intersect(cal[x], cal[y])
        idx_iv = {c: ([s for s, _ in v], [e for _, e in v]) for c, v in shared.items()}
        both_all = diff_all = both_in = diff_in = 0
        for c, pos, ref, alt, q, info, g in sites:
            ax, ay = g[x][0], g[y][0]
            if ax is None or ay is None:
                continue
            inside = c in idx_iv and in_iv(idx_iv[c][0], idx_iv[c][1], pos)
            both_all += 1
            both_in += inside
            if ax != ay:
                diff_all += 1
                diff_in += inside
                if a.st5_case not in (x, y):
                    diff_rows.append([x, y, c, str(pos), ref, alt, ax, ay, g[x][1], g[y][1], g[x][2], g[y][2], g[x][3], g[y][3], "1" if inside else "0"])
        shared_bp = total_bp(shared)
        d_norm = diff_in / shared_bp * analyzable_bp if shared_bp else None
        pair_rows.append([x, y, str(shared_bp), str(both_in), str(diff_in), f"{d_norm:.1f}" if d_norm is not None else "NA", str(both_all), str(diff_all)])
    with open(os.path.join(a.out_dir, "pairwise_distance.tsv"), "w") as fh:
        fh.write("sample_a\tsample_b\tshared_callable_bp\tsites_both_called_in_shared\tD_raw_in_shared\tD_per_analyzable_genome\tsites_both_called_all\tD_raw_all_called\n")
        for r in pair_rows:
            fh.write("\t".join(r) + "\n")
    with open(os.path.join(a.out_dir, "pairwise_diff_sites.tsv"), "w") as fh:
        fh.write("sample_a\tsample_b\tchrom\tpos\tref\talt\tallele_a\tallele_b\tad_a\tad_b\tsrc_a\tsrc_b\tadv_a\tadv_b\tin_shared_callable\n")
        for r in diff_rows:
            fh.write("\t".join(r) + "\n")

    # exact-allele carrier patterns
    patterns = {}
    for c, pos, ref, alt, q, info, g in sites:
        carried = []
        for al in alt.split(","):
            who = [n for n in names if g[n][0] == al]
            if who:
                carried.append(f"{al}:{','.join(who)}")
        missing = ",".join(n for n in names if g[n][0] is None)
        key = (";".join(carried) or "-", missing or "-")
        patterns[key] = patterns.get(key, 0) + 1
    with open(os.path.join(a.out_dir, "site_patterns.tsv"), "w") as fh:
        fh.write("allele_carriers\tN_samples\tsites\n")
        for (carriers, missing), k in sorted(patterns.items(), key=lambda kv: -kv[1]):
            fh.write(f"{carriers}\t{missing}\t{k}\n")

    # alleles carried by any H99 stock or by Cn76
    def state(g, n, ref, allele):
        al = g[n][0]
        if al is None:
            return "N"
        if al == allele:
            return "same"
        if al == ref:
            return "REF"
        return f"other:{al}"
    rows = []
    for c, pos, ref, alt, q, info, g in sites:
        for al in alt.split(","):
            h_carriers = [n for n in h99 if g[n][0] == al]
            cn76 = "Cn76" in names and g["Cn76"][0] == al
            if not h_carriers and not cn76:
                continue
            rows.append((c, pos, ref, al, q, info, g, h_carriers))
    want = {(c, pos) for c, pos, *_ in rows}
    main_raw = {}
    with open_text(a.main_raw) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 10)
            key = (f[0], int(f[1]))
            if key not in want:
                continue
            f = line.rstrip("\n").split("\t")
            alle = [f[3]] + f[4].split(",")
            gts = [v.split(":")[0] for v in f[9:]]
            main_raw.setdefault(key, []).append((f[3], alle, gts))
    mask = read_bed(a.mask)
    out_rows = []
    for c, pos, ref, al, q, info, g, h_carriers in rows:
        recs = main_raw.get((c, pos))
        if not recs:
            m_rec, m_car, m_called = "no_record", "NA", "NA"
        else:
            car = called = 0
            for mref, alle, gts in recs:
                same = {i for i, x in enumerate(alle) if i > 0 and trim(mref, x) == trim(ref, al)}
                for gt in gts:
                    if gt in (".", "./.", ".|."):
                        continue
                    called += 1
                    if gt.isdigit() and int(gt) in same:
                        car += 1
            m_rec, m_car, m_called = "present", str(car), str(called)
        h_n = [n for n in h99 if g[n][0] is None]
        out_rows.append([c, str(pos), ref, al, q, info.get("MQ", "."), info.get("QD", "."), info.get("FS", "."), nearest_distance(mask, c, pos), ",".join(h_carriers) or "-", ",".join(h_n) or "-"]
                        + [f"{state(g, n, ref, al)}({g[n][1]};{g[n][2]};{g[n][3]})" for n in names] + [m_rec, m_car, m_called])
    with open(os.path.join(a.out_dir, "h99_allele_carriers.tsv"), "w") as fh:
        fh.write("chrom\tpos\tref\tallele\tqual\tMQ\tQD\tFS\tmask_distance_bp\tH99_carriers_same_allele\tH99_N\t" + "\t".join(names) + "\tmain_raw_record\tmain_raw_carriers_same_allele\tmain_raw_called\n")
        for r in out_rows:
            fh.write("\t".join(r) + "\n")

    # sequencing metadata from FASTQ headers
    seq_rows = []
    with open(a.samples_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            path = os.path.join(a.raw_root, r["fastq_R1"])
            try:
                with gzip.open(path, "rt") as gz:
                    hdr = gz.readline().strip()
            except OSError as e:
                hdr = f"unreadable: {e}"
            kind, inst, run, cell, lane, idx = parse_header(hdr)
            seq_rows.append([r["sample"], r["role"], kind, inst, run, cell, lane, idx, hdr])
    with open(os.path.join(a.out_dir, "sequencing_metadata.tsv"), "w") as fh:
        fh.write("sample\trole\tstyle\tinstrument\trun\tflowcell\tlane\tindex\tfirst_R1_header\n")
        for r in seq_rows:
            fh.write("\t".join(r) + "\n")

    # summary
    L = ["X1 side check: Cn76, Cn84 and the H99 controls on one set of reliable sites", "=" * 72]
    L.append(f"sites (filtered SNPs of the 7-sample joint call, static mask removed): {len(sites)}")
    L.append("Pairwise SNP distance inside the region callable in both samples (base 2 of each, intersected with the analyzable region); D_norm = D_raw_in_shared / shared callable bp x analyzable bp. Last two columns: all sites called in both, for reference.")
    for r in pair_rows:
        L.append(f"  {r[0]:<8} {r[1]:<8} shared_callable_bp={r[2]:<10} both_called_in_shared={r[3]:<7} D_raw_in_shared={r[4]:<7} D_norm={r[5]:<8} (all called: {r[6]} sites, D_raw {r[7]})")
    L.append("Alleles carried by at least one H99 stock or by Cn76 (exact allele), by H99 carrier set and Cn76 state; main = carriers of the same allele in the main raw callset (no_record = the main raw callset has no record at that position):")
    tally = {}
    for r in out_rows:
        cn76 = r[11 + names.index("Cn76")].split("(")[0] if "Cn76" in names else "NA"
        main = "no_record" if r[-3] == "no_record" else ("0" if r[-2] == "0" else ">=1")
        key = (r[9], cn76, main)
        tally[key] = tally.get(key, 0) + 1
    for (hc, cn76, main), k in sorted(tally.items(), key=lambda kv: (-kv[1], kv[0])):
        L.append(f"  H99={hc:<40} Cn76={cn76:<10} main={main:<9} alleles={k}")
    if not tally:
        L.append("  none")
    L.append("Sequencing metadata of the 7 samples (first R1 header):")
    for r in seq_rows:
        if r[0] in names:
            L.append(f"  {r[0]:<8} style={r[2]} instrument={r[3]} run={r[4]} flowcell={r[5]} lane={r[6]} index={r[7]}")
    for n in names:
        me = next((r for r in seq_rows if r[0] == n), None)
        if me and me[5]:
            same = [r[0] for r in seq_rows if r[5] == me[5] and r[6] == me[6] and r[0] != n]
            L.append(f"  same flowcell+lane as {n}: {len(same)} samples" + (f" ({', '.join(same[:15])}{' ...' if len(same) > 15 else ''})" if same else ""))
    L.append("Notes: carrier patterns are reported as observed; no inference on ancestral state, reference-sequence errors or laboratory origin is drawn here.")
    with open(os.path.join(a.out_dir, "summary.txt"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
