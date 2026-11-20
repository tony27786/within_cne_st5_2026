#!/usr/bin/env python3
"""qc_gate_table.py -- strain-level QC gate table, one row per sample. Standard library only.

Gate items: callable fraction of the analyzable region (PASS if >= 0.90 x the median of the reference role), median depth (PASS if >= 30), distance to the reference (robust z of D_norm, PASS if |z| <= 3.5, ST5 only) and the MLST group from samples.tsv.
The intermediate allele fraction (base-1 sites with major-allele AD fraction in [0.2, 0.8] over sites with confirmed AD and sum(AD) >= 10) is reported for every sample without any threshold. Nothing is excluded here; the table only records the gate.
"""
import argparse
import gzip
import math
import os
import statistics
import sys
from fractions import Fraction

NA = "NA"


def read_tsv(path):
    """Return (header, rows as dicts) of a TSV with a header line."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return header, rows


def num(x):
    if x is None or x in ("", NA, "."):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def fmt(x, nd=6):
    if x is None:
        return NA
    if isinstance(x, bool):
        return "PASS" if x else "FAIL"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        if math.isinf(x):
            return "inf" if x > 0 else "-inf"
        return format(x, f".{nd}g")
    return str(x)


def median(vals):
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else None


def parse_flagstat(path):
    """pct_mapped, pct_dup over primary reads (falls back to totals)."""
    if not os.path.exists(path):
        return None, None
    vals = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 4 or parts[1] != "+":
                continue
            try:
                n = int(parts[0]) + int(parts[2])
            except ValueError:
                continue
            label = " ".join(parts[3:]).split(" (")[0]
            vals.setdefault(label, n)
    tot = vals.get("primary") or vals.get("in total")
    mapped = vals.get("primary mapped", vals.get("mapped"))
    dup = vals.get("primary duplicates", vals.get("duplicates"))
    if not tot:
        return None, None
    return (100.0 * mapped / tot if mapped is not None else None, 100.0 * dup / tot if dup is not None else None)


def scan_base1(path, iaf_lo, iaf_hi, min_ad, detail=None):
    """Per-sample counts from a base-1 VCF (GT already masked, AD intact). The intermediate allele fraction uses only AD confirmed against the sample gVCF (FORMAT/ADV = M), whatever the GT. If detail is a list, every intermediate-fraction genotype is appended to it."""
    lo_n, lo_d = iaf_lo.numerator, iaf_lo.denominator
    hi_n, hi_d = iaf_hi.numerator, iaf_hi.denominator
    counts = {}
    samples = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                counts = {s: {"n_sites": 0, "n_called": 0, "n_ref": 0, "n_alt": 0, "n_N": 0, "n_ad10": 0, "n_iaf": 0, "n_ad10_synth": 0, "n_ad10_unverified": 0, "n_no_adv": 0} for s in samples}
                continue
            f = line.rstrip("\n").split("\t")
            fmt_keys = f[8].split(":")
            i_ad = fmt_keys.index("AD") if "AD" in fmt_keys else -1
            i_adv = fmt_keys.index("ADV") if "ADV" in fmt_keys else -1
            i_src = fmt_keys.index("SRC") if "SRC" in fmt_keys else -1
            for j, s in enumerate(samples):
                c = counts[s]
                c["n_sites"] += 1
                v = f[9 + j].split(":")
                gt = v[0]
                if gt in (".", "./.", ".|."):
                    c["n_N"] += 1
                else:
                    c["n_called"] += 1
                    if gt == "0":
                        c["n_ref"] += 1
                    else:
                        c["n_alt"] += 1
                if i_adv < 0:
                    c["n_no_adv"] += 1          # base 1 without the AD evidence check: not evaluable
                    continue
                adv = v[i_adv] if i_adv < len(v) else "."
                ad = None
                if 0 <= i_ad < len(v):
                    parts = v[i_ad].split(",")
                    if all(p not in ("", ".") for p in parts):
                        ad = [int(p) for p in parts]
                if ad is None or sum(ad) < min_ad:
                    continue
                tot = sum(ad)
                if adv == "S":
                    c["n_ad10_synth"] += 1
                    continue
                if adv != "M":
                    c["n_ad10_unverified"] += 1
                    continue
                c["n_ad10"] += 1
                m = max(ad)
                if m * lo_d >= lo_n * tot and m * hi_d <= hi_n * tot:
                    c["n_iaf"] += 1
                    if detail is not None:
                        src = v[i_src] if 0 <= i_src < len(v) else "."
                        detail.append([s, f[0], f[1], f[3], f[4], gt, v[i_ad], str(tot), f"{m / tot:.4f}", src, adv, f[5]])
    return counts


def main():
    ap = argparse.ArgumentParser(description="G0 strain-level QC gate")
    ap.add_argument("--samples", required=True, help="sample, group, role, callset TSV")
    ap.add_argument("--analyzable-bp", type=int, required=True)
    ap.add_argument("--callable-dir", required=True, help="03_callable/summary")
    ap.add_argument("--union-callable", required=True, help="sample, callable_bp_analyzable_union")
    ap.add_argument("--flagstat-dir", required=True)
    ap.add_argument("--contigs", required=True, help="chroms.list (.fai order)")
    ap.add_argument("--mito", required=True)
    ap.add_argument("--base1", action="append", required=True, help="callset=path (repeat)")
    ap.add_argument("--gtmask", action="append", default=[], help="callset=path (repeat)")
    ap.add_argument("--mask-summary", default=None)
    ap.add_argument("--ref-role", required=True)
    ap.add_argument("--callable-ratio-min", required=True)
    ap.add_argument("--depth-median-min", type=float, required=True)
    ap.add_argument("--iaf-lo", required=True, help="lower bound of the intermediate major-allele fraction, e.g. 0.2")
    ap.add_argument("--iaf-hi", required=True, help="upper bound, e.g. 0.8")
    ap.add_argument("--iaf-min-ad", type=int, required=True, help="denominator: sum(AD) >= this")
    ap.add_argument("--refdist-z-max", type=float, required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-depth", required=True)
    ap.add_argument("--out-iaf-sites", required=True, help="every intermediate-fraction genotype of the h99 and case callsets, site by site")
    ap.add_argument("--detail-callsets", default="h99,case")
    a = ap.parse_args()

    ratio_min = float(Fraction(a.callable_ratio_min))
    iaf_lo, iaf_hi = Fraction(a.iaf_lo), Fraction(a.iaf_hi)

    _, manifest = read_tsv(a.samples)
    names = [r["sample"] for r in manifest]
    warnings = []

    # base 2 summaries (stage 03)
    callable_info = {}
    depth_chrom = {}
    for s in names:
        p = os.path.join(a.callable_dir, f"{s}.callable.tsv")
        if os.path.exists(p):
            _, rows = read_tsv(p)
            callable_info[s] = rows[0] if rows else {}
        else:
            callable_info[s] = {}
            warnings.append(f"{s}: no callable summary (stage 03)")
        q = os.path.join(a.callable_dir, f"{s}.depth_by_chrom.tsv")
        depth_chrom[s] = {}
        if os.path.exists(q):
            _, rows = read_tsv(q)
            depth_chrom[s] = {r["chrom"]: r["median_depth"] for r in rows}
    _, urows = read_tsv(a.union_callable)
    union_bp = {r["sample"]: num(r.get("callable_bp_analyzable_union")) for r in urows}

    # base 1 scans (stage 06), one per callset
    b1 = {}
    detail_cs = set(x for x in a.detail_callsets.split(",") if x)
    iaf_detail = []
    for item in a.base1:
        cs, path = item.split("=", 1)
        if os.path.exists(path):
            det = [] if cs in detail_cs else None
            counts = scan_base1(path, iaf_lo, iaf_hi, a.iaf_min_ad, det)
            for s, c in counts.items():
                b1[s] = (cs, c)
            n_no_adv = sum(c["n_no_adv"] for c in counts.values())
            if n_no_adv:
                warnings.append(f"base 1 of callset {cs} has {n_no_adv} genotypes without FORMAT/ADV (written by an older 06?): intermediate allele fraction set to NA for its samples")
            if det:
                iaf_detail.extend([cs] + row for row in det)
        else:
            warnings.append(f"base 1 missing for callset {cs}: {path}")
    gtm = {}
    for item in a.gtmask:
        cs, path = item.split("=", 1)
        if os.path.exists(path):
            _, rows = read_tsv(path)
            for r in rows:
                gtm[r["sample"]] = r

    # per-sample raw metrics
    rec = {}
    for r in manifest:
        s = r["sample"]
        ci = callable_info.get(s, {})
        cb = num(ci.get("callable_bp_analyzable"))
        c = b1.get(s, (None, None))[1]
        if r["callset"] != b1.get(s, (r["callset"], None))[0]:
            warnings.append(f"{s}: found in base 1 of callset {b1[s][0]}, manifest says {r['callset']}")
        pm, pd = parse_flagstat(os.path.join(a.flagstat_dir, f"{s}.flagstat.txt"))
        rec[s] = {
            "callable_bp": cb,
            "callable_frac": (cb / a.analyzable_bp) if cb is not None else None,
            "callable_bp_genome": num(ci.get("callable_bp_genome")),
            "callable_frac_union": None,
            "median_depth": num(ci.get("median_depth_analyzable")),
            "mean_depth": num(ci.get("mean_depth_analyzable")),
            "pct_mapped": pm,
            "pct_dup": pd,
            "c": c,
        }
        ub = union_bp.get(s)
        if ub is not None:
            rec[s]["callable_bp_union"] = ub
        if c is None:
            warnings.append(f"{s}: not found in any base-1 VCF")

    # analyzable bp after the union mask = denominator of callable_frac_union
    union_analyzable = None
    if a.mask_summary and os.path.exists(a.mask_summary):
        _, ms = read_tsv(a.mask_summary)
        for row in ms:
            if row.get("metric") == "analyzable_union_bp":
                union_analyzable = num(row.get("value"))
    for s in names:
        ub = rec[s].get("callable_bp_union")
        if ub is not None and union_analyzable:
            rec[s]["callable_frac_union"] = ub / union_analyzable

    # reference distributions (reference role)
    reps = [r["sample"] for r in manifest if r["role"] == a.ref_role]
    med_callable_frac = median([rec[s]["callable_frac"] for s in reps])
    med_callable_bp = median([rec[s]["callable_bp"] for s in reps])
    for s in names:
        c = rec[s]["c"]
        d_raw = c["n_alt"] if c else None
        rec[s]["d_raw"] = d_raw
        cb = rec[s]["callable_bp"]
        rec[s]["d_norm"] = (d_raw / cb * med_callable_bp) if (d_raw is not None and cb and med_callable_bp) else None
    rep_dnorm = [rec[s]["d_norm"] for s in reps if rec[s]["d_norm"] is not None]
    med_dnorm = median(rep_dnorm)
    mad_dnorm = median([abs(x - med_dnorm) for x in rep_dnorm]) if rep_dnorm else None
    rep_draw = [rec[s]["d_raw"] for s in reps if rec[s]["d_raw"] is not None]
    med_draw = median(rep_draw)
    mad_draw = median([abs(x - med_draw) for x in rep_draw]) if rep_draw else None

    # gate
    group_of = {r["sample"]: r["group"] for r in manifest}
    out_cols = [
        "sample", "group", "role", "callset",
        "callable_bp_analyzable", "callable_frac", "callable_ratio_to_median", "pass_callable",
        "median_depth", "pass_depth",
        "n_sites_ad10_measured", "n_intermediate_af_sites", "intermediate_af_frac", "n_ad10_excluded_synthesized", "n_ad10_excluded_unverified",
        "d_raw_to_ref", "d_raw_minus_median", "d_norm_to_ref", "d_norm_minus_median", "refdist_robust_z", "pass_refdist",
        "mlst_st_stringmlst", "pass_mlst",
        "l1_discovery_eligible", "gate_overall", "gate_flags", "gate_actions",
        "callable_bp_genome", "callable_frac_union_mask", "mean_depth",
        "pct_mapped_primary", "pct_dup_primary",
        "n_base1_sites", "n_called", "n_ref", "n_alt", "n_N",
        "n_mask_gq", "n_mask_dp", "n_mask_ad", "n_mask_star", "n_ruleA_kept",
    ]
    rows_out = []
    tally = {"callable": 0, "depth": 0, "refdist": 0, "mlst": 0, "na": 0}
    for r in manifest:
        s = r["sample"]
        x = rec[s]
        c = x["c"]
        flags, actions = [], []

        ratio = (x["callable_frac"] / med_callable_frac) if (x["callable_frac"] is not None and med_callable_frac) else None
        p_callable = (ratio >= ratio_min) if ratio is not None else None
        p_depth = (x["median_depth"] >= a.depth_median_min) if x["median_depth"] is not None else None

        # Intermediate allele fraction: reported only (no threshold, no flag); NA when the denominator is 0 or base 1 lacks FORMAT/ADV.
        iaf_ok = bool(c) and c["n_no_adv"] == 0
        n_ad10 = c["n_ad10"] if iaf_ok else None
        n_iaf = c["n_iaf"] if iaf_ok else None
        iaf_frac = (n_iaf / n_ad10) if (iaf_ok and n_ad10) else None

        z = None
        p_ref = None
        if group_of[s] == "ST5" and x["d_norm"] is not None and med_dnorm is not None:
            if mad_dnorm:
                z = (x["d_norm"] - med_dnorm) / (1.4826 * mad_dnorm)
            else:
                z = 0.0 if x["d_norm"] == med_dnorm else math.copysign(math.inf, x["d_norm"] - med_dnorm)
            p_ref = abs(z) <= a.refdist_z_max

        # MLST = pre-inclusion stringMLST (samples.tsv group). Not applicable to the H99 controls, nor to Cn76, which is ST2 by design (case report).
        mlst = group_of[s]
        if r["role"] == "technical_replicate":
            p_mlst = None
        elif r["role"] == "case_report_ST2_to_ST5":
            p_mlst = True if mlst == "ST5" else None
        else:
            p_mlst = (mlst == "ST5")

        for item, ok in (("callable", p_callable), ("depth", p_depth), ("refdist", p_ref), ("mlst", p_mlst)):
            if ok is False:
                flags.append(f"{item}_fail")
                tally[item] += 1
        missing = []
        if p_callable is None:
            missing.append("callable")
        if p_depth is None:
            missing.append("depth")
        if group_of[s] == "ST5" and p_ref is None:
            missing.append("refdist")
        if missing:
            flags.append("not_evaluable:" + ",".join(missing))
            tally["na"] += 1

        if p_callable is False or p_depth is False:
            actions.append("exclude_from_L1_site_discovery(distances_kept_flagged)")
        if p_ref is False:
            actions.append("check_species_or_contamination")
        if p_mlst is False:
            actions.append("exclude_and_record")

        if r["role"] == a.ref_role:
            l1 = "YES" if (p_callable and p_depth) else "NO"
        else:
            l1 = NA
        overall = "PASS" if not flags else "FLAG"
        g = gtm.get(s, {})

        rows_out.append({
            "sample": s,
            "group": r["group"],
            "role": r["role"],
            "callset": r["callset"],
            "callable_bp_analyzable": fmt(int(x["callable_bp"]) if x["callable_bp"] is not None else None),
            "callable_frac": fmt(x["callable_frac"]),
            "callable_ratio_to_median": fmt(ratio),
            "pass_callable": fmt(p_callable),
            "median_depth": fmt(x["median_depth"]),
            "pass_depth": fmt(p_depth),
            "n_sites_ad10_measured": fmt(n_ad10),
            "n_intermediate_af_sites": fmt(n_iaf),
            "intermediate_af_frac": fmt(iaf_frac),
            "n_ad10_excluded_synthesized": fmt(c["n_ad10_synth"] if iaf_ok else None),
            "n_ad10_excluded_unverified": fmt(c["n_ad10_unverified"] if iaf_ok else None),
            "d_raw_to_ref": fmt(x["d_raw"]),
            "d_raw_minus_median": fmt((x["d_raw"] - med_draw) if (x["d_raw"] is not None and med_draw is not None and group_of[s] == "ST5") else None),
            "d_norm_to_ref": fmt(x["d_norm"]),
            "d_norm_minus_median": fmt((x["d_norm"] - med_dnorm) if (x["d_norm"] is not None and med_dnorm is not None and group_of[s] == "ST5") else None),
            "refdist_robust_z": fmt(z, 4),
            "pass_refdist": fmt(p_ref),
            "mlst_st_stringmlst": mlst,
            "pass_mlst": fmt(p_mlst),
            "l1_discovery_eligible": l1,
            "gate_overall": overall,
            "gate_flags": ";".join(flags) if flags else "-",
            "gate_actions": ";".join(actions) if actions else "-",
            "callable_bp_genome": fmt(int(x["callable_bp_genome"]) if x["callable_bp_genome"] is not None else None),
            "callable_frac_union_mask": fmt(x["callable_frac_union"]),
            "mean_depth": fmt(x["mean_depth"]),
            "pct_mapped_primary": fmt(x["pct_mapped"], 5),
            "pct_dup_primary": fmt(x["pct_dup"], 4),
            "n_base1_sites": fmt(c["n_sites"] if c else None),
            "n_called": fmt(c["n_called"] if c else None),
            "n_ref": fmt(c["n_ref"] if c else None),
            "n_alt": fmt(c["n_alt"] if c else None),
            "n_N": fmt(c["n_N"] if c else None),
            "n_mask_gq": g.get("n_mask_gq", NA),
            "n_mask_dp": g.get("n_mask_dp", NA),
            "n_mask_ad": g.get("n_mask_ad", NA),
            "n_mask_star": g.get("n_mask_star", NA),
            "n_ruleA_kept": g.get("n_ruleA_kept", NA),
        })

    with open(a.out_table, "w") as fh:
        fh.write("\t".join(out_cols) + "\n")
        for row in rows_out:
            fh.write("\t".join(row[k] for k in out_cols) + "\n")

    with open(a.out_iaf_sites, "w") as fh:
        fh.write("callset\tsample\tchrom\tpos\tref\talt\tgt_after_mask\tad\tsum_ad\tmajor_frac\tsrc\tadv\tqual\n")
        for row in iaf_detail:
            fh.write("\t".join(row) + "\n")

    # per-chromosome median depth
    with open(a.contigs) as fh:
        chroms = [l.strip() for l in fh if l.strip() and l.strip() != a.mito]
    with open(a.out_depth, "w") as fh:
        fh.write("sample\tgenome_median_depth\t" + "\t".join(chroms) + "\t" + "\t".join(f"{c}_ratio" for c in chroms) + "\n")
        for s in names:
            gmed = rec[s]["median_depth"]
            meds = [num(depth_chrom[s].get(c)) for c in chroms]
            ratios = [(m / gmed) if (m is not None and gmed) else None for m in meds]
            fh.write(s + "\t" + fmt(gmed) + "\t" + "\t".join(fmt(m) for m in meds) + "\t" + "\t".join(fmt(q, 4) for q in ratios) + "\n")

    # human-readable summary
    L = []
    L.append("G0 strain-level QC gate")
    L.append("=" * 60)
    L.append(f"samples in table          : {len(rows_out)}")
    L.append(f"reference set             : role == {a.ref_role}, n = {len(reps)}")
    L.append(f"analyzable region         : {a.analyzable_bp} bp (nuclear minus static mask)")
    L.append("")
    L.append("Thresholds")
    L.append(f"  callable fraction  >= {a.callable_ratio_min} x median(reps) = {fmt(ratio_min * med_callable_frac) if med_callable_frac else NA}  (median {fmt(med_callable_frac)})")
    L.append(f"  median depth       >= {a.depth_median_min:g}x")
    L.append(f"  distance to H99    |robust z| <= {a.refdist_z_max:g} on D_norm (ST5 only)")
    L.append("  MLST               stringMLST ST5 (pre-inclusion typing, samples.tsv group)")
    L.append("")
    L.append("Intermediate allele fraction (reported only: no threshold, no pass/fail, no flag)")
    L.append(f"  numerator   base-1 sites whose major-allele AD fraction is in [{a.iaf_lo}, {a.iaf_hi}]")
    L.append(f"  denominator base-1 sites with AD confirmed against the sample gVCF (FORMAT/ADV = M) and sum(AD) >= {a.iaf_min_ad}, whatever the GT; fraction NA when the denominator is 0")
    L.append("  Synthesized AD (reference blocks) and unverified AD are left out and counted in n_ad10_excluded_synthesized / n_ad10_excluded_unverified.")
    L.append("  No mixed-infection classification is made from these numbers.")
    rep_frac = [num(row["intermediate_af_frac"]) for row in rows_out if row["role"] == a.ref_role and num(row["intermediate_af_frac"]) is not None]
    if rep_frac:
        L.append(f"  representatives (n = {len(rep_frac)} with a denominator): median {fmt(median(rep_frac))}, min {fmt(min(rep_frac))}, max {fmt(max(rep_frac))}")
    L.append("")
    L.append("Distance to H99 (CNA3) among the representatives")
    L.append(f"  D_raw  median {fmt(med_draw)}  MAD {fmt(mad_draw)}")
    L.append(f"  D_norm median {fmt(med_dnorm)}  MAD {fmt(mad_dnorm)}  (median callable bp {fmt(med_callable_bp)})")
    L.append("  Read as clustering / outliers; absolute values are not expected to match 46,218.")
    L.append("  Flags are kept and no sample is excluded for them; the absolute differences below are reported with z. A flag is not read as sub-lineage structure without review.")
    L.append("")
    L.append("Failures per item (any sample): " + ", ".join(f"{k}={v}" for k, v in tally.items()))
    flagged = [row for row in rows_out if row["gate_overall"] != "PASS"]
    L.append(f"Flagged samples: {len(flagged)} / {len(rows_out)}")
    for row in flagged:
        extra = f"  D_raw-med={row['d_raw_minus_median']} D_norm-med={row['d_norm_minus_median']} z={row['refdist_robust_z']}" if row["pass_refdist"] == "FAIL" else ""
        L.append(f"  {row['sample']:<10} {row['role']:<28} {row['gate_flags']}  -> {row['gate_actions']}{extra}")
    L.append("")
    for cs in ("h99", "case"):
        L.append(f"Callset {cs}: " + ", ".join(f"{row['sample']} D_raw={row['d_raw_to_ref']} intermediate={row['n_intermediate_af_sites']}/{row['n_sites_ad10_measured']}={row['intermediate_af_frac']}" for row in rows_out if row["callset"] == cs))
    L.append(f"Intermediate-fraction genotypes of the h99 / case callsets, listed site by site: {len(iaf_detail)} (intermediate_af_sites_h99_case.tsv)")
    n_kept = sum(int(row["n_ruleA_kept"]) for row in rows_out if row["n_ruleA_kept"] not in (NA, ""))
    L.append(f"Rule A: {n_kept} reference-block REF calls kept because base 2 is callable at the base (see 06_filter/<callset>/refblock_kept.tsv.gz).")
    if a.mask_summary and os.path.exists(a.mask_summary):
        L.append("")
        L.append("Static mask + empirical inaccessible region (mask_union_summary.tsv)")
        with open(a.mask_summary) as fh:
            for line in fh:
                L.append("  " + line.rstrip("\n"))
    if warnings:
        L.append("")
        L.append(f"WARNINGS ({len(warnings)})")
        L.extend("  " + w for w in warnings)
    with open(a.out_summary, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
