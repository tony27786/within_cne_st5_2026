#!/usr/bin/env python3
"""gt_mask.py -- genotype-level N for base 1, with rule A. Standard library only.

Reads a haploid joint-called VCF on stdin and writes VCF on stdout. A called genotype becomes '.' when GQ < --gq-min, DP < --dp-min, AD[called allele] / sum(AD) < --ad-frac-min, or the allele is '*'.
Rule A (needs --src-dir from gt_source.py): a REF genotype from a gVCF reference block with GQ >= --gq-min and DP equal to the block MIN_DP is judged by base 2 at that base instead of by DP (FORMAT/SRC=BC).
FORMAT/SRC (V, B, BC, D, -) and FORMAT/ADV (M, S, A, X, U) record the genotype source and whether AD is confirmed against the sample gVCF; INFO/AC, AN and AF are recomputed.
"""
import argparse
import gzip
import os
import sys
from fractions import Fraction

SAMPLE_KEYS = (
    "n_sites",             # records seen
    "n_missing_in",        # GT already missing on input
    "n_called_in",         # GT called on input
    "n_masked",            # called on input -> set to N (any reason)
    "n_mask_gq",           # reason counts (non-exclusive)
    "n_mask_dp",
    "n_mask_ad",
    "n_mask_star",
    "n_ad_missing",        # AD absent/incomplete on a called genotype
    "n_ad_missing_dp_ok",  # ... while DP >= dp-min (expected 0)
    "n_called_not_major",  # called-allele rule masks, max(AD) rule would pass
    "n_src_V",             # genotype source counts (all genotypes)
    "n_src_B",
    "n_src_D",
    "n_src_uncovered",
    "n_ruleA_candidates",  # B, GT=0, GQ ok, DP < dp-min, DP == block MIN_DP
    "n_ruleA_kept",        # ... and base 2 callable -> DP rule not applied (SRC=BC)
    "n_ruleA_not_callable",  # ... base 2 not callable -> still N by DP
    "n_ruleA_kept_but_masked_other",  # kept by rule A for DP but masked by another rule
    "n_adv_M",             # AD evidence check (all genotypes): measured, matches the gVCF variant record
    "n_adv_S",             # synthesized (reference block)
    "n_adv_A",             # AD absent / incomplete
    "n_adv_X",             # variant record whose AD does not match
    "n_adv_U",             # not checkable (deletion-spanned, uncovered, different REF)
    "n_ref_out",           # remaining calls = REF
    "n_alt_out",           # remaining calls = non-REF
)

HEADER_LINES = {
    "AC": '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">',
    "AN": '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
    "AF": '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">',
}
SRC_HEADER = '##FORMAT=<ID=SRC,Number=1,Type=String,Description="G0 genotype source in the sample gVCF: V = variant record at this site (AD measured); B = reference block (DP = block MIN_DP, AD synthesized from MIN_DP, not a per-site allele count); BC = reference block kept by rule A (block MIN_DP < 10 but the base is callable in base 2: depth >= 10 from MAPQ >= 30 reads); D = spanned by an upstream deletion record; - = not covered by the gVCF; . = unknown">'
ADV_HEADER = '##FORMAT=<ID=ADV,Number=1,Type=String,Description="G0 AD evidence check against the sample gVCF: M = measured (variant record at this base with the same REF and an AD equal to the AD here, allele by allele); S = synthesized from a reference block MIN_DP (not a per-site allele count); A = AD absent or incomplete; X = variant record whose AD does not match; U = not checkable (deletion-spanned, not covered, different REF); . = unknown">'


def to_int(x):
    """Parse a VCF integer field; None for missing/unparsable."""
    if x is None or x == "" or x == ".":
        return None
    try:
        return int(x)
    except ValueError:
        try:
            return int(float(x))
        except ValueError:
            return None


def ad_check(src, alleles, ad, rec_ref, rec_alt, rec_ad):
    """FORMAT/ADV of one genotype (see the module docstring). ad = joint AD as ints, or None."""
    if src == ".":
        return "."
    if src == "B":
        return "S"
    if src != "V":
        return "U"
    if rec_ref != alleles[0]:
        return "U"
    rparts = [] if rec_ad in ("", ".") else rec_ad.split(",")
    if not rparts or any(to_int(x) is None for x in rparts) or ad is None:
        return "A"
    ralleles = [rec_ref] + rec_alt.split(",")
    if len(rparts) != len(ralleles):
        return "X"
    rad = dict(zip(ralleles, (to_int(x) for x in rparts)))
    nonref = rad.get("<NON_REF>", 0)
    for allele, n in zip(alleles, ad):
        if allele in rad:
            if n != rad[allele]:
                return "X"
        elif n > nonref:
            return "X"
    return "M"


def parse_ad(v, i_ad, n_alleles):
    """Joint AD of one genotype as ints, or None when absent / incomplete."""
    if i_ad < 0:
        return None
    parts = v[i_ad].split(",")
    vals = [to_int(x) for x in parts]
    if len(vals) == n_alleles and all(x is not None for x in vals):
        return vals
    return None


def set_info(info, updates):
    """Replace or append key=value pairs in a VCF INFO string, keeping order."""
    items = [] if info in ("", ".") else info.split(";")
    done = set()
    out = []
    for it in items:
        key = it.split("=", 1)[0]
        if key in updates:
            out.append(f"{key}={updates[key]}")
            done.add(key)
        else:
            out.append(it)
    for key, val in updates.items():
        if key not in done:
            out.append(f"{key}={val}")
    return ";".join(out) if out else "."


def main():
    ap = argparse.ArgumentParser(description="Genotype-level N for base 1, with rule A.")
    ap.add_argument("--gq-min", type=int, required=True)
    ap.add_argument("--dp-min", type=int, required=True)
    ap.add_argument("--ad-frac-min", required=True, help="parsed exactly, e.g. 0.90 -> 9/10")
    ap.add_argument("--src-dir", default=None, help="gt_source.py outputs <sample>.tsv.gz (enables SRC tags and rule A)")
    ap.add_argument("--kept-out", default=None, help="gz TSV of genotypes kept by rule A (required with --src-dir)")
    ap.add_argument("--stats", required=True, help="per-sample TSV")
    ap.add_argument("--site-stats", required=True, help="site-level TSV (metric, value)")
    args = ap.parse_args()
    if args.src_dir and not args.kept_out:
        sys.exit("gt_mask.py: --kept-out is required with --src-dir")

    thr = Fraction(args.ad_frac_min)
    num, den = thr.numerator, thr.denominator

    out = sys.stdout
    samples = []
    srcs = []
    stats = {}
    site = {"n_records": 0, "n_multiallelic": 0, "n_with_star_allele": 0, "n_no_alt_call_after_mask": 0}
    seen_info = set()
    kept = gzip.open(args.kept_out, "wt") if args.kept_out else None
    if kept:
        kept.write("chrom\tpos\tsample\tvcf_dp\tgq\tblock_start\tblock_end\tblock_min_dp\tblock_dp\n")

    for line in sys.stdin:
        if line.startswith("##"):
            if line.startswith("##INFO=<ID="):
                seen_info.add(line[len("##INFO=<ID="):].split(",", 1)[0])
            if line.startswith("##FORMAT=<ID=SRC,") or line.startswith("##FORMAT=<ID=ADV,"):
                continue
            out.write(line)
            continue
        if line.startswith("#CHROM"):
            for key in ("AC", "AN", "AF"):
                if key not in seen_info:
                    out.write(HEADER_LINES[key] + "\n")
            out.write(SRC_HEADER + "\n")
            out.write(ADV_HEADER + "\n")
            out.write("##G0_gt_mask=<Description=\"GT set to missing if GQ<%d, DP<%d, AD missing or sum(AD)==0, AD[GT]/sum(AD)<%s, or GT allele is '*'. Rule A: REF calls from a gVCF reference block with GQ>=%d and DP equal to the block MIN_DP are judged for depth by base 2 instead (SRC=BC). FORMAT/ADV = AD evidence check against the sample gVCF. AD/DP/GQ/PL unchanged; INFO AC/AN/AF recomputed from masked GT.\">\n" % (args.gq_min, args.dp_min, args.ad_frac_min, args.gq_min))
            out.write(line)
            samples = line.rstrip("\n").split("\t")[9:]
            stats = {s: dict.fromkeys(SAMPLE_KEYS, 0) for s in samples}
            if args.src_dir:
                for s in samples:
                    h = gzip.open(os.path.join(args.src_dir, f"{s}.tsv.gz"), "rt")
                    h.readline()
                    srcs.append(h)
            continue

        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            sys.exit(f"gt_mask.py: record without samples: {line[:80]}")
        alleles = [f[3]] + f[4].split(",")
        n_alleles = len(alleles)
        fmt = f[8].split(":")
        if fmt[0] != "GT":
            sys.exit(f"gt_mask.py: FORMAT must start with GT at {f[0]}:{f[1]}")
        if "SRC" in fmt or "ADV" in fmt:
            sys.exit(f"gt_mask.py: input already has FORMAT/SRC or ADV at {f[0]}:{f[1]}")
        i_ad = fmt.index("AD") if "AD" in fmt else -1
        i_dp = fmt.index("DP") if "DP" in fmt else -1
        i_gq = fmt.index("GQ") if "GQ" in fmt else -1
        n_fmt = len(fmt)
        f[8] = f[8] + ":SRC:ADV"

        ac = [0] * n_alleles
        an = 0
        for j, s in enumerate(samples):
            st = stats[s]
            st["n_sites"] += 1
            col = 9 + j
            v = f[col].split(":")
            if len(v) < n_fmt:
                v += ["."] * (n_fmt - len(v))
            if srcs:
                sl = srcs[j].readline().rstrip("\n").split("\t")
                if len(sl) < 12 or sl[0] != f[0] or sl[1] != f[1]:
                    sys.exit(f"gt_mask.py: source file of {s} out of step at {f[0]}:{f[1]} (or written by an older gt_source.py without rec_ref/rec_alt/rec_ad)")
                src, bstart, bend, bmin, bdp, callable_ = sl[2], sl[3], sl[4], sl[5], sl[6], sl[8]
                rec_ref, rec_alt, rec_ad = sl[9], sl[10], sl[11]
                st["n_src_" + {"V": "V", "B": "B", "D": "D"}.get(src, "uncovered")] += 1
            else:
                src, bstart, bend, bmin, bdp, callable_ = ".", ".", ".", ".", ".", "0"
                rec_ref, rec_alt, rec_ad = ".", ".", "."
            tag = src
            ad = parse_ad(v, i_ad, n_alleles)
            adv = ad_check(src, alleles, ad, rec_ref, rec_alt, rec_ad)
            if adv != ".":
                st["n_adv_" + adv] += 1
            gt = v[0]
            if gt in (".", "./.", ".|."):
                st["n_missing_in"] += 1
                v += [tag, adv]
                f[col] = ":".join(v)
                continue
            if "/" in gt or "|" in gt:
                sys.exit(f"gt_mask.py: non-haploid GT '{gt}' for {s} at {f[0]}:{f[1]} (G0 expects HaplotypeCaller --sample-ploidy 1)")
            a = to_int(gt)
            if a is None or a < 0 or a >= n_alleles:
                sys.exit(f"gt_mask.py: bad GT '{gt}' for {s} at {f[0]}:{f[1]}")
            st["n_called_in"] += 1

            gq = to_int(v[i_gq]) if i_gq >= 0 else None
            dp = to_int(v[i_dp]) if i_dp >= 0 else None

            reasons = []
            if gq is None or gq < args.gq_min:
                reasons.append("gq")
            dp_low = dp is None or dp < args.dp_min
            rule_a_kept = False
            if dp_low and src == "B" and a == 0 and gq is not None and gq >= args.gq_min and dp is not None and to_int(bmin) == dp:
                st["n_ruleA_candidates"] += 1
                if callable_ == "1":
                    rule_a_kept = True
                    dp_low = False
                    tag = "BC"
                    st["n_ruleA_kept"] += 1
                else:
                    st["n_ruleA_not_callable"] += 1
            if dp_low:
                reasons.append("dp")
            if ad is None:
                reasons.append("ad")
                st["n_ad_missing"] += 1
                if dp is not None and dp >= args.dp_min:
                    st["n_ad_missing_dp_ok"] += 1
            else:
                tot = sum(ad)
                if tot == 0 or ad[a] * den < num * tot:
                    reasons.append("ad")
                    if tot > 0 and max(ad) * den >= num * tot:
                        st["n_called_not_major"] += 1
            if alleles[a] == "*":
                reasons.append("star")

            if rule_a_kept:
                if reasons:
                    st["n_ruleA_kept_but_masked_other"] += 1
                elif kept:
                    kept.write(f"{f[0]}\t{f[1]}\t{s}\t{dp}\t{gq}\t{bstart}\t{bend}\t{bmin}\t{bdp}\n")
            if reasons:
                st["n_masked"] += 1
                for r in reasons:
                    st["n_mask_" + r] += 1
                v[0] = "."
            else:
                an += 1
                ac[a] += 1
                if a == 0:
                    st["n_ref_out"] += 1
                else:
                    st["n_alt_out"] += 1
            v += [tag, adv]
            f[col] = ":".join(v)

        site["n_records"] += 1
        if n_alleles > 2:
            site["n_multiallelic"] += 1
        if "*" in alleles[1:]:
            site["n_with_star_allele"] += 1
        if sum(ac[1:]) == 0:
            site["n_no_alt_call_after_mask"] += 1
        f[7] = set_info(f[7], {
            "AC": ",".join(str(x) for x in ac[1:]),
            "AN": str(an),
            "AF": ",".join((format(x / an, ".6g") if an else ".") for x in ac[1:]),
        })
        out.write("\t".join(f) + "\n")

    if not samples:
        sys.exit("gt_mask.py: no #CHROM header line on input")
    for h in srcs:
        if h.readline().strip():
            sys.exit("gt_mask.py: a source file has more lines than the VCF has records")
        h.close()
    if kept:
        kept.close()

    with open(args.stats, "w") as fh:
        fh.write("sample\t" + "\t".join(SAMPLE_KEYS) + "\n")
        for s in samples:
            fh.write(s + "\t" + "\t".join(str(stats[s][k]) for k in SAMPLE_KEYS) + "\n")
    with open(args.site_stats, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, val in site.items():
            fh.write(f"{k}\t{val}\n")
        tot = {k: sum(stats[s][k] for s in samples) for k in SAMPLE_KEYS}
        for k in ("n_called_in", "n_masked", "n_mask_gq", "n_mask_dp", "n_mask_ad", "n_mask_star", "n_ad_missing", "n_ad_missing_dp_ok", "n_called_not_major", "n_src_V", "n_src_B", "n_src_D", "n_src_uncovered", "n_ruleA_candidates", "n_ruleA_kept", "n_ruleA_not_callable", "n_ruleA_kept_but_masked_other", "n_adv_M", "n_adv_S", "n_adv_A", "n_adv_X", "n_adv_U"):
            fh.write(f"genotypes_{k}\t{tot[k]}\n")


if __name__ == "__main__":
    main()
