# Table 2 (the five atypical isolate pairs: two inter-patient pairs at 0 SNPs, three intra-patient pairs > 5 SNPs) and the supplementary tables with per-isolate checks of the nine isolates involved (TableS_exceptions) and the Cn13 allele depths (TableS_cn13).
# Inputs: G1 output (02_l0 distances and ruler B pairs, SNP address, 01_prep stats, supported clades), G0 07_qc_gate (qc_gate, depth_by_chrom), data/strain_hiv.tsv, isolates table, errata, and the Cn13 allele-depth check results/misscheck/cn13_mix_check.tsv.
# Outputs: tables/out/{Table2, TableS_exceptions, TableS_cn13}.tsv. Run: Rscript tables/table2.R

if (!l10n_info()[["UTF-8"]]) invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
TAB_DIR <- get_script_dir()
PROJ <- normalizePath(file.path(TAB_DIR, ".."), mustWork = TRUE)
OUT_DIR <- Sys.getenv("TAB_OUTPUT_DIR", file.path(TAB_DIR, "out"))
abs_path <- function(rel) file.path(PROJ, rel)

G1 <- "results/g1/"
G0 <- "results/g0/07_qc_gate/"
MC <- "results/misscheck/"
IN <- c(
  molten = paste0(G1, "02_l0/dists/main.molten.tsv"),
  ruler_b = paste0(G1, "02_l0/main/ruler_B_pairs.tsv"),
  address = paste0(G1, "02_l0/main/snp_address.tsv"),
  stats = paste0(G1, "01_prep/sets/main/stats.tsv"),
  clades = paste0(G1, "04_tree/supported_clades.tsv"),
  qc = paste0(G0, "qc_gate.tsv"),
  depth = paste0(G0, "depth_by_chrom.tsv"),
  roles = "data/strain_hiv.tsv",
  isolates = "data/isolates.csv",
  cn13 = paste0(MC, "cn13_mix_check.tsv"),
  errata = "data/errata_clinical.tsv"
)

rd <- function(k, ...) read.delim(abs_path(IN[[k]]), check.names = FALSE, stringsAsFactors = FALSE, ...)

molten <- read.delim(abs_path(IN[["molten"]]), header = FALSE, col.names = c("a", "b", "d"), stringsAsFactors = FALSE)
rb <- rd("ruler_b")
address <- rd("address")
stats <- rd("stats")
clades <- rd("clades")
qc <- rd("qc")
depth <- rd("depth")
roles <- rd("roles")
iso <- read.csv(abs_path(IN[["isolates"]]), colClasses = "character", encoding = "UTF-8", check.names = FALSE)
names(iso)[1] <- sub("^\ufeff", "", names(iso)[1])
iso <- iso[, c("sample", "patient_id", "sampling_date", "specimen_type", "wgs_role")]
# Errata to the isolate table are applied only if the archived value matches.
errata <- read.delim(abs_path(IN[["errata"]]), colClasses = "character", na.strings = character())
ei <- errata[errata$table == "isolates", ]
for (k in seq_len(nrow(ei))) {
  r <- which(iso[[ei$key_field[k]]] == ei$key_value[k])
  stopifnot(length(r) == 1, ei$field[k] %in% names(iso), iso[r, ei$field[k]] == ei$archived_value[k])
  iso[r, ei$field[k]] <- ei$corrected_value[k]
}

# Cohort, distance matrix and helpers
cohort <- roles$sample[roles$wgs_role %in% c("representative", "within_patient_extra")]   # 191 ST5 isolates, no Cn169b
reps <- roles$sample[roles$wgs_role == "representative"]
pat <- setNames(roles$patient_id, roles$sample)

smp <- sort(unique(molten$a))
D <- matrix(NA_real_, length(smp), length(smp), dimnames = list(smp, smp))
D[cbind(molten$a, molten$b)] <- molten$d
dist <- function(x, y) D[x, y]

# Nearest isolates from another patient, among the 191; Cn13 is not a candidate because its mixed sites are
# called as N, which can underestimate its pairwise distances. Ties are all kept.
nearest_other <- function(s) {
  cand <- setdiff(cohort[pat[cohort] != pat[[s]]], "Cn13")
  v <- D[s, cand]
  hit <- names(v)[v == min(v)]
  paste(sprintf("%s (%s), %d", hit, pat[hit], as.integer(min(v))), collapse = "; ")
}

pair_tab <- data.frame(
  a = c("Cn13", "Cn101", "Cn103", "Cn212", "Cn174"),
  b = c("Cn14", "Cn103", "Cn111", "Cn211", "Cn191"),
  type = c("between", "between", "within", "within", "within"),
  stringsAsFactors = FALSE
)
pair_tab$d <- mapply(dist, pair_tab$a, pair_tab$b)

others_174 <- setdiff(cohort, c("Cn174", "Cn191"))

gt <- setNames(stats$n_gt_missing, stats$sample)
gt_med <- median(gt[cohort]); gt_p95 <- unname(quantile(gt[cohort], 0.95))
iaf <- setNames(qc$intermediate_af_frac, qc$sample)
iaf_pct <- function(s) 100 * mean(iaf[cohort] < iaf[[s]])   # percentage of the 191 isolates with a lower value

ratio_cols <- grep("_ratio$", names(depth), value = TRUE)
depth_flag <- function(s) {
  r <- unlist(depth[depth$sample == s, ratio_cols])
  out <- sub("_ratio$", "", ratio_cols[r < 0.85 | r > 1.15])
  if (length(out)) paste(out, collapse = ", ") else "None"
}
iso9 <- c("Cn13", "Cn14", "Cn101", "Cn103", "Cn111", "Cn212", "Cn211", "Cn174", "Cn191")
# QC-gate flag for the distance to H99 (|robust z| > 3.5); FAIL only means outside the range, so the direction
# comes from the deviation from the cohort median.
refflag <- setNames(qc$pass_refdist, qc$sample)
refdev <- setNames(qc$d_raw_minus_median, qc$sample)
ref_label <- function(s) ifelse(refflag[s] != "FAIL", "Within range", ifelse(refdev[s] > 0, "High (flagged)", "Low (flagged)"))

clade_tips <- strsplit(clades$tips_below, ",")
smallest_clade <- function(s) {
  if (!s %in% reps) return("Not in the tree (additional isolate)")
  n <- clades$n_tips_below[vapply(clade_tips, function(t) s %in% t, logical(1))]
  if (length(n)) as.character(min(n)) else "None"
}

# Laboratory-record check of the nine isolates in the laboratory information system (manual), with remarks from the original reports.
lis_checked <- c("Cn13", "Cn14", "Cn101", "Cn103", "Cn111", "Cn212", "Cn211", "Cn174", "Cn191")
lis_note <- c(Cn212 = "positive in blood cultures from both arms", Cn211 = "sputum graded unsatisfactory by the laboratory",
              Cn174 = "7 colonies")

# Cn13 allele depths at sites where Cn14 and a comparator (Cn95 or Cn211) are both called and differ
mix <- read.delim(abs_path(IN[["cn13"]]), header = FALSE, colClasses = "character",
                  col.names = c("comp", "chrom", "pos", "gt14", "gtcomp", "gt13", "ad13"))
ad <- lapply(strsplit(mix$ad13, ","), as.numeric)
mix$depth <- vapply(ad, sum, numeric(1))
mix$f14 <- mapply(function(a, g) a[as.integer(g) + 1], ad, mix$gt14) / mix$depth   # Cn14-allele read fraction
mix$minor <- vapply(ad, function(a) 1 - max(a) / sum(a), numeric(1))
mix_sum <- do.call(rbind, lapply(c("Cn95", "Cn211"), function(r) {
  m <- mix[mix$comp == r, ]
  nn <- m$gt13 == "."
  data.frame(comp = r, sites = nrow(m), same14 = sum(!nn & m$gt13 == m$gt14), samecomp = sum(!nn & m$gt13 == m$gtcomp),
             n_N = sum(nn), f14_med = median(m$f14[nn]), f14_q1 = unname(quantile(m$f14[nn], 0.25)),
             f14_q3 = unname(quantile(m$f14[nn], 0.75)), maj14 = sum(m$f14[nn] > 0.5),
             minor_mid = sum(m$minor[nn] >= 0.3 & m$minor[nn] <= 0.5), min_depth = min(m$depth[nn]),
             d13 = dist("Cn13", r), d14 = dist("Cn14", r))
}))

# Table 2
spec <- function(s) {
  x <- iso$specimen_type[match(s, iso$sample)]
  out <- c("Cerebrospinal fluid" = "CSF", "Whole blood" = "blood", "Drainage fluid" = "drainage fluid", "Sputum" = "sputum")[x]
  ifelse(is.na(out), tolower(x), out)
}
sdate <- function(s) iso$sampling_date[match(s, iso$sample)]
days <- function(a, b) as.integer(abs(as.Date(sdate(b)) - as.Date(sdate(a))))
admission <- function(a, b) {
  r <- rb[(rb$isolate_1 == a & rb$isolate_2 == b) | (rb$isolate_1 == b & rb$isolate_2 == a), ]
  if (r$admission_axis == "same") "same admission" else "two admissions"
}
sampling <- function(i) {
  a <- pair_tab$a[i]; b <- pair_tab$b[i]
  sp <- if (spec(a) == spec(b)) sprintf("both %s", spec(a)) else sprintf("%s and %s", spec(a), spec(b))
  if (pair_tab$type[i] == "between") sprintf("%d days apart, %s", days(a, b), sp)   # full dates are not published
  else sprintf("%d days apart, %s, %s", days(a, b), admission(a, b), sp)
}

evidence <- c(
  sprintf("Asymmetric distances to Cn95 (%d vs. %d SNPs). Cn13 has the highest intermediate allele fraction among the 191 ST5 isolates",
          dist("Cn13", "Cn95"), dist("Cn14", "Cn95")),
  sprintf("Both %d SNPs from Cn95. Intermediate allele fraction and missing genotype calls within the cohort range",
          dist("Cn101", "Cn95")),
  sprintf("Cn103 is %d SNPs from Cn101 (P029). The nearest inter-patient neighbor of Cn111 is Cn211 (%d SNPs)",
          dist("Cn103", "Cn101"), dist("Cn111", "Cn211")),
  sprintf("Cn211 is %d SNPs from Cn95 (P020, sampled in %s) and has %d missing genotype calls (cohort 95th percentile, %.1f)",
          dist("Cn211", "Cn95"), substr(sdate("Cn95"), 1, 4), gt[["Cn211"]], gt_p95),
  sprintf("Each other's nearest neighbor, \u2265 %d SNPs from all other patients' isolates. Cn174 is in the 15-isolate clade",
          as.integer(min(D["Cn174", setdiff(others_174, "Cn13")])))
)
interpretation <- c(
  "Allele depths consistent with a two-genotype mixture in Cn13 (about 6:4), the major genotype matching Cn14 at the sites examined, with mixed sites called as N",
  "Indistinguishable at mutually callable sites. The same strain in two patients cannot be distinguished from a post-archiving error",
  "Undetermined", "Undetermined", "Undetermined"
)

tab2 <- data.frame(
  Pair = sprintf("%s\u2013%s (%s)", pair_tab$a, pair_tab$b,
                 ifelse(pair_tab$type == "between", paste(pat[pair_tab$a], pat[pair_tab$b], sep = ", "), pat[pair_tab$a])),
  Sampling = vapply(seq_len(nrow(pair_tab)), sampling, character(1)),
  SNPs = as.character(as.integer(pair_tab$d)),
  `Key genomic evidence` = evidence,
  Interpretation = interpretation,
  check.names = FALSE, stringsAsFactors = FALSE
)
grp_row <- function(label) data.frame(Pair = label, Sampling = "", SNPs = "", `Key genomic evidence` = "", Interpretation = "",
                                      check.names = FALSE)
tab2 <- rbind(grp_row("Inter-patient pairs (0 SNPs)"), tab2[1:2, ], grp_row("Intra-patient pairs (> 5 SNPs)"), tab2[3:5, ])
rownames(tab2) <- NULL

# Supplementary table: per-isolate checks
role_lab <- c(representative = "Representative", within_patient_extra = "Additional")
qrow <- function(s) qc[qc$sample == s, ]
tabS <- data.frame(
  Isolate = iso9,
  Patient = pat[iso9],
  Role = role_lab[roles$wgs_role[match(iso9, roles$sample)]],
  `Sampling (year-month)` = substr(sdate(iso9), 1, 7),
  Specimen = spec(iso9),
  `Intermediate allele sites, n/N (%)` = vapply(iso9, function(s) {
    r <- qrow(s); sprintf("%d/%s (%.2f)", r$n_intermediate_af_sites, format(r$n_sites_ad10_measured, big.mark = ","), 100 * r$intermediate_af_frac)
  }, character(1)),
  `Intermediate allele fraction, percentile` = sprintf("%.1f", vapply(iso9, iaf_pct, numeric(1))),
  `Missing genotype calls` = as.integer(gt[iso9]),
  `Relative chromosomal depth outside 0.85–1.15` = vapply(iso9, depth_flag, character(1)),
  `Nearest inter-patient isolate (SNPs)` = vapply(iso9, nearest_other, character(1)),
  `Smallest supported clade (isolates)` = vapply(iso9, smallest_clade, character(1)),
  `SNP address` = ifelse(iso9 %in% reps, address$address[match(iso9, address$sample)], "\u2014"),
  `Normalized distance to H99` = unname(ref_label(iso9)),
  `Laboratory record` = ifelse(!iso9 %in% lis_checked, "Not checked", ifelse(iso9 %in% names(lis_note), paste0("Checked, ", lis_note[iso9]), "Checked")),
  check.names = FALSE, stringsAsFactors = FALSE
)
rownames(tabS) <- NULL

# Supplementary table: allele depths of Cn13 at sites that separate Cn14 from a comparator
tabC <- with(mix_sum, data.frame(
  Comparator = sprintf("%s (%s)", comp, pat[comp]),
  `Distance to comparator (SNPs), Cn13 / Cn14` = sprintf("%d / %d", as.integer(d13), as.integer(d14)),
  `Sites where Cn14 and comparator differ` = sites,
  `Cn13 called as Cn14` = same14,
  `Cn13 called as comparator` = samecomp,
  `Cn13 missing (N)` = n_N,
  `N sites, Cn14 allele read fraction, median (IQR)` = sprintf("%.2f (%.2f\u2013%.2f)", f14_med, f14_q1, f14_q3),
  `N sites, Cn14 allele in the majority` = sprintf("%d/%d", maj14, n_N),
  `N sites, minor allele fraction 0.3–0.5` = sprintf("%d/%d", minor_mid, n_N),
  check.names = FALSE, stringsAsFactors = FALSE
))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.table(tab2, file.path(OUT_DIR, "Table2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tabS, file.path(OUT_DIR, "TableS_exceptions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tabC, file.path(OUT_DIR, "TableS_cn13.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Table 2 done")
