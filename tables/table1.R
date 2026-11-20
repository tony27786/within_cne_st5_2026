# Table 1: characteristics of the 160 patients (one representative isolate each) by cryptococcal meningitis (CM) and HIV status; descriptive only.
# Inputs: association/design.tsv (build_design.py) and the additional columns of data/admissions.csv, joined by admission_id.
# Output: tables/out/Table1.tsv. Run: Rscript tables/table1.R

# Labels contain non-ASCII characters (written as \u escapes); Rscript may start in the C locale.
if (!l10n_info()[["UTF-8"]]) invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}
TAB_DIR <- get_script_dir()
PROJ <- normalizePath(file.path(TAB_DIR, ".."), mustWork = TRUE)
OUT_DIR <- Sys.getenv("TAB_OUTPUT_DIR", file.path(TAB_DIR, "out"))

P_DESIGN <- "association/design.tsv"
P_ADM <- "data/admissions.csv"
ADD_COLS <- c("sex", "headache", "fever", "vomiting", "neck_stiffness", "seizure", "visual_change",
              "hiv_viral_load", "hiv_viral_load_censor", "csf_chloride", "blood_culture", "csf_culture",
              "admission_date", "discharge_date")

abs_path <- function(rel) file.path(PROJ, rel)

d <- read.delim(abs_path(P_DESIGN), colClasses = "character", na.strings = character())
adm <- read.csv(abs_path(P_ADM), colClasses = "character", encoding = "UTF-8", na.strings = character(),
                check.names = FALSE)
names(adm)[1] <- sub("^\ufeff", "", names(adm)[1])
adm <- adm[, c("admission_id", "patient_id", "altered_consciousness", ADD_COLS)]

m <- adm[match(d$admission_id, adm$admission_id), ]
d <- cbind(d, m[, ADD_COLS])

num <- function(x) ifelse(x == "", NA_real_, suppressWarnings(as.numeric(x)))
is1 <- function(x) x == "1"

# Groups
grp_lab <- c(pos = "CM, HIV-positive", neg = "CM, HIV-negative", noncm = "Non-CM", all = "All")
d$grp <- ifelse(d$cm == "1" & d$hiv_positive == "1", "pos", ifelse(d$cm == "1" & d$hiv_positive == "0", "neg", "noncm"))
symptoms <- c("headache", "fever", "vomiting", "altered_consciousness", "neck_stiffness", "seizure", "visual_change")

# CSF culture: 0 is negative in CM patients and in non-CM patients with CSF chemistry; 0 without chloride, glucose or protein means not done.
chem <- !is.na(num(d$csf_chloride)) | !is.na(num(d$csf_glucose)) | !is.na(num(d$csf_protein))
d$csf_cult <- ifelse(is1(d$csf_culture), "positive", ifelse(d$cm == "1" | chem, "negative", "not done"))

hiv <- d$hiv_positive == "1"
vl_state <- ifelse(d$hiv_viral_load_censor == "", "not tested", d$hiv_viral_load_censor)
op <- num(d$opening_pressure)
op_cens <- d$opening_pressure_censored == "1"
los <- as.integer(as.Date(d$discharge_date) - as.Date(d$admission_date))
fu10 <- d$alive_10w; fu6 <- d$alive_6m

# Cell formats
groups <- list(pos = d$grp == "pos", neg = d$grp == "neg", noncm = d$grp == "noncm", all = rep(TRUE, nrow(d)))
# Round half up (away from zero) to d decimals; sprintf and formatC round the stored double (ties to even at exact halves: a median of 38.5 was shown as 38). "+ 0" turns -0 into 0.
fmt_num <- function(x, d = 0) formatC(sign(x) * floor(abs(x) * 10^d + 0.5 + 1e-9) / 10^d + 0, format = "f", digits = d)
pct <- function(n, N) { p <- 100 * n / N; ifelse(n > 0 & p < 0.5, "<1", fmt_num(p)) }

# n (%) among 'eligible' patients of the column; n/N (%) when N differs from the column size; "—" if none eligible.
cell_n <- function(event, eligible = rep(TRUE, nrow(d)), known = rep(TRUE, nrow(d))) {
  vapply(groups, function(g) {
    N <- sum(g & eligible & known)
    if (sum(g & eligible) == 0) return("\u2014")
    n <- sum(g & eligible & known & event)
    if (N == sum(g)) sprintf("%d (%s)", n, pct(n, N)) else sprintf("%d/%d (%s)", n, N, pct(n, N))
  }, character(1))
}
# count only (no percentage), for states reported as numbers of patients
cell_count <- function(event, eligible = rep(TRUE, nrow(d))) {
  vapply(groups, function(g) if (sum(g & eligible) == 0) "\u2014" else as.character(sum(g & eligible & event)), character(1))
}
# median (IQR); ", n = k" when not every eligible patient has a value
cell_med <- function(x, digits = 0, eligible = rep(TRUE, nrow(d))) {
  vapply(groups, function(g) {
    v <- x[g & eligible & !is.na(x)]
    if (!length(v)) return("\u2014")
    q <- quantile(v, c(0.5, 0.25, 0.75), names = FALSE)
    out <- sprintf("%s (%s\u2013%s)", fmt_num(q[1], digits), fmt_num(q[2], digits), fmt_num(q[3], digits))
    if (length(v) < sum(g & eligible)) out <- sprintf("%s, n = %d", out, length(v))
    out
  }, character(1))
}

rows <- list()
add <- function(label, cells = rep("", 4)) rows[[length(rows) + 1]] <<- c(label, cells)
sec <- function(label) add(label)

sec("Demographics")
add("Age, years, median (IQR)", cell_med(num(d$age_years)))
add("Sex, n (%)")
add("\u2003Male", cell_n(d$sex == "male"))
add("\u2003Female", cell_n(d$sex == "female"))

sec("HIV infection")
add("HIV positive, n (%)", cell_n(hiv))
add("On ART at admission, n (%)", cell_n(is1(d$art_at_admission), eligible = hiv))
add("HIV RNA tested, n (%)", cell_n(vl_state != "not tested", eligible = hiv))
add("\u2003Quantifiable, log\u2081\u2080 copies/mL, median (IQR)", cell_med(log10(num(d$hiv_viral_load)), 1, eligible = hiv & vl_state == "none"))
add("\u2003Quantifiable, n", cell_count(vl_state == "none", eligible = hiv))
add("\u2003Above the upper limit of quantification, n", cell_count(vl_state == "right", eligible = hiv))
add("\u2003Below the lower limit of quantification, n", cell_count(vl_state == "left", eligible = hiv))
add("\u2003Not detected, n", cell_count(vl_state == "undetectable", eligible = hiv))
add("CD4 count, cells/\u00b5L, median (IQR)", cell_med(num(d$cd4_count)))

sec("Signs and symptoms at admission")
sym_lab <- c(headache = "Headache", fever = "Fever", vomiting = "Vomiting", altered_consciousness = "Altered consciousness",
             neck_stiffness = "Neck stiffness", seizure = "Seizure", visual_change = "Visual change")
for (s in symptoms) add(sprintf("%s, n (%%)", sym_lab[[s]]), cell_n(is1(d[[s]])))

sec("Representative isolate")
add("Specimen, n (%)")
add("\u2003CSF", cell_n(d$specimen_group == "csf"))
add("\u2003Blood", cell_n(d$specimen_group == "blood"))
add("\u2003Other", cell_n(d$specimen_group == "other"))

sec("Microbiology")
add("Blood culture, positive record, n", cell_count(is1(d$blood_culture)))
add("CSF culture positive, n/N (%)", cell_n(d$csf_cult == "positive", known = d$csf_cult != "not done"))
# Only positive blood CrAg results were recorded (none negative), so the row counts positive records, as for blood culture
add("Blood CrAg, positive record, n", cell_count(d$blood_crag_result == "positive"))
add("CSF CrAg positive, n/N (%)", cell_n(d$csf_crag_result == "positive", known = d$csf_crag_result != ""))

sec("CSF at the first lumbar puncture")
add("White cells, \u00d710\u2076/L, median (IQR)", cell_med(num(d$csf_wbc)))
add("Protein, mg/L, median (IQR)", cell_med(num(d$csf_protein)))
add("Glucose, mmol/L, median (IQR)", cell_med(num(d$csf_glucose), 1))
add("Chloride, mmol/L, median (IQR)", cell_med(num(d$csf_chloride)))
add("Opening pressure, mmH\u2082O, median (IQR)", cell_med(op))
add("Opening pressure \u2265 250 mmH\u2082O, n/N (%)", cell_n(!is.na(op) & op >= 250, known = !is.na(op)))

sec("Hospital stay and outcome at discharge")
add("Length of hospital stay, days, median (IQR)", cell_med(los))
add("Outcome at discharge, n (%)")
out_lab <- c(improved = "Improved", died = "Died", discharged_against_advice = "Discharged against medical advice",
             transferred = "Transferred")
for (o in names(out_lab)) add(sprintf("\u2003%s", out_lab[[o]]), cell_n(d$discharge_outcome == o))

sec("Follow-up after admission")
add("All-cause death within 10 weeks, n/N (%)", cell_n(fu10 == "0", known = fu10 != ""))
add("All-cause death within 6 months, n/N (%)", cell_n(fu6 == "0", known = fu6 != ""))

tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(tab) <- c("Characteristic", sprintf("%s (n = %d)", grp_lab, vapply(groups, sum, integer(1))))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.table(tab, file.path(OUT_DIR, "Table1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Table 1 done: ", nrow(tab), " rows")
