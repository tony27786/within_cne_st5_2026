# Table 3 (primary and prognostic outcomes), Table S-A1a (exploratory outcomes) and Table S-A1b (sensitivity analyses) of the genome-phenotype analysis, generated from the A1 run outputs.
# Inputs: association/run_{main,sensitivity}.tsv (a1_run.R). Outputs: tables/out/{Table3, TableS_A1a, TableS_A1b}.tsv.
# Number format: P and q with two significant digits (two decimals when >= 0.1); h2 and interval bounds with two decimals, an exact 0 as "0", and 0 < h2 < 0.01 with one significant digit. Run: Rscript tables/table_a1.R

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

IN <- c(main = "association/run_main.tsv", sens = "association/run_sensitivity.tsv")

rd <- function(k) read.delim(abs_path(IN[[k]]), check.names = FALSE, stringsAsFactors = FALSE)
main <- rd("main")
sens <- rd("sens")

fmt_p <- function(p) ifelse(is.na(p), "\u2014", ifelse(p >= 0.1, sprintf("%.2f", p),
                                                       formatC(signif(p, 2), format = "fg", digits = 2, flag = "#")))
fmt_h2 <- function(h) ifelse(h == 0, "0", ifelse(h < 0.01, formatC(signif(h, 1), format = "fg", digits = 1), sprintf("%.2f", h)))
fmt_b <- function(b) ifelse(b == 0, "0", sprintf("%.2f", b))
h2_cell <- function(r) {
  if (r$type == "continuous") sprintf("%s (%s–%s)", fmt_h2(r$h2), fmt_b(r$ci_lo), fmt_b(r$ci_hi))
  else fmt_h2(r$h2)
}
n_cell <- function(r) if (r$type == "binary") sprintf("%d (%d)", r$n, as.integer(r$events)) else as.character(r$n)
test_lab <- c(binary = "Score, bootstrap", continuous = "RLRT", ordinal = "RLRT, approximate")

# Table 3: primary outcomes (BH over the four) and prognostic outcomes (prespecified, uncorrected)
m <- function(u) main[main$unit == u, ]

lab3 <- c(A1 = "CSF white cell count, \u00d710\u2076/L, ln(count + 1)", A2 = "CSF protein, mg/L, ln", A3 = "CSF glucose, mmol/L",
          B = "Opening pressure \u2265 250 mmH\u2082O", C = "10-week all-cause mortality", D = "6-month all-cause mortality")
row3 <- function(u) data.frame(Outcome = lab3[[u]], n = n_cell(m(u)), h2ci = h2_cell(m(u)), P = fmt_p(m(u)$p),
                               qbh = if (u %in% c("C", "D")) "\u2014" else fmt_p(m(u)$q), stringsAsFactors = FALSE)
hdr3 <- function(t) data.frame(Outcome = t, n = "", h2ci = "", P = "", qbh = "", stringsAsFactors = FALSE)
tab3 <- rbind(hdr3("Primary outcomes"), row3("A1"), row3("A2"), row3("A3"), row3("B"),
              hdr3("Prognostic outcomes"), row3("C"), row3("D"))
rownames(tab3) <- NULL
names(tab3) <- c("Outcome", "n (events)", "h\u00b2 (approximate 95% CI)", "P", "q")

# Table S-A1a: exploratory outcomes
labA <- c(F = "Altered consciousness",
          G2 = "CSF CrAg titer, log\u2082, 1:5 series", G3 = "CSF CrAg titer, log\u2082, 1:2 series",
          G4 = "CSF CrAg merged grade (1/2/3)", H1 = "Blood CrAg positive record",
          H2 = "Blood CrAg titer, log\u2082, 1:5 series", H3 = "Blood CrAg titer, log\u2082, 1:2 series",
          H4 = "Blood CrAg merged grade (1/2/3)", I = "ln(CD4 count), HIV-positive patients")
popA <- c(F = "CM", G2 = "CM", G3 = "CM", G4 = "CM", H1 = "All", H2 = "All", H3 = "All", H4 = "All", I = "HIV+")
tabA <- do.call(rbind, lapply(names(labA), function(u) {
  r <- m(u)
  data.frame(Outcome = labA[[u]], Population = popA[[u]],
             `n (events)` = n_cell(r), h2ci = h2_cell(r),
             P = fmt_p(r$p), q = fmt_p(r$q), Test = test_lab[[r$type]],
             check.names = FALSE, stringsAsFactors = FALSE)
}))
names(tabA)[names(tabA) == "h2ci"] <- "h\u00b2 (approximate 95% CI)"

# Table S-A1b: sensitivity analyses for the primary outcomes
lab_an <- c(S1b_same_rows_no_cd4 = "Same patients, no CD4 adjustment", S1c_all_rows_no_cd4 = "All eligible CM patients, no CD4 adjustment",
            S2_rank_inverse_normal = "Rank-based inverse normal transformation", "S3a_without_Cn101_Cn103 (1 removed)" = "Without Cn103")
lab_out <- c(A1 = "CSF white cell count", A2 = "CSF protein", A3 = "CSF glucose", B = "Opening pressure \u2265 250 mmH\u2082O")
tabB <- do.call(rbind, lapply(seq_len(nrow(sens)), function(i) {
  r <- sens[i, ]
  first <- i == 1 || sens$unit[i - 1] != r$unit
  data.frame(Outcome = if (first) lab_out[[r$unit]] else "", Analysis = lab_an[[r$analysis]], n = n_cell(r),
             h2ci = h2_cell(r), P = fmt_p(r$p), check.names = FALSE, stringsAsFactors = FALSE)
}))
names(tabB)[names(tabB) == "h2ci"] <- "h\u00b2 (approximate 95% CI)"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write_tab <- function(tab, file) write.table(tab, file.path(OUT_DIR, paste0(file, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write_tab(tab3, "Table3")
write_tab(tabA, "TableS_A1a")
write_tab(tabB, "TableS_A1b")
message("A1 tables done")
