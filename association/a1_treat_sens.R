#!/usr/bin/env Rscript
# Post-hoc treatment-adjustment sensitivity analysis of the A1 mortality units C (10-week death) and D (6-month death): reproduces the Table 3 rows C and D from run_main.tsv (stops on mismatch), then refits the K-score bootstrap with the initial induction regimen in the null model (ridge penalty sigma 1.08 main, 0.58 and 1.41) and the h2 with the regimen as unpenalised fixed effects.
# Inputs: design table, K, admissions table (induction_regimen) and run_main.tsv (--ref). Outputs: treat_sens_{checks,models,desc,coef}.tsv in association/ (treat_sens_<tag>_{...}.tsv with --tag=<tag>).
# Usage (from the project root): Rscript association/a1_treat_sens.R   (optional: --design=, --kin=, --adm=, --ref=, --out=, --tag=<label> to add a label to the output names)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(root, "association/a1_lmm.R"))
source(file.path(root, "association/a1_treat_lib.R"))

design_path <- get_arg("design", file.path(root, "association/design.tsv"))
kin_path <- get_arg("kin", file.path(root, K_REL))
adm_path <- get_arg("adm", file.path(root, "data/admissions.csv"))
ref_path <- get_arg("ref", file.path(root, "association/run_main.tsv"))
out_dir <- get_arg("out", file.path(root, "association"))
tag <- get_arg("tag", "")
prefix <- file.path(out_dir, if (nzchar(tag)) paste0("treat_sens_", tag) else "treat_sens")

d0 <- read_design(design_path)
K <- read_kinship(kin_path)
adm <- read_admissions(adm_path)
ref <- read.delim(ref_path, colClasses = "character", na.strings = "", check.names = FALSE)
d <- derive_vars(d0)
units <- a1_units()
B <- 10000L
pus <- lapply(c(C = "C", D = "D"), prep_unit, units = units, d = d, K = K, adm = adm)

# 1 Reproduction of the Table 3 rows C and D (stops on mismatch)
# The run output stores numbers to 6 significant digits (write_tsv), so numbers are compared at that precision.
same6 <- function(a, b) isTRUE(is.finite(a) && is.finite(b) && abs(signif(a, 6) - b) <= 1e-12 * max(1, abs(b)))
base <- list()
checks <- list()
for (uid in names(pus)) {
  pu <- pus[[uid]]
  r <- ref[ref$unit == uid & ref$analysis == "main", , drop = FALSE]
  if (nrow(r) != 1) stop(sprintf("the frozen output has no single main row for unit %s", uid))
  rn <- function(col) as.numeric(r[[col]])
  seed0 <- A1_SEED + 1000L * pu$index + 0L
  f0 <- fit_lmm(pu$y, pu$X, pu$kin)
  sb0 <- score_boot(pu$y, pu$X, pu$kin$K, B = B, seed = seed0)
  gb0 <- boot_score(pu$y, glm_fitter(pu$X), pu$kin$K, B = B, seed = seed0)
  covs <- paste(attr(pu$x, "covs"), collapse = "+")
  counts0 <- c(sb0$B_valid, sb0$n_error, sb0$n_nonconv, sb0$n_sep_warn)
  ck <- function(item, frozen, new, ok) data.frame(unit = uid, item = item, frozen = as.character(frozen), reproduced = as.character(new),
                                                   ok = isTRUE(ok), stringsAsFactors = FALSE)
  checks[[uid]] <- rbind(
    ck("n", r$n, nrow(pu$x), nrow(pu$x) == rn("n")),
    ck("events", r$events, sum(pu$y), sum(pu$y) == rn("events")),
    ck("covariates", r$covariates, covs, identical(covs, r$covariates)),
    ck("c", r$c, signif(pu$kin$c, 6), same6(pu$kin$c, rn("c"))),
    ck("tau", r$tau, signif(f0$tau, 6), same6(f0$tau, rn("tau"))),
    ck("s2", r$s2, signif(f0$s2, 6), same6(f0$s2, rn("s2"))),
    ck("h2", r$h2, signif(f0$h2, 6), same6(f0$h2, rn("h2"))),
    ck("Q", r$stat, signif(sb0$score, 6), same6(sb0$score, rn("stat"))),
    ck("P (bootstrap)", r$p, signif(sb0$p_boot, 6), same6(sb0$p_boot, rn("p"))),
    ck("P (Davies)", r$p_davies, signif(sb0$p_davies, 6), same6(sb0$p_davies, rn("p_davies"))),
    ck("valid / error / non-converged / separation", paste(r$boot_valid, r$boot_error, r$boot_nonconv, r$boot_sep_warn, sep = " / "),
       paste(counts0, collapse = " / "),
       identical(counts0, as.integer(c(rn("boot_valid"), rn("boot_error"), rn("boot_nonconv"), rn("boot_sep_warn"))))),
    ck("generic bootstrap = frozen score_boot (Q, P, counts)", sprintf("%.10g / %.10g", sb0$score, sb0$p_boot),
       sprintf("%.10g / %.10g", gb0$score, gb0$p),
       identical(gb0$score, sb0$score) && identical(gb0$p, sb0$p_boot) &&
         identical(c(gb0$B_valid, gb0$n_error, gb0$n_nonconv, gb0$n_sep_warn), counts0)))
  base[[uid]] <- list(f0 = f0, sb0 = sb0, gb0 = gb0, seed = seed0)
}
checks <- do.call(rbind, checks)
rownames(checks) <- NULL
if (!all(checks$ok)) {
  print(checks[!checks$ok, ], row.names = FALSE)
  stop("reproduction of the frozen Table 3 rows failed; no treatment-adjusted model was fitted and nothing was written")
}
cat("reproduction of the frozen rows C and D: all checks passed\n")

# 2-3 Treatment-adjusted models and h2
models <- list()
coefs <- list()
h2s <- list()
for (uid in names(pus)) {
  pu <- pus[[uid]]
  b <- base[[uid]]
  models[[paste(uid, "baseline")]] <- data.frame(
    unit = uid, model = "baseline", sigma = NA_real_, seed = b$seed, n = nrow(pu$x), events = sum(pu$y), stat = b$sb0$score,
    p = b$sb0$p_boot, p_lower = b$gb0$p_lower, p_upper = b$gb0$p_upper, boot_valid = b$sb0$B_valid, boot_error = b$sb0$n_error,
    boot_nonconv = b$sb0$n_nonconv, boot_sep_warn = b$sb0$n_sep_warn, reportable = b$sb0$reportable, obs_status = "converged",
    obs_iter = NA_integer_, obs_capped = NA_integer_, obs_halvings = NA_integer_, obs_scaled_grad = NA_real_,
    obs_sep_warn = isTRUE(b$sb0$obs_sep_warn), sum_regimen_effects = NA_real_, tau = b$f0$tau, s2 = b$f0$s2, h2 = b$f0$h2,
    note = b$sb0$note, stringsAsFactors = FALSE)
  # h2: regimen as unpenalised fixed effects, largest group as reference; a second reference level as a check
  Xa <- cbind(pu$X, pu$Z[, -1, drop = FALSE])
  Xb <- cbind(pu$X, pu$Z[, -2, drop = FALSE])
  if (qr(Xa)$rank < ncol(Xa)) stop(sprintf("unit %s: host terms plus regimen indicators are not of full rank", uid))
  fa <- fit_lmm(pu$y, Xa, pu$kin)
  fb <- fit_lmm(pu$y, Xb, pu$kin)
  h2s[[uid]] <- data.frame(unit = uid, reference = pu$levels[1], tau = fa$tau, s2 = fa$s2, h2 = fa$h2, second_reference = pu$levels[2],
                           h2_second_reference = fb$h2, abs_difference = abs(fa$h2 - fb$h2), stringsAsFactors = FALSE)
  for (k in names(SIGMAS)) {
    sig <- SIGMAS[[k]]
    seed <- ts_seed(pu$index, SIGMA_CODE[[k]])
    fo <- pen_logit_fit(pu$y, pu$X, pu$Z, sig)
    gb <- boot_score(pu$y, pen_fitter(pu$X, pu$Z, sig), pu$kin$K, B = B, seed = seed)
    models[[paste(uid, k)]] <- data.frame(
      unit = uid, model = "treatment", sigma = sig, seed = seed, n = nrow(pu$x), events = sum(pu$y), stat = gb$score, p = gb$p,
      p_lower = gb$p_lower, p_upper = gb$p_upper, boot_valid = gb$B_valid, boot_error = gb$n_error, boot_nonconv = gb$n_nonconv,
      boot_sep_warn = gb$n_sep_warn, reportable = gb$reportable, obs_status = fo$status, obs_iter = fo$iter, obs_capped = fo$capped,
      obs_halvings = fo$halvings, obs_scaled_grad = fo$max_scaled_grad, obs_sep_warn = fo$sep,
      sum_regimen_effects = if (fo$converged) sum(fo$coef[colnames(pu$Z)]) else NA_real_, tau = fa$tau, s2 = fa$s2, h2 = fa$h2,
      note = gb$note, stringsAsFactors = FALSE)
    if (!is.null(fo$coef)) {
      coefs[[paste(uid, k)]] <- data.frame(unit = uid, sigma = sig, term = names(fo$coef), estimate = unname(fo$coef), stringsAsFactors = FALSE)
    }
    cat(sprintf("unit %s, sigma %.2f: P %s\n", uid, sig, format(signif(gb$p, 4))))
  }
}
models <- do.call(rbind, models)
coefs <- do.call(rbind, coefs)
h2s <- do.call(rbind, h2s)
rownames(models) <- rownames(coefs) <- rownames(h2s) <- NULL

# 4 Patients by regimen (unit C rows; unit D is nested)
xc <- pus$C$x
xd <- pus$D$x
x_all <- rbind(xc, xd[!(xd$sample %in% xc$sample), names(xc), drop = FALSE])
q3 <- function(v) quantile(v, c(0.5, 0.25, 0.75), names = FALSE, type = 7)
desc_row <- function(label, in_x, in_c, in_d) {
  z <- x_all[in_x, , drop = FALSE]
  show <- nrow(z) > 2
  cd4 <- if (show) q3(z$cd4) else rep(NA_real_, 3)
  age <- if (show) q3(z$age) else rep(NA_real_, 3)
  data.frame(regimen = label, n_10w = sum(in_c), deaths_10w = sum(xc$death10[in_c]), n_6m = sum(in_d), deaths_6m = sum(xd$death6[in_d]),
             n_described = nrow(z), hiv_pos = sum(z$hiv), cd4_median = cd4[1], cd4_q1 = cd4[2], cd4_q3 = cd4[3], age_median = age[1],
             age_q1 = age[2], age_q3 = age[3], year_min = if (show) min(z$year) else NA_real_, year_max = if (show) max(z$year) else NA_real_,
             icp_ge250 = sum(z$op_high %in% 1), icp_valid = sum(!is.na(z$op_high)), icp_missing = sum(is.na(z$op_high)),
             consc = sum(z$consc %in% 1), consc_valid = sum(!is.na(z$consc)), consc_missing = sum(is.na(z$consc)), stringsAsFactors = FALSE)
}
lv_all <- regimen_levels(x_all$regimen)
desc <- rbind(do.call(rbind, lapply(lv_all, function(l) desc_row(l, x_all$regimen == l, xc$regimen == l, xd$regimen == l))),
              desc_row("Total", rep(TRUE, nrow(x_all)), rep(TRUE, nrow(xc)), rep(TRUE, nrow(xd))))

write_tsv(checks, paste0(prefix, "_checks.tsv"))
write_tsv(models, paste0(prefix, "_models.tsv"))
write_tsv(desc, paste0(prefix, "_desc.tsv"))
write_tsv(coefs, paste0(prefix, "_coef.tsv"))
cat(sprintf("a1_treat_sens done: %s_{checks,models,desc,coef}.tsv\n", prefix))
