#!/usr/bin/env Rscript
# Host-model preparation for the A1 analysis (host covariates only; K and the genomic columns are not read): VIF per unit, linear / logistic host models for the four primary outcomes with influence diagnostics, functional form against natural splines and a patient-level bootstrap (2,000 resamples, seed 20261026), and fit checks for the other units.
# Input: design table (build_design.py). Outputs: host_prep_{vif,coef,boot,influence,feasibility,funcform}.tsv and host_prep_figs.pdf in association/ (host_prep_<tag>_{...} with --tag=<tag>).
# Usage (from the project root): Rscript association/host_prep.R   (optional: --design=, --out=, --tag=<label> to add a label to the output names)

suppressPackageStartupMessages({
  library(car)
  library(splines)
  library(boot)
})

SEED <- 20261026L
N_BOOT <- 2000L

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_path), ".."))
design_path <- get_arg("design", file.path(root, "association/design.tsv"))
out_dir <- get_arg("out", file.path(root, "association"))
tag <- get_arg("tag", "")
prefix <- file.path(out_dir, if (nzchar(tag)) paste0("host_prep_", tag) else "host_prep")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

d0 <- read.delim(design_path, colClasses = "character", na.strings = "", check.names = FALSE)
needed <- c("sample", "admission_id", "a1_main", "cm", "hiv_positive", "age_years", "sampling_year", "cd4_count",
            "art_at_admission", "csf_wbc", "csf_protein", "csf_glucose", "pop_csf", "opening_pressure",
            "opening_pressure_censored", "pop_op", "altered_consciousness", "alive_10w", "alive_6m",
            "csf_crag_pos_record", "csf_crag_titer_log2", "csf_crag_platform", "csf_crag_grade",
            "blood_crag_pos_record", "blood_crag_titer_log2", "blood_crag_platform", "blood_crag_grade")
miss_cols <- setdiff(needed, names(d0))
if (length(miss_cols)) stop("design table lacks columns: ", paste(miss_cols, collapse = ", "))

# Derived variables (A1 main set only)
num <- function(x) suppressWarnings(as.numeric(x))
d <- d0[d0$a1_main %in% "1", ]
d$cm <- as.integer(d$cm)
d$hiv <- as.integer(d$hiv_positive)
d$age <- num(d$age_years)
d$year <- num(d$sampling_year)
year_center <- mean(d$year)
d$year_c <- d$year - year_center
d$cd4 <- num(d$cd4_count)
if (any(d$cd4 <= 0, na.rm = TRUE)) stop("CD4 <= 0 present; ln(CD4) is undefined, decide the transform before running")
d$ln_cd4 <- log(d$cd4)
d$art <- as.integer(d$art_at_admission)
d$ln_wbc <- log(num(d$csf_wbc) + 1)
d$ln_prot <- log(num(d$csf_protein))
d$glu <- num(d$csf_glucose)
op <- num(d$opening_pressure)
d$op_high <- ifelse(is.na(op), NA, as.integer(d$opening_pressure_censored %in% "1" | op >= 250))
d$death10 <- 1L - as.integer(d$alive_10w)
d$death6 <- 1L - as.integer(d$alive_6m)
d$consc <- as.integer(d$altered_consciousness)
d$csf_pos <- as.integer(d$csf_crag_pos_record)
d$blood_pos <- as.integer(d$blood_crag_pos_record)
d$csf_lt <- num(d$csf_crag_titer_log2)
d$blood_lt <- num(d$blood_crag_titer_log2)
grade_score <- c(low = 1, mid = 2, high = 3)
d$csf_gs <- unname(grade_score[d$csf_crag_grade])
d$blood_gs <- unname(grade_score[d$blood_crag_grade])
d$csf_plat <- d$csf_crag_platform
d$blood_plat <- d$blood_crag_platform

# Analysis units
cov_cm <- c("hiv", "ln_cd4", "age", "year_c")
U <- function(id, role, y, family, covs, filt) list(id = id, role = role, y = y, family = family, covs = covs, filt = filt)
units <- list(
  U("A1", "primary", "ln_wbc", "gaussian", cov_cm, quote(cm == 1 & pop_csf %in% "1" & !is.na(ln_cd4))),
  U("A2", "primary", "ln_prot", "gaussian", cov_cm, quote(cm == 1 & pop_csf %in% "1" & !is.na(ln_cd4))),
  U("A3", "primary", "glu", "gaussian", cov_cm, quote(cm == 1 & pop_csf %in% "1" & !is.na(ln_cd4))),
  U("B", "primary", "op_high", "binomial", cov_cm, quote(cm == 1 & pop_op %in% "1" & !is.na(ln_cd4))),
  U("C", "prognostic", "death10", "binomial", cov_cm, quote(cm == 1 & !is.na(death10) & !is.na(ln_cd4))),
  U("D", "prognostic", "death6", "binomial", cov_cm, quote(cm == 1 & !is.na(death6) & !is.na(ln_cd4))),
  U("F", "exploratory", "consc", "binomial", cov_cm, quote(cm == 1 & !is.na(consc) & !is.na(ln_cd4))),
  U("G1", "exploratory", "csf_pos", "binomial", cov_cm, quote(cm == 1 & !is.na(csf_pos) & !is.na(ln_cd4))),
  U("G2", "exploratory", "csf_lt", "gaussian", cov_cm, quote(cm == 1 & csf_plat %in% "1:5" & !is.na(csf_lt) & !is.na(ln_cd4))),
  U("G3", "exploratory", "csf_lt", "gaussian", cov_cm, quote(cm == 1 & csf_plat %in% "1:2" & !is.na(csf_lt) & !is.na(ln_cd4))),
  U("G4", "exploratory", "csf_gs", "gaussian", c(cov_cm, "csf_plat"), quote(cm == 1 & !is.na(csf_gs) & !is.na(ln_cd4))),
  U("H1", "exploratory", "blood_pos", "binomial", c(cov_cm, "cm"), quote(!is.na(blood_pos) & !is.na(ln_cd4))),
  U("H2", "exploratory", "blood_lt", "gaussian", c(cov_cm, "cm"), quote(blood_plat %in% "1:5" & !is.na(blood_lt) & !is.na(ln_cd4))),
  U("H3", "exploratory", "blood_lt", "gaussian", c(cov_cm, "cm"), quote(blood_plat %in% "1:2" & !is.na(blood_lt) & !is.na(ln_cd4))),
  U("H4", "exploratory", "blood_gs", "gaussian", c(cov_cm, "cm", "blood_plat"), quote(!is.na(blood_gs) & !is.na(ln_cd4))),
  U("I", "exploratory", "ln_cd4", "gaussian", c("age", "year_c", "cm", "art"), quote(hiv == 1 & !is.na(ln_cd4) & !is.na(art)))
)
names(units) <- vapply(units, `[[`, "", "id")
unit_data <- function(u) {
  x <- d[with(d, eval(u$filt)), , drop = FALSE]
  x <- x[stats::complete.cases(x[, c(u$y, u$covs)]), , drop = FALSE]
  rownames(x) <- x$sample  # diagnostic plots label points by sample
  x
}
# Factor covariates with a single level in a subset are dropped and noted
usable_covs <- function(u, x) {
  keep <- vapply(u$covs, function(v) length(unique(x[[v]])) > 1, TRUE)
  u$covs[keep]
}

fit_model <- function(y, covs, family, x) {
  f <- reformulate(covs, response = y)
  warns <- character(0)
  fit <- withCallingHandlers(
    if (family == "gaussian") lm(f, data = x) else glm(f, data = x, family = binomial()),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  attr(fit, "warns") <- warns
  fit
}
glm_flags <- function(fit) {
  if (!inherits(fit, "glm")) return(list(converged = TRUE, separation = FALSE, max_abs_coef = NA, extreme_fitted = NA))
  p <- fitted(fit)
  mac <- max(abs(coef(fit)), na.rm = TRUE)
  ext <- sum(p < 1e-6 | p > 1 - 1e-6)
  list(converged = fit$converged, separation = (!fit$converged) || mac > 10 || ext > 0, max_abs_coef = mac, extreme_fitted = ext)
}
# Aliased (inestimable) terms are kept as rows with NA so that they stay visible
coef_table <- function(fit, unit, model) {
  s <- summary(fit)$coefficients
  ci <- if (inherits(fit, "glm")) suppressMessages(confint.default(fit)) else confint(fit)
  tm <- names(coef(fit))
  i <- match(tm, rownames(s))
  j <- match(tm, rownames(ci))
  data.frame(unit = unit, model = model, term = tm, estimate = s[i, 1], se = s[i, 2],
             ci_lo = ci[j, 1], ci_hi = ci[j, 2], p = s[i, 4], row.names = NULL)
}
write_tsv <- function(df, path) {
  df[] <- lapply(df, function(v) if (is.numeric(v)) signif(v, 6) else v)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

vif_rows <- list()
for (u in units) {
  x <- unit_data(u)
  covs <- usable_covs(u, x)
  if (length(covs) < 2) next
  x$.pseudo <- seq_len(nrow(x))
  v <- car::vif(lm(reformulate(covs, response = ".pseudo"), data = x))
  if (is.matrix(v)) v <- v[, ncol(v)]^2
  vif_rows[[u$id]] <- data.frame(unit = u$id, term = names(v), vif = unname(v))
}
vif_df <- do.call(rbind, vif_rows)
write_tsv(vif_df, paste0(prefix, "_vif.tsv"))

# Host models for the primary outcomes: coefficients, functional form, influence diagnostics, bootstrap
coef_rows <- list(); ff_rows <- list(); infl_rows <- list(); boot_rows <- list()
pdf(paste0(prefix, "_figs.pdf"), width = 8.5, height = 8.5)
for (k in seq_along(units)) {
  u <- units[[k]]
  if (u$role != "primary") next
  x <- unit_data(u)
  n <- nrow(x)
  fit <- fit_model(u$y, u$covs, u$family, x)
  coef_rows[[paste(u$id, "main")]] <- coef_table(fit, u$id, "main")
  fit_nocd4 <- fit_model(u$y, setdiff(u$covs, "ln_cd4"), u$family, x)
  coef_rows[[paste(u$id, "noCD4")]] <- coef_table(fit_nocd4, u$id, "same_rows_without_ln_cd4")

  # Functional form against natural splines
  for (term in c("age", "year_c", "ln_cd4")) {
    covs_alt <- c(setdiff(u$covs, term), sprintf("splines::ns(%s, df = 3)", term))
    alt <- fit_model(u$y, covs_alt, u$family, x)
    a <- anova(fit, alt, test = if (u$family == "gaussian") "F" else "Chisq")
    ff_rows[[paste(u$id, term)]] <- data.frame(unit = u$id, term = term, p_spline_vs_linear = a[2, ncol(a)], delta_aic = AIC(alt) - AIC(fit))
  }

  # Influence diagnostics (reported, never dropped)
  cooks <- cooks.distance(fit)
  dfb <- dfbetas(fit)[, -1, drop = FALSE]
  maxdfb <- apply(abs(dfb), 1, max)
  which_term <- colnames(dfb)[apply(abs(dfb), 1, which.max)]
  flag <- cooks > 4 / n | maxdfb > 2 / sqrt(n)
  if (any(flag)) {
    infl <- data.frame(unit = u$id, sample = x$sample[flag], admission_id = x$admission_id[flag], cooks_d = cooks[flag],
                       max_abs_dfbeta = maxdfb[flag], dfbeta_term = which_term[flag], row.names = NULL)
    infl <- infl[order(-infl$cooks_d), ]
    infl_rows[[u$id]] <- infl
    # descriptive refit without every flagged case; inestimable terms are shown, not dropped
    x_drop <- x[!flag, ]
    fit_drop <- fit_model(u$y, u$covs, u$family, x_drop)
    coef_rows[[paste(u$id, "drop")]] <- coef_table(fit_drop, u$id, "without_flagged_cases")
  }

  # Patient-level bootstrap of the host coefficients (percentile intervals)
  set.seed(SEED + k)
  terms0 <- names(coef(fit))
  stat <- function(dat, idx) {
    f2 <- fit_model(u$y, u$covs, u$family, dat[idx, ])
    cf <- coef(f2)[terms0]
    c(cf, flag = as.numeric(glm_flags(f2)$separation))
  }
  b <- boot::boot(x, stat, R = N_BOOT)
  bt <- b$t[, seq_along(terms0), drop = FALSE]
  bs <- data.frame(unit = u$id, term = terms0, estimate = b$t0[seq_along(terms0)],
                   boot_median = apply(bt, 2, median, na.rm = TRUE),
                   boot_lo = apply(bt, 2, quantile, 0.025, na.rm = TRUE), boot_hi = apply(bt, 2, quantile, 0.975, na.rm = TRUE),
                   same_sign = apply(sweep(sign(bt), 2, sign(b$t0[seq_along(terms0)]), "=="), 2, mean, na.rm = TRUE),
                   n_na = colSums(is.na(bt)), row.names = NULL)
  boot_rows[[u$id]] <- bs

  # figures
  op_par <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  plot(fit, which = c(1, 2, 4, 5), caption = paste(u$id, c("Residuals vs fitted", "Q-Q", "Scale-location", "Cook's distance", "Residuals vs leverage", "Cook's distance vs leverage")))
  par(op_par)
  car::crPlots(fit, terms = ~ age + year_c + ln_cd4, layout = c(2, 2), main = paste(u$id, "partial residuals"))
}
invisible(dev.off())
ff <- do.call(rbind, ff_rows)

# Fit checks for the other units
feas <- list()
for (u in units) {
  if (u$role == "primary") next
  x <- unit_data(u)
  covs <- usable_covs(u, x)
  note <- if (length(covs) < length(u$covs)) paste("single-level covariate dropped:", paste(setdiff(u$covs, covs), collapse = ",")) else ""
  if (nrow(x) <= length(covs) + 1) {
    feas[[u$id]] <- data.frame(unit = u$id, n = nrow(x), events = NA, converged = NA, separation_flag = NA, full_rank = NA, warnings = NA, note = paste(note, "sample size too small to fit"))
    next
  }
  fit <- fit_model(u$y, covs, u$family, x)
  fl <- glm_flags(fit)
  feas[[u$id]] <- data.frame(unit = u$id, n = nrow(x), events = if (u$family == "binomial") sum(x[[u$y]]) else NA,
                             converged = fl$converged, separation_flag = fl$separation, full_rank = !any(is.na(coef(fit))),
                             warnings = length(attr(fit, "warns")), note = note)
}
feas_df <- do.call(rbind, feas)

# Outputs
coef_df <- do.call(rbind, coef_rows)
write_tsv(coef_df, paste0(prefix, "_coef.tsv"))
write_tsv(do.call(rbind, boot_rows), paste0(prefix, "_boot.tsv"))
write_tsv(if (length(infl_rows)) do.call(rbind, infl_rows) else data.frame(unit = character(0)), paste0(prefix, "_influence.tsv"))
write_tsv(feas_df, paste0(prefix, "_feasibility.tsv"))
write_tsv(ff, paste0(prefix, "_funcform.tsv"))
cat(sprintf("host_prep done: %s_figs.pdf + %s_{vif,coef,boot,influence,feasibility,funcform}.tsv\n", prefix, prefix))
