# Shared functions for the A1 analysis: variable derivation, analysis units, LMM fit, RLRT null, profile interval and parametric-bootstrap score test.
# Model: y = Xb + g + e, g ~ N(0, tau K_S), e ~ N(0, s2 I); h2 = c tau / (c tau + s2) with c = tr(K_S) / n_S and K unscaled (gaston's h2 = tau / (tau + s2); to_gaston_h2() converts).
# Sourced by a1_run.R and a1_treat_sens.R; derive_vars() and the unit filters mirror host_prep.R so that the populations match.

suppressPackageStartupMessages({
  library(gaston)
  library(RLRsim)
})

A1_SEED <- 20261026L
G1_OUT <- "results/g1"   # G1 output directory, relative to the project root
K_REL <- file.path(G1_OUT, "03_l1", "K.tsv")
EIG_MIN <- 1e-6  # lmm.diago floors eigenvalues below 1e-6 and lmm.diago.likelihood does not; such K_S are refused

read_design <- function(design_path) read.delim(design_path, colClasses = "character", na.strings = "", check.names = FALSE)

read_kinship <- function(path) {
  K <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
  if (!identical(rownames(K), colnames(K)) || !isSymmetric(unname(K))) stop("K must be square, symmetric and named identically on both margins")
  K
}

derive_vars <- function(d0) {
  num <- function(x) suppressWarnings(as.numeric(x))
  d <- d0[d0$a1_main %in% "1", ]
  d$cm <- as.integer(d$cm)
  d$hiv <- as.integer(d$hiv_positive)
  d$age <- num(d$age_years)
  d$year <- num(d$sampling_year)
  year_center <- mean(d$year)
  d$year_c <- d$year - year_center
  d$cd4 <- num(d$cd4_count)
  if (any(d$cd4 <= 0, na.rm = TRUE)) stop("CD4 <= 0 present; ln(CD4) is undefined")
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
  d$blood_pos <- as.integer(d$blood_crag_pos_record)
  d$csf_lt <- num(d$csf_crag_titer_log2)
  d$blood_lt <- num(d$blood_crag_titer_log2)
  grade_score <- c(low = 1, mid = 2, high = 3)
  d$csf_gs <- unname(grade_score[d$csf_crag_grade])
  d$blood_gs <- unname(grade_score[d$blood_crag_grade])
  d$csf_plat <- d$csf_crag_platform
  d$blood_plat <- d$blood_crag_platform
  attr(d, "year_center") <- year_center
  d
}

# type: continuous (RLRT), ordinal (RLRT, approximate: scores 1/2/3), binary (parametric bootstrap score test)
a1_units <- function() {
  cov_cm <- c("hiv", "ln_cd4", "age", "year_c")
  U <- function(id, role, y, type, covs, filt) list(id = id, role = role, y = y, type = type, covs = covs, filt = filt)
  units <- list(
    U("A1", "primary", "ln_wbc", "continuous", cov_cm, quote(cm == 1 & pop_csf %in% "1" & !is.na(ln_cd4))),
    U("A2", "primary", "ln_prot", "continuous", cov_cm, quote(cm == 1 & pop_csf %in% "1" & !is.na(ln_cd4))),
    U("A3", "primary", "glu", "continuous", cov_cm, quote(cm == 1 & pop_csf %in% "1" & !is.na(ln_cd4))),
    U("B", "primary", "op_high", "binary", cov_cm, quote(cm == 1 & pop_op %in% "1" & !is.na(ln_cd4))),
    U("C", "prognostic", "death10", "binary", cov_cm, quote(cm == 1 & !is.na(death10) & !is.na(ln_cd4))),
    U("D", "prognostic", "death6", "binary", cov_cm, quote(cm == 1 & !is.na(death6) & !is.na(ln_cd4))),
    U("F", "exploratory", "consc", "binary", cov_cm, quote(cm == 1 & !is.na(consc) & !is.na(ln_cd4))),
    U("G2", "exploratory", "csf_lt", "continuous", cov_cm, quote(cm == 1 & csf_plat %in% "1:5" & !is.na(csf_lt) & !is.na(ln_cd4))),
    U("G3", "exploratory", "csf_lt", "continuous", cov_cm, quote(cm == 1 & csf_plat %in% "1:2" & !is.na(csf_lt) & !is.na(ln_cd4))),
    U("G4", "exploratory", "csf_gs", "ordinal", c(cov_cm, "csf_plat"), quote(cm == 1 & !is.na(csf_gs) & !is.na(ln_cd4))),
    U("H1", "exploratory", "blood_pos", "binary", c(cov_cm, "cm"), quote(!is.na(blood_pos) & !is.na(ln_cd4))),
    U("H2", "exploratory", "blood_lt", "continuous", c(cov_cm, "cm"), quote(blood_plat %in% "1:5" & !is.na(blood_lt) & !is.na(ln_cd4))),
    U("H3", "exploratory", "blood_lt", "continuous", c(cov_cm, "cm"), quote(blood_plat %in% "1:2" & !is.na(blood_lt) & !is.na(ln_cd4))),
    U("H4", "exploratory", "blood_gs", "ordinal", c(cov_cm, "cm", "blood_plat"), quote(!is.na(blood_gs) & !is.na(ln_cd4))),
    U("I", "exploratory", "ln_cd4", "continuous", c("age", "year_c", "cm", "art"), quote(hiv == 1 & !is.na(ln_cd4) & !is.na(art)))
  )
  names(units) <- vapply(units, `[[`, "", "id")
  units
}

# Rows of a unit; covariates with a single level in the subset are dropped and returned in attr "dropped"
unit_data <- function(u, d) {
  x <- d[with(d, eval(u$filt)), , drop = FALSE]
  x <- x[stats::complete.cases(x[, c(u$y, u$covs)]), , drop = FALSE]
  rownames(x) <- x$sample
  keep <- vapply(u$covs, function(v) length(unique(x[[v]])) > 1, TRUE)
  attr(x, "covs") <- u$covs[keep]
  attr(x, "dropped") <- u$covs[!keep]
  x
}

unit_X <- function(x) {
  X <- model.matrix(reformulate(attr(x, "covs")), data = x)
  colnames(X) <- make.names(colnames(X))
  X
}

# K_S matched to the unit rows by sample name, with its eigen decomposition and scale c
unit_kin <- function(K, x) {
  s <- x$sample
  miss <- setdiff(s, rownames(K))
  if (length(miss)) stop("samples missing from K: ", paste(head(miss), collapse = ", "))
  Ks <- K[s, s, drop = FALSE]
  eig <- eigen(Ks, symmetric = TRUE)
  if (min(eig$values) <= EIG_MIN) stop(sprintf("K_S smallest eigenvalue %.3g <= %g; not positive definite enough for lmm.diago", min(eig$values), EIG_MIN))
  list(K = Ks, eig = eig, c = mean(diag(Ks)), min_eig = min(eig$values))
}

to_gaston_h2 <- function(h2, c) {
  r <- h2 / (c * (1 - h2))
  ifelse(h2 >= 1, 1, r / (1 + r))
}

fit_lmm <- function(y, X, kin) {
  f <- lmm.diago(Y = y, X = X, eigenK = kin$eig, verbose = FALSE)
  l1 <- as.numeric(lmm.diago.likelihood(tau = f$tau, s2 = f$sigma2, Y = y, X = X, eigenK = kin$eig))
  l0 <- as.numeric(lmm.diago.likelihood(h2 = 0, Y = y, X = X, eigenK = kin$eig)$likelihood)
  list(tau = f$tau, s2 = f$sigma2, h2 = kin$c * f$tau / (kin$c * f$tau + f$sigma2),
       beta = setNames(as.numeric(f$BLUP_beta), colnames(X)), logL1 = l1, logL0 = l0, rlrt = max(0, 2 * (l1 - l0)))
}

# Profile restricted-likelihood interval for h2 (2 * drop <= crit), truncated to [0, 1], on the grid 0, 0.01, ..., 1 refined by uniroot; max_gap > 0 would mean the grid beat the fitted maximum.
profile_ci <- function(y, X, kin, fit, crit = qchisq(0.95, 1)) {
  prof <- function(h) as.numeric(lmm.diago.likelihood(h2 = to_gaston_h2(h, kin$c), Y = y, X = X, eigenK = kin$eig)$likelihood)
  dev <- function(h) 2 * (fit$logL1 - prof(h)) - crit
  grid <- seq(0, 1, by = 0.01)
  dv <- dev(grid)
  h <- fit$h2
  below <- which(grid < h & dv > 0)
  lower <- if (!length(below)) 0 else {
    i <- max(below)
    uniroot(dev, c(grid[i], min(grid[i + 1], h)), tol = 1e-6)$root
  }
  above <- which(grid > h & dv > 0)
  upper <- if (!length(above)) 1 else {
    j <- min(above)
    uniroot(dev, c(max(grid[j - 1], h), grid[j]), tol = 1e-6)$root
  }
  list(lower = lower, upper = upper, max_gap = max(0, max(prof(grid)) - fit$logL1))
}

# RLRTSim null with Z built from the n - p non-zero eigenpairs of (I - H) K (I - H), so that the p near-zero directions inside span(X) do not add chi-square terms.
rlrt_null <- function(X, kin, nsim = 10000L, seed = A1_SEED) {
  n <- nrow(X)
  qrX <- qr(X)
  PK <- qr.resid(qrX, t(qr.resid(qrX, kin$K)))
  e <- eigen((PK + t(PK)) / 2, symmetric = TRUE)
  k <- n - qrX$rank
  Z <- e$vectors[, seq_len(k), drop = FALSE] %*% diag(sqrt(pmax(e$values[seq_len(k)], 0)), k)
  as.numeric(RLRTSim(X = X, Z = Z, qrX = qrX, sqrt.Sigma = diag(k), nsim = nsim, seed = seed))
}

p_from_null <- function(stat, null, tol = 1e-8) (1 + sum(null >= stat - tol)) / (1 + length(null))

# Logistic variance-component score statistic s = (y - pi0)' K (y - pi0) / 2 with pi0 from the host null model (as gaston::score.variance.logistic); P from a parametric bootstrap that refits the null model in every replicate.
# Failed or non-converged replicates are excluded and counted, separation warnings are kept and counted; the unit is not reportable (P = NA) if the observed null fit fails or fewer than min_valid * B replicates remain.
score_boot <- function(y, X, Ks, B = 10000L, seed = A1_SEED, min_valid = 0.99) {
  null_fit <- function(yy) {
    warns <- character(0)
    fit <- withCallingHandlers(tryCatch(glm.fit(X, yy, family = binomial()), error = function(e) NULL),
                               warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(fit = fit, sep = any(grepl("fitted probabilities numerically 0 or 1", warns)), other_warn = sum(!grepl("fitted probabilities numerically 0 or 1", warns)))
  }
  score_of <- function(yy, fit) {
    r <- yy - fit$fitted.values
    drop(crossprod(r, Ks %*% r)) / 2
  }
  obs <- null_fit(y)
  if (is.null(obs$fit) || !obs$fit$converged) {
    return(list(score = NA, p_boot = NA, p_davies = NA, B = B, B_valid = 0L, n_error = NA, n_nonconv = NA, n_sep_warn = NA,
                n_other_warn = NA, obs_converged = FALSE, obs_sep_warn = obs$sep, reportable = FALSE,
                note = "observed null model failed or did not converge"))
  }
  s_obs <- score_of(y, obs$fit)
  dav <- score.variance.logistic(K0 = unname(Ks), Y = y, X = X)
  if (abs(as.numeric(dav$score) - s_obs) > 1e-6 * max(1, s_obs)) stop("score statistic differs from gaston::score.variance.logistic")
  pi0 <- obs$fit$fitted.values
  set.seed(seed)
  sb <- numeric(B); n_error <- 0L; n_nonconv <- 0L; n_sep <- 0L; n_warn <- 0L
  for (b in seq_len(B)) {
    yb <- rbinom(length(y), 1, pi0)
    nf <- null_fit(yb)
    if (is.null(nf$fit)) { sb[b] <- NA; n_error <- n_error + 1L; next }
    if (!nf$fit$converged) { sb[b] <- NA; n_nonconv <- n_nonconv + 1L; next }
    n_sep <- n_sep + nf$sep
    n_warn <- n_warn + (nf$other_warn > 0)
    sb[b] <- score_of(yb, nf$fit)
  }
  ok <- !is.na(sb)
  reportable <- sum(ok) >= min_valid * B
  list(score = s_obs, p_boot = if (reportable) (1 + sum(sb[ok] >= s_obs - 1e-10)) / (1 + sum(ok)) else NA,
       p_davies = as.numeric(dav$p), B = B, B_valid = sum(ok), n_error = n_error, n_nonconv = n_nonconv, n_sep_warn = n_sep,
       n_other_warn = n_warn, obs_converged = TRUE, obs_sep_warn = obs$sep, reportable = reportable,
       note = if (reportable) "" else sprintf("valid replicates %d < %d", sum(ok), ceiling(min_valid * B)))
}

write_tsv <- function(df, path) {
  df[] <- lapply(df, function(v) if (is.numeric(v)) signif(v, 6) else v)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
