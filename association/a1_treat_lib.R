# Functions for the treatment-adjustment sensitivity analysis of the A1 mortality units: regimen coding, ridge-penalised logistic null model, generic K-score parametric bootstrap and the admissions reader.
# Sourced by a1_treat_sens.R after a1_lmm.R.

TS_SEED <- 20261109L                    # new analyses: TS_SEED + 1000 x unit index + sigma code
SIGMAS <- c(main = 1.08, strong = 0.58, weak = 1.41)
SIGMA_CODE <- c(main = 1L, strong = 2L, weak = 3L)
GLM_EPS <- 10 * .Machine$double.eps     # glm.fit's threshold for "fitted probabilities numerically 0 or 1"
SEP_MSG <- "fitted probabilities numerically 0 or 1"
REGIMENS <- c("AmB+5FC+FLU", "AmB+5FC+VRC", "AmB+5FC", "5FC+VRC", "5FC+FLU", "AmB+FLU", "AmB")  # labels in Table 3 rows

ts_seed <- function(unit_index, code) TS_SEED + 1000L * unit_index + code

# Levels ordered by group size (largest first), ties by label, so the order is deterministic
regimen_levels <- function(reg) {
  tab <- table(reg)
  names(tab)[order(-as.integer(tab), names(tab))]
}

# Full dummy coding: one 0/1 column per regimen level
regimen_Z <- function(reg, levels) {
  Z <- matrix(0, length(reg), length(levels), dimnames = list(NULL, paste0("reg:", levels)))
  for (j in seq_along(levels)) Z[, j] <- as.numeric(reg == levels[j])
  if (any(rowSums(Z) != 1)) stop("every row needs exactly one regimen level")
  Z
}

# y 0/1; X unpenalised columns with the intercept first; Z regimen indicators penalised with 1 / sigma^2 (ridge, an N(0, sigma^2) prior on each regimen effect; the effects sum to zero at the optimum).
# status: "converged", "nonconv" (iteration limit, or step halving could not improve) or "error"; sep flags fitted probabilities numerically 0 or 1 as glm.fit.
pen_logit_fit <- function(y, X, Z, sigma, start = NULL, tol_dev = 1e-8, tol_grad = 1e-6, maxit = 50L, max_halve = 30L,
                          max_eta_step = 10) {
  M <- cbind(X, Z)
  p <- ncol(M)
  lam <- c(rep(0, ncol(X)), rep(1 / sigma^2, ncol(Z)))
  capped <- 0L
  result <- function(status, th = NULL, it = 0L, halv = 0L, sg = NA_real_, msg = "") {
    fitted <- if (is.null(th)) NULL else plogis(drop(M %*% th))
    list(status = status, converged = identical(status, "converged"), coef = if (is.null(th)) NULL else setNames(th, colnames(M)),
         fitted = fitted, iter = it, halvings = halv, capped = capped, max_scaled_grad = sg,
         sep = !is.null(fitted) && any(fitted < GLM_EPS | fitted > 1 - GLM_EPS), message = msg)
  }
  if (length(y) != nrow(M) || !all(is.finite(M)) || !all(y %in% c(0, 1))) {
    return(result("error", msg = "non-finite design, or outcome not 0/1"))
  }
  if (!all(X[, 1] == 1)) return(result("error", msg = "the first column of X must be the intercept"))
  # penalised log-likelihood, with log(1 + exp(eta)) computed stably
  objective <- function(th) {
    eta <- drop(M %*% th)
    if (!all(is.finite(eta))) return(-Inf)
    sum(y * eta - (pmax(eta, 0) + log1p(exp(-abs(eta))))) - 0.5 * sum(lam * th^2)
  }
  grad_hess <- function(th) {
    pr <- plogis(drop(M %*% th))
    w <- pr * (1 - pr)
    list(g = drop(crossprod(M, y - pr)) - lam * th, H = crossprod(M * sqrt(w)) + diag(lam, nrow = p))
  }
  scaled_grad <- function(gh) max(abs(gh$g) / sqrt(diag(gh$H)))
  th <- if (is.null(start)) c(qlogis((sum(y) + 0.5) / (length(y) + 1)), rep(0, p - 1)) else start
  f <- objective(th)
  if (!is.finite(f)) return(result("error", msg = "non-finite objective at the start"))
  dev_old <- -2 * f
  halv <- 0L
  for (it in seq_len(maxit)) {
    gh <- grad_hess(th)
    if (!all(is.finite(gh$g)) || !all(is.finite(gh$H))) return(result("error", th, it, halv, msg = "non-finite gradient or Hessian"))
    R <- tryCatch(chol(gh$H), error = function(e) NULL)
    if (is.null(R)) return(result("error", th, it, halv, msg = "Hessian not positive definite"))
    step <- backsolve(R, forwardsolve(t(R), gh$g))
    # Shorten (never turn) a step that moves any linear predictor by more than max_eta_step logit units.
    d_eta <- max(abs(M %*% step))
    if (is.finite(d_eta) && d_eta > max_eta_step) {
      step <- step * (max_eta_step / d_eta)
      capped <- capped + 1L
    }
    improves <- function(fn) is.finite(fn) && fn >= f - 1e-12 * abs(f)
    s <- 1
    f_new <- objective(th + step)
    k <- 0L
    while (!improves(f_new) && k < max_halve) {
      s <- s / 2
      f_new <- objective(th + s * step)
      k <- k + 1L
    }
    halv <- halv + k
    if (!improves(f_new)) return(result("nonconv", th, it, halv, scaled_grad(gh), "step halving did not improve the objective"))
    th <- th + s * step
    f <- f_new
    dev <- -2 * f
    if (abs(dev - dev_old) / (abs(dev) + 0.1) < tol_dev) {
      sg <- scaled_grad(grad_hess(th))
      if (is.finite(sg) && sg < tol_grad) return(result("converged", th, it, halv, sg))
    }
    dev_old <- dev
  }
  result("nonconv", th, maxit, halv, scaled_grad(grad_hess(th)), "iteration limit reached")
}

# Null-model fitters for the bootstrap; each returns list(status = "converged" / "nonconv" / "error", fitted, sep). glm_fitter is the null fit of score_boot (glm.fit with default control, warnings captured).
glm_fitter <- function(X) {
  function(yy) {
    warns <- character(0)
    fit <- withCallingHandlers(tryCatch(glm.fit(X, yy, family = binomial()), error = function(e) NULL),
                               warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    if (is.null(fit)) return(list(status = "error", fitted = NULL, sep = FALSE))
    list(status = if (fit$converged) "converged" else "nonconv", fitted = fit$fitted.values, sep = any(grepl(SEP_MSG, warns)))
  }
}

pen_fitter <- function(X, Z, sigma) {
  function(yy) tryCatch(pen_logit_fit(yy, X, Z, sigma),
                        error = function(e) list(status = "error", fitted = NULL, sep = FALSE, message = conditionMessage(e)))
}

# K-score parametric bootstrap: Q = (y - pi)' Ks (y - pi) / 2 with pi from the fitted null model; replicates are simulated from the observed fit and refitted with the same fitter.
# With the glm fitter and the same seed this reproduces score_boot draw for draw; seed = NULL continues from the current random-number state.
boot_score <- function(y, fitter, Ks, B, seed = NULL, min_valid = 0.99, keep_draws = FALSE) {
  q_of <- function(yy, fitted) {
    r <- yy - fitted
    drop(crossprod(r, Ks %*% r)) / 2
  }
  obs <- fitter(y)
  if (!identical(obs$status, "converged")) {
    return(list(score = NA_real_, p = NA_real_, p_lower = NA_real_, p_upper = NA_real_, B = B, B_valid = 0L, n_error = NA_integer_,
                n_nonconv = NA_integer_, n_sep_warn = NA_integer_, n_exceed = NA_integer_, obs_status = obs$status,
                obs_sep_warn = isTRUE(obs$sep), reportable = FALSE, note = "observed null model failed or did not converge"))
  }
  s_obs <- q_of(y, obs$fitted)
  pi0 <- obs$fitted
  if (!is.null(seed)) set.seed(seed)
  sb <- rep(NA_real_, B)
  n_error <- 0L
  n_nonconv <- 0L
  n_sep <- 0L
  for (b in seq_len(B)) {
    yb <- rbinom(length(y), 1, pi0)
    f <- fitter(yb)
    if (identical(f$status, "error")) { n_error <- n_error + 1L; next }
    if (!identical(f$status, "converged")) { n_nonconv <- n_nonconv + 1L; next }
    n_sep <- n_sep + isTRUE(f$sep)
    sb[b] <- q_of(yb, f$fitted)
  }
  ok <- !is.na(sb)
  V <- sum(ok)
  E <- sum(sb[ok] >= s_obs - 1e-10)
  reportable <- V >= min_valid * B
  res <- list(score = s_obs, p = if (reportable) (1 + E) / (1 + V) else NA_real_, p_lower = (1 + E) / (1 + B),
              p_upper = (1 + E + (B - V)) / (1 + B), B = B, B_valid = V, n_error = n_error, n_nonconv = n_nonconv,
              n_sep_warn = n_sep, n_exceed = E, obs_status = "converged", obs_sep_warn = isTRUE(obs$sep), reportable = reportable,
              note = if (reportable) "" else sprintf("valid replicates %d < %d", V, ceiling(min_valid * B)))
  if (keep_draws) res$draws <- sb
  res
}

# Admissions table reduced to admission_id and induction_regimen; the byte-order mark is stripped from the first column name.
read_admissions <- function(path) {
  a <- read.csv(path, colClasses = "character", check.names = FALSE, na.strings = "", encoding = "UTF-8")
  names(a)[1] <- sub("^\xef\xbb\xbf", "", names(a)[1], useBytes = TRUE)
  need <- c("admission_id", "induction_regimen")
  if (!all(need %in% names(a))) stop("admissions table lacks admission_id or induction_regimen")
  a <- a[, need]
  if (anyNA(a$admission_id) || anyDuplicated(a$admission_id)) stop("admission_id missing or duplicated in the admissions table")
  a
}

# Unit rows from the a1_lmm.R functions, with the regimen of the representative strain's own admission attached
prep_unit <- function(uid, units, d, K, adm) {
  u <- units[[uid]]
  x <- unit_data(u, d)
  reg <- adm$induction_regimen[match(x$admission_id, adm$admission_id)]
  bad <- is.na(reg) | !nzchar(reg)
  if (any(bad)) stop(sprintf("unit %s: %d rows without an initial regimen", uid, sum(bad)))
  x$regimen <- reg
  lv <- regimen_levels(reg)
  list(id = uid, index = match(uid, names(units)), u = u, x = x, y = x[[u$y]], X = unit_X(x), kin = unit_kin(K, x), reg = reg,
       levels = lv, Z = regimen_Z(reg, lv))
}
