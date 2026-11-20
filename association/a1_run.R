#!/usr/bin/env Rscript
# A1 tests of the 15 modelled units (RLRT against the RLRTSim null with a profile interval for continuous units, the same RLRT without interval for ordinal units, parametric-bootstrap score test for binary units) and sensitivity analyses of the four primary units.
# Inputs: design table (build_design.py) and the kinship matrix K of the G1 output. Outputs: run_{main,sensitivity,loo}.tsv in association/ (run_<tag>_{...}.tsv with --tag=<tag>).
# Usage (from the project root): Rscript association/a1_run.R   (optional: --design=, --kin=, --out=, --tag=<label> to add a label to the output names)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(root, "association/a1_lmm.R"))

design_path <- get_arg("design", file.path(root, "association/design.tsv"))
kin_path <- get_arg("kin", file.path(root, K_REL))
out_dir <- get_arg("out", file.path(root, "association"))
tag <- get_arg("tag", "")
prefix <- file.path(out_dir, if (nzchar(tag)) paste0("run_", tag) else "run")

d0 <- read_design(design_path)
K <- read_kinship(kin_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
d <- derive_vars(d0)
units <- a1_units()
NSIM <- 10000L
B <- 10000L

with_covs <- function(x, covs) {
  attr(x, "covs") <- covs
  attr(x, "dropped") <- character(0)
  x
}
seed_of <- function(uid, code) A1_SEED + 1000L * match(uid, names(units)) + code

# One analysis of one unit: h2, its interval (continuous only) and the P value by the unit's test
analyze <- function(u, x, analysis, code, y = x[[u$y]]) {
  X <- unit_X(x)
  kin <- unit_kin(K, x)
  f <- fit_lmm(y, X, kin)
  out <- data.frame(unit = u$id, analysis = analysis, type = u$type, n = nrow(x),
                    events = if (u$type == "binary") sum(y) else NA, covariates = paste(attr(x, "covs"), collapse = "+"),
                    dropped = paste(attr(x, "dropped"), collapse = "+"), c = kin$c, min_eig = kin$min_eig,
                    tau = f$tau, s2 = f$s2, h2 = f$h2, ci_lo = NA, ci_hi = NA, stat = NA, p = NA, p_davies = NA,
                    boot_valid = NA, boot_error = NA, boot_nonconv = NA, boot_sep_warn = NA, ci_max_gap = NA,
                    reportable = TRUE, note = "", stringsAsFactors = FALSE)
  if (u$type == "binary") {
    sb <- score_boot(y, X, kin$K, B = B, seed = seed_of(u$id, code))
    out[c("stat", "p", "p_davies", "boot_valid", "boot_error", "boot_nonconv", "boot_sep_warn", "reportable")] <-
      list(sb$score, sb$p_boot, sb$p_davies, sb$B_valid, sb$n_error, sb$n_nonconv, sb$n_sep_warn, sb$reportable)
    out$note <- paste(c(sb$note, if (isTRUE(sb$obs_sep_warn)) "separation warning in the observed null model"), collapse = "; ")
    out$note <- sub("^; ", "", out$note)
  } else {
    out$stat <- f$rlrt
    out$p <- p_from_null(f$rlrt, rlrt_null(X, kin, nsim = NSIM, seed = seed_of(u$id, code)))
    if (u$type == "continuous") {
      ci <- profile_ci(y, X, kin, f)
      out[c("ci_lo", "ci_hi", "ci_max_gap")] <- list(ci$lower, ci$upper, ci$max_gap)
    }
  }
  out
}

# Main analysis
main <- do.call(rbind, lapply(units, function(u) analyze(u, unit_data(u, d), "main", 0L)))
grp <- vapply(units, `[[`, "", "role")[main$unit]
main$q <- NA
for (g in c("primary", "exploratory")) main$q[grp == g] <- p.adjust(main$p[grp == g], method = "BH")

# Sensitivity analyses of the primary units
full_filt <- list(A1 = quote(cm == 1 & pop_csf %in% "1"), A2 = quote(cm == 1 & pop_csf %in% "1"),
                  A3 = quote(cm == 1 & pop_csf %in% "1"), B = quote(cm == 1 & pop_op %in% "1"))
blom <- function(v) qnorm((rank(v) - 3 / 8) / (length(v) + 1 / 4))
sens <- list()
loo <- list()
for (uid in c("A1", "A2", "A3", "B")) {
  u <- units[[uid]]
  x <- unit_data(u, d)
  sens[[length(sens) + 1]] <- analyze(u, with_covs(x, setdiff(attr(x, "covs"), "ln_cd4")), "S1b_same_rows_no_cd4", 1L)
  uf <- u
  uf$filt <- full_filt[[uid]]
  uf$covs <- setdiff(u$covs, "ln_cd4")
  sens[[length(sens) + 1]] <- analyze(uf, unit_data(uf, d), "S1c_all_rows_no_cd4", 2L)
  if (u$type == "continuous") {
    r <- residuals(lm(x[[u$y]] ~ unit_X(x) - 1))
    sens[[length(sens) + 1]] <- analyze(u, x, "S2_rank_inverse_normal", 3L, y = blom(r))
  }
  keep <- !(x$sample %in% c("Cn101", "Cn103"))
  sens[[length(sens) + 1]] <- analyze(u, with_covs(x[keep, ], attr(x, "covs")), sprintf("S3a_without_Cn101_Cn103 (%d removed)", sum(!keep)), 4L)
  h2_loo <- vapply(seq_len(nrow(x)), function(i) {
    xi <- with_covs(x[-i, ], attr(x, "covs"))
    fit_lmm(xi[[u$y]], unit_X(xi), unit_kin(K, xi))$h2
  }, 0)
  h2_main <- main$h2[main$unit == uid]
  loo[[uid]] <- data.frame(unit = uid, sample = x$sample, h2_without = h2_loo, change = h2_loo - h2_main)
}
sens <- do.call(rbind, sens)
loo_df <- do.call(rbind, loo)

write_tsv(main, paste0(prefix, "_main.tsv"))
write_tsv(sens, paste0(prefix, "_sensitivity.tsv"))
write_tsv(loo_df, paste0(prefix, "_loo.tsv"))
cat(sprintf("a1_run done: %s_{main,sensitivity,loo}.tsv\n", prefix))
