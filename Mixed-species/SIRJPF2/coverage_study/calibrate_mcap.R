#!/usr/bin/env Rscript
## Option-B sensitivity check: use the existing parametric-bootstrap profiles to
## CALIBRATE the MCAP cutoff so the interval attains nominal 95% coverage.
## Local, no cluster: re-reads coverage_results/profile_<param>_<b>.rds.
##
## For each simulated dataset b we form the smoothed-profile half-deviance at the
## TRUE value, D_b = max(smoothed_ll) - smoothed_ll(true).  The interval covers iff
## D_b <= cutoff.  The asymptotic cutoff is qchisq(.95,1)/2 = 1.92; the CALIBRATED
## cutoff cstar is the 95th percentile of {D_b}.
##
## CAVEATS (verified by rigor audit 2026-06-03):
##  - calibrated_cov is IN-SAMPLE / tautological (mean(D_b <= quantile(D_b,.95)) ~ .95
##    for any sample); it is NOT out-of-sample validation. The load-bearing numbers are
##    cstar (the calibrated cutoff) and the empirical CI-width ratio.
##  - nominal_cov uses each dataset's OWN mcap delta (per-dataset MC-adjusted cutoff);
##    the D_b<=1.92 fixed-rule coverage can differ by ~1 pt. Both reported below.
##  - width ratio is computed EMPIRICALLY from the smoothed profile (no quadratic
##    assumption); a quadratic-approx factor sqrt(cstar/1.92) is also shown for reference.
suppressMessages({library(pomp); library(dplyr)})

tp     <- readRDS("simulated_data/true_params.rds")
params <- c("rn", "sigP", "sigF")          # ri excluded: non-identifiable (flat profile)
nominal_cutoff <- qchisq(0.95, 1) / 2      # 1.9207

## width of the {smoothed >= smax - cutoff} region of a profile, on the log scale
ci_width <- function(par, smoothed, cutoff) {
  keep <- smoothed >= (max(smoothed) - cutoff)
  if (!any(keep)) return(NA_real_)
  diff(range(par[keep]))
}

rows <- list()
for (p in params) {
  log_true <- log(as.numeric(tp[p]))
  lc <- paste0("log_", p)
  files <- list.files("coverage_results",
                      pattern = sprintf("^profile_%s_\\d+\\.rds$", p),
                      full.names = TRUE)
  Db <- numeric(0); cov_nom <- integer(0); fits <- list()
  for (f in files) {
    d <- readRDS(f)
    m <- tryCatch(suppressWarnings(mcap(d$loglik, d[[lc]], level = 0.95, span = 0.95, Ngrid = 1000)),
                  error = function(e) NULL)
    if (is.null(m)) next
    smax      <- max(m$fit$smoothed)
    s_at_true <- approx(m$fit$parameter, m$fit$smoothed, xout = log_true)$y  # NA if truth beyond grid
    Db        <- c(Db, if (is.na(s_at_true)) Inf else smax - s_at_true)
    cov_nom   <- c(cov_nom, as.integer(log_true >= m$ci[1] && log_true <= m$ci[2]))
    fits[[length(fits) + 1]] <- list(par = m$fit$parameter, sm = m$fit$smoothed)
  }
  n     <- length(Db)
  cstar <- as.numeric(quantile(Db, 0.95, type = 7))
  ## empirical CI-width ratio (calibrated cutoff vs fixed 1.92), median over datasets
  wr <- vapply(fits, function(g) ci_width(g$par, g$sm, cstar) / ci_width(g$par, g$sm, nominal_cutoff),
               numeric(1))
  rows[[p]] <- data.frame(
    param              = p,
    n                  = n,
    nominal_cov        = mean(cov_nom),                 # mcap per-dataset-delta CI membership
    cov_at_192         = mean(Db <= nominal_cutoff),    # fixed 1.92 half-deviance rule
    calibrated_cutoff  = cstar,
    calibrated_cov     = mean(Db <= cstar),             # IN-SAMPLE / by construction ~0.95
    width_ratio_emp    = median(wr, na.rm = TRUE),      # empirical, no quadratic assumption
    width_ratio_quad   = sqrt(cstar / nominal_cutoff),  # quadratic-approx reference
    truth_beyond_grid  = mean(!is.finite(Db))
  )
}
res <- bind_rows(rows)
cat("=== MCAP cutoff calibration (option B sensitivity) ===\n")
cat("calibrated_cov is IN-SAMPLE (tautological); load-bearing numbers are cutoff + width_ratio_emp.\n\n")
print(res, row.names = FALSE, digits = 3)
saveRDS(res, "mcap_calibration.rds")
cat("\nSaved to mcap_calibration.rds\n")
