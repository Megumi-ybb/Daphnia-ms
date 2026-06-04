#!/usr/bin/env Rscript
## Ed's option-B table: for each parameter,
##   (1) the basic MCAP 95% interval (real data),
##   (2) its empirical coverage from the parametric-bootstrap simulation, and
##   (3) the cutoff-corrected interval that targets 95% coverage.
## Computed entirely from SAVED results:
##   - simulation profiles : coverage_results/profile_<param>_<b>.rds
##   - real-data profiles  : ../profile/<param>/<param>.RData  (object final_params)
## The corrected cutoff cstar is the 95th percentile of the simulated half-deviance
## at the truth, D_b = max(smoothed_ll) - smoothed_ll(true); re-cutting the real-data
## profile at cstar gives ~95% coverage by construction.
suppressMessages({library(pomp); library(dplyr)})

true_params    <- readRDS("simulated_data/true_params.rds")
nominal_cutoff <- qchisq(0.95, 1) / 2          # 1.92

## endpoints of {smoothed >= max(smoothed) - cutoff} on a profile grid
ci_from_smoothed <- function(par, sm, cutoff) {
  keep <- sm >= (max(sm) - cutoff)
  if (!any(keep)) return(c(NA_real_, NA_real_))
  range(par[keep])
}

## calibrated cutoff + basic (mcap) coverage from the SIMULATION profiles
sim_calibration <- function(param) {
  log_true <- log(as.numeric(true_params[param]))
  lc <- paste0("log_", param)
  files <- list.files("coverage_results",
                      pattern = sprintf("^profile_%s_[0-9]+\\.rds$", param), full.names = TRUE)
  Db <- numeric(0); cov <- integer(0)
  for (f in files) {
    d <- readRDS(f)
    m <- tryCatch(suppressWarnings(mcap(d$loglik, d[[lc]], level = 0.95, span = 0.95, Ngrid = 1000)),
                  error = function(e) NULL)
    if (is.null(m)) next
    s_true <- approx(m$fit$parameter, m$fit$smoothed, xout = log_true)$y
    Db  <- c(Db, if (is.na(s_true)) Inf else max(m$fit$smoothed) - s_true)
    cov <- c(cov, as.integer(log_true >= m$ci[1] && log_true <= m$ci[2]))
  }
  list(cstar = as.numeric(quantile(Db, 0.95, type = 7)),
       coverage = mean(cov), corrected_cov = mean(Db <= as.numeric(quantile(Db, 0.95, type = 7))),
       n = length(Db))
}

## basic mcap interval + corrected (cstar) interval from the REAL-DATA profile
realdata_intervals <- function(param, cstar) {
  e <- new.env(); load(sprintf("../profile/%s/%s.RData", param, param), e)
  sub   <- e$final_params %>% group_by(.data[[param]]) %>% filter(loglik == max(loglik)) %>% ungroup()
  m     <- suppressWarnings(mcap(sub$loglik, log(sub[[param]]), level = 0.95, span = 0.95, Ngrid = 1000))
  list(basic = as.numeric(m$ci),
       corrected = ci_from_smoothed(m$fit$parameter, m$fit$smoothed, cstar),
       mle = m$mle)
}

## auto-detect: every param with BOTH a simulation profile and a real-data profile.
## (After the #3 cluster run this produces the full table automatically.)
sim_params <- unique(sub("^profile_([A-Za-z_]+)_[0-9]+\\.rds$", "\\1",
                         list.files("coverage_results", pattern = "^profile_[A-Za-z_]+_[0-9]+\\.rds$")))
params <- sort(sim_params[sapply(sim_params, function(p) file.exists(sprintf("../profile/%s/%s.RData", p, p)))])
params <- setdiff(params, "ri")       # ri non-identifiable (degenerate, ~95.8%); calibration N/A
cat("Parameters with both simulation + real-data profiles:", paste(params, collapse = ", "), "\n\n")
rows <- list()
for (p in params) {
  cal <- sim_calibration(p)
  iv  <- realdata_intervals(p, cal$cstar)
  rows[[p]] <- data.frame(
    param            = p,
    MLE_log          = iv$mle,
    basic_lo_log     = iv$basic[1],   basic_hi_log = iv$basic[2],
    basic_coverage   = cal$coverage,
    cutoff_basic     = nominal_cutoff, cutoff_corrected = cal$cstar,
    corrected_lo_log = iv$corrected[1], corrected_hi_log = iv$corrected[2],
    corrected_coverage = cal$corrected_cov,
    width_ratio      = (iv$corrected[2] - iv$corrected[1]) / (iv$basic[2] - iv$basic[1]),
    n_sim            = cal$n
  )
}
res <- bind_rows(rows)

cat("=== Option-B coverage-corrected interval table (LOG scale) ===\n")
print(res[, c("param","basic_lo_log","basic_hi_log","basic_coverage",
              "cutoff_corrected","corrected_lo_log","corrected_hi_log",
              "corrected_coverage","width_ratio")], row.names = FALSE, digits = 4)

nat <- transform(res,
  basic_lo = exp(basic_lo_log), basic_hi = exp(basic_hi_log),
  corr_lo  = exp(corrected_lo_log), corr_hi = exp(corrected_hi_log))
cat("\n=== natural scale (exp of log endpoints) ===\n")
print(nat[, c("param","basic_lo","basic_hi","basic_coverage","corr_lo","corr_hi","width_ratio")],
      row.names = FALSE, digits = 4)

saveRDS(res, "corrected_ci_table.rds")
cat("\nSaved corrected_ci_table.rds\n")
