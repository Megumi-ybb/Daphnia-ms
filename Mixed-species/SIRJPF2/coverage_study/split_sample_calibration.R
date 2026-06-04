#!/usr/bin/env Rscript
## Out-of-sample check of the cutoff calibration: does a cutoff calibrated on
## one half of the simulated datasets achieve ~95% coverage on the OTHER half?
## If yes, the calibration generalizes (not just refits its own sims). If the
## out-of-sample coverage falls below 95%, the in-sample 95% is optimistic at B=100.
suppressMessages({library(pomp); library(dplyr)})

true_params <- readRDS("simulated_data/true_params.rds")
params <- c("rn", "sigP", "sigF")

## half-deviance of the smoothed profile at the truth, per simulated dataset
compute_Db <- function(param) {
  log_true <- log(as.numeric(true_params[param])); lc <- paste0("log_", param)
  files <- list.files("coverage_results",
                      pattern = sprintf("^profile_%s_[0-9]+\\.rds$", param), full.names = TRUE)
  Db <- vapply(files, function(f) {
    d <- readRDS(f)
    m <- tryCatch(suppressWarnings(mcap(d$loglik, d[[lc]], level = 0.95, span = 0.95, Ngrid = 1000)),
                  error = function(e) NULL)
    if (is.null(m)) return(NA_real_)
    s <- approx(m$fit$parameter, m$fit$smoothed, xout = log_true)$y
    if (is.na(s)) Inf else max(m$fit$smoothed) - s
  }, numeric(1))
  Db[!is.na(Db)]
}

set.seed(1)
nrep <- 500
cat("param  n  in_sample  out_of_sample (mean, SD over", nrep, "random 50/50 splits)\n")
for (p in params) {
  Db <- compute_Db(p); n <- length(Db)
  insample <- mean(Db <= quantile(Db, 0.95, type = 7))
  oos <- replicate(nrep, {
    idx   <- sample(n, n %/% 2)
    cstar <- quantile(Db[idx], 0.95, type = 7)   # calibrate on one half
    mean(Db[-idx] <= cstar)                       # coverage on the held-out half
  })
  cat(sprintf("%-5s  %d  %.3f      %.3f  (SD %.3f)\n", p, n, insample, mean(oos), sd(oos)))
}
