#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(pomp)
library(dplyr)
library(ggplot2)

## ---------------------------------------------------------------------------
## 1. Settings
## ---------------------------------------------------------------------------
level <- 0.95

## True parameter values
true_params <- readRDS("simulated_data/true_params.rds")

## Parameters to process: auto-detect every single-parameter profile on disk.
## (To restrict, set param_names manually, e.g. c("rn","sigP","sigF").)
profile_files_all <- list.files("coverage_results", pattern = "^profile_[A-Za-z_]+_[0-9]+\\.rds$")
param_names <- sort(unique(sub("^profile_([A-Za-z_]+)_[0-9]+\\.rds$", "\\1", profile_files_all)))
if (length(param_names) == 0) stop("No profile_<param>_<b>.rds files found in coverage_results/.")
cat("Single-parameter profiles detected:", paste(param_names, collapse = ", "), "\n\n")
## Wilson score confidence interval for the coverage proportion
wilson_ci <- function(p, n, z = 1.96) {
  denom <- 1 + z^2/n
  center <- (p + z^2/(2*n)) / denom
  margin <- z * sqrt((p*(1-p)/n + z^2/(4*n^2))) / denom
  c(lower = center - margin, upper = center + margin)
}

## ---------------------------------------------------------------------------
## 2. Loop over parameters
## ---------------------------------------------------------------------------

all_results <- list()

for (param in param_names) {

  true_val     <- true_params[param]
  true_log_val <- log(true_val)
  log_col      <- paste0("log_", param)

  cat("============================================================\n")
  cat("Processing parameter:", param, "\n")
  cat(sprintf("  True %s      = %.4f\n", param, true_val))
  cat(sprintf("  True log(%s) = %.4f\n", param, true_log_val))
  cat("============================================================\n\n")

  ## Discover profile files for this parameter
  pattern <- sprintf("^profile_%s_\\d+\\.rds$", param)
  profile_files <- list.files("coverage_results", pattern = pattern,
                              full.names = TRUE)
  cat("Found", length(profile_files), "completed profile files for", param, "\n\n")

  if (length(profile_files) == 0) {
    cat("  Skipping", param, "— no profile files found.\n\n")
    next
  }

  ## Extract dataset indices from filenames
  idx_pattern <- sprintf(".*profile_%s_(\\d+)\\.rds$", param)
  dataset_indices <- as.integer(gsub(idx_pattern, "\\1", profile_files))

  results <- data.frame(
    b      = integer(0),
    ci_lo  = numeric(0),
    ci_hi  = numeric(0),
    mle    = numeric(0),
    covers = logical(0)
  )

  failed <- c()

  for (k in seq_along(profile_files)) {
    b     <- dataset_indices[k]
    fname <- profile_files[k]

    subset_data <- readRDS(fname)

    ## Compute MCAP CI
    mcap_obj <- tryCatch(
      mcap(subset_data$loglik, subset_data[[log_col]],
           level = level, span = 0.95, Ngrid = 1000),
      error = function(e) {
        cat("  WARNING: mcap failed for dataset", b, ":", e$message, "\n")
        return(NULL)
      }
    )

    if (is.null(mcap_obj)) {
      failed <- c(failed, b)
      next
    }

    ci_lo  <- mcap_obj$ci[1]
    ci_hi  <- mcap_obj$ci[2]
    covers <- (true_log_val >= ci_lo) & (true_log_val <= ci_hi)

    results <- rbind(results, data.frame(
      b      = b,
      ci_lo  = ci_lo,
      ci_hi  = ci_hi,
      mle    = mcap_obj$mle,
      covers = covers
    ))
  }

  ## -------------------------------------------------------------------------
  ## 3. Compute coverage
  ## -------------------------------------------------------------------------

  n_completed <- nrow(results)
  n_covering  <- sum(results$covers)
  coverage    <- n_covering / n_completed

  cov_ci <- wilson_ci(coverage, n_completed)

  ## -------------------------------------------------------------------------
  ## 4. Print report
  ## -------------------------------------------------------------------------

  report <- paste0(
    sprintf("=== MCAP Coverage Simulation Study: %s ===\n", param),
    "Nominal level: ", level*100, "%\n",
    "Profile files found: ", length(profile_files), "\n",
    "Successfully processed: ", n_completed, "\n",
    "mcap() failures: ", length(failed), "\n",
    if (length(failed) > 0) paste0("  Failed datasets: ", paste(failed, collapse=", "), "\n") else "",
    "\n",
    "--- Results ---\n",
    sprintf("CIs covering true log(%s): %d / %d\n", param, n_covering, n_completed),
    sprintf("Empirical coverage: %.1f%%\n", coverage*100),
    sprintf("True log(%s) = %.4f\n", param, true_log_val),
    sprintf("Mean MCAP MLE = %.4f\n", mean(results$mle)),
    sprintf("Mean CI width = %.4f\n", mean(results$ci_hi - results$ci_lo))
  )

  cat(report)
  # writeLines(report, sprintf("coverage_results/coverage_report_%s.txt", param))

  ## -------------------------------------------------------------------------
  ## 5. Save full results
  ## -------------------------------------------------------------------------

  # saveRDS(results, file = sprintf("coverage_results/coverage_summary_%s.rds", param))

  ## -------------------------------------------------------------------------
  ## 6. Coverage plot
  ## -------------------------------------------------------------------------

  results$b_ordered <- rank(results$mle)

  y_label <- paste0("log(", param, ")")

  p <- ggplot(results, aes(x = b_ordered)) +
    geom_segment(aes(xend = b_ordered, y = ci_lo, yend = ci_hi,
                     color = covers), linewidth = 0.8) +
    geom_point(aes(y = mle), size = 1.2) +
    geom_hline(yintercept = true_log_val, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "red"),
                       labels = c("TRUE" = "Covers", "FALSE" = "Misses"),
                       name   = "") +
    labs(
      x     = "Simulated dataset (ordered by MLE)",
      y     = y_label,
      title = sprintf("MCAP %d%% CI Coverage for %s  (%.0f/%d = %.1f%%)",
                      level*100, param, n_covering, n_completed, coverage*100)
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  plot_file <- sprintf("coverage_results/coverage_plot_%s.pdf", param)
  ggsave(plot_file, p, width = 8, height = 5)

  cat(sprintf("\nPlot saved to %s\n", plot_file))
  cat(sprintf("Summary saved to coverage_results/coverage_summary_%s.rds\n", param))
  cat(sprintf("Report saved to coverage_results/coverage_report_%s.txt\n\n", param))

  all_results[[param]] <- results
}

## ---------------------------------------------------------------------------
## 7. Combined summary
## ---------------------------------------------------------------------------

cat("\n========== COMBINED SUMMARY ==========\n")
for (param in names(all_results)) {
  res <- all_results[[param]]
  n   <- nrow(res)
  cov <- sum(res$covers) / n
  cat(sprintf("  %s : %d/%d = %.1f%% coverage\n", param, sum(res$covers), n, cov*100))
}
cat("======================================\n")


## ---------------------------------------------------------------------------
## 8. Composite ridge targets (ri*f_Si, probi*f_Si)  --  METHOD NOT YET CONFIRMED
## ---------------------------------------------------------------------------
## ri and probi are individually non-identifiable; the identified quantities are the
## products ri*f_Si (recruitment) and probi*f_Si (infection), reported in Table 1.
## To assess composite coverage we read the profile of its DRIVING parameter (the
## generalized coverage_profile.R persists ALL coef columns), form the composite at
## each profile point, and run mcap on (loglik, log(composite)).
##
## REQUIRES profiles produced by the generalized coverage_profile.R (the old
## profile_ri_*.rds carry only (ri, loglik, log_ri), no f_Si column -> skipped).
## !! VERIFY before reporting: confirm this profile-derived composite CI matches the
##    construction of the Table-1 composite CIs (data/.../profile_ci_table.rda).
## Disabled by default so an unconfirmed method never enters the summary.
do_composite <- FALSE
composites <- list(
  ri_f_Si    = c("ri",    "f_Si"),
  probi_f_Si = c("probi", "f_Si")
)
composite_coverage <- function(name, factors, driver = factors[1], level = 0.95) {
  true_comp_log <- log(prod(sapply(factors, function(p) as.numeric(true_params[p]))))
  files  <- list.files("coverage_results",
                       pattern = sprintf("^profile_%s_[0-9]+\\.rds$", driver), full.names = TRUE)
  covers <- logical(0)
  for (f in files) {
    d <- readRDS(f)
    if (!all(factors %in% names(d))) next            # need partner coef columns
    log_comp <- log(Reduce(`*`, lapply(factors, function(p) d[[p]])))
    m <- tryCatch(mcap(d$loglik, log_comp, level = level, span = 0.95, Ngrid = 1000),
                  error = function(e) NULL)
    if (is.null(m)) next
    covers <- c(covers, true_comp_log >= m$ci[1] && true_comp_log <= m$ci[2])
  }
  n <- length(covers)
  cat(sprintf("  %-12s : %d/%d = %.1f%% (true log = %.3f)\n",
              name, sum(covers), n, if (n) 100*mean(covers) else NA, true_comp_log))
}
if (do_composite) {
  cat("\n=== Composite ridge coverage (METHOD UNCONFIRMED -- verify vs Table 1) ===\n")
  for (nm in names(composites)) composite_coverage(nm, composites[[nm]], level = level)
}

## ---------------------------------------------------------------------------
## 9. One-sided / boundary parameters
## ---------------------------------------------------------------------------
## For a param whose MLE sits at a bound (e.g. theta_P, sigIi) the two-sided MCAP CI
## can be degenerate; score coverage against the single finite bound instead.
## side="lower": CI is [ci_lo, Inf), covered iff true_log >= ci_lo;
## side="upper": CI is (-Inf, ci_hi], covered iff true_log <= ci_hi.
## Disabled by default; set the sides from the actual boundary before reporting.
do_onesided  <- FALSE
onesided_side <- c(theta_P = "lower", sigIi = "lower")   # EDIT to match each param's bound
onesided_coverage <- function(param, side, level = 0.95) {
  true_log <- log(as.numeric(true_params[param]))
  log_col  <- paste0("log_", param)
  files    <- list.files("coverage_results",
                         pattern = sprintf("^profile_%s_[0-9]+\\.rds$", param), full.names = TRUE)
  covers <- logical(0)
  for (f in files) {
    d <- readRDS(f)
    m <- tryCatch(mcap(d$loglik, d[[log_col]], level = level, span = 0.95, Ngrid = 1000),
                  error = function(e) NULL)
    if (is.null(m)) next
    covers <- c(covers, if (side == "lower") true_log >= m$ci[1] else true_log <= m$ci[2])
  }
  n <- length(covers)
  cat(sprintf("  %-10s (%s-bound): %d/%d = %.1f%%\n",
              param, side, sum(covers), n, if (n) 100*mean(covers) else NA))
}
if (do_onesided) {
  cat("\n=== One-sided coverage (CONFIRM each side) ===\n")
  for (p in names(onesided_side)) onesided_coverage(p, onesided_side[[p]], level = level)
}
