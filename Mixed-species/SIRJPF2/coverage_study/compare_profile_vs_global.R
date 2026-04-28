#!/usr/bin/env Rscript
library(dplyr)

## ---------------------------------------------------------------------------
## For each (dataset b, profile param p), compare three loglik values:
##   ll_truth         : pfilter at theta_true on dataset b  (diagnose_<b>.rds)
##   ll_profile_max   : best loglik found by the coverage profile for p on b
##                      (coverage_results/profile_<p>_<b>.rds, max of loglik)
##   ll_global_best   : best loglik from unconstrained mif2 warm-started at
##                      truth on b (diagnose_<b>.rds, max of fits$loglik)
##
## Question: is the coverage profile under-converged?
##   - profile_max ~ global_best   -> profile is converged; gap from truth is real
##   - profile_max < global_best   -> profile mif2 is leaving ll on the table
##   - profile_max < truth         -> profile severely under-converged
## ---------------------------------------------------------------------------

profile_params <- c("ri","rn","sigP","sigF")

diag_files <- list.files("diagnose_results", pattern = "^diagnose_\\d+\\.rds$",
                         full.names = TRUE)
b_indices  <- as.integer(gsub(".*diagnose_(\\d+)\\.rds$", "\\1", diag_files))

cat("Comparing profile-max loglik vs global-best loglik vs ll(truth)\n")
cat("on", length(b_indices), "datasets x", length(profile_params), "profile params.\n\n")

rows <- list()

for (i in seq_along(diag_files)) {
  b <- b_indices[i]
  diag <- readRDS(diag_files[i])

  ll_truth       <- diag$ll_truth[1]
  ll_truth_se    <- diag$ll_truth[2]
  ll_global_best <- max(diag$fits$loglik)

  for (p in profile_params) {
    f_prof <- sprintf("coverage_results/profile_%s_%03d.rds", p, b)
    if (!file.exists(f_prof)) {
      cat(sprintf("  [skip] %s missing for dataset %d\n", p, b))
      next
    }
    sub <- readRDS(f_prof)
    ll_profile_max <- max(sub$loglik)

    rows[[length(rows) + 1]] <- data.frame(
      b               = b,
      param           = p,
      ll_truth        = ll_truth,
      ll_truth_se     = ll_truth_se,
      ll_profile_max  = ll_profile_max,
      ll_global_best  = ll_global_best,
      gap_prof_truth  = ll_profile_max - ll_truth,
      gap_glob_truth  = ll_global_best - ll_truth,
      gap_glob_prof   = ll_global_best - ll_profile_max
    )
  }
}

tab <- bind_rows(rows) %>% arrange(param, b)

cat("=== Loglik comparison ===\n")
cat("gap_prof_truth = ll_profile_max - ll_truth   (positive = profile beats truth)\n")
cat("gap_glob_truth = ll_global_best - ll_truth   (the diagnostic gain)\n")
cat("gap_glob_prof  = ll_global_best - ll_profile_max  (positive = profile leaves ll on the table)\n\n")

print(tab %>%
        select(b, param, ll_truth, ll_profile_max, ll_global_best,
               gap_prof_truth, gap_glob_truth, gap_glob_prof),
      row.names = FALSE, digits = 4)

cat("\n=== Per-parameter means ===\n")
summ <- tab %>%
  group_by(param) %>%
  summarise(
    mean_gap_prof_truth = mean(gap_prof_truth),
    mean_gap_glob_truth = mean(gap_glob_truth),
    mean_gap_glob_prof  = mean(gap_glob_prof),
    n = n(),
    .groups = "drop"
  )
print(summ, row.names = FALSE, digits = 3)

cat("\nInterpretation key:\n")
cat("  gap_glob_prof ~ 0   -> profile mif2 is converged; coverage failures are not a compute problem\n")
cat("  gap_glob_prof > ~2  -> profile mif2 is under-converged; more compute would help\n")
cat("  gap_prof_truth < 0  -> profile didn't even reach truth's ll; severely under-converged\n")
