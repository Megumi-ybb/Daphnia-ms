#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(pomp)
library(dplyr)
library(ggplot2)

## ---------------------------------------------------------------------------
## 1. Settings
## ---------------------------------------------------------------------------
level <- 0.95

## True parameter value (on log scale, since mcap is done on log(ri))
true_params <- readRDS("simulated_data/true_params.rds")
true_log_ri <- log(true_params["ri"])

cat("True log(ri) =", true_log_ri, "\n")
cat("True ri      =", true_params["ri"], "\n\n")

## ---------------------------------------------------------------------------
## 2. Auto-discover all completed profile results
## ---------------------------------------------------------------------------

profile_files <- list.files("coverage_results", pattern = "^profile_ri_\\d+\\.rds$",
                            full.names = TRUE)
cat("Found", length(profile_files), "completed profile files\n\n")

if (length(profile_files) == 0) {
  stop("No profile result files found in coverage_results/")
}

## Extract dataset indices from filenames
dataset_indices <- as.integer(gsub(".*profile_ri_(\\d+)\\.rds$", "\\1", profile_files))

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
  
  subset_data_ri <- readRDS(fname)
  
  # Compute MCAP CI (same as Target_profiling_plots.R)
  mcap_obj <- tryCatch(
    mcap(subset_data_ri$loglik, subset_data_ri$log_ri,
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
  covers <- (true_log_ri >= ci_lo) & (true_log_ri <= ci_hi)
  
  results <- rbind(results, data.frame(
    b      = b,
    ci_lo  = ci_lo,
    ci_hi  = ci_hi,
    mle    = mcap_obj$mle,
    covers = covers
  ))
}

## ---------------------------------------------------------------------------
## 3. Compute coverage
## ---------------------------------------------------------------------------

n_completed <- nrow(results)
n_covering  <- sum(results$covers)
coverage    <- n_covering / n_completed

# Wilson score confidence interval for the coverage proportion
wilson_ci <- function(p, n, z = 1.96) {
  denom <- 1 + z^2/n
  center <- (p + z^2/(2*n)) / denom
  margin <- z * sqrt((p*(1-p)/n + z^2/(4*n^2))) / denom
  c(lower = center - margin, upper = center + margin)
}
cov_ci <- wilson_ci(coverage, n_completed)

## ---------------------------------------------------------------------------
## 4. Print report
## ---------------------------------------------------------------------------

report <- paste0(
  "=== MCAP Coverage Simulation Study: ri ===\n",
  "Nominal level: ", level*100, "%\n",
  "Profile files found: ", length(profile_files), "\n",
  "Successfully processed: ", n_completed, "\n",
  "mcap() failures: ", length(failed), "\n",
  if (length(failed) > 0) paste0("  Failed datasets: ", paste(failed, collapse=", "), "\n") else "",
  "\n",
  "--- Results ---\n",
  "CIs covering true log(ri): ", n_covering, " / ", n_completed, "\n",
  "Empirical coverage: ", sprintf("%.1f%%", coverage*100), "\n",
  # "95% Wilson CI for coverage: [", sprintf("%.1f%%", cov_ci[1]*100), ", ",
  # sprintf("%.1f%%", cov_ci[2]*100), "]\n",
  # "\n",
  "True log(ri) = ", sprintf("%.4f", true_log_ri), "\n",
  "Mean MCAP MLE = ", sprintf("%.4f", mean(results$mle)), "\n",
  "Mean CI width = ", sprintf("%.4f", mean(results$ci_hi - results$ci_lo)), "\n"
)

cat(report)
writeLines(report, "coverage_results/coverage_report.txt")

## ---------------------------------------------------------------------------
## 5. Save full results
## ---------------------------------------------------------------------------

saveRDS(results, file = "coverage_results/coverage_summary.rds")

## ---------------------------------------------------------------------------
## 6. Coverage plot
## ---------------------------------------------------------------------------

results$b_ordered <- rank(results$mle)

p <- ggplot(results, aes(x = b_ordered)) +
  geom_segment(aes(xend = b_ordered, y = ci_lo, yend = ci_hi,
                   color = covers), linewidth = 0.8) +
  geom_point(aes(y = mle), size = 1.2) +
  geom_hline(yintercept = true_log_ri, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "red"),
                     labels = c("TRUE" = "Covers", "FALSE" = "Misses"),
                     name   = "") +
  labs(
    x     = "Simulated dataset (ordered by MLE)",
    y     = expression(log(r[i])),
    title = sprintf("MCAP %d%% CI Coverage for ri  (%.0f/%d = %.1f%%)",
                    level*100, n_covering, n_completed, coverage*100)
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("coverage_results/coverage_plot.pdf", p, width = 8, height = 5)

cat("\nPlot saved to coverage_results/coverage_plot.pdf\n")
cat("Summary saved to coverage_results/coverage_summary.rds\n")
cat("Report saved to coverage_results/coverage_report.txt\n")
