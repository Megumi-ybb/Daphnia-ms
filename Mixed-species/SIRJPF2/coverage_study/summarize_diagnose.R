#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(dplyr)
library(tidyr)
library(ggplot2)

## ---------------------------------------------------------------------------
## Summarize results from diagnose_bias.R across datasets.
##
## Two questions:
##   (Q1) Does the BEST mif2 fit (warm-started at truth) end up near truth,
##        or systematically off?  -> coef drift on log scale.
##   (Q2) Does the best mif2 fit reach a HIGHER loglik than pfilter(truth)?
##        If yes,  the simulated dataset's true MLE is itself off-truth
##                 (genuine finite-sample bias).
##        If no,   mif2 stayed at truth (coverage misses are caused by the
##                 random-start global search not reaching truth).
## ---------------------------------------------------------------------------

results_dir <- "diagnose_results"
focus_params <- c("ri","rn","sigP","sigF","f_Si","f_Sn","probi","probn",
                  "sigIn","sigIi","sigJi","sigJn")

files <- list.files(results_dir, pattern = "^diagnose_\\d+\\.rds$",
                    full.names = TRUE)

if (length(files) == 0) {
  stop("No diagnose_*.rds files found in ", results_dir)
}

cat("Found", length(files), "diagnostic dataset(s):\n  ",
    paste(basename(files), collapse = ", "), "\n\n")

## ---------------------------------------------------------------------------
## 1. Per-dataset best fit and loglik comparison
## ---------------------------------------------------------------------------

best_rows <- list()
ll_compare <- list()

for (f in files) {
  res <- readRDS(f)
  b   <- res$b
  best_idx <- which.max(res$fits$loglik)

  best_rows[[as.character(b)]] <- res$fits[best_idx, , drop = FALSE]
  ll_compare[[as.character(b)]] <- data.frame(
    b            = b,
    ll_truth     = res$ll_truth[1],
    ll_truth_se  = res$ll_truth[2],
    ll_best      = res$fits$loglik[best_idx],
    ll_best_se   = res$fits$loglik_se[best_idx],
    ll_gain      = res$fits$loglik[best_idx] - res$ll_truth[1]
  )
}

best_df <- bind_rows(lapply(seq_along(best_rows), function(i) {
  row <- best_rows[[i]]
  row$b <- as.integer(names(best_rows)[i])
  row
}))

ll_df <- bind_rows(ll_compare)

## ---------------------------------------------------------------------------
## 2. Q2 — loglik gain from truth -> mif2-best
## ---------------------------------------------------------------------------

cat("=== Q2: ll(mif2-best) vs ll(truth) per dataset ===\n")
cat("Positive gain = mif2 found a higher-likelihood point than truth.\n")
cat("If gain ~ 0 (within ~2 SE), mif2 stayed at truth.\n\n")

ll_print <- ll_df %>%
  arrange(b) %>%
  mutate(
    pooled_se = sqrt(ll_truth_se^2 + ll_best_se^2),
    z         = ll_gain / pooled_se,
    flag      = ifelse(abs(z) < 2, "stayed",
                ifelse(z >= 2, "drift+", "lower!"))
  )

print(ll_print %>%
        select(b, ll_truth, ll_best, ll_gain, pooled_se, z, flag),
      row.names = FALSE)

cat(sprintf("\nMean ll_gain across datasets: %+.3f  (SD %.3f)\n",
            mean(ll_df$ll_gain), sd(ll_df$ll_gain)))

## ---------------------------------------------------------------------------
## 3. Q1 — per-parameter drift on log scale (best mif2 vs truth)
## ---------------------------------------------------------------------------

cat("\n=== Q1: log-scale drift of best mif2 fit from truth ===\n")
cat("drift = log(best_param) - log(true_param)\n")
cat("Mean ~ 0  -> mif2 returns to truth.\n")
cat("Mean << 0 -> systematic downward bias (mif2 prefers smaller values).\n")
cat("Mean >> 0 -> systematic upward bias.\n\n")

true_params <- readRDS(files[1]) |>
  (\(r) r$true_params)()

drift_summary <- data.frame(
  param      = focus_params,
  truth      = sapply(focus_params, function(p) as.numeric(true_params[p])),
  log_truth  = sapply(focus_params, function(p) log(as.numeric(true_params[p]))),
  mean_drift = NA_real_,
  sd_drift   = NA_real_,
  min_drift  = NA_real_,
  max_drift  = NA_real_
)

for (i in seq_along(focus_params)) {
  p <- focus_params[i]
  drift <- log(best_df[[p]]) - log(as.numeric(true_params[p]))
  drift_summary$mean_drift[i] <- mean(drift)
  drift_summary$sd_drift[i]   <- sd(drift)
  drift_summary$min_drift[i]  <- min(drift)
  drift_summary$max_drift[i]  <- max(drift)
}

print(drift_summary, row.names = FALSE, digits = 3)

## ---------------------------------------------------------------------------
## 4. Within-dataset spread across the n_reps mif2 chains (consistency check)
## ---------------------------------------------------------------------------

cat("\n=== Within-dataset chain spread (log scale) ===\n")
cat("Big spread -> mif2 chains end up in different places on the same dataset\n")
cat("             (algorithm noisy, even from a single starting point).\n\n")

within_spread <- data.frame(param = focus_params, mean_within_sd = NA_real_)
for (i in seq_along(focus_params)) {
  p <- focus_params[i]
  sds <- sapply(files, function(f) {
    res <- readRDS(f)
    sd(log(res$fits[[p]]))
  })
  within_spread$mean_within_sd[i] <- mean(sds)
}
print(within_spread, row.names = FALSE, digits = 3)

## ---------------------------------------------------------------------------
## 5. Save and plot
## ---------------------------------------------------------------------------

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(ll_df = ll_df, best_df = best_df,
             drift_summary = drift_summary,
             within_spread = within_spread),
        file = file.path(results_dir, "diagnose_summary.rds"))

# Plot: per-parameter drift, one boxplot per parameter
plot_df <- bind_rows(lapply(focus_params, function(p) {
  data.frame(
    param = p,
    drift = log(best_df[[p]]) - log(as.numeric(true_params[p]))
  )
}))

p_drift <- ggplot(plot_df, aes(x = param, y = drift)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  labs(
    x = "parameter",
    y = "log(best mif2 fit) - log(truth)",
    title = "Drift of best mif2 fit from truth, warm-started at truth",
    subtitle = sprintf("%d datasets, %d mif2 chains each", nrow(best_df),
                       nrow(readRDS(files[1])$fits))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(results_dir, "diagnose_drift.pdf"), p_drift,
       width = 8, height = 5)
cat(sprintf("\nDrift plot saved to %s/diagnose_drift.pdf\n", results_dir))
cat(sprintf("Summary saved to %s/diagnose_summary.rds\n", results_dir))
