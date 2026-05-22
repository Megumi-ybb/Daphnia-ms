#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
# Aggregate paired (null, alt) fits across all completed (truth, b) combinations,
# tabulate AIC-selection rates. Output to ../../../daphnia-article/data/.

B <- 50
truths <- c("shared", "unit_specific")

read_fit <- function(model, truth, b) {
  f <- sprintf("results/fit_%s_%s_%d.rds", model, truth, b)
  if (!file.exists(f)) return(NULL)
  readRDS(f)
}

rows <- list()
for (truth in truths) {
  for (b in 1:B) {
    null_res <- read_fit("null", truth, b)
    alt_res  <- read_fit("alt",  truth, b)
    if (is.null(null_res) || is.null(alt_res)) next
    rows[[length(rows) + 1]] <- data.frame(
      truth     = truth,
      b         = b,
      ll_null   = null_res$ll,
      AIC_null  = null_res$AIC,
      k_null    = null_res$k,
      ll_alt    = alt_res$ll,
      AIC_alt   = alt_res$AIC,
      k_alt     = alt_res$k,
      Delta_AIC = alt_res$AIC - null_res$AIC,  # negative => alt selected
      selected  = ifelse(alt_res$AIC < null_res$AIC, "alt", "null"),
      stringsAsFactors = FALSE
    )
  }
}

if (length(rows) == 0) {
  stop("No paired (null, alt) fits found. Run fit_null.R and fit_alt.R first.")
}

paired <- do.call(rbind, rows)
cat(sprintf("Paired fits collected: %d (truth=shared: %d, truth=unit_specific: %d)\n",
            nrow(paired),
            sum(paired$truth == "shared"),
            sum(paired$truth == "unit_specific")))

# ---- Per-truth selection-rate summary ----
summarize_truth <- function(sub) {
  data.frame(
    truth          = sub$truth[1],
    n              = nrow(sub),
    selects_null   = sum(sub$selected == "null"),
    selects_alt    = sum(sub$selected == "alt"),
    pct_select_alt = mean(sub$selected == "alt") * 100,
    mean_Delta_AIC = mean(sub$Delta_AIC),
    median_Delta_AIC = median(sub$Delta_AIC),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(split(paired, paired$truth), summarize_truth))
rownames(summary_df) <- NULL

cat("\n=== AIC selection-rate summary ===\n")
print(summary_df, row.names = FALSE)

# ---- Save ----
proj_root <- normalizePath("../../..")
out_dir <- file.path(proj_root, "daphnia-article", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(list(paired = paired, summary = summary_df),
        file.path(out_dir, "aic_power_summary.rds"))
write.csv(paired,     file.path(out_dir, "aic_power_paired.csv"),  row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "aic_power_summary.csv"), row.names = FALSE)

# ---- Human-readable text ----
txt <- character()
push <- function(...) txt <<- c(txt, sprintf(...))

push("AIC selection-rate simulation for SIRJPF2 (Task #2)")
push("===================================================")
push("")
push("Design: parametric simulation following Stocks 2020 Table D3.")
push("  Truth A (shared):        theta_In is the same across 8 units (all-shared MLE).")
push("  Truth B (unit_specific): theta_In varies +/- 0.5 log units across units")
push("                           (factor 2.7 spread, comparable to Dent_SIRJPF/theta_In)")
push("  For each simulated panel, fit both models and select the lower-AIC one.")
push("")
push("Per-truth selection rates:")
push("")
push("%-15s %4s %12s %11s %14s %15s",
     "truth", "n", "selects_null", "selects_alt", "pct_select_alt", "mean_Delta_AIC")
for (i in seq_len(nrow(summary_df))) {
  push("%-15s %4d %12d %11d %13.1f%% %15.2f",
       summary_df$truth[i],
       summary_df$n[i],
       summary_df$selects_null[i],
       summary_df$selects_alt[i],
       summary_df$pct_select_alt[i],
       summary_df$mean_Delta_AIC[i])
}
push("")
push("Delta_AIC = AIC_alt - AIC_null; negative => alt (unit-specific) selected.")
push("")
push("Interpretation guide:")
push("  - Under truth=shared:        AIC should select null on most panels (low false-positive).")
push("  - Under truth=unit_specific: AIC should select alt on most panels (high power).")

writeLines(txt, file.path(out_dir, "aic_power_summary.txt"))

cat(sprintf("\nSaved:\n  %s\n  %s\n  %s\n  %s\n",
            file.path(out_dir, "aic_power_summary.rds"),
            file.path(out_dir, "aic_power_summary.csv"),
            file.path(out_dir, "aic_power_paired.csv"),
            file.path(out_dir, "aic_power_summary.txt")))
