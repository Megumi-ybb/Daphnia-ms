#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
# Collect bootstrap LRT results and compute p-values
# Run locally after all HPC jobs complete

library(ggplot2)

B <- 100

alt_names <- c("xi", "theta_Sn", "theta_In", "theta_P", "probn", "rn", "f_Sn")

# ---- Load observed LRT statistics from parameter table ----
load("../../../../data/Simple_dynamics/Dent/para/Dent_para_loglik_df.rds")
null_ll_obs <- dent_para_parameter_table["all_shared", "ll"]
null_AIC_obs <- dent_para_parameter_table["all_shared", "AIC"]

cat("Observed null log-likelihood:", null_ll_obs, "\n\n")

# ---- Load bootstrap null results ----
null_lls <- numeric(B)
null_completed <- 0
for (b in 1:B) {
  f <- sprintf("results_null/lrt_null_%d.rds", b)
  if (file.exists(f)) {
    res <- readRDS(f)
    null_lls[b] <- res$ll
    null_completed <- null_completed + 1
  } else {
    null_lls[b] <- NA
  }
}
cat("Null fits completed:", null_completed, "/", B, "\n")

# ---- For each alternative, compute bootstrap p-value ----
results <- data.frame(
  alt_name = character(),
  Lambda_obs = numeric(),
  p_chisq = numeric(),
  p_boot = numeric(),
  B_completed = integer(),
  stringsAsFactors = FALSE
)

# Row name mapping: alt_name directly matches row names in parameter table
# Row names are: all_shared, xi, theta_Sn, theta_P, theta_In, rn, probn, f_Sn
# (sorted by AIC, but we index by name)
row_map <- list(
  xi       = "xi",
  theta_Sn = "theta_Sn",
  theta_In = "theta_In",
  theta_P  = "theta_P",
  probn    = "probn",
  rn       = "rn",
  f_Sn     = "f_Sn"
)

for (alt in alt_names) {
  cat("\n--- Alternative:", alt, "---\n")

  # Observed LRT statistic
  row_idx <- row_map[[alt]]
  alt_ll_obs <- dent_para_parameter_table[row_idx, "ll"]
  Lambda_obs <- 2 * (alt_ll_obs - null_ll_obs)
  cat("  Observed Lambda:", Lambda_obs, "\n")

  # Chi-square p-value (from existing table)
  p_chisq <- dent_para_parameter_table[row_idx, "lrt_pval"]

  # Load bootstrap alt results
  alt_lls <- numeric(B)
  alt_completed <- 0
  for (b in 1:B) {
    f <- sprintf("results_alt/lrt_%s_%d.rds", alt, b)
    if (file.exists(f)) {
      res <- readRDS(f)
      alt_lls[b] <- res$ll
      alt_completed <- alt_completed + 1
    } else {
      alt_lls[b] <- NA
    }
  }
  cat("  Alt fits completed:", alt_completed, "/", B, "\n")

  # Compute bootstrap LRT statistics (only for complete pairs)
  valid <- !is.na(null_lls) & !is.na(alt_lls)
  Lambda_boot <- 2 * (alt_lls[valid] - null_lls[valid])
  B_valid <- sum(valid)

  # Bootstrap p-value with correction
  p_boot <- (1 + sum(Lambda_boot >= Lambda_obs)) / (1 + B_valid)
  cat("  Bootstrap p-value:", p_boot, "(B =", B_valid, ")\n")
  cat("  Chi-square p-value:", p_chisq, "\n")

  results <- rbind(results, data.frame(
    alt_name = alt,
    Lambda_obs = Lambda_obs,
    p_chisq = ifelse(is.null(p_chisq), NA, p_chisq),
    p_boot = p_boot,
    B_completed = B_valid,
    stringsAsFactors = FALSE
  ))

  # Diagnostic histogram
  if (B_valid > 10) {
    # Compute df for chi-square overlay
    alt_AIC_obs <- dent_para_parameter_table[row_idx, "AIC"]
    k_alt <- (alt_AIC_obs + 2 * alt_ll_obs) / 2
    k_null <- (null_AIC_obs + 2 * null_ll_obs) / 2
    df <- round(k_alt - k_null)

    p <- ggplot(data.frame(Lambda = Lambda_boot), aes(x = Lambda)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", alpha = 0.7) +
      stat_function(fun = dchisq, args = list(df = max(df, 1)), color = "red", linewidth = 1) +
      geom_vline(xintercept = Lambda_obs, color = "darkgreen", linetype = "dashed", linewidth = 1) +
      labs(title = paste0("Bootstrap LRT: ", alt, " (df=", df, ")"),
           subtitle = paste0("p_boot=", round(p_boot, 3), ", p_chisq=", round(p_chisq, 3)),
           x = expression(Lambda), y = "Density") +
      theme_minimal()
    ggsave(sprintf("diagnostic_%s.pdf", alt), p, width = 6, height = 4)
  }
}

# ---- Save results ----
saveRDS(results, file = "bootstrap_pvalues.rds")
cat("\n\n=== Summary ===\n")
print(results)
cat("\nResults saved to bootstrap_pvalues.rds\n")
