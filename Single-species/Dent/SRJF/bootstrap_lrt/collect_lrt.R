#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
# Collect bootstrap LRT results and compute p-values
# Run locally after all HPC jobs complete

library(ggplot2)

B <- 100

alt_names <- c("theta_Sn", "rn", "f_Sn")

# ---- LRT helper functions (from si.Rnw) ----
lrt_pvalue <- function(ll_alt, AIC_alt, ll_null, AIC_null) {
  k_alt <- (AIC_alt + 2 * ll_alt) / 2
  k_null <- (AIC_null + 2 * ll_null) / 2
  df <- round(k_alt - k_null)
  if (df <= 0) return(NA)
  LRT <- 2 * (ll_alt - ll_null)
  if (LRT < 0) return(1)
  pchisq(LRT, df = df, lower.tail = FALSE)
}

# ---- Load observed LRT statistics from parameter table ----
# The all-shared null has no unit-specific parameters, so blocking is a no-op
# for it; its plain ll / AIC ARE its block ll / AIC. The alt observed
# log-likelihoods are read from block_ll / block_AIC inside the loop below.
load("../../../../data/Simple_dynamics/Dent/no_para/Dent_no_para_loglik_df.rds")
null_ll_obs <- dent_no_para_parameter_table["all_shared", "ll"]
null_AIC_obs <- dent_no_para_parameter_table["all_shared", "AIC"]

cat("Observed null log-likelihood:", null_ll_obs, "\n\n")

# ---- Load bootstrap null results ----
null_lls <- numeric(B)
null_completed <- 0
null_meta <- NULL
for (b in 1:B) {
  f <- sprintf("results_null/lrt_null_%d.rds", b)
  if (file.exists(f)) {
    res <- readRDS(f)
    null_lls[b] <- res$ll
    null_completed <- null_completed + 1
    if (is.null(null_meta) && !is.null(res$Np)) null_meta <- res[c("Np", "Mp", "Np_rep", "block")]
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

# Row name mapping from parameter table to alt_name
# Row 1: all_shared (null)
# Row 2: specific_theta_Sn
# Row 3: specific_rn
# Row 4: specific_f_Sn
row_map <- list(
  theta_Sn = 2,
  rn       = 3,
  f_Sn     = 4
)

for (alt in alt_names) {
  cat("\n--- Alternative:", alt, "---\n")

  # Observed LRT statistic.
  # The bootstrap alt fits use block = TRUE, so the OBSERVED alt log-likelihood
  # must come from the BLOCK fit (block_ll / block_AIC) for the observed and
  # bootstrap Lambda to be on the same footing. Mirrors the corrected SIRJPF2
  # collect_lrt.R, which loads block_ll / block_AIC for the alt arm.
  row_idx <- row_map[[alt]]
  alt_ll_obs  <- dent_no_para_parameter_table[row_idx, "block_ll"]
  alt_AIC_obs <- dent_no_para_parameter_table[row_idx, "block_AIC"]
  Lambda_obs <- 2 * (alt_ll_obs - null_ll_obs)
  cat("  Observed Lambda:", Lambda_obs, "\n")

  # Chi-square p-value, computed from the block-fit alt vs the all-shared null
  # (the all-shared null has no unit-specific params, so block is a no-op and
  # its plain ll / AIC are the block ll / AIC).
  p_chisq <- lrt_pvalue(alt_ll_obs, alt_AIC_obs, null_ll_obs, null_AIC_obs)

  # Load bootstrap alt results
  alt_lls <- numeric(B)
  alt_completed <- 0
  alt_meta <- NULL
  for (b in 1:B) {
    f <- sprintf("results_alt/lrt_%s_%d.rds", alt, b)
    if (file.exists(f)) {
      res <- readRDS(f)
      alt_lls[b] <- res$ll
      alt_completed <- alt_completed + 1
      if (is.null(alt_meta) && !is.null(res$Np)) alt_meta <- res[c("Np", "Mp", "Np_rep", "block")]
    } else {
      alt_lls[b] <- NA
    }
  }
  cat("  Alt fits completed:", alt_completed, "/", B, "\n")

  # ---- Guard: null and alt must be fit at identical particle settings ----
  if (!is.null(null_meta) && !is.null(alt_meta)) {
    if (!isTRUE(all.equal(null_meta[c("Np", "Mp", "Np_rep")], alt_meta[c("Np", "Mp", "Np_rep")]))) {
      stop(sprintf("Particle-setting mismatch for %s: null (Np=%s, Mp=%s) vs alt (Np=%s, Mp=%s). Re-fit both arms at equal settings before collecting.",
                   alt, null_meta$Np, null_meta$Mp, alt_meta$Np, alt_meta$Mp))
    }
  } else {
    warning(sprintf("Cannot verify null/alt particle symmetry for %s (results predate metadata).", alt))
  }

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
    # Compute df for chi-square overlay (block fit, consistent with above)
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
