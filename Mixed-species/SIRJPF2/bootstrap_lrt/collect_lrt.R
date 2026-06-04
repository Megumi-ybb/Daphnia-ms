#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
# Collect bootstrap LRT results and compute p-values
# Run locally after all HPC jobs complete

library(ggplot2)

B <- 100

alt_names <- c("xi", "theta_Si_theta_Sn", "theta_P",
               "theta_Ii_theta_In", "ri_rn", "probi_probn", "f_Si_f_Sn")

# ---- Classical chi-square LRT p-value (mirrors si.Rnw lrt_pvalue helper) ----
lrt_pvalue <- function(ll_alt, AIC_alt, ll_null, AIC_null) {
  k_alt  <- (AIC_alt  + 2 * ll_alt)  / 2
  k_null <- (AIC_null + 2 * ll_null) / 2
  df <- round(k_alt - k_null)
  if (df <= 0) return(NA_real_)
  Lambda <- 2 * (ll_alt - ll_null)
  if (Lambda < 0) return(1)
  pchisq(Lambda, df = df, lower.tail = FALSE)
}

# ---- Load observed log-likelihoods from the SIRJPF2 parameter table ----
load("../../../data/Target_dynamics/para/Target_para_loglik_df.rds")
null_ll_obs  <- target_para_parameter_table["all_shared", "ll"]
null_AIC_obs <- target_para_parameter_table["all_shared", "AIC"]

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

# Row name mapping from parameter table to alt_name
row_map <- list(
  xi                = 5,
  theta_Si_theta_Sn = 6,
  theta_P           = 4,
  theta_Ii_theta_In = 2,
  ri_rn             = 8,
  probi_probn       = 3,
  f_Si_f_Sn         = 7
)

results <- data.frame(
  family       = character(),
  alt_name     = character(),
  k_null       = integer(),
  k_alt        = integer(),
  df           = integer(),
  ll_null_obs  = numeric(),
  ll_alt_obs   = numeric(),
  Lambda_obs   = numeric(),
  p_chisq      = numeric(),
  p_boot       = numeric(),
  B_completed  = integer(),
  Lambda_boot_mean = numeric(),
  Lambda_boot_sd   = numeric(),
  stringsAsFactors = FALSE
)

for (alt in alt_names) {
  cat("\n--- Alternative:", alt, "---\n")

  row_idx     <- row_map[[alt]]
  alt_ll_obs  <- target_para_parameter_table[row_idx, "ll"]
  alt_AIC_obs <- target_para_parameter_table[row_idx, "AIC"]
  Lambda_obs  <- 2 * (alt_ll_obs - null_ll_obs)
  k_alt       <- (alt_AIC_obs  + 2 * alt_ll_obs)  / 2
  k_null      <- (null_AIC_obs + 2 * null_ll_obs) / 2
  df          <- round(k_alt - k_null)
  p_chisq     <- lrt_pvalue(alt_ll_obs, alt_AIC_obs, null_ll_obs, null_AIC_obs)

  cat(sprintf("  k_null=%d  k_alt=%d  df=%d  Lambda_obs=%.3f\n",
              as.integer(k_null), as.integer(k_alt), df, Lambda_obs))

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
  # Lambda_boot = 2*(ll_alt - ll_null); the particle-filter loglik is biased by
  # ~ -O(1/Np), so unequal Np would bias Lambda (and p_boot) in the wrong
  # direction. Refuse to report a p-value built from mismatched arms.
  if (!is.null(null_meta) && !is.null(alt_meta)) {
    if (!isTRUE(all.equal(null_meta[c("Np", "Mp", "Np_rep")], alt_meta[c("Np", "Mp", "Np_rep")]))) {
      stop(sprintf(
        "Particle-setting mismatch for %s: null (Np=%s, Mp=%s) vs alt (Np=%s, Mp=%s). Re-fit both arms at equal settings before collecting.",
        alt, null_meta$Np, null_meta$Mp, alt_meta$Np, alt_meta$Mp))
    }
  } else {
    warning(sprintf(
      "Cannot verify null/alt particle symmetry for %s (results predate the Np/Mp metadata); ensure both arms used the same settings.", alt))
  }

  valid <- !is.na(null_lls) & !is.na(alt_lls)
  Lambda_boot <- 2 * (alt_lls[valid] - null_lls[valid])
  B_valid <- sum(valid)

  p_boot <- (1 + sum(Lambda_boot >= Lambda_obs)) / (1 + B_valid)
  cat(sprintf("  p_chisq = %.4f   p_boot = %.4f   (B = %d)\n", p_chisq, p_boot, B_valid))

  results <- rbind(results, data.frame(
    family       = "SIRJPF2",
    alt_name     = alt,
    k_null       = as.integer(k_null),
    k_alt        = as.integer(k_alt),
    df           = df,
    ll_null_obs  = null_ll_obs,
    ll_alt_obs   = alt_ll_obs,
    Lambda_obs   = Lambda_obs,
    p_chisq      = p_chisq,
    p_boot       = p_boot,
    B_completed  = B_valid,
    Lambda_boot_mean = mean(Lambda_boot),
    Lambda_boot_sd   = sd(Lambda_boot),
    stringsAsFactors = FALSE
  ))

  if (B_valid > 10) {
    p <- ggplot(data.frame(Lambda = Lambda_boot), aes(x = Lambda)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", alpha = 0.7) +
      stat_function(fun = dchisq, args = list(df = max(df, 1)), color = "red", linewidth = 1) +
      geom_vline(xintercept = Lambda_obs, color = "darkgreen", linetype = "dashed", linewidth = 1) +
      labs(title = paste0("SIRJPF2 bootstrap LRT: ", alt, " (df=", df, ")"),
           subtitle = sprintf("p_boot=%.3f, p_chisq=%.3f, B=%d",
                              p_boot, p_chisq, B_valid),
           x = expression(Lambda), y = "Density") +
      theme_minimal()
    ggsave(sprintf("diagnostic_%s.pdf", alt), p, width = 6, height = 4)
  }
}

# ---- Save ----
saveRDS(results, file = "bootstrap_pvalues.rds")
write.csv(results, file = "bootstrap_pvalues.csv", row.names = FALSE)
cat("\n\n=== Summary ===\n")
print(results, row.names = FALSE)
cat("\nResults saved to bootstrap_pvalues.{rds,csv}\n")
