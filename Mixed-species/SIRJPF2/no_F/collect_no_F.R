#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
# Collect the no-F ablation fit(s), compute AIC, and compare against the
# all-shared SIRJPF2 model via Delta AIC. Run locally after the HPC job(s).
#
# Model accounting:
#   all-shared:  k = 24 free parameters  (sigSn = sigSi = 0 boundary-fixed, excluded)
#   no-F:        k = 23                  (sigF additionally dropped)
# The two models are non-nested (freezing F dynamics is structural), so AIC --
# not a likelihood-ratio test -- is the appropriate comparison.

K_NO_F <- 23

# ---- Best no-F fit across whatever rep files exist ----
rep_files <- list.files("results", pattern = "^no_F_[0-9]+\\.rds$", full.names = TRUE)
if (length(rep_files) == 0) stop("No results/no_F_*.rds files found. Run fit_no_F.R first.")

lls <- sapply(rep_files, function(f) readRDS(f)$ll)
cat("Found", length(rep_files), "no-F fit(s):\n")
for (i in seq_along(rep_files)) cat("  ", basename(rep_files[i]), "ll =", lls[i], "\n")

best_idx  <- which.max(lls)
best_fit  <- readRDS(rep_files[best_idx])
ll_no_F   <- best_fit$ll
se_no_F   <- best_fit$se
AIC_no_F  <- -2 * ll_no_F + 2 * K_NO_F

# ---- All-shared baseline from the SIRJPF2 parameter table ----
load("../../../data/Target_dynamics/para/Target_para_loglik_df.rds")
ll_all_shared  <- target_para_parameter_table["all_shared", "ll"]
AIC_all_shared <- target_para_parameter_table["all_shared", "AIC"]
k_all_shared   <- round((AIC_all_shared + 2 * ll_all_shared) / 2)   # = 24

delta_AIC <- AIC_no_F - AIC_all_shared   # > 0  =>  all-shared (with F) preferred

summary_df <- data.frame(
  model     = c("all_shared", "no_F"),
  k         = c(k_all_shared, K_NO_F),
  ll        = c(ll_all_shared, ll_no_F),
  ll_se     = c(NA_real_, se_no_F),
  AIC       = c(AIC_all_shared, AIC_no_F),
  delta_AIC = c(0, delta_AIC)
)

cat("\n==== No-F ablation vs all-shared SIRJPF2 ====\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nDelta AIC (no_F - all_shared) = %.2f\n", delta_AIC))
cat(if (delta_AIC > 0)
      sprintf("=> all-shared model (with depletable F) preferred by %.1f AIC units.\n", delta_AIC)
    else
      "=> no-F model preferred (unexpected); inspect the fit.\n")

# ---- Persist for the SI / response letter ----
out_dir <- "../../../daphnia-article/data"
saveRDS(list(summary = summary_df, delta_AIC = delta_AIC, best_coef = best_fit$coef),
        file = file.path(out_dir, "no_F_ablation_summary.rds"))
write.csv(summary_df, file = file.path(out_dir, "no_F_ablation_summary.csv"), row.names = FALSE)

sink(file.path(out_dir, "no_F_ablation_summary.txt"))
cat("No-F ablation of the all-shared SIRJPF2 model\n")
cat("=============================================\n\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nDelta AIC (no_F - all_shared) = %.2f\n", delta_AIC))
sink()

cat("\nWrote no_F_ablation_summary.{rds,csv,txt} to", out_dir, "\n")
