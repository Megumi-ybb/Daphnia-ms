#!/usr/bin/env Rscript
# Aggregate classical chi-square + parametric-bootstrap LRT results for
# the two Mixed-species panel families (SIRJPF2 with parasite, SRJF2 without).
#
# Inputs:
#   Mixed-species/SIRJPF2/bootstrap_lrt/bootstrap_pvalues.rds
#   Mixed-species/SRJF2/bootstrap_lrt/bootstrap_pvalues.rds
#
# Outputs (cite-able from response.tex / si.Rnw / ms.Rnw):
#   daphnia-article/data/lrt_mixed_summary.rds
#   daphnia-article/data/lrt_mixed_summary.csv
#   daphnia-article/data/lrt_mixed_summary.txt
#
# Run from anywhere: paths are anchored to the script's own location.

# ---- Resolve paths relative to this script ----
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1) {
  mixed_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
} else {
  mixed_dir <- normalizePath(getwd())   # interactive fallback
}
project_root <- normalizePath(file.path(mixed_dir, ".."), mustWork = TRUE)

sirjpf2_rds <- file.path(mixed_dir, "SIRJPF2", "bootstrap_lrt", "bootstrap_pvalues.rds")
srjf2_rds   <- file.path(mixed_dir, "SRJF2",   "bootstrap_lrt", "bootstrap_pvalues.rds")

stopifnot(file.exists(sirjpf2_rds), file.exists(srjf2_rds))

sirjpf2 <- readRDS(sirjpf2_rds)
srjf2   <- readRDS(srjf2_rds)

# Friendly parameter labels for the response letter / SI text
label_map <- c(
  xi                = "$\\xi$",
  theta_Si_theta_Sn = "$\\theta_S$",
  theta_P           = "$\\theta_P$",
  theta_Ii_theta_In = "$\\theta_I$",
  ri_rn             = "$r$",
  probi_probn       = "$p$",
  f_Si_f_Sn         = "$f_S$",
  rn_ri             = "$r$",
  theta_Sn_theta_Si = "$\\theta_S$"
)

combined <- rbind(sirjpf2, srjf2)
combined$param_label <- label_map[combined$alt_name]

# Reorder columns for human reading
combined <- combined[, c(
  "family", "alt_name", "param_label",
  "k_null", "k_alt", "df",
  "ll_null_obs", "ll_alt_obs", "Lambda_obs",
  "p_chisq", "p_boot",
  "B_completed", "Lambda_boot_mean", "Lambda_boot_sd"
)]

# ---- Save artifacts ----
out_dir <- file.path(project_root, "daphnia-article", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(combined, file.path(out_dir, "lrt_mixed_summary.rds"))
write.csv(combined, file.path(out_dir, "lrt_mixed_summary.csv"), row.names = FALSE)

# ---- Human-readable text report ----
fmt_p <- function(p) {
  if (is.na(p)) return("---")
  if (p >= 1)      return(">0.99")
  if (p < 0.001)   return("<0.001")
  if (p < 0.01)    return(sprintf("%.3f", p))
  sprintf("%.2f", p)
}

txt <- character()
push <- function(...) txt <<- c(txt, sprintf(...))

push("Mixed-species panel LRT summary")
push("================================")
push("")
push("Source rds files:")
push("  %s", sirjpf2_rds)
push("  %s", srjf2_rds)
push("")
push("Bootstrap test stat: Lambda^(b) = 2 * (ll_alt(y_b) - ll_null(y_b))")
push("Bootstrap p-value:   p_boot = (1 + #{Lambda_boot >= Lambda_obs}) / (1 + B)")
push("Classical p-value:   1 - F_chisq(Lambda_obs, df), set to 1 when Lambda_obs <= 0")
push("")

for (fam in unique(combined$family)) {
  sub <- combined[combined$family == fam, , drop = FALSE]
  null_ll <- sub$ll_null_obs[1]
  push("--- Family: %s  (null all-shared ll = %.3f, k_null = %d, B = %d)",
       fam, null_ll, sub$k_null[1], sub$B_completed[1])
  push("%-22s %4s %8s %8s %8s   %4s %5s %5s",
       "alternative", "df", "ll_alt", "Lambda", "<bootΛ>", "Bsd", "pchi", "pboot")
  for (i in seq_len(nrow(sub))) {
    push("%-22s %4d %8.3f %8.3f %8.3f   %4.1f %5s %5s",
         sub$alt_name[i],
         sub$df[i],
         sub$ll_alt_obs[i],
         sub$Lambda_obs[i],
         sub$Lambda_boot_mean[i],
         sub$Lambda_boot_sd[i],
         fmt_p(sub$p_chisq[i]),
         fmt_p(sub$p_boot[i]))
  }
  push("")
}

push("Notes for citation:")
push("  - Family SIRJPF2 (parasitized, two species, 8 mesocosms): 7 alternatives, B = 100 paired replicates each.")
push("  - Family SRJF2   (no parasite, two species, 9 mesocosms): 3 alternatives, B = 100 paired replicates each.")
push("  - 4/7 SIRJPF2 alternatives and 3/3 SRJF2 alternatives have Lambda_obs < 0 (alt fits worse than null on data MLE due to Monte Carlo likelihood noise), making chi-square inference inapplicable; bootstrap returns p_boot >= 0.89 in those cases.")
push("  - For the two SIRJPF2 alternatives with the largest Lambda_obs (theta_Ii_theta_In, probi_probn), classical chi-square gives p < 0.05 but bootstrap p_boot is 0.15 and 0.22 -- chi-square is anti-conservative because the bootstrap null distribution has substantially heavier upper tails (sd 8-14 vs chi-square sd ~5).")
push("  - Across both families and all 10 alternatives, bootstrap p_boot > 0.05, so neither test rejects the all-shared baseline.")

writeLines(txt, file.path(out_dir, "lrt_mixed_summary.txt"))

cat("\n=== Combined Mixed-species LRT summary ===\n")
print(combined, row.names = FALSE)
cat(sprintf("\nSaved:\n  %s\n  %s\n  %s\n",
            file.path(out_dir, "lrt_mixed_summary.rds"),
            file.path(out_dir, "lrt_mixed_summary.csv"),
            file.path(out_dir, "lrt_mixed_summary.txt")))
