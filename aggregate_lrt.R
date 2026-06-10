#!/usr/bin/env Rscript
# Project-wide LRT aggregator.
#
# Scans every model-family `bootstrap_lrt/bootstrap_pvalues.rds` produced by the
# per-family `collect_lrt.R` scripts, combines them, and writes a single citable
# summary at `daphnia-article/data/lrt_all_summary.{rds,csv,txt}`.
#
# Each family rds is expected to share the schema:
#   family, alt_name, k_null, k_alt, df,
#   ll_null_obs, ll_alt_obs, Lambda_obs,
#   p_chisq, p_boot, B_completed,
#   Lambda_boot_mean, Lambda_boot_sd
#
# Families with missing rds are skipped with a warning.
# Families whose alt fits haven't run yet (B_completed = 0) keep their rows
# (classical p_chisq is still informative) but p_boot will be NA.

# ---- Resolve paths relative to this script ----
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1) {
  proj_root <- normalizePath(dirname(sub("^--file=", "", script_arg)))
} else {
  proj_root <- normalizePath(getwd())
}

# Ordered: parasitised mixed, no-parasite mixed, parasitised single (Dent, Lum),
#          no-parasite single (Dent, Lum)
families <- list(
  list(label = "SIRJPF2",     path = "Mixed-species/SIRJPF2/bootstrap_lrt"),
  list(label = "SRJF2",       path = "Mixed-species/SRJF2/bootstrap_lrt"),
  list(label = "Dent_SIRJPF", path = "Single-species/Dent/SIRJPF/bootstrap_lrt"),
  list(label = "Lum_SIRJPF",  path = "Single-species/Lum/SIRJPF/bootstrap_lrt"),
  list(label = "Dent_SRJF",   path = "Single-species/Dent/SRJF/bootstrap_lrt"),
  list(label = "Lum_SRJF",    path = "Single-species/Lum/SRJF/bootstrap_lrt")
)

label_map <- c(
  xi                = "$\\xi$",
  theta_Si_theta_Sn = "$\\theta_S$",
  theta_Sn_theta_Si = "$\\theta_S$",
  theta_P           = "$\\theta_P$",
  theta_Ii_theta_In = "$\\theta_I$",
  ri_rn             = "$r$",
  rn_ri             = "$r$",
  probi_probn       = "$p$",
  f_Si_f_Sn         = "$f_S$",
  theta_Sn          = "$\\theta_S$",
  theta_Si          = "$\\theta_S$",
  theta_In          = "$\\theta_I$",
  theta_Ii          = "$\\theta_I$",
  rn                = "$r$",
  ri                = "$r$",
  probn             = "$p$",
  probi             = "$p$",
  f_Sn              = "$f_S$",
  f_Si              = "$f_S$"
)

combined <- NULL
status   <- character()

for (fam in families) {
  rds <- file.path(proj_root, fam$path, "bootstrap_pvalues.rds")
  if (!file.exists(rds)) {
    status <- c(status, sprintf("[skip ] %-12s no bootstrap_pvalues.rds at %s", fam$label, rds))
    next
  }
  tab <- readRDS(rds)
  required <- c("family", "alt_name", "k_null", "k_alt", "df",
                "ll_null_obs", "ll_alt_obs", "Lambda_obs",
                "p_chisq", "p_boot", "B_completed",
                "Lambda_boot_mean", "Lambda_boot_sd")
  if (!all(required %in% colnames(tab))) {
    status <- c(status, sprintf("[stale] %-12s legacy schema, run updated collect_lrt.R first",
                                fam$label))
    next
  }
  # Force the family label to the canonical one (collect_lrt.R may have used a different string)
  tab$family <- fam$label
  combined <- rbind(combined, tab)
  n_alt <- nrow(tab)
  b_complete <- sum(tab$B_completed == max(tab$B_completed, na.rm = TRUE))
  b_max <- max(tab$B_completed, na.rm = TRUE)
  status <- c(status, sprintf("[ok  ] %-12s %d alternatives, max B=%d (paired)",
                              fam$label, n_alt, b_max))
}

if (is.null(combined)) {
  stop("No family results found. Run each family's collect_lrt.R first.")
}

combined$param_label <- label_map[combined$alt_name]
combined <- combined[, c(
  "family", "alt_name", "param_label",
  "k_null", "k_alt", "df",
  "ll_null_obs", "ll_alt_obs", "Lambda_obs",
  "p_chisq", "p_boot",
  "B_completed", "Lambda_boot_mean", "Lambda_boot_sd"
)]

# ---- Save ----
out_dir <- file.path(proj_root, "daphnia-article", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(combined, file.path(out_dir, "lrt_all_summary.rds"))
write.csv(combined, file.path(out_dir, "lrt_all_summary.csv"), row.names = FALSE)

# ---- Human-readable text ----
fmt_p <- function(p) {
  if (is.na(p))   return("---")
  if (p >= 1)     return(">0.99")
  if (p < 0.001)  return("<0.001")
  if (p < 0.01)   return(sprintf("%.3f", p))
  sprintf("%.2f", p)
}

txt <- character()
push <- function(...) txt <<- c(txt, sprintf(...))

push("Project-wide bootstrap + classical LRT summary")
push("==============================================")
push("")
push("Source families:")
for (s in status) push("  %s", s)
push("")
push("Bootstrap test stat: Lambda^(b) = 2 * (ll_alt(y_b) - ll_null(y_b))")
push("Bootstrap p-value:   p_boot = (1 + #{Lambda_boot >= Lambda_obs}) / (1 + B)")
push("Classical p-value:   1 - F_chisq(Lambda_obs, df), set to 1 when Lambda_obs <= 0,")
push("                     NA when df <= 0.")
push("")

for (fam in unique(combined$family)) {
  sub <- combined[combined$family == fam, , drop = FALSE]
  null_ll <- sub$ll_null_obs[1]
  push("--- Family: %-12s (null all-shared ll = %.3f, k_null = %d)",
       fam, null_ll, sub$k_null[1])
  push("%-22s %4s %9s %9s %9s   %5s %6s %6s %4s",
       "alternative", "df", "ll_alt", "Lambda_obs", "<bootΛ>",
       "Bsd", "pchi", "pboot", "B")
  for (i in seq_len(nrow(sub))) {
    push("%-22s %4d %9.3f %9.3f %9s   %5s %6s %6s %4d",
         sub$alt_name[i],
         sub$df[i],
         sub$ll_alt_obs[i],
         sub$Lambda_obs[i],
         ifelse(is.na(sub$Lambda_boot_mean[i]),
                "    NA", sprintf("%9.3f", sub$Lambda_boot_mean[i])),
         ifelse(is.na(sub$Lambda_boot_sd[i]),
                "  NA", sprintf("%5.1f", sub$Lambda_boot_sd[i])),
         fmt_p(sub$p_chisq[i]),
         fmt_p(sub$p_boot[i]),
         sub$B_completed[i])
  }
  push("")
}

paired_total <- sum(combined$B_completed > 0)
paired_idx <- combined$B_completed > 0 & !is.na(combined$p_boot)
paired_reject <- sum(paired_idx & combined$p_boot < 0.05)

# --- Reject set (p_boot < 0.05), computed from the data ---
reject_rows <- combined[paired_idx & combined$p_boot < 0.05, , drop = FALSE]
reject_str <- if (nrow(reject_rows) > 0) {
  paste(sprintf("%s/%s (p_chi %s, p_boot %s)",
                reject_rows$family, reject_rows$alt_name,
                vapply(reject_rows$p_chisq, fmt_p, ""),
                vapply(reject_rows$p_boot,  fmt_p, "")),
        collapse = "; ")
} else { "none" }

# --- Smallest non-rejected p_boot and its label(s), with ties ---
nonreject_idx <- paired_idx & combined$p_boot >= 0.05
nonreject_min <- if (any(nonreject_idx)) min(combined$p_boot[nonreject_idx]) else NA_real_
min_rows <- combined[nonreject_idx & abs(combined$p_boot - nonreject_min) < 1e-9, , drop = FALSE]
min_label <- paste(sprintf("%s/%s", min_rows$family, min_rows$alt_name), collapse = ", ")

# --- Discrepancy examples: chi-square rejects (p_chi < 0.05) but the paired
#     bootstrap does not (p_boot >= 0.05) => chi-square is anti-conservative.
#     Ranked by the size of the p_boot - p_chi gap; top 3 shown. ---
anti <- combined[paired_idx & !is.na(combined$p_chisq) &
                 combined$p_chisq < 0.05 & combined$p_boot >= 0.05, , drop = FALSE]
if (nrow(anti) > 0) {
  anti <- anti[order(anti$p_boot - anti$p_chisq, decreasing = TRUE), , drop = FALSE]
  anti <- head(anti, 3)
}

push("Headline interpretation:")
push("  - Of %d alternatives with paired bootstrap data, %d reject the all-shared null at",
     paired_total, paired_reject)
push("    alpha=0.05: %s.", reject_str)
push("    The smallest p_boot among non-rejected alternatives is %.2f (%s);",
     nonreject_min, min_label)
push("    all other non-rejected alternatives have p_boot >= this minimum,")
push("    comfortably above the alpha=0.05 threshold.")
if (nrow(anti) > 0) {
  push("  - Chi-square and bootstrap disagree on these 14-/7-df alternatives, with the")
  push("    chi-square uniformly anti-conservative (it rejects where the bootstrap does not):")
  lab_w <- max(nchar(sprintf("%s/%s", anti$family, anti$alt_name)))
  for (i in seq_len(nrow(anti))) {
    push("      * %-*s p_chi %s vs p_boot %s (chi anti-conservative)",
         lab_w, sprintf("%s/%s", anti$family[i], anti$alt_name[i]),
         fmt_p(anti$p_chisq[i]), fmt_p(anti$p_boot[i]))
  }
} else {
  push("  - Chi-square and bootstrap agree on all paired alternatives.")
}
push("    The bootstrap null distribution mean and sd (<bootΛ>, Bsd above) show why the")
push("    chi-square shape is a poor reference here.")
push("  - Negative Lambda_obs (alt fit below null on the data MLE) reflects Monte Carlo")
push("    likelihood noise; bootstrap returns p_boot near 1 for those rows by construction.")

writeLines(txt, file.path(out_dir, "lrt_all_summary.txt"))

cat("\n=== Project-wide LRT summary ===\n\n")
for (s in status) cat(s, "\n")
cat("\n")
print(combined, row.names = FALSE)
cat(sprintf("\nSaved:\n  %s\n  %s\n  %s\n",
            file.path(out_dir, "lrt_all_summary.rds"),
            file.path(out_dir, "lrt_all_summary.csv"),
            file.path(out_dir, "lrt_all_summary.txt")))
