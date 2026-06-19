library(stats)

set.seed(20260615)


################################# Toy exact-MLE functions #################################

# This function is the likelihood-ratio statistic for a simple Gaussian toy model.
# It compares:
#   H1: each mesocosm/block can have its own mean
#   H0: a random-effect model with a between-block variance component
#
# SSW = within-block sum of squares
# SSB = between-block sum of squares
# I   = number of blocks/mesocosms
# J   = number of observations inside each block
#
# The important point is the boundary condition:
# when the estimated random-effect variance would be negative, H0 collapses to
# the pooled-variance model.  That boundary behavior is one reason the usual
# chi-square reference can be too small.
calculate_Lambda_from_SS = function(SSW, SSB, I, J) {
  N = I * J
  
  H1_loglik = -N / 2 * (log(2 * pi) + log(SSW / N) + 1)
  
  within_variance = SSW / (I * (J - 1))
  between_variance = SSB / I
  interior_variance = between_variance >= within_variance
  
  H0_loglik = numeric(length(SSW))
  
  H0_loglik[interior_variance] =
    -N / 2 * log(2 * pi) -
    I / 2 * log(between_variance[interior_variance]) -
    I * (J - 1) / 2 * log(within_variance[interior_variance]) -
    N / 2
  
  pooled_variance = (SSW + SSB) / N
  H0_loglik[!interior_variance] =
    -N / 2 * (log(2 * pi) + log(pooled_variance[!interior_variance]) + 1)
  
  Lambda = 2 * (H1_loglik - H0_loglik)
  return(Lambda)
}


# Simulate one independent block of the toy model.
#
# phi2 controls how strong the block/mesocosm effect is.  Larger phi2 gives
# larger between-block variation, so the fixed-effect model gains more loglik
# over the random-effect/null model.
simulate_one_toy_block = function(Bn, I, J, phi2) {
  sum_square_matrix = vapply(seq_len(Bn), function(rep_idx) {
    Y = matrix(rnorm(I * J), I, J) + rnorm(I, 0, sqrt(phi2))
    row_mean = rowMeans(Y)
    
    SSW = sum((Y - row_mean)^2)
    SSB = J * sum((row_mean - mean(Y))^2)
    
    c(SSW, SSB)
  }, numeric(2))
  
  toy_Lambda = calculate_Lambda_from_SS(sum_square_matrix[1, ],
                                        sum_square_matrix[2, ],
                                        I,
                                        J)
  return(toy_Lambda)
}


# Sum independent toy blocks.  This matches the LRT df structure:
# one block gives roughly 7 df, two blocks give roughly 14 df.
simulate_total_toy_Lambda = function(Bn, I, J, phi2, number_of_blocks) {
  toy_Lambda = numeric(Bn)
  
  for (block_idx in seq_len(number_of_blocks)) {
    toy_Lambda = toy_Lambda + simulate_one_toy_block(Bn, I, J, phi2)
  }
  
  return(toy_Lambda)
}


# rho is an easier scale to read than phi2:
#   rho = 1 + J * phi2
#
# rho = 1 is the calibration control with no extra block effect.
# rho > 1 means stronger block heterogeneity.
estimate_toy_mean_Lambda = function(rho, number_of_blocks) {
  phi2 = (rho - 1) / J
  mean_toy_Lambda = mean(simulate_total_toy_Lambda(4000, I, J, phi2, number_of_blocks))
  return(mean_toy_Lambda)
}


# For each real bootstrap mean, find the toy rho that gives the same mean LRT.
# This gives a rough "how much random-effect heterogeneity would explain this?"
# scale.
find_implied_rho = function(target_mean_Lambda, number_of_blocks) {
  toy_floor = estimate_toy_mean_Lambda(1.001, number_of_blocks)
  
  if (target_mean_Lambda <= toy_floor) {
    return(NA_real_)
  }
  
  rho_root = tryCatch(
    uniroot(function(rho) {
      estimate_toy_mean_Lambda(rho, number_of_blocks) - target_mean_Lambda
    }, c(1.001, 80))$root,
    error = function(e) NA_real_
  )
  
  return(rho_root)
}


################################# Analysis setup #################################

B = 100
I = 8
J = 10

alternative_names = c("xi",
                      "theta_Si_theta_Sn",
                      "theta_P",
                      "theta_Ii_theta_In",
                      "ri_rn",
                      "probi_probn",
                      "f_Si_f_Sn")

target_table_row = list(xi = 5,
                        theta_Si_theta_Sn = 6,
                        theta_P = 4,
                        theta_Ii_theta_In = 2,
                        ri_rn = 8,
                        probi_probn = 3,
                        f_Si_f_Sn = 7)


################################# Observed model log likelihoods #################################

# This table is produced by util/paramter_table.R.
# It contains the observed shared model and observed block-specific alternatives.
load("../../../data/Target_dynamics/para/Target_para_loglik_df.rds")
target_loglik_table = target_para_parameter_table

null_observed_loglik = target_loglik_table["all_shared", "ll"]
null_number_parameter = (target_loglik_table["all_shared", "AIC"] +
                           2 * null_observed_loglik) / 2


################################# Read bootstrap fits #################################

read_one_bootstrap_fit = function(file_name) {
  if (file.exists(file_name)) {
    return(readRDS(file_name))
  } else {
    return(NULL)
  }
}


extract_from_fit_list = function(fit_list, variable_name) {
  value_vector = sapply(fit_list, function(fit_object) {
    if (is.null(fit_object)) {
      return(NA)
    }
    
    if (variable_name %in% c("ll", "se")) {
      return(fit_object[[variable_name]])
    }
    
    return(unname(fit_object$coef[variable_name]))
  })
  
  return(value_vector)
}


null_fit_list = lapply(seq_len(B), function(bootstrap_idx) {
  file_name = paste0("results_null/lrt_null_", bootstrap_idx, ".rds")
  read_one_bootstrap_fit(file_name)
})

null_bootstrap_loglik = extract_from_fit_list(null_fit_list, "ll")
null_bootstrap_se = extract_from_fit_list(null_fit_list, "se")
null_sigP = extract_from_fit_list(null_fit_list, "sigP")
null_sigF = extract_from_fit_list(null_fit_list, "sigF")


alternative_bootstrap = lapply(alternative_names, function(alternative_name) {
  alt_fit_list = lapply(seq_len(B), function(bootstrap_idx) {
    file_name = paste0("results_alt/lrt_", alternative_name, "_", bootstrap_idx, ".rds")
    read_one_bootstrap_fit(file_name)
  })
  
  good_fit_idx = which(!sapply(alt_fit_list, is.null))
  first_good_fit = if (length(good_fit_idx) == 0) NULL else alt_fit_list[[good_fit_idx[1]]]
  
  list(
    loglik = extract_from_fit_list(alt_fit_list, "ll"),
    se = extract_from_fit_list(alt_fit_list, "se"),
    sigP = extract_from_fit_list(alt_fit_list, "sigP"),
    sigF = extract_from_fit_list(alt_fit_list, "sigF"),
    ncoef = if (is.null(first_good_fit)) NA else length(first_good_fit$coef),
    block = if (is.null(first_good_fit)) NA else isTRUE(first_good_fit$block)
  )
})

names(alternative_bootstrap) = alternative_names


################################# Main LRT diagnostic table #################################

LRT_table = do.call(rbind, lapply(alternative_names, function(alternative_name) {
  row_idx = target_table_row[[alternative_name]]
  
  alt_observed_loglik = target_loglik_table[row_idx, "block_ll"]
  alt_number_parameter = (target_loglik_table[row_idx, "block_AIC"] +
                            2 * alt_observed_loglik) / 2
  
  LRT_df = round(alt_number_parameter - null_number_parameter)
  observed_Lambda = 2 * (alt_observed_loglik - null_observed_loglik)
  
  paired_success = !is.na(null_bootstrap_loglik) &
    !is.na(alternative_bootstrap[[alternative_name]]$loglik)
  
  bootstrap_Lambda = 2 * (alternative_bootstrap[[alternative_name]]$loglik[paired_success] -
                            null_bootstrap_loglik[paired_success])
  
  exceedance_count = sum(bootstrap_Lambda >= observed_Lambda)
  number_success = sum(paired_success)
  
  # Add-one smoothing avoids reporting p = 0 when B is only 100.
  bootstrap_p_value = (1 + exceedance_count) / (1 + number_success)
  bootstrap_p_value_CI = binom.test(exceedance_count, number_success)$conf.int
  
  data.frame(
    alt = alternative_name,
    df = LRT_df,
    Lambda_obs = round(observed_Lambda, 2),
    boot_mean = round(mean(bootstrap_Lambda), 2),
    infl = round(mean(bootstrap_Lambda) / LRT_df, 2),
    var_2df = round(var(bootstrap_Lambda) / (2 * LRT_df), 2),
    altse_nullse = round(mean(alternative_bootstrap[[alternative_name]]$se[paired_success]) /
                           mean(null_bootstrap_se[paired_success]), 2),
    negfrac = round(mean(bootstrap_Lambda < 0), 2),
    p_chisq = signif(pchisq(observed_Lambda, LRT_df, lower.tail = FALSE), 3),
    p_boot = round(bootstrap_p_value, 3),
    p_boot_CI = paste0("[",
                       round(bootstrap_p_value_CI[1], 3),
                       ",",
                       round(bootstrap_p_value_CI[2], 3),
                       "]"),
    rho_star = round(find_implied_rho(mean(bootstrap_Lambda), LRT_df / 7), 2),
    stringsAsFactors = FALSE
  )
}))

cat("\n################################# S1: observed LRT vs bootstrap null #################################\n")
cat("If Wilks' theorem worked here, infl should be close to 1 and var_2df close to 1.\n")
cat("Large infl means the chi-square null is too small, so p_chisq is anti-conservative.\n\n")
print(LRT_table, row.names = FALSE)


################################# S2: toy calibration #################################

toy_calibration_table = do.call(rbind, lapply(c(1, 2, 3, 5, 10), function(rho) {
  toy_Lambda = simulate_total_toy_Lambda(12000, I, J, (rho - 1) / J, 2)
  
  data.frame(
    rho = rho,
    phi2 = round((rho - 1) / J, 3),
    meanL14 = round(mean(toy_Lambda), 2),
    infl = round(mean(toy_Lambda) / 14, 2),
    size_chisq14 = round(mean(toy_Lambda > qchisq(0.95, 14)), 3)
  )
}))

cat("\n################################# S2: exact-MLE toy calibration #################################\n")
cat("rho = 1 is the no-extra-block-effect control.  As rho increases, the nominal chi-square test rejects too often.\n\n")
print(toy_calibration_table, row.names = FALSE)


theta_I_target_mean = LRT_table$boot_mean[LRT_table$alt == "theta_Ii_theta_In"]
theta_I_rho_star = find_implied_rho(theta_I_target_mean, 2)
theta_I_toy_match = simulate_total_toy_Lambda(15000, I, J, (theta_I_rho_star - 1) / J, 2)
theta_I_observed_Lambda = LRT_table$Lambda_obs[LRT_table$alt == "theta_Ii_theta_In"]

cat("\n################################# S2: theta_I matched toy mechanism #################################\n")
cat("This asks: if the toy model has the same mean LRT as theta_I, how bad is the chi-square reference?\n")
cat(paste0("rho_star = ", round(theta_I_rho_star, 2),
           "; toy_mean = ", round(mean(theta_I_toy_match), 2),
           "; real_boot_mean = ", theta_I_target_mean,
           "; nominal_5pct_test_size = ", round(mean(theta_I_toy_match > qchisq(0.95, 14)), 2),
           "; nominal_1pct_test_size = ", round(mean(theta_I_toy_match > qchisq(0.99, 14)), 3),
           "; observed_Lambda_quantile = ", round(mean(theta_I_toy_match <= theta_I_observed_Lambda), 2),
           "\n"))


################################# S3: process-noise collapse #################################

noise_collapse_table = do.call(rbind, lapply(alternative_names, function(alternative_name) {
  paired_success = !is.na(null_bootstrap_loglik) &
    !is.na(alternative_bootstrap[[alternative_name]]$loglik)
  
  data.frame(
    alt = alternative_name,
    sigP_alt_null = round(mean(alternative_bootstrap[[alternative_name]]$sigP[paired_success]) /
                            mean(null_sigP[paired_success]), 2),
    sigF_alt_null = round(mean(alternative_bootstrap[[alternative_name]]$sigF[paired_success]) /
                            mean(null_sigF[paired_success]), 2),
    infl = round(LRT_table$boot_mean[LRT_table$alt == alternative_name] /
                   LRT_table$df[LRT_table$alt == alternative_name], 2),
    stringsAsFactors = FALSE
  )
}))

cat("\n################################# S3: process-noise collapse under H1 #################################\n")
cat("A ratio below 1 means the alternative model absorbed variation that the null model left as process noise.\n")
cat("This is the model-specific fingerprint for why the LR null is inflated.\n\n")
print(noise_collapse_table, row.names = FALSE)

minimum_noise_ratio = pmin(noise_collapse_table$sigP_alt_null,
                           noise_collapse_table$sigF_alt_null)

cat(paste0("Inflation vs minimum noise ratio: Pearson = ",
           round(cor(noise_collapse_table$infl, minimum_noise_ratio), 2),
           "; Spearman = ",
           round(cor(noise_collapse_table$infl, minimum_noise_ratio, method = "spearman"), 2),
           ". Smaller noise ratio means stronger collapse; stronger collapse corresponds to larger inflation if the correlation is negative.\n"))


################################# S4: particle-filter and optimizer checks #################################

pf_noise_table = do.call(rbind, lapply(alternative_names, function(alternative_name) {
  paired_success = !is.na(null_bootstrap_loglik) &
    !is.na(alternative_bootstrap[[alternative_name]]$loglik)
  
  bootstrap_Lambda = 2 * (alternative_bootstrap[[alternative_name]]$loglik[paired_success] -
                            null_bootstrap_loglik[paired_success])
  
  pf_sd = mean(sqrt(4 * (alternative_bootstrap[[alternative_name]]$se[paired_success]^2 +
                           null_bootstrap_se[paired_success]^2)))
  
  data.frame(
    alt = alternative_name,
    sd_Lambda = round(sd(bootstrap_Lambda), 2),
    pf_sd = round(pf_sd, 2),
    pf_var_pct = round(100 * pf_sd^2 / var(bootstrap_Lambda), 2),
    altse_nullse = round(mean(alternative_bootstrap[[alternative_name]]$se[paired_success]) /
                           mean(null_bootstrap_se[paired_success]), 2),
    stringsAsFactors = FALSE
  )
}))

cat("\n################################# S4a: particle-filter noise contribution #################################\n")
cat("pf_var_pct is the approximate percent of bootstrap LR variance explained by particle-filter loglik SE.\n")
cat("Small values mean filter noise is not the main reason the bootstrap null is wide.\n\n")
print(pf_noise_table, row.names = FALSE)


theta_P_success = !is.na(null_bootstrap_loglik) &
  !is.na(alternative_bootstrap[["theta_P"]]$loglik)
theta_P_Lambda = 2 * (alternative_bootstrap[["theta_P"]]$loglik[theta_P_success] -
                        null_bootstrap_loglik[theta_P_success])
toy_negative_count = sum(simulate_total_toy_Lambda(50000, I, J, (1.93 - 1) / J, 1) < 0)

cat("\n################################# S4b: optimizer slack check #################################\n")
cat("Negative Lambda values are possible from numerical optimizer/filter slack, because the fitted alternative can be worse than the fitted null.\n")
cat("But the exact toy has essentially no negatives, so negatives explain local noise, not the large positive mean inflation.\n")
cat(paste0("theta_P negative Lambda count = ", sum(theta_P_Lambda < 0),
           "/", sum(theta_P_success),
           "; minimum theta_P Lambda = ", round(min(theta_P_Lambda), 2),
           "; exact toy negatives in 50000 = ", toy_negative_count,
           "\n"))

cat("\n################################# S4c: inflation is not only df #################################\n")
cat("xi and theta_P both have df = 7, but different bootstrap inflation.\n")
cat(paste0("xi infl = ", LRT_table$infl[LRT_table$alt == "xi"],
           "; theta_P infl = ", LRT_table$infl[LRT_table$alt == "theta_P"],
           "\n"))


################################# S5: bootstrap state #################################

cat("\n################################# S5: bootstrap fit state #################################\n")
cat("This checks that theta_I uses the intended block-vs-block paired bootstrap comparison.\n")
cat(paste0("Null fits = ", sum(!is.na(null_bootstrap_loglik)), "/", B,
           "; theta_I alt block = ", alternative_bootstrap[["theta_Ii_theta_In"]]$block,
           "; theta_I ncoef = ", alternative_bootstrap[["theta_Ii_theta_In"]]$ncoef,
           "; matched p_boot = ",
           LRT_table$p_boot[LRT_table$alt == "theta_Ii_theta_In"],
           "\n"))


################################# S6: theta_I uncertainty and adverse checks #################################

theta_I_success = !is.na(null_bootstrap_loglik) &
  !is.na(alternative_bootstrap[["theta_Ii_theta_In"]]$loglik)
theta_I_bootstrap_Lambda = 2 * (alternative_bootstrap[["theta_Ii_theta_In"]]$loglik[theta_I_success] -
                                  null_bootstrap_loglik[theta_I_success])

theta_I_exceedance_count = sum(theta_I_bootstrap_Lambda >= theta_I_observed_Lambda)
theta_I_boot_mean = mean(theta_I_bootstrap_Lambda)
theta_I_CI = binom.test(theta_I_exceedance_count, length(theta_I_bootstrap_Lambda))$conf.int

# Adverse check: keep the mean but shrink the variance to a toy ceiling.
# If the conclusion changes under this artificial variance channel, the result is thin.
theta_I_adverse_Lambda =
  theta_I_boot_mean +
  (theta_I_bootstrap_Lambda - theta_I_boot_mean) *
  sqrt(1.34 * 2 * 14 / var(theta_I_bootstrap_Lambda))

theta_I_adverse_p =
  (1 + sum(theta_I_adverse_Lambda >= theta_I_observed_Lambda)) /
  (1 + length(theta_I_adverse_Lambda))

theta_I_bootstrap_p =
  (1 + theta_I_exceedance_count) / (1 + length(theta_I_bootstrap_Lambda))

cat("\n################################# S6a: variance-channel adverse check #################################\n")
cat("This keeps the theta_I bootstrap mean but shrinks variance toward the toy ceiling.\n")
cat("If p_boot moves close to 0.05, the result is not comfortably separated from rejection.\n")
cat(paste0("var/2df = ", round(var(theta_I_bootstrap_Lambda) / (2 * 14), 2),
           "; p_boot = ", round(theta_I_bootstrap_p, 3),
           "; adverse p_boot = ", round(theta_I_adverse_p, 3),
           "\n"))

cat("\n################################# S6b: depth asymmetry #################################\n")
cat("The observed numerator used deeper settings than the bootstrap fits: Np = 1e4 with 360 starts vs bootstrap Np = 1500 with 300 starts.\n")
cat("The bias sign is ambiguous, so the clean fix is matched-depth refitting, not just more bootstrap replicates.\n")

cat("\n################################# S6c: Monte Carlo uncertainty from B = 100 #################################\n")
cat("With B = 100, the bootstrap p-value has visible binomial uncertainty.\n")
cat(paste0("theta_I exceedance k = ", theta_I_exceedance_count,
           "/", length(theta_I_bootstrap_Lambda),
           "; p_boot = ", round(theta_I_bootstrap_p, 3),
           "; CI = [", round(theta_I_CI[1], 3), ",", round(theta_I_CI[2], 3), "]",
           "; Jeffreys_P(p < 0.05) = ",
           round(pbeta(0.05,
                       theta_I_exceedance_count + 0.5,
                       length(theta_I_bootstrap_Lambda) - theta_I_exceedance_count + 0.5), 3),
           "; 5th largest Lambda = ",
           round(sort(theta_I_bootstrap_Lambda, decreasing = TRUE)[5], 2),
           "\n"))

cat("\n################################# S6d: combined theta_I conclusion #################################\n")
cat(paste0("theta_I is not rejected by the paired bootstrap: p_boot = ",
           round(theta_I_bootstrap_p, 3),
           ". The margin is thin because the adverse check is about ",
           round(theta_I_adverse_p, 2),
           ". Certify with B around 500 and matched-depth refits.\n"))


################################# S7: final interpretation #################################

cat("\n################################# S7: final interpretation #################################\n")
cat("Use the paired bootstrap as the primary LRT reference.  The chi-square reference is anti-conservative here.\n")
cat(paste0("The smallest bootstrap p-value is theta_I p_boot = ",
           round(theta_I_bootstrap_p, 3),
           "; adverse p_boot is about ",
           round(theta_I_adverse_p, 2),
           "; toy chi-square true size at rho_star is about ",
           round(mean(theta_I_toy_match > qchisq(0.95, 14)), 2),
           ". Nothing rejects under the bootstrap reference.\n"))


################################# Figure: real bootstrap vs toy vs chi-square #################################

pdf("re_vs_fe_mechanism.pdf", width = 7.2, height = 4.6)
old_par = par(mar = c(4.3, 4.3, 3, 1))

real_boot_density = density(theta_I_bootstrap_Lambda)
toy_match_density = density(theta_I_toy_match)
x_grid = seq(0, max(60, theta_I_observed_Lambda + 5), length.out = 400)

plot(toy_match_density,
     col = "firebrick",
     lwd = 2,
     xlim = range(x_grid),
     ylim = c(0, max(real_boot_density$y,
                     toy_match_density$y,
                     dchisq(12, 14)) * 1.05),
     main = "SIRJPF2 theta_I: null log-LR vs chi-square",
     xlab = expression(Lambda == 2 * (hat(l)[1] - hat(l)[0])),
     ylab = "density")

lines(real_boot_density, col = "steelblue", lwd = 2)
lines(x_grid, dchisq(x_grid, 14), col = "black", lty = 2, lwd = 2)
abline(v = 14, col = "grey50", lty = 3)
abline(v = theta_I_observed_Lambda, col = "darkgreen", lwd = 2, lty = 2)

legend("topright",
       bty = "n",
       lwd = 2,
       col = c("firebrick", "steelblue", "black", "darkgreen"),
       lty = c(1, 1, 2, 2),
       legend = c(paste0("toy rho* = ", round(theta_I_rho_star, 1)),
                  "bootstrap (real)",
                  "chi-square df = 14",
                  paste0("observed ", round(theta_I_observed_Lambda, 1))))

par(old_par)
invisible(dev.off())

cat("\nFigure written: re_vs_fe_mechanism.pdf\n")
