#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))

# Effect-size sanity check before launching the full 100-fit AIC sim.
#
# The AIC selection-rate study compares an all-shared SIRJPF2 (k=26) against
# a theta_In-unit-specific SIRJPF2 (k=33). AIC of alt is lower iff
#     2*(ll_alt - ll_null) > 2*(33 - 26) = 14.
# If the *oracle* log-likelihood gain (i.e. true alt params vs true shared
# geomean) does not comfortably exceed 7 nats on the data-generating panels,
# the simulation will mechanically report "AIC fails to detect unit-specific
# theta_In" when in fact AIC is working correctly and the effect is too
# small relative to the parameter penalty.
#
# This script evaluates the oracle log-likelihood gain on 5 of the
# unit-specific panels at HIGH Np (no mif2 search needed). If the average
# 2 * ll-gain >> 14 the chosen perturbation (+/- 0.5 log units) is fine;
# if it is comparable to or smaller than 14 the perturbation should be
# increased before the full study is submitted.
#
# Usage:
#   Rscript pilot_oracle.R
# Requires:
#   simulate_unit_specific.R has already produced simulated_data_unit_specific/

library(foreach)
library(doParallel)
registerDoParallel(cores = min(parallel::detectCores(), 20))
library(pomp)
library(panelPomp)

set.seed(802 * 1000 + 99999)  # Distinct from all fit_*.R seeds.

# ---- Load truth and the 5 pilot panels ----
truth_obj <- readRDS("simulated_data_unit_specific/true_params.rds")
shared_params <- truth_obj$shared
theta_In_unit_values <- truth_obj$theta_In_unit_values
b_pilot <- c(1, 13, 25, 37, 49)
cat("Pilot panels:", paste(b_pilot, collapse = ", "), "\n")
cat("Truth (unit-specific theta_In):\n"); print(theta_In_unit_values)
cat("Null geomean theta_In:", exp(mean(log(theta_In_unit_values))), "\n\n")

# ---- Model spec (must match fit_*.R exactly) ----
dyn_rpro <- Csnippet("
  double Sn_term, In_term, F_term, P_term, Si_term, Ii_term, Jn_term, Ji_term;
  double noiSn, noiIn, noiSi, noiIi, noiF, noiP, noiJn, noiJi;
  double delta = 0.013;

  noiSn = rnorm(0, sigSn * sqrt(dt));
  noiIn = rnorm(0, sigIn * sqrt(dt));
  noiSi = rnorm(0, sigSi * sqrt(dt));
  noiIi = rnorm(0, sigIi * sqrt(dt));
  noiJn = rnorm(0, sigJn * sqrt(dt));
  noiJi = rnorm(0, sigJi * sqrt(dt));
  noiF  = rnorm(0, sigF  * sqrt(dt));
  noiP  = rnorm(0, sigP  * sqrt(dt));

  Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt - probn * f_Sn * Sn * P * dt - delta * Sn * dt + Sn * noiSn;
  Jn_term = rn * f_Sn * F * Sn * dt - 0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;
  In_term = probn * f_Sn * Sn * P * dt - theta_In * In * dt - delta * In * dt + In * noiIn;
  Si_term = 0.1 * Ji * dt - theta_Si * Si * dt - probi * f_Si * Si * P * dt - delta * Si * dt + Si * noiSi;
  Ji_term = ri * f_Si * F * Si * dt - 0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi;
  Ii_term = probi * f_Si * Si * P * dt - theta_Ii * Ii * dt - delta * Ii * dt + Ii * noiIi;
  F_term  = F * noiF - f_Sn * F * (Sn + xi * In + 1 * Jn) * dt - f_Si * F * (Si + xi * Ii + 1 * Ji) * dt - delta * F * dt + 0.37 * dt;
  P_term  = 30 * theta_In * In * dt + 30 * theta_Ii * Ii * dt - f_Sn * (Sn + xi * In) * P * dt - f_Si * (Si + xi * Ii) * P * dt - theta_P * P * dt - delta * P * dt + P * noiP;

  F  += F_term;
  Sn += Sn_term;
  In += In_term;
  Si += Si_term;
  Ji += Ji_term;
  Jn += Jn_term;
  Ii += Ii_term;
  P  += P_term;

  if (t - 4.0 < 0.001 && t - 4.0 > -0.001) { P += 25; }

  if (Sn < 0.0 || Sn > 1e5) { Sn = 0.0; error_count += 1; }
  if (Si < 0.0 || Si > 1e5) { Si = 0.0; error_count += 1000000; }
  if (F  < 0.0 || F  > 1e20){ F  = 0.0; error_count += 1000; }
  if (In < 0.0 || In > 1e5) { In = 0.0; error_count += 0.001; }
  if (Ii < 0.0 || Ii > 1e5) { Ii = 0.0; error_count += 0.000000001; }
  if (Jn < 0.0 || Jn > 1e5) { Jn = 0.0; error_count += 0.001; }
  if (Ji < 0.0 || Ji > 1e5) { Ji = 0.0; error_count += 0.000000001; }
  if ((P < 0.0 || P > 1e20) && t > 3.9) { P = 0.0; error_count += 0.000001; }

  T_Sn = fabs(Sn); T_In = fabs(In); T_Si = fabs(Si); T_Ii = fabs(Ii);
")

dyn_init <- Csnippet("
  Sn = 2.333; Si = 0.667; F = 16.667;
  Jn = 0; Ji = 0;
  T_Sn = 0.0; T_Si = 0.0; T_In = 0.0; T_Ii = 0.0;
  In = 0.0; Ii = 0.0;
  error_count = 0.0; P = 0;
")

dmeas <- Csnippet("
  if (error_count > 0.0) {
    lik = -150;
  } else {
    if (give_log) {
      lik = dnbinom_mu(lumadult,k_Si,T_Si,give_log) + dnbinom_mu(luminf,k_Ii,T_Ii,give_log) + dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) + dnbinom_mu(dentinf,k_In,T_In,give_log);
    } else {
      lik = dnbinom_mu(lumadult,k_Si,T_Si,give_log) * dnbinom_mu(luminf,k_Ii,T_Ii,give_log) * dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) * dnbinom_mu(dentinf,k_In,T_In,give_log);
    }
  }
")

pt <- parameter_trans(
  log = c("sigSn","sigIn","sigSi","sigIi","sigF","sigP","f_Sn","f_Si",
          "rn","ri","k_Ii","k_In","k_Sn","k_Si","sigJi","sigJn","probn","probi",
          "theta_Sn","theta_In","theta_Si","theta_Ii","theta_P","theta_Jn","theta_Ji","xi")
)

param_names <- c("xi","sigSn","sigIn","sigSi","sigIi","sigF","sigP",
                 "theta_Sn","theta_In","theta_Si","theta_P","theta_Ii",
                 "f_Sn","f_Si","rn","ri","probn","probi",
                 "k_Ii","k_In","k_Sn","k_Si","sigJi","sigJn","theta_Jn","theta_Ji")
state_names <- c("Sn","In","Si","Jn","Ji","Ii","error_count","F","T_Sn","T_In","T_Si","T_Ii","P")

# ---- Build panel object for a given parameter setting (alt or null) ----
build_panel <- function(sim_data_list, theta_In_value_or_vec) {
  unit_specific <- length(theta_In_value_or_vec) == 8
  pomplist <- list()
  for (i in 1:8) {
    u_name <- paste0("u", i)
    dat <- sim_data_list[[u_name]]
    colnames(dat) <- c("day","dentadult","dentinf","lumadult","luminf")
    pomplist[[i]] <- pomp(
      data = dat, times = "day", t0 = 1,
      rprocess = euler(dyn_rpro, delta.t = 1/4),
      rinit = dyn_init, dmeasure = dmeas,
      partrans = pt,
      obsnames = c("dentadult","dentinf","lumadult","luminf"),
      accumvars = c("error_count"),
      paramnames = param_names, statenames = state_names
    )
    if (unit_specific) {
      pars_i <- c(shared_params, theta_In = theta_In_value_or_vec[i])
    } else {
      pars_i <- c(shared_params, theta_In = theta_In_value_or_vec)
    }
    names(pars_i)[length(pars_i)] <- "theta_In"
    coef(pomplist[[i]]) <- pars_i[param_names]
  }
  names(pomplist) <- paste0("u", 1:8)
  if (unit_specific) {
    specific_mat <- matrix(theta_In_value_or_vec, nrow = 1, ncol = 8, byrow = TRUE,
                          dimnames = list("theta_In", paste0("u", 1:8)))
    panelPomp(pomplist, shared = shared_params, specific = specific_mat)
  } else {
    full_shared <- c(shared_params, theta_In = theta_In_value_or_vec)
    names(full_shared)[length(full_shared)] <- "theta_In"
    panelPomp(pomplist, shared = full_shared[param_names])
  }
}

# ---- Oracle pfilter at the true alt and the null geomean ----
Np_pilot <- 5000   # large to suppress MC noise
Np_rep   <- 5     # average over a few pfilter replications

theta_In_null <- exp(mean(log(theta_In_unit_values)))

oracle <- function(panel_obj) {
  ll_reps <- replicate(Np_rep, unitLogLik(pfilter(panel_obj, Np = Np_pilot)))
  panel_logmeanexp(ll_reps, MARGIN = 1, se = TRUE)
}

rows <- list()
for (b in b_pilot) {
  cat("--- b =", b, "---\n")
  sim_data_list <- readRDS(sprintf("simulated_data_unit_specific/sim_data_%03d.rds", b))

  alt_panel  <- build_panel(sim_data_list, theta_In_unit_values)
  null_panel <- build_panel(sim_data_list, theta_In_null)

  ll_alt  <- oracle(alt_panel)
  ll_null <- oracle(null_panel)

  delta_ll <- ll_alt[1] - ll_null[1]
  Lambda_oracle <- 2 * delta_ll
  Delta_AIC_oracle <- 2 * (33 - 26) - 2 * delta_ll   # = 14 - 2*delta_ll

  cat(sprintf("  ll_alt  (true unit-specific theta_In) = %.2f (se %.2f)\n", ll_alt[1], ll_alt[2]))
  cat(sprintf("  ll_null (geomean theta_In)            = %.2f (se %.2f)\n", ll_null[1], ll_null[2]))
  cat(sprintf("  delta_ll (alt - null)                 = %.2f\n", delta_ll))
  cat(sprintf("  2*delta_ll (oracle Lambda)            = %.2f\n", Lambda_oracle))
  cat(sprintf("  Delta_AIC_oracle (alt - null)         = %.2f  (negative => alt preferred)\n\n",
              Delta_AIC_oracle))

  rows[[length(rows) + 1]] <- data.frame(
    b = b, ll_alt = ll_alt[1], se_alt = ll_alt[2],
    ll_null = ll_null[1], se_null = ll_null[2],
    delta_ll = delta_ll, Lambda_oracle = Lambda_oracle,
    Delta_AIC_oracle = Delta_AIC_oracle
  )
}

pilot <- do.call(rbind, rows)
cat("=== Summary across", length(b_pilot), "pilot panels ===\n")
print(pilot, row.names = FALSE)

cat(sprintf("\nMean oracle Lambda      = %.2f   (must clear 14 for AIC to favour alt)\n",
            mean(pilot$Lambda_oracle)))
cat(sprintf("Mean oracle Delta_AIC   = %.2f   (negative => unit-specific selected)\n",
            mean(pilot$Delta_AIC_oracle)))

decision <- if (mean(pilot$Lambda_oracle) > 14 + 2 * sd(pilot$Lambda_oracle)) {
  "OK: perturbation comfortably clears the 14-unit AIC penalty. Proceed with the full study."
} else if (mean(pilot$Lambda_oracle) > 14) {
  "BORDERLINE: oracle gain only marginally exceeds 14. Consider increasing log-offset spread to +/- 0.7 or +/- 0.8 before full submission."
} else {
  "BLOCK: oracle gain does NOT clear 14. The current +/- 0.5 log-unit perturbation is too small; the full study will mechanically conclude AIC fails to detect when the effect is simply too small. Re-run simulate_unit_specific.R with log_offsets = seq(-0.8, 0.8, length.out = 8) or wider, then re-run this pilot."
}
cat("\nDecision:", decision, "\n")

saveRDS(pilot, file = "pilot_oracle.rds")
write.csv(pilot, file = "pilot_oracle.csv", row.names = FALSE)
