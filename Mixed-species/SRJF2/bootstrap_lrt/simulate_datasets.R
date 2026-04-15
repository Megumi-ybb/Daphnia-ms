#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(pomp)
library(panelPomp)

set.seed(801)

B <- 100

# ---- Load MLE from best_result.rda ----
load("../best_result.rda")
# mif.estimate is loaded from the .rda file

# ---- Model specification (same C snippets as all_shared.R) ----
dyn_rpro <- Csnippet("
  double Sn_term, F_term, Si_term, Jn_term, Ji_term;
  double noiSn, noiSi, noiF, noiJn, noiJi;
  double delta = 0.013;

  noiSn = rnorm(0, sigSn * sqrt(dt));
  noiSi = rnorm(0, sigSi * sqrt(dt));
  noiJn = rnorm(0, sigJn * sqrt(dt));
  noiJi = rnorm(0, sigJi * sqrt(dt));
  noiF  = rnorm(0, sigF  * sqrt(dt));

  Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt - delta * Sn * dt + Sn * noiSn;
  Jn_term = rn * f_Sn * F * Sn * dt - 0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;
  Si_term = 0.1 * Ji * dt - theta_Si * Si * dt - delta * Si * dt + Si * noiSi;
  Ji_term = ri * f_Si * F * Si * dt - 0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi;
  F_term  = F * noiF - f_Sn * F * (Sn + 1 * Jn) * dt - f_Si * F * (Si + 1 * Ji) * dt - delta * F * dt + 0.37 * dt;

  F  += F_term;
  Sn += Sn_term;
  Si += Si_term;
  Ji += Ji_term;
  Jn += Jn_term;

  if (Sn < 0.0 || Sn > 1e5) { Sn = 0.0; error_count += 1; }
  if (Si < 0.0 || Si > 1e5) { Si = 0.0; error_count += 1000000; }
  if (F  < 0.0 || F  > 1e20){ F  = 0.0; error_count += 1000; }
  if (Jn < 0.0 || Jn > 1e5) { Jn = 0.0; error_count += 0.001; }
  if (Ji < 0.0 || Ji > 1e5) { Ji = 0.0; error_count += 0.000000001; }

  T_Sn = fabs(Sn);
  T_Si = fabs(Si);
")

dyn_init <- Csnippet("
  Sn = 2.333;
  Si = 0.667;
  F  = 16.667;
  Jn = 0;
  Ji = 0;
  T_Sn = 0.0;
  T_Si = 0.0;
  error_count = 0.0;
")

dmeas <- Csnippet("
  if (error_count > 0.0) {
    lik = -150;
  } else {
    if (give_log) {
      lik = dnbinom_mu(lumadult,k_Si,T_Si,give_log) + dnbinom_mu(dentadult,k_Sn,T_Sn,give_log);
    } else {
      lik = dnbinom_mu(lumadult,k_Si,T_Si,give_log) * dnbinom_mu(dentadult,k_Sn,T_Sn,give_log);
    }
  }
")

rmeas <- Csnippet("
  dentadult = rnbinom_mu(k_Sn, T_Sn);
  lumadult  = rnbinom_mu(k_Si, T_Si);
")

pt <- parameter_trans(
  log = c("sigSn","sigSi","sigF","f_Sn","f_Si",
          "rn","ri","k_Sn","k_Si","sigJi","sigJn",
          "theta_Sn","theta_Si","theta_Jn","theta_Ji")
)

param_names <- c("sigSn","sigSi","sigF","f_Sn","f_Si",
                 "rn","ri","k_Sn","k_Si","sigJi","sigJn",
                 "theta_Sn","theta_Si","theta_Jn","theta_Ji")
state_names <- c("Sn","Si","Jn","Ji","error_count","F","T_Sn","T_Si")

# ---- Observation times (same as original data) ----
obs_times <- (0:9) * 5 + 7

template_data <- data.frame(
  day       = obs_times,
  dentadult = rep(0, length(obs_times)),
  lumadult  = rep(0, length(obs_times))
)

parameters <- mif.estimate[param_names]

# ---- Build panelPomp with 9 units ----
pomplist <- list()
for (i in 1:9) {
  pomplist[[i]] <- pomp(
    data      = template_data,
    times     = "day",
    t0        = 1,
    rprocess  = euler(dyn_rpro, delta.t = 1/4),
    rinit     = dyn_init,
    dmeasure  = dmeas,
    rmeasure  = rmeas,
    partrans  = pt,
    obsnames  = c("dentadult","lumadult"),
    accumvars = c("error_count"),
    paramnames = param_names,
    statenames = state_names
  )
  coef(pomplist[[i]]) <- parameters
}
names(pomplist) <- paste0("u", 1:9)

panelfood <- panelPomp(pomplist, shared = mif.estimate)

dir.create("simulated_data", showWarnings = FALSE, recursive = TRUE)

cat("Simulating", B, "datasets from MLE...\n")

for (b in 1:B) {
  sim_data_list <- list()
  for (u in names(panelfood)) {
    unit_model <- unit_objects(panelfood)[[u]]
    sim <- pomp::simulate(unit_model, nsim = 1, format = "data.frame",
                          params = mif.estimate)
    sim_data_list[[u]] <- sim[, c("day","dentadult","lumadult")]
  }
  saveRDS(sim_data_list, file = sprintf("simulated_data/sim_data_%03d.rds", b))
  cat("  Dataset", b, "saved.\n")
}

saveRDS(mif.estimate, file = "simulated_data/true_params.rds")

cat("Done. All", B, "simulated datasets saved in simulated_data/\n")
cat("True log(ri) =", log(mif.estimate["ri"]), "\n")
