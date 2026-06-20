#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(pomp)
library(panelPomp)

set.seed(801)

B <- 100

# ---- Load MLE from best_result.rda ----
load("../model/best_result.rda")
cat("MLE parameters loaded from best_result.rda\n")
cat("Parameters:", paste(names(mif.estimate), collapse = ", "), "\n")

# ---- Model specification (from model/all_shared.R, with probi naming) ----
dyn_rpro <- Csnippet("
  double Si_term, Ii_term, F_term, P_term, Ji_term;
  double noiSi, noiIi, noiF, noiP, noiJi;
  double delta = 0.013;

  noiSi = rnorm(0, sigSi * sqrt(dt));
  noiIi = rnorm(0, sigIi * sqrt(dt));
  noiJi = rnorm(0, sigJi * sqrt(dt));
  noiF  = rnorm(0, sigF  * sqrt(dt));
  noiP  = rnorm(0, sigP  * sqrt(dt));

  Si_term = 0.1 * Ji * dt - theta_Si * Si * dt - probi * f_Si * Si * P * dt - delta * Si * dt + Si * noiSi;
  Ji_term = ri * f_Si * F * Si * dt - 0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi;
  Ii_term = probi * f_Si * Si * P * dt - theta_Ii * Ii * dt - delta * Ii * dt + Ii * noiIi;
  F_term  = F * noiF - f_Si * F * (Si + xi * Ii + 1 * Ji) * dt - delta * F * dt + 0.37 * dt;
  P_term  = 30 * theta_Ii * Ii * dt - f_Si * (Si + xi * Ii) * P * dt - theta_P * P * dt - delta * P * dt + P * noiP;

  F  += F_term;
  Si += Si_term;
  Ii += Ii_term;
  Ji += Ji_term;
  P  += P_term;

  if (t - 4.0 < 0.001 && t - 4.0 > -0.001) { P += 25; }

  if (Si < 0.0 || Si > 1e5) { Si = 0.0; error_count += 1; }
  if (F  < 0.0 || F  > 1e20){ F  = 0.0; error_count += 1000; }
  if (Ii < 0.0 || Ii > 1e5) { Ii = 0.0; error_count += 0.001; }
  if (Ji < 0.0 || Ji > 1e5) { Ji = 0.0; error_count += 0.001; }
  if ((P < 0.0 || P > 1e20) && t > 3.9) { P = 0.0; error_count += 0.000001; }

  T_Si = fabs(Si);
  T_Ii = fabs(Ii);
")

dyn_init <- Csnippet("
  Si = 3;
  F  = 16.667;
  Ji = 0;
  T_Si = 0.0;
  T_Ii = 0.0;
  Ii = 0.0;
  error_count = 0.0;
  P = 0;
")

dmeas <- Csnippet("
  if (error_count > 0.0) {
    lik = -150;
  } else {
    if (give_log) {
      lik = dnbinom_mu(dentadult,k_Si,T_Si,give_log) + dnbinom_mu(dentinf,k_Ii,T_Ii,give_log);
    } else {
      lik = dnbinom_mu(dentadult,k_Si,T_Si,give_log) * dnbinom_mu(dentinf,k_Ii,T_Ii,give_log);
    }
  }
")

rmeas <- Csnippet("
  dentadult = rnbinom_mu(k_Si, T_Si);
  dentinf   = rnbinom_mu(k_Ii, T_Ii);
")

pt <- parameter_trans(
  log = c("sigSi","sigIi","sigF","sigP","f_Si",
          "ri","k_Ii","k_Si","sigJi","probi",
          "theta_Si","theta_Ii","theta_P","theta_Ji","xi")
)

param_names <- c("sigSi","sigIi","sigF","sigP","f_Si",
                 "ri","k_Ii","k_Si","sigJi","probi",
                 "theta_Si","theta_Ii","theta_P","theta_Ji","xi")
state_names <- c("Si","Ii","Ji","error_count","F","T_Si","T_Ii","P")

obs_times <- (0:9) * 5 + 7

template_data <- data.frame(
  day      = obs_times,
  dentadult = rep(0, length(obs_times)),
  dentinf   = rep(0, length(obs_times))
)

parameters <- mif.estimate[param_names]

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
    obsnames  = c("dentadult","dentinf"),
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
    sim_data_list[[u]] <- sim[, c("day","dentadult","dentinf")]
  }
  saveRDS(sim_data_list, file = sprintf("simulated_data/sim_data_%03d.rds", b))
  cat("  Dataset", b, "saved.\n")
}

saveRDS(mif.estimate, file = "simulated_data/true_params.rds")

cat("Done. All", B, "simulated datasets saved in simulated_data/\n")
cat("True log(ri) =", log(mif.estimate["ri"]), "\n")
