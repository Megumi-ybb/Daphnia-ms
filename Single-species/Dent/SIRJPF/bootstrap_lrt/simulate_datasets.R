#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(pomp)
library(panelPomp)

set.seed(801)

B <- 100

# ---- Load MLE from model fitting ----
load("../model/all_shared.RData")

mif.estimate <- c(
  rn = 69.34499, f_Sn = 0.0007281269, probn = 0.3193476, xi = 14.43326,
  theta_Sn = 0.01660319, theta_In = 0.542057, theta_P = 7.917627e-05, theta_Jn = 0.002125109,
  sigSn = 0, sigIn = 0.5887923, sigJn = 0.3089306, sigF = 0.07980189,
  sigP = 0.4788064, k_In = 1.117379, k_Sn = 12.36433
)

# ---- Model specification (same C snippets as model/all_shared.R) ----
dyn_rpro <- Csnippet("
  double Sn_term, In_term, F_term, P_term, Jn_term;
  double noiSn, noiIn, noiF, noiP, noiJn;
  double delta = 0.013;

  noiSn = rnorm(0, sigSn * sqrt(dt));
  noiIn = rnorm(0, sigIn * sqrt(dt));
  noiJn = rnorm(0, sigJn * sqrt(dt));
  noiF  = rnorm(0, sigF  * sqrt(dt));
  noiP  = rnorm(0, sigP  * sqrt(dt));

  Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt - probn * f_Sn * Sn * P * dt - delta * Sn * dt + Sn * noiSn;
  Jn_term = rn * f_Sn * F * Sn * dt - 0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;
  In_term = probn * f_Sn * Sn * P * dt - theta_In * In * dt - delta * In * dt + In * noiIn;
  F_term  = F * noiF - f_Sn * F * (Sn + xi * In + 1 * Jn) * dt - delta * F * dt + 0.37 * dt;
  P_term  = 30 * theta_In * In * dt - f_Sn * (Sn + xi * In) * P * dt - theta_P * P * dt - delta * P * dt + P * noiP;

  F  += F_term;
  Sn += Sn_term;
  In += In_term;
  Jn += Jn_term;
  P  += P_term;

  if (t - 4.0 < 0.001 && t - 4.0 > -0.001) { P += 25; }

  if (Sn < 0.0 || Sn > 1e5) { Sn = 0.0; error_count += 1; }
  if (F  < 0.0 || F  > 1e20){ F  = 0.0; error_count += 1000; }
  if (In < 0.0 || In > 1e5) { In = 0.0; error_count += 0.001; }
  if (Jn < 0.0 || Jn > 1e5) { Jn = 0.0; error_count += 0.001; }
  if ((P < 0.0 || P > 1e20) && t > 3.9) { P = 0.0; error_count += 0.000001; }

  T_Sn = fabs(Sn);
  T_In = fabs(In);
")

dyn_init <- Csnippet("
  Sn = 3;
  F  = 16.667;
  Jn = 0;
  T_Sn = 0.0;
  T_In = 0.0;
  In = 0.0;
  error_count = 0.0;
  P = 0;
")

dmeas <- Csnippet("
  if (error_count > 0.0) {
    lik = -150;
  } else {
    if (give_log) {
      lik = dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) + dnbinom_mu(dentinf,k_In,T_In,give_log);
    } else {
      lik = dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) * dnbinom_mu(dentinf,k_In,T_In,give_log);
    }
  }
")

rmeas <- Csnippet("
  dentadult = rnbinom_mu(k_Sn, T_Sn);
  dentinf   = rnbinom_mu(k_In, T_In);
")

pt <- parameter_trans(
  log = c("sigSn","sigIn","sigF","sigP","f_Sn",
          "rn","k_In","k_Sn","sigJn","probn",
          "theta_Sn","theta_In","theta_P","theta_Jn","xi")
)

obs_times <- (0:9) * 5 + 7

template_data <- data.frame(
  day       = obs_times,
  dentadult = rep(0, length(obs_times)),
  dentinf   = rep(0, length(obs_times))
)

param_names <- c("sigSn","sigIn","sigF","sigP","f_Sn","rn","k_In","k_Sn",
                 "sigJn","probn","theta_Sn","theta_In","theta_P","theta_Jn","xi")
state_names <- c("Sn","In","Jn","error_count","F","T_Sn","T_In","P")

parameters <- mif.estimate[param_names]

pomplist <- list()
for (i in 1:8) {
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
names(pomplist) <- paste0("u", 1:8)

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
