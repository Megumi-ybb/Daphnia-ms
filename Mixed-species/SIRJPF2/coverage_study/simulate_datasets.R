#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(pomp)
library(panelPomp)

set.seed(801)

B <- 100

mif.estimate <- c(
  ri        = 2.152877e+05,
  rn        = 4.082784e+01,
  f_Si      = 2.419416e-07,
  f_Sn      = 1.100347e-03,
  probi     = 1.341780e+03,
  probn     = 2.715528e-01,
  xi        = 2.223829e+01,
  theta_Sn  = 8.478459e-04,
  theta_Si  = 2.524539e-03,
  theta_Ii  = 3.854897e-01,
  theta_In  = 5.837784e-01,
  theta_P   = 9.477111e-04,
  theta_Ji  = 5.620201e-04,
  theta_Jn  = 1.868651e-05,
  sigSn     = 0.000000e+00,
  sigSi     = 0.000000e+00,
  sigIn     = 2.930410e-04,
  sigIi     = 1.727540e-07,
  sigJi     = 3.021246e-01,
  sigJn     = 2.840041e-01,
  sigF      = 1.436943e-01,
  sigP      = 2.714232e-01,
  k_Ii      = 1.387153e+00,
  k_In      = 9.020138e-01,
  k_Si      = 5.262009e+00,
  k_Sn      = 4.103463e+00
)

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

  T_Sn = fabs(Sn);
  T_In = fabs(In);
  T_Si = fabs(Si);
  T_Ii = fabs(Ii);
")

dyn_init <- Csnippet("
  Sn = 2.333;
  Si = 0.667;
  F  = 16.667;
  Jn = 0; Ji = 0;
  T_Sn = 0.0; T_Si = 0.0; T_In = 0.0; T_Ii = 0.0;
  In = 0.0; Ii = 0.0;
  error_count = 0.0;
  P = 0;
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

rmeas <- Csnippet("
  dentadult = rnbinom_mu(k_Sn, T_Sn);
  dentinf   = rnbinom_mu(k_In, T_In);
  lumadult  = rnbinom_mu(k_Si, T_Si);
  luminf    = rnbinom_mu(k_Ii, T_Ii);
")

pt <- parameter_trans(
  log = c("sigSn","sigIn","sigSi","sigIi","sigF","sigP","f_Sn","f_Si",
          "rn","ri","k_Ii","k_In","k_Sn","k_Si","sigJi","sigJn","probn","probi",
          "theta_Sn","theta_In","theta_Si","theta_Ii","theta_P","theta_Jn","theta_Ji","xi")
)

obs_times <- (0:9) * 5 + 7

template_data <- data.frame(
  day       = obs_times,
  dentadult = rep(0, length(obs_times)),
  dentinf   = rep(0, length(obs_times)),
  lumadult  = rep(0, length(obs_times)),
  luminf    = rep(0, length(obs_times))
)

param_names <- c("xi","sigSn","sigIn","sigSi","sigIi","sigF","sigP",
                 "theta_Sn","theta_In","theta_Si","theta_P","theta_Ii",
                 "f_Sn","f_Si","rn","ri","probn","probi",
                 "k_Ii","k_In","k_Sn","k_Si","sigJi","sigJn","theta_Jn","theta_Ji")
state_names <- c("Sn","In","Si","Jn","Ji","Ii","error_count","F","T_Sn","T_In","T_Si","T_Ii","P")

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
    obsnames  = c("dentadult","dentinf","lumadult","luminf"),
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
    sim_data_list[[u]] <- sim[, c("day","dentadult","dentinf","lumadult","luminf")]
  }
  saveRDS(sim_data_list, file = sprintf("simulated_data/sim_data_%03d.rds", b))
  cat("  Dataset", b, "saved.\n")
}

saveRDS(mif.estimate, file = "simulated_data/true_params.rds")

cat("Done. All", B, "simulated datasets saved in simulated_data/\n")
cat("True log(ri) =", log(mif.estimate["ri"]), "\n")
