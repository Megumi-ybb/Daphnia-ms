#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))

library(foreach)
library(doParallel)
registerDoParallel(cores = 100)
library(pomp)
library(panelPomp)

args <- commandArgs(trailingOnly = TRUE)
b <- as.integer(args[1])
alt_name <- args[2]
cat("Bootstrap LRT alt fit: dataset", b, ", alternative:", alt_name, "\n")

set.seed(801 * 1000 + b + 100000)

# ---- Map alternative name to parameter vector ----
alt_params_map <- list(
  xi       = c("xi"),
  theta_Sn = c("theta_Sn"),
  theta_In = c("theta_In"),
  theta_P  = c("theta_P"),
  probn    = c("probn"),
  rn       = c("rn"),
  f_Sn     = c("f_Sn")
)

if (!(alt_name %in% names(alt_params_map))) {
  stop("Unknown alternative: ", alt_name,
       ". Valid options: ", paste(names(alt_params_map), collapse = ", "))
}

specific_names <- alt_params_map[[alt_name]]
cat("Unit-specific parameters:", paste(specific_names, collapse = ", "), "\n")

# ---- Helper function ----
create_parameters <- function(parameter_names, parameters) {
  shared_parameter <- parameters[!(names(parameters) %in% parameter_names)]
  specific_values <- lapply(parameter_names, function(p) rep(parameters[[p]], 8))
  specific_mat <- matrix(
    data = unlist(specific_values),
    nrow = length(parameter_names),
    ncol = 8,
    byrow = TRUE,
    dimnames = list(parameter_names, paste0("u", 1:8))
  )
  return(list(shared_parameter = shared_parameter, specific_mat = specific_mat))
}

# ---- Load simulated data and true params ----
sim_data_list <- readRDS(sprintf("simulated_data/sim_data_%03d.rds", b))
true_params <- readRDS("simulated_data/true_params.rds")

# ---- Model specification (same C snippets as all_shared) ----
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

param_names <- c("sigSn","sigIn","sigF","sigP","f_Sn","rn","k_In","k_Sn",
                 "sigJn","probn","theta_Sn","theta_In","theta_P","theta_Jn","xi")
state_names <- c("Sn","In","Jn","error_count","F","T_Sn","T_In","P")

# ---- Build panelPomp from simulated data with unit-specific params ----
parameters <- true_params[param_names]

pomplist <- list()
for (i in 1:8) {
  u_name <- paste0("u", i)
  dat <- sim_data_list[[u_name]]
  colnames(dat) <- c("day","dentadult","dentinf")

  pomplist[[i]] <- pomp(
    data      = dat,
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

result_params <- create_parameters(specific_names, true_params)
panelfood <- panelPomp(pomplist, shared = result_params$shared_parameter,
                       specific = result_params$specific_mat)

# ---- Algorithmic parameters ----
Np <- 1000
Np_rep <- 10
Mp <- 1000

# ---- Round 1: mif2 from warm start ----
dent_rw.sd <- 0.05

parameter_candidates <- list(shared = result_params$shared_parameter,
                             specific = result_params$specific_mat)

mf1 <- foreach(
  i = 1:(3 * getDoParWorkers()),
  .packages = c("pomp", "panelPomp"),
  .inorder = FALSE,
  .options.multicore = list(set.seed = TRUE)
) %dopar% {
  mif2(
    panelfood,
    block = TRUE,   # block-resample unit-specific params; matches specific_model_block real-data fits
    Nmif = 150,
    shared.start = parameter_candidates$shared,
    specific.start = parameter_candidates$specific,
    rw.sd = rw_sd(xi=dent_rw.sd, sigSn=0, sigIn=dent_rw.sd,
                  sigF=dent_rw.sd, theta_Sn=dent_rw.sd, theta_In=dent_rw.sd,
                  theta_P=dent_rw.sd, k_Sn=dent_rw.sd, k_In=dent_rw.sd,
                  f_Sn=dent_rw.sd, rn=dent_rw.sd, probn=dent_rw.sd,
                  sigP=dent_rw.sd, sigJn=dent_rw.sd, theta_Jn=dent_rw.sd),
    cooling.type = "geometric",
    cooling.fraction.50 = 0.7,
    Np = Mp
  ) -> m1

  ll <- replicate(n = Np_rep, unitLogLik(pfilter(m1, Np = Np)))
  list(mif = m1, ll = panel_logmeanexp(x = ll, MARGIN = 1, se = TRUE))
}

# ---- Select top 25% and expand for round 2 ----
log_list <- sapply(mf1, function(x) x$ll[1])
select <- order(log_list, decreasing = TRUE)[1:ceiling(length(mf1)/4)]

shared_dataframe <- data.frame(t(mf1[[select[1]]]$mif@shared))
for (i in 1:length(select)) {
  shared_dataframe[i,] <- t(mf1[[select[i]]]$mif@shared)
}
shared_dataframe <- shared_dataframe[rep(1:nrow(shared_dataframe), each = 4), ]

specific_list <- lapply(select, function(idx) mf1[[idx]]$mif@specific)
replicated_specific_list <- list()
for (i in 1:length(specific_list)) {
  replicated_specific_list <- c(replicated_specific_list,
                                replicate(4, specific_list[[i]], simplify = FALSE))
}

# ---- Round 2 ----
dent_rw.sd <- 0.04

mf <- foreach(
  i = 1:(3 * getDoParWorkers()),
  .packages = c("pomp", "panelPomp"),
  .inorder = FALSE,
  .options.multicore = list(set.seed = TRUE)
) %dopar% {
  share_para_temp <- as.numeric(shared_dataframe[i,])
  names(share_para_temp) <- colnames(shared_dataframe)
  specific_para_temp <- as.matrix(replicated_specific_list[[i]])

  mif2(
    panelfood,
    block = TRUE,   # block-resample unit-specific params; matches specific_model_block real-data fits
    Nmif = 150,
    shared.start = share_para_temp,
    specific.start = specific_para_temp,
    rw.sd = rw_sd(xi=dent_rw.sd, sigSn=0, sigIn=dent_rw.sd,
                  sigF=dent_rw.sd, theta_Sn=dent_rw.sd, theta_In=dent_rw.sd,
                  theta_P=dent_rw.sd, k_Sn=dent_rw.sd, k_In=dent_rw.sd,
                  f_Sn=dent_rw.sd, rn=dent_rw.sd, probn=dent_rw.sd,
                  sigP=dent_rw.sd, sigJn=dent_rw.sd, theta_Jn=dent_rw.sd),
    cooling.type = "geometric",
    cooling.fraction.50 = 0.7,
    Np = Mp
  ) -> m1

  ll <- replicate(n = Np_rep, unitLogLik(pfilter(m1, Np = Np)))
  list(mif = m1, ll = panel_logmeanexp(x = ll, MARGIN = 1, se = TRUE))
}

# ---- Extract best result ----
lls <- matrix(unlist(sapply(mf, getElement, "ll")), nrow = 2)
best <- which.max(lls[1,])
result <- list(
  ll   = unname(mf[[best]]$ll[1]),
  se   = unname(mf[[best]]$ll[2]),
  coef = coef(mf[[best]]$mif),
  # provenance for the paired-settings guard in collect_lrt.R
  Np      = Np,
  Mp      = Mp,
  Np_rep  = Np_rep,
  block   = TRUE,
  Nmif    = c(round1 = 150, round2 = 150)
)

saveRDS(result, file = paste0('results_alt/lrt_',alt_name,"_",b,'.rds'))
cat("Done. ll =", result$ll, "se =", result$se, "\n")