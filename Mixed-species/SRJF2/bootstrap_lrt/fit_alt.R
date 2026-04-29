#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))

library(foreach)
library(doParallel)
registerDoParallel(cores = detectCores())
library(pomp)
library(panelPomp)

args <- commandArgs(trailingOnly = TRUE)
b <- as.integer(args[1])
alt_name <- args[2]
cat("Bootstrap LRT alt fit: dataset", b, ", alternative:", alt_name, "\n")

set.seed(801 * 1000 + b + 100000)

# ---- Map alternative name to parameter vector ----
alt_params_map <- list(
  f_Si_f_Sn         = c("f_Si", "f_Sn"),
  rn_ri             = c("rn", "ri"),
  theta_Sn_theta_Si = c("theta_Sn", "theta_Si")
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
  specific_values <- lapply(parameter_names, function(p) rep(parameters[[p]], 9))
  specific_mat <- matrix(
    data = unlist(specific_values),
    nrow = length(parameter_names),
    ncol = 9,
    byrow = TRUE,
    dimnames = list(parameter_names, paste0("u", 1:9))
  )
  return(list(shared_parameter = shared_parameter, specific_mat = specific_mat))
}

# ---- Load simulated data and true params ----
sim_data_list <- readRDS(sprintf("simulated_data/sim_data_%03d.rds", b))
true_params <- readRDS("simulated_data/true_params.rds")

# ---- Model specification (same C snippets as all_shared) ----
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

# ---- Build panelPomp from simulated data with unit-specific params ----
parameters <- true_params[param_names]

pomplist <- list()
for (i in 1:9) {
  u_name <- paste0("u", i)
  dat <- sim_data_list[[u_name]]
  colnames(dat) <- c("day","dentadult","lumadult")

  pomplist[[i]] <- pomp(
    data      = dat,
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
    Nmif = 150,
    shared.start = parameter_candidates$shared,
    specific.start = parameter_candidates$specific,
    rw.sd = rw_sd(sigSn=0, sigSi=0, sigF=dent_rw.sd,
                  theta_Sn=dent_rw.sd, theta_Si=dent_rw.sd,
                  k_Sn=dent_rw.sd, k_Si=dent_rw.sd,
                  f_Sn=dent_rw.sd, f_Si=dent_rw.sd,
                  rn=dent_rw.sd, ri=dent_rw.sd,
                  sigJi=dent_rw.sd, sigJn=dent_rw.sd,
                  theta_Jn=dent_rw.sd, theta_Ji=dent_rw.sd),
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
    Nmif = 150,
    shared.start = share_para_temp,
    specific.start = specific_para_temp,
    rw.sd = rw_sd(sigSn=0, sigSi=0, sigF=dent_rw.sd,
                  theta_Sn=dent_rw.sd, theta_Si=dent_rw.sd,
                  k_Sn=dent_rw.sd, k_Si=dent_rw.sd,
                  f_Sn=dent_rw.sd, f_Si=dent_rw.sd,
                  rn=dent_rw.sd, ri=dent_rw.sd,
                  sigJi=dent_rw.sd, sigJn=dent_rw.sd,
                  theta_Jn=dent_rw.sd, theta_Ji=dent_rw.sd),
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
  coef = coef(mf[[best]]$mif)
)

saveRDS(result, file = paste0('results_alt/lrt_',alt_name,"_",b,'.rds'))
cat("Done. ll =", result$ll, "se =", result$se, "\n")
