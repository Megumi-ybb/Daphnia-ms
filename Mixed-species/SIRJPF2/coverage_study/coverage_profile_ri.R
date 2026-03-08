#!/apps/bin/Rscript
library(reshape2)
library(magrittr)
library(foreach)
library(doParallel)
library(pomp)
library(panelPomp)
library(tidyverse)

## ---------------------------------------------------------------------------
## 0. Read dataset index from command-line argument
## ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: ./coverage_profile_ri.R <dataset_index>\n  e.g. ./coverage_profile_ri.R 42")
}
b <- as.integer(args[1])
if (is.na(b) || b < 1 || b > 100) {
  stop("Dataset index must be an integer between 1 and 100.")
}

cat("=== Coverage profile ri: dataset", b, "===\n")
cat("Start time:", format(Sys.time()), "\n")

## Register all available cores
n_cores <- parallel::detectCores(logical = TRUE)
registerDoParallel(cores = n_cores)
cat("Using", n_cores, "cores\n")

## Random seed: reproducible but different per dataset
set.seed(801 * 1000 + b)

name_str <- "ri"
run_level <- 3

## ---------------------------------------------------------------------------
## 1. Load simulated dataset
## ---------------------------------------------------------------------------
sim_data_list <- readRDS(sprintf("simulated_data/sim_data_%03d.rds", b))
true_params   <- readRDS("simulated_data/true_params.rds")

## ---------------------------------------------------------------------------
## 2. Model specification (identical C snippets)
## ---------------------------------------------------------------------------

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

param_names <- c("xi","sigSn","sigIn","sigSi","sigIi","sigF","sigP",
                 "theta_Sn","theta_In","theta_Si","theta_P","theta_Ii",
                 "f_Sn","f_Si","rn","ri","probn","probi",
                 "k_Ii","k_In","k_Sn","k_Si","sigJi","sigJn","theta_Jn","theta_Ji")
state_names <- c("Sn","In","Si","Jn","Ji","Ii","error_count","F","T_Sn","T_In","T_Si","T_Ii","P")

## ---------------------------------------------------------------------------
## 3. Build panelPomp with simulated data
## ---------------------------------------------------------------------------

pomplist <- list()
for (i in 1:8) {
  unit_name <- paste0("u", i)
  dat <- sim_data_list[[unit_name]]
  colnames(dat) <- c("day","dentadult","dentinf","lumadult","luminf")

  pomplist[[i]] <- pomp(
    data      = dat,
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
  coef(pomplist[[i]]) <- true_params[param_names]
}
names(pomplist) <- paste0("u", 1:8)

panelfood <- panelPomp(pomplist, shared = true_params)

## ---------------------------------------------------------------------------
## 4. Profile design for ri (same structure as profile_ri.R)
## ---------------------------------------------------------------------------

generate_parameter_profile <- function(prof_name, nprof = 80) {
  shared_ub <- true_params * 10
  shared_lb <- shared_ub / 100

  ub_unit <- log(shared_ub[prof_name])
  lb_unit <- log(shared_lb[prof_name])

  shared_lb <- shared_lb[!(names(shared_lb) %in% c("sigSn","sigSi", prof_name))]
  shared_ub <- shared_ub[!(names(shared_ub) %in% c("sigSn","sigSi", prof_name))]

  parameter_shared <- pomp::profile_design(
    temp  = seq(lb_unit, ub_unit, length.out = nprof),
    lower = log(shared_lb),
    upper = log(shared_ub),
    type  = "runif",
    nprof = nprof
  )
  parameter_shared <- parameter_shared %>% rename(!!prof_name := temp)
  parameter_shared <- exp(parameter_shared)
  parameter_shared$sigSn <- 0
  parameter_shared$sigSi <- 0
  return(parameter_shared)
}

generate_sd <- function(x = 0.05, profile_name) {
  sd_list <- c(
    ri = x, rn = x, f_Si = x, f_Sn = x, probi = x, probn = x,
    xi = x, theta_Sn = x, theta_Si = x, theta_Ii = x, theta_In = x,
    theta_P = x, theta_Ji = x, theta_Jn = x,
    sigSn = 0, sigSi = 0, sigIn = x, sigIi = x,
    sigJi = x, sigJn = x, sigF = x, sigP = x,
    k_Ii = x, k_In = x, k_Si = x, k_Sn = x
  )
  sd_list[profile_name] <- 0
  return(sd_list)
}

parameter_shared <- generate_parameter_profile(name_str)

algorithmic.params <- list(
  Np     = c(50, 320, 1000),
  Np_rep = c( 2,  10,   20),
  Mp     = c(50, 400,  500),
  Nmif   = c( 2, 320,  250)
)

## ---------------------------------------------------------------------------
## 5. Round 1: initial mif2 with random starts
## ---------------------------------------------------------------------------

cat("Round 1 starting...\n")
dent_rw_sd_first <- generate_sd(x = 0.05, profile_name = name_str)

{
  foreach(
    i = 1:nrow(parameter_shared),
    .packages = c("pomp","panelPomp"),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE)
  ) %dopar% {
    guessed.parameter.values <- as.numeric(parameter_shared[i,])
    names(guessed.parameter.values) <- colnames(parameter_shared)

    m1 <- mif2(
      panelfood,
      Nmif = 200,
      shared.start = guessed.parameter.values,
      rw.sd = rw_sd(
        xi=dent_rw_sd_first['xi'],
        sigSn=dent_rw_sd_first['sigSn'], sigIn=dent_rw_sd_first['sigIn'],
        sigSi=dent_rw_sd_first['sigSi'], sigIi=dent_rw_sd_first['sigIi'],
        sigF=dent_rw_sd_first['sigF'],
        theta_Sn=dent_rw_sd_first['theta_Sn'], theta_In=dent_rw_sd_first['theta_In'],
        theta_Si=dent_rw_sd_first['theta_Si'], theta_P=dent_rw_sd_first['theta_P'],
        theta_Ii=dent_rw_sd_first['theta_Ii'],
        k_Sn=dent_rw_sd_first['k_Sn'], k_In=dent_rw_sd_first['k_In'],
        k_Si=dent_rw_sd_first['k_Si'], k_Ii=dent_rw_sd_first['k_Ii'],
        f_Sn=dent_rw_sd_first['f_Sn'], f_Si=dent_rw_sd_first['f_Si'],
        rn=dent_rw_sd_first['rn'], ri=dent_rw_sd_first['ri'],
        probn=dent_rw_sd_first['probn'], probi=dent_rw_sd_first['probi'],
        sigP=dent_rw_sd_first['sigP'],
        sigJi=dent_rw_sd_first['sigJi'], sigJn=dent_rw_sd_first['sigJn'],
        theta_Jn=dent_rw_sd_first['theta_Jn'], theta_Ji=dent_rw_sd_first['theta_Ji']
      ),
      cooling.type = "geometric",
      cooling.fraction.50 = 0.7,
      Np = algorithmic.params$Mp[run_level]
    )

    ll <- replicate(
      n = algorithmic.params$Np_rep[run_level],
      unitlogLik(pfilter(m1, Np = algorithmic.params$Np[run_level]))
    )

    list(
      mif_coef = coef(m1),
      ll       = panel_logmeanexp(x = ll, MARGIN = 1, se = TRUE)
    )
  }
} -> mf1

cat("Round 1 done.\n")

## ---------------------------------------------------------------------------
## 6. Select top 25% for Round 2
## ---------------------------------------------------------------------------

log_list <- sapply(mf1, function(x) x$ll[1])
select   <- order(log_list, decreasing = TRUE)[1:ceiling(length(mf1)/4)]

shared_dataframe <- data.frame(t(mf1[[select[1]]]$mif_coef))
for (i in seq_along(select)) {
  shared_dataframe[i,] <- t(mf1[[select[i]]]$mif_coef)
}
# Replicate each row 4 times to get back to nrow(parameter_shared) rows
shared_dataframe <- shared_dataframe[rep(1:nrow(shared_dataframe), each = 4), ]

## ---------------------------------------------------------------------------
## 7. Round 2: refined mif2
## ---------------------------------------------------------------------------

cat("Round 2 starting...\n")
dent_rw_sd_second <- generate_sd(x = 0.04, profile_name = name_str)

{
  foreach(
    i = 1:nrow(parameter_shared),
    .packages = c("pomp","panelPomp"),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE)
  ) %dopar% {
    share_para_temp <- as.numeric(shared_dataframe[i,])
    names(share_para_temp) <- colnames(shared_dataframe)

    m1 <- mif2(
      panelfood,
      Nmif = 300,
      shared.start = share_para_temp,
      rw.sd = rw_sd(
        xi=dent_rw_sd_second['xi'],
        sigSn=dent_rw_sd_second['sigSn'], sigIn=dent_rw_sd_second['sigIn'],
        sigSi=dent_rw_sd_second['sigSi'], sigIi=dent_rw_sd_second['sigIi'],
        sigF=dent_rw_sd_second['sigF'],
        theta_Sn=dent_rw_sd_second['theta_Sn'], theta_In=dent_rw_sd_second['theta_In'],
        theta_Si=dent_rw_sd_second['theta_Si'], theta_P=dent_rw_sd_second['theta_P'],
        theta_Ii=dent_rw_sd_second['theta_Ii'],
        k_Sn=dent_rw_sd_second['k_Sn'], k_In=dent_rw_sd_second['k_In'],
        k_Si=dent_rw_sd_second['k_Si'], k_Ii=dent_rw_sd_second['k_Ii'],
        f_Sn=dent_rw_sd_second['f_Sn'], f_Si=dent_rw_sd_second['f_Si'],
        rn=dent_rw_sd_second['rn'], ri=dent_rw_sd_second['ri'],
        probn=dent_rw_sd_second['probn'], probi=dent_rw_sd_second['probi'],
        sigP=dent_rw_sd_second['sigP'],
        sigJi=dent_rw_sd_second['sigJi'], sigJn=dent_rw_sd_second['sigJn'],
        theta_Jn=dent_rw_sd_second['theta_Jn'], theta_Ji=dent_rw_sd_second['theta_Ji']
      ),
      cooling.type = "geometric",
      cooling.fraction.50 = 0.7,
      Np = algorithmic.params$Mp[run_level]
    )

    ll <- replicate(
      n = algorithmic.params$Np_rep[run_level],
      unitlogLik(pfilter(m1, Np = algorithmic.params$Np[run_level]))
    )

    ## STORAGE OPTIMIZATION: save only coefficients + loglik, not the full mif object
    list(
      mif_coef = coef(m1),
      ll       = panel_logmeanexp(x = ll, MARGIN = 1, se = TRUE)
    )
  }
} -> mf

cat("Round 2 done.\n")

## ---------------------------------------------------------------------------
## 8. Extract subset_data (ri + loglik only) — the minimal data for mcap()
## ---------------------------------------------------------------------------

# Get final log-likelihoods
final_likes <- sapply(mf, function(x) x$ll[1])

# Build final_params: each row = final coefficients of one profile point
final_params <- data.frame(t(mf[[1]]$mif_coef))
for (i in 2:length(mf)) {
  final_params[i,] <- t(mf[[i]]$mif_coef)
}
final_params$loglik <- final_likes

# Group by ri, keep only the best loglik at each profile point
# (same pattern as Target_profiling_plots.R)
subset_data_ri <- final_params %>%
  group_by(ri) %>%
  filter(loglik == max(loglik)) %>%
  ungroup() %>%
  select(ri, loglik)

# Add log-transformed ri (needed by mcap)
subset_data_ri$log_ri <- log(subset_data_ri$ri)

## ---------------------------------------------------------------------------
## 9. Save minimal output
## ---------------------------------------------------------------------------

dir.create("coverage_results", showWarnings = FALSE, recursive = TRUE)

saveRDS(
  subset_data_ri,
  file = sprintf("coverage_results/profile_ri_%03d.rds", b)
)

cat("Saved subset_data_ri for dataset", b, "\n")
cat("  Rows:", nrow(subset_data_ri), "\n")
cat("  File size:", file.size(sprintf("coverage_results/profile_ri_%03d.rds", b)), "bytes\n")
cat("End time:", format(Sys.time()), "\n")
