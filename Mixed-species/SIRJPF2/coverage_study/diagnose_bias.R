#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(magrittr)
library(foreach)
library(doParallel)
library(pomp)
library(panelPomp)
library(tidyverse)

## ---------------------------------------------------------------------------
## Diagnostic: warm-start mif2 at the truth on a single simulated dataset.
##
## Question: does mif2 stay near the truth, or does it drift toward a ridge?
##   - If it stays  -> coverage misses are caused by poor global search
##                     (random starts in the profile don't reach the truth).
##   - If it drifts -> either the simulated dataset's actual MLE is genuinely
##                     off-truth (finite-sample bias) or mif2 is descending a
##                     ridge. The pfilter(truth) loglik tells these apart.
##
## Usage: ./diagnose_bias.R <dataset_index>
## Output: diagnose_results/diagnose_<b>.rds
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: ./diagnose_bias.R <dataset_index>\n  e.g. ./diagnose_bias.R 42")
}
b <- as.integer(args[1])
if (is.na(b) || b < 1 || b > 100) {
  stop("Dataset index must be an integer between 1 and 100.")
}

cat("=== Diagnose bias: dataset", b, "===\n")
cat("Start time:", format(Sys.time()), "\n")

n_cores <- parallel::detectCores(logical = TRUE)
registerDoParallel(cores = n_cores)
cat("Using", n_cores, "cores\n")

set.seed(900 * 1000 + b)

run_level <- 3
n_reps    <- 20  # mif2 chains warm-started at truth on this dataset

## ---------------------------------------------------------------------------
## 1. Load simulated dataset and truth
## ---------------------------------------------------------------------------
sim_data_list <- readRDS(sprintf("simulated_data/sim_data_%03d.rds", b))
true_params   <- readRDS("simulated_data/true_params.rds")

## ---------------------------------------------------------------------------
## 2. Model specification (identical C snippets to the rest of coverage_study)
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
## 3. Build panelPomp with the simulated data, parameters set to truth
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

algorithmic.params <- list(
  Np     = c(50, 320, 1000),
  Np_rep = c( 2,  10,   20),
  Mp     = c(50, 400,  500),
  Nmif   = c( 2, 320,  250)
)

## ---------------------------------------------------------------------------
## 4. Reference loglik at truth (no mif2): tells us where the truth sits
##    on the likelihood surface for THIS simulated dataset.
## ---------------------------------------------------------------------------

cat("Computing pfilter loglik at truth (no mif2)...\n")

ll_truth_mat <- foreach(
  rep = 1:algorithmic.params$Np_rep[run_level],
  .packages = c("pomp","panelPomp"),
  .combine  = cbind,
  .inorder  = FALSE,
  .options.multicore = list(set.seed = TRUE)
) %dopar% {
  unitLogLik(pfilter(panelfood, Np = algorithmic.params$Np[run_level]))
}

ll_truth <- panel_logmeanexp(x = ll_truth_mat, MARGIN = 1, se = TRUE)
cat(sprintf("  ll(truth) = %.3f  (SE %.3f)\n", ll_truth[1], ll_truth[2]))

## ---------------------------------------------------------------------------
## 5. mif2 chains warm-started at truth (full parameter space free)
##
## Same Round 1 + Round 2 schedule used by coverage_profile_*.R, except
## NO parameter is pinned. Starting point is the truth, perturbed by rw_sd.
## ---------------------------------------------------------------------------

build_rw <- function(x) {
  c(
    ri = x, rn = x, f_Si = x, f_Sn = x, probi = x, probn = x,
    xi = x, theta_Sn = x, theta_Si = x, theta_Ii = x, theta_In = x,
    theta_P = x, theta_Ji = x, theta_Jn = x,
    sigSn = 0, sigSi = 0, sigIn = x, sigIi = x,
    sigJi = x, sigJn = x, sigF = x, sigP = x,
    k_Ii = x, k_In = x, k_Si = x, k_Sn = x
  )
}

rw1 <- build_rw(0.05)
rw2 <- build_rw(0.04)

cat("Running", n_reps, "mif2 chains warm-started at truth...\n")

mf <- foreach(
  rep = 1:n_reps,
  .packages = c("pomp","panelPomp"),
  .inorder  = FALSE,
  .options.multicore = list(set.seed = TRUE)
) %dopar% {

  m1 <- mif2(
    panelfood,
    Nmif = 200,
    shared.start = true_params,
    rw.sd = rw_sd(
      xi=rw1['xi'],
      sigSn=rw1['sigSn'], sigIn=rw1['sigIn'],
      sigSi=rw1['sigSi'], sigIi=rw1['sigIi'],
      sigF=rw1['sigF'],
      theta_Sn=rw1['theta_Sn'], theta_In=rw1['theta_In'],
      theta_Si=rw1['theta_Si'], theta_P=rw1['theta_P'],
      theta_Ii=rw1['theta_Ii'],
      k_Sn=rw1['k_Sn'], k_In=rw1['k_In'],
      k_Si=rw1['k_Si'], k_Ii=rw1['k_Ii'],
      f_Sn=rw1['f_Sn'], f_Si=rw1['f_Si'],
      rn=rw1['rn'], ri=rw1['ri'],
      probn=rw1['probn'], probi=rw1['probi'],
      sigP=rw1['sigP'],
      sigJi=rw1['sigJi'], sigJn=rw1['sigJn'],
      theta_Jn=rw1['theta_Jn'], theta_Ji=rw1['theta_Ji']
    ),
    cooling.type = "geometric",
    cooling.fraction.50 = 0.7,
    Np = algorithmic.params$Mp[run_level]
  )

  m2 <- mif2(
    m1,
    Nmif = 300,
    rw.sd = rw_sd(
      xi=rw2['xi'],
      sigSn=rw2['sigSn'], sigIn=rw2['sigIn'],
      sigSi=rw2['sigSi'], sigIi=rw2['sigIi'],
      sigF=rw2['sigF'],
      theta_Sn=rw2['theta_Sn'], theta_In=rw2['theta_In'],
      theta_Si=rw2['theta_Si'], theta_P=rw2['theta_P'],
      theta_Ii=rw2['theta_Ii'],
      k_Sn=rw2['k_Sn'], k_In=rw2['k_In'],
      k_Si=rw2['k_Si'], k_Ii=rw2['k_Ii'],
      f_Sn=rw2['f_Sn'], f_Si=rw2['f_Si'],
      rn=rw2['rn'], ri=rw2['ri'],
      probn=rw2['probn'], probi=rw2['probi'],
      sigP=rw2['sigP'],
      sigJi=rw2['sigJi'], sigJn=rw2['sigJn'],
      theta_Jn=rw2['theta_Jn'], theta_Ji=rw2['theta_Ji']
    ),
    cooling.type = "geometric",
    cooling.fraction.50 = 0.7
  )

  ll <- replicate(
    n = algorithmic.params$Np_rep[run_level],
    unitLogLik(pfilter(m2, Np = algorithmic.params$Np[run_level]))
  )

  list(
    coef = coef(m2),
    ll   = panel_logmeanexp(x = ll, MARGIN = 1, se = TRUE)
  )
}

cat("All", n_reps, "chains done.\n")

## ---------------------------------------------------------------------------
## 6. Save results
## ---------------------------------------------------------------------------

all_coefs <- do.call(rbind, lapply(mf, function(x) x$coef))
all_lls   <- sapply(mf, function(x) x$ll[1])
all_ses   <- sapply(mf, function(x) x$ll[2])

result <- list(
  b           = b,
  true_params = true_params,
  ll_truth    = ll_truth,           # c(estimate, se)
  fits        = data.frame(all_coefs, loglik = all_lls, loglik_se = all_ses)
)

dir.create("diagnose_results", showWarnings = FALSE, recursive = TRUE)
saveRDS(result, file = sprintf("diagnose_results/diagnose_%03d.rds", b))

cat(sprintf("Saved %d mif2 fits + ll(truth)=%.3f for dataset %d\n",
            n_reps, ll_truth[1], b))
cat("End time:", format(Sys.time()), "\n")
