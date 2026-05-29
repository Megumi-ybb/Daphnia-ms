#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))

# ============================================================================
# No-F ablation of the all-shared SIRJPF2 model.
#
# Addresses the editor's request for "results without F" -- i.e. evidence that
# the depletable food compartment is load-bearing for the observed
# peak-then-collapse dynamics.
#
# Ablation relative to model/all_shared.R:
#   * F is HELD at its initial value (16.667) for the whole trajectory:
#     - F_term and the `F += F_term` update are removed,
#     - the F process noise (noiF / sigF) is removed,
#     - the F positivity/upper-bound check is removed (F is now constant).
#   * Daphnia recruitment still reads F (`rn*f_Sn*F*Sn`, `ri*f_Si*F*Si`) but F
#     no longer depletes, so birth rates no longer respond to a depletable
#     resource. F is kept explicit (not folded into rn) so the all-shared MLE
#     transfers directly as a warm start and rn stays on the same scale.
#
# Everything else -- data, measurement model, iterated-filtering protocol
# (3 rounds, Nmif 150/150/400, Np=Mp=1000, 10x workers chains, cooling 0.7) --
# is identical to model/all_shared.R, so the AIC comparison is directly
# interpretable. The two models are NOT nested (freezing F dynamics is a
# structural change, not a parameter restriction), so we compare via AIC, not
# a likelihood-ratio test.
#
# Parameter count: all-shared has k = 24 free parameters (sigSn = sigSi = 0 are
# boundary-fixed and not counted). Dropping sigF gives k = 23 here.
#
# Usage:  ./fit_no_F.R [rep]        (rep defaults to 1)
#   Matches the original single-job all-shared fit. Pass distinct rep indices
#   if you want independent global searches; collect_no_F.R takes the max ll.
# ============================================================================

library(reshape2)
library(magrittr)
library(foreach)
library(readxl)
library(doParallel)
registerDoParallel(cores = 100)
library(pomp)
library(panelPomp)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
cat("No-F ablation fit, rep", rep_id, "\n")

set.seed(917 * 1000 + rep_id)

run_level <- 3

# ---- Real data (identical ingestion to model/all_shared.R) ----
# NOTE: model/all_shared.R reads "./Mesocosmdata.xlsx" sheet 3. The repo root
# currently carries Mesocosmdata.xls; ensure the .xlsx used for the production
# all-shared fit is reachable at the path below on the HPC node.
Mesocosm_data <- read_excel("../../../Mesocosmdata.xls", 3)

dentNoPara <- Mesocosm_data[91:170, ]
dentNoPara <- subset(dentNoPara, select = c(rep, day, dent.adult, dent.inf, lum.adult, lum.adult.inf))
dentNoPara <- dentNoPara[80:1, ]
dentNoPara$day <- (dentNoPara$day - 1) * 5 + 7

data <- list()
trails <- c("K", "L", "M", "N", "O", "P", "Q", "S")
for (i in 1:length(trails)) {
  data[[i]] <- subset(dentNoPara,
                      select = c("day", "dent.adult", "dent.inf", "lum.adult", "lum.adult.inf"),
                      dentNoPara$rep == trails[i])
}

# ---- Process model: F frozen at its initial value ----
dyn_rpro <- Csnippet("
                      double Sn_term, In_term, P_term , Si_term, Ii_term, Jn_term, Ji_term;
                      double noiSn, noiIn, noiSi , noiIi , noiP, noiJn, noiJi;
                      double delta = 0.013; //fraction of volume replaced day-1

                      noiSn = rnorm(0, sigSn * sqrt(dt));
                      noiIn = rnorm(0, sigIn * sqrt(dt));
                      noiSi = rnorm(0, sigSi * sqrt(dt));
                      noiIi = rnorm(0, sigIi * sqrt(dt));
                      noiJn = rnorm(0, sigJn * sqrt(dt));
                      noiJi = rnorm(0, sigJi * sqrt(dt));
                      noiP = rnorm(0, sigP * sqrt(dt));

                      //------------Sn-------------
                      Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt -  probn * f_Sn * Sn * P * dt - delta * Sn * dt + Sn * noiSn;

                      //-----------Jn------------- (F held constant: no resource depletion feedback)
                      Jn_term = rn * f_Sn * F * Sn  * dt  -  0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;

                      //------------In--------------
                      In_term = probn * f_Sn * Sn * P * dt - theta_In * In *dt - delta * In * dt + In * noiIn;

                      //------------Si-------------
                      Si_term = 0.1 * Ji * dt - theta_Si * Si * dt -  probi * f_Si * Si * P * dt - delta * Si * dt + Si * noiSi;

                      //-----------Ji------------- (F held constant: no resource depletion feedback)
                      Ji_term = ri * f_Si *F * Si  * dt - 0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi;

                      //------------Ii--------------
                      Ii_term = probi * f_Si * Si * P * dt - theta_Ii * Ii *dt - delta * Ii * dt + Ii * noiIi;

                      //-----------F--------------- ablated: F is NOT updated, stays at its initial value

                      //----------P---------------
                      P_term = 30 * theta_In * In * dt + 30 * theta_Ii * Ii * dt - f_Sn * (Sn + xi * In) * P * dt - f_Si * (Si + xi * Ii) * P * dt- theta_P * P * dt - delta * P * dt + P * noiP;

                      Sn += Sn_term;
                      In += In_term;
                      Si += Si_term;
                      Ji += Ji_term;
                      Jn += Jn_term;
                      Ii += Ii_term;
                      P += P_term;

                      //Trace time
                     if (t - 4.0 < 0.001 && t - 4.0 > -0.001){
                      //Initial statement
                        P += 25;
                      }

                      if (Sn < 0.0 || Sn > 1e5) {
                        Sn = 0.0;
                        error_count += 1;
                      }
                      if (Si < 0.0 || Si > 1e5) {
                        Si = 0.0;
                        error_count += 1000000;
                      }
                      if (In < 0.0 || In > 1e5) {
                        In = 0.0;
                        error_count += 0.001;
                      }
                      if (Ii < 0.0 || Ii > 1e5) {
                        Ii = 0.0;
                        error_count += 0.000000001;
                      }
                      if (Jn < 0.0 || Jn > 1e5) {
                        Jn = 0.0;
                        error_count += 0.001;
                      }
                      if (Ji < 0.0 || Ji > 1e5) {
                        Ji = 0.0;
                        error_count += 0.000000001;
                      }
                      if ((P < 0.0 || P > 1e20)&& t > 3.9) {
                        P = 0.0;
                        error_count += 0.000001;
                      }

                      T_Sn = fabs(Sn);
                      T_In = fabs(In);
                      T_Si = fabs(Si);
                      T_Ii = fabs(Ii);
                      ")

# Initial state. Assume t0 = day 4. F set here and never changed thereafter.
dyn_init <- Csnippet("
                     Sn = 2.333; //2.3333 = 35/15L
                     Si = 0.667; //0.667 = 10/15L
                     F = 16.667;
                     Jn = 0;
                     Ji = 0;

                     T_Sn = 0.0;
                     T_Si = 0.0;
                     T_In = 0.0;
                     T_Ii = 0.0;
                     In = 0.0;
                     Ii = 0.0;

                     error_count = 0.0;
                     P = 0;
                     ")

dmeas <- Csnippet("
                 if (error_count > 0.0) {
                   lik = -150;
                  }
                 else{
                    if(give_log){
                    lik  = dnbinom_mu(lumadult,k_Si,T_Si,give_log) +  dnbinom_mu(luminf,k_Ii,T_Ii,give_log) + dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) +  dnbinom_mu(dentinf,k_In,T_In,give_log);
                    }
                    else{
                    lik  = dnbinom_mu(lumadult,k_Si,T_Si,give_log) *  dnbinom_mu(luminf,k_Ii,T_Ii,give_log) * dnbinom_mu(dentadult,k_Sn,T_Sn,give_log) *  dnbinom_mu(dentinf,k_In,T_In,give_log);
                    }
                      }
                 ")

rmeas <- Csnippet("
                 dentadult = rnbinom_mu(k_Sn,T_Sn);
                 dentinf = rnbinom_mu(k_In,T_In);
                 lumadult = rnbinom_mu(k_Si,T_Si);
                 luminf = rnbinom_mu(k_Ii,T_Ii);
                 ")

# sigF dropped from the estimation-scale log transform.
pt <- parameter_trans(
  log = c( "sigSn", "sigIn", "sigSi", "sigIi", "sigP", "f_Sn", "f_Si",
           "rn", "ri", "k_Ii", "k_In", "k_Sn", "k_Si", "sigJi", "sigJn", "probn", "probi",
           "theta_Sn", "theta_In", "theta_Si", "theta_Ii", "theta_P", "theta_Jn", "theta_Ji", "xi")
)

# Parameter vector: identical to all_shared minus sigF (25 names; 23 free once
# the boundary-fixed sigSn = sigSi = 0 are excluded).
parameters <- c(28.6562, 0, 0.0003063207, 0, 0.02208698, 0.238589, 0.1479834,
                0.5489315, 0.0318604, 0.02024991, 0.3531879, 0.001105668, 1.838259e-05,
                59.04676, 13076.83, 0.2565626, 31.10083, 1.241092, 1.005756, 4.282648,
                4.715556, 0.2727418, 0.2836891, 0.0001532613, 0.0001299562)
names(parameters) <- c("xi", "sigSn", "sigIn", "sigSi", "sigIi", "sigP", "theta_Sn",
                       "theta_In", "theta_Si", "theta_P", "theta_Ii", "f_Sn", "f_Si", "rn", "ri", "probn",
                       "probi", "k_Ii", "k_In", "k_Sn", "k_Si", "sigJi", "sigJn", "theta_Jn", "theta_Ji")

param_names <- c("xi", "sigSn", "sigIn", "sigSi", "sigIi", "sigP", "theta_Sn", "theta_In", "theta_Si", "theta_P",
                 "theta_Ii", "f_Sn", "f_Si", "rn", "ri", "probn", "probi", "k_Ii", "k_In", "k_Sn", "k_Si", "sigJi", "sigJn", "theta_Jn", "theta_Ji")
# F remains a state variable so it is initialized to 16.667 and read by the
# recruitment terms; it is simply never updated.
state_names <- c("Sn", "In", "Si", "Jn", "Ji", "Ii", "error_count", "F", "T_Sn", "T_In", "T_Si", "T_Ii", "P")

pomplist <- list()
for (i in 1:8) {
  colnames(data[[i]]) <- c('day', 'dentadult', 'dentinf', 'lumadult', 'luminf')
  pomp(data = data[[i]],
       times = "day",
       t0 = 1,
       rprocess = euler(dyn_rpro, delta.t = 1/4),
       rinit = dyn_init,
       dmeasure = dmeas,
       rmeasure = rmeas,
       partrans = pt,
       obsnames = c("dentadult", "dentinf", "lumadult", "luminf"),
       accumvars = c("error_count"),
       paramnames = param_names,
       statenames = state_names
  ) -> pomplist[[i]]
  coef(pomplist[[i]]) <- parameters
}
names(pomplist) <- paste("u", 1:8, sep = "")

# All-shared MLE used as the warm start, with sigF removed.
shared_parameter <- c(
  ri = 13076.83, rn = 59.04676, f_Si = 1.838259e-05, f_Sn = 0.001105668, probi = 31.10083,
  probn = 0.2565626, xi = 28.6562, theta_Sn = 0.1479834, theta_Si = 0.0318604,
  theta_Ii = 0.3531879, theta_In = 0.5489315, theta_P = 0.02024991, theta_Ji = 0.0001299562,
  theta_Jn = 0.0001532613, sigSn = 0, sigSi = 0, sigIn = 0.0003063207,
  sigIi = 0.02208698, sigJi = 0.2727418, sigJn = 0.2836891, sigP = 0.238589,
  k_Ii = 1.241092, k_In = 1.005756, k_Si = 4.715556, k_Sn = 4.282648
)

panelfood <- panelPomp(pomplist, shared = shared_parameter)

algorithmic.params <- list(
  Np =     c(50, 500, 1000),
  Np_rep = c( 2,  10,  10),
  Mp =     c(50, 500, 1000),
  Nmif =   c( 2,  320, 250)
)

parameter_candidates <- list(shared_parameter)
names(parameter_candidates) <- c("shared")

# rw.sd: identical to all_shared minus sigF. sigSn = sigSi = 0 stay boundary-fixed.
no_F_rw_sd <- function(s) {
  rw_sd(xi = s, sigSn = 0, sigIn = s, sigSi = 0, sigIi = s,
        theta_Sn = s, theta_In = s, theta_Si = s, theta_P = s, theta_Ii = s,
        k_Sn = s, k_In = s, k_Si = s, k_Ii = s,
        f_Sn = s, f_Si = s, rn = s, ri = s, probn = s, probi = s,
        sigP = s, sigJi = s, sigJn = s, theta_Jn = s, theta_Ji = s)
}

# ---- Round 1: warm start from the all-shared MLE ----
dent_rw.sd <- 0.05
{
  foreach(
    i = 1:(10 * getDoParWorkers()),
    .packages = c("pomp", "panelPomp"),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE)
  ) %dopar%
    {
      guessed.parameter.values <- parameter_candidates
      mif2(
        panelfood,
        Nmif = 150,
        shared.start = guessed.parameter.values$shared,
        rw.sd = no_F_rw_sd(dent_rw.sd),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
      ) -> m1

      ll <- replicate(n = algorithmic.params$Np_rep[run_level],
                      unitLogLik(pfilter(m1, Np = algorithmic.params$Np[run_level])))

      list(mif = m1, ll = panel_logmeanexp(x = ll, MARGIN = 1, se = TRUE))
    }
} -> mf1

log_list <- c()
for (i in 1:length(mf1)) log_list <- c(log_list, mf1[[i]]$ll[1])
select <- order(log_list, decreasing = TRUE)[1:ceiling(length(mf1)/4)]
shared_dataframe <- data.frame(t(mf1[[select[1]]]$mif@shared))
for (i in 1:length(select)) shared_dataframe[i, ] <- t(mf1[[select[i]]]$mif@shared)
shared_dataframe <- shared_dataframe[rep(1:nrow(shared_dataframe), each = 4), ]

# ---- Round 2 ----
dent_rw.sd <- 0.04
{
  foreach(
    i = 1:(10 * getDoParWorkers()),
    .packages = c("pomp", "panelPomp"),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE)
  ) %dopar%
    {
      share_para_temp <- as.numeric(shared_dataframe[i, ])
      names(share_para_temp) <- colnames(shared_dataframe)

      mif2(
        panelfood,
        Nmif = 150,
        shared.start = share_para_temp,
        rw.sd = no_F_rw_sd(dent_rw.sd),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
      ) -> m1

      ll <- replicate(n = algorithmic.params$Np_rep[run_level],
                      unitLogLik(pfilter(m1, Np = algorithmic.params$Np[run_level])))

      list(mif = m1, ll = panel_logmeanexp(x = ll, MARGIN = 1, se = TRUE))
    }
} -> mf2

log_list <- c()
for (i in 1:length(mf2)) log_list <- c(log_list, mf2[[i]]$ll[1])
select <- order(log_list, decreasing = TRUE)[1:ceiling(length(mf2)/4)]
shared_dataframe <- data.frame(t(mf2[[select[1]]]$mif@shared))
for (i in 1:length(select)) shared_dataframe[i, ] <- t(mf2[[select[i]]]$mif@shared)
shared_dataframe <- shared_dataframe[rep(1:nrow(shared_dataframe), each = 4), ]

# ---- Round 3 ----
dent_rw.sd <- 0.04
{
  foreach(
    i = 1:(10 * getDoParWorkers()),
    .packages = c("pomp", "panelPomp"),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE)
  ) %dopar%
    {
      share_para_temp <- as.numeric(shared_dataframe[i, ])
      names(share_para_temp) <- colnames(shared_dataframe)

      mif2(
        panelfood,
        Nmif = 400,
        shared.start = share_para_temp,
        rw.sd = no_F_rw_sd(dent_rw.sd),
        cooling.type = "geometric",
        cooling.fraction.50 = 0.7,
        Np = algorithmic.params$Mp[run_level]
      ) -> m1

      ll <- replicate(n = algorithmic.params$Np_rep[run_level],
                      unitLogLik(pfilter(m1, Np = algorithmic.params$Np[run_level])))

      list(mif = m1, ll = panel_logmeanexp(x = ll, MARGIN = 1, se = TRUE))
    }
} -> mf

lls <- matrix(unlist(sapply(mf, getElement, "ll")), nrow = 2)
best <- which.max(lls[1, ])
mif.estimate <- coef(mf[[best]]$mif)
pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[1])
s.e.of.pf.loglik.of.mif.estimate <- unname(mf[[best]]$ll[2])

cat("Best ll =", pf.loglik.of.mif.estimate, "se =", s.e.of.pf.loglik.of.mif.estimate, "\n")

# Flat-save convention (matches bootstrap_lrt template): one result list per rep.
result <- list(
  ll   = pf.loglik.of.mif.estimate,
  se   = s.e.of.pf.loglik.of.mif.estimate,
  coef = mif.estimate
)
saveRDS(result, file = paste0('results/no_F_', rep_id, '.rds'))
cat("Saved results/no_F_", rep_id, ".rds\n", sep = "")
