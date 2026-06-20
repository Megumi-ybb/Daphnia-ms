#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(pomp)
library(panelPomp)

set.seed(801)

B <- 100
max_attempts <- 500   # draw this many candidates; keep the first B that pass the abnormality bounds

# ---- Load MLE from model fit ----
load("../model/best_result.rda")
cat("Loaded MLE parameters:\n")
print(mif.estimate)

# ---- Model specification (from model/all_shared.R) ----
dyn_rpro <- Csnippet("
  double Sn_term, F_term, Jn_term;
  double noiSn, noiF, noiJn;
  double delta = 0.013;

  noiSn = rnorm(0, sigSn * sqrt(dt));
  noiJn = rnorm(0, sigJn * sqrt(dt));
  noiF  = rnorm(0, sigF  * sqrt(dt));

  Sn_term = 0.1 * Jn * dt - theta_Sn * Sn * dt - delta * Sn * dt + Sn * noiSn;
  Jn_term = rn * f_Sn * F * Sn * dt - 0.1 * Jn * dt - theta_Jn * Jn * dt - delta * Jn * dt + Jn * noiJn;
  F_term  = F * noiF - f_Sn * F * (Sn + 1 * Jn) * dt - delta * F * dt + 0.37 * dt;

  F  += F_term;
  Sn += Sn_term;
  Jn += Jn_term;

  if (Sn < 0.0 || Sn > 1e5) { Sn = 0.0; error_count += 1; }
  if (F  < 0.0 || F  > 1e20){ F  = 0.0; error_count += 1000; }
  if (Jn < 0.0 || Jn > 1e5) { Jn = 0.0; error_count += 0.001; }

  T_Sn = fabs(Sn);
")

dyn_init <- Csnippet("
  Sn = 3;
  F  = 16.667;
  Jn = 0;
  T_Sn = 0.0;
  error_count = 0.0;
")

dmeas <- Csnippet("
  if (error_count > 0.0) {
    lik = -150;
  } else {
    if (give_log) {
      lik = dnbinom_mu(dentadult, k_Sn, T_Sn, give_log);
    } else {
      lik = dnbinom_mu(dentadult, k_Sn, T_Sn, give_log);
    }
  }
")

rmeas <- Csnippet("
  dentadult = rnbinom_mu(k_Sn, T_Sn);
")

pt <- parameter_trans(
  log = c("sigSn", "sigF", "f_Sn", "rn", "k_Sn", "sigJn",
          "theta_Sn", "theta_Jn")
)

param_names <- c("f_Sn", "sigF", "sigSn", "rn", "k_Sn", "sigJn",
                 "theta_Sn", "theta_Jn")
state_names <- c("Sn", "Jn", "error_count", "F", "T_Sn")

# ---- Observation times (same as data: day 7, 12, 17, ..., 52) ----
obs_times <- (0:9) * 5 + 7

template_data <- data.frame(
  day       = obs_times,
  dentadult = rep(0, length(obs_times))
)

# ---- Build panelPomp with template data ----
parameters <- mif.estimate[param_names]

pomplist <- list()
for (i in 1:10) {
  pomplist[[i]] <- pomp(
    data      = template_data,
    times     = "day",
    t0        = 1,
    rprocess  = euler(dyn_rpro, delta.t = 1/4),
    rinit     = dyn_init,
    dmeasure  = dmeas,
    rmeasure  = rmeas,
    partrans  = pt,
    obsnames  = c("dentadult"),
    accumvars = c("error_count"),
    paramnames = param_names,
    statenames = state_names
  )
  coef(pomplist[[i]]) <- parameters
}
names(pomplist) <- paste0("u", 1:10)

panelfood <- panelPomp(pomplist, shared = mif.estimate)

# ---- Simulate B datasets ----
dir.create("simulated_data", showWarnings = FALSE, recursive = TRUE)

# ---- Abnormality filter (mirrors SIRJPF2 passes_bounds): reject a candidate if
#      ANY unit's ANY observation exceeds 4x the real-data max.
#      dent.adult real max 404 -> 1616. ----
passes_bounds <- function(sim_data_list) {
  all(vapply(sim_data_list, function(sim_data) {
    all(
      sim_data$dentadult <= 1616,
      na.rm = TRUE
    )
  }, logical(1)))
}

cat("Simulating up to", max_attempts, "candidate datasets from MLE...\n")
candidate_data <- vector("list", max_attempts)
for (candidate_id in seq_len(max_attempts)) {
  sim_data_list <- list()
  for (u in names(panelfood)) {
    unit_model <- unit_objects(panelfood)[[u]]
    sim <- pomp::simulate(unit_model, nsim = 1, format = "data.frame",
                          params = mif.estimate)
    sim_data_list[[u]] <- sim[, c("day", "dentadult")]
  }
  candidate_data[[candidate_id]] <- sim_data_list
}

valid_candidates <- which(vapply(candidate_data, passes_bounds, logical(1)))
cat(" ", length(valid_candidates), "of", max_attempts, "candidates passed the abnormality bounds.\n")
if (length(valid_candidates) < B) {
  stop("Only ", length(valid_candidates), " of ", max_attempts,
       " candidates satisfy the abnormality bounds; raise max_attempts.")
}
selected_candidate_ids <- valid_candidates[seq_len(B)]
for (b in seq_len(B)) {
  saveRDS(candidate_data[[selected_candidate_ids[[b]]]],
          file = sprintf("simulated_data/sim_data_%03d.rds", b))
}

saveRDS(mif.estimate, file = "simulated_data/true_params.rds")

cat("Done. All", B, "simulated datasets saved in simulated_data/\n")
