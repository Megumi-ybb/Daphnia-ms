#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))
library(pomp)
library(panelPomp)

set.seed(801)

B <- 100
max_attempts <- 500   # draw this many candidates; keep the first B that pass the abnormality bounds

mif.estimate <- c(
  ri       = 708.6439,
  f_Si     = 0.0002316239,
  theta_Si = 0.5255919,
  theta_Ji = 0.07714973,
  sigSi    = 0,
  sigJi    = 0.3588666,
  sigF     = 0.0458815,
  k_Si     = 2.288263
)

dyn_rpro <- Csnippet("
  double Si_term, F_term, Ji_term;
  double noiSi, noiF, noiJi;
  double delta = 0.013;

  noiSi = rnorm(0, sigSi * sqrt(dt));
  noiJi = rnorm(0, sigJi * sqrt(dt));
  noiF  = rnorm(0, sigF  * sqrt(dt));

  Si_term = 0.1 * Ji * dt - theta_Si * Si * dt - delta * Si * dt + Si * noiSi;
  Ji_term = ri * f_Si * F * Si * dt - 0.1 * Ji * dt - theta_Ji * Ji * dt - delta * Ji * dt + Ji * noiJi;
  F_term  = F * noiF - f_Si * F * (Si + 1 * Ji) * dt - delta * F * dt + 0.37 * dt;

  F  += F_term;
  Si += Si_term;
  Ji += Ji_term;

  if (Si < 0.0 || Si > 1e5) { Si = 0.0; error_count += 1; }
  if (F  < 0.0 || F  > 1e20){ F  = 0.0; error_count += 1000; }
  if (Ji < 0.0 || Ji > 1e5) { Ji = 0.0; error_count += 0.001; }

  T_Si = fabs(Si);
")

dyn_init <- Csnippet("
  Si = 3;
  F  = 16.667;
  Ji = 0;
  T_Si = 0.0;
  error_count = 0.0;
")

dmeas <- Csnippet("
  if (error_count > 0.0) {
    lik = -150;
  } else {
    if (give_log) {
      lik = dnbinom_mu(lumadult, k_Si, T_Si, give_log);
    } else {
      lik = dnbinom_mu(lumadult, k_Si, T_Si, give_log);
    }
  }
")

rmeas <- Csnippet("
  lumadult = rnbinom_mu(k_Si, T_Si);
")

pt <- parameter_trans(
  log = c("sigSi","sigF","f_Si","ri","k_Si","sigJi","theta_Si","theta_Ji")
)

obs_times <- (0:9) * 5 + 7

template_data <- data.frame(
  day      = obs_times,
  lumadult = rep(0, length(obs_times))
)

param_names <- c("sigSi","sigF","f_Si","ri","k_Si","sigJi","theta_Si","theta_Ji")
state_names <- c("Si","Ji","error_count","F","T_Si")

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
    obsnames  = c("lumadult"),
    accumvars = c("error_count"),
    paramnames = param_names,
    statenames = state_names
  )
  coef(pomplist[[i]]) <- parameters
}
names(pomplist) <- paste0("u", 1:10)

panelfood <- panelPomp(pomplist, shared = mif.estimate)

dir.create("simulated_data", showWarnings = FALSE, recursive = TRUE)

# ---- Abnormality filter (mirrors SIRJPF2 passes_bounds): reject a candidate if
#      ANY unit's ANY observation exceeds 4x the real-data max.
#      lum.adult real max 104 -> 416. ----
passes_bounds <- function(sim_data_list) {
  all(vapply(sim_data_list, function(sim_data) {
    all(
      sim_data$lumadult <= 416,
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
    sim_data_list[[u]] <- sim[, c("day","lumadult")]
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
cat("True log(ri) =", log(mif.estimate["ri"]), "\n")
