#!/user/by2418/.conda/envs/r-pomp/bin/Rscript
Sys.setenv(PATH = paste("/user/by2418/.conda/envs/r-pomp/bin", Sys.getenv("PATH"), sep = ":"))

# ============================================================================
# AIC selection-rate study driver -- ONE job of the 10-job x 100-core scheme.
# ----------------------------------------------------------------------------
# The full study is 200 fits:
#   50 panels  x  2 truths {shared, unit_specific}  x  2 models {null, alt}.
# This driver partitions those 200 fits into 10 jobs of ~20 fits each. Each job
# runs its fits SEQUENTIALLY by shelling out to the existing, validated
# single-fit scripts fit_null.R / fit_alt.R, which themselves each call
# registerDoParallel(cores = 100) and run the 2-round foreach mif2. So each job
# wants a 100-core node and uses the whole node one fit at a time.
#
#   grid_run ... --grid_ncpus=100 ./run_aic_job.R <job>     # <job> in 1..10
#
# RESUME-SAFE: a (truth, b, model) whose results/fit_<model>_<truth>_<b>.rds
# already exists is skipped. Re-submitting a job continues where it left off,
# and re-submitting the whole sweep only fills the gaps.
#
# Job -> fit mapping (deterministic, balanced):
#   The 200 fits are enumerated in a fixed order (truth, model, b) and dealt
#   round-robin into 10 jobs, so each job gets exactly 20 fits and the (slower)
#   alt fits are spread evenly across jobs rather than piling onto two of them.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: run_aic_job.R <job index 1..10>")
job <- as.integer(args[1])
N_JOBS <- 10
if (is.na(job) || job < 1 || job > N_JOBS) {
  stop(sprintf("job index must be 1..%d, got: %s", N_JOBS, args[1]))
}

B <- 50
truths <- c("shared", "unit_specific")
models <- c("null", "alt")

# ---- Enumerate all 200 fits in a fixed order, then deal round-robin ----
grid <- expand.grid(b = 1:B, model = models, truth = truths,
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$truth, grid$model, grid$b), ]
grid$fit_id <- seq_len(nrow(grid))            # 1..200
grid$job <- ((grid$fit_id - 1) %% N_JOBS) + 1 # round-robin -> 20 per job

my_fits <- grid[grid$job == job, ]
cat(sprintf("=== AIC job %d / %d : %d fits assigned ===\n",
            job, N_JOBS, nrow(my_fits)))

# ---- Locate the single-fit scripts and the Rscript that should run them ----
this_dir   <- tryCatch(dirname(normalizePath(sub("--file=", "",
                grep("--file=", commandArgs(FALSE), value = TRUE)[1]))),
                error = function(e) getwd())
if (is.na(this_dir) || !nzchar(this_dir)) this_dir <- getwd()
fit_script <- function(model) file.path(this_dir, sprintf("fit_%s.R", model))
rscript    <- file.path(R.home("bin"), "Rscript")
dir.create(file.path(this_dir, "results"), showWarnings = FALSE)

result_path <- function(model, truth, b)
  file.path(this_dir, "results", sprintf("fit_%s_%s_%d.rds", model, truth, b))

# ---- Run each assigned fit sequentially, skipping finished ones ----
n_done <- 0L; n_skip <- 0L; n_fail <- 0L
for (i in seq_len(nrow(my_fits))) {
  b     <- my_fits$b[i]
  model <- my_fits$model[i]
  truth <- my_fits$truth[i]
  out   <- result_path(model, truth, b)

  if (file.exists(out)) {
    cat(sprintf("[job %d] skip (exists): %s\n", job, basename(out)))
    n_skip <- n_skip + 1L
    next
  }

  cat(sprintf("[job %d] %d/%d  fit_%s.R  b=%d truth=%s  -> %s\n",
              job, i, nrow(my_fits), model, b, truth, basename(out)))
  t0 <- Sys.time()
  status <- system2(rscript, args = c(fit_script(model), b, truth),
                    stdout = "", stderr = "")
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

  if (status == 0 && file.exists(out)) {
    cat(sprintf("[job %d] done  b=%d truth=%s model=%s  (%.1f min)\n",
                job, b, truth, model, dt))
    n_done <- n_done + 1L
  } else {
    cat(sprintf("[job %d] !!! FAILED b=%d truth=%s model=%s (exit %s, %.1f min) !!!\n",
                job, b, truth, model, status, dt))
    n_fail <- n_fail + 1L
  }
}

cat(sprintf("=== AIC job %d complete: %d run, %d skipped, %d failed ===\n",
            job, n_done, n_skip, n_fail))
if (n_fail > 0) quit(status = 1)
