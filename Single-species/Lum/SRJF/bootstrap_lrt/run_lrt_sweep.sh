#!/bin/bash
# run_lrt_sweep.sh -- drive the pypomp/GPU bootstrap-LRT sweep for ONE family on ONE
# GPU, unattended.  One python process at a time; stack several of these (disjoint b
# ranges) on the same GPU via launch_lrt.sh to fill it.
#
# Usage:
#   ./run_lrt_sweep.sh <tag> <b_start> <b_end> <target> [target ...]
#     <tag>            : label for the log + per-worker done-list (e.g. g1w1)
#     <b_start..b_end> : bootstrap-replicate range (1..100)
#     <target ...>     : "null" plus this family's alternative names; put NULL FIRST
#                        so each b gets a complete (null + all alts) LRT set in order.
#
# Pin the GPU + budget via the environment, e.g.:
#   CUDA_VISIBLE_DEVICES=1 PYPOMP_NSTARTS=300 PYPOMP_BATCH=300 \
#     ./run_lrt_sweep.sh g1w1 1 17 null theta_In ...
#
# WRITES THE CANONICAL lrt_*.rds (PYPOMP_OVERWRITE=1) that collect_lrt.R reads --
# this INTENTIONALLY overwrites the stale R results for this family (that is the
# whole point of the re-run).  The originals are recoverable via `git checkout`.
#
# RESUME: every completed (target,b) is appended to lrt_done_<tag>.log and SKIPPED on
# re-run.  We resume on the done-list, NOT on output existence, precisely because the
# stale lrt_*.rds already exist.  Delete lrt_done_<tag>.log to force a full recompute.
set -u

cd "$(dirname "$0")" || exit 1          # always run from the family's bootstrap_lrt/ dir

TAG="${1:?need <tag> (e.g. g1w1)}"
B0="${2:?need <b_start>}"
B1="${3:?need <b_end>}"
shift 3
TARGETS=("$@"); [ "${#TARGETS[@]}" -eq 0 ] && { echo "need >=1 target (null + alts)"; exit 2; }

# A fresh tmux shell may not have the conda env active -> activate py313
# (the env with jax[cuda12] + pypomp(-e) + pyreadr).  Edit the path if anaconda
# lives elsewhere on this node.
if [ "${CONDA_DEFAULT_ENV:-}" != "py313" ]; then
  source /apps/anaconda3/etc/profile.d/conda.sh 2>/dev/null && conda activate py313 2>/dev/null || true
fi

export XLA_PYTHON_CLIENT_PREALLOCATE="${XLA_PYTHON_CLIENT_PREALLOCATE:-false}"  # let procs share a GPU
# Persistent XLA compile cache shared across ALL jobs/GPUs.  Each (target,b) runs in a
# fresh python process; without this every fit recompiles the mif/pfilter graph (~min).
# All fits in a family have identical array shapes (same units/obs/batch), so after the
# first compile they all hit the cache.  Safe for concurrent procs (atomic writes).
export JAX_COMPILATION_CACHE_DIR="${JAX_COMPILATION_CACHE_DIR:-$HOME/.jax_cache}"
export PYPOMP_USE_GPU="${PYPOMP_USE_GPU:-1}"
export PYPOMP_X64="${PYPOMP_X64:-1}"            # double precision (states span 0..1e20)
export PYPOMP_RUN_LEVEL="${PYPOMP_RUN_LEVEL:-3}"  # 3 = real R settings (Np=Mp=1500, Nmif 150/250)
export PYPOMP_NSTARTS="${PYPOMP_NSTARTS:-300}"  # restarts per fit
export PYPOMP_BATCH="${PYPOMP_BATCH:-300}"      # one even batch over the 300 starts (no ragged recompile)
export PYPOMP_OVERWRITE=1                        # write canonical lrt_*.rds (feeds collect_lrt.R)
PY="${PYPOMP_PY:-python}"
LOG="lrt_sweep_${TAG}.log"
DONE="lrt_done_${TAG}.log"; touch "$DONE"

echo "=== LRT sweep $TAG | GPU=${CUDA_VISIBLE_DEVICES:-?} | b ${B0}..${B1} | targets=${TARGETS[*]} | nstarts=${PYPOMP_NSTARTS} batch=${PYPOMP_BATCH} run_level=${PYPOMP_RUN_LEVEL} | start $(date) ===" | tee -a "$LOG"

for b in $(seq "$B0" "$B1"); do
  for t in "${TARGETS[@]}"; do
    # Resume on canonical OUTPUT existence too (not just the per-worker done-list),
    # so the run can be re-sharded across GPUs without redoing finished fits.
    # Safe here: this run's results dirs started empty, so any lrt_*.rds present is fresh.
    if [ "$t" = "null" ]; then outf="results_null/lrt_null_${b}.rds"; else outf="results_alt/lrt_${t}_${b}.rds"; fi
    if [ -f "$outf" ]; then echo "skip (output exists) target=$t b=$b" | tee -a "$LOG"; continue; fi
    if grep -qx "$t $b" "$DONE"; then
      echo "skip (done) target=$t b=$b" | tee -a "$LOG"; continue
    fi
    echo ">>> $(date) target=$t b=$b" | tee -a "$LOG"
    if "$PY" pypomp_bootstrap_lrt.py "$b" "$t" >> "$LOG" 2>&1; then
      echo "$t $b" >> "$DONE"
      echo "    done target=$t b=$b" | tee -a "$LOG"
    else
      echo "    !!! FAILED target=$t b=$b (see $LOG) !!!" | tee -a "$LOG"
    fi
  done
done
echo "=== sweep $TAG COMPLETE $(date) ===" | tee -a "$LOG"
