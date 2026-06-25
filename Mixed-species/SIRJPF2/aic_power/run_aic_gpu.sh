#!/bin/bash
# run_aic_gpu.sh -- drive the GPU/pypomp AIC selection-rate fits on ONE GPU.
# ---------------------------------------------------------------------------
# Mirrors ../coverage_study/run_sweep.sh: a single unattended loop over the
# study's (truth, b, model) grid, one fresh `python pypomp_aic.py` per fit, with
# FILE-SKIP RESUME (a fit whose results/fit_<model>_<truth>_<b>.rds exists is
# skipped). Pin the GPU + batch via the environment.
#
# Usage:
#   ./run_aic_gpu.sh <tag> <b_start> <b_end> [truth ...] [-- model ...]
#     <tag>            : log label (e.g. A, B)
#     <b_start..b_end> : panel range within 1..50
#     [truth ...]      : shared unit_specific   (default: both)
#   models default to both (null alt); to restrict, append e.g.  -- null
#
# Two GPUs (split the 50 panels), in tmux:
#   tmux new -d -s aicA 'CUDA_VISIBLE_DEVICES=5 PYPOMP_BATCH=300 ./run_aic_gpu.sh A 1  25'
#   tmux new -d -s aicB 'CUDA_VISIBLE_DEVICES=6 PYPOMP_BATCH=300 ./run_aic_gpu.sh B 26 50'
#   tail -f aic_A.log
#
# Full sweep on one GPU:
#   CUDA_VISIBLE_DEVICES=5 PYPOMP_BATCH=300 ./run_aic_gpu.sh A 1 50
set -u

TAG="${1:?need a <tag> (e.g. A)}"
B0="${2:?need <b_start>}"
B1="${3:?need <b_end>}"
shift 3

TRUTHS=(); MODELS=(null alt)
while [ "$#" -gt 0 ]; do
  if [ "$1" = "--" ]; then shift; MODELS=("$@"); break; fi
  TRUTHS+=("$1"); shift
done
[ "${#TRUTHS[@]}" -eq 0 ] && TRUTHS=(shared unit_specific)

# Activate the GPU pypomp env (py313) if a fresh tmux shell lacks it.
if [ "${CONDA_DEFAULT_ENV:-}" != "py313" ]; then
  source /apps/anaconda3/etc/profile.d/conda.sh 2>/dev/null && conda activate py313 2>/dev/null || true
fi

export PYPOMP_USE_GPU="${PYPOMP_USE_GPU:-1}"
export XLA_PYTHON_CLIENT_PREALLOCATE="${XLA_PYTHON_CLIENT_PREALLOCATE:-false}"
# Shared XLA compile cache across all fits (identical array shapes -> one compile).
export JAX_COMPILATION_CACHE_DIR="${JAX_COMPILATION_CACHE_DIR:-$HOME/.jax_cache}"
export PYPOMP_BATCH="${PYPOMP_BATCH:-300}"
PY="${PYPOMP_PY:-python}"
LOG="aic_${TAG}.log"

echo "=== aic-gpu $TAG | GPU=${CUDA_VISIBLE_DEVICES:-?} | panels ${B0}..${B1} | truths=${TRUTHS[*]} | models=${MODELS[*]} | batch=${PYPOMP_BATCH} | start $(date) ===" | tee -a "$LOG"

for truth in "${TRUTHS[@]}"; do
  for model in "${MODELS[@]}"; do
    for b in $(seq "$B0" "$B1"); do
      out="results/fit_${model}_${truth}_${b}.rds"
      if [ -f "$out" ]; then
        echo "skip (exists) $out" | tee -a "$LOG"
        continue
      fi
      echo ">>> $(date) truth=$truth b=$b model=$model -> $out" | tee -a "$LOG"
      if "$PY" pypomp_aic.py "$b" "$model" --truth "$truth" >> "$LOG" 2>&1; then
        echo "    done $out" | tee -a "$LOG"
      else
        echo "    !!! FAILED truth=$truth b=$b model=$model (see $LOG) !!!" | tee -a "$LOG"
      fi
    done
  done
done

echo "=== aic-gpu $TAG done $(date) ===" | tee -a "$LOG"
