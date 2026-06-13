#!/bin/bash
# run_sweep.sh -- drive the pypomp coverage-profile sweep on ONE GPU, unattended.
#
# Usage:
#   ./run_sweep.sh <tag> <b_start> <b_end> [param ...]
#     <tag>        : a label for the log file (e.g. A, B)
#     <b_start..b_end> : dataset range (1..100)
#     [param ...]  : profiled params (default: rn ri sigF sigP)
#
# Pin the GPU + memory via the environment, e.g.:
#   CUDA_VISIBLE_DEVICES=6 PYPOMP_BATCH=400 ./run_sweep.sh A 1 50 rn ri sigF sigP
#
# Two GPUs (clean ~2x on the 400-job sweep) -- one process per GPU, in tmux:
#   tmux new -d -s sweepA 'CUDA_VISIBLE_DEVICES=6 PYPOMP_BATCH=400 ./run_sweep.sh A 1 50  rn ri sigF sigP'
#   tmux new -d -s sweepB 'CUDA_VISIBLE_DEVICES=7 PYPOMP_BATCH=400 ./run_sweep.sh B 51 100 rn ri sigF sigP'
#   tail -f sweep_A.log        # watch;  tmux attach -t sweepA  (Ctrl-b d to detach)
#
# RESUMABLE: a (param,b) whose coverage_results/profile_<param>_<bbb>.rds already
# exists is SKIPPED, so you can re-run after an interruption and it continues.
set -u

TAG="${1:?need a <tag> (e.g. A)}"
B0="${2:?need <b_start>}"
B1="${3:?need <b_end>}"
shift 3
PARAMS=("$@"); [ "${#PARAMS[@]}" -eq 0 ] && PARAMS=(rn ri sigF sigP)

# A fresh tmux shell may not have the conda env active -> activate py313.
# (edit the path if your anaconda lives elsewhere)
if [ "${CONDA_DEFAULT_ENV:-}" != "py313" ]; then
  source /apps/anaconda3/etc/profile.d/conda.sh 2>/dev/null && conda activate py313 2>/dev/null || true
fi

export XLA_PYTHON_CLIENT_PREALLOCATE="${XLA_PYTHON_CLIENT_PREALLOCATE:-false}"
export PYPOMP_BATCH="${PYPOMP_BATCH:-400}"
PY="${PYPOMP_PY:-python}"
LOG="sweep_${TAG}.log"

echo "=== sweep $TAG | GPU=${CUDA_VISIBLE_DEVICES:-?} | datasets ${B0}..${B1} | params=${PARAMS[*]} | batch=${PYPOMP_BATCH} | start $(date) ===" | tee -a "$LOG"

for p in "${PARAMS[@]}"; do
  for b in $(seq "$B0" "$B1"); do
    out=$(printf "coverage_results/profile_%s_%03d.rds" "$p" "$b")
    if [ -f "$out" ]; then
      echo "skip (exists) $out" | tee -a "$LOG"
      continue
    fi
    echo ">>> $(date) param=$p b=$b -> $out" | tee -a "$LOG"
    if "$PY" pypomp_coverage_profile.py "$b" "$p" >> "$LOG" 2>&1; then
      echo "    done $out" | tee -a "$LOG"
    else
      echo "    !!! FAILED param=$p b=$b (see $LOG) !!!" | tee -a "$LOG"
    fi
  done
done

echo "=== sweep $TAG done $(date) ===" | tee -a "$LOG"
