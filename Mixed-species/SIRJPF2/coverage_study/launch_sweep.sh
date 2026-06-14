#!/bin/bash
# launch_sweep.sh <procs_per_gpu> [NSTARTS] [GPU_A] [GPU_B]
# ---------------------------------------------------------------------------
# Launch the coverage sweep as <procs_per_gpu> tmux sessions PER GPU, splitting
# datasets 1..100 into disjoint contiguous ranges (GPU_A gets 1..50, GPU_B 51..100).
# Each session runs all 4 study params over its range. BATCH = 60 * NSTARTS (single
# pass; nprof=60 grid). Resumable: finished profile_<param>_<b>.rds are skipped.
#
#   ./launch_sweep.sh 5 30           # 5/GPU, NSTARTS=30, all 100 datasets
#   ./launch_sweep.sh 5 30 30        # ...but only datasets 1..30 (faster, fewer coverage reps)
#   ./launch_sweep.sh 5 30 25 5 6    # datasets 1..25 on GPUs 5 and 6
#
# Stop a run first:  for s in $(tmux ls -F '#S' | grep -E '^[AB][0-9]+$'); do tmux kill-session -t $s; done
set -u
PPG="${1:?usage: launch_sweep.sh <procs_per_gpu> [NSTARTS=50] [N_DATASETS=100] [GPU_A=5] [GPU_B=6]}"
NSTARTS="${2:-50}"
NDS="${3:-100}"                      # number of datasets (1..NDS), split across the two GPUs
GA="${4:-5}"; GB="${5:-6}"
BATCH=$(( 60 * NSTARTS ))
PARAMS="rn ri sigF sigP"
HALF=$(( (NDS + 1) / 2 ))            # GPU_A gets 1..HALF, GPU_B gets HALF+1..NDS

launch_gpu () {                      # <gpu> <lo> <hi> <tag_prefix> <n_chunks>
  local gpu="$1" lo="$2" hi="$3" tagp="$4" n="$5"
  local total=$(( hi - lo + 1 ))
  local per=$(( (total + n - 1) / n ))     # ceil(total / n)
  local s="$lo" i=1
  while [ "$s" -le "$hi" ]; do
    local e=$(( s + per - 1 )); [ "$e" -gt "$hi" ] && e="$hi"
    local tag="${tagp}${i}"
    tmux new -d -s "$tag" \
      "CUDA_VISIBLE_DEVICES=$gpu PYPOMP_NSTARTS=$NSTARTS PYPOMP_BATCH=$BATCH ./run_sweep.sh $tag $s $e $PARAMS"
    echo "  $tag : GPU $gpu  datasets ${s}..${e}"
    s=$(( e + 1 )); i=$(( i + 1 ))
  done
}

echo "=== launch | ${PPG} procs/GPU | NSTARTS=${NSTARTS} BATCH=${BATCH} | datasets 1..${NDS} | GPU ${GA}=(1..${HALF}) ${GB}=(${HALF}+1..${NDS}) ==="
launch_gpu "$GA" 1            "$HALF"  A "$PPG"
launch_gpu "$GB" $((HALF+1))  "$NDS"   B "$PPG"
echo "=== launched $(( 2 * PPG )) sessions. watch: tail -f sweep_A1.log | view: tmux ls ==="
tmux ls 2>/dev/null | grep -E '^[AB][0-9]+:' || true
