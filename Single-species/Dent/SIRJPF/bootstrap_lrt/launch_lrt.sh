#!/bin/bash
# launch_lrt.sh -- spawn K parallel tmux workers on ONE GPU, splitting b=1..100 into K
# contiguous ranges.  Run this from inside a family's bootstrap_lrt/ directory.
#
# Usage:
#   ./launch_lrt.sh <gpu> <K> <target> [target ...]
#     <gpu>        : CUDA device id (e.g. 0..4)
#     <K>          : number of parallel worker processes to stack on this GPU
#     <target ...> : "null" + this family's alternative names (put NULL FIRST)
#
# Example -- Dent/SIRJPF on GPU 1, 6 workers:
#   ./launch_lrt.sh 1 6 null xi theta_Sn theta_In theta_P probn rn f_Sn
#
# WHY STACK: one bootstrap-LRT fit is a ~300-wide vmap -- ~10x narrower than the
# coverage sweep's 3000-wide batches, so a single fit under-fills an A100.  Stack
# K workers to saturate it.  Check `nvidia-smi pmon -s u` (or watch nvidia-smi):
# bump K until the per-process SM% on this GPU SUMS to ~90-100%.  Memory is a
# non-issue (~3-4 GB/proc of 80 GB).
set -u
GPU="${1:?need <gpu> (CUDA device id)}"
K="${2:?need <K> (workers on this GPU)}"
shift 2
TARGETS=("$@"); [ "${#TARGETS[@]}" -eq 0 ] && { echo "need targets (null + alts)"; exit 2; }
HERE="$(cd "$(dirname "$0")" && pwd)"
# b sub-range (default full 1..100). Set LRT_B0/LRT_B1 to split ONE family across
# GPUs, e.g. GPU 1 does LRT_B0=1 LRT_B1=50, GPU 2 does LRT_B0=51 LRT_B1=100.
N0="${LRT_B0:-1}"; N1="${LRT_B1:-100}"
span=$(( N1 - N0 + 1 ))
per=$(( (span + K - 1) / K ))        # ceil(span/K) -> K contiguous chunks within [N0,N1]

i=0; b0=$N0
while [ "$b0" -le "$N1" ]; do
  i=$((i+1)); b1=$(( b0 + per - 1 )); [ "$b1" -gt "$N1" ] && b1=$N1
  tag="g${GPU}w${i}"; sess="lrt_${tag}"
  tmux new -d -s "$sess" \
    "cd '$HERE' && CUDA_VISIBLE_DEVICES=$GPU PYPOMP_BATCH=${PYPOMP_BATCH:-300} PYPOMP_NSTARTS=${PYPOMP_NSTARTS:-300} ./run_lrt_sweep.sh $tag $b0 $b1 ${TARGETS[*]}"
  echo "launched tmux '$sess' : GPU $GPU  b ${b0}..${b1}  ($((b1-b0+1)) reps x ${#TARGETS[@]} targets)"
  b0=$(( b1 + 1 ))
done
echo
echo "GPU $GPU: $i workers launched."
echo "  watch:   tmux ls   |   tail -f ${HERE}/lrt_sweep_g${GPU}w1.log   |   nvidia-smi"
echo "  attach:  tmux attach -t lrt_g${GPU}w1   (Ctrl-b d to detach)"
echo "  stop:    for s in \$(tmux ls -F '#S' | grep ^lrt_g${GPU}w); do tmux kill-session -t \$s; done"
