#!/bin/bash
# lrt_progress.sh -- show GPU bootstrap-LRT progress (finished canonical jobs per
# family/target). Counts results_{null,alt}/lrt_*.rds (excludes *_pypomp sidecars).
# Run from anywhere:
#   bash lrt_progress.sh
#   watch -n 60 bash ~/Documents/Daphnia-ms/lrt_progress.sh
ROOT="$(cd "$(dirname "$0")" && pwd)"
FAMS=(
  "SRJF2|Mixed-species/SRJF2|null f_Si_f_Sn rn_ri theta_Sn_theta_Si"
  "Dent/SIRJPF|Single-species/Dent/SIRJPF|null xi theta_Sn theta_In theta_P probn rn f_Sn"
  "Dent/SRJF|Single-species/Dent/SRJF|null theta_Sn rn f_Sn"
  "Lum/SIRJPF|Single-species/Lum/SIRJPF|null xi theta_Si theta_Ii theta_P probi ri f_Si"
  "Lum/SRJF|Single-species/Lum/SRJF|null theta_Si ri f_Si"
)
gtot=0; gdone=0
printf "%-12s %9s   %s\n" "FAMILY" "DONE/TOT" "per-target (done/100)"
printf -- "----------------------------------------------------------------------------\n"
for spec in "${FAMS[@]}"; do
  lbl="${spec%%|*}"; rest="${spec#*|}"; f="${rest%%|*}"; tgts="${rest#*|}"; bl="$ROOT/$f/bootstrap_lrt"
  tot=0; done=0; line=""
  for t in $tgts; do
    if [ "$t" = null ]; then
      n=$(ls "$bl"/results_null/lrt_null_[0-9]*.rds 2>/dev/null | grep -vc _pypomp)
    else
      n=$(ls "$bl"/results_alt/lrt_${t}_[0-9]*.rds 2>/dev/null | grep -vc _pypomp)
    fi
    tot=$((tot + 100)); done=$((done + n)); line="$line ${t}:${n}"
  done
  gtot=$((gtot + tot)); gdone=$((gdone + done))
  printf "%-12s %4d/%-4d  %s\n" "$lbl" "$done" "$tot" "$line"
done
printf -- "----------------------------------------------------------------------------\n"
pct=$(( gtot > 0 ? 100 * gdone / gtot : 0 ))
workers=$(tmux ls 2>/dev/null | grep -c '^lrt_g')
printf "%-12s %4d/%-4d  (%d%%)   live tmux workers: %s\n" "TOTAL" "$gdone" "$gtot" "$pct" "$workers"
