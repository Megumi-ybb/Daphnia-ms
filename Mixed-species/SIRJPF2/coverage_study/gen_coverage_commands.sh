#!/bin/bash
# Generate GRID@CBS command files for the generalized coverage profiler
# (coverage_profile.R <b> <param>). One file per parameter, B datasets each.
#
# Edit PARAMS to match the agreed "all params of interest" scope (#3 / Ed).
# Default below = the parameters carrying a reported CI, minus the 4 already
# run (rn, ri, sigF, sigP). ri/f_Si/probi are profiled so that the composites
# ri*f_Si and probi*f_Si can be formed in collect_coverage.R.
# NOTE: theta_Jn/theta_Ji (juvenile mortality) are profilable but are intentionally
#       NOT in the default PARAMS below -- add them iff they carry a reported CI
#       in the agreed scope (confirm with Ed before the full submission).
set -euo pipefail

PARAMS="xi theta_Sn theta_Si theta_In theta_Ii theta_P f_Sn f_Si probn probi k_Sn k_Si k_In k_Ii sigIn sigIi sigJn sigJi"
B=100

for p in $PARAMS; do
  out="commands_${p}.txt"
  : > "$out"
  for b in $(seq 1 "$B"); do
    printf 'grid_run --grid_submit=batch --grid_ncpus=100 --grid_mem=150G ./coverage_profile.R %d %s  # -> coverage_results/profile_%s_%03d.rds\n' \
      "$b" "$p" "$p" "$b" >> "$out"
  done
  echo "wrote $out (${B} lines)"
done

echo
echo "On the cluster:  chmod +x coverage_profile.R ; for f in commands_*.txt; do bash \"\$f\"; done"
echo "Then locally:    Rscript collect_coverage.R   (after extending its param_names + composite handling)"
