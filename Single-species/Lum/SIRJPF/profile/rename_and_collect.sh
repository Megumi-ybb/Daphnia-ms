#!/usr/bin/env bash
#This is a script to help collet the output of profile results from each folder so that we can validate if the profile is successful.
#The collect file scripts works under specific environments, please adjust the code as you needed.
set -euo pipefail
shopt -s nullglob

for out in */*.out; do
  mv "$out" "${out%.out}.txt"
done

summary="last_two_rows_summary.txt"
: > "$summary"

for txt in */*.txt; do
  folder="${txt%%/*}"
  printf '%s :\n' "$folder" >> "$summary"
  tail -n 2 "$txt"          >> "$summary"
  printf '\n'               >> "$summary"
done
echo "Done! See $summary"
