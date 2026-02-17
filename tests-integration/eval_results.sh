#!/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cases_folder="$script_dir/cases"
cmp_folder="$script_dir/cases_solved"
manifest_file="$script_dir/cases_manifest.txt"
res_file="/RESULTS/SUMMARY.csv"
diff_file="$script_dir/test_results"

# Truncate previous results
: >"$diff_file"

while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "${line//[[:space:]]/}" || "$line" =~ ^[[:space:]]*# ]] && continue
  dir="${line%%|*}"

  echo "▶ Checking: $dir"
  echo "FLUNED TEST FOLDER: $dir" >>"$diff_file"

  if diff -q -i -w \
    "$cases_folder$dir$res_file" \
    "$cmp_folder$dir$res_file" >/dev/null; then

    echo "  ✔ PASSED"
    echo "PASSED" >>"$diff_file"

  else
    echo "  ✘ NOT PASSED (differences found)"
    echo "NOT PASSED - SUMMARY.csv differences below" >>"$diff_file"

    diff --unified=0 -i -w \
      "$cases_folder$dir$res_file" \
      "$cmp_folder$dir$res_file" >>"$diff_file"
  fi

  echo >>"$diff_file"
  echo
done <"$manifest_file"

echo "Done. Full diff log: $diff_file"
