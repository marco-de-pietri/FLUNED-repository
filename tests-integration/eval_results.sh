#!/bin/bash

cases_folder="./cases"
cmp_folder="./cases_solved"
res_file="/RESULTS/SUMMARY.csv"
diff_file="./test_results"

# Truncate previous results
: >"$diff_file"

directories=(
  "/01_ACTIVATION/06_mockup2-Incompressible-fixGeom3_kEpsilon/FLUNED_01_DEFAULT/"
  "/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_01_DEFAULT_N16/"
  "/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_02_DEFAULT_N17/"
  "/03_FLUENT/01/FLUNED_01_DEFAULT/"
  "/04_LAMINAR_FLOW/01_LAMINAR_005/FLUNED_01_DEFAULT_N16/"
)

for dir in "${directories[@]}"; do
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
done

echo "Done. Full diff log: $diff_file"
