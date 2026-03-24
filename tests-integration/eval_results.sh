#!/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cases_folder="$script_dir/cases"
cmp_folder="$script_dir/cases_solved"
manifest_file="$script_dir/cases_manifest.txt"
res_xml_file="/RESULTS/fluned_summary.xml"
diff_file="$script_dir/test_results"
xml_cmp_script="$script_dir/compare_fluned_summary_xml.py"
xml_rel_tol="${XML_REL_TOL:-1e-5}"
xml_abs_tol="${XML_ABS_TOL:-1e-12}"

# Truncate previous results
: >"$diff_file"

while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "${line//[[:space:]]/}" || "$line" =~ ^[[:space:]]*# ]] && continue
  dir="${line%%|*}"

  echo "▶ Checking: $dir"
  echo "FLUNED TEST FOLDER: $dir" >>"$diff_file"

  src_xml="$cases_folder$dir$res_xml_file"
  ref_xml="$cmp_folder$dir$res_xml_file"

  if [[ ! -f "$src_xml" ]]; then
    echo "  ✘ NOT PASSED (missing source XML)"
    echo "NOT PASSED - missing source XML: $src_xml" >>"$diff_file"
  elif [[ ! -f "$ref_xml" ]]; then
    echo "  ✘ NOT PASSED (missing reference XML)"
    echo "NOT PASSED - missing reference XML: $ref_xml" >>"$diff_file"
  elif python3 "$xml_cmp_script" \
    --rel-tol "$xml_rel_tol" \
    --abs-tol "$xml_abs_tol" \
    "$ref_xml" \
    "$src_xml" >>"$diff_file"; then

    echo "  ✔ PASSED"
    echo "PASSED" >>"$diff_file"

  else
    echo "  ✘ NOT PASSED (XML differences found)"
  fi

  echo >>"$diff_file"
  echo
done <"$manifest_file"

echo "Done. Full diff log: $diff_file"
