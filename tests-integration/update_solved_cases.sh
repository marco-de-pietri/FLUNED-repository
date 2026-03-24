#!/bin/bash
#
# WARNING - USE ONLY WHEN THE REFERENCE RESULT FILES OR THEIR FORMAT CHANGE
#
#
# Replace each RESULTS/SUMMARY.csv and RESULTS/fluned_summary.xml in the
# comparison tree with the files from the corresponding case directory.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cases_folder="$script_dir/cases"
cmp_folder="$script_dir/cases_solved" # destination tree to update
manifest_file="$script_dir/cases_manifest.txt"
res_file="/RESULTS/SUMMARY.csv"
res_xml_file="/RESULTS/fluned_summary.xml"
log_file="$script_dir/replace_log" # optional: keeps a record

while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "${line//[[:space:]]/}" || "$line" =~ ^[[:space:]]*# ]] && continue
  dir="${line%%|*}"

  src="$cases_folder$dir$res_file"
  dst="$cmp_folder$dir$res_file"

  src_xml="$cases_folder$dir$res_xml_file"
  dst_xml="$cmp_folder$dir$res_xml_file"

  if [[ ! -f "$src" ]]; then
    echo "MISSING SOURCE: $src" | tee -a "$log_file"
    continue
  fi

  # Ensure destination directory exists
  mkdir -p "$(dirname "$dst")"

  # Copy (overwrite) the reference summary files
  cp -f "$src" "$dst"
  cp -f "$src_xml" "$dst_xml"
  echo "UPDATED: $dst" | tee -a "$log_file"
done <"$manifest_file"
