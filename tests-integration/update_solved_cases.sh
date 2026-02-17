#!/bin/bash
#
# WARNING - USE ONLY WHEN THERE ARE BREAKING CHANGES IN THE SUMMARY.csv FILE
# FORMAT
#
#
# Replace each RESULTS/SUMMARY.csv in the comparison tree with the one
# from the corresponding case directory.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cases_folder="$script_dir/cases"
cmp_folder="$script_dir/cases_solved"          # destination tree to update
manifest_file="$script_dir/cases_manifest.txt"
res_file="/RESULTS/SUMMARY.csv"
log_file="$script_dir/replace_log"             # optional: keeps a record

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "${line//[[:space:]]/}" || "$line" =~ ^[[:space:]]*# ]] && continue
    dir="${line%%|*}"

    src="$cases_folder$dir$res_file"
    dst="$cmp_folder$dir$res_file"

    if [[ ! -f "$src" ]]; then
        echo "MISSING SOURCE: $src" | tee -a "$log_file"
        continue
    fi

    # Ensure destination directory exists
    mkdir -p "$(dirname "$dst")"

    # Copy (overwrite) the SUMMARY.csv
    cp -f "$src" "$dst"
    echo "UPDATED: $dst" | tee -a "$log_file"
done <"$manifest_file"
