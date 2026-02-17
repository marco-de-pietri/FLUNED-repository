#!/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cases_folder="$script_dir/cases"
manifest_file="$script_dir/cases_manifest.txt"

while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "${line//[[:space:]]/}" || "$line" =~ ^[[:space:]]*# ]] && continue
  dir="${line%%|*}"
  rm -rf "$cases_folder$dir"
done <"$manifest_file"

rm -f "$script_dir/test_results"
