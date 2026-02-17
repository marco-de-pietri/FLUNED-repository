#!/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cases_folder="$script_dir/cases"
manifest_file="$script_dir/cases_manifest.txt"

while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "${line//[[:space:]]/}" || "$line" =~ ^[[:space:]]*# ]] && continue

  dir="${line%%|*}"
  options=""
  if [[ "$line" == *"|"* ]]; then
    options="${line#*|}"
  fi

  case_path="$cases_folder$dir"

  if [[ -n "$options" ]]; then
    (cd "$case_path" && fluned-post $options)
  else
    (cd "$case_path" && fluned-post)
  fi
done <"$manifest_file"
