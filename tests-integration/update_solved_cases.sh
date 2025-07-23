#!/bin/bash
#
# WARNING - USE ONLY WHEN THERE ARE BREAKING CHANGES IN THE SUMMARY.csv FILE
# FORMAT
#
#
# Replace each RESULTS/SUMMARY.csv in the comparison tree with the one
# from the corresponding case directory.

cases_folder="./cases"
cmp_folder="./cases_solved"          # destination tree to update
res_file="/RESULTS/SUMMARY.csv"
log_file="./replace_log"             # optional: keeps a record

# Array of relative test-case directories
directories=(
    "/01_ACTIVATION/06_mockup2-Incompressible-fixGeom3_kEpsilon/FLUNED_01_DEFAULT/"
    "/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_01_DEFAULT_N16/"
    "/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_02_DEFAULT_N17/"
    "/03_FLUENT/01/FLUNED_01_DEFAULT/"
    "/04_LAMINAR_FLOW/01_LAMINAR_005/FLUNED_01_DEFAULT_N16/"
)


for dir in "${directories[@]}"; do
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
done
