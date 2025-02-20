#!/bin/bash

cases_folder="./cases"
cmp_folder="./cases_solved"
res_file="/RESULTS/SUMMARY.csv"
diff_file="./test_results"


# Array of directories
directories=(
    "/01_ACTIVATION/06_mockup2-Incompressible-fixGeom3_kEpsilon/FLUNED_01_DEFAULT/"
    "/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_01_DEFAULT_N16/"
    "/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_02_DEFAULT_N17/"
    "/03_FLUENT/01/FLUNED_01_DEFAULT/"
    "/04_LAMINAR_FLOW/01_LAMINAR_005/FLUNED_01_DEFAULT_N16/"
)

length=${#directories[@]}

# Loop through each directory and compare files
for (( i=0; i<$length; i++ )); do

      echo "FLUNED TEST FOLDER: ${directories[$i]}" >> $diff_file 
	  
      if cmp -s "$cases_folder${directories[$i]}$res_file" "$cmp_folder${directories[$i]}$res_file"; then
          echo "PASSED" >> $diff_file 
      else
          echo "NOT PASSED - SUMMARY.csv differences below" >> $diff_file 
		  diff --unified=0 -i -w  "$cases_folder${directories[$i]}$res_file" "$cmp_folder${directories[$i]}$res_file" >> $diff_file
      fi
	  
	  echo "" >> $diff_file
	  
done
