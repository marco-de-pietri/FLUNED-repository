#!/bin/bash

cases_folder="./cases"

# Array of directories
directories=(
    "/01_ACTIVATION/06_mockup2-Incompressible-fixGeom3_kEpsilon/FLUNED_01_DEFAULT/"
    "/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_01_DEFAULT_N16/"
    "/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_02_DEFAULT_N17/"
    "/03_FLUENT/01/FLUNED_01_DEFAULT/"
    "/04_LAMINAR_FLOW/01_LAMINAR_005/FLUNED_01_DEFAULT_N16/"
)


# flunedPost options
options=(
	"-s"
	""
	""
	""
	""

)

length=${#directories[@]}

# Loop through each directory and execute flunedPost
for (( i=0; i<$length; i++ )); do

      (cd $cases_folder${directories[$i]} && exec flunedPost ${options[$i]})
done


