#!/bin/bash

casesfolder=./cases
cmpfolder=./cases_solved
resfolder=/RESULTS/SUMMARY.csv
diff_file=./test_results

folder1=/01_ACTIVATION/06_mockup2-Incompressible-fixGeom3_kEpsilon/FLUNED_01_DEFAULT 
echo $folder1 >> $diff_file
(diff $casesfolder$folder1$resfolder $cmpfolder$folder1$resfolder >> $diff_file)

folder2=/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_01_DEFAULT_N16
echo $folder2 >> $diff_file
(diff $casesfolder$folder2$resfolder $cmpfolder$folder2$resfolder >> $diff_file)

folder3=/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_02_DEFAULT_N17
echo $folder3 >> $diff_file
(diff $casesfolder$folder3$resfolder $cmpfolder$folder3$resfolder >> $diff_file)

folder4=/03_FLUENT/01/FLUNED_01_DEFAULT
echo $folder4 >> $diff_file
(diff $casesfolder$folder4$resfolder $cmpfolder$folder4$resfolder >> $diff_file)

folder5=/04_LAMINAR_FLOW/01_LAMINAR_005/FLUNED_01_DEFAULT_N16
echo $folder5 >> $diff_file
(diff $casesfolder$folder5$resfolder $cmpfolder$folder5$resfolder >> $diff_file)
