#!/bin/bash

casesfolder=./cases

folder1=/01_ACTIVATION/06_mockup2-Incompressible-fixGeom3_kEpsilon/FLUNED_01_DEFAULT/
(rm -rf $casesfolder$folder1)

folder2=/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_01_DEFAULT_N16/
(rm -rf $casesfolder$folder2)

folder3=/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_02_DEFAULT_N17/
(rm -rf $casesfolder$folder3)

folder4=/03_FLUENT/01/FLUNED_01_DEFAULT/
(rm -rf $casesfolder$folder4)

folder5=/04_LAMINAR_FLOW/01_LAMINAR_005/FLUNED_01_DEFAULT_N16/
(rm -rf $casesfolder$folder5)

(rm ./test_results)