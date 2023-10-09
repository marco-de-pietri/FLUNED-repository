#!/bin/bash

casesfolder=./cases

(cd $casesfolder/01_ACTIVATION/06_mockup2-Incompressible-fixGeom3_kEpsilon/FLUNED_01_DEFAULT && exec flunedPost -s)
(cd $casesfolder/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_01_DEFAULT_N16/ && exec flunedPost)
(cd $casesfolder/02_DECAY/01_tewn_2022-Incompressible-tewn_2022_ref33/FLUNED_02_DEFAULT_N17/ && exec flunedPost)
(cd $casesfolder/03_FLUENT/01/FLUNED_01_DEFAULT/ && exec flunedPost)
(cd $casesfolder/04_LAMINAR_FLOW/01_LAMINAR_005/FLUNED_01_DEFAULT_N16/ && exec flunedPost)

