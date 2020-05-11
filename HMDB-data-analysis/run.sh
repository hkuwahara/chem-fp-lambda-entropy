#! /bin/bash

Rscript --vanilla ./reduce-features2.R maccs
Rscript --vanilla ./reduce-features2.R pubchem
./parallel-sim-matrix-gen2.sh
Rscript --vanilla ./sim_results_gen.R
./get-diff.sh
./find-different-samples.sh

