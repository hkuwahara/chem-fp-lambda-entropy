#! /bin/bash

mkdir -p maccs
mkdir -p pubchem

Rscript --vanilla ./gen-sim-results2.R
./combine-sim-results2.sh
Rscript --vanilla ./analyze-sim-results3.R
#./perform-random-test2.sh
Rscript --vanilla ./gen-pval-random.R

