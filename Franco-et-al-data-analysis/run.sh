#! /bin/bash


./reduced-fp-list.sh
Rscript --vanilla ./gen-sim-results.R
./combine-sim-results.sh
Rscript --vanilla ./analyze-sim-results.R

