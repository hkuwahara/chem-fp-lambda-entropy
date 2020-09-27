#! /bin/bash

mkdir -p maccs
mkdir -p pubchem

#./seperate-data.sh

Rscript --vanilla ./reduce-features2.R maccs
Rscript --vanilla ./reduce-features2.R pubchem

./parallel-sim-matrix-gen3.sh
