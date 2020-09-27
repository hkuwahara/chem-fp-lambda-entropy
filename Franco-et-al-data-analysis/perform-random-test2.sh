#! /bin/bash



parallel --link Rscript --vanilla randomized-sim-test2.R {1} {2} 10000 ::: maccs pubchem ::: 132 411

./convert-rows-and-cols.sh random-results-maccs-keys-132.tsv random-results-maccs-keys-132-t.tsv
./convert-rows-and-cols.sh random-results-pubchem-keys-411.tsv random-results-pubchem-keys-411-t.tsv


