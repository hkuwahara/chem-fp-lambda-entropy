#! /bin/bash

parallel -a <(echo -e "maccs\npubchem") -a query_id.txt Rscript --vanilla ./sim-matrix-gen3.R {1} {2}  


