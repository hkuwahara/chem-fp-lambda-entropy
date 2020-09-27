#! /bin/bash

ref_file=sim_testset.tsv
maccs_file=sim-maccs.tsv
pubchem_file=sim-pubchem.tsv
outfile=sim-results.tsv

awk 'BEGIN{FS="\t"; OFS="\t";} NR==1{header="id" OFS "ref"; next;} NR==FNR{a[$1] = $1 OFS ($4/100); next;} FNR==1{header = header OFS $2 OFS $3; next;} {a[$1] = a[$1] OFS  $2 OFS $3; next;} END{print header; for( i = 1; i <= 100; i++ ) {print a[i "a"];} }' ${ref_file} ${maccs_file} ${pubchem_file} > ${outfile}

