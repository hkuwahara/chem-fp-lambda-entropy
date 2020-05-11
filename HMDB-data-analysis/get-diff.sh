#! /bin/bash

toplevel=.

levels=("0.1" "0.2" "0.3")
for level in ${levels[@]};  
do
	paste sim-ce-reduced-0-blood-maccs-melted.tsv sim-ce-reduced-${level}-blood-maccs-melted.tsv > sim-merge-${level}.tsv 
	awk 'BEGIN{FS="\t"; OFS="\t";} NR==FNR{a[$1] = 1; next;}  ($1 in a) && ($2 in a){t=$3-$6; t2=sqrt(t*t)/($3+0.01); print $1, $2, $3, $6, t, t2;}' filtered-blood-compounds.list  sim-merge-${level}.tsv  | awk 'BEGIN{FS="\t";} ($1 != $2){print $0}' | sort -t $'\t' -k6,6rn > blood-maccs-diff-sim-${level}.tsv 
	rm -f sim-merge-${level}.tsv 
done

for level in ${levels[@]};  
do
	paste sim-ce-reduced-0-blood-pubchem-melted.tsv sim-ce-reduced-${level}-blood-pubchem-melted.tsv > sim-merge-${level}.tsv 
	awk 'BEGIN{FS="\t"; OFS="\t";} NR==FNR{a[$1] = 1; next;}  ($1 in a) && ($2 in a){t=$3-$6; t2=sqrt(t*t)/($3+0.01); print $1, $2, $3, $6, t, t2;}' filtered-blood-compounds.list  sim-merge-${level}.tsv  | awk 'BEGIN{FS="\t";} ($1 != $2){print $0}' | sort -t $'\t' -k6,6rn > blood-pubchem-diff-sim-${level}.tsv 
	rm -f sim-merge-${level}.tsv 
done

