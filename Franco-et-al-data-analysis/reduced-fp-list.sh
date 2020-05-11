#! /bin/bash

analysis_toplevel=$(cd ..; pwd)
fp_types=("maccs" "pubchem")
reduced_levels=("0.1" "0.2" "0.3")

for fp_type in ${fp_types[@]};
do
	for level in ${reduced_levels[@]};
	do
		input_file=${analysis_toplevel}/HMDB-data-analysis/ce-reduced-${level}-blood-${fp_type}.tsv
		head -1 ${input_file} > blood-${fp_type}-${level}-fp-list.tsv 
	done 	
done

