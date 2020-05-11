#! /bin/bash

input_file=/home/hiro/research/kaust/chem-fp/analysis/blood-maccs/diff-sim-0.3.tsv
ref_val="0.20"
output_file=/home/hiro/research/kaust/chem-fp/analysis/blood-maccs/diverged-samples.tsv
awk -v cutoff="${ref_val}" 'function abs(x){return x >= 0 ? x : -x;} BEGIN{FS="\t"; OFS="\t";} {if(abs($5) >= (cutoff+0)){print $1, $2, $3, $4, abs($5);}}' ${input_file} | sort -t $'\t' -k5,5rn | head -50 > ${output_file}


input_file=/home/hiro/research/kaust/chem-fp/analysis/blood-pubchem/diff-sim-0.3.tsv
ref_val="0.20"
output_file=/home/hiro/research/kaust/chem-fp/analysis/blood-pubchem/diverged-samples.tsv
awk -v cutoff="${ref_val}" 'function abs(x){return x >= 0 ? x : -x;} BEGIN{FS="\t"; OFS="\t";} {if(abs($5) >= (cutoff+0)){print $1, $2, $3, $4, abs($5);}}' ${input_file} | sort -t $'\t' -k5,5rn | head -50 > ${output_file}

