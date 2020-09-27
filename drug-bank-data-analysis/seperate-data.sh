#! /bin/bash

input_file=approved-drug-maccs.tsv
id_file=drug-id.txt
cut -f1 ${input_file} > ${id_file}

shuf ${id_file} > shuffled_id.txt

awk 'BEGIN{query_file="query_id.txt"; candidate_id="candidate_id.txt";} (FNR <= 10){ print $1 > query_file; next;} (FNR > 10){ print $1 > candidate_id; next;}' shuffled_id.txt   

rm -f shuffled_id.txt


awk 'BEGIN{FS="\t"; OFS="\t";} NR==FNR{a[$1] = 1; next;} ($1 in a){print $0;}' query_id.txt approved-drug-maccs.tsv > approved-drug-maccs-query.tsv
awk 'BEGIN{FS="\t"; OFS="\t";} NR==FNR{a[$1] = 1; next;} ($1 in a){print $0;}' candidate_id.txt approved-drug-maccs.tsv > approved-drug-maccs-candidate.tsv

awk 'BEGIN{FS="\t"; OFS="\t";} NR==FNR{a[$1] = 1; next;} ($1 in a){print $0;}' query_id.txt approved-drug-pubchem.tsv > approved-drug-pubchem-query.tsv
awk 'BEGIN{FS="\t"; OFS="\t";} NR==FNR{a[$1] = 1; next;} ($1 in a){print $0;}' candidate_id.txt approved-drug-pubchem.tsv > approved-drug-pubchem-candidate.tsv

