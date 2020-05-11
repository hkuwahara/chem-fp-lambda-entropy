#! /bin/bash

toplevel_dir=$PWD
basename -s ".tsv" ${toplevel_dir}/ce-reduced*.tsv | parallel Rscript --vanilla ./sim-matrix-gen2.R ${toplevel_dir}/${exp_id} {}  


