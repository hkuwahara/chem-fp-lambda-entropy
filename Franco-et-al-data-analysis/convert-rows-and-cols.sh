#! /bin/bash

input_file=$1
output_file=$2
datamash transpose < ${input_file} > ${output_file}

