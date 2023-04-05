#!/bin/bash

while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

for i in 1 2 3 4 5 6
    do
      ./stan_models/${modelname} sample num_warmup=500 num_samples=2500 data file=datafiles/Bcell_Imm.Rdump \
      output file=save_csv/${modelname}_${i}.csv &
    done
