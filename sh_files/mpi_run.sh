#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test_sanket
#SBATCH -D /opt/mesh/eigg/sanket/cmdstan/MZ_New_dynamics
#SBATCH --exclude=raasay
#SBATCH --nodes=3
#SBATCH --ntasks=120


echo "models/MAP_ki67_SHM_MZ sample num_samples=500 num_warmup=300 data file=datafiles/MZ_data.Rdump output file=output_csv/SHM_T1.csv";

mpirun -np 4 stan_models/MAP_ki67_SHM_MZ sample num_samples=500 num_warmup=300 data file=datafiles/MZ_data.Rdump output file=save_csv/SHM_T1.csv


