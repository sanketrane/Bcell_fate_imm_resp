#!/bin/sh

#SBATCH --account=iicd           # Replace ACCOUNT with your group account name
#SBATCH --job-name=firstjob      # The job name
#SBATCH -c 25                     # The number of cpu cores to use (up to 32 cores per server)

Rscript Parallel_clone.R 1 1 > job.out 2> erjob.err &
echo "DONE!"
date

