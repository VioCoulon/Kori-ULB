#!/bin/bash
# Submission script for RunASE
#SBATCH --time=1-23:59:00
#SBATCH --job-name=Kori
#SBATCH -o kori.out
#SBATCH -e kori.err
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2500 #2.5Gb
#SBATCH --partition=batch

module --force purge
module load tis/2018.01
module load MCR/R2018b

path=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/

cd $path
srun ./RunASE_cluster_long 2>&1 > /dev/null
