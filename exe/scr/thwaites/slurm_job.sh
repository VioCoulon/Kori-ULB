#!/bin/bash
# Submission script for RunASE
#SBATCH --time=1-23:59:00
#SBATCH --job-name=Kori
#SBATCH -o kori.out
#SBATCH -e kori.err
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000 #2.5Gb
#SBATCH --partition=batch

module --force purge
module load tis/2018.01
module load MCR/R2018b


# Run the executable.
#echo "Running $EXECUTABLE"
"$EXECUTABLE"     # It works, but it leaves the job running in background after the sim finishes.
#"$EXECUTABLE" > /dev/null 2>&1

# Try if this avoids background running.
#srun ."$EXECUTABLE" # Not working.
#srun "$EXECUTABLE"
