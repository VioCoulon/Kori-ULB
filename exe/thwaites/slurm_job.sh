#!/bin/bash
# Submission script for RunASE
#SBATCH --time=1-23:59:00
#SBATCH --job-name=Kori
#SBATCH -o kori.out
#SBATCH -e kori.err
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000 # 2.5Gb, 3000
#SBATCH --partition=batch

module --force purge
module load tis/2018.01
module load MCR/R2018b


# Ensure the required library paths are included
# ldd RunASE_ceci | grep mwlaunchermain
#export LD_LIBRARY_PATH=/opt/cecisw/noarch/easybuild/tis-2018.01/software/MCR/R2018b/v95/bin/glnxa64:/opt/cecisw/noarch/easybuild/tis-2018.01/software/MCR/R2018b/v95/runtime/glnxa64:/opt/cecisw/noarch/easybuild/tis-2018.01/software/MCR/R2018b/v95/sys/os/glnxa64:$LD_LIBRARY_PATH

# Optional: confirm inside job log
#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

# Run the executable.
#echo "Running $EXECUTABLE"
"$EXECUTABLE"     # It works, but it leaves the job running in background after the sim finishes.
#srun "$EXECUTABLE"

#"$EXECUTABLE" > /dev/null 2>&1

# Try if this avoids background running.
#srun ."$EXECUTABLE" # Not working.
#srun "$EXECUTABLE"
