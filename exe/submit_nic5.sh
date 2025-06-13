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

#path=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/
#cd $path
#srun ./RunASE_nic5 2>&1 > /dev/null


# Nic5.
# /scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/sigma_oce400/meltfac3000_seed7
# /scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/sigma_oce200/meltfac0250_seed2
path=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/sigma_oce200/meltfac0250_seed2
path_out=/scratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/stoch/sigma_oce400/meltfac3000_seed7

cd $path

# OPTION 1.
# Submit just one job. Execute with: "sbatch submit_nic5.sh" from terminal.
srun ./RunASE_nic5 2>&1 > /dev/null

# OPTION 2.
# Submit several jobs by looping through each file in the directory.
#for exe in ./*; do
    # Check if the file is executable
#    if [[ -x "$exe" && ! -d "$exe" ]]; then

        # Create a directory with the same name as the executable if it doesn't exist
        # Extract the filename without the path
#        exe_name=$(basename "$exe")
#        echo "Creating directory $path_out/$exe_name"
#        mkdir -p "$path_out/$exe_name"

        # Submit the executable as a job with srun
#        echo "Submitting job for $exe"

#        srun -b ./"$exe" 2>&1 > "$path_out/$exe_name/output.log" 
#    fi
#done

#wait # Ensures script waits for all jobs to finish before exiting if using srun with "&"