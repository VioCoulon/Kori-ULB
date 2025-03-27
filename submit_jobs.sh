#!/bin/bash

# Test without the submitting file settings because this is run directly.
module --force purge
module load tis/2018.01
module load MCR/R2018b

#path=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/
#cd $path
#srun ./RunASE_nic5 2>&1 > /dev/null


# Execute file with ./submit_job.sh

# Desired experiment.
#exp=seed100

# Nic5.
#path=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/sigma_oce100/$exp
#path_out=/scratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/stoch/sigma_oce100/$exp


# Test for all subfolders.
#exp=sigma_oce050
exp=HR

path=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/deter/$exp
path_out=/scratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/deter/$exp



# Slurm job script.
SBATCH_SCRIPT="/scratch/ulb/glaciol/dmoreno/Kori-ULB/slurm_job.sh"

cd $path


##########################################################################
# OPTION 1.
# Submit just one job.
#srun ./RunASE_lemaitre4_stoch_gamma10e5_sigma_oce1_seed2 2>&1 > /dev/null
##########################################################################



##########################################################################
# OPTION 2.
# Submit several jobs by looping through each file in the directory.

find . -type f -executable | while read -r exe; do
    # Extract the filename without the path
    exe_name=$(basename "$exe")
    
    # Create a corresponding output directory
    exe_dir=$(dirname "$exe")

    # Get the absolute path of the executable (including name).
    exe_abs_path=$(realpath "$exe")

    # Extract the directory where the executable is located
    exe_path=$(dirname "$exe_abs_path")
    exe_full_path="$exe_path/$exe_name"
    
    # Remove the executable name to get the path below the executable
    output_dir="$path_out/${exe_dir#./}"

    mkdir -p "$output_dir"

    # Submit the executable as a job with sbatch
    
    #sbatch --export=EXECUTABLE="$exe_dir" "$SBATCH_SCRIPT"

    echo "Working directory of the job : $exe_path"
    echo "Output directory             : $output_dir"
    echo "Submitted job                : $exe_full_path"
    
    sbatch --export=EXECUTABLE="$exe_full_path" --chdir="$exe_path" "$SBATCH_SCRIPT"
    #srun "$exe" "2>&1"
done
##########################################################################
