#!/bin/bash

# Test without the submitting file settings because this is run directly.
#module --force purge
#module load tis/2018.01
#module load MCR/R2018b


# EXPERIMENT.
#exp=deter
#exp=sigma_oce400
exp=revert_t_m25     # revert_t_m20_gammas
#exp=DIVA
#exp=CaMIP_Exp3_CECI

# PATHS.
# Stoch.
#path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/tau_To_10/$exp
#path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/stoch/tau_To_10/$exp

# Deter.
#path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/deter/$exp
#path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/deter/$exp

path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/deter/revert_basal_friction/coulomb/m_10/$exp
path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/deter/revert_basal_friction/coulomb/m_10/$exp

# Init.
#path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/init/$exp
#path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/ice_data/eta1e7/ground_melt_0/$exp

# CalvingMIP.
#path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/calvingMIP/Exp3-4/dx_2km/OceanVisc_1e10
#path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/calvingMIP/Exp3-4/dx_2km/OceanVisc_1e10


    
# Slurm job script.
SBATCH_SCRIPT="/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/slurm_job.sh"

cd $path

# COMPILING OPTIONS: individual_file, ensemble.
option="ensemble"  


##########################################################################
# OPTION 1.
# Submit just one job.
if [ "$option" = "individual_file" ]; then

#srun ./RunASE_lemaitre4_stoch_gamma10e5_sigma_oce1_seed2 2>&1 > /dev/null
#srun ./$exp 2>&1 > /dev/null

echo "Working directory of the job : $path"
echo "Submited job                 : $path/$exp"

#mkdir -p "$path_out"

sbatch --export=EXECUTABLE="$path/$exp" --chdir="$path" "$SBATCH_SCRIPT"
##########################################################################



##########################################################################
# OPTION 2.
# Submit several jobs by looping through each file in the directory.
elif [ "$option" = "ensemble" ]; then

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


fi
