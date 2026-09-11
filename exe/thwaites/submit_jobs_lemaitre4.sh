#!/bin/bash

# Test without the submitting file settings because this is run directly.
module --force purge
module load tis/2018.01
module load MCR/R2018b


# EXPERIMENT.
#exp=deter
#exp=sigma_oce400
#exp=revert_t_m25
# dutrieux2009, dutrieux2012, mathiot_cold, naughten_cold, zhou.
exp=zhou

# PATHS.
# Stoch.
#path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/tau_To_10/$exp
#path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/stoch/tau_To_10/$exp

#path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/extra_runs/tau_To_70/$exp
#path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/stoch/extra_runs/tau_To_70/$exp

# Deter.
#path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/$exp
#path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/$exp

# Reversibility tests.
#path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/deter/revert_basal_friction/weertman/m_10/$exp
#path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/thwaites/deter/revert_basal_friction/weertman/m_10/$exp

# Calibration.
meltname=quad_semi_local_local_slope
p="95"

# No dT corrections.
#path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/ensembles/calibration/dT/$meltname/$exp
#path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/calibration/dT/$meltname/$exp

# dT corrections.
path=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/ensembles/calibration/dT/$meltname/percentile_$p/$exp
path_out=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/output/calibration/dT/$meltname/percentile_$p/$exp


# Slurm job script.
SBATCH_SCRIPT="/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/slurm_job.sh"

# Submitting options: individual_file, ensemble.
option="ensemble"


cd $path


##########################################################################
# OPTION 1.
# Submit just one job.
if [ "$option" = "individual_file" ]; then

file=KoriCalibration_CECI


#srun ./KoriCalibration_CECI 2>&1 > /dev/null
sbatch --export=EXECUTABLE="$path/$file" --chdir="$path" "$SBATCH_SCRIPT"
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