#!/bin/bash
# Bash file to compile Kori-ULB and copy it to CECI cluster.
# Execute it with: bash compile_Kori.sh
path_kori=/home/daniel/models/Kori-ULB
path_kori_subroutines=$path_kori/subroutines



# Define the remote server.
REMOTE_HOST=nic5

# Experiment name.
exp=deter
#exp=sigma_oce400

# LOCAL PATHS.
path_exe=$path_kori/exe/thwaites
#path_param=$path_exe/nic5/$exp          # Stochastic ensemble.

#path_param=$path_exe/nic5/stoch/tau_To_20/$exp      # Stochastic.
path_param=$path_exe/$REMOTE_HOST/$exp      # Deter.



# CLUSTER PATHS.
# Nic5.
#path_cluster=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/deter/$exp
#path_cluster=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/$exp
path_cluster=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/


# Enter path with matlab scripts to be compiled.
cd $path_exe


# OPTION 1.
# Compilation for single files.
#file_name=RunASE_lemaitre4.m
#exe_name=RunASE_lemaitre4_stoch_gamma10e5_sigma_oce1_seed2
#echo "Compiling: $file_name"
#echo "Exe_name : $exe_name"
#mcc -m $file_name -a $path_kori/KoriModel.m -a $path_kori_subroutines -o $exe_name
#scp $exe_name lemaitre4:$path_lemaitre4


# OPTION 2.
# Loop over each file in the directory.
file=RunASE_ceci.m
exe_name=RunASE_ceci

# Define the remote server.
#REMOTE_HOST="nic5"
REMOTE_HOST=nic5

# Create an array with subfolder names.
SUBFOLDERS=($(find $path_param -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort))

printf "%s\n" "${SUBFOLDERS[@]}"


# Compile Kori only once.
echo "Path_par  : $path_param"
echo "Path_exe  : $path_exe"
echo "Compiling : $file"
mcc -m "$file" -a "$path_kori/KoriModel.m" -a "$path_kori_subroutines" -o "$exe_name"



# Keep all the ensemble under the same directory.
for folder in "${SUBFOLDERS[@]}"; do
        
    # Define the source and destination paths
    LOCAL_FILE="$path_param/$folder"
    echo "Local  : $LOCAL_FILE"

    # Try copying the exe to each folder to copy everythin at once.
    cp "$exe_name" "$path_param/$folder/"

done


echo "Creating remote directories in a single SSH connection"
ssh $REMOTE_HOST "mkdir -p $path_cluster"
rsync -avz --progress "$path_param" "$REMOTE_HOST:$path_cluster/"

echo "All files copied!"