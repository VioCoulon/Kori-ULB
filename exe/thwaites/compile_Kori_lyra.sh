#!/bin/bash
# Bash file to compile Kori-ULB and copy it to CECI cluster.
# Execute it with: bash compile_Kori.sh
path_kori=/home/daniel/models/Kori-ULB
path_kori_subroutines=$path_kori/subroutines


# Define the remote server.
REMOTE_HOST=lyra

# Experiment name.
#exp=sigma_oce050
#exp=HR
exp=sigma_oce010

# LOCAL PATHS.
path_exe=$path_kori/exe/thwaites
#path_param=$path_exe/nic5/$exp          # Stochastic ensemble.

#path_param=$path_exe/nic5/stoch/$exp      # Stochastic.
path_param=$path_exe/$REMOTE_HOST/stoch/tau_To_01/$exp      # Stochastic.


# CLUSTER PATHS.
# Lemaitre.
#path_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/sigma_oce1

# Lyra.
#cluster=lyra
#path_cluster=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/deter/$exp
#path_cluster=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/$exp
path_cluster=/globalsc/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/tau_To_01/


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