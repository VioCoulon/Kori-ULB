#!/bin/bash
# Bash file to compile Kori-ULB and copy it to CECI cluster.
# Execute it with: bash compile_Kori.sh
path_kori=/home/daniel/models/Kori-ULB
path_kori_subroutines=$path_kori/subroutines


# Define the remote server.
REMOTE_HOST=lyra

# Experiment name.
exp=revert_t_m25
  

# LOCAL PATHS.
path_exe=$path_kori/exe/thwaites
#path_param=$path_exe/nic5/$exp          # Stochastic ensemble.

#path_param=$path_exe/$REMOTE_HOST/stoch/tau_To_140/$exp      # Stochastic.
#path_param=$path_exe/$REMOTE_HOST/$exp      # Deter.
path_param=$path_exe/$REMOTE_HOST/deter/revert_basal_friction/weertman/m_10/$exp       # Deter, revert, calving.


# CLUSTER PATHS.
# Lyra.
#path_cluster=/globalsc/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/tau_To_140
#path_exe_cluster=/globalsc/ulb/glaciol/dmoreno/Kori-ULB/exe/precompiled
#path_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/deter/

path_cluster=/globalsc/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/deter/revert_basal_friction/weertman/m_10
path_exe_cluster=/globalsc/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/precompiled


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
#SUBFOLDERS=($(find $path_param -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort))
#printf "%s\n" "${SUBFOLDERS[@]}"


# Compile Kori only once.
echo "Path_par  : $path_param"
echo "Path_exe  : $path_exe"
echo "Compiling : $file"

# RECOMPILE IN CASE OF CHANGES IN THE CODE!
#mcc -m "$file" -a "$path_kori/KoriModel.m" -a "$path_kori_subroutines" -o "$exe_name"
#/usr/local/MATLAB/R2018b/bin/mcc -m "$file" -a "$path_kori/KoriModel.m" -a "$path_kori_subroutines" -o "$exe_name"

ssh $REMOTE_HOST "mkdir -p $path_cluster"

rsync -avz --progress "$exe_name" "$REMOTE_HOST:$path_exe_cluster/"

# Copy folder with param files.
rsync -avz --progress "$path_param" "$REMOTE_HOST:$path_cluster/"



echo "Copying executable to all directories in a single SSH connection"


# Copy via ssh only once and then copy within the cluster.
# File(s) to copy (you can use wildcards like *.sh)
# SSH command to copy files to all subdirectories
# cp "$path_exe_cluster/$exe_name" "\$dir"
ssh $REMOTE_HOST << EOF
for dir in "$path_cluster/$exp"/*/; do
    if [ -d "\$dir" ]; then
        echo "Exp  : \$dir"
        cp "$path_exe_cluster/$exe_name" "\$dir"
    fi
done
EOF


echo "All files copied!"