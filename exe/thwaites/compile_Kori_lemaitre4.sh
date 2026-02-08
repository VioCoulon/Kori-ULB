#!/bin/bash
# Bash file to compile Kori-ULB and copy it to CECI cluster.
# Execute it with: bash compile_Kori.sh
path_kori=/home/daniel/models/Kori-ULB
path_kori_subroutines=$path_kori/subroutines


# Define the remote server.
REMOTE_HOST=lemaitre4

# EXPERIMENT NAME.
#exp=revert_t_m05_gammas
exp=INIT_DIVA_2
#exp=sigma_oce400  


# LOCAL PATHS.
path_exe=$path_kori/exe/thwaites
#path_param=$path_exe/nic5/$exp          # Stochastic ensemble.

# Deterministic.
#path_param=$path_exe/$REMOTE_HOST/deter/$exp      # Stochastic.

# Stochastic.
#path_param=$path_exe/$REMOTE_HOST/stoch/smb_0/tau_To_001/$exp      # Stochastic.

# Initialization.
path_param=$path_exe/$REMOTE_HOST/init/$exp      # Deter.


# CLUSTER PATHS.
# Lemaitre4.
# Deterministic.
#path_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/deter

# Stochastic.
#path_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/stoch

# Initialization.
path_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/init

# Path to precompile executable.
path_exe_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/thwaites/precompiled



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
echo "Path_par     : $path_param"
echo "Path_exe     : $path_exe"
echo "Path_cluster : $path_cluster"
echo "Compiling    : $file"

# RECOMPILE IN CASE OF CHANGES IN THE CODE!
mcc -m "$file" -a "$path_kori/KoriModel.m" -a "$path_kori_subroutines" -o "$exe_name"
rsync -avz --progress "$exe_name" "$REMOTE_HOST:$path_exe_cluster/"


# Copy folder with param files.
ssh $REMOTE_HOST "mkdir -p $path_cluster"
rsync -avz --progress "$path_param" "$REMOTE_HOST:$path_cluster/"



echo "Copying executable to all directories in a single SSH connection"
#ssh $REMOTE_HOST "mkdir -p $path_cluster"

# Copy via ssh only once and then copy within the cluster.
# File(s) to copy (you can use wildcards like *.sh)
# SSH command to copy files to all subdirectories
ssh $REMOTE_HOST << EOF
for dir in "$path_cluster/$exp"/*/; do
    if [ -d "\$dir" ]; then
        echo "Exp  : \$dir"
        cp "$path_exe_cluster/$exe_name" "\$dir"
    fi
done
EOF



# Original.
#rsync -avz --progress "$path_param" "$REMOTE_HOST:$path_cluster/"

echo "All files copied!"