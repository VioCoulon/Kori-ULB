#!/bin/bash
# Bash file to compile Kori-ULB and copy it to CECI cluster.
# Execute it with: bash compile_Kori.sh
path_kori=/home/daniel/models/Kori-ULB
path_kori_subroutines=$path_kori/subroutines


# Define the remote server.
REMOTE_HOST=lemaitre4

# Experiment name.
exp=sigma_oce400
#exp=deter
#exp=tau_To_05    

# LOCAL PATHS.
path_exe=$path_kori/exe/thwaites
#path_param=$path_exe/nic5/$exp          # Stochastic ensemble.

path_param=$path_exe/$REMOTE_HOST/stoch/extra_runs/tau_To_70/$exp      # Stochastic.
#path_param=$path_exe/$REMOTE_HOST/$exp      # Deter.


# CLUSTER PATHS.
# Lemaitre4.
path_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/extra_runs/tau_To_70
path_exe_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/precompiled
#path_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/deter/


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
#rsync -avz --progress "$exe_name" "$REMOTE_HOST:$path_exe_cluster/"

# Copy folder with param files.
rsync -avz --progress "$path_param" "$REMOTE_HOST:$path_cluster/"

# Keep all the ensemble under the same directory.
#for folder in "${SUBFOLDERS[@]}"; do
        
    # Define the source and destination paths
#    LOCAL_FILE="$path_param/$folder"
#    echo "Local  : $LOCAL_FILE"

    # Try copying the exe to each folder to copy everything at once.
    #cp "$exe_name" "$path_param/$folder/"

#done


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

#for dir in "$path_cluster/$exp"/*/; do
#    if [ -d "\$dir" ]; then
#        echo "Exp  : $dir"
#        cp "$path_exe_cluster/$exe_name" "\$dir"
#    fi
#done


#find "$path_cluster/$exp" -type d | while read -r dir; do
#    # Check if directory contains any subdirectories
#    if [ -z "$(find "$dir" -mindepth 1 -type d 2>/dev/null)" ]; then
#        cp "$path_exe_cluster/$exe_name" "$dir"
#    fi
#done


# Original.
#rsync -avz --progress "$path_param" "$REMOTE_HOST:$path_cluster/"

echo "All files copied!"