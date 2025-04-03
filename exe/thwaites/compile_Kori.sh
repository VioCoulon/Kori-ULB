#!/bin/bash
# Bash file to compile Kori-ULB and copy it to CECI cluster.
# Execute it with: bash compile_Kori.sh
path_kori=/home/daniel/models/Kori-ULB
path_kori_subroutines=$path_kori/subroutines

# Experiment name.
#exp=sigma_oce050
#exp=HR
exp=sigma_oce100

# LOCAL PATHS.
path_exe=$path_kori/exe/thwaites
#path_param=$path_exe/nic5/$exp          # Stochastic ensemble.

#path_param=$path_exe/nic5/stoch/HR/$exp      # Stochastic.
path_param=$path_exe/nic5/stoch/$exp      # Stochastic.


# CLUSTER PATHS.
# Lemaitre.
#path_cluster=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/sigma_oce1

# Nic5.
cluster=nic5
#path_cluster=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/deter/$exp
path_cluster=/scratch/ulb/glaciol/dmoreno/Kori-ULB/exe/stoch/$exp


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
file=RunASE_nic5.m
exe_name=RunASE_nic5

# Define the remote server.
REMOTE_HOST="nic5"

# Create an array with subfolder names.
SUBFOLDERS=($(find $path_param -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort))

printf "%s\n" "${SUBFOLDERS[@]}"


# Compile Kori only once.
echo "Path_exe  : $path_exe"
echo "Compiling : $file"
mcc -m "$file" -a "$path_kori/KoriModel.m" -a "$path_kori_subroutines" -o "$exe_name"



# Collect all directories that need to be created.
REMOTE_DIRS=()
FILES_TO_COPY=()

#for main in "${MAIN_FOLDERS[@]}"; do
#    for sub in "${SUBFOLDERS[@]}"; do
        
        # Define the source and destination paths
#        LOCAL_FILE="$path_param/$main/${main}_$sub/params.mat"
#        REMOTE_DIR="$path_cluster/$main/${main}_$sub"

#        echo "Local  : $LOCAL_FILE"
#        echo "Remote : $REMOTE_DIR"

        # Check if the file exists before copying
#        if [[ -f "$LOCAL_FILE" ]]; then
#            REMOTE_DIRS+=("$REMOTE_DIR")
#            FILES_TO_COPY+=("$LOCAL_FILE|$REMOTE_DIR")
#        else
#            echo "WARNING: File $LOCAL_FILE does not exist, skipping..."
#        fi
#    done
#done

# Keep all the ensemble under the same directory.
for folder in "${SUBFOLDERS[@]}"; do
        
    # Define the source and destination paths
    LOCAL_FILE="$path_param/$folder/params.mat"
    REMOTE_DIR="$path_cluster/$folder"

    echo "Local  : $LOCAL_FILE"
    echo "Remote : $REMOTE_DIR"

    # Check if the file exists before copying
    if [[ -f "$LOCAL_FILE" ]]; then
        REMOTE_DIRS+=("$REMOTE_DIR")
        FILES_TO_COPY+=("$LOCAL_FILE|$REMOTE_DIR")
    else
        echo "WARNING: File $LOCAL_FILE does not exist, skipping..."
    fi
done


# Create all necessary directories in a single SSH connection
if [[ ${#REMOTE_DIRS[@]} -gt 0 ]]; then
    echo "Creating remote directories in a single SSH connection"
    ssh $REMOTE_HOST "mkdir -p ${REMOTE_DIRS[*]}"
fi

# Copy all files using rsync
for entry in "${FILES_TO_COPY[@]}"; do
    IFS='|' read -r LOCAL_FILE REMOTE_DIR <<< "$entry"
    
    echo "Copying files to $REMOTE_DIR"
    
    rsync -avz --progress "$LOCAL_FILE" "$REMOTE_HOST:$REMOTE_DIR/"
    rsync -avz --progress "$path_exe/$exe_name" "$REMOTE_HOST:$REMOTE_DIR/"
    
    echo "File copied successfully!"
done

echo "All files copied!"