#!/bin/bash
# Bash file to compile Kori-ULB and copy it to CECI cluster.
# Execute it with: bash compile_Kori.sh
path_kori=/home/daniel/models/Kori-ULB
path_kori_subroutines=$path_kori/subroutines
path_exe=$path_kori/exe/thwaites/sigma_oce1

path_lemaitre4=/globalscratch/ulb/glaciol/dmoreno/Kori-ULB/exe/sigma_oce1


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
# Loop over each file in the directory
for file in "$path_exe"/*.m; do
    # Check if it's a file (not a directory)
    if [ -f "$file" ]; then
        # Extract the base name of the file (without path)
        file_name=$(basename "$file")
        
        # Set the exe_name based on your naming conventions
        exe_name="${file_name%.*}" # removes file extension
        
        # Print the status for debugging
        echo "Compiling: $file_name"
        echo "Exe_name : $exe_name"

        # Run the command
        mcc -m "$file" -a "$path_kori/KoriModel.m" -a "$path_kori_subroutines" -o "$exe_name"
        
        # Securely copy the executable to the remote server
        scp "$exe_name" lemaitre4:"$path_lemaitre4"
    fi
done
