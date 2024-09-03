#!/bin/bash
#SBATCH --nodes 7
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem="32GB"
#SBATCH --output=/share/nas/kelvinw/sbatch/logs/out-slurm_%j.out
#SBATCH --error=/share/nas/kelvinw/sbatch/logs/out-slurm_%j.err
#SBATCH --job-name tclean
#SBATCH --mail-user=kelvin.wandia@manchester.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE


# Start time logging
echo "Job started on:"
date +'%Y-%m-%d %H:%M:%S'
DATE_START=$(date +%s)

# Load necessary modules
echo "Loading openmpi module..."
module load openmpi-2.1.1 || { echo "Failed to load openmpi module"; exit 1; }
echo "MPI Environment Variables:"
env | grep MPI



# # Paths and files
# directory="/state/partition1/kelvinw/m15psrc_working_dir"
# fileToCopy="/share/nas/kelvinw/trial/M15PSRC.ms"
# fileToCheck="$directory/M15PSRC.ms"
# fileToCopy1="/share/nas/kelvinw/casa_vlbi"
# fileToCheck1="$directory/casa_vlbi"

# # Function to copy files and check for success
# copy_file_if_not_exists() {
#     local src=$1
#     local dest=$2

#     if [ -e "$dest" ]; then
#         echo "File or directory $dest already exists. Skipping copy."
#     else
#         mkdir -p "$(dirname "$dest")"
#         echo "Copying from $src to $dest"
#         cp -r "$src" "$dest" || { echo "Failed to copy $src to $dest"; exit 1; }
#         echo "Finished copying $src to $dest"
#     fi
# }

# # Copy files
# copy_file_if_not_exists "$fileToCopy1" "$fileToCheck1"
# copy_file_if_not_exists "$fileToCopy" "$fileToCheck"

# # Change to working directory
# cd "$directory" || { echo "Failed to change directory to $directory"; exit 1; }

### Use this when using ssh-compute-0-XXX ie when in /state/partition1
# FILENAME="/casa_vlbi/data/white_dwarfs_propagated_coords.txt"
# VIS="M15PSRC.ms"
# FIELD="M15PSRC"  
# FILEPATH="trial/"

# Filenames and field
FILENAME="/mnt/casa_vlbi/data/white_dwarfs_propagated_coords.txt"
VIS="/mnt/trial/M15PSRC.ms"
FIELD="M15PSRC"  
FILEPATH="/mnt/trial/"

# Set environment variables for the container
export CONTAINER=/share/nas/kelvinw/singularity_images/casa6_openmpi_2.1.1.sif
export SCRIPT=/mnt/casa_vlbi/imaging.py

# Run the main processing command using Singularity and MPI
echo "Running singularity with mpirun..."
mpirun singularity exec -B /share/nas/kelvinw:/mnt,/state/partition1 \
    "$CONTAINER" python "$SCRIPT" \
    --filename "$FILENAME" \
    --vis "$VIS" \
    --field "$FIELD" \
    --filepath "$FILEPATH" || { echo "Singularity execution failed"; exit 1; } 

# End time logging
echo "Job finished on:"
date +'%Y-%m-%d %H:%M:%S'
DATE_END=$(date +%s)

# Calculate and display total runtime
TOTAL_SECONDS=$((DATE_END - DATE_START))
HOURS=$((TOTAL_SECONDS / 3600))
MINUTES=$(( (TOTAL_SECONDS % 3600) / 60 ))
SECONDS=$((TOTAL_SECONDS % 60))
echo "Total run time: $HOURS Hours $MINUTES Minutes $SECONDS Seconds"
