#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH --mem="16GB"
#SBATCH --output=/share/nas/kelvinw/sbatch/logs/out-slurm_%j.out
#SBATCH --error=/share/nas/kelvinw/sbatch/logs/out-slurm_%j.err
#SBATCH --job-name tclean
#SBATCH --mail-user=kelvin.wandia@manchester.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE


# Start time logging
echo "Job started on:"
date +'%Y-%m-%d %H:%M:%S'
DATE_START=$(date +%s)

# Load necessary modules #SBATCH -w compute-0-116
echo "Loading openmpi module..."
module load openmpi-4.0.4 || { echo "Failed to load openmpi module"; exit 1; }
echo "MPI Environment Variables:"
env | grep MPI

export OMPI_MCA_btl=vader,self
export OMPI_MCA_orte_base_help_aggregate=0


# Set environment variables for the container
export CONTAINER=/share/nas/kelvinw/singularity_images/casa6_openmpi_4.0.4.sif
export SCRIPT=/mnt/casa_vlbi/main.py


# Run the main processing command using Singularity and MPI
echo "Running singularity with mpirun..."
mpirun singularity exec -B /share/nas/kelvinw:/mnt,/state/partition1 \
    "$CONTAINER"  python "$SCRIPT" \

    
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
