#!/bin/bash
#SBATCH --job-name=AIMD 
#SBATCH --output=%x_id_%j.out  # Output file based on job name and job ID
#SBATCH --error=%x_id_%j.err
#SBATCH --time=5-00:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --mem=30G 

# Load environment
module load openmpi/4.1.6
module load cp2k/2023.2

export OMP_NUM_THREADS=1
# Calculate total number of processes
TOTAL_TASKS=$((SLURM_NNODES * SLURM_NTASKS_PER_NODE))
echo "Running job with name ${SLURM_JOB_NAME} and ID ${SLURM_JOB_ID}"


start_time=$(date +%s)
echo "Job started at $(date)"

# --- Your code here ---
srun -n $TOTAL_TASKS cp2k.psmp -i input.inp -o ${SLURM_JOB_NAME}_id_${SLURM_JOB_ID}.log
# --- End of code ---

end_time=$(date +%s)
echo "Job ended at $(date)"
# Calculate and log the total time used in hours
total_time=$((end_time - start_time))
total_hours=$(echo "scale=2; $total_time / 3600" | bc)
echo "Total time used: $total_hours hours"
