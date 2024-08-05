#!/bin/bash
#SBATCH --job-name=SurfaceTension  # Job name 
#SBATCH --output=%x_id_%j.out  # Output file based on job name and job ID
#SBATCH --error=%x_id_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G 

# Load environment
module load mamba gcc
source activate ihb


start_time=$(date +%s)
echo "Job started at $(date)"

# --- Your code here ---
./2_batchSurfaceTension.sh
# --- End of code ---

end_time=$(date +%s)
echo "Job ended at $(date)"
# Calculate and log the total time used in hours
total_time=$((end_time - start_time))
total_hours=$(echo "scale=2; $total_time / 3600" | bc)
echo "Total time used: $total_hours hours"
