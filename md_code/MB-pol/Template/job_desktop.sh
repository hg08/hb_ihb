#!/bin/bash

#SBATCH --job-name=125h2o                  # Job name
#SBATCH --output=%x_id_%j.out                # Output file based on job name and job ID
#SBATCH --error=%x_id_%j.err                 # Error file based on job name and job ID
#SBATCH --time=INFINITE                   # Time limit
#SBATCH --partition=debug                 # Partition name

echo "Running job with name ${SLURM_JOB_NAME} and ID ${SLURM_JOB_ID}"
start_time=$(date +%s)
echo "Job started at $(date)"

# --- Your code here ---

LAMMPS_HOME=${LAMMPS_HOME}

MPI_PROC=1        # number of mpi processes
OMP_NUM_THREADS=4 # number of openmp threads per mpi process

# optional command line arguments to override the default mpi and openmp settings
if [ $# -eq 2 ]; then
    MPI_PROC=$1
    OMP_NUM_THREADS=$2
fi


echo "Running MBX+LAMMPS using MPI=$MPI_PROC and OMP=$OMP_NUM_THREADS"

export OMP_NUM_THREADS
mpirun -np ${MPI_PROC} $LAMMPS_HOME/src/lmp_mpi_mbx -in lammps.in

if [ -s log.lammps ] && grep -q "Total wall time" log.lammps; then
    date | tee -a performance.txt
    echo "<$(hostname -s)>" | tee -a performance.txt
    grep 'ns/day' log.lammps | tee -a performance.txt
    grep '% CPU use' log.lammps | tee -a performance.txt
    grep 'Total wall time' log.lammps | tee -a performance.txt
    echo "-----------------------------" | tee -a performance.txt
fi
# --- End of your code ---

end_time=$(date +%s)
echo "Job ended at $(date)"
# Calculate and log the total time used in hours
total_time=$((end_time - start_time))
total_hours=$(echo "scale=2; $total_time / 3600" | bc)
echo "Total time used: $total_hours hours"
