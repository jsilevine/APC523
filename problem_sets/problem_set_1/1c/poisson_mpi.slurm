#!/bin/bash
#SBATCH --job-name=poisson_omp     # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=2               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1G         # memory per cpu-core (4G is default)
#SBATCH --time=00:00:20          # total run time limit (HH:MM:SS)

module purge
module load openmpi/gcc/4.1.0
module load R/4.0.5

## cleanup previous run
rm time.csv
rm final_phi.csv
rm phi_128.png

mpirun --use-hwthread-cpus -np 2 ./a.out 16 1e-6 2
Rscript p1a_fig_gen.R
