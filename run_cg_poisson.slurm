#!/bin/bash

#SBATCH --exclusive
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH --mem=0
#SBATCH --account=Education-EEMCS-Courses-WI4450

# run CG for a n^3 Poisson problem
n=128
export OMP_NUM_THREADS=8
srun ./main_cg_poisson.x 1024 128 128
# srun ./main_benchmarks.x $n >> benchmarks-weak-${OMP_NUM_THREADS}.out
