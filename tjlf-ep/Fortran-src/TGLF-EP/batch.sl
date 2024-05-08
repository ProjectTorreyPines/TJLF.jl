#!/bin/bash -l

#SBATCH -p debug
#SBATCH -A m808
#SBATCH -o ./run.out
#SBATCH -e ./run.err
#SBATCH -N 42
#SBATCH -t 00:30:00
#SBATCH -J TGLFEP

cd $SLURM_SUBMIT_DIR
srun -n 1000 ./TGLFEP_driver
