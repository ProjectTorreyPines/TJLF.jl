#!/bin/bash -l

#SBATCH -p debug
#SBATCH -A m808
#SBATCH -o ./run.out
#SBATCH -e ./run.err
#SBATCH -N 2
#SBATCH -t 00:30:00
#SBATCH -J trace_EP
#SBATCH --constraint=haswell

srun -n  50 $TGLFEP_DIR/TGLFEP_driver
