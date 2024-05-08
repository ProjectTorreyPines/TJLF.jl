#!/bin/bash -l

#SBATCH -p debug
#SBATCH -A m808
#SBATCH -o ./run.out
#SBATCH -e ./run.err
#SBATCH -N 37
#SBATCH -t 00:30:00
#SBATCH -J qmin2_D3D
#SBATCH --constraint=knl

srun -n  2500 $TGLFEP_DIR/TGLFEP_driver
