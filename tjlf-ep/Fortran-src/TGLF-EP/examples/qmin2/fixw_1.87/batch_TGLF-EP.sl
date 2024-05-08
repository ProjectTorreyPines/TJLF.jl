#!/bin/bash -l

#SBATCH -p debug
#SBATCH -A m808
#SBATCH -o ./run.out
#SBATCH -e ./run.err
#SBATCH -N 40
#SBATCH -t 00:30:00
#SBATCH -J qmin2_87
#SBATCH --constraint=haswell

srun -n  1250 /global/u2/b/bassem/gacode_cori_7-2018/TGLFEP/TGLFEP_driver
