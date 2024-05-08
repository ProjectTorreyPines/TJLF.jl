#!/bin/bash

#SBATCH -p medium
#SBATCH -n 3
#SBATCH -t 24:00:00

#SBATCH -o ./run.out
#SBATCH -e ./run.err
#SBATCH -J TGLF-EP_OMFIT
#SBATCH --export=all

fc
env > env.txt

mpirun -n 3 /home/agnewb/gacode_add/TGLF-EP/TGLFEP_driver
