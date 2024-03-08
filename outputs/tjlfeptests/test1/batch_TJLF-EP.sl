#!/bin/bash

#SBATCH -p long
#SBATCH -n 2
#SBATCH -t 24:00:00

#SBATCH -o ./run.out
#SBATCH -e ./run.err
#SBATCH -J
#SBATCH --export=all

env > env.txt

mpiexec -n 2 