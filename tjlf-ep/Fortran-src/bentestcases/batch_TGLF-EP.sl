#!/bin/bash

#SBATCH -p medium
#SBATCH -n 6
#SBATCH -t 24:00:00
#SBATCH -o ./run.out
#SBATCH -e ./run.err
#SBATCH -J TGLF-EP_OMFIT
#SBATCH --export=all


export PATH="/home/towlej/gacode/tglf/bin/tglf:$PATH"
export GACODE_PLATFORM=OMEGA_INTEL
export GACODE_ROOT=/home/towlej/gacode
export GACODE_ADD_ROOT=/home/towlej/.julia/dev/TJLF.jl/tjlf-ep/Fortran-src/TGLF-EP
. ${GACODE_ROOT}/shared/bin/gacode_setup
. ${GACODE_ROOT}/platform/env/env.$GACODE_PLATFORM

# export PATH="/home/towlej/gacode/tglf/bin/tglf:$PATH"
# export GACODE_PLATFORM=OMEGA_INTEL
# export GACODE_ROOT=/home/towlej/gacode
# export GACODE_ADD_ROOT=/home/towlej/gacode_add/TGLF-EP
# . ${GACODE_ROOT}/shared/bin/gacode_setup
# . ${GACODE_ROOT}/platform/env/env.$GACODE_PLATFORM


#module purge

env > env.txt

mpirun -n 6 /home/towlej/.julia/dev/TJLF.jl/tjlf-ep/Fortran-src/TGLF-EP/TGLFEP_driver

