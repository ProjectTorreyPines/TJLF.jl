using Pkg
path = joinpath(@__DIR__, "..")
Pkg.activate(path)
Pkg.build("TJLF")
using TJLF
using BenchmarkTools
using LinearAlgebra
BLAS.set_num_threads(1)

N_cases = 7
inputTJLFVector = TJLF.InputTJLF[readInput("$path/outputs/TIM_case/case$(case+1)/input.tglf") for case in 1:N_cases]
TJLF.run_tjlf(inputTJLFVector)

@btime TJLF.run_tjlf(inputTJLFVector) |> x -> nothing