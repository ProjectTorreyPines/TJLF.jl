using Pkg
Pkg.activate("..")
Pkg.build("TJLF")
using TJLF
using BenchmarkTools
using LinearAlgebra
BLAS.set_num_threads(1)

inputTJLFVector = Vector{Main.TJLF.InputTJLF}(undef, 7)

baseDirectory = "../outputs/TIM_case/case2/"
inputTJLFVector[1] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case3/"
inputTJLFVector[2] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case4/"
inputTJLFVector[3] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case5/"
inputTJLFVector[4] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case6/"
inputTJLFVector[5] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case7/"
inputTJLFVector[6] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case8/"
inputTJLFVector[7] = readInput(baseDirectory*"input.tglf")

TJLF.run_tjlf(inputTJLFVector)

@btime TJLF.run_tjlf(inputTJLFVector) |> x -> nothing
