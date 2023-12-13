using Pkg
Pkg.activate("..")
Pkg.build("TJLF")
using TJLF
using BenchmarkTools
using LinearAlgebra
BLAS.set_num_threads(1)

N_cases = 7
inputTJLFVector = Vector{Main.TJLF.InputTJLF}(undef, N_cases)

for case in 1:N_cases
    inputTJLFVector[case] = "../outputs/TIM_case/case$(case+1)/input.tglf"
end
TJLF.run_tjlf(inputTJLFVector)

@btime TJLF.run_tjlf(inputTJLFVector) |> x -> nothing
