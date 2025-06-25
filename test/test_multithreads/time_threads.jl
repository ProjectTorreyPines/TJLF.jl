using Test
using TJLF

n_cases = 8
@show Threads.nthreads()

# this case was picked because it takes the longest
case = readInput(joinpath(dirname(@__DIR__), "tglf_regression", "tglf10", "input.tglf"))
inputTJLFVector = TJLF.InputTJLF[case for _ in 1:n_cases]

# single run just to ensure things get precompiled
TJLF.run_tjlf(inputTJLFVector[1])

# serial run
@time TJLF.run_tjlf(inputTJLFVector[1])

# parallel run
@time TJLF.run_tjlf(inputTJLFVector)

println("")