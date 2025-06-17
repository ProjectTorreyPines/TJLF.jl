using Test
using TJLF
#using BenchmarkTools  #for using @btime 


@show Threads.nthreads()


inputTJLFVector = Vector{Main.TJLF.InputTJLF}(undef, 7)
directory = joinpath(dirname(@__DIR__), "tglf_regression")

baseDirectory = joinpath(directory, "tglf01/")
inputTJLFVector[1] = readInput(baseDirectory*"input.tglf")
baseDirectory = joinpath(directory, "tglf02/")
inputTJLFVector[2] = readInput(baseDirectory*"input.tglf")
baseDirectory = joinpath(directory, "tglf04/")
inputTJLFVector[3] = readInput(baseDirectory*"input.tglf")
baseDirectory = joinpath(directory, "tglf05/")
inputTJLFVector[4] = readInput(baseDirectory*"input.tglf")
baseDirectory = joinpath(directory, "tglf06/")
inputTJLFVector[5] = readInput(baseDirectory*"input.tglf")
baseDirectory = joinpath(directory, "tglf08/")
inputTJLFVector[6] = readInput(baseDirectory*"input.tglf")
baseDirectory = joinpath(directory, "tglf10/")
inputTJLFVector[7] = readInput(baseDirectory*"input.tglf")


# run several times to avoid confusion with a compilation time
for i in 1:3

@time   TJLF.run_tjlf(inputTJLFVector) 

@time   TJLF.run_tjlf(inputTJLFVector[7]) 
println("Iteration $i completed -------")
end

println("")