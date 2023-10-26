using Test
using Base.Filesystem
include("../src/TJLF.jl")
using ..TJLF


# saturation rule test
directory = "../outputs/tglf_regression/"
tests = readdir(directory)

excludeFolders = ["tglf05"]
testFolders = [joinpath(directory,item) for item in readdir(directory) if isdir(joinpath(directory,item)) && item ∉ excludeFolders]

for baseDirectory in testFolders
    
    inputTJLF = readInput(baseDirectory)
    println(dir_name)
    fluxes = TJLF.run(inputTJLF)
    
    @assert !isnan(sum(fluxes))

end

println("REGRESSION SUCCESS")