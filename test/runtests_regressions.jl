using Test
using Base.Filesystem
include("../src/TJLF.jl")
using ..TJLF


# saturation rule test
directory = "../outputs/tglf_regression/"
tests = readdir(directory)

# 03 is s-alpha geometry
# 05 is ns = 3
excludeFolders = ["tglf03", "tglf05"]
testFolders = [joinpath(directory,item) for item in readdir(directory) if isdir(joinpath(directory,item)) && item ∉ excludeFolders]

for baseDirectory in testFolders
    @testset "$baseDirectory" begin
        
        fileDirectory = joinpath(baseDirectory, "out.tglf.QL_flux_spectrum")
        lines = readlines(fileDirectory)
        (ntype, nspecies, nfield, nky, nmodes) = parse.(Int32, split(lines[4]))
        fluxesFortran = Vector{Float64}()
        for line in lines[7:length(lines)]
            line = split(line)
            if any(occursin.(["m","s"],string(line))) continue end
            for x in line
                push!(fluxesFortran,parse(Float64, string(x)))
            end
        end
        
        inputTJLF = readInput(baseDirectory)
        fluxesJulia = TJLF.run(inputTJLF)
        
        @test isapprox(sum(fluxesJulia),sum(fluxesFortran), rtol=1e-6)
    end

end

println("REGRESSION SUCCESS")