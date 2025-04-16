using Test
using TJLF

# tglf regression test
directory = joinpath(@__DIR__, "tglf_regression")
tests = readdir(directory)

# 03 is s-alpha geometry, not implemented
excludeFolders = ["tglf03"]
testFolders = [joinpath(directory, item) for item in readdir(directory) if isdir(joinpath(directory, item)) && item ∉ excludeFolders]

for baseDirectory in testFolders
    @testset "$baseDirectory" begin

        fileDirectory = joinpath(baseDirectory, "out.tglf.gbflux")
        lines = readlines(fileDirectory)
        fluxesFortran = parse.(Float64, split(lines[1]))

        inputTJLF = readInput(joinpath(baseDirectory,"input.tglf"))
        fluxesJulia = sum(TJLF.run_tjlf(inputTJLF); dims=1)[1, :, :]
        
        # Test saving InputTJLF and running on it
        temp_dir = mktempdir(baseDirectory)
        TJLF.save(inputTJLF, joinpath(temp_dir, "test_input_tjlf"))
        inputTJLF = readInput(joinpath(temp_dir, "test_input_tjlf"))
        fluxesJulia = sum(TJLF.run_tjlf(inputTJLF); dims=1)[1, :, :]

        for i in 1:3*inputTJLF.NS
            @test isapprox(sum(fluxesJulia[i]), sum(fluxesFortran[i]), atol=1e-2)
        end
        for i in 3*inputTJLF.NS+1:4*inputTJLF.NS
            @test isapprox(sum(fluxesJulia[i+inputTJLF.NS]), sum(fluxesFortran[i]), atol=1e-2)
        end
    end
end
