using Test
using TJLF

kygridDirectory = joinpath(@__DIR__, "test_kygrid")
tests = readdir(kygridDirectory)
for dir_name in tests
    if dir_name == ".DS_Store"
        continue
    end
    baseDirectory = joinpath(kygridDirectory, dir_name)
    @testset "$baseDirectory" begin
        # Get ky spectrum
        fileDirectory = joinpath(baseDirectory, "out.tglf.ky_spectrum")
        lines = readlines(fileDirectory)
        ky_spect = Array{Float64}(undef, 0)
        nky = parse(Int, strip(lines[2]))
        for line in lines[3:length(lines)]
            push!(ky_spect, parse(Float64, line))
        end

        #******************************************************************************#************************
        # Read input.tglf
        #******************************************************************************#************************

        inputTJLF = readInput(joinpath(baseDirectory, "input.tglf"))

        #*******************************************************************************************************
        #   start running stuff
        #*******************************************************************************************************

        satParams = get_sat_params(inputTJLF)
        Julia_ky_spect = get_ky_spectrum(inputTJLF, satParams.grad_r0)
        @test isapprox(Julia_ky_spect, ky_spect, rtol=1e-6)

    end
end