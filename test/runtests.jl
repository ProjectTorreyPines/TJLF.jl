using Test
using Base.Filesystem
include("../src/TJLF.jl")
using ..TJLF


# saturation rule test
satRuleDirectory = "../outputs/test_SAT_rules/"
tests = readdir(satRuleDirectory)
for dir_name in tests
    if dir_name == ".DS_Store" continue end
    baseDirectory = satRuleDirectory*dir_name*"/"

    @testset "$baseDirectory" begin

        fileDirectory = baseDirectory * "out.tglf.QL_flux_spectrum"
        lines = readlines(fileDirectory)
        # 5, 2, 3, 21, 2
        (ntype, nspecies, nfield, nky, nmodes) = parse.(Int32, split(lines[4]))
        ql = Vector{Float64}()
        for line in lines[7:length(lines)]
            line = split(line)
            if any(occursin.(["m","s"],string(line))) continue end

            for x in line
                push!(ql,parse(Float64, string(x)))
            end

        end
        # 5, 21, 2, 3, 2
        QLw = reshape(ql, (ntype, nky, nmodes, nfield, nspecies))
        QL_data = permutedims(QLw,(4,5,3,2,1))
        # QL = AxisArray(QL_data, Axis{:nky}(1:nky), Axis{:nmodes}(1:nmodes), Axis{:nspecies}(1:nspecies), Axis{:nfield}(1:nfield), Axis{:ntype}(("particle","energy","toroidal_stress","parrllel_stress","exchange_stress")))

        particle_QL = QL_data[:, :, :, :, 1]
        energy_QL = QL_data[:, :, :, :, 2]
        toroidal_stress_QL = QL_data[:, :, :, :, 3]
        parallel_stress_QL = QL_data[:, :, :, :, 4]
        exchange_QL = QL_data[:, :, :, :, 5]

        # Read spectral shift and ave_p0 (only needed for SAT0)
        fileDirectory = baseDirectory * "out.tglf.spectral_shift_spectrum"
        lines = readlines(fileDirectory)
        kx0_e = Vector{Float64}()
        for line in lines[6:length(lines)]
                push!(kx0_e,parse(Float64, line))
        end

        fileDirectory = baseDirectory * "out.tglf.ave_p0_spectrum"
        lines = readlines(fileDirectory)
        ave_p0 = Vector{Float64}()
        for line in lines[4:length(lines)]
                push!(ave_p0,parse(Float64, line))
        end

        # Read scalar saturation parameters
        fileDirectory = baseDirectory * "out.tglf.scalar_saturation_parameters"
        lines = readlines(fileDirectory)
        inputComparison = Dict()
        for line in lines[1:length(lines)]
            line = split(line, "\n")
            #### no idea why this is here
            if any(occursin.(["!", "UNITS", "SAT_RULE", "XNU_MODEL", "ETG_FACTOR", "R_unit", "ALPHA_ZF", "RULE"],string(line)))
                if contains(line[1],"R_unit")
                    line = split(line[1]," = ")
                    R_unit = parse(Float64, strip(line[2]))
                end
                continue
            end
                
            line = split(line[1]," = ")
            line .= strip.(line)
            inputComparison[string(line[1])] = parse(Float64, line[2])
        end

        # Get ky spectrum
        fileDirectory = baseDirectory * "out.tglf.ky_spectrum"
        lines = readlines(fileDirectory)
        ky_spect = Array{Float64}(undef, 0)
        for line in lines[3:length(lines)]
                push!(ky_spect,parse(Float64, line))
        end

        # Get eigenvalue spectrum
        fileDirectory = baseDirectory * "out.tglf.eigenvalue_spectrum"
        lines = readlines(fileDirectory)
        lines = split(join(lines[3:length(lines)]))
        lines = [parse(Float64, l) for l in lines]

        gamma = []
        freq = []
        for k in 1:nmodes
            push!(gamma, lines[2k-1:2*nmodes:end])
            push!(freq, lines[2k-1:2*nmodes:end])
        end


        gammas = hcat(gamma...)
        R_unit = ones(size(gammas)) * R_unit

        # Get potential spectrum
        fileDirectory = baseDirectory * "out.tglf.field_spectrum"
        lines = readlines(fileDirectory)

        columns = split.(lines[2],",")
        nc = length(columns)

        lines = split(join(lines[6:length(lines)]))

        tmpdict = Dict()
        for (ik, k) in enumerate(columns)
            tmp = []
            for nm in 1:nmodes
                push!(tmp, parse.(Float64,lines[ik - 3 + nm * nc:nc*nmodes:end]))
            end
            tmpdict[k] = tmp
        end
        potentialTmp = tmpdict[columns[length(columns)]]
        potential = hcat(potentialTmp...)


        fileDirectory = baseDirectory * "out.tglf.gbflux"
        lines = readlines(fileDirectory)
        width::Integer = round(length(split(lines[1]))/4)
        fluxes = transpose(reshape(parse.(Float64,split(lines[1])), (width,4)))

        #******************************************************************************#************************
        # Read input.tglf
        #******************************************************************************#************************
        
        inputTJLF = readInput(baseDirectory)

        #*******************************************************************************************************
        #   start running stuff
        #*******************************************************************************************************

        satParams = get_sat_params(inputTJLF)

        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
        @test isapprox(ky_spect, inputTJLF.KY_SPECTRUM, rtol=1e-6)

        kx0epy = xgrid_functions_geo(inputTJLF, satParams, Matrix(gammas'))
        @test isapprox(kx0epy, kx0_e, rtol=1e-3)
        @test isapprox(inputComparison["SAT_geo0_out"], satParams.SAT_geo0, rtol=1e-6)
        @test isapprox(inputComparison["SAT_geo1_out"], satParams.SAT_geo1, rtol=1e-6)
        @test isapprox(inputComparison["SAT_geo2_out"], satParams.SAT_geo2, rtol=1e-6)
        @test isapprox(R_unit[1, 1], satParams.R_unit,  rtol=1e-6)
        @test isapprox(inputComparison["Bt0_out"], satParams.Bt0, rtol=1e-6)
        @test isapprox(inputComparison["grad_r0_out"], satParams.grad_r0, rtol=1e-6)

        if inputTJLF.VEXB_SHEAR != 0.0
            @test isapprox(inputComparison["B_geo0_out"], satParams.B_geo[1], rtol=1e-6)
        end

        QL_flux_out, flux_out = sum_ky_spectrum(inputTJLF, satParams, Matrix(gammas'), QL_data)
        # julia_sat1 = sum(sum(sat_1["energy_flux_integral"], dims=3)[:,:,1], dims=1)[1,:]
        julia_sat1 = sum(QL_flux_out, dims = 1)[1,:,2]
        expected_sat1 = fluxes[2,:]
        @test isapprox(julia_sat1, expected_sat1, rtol=1e-3)
    end
end














satRuleDirectory = "../outputs/test_kygrid/"
tests = readdir(satRuleDirectory)
for dir_name in tests
    if dir_name == ".DS_Store" continue end
    baseDirectory = satRuleDirectory*dir_name*"/"
    @testset "$baseDirectory" begin

        # Get ky spectrum
        fileDirectory = baseDirectory * "out.tglf.ky_spectrum"
        lines = readlines(fileDirectory)
        ky_spect = Array{Float64}(undef, 0)
        nky = parse(Int, strip(lines[2]))
        for line in lines[3:length(lines)]
            push!(ky_spect,parse(Float64, line))
        end

        #******************************************************************************#************************
        # Read input.tglf
        #******************************************************************************#************************

        inputTJLF = readInput(baseDirectory)

        #*******************************************************************************************************
        #   start running stuff
        #*******************************************************************************************************
        
        satParams = get_sat_params(inputTJLF)
        Julia_ky_spect = get_ky_spectrum(inputTJLF, satParams.grad_r0)
        @test isapprox(Julia_ky_spect, ky_spect, rtol=1e-6)

    end

end












testTMDirectory = "../outputs/test_TM/"
tests = readdir(testTMDirectory)
for dir_name in tests
    if dir_name == ".DS_Store" continue end
    baseDirectory = testTMDirectory*dir_name*"/"

    @testset "$baseDirectory" begin

        # get particle, energy, exchange fluxes
        fileDirectory = baseDirectory * "out.tglf.QL_flux_spectrum"
        lines = readlines(fileDirectory)
        (ntype, nspecies, nfield, nky, nmodes) = parse.(Int32, split(lines[4]))
        ql = Vector{Float64}()
        for line in lines[7:length(lines)]
            line = split(line)
            if any(occursin.(["m","s"],string(line))) continue end
        
            for x in line
                push!(ql,parse(Float64, string(x)))
            end
        
        end
        QLw = reshape(ql, (ntype, nky, nmodes, nfield, nspecies))
        QL_data = permutedims(QLw,(4,5,3,2,1))

        # Get eigenvalue spectrum
        fileDirectory = baseDirectory * "out.tglf.eigenvalue_spectrum"
        lines = readlines(fileDirectory)
        lines = split(join(lines[3:length(lines)]))
        lines = [parse(Float64, l) for l in lines]

        gamma = Vector{Vector{Float64}}()
        freq = Vector{Vector{Float64}}()
        for k in 1:nmodes
            push!(gamma, lines[2k-1:2*nmodes:end])
            push!(freq, lines[2k:2*nmodes:end])
        end
        gamma = vcat(gamma...)
        freq = vcat(freq...)

        # Get the integral flux
        fileDirectory = baseDirectory * "out.tglf.gbflux"
        lines = readlines(fileDirectory)
        width::Integer = round(length(split(lines[1]))/4)
        fluxesFortran = transpose(reshape(parse.(Float64,split(lines[1])), (width,4)))

        #******************************************************************************************************
        # Read input.tglf
        #******************************************************************************************************

        inputTJLF = readInput(baseDirectory)
        inputTJLF2 = readInput(baseDirectory)

        #*******************************************************************************************************
        #   start running stuff
        #*******************************************************************************************************

        outputHermite = gauss_hermite(inputTJLF)
        satParams = get_sat_params(inputTJLF)
        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)

        fluxes, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
        gammaJulia = eigenvalue[:,:,1]
        freqJulia = eigenvalue[:,:,2]

        for i in eachindex(fluxes[1,:,:,:,:])
            @test isapprox(QL_data[i], (fluxes[1,:,:,:,:])[i], atol=1e-6)
        end
        for i in eachindex(gammaJulia)
            @test isapprox(gamma[i],gammaJulia'[i], atol=1e-6)
            @test isapprox(freq[i],freqJulia'[i], atol=1e-6)
        end

        QL_flux_out, flux_out = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], fluxes)
        juliaEnergy = sum(QL_flux_out, dims = 1)[1,:,2]
        fortranEnergy = fluxesFortran[2,:]
        @test isapprox(juliaEnergy, fortranEnergy, rtol=1e-3)


        #*******************************************************************************************************0
        #   test that providing widths yields same results (NO LONGER VALID AFTER NEW QUICK/ITERATIVE EIGENSOLVER METHOD)
        #*******************************************************************************************************

        # inputTJLF2.KY_SPECTRUM .= inputTJLF.KY_SPECTRUM
        # inputTJLF2.WIDTH_SPECTRUM .= inputTJLF.WIDTH_SPECTRUM
        # inputTJLF2.GAMMA_SPECTRUM .= eigenvalue[inputTJLF.NMODES,:,1]
        # inputTJLF2.FIND_WIDTH = false

        # fluxes2, eigenvalue2 = tjlf_TM(inputTJLF2, satParams, outputHermite)

        # for i in eachindex(fluxes)
        #     @test isapprox(fluxes[i], fluxes2[i],atol=1e-8)
        # end

    end

end



println("SUCCESS")