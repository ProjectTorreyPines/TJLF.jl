using Test
using TJLF

# Tolerance for comparing TJLF vs TGLF results
const RTOL = 0.00001  # 1% tolerance for EM case

@testset "TJLF vs TGLF electromagnetic comparison" begin
    # Test directory with input.tglf
    test_dir = joinpath(@__DIR__, "test_EM")
    
    # Get particle, energy, exchange fluxes from TGLF output
    fileDirectory = joinpath(test_dir, "out.tglf.QL_flux_spectrum")
    lines = readlines(fileDirectory)
    (ntype, nspecies, nfield, nky, nmodes) = parse.(Int32, split(lines[4]))
    ql = Vector{Float64}()
    for line in lines[7:length(lines)]
        line = split(line)
        if any(occursin.(["m", "s"], string(line)))
            continue
        end
        for x in line
            push!(ql, parse(Float64, string(x)))
        end
    end
    QLw = reshape(ql, (ntype, nky, nmodes, nfield, nspecies))
    QL_data = permutedims(QLw, (4, 5, 3, 2, 1))

    # Get eigenvalue spectrum from TGLF output
    fileDirectory = joinpath(test_dir, "out.tglf.eigenvalue_spectrum")
    lines = readlines(fileDirectory)
    lines = lines[3:end]  # Skip header lines
    gamma = Float64[]
    freq = Float64[]
    for line in lines
        if !isempty(strip(line))
            nums = parse.(Float64, split(line))
            # Each line contains gamma1, freq1, gamma2, freq2
            # We only want the first pair (gamma1, freq1)
            push!(gamma, nums[1])
            push!(freq, nums[2])
        end
    end

    # Get the integral flux from TGLF output
    fileDirectory = joinpath(test_dir, "out.tglf.gbflux")
    lines = readlines(fileDirectory)
    width::Integer = round(length(split(lines[1])) / 4)
    fluxesFortran = transpose(reshape(parse.(Float64, split(lines[1])), (width, 4)))

    # Run TJLF
    input_tjlf = readInput(joinpath(test_dir, "input.tglf"))
    satParams = get_sat_params(input_tjlf)
    input_tjlf.KY_SPECTRUM .= get_ky_spectrum(input_tjlf, satParams.grad_r0)
    
    outputHermite = gauss_hermite(input_tjlf)
    fluxes, eigenvalue = tjlf_TM(input_tjlf, satParams, outputHermite)
    gammaJulia = eigenvalue[:, :, 1]
    freqJulia = eigenvalue[:, :, 2]

    # Compare eigenvalues
    println("\nEigenvalue Comparison:")
    for i in 1:length(gamma)  # Use length of TGLF gamma vector
        println("Mode $i: TGLF (γ=$(gamma[i]), ω=$(freq[i])) vs TJLF (γ=$(gammaJulia'[i]), ω=$(freqJulia'[i]))")
        @test isapprox(gamma[i], gammaJulia'[i], rtol=RTOL)  # Using global tolerance
        @test isapprox(freq[i], freqJulia'[i], rtol=RTOL)    # Using global tolerance
    end

    # Compare fluxes
    QL_flux_out, flux_out = sum_ky_spectrum(input_tjlf, satParams, eigenvalue[:, :, 1], fluxes)
    juliaEnergy = sum(QL_flux_out; dims=1)[1, :, 2]
    fortranEnergy = fluxesFortran[2, :]

    println("\nFlux Comparison:")
    println("TGLF Energy Fluxes: ", fortranEnergy)
    println("TJLF Energy Fluxes: ", juliaEnergy)
    @test isapprox(juliaEnergy, fortranEnergy, rtol=RTOL)  # Using global tolerance
end 