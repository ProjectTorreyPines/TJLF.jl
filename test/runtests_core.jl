using Test
using TJLF

# Tolerance for comparing TJLF vs TGLF results
const RTOL = 0.0005  # 0.05% tolerance for NT core case

@testset "TJLF vs TGLF NT core comparison" begin
    # Test directory with input.tglf
    test_dir = joinpath(@__DIR__, "test_core")
    
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
    
    QL_weights, firstPass_eigenvalue, secondPass_eigenvalue = tjlf_TM(input_tjlf, satParams, outputHermite; return_both_eigenvalues=true)
    
    # For eigenvalue comparison, use first pass eigenvalues (matches TGLF printed output)
    gammaJulia = firstPass_eigenvalue[:, :, 1]
    freqJulia = firstPass_eigenvalue[:, :, 2]

    # Compare eigenvalues
    println("\nEigenvalue Comparison:")
    for i in 1:length(gamma)  # Use length of TGLF gamma vector
        println("Mode $i: TGLF (γ=$(gamma[i]), ω=$(freq[i])) vs TJLF (γ=$(gammaJulia'[i]), ω=$(freqJulia'[i]))")
        @test isapprox(gamma[i], gammaJulia'[i], rtol=RTOL)  # Using global tolerance
        @test isapprox(freq[i], freqJulia'[i], rtol=RTOL)    # Using global tolerance
    end

    # Compare fluxes
    # For flux calculation, use the appropriate eigenvalues based on the input parameters
    eigenvalue_for_flux = if input_tjlf.ALPHA_QUENCH != 0.0 || input_tjlf.VEXB_SHEAR * input_tjlf.SIGN_IT == 0.0
        # Single pass case: use first pass eigenvalues
        firstPass_eigenvalue[:, :, 1]
    else
        # Spectral shift case: use second pass eigenvalues for flux calculation
        secondPass_eigenvalue[:, :, 1]
    end
    
    if input_tjlf.SAT_RULE == 2 || input_tjlf.SAT_RULE == 3
        most_unstable_gamma_first_pass = firstPass_eigenvalue[1, :, 1]
        vzf_out_first_pass, kymax_out_first_pass, jmax_out_first_pass = TJLF.get_zonal_mixing(input_tjlf, satParams, most_unstable_gamma_first_pass)
        
        println("\nZonal mixing parameters for flux calculation:")
        println("  vzf_out_first_pass = ", vzf_out_first_pass)
        println("  kymax_out_first_pass = ", kymax_out_first_pass)
        println("  jmax_out_first_pass = ", jmax_out_first_pass)
        
        QL_flux_out, flux_out = sum_ky_spectrum(input_tjlf, satParams, eigenvalue_for_flux, QL_weights; 
                                                vzf_out_param=vzf_out_first_pass, 
                                                kymax_out_param=kymax_out_first_pass, 
                                                jmax_out_param=jmax_out_first_pass)
    else
        QL_flux_out, flux_out = sum_ky_spectrum(input_tjlf, satParams, eigenvalue_for_flux, QL_weights)
    end
    
    juliaEnergy = sum(QL_flux_out; dims=1)[1, :, 2]
    fortranEnergy = fluxesFortran[2, :]

    println("\nFlux Comparison:")
    println("TGLF Energy Fluxes: ", fortranEnergy)
    println("TJLF Energy Fluxes: ", juliaEnergy)
    @test isapprox(juliaEnergy, fortranEnergy, rtol=RTOL)  # Using global tolerance
end 