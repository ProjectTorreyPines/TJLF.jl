function run(inputTJLF::InputTJLF)
    checkInput(inputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    QL_weights, firstPass_eigenvalue, secondPass_eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite; return_both_eigenvalues=true)
    
    # Use appropriate eigenvalues for flux calculation
    # For two-pass case (ALPHA_QUENCH=0 and VEXB_SHEAR≠0), use second pass eigenvalues
    # For single-pass case, first and second pass are identical
    vexb_shear_s = inputTJLF.VEXB_SHEAR * inputTJLF.SIGN_IT
    if inputTJLF.ALPHA_QUENCH == 0.0 && vexb_shear_s != 0.0
        # Two-pass case: use second pass eigenvalues for flux calculation
        eigenvalue_for_flux = secondPass_eigenvalue
    else
        # Single-pass case: use first pass eigenvalues
        eigenvalue_for_flux = firstPass_eigenvalue
    end
    
    # Calculate fluxes with proper zonal mixing for SAT_RULE 2/3
    if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3
        most_unstable_gamma_first_pass = firstPass_eigenvalue[1, :, 1]
        vzf_out_first_pass, kymax_out_first_pass, jmax_out_first_pass = get_zonal_mixing(inputTJLF, satParams, most_unstable_gamma_first_pass)
        QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue_for_flux[:, :, 1], QL_weights; 
                                                    vzf_out_param=vzf_out_first_pass, 
                                                    kymax_out_param=kymax_out_first_pass, 
                                                    jmax_out_param=jmax_out_first_pass)
    else
        QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue_for_flux[:, :, 1], QL_weights)
    end
    
    # Return first pass eigenvalues for output (matching TGLF behavior)
    return (QL_weights=QL_weights, eigenvalue=firstPass_eigenvalue, QL_flux_out=QL_flux_out, flux_spectrum=flux_spectrum)
end

function run(inputTGLF::InputTGLF)
    inputTJLF = InputTJLF{Float64}(inputTGLF)
    return run(inputTJLF)
end

"""
    run_tjlf(inputTJLF::InputTJLF)

parameters:
    input_tjlfs::InputTJLF                  - InputTJLF struct

outputs:
    outputs                                 - fluxes (field, species, type)

description:
    Runs TJLF on a single InputTJLF struct, during the run, will save the width spectrum and eigenvalue spectrum to the InputTJLF struct.
    If you want to use these widths and eigenvalues in future runs, there is a flag: FIND_WIDTH and FIND_EIGEN that you set to false
"""
function run_tjlf(inputTJLF::InputTJLF)
    return run(inputTJLF).QL_flux_out
end

"""
    run_tjlf(input_tjlfs::Vector{InputTJLF{T}}) where {T<:Real}

parameters:
    input_tjlfs::Vector{InputTJLF{T}}       - vector of InputTJLF structs

outputs:
    outputs                                 - vector of fluxes (field, species, type)

description:
    Runs TJLF on a vector of InputTJLF structs, during the run, will save the width spectrum and eigenvalue spectrum to the InputTJLF struct.
    If you want to use these widths and eigenvalues in future runs, there is a flag: FIND_WIDTH and FIND_EIGEN that you set to false
"""
function run_tjlf(input_tjlfs::Vector{InputTJLF{T}}) where {T<:Real}
    checkInput(input_tjlfs)
    outputs = Vector{Array{Float64,3}}(undef, length(input_tjlfs))
    Threads.@threads for idx in eachindex(input_tjlfs)
        outputs[idx] = TJLF.run_tjlf(input_tjlfs[idx])
    end
    return outputs
end

# (field, species, type)
# type: (particle, energy, torodial stress, parallel stress, exchange)

Qe(QL_flux_out::Array{<:Real}) = sum(QL_flux_out[:, 1, 2])

Qi(QL_flux_out::Array{<:Real}) = sum(QL_flux_out[:, 2:end, 2])

Πi(QL_flux_out::Array{<:Real}) = sum(QL_flux_out[:, 2:end, 3])

Γe(QL_flux_out::Array{<:Real}) = sum(QL_flux_out[:, 1, 1])

Γi(QL_flux_out::Array{<:Real}) = [sum(QL_flux_out[:, k, 1]) for k in 2:size(QL_flux_out)[2]]
