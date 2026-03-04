function run(inputTJLF::InputTJLF)
    use_tm = ismissing(inputTJLF.USE_TRANSPORT_MODEL) ? true : inputTJLF.USE_TRANSPORT_MODEL
    if use_tm
        checkInput(inputTJLF)
        outputHermite = gauss_hermite(inputTJLF)
        satParams = get_sat_params(inputTJLF)
        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
        tm = tjlf_TM(inputTJLF, satParams, outputHermite)
        QL_weights            = tm.QL_weights
        firstPass_eigenvalue  = tm.firstPass_eigenvalue
        secondPass_eigenvalue = tm.secondPass_eigenvalue
        field_weight_out      = tm.field_weight_out
        phi_bar_matrix        = tm.phi_bar_matrix

        vexb_shear_s = inputTJLF.VEXB_SHEAR * inputTJLF.SIGN_IT
        if inputTJLF.ALPHA_QUENCH == 0.0 && vexb_shear_s != 0.0
            eigenvalue_for_flux = secondPass_eigenvalue
        else
            eigenvalue_for_flux = firstPass_eigenvalue
        end

        if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3
            most_unstable_gamma_first_pass = firstPass_eigenvalue[1, :, 1]
            vzf_out_first_pass, kymax_out_first_pass, jmax_out_first_pass = get_zonal_mixing(inputTJLF, satParams, most_unstable_gamma_first_pass)
            QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue_for_flux[:, :, 1], QL_weights;
                                                         vzf_out_param=vzf_out_first_pass,
                                                         kymax_out_param=kymax_out_first_pass,
                                                         jmax_out_param=jmax_out_first_pass)
        else
            QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue_for_flux[:, :, 1], QL_weights;
                                                         phi_bar_matrix=phi_bar_matrix)
        end
        return (QL_weights=QL_weights, eigenvalue=firstPass_eigenvalue, QL_flux_out=QL_flux_out,
                flux_spectrum=flux_spectrum, field_weight_out=field_weight_out)
    else
        checkInput(inputTJLF)
        outputHermite = gauss_hermite(inputTJLF)
        satParams = get_sat_params(inputTJLF)
        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
        ns     = inputTJLF.NS
        nmodes = inputTJLF.NMODES
        nbasis = inputTJLF.NBASIS_MAX
        nmodes_out, gamma_out, freq_out,
            particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out,
            _ft, field_weight_out_3d, _phi = tjlf_LS(inputTJLF, satParams, outputHermite,
                                                      inputTJLF.KY, nbasis, inputTJLF.VEXB_SHEAR, 1)

        eigenvalue       = zeros(Float64, nmodes, 1, 2)
        QL_weights       = zeros(Float64, 3, ns, nmodes, 1, 5)
        field_weight_out = zeros(ComplexF64, 3, nbasis, nmodes, 1)

        eigenvalue[:, 1, 1]                      .= gamma_out
        eigenvalue[:, 1, 2]                      .= freq_out
        QL_weights[:, :, 1:nmodes_out, 1, 1]    .= particle_QL_out[:, :, 1:nmodes_out]
        QL_weights[:, :, 1:nmodes_out, 1, 2]    .= energy_QL_out[:, :, 1:nmodes_out]
        QL_weights[:, :, 1:nmodes_out, 1, 3]    .= stress_tor_QL_out[:, :, 1:nmodes_out]
        QL_weights[:, :, 1:nmodes_out, 1, 4]    .= stress_par_QL_out[:, :, 1:nmodes_out]
        QL_weights[:, :, 1:nmodes_out, 1, 5]    .= exchange_QL_out[:, :, 1:nmodes_out]
        field_weight_out[:, :, 1:nmodes_out, 1] .= field_weight_out_3d[:, :, 1:nmodes_out]

        QL_flux_out   = Array{Float64}(undef, 0, 0, 0)
        flux_spectrum = Array{Float64}(undef, 0, 0, 0)
        return (QL_weights=QL_weights, eigenvalue=eigenvalue, QL_flux_out=QL_flux_out,
                flux_spectrum=flux_spectrum, field_weight_out=field_weight_out)
    end
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
    outputs = Vector{Array{T,3}}(undef, length(input_tjlfs))
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
