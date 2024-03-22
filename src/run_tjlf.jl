function run(inputTJLF::InputTJLF)
    checkInput(inputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    QL_weights, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
    QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:, :, 1], QL_weights)
    return (QL_weights=QL_weights, eigenvalue=eigenvalue, QL_flux_out=QL_flux_out, flux_spectrum=flux_spectrum)
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
    run_tjlf(input_tjlfs::Vector{InputTJLF})

parameters:
    input_tjlfs::Vector{InputTJLF}          - vector of InputTJLF structs
    
outputs:
    outputs                                 - vector of fluxes (field, species, type)
    
description:
    Runs TJLF on a vector of InputTJLF structs, during the run, will save the width spectrum and eigenvalue spectrum to the InputTJLF struct.
    If you want to use these widths and eigenvalues in future runs, there is a flag: FIND_WIDTH and FIND_EIGEN that you set to false
"""
function run_tjlf(input_tjlfs::Vector{InputTJLF})
    checkInput(input_tjlfs)
    outputs = Vector{Array{Float64,3}}(undef, length(input_tjlfs))
    Threads.@threads for idx in eachindex(input_tjlfs)
        outputs[idx] = TJLF.run_tjlf(input_tjlfs[idx])
    end
    return outputs
end

# (field, species, type)
# type: (particle, energy, torodial stress, parallel stress, exchange)

Qe(QL_flux_out::Array{Float64}) = sum(QL_flux_out[:, 1, 2])

Qi(QL_flux_out::Array{Float64}) = sum(QL_flux_out[:, 2:end, 2])

Πi(QL_flux_out::Array{Float64}) = sum(QL_flux_out[:, 2:end, 3])

Γe(QL_flux_out::Array{Float64}) = sum(QL_flux_out[:, 1, 1])

Γi(QL_flux_out::Array{Float64}) = [sum(QL_flux_out[:, k, 1]) for k in 1:size(QL_flux_out)[2]-1]
