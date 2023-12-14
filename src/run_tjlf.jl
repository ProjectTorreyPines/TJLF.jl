function run(inputTJLF::InputTJLF)
    checkInput(inputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    QL_weights, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
    # for i in 1:size(eigenvalue,1)
    #     eigenvalue[i,1,1] = min(eigenvalue[i,1,1]/inputTJLF.KY_SPECTRUM[1], eigenvalue[i,2,1]/inputTJLF.KY_SPECTRUM[2]*.99)*inputTJLF.KY_SPECTRUM[1]
    # end -ORSO hack
    QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weights)
    return QL_weights, eigenvalue, QL_flux_out, flux_spectrum
end
function run(inputTGLF::InputTGLF)
    inputTJLF = InputTJLF{Float64}(inputTGLF)
    checkInput(inputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    QL_weights, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
    QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weights)
    return QL_weights, eigenvalue, QL_flux_out, flux_spectrum
end


"""
    function run_tjlf(input_tjlfs::Vector{InputTJLF})

parameters:
    input_tjlfs::InputTJLF                  - InputTJLF struct
    
description:
    Runs TJLF on a single InputTJLF struct, during the run, will save the width spectrum and eigenvalue spectrum to the InputTJLF struct.
    If you want to use these widths and eigenvalues in future runs, there is a flag: FIND_WIDTH and FIND_EIGEN that you set to false
"""
function run_tjlf(inputTJLF::InputTJLF)
    checkInput(inputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    QL_weights, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
    return  sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weights)[1]
end
"""
    function run_tjlf(input_tjlfs::Vector{InputTJLF})

parameters:
    input_tjlfs::Vector{InputTJLF}          - vector of InputTJLF structs
    
outputs:
    outputs                                 - vector of the fluxes
    
description:
    Runs TJLF on a vector of InputTJLF structs, during the run, will save the width spectrum and eigenvalue spectrum to the InputTJLF struct.
    If you want to use these widths and eigenvalues in future runs, there is a flag: FIND_WIDTH and FIND_EIGEN that you set to false
"""
function run_tjlf(input_tjlfs::Vector{InputTJLF})
    checkInput(input_tjlfs)
    outputs = Vector{Array{Float64, 3}}(undef,length(input_tjlfs))
    Threads.@threads for idx in eachindex(input_tjlfs)
        outputs[idx] = TJLF.run_tjlf(input_tjlfs[idx])
    end

    return outputs
end