function run(inputTJLF::InputTJLF{Float64})
    #println(inputTJLF.USE_TRANSPORT_MODEL)
    if (inputTJLF.USE_TRANSPORT_MODEL)
        checkInput(inputTJLF)
        outputHermite = gauss_hermite(inputTJLF)
        satParams = get_sat_params(inputTJLF)
        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
        QL_weights, eigenvalue, field_weight_out = tjlf_TM(inputTJLF, satParams, outputHermite)
        # for i in 1:size(eigenvalue,1)
        #     eigenvalue[i,1,1] = min(eigenvalue[i,1,1]/inputTJLF.KY_SPECTRUM[1], eigenvalue[i,2,1]/inputTJLF.KY_SPECTRUM[2]*.99)*inputTJLF.KY_SPECTRUM[1]
        # end -ORSO hack
        QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weights)
        return QL_weights, eigenvalue, QL_flux_out, flux_spectrum, field_weight_out, satParams
    else
        #println(inputTJLF.FIND_WIDTH)
        if (inputTJLF.FIND_WIDTH)
            # Do tjlf_max
            #println("max route")
            checkInput(inputTJLF)
            outputHermite = gauss_hermite(inputTJLF)
            satParams = get_sat_params(inputTJLF)
            inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0) # look more
            for i in 1:length(inputTJLF.KY_SPECTRUM)
                if (inputTJLF.KY == inputTJLF.KY_SPECTRUM[i])
                    kyIndex = i
                    break
                else
                    throw(error("Ky didn't match any of the ky_spectrum"))
                end
            end
            nmodes_out, gamma_nb_min_out, gamma_out, freq_out, particle_QL_out,
            energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out, field_weight_out = tjlf_max(inputTJLF, 
            satParams, outputHermite, inputTJLF.KY, inputTJLF.VEXB_SHEAR, 1)
            return gamma_out, freq_out, particle_QL_out, energy_QL_out, stress_par_QL_out, exchange_QL_out, field_weight_out, satParams, nmodes_out
        else
            # Do tjlf_ls
            #println("LS route")
            checkInput(inputTJLF)
            outputHermite = gauss_hermite(inputTJLF)
            satParams = get_sat_params(inputTJLF)
            inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0) # look more
            nmodes_out, gamma_out, freq_out, particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out,
            NaN, field_weight_out = tjlf_LS(inputTJLF, satParams, outputHermite, inputTJLF.KY, inputTJLF.NBASIS_MAX, inputTJLF.VEXB_SHEAR,
            1)
            return gamma_out, freq_out, particle_QL_out, energy_QL_out, stress_par_QL_out, exchange_QL_out, field_weight_out, satParams, nmodes_out
        end
    end
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
    function run(inputTJLF::InputTJLF)

parameters:
    input_tjlfs::InputTJLF                  - InputTJLF struct
    
description:
    Runs TJLF on a single InputTJLF struct, during the run, will save the width spectrum and eigenvalue spectrum to the InputTJLF struct.
    If you want to use these widths and eigenvalues in future runs, there is a flag: FIND_WIDTH and FIND_EIGEN that you set to false
"""
function run(inputTJLF::InputTJLF)
    checkInput(inputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    QL_weights, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
    return  sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weights)[1]
end

"""
    function run(input_tjlfs::Vector{InputTJLF})

parameters:
    input_tjlfs::Vector{InputTJLF}          - vector of InputTJLF structs
    
outputs:
    outputs                                 - vector of the fluxes
    
description:
    Runs TJLF on a vector of InputTJLF structs, during the run, will save the width spectrum and eigenvalue spectrum to the InputTJLF struct.
    If you want to use these widths and eigenvalues in future runs, there is a flag: FIND_WIDTH and FIND_EIGEN that you set to false
"""
function run(input_tjlfs::Vector{InputTJLF})
    checkInput(input_tjlfs)
    outputs = Vector{Array{Float64, 3}}(undef,length(input_tjlfs))
    Threads.@threads for idx in eachindex(input_tjlfs)
        outputs[idx] = TJLF.run(input_tjlfs[idx])
    end

    return outputs
end