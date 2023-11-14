function run(inputTJLF::InputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    QL_weights, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
    QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weights)
    return QL_weights, eigenvalue, QL_flux_out, flux_spectrum
end


function run(inputTGLF::InputTGLF)
    inputTJLF = InputTJLF{Float64}(inputTGLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    QL_weights, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
    QL_flux_out, flux_spectrum = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weights)
    return QL_weights, eigenvalue, QL_flux_out, flux_spectrum
end