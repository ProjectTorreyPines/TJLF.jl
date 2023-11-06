function run(inputTJLF::InputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    if isnan(inputTJLF.KY_SPECTRUM[1])
        inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    end
    fluxes, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
    return fluxes
end