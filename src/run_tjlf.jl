function run(inputTJLF::InputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    fluxes, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
    QL_flux_out = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], fluxes)
    return fluxes, eigenvalue, QL_flux_out
end