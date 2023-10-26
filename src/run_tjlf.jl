function run(inputTJLF::InputTJLF)
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    ky_spect, nky = get_ky_spectrum(inputTJLF, satParams.grad_r0)
    fluxes, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite, ky_spect)
    return fluxes
end