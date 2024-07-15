using Test
using Base.Filesystem
include("../../src/TJLF.jl")
using ..TJLF
 #******************************************************************************#************************
 # Read input files
 #******************************************************************************#************************
        
input_miller = readInput("input_miller.tglf")
input_mxh = readInput("input_mxh.tglf")

satParams_miller = get_sat_params(input_miller)
satParams_mxh = get_sat_params(input_mxh)

outputHermite_miller = gauss_hermite(input_miller)
outputHermite_mxh = gauss_hermite(input_mxh)
input_miller.KY_SPECTRUM .= get_ky_spectrum(input_miller, satParams_miller.grad_r0)
input_mxh.KY_SPECTRUM .= get_ky_spectrum(input_mxh, satParams_mxh.grad_r0)
fluxes_miller, eigenvalue_miller = tjlf_TM(input_miller, satParams_miller, outputHermite_miller)
fluxes_mxh, eigenvalue_mxh = tjlf_TM(input_mxh, satParams_mxh, outputHermite_mxh)
QL_flux_out_miller, flux_out_miller = sum_ky_spectrum(input_miller, satParams_miller, eigenvalue_miller[:,:,1], fluxes_miller)
QL_flux_out_mxh, flux_out_mxh = sum_ky_spectrum(input_mxh, satParams_mxh, eigenvalue_mxh[:,:,1], fluxes_mxh)
using Plots
plot(1:21,eigenvalue_miller[2,:,1], xlabel = "k_y", ylabel = "omega", label = "miller")
plot!(1:21,eigenvalue_mxh[2,:,1], label = "mxh")



