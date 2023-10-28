# calls the fortran code as a shared library
# ccall((:main, "./src/Fortran/tglf.so"), Cvoid, () ,)
include("../src/TJLF.jl")
using .TJLF

#******************************************************************************************************
# Read input.tglf
#******************************************************************************************************
# location for the input.tglf file
baseDirectory = "../outputs/test_TM/simple_test/"
# baseDirectory = "../outputs/test_TM/exb_shear0/"

inputTJLF = readInput(baseDirectory)

#*******************************************************************************************************
#   start running stuff
#*******************************************************************************************************

outputHermite = gauss_hermite(inputTJLF)
satParams = get_sat_params(inputTJLF)
ky_spect, nky = get_ky_spectrum(inputTJLF, satParams.grad_r0)
fluxes, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite, ky_spect)

sum_ky_spectrum(inputTJLF, satParams, ky_spect, eigenvalue[1,:,:], fluxes)


@profview tjlf_TM(inputTJLF, satParams, outputHermite, ky_spect)
@profview_allocs tjlf_TM(inputTJLF, satParams, outputHermite, ky_spect)
# using Profile
# @profile tjlf_TM(inputTJLF, satParams, outputHermite, ky_spect)
# @codewarning
using BenchmarkTools
@btime tjlf_TM(inputTJLF, satParams, outputHermite, ky_spect)


#*******************************************************************************************************
#   plot results
#*******************************************************************************************************

using Plots

# calls the fortran code as an executable
# path = "./Fortran/tglf"
# run(`$(path) $baseDirectory`)

fileDirectory = baseDirectory * "out.tglf.QL_flux_spectrum"
lines = readlines(fileDirectory)
(ntype, nspecies, nfield, nky, nmodes) = parse.(Int32, split(lines[4]))
ql = Vector{Float64}()
for line in lines[7:length(lines)]
    line = split(line)
    if any(occursin.(["m","s"],string(line))) continue end

    for x in line
        push!(ql,parse(Float64, string(x)))
    end

end
QLw = reshape(ql, (ntype, nky, nmodes, nfield, nspecies))
QL_data = permutedims(QLw,(2,3,5,4,1))
particle_QL = QL_data[:, :, :, :, 1]
energy_QL = QL_data[:, :, :, :, 2]
toroidal_stress_QL = QL_data[:, :, :, :, 3]
parallel_stress_QL = QL_data[:, :, :, :, 4]
exchange_QL = QL_data[:, :, :, :, 5]

plot(ky_spect, particle_QL[:,1,1,1], label="Fortran")
plot!(ky_spect, fluxes[1,1,1,:,1], label="Julia", title="particle flux")

plot(ky_spect, energy_QL[:,1,1,1], label="Fortran")
plot!(ky_spect, fluxes[1,1,1,:,2], label="Julia", title="energy flux")

plot(ky_spect, exchange_QL[:,1,1,1], label="Fortran")
plot!(ky_spect, fluxes[1,1,1,:,5], label="Julia", title="exchange flux")

plot(ky_spect, toroidal_stress_QL[:,1,1,1], label="Fortran", title="toroidal stress")
plot!(ky_spect, fluxes[1,1,1,:,3], label="Julia", title="toroidal stress",linestyle=:dash)

plot(ky_spect, parallel_stress_QL[:,1,1,1], label="Fortran", title="parallel stress")
plot!(ky_spect, fluxes[1,1,1,:,4], label="Julia", title="parallel stress",linestyle=:dash)


# Get eigenvalue spectrum
nmodes = inputTJLF.NMODES
fileDirectory = baseDirectory * "out.tglf.eigenvalue_spectrum"
lines = readlines(fileDirectory)
lines = split(join(lines[3:length(lines)]))
lines = [parse(Float64, l) for l in lines]

gamma = []
freq = []
for k in 1:nmodes
    push!(gamma, lines[2k-1:2*nmodes:end])
    push!(freq, lines[2k:2*nmodes:end])
end

gammaJulia = eigenvalue[1,:,1]
freqJulia = eigenvalue[2,:,1]

plot(ky_spect, freq, label="Fortran")
plot!(ky_spect, freqJulia, label="Julia", title="Frequency",linestyle=:dash)
plot(ky_spect, gamma, label="Fortran")
plot!(ky_spect, gammaJulia, label="Julia", title="Growth Rate",linestyle=:dash)