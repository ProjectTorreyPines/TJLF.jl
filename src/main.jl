# calls the fortran code as a shared library
# ccall((:main, "./src/Fortran/tglf.so"), Cvoid, () ,)
include("../src/TJLF.jl")
using .TJLF

#******************************************************************************************************
# Read input.tglf
#******************************************************************************************************
# location for the input.tglf file
baseDirectory = "../outputs/test_TM/nmodes2/"

inputTJLF = readInput(baseDirectory)
inputTJLF2 = readInput(baseDirectory)

#*******************************************************************************************************
#   start running stuff
#*******************************************************************************************************

outputHermite = gauss_hermite(inputTJLF)
satParams = get_sat_params(inputTJLF)
inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
fluxes, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
QL_flux_out = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], fluxes)

inputTJLF2.KY_SPECTRUM .= inputTJLF.KY_SPECTRUM
inputTJLF2.WIDTH_SPECTRUM .= inputTJLF.WIDTH_SPECTRUM
inputTJLF2.FIND_WIDTH = false

fluxes2, eigenvalue2 = tjlf_TM(inputTJLF2, satParams, outputHermite)

for i in eachindex(fluxes)
    if(!isapprox(fluxes[i],fluxes2[i],rtol=1e-6))
        println(i)
    end
end

#*******************************************************************************************************
#   plot energy flux with varied RLTS1
#*******************************************************************************************************

using Plots

inputTJLF2.WIDTH_SPECTRUM .= inputTJLF2.WIDTH 
inputTJLF2.WIDTH_SPECTRUM .= inputTJLF.WIDTH_SPECTRUM

rltsGrid = []
electronEnergy = []
ionEnergy = []

for rlts in 2:0.2:4
    inputTJLF2.RLTS[1] = rlts
    fluxes2, eigenvalue2 = tjlf_TM(inputTJLF2, satParams, outputHermite)
    QL_flux_out = sum_ky_spectrum(inputTJLF2, satParams, eigenvalue2[:,:,1], fluxes2)

    push!(rltsGrid,rlts)
    push!(electronEnergy,sum(QL_flux_out[1,1,2]))
    push!(ionEnergy,sum(QL_flux_out[1,2,2]))
end
plot(rltsGrid, electronEnergy, title="Q_e vs RLTS_1; varied width", xlabel="RLTS_1", ylabel="Q_e", label = nothing)
plot(rltsGrid, ionEnergy, title="Q_i vs RLTS_1; varied width", xlabel="RLTS_1", ylabel="Q_i", label = nothing)


plot(inputTJLF2.KY_SPECTRUM, inputTJLF2.WIDTH_SPECTRUM, title="widths", label = nothing)

#*******************************************************************************************************
#   profiling
#*******************************************************************************************************

@profview tjlf_TM(inputTJLF2, satParams, outputHermite)
# @profview_allocs tjlf_TM(inputTJLF, satParams, outputHermite)
# using Profile
# @profile tjlf_TM(inputTJLF, satParams, outputHermite)
# @codewarning
using BenchmarkTools
@btime tjlf_TM(inputTJLF2, satParams, outputHermite)


#*******************************************************************************************************
#   plot results
#*******************************************************************************************************

using Plots

# calls the fortran code as an executable
path = "./Fortran/tglf"
run(`$(path) $baseDirectory`)

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
QL_data = permutedims(QLw,(4,5,3,2,1)) # (nf,ns,nm,nky,ntype)

# (nf,ns,nm,nky,ntype)
plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,1], label="Fortran")
plot!(inputTJLF.KY_SPECTRUM, fluxes[1,1,1,:,1], label="Julia", title="particle flux")

plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,2], label="Fortran")
plot!(inputTJLF.KY_SPECTRUM, fluxes[1,1,1,:,2], label="Julia", title="energy flux")

plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,5], label="Fortran")
plot!(inputTJLF.KY_SPECTRUM, fluxes[1,1,1,:,5], label="Julia", title="exchange flux")

plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,3], label="Fortran", title="toroidal stress")
plot!(inputTJLF.KY_SPECTRUM, fluxes[1,1,1,:,3], label="Julia", title="toroidal stress",linestyle=:dash)

plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,4], label="Fortran", title="parallel stress")
plot!(inputTJLF.KY_SPECTRUM, fluxes[1,1,1,:,4], label="Julia", title="parallel stress",linestyle=:dash)


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
freqJulia = eigenvalue[1,:,2]

plot(inputTJLF.KY_SPECTRUM, freq, label="Fortran")
plot!(inputTJLF.KY_SPECTRUM, freqJulia, label="Julia", title="Frequency",linestyle=:dash)
plot(inputTJLF.KY_SPECTRUM, gamma, label="Fortran")
plot!(inputTJLF.KY_SPECTRUM, gammaJulia, label="Julia", title="Growth Rate",linestyle=:dash)