# calls the fortran code as a shared library
# ccall((:main, "./src/Fortran/tglf.so"), Cvoid, () ,)
# calls the fortran code as an executable
# path = "./Fortran/tglf"
# run(`$(path) $baseDirectory`)
include("../src/TJLF.jl")
using .TJLF
using Plots
using Lowess

#******************************************************************************************************
# Read input.tglf
#******************************************************************************************************
# location for the input.tglf file
# baseDirectory = "../outputs/test_TM/simple_test/"
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
#   plot fluxes with varied RLTS2
#*******************************************************************************************************

rlts = []
xGrid = []
particleFlux = []
energyFlux = []
xGrid2 = []
particleFlux2 = []
energyFlux2 = []

for val in 0.1:0.1:3.0
    inputTJLF.RLTS[2] = val
    inputTJLF2.RLTS[2] = val
    push!(rlts,val)

    flux, _, _ = TJLF.run(inputTJLF)
    push!(xGrid,inputTJLF.KY_SPECTRUM)
    push!(particleFlux,flux[1,1,1,:,1])
    push!(energyFlux,flux[1,1,1,:,2])

    flux2, _, _ = TJLF.run(inputTJLF2)
    push!(xGrid2,inputTJLF2.KY_SPECTRUM)
    push!(particleFlux2,flux2[1,1,1,:,1])
    push!(energyFlux2,flux2[1,1,1,:,2])

end
for i in 1:4
    plot()
    for r in rlts
        if i==1
            plot!(xGrid[i],particleFlux[i],title="particle flux; ITG; find width", xlabel="ky", ylabel="flux", label = "$RLTS(r)")
        elseif i==2
            plot!(xGrid[i],energyFlux[i],title="energy flux; ITG; find width", xlabel="ky", ylabel="flux", label = "$(r)")
        elseif i==3
            plot!(xGrid[i],particleFlux2[i],title="particle flux; ITG; fixed width", xlabel="ky", ylabel="flux", label = "$(r)")
        elseif i==4
            plot!(xGrid[i],energyFlux2[i],title="particle flux; ITG; fixed width", xlabel="ky", ylabel="flux", label = "$(r)")
        end
    end
    plot!(legendfont=6)
end

#*******************************************************************************************************
#   plot smoothed fluxes
#*******************************************************************************************************

plot(inputTJLF.KY_SPECTRUM, fluxes[1,2,1,:,1], title="particle flux", xlabel="ky", ylabel="flux", label = "original data")
plot!(inputTJLF.KY_SPECTRUM, lowess(inputTJLF.KY_SPECTRUM, fluxes[1,2,1,:,1]), title="particle flux", xlabel="ky", ylabel="flux", label = "smooth data")

#*******************************************************************************************************
#   plot energy fluxes with varied input variable
#*******************************************************************************************************

xGrid = []
electronEnergy = []
ionEnergy = []
electronEnergy2 = []
ionEnergy2 = []
for val in 0.1:0.1:6.0
    inputTJLF.RLTS[2] = val
    inputTJLF2.RLTS[2] = val
    push!(xGrid,val)

    _, eigenvalue, QL_flux_out = TJLF.run(inputTJLF)
    # plot!(inputTJLF.KY_SPECTRUM, eigenvalue[2,:,1], title="second mode; find width; ITG", xlabel="ky", ylabel="gamma", label = "RLTS_2 = $val")
    push!(electronEnergy,QL_flux_out[1,1,2])
    push!(ionEnergy,QL_flux_out[1,2,2])

    _, eigenvalue2, QL_flux_out2 = TJLF.run(inputTJLF2)
    # plot!(inputTJLF2.KY_SPECTRUM, eigenvalue2[2,:,1], title="second mode; fixed width", xlabel="ky", ylabel="gamma", label = "RLTS_2 = $val")
    push!(electronEnergy2,QL_flux_out2[1,1,2])
    push!(ionEnergy2,QL_flux_out2[1,2,2])

end
plot(xGrid, electronEnergy, title="Q_e vs RLTS_2", xlabel="RLTS_2", ylabel="Q_e", label = "find widths")
plot!(xGrid, electronEnergy2, title="Q_e vs RLTS_2; TEM", xlabel="RLTS_2", ylabel="Q_e", label = "fixed widths")
plot(xGrid, ionEnergy, title="Q_i vs RLTS_2", xlabel="RLTS_2", ylabel="Q_i", label = "find widths")
plot!(xGrid, ionEnergy2, title="Q_i vs RLTS_2; TEM", xlabel="RLTS_2", ylabel="Q_i", label = "fixed widths")

#*******************************************************************************************************
#   compare flux/eigenvalues to Fortran
#*******************************************************************************************************

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