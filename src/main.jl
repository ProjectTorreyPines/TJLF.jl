# calls the fortran code as a shared library
# ccall((:main, "./src/Fortran/tglf.so"), Cvoid, () ,)
# calls the fortran code as an executable
# path = "./Fortran/tglf"
# run(`$(path) $baseDirectory`)
using Revise
include("../src/TJLF.jl")
using .TJLF
using Plots

#******************************************************************************************************
# Read input.tglf
#******************************************************************************************************
# location for the input.tglf file
baseDirectory = "../outputs/test_TM/simple_test/"
# baseDirectory = "../outputs/test_TM/nmodes2/"

inputTJLF = readInput(baseDirectory)
inputTJLF2 = readInput(baseDirectory)

#*******************************************************************************************************
#   start running stuff
#*******************************************************************************************************
outputHermite = gauss_hermite(inputTJLF)
satParams = get_sat_params(inputTJLF)
inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
QL_weight, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
QL_flux_out, flux_out = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weight)

inputTJLF2.KY_SPECTRUM .= inputTJLF.KY_SPECTRUM
inputTJLF2.WIDTH_SPECTRUM .= inputTJLF.WIDTH_SPECTRUM
inputTJLF2.GAMMA_SPECTRUM .= eigenvalue[inputTJLF.NMODES,:,1]
inputTJLF2.FIND_WIDTH = false

QL_weight2, eigenvalue2 = tjlf_TM(inputTJLF2, satParams, outputHermite)

for i in eachindex(QL_weight)
    if(!isapprox(QL_weight[i],QL_weight2[i],rtol=1e-6))
        println("$i: $(QL_weight[i]) and $(QL_weight2[i])")
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
#   plot QL weights with varied RLTS2
#*******************************************************************************************************

rlts = []
xGrid = []
particleFlux = []
energyFlux = []
xGrid2 = []
particleFlux2 = []
energyFlux2 = []

for val in 3.0:0.1:6.0
    inputTJLF.RLTS[2] = val
    inputTJLF2.RLTS[2] = val
    push!(rlts,val)

    QL_weight, _, _, _ = TJLF.run(inputTJLF)
    push!(xGrid,inputTJLF.KY_SPECTRUM)
    # push!(particleFlux,flux[1,1,1,:,1])
    push!(energyFlux,QL_weight[1,2,1,:,2])

    QL_weight2, eigenvalue2, _, _ = TJLF.run(inputTJLF2)
    inputTJLF2.GAMMA_SPECTRUM .= eigenvalue2[:,:,1]
    push!(xGrid2,inputTJLF2.KY_SPECTRUM)
    # push!(particleFlux2,flux2[1,1,1,:,1])
    push!(energyFlux2,QL_weight2[1,2,1,:,2])

end

for i in 1:2
    plot()
    for r in eachindex(rlts)
        if i==0
            # plot!(xGrid[r],particleFlux[r],title="electron particle flux; ITG; find width", xlabel="ky", ylabel="flux", label = "RLTS_2 = $(rlts[r])")
        elseif i==1
            plot!(xGrid[r],energyFlux[r],title="ion energy flux; ITG; find width", xlabel="ky", ylabel="flux", label = "RLTS_2 = $(rlts[r])")
        elseif i==2
            plot!(xGrid2[r],energyFlux2[r],title="ion particle flux; ITG; fixed width", xlabel="ky", ylabel="flux", label = "RLTS_2 = $(rlts[r])")
        elseif i==3
            # plot!(xGrid2[r],energyFlux2[r],title="ion energy flux; ITG; fixed width", xlabel="ky", ylabel="flux", label = "RLTS_2 = $(rlts[r])")
        end
    end
    display(plot!(legendfontsize=3))
end

#*******************************************************************************************************
#   plot fluxes with varied RLTS2
#*******************************************************************************************************

rlts = []
xGrid = []
particleFlux = []
energyFlux = []


for val in 0.1:0.1:3.0
    inputTJLF.RLTS[2] = val
    inputTJLF2.RLTS[2] = val
    push!(rlts,val)

    _, _, _, flux_spect = TJLF.run(inputTJLF)
    push!(xGrid,inputTJLF.KY_SPECTRUM)
    push!(energyFlux,flux_spect[1,1,1,:,2])
end

plot()
for r in eachindex(rlts)
    plot!(xGrid[r],energyFlux[r],title="electron energy integral flux; ITG", xlabel="ky", ylabel="flux", label = "RLTS_2 = $(rlts[r])")
end
display(plot!(legendfontsize=3))
xlims!(0, .5)

#*******************************************************************************************************
#   plot integral energy fluxes with varied input variable
#*******************************************************************************************************

begin
xGrid = []
electronEnergy = []
electronEnergyZF = []
# ionEnergy = []
electronEnergy2 = []
electronEnergy2ZF = []
# ionEnergy2 = []
# plot()
@time for val in 0.1:0.5:6.1
    inputTJLF.RLTS[2] = val
    inputTJLF2.RLTS[2] = val
    push!(xGrid,val)

    # inputTJLF.ALPHA_ZF = 1.0
    # weights_out, eigenvalue, QL_flux_out, flux_out = TJLF.run(inputTJLF)
    # push!(electronEnergy,QL_flux_out[1,1,1])
    # inputTJLF.ALPHA_ZF = -1.0
    # weights_out, eigenvalue, QL_flux_out, flux_out = TJLF.run(inputTJLF)
    # push!(electronEnergyZF,QL_flux_out[1,1,1])
    # display(scatter!(inputTJLF.KY_SPECTRUM, eigenvalue[1,:,1]./inputTJLF.KY_SPECTRUM, title="gamma spectrum; ITG", xlabel="ky", ylabel="gamma/ky", label = "RLTS_2 = $val"))
    # plot!(inputTJLF.KY_SPECTRUM, weights_out[1,1,1,:,2], title="energy weight before integral", xlabel="ky", ylabel="flux", label = "RLTS_2 = $val")
    # plot!(inputTJLF.KY_SPECTRUM, intensity[:,1], title="phinorm", xlabel="ky", ylabel="phinorm", label = "RLTS_2 = $val")
    # push!(ionEnergy,QL_flux_out[1,2,2])

    inputTJLF2.ALPHA_ZF = -1.0
    weights_out2, eigenvalue2, QL_flux_out2, _ = TJLF.run(inputTJLF2)
    push!(electronEnergy2,QL_flux_out2[1,1,1])
    # inputTJLF2.ALPHA_ZF = -1.0
    # weights_out2, eigenvalue2, QL_flux_out2, _, = TJLF.run(inputTJLF2)
    # push!(electronEnergy2ZF,QL_flux_out2[1,1,1])
end
end
# vline!([0.24465894629054544], linestyle=:dash, color=:black, label = "kymin")
# plot!(xscale=:log10, yscale=:log10)
# ylims!(0.1,1)
# ylims!(0.9,2)
# xlims!(0.1, 0.3)
# ylims!(0,.5)
# display(plot!(legendfontsize=3))
plot(xGrid, electronEnergy, title="particle flux vs RLTS_2; ITG", xlabel="RLTS_2", ylabel="particle flux", label = "find widths,α_zf = 1")
plot!(xGrid, electronEnergy2, title="particle flux vs RLTS_2; ITG", xlabel="RLTS_2", ylabel="particle flux", label = "fixed widths; α_zf = 1; tol = 1e-2")
plot!(xGrid, electronEnergyZF, title="particle flux vs RLTS_2; ITG", xlabel="RLTS_2", ylabel="particle flux", label = "find widths; α_zf = -1")
plot!(xGrid, electronEnergy2ZF, title="particle flux vs RLTS_2; ITG", xlabel="RLTS_2", ylabel="particle flux", label = "fixed widths; α_zf = -1")
vline!([3], linestyle=:dash, color=:black, label = "width value")

#*******************************************************************************************************
#   compare weights/eigenvalues to Fortran
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
begin
plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,1]; label="Fortran")
display(plot!(inputTJLF.KY_SPECTRUM, QL_weight[1,1,1,:,1]; label="Julia", title="particle flux"))

plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,2]; label="Fortran")
display(plot!(inputTJLF.KY_SPECTRUM, QL_weight[1,1,1,:,2]; label="Julia", title="energy flux"))

plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,5]; label="Fortran")
display(plot!(inputTJLF.KY_SPECTRUM, QL_weight[1,1,1,:,5]; label="Julia", title="exchange flux"))

plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,3]; label="Fortran", title="toroidal stress")
display(plot!(inputTJLF.KY_SPECTRUM, QL_weight[1,1,1,:,3]; label="Julia", title="toroidal stress",linestyle=:dash))

plot(inputTJLF.KY_SPECTRUM, QL_data[1,1,1,:,4]; label="Fortran", title="parallel stress")
display(plot!(inputTJLF.KY_SPECTRUM, QL_weight[1,1,1,:,4]; label="Julia", title="parallel stress",linestyle=:dash))
end

# Get eigenvalue spectrum
nmodes = inputTJLF.NMODES
fileDirectory = baseDirectory * "out.tglf.eigenvalue_spectrum"
lines = readlines(fileDirectory)
lines = split(join(lines[3:length(lines)]))
lines = [parse(Float64, l) for l in lines]

begin
gamma = []
freq = []
for k in 1:nmodes
    push!(gamma, lines[2k-1:2*nmodes:end])
    push!(freq, lines[2k:2*nmodes:end])
end
gammaJulia = eigenvalue[:,:,1]
freqJulia = eigenvalue[:,:,2]

plot(inputTJLF.KY_SPECTRUM, freq; label="Fortran")
display(plot!(inputTJLF.KY_SPECTRUM, freqJulia[1,:]; label="Julia", title="Frequency",linestyle=:dash))
# display(plot!(inputTJLF.KY_SPECTRUM, freqJulia[2,:], label="Julia", title="Frequency",linestyle=:dash))
plot(inputTJLF.KY_SPECTRUM, gamma; label="Fortran")
display(plot!(inputTJLF.KY_SPECTRUM, gammaJulia[1,:]; label="Julia", title="Growth Rate",linestyle=:dash))
end