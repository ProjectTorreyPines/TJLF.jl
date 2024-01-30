# calls the fortran code as an executable
path = "/Users/benagnew/TJLF.jl/src/Fortran/tglf"
baseDirectory = "/Users/benagnew/TJLF.jl/outputs/testB/tglf01"
run(`$(path) $(baseDirectory)`) 

using Revise
include("../src/TJLF.jl")
using .TJLF
using Plots
using Base.Threads
Threads.nthreads()
using LinearAlgebra
BLAS.set_num_threads(1)


#******************************************************************************************************
# Read input.tglf
#******************************************************************************************************
# location for the input.tglf file
inputTJLFVector = Vector{Main.TJLF.InputTJLF}(undef, 7)
begin
baseDirectory = "../outputs/TIM_case/case2/"
inputTJLFVector[1] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case3/"
inputTJLFVector[2] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case4/"
inputTJLFVector[3] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case5/"
inputTJLFVector[4] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case6/"
inputTJLFVector[5] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case7/"
inputTJLFVector[6] = readInput(baseDirectory*"input.tglf")
baseDirectory = "../outputs/TIM_case/case8/"
inputTJLFVector[7] = readInput(baseDirectory*"input.tglf")
end


baseDirectory = "/Users/benagnew/TJLF.jl/outputs/testB/tglf01/"
inputTJLF = readInput(baseDirectory*"input.tglf")
baseDirectory = "/Users/benagnew/TJLF.jl/outputs/tglf_regression/tglf01/"
inputTJLF = readInput(baseDirectory*"input.tglf")
inputTJLF2 = inputTJLF



inputTJLF.ALPHA_ZF = 1.0
inputTJLF.AS_1
inputTJLF.AS_2
inputTJLF.AS_3
inputTJLF.MASS_1
inputTJLF.MASS_2
inputTJLF.MASS_3
inputTJLF.RLNS_1
inputTJLF.RLNS_2
inputTJLF.RLNS_3
inputTJLF.RLTS_1
inputTJLF.RLTS_2
inputTJLF.RLTS_3
inputTJLF.TAUS_1
inputTJLF.TAUS_2
inputTJLF.TAUS_3
inputTJLF.USE_BPER = true
inputTJLF.VPAR_1
inputTJLF.VPAR_2
inputTJLF.VPAR_3
inputTJLF.VPAR_SHEAR_1
inputTJLF.VPAR_SHEAR_2
inputTJLF.VPAR_SHEAR_3
inputTJLF.ZS_1
inputTJLF.ZS_2
inputTJLF.ZS_3
inputTJLF.SMALL = 1.0e-13
tDirectory = "/Users/benagnew/TJLF.jl/outputs/testB/tglf01/"
#CHECK TJLF_READ_INPUT.jl FOR inP status.

#*******************************************************************************************************
#   start running stuff
#*******************************************************************************************************
begin
    outputHermite = gauss_hermite(inputTJLF)
    satParams = get_sat_params(inputTJLF)
    inputTJLF.KY_SPECTRUM .= get_ky_spectrum(inputTJLF, satParams.grad_r0)
    QL_weight, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite)
end

QL_flux_out, flux_out = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weight)
final_flux = sum(QL_flux_out,dims=1)[1,:,:]

#The spectra don't necessarily need to be set identically because it is very likely that rho_i or rho_e, which are used to define ky grid points of equidistance and logarithmic spacing.
inputTJLF2.KY_SPECTRUM .= inputTJLF.KY_SPECTRUM
inputTJLF2.WIDTH_SPECTRUM .= inputTJLF.WIDTH_SPECTRUM
inputTJLF2.FIND_WIDTH = false

#Testing purposes:
begin
satParams2 = get_sat_params(inputTJLF2) #--
outputHermite2 = gauss_hermite(inputTJLF2) #--
inputTJLF2.KY_SPECTRUM .= get_ky_spectrum(inputTJLF2, satParams2.grad_r0) #--
QL_weight2, eigenvalue2 = tjlf_TM(inputTJLF2, satParams2, outputHermite2) #Prev. def. w/ (inputTJLF2, satParams, outputHermite)
end
#Why would satParams be used here? We just defined it as being derived from only inputTJLF but we apply it
#to inputTJLF2 which right now is /TIM_case2/case2. Are we defining these two to have exactly the same satParams?
#and outputHermite? -- This is because the ky spectrum is defined already for each of the cases: check out get_ky_spectrum function 

begin
    plot(inputTJLF.KY_SPECTRUM, eigenvalue[1,:,1]; label="correct")
    #display(plot!(inputTJLF2.KY_SPECTRUM, eigenvalue2[1,:,1]; label="fast", title="gamma"))
end
#This plot is showing the ky-spectrum vs. the eigenvalue growth rate gamma. There may be some scaling things to consider for gamma.

for i in eachindex(QL_weight)
    if(!isapprox(QL_weight[i],QL_weight2[i],rtol=1e-6))
        println("$i: $(QL_weight[i]) and $(QL_weight2[i])")
    end
end

#*******************************************************************************************************
#   profiling
#*******************************************************************************************************

@profview tjlf_TM(inputTJLF, satParams, outputHermite)
# @profview_allocs tjlf_TM(inputTJLF, satParams, outputHermite)
# using Profile
# @profile tjlf_TM(inputTJLF, satParams, outputHermite)
using BenchmarkTools
@btime tjlf_TM(inputTJLF2, satParams, outputHermite)

#*******************************************************************************************************
#   plot QL weights with varied RLTS2
#*******************************************************************************************************
begin
rlts = []
xGrid = []
particleFlux = []
energyFlux = []
xGrid2 = []
particleFlux2 = []
energyFlux2 = []

#This does involve KY_SPECTRUM which we would have just set to be equal to the first input.
for val in 3.0:0.1:6.0
    inputTJLF.RLTS[2] = val
    inputTJLF2.RLTS[2] = val
    push!(rlts,val)

    QL_weight, _, _, _ = TJLF.run(inputTJLF)
    push!(xGrid,inputTJLF.KY_SPECTRUM)
    # push!(particleFlux,flux[1,1,1,:,1])
    push!(energyFlux,QL_weight[1,2,1,:,2])

    QL_weight2, eigenvalue2, _, _ = TJLF.run(inputTJLF2)
    #inputTJLF2.EIGEN_SPECTRUM .= eigenvalue2[:,:,1]
    #What is the above line accomplishing?
    push!(xGrid2,inputTJLF2.KY_SPECTRUM)
    # push!(particleFlux2,flux2[1,1,1,:,1])
    push!(energyFlux2,QL_weight2[1,2,1,:,2])

end
end

#But the eigenvalue2 is reset here so what exactly was the point of the beginning needing
#eigenvalue2 to be set with inputTJLF's satParams and outputHermite?

#inputTJLF2.EIGEN_SPECTRUM and
#eigenvalue2[1,:,1] are identical except evalue2 has an imaginary part. Is the point to 
#just take the real part of the eigenvalue of TJLF2 here?

#I added an encompassing begin-end statement so as to not create any issues with xGrid length

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

begin
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

begin
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
fileDirectory = tDirectory * "out.tglf.eigenvalue_spectrum" #baseDirectory * "out.tglf.eigenvalue_spectrum"
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
gammaJulia = eigenvalue[:,:,1] #Was listed as eigenvalue2 for this and below line
freqJulia = eigenvalue[:,:,2]

plot(inputTJLF.KY_SPECTRUM, freq[1]; label="Fortran")
display(plot!(inputTJLF.KY_SPECTRUM, freqJulia[1,:]; label="Julia", title="Frequency",linestyle=:dash))
# display(plot!(inputTJLF.KY_SPECTRUM, freqJulia[2,:], label="Julia", title="Frequency",linestyle=:dash))
plot(inputTJLF.KY_SPECTRUM, gamma[1]; label="Fortran")
display(plot!(inputTJLF.KY_SPECTRUM, gammaJulia[1,:]; label="Julia", title="Growth Rate",linestyle=:dash))
end