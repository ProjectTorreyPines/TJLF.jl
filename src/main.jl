# calls the fortran code as an executable
##path = "/Users/benagnew/TJLF.jl/src/Fortran/tglf"
##baseDirectory = "-e /Users/benagnew/TJLF.jl/outputs/testB/tglf01"
##run(`$(path) $(baseDirectory)`) 


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

baseDirectory = "/Users/benagnew/TJLF.jl/outputs/testB/mxh_tglf/" 
inputTJLF = readInput(baseDirectory*"input.tglf")
baseDirectory = "/Users/benagnew/TJLF.jl/outputs/testB/tglf_diff/"
inputTJLF2 = readInput(baseDirectory*"input.tglf") # Creates the inputTJLF struct from the input.tglf file
inputTJLF2 = inputTJLF
inputTJLF.USE_TRANSPORT_MODEL
inputTJLF.WIDTH
inputTJLF.SMALL
inputTJLF.ALPHA_ZF #= 1.0=#
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
inputTJLF.USE_BPER #= true=#
inputTJLF.VPAR_1
inputTJLF.VPAR_2
inputTJLF.VPAR_3
inputTJLF.VPAR_SHEAR_1
inputTJLF.VPAR_SHEAR_2
inputTJLF.VPAR_SHEAR_3
inputTJLF.ZS_1
inputTJLF.ZS_2
inputTJLF.ZS_3
inputTJLF.SMALL #= 1.0e-13 =#
tDirectory = "/Users/benagnew/TJLF.jl/outputs/testB/mxh_tglf/"
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

# eigenvalue
# length(inputTJLF.KY_SPECTRUM)

QL_flux_out, flux_out = sum_ky_spectrum(inputTJLF, satParams, eigenvalue[:,:,1], QL_weight)
final_flux = sum(QL_flux_out,dims=1)[1,:,:]

inputTJLF2.KY_SPECTRUM .= inputTJLF.KY_SPECTRUM
inputTJLF2.WIDTH_SPECTRUM .= inputTJLF.WIDTH_SPECTRUM
inputTJLF2.FIND_WIDTH = false
#println(eigenvalue[1,:,1])
begin
    plot(inputTJLF.KY_SPECTRUM, eigenvalue[1,:,1]; label="mxh")
    display(plot!(inputTJLF.KY_SPECTRUM, [0.048818761241213, 0.0, 0.2305094014790693, 0.019255098480247684, 0.0, 0.0708815920453197, 0.07220836901415158, 0.02321552087525964, 0.029154359286357928, 0.006948842338220094, 0.07837577272431888, 0.03892311389470833, 0.14061432243514374, 0.23279705433084796, 0.44403088011383474, 0.6195920787404735, 0.7898350680943304, 2.081575910068458, 2.6331660463152904, 3.388316243354535, 3.024347207799407]; label="0s", title="gamma"))
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
@btime tjlf_TM(inputTJLF, satParams, outputHermite)

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

# Plotting a test MXH case for geometry:

# MXH:
ts = [0.0, -0.044219048771875484, -0.08843850764747149, -0.13265933572566832, -0.1768835730404605, -0.22111483984976155, -0.2653587904324358, -0.3096235122300335, -0.3539198652751106, -0.3982617608874229, -0.4426663821506042, -0.4871543513706995, -0.5317498513529563, -0.5764807078157208, -0.6213784395703039, -0.6664782812473689, -0.7118191803438845, -0.7574437661425031, -0.803398282467054, -0.849732469017662, -0.8964993667899864, -0.9437550113428812, -0.9915579629398543, -1.039968604520359, -1.0890481172489859, -1.1388570204054267, -1.1894531412175287, -1.2408888683398043, -1.2932075533087843, -1.346438978254564, -1.400593933175975, -1.4556581703934668, -1.511586338420511, -1.5682969050600049, -1.6256694383718626, -1.683545709673742, -1.741735675733485, -1.800028390883123, -1.8582065224498905, -1.916061957392908, -1.973409600163563, -2.0300970939702094, -2.0860095095300815, -2.141069388643852, -2.19523339037043, -2.2484870281689653, -2.3008387751049537, -2.3523144116030728, -2.4029520936119955, -2.452798320731159, -2.501904798344909, -2.550326092833351, -2.5981179438669946, -2.645336096795927, -2.692035533462986, -2.7382700007755107, -2.784091757441902, -2.8295514778871422, -2.874698267693338, -2.919579757011117, -2.964242247661221, -3.0087308966004516, -3.0530899235473314, -3.0973628342630595, -3.141592653589793, -3.1858224729165268, -3.230095383632255, -3.2744544105791347, -3.3189430595183653, -3.363605550168469, -3.408487039486248, -3.453633829292444, -3.499093549737684, -3.5449153064040755, -3.5911497737166003, -3.637849210383659, -3.6850673633125917, -3.7328592143462354, -3.7812805088346773, -3.8303869864484272, -3.8802332135675908, -3.9308708955765135, -3.9823465320746325, -4.034698279010621, -4.087951916809156, -4.142115918535734, -4.197175797649505, -4.253088213209377, -4.309775707016023, -4.367123349786678, -4.424978784729696, -4.483156916296464, -4.541449631446101, -4.599639597505844, -4.657515868807724, -4.714888402119581, -4.771598968759076, -4.827527136786119, -4.882591374003612, -4.936746328925022, -4.989977753870802, -5.042296438839782, -5.093732165962058, -5.14432828677416, -5.194137189930601, -5.243216702659227, -5.291627344239732, -5.339430295836705, -5.3866859403896, -5.433452838161925, -5.479787024712532, -5.525741541037084, -5.571366126835702, -5.616707025932217, -5.6618068676092825, -5.706704599363865, -5.75143545582663, -5.796030955808886, -5.840518925028982, -5.884923546292163, -5.929265441904476, -5.9735617949495525, -6.01782651674715, -6.062070467329825, -6.106301734139126, -6.150525971453918, -6.194746799532115, -6.238966258407711, -6.283185307179586]
Rs = [3.4727160921602036, 3.4731343730035356, 3.472138199764138, 3.4697347138408614, 3.465937403562048, 3.4607655078449566, 3.454243227662104, 3.446398792758445, 3.4372634368268273, 3.4268703355862016, 3.4152535594436313, 3.4024470866028307, 3.388483914771467, 3.37339530120203, 3.3572101527335, 3.3399545806218995, 3.3216516298474863, 3.3023211896235107, 3.2819800911766475, 3.2606424005531194, 3.238319918152206, 3.2150229027164094, 3.1907610452595208, 3.16554472720662, 3.139386605551134, 3.11230357370762, 3.0843191457781423, 3.055466297352693, 3.025790757585095, 2.9953546718159823, 2.964240427600195, 2.932554252667275, 2.9004289653668027, 2.8680250417442354, 2.8355290685730576, 2.8031488265771975, 2.7711048053611047, 2.739618847313436, 2.7089015678210346, 2.679140758889292, 2.650492814299829, 2.6230783500224577, 2.596982044452761, 2.572255794597244, 2.5489238633548217, 2.526988755219752, 2.5064368950257796, 2.487243580287351, 2.469377001488903, 2.4528013351833167, 2.437479024346714, 2.4233723999848604, 2.4104447970870875, 2.3986612980355804, 2.387989210487063, 2.3783983613448982, 2.369861266655029, 2.362353219878841, 2.355852327782143, 2.3503395134931346, 2.3457984994033794, 2.3422157778558343, 2.339580574435246, 2.337884806715703, 2.3371230401949417, 2.3372924470317438, 2.3383927813105627, 2.340426315926369, 2.3433977902093934, 2.3473143693434526, 2.3521856084708745, 2.3580234226868546, 2.3648420634019494, 2.3726581001345224, 2.381490404381117, 2.3913601284009753, 2.402290666022835, 2.4143075742929665, 2.427438423186775, 2.4417125249523925, 2.457160474436678, 2.473813407209228, 2.49170185539501, 2.510854057049466, 2.5312935643006584, 2.553036016732805, 2.5760850263817923, 2.600427289491046, 2.6260273152203704, 2.6528225165299286, 2.680719738664939, 2.709594418986905, 2.739293283215713, 2.7696407362607838, 2.800448117563773, 2.8315241876663384, 2.862684970727389, 2.893761466642543, 2.9246045086758463, 2.9550868047934338, 2.985102706137461, 3.014566436046151, 3.043409468386103, 3.0715775844312363, 3.0990279568489876, 3.125726457533016, 3.151645278539876, 3.176760888889569, 3.201052314302989, 3.224499711600302, 3.2470832062024435, 3.2687819642629448, 3.2895734768097578, 3.3094330397581606, 3.328333419541963, 3.3462446987045262, 3.3631342986825663, 3.378967177948448, 3.3937062024985956, 3.4073126823406197, 3.4197470622355253, 3.43096974778227, 3.4409420395178354, 3.449627138835527, 3.4569911812234158, 3.463004245773494, 3.4676412862970643, 3.4708829297031936, 3.4727160921602036]
Zs = [0.0, -0.03227327091250318, -0.06449058396410673, -0.09659636591624565, -0.12853578779384095, -0.16025507634030345, -0.19170175741753354, -0.22282481592865952, -0.2535747618179175, -0.28390359662528514, -0.31376467943394476, -0.343112494480476, -0.3719023250062661, -0.400089839087936, -0.42763059331185665, -0.45447945947910723, -0.4805899783465567, -0.5059136430922742, -0.5303991141689536, -0.553991367000218, -0.5766307752524229, -0.5982521360798302, -0.618783651029092, -0.6381458888795799, -0.6562507767870748, -0.6730006963478824, -0.6882878043505899, -0.7019937558363182, -0.7139900783542609, -0.7241395227800669, -0.7322987770932069, -0.7383229360966618, -0.7420720135766691, -0.7434194982320302, -0.7422624558948383, -0.7385320259800202, -0.7322025540414645, -0.723297369576076, -0.7118896419167491, -0.6980978305559539, -0.6820766032493347, -0.6640051368684801, -0.6440750292176369, -0.6224796410981261, -0.5994058964883371, -0.5750287880351049, -0.5495082992077707, -0.5229882047554762, -0.4955961753166919, -0.46744469223370494, -0.4386323987577114, -0.40924563025866895, -0.3793599604784938, -0.34904166979175255, -0.3183490881007197, -0.2873337946416722, -0.25604167465280464, -0.2245138425854917, -0.1927874462661857, -0.1608963681461474, -0.12887183982607223, -0.09674298524274379, -0.06453730676171471, -0.03228112723245714, -9.104292558354331e-17, 0.032281127232456956, 0.06453730676171451, 0.09674298524274363, 0.12887183982607206, 0.16089636814614722, 0.1927874462661855, 0.22451384258549154, 0.2560416746528045, 0.287333794641672, 0.31834908810071955, 0.3490416697917524, 0.37935996047849363, 0.4092456302586688, 0.4386323987577112, 0.4674446922337048, 0.49559617531669176, 0.5229882047554761, 0.5495082992077707, 0.5750287880351052, 0.5994058964883372, 0.6224796410981259, 0.6440750292176369, 0.6640051368684801, 0.6820766032493347, 0.6980978305559538, 0.7118896419167491, 0.7232973695760762, 0.7322025540414643, 0.7385320259800202, 0.7422624558948383, 0.7434194982320302, 0.7420720135766691, 0.7383229360966619, 0.7322987770932068, 0.7241395227800672, 0.713990078354261, 0.7019937558363182, 0.68828780435059, 0.6730006963478824, 0.6562507767870748, 0.6381458888795801, 0.6187836510290922, 0.5982521360798306, 0.5766307752524227, 0.5539913670002179, 0.5303991141689535, 0.505913643092274, 0.48058997834655665, 0.4544794594791078, 0.4276305933118565, 0.4000898390879359, 0.3719023250062667, 0.3431124944804763, 0.31376467943394515, 0.2839035966252857, 0.2535747618179174, 0.22282481592866013, 0.19170175741753379, 0.1602550763403034, 0.12853578779384103, 0.09659636591624551, 0.06449058396410723, 0.03227327091250338, 1.8208585116708663e-16]

# Standard:
RsM = [3.473164, 3.4724692476602654, 3.4703875004760714, 3.466926216426579, 3.462097591220862, 3.45591822106341, 3.448408657368658, 3.4395928771004702, 3.4294976952166696, 3.418152146396027, 3.405586862080865, 3.39183346634186, 3.376924010723733, 3.360890464636369, 3.343764274538818, 3.3255760025575047, 3.306355053619644, 3.286129499914885, 3.26492601270262, 3.242769914293686, 3.21968536754987, 3.195695726500182, 3.1708240795760076, 3.14509402614083, 3.1185307364940327, 3.091162353400151, 3.06302179576902, 3.0341490161437576, 3.0045937334871082, 2.9744185982307547, 2.943702633101942, 2.9125446215353574, 2.8810658948990726, 2.849411747843275, 2.817750588917168, 2.7862700539318443, 2.755169790567553, 2.7246514425357504, 2.6949072702008863, 2.6661094255372806, 2.6384018145462966, 2.611895729384635, 2.5866693779793097, 2.5627705535916894, 2.5402212556795285, 2.5190230933921676, 2.499162591669212, 2.4806158774833986, 2.463352524897231, 2.4473385390937277, 2.4325385690222983, 2.4189174822270503, 2.4064414403369123, 2.395078599264097, 2.384799536638088, 2.3755774870194517, 2.3673884460092482, 2.3602111884762933, 2.3540272337284462, 2.3488207810897967, 2.344578632429712, 2.3412901131655435, 2.3389469996539964, 2.3375434583151256, 2.337076, 2.3375434583151256, 2.338946999653996, 2.3412901131655435, 2.344578632429712, 2.3488207810897967, 2.3540272337284462, 2.3602111884762933, 2.3673884460092482, 2.3755774870194517, 2.384799536638088, 2.395078599264097, 2.4064414403369123, 2.4189174822270503, 2.4325385690222983, 2.4473385390937277, 2.4633525248972306, 2.4806158774833986, 2.499162591669212, 2.5190230933921676, 2.5402212556795285, 2.5627705535916894, 2.5866693779793097, 2.6118957293846354, 2.638401814546296, 2.66610942553728, 2.6949072702008863, 2.7246514425357504, 2.755169790567553, 2.786270053931844, 2.8177505889171677, 2.8494117478432743, 2.881065894899073, 2.912544621535357, 2.9437026331019425, 2.9744185982307543, 3.0045937334871082, 3.0341490161437568, 3.0630217957690196, 3.091162353400151, 3.1185307364940322, 3.14509402614083, 3.1708240795760076, 3.195695726500182, 3.21968536754987, 3.242769914293686, 3.2649260127026194, 3.286129499914885, 3.306355053619644, 3.3255760025575047, 3.3437642745388176, 3.360890464636369, 3.376924010723733, 3.3918334663418594, 3.4055868620808645, 3.418152146396027, 3.4294976952166696, 3.4395928771004707, 3.448408657368658, 3.45591822106341, 3.462097591220862, 3.466926216426579, 3.4703875004760714, 3.4724692476602654, 3.473164]
ZsM = [0.0, -0.03227327091250318, -0.06449058396410673, -0.09659636591624565, -0.12853578779384095, -0.16025507634030345, -0.19170175741753354, -0.22282481592865952, -0.2535747618179175, -0.28390359662528514, -0.31376467943394476, -0.343112494480476, -0.3719023250062661, -0.400089839087936, -0.42763059331185665, -0.45447945947910723, -0.4805899783465567, -0.5059136430922742, -0.5303991141689536, -0.553991367000218, -0.5766307752524229, -0.5982521360798302, -0.618783651029092, -0.6381458888795799, -0.6562507767870748, -0.6730006963478824, -0.6882878043505899, -0.7019937558363182, -0.7139900783542609, -0.7241395227800669, -0.7322987770932069, -0.7383229360966618, -0.7420720135766691, -0.7434194982320302, -0.7422624558948383, -0.7385320259800202, -0.7322025540414645, -0.723297369576076, -0.7118896419167491, -0.6980978305559539, -0.6820766032493347, -0.6640051368684801, -0.6440750292176369, -0.6224796410981261, -0.5994058964883371, -0.5750287880351049, -0.5495082992077707, -0.5229882047554762, -0.4955961753166919, -0.46744469223370494, -0.4386323987577114, -0.40924563025866895, -0.3793599604784938, -0.34904166979175255, -0.3183490881007197, -0.2873337946416722, -0.25604167465280464, -0.2245138425854917, -0.1927874462661857, -0.1608963681461474, -0.12887183982607223, -0.09674298524274379, -0.06453730676171471, -0.03228112723245714, -9.104292558354331e-17, 0.032281127232456956, 0.06453730676171451, 0.09674298524274363, 0.12887183982607206, 0.16089636814614722, 0.1927874462661855, 0.22451384258549154, 0.2560416746528045, 0.287333794641672, 0.31834908810071955, 0.3490416697917524, 0.37935996047849363, 0.4092456302586688, 0.4386323987577112, 0.4674446922337048, 0.49559617531669176, 0.5229882047554761, 0.5495082992077707, 0.5750287880351052, 0.5994058964883372, 0.6224796410981259, 0.6440750292176369, 0.6640051368684801, 0.6820766032493347, 0.6980978305559538, 0.7118896419167491, 0.7232973695760762, 0.7322025540414643, 0.7385320259800202, 0.7422624558948383, 0.7434194982320302, 0.7420720135766691, 0.7383229360966619, 0.7322987770932068, 0.7241395227800672, 0.713990078354261, 0.7019937558363182, 0.68828780435059, 0.6730006963478824, 0.6562507767870748, 0.6381458888795801, 0.6187836510290922, 0.5982521360798306, 0.5766307752524227, 0.5539913670002179, 0.5303991141689535, 0.505913643092274, 0.48058997834655665, 0.4544794594791078, 0.4276305933118565, 0.4000898390879359, 0.3719023250062667, 0.3431124944804763, 0.31376467943394515, 0.2839035966252857, 0.2535747618179174, 0.22282481592866013, 0.19170175741753379, 0.1602550763403034, 0.12853578779384103, 0.09659636591624551, 0.06449058396410723, 0.03227327091250338, 1.8208585116708663e-16]

plot(Rs, Zs; aspect_ratio=:equal, linestyle=:dash)
display(plot!(RsM, ZsM; linestyle=:dash))