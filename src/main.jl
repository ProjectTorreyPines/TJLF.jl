# calls the fortran code as a shared library
# ccall((:main, "./src/Fortran/tglf.so"), Cvoid, () ,)

# location for the input.tglf.gen file
baseDirectory = "../../simple_test/"
# calls the fortran code as an executable
path = "./Fortran/tglf"
# run(`$(path) $baseDirectory`)

using Plots


include("tjlf_modules.jl")
include("tjlf_hermite.jl")
include("tjlf_kygrid.jl")
include("tjlf_geometry.jl")
include("tjlf_TRANSPORT_MODEL.jl")
include("tjlf_LINEAR_SOLUTION.jl")

#******************************************************************************************************
# Read input.tglf
#******************************************************************************************************
fileDirectory = baseDirectory * "input.tglf"
lines = readlines(fileDirectory)
inputTJLF = InputTJLF{Float64}()

for line in lines[1:length(lines)]
    line = split(line, "\n")
    line = split(line[1],"=")
    if line[1] == "NS"
        inputTJLF.ZS = Vector{Float64}(undef, parse(Int, strip(line[2])))
        inputTJLF.AS = Vector{Float64}(undef, parse(Int, strip(line[2])))
        inputTJLF.MASS = Vector{Float64}(undef, parse(Int, strip(line[2])))
        inputTJLF.TAUS = Vector{Float64}(undef, parse(Int, strip(line[2])))
        inputTJLF.RLNS = Vector{Float64}(undef, parse(Int, strip(line[2])))
        inputTJLF.RLTS = Vector{Float64}(undef, parse(Int, strip(line[2])))
        inputTJLF.VPAR = Vector{Float64}(undef, parse(Int, strip(line[2])))
        inputTJLF.VPAR_SHEAR = Vector{Float64}(undef, parse(Int, strip(line[2])))
        inputTJLF.VNS_SHEAR = Vector{Float64}(undef, parse(Int, strip(line[2])))
        inputTJLF.VTS_SHEAR = Vector{Float64}(undef, parse(Int, strip(line[2])))
    end
end

for line in lines[1:length(lines)]
    line = split(line, "\n")
    line = split(line[1],"=")

    #### for the species vector
    check = match(r"_\d",line[1])
    if check !== nothing
        temp = split(line[1],"_")
        speciesField = Symbol(replace(line[1], r"_\d"=>""))
        speciesIndex = check.match[2:end]
        if parse(Int,speciesIndex) > length(inputTJLF.ZS) continue end
        getfield(inputTJLF, speciesField)[parse(Int,speciesIndex)] = parse(Float64,strip(line[2], ['\'','.',' ']))
        # setfield!(inputSpecies[parse(Int,speciesIndex)],    speciesField,     parse(Float64,strip(line[2], ['\'','.',' '])))
    else # if not for the species vector
        field = Symbol(line[1])
        
        # string
        if line[2][1] == '\''
            val = string(strip(line[2], ['\'']))
        # bool
        elseif line[2][1] == '.'
            val = strip(line[2], ['\'','.']) == "true"
        # int
        elseif !contains(line[2],'.')
            val = parse(Int, strip(line[2], ['\'','.',' ']))
        # float
        else
            val = parse(Float64,strip(line[2], ['\'','.',' ']))
        end

        try
            setfield!(inputTJLF,field,val)
        catch
            throw(error(field))
        end

    end 
end

if inputTJLF.SAT_RULE == 2 || inputTJLF.SAT_RULE == 3
    inputTJLF.UNITS = "CGYRO"
    ####### WTF
    inputTJLF.XNU_MODEL = 3
    inputTJLF.WDIA_TRAPPED = 1.0
end


#******************************************************************************#************************
#   start running stuff
#******************************************************************************#************************

outputHermite = gauss_hermite(inputTJLF)
satParams = get_sat_params(inputTJLF)
ky_spect, nky = get_ky_spectrum(inputTJLF, satParams.grad_r0)
fluxes, eigenvalue = tjlf_TM(inputTJLF, satParams, outputHermite, ky_spect)

# gamma = eigenvalue[1,:,:]
# freq = eigenvalue[2,:,:]



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
plot!(ky_spect, fluxes[2,1,1,:,1], label="Julia", title="energy flux")

plot(ky_spect, exchange_QL[:,1,1,1], label="Fortran")
plot!(ky_spect, fluxes[5,1,1,:,1], label="Julia", title="exchange flux")


plot(ky_spect, toroidal_stress_QL[:,1,1,1], label="Fortran toroidal", title="Stresses")
plot!(ky_spect, fluxes[3,1,1,:,1], label="Julia toroidal", title="Stresses",linestyle=:dash)
plot!(ky_spect, parallel_stress_QL[:,1,1,1], label="Fortran parallel", title="Stresses")
plot!(ky_spect, fluxes[4,1,1,:,1], label="Julia parallel", title="Stresses",linestyle=:dash)












# ################################# Needed from Fortran output files ############################################

# Get ky spectrum
# fileDirectory = baseDirectory * "out.tglf.ky_spectrum"
# lines = readlines(fileDirectory)
# ky_spect = Array{Float64}(undef, 0)
# for line in lines[3:length(lines)]
#         push!(ky_spect,parse(Float64, line))
# end
# # Get eigenvalue spectrum
# fileDirectory = baseDirectory * "out.tglf.eigenvalue_spectrum"
# lines = readlines(fileDirectory)
# lines = split(join(lines[3:length(lines)]))
# lines = [parse(Float64, l) for l in lines]

# gamma = []
# freq = []
# for k in 1:nmodes
#     push!(gamma, lines[2k-1:2*nmodes:end])
#     push!(freq, lines[2k-1:2*nmodes:end])
# end
# gammas = hcat(gamma...)

# fileDirectory = baseDirectory * "out.tglf.QL_flux_spectrum"
# lines = readlines(fileDirectory)
# (ntype, nspecies, nfield, nky, nmodes) = parse.(Int32, split(lines[4]))
# ql = Vector{Float64}()
# for line in lines[7:length(lines)]
#     line = split(line)
#     if any(occursin.(["m","s"],string(line))) continue end

#     for x in line
#         push!(ql,parse(Float64, string(x)))
#     end

# end
# QLw = reshape(ql, (ntype, nky, nmodes, nfield, nspecies))
# QL_data = permutedims(QLw,(2,3,5,4,1))
# particle_QL = QL_data[:, :, :, :, 1]
# energy_QL = QL_data[:, :, :, :, 2]
# toroidal_stress_QL = QL_data[:, :, :, :, 3]
# parallel_stress_QL = QL_data[:, :, :, :, 4]
# exchange_QL = QL_data[:, :, :, :, 5]
# #Read ave_p0 (only needed for SAT0)
# fileDirectory = baseDirectory * "out.tglf.ave_p0_spectrum"
# lines = readlines(fileDirectory)
# ave_p0 = Vector{Float64}()
# for line in lines[4:length(lines)]
#         push!(ave_p0,parse(Float64, line))
# end
# # Get potential spectrum
# fileDirectory = baseDirectory * "out.tglf.field_spectrum"
# lines = readlines(fileDirectory)

# columns = split.(lines[2],",")
# nc = length(columns)

# lines = split(join(lines[6:length(lines)]))

# tmpdict = Dict()
# for (ik, k) in enumerate(columns)
#     tmp = []
#     for nm in 1:nmodes
#         push!(tmp, parse.(Float64,lines[ik - 3 + nm * nc:nc*nmodes:end]))
#     end
#     tmpdict[k] = tmp
# end
# potentialTmp = tmpdict[columns[length(columns)]]
# potential = hcat(potentialTmp...)
