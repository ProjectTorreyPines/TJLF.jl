using Test
using Revise
# using AxisArrays

include("../src/tjlf_modules.jl")
include("../src/tjlf_multiscale_spectrum.jl")
include("../src/tjlf_geometry.jl")

baseDirectory = "./outputs/test_SAT3/"





fileDirectory = baseDirectory * "out.tglf.QL_flux_spectrum"
lines = readlines(fileDirectory)
# 5, 2, 3, 21, 2
(ntype, nspecies, nfield, nky, nmodes) = parse.(Int32, split(lines[4]))
ql = Vector{Float64}()
for line in lines[7:length(lines)]
    line = split(line)
    if any(occursin.(["m","s"],string(line))) continue end

    for x in line
        push!(ql,parse(Float64, string(x)))
    end

end
# 5, 21, 2, 3, 2
QLw = reshape(ql, (ntype, nky, nmodes, nfield, nspecies))
QL_data = permutedims(QLw,(2,3,5,4,1))
# QL = AxisArray(QL_data, Axis{:nky}(1:nky), Axis{:nmodes}(1:nmodes), Axis{:nspecies}(1:nspecies), Axis{:nfield}(1:nfield), Axis{:ntype}(("particle","energy","toroidal_stress","parrllel_stress","exchange_stress")))

particle_QL = QL_data[:, :, :, :, 1]
energy_QL = QL_data[:, :, :, :, 2]
toroidal_stress_QL = QL_data[:, :, :, :, 3]
parallel_stress_QL = QL_data[:, :, :, :, 4]
exchange_QL = QL_data[:, :, :, :, 5]

# Read spectral shift and ave_p0 (only needed for SAT0)
fileDirectory = baseDirectory * "out.tglf.spectral_shift_spectrum"
lines = readlines(fileDirectory)
kx0_e = Vector{Float64}()
for line in lines[6:length(lines)]
        push!(kx0_e,parse(Float64, line))
end

fileDirectory = baseDirectory * "out.tglf.ave_p0_spectrum"
lines = readlines(fileDirectory)
ave_p0 = Vector{Float64}()
for line in lines[4:length(lines)]
        push!(ave_p0,parse(Float64, line))
end

# Read scalar saturation parameters
fileDirectory = baseDirectory * "out.tglf.scalar_saturation_parameters"
lines = readlines(fileDirectory)
inputComparison = Dict()
for line in lines[1:length(lines)]
    line = split(line, "\n")
    #### no idea why this is here
    if any(occursin.(["!", "UNITS", "SAT_RULE", "XNU_MODEL", "ETG_FACTOR", "R_unit", "ALPHA_ZF", "RULE"],string(line)))
        if contains(line[1],"R_unit")
            line = split(line[1]," = ")
            global R_unit = parse(Float64, strip(line[2]))
        end
        continue
    end
        
    line = split(line[1]," = ")
    line .= strip.(line)
    inputComparison[string(line[1])] = parse(Float64, line[2])
end

# Get ky spectrum
fileDirectory = baseDirectory * "out.tglf.ky_spectrum"
lines = readlines(fileDirectory)
ky_spect = Array{Float64}(undef, 0)
for line in lines[3:length(lines)]
        push!(ky_spect,parse(Float64, line))
end

# Get eigenvalue spectrum
fileDirectory = baseDirectory * "out.tglf.eigenvalue_spectrum"
lines = readlines(fileDirectory)
lines = split(join(lines[3:length(lines)]))
lines = [parse(Float64, l) for l in lines]

gamma = []
freq = []
for k in 1:nmodes
    push!(gamma, lines[2k-1:2*nmodes:end])
    push!(freq, lines[2k-1:2*nmodes:end])
end
# gamma = xr.DataArray(gamma, dims=('mode_num', 'ky'), coords={'ky': ky_spect, 'mode_num': np.arange(nmodes) + 1})
# freq = xr.DataArray(freq, dims=('mode_num', 'ky'), coords={'ky': ky_spect, 'mode_num': np.arange(nmodes) + 1})


gammas = hcat(gamma...)
R_unit = ones(size(gammas)) * R_unit

# Get potential spectrum
fileDirectory = baseDirectory * "out.tglf.field_spectrum"
lines = readlines(fileDirectory)

columns = split.(lines[2],",")
nc = length(columns)

lines = split(join(lines[6:length(lines)]))

tmpdict = Dict()
for (ik, k) in enumerate(columns)
    tmp = []
    for nm in 1:nmodes
        push!(tmp, parse.(Float64,lines[ik - 3 + nm * nc:nc*nmodes:end]))
    end
    tmpdict[k] = tmp
end
potentialTmp = tmpdict[columns[length(columns)]]
potential = hcat(potentialTmp...)


fileDirectory = baseDirectory * "out.tglf.gbflux"
lines = readlines(fileDirectory)
width::Integer = round(length(split(lines[1]))/4)
fluxes = transpose(reshape(parse.(Float64,split(lines[1])), (width,4)))




# Read input.tglf
fileDirectory = baseDirectory * "input.tglf"
lines = readlines(fileDirectory)
# struct attempt
inputTJLF = InputTJLF{Float64}()

for line in lines[1:length(lines)]
    line = split(line, "\n")
    line = split(line[1],"=")
    if line[1] == "NS"
        global inputSpecies = Vector{Species{Float64}}(undef, parse(Int, strip(line[2])))
        for i in 1:length(inputSpecies)
            inputSpecies[i] = Species{Float64}()
        end
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
        if parse(Int,speciesIndex) > length(inputSpecies) continue end
        setfield!(inputSpecies[parse(Int,speciesIndex)],    speciesField,     parse(Float64,strip(line[2], ['\'','.',' '])))
    
        # if not for the species vector
    else
        field = Symbol(line[1])
        # string
        if line[2][1] == '\''
            val = string(strip(line[2], ['\'']))
        # bool
        elseif line[2][1] == '.'
            val = strip(line[2], ['\'']) == "true"
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
setfield!(inputTJLF,:SPECIES,inputSpecies)

inputTJLF.UNITS = "CGYRO"

kx0epy, satgeo0, satgeo1, satgeo2, runit, bt0, bgeo0, gradr0, _, _, _, _ = get_sat_params(inputTJLF, ky_spect, Matrix(gammas))
@assert isapprox(kx0epy, kx0_e, rtol=1e-3)
@assert isapprox(inputComparison["SAT_geo0_out"], satgeo0, rtol=1e-6)
@assert isapprox(inputComparison["SAT_geo1_out"], satgeo1, rtol=1e-6)
@assert isapprox(inputComparison["SAT_geo2_out"], satgeo2, rtol=1e-6)
@assert isapprox(R_unit[1, 1], runit,  rtol=1e-6)
@assert isapprox(inputComparison["Bt0_out"], bt0, rtol=1e-6)
@assert isapprox(inputComparison["grad_r0_out"], gradr0, rtol=1e-6)

if inputTJLF.VEXB_SHEAR != 0.0
    @assert isapprox(inputComparison["B_geo0_out"], bgeo0, rtol=1e-6)
end


sat_1 = sum_ky_spectrum(inputTJLF, ky_spect, gammas, ave_p0, R_unit, kx0_e, potential, particle_QL, energy_QL, toroidal_stress_QL, parallel_stress_QL, exchange_QL)
julia_sat1 = sum(sum(sat_1["energy_flux_integral"], dims=3)[:,:,1], dims=1)[1,:]
expected_sat1 = fluxes[2,:]
@assert isapprox(julia_sat1, expected_sat1, rtol=1e-3)


println("SUCCESS")