module test
include("../src/tjlf_geometry.jl")
using .geometry

baseDirectory = "../outputs/test2"
fileDirectory = baseDirectory * "out.tglf.QL_flux_spectrum"
lines = readlines(fileDirectory)

# 5, 2, 3, 21, 2
(ntype, nspecies, nfield, nky, nmodes) = parse.(Int32, split(lines[4]))

ql = []
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

# QL_data = QL_data_array.transpose('ky', 'mode', 'species', 'field', 'type').data
particle_QL = QL_data[:, :, :, :, 1]
energy_QL = QL_data[:, :, :, :, 2]
toroidal_stress_QL = QL_data[:, :, :, :, 3]
parallel_stress_QL = QL_data[:, :, :, :, 4]
exchange_QL = QL_data[:, :, :, :, 5]

# Read spectral shift and ave_p0 (only needed for SAT0)
fileDirectory = baseDirectory * "out.tglf.spectral_shift_spectrum"
lines = readlines(fileDirectory)
kx0_e = []
for line in lines[6:length(lines)]
        push!(kx0_e,parse(Float64, line))
end

fileDirectory = baseDirectory * "out.tglf.ave_p0_spectrum"
lines = readlines(fileDirectory)
ave_p0 = []
for line in lines[4:length(lines)]
        push!(ave_p0,parse(Float64, line))
end

# Read scalar saturation parameters
fileDirectory = baseDirectory * "out.tglf.scalar_saturation_parameters"
lines = readlines(fileDirectory)
inputs = Dict()
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
    inputs[string(line[1])] = parse(Float64, line[2])
end

# Read input.tglf
fileDirectory = baseDirectory * "input.tglf"
lines = readlines(fileDirectory)
for line in lines[2:length(lines)]
    line = split(line, "\n")
    line = strip.(split(line[1],"="))
    try
        inputs[string(line[1])] = parse(Float64, line[2])
    catch ValueError
        continue
    end
        
end

# Added inputs
inputs["UNITS"] = "GYRO"
inputs["ALPHA_ZF"] = 1.0
inputs["RLNP_CUTOFF"] = 18.0
inputs["NS"] = ceil(inputs["NS"])
inputs["ALPHA_QUENCH"] = 0.0

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
# for k, v in list(tmpdict.items()):
#     potential = xr.DataArray(v, dims=('mode_num', 'ky'), coords={'ky': ky_spect, 'mode_num': np.arange(nmodes) + 1})
# potential = potential.T
potentialTmp = tmpdict[columns[length(columns)]]
potential = hcat(potentialTmp...)


fileDirectory = baseDirectory * "out.tglf.gbflux"
lines = readlines(fileDirectory)
width::Integer = length(lines)/2
fluxes = transpose(reshape(parse.(Float64,split(lines[1])), (2,width)))

# sat_1 = sum_ky_spectrum(inputs['SAT_RULE'], ky_spect, gammas, ave_p0, R_unit, kx0_e, potential,
#                         particle_QL, energy_QL, toroidal_stress_QL, parallel_stress_QL, exchange_QL, **inputs)

# expected_sat1 = fluxes[1]
# python_sat1 = np.sum(np.sum(sat_1['energy_flux_integral'], axis=2), axis=0)

# assert_allclose(python_sat1, expected_sat1, rtol=1e-3)

inputs["DRMINDX_LOC"] = 1.0
inputs["ALPHA_E"] = 1.0
inputs["VEXB_SHEAR"] = 0.080234
inputs["SIGN_IT"] = 1.0
kx0epy, satgeo1, satgeo2, runit, bt0, bgeo0, gradr0, _, _, _, _ = get_sat_params(1, ky_spect, Matrix(gammas'), inputs)

@assert isapprox(kx0epy, kx0_e, rtol=1e-3)
@assert isapprox(inputs["SAT_geo1_out"], satgeo1, rtol=1e-6)
@assert isapprox(inputs["SAT_geo2_out"], satgeo2, rtol=1e-6)
@assert isapprox(R_unit[1, 1], runit,  rtol=1e-6)
@assert isapprox(inputs["Bt0_out"], bt0, rtol=1e-6)
@assert isapprox(inputs["grad_r0_out"], gradr0, rtol=1e-6)

if inputs["VEXB_SHEAR"] != 0.0
    @assert isapprox(inputs["B_geo0_out"], bgeo0, rtol=1e-6)
end

end