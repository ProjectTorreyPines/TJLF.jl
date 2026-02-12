
ENV["GKSwstype"] = "100"

using TJLF
using TurbulentTransport
using Plots

tmpdir = mktempdir()
filepath = joinpath(tmpdir, "input.tglf")
input_tglf_lines = """
ADIABATIC_ELEC = .false.
ALPHA_E = 1.0
ALPHA_MACH = 0.0
ALPHA_P = 1.0
ALPHA_QUENCH = 0
ALPHA_ZF = -1.0
AS_1 = 1.0
AS_2 = 0.784867
AS_3 = 0.0302081
BETAE = 0.00362972
BETA_LOC = 0.0
DAMP_PSI = 0.0
DAMP_SIG = 0.0
DEBYE = 0.0217677
DEBYE_FACTOR = 1.0
DELTA_LOC = 0.0681444
DRMAJDX_LOC = -0.189065
DRMINDX_LOC = 1.0
DZMAJDX_LOC = 0.00278328
ETG_FACTOR = 1.25
FILTER = 2.0
FIND_WIDTH = .true.
GCHAT = 1.0
GHAT = 1.0
GRADB_FACTOR = 0.0
IBRANCH = -1
IFLUX = .true.
KAPPA_LOC = 1.40438
KX0_LOC = 0.0
KY = 0.3
KYGRID_MODEL = 4
LINSKER_FACTOR = 0.0
MASS_1 = 0.000272445
MASS_2 = 1.0
MASS_3 = 6.0
NBASIS_MAX = 6
NBASIS_MIN = 2
NEW_EIKONAL = .true.
NKY = 12
NMODES = 2
NS = 3
NWIDTH = 21
NXGRID = 16
PARK = 1.0
P_PRIME_LOC = -0.00355359
Q_LOC = 2.00545
Q_PRIME_LOC = 14.7947
RLNP_CUTOFF = 18.0
RLNS_1 = 0.513787
RLNS_2 = 0.758616
RLNS_3 = -0.872531
RLTS_1 = 2.03987
RLTS_2 = 2.20153
RLTS_3 = 2.20153
RMAJ_LOC = 2.86212
RMIN_LOC = 0.573129
SAT_RULE = 3
SIGN_BT = -1
SIGN_IT = 1
S_DELTA_LOC = 0.116297
S_KAPPA_LOC = 0.125574
S_ZETA_LOC = -0.0258657
TAUS_1 = 1.0
TAUS_2 = 1.39296
TAUS_3 = 1.39296
THETA_TRAPPED = 0.7
UNITS = 'GYRO'
USE_AVE_ION_GRID = .false.
USE_BISECTION = .true.
USE_BPAR = .true.
USE_BPER = .true.
USE_INBOARD_DETRAPPED = .false.
USE_MHD_RULE = .false.
VEXB_SHEAR = 0.080234
VPAR_1 = 0.419061
VPAR_2 = 0.419061
VPAR_3 = 0.419061
VPAR_MODEL = 0
VPAR_SHEAR_1 = 0.803536
VPAR_SHEAR_2 = 0.803536
VPAR_SHEAR_3 = 0.803536
VPAR_SHEAR_MODEL = 1
WDIA_TRAPPED = 1.0
WD_ZERO = 0.1
WIDTH = 1.65
WIDTH_MIN = 0.3
XNUE = 0.0948099
XNU_FACTOR = 1.0
XNU_MODEL = 3
ZEFF = 1.90624
ZETA_LOC = -0.0148888
ZMAJ_LOC = -0.0576768
ZS_1 = -1.0
ZS_2 = 1.0
ZS_3 = 6.0
""";
open(filepath, "w") do f
    write(f, input_tglf_lines)
end;

# lets load the input.tglf file
input_tglf = TurbulentTransport.load(InputTGLF(),filepath);
# Ensure TJLF runs with the same parameters in TGLF for USE_PRESETS=.true.
TurbulentTransport.apply_presets!(input_tglf)

input_tjlf = InputTJLF{Float64}(input_tglf)

result = TJLF.run(input_tjlf)
eig = result.eigenvalue
println("eig has size", size(eig))
KY_spec = input_tjlf.KY_SPECTRUM

gamma = eig[1,:,1]
println(size(gamma))
omega = eig[1,:,2]
shape = size(eig)
println(eig)
println("eig shape, $shape")

# quit()

println(size(omega))
QL_spec = result.QL_weights

gr()

plot(KY_spec, gamma, label="gamma")
plot!(KY_spec, omega, label="omega")

println("gamma: ", gamma)
println("omega: ", omega)

savefig("outputs/tjlf_gammaOmega_vs_ky.png")

# println("KY_spec shape: ", size(KY_spec))
# println("QL_spec shape: ", size(QL_spec))

# flux_types = ["particle", "momentum", "energy"]
# species = ["electrons", "ions", "impurities"]
# plot()
# for i in 1:3
#     for j in 1:3
#         flux = QL_spec[i, j, 1, :, 1]
#         plot!(KY_spec, flux, label="$(flux_types[i]) $(species[j])")
#     end
# end
# savefig("outputs/tjlf_QL_vs_ky.png")
