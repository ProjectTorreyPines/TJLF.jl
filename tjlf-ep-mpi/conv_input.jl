using Revise
include("../src/TJLF.jl")
using .TJLF
using Plots
using Base.Threads
Threads.nthreads()
using LinearAlgebra
BLAS.set_num_threads(1)

# Why is InputTJLF not defined when I call these? What the hell?
#baseDirectory = "/Users/benagnew/TJLF.jl/outputs/TIM_case/case2/" 
#inputTJLF = readInput(baseDirectory*"input.tglf")

"""
convert_input is meant for corrected for the type-discrepancy in TJLF and TJLFEP. They are identical, but Julia cannot see that yet.

inputs: TJLFEP.InputTJLF{Float64} struct

outputs: Main.TJLF.InputTJLF{Float64} struct of the same values
"""
function convert_input(input::TJLFEP.InputTJLF{Float64}, ns::Int64, nky::Int64)
    # Extract relevant fields from InputTJLF{Float64} and construct Main.TJLF.InputTJLF{Float64}
    new_input = Main.TJLFEP.TJLF.InputTJLF{Float64}(ns, nky)
    new_input.UNITS = input.UNITS
    new_input.USE_BPER = input.USE_BPER
    new_input.USE_BPAR = input.USE_BPAR
    new_input.USE_MHD_RULE = input.USE_MHD_RULE
    new_input.USE_BISECTION = input.USE_BISECTION
    new_input.USE_INBOARD_DETRAPPED = input.USE_INBOARD_DETRAPPED
    new_input.USE_AVE_ION_GRID = input.USE_AVE_ION_GRID
    new_input.NEW_EIKONAL = input.NEW_EIKONAL
    new_input.FIND_WIDTH = input.FIND_WIDTH
    new_input.IFLUX = input.IFLUX
    new_input.ADIABATIC_ELEC = input.ADIABATIC_ELEC
    new_input.SAT_RULE = input.SAT_RULE
    new_input.NS = input.NS
    new_input.NMODES = input.NMODES
    new_input.NWIDTH = input.NWIDTH
    new_input.NBASIS_MAX = input.NBASIS_MAX
    new_input.NBASIS_MIN = input.NBASIS_MIN
    new_input.NXGRID = input.NXGRID
    new_input.NKY = input.NKY
    new_input.KYGRID_MODEL = input.KYGRID_MODEL
    new_input.XNU_MODEL = input.XNU_MODEL
    new_input.VPAR_MODEL = input.VPAR_MODEL
    new_input.IBRANCH = input.IBRANCH
    new_input.ZS = input.ZS
    new_input.MASS = input.MASS
    new_input.RLNS = input.RLNS
    new_input.RLTS = input.RLTS
    new_input.TAUS = input.TAUS
    new_input.AS = input.AS
    new_input.VPAR = input.VPAR
    new_input.VPAR_SHEAR = input.VPAR_SHEAR
    new_input.WIDTH_SPECTRUM = input.WIDTH_SPECTRUM
    new_input.KY_SPECTRUM = input.KY_SPECTRUM
    new_input.EIGEN_SPECTRUM = input.EIGEN_SPECTRUM
    new_input.FIND_EIGEN = input.FIND_EIGEN
    new_input.SIGN_BT = input.SIGN_BT
    new_input.SIGN_IT = input.SIGN_IT
    new_input.KY = input.KY
    new_input.VEXB_SHEAR = input.VEXB_SHEAR
    new_input.BETAE = input.BETAE
    new_input.XNUE = input.XNUE
    new_input.ZEFF = input.ZEFF
    new_input.DEBYE = input.DEBYE
    new_input.ALPHA_MACH = input.ALPHA_MACH
    new_input.ALPHA_E = input.ALPHA_E
    new_input.ALPHA_P = input.ALPHA_P
    new_input.ALPHA_QUENCH = input.ALPHA_QUENCH
    new_input.ALPHA_ZF = input.ALPHA_ZF
    new_input.XNU_FACTOR = input.XNU_FACTOR
    new_input.DEBYE_FACTOR = input.DEBYE_FACTOR
    new_input.ETG_FACTOR = input.ETG_FACTOR
    new_input.RLNP_CUTOFF = input.RLNP_CUTOFF
    new_input.WIDTH = input.WIDTH
    new_input.WIDTH_MIN = input.WIDTH_MIN
    new_input.RMIN_LOC = input.RMIN_LOC
    new_input.RMAJ_LOC = input.RMAJ_LOC
    new_input.ZMAJ_LOC = input.ZMAJ_LOC
    new_input.DRMINDX_LOC = input.DRMINDX_LOC
    new_input.DRMAJDX_LOC = input.DRMAJDX_LOC
    new_input.DZMAJDX_LOC = input.DZMAJDX_LOC
    new_input.Q_LOC = input.Q_LOC
    new_input.KAPPA_LOC = input.KAPPA_LOC
    new_input.S_KAPPA_LOC = input.S_KAPPA_LOC
    new_input.DELTA_LOC = input.DELTA_LOC
    new_input.S_DELTA_LOC = input.S_DELTA_LOC
    new_input.ZETA_LOC = input.ZETA_LOC
    new_input.S_ZETA_LOC = input.S_ZETA_LOC
    new_input.P_PRIME_LOC = input.P_PRIME_LOC
    new_input.Q_PRIME_LOC = input.Q_PRIME_LOC
    new_input.BETA_LOC = input.BETA_LOC
    new_input.KX0_LOC = input.KX0_LOC
    new_input.DAMP_PSI = input.DAMP_PSI
    new_input.DAMP_SIG = input.DAMP_SIG
    new_input.WDIA_TRAPPED = input.WDIA_TRAPPED
    new_input.PARK = input.PARK
    new_input.GHAT = input.GHAT
    new_input.GCHAT = input.GCHAT
    new_input.WD_ZERO = input.WD_ZERO
    new_input.LINSKER_FACTOR = input.LINSKER_FACTOR
    new_input.GRADB_FACTOR = input.GRADB_FACTOR
    new_input.FILTER = input.FILTER
    new_input.THETA_TRAPPED = input.THETA_TRAPPED
    new_input.SMALL = input.SMALL
    new_input.USE_TRANSPORT_MODEL = input.USE_TRANSPORT_MODEL

    return new_input
end