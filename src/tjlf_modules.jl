Base.@kwdef mutable struct Species{T<:Real}
    ZS::Union{T,Missing} = missing
    MASS::Union{T,Missing} = missing
    RLNS::Union{T,Missing} = missing
    RLTS::Union{T,Missing} = missing
    TAUS::Union{T,Missing} = missing
    AS::Union{T,Missing} = missing
    VPAR::Union{T,Missing} = missing
    VPAR_SHEAR::Union{T,Missing} = missing
    VNS_SHEAR::Union{T,Missing} = missing
    VTS_SHEAR::Union{T,Missing} = missing
end


Base.@kwdef mutable struct InputTJLF{T<:Real}

    UNITS::Union{String,Missing} = missing
    NS::Union{Int,Missing} = missing

    USE_TRANSPORT_MODEL::Union{Bool,Missing} = missing
    GEOMETRY_FLAG::Union{Int,Missing} = missing
    USE_BPER::Union{Bool,Missing} = missing
    USE_BPAR::Union{Bool,Missing} = missing
    USE_MHD_RULE::Union{Bool,Missing} = missing
    USE_BISECTION::Union{Bool,Missing} = missing
    USE_INBOARD_DETRAPPED::Union{Bool,Missing} = missing
    USE_AVE_ION_GRID::Union{Bool,Missing} = missing

    SAT_RULE::Union{Int,Missing} = missing
    KYGRID_MODEL::Union{Int,Missing} = missing
    XNU_MODEL::Union{Int,Missing} = missing
    VPAR_MODEL::Union{Int,Missing} = missing
    VPAR_SHEAR_MODEL::Union{Int,Missing} = missing



    SPECIES::Vector{Species{T}}=Vector{Species{T}}()



    SIGN_BT::Union{T,Missing} = missing
    SIGN_IT::Union{T,Missing} = missing
    KY::Union{T,Missing} = missing
    NEW_EIKONAL::Union{Bool,Missing} = missing
    VEXB::Union{T,Missing} = missing
    VEXB_SHEAR::Union{T,Missing} = missing
    BETAE::Union{T,Missing} = missing
    XNUE::Union{T,Missing} = missing
    ZEFF::Union{T,Missing} = missing
    DEBYE::Union{T,Missing} = missing

    IFLUX::Union{Bool,Missing} = missing
    IBRANCH::Union{Int,Missing} = missing

    NMODES::Union{Int,Missing} = missing
    NBASIS_MAX::Union{Int,Missing} = missing
    NBASIS_MIN::Union{Int,Missing} = missing
    NXGRID::Union{Int,Missing} = missing
    NKY::Union{Int,Missing} = missing

    ADIABATIC_ELEC::Union{Bool,Missing} = missing
    ALPHA_MACH::Union{T,Missing} = missing
    ALPHA_E::Union{T,Missing} = missing
    ALPHA_P::Union{T,Missing} = missing
    ALPHA_QUENCH::Union{T,Missing} = missing
    ALPHA_ZF::Union{T,Missing} = missing
    XNU_FACTOR::Union{T,Missing} = missing
    DEBYE_FACTOR::Union{T,Missing} = missing
    ETG_FACTOR::Union{T,Missing} = missing
    RLNP_CUTOFF::Union{T,Missing} = missing

    WRITE_WAVEFUNCTION_FLAG::Union{Int,Missing} = missing

    WIDTH::Union{T,Missing} = missing
    WIDTH_MIN::Union{T,Missing} = missing
    FIND_WIDTH::Union{Bool,Missing} = missing
    NWIDTH::Union{Int,Missing} = missing
    RMIN_LOC::Union{T,Missing} = missing
    RMAJ_LOC::Union{T,Missing} = missing
    ZMAJ_LOC::Union{T,Missing} = missing
    DRMINDX_LOC::Union{T,Missing} = missing
    DRMAJDX_LOC::Union{T,Missing} = missing
    DZMAJDX_LOC::Union{T,Missing} = missing
    Q_LOC::Union{T,Missing} = missing
    KAPPA_LOC::Union{T,Missing} = missing
    S_KAPPA_LOC::Union{T,Missing} = missing
    DELTA_LOC::Union{T,Missing} = missing
    S_DELTA_LOC::Union{T,Missing} = missing
    ZETA_LOC::Union{T,Missing} = missing
    S_ZETA_LOC::Union{T,Missing} = missing
    P_PRIME_LOC::Union{T,Missing} = missing
    Q_PRIME_LOC::Union{T,Missing} = missing
    BETA_LOC::Union{T,Missing} = missing
    KX0_LOC::Union{T,Missing} = missing
    RMIN_SA::Union{T,Missing} = missing
    RMAJ_SA::Union{T,Missing} = missing
    Q_SA::Union{T,Missing} = missing
    SHAT_SA::Union{T,Missing} = missing
    ALPHA_SA::Union{T,Missing} = missing
    XWELL_SA::Union{T,Missing} = missing
    THETA0_SA::Union{T,Missing} = missing
    B_MODEL_SA::Union{Int,Missing} = missing
    FT_MODEL_SA::Union{Int,Missing} = missing
    
    DAMP_PSI::Union{T,Missing} = missing
    DAMP_SIG::Union{T,Missing} = missing
    PARK::Union{T,Missing} = missing
    GHAT::Union{T,Missing} = missing
    GCHAT::Union{T,Missing} = missing
    WD_ZERO::Union{T,Missing} = missing
    LINSKER_FACTOR::Union{T,Missing} = missing
    GRADB_FACTOR::Union{T,Missing} = missing
    FILTER::Union{T,Missing} = missing
    THETA_TRAPPED::Union{T,Missing} = missing
    NN_MAX_ERROR::Union{T,Missing} = missing
end

# # struct RunTGLF
# #     input::InputTGLF
# #     geometry::GeoTGLF
# #     ....
# # end

# # run_tglf = RunTGLF()
# # run_tglf.input.specie[1].AS
# # RMAJ_LOC  = run_tglf.input.RMAJ_LOC
# # run_tglf.geometry.B 


# # function maj_s(input::InputTGLF)
# #     if input.geometry...
# #         return input.RMAJOR
# #     else
# #         ..
# #     end
# # end



