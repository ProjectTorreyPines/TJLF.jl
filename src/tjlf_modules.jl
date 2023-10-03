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



    SPECIES::Vector{Species{T}} = Vector{Species{T}}()



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
    WDIA_TRAPPED::Union{T,Missing} = missing
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


##########################################################
#       Get from xgrid_functions_geo
##########################################################

mutable struct OutputGeometry{T<:Real}

    kxx::Vector{T}
    
    wdx::Vector{T}
    wdpx::Vector{T}

    b0x::Vector{T}
    b2x::Vector{T}

    cx_tor_par::Vector{T}
    cx_tor_per::Vector{T}
    cx_par_par::Vector{T}
    
end

##########################################################
#       Get from tjlf_hermite
##########################################################

mutable struct OutputHermite{T<:Real}

    x::Vector{T}
    wx::Vector{T}
    h::Matrix{T}
    
end


##########################################################
#       Used for get_matrx()
##########################################################

mutable struct Ave{T<:Real}
    kx::Matrix{T}

    wdh::Matrix{T}
    modwdh::Matrix{T}
    wdg::Matrix{T}
    modwdg::Matrix{T}

    b0::Matrix{T}
    b0inv::Matrix{T}
    lnB::Matrix{T}
    p0::Matrix{T}
    p0inv::Matrix{T}
    bp::Matrix{T}
    bpinv::Matrix{T}

    c_par_par::Matrix{T}
    c_tor_par::Matrix{T}
    c_tor_per::Matrix{T}

    kpar::Matrix{T}
    modkpar::Matrix{T}
    kpar_eff::Array{ComplexF64}
    modkpar_eff::Array{ComplexF64}

    function Ave{T}(ns::Int, nbasis::Int) where T<:Real
        kx = zeros(T, nbasis, nbasis)
        wdh = zeros(T, nbasis, nbasis)
        modwdh = zeros(T, nbasis, nbasis)
        wdg = zeros(T, nbasis, nbasis)
        modwdg = zeros(T, nbasis, nbasis)
        b0 = zeros(T, nbasis, nbasis)
        b0inv = zeros(T, nbasis, nbasis)
        lnB = zeros(T, nbasis, nbasis)
        p0 = zeros(T, nbasis, nbasis)
        p0inv = zeros(T, nbasis, nbasis)
        bp = zeros(T, nbasis, nbasis)
        bpinv = zeros(T, nbasis, nbasis)
        c_par_par = zeros(T, nbasis, nbasis)
        c_tor_par = zeros(T, nbasis, nbasis)
        c_tor_per = zeros(T, nbasis, nbasis)
        kpar = zeros(T, nbasis, nbasis)
        modkpar = zeros(T, nbasis, nbasis)
        kpar_eff = Array{ComplexF64}(undef, ns, nbasis, nbasis)
        modkpar_eff = Array{ComplexF64}(undef, ns, nbasis, nbasis)


        new(kx,
        wdh,modwdh,wdg,modwdg,
        b0,b0inv,lnB,p0,p0inv,bp,bpinv,
        c_par_par,c_tor_par,c_tor_per,
        kpar,modkpar,kpar_eff,modkpar_eff)
    end
end


mutable struct AveH{T<:Real}

    hn::Array{T}
    hp1::Array{T}
    hp3::Array{T}
    hr11::Array{T}
    hr13::Array{T}
    hr33::Array{T}
    hw113::Array{T}
    hw133::Array{T}
    hw333::Array{T}

    ht1::Array{T}
    ht3::Array{T}
    hu1::Array{T}
    hu3::Array{T}
    hu33::Array{T}
    hu3ht1::Array{T}
    hu33ht1::Array{T}
    hu3ht3::Array{T}
    hu33ht3::Array{T}

    hnp0::Array{T}
    hp1p0::Array{T}
    hp3p0::Array{T}
    hr11p0::Array{T}
    hr13p0::Array{T}
    hr33p0::Array{T}
    c_tor_par_hp1p0::Array{T}
    c_tor_par_hr11p0::Array{T}
    c_tor_par_hr13p0::Array{T}

    hnb0::Array{T}
    hp1b0::Array{T}
    hp3b0::Array{T}
    hr11b0::Array{T}
    hr13b0::Array{T}
    hr33b0::Array{T}
    hw113b0::Array{T}
    hw133b0::Array{T}
    hw333b0::Array{T}

    hnbp::Array{T}
    hp1bp::Array{T}
    hp3bp::Array{T}
    hr11bp::Array{T}
    hr13bp::Array{T}
    hr33bp::Array{T}
    hw113bp::Array{T}
    hw133bp::Array{T}
    hw333bp::Array{T}

    function AveH{T}(ns::Int, nbasis::Int) where T<:Real
        hn = zeros(T, ns, nbasis, nbasis)
        hp1 = zeros(T, ns, nbasis, nbasis)
        hp3 = zeros(T, ns, nbasis, nbasis)
        hr11 = zeros(T, ns, nbasis, nbasis)
        hr13 = zeros(T, ns, nbasis, nbasis)
        hr33 = zeros(T, ns, nbasis, nbasis)
        hw113 = zeros(T, ns, nbasis, nbasis)
        hw133 = zeros(T, ns, nbasis, nbasis)
        hw333 = zeros(T, ns, nbasis, nbasis)
        ht1 = zeros(T, ns, nbasis, nbasis)
        ht3 = zeros(T, ns, nbasis, nbasis)
        hu1 = zeros(T, ns, nbasis, nbasis)
        hu3 = zeros(T, ns, nbasis, nbasis)
        hu33 = zeros(T, ns, nbasis, nbasis)
        hu3ht1 = zeros(T, ns, nbasis, nbasis)
        hu33ht1 = zeros(T, ns, nbasis, nbasis)
        hu3ht3 = zeros(T, ns, nbasis, nbasis)
        hu33ht3 = zeros(T, ns, nbasis, nbasis)
        hnp0 = zeros(T, ns, nbasis, nbasis)
        hp1p0 = zeros(T, ns, nbasis, nbasis)
        hp3p0 = zeros(T, ns, nbasis, nbasis)
        hr11p0 = zeros(T, ns, nbasis, nbasis)
        hr13p0 = zeros(T, ns, nbasis, nbasis)
        hr33p0 = zeros(T, ns, nbasis, nbasis)
        c_tor_par_hp1p0 = zeros(T, ns, nbasis, nbasis)
        c_tor_par_hr11p0 = zeros(T, ns, nbasis, nbasis)
        c_tor_par_hr13p0 = zeros(T, ns, nbasis, nbasis)
        hnb0 = zeros(T, ns, nbasis, nbasis)
        hp1b0 = zeros(T, ns, nbasis, nbasis)
        hp3b0 = zeros(T, ns, nbasis, nbasis)
        hr11b0 = zeros(T, ns, nbasis, nbasis)
        hr13b0 = zeros(T, ns, nbasis, nbasis)
        hr33b0 = zeros(T, ns, nbasis, nbasis)
        hw113b0 = zeros(T, ns, nbasis, nbasis)
        hw133b0 = zeros(T, ns, nbasis, nbasis)
        hw333b0 = zeros(T, ns, nbasis, nbasis)
        hnbp = zeros(T, ns, nbasis, nbasis)
        hp1bp = zeros(T, ns, nbasis, nbasis)
        hp3bp = zeros(T, ns, nbasis, nbasis)
        hr11bp = zeros(T, ns, nbasis, nbasis)
        hr13bp = zeros(T, ns, nbasis, nbasis)
        hr33bp = zeros(T, ns, nbasis, nbasis)
        hw113bp = zeros(T, ns, nbasis, nbasis)
        hw133bp = zeros(T, ns, nbasis, nbasis)
        hw333bp = zeros(T, ns, nbasis, nbasis)

        new(
            hn,hp1,hp3,hr11,hr13,hr33,hw113,hw133,hw333,
            ht1,ht3,hu1,hu3,hu33,hu3ht1,hu33ht1,hu3ht3,hu33ht3,
            hnp0,hp1p0,hp3p0,hr11p0,hr13p0,hr33p0,c_tor_par_hp1p0,c_tor_par_hr11p0,c_tor_par_hr13p0,
            hnb0,hp1b0,hp3b0,hr11b0,hr13b0,hr33b0,hw113b0,hw133b0,hw333b0,
            hnbp,hp1bp,hp3bp,hr11bp,hr13bp,hr33bp,hw113bp,hw133bp,hw333bp
            )
    end

end



mutable struct AveW{T<:Real}

    wdhp1p0::Array{T}
    wdhr11p0::Array{T}
    wdhr13p0::Array{T}
    wdht1::Array{T}
    wdht3::Array{T}
    wdhu1::Array{T}
    wdhu3::Array{T}
    wdhu3ht1::Array{T}
    wdhu3ht3::Array{T}
    wdhu33::Array{T}
    wdhu33ht1::Array{T}
    wdhu33ht3::Array{T}

    modwdht1::Array{T}
    modwdht3::Array{T}
    modwdhu1::Array{T}
    modwdhu3::Array{T}
    modwdhu3ht1::Array{T}
    modwdhu3ht3::Array{T}
    modwdhu33::Array{T}
    modwdhu33ht1::Array{T}
    modwdhu33ht3::Array{T}

    wdhp1b0::Array{T}
    wdhr11b0::Array{T}
    wdhr13b0::Array{T}

    wdhp1bp::Array{T}
    wdhr11bp::Array{T}
    wdhr13bp::Array{T}

    function AveW{T}(ns::Int, nbasis::Int) where T<:Real
        wdhp1p0 = zeros(T, ns, nbasis, nbasis)
        wdhr11p0 = zeros(T, ns, nbasis, nbasis)
        wdhr13p0 = zeros(T, ns, nbasis, nbasis)
        wdht1 = zeros(T, ns, nbasis, nbasis)
        wdht3 = zeros(T, ns, nbasis, nbasis)
        wdhu1 = zeros(T, ns, nbasis, nbasis)
        wdhu3 = zeros(T, ns, nbasis, nbasis)
        wdhu3ht1 = zeros(T, ns, nbasis, nbasis)
        wdhu3ht3 = zeros(T, ns, nbasis, nbasis)
        wdhu33 = zeros(T, ns, nbasis, nbasis)
        wdhu33ht1 = zeros(T, ns, nbasis, nbasis)
        wdhu33ht3 = zeros(T, ns, nbasis, nbasis)

        modwdht1 = zeros(T, ns, nbasis, nbasis)
        modwdht3 = zeros(T, ns, nbasis, nbasis)
        modwdhu1 = zeros(T, ns, nbasis, nbasis)
        modwdhu3 = zeros(T, ns, nbasis, nbasis)
        modwdhu3ht1 = zeros(T, ns, nbasis, nbasis)
        modwdhu3ht3 = zeros(T, ns, nbasis, nbasis)
        modwdhu33 = zeros(T, ns, nbasis, nbasis)
        modwdhu33ht1 = zeros(T, ns, nbasis, nbasis)
        modwdhu33ht3 = zeros(T, ns, nbasis, nbasis)

        wdhp1b0 = zeros(T, ns, nbasis, nbasis)
        wdhr11b0 = zeros(T, ns, nbasis, nbasis)
        wdhr13b0 = zeros(T, ns, nbasis, nbasis)

        wdhp1bp = zeros(T, ns, nbasis, nbasis)
        wdhr11bp = zeros(T, ns, nbasis, nbasis)
        wdhr13bp = zeros(T, ns, nbasis, nbasis)

        new(
            wdhp1p0,wdhr11p0,wdhr13p0,wdht1,wdht3,wdhu1,wdhu3,wdhu3ht1,wdhu3ht3,wdhu33,wdhu33ht1,wdhu33ht3,
            modwdht1,modwdht3,modwdhu1,modwdhu3,modwdhu3ht1,modwdhu3ht3,modwdhu33,modwdhu33ht1,modwdhu33ht3,
            wdhp1b0,wdhr11b0,wdhr13b0,
            wdhp1bp,wdhr11bp,wdhr13bp
            )
    end

end



mutable struct AveK

    kparhnp0::Array{ComplexF64}
    kparhp1p0::Array{ComplexF64}
    kparhp3p0::Array{ComplexF64}
    kparhu1::Array{ComplexF64}
    kparhu3::Array{ComplexF64}
    kparht1::Array{ComplexF64}
    kparht3::Array{ComplexF64}
    modkparhu1::Array{ComplexF64}
    modkparhu3::Array{ComplexF64}

    kparhp1b0::Array{ComplexF64}
    kparhr11b0::Array{ComplexF64}
    kparhr13b0::Array{ComplexF64}

    kparhnbp::Array{ComplexF64}
    kparhp3bp::Array{ComplexF64}
    kparhp1bp::Array{ComplexF64}
    kparhr11bp::Array{ComplexF64}
    kparhr13bp::Array{ComplexF64}

    function AveK(ns::Int, nbasis::Int)
        kparhnp0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhp1p0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhp3p0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhu1 = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhu3 = zeros(ComplexF64, ns, nbasis, nbasis)
        kparht1 = zeros(ComplexF64, ns, nbasis, nbasis)
        kparht3 = zeros(ComplexF64, ns, nbasis, nbasis)
        modkparhu1 = zeros(ComplexF64, ns, nbasis, nbasis)
        modkparhu3 = zeros(ComplexF64, ns, nbasis, nbasis)

        kparhp1b0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhr11b0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhr13b0 = zeros(ComplexF64, ns, nbasis, nbasis)
                
        kparhnbp = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhp3bp = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhp1bp = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhr11bp = zeros(ComplexF64, ns, nbasis, nbasis)
        kparhr13bp = zeros(ComplexF64, ns, nbasis, nbasis)

        new(
            kparhnp0,kparhp1p0,kparhp3p0,kparhu1,kparhu3,kparht1,kparht3,modkparhu1,modkparhu3,
            kparhp1b0,kparhr11b0kparhr13b0,
            kparhnbp,kparhp3bp,kparhp1bp,kparhr11bp,kparhr13bp
            )
    end

end