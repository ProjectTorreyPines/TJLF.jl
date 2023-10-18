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


    # SPECIES::Vector{Species{T}} = Vector{Species{T}}()
    ZS::Union{Vector{T},Missing} = missing
    MASS::Union{Vector{T},Missing} = missing
    RLNS::Union{Vector{T},Missing} = missing
    RLTS::Union{Vector{T},Missing} = missing
    TAUS::Union{Vector{T},Missing} = missing
    AS::Union{Vector{T},Missing} = missing
    VPAR::Union{Vector{T},Missing} = missing
    VPAR_SHEAR::Union{Vector{T},Missing} = missing
    VNS_SHEAR::Union{Vector{T},Missing} = missing
    VTS_SHEAR::Union{Vector{T},Missing} = missing


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
#       Get from tjlf_hermite
##########################################################

mutable struct OutputHermite{T<:Real}

    x::Vector{T}
    wx::Vector{T}
    h::Matrix{T}
    
end

##########################################################
#       Get from xgrid_functions_geo
##########################################################

mutable struct OutputGeometry{T<:Real}

    kx0_e::T

    fts::Vector{T}

    kxx::Vector{T}
    
    wdx::Vector{T}
    wdpx::Vector{T}

    b0x::Vector{T}
    b2x::Vector{T}

    cx_tor_par::Vector{T}
    cx_tor_per::Vector{T}
    cx_par_par::Vector{T}
    
end

mutable struct SaturationParameters{T<:Real}

    SAT_geo0::T
    SAT_geo1::T
    SAT_geo2::T

    y::Vector{T}
    R_unit::T
    B_unit::T
    q_unit::T

    R::Vector{T}
    Bp::Vector{T}
    Bt::Vector{T}

    Bt0::T
    grad_r0::T

    S_prime::Vector{T}
    kx_factor::Vector{T}

    B_geo::Vector{T}
    qrat_geo::Vector{T}

    sintheta_geo::Vector{T}
    costheta_geo::Vector{T}
    costheta_p_geo::Vector{T}

    theta::Vector{T}
    
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

    gradB::Matrix{T}
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
        gradB = zeros(T, nbasis, nbasis)
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
        gradB,b0,b0inv,lnB,p0,p0inv,bp,bpinv,
        c_par_par,c_tor_par,c_tor_per,
        kpar,modkpar,kpar_eff,modkpar_eff)
    end
end


mutable struct AveH{T<:Real}

    # average h-bessel functions
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



mutable struct AveWH{T<:Real}

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

    function AveWH{T}(ns::Int, nbasis::Int) where T<:Real
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



mutable struct AveKH

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

    function AveKH(ns::Int, nbasis::Int)
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
            kparhp1b0,kparhr11b0,kparhr13b0,
            kparhnbp,kparhp3bp,kparhp1bp,kparhr11bp,kparhr13bp
            )
    end
end


mutable struct AveG{T<:Real}

    # average g-bessel functions
    gn::Array{T}
    gp1::Array{T}
    gp3::Array{T}
    gr11::Array{T}
    gr13::Array{T}
    gr33::Array{T}
    gw113::Array{T}
    gw133::Array{T}
    gw333::Array{T}

    gt1::Array{T}
    gt3::Array{T}
    gu1::Array{T}
    gu3::Array{T}
    gu33::Array{T}
    gu3gt1::Array{T}
    gu3gt3::Array{T}
    gu33gt1::Array{T}
    gu33gt3::Array{T}

    gnp0::Array{T}
    gp1p0::Array{T}
    gp3p0::Array{T}
    gr11p0::Array{T}
    gr13p0::Array{T}
    gr33p0::Array{T}
    c_tor_par_gp1p0::Array{T}
    c_tor_par_gr11p0::Array{T}   
    c_tor_par_gr13p0::Array{T}

    gnb0::Array{T}
    gp1b0::Array{T}
    gp3b0::Array{T}
    gr11b0::Array{T}
    gr13b0::Array{T}
    gr33b0::Array{T}
    gw113b0::Array{T}
    gw133b0::Array{T}
    gw333b0::Array{T}

    gnbp::Array{T}
    gp1bp::Array{T}
    gp3bp::Array{T}
    gr11bp::Array{T}
    gr13bp::Array{T}
    gr33bp::Array{T}
    gw113bp::Array{T}
    gw133bp::Array{T}
    gw333bp::Array{T}

    function AveG{T}(ns::Int, nbasis::Int) where T<:Real
        gn = zeros(T, ns, nbasis, nbasis)
        gp1 = zeros(T, ns, nbasis, nbasis)
        gp3 = zeros(T, ns, nbasis, nbasis)
        gr11 = zeros(T, ns, nbasis, nbasis)
        gr13 = zeros(T, ns, nbasis, nbasis)
        gr33 = zeros(T, ns, nbasis, nbasis)
        gw113 = zeros(T, ns, nbasis, nbasis)
        gw133 = zeros(T, ns, nbasis, nbasis)
        gw333 = zeros(T, ns, nbasis, nbasis)

        gt1 = zeros(T, ns, nbasis, nbasis)
        gt3 = zeros(T, ns, nbasis, nbasis)
        gu1 = zeros(T, ns, nbasis, nbasis)
        gu3 = zeros(T, ns, nbasis, nbasis)
        gu33 = zeros(T, ns, nbasis, nbasis)
        gu3gt1 = zeros(T, ns, nbasis, nbasis)
        gu3gt3 = zeros(T, ns, nbasis, nbasis)
        gu33gt1 = zeros(T, ns, nbasis, nbasis)
        gu33gt3 = zeros(T, ns, nbasis, nbasis)

        gnp0 = zeros(T, ns, nbasis, nbasis)
        gp1p0 = zeros(T, ns, nbasis, nbasis)
        gp3p0 = zeros(T, ns, nbasis, nbasis)
        gr11p0 = zeros(T, ns, nbasis, nbasis)
        gr13p0 = zeros(T, ns, nbasis, nbasis)
        gr33p0 = zeros(T, ns, nbasis, nbasis)
        c_tor_par_gp1p0 = zeros(T, ns, nbasis, nbasis)
        c_tor_par_gr11p0 = zeros(T, ns, nbasis, nbasis)   
        c_tor_par_gr13p0 = zeros(T, ns, nbasis, nbasis)

        gnb0 = zeros(T, ns, nbasis, nbasis)
        gp1b0 = zeros(T, ns, nbasis, nbasis)
        gp3b0 = zeros(T, ns, nbasis, nbasis)
        gr11b0 = zeros(T, ns, nbasis, nbasis)
        gr13b0 = zeros(T, ns, nbasis, nbasis)
        gr33b0 = zeros(T, ns, nbasis, nbasis)
        gw113b0 = zeros(T, ns, nbasis, nbasis)
        gw133b0 = zeros(T, ns, nbasis, nbasis)
        gw333b0 = zeros(T, ns, nbasis, nbasis)

        gnbp = zeros(T, ns, nbasis, nbasis)
        gp1bp = zeros(T, ns, nbasis, nbasis)
        gp3bp = zeros(T, ns, nbasis, nbasis)
        gr11bp = zeros(T, ns, nbasis, nbasis)
        gr13bp = zeros(T, ns, nbasis, nbasis)
        gr33bp = zeros(T, ns, nbasis, nbasis)
        gw113bp = zeros(T, ns, nbasis, nbasis)
        gw133bp = zeros(T, ns, nbasis, nbasis)
        gw333bp = zeros(T, ns, nbasis, nbasis)

        new(
            gn,gp1,gp3,gr11,gr13,gr33,gw113,gw133,gw333,
            gt1,gt3,gu1,gu3,gu33,gu3gt1,gu3gt3,gu33gt1,gu33gt3,
            gnp0,gp1p0,gp3p0,gr11p0,gr13p0,gr33p0,c_tor_par_gp1p0,c_tor_par_gr11p0,c_tor_par_gr13p0,
            gnb0,gp1b0,gp3b0,gr11b0,gr13b0,gr33b0,gw113b0,gw133b0,gw333b0,
            gnbp,gp1bp,gp3bp,gr11bp,gr13bp,gr33bp,gw113bp,gw133bp,gw333bp
            )
    end
end


mutable struct AveWG{T<:Real}

    wdgp1p0::Array{T}
    wdgr11p0::Array{T}
    wdgr13p0::Array{T}
    wdgu1::Array{T}
    wdgu3::Array{T}
    wdgu33::Array{T}
    wdgt1::Array{T}
    wdgt3::Array{T}
    wdgu3gt1::Array{T}
    wdgu3gt3::Array{T}
    wdgu33gt1::Array{T}
    wdgu33gt3::Array{T}

    modwdgu1::Array{T}
    modwdgu3::Array{T}
    modwdgu33::Array{T}
    modwdgt1::Array{T}
    modwdgt3::Array{T}
    modwdgu3gt1::Array{T}
    modwdgu3gt3::Array{T}
    modwdgu33gt1::Array{T}
    modwdgu33gt3::Array{T}

    wdgp1b0::Array{T}
    wdgr11b0::Array{T}
    wdgr13b0::Array{T}

    wdgp1bp::Array{T}
    wdgr11bp::Array{T}
    wdgr13bp::Array{T}

    function AveWG{T}(ns::Int, nbasis::Int) where T<:Real
        wdgp1p0 = zeros(T, ns, nbasis, nbasis)
        wdgr11p0 = zeros(T, ns, nbasis, nbasis)
        wdgr13p0 = zeros(T, ns, nbasis, nbasis)
        wdgu1 = zeros(T, ns, nbasis, nbasis)
        wdgu3 = zeros(T, ns, nbasis, nbasis)
        wdgu33 = zeros(T, ns, nbasis, nbasis)
        wdgt1 = zeros(T, ns, nbasis, nbasis)
        wdgt3 = zeros(T, ns, nbasis, nbasis)
        wdgu3gt1 = zeros(T, ns, nbasis, nbasis)
        wdgu3gt3 = zeros(T, ns, nbasis, nbasis)
        wdgu33gt1 = zeros(T, ns, nbasis, nbasis)
        wdgu33gt3 = zeros(T, ns, nbasis, nbasis)

        modwdgu1 = zeros(T, ns, nbasis, nbasis)
        modwdgu3 = zeros(T, ns, nbasis, nbasis)
        modwdgu33 = zeros(T, ns, nbasis, nbasis)
        modwdgt1 = zeros(T, ns, nbasis, nbasis)
        modwdgt3 = zeros(T, ns, nbasis, nbasis)
        modwdgu3gt1 = zeros(T, ns, nbasis, nbasis)
        modwdgu3gt3 = zeros(T, ns, nbasis, nbasis)
        modwdgu33gt1 = zeros(T, ns, nbasis, nbasis)
        modwdgu33gt3 = zeros(T, ns, nbasis, nbasis)

        wdgp1b0 = zeros(T, ns, nbasis, nbasis)
        wdgr11b0 = zeros(T, ns, nbasis, nbasis)
        wdgr13b0 = zeros(T, ns, nbasis, nbasis)

        wdgp1bp = zeros(T, ns, nbasis, nbasis)
        wdgr11bp = zeros(T, ns, nbasis, nbasis)
        wdgr13bp = zeros(T, ns, nbasis, nbasis)

        new(
            wdgp1p0,wdgr11p0,wdgr13p0,wdgu1,wdgu3,wdgu33,wdgt1,wdgt3,wdgu3gt1,wdgu3gt3,wdgu33gt1,wdgu33gt3,
            modwdgu1,modwdgu3,modwdgu33,modwdgt1,modwdgt3,modwdgu3gt1,modwdgu3gt3,modwdgu33gt1,modwdgu33gt3,
            wdgp1b0,wdgr11b0,wdgr13b0,
            wdgp1bp,wdgr11bp,wdgr13bp
            )
    end
end

mutable struct AveKG

    kpargnp0::Array{ComplexF64}
    kpargp1p0::Array{ComplexF64}
    kpargp3p0::Array{ComplexF64}
    kpargu1::Array{ComplexF64}
    kpargu3::Array{ComplexF64}
    kpargt1::Array{ComplexF64}
    kpargt3::Array{ComplexF64}
    modkpargu1::Array{ComplexF64}
    modkpargu3::Array{ComplexF64}

    kpargp1b0::Array{ComplexF64}
    kpargr11b0::Array{ComplexF64}
    kpargr13b0::Array{ComplexF64}

    kpargnbp::Array{ComplexF64}
    kpargp3bp::Array{ComplexF64}
    kpargp1bp::Array{ComplexF64}
    kpargr11bp::Array{ComplexF64}
    kpargr13bp::Array{ComplexF64}

    function AveKG(ns::Int, nbasis::Int)
        kpargnp0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargp1p0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargp3p0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargu1 = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargu3 = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargt1 = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargt3 = zeros(ComplexF64, ns, nbasis, nbasis)
        modkpargu1 = zeros(ComplexF64, ns, nbasis, nbasis)
        modkpargu3 = zeros(ComplexF64, ns, nbasis, nbasis)

        kpargp1b0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargr11b0 = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargr13b0 = zeros(ComplexF64, ns, nbasis, nbasis)
                
        kpargnbp = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargp3bp = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargp1bp = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargr11bp = zeros(ComplexF64, ns, nbasis, nbasis)
        kpargr13bp = zeros(ComplexF64, ns, nbasis, nbasis)

        new(
            kpargnp0,kpargp1p0,kpargp3p0,kpargu1,kpargu3,kpargt1,kpargt3,modkpargu1,modkpargu3,
            kpargp1b0,kpargr11b0,kpargr13b0,
            kpargnbp,kpargp3bp,kpargp1bp,kpargr11bp,kpargr13bp
            )
    end
end