Base.@kwdef mutable struct InputTGLF
    SIGN_BT::Union{Int,Missing} = missing
    SIGN_IT::Union{Int,Missing} = missing
    NS::Union{Int,Missing} = missing
    ZMAJ_LOC::Union{Float64,Missing} = missing
    DRMINDX_LOC::Union{Float64,Missing} = missing
    DZMAJDX_LOC::Union{Float64,Missing} = missing
    S_DELTA_LOC::Union{Float64,Missing} = missing
    ZETA_LOC::Union{Float64,Missing} = missing
    S_ZETA_LOC::Union{Float64,Missing} = missing

    MASS_1::Union{Float64,Missing} = missing
    ZS_1::Union{Float64,Missing} = missing
    AS_1::Union{Float64,Missing} = missing
    TAUS_1::Union{Float64,Missing} = missing

    MASS_2::Union{Float64,Missing} = missing
    ZS_2::Union{Float64,Missing} = missing
    VPAR_2::Union{Float64,Missing} = missing
    VPAR_SHEAR_2::Union{Float64,Missing} = missing

    MASS_3::Union{Float64,Missing} = missing
    ZS_3::Union{Float64,Missing} = missing
    RLTS_3::Union{Float64,Missing} = missing
    TAUS_3::Union{Float64,Missing} = missing
    VPAR_3::Union{Float64,Missing} = missing
    VPAR_SHEAR_3::Union{Float64,Missing} = missing

    # TGLF-NN uses 3 species
    # This is why parameters for species 1:3 are sorted differently than 4:10
    MASS_4::Union{Float64,Missing} = missing
    AS_4::Union{Float64,Missing} = missing
    ZS_4::Union{Float64,Missing} = missing
    RLNS_4::Union{Float64,Missing} = missing
    RLTS_4::Union{Float64,Missing} = missing
    TAUS_4::Union{Float64,Missing} = missing
    VPAR_4::Union{Float64,Missing} = missing
    VPAR_SHEAR_4::Union{Float64,Missing} = missing

    MASS_5::Union{Float64,Missing} = missing
    AS_5::Union{Float64,Missing} = missing
    ZS_5::Union{Float64,Missing} = missing
    RLNS_5::Union{Float64,Missing} = missing
    RLTS_5::Union{Float64,Missing} = missing
    TAUS_5::Union{Float64,Missing} = missing
    VPAR_5::Union{Float64,Missing} = missing
    VPAR_SHEAR_5::Union{Float64,Missing} = missing

    MASS_6::Union{Float64,Missing} = missing
    AS_6::Union{Float64,Missing} = missing
    ZS_6::Union{Float64,Missing} = missing
    RLNS_6::Union{Float64,Missing} = missing
    RLTS_6::Union{Float64,Missing} = missing
    TAUS_6::Union{Float64,Missing} = missing
    VPAR_6::Union{Float64,Missing} = missing
    VPAR_SHEAR_6::Union{Float64,Missing} = missing

    MASS_7::Union{Float64,Missing} = missing
    AS_7::Union{Float64,Missing} = missing
    ZS_7::Union{Float64,Missing} = missing
    RLNS_7::Union{Float64,Missing} = missing
    RLTS_7::Union{Float64,Missing} = missing
    TAUS_7::Union{Float64,Missing} = missing
    VPAR_7::Union{Float64,Missing} = missing
    VPAR_SHEAR_7::Union{Float64,Missing} = missing

    MASS_8::Union{Float64,Missing} = missing
    AS_8::Union{Float64,Missing} = missing
    ZS_8::Union{Float64,Missing} = missing
    RLNS_8::Union{Float64,Missing} = missing
    RLTS_8::Union{Float64,Missing} = missing
    TAUS_8::Union{Float64,Missing} = missing
    VPAR_8::Union{Float64,Missing} = missing
    VPAR_SHEAR_8::Union{Float64,Missing} = missing

    MASS_9::Union{Float64,Missing} = missing
    AS_9::Union{Float64,Missing} = missing
    ZS_9::Union{Float64,Missing} = missing
    RLNS_9::Union{Float64,Missing} = missing
    RLTS_9::Union{Float64,Missing} = missing
    TAUS_9::Union{Float64,Missing} = missing
    VPAR_9::Union{Float64,Missing} = missing
    VPAR_SHEAR_9::Union{Float64,Missing} = missing

    MASS_10::Union{Float64,Missing} = missing
    AS_10::Union{Float64,Missing} = missing
    ZS_10::Union{Float64,Missing} = missing
    RLNS_10::Union{Float64,Missing} = missing
    RLTS_10::Union{Float64,Missing} = missing
    TAUS_10::Union{Float64,Missing} = missing
    VPAR_10::Union{Float64,Missing} = missing
    VPAR_SHEAR_10::Union{Float64,Missing} = missing

    AS_2::Union{Float64,Missing} = missing
    AS_3::Union{Float64,Missing} = missing
    BETAE::Union{Float64,Missing} = missing
    DEBYE::Union{Float64,Missing} = missing
    DELTA_LOC::Union{Float64,Missing} = missing
    DRMAJDX_LOC::Union{Float64,Missing} = missing
    KAPPA_LOC::Union{Float64,Missing} = missing
    P_PRIME_LOC::Union{Float64,Missing} = missing
    Q_LOC::Union{Float64,Missing} = missing
    Q_PRIME_LOC::Union{Float64,Missing} = missing
    RLNS_1::Union{Float64,Missing} = missing
    RLNS_2::Union{Float64,Missing} = missing
    RLNS_3::Union{Float64,Missing} = missing
    RLTS_1::Union{Float64,Missing} = missing
    RLTS_2::Union{Float64,Missing} = missing
    RMAJ_LOC::Union{Float64,Missing} = missing
    RMIN_LOC::Union{Float64,Missing} = missing
    S_KAPPA_LOC::Union{Float64,Missing} = missing
    TAUS_2::Union{Float64,Missing} = missing
    VEXB_SHEAR::Union{Float64,Missing} = missing
    VPAR_1::Union{Float64,Missing} = missing
    VPAR_SHEAR_1::Union{Float64,Missing} = missing
    XNUE::Union{Float64,Missing} = missing
    ZEFF::Union{Float64,Missing} = missing

    # switches
    UNITS::Union{String,Missing} = missing
    ALPHA_ZF::Union{Float64,Missing} = missing
    USE_MHD_RULE::Union{Bool,Missing} = missing
    NKY::Union{Int,Missing} = missing
    SAT_RULE::Union{Int,Missing} = missing
    KYGRID_MODEL::Union{Int,Missing} = missing
    NMODES::Union{Int,Missing} = missing
    NBASIS_MIN::Union{Int,Missing} = missing
    NBASIS_MAX::Union{Int,Missing} = missing
    XNU_MODEL::Union{Int,Missing} = missing
    USE_AVE_ION_GRID::Union{Bool,Missing} = missing
    ALPHA_QUENCH::Union{Int,Missing} = missing
    ALPHA_MACH::Union{Float64,Missing} = missing
    WDIA_TRAPPED::Union{Float64,Missing} = missing
    USE_BPAR::Union{Bool,Missing} = missing
    USE_BPER::Union{Bool,Missing} = missing

    _Qgb::Union{Float64,Missing} = missing

    # missing
    USE_TRANSPORT_MODEL::Bool = true
    USE_BISECTION::Bool = true
    USE_INBOARD_DETRAPPED::Bool = false
    NEW_EIKONAL::Bool = true
    FIND_WIDTH::Bool = true
    IFLUX::Bool = true
    ADIABATIC_ELEC::Bool = false

    GEOMETRY_FLAG::Int = 1
    NWIDTH::Int = 21
    NXGRID::Int = 16
    VPAR_MODEL::Int = 0
    FT_MODEL_SA::Int = 1
    VPAR_SHEAR_MODEL::Int = 1
    IBRANCH::Int = -1
    WRITE_WAVEFUNCTION_FLAG::Int = 0

    VNS_SHEAR_1::Float64 = 0.0
    VNS_SHEAR_2::Float64 = 0.0
    VNS_SHEAR_3::Float64 = 0.0
    VNS_SHEAR_4::Float64 = 0.0
    VNS_SHEAR_5::Float64 = 0.0
    VNS_SHEAR_6::Float64 = 0.0
    VNS_SHEAR_7::Float64 = 0.0
    VNS_SHEAR_8::Float64 = 0.0
    VNS_SHEAR_9::Float64 = 0.0
    VNS_SHEAR_10::Float64 = 0.0
    VTS_SHEAR_1::Float64 = 0.0
    VTS_SHEAR_2::Float64 = 0.0
    VTS_SHEAR_3::Float64 = 0.0
    VTS_SHEAR_4::Float64 = 0.0
    VTS_SHEAR_5::Float64 = 0.0
    VTS_SHEAR_6::Float64 = 0.0
    VTS_SHEAR_7::Float64 = 0.0
    VTS_SHEAR_8::Float64 = 0.0
    VTS_SHEAR_9::Float64 = 0.0
    VTS_SHEAR_10::Float64 = 0.0

    KY::Float64 = 0.3
    VEXB::Float64 = 0.0
    ALPHA_E::Float64 = 1.0
    ALPHA_P::Float64 = 1.0
    XNU_FACTOR::Float64 = 1.0
    DEBYE_FACTOR::Float64 = 1.0
    RLNP_CUTOFF::Float64 = 18.0
    WIDTH::Float64 = 1.65
    WIDTH_MIN::Float64 = 0.3
    BETA_LOC::Float64 = 1.0
    KX0_LOC::Float64 = 1.0
    RMIN_SA::Float64 = 0.5
    RMAJ_SA::Float64 = 3.0
    Q_SA::Float64 = 2.0
    SHAT_SA::Float64 = 1.0
    ALPHA_SA::Float64 = 0.0
    PARK::Float64 = 1.0
    GHAT::Float64 = 1.0
    GCHAT::Float64 = 1.0
    WD_ZERO::Float64 = 0.1
    LINSKER_FACTOR::Float64 = 0.0
    GRADB_FACTOR::Float64 = 0.0
    FILTER::Float64 = 2.0
    THETA_TRAPPED::Float64 = 0.7
    NN_MAX_ERROR::Float64 = -1.0

end












mutable struct InputTJLF{T<:Real}

    UNITS::String

    USE_TRANSPORT_MODEL::Union{Bool,Missing}
    USE_BPER::Union{Bool,Missing}
    USE_BPAR::Union{Bool,Missing}
    USE_MHD_RULE::Union{Bool,Missing}
    USE_BISECTION::Union{Bool,Missing}
    USE_INBOARD_DETRAPPED::Union{Bool,Missing}
    USE_AVE_ION_GRID::Union{Bool,Missing} # not used?
    NEW_EIKONAL::Union{Bool,Missing}
    FIND_WIDTH::Union{Bool,Missing}
    IFLUX::Union{Bool,Missing}
    ADIABATIC_ELEC::Union{Bool,Missing}

    GEOMETRY_FLAG::Union{Int,Missing} # used to specify geometry type
    SAT_RULE::Union{Int,Missing}
    NS::Union{Int,Missing}
    NMODES::Union{Int,Missing}
    NWIDTH::Union{Int,Missing}
    NBASIS_MAX::Union{Int,Missing}
    NBASIS_MIN::Union{Int,Missing}
    NXGRID::Union{Int,Missing}
    NKY::Union{Int,Missing}
    KYGRID_MODEL::Union{Int,Missing}
    XNU_MODEL::Union{Int,Missing}
    VPAR_MODEL::Union{Int,Missing}
    B_MODEL_SA::Union{Int,Missing} # used for SA geometry
    FT_MODEL_SA::Union{Int,Missing} # used for SA geometry
    VPAR_SHEAR_MODEL::Union{Int,Missing} # commented out in TGLF
    IBRANCH::Union{Int,Missing}
    WRITE_WAVEFUNCTION_FLAG::Union{Int,Missing} 

    ZS::Vector{T}
    MASS::Vector{T}
    RLNS::Vector{T}
    RLTS::Vector{T}
    TAUS::Vector{T}
    AS::Vector{T}
    VPAR::Vector{T}
    VPAR_SHEAR::Vector{T}
    VNS_SHEAR::Vector{T}
    VTS_SHEAR::Vector{T}

    WIDTH_SPECTRUM::Vector{T}
    KY_SPECTRUM::Vector{T}
    GAMMA_SPECTRUM::Vector{T}

    SIGN_BT::Int
    SIGN_IT::Int
    KY::T

    VEXB::T
    VEXB_SHEAR::T
    BETAE::T
    XNUE::T
    ZEFF::T
    DEBYE::T
    
    ALPHA_MACH::T
    ALPHA_E::T
    ALPHA_P::T
    ALPHA_QUENCH::Int
    ALPHA_ZF::T
    XNU_FACTOR::T
    DEBYE_FACTOR::T
    ETG_FACTOR::T
    RLNP_CUTOFF::T
    
    WIDTH::T
    WIDTH_MIN::T

    RMIN_LOC::T
    RMAJ_LOC::T
    ZMAJ_LOC::T
    DRMINDX_LOC::T
    DRMAJDX_LOC::T
    DZMAJDX_LOC::T
    Q_LOC::T
    KAPPA_LOC::T
    S_KAPPA_LOC::T
    DELTA_LOC::T
    S_DELTA_LOC::T
    ZETA_LOC::T
    S_ZETA_LOC::T
    P_PRIME_LOC::T
    Q_PRIME_LOC::T
    BETA_LOC::T
    KX0_LOC::T
    RMIN_SA::T
    RMAJ_SA::T
    Q_SA::T
    SHAT_SA::T
    ALPHA_SA::T
    XWELL_SA::T
    THETA0_SA::T

    DAMP_PSI::T
    DAMP_SIG::T
    WDIA_TRAPPED::T
    PARK::T
    GHAT::T
    GCHAT::T
    WD_ZERO::T
    LINSKER_FACTOR::T
    GRADB_FACTOR::T
    FILTER::T
    THETA_TRAPPED::T
    NN_MAX_ERROR::T


    function InputTJLF{T}(ns::Int,nky::Int) where T<:Real
        new("",
        missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
        missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
        fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),
        fill(NaN,(nky)),fill(NaN,(nky)),fill(NaN,(nky)),
        0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
        NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
        NaN,NaN,NaN,NaN,NaN,NaN,)
    end

    # create InputTJLF struct given a InputTGLF struct
    function InputTJLF{T}(inputTGLF::InputTGLF) where T<:Real

        inputTJLF = InputTJLF{T}(inputTGLF.NS,inputTGLF.NWIDTH)

        for fieldname in fieldnames(inputTGLF)
            if occursin(r"\d",String(fieldname)) || fieldname==:_Qgb # species parameter
                continue
            end
            setfield!(inputTJLF,fieldname,getfield(inputTGLF,fieldname))
        end
        for i in 1:inputTGLF.NS
            inputTJLF.ZS[i] = getfield(inputTGLF,Symbol("ZS_",i))
            inputTJLF.AS[i] = getfield(inputTGLF,Symbol("AS_",i))
            inputTJLF.MASS[i] = getfield(inputTGLF,Symbol("MASS_",i))
            inputTJLF.RLNS[i] = getfield(inputTGLF,Symbol("RLNS_",i))
            inputTJLF.RLTS[i] = getfield(inputTGLF,Symbol("RLTS_",i))
            inputTJLF.TAUS[i] = getfield(inputTGLF,Symbol("TAUS_",i))
            inputTJLF.VPAR[i] = getfield(inputTGLF,Symbol("VPAR_",i))
            inputTJLF.VPAR_SHEAR[i] = getfield(inputTGLF,Symbol("VPAR_SHEAR_",i))
            inputTJLF.VNS_SHEAR[i] = getfield(inputTGLF,Symbol("VNS_SHEAR_",i))
            inputTJLF.VTS_SHEAR[i] = getfield(inputTGLF,Symbol("VTS_SHEAR_",i))
        end
        inputTJLF.WIDTH_SPECTRUM .= inputTJLF.WIDTH

        field_names = fieldnames(InputTJLF)
        for field_name in field_names
            field_value = getfield(inputTJLF, field_name)
            if typeof(field_value)<:Real
                @assert !isnan(field_value) && !ismissing(field_value) "Did not properly populate inputTJLF for $field_name"
            end
            if typeof(field_value)<:Vector && field_name!=:KY_SPECTRUM && field_name!=:GAMMA_SPECTRUM
                for val in field_value
                    @assert !isnan(val) "Did not properly populate inputTJLF for array $field_name"
                end
            end
        end
    
        return inputTJLF
    end
end


##########################################################
#       Get from tjlf_hermite
##########################################################

struct OutputHermite{T<:Real}

    x::Vector{T}
    wx::Vector{T}
    h::Matrix{T}
    _dvec::Vector{T}
end

function OutputHermite(x, wx, h)
    _dvec = similar(wx)
    return OutputHermite(x, wx, h, _dvec)
end

##########################################################
#       Get from xgrid_functions_geo
##########################################################

struct OutputGeometry{T<:Real}

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

struct SaturationParameters{T<:Real}

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
    kpar_eff::Array{ComplexF64,3}
    modkpar_eff::Array{ComplexF64,3}

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
        kpar_eff = zeros(ComplexF64, ns, nbasis, nbasis)
        modkpar_eff = zeros(ComplexF64, ns, nbasis, nbasis)


        new(kx,
        wdh,modwdh,wdg,modwdg,
        gradB,b0,b0inv,lnB,p0,p0inv,bp,bpinv,
        c_par_par,c_tor_par,c_tor_per,
        kpar,modkpar,kpar_eff,modkpar_eff)
    end
end


mutable struct AveH{T<:Real}

    # average h-bessel functions
    hn::Array{T,3}
    hp1::Array{T,3}
    hp3::Array{T,3}
    hr11::Array{T,3}
    hr13::Array{T,3}
    hr33::Array{T,3}
    hw113::Array{T,3}
    hw133::Array{T,3}
    hw333::Array{T,3}

    ht1::Array{T,3}
    ht3::Array{T,3}
    hu1::Array{T,3}
    hu3::Array{T,3}
    hu33::Array{T,3}
    hu3ht1::Array{T,3}
    hu33ht1::Array{T,3}
    hu3ht3::Array{T,3}
    hu33ht3::Array{T,3}

    hnp0::Array{T,3}
    hp1p0::Array{T,3}
    hp3p0::Array{T,3}
    hr11p0::Array{T,3}
    hr13p0::Array{T,3}
    hr33p0::Array{T,3}
    c_tor_par_hp1p0::Array{T,3}
    c_tor_par_hr11p0::Array{T,3}
    c_tor_par_hr13p0::Array{T,3}

    hnb0::Array{T,3}
    hp1b0::Array{T,3}
    hp3b0::Array{T,3}
    hr11b0::Array{T,3}
    hr13b0::Array{T,3}
    hr33b0::Array{T,3}
    hw113b0::Array{T,3}
    hw133b0::Array{T,3}
    hw333b0::Array{T,3}

    hnbp::Array{T,3}
    hp1bp::Array{T,3}
    hp3bp::Array{T,3}
    hr11bp::Array{T,3}
    hr13bp::Array{T,3}
    hr33bp::Array{T,3}
    hw113bp::Array{T,3}
    hw133bp::Array{T,3}
    hw333bp::Array{T,3}

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

    wdhp1p0::Array{T,3}
    wdhr11p0::Array{T,3}
    wdhr13p0::Array{T,3}
    wdht1::Array{T,3}
    wdht3::Array{T,3}
    wdhu1::Array{T,3}
    wdhu3::Array{T,3}
    wdhu3ht1::Array{T,3}
    wdhu3ht3::Array{T,3}
    wdhu33::Array{T,3}
    wdhu33ht1::Array{T,3}
    wdhu33ht3::Array{T,3}

    modwdht1::Array{T,3}
    modwdht3::Array{T,3}
    modwdhu1::Array{T,3}
    modwdhu3::Array{T,3}
    modwdhu3ht1::Array{T,3}
    modwdhu3ht3::Array{T,3}
    modwdhu33::Array{T,3}
    modwdhu33ht1::Array{T,3}
    modwdhu33ht3::Array{T,3}

    wdhp1b0::Array{T,3}
    wdhr11b0::Array{T,3}
    wdhr13b0::Array{T,3}

    wdhp1bp::Array{T,3}
    wdhr11bp::Array{T,3}
    wdhr13bp::Array{T,3}

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

    kparhnp0::Array{ComplexF64,3}
    kparhp1p0::Array{ComplexF64,3}
    kparhp3p0::Array{ComplexF64,3}
    kparhu1::Array{ComplexF64,3}
    kparhu3::Array{ComplexF64,3}
    kparht1::Array{ComplexF64,3}
    kparht3::Array{ComplexF64,3}
    modkparhu1::Array{ComplexF64,3}
    modkparhu3::Array{ComplexF64,3}

    kparhp1b0::Array{ComplexF64,3}
    kparhr11b0::Array{ComplexF64,3}
    kparhr13b0::Array{ComplexF64,3}

    kparhnbp::Array{ComplexF64,3}
    kparhp3bp::Array{ComplexF64,3}
    kparhp1bp::Array{ComplexF64,3}
    kparhr11bp::Array{ComplexF64,3}
    kparhr13bp::Array{ComplexF64,3}

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
    gn::Array{T,3}
    gp1::Array{T,3}
    gp3::Array{T,3}
    gr11::Array{T,3}
    gr13::Array{T,3}
    gr33::Array{T,3}
    gw113::Array{T,3}
    gw133::Array{T,3}
    gw333::Array{T,3}

    gt1::Array{T,3}
    gt3::Array{T,3}
    gu1::Array{T,3}
    gu3::Array{T,3}
    gu33::Array{T,3}
    gu3gt1::Array{T,3}
    gu3gt3::Array{T,3}
    gu33gt1::Array{T,3}
    gu33gt3::Array{T,3}

    gnp0::Array{T,3}
    gp1p0::Array{T,3}
    gp3p0::Array{T,3}
    gr11p0::Array{T,3}
    gr13p0::Array{T,3}
    gr33p0::Array{T,3}
    c_tor_par_gp1p0::Array{T,3}
    c_tor_par_gr11p0::Array{T,3}
    c_tor_par_gr13p0::Array{T,3}

    gnb0::Array{T,3}
    gp1b0::Array{T,3}
    gp3b0::Array{T,3}
    gr11b0::Array{T,3}
    gr13b0::Array{T,3}
    gr33b0::Array{T,3}
    gw113b0::Array{T,3}
    gw133b0::Array{T,3}
    gw333b0::Array{T,3}

    gnbp::Array{T,3}
    gp1bp::Array{T,3}
    gp3bp::Array{T,3}
    gr11bp::Array{T,3}
    gr13bp::Array{T,3}
    gr33bp::Array{T,3}
    gw113bp::Array{T,3}
    gw133bp::Array{T,3}
    gw333bp::Array{T,3}

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

    wdgp1p0::Array{T,3}
    wdgr11p0::Array{T,3}
    wdgr13p0::Array{T,3}
    wdgu1::Array{T,3}
    wdgu3::Array{T,3}
    wdgu33::Array{T,3}
    wdgt1::Array{T,3}
    wdgt3::Array{T,3}
    wdgu3gt1::Array{T,3}
    wdgu3gt3::Array{T,3}
    wdgu33gt1::Array{T,3}
    wdgu33gt3::Array{T,3}

    modwdgu1::Array{T,3}
    modwdgu3::Array{T,3}
    modwdgu33::Array{T,3}
    modwdgt1::Array{T,3}
    modwdgt3::Array{T,3}
    modwdgu3gt1::Array{T,3}
    modwdgu3gt3::Array{T,3}
    modwdgu33gt1::Array{T,3}
    modwdgu33gt3::Array{T,3}

    wdgp1b0::Array{T,3}
    wdgr11b0::Array{T,3}
    wdgr13b0::Array{T,3}

    wdgp1bp::Array{T,3}
    wdgr11bp::Array{T,3}
    wdgr13bp::Array{T,3}

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

    kpargnp0::Array{ComplexF64,3}
    kpargp1p0::Array{ComplexF64,3}
    kpargp3p0::Array{ComplexF64,3}
    kpargu1::Array{ComplexF64,3}
    kpargu3::Array{ComplexF64,3}
    kpargt1::Array{ComplexF64,3}
    kpargt3::Array{ComplexF64,3}
    modkpargu1::Array{ComplexF64,3}
    modkpargu3::Array{ComplexF64,3}

    kpargp1b0::Array{ComplexF64,3}
    kpargr11b0::Array{ComplexF64,3}
    kpargr13b0::Array{ComplexF64,3}

    kpargnbp::Array{ComplexF64,3}
    kpargp3bp::Array{ComplexF64,3}
    kpargp1bp::Array{ComplexF64,3}
    kpargr11bp::Array{ComplexF64,3}
    kpargr13bp::Array{ComplexF64,3}

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