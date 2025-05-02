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
    USE_BISECTION::Bool = true
    USE_INBOARD_DETRAPPED::Bool = false
    NEW_EIKONAL::Bool = true
    FIND_WIDTH::Bool = true
    IFLUX::Bool = true
    ADIABATIC_ELEC::Bool = false

    NWIDTH::Int = 21
    NXGRID::Int = 16
    VPAR_MODEL::Int = 0
    VPAR_SHEAR_MODEL::Int = 1
    IBRANCH::Int = -1

    KY::Float64 = 0.3
    ALPHA_E::Float64 = 1.0
    ALPHA_P::Float64 = 1.0
    XNU_FACTOR::Float64 = 1.0
    DEBYE_FACTOR::Float64 = 1.0
    RLNP_CUTOFF::Float64 = 18.0
    WIDTH::Float64 = 1.65
    WIDTH_MIN::Float64 = 0.3
    BETA_LOC::Float64 = 1.0
    KX0_LOC::Float64 = 1.0
    PARK::Float64 = 1.0
    GHAT::Float64 = 1.0
    GCHAT::Float64 = 1.0
    WD_ZERO::Float64 = 0.1
    LINSKER_FACTOR::Float64 = 0.0
    GRADB_FACTOR::Float64 = 0.0
    FILTER::Float64 = 2.0
    THETA_TRAPPED::Float64 = 0.7
    ETG_FACTOR::Float64 = 1.25
    DAMP_PSI::Float64 = 0.0
    DAMP_SIG::Float64 = 0.0

end










mutable struct InputTJLF{T<:Real}

    UNITS::Union{String,Missing}

    USE_BPER::Union{Bool,Missing}
    USE_BPAR::Union{Bool,Missing}
    USE_MHD_RULE::Union{Bool,Missing}
    USE_BISECTION::Union{Bool,Missing}
    USE_INBOARD_DETRAPPED::Union{Bool,Missing}
    USE_AVE_ION_GRID::Union{Bool,Missing}
    NEW_EIKONAL::Union{Bool,Missing}
    FIND_WIDTH::Union{Bool,Missing}
    IFLUX::Union{Bool,Missing}
    ADIABATIC_ELEC::Union{Bool,Missing}

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
    IBRANCH::Union{Int,Missing}

    ZS::Union{Vector{T},Missing}
    MASS::Union{Vector{T},Missing}
    RLNS::Union{Vector{T},Missing}
    RLTS::Union{Vector{T},Missing}
    TAUS::Union{Vector{T},Missing}
    AS::Union{Vector{T},Missing}
    VPAR::Union{Vector{T},Missing}
    VPAR_SHEAR::Union{Vector{T},Missing}

    # NOT IN TGLF
    WIDTH_SPECTRUM::Union{Vector{T},Missing}
    KY_SPECTRUM::Union{Vector{T},Missing}
    EIGEN_SPECTRUM::Union{Vector{ComplexF64},Missing}
    FIND_EIGEN::Union{Bool,Missing}
    # NOT IN TGLF

    SIGN_BT::Union{Int,Missing}
    SIGN_IT::Union{Int,Missing}
    KY::Union{T,Missing}

    VEXB_SHEAR::Union{T,Missing}
    BETAE::Union{T,Missing}
    XNUE::Union{T,Missing}
    ZEFF::Union{T,Missing}
    DEBYE::Union{T,Missing}

    ALPHA_MACH::Union{T,Missing}
    ALPHA_E::Union{T,Missing}
    ALPHA_P::Union{T,Missing}
    ALPHA_QUENCH::Union{Int,Missing}
    ALPHA_ZF::Union{T,Missing}
    XNU_FACTOR::Union{T,Missing}
    DEBYE_FACTOR::Union{T,Missing}
    ETG_FACTOR::Union{T,Missing}
    RLNP_CUTOFF::Union{T,Missing}

    WIDTH::Union{T,Missing}
    WIDTH_MIN::Union{T,Missing}

    RMIN_LOC::Union{T,Missing}
    RMAJ_LOC::Union{T,Missing}
    ZMAJ_LOC::Union{T,Missing}
    DRMINDX_LOC::Union{T,Missing}
    DRMAJDX_LOC::Union{T,Missing}
    DZMAJDX_LOC::Union{T,Missing}
    Q_LOC::Union{T,Missing}
    KAPPA_LOC::Union{T,Missing}
    S_KAPPA_LOC::Union{T,Missing}
    DELTA_LOC::Union{T,Missing}
    S_DELTA_LOC::Union{T,Missing}
    ZETA_LOC::Union{T,Missing}
    S_ZETA_LOC::Union{T,Missing}
    P_PRIME_LOC::Union{T,Missing}
    Q_PRIME_LOC::Union{T,Missing}
    BETA_LOC::Union{T,Missing}
    KX0_LOC::Union{T,Missing}

    DAMP_PSI::Union{T,Missing}
    DAMP_SIG::Union{T,Missing}
    WDIA_TRAPPED::Union{T,Missing}
    PARK::Union{T,Missing}
    GHAT::Union{T,Missing}
    GCHAT::Union{T,Missing}
    WD_ZERO::Union{T,Missing}
    LINSKER_FACTOR::Union{T,Missing}
    GRADB_FACTOR::Union{T,Missing}
    FILTER::Union{T,Missing}
    THETA_TRAPPED::Union{T,Missing}
    SMALL::Union{T,Missing}

    #MXH params
    SHAPE_COS0::Union{T,Missing}
    SHAPE_COS1::Union{T,Missing}
    SHAPE_COS2::Union{T,Missing}
    SHAPE_COS3::Union{T,Missing}
    SHAPE_COS4::Union{T,Missing}
    SHAPE_COS5::Union{T,Missing}
    SHAPE_COS6::Union{T,Missing}

    SHAPE_SIN3::Union{T,Missing}
    SHAPE_SIN4::Union{T,Missing}
    SHAPE_SIN5::Union{T,Missing}
    SHAPE_SIN6::Union{T,Missing}

    SHAPE_S_COS0::Union{T,Missing}
    SHAPE_S_COS1::Union{T,Missing}
    SHAPE_S_COS2::Union{T,Missing}
    SHAPE_S_COS3::Union{T,Missing}
    SHAPE_S_COS4::Union{T,Missing}
    SHAPE_S_COS5::Union{T,Missing}
    SHAPE_S_COS6::Union{T,Missing}

    SHAPE_S_SIN3::Union{T,Missing}
    SHAPE_S_SIN4::Union{T,Missing}
    SHAPE_S_SIN5::Union{T,Missing}
    SHAPE_S_SIN6::Union{T,Missing}


    function InputTJLF()
        return InputTJLF{Float64}()
    end
    function InputTJLF{T}() where {T<:Real}
        new(
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing)
    end

    #For list-format inputs:
    function InputTJLF{T}(inP::Bool) where {T<:Real}
        if inP
            new("CGYRO", false, false, false, true, false, true, true, true, true, false, 2, 3, 5, 21, 6, 4, 16, 12, 4, 3, 0, -1, [-1.0, 1.0, 6.0], [0.0002723125672605524, 1.0, 6.0], [0.9691383387573976, 1.078021414201318, 0.0733427635614379], [3.332037619158914, 2.0626412607995435, 2.0626412607995435], [1.0, 1.3661261082028286, 1.3661261082028286], [1.0, 0.8075398023805694, 0.030988644410732645], [0.30611236015079274, 0.30611236015079274, 0.30611236015079274], [1.5491649356389778, 1.5491649356389778, 1.5491649356389778], [1.65, 1.65, 1.65, 1.65, 0.9685467847385054, 0.7035623639735143, 0.6324554384084142, 0.591251466897806, 0.5250292077902518, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.551199191728457, 1.5603157485179238, 1.5643992899672103, 1.5643992899672103, 1.4832394778484315, 0.6324554384084142, 0.5390399200630244, 0.4845607815598427], [0.05994688615887238, 0.11989377231774476, 0.17984065847661712, 0.2397875446354895, 0.2997344307943619, 0.5395219754298515, 0.6594157477475961, 0.7793095200653409, 0.8992032923830857, 1.0190970647008304, 1.138990837018575, 1.25888460933632, 1.1989377231774476, 1.60763400555729, 2.1556474918269477, 2.890468908319081, 3.875777714879744, 5.196960552637156, 6.968510831252536, 9.343950702231606, 12.529135254288084, 16.80009186935815, 22.52694069387401, 30.205969167995068], ComplexF64[0.01903754432811539 - 0.03822460700618263im, 0.066744785730153 - 0.08638900959186772im, 0.12700366682079575 - 0.13466694388477374im, 0.17881630543473276 - 0.16774618757919235im, 0.21216108783162443 - 0.18009019917758487im, 0.33316817373508417 - 0.3198305942745321im, 0.34745089281285046 - 0.3955067169459169im, 0.3399817917167648 - 0.4616241685193869im, 0.3137572066001085 - 0.5077973469905399im, 0.28555930064941276 + 0.45333962118452586im, 0.34112248363191056 + 0.501904780504562im, 0.3914434522233301 + 0.5529917255055463im, 0.36420457375571014 + 0.5235784686127289im, 0.531883979920564 + 0.7034735644682284im, 0.7500443094115281 + 0.9612053548664551im, 1.0612410507865202 + 1.3223708279218578im, 1.2154632046734721 + 1.7048988962925355im, 1.8556234007580341 + 2.2437573929973844im, 2.891191856550966 + 3.1140102214419576im, 4.213253106128079 + 4.483889088628243im, 5.747152586787458 + 6.01284811393025im, 7.352163708936412 + 7.770351222625782im, 9.652984589008142 + 10.27500211396812im, 11.970055852084016 + 14.27560461246188im], true, -1, 1, 0.3, 0.148365431821359, 0.0009809454014984833, 0.2658337070903717, 1.9296593, 0.029821537734289975, 0.0, 1.0, 1.0, 0, -1.0, 1.0, 1.0, 1.25, 18.0, 1.65, 0.3, 0.8896452200962354, 2.8058920740841784, 0.0, 1.0, -0.19752155788650919, 0.0, 3.3106313319155714, 1.6054967596315595, 0.39307251195418547, 0.21740011375976812, 0.7746695322421236, -0.05113765116526302, -0.2388377806334241, -0.0011605489188390146, 35.13388509054382, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.0, 0.0, 2.0, 0.7, 1.0e-12,  0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0)        else
        end
    end

    function InputTJLF{T}(ns::Int, nky::Int) where {T<:Real}
        new("",
        missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
        missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
        fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),
        fill(NaN,(nky)),fill(NaN,(nky)),fill(NaN*im,(nky)),missing,
        0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
        NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,1.0e-13,  0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0)
    end

    # create InputTJLF struct given a InputTGLF struct
    function InputTJLF{T}(inputTGLF::InputTGLF) where {T<:Real}

        nky = get_ky_spectrum_size(inputTGLF.NKY,inputTGLF.KYGRID_MODEL)
        inputTJLF = InputTJLF{T}(inputTGLF.NS, nky)

        for fieldname in fieldnames(inputTGLF)
            if occursin(r"\d", String(fieldname)) || fieldname == :_Qgb # species parameter
                continue
            end
            setfield!(inputTJLF, fieldname, getfield(inputTGLF, fieldname))
        end
        for i in 1:inputTGLF.NS
            inputTJLF.ZS[i] = getfield(inputTGLF, Symbol("ZS_", i))
            inputTJLF.AS[i] = getfield(inputTGLF, Symbol("AS_", i))
            inputTJLF.MASS[i] = getfield(inputTGLF, Symbol("MASS_", i))
            inputTJLF.RLNS[i] = getfield(inputTGLF, Symbol("RLNS_", i))
            inputTJLF.RLTS[i] = getfield(inputTGLF, Symbol("RLTS_", i))
            inputTJLF.TAUS[i] = getfield(inputTGLF, Symbol("TAUS_", i))
            inputTJLF.VPAR[i] = getfield(inputTGLF, Symbol("VPAR_", i))
            inputTJLF.VPAR_SHEAR[i] = getfield(inputTGLF, Symbol("VPAR_SHEAR_", i))
        end
        inputTJLF.WIDTH_SPECTRUM .= inputTJLF.WIDTH

        checkInput(inputTJLF)

        return inputTJLF
    end
end

function minimal_scalar_copy(inputs::TJLF.InputTJLF{T}) where T<:Real
    # Create a new instance using the constructor that sets up NS and the number of ky points.
    local_inputs = TJLF.InputTJLF{T}(inputs.NS, length(inputs.KY_SPECTRUM))
    # Loop over all fields and assign the value from the original.
    for f in fieldnames(typeof(inputs))
        setfield!(local_inputs, f, getfield(inputs, f))
    end
    return local_inputs
end

##########################################################
#       Get from tjlf_hermite
##########################################################

struct OutputHermite{T<:Real}
    x::Vector{T}
    wx::Vector{T}
    h::Matrix{T}
    _dvec::Matrix{T}
end

function OutputHermite(x, wx, h, nky::Int)
    _dvec = Matrix{typeof(wx[1])}(undef, size(wx, 1), nky)
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

    function Ave{T}(ns::Int, nbasis::Int) where {T<:Real}
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
            wdh, modwdh, wdg, modwdg,
            gradB, b0, b0inv, lnB, p0, p0inv, bp, bpinv,
            c_par_par, c_tor_par, c_tor_per,
            kpar, modkpar, kpar_eff, modkpar_eff)
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

    function AveH{T}(ns::Int, nbasis::Int) where {T<:Real}
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
            hn, hp1, hp3, hr11, hr13, hr33, hw113, hw133, hw333,
            ht1, ht3, hu1, hu3, hu33, hu3ht1, hu33ht1, hu3ht3, hu33ht3,
            hnp0, hp1p0, hp3p0, hr11p0, hr13p0, hr33p0, c_tor_par_hp1p0, c_tor_par_hr11p0, c_tor_par_hr13p0,
            hnb0, hp1b0, hp3b0, hr11b0, hr13b0, hr33b0, hw113b0, hw133b0, hw333b0,
            hnbp, hp1bp, hp3bp, hr11bp, hr13bp, hr33bp, hw113bp, hw133bp, hw333bp
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

    function AveWH{T}(ns::Int, nbasis::Int) where {T<:Real}
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
            wdhp1p0, wdhr11p0, wdhr13p0, wdht1, wdht3, wdhu1, wdhu3, wdhu3ht1, wdhu3ht3, wdhu33, wdhu33ht1, wdhu33ht3,
            modwdht1, modwdht3, modwdhu1, modwdhu3, modwdhu3ht1, modwdhu3ht3, modwdhu33, modwdhu33ht1, modwdhu33ht3,
            wdhp1b0, wdhr11b0, wdhr13b0,
            wdhp1bp, wdhr11bp, wdhr13bp
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
            kparhnp0, kparhp1p0, kparhp3p0, kparhu1, kparhu3, kparht1, kparht3, modkparhu1, modkparhu3,
            kparhp1b0, kparhr11b0, kparhr13b0,
            kparhnbp, kparhp3bp, kparhp1bp, kparhr11bp, kparhr13bp
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

    function AveG{T}(ns::Int, nbasis::Int) where {T<:Real}
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
            gn, gp1, gp3, gr11, gr13, gr33, gw113, gw133, gw333,
            gt1, gt3, gu1, gu3, gu33, gu3gt1, gu3gt3, gu33gt1, gu33gt3,
            gnp0, gp1p0, gp3p0, gr11p0, gr13p0, gr33p0, c_tor_par_gp1p0, c_tor_par_gr11p0, c_tor_par_gr13p0,
            gnb0, gp1b0, gp3b0, gr11b0, gr13b0, gr33b0, gw113b0, gw133b0, gw333b0,
            gnbp, gp1bp, gp3bp, gr11bp, gr13bp, gr33bp, gw113bp, gw133bp, gw333bp
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

    function AveWG{T}(ns::Int, nbasis::Int) where {T<:Real}
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
            wdgp1p0, wdgr11p0, wdgr13p0, wdgu1, wdgu3, wdgu33, wdgt1, wdgt3, wdgu3gt1, wdgu3gt3, wdgu33gt1, wdgu33gt3,
            modwdgu1, modwdgu3, modwdgu33, modwdgt1, modwdgt3, modwdgu3gt1, modwdgu3gt3, modwdgu33gt1, modwdgu33gt3,
            wdgp1b0, wdgr11b0, wdgr13b0,
            wdgp1bp, wdgr11bp, wdgr13bp
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
            kpargnp0, kpargp1p0, kpargp3p0, kpargu1, kpargu3, kpargt1, kpargt3, modkpargu1, modkpargu3,
            kpargp1b0, kpargr11b0, kpargr13b0,
            kpargnbp, kpargp3bp, kpargp1bp, kpargr11bp, kpargr13bp
        )
    end
end




mutable struct AveGrad{T<:Real}
    gradhp1::Array{T,3}
    gradhr11::Array{T,3}
    gradhr13::Array{T,3}
    gradhp1p1::Array{T,3}
    gradhr11p1::Array{T,3}
    gradhr13p1::Array{T,3}
    gradhp1p0::Array{T,3}
    gradhr11p0::Array{T,3}
    gradhr13p0::Array{T,3}

    gradgp1::Array{T,3}
    gradgr11::Array{T,3}
    gradgr13::Array{T,3}
    gradgp1p1::Array{T,3}
    gradgr11p1::Array{T,3}
    gradgr13p1::Array{T,3}
    gradgp1p0::Array{T,3}
    gradgr11p0::Array{T,3}
    gradgr13p0::Array{T,3}

    function AveGrad{T}(ns::Int, nbasis::Int) where {T<:Real}
        gradhp1 = zeros(T, ns, nbasis, nbasis)
        gradhr11 = zeros(T, ns, nbasis, nbasis)
        gradhr13 = zeros(T, ns, nbasis, nbasis)
        gradhp1p1 = zeros(T, ns, nbasis, nbasis)
        gradhr11p1 = zeros(T, ns, nbasis, nbasis)
        gradhr13p1 = zeros(T, ns, nbasis, nbasis)
        gradhp1p0 = zeros(T, ns, nbasis, nbasis)
        gradhr11p0 = zeros(T, ns, nbasis, nbasis)
        gradhr13p0 = zeros(T, ns, nbasis, nbasis)

        gradgp1 = zeros(T, ns, nbasis, nbasis)
        gradgr11 = zeros(T, ns, nbasis, nbasis)
        gradgr13 = zeros(T, ns, nbasis, nbasis)
        gradgp1p1 = zeros(T, ns, nbasis, nbasis)
        gradgr11p1 = zeros(T, ns, nbasis, nbasis)
        gradgr13p1 = zeros(T, ns, nbasis, nbasis)
        gradgp1p0 = zeros(T, ns, nbasis, nbasis)
        gradgr11p0 = zeros(T, ns, nbasis, nbasis)
        gradgr13p0 = zeros(T, ns, nbasis, nbasis)

        new(
            gradhp1, gradhr11, gradhr13, gradhp1p1, gradhr11p1, gradhr13p1, gradhp1p0, gradhr11p0, gradhr13p0,
            gradgp1, gradgr11, gradgr13, gradgp1p1, gradgr11p1, gradgr13p1, gradgp1p0, gradgr11p0, gradgr13p0
        )
    end
end

mutable struct AveGradB{T<:Real}
    gradBhp1::Array{T,3}
    gradBhp3::Array{T,3}
    gradBhr11::Array{T,3}
    gradBhr13::Array{T,3}
    gradBhr33::Array{T,3}
    gradBhu1::Array{T,3}
    gradBhu3::Array{T,3}
    gradBhu33::Array{T,3}

    gradBgp1::Array{T,3}
    gradBgp3::Array{T,3}
    gradBgr11::Array{T,3}
    gradBgr13::Array{T,3}
    gradBgr33::Array{T,3}
    gradBgu1::Array{T,3}
    gradBgu3::Array{T,3}
    gradBgu33::Array{T,3}

    function AveGradB{T}(ns::Int, nbasis::Int) where {T<:Real}
        gradBhp1 = zeros(T, ns, nbasis, nbasis)
        gradBhp3 = zeros(T, ns, nbasis, nbasis)
        gradBhr11 = zeros(T, ns, nbasis, nbasis)
        gradBhr13 = zeros(T, ns, nbasis, nbasis)
        gradBhr33 = zeros(T, ns, nbasis, nbasis)
        gradBhu1 = zeros(T, ns, nbasis, nbasis)
        gradBhu3 = zeros(T, ns, nbasis, nbasis)
        gradBhu33 = zeros(T, ns, nbasis, nbasis)

        gradBgp1 = zeros(T, ns, nbasis, nbasis)
        gradBgp3 = zeros(T, ns, nbasis, nbasis)
        gradBgr11 = zeros(T, ns, nbasis, nbasis)
        gradBgr13 = zeros(T, ns, nbasis, nbasis)
        gradBgr33 = zeros(T, ns, nbasis, nbasis)
        gradBgu1 = zeros(T, ns, nbasis, nbasis)
        gradBgu3 = zeros(T, ns, nbasis, nbasis)
        gradBgu33 = zeros(T, ns, nbasis, nbasis)

        new(
            gradBhp1, gradBhp3, gradBhr11, gradBhr13, gradBhr33, gradBhu1, gradBhu3, gradBhu33,
            gradBgp1, gradBgp3, gradBgr11, gradBgr13, gradBgr33, gradBgu1, gradBgu3, gradBgu33
        )
    end


#----------------------------------------------
struct ShiftAndInvert{TA,TB,TT}
    A_lu::TA
    B::TB
    temp::TT
end

function (M::ShiftAndInvert)(y, x)
    mul!(M.temp, M.B, x)
    ldiv!(y, M.A_lu, M.temp)
end


end
