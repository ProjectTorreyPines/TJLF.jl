# I need to create a struct for the TGLF input so I can then call it
# from TJLF.
# I also need a struct for the inputs .profile and .TGLFEP.

# In order to call tjlf_read_input, I actually am going to need to 
# produce an input.tglf file for a given input .profile and .TGLFEP.
# Then I can directly call tjlf_read_input function on the new file and
# I won't even have to transfer any of that. 

# I could also try to devise a way to transfer the information
# directly from the gacode and TGLFEP files.

# Lastly, I could also try to follow more in the exact footsteps
# of tglf-ep, but I honestly think that will take more time in execution
# of the code.

# I need to determine what needs to go into the tglf file first and then
# print the inputs that I get from either the TGLFEP or gacode (or profile)
# 

# Would it be easier to try to transfer direct information into an InputTJLF
# struct instead? I would need to make a new readGAInput function or 
# readEPInput or something. Then it could be called as in TJLF.
# The issue I'm thinking of with trying to write a new input.tglf is that
# I'm going to have to scan through the gacode and/or TGLFEP files anyways
# so why not directly translate them into TJLF? That is what it seems
# TJLF is doing from the input.tglf files. It is a direct translation
# and wouldn't require any extra calls to TJLF, even if it wouldn't be too
# much of a pain.

# In that case, I do want to create both a read__Input function and a TJLF 
# input struct:


mutable struct InputTJLF{T<:Real}

    UNITS::Union{String,Missing}

    USE_BPER::Union{Bool,Missing}
    USE_BPAR::Union{Bool,Missing}
    USE_MHD_RULE::Union{Bool,Missing}
    USE_BISECTION::Union{Bool,Missing} # 5
    USE_INBOARD_DETRAPPED::Union{Bool,Missing}
    USE_AVE_ION_GRID::Union{Bool,Missing} # don't see this in TGLF
    NEW_EIKONAL::Union{Bool,Missing} ## this seems useless, the flag has to be both true and false to do anything
    FIND_WIDTH::Union{Bool,Missing}
    IFLUX::Union{Bool,Missing} # 10
    ADIABATIC_ELEC::Union{Bool,Missing}

    SAT_RULE::Union{Int,Missing}
    NS::Union{Int,Missing}
    NMODES::Union{Int,Missing}
    NWIDTH::Union{Int,Missing} # 15
    NBASIS_MAX::Union{Int,Missing}
    NBASIS_MIN::Union{Int,Missing}
    NXGRID::Union{Int,Missing}
    NKY::Union{Int,Missing}
    KYGRID_MODEL::Union{Int,Missing} # 20
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

    SIGN_BT::Union{T,Missing}
    SIGN_IT::Union{T,Missing}
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
    RLNP_CUTOFF::Union{T,Missing} # not in tglf documentation

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
    BETA_LOC::Union{T,Missing} # not in tglf documentation
    KX0_LOC::Union{T,Missing}

    DAMP_PSI::Union{T,Missing} # not in tglf documentation
    DAMP_SIG::Union{T,Missing} # not in tglf documentation
    WDIA_TRAPPED::Union{T,Missing} # not in tglf documentation
    PARK::Union{T,Missing}
    GHAT::Union{T,Missing}
    GCHAT::Union{T,Missing}
    WD_ZERO::Union{T,Missing}
    LINSKER_FACTOR::Union{T,Missing}
    GRADB_FACTOR::Union{T,Missing}
    FILTER::Union{T,Missing}
    THETA_TRAPPED::Union{T,Missing}
    SMALL::Union{T,Missing}

    USE_TRANSPORT_MODEL::Union{Bool, Missing}
    
    # Any needed constructions of InputTJLF will be put here:
    
    # One of the TJLF methods is for an InputTGLF struct, but that isn't
    # used at all to my knowledge in TJLF, so I will not replicate it for
    # now. I will be creating the blank method, though:
    
    function InputTJLF{T}(ns::Int, nky::Int, dflt::Bool) where {T<:Real}
        if dflt
            new("GYRO",
            false,false,true,true,false,missing,true,false,true,false,
            0,2,2,21,4,4,32,12,0,2,0,-1,
            fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),
            fill(NaN,(nky)),fill(NaN,(nky)),fill(NaN*im,(nky)),missing,
            1.0,1.0,0.3,0.0,0.0,0.0,1.0,0.0,0.0,1.0,
            1.0,0.0,1.0,1.0,1.0,1.25,18.0,1.65,0.3,0.5,
            3.0,0.0,1.0,0.0,0.0,2.0,1.0,16.0,0.0,0.0,
            0.0,0.0,0.0,16.0,0.0,0.0,0.0,0.0,0.0,1.0,
            1.0,1.0,0.1,0.0,0.0,0.0,0.7,1.0e-13,true)
        else
            new("",
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,missing,
            fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),fill(NaN,(ns)),
            fill(NaN,(nky)),fill(NaN,(nky)),fill(NaN*im,(nky)),missing,
            0.0,0.0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
            NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,1.0e-13,true)
        end
    end

end

# When adjusting a module, you need to restart your REPL even though it is mutable.

mutable struct InputTJLFEP{T<:Real} # This acts as the interface module of Fortran, essentially. It reads the TGLFEP file

    PROCESS_IN::Union{Int, Missing} # May need to be a MODE_IN::Union{T, Missing} for PROCESS_IN <= 1

    THRESHOLD_FLAG::Union{Int, Missing}
    N_BASIS::Union{Int, Missing}
    SCAN_METHOD::Union{Int, Missing}

    REJECT_I_PINCH_FLAG::Union{Int, Missing} # 5
    REJECT_E_PINCH_FLAG::Union{Int, Missing}
    REJECT_TH_PINCH_FLAG::Union{Int, Missing}
    REJECT_EP_PINCH_FLAG::Union{Int, Missing}
    REJECT_TEARING_FLAG::Union{Int, Missing}
    REJECT_MAX_OUTER_PANEL_FLAG::Union{Int, Missing} # 10

    QL_THRESH_RATIO::Union{T, Missing}
    Q_SCALE::Union{T, Missing}

    WRITE_WAVEFUNCTION::Union{Int, Missing} # Maybe not needed for my purposes? Not sure yet
    KY_MODEL::Union{Int, Missing}

    SCAN_N::Union{Int} # 15
    IRS::Union{Int, Missing}

    FACTOR_IN_PROFILE::Union{Bool, Missing}
    FACTOR_IN::Union{T, Missing}

    WIDTH_IN_FLAG::Union{Bool}
    WIDTH_IN::Union{T, Missing} # 20
    WIDTH_MIN::Union{T, Missing}
    WIDTH_MAX::Union{T, Missing}

    INPUT_PROFILE_METHOD::Union{Int, Missing} # I very likely will not use this (?). It seems that INPUT_PROFILE_METHOD is meant for 
    # determining the type of input that will be used: input.profile or input.profiles. EXPRO is the method used to do the latter
    # but I cannot find the runction in Fortran for expro_read because it is imported from somewhere else. 
    
    #NION
    #

    IS_EP::Union{Int, Missing}
    M_I::Union{T, Missing} # 25
    Q_I::Union{T, Missing}
    M_EP::Union{T, Missing}
    Q_EP::Union{T, Missing}

    IR::Union{Int, Missing}
    NMODES::Union{Int} # 30 

    MODE_IN::Union{Int, Missing} # Likely will not be used for now. There are some default settings dependent on the choice of process-in. = 2 for process_in = 5: the case we are interested in.
    FACTOR_MAX_PROFILE::Union{Vector{T}, Missing}
    FACTOR_MAX::Union{T, Missing}
    KY_IN::Union{T, Missing}
    KYMARK::Union{T, Missing} # 35

    # All of the next 7 vectors require the value nr to be recognized for the array. This is already done in the profile struct.
    # I need to decide what I want to do with this then. These are for process_in = 6 apparently.
    DEP_TRACE_LOC::Union{T, Missing}
    DEP_TRACE::Union{Vector{T}, Missing}
    DEP_TRACE_COMPLETE::Union{Vector{T}, Missing}
    DEP_TRACE_SCALE::Union{Matrix{T}, Missing}

    QL_RATIO_LOC::Union{T, Missing} # 40
    QL_RATIO::Union{Vector{T}, Missing}
    QL_RATIO_COMPLETE::Union{Vector{T}, Missing}
    QL_RATIO_SCAN::Union{Matrix{T}, Missing}

    CHI_GB::Union{Vector{T}, Missing}
    IR_EXP::Union{Vector{Int64}, Missing} # 45 

    NBASIS::Union{Int, Missing}
    NTOROIDAL::Union{Int, Missing}
    KYHAT_IN::Union{T, Missing}

    SCAN_FACTOR::Union{T, Missing} # Default to 1
    JTSCALE::Union{Int, Missing} # 50
    JTSCALE_MAX::Union{Int, Missing} # Default to 1
    TSCALE_INTERVAL::Union{T, Missing} # Default to 1.0

    TGLFEP_NION::Union{Int, Missing} # Equal to NS-1?
    
    FREQ_CUTOFF::Union{T, Missing} # Default to -0.2
    FREQ_AE_UPPER::Union{T, Missing} # 55

    NN::Union{Int, Missing}
    FACTOR_OUT::Union{Vector{T}, Missing} # vector of nn elements
    
    # May or may not need id_2, np_2, id_3, np_3 ?

    # Ignoring suffix and str_wf_file right now.
    # As with l_print, l-WRITE_WAVEFUNCTION, l_wavefunction_out

    LKEEP::Vector{Bool} # All the L-starting parts have 4 dimensions. See kwscale_scan for that.
    LTEARING::Vector{Bool} # This needs to be defined somewhat differently then.
    L_TH_PINCH::Vector{Bool} # 60
    L_I_PINCH::Vector{Bool}
    L_E_PINCH::Vector{Bool}
    L_EP_PINCH::Vector{Bool}
    L_QL_RATIO::Vector{Bool}
    L_MAX_OUTER_PANEL::Vector{Bool} # 65

    FACTOR::Union{Vector{T}, Missing}
    SUFFIX::Union{String, Missing}
    F_REAL::Union{Vector{T}, Missing}
    REAL_FREQ::Int

    WIDTH::Union{Vector{T}, Missing}

    # I might want a default constructor for default values. I'm going to make one (no paramters)
    # which follows the 'default' case that TGLFEP used in my first test of TGLFEP (OMFIT).

    # There are a bunch of cases where the dimension of the fills or the existence of a value (WIDTH_IN and WIDTH_MIN/MAX)
    # depends on other input variables:

    # If WIDTH_IN_FLAG is true the fill needs (?) to be SCAN_N long of that same value and the min/max are missing
    # There is a bit of confusion on if the same goes for the FACTOR_IN_PROFILE. 

    # There are few things that depend on the profile reading: nr
    function InputTJLFEP{T}(nscan_in::Int, widthin::Bool, nn::Int, nr::Int, jtscale_max::Int, nmodes::Int) where {T<:Real}
        if(widthin)
            new(missing, missing, missing, missing, missing, missing, missing, missing, missing, missing,
            NaN, NaN, missing, missing, nscan_in, missing, missing, NaN, widthin, NaN,
            NaN, NaN, missing, missing, NaN, NaN, NaN, NaN, missing, nmodes, missing, fill(NaN, nscan_in),
            NaN, NaN, NaN, NaN, fill(NaN, nr), fill(NaN, nr), fill(NaN, (jtscale_max, nr)), NaN, fill(NaN, nr), fill(NaN, nr),
            fill(NaN, (jtscale_max, nr)), fill(NaN, nr), fill(0, nscan_in), missing, missing, NaN, 1, missing, 1, 
            1.0, missing, -0.2, NaN, nn, fill(NaN, nn), fill(false, 4), fill(false, 4),
            fill(false, 4), fill(false, 4), fill(false, 4),
            fill(false, 4), fill(false, 4), fill(false, 4),
            fill(NaN, nscan_in), missing, fill(NaN, nr), 0, fill(NaN, nscan_in))
        else
            new(missing, missing, missing, missing, missing, missing, missing, missing, missing, missing,
            NaN, NaN, missing, missing, nscan_in, missing, missing, NaN, widthin, 0.0,
            NaN, NaN, missing, missing, NaN, NaN, NaN, NaN, missing, nmodes, missing, fill(NaN, nscan_in),
            NaN, NaN, NaN, NaN, fill(NaN, nr), fill(NaN, nr), fill(NaN, (jtscale_max, nr)), NaN, fill(NaN, nr), fill(NaN, nr),
            fill(NaN, (jtscale_max, nr)), fill(NaN, nr), fill(0, nscan_in), missing, missing, NaN, 1, missing, 1, 
            1.0, missing, -0.2, NaN, nn, fill(NaN, nn), fill(false, 4), fill(false, 4),
            fill(false, 4), fill(false, 4), fill(false, 4),
            fill(false, 4), fill(false, 4), fill(false, 4),
            fill(NaN, nscan_in), missing, fill(NaN, nr), 0, fill(NaN, nscan_in))
        end
    end
end

mutable struct profile{T<:Real}
    # This struct is an intermediary between the processes in read_inputs
    
    SIGN_BT::Union{T, Missing} #
    SIGN_IT::Union{T, Missing} #

    NR::Union{Int, Missing} #
    NS::Union{Int, Missing} #
    GEOMETRY_FLAG::Union{Int, Missing} #
    ROTATION_FLAG::Union{Int, Missing} #

    ZS::Union{Vector{T}, Missing} #
    MASS::Union{Vector{T}, Missing} #

    # This next section is especially why the profile struct must exist
    AS::Union{Matrix{T}, Missing} #
    TAUS::Union{Matrix{T}, Missing} #
    RLNS::Union{Matrix{T}, Missing} #
    RLTS::Union{Matrix{T}, Missing} #
    VPAR::Union{Matrix{T}, Missing} #
    VPAR_SHEAR::Union{Matrix{T}, Missing} #

    RMIN::Union{Vector{T}, Missing}
    RMAJ::Union{Vector{T}, Missing}
    SHIFT::Union{Vector{T}, Missing}
    Q::Union{Vector{T}, Missing}
    SHEAR::Union{Vector{T}, Missing}
    ALPHA::Union{Vector{T}, Missing}
    Q_PRIME::Union{Vector{T}, Missing}
    P_PRIME::Union{Vector{T}, Missing}
    KAPPA::Union{Vector{T}, Missing}
    S_KAPPA::Union{Vector{T}, Missing}
    DELTA::Union{Vector{T}, Missing}
    S_DELTA::Union{Vector{T}, Missing}
    ZETA::Union{Vector{T}, Missing}
    S_ZETA::Union{Vector{T}, Missing}
    ZEFF::Union{Vector{T}, Missing}
    BETAE::Union{Vector{T}, Missing}
    RHO_STAR::Union{Vector{T}, Missing}
    OMEGA_TAE::Union{Vector{T}, Missing}
    B_UNIT::Union{Vector{T}, Missing}

    IS::Union{Int, Missing}
    IRS::Union{Int, Missing}

    A_QN::Union{T, Missing} # quasineutrality scale factor for non-EP species
    TGLFEP_NION::Union{Int, Missing} #

    # As of right now, I don't believe there needs to be parameters, but the vectors
    # are probably the most of concern there. 
    function profile{T}(nr::Int, ns::Int) where (T<:Real)
        new(NaN, NaN, missing, missing, missing, missing, fill(NaN, ns), fill(NaN, ns),
        fill(NaN, (nr, ns)), fill(NaN, (nr, ns)), fill(NaN, (nr, ns)), fill(NaN, (nr, ns)), 
        fill(NaN, (nr, ns)), fill(NaN, (nr, ns)), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), 
        fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), 
        fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), 
        fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), fill(NaN, nr), missing, missing, NaN, missing)
    end
end

