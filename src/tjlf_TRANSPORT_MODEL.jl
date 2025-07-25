"""
    tjlf_TM(inputs::InputTJLF{T},satParams::SaturationParameters{T},outputHermite::OutputHermite{T}; return_both_eigenvalues::Bool=false) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl
    return_both_eigenvalues::Bool       - If true, returns both first and second pass eigenvalues (default: false)

outputs:
    QL_weights                          - 5d array of QL weights (field, species, mode, ky, type),
                                          type: (particle, energy, torodial stress, parallel stress, exchange)
    firstPass_eigenvalue                - 3d array of eigenvalues (mode, ky, type)
                                          type: (gamma, frequency)
    secondPass_eigenvalue               - 3d array of second pass eigenvalues (mode, ky, type) [only if return_both_eigenvalues=true]
                                          type: (gamma, frequency)

description:
    Main transport model function.
    Calls linear TGLF over a spectrum of ky's and computes spectral integrals of field, intensity, and fluxes.
    For backward compatibility, by default only returns first pass eigenvalues.
"""
function tjlf_TM(inputs::TJLF.InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}; return_both_eigenvalues::Bool=false) where T<:Real

    alpha_quench_in = inputs.ALPHA_QUENCH
    vexb_shear_s = inputs.VEXB_SHEAR * inputs.SIGN_IT
    ns = inputs.NS
    nmodes = inputs.NMODES
    ky_spect = inputs.KY_SPECTRUM
    nky = length(ky_spect)

    original_iflux = inputs.IFLUX

    # Output arrays
    firstPass_eigenvalue = zeros(Float64, nmodes, nky, 2)
    secondPass_eigenvalue = zeros(Float64, nmodes, nky, 2)  # Separate array for second pass results
    QL_weights = zeros(Float64, 3, ns, nmodes, nky, 5)

    if alpha_quench_in != 0.0 || vexb_shear_s == 0.0
        Threads.@threads for ky_index in eachindex(ky_spect)
            local_inputs = minimal_scalar_copy(inputs)
            onePass!(local_inputs, satParams, outputHermite, vexb_shear_s,
                     firstPass_eigenvalue, QL_weights, ky_index)
        end
        # No second pass needed for single-pass cases
        secondPass_eigenvalue .= firstPass_eigenvalue

    elseif !inputs.FIND_WIDTH
        inputs.IFLUX = false
        Threads.@threads for ky_index in eachindex(ky_spect)
            local_inputs = minimal_scalar_copy(inputs)
            widthPass!(local_inputs, satParams, outputHermite, firstPass_eigenvalue, ky_index)
        end
        # Initialize secondPass_eigenvalue with firstPass results before spectral shift
        secondPass_eigenvalue .= firstPass_eigenvalue
        
        # Calculate spectral shift and run second pass
        kx0_e = calculate_spectral_shift(inputs, satParams, firstPass_eigenvalue)
        inputs.IFLUX = true
        Threads.@threads for ky_index in eachindex(ky_spect)
            local_inputs = minimal_scalar_copy(inputs)
            secondPass!(local_inputs, satParams, outputHermite, kx0_e[ky_index],
                        firstPass_eigenvalue, secondPass_eigenvalue, QL_weights, ky_index)
        end

    else
        inputs.WIDTH_SPECTRUM .= inputs.WIDTH
        inputs.IFLUX = false
        Threads.@threads for ky_index in eachindex(ky_spect)
            local_inputs = minimal_scalar_copy(inputs)
            firstPass!(local_inputs, satParams, outputHermite, firstPass_eigenvalue, ky_index)
        end
        # Initialize secondPass_eigenvalue with firstPass results before spectral shift
        secondPass_eigenvalue .= firstPass_eigenvalue
        
        # Calculate spectral shift and run second pass
        kx0_e = calculate_spectral_shift(inputs, satParams, firstPass_eigenvalue)
        inputs.IFLUX = true
        Threads.@threads for ky_index in eachindex(ky_spect)
            local_inputs = minimal_scalar_copy(inputs)
            secondPass!(local_inputs, satParams, outputHermite, kx0_e[ky_index],
                        firstPass_eigenvalue, secondPass_eigenvalue, QL_weights, ky_index)
        end
    end

    inputs.IFLUX = original_iflux
    inputs.EIGEN_SPECTRUM .= firstPass_eigenvalue[1,:,1] .+ firstPass_eigenvalue[1,:,2] .* im

    if inputs.SAT_RULE == 0
        @debug "Using SAT0 means the return value of TM is not the QL weights, but actually flux_spectrum_out = QL_weights * phi_bar_out. Notice the difference near LS return statement."
    end

    # Return eigenvalue arrays based on flag
    if return_both_eigenvalues
        return QL_weights, firstPass_eigenvalue, secondPass_eigenvalue
    else
        return QL_weights, firstPass_eigenvalue
    end
end

"""
    calculate_spectral_shift(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, firstPass_eigenvalue::Array{T,3}) where T<:Real

Helper function to calculate spectral shift parameters. 
For SAT_RULE=2&3, calculates zonal mixing parameters from first pass eigenvalues and passes them to xgrid_functions_geo.
For SAT_RULE=1, uses the original behavior without zonal mixing parameters.

Returns:
    kx0_e - spectral shift array for all ky values
"""
function calculate_spectral_shift(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, firstPass_eigenvalue::Array{T,3}) where T<:Real
    # CRITICAL FIX: For SAT_RULE=2 and SAT_RULE=3, calculate zonal mixing parameters from first pass
    # and pass them to xgrid_functions_geo for proper spectral shift calculation
    if inputs.SAT_RULE == 2 || inputs.SAT_RULE == 3
        # Calculate zonal mixing parameters from first pass eigenvalues
        most_unstable_gamma_first_pass = firstPass_eigenvalue[1, :, 1]
        vzf_out_first_pass, kymax_out_first_pass, jmax_out_first_pass = get_zonal_mixing(inputs, satParams, most_unstable_gamma_first_pass)            
        # Use these parameters in spectral shift calculation
        kx0_e = xgrid_functions_geo(inputs, satParams, firstPass_eigenvalue[:,:,1]; 
                                    vzf_out_param=vzf_out_first_pass, 
                                    kymax_out_param=kymax_out_first_pass)
    else
        # For SAT_RULE=1, use the original behavior
        kx0_e = xgrid_functions_geo(inputs, satParams, firstPass_eigenvalue[:,:,1])
    end
    
    return kx0_e
end

#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

"""
    onePass!(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T},vexb_shear_s::T,eigenvalue_spectrum_out::Array{T,3}, QL_weights::Array{T,5}, ky_index::Int) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl
    vexb_shear_s::T                     - exb value, VEXB_SHEAR*SIGN_IT from input.tglf file
    eigenvalue_spectrum_out::Array{T,3} - eigenvalue output array (mode, ky, type)
    QL_weights::Array{T,5}              - QL weights output array (field, species, mode, ky, type)
    ky_index::Int                       - index used for multithreading

outputs (writen to, not outputed):
    QL_weights                          - 5d array of QL weights (field, species, mode, ky, type),
                                          type: (particle, energy, torodial stress, parallel stress, exchange)
    eigenvalue_spectrum_out             - 3d array of eigenvalues (mode, ky, type)
                                          type: (gamma, frequency)

description:
    calculate both eigenvalues and fluxes on firstPass if vexb_shear = 0 or alpha_quench != 0, no secondPass
"""
function onePass!(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T},
    vexb_shear_s::T,
    eigenvalue_spectrum_out::Array{T,3}, QL_weights::Array{T,5}, ky_index::Int) where T<:Real

    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    # increment through the ky_spectrum and find the width/eigenvalues of each ky
    ky = inputs.KY_SPECTRUM[ky_index]

    if(inputs.FIND_WIDTH)
        inputs.IFLUX = false
        if(inputs.NEW_EIKONAL)
            # println("this is 1")
            nmodes_out, gamma_nb_min_out,
            gamma_out, freq_out,
            particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out = tjlf_max(inputs, satParams, outputHermite, ky, vexb_shear_s, ky_index)
        else
            error("NOT IMPLEMENTED YET -DSUN")
        end
    else
        inputs.IFLUX = true
        # use new nbasis value
        nbasis = inputs.NBASIS_MAX
        gamma_nb_min_out = NaN
        # println("this is XII")
        nmodes_out, gamma_out, freq_out,
        particle_QL_out,energy_QL_out,stress_tor_QL_out,stress_par_QL_out,exchange_QL_out,
        ft_test = tjlf_LS(inputs, satParams, outputHermite, ky, nbasis, vexb_shear_s, ky_index)

        if(inputs.USE_INBOARD_DETRAPPED && inputs.IBRANCH==-1) # check for inward ballooning modes
            b_geo = satParams.B_geo
            Bmax,_ = findmax(b_geo)
            Bmin,_ = findmin(b_geo)
            modB_test = 0.5*(Bmax + Bmin)/Bmin
            if(ft_test > modB_test)
                ft_min = 0.01
                outputGeo = xgrid_functions_geo(inputs, satParams, outputHermite, ky, ky_index; kx0_e)
                outputGeo.fts .= ft_min
                @warn "NOT TESTED TM.jl ln 144"
                # println("this is XIII")
                nmodes_out, gamma_out, freq_out,
                particle_QL_out,energy_QL_out,stress_tor_QL_out,stress_par_QL_out,exchange_QL_out,
                _ = tjlf_LS(inputs, satParams, outputHermite, ky, nbasis, vexb_shear_s, ky_index;outputGeo)
            end
        end
    end

    unstable = true
    gamma_max = findmax(gamma_out)[1] # this covers ibranch=-1,0
    if(gamma_max == 0.0 || gamma_nb_min_out == 0.0) unstable = false end

    if(unstable)
        eigenvalue_spectrum_out[1:nmodes_out,ky_index,1] .= gamma_out[1:nmodes_out]
        eigenvalue_spectrum_out[1:nmodes_out,ky_index,2] .= freq_out[1:nmodes_out]
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,1] .= particle_QL_out[:,ns0:ns,1:nmodes_out]
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,2] .= energy_QL_out[:,ns0:ns,1:nmodes_out]
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,3] .= stress_tor_QL_out[:,ns0:ns,1:nmodes_out]
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,4] .= stress_par_QL_out[:,ns0:ns,1:nmodes_out]
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,5] .= exchange_QL_out[:,ns0:ns,1:nmodes_out]
    end

end

#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
"""
    firstPass(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}, firstPass_eigenvalue::Array{T,3}, ky_index::Int) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl
    firstPass_eigenvalue::Array{T,3}    - eigenvalue output array (mode, ky, type)
    ky_index::Int                       - index used for multithreading

outputs (writen to, not outputed):
    firstPass_eigenvalue                - 3d array of eigenvalues (mode, ky, type)
                                          type: (gamma, frequency)

description:
    calculate the widths and eigenvalues with vexb_shear = 0.0
"""
function firstPass!(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}, firstPass_eigenvalue::Array{T,3}, ky_index::Int) where T<:Real

    # increment through the ky_spectrum and find the width/eigenvalues of each k
    ky = inputs.KY_SPECTRUM[ky_index]

    if(inputs.NEW_EIKONAL) # not sure what this is -DSUN
        # println("this is 1")
        nmodes_out, gamma_nb_min_out,
        gamma_out, freq_out,
        _,_,_,_,_ = tjlf_max(inputs, satParams, outputHermite, ky, 0.0, ky_index)
    else
        error("NOT IMPLEMENTED YET -DSUN")
    end

    unstable = true
    gamma_max = findmax(gamma_out)[1] # this covers ibranch=-1,0
    if(gamma_max == 0.0 || gamma_nb_min_out == 0.0) unstable = false end

    if(unstable)
        firstPass_eigenvalue[1:nmodes_out,ky_index,1] .= gamma_out[1:nmodes_out]
        firstPass_eigenvalue[1:nmodes_out,ky_index,2] .= freq_out[1:nmodes_out]
    end

end

#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
"""
    widthPass(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}, firstPass_eigenvalue::Array{T,3}, ky_index::Int) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl
    firstPass_eigenvalue::Array{T,3}    - eigenvalue output array (mode, ky, type)
    ky_index::Int                       - index used for multithreading

outputs (writen to, not outputed):
    firstPass_eigenvalue                - 3d array of eigenvalues (mode, ky, type)
                                          type: (gamma, frequency)

description:
    calculate the eigenvalues to calculate the spectral shift (kx0_e) on secondPass
"""
function widthPass!(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}, firstPass_eigenvalue::Array{T,3}, ky_index::Int) where T<:Real

    if(!inputs.NEW_EIKONAL)
        error("not sure what this flag is -DSUN")
    end

    # increment through the ky/width spectrum
    ky = inputs.KY_SPECTRUM[ky_index]
    if inputs.WIDTH_SPECTRUM[ky_index]<0.0 # the abs() was an idea to set the width to be <0 if no growth rate was found
        return
    end

    # calculate the eigenvalues with no shear
    nbasis = inputs.NBASIS_MAX
    nmodes_out, gamma_out, freq_out,
    _,_,_,_,_,ft_test = tjlf_LS(inputs, satParams, outputHermite, ky, nbasis, 0.0,ky_index)

    if(inputs.IBRANCH==-1) # check for inward ballooning modes
        if(inputs.USE_INBOARD_DETRAPPED) ####### find ft_test and modB_test
            if(inputs.USE_INBOARD_DETRAPPED)
                b_geo = satParams.B_geo
                Bmax,_ = findmax(b_geo)
                Bmin,_ = findmin(b_geo)
                modB_test = 0.5*(Bmax + Bmin)/Bmin
                if(ft_test > modB_test)
                    ft_min = 0.01
                    outputGeo = xgrid_functions_geo(inputs, satParams, outputHermite, ky, ky_index; kx0_e)
                    outputGeo.fts .= ft_min
                    @warn "NOT TESTED max.jl ln 302"
                    # println("this is XIII")
                    nmodes_out, gamma_out, freq_out,
                    _,_,_,_,_,_ = tjlf_LS(inputs, satParams, outputHermite, ky, nbasis, vexb_shear_s, ky_index;outputGeo)
                end
            end
        end
    end

    unstable = true
    gamma_max = findmax(gamma_out)[1] # this covers ibranch=-1,0
    if(gamma_max == 0.0) unstable = false end ############################ deleted gammma_nb_min_out ############################

    if(unstable)
        firstPass_eigenvalue[1:nmodes_out,ky_index,1] .= gamma_out[1:nmodes_out]
        firstPass_eigenvalue[1:nmodes_out,ky_index,2] .= freq_out[1:nmodes_out]
    end

end
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

"""
    secondPass!(inputs::InputTJLF{T},satParams::SaturationParameters{T},outputHermite::OutputHermite{T},kx0_e::T,firstPass_eigenvalue::Array{T,3}, QL_weights::Array{T,5}, ky_index::Int) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl
    kx0_e::T                            - spectral shift calculated using first/width pass eigenvalues
    firstPass_eigenvalue::Array{T,3}    - array of eigenvalues found in first pass
    QL_weights::Array{T,5}              - QL weights output array (field, species, mode, ky, type)
    ky_index::Int                       - index used for multithreading

outputs (writen to, not outputed):
    QL_weights                          - 5d array of QL weights (field, species, mode, ky, type),
                                          type: (particle, energy, torodial stress, parallel stress, exchange)

description:
    calculate the QL weights using the width and eigenvalues from the first pass as reference
"""
function secondPass!(inputs::InputTJLF{T}, satParams::SaturationParameters{T},outputHermite::OutputHermite{T},
    kx0_e::T,
    firstPass_eigenvalue::Array{T,3}, secondPass_eigenvalue::Array{T,3}, QL_weights::Array{T,5}, ky_index::Int) where T<:Real

    # input values
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    nbasis_max_in = inputs.NBASIS_MAX
    nbasis_min_in = inputs.NBASIS_MIN
    vexb_shear_s = inputs.VEXB_SHEAR*inputs.SIGN_IT

    # saturation values
    R_unit = satParams.R_unit

    # initialize reference eigenvalue arrays
    gamma_reference_kx0 = similar(firstPass_eigenvalue[:,1,1])
    freq_reference_kx0 = similar(firstPass_eigenvalue[:,1,2])

    ky = inputs.KY_SPECTRUM[ky_index]

    if(inputs.NEW_EIKONAL)
        gamma_reference_kx0 .= firstPass_eigenvalue[:,ky_index,1]
        freq_reference_kx0 .= firstPass_eigenvalue[:,ky_index,2]
        nbasis = nbasis_max_in

        nmodes_out, gamma_out, freq_out,
        particle_QL_out,energy_QL_out,stress_tor_QL_out,stress_par_QL_out,exchange_QL_out,
        _ = tjlf_LS(inputs, satParams, outputHermite, ky, nbasis, vexb_shear_s, ky_index;
                            kx0_e, gamma_reference_kx0, freq_reference_kx0)

        gamma_nb_min_out = gamma_out[1]
    else
        error("NOT IMPLEMENTED YET -DSUN")
    end

    unstable = true
    if inputs.IBRANCH==0
        gamma_max = max(gamma_out[1],gamma_out[2])
    else
        gamma_max = max(gamma_out[1],0.0)
    end
    if(gamma_max == 0.0 || gamma_nb_min_out == 0.0) unstable = false end

    # CRITICAL FIX: Store eigenvalues returned by tjlf_LS (which are first pass eigenvalues in spectral shift mode)
    # This matches TGLF behavior where second pass returns first pass eigenvalues for spectral shift
    if nmodes_out > 0
        secondPass_eigenvalue[1:nmodes_out,ky_index,1] .= gamma_out[1:nmodes_out]
        secondPass_eigenvalue[1:nmodes_out,ky_index,2] .= freq_out[1:nmodes_out]
    else
        # When eigensolver completely fails (nmodes_out = 0), store zeros in first mode
        secondPass_eigenvalue[1,ky_index,1] = 0.0
        secondPass_eigenvalue[1,ky_index,2] = 0.0
    end

    if(unstable)
        # println("DSUN")
        gamma_cutoff = 0.1*ky/R_unit
        rexp = 1.0
        reduce = 1.0
        if(nbasis_max_in != nbasis_min_in)
            if(gamma_nb_min_out < gamma_out[1] && gamma_nb_min_out < gamma_cutoff)
                reduce = (gamma_nb_min_out/gamma_cutoff)^rexp
            end
        end

        @views QL_weights[:,ns0:ns,1:nmodes_out,ky_index,1] .= particle_QL_out[:,ns0:ns,1:nmodes_out]
        @views QL_weights[:,ns0:ns,1:nmodes_out,ky_index,2] .= energy_QL_out[:,ns0:ns,1:nmodes_out]
        @views QL_weights[:,ns0:ns,1:nmodes_out,ky_index,3] .= stress_tor_QL_out[:,ns0:ns,1:nmodes_out]
        @views QL_weights[:,ns0:ns,1:nmodes_out,ky_index,4] .= stress_par_QL_out[:,ns0:ns,1:nmodes_out]
        @views QL_weights[:,ns0:ns,1:nmodes_out,ky_index,5] .= exchange_QL_out[:,ns0:ns,1:nmodes_out]

    end

    # recompute spectrum using non-local in ky multiscale model
    # if(sat_rule_in >= 1) get_multiscale_spectrum(inputs) end
    # if(new_eikonal_in) eikonal_unsaved = false end
    # gamma_out[1] = gmax
    # freq_out[1] = fmax
    # inputs.NEW_EIKONAL = true  # reset default for next call to tjlf_TM
end