"""
    function tjlf_TM(inputs::InputTJLF{T},satParams::SaturationParameters{T},outputHermite::OutputHermite{T}) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl

outputs:
    QL_weights                          - 5d array of QL weights (field, species, mode, ky, type),
                                          type: (particle, energy, torodial stress, parallel stress, exchange)
    firstPass_eigenvalue                - 3d array of eigenvalues (mode, ky, type)
                                          type: (gamma, frequency)

description:
    Main transport model subroutine.
    Calls linear TGLF over a spectrum of ky's and computes spectral integrals of field, intensity and fluxes.
"""

function tjlf_TM(inputs::InputTJLF{T},satParams::SaturationParameters{T},outputHermite::OutputHermite{T}) where T<:Real

    alpha_quench_in = inputs.ALPHA_QUENCH
    vexb_shear_s = inputs.VEXB_SHEAR*inputs.SIGN_IT
    ns = inputs.NS
    nmodes = inputs.NMODES
    ky_spect = inputs.KY_SPECTRUM
    nky = length(ky_spect)
    
    original_iflux = inputs.IFLUX

    ### output values
    firstPass_eigenvalue = zeros(Float64, nmodes, nky, 2)
    QL_weights::Array{Float64,5} = zeros(Float64, 3, ns, nmodes, nky, 5)

    # compute the flux spectrum and eigenvalues
    if alpha_quench_in!=0.0 || vexb_shear_s==0.0 # do not calculate spectral shift

        #Threads.@threads
        for ky_index = eachindex(ky_spect)
            # print("this is c")
            onePass!(inputs, satParams, outputHermite, vexb_shear_s, firstPass_eigenvalue, QL_weights, ky_index)
        end

    elseif !inputs.FIND_WIDTH # calculate spectral shift with predetermined width spectrum

        # get the gammas to calculate the spectral shift on second pass
        inputs.IFLUX = false # do not compute QL on first pass
        #Threads.@threads 
        for ky_index = eachindex(ky_spect)
            inputs.IFLUX = false
            widthPass!(inputs, satParams, outputHermite, firstPass_eigenvalue, ky_index)
        end
        kx0_e = xgrid_functions_geo(inputs,satParams,firstPass_eigenvalue[:,:,1]) # spectral shift
        
        inputs.IFLUX = true # do not compute QL on first pass
        #Threads.@threads 
        for ky_index = eachindex(ky_spect)
            inputs.IFLUX = original_iflux
            secondPass!(inputs, satParams, outputHermite, kx0_e[ky_index], firstPass_eigenvalue, QL_weights, ky_index)
        end

    else # calculate the spectral shift and the width spectrum
        # initial guess if finding width
        inputs.WIDTH_SPECTRUM .= inputs.WIDTH
        inputs.IFLUX = false # do not compute QL on first pass
        #Threads.@threads 
        for ky_index = eachindex(ky_spect)
            # println("this is a")
            firstPass!(inputs, satParams, outputHermite, firstPass_eigenvalue, ky_index) #  spectral shift model double pass
        end
        kx0_e = xgrid_functions_geo(inputs,satParams,firstPass_eigenvalue[:,:,1]) # spectral shift

        inputs.IFLUX = true # do not compute QL on first pass
        #Threads.@threads 
        for ky_index = eachindex(ky_spect)
            # println("this is b")
            inputs.IFLUX = original_iflux # compute QL on second pass
            secondPass!(inputs, satParams, outputHermite, kx0_e[ky_index], firstPass_eigenvalue, QL_weights, ky_index)
        end

    end
    inputs.IFLUX = original_iflux

    return QL_weights, firstPass_eigenvalue

end

#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

"""
    function onePass!(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl

outputs:
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
        if(inputs.NEW_EIKONAL)
            # println("this is 1")
            nmodes_out, gamma_nb_min_out,
            gamma_out, freq_out,
            particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out = tjlf_max(inputs, satParams, outputHermite, ky, vexb_shear_s, ky_index)
        else
            error("NOT IMPLEMENTED YET -DSUN")
        end
    else
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
    function firstPass(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl

outputs:
    QL_weights                          - 5d array of QL weights (field, species, mode, ky, type),
                                          type: (particle, energy, torodial stress, parallel stress, exchange)
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
    function widthPass(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl

outputs:
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
    function secondPass(inputs::InputTJLF{T},satParams::SaturationParameters{T},outputHermite::OutputHermite{T},firstPass_eigenvalue::Array{T,3}) 

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl
    firstPass_eigenvalue::Array{T,3}    - array of eigenvalues found in first pass

outputs:
    QL_weights                          - 5d array of QL weights (field, species, mode, ky, type),
                                          type: (particle, energy, torodial stress, parallel stress, exchange)

description:
    calculate the QL weights using the width and eigenvalues from the first pass as reference
"""
function secondPass!(inputs::InputTJLF{T}, satParams::SaturationParameters{T},outputHermite::OutputHermite{T},
    kx0_e::T,
    firstPass_eigenvalue::Array{T,3}, QL_weights::Array{T,5}, ky_index::Int) where T<:Real

    # input values
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    nbasis_max_in = inputs.NBASIS_MAX
    nbasis_min_in = inputs.NBASIS_MIN
    nx = 2*inputs.NXGRID - 1
    vexb_shear_s = inputs.VEXB_SHEAR*inputs.SIGN_IT

    # saturation values
    R_unit = satParams.R_unit
    q_unit = satParams.q_unit

    # mask_save::Vector{Int} = zeros(Int, nky)

    # initialize reference eigenvalue arrays
    gamma_reference_kx0 = similar(firstPass_eigenvalue[:,1,1])
    freq_reference_kx0 = similar(firstPass_eigenvalue[:,1,2])

    ky = inputs.KY_SPECTRUM[ky_index]

    if(inputs.NEW_EIKONAL)
        gamma_reference_kx0 .= firstPass_eigenvalue[:,ky_index,1]
        freq_reference_kx0 .= firstPass_eigenvalue[:,ky_index,2]
        nbasis = nbasis_max_in
        # println("this is 3")

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
        
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,1] .= particle_QL_out[:,ns0:ns,1:nmodes_out]
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,2] .= energy_QL_out[:,ns0:ns,1:nmodes_out]
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,3] .= stress_tor_QL_out[:,ns0:ns,1:nmodes_out]
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,4] .= stress_par_QL_out[:,ns0:ns,1:nmodes_out]
        QL_weights[:,ns0:ns,1:nmodes_out,ky_index,5] .= exchange_QL_out[:,ns0:ns,1:nmodes_out]
        
    end

    # recompute spectrum using non-local in ky multiscale model
    # if(sat_rule_in >= 1) get_multiscale_spectrum(inputs) end
    # if(new_eikonal_in) eikonal_unsaved = false end
    # gamma_out[1] = gmax
    # freq_out[1] = fmax
    # inputs.NEW_EIKONAL = true  # reset default for next call to tjlf_TM

end