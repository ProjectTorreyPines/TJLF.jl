"""
    function tjlf_TM(inputs::InputTJLF{T},satParams::SaturationParameters{T},outputHermite::OutputHermite{T}) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl

outputs:
    fluxes                              - 5d array of QL fluxes (field, species, mode, ky, type),
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

    original_iflux = inputs.IFLUX
    original_width = inputs.WIDTH

    # compute the flux spectrum

    if alpha_quench_in!=0.0 || vexb_shear_s==0.0 # do not calculate spectral shift
        # print("this is c")
        fluxes, firstPass_eigenvalue = onePass(inputs, satParams, outputHermite, vexb_shear_s)

    elseif !inputs.FIND_WIDTH # calculate spectral shift with predetermined width spectrum
        # get the gammas to calculate the spectral shift on second pass
        inputs.IFLUX = false
        firstPass_eigenvalue = widthPass(inputs, satParams, outputHermite)

        inputs.IFLUX = original_iflux
        fluxes = secondPass(inputs, satParams, outputHermite, firstPass_eigenvalue)

    else # calculate the spectral shift and the width spectrum
        #  spectral shift model double pass
        # println("this is a")
        inputs.IFLUX = false # do not compute QL on first pass
        firstPass_eigenvalue = firstPass(inputs, satParams, outputHermite)

        # println("this is b")
        inputs.IFLUX = original_iflux # compute QL on second pass
        fluxes = secondPass(inputs, satParams, outputHermite, firstPass_eigenvalue)
    end

    inputs.WIDTH = original_width

    return fluxes, firstPass_eigenvalue

    # sum_ky_spectrum(inputs, eigenvalue_spectrum_out[1,:,:],
    #             ave_p0,potential,
    #             particle_QL,
    #             energy_QL,
    #             toroidal_stress_QL,
    #             parallel_stress_QL,
    #             exchange_QL)

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
    fluxes                              - 5d array of QL fluxes (field, species, mode, ky, type),
                                          type: (particle, energy, torodial stress, parallel stress, exchange)
    firstPass_eigenvalue                - 3d array of eigenvalues (mode, ky, type)
                                          type: (gamma, frequency)

description:
    calculate the widths and eigenvalues with vexb_shear = 0.0
"""
function firstPass(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}) where T<:Real

    nmodes = inputs.NMODES
    new_eikonal_in = inputs.NEW_EIKONAL
    ky_spect = inputs.KY_SPECTRUM
    nky = length(ky_spect)

    original_width = inputs.WIDTH

    ### output values
    eigenvalue_spectrum_out = zeros(Float64, nmodes, nky, 2)

    # increment through the ky_spectrum and find the width/eigenvalues of each ky
    for i = eachindex(ky_spect)
        ky = ky_spect[i]

        if(new_eikonal_in) # not sure what this is -DSUN
            # println("this is 1")
            nmodes_out, gamma_nb_min_out,
            gamma_out, freq_out,
            _,_,_,_,_ = tjlf_max2(inputs, satParams, outputHermite, ky, 0.0)

            inputs.WIDTH_SPECTRUM[i] = inputs.WIDTH
            ### reset value
            inputs.WIDTH = original_width
        else
            error("NOT IMPLEMENTED YET -DSUN")
        end

        unstable = true
        gamma_max = max(gamma_out[1],gamma_out[2]) # this covers ibranch=-1,0
        if(gamma_max == 0.0 || gamma_nb_min_out == 0.0) unstable = false end

        if(unstable)
            eigenvalue_spectrum_out[1:nmodes_out,i,1] .= gamma_out[1:nmodes_out]
            eigenvalue_spectrum_out[1:nmodes_out,i,2] .= freq_out[1:nmodes_out]
        end
    end

    return eigenvalue_spectrum_out

end

#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

"""
    function onePass(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}) where T<:Real

parameters:
    inputs::InputTJLF{T}                - InputTJLF struct constructed in tjlf_read_input.jl
    satParams::SaturationParameters{T}  - SaturationParameters struct constructed in tjlf_geometry.jl
    outputHermite::OutputHermite{T}     - OutputHermite struct constructed in tjlf_hermite.jl

outputs:
    QL_flux_spectrum_out                - 5d array of QL fluxes (field, species, mode, ky, type),
                                          type: (particle, energy, torodial stress, parallel stress, exchange)
    eigenvalue_spectrum_out              - 3d array of eigenvalues (mode, ky, type)
                                          type: (gamma, frequency)

description:
    calculate both eigenvalues and fluxes on firstPass if vexb_shear = 0 or alpha_quench != 0, no secondPass
"""
function onePass(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}, vexb_shear_s::T) where T<:Real

    nmodes = inputs.NMODES
    new_eikonal_in = inputs.NEW_EIKONAL
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    ky_spect = inputs.KY_SPECTRUM
    nky = length(ky_spect)

    original_width = inputs.WIDTH

    ### output values
    eigenvalue_spectrum_out = zeros(Float64, nmodes, nky, 2)
    QL_flux_spectrum_out::Array{Float64,5} = zeros(Float64, 3, ns, nmodes, nky, 5)

    # increment through the ky_spectrum and find the width/eigenvalues of each ky
    for i = eachindex(ky_spect)
        ky = ky_spect[i]
        inputs.WIDTH = inputs.WIDTH_SPECTRUM[i]

        if(new_eikonal_in) # not sure what this is -DSUN
            # println("this is 1")
            nmodes_out, gamma_nb_min_out,
            gamma_out, freq_out,
            particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out = tjlf_max2(inputs, satParams, outputHermite, ky, vexb_shear_s)
        else
            error("NOT IMPLEMENTED YET -DSUN")
        end

        unstable = true
        gamma_max = max(gamma_out[1],gamma_out[2]) # this covers ibranch=-1,0
        if(gamma_max == 0.0 || gamma_nb_min_out == 0.0) unstable = false end

        if(unstable)
            eigenvalue_spectrum_out[1:nmodes_out,i,1] .= gamma_out[1:nmodes_out]
            eigenvalue_spectrum_out[1:nmodes_out,i,2] .= freq_out[1:nmodes_out]
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,1] .= particle_QL_out[:,ns0:ns,1:nmodes_out]
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,2] .= energy_QL_out[:,ns0:ns,1:nmodes_out]
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,3] .= stress_tor_QL_out[:,ns0:ns,1:nmodes_out]
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,4] .= stress_par_QL_out[:,ns0:ns,1:nmodes_out]
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,5] .= exchange_QL_out[:,ns0:ns,1:nmodes_out]
        end
    end

    inputs.WIDTH = original_width

    return QL_flux_spectrum_out, eigenvalue_spectrum_out

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
function widthPass(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}) where T<:Real

    nmodes = inputs.NMODES
    ky_spect = inputs.KY_SPECTRUM
    nky = length(ky_spect)

    original_width = inputs.WIDTH

    if(!inputs.NEW_EIKONAL)
        error("not sure what this flag is -DSUN")
    end

    ### output values
    eigenvalue_spectrum_out = zeros(Float64, nmodes, nky, 2)

    # increment through the ky/width spectrum
    for i = eachindex(ky_spect)
        ky = ky_spect[i]
        inputs.WIDTH = inputs.WIDTH_SPECTRUM[i]
        if inputs.WIDTH<0.0 # invalid value found in max()
            continue
        end

        # calculate the eigenvalues with no shear
        nbasis = inputs.NBASIS_MAX
        nmodes_out, gamma_out, freq_out,
        _,_,_,_,_ = tjlf_LS(inputs, satParams, outputHermite, ky, nbasis, 0.0)

        if(inputs.IBRANCH==-1) # check for inward ballooning modes
            if(inputs.USE_INBOARD_DETRAPPED && ft_test > modB_test) ####### find ft_test and modB_test
                error("not implemented XIII")
            end
        end

        unstable = true
        gamma_max = max(gamma_out[1],gamma_out[2]) # this covers ibranch=-1,0
        if(gamma_max == 0.0) unstable = false end ############################ deleted gammma_nb_min_out ############################

        if(unstable)
            eigenvalue_spectrum_out[1:nmodes_out,i,1] .= gamma_out[1:nmodes_out]
            eigenvalue_spectrum_out[1:nmodes_out,i,2] .= freq_out[1:nmodes_out]
        end
    end
    inputs.WIDTH = original_width

    return eigenvalue_spectrum_out

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
    QL_flux_spectrum_out                - 5d array of QL fluxes (field, species, mode, ky, type),
                                          type: (particle, energy, torodial stress, parallel stress, exchange)

description:
    calculate the fluxes using the width and eigenvalues from the first pass as reference
"""
function secondPass(inputs::InputTJLF{T}, satParams::SaturationParameters{T},outputHermite::OutputHermite{T},
    firstPass_eigenvalue::Array{T,3}) where T<:Real

    ### input values
    sat_rule_in = inputs.SAT_RULE
    nmodes = inputs.NMODES
    ns = inputs.NS
    new_eikonal_in = inputs.NEW_EIKONAL
    nbasis_max_in = inputs.NBASIS_MAX
    nbasis_min_in = inputs.NBASIS_MIN
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    nx = 2*inputs.NXGRID - 1
    ky_spect = inputs.KY_SPECTRUM
    nky = length(ky_spect)
    vexb_shear_s = inputs.VEXB_SHEAR*inputs.SIGN_IT
    ### saturation values
    R_unit = satParams.R_unit
    q_unit = satParams.q_unit
    ### change input value
    inputs.IFLUX = true

    mask_save::Vector{Int} = zeros(Int, nky)

    # initialize output arrays
    QL_flux_spectrum_out::Array{Float64,5} = zeros(Float64, 3, ns, nmodes, nky, 5)
    gamma_reference_kx0 = similar(firstPass_eigenvalue[:,1,1])
    freq_reference_kx0 = similar(firstPass_eigenvalue[:,1,2])
    kx0_e = xgrid_functions_geo(inputs,satParams,firstPass_eigenvalue[1,:,:])

    for i = eachindex(ky_spect)
        ky = ky_spect[i]

        if(new_eikonal_in)
            gamma_reference_kx0 .= firstPass_eigenvalue[:,i,1]
            freq_reference_kx0 .= firstPass_eigenvalue[:,i,2]
            inputs.WIDTH = abs(inputs.WIDTH_SPECTRUM[i])
            nbasis = nbasis_max_in
            # println("this is 3")

            nmodes_out, gamma_out, freq_out,
            particle_QL_out,
            energy_QL_out,
            stress_tor_QL_out,
            stress_par_QL_out,
            exchange_QL_out = tjlf_LS(inputs, satParams, outputHermite, ky, nbasis, vexb_shear_s,
                                kx0_e[i], gamma_reference_kx0, freq_reference_kx0)

            gamma_nb_min_out = gamma_out[1]

            mask_save[i] = 1
            if(gamma_out[1] == 0.0) mask_save[i] = 0 end
        else
        ############### need these values from the function calls earlier probably
            error("NOT IMPLEMENTED YET -DSUN")
            gamma_nb_min_out = gamma_nb_min_save[i]
            inputs.WIDTH = width_save[i]
            ft = ft_save[i]
            R_unit = R_unit_save[i]
            q_unit = q_unit_save[i]
            for j = 1:nx
                wdx[j] = wdx_save[i,j]
                b0x[j] = b0x_save[i,j]
                b2x[j] = b2x_save[i,j]
                cx_par_par[j] = cx_par_par_save[i,j]
                cx_tor_par[j] = cx_tor_par_save[i,j]
                cx_tor_per[j] = cx_tor_per_save[i,j]
                kxx[j] = kxx_save[i,j]
            end

            if(mask_save[i] == 1)
                # println("this is 4")
                tjlf_LS() ############### have to create this ###############
            else
                gamma_out[1] = 0.0
            end
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
            if(sat_rule_in >= 1) reduce = 1.0 end
            # # save the spectral shift of the radial wavenumber due to VEXB_SHEAR
            # spectral_shift_out[i] = kx0_e
            # ave_p0_spectrum_out[i] = ave_p0_out
            # # save field_spectrum_out and eigenvalue_spectrum_out

            # for imax = 1:nmodes_out
            #     QL_field_spectrum_out[1,i,imax] = v_QL_out[imax]
            #     QL_field_spectrum_out[2,i,imax] = 1.0
            #     QL_field_spectrum_out[3,i,imax] = a_par_QL_out[imax]
            #     QL_field_spectrum_out[4,i,imax] = b_par_QL_out[imax]
            #     field_spectrum_out[1,i,imax] = reduce*v_bar_out[imax]
            #     field_spectrum_out[2,i,imax] = reduce*phi_bar_out[imax]
            #     field_spectrum_out[3,i,imax] = reduce*a_par_bar_out[imax]
            #     field_spectrum_out[4,i,imax] = reduce*b_par_bar_out[imax]
            #     if(ky <= 1.0 && gamma_out[imax] > gmax)
            #         gmax=gamma_out[imax]
            #         fmax=freq_out[imax]
            #     end
            # end
            # save intensity_spectrum_out
            # for is = ns0:ns
            #     for imax = 1:nmodes_out
            #         QL_intensity_spectrum_out[1,is,i,imax] = N_QL_out[imax,is]
            #         QL_intensity_spectrum_out[2,is,i,imax] = T_QL_out[imax,is]
            #         QL_intensity_spectrum_out[3,is,i,imax] = U_QL_out[imax,is]
            #         QL_intensity_spectrum_out[4,is,i,imax] = Q_QL_out[imax,is]
            #         intensity_spectrum_out[1,is,i,imax] = N_bar_out[imax,is]
            #         intensity_spectrum_out[2,is,i,imax] = T_bar_out[imax,is]
            #         intensity_spectrum_out[3,is,i,imax] = U_bar_out[imax,is]
            #         intensity_spectrum_out[4,is,i,imax] = Q_bar_out[imax,is]
            #     end #imax
            # end  # is
            # save flux_spectrum_out
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,1] .= particle_QL_out[:,ns0:ns,1:nmodes_out]
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,2] .= energy_QL_out[:,ns0:ns,1:nmodes_out]
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,3] .= stress_tor_QL_out[:,ns0:ns,1:nmodes_out]
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,4] .= stress_par_QL_out[:,ns0:ns,1:nmodes_out]
            QL_flux_spectrum_out[:,ns0:ns,1:nmodes_out,i,5] .= exchange_QL_out[:,ns0:ns,1:nmodes_out]
            # for is = ns0:ns
            #     for j = 1:3
            #         for imax = 1:nmodes_out
                        # phi_bar = reduce*phi_bar_out[imax]
                        # pflux1 = particle_QL_out[imax,is,j]
                        # eflux1 = energy_QL_out[imax,is,j]
                        # stress_tor1 = stress_tor_QL_out[imax,is,j]
                        # stress_par1 = stress_par_QL_out[imax,is,j]
                        # exch1 = exchange_QL_out[imax,is,j]
                        # QL_flux_spectrum_out[j,is,imax,i,1] = pflux1
                        # QL_flux_spectrum_out[j,is,imax,i,2] = eflux1
                        # QL_flux_spectrum_out[j,is,imax,i,3] = stress_tor1
                        # QL_flux_spectrum_out[j,is,imax,i,4] = stress_par1
                        # QL_flux_spectrum_out[j,is,imax,i,5] = exch1
                        # flux_spectrum_out[1,is,j,i,imax] = phi_bar*pflux1
                        # flux_spectrum_out[2,is,j,i,imax] = phi_bar*eflux1
                        # flux_spectrum_out[3,is,j,i,imax] = phi_bar*stress_tor1
                        # flux_spectrum_out[4,is,j,i,imax] = phi_bar*stress_par1
                        # flux_spectrum_out[5,is,j,i,imax] = phi_bar*exch1
            #         end #imax
            #     end # j
            # end  # is

        end

    end  # i

    # recompute spectrum using non-local in ky multiscale model
    # if(sat_rule_in >= 1) get_multiscale_spectrum(inputs) end
    # if(new_eikonal_in) eikonal_unsaved = false end
    # gamma_out[1] = gmax
    # freq_out[1] = fmax
    # inputs.NEW_EIKONAL = true  # reset default for next call to tjlf_TM

    return QL_flux_spectrum_out

end