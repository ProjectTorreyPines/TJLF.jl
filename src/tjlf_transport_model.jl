include("tjlf_modules.jl")
include("tjlf_max.jl")


function get_bilinear_spectrum(inputs::InputTJLF, outputGeo::OutputGeometry, outputHermite::OutputHermite, ky_spect::Vector{T}, vexb_shear_s::T, jmax_out::Int) where T<:Real
#
# computes the bilinear fluctuation moments 
# and saves them in flux_spectrum_out, intensity_spectrum_out
# and field_spectrum_out
#
    ### DELETE THIS EVENTUALLY
    R_unit, q_unit= get_sat_params(:rq_units, inputs)

    sat_rule_in = inputs.SAT_RULE
    nmodes = inputs.NMODES
    ns = inputs.NS
    new_eikonal_in = inputs.NEW_EIKONAL
    find_width_in = inputs.FIND_WIDTH
    nbasis_max_in = inputs.NBASIS_MAX
    nbasis_min_in = inputs.NBASIS_MIN
    nxgrid_in = inputs.NXGRID
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end

    nx = 2*nxgrid_in - 1
    nky = length(ky_spect)

    # initialize output arrays
    spectral_shift_out = zeros(nky)
    ave_p0_spectrum_out = zeros(nky)
    eigenvalue_spectrum_out = zeros(2, nky, nmodes)
    field_spectrum_out = zeros(4, nky, nmodes)
    ### apparently dimension 2(labeled ns) may not be zero at 1
    ### if adiabatic_electron_in in inputs is true
    intensity_spectrum_out = zeros(4, ns, nky, nmodes)
    flux_spectrum_out = zeros(5, ns, 3, nky, nmodes)
    nsts_phase_spectrum_out = zeros(ns, nky, nmodes)
    ne_te_phase_spectrum_out = zeros(nky, nmodes)

    # save maximum width
    ##### this is in tjlf_max.f90 but its too much to decipher rn COME BACK!!!!!!
    width_in = 10^10
    save_width = width_in

    iflux_in = true 
    gmax = 0.0
    fmax = 0.0
    mask_save = zeros(Int, nky)
    QL_field_spectrum_out = Array{Float64, 3}(undef, 4, nky, nmodes)
    field_spectrum_out = Array{Float64, 3}(undef, 4, nky, nmodes)
    eigenvalue_spectrum_out = Array{Float64, 3}(undef, 2, nky, nmodes)
    QL_intensity_spectrum_out = Array{Float64, 4}(undef, 4, ns, nky, nmodes)
    intensity_spectrum_out = Array{Float64, 4}(undef, 4, ns, nky, nmodes)
    QL_flux_spectrum_out = Array{Float64, 5}(undef, 5, ns, 3, nky, nmodes)
    flux_spectrum_out = Array{Float64, 5}(undef, 5, ns, 3, nky, nmodes)
    for i = 1:nky
        ky_s = ky_spect[i]
        new_width = true

        if(new_eikonal_in)
            if(jmax_out == 0) ################ first pass? ###############
                if(find_width_in)
                    println("this is 1")
                    gamma_nb_min_out,
                    gamma_out,
                    freq_out,
                    v_QL_out, 
                    a_par_QL_out,
                    b_par_QL_out,
                    phi_bar_out,
                    a_par_bar_out,
                    b_par_bar_out,
                    v_bar_out,
                    ne_te_phase_out,
                    field_weight_out,
                    particle_QL_out,
                    energy_QL_out,
                    stress_par_QL_out,
                    stress_tor_QL_out,
                    exchange_QL_out,
                    N_QL_out,
                    T_QL_out,
                    U_QL_out,
                    Q_QL_out,
                    N_bar_out,
                    T_bar_out,
                    U_bar_out,
                    Q_bar_out,
                    Ns_Ts_phase_out = tjlf_max(inputs, outputGeo, outputHermite, ky_s, vexb_shear_s) ############### have to create this ###############
                else
                    nbasis = nbasis_max_in
                    new_width = true
                    tjlf_LS() ############### have to create this ###############
                    gamma_nb_min_out = gamma_out[1]
                end
            else   # second pass
                gamma_reference_kx0 .= eigenvalue_first_pass[1,i,:]
                freq_reference_kx0 .= eigenvalue_first_pass[2,i,:]
                width_in = width_out[i]
                nbasis = nbasis_max_in
                new_width = true
                tjlf_LS() ############### have to create this ###############
                gamma_nb_min_out = gamma_out[1]
            end
            
            ############### need these values from the function calls earlier probably
            mask_save[i] = 1
            if(gamma_out[1] == 0.0) mask_save[i] = 0 end
            #### not used yet
            gamma_nb_min_save[i] = gamma_nb_min_out
            width_save[i] = width_in
            ft_save[i] = ft
            R_unit_save[i] = R_unit
            q_unit_save[i] = q_unit
                ############ need these values from the function calls earlier probably
            for j = 1:nx
                wdx_save[i,j] = wdx[j]
                b0x_save[i,j] = b0x[j]
                b2x_save[i,j] = b2x[j]
                cx_par_par_save[i,j] = cx_par_par[j]
                cx_tor_par_save[i,j] = cx_tor_par[j]
                cx_tor_per_save[i,j] = cx_tor_per[j]
                kxx_save[i,j] = kxx[j]
            end
        else
            ############### need these values from the function calls earlier probably
            gamma_nb_min_out = gamma_nb_min_save[i]
            width_in = width_save[i]
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
                tjlf_LS() ############### have to create this ###############
            else
                gamma_out[1] = 0.0
            end
        end

        unstable = true
        gamma_max = max(gamma_out[1],gamma_out[2]) # this covers ibranch=-1,0
        if(gamma_max == 0.0 || gamma_nb_min_out == 0.0) unstable = false end     
         
        gamma_net_1 = gamma_nb_min_out 
        gamma_cutoff = 0.1*ky_s/R_unit
        rexp = 1.0
        reduce = 1.0
        if(nbasis_max_in != nbasis_min_in)
            if(gamma_net_1 < gamma_out[1] && gamma_net_1 < gamma_cutoff)
                reduce = (gamma_net_1/gamma_cutoff)^rexp
            end
        end
        if(sat_rule_in > 1) reduce = 1.0 end


        # width_out[i] = width_in
        if(unstable)
            # save the spectral shift of the radial wavenumber due to VEXB_SHEAR
            # spectral_shift_out[i] = kx0_e
            # ave_p0_spectrum_out[i] = ave_p0_out
            # save field_spectrum_out and eigenvalue_spectrum_out
            nmodes_out = nmodes #### this might be wrong?
            for imax = 1:nmodes_out
                QL_field_spectrum_out[1,i,imax] = v_QL_out[imax]
                QL_field_spectrum_out[2,i,imax] = 1.0
                QL_field_spectrum_out[3,i,imax] = a_par_QL_out[imax]
                QL_field_spectrum_out[4,i,imax] = b_par_QL_out[imax]
                field_spectrum_out[1,i,imax] = reduce*v_bar_out[imax]
                field_spectrum_out[2,i,imax] = reduce*phi_bar_out[imax]
                field_spectrum_out[3,i,imax] = reduce*a_par_bar_out[imax]
                field_spectrum_out[4,i,imax] = reduce*b_par_bar_out[imax]
                eigenvalue_spectrum_out[1,i,imax]=gamma_out[imax]
                eigenvalue_spectrum_out[2,i,imax]=freq_out[imax]

                if(ky_s <= 1.0 && gamma_out[imax] > gmax)
                    gmax=gamma_out[imax]
                    fmax=freq_out[imax]
                end
            end
            # save intensity_spectrum_out
            for is = ns0:ns
                for imax = 1:nmodes_out
                    QL_intensity_spectrum_out[1,is,i,imax] = N_QL_out[imax,is]
                    QL_intensity_spectrum_out[2,is,i,imax] = T_QL_out[imax,is]
                    QL_intensity_spectrum_out[3,is,i,imax] = U_QL_out[imax,is]
                    QL_intensity_spectrum_out[4,is,i,imax] = Q_QL_out[imax,is]
                    intensity_spectrum_out[1,is,i,imax] = N_bar_out[imax,is]
                    intensity_spectrum_out[2,is,i,imax] = T_bar_out[imax,is]
                    intensity_spectrum_out[3,is,i,imax] = U_bar_out[imax,is]
                    intensity_spectrum_out[4,is,i,imax] = Q_bar_out[imax,is]
                end #imax
            end  # is
            # save flux_spectrum_out 
            for is = ns0:ns
                for j = 1:3
                    for imax = 1:nmodes_out
                        phi_bar = reduce*phi_bar_out[imax]
                        pflux1 = particle_QL_out[imax,is,j]
                        eflux1 = energy_QL_out[imax,is,j]
                        stress_tor1 = stress_tor_QL_out[imax,is,j]
                        stress_par1 = stress_par_QL_out[imax,is,j]
                        exch1 = exchange_QL_out[imax,is,j]
                        QL_flux_spectrum_out[1,is,j,i,imax] = pflux1
                        QL_flux_spectrum_out[2,is,j,i,imax] = eflux1
                        QL_flux_spectrum_out[3,is,j,i,imax] = stress_tor1
                        QL_flux_spectrum_out[4,is,j,i,imax] = stress_par1
                        QL_flux_spectrum_out[5,is,j,i,imax] = exch1
                        flux_spectrum_out[1,is,j,i,imax] = phi_bar*pflux1
                        flux_spectrum_out[2,is,j,i,imax] = phi_bar*eflux1
                        flux_spectrum_out[3,is,j,i,imax] = phi_bar*stress_tor1
                        flux_spectrum_out[4,is,j,i,imax] = phi_bar*stress_par1
                        flux_spectrum_out[5,is,j,i,imax] = phi_bar*exch1
                    end #imax
                end # j
            end  # is 
            # save ne_te crossphase
            for imax = 1:nmodes_out
                ne_te_phase_spectrum_out[i,imax] = ne_te_phase_out[imax]
            end  #imax
            # save ns_ts crossphase
            for is = 1:ns
                for imax = 1:nmodes_out
                    nsts_phase_spectrum_out[is,i,imax] = Ns_Ts_phase_out[imax,is]
                end  #imax
            end    # is
        end #unstable .T.

        # reset width to maximum if used tjlf_max
        if(find_width_in) width_in=save_width end

    end  # i 

    # recompute spectrum using non-local in ky multiscale model
    # if(sat_rule_in >= 1) get_multiscale_spectrum(inputs) end

    # if(new_eikonal_in) eikonal_unsaved = false end
    # gamma_out[1] = gmax
    # freq_out[1] = fmax
    # inputs.NEW_EIKONAL = true  # reset default for next call to tjlf_TM


    # for is=ns0:ns
    # for j=1:3
    # for m=1:nmodes
    # for i=1:nky
    # for k=1:5
    #        print(QL_flux_spectrum_out[k,is,j,i,m])
    #        print(" ")
    # end
    # println()
    # end
    # end
    # end
    # end

    return eigenvalue_spectrum_out
      
end
