function intensity_sat1(inputs::InputTJLF, ky_spect::Vector{T}, gp::Array{T}, QL_data::Array{T}, expsub::T=2.0, return_phi_params::Bool=false) where T<: Real
    
    # use the growth rates of the most unstable modes
    length(size(gp)) > 1 ? gammas1 = gp[:, 1] : gammas1 = gp

    kx0_e, SAT_geo0_out, _, _, _, _, _, _, _, _, _, _ = get_sat_params(inputs,ky_spect,gp)
    
    rlnp_cutoff = inputs.RLNP_CUTOFF
    beta_loc = inputs.BETA_LOC
    ns = inputs.NS
    rmaj_loc = inputs.RMAJ_LOC
    p_prime_loc = inputs.P_PRIME_LOC
    alpha_zf = inputs.ALPHA_ZF
    alpha_quench = inputs.ALPHA_QUENCH

    zs_2 = inputs.SPECIES[2].ZS
    taus_2 = inputs.SPECIES[2].TAUS
    mass_2 = inputs.SPECIES[2].MASS
    taus_1 = inputs.SPECIES[1].TAUS
    rho_ion = abs(zs_2) / √(taus_2*mass_2)

    small = 10^-10
    nky = length(ky_spect)
    nmodes = length(QL_data[1, :, 1, 1, 1])
    gamma_net = zeros(nky)

    if(rlnp_cutoff > 0.0)
        if(beta_loc == 0.0)
            dlnpdr = 0.0
            ptot = 0.0
            for i in 1:ns
                ptot = ptot + inputs.SPECIES[i].AS*inputs.SPECIES[i].TAUS
                dlnpdr = (dlnpdr + 
                    inputs.SPECIES[i].AS*inputs.SPECIES[i].TAUS*
                    (inputs.SPECIES[i].RLNS+inputs.SPECIES[i].RLTS))
            end
            ### kwargs["RMAJ_LOC"] used for rmaj_input
            dlnpdr = rmaj_loc*dlnpdr/max(ptot,0.01)
        else
            dlnpdr = (-p_prime_loc*(8.0π/beta_loc)*
                (rmin_input/q_loc)*rmaj_input) ##### need to define rmin_input
        end

        if(dlnpdr > rlnp_cutoff) dlnpdr = rlnp_cutoff_in end
        if(dlnpdr < 4.0) dlnpdr = 4.0 end

    else
        dlnpdr = 12.0
    end

    ### coefficients for SAT_RULE = 1
    czf = abs(alpha_zf)
    cnorm = 14.29
    cz1=0.48*czf
    cz2=1.0*czf
    cky=3.0
    sqcky=√(cky)
    etg_streamer=1.05
    kyetg = etg_streamer * abs(zs_2) / √(taus_2 * mass_2)
    measure = √(taus_1 * mass_2)
    if(alpha_quench != 0.0) etg_streamer=2.1 end

    # coefficients for spectral shift model for ExB shear
    ax=0.0
    ay=0.0
    exp_ax = 1
    if(alpha_quench == 0.0)
        #spectral shift model parameters
        ax = 1.15
        ay = 0.56
        exp_ax = 4
    end

    gamma_net .= gammas1 ./ (1.0 .+ abs.(ax.*kx0_e).^exp_ax)
    vzf_out, kymax_out, jmax_out = get_zonal_mixing(inputs, ky_spect, gamma_net)

    # compute multi-scale phi-intensity spectrum field_spectrum(2,,) = phi_bar_out
    # note that the field_spectrum(1,,) = v_bar_out = 1.0 for sat_rule_in = 1
    gammamax1= vzf_out*kymax_out
    kymax1 = kymax_out
    jmax1 = jmax_out
    vzf1 = vzf_out

    # include zonal flow effects on growth rate model:
    gamma_mix1 = zeros(nky)
    gamma = zeros(nky)

    gamma .= ifelse.(ky_spect.<kymax1, 
                max.(gamma_net .- (cz1.*(kymax1 .- ky_spect).*vzf1),0.0),
                (cz2*gammamax1) .+ max.(gamma_net .- (cz2.*vzf1.*ky_spect),0.0))

    gamma_mix1 .= gamma

    mixnorm = 0.0
    for j in jmax1+2:nky
        gamma_ave = 0.0
        ky0 = ky_spect[j]
        mixnorm = (ky0*
                    (atan(sqcky*(ky_spect[nky]/ky0 - 1.0))
                    - atan(sqcky*(ky_spect[jmax1+1]/ky0 - 1.0)))
        )
        for i in jmax1+1:nky-1
            ky1 = ky_spect[i]
            ky2 = ky_spect[i+1]
            mix1 = (ky0*
                    (atan(sqcky*(ky2/ky0 - 1.0))
                    - atan(sqcky*(ky1/ky0-1.0)))
            )
            delta = (gamma[i+1]-gamma[i])/(ky2-ky1)
            mix2 = (ky0*mix1 
                + (ky0*ky0/(2.0*sqcky))
                * (log(cky*(ky2-ky0)^2 + ky0^2)
                -  log(cky*(ky1-ky0)^2 + ky0^2))
            )
            gamma_ave = gamma_ave + (gamma[i] - ky1*delta)*mix1 + delta*mix2
        end
        gamma_mix1[j] = gamma_ave/mixnorm
    end

    # intensity model

    field_spectrum_out = zeros((nky,nmodes))
    gammaeff_out = zeros((nky,nmodes))
    kx_width_out = zeros(nky)
    sat_geo_factor_out = zeros(nky)

    for j in 1:nky 
        gamma0 = gammas1[j]
        ky0 = ky_spect[j]
        kx = kx0_e[j]

        sat_geo_factor = SAT_geo0_out
        kx_width = ky0

        kx_width_out[j] = kx_width
        sat_geo_factor_out[j] = sat_geo_factor

        
        for i in 1:nmodes
            gammaeff = 0.0
            if(gamma0>small)
                gammaeff = gamma_mix1[j]*(gp[j,i]/gamma0)^expsub
            end
            if(ky0>kyetg) gammaeff = gammaeff*√(ky0/kyetg) end

            field_spectrum_out[j,i] = (measure*cnorm*
                            ((gammaeff/(kx_width*ky0))/(1.0+ay*kx^2))^2)

            if(inputs.UNITS != "GYRO")
                field_spectrum_out[j,i] = sat_geo_factor*field_spectrum_out[j,i]
            end

            gammaeff_out[j, i] = gammaeff
        end    
    end
    
        
    #SAT3 QLA here
    QLA_P = 1.0
    QLA_E = 1.0
    QLA_O = 1.0

    phinorm = field_spectrum_out
    # so the normal behavior doesn't change,
    if return_phi_params
        out = Dict(
            (kx_width=kx_width_out),  # [nky] kx_model (kx rms width)
            (gammaeff=gammaeff_out),  # [nky, nmodes] effective growthrate
            (kx0_e=kx0_e),  # [nky] spectral shift in kx
            (ax=ax),  # SAT1 (cx), SAT2 (alpha_x)
            (ay=ay),  # SAT1 (cy)
            (exp_ax=exp_ax),  # SAT2 (sigma_x)
        )
    else
        out = phinorm, QLA_P, QLA_E, QLA_O  # SAT123 intensity and QLA params
    end

    return out
       
end


#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################




function intensity_sat2(inputs::InputTJLF, ky_spect::Vector{T}, gp::Array{T}, QL_data::Array{T}, expsub::T=2.0, return_phi_params::Bool=false) where T<: Real
    
    # use the growth rates of the most unstable modes
    length(size(gp)) > 1 ? gammas1 = gp[:, 1] : gammas1 = gp

    kx0_e, SAT_geo0_out, SAT_geo1_out, SAT_geo2_out, _, Bt0_out, B_geo0_out, grad_r0_out, _, _, _, _ = get_sat_params(inputs,ky_spect,gp)
    
    rlnp_cutoff = inputs.RLNP_CUTOFF
    beta_loc = inputs.BETA_LOC
    ns = inputs.NS
    rmaj_loc = inputs.RMAJ_LOC
    p_prime_loc = inputs.P_PRIME_LOC
    alpha_zf = inputs.ALPHA_ZF
    alpha_quench = inputs.ALPHA_QUENCH

    zs_2 = inputs.SPECIES[2].ZS
    taus_2 = inputs.SPECIES[2].TAUS
    mass_2 = inputs.SPECIES[2].MASS
    taus_1 = inputs.SPECIES[1].TAUS
    rho_ion = abs(zs_2) / √(taus_2*mass_2)

    small = 10^-10
    nky = length(ky_spect)
    nmodes = length(QL_data[1, :, 1, 1, 1])
    gamma_net = zeros(nky)

    if(rlnp_cutoff > 0.0)
        if(beta_loc == 0.0)
            dlnpdr = 0.0
            ptot = 0.0
            for i::Integer in 1:ns
                ptot = ptot + inputs.SPECIES[i].AS*inputs.SPECIES[i].TAUS
                dlnpdr = (dlnpdr + 
                    inputs.SPECIES[i].AS*inputs.SPECIES[i].TAUS*
                    (inputs.SPECIES[i].RLNS+inputs.SPECIES[i].RLTS))
            end
            ### kwargs["RMAJ_LOC"] used for rmaj_input
            dlnpdr = rmaj_loc*dlnpdr/max(ptot,0.01)
        else
            ##### need to define rmin_input
            dlnpdr = (-p_prime_loc*(8.0π/beta_loc)*
                (rmin_input/q_loc)*rmaj_input)
        end

        if(dlnpdr > rlnp_cutoff) dlnpdr = rlnp_cutoff_in end
        if(dlnpdr < 4.0) dlnpdr = 4.0 end

    else
        dlnpdr = 12.0
    end


    etg_streamer=1.05

    if(alpha_quench != 0.0) etg_streamer=2.1 end
    
    b0 = 0.76
    b1 = 1.22
    b2 = ifelse(nmodes > 1, 3.55, 3.74)
    b3 = 1.0

    d1 = (Bt0_out/b_geo0_out)^4/grad_r0_out   # PPCF paper 2020
    Gq = b_geo0_out/grad_r0_out
    d2 = b3/(Gq^2)

    czf = abs(alpha_zf)
    cnorm = b2*(12.0/dlnpdr)
    kyetg = 1000.0   # does not impact SAT2
    cky = 3.0
    sqcky = √(cky)
    kycut = b0*kymax_out
    cz1 = 0.0
    cz2 = 1.05*czf
    measure = 1.0/kymax_out


    # coefficients for spectral shift model for ExB shear
    ax=0.0
    ay=0.0
    exp_ax = 1
    if(alpha_quench == 0.0)
        ax = 1.21
        ay = 1.0
        exp_ax = 2
        units_in = "CGYRO"
    end

    for i in 1:nky
        kx = kx0_e[i]
        ky0 = ky_spect[i]
        if(ky0<kycut)
            kx_width = kycut/grad_r0_out
        else
            kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
        end
        kx = kx*ky0/kx_width
        gamma_net[i] = gammas1[i]/(1.0 + abs(ax*kx)^exp_ax)
    end

    vzf_out, kymax_out, jmax_out = get_zonal_mixing(inputs, ky_spect, gammas1)
    vzf_out_fp = vzf_out # used for recreation of first pass for SAT3
    vzf_out = vzf_out*gamma_net[jmax_out]/max(gammas1[jmax_out], small)


    # compute multi-scale phi-intensity spectrum field_spectrum(2,,) = phi_bar_out
    # note that the field_spectrum(1,,) = v_bar_out = 1.0 for sat_rule_in = 1
    gammamax1= vzf_out*kymax_out
    kymax1 = kymax_out
    jmax1 = jmax_out
    vzf1 = vzf_out

    # include zonal flow effects on growth rate model:
    gamma_mix1 = zeros(nky)
    gamma = zeros(nky)

    gamma .= ifelse.(ky_spect.<kymax1,
                gamma_net,
                gammamax1 .+ max.(gamma_net .- (cz2*vzf1.*ky_spect), 0.0))
    gamma_mix1 .= gamma




    mixnorm = 0.0
    for j in jmax1+2:nky
        gamma_ave = 0.0
        ky0 = ky_spect[j]
        mixnorm = (ky0*
                    (atan(sqcky*(ky_spect[nky]/ky0 - 1.0))
                    - atan(sqcky*(ky_spect[jmax1+1]/ky0 - 1.0)))
        )
        for i in jmax1+1:nky-1
            ky1 = ky_spect[i]
            ky2 = ky_spect[i+1]
            mix1 = (ky0*
                    (atan(sqcky*(ky2/ky0 - 1.0))
                    - atan(sqcky*(ky1/ky0-1.0)))
            )
            delta = (gamma[i+1]-gamma[i])/(ky2-ky1)
            mix2 = (ky0*mix1 
                + (ky0*ky0/(2.0*sqcky))
                * (log(cky*(ky2-ky0)^2 + ky0^2)
                -  log(cky*(ky1-ky0)^2 + ky0^2))
            )
            gamma_ave = gamma_ave + (gamma[i] - ky1*delta)*mix1 + delta*mix2
        end
        gamma_mix1[j] = gamma_ave/mixnorm
    end
  
    
    
    # intensity model

    field_spectrum_out = zeros((nky,nmodes))
    gammaeff_out = zeros((nky,nmodes))
    kx_width_out = zeros(nky)
    sat_geo_factor_out = zeros(nky)

    for j in 1:nky 
        gamma0 = gammas1[j]
        ky0 = ky_spect[j]
        kx = kx0_e[j]
        
        if(ky0 < kycut)
            kx_width = kycut/ grad_r0_out
            sat_geo_factor = SAT_geo0_out * d1 * SAT_geo1_out
        else
            kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
            sat_geo_factor = (SAT_geo0_out*
                                (d1*SAT_geo1_out*kycut +
                                (ky0 - kycut)*d2*SAT_geo2_out)/ky0)
        end
        kx = kx*ky0/kx_width

        kx_width_out[j] = kx_width
        sat_geo_factor_out[j] = sat_geo_factor

        
        for i in 1:nmodes
            gammaeff = 0.0
            if(gamma0>small)
                gammaeff = gamma_mix1[j]*(gp[j,i]/gamma0)^expsub
            end
            if(ky0>kyetg) gammaeff = gammaeff*√(ky0/kyetg) end

            field_spectrum_out[j,i] = (measure*cnorm*
                            ((gammaeff/(kx_width*ky0))/(1.0+ay*kx^2))^2)

            if(inputs.UNITS != "GYRO")
                field_spectrum_out[j,i] = sat_geo_factor*field_spectrum_out[j,i]
            end

            gammaeff_out[j, i] = gammaeff
        end    
        
    end
    
       
    QLA_P = 1.0
    QLA_E = 1.0
    QLA_O = 1.0 

    phinorm = field_spectrum_out
    # so the normal behavior doesn't change,
    if return_phi_params
        out = Dict(
            (kx_width=kx_width_out),  # [nky] kx_model (kx rms width)
            (gammaeff=gammaeff_out),  # [nky, nmodes] effective growthrate
            (kx0_e=kx0_e),  # [nky] spectral shift in kx
            (ax=ax),  # SAT1 (cx), SAT2 (alpha_x)
            (ay=ay),  # SAT1 (cy)
            (exp_ax=exp_ax),  # SAT2 (sigma_x)
        )
        merge!(out,Dict((d1=d1), (d2=d2), (kycut=kycut), (b3=b3)))
    else
        out = phinorm, QLA_P, QLA_E, QLA_O  # SAT123 intensity and QLA params
    end

    return out
end