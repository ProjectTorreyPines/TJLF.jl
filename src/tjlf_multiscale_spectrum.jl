using Revise

include("tjlf_modules.jl")
include("intensity_sat_rules.jl")

#
#     parameters: 
#     inputs - InputTJLF struct constructed using the input.TGLF file
#     ky_mix - array of ky (mode number)
#     gamma_mix - array of gamma (net growth rate)
#     
#     outputs:
#     vzf_mix - zonal flow mixing rate
#     kymax_mix - ky value at the calculated vzf_mix
#     jmax_mix - index of ky_mix array where vzf_mix is calculated
#
#     finds the maximum of gamma/ky spectrum at low-k values by going through the ky_mix
#     array, then after finding the max, interpolate the value to improve accuracy
#   
function get_zonal_mixing(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, ky_mix::AbstractVector{T}, gamma_mix::AbstractArray{T}) where T<:Real


    sat_rule_in = inputs.SAT_RULE
    alpha_zf = inputs.ALPHA_ZF

    zs_2 = inputs.ZS[2]
    taus_2 = inputs.TAUS[2]
    mass_2 = inputs.MASS[2]

    rho_ion =  √(taus_2*mass_2) / abs(zs_2)

    # find the local maximum of gamma_mix/ky_mix with the largest gamma_mix/ky_mix^2
    kycut = 0.8/rho_ion
    kymin = 0.0  

    testmax = 0.0
    jmax_mix = 1

    if(alpha_zf < 0.0) kymin = 0.173 * √(2.0) / rho_ion end
    # saturation rules
    if sat_rule_in==2 || sat_rule_in==3
        grad_r0 = satParams.grad_r0
        kycut = grad_r0 * kycut
        kymin = grad_r0 * kymin
    end

    # find the low and high ky peaks of gamma/ky

    # go through ky_mix except for last value
    # if next val is >= kymin and curr val is less than kycut save it
    # update testmax and max index if necessary

    ##### What should i initialize j1 to???
    j1 = nothing
    for j in 1:length(ky_mix)-1
        if((ky_mix[j+1] >= kymin) && (ky_mix[j] <= kycut))
            # save index in case no max
            j1=j
            testmax1 = gamma_mix[j]/ky_mix[j]
            # find maxes
            if(testmax1 > testmax)
                testmax = testmax1
                jmax_mix = j
            end
        end
    end

    ###### What do you do if jmax_mix = NOTHING
    if isnothing(j1)
        error("There is nothing in ky and gamma that meet saturation requirements, IDK what to do -DSun")
    end
    # if testmax is not updated a single time
    if(testmax==0.0) jmax_mix=j1 end

    # no unstable modes in range set kymax index to end of range
    # this is cut of at j1 since a maximum may not exist in the low-k range
    kymax1 = ky_mix[jmax_mix]
    gammamax1 = gamma_mix[jmax_mix]




    # linearly interpolate to find a more accurate low-k maximum gamma/ky
    if(kymax1<kymin)
        kymax1 = kymin
        gammamax1 = (gamma_mix[1] 
            + (gamma_mix[2]-gamma_mix[1])*(kymin-ky_mix[1])/(ky_mix[2]-ky_mix[1]))
    end
    # determine kymax1 and gammamax1 bounded by the tree points f0,f1,f2
    # use a quadratic fit: f = a + b x + c x^2    to f = gamma/ky centered at jmax1
    # scale it to be quadratic where x goes from 0 to 1
    if(jmax_mix>1 && jmax_mix < j1)
        jmax1 = jmax_mix
        f0 = gamma_mix[jmax1-1] / ky_mix[jmax1-1]
        f1 = gamma_mix[jmax1] / ky_mix[jmax1]
        f2 = gamma_mix[jmax1+1] / ky_mix[jmax1+1]
        deltaky = ky_mix[jmax1+1] - ky_mix[jmax1-1]
        x1 = (ky_mix[jmax1] - ky_mix[jmax1-1])/deltaky
        a = f0
        b = (f1 - f0*(1-x1*x1)-f2*x1*x1)/(x1-x1*x1)
        c = f2 - f0 - b

        # if f0>f1 then f1 is not a local maximum
        if(f0 > f1)
            kymax1 = ky_mix[jmax1-1]
            gammamax1 = f0*kymax1

            #interpolate to find the value of gammamax1 at kymin
            if(kymax1 < kymin)
                kymax1 = kymin
                xmin = (kymin - ky_mix[jmax1-1])/deltaky
                gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
            end
        end

        # if f0<f1 then f1>f2 due to the maximum search
        # use the quadratic fit to refine the local maximum:
        if(f0<f1)
            xmax = -b/(2.0*c)
            xmin = 0.0

            if(ky_mix[jmax1-1]<kymin)
                xmin = (kymin - ky_mix[jmax1-1])/deltaky
            end

            # if xmax >= 1    use f2 as the maximum
            if(xmax >= 1.0)
                kymax1 = ky_mix[jmax1+1]
                gammamax1 = f2*kymax1
            elseif(xmax < xmin)   
                # use the quadratic fit to determine gammamax1 at kymin 
                if(xmin > 0.0)
                    kymax1 = kymin
                    gammamax1 = (a + b*xmin + c*xmin*xmin)*kymin
                # if xmax<=0 use f0 as the maximum
                else
                    kymax1 = ky_mix[jmax1-1]
                    gammamax1 = f0*kymax1
                end

            # the conditions f0<f1<f2 and xmin<xmax<1 are satisfied
            # use the quadratic fit to determine gammamax1 and kymax1
            else
                kymax1 = ky_mix[jmax1-1] + deltaky*xmax
                gammamax1 = (a + b*xmax + c*xmax*xmax)*kymax1
            end #xmax tests
        
        end #f0 < f1

    end    # jmax_mix > 1
    vzf_mix = gammamax1/kymax1
    kymax_mix = kymax1
    ### commented out in original Fortran code
    # jmax_mix = jmax1

    return vzf_mix, kymax_mix, jmax_mix

end









#
#     helper function that returns either the smaller value if x is
#     less than, larger value if x is greater than, or a linear 
#     interpolation of the two if x is within the bounds
# 
function mode_transition_function(x::T, y1::T, y2::T, x_ITG::T, x_TEM::T) where T<:Real

    if (x < x_ITG)
        y = y1
    elseif (x > x_TEM)
        y = y2
    else
        y = y1*((x_TEM - x) / (x_TEM - x_ITG)) + y2*((x - x_ITG) / (x_TEM - x_ITG))
    end

    return y
end

#
#     helper function that returns either the linear interpolated y value
#     at a x0 value given the x and y arrays
# 
function linear_interpolation(x::Array{T}, y::Array{T}, x0::T) where T<: Real
    i = 1
    while (x[i] <= x0)
        i += 1
    end
    return ((y[i] - y[i-1]) * x0 + (x[i] * y[i-1] - x[i-1]*y[i])) / (x[i] - x[i-1])

end


######### DSUN thoughts of splitting intensity_sat into 3 functions to help with readability?
######### abandoned for now, but can be seen in the intensity.jl file
# function intensity_sat(inputs::InputTJLF{T}, ky_spect::Vector{T}, gp::Array{T}, QL_data::Array{T}, expsub::T=2.0, return_phi_params::Bool=false) where T<: Real
#     if inputs.SAT_RULE == 1
#         return intensity_sat1(inputs, ky_spect, gp, QL_data)
#     elseif inputs.SAT_RULE == 2
#         return intensity_sat2(inputs, ky_spect, gp, QL_data)
#     elseif inputs.SAT_RULE == 3
#         return intensity_sat3(inputs, ky_spect, gp, QL_data)
#     end
# end


#
#     parameters: 
#     inputs                  - InputTJLF struct constructed using the input.TGLF file
#     ky_spect                - array of ky (mode number)
#     gp                      - array of gamma (net growth rate)
#     QL_data                 - array of quasilinear weights
#     expsub (opt)            - float for the exponent
#     return_phi_params (opt) - boolean flag to change how much is output
#     
#     outputs:
#     phinorm           - intensity of the saturation rule
#     QLA_P             - QLA particle value
#     QLA_E             - QLA energy value
#     QLA_O             - QLA value
#
#     takes in the input.tglf file, ky values, gammas, quasilinear weights, and
#     some other optional parameters and returns the intensity value and QLA
#     parameter values dependent on the 3 saturation rules
#     called get_multiscale_spectrum in the Fortran
# 
function intensity_sat(
    inputs::InputTJLF{T},
    satParams::SaturationParameters{T},
    ky_spect::Vector{T},
    gp::Array{T},
    QL_data::Array{T}, ### taken from the output file
    expsub::T=2.0,
    return_phi_params::Bool=false) where T<:Real

    ############ figure out how to make this prettier
    units_in = inputs.UNITS
    if inputs.SAT_RULE == 2 || inputs.SAT_RULE == 3
        inputs.UNITS = "CGYRO"
        units_in = "CGYRO"
    end

    sat_rule_in = inputs.SAT_RULE
    rlnp_cutoff = inputs.RLNP_CUTOFF
    beta_loc = inputs.BETA_LOC
    ns = inputs.NS
    rmaj_loc = inputs.RMAJ_LOC
    p_prime_loc = inputs.P_PRIME_LOC
    alpha_zf = inputs.ALPHA_ZF
    alpha_quench = inputs.ALPHA_QUENCH

    zs_2 = inputs.ZS[2]
    taus_2 = inputs.TAUS[2]
    mass_2 = inputs.MASS[2]
    taus_1 = inputs.TAUS[1]
    nmodes = inputs.NMODES
    rho_ion = √(taus_2*mass_2)/ abs(zs_2)
    small = 10^-10
    nky = length(ky_spect)

    SAT_geo0_out = satParams.SAT_geo0
    SAT_geo1_out = satParams.SAT_geo1
    SAT_geo2_out = satParams.SAT_geo2
    Bt0_out = satParams.Bt0
    b_geo0_out = satParams.B_geo0
    grad_r0_out = satParams.grad_r0
    kx0_e = xgrid_functions_geo(inputs, satParams, ky_spect,  gp)


    if length(size(gp)) > 1
        gammas1 = gp[:, 1] # SAT1 and SAT2 use the growth rates of the most unstable modes
    else
        gammas1 = gp
    end

    vzf_out, kymax_out, jmax_out = get_zonal_mixing(inputs, satParams, ky_spect, gammas1)
    
    # model fit parameters
    # Miller geometry values igeo=1
    if(rlnp_cutoff > 0.0)
        if(beta_loc == 0.0)
            dlnpdr = 0.0
            ptot = 0.0
            for i in 1:ns
                ptot = ptot + inputs.AS[i]*inputs.TAUS[i]
                dlnpdr = (dlnpdr + 
                    inputs.AS[i]*inputs.TAUS[i]*
                    (inputs.RLNS[i]+inputs.RLTS[i]))
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


    ### coefficients for SAT_RULE = 1
    czf = abs(alpha_zf)
    cnorm = 14.29
    cz1=0.48*czf
    cz2=1.0*czf
    cky=3.0
    sqcky=√(cky)
    cnorm = 14.29
    etg_streamer=1.05
    kyetg = etg_streamer * abs(zs_2) / √(taus_2 * mass_2)
    
    measure = √(taus_1 * mass_2)
    # if(USE_SUB1)
    #     cnorm=12.12
    #     expsub=1
    # end
    if(alpha_quench != 0.0) etg_streamer=2.1 end
    # if(igeo == 0)then # s-alpha
    #     cnorm=14.63
    #     cz1=0.90*czf
    #     cz2=1.0*czf
    # end


    ### coefficents for SAT_RULE = 2
    if(sat_rule_in == 2 || sat_rule_in == 3)
        b0 = 0.76
        b1 = 1.22
        b2 = 3.74
        if(nmodes > 1) b2 = 3.55 end
        b3 = 1.0

        d1 = (Bt0_out/b_geo0_out)^4 
        d1 = d1/grad_r0_out
        Gq = b_geo0_out/grad_r0_out
        d2 = b3/(Gq^2)
        cnorm = b2*(12.0/dlnpdr)
        kyetg = 1000.0
        cky = 3.0
        sqcky = √(cky)
        kycut = b0*kymax_out
        cz1 = 0.0
        cz2 = 1.05*czf
        measure = 1.0/kymax_out
    end
    ### coefficents for SAT_RULE = 3
    if(sat_rule_in==3)
        kmax = kymax_out
        gmax = vzf_out * kymax_out
        kmin = 0.685 * kmax
        aoverb = - 1.0 / (2 * kmin)
        coverb = - 0.751 * kmax
        kT = 1.0/rho_ion
        k0 = 0.6 * kmin
        kP = 2.0 * kmin
        c_1 = - 2.42
        x_ITG = 0.8
        x_TEM = 1.0
        Y_ITG = 3.3 * (gmax^2) / (kmax^5)
        Y_TEM = 12.7 * (gmax^2) / (kmax^4)
        scal = 0.82 # Q(SAT3 GA D) / (2 * QLA(ITG,Q) * Q(SAT2 GA D))
        
        Ys = zeros(nmodes)
        xs = zeros(nmodes)

        for k in 1:nmodes
            sum_W_i = zeros(length(QL_data[:, k, 1, 1, 2]))
            for is in 2:size(QL_data)[3] # sum over ion species, requires electrons to be species 1
                sum_W_i .= sum_W_i .+ QL_data[:, k, is, 1, 2]
            end
            # check for singularities in weight ratio near kmax
            ### isn't i == jmax_out? I don't think you need the loop
            i = 1
            while (ky_spect[i] < kmax)
	            i += 1
	        end

            ### why is sum_W_i an array now???
            if(sum_W_i[i]==0.0 || sum_W_i[i-1]==0.0)
                x = 0.5
            else
                abs_W_ratio = abs.(QL_data[:,k,1,1,2]./sum_W_i)
                x = linear_interpolation(ky_spect, abs_W_ratio, kmax)
            end
            Y = mode_transition_function(x, Y_ITG, Y_TEM, x_ITG, x_TEM)

            xs[k] = x
            Ys[k] = Y
        end
    end


    # coefficients for spectral shift model for ExB shear
    ax=0.0
    ay=0.0
    exp_ax = 1
    if(alpha_quench == 0.0)
        #spectral shift model parameters
        ax = 1.15
        ay = 0.56
        exp_ax = 4

        if(sat_rule_in==2 || sat_rule_in==3)
            ax = 1.21
            ay = 1.0
            exp_ax = 2
        end
    end

    # if(!first_pass)  # second pass for spectral shift model
    gamma_net = zeros(nky)
    for i in 1:nky
        kx = kx0_e[i]
        if(sat_rule_in==2 || sat_rule_in==3)
            ky0 = ky_spect[i]
            if(ky0<kycut)
                kx_width = kycut/grad_r0_out
            else
                kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
            end
            kx = kx*ky0/kx_width
        end
        gamma_net[i] = gammas1[i]/(1.0 + abs(ax*kx)^exp_ax)
    end

    if(sat_rule_in==1)
        vzf_out, kymax_out, jmax_out = get_zonal_mixing(inputs, satParams, ky_spect, gamma_net)
    else
        vzf_out_fp = vzf_out # used for recreation of first pass for SAT3
        vzf_out = vzf_out*gamma_net[jmax_out]/max(gammas1[jmax_out], small)
    end
    # end   # second pass complete

    gammamax1= vzf_out*kymax_out
    kymax1 = kymax_out
    jmax1 = jmax_out
    vzf1 = vzf_out



    gamma_mix1 = zeros(nky)
    gamma = zeros(nky)

    if(sat_rule_in==1)
        gamma .= ifelse.(ky_spect.<kymax1, 
                    max.(gamma_net .- (cz1.*(kymax1 .- ky_spect).*vzf1),0.0),
                    (cz2*gammamax1) .+ max.(gamma_net .- (cz2.*vzf1.*ky_spect),0.0))
    elseif(sat_rule_in==2 || sat_rule_in==3)
        gamma .= ifelse.(ky_spect.<kymax1,
                    gamma_net,
                    gammamax1 .+ max.(gamma_net .- cz2*vzf1.*ky_spect, 0.0))
    end
    gamma_mix1 .= gamma


    # if(USE_MIX)
        # mix over ky > kymax with integration weight = sqcky*ky0**2/(ky0**2 + cky*(ky-ky0)**2)
        #### unsure of this index
        mixnorm1 = 0.0
        for j in jmax1+2:nky
            gamma_ave = 0.0
            ky0 = ky_spect[j]
            mixnorm1 = (ky0*
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
            gamma_mix1[j] = gamma_ave/mixnorm1
        end
    # end
    
    # recreate first pass for ExB rule in SAT3
    if(sat_rule_in==3)
        # if(!first_pass)
            gamma_fp = zeros(size(ky_spect))
            gamma = zeros(size(ky_spect))

            gamma .= ifelse.(ky_spect.<kymax1,
                            gammas1, 
                            ((gammamax1 * (vzf_out_fp/vzf_out)) .+ max.(gammas1 .- (cz2.*vzf_out_fp.*ky_spect),0.0)))
        
            gamma_fp .= gamma
        # end


        # if(USE_MIX)
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
                            (atan(sqcky*(ky2/ky0-1.0))
                            - atan(sqcky*(ky1/ky0-1.0)))
                    )
                    delta = (gamma[i+1]-gamma[i])/(ky2-ky1)
                    mix2 = (ky0*mix1 
                            + (ky0*ky0/(2.0*sqcky))*(log(cky*(ky2-ky0)^2+ky0^2)
                            - log(cky*(ky1-ky0)^2+ky0^2))
                    )
                    gamma_ave = gamma_ave + (gamma[i] - ky1*delta)*mix1 + delta*mix2
                end
                gamma_fp[j] = gamma_ave/mixnorm
            end  
        # end
    # else
    #     gamma_fp = gamma_mix1
    end
    


    # generate SAT2 potential for SAT3 to connect to for electron scale

    YTs = zeros(nmodes)
    if(sat_rule_in==3) 
        if(ky_spect[nky] >= kT)
            dummy_interp = zeros(size(ky_spect))
            k = 1
            while ky_spect[k] < kT
                k += 1
            end
            for l in k-1:k
                gamma0 = gammas1[l]
	            ky0 = ky_spect[l]
	            kx = kx0_e[l]

                if(ky0 < kycut)
                    kx_width = kycut/grad_r0_out
                    sat_geo_factor = SAT_geo0_out*d1*SAT_geo1_out
                else
                    kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
                    sat_geo_factor = (SAT_geo0_out*
                                (d1*SAT_geo1_out*kycut +
                                (ky0 - kycut)*d2*SAT_geo2_out)/ky0
                    )
                end
	            kx = kx*ky0/kx_width
	            gammaeff = 0.0
	            if(gamma0 > small) gammaeff = gamma_fp[l] end
	            # potentials without multimode and ExB effects, added later
	            dummy_interp[l] = scal*measure*cnorm*(gammaeff/(kx_width*ky0))^2
                if(units_in != "GYRO") dummy_interp[l] = sat_geo_factor*dummy_interp[l] end
            end
	        YT = linear_interpolation(ky_spect, dummy_interp, kT)
	        YTs .= YT
	    else
	        if(aoverb*(kP^2)+kP+coverb-((kP-kT)*(2*aoverb*kP+1))!=0)
	            for l in 1:nmodes
	                YTs[l] = (Ys[l]*
                            (
                                ((aoverb*(k0^2) + k0+coverb)
                                /(
                                    aoverb*(kP^2) 
                                    + kP
                                    + coverb
                                    -((kP - kT)*(2*aoverb*kP+1))
                                    )
                                )^abs(c_1))
                    )
                end
            end
        end
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
        if(sat_rule_in==1)
            sat_geo_factor = SAT_geo0_out
            kx_width = ky0
        elseif(sat_rule_in==2 || sat_rule_in==3)
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
        end
        kx_width_out[j] = kx_width
        sat_geo_factor_out[j] = sat_geo_factor

        if(sat_rule_in==1 || sat_rule_in==2)
            for i in 1:nmodes
                gammaeff = 0.0
                if(gamma0>small)
                    gammaeff = gamma_mix1[j]*(gp[j,i]/gamma0)^expsub
                end
                if(ky0>kyetg) gammaeff = gammaeff * √(ky0/kyetg) end

                field_spectrum_out[j,i] = (measure*cnorm*
                                ((gammaeff/(kx_width*ky0))/(1.0+ay*kx^2))^2)
                if(units_in != "GYRO")
                    field_spectrum_out[j,i] = sat_geo_factor*field_spectrum_out[j,i]
                end
                gammaeff_out[j, i] = gammaeff
            end
    
        elseif(sat_rule_in==3) # SAT3
            if(gamma_fp[j]==0)
                Fky=0.0
            else
		        Fky =  ( (gamma_mix1[j] / gamma_fp[j])^2 ) / ((1.0 + ay*(kx^2))^2)
            end
            for i in 1:nmodes
                field_spectrum_out[j,i] = 0.0  
                gammaeff = 0.0
                if(gamma0 > small)
                    if (ky0 <= kP) # initial quadratic

			            sig_ratio = (aoverb * (ky0^2) + ky0 + coverb) / (aoverb * (k0^2) + k0 + coverb)	
			            field_spectrum_out[j,i] = Ys[i] * (sig_ratio^c_1) * Fky * (gp[j,i]/gamma0)^(2 * expsub)

                    elseif (ky0 <= kT) # connecting quadratic

		                if(YTs[i]==0.0 || kP==kT) 
                            field_spectrum_out[j,i] = 0.0
		                else
                            doversig0 = ( (Ys[i] / YTs[i]) ^ (1.0/abs(c_1))
                                        - (aoverb*kP^2 + kP + coverb - (kP-kT)*(2*aoverb*kP + 1))
                                        /(aoverb*k0^2 + k0 + coverb))
                            doversig0 = doversig0 * (1.0/((kP-kT)^2))
                            eoversig0 = -2 * doversig0 * kP + ((2 * aoverb * kP + 1)/(aoverb * (k0^2) + k0 + coverb))
                            foversig0 = ((Ys[i] / YTs[i])^(1.0/abs(c_1))) - eoversig0 * kT - doversig0 * (kT^2)
                            sig_ratio = doversig0*(ky0^2) + eoversig0*ky0 + foversig0
                            field_spectrum_out[j,i] = Ys[i] * (sig_ratio ^ c_1) * Fky * (gp[j,i]/gamma0)^(2*expsub)
                        end

                    else # SAT2 for electron scale

                        gammaeff = gamma_mix1[j]*(gp[j,i]/gamma0)^expsub
                        if(ky0 > kyetg) gammaeff = gammaeff*√(ky0/kyetg) end

                        field_spectrum_out[j,i] = (scal*measure*cnorm*
                                        ((gammaeff/(kx_width*ky0))/(1.0+ay*kx^2))^2)
                        if(units_in != "GYRO")
                            field_spectrum_out[j,i] = sat_geo_factor*field_spectrum_out[j,i]
                        end
                        
                    end
                end
                gammaeff_out[j, i] = gammaeff
            end
        end
    end
    
	#SAT3 QLA here
    if(sat_rule_in==3)
        QLA_P = zeros(nmodes)
	    QLA_E = zeros(nmodes)
	    for k in 1:nmodes
	        # factor of 2 included for real symmetry
            QLA_P[k] = 2 * mode_transition_function(xs[k], 1.1, 0.6, x_ITG, x_TEM)
            QLA_E[k] = 2 * mode_transition_function(xs[k], .75, 0.6, x_ITG, x_TEM)
        end
	    QLA_O = 2 * 0.8
    else
	   QLA_P = 1.0
	   QLA_E = 1.0
	   QLA_O = 1.0 
    end

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
        if sat_rule_in == 2
            # add bonus geometry params,
            merge!(out,Dict((d1=d1), (d2=d2), (kycut=kycut), (b3=b3)))
        end
    else
        out = phinorm, QLA_P, QLA_E, QLA_O  # SAT123 intensity and QLA params
    end

    return out
		
end

#
#     helper function that calculates the flux integrals given the index,
#     and ky/QL values among other things
# 
######### i think dky0 is sometimes a different type for some reason
function flux_integrals(
    i::Int,
    ky::T,
    dky0::K,
    dky1::T,
    particle::Array{T},
    energy::Array{T},
    toroidal_stress::Array{T},
    parallel_stress::Array{T},
    exchange::Array{T},
    particle_flux_out,
    energy_flux_out,
    stress_tor_out,
    stress_par_out,
    exchange_out,
    q_low_out,
    taus_1=1.0,
    mass_2=1.0,
) where {T <: Real, K<: Real}
    """
    Compute the flux integrals
    """
    particle_flux_out .+= (
        dky0 .* (i==1 ? 0 : particle[i - 1,:,:,:])
        .+ dky1 .* particle[i,:,:,:]
    )
    energy_flux_out .+= (
        dky0 .* (i==1 ? 0 : energy[i - 1,:,:,:])
        .+ dky1 .* energy[i,:,:,:]
    )
    stress_tor_out .+= (
        dky0 .* (i==1 ? 0 : toroidal_stress[i - 1,:,:,:])
        .+ dky1 .* toroidal_stress[i,:,:,:]
    )
    stress_par_out .+= (
        dky0 .* (i==1 ? 0 : parallel_stress[i - 1,:,:,:])
        .+ dky1 .* parallel_stress[i,:,:,:]
    )
    exchange_out .+= (
        dky0 .* (i==1 ? 0 : exchange[i - 1,:,:,:])
        .+ dky1 .* exchange[i,:,:,:]
    )

    if ky * taus_1 * mass_2 <= 1
        q_low_out .= (
            energy_flux_out[:,:,1] + energy_flux_out[:,:,2]
        )
    end
    return (
        particle_flux_out,
        energy_flux_out,
        stress_tor_out,
        stress_par_out,
        exchange_out,
        q_low_out,
    )
end




#
#     parameters: 
#     inputs            - InputTJLF struct constructed using the input.TGLF file
#     ky_spect          - array of ky (mode number)
#     gp                - array of gamma (net growth rate)
#     QL_data           - split into separate types of QL weights
#     (optional)        - a lot of optional parameters that I don't use -DSUN
#     
#     outputs:
#     a dictionary with the flux integral values:
#           "particle_flux_integral" =>particle_flux_out,
#           "energy_flux_integral" => energy_flux_out,
#           "toroidal_stresses_integral" => stress_tor_out,
#           "parallel_stresses_integral" => stress_par_out,
#           "exchange_flux_integral" => exchange_out,
#
#     takes in the input.tglf file, ky values, gammas, quasilinear weights, calls
#     intensity_sat() which returns the QL intensity values, then calls
#     flux_integrals() which numerically integrates the fluxes
# 

function sum_ky_spectrum(
    inputs::InputTJLF{T},
    satParams::SaturationParameters{T},
    ky_spect::Vector{T},
    gp::Array{T},
    ave_p0::Vector{T},
    potential::Array{T},
    particle_QL::Array{T},
    energy_QL::Array{T},
    toroidal_stress_QL::Array{T},
    parallel_stress_QL::Array{T},
    exchange_QL::Array{T},
    etg_fact::T=1.25,
    c0::T=32.48,
    c1::T=0.534,
    exp1::T=1.547,
    cx_cy::T=0.56,
    alpha_x::T=1.15,
)where T <: Real

    sat_rule_in = inputs.SAT_RULE

    NM = length(energy_QL[1, :, 1, 1])  # get the number of modes
    NS = length(energy_QL[1, 1, :, 1])  # get the number of species
    NF = length(energy_QL[1, 1, 1, :])  # get the number of fields
    particle_flux_out = zeros((NM, NS, NF))
    energy_flux_out = zeros((NM, NS, NF))
    stress_tor_out = zeros((NM, NS, NF))
    stress_par_out = zeros((NM, NS, NF))
    exchange_out = zeros((NM, NS, NF))
    q_low_out = zeros((NM, NS))

    
    QL_data = cat(
        dims=5, particle_QL, energy_QL, toroidal_stress_QL, parallel_stress_QL, exchange_QL
    )
    # Multiply QL weights with desired intensity
    if sat_rule_in in [1.0, 1, "SAT1", 2.0, 2, "SAT2", 3.0, 3, "SAT3"]
        intensity_factor, QLA_P, QLA_E, QLA_O = intensity_sat(
            inputs, satParams, ky_spect, gp, QL_data
        )
    else
        throw(error(
            "sat_rule_in must be [1.0, 1, 'SAT1', 2.0, 2, 'SAT2', 3.0, 3, 'SAT3], not $sat_rule_in"
            )
        )
    end

    if size(particle_QL) == size(energy_QL) == size(toroidal_stress_QL) == size(parallel_stress_QL) == size(exchange_QL)
        shapes = size(particle_QL)
    else
        error("QL matrices are different sizes!")
    end

    particle = zeros(shapes)
    energy = zeros(shapes)
    toroidal_stress = zeros(shapes)
    parallel_stress = zeros(shapes)
    exchange = zeros(shapes)


    particle .= particle_QL .* (intensity_factor .* QLA_P')
    energy .= energy_QL .* (intensity_factor .* QLA_E')
    toroidal_stress .= toroidal_stress_QL .* (intensity_factor .* QLA_O)
    parallel_stress .= parallel_stress_QL .* (intensity_factor .* QLA_O)
    exchange .= exchange_QL .* (intensity_factor .* QLA_O)

    
    dky0 = 0
    ky0 = 0
    for i in eachindex(ky_spect)
        ky = ky_spect[i]
        ky1 = ky
        if i == 1
            dky1 = ky1
        else
            dky = log(ky1 / ky0) / (ky1 - ky0)
            dky1 = ky1 * (1.0 - ky0 * dky)
            dky0 = ky0 * (ky1 * dky - 1.0)
        end

        (
            particle_flux_out,
            energy_flux_out,
            stress_tor_out,
            stress_par_out,
            exchange_out,
            q_low_out,
        ) = flux_integrals(
            i,
            ky,
            dky0,
            dky1,
            particle,
            energy,
            toroidal_stress,
            parallel_stress,
            exchange,
            particle_flux_out,
            energy_flux_out,
            stress_tor_out,
            stress_par_out,
            exchange_out,
            q_low_out,
        )
        ky0 = ky1
    end

    results = Dict(
            "particle_flux_integral" =>particle_flux_out,
            "energy_flux_integral" => energy_flux_out,
            "toroidal_stresses_integral" => stress_tor_out,
            "parallel_stresses_integral" => stress_par_out,
            "exchange_flux_integral" => exchange_out,
    )

    return results

end