function tjlf_max2(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}, ky_s::T, vexb_shear_s::T) where T<:Real

    ### saturation parameters
    R_unit = satParams.R_unit
    q_unit = satParams.q_unit
    ### input parameters
    alpha_p_in = inputs.ALPHA_P
    sign_It_in = inputs.SIGN_IT
    use_bisection_in = inputs.USE_BISECTION
    sat_rule_in = inputs.SAT_RULE
    width_min_in = inputs.WIDTH_MIN
    nwidth_in = inputs.NWIDTH
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)

    width_min = inputs.WIDTH_MIN
    width_max = abs(inputs.WIDTH)


    ### original values
    original_iflux = inputs.IFLUX
    original_ibranch = inputs.IBRANCH
    original_bper = inputs.USE_BPER
    original_bpar = inputs.USE_BPAR
    original_width = inputs.WIDTH
    ### change input values
    inputs.IBRANCH = -1
    inputs.IFLUX = false
    if(sat_rule_in==2 || sat_rule_in==3)
        inputs.USE_BPER = false
        inputs.USE_BPAR = false
    end
    nbasis = ifelse(inputs.NBASIS_MIN!=0, inputs.NBASIS_MIN, inputs.NBASIS_MAX)

    if(alpha_p_in > 0.0)
        for is = ns0:ns
            mass = inputs.MASS[is]
            taus = inputs.TAUS[is]
            vpar_shear = inputs.VPAR_SHEAR[is]
            zs = inputs.ZS[is]
            vs = √(taus/mass)

            kyi = ky_s* √(taus*mass)/abs(zs)
            wgp_max = abs((taus/zs)*alpha_p_in*sign_It_in*vpar_shear
                            /vs)*ky_s/(1+kyi^2)

            width_p_max = 3.6*vs/(√(2)*R_unit*q_unit*max(wgp_max,0.001))
            width_p_max=max(width_p_max,0.1)
            if(width_p_max < width_min_in)
                width_min = width_p_max
            end
        end
    end


    tmax = log10(width_max)
    tmin = log10(width_min)
    nt = nwidth_in
    dtmin = (tmax-tmin)/(nt-1)
    if(use_bisection_in) nt=5 end

    dt = (tmax-tmin)/(nt-1)
    tp = tmin

    gamma_n = zeros(Float64, nt)
    # freq_n = zeros(Float64, nt)
    width_n = zeros(Float64, nt)
    for i = 1:nt
        tp = tmin + (i-1)*dt
        inputs.WIDTH = 10.0^tp
        new_width = true
        # println("this is I")
        nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)

        width_n[i] = inputs.WIDTH
        gamma_n[i] = gamma_out[1]
    end

    ### find the most unstable and save that width
    (gamma_max, imax) = findmax(gamma_n)
    inputs.WIDTH = width_n[imax]

    # use bounded bisection search to refine width
    if(use_bisection_in && gamma_max > 0.0)
        # maximum is against bottom width
        if(imax==1)

            g1 = gamma_n[1]
            t1 = tmin
            g2 = gamma_n[2]
            t2 = log10(width_n[2])
            tp = (t2+t1)/2.0
            inputs.WIDTH = 10.0^tp
            new_width = true
            # println("this is II")
            nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)
            gm = gamma_out[1]
            tm = tp

        # maximum is against top width
        elseif(imax==nt)
            g1 = gamma_n[nt-1]
            t1 = log10(width_n[nt-1])
            g2 = gamma_n[nt]
            t2 = tmax
            tp = (t2+t1)/2.0
            inputs.WIDTH = 10.0^tp
            new_width = true
            # println("this is III")
            nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)

            gm = gamma_out[1]
            tm = tp
        # maximum is away from boundaries
        else
            g1 = gamma_n[imax-1]
            t1 = log10(width_n[imax-1])
            g2 = gamma_n[imax+1]
            t2 = log10(width_n[imax+1])
            gm = gamma_n[imax]
            tm = log10(width_n[imax])
        end

        # start bisection search
        dt=(t2-t1)/2.0
        while(dt > dtmin)
            dt=dt/2.0

            gmax = max(gm,g1,g2)

            if(g1==gmax)
                if(t1>tmin)
                # shift past t1 and compute new g1,t1
                    tp = t1 - dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    # println("this is IV")
                    nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)

                    tm = t1
                    gm = g1
                    g1 = gamma_out[1]
                    t1 = tp

                    # compute new g2,t2
                    tp = tm + dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    # println("this is V")
                    nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)

                    g2 = gamma_out[1]
                    t2 = tp

                else   # t1 at tmin; shrink towards t1
                    tp = t1 + dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    # println("this is VI")
                    tjlf_LS() ############### have to create this ###############
                    g2 = gm
                    t2 = tm
                    gm = gamma_out[1]
                    tm = tp
                end


            elseif(g2==gmax)
                if(t2<tmax) # shift past t2 and compute new g2,t2
                    tp = t2 + dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    # println("this is VII")
                    nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)

                    gm = g2
                    tm = t2
                    g2 = gamma_out[1]
                    t2 = tp

                    # compute new g1,t1
                    tp = tm - dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    # println("this is VIII")
                    nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)

                    g1 = gamma_out[1]
                    t1 = tp
                else  # t2 at tmax shrink towards t2
                    tp = t2 - dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    # println("this is IX")
                    nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)

                    g1 = gm
                    t1 = tm
                    gm = gamma_out[1]
                    tm = tp
                end


            else  # gm==gmax
                # compute new g1,t1 and g2,t2 closer to gm,tm
                tp = tm - dt
                inputs.WIDTH = 10.0^tp
                new_width = true
                # println("this is X")
                nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)

                g1 = gamma_out[1]
                t1 = tp

                tp = tm + dt
                inputs.WIDTH = 10.0^tp
                new_width = true
                # println("this is XI")
                nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)

                g2 = gamma_out[1]
                t2 = tp
            end
        end  # end of bisection search main loop

        # find final maximum
        gmax = gm
        tp = tm
        if(g1>gmax)
            gmax = g1
            tp = t1
        end
        if(g2>gmax)
            gmax = g2
            tp = t2
        end
        gamma_max = gmax
        inputs.WIDTH = 10.0^tp
    end # done with bisection search

    gamma_nb_min_out = gamma_max

    ### reset values (IFLUX determines if QL was calculated)
    inputs.IBRANCH = original_ibranch
    inputs.IFLUX = original_iflux
    if(gamma_max!=0.0) # refine eigenvalue with more basis functions
        ### reset values
        nbasis = inputs.NBASIS_MAX
        if(sat_rule_in==2 || sat_rule_in==3)
            inputs.USE_BPER = original_bper
            inputs.USE_BPAR = original_bpar
        end
        # println("this is XII")
        if !inputs.IFLUX
            nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)
        else
            nmodes_out, gamma_out, freq_out,
            particle_QL_out, 
            energy_QL_out, 
            stress_tor_QL_out, 
            stress_par_QL_out, 
            exchange_QL_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)
        end

        if(inputs.IBRANCH==-1) # check for inward ballooning modes
            if(inputs.USE_INBOARD_DETRAPPED && ft_test > modB_test) ####### find ft_test and modB_test
                error("not implemented XIII")
                fts(:) =  ft_min
                new_geometry = false
                new_width = false
                new_matrix = true
                # println("this is XIII")
                if !inputs.IFLUX
                    nmodes_out, gamma_out, freq_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)
                else
                    nmodes_out, gamma_out, freq_out,
                    particle_QL_out, 
                    energy_QL_out, 
                    stress_tor_QL_out, 
                    stress_par_QL_out, 
                    exchange_QL_out = tjlf_LS(inputs, satParams, outputHermite, ky_s, nbasis, vexb_shear_s)
                end
           end
        end

        gamma_max = max(gamma_out[1],gamma_out[2])  # works for both ibranch_in cases

    end

    if(gamma_max==0.0)
        inputs.WIDTH = original_width
        if(sat_rule_in==2 || sat_rule_in==3)
            inputs.USE_BPER = original_bper
            inputs.USE_BPAR = original_bpar
        end
        maxmodes = 16 #### from tglf_modules
        gamma_out = zeros(Float64,maxmodes)
        freq_out = zeros(Float64,maxmodes)
        particle_QL_out = nothing
        energy_QL_out = nothing
        stress_tor_QL_out = nothing
        stress_par_QL_out = nothing
        exchange_QL_out = nothing
    end


    if !inputs.IFLUX
        return nmodes_out, gamma_nb_min_out, gamma_out, freq_out
    else
        return nmodes_out, gamma_nb_min_out, gamma_out, freq_out,
        particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out
    end
    
    

end