include("tjlf_modules.jl")
include("tjlf_geometry.jl")
include("tjlf_LS.jl")

function tjlf_max(inputs::InputTJLF, ky_s::T, vexb_shear_s::T) where T<:Real

    R_unit, q_unit= get_sat_params(:rq_units, inputs)
    ### basically calls xgrid_functions_geo() to get R_unit, q_unit, make sure this is not redundant! maybe create output struct

    #### organize this better
    save_iflux = inputs.IFLUX
    save_ibranch = inputs.IBRANCH
    save_nbasis = inputs.NBASIS_MAX
    save_bper = inputs.USE_BPER
    save_bpar = inputs.USE_BPAR
    alpha_p_in = inputs.ALPHA_P
    sign_It_in = inputs.SIGN_IT
    use_bisection_in = inputs.USE_BISECTION

    sat_rule_in = inputs.SAT_RULE
    width_min_in = inputs.WIDTH_MIN
    nbasis_min_in = inputs.NBASIS_MIN
    nwidth_in = inputs.NWIDTH
    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end

    save_width = inputs.WIDTH
    width_min = inputs.WIDTH_MIN
    width_max = abs(inputs.WIDTH)

    inputs.IBRANCH = -1

    if(alpha_p_in > 0.0)
        for is = ns0:ns
            mass = inputs.SPECIES[is].MASS
            taus = inputs.SPECIES[is].TAUS
            vpar_shear = inputs.SPECIES[is].VPAR_SHEAR
            zs = inputs.SPECIES[is].ZS
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

    inputs.IFLUX = false
    if(nbasis_min_in!=0)
        nbasis = nbasis_min_in
    end
    if(sat_rule_in==2 || sat_rule_in==3)
        inputs.USE_BPER = false
        inputs.USE_BPAR = false
    end

    tmax = log10(width_max)
    tmin = log10(width_min)
    nt = nwidth_in
    dtmin = (tmax-tmin)/(nt-1)
    if(use_bisection_in) nt=5 end

    dt = (tmax-tmin)/(nt-1)
    tp = tmin

    gamma_n = zeros(nt)
    freq_n = zeros(nt)
    width_n = zeros(nt)
    for i = 1:nt
        tp = tmin + (i-1)*dt
        inputs.WIDTH = 10.0^tp
        new_width = true
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
        Ns_Ts_phase_out = tjlf_LS(inputs, ky_s, vexb_shear_s) ############### have to create this ###############
        ############### these are prob outputs ###############
        width_n[i] = inputs.WIDTH
        gamma_n[i] = gamma_out[1]
        freq_n[i] = freq_out[1]
    end

    ### this might be off by one if there are repeats, og code finds last max value
    ### this function one finds first max
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
            tglf_LS() ############### have to create this ###############
            ############### these are prob outputs ###############
            gm = gamma_out[1]
            tm = tp

        # maximum is against top width
        elseif(imax==nt)
            g1 = gamma_n[nt-1]
            t1 = log10(width_n[nt-1])
            g2 = gamma_n[nt]
            t2 = tmax           
            tp = (t2+t1)/2.0    
            width_in = 10.0^tp
            new_width = true
            tglf_LS() ############### have to create this ###############
            ############### these are prob outputs ###############
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
            ###### rewrite this, it is so ugly
            if(g1==gmax)
                if(t1>tmin)
                # shift past t1 and compute new g1,t1 
                    tp = t1 - dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    tglf_LS() ############### have to create this ###############
                    ############### these are prob outputs ###############
                    tm = t1
                    gm = g1
                    g1 = gamma_out[1]
                    t1 = tp

                    # compute new g2,t2
                    tp = tm + dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    tglf_LS() ############### have to create this ###############
                    ############### these are prob outputs ###############
                    g2 = gamma_out[1]
                    t2 = tp

                else   # t1 at tmin; shrink towards t1
                    tp = t1 + dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    tglf_LS() ############### have to create this ###############
                    ############### these are prob outputs ###############
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
                    tglf_LS() ############### have to create this ###############
                    ############### these are prob outputs ###############
                    gm = g2
                    tm = t2
                    g2 = gamma_out[1]
                    t2 = tp

                    # compute new g1,t1
                    tp = tm - dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    tglf_LS() ############### have to create this ###############
                    ############### these are prob outputs ###############
                    g1 = gamma_out[1]
                    t1 = tp
                else  # t2 at tmax shrink towards t2
                    tp = t2 - dt
                    inputs.WIDTH = 10.0^tp
                    new_width = true
                    tglf_LS() ############### have to create this ###############
                    ############### these are prob outputs ###############
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
                tglf_LS() ############### have to create this ###############
                ############### these are prob outputs ###############
                g1 = gamma_out[1]
                t1 = tp

                tp = tm + dt
                inputs.WIDTH = 10.0^tp
                new_width = true
                tglf_LS() ############### have to create this ###############
                ############### these are prob outputs ###############
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



    ### output value??
    gamma_nb_min_out = gamma_max
    #  reset ibranch_in
    inputs.IBRANCH = save_ibranch

    if(gamma_max!=0.0) # refine eigenvalue with more basis functions
        
        nbasis = save_nbasis
        if(sat_rule_in==2 || sat_rule_in==3)
            inputs.USE_BPER = save_bper
            inputs.USE_BPAR = save_bpar
        end
         
        inputs.IFLUX = save_iflux
        new_width = true
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
        Ns_Ts_phase_out = tjlf_LS(inputs, ky_s, vexb_shear_s) ############### have to create this ###############

        if(inputs.IBRANCH==-1) # check for inward ballooning modes
            if(inputs.USE_INBOARD_DETRAPPED && ft_test > modB_test) ####### find ft_test and modB_test
                fts(:) =  ft_min
                new_geometry = false
                new_width = false
                new_matrix = true
                tjlf_LS() ############### have to create this ###############
           end
        end
         
        gamma_max = max(gamma_out[1],gamma_out[2])  # works for both ibranch_in cases

    end

    if(gamma_max==0.0)
        inputs.WIDTH = save_width
        if(sat_rule_in==2 || sat_rule_in==3)
            inputs.USE_BPER = save_bper
            inputs.USE_BPAR = save_bpar
        end
        gamma_out = zeros(Float64,inputs.NMODES)
        freq_out = zeros(Float64,inputs.NMODES)
    end

    nbasis = save_nbasis ####### worrisome, is this an output? 
    inputs.IFLUX = save_iflux

    return gamma_nb_min_out,
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
            Ns_Ts_phase_out

end