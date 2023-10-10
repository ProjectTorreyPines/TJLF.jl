import LinearAlgebra.LAPACK.gesv!

include("tjlf_modules.jl")
include("tjlf_geometry.jl")
include("tjlf_eigensolver.jl")
include("tjlf_matrix.jl")
#***********************************************************************
#  TGLF Linear Stability driver computes all quasilinear quantities 
#  for a single ky
#
#***********************************************************************
function tjlf_LS(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}, 
    ky::T) where T <: Real
    nbasis = ifelse(inputs.NBASIS_MIN!=0, inputs.NBASIS_MIN, inputs.NBASIS_MAX)
    width_in = inputs.WIDTH
    iflux_in = inputs.IFLUX
    ibranch_in = inputs.IBRANCH
    use_bper_in = inputs.USE_BPER
    use_bpar_in = inputs.USE_BPAR
    vexb_shear_s = inputs.VEXB_SHEAR*inputs.SIGN_IT
    return tjlf_LS(inputs,satParams,outputHermite,
                    ky, nbasis,
                    width_in,iflux_in,ibranch_in,use_bper_in,use_bpar_in,vexb_shear_s)
end
function tjlf_LS(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outputHermite::OutputHermite{T}, 
            ky::T, 
            nbasis::Int, width_parameter::T,iflux_parameter::Bool,ibranch_parameter::Int,
            use_bper_parameter::Bool,use_bpar_parameter::Bool,
            vexb_shear_s::T) where T <: Real

    small = 1.0e-13
    epsilon1 = 1.0e-12
    nmodes_in = inputs.NMODES
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    alpha_quench_in = inputs.ALPHA_QUENCH
    new_eikonal_in = inputs.NEW_EIKONAL
    

    R_unit = satParams.R_unit
    q_unit = satParams.q_unit

    new_geometry = true ###### hardcoded for now
    # check co-dependencies
    if(new_geometry) new_width= true end 
    if(new_width) new_matrix = true end

    if(new_eikonal_in)
        if(new_geometry)
            # MILLER geometry
            igeo = 1
            if(igeo==1)
                #### i think this is for like debugging? like you can look at trace_path to help figure out what ran
                # trace_path[4]=1
                # miller_geo(inputs)
                # mercier_luc(inputs)

            # FOURIER geometry
            elseif(igeo==2)
                error("sorry fourier geometry not implemented yet :(")
                # trace_path[5]=1
                fourier_geo(inputs)

            # ELITE geometry
            elseif(igeo==3)
                error("sorry fourier geometry not implemented yet :(")
                #trace_path[8]=1
                ELITE_geo(inputs)
            end
        
            # compute the eikonal functions for general geometry (igeo>0)
            # if(igeo > 0) mercier_luc(inputs) end

        end
        new_geometry = false


        #  load the x-grid eikonal functions v_QL_out,b0x
        if(new_width)
            # trace_path[6]=1
            # get_xgrid_functions() ############### have to create this ###############
            outputGeo = xgrid_functions_geo(inputs, outputHermite, ky, width_parameter)
            R_unit, q_unit = get_sat_params(:rq_units,inputs) ### no work if igeo is 0
        end
    end  #new_eikonal_in


    # compute the rank of the eigenmatrix iur
    nroot=15
    iur = (ns-ns0+1)*nroot*nbasis
    new_matrix = true ####### hard code for now
    if(new_matrix)
        # trace_path[7]=1
        ave, aveH, aveWH, aveKH, aveG, aveWG, aveKG = get_matrix(inputs, outputGeo, outputHermite, ky, width_parameter) ############### have to create this ###############
    end

    #  solver for linear eigenmodes of tglf equations
    amat, bmat, alpha, beta, rr, ri = tjlf_eigensolver(inputs,outputGeo,satParams,ave,aveH,aveWH,aveKH,aveG,aveWG,aveKG,ky,width_parameter)
    # println(rr)
    # println(ri)
    # println(alpha)
    # println(beta)
    # println(amat)
    # println(bmat)

    #  initalize output to zero
    maxmodes = 16 #### from tglf_modules
    jmax = zeros(Int, maxmodes)
    gamma_out = zeros(Float64, maxmodes)
    freq_out = zeros(Float64, maxmodes)
    v_QL_out = zeros(Float64, maxmodes)
    a_par_QL_out = zeros(Float64, maxmodes)
    b_par_QL_out = zeros(Float64, maxmodes)
    phi_bar_out = zeros(Float64, maxmodes)
    a_par_bar_out = zeros(Float64, maxmodes)
    b_par_bar_out = zeros(Float64, maxmodes)
    v_bar_out = zeros(Float64, maxmodes)
    ne_te_phase_out = zeros(Float64, maxmodes)

    field_weight_out = zeros(ComplexF64, maxmodes, 3, nbasis)
    particle_QL_out = zeros(Float64, maxmodes, ns, 3)
    energy_QL_out = zeros(Float64, maxmodes, ns, 3)
    stress_par_QL_out = zeros(Float64, maxmodes, ns, 3)
    stress_tor_QL_out = zeros(Float64, maxmodes, ns, 3)
    exchange_QL_out = zeros(Float64, maxmodes, ns, 3)

    N_QL_out = zeros(Float64, maxmodes, ns)
    T_QL_out = zeros(Float64, maxmodes, ns)
    U_QL_out = zeros(Float64, maxmodes, ns)
    Q_QL_out = zeros(Float64, maxmodes, ns)
    N_bar_out = zeros(Float64, maxmodes, ns)
    T_bar_out = zeros(Float64, maxmodes, ns)
    U_bar_out = zeros(Float64, maxmodes, ns)
    Q_bar_out = zeros(Float64, maxmodes, ns)
    Ns_Ts_phase_out = zeros(Float64, maxmodes, ns)


    if(ibranch_parameter==0)
        di = zeros(Int, iur)
        de = zeros(Int, iur)
        nmodes_out = nmodes_in
        # sort the unstable modes into electron and ion frequencies
        mi = 0
        me = 0
        for j1 = 1:iur
            if(rr[j1]>epsilon1)
                if(ri[j1]>0.0)
                    # note that ri = -freq, rr = gamma
                    mi = mi+1
                    di[mi] = j1
                else
                    me = me+1
                    de[me] = j1
                end
            end
        end


        # find the most unstable mode for each branch
        if(me>0)
            zgamax = 0.0
            for iroot = 1:me
                if(rr[de[iroot]]>zgamax)then
                    zgamax = rr[de[iroot]]
                    jmax[1] = de[iroot]
                end
            end
            gamma_out[1] = rr[jmax[1]]
            freq_out[1] = -ri[jmax[1]]
        end
        if(mi>0)
            zgamax = 0.0
            for iroot = 1:mi
                if(rr[di[iroot]]>zgamax)
                    zgamax = rr[di[iroot]]
                    jmax[2] = di[iroot]
                end
            end
            gamma_out[2] = rr[jmax[2]]
            freq_out[2] = -ri[jmax[2]]
        end


    elseif(ibranch_parameter==-1)
        # find the top nmodes most unstable modes
        ### put the unstable modes in ascending order by growthrate
        sortperm
        jmax = sortperm(rr) #### sort_eigenvalues(nmodes_in,jmax), i believe this does the same thing
        jmax .= ifelse.(rr[jmax].>epsilon1, jmax, 0)
        reverse!(jmax)
        nmodes_out = 0
        for j1 = 1:nmodes_in 
            if(jmax[j1]!=0)
                nmodes_out = nmodes_out + 1
                gamma_out[j1] = rr[jmax[j1]]
                freq_out[j1] = -ri[jmax[j1]]
            end
        end
    end



    # apply quench rule
    if(alpha_quench_in!=0.0)
        for j1 = 1:nmodes_in
            gamma_out[j1] = get_gamma_net(inputs, gamma_out[j1],vexb_shear_s)
        end
        
    # use spectral shift model for second pass
    elseif(vexb_shear_s!=0.0)
        for j1 = 1:nmodes_in
            gamma_out[j1] = gamma_reference_kx0[j1] ####### pass this in as a parameter ########
            freq_out[j1] = freq_reference_kx0[j1] ####### pass this in as a parameter ########
        end
    end
    
    # get the fluxes for the most unstable modes
    v = zeros(ComplexF64, iur)
    zmat = Matrix{ComplexF64}(undef, iur, iur)
    if(iflux_parameter)
        wd_bar_out = Vector{Float64}(undef, nmodes_out)
        b0_bar_out = Vector{Float64}(undef, nmodes_out)
        modB_bar_out = Vector{Float64}(undef, nmodes_out)
        v_QL_out = Vector{Float64}(undef, nmodes_out)
        a_par_QL_out = Vector{Float64}(undef, nmodes_out)
        b_par_QL_out = Vector{Float64}(undef, nmodes_out)
        kx_bar_out = Vector{Float64}(undef, nmodes_out)
        kpar_bar_out = Vector{Float64}(undef, nmodes_out)
        for imax = 1:nmodes_out
            if(jmax[imax]>0)
                for i = 1:iur
                    v[i] = small ##### what is v
                    for j = 1:iur
                        zmat[i,j] = beta[jmax[imax]]*amat[i,j] - (small +alpha[jmax[imax]])*bmat[i,j]
                    end
                end
                ### gesv!(A,B) solves Ax = B, A becomes LU factor, and B becomes solution
                gesv!(zmat,v) 
                # println(v)

                eigenvalue = im*alpha[jmax[imax]]/beta[jmax[imax]]
                
                Ns_Ts_phase,
                Ne_Te_phase,
                N_weight, 
                T_weight, 
                U_weight, 
                Q_weight, 
                wd_bar, 
                b0_bar, 
                modB_bar, 
                v_weight, 
                a_par_weight, 
                b_par_weight, 
                kx_bar, 
                kpar_bar, 
                field_weight_QL_out, 
                particle_weight, 
                energy_weight, 
                stress_par_weight, 
                stress_tor_weight, 
                exchange_weight = get_QL_weights(inputs, ky, v, eigenvalue, ave_p0inv, ave_b0inv, ave_bpinv, ave_wdh, ave_b0, ave_kx, ave_c_par_par, ave_kpar, ave_c_tor_par, ave_c_tor_per, ave_hp1)############### have to create this ###############

                #### probably outputs
                wd_bar_out[imax] = wd_bar
                b0_bar_out[imax] = b0_bar
                modB_bar_out[imax] = modB_bar
                v_QL_out[imax] = v_weight
                a_par_QL_out[imax] = a_par_weight
                b_par_QL_out[imax] = b_par_weight
                kx_bar_out[imax] = kx_bar
                kpar_bar_out[imax] = kpar_bar/(R_unit*q_unit*width_parameter)

                for i = 1:nbasis
                    for j = 1:3
                        field_weight_out[imax,j,i] = field_weight_QL_out[j,i]
                    end
                end

                for is = ns0:ns
                    for j = 1:3
                        particle_QL_out[imax,is,j] = particle_weight[is,j]
                        energy_QL_out[imax,is,j] = energy_weight[is,j]
                        stress_par_QL_out[imax,is,j] = stress_par_weight[is,j]
                        stress_tor_QL_out[imax,is,j] = stress_tor_weight[is,j]
                        exchange_QL_out[imax,is,j] = exchange_weight[is,j]
                    end
                    N_QL_out[imax,is] = N_weight[is]
                    T_QL_out[imax,is] = T_weight[is]
                    U_QL_out[imax,is] = U_weight[is]
                    Q_QL_out[imax,is] = Q_weight[is]
                    Ns_Ts_phase_out[imax,is] = Ns_Ts_phase[is]
                end
                ne_te_phase_out[imax] = Ne_Te_phase
          
                kyi = ky
                if(abs(v_QL_out[imax])<epsilon1)
                    v_bar_out[imax] = 0.0
                    phi2_bar = 0.0
                else
                    v_bar_out[imax] = get_intensity(inputs, ave, R_unit, kyi,gamma_out[imax]) ############### can use some cleaning
                    phi2_bar = v_bar_out[imax]/v_QL_out[imax]
                end

                phi_bar_out[imax] = phi2_bar
                a_par_bar_out[imax] = phi2_bar*a_par_QL_out[imax]
                b_par_bar_out[imax] = phi2_bar*b_par_QL_out[imax]

                for is = ns0:ns
                    N_bar_out[imax,is] = phi2_bar*N_QL_out[imax,is]
                    T_bar_out[imax,is] = phi2_bar*T_QL_out[imax,is]
                    U_bar_out[imax,is] = phi2_bar*U_QL_out[imax,is]
                    Q_bar_out[imax,is] = phi2_bar*Q_QL_out[imax,is]
                end
            end
        end

        # check for inward ballooing 
        ft_test = 0.0
        sum_modB_bar = 0.0
        sum_v_bar = 0.0
        for i = 1:nmodes_out
            sum_modB_bar = sum_modB_bar + v_bar_out[i]*modB_bar_out[i]
            sum_v_bar = sum_v_bar + v_bar_out[i]
        end
        if(sum_v_bar>epsilon1) ft_test = sum_modB_bar/sum_v_bar end
        modB_min = get_sat_params(:minB, inputs) ####### temporary value for now
        ft_test = ft_test/modB_min
    end

    return nmodes_out, gamma_out, freq_out,
    particle_QL_out, energy_QL_out, stress_tor_QL_out, stress_par_QL_out, exchange_QL_out

    # return  gamma_out, 
    #         freq_out, 
    #         v_QL_out, 
    #         a_par_QL_out,
    #         b_par_QL_out,
    #         phi_bar_out,
    #         a_par_bar_out,
    #         b_par_bar_out,
    #         v_bar_out,
    #         ne_te_phase_out,
    #         field_weight_out,
    #         particle_QL_out,
    #         energy_QL_out,
    #         stress_par_QL_out,
    #         stress_tor_QL_out,
    #         exchange_QL_out,
    #         N_QL_out,
    #         T_QL_out,
    #         U_QL_out,
    #         Q_QL_out,
    #         N_bar_out,
    #         T_bar_out,
    #         U_bar_out,
    #         Q_bar_out,
    #         Ns_Ts_phase_out

end

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------


function get_intensity(inputs::InputTJLF{T}, ave::Ave{T}, R_unit::T, kp::T, gp::T) where T<:Real

    nmodes_in = inputs.NMODES
    sat_rule_in = inputs.SAT_RULE
    etg_factor_in = inputs.ETG_FACTOR
    alpha_quench_in = inputs.ALPHA_QUENCH

    pols = (ave.p0[1,1]/abs(inputs.SPECIES[1].AS*inputs.SPECIES[1].ZS^2))^2
    ks = kp* √(inputs.SPECIES[1].TAUS*inputs.SPECIES[2].MASS)/abs(inputs.SPECIES[1].ZS)
    measure = √(inputs.SPECIES[1].TAUS*inputs.SPECIES[2].MASS)
    
    if(sat_rule_in==0)
        igeo = 1 ################### MILLER GEO for now
        if(igeo==0)
            if (nmodes_in<=2) #this fit is for nmodes_in=2
                cnorm = 30.40*pols
                exponent1 = 1.657
            else #this fit is for nmodes_in=4
                cnorm = 24.58*pols
                exponent1 = 1.761
            end
        
            if(ks>1.0) cnorm=cnorm/(ks)^etg_factor_in end
            c1 = 0.0

        elseif(igeo>=1)
            if(nmodes_in<=2) #this fit is for nmodes_in=2
                cnorm = 32.48*pols
                exponent1 = 1.547
                c1 = 0.534
            else #this fit is for nmodes_in=4
                cnorm = 30.03*pols
                exponent1 = 1.66
                c1 = 0.234
            end
        
            if(ks>1.0) cnorm=cnorm/(ks)^etg_factor_in end
        end

        wd0 = ks*√(inputs.SPECIES[1].TAUS*inputs.SPECIES[2].MASS)/R_unit  #renomalized for scale invariance
        gnet = gp/wd0
        intensity = cnorm*(wd0^2)*(gnet^exponent1 + c1*gnet)/(kp^4)
       
        if(alpha_quench_in==0.0 && abs(kx0_e)>0.0)
            intensity = intensity/(1.0+0.56*kx0_e^2)^2
            intensity = intensity/(1.0+(1.15*kx0_e)^4)^2
        end
        
        SAT_geo0_out = 1.0 ############# BRUH #################
        intensity = intensity*SAT_geo0_out*measure
    elseif(sat_rule_in>=1)
        intensity = 1.0
    end
       
    return intensity
end

#--------------------------------------------------------------

function get_gamma_net(inputs::InputTJLF{T}, vexb_shear_s::T, gp::T) where T<:Real
    alpha_quench_in = inputs.ALPHA_QUENCH
    kappa_loc = inputs.KAPPA_LOC ##### only true for MILLER
    alpha_exb = 0.3
    igeo = 1 ####### MILLER GEOMETRY for now

    if(igeo==1) alpha_exb=0.3*√(kappa_loc) end
    get_gamma_net =  max(gp - abs(alpha_exb*alpha_quench_in*vexb_shear_s),0.0)

    return get_gamma_net
end

#--------------------------------------------------------------

#compute the quasilinear weights for a single eigenmode
#with eigenvector v. All of the QL weights are normalized to phi_norm
function get_QL_weights(inputs::InputTJLF{T}, ave::Ave{T}, 
    ky::T, nbasis::T,
    use_bper_parameter::T,use_bpar_parameter::T) where T<:Real
    # , ky, v, eigenvalue,
    # ave_p0inv, ave_b0inv, ave_bpinv, ave_wdh, ave_b0, ave_kx, ave_c_par_par, ave_kpar, ave_c_tor_par, ave_c_tor_per, ave_hp1)

    epsilon1 = 1.e-12
    sat_rule_in = inputs.SAT_RULE
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    nroot=15 ### hardcoded
    iur = (ns-ns0+1)*nroot*nbasis

    vpar_model_in = inputs.VPAR_MODEL
    alpha_mach_in = inputs.ALPHA_MACH
    sign_It_in = inputs.SIGN_IT
    alpha_p_in = inputs.ALPHA_P
    freq_QL = eigenvalue

    n = Matrix{ComplexF64}(undef, ns,nbasis)
    u_par = Matrix{ComplexF64}(undef, ns,nbasis)
    p_par = Matrix{ComplexF64}(undef, ns,nbasis)
    p_tot = Matrix{ComplexF64}(undef, ns,nbasis)
    q_par = Matrix{ComplexF64}(undef, ns,nbasis)
    q_tot = Matrix{ComplexF64}(undef, ns,nbasis)


    vnorm = 0.0
    for is = ns0:ns
        j = (is-ns0)*nroot*nbasis
        for i = 1:nbasis          
            n[is,i] = v[j+i]
            u_par[is,i] = v[j+nbasis+i]
            p_par[is,i] = v[j+nbasis*2+i]
            p_tot[is,i] = v[j+nbasis*3+i]
            q_par[is,i] = v[j+nbasis*4+1] #### why is this just +1??? -DSUN
            q_tot[is,i] = v[j+nbasis*5+i]
            if(nroot>6)
                n[is,i] = n[is,i] -v[j+nbasis*6+i]+v[j+nbasis*12+i]
                u_par[is,i] = u_par[is,i] -v[j+nbasis*7+i]
                p_par[is,i] = p_par[is,i] -v[j+nbasis*8+i]+v[j+nbasis*13+i]
                p_tot[is,i] = p_tot[is,i] -v[j+nbasis*9+i]+v[j+nbasis*14+i]
                q_par[is,i] = q_par[is,i] -v[j+nbasis*10+i]
                q_tot[is,i] = q_tot[is,i] -v[j+nbasis*11+i]
            end
        end
    end


    # compute vnorm
    vnorm = 0.0
    if(ns<=2)
        vnorm = norm(v)
    else #weight vnorm by equililibrium densities
        j = 1
        for i = 1:iur
            if(i>j*nbasis*nroot) j=j+1 end
            as = inputs.SPECIES[j].AS
            zs = inputs.SPECIES[j].ZS
            vnorm = vnorm + real(v[i]*conj(v[i]))*abs(as*zs)
        end
        vnorm = vnorm/abs(inputs.SPECIES[1].AS*inputs.SPECIES[1].ZS)   #normalize to electron charge density
    end
    
    # compute the electromagnetic potentials
    betae_s = inputs.BETAE ##### not true for 'GENE' units
    betae_psi = 0.0
    if(use_bper_parameter) betae_psi = 0.5*betae_s/(ky*ky) end
    betae_sig = 0.0
    if(use_bpar_parameter) betae_sig = 0.5*betae_s end

    phi = zeros(ComplexF64, nbasis)
    psi = zeros(ComplexF64, nbasis)
    bsig = zeros(ComplexF64, nbasis)
    U0 = 0.0 ### defined in startup.f90
    for is = 1:ns
        vpar_s = inputs.ALPHA_MACH*inputs.SIGN_IT*inputs.SPECIES[is].VPAR
        taus = inputs.SPECIES[is].TAUS
        as = inputs.SPECIES[is].AS
        zs = inputs.SPECIES[is].ZS
        
        U0 = U0 + as*vpar_s*zs^2/taus
    end


    ### there must be some lienar algebra to clean this up
    for i = 1:nbasis
        for is = ns0:ns
            taus = inputs.SPECIES[is].TAUS
            mass = inputs.SPECIES[is].MASS
            as = inputs.SPECIES[is].AS
            zs = inputs.SPECIES[is].ZS
            vs = √(taus/mass)

            phi[i] += (ave.p0inv[i,:].*(as*zs))*transpose(n[is,:])

            if(use_bper_parameter)
                psi[i] +=( ave.b0inv[i,:].*(as*vs*zs*betae_psi))*transpose(u_par[is,:])
                
                if(vpar_model_in==0)
                    phi[i] += (ave.bpinv[i,:].*(U0*betae_psi*as*zs*vs)) *transpose(u_par[is,:])
                    psi[i] -= (ave.bpinv[i,:].*(U0*betae_psi*as*zs))    *transpose(n[is,:])
                end
            end
          
            if(use_bpar_parameter)
                bsig[i] = bsig[i] - betae_sig* as * taus * (1.5*p_tot[is,i]-0.5*p_par[is,i])
            end
        end
    end

    # add the adiabatic terms to the total moments
    for is = ns0:ns
        zs = inputs.SPECIES[is].ZS
        taus = inputs.SPECIES[is].TAUS
        n[is,:] .= n[is,:] .- (phi.*(zs/taus))
        p_par[is,:] .= p_par[is,:] .- (phi.*(zs/taus))
        p_tot[is,:] .= p_tot[is,:] .- (phi.*(zs/taus))
    end

    # compute phi_norm, psi_norm, sig_norm
    phi_norm = norm(phi)
    psi_norm = norm(psi)
    bsig_norm = norm(bsig)
    if(phi_norm<epsilon1) phi_norm = epsilon1 end

    #save the field weights
    field_weight_QL_out = Matrix{ComplexF64}(undef, 3,nbasis)
    field_weight_QL_out[1,:] .= phi .* (im/√(phi_norm))
    field_weight_QL_out[2,:] .= psi .* (im/√(phi_norm))
    field_weight_QL_out[3,:] .= bsig.* (im/√(phi_norm))

    #compute <phi|*|phi> averages
    phi_wd_phi = conj(phi) * ave_wdh * transpose(phi)
    phi_b0_phi = conj(phi) * ave_b0 * transpose(phi)
    phi_modB_phi =  conj(phi) * ave_c_par_par * transpose(phi)
    phi_kx_phi = conj(phi)* ave_kx * transpose(phi)
    phi_kpar_phi = conj(phi) * (im.*ave_kpar) * transpose(phi)

    wd_bar = real(phi_wd_phi)/phi_norm
    b0_bar = real(phi_b0_phi)/phi_norm
    modB_bar = abs(real(phi_modB_phi)/phi_norm)
    kx_bar = real(phi_kx_phi)/phi_norm
    kpar_bar = real(phi_kpar_phi)/phi_norm

    stress_correction = 1.0

    stress_par = Array{ComplexF64, 3}(undef, ns,nbasis,3)
    stress_per = Array{ComplexF64, 3}(undef, ns,nbasis,3)
    for is = ns0:ns
        vpar_shear_in = inputs.SPECIES[is].VPAR_SHEAR
        taus = inputs.SPECIES[is].TAUS
        mass = inputs.SPECIES[is].MASS
        vs = √(taus/mass)

        wp = ky*ave_hp1[is,1,1]*abs(alpha_p_in*vpar_shear_in)/vs
        if(sat_rule_in==0) stress_correction = (imag(freq_QL)+2.0*wp)/(imag(freq_QL)+wp) end

        stress_par[is,:,1] .= u_par[is,:]*stress_correction
        stress_par[is,:,2] .= p_par[is,:]*stress_correction
        stress_per[is,:,1] .= (im*ky) .* (1.5.*p_tot[is,:] .- 0.5.*p_par[is,:])*transpose(ave_kx) # (is,j) x (j,i)
        stress_per[is,:,2] .= (im*ky) .* (1.5.*q_tot[is,:] .- 0.5.*q_par[is,:])*transpose(ave_kx)
    end


    # compute the quasilinear weights for the fluxes
    particle_weight = zeros(Float64, ns, 3)
    energy_weight = zeros(Float64, ns, 3)
    stress_par_weight = zeros(Float64, ns, 3)
    stress_tor_weight = zeros(Float64, ns, 3)
    exchange_weight = zeros(Float64, ns, 3)

    for is = ns0:ns
        mass = inputs.SPECIES[is].MASS
        taus = inputs.SPECIES[is].TAUS
        vs = √(taus/mass)
        zs = inputs.SPECIES[is].ZS
        as = inputs.SPECIES[is].AS
        for i = 1:nbasis
            particle_weight[is,1] = particle_weight[is,1] + real(im*conj(phi[i])*n[is,i])
            energy_weight[is,1] = energy_weight[is,1] + real(im*conj(phi[i])*p_tot[is,i])
            
            for j = 1:nbasis
                stress_par_weight[is,1] = (stress_par_weight[is,1] 
                        + real(im*conj(phi[i])*ave_c_par_par[i,j]*stress_par[is,j,1]))         
                stress_tor_weight[is,1] = (stress_tor_weight[is,1] 
                        + real(im*conj(phi[i])
                        * (ave_c_tor_par[i,j]*stress_par[is,j,1]+ave_c_tor_per[i,j]*stress_per[is,j,1])))
            end

            exchange_weight[is,1] = (exchange_weight[is,1] 
                    + zs*real(im*freq_QL*conj(phi[i])*n[is,i]))
          
            if(use_bper_parameter)
                particle_weight[is,2] = (particle_weight[is,2]
                    - vs*real(im*conj(psi[i])*u_par[is,i]))
                
                energy_weight[is,2] = (energy_weight[is,2]
                    - vs*real(im*conj(psi[i])*q_tot[is,i]))
            
                exchange_weight[is,2] = (exchange_weight[is,2]
                    - zs*vs
                    *real(im*freq_QL*conj(psi[i])*u_par[is,i]))

                for j=1:nbasis
                    stress_par_weight[is,2] = (stress_par_weight[is,2]
                        - real(im*conj(psi[i])*ave_c_par_par[i,j]*stress_par[is,j,2]))        
                    stress_tor_weight[is,2] = (stress_tor_weight[is,2]
                        - real(im*conj(psi[i])*(ave_c_tor_par[i,j]*stress_par[is,j,2]
                        + ave_c_tor_per[i,j]*stress_per[is,j,2])))
                end
            end
            if(use_bpar_parameter)
                particle_weight[is,3] = (particle_weight[is,3]
                    + real(im*conj(bsig[i])*(1.5*p_tot[is,i]-0.5*p_par[is,i]))
                    * taus/zs)
                exchange_weight[is,3] = (exchange_weight[is,3]
                    + taus*real(conj(-im*freq_QL*bsig[i])
                    *(1.5*p_tot[is,i]-0.5*p_par[is,i])))
            end
        end

        for j = 1:3
            particle_weight[is,j] = as*ky*particle_weight[is,j]/phi_norm
            energy_weight[is,j] = as*taus*1.5*ky*energy_weight[is,j]/phi_norm
            stress_par_weight[is,j] = mass*as*vs*ky*stress_par_weight[is,j]/phi_norm
            stress_tor_weight[is,j] = sign_It_in*mass*as*vs*ky*stress_tor_weight[is,j]/phi_norm
            exchange_weight[is,j] = as*exchange_weight[is,j]/phi_norm
        end
    end


    #  add the vpar shifts to the total  moments
    if(vpar_model_in==0)
        for is = ns0:ns
            taus = inputs.SPECIES[is].TAUS
            mass = inputs.SPECIES[is].MASS
            vpar_s = alpha_mach_in*sign_It_in*inputs.SPECIES[is].VPAR
            zs = inputs.SPECIES[is].ZS
            vs = √(taus/mass)
            
            for j = 1:nbasis
                n[is,j] = n[is,j] + vpar_s*(zs/taus)*psi[j]
                u_par[is,j] = u_par[is,j] -(vpar_s/vs)*(zs/taus)*phi[j]
                p_par[is,j] = p_par[is,j] + vpar_s*(zs/taus)*psi[j]
                p_tot[is,j] = p_tot[is,j] + vpar_s*(zs/taus)*psi[j]
                q_par[is,j] = q_par[is,j] - 3.0*(vpar_s/vs)*(zs/taus)*phi[j]
                q_tot[is,j] = q_tot[is,j] -(5.0/3.0)*(vpar_s/vs)*(zs/taus)*phi[j]
            end
        end
    end


    #### outputs
    # compute the density and temperature amplitude weights
    N_weight = zeros(Float64, ns)
    T_weight = zeros(Float64, ns)
    U_weight = zeros(Float64, ns)
    Q_weight = zeros(Float64, ns)
    temp = Matrix{ComplexF64}(undef, ns,nbasis)
    
    for is = ns0:ns
        for i = 1:nbasis
            temp[is,i] = p_tot[is,i] - n[is,i]
            N_weight[is] = N_weight[is] + real(n[is,i]*conj(n[is,i]))
            T_weight[is] = T_weight[is] + real(temp[is,i]*conj(temp[is,i]))
            U_weight[is] = U_weight[is] + real(u_par[is,i]*conj(u_par[is,i]))
            Q_weight[is] = Q_weight[is] + real(q_tot[is,i]*conj(q_tot[is,i]))
        end
        N_weight[is] = N_weight[is]/phi_norm
        T_weight[is] = T_weight[is]/phi_norm
        U_weight[is] = U_weight[is]/phi_norm
        Q_weight[is] = Q_weight[is]/phi_norm
    end
    v_weight = vnorm/phi_norm
    a_par_weight = psi_norm/phi_norm
    b_par_weight = bsig_norm/phi_norm


    #compute electron density-temperature phase 
    Ne_Te_phase = 0.0
    Ne_Te_cos = 0.0
    Ne_Te_sin = 0.0
    for i = 1:nbasis
         Ne_Te_cos = Ne_Te_cos + real(conj(n[1,i])*temp[1,i])
         Ne_Te_sin = Ne_Te_sin + imag(conj(n[1,i])*temp[1,i])
    end
    Ne_Te_phase = atan(Ne_Te_sin,Ne_Te_cos)

    #compute species density-temperature phase
    Ns_Ts_phase = zeros(Float64, ns)
    for is = ns0:ns
        Ns_Ts_cos = 0.0
        Ns_Ts_sin = 0.0
        for i = 1:nbasis
            Ns_Ts_cos = Ns_Ts_cos + real(conj(n[is,i])*temp[is,i])
            Ns_Ts_sin = Ns_Ts_sin + imag(conj(n[is,i])*temp[is,i])
        end
        Ns_Ts_phase[is] = atan(Ns_Ts_sin,Ns_Ts_cos)
    end

    return Ns_Ts_phase, Ne_Te_phase,
    N_weight, T_weight, U_weight, Q_weight, 
    wd_bar, b0_bar, modB_bar, v_weight, a_par_weight, b_par_weight, kx_bar, kpar_bar, 
    field_weight_QL_out, particle_weight, energy_weight, stress_par_weight, stress_tor_weight, exchange_weight

end