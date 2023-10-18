import LinearAlgebra.LAPACK.gesv!
using Revise

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
            ky::T, 
            nbasis::Int,
            vexb_shear_s::T,
            gamma_reference_kx0::Union{Vector{T},Missing} = missing,
            freq_reference_kx0::Union{Vector{T},Missing} = missing) where T <: Real

    small = 1.0e-13
    epsilon1 = 1.0e-12
    nmodes_in = inputs.NMODES
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    # compute the rank of the eigenmatrix iur
    nroot=15
    iur = (ns-ns0+1)*nroot*nbasis

    alpha_quench_in = inputs.ALPHA_QUENCH
    new_eikonal_in = inputs.NEW_EIKONAL
    filter_in = inputs.FILTER
    
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
            outputGeo = xgrid_functions_geo(inputs, outputHermite, ky)
            R_unit = satParams.R_unit
            q_unit = satParams.q_unit
        end
    end  #new_eikonal_in

    new_matrix = true ####### hardcode for now
    if(new_matrix)
        # trace_path[7]=1
        ave, aveH, aveWH, aveKH, 
        aveG, aveWG, aveKG = get_matrix(inputs, outputGeo, outputHermite, ky, nbasis)
    end

    #  solver for linear eigenmodes of tglf equations
    eigenvalues, v = tjlf_eigensolver(inputs,outputGeo,satParams,ave,aveH,aveWH,aveKH,aveG,aveWG,aveKG,nbasis,ky)
    rr = real.(eigenvalues)
    ri = imag.(eigenvalues)

    # filter out numerical instabilities that sometimes occur with high mode frequency
    if filter_in>0.0
        max_freq = 2*abs(ave.wdh[1,1])/R_unit
        for is = ns0:ns
            rlts = inputs.RLTS[is]
            rlns = inputs.RLNS[is]
            as = inputs.AS[is]
            zs = inputs.ZS[is]
            test = abs(as*zs*(aveH.hp3p0[is,1,1]*rlns + 1.5*(aveH.hr13p0[is,1,1] - aveH.hp3p0[is,1,1])*rlts))
            max_freq = max(max_freq,test)
        end
        max_freq *= filter_in*abs(ky)

        rr .*= ifelse.((rr.>0.0)  .&  (abs.(ri).>max_freq), -1 , 1)
    end

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


    if(inputs.IBRANCH==0)
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
                if(rr[de[iroot]]>zgamax)
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


    elseif(inputs.IBRANCH==-1)
        # find the top nmodes most unstable modes
        ### put the unstable modes in ascending order by growthrate
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
        gamma_out .= gamma_reference_kx0 
        freq_out .= freq_reference_kx0
    end
    
    # get the fluxes for the most unstable modes
    # v = zeros(ComplexF64, iur)
    # zmat = Matrix{ComplexF64}(undef, iur, iur)
    if(inputs.IFLUX)
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
                # v .= small
                # for i = 1:iur
                #     for j = 1:iur
                #         # zmat[i,j] = (beta[jmax[imax]]*amat[i,j] - (small + alpha[jmax[imax]])*bmat[i,j])
                #         zmat[i,j] = (beta[jmax[imax]]*amat[i,j] - (alpha[jmax[imax]])*bmat[i,j])
                #         if i==j
                #             zmat[i,j] -= small
                #         end
                #     end
                # end
                # ### gesv!(A,B) solves Ax = B, A becomes LU factor, and B becomes solution
                # # gesv!(zmat,v)
                # v = zmat \ v

                # eigenvalue = im*alpha[jmax[imax]]/beta[jmax[imax]]
                
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
                exchange_weight = get_QL_weights(inputs, ave, aveH, ky, nbasis, eigenvalues[jmax[imax]], v[:,jmax[imax]])
                #### probably outputs
                wd_bar_out[imax] = wd_bar
                b0_bar_out[imax] = b0_bar
                modB_bar_out[imax] = modB_bar
                v_QL_out[imax] = v_weight
                a_par_QL_out[imax] = a_par_weight
                b_par_QL_out[imax] = b_par_weight
                kx_bar_out[imax] = kx_bar
                kpar_bar_out[imax] = kpar_bar/(R_unit*q_unit*inputs.WIDTH)

                field_weight_out[imax,:,:]  .= field_weight_QL_out
                particle_QL_out[imax,:,:]   .= particle_weight
                energy_QL_out[imax,:,:]     .= energy_weight
                stress_par_QL_out[imax,:,:] .= stress_par_weight
                stress_tor_QL_out[imax,:,:] .= stress_tor_weight
                exchange_QL_out[imax,:,:]   .= exchange_weight

                N_QL_out[imax,:] .= N_weight
                T_QL_out[imax,:] .= T_weight
                U_QL_out[imax,:] .= U_weight
                Q_QL_out[imax,:] .= Q_weight
                Ns_Ts_phase_out[imax,:] .= Ns_Ts_phase

                ne_te_phase_out[imax] = Ne_Te_phase
          
                if(abs(v_QL_out[imax])<epsilon1)
                    v_bar_out[imax] = 0.0
                    phi2_bar = 0.0
                else
                    kyi = ky
                    v_bar_out[imax] = get_intensity(inputs, ave, outputGeo.kx0_e, R_unit, kyi,gamma_out[imax]) ############### can use some cleaning
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
        modB_min = abs(satParams.minB_geo)
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


function get_intensity(inputs::InputTJLF{T}, ave::Ave{T}, kx0_e::T, R_unit::T, kp::T, gp::T) where T<:Real

    nmodes_in = inputs.NMODES
    sat_rule_in = inputs.SAT_RULE
    etg_factor_in = inputs.ETG_FACTOR
    alpha_quench_in = inputs.ALPHA_QUENCH

    pols = (ave.p0[1,1]/abs(inputs.AS[1]*inputs.ZS[1]^2))^2
    ks = kp* √(inputs.TAUS[1]*inputs.MASS[2])/abs(inputs.ZS[1])
    measure = √(inputs.TAUS[1]*inputs.MASS[2])
    
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

        wd0 = ks*√(inputs.TAUS[1]*inputs.MASS[2])/R_unit  #renomalized for scale invariance
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
function get_QL_weights(inputs::InputTJLF{T}, ave::Ave{T}, aveH::AveH{T},
    ky::T, nbasis::Int, eigenvalue::K, v::Vector{K}) where T<:Real where K<:Complex

    epsilon1 = 1.e-12
    sat_rule_in = inputs.SAT_RULE
    ns = inputs.NS
    ns0 = ifelse(inputs.ADIABATIC_ELEC, 2, 1)
    nroot=15 ### hardcoded
    iur = (ns-ns0+1)*nroot*nbasis


    taus = inputs.TAUS
    mass = inputs.MASS
    as = inputs.AS
    zs = inputs.ZS
    vs = .√(taus./mass)
    vpar = inputs.VPAR
    vpar_shear = inputs.VPAR_SHEAR
    
    vpar_model_in = inputs.VPAR_MODEL
    alpha_mach_in = inputs.ALPHA_MACH
    sign_It_in = inputs.SIGN_IT
    alpha_p_in = inputs.ALPHA_P
    freq_QL = im*eigenvalue

    n = Matrix{ComplexF64}(undef, ns,nbasis)
    u_par = Matrix{ComplexF64}(undef, ns,nbasis)
    p_par = Matrix{ComplexF64}(undef, ns,nbasis)
    p_tot = Matrix{ComplexF64}(undef, ns,nbasis)
    q_par = Matrix{ComplexF64}(undef, ns,nbasis)
    q_tot = Matrix{ComplexF64}(undef, ns,nbasis)

    for is = ns0:ns
        j = (is-ns0)*nroot*nbasis
        for i = 1:nbasis          
            n[is,i] = v[j+i]
            u_par[is,i] = v[j+nbasis+i]
            p_par[is,i] = v[j+nbasis*2+i]
            p_tot[is,i] = v[j+nbasis*3+i]
            q_par[is,i] = v[j+nbasis*4+i]
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
        vnorm = real((adjoint(v) * v))
    else #weight vnorm by equililibrium densities
        j = 1
        for i = 1:iur
            if(i>j*nbasis*nroot) j=j+1 end
            vnorm = vnorm + real((adjoint(v) * v))*abs(inputs.AS[j]*inputs.ZS[j])
        end
        vnorm = vnorm/abs(inputs.AS[1]*inputs.ZS[1])   #normalize to electron charge density
    end
    
    # compute the electromagnetic potentials
    betae_s = inputs.BETAE ##### not true for 'GENE' units
    betae_psi = 0.0
    if(inputs.USE_BPER) betae_psi = 0.5*betae_s/(ky*ky) end
    betae_sig = 0.0
    if(inputs.USE_BPAR) betae_sig = 0.5*betae_s end

    phi = zeros(ComplexF64, nbasis)
    psi = zeros(ComplexF64, nbasis)
    bsig = zeros(ComplexF64, nbasis)
    U0 = sum((alpha_mach_in*sign_It_in).*vpar.*zs.^2 .* as./taus) ### defined in startup.f90
    phi .= sum(ave.p0inv * transpose(Diagonal((as.*zs)[ns0:ns]) * n[ns0:ns,:]), dims=2)
    if(inputs.USE_BPER)
        psi .= (betae_psi .*     sum(ave.b0inv * transpose(Diagonal((as.*zs.*vs)[ns0:ns])*u_par[ns0:ns,:]), dims=2))
        if(vpar_model_in==0)
            phi .= phi .+   (U0*betae_psi)  .*sum(ave.bpinv * transpose(Diagonal((as.*zs.*vs)[ns0:ns])*u_par[ns0:ns,:]), dims=2)
            psi .= psi .-   (U0*betae_psi)  .*sum(ave.bpinv * transpose(Diagonal((as.*zs)[ns0:ns]) * n[ns0:ns,:]), dims=2)
        end
    end
    if(inputs.USE_BPAR)
        bsig .=  transpose((-betae_sig).*(as.*taus)[ns0:ns]) * (1.5.*p_tot[ns0:ns,:] .- 0.5.*p_par[ns0:ns,:])
    end

    # add the adiabatic terms to the total moments
    n[ns0:ns,:] .= n[ns0:ns,:] .- ((zs./taus)[ns0:ns] .* transpose(phi)) #### outer product to make matrix, idk why its this order tbh -DSUN
    p_par[ns0:ns,:] .= p_par[ns0:ns,:] .- ((zs./taus)[ns0:ns] .* transpose(phi))
    p_tot[ns0:ns,:] .= p_tot[ns0:ns,:] .- ((zs./taus)[ns0:ns] .* transpose(phi))

    # compute phi_norm, psi_norm, sig_norm
    phi_norm = real(adjoint(phi) * phi)
    psi_norm = real(adjoint(psi) * psi)
    bsig_norm = real(adjoint(bsig) * bsig)
    if(phi_norm<epsilon1) phi_norm = epsilon1 end

    #save the field weights
    field_weight_QL_out = Matrix{ComplexF64}(undef, 3,nbasis)
    field_weight_QL_out[1,:] .= phi .* (im/√(phi_norm))
    field_weight_QL_out[2,:] .= psi .* (im/√(phi_norm))
    field_weight_QL_out[3,:] .= bsig.* (im/√(phi_norm))

    #compute <phi|*|phi> averages
    phi_wd_phi = adjoint(phi) * ave.wdh * phi
    phi_b0_phi = adjoint(phi) * ave.b0 * phi
    phi_modB_phi =  adjoint(phi) * ave.c_par_par * phi
    phi_kx_phi = adjoint(phi)* ave.kx * phi
    phi_kpar_phi = adjoint(phi) * (im.*ave.kpar) * phi

    wd_bar = real(phi_wd_phi)/phi_norm
    b0_bar = real(phi_b0_phi)/phi_norm
    modB_bar = abs(real(phi_modB_phi)/phi_norm)
    kx_bar = real(phi_kx_phi)/phi_norm
    kpar_bar = real(phi_kpar_phi)/phi_norm

    
    # fill the stress moments
    stress_par = Array{ComplexF64, 3}(undef, ns,nbasis,3)
    stress_per = Array{ComplexF64, 3}(undef, ns,nbasis,3)

    stress_correction = 1.0
    if(sat_rule_in==0) 
        wp = (ky*abs(alpha_p_in)) .* aveH.hp1[ns0:ns,1,1].*vpar_shear./vs
        stress_correction = (imag(freq_QL).+2.0.*wp)./(imag(freq_QL).+wp)
    end

    stress_par[ns0:ns,:,1] .= u_par[ns0:ns,:].*stress_correction
    stress_par[ns0:ns,:,2] .= p_par[ns0:ns,:].*stress_correction
    stress_per[ns0:ns,:,1] .= (im*ky) .* (1.5 .*p_tot[ns0:ns,:] .- 0.5 .*p_par[ns0:ns,:]) * (ave.kx)' # (is,j) x (j,i)
    stress_per[ns0:ns,:,2] .= (im*ky) .* (1.5 .*q_tot[ns0:ns,:] .- 0.5 .*q_par[ns0:ns,:]) * (ave.kx)'


    # compute the quasilinear weights for the fluxes
    particle_weight = zeros(Float64, ns, 3)
    energy_weight = zeros(Float64, ns, 3)
    stress_par_weight = zeros(Float64, ns, 3)
    stress_tor_weight = zeros(Float64, ns, 3)
    exchange_weight = zeros(Float64, ns, 3)

    ### real() with an im in it is funky
    #### CHECK THIS PLEASE
    particle_weight[ns0:ns,1]   .= real.(im.*n[ns0:ns,:]*conj(phi))
    energy_weight[ns0:ns,1]     .= real.(im.*p_tot[ns0:ns,:]*conj(phi))
    stress_par_weight[ns0:ns,1] .= real.(im.* (stress_par[ns0:ns,:,1]*ave.c_par_par')*conj(phi))
    stress_tor_weight[ns0:ns,1] .= real.(im.* (stress_par[ns0:ns,:,1]*ave.c_tor_par' 
                                            .+ stress_per[ns0:ns,:,1]*ave.c_tor_per')*conj(phi))
    exchange_weight[ns0:ns,1]   .= Diagonal(zs[ns0:ns]) * real.((im*freq_QL).* n[ns0:ns,:]*conj(phi))
    

    if(inputs.USE_BPER)
        particle_weight[ns0:ns,2]   .= - Diagonal(vs[ns0:ns])      * real.(im.* u_par[ns0:ns,:]*conj(psi))
        energy_weight[ns0:ns,2]     .= - (Diagonal(vs[ns0:ns])     * real.(im.* q_tot[ns0:ns,:]*conj(psi)) 
                                      .+ Diagonal((zs.*vs)[ns0:ns])* real.((im*freq_QL).* u_par[ns0:ns,:]*conj(psi)))
        stress_par_weight[ns0:ns,2] .= - real.(im.* (stress_par[ns0:ns,:,2]*ave.c_par_par')*conj(psi))      
        stress_tor_weight[is,2]     .= - real.(im.* (stress_par[ns0:ns,:,2]*ave.c_tor_par'
                                                  .+ stress_per[ns0:ns,:,2]*ave.c_tor_per')*conj(psi))
    end
    if(inputs.USE_BPAR)
        particle_weight[ns0:ns,3]   .= real.(im.* Diagonal((taus./zs)[ns0:ns]) *(1.5*p_tot[ns0:ns,:] .- 0.5*p_par[ns0:ns,:])*conj(bsig))
        exchange_weight[ns0:ns,3]   .= Diagonal(taus[ns0:ns]) * real.(          (1.5*p_tot[ns0:ns,:] .- 0.5*p_par[ns0:ns,:])*conj((-im*freq_QL).*bsig))
    end

    particle_weight[ns0:ns,:]   .= ky.*       (Diagonal(as[ns0:ns])*particle_weight[ns0:ns,:])./phi_norm
    energy_weight[ns0:ns,:]     .= (1.5*ky).* (Diagonal((as.*taus)[ns0:ns])*energy_weight[ns0:ns,:])./phi_norm
    stress_par_weight[ns0:ns,:] .= ky.*       (Diagonal((as.*mass.*vs)[ns0:ns])*stress_par_weight[ns0:ns,:])./phi_norm
    stress_tor_weight[ns0:ns,:] .= ky.*sign_It_in.*(Diagonal((mass.*as.*vs)[ns0:ns])*stress_tor_weight[ns0:ns,:])./phi_norm
    exchange_weight[ns0:ns,:]   .= (Diagonal(as[ns0:ns])*exchange_weight[ns0:ns,:])./phi_norm


    #  add the vpar shifts to the total  moments
    if(vpar_model_in==0)
        vpar_s = (alpha_mach_in*sign_It_in).*vpar
        n[ns0:ns,:] = (vpar_s.*zs./taus)[ns0:ns] .* transpose(psi) ### outer product is slightly weird
        u_par[ns0:ns,:] = (-(vpar_s./vs).*(zs./taus))[ns0:ns] .* transpose(phi)
        p_par[ns0:ns,:] = (vpar_s.*(zs./taus))[ns0:ns] .* transpose(psi)
        p_tot[ns0:ns,:] = (vpar_s.*(zs./taus))[ns0:ns] .* transpose(psi)
        q_par[ns0:ns,:] = (-3 .*(vpar_s./vs).*(zs./taus))[ns0:ns] .* transpose(phi)
        q_tot[ns0:ns,:] = (-(5/3).*(vpar_s./vs).*(zs./taus))[ns0:ns] .* transpose(phi)
    end



    #### outputs
    # compute the density and temperature amplitude weights
    N_weight = zeros(Float64, ns)
    T_weight = zeros(Float64, ns)
    U_weight = zeros(Float64, ns)
    Q_weight = zeros(Float64, ns)
    temp = Matrix{ComplexF64}(undef, ns, nbasis)
    
    temp .= p_tot .- n
    N_weight .= sum(real(n     .* conj(n)), dims=2)./phi_norm
    T_weight .= sum(real(temp  .* conj(temp)), dims=2)./phi_norm
    U_weight .= sum(real(u_par .* conj(u_par)), dims=2)./phi_norm
    Q_weight .= sum(real(q_tot .* conj(q_tot)), dims=2)./phi_norm
    v_weight = vnorm/phi_norm
    a_par_weight = psi_norm/phi_norm
    b_par_weight = bsig_norm/phi_norm


    #compute electron density-temperature phase 
    Ne_Te_cos = real(adjoint(n[1,:])*temp[1,:])
    Ne_Te_sin = imag(adjoint(n[1,:])*temp[1,:])
    Ne_Te_phase = atan(Ne_Te_sin,Ne_Te_cos)

    #compute species density-temperature phase
    Ns_Ts_phase = zeros(Float64, ns)
    Ns_Ts_cos = zeros(Float64, ns)
    Ns_Ts_sin = zeros(Float64, ns)
    Ns_Ts_cos .= sum(real(conj(n[ns0:ns,:]).*temp[ns0:ns,:]), dims=2)
    Ns_Ts_sin .= sum(imag(conj(n[ns0:ns,:]).*temp[ns0:ns,:]), dims=2)
    Ns_Ts_phase .= atan.(Ns_Ts_sin,Ns_Ts_cos)
    atan

    return Ns_Ts_phase, Ne_Te_phase,
    N_weight, T_weight, U_weight, Q_weight, 
    wd_bar, b0_bar, modB_bar, v_weight, a_par_weight, b_par_weight, kx_bar, kpar_bar, 
    field_weight_QL_out, particle_weight, energy_weight, stress_par_weight, stress_tor_weight, exchange_weight

end