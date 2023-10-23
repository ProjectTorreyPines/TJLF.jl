# export xgrid_functions_geo, get_sat_params, mercier_luc, miller_geo

include("tjlf_modules.jl")
include("tjlf_multiscale_spectrum.jl")

function xgrid_functions_geo(inputs::InputTJLF, satParams::SaturationParameters{T}, ky::Vector{T}, gammas::Matrix{T},
    small::T=0.00000001) where T<:Real
    sign_IT = inputs.SIGN_IT
    vexb_shear = inputs.VEXB_SHEAR
    sat_rule_in = inputs.SAT_RULE
    alpha_quench_in = inputs.ALPHA_QUENCH
    alpha_e_in = inputs.ALPHA_E
    
    mass_2 = inputs.MASS[2]
    taus_2 = inputs.TAUS[2]
    zs_2 = inputs.ZS[2]
    units_in = inputs.UNITS
    nx = 2*inputs.NXGRID - 1
    vs_2 = √(taus_2 / mass_2)

    grad_r0 = satParams.grad_r0
    B_geo0 = satParams.B_geo[1]

    # kx0 = kx0_loc/ky # note that kx0 is kx/ky

    # generalized quench rule kx0 shift
    gamma_reference_kx0 = gammas[:, 1]
    if(alpha_quench_in==0.0 && gamma_reference_kx0[1]!=0.0)
        vexb_shear_s = vexb_shear * sign_IT
        vexb_shear_kx0 = alpha_e_in*vexb_shear_s

        kyi = ky.*(vs_2*mass_2/abs(zs_2))
        wE = zeros(Float64, length(ky))

        if(units_in=="GYRO")
            kx0_factor = abs(grad_r0^2/B_geo0)
            kx0_factor = 1.0 + 0.40*(kx0_factor-1.0)^2
            wE .= (kx0_factor*vexb_shear_kx0) .* (ifelse.(kyi.<0.3,kyi./0.3,1.0))./gamma_reference_kx0
        else
            kx0_factor = 1.0
        end

        kx0_e = -(0.36*vexb_shear_kx0./gamma_reference_kx0
                .+ (0.38.*wE.*tanh.((0.69.*wE).^6)))
        
        a0 = 1.3
        if(sat_rule_in==1)
            a0 = 1.45
            kx0_e = -(0.53*vexb_shear_kx0./gamma_reference_kx0 .+
                    (0.25.*wE.*tanh.((0.69.*wE).^6)))
        elseif(sat_rule_in==2 || sat_rule_in==3)
            a0=1.6
            vzf_out, kymax_out, _ = get_zonal_mixing(inputs, satParams, ky, gamma_reference_kx0)
            if(abs(kymax_out*vzf_out*vexb_shear_kx0) > small)
                kx0_e = -(0.32*vexb_shear_kx0).*((ky./kymax_out).^0.3)./(ky.*vzf_out)
            else
                kx0_e = zeros(Float64,length(ky))
            end
        end

        kx0_e = ifelse.(abs.(kx0_e) .> a0 , a0.*kx0_e./abs.(kx0_e) , kx0_e)
        kx0_e = ifelse.(isnan.(kx0_e) , 0 , kx0_e)

        #### not in the Python
        # if(units_in=="GYRO")
        #     kx0 = sign_Bt_in*kx0_e # cancel the sign_Bt_in factor in kxx below
        # else
        #     if(sat_rule_in.eq.1)kx0 = sign_Bt_in*kx0_e/(2.1)end  # goes with xnu_model=2
        #     if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)kx0 = sign_Bt_in*kx0_e*0.7/grad_r0_out^2 end     # goes with xnu_model=3, the factor 0.7/grad_r0_out^2 is needed for stress_tor
        #     # note kx0 = alpha_e*gamma_ExB_HB/gamma Hahm - Burrell form of gamma_ExB
        #     # The 2.1 effectively increases ay0 & ax0 and reduces toroidal stress to agree with CGYRO
        # end
    end

    return kx0_e
end


function xgrid_functions_geo(inputs::InputTJLF{T}, satParams::SaturationParameters{T}, outHermite::OutputHermite{T},
    ky::T, kx0_e::T=0.0,
    mts::T=5.0, ms::Int=128, small::T=0.00000001) where T<:Real

    width_in = inputs.WIDTH
    sat_rule_in = inputs.SAT_RULE

    ### different for different geometries!!!
    rmaj_s = inputs.RMAJ_LOC
    rmin_s = inputs.RMIN_LOC
    q_s = inputs.Q_LOC

    sign_Bt_in = inputs.SIGN_BT
    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end
    nx = 2*inputs.NXGRID - 1
    
    #*************************************************************
    # begin calculation of wdx and b0x
    #*************************************************************    
    y = satParams.y
    Ly = y[ms+1]
    R = satParams.R
    Bp = satParams.Bp
    R_unit = satParams.R_unit
    q_unit = satParams.q_unit
    b_geo = satParams.B_geo
    qrat_geo = satParams.qrat_geo
    S_prime = satParams.S_prime
    kx_factor = satParams.kx_factor
    sintheta_geo = satParams.sintheta_geo
    costheta_geo = satParams.costheta_geo
    costheta_p_geo = satParams.costheta_p_geo
    grad_r0_out = satParams.grad_r0

    f = satParams.Bt0 * inputs.RMAJ_LOC # Bt0_out = f/rmaj_input defined

    kx0 = inputs.KX0_LOC/ky
    if(inputs.UNITS=="GYRO")
        kx0 = sign_Bt_in*kx0_e 
    else
       if(sat_rule_in==1) kx0 = sign_Bt_in*kx0_e/(2.1) end
       if(sat_rule_in==2 || sat_rule_in==3) kx0 = sign_Bt_in*kx0_e*0.7/grad_r0_out^2 end
    end

    x = outHermite.x
    kxx = Vector{Float64}(undef,nx)
    wdx = Vector{Float64}(undef,nx)
    wdpx = Vector{Float64}(undef,nx)
    b0x = Vector{Float64}(undef,nx)
    b2x = Vector{Float64}(undef,nx)
    cx_tor_par = Vector{Float64}(undef,nx)
    cx_tor_per = Vector{Float64}(undef,nx)
    cx_par_par = Vector{Float64}(undef,nx)
    for i = 1:nx
        thx = width_in * x[i]
        sign_theta = ifelse(thx>=0, 1.0, -1.0)

        loops = Int(floor(abs(thx/(2π))))
        y_x = Ly*(abs(thx) - loops*2π)/ (2π)
        if(thx<0.0)
            y_x = Ly - y_x
            loops = loops + 1
        end
        lastIndex = 2
        for m = 2:ms+1
            lastIndex = m
            if(y[m]>=y_x) break end
        end
        m1 = lastIndex - 1
        m2 = lastIndex

        dkxky1 = sign_theta*loops*S_prime[ms+1]
        y1 = y[m1]
        dkxky2 = sign_theta*loops*S_prime[ms+1]
        y2 = y[m2]

        ### kxx
        kxx1 = (kx_factor[m1]*(S_prime[m1]+dkxky1)  -   kx0*b_geo[m1]/qrat_geo[m1]^2)    *qrat_geo[m1]/b_geo[m1]
        kxx2 = (kx_factor[m2]*(S_prime[m2]+dkxky2)  -   kx0*b_geo[m2]/qrat_geo[m2]^2)    *qrat_geo[m2]/b_geo[m2]
        kxx[i] = kxx1 + (kxx2 - kxx1)*(y_x - y1)/(y2 - y1)
        kxx[i] = sign_Bt_in*kxx[i]
        
        ### wdx
        wd1 = (qrat_geo[m1]/b_geo[m1])*     (costheta_geo[m1] + (kx_factor[m1]*(S_prime[m1]+dkxky1) - kx0*b_geo[m1]/qrat_geo[m1]^2)*sintheta_geo[m1])
        wd2 = (qrat_geo[m2]/b_geo[m2])*     (costheta_geo[m2] + (kx_factor[m2]*(S_prime[m2]+dkxky2) - kx0*b_geo[m2]/qrat_geo[m2]^2)*sintheta_geo[m2])
        wdx[i] = wd1 + (wd2 - wd1)*(y_x - y1)/(y2 - y1)
        wdx[i] = (R_unit/rmaj_s)*wdx[i]

        ### wdpx
        wdp1 = (qrat_geo[m1]/b_geo[m1])*costheta_p_geo[m1]
        wdp2 = (qrat_geo[m2]/b_geo[m2])*costheta_p_geo[m2]
        wdpx[i] = wdp1 + (wdp2 - wdp1)*(y_x - y1)/(y2 - y1)
        wdpx[i] = (R_unit/rmaj_s)*wdpx[i]
        
        ### b2x
        b1 = b_geo[m1]^2
        b2 = b_geo[m2]^2
        b2x[i] =  b1 +(b2-b1)*(y_x-y1)/(y2-y1)

        ### b0x
        b1 = (1.0 + (kx_factor[m1]*(S_prime[m1]+dkxky1) - kx0*b_geo[m1]/qrat_geo[m1]^2)^2)*qrat_geo[m1]^2
        b2 = (1.0 + (kx_factor[m2]*(S_prime[m2]+dkxky2) - kx0*b_geo[m2]/qrat_geo[m2]^2)^2)*qrat_geo[m2]^2
        b0x[i] =  b1 +(b2 - b1)*(y_x - y1)/(y2 - y1)

        if(b0x[i]<0.0)
            error("interpolation error b0x < 0")
        end

        ### viscous stress projection coefficients
        cxtorper1 = -R[m1]*Bp[m1]/b_geo[m1]
        cxtorper2 = -R[m2]*Bp[m2]/b_geo[m2]
        cx_tor_par[i] = f/b_geo[m1] + (f/b_geo[m2]-f/b_geo[m1])*(y_x-y1)/(y2-y1)
        cx_tor_par[i] = sign_Bt_in*cx_tor_par[i]
        cx_tor_per[i] = cxtorper1 + (cxtorper2-cxtorper1)*(y_x-y1)/(y2-y1)
        cx_par_par[i] = b_geo[m1] + (b_geo[m2]-b_geo[m1])*(y_x-y1)/(y2-y1)
    end

    #*************************************************************
    # calculate fts (fraction of trapped electrons)
    #*************************************************************

    nb_grid = 25
    By = Vector{Float64}(undef, nb_grid+1)
    pm = Matrix{Int}(undef, 2, nb_grid+1)

    Bmax, m_max = findmax(b_geo)
    Bmin, m_min = findmin(b_geo)

    By[1] = Bmin
    db = (Bmax - Bmin)/(nb_grid)
    for i = 2:nb_grid+1
        By[i] = By[i-1] + db
    end

    # find pairs of m's at the same B starting at Bmin and taking the farthest pair
    pm[1,1] = m_min
    pm[2,1] = m_min
    qm = 1
    for i = 2:nb_grid
        ### get the number of values between m_max and m_min
        j_max = m_max - m_min
        if(j_max<0) j_max = j_max + ms end
        ### increment through them
        for m = 1:j_max
            # find the farthest gridpoint where b_geo=By
            j = m_min + m
            ### remember 1 = ms+1
            if(j>ms+1) j = j - ms end

            test1 = b_geo[j-1] - By[i]
            test2 = b_geo[j] - By[i]

            if(test1*test2<=0)
                if(abs(test1)<abs(test2))
                    qm = j-1
                else
                    qm = j
                end
            end
        end
        pm[1,i] = qm
        
        #### reverse direction
        j_max = m_min - m_max
        if(j_max<0) j_max = j_max + ms end
        qm = 1
        for m = 1:j_max
            # find the farthest gridpoint where b_geo=By
            j = m_min - m
            if(j<1) j = j + ms end

            test1 = b_geo[j+1] - By[i]
            test2 = b_geo[j] - By[i]

            if(test1*test2<=0.0)
                if(abs(test1)<abs(test2))
                    qm = j+1
                else
                    qm = j
                end
            end
        end
        pm[2,i] = qm
    end
    pm[1,nb_grid+1] = m_max
    pm[2,nb_grid+1] = m_max
    for i = 1:nb_grid+1
        if(pm[1,i]>ms+1 || pm[2,i]>ms+1)
            error("error in get_ft_geo: pm out of bounds")
        end
    end

    delta_y = Vector{Float64}(undef, nb_grid+1)
    delta_y[1] = 0
    for i = 2:nb_grid+1
        if(y[pm[1,i]]>y[pm[2,i]])
            delta_y[i] =      y[pm[1,i]] - y[pm[2,i]]
        else
            delta_y[i] = Ly + y[pm[1,i]] - y[pm[2,i]]
        end
    end
    #*************************************************************
    # compute trapped fraction
    #*************************************************************
    ### This can be made neater for sure -DSUN
    kpar = 2π/(Ly*√(2)*width_in)
    bounce_y = min(Ly, π*inputs.THETA_TRAPPED/kpar)

    B_bounce = Bmax
    if(bounce_y<Ly)
        index = 1
        for i = 2:nb_grid+1
            index = i
            if(delta_y[i]>bounce_y) break end
        end
        B_bounce = By[index-1] + (By[index] - By[index-1])*(bounce_y-delta_y[index-1])/(delta_y[index]-delta_y[index-1])
    end
    ft = √(1.0 - Bmin/B_bounce)
    modB_min = abs(Bmin)
    modB_test = 0.5*(Bmax + Bmin)/Bmin
    fts = Vector{Float64}(undef, ns)
    for is = ns0:ns
        fts[is] = max(ft,0.01)
    end

    xnu_model_in = inputs.XNU_MODEL
    wdia_trapped_in = inputs.WDIA_TRAPPED
    if(xnu_model_in==3 && wdia_trapped_in>0.0) 
        theta_trapped_in = inputs.THETA_TRAPPED
        for is = ns0:ns
            taus = inputs.TAUS[is]
            mass = inputs.MASS[is]
            rlns = inputs.RLNS[is]
            vs = √(taus/mass)

            wdia = abs(ky*rlns)/vs
            kpar = 2π/(Ly*√(2)*width_in)
            ft0 = √(1.0 - Bmin/Bmax)
            cdt = 3*wdia_trapped_in*(1-ft0^2)
            kpar = kpar/max(theta_trapped_in,0.0001) + wdia*cdt
            bounce_y = min(Ly,π/kpar)
            B_bounce = Bmax
            if(bounce_y<Ly)
                index = 1
                for i = 2:nb_grid+1
                    index = i
                    if(delta_y[i]>bounce_y) break end
                end
                B_bounce = By[index-1] + (By[index] - By[index-1])*(bounce_y-delta_y[index-1])/(delta_y[index]-delta_y[index-1])
            end
            fts[is] = max(√(1 - Bmin/B_bounce),0.01)
        end
    end

    return OutputGeometry{Float64}(0, fts, kxx,wdx,wdpx,b0x,b2x,cx_tor_par,cx_tor_per,cx_par_par)

end








#### LINES 220-326, 478 in tglf_geometry.f90
function get_sat_params(inputs::InputTJLF, mts::T=5.0, ms::Int=128, small::T=0.00000001) where T<:Real
    #******************************************************************************#************************
    # PURPOSE: compute the geometric coefficients on the x-grid
    #******************************************************************************#************************
    
    ### different for different geometries!!!
    rmaj_s = inputs.RMAJ_LOC
    rmin_s = inputs.RMIN_LOC
    q_s = inputs.Q_LOC

    
    sat_rule_in = inputs.SAT_RULE
    units_in = inputs.UNITS
    sign_Bt_in = inputs.SIGN_BT
    ns = inputs.NS
    ns0 = 1
    if(inputs.ADIABATIC_ELEC) ns0 = 2 end
    nx = 2*inputs.NXGRID - 1

    s_p, Bp, b_geo, pk_geo, qrat_geo, sintheta_geo, costheta_geo, costheta_p_geo, Bt0_out, B_unit_out, ds, t_s, f, R, B, S_prime, kx_factor = mercier_luc(inputs)
    

    #*************************************************************
    # find length along magnetic field y
    #*************************************************************
    y = zeros(Float64,ms+1)
    y[1]=0.0
    for m in 2:ms+1
        y[m] = y[m-1]+s_p[m]*ds*4.0/(pk_geo[m]+pk_geo[m-1])
    end
    # set the global units
    Ly = y[ms+1]
    R_unit = rmaj_s*b_geo[1]/(qrat_geo[1]*costheta_geo[1])
    q_unit = Ly/(2π*R_unit)

    # # # save f for output
    # # RBt_ave_out = f/B_unit


    #*************************************************************
    #  compute flux surface averages
    #*************************************************************
    norm_ave = 0.0
    SAT_geo1_out = 0.0
    SAT_geo2_out = 0.0
    for m in 2:ms+1
        dlp = s_p[m]*ds*(0.5/Bp[m]+0.5/Bp[m-1])
        norm_ave += dlp
        # B2x1 = b_geo[i-1]^2
        # B2x2 = b_geo[i]^2
        # B2_ave_out = B2_ave_out + dlp*(B2x1+B2x2)/2.0
        # R2x1 = R(i-1)**2
        # R2x2 = R[i]**2
        # R2_ave_out = R2_ave_out + dlp*(R2x1+R2x2)/2.0
        # B_ave_out = B_ave_out + dlp*(b_geo(i-1)+b_geo[i])/2.0
        # Bt_ave_out = Bt_ave_out + dlp*(f/b_geo(i-1)+f/b_geo[i])/(2.0*Rmaj_s)
        # Grad_r_ave_out = Grad_r_ave_out + dlp*0.5*((R(i-1)*Bp(i-1))**2+(R[i]*Bp[i])**2)*(q_s/rmin_s)**2
        # kykx_geo_ave = kykx_geo_ave + dlp*0.5*(B_geo(i-1)**2/qrat_geo(i-1)**4+B_geo[i]**2/qrat_geo[i]**4)
        SAT_geo1_out += dlp*((b_geo[1]/b_geo[m-1])^4 +(b_geo[1]/b_geo[m])^4)/2.0
        SAT_geo2_out += dlp*((qrat_geo[1]/qrat_geo[m-1])^4 +(qrat_geo[1]/qrat_geo[m])^4)/2.0
    end
    # R2_ave_out = R2_ave_out/norm_ave
    # B2_ave_out = B2_ave_out/norm_ave
    # B2_ave_out = B2_ave_out/B_unit**2
    # B_ave_out = B_ave_out/norm_ave
    # Bt_ave_out = Bt_ave_out/norm_ave
    # Grad_r_ave_out = Grad_r_ave_out/norm_ave
    # kykx_geo_ave = kykx_geo_ave/norm_ave
    SAT_geo1_out = SAT_geo1_out/norm_ave
    SAT_geo2_out = SAT_geo2_out/norm_ave
    #
    # poloidal magnetic field on outboard midplane
    #
    # Bp0_out = Bp[1]/B_unit
    if(units_in=="GYRO")
        SAT_geo0_out = 1.0
        SAT_geo1_out = 1.0
        SAT_geo2_out = 1.0
    else
        SAT_geo0_out = 0.946/qrat_geo[1]          # normed to GASTD with CGYRO
        if(sat_rule_in==2 || sat_rule_in==3) SAT_geo0_out = 1.0 end
    end
    ### line 276-277
    grad_r0_out = b_geo[1]/qrat_geo[1]
    # Additional outputs for SAT2 G1(theta), Gq(theta)
    theta_out = t_s  # theta grid over which everything is calculated.

    return SaturationParameters{Float64}(SAT_geo0_out,SAT_geo1_out,SAT_geo2_out,
                                    y, R_unit, B_unit_out[end], q_unit,
                                    R, Bp, B,
                                    Bt0_out, grad_r0_out,
                                    S_prime,kx_factor,
                                    b_geo, qrat_geo, 
                                    sintheta_geo, costheta_geo, costheta_p_geo,
                                    theta_out)

end









function mercier_luc(inputs::InputTJLF,
    mts::Float64=5.0, ms::Int=128, small::Float64=0.00000001)
    #-------------------------------------------
    # the following must be defined from a previous call to one of the
    # geometry routines miller_geo, fourier_geo,ELITE_geo and stored in tglf_sgrid:
    # ms # the number of points in the s-grid (flux surface contour)
    # ds # the arc length differential on a flux surface
    # R(ms) # the major radius on the s-grid
    # Z(ms) # the vertical coordinate on the s-grid
    # Bp(ms) # the poloidal magnetic field on the s-grid normalized to B_unit
    # q_s = local flux surface safety factor
    # q_prime_s = dq/dpsi
    # p_prime_s = dp/dpsi
    #-----------------------
    #
    # compute the first and second derivatives of R,Z on the s-grid
    # and the local radius of curvature.
    # Note that for the Mercier-Luc coordinate dR/ds = cos(u), dZ/ds = -sin(u)
    # so (dR/ds)**2+(dZ/ds)**2 = 1, error_check compute the error in this relation
    # to make sure that the input flux surface coordinates R(s), Z(s) are ok.
    #


    ### geometry dependent!!!!!!, maybe save from output of miller_geo()?
    q_s = inputs.Q_LOC
    rmaj_input = inputs.RMAJ_LOC
    rmin_s = inputs.RMIN_LOC
    rmaj_s = inputs.RMAJ_LOC
    b_geo = zeros(Float64, ms+1)

    ### technically don't have to initialize here, but maybe better looking?
    ### if remove, switch .= to = for the below functions
    Bt = zeros(Float64, ms+1)
    B = zeros(Float64, ms+1)
    pk_geo = zeros(Float64, ms+1)
    qrat_geo = zeros(Float64, ms+1)
    
    

    R, Bp, Z, q_prime_s, p_prime_s, B_unit_out, ds, t_s = miller_geo(inputs)
    
    
    psi_x = R.*Bp
    delta_s = 12.0*ds
    ds2 = 12.0*ds^2


    # note that the point 1 and ms+1 are the same so m+1->1 and m-1->ms-1 at m=0
    s_p = zeros(Float64,ms+1)
    r_curv = zeros(Float64,ms+1)
    sin_u = zeros(Float64,ms+1)
    for m in 1:ms + 1
        m1 = ((ms + m - 2) % (ms+1)) + 1
        m2 = ((ms + m - 1) % (ms+1)) + 1
        m3 = (m % (ms+1)) + 1
        m4 = ((m + 1) % (ms+1)) + 1
        ### had to do this weird short-circuiting bc of weird indexing
        m1 < ms || (m1 = ((ms + m - 3) % (ms+1)) + 1)
        m2 < ms+1 || (m2 = ((ms + m - 2) % (ms+1)) + 1)
        m3 > 1 || (m3 = ((m + 1) % (ms+1)) + 1)
        m4 > 2 || (m4 = ((m + 2) % (ms+1)) + 1)

        R_s = (R[m1] - 8.0*R[m2] + 8.0*R[m3] - R[m4])/delta_s
        Z_s = (Z[m1] - 8.0*Z[m2] + 8.0*Z[m3] - Z[m4])/delta_s
        s_p[m] = √(R_s^2 + Z_s^2)

        R_ss = (-R[m1] + 16.0*R[m2] - 30.0*R[m] + 16.0*R[m3] - R[m4])/ds2
        Z_ss = (-Z[m1] + 16.0*Z[m2] - 30.0*Z[m] + 16.0*Z[m3] - Z[m4])/ds2
        r_curv[m] = (s_p[m]^3)/(R_s*Z_ss - Z_s*R_ss)
        sin_u[m] = -Z_s/s_p[m]
    end
    


    #---------------------------------------------------------------
    # Compute f=R*Bt such that the eikonal S which solves
    # B*Grad(S)=0 has the correct quasi-periodicity S(s+Ls)=S(s)-2*pi*q_s
    f = 0.0
    for m in 2:ms+1
        f = f + 0.5*ds*(s_p[m-1]/(R[m-1]*psi_x[m-1]) + s_p[m]/(R[m]*psi_x[m]))
    end
    f = 2π*q_s/f

    ### return value
    Bt0_out = f/rmaj_input

    ###### not in the python
    # Bref_out = 1.0
    # betae_s = betae_in
    # debye_s = debye_in
    # if(units_in .eq. 'GENE')then
    # # convert inputs from GENE reference magnetic field to Bunit
    #     Bref_out = f/Rmaj_input # Bref/Bunit
    #     # write(*,*)"Bref/Bunit = ",Bref_out
    #     betae_s = betae_in*Bref_out**2
    #     p_prime_s = p_prime_loc*Bref_out**2
    #     debye_s = debye_in/Bref_out
    # endif



    #*************************************************************
    # Compute toroidal and total fields:
    #*************************************************************
    Bt .= f ./R
    B .= .√(Bt.^2 .+ Bp.^2)


    #*************************************************************
    # # Compute Miller's D0 , Dp and Dff'p needed for kx.
    #*************************************************************
    d_0 = zeros(Float64, ms+1)
    d_p = zeros(Float64, ms+1)
    d_ffp = zeros(Float64, ms+1)
    
    dq1 = ds*s_p[1]*f/(R[1]*psi_x[1]^2)
    d0_s1 = -dq1*(2.0/r_curv[1] + 2.0*sin_u[1]/R[1])
    dp_s1 = dq1*4.0*pi*R[1] /Bp[1]
    dffp_s1 = dq1*(R[1]/Bp[1])*(B[1]/f)^2
    
    for m = 2:ms+1
        dq2 = ds*s_p[m]*f/(R[m]*psi_x[m]^2)
        d0_s2 = -dq2*(2.0/r_curv[m]+2.0*sin_u[m]/R[m])
        dp_s2 = dq2*4π*R[m]/Bp[m]
        dffp_s2 = dq2*(R[m]/Bp[m])*(B[m]/f)^2
        
        d_0[m] = d_0[m-1]+(d0_s1+d0_s2)/2
        d_p[m] = d_p[m-1]+(dp_s1+dp_s2)/2
        d_ffp[m] = d_ffp[m-1]+(dffp_s1+dffp_s2)/2
        
        d0_s1 = d0_s2
        dp_s1 = dp_s2
        dffp_s1 = dffp_s2
    end


    #*************************************************************
    # Begin computing geometric quantities required for solution
    # of gyrokinetic equation:
    #*************************************************************
    b_geo .= B
    pk_geo .= 2.0 .*Bp./B
    qrat_geo .= (rmin_s./R).*(B./Bp)/q_s


    #### not done in the python
    #*************************************************************
    # Determine ff_prime from:
    #*************************************************************
    ff_prime = (2π*q_prime_s-d_0[ms+1]-d_p[ms+1]*p_prime_s)/d_ffp[ms+1]

    #*************************************************************
    # Compute [[kx/ky]] (herein, kxoky_geo) from Waltz-Miller [2] paper.
    #*************************************************************
   
    S_prime = zeros(Float64, ms+1)
    kx_factor = zeros(Float64, ms+1)
    kxoky_geo = zeros(Float64, ms+1)
    for m = 1:ms+1
        S_prime[m] = -(d_0[m]+d_p[m]*p_prime_s+d_ffp[m]*ff_prime)
        kx_factor[m] = (psi_x[m]^2)/B[m]
        kxoky_geo[m] = S_prime[m]*kx_factor[m]
    end


    costheta_geo = zeros(Float64, ms+1)
    #*************************************************************
    # Compute drift coefficients:
    #*************************************************************
    p_prime_zero_s = 1.0
    if(inputs.USE_MHD_RULE) p_prime_zero_s = 0.0 end

    # ### not used
    # epsl_geo = (2.0/rmaj_s).*qrat_geo./b_geo

    costheta_p_geo = (4π*p_prime_s*p_prime_zero_s*rmaj_s) .* (Bp.*R./B.^2)
    costheta_geo .= (
            -rmaj_s.* (Bp./(B.^2)) .*
            ((Bp./r_curv) .- (f^2 ./(Bp.*R.^3) ).*sin_u) )
    

    #*************************************************************
    # Functions which require theta-derivatives:
    #*************************************************************
    sintheta_geo = zeros(Float64, ms+1)
    for m = 1:ms+1

        m1 = ((ms + m - 2) % (ms+1)) + 1
        m2 = ((ms + m - 1) % (ms+1)) + 1
        m3 = (m % (ms+1)) + 1
        m4 = ((m + 1) % (ms+1)) + 1
        ### had to do this weird short-circuiting bc of weird indexing
        m1 < ms || (m1 = ((ms + m - 3) % (ms+1)) + 1)
        m2 < ms+1 || (m2 = ((ms + m - 2) % (ms+1)) + 1)
        m3 > 1 || (m3 = ((m + 1) % (ms+1)) + 1)
        m4 > 2 || (m4 = ((m + 2) % (ms+1)) + 1)

        sintheta_geo[m] = -rmaj_s*(f/(R[m]*B[m]^2)) * (B[m1]-8*B[m2]+8*B[m3]-B[m4])/(delta_s*s_p[m])

    end

    #*************************************************************
    # compute vprime and vpp
    #*************************************************************
    # vprime = 0.0
    # vpp = 0.0
    # dvpp1 = (s_p(0)*R(0)/psi_x(0)**3) &
    #     *(4.0*pi*p_prime_s*R(0)**2 + ff_prime - 2.0*psi_x(0)/r_curv(0))
    # do m=1,ms
    #     dvpp2 = (s_p(m)*R(m)/psi_x(m)**3) &
    #         *(4.0*pi*p_prime_s*R(m)**2 + ff_prime - 2.0*psi_x(m)/r_curv(m))
    #     vprime = vprime + 0.5*ds*(s_p(m-1)/Bp(m-1) + s_p(m)/Bp(m))
    #     vpp = vpp + 0.5*ds*(dvpp1 + dvpp2)
    #     dvpp1 = dvpp2
    # enddo
    # vprime = pi_2*vprime
    # vpp = pi_2*vpp
    # # write(*,*)"vprime = ",vprime,"vpp = ",vpp
    # #
    # ave_M1 = 0.0
    # ave_M2 = 0.0
    # ave_M3 = 0.0
    # ave_M4 = 0.0
    # do m=1,ms
    #     ave_M1 = ave_M1 + 0.5*ds*(s_p(m-1)/Bp(m-1)**3 + s_p(m)/Bp(m)**3)
    #     ave_M2 = ave_M2 + 0.5*ds*(s_p(m-1)*R(m-1)/psi_x(m-1)**3 + s_p(m)*R(m)/psi_x(m)**3)
    #     ave_M3 = ave_M3 + 0.5*ds*(s_p(m-1)*(B(m-1)/psi_x(m-1))**2/Bp(m-1) + s_p(m)*(B(m)/psi_x(m))**2/Bp(m))
    #     ave_M4 = ave_M4 + 0.5*ds*(s_p(m-1)*B(m-1)**2/Bp(m-1) + s_p(m)*B(m)**2/Bp(m))
    # enddo
    # # write(*,*)"M1=",ave_M1,"M2=",ave_M2,"M3=",ave_M3,"M4=",ave_M4
    # p_prime_M = 4.0*pi*p_prime_s
    # q_prime_M = pi_2*q_prime_s
    # DM_out = 0.25 + (p_prime_M/MAX(q_prime_M**2,1e-12))*((vpp/pi_2 - p_prime_M*ave_M1)*ave_M3 &
    #     + (f**2*p_prime_M*ave_M2 - q_prime_M*f)*ave_M2)
    # H = (f*p_prime_M*q_prime_M/MAX(q_prime_M**2,1e-12))*ave_M3*(ave_M2/ave_M3 - vprime/(pi_2*ave_M4))
    # # write(*,*)"H = ",H
    # DR_out = DM_out - (0.5 - H)**2

    return s_p, Bp, b_geo, pk_geo, qrat_geo, sintheta_geo, costheta_geo, costheta_p_geo, Bt0_out, B_unit_out, ds, t_s, f, R, B, S_prime, kx_factor
    
end


function miller_geo(inputs::InputTJLF, mts::Float64=5.0, ms::Int=128)

    rmin_loc = inputs.RMIN_LOC
    rmaj_loc = inputs.RMAJ_LOC
    
    delta_loc = inputs.DELTA_LOC
    kappa_loc = inputs.KAPPA_LOC
    zeta_loc = inputs.ZETA_LOC
    q_loc = inputs.Q_LOC
    #### these might be global variables
    p_prime_s = inputs.P_PRIME_LOC
    q_prime_s = inputs.Q_PRIME_LOC

    drmajdx_loc = inputs.DRMAJDX_LOC
    drmindx_loc = inputs.DRMINDX_LOC


    ### these are set to 0.0 despite having definitions in the inputs
    zmaj_loc = 0.0
    dzmajdx_loc = 0.0
    s_zeta_loc = inputs.S_ZETA_LOC
    s_delta_loc = inputs.S_DELTA_LOC
    s_kappa_loc = inputs.S_KAPPA_LOC

    R = zeros(Float64,ms+1)
    Z = zeros(Float64,ms+1)
    Bp = zeros(Float64,ms+1)

    
    if(rmin_loc<0.00001) rmin_loc=0.00001 end
    #
    # compute the arclength around the flux surface
    #
    x_delta = asin(delta_loc)
    theta = 0.0
    arg_r = theta + x_delta*sin(theta)
    darg_r = 1.0 + x_delta*cos(theta)
    arg_z = theta + zeta_loc*sin(2.0*theta)
    darg_z = 1.0 + zeta_loc*2.0*cos(2.0*theta)
    r_t = -rmin_loc*sin(arg_r)*darg_r
    z_t = kappa_loc*rmin_loc*cos(arg_z)*darg_z
    l_t = √(r_t^2 + z_t^2)
    # scale dtheta by l_t to keep mts points in each ds interval of size pi_2/ms
    dtheta = 2π/(mts*ms*l_t)
    l_t1 = l_t
    scale_max = l_t
    arclength = 0.0
    
    while(theta<2π)
        theta = theta + dtheta
        if(theta>2π)
            theta = theta-dtheta
            dtheta = 2π-theta
            theta = 2π
        end

        arg_r = theta + x_delta*sin(theta)
        darg_r = 1.0 + x_delta*cos(theta) # d(arg_r)/dtheta
        r_t = -rmin_loc*sin(arg_r)*darg_r # dR/dtheta
        
        arg_z = theta + zeta_loc*sin(2.0*theta)
        darg_z = 1.0 + zeta_loc*2.0*cos(2.0*theta) # d(arg_z)/dtheta
        
        z_t = kappa_loc*rmin_loc*cos(arg_z)*darg_z # dZ/dtheta
        l_t = √(r_t^2 + z_t^2) # dl/dtheta
        
        # arclength along flux surface in poloidal direction
        arclength = arclength + 0.50*(l_t + l_t1)*dtheta
        # save maximum expansion scale for later
        if(l_t>scale_max)  scale_max = l_t end

        l_t1 = l_t
    end

    ds = arclength/ms
    
    # Find the theta points which map to an equally spaced s-grid of ms points along the arclength
    # going clockwise from the outboard midplane around the flux surface
    # by searching for the theta where dR**2 + dZ**2 >= ds**2 for a centered difference df=f(m+1)-f(m-1).
    # This keeps the finite difference error of dR/ds, dZ/ds on the s-grid small
    
    t_s = zeros(Float64,ms+1)
    t_s[ms+1]=-2π
    # make a first guess based on theta=0.0
    theta = 0.0

    for m in 1:round(Integer, ms/2)+1
        arg_r = theta + x_delta*sin(theta)
        darg_r = 1.0 + x_delta*cos(theta)
        arg_z = theta + zeta_loc*sin(2.0*theta)
        darg_z = 1.0 + zeta_loc*2.0*cos(2.0*theta)
        r_t = -rmin_loc*sin(arg_r)*darg_r
        z_t = kappa_loc*rmin_loc*cos(arg_z)*darg_z
        l_t = √(r_t^2 + z_t^2)
        if m == 1
            dtheta = -ds/l_t
            theta = dtheta
        else
            dtheta = -ds/(0.5*(l_t+l_t1))
            t_s[m] = t_s[m-1] + dtheta
            theta = t_s[m] + dtheta
        end
        l_t1 = l_t
    end

    # distribute endpoint error over interior points
    dtheta = (t_s[round(Int, ms/2+1)]+π) / (ms/2)

    for m in 2:round(Integer, ms/2)+1
        t_s[m] = t_s[m] - (m-1)*dtheta
        t_s[ms-m+2] = -2π - t_s[m]
    end

    # Loop to compute most geometrical quantities needed for Mercie-Luc expansion
    # R, Z, R*Bp on flux surface s-grid
    B_unit_out = zeros(Float64, ms + 1)
    # grad_r_out = zeros(Float64, ms + 1)
    for m in 1:ms+1
        theta = t_s[m]
        arg_r = theta + x_delta*sin(theta)
        darg_r = 1.0 + x_delta*cos(theta)
        arg_z = theta + zeta_loc*sin(2.0*theta)
        darg_z = 1.0 + zeta_loc*2.0*cos(2.0*theta)

        R[m] = rmaj_loc + rmin_loc*cos(arg_r) # R(theta)
        Z[m] = zmaj_loc + kappa_loc*rmin_loc*sin(arg_z) # Z(theta)
        
        R_t = -rmin_loc*sin(arg_r)*darg_r # dR/dtheta
        Z_t = kappa_loc*rmin_loc*cos(arg_z)*darg_z # dZ/dtheta
        l_t = √(R_t^2 + Z_t^2) # dl/dtheta

        # dR/dr
        R_r = drmajdx_loc + drmindx_loc*cos(arg_r) - sin(arg_r)*s_delta_loc*sin(theta)/√(1.0 - delta_loc^2)
        # dZ/dr
        Z_r = dzmajdx_loc + kappa_loc*sin(arg_z)*(drmindx_loc +s_kappa_loc) + kappa_loc*cos(arg_z)*s_zeta_loc*sin(2.0*theta)
        # Jacobian
        det = R_r*Z_t - R_t*Z_r

        grad_r = abs(l_t/det)
        if m==1
            global B_unit = 1.0/grad_r # B_unit choosen to make bx(0)=ky**2 i.e. qrat_geo(0)/b_geo(0)=1.0
            if(drmindx_loc==1.0) global B_unit=1.0 end # Waltz-Miller convention
        end
        B_unit_out[m] = B_unit
        # grad_r_out[m] = grad_r

        # changes q_s to q_loc
        Bp[m] = (rmin_loc/(q_loc*R[m]))*grad_r*B_unit
        p_prime_s = p_prime_s * B_unit
        q_prime_s = q_prime_s / B_unit
        
    end
    
    return R, Bp, Z, q_prime_s, p_prime_s, B_unit_out, ds, t_s
end