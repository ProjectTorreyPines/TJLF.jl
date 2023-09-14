# export get_sat_params, xgrid_functions_geo, mercier_luc, miller_geo

### get_zonal_mixing()
include("tjlf_modules.jl")
include("tjlf_multiscale_spectrum.jl")

function get_sat_params(inputs::InputTJLF, ky::AbstractVector{T}, gammas::AbstractMatrix{T}) where T<:Real

    kx0_e,
    SAT_geo0_out,
    SAT_geo1_out,
    SAT_geo2_out,
    R_unit,
    Bt0_out,
    B_geo0_out,
    grad_r0_out,
    theta_out,
    Bt_out,
    grad_r_out,
    B_unit_out = xgrid_functions_geo(inputs, ky, gammas)

    return (
        kx0_e, # xgrid_functions_geo
        SAT_geo0_out, # xgrid_functions_geo
        SAT_geo1_out, # xgrid_functions_geo
        SAT_geo2_out, # xgrid_functions_geo
        R_unit, # xgrid_functions_geo

        Bt0_out, # mercier_luc
        B_geo0_out, # mercier_luc
        grad_r0_out, # mercier_luc

        theta_out, # miller_geo
        Bt_out, # miller_geo
        grad_r_out, # miller_geo
        B_unit_out, # miller_geo
    )
end


function get_sat_params(param::Type{Val{:grad_r0}}, inputs::InputTJLF)
    _, _, b_geo, _, qrat_geo, _, _, _, _, _, _, _ = mercier_luc(inputs)
    return b_geo[1]/qrat_geo[1]
end



















#### LINES 220-326, 439-478 in tglf_geometry.f90
function xgrid_functions_geo(inputs::InputTJLF, ky::AbstractVector{T}, gammas::AbstractMatrix{T}, 
    mts::T=5.0, ms::Int=128, small::T=0.00000001) where T<:Real
    #******************************************************************************#************************
    #
    # PURPOSE: compute the geometric coefficients on the x-grid
    #
    #
    #******************************************************************************#************************
    #
    ### different geometries!!!
    rmaj_s = inputs.RMAJ_LOC
    rmin_s = inputs.RMIN_LOC
    q_s = inputs.Q_LOC

    ### needed for this function
    sat_rule_in = inputs.SAT_RULE
    alpha_quench_in = inputs.ALPHA_QUENCH
    alpha_e_in = inputs.ALPHA_E
    sign_IT = inputs.SIGN_IT
    vexb_shear = inputs.VEXB_SHEAR
    mass_2 = inputs.SPECIES[2].MASS
    taus_2 = inputs.SPECIES[2].TAUS
    zs_2 = inputs.SPECIES[2].ZS
    units_in = inputs.UNITS

    vs_2 = √(taus_2 / mass_2)
    gamma_reference_kx0 = gammas[:, 1]

    s_p, Bp, b_geo, pk_geo, qrat_geo, costheta_geo, Bt0_out, B_unit_out, grad_r_out, ds, t_s, B = mercier_luc(inputs)
    
    y = zeros(Float64,ms+1)

    #
    # find length along magnetic field y
    #
    y[1]=0.0
    for m in 2:ms+1
        
        y[m] = y[m-1]+s_p[m]*ds*4.0/(pk_geo[m]+pk_geo[m-1])
    end
    # set the global units
    Ly=y[ms+1]
    R_unit = rmaj_s*b_geo[1]/(qrat_geo[1]*costheta_geo[1])
    q_unit = Ly/(2π*R_unit)

    #### not used in Python
    # # midplane effective shear: reduces to s-alpha in shifted circle
    # # note: S_prime(0)=0.0, S_prime(ms)=-2 pi q_prime, y(0)=0.0, y(ms)=Ly
    # # midplane shear is average of left and right y-derivatives at midplane for general geometry
    # midplane_shear = -(Ly/(2π))*((rmin_s/q_s)^2)*0.5
    #         *(S_prime[2]/y[2]+(S_prime[ms+1]-S_prime[ms])/(y[ms+1]-y[ms]))
    # midplane_shear += 0.11
    # # # save f for output
    # # RBt_ave_out = f/B_unit


    # #
    # # generalized quench rule kx0 shift
    # #
    # kx0 = kx0_loc/ky # note that kx0 is kx/ky
    # kx0_e=0.0
    # kx0_p=0.0
    # #EPS2011 sign_kx0=1.0

    

    if(alpha_quench_in==0.0 && gamma_reference_kx0[1]!=0.0)
        vexb_shear_s = vexb_shear * sign_IT
        vexb_shear_kx0 = alpha_e_in*vexb_shear_s

        kyi = ky.*(vs_2*mass_2/abs(zs_2))
        wE = zeros(Float64, length(ky))
        ### not used 
        wd0 = abs.(ky./rmaj_s)
        
        ### ordering of the conditions is different between fortran and Python
        if(units_in=="GYRO")
            kx0_factor = abs(b_geo[1]/qrat_geo[1]^2)
            kx0_factor = 1.0 + 0.40*(kx0_factor-1.0)^2
            ### wE is an array for Python but not Fotran
            wE .= (kx0_factor*vexb_shear_kx0) .*
                (ifelse.(kyi.<0.3,kyi./0.3,1.0))./gamma_reference_kx0
        else
            kx0_factor = 1.0
        end
        grad_r0_out = b_geo[1]/qrat_geo[1]
        B_geo0_out = b_geo[1]
        kx_geo0_out= 1.0/qrat_geo[1]

        kx0_e = -(0.36*vexb_shear_kx0./gamma_reference_kx0
                .+ (0.38.*wE.*tanh.((0.69.*wE).^6)))



        a0 = 1.3
        if(sat_rule_in==1)
            a0 = 1.45
            kx0_e = -(0.53*vexb_shear_kx0./gamma_reference_kx0 .+
                    (0.25.*wE.*tanh.((0.69.*wE).^6)))
        elseif(sat_rule_in==2 || sat_rule_in==3)
            a0=1.6
            vzf_out, kymax_out, _ = get_zonal_mixing(inputs, ky, gamma_reference_kx0)
            if(abs(kymax_out*vzf_out*vexb_shear_kx0) > small)
                kx0_e = -(0.32*vexb_shear_kx0).*((ky./kymax_out).^0.3)./(ky.*vzf_out)
            else
                kx0_e = zeros(Float64,length(ky))
            end
        end

        kx0_e = ifelse.(abs.(kx0_e) .> a0 , a0.*kx0_e./abs.(kx0_e) , kx0_e)
        ### copied from the Python, why is this a thing
        kx0_e = ifelse.(isnan.(kx0_e) , 0 , kx0_e)

        #### not in the Python
        # if(units_in=="GYRO")
        #     kx0 = sign_Bt_in*kx0_e # cancel the sign_Bt_in factor in kxx below
        # else
        #     if(sat_rule_in.eq.1)kx0 = sign_Bt_in*kx0_e/(2.1)end  # goes with xnu_model=2
        #     if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)kx0 = sign_Bt_in*kx0_e*0.7/grad_r0_out**2 end     # goes with xnu_model=3, the factor 0.7/grad_r0_out**2 is needed for stress_tor
        #     # note kx0 = alpha_e*gamma_ExB_HB/gamma Hahm - Burrell form of gamma_ExB
        #     # The 2.1 effectively increases ay0 & ax0 and reduces toroidal stress to agree with CGYRO
        # end
        
    end

    #####################################################################################################################

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
        # R2x2 = R(i)**2
        # R2_ave_out = R2_ave_out + dlp*(R2x1+R2x2)/2.0
        # B_ave_out = B_ave_out + dlp*(b_geo(i-1)+b_geo(i))/2.0
        # Bt_ave_out = Bt_ave_out + dlp*(f/b_geo(i-1)+f/b_geo(i))/(2.0*Rmaj_s)
        # Grad_r_ave_out = Grad_r_ave_out + dlp*0.5*((R(i-1)*Bp(i-1))**2+(R(i)*Bp(i))**2)*(q_s/rmin_s)**2
        # kykx_geo_ave = kykx_geo_ave + dlp*0.5*(B_geo(i-1)**2/qrat_geo(i-1)**4+B_geo(i)**2/qrat_geo(i)**4)
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
    B_geo0_out = b_geo[1]
    ### so this value is defined in mercier_luc
    Bt0_out = Bt0_out
    B_unit_out = B_unit_out
    grad_r_out = grad_r0_out

    # Additional outputs for SAT2 G1(theta), Gq(theta)
    theta_out = t_s  # theta grid over which everything is calculated.
    Bt_out = B  # total magnetic field matching theta_out grid.

    return kx0_e,
            SAT_geo0_out,
            SAT_geo1_out,
            SAT_geo2_out,
            R_unit,
            Bt0_out,
            B_geo0_out,
            grad_r0_out,
            theta_out,
            Bt_out,
            grad_r_out,
            B_unit_out

end









function mercier_luc(inputs::InputTJLF, mts::AbstractFloat=5.0, ms::Integer=128, small::AbstractFloat=0.00000001)
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
    b_geo = zeros(Real, ms+1)

    ### technically don't have to initialize here, but maybe better looking?
    ### if remove, switch .= to = for the below functions
    Bt = zeros(Real, ms+1)
    B = zeros(Real, ms+1)
    pk_geo = zeros(Real, ms+1)
    qrat_geo = zeros(Real, ms+1)
    costheta_geo = zeros(Real, ms+1)
    

    R, Bp, Z, q_prime_s, p_prime_s, B_unit_out, grad_r_out, ds, t_s = miller_geo(inputs)
    
    
    psi_x = R.*Bp
    delta_s = 12.0*ds
    ds2 = 12.0*ds^2


    # note that the point 1 and ms+1 are the same so m+1->1 and m-1->ms-1 at m=0
    s_p = zeros(Real,ms+1)
    r_curv = zeros(Real,ms+1)
    sin_u = zeros(Real,ms+1)
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



    #-----------------------------------------------------------
    #-----------------------------------------------------------
    # Compute toroidal and total fields:
    Bt .= f ./R
    B .= .√(Bt.^2 .+ Bp.^2)

    #### not done in the python
    #-----------------------------------------------------------
    #-----------------------------------------------------------
    # # Compute Miller's D0 , Dp and Dff'p needed for kx.
    # d_0(0) = 0.0
    # d_p(0) = 0.0
    # d_ffp(0) = 0.0
    # #
    # dq1 = ds*s_p(0)*f/(R(0)*psi_x(0)**2)
    # d0_s1 = -dq1*(2.0/r_curv(0)+2.0*sin_u(0)/R(0))
    # dp_s1 = dq1*4.0*pi*R(0)/Bp(0)
    # dffp_s1 = dq1*(R(0)/Bp(0))*(B(0)/f)**2
    # #
    # do m=1,ms
    #     dq2 = ds*s_p(m)*f/(R(m)*psi_x(m)**2)
    #     d0_s2 = -dq2*(2.0/r_curv(m)+2.0*sin_u(m)/R(m))
    #     dp_s2 = dq2*4.0*pi*R(m)/Bp(m)
    #     dffp_s2 = dq2*(R(m)/Bp(m))*(B(m)/f)**2
    #     #
    #     d_0(m) = d_0(m-1)+0.5*(d0_s1+d0_s2)
    #     d_p(m) = d_p(m-1)+0.5*(dp_s1+dp_s2)
    #     d_ffp(m) = d_ffp(m-1)+0.5*(dffp_s1+dffp_s2)
    #     #
    #     d0_s1 = d0_s2
    #     dp_s1 = dp_s2
    #     dffp_s1 = dffp_s2
    #     #
    # enddo


    #-----------------------------------------------------------
    #-----------------------------------------------------------
    # Begin computing geometric quantities required for solution
    # of gyrokinetic equation:
    #
    # - b_geo replaces bmaj(j)=b_theta(j)/b_unit
    #
    # - pk_geo is close to pk=2*rmin/(rmaj*q), the coefficient
    # of d/dtheta
    #
    # - qrat_geo -> 1 in a circle
    #
    # Note that for the physical quantity
    #
    # k_theta = nq/r
    #
    # we use
    #
    # kyrhos_s = n*q_s/rmin_s*rhos_unit_s
    #
    # which is exactly the same as for the circle.
    #
    # Also, "omega_star" remains unchanged from circle with
    # logarithmic density gradients along minor axis.
    #
    # - "ky*rhos" in Bessel function is kyrhos_s*qrat_geo(j)/b_geo(j)
    #
    # - "kx*rhos" is kxoky_geo(j)*ky*rhos
    #
    
    b_geo .= B
    pk_geo .= 2.0 .*Bp./B
    qrat_geo .= (rmin_s./R).*(B./Bp)/q_s


    #### not done in the python
    #-----------------------------------------------------------
    #---------------------------------------------------------------
    # # Determine ff_prime from:
    # #
    # # 2 pi q_prime = d_0(ms)
    # # +d_p(ms)*p_prime
    # # +d_ffp(ms)*ff_prime
    # #
    # ff_prime = (2π*q_prime_s-d_0(ms)-d_p(ms)*p_prime_s) &
    #     /d_ffp(ms)

    #### not done in the python
    #---------------------------------------------------------------
    #--------------------------------------------------------------
    # # Compute [[kx/ky]] (herein, kxoky_geo) from Waltz-Miller [2] paper.
    # # 2
    # # (R B_p) S1
    # # kxoky_geo = ---------- ------
    # # B R B_p
    # #
    # # S1
    # # ------ = -(d_0(theta)+d_p(theta)*p_prime+d_ffp(theta)*ff_prime)
    # # R B_p
    # #
    # do m=0,ms

    #     S_prime(m) = -(d_0(m)+d_p(m)*p_prime_s+d_ffp(m)*ff_prime)
    #     kx_factor(m) = (psi_x(m)**2)/B(m)
    #     kxoky_geo(m) = S_prime(m)*kx_factor(m)
    #     # write(*,*)"check s_prime",S_prime(m)

    # enddo


    #---------------------------------------------------------------
    #---------------------------------------------------------------
    # # Compute drift coefficients:
    # # p_prime_zero forces grad-B-curvature to zero to compensates
    # # for b_par =0
    # #
    # p_prime_zero_s = 1.0
    # if(use_mhd_rule_in) p_prime_zero_s = 0.0 end

    # ### not used
    # epsl_geo = (2.0/rmaj_s).*qrat_geo./b_geo
    # ### not used
    # costheta_p_geo = (4.0*π*p_prime_s*p_prime_zero_s*rmaj_s) .* (Bp.*R./B.^2)
    costheta_geo .= (
            -rmaj_s.* (Bp./(B.^2)) .*
            ((Bp./r_curv) .- (f^2 ./(Bp.*R.^3) ).*sin_u) )
    

    #### not done in the python
    #---------------------------------------------------------------
    #-------------------------------------------------------------
    # Functions which require theta-derivatives:
    #
    # do m=0,ms

    #     # Waltz/Miller [[sin]]
    #     m1=MOD(ms+m-2,ms)
    #     m2=MOD(ms+m-1,ms)
    #     m3=MOD(m+1,ms)
    #     m4=MOD(m+2,ms)
    #     sintheta_geo(m) = -rmaj_s*(f/(R(m)*B(m)**2))* &
    #         (B(m1)-8.0*B(m2)+8.0*B(m3)-B(m4))/(delta_s*s_p(m))
    #     # write(*,*)m,m1,m2,m3,m4,"sintheta_geo=",sintheta_geo(m)

    # enddo
    # #
    # # compute vprime and vpp
    # #
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

    return s_p, Bp, b_geo, pk_geo, qrat_geo, costheta_geo, Bt0_out, B_unit_out, grad_r_out, ds, t_s, B
    
end


function miller_geo(inputs::InputTJLF, mts::AbstractFloat=5.0, ms::Integer=128)

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

    R = zeros(Real,ms+1)
    Z = zeros(Real,ms+1)
    Bp = zeros(Real,ms+1)

    
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
    
    t_s = zeros(Real,ms+1)
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
    B_unit_out = zeros(Real, ms + 1)
    grad_r_out = zeros(Real, ms + 1)
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
        grad_r_out[m] = grad_r

        # changes q_s to q_loc
        Bp[m] = (rmin_loc/(q_loc*R[m]))*grad_r*B_unit
        p_prime_s = p_prime_s * B_unit
        q_prime_s = q_prime_s / B_unit
        
    end
    
    return R, Bp, Z, q_prime_s, p_prime_s, B_unit_out, grad_r_out, ds, t_s
end