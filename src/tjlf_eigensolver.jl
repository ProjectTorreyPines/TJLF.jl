import LinearAlgebra.LAPACK.ggev!
include("tjlf_modules.jl")

function tjlf_eigensolver(inputs::InputTJLF)

    #   CHARACTER(1) :: rightvectors
    #   INTEGER :: is,js
    #   INTEGER ::  j1, j2, j, i
    #   INTEGER :: ifail,ib, jb, ia, ja, ia0, ja0
    #   INTEGER :: info,lwork
    #   INTEGER :: xnu_model
    #   INTEGER :: idum=0
    #   REAL :: ft2,ft3,ft4,ft5
    #   REAL :: Linsker,am,bm
    #   REAL :: w_s, w_dh, w_dg, w_d1, modw_d1,w_d0, w_cd
    #   REAl :: wd_psi,modwd_psi
    #   REAL :: E_i,M_i,N_j,J_j
    #   REAL :: k_par0,k_par1,modk_par0,modk_par1
    #   REAL :: h10n,h10p1,h10p3,h10r13,h10r33
    #   REAL :: hnb0,hp1b0,hp3b0,hr11b0,hr13b0,hr33b0
    #   REAL :: hw113b0,hw133b0,hw333b0
    #   REAL :: hnbp,hp1bp,hp3bp,hr11bp,hr13bp,hr33bp
    #   REAL :: hw113bp,hw133bp,hw333bp
    #   REAL :: c_tor_par_hp1,c_tor_par_hr11
    #   REAL :: c_tor_par_hr13
    #   REAL :: g10n,g10p1,g10p3,g10r13,g10r33
    #   REAL :: gnb0,gp1b0,gp3b0,gr11b0,gr13b0,gr33b0
    #   REAL :: gw113b0,gw133b0,gw333b0
    #   REAL :: gnbp,gp1bp,gp3bp,gr11bp,gr13bp,gr33bp
    #   REAL :: gw113bp,gw133bp,gw333bp
    #   REAL :: c_tor_par_gp1,c_tor_par_gr11
    #   REAL :: c_tor_par_gr13
    #   REAL :: betae_psi,betae_sig
    #   REAL :: max_freq,test
    #   REAL :: gradB1
    #   REAL :: c35
    #   REAL :: xnuei,d_ab,d_ij,d_ee,d_ab_psi,d_11,d_1
    #   REAL :: xnuion
    #   REAL :: xnu_bndry,xnu_hat,xnu_a,xnu_b,xnu_c
    #   REAL :: xnu_phi_b
    #   REAL :: xnu_n_b
    #   REAL :: xnu_p1_b
    #   REAL :: xnu_p3_b
    #   REAL :: xnu_u_b
    #   REAL :: xnu_q1_b
    #   REAL :: xnu_q3_b
    #   REAL :: xnu_p1_1
    #   REAL :: xnu_u_u_1,xnu_u_q3_1
    #   REAL :: xnu_q1_u_1,xnu_q1_q1_1,xnu_q1_q3_1
    #   REAL :: xnu_q3_u_1,xnu_q3_q3_1
    #   REAL :: c01,c02,c03,c04,c05
    #   REAL :: c06,c07,c08,c09,c010
    #   REAL :: cb1,cb2,cb3,cb4,cb5,cb6,cb7,cb8
    #   REAL :: an,ap3,ap1,bn,bp3,bp1
    #   REAL :: cnuei,kparvthe,cxf
    #   REAL :: nuei_p1_p1_1,nuei_p1_p3_1
    #   REAL :: nuei_u_u_1,nuei_u_q3_1
    #   REAL :: nuei_q1_u_1,nuei_q1_q1_1,nuei_q1_q3_1
    #   REAL :: nuei_q3_u_1,nuei_q3_q3_1
    #   REAL :: nuei_p1_p1_t,nuei_p1_p3_t
    #   REAL :: nuei_u_u_t,nuei_u_q3_t
    #   REAL :: nuei_q1_u_t,nuei_q1_q1_t,nuei_q1_q3_t
    #   REAL :: nuei_q3_u_t,nuei_q3_q3_t
    #   REAL :: nuei_n_n,nuei_n_p1,nuei_n_p3
    #   REAL :: nuei_p3_n,nuei_p3_p1,nuei_p3_p3
    #   REAL :: nuei_p1_n,nuei_p1_p3,nuei_p1_p1
    #   REAL :: nuei_u_u,nuei_u_q1,nuei_u_q3
    #   REAL :: nuei_q1_u,nuei_q1_q1,nuei_q1_q3
    #   REAL :: nuei_q3_u,nuei_q3_q1,nuei_q3_q3
    #   REAL :: k1,k2,k3,k4,k5
    #   REAL :: ki,ks0,charge_tot,gradne,gradne_s
    #   REAL :: beta2,bs
    #   REAL :: damp_psi,damp_sig
    #   REAL,DIMENSION(nsm) :: vpar, vpar_shear
    #   COMPLEX :: phi_A,phi_B,phi_AU,phi_BU
    #   COMPLEX :: psi_A,psi_B,psi_AN,psi_BN  
    #   COMPLEX :: sig_A,sig_B
    #   COMPLEX :: k_par,k_par_psi
    #   REAL,ALLOCATABLE,DIMENSION(:) :: rwork
    #   COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: at,bt
    #   COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: vleft,vright
    #   COMPLEX,ALLOCATABLE,DIMENSION(:) :: work
    #   COMPLEX,ALLOCATABLE,DIMENSION(:) :: zomega
#
#     ifail = 0
#     lwork = 8*iur
#     ALLOCATE(rwork(lwork))
#     ALLOCATE(at(iur,iur))
#     ALLOCATE(bt(iur,iur))
#     ALLOCATE(vleft(iur,iur))
#     ALLOCATE(vright(iur,iur))
#     lwork = 33*iur
#     ALLOCATE(work(lwork))
#     ALLOCATE(zomega(iur))


#     c35 = 3.0/5.0
#     ft = fts(1)  # electrons
#     ft2 = ft*ft
#     ft3 = ft*ft2
#     ft4 = ft*ft3
#     ft5 = ft*ft4

#     ############################################
#     width_in = inputs.WIDTH
#     use_bpar_in = inputs.USE_BPAR
#     use_bper_in = inputs.USE_BPER
#     vpar_model_in = inputs.VPAR_MODEL
#     damp_psi_in = inputs.DAMP_PSI

#     k_par0 = park_in/(R_unit*q_unit*width_in)
#     w_d0 = ky/R_unit
#     w_cd = -gchat_in*w_d0
#     w_s = -ky/B_unit
#     wd_psi = w_cd
#     modwd_psi = abs(wd_psi)
#     betae_psi = 0.0
#     damp_psi = 0.0

#     if(!use_bper_in)
#         hnb0 = 0.0
#         hp1b0 = 0.0
#         hp3b0 = 0.0
#         hr11b0 = 0.0
#         hr13b0 = 0.0
#         hr33b0 = 0.0
#         hw113b0 = 0.0
#         hw133b0 = 0.0
#         hw333b0 = 0.0
#         kpar_hp1b0 = 0.0
#         kpar_hr11b0 = 0.0
#         kpar_hr13b0 = 0.0
#         wdhp1b0 = 0.0
#         wdhr11b0 = 0.0
#         wdhr13b0 = 0.0
#         gnb0 = 0.0
#         gp1b0 = 0.0
#         gp3b0 = 0.0
#         gr11b0 = 0.0
#         gr13b0 = 0.0
#         gr33b0 = 0.0
#         gw113b0 = 0.0
#         gw133b0 = 0.0
#         gw333b0 = 0.0
#         kpar_gp1b0 = 0.0
#         kpar_gr11b0 = 0.0
#         kpar_gr13b0 = 0.0
#         wdgp1b0 = 0.0
#         wdgr11b0 = 0.0
#         wdgr13b0 = 0.0
#     end
#     if(!use_bpar_in)
#         h10n = 0.0
#         h10p1 = 0.0
#         h10p3 = 0.0
#         h10r13 = 0.0
#         h10r33 = 0.0
#         g10n = 0.0
#         g10p1 = 0.0
#         g10p3 = 0.0
#         g10r13 = 0.0
#         g10r33 = 0.0
#     end
#     if(use_bper_in)
#         if(nbasis==2)
#             betae_psi = 0.5*betae_s/(ky^2+(damp_psi_in*vs(2)/(q_unit*width_in))^2)
#         else
#            betae_psi = 0.5*betae_s/ky^2
#         end
#     endif
#     if(vpar_model_in!=0)
#         hnbp = 0.0
#         hp1bp = 0.0
#         hp3bp = 0.0
#         hr11bp = 0.0
#         hr13bp = 0.0
#         hr33bp = 0.0
#         hw113bp = 0.0
#         hw133bp = 0.0
#         hw333bp = 0.0
#         wdhp1bp = 0.0
#         wdhr11bp = 0.0
#         wdhr13bp = 0.0
#         kpar_hnbp = 0.0
#         kpar_hp1bp = 0.0
#         kpar_hp3bp = 0.0
#         kpar_hr11bp = 0.0
#         kpar_hr13bp = 0.0
#         gnbp = 0.0
#         gp1bp = 0.0
#         gp3bp = 0.0
#         gr11bp = 0.0
#         gr13bp = 0.0
#         gr33bp = 0.0
#         gw113bp = 0.0
#         gw133bp = 0.0
#         gw333bp = 0.0
#         wdgp1bp = 0.0
#         wdgr11bp = 0.0
#         wdgr13bp = 0.0
#         kpar_gnbp = 0.0
#         kpar_gp1bp = 0.0
#         kpar_gp3bp = 0.0
#         kpar_gr11bp = 0.0
#         kpar_gr13bp = 0.0
#     end

#     max_freq=2.0*ABS(ave_wdh(1,1))/R_unit
#     do is=ns0,ns
#         test = ABS(as(is)*zs(is)*(ave_hp3p0(is,1,1)*rlns(is)+1.5*(ave_hr13p0(is,1,1)-ave_hp3p0(is,1,1))*rlts(is)))
#         max_freq = MAX(max_freq,test)
#     enddo 
#     max_freq = filter_in*ABS(ky)*max_freq


#     betae_sig = 0.0
#     damp_sig = 0.0
#     if(use_bpar_in)then
#         if(nbasis.eq.2)then
#             betae_psi = 0.5*betae_s/(ky*ky+(damp_sig_in*vs(2)/(q_unit*width_in))**2)
#         else
#             betae_psi = 0.5*betae_s/(ky*ky)
#         endif
#     endif

#     Linsker = 0.5*Linsker_factor_in
#     if(nbasis.eq.1)Linsker=0.0
#         am = 1.0
#         bm = 0.0
#         if(Linsker.ne.0.0)then
#             am = am/2.0
#             bm = bm/2.0
#         endif
# #
# #  GLF toroidal closure coefficients
# #
#        call get_v
# #
# # GLF parallel closure coefficients
# #

# ############ top of SECTION III and between eqn 22/23 ###########
#         bpar_HP = 3.0+(32.0 - 9.0*pi)/(3.0*pi - 8.0)
#         bper_DH = 1.0
#     # note sqrt(2) included in definition of d_per,d_par
#         dper_DH = sqrt_two*sqrt_pi/2.0
#         dpar_HP = 2.0*sqrt_two*sqrt_pi/(3.0*pi - 8.0)
#         b1 = bpar_HP
#         d1 = dpar_HP
#         b3 = bper_DH
#         d3 = dper_DH
#         b33 = (b1 - b3)/(3.0)
#         d33 = (d1 - d3)/(3.0)

#         # include R(theta)/R0 factor like gyro convetions. Note that sign_Bt_in is in ave_c_tor_par
#         do is=1,ns
#             vpar_shear(is)=vpar_shear_s(is)*ave_c_tor_par(1,1)/Rmaj_input
#             if(vpar_model_in.eq.0)vpar(is) = vpar_s(is)*ave_c_tor_par(1,1)/Rmaj_input
#         enddo
#         #
#         #  compute electron-ion collsion model coefficients
#         #
#         xnu_model = xnu_model_in
#         # different trapped/passing boundary layer models
#         # xnu_model = 0  version 1.80 large xnu limit = 0.0
#         # xnu_model = 1  version 1.81 used for APS07 , large xnu limit = adiabatic
#         # xnu_model = 2  version 1.85 large xnu_limit = circulating response 
#         # xnu_model = 3  retuned trapped boundary term to fit CGYRO with Lorentz operator2/8/2017 & with wdia_trapped 9/15/20
#         # xnu_model = 4  best fit to response function Phys. Plasmas 17, (2010) 122309.
#         k1=0.0
#         k2=0.0
#         k3=0.0
#         k4=0.0
#         k5=0.0
#         if(xnu_model.le.1)then
#             k1 = 0.601248 + 1.12838*zeff_in
#             k2 = 0.797885 + 1.12838*zeff_in
#             k3 = 1.795240 + 2.25676*zeff_in
#         endif
#         xnu_p1_1 = -(4.0/5.0)*k2
#         xnu_u_u_1 = -(2.0/3.0)*k1
#         xnu_u_q3_1 = k1 - (2.0/5.0)*k2
#         xnu_q1_u_1 = -(4.0/5.0)*k2
#         xnu_q1_q1_1 = -(16.0/35.0)*k3
#         xnu_q1_q3_1 = (6.0/5.0)*k2 + (12.0/35.0)*k3
#         xnu_q3_u_1 = -(4.0/9.0)*k2    
#         xnu_q3_q3_1 = (2.0/3.0)*k2 - (4.0/15.0)*k3

#         if(park_in.eq.0.0)then
#             xnu_u_u_1 = 0.0
#             xnu_u_q3_1 =0.0
#             xnu_q1_u_1=0.0
#             xnu_q1_q1_1=0.0
#             xnu_q1_q3_1=0.0
#             xnu_q3_u_1=0.0
#             xnu_q3_q3_1=0.0
#         endif
#         if(nroot.gt.6)then
#             xnu_hat = xnue_s/(ky*taus(1)/R_unit)

#             # model for effective ion wavenumber averaged with ion charge densities zs*as
#             # ki = k_theta*sqrt(sum_s(rho_s**2*as*zs))

#             ki=0.0
#             charge_tot=0.0
#             do is=2,ns
#                 ki = ki + taus(is)*mass(is)*as(is)*zs(is)
#                 charge_tot = charge_tot + as(is)*zs(is)
#             enddo
#             ki = SQRT(ki/charge_tot)*ky_s
#             ks0 = ky*SQRT(taus(1)*mass(2))


#             xnu_phi_b=0.0
#             if(xnu_model.eq.1)xnu_phi_b=1.0
#             if(xnu_phi_b.eq.0.0)then
#                 # boundary collision model without phi terms
#                 gradne = rlns(1)*R_unit
#                 gradne_s = MAX(gradne+10.8,1.8)
#                 xnu_c = gradne_s*(1.5*(1.0-TANH((gradne_s/12.6)**2))+0.13)

#                 xnu_a = ki*(MAX(0.36+0.10*gradne,0.0) + xnu_c*(ki/ks0)*(1.0-TANH(ki/0.55)))
#                 xnu_b = 3.1/(1.0+(2.1*ki + 8.0*ki*ki)*xnu_hat)
#             else
#                 # boundary collsion model with phi terms
#                 xnu_a = 0.41 + ((ki/ks0)**1.7)*0.70*(1.0 + 1.4*ki/0.38)/(1.0 + (ki/0.38)**4)
#                 xnu_b = 3.79/(1.0 + 4.63*ki*xnu_hat)
#             endif
#             xnu_bndry = (1.0 - ft2)*xnu_a*xnu_b

#             xnu_n_b = xnu_factor_in*xnu_bndry*xnu_p1_1
#             xnu_p3_b =  xnu_n_b
#             xnu_p1_b = xnu_n_b
#             xnu_u_b = xnu_factor_in*xnu_bndry*xnu_q1_q1_1 
#             xnu_q3_b = xnu_u_b
#             xnu_q1_b = xnu_u_b
#             if(park_in.eq.0.0)then
#                 xnu_u_b = 0.0
#                 xnu_q1_b=0.0
#                 xnu_q3_b =0.0
#             endif
#         endif
        
#         cnuei = 0.0
#         if(xnu_model.ge.2)cnuei = xnue_s
#         kparvthe = ABS(k_par0)*vs(1)/sqrt_two
#         kparvthe=MAX(kparvthe,1.0E-10)

# #
# # full velocity space terms 
# #
#         k1 = cnuei*(nuei_c1(1) +zeff_in*nuei_c1(2))
#         k2 = cnuei*(nuei_c1(3) +zeff_in*nuei_c1(4))
#         k3 = cnuei*(nuei_c1(5) +zeff_in*nuei_c1(6))
#         k4 = cnuei*(nuei_c1(7) +zeff_in*nuei_c1(8))
#         k5 = cnuei*(nuei_c1(9) +zeff_in*nuei_c1(10))
#         c01 = (4.0/5.0)*k4
#         c02 = (2.0/5.0)*k2 - k1
#         c03 = (2.0/3.0)*k1
#         c04 = (4.0/15.0)*k3-(4.0/3.0)*k2+(5.0/3.0)*k1
#         c05 = (16.0/35.0)*k5
#         c06 = 0.0
#         c07 = 0.0
#         c08 = 0.0
#         c09 = 0.0
#         c010 = 0.0

#         nuei_p1_p1_1 = c01
#         nuei_p1_p3_1 = -nuei_p1_p1_1
#         nuei_u_q3_1 = c02
#         nuei_u_u_1 = c03 -(5.0/3.0)*nuei_u_q3_1
#         nuei_q3_q3_1 = c04
#         nuei_q3_u_1 = (10.0/9.0)*c02 -(5.0/3.0)*nuei_q3_q3_1
#         nuei_q3_u_1 = nuei_q3_u_1 + (5.0/3.0)*nuei_u_u_1
#         nuei_q3_q3_1 = nuei_q3_q3_1 + (5.0/3.0)*nuei_u_q3_1
#         nuei_q1_q1_1 = c05
#         nuei_q1_q3_1 = -(9.0/5.0)*nuei_q1_q1_1
#         nuei_q1_u_1 = (9.0/5.0)*nuei_q3_u_1
#         nuei_q1_q3_1 = nuei_q1_q3_1 + (9.0/5.0)*nuei_q3_q3_1
        
#         nuei_p1_p1_t = c01
#         nuei_p1_p3_t = -ft2*nuei_p1_p1_t
#         nuei_u_q3_t = c02
#         nuei_u_u_t = c03 -(5.0/3.0)*nuei_u_q3_t
#         nuei_q3_q3_t = c04
#         nuei_q3_u_t = (10.0/9.0)*c02 -(5.0/3.0)*nuei_q3_q3_t
#         nuei_q3_u_t = nuei_q3_u_t + (5.0/3.0)*nuei_u_u_t
#         nuei_q3_q3_t = nuei_q3_q3_t + (5.0/3.0)*nuei_u_q3_t
#         nuei_q1_q1_t = c05
#         nuei_q1_q3_t = -(9.0/5.0)*ft2*nuei_q1_q1_t
#         nuei_q1_u_t = (9.0/5.0)*ft2*nuei_q3_u_t
#         nuei_q1_q3_t = nuei_q1_q3_t + (9.0/5.0)*ft2*nuei_q3_q3_t
        
#         an = 0.75 
#         ap3 = 1.25 
#         ap1 = 2.25 
#         bn = ft
#         bp3 = ft
#         bp1 = ft3
#         cb3=0.0
#         cb5 =0.0
# #recalibrated 8/20/14      
#         cb1 = 0.163*SQRT(kparvthe*cnuei*(1.0 + 0.82*zeff_in))
#         if(xnu_model_in.eq.3)then
#             if(wdia_trapped_in.eq.0.0)then
#                 cb1 = 0.50*(kparvthe**0.34)*(cnuei*(1.0 + 0.82*zeff_in))**0.66
#             else
#                 cb1 = 0.315*(kparvthe**0.34)*(cnuei*(1.0 + 0.82*zeff_in))**0.66
#             endif
#         endif
#         cb1 = cb1*xnu_factor_in
#         cb2 = cb1
#         cb4 = cb1
        
#         # even trapped region terms
#         nuei_n_n = (1.0 -ft2)*cb1
#         nuei_n_p3 = (1.0 - ft2)*cnuei*cb3
#         nuei_n_p1 = 0.0
#         nuei_n_p3 = nuei_n_p3 - ft2*nuei_n_p1
#         nuei_n_n = nuei_n_n -nuei_n_p3 - ft2*nuei_n_p1
#         nuei_p3_n = (2.0/3.0)*(1.0 - ft2)*cnuei*cb3
#         nuei_p3_p3 = (1.0 -ft2)*cb2
#         nuei_p3_p1 = 0.0
#         nuei_p3_p3 = nuei_p3_p3 - ft2*nuei_p3_p1
#         nuei_p3_n = nuei_p3_n -nuei_p3_p3 -ft2*nuei_p3_p1
#         nuei_p3_n = nuei_p3_n +nuei_n_n
#         nuei_p3_p3 = nuei_p3_p3 +nuei_n_p3
#         nuei_p3_p1 = nuei_p3_p1 +nuei_n_p1
#         nuei_p1_n = 0.0
#         nuei_p1_p3 = 0.0
#         nuei_p1_p1 = (1.0 - ft2)*cb4
#         nuei_p1_p3 = nuei_p1_p3 -ft2*nuei_p1_p1
#         nuei_p1_n = nuei_p1_n -nuei_p1_p3 -ft2*nuei_p1_p1
#         nuei_p1_n = nuei_p1_n + ft2*nuei_p3_n
#         nuei_p1_p3 = nuei_p1_p3 +ft2*nuei_p3_p3
#         nuei_p1_p1 = nuei_p1_p1 +ft2*nuei_p3_p1      

#         cb5=0.0
#         cb6=0.0
#         cb7=0.0
#         cb8=0.0
#         nuei_u_u = (1.0 -ft2)*cnuei*cb5
#         nuei_u_q3 = (1.0 -ft2)*cnuei*cb7
        
#         nuei_u_q1 = 0.0
#         nuei_u_q3 = nuei_u_q3 - (9.0/5.0)*ft2*nuei_u_q1
#         nuei_u_u = nuei_u_u -(5.0/3.0)*nuei_u_q3 -3.0*ft2*nuei_u_q1
#         nuei_q3_u = (10.0/9.0)*(1.0 -ft2)*cnuei*cb7
#         nuei_q3_q3 = (1.0 - ft2)*cnuei*cb6
        
#         nuei_q3_q1 = 0.0
#         nuei_q3_q3 = nuei_q3_q3 -(9.0/5.0)*ft2*nuei_q3_q1
#         nuei_q3_u = nuei_q3_u - (5.0/3.0)*nuei_q3_q3 -3.0*ft2*nuei_q3_q1
#         nuei_q3_u = nuei_q3_u + (5.0/3.0)*nuei_u_u
#         nuei_q3_q3 = nuei_q3_q3 +(5.0/3.0)*nuei_u_q3
#         nuei_q3_q1 = nuei_q3_q1 +(5.0/3.0)*nuei_u_q1
#         nuei_q1_u = 0.0
#         nuei_q1_q3 = 0.0
#         nuei_q1_q1 = (1.0 -ft2)*cnuei*cb8
#         nuei_q1_q3 = nuei_q1_q3 - (9.0/5.0)*ft2*nuei_q1_q1
#         nuei_q1_u = nuei_q1_u - (5.0/3.0)*nuei_q1_q3 - 3.0*ft2*nuei_q1_q1
#         nuei_q1_u = nuei_q1_u +(9.0/5.0)*ft2*nuei_q3_u
#         nuei_q1_q3 = nuei_q1_q3 +(9.0/5.0)*ft2*nuei_q3_q3
#         nuei_q1_q1 = nuei_q1_q1 + (9.0/5.0)*ft2*nuei_q3_q1


# # done with electron-ion collision model
# #
# #
# # start of loop over species is,js for amat
# #
#         do is = ns0,ns

#             ft = fts(is)
#             ft2 = ft*ft
#             ft3 = ft*ft2
#             ft4 = ft*ft3
#             ft5 = ft*ft4
#             if(nroot.gt.6)then
#                 call get_u(ft)
#                 u2_r = u2_r*ft2
#                 u2_i = u2_i*ft2
#                 u3_r = u3_r/ft2
#                 u3_i = u3_i/ft2
#                 u5_r = u5_r*ft2
#                 u5_i = u5_i*ft2
#                 u7_r = u7_r*ft2
#                 u7_i = u7_i*ft2
#                 u9_r = u9_r/ft2
#                 u9_i = u9_i/ft2
#                 ub2_r = ub2_r*ft2
#                 ub2_i = ub2_i*ft2
#                 ub3_r = ub3_r/ft2
#                 ub3_i = ub3_i/ft2
#                 ub5_r = ub5_r*ft2
#                 ub5_i = ub5_i*ft2
#                 ub7_r = ub7_r*ft2
#                 ub7_i = ub7_i*ft2
#                 ub9_r = ub9_r/ft2
#                 ub9_i = ub9_i/ft2
#             endif
#             do js = ns0,ns
# #
# # start of loop over basis ib,jb for amat
# #
#                 do ib = 1,nbasis
#                     do jb = 1,nbasis
#                         if(idum.eq.1)write(*,*)"fooled you"
# #
# # collision  model
# #
#                         xnuei = 0.0
#                         if(is.eq.1)then
#                         xnuei = xnue_s
#                         endif
#                         xnuion = 0.0
#                         d_ab=0.0
#                         d_ij=0.0
#                         d_ee=0.0
#                         d_ab_psi=0.0
#                         d_11 = 0.0
#                         d_1 = 0.0
#                         if(is.eq.1)d_1=1.0
#                         if(is.eq.1.and.js.eq.1)d_11=1.0
#                         if(is.eq.js)d_ij=1.0
#                         if(ib.eq.jb.and.is.eq.js)d_ab=1.0
#                         if(is.eq.1.and.d_ab.eq.1.0)d_ee=1.0
#                         if(ib.eq.jb)d_ab_psi=1.0
#                         b1 = bpar_HP
#                         d1 = dpar_HP
#                         b3 = bper_DH
#                         d3 = dper_DH
#                         bs = 20.0*SQRT(taus(is)*mass(is))*ky/ABS(zs(is))
#                         b33 = (b1 - b3)/(3.0)
#                         d33 = (d1 - d3)/(3.0)
#                 #
#                         hn = ave_hnp0(is,ib,jb)
#                         hp1 = ave_hp1p0(is,ib,jb)
#                         hp3 = ave_hp3p0(is,ib,jb)
#                         hr11 = ave_hr11p0(is,ib,jb)
#                         hr13 = ave_hr13p0(is,ib,jb)
#                         hr33 = ave_hr33p0(is,ib,jb)
#                         c_tor_par_hp1 = ave_c_tor_par_hp1p0(is,ib,jb)/Rmaj_input
#                         c_tor_par_hr11 = ave_c_tor_par_hr11p0(is,ib,jb)/Rmaj_input
#                         c_tor_par_hr13 = ave_c_tor_par_hr13p0(is,ib,jb)/Rmaj_input

#                         if(use_bper_in)then
#                             hnb0 = ave_hnb0(is,ib,jb)
#                             hp1b0 = ave_hp1b0(is,ib,jb)
#                             hp3b0 = ave_hp3b0(is,ib,jb)
#                             hr11b0 = ave_hr11b0(is,ib,jb)
#                             hr13b0 = ave_hr13b0(is,ib,jb)
#                             hr33b0 = ave_hr33b0(is,ib,jb)
#                             hw113b0 = ave_hw113b0(is,ib,jb)
#                             hw133b0 = ave_hw133b0(is,ib,jb)
#                             hw333b0 = ave_hw333b0(is,ib,jb)
#                             if(vpar_model_in.eq.0)then
#                                 hnbp = ave_hnbp(is,ib,jb)
#                                 hp1bp = ave_hp1bp(is,ib,jb)
#                                 hp3bp = ave_hp3bp(is,ib,jb)
#                                 hr11bp = ave_hr11bp(is,ib,jb)
#                                 hr13bp = ave_hr13bp(is,ib,jb)
#                                 hr33bp = ave_hr33bp(is,ib,jb)
#                                 hw113bp = ave_hw113bp(is,ib,jb)
#                                 hw133bp = ave_hw133bp(is,ib,jb)
#                                 hw333bp = ave_hw333bp(is,ib,jb)
#                             endif
#                         endif
#                         if(use_bpar_in)then
#                             h10n = 1.5*(hnb0-hp3b0)
#                             h10p1 = 2.5*hp1b0 - 1.5*hr13b0
#                             h10p3 = 2.5*hp3b0 - 1.5*hr33b0
#                             h10r13 = 3.5*hr13b0 - 1.5*hw133b0
#                             h10r33 = 3.5*hr33b0 - 1.5*hw333b0
#                         endif
#                         hu1 = ave_hu1(is,ib,jb)
#                         hu3 = ave_hu3(is,ib,jb)
#                         ht1 = ave_ht1(is,ib,jb)
#                         ht3 = ave_ht3(is,ib,jb)
#                         wdhp1p0 = ave_wdhp1p0(is,ib,jb)
#                         wdhr11p0 = ave_wdhr11p0(is,ib,jb)
#                         wdhr13p0 = ave_wdhr13p0(is,ib,jb)
#                         if(use_bper_in)then
#                             wdhp1b0 = ave_wdhp1b0(is,ib,jb)
#                             wdhr11b0 = ave_wdhr11b0(is,ib,jb)
#                             wdhr13b0 = ave_wdhr13b0(is,ib,jb)
#                             if(vpar_model_in.eq.0)then
#                                 wdhp1bp = ave_wdhp1bp(is,ib,jb)
#                                 wdhr11bp = ave_wdhr11bp(is,ib,jb)
#                                 wdhr13bp = ave_wdhr13bp(is,ib,jb)
#                             endif
#                         endif
#                         wdhu1 = ave_wdhu1(is,ib,jb)
#                         wdhu3 = ave_wdhu3(is,ib,jb)
#                         wdhu3ht1 = ave_wdhu3ht1(is,ib,jb)
#                         wdhu3ht3 = ave_wdhu3ht3(is,ib,jb)
#                         wdhu33 = ave_wdhu33(is,ib,jb)
#                         wdhu33ht1 = ave_wdhu33ht1(is,ib,jb)
#                         wdhu33ht3 = ave_wdhu33ht3(is,ib,jb)
#                         modwdhu1 = ave_modwdhu1(is,ib,jb)
#                         modwdhu3 = ave_modwdhu3(is,ib,jb)
#                         modwdhu3ht1 = ave_modwdhu3ht1(is,ib,jb)
#                         modwdhu3ht3 = ave_modwdhu3ht3(is,ib,jb)
#                         modwdhu33 = ave_modwdhu33(is,ib,jb)
#                         modwdhu33ht1 = ave_modwdhu33ht1(is,ib,jb)
#                         modwdhu33ht3 = ave_modwdhu33ht3(is,ib,jb)
#                         hv1r = (v1_r-vb1_r)*ave_modwdh(ib,jb)+vb1_r*modwdhu3*c35
#                         hv2r = (v2_r-vb2_r)*ave_modwdh(ib,jb)+vb2_r*modwdhu3*c35
#                         hv3r = (v3_r-vb3_r)*ave_modwdh(ib,jb)+vb3_r*modwdhu33*c35
#                         hv4r = (v4_r-vb4_r)*ave_modwdh(ib,jb)+vb4_r*modwdhu33*c35
#                         hv5r = (v5_r-vb5_r)*ave_modwdh(ib,jb)+vb5_r*modwdhu3*c35
#                         hv6r = (v6_r-vb6_r)*ave_modwdh(ib,jb)+vb6_r*modwdhu3*c35
#                         hv7r = (v7_r-vb7_r)*ave_modwdh(ib,jb)+vb7_r*modwdhu3*c35
#                         hv8r = (v8_r-vb8_r)*ave_modwdh(ib,jb)+vb8_r*modwdhu33*c35
#                         hv9r = (v9_r-vb9_r)*ave_modwdh(ib,jb)+vb9_r*modwdhu33*c35
#                         hv10r = (v10_r-vb10_r)*ave_modwdh(ib,jb) &
#                         +vb10_r*modwdhu33*c35
#                         hv1rht1 = (v1_r-vb1_r)*ave_modwdht1(is,ib,jb) &
#                         +vb1_r*modwdhu3ht1*c35
#                         hv2rht3 = (v2_r-vb2_r)*ave_modwdht3(is,ib,jb) &
#                         +vb2_r*modwdhu3ht3*c35
#                         hv3rht1 = (v3_r-vb3_r)*ave_modwdht1(is,ib,jb) &
#                         +vb3_r*modwdhu33ht1*c35
#                         hv4rht3 = (v4_r-vb4_r)*ave_modwdht3(is,ib,jb) &
#                         +vb4_r*modwdhu33ht3*c35
#                         hv6rhu1 = (v6_r*modwdhu1)
#                         hv7rhu3 = (v7_r*modwdhu3)
#                         hv9rhu1 = (v9_r*modwdhu1)
#                         hv10rhu3= (v10_r*modwdhu3)
#                         hv1i = (v1_i-vb1_i)*ave_wdh(ib,jb)+vb1_i*wdhu3*c35
#                         hv2i = (v2_i-vb2_i)*ave_wdh(ib,jb)+vb2_i*wdhu3*c35
#                         hv3i = (v3_i-vb3_i)*ave_wdh(ib,jb)+vb3_i*wdhu33*c35
#                         hv4i = (v4_i-vb4_i)*ave_wdh(ib,jb)+vb4_i*wdhu33*c35
#                         hv5i = (v5_i-vb5_i)*ave_wdh(ib,jb)+vb5_i*wdhu3*c35
#                         hv6i = (v6_i-vb6_i)*ave_wdh(ib,jb)+vb6_i*wdhu3*c35
#                         hv7i = (v7_i-vb7_i)*ave_wdh(ib,jb)+vb7_i*wdhu3*c35
#                         hv8i = (v8_i-vb8_i)*ave_wdh(ib,jb)+vb8_i*wdhu33*c35
#                         hv9i = (v9_i-vb9_i)*ave_wdh(ib,jb)+vb9_i*wdhu33*c35
#                         hv10i = (v10_i-vb10_i)*ave_wdh(ib,jb)+vb10_i*wdhu33*c35
#                         hv1iht1 = (v1_i-vb1_i)*ave_wdht1(is,ib,jb) &
#                         +vb1_i*wdhu3ht1*c35
#                         hv2iht3 = (v2_i-vb2_i)*ave_wdht3(is,ib,jb) &
#                         +vb2_i*wdhu3ht3*c35
#                         hv3iht1 = (v3_i-vb3_i)*ave_wdht1(is,ib,jb) &
#                         +vb3_i*wdhu33ht1*c35
#                         hv4iht3 = (v4_i-vb4_i)*ave_wdht3(is,ib,jb) &
#                         +vb4_i*wdhu33ht3*c35
#                         hv6ihu1 = (v6_i*wdhu1)
#                         hv7ihu3 = (v7_i*wdhu3)
#                         hv9ihu1 = (v9_i*wdhu1)
#                         hv10ihu3= (v10_i*wdhu3)

#                         kpar_hnp0 = k_par0*ave_kparhnp0(is,ib,jb)
#                         kpar_hp1p0 = k_par0*ave_kparhp1p0(is,ib,jb)
#                         kpar_hp3p0 = k_par0*ave_kparhp3p0(is,ib,jb)
#                         if(use_bper_in)then
#                             kpar_hp1b0 = k_par0*ave_kparhp1b0(is,ib,jb)
#                             kpar_hr11b0 = k_par0*ave_kparhr11b0(is,ib,jb)
#                             kpar_hr13b0 = k_par0*ave_kparhr13b0(is,ib,jb)
#                             if(vpar_model_in.eq.0)then
#                                 kpar_hnbp = k_par0*ave_kparhnbp(is,ib,jb)
#                                 kpar_hp1bp = k_par0*ave_kparhp1bp(is,ib,jb)
#                                 kpar_hp3bp = k_par0*ave_kparhp3bp(is,ib,jb)
#                                 kpar_hr11bp = k_par0*ave_kparhr11bp(is,ib,jb)
#                                 kpar_hr13bp = k_par0*ave_kparhr13bp(is,ib,jb)
#                             endif
#                         endif
#                         kpar_hu1 = ave_kparhu1(is,ib,jb)
#                         kpar_hu3 = ave_kparhu3(is,ib,jb)
#                         kpar_hb1 = b1*ave_kpar_eff(is,ib,jb)
#                         kpar_hb3 = b3*ave_kpar_eff(is,ib,jb)
#                         kpar_hb33 = b33*ave_kpar_eff(is,ib,jb)
#                         kpar_hb1ht1 = b1*ave_kparht1(is,ib,jb)
#                         kpar_hb3ht3 = b3*ave_kparht3(is,ib,jb)
#                         kpar_hb33ht1 = b33*ave_kparht1(is,ib,jb)
#                         modkpar_hd1 = d1*ave_modkpar_eff(is,ib,jb)
#                         modkpar_hd3 = d3*ave_modkpar_eff(is,ib,jb)
#                         modkpar_hd33 = d33*ave_modkpar_eff(is,ib,jb)
#                         modkpar_hd1hu1 = d1*ave_modkparhu1(is,ib,jb)
#                         modkpar_hd3hu3 = d3*ave_modkparhu3(is,ib,jb)
#                         modkpar_hd33hu1 = d33*ave_modkparhu1(is,ib,jb)
#                         grad_hu1 = 0.0
#                         grad_hu3 = 0.0
#                         dhr13 = -b3*(ave_kparht3(is,ib,jb)-ave_kparht1(is,ib,jb)/3.0 &
#                                 -ave_kparhu3(is,ib,jb) + ave_kparhu1(is,ib,jb)/3.0)
            
#                         if(Linsker.eq.0.0)then
#                             gradhp1=0.0
#                             gradhr11=0.0
#                             gradhr13=0.0
#                             gradhp1p1=0.0
#                             gradhr11p1=0.0
#                             gradhr13p1=0.0
#                         else
#                             gradhp1 = Linsker*ave_gradhp1p0(is,ib,jb)
#                             gradhr11 = Linsker*ave_gradhr11p0(is,ib,jb)
#                             gradhr13 = Linsker*ave_gradhr13p0(is,ib,jb)
#                             gradhp1p1 = Linsker*ave_gradhp1p1(is,ib,jb)
#                             gradhr11p1 = Linsker*ave_gradhr11p1(is,ib,jb)
#                             gradhr13p1 = Linsker*ave_gradhr13p1(is,ib,jb) 
#                         endif

#                         if(nroot.gt.6)then
#                             gn = ave_gnp0(is,ib,jb)
#                             gp1 = ave_gp1p0(is,ib,jb)
#                             gp3 = ave_gp3p0(is,ib,jb)
#                             gr11 = ave_gr11p0(is,ib,jb)
#                             gr13 = ave_gr13p0(is,ib,jb)
#                             gr33 = ave_gr33p0(is,ib,jb)
#                             c_tor_par_gp1 = ave_c_tor_par_gp1p0(is,ib,jb)/Rmaj_input
#                             c_tor_par_gr11 = ave_c_tor_par_gr11p0(is,ib,jb)/Rmaj_input
#                             c_tor_par_gr13 = ave_c_tor_par_gr13p0(is,ib,jb)/Rmaj_input
#                             if(use_bper_in)then
#                                 gnb0 = ave_gnb0(is,ib,jb)
#                                 gp1b0 = ave_gp1b0(is,ib,jb)
#                                 gp3b0 = ave_gp3b0(is,ib,jb)
#                                 gr11b0 = ave_gr11b0(is,ib,jb)
#                                 gr13b0 = ave_gr13b0(is,ib,jb)
#                                 gr33b0 = ave_gr33b0(is,ib,jb)
#                                 gw113b0 = ave_gw113b0(is,ib,jb)
#                                 gw133b0 = ave_gw133b0(is,ib,jb)
#                                 gw333b0 = ave_gw333b0(is,ib,jb)
#                                 if(vpar_model_in.eq.0)then
#                                     gnbp = ave_gnbp(is,ib,jb)
#                                     gp1bp = ave_gp1bp(is,ib,jb)
#                                     gp3bp = ave_gp3bp(is,ib,jb)
#                                     gr11bp = ave_gr11bp(is,ib,jb)
#                                     gr13bp = ave_gr13bp(is,ib,jb)
#                                     gr33bp = ave_gr33bp(is,ib,jb)
#                                     gw113bp = ave_gw113bp(is,ib,jb)
#                                     gw133bp = ave_gw133bp(is,ib,jb)
#                                     gw333bp = ave_gw333bp(is,ib,jb)
#                                 endif
#                             endif
#                             if(use_bpar_in)then
#                                 g10n = 1.5*(gnb0-gp3b0)
#                                 g10p1 = 2.5*gp1b0 - 1.5*gr13b0
#                                 g10p3 = 2.5*gp3b0 - 1.5*gr33b0
#                                 g10r13 = 3.5*gr13b0 - 1.5*gw133b0
#                                 g10r33 = 3.5*gr33b0 - 1.5*gw333b0
#                             endif
                            
#                             gu1 = ave_gu1(is,ib,jb)
#                             gu3 = ave_gu3(is,ib,jb)
#                             gt1 = ave_gt1(is,ib,jb)
#                             gt3 = ave_gt3(is,ib,jb)

#                             wdgp1p0 = ave_wdgp1p0(is,ib,jb)
#                             wdgr11p0 = ave_wdgr11p0(is,ib,jb)
#                             wdgr13p0 = ave_wdgr13p0(is,ib,jb)
#                             if(use_bper_in)then
#                                 wdgp1b0 = ave_wdgp1b0(is,ib,jb)
#                                 wdgr11b0 = ave_wdgr11b0(is,ib,jb)
#                                 wdgr13b0 = ave_wdgr13b0(is,ib,jb)
#                                 if(vpar_model_in.eq.0)then
#                                     wdgp1bp = ave_wdgp1bp(is,ib,jb)
#                                     wdgr11bp = ave_wdgr11bp(is,ib,jb)
#                                     wdgr13bp = ave_wdgr13bp(is,ib,jb)
#                                 endif
#                             endif

#                             wdgu1 = ave_wdgu1(is,ib,jb)
#                             wdgu3 = ave_wdgu3(is,ib,jb)
#                             wdgu33 = ave_wdgu33(is,ib,jb)
#                             wdgu3gt1 = ave_wdgu3gt1(is,ib,jb)
#                             wdgu3gt3 = ave_wdgu3gt3(is,ib,jb)
#                             wdgu33gt1 = ave_wdgu33gt1(is,ib,jb)
#                             wdgu33gt3 = ave_wdgu33gt3(is,ib,jb)
#                             modwdgu1 = ave_modwdgu1(is,ib,jb)
#                             modwdgu3 = ave_modwdgu3(is,ib,jb)
#                             modwdgu33 = ave_modwdgu33(is,ib,jb)
#                             modwdgu3gt1 = ave_modwdgu3gt1(is,ib,jb)
#                             modwdgu3gt3 = ave_modwdgu3gt3(is,ib,jb)
#                             modwdgu33gt1 = ave_modwdgu33gt1(is,ib,jb)
#                             modwdgu33gt3 = ave_modwdgu33gt3(is,ib,jb)
#                             gu1r = (u1_r-ub1_r)*ave_modwdg(ib,jb)+ub1_r*modwdgu3*c35
#                             gu2r = (u2_r-ub2_r)*ave_modwdg(ib,jb)+ub2_r*modwdgu3*c35
#                             gu3r = (u3_r-ub3_r)*ave_modwdg(ib,jb)+ub3_r*modwdgu33*c35
#                             gu4r = (u4_r-ub4_r)*ave_modwdg(ib,jb)+ub4_r*modwdgu33*c35
#                             gu5r = (u5_r-ub5_r)*ave_modwdg(ib,jb)+ub5_r*modwdgu3*c35
#                             gu6r = (u6_r-ub6_r)*ave_modwdg(ib,jb)+ub6_r*modwdgu3*c35
#                             gu7r = (u7_r-ub7_r)*ave_modwdg(ib,jb)+ub7_r*modwdgu3*c35
#                             gu8r = (u8_r-ub8_r)*ave_modwdg(ib,jb)+ub8_r*modwdgu33*c35
#                             gu9r = (u9_r-ub9_r)*ave_modwdg(ib,jb)+ub9_r*modwdgu33*c35
#                             gu10r = (u10_r-ub10_r)*ave_modwdg(ib,jb) &
#                                     +ub10_r*modwdgu33*c35
#                             gu1rgt1 = (u1_r-ub1_r)*ave_modwdgt1(is,ib,jb) &
#                                     +ub1_r*modwdgu3gt1*c35
#                             gu2rgt3 = (u2_r-ub2_r)*ave_modwdgt3(is,ib,jb) &
#                                     +ub2_r*modwdgu3gt3*c35
#                             gu3rgt1 = (u3_r-ub3_r)*ave_modwdgt1(is,ib,jb) &
#                                     +ub3_r*modwdgu33gt1*c35
#                             gu4rgt3 = (u4_r-ub4_r)*ave_modwdgt3(is,ib,jb) &
#                                     +ub4_r*modwdgu33gt3*c35
#                             gu6rgu1 = (u6_r*modwdgu1)
#                             gu7rgu3 = (u7_r*modwdgu3)
#                             gu9rgu1 = (u9_r*modwdgu1)
#                             gu10rgu3= (u10_r*modwdgu3)
#                             gu1i = (u1_i-ub1_i)*ave_wdg(ib,jb)+ub1_i*wdgu3*c35
#                             gu2i = (u2_i-ub2_i)*ave_wdg(ib,jb)+ub2_i*wdgu3*c35
#                             gu3i = (u3_i-ub3_i)*ave_wdg(ib,jb)+ub3_i*wdgu33*c35
#                             gu4i = (u4_i-ub4_i)*ave_wdg(ib,jb)+ub4_i*wdgu33*c35
#                             gu5i = (u5_i-ub5_i)*ave_wdg(ib,jb)+ub5_i*wdgu3*c35
#                             gu6i = (u6_i-ub6_i)*ave_wdg(ib,jb)+ub6_i*wdgu3*c35
#                             gu7i = (u7_i-ub7_i)*ave_wdg(ib,jb)+ub7_i*wdgu3*c35
#                             gu8i = (u8_i-ub8_i)*ave_wdg(ib,jb)+ub8_i*wdgu33*c35
#                             gu9i = (u9_i-ub9_i)*ave_wdg(ib,jb)+ub9_i*wdgu33*c35
#                             gu10i = (u10_i-ub10_i)*ave_wdg(ib,jb) &
#                                     +ub10_i*wdgu33*c35
#                             gu1igt1 = (u1_i-ub1_i)*ave_wdgt1(is,ib,jb) &
#                                     +ub1_i*wdgu3gt1*c35
#                             gu2igt3 = (u2_i-ub2_i)*ave_wdgt3(is,ib,jb) &
#                                     +ub2_i*wdgu3gt3*c35
#                             gu3igt1 = (u3_i-ub3_i)*ave_wdgt1(is,ib,jb) &
#                                     +ub3_i*wdgu33gt1*c35
#                             gu4igt3 = (u4_i-ub4_i)*ave_wdgt3(is,ib,jb) &
#                                     +ub4_i*wdgu33gt3*c35
#                             gu6igu1 = (u6_i*wdgu1)
#                             gu7igu3 = (u7_i*wdgu3)
#                             gu9igu1 = (u9_i*wdgu1)
#                             gu10igu3= (u10_i*wdgu3)

#                             kpar_gnp0 = k_par0*ave_kpargnp0(is,ib,jb)
#                             kpar_gp1p0 = k_par0*ave_kpargp1p0(is,ib,jb)
#                             kpar_gp3p0 = k_par0*ave_kpargp3p0(is,ib,jb)
#                             if(use_bper_in)then
#                                 kpar_gp1b0 = k_par0*ave_kpargp1p0(is,ib,jb)
#                                 kpar_gr11b0 = k_par0*ave_kpargr11b0(is,ib,jb)
#                                 kpar_gr13b0 = k_par0*ave_kpargr13b0(is,ib,jb)
#                                 if(vpar_model_in.eq.0)then
#                                     kpar_gnbp = k_par0*ave_kpargnbp(is,ib,jb)
#                                     kpar_gp1bp = k_par0*ave_kpargp1bp(is,ib,jb)
#                                     kpar_gp3bp = k_par0*ave_kpargp3bp(is,ib,jb)
#                                     kpar_gr11bp = k_par0*ave_kpargr11bp(is,ib,jb)
#                                     kpar_gr13bp = k_par0*ave_kpargr13bp(is,ib,jb)
#                                 endif
#                             endif
#                             kpar_gu1 = ave_kpargu1(is,ib,jb)
#                             kpar_gu3 = ave_kpargu3(is,ib,jb)
#                             kpar_gb1 = ft2*b1*ave_kpar_eff(is,ib,jb)
#                             kpar_gb3 = ft2*b3*ave_kpar_eff(is,ib,jb)
#                             kpar_gb33 = b33*ave_kpar_eff(is,ib,jb)
#                             kpar_gb1gt1 = ft2*b1*ave_kpargt1(is,ib,jb)
#                             kpar_gb3gt3 = ft2*b3*ave_kpargt3(is,ib,jb)
#                             kpar_gb33gt1 = b33*ave_kpargt1(is,ib,jb)
#                             modkpar_gd1 = ft*d1*ave_modkpar_eff(is,ib,jb)
#                             modkpar_gd3 = ft*d3*ave_modkpar_eff(is,ib,jb)
#                             modkpar_gd33 = (d33/ft)*ave_modkpar_eff(is,ib,jb)
#                             modkpar_gd1gu1 = ft*d1*ave_modkpargu1(is,ib,jb)
#                             modkpar_gd3gu3 = ft*d3*ave_modkpargu3(is,ib,jb)
#                             modkpar_gd33gu1 = (d33/ft)*ave_modkpargu1(is,ib,jb)

#                             grad_gu1 = 0.0
#                             grad_gu3 = 0.0
#                             dgr13 = -b3*(ft2*ave_kpargt3(is,ib,jb)-ave_kpargt1(is,ib,jb)/3.0 &
#                                 -ft2*ave_kpargu3(is,ib,jb) + ave_kpargu1(is,ib,jb)/3.0)

#                             if(nbasis.eq.1.or.Linsker.eq.0.0)then
#                                 gradgp1=0.0
#                                 gradgr11=0.0
#                                 gradgr13=0.0
#                                 gradgp1p1=0.0
#                                 gradgr11p1=0.0
#                                 gradgr13p1=0.0
#                             else
#                                 gradgp1 = Linsker*ave_gradgp1p0(is,ib,jb)
#                                 gradgr11 = Linsker*ave_gradgr11p0(is,ib,jb)
#                                 gradgr13 = Linsker*ave_gradgr13p0(is,ib,jb)
#                                 gradgp1p1 = Linsker*ave_gradgp1p1(is,ib,jb)
#                                 gradgr11p1 = Linsker*ave_gradgr11p1(is,ib,jb)
#                                 gradgr13p1 = Linsker*ave_gradgr13p1(is,ib,jb) 
#                             endif    
#                         endif  # nroot>6

#                         w_d1 = -ghat_in*w_d0
#                         k_par1 = k_par0
#                         if(js.ne.is)then
#                             w_d1 = 0.0
#                             k_par1 = 0.0
#                         endif
#                         modw_d1 = ABS(w_d1)
#                         modk_par1 = ABS(k_par1)
#                         modk_par0 = ABS(k_par0)
#                         w_dh= w_d1*ave_wdh(ib,jb)
#                         w_dg= w_d1*ave_wdg(ib,jb)
#                         k_par = k_par1*ave_kpar_eff(is,ib,jb)
#                         k_par_psi = k_par0*ave_kpar_eff(is,ib,jb)
#                         gradB1=k_par1*gradB_factor_in
#                         if(nbasis.eq.1.or.gradB1.eq.0.0)then
#                             gradB=0.0
#                             gradBhp1=0.0
#                             gradBhp3=0.0
#                             gradBhr11=0.0
#                             gradBhr13=0.0
#                             gradBhr33=0.0
#                             gradBhu1=0.0
#                             gradBhu3=0.0
#                             gradBhu33=0.0
#                             if(nroot.gt.6)then
#                                 gradBgp1=0.0
#                                 gradBgp3=0.0
#                                 gradBgr11=0.0
#                                 gradBgr13=0.0
#                                 gradBgr33=0.0
#                                 gradBgu1=0.0
#                                 gradBgu3=0.0
#                                 gradBgu33=0.0
#                             endif
#                         else
#                         gradB = gradB1*ave_gradB(ib,jb)
#                         gradBhp1=gradB1*ave_gradBhp1(is,ib,jb)
#                         gradBhp3=gradB1*ave_gradBhp3(is,ib,jb)
#                         gradBhr11=gradB1*ave_gradBhr11(is,ib,jb)
#                         gradBhr13=gradB1*ave_gradBhr13(is,ib,jb)
#                         gradBhr33=gradB1*ave_gradBhr33(is,ib,jb)
#                         gradBhu1=gradB1*ave_gradBhu1(is,ib,jb)
#                         gradBhu3=gradB1*ave_gradBhu3(is,ib,jb)
#                         gradBhu33=gradB1*ave_gradBhu33(is,ib,jb)
#                         if(nroot.gt.6)then
#                             gradBgp1=gradB1*ave_gradBgp1(is,ib,jb)
#                             gradBgp3=gradB1*ave_gradBgp3(is,ib,jb)
#                             gradBgr11=gradB1*ave_gradBgr11(is,ib,jb)
#                             gradBgr13=gradB1*ave_gradBgr13(is,ib,jb)
#                             gradBgr33=gradB1*ave_gradBgr33(is,ib,jb)
#                             gradBgu1=gradB1*ave_gradBgu1(is,ib,jb)
#                             gradBgu3=gradB1*ave_gradBgu3(is,ib,jb)
#                             gradBgu33=gradB1*ave_gradBgu33(is,ib,jb)
#                         endif
#                     endif

#                     M_i = zs(is)*vs(is)/taus(is)
#                     J_j = zs(js)*as(js)*vs(js)
#                     E_i = zs(is)/taus(is)
#                     N_j = zs(js)*as(js)
# #
# # matrix in order n, u, p1, p3, q1,q3
# #
# # n_u equ #1
# #
#                     ia0 = (is-ns0)*nroot*nbasis
#                     ja0 = (js-ns0)*nroot*nbasis
      
#                     ia = ib + ia0
# #
# #  untrapped terms
# #
#                     ja = jb + ja0

#                     phi_A = N_j*xi*w_s*(rlns(is)*hn + rlts(is)*1.5*(hp3-hn)) 
#                     if(vpar_model_in.eq.0)then
#                         phi_A = phi_A + N_j*E_i*kpar_hnp0*vpar(is)
#                     endif
#                     phi_B = -hn*E_i*N_j
#                     sig_A = 0.0
#                     sig_B = 0.0
#                     psi_A = 0.0
#                     psi_B = 0.0
#                     phi_AU = 0.0
#                     phi_BU = 0.0
#                     psi_AN = 0.0
#                     psi_BN = 0.0
#                     if(use_bpar_in)then
#                         sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
#                         (xi*w_s*(rlns(is)*h10n + rlts(is)*1.5*(h10p3-h10n)))
#                         sig_B = betae_sig*h10n*as(js)*taus(js)*zs(is)*zs(is) &
#                         /(taus(is)*mass(is))
#                         sig_A = sig_A - damp_sig*sig_B
#                     endif
#                     if(use_bper_in)then
#                         psi_A = -betae_psi*J_j*xi*w_s*vpar_shear(is)*hp1b0
#                         if(vpar_model_in.eq.0)then
#                             psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar(is)*wdhp1b0
#                             psi_B = betae_psi*M_i*J_j*vpar(is)*hp1b0/vs(is)
#                             phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*hnbp + rlts(is)*1.5*(hp3bp-hnbp)) &
#                             + E_i*kpar_hnbp*vpar(is))
#                             phi_BU = -betae_psi*U0*E_i*J_j*vpar(is)*hp1bp
#                             psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar(is)*wdhp1bp+w_s*vpar_shear(is)*hp1bp)
#                             psi_BN = -betae_psi*U0*M_i*N_j*vpar(is)*hp1bp/vs(is)
#                         endif
#                     endif

#                     amat(ia,ja) = phi_A + psi_AN  
#                     bmat(ia,ja) = d_ab + phi_B + psi_BN

#                     ja = nbasis+jb + ja0
#                     amat(ia,ja) = psi_A + phi_AU -k_par*vs(is) + am*gradB*vs(is)
#                     bmat(ia,ja) = psi_B + phi_BU

#                     ja = 2*nbasis+jb + ja0
#                     amat(ia,ja) = -0.5*sig_A -0.5*xi*w_dh*taus(is)/zs(is)
#                     bmat(ia,ja) = -0.5*sig_B

#                     ja = 3*nbasis+jb + ja0
#                     amat(ia,ja) = 1.5*sig_A -1.5*xi*w_dh*taus(is)/zs(is)
#                     bmat(ia,ja) = 1.5*sig_B

#                     ja = 4*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0

#                     ja = 5*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0

#                     if(nroot.gt.6)then
# #
# #  n_u ghost terms
# #
#                         ja = 6*nbasis+jb + ja0
#                         amat(ia,ja) = -1.0*(phi_A + psi_AN)
#                         bmat(ia,ja) = -1.0*(phi_B + psi_BN)

#                         ja = 7*nbasis+jb + ja0
#                         amat(ia,ja) = -1.0*(psi_A + phi_AU)
#                         bmat(ia,ja) = -1.0*(psi_B + phi_BU)

#                         ja = 8*nbasis+jb + ja0
#                         amat(ia,ja) = -0.5*(-1.0*sig_A)
#                         bmat(ia,ja) = -0.5*(-1.0*sig_B)

#                         ja = 9*nbasis+jb + ja0
#                         amat(ia,ja) = 1.5*(-1.0*sig_A)
#                         bmat(ia,ja) = 1.5*(-1.0*sig_B)

#                         ja = 10*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0
#                         bmat(ia,ja) = 0.0

#                         ja = 11*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0
#                         bmat(ia,ja) = 0.0
# #
# #  n_u trapped particle terms
# #
#                         ja = 12*nbasis+jb + ja0
#                         amat(ia,ja) = phi_A + psi_AN
#                         bmat(ia,ja) = phi_B + psi_BN

#                         ja = 13*nbasis+jb + ja0
#                         amat(ia,ja) = -0.5*sig_A
#                         bmat(ia,ja) = -0.5*sig_B

#                         ja = 14*nbasis+jb + ja0
#                         amat(ia,ja) = 1.5*sig_A
#                         bmat(ia,ja) = 1.5*sig_B

#                     endif #  nroot.gt.6
# #
# # u_par_u equ #2
# #
#                     ia = nbasis+ib + ia0
# #
#                     phi_A = N_j*xi*w_s*vpar_shear(is)*hp1/vs(is)
#                     phi_B = 0.0 
#                     if(vpar_model_in.eq.0)then
#                         phi_A = phi_A  +N_j*xi*w_cd*wdhp1p0*vpar(is)/vs(is)  &
#                             + d_1*(nuei_u_u_1+nuei_u_q3_1*5.0/3.0)*hp1*E_i*N_j*vpar(is)/vs(is)
#                         phi_B = -E_i*N_j*hp1*vpar(is)/vs(is)
#                     endif
#                     sig_A = 0.0
#                     sig_B = 0.0
#                     psi_A = 0.0
#                     psi_B = 0.0
#                     phi_AU = 0.0
#                     phi_BU = 0.0
#                     psi_AN = 0.0
#                     psi_BN = 0.0
#                     if(use_bper_in)then
#                         psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*hp1b0+1.5*rlts(is)*(hr13b0-hp1b0))  
#                         psi_B = betae_psi*M_i*J_j*hp1b0
#                         psi_A = psi_A - damp_psi*psi_B
#                         if(vpar_model_in.eq.0)then
#                         psi_A = psi_A  -betae_psi*J_j*M_i*kpar_hp1b0*vpar(is)
#                         phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear(is)*hp1bp +w_cd*wdhp1bp*vpar(is))/vs(is) &
#                             + d_1*(nuei_u_u_1+nuei_u_q3_1*5.0/3.0)*hp1*E_i*betae_psi*U0*J_j*vpar(is)/vs(is)
#                         phi_BU = -betae_psi*U0*E_i*J_j*hp1bp*vpar(is)/vs(is)
#                         psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*hp1bp+1.5*rlts(is)*(hr13bp-hp1bp))  &
#                             + M_i*kpar_hp1bp*vpar(is))
#                         psi_BN = -betae_psi*U0*M_i*N_j*hp1bp
#                         endif
#                     endif
# #
# # u_par_u  untrapped terms
# #
# #
#                     ja = jb + ja0
#                     amat(ia,ja) = phi_A + psi_AN
#                     bmat(ia,ja) = phi_B + psi_BN

#                     ja = nbasis+jb + ja0
#                     amat(ia,ja) =  psi_A +phi_AU  &
#                     -d_ee*nuei_u_u_1  &
#                     +xnuei*(d_ab*xnu_u_u_1 - d_ij*hu3*xnu_u_q3_1)

#                     bmat(ia,ja) = d_ab + psi_B +phi_BU

#                     ja = 2*nbasis+jb + ja0
#                     amat(ia,ja) =  -(k_par - k_par1*gradhp1p1)*vs(is) &
#                     + (am+bm*0.5)*gradB*vs(is) 
#                     bmat(ia,ja) = 0.0

#                     ja = 3*nbasis+jb + ja0
#                     amat(ia,ja) =  - bm*1.5*gradB*vs(is)
#                     bmat(ia,ja) = 0.0

#                     ja = 4*nbasis+jb + ja0
#                     amat(ia,ja) = -0.5*xi*w_dh*taus(is)/zs(is)
#                     bmat(ia,ja) = 0.0

#                     ja = 5*nbasis+jb + ja0
#                     amat(ia,ja) = -1.5*xi*w_dh*taus(is)/zs(is) &
#                     -d_ee*nuei_u_q3_1 &
#                     +xnuei*d_ab*xnu_u_q3_1
#                     bmat(ia,ja) = 0.0

#                     if(nroot.gt.6)then
# #
# #  u_par_u ghost terms
# #
#                         ja = 6*nbasis+jb + ja0
#                         amat(ia,ja) = -1.0*(phi_A + psi_AN)
#                         bmat(ia,ja) = -1.0*(phi_B + psi_BN)
    
#                         ja = 7*nbasis+jb + ja0
#                         amat(ia,ja) = -1.0*(psi_A + phi_AU)
#                         bmat(ia,ja) = -1.0*(psi_B + phi_BU)
    
#                         ja = 8*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0
#                         bmat(ia,ja) = 0.0
    
#                         ja = 9*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0 
#                         bmat(ia,ja) = 0.0
    
#                         ja = 10*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0
#                         bmat(ia,ja) = 0.0
    
#                         ja = 11*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0
#                         bmat(ia,ja) = 0.0
#                     #
#                     # u_par_u trapped particle terms
#                     #
#                         ja = 12*nbasis+jb + ja0
#                         amat(ia,ja) = phi_A + psi_AN
#                         bmat(ia,ja) = phi_B + psi_BN

#                         ja = 13*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0
#                         bmat(ia,ja) = 0.0

#                         ja = 14*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0
#                         bmat(ia,ja) = 0.0
# #
#                     endif  #  nroot.gt.6
# #
# # p_par_u equ #3
# #
#                     ia = 2*nbasis+ib + ia0

#                     phi_A = N_j*xi*w_s*(rlns(is)*hp1 + rlts(is)*1.5*(hr13-hp1))  
#                     if(vpar_model_in.eq.0)then
#                         phi_A = phi_A +N_j*E_i*kpar_hp1p0*vpar(is)
#                     endif
#                     phi_B = -hp1*E_i*N_j
#                     sig_A = 0.0
#                     sig_B = 0.0
#                     psi_A = 0.0
#                     psi_B = 0.0
#                     phi_AU = 0.0
#                     phi_BU = 0.0
#                     psi_AN = 0.0
#                     psi_BN = 0.0
#                     if(use_bpar_in)then
#                         sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
#                         (xi*w_s*(rlns(is)*h10p1 + rlts(is)*1.5*(h10r13-h10p1)))
#                         sig_B = betae_sig*h10p1*as(js)*taus(js)*zs(is)*zs(is) &
#                         /(taus(is)*mass(is))
#                         sig_A = sig_A - damp_sig*sig_B
#                     endif
#                     if(use_bper_in)then
#                         psi_A = -betae_psi*J_j*xi*w_s*vpar_shear(is)*hr11b0
#                         if(vpar_model_in.eq.0)then
#                             psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar(is)*wdhr11b0
#                             psi_B = betae_psi*M_i*J_j*vpar(is)*hr11b0/vs(is)
#                             phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*hp1bp + rlts(is)*1.5*(hr13bp-hp1bp)) &
#                             +E_i*kpar_hp1bp*vpar(is))
#                             phi_BU = -betae_psi*U0*hp1bp*E_i*J_j
#                             psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar(is)*wdhr11bp+w_s*vpar_shear(is)*hr11bp)
#                             psi_BN = -betae_psi*U0*M_i*N_j*vpar(is)*hr11bp/vs(is)
#                         endif
#                     endif
# #
# #  p_par_u  untrapped terms
# #
#                     ja = jb + ja0
#                     amat(ia,ja) = phi_A + psi_AN  &
#                     +2.0*taus(is)*(modw_d1*hv1rht1/ABS(zs(is)) +w_d1*xi*hv1iht1/zs(is)) &
#                     +2.0*taus(is)*(modw_d1*hv2rht3/ABS(zs(is)) +w_d1*xi*hv2iht3/zs(is)) &
#                     -xnuei*d_ij*xnu_p1_1*(ht1-ht3)
#                     bmat(ia,ja) = phi_B + psi_BN
#                 #
#                     ja = nbasis+jb + ja0
#                     amat(ia,ja) = k_par1*grad_hu1*vs(is) + psi_A + phi_AU
#                     bmat(ia,ja) = psi_B + phi_BU
#                 #
#                     ja = 2*nbasis+jb + ja0
#                     amat(ia,ja) =  -0.5*sig_A   & 
#                     -xi*w_d1*(taus(is)/zs(is))*(0.5*wdhu1+1.5*wdhu3) &
#                     -2.0*taus(is)*(modw_d1*hv1r/ABS(zs(is)) +w_d1*xi*hv1i/zs(is)) &
#                     -d_ee*nuei_p1_p1_1  &
#                     +xnuei*d_ab*xnu_p1_1
#                     bmat(ia,ja) = d_ab -0.5*sig_B
#                 #
#                     ja = 3*nbasis+jb + ja0
#                     amat(ia,ja) = 1.5*sig_A -2.0*taus(is)* &
#                         (modw_d1*hv2r/ABS(zs(is)) +w_d1*xi*hv2i/zs(is)) &
#                     -d_ee*nuei_p1_p3_1  &
#                     -xnuei*d_ab*xnu_p1_1
#                     bmat(ia,ja) = 1.5*sig_B
#                 #
#                     ja = 4*nbasis+jb + ja0
#                     amat(ia,ja) = -k_par*vs(is) + (am +bm)*gradB*vs(is)
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 5*nbasis+jb + ja0
#                     amat(ia,ja) = -bm*3.0*gradB*vs(is)
#                     bmat(ia,ja) = 0.0
#                 #
#                     if(nroot.gt.6)then
# #
# #   p_par_u ghost terms
# #
#                         ja = 6*nbasis+jb + ja0
#                         amat(ia,ja) = -1.0*(phi_A + psi_AN)
#                         bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#                     #
#                         ja = 7*nbasis+jb + ja0
#                         amat(ia,ja) = -1.0*(psi_A + phi_AU)
#                         bmat(ia,ja) = -1.0*(psi_B + phi_BU)
#                     #
#                         ja = 8*nbasis+jb + ja0
#                         amat(ia,ja) = -0.5*(-1.0*sig_A)
#                         bmat(ia,ja) = -0.5*(-1.0*sig_B)
#                     #
#                         ja = 9*nbasis+jb + ja0
#                         amat(ia,ja) = 1.5*(-1.0*sig_A)
#                         bmat(ia,ja) = 1.5*(-1.0*sig_B)
#                     #
#                         ja = 10*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0
#                         bmat(ia,ja) = 0.0
#                     #
#                         ja = 11*nbasis+jb + ja0
#                         amat(ia,ja) = 0.0
#                         bmat(ia,ja) = 0.0
# #
# #  p_par_u trapped particle terms
# #
#                     ja = 12*nbasis+jb + ja0
#                     amat(ia,ja) = phi_A + psi_AN
#                     bmat(ia,ja) = phi_B + psi_BN
#                 #
#                     ja = 13*nbasis+jb + ja0
#                     amat(ia,ja) = -0.5*sig_A
#                     bmat(ia,ja) = -0.5*sig_B
#                 #
#                     ja = 14*nbasis+jb + ja0
#                     amat(ia,ja) = 1.5*sig_A
#                     bmat(ia,ja) = 1.5*sig_B
#                 #
#                 endif   # nroot >6
# #
# # p_tot_u equ #4
# #
#                 ia = 3*nbasis+ib + ia0
# #
#                 phi_A = N_j*xi*w_s*(rlns(is)*hp3+rlts(is)*1.5*(hr33-hp3))  
#                 if(vpar_model_in.eq.0)then
#                     phi_A = phi_A + N_j*E_i*kpar_hp3p0*vpar(is)
#                 endif
#                 phi_B = -hp3*E_i*N_j
#                 sig_A = 0.0
#                 sig_B = 0.0
#                 psi_A = 0.0
#                 psi_B = 0.0
#                 phi_AU = 0.0
#                 phi_BU = 0.0
#                 psi_AN = 0.0
#                 psi_BN = 0.0
#                 if(use_bpar_in)then
#                     sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
#                     (xi*w_s*(rlns(is)*h10p3 + rlts(is)*1.5*(h10r33-h10p3)))
#                     sig_B = betae_sig*h10p3*as(js)*taus(js)*zs(is)*zs(is) &
#                     /(taus(is)*mass(is))
#                     sig_A = sig_A - damp_sig*sig_B
#                 endif
#                 if(use_bper_in)then
#                     psi_A = -betae_psi*J_j*xi*w_s*vpar_shear(is)*hr13b0
#                     if(vpar_model_in.eq.0)then
#                         psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar(is)*wdhr13b0
#                         psi_B = betae_psi*M_i*J_j*vpar(is)*hr13b0/vs(is)
#                         phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*hp3bp+rlts(is)*1.5*(hr33bp-hp3bp))    &
#                         + E_i*kpar_hp3bp*vpar(is))
#                         phi_BU = -betae_psi*U0*hp3bp*E_i*J_j
#                         psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar(is)*wdhr13bp +w_s*vpar_shear(is)*hr13bp)
#                         psi_BN = -betae_psi*U0*M_i*N_j*vpar(is)*hr13bp/vs(is)
#                     endif
#                 endif
# #
# #   p_tot_u untrapped terms
# #
#                 ja = jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN  &
#                 +2.0*taus(is)*(modw_d1*hv3rht1/ABS(zs(is)) +w_d1*xi*hv3iht1/zs(is)) &
#                 +2.0*taus(is)*(modw_d1*hv4rht3/ABS(zs(is)) +w_d1*xi*hv4iht3/zs(is))
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = nbasis+jb + ja0
#                 amat(ia,ja) = k_par1*grad_hu3*vs(is) + psi_A + phi_AU
#                 bmat(ia,ja) = psi_B + phi_BU
#             #
#                 ja = 2*nbasis+jb + ja0
#                 amat(ia,ja) = -0.5*sig_A  &
#                 -xi*w_d1*(taus(is)/zs(is))*0.5*wdhu3 &
#                 -2.0*taus(is)*(modw_d1*hv3r/ABS(zs(is)) +w_d1*xi*hv3i/zs(is))
#                 bmat(ia,ja) = -0.5*sig_B
#             #
#                 ja = 3*nbasis+jb + ja0
#                 amat(ia,ja) =  1.5*sig_A  &
#                 -xi*w_d1*(taus(is)/zs(is))*1.5*wdhu33 &
#                 -2.0*taus(is)*(modw_d1*hv4r/ABS(zs(is)) +w_d1*xi*hv4i/zs(is))
#                 bmat(ia,ja) = d_ab + 1.5*sig_B
#             #
#                 ja = 4*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 5*nbasis+jb + ja0
#                 amat(ia,ja) = -k_par*vs(is) + am*gradB*vs(is)
#                 bmat(ia,ja) = 0.0
#             #
#                 if(nroot.gt.6)then
#             #
#             #   p_tot_u ghost terms
#             #
#                 ja = 6*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(phi_A + psi_AN)
#                 bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#             #
#                 ja = 7*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(psi_A + phi_AU)
#                 bmat(ia,ja) = -1.0*(psi_B + phi_BU)
#             #
#                 ja = 8*nbasis+jb + ja0
#                 amat(ia,ja) = -0.5*(-1.0*sig_A)
#                 bmat(ia,ja) = -0.5*(-1.0*sig_B)
#             #
#                 ja = 9*nbasis+jb + ja0
#                 amat(ia,ja) = 1.5*(-1.0*sig_A)
#                 bmat(ia,ja) = 1.5*(-1.0*sig_B)
#             #
#                 ja = 10*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 11*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#             #   p_tot_u trapped terms
#             #
#                 ja = 12*nbasis+jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = 13*nbasis+jb + ja0
#                 amat(ia,ja) = -0.5*sig_A
#                 bmat(ia,ja) = -0.5*sig_B
#             #
#                 ja = 14*nbasis+jb + ja0
#                 amat(ia,ja) = 1.5*sig_A
#                 bmat(ia,ja) = 1.5*sig_B
# #
#             endif  # nroot>6
# #
# # q_par_u equ #5
# #
#             ia = 4*nbasis+ib + ia0
#         #
#             phi_A = N_j*xi*w_s*hr11*vpar_shear(is)/vs(is)
#             phi_B = 0.0
#             if(vpar_model_in.eq.0)then
#                 phi_A = phi_A +N_j*xi*w_cd*wdhr11p0*vpar(is)/vs(is) &
#                 + d_1*(nuei_q1_u_1 + (5.0/3.0)*nuei_q1_q3_1+3.0*nuei_q1_q1_1)*hr11*E_i*N_j*vpar(is)/vs(is)
#                 phi_B = -E_i*N_j*hr11*vpar(is)/vs(is)
#             endif
#             sig_A = 0.0
#             sig_B = 0.0
#             psi_A = 0.0
#             psi_B = 0.0
#             phi_AU = 0.0
#             phi_BU = 0.0
#             psi_AN = 0.0
#             psi_BN = 0.0
#             if(use_bper_in)then
#                 psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*hr11b0+1.5*rlts(is)*(hw113b0-hr11b0)) 
#                 psi_B =betae_psi*M_i*J_j*hr11b0
#                 psi_A = psi_A - damp_psi*psi_B
#                 if(vpar_model_in.eq.0)then
#                     psi_A = psi_A  -betae_psi*J_j*M_i*kpar_hr11b0*vpar(is)
#                     phi_AU = betae_psi*U0*J_j*xi*(w_s*hr11bp*vpar_shear(is) +w_cd*wdhr11bp*vpar(is))/vs(is) &
#                     + d_1*(nuei_q1_u_1 + (5.0/3.0)*nuei_q1_q3_1+3.0*nuei_q1_q1_1)*hr11bp*E_i*betae_psi*U0*J_j*vpar(is)/vs(is)
#                     phi_BU = -betae_psi*U0*E_i*J_j*hr11bp*vpar(is)/vs(is)
#                     psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*hr11bp+1.5*rlts(is)*(hw113bp-hr11bp)) &
#                     +M_i*kpar_hr11bp*vpar(is))
#                     psi_BN =-betae_psi*U0*M_i*N_j*hr11bp
#                 endif
#             endif
# #
# #  q_par_u untrapped terms
# #
#             ja =          jb + ja0
#             amat(ia,ja) = k_par1*kpar_hb1ht1*vs(is) + phi_A + psi_AN
#             bmat(ia,ja) = phi_B + psi_BN
#         #
#             ja =   nbasis+jb + ja0
#             amat(ia,ja) = psi_A + phi_AU +modk_par1*modkpar_hd1hu1*vs(is) &
#             - taus(is)*(modw_d1*hv5r/ABS(zs(is)) + w_d1*xi*hv5i/zs(is)) &
#             -d_ee*nuei_q1_u_1  &
#             +xnuei*(d_ab*xnu_q1_u_1 - d_ij*hu1*xnu_q1_q1_1 &
#             - d_ij*hu3*xnu_q1_q3_1)
#             bmat(ia,ja) = psi_B + phi_BU
#         #
#             ja = 2*nbasis+jb + ja0
#             amat(ia,ja) = -k_par1*(kpar_hu1 -grad_hu1 - gradhr11p1)*vs(is) &
#             - k_par1*kpar_hb1*vs(is) + (am+bm*1.5)*gradBhu1*vs(is)      &
#             -d_11*k_par*vs(is)*c06
#             bmat(ia,ja) = 0.0
#         #
#             ja = 3*nbasis+jb + ja0
#             amat(ia,ja) = -bm*4.5*gradBhu3*vs(is)    &
#                 +d_11*k_par*vs(is)*c06 - d_11*k_par*vs(is)*c08
#             bmat(ia,ja) = 0.0
#         #
#             ja = 4*nbasis+jb + ja0
#             amat(ia,ja) =   &
#                 -modk_par1*modkpar_hd1*vs(is) &
#                 - taus(is)*(modw_d1*hv6r/ABS(zs(is)) + w_d1*xi*hv6i/zs(is)) &
#             -d_ee*nuei_q1_q1_1  &
#             +xnuei*d_ab*xnu_q1_q1_1
#             bmat(ia,ja) = d_ab
#         #
#             ja = 5*nbasis+jb + ja0
#             amat(ia,ja) = - taus(is)*(modw_d1*hv7r/ABS(zs(is)) + w_d1*xi*hv7i/zs(is)) &
#             -d_ee*nuei_q1_q3_1   &
#             +xnuei*d_ab*xnu_q1_q3_1
#             bmat(ia,ja) = 0.0
#         #
#             if(nroot.gt.6)then
# #
# #  q_par_u ghost terms
# #
#                 ja = 6*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(phi_A + psi_AN)
#                 bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#             #
#                 ja = 7*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(psi_A + phi_AU)
#                 bmat(ia,ja) = -1.0*(psi_B + phi_BU)
#             #
#                 ja = 8*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 9*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 10*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 11*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#             #  q_par_u trapped terms
#             #
#                 ja = 12*nbasis+jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = 13*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 14*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#             endif   # nroot >6
# #
# # q_tot_u equ #6
# #
#             ia = 5*nbasis+ib + ia0
# #
#             phi_A = N_j*xi*w_s*vpar_shear(is)*hr13/vs(is)
#             phi_B = 0.0
#             if(vpar_model_in.eq.0)then
#                 phi_A = phi_A +N_j*xi*w_cd*wdhr13p0*vpar(is)/vs(is) &
#                 + d_1*(nuei_q3_u_1 + (5.0/3.0)*nuei_q3_q3_1)*hr13*E_i*N_j*vpar(is)/vs(is)
#                 phi_B = -E_i*N_j*hr13*vpar(is)/vs(is)
#             endif
#             sig_A = 0.0
#             sig_B = 0.0
#             psi_A = 0.0
#             psi_B = 0.0
#             phi_AU = 0.0
#             phi_BU = 0.0
#             psi_AN = 0.0
#             psi_BN = 0.0
#             if(use_bper_in)then
#                 psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*hr13b0+1.5*rlts(is)*(hw133b0-hr13b0)) 
#                 psi_B = hr13b0*betae_psi*M_i*J_j
#                 psi_A = psi_A - damp_psi*psi_B
#                 if(vpar_model_in.eq.0)then
#                     psi_A = psi_A -betae_psi*J_j*M_i*kpar_hr13b0*vpar(is)
#                     phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear(is)*hr13bp +w_cd*wdhr13bp*vpar(is))/vs(is) &
#                     + d_1*(nuei_q3_u_1 + (5.0/3.0)*nuei_q3_q3_1)*hr13bp*E_i*betae_psi*U0*J_j*vpar(is)/vs(is)
#                     phi_BU = -betae_psi*U0*E_i*J_j*hr13bp*vpar(is)/vs(is)
#                     psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*hr13bp+1.5*rlts(is)*(hw133bp-hr13bp)) &
#                     + M_i*kpar_hr13bp*vpar(is))
#                     psi_BN = -betae_psi*U0*hr13bp*M_i*N_j
#                 endif
#             endif
# #
# #  q_tot_u untrapped terms
# #
#             ja = jb + ja0
#             amat(ia,ja) = phi_A + psi_AN  &
#             +k_par1*(kpar_hb3ht3 -dhr13+ kpar_hb33ht1)*vs(is)
#             bmat(ia,ja) = phi_B + psi_BN
#         #
#             ja = nbasis+jb + ja0
#             amat(ia,ja) =  psi_A  + phi_AU + &
#             modk_par1*(modkpar_hd3hu3 + modkpar_hd33hu1)*vs(is) &
#                 - taus(is)*(modw_d1*hv8r/ABS(zs(is)) + w_d1*xi*hv8i/zs(is)) &
#             -d_ee*nuei_q3_u_1  &
#             +xnuei*(d_ab*xnu_q3_u_1 - d_ij*hu3*xnu_q3_q3_1)
#             bmat(ia,ja) = psi_B + phi_BU
#         #
#             ja = 2*nbasis+jb + ja0
#             amat(ia,ja) = - k_par1*(kpar_hu3 -grad_hu3 - gradhr13p1)*vs(is) & 
#             - k_par1*kpar_hb33*vs(is) + (am+bm*0.5)*gradBhu3*vs(is)      &
#             -d_11*k_par*vs(is)*c07
#             bmat(ia,ja) = 0.0
#         #
#             ja = 3*nbasis+jb + ja0
#             amat(ia,ja) = -k_par1*kpar_hb3*vs(is) - bm*1.5*gradBhu33*vs(is) &
#             +d_11*k_par*vs(is)*c07 
#             bmat(ia,ja) = 0.0
#         #
#             ja = 4*nbasis+jb + ja0
#             amat(ia,ja) = -modk_par1*modkpar_hd33*vs(is) &
#             - taus(is)*(modw_d1*hv9r/ABS(zs(is)) + w_d1*xi*hv9i/zs(is))
#             bmat(ia,ja) = 0.0
#         #
#             ja = 5*nbasis+jb + ja0
#             amat(ia,ja) =   &
#             -modk_par1*modkpar_hd3*vs(is) &
#             - taus(is)*(modw_d1*hv10r/ABS(zs(is)) + w_d1*xi*hv10i/zs(is)) &
#             -d_ee*nuei_q3_q3_1  &
#             +xnuei*d_ab*xnu_q3_q3_1
#             bmat(ia,ja) = d_ab
#         #
#             if(nroot.gt.6)then
# #
# #  q_tot_u ghost terms
# #
#                 ja = 6*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(phi_A + psi_AN)
#                 bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#             #
#                 ja = 7*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(psi_A + phi_AU)
#                 bmat(ia,ja) = -1.0*(psi_B + phi_BU)
#             #
#                 ja = 8*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 9*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 10*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 11*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#             #  q_tot_u trapped terms
#             #
#                 ja = 12*nbasis+jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = 13*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 14*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0

#                 endif    # nroot>6
#             #
#                 if(nroot.gt.6)then
# #
# # n_g equ #7
# #
#                     ia = 6*nbasis+ib + ia0
#                 #
#                     phi_A = N_j*xi*w_s*(rlns(is)*gn + rlts(is)*1.5*(gp3-gn))  
#                     if(vpar_model_in.eq.0)then
#                         phi_A = phi_A + N_j*E_i*kpar_gnp0*vpar(is)
#                     endif
#                     phi_B = -E_i*N_j*gn
#                     phi_A = phi_A +xnu_phi_b*xnuei*xnu_n_b*phi_B
#                     sig_A = 0.0
#                     sig_B = 0.0
#                     psi_A = 0.0
#                     psi_B = 0.0
#                     phi_AU = 0.0
#                     phi_BU = 0.0
#                     psi_AN = 0.0
#                     psi_BN = 0.0
#                     if(use_bpar_in)then
#                         sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
#                         (xi*w_s*(rlns(is)*g10n + rlts(is)*1.5*(g10p3-g10n)))
#                         sig_B = betae_sig*g10n*as(js)*taus(js)*zs(is)*zs(is) &
#                         /(taus(is)*mass(is))
#                         sig_A = sig_A - damp_sig*sig_B
#                     endif
#                     if(use_bper_in)then
#                         psi_A = -betae_psi*J_j*xi*w_s*vpar_shear(is)*gp1b0
#                         if(vpar_model_in.eq.0)then
#                             psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar(is)*wdgp1b0
#                             psi_B = betae_psi*M_i*J_j*vpar(is)*gp1b0/vs(is)
#                             phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*gnbp + rlts(is)*1.5*(gp3bp-gnbp))  &
#                             + E_i*kpar_gnbp*vpar(is))
#                             phi_BU = -betae_psi*U0*E_i*J_j*gnbp
#                             psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar(is)*wdgp1bp+w_s*vpar_shear(is)*gp1bp)
#                             psi_BN = -betae_psi*U0*M_i*N_j*vpar(is)*gp1bp/vs(is)
#                         endif
#                     endif
# #
# #  n_g untrapped terms
# #
#                     ja = jb + ja0
#                     amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_n_n*bn
#                     bmat(ia,ja) = phi_B + psi_BN
#                 #
#                     ja = nbasis+jb + ja0
#                     amat(ia,ja) = psi_A + phi_AU
#                     bmat(ia,ja) = psi_B + phi_BU
#                 #
#                     ja = 2*nbasis+jb + ja0
#                     amat(ia,ja) = -0.5*sig_A +d_ee*nuei_n_p1*bp1  &
#                     -d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1)
#                     bmat(ia,ja) = -0.5*sig_B
#                 #
#                     ja = 3*nbasis+jb + ja0
#                     amat(ia,ja) = 1.5*sig_A +d_ee*nuei_n_p3*bp3  &
#                     +d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1)
#                     bmat(ia,ja) = 1.5*sig_B
#                 #
#                     ja = 4*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 5*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0
#                 #
#                 #  n_g ghost terms
#                 #
#                     ja = 6*nbasis+jb + ja0
#                     amat(ia,ja) = -1.0*(phi_A +psi_AN)  &
#                     -d_ee*nuei_n_n  &
#                     +xnuei*d_ab*xnu_n_b 
#                     bmat(ia,ja) = d_ab - 1.0*(phi_B + psi_BN)
#                 #
#                     ja = 7*nbasis+jb + ja0
#                     amat(ia,ja) = -k_par*vs(is) + am*gradB*vs(is) -1.0*(psi_A + phi_AU)
#                     bmat(ia,ja) = -1.0*(psi_B + phi_BU)
#                 #
#                     ja = 8*nbasis+jb + ja0
#                     amat(ia,ja) = -0.5*(-1.0*sig_A) &
#                         -0.5*xi*w_dg*taus(is)/zs(is)   &
#                         -d_ee*nuei_n_p1
#                     bmat(ia,ja) = -0.5*(-1.0*sig_B)
#                 #
#                     ja = 9*nbasis+jb + ja0
#                     amat(ia,ja) = 1.5*(-1.0*sig_A) &
#                         -1.5*xi*w_dg*taus(is)/zs(is)  &
#                         -d_ee*nuei_n_p3
#                     bmat(ia,ja) = 1.5*(-1.0*sig_B)
#                 #
#                     ja = 10*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 11*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0
#                 #
#                 # n_g trapped terms
#                 #
#                     ja = 12*nbasis+jb + ja0
#                     amat(ia,ja) = phi_A + psi_AN
#                     bmat(ia,ja) = phi_B + psi_BN
#                 #
#                     ja = 13*nbasis+jb + ja0
#                     amat(ia,ja) = -0.5*sig_A 
#                     bmat(ia,ja) = -0.5*sig_B
#                 #
#                     ja = 14*nbasis+jb + ja0
#                     amat(ia,ja) = 1.5*sig_A  
#                     bmat(ia,ja) = 1.5*sig_B
# #
# # u_par_g equ #8
# #
#                     ia = 7*nbasis+ib + ia0
#                 #
#                     phi_A = N_j*xi*w_s*vpar_shear(is)*gp1/vs(is)
#                     if(vpar_model_in.eq.0)then
#                         phi_A = phi_A  + N_j*xi*w_cd*wdgp1p0*vpar(is)/vs(is)
#                         phi_B = -E_i*N_j*gp1*vpar(is)/vs(is)
#                     endif
#                     sig_A = 0.0
#                     sig_B = 0.0
#                     psi_A = 0.0
#                     psi_B = 0.0
#                     phi_AU = 0.0
#                     phi_BU = 0.0
#                     psi_AN = 0.0
#                     psi_BN = 0.0
#                     if(use_bper_in)then
#                         psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*gp1b0+1.5*rlts(is)*(gr13b0-gp1b0))
#                         psi_B =betae_psi*M_i*J_j*gp1b0
#                         psi_A = psi_A - damp_psi*psi_B
#                         if(vpar_model_in.eq.0)then
#                             psi_A = psi_A  -betae_psi*J_j*M_i*kpar_gp1b0*vpar(is)
#                             phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear(is)*gp1bp +w_cd*wdgp1bp*vpar(is))/vs(is)
#                             phi_BU = -betae_psi*U0*E_i*J_j*gp1bp*vpar(is)/vs(is)
#                             psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*gp1bp+1.5*rlts(is)*(gr13bp-gp1bp)) &
#                             + M_i*kpar_gp1bp*vpar(is))
#                             psi_BN = -betae_psi*U0*M_i*N_j*gp1bp
#                         endif
#                     endif
# #
# # u_par_g untrapped terms
# #
#                     ja = jb + ja0
#                     amat(ia,ja) = phi_A + psi_AN
#                     bmat(ia,ja) = phi_B + psi_BN
#                 #
#                     ja = nbasis+jb + ja0
#                     amat(ia,ja) = psi_A + phi_AU + d_ee*ft3*nuei_u_u 
#                     bmat(ia,ja) = psi_B + phi_BU
#                 #
#                     ja = 2*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 3*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 4*nbasis+jb + ja0
#                     amat(ia,ja) = +d_ee*ft5*nuei_u_q1  &
#                             -d_ee*(1.0 -ft2)*(ft3*nuei_u_u*1.25 +ft3*nuei_u_q3*35.0/12.0 +ft5*nuei_u_q1*6.25) 
                                
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 5*nbasis+jb + ja0
#                     amat(ia,ja) = d_ee*ft3*nuei_u_q3  &
#                         +d_ee*(1.0-ft2)*(ft3*nuei_u_u*2.25 +ft3*nuei_u_q3*5.25 +ft5*nuei_u_q1*11.25)  
                                
#                     bmat(ia,ja) = 0.0
#                 #
#                 # u_par_g ghost terms
#                 #
#                     ja = 6*nbasis+jb + ja0
#                     amat(ia,ja) = -1.0*(phi_A + psi_AN)
#                     bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#                 #
#                     ja = 7*nbasis+jb + ja0
#                     amat(ia,ja) = -1.0*(psi_A + phi_AU)  &
#                     -d_ee*nuei_u_u_t -d_ee*nuei_u_u  &
#                     +xnuei*(d_ab*xnu_u_u_1 - d_ij*gu3*xnu_u_q3_1) &
#                     +xnuei*d_ab*xnu_u_b 
#                 #    &   +resist(is,js)
#                     bmat(ia,ja) = d_ab -1.0*(psi_B + phi_BU)
#                 #
#                     ja = 8*nbasis+jb + ja0
#                     amat(ia,ja) = -(k_par - k_par1*gradgp1p1)*vs(is) &
#                     +(am + bm*0.5)*gradB*vs(is) 
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 9*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0 &
#                     -bm*1.5*gradB*vs(is)
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 10*nbasis+jb + ja0
#                     amat(ia,ja) = -0.5*xi*w_dg*taus(is)/zs(is)  &
#                         -d_ee*nuei_u_q1
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 11*nbasis+jb + ja0
#                     amat(ia,ja) = -1.5*xi*w_dg*taus(is)/zs(is) &
#                     -d_ee*nuei_u_q3_t -d_ee*nuei_u_q3  &
#                     +xnuei*d_ab*xnu_u_q3_1
#                     bmat(ia,ja) = 0.0
#                 #
#                 # u_par_g trapped terms
#                 #
#                     ja = 12*nbasis+jb + ja0
#                     amat(ia,ja) = phi_A + psi_AN
#                     bmat(ia,ja) = phi_B + psi_BN
#                 #
#                     ja = 13*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0
#                 #
#                     ja = 14*nbasis+jb + ja0
#                     amat(ia,ja) = 0.0
#                     bmat(ia,ja) = 0.0
#                 #
#                 # p_par_g equ #9
#                 #
#                     ia = 8*nbasis+ib + ia0
#                 #
#                     phi_A = N_j*xi*w_s*(rlns(is)*gp1 + rlts(is)*1.5*(gr13-gp1))  
#                 if(vpar_model_in.eq.0)then
#                     phi_A = phi_A  +N_j*E_i*kpar_gp1p0*vpar(is)
#                 endif
#                 phi_B = -E_i*N_j*gp1
#                 phi_A = phi_A +xnu_phi_b*xnuei*xnu_p1_b*phi_B
#                 sig_A = 0.0
#                 sig_B = 0.0
#                 psi_A = 0.0
#                 psi_B = 0.0
#                 phi_AU = 0.0
#                 phi_BU = 0.0
#                 psi_AN = 0.0
#                 psi_BN = 0.0
#                 if(use_bpar_in)then
#                     sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
#                     (xi*w_s*(rlns(is)*g10p1 + rlts(is)*1.5*(g10r13-g10p1)))
#                     sig_B = betae_sig*g10p1*as(js)*taus(js)*zs(is)*zs(is) &
#                     /(taus(is)*mass(is))
#                     sig_A = sig_A - damp_sig*sig_B
#                 endif
#                 if(use_bper_in)then
#                     psi_A = -betae_psi*J_j*xi*w_s*vpar_shear(is)*gr11b0
#                     if(vpar_model_in.eq.0)then
#                         psi_A = psi_A  -betae_psi*J_j*xi*w_cd*vpar(is)*wdgr11b0
#                         psi_B = betae_psi*M_i*J_j*vpar(is)*gr11b0/vs(is)
#                         phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*gp1bp + rlts(is)*1.5*(gr13bp-gp1bp)) &
#                         + E_i*kpar_gp1bp*vpar(is))
#                         phi_BU = -betae_psi*U0*E_i*J_j*gp1bp
#                         psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar(is)*wdgr11bp+w_s*vpar_shear(is)*gr11bp)
#                         psi_BN = -betae_psi*U0*M_i*N_j*vpar(is)*gr11bp/vs(is)
#                     endif
#                 endif
# #
# # p_par_g untrapped terms
# #
#                 ja = jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_p1_n*bn
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = nbasis+jb + ja0
#                 amat(ia,ja) = psi_A + phi_AU
#                 bmat(ia,ja) = psi_B + phi_BU
#             #
#                 ja = 2*nbasis+jb + ja0
#                 amat(ia,ja) = -0.5*sig_A +d_ee*nuei_p1_p1*bp1  &
#                 -d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1)
#                 bmat(ia,ja) = -0.5*sig_B
#             #
#                 ja = 3*nbasis+jb + ja0
#                 amat(ia,ja) = 1.5*sig_A +d_ee*nuei_p1_p3*bp3    &
#                 +d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1)
#                 bmat(ia,ja) = 1.5*sig_B
#             #
#                 ja = 4*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 5*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#             # p_par_g ghost terms
#             #
#             #
#                 ja = 6*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(phi_A + psi_AN) &
#                     +2.0*taus(is)*(modw_d1*gu1rgt1/ABS(zs(is)) +w_d1*xi*gu1igt1/zs(is)) &
#                     +2.0*taus(is)*(modw_d1*gu2rgt3/ABS(zs(is)) +w_d1*xi*gu2igt3/zs(is)) &
#                     -d_ee*nuei_p1_n  &
#                     -xnuei*d_ij*xnu_p1_1*(gt1-ft2*gt3) 
#                 bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#             #
#                 ja = 7*nbasis+jb + ja0
#                 amat(ia,ja) = k_par1*grad_gu1*vs(is) -1.0*(psi_A + phi_AU)
#                 bmat(ia,ja) = -1.0*(psi_B + phi_BU)
#             #
#                 ja = 8*nbasis+jb + ja0
#                 amat(ia,ja) = -0.5*(-1.0*sig_A)  &
#                     -xi*w_d1*(taus(is)/zs(is))*(0.5*wdgu1+1.5*wdgu3) &
#                     -2.0*taus(is)*(modw_d1*gu1r/ABS(zs(is)) +w_d1*xi*gu1i/zs(is)) & 
#                 -d_ee*nuei_p1_p1_t -d_ee*nuei_p1_p1  &
#                 +xnuei*d_ab*(xnu_p1_1 + xnu_p1_b)
#                 bmat(ia,ja) = d_ab -0.5*(-1.0*sig_B)
#             #
#                 ja = 9*nbasis+jb + ja0
#                 amat(ia,ja) = 1.5*(-1.0*sig_A) &
#                 -2.0*taus(is)*(modw_d1*gu2r/ABS(zs(is)) +w_d1*xi*gu2i/zs(is)) &
#                 -d_ee*nuei_p1_p3_t -d_ee*nuei_p1_p3  &
#                 -xnuei*d_ab*xnu_p1_1*ft2
#                 bmat(ia,ja) = 1.5*(-1.0*sig_B)
#             #
#                 ja = 10*nbasis+jb + ja0
#                 amat(ia,ja) = -k_par*vs(is) +(am +bm)*gradB*vs(is)
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 11*nbasis+jb + ja0
#                 amat(ia,ja) = -bm*3.0*gradB*vs(is)
#                 bmat(ia,ja) = 0.0
#             #
#             # p_par_g trapped terms
#             #
#                 ja = 12*nbasis+jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = 13*nbasis+jb + ja0
#                 amat(ia,ja) = -0.5*sig_A 
#                 bmat(ia,ja) = -0.5*sig_B
#             #
#                 ja = 14*nbasis+jb + ja0
#                 amat(ia,ja) = 1.5*sig_A 
#                 bmat(ia,ja) = 1.5*sig_B
#             #
#             # p_tot_g equ #10
#             #
#                 ia = 9*nbasis+ib + ia0
#             #
#                 phi_A = N_j*xi*w_s*(rlns(is)*gp3+rlts(is)*1.5*(gr33-gp3))
#                 if(vpar_model_in.eq.0)then 
#                     phi_A = phi_A + N_j*E_i*kpar_gp3p0*vpar(is)
#                 endif
#                 phi_B = -gp3*E_i*N_j
#                 phi_A = phi_A +xnu_phi_b*xnuei*xnu_p3_b*phi_B
#                 sig_A = 0.0
#                 sig_B = 0.0
#                 psi_A = 0.0
#                 psi_B = 0.0
#                 phi_AU = 0.0
#                 phi_BU = 0.0
#                 psi_AN = 0.0
#                 psi_BN = 0.0
#                 if(use_bpar_in)then
#                     sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
#                     (xi*w_s*(rlns(is)*g10p3 + rlts(is)*1.5*(g10r33-g10p3)))
#                     sig_B = betae_sig*g10p3*as(js)*taus(js)*zs(is)*zs(is) &
#                     /(taus(is)*mass(is))
#                     sig_A = sig_A - damp_sig*sig_B
#                 endif
#                 if(use_bper_in)then
#                     psi_A = -betae_psi*J_j*xi*w_s*vpar_shear(is)*gr13b0
#                     if(vpar_model_in.eq.0)then
#                         psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar(is)*wdgr13b0
#                         psi_B = betae_psi*M_i*J_j*vpar(is)*gr13b0/vs(is)
#                         phi_AU = betae_psi*U0*J_j*(xi*w_s*(rlns(is)*gp3bp+rlts(is)*1.5*(gr33bp-gp3bp))  &
#                         + E_i*kpar_gp3bp*vpar(is))
#                         phi_BU = -betae_psi*U0*gp3bp*E_i*J_j
#                         psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar(is)*wdgr13bp+w_s*vpar_shear(is)*gr13bp)
#                         psi_BN = -betae_psi*U0*M_i*N_j*vpar(is)*gr13bp/vs(is)
#                     endif
#                 endif
# #
# #  p_tot_g untrapped terms
# #
#                 ja = jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_p3_n*bn
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = nbasis+jb + ja0
#                 amat(ia,ja) = psi_A + phi_AU
#                 bmat(ia,ja) = psi_B + phi_BU
#             #
#                 ja = 2*nbasis+jb + ja0
#                 amat(ia,ja) = -0.5*sig_A +d_ee*nuei_p3_p1*bp1  &
#                 -d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1)
#                 bmat(ia,ja) = -0.5*sig_B
#             #
#                 ja = 3*nbasis+jb + ja0
#                 amat(ia,ja) = 1.5*sig_A +d_ee*nuei_p3_p3*bp3    &
#                 +d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1)
#                 bmat(ia,ja) = 1.5*sig_B
#             #
#                 ja = 4*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 5*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#             #  p_tot_g ghost terms
#             #
#                 ja = 6*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(phi_A + psi_AN) &
#                 +2.0*taus(is)*(modw_d1*gu3rgt1/ABS(zs(is)) +w_d1*xi*gu3igt1/zs(is)) &
#                 +2.0*taus(is)*(modw_d1*gu4rgt3/ABS(zs(is)) +w_d1*xi*gu4igt3/zs(is)) &
#                 -d_ee*nuei_p3_n
#                 bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#             #
#                 ja = 7*nbasis+jb + ja0
#                 amat(ia,ja) = k_par1*grad_gu3*vs(is) -1.0*(psi_A + phi_AU)
#                 bmat(ia,ja) = -1.0*(psi_B + phi_BU)
#             #
#                 ja = 8*nbasis+jb + ja0
#                 amat(ia,ja) = -0.5*(-1.0*sig_A) &
#                     -xi*w_d1*(taus(is)/zs(is))*0.5*wdgu3 & 
#                     -2.0*taus(is)*(modw_d1*gu3r/ABS(zs(is)) +w_d1*xi*gu3i/zs(is))  &
#                     -d_ee*nuei_p3_p1
#                 bmat(ia,ja) = -0.5*(-1.0*sig_B)
#             #
#                 ja = 9*nbasis+jb + ja0
#                 amat(ia,ja) = 1.5*(-1.0*sig_A)  &
#                 -xi*w_d1*(taus(is)/zs(is))*1.5*wdgu33 &
#                 -2.0*taus(is)*(modw_d1*gu4r/ABS(zs(is)) +w_d1*xi*gu4i/zs(is)) &
#                 -d_ee*nuei_p3_p3  &
#                 +xnuei*d_ab*xnu_p3_b
#                 bmat(ia,ja) = d_ab + 1.5*(-1.0*sig_B)
#             #
#                 ja = 10*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 11*nbasis+jb + ja0
#                 amat(ia,ja) = -k_par*vs(is) + am*gradB*vs(is)
#                 bmat(ia,ja) = 0.0
#             #
#             #  p_tot_g trapped terms
#             #
#                 ja = 12*nbasis+jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = 13*nbasis+jb + ja0
#                 amat(ia,ja) = -0.5*sig_A 
#                 bmat(ia,ja) = -0.5*sig_B
#             #
#                 ja = 14*nbasis+jb + ja0
#                 amat(ia,ja) = 1.5*sig_A 
#                 bmat(ia,ja) = 1.5*sig_B
#             #
#             #  q_par_g equ #11
#             #
#                 ia = 10*nbasis+ib + ia0
#             #
#                 phi_A = N_j*xi*w_s*vpar_shear(is)*gr11/vs(is)
#                 if(vpar_model_in.eq.0)then
#                     phi_A = phi_A + N_j*xi*w_cd*wdgr11p0*vpar(is)/vs(is)
#                     phi_B = -E_i*N_j*gr11*vpar(is)/vs(is)
#                 endif
#                 sig_A = 0.0
#                 sig_B = 0.0
#                 psi_A = 0.0
#                 psi_B = 0.0
#                 phi_AU = 0.0
#                 phi_BU = 0.0
#                 psi_AN = 0.0
#                 psi_BN = 0.0
#                 if(use_bper_in)then
#                     psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*gr11b0+1.5*rlts(is)*(gw113b0-gr11b0))
#                     psi_B = gr11b0*betae_psi*M_i*J_j
#                     psi_A = psi_A - damp_psi*psi_B
#                     if(vpar_model_in.eq.0)then
#                         psi_A = psi_A  -betae_psi*J_j*M_i*kpar_gr11b0*vpar(is)
#                         phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear(is)*gr11bp +w_cd*wdgr11bp*vpar(is))/vs(is)
#                         phi_BU = -betae_psi*U0*E_i*J_j*gr11bp*vpar(is)/vs(is)
#                         psi_AN = -betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*gr11bp+1.5*rlts(is)*(gw113bp-gr11bp)) &
#                         +M_i*kpar_gr11bp*vpar(is))
#                         psi_BN = -gr11bp*betae_psi*U0*M_i*N_j
#                     endif
#                 endif
# #
# # q_par_g untrapped terms
# #
#                 ja = jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = nbasis+jb + ja0
#                 amat(ia,ja) = psi_A + phi_AU +d_ee*ft3*nuei_q1_u
#                 bmat(ia,ja) = psi_B + phi_BU
#             #
#                 ja = 2*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 3*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 4*nbasis+jb + ja0
#                 amat(ia,ja) = +d_ee*ft5*nuei_q1_q1 &
#                     -d_ee*(1.0 -ft2)*(ft3*nuei_q1_u*1.25 +ft3*nuei_q1_q3*35.0/12.0+ft5*nuei_q1_q1*6.25)
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 5*nbasis+jb + ja0
#                 amat(ia,ja) = d_ee*ft3*nuei_q1_q3  &
#                 +d_ee*(1.0 -ft2)*(ft3*nuei_q1_u*2.25 +ft3*nuei_q1_q3*5.25 +ft5*nuei_q1_q1*11.25)
#                 bmat(ia,ja) = 0.0
#             #
#             # q_par_g ghost terms
#             #
#                 ja = 6*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(phi_A + psi_AN) + k_par1*kpar_gb1gt1*vs(is) 
#                 bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#             #
#                 ja = 7*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(psi_A + phi_AU) +modk_par1*modkpar_gd1gu1*vs(is)  &
#                     - taus(is)*(modw_d1*gu5r/ABS(zs(is)) + w_d1*xi*gu5i/zs(is)) &
#                     -d_ee*nuei_q1_u_t -d_ee*nuei_q1_u &
#                     +xnuei*(d_ab*xnu_q1_u_1*ft2 - d_ij*gu1*xnu_q1_q1_1 &
#                     - d_ij*gu3*xnu_q1_q3_1*ft2)
#                 bmat(ia,ja) = -1.0*(psi_B + phi_BU)
#             #
#                 ja = 8*nbasis+jb + ja0
#                 amat(ia,ja) =  -k_par1*(kpar_gu1 -grad_gu1 - gradgr11p1)*vs(is) &
#                     - k_par1*kpar_gb1*vs(is) + (am+bm*1.5)*gradBgu1*vs(is)
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 9*nbasis+jb + ja0
#                 amat(ia,ja) = -bm*4.5*gradBgu3*vs(is) 
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 10*nbasis+jb + ja0
#                 amat(ia,ja) =    &
#                     -modk_par1*modkpar_gd1*vs(is) &
#                     - taus(is)*(modw_d1*gu6r/ABS(zs(is)) + w_d1*xi*gu6i/zs(is)) &
#                     -d_ee*nuei_q1_q1_t -d_ee*nuei_q1_q1  &
#                     +xnuei*d_ab*(xnu_q1_q1_1 + xnu_q1_b)
#                 bmat(ia,ja) = d_ab
#             #
#                 ja = 11*nbasis+jb + ja0
#                 amat(ia,ja) = - taus(is)*(modw_d1*gu7r/ABS(zs(is)) + w_d1*xi*gu7i/zs(is)) &
#                     -d_ee*nuei_q1_q3_t -d_ee*nuei_q1_q3  &
#                     +xnuei*d_ab*xnu_q1_q3_1*ft2
#                 bmat(ia,ja) = 0.0
#             #
#             # q_par_g trapped terms
#             #
#                 ja = 12*nbasis+jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = 13*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 14*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#             # q_tot_g equ #12
#             #
#                 ia = 11*nbasis+ib + ia0
#             #
#                 phi_A = N_j*xi*w_s*vpar_shear(is)*gr13/vs(is)
#                 if(vpar_model_in.eq.0)then
#                     phi_A = phi_A + N_j*xi*w_cd*wdgr13p0*vpar(is)/vs(is)
#                     phi_B = -E_i*N_j*gr13*vpar(is)/vs(is)
#                 endif
#                 sig_A = 0.0
#                 sig_B = 0.0
#                 psi_A = 0.0
#                 psi_B = 0.0
#                 phi_AU = 0.0
#                 phi_BU = 0.0
#                 psi_AN = 0.0
#                 psi_BN = 0.0
#                 if(use_bper_in)then
#                     psi_A = -betae_psi*J_j*vs(is)*xi*w_s*(rlns(is)*gr13b0+1.5*rlts(is)*(gw133b0-gr13b0))
#                     psi_B = gr13b0*betae_psi*M_i*J_j
#                     psi_A = psi_A - damp_psi*psi_B
#                     if(vpar_model_in.eq.0)then
#                         psi_A = psi_A  -betae_psi*J_j*M_i*kpar_gr13b0*vpar(is)
#                         phi_AU = betae_psi*U0*J_j*xi*(w_s*vpar_shear(is)*gr13bp +w_cd*wdgr13bp*vpar(is))/vs(is)
#                         phi_BU = -betae_psi*U0*E_i*J_j*gr13bp*vpar(is)/vs(is)
#                         psi_AN = betae_psi*U0*N_j*(vs(is)*xi*w_s*(rlns(is)*gr13bp+1.5*rlts(is)*(gw133bp-gr13bp)) &
#                         +M_i*kpar_gr13bp*vpar(is))
#                         psi_BN = -gr13bp*betae_psi*U0*M_i*N_j
#                     endif
#                 endif
# #
# # q_tot_g untrapped terms
# #
#                 ja = jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = nbasis+jb + ja0
#                 amat(ia,ja) = psi_A + phi_AU + d_ee*ft3*nuei_q3_u
#                 bmat(ia,ja) = psi_B + phi_BU
#             #
#                 ja = 2*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 3*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 4*nbasis+jb + ja0
#                 amat(ia,ja) = d_ee*ft5*nuei_q3_q1   &
#                 -d_ee*(1.0 -ft2)*(ft3*nuei_q3_u*1.25 +ft3*nuei_q3_q3*35.0/12.0 +ft5*nuei_q3_q1*6.25)
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 5*nbasis+jb + ja0
#                 amat(ia,ja) = d_ee*ft3*nuei_q3_q3  &
#                 +d_ee*(1.0 -ft2)*(ft3*nuei_q3_u*2.25 + ft3*nuei_q3_q3*5.25 + ft5*nuei_q3_q1*11.25)
#                 bmat(ia,ja) = 0.0
#             #
#             # q_tot_g ghost terms
#             #
#                 ja = 6*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(phi_A + psi_AN) + k_par1*(kpar_gb3gt3 -dgr13 + kpar_gb33gt1)*vs(is)
#                 bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#             #
#                 ja = 7*nbasis+jb + ja0
#                 amat(ia,ja) = -1.0*(psi_A + phi_AU) +  &
#                 modk_par1*(modkpar_gd3gu3 + modkpar_gd33gu1)*vs(is) &
#                     - taus(is)*(modw_d1*gu8r/ABS(zs(is)) + w_d1*xi*gu8i/zs(is)) &
#                 -d_ee*nuei_q3_u_t -d_ee*nuei_q3_u  &
#                 +xnuei*(d_ab*xnu_q3_u_1 - d_ij*gu3*xnu_q3_q3_1)
#                 bmat(ia,ja) = -1.0*(psi_B + phi_BU)
#             #
#                 ja = 8*nbasis+jb + ja0
#                 amat(ia,ja) = - k_par1*(kpar_gu3 -grad_gu3 - gradgr13p1)*vs(is) &
#                 - k_par1*kpar_gb33*vs(is) + (am+bm*0.5)*gradBgu3*vs(is)
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 9*nbasis+jb + ja0
#                 amat(ia,ja) = &
#                     - k_par1*kpar_gb3*vs(is) - bm*1.5*gradBgu33*vs(is) 
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 10*nbasis+jb + ja0
#                 amat(ia,ja) = -modk_par1*modkpar_gd33*vs(is) &
#                 - taus(is)*(modw_d1*gu9r/ABS(zs(is)) + w_d1*xi*gu9i/zs(is))  &
#                 -d_ee*nuei_q3_q1
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 11*nbasis+jb + ja0
#                 amat(ia,ja) =   &
#                     -modk_par1*modkpar_gd3*vs(is) &
#                     - taus(is)*(modw_d1*gu10r/ABS(zs(is)) +w_d1* xi*gu10i/zs(is)) &
#                     -d_ee*nuei_q3_q3_t -d_ee*nuei_q3_q3  &
#                     +xnuei*d_ab*(xnu_q3_q3_1 + xnu_q3_b)
#                 bmat(ia,ja) = d_ab
#             #
#             # q_tot_g trapped terms
#             #
#                 ja = 12*nbasis+jb + ja0
#                 amat(ia,ja) = phi_A + psi_AN
#                 bmat(ia,ja) = phi_B + psi_BN
#             #
#                 ja = 13*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#                 ja = 14*nbasis+jb + ja0
#                 amat(ia,ja) = 0.0
#                 bmat(ia,ja) = 0.0
#             #
#             # n_t equ #13
#             #
#                 ia = 12*nbasis+ib + ia0
#             #
#                 phi_A = N_j*xi*w_s*(rlns(is)*gn + rlts(is)*1.5*(gp3-gn))       
#                 phi_B = -gn*E_i*N_j
#                 phi_A = phi_A +xnu_phi_b*xnuei*xnu_n_b*phi_B
#                 sig_A = 0.0
#                 sig_B = 0.0
#                 psi_A = 0.0
#                 psi_B = 0.0
#                 phi_AU = 0.0
#                 phi_BU = 0.0
#                 psi_AN = 0.0
#                 psi_BN = 0.0
#             if(use_bpar_in)then
#                 sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
#                 (xi*w_s*(rlns(is)*g10n + rlts(is)*1.5*(g10p3-g10n)))
#                 sig_B = betae_sig*g10n*as(js)*taus(js)*zs(is)*zs(is) &
#                 /(taus(is)*mass(is))
#                 sig_A = sig_A - damp_sig*sig_B
#             endif
#             if(use_bper_in)then
#                 psi_A = -betae_psi*J_j*xi*w_s*vpar_shear(is)*gp1b0
#                 if(vpar_model_in.eq.0)then
#                     psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar(is)*wdgp1b0
#                     psi_B = betae_psi*M_i*J_j*vpar(is)*gp1b0/vs(is)
#                     phi_AU = betae_psi*U0*J_j*xi*w_s*(rlns(is)*gnbp + rlts(is)*1.5*(gp3bp-gnbp))
#                     phi_BU = -betae_psi*U0*gnbp*E_i*J_j
#                     psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar(is)*wdgp1bp+w_s*vpar_shear(is)*gp1bp)
#                     psi_BN = -betae_psi*U0*M_i*N_j*vpar(is)*gp1bp/vs(is)
#                 endif
#             endif
# #
# # n_t untrapped terms
# #
#             ja = jb + ja0
#             amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_n_n*bn
#             bmat(ia,ja) = phi_B + psi_BN
#         #
#             ja = nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 2*nbasis+jb + ja0
#             amat(ia,ja) = -0.5*sig_A +d_ee*nuei_n_p1*bp1  &
#             -d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1)
#             bmat(ia,ja) = -0.5*sig_B
#         #
#             ja = 3*nbasis+jb + ja0
#             amat(ia,ja) = 1.5*sig_A +d_ee*nuei_n_p3*bp3    &
#             +d_ee*(1.0 - ft2)*(bn*an*nuei_n_n + bp3*ap3*nuei_n_p3 + bp1*ap1*nuei_n_p1)
#             bmat(ia,ja) = 1.5*sig_B
#         #
#             ja = 4*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 5*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#         # n_t ghost terms
#         #
#             ja = 6*nbasis+jb + ja0
#             amat(ia,ja) = -1.0*(phi_A + psi_AN)
#             bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#         #
#             ja = 7*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 8*nbasis+jb + ja0
#             amat(ia,ja) = -0.5*(-1.0*sig_A) 
#             bmat(ia,ja) = -0.5*(-1.0*sig_B)
#         #
#             ja = 9*nbasis+jb + ja0
#             amat(ia,ja) = 1.5*(-1.0*sig_A) 
#             bmat(ia,ja) = 1.5*(-1.0*sig_B)
#         #
#             ja = 10*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 11*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#         # n_t trapped terms
#         #
#             ja = 12*nbasis+jb + ja0
#             amat(ia,ja) = phi_A + psi_AN   &
#                 -d_ee*nuei_n_n  &
#                 +xnuei*d_ab*xnu_n_b 
#             bmat(ia,ja) = d_ab + phi_B + psi_BN
#         #
#             ja = 13*nbasis+jb + ja0
#             amat(ia,ja) = -0.5*sig_A -0.5*xi*w_dg*taus(is)/zs(is)  &
#                 -d_ee*nuei_n_p1
#             bmat(ia,ja) = -0.5*sig_B
#         #
#             ja = 14*nbasis+jb + ja0
#             amat(ia,ja) = 1.5*sig_A -1.5*xi*w_dg*taus(is)/zs(is)  &
#             -d_ee*nuei_n_p3
#             bmat(ia,ja) = 1.5*sig_B
#         #
#         # p_par_t equ #14
#         #
#             ia = 13*nbasis+ib + ia0
#         #
#             phi_A = N_j*xi*w_s*(rlns(is)*gp1 + rlts(is)*1.5*(gr13-gp1))     
#             phi_B = -gp1*E_i*N_j
#             phi_A = phi_A +xnu_phi_b*xnuei*xnu_p1_b*phi_B
#             sig_A = 0.0
#             sig_B = 0.0
#             psi_A = 0.0
#             psi_B = 0.0
#             phi_AU = 0.0
#             phi_BU = 0.0
#             psi_AN = 0.0
#             psi_BN = 0.0
#             if(use_bpar_in)then
#                 sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
#                 (xi*w_s*(rlns(is)*g10p1 + rlts(is)*1.5*(g10r13-g10p1)))
#                 sig_B = betae_sig*g10p1*as(js)*taus(js)*zs(is)*zs(is) &
#                 /(taus(is)*mass(is))
#                 sig_A = sig_A - damp_sig*sig_B
#             endif
#             if(use_bper_in)then
#                 psi_A = -betae_psi*J_j*xi*w_s*vpar_shear(is)*gr11b0
#                 if(vpar_model_in.eq.0)then
#                     psi_A = psi_A  -betae_psi*J_j*xi*w_cd*vpar(is)*wdgr11b0
#                     psi_B = betae_psi*M_i*J_j*vpar(is)*gr11b0/vs(is)
#                     phi_AU = betae_psi*U0*J_j*xi*w_s*(rlns(is)*gp1bp + rlts(is)*1.5*(gr13bp-gp1bp))
#                     phi_BU = -betae_psi*U0*gp1bp*E_i*J_j
#                     psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar(is)*wdgr11bp+ w_s*vpar_shear(is)*gr11bp)
#                     psi_BN = -betae_psi*U0*M_i*N_j*vpar(is)*gr11bp/vs(is)
#                 endif
#             endif
# #
# # p_par_t  untrapped terms
# #
#             ja = jb + ja0
#             amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_p1_n*bn
#             bmat(ia,ja) = phi_B + psi_BN
#         #
#             ja = nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 2*nbasis+jb + ja0
#             amat(ia,ja) = -0.5*sig_A +d_ee*nuei_p1_p1*bp1  &
#             -d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1)
#             bmat(ia,ja) = -0.5*sig_B
#         #
#             ja = 3*nbasis+jb + ja0
#             amat(ia,ja) = 1.5*sig_A
#             bmat(ia,ja) = 1.5*sig_B +d_ee*nuei_p1_p3*bp3   &
#             +d_ee*(1.0 - ft2)*(bn*an*nuei_p1_n + bp3*ap3*nuei_p1_p3 + bp1*ap1*nuei_p1_p1)
#         #
#             ja = 4*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 5*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#         # p_par_t  ghost terms
#         #
#             ja = 6*nbasis+jb + ja0
#             amat(ia,ja) = -1.0*(phi_A + psi_AN) 
#             bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#         #
#             ja = 7*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 8*nbasis+jb + ja0
#             amat(ia,ja) = -0.5*(-1.0*sig_A)
#             bmat(ia,ja) = -0.5*(-1.0*sig_B)
#         #
#             ja = 9*nbasis+jb + ja0
#             amat(ia,ja) = 1.5*(-1.0*sig_A)
#             bmat(ia,ja) = 1.5*(-1.0*sig_B)
#         #
#             ja = 10*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 11*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#         # p_par_t  trapped terms
#         #
#         #
#             ja = 12*nbasis+jb + ja0
#             amat(ia,ja) = phi_A +psi_AN &
#                 +2.0*taus(is)*(modw_d1*gu1rgt1/ABS(zs(is)) +w_d1*xi*gu1igt1/zs(is)) &
#                 +2.0*taus(is)*(modw_d1*gu2rgt3/ABS(zs(is)) +w_d1*xi*gu2igt3/zs(is)) &
#                 -d_ee*nuei_p1_n  &
#                 -xnuei*d_ij*xnu_p1_1*(gt1-ft2*gt3) 
#             bmat(ia,ja) = phi_B + psi_BN
#         #
#             ja = 13*nbasis+jb + ja0
#             amat(ia,ja) =  -0.5*sig_A   &
#                 -xi*w_d1*(taus(is)/zs(is))*(0.5*wdgu1+1.5*wdgu3) &
#                 -2.0*taus(is)*(modw_d1*gu1r/ABS(zs(is)) +w_d1*xi*gu1i/zs(is)) &
#                 -d_ee*nuei_p1_p1_t -d_ee*nuei_p1_p1  &
#                 +xnuei*d_ab*(xnu_p1_1 + xnu_p1_b)
#             bmat(ia,ja) = d_ab -0.5*sig_B
#         #
#             ja = 14*nbasis+jb + ja0
#             amat(ia,ja) = 1.5*sig_A  &   
#             -2.0*taus(is)*(modw_d1*gu2r/ABS(zs(is)) +w_d1*xi*gu2i/zs(is)) &
#             -d_ee*nuei_p1_p3_t -d_ee*nuei_p1_p3  &
#             -xnuei*d_ab*xnu_p1_1*ft2
#             bmat(ia,ja) = 1.5*sig_B
#         #
#         # p_tot_t equ #15
#         #
#             ia = 14*nbasis+ib + ia0
#         #
#             phi_A = N_j*xi*w_s*(rlns(is)*gp3+rlts(is)*1.5*(gr33-gp3)) 
#             phi_B = -gp3*E_i*N_j
#             phi_A = phi_A +xnu_phi_b*xnuei*xnu_p3_b*phi_B
#             sig_A = 0.0
#             sig_B = 0.0
#             psi_A = 0.0
#             psi_B = 0.0
#             phi_AU = 0.0
#             phi_BU = 0.0
#             psi_AN = 0.0
#             psi_BN = 0.0
#             if(use_bpar_in)then
#                 sig_A = -betae_sig*(as(js)*taus(js)*zs(is)/mass(is))* &
#                 (xi*w_s*(rlns(is)*g10p3 + rlts(is)*1.5*(g10r33-g10p3)))
#                 sig_B = betae_sig*g10p3*as(js)*taus(js)*zs(is)*zs(is) &
#                 /(taus(is)*mass(is))
#                 sig_A = sig_A - damp_sig*sig_B
#             endif
#             if(use_bper_in)then
#                 psi_A = -betae_psi*J_j*xi*w_s*vpar_shear(is)*gr13b0
#                 if(vpar_model_in.eq.0)then
#                     psi_A = psi_A -betae_psi*J_j*xi*w_cd*vpar(is)*wdgr13b0
#                     psi_B = betae_psi*M_i*J_j*vpar(is)*gr13b0/vs(is)
#                     phi_AU = betae_psi*U0*J_j*xi*w_s*(rlns(is)*gp3bp+rlts(is)*1.5*(gr33bp-gp3bp)) 
#                     phi_BU = -betae_psi*U0*gp3bp*E_i*J_j
#                     psi_AN = betae_psi*U0*N_j*xi*(w_cd*vpar(is)*wdgr13bp+w_s*vpar_shear(is)*gr13bp)
#                     psi_BN = -betae_psi*U0*M_i*N_j*vpar(is)*gr13bp/vs(is)
#                 endif
#             endif
# #
# # p_tot_t untrapped terms
# #
#             ja = jb + ja0
#             amat(ia,ja) = phi_A + psi_AN +d_ee*nuei_p3_n*bn
#             bmat(ia,ja) = phi_B + psi_BN
#         #
#             ja = nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 2*nbasis+jb + ja0
#             amat(ia,ja) = -0.5*sig_A +d_ee*nuei_p3_p1*bp1  &
#             -d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1)
#             bmat(ia,ja) = -0.5*sig_B
#         #
#             ja = 3*nbasis+jb + ja0
#             amat(ia,ja) = 1.5*sig_A +d_ee*nuei_p3_p3*bp3    &
#             +d_ee*(1.0 - ft2)*(bn*an*nuei_p3_n + bp3*ap3*nuei_p3_p3 + bp1*ap1*nuei_p3_p1)
#             bmat(ia,ja) = 1.5*sig_B
#         #
#             ja = 4*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 5*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#         # p_tot_t ghost terms
#         #
#             ja = 6*nbasis+jb + ja0
#             amat(ia,ja) = -1.0*(phi_A + psi_AN) 
#             bmat(ia,ja) = -1.0*(phi_B + psi_BN)
#         #
#             ja = 7*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 8*nbasis+jb + ja0
#             amat(ia,ja) = -0.5*(-1.0*sig_A) 
#             bmat(ia,ja) = -0.5*(-1.0*sig_B)
#         #
#             ja = 9*nbasis+jb + ja0
#             amat(ia,ja) = 1.5*(-1.0*sig_A) 
#             bmat(ia,ja) = 1.5*(-1.0*sig_B)
#         #
#             ja = 10*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#             ja = 11*nbasis+jb + ja0
#             amat(ia,ja) = 0.0
#             bmat(ia,ja) = 0.0
#         #
#         # p_tot_t trapped terms
#         #
#             ja = 12*nbasis+jb + ja0
#             amat(ia,ja) = phi_A + psi_AN &
#             +2.0*taus(is)*(modw_d1*gu3rgt1/ABS(zs(is)) +w_d1*xi*gu3igt1/zs(is)) &
#             +2.0*taus(is)*(modw_d1*gu4rgt3/ABS(zs(is)) +w_d1*xi*gu4igt3/zs(is)) &
#             -d_ee*nuei_p3_n 
#             bmat(ia,ja) = phi_B + psi_BN
#         #
#             ja = 13*nbasis+jb + ja0
#             amat(ia,ja) = -0.5*sig_A -xi*w_d1*(taus(is)/zs(is))*0.5*wdgu3   &
#                 -2.0*taus(is)*(modw_d1*gu3r/ABS(zs(is)) +w_d1*xi*gu3i/zs(is)) &
#                 -d_ee*nuei_p3_p1
#             bmat(ia,ja) = -0.5*sig_B
#         #
#             ja = 14*nbasis+jb + ja0
#             amat(ia,ja) =    &
#             +1.5*sig_A -xi*w_d1*(taus(is)/zs(is))*1.5*wdgu33 &
#             -2.0*taus(is)*(modw_d1*gu4r/ABS(zs(is)) +w_d1*xi*gu4i/zs(is)) &
#             -d_ee*nuei_p3_p3  &
#             +xnuei*d_ab*xnu_p3_b
#             bmat(ia,ja) = d_ab + 1.5*sig_B
# #
#                 endif # nroot .gt. 6
# #
#             enddo # end of js loop 
#         enddo # end of ib loop
#     enddo # end of jb loop
# enddo # end of is loop
# #
# #
# #..find the eigenvalues and eigenvectors
# #
# #
#       lwork=33*iur
#       rightvectors="N"
# #      if(iflux_in)rightvectors="V"
#       do i=1,iur
#       do j=1,iur
#         at(i,j) = amat(i,j)
#         bt(i,j) = bmat(i,j)
# #        write(*,*)i,j,at(i,j),bt(i,j)
#       enddo
#       enddo
    
    fileDirectory = "./matrixA"
    lines = readlines(fileDirectory)
    iur = parse(Int,split(lines[1])[3])
    matrixA = Vector{ComplexF64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        realPart = parse(Float64,split(strip(line, [' ', '(',')']),",")[1])
        imagPart = parse(Float64,split(strip(line, [' ', '(',')']),",")[2])
        push!(matrixA,realPart + imagPart*im)
    end
    matrixA = Matrix(reshape(matrixA, iur, iur)')

    fileDirectory = "./matrixB"
    lines = readlines(fileDirectory)
    iur = parse(Int,split(lines[1])[3])
    matrixB = Vector{ComplexF64}()
    for line in lines[2:length(lines)]
        if contains(line,"=") continue end
        realPart = parse(Float64,split(strip(line, [' ', '(',')']),",")[1])
        imagPart = parse(Float64,split(strip(line, [' ', '(',')']),",")[2])
        push!(matrixB,realPart + imagPart*im)
    end
    matrixB = Matrix(reshape(matrixB, iur, iur)')

    (alpha, beta, vl, vr) = ggev!('N','N',matrixA,matrixB)

    zomega = zeros(ComplexF64, iur)
    rr = zeros(Float64, iur)
    ri = zeros(Float64, iur)
    #### not sure what vr and vi are for
    vr = zeros(iur, iur)
    vi = zeros(iur, iur)
    for j1 = 1:iur
        
        beta2 = real(conj(beta[j1])*beta[j1])
        if(beta2!=0.0)
            zomega[j1] = alpha[j1]*conj(beta[j1])/beta2
        else
            zomega[j1]=0.0
        end

        rr[j1] = real(zomega[j1])
        ri[j1] = imag(zomega[j1])
        # filter out numerical instabilities that sometimes occur 
        # with high mode frequency
        if(inputs.FILTER>0.0)
            max_freq = 47.280929835746996####### this is hard coded for now
            if(rr[j1]>0.0 && abs(ri[j1])>max_freq)
                rr[j1]=-rr[j1]
            end
        end
        #### not sure what this is used for
        for j2 = 1:iur
            vr[j1,j2] = 0.0
            vi[j1,j2] = 0.0
        end
    end

    return matrixA, matrixB, alpha, beta, rr, ri
end