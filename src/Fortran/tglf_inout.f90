!-----------------------------------------------------------------
!  input routines
!-----------------------------------------------------------------
SUBROUTINE put_units(units)
!
USE tglf_global
!
IMPLICIT NONE
CHARACTER (len=*) :: units

if(units .eq. 'GYRO' .or. units .eq. 'CGYRO' .or. units .eq. 'GENE')then
   units_in = units
else
   call tglf_error(1,"Invalid units selected")
endif
!
END SUBROUTINE put_units
!
!-----------------------------------------------------------------
!
SUBROUTINE put_species(nsp,zsp,msp)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN):: nsp
  REAL,INTENT(IN) :: zsp(nsm),msp(nsm)
  INTEGER :: is
  !
  use_default_species=.FALSE.
! check for valid input
  if(nsp.lt.2.or.nsp.gt.nsm)call tglf_error(1,"input number of species invlaid")
  do is=1,nsp
    if(tglf_isnan(zsp(is)))call tglf_error(1,"input zs_in is NAN")
    if(tglf_isinf(zsp(is)))call tglf_error(1,"input zs_in is INF") 
    if(tglf_isnan(msp(is)))call tglf_error(1,"input mass_in is NAN")
    if(tglf_isinf(msp(is)))call tglf_error(1,"input mass_in is INF")
    if(zsp(is).eq.0.0)call tglf_error(1,"input zs_in is = 0")
    if(msp(is).le.0.0)call tglf_error(1,"input mass_in is <= 0") 
  enddo
  ! transfer values
  ns_in = nsp
  do is=1,nsp
      zs_in(is)=zsp(is)
      mass_in(is)=msp(is)
  enddo
  new_matrix = .TRUE.
  !
END SUBROUTINE put_species
!
!-----------------------------------------------------------------
!
SUBROUTINE put_kys(kys)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: kys

  if(tglf_isnan(kys))call tglf_error(1,"input ky_in is NAN")
  if(tglf_isinf(kys))call tglf_error(1,"input ky_in is INF")
  if(kys.le.0)call tglf_error(1,"input kys is <= 0")
  !
  ! check for changes and update flow controls
  ! 
  if(kys.ne.ky_in)new_matrix = .TRUE.
  !
  ! transfer values
  !
  ky_in = kys
  !
END SUBROUTINE put_kys
!
!-----------------------------------------------------------------
!
SUBROUTINE put_signs(sign_Bt,sign_It)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: sign_Bt,sign_It

  if(tglf_isnan(sign_Bt))call tglf_error(1,"input sign_Bt_in is NAN")
  if(tglf_isnan(sign_It))call tglf_error(1,"input sign_It_in is NAN")
  !
  ! transfer values
  !
  sign_Bt_in = sign_Bt
  sign_It_in = sign_It
  !
END SUBROUTINE put_signs
!
!-----------------------------------------------------------------
!
SUBROUTINE put_gaussian_width(width,width_min,nwidth,find_width)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  LOGICAL,INTENT(IN) :: find_width
  INTEGER,INTENT(IN) :: nwidth
  REAL,INTENT(IN) :: width,width_min

  if(tglf_isnan(width))call tglf_error(1,"input width_in is NAN")
  if(tglf_isinf(width))call tglf_error(1,"input width_in is INF")
  if(tglf_isnan(width_min))call tglf_error(1,"input width_min_in is NAN")
  if(tglf_isinf(width_min))call tglf_error(1,"input width_min_in is INF")
  if(nwidth.le.0)call tglf_error(1,"input nwidth_in <= 0")

  !
  ! check for changes and update flow controls
  ! 
  if(width.ne.width_in)new_width = .TRUE.
  !
  ! transfer values
  !
  width_in = width
  width_min_in=width_min
  find_width_in = find_width
  nwidth_in=MIN(nwidth,nt0)
  !
END SUBROUTINE put_gaussian_width
!
!
SUBROUTINE put_eikonal(new_eikonal)
  !*********************************************
  !
  !*********************************************
  USE tglf_global
  !
  IMPLICIT NONE
  LOGICAL,INTENT(IN) :: new_eikonal
  !
  !  set flow control switch
  !
  new_eikonal_in = new_eikonal
  ! check consistency
  if(new_eikonal_in)then
     eikonal_unsaved=.TRUE.
  else
     if(eikonal_unsaved)then
        write(*,*)"warning put_eikonal:"
        write(*,*)"new_eikonal = .FALSE.attempted before call with new_eikonal=.TRUE."
        new_eikonal_in = .TRUE.
     endif
  endif
  !
END SUBROUTINE put_eikonal
!
!-----------------------------------------------------------------
!
SUBROUTINE put_gradients(rln,rlt,vpar_shear,vexb_shear)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: rln(nsm),rlt(nsm),vpar_shear(nsm)
  REAL,INTENT(IN) :: vexb_shear
  INTEGER :: is

  do is=1,nsm
   if(tglf_isnan(rln(is)))call tglf_error(1,"input rlns_in is NAN")
   if(tglf_isinf(rln(is)))call tglf_error(1,"input rlns_in is INF")
   if(tglf_isnan(rlt(is)))call tglf_error(1,"input rlts_in is NAN")
   if(tglf_isnan(rlt(is)))call tglf_error(1,"input rlts_in is INF")
   if(tglf_isnan(vpar_shear(is)))call tglf_error(1,"input vpar_shear_in is NAN")
   if(tglf_isnan(vpar_shear(is)))call tglf_error(1,"input vpar_shear_in is INF")
  enddo
  if(tglf_isnan(vexb_shear))call tglf_error(1,"input vexb_shear_in is NAN")
  if(tglf_isinf(vexb_shear))call tglf_error(1,"input vexb_shear_in is INF")

  !
  ! transfer values
  !
  do is=1,nsm
     rlns_in(is) = rln(is)
     rlts_in(is) = rlt(is)
     vpar_shear_in(is) = vpar_shear(is)
  enddo
  vexb_shear_in = vexb_shear
  !    
END SUBROUTINE put_gradients
!
!
!-----------------------------------------------------------------
!
SUBROUTINE put_profile_shear(vns_shear,vts_shear)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: vns_shear(nsm),vts_shear(nsm)
  INTEGER :: is
  !
  ! transfer values
  !
  do is=1,nsm
     vns_shear_in(is) = vns_shear(is)
     vts_shear_in(is) = vts_shear(is)
  enddo
  !    
END SUBROUTINE put_profile_shear
!
!-----------------------------------------------------------------
!
SUBROUTINE put_averages(tsp,asp,vpar,vexb,betae,xnue,zeff,debye)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: tsp(nsm),asp(nsm),vpar(nsm)
  REAL,INTENT(IN) :: vexb,betae,xnue,zeff,debye
  INTEGER :: is, count

  do is=1,nsm
    if(tglf_isnan(tsp(is)))call tglf_error(1,"input taus_in is NAN")
    if(tglf_isinf(tsp(is)))call tglf_error(1,"input taus_in is INF")
    if(tsp(is).lt.0.0)call tglf_error(1,"input taus_in is < 0")
    if(tglf_isnan(asp(is)))call tglf_error(1,"input as_in is NAN")
    if(tglf_isinf(asp(is)))call tglf_error(1,"input as_in is INF")
    if(asp(is).lt.0.0)call tglf_error(1,"input as_in is < 0")
    if(tglf_isnan(vpar(is)))call tglf_error(1,"input vpar_in is NAN")
    if(tglf_isinf(vpar(is)))call tglf_error(1,"input vpar_in is INF")
  enddo
    if(tglf_isnan(vexb))call tglf_error(1,"input vexb_in is NAN")
    if(tglf_isinf(vexb))call tglf_error(1,"input vexb_in is INF")
    if(tglf_isnan(betae))call tglf_error(1,"input betae_in is NAN")
    if(tglf_isinf(betae))call tglf_error(1,"input betae_in is INF")
    if(betae.lt.0.0)call tglf_error(1,"input betae_in is < 0")
    if(tglf_isnan(xnue))call tglf_error(1,"input xnue_in is NAN")
    if(tglf_isinf(xnue))call tglf_error(1,"input xnue_in is INF")
    if(xnue.lt.0.0)call tglf_error(1,"input xnue_in is < 0")
    if(tglf_isnan(zeff))call tglf_error(1,"input zeff_in is NAN")
    if(tglf_isinf(zeff))call tglf_error(1,"input zeff_in is INF")
    if(zeff.lt.0.0)call tglf_error(1,"input zeff_in is < 0")
    if(tglf_isnan(debye))call tglf_error(1,"input debye_in is NAN")
    if(tglf_isinf(debye))call tglf_error(1,"input debye_in is INF")
    if(debye.lt.0.0)call tglf_error(1,"input debye_in is < 0")
  !
  ! set flow control switch
  new_matrix = .TRUE.
  ! transfer values
  count = 0
  do is=1,nsm
     taus_in(is) = tsp(is)
     as_in(is) = asp(is)
     vpar_in(is) = vpar(is)
      if(as_in(is) .gt. 0.0) count = count +1
  enddo
  nstotal_in = count
  if(as_in(1) .le. 0.0 .or. as_in(2) .le. 0.0)then
     write(*,*)"error: electron or ion density <= zero. Setting them both to 1.0"
      as_in(1) = 1.0
      as_in(2) = 1.0
  endif
  if(taus_in(1) .le. 0.0 .or. taus_in(2) .le. 0.0)then
     write(*,*)"error: electron or ion temperature <= zero. Setting them both to 1.0"
      taus_in(1) = 1.0
      taus_in(2) = 1.0
  endif
  vexb_in = vexb
  betae_in = betae
  xnue_in = xnue
  zeff_in = zeff
  debye_in = debye
  !      
END SUBROUTINE put_averages
!
!-----------------------------------------------------------------
!
SUBROUTINE put_switches(iflux,use_bper,use_bpar,use_mhd_rule,use_bisection, &
     use_inboard_detrapped,ibranch,nmodes,nb_max,nb_min,nxgrid,nkys,use_ave_ion_grid)
  !
  USE tglf_global
  USE tglf_dimensions
  !
  IMPLICIT NONE
  LOGICAL :: iflux,use_bper,use_bpar,use_mhd_rule,use_bisection
  LOGICAL :: use_inboard_detrapped, use_ave_ion_grid
  INTEGER :: ibranch,nmodes,nb_max,nb_min
  INTEGER :: nxgrid,nkys
  !
  ! validaty checks
  ! reset to defaults if invlaid
  !
  if(nb_max.lt.2.or.nb_max.gt.nb)nb_max=nbasis_max_in
  if(nb_min.lt.2.or.nb_min.gt.nb)nb_min=nbasis_min_in
  if(nb_max.lt.nb_min)nb_max=nb_min
  if(2*(nb_max/2).ne.nb_max)nb_max = 2*(nb_max/2)  ! must be even
  if(2*(nb_min/2).ne.nb_min)nb_min = 2*(nb_min/2)  ! must be even
  if(ibranch.lt.-1.or.ibranch.gt.0)ibranch=ibranch_in
  if(nxgrid.lt.1.or.2*nxgrid-1.gt.nxm)nxgrid=MIN((nxm+1)/2,nxgrid_in)
  if(nmodes.lt.1.or.nmodes.gt.maxmodes)nmodes=nmodes_in
  if(nkys.lt.2.or.nkys.gt.nkym)nkys=nky_in
  !      write(*,*)nb_max,nb_min,ibranch,nxgrid,nmodes,nkys
  !
  ! check for changes and update flow controls
  !
  !      if(nxgrid.ne.nxgrid_in)new_start = .TRUE.
  !      if(nb_max.ne.nbasis_max_in)new_start = .TRUE.
  !
  ! transfer values
  !
  iflux_in = iflux
  use_bper_in = use_bper
  use_bpar_in = use_bpar
  use_mhd_rule_in = use_mhd_rule
  use_bisection_in = use_bisection
  ibranch_in = ibranch
  nmodes_in = nmodes
  if(ibranch_in.eq.0)nmodes_in=2
  nbasis_max_in = nb_max
  nbasis_min_in = nb_min
  nxgrid_in = nxgrid
  nky_in = nkys
  use_inboard_detrapped_in = use_inboard_detrapped
  use_ave_ion_grid_in = use_ave_ion_grid
  !
END SUBROUTINE put_switches
!
!-----------------------------------------------------------------
!
SUBROUTINE put_rare_switches(rtheta_trap,rwdia_trap,rpark,rghat,rgchat, &
     rwd_zero,rLinsker,rgradB,rfilter,rdamp_psi,rdamp_sig)
  !
  USE tglf_global
  USE tglf_dimensions
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: rtheta_trap,rwdia_trap,rpark,rghat,rgchat,rwd_zero
  REAL,INTENT(IN) :: rLinsker,rgradB,rfilter
  REAL,INTENT(IN) :: rdamp_psi,rdamp_sig
  !
  ! transfer values
  !
  theta_trapped_in = rtheta_trap
  wdia_trapped_in = rwdia_trap
  park_in = rpark
  ghat_in = rghat
  gchat_in = rgchat
  wd_zero_in = rwd_zero
  Linsker_factor_in = rLinsker
  gradB_factor_in = rgradB
  filter_in = rfilter
  damp_psi_in = rdamp_psi
  damp_sig_in = rdamp_sig
  !
END SUBROUTINE put_rare_switches
!
!-----------------------------------------------------------------
!
SUBROUTINE put_model_parameters(adi_elec,alpha_e,alpha_p,alpha_mach,  &
     alpha_quench,alpha_zf,xnu_fac,debye_fac,etg_fac,rlnp_cut,        &
     sat_rule,kygrid_model,xnu_model,vpar_model,vpar_shear_model)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  LOGICAL,INTENT(IN) :: adi_elec
  INTEGER :: sat_rule,xnu_model,kygrid_model
  INTEGER :: vpar_model,vpar_shear_model
  REAL,INTENT(IN) :: alpha_e,alpha_p,alpha_mach,alpha_zf
  REAL,INTENT(IN) :: alpha_quench,xnu_fac,debye_fac,etg_fac,rlnp_cut
  !
  ! check for changes and update flow controls
  !
  if(adi_elec .NEQV. adiabatic_elec_in)new_matrix = .TRUE.
  if(kygrid_model.lt.0.or.kygrid_model.gt.5)kygrid_model = kygrid_model_in
  if(xnu_model.lt.0.or.xnu_model.gt.3)xnu_model = xnu_model_in
  if(sat_rule.lt.0.or.sat_rule.gt.3)sat_rule=sat_rule_in
  !if(vpar_model.lt.-1.or.vpar_model.gt.1)vpar_model=vpar_model_in
  if(vpar_shear_model.lt.0.or.vpar_shear_model.gt.1)vpar_shear_model=vpar_shear_model_in
  !
  ! transfer values
  !
  adiabatic_elec_in = adi_elec
  alpha_mach_in = alpha_mach
  alpha_p_in = alpha_p
  alpha_e_in = alpha_e
  alpha_quench_in = alpha_quench
  alpha_zf_in = alpha_zf
  xnu_factor_in = xnu_fac
  debye_factor_in = debye_fac
  etg_factor_in = etg_fac
  rlnp_cutoff_in = rlnp_cut
  sat_rule_in = sat_rule
  xnu_model_in = xnu_model
  kygrid_model_in = kygrid_model
  !vpar_model_in = vpar_model      depreciated input switch
  vpar_shear_model_in = vpar_shear_model
  !
  if(alpha_quench_in .ne.0.0)then
     ! turn off spectral shift and model and only use Waltz quench rule
     alpha_e_in = 0.0
  endif
  !
END SUBROUTINE put_model_parameters
!
!-----------------------------------------------------------------
!
SUBROUTINE put_s_alpha_geometry(rmin,rmaj,q,shat,alpha,xwell, &
     theta0,b_model,ft_model)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER:: b_model,ft_model
  REAL,INTENT(IN) :: rmin,rmaj,q,shat,alpha,theta0,xwell

  if(tglf_isnan(rmin))call tglf_error(1,"input rmin_sa is NAN")
  if(tglf_isinf(rmin))call tglf_error(1,"input rmin_sa is INF")
  if(tglf_isnan(rmaj))call tglf_error(1,"input rmaj_sa is NAN")
  if(tglf_isinf(rmaj))call tglf_error(1,"input rmaj_sa is INF")
  if(tglf_isnan(q))call tglf_error(1,"input q_sa is NAN")
  if(tglf_isinf(q))call tglf_error(1,"input q_sa is INF")
  if(tglf_isnan(shat))call tglf_error(1,"input shat_sa is NAN")
  if(tglf_isinf(shat))call tglf_error(1,"input shat_sa is INF")
  if(tglf_isnan(alpha))call tglf_error(1,"input alpha_sa is NAN")
  if(tglf_isinf(alpha))call tglf_error(1,"input alpha_sa is INF")
  if(tglf_isnan(xwell))call tglf_error(1,"input xwell_sa is NAN")
  if(tglf_isinf(xwell))call tglf_error(1,"input xwell_sa is INF")
  if(tglf_isnan(theta0))call tglf_error(1,"input theta0_sa is NAN")
  if(tglf_isinf(theta0))call tglf_error(1,"input theta0_sa is INF")
  !
  ! set geometry type flag for shifted circle
  igeo = 0
  ! set flow control switch
  new_geometry = .TRUE.
  !
  ! transfer values
  !
  rmin_sa = rmin
  rmaj_sa = rmaj
  q_sa = ABS(q)
  q_in = q_sa   ! needed for kygrid_model_in = 3
  shat_sa = shat
  alpha_sa = alpha
  xwell_sa = xwell
  theta0_sa = theta0
  b_model_sa = b_model
  ft_model_sa = ft_model
  !
  ! validatiy checks
  !
  if(rmin_sa.ge.rmaj_sa)rmin_sa=0.999*rmaj_sa   
  if(ft_model_sa.lt.0.or.ft_model_sa.gt.3)then
     write(*,*)"******* ERROR ft_model_sa invalid *******"  
     ft_model_sa = 1
  endif
  !
END SUBROUTINE put_s_alpha_geometry
!
!-----------------------------------------------------------------
!

SUBROUTINE put_Miller_geometry(rmin,rmaj,zmaj,drmindx,drmajdx,dzmajdx, &
  kappa,s_kappa,delta,s_delta,zeta,s_zeta,q,q_prime,p_prime,beta,kx0_m)
  !
  ! This routine eliminates the need for subroutine miller_init 
  ! and the miller.dat input file.
  !
  USE tglf_global
  !
  IMPLICIT NONE
  REAL,INTENT(IN) :: rmin,rmaj,zmaj,q,q_prime,p_prime,kx0_m,beta
  REAL,INTENT(IN) :: drmindx,drmajdx,dzmajdx
  REAL,INTENT(IN) :: kappa,s_kappa,delta,s_delta,zeta,s_zeta

  if(tglf_isnan(rmin))call tglf_error(1,"input rmin_loc is NAN")
  if(tglf_isinf(rmin))call tglf_error(1,"input rmin_loc is INF")
  if(tglf_isnan(rmaj))call tglf_error(1,"input rmaj_loc is NAN")
  if(tglf_isinf(rmaj))call tglf_error(1,"input rmaj_loc is INF")
  if(tglf_isnan(drmindx))call tglf_error(1,"input drmindx_loc is NAN")
  if(tglf_isinf(drmindx))call tglf_error(1,"input drmindx_loc is INF")
  if(tglf_isnan(drmajdx))call tglf_error(1,"input drmajdx_loc is NAN")
  if(tglf_isinf(drmajdx))call tglf_error(1,"input drmajdx_loc is INF")
  if(tglf_isnan(dzmajdx))call tglf_error(1,"input dzmajdx_loc is NAN")
  if(tglf_isinf(dzmajdx))call tglf_error(1,"input dzmajdx_loc is INF")
  if(tglf_isnan(kappa))call tglf_error(1,"input kappa_loc is NAN")
  if(tglf_isinf(kappa))call tglf_error(1,"input kappa_loc is INF")
  if(tglf_isnan(s_kappa))call tglf_error(1,"input s_kappa_loc is NAN")
  if(tglf_isinf(s_kappa))call tglf_error(1,"input s_kappa_loc is INF")
  if(tglf_isnan(delta))call tglf_error(1,"input delta_loc is NAN")
  if(tglf_isinf(delta))call tglf_error(1,"input delta_loc is INF")
  if(tglf_isnan(s_delta))call tglf_error(1,"input s_delta_loc is NAN")
  if(tglf_isinf(s_delta))call tglf_error(1,"input s_delta_loc is INF")
  if(tglf_isnan(zeta))call tglf_error(1,"input zeta_loc is NAN")
  if(tglf_isinf(zeta))call tglf_error(1,"input zeta_loc is INF")
  if(tglf_isnan(s_zeta))call tglf_error(1,"input s_zeta_loc is NAN")
  if(tglf_isinf(s_zeta))call tglf_error(1,"input s_zeta_loc is INF")
  if(tglf_isnan(q))call tglf_error(1,"input q_loc is NAN")
  if(tglf_isinf(q))call tglf_error(1,"input q_loc is INF")
  if(tglf_isnan(q_prime))call tglf_error(1,"input q_prime_loc is NAN")
  if(tglf_isinf(q_prime))call tglf_error(1,"input q_prime_loc is INF")
  if(tglf_isnan(p_prime))call tglf_error(1,"input p_prime_loc is NAN")
  if(tglf_isinf(p_prime))call tglf_error(1,"input p_prime_loc is INF")
  if(tglf_isnan(beta))call tglf_error(1,"input beta_loc is NAN")
  if(tglf_isinf(beta))call tglf_error(1,"input beta_loc is INF")
  if(tglf_isnan(kx0_m))call tglf_error(1,"input kx0_loc is NAN")
  if(tglf_isinf(kx0_m))call tglf_error(1,"input kx0_loc is INF")

  !
  ! set geometry type flag for Miller
  !
  igeo = 1
  !
  ! set flow control switch
  !
  new_geometry = .TRUE.
  !
  ! transfer values
  !
  rmin_loc = rmin
  rmaj_loc = rmaj
  zmaj_loc = zmaj
  q_loc = ABS(q)
  q_in = q_loc   ! needed for kygrid_model_in = 3
  p_prime_loc = p_prime
  q_prime_loc = q_prime
  drmindx_loc = drmindx
  drmajdx_loc = drmajdx
  dzmajdx_loc = dzmajdx
  kappa_loc = kappa
  s_kappa_loc = s_kappa
  delta_loc = delta
  s_delta_loc = s_delta
  zeta_loc = zeta
  s_zeta_loc = s_zeta
  beta_loc = beta
  kx0_loc = kx0_m
  !
  ! validatiy checks
  !
  if(rmin_loc.ge.rmaj_loc)rmin_loc = 0.999*rmaj_loc
  if(drmajdx_loc.le.-1.0)drmajdx_loc=-0.999
  !
END SUBROUTINE put_Miller_geometry
!
!-----------------------------------------------------------------
!
SUBROUTINE put_Fourier_geometry(q,q_prime,p_prime,nf,f)
  !
  ! This routine transfers the inputs for the fourier_geo model 
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nf
  REAL,INTENT(IN) :: q,q_prime,p_prime
  REAL,INTENT(IN) :: f(8,0:max_fourier)
  INTEGER :: i,j
  !
  ! validatiy checks
  !
  if(nf.gt.max_fourier)call tglf_error(1,"input nfourier_in exceeds maximum")
  if(tglf_isnan(q))call tglf_error(1,"input q_fourier_in is NAN")
  if(tglf_isinf(q))call tglf_error(1,"input q_fourier_in is INF")
  if(tglf_isnan(q_prime))call tglf_error(1,"input q_prime_fourier_in is NAN")
  if(tglf_isinf(q_prime))call tglf_error(1,"input q_prime_fourier_in is INF")
  if(tglf_isnan(p_prime))call tglf_error(1,"input p_prime_fourier_in is NAN")
  if(tglf_isinf(p_prime))call tglf_error(1,"input p_prime_fourier_in is INF")
  do i=1,nf
    do j=1,8
      if(tglf_isnan(f(j,i)))call tglf_error(1,"input fourrier_in is NAN")
      if(tglf_isinf(f(j,i)))call tglf_error(1,"input fourrier_in is INF")
    enddo
  enddo
  !
  ! set geometry type flag for Fourier
  !
  igeo = 2
  !
  ! set flow control switch
  !
  new_geometry = .TRUE.
  !
  ! transfer values
  !
  q_fourier_in = ABS(q)
  q_in = q_loc   ! needed for kygrid_model_in = 3
  p_prime_fourier_in = p_prime
  q_prime_fourier_in = q_prime
  nfourier_in = nf
  fourier_in(:,:)=f(:,:)
END SUBROUTINE put_Fourier_geometry
!
!-----------------------------------------------------------------
!
SUBROUTINE put_ELITE_geometry(nc,q,q_prime,p_prime,r_c,z_c,bp_c)
  !
  ! This routine requires having read the data from an ELITE geometry file
  ! giving R,Z,Bp on a flux surface contour with nc points.
  !
  USE tglf_global
  USE tglf_sgrid
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nc
  REAL,INTENT(IN) :: q,q_prime,p_prime
  REAL,INTENT(IN) :: r_c(0:max_ELITE),z_c(0:max_ELITE),bp_c(0:max_ELITE)
  !
  INTEGER i
  !
  ! set geometry type flag for ELITE
  !
  igeo = 2
  !
  ! set flow control switch
  !
  new_geometry = .TRUE.
  !
  ! validatiy checks
  !
  if(nc.gt.max_ELITE)call tglf_error(1,"input n_ELITE exceeds limit")
  if(nc.lt.ms)call tglf_error(1,"input n_ELITE is too small")
  if(tglf_isnan(q_prime))call tglf_error(1,"input q_ELITE is NAN")
  if(tglf_isinf(q_prime))call tglf_error(1,"input q_ELITE is INF")
  do i=0,nc
    if(tglf_isnan(r_c(i)))call tglf_error(1,"input R_ELITE is NAN")
    if(tglf_isinf(r_c(i)))call tglf_error(1,"input R_ELITE is INF")
    if(tglf_isnan(z_c(i)))call tglf_error(1,"input Z_ELITE is NAN")
    if(tglf_isinf(z_c(i)))call tglf_error(1,"input Z_ELITE is INF")
    if(tglf_isnan(bp_c(i)))call tglf_error(1,"input Bp_ELITE is NAN")
    if(tglf_isinf(bp_c(i)))call tglf_error(1,"input Bp_ELITE is INF")
  enddo
  !
  ! transfer values
  !
  n_ELITE = nc
  q_ELITE = q
  q_in = q   ! needed for kygrid=3
  q_prime_ELITE = q_prime
  p_prime_ELITE = p_prime
  ! note direction change for contour angle
  do i=0,n_ELITE
     R_ELITE(n_ELITE-i) = r_c(i)
     Z_ELITE(n_ELITE-i) = z_c(i)
     Bp_ELITE(n_ELITE-i) = bp_c(i)
  enddo
  !
END SUBROUTINE put_ELITE_geometry
!
!-----------------------------------------------------------------
!  output routines
!-----------------------------------------------------------------
!
REAL FUNCTION get_growthrate(index1)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: index1
  INTEGER :: i3
  !
  i3 = SIZE(gamma_out)
  if(index1.gt.i3)then
     write(*,*)"requested growthrate index out of bounds",i3
     get_growthrate = 0.0
  else
     get_growthrate = gamma_out(index1)
  endif
  !
END FUNCTION get_growthrate
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_frequency(index1)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: index1
  INTEGER :: i3
  !
  i3 = SIZE(freq_out)
  if(index1.gt.i3)then
     write(*,*)"requested frequency index is of bounds",i3
     get_frequency = 0.0
  else
     get_frequency = freq_out(index1)
  endif
  !
END FUNCTION get_frequency
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_QL_particle_flux(i1,i2,i3)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2,i3
  INTEGER :: i4
  !
  i4=SIZE(particle_QL_out)
  if(i1*i2*i3.gt.i4)then
     write(*,*)"requested QL particle flux index is of bounds",i4
     get_QL_particle_flux = 0.0
  else
     get_QL_particle_flux = particle_QL_out(i1,i2,i3)
  endif
  !
END FUNCTION get_QL_particle_flux
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_QL_energy_flux(i1,i2,i3)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2,i3
  INTEGER :: i4
  !
  i4=SIZE(energy_QL_out)
  if(i1*i2*i3.gt.i4)then
     write(*,*)"requested QL energy flux index is of bounds",i3
     get_QL_energy_flux = 0.0
  else
     get_QL_energy_flux = energy_QL_out(i1,i2,i3)
  endif
  !
END FUNCTION get_QL_energy_flux
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_QL_stress_par(i1,i2,i3)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2,i3
  INTEGER :: i4
  !
  i4=SIZE(stress_par_QL_out)
  if(i1*i2*i3.gt.i4)then
     write(*,*)"requested QL stress_par index is of bounds",i3
     get_QL_stress_par = 0.0
  else
     get_QL_stress_par = stress_par_QL_out(i1,i2,i3)
  endif
  !
END FUNCTION get_QL_stress_par
!
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_QL_stress_tor(i1,i2,i3)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2,i3
  INTEGER :: i4
  !
  i4=SIZE(stress_tor_QL_out)
  if(i1*i2*i3.gt.i4)then
     write(*,*)"requested QL stress_tor index is of bounds",i3
     get_QL_stress_tor = 0.0
  else
     get_QL_stress_tor = stress_tor_QL_out(i1,i2,i3)
  endif
  !
END FUNCTION get_QL_stress_tor
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_QL_exchange(i1,i2,i3)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2,i3
  INTEGER :: i4
  !
  i4=SIZE(exchange_QL_out)
  if(i1*i2*i3.gt.i4)then
     write(*,*)"requested QL exchange index is of bounds",i3
     get_QL_exchange = 0.0
  else
     get_QL_exchange = exchange_QL_out(i1,i2,i3)
  endif
  !
END FUNCTION get_QL_exchange
!-----------------------------------------------------------------
!
REAL FUNCTION get_QL_density(i1,i2)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  i3=SIZE(N_QL_out)
  if(i1*i2.gt.i3)then
     write(*,*)"requested QL_density index is of bounds",i3
     get_QL_density = 0.0
  else
     get_QL_density = N_QL_out(i1,i2)
  endif
  !
END FUNCTION get_QL_density
!-----------------------------------------------------------------
!
REAL FUNCTION get_QL_temperature(i1,i2)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  i3=SIZE(T_QL_out)
  if(i1*i2.gt.i3)then
     write(*,*)"requested QL_temperature index is of bounds",i3
     get_QL_temperature = 0.0
  else
     get_QL_temperature = T_QL_out(i1,i2)
  endif
  !
END FUNCTION get_QL_temperature
!-----------------------------------------------------------------
!
REAL FUNCTION get_ne_te_phase(i1)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1
  INTEGER :: i3
  !
  i3=SIZE(ne_te_phase_out)
  if(i1.gt.i3)then
     write(*,*)"requested ne_te_phase index is of bounds",i3
     get_ne_te_phase = 0.0
  else
     get_ne_te_phase = ne_te_phase_out(i1)
  endif
  !
END FUNCTION get_ne_te_phase
!-----------------------------------------------------------------
!
REAL FUNCTION get_phi_bar(i1)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1
  INTEGER :: i3
  !
  i3=SIZE(phi_bar_out)
  if(i1.gt.i3)then
     write(*,*)"requested phi_bar index is of bounds",i3
     get_phi_bar = 0.0
  else
     get_phi_bar = phi_bar_out(i1)
  endif
  !
END FUNCTION get_phi_bar
!-----------------------------------------------------------------
!
REAL FUNCTION get_kpar_bar(i1)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1
  INTEGER :: i3
  !
  i3=SIZE(kpar_bar_out)
  if(i1.gt.i3)then
     write(*,*)"requested kpar_bar index is of bounds",i3
     get_kpar_bar = 0.0
  else
     get_kpar_bar = kpar_bar_out(i1)
  endif
  !
END FUNCTION get_kpar_bar
!-----------------------------------------------------------------
!
REAL FUNCTION get_v_bar(i1)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1
  INTEGER :: i3
  !
  i3=SIZE(v_bar_out)
  if(i1.gt.i3)then
     write(*,*)"requested v_bar index is of bounds",i3
     get_v_bar = 0.0
  else
     get_v_bar = v_bar_out(i1)
  endif
  !
END FUNCTION get_v_bar
!-----------------------------------------------------------------
!
REAL FUNCTION get_n_bar(i1,i2)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  i3=SIZE(n_bar_out)
  if(i1*i2.gt.i3)then
     write(*,*)"requested n_bar index is of bounds",i3
     get_n_bar = 0.0
  else
     get_n_bar = n_bar_out(i1,i2)
  endif
  !
END FUNCTION get_n_bar
!-----------------------------------------------------------------
!
REAL FUNCTION get_t_bar(i1,i2)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  i3=SIZE(t_bar_out)
  if(i1*i2.gt.i3)then
     write(*,*)"requested t_bar index is of bounds",i3
     get_t_bar = 0.0
  else
     get_t_bar = t_bar_out(i1,i2)
  endif
  !
END FUNCTION get_t_bar
!-----------------------------------------------------------------
!
REAL FUNCTION get_gaussian_width()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_gaussian_width = width_in
  !
END FUNCTION get_gaussian_width
!-----------------------------------------------------------------
!
REAL FUNCTION get_R_unit()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_R_unit = R_unit
  !
END FUNCTION get_R_unit
!-----------------------------------------------------------------
!
REAL FUNCTION get_ft()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_ft = ft
  !
END FUNCTION get_ft
!-----------------------------------------------------------------
REAL FUNCTION get_B_unit()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_B_unit = B_unit
  !
END FUNCTION get_B_unit
!-----------------------------------------------------------------
!
REAL FUNCTION get_q_unit()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_q_unit = q_unit
  !
END FUNCTION get_q_unit
!-----------------------------------------------------------------
!
REAL FUNCTION get_B2_ave()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_B2_ave = B2_ave_out
  !
END FUNCTION get_B2_ave
!-----------------------------------------------------------------
!
REAL FUNCTION get_R2_ave()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_R2_ave = R2_ave_out
  !
END FUNCTION get_R2_ave
!-----------------------------------------------------------------
!
REAL FUNCTION get_B_ave()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_B_ave = B_ave_out
  !
END FUNCTION get_B_ave
!-----------------------------------------------------------------
!
REAL FUNCTION get_Bt_ave()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_Bt_ave = Bt_ave_out
  !
END FUNCTION get_Bt_ave
!-----------------------------------------------------------------
!
REAL FUNCTION get_Bp0()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_Bp0 = Bp0_out
  !
END FUNCTION get_Bp0
!-----------------------------------------------------------------
!
REAL FUNCTION get_RBt_ave()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_RBt_ave = RBt_ave_out
  !
END FUNCTION get_RBt_ave
!-----------------------------------------------------------------
!
REAL FUNCTION get_b0_bar(n1)
  !
  USE tglf_global
  USE tglf_coeff
  !
  IMPLICIT NONE
  ! 
  INTEGER,INTENT(IN):: n1
  !
  get_b0_bar = b0_bar_out(n1)
  !
END FUNCTION get_b0_bar
!-----------------------------------------------------------------
!
REAL FUNCTION get_wd_bar(n1)
  !
  USE tglf_global
  USE tglf_coeff
  !
  IMPLICIT NONE
  ! 
  INTEGER,INTENT(IN):: n1
  !
  get_wd_bar = wd_bar_out(n1)
  !
END FUNCTION get_wd_bar
!-----------------------------------------------------------------
!
REAL FUNCTION get_particle_flux(i1,i2)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  i3=SIZE(particle_flux_out)
  if(i1*i2.gt.i3)then
     write(*,*)"requested particle flux index is of bounds",i3
     get_particle_flux = 0.0
  else
     get_particle_flux = particle_flux_out(i1,i2)
  endif
  !
END FUNCTION get_particle_flux
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_energy_flux(i1,i2)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  i3=SIZE(energy_flux_out)
  if(i1*i2.gt.i3)then
     write(*,*)"requested energy flux index is of bounds",i3
     get_energy_flux = 0.0
  else
     get_energy_flux = energy_flux_out(i1,i2)
  endif
  !
END FUNCTION get_energy_flux
!-----------------------------------------------------------------
!
REAL FUNCTION get_stress_par(i1,i2)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  i3=SIZE(stress_par_out)
  if(i1*i2.gt.i3)then
     write(*,*)"requested stress_par index is of bounds",i3
     get_stress_par = 0.0
  else
     get_stress_par = stress_par_out(i1,i2)
  endif
  !
END FUNCTION get_stress_par
!!-----------------------------------------------------------------
!
REAL FUNCTION get_stress_tor(i1,i2)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  i3=SIZE(stress_tor_out)
  if(i1*i2.gt.i3)then
     write(*,*)"requested stress_tor index is of bounds",i3
     get_stress_tor = 0.0
  else
     get_stress_tor = stress_tor_out(i1,i2)
  endif
  !
END FUNCTION get_stress_tor
!-----------------------------------------------------------------
!
REAL FUNCTION get_exchange(i1,i2)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1,i2
  INTEGER :: i3
  !
  i3=SIZE(exchange_out)
  if(i1*i2.gt.i3)then
     write(*,*)"requested exchange index is of bounds",i3
     get_exchange = 0.0
  else
     get_exchange = exchange_out(i1,i2)
  endif
  !
END FUNCTION get_exchange
!
!-----------------------------------------------------------------
!
REAL FUNCTION get_n_bar_sum(i1)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1
  INTEGER :: i3
  !
  i3=SIZE(n_bar_sum_out)
  if(i1.gt.i3)then
     write(*,*)"requested n_bar_sum index is of bounds",i3
     get_n_bar_sum = 0.0
  else
     get_n_bar_sum = n_bar_sum_out(i1)
  endif
  !
END FUNCTION get_n_bar_sum
!-----------------------------------------------------------------
!
REAL FUNCTION get_q_low(i1)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1
  INTEGER :: i3
  !
  i3=SIZE(q_low_out)
  if(i1.gt.i3)then
     write(*,*)"requested energy flux index is of bounds",i3
     get_q_low = 0.0
  else
     get_q_low = q_low_out(i1)
  endif
  !
END FUNCTION get_q_low
!-----------------------------------------------------------------
!
REAL FUNCTION get_t_bar_sum(i1)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i1
  INTEGER :: i3
  !
  i3=SIZE(t_bar_sum_out)
  if(i1.gt.i3)then
     write(*,*)"requested t_bar_sum index is of bounds",i3
     get_t_bar_sum = 0.0
  else
     get_t_bar_sum = t_bar_sum_out(i1)
  endif
  !
END FUNCTION get_t_bar_sum
!-----------------------------------------------------------------
!
REAL FUNCTION get_phi_bar_sum()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_phi_bar_sum = phi_bar_sum_out
  !
END FUNCTION get_phi_bar_sum
!-----------------------------------------------------------------
!
REAL FUNCTION get_v_bar_sum()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_v_bar_sum = v_bar_sum_out
  !
END FUNCTION get_v_bar_sum
!-----------------------------------------------------------------
!
INTEGER FUNCTION get_nky_out()
  !
  USE tglf_global
  USE tglf_kyspectrum
  !
  IMPLICIT NONE
  !
  get_nky_out = nky
  !
END FUNCTION get_nky_out
!----------------------------------------------------------------
!
REAL FUNCTION get_flux_spectrum_out(itype,ispec,ifield,iky,imode)
  !
  USE tglf_global
  !   
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: itype,ispec,ifield,iky,imode
  INTEGER :: error
  !
  error=0
  get_flux_spectrum_out = 0.0
  if(itype.lt.1.or.itype.gt.5)then
     write(*,*)"itype out of bounds",1,5
     error=1
  elseif(ispec.lt.1.or.ispec.gt.nsm)then
     write(*,*)"ispec out of bounds",1,nsm
     error=1        
  elseif(ifield.lt.1.or.ifield.gt.3)then
     write(*,*)"ifield out of bounds",1,3
     error=1        
  elseif(iky.lt.1.or.iky.gt.nkym)then
     write(*,*)"iky out of bounds",1,nkym
     error=1
  elseif(imode.lt.1.or.imode.gt.maxmodes)then
     write(*,*)"imode out of bounds",1,maxmodes
     error=1
  endif
  !
  if(error.eq.0)then
     get_flux_spectrum_out=flux_spectrum_out(itype,ispec,ifield,iky,imode)
  endif
  !
END FUNCTION get_flux_spectrum_out
!----------------------------------------------------------------
!
REAL FUNCTION get_eigenvalue_spectrum_out(itype,iky,imode)
  !
  USE tglf_global
  !   
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: itype,iky,imode
  INTEGER :: error
  !
  error=0
  get_eigenvalue_spectrum_out = 0.0
  if(itype.lt.1.or.itype.gt.2)then
     write(*,*)"ntype out of bounds",1,2
     error=1
  elseif(iky.lt.1.or.iky.gt.nkym)then
     write(*,*)"iky out of bounds",1,nkym
     error=1
  elseif(imode.lt.1.or.imode.gt.maxmodes)then
     write(*,*)"imode out of bounds",1,maxmodes
     error=1
  endif
  !
  if(error.eq.0)then
     get_eigenvalue_spectrum_out=eigenvalue_spectrum_out(itype,iky,imode)
  endif
  !
END FUNCTION get_eigenvalue_spectrum_out
!----------------------------------------------------------------
!
REAL FUNCTION get_intensity_spectrum_out(itype,ispec,iky,imode)
  !
  USE tglf_global
  !   
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: itype,ispec,iky,imode
  INTEGER :: error
  !
  error=0
  get_intensity_spectrum_out = 0.0
  if(itype.lt.1.or.itype.gt.4)then
     write(*,*)"ntype out of bounds",1,4
     error=1
  elseif(ispec.lt.1.or.ispec.gt.nsm)then
     write(*,*)"ispec out of bounds",1,nsm
     error=1        
  elseif(iky.lt.1.or.iky.gt.nkym)then
     write(*,*)"iky out of bounds",1,nkym
     error=1
  elseif(imode.lt.1.or.imode.gt.maxmodes)then
     write(*,*)"imode out of bounds",1,maxmodes
     error=1
  endif
  !
  if(error.eq.0)then
     get_intensity_spectrum_out=intensity_spectrum_out(itype,ispec,iky,imode)
  endif
  !
END FUNCTION get_intensity_spectrum_out
!----------------------------------------------------------------
!
REAL FUNCTION get_field_spectrum_out(itype,iky,imode)
  !
  USE tglf_global
  !   
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: itype,iky,imode
  INTEGER :: error
  !
  error = 0
  get_field_spectrum_out = 0.0
  if(itype.lt.1.or.itype.gt.4)then
     write(*,*)"itype out of bounds",1,2
     error=1
  elseif(iky.lt.1.or.iky.gt.nkym)then
     write(*,*)"nky out of bounds",1,nkym
     error=1
  elseif(imode.lt.1.or.imode.gt.maxmodes)then
     write(*,*)"imode out of bounds",1,maxmodes
     error=1
  endif
  !
  if(error.eq.0)then
     get_field_spectrum_out=field_spectrum_out(itype,iky,imode)
  endif
  !
END FUNCTION get_field_spectrum_out
!-----------------------------------------------------------------
!
REAL FUNCTION get_ky_spectrum_out(iky)
  !
  USE tglf_global
  USE tglf_kyspectrum
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: iky
  !
  if(iky.lt.1.or.iky.gt.nkym)then
     write(*,*)"iky out of bounds",1,nkym
     get_ky_spectrum_out = 0.0
  else
     get_ky_spectrum_out = ky_spectrum(iky)
  endif
  !
END FUNCTION get_ky_spectrum_out
!-----------------------------------------------------------------
!

REAL FUNCTION get_DM()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_DM = DM_out
  !
END FUNCTION get_DM
!-----------------------------------------------------------------
!
REAL FUNCTION get_DR()
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  get_DR = DR_out
  !
END FUNCTION get_DR
!-----------------------------------------------------------------
!
SUBROUTINE get_DEP_parameters(r_dep,rmaj_dep,q_dep,taui_dep,rlni_dep,rlti_dep,ni_dep)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  !
  REAL,INTENT(OUT) :: r_dep,rmaj_dep,q_dep,taui_dep,rlni_dep,rlti_dep,ni_dep
  !
  ! warning this routine assumes that the call put_species set the main ion species to be index 2
  !
  taui_dep = taus_in(2)
  rlni_dep = rlns_in(2)
  rlti_dep = rlts_in(2)
  ni_dep = as_in(2)
  ! note that rmin_input,Rmaj_input,q_input are set the different geometry routines 
  r_dep = rmin_input
  rmaj_dep = Rmaj_input
  q_dep = q_in
  !
END SUBROUTINE get_DEP_parameters
!-----------------------------------------------------------------
!
SUBROUTINE get_wavefunction_out(nmodes,nfields,nplot,angle,wavefunction)
  !
  USE tglf_global
  !
  IMPLICIT NONE
  INTEGER,INTENT(OUT) :: nmodes,nfields,nplot
  REAL,INTENT(OUT) :: angle(max_plot)
  COMPLEX,INTENT(OUT) :: wavefunction(maxmodes,3,max_plot)
  !
  if(new_start)then
     write(*,*)"error: tglf must be called before get_wavefunction_out"
  else
     call get_wavefunction
     nmodes = nmodes_out
     nfields = nfields_out
     nplot = max_plot
     angle(:) = plot_angle_out(:)
     wavefunction(:,:,:) = plot_field_out(:,:,:)
  endif
  !
END SUBROUTINE get_wavefunction_out
!-----------------------------------------------------------------



SUBROUTINE get_error_status(message,status)
!-------------------------------------------------------------------
! Determine internal TGLF error status and return string and integer
!
! TGLF success:
!  status=0
!  message='null'
!
! TGLF failure
!  status=1
!  message=<error description>
!
! TGLF warning
!  status=2
!  message=<warning description>
!-------------------------------------------------------------------

  use tglf_global

  implicit none

  character (len=*), intent(inout) :: message
  integer, intent(inout) :: status

  message = trim(error_msg)

  if (message == 'null') then
     status = 0
  else
     status = 1
  endif

END SUBROUTINE get_error_status

!
!-----------------------------------------------------------------
!
SUBROUTINE write_tglf_input(tglf_path_in)
  !
  USE tglf_pkg
  USE tglf_global
  !
  IMPLICIT NONE
  !
  INTEGER :: is
  character (len=256)  :: tglf_path_in
  !
  open(unit=11,file=trim(tglf_path_in)//"tglf_input",status="unknown")
  write(11,*)"&tglfin  ! version 1.93"
  write(11,*)"      adiabatic_elec_tg=",adiabatic_elec_in
  write(11,*)"      kygrid_model_tg = ",kygrid_model_in
  write(11,*)"      xnu_model_tg= ",xnu_model_in
  write(11,*)"      sat_rule_tg= ",sat_rule_in
  write(11,*)"      vpar_model_tg= ",vpar_model_in
  write(11,*)"      vpar_shear_model_tg= ",vpar_shear_model_in
  write(11,*)"      nbasis_max_tg= ",nbasis_max_in
  write(11,*)"      nbasis_min_tg= ",nbasis_min_in
  write(11,*)"      nxgrid_tg= ",nxgrid_in
  write(11,*)"      ibranch_tg= ",ibranch_in
  write(11,*)"      nmodes_tg= ",nmodes_in
  write(11,*)"      ns_tg= ",ns_in
  write(11,*)"      nky_tg = ",nky_in
  write(11,*)"      ky_tg= ",ky_in
  write(11,*)"      sign_Bt_tg= ",sign_Bt_in
  write(11,*)"      iflux_tg=",iflux_in
  write(11,*)"      igeo_tg= ",igeo
  write(11,*)"      width_max_tg= ",width_in
  write(11,*)"      width_min_tg= ",width_min_in
  write(11,*)"      find_width_tg=",find_width_in
  write(11,*)"      nwidth_tg= ",nwidth_in
  write(11,*)"      park_tg= ",park_in
  write(11,*)"      ghat_tg= ",ghat_in
  write(11,*)"      gchat_tg= ",gchat_in
  if(igeo.eq.0)then
     write(11,*)"      rmin_tg= ",rmin_sa
     write(11,*)"      rmaj_tg= ",rmaj_sa
     write(11,*)"      alpha_tg= ",alpha_sa
     write(11,*)"      xwell_tg= ",xwell_sa
     write(11,*)"      q_tg= ",q_sa
     write(11,*)"      shat_tg= ",shat_sa
     write(11,*)"      theta0_tg= ",theta0_sa
     write(11,*)"      b_model_tg= ",b_model_sa
     write(11,*)"      ft_model_tg= ",ft_model_sa
  elseif(igeo.eq.1)then
     write(11,*)"      rmin_tg= ",rmin_loc
     write(11,*)"      rmaj_tg= ",rmaj_loc
     write(11,*)"      zmaj_tg= ",zmaj_loc
     write(11,*)"      drmajdx_tg= ",drmajdx_loc
     write(11,*)"      dzmajdx_tg= ",dzmajdx_loc
     write(11,*)"      drmindx_tg= ",drmindx_loc
     write(11,*)"      kappa_tg= ",kappa_loc
     write(11,*)"      s_kappa_tg= ",s_kappa_loc
     write(11,*)"      delta_tg= ",delta_loc
     write(11,*)"      s_delta_tg= ",s_delta_loc
     write(11,*)"      zeta_tg= ",zeta_loc
     write(11,*)"      s_zeta_tg= ",s_zeta_loc
     write(11,*)"      q_tg= ",q_loc
     write(11,*)"      p_prime_tg= ",p_prime_loc
     write(11,*)"      q_prime_tg= ",q_prime_loc
  endif
  write(11,*)"      betae_tg = ",betae_in
  write(11,*)"      damp_psi_tg= ",damp_psi_in
  write(11,*)"      damp_sig_tg= ",damp_sig_in
  write(11,*)"      filter_tg= ",filter_in
  write(11,*)"      debye_tg= ",debye_in
  write(11,*)"      use_bper_tg=",use_bper_in
  write(11,*)"      use_bpar_tg=",use_bpar_in
  write(11,*)"      use_mhd_rule_tg=",use_mhd_rule_in
  write(11,*)"      use_bisection_tg=",use_bisection_in
  write(11,*)"      zeff_tg= ",zeff_in
  write(11,*)"      xnue_tg= ",xnue_in
  write(11,*)"      theta_trapped_tg= ",theta_trapped_in
  write(11,*)"      wdia_trapped_tg = ",wdia_trapped_in
  write(11,*)"      Linsker_factor_tg= ",Linsker_factor_in
  write(11,*)"      gradB_factor_tg= ", gradB_factor_in
  write(11,*)"      wd_zero_tg= ",wd_zero_in
  write(11,*)"      vexb_shear_tg = ",vexb_shear_in
  write(11,*)"      alpha_e_tg= ",alpha_e_in
  write(11,*)"      alpha_p_tg= ",alpha_p_in
  write(11,*)"      alpha_mach_tg= ",alpha_mach_in
  write(11,*)"      alpha_n_tg= ",alpha_n_in
  write(11,*)"      alpha_t_tg= ",alpha_t_in
  write(11,*)"      alpha_quench_tg= ",alpha_quench_in
  write(11,*)"      xnu_factor_tg= ",xnu_factor_in
  write(11,*)"      debye_factor_tg= ",debye_factor_in
  write(11,*)"      etg_factor_tg = ",etg_factor_in
  do is = 1,ns_in
     write(11,*)"      mass_tg(",is,")= ",mass_in(is)
     write(11,*)"      zs_tg(",is,")= ",zs_in(is)
     write(11,*)"      taus_tg(",is,")= ",taus_in(is)
     write(11,*)"      as_tg(",is,")= ",as_in(is)
     write(11,*)"      vpar_tg(",is,")= ",vpar_in(is)
     write(11,*)"      rlns_tg(",is,")= ",rlns_in(is)
     write(11,*)"      rlts_tg(",is,")= ",rlts_in(is)
     write(11,*)"      vpar_shear_tg(",is,")= ",vpar_shear_in(is)
     write(11,*)"      vns_shear_tg(",is,")= ",vns_shear_in(is)
     write(11,*)"      vts_shear_tg(",is,")= ",vts_shear_in(is)
  enddo
  write(11,*)"/"
  do is=1,ns_in
     write(11,*)"! particle_flux(",is,",1) = ",get_particle_flux(is,1)
     write(11,*)"! energy_flux(",is,",1) = ",get_energy_flux(is,1)
     write(11,*)"! stress_tor(",is,",1) = ",get_stress_tor(is,1)
  enddo
  write(11,*)"! trace_path=",(trace_path(is),is=1,7)
  write(11,*)"      nn_max_error_tg=",nn_max_error_in
  close(11)
  !
END SUBROUTINE write_tglf_input
!-----------------------------------------------------------------
!
SUBROUTINE write_wavefunction_out(datafile, tglf_path_in)
  !
  USE tglf_global
  !
  IMPLICIT NONE

  character (len=*), intent(in) :: datafile

  INTEGER :: i,n,k,noff
  REAL :: wave(maxmodes*6)
  CHARACTER(len=81) :: header
  CHARACTER(len=11) :: theta="    theta  "
  CHARACTER(len=22) :: phi="  RE(phi)    IM(phi)  "
  CHARACTER(len=24) :: Bper="  RE(Bper)    IM(Bper)  "
  CHARACTER(len=24) :: Bpar="  RE(Bpar)    IM(Bpar)  "
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf must be called before write_wavefunction_out"
  else
     call get_wavefunction
     !
     open(unit=33,file=trim(tglf_path_in)//datafile,status='replace')
     header = theta//phi
     if(use_bper_in)header = theta//phi//Bper
     if(use_bpar_in)header = theta//phi//Bpar
     if(use_bper_in.and.use_bpar_in)header = theta//phi//Bper//Bpar
     !
     write(33,*)nmodes_out,nfields_out,max_plot
     write(33,*)header
     do i = 1,max_plot
        do n=1,nmodes_out
           noff=2*nfields_out*(n-1)
           wave(noff+1) = REAL(plot_field_out(n,1,i))
           wave(noff+2) = AIMAG(plot_field_out(n,1,i))
           if(use_bper_in)then
              wave(noff+3) = REAL(plot_field_out(n,2,i))
              wave(noff+4) = AIMAG(plot_field_out(n,2,i))
           endif
           if(use_bpar_in)then
              wave(noff+5) = REAL(plot_field_out(n,3,i))
              wave(noff+6) = AIMAG(plot_field_out(n,3,i))
           endif
        enddo
        write(33,*)plot_angle_out(i),(wave(k),k=1,nmodes_out*nfields_out*2)
     enddo
     close(33)
  endif
  !
END SUBROUTINE write_wavefunction_out
!-----------------------------------------------------------------
!
SUBROUTINE write_tglf_sum_flux_spectrum(tglf_path_in)
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
  !   
  IMPLICIT NONE
  CHARACTER(26) :: fluxfile="out.tglf.sum_flux_spectrum"
  INTEGER :: i,j,is,imax,jmax
  REAL :: dky
  REAL :: dky0,dky1,ky0,ky1
  REAL :: pflux0(nsm,3),eflux0(nsm,3)
  REAL :: stress_par0(nsm,3),stress_tor0(nsm,3)
  REAL :: exch0(nsm,3)
  REAL :: pflux1(nsm,3),eflux1(nsm,3)
  REAL :: stress_par1(nsm,3),stress_tor1(nsm,3)
  REAL :: exch1(nsm,3)
  REAL :: pflux_out,eflux_out,mflux_tor_out,mflux_par_out,exch_out
  character (len=256)  :: tglf_path_in

  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_flux_spectrum"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
  !
  !
  ! initialize fluxes
  !
!  write(*,*)"ns0,ns,nky,nmodes",ns0,ns,nky,nmodes_in
!
  do is=ns0,ns
     do j=1,3 
        pflux0(is,j) = 0.0
        eflux0(is,j) = 0.0
        stress_par0(is,j) = 0.0
        stress_tor0(is,j) = 0.0
        exch0(is,j) = 0.0
     enddo
  enddo
  ! set number of field
  jmax = 1
  if(use_Bper_in)jmax=2
  if(use_Bpar_in)jmax=3
  ! loop over species and fields
  do is=ns0,ns
     do j=1,jmax
        write(33,*)"species = ",is,"field =",j
        write(33,*)" particle flux,energy flux,toroidal stress,parallel stress,exchange"

        !
        ! loop over ky spectrum
        !
        iflux_in=.TRUE. 
        dky0=0.0
        ky0=0.0 
        do i=1,nky
           ky_in = ky_spectrum(i)
           dky = dky_spectrum(i)
           ky1=ky_in
           if(i.eq.1)then
              dky1=ky1
           elseif(kygrid_model_in.ne.0)then
              dky = LOG(ky1/ky0)/(ky1-ky0)
              dky1 = ky1*(1.0 - ky0*dky)
              dky0 = ky0*(ky1*dky - 1.0)
           endif
           ! normalize the ky integral to make it independent of the 
           ! choice of temperature and mass scales 
!           dky0 = dky0*SQRT(taus_in(1)*mass_in(2))
!           dky1 = dky1*SQRT(taus_in(1)*mass_in(2))
           !
           ! compute the fluxes in the same way as tglf_TM.f90 
           !
           pflux1(is,j) = 0.0
           eflux1(is,j) = 0.0
           stress_tor1(is,j) = 0.0
           stress_par1(is,j) = 0.0
           exch1(is,j) = 0.0
           do imax = 1,nmodes_in
              pflux1(is,j) = pflux1(is,j) + flux_spectrum_out(1,is,j,i,imax)
              eflux1(is,j) = eflux1(is,j) + flux_spectrum_out(2,is,j,i,imax)
              stress_tor1(is,j) = stress_tor1(is,j) + &
                   flux_spectrum_out(3,is,j,i,imax)
              stress_par1(is,j) = stress_par1(is,j) + &
                   flux_spectrum_out(4,is,j,i,imax)
              exch1(is,j) = exch1(is,j) + flux_spectrum_out(5,is,j,i,imax)
           enddo !imax
           pflux_out = dky0*pflux0(is,j) + dky1*pflux1(is,j)
           eflux_out = dky0*eflux0(is,j) + dky1*eflux1(is,j)
           mflux_tor_out = dky0*stress_tor0(is,j) + dky1*stress_tor1(is,j)
           mflux_par_out = dky0*stress_par0(is,j) + dky1*stress_par1(is,j)
           exch_out = dky0*exch0(is,j) + dky1*exch1(is,j)
           write(33,*)pflux_out,eflux_out,mflux_tor_out,mflux_par_out,exch_out
           pflux0(is,j) = pflux1(is,j)
           eflux0(is,j) = eflux1(is,j)
           stress_par0(is,j) = stress_par1(is,j)
           stress_tor0(is,j) = stress_tor1(is,j)
           exch0(is,j) = exch1(is,j)
           ky0 = ky1
        enddo  ! i
     enddo  ! j
  enddo  ! is 
  !
  CLOSE(33)
  !
END SUBROUTINE write_tglf_sum_flux_spectrum
!-----------------------------------------------------------------

SUBROUTINE write_tglf_density_spectrum(tglf_path_in)
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
  !   
  IMPLICIT NONE
  CHARACTER(25) :: fluxfile="out.tglf.density_spectrum"
  INTEGER :: i,is,n
  REAL :: density(nsm)
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_density_spectrum"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"gyro-bohm normalized density fluctuation amplitude spectra"
  write(33,*)"(density(is),is=1,ns_in)"
  do i=1,nky
    do is=1,ns
      density(is) = 0.0
      do n = 1,nmodes_in
        density(is) = density(is) + intensity_spectrum_out(1,is,i,n)
      enddo
      density(is) = SQRT(density(is))
    enddo
    write(33,*)(density(is),is=1,ns)
  enddo
!
  CLOSE(33)
!
END SUBROUTINE write_tglf_density_spectrum
!-----------------------------------------------------------------

SUBROUTINE write_tglf_temperature_spectrum(tglf_path_in)
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
  !   
  IMPLICIT NONE
  CHARACTER(29) :: fluxfile="out.tglf.temperature_spectrum"
  INTEGER :: i,is,n
  REAL :: temp(nsm)
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_temperature_spectrum"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"gyro-bohm normalized temperature fluctuation amplitude spectra"
  write(33,*)"(temperature(is),is=1,ns_in)"
  do i=1,nky
    do is=1,ns
      temp(is) = 0.0
      do n = 1,nmodes_in
        temp(is) = temp(is) + intensity_spectrum_out(2,is,i,n)
      enddo
      temp(is) = SQRT(temp(is))
    enddo
    write(33,*)(temp(is),is=1,ns)
  enddo
!
  CLOSE(33)
!
END SUBROUTINE write_tglf_temperature_spectrum
!-----------------------------------------------------------------

SUBROUTINE write_tglf_intensity_spectrum(tglf_path_in)
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
  !   
  IMPLICIT NONE
  CHARACTER(36) :: fluxfile="out.tglf.intensity_spectrum"
  INTEGER :: i,is,n
  REAL :: den,tem, u_par, q_tot
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_intensity_spectrum_per_mode"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"gyro-bohm normalized fluctuation intensity spectra per mode"
  write(33,*)"density,temperature,parallel velocity,parallel energy"
  write(33,*)"index limits: ns,nky,nmodes"
  write(33,*)ns,nky,nmodes_in
  do is=ns0,ns
    do i=1,nky
      do n = 1,nmodes_in
        den = intensity_spectrum_out(1,is,i,n)
        tem = intensity_spectrum_out(2,is,i,n)
        u_par = intensity_spectrum_out(3,is,i,n)
        q_tot = intensity_spectrum_out(4,is,i,n)
        write(33,*) den,tem,u_par,q_tot
      enddo
    enddo
  enddo
!
  CLOSE(33)
!
END SUBROUTINE write_tglf_intensity_spectrum

!-----------------------------------------------------------------

SUBROUTINE write_tglf_field_spectrum(tglf_path_in)
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
  !   
  IMPLICIT NONE
  CHARACTER(36) :: fluxfile="out.tglf.field_spectrum"
  INTEGER :: i,n
  REAL :: v,phi,a_par,b_par
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_field_spectrum_per_mode"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"gyro-bohm normalized field fluctuation intensity spectra per mode"
  write(33,*)"vector,potential,A_par,B_par"
  write(33,*)"index limits: nky,nmodes"
  write(33,*)nky,nmodes_in
  if(use_bper_in)then
     write(33,*)"a_par_yes"
  else
     write(33,*)"a_par_no"
  endif
  if(use_bpar_in) then
     write(33,*)"b_par_yes"
  else
     write(33,*)"b_par_no"
  endif
  do i=1,nky
    phi = 0.0
    a_par=0.0
    b_par=0.0
    do n = 1,nmodes_in
      v = field_spectrum_out(1,i,n)
      phi = field_spectrum_out(2,i,n)
      a_par = field_spectrum_out(3,i,n)
      b_par = field_spectrum_out(4,i,n)
      write(33,*)v,phi,a_par,b_par
    enddo
  enddo
!
  CLOSE(33)
!
END SUBROUTINE write_tglf_field_spectrum
!-----------------------------------------------------------------

SUBROUTINE write_tglf_eigenvalue_spectrum(tglf_path_in)
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
  !   
  IMPLICIT NONE
  CHARACTER(28) :: fluxfile="out.tglf.eigenvalue_spectrum"
  INTEGER :: i,n
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_eigenvalue_spectrum"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"gyro-bohm normalized eigenvalue spectra"
  write(33,*)"(gamma(n),freq(n),n=1,nmodes_in)"
  do i=1,nky
    write(33,*)(eigenvalue_spectrum_out(1,i,n), &
              eigenvalue_spectrum_out(2,i,n),n=1,nmodes_in)
  enddo
!
  CLOSE(33)
!
END SUBROUTINE write_tglf_eigenvalue_spectrum
!-----------------------------------------------------------------

SUBROUTINE write_tglf_nete_crossphase_spectrum(tglf_path_in)
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
  !   
  IMPLICIT NONE
  CHARACTER(33) :: fluxfile="out.tglf.nete_crossphase_spectrum"
  INTEGER :: i,n
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_nete_crossphase_spectrum"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"electron density-temperature cross phase spectra per mode"
  write(33,*)"(ne_te_phase,n=1,nmodes_in)"
  do i=1,nky
    write(33,*)(ne_te_phase_spectrum_out(i,n),n=1,nmodes_in)
  enddo
!
  CLOSE(33)
!
END SUBROUTINE write_tglf_nete_crossphase_spectrum
!-----------------------------------------------------------------
!
SUBROUTINE write_tglf_nsts_crossphase_spectrum(tglf_path_in)
  !
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
  !
  IMPLICIT NONE
  CHARACTER(33) :: fluxfile="out.tglf.nsts_crossphase_spectrum"
  INTEGER :: is,i,n
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_nsts_crossphase_spectrum"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"density-temperature cross phase spectra per mode for ",ns_in," species"
   do is=1,ns_in
    write(33,*)"species index = ",is
    write(33,*)"(nsts_phase_spectrum_out,n=1,nmodes_in)"
    do i=1,nky
      write(33,*)(nsts_phase_spectrum_out(is,i,n),n=1,nmodes_in)
    enddo
  enddo
!
  CLOSE(33)
!
END SUBROUTINE write_tglf_nsts_crossphase_spectrum
!-----------------------------------------------------------------
!
 SUBROUTINE write_tglf_QL_flux_spectrum(tglf_path_in)
!
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
! 
  IMPLICIT NONE
  CHARACTER(33) :: fluxfile="out.tglf.QL_flux_spectrum"
  INTEGER :: i,j,k,is,m,jflds
  REAL,PARAMETER :: small=1.0E-10
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_QL_weight_spectrum"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  jflds=1
  if(use_bper_in)jflds = jflds+1
  if(use_bpar_in)jflds = jflds+1
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"QL weights per mode:"
!  write(33,*)"type: 1=particle,2=energy,3=toroidal stress,4=parallel stress,5=exchange"
  write(33,*)"QL_flux_spectrum_out(type,nspecies,field,ky,mode)"
  write(33,*)"index limits: type,ns,field,nky,nmodes"
  write(33,*)5,ns,jflds,nky,nmodes_in
!
!      do j=1,nky
!         do i=1,nmodes_in
!            do is=1,ns
!             do k=1,3
!                  write(33,*)QL_flux_spectrum_out(1,is,k,j,i)
!                  write(33,*)QL_flux_spectrum_out(2,is,k,j,i)
!                  write(33,*)QL_flux_spectrum_out(3,is,k,j,i)
!                  write(33,*)QL_flux_spectrum_out(4,is,k,j,i)
!                  write(33,*)QL_flux_spectrum_out(5,is,k,j,i)
!              enddo
!           enddo
!        enddo
!      enddo
do is=ns0,ns
    do j=1,jflds
       write(33,*)"species = ",is,"field =",j
       do m=1,nmodes_in
         write(33,*)"mode = ",m
         do i=1,nky
           write(33,*)(QL_flux_spectrum_out(k,is,j,i,m),k=1,5)
         enddo  ! i
      enddo ! m
    enddo  ! j
 enddo  ! is

  CLOSE(33)
!
 END SUBROUTINE write_tglf_QL_flux_spectrum
!-----------------------------------------------------------------
!
 SUBROUTINE write_tglf_ky_spectrum(tglf_path_in)
!
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
!
  IMPLICIT NONE
  CHARACTER(20) :: fluxfile="out.tglf.ky_spectrum"
  INTEGER :: i
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_spectral_shift"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"index limits: nky"
  write(33,*)nky
!
  do i=1,nky
    write(33,*)ky_spectrum(i)
  enddo
  CLOSE(33)
!
 END SUBROUTINE write_tglf_ky_spectrum
!-----------------------------------------------------------------
!
 SUBROUTINE write_tglf_spectral_shift_spectrum(tglf_path_in)
!
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
!
  IMPLICIT NONE
  CHARACTER(33) :: fluxfile="out.tglf.spectral_shift_spectrum"
  INTEGER :: i
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_spectral_shift"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"kx spectral shift model is used when ALPHA_QUENCH=0 and ALPHA_E=1.0"
  write(33,*)"note: the model for the spectral shift (kx_e) = <phi|kx/ky|phi>/<phi|phi>"
  write(33,*)"depends on which staturation model is being used: SAT_RULE and UNITS settings"
  write(33,*)"index limits: nky"
  write(33,*)nky
!
  do i=1,nky
    write(33,*)spectral_shift_out(i)
  enddo
  CLOSE(33)
!
 END SUBROUTINE write_tglf_spectral_shift_spectrum
!
!______________________________________________________________
!
SUBROUTINE write_tglf_ave_p0_spectrum(tglf_path_in)
!
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
!
  IMPLICIT NONE
  CHARACTER(33) :: fluxfile="out.tglf.ave_p0_spectrum"
  INTEGER :: i
  character (len=256)  :: tglf_path_in
!
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_ave_p0_spectrum"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
!
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"ave_p0 is used for normalization of SAT0"
  write(33,*)"index limits: nky"
  write(33,*)nky
!
  do i=1,nky
    write(33,*)ave_p0_spectrum_out(i)
  enddo
  CLOSE(33)
!
 END SUBROUTINE write_tglf_ave_p0_spectrum
!
!______________________________________________________________
!
 SUBROUTINE write_tglf_width_spectrum(tglf_path_in)
!
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
!
  IMPLICIT NONE
  CHARACTER(33) :: fluxfile="out.tglf.width_spectrum"
  INTEGER :: i
  character (len=256)  :: tglf_path_in
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_width_spectrum"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=trim(tglf_path_in)//fluxfile,status='replace')
!
  write(33,*)"width of the Gaussian envelope for the Hermite basis functions"
  write(33,*)"index limits: nky"
  write(33,*)nky
!
  do i=1,nky
    write(33,*)width_out(i)
  enddo
  CLOSE(33)
!
 END SUBROUTINE write_tglf_width_spectrum
!-----------------------------------------------------------------
!
 SUBROUTINE write_tglf_scalar_saturation_parameters
!
  USE tglf_dimensions
  USE tglf_global
  USE tglf_species
  USE tglf_kyspectrum
!
  IMPLICIT NONE
  CHARACTER(37) :: fluxfile="out.tglf.scalar_saturation_parameters"
  !
  if(new_start)then
     write(*,*)"error: tglf_TM must be called before write_tglf_saturation_parameters"
     write(*,*)"       NN doesn't compute spectra -> if needed set tglf_nn_max_error_in=-1"
  endif
  !
  OPEN(unit=33,file=fluxfile,status='replace')
!
  write(33,*)"!This file has all of the scalar staturation parameters used for different SAT_RULE options"
  write(33,*)"SAT_RULE = ", sat_rule_in
  write(33,*)"UNITS = ", units_in
  write(33,*)"XNU_MODEL = ",xnu_model_in
  write(33,*)"!   SAT0 model "
  write(33,*)"ETG_FACTOR = ",etg_factor_in
  write(33,*)"B_unit = ",B_unit
  write(33,*)"R_unit = ",R_unit
  write(33,*)"q_unit = ",q_unit
  write(33,*)"!   SAT1 & SAT2 models "
  write(33,*)"ALPHA_ZF = ",alpha_zf_in
  write(33,*)"SAT_geo0_out = ",SAT_geo0_out
  write(33,*)"SAT_geo1_out = ",SAT_geo1_out
  write(33,*)"SAT_geo2_out = ",SAT_geo2_out
  write(33,*)"Bt0_out = ",Bt0_out
  write(33,*)"B_geo0_out = ",B_geo0_out
  write(33,*)"grad_r0_out = ",grad_r0_out
  write(33,*)"rho_ion = ",rho_ion
  write(33,*)"rho_e = ",rho_e
  write(33,*)"kymax_out = ",kymax_out
  write(33,*)"vzf_out = ",vzf_out
!
  CLOSE(33)
!
 END SUBROUTINE write_tglf_scalar_saturation_parameters




