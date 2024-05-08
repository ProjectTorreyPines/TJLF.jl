!------ TGLFEP_read_EXPRO.f90 -----------------
! PURPOSE: Uses EXPRO subroutines to read
!   profile data from input.profiles (not the
!   "input.profile" TGLFEP native format) into
!   TGLFEP profile variables
!----------------------------------------------

subroutine TGLFEP_read_EXPRO

  use mpi
  use TGLFEP_interface
  use TGLFEP_profile

  use expro

  real, parameter :: pi=3.14159265359
  real :: a_meters
  real :: eps_grad
  integer :: i, j, id, ierr, jr_exp
  real, dimension(:), allocatable :: sum0
  real, dimension(:), allocatable :: a_qn_rho

  integer :: values(8)

  real, dimension(:), allocatable :: beta_unit
  real, dimension(:), allocatable :: B_unit

  call MPI_COMM_RANK(TGLFEP_COMM_WORLD,id,ierr)

!  call EXPRO_palloc(TGLFEP_COMM_WORLD,'./',1)

!  EXPRO_ctrl_density_method = 2  ! Enforce quasineutrality
  EXPRO_ctrl_signb = -1.0
  EXPRO_ctrl_signq = 1.0
  EXPRO_ctrl_rotation_method = 1
  EXPRO_ctrl_numeq_flag = 0
!  EXPRO_ctrl_z(:) = 1.0    ! Hardcoded ion charges, D.
  EXPRO_ctrl_quasineutral_flag = 0
!  EXPRO_ctrl_n_ion = 3     ! Thermal+Carbon+beam, not used with EXPRO_ctrl_density_method=2

  rotation_flag = 0 !
!  rotation_flag = 1

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call date_and_time(values=values)
  if (id==0) print '(A21,I2,A1,I2,A1,I2)', 'Before expro_read at ', values(5),':',values(6),':',values(7)
!  print *, 'TGLFEP_COMM_WORLD is', TGLFEP_COMM_WORLD
  call expro_read('input.gacode',TGLFEP_COMM_WORLD)
  call date_and_time(values=values)
  if (id==0) print '(A47,I2,A1,I2,A1,I2)', 'After expro_read, before TGLFEP_scalefactor at ', values(5),':',values(6),':',values(7)
  if (id==0) call expro_write('dump.gacode')
!  print *, 'After calling EXPRO_pread.'

  sign_bt = -1.0 !
  sign_it = -1.0 !
  nr = EXPRO_n_exp !
!  ns = EXPRO_n_ion - 1 + 1   ! Less carbon and add elctron
!  ns = 3     ! Harcoded ions, electrons, EPs
  if (TGLFEP_nion .gt. expro_n_ion) then
    TGLFEP_nion = EXPRO_n_ion
    if (id==0) print *, 'WARNING: TGLFEP_nion reduced to number of ions in input.gacode.'
  endif
  ns = TGLFEP_nion + 1  ! Includes electrons
  if (is_EP .gt. TGLFEP_nion) then
    if (id==0) print *, 'ERROR: EP species index cannot exceed TGLFEP_nion.'
    STOP
  endif
  geometry_flag = 1

  allocate(f_real(nr))
  f_real(:) = 1.0

  EXPRO_rmin(:) = EXPRO_rmin(:) + 1.0E-6   ! Correct for division by zero at origin

  if (l_real_units .eq. 1) f_real(:) = (EXPRO_cs(:)/(EXPRO_rmin(nr)))/(2*3.14159265*1.0E3)

  allocate(beta_unit(nr))

  allocate(zs(ns))
  allocate(mass(ns))

  allocate(as(nr,ns))
  allocate(taus(nr,ns))
  allocate(rlns(nr,ns))
  allocate(rlts(nr,ns))

  allocate(vpar(nr,ns))
  allocate(vpar_shear(nr,ns))
  allocate(a_qn_rho(nr))
  allocate(sum0(nr))

! Quasineutrality logic is added here, enforced among ions i<=TGLFEP_nion  
  sum0(:) = 0.0
  do i = 1, TGLFEP_nion
    if (i .ne. is_EP) sum0(:) = sum0(:) + EXPRO_z(i)*EXPRO_ni(i,:)/EXPRO_ne(:)
  enddo
  a_qn_rho(:) = (1.0 - EXPRO_z(is_EP)*EXPRO_ni(is_EP,:)/EXPRO_ne(:)) / sum0(:)
!  if (id_2==0 .and. id_3==0) print *, 'a_qn_rho', a_qn_rho
  do i = 1, TGLFEP_nion
    if (i .ne. is_EP) EXPRO_ni(i,:) = a_qn_rho(:)*EXPRO_ni(i,:)
  enddo

  do i = 1, scan_n
    if (scan_n .ne. 1) then
      jr_exp = irs + floor(1.0*(i-1)*(nr-irs)/(scan_n-1))
    else
      jr_exp = irs
    endif
  !    if (id==0) print *, 'grid ', i, ' exp_grid ', jr_exp
    ir_exp(i) = jr_exp
  enddo
!  zs(1) = -1.0
  zs(1) = expro_ze
!  mass(1) = 2.7233115E-04    ! me / mD
  mass(1) = expro_masse / 2.0 !
  as(:,1) = 1.0               ! electron density
  taus(:,1) = 1.0             ! and temperature normalized to one
  a_meters = EXPRO_rmin(nr)
  eps_grad = 1.0
  rlns(:,1) = eps_grad*EXPRO_dlnnedr*a_meters
  rlts(:,1) = eps_grad*EXPRO_dlntedr*a_meters

  vpar(:,:) = 0.0
  vpar_shear(:,:) = 0.0

  if (rotation_flag .eq. 1) then
!    print *, 'EXPRO_cs', EXPRO_cs(:)
    do i=1,ns
      vpar(:,i) = EXPRO_w0*EXPRO_rmaj/EXPRO_cs
      vpar_shear(:,i) = -EXPRO_rmaj*EXPRO_w0p*a_meters/EXPRO_cs
!      print *, 'species', i
!      print *, 'vpar', vpar(:,i)
!      print *, 'vpar shear', vpar_shear(:,i)
    enddo
  endif

!  j = 2
!      zs(2) = q_i     ! Read from input.TGLFEP
!      mass(2) = m_i   !
!!      as(:,2) = EXPRO_ni(1,:)/EXPRO_ne(:)
!      as(:,2) = (EXPRO_ne(:)-q_EP*EXPRO_ni(is_EP,:))/EXPRO_ne(:) ! Harcoded quasineutral with D-D
!      taus(:,2) = EXPRO_ti(1,:)/EXPRO_te(:)
!      rlns(:,2) = eps_grad*EXPRO_dlnnidr(1,:)*a_meters    !  These are wrong with quasineutrality on,
!      rlts(:,2) = eps_grad*EXPRO_dlntidr(1,:)*a_meters    !  but they're turned off in TGLFEP_tglf_map, mode_in=2.
!
!      zs(3) = q_EP    ! Read from input.TGLFEP
!      mass(3) = m_EP  !
!      as(:,3) = EXPRO_ni(is_EP,:)/EXPRO_ne(:)
!!      print *, 'as(:,3) in read_EXPRO', as(:,3)
!      taus(:,3) = EXPRO_ti(is_EP,:)/EXPRO_te(:)
  do i = 1, TGLFEP_nion
    zs(i+1) = expro_z(i)    ! Read from input.TGLFEP
    mass(i+1) = expro_mass(i)/2.0  ! Normalized by deuterium mass
    as(:,i+1) = EXPRO_ni(i,:)/EXPRO_ne(:)
!    print *, 'as(:,3) in read_EXPRO', as(:,3)
    taus(:,i+1) = EXPRO_ti(i,:)/EXPRO_te(:)
    rlns(:,i+1) = eps_grad*EXPRO_dlnnidr(i,:)*a_meters    !  These are wrong with quasineutrality on,
    rlts(:,i+1) = eps_grad*EXPRO_dlntidr(i,:)*a_meters    !  but they're turned off in TGLFEP_tglf_map,
                                                          !  mode_in=2.
  enddo

  if (id==0) print *, 'masses:', mass(:)
  if (id==0) print *, 'charges:', zs(:)

  do j = 1, nr                                                ! EMB 10-25-2021
    EXPRO_dlnnidr(is_EP,j) = max(EXPRO_dlnnidr(is_EP,j),1.0)  ! assure (dn_EP/dr)/n_EP >= 1.0
  enddo                                                       !
  rlns(:,is_EP+1) = EXPRO_dlnnidr(is_EP,:)*a_meters
  rlts(:,is_EP+1) = EXPRO_dlntidr(is_EP,:)*a_meters

  allocate(rmin(nr))
  allocate(rmaj(nr))
  allocate(q(nr))
  allocate(shear(nr))

  allocate(q_prime(nr))
  allocate(p_prime(nr))
  allocate(shift(nr))
  allocate(kappa(nr))
  allocate(s_kappa(nr))
  allocate(delta(nr))
  allocate(s_delta(nr))
  allocate(zeta(nr))
  allocate(s_zeta(nr))

  allocate(zeff(nr))
  allocate(betae(nr))

  allocate(rho_star(nr))
  allocate(omega_TAE(nr))
  allocate(B_unit(nr))

  rmin(:) = EXPRO_rmin(:)/a_meters
  rmaj(:) = EXPRO_rmaj(:)/a_meters
  q(:) = q_scale*EXPRO_q(:)
  shear(:) = EXPRO_s(:)
  q_prime(:) = (q(:)/rmin(:))**2 * shear(:)
!  B_unit(:) = abs(EXPRO_b_ref)*(EXPRO_kappa(:)+0.5*EXPRO_kappa(:)*EXPRO_skappa(:))
  B_unit(:) = expro_bunit(:)
  beta_unit(:) = 2.0*4*pi*1.0E-7*EXPRO_ptot(:)/(B_unit(:)**2)
  p_prime(:) = -1.0*(abs(q(:))/rmin(:))*beta_unit(:)/(8*pi)*EXPRO_dlnptotdr(:)*a_meters
  shift(:) = EXPRO_drmaj(:)
  kappa(:) = EXPRO_kappa(:)
  s_kappa(:) = EXPRO_skappa(:)
  delta(:) = EXPRO_delta(:)
  s_delta(:) = EXPRO_sdelta(:)
  zeta(:) = EXPRO_zeta(:)
  s_zeta(:) = EXPRO_szeta(:)
  zeff(:) = EXPRO_z_eff(:)
  betae(:) = beta_unit(:)*(1.6022*1.0E3*EXPRO_ne(:)*EXPRO_te(:)/EXPRO_ptot(:))
  rho_star(:) = EXPRO_rhos(:)/a_meters
  omega_TAE(:) = (2./betae(:))**0.5/2./q(:)/rmaj(:)
end
