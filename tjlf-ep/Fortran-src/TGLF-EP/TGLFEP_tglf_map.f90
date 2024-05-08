module TGLFEP_profile



  implicit none



  real :: sign_bt,sign_it

  integer :: nr,ns,geometry_flag,  &

             rotation_flag

  real,allocatable,dimension(:) ::  zs,mass

  real,allocatable,dimension(:,:) :: as,taus,rlns,rlts, &

                                     vpar,              &

                                     vpar_shear

  real,allocatable,dimension(:) :: rmin,rmaj,shift,&

                                   q,shear,alpha,  &

                                   q_prime,p_prime,&

                                   kappa,s_kappa,  &

                                   delta,s_delta,  &

                                   zeta,s_zeta,    &

                                   zeff,betae,     &

                                   rho_star,omega_TAE



  integer :: is !energetic particle species (density scan)



  integer :: scan_n = 1  ! How many points are actually run.

  integer :: irs         ! First run point in profile.

  real :: a_qn  ! quasineutrality scale factor for non-EP species



end module TGLFEP_profile



subroutine TGLFEP_tglf_map 



  use tglf_interface

  use TGLFEP_interface

  use TGLFEP_profile



  implicit none



  integer :: i

  real :: sum0,sum1

 

  if(ir .lt. 1 .or. ir .gt. nr) then

    print *, 'ir=', ir, ' is beyond the boundary with nr=', nr

    stop

  endif



  !if(mode_in .lt. 1 .or. mode_in .gt. 4) then

  !  write(*,*) 'mode_flag must be 1. EP+ITG/TEM drive'

  !  write(*,*) '                  2. EP drive only'

  !  write(*,*) '                  3. ITG/TEM drive only'

  !  write(*,*) '                  4. ITG/TEM in EP+ITG/TEM drive'

  !  write(*,*) 'mode_flag has been set as default 1.'

  !endif



  tglf_sign_bt_in = sign_bt

  tglf_sign_it_in = sign_it



  tglf_sat_rule_in = 0

!  tglf_sat_rule_in = 1



  tglf_ns_in = ns



!  is = ns  !the density scan species

  is = is_EP + 1



  tglf_geometry_flag_in = geometry_flag



  do i = 1,ns

    tglf_zs_in(i)   = zs(i)

    tglf_mass_in(i) = mass(i)

    tglf_as_in(i)   = as(ir,i)

!    print *, 'as(',ir,',',i,')=',as(ir,i)

!    print *, '(1)', 'tglf_as_in(',i,')=', tglf_as_in(i), 'at ir=',ir

    tglf_taus_in(i) = taus(ir,i)

  enddo

!  print *, 'In tglf_map (1) factor_in=', factor_in

  tglf_zs_in(1) = -1.0 !electron



  factor_max = 0.5*1.0/(zs(is)*as(ir,is))

  if (scan_method .eq. 2) factor_max = 1.0E3

  if(factor_in .lt. 0.) factor_in = 0.0

  if(factor_in .gt. factor_max) factor_in = factor_max

  !factor_in = nint(100000*factor_in)/100000

  if (scan_method .eq. 1) then
    tglf_as_in(is) = as(ir,is)*factor_in   !  if logic EMB 2/2023
  endif


!  print *, '(2)', 'tglf_as_in(',is,')=', tglf_as_in(is), 'at ir=',ir, 'factor_in=',factor_in,'factor_max=',factor_max



                 ! EMB changed quasineutrality enforcement to

  sum0 = 0.      ! include all species but EP 7/25/2022

  do i = 2,ns

    if (i .ne. is) sum0 = sum0 + tglf_zs_in(i)*tglf_as_in(i)

  enddo
  if (ir .eq. 3) then
    print *, "=============="
    print *, a_qn
    print *, tglf_zs_in
  endif
  a_qn = (1.0 - tglf_zs_in(is)*tglf_as_in(is)) / sum0

!  tglf_as_in(2) = (1.0 - sum0)/tglf_zs_in(2) !electric neutrality
  if (ir .eq. 3) then
    print *, sum0
    print *, a_qn
    print *, tglf_as_in
    print *, "=============="
  endif

  do i = 2,ns
    if (i .ne. is) tglf_as_in(i) = a_qn*tglf_as_in(i)
  enddo

  if (ir .eq. 3) then
    print *, tglf_as_in
  endif
  !print *, "is, factor_max, factor_in, sum0, AS:"
  !print *, is, "/", factor_max, "/", factor_in, "/", sum0, "/", tglf_as_in

  if(mode_in .eq. 2) then !EP drive only

    do i = 1,ns  ! EMB changed 2 to ns and added "if" logic 6/17/2022

      if (i .ne. is) then

!        tglf_rlns_in(i) = 0.

!        tglf_rlts_in(i) = 0.

        tglf_rlns_in(i) = 1.0E-6

        tglf_rlts_in(i) = 1.0E-6

      else

        tglf_rlns_in(i) = rlns(ir,i)

        tglf_rlts_in(i) = rlts(ir,i)

      endif

!      if (id_2==0) print *, 'species:', i, '   grid:', ir, '  rlns:', tglf_rlns_in(i), '  rlts:', tglf_rlts_in(i)

    enddo

  else

    do i = 1,ns      ! EMB changed 2 to ns 6/17/2022

      tglf_rlns_in(i) = rlns(ir,i)

      tglf_rlts_in(i) = rlts(ir,i)

    enddo

  endif



  !EP species

  ! For multiple EP species, mode_in 3 will turn off

  ! drives for species index >= is_EP (is in TGLF-EP)

  if(mode_in .eq. 3) then !ITG/TEM drive only  

    do i = is, ns   ! EMB changed 3 to is 6/17/2022

!      tglf_rlns_in(i) = 1.0E-6

!      tglf_rlts_in(i) = 1.0E-6

      tglf_rlns_in(i) = 1.0E-5

      tglf_rlts_in(i) = 1.0E-5

    enddo



    !tglf_filter_in = 2.0

  else

!    do i = is, ns                       !! EMB commented out 6/18/2022 and 8/24/2022

!      tglf_rlns_in(i) = rlns(ir,i)       ! Changed 3 to is 8/18/2022

!      tglf_rlts_in(i) = rlts(ir,i)

!      if (id_2==0) print *, 'species:', i, '   grid:', ir, '  rlns:', tglf_rlns_in(i), '  rlts:', tglf_rlts_in(i) 

!    enddo



    tglf_rlns_in(is) = rlns(ir,is)*scan_factor

!    if (id_2==0) print *, 'grid:', ir, '    scan factor:', scan_factor



    tglf_filter_in = 0.0   !The frequency threshold off

  endif



  if (scan_method .eq. 2) tglf_rlns_in(is) = factor_in*rlns(ir,is)    ! EMB 2/2023



  !----------------------------------------------------------------

  ! Geometry parameters:

  if(geometry_flag .eq. 0) then

    ! s-alpha

    tglf_rmin_sa_in     = rmin(ir)

    tglf_rmaj_sa_in     = rmaj(ir)

    tglf_q_sa_in        = q(ir)

    tglf_shat_sa_in     = shear(ir)

    tglf_alpha_sa_in    = alpha(ir)

    tglf_xwell_sa_in    = 0.0

    tglf_theta0_sa_in   = 0.0

    tglf_b_model_sa_in  = 1

    tglf_ft_model_sa_in = 1

  endif

  if(geometry_flag .eq. 1) then

    ! Miller

    tglf_rmin_loc_in    = rmin(ir)

    tglf_rmaj_loc_in    = rmaj(ir)

    tglf_zmaj_loc_in    = 0.0

    tglf_drmajdx_loc_in = shift(ir)

    tglf_dzmajdx_loc_in = 0.0

    tglf_kappa_loc_in   = kappa(ir)

    tglf_s_kappa_loc_in = s_kappa(ir)

    tglf_delta_loc_in   = delta(ir)

    tglf_s_delta_loc_in = s_delta(ir)

    tglf_zeta_loc_in    = zeta(ir)

    tglf_s_zeta_loc_in  = s_zeta(ir)

    tglf_q_loc_in       = abs(q(ir))

    tglf_q_prime_loc_in = q_prime(ir)



    sum0 = 0.

    do i = 1,ns

      sum0 = sum0 + tglf_as_in(i)*tglf_taus_in(i)*(tglf_rlns_in(i)+tglf_rlts_in(i))

    enddo

    sum1 = 0.

    do i = 1,ns

      sum1 = sum1 + as(ir,i)*taus(ir,i)*(rlns(ir,i)+rlts(ir,i))

    enddo

    tglf_p_prime_loc_in = p_prime(ir)*sum0/sum1

    ! It is usually negative comfirmed by Gary 7.13.2016

  endif

  !----------------------------------------------------------------



  if (rotation_flag .eq. 1) then

    tglf_vpar_in(1:ns) = vpar(ir,1:ns)

!    print *, 'radius', ir, 'tglf_vpar_in', tglf_vpar_in(:)

    tglf_vpar_shear_in(1:ns) = vpar_shear(ir,1:ns)

  else

    tglf_vpar_in(1:ns) = 0.0

    tglf_vpar_shear_in(1:ns) = 0.0

  endif



  tglf_use_bper_in = .true.

!  tglf_use_bper_in = .false.

  tglf_betae_in = betae(ir)



  tglf_xnue_in = 0.

  tglf_zeff_in = zeff(ir)



  if(mode_in .eq. 4) then !ITG/TEM in EP+ITG/TEM drive

    tglf_filter_in = 2.0

  endif



  select case(ky_model)

  case (0)

    ky_in = 0.01*n_toroidal

  case (1)

    ky_in = n_toroidal*tglf_q_loc_in/tglf_rmin_loc_in*rho_star(ir) !ky = n*q/(r/a)*rho_star

  case (2)

    ky_in = n_toroidal*0.1*tglf_zs_in(is)/sqrt(tglf_mass_in(is)*tglf_taus_in(is)) !ky_ep = 0.1*n

  case (3)

    ky_in = kyhat_in*tglf_zs_in(is)/sqrt(tglf_mass_in(is)*tglf_taus_in(is))

  end select



  freq_AE_upper = freq_cutoff*abs(omega_TAE(ir)) !freq_cutoff*omega_TAE

  
  !!print *, tglf_ibranch_in
  !print *, tglf_width_in

!  if (l_print) print *, 'omega_TAE ', omega_TAE, ' freq_cutoff ', freq_AE_upper



!  print *, 'In tglf_map (3) factor_in=', factor_in

  !print *, a_qn, " ", tglf_as_in, " ", ky_in

end subroutine TGLFEP_tglf_map



subroutine read_input_profile



  use TGLFEP_profile



  implicit none

  logical :: iexist

  integer :: i,j



  inquire(file='input.profile',exist=iexist)

  if(iexist) then

    open(unit=55,file='input.profile',status='old')

    read(55,*) sign_bt

    read(55,*) sign_it

    read(55,*) nr

    read(55,*) ns

    read(55,*) geometry_flag



    allocate(zs(ns))

    allocate(mass(ns))



    allocate(as(nr,ns))

    allocate(taus(nr,ns))

    allocate(rlns(nr,ns))

    allocate(rlts(nr,ns))



    do i = 1,ns

      read(55,*)

      read(55,*)

      read(55,*) zs(i)

      read(55,*) mass(i)



      if(i .ne. 1) then

        read(55,*)

        do j = 1,nr

          read(55,*) as(j,i)

        enddo

        read(55,*)

        do j = 1,nr

          read(55,*) taus(j,i)

        enddo

      else

        as(:,i)   = 1.0  !electron density

        taus(:,i) = 1.0  !electron temperature

      endif



      read(55,*)

      do j = 1,nr

        read(55,*) rlns(j,i)

      enddo

      read(55,*)

      do j = 1,nr

        read(55,*) rlts(j,i)

      enddo

    enddo



    allocate(rmin(nr))

    allocate(rmaj(nr))

    allocate(q(nr))

    allocate(shear(nr))



    if(geometry_flag .eq. 0) then

      allocate(alpha(nr))

    endif



    if(geometry_flag .eq. 1) then

      allocate(q_prime(nr))

      allocate(p_prime(nr))

      allocate(shift(nr))

      allocate(kappa(nr))

      allocate(s_kappa(nr))

      allocate(delta(nr))

      allocate(s_delta(nr))

      allocate(zeta(nr))

      allocate(s_zeta(nr))

    endif



    allocate(zeff(nr))

    allocate(betae(nr))



    allocate(rho_star(nr))

    allocate(omega_TAE(nr))



    read(55,*)

    read(55,*)

    if(geometry_flag .eq. 0) then



      read(55,*)

      do i = 1,nr

        read(55,*) rmin(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) rmaj(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) q(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) shear(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) alpha(i)

      enddo

    

    endif

    

    if(geometry_flag .eq. 1) then



      read(55,*)

      do i = 1,nr

        read(55,*) rmin(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) rmaj(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) q(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) shear(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) q_prime(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) p_prime(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) shift(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) kappa(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) s_kappa(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) delta(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) s_delta(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) zeta(i)

      enddo

      read(55,*)

      do i = 1,nr

        read(55,*) s_zeta(i)

      enddo



    endif



    read(55,*)

    read(55,*)

    do i = 1,nr

      read(55,*) zeff(i)

    enddo

    read(55,*)

    do i = 1,nr

      read(55,*) betae(i)

    enddo

    read(55,*)

    !As Gary Staebler suggested, change rho_star to ky1 (n=1)

    do i = 1,nr

      read(55,*) rho_star(i)

    enddo

    !ky = (nq/r)*rhos

    do i = 1,nr

      rho_star(i) = rho_star(i)/q(i)*rmin(i)

    enddo

    !read(21,*)

    ! do i = 1,nr

    !   read(21,*) omega_TAE(i)

    ! enddo

    !omega_TAE don't need to read from file, 4.24.2017

    !since omega_TAE/(c_s/a) = sqrt(2/beta_e)/2/q/(R/a)

    do i = 1,nr

      omega_TAE(i) = (2./betae(i))**0.5/2./q(i)/rmaj(i)

    enddo



    close(55)



  else



    print *, 'input.profile file not found'

    stop



  endif



end subroutine read_input_profile



subroutine dump_profile



  use TGLFEP_profile



  implicit none

  integer :: i,j



  open(unit=21,file='dump.profile',status='replace')

  write(21,*) sign_bt,'  SIGN_BT'

  write(21,*) sign_it,'  SIGN_IT'

  write(21,*) nr,'  NR'

  write(21,*) ns,'  NS'

  write(21,*) geometry_flag, '  GEOMETRY_FLAG'



  do i = 1,ns

    if(i .eq. 1) then

      write(21,*) '--------------------------------------------------------------'

      write(21,*) '# electron species '

    else

      write(21,*) '--------------------------------------------------------------'

      write(21,*) '# ion species ',i-1

    endif



    write(21,*) zs(i),'  ZS'

    write(21,*) mass(i),'  MASS m/m_D'



    if(i .ne. 1) then

      write(21,*) '# normalized density: as'

      write(21,10) as(:,i)

      write(21,*) '# normalized temperature: taus'

      write(21,10) taus(:,i)

    endif



    write(21,*) '# normalized density gradients: rlns'

    write(21,10) rlns(:,i)

    write(21,*) '# normalized temperature gradients: rlts'

    write(21,10) rlts(:,i)

  enddo



  write(21,*) '--------------------------------------------------------------'

  write(21,*) '# Geometry '

  if(geometry_flag .eq. 0) then

    write(21,*) '# minor radius: rmin'

    write(21,10) rmin

    write(21,*) '# major radius: rmaj'

    write(21,10) rmaj

    write(21,*) '# safety factor: q'

    write(21,10) q

    write(21,*) '# magnetic shear: shear'

    write(21,10) shear

    write(21,*) '# normalized pressure gradient: alpha'

    write(21,10) alpha

  else if(geometry_flag .eq. 1) then

    write(21,*) '# minor radius: rmin'

    write(21,10) rmin

    write(21,*) '# major radius: rmaj'

    write(21,10) rmaj

    write(21,*) '# safety factor: q'

    write(21,10) q

    write(21,*) '# magnetic shear: shear'

    write(21,10) shear

    write(21,*) '# q_prime'

    write(21,10) q_prime

    write(21,*) '# p_prime'

    write(21,10) p_prime

    write(21,*) '# shift'

    write(21,10) shift

    write(21,*) '# elogation: kappa'

    write(21,10) kappa

    write(21,*) '# shear in elogation: s_kappa'

    write(21,10) s_kappa

    write(21,*) '# triangularity: delta'

    write(21,10) delta

    write(21,*) '# shear in triangularity: s_delta'

    write(21,10) s_delta

    write(21,*) '# squareness: zeta'

    write(21,10) zeta

    write(21,*) '# shear in squareness: s_zeta'

    write(21,10) s_zeta

  endif



  write(21,*) '--------------------------------------------------------------'

  write(21,*) '# effective ion charge: zeff'

  write(21,10) zeff

  write(21,*) '# betae'

  write(21,10) betae

  !write(21,*) '# rho_star = rho_s/a'

  write(21,*) '# ky = (nq/r)*rho_s for n = 1'

  write(21,10) rho_star*q/rmin

  !write(21,*) '# omega_TAE / (c_s/a)'

  !write(21,10) omega_TAE



  close(21)



10  format(F11.6)



end subroutine dump_profile



subroutine deallocate_profile



  use TGLFEP_profile



  implicit none



  deallocate(zs)

  deallocate(mass)



  deallocate(as)

  deallocate(taus)

  deallocate(rlns)

  deallocate(rlts)



  deallocate(rmin)

  deallocate(rmaj)

  deallocate(q)

  deallocate(shear)



  if(geometry_flag .eq. 0) then

    deallocate(alpha)

  endif



  if(geometry_flag .eq. 1) then

    deallocate(q_prime)

    deallocate(p_prime)

    deallocate(shift)



    deallocate(kappa)

    deallocate(s_kappa)

    deallocate(delta)

    deallocate(s_delta)

    deallocate(zeta)

    deallocate(s_zeta)

  endif



  deallocate(zeff)

  deallocate(betae)



  deallocate(rho_star)

  deallocate(omega_TAE)



end subroutine deallocate_profile

