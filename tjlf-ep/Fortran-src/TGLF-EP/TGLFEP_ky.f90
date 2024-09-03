!------------------------------------------------------------ TGLFEP_ky.f90
!
! PURPOSE: Calculate the growth rate and frequency at a single ky
!------------------------------------------------------------

subroutine TGLFEP_ky

  use tglf_interface 
  use tglf_pkg 
  use TGLFEP_interface
  use tglf_max_dimensions

  implicit none 
  logical :: iexist

  integer :: i,n, jfields, kcomp, j_ion

  real :: g(nmodes),f(nmodes)

  integer :: nmodes_out 
  integer :: nfields_out 
  integer :: max_plot_out 
  real, dimension(max_plot) :: angle_out
  complex, dimension(maxmodes,3,max_plot) :: wavefunction_out  
  real, dimension(20) :: wave_write 
  real, dimension(4) :: x_tear_test 
  real :: wave_max 
  integer, dimension(1) :: wave_max_loc
  integer :: n_balloon_pi
  integer :: i_mid_plot
  integer :: n_out
  real :: max_phi 
  real :: max_apar 
  real :: max_field 
  real :: phase 
  real, dimension(4) :: DEP
  real, dimension(4) :: chi_th 
  real, dimension(4) :: chi_i 
  real, dimension(4) :: chi_i_cond 
  real, dimension(4) :: chi_e 
  real, dimension(4) :: chi_e_cond

  real :: EP_QL_flux
  real :: EP_QL_e_flux
  real :: th_QL_flux
  real :: i_QL_flux
  real, dimension(4) :: i_QL_cond_flux
  real :: e_QL_flux
  real, dimension(4) :: e_QL_cond_flux 
  real :: th_eff_grad
  real :: i_eff_grad
  real, dimension(4) :: QL_flux_ratio
  real, dimension(4) :: EP_conv_frac
  real, dimension(4) :: theta_2_moment
  real :: ef_phi_norm

  call TGLFEP_tglf_map

!  print *, 'In TGLFEP_ky (1) factor_in=', factor_in

  tglf_use_transport_model_in = .false.

  if(process_in .eq. 0) then
    inquire(file='input.ky',exist=iexist)
    if(iexist) then
      open(unit=44,file='input.ky',status='old')
      read(44,*) ky_in
      close(44)
    else
      print *, 'input.ky file not found'
      stop
    endif
    
!    tglf_write_wavefunction_flag_in = 1
  endif

!  tglf_write_wavefunction_flag_in = 0  ! Test to write wavefunction

  tglf_kygrid_model_in = 0
  tglf_ky_in           = ky_in

  tglf_nmodes_in = nmodes

!  tglf_nbasis_min_in = 16
!  tglf_nbasis_max_in = 16
!  tglf_nbasis_min_in = 32
!  tglf_nbasis_max_in = 32
  tglf_nbasis_min_in = n_basis
  tglf_nbasis_max_in = n_basis
!  tglf_nbasis_min_in = 4
!  tglf_nbasis_max_in = 4
!  tglf_nbasis_min_in = 2
!  tglf_nbasis_max_in = 2
!  tglf_nbasis_min_in = 8
!  tglf_nbasis_max_in = 8


  tglf_nxgrid_in     = 32
!  tglf_nxgrid_in     = 16
  
  tglf_width_in      = width_in
  tglf_find_width_in = .false.

!  print *, 'In TGLFEP_ky (2) factor_in=', factor_in
  call tglf_run
!  print *, 'In TGLFEP_ky (3) factor_in=', factor_in

  do n=1,nmodes
    g(n) = get_growthrate(n)
    f(n) = get_frequency(n)
  enddo

  do n = 1, nmodes
    lkeep(n) = f(n) .lt. freq_AE_upper
    lkeep(n) = lkeep(n) .and. (g(n) .gt. gamma_thresh)
  enddo

  call get_wavefunction_out(nmodes_out,nfields_out,max_plot_out,angle_out,wavefunction_out)
!  print *, 'In TGLFEP_ky (4) factor_in=', factor_in
  ltearing(:) = .false.
  l_i_pinch(:) = .false.
  l_e_pinch(:) = .false.
  l_EP_pinch(:) = .false.
  l_th_pinch(:) = .false.
  l_QL_ratio(:) = .false.
!  l_max_outer_panel(:) = .false.
  do n = 1, nmodes
    x_tear_test(n) = 0.0
    wave_max=maxval(abs(wavefunction_out(n,1,:)))+1.0E-3
    wave_max_loc=maxloc(abs(wavefunction_out(n,1,:)))
    n_balloon_pi = (max_plot-1)/9
    i_mid_plot = (max_plot-1)/2 + 1
!    l_max_outer_panel(n) = (wave_max_loc(1) .lt. (i_mid_plot-n_balloon_pi)) .or. (wave_max_loc(1) .gt. (i_mid_plot+n_balloon_pi))
    theta_2_moment(n) = 0.0
    ef_phi_norm = 0.0
    do i = 1, max_plot
       x_tear_test(n) = max(x_tear_test(n),abs(wavefunction_out(n,1,i)-wavefunction_out(n,1,max_plot+1-i))/wave_max)
       ef_phi_norm = ef_phi_norm + abs(wavefunction_out(n,1,i))
       theta_2_moment(n) = theta_2_moment(n) + (9*3.14159265*(-1.0+(2.0*(i-1))/(max_plot-1)))**2*abs(wavefunction_out(n,1,i))
    enddo
    theta_2_moment(n) = theta_2_moment(n) / ef_phi_norm
    if (x_tear_test(n) .GT. 1.0E-1) ltearing(n) = .true.
    EP_QL_flux = 0.0    ! Particle flux for EPs
    EP_QL_e_flux = 0.0  ! Energy flux for EPs
    i_QL_flux = 0.0     ! Energy flux for e- and ions
    i_QL_cond_flux(n) = 0.0
    i_eff_grad = 0.0
    e_QL_flux = 0.0
    e_QL_cond_flux(n) = 0.0
    th_QL_flux = 0.0
    th_eff_grad = 0.0
    do jfields = 1, 3
      EP_QL_flux = EP_QL_flux + get_QL_particle_flux(n,is_EP+1,jfields)
      EP_QL_e_flux = EP_QL_e_flux + get_QL_energy_flux(n,is_EP+1,jfields)
      e_QL_flux  = e_QL_flux + get_QL_energy_flux(n,1,jfields)
      e_QL_cond_flux(n) = e_QL_cond_flux(n) +  get_QL_energy_flux(n,1,jfields) -  &
                      1.5*tglf_taus_in(1)*get_QL_particle_flux(n,1,jfields)
      do j_ion = 2, is_EP
        i_QL_flux = i_QL_flux + get_QL_energy_flux(n,j_ion,jfields)
        i_QL_cond_flux(n) = i_QL_cond_flux(n) + get_QL_energy_flux(n,j_ion,jfields) - &
                          1.5*tglf_taus_in(j_ion)*get_QL_particle_flux(n,j_ion,jfields)
      enddo
    enddo
    do j_ion = 2, is_EP
      i_eff_grad = i_eff_grad + tglf_rlts_in(j_ion)*tglf_as_in(j_ion)*tglf_taus_in(j_ion)
    enddo
    th_eff_grad = i_eff_grad + tglf_rlts_in(1)*tglf_as_in(1)
!    th_QL_flux = i_QL_flux + e_QL_flux
    th_QL_flux = i_QL_cond_flux(n) + e_QL_cond_flux(n)                     ! Thermal species inverse gradient
    DEP(n) = EP_QL_flux / (tglf_as_in(is_EP+1)*tglf_rlns_in(is_EP+1))      ! lengths are set to 1.0E-6
    chi_e(n) = e_QL_flux / (tglf_rlts_in(1)*tglf_as_in(1))                 ! in usual operation.
    chi_e_cond(n) = e_QL_cond_flux(n) / (tglf_rlts_in(1)*tglf_as_in(1))
    chi_i(n) = i_QL_flux / i_eff_grad
    chi_i_cond(n) = i_QL_cond_flux(n) / i_eff_grad
    chi_th(n) = th_QL_flux / th_eff_grad
!    QL_flux_ratio(n) = (EP_QL_flux/tglf_as_in(is_EP+1))/(abs(i_QL_cond_flux(n))/(tglf_as_in(1)-tglf_as_in(is_EP+1)))
!    QL_flux_ratio(n) = DEP(n)/chi_i_cond(n)
    QL_flux_ratio(n) = EP_QL_e_flux / abs(i_QL_cond_flux(n))
    EP_conv_frac(n) = EP_QL_flux * 1.5 * tglf_taus_in(is_EP+1) / EP_QL_e_flux
    if (chi_i(n) .lt. 0.0) l_i_pinch(n) = .true.
    if (chi_e(n) .lt. 0.0) l_e_pinch(n) = .true.
    if (chi_th(n) .lt. 0.0) l_th_pinch(n) = .true.
    if (DEP(n) .lt. 0.0) l_EP_pinch(n) = .true.
    if (QL_flux_ratio(n) .lt. QL_ratio_thresh) l_QL_ratio(n) = .true. 
    if (theta_2_moment(n) .gt. theta_sq_thresh) l_theta_sq(n) = .true.
  enddo

  do n = 1, nmodes
    if (reject_tearing_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. ltearing(n)
    if (reject_i_pinch_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_i_pinch(n)
    if (reject_e_pinch_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_e_pinch(n)
    if (reject_th_pinch_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_th_pinch(n)
    if (reject_EP_pinch_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_EP_pinch(n)
!    if (reject_max_outer_panel_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_max_outer_panel(n)
    lkeep(n) = lkeep(n) .and. .not. l_QL_ratio(n)
    lkeep(n) = lkeep(n) .and. .not. l_theta_sq(n)
  enddo

  if (l_wavefunction_out .eq. 1) then
    open (unit=22,file=trim(str_wf_file),status='replace')
    write (22,*) 'nmodes=', nmodes_out, ' nfields=', nfields_out, ' max_plot=', max_plot_out
    write (22,*) 'ky=', ky_in, '  width=', width_in
    write (22,*) 'theta     ', '((Re(field_i), Im(field_i),i=(1,nfields)),j=1,nmodes)'
    write (22,'(A,4E16.6E3,A)') 'Tearing metric:  [', x_tear_test(:), ' ]'
    write (22,'(A,4E16.6E3,A)') 'DEP: [' , DEP(:) , ']'
    write (22,'(A,4E16.6E3,A)') 'chi_th: [', chi_th(:), ']'
    write (22,'(A,4E16.6E3,A)') 'chi_i: [', chi_i(:), ']'
    write (22,'(A,4E16.6E3,A)') 'chi_i_cond: [', chi_i_cond(:), ']'
    write (22,'(A,4E16.6E3,A)') 'chi_e: [', chi_e(:), ']'
    write (22,'(A,4E16.6E3,A)') 'chi_e_cond: [', chi_e_cond(:), ']'
    write (22,'(A,4E16.6E3,A)') 'i_QL_cond_flux: [', i_QL_cond_flux(:), ']'
    write (22,'(A,4E16.6E3,A)') 'e_QL_cond_flux: [', e_QL_cond_flux(:), ']'
    write (22,'(A,4E16.6E3,A)') 'QL_ratio: [', QL_flux_ratio(:), ']'
    write (22,'(A,4E16.6E3,A)') 'EP QL convective fraction: [', EP_conv_frac(:), ']'
    write (22,'(A,4E16.6E3,A)') '<theta^2>: [', theta_2_moment(:), ']'
    write (22,*) lkeep(:)
! First, renormalize and adjust phases
    n_out = 0
    do n = 1, nmodes_out
      max_phi = maxval(abs(wavefunction_out(n,1,:)))
      max_apar = maxval(abs(wavefunction_out(n,2,:)))
      max_field = max(max_phi,max_apar)
      phase = atan2(aimag(wavefunction_out(n,1,(max_plot_out+1)/2)),real(wavefunction_out(n,1,(max_plot_out+1)/2)))
      do jfields = 1, nfields_out
        wavefunction_out(n,jfields,:) = wavefunction_out(n,jfields,:)/(max_field*exp(cmplx(0,1)*phase))
      enddo
      if ((n_out.eq.0) .and. lkeep(n)) n_out = n
    enddo

    if (n_out .eq. 0) then
      n_out = 1
      write (22,*) "No kept modes at nominal write parameters. Showing leading mode."   
    endif

! Write renormalized, re-phased eigenfunctions out to str_wf_file.
    do i = 1, max_plot
!      kcomp = 0
!      wave_write(:) = 0.0
!      do n = 1, nmodes_out
!        do jfields = 1, nfields_out
!          kcomp = kcomp + 1
!          wave_write(kcomp) = real(wavefunction_out(n,jfields,i))
!          kcomp = kcomp + 1
!          wave_write(kcomp) = aimag(wavefunction_out(n,jfields,i))
!        enddo
!      enddo
      write (22,'(F10.6,4E16.6E3)') angle_out(i), real(wavefunction_out(n_out,1,i)), aimag(wavefunction_out(n_out,1,i)),&
                                   real(wavefunction_out(n_out,2,i)), aimag(wavefunction_out(n_out,2,i)) 
                ! Just first kept mode with two fields hardcoded for now
    enddo
    close (unit=22)
  endif 
end subroutine TGLFEP_ky
