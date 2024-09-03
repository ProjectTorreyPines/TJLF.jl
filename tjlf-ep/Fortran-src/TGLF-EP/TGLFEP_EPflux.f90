!------------------------------------------------------------
! TGLFEP_EPflux.f90
!
! PURPOSE:
!  Calculate the quasilinear flux of trace EPs
!------------------------------------------------------------

subroutine TGLFEP_EPflux

  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface
  use TGLFEP_profile

!  use EXPRO_interface
  use expro

  implicit none
  logical :: iexist

  integer :: i,n

  real :: g(nmodes)
  real :: f(nmodes)

  real :: chi_eff

!  allocate(DEP_trace(nr))
!  allocate(QL_ratio(nr))
!  allocate(DEP_trace_complete(nr))
!  allocate(QL_ratio_complete(nr))
!  allocate(chi_gB(nr))

  DEP_trace(:) = 0.0
  QL_ratio(:)  = 0.0
  chi_gB(:)    = 0.0

  call TGLFEP_tglf_map

  tglf_use_transport_model_in = .true.

  tglf_kygrid_model_in = 1    ! standard TGLF ky spectrum
!  tglf_kygrid_model_in = 0     ! Dense spectrum
!  tglf_nky_in = 80             ! At low-k
!  tglf_ky_in = 4.0

!  tglf_nmodes_in = 2    ! TGLF default
  tglf_nmodes_in = 4    

  tglf_nbasis_min_in = 2
  tglf_nbasis_max_in = 4

  tglf_nxgrid_in     = 16
  
  tglf_width_in = 0.3   ! TGLF default
  tglf_find_width_in = .true.  ! use TGLF width-searching algorithm

  tglf_as_in(3) = 1.0E-6*tglf_as_in(3)  ! Set EP density to trace
  tglf_rlts_in(1) = (jtscale*tscale_interval)*tglf_rlts_in(1)  ! Scan thermal T gradient scale in increments of tscale_interval
  tglf_rlts_in(2) = (jtscale*tscale_interval)*tglf_rlts_in(2)  ! up to jtscale_max*tscale_interval 

  call tglf_run

  do n=1,nmodes
    g(n) = get_growthrate(n)
    f(n) = get_frequency(n)
  enddo

  chi_gB = EXPRO_rhos(:)**2 * (EXPRO_cs(:)/EXPRO_rmin(nr))

  DEP_trace_loc = chi_gb(ir)*tglf_ion_pflux_out(2)*EXPRO_rmin(nr)/(tglf_as_in(3)*EXPRO_ne(ir)*tglf_rlns_in(3))
  chi_eff = chi_gB(ir)*(tglf_elec_eflux_out+tglf_ion_eflux_out(1))*EXPRO_rmin(nr) / &
             (0.5*EXPRO_ne(ir)*( tglf_rlts_in(1)+tglf_as_in(2)*tglf_taus_in(2)*tglf_rlts_in(3) ))
  QL_ratio_loc = DEP_trace_loc / (chi_eff+1.0E-8)
  print *, 'ir =', ir, ' DEP =', DEP_trace_loc, 'QL_ratio =', QL_ratio_loc 

end subroutine TGLFEP_EPflux
