!---------------------------------------------------------
! TGLFEP_interface.f90
! 7.13.2016
! PURPOSE:
!  The parameters for TGLFEP
!---------------------------------------------------------

module TGLFEP_interface
  
  implicit none

  integer :: TGLFEP_COMM
  integer :: TGLFEP_COMM_WORLD

  integer,parameter :: nmodes=4

  integer :: process_in
  integer :: input_profile_method=1      ! EMB 10/2018
  integer :: threshold_flag
  integer :: scan_method

  integer :: mode_in
  logical :: width_in_flag
  real :: width_in, width_min, width_max
  real, allocatable, dimension(:) :: factor_max_profile
  real :: factor_in, factor_max, ky_in, kymark

  real :: DEP_trace_loc
  real, allocatable, dimension(:) :: DEP_trace
  real, allocatable, dimension(:) :: DEP_trace_complete
  real, allocatable, dimension(:,:) :: DEP_trace_scan
  real :: QL_ratio_loc
  real, allocatable, dimension(:) :: QL_ratio
  real, allocatable, dimension(:) :: QL_ratio_complete
  real, allocatable, dimension(:,:) :: QL_ratio_scan
  real, allocatable, dimension(:) :: chi_gB
  integer, allocatable, dimension(:) :: ir_exp

  integer :: n_basis
  integer :: ir, n_toroidal, ky_model
  real :: kyhat_in
  real :: scan_factor = 1.0
  integer :: is_EP
  integer :: jtscale
  integer, parameter :: jtscale_max=1
  real, parameter :: tscale_interval=1.0
  integer :: TGLFEP_nion
  real, dimension(10) :: m_i  = 1.0
!  real :: m_EP = 1.0
  real, dimension(10) :: q_i  = 1.0
!  real :: q_EP = 1.0

  real,parameter :: freq_cutoff = -0.2
  real :: freq_AE_upper

  integer :: nn
  real,allocatable,dimension(:) :: factor_out

  integer :: id_2, np_2, id_3, np_3

  character(16) :: suffix = ''
  character(50) :: str_wf_file

  logical :: l_print
  integer :: l_write_wavefunction
  integer :: l_wavefunction_out

  integer :: reject_i_pinch_flag
  integer :: reject_e_pinch_flag
  integer :: reject_th_pinch_flag
  integer :: reject_EP_pinch_flag
  integer :: reject_tearing_flag
  integer :: reject_max_outer_panel_flag
  real :: QL_ratio_thresh
  real :: q_scale

  logical, dimension(4) :: lkeep
  logical, dimension(4) :: ltearing
  logical, dimension(4) :: l_th_pinch
  logical, dimension(4) :: l_i_pinch
  logical, dimension(4) :: l_e_pinch
  logical, dimension(4) :: l_EP_pinch
  logical, dimension(4) :: l_QL_ratio
  logical, dimension(4) :: l_max_outer_panel

  integer :: l_real_units = 0
  real, allocatable, dimension(:) :: f_real

end module TGLFEP_interface
