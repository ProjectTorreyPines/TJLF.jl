!------------------------------------------------------------ TGLFEP_ky.f90
!
! PURPOSE: Calculate the growth rate and frequency at a single ky
!------------------------------------------------------------

subroutine TGLFEP_ky

  use tglf_interface !module within TGLF
  use tglf_pkg !module within TGLF
  use TGLFEP_interface !module within TGLFEP
  use tglf_max_dimensions !module within TGLF

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
  real :: th_QL_flux
  real :: i_QL_flux
  real, dimension(4) :: i_QL_cond_flux
  real :: e_QL_flux
  real, dimension(4) :: e_QL_cond_flux 
  real :: th_eff_grad
  real :: i_eff_grad
  real, dimension(4) :: QL_flux_ratio

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

  tglf_kx0_loc_in = 0.0

!  printing out the tglf struct:
  if (ir .eq. 102) then ! SET BACK TO A VALID IR TO CHECK VALUES
!    open(1000, file='input.tglf', status='replace')
    print *, "tglf_units_in", tglf_units_in
    print *, "tglf_use_bper_in", tglf_use_bper_in
    print *, "tglf_use_bpar_in", tglf_use_bpar_in
    print *, "tglf_use_mhd_rule_in", tglf_use_mhd_rule_in
    print *, "tglf_use_bisection_in", tglf_use_bisection_in
    print *, "tglf_use_inboard_detrapped_in", tglf_use_inboard_detrapped_in
    print *, "tglf_use_ave_ion_grid_in", tglf_use_ave_ion_grid_in
    print *, "tglf_new_eikonal_in", tglf_new_eikonal_in
    print *, "tglf_find_width_in", tglf_find_width_in
    print *, "tglf_iflux_in", tglf_iflux_in
    print *, "tglf_adiabatic_elec_in", tglf_adiabatic_elec_in

    print *, "tglf_sat_rule_in", tglf_sat_rule_in
    print *, "tglf_ns_in", tglf_ns_in
    print *, "tglf_nmodes_in", tglf_nmodes_in
    print *, "tglf_nwidth_in", tglf_nwidth_in 
    print *, "tglf_nbasis_max_in", tglf_nbasis_max_in
    print *, "tglf_nbasis_min_in", tglf_nbasis_min_in
    print *, "tglf_nxgrid_in", tglf_nxgrid_in
    print *, "tglf_nky_in", tglf_nky_in
    print *, "tglf_kygrid_model_in", tglf_kygrid_model_in
    print *, "tglf_xnu_model_in", tglf_xnu_model_in
    print *, "tglf_vpar_model_in", tglf_vpar_model_in
    print *, "tglf_ibranch_in", tglf_ibranch_in

    print *, "tglf_zs_in(1)", tglf_zs_in(1)
    print *, "tglf_zs_in(2)", tglf_zs_in(2)
    print *, "tglf_zs_in(3)", tglf_zs_in(3)
    print *, "tglf_zs_in(4)", tglf_zs_in(4)
!    print *, tglf_zs_in(5)
!    print *, tglf_zs_in(6)
!    print *, tglf_zs_in(7)

    print *, "tglf_mass_in(1)", tglf_mass_in(1)
    print *, "tglf_mass_in(2)", tglf_mass_in(2)
    print *, "tglf_mass_in(3)", tglf_mass_in(3)
    print *, "tglf_mass_in(4)", tglf_mass_in(4)
    !    print *, tglf_mass_in(5)
    !    print *, tglf_mass_in(6)
    !    print *, tglf_mass_in(7)
    
    print *, "tglf_rlns_in(1)", tglf_rlns_in(1)
    print *, "tglf_rlns_in(2)", tglf_rlns_in(2)
    print *, "tglf_rlns_in(3)", tglf_rlns_in(3)
    print *, "tglf_rlns_in(4)", tglf_rlns_in(4)
    !    print *, tglf_rlns_in(5)
    !    print *, tglf_rlns_in(6)
    !    print *, tglf_rlns_in(7)
    
    print *, "tglf_rlts_in(1)", tglf_rlts_in(1)
    print *, "tglf_rlts_in(2)", tglf_rlts_in(2)
    print *, "tglf_rlts_in(3)", tglf_rlts_in(3)
    print *, "tglf_rlts_in(4)", tglf_rlts_in(4)
    !    print *, tglf_rlts_in(5)
    !    print *, tglf_rlts_in(6)
    !    print *, tglf_rlts_in(7)
    
    print *, "tglf_taus_in(1)", tglf_taus_in(1)
    print *, "tglf_taus_in(2)", tglf_taus_in(2)
    print *, "tglf_taus_in(3)", tglf_taus_in(3)
    print *, "tglf_taus_in(4)", tglf_taus_in(4)
    !    print *, tglf_taus_in(5)
    !    print *, tglf_taus_in(6)
    !    print *, tglf_taus_in(7)
    
    print *, "tglf_as_in(1)", tglf_as_in(1)
    print *, "tglf_as_in(2)", tglf_as_in(2)
    print *, "tglf_as_in(3)", tglf_as_in(3)
    print *, "tglf_as_in(4)", tglf_as_in(4)
    !    print *, tglf_as_in(5)
    !    print *, tglf_as_in(6)
    !    print *, tglf_as_in(7)
    
    print *, "tglf_vpar_in(1)", tglf_vpar_in(1)
    print *, "tglf_vpar_in(2)", tglf_vpar_in(2)
    print *, "tglf_vpar_in(3)", tglf_vpar_in(3)
    print *, "tglf_vpar_in(4)", tglf_vpar_in(4)
    !    print *, tglf_vpar_in(5)
    !    print *, tglf_vpar_in(6)
    !    print *, tglf_vpar_in(7)
    
    print *, "tglf_vpar_shear_in(1)", tglf_vpar_shear_in(1)
    print *, "tglf_vpar_shear_in(2)", tglf_vpar_shear_in(2)
    print *, "tglf_vpar_shear_in(3)", tglf_vpar_shear_in(3)
    print *, "tglf_vpar_shear_in(4)", tglf_vpar_shear_in(4)
    !    print *, tglf_vpar_shear_in(5)
    !    print *, tglf_vpar_shear_in(6)
    !    print *, tglf_vpar_shear_in(7)
    
    print *, "tglf_sign_bt_in", tglf_sign_bt_in
    print *, "tglf_sign_it_in", tglf_sign_it_in
    print *, "tglf_ky_in", tglf_ky_in
    
    print *, "tglf_vexb_shear_in", tglf_vexb_shear_in
    print *, "tglf_betae_in", tglf_betae_in
    print *, "tglf_xnue_in", tglf_xnue_in
    print *, "tglf_zeff_in", tglf_zeff_in
    print *, "tglf_debye_in", tglf_debye_in
    
    print *, "tglf_alpha_mach_in", tglf_alpha_mach_in
    print *, "tglf_alpha_e_in", tglf_alpha_e_in
    print *, "tglf_alpha_p_in", tglf_alpha_p_in
    print *, "tglf_alpha_quench_in", tglf_alpha_quench_in
    print *, "tglf_alpha_zf_in", tglf_alpha_zf_in
    print *, "tglf_xnu_factor_in", tglf_xnu_factor_in
    print *, "tglf_debye_factor_in", tglf_debye_factor_in
    print *, "tglf_etg_factor_in", tglf_etg_factor_in
    print *, "tglf_rlnp_cutoff_in", tglf_rlnp_cutoff_in
    
    print *, "tglf_width_in", tglf_width_in
    print *, "tglf_width_min_in", tglf_width_min_in
    
    print *, "tglf_rmin_loc_in", tglf_rmin_loc_in
    print *, "tglf_rmaj_loc_in", tglf_rmaj_loc_in
    print *, "tglf_zmaj_loc_in", tglf_zmaj_loc_in
    print *, "tglf_drmindx_loc_in", tglf_drmindx_loc_in
    print *, "tglf_drmajdx_loc_in", tglf_drmajdx_loc_in
    print *, "tglf_dzmajdx_loc_in", tglf_dzmajdx_loc_in
    print *, "tglf_q_loc_in", tglf_q_loc_in
    print *, "tglf_kappa_loc_in", tglf_kappa_loc_in
    print *, "tglf_s_kappa_loc_in", tglf_s_kappa_loc_in
    print *, "tglf_delta_loc_in", tglf_delta_loc_in
    print *, "tglf_s_delta_loc_in", tglf_s_delta_loc_in
    print *, "tglf_zeta_loc_in", tglf_zeta_loc_in
    print *, "tglf_s_zeta_loc_in", tglf_s_zeta_loc_in
    print *, "tglf_p_prime_loc_in", tglf_p_prime_loc_in
    print *, "tglf_q_prime_loc_in", tglf_q_prime_loc_in
    print *, "tglf_beta_loc_in", tglf_beta_loc_in
    print *, "tglf_kx0_loc_in", tglf_kx0_loc_in
    
    print *, "tglf_damp_psi_in", tglf_damp_psi_in
    print *, "tglf_damp_sig_in", tglf_damp_sig_in
    print *, "tglf_wdia_trapped_in", tglf_wdia_trapped_in
    print *, "tglf_park_in", tglf_park_in
    print *, "tglf_ghat_in", tglf_ghat_in
    print *, "tglf_gchat_in", tglf_gchat_in
    print *, "tglf_wd_zero_in", tglf_wd_zero_in
    print *, "tglf_linsker_factor_in", tglf_linsker_factor_in
    print *, "tglf_gradb_factor_in", tglf_gradb_factor_in
    print *, "tglf_filter_in", tglf_filter_in
    print *, "tglf_theta_trapped_in", tglf_theta_trapped_in
    print *, "tglf_use_transport_model_in", tglf_use_transport_model_in
    print *, "------------------"
  endif


  if (ir .eq. 3) then
    !print *, tglf_ibranch_in
    !print *, tglf_width_in
    !print *, tglf_as_in(1)
    !print *, tglf_as_in(2)
    !print *, tglf_as_in(3)
    !print *, tglf_as_in(4)
    
    print *, tglf_ky_in
    print *, tglf_nky_in
    print *, tglf_kygrid_model_in
    print *, " "
  endif
!  print *, 'In TGLFEP_ky (2) factor_in=', factor_in
  call tglf_run
  !print *, "tglf_run called: TGLFEP_ky.f90"
!  print *, 'In TGLFEP_ky (3) factor_in=', factor_in

  do n=1,nmodes
    g(n) = get_growthrate(n)
    f(n) = get_frequency(n)
  enddo

  if (ir .eq. 101) then
    do n = 1, nmodes
      print *, g(n)
    end do
    print *, "------------------"
    do n = 1, nmodes
      print *, f(n)
    enddo
    print *, "=================="

  endif

  do n = 1, nmodes
    lkeep(n) = f(n) .lt. freq_AE_upper
    lkeep(n) = lkeep(n) .and. g(n) .gt. 1.0E-7
  enddo

  if (ir .eq. 101) then ! SET BACK TO A VALID IR TO CHECK VALUES
    !    open(1000, file='input.tglf', status='replace')
        print *, "tglf_units_in", tglf_units_in
        print *, "tglf_use_bper_in", tglf_use_bper_in
        print *, "tglf_use_bpar_in", tglf_use_bpar_in
        print *, "tglf_use_mhd_rule_in", tglf_use_mhd_rule_in
        print *, "tglf_use_bisection_in", tglf_use_bisection_in
        print *, "tglf_use_inboard_detrapped_in", tglf_use_inboard_detrapped_in
        print *, "tglf_use_ave_ion_grid_in", tglf_use_ave_ion_grid_in
        print *, "tglf_new_eikonal_in", tglf_new_eikonal_in
        print *, "tglf_find_width_in", tglf_find_width_in
        print *, "tglf_iflux_in", tglf_iflux_in
        print *, "tglf_adiabatic_elec_in", tglf_adiabatic_elec_in
    
        print *, "tglf_sat_rule_in", tglf_sat_rule_in
        print *, "tglf_ns_in", tglf_ns_in
        print *, "tglf_nmodes_in", tglf_nmodes_in
        print *, "tglf_nwidth_in", tglf_nwidth_in 
        print *, "tglf_nbasis_max_in", tglf_nbasis_max_in
        print *, "tglf_nbasis_min_in", tglf_nbasis_min_in
        print *, "tglf_nxgrid_in", tglf_nxgrid_in
        print *, "tglf_nky_in", tglf_nky_in
        print *, "tglf_kygrid_model_in", tglf_kygrid_model_in
        print *, "tglf_xnu_model_in", tglf_xnu_model_in
        print *, "tglf_vpar_model_in", tglf_vpar_model_in
        print *, "tglf_ibranch_in", tglf_ibranch_in
    
        print *, "tglf_zs_in(1)", tglf_zs_in(1)
        print *, "tglf_zs_in(2)", tglf_zs_in(2)
        print *, "tglf_zs_in(3)", tglf_zs_in(3)
        print *, "tglf_zs_in(4)", tglf_zs_in(4)
    !    print *, tglf_zs_in(5)
    !    print *, tglf_zs_in(6)
    !    print *, tglf_zs_in(7)
    
        print *, "tglf_mass_in(1)", tglf_mass_in(1)
        print *, "tglf_mass_in(2)", tglf_mass_in(2)
        print *, "tglf_mass_in(3)", tglf_mass_in(3)
        print *, "tglf_mass_in(4)", tglf_mass_in(4)
        !    print *, tglf_mass_in(5)
        !    print *, tglf_mass_in(6)
        !    print *, tglf_mass_in(7)
        
        print *, "tglf_rlns_in(1)", tglf_rlns_in(1)
        print *, "tglf_rlns_in(2)", tglf_rlns_in(2)
        print *, "tglf_rlns_in(3)", tglf_rlns_in(3)
        print *, "tglf_rlns_in(4)", tglf_rlns_in(4)
        !    print *, tglf_rlns_in(5)
        !    print *, tglf_rlns_in(6)
        !    print *, tglf_rlns_in(7)
        
        print *, "tglf_rlts_in(1)", tglf_rlts_in(1)
        print *, "tglf_rlts_in(2)", tglf_rlts_in(2)
        print *, "tglf_rlts_in(3)", tglf_rlts_in(3)
        print *, "tglf_rlts_in(4)", tglf_rlts_in(4)
        !    print *, tglf_rlts_in(5)
        !    print *, tglf_rlts_in(6)
        !    print *, tglf_rlts_in(7)
        
        print *, "tglf_taus_in(1)", tglf_taus_in(1)
        print *, "tglf_taus_in(2)", tglf_taus_in(2)
        print *, "tglf_taus_in(3)", tglf_taus_in(3)
        print *, "tglf_taus_in(4)", tglf_taus_in(4)
        !    print *, tglf_taus_in(5)
        !    print *, tglf_taus_in(6)
        !    print *, tglf_taus_in(7)
        
        print *, "tglf_as_in(1)", tglf_as_in(1)
        print *, "tglf_as_in(2)", tglf_as_in(2)
        print *, "tglf_as_in(3)", tglf_as_in(3)
        print *, "tglf_as_in(4)", tglf_as_in(4)
        !    print *, tglf_as_in(5)
        !    print *, tglf_as_in(6)
        !    print *, tglf_as_in(7)
        
        print *, "tglf_vpar_in(1)", tglf_vpar_in(1)
        print *, "tglf_vpar_in(2)", tglf_vpar_in(2)
        print *, "tglf_vpar_in(3)", tglf_vpar_in(3)
        print *, "tglf_vpar_in(4)", tglf_vpar_in(4)
        !    print *, tglf_vpar_in(5)
        !    print *, tglf_vpar_in(6)
        !    print *, tglf_vpar_in(7)
        
        print *, "tglf_vpar_shear_in(1)", tglf_vpar_shear_in(1)
        print *, "tglf_vpar_shear_in(2)", tglf_vpar_shear_in(2)
        print *, "tglf_vpar_shear_in(3)", tglf_vpar_shear_in(3)
        print *, "tglf_vpar_shear_in(4)", tglf_vpar_shear_in(4)
        !    print *, tglf_vpar_shear_in(5)
        !    print *, tglf_vpar_shear_in(6)
        !    print *, tglf_vpar_shear_in(7)
        
        print *, "tglf_sign_bt_in", tglf_sign_bt_in
        print *, "tglf_sign_it_in", tglf_sign_it_in
        print *, "tglf_ky_in", tglf_ky_in
        
        print *, "tglf_vexb_shear_in", tglf_vexb_shear_in
        print *, "tglf_betae_in", tglf_betae_in
        print *, "tglf_xnue_in", tglf_xnue_in
        print *, "tglf_zeff_in", tglf_zeff_in
        print *, "tglf_debye_in", tglf_debye_in
        
        print *, "tglf_alpha_mach_in", tglf_alpha_mach_in
        print *, "tglf_alpha_e_in", tglf_alpha_e_in
        print *, "tglf_alpha_p_in", tglf_alpha_p_in
        print *, "tglf_alpha_quench_in", tglf_alpha_quench_in
        print *, "tglf_alpha_zf_in", tglf_alpha_zf_in
        print *, "tglf_xnu_factor_in", tglf_xnu_factor_in
        print *, "tglf_debye_factor_in", tglf_debye_factor_in
        print *, "tglf_etg_factor_in", tglf_etg_factor_in
        print *, "tglf_rlnp_cutoff_in", tglf_rlnp_cutoff_in
        
        print *, "tglf_width_in", tglf_width_in
        print *, "tglf_width_min_in", tglf_width_min_in
        
        print *, "tglf_rmin_loc_in", tglf_rmin_loc_in
        print *, "tglf_rmaj_loc_in", tglf_rmaj_loc_in
        print *, "tglf_zmaj_loc_in", tglf_zmaj_loc_in
        print *, "tglf_drmindx_loc_in", tglf_drmindx_loc_in
        print *, "tglf_drmajdx_loc_in", tglf_drmajdx_loc_in
        print *, "tglf_dzmajdx_loc_in", tglf_dzmajdx_loc_in
        print *, "tglf_q_loc_in", tglf_q_loc_in
        print *, "tglf_kappa_loc_in", tglf_kappa_loc_in
        print *, "tglf_s_kappa_loc_in", tglf_s_kappa_loc_in
        print *, "tglf_delta_loc_in", tglf_delta_loc_in
        print *, "tglf_s_delta_loc_in", tglf_s_delta_loc_in
        print *, "tglf_zeta_loc_in", tglf_zeta_loc_in
        print *, "tglf_s_zeta_loc_in", tglf_s_zeta_loc_in
        print *, "tglf_p_prime_loc_in", tglf_p_prime_loc_in
        print *, "tglf_q_prime_loc_in", tglf_q_prime_loc_in
        print *, "tglf_beta_loc_in", tglf_beta_loc_in
        print *, "tglf_kx0_loc_in", tglf_kx0_loc_in
        
        print *, "tglf_damp_psi_in", tglf_damp_psi_in
        print *, "tglf_damp_sig_in", tglf_damp_sig_in
        print *, "tglf_wdia_trapped_in", tglf_wdia_trapped_in
        print *, "tglf_park_in", tglf_park_in
        print *, "tglf_ghat_in", tglf_ghat_in
        print *, "tglf_gchat_in", tglf_gchat_in
        print *, "tglf_wd_zero_in", tglf_wd_zero_in
        print *, "tglf_linsker_factor_in", tglf_linsker_factor_in
        print *, "tglf_gradb_factor_in", tglf_gradb_factor_in
        print *, "tglf_filter_in", tglf_filter_in
        print *, "tglf_theta_trapped_in", tglf_theta_trapped_in
        print *, "tglf_use_transport_model_in", tglf_use_transport_model_in
        print *, "------------------"

        !print *, "satParams y: ", y
        !print *, "satParams t_s: ", t_s
        !print *, "------------------"
        print *, "nmodes_out: ", nmodes_out
        print *, "nbasis: ", n_basis
        print *, "=================="


  endif

  call get_wavefunction_out(nmodes_out,nfields_out,max_plot_out,angle_out,wavefunction_out)
!  print *, 'In TGLFEP_ky (4) factor_in=', factor_in
  ltearing(:) = .false.
  l_i_pinch(:) = .false.
  l_e_pinch(:) = .false.
  l_EP_pinch(:) = .false.
  l_th_pinch(:) = .false.
  l_QL_ratio(:) = .false.
  l_max_outer_panel(:) = .false.
  if (ir .eq. 101) then
    print *, "======="
  endif
  do n = 1, nmodes
    x_tear_test(n) = 0.0
    
    wave_max=maxval(abs(wavefunction_out(n,1,:)))+1.0E-3
    wave_max_loc=maxloc(abs(wavefunction_out(n,1,:)))
!    shouldn't wavefunction be 0 if growthrate is 0?
!    if (g(n) .eq. 0.0) then
!      wave_max=0.001
!      wave_max_loc=1
!    endif
    n_balloon_pi = (max_plot-1)/9
    i_mid_plot = (max_plot-1)/2 + 1
    l_max_outer_panel(n) = (wave_max_loc(1) .lt. (i_mid_plot-n_balloon_pi)) .or. (wave_max_loc(1) .gt. (i_mid_plot+n_balloon_pi))

    if (ir .eq. 101 .and. n .eq. 4) then
      print *, "wave_max: ", wave_max
      print *, "wave_max_loc: ", wave_max_loc(1)
      print *, "mode: ", n
      print *, "l_max_outer_panel(n): ", l_max_outer_panel(n)
      print *, "-------"
      print *, "wavefunction: "
      do i = wave_max_loc(1), wave_max_loc(1) + 10
          print '(F12.6, " ")', abs(wavefunction_out(n, 1, i))
      end do
    endif

    do i = 1, max_plot
       x_tear_test(n) = max(x_tear_test(n),abs(wavefunction_out(n,1,i)-wavefunction_out(n,1,max_plot+1-i))/wave_max)
    enddo
    if (x_tear_test(n) .GT. 1.0E-1) ltearing(n) = .true.
    EP_QL_flux = 0.0    ! Particle flux for EPs
    i_QL_flux = 0.0     ! Energy flux for e- and ions
    i_QL_cond_flux(n) = 0.0
    i_eff_grad = 0.0
    e_QL_flux = 0.0
    e_QL_cond_flux(n) = 0.0
    th_QL_flux = 0.0
    th_eff_grad = 0.0
    do jfields = 1, 3
      EP_QL_flux = EP_QL_flux + get_QL_particle_flux(n,is_EP+1,jfields)
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
    QL_flux_ratio(n) = (EP_QL_flux/tglf_as_in(is_EP+1))/(abs(i_QL_cond_flux(n))/(tglf_as_in(1)-tglf_as_in(is_EP+1)))
!    QL_flux_ratio(n) = DEP(n)/chi_i_cond(n)
    if (chi_i(n) .lt. 0.0) l_i_pinch(n) = .true.
    if (chi_e(n) .lt. 0.0) l_e_pinch(n) = .true.
    if (chi_th(n) .lt. 0.0) l_th_pinch(n) = .true.
    if (DEP(n) .lt. 0.0) l_EP_pinch(n) = .true.
    if (QL_flux_ratio(n) .lt. QL_ratio_thresh) then
      l_QL_ratio = .true.
    endif
  enddo

  do n = 1, nmodes
    if (reject_tearing_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. ltearing(n)
    if (reject_i_pinch_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_i_pinch(n)
    if (reject_e_pinch_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_e_pinch(n)
    if (reject_th_pinch_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_th_pinch(n)
    if (reject_EP_pinch_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_EP_pinch(n)
    if (reject_max_outer_panel_flag .eq. 1) lkeep(n) = lkeep(n) .and. .not. l_max_outer_panel(n)
    lkeep(n) = lkeep(n) .and. .not. l_QL_ratio(n)
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
    write (22,*) lkeep(:)
! First, renormalize and adjust phases
    do n = 1, nmodes_out
      max_phi = maxval(abs(wavefunction_out(n,1,:)))
      max_apar = maxval(abs(wavefunction_out(n,2,:)))
      max_field = max(max_phi,max_apar)
      phase = atan2(aimag(wavefunction_out(n,1,(max_plot_out+1)/2)),real(wavefunction_out(n,1,(max_plot_out+1)/2)))
      do jfields = 1, nfields_out
        wavefunction_out(n,jfields,:) = wavefunction_out(n,jfields,:)/(max_field*exp(cmplx(0,1)*phase))
      enddo
    enddo
! Write renormalized, re-phased eigenfunctions out to str_wf_file.
    do i = 1, max_plot
      kcomp = 0
      wave_write(:) = 0.0
      do n = 1, nmodes_out
        do jfields = 1, nfields_out
          kcomp = kcomp + 1
          wave_write(kcomp) = real(wavefunction_out(n,jfields,i))
          kcomp = kcomp + 1
          wave_write(kcomp) = aimag(wavefunction_out(n,jfields,i))
        enddo
      enddo
      write (22,'(F10.6,4E16.6E3)') angle_out(i), real(wavefunction_out(1,1,i)), aimag(wavefunction_out(1,1,i)),&
                                   real(wavefunction_out(1,2,i)), aimag(wavefunction_out(1,2,i)) 
                ! Just leading mode with two fields hardcoded for now
    enddo
    close (unit=22)
  endif 
end subroutine TGLFEP_ky
