      SUBROUTINE tglf_setup_geometry
!***********************************************************
!  This routine is used to compute the geometric quantities
!  needed by a transport code. It does not compute the linear
!  eigenmodes. As a side effect the species and basis functions
!  are defined so be sure that put_species and put_switches
!  have been called if you don't want to use the defaults.
!  R2_ave_out
!  B2_ave_out
!  others to be added later
!***********************************************************
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
!
      new_geometry = .TRUE.
!
!
! set global constants
!
      pi = atan2( 0.0, -1.0 )
      pi_2 = 2.0*pi
      sqrt_pi = SQRT(pi)
      sqrt_two = SQRT(2.0)
!
      new_width=.TRUE.
!
!   load the xgrid and gauss-hermite weights
!
       if(gauher_uncalled)CALL gauher
!
! set ky in units of k_theta*rho_s  rho_s=C_s/omega_s
! C_s=sqrt(Te/mi), omega_s=eB/(mi c)
!
     ky = ky_s
!      write(*,*)"ky = ",ky
!
! check co-dependencies 
!
      if(new_geometry)new_width=.TRUE.
      if(new_width)new_matrix=.TRUE.
!      write(*,*)"new_start=",new_start
!      write(*,*)"new_geometry=",new_geometry
!      write(*,*)"new_width=",new_width
!      write(*,*)"new_matrix=",new_matrix
!      write(*,*)"width_in=",width_in
!      write(*,*)"nbasis_max_in=",nbasis_max_in
!      write(*,*)"new_eikonal_in=",new_eikonal_in
!      write(*,*)"betae_in=",betae_in
!      write(*,*)"debye_in=",debye_in
!      write(*,*)"xnuei_in=",xnuei_in
!      write(*,*)"zeff_in=",zeff_in
!      write(*,*)"filter_in=",filter_in
!      write(*,*)"park_in=",park_in
!      write(*,*)"ghat_in=",ghat_in
!      write(*,*)"wd_zero_in=",wd_zero_in
!      write(*,*)"Linsker_factor_in=",Linsker_factor_in
!      write(*,*)"gradB_factor_in=",gradB_factor_in
!      write(*,*)"x_psi_in=",x_psi_in
!      write(*,*)"xnu_factor_in=",xnu_factor_in
!      write(*,*)"debye_factor_in=",debye_factor_in
!      write(*,*)"theta_trapped_in=",theta_trapped_in
!      write(*,*)"alpha_p_in=",alpha_p_in
!      write(*,*)"alpha_e_in=",alpha_e_in
!
!
       if(igeo.eq.1)then
! set up MILLER geometry
         call miller_geo
       elseif(igeo.eq.2)then
! set up Fourier geometry
         call fourier_geo
       elseif(igeo.eq.3)then
! set up ELITE geometry
         call ELITE_geo
       endif
!
! compute the eikonal functions for general geometry (igeo > 0)
!
       if(igeo.gt.0)call mercier_luc 
!
       call get_xgrid_functions
!
       new_geometry=.FALSE.
!
      END SUBROUTINE tglf_setup_geometry
