!-----------------------------------------------------------------
!
      SUBROUTINE tglf_TM_mpi
!
!  Main transport model subroutine.
!  Calls linear TGLF over a spectrum of ky's and computes spectral integrals of 
!  field, intensity and fluxes.
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_kyspectrum
      IMPLICIT NONE
!
      INTEGER :: i,j,is,imax, jmax_mix
      LOGICAL :: save_iflux, save_find_width
      REAL :: dky
      REAL :: phi_bar0,phi_bar1
      REAl :: v_bar0,v_bar1
      REAL :: dky0,dky1,ky0,ky1
      REAL :: pflux0(nsm,3),eflux0(nsm,3)
      REAL :: stress_par0(nsm,3),stress_tor0(nsm,3)
      REAL :: exch0(nsm,3)
      REAL :: nsum0(nsm),tsum0(nsm)
      REAL :: pflux1(nsm,3),eflux1(nsm,3)
      REAL :: stress_par1(nsm,3),stress_tor1(nsm,3)
      REAL :: exch1(nsm,3)
      REAL :: nsum1(nsm),tsum1(nsm)
      REAL :: save_vexb_shear
!
      CALL tglf_startup
!
! initialize fluxes
!
      do is=1,ns
        do j=1,3 
          particle_flux_out(is,j) = 0.0
          energy_flux_out(is,j) = 0.0
          stress_par_out(is,j) = 0.0
          stress_tor_out(is,j) = 0.0
          exchange_out(is,j) = 0.0
          pflux0(is,j) = 0.0
          eflux0(is,j) = 0.0
          stress_par0(is,j) = 0.0
          stress_tor0(is,j) = 0.0
          exch0(is,j) = 0.0
        enddo
        n_bar_sum_out(is) = 0.0
        t_bar_sum_out(is) = 0.0
        q_low_out(is) = 0.0
        nsum0(is) = 0.0
        tsum0(is) = 0.0
      enddo
      phi_bar_sum_out = 0.0
      v_bar_sum_out = 0.0
      v_bar0 = 0.0
      phi_bar0 = 0.0
      do i=1,nky
       width_out(i) = width_in ! needed for spectral shift model double pass
      enddo
!
! compute the flux spectrum
!
    if(alpha_quench_in .eq. 0.0 .and. vexb_shear_s .ne.0.0)then
!  spectral shift model double pass
       save_vexb_shear = vexb_shear_s
       save_find_width = find_width_in
       save_iflux = iflux_in
       iflux_in=.FALSE.     ! do not compute eigenvectors on first pass
       vexb_shear_s = 0.0  ! do not use spectral shift on first pass
       jmax_out = 0         ! ID for first pass
       CALL get_bilinear_spectrum_mpi
       eigenvalue_first_pass(:,:,:) = eigenvalue_spectrum_out(:,:,:)
       vexb_shear_s = save_vexb_shear
       find_width_in = .FALSE.
       iflux_in = save_iflux
       if(sat_rule_in.eq.0)jmax_out = 1   ! flag for second pass
       CALL get_bilinear_spectrum_mpi
!  reset eigenvalues to the values with vexb_shear=0.
!  note ql weights are with vexb_shear
       eigenvalue_spectrum_out(:,:,:) = eigenvalue_first_pass(:,:,:)
       find_width_in = save_find_width
      else
        ! only one pass for alpha_quench_in .ne. 0.0 or vexb_shear_s .eq. 0.0
        jmax_out = 0
        CALL get_bilinear_spectrum_mpi
      endif
!
! sum over ky spectrum
!
      iflux_in=.TRUE. 
      dky0=0.0
      ky0=0.0 
      do i=1,nky
        ky_s = ky_spectrum(i)
        dky = dky_spectrum(i)
        ky1=ky_s
        if(i.eq.1)then
          dky1=ky1
        elseif(kygrid_model_in.ne.0)then
            dky = LOG(ky1/ky0)/(ky1-ky0)
            dky1 = ky1*(1.0 - ky0*dky)
            dky0 = ky0*(ky1*dky - 1.0)
        endif
! backward comatiblily for UNITS=GYRO
!        if(units_in.eq.'GYRO')then
! this is an error in the original integrator
!           dky0 = dky0*SQRT(taus_in(1)*mass_in(2))
!           dky1 = dky1*SQRT(taus_in(1)*mass_in(2))
!        endif
!
! compute the field integrals
!
        v_bar1 = 0.0
        phi_bar1 = 0.0
        do imax = 1,nmodes_in
          v_bar1 = v_bar1 + field_spectrum_out(1,i,imax)
          phi_bar1 = phi_bar1 + field_spectrum_out(2,i,imax)
        enddo
        phi_bar_sum_out = phi_bar_sum_out + dky0*phi_bar0 + dky1*phi_bar1
        v_bar_sum_out = v_bar_sum_out + dky0*v_bar0 + dky1*v_bar1
        phi_bar0 = phi_bar1
        v_bar0 = v_bar1
!
! compute the intensity integrals
!
        do is=ns0,ns
          nsum1(is) = 0.0
          tsum1(is) = 0.0
          do imax = 1,nmodes_in
            nsum1(is) = nsum1(is) + intensity_spectrum_out(1,is,i,imax)
            tsum1(is) = tsum1(is) + intensity_spectrum_out(2,is,i,imax)
          enddo
          n_bar_sum_out(is) = n_bar_sum_out(is) &
              + dky0*nsum0(is) + dky1*nsum1(is)
          t_bar_sum_out(is) = t_bar_sum_out(is) &
              + dky0*tsum0(is) + dky1*tsum1(is)
          nsum0(is) = nsum1(is)
          tsum0(is) = tsum1(is)
        enddo
!
! compute the flux integrals
!
        do is=ns0,ns
          do j=1,3
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
            particle_flux_out(is,j) = particle_flux_out(is,j) &
              + dky0*pflux0(is,j) + dky1*pflux1(is,j)
            energy_flux_out(is,j) = energy_flux_out(is,j) &
              + dky0*eflux0(is,j) + dky1*eflux1(is,j)
            stress_tor_out(is,j) = stress_tor_out(is,j) &
              + dky0*stress_tor0(is,j) + dky1*stress_tor1(is,j)
            stress_par_out(is,j) = stress_par_out(is,j) &
              + dky0*stress_par0(is,j) + dky1*stress_par1(is,j)
            exchange_out(is,j) = exchange_out(is,j) &
              + dky0*exch0(is,j) + dky1*exch1(is,j)
!            write(*,*)is,j,i
!            write(*,*)"ky0=",ky0,"ky1=",ky1
!            write(*,*)"pflux0=",pflux0,"pflux1=",pflux1
!            write(*,*)"eflux0=",eflux0,"eflux1=",eflux1
!            write(*,*)dky0*pflux0+dky1*pflux1
!            write(*,*)dky0*eflux0+dky1*eflux1
!            write(*,*)"stress_tor_out=",stress_tor_out(is,1)
             pflux0(is,j) = pflux1(is,j)
             eflux0(is,j) = eflux1(is,j)
             stress_par0(is,j) = stress_par1(is,j)
             stress_tor0(is,j) = stress_tor1(is,j)
             exch0(is,j) = exch1(is,j)
           enddo  ! j
           if(ky_s.le.1.0)then
             q_low_out(is) = energy_flux_out(is,1)+energy_flux_out(is,2)
           endif
         enddo  ! is 
!
        ky0 = ky1
      enddo  ! i
!
      CALL tglf_shutdown
!
      END SUBROUTINE tglf_TM_mpi
!
!-----------------------------------------------------------------
!
      SUBROUTINE get_bilinear_spectrum_mpi
!
! computes the bilinear fluctuation moments 
! and saves them in flux_spectrum_out, intensity_spectrum_out
! and field_spectrum_out
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_mpi
      USE tglf_species
      USE tglf_kyspectrum
      USE tglf_xgrid
      IMPLICIT NONE
!
      LOGICAL :: unstable
      INTEGER :: i,j,k,is,imax,t
      REAL :: save_width
      REAL :: gmax,fmax
      REAL :: phi_bar
      REAL :: gamma_cutoff,reduce,rexp
      REAL :: gamma_net_1
      REAL :: pflux1,eflux1
      REAL :: stress_tor1,stress_par1
      REAL :: exch1, gamma_max
      ! mpi
      REAL :: ne_te_phase_spectrum_save(nkym,maxmodes)
      REAL :: nsts_phase_spectrum_save(nsm,nkym,maxmodes)
      REAL :: eigenvalue_spectrum_save(2,nkym,maxmodes)
      REAL :: field_spectrum_save(4,nkym,maxmodes)
      REAL :: QL_field_spectrum_save(4,nkym,maxmodes)
      REAL :: intensity_spectrum_save(4,nsm,nkym,maxmodes)
      REAL :: QL_intensity_spectrum_save(4,nsm,nkym,maxmodes)
      REAL :: flux_spectrum_save(5,nsm,3,nkym,maxmodes)
      REAL :: QL_flux_spectrum_save(5,nsm,3,nkym,maxmodes)
      REAL :: spectral_shift_save(nkym)
      REAL :: ave_p0_spectrum_save(nkym)
      REAL :: width_out_save(nkym)
      INTEGER :: ierr
!
!
! setup the ky-spectrum
!
      CALL get_ky_spectrum
!
! initialize mpi buffers for output arrays 
!
      do i=1,nky
       spectral_shift_save(i) = 0.0
       ave_p0_spectrum_save(i) = 0.0
       width_out_save(i)=0.0
       do k=1,nmodes_in
        do t = 1,2
          eigenvalue_spectrum_save(t,i,k) = 0.0
        enddo
        do t=1,4
          field_spectrum_save(t,i,k) = 0.0
          QL_field_spectrum_save(t,i,k) = 0.0
        enddo
        do is=ns0,ns
          do t=1,4
            intensity_spectrum_save(t,is,i,k) = 0.0
            QL_intensity_spectrum_save(t,is,i,k) = 0.0
          enddo
          do j=1,3
            do t=1,5
              flux_spectrum_save(t,is,j,i,k) = 0.0
              QL_flux_spectrum_save(t,is,j,i,k) = 0.0
            enddo
          enddo ! j
          nsts_phase_spectrum_save(is,i,k)=0.0
        enddo ! is
        ne_te_phase_spectrum_save(i,k)=0.0
       enddo  !k
      enddo  !i
!
! loop over ky spectrum
!
! save maximum width
      save_width = width_in
      iflux_in=.TRUE.
      gmax = 0.0
      fmax = 0.0
      do i=1+iProcTglf,nky,nProcTglf
        ky_s = ky_spectrum(i)
!
        new_width=.TRUE.
!
        if(new_eikonal_in)then
          if(jmax_out.eq.0)then   ! first pass
            if(find_width_in)then
              CALL tglf_max
            else
              nbasis = nbasis_max_in
              new_width = .TRUE.
              CALL tglf_LS
              gamma_nb_min_out = gamma_out(1)
            endif
          else   ! second pass
            gamma_reference_kx0(:) = eigenvalue_first_pass(1,i,:)
            freq_reference_kx0(:) = eigenvalue_first_pass(2,i,:)
            width_in = width_out(i)
            nbasis = nbasis_max_in
            new_width=.TRUE.
            CALL tglf_LS
            gamma_nb_min_out = gamma_out(1)
          endif
          mask_save(i) = 1
          if(gamma_out(1).eq.0.0)mask_save(i)=0
!          write(*,*)i,"ky=",ky_s,mask_save(i),gamma_out(1)
          gamma_nb_min_save(i) = gamma_nb_min_out
          width_save(i) = width_in
          ft_save(i) = ft
          R_unit_save(i) = R_unit
          q_unit_save(i) = q_unit
          do j=1,nx
            wdx_save(i,j) = wdx(j)
            b0x_save(i,j) = b0x(j)
            cx_par_par_save(i,j) = cx_par_par(j)
            cx_tor_par_save(i,j) = cx_tor_par(j)
            cx_tor_per_save(i,j) = cx_tor_per(j)
            kxx_save(i,j) = kxx(j)
          enddo
        else
          gamma_nb_min_out = gamma_nb_min_save(i)
          width_in = width_save(i)
          ft = ft_save(i)
          R_unit = R_unit_save(i)
          q_unit = q_unit_save(i)
          do j=1,nx
             wdx(j) = wdx_save(i,j)
             b0x(j) = b0x_save(i,j)
             cx_par_par(j) = cx_par_par_save(i,j)
             cx_tor_par(j) = cx_tor_par_save(i,j)
             cx_tor_per(j) = cx_tor_per_save(i,j)
             kxx(j) = kxx_save(i,j)
          enddo
          if(mask_save(i).eq.1)then
            CALL tglf_LS
          else
            gamma_out(1)=0.0
          endif
        endif
!        write(*,*)i,"ky=",ky_s,"width=",width_in
!        write(*,*)"nbasis=",nbasis_max_in,nbasis_min_in
!        write(*,*)"ft=",ft,"R=",R_unit,"q=",q_unit
!        write(*,*)"wdx=",wdx(1),"b0x=",b0x(1)
!
        unstable=.TRUE.
        gamma_max = MAX(gamma_out(1),gamma_out(2)) ! this covers ibranch=-1,0
        if(gamma_max.eq.0.0.or.gamma_nb_min_out.eq.0.0)unstable=.FALSE.      
        gamma_net_1 = gamma_nb_min_out 
!       gamma_cutoff = (0.1*ky_s/R_unit)*SQRT(taus(1)*mass(2))  ! scaled like gamma
        gamma_cutoff = 0.1*ky_s/R_unit
        rexp = 1.0
        reduce = 1.0
        if(nbasis_max_in.ne.nbasis_min_in)then
          if(gamma_net_1.lt.gamma_out(1) &
             .and.gamma_net_1.lt.gamma_cutoff)then
          reduce = (gamma_net_1/gamma_cutoff)**rexp
!            write(*,*)"phi reduced",ky_s,gamma_nb_min_out,gamma_out(1)
          endif
       endif
       if(sat_rule_in.ge.1)reduce=1.0
!
       width_out_save(i) = width_in
!
       if(unstable)then
! save the spectral shift of the radial wavenumber due to VEXB_SHEAR
         spectral_shift_save(i) = kx0_e
         ave_p0_spectrum_save(i) = ave_p0_out
! save field_spectrum_out and eigenvalue_spectrum_out
         do imax=1,nmodes_out
           QL_field_spectrum_save(1,i,imax) = v_QL_out(imax)
           QL_field_spectrum_save(2,i,imax) = 1.0
           QL_field_spectrum_save(3,i,imax) = a_par_QL_out(imax)
           QL_field_spectrum_save(4,i,imax) = b_par_QL_out(imax)
           field_spectrum_save(1,i,imax) = reduce*v_bar_out(imax)
           field_spectrum_save(2,i,imax) = reduce*phi_bar_out(imax)
           field_spectrum_save(3,i,imax) = reduce*a_par_bar_out(imax)
           field_spectrum_save(4,i,imax) = reduce*b_par_bar_out(imax)
           eigenvalue_spectrum_save(1,i,imax)=gamma_out(imax)
           eigenvalue_spectrum_save(2,i,imax)=freq_out(imax)
           if(ky_s.le.1.0.and.gamma_out(imax).gt.gmax)then
             gmax=gamma_out(imax)
             fmax=freq_out(imax)
           endif
!          write(*,*)ky_s,width_in
!          write(*,*)"modes",imax,phi_QL_out(imax)
!          write(*,*)gamma_out(imax),freq_out(imax)
         enddo
! save intensity_spectrum_out
         do is=ns0,ns
          do imax=1,nmodes_out
            QL_intensity_spectrum_save(1,is,i,imax) = N_QL_out(imax,is)
            QL_intensity_spectrum_save(2,is,i,imax) = T_QL_out(imax,is)
            QL_intensity_spectrum_save(3,is,i,imax) = U_QL_out(imax,is)
            QL_intensity_spectrum_save(4,is,i,imax) = Q_QL_out(imax,is)
            intensity_spectrum_save(1,is,i,imax) = N_bar_out(imax,is)
            intensity_spectrum_save(2,is,i,imax) = T_bar_out(imax,is)
            intensity_spectrum_save(3,is,i,imax) = U_bar_out(imax,is)
            intensity_spectrum_save(4,is,i,imax) = Q_bar_out(imax,is)
           enddo !imax
         enddo  ! is
! save flux_spectrum_out 
         do is=ns0,ns
          do j=1,3
            do imax=1,nmodes_out
              phi_bar = reduce*phi_bar_out(imax)
              pflux1 = particle_QL_out(imax,is,j)
              eflux1 = energy_QL_out(imax,is,j)
              stress_tor1 = stress_tor_QL_out(imax,is,j)
              stress_par1 = stress_par_QL_out(imax,is,j)
              exch1 = exchange_QL_out(imax,is,j)
              QL_flux_spectrum_save(1,is,j,i,imax) = pflux1
              QL_flux_spectrum_save(2,is,j,i,imax) = eflux1
              QL_flux_spectrum_save(3,is,j,i,imax) = stress_tor1
              QL_flux_spectrum_save(4,is,j,i,imax) = stress_par1
              QL_flux_spectrum_save(5,is,j,i,imax) = exch1
              flux_spectrum_save(1,is,j,i,imax) = phi_bar*pflux1
              flux_spectrum_save(2,is,j,i,imax) = phi_bar*eflux1
              flux_spectrum_save(3,is,j,i,imax) = phi_bar*stress_tor1
              flux_spectrum_save(4,is,j,i,imax) = phi_bar*stress_par1
              flux_spectrum_save(5,is,j,i,imax) = phi_bar*exch1
            enddo !imax
           enddo ! j
         enddo  ! is 
! save ne_te crossphase
         do imax=1,nmodes_out
           ne_te_phase_spectrum_save(i,imax) = ne_te_phase_out(imax)
         enddo  !imax
! save ns_ts crossphase
         do is=1,ns_in
           do imax=1,nmodes_out
             nsts_phase_spectrum_save(is,i,imax) = Ns_Ts_phase_out(imax,is)
           enddo  !imax
          enddo    ! is

        endif !unstable .T.
!
! reset width to maximum if used tglf_max
        if(find_width_in)width_in=save_width
!
      enddo  ! i 
!
! collect and broadcast the results 
      call MPI_BARRIER(iCommTglf,ierr)

      call MPI_ALLREDUCE(ne_te_phase_spectrum_save   &
                        ,ne_te_phase_spectrum_out    &
                        ,nkym*maxmodes               &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
       call MPI_ALLREDUCE(nsts_phase_spectrum_save   &
                        ,nsts_phase_spectrum_out     &
                        ,nsm*nkym*maxmodes           &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
      call MPI_ALLREDUCE(eigenvalue_spectrum_save    &
                        ,eigenvalue_spectrum_out     &
                        ,2*nkym*maxmodes             &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
      call MPI_ALLREDUCE(field_spectrum_save         &
                        ,field_spectrum_out          &
                        ,4*nkym*maxmodes             &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
      call MPI_ALLREDUCE(QL_field_spectrum_save      &
                        ,QL_field_spectrum_out       &
                        ,4*nkym*maxmodes             &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
      call MPI_ALLREDUCE(intensity_spectrum_save     &
                        ,intensity_spectrum_out      &
                        ,4*nsm*nkym*maxmodes         &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
      call MPI_ALLREDUCE(QL_intensity_spectrum_save  &
                        ,QL_intensity_spectrum_out   &
                        ,4*nsm*nkym*maxmodes         &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
      call MPI_ALLREDUCE(flux_spectrum_save          &
                        ,flux_spectrum_out           &
                        ,5*nsm*3*nkym*maxmodes       &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
      call MPI_ALLREDUCE(QL_flux_spectrum_save       &
                        ,QL_flux_spectrum_out        &
                        ,5*nsm*3*nkym*maxmodes       &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
      call MPI_ALLREDUCE(spectral_shift_save         &
                        ,spectral_shift_out          &
                        ,nkym                        &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
      call MPI_ALLREDUCE(ave_p0_spectrum_save         &
                        ,ave_p0_spectrum_out          &
                        ,nkym                        &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)
     call MPI_ALLREDUCE(width_out_save               &
                        ,width_out                   &
                        ,nkym                        &
                        ,MPI_DOUBLE_PRECISION        &
                        ,MPI_SUM                     &
                        ,iCommTglf                   &
                        ,ierr)

      ! recompute spectrum using non-local in ky multiscale model
      
      if(sat_rule_in.ge.1)call get_multiscale_spectrum
      
      !
      if(new_eikonal_in)eikonal_unsaved=.FALSE.
      gamma_out(1) = gmax
      freq_out(1) = fmax
      new_eikonal_in = .TRUE.  ! reset default for next call to tglf_TM
      
      END SUBROUTINE get_bilinear_spectrum_mpi
