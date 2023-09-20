!
!--------------------------------------------------------------
!
      SUBROUTINE tglf_LS
!
!***********************************************************************
! questions  should be addressed to
!  Gary Staebler 858-455-3466  or email: gary.staebler@gat.com
!***********************************************************************
!  TGLF Linear Stability driver computes all quasilinear quantities 
!  for a single ky
!
!***********************************************************************
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_eigen
      USE tglf_weight
!
      IMPLICIT NONE
      REAL,PARAMETER :: epsilon1 = 1.E-12
      LOGICAL :: NO_TRAPPED = .false.
      INTEGER :: j1, j, i, jmax(maxmodes), iroot
      INTEGER :: imax,is
      INTEGER :: mi,me
      REAL :: zgamax
      REAL :: phi2_bar,kyi
      REAL :: sum_v_bar, sum_modB_bar
      REAL :: get_intensity, get_gamma_net
!  ZGESV storage
      REAL :: small = 1.0E-13
      INTEGER :: info
      INTEGER,ALLOCATABLE,DIMENSION(:) :: di,de
      INTEGER,ALLOCATABLE,DIMENSION(:) :: ipiv
      COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: zmat
! 
!      cputime0=MPI_WTIME()
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
!      write(*,*)"xnue_in=",xnue_in
!      write(*,*)"zeff_in=",zeff_in
!      write(*,*)"filter_in=",filter_in
!      write(*,*)"park_in=",park_in
!      write(*,*)"ghat_in=",ghat_in
!      write(*,*)"wd_zero_in=",wd_zero_in
!      write(*,*)"Linsker_factor_in=",Linsker_factor_in
!      write(*,*)"gradB_factor_in=",gradB_factor_in
!      write(*,*)"xnu_factor_in=",xnu_factor_in
!      write(*,*)"debye_factor_in=",debye_factor_in
!      write(*,*)"theta_trapped_in=",theta_trapped_in
!      write(*,*)"alpha_p_in=",alpha_p_in
!      write(*,*)"alpha_e_in=",alpha_e_in
!
      if(new_eikonal_in)then
        if(new_geometry)then
          if(igeo.eq.1)then
! set up MILLER geometry
            trace_path(4)=1
            call miller_geo
          elseif(igeo.eq.2)then
! set up FOURIER geometry
            trace_path(5)=1
            call fourier_geo
          elseif(igeo.eq.3)then
! set up ELITE geometry
            trace_path(8)=1
            call ELITE_geo
          endif
!        write(*,*)"geometry done"
!
! compute the eikonal functions for general geometry (igeo>0)
!
          if(igeo.gt.0)call mercier_luc 
!
        endif ! new_geometry
        new_geometry=.FALSE.
!
!  load the x-grid eikonal functions wdx,b0x
!
        if(new_width)then
          trace_path(6)=1
          call get_xgrid_functions
        endif
!       write(*,*)"ft=",ft
      endif  !new_eikonal_in
!      write(*,*)"eikonal done"
!
      !
      !  compute the rank of the eigenmatrix iur
      !
      nroot=15
      if(NO_TRAPPED)nroot=6
      iur  = (ns-ns0+1)*nroot*nbasis
      !      write(*,*)"iur = ",iur,"nroot=",nroot,"ns0=",ns0,"ft=",ft
!
if(new_matrix)then
        trace_path(7)=1
        call get_matrix
      endif
!      write(*,*)"matrix done"
!
! allocate eigenvalues, eigenvectors
!
!   tglf_LS local 
      ALLOCATE(di(iur))
      ALLOCATE(de(iur))
      ALLOCATE(ipiv(iur))
      ALLOCATE(zmat(iur,iur))
!   tglf_eigen module
      ALLOCATE(fv1(iur))
      ALLOCATE(fv2(iur))
      ALLOCATE(fv3(iur))            
      ALLOCATE(rr(iur))
      ALLOCATE(ri(iur))
      ALLOCATE(ar(iur,iur))
      ALLOCATE(ai(iur,iur))          
      ALLOCATE(vr(iur,iur))
      ALLOCATE(vi(iur,iur))
      ALLOCATE(amat(iur,iur))
      ALLOCATE(bmat(iur,iur))
      ALLOCATE(v(iur))
      ALLOCATE(alpha(iur))
      ALLOCATE(beta(iur))
!      write(*,*)"eigen allocation done"
!
!  solver for linear eigenmodes of tglf equations
!
      call tglf_eigensolver
!      write(*,*)"eigensolver done"
!
!      initalize output to zero
!
      do j1=1,maxmodes
        jmax(j1) = 0
        gamma_out(j1) = 0.0
        freq_out(j1) = 0.0
        v_QL_out(j1) = 0.0
        a_par_QL_out(j1) = 0.0
        b_par_QL_out(j1) = 0.0
        phi_bar_out(j1) = 0.0
        a_par_bar_out(j1) = 0.0
        b_par_bar_out(j1) = 0.0
        v_bar_out(j1) = 0.0
        do is=1,nsm
          do j=1,3
            particle_QL_out(j1,is,j) = 0.0
            energy_QL_out(j1,is,j) = 0.0
            stress_par_QL_out(j1,is,j) = 0.0
            stress_tor_QL_out(j1,is,j) = 0.0
            exchange_QL_out(j1,is,j) = 0.0
          enddo
          N_QL_out(j1,is) = 0.0
          T_QL_out(j1,is) = 0.0
          U_QL_out(j1,is) = 0.0
          Q_QL_out(j1,is) = 0.0
          N_bar_out(j1,is) = 0.0
          T_bar_out(j1,is) = 0.0
          U_bar_out(j1,is) = 0.0
          Q_bar_out(j1,is) = 0.0
          Ns_Ts_phase_out(j1,is) = 0.0
        enddo
        ne_te_phase_out(j1) = 0.0
      enddo
!
      if(ibranch_in.eq.0)then
!       write(*,*)"nmodes_in=",nmodes_in
      nmodes_out = nmodes_in
!
! sort the unstable modes into electron and ion frequencies
      mi = 0
      me = 0
      do j1=1,iur
        if(rr(j1).gt.epsilon1)then
         if(ri(j1).gt.0.0)then
! note that ri = -freq, rr = gamma
           mi = mi+1
           di(mi)=j1
         else
           me = me+1
           de(me)=j1
         endif
!          write(*,*)"debug sort",j1,rr(j1),ri(j1),me,mi
        endif
      enddo
! find the most unstable mode for each branch
      if(me.gt.0)then
        zgamax = 0.0
        do iroot=1,me
         if(rr(de(iroot)).gt.zgamax)then
           zgamax = rr(de(iroot))
           jmax(1) = de(iroot)
         endif
        enddo
        gamma_out(1)=rr(jmax(1))
        freq_out(1)=-ri(jmax(1))
      endif
      if(mi.gt.0)then
        zgamax = 0.0
        do iroot=1,mi
          if(rr(di(iroot)).gt.zgamax)then
           zgamax = rr(di(iroot))
           jmax(2) = di(iroot)
         endif
        enddo
        gamma_out(2)=rr(jmax(2))
        freq_out(2)=-ri(jmax(2))
      endif
!        write(*,*)"debug sort", ky
!        write(*,*)"nmodes_out = ",nmodes_out,nmodes_in
!        write(*,*)"me = ",me,"mi = ",mi
!        write(*,*)"gamma_out=",gamma_out(1),gamma_out(2)
!        write(*,*)"freq_out=",freq_out(1),freq_out(2)
      endif ! ibranch_in .eq.0 
      if(ibranch_in.eq.-1)then
!
!  find the top nmodes most unstable modes
       CALL sort_eigenvalues(nmodes_in,jmax)
        nmodes_out = 0
         do j1=1,nmodes_in 
          if(jmax(j1).ne.0)then
            nmodes_out = nmodes_out + 1
            gamma_out(j1)=rr(jmax(j1))
            freq_out(j1)=-ri(jmax(j1))
          endif
        enddo
!       write(*,*)"gamma_out = ",gamma_out(1)
!      write(*,*)"debug jmax =",jmax(1),jmax(2)
!      write(*,*)"debug nmodes_out = ",nmodes_out
      endif ! ibranch_in .eq. -1
!
      if(alpha_quench_in.ne.0.0)then
! apply quench rule
        do j1=1,nmodes_in
          gamma_out(j1) = get_gamma_net(gamma_out(j1))
        enddo
      elseif(vexb_shear_s.ne.0.0)then
! use spectral shift model for second pass
       do j1=1,nmodes_in
          gamma_out(j1) = gamma_reference_kx0(j1)
          freq_out(j1) = freq_reference_kx0(j1)
        enddo
      endif
!      
!  get the fluxes for the most unstable modes
      if(iflux_in)then
        do imax=1,nmodes_out
         if(jmax(imax).gt.0)then
          do i=1,iur
            v(i) = small
            do j=1,iur
              zmat(i,j) = beta(jmax(imax))*amat(i,j)- &
                         (small +alpha(jmax(imax)))*bmat(i,j)
            enddo
          enddo
          call zgesv(iur,1,zmat,iur,ipiv,v,iur,info)

          ! Check for clean exit from ZGESV

          if (info /= 0)CALL tglf_error(1,"ZGESV failed in tglf_LS")

!  alpha/beta=-xi*(frequency+xi*growthrate)
          eigenvalue = xi*alpha(jmax(imax))/beta(jmax(imax))  
!          write(*,*)"eigenvalue=",eigenvalue,imax
          call get_QL_weights
!
          wd_bar_out(imax)=wd_bar
          b0_bar_out(imax)=b0_bar
          modB_bar_out(imax)=modB_bar
          v_QL_out(imax)=v_weight
          a_par_QL_out(imax)=a_par_weight
          b_par_QL_out(imax)=b_par_weight
          kx_bar_out(imax)=kx_bar
          kpar_bar_out(imax)=kpar_bar/(R_unit*q_unit*width_in)
          do i=1,nbasis
            do j=1,3
             field_weight_out(imax,j,i)=field_weight_QL_out(j,i)
            enddo
          enddo
!
          do is=ns0,ns
            do j=1,3
              particle_QL_out(imax,is,j)=particle_weight(is,j)
              energy_QL_out(imax,is,j)=energy_weight(is,j)
              stress_par_QL_out(imax,is,j)=stress_par_weight(is,j)
              stress_tor_QL_out(imax,is,j)=stress_tor_weight(is,j)
              exchange_QL_out(imax,is,j)=exchange_weight(is,j)
            enddo 
            N_QL_out(imax,is)=N_weight(is)
            T_QL_out(imax,is)=T_weight(is)
            U_QL_out(imax,is)=U_weight(is)
            Q_QL_out(imax,is)=Q_weight(is)
            Ns_Ts_phase_out(imax,is)=Ns_Ts_phase(is)
          enddo
          ne_te_phase_out(imax) = Ne_Te_phase
          ! SAT_RULE = 0 intensity model assumes the vector norm v_bar
          kyi=ky
          if(ABS(v_QL_out(imax)).lt.epsilon1)then
            v_bar_out(imax)=0.0
            phi2_bar = 0.0
          else
            v_bar_out(imax) = get_intensity(kyi,gamma_out(imax))
            phi2_bar = v_bar_out(imax)/v_QL_out(imax)
          endif
!
          phi_bar_out(imax) = phi2_bar
          a_par_bar_out(imax) = phi2_bar*a_par_QL_out(imax)
          b_par_bar_out(imax) = phi2_bar*b_par_QL_out(imax)
          do is=ns0,ns
            N_bar_out(imax,is)=phi2_bar*N_QL_out(imax,is)
            T_bar_out(imax,is)=phi2_bar*T_QL_out(imax,is)
            U_bar_out(imax,is)=phi2_bar*U_QL_out(imax,is)
            Q_bar_out(imax,is)=phi2_bar*Q_QL_out(imax,is)
          enddo
         endif
        enddo !imax
! check for inward ballooing 
     ft_test = 0.0
     sum_modB_bar=0.0
     sum_v_bar = 0.0
     do i=1,nmodes_out
       sum_modB_bar = sum_modB_bar + v_bar_out(i)*modB_bar_out(i)
       sum_v_bar = sum_v_bar + v_bar_out(i)
     enddo
     if(sum_v_bar.gt.epsilon1)ft_test = sum_modB_bar/sum_v_bar
     ft_test = ft_test/modB_min
!     write(*,*)ky,"ft_test=",ft_test
!           
      endif  ! iflux_in
!
!
!
! dealocate eigenvalues, eigenvectors 
!
!    write(*,*)"ready to deallocate tglf_LS"
!  tglf_LS local
      if(ALLOCATED(di))DEALLOCATE(di)
      if(ALLOCATED(de))DEALLOCATE(de)
      if(ALLOCATED(ipiv))DEALLOCATE(ipiv)
      if(ALLOCATED(zmat))DEALLOCATE(zmat)
!      write(*,*)"deallocated tglf_LS local"
!  tglf_eigen module
      if(ALLOCATED(fv1))DEALLOCATE(fv1)
      if(ALLOCATED(fv2))DEALLOCATE(fv2)
      if(ALLOCATED(fv3))DEALLOCATE(fv3)
      if(ALLOCATED(rr))DEALLOCATE(rr)
      if(ALLOCATED(ri))DEALLOCATE(ri)
      if(ALLOCATED(ar))DEALLOCATE(ar)
      if(ALLOCATED(ai))DEALLOCATE(ai)
      if(ALLOCATED(vr))DEALLOCATE(vr)
      if(ALLOCATED(vi))DEALLOCATE(vi)
      if(ALLOCATED(amat))DEALLOCATE(amat)
      if(ALLOCATED(bmat))DEALLOCATE(bmat)
      if(ALLOCATED(v))DEALLOCATE(v)
      if(ALLOCATED(alpha))DEALLOCATE(alpha)
      if(ALLOCATED(beta))DEALLOCATE(beta)
!      write(*,*)"deallocated eigen"
!

!
      END SUBROUTINE tglf_LS
!
!
!--------------------------------------------------------------
!
      REAL FUNCTION get_intensity(kp,gp)
!
      USE tglf_species
      USE tglf_global
      USE tglf_coeff
      IMPLICIT NONE
!
      REAL,INTENT(IN) :: kp,gp
      REAL :: cnorm,exponent1
      REAL :: wd0,gnet
      REAL :: c1,pols,ks
      REAL :: measure
      REAL :: get_GAM_freq
      REAL :: intensity
!
      pols = (ave_p0(1,1)/ABS(as(1)*zs(1)*zs(1)))**2 ! scale invariant pol
      ks = kp*SQRT(taus(1)*mass(2))/ABS(zs(1))   ! scale invariant gyroradius * poloidal wavenumber
      measure = SQRT(taus(1)*mass(2))
      if(sat_rule_in.eq.0)then
       if(igeo.eq.0)then
        if(nmodes_in.le.2)then
! this fit is for nmodes_in=2
          cnorm = 30.40*pols
          exponent1 = 1.657
        else
! this fit is for nmodes_in=4
          cnorm = 24.58*pols
          exponent1 = 1.761
        endif
        if(ks.gt.1.0)cnorm=cnorm/(ks)**etg_factor_in
        c1 = 0.0
      elseif(igeo.ge.1)then
        if(nmodes_in.le.2)then
! this fit is for nmodes_in=2
          cnorm = 32.48*pols
          exponent1 = 1.547
          c1 = 0.534
        else
! this fit is for nmodes_in=4
          cnorm = 30.03*pols
          exponent1 = 1.66
          c1 = 0.234
        endif
        if(ks.gt.1.0)cnorm=cnorm/(ks)**etg_factor_in
       endif
       wd0 =ks*SQRT(taus(1)/mass(2))/R_unit  ! renomalized for scale invariance
       gnet = gp/wd0
       intensity = cnorm*(wd0**2)*(gnet**exponent1 &
        + c1*gnet)/(kp**4)
       if(alpha_quench_in.eq.0.0.and.ABS(kx0_e).gt.0.0)then
         intensity = intensity/(1.0+0.56*kx0_e**2)**2
         intensity = intensity/(1.0+(1.15*kx0_e)**4)**2
       endif
         intensity = intensity*SAT_geo0_out*measure
!       write(*,*)" sat0 ",SAT_geo0_out
      elseif(sat_rule_in.ge.1)then
!
!   will be computed later by get_multiscale_spectrum
!
       intensity = 1.0
      endif
!     
      get_intensity = intensity
!
      END FUNCTION get_intensity
!
!--------------------------------------------------------------
!
      REAL FUNCTION get_gamma_net(gp)
!
      USE tglf_global
      USE tglf_species
!
      REAL,INTENT(IN) :: gp
      REAL :: alpha_exb
!
      alpha_exb = 0.3
      if(igeo.eq.1)alpha_exb=0.3*SQRT(kappa_loc)
      get_gamma_net =  MAX(gp - ABS(alpha_exb*alpha_quench_in*vexb_shear_s),0.0)
!
      END FUNCTION get_gamma_net
!
!--------------------------------------------------------------
!
      REAL FUNCTION get_GAM_freq()
!
      USE tglf_global
      USE tglf_species
!
      REAL :: elong=1.0
!
      if(igeo.eq.1)elong = kappa_loc
!      get_GAM_freq = SQRT(3.5*taus(2)+2.0*taus(1))/R_unit
      get_GAM_freq = (2.0/(1.0+elong))*SQRT(taus(2)+taus(1))/R_unit
!
      END FUNCTION get_GAM_freq
!
!--------------------------------------------------------------
!
!
      SUBROUTINE sort_eigenvalues(nsorted,sorted_index)
!
      USE tglf_dimensions
      USE tglf_eigen
!
      IMPLICIT NONE
      INTEGER :: i, j, ngrow,imid, iend, it, nsorted
      INTEGER :: grow_index(iur)
      INTEGER :: sorted_index(maxmodes)
      REAL,PARAMETER ::  epsilon1=1.0E-12
!
!  find all of the unstable modes and the most unstable
!
      ngrow = 0
!      write(*,*)iur,"eigenvalues"
      do i=1,iur
!        growthrate(i)=0.0
!        frequency(i)=0.0
        grow_index(i)=0
        if(rr(i).gt.epsilon1)then
          ngrow = ngrow + 1
          grow_index(ngrow) = i
!          write(*,*)
!     >    ngrow,grow_index(ngrow),rr(grow_index(ngrow))
        endif 
      enddo
!
      if(ngrow.le.1)go to 30
!
! put the unstable modes in ascending order by growthrate
! using a heapsort algorithm 
!
      imid = ngrow/2 + 1
      iend = ngrow
!      write(*,*)"imid = ",imid,"iend = ",iend
 10   continue
        if(imid.gt.1)then
          imid = imid-1
          it = grow_index(imid)
        else
          it = grow_index(iend)
          grow_index(iend) = grow_index(1)
          iend = iend-1
          if(iend.eq.1)then
            grow_index(1)=it
          do i=1,ngrow
!           write(*,*)i,grow_index(i),rr(grow_index(i))
          enddo
            go to 30
         endif
        endif
        i=imid
        j = imid+imid
!        write(*,*)"i = ",i,"j = ",j
 20     if(j.le.iend)then
          if(j.lt.iend)then
           if(rr(grow_index(j)).gt.rr(grow_index(j+1)))j=j+1
          endif
          if(rr(it).gt.rr(grow_index(j)))then
            grow_index(i)=grow_index(j)
            i=j
            j=j+j
          else
            j = iend+1
          endif
        go to 20
        endif
        grow_index(i)=it
      go to 10
 30   continue
      do i=1,nsorted
        sorted_index(i)=0
        if(i.le.ngrow)sorted_index(i)=grow_index(i)
      enddo
! debug
!      do i=1,ngrow
!        growthrate(i)=rr(grow_index(i))
!        frequency(i)=-ri(grow_index(i))
!      enddo
!      do i=1,ngrow
!      write(*,*)i,growthrate(i),frequency(i),grow_index(i)
!      enddo
!
      END SUBROUTINE sort_eigenvalues
!
      SUBROUTINE get_QL_weights
! **************************************************************
!
! compute the quasilinear weights for a single eigenmode
! with eigenvector v. All of the QL weights are normalized to phi_norm
!
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
!      USE tglf_hermite
      USE tglf_eigen
      USE tglf_coeff
      USE tglf_weight
      USE tglf_sgrid
! 
      IMPLICIT NONE
!
      INTEGER :: i,j,is
      COMPLEX :: n(nsm,nb)
      COMPLEX :: u_par(nsm,nb)
      COMPLEX :: p_par(nsm,nb)
      COMPLEX :: p_tot(nsm,nb)
      COMPLEX :: q_par(nsm,nb)
      COMPLEX :: q_tot(nsm,nb)
!      COMPLEX :: ng(nsm,nb)
!      COMPLEX :: ug_par(nsm,nb)
!      COMPLEX :: pg_par(nsm,nb)
!      COMPLEX :: pg_tot(nsm,nb)
!      COMPLEX :: qg_par(nsm,nb)
!      COMPLEX :: qg_tot(nsm,nb)
!      COMPLEX :: nt(nsm,nb)
!      COMPLEX :: pt_par(nsm,nb)
!      COMPLEX :: pt_tot(nsm,nb)
      COMPLEX :: temp(nsm,nb)
      COMPLEX :: stress_par(nsm,nb,3),stress_per(nsm,nb,3)
      COMPLEX,DIMENSION(nb) :: phi=0.0,psi=0.0,bsig=0.0
      COMPLEX :: phi_wd_phi,wd_phi
      COMPLEX :: phi_b0_phi,b0_phi
      COMPLEX :: phi_modB_phi,modB_phi
      COMPLEX :: phi_kx_phi,kx_phi
      COMPLEX :: phi_kpar_phi,kpar_phi
      COMPLEX :: freq_QL
      REAL :: betae_psi,betae_sig
      REAL :: phi_norm,psi_norm,bsig_norm,vnorm
      REAL :: epsilon1
      REAL :: stress_correction,wp
!
!      xi=(0.0,1.0)
      epsilon1 = 1.E-12
      freq_QL = eigenvalue
!      write(*,*)"freq_QL = ",freq_QL
!
!  fill the density and total pressure vectors
!
      vnorm = 0.0
      do is = ns0,ns
        j = (is-ns0)*nroot*nbasis
        do i=1,nbasis          
          n(is,i) = v(j+i)
          u_par(is,i) = v(j+nbasis+i)
          p_par(is,i) = v(j+nbasis*2+i)
          p_tot(is,i) = v(j+nbasis*3+i)
          q_par(is,i) = v(j+nbasis*4+1)
          q_tot(is,i) = v(j+nbasis*5+i)
          if(nroot.gt.6)then
            n(is,i) = n(is,i) -v(j+nbasis*6+i)+v(j+nbasis*12+i)
            u_par(is,i) = u_par(is,i) -v(j+nbasis*7+i)
            p_par(is,i) = p_par(is,i) -v(j+nbasis*8+i)+v(j+nbasis*13+i)
            p_tot(is,i) = p_tot(is,i) -v(j+nbasis*9+i)+v(j+nbasis*14+i)
            q_par(is,i) = q_par(is,i) -v(j+nbasis*10+i)
            q_tot(is,i) = q_tot(is,i) -v(j+nbasis*11+i)
          endif
        enddo
      enddo
! 
!  compute vnorm
!
      vnorm = 0.0
      if(ns.le.2)then
        do i=1,iur
          vnorm = vnorm + REAL(v(i)*CONJG(v(i)))
        enddo
      else ! weight vnorm by equililibrium densities
        j=1
        do i=1,iur
          if(i.gt.j*nbasis*nroot)j=j+1
          vnorm = vnorm + REAL(v(i)*CONJG(v(i)))*ABS(as(j)*zs(j))
        enddo
        vnorm = vnorm/ABS(as(1)*zs(1))   ! normalize to electron charge density
      endif
!      write(*,*)"vnorm =",vnorm
!
!  compute the electromagnetic potentials
!
      betae_psi = 0.0
      if(use_bper_in)betae_psi = 0.5*betae_s/(ky*ky)
      betae_sig = 0.0
      if(use_bpar_in)betae_sig = 0.5*betae_s
      do i=1,nbasis
        phi(i)=0.0
        psi(i)=0.0
        bsig(i)=0.0
        do is=ns0,ns
          do j=1,nbasis
            phi(i) = phi(i) +ave_p0inv(i,j)*as(is)*zs(is)*n(is,j)  
          enddo
          if(use_bper_in)then
            do j=1,nbasis
              psi(i) = psi(i) + &
              betae_psi*ave_b0inv(i,j)*as(is)*zs(is)*vs(is)*u_par(is,j)
            enddo
            if(vpar_model_in.eq.0)then
              do j=1,nbasis
               phi(i) = phi(i) + U0*betae_psi*ave_bpinv(i,j)*as(is)*zs(is)*vs(is)*u_par(is,j)
               psi(i) = psi(i) - U0*betae_psi*ave_bpinv(i,j)*as(is)*zs(is)*n(is,j)
              enddo
            endif
          endif
          if(use_bpar_in)then
              bsig(i) = bsig(i) - betae_sig*as(is)*taus(is)* &
                   (1.5*p_tot(is,i)-0.5*p_par(is,i))
          endif
        enddo
      enddo
! 
!  add the adiabatic terms to the total moments
!    
      do is=ns0,ns
        do i=1,nbasis
          n(is,i) = n(is,i) - zs(is)*phi(i)/taus(is)
          p_par(is,i) = p_par(is,i) - zs(is)*phi(i)/taus(is)
          p_tot(is,i) = p_tot(is,i) - zs(is)*phi(i)/taus(is)
        enddo
      enddo
!
!  compute phi_norm, psi_norm, sig_norm
!
      phi_norm = 0.0
      psi_norm = 0.0
      bsig_norm = 0.0
      do i=1,nbasis
        phi_norm = phi_norm + REAL(phi(i)*CONJG(phi(i)))
        psi_norm = psi_norm + REAL(psi(i)*CONJG(psi(i)))
        bsig_norm = bsig_norm + REAL(bsig(i)*CONJG(bsig(i)))
      enddo
      if(phi_norm.lt.epsilon1)phi_norm = epsilon1
!      write(*,*)"phi_norm =",phi_norm
!
! save the field weights
      do i=1,nbasis
        field_weight_QL_out(1,i) = xi*phi(i)/SQRT(phi_norm)
        field_weight_QL_out(2,i) = xi*psi(i)/SQRT(phi_norm)
        field_weight_QL_out(3,i) = xi*bsig(i)/SQRT(phi_norm)
!        write(*,*)i,"field_weight=",field_weight_QL_out(1,i),field_weight_QL_out(2,i),field_weight_QL_out(3,i)
       enddo
!
! compute <phi|*|phi> averages
      phi_wd_phi = 0.0
      phi_b0_phi = 0.0
      phi_modB_phi = 0.0
      phi_kx_phi = 0.0
      phi_kpar_phi = 0.0
      do i=1,nbasis
         wd_phi = 0.0
         b0_phi = 0.0
         modB_phi = 0.0
         kx_phi = 0.0
         kpar_phi = 0.0
         do j=1,nbasis
           wd_phi = wd_phi +ave_wdh(i,j)*phi(j)
           b0_phi = b0_phi +ave_b0(i,j)*phi(j)
           modB_phi = modB_phi +ave_c_par_par(i,j)*phi(j)
           kx_phi = kx_phi +ave_kx(i,j)*phi(j)
           kpar_phi = kpar_phi +xi*ave_kpar(i,j)*phi(j)
         enddo
         phi_wd_phi = phi_wd_phi + CONJG(phi(i))*wd_phi
         phi_b0_phi = phi_b0_phi + CONJG(phi(i))*b0_phi
         phi_modB_phi = phi_modB_phi + CONJG(phi(i))*modB_phi
         phi_kx_phi = phi_kx_phi + CONJG(phi(i))*kx_phi
         phi_kpar_phi = phi_kpar_phi + CONJG(phi(i))*kpar_phi
      enddo
      wd_bar = REAL(phi_wd_phi)/phi_norm
      b0_bar = REAL(phi_b0_phi)/phi_norm
      modB_bar = ABS(REAL(phi_modB_phi)/phi_norm)
      kx_bar = REAL(phi_kx_phi)/phi_norm
      kpar_bar = REAL(phi_kpar_phi)/phi_norm
!      write(*,*)"wd_bar = ",wd_bar
!      write(*,*)"b0_bar = ",b0_bar
!       write(*,*)"modB_bar = ",modB_bar
!      write(*,*)"kx_bar = ",kx_bar
!      write(*,*)"kpar_bar = ",kpar_bar
!
! fill the stress moments
!
      stress_correction = 1.0
!
      do is=ns0,ns
        wp = ky*ave_hp1(is,1,1)*ABS(alpha_p_in*vpar_shear_in(is))/vs(is)
!        wp = wp*taus(is)/(ABS(zs(is))*ave_hn(is,1,1))
!        wp = ABS(vpar_shear_in(is))*R_unit/vs(is)
!        wp =(wp*sqrt_two*R_unit*q_unit*width_in/vs(is))/3.6
!        stress_correction = (1.0 + alpha_p_in*wp)/(1.0 + wp)
!        stress_correction = (1.0 + 2.3*wp)/(1.0 + wp)
        if(sat_rule_in.eq.0)stress_correction = (AIMAG(freq_QL)+2.0*wp)/(AIMAG(freq_QL)+wp)
!            write(*,*)is,"stress_corr=",stress_correction
        do i=1,nbasis
          stress_par(is,i,1) = u_par(is,i)*stress_correction
          stress_par(is,i,2) = p_par(is,i)*stress_correction
          stress_per(is,i,1) = 0.0
          stress_per(is,i,2) = 0.0
          do j=1,nbasis
              stress_per(is,i,1) = stress_per(is,i,1)  &
            + xi*ky*ave_kx(i,j)*(1.5*p_tot(is,j)-0.5*p_par(is,j)) 
              stress_per(is,i,2) = stress_per(is,i,2)  &
            + xi*ky*ave_kx(i,j)*(1.5*q_tot(is,j)-0.5*q_par(is,j)) 
          enddo
        enddo
      enddo
!
!  compute the quasilinear weights for the fluxes
!
      do is=ns0,ns
        do j=1,3
          particle_weight(is,j)=0.0
          energy_weight(is,j)=0.0
          stress_par_weight(is,j)=0.0
          stress_tor_weight(is,j)=0.0
          exchange_weight(is,j)=0.0
        enddo
        do i=1,nbasis
          particle_weight(is,1) = particle_weight(is,1) &
          + REAL(xi*CONJG(phi(i))*n(is,i))
          energy_weight(is,1) = energy_weight(is,1) &
          + REAL(xi*CONJG(phi(i))*p_tot(is,i))
          do j=1,nbasis
              stress_par_weight(is,1) = stress_par_weight(is,1)  &
            + REAL(xi*CONJG(phi(i))*ave_c_par_par(i,j)*stress_par(is,j,1))           
              stress_tor_weight(is,1) = stress_tor_weight(is,1)  &
            + REAL(xi*CONJG(phi(i))*(ave_c_tor_par(i,j)*stress_par(is,j,1)+ave_c_tor_per(i,j)*stress_per(is,j,1)))
          enddo
!          stress_par_weight(is,1) = stress_par_weight(is,1)  &
!           + REAL(xi*CONJG(phi(i))*u_par(is,i))
!          stress_par_weight(is,1) = stress_par_weight(is,1)  &
!            + REAL(xi*CONJG(phi(i))*stress_par(is,i,1)*ave_c_par_par(1,1))           
!          stress_tor_weight(is,1) = stress_tor_weight(is,1)  &
!            + REAL(xi*CONJG(phi(i))*(ave_c_tor_par(1,1)*stress_par(is,i,1)+ave_c_tor_per(1,1)*stress_per(is,i,1)))
          exchange_weight(is,1) = exchange_weight(is,1) &
          + zs(is)*REAL(xi*freq_QL*CONJG(phi(i))*n(is,i))
          if(use_bper_in)then
            particle_weight(is,2) = particle_weight(is,2) &
            - vs(is)*REAL(xi*CONJG(psi(i))*u_par(is,i))
            energy_weight(is,2) = energy_weight(is,2) &
            - vs(is)*REAL(xi*CONJG(psi(i))*q_tot(is,i))
            exchange_weight(is,2) = exchange_weight(is,2) &
            - zs(is)*vs(is)*REAL(xi*freq_QL*CONJG(psi(i))*u_par(is,i))
            do j=1,nbasis
              stress_par_weight(is,2) = stress_par_weight(is,2)  &
              - REAL(xi*CONJG(psi(i))*ave_c_par_par(i,j)*stress_par(is,j,2))           
              stress_tor_weight(is,2) = stress_tor_weight(is,2)  &
              - REAL(xi*CONJG(psi(i))*(ave_c_tor_par(i,j)*stress_par(is,j,2)+ave_c_tor_per(i,j)*stress_per(is,j,2)))
            enddo
          endif
          if(use_bpar_in)then
            particle_weight(is,3) = particle_weight(is,3) &
            + REAL(xi*CONJG(bsig(i))*(1.5*p_tot(is,i)-0.5*p_par(is,i)))*taus(is)/zs(is)
            exchange_weight(is,3) = exchange_weight(is,3) &
            + taus(is)*REAL(CONJG(-xi*freq_QL*bsig(i))*(1.5*p_tot(is,i)-0.5*p_par(is,i)))
          endif
        enddo
!
        do j=1,3
          particle_weight(is,j) = as(is)*ky*particle_weight(is,j)/phi_norm
          energy_weight(is,j) = as(is)*taus(is)*1.5*ky*energy_weight(is,j)/phi_norm
          stress_par_weight(is,j) = mass(is)*as(is)*vs(is)*ky*stress_par_weight(is,j)/phi_norm
          stress_tor_weight(is,j) = sign_It_in*mass(is)*as(is)*vs(is)*ky*stress_tor_weight(is,j)/phi_norm
          exchange_weight(is,j) = as(is)*exchange_weight(is,j)/phi_norm
        enddo
      enddo
!
!   add the vpar shifts to the total  moments
!
      if(vpar_model_in.eq.0)then
        do is=ns0,ns
        do j=1,nbasis
          n(is,j) = n(is,j) + vpar_s(is)*(zs(is)/taus(is))*psi(j)
          u_par(is,j) = u_par(is,j) -(vpar_s(is)/vs(is))*(zs(is)/taus(is))*phi(j)
          p_par(is,j) = p_par(is,j) + vpar_s(is)*(zs(is)/taus(is))*psi(j)
          p_tot(is,j) = p_tot(is,j) + vpar_s(is)*(zs(is)/taus(is))*psi(j)
          q_par(is,j) = q_par(is,j) - 3.0*(vpar_s(is)/vs(is))*(zs(is)/taus(is))*phi(j)
          q_tot(is,j) = q_tot(is,j) -(5.0/3.0)*(vpar_s(is)/vs(is))*(zs(is)/taus(is))*phi(j)
        enddo
        enddo
      endif
! 
!  compute the density and temperature amplitude weights 
!    
      do is=ns0,ns
        N_weight(is)=0.0
        T_weight(is)=0.0
        U_weight(is)=0.0
        Q_weight(is)=0.0
        do i=1,nbasis
          temp(is,i) = p_tot(is,i) - n(is,i)
          N_weight(is) = N_weight(is) + REAL(n(is,i)*CONJG(n(is,i)))
          T_weight(is) = T_weight(is) + REAL(temp(is,i)*CONJG(temp(is,i)))
          U_weight(is) = U_weight(is) + REAL(u_par(is,i)*CONJG(u_par(is,i)))
          Q_weight(is) = Q_weight(is) + REAL(q_tot(is,i)*CONJG(q_tot(is,i)))
        enddo
        N_weight(is) = N_weight(is)/phi_norm
        T_weight(is) = T_weight(is)/phi_norm
        U_weight(is) = U_weight(is)/phi_norm
        Q_weight(is) = Q_weight(is)/phi_norm
      enddo
      v_weight = vnorm/phi_norm
      a_par_weight = psi_norm/phi_norm
      b_par_weight = bsig_norm/phi_norm
!
!      write(*,*) 'ns = ',ns
!      write(*,*) 'is  particle_weight   energy_weight'
!      do is=ns0,ns
!        write(*,100)is,particle_weight(is),energy_weight(is)
!      enddo
!      write(*,*) 'Q_i/Q_e = ',energy_weight(2)/energy_weight(1)
!      write(*,*) 'G_i/Q_i = ',particle_weight(2)/energy_weight(2)
!      write(*,*)
!
! 
! compute electron density-temperature phase 
      Ne_Te_phase = 0.0
      Ne_Te_cos = 0.0
      Ne_Te_sin = 0.0
      do i=1,nbasis
         Ne_Te_cos = Ne_Te_cos + REAL(CONJG(n(1,i))*temp(1,i))
         Ne_Te_sin = Ne_Te_sin + AIMAG(CONJG(n(1,i))*temp(1,i))
      enddo
      Ne_Te_phase = ATAN2(Ne_Te_sin,Ne_Te_cos)
!
! compute species density-temperature phase
      do is=ns0,ns
        Ns_Ts_phase(is) = 0.0
        Ns_Ts_cos = 0.0
        Ns_Ts_sin = 0.0
        do i=1,nbasis
           Ns_Ts_cos = Ns_Ts_cos + REAL(CONJG(n(is,i))*temp(is,i))
           Ns_Ts_sin = Ns_Ts_sin + AIMAG(CONJG(n(is,i))*temp(is,i))
        enddo
        Ns_Ts_phase(is) = ATAN2(Ns_Ts_sin,Ns_Ts_cos)
      enddo
!
      END SUBROUTINE get_QL_weights
!
!    
!***********************
! start of get_wavefunction
      SUBROUTINE get_wavefunction
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_sgrid
!
      IMPLICIT NONE
!
      INTEGER :: n,i,j,k,np,npi,j0,imax
      REAL :: dx,hp0
      REAL :: hp(nb,max_plot)
      REAL :: xp(max_plot)
!
! set up the theta-grid
! npi is the number of pi intervals for the plot from -npi Pi to +npi Pi
     npi=9  !npi <= 9 limited by max_plot = 2*npi*ms/8+1
     np = ms/8   ! number of points per 1/2 period = 16 for ms=12
     if(igeo.eq.0)then
       dx = REAL(npi)*2.0*pi/REAL(max_plot-1)
       do i=1,max_plot
         xp(i) = -REAL(npi)*pi + REAL(i-1)*dx
         plot_angle_out(i) = xp(i)        
       enddo
     else
! general geometry case 0<y<Ly one loop counterclockwise 
! y is the the straight field line coordinant of the Hermite basis
! t_s is the mapping of the original theta coordinate 0 < t_s < -2Pi with t_s(0)=0, t_s(ms)=-2 Pi
! this will be used for plotting, Note that t_s has the opposite sign to y.
       dx = 2.0*pi/(y(ms)*width_in)
       j0 = npi*np+1 ! the index of the midpoint
       xp(j0) = 0.0
       plot_angle_out(j0) = 0.0
       j = 0  ! j is the local index for one circuit poloidally
       k = 0  ! counts the number of 2 pi  intervals
       imax=np*npi
       do i=1,imax
           j = j+1
           if(j.gt.2*np)then
             j = j - 2*np
             k = k + 1
           endif
           xp(j0+i) = (REAL(k)*y(ms) + y(4*j))*dx  ! remember y is positive 0 <= y <= Ly
           xp(j0-i) = -(REAL(k+1)*y(ms) - y(ms-4*j))*dx ! ok for up/down assymetric cases
           plot_angle_out(j0+i) = -(REAL(k)*t_s(ms) + t_s(4*j))  ! remember t_s is negative  0 <= t_s <= -2 Pi
           plot_angle_out(j0-i) = REAL(k+1)*t_s(ms) - t_s(ms-4*j)
       enddo
     endif
 !     write(*,*)"t_s = ",(t_s(i),i=0,ms)
 !   do i=1,imax
 !     write(*,*)i,"xp=",xp(i),"tp=",plot_angle_out(i)
 !   enddo
! compute the hermite polynomials on the theta-grid using recursion
      hp0 = sqrt_two/pi**0.25
      do i=1,max_plot
       hp(1,i) = hp0*EXP(-xp(i)*xp(i)/2.0)
       hp(2,i) = xp(i)*sqrt_two*hp0*EXP(-xp(i)*xp(i)/2.0)
       if(nbasis.gt.2)then
         do j=3,nbasis
          hp(j,i) = xp(i)*SQRT(2.0/REAL(j-1))*hp(j-1,i) &
           - SQRT(REAL(j-2)/REAL(j-1))*hp(j-2,i)
         enddo
       endif
      enddo
!      compute the fields on the theta-grid
      do n=1,nmodes_out
       do k=1,3
!        write(*,*)"field_weight_out=",(field_weight_out(n,k,j),j=1,nbasis)
        do i=1,max_plot
         plot_field_out(n,k,i) = 0.0
         if((k.eq.1).or.(k.eq.2.and.use_bper_in).or.(k.eq.3.and.use_bpar_in))then
          do j=1,nbasis
           plot_field_out(n,k,i) = plot_field_out(n,k,i) &
           + field_weight_out(n,k,j)*hp(j,i)
          enddo
         endif
        enddo
       enddo
      enddo
!
      END SUBROUTINE get_wavefunction
!
