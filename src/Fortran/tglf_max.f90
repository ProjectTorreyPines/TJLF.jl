!
      SUBROUTINE tglf_max
!
      USE tglf_global
      USE tglf_dimensions
      USE tglf_species
      USE tglf_pkg
!
      IMPLICIT NONE
      LOGICAL :: save_iflux
      LOGICAL :: save_bper
      LOGICAL :: save_bpar
      INTEGER :: nt,i,is,imax
      INTEGER :: save_nbasis
      INTEGER :: save_ibranch
      REAL :: width_max=1.65
      REAL :: width_min=0.15
      REAL :: tp,dt
      REAL :: t1,t2,tm,g1,g2,gm
      REAL :: tmax,tmin,gamma_max,gmax
      REAL :: save_width,dtmin
      REAL :: gamma_n(nt0),freq_n(nt0),width_n(nt0)
!      REAL :: save_vexb_shear
      REAL :: wgp_max,width_p_max
      REAL :: kyi,ft2
!
      CALL tglf_setup_geometry
!
!      do i=1,nmodes_in
!        gamma_reference_kx0(i)=0.0
!        freq_reference_kx0(i)=0.0
!      enddo
      save_iflux = iflux_in
      save_ibranch = ibranch_in
      save_nbasis = nbasis_max_in
      save_width = width_in
      save_bper = use_bper_in
      save_bpar = use_bpar_in
!      save_vexb_shear = vexb_shear_s
!      if(alpha_quench_in.eq.0.0)vexb_shear_s = 0.0
      ibranch_in = -1
      width_min = width_min_in
      width_max = ABS(width_in)
!
!      write(*,*)"R_unit=",R_unit,"q_unit=",q_unit
!      write(*,*)ns0,ns,ky_s
    if(alpha_p_in.gt.0.0)then
      do is=ns0,ns
        kyi = ky_s*SQRT(taus(is)*mass(is))/ABS(zs(is))
        wgp_max = ABS((taus(is)/zs(is))*vpar_shear_s(is)/vs(is))*ky_s/(1+kyi**2)
        width_p_max = 3.6*vs(is)/(sqrt_two*R_unit*q_unit*MAX(wgp_max,0.001))
        width_p_max=MAX(width_p_max,0.1)
         if(width_p_max.lt.width_min_in)then
          width_min = width_p_max
        endif
      enddo
    endif
!        kyi = ky_s*SQRT(taus(2)*mass(2))/ABS(zs(2))
!        wgp_max = ABS(vpar_shear_in(2)/vs(2))*kyi/(1+kyi**2)
!        width_p_max = 3.6/(sqrt_two*R_unit*q_unit*MAX(wgp_max,0.001))
!        width_p_max=MAX(width_p_max,0.01)
!         if(width_p_max.lt.width_min_in)then
!          width_min = width_p_max
!        endif
!      write(*,*)ky," width_p_max = ", width_p_max,width_min
!
! for ibranch_in > 0 the most unstable positive frequency mode is stored 
! in gamma_out(1) and the most unstable negative frequency mode 
! (ion diamagnetic drift direction) is stored in gamma_out(2)
! for ibranch_in = -1 the unstable modes in rank order are stored in
! gamma_out(i) i=1,nmodes_in
!
      iflux_in=.FALSE.
      if(nbasis_min_in.ne.0)then
        nbasis = nbasis_min_in
      endif
      if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)then
        use_bper_in = .false.
        use_bpar_in = .false.
      endif
!       write(*,*)"nbasis = ",nbasis
      tmin=LOG10(width_min)
      tmax=LOG10(width_max)
      nt=nwidth_in
      dtmin=(tmax-tmin)/REAL(nt-1)
      if(use_bisection_in)nt=5
!
!      open(2,file = 'width.dat',status='unknown')
      dt = (tmax-tmin)/REAL(nt-1)
      tp = tmin
      do i=1,nt
       gamma_n(i)=0.0
       freq_n(i)=0.0
       width_n(i)=0.0
       tp = tmin + REAL(i-1)*dt
       width_in = 10.0**tp
!       write(*,*)"width_in = ",width_in
       new_width = .TRUE.
       call tglf_LS
!       write(*,*)i,width_in,gamma_out(1),freq_out(1),ft
       width_n(i)=width_in
       gamma_n(i) = gamma_out(1)
       freq_n(i) = freq_out(1)
      enddo
!      close(2)
! find the global maximum
      gamma_max=gamma_n(nt)
      imax=nt
      do i=nt-1,1,-1
        if(gamma_n(i).gt.gamma_max)then
          gamma_max=gamma_n(i)
          imax=i
        endif
!        write(*,*)"debug",i,imax,gamma_max,gamma_n(i)
      enddo
      width_in = width_n(imax)
!
      if(use_bisection_in.and.gamma_max.gt.0.0)then
!
! use bounded bisection search to refine width
!
         if(imax.eq.1)then
! maximum is against bottom width
           g1 = gamma_n(1)
           t1 = tmin
           g2 = gamma_n(2)
           t2 = LOG10(width_n(2))
           tp = (t2+t1)/2.0    
           width_in = 10.0**tp
           new_width = .TRUE.
           call tglf_LS
           gm = gamma_out(1)
           tm = tp
         elseif(imax.eq.nt)then
! maximum is against top width
           g1 = gamma_n(nt-1)
           t1 = LOG10(width_n(nt-1))
           g2 = gamma_n(nt)
           t2 = tmax           
           tp = (t2+t1)/2.0    
           width_in = 10.0**tp
           new_width = .TRUE.
           call tglf_LS
           gm = gamma_out(1)
           tm = tp
         else
! maximum is away from boundaries
           g1 = gamma_n(imax-1)
           t1 = LOG10(width_n(imax-1))
           g2 = gamma_n(imax+1)
           t2 = LOG10(width_n(imax+1))
           gm = gamma_n(imax)
           tm = LOG10(width_n(imax))
         endif
! start bisection search
         dt=(t2-t1)/2.0
!         write(*,*)"dtmin=",dtmin,"tmin=",tmin,"tmax=",tmax
         do while(dt.gt.dtmin)
           dt=dt/2.0
           gmax = MAX(gm,MAX(g1,g2))
!           write(*,*)"dt=",dt,gmax
!           write(*,*)g1,gm,g2
!           write(*,*)t1,tm,t2
!
           if(g1.eq.gmax)then  
             if(t1.gt.tmin)then
! shift past t1 and compute new g1,t1 
               tp = t1-dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               tm = t1
               gm = g1
               g1 = gamma_out(1)
               t1 = tp
! compute new g2,t2
               tp = tm+dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               g2 = gamma_out(1)
               t2 = tp
             else   ! t1 at tmin
! shrink towards t1
               tp = t1 + dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               g2 = gm
               t2 = tm
               gm = gamma_out(1)
               tm = tp               
             endif
           elseif(g2.eq.gmax)then
             if(t2.lt.tmax)then
! shift past t2 and compute new g2,t2
               tp = t2 + dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               gm = g2
               tm = t2
               g2 = gamma_out(1)
               t2 = tp
! compute new g1,t1
               tp = tm - dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               g1 = gamma_out(1)
               t1 = tp
             else  ! t2 at tmax 
! shrink towards t2
               tp = t2 - dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               g1 = gm
               t1 = tm
               gm = gamma_out(1)
               tm = tp               
             endif
           else  ! gm.eq.gmax
! compute new g1,t1 and g2,t2 closer to gm,tm
             tp = tm - dt
             width_in = 10.0**tp
             new_width = .TRUE.
             call tglf_LS
             g1 = gamma_out(1)
             t1 = tp
!
             tp = tm + dt
             width_in = 10.0**tp
             new_width = .TRUE.
             call tglf_LS
             g2 = gamma_out(1)
             t2 = tp
           endif           
         enddo  ! end of bisection search main loop
! find final maximum
         gmax=gm
         tp=tm
         if(g1.gt.gmax)then
          gmax=g1
          tp=t1
         endif
         if(g2.gt.gmax)then
          gmax=g2
          tp=t2
         endif
         gamma_max = gmax
         width_in = 10.0**tp        
      endif ! done with bisection search
!      write(*,*)"gamma_max=",gamma_max,"width_in=",width_in
       gamma_nb_min_out = gamma_max
!  reset ibranch_in
       ibranch_in = save_ibranch
!
       if(gamma_max.ne.0.0)then
! refine eigenvalue with more basis functions
         nbasis = save_nbasis
         if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)then
          use_bper_in = save_bper
          use_bpar_in = save_bpar
         endif
!         write(*,*)"nbasis=",nbasis
         iflux_in=save_iflux
         new_width=.TRUE.
         call tglf_LS
         if(ibranch_in.eq.-1)then
! check for inward ballooning modes
!         write(*,*) ky_s,use_inboard_detrapped_in
!         write(*,*)"modB_test = ",modB_test
!         write(*,*)"ft_test = ",ft_test
           if(use_inboard_detrapped_in.and.ft_test.gt.modB_test)then
             fts(:) =  ft_min
!           write(*,*)"changed ft",ft
             new_geometry = .FALSE.
             new_width = .FALSE.
             new_matrix = .TRUE.
             call tglf_LS
           endif
!           if(alpha_quench_in.eq.0.0)then
!             if(save_vexb_shear.ne.0.0.or.wgp_max.ne.0.0)then
!               do i=1,nmodes_out
!                 gamma_reference_kx0(i) = gamma_out(i)
!                 freq_reference_kx0(i) = freq_out(i)
!               enddo
!               vexb_shear_s = save_vexb_shear
!               iflux_in=save_iflux
!               new_width=.TRUE.
!               call tglf_LS
!             endif
!           endif
         endif  ! ibranch_in == -1
         gamma_max = MAX(gamma_out(1),gamma_out(2))  ! works for both ibranch_in cases
!
!        write(*,*)"width_p_max = ", width_p_max,width_in,gamma_max
!        write(*,*)ky,width_in,gamma_out(1),freq_out(1)
!        write(*,*)" maximum gamma for nbasis = ",nbasis_max_in
!        write(*,*)"gamma_out(1) = ",gamma_out(1)
!        write(*,*)"freq_out(1) = ",freq_out(1)
!        write(*,*)"gamma_out(2) = ",gamma_out(2)
!        write(*,*)"freq_out(2) = ",freq_out(2)
!        gamma_n(igamma)=gamma_out(branch)
!        freq_n(igamma) = freq_out(branch)
       endif
!
       if(gamma_max.eq.0.0)then
         width_in=save_width
         if(sat_rule_in.eq.2 .OR. sat_rule_in.eq.3)then
          use_bper_in = save_bper
          use_bpar_in = save_bpar
         endif
         do i=1,nmodes_in
          gamma_out(i)=0.0
          freq_out(i)=0.0
 !         gamma_reference_kx0(i)=0.0
 !         freq_reference_kx0(i)=0.0
         enddo
       endif
!
       nbasis=save_nbasis
       iflux_in=save_iflux
!       vexb_shear_s = save_vexb_shear
!
      END SUBROUTINE tglf_max   
