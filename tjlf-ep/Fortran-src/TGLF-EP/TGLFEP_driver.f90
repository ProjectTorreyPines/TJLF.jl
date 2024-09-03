!------------------------------------------------------------------
! TGLFEP_driver.f90
! 7.13.2016
! PURPOSE:
!  TGLFEP_driver calls TGLFEP_mainsub
!------------------------------------------------------------------

program TGLFEP_driver

  use mpi
  use TGLFEP_interface
  use TGLFEP_profile

!  use EXPRO_interface
  use expro

  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------
   
  integer :: id, np, color, key, ierr, STATUS(MPI_STATUS_SIZE)
  integer :: TGLFEP_COMM_IN

  logical :: iexist, factor_in_profile
  integer :: i,ii,i_r
  real,allocatable,dimension(:) :: factor, width, kymark_out, SFmin, SFmin_out
  real, allocatable, dimension(:) :: dndr_crit, dndr_crit_out, dpdr_crit, dpdr_crit_out
  real :: dpdr_scale
  real,allocatable,dimension(:) :: SFmin_out_n
  real,allocatable,dimension(:) :: SFmin_out_p
  character(3) :: str_r

  integer :: ir_min
  integer :: ir_max
  integer :: ir_dum_1
  integer :: ir_dum_2

  real, allocatable, dimension(:) :: dpdr_EP
  real :: dpdr_EP_max
  integer :: dpdr_EP_max_loc
  real :: n_at_max
  real :: dlnndr_at_max
  
  integer :: values(8)
  integer :: values_begin(8)
  integer :: values_end(8)

  logical, allocatable, dimension(:) :: l_accept_profile

  call MPI_INIT(ierr)
  TGLFEP_COMM_WORLD = MPI_COMM_WORLD

  call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

  np_3 = np
  id_3 = id

  l_print = .false.
  if (id .eq. 0) l_print = .true.

!  if (id .EQ. 0) print *, 'MPI environment initialized'

  !read input
  inquire(file='input.TGLFEP',exist=iexist)
  if(iexist) then
    open(unit=22,file='input.TGLFEP',status='old')
    read(22,'(I1)') process_in
    if (id .eq. 0) print *, 'Read process_in = ', process_in

    if(process_in .le. 1) read(22,'(I1)') mode_in
    if ((process_in .eq. 4).or.(process_in.eq.5).or.(process_in.eq.6)) then
      read(22,'(I1)') threshold_flag
      read(22,'(I2)') n_basis
      read(22,'(I1)') scan_method
      read(22,'(I1)') reject_i_pinch_flag
      read(22,'(I1)') reject_e_pinch_flag
      read(22,'(I1)') reject_th_pinch_flag
      read(22,'(I1)') reject_EP_pinch_flag
      read(22,'(I1)') reject_tearing_flag
      read(22,'(I1)') rotational_suppression_flag
!      read(22,'(I1)') reject_max_outer_panel_flag
      read(22,*) QL_ratio_thresh
      read(22,*) theta_sq_thresh
      read(22,*) q_scale
    endif

    read(22,'(I1)') l_write_wavefunction
    read(22,'(I1)') ky_model
    if (id .eq. 0) print *, 'Read ky_model = ', ky_model

    read(22,*) scan_n
    if (id .eq. 0) print *, 'Read scan_n = ', scan_n
    read(22,*) irs
    if (id .eq. 0) print *, 'Read irs = ', irs

    read(22,*) factor_in_profile
    if (id .eq. 0) print *, 'Read factor_in_profile = ', factor_in_profile
    allocate(factor(scan_n))
    allocate(factor_max_profile(scan_n))
    allocate(ir_exp(scan_n))
    if(factor_in_profile) then
      do i = 1,scan_n
        read(22,*) factor(i)
      enddo
      if (id .eq. 0) print *, 'Read factor(scan_n) = ', factor(scan_n)
    else
      read(22,*) factor_in
      if (id .eq. 0) print *, 'Read factor_in = ', factor_in
      factor(:) = factor_in
    endif
    factor_max_profile(:) = factor(:)

    read(22,*) width_in_flag
    if (id .eq. 0) print *, 'Read width_in_flag = ', width_in_flag
    allocate(width(scan_n))
    allocate(kymark_out(scan_n))
    if(width_in_flag) then
      do i = 1,scan_n
        read(22,*) width(i)
      enddo
      if (id .eq. 0) print *, 'Read width(scan_n) = ', width(scan_n)
    else
      width = 0.00
      read(22,*) width_min
      if (id .eq. 0) print *, 'Read width_min = ', width_min
      read(22,*) width_max
      if (id .eq. 0) print *, 'Read width_max = ', width_max
    endif

    read(22,*,end=100) input_profile_method                                       ! EMB 10/2018
    if (id .eq. 0) print *, 'Read input_profile_method = ', input_profile_method  !
    read (22,*) TGLFEP_nion
    read (22,*) is_EP
!    read (22,*) m_i        ! Main ion mass over deuterium mass
!    read (22,*) q_i        ! Main ion charge over e
!    read (22,*) m_EP       ! EP mass over deuterium mass
!    read (22,*) q_EP       ! EP charge over e
  ! EMB 10/2020
    read (22,*) l_real_units  
100 continue

    close(22)
  else
    print *, 'input.TGLFEP file not found'
    stop
  endif

!  if (id .EQ. 0) print *, 'input.TGLFEP read and variables allocated'

!  call read_input_profile                        !
  select case (input_profile_method)              ! EMB 10/2018
    case (1)                                      !
      call read_input_profile                     !
    case (2)                                      !
      call TGLFEP_read_EXPRO
    case (3)                                      !
      call TGLFEP_generate_input                  !
    case default                                  !
      call read_input_profile                     !
  end select                                      !

! EMB 5-26-2021
! For experimental profiles, adjust max scale factor to be constant across the domain. 
  allocate(dpdr_EP(nr))
  if (input_profile_method .eq. 2) then
    dpdr_EP = EXPRO_ni(is_EP,:)*EXPRO_Ti(is_EP,:)*(EXPRO_dlnnidr(is_EP,:)+EXPRO_dlntidr(is_EP,:))  ! A.U.
    dpdr_EP_max = maxval(abs(dpdr_EP))
    dpdr_EP_max_loc = maxloc(abs(dpdr_EP),1)
    n_at_max = EXPRO_ni(is_EP,dpdr_EP_max_loc)
    dlnndr_at_max = EXPRO_dlnnidr(is_EP,dpdr_EP_max_loc)
    if (process_in .ne. 5) then           ! Factor rescale disabled for process_in=5
      do ir = 1, scan_n
        factor(ir) = factor(ir)*dpdr_EP_max/abs(dpdr_EP(ir_exp(ir)))
      enddo
!    else
!      EXPRO_ni(is_EP,:) = n_at_max
!      EXPRO_dlnnidr(is_EP,:) = dlnndr_at_max
    endif
    factor_max_profile(:) = factor(:)
!    print *, 'factor profile', factor(:), 'dpdr_EP profile', dpdr_EP(:)
  endif

!  allocate(f_real(nr))
!  f_real(:) = 1.0

!  if (id .eq. 0) print *, zs 

!  if (id .EQ. 0) print *, 'read_input_profile executed'

  if(ky_model .eq. 0) then
    n_toroidal = 4
  else
    n_toroidal = 3
  endif

  if ((process_in .eq. 4).or.(process_in .eq. 5).or.(process_in .eq. 6)) then
    if(threshold_flag .eq. 0) then
      nn = 5
!      nn = 10
!      nn = 1
    else
      nn = 15
    endif
    allocate(factor_out(nn))
  endif

  !run mainsub
  key = id / (scan_n)
  color = id - key*(scan_n)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,TGLFEP_COMM_IN,ierr)

!  if (id .EQ. 0) print *, 'MPI communicator split.'

  if (input_profile_method .eq. 2) then
    ir = ir_exp(color+1)
!    print *, 'id: ', id, ' color: ', color, ' ir: ', ir 
  else
    ir = color + irs
  endif
  str_r = achar(ir/100+iachar("0"))//achar(mod(ir,100)/10+iachar("0"))//achar(mod(ir,10)+iachar("0"))
  write (suffix,'(A2,A3)') '_r', str_r

  factor_in = factor(color+1)
  if(width_in_flag) width_in = width(color+1)

  call date_and_time(values=values)  
  values_begin(:) = values(:)
  call TGLFEP_mainsub(TGLFEP_COMM_IN)
  if (id .eq. 0) print *, 'Finished TGLFEP_mainsub.'
  call date_and_time(values=values)
  values_end(:) = values(:)
  values(7) = modulo(values_end(7)-values_begin(7),60)
  values(6) = modulo(values_end(6)-values_begin(6),60)
  if (values_end(7) .lt. values_begin(7)) values(6) = values(6) - 1
  values(5) = modulo(values_end(5)-values_begin(5),24)
  if (values_end(6) .lt. values_begin(6)) values(5) = values(5) - 1

  if (id .eq. 0) print *, "Run time: ", values(5), ':', values(6), ':', values(7)

!  if (id .EQ. 0) print *, 'TGLFEP_mainsub executed.'
!  goto 1000

  !write output
  if(.not. width_in_flag) then !get width_out and kymark_out
    if(id .eq. 0) then
      width(1) = width_in
      kymark_out(1) = kymark
      do i = 1,scan_n-1
        call MPI_RECV(width(i+1),1,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,STATUS,ierr)
        call MPI_RECV(kymark_out(i+1),1,MPI_DOUBLE_PRECISION,i,i+scan_n,MPI_COMM_WORLD,STATUS,ierr)
      enddo
    else if(id .lt. scan_n) then
      call MPI_SEND(width_in,1,MPI_DOUBLE_PRECISION,0,id,MPI_COMM_WORLD,STATUS,ierr)
      call MPI_SEND(kymark,1,MPI_DOUBLE_PRECISION,0,id+scan_n,MPI_COMM_WORLD,STATUS,ierr)
    endif
  endif

  if(id .eq. 0) then

    open(unit=18,file='out.TGLFEP',status='replace')
    write(18,*) 'process_in = ',process_in

    if(process_in .le. 1) write(18,*) 'mode_in = ',mode_in
    if ((process_in .eq. 4) .or. (process_in .eq. 5) .or. (process_in .eq. 6)) write(18,*) 'threshold_flag = ',threshold_flag

    write(18,*) 'ky_model = ',ky_model

    write(18,*) '--------------------------------------------------------------'
    write(18,*) 'scan_n = ',scan_n
    write(18,*) 'irs = ',irs
    write(18,*) 'n_basis= ', n_basis
    write(18,*) 'scan_method= ', scan_method

    if(width_in_flag) then
      write(18,*) 'ir,  width'
      do i = 1,scan_n
        write(18,'(I3,F8.2)') irs+i-1,width(i)
      enddo
    else
      write(18,*) 'ir,  width,  kymark'
      do i = 1,scan_n
        write(18,'(I3,F8.2,F9.3)') irs+i-1,width(i),kymark_out(i)
      enddo
    endif

    write(18,*) '--------------------------------------------------------------'
    write(18,*) 'factor_in_profile = ',factor_in_profile
    if(factor_in_profile) then
      do i = 1,scan_n
        write(18,'(F7.2)') factor(i)
      enddo
    else
      write(18,*) 'factor_in = ',factor(1)
    endif

    write(18,*) 'width_in_flag = ',width_in_flag
    if(.not. width_in_flag) write(18,*) 'width_min = ',width_min,'width_max = ',width_max

    close(18)

  endif

  if ((process_in .eq. 4).or.(process_in .eq. 5).or.(process_in .eq. 6)) then !print out 'density threshold'
    if(id .eq. 0) then
      open(unit=11,file='out.TGLFEP',status='old',position='append')
      write(11,*) '**************************************************************'
      write(11,*) '************** The critical EP density gradient **************'
      write(11,*) '**************************************************************'

      allocate(SFmin(scan_n))    ! Only points run in simulations.
      allocate(SFmin_out(nr))    ! All points in profile input, with interpolations.
      allocate(dndr_crit(scan_n))
      allocate(dndr_crit_out(nr))
      allocate(dpdr_crit(scan_n))
      allocate(dpdr_crit_out(nr))
!      allocate(SFmin_out_n(nr))
!      allocate(SFmin_out_p(nr))
      allocate(l_accept_profile(scan_n))

      print *, "After second allocation series."
      if(threshold_flag .eq. 0) then

        SFmin(1) = factor_in
        print *, "Before MPI_RECV for factor_in."
        do i = 1,scan_n-1
          call MPI_RECV(factor_in,1,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,STATUS,ierr)
          SFmin(i+1) = factor_in
        enddo
        print *, "After MPI_RECV for factor_in."

        write(11,*) '--------------------------------------------------------------'
        write(11,*) 'SFmin'
!        write(11,10) SFmin

! Fill in all interior blanks with interpolation and determine first and last solution coordinate.
        call TGLFEP_complete_output(SFmin,SFmin_out,ir_min,ir_max,l_accept_profile)
! Fill in edge regions with maximum value tested (with still no mode found) for SFmin and SFmin_out.
        if ((ir_min-irs+1) .gt. 1) then
          SFmin(1:ir_min-irs) = factor_max_profile(1:ir_min-irs)
          if (irs .gt. 1) SFmin_out(1:irs-1) = factor_max_profile(1)
          SFmin_out(irs:ir_min-1) = factor_max_profile(1:ir_min-irs)
        endif
        if ((ir_max-irs+1) .lt. scan_n) then
          SFmin(ir_max-irs+2:scan_n) = factor_max_profile(ir_max-irs+2:scan_n)
          if (irs+scan_n-1 .lt. nr)  SFmin_out(irs+scan_n:nr) = factor_max_profile(scan_n)
          SFmin_out(ir_max+1:irs+scan_n-1) =factor_max_profile(ir_max-irs+2:scan_n)
        endif

! These extra lines make the absolute critical gradient in the extrapolated region fixed
! instead of making the scale factor fixed.

!        SFmin_out_n(:) = SFmin_out(:)
!        SFmin_out_n(1:ir_min)  = SFmin_out(1:ir_min)*EXPRO_ni(is_EP,ir_min)*EXPRO_dlnnidr(is_EP,ir_min) / &
!                                   (EXPRO_ni(is_EP,1:ir_min)*EXPRO_dlnnidr(is_EP,1:ir_min))
!        SFmin_out_n(ir_max:nr) = SFmin_out(ir_max:nr)*EXPRO_ni(is_EP,ir_max)*EXPRO_dlnnidr(is_EP,ir_max) / &
!                                   (EXPRO_ni(is_EP,ir_max:nr)*EXPRO_dlnnidr(is_EP,ir_max:nr))

        do i = 1, scan_n
          if (l_accept_profile(i)) then
            write(11,10) SFmin(i)
          else
            write(11,20) SFmin_out(i+irs-1), '   (*)'
          endif
        enddo

        if (input_profile_method .eq. 2) then
          dndr_crit(:) = 10000.
          do i = 1, scan_n
            if (SFmin(i) .lt. 9000.) then
!              dndr_crit(i) = SFmin(i)*EXPRO_ni(is_EP,i+irs-1)*EXPRO_dlnnidr(is_EP,i+irs-1)
              dndr_crit(i) = SFmin(i)*EXPRO_ni(is_EP,ir_exp(i))*EXPRO_dlnnidr(is_EP,ir_exp(i))
            else if ((i .lt. ir_min-irs+1) .or. (i .gt. ir_max-irs+1)) then
!              dndr_crit(i) = factor_max_profile(i)*EXPRO_ni(is_EP,i+irs-1)*EXPRO_dlnnidr(is_EP,i+irs-1)
              dndr_crit(i) = factor_max_profile(i)*EXPRO_ni(is_EP,ir_exp(i))*EXPRO_dlnnidr(is_EP,ir_exp(i))
            endif
          enddo
          call TGLFEP_complete_output(dndr_crit,dndr_crit_out,ir_dum_1,ir_dum_2,l_accept_profile)
          open(unit=22,file='alpha_dndr_crit.input',status='replace')
          write(22,'(A33)') 'Density critical gradient (10^19/m^4)'
          write(22,10) dndr_crit_out
          close(22)
        endif

! Generate here an alternatrive file where the pressure gradient is held fixed in the
! edge regions instead of the density gradient.

!        SFmin_out_p(:) = SFmin_out(:)
!        allocate(dpdr_EP(nr))
!        dpdr_EP(:) = EXPRO_ni(is_EP,:)*EXPRO_Ti(is_EP,:)*(EXPRO_dlnnidr(is_EP,:)+EXPRO_dlntidr(is_EP,:))*0.16022  ! 10 kPa/m
!        SFmin_out_p(1:ir_min)  = SFmin_out(1:ir_min)*dpdr_EP(ir_min) / dpdr_EP(1:ir_min)
!        SFmin_out_p(ir_max:nr) = SFmin_out(ir_max:nr)*dpdr_EP(ir_max) / dpdr_EP(ir_max:nr)

        if (input_profile_method .eq. 2) then
          dpdr_crit = 10000.
          dpdr_EP(:) = EXPRO_ni(is_EP,:)*EXPRO_Ti(is_EP,:)*(EXPRO_dlnnidr(is_EP,:)+EXPRO_dlntidr(is_EP,:))*0.16022  ! 10 kPa/m
          do i = 1, scan_n
            if (SFmin(i) .lt. 9000.) then
              if ((process_in .eq. 4).or.(process_in.eq.5)) then
                select case(scan_method)
                case(1)
                  dpdr_scale = SFmin(i)
                case(2)
                  dpdr_scale = (SFmin(i)*EXPRO_dlnnidr(is_EP,ir_exp(i))+EXPRO_dlntidr(is_EP,ir_exp(i))) / &
                                                     (EXPRO_dlnnidr(is_EP,ir_exp(i))+EXPRO_dlntidr(is_EP,ir_exp(i)))
                end select
!              dpdr_crit(i) = SFmin(i)*dpdr_EP(i+irs-1)
              dpdr_crit(i) = dpdr_scale*dpdr_EP(ir_exp(i))
!            else if ((i .lt. ir_min-irs+1) .or. (i .gt. ir_max-irs+1)) then
!              dpdr_crit(i) = factor_max_profile(i)*dpdr_EP(i+irs-1)
              endif
            endif
          enddo
          call TGLFEP_complete_output(dpdr_crit,dpdr_crit_out,ir_dum_1,ir_dum_2,l_accept_profile)
          open(unit=22,file='alpha_dpdr_crit.input',status='replace')
          write(22,'(A37)') 'Pressure critical gradient (10 kPa/m)'
          write(22,10) dpdr_crit_out
          close(22)
        endif

        if (input_profile_method .eq. 3) then
          open(unit=22,file='out.TGLFEP.generate',status='replace')
          write(22,*) 'a/Ln_j, a/LT_j, j=1,ns ;    ns=', ns
          write(22,*) 'rmin/a, rmaj/a, q, s_hat, dR/dr, kappa, s_kappa, delta, s_delta, zeta, s_zeta, betae:  '
          write(22,*) '(a/ne)*dn_EP/dr_crit'
          do i_r = 1, scan_n
            do i = 1, ns
              write(22,'(2E13.6)') rlns(i_r,i), rlts(i_r,i)
            enddo
            write(22,'(12E13.6,A3)') rmin(i_r), rmaj(i_r), q(i_r), shear(i_r), shift(i_r), kappa(i_r), s_kappa(i_r), delta(i_r), s_delta(i_r), zeta(i_r), s_zeta(i_r), betae(i_r), ' : '
            write(22,'(E13.6)') SFmin(i_r)*as(i_r,ns)*rlns(i_r,ns)
          enddo
          close(22)
        endif

        deallocate(dpdr_EP)

        write(11,*) '--------------------------------------------------------------'
        write(11,*) 'The EP density threshold n_EP/n_e (%) for gamma_AE = 0'
        do i = 1,scan_n
          write(11,10) SFmin(i)*as(irs+i-1,is)*100. !percent
        enddo

        write(11,*) '--------------------------------------------------------------'
        write(11,*) 'The EP beta crit (%) = beta_e*(n_EP_th/n_e)*(T_EP/T_e)'
        do i = 1,scan_n
          if(geometry_flag .eq. 0) then
            write(11,10) SFmin(i)*betae(irs+i-1)*100.*as(irs+i-1,is)*taus(irs+i-1,is) !percent
          else
            write(11,10) SFmin(i)*betae(irs+i-1)*100.*as(irs+i-1,is)*taus(irs+i-1,is)*kappa(irs+i-1)**2 !percent
          endif
        enddo

        deallocate(SFmin)
        deallocate(SFmin_out)
        deallocate(dndr_crit)
        deallocate(dndr_crit_out)
        deallocate(dpdr_crit)
        deallocate(dpdr_crit_out)
!        deallocate(SFmin_out_n)
!        deallocate(SFmin_out_p)
        deallocate(l_accept_profile)

        if (process_in .eq. 4) then
          write(11,*) '--------------------------------------------------------------'
          write(11,*) 'The scale factor for EP density threshold at each n:'
          write(11,*) (factor_out(ii),ii=1,nn)
          do i = 1,scan_n-1
            call MPI_RECV(factor_out,nn,MPI_DOUBLE_PRECISION,i,i+scan_n,MPI_COMM_WORLD,STATUS,ierr)
            write(11,*) (factor_out(ii),ii=1,nn)
          enddo
        endif

      else    ! threshold_flag /= 0

        write(11,*) '--------------------------------------------------------------'
        write(11,*) 'a/Ln_EP multiplied by: 0.1, 0.2, ..., 1.5'
        write(11,*) 'The scale factor for EP density threshold:'
        write(11,*) (factor_out(ii),ii=1,nn)
        do i = 1,scan_n-1
          call MPI_RECV(factor_out,nn,MPI_DOUBLE_PRECISION,i,i+scan_n,MPI_COMM_WORLD,STATUS,ierr)
          write(11,*) (factor_out(ii),ii=1,nn)
        enddo

      endif  ! threshold flag = 0 or /=0

      close(11)
    else if(id .lt. scan_n) then     ! id /=0
      if(threshold_flag .eq. 0) then
        call MPI_SEND(factor_in,1,MPI_DOUBLE_PRECISION,0,id,MPI_COMM_WORLD,STATUS,ierr)
      endif

      call MPI_SEND(factor_out,nn,MPI_DOUBLE_PRECISION,0,id+scan_n,MPI_COMM_WORLD,STATUS,ierr)
    endif     ! id =0 or /=0

  endif       ! process_in = 4 or 5

  if ((process_in .eq. 7) .and. (id .eq. 0)) then
    open (unit=6,file='out.DEP_trace',status='replace')
    write(6,*) '( Trace DEP [m^2/s], chi_eff [m^2/s], DEP/chi_eff )'
    do jtscale = 1, jtscale_max
      write(6,*) 'Temeprature scale factor =', jtscale*tscale_interval
      do i = 1, nr
        write(6,*) DEP_trace_scan(jtscale,i), DEP_trace_scan(jtscale,i)/(QL_ratio_scan(jtscale,i)+1.0E-10), QL_ratio_scan(jtscale,i)
      enddo  ! i
    enddo  !  jtscale
    close (unit=6)
    open (unit=6,file='out.DEP_trace_interp',status='replace')
    write(6,*) '( Trace DEP [m^2/s], chi_eff [m^2/s], DEP/chi_eff )'
    do jtscale = 1, jtscale_max
      DEP_trace(:) = DEP_trace_scan(jtscale,:)
      QL_ratio(:) = QL_ratio_scan(jtscale,:)
      call TGLFEP_interpolate(DEP_trace,DEP_trace_complete)
      call TGLFEP_interpolate(QL_ratio,QL_ratio_complete) 
      do i = 1, nr
        write(6,*) DEP_trace_complete(i), DEP_trace_complete(i)/(QL_ratio_complete(i)+1.0E-10), QL_ratio_complete(i)
      enddo  ! i
    enddo  ! jtscale
    close (unit=6)
  endif

  if ((process_in .eq. 4) .or. (process_in .eq. 5) .or. (process_in .eq. 6)) deallocate(factor_out)
  deallocate(width)
  deallocate(kymark_out)
  deallocate(factor)
  deallocate(factor_max_profile)

  if (process_in .eq. 7) then
    deallocate(DEP_trace)
    deallocate(QL_ratio)
    deallocate(DEP_trace_complete)
    deallocate(QL_ratio_complete)
    deallocate(chi_gB)
  endif

  if(id .eq. 0) call dump_profile

  call deallocate_profile
  deallocate(f_real)
  deallocate(ir_exp)
  call MPI_FINALIZE(ierr)

10 format(F12.4)
20 format(F12.4,A6)

1000 continue

end program TGLFEP_driver
