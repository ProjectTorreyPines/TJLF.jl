!------------------------------------------------------------
! TGLFEP_kwscale_scan.f90
!
! PURPOSE:
!  Calculate the EP density threshold gradient
!  but with floating values of scalefactor,
!  ky, and mode width.
!------------------------------------------------------------

subroutine TGLFEP_kwscale_scan

  use mpi
  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface 
  use TGLFEP_profile

  implicit none

  integer :: id,np,ierr,STATUS(MPI_STATUS_SIZE)
  integer :: id_world, np_world
  integer :: i,n,k
  integer :: k_max = 4
  logical :: iexist
  logical :: l_write_out = .true.
  real :: f0,f1,ft
  real :: w0, w1
  real :: kyhat0, kyhat1            ! kyhat = ky / rho_EP (ky_model=2) or 
                                    ! kyhat = n_toroidal (ky_model=1)
  integer, parameter :: nfactor = 10
  integer, parameter :: nefwid  = 10    ! Distinct from nwidth used in alternate width finder
  integer, parameter :: nkyhat  = 5    ! Distinct from nn used when process_in/=5
  integer :: nkwf = nfactor*nefwid*nkyhat
  real, dimension(nfactor) :: factor
  real, dimension(nefwid) :: efwid
  real, dimension(nkyhat) :: kyhat 

  integer, dimension(nkyhat,nefwid) :: imark
  integer :: imark_min
  integer :: imark_ref

  integer :: ifactor
  integer :: iefwid
  integer :: ikyhat
  integer :: ifactor_write
  integer :: iefwid_write
  integer :: ikyhat_write
  integer :: iefwid_mark
  integer :: ikyhat_mark

  integer :: ikymin
  integer :: ikymax
  integer :: iefwmin
  integer :: iefwmax

  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: lkeep_i
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: lkeep_i_out
  logical, dimension(nkyhat,nefwid) :: lkeep_ref
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: ltearing_i
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: ltearing_i_out
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_th_pinch_i
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_th_pinch_i_out
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_i_pinch_i
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_i_pinch_i_out
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_e_pinch_i
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_e_pinch_i_out
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_EP_pinch_i
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_EP_pinch_i_out
!  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_max_outer_panel_i
!  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_max_outer_panel_i_out
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_QL_ratio_i
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_QL_ratio_i_out
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_theta_sq_i
  logical, dimension(nkyhat,nefwid,nfactor,nmodes) :: l_theta_sq_i_out

  real,dimension(nkyhat,nefwid,nfactor,nmodes) :: growthrate,growthrate_out, &
                                                  frequency,frequency_out
  real,dimension(nmodes) :: g,f
  real :: gmark,fmark

  real :: delw
  real :: delky
  real :: delf
  real, dimension(nkyhat,nefwid) :: f_guess
  real :: f_guess_mark
  real :: f_g1
  real :: f_g2
  real :: gamma_g1
  real :: gamma_g2
  real :: gamma_mark_i_1(nkyhat,nefwid)
  real :: gamma_mark_i_2(nkyhat,nefwid)
  real :: f_mark_i(nkyhat,nefwid)
  real :: kyhat_max
  real :: kyhat_min

  real :: ky_write

  character(7)  :: str_sf
  character(4), dimension(nmodes)  :: keep_label

  integer :: values(8)

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)
  call MPI_COMM_RANK(TGLFEP_COMM_WORLD,id_world,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM_WORLD,np_world,ierr)

  kyhat_min = 0.0
  kyhat_max = 1.0

  growthrate       = 0.0
  growthrate_out   = 0.0
  frequency        = 0.0
  frequency_out    = 0.0

  lkeep_i(:,:,:,:)     = .true.
  lkeep_i_out(:,:,:,:) = .true.
  ltearing_i(:,:,:,:)     = .false.
  ltearing_i_out(:,:,:,:) = .false.
  l_th_pinch_i(:,:,:,:)     = .false.
  l_th_pinch_i_out(:,:,:,:) = .false.
  l_i_pinch_i(:,:,:,:)     = .false.             
  l_i_pinch_i_out(:,:,:,:) = .false.
  l_e_pinch_i(:,:,:,:)     = .false.             
  l_e_pinch_i_out(:,:,:,:) = .false.
  l_EP_pinch_i(:,:,:,:) = .false.
  l_EP_pinch_i_out(:,:,:,:) = .false.
!  l_max_outer_panel_i(:,:,:,:) = .false.
!  l_max_outer_panel_i_out(:,:,:,:) = .false.
  l_QL_ratio_i(:,:,:,:) = .false.
  l_QL_ratio_i_out(:,:,:,:) = .false.
  l_theta_sq_i(:,:,:,:) = .false.
  l_theta_sq_i_out(:,:,:,:) = .false.

  f0 = 0.0
  f1 = factor_in
  w0 = width_min
  w1 = width_max
  kyhat0 = kyhat_min
  kyhat1 = kyhat_max
!  print *, "factor_in in scalefactor (1)", factor_in

  ikyhat_write = floor(nkyhat/2.)
  iefwid_write = floor(nefwid/2.)
  ifactor_write = nfactor
  do k = 1,k_max

    do i = 1,nfactor
      factor(i) = (f1-f0)/nfactor*i+f0
    enddo
    do i = 1, nefwid
      efwid(i) = (w1-w0)/(nefwid-1)*(i-1) + w0
    enddo
    do i = 1, nkyhat
      if (kyhat0 .eq. 0.) then
        kyhat(i) = (kyhat1-kyhat0)/nkyhat*i + kyhat0
      else
        kyhat(i) = (kyhat1-kyhat0)/(nkyhat-1)*(i-1) + kyhat0  ! include lower bound if /= 0
      endif
    enddo

    do i = 1+id, nkwf, np
      l_wavefunction_out = 0
      
      ikyhat = floor((i-1.)/(nefwid*nfactor))+1
      iefwid = floor(1.0*mod(i-1,nefwid*nfactor)/nfactor)+1
      ifactor = mod(i-1,nfactor)+1

      factor_in  = factor(ifactor)
      kyhat_in = kyhat(ikyhat)
      width_in = efwid(iefwid)
!      if (id .eq. 0) print *, ir, ikyhat, iefwid, ifactor, factor_in, n_toroidal, width_in 
      str_sf =  achar(mod(floor(factor_in/100.),10)+iachar("0"))  //  &
                achar(mod(floor(factor_in/10.),10)+iachar("0"))   //  &
                achar(mod(floor(factor_in),10)+iachar("0"))       //  &
                '.'                                               //  &
                achar(mod(floor(10*factor_in),10)+iachar("0"))    //  &
                achar(mod(floor(100*factor_in),10)+iachar("0"))   //  &
                achar(mod(floor(1000*factor_in),10)+iachar("0"))

      write(str_wf_file,'(A16,A5,A3,A7)') 'out.wavefunction', trim(suffix), '_sf', str_sf

      if ( (l_write_wavefunction .eq. 1) .and.                     &
           (ikyhat .eq. ikyhat_write)    .and.                     &
           (iefwid .eq. iefwid_write)    .and.                     &
           (ifactor .eq. ifactor_write)  .and.                     &
           (k .eq. k_max) )   l_wavefunction_out = 1

      call TGLFEP_ky
      call date_and_time(values=values)
!   if (id==0) print 30, id,id_2,id_3, 'after TGLFEP_ky at i=',i, 'time  ', values(5),values(6),values(7)
!    print *, 'In scalefactor (2.25) factor=', factor
!    print *, 'proc_world: ', id_world, 'proc_loc: ', id, '  suffix:', suffix
      do n=1,nmodes
        growthrate(ikyhat,iefwid,ifactor,n) = get_growthrate(n)
        frequency(ikyhat,iefwid,ifactor,n)  = get_frequency(n)
        lkeep_i(ikyhat,iefwid,ifactor,n) = lkeep(n)
        ltearing_i(ikyhat,iefwid,ifactor,n) = ltearing(n)
        l_th_pinch_i(ikyhat,iefwid,ifactor,n) = l_th_pinch(n)
        l_i_pinch_i(ikyhat,iefwid,ifactor,n) = l_i_pinch(n)
        l_e_pinch_i(ikyhat,iefwid,ifactor,n) = l_e_pinch(n)
        l_EP_pinch_i(ikyhat,iefwid,ifactor,n) = l_EP_pinch(n)
!        l_max_outer_panel_i(ikyhat,iefwid,ifactor,n) = l_max_outer_panel(n)
        l_QL_ratio_i(ikyhat,iefwid,ifactor,n) = l_QL_ratio(n)
        l_theta_sq_i(ikyhat,iefwid,ifactor,n) = l_theta_sq(n)
      enddo

    enddo  !  i (kyhat, efwid, factor loop)

    call MPI_BARRIER(TGLFEP_COMM,ierr)
    call date_and_time(values=values)
!    if (id_3==0) print 50, 'All after MPI_BARRIER at ', values(5),values(6),values(7)

    call MPI_ALLREDUCE(growthrate,                     &
                       growthrate_out,                 &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_DOUBLE_PRECISION,           &
                       MPI_SUM,                        &
                       TGLFEP_COMM,                    &
                       ierr)
    
    call MPI_ALLREDUCE(frequency,                      &
                       frequency_out,                  &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_DOUBLE_PRECISION,           &
                       MPI_SUM,                        &
                       TGLFEP_COMM,                    &
                       ierr)

    call MPI_ALLREDUCE(lkeep_i,                        &
                       lkeep_i_out,                    &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_LOGICAL,                    &
                       MPI_LAND,                       &
                       TGLFEP_COMM,                    &
                       ierr)

    call MPI_ALLREDUCE(ltearing_i,                     &
                       ltearing_i_out,                 &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_LOGICAL,                    &
                       MPI_LOR,                        &
                       TGLFEP_COMM,                    &
                       ierr)

    call MPI_ALLREDUCE(l_th_pinch_i,                   &
                       l_th_pinch_i_out,               &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_LOGICAL,                    &
                       MPI_LOR,                        &
                       TGLFEP_COMM,                    &
                       ierr)

    call MPI_ALLREDUCE(l_i_pinch_i,                    &
                       l_i_pinch_i_out,                &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_LOGICAL,                    &
                       MPI_LOR,                        &
                       TGLFEP_COMM,                    &
                       ierr)

    call MPI_ALLREDUCE(l_e_pinch_i,                    &
                       l_e_pinch_i_out,                &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_LOGICAL,                    &
                       MPI_LOR,                        &
                       TGLFEP_COMM,                    &
                       ierr)

    call MPI_ALLREDUCE(l_EP_pinch_i,                   &
                       l_EP_pinch_i_out,               &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_LOGICAL,                    &
                       MPI_LOR,                        &
                       TGLFEP_COMM,                    &
                       ierr)

!    call MPI_ALLREDUCE(l_max_outer_panel_i,            &
!                       l_max_outer_panel_i_out,        &
!                       nkyhat*nefwid*nfactor*nmodes,   &
!                       MPI_LOGICAL,                    &
!                       MPI_LOR,                        &
!                       TGLFEP_COMM,                    &
!                       ierr)

    call MPI_ALLREDUCE(l_QL_ratio_i,                   &
                       l_QL_ratio_i_out,               &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_LOGICAL,                    &
                       MPI_LOR,                        &
                       TGLFEP_COMM,                    &
                       ierr)

    call MPI_ALLREDUCE(l_theta_sq_i,                   &
                       l_theta_sq_i_out,               &
                       nkyhat*nefwid*nfactor*nmodes,   &
                       MPI_LOGICAL,                    &
                       MPI_LOR,                        &
                       TGLFEP_COMM,                    &
                       ierr)

!    call date_and_time(values=values)
!    if (id==0) print 40, id, id_2, id_3, 'after MPI_ALLREDUCE at', values(5),values(6),values(7)

!    print *, 'In scalefactor (2.5) factor=', factor
    imark(:,:) = nfactor+1
    do ikyhat = 1, nkyhat
      do iefwid = 1, nefwid
        do ifactor = 1, nfactor
          do n = 1,nmodes
            if (lkeep_i_out(ikyhat,iefwid,ifactor,n)) then
              imark(ikyhat,iefwid) = ifactor
              exit
            endif
          enddo ! n (modes iterator)
          if(imark(ikyhat,iefwid) .le. nfactor) exit
        enddo ! ifactor
      enddo ! iefwid
    enddo ! ikyhat

    imark_min = nfactor+1
    do ikyhat = 1, nkyhat
      do iefwid = 1, nefwid
        imark_min = min(imark(ikyhat,iefwid),imark_min)
      enddo
    enddo

    fmark = 1.0E20
    gmark = 0.0
    f_guess_mark = 1.0E20
    if (imark_min .le. nfactor) then   ! (mode found)
      ikyhat_mark = nkyhat
      iefwid_mark = nefwid
      ikyhat_write = floor(nkyhat/2.)
      iefwid_write = floor(nefwid/2.)
      do ikyhat = 1, nkyhat
        do iefwid = 1, nefwid
          lkeep_ref(ikyhat,iefwid) = .false.
          imark_ref = nfactor-1
          if (imark(ikyhat,iefwid) .lt. nfactor) then
            imark_ref = imark(ikyhat,iefwid)
            lkeep_ref(ikyhat,iefwid) = .true.
          else
!            imark_ref = imark(ikyhat,iefwid) - 1
            if (imark(ikyhat,iefwid) .eq. nfactor) lkeep_ref(ikyhat,iefwid) = .true.
          endif
          f_g1 = factor(imark_ref)
          f_g2 = factor(imark_ref+1)
!          gamma_g1 = -2.0
!          gamma_g2 = -1.0
          gamma_g1 = 0.0
          gamma_g2 = 100.0
          do n = 1, nmodes
            if (lkeep_i_out(ikyhat,iefwid,imark_ref,n)) gamma_g1 = max(growthrate_out(ikyhat,iefwid,imark_ref,n),gamma_g1)
            if (lkeep_i_out(ikyhat,iefwid,imark_ref+1,n)) gamma_g2 = max(growthrate_out(ikyhat,iefwid,imark_ref+1,n),gamma_g2)
          enddo
!          f_guess(ikyhat,iefwid) = (gamma_g1*f_g2-gamma_g2*f_g1)/(gamma_g1-gamma_g2)  ! For gamma=0
          f_guess(ikyhat,iefwid) = f_g1 + (gamma_thresh-gamma_g1)*(f_g2-f_g1)/(gamma_g2-gamma_g1)  ! For gamma=gamma_thresh
          gamma_mark_i_1(ikyhat,iefwid) = gamma_g1
          gamma_mark_i_2(ikyhat,iefwid) = gamma_g2
          if (imark(ikyhat,iefwid) .lt. nfactor) then
            f_mark_i(ikyhat,iefwid) = f_g1
          else
            f_mark_i(ikyhat,iefwid) = f_g2
          endif
        enddo ! iefwid
      enddo  ! ikyhat
!      f_guess_mark = 10000.0
!      gmark = 0.
!      fmark = 1.0E20
      ikymin = nkyhat
      ikymax = 1
      iefwmin = nefwid
      iefwmax = 1
      do ikyhat = 1, nkyhat
        do iefwid = 1, nefwid
          if ( (lkeep_ref(ikyhat,iefwid)) .and.                                              &
!               (f_guess(ikyhat,iefwid) .gt. 0.0) .and.       &
!               (f_guess(ikyhat,iefwid) .lt. f_guess_mark) )  &
!               (gamma_mark_i_1(ikyhat,iefwid) .gt. gmark) .and.                              &
!               (gamma_mark_i_1(ikyhat,iefwid) .lt. 0.95*gamma_mark_i_2(ikyhat,iefwid)) .and.  &   ! Looking for (dgamma/d(dnEP/dr))>0
               (f_mark_i(ikyhat,iefwid) .le. fmark) )                                        &
          then
!              ikyhat_mark = ikyhat
!              iefwid_mark = iefwid
!              gmark = gamma_mark_i_1(ikyhat,iefwid)
            if (f_mark_i(ikyhat,iefwid) .lt. fmark) then  ! Reset bounds in ky, efwid and f_guess_mark
              fmark = f_mark_i(ikyhat,iefwid)
              ikymin = nkyhat
              ikymax = 1
              iefwmin = nefwid
              iefwmax = 1
              f_guess_mark = 1.0E20
              ikyhat_write = ikyhat
              iefwid_write = iefwid
              ikyhat_mark = ikyhat
              iefwid_mark = iefwid
            endif
!            f_guess_mark = f_guess(ikyhat,iefwid)
            ikymin = min(ikyhat,ikymin)
            ikymax = max(ikyhat,ikymax)
            iefwmin = min(iefwid,iefwmin)
            iefwmax = max(iefwid,iefwmax)
            f_guess_mark = min(f_guess(ikyhat,iefwid),f_guess_mark)
!            if ( (gamma_mark_i_1(ikyhat,iefwid) .gt. gmark) .or. (f_mark_i(ikyhat,iefwid) .lt. fmark) ) then
!              ikyhat_write = ikyhat
!              iefwid_write = iefwid
!            endif
          endif
        enddo ! iefwid
      enddo ! ikyhat
    endif  ! mode found

    if (l_write_out .and. id .eq. 0) then
      inquire(file=trim('out.scalefactor'//suffix),exist=iexist)
      if(iexist) then
        open(unit=88,file=trim('out.scalefactor'//suffix),status='old',position='append')
      else
        open(unit=88,file=trim('out.scalefactor'//suffix),status='new')
        write(88,*)"factor,(gamma(n),freq(n),flag,n=1,nmodes_in)"
        write(88,*) "flag key:  'K'  mode is kept"
        if (reject_tearing_flag .eq. 1)  write(88,*) "           'T' rejected for tearing parity"
        if (reject_i_pinch_flag .eq. 1)  write(88,*) "           'Pi' rejected for ion pinch"
        if (reject_e_pinch_flag .eq. 1)  write(88,*) "           'Pe' rejected for electron pinch"
        if (reject_th_pinch_flag .eq. 1) write(88,*) "           'Pth' rejected for thermal pinch"
        if (reject_EP_pinch_flag .eq. 1) write(88,*) "           'PEP' rejected for EP pinch"
!        if (reject_max_outer_panel_flag .eq. 1) write(88,*) "           'OP' rejected for ballooning-space max outside first panel"
        write(88,*) "           'QLR' rejected for QL ratio EP/|chi_i| < ", QL_ratio_thresh
        write(88,*) "           'TH2' rejected for <theta^2>  > ", theta_sq_thresh
        write(88,*) "           'F'   rejected for non-AE frequency > ", f_real(ir)*freq_AE_upper, " kHz"
        write(88,*) "           'BT'  mode growth rate is below threshold gamma_thresh = ", f_real(ir)*gamma_thresh, " kHz"
        write(88,*) "omega_TAE = ", f_real(ir)*omega_TAE(ir), " ;  omega_GAM = ", -f_real(ir)*omega_GAM(ir)
        if (l_real_units .eq. 1) write (88,*) "Frequencies in real units, plasma frame [kHz] ; (c_s/a)/(2*pi) = ", f_real(ir), " kHz"
      endif

      ky_write = kyhat(ikyhat_write)*tglf_zs_in(is_EP+1)/sqrt(tglf_mass_in(is_EP+1)*tglf_taus_in(is_EP+1))

      write(88,60) '--------------- ky*rho_EP= ', kyhat(ikyhat_write), '  (ky*rho_s=', ky_write, ')', &
                   '   width= ', efwid(iefwid_write), &
                   '   scalefactor= ', fmark, ' -------------------'
!      write (88,*) "np_1:", np, "np_2:", np_2, "np_3", np_3
!      write (88,*) "id_1:", id, "id_2:", id_2, "id_3", id_3
      do ifactor = 1,nfactor
        ikyhat = ikyhat_write
        iefwid = iefwid_write
        do n = 1,nmodes
          g(n) = f_real(ir)*growthrate_out(ikyhat,iefwid,ifactor,n)
          f(n) = f_real(ir)*frequency_out(ikyhat,iefwid,ifactor,n)
          if (lkeep_i_out(ikyhat,iefwid,ifactor,n)) then
            keep_label(n) = ' K  '
          else if (growthrate_out(ikyhat,iefwid,ifactor,n) .lt. gamma_thresh) then
            keep_label(n) = ' BT  '
          else if ((ltearing_i_out(ikyhat,iefwid,ifactor,n)).and.(reject_tearing_flag .eq. 1)) then
            keep_label(n) = ' T  '
          else if ((l_i_pinch_i_out(ikyhat,iefwid,ifactor,n)).and.(reject_i_pinch_flag .eq. 1)) then
            keep_label(n) = ' Pi '
          else if ((l_e_pinch_i_out(ikyhat,iefwid,ifactor,n)).and.(reject_e_pinch_flag .eq. 1)) then
            keep_label(n) = ' Pe '
          else if ((l_th_pinch_i_out(ikyhat,iefwid,ifactor,n)).and.(reject_th_pinch_flag .eq. 1)) then
            keep_label(n) = ' Pth'
          else if ((l_EP_pinch_i_out(ikyhat,iefwid,ifactor,n)).and.(reject_EP_pinch_flag .eq. 1)) then
            keep_label(n) = ' PEP'
!          else if ((l_max_outer_panel_i_out(ikyhat,iefwid,ifactor,n)).and.(reject_max_outer_panel_flag .eq. 1)) then
!            keep_label(n) = ' OP'
          else  if (l_QL_ratio_i_out(ikyhat,iefwid,ifactor,n)) then
            keep_label(n) = ' QLR'
          else  if (l_theta_sq_i_out(ikyhat,iefwid,ifactor,n)) then
            keep_label(n) = ' TH2'
          else if (frequency_out(ikyhat,iefwid,ifactor,n) .gt. freq_AE_upper) then 
            keep_label(n) = ' F  '
          else
            keep_label(n) = ' ?  '
          endif
        enddo
        if (l_real_units .eq. 0) then
          write(88,10)factor(ifactor),(g(n),f(n),keep_label(n),n=1,nmodes)
        else
          write(88,20)factor(ifactor),(g(n),f(n),keep_label(n),n=1,nmodes)
!          write(88,*) 'TGLF densities:', tglf_as_in(:)
!          write(88,*) 'TGLF a/L_n vals:', tglf_rlns_in(:)
        endif
      enddo  ! ifactor
      close(88)
    endif  ! l_write_out and id = 0 

!    if (imark_min .le. nfactor) then     ! mode found
    if (fmark .lt. 1.0E10) then     ! accepted mode with all constraints

      f1 = 2.0 * fmark
      f0 = 0.0
      delw = (w1-w0) / (nefwid-1)

      w1 = efwid(iefwmax) + delw
      w0 = efwid(iefwmin) - delw
      if (w1 .gt. width_max) w1 = width_max
      if (w0 .lt. width_min) w0 = width_min

      delky = (kyhat1-kyhat0) / nkyhat
      kyhat1 = kyhat(ikymax) + delky
      kyhat0 = kyhat(ikymin) - delky
      if (kyhat1 .gt. kyhat_max) kyhat1 = kyhat_max
      if (kyhat0 .lt. kyhat_min) kyhat0 = kyhat_min

!      if (k .eq. 1) then
!        f1 = fmark
!        f0 = 0.0
!      else ! (iterations after the first)
!        delf = (f1-f0) / 3.
!!        f1 = f_guess(ikyhat_mark,iefwid_mark) + delf
!!        f0 = f_guess(ikyhat_mark,iefwid_mark) - delf
!        f1 = fmark + delf
!        f0 = fmark - delf
!        if (f0 .lt. 0.0) f0 = 0.0
!        delw = (w1-w0) / 4.
!        w1 = efwid(iefwid_mark) + delw
!        w0 = efwid(iefwid_mark) - delw
!        if (w1 .gt. width_max) then
!          w0 = w0 - (w1-width_max)
!          w1 = width_max
!        else if (w0 .lt. width_min) then
!          w1 = w1 + (width_min-w0)
!          w0 = width_min
!        endif
!
!        delky = (kyhat1-kyhat0) / 4.
!        kyhat1 = kyhat(ikyhat_mark) + delky
!        kyhat0 = kyhat(ikyhat_mark) - delky
!        if (kyhat1 .gt. kyhat_max) then
!          kyhat0 = kyhat0 - (kyhat1-kyhat_max)
!          kyhat1 = kyhat_max
!        else if (kyhat0 .lt. kyhat_min) then
!          kyhat1 = kyhat1 + (kyhat_min-kyhat0)
!          kyhat0 = kyhat_min
!        endif
!
!      endif !  k > 1

    else  
      f0 = f1
      f1 = 10.*f1
    endif  ! Mode found

                   

  enddo ! k (recursive subdivision iterator)

  if(imark_min .gt. nfactor) then
    factor_in = 10000. !NaN, no threshold found
    width_in = efwid(1)
    kymark = kyhat(1)
  else
!    factor_in = factor(imark_min) !find the density threshold
!    factor_in = f_guess(ikyhat_mark,iefwid_mark)
    factor_in = f_guess_mark
    width_in = efwid(iefwid_mark)
    kymark = kyhat(ikyhat_mark)
  endif
!  print *, 'In scalefactor (3) factor_in=', factor_in, 'imark=', imark
!  print *, 'In scalefactor (4) factor=', factor

  ! if(write_out_flag .and. id .eq. 0) then
  !   print *, 'At',trim(suffix),' Scale Factor is ',factor_in,'with width = ',width_in,'ky = ',ky_in,'and freq = ',fmark
  ! endif

10 format(F11.4,4(2F14.7,A4))
15 format(E16.4,4(2ES14.5,A4))
20 format(F11.4,4(2ES14.5,A4))
25 format(I4,A30,I2,A10,I2,':',I2,':',I2)
30 format(I4,I4,I4,A30,I2,A10,I2,':',I2,':',I2)
40 format(I4,I4,I4,A30,I2,':',I2,':',I2)
50 format(A30,I2,':',I2,':',I2)
60 format(A27,F8.4,A12,F8.4,A1,A10,F8.4,A16,F10.4,A20)

end subroutine TGLFEP_kwscale_scan
