!------------------------------------------------------------
! TGLFEP_scalefactor.f90
!
! PURPOSE:
!  Calculate the EP density threshold 
!  usually called after TGLFEP_ky_widthscan
!------------------------------------------------------------

subroutine TGLFEP_scalefactor

  use mpi
  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface 

  implicit none
  integer :: id,np,ierr,STATUS(MPI_STATUS_SIZE)
  integer :: id_world, np_world
  integer :: i,n,k,imark
  integer :: k_max = 4
  logical :: iexist,write_out_flag = .false.
  real :: f0,f1,ft
  integer,parameter :: nfactor = 10
  real,dimension(nfactor) :: factor
  real,dimension(nfactor,nmodes) :: growthrate,growthrate_out &
                                   ,frequency,frequency_out
  real,dimension(nmodes) :: g,f
  real :: gmark,fmark

  character(7)  :: str_sf

  logical :: lkeep_i(nfactor,nmodes)
  logical :: lkeep_i_out(nfactor,nmodes)
  logical :: ltearing_i(nfactor,nmodes)
  logical :: ltearing_i_out(nfactor,nmodes)
  logical :: l_th_pinch_i(nfactor,nmodes)
  logical :: l_th_pinch_i_out(nfactor,nmodes)
  logical :: l_i_pinch_i(nfactor,nmodes)
  logical :: l_i_pinch_i_out(nfactor,nmodes)
  logical :: l_e_pinch_i(nfactor,nmodes)
  logical :: l_e_pinch_i_out(nfactor,nmodes)
  logical :: l_EP_pinch_i(nfactor,nmodes)
  logical :: l_EP_pinch_i_out(nfactor,nmodes)
  logical :: l_QL_ratio_i(nfactor,nmodes)
  logical :: l_QL_ratio_i_out(nfactor,nmodes)

  integer :: values(8)

  character(4), dimension(nmodes) :: keep_label

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)
  call MPI_COMM_RANK(TGLFEP_COMM_WORLD,id_world,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM_WORLD,np_world,ierr)

  growthrate       = 0.0
  growthrate_out   = 0.0
  frequency        = 0.0
  frequency_out    = 0.0
  lkeep_i(:,:)     = .true.
  lkeep_i_out(:,:) = .true.
  ltearing_i(:,:)     = .false.
  ltearing_i_out(:,:) = .false.
  l_th_pinch_i(:,:)     = .false.
  l_th_pinch_i_out(:,:) = .false.
  l_i_pinch_i(:,:)     = .false.
  l_i_pinch_i_out(:,:) = .false.
  l_e_pinch_i(:,:)     = .false.
  l_e_pinch_i_out(:,:) = .false.
  l_EP_pinch_i(:,:)     = .false.
  l_EP_pinch_i_out(:,:) = .false.
  l_QL_ratio_i(:,:)     = .false.
  l_QL_ratio_i_out(:,:) = .false.


  f0 = 0.0
  f1 = factor_in
!  print *, "factor_in in scalefactor (1)", factor_in

  tglf_write_wavefunction_flag_in = 0
  do k = 1,k_max

    do i = 1,nfactor
      factor(i) = (f1-f0)/nfactor*i+f0
!      print *, 'factor', i, factor(i)
    enddo

    do i = 1+id,nfactor,np
    l_wavefunction_out = 0
    factor_in = factor(i)
!      print *, "factor_in in scalefactor (2)", factor_in, 'at i=', i
    str_sf =  achar(mod(floor(factor_in/100.),10)+iachar("0"))  //  &
              achar(mod(floor(factor_in/10.),10)+iachar("0"))   //  &
              achar(mod(floor(factor_in),10)+iachar("0"))       //  &
              '.'                                               //  &
              achar(mod(floor(10*factor_in),10)+iachar("0"))    //  &
              achar(mod(floor(100*factor_in),10)+iachar("0"))   //  &
              achar(mod(floor(1000*factor_in),10)+iachar("0"))
    write(str_wf_file,'(A16,A8,A3,A7)') 'out.wavefunction', suffix, '_sf', str_sf

!    tglf_ef_file_in = trim(str_file)//suffix
!    tglf_ef_file_in = 'out.wavefunction.scratch'

!    tglf_dump_flag_in = .true.
!    file_dump_local = 'tglf_dump'//suffix

!    print *, 'In scalefactor (2.1) factor=', factor
!    call date_and_time(values=values)
!    if (id==0) print 25, id_3, 'before TGLFEP_ky at i=',i, 'time  ', values(5),values(6),values(7)

    if ((k .eq. 4) .and. (i .eq. nfactor) .and. (l_write_wavefunction .eq. 1)) l_wavefunction_out = 1
!    if ((k .eq. 1) .and. (i .eq. 1) .and. (l_write_wavefunction .eq. 1)) l_wavefunction_out = 1
    call TGLFEP_ky
    call date_and_time(values=values)
!   if (id==0) print 30, id,id_2,id_3, 'after TGLFEP_ky at i=',i, 'time  ', values(5),values(6),values(7)
!    print *, 'In scalefactor (2.25) factor=', factor
!    STOP
!    print *, 'proc_world: ', id_world, 'proc_loc: ', id, '  suffix:', suffix
      do n=1,nmodes
        growthrate(i,n) = get_growthrate(n)
        frequency(i,n)  = get_frequency(n)
        lkeep_i(i,n) = lkeep(n)
        ltearing_i(i,n) = ltearing(n)
        l_th_pinch_i(i,n) = l_th_pinch(n)
        l_i_pinch_i(i,n) = l_i_pinch(n)
        l_e_pinch_i(i,n) = l_e_pinch(n)
        l_EP_pinch_i(i,n) = l_EP_pinch(n)
        l_QL_ratio_i(i,n) = l_QL_ratio(n)
      enddo

    enddo  !  i (factor_in iterator)

    call MPI_BARRIER(TGLFEP_COMM,ierr)
    call date_and_time(values=values)
!    if (id_3==0) print 50, 'All after MPI_BARRIER at ', values(5),values(6),values(7)

    call MPI_ALLREDUCE(growthrate                      &
                      ,growthrate_out                  &
                      ,nfactor*nmodes                  &
                      ,MPI_DOUBLE_PRECISION            &
                      ,MPI_SUM                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)
    
    call MPI_ALLREDUCE(frequency                       &
                      ,frequency_out                   &
                      ,nfactor*nmodes                  &
                      ,MPI_DOUBLE_PRECISION            &
                      ,MPI_SUM                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)

    call MPI_ALLREDUCE(lkeep_i                         &
                      ,lkeep_i_out                     &
                      ,nfactor*nmodes                  &
                      ,MPI_LOGICAL                     &
                      ,MPI_LAND                        &
                      ,TGLFEP_COMM                     &
                      ,ierr)
    
    call MPI_ALLREDUCE(ltearing_i                      &
                      ,ltearing_i_out                  &
                      ,nfactor*nmodes                  &
                      ,MPI_LOGICAL                     &
                      ,MPI_LOR                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)

    call MPI_ALLREDUCE(l_th_pinch_i                    &
                      ,l_th_pinch_i_out                &
                      ,nfactor*nmodes                  &
                      ,MPI_LOGICAL                     &
                      ,MPI_LOR                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)

    call MPI_ALLREDUCE(l_i_pinch_i                     &
                      ,l_i_pinch_i_out                 &
                      ,nfactor*nmodes                  &
                      ,MPI_LOGICAL                     &
                      ,MPI_LOR                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)

    call MPI_ALLREDUCE(l_e_pinch_i                     &
                      ,l_e_pinch_i_out                 &
                      ,nfactor*nmodes                  &
                      ,MPI_LOGICAL                     &
                      ,MPI_LOR                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)

    call MPI_ALLREDUCE(l_EP_pinch_i                    &
                      ,l_EP_pinch_i_out                &
                      ,nfactor*nmodes                  &
                      ,MPI_LOGICAL                     &
                      ,MPI_LOR                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)

    call MPI_ALLREDUCE(l_QL_ratio_i                    &
                      ,l_QL_ratio_i_out                &
                      ,nfactor*nmodes                  &
                      ,MPI_LOGICAL                     &
                      ,MPI_LOR                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)

!    call date_and_time(values=values)
!    if (id==0) print 40, id, id_2, id_3, 'after MPI_ALLREDUCE at', values(5),values(6),values(7)

    write_out_flag = .true.
    if(write_out_flag .and. id .eq. 0) then
      inquire(file=trim('out.scalefactor'//suffix),exist=iexist)
      if(iexist) then
        open(unit=88,file=trim('out.scalefactor'//suffix),status='old',position='append')
      else
        open(unit=88,file=trim('out.scalefactor'//suffix),status='new')
        write(88,*) 'mode_flag ',mode_in,'ky ',ky_in,'width ',width_in
        write(88,*)"factor,(gamma(n),freq(n),flag,n=1,nmodes_in)"
        write(88,*) "flag key:  'K' mode is kept"
        if (reject_tearing_flag .eq. 1)  write(88,*) "           'T' rejected for tearing parity"
        if (reject_i_pinch_flag .eq. 1)  write(88,*) "           'Pi' rejected for ion pinch"
        if (reject_e_pinch_flag .eq. 1)  write(88,*) "           'Pe' rejected for electron pinch"
        if (reject_th_pinch_flag .eq. 1) write(88,*) "           'Pth' rejected for thermal pinch"
        if (reject_EP_pinch_flag .eq. 1) write(88,*) "           'PEP' rejected for EP pinch"
        write(88,*) "           'QLR' rejected for QL ratio EP/|chi_i| < ", QL_ratio_thresh
        write(88,'(A47,ES14.5,A4)') "             'F' rejected for non-AE frequency > ", f_real(ir)*freq_AE_upper, " kHz"
        if (l_real_units .eq. 1) write (88,*) "Frequencies in real units, plasma frame [kHz]"
      endif     

      write (88,*) "--------------------------------------------------"       
!      write (88,*) "np_1:", np, "np_2:", np_2, "np_3", np_3
!      write (88,*) "id_1:", id, "id_2:", id_2, "id_3", id_3
      do i = 1,nfactor
        do n = 1,nmodes
          g(n) = f_real(ir)*growthrate_out(i,n)
          f(n) = f_real(ir)*frequency_out(i,n)
          if (lkeep_i_out(i,n)) then
            keep_label(n) = ' K  '
          else if ((ltearing_i_out(i,n)).and.(reject_tearing_flag .eq. 1)) then 
            keep_label(n) = ' T  '
          else if ((l_i_pinch_i_out(i,n)).and.(reject_i_pinch_flag .eq. 1)) then
            keep_label(n) = ' Pi '
          else if ((l_e_pinch_i_out(i,n)).and.(reject_e_pinch_flag .eq. 1)) then
            keep_label(n) = ' Pe '
          else if ((l_th_pinch_i_out(i,n)).and.(reject_th_pinch_flag .eq. 1)) then
            keep_label(n) = ' Pth'
          else if ((l_EP_pinch_i_out(i,n)).and.(reject_EP_pinch_flag .eq. 1)) then
            keep_label(n) = ' PEP'
          else  if (l_QL_ratio_i_out(i,n)) then
            keep_label(n) = ' QLR'
          else
            keep_label(n) = ' F  '
          endif
        enddo
        if (l_real_units .eq. 0) then
          write(88,10)factor(i),(g(n),f(n),keep_label(n),n=1,nmodes)
        else
          write(88,20)factor(i),(g(n),f(n),keep_label(n),n=1,nmodes)
!          write(88,*) 'TGLF densities:', tglf_as_in(:)
!          write(88,*) 'TGLF a/L_n vals:', tglf_rlns_in(:)
        endif
      enddo

      close(88)
    endif

!    print *, 'In scalefactor (2.5) factor=', factor
    imark = 0
    do i = 1,nfactor
      do n = 1,nmodes
        gmark = growthrate_out(i,n)
        fmark = frequency_out(i,n)
        if ( lkeep_i_out(i,n) ) then
          f1 = factor(i)
          imark = i
          exit
        endif
      enddo
      if(imark .gt. 0) exit
    enddo

    if(imark .gt. 1) f0 = factor(imark-1)

  enddo !k

  if(imark .eq. 0) then
    factor_in = 10000. !NaN, no threshold found
  else
    factor_in = factor(imark) !find the density threshold
  endif
!  print *, 'In scalefactor (3) factor_in=', factor_in, 'imark=', imark
!  print *, 'In scalefactor (4) factor=', factor

  ! if(write_out_flag .and. id .eq. 0) then
  !   print *, 'At',trim(suffix),' Scale Factor is ',factor_in,'with width = ',width_in,'ky = ',ky_in,'and freq = ',fmark
  ! endif

10 format(F8.4,4(2F14.7,A4))
20 format(F8.4,4(2ES14.5,A4))
25 format(I4,A30,I2,A10,I2,':',I2,':',I2)
30 format(I4,I4,I4,A30,I2,A10,I2,':',I2,':',I2)
40 format(I4,I4,I4,A30,I2,':',I2,':',I2)
50 format(A30,I2,':',I2,':',I2)

end subroutine TGLFEP_scalefactor
