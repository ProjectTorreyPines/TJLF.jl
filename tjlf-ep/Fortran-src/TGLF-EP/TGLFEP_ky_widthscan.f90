!------------------------------------------------------------
! TGLFEP_ky_widthscan.f90
!
! PURPOSE:
!  Calculate the growth rate and frequency vs width
!  Find the suitable width for max gamma_AE
!------------------------------------------------------------

subroutine TGLFEP_ky_widthscan

  use mpi
  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface
  use tglf_dimensions
  use tglf_weight
  use tglf_pkg

  implicit none
  integer :: id,np,ierr,STATUS(MPI_STATUS_SIZE)
  logical :: iexist, write_out_flag = .true.
  character(30) :: str_file
  character(4)  :: str_width
  integer :: i,n,k,noff,i_plot
  integer :: nwidth
  integer :: i_spec
  integer :: i_flux
  integer :: i_field
  real,allocatable,dimension(:) :: width
  real,allocatable,dimension(:,:) :: growthrate,growthrate_out &
                                    ,frequency,frequency_out
  real, allocatable, dimension(:,:,:,:,:) :: QL_weights
  real, allocatable, dimension(:,:,:,:,:) :: QL_weights_out
  real :: delta_width
  real :: w,g(nmodes),f(nmodes)
  logical :: find_max
  real :: gmark,fmark

  logical, allocatable, dimension(:,:) :: lkeep_i, lkeep_i_out

!  real :: wave(maxmodes*4)
!  character(len=81) :: header
!  character(len=11) :: theta="    theta  "
!  character(len=22) :: phi="  RE(phi)    IM(phi)  "
!  character(len=24) :: Bper="  RE(Bper)    IM(Bper)  "


  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)

  delta_width = 0.01
!  delta_width = 0.02
!  delta_width = 0.03
!  delta_width = 0.05

  nwidth = nint((width_max - width_min)/delta_width)

  allocate(width(nwidth))
  allocate(growthrate(nwidth,nmodes))
  allocate(growthrate_out(nwidth,nmodes))
  allocate(frequency(nwidth,nmodes))
  allocate(frequency_out(nwidth,nmodes))
  allocate(lkeep_i(nwidth,nmodes))
  allocate(lkeep_i_out(nwidth,nmodes))
  allocate(QL_weights(3,2,3,nwidth,nmodes))     ! 3=(e-,ion1,ion2),2=(energy,particle),3=(phi,B_perp,B_par)
  allocate(QL_weights_out(3,2,3,nwidth,nmodes)) !
    
  growthrate     = 0.0
  growthrate_out = 0.0
  frequency      = 0.0
  frequency_out  = 0.0
  QL_weights = 0.0
  QL_weights_out = 0.0
  lkeep_i(:,:)     = .true.
  lkeep_i_out(:,:) = .true.

  do i = 1,nwidth
    width(i) = width_min + delta_width*(i-1)
  enddo

  do i = 1+id,nwidth,np

    width_in = width(i)

    str_width = achar(floor(width_in)+iachar("0"))//'.'//achar(mod(floor(10*width_in),10)+iachar("0"))&
                //achar(mod(floor(100*width_in),10)+iachar("0"))
    write(str_file,'(A18,A4)')'out.wavefunction_w',str_width

    call TGLFEP_ky

    call write_wavefunction_out(trim(str_file)//suffix)
    print *, 'Wavefunction write called by ', id, ' for ', suffix, '.'


!!!! Mostly copied from tglf_inout.f90 subroutine write_wavefunction_out.
!     open(unit=33,file=trim(str_file)//suffix,status='replace')
!     header = theta//phi//Bper
!     !
!     write(33,*) nmodes_out,nfields_out,max_plot_out
!     write(33,*) header
!     do i_plot = 1,max_plot_out
!        do n=1,nmodes_out
!           noff=2*nfields_out*(n-1)
!           wave(noff+1) = REAL(wavefunction_out(n,1,i_plot))
!           wave(noff+2) = AIMAG(wavefunction_out(n,1,i_plot))
!           wave(noff+3) = REAL(wavefunction_out(n,2,i_plot))
!           wave(noff+4) = AIMAG(wavefunction_out(n,2,i_plot))
!        enddo
!        write(33,*) angle_out(i_plot),(wave(k),k=1,nmodes_out*nfields_out*2)
!     enddo
!     close(33)
!!!!!!!!!! End copied code.
    do n = 1,nmodes
      growthrate(i, n) = get_growthrate(n)
      frequency(i, n)  = get_frequency(n)
      lkeep_i(i,n) = lkeep(n)
      do i_field = 1, 3
        QL_weights(1,2,i_field,i,n) = get_QL_particle_flux(n,1,i_field)
        QL_weights(1,1,i_field,i,n) = get_QL_energy_flux(n,1,i_field)
        QL_weights(2,2,i_field,i,n) = get_QL_particle_flux(n,2,i_field)
        QL_weights(2,1,i_field,i,n) = get_QL_energy_flux(n,2,i_field)
        QL_weights(3,2,i_field,i,n) = get_QL_particle_flux(n,3,i_field)
        QL_weights(3,1,i_field,i,n) = get_QL_energy_flux(n,3,i_field)
      enddo
    enddo

  enddo

  call MPI_BARRIER(TGLFEP_COMM,ierr)
  call MPI_ALLREDUCE(growthrate                      &
                    ,growthrate_out                  &
                    ,nwidth*nmodes                   &
                    ,MPI_DOUBLE_PRECISION            &
                    ,MPI_SUM                         &
                    ,TGLFEP_COMM                     &
                    ,ierr)
  
  call MPI_ALLREDUCE(frequency                       &
                    ,frequency_out                   &
                    ,nwidth*nmodes                   &
                    ,MPI_DOUBLE_PRECISION            &
                    ,MPI_SUM                         &
                    ,TGLFEP_COMM                     &
                    ,ierr)

  call MPI_ALLREDUCE(lkeep_i                         &
                    ,lkeep_i_out                     &
                    ,nwidth*nmodes                   &
                    ,MPI_LOGICAL                     &
                    ,MPI_LAND                        &
                    ,TGLFEP_COMM                     &
                    ,ierr)

  call MPI_ALLREDUCE(QL_weights                      &
                    ,QL_weights_out                  &
                    ,nwidth*nmodes*18                &
                    ,MPI_DOUBLE_PRECISION            &
                    ,MPI_SUM                         &
                    ,TGLFEP_COMM                     &
                    ,ierr)

  if(write_out_flag .and. id .eq. 0) then
    write(str_file,'(A18,I1)')'out.ky_widthscan_m',mode_in
    open(unit=9,file=trim(str_file)//suffix,status='replace')
    ! inquire(file=trim(str_file//suffix),exist=iexist)
    ! if(iexist) then
    !   open(unit=33,file=trim(str_file//suffix),status='old',position='append')
    ! else
    !   open(unit=33,file=trim(str_file//suffix),status='new')
    ! endif

    write(9,*)"widthscan at ky =",ky_in,'mode_flag ',mode_in,'factor ',factor_in
    write(9,*)"width,(gamma(n),freq(n),n=1,nmodes_in)"
    if (l_real_units .eq. 1) write (9,*) "Frequency in real units, plasma frame [kHz]"
    do i=1,nwidth
      w = width(i)
      do n=1,nmodes
        g(n) = f_real(ir)*growthrate_out(i,n)
        f(n) = f_real(ir)*frequency_out(i,n)
      enddo
      if (l_real_units .eq. 0) then
        write(9,10)w,(g(n),f(n),n=1,nmodes)
      else
        write(9,20)w,(g(n),f(n),n=1,nmodes)
      endif
    enddo
    close(9)

    write(str_file,'(A21,I1)')'out.wscan_ql_weight_m',mode_in
    open(unit=9,file=trim(str_file)//suffix,status='replace')
    write(9,*)"widthscan at ky =",ky_in,'mode_flag ',mode_in,'factor ',factor_in
    write(9,*)"width, [fluxes (elec_E(1:3), ion_E(1:3), EP_E(1:3), elec_n(1:3), ion_n(1:3), EP_n(1:3)),i=1,nmodes_in]"
    do i=1,nwidth
      w = width(i)
      write(9,30) w, ((((QL_weights_out(i_spec,i_flux,i_field,i,n),i_field=1,3),i_spec=1,3),i_flux=1,2),n=1,nmodes)
    enddo
    close(9)
  endif

  find_max = .true. 
  !find_max = .false.
  if(find_max) then
    width_in = width_min !- 0.01
    gmark = 0.0
    fmark = 0.0
    do i = 1,nwidth
      do n = 1,nmodes
        if(lkeep_i_out(i,n) .and. growthrate_out(i,n) .gt. gmark) then
          gmark = growthrate_out(i,n)
          fmark = frequency_out(i,n)
          width_in = width(i)
        endif
      enddo
    enddo
    !if(id .eq. 0) print *,ir,'Find the width = ',width_in,'for mode_flag ',mode_in,', factor ',factor_in,', ky ',ky_in
  endif

  deallocate(width)
  deallocate(growthrate)
  deallocate(growthrate_out)
  deallocate(frequency)
  deallocate(frequency_out)
  deallocate(lkeep_i)
  deallocate(lkeep_i_out)

10 format(F5.2,8F12.7)
20 format(F5.2,8ES12.5)
30 format(F5.2,72F12.7)

end subroutine TGLFEP_ky_widthscan
