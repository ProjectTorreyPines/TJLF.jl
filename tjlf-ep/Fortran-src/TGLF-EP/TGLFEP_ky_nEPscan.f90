!------------------------------------------------------------
! TGLFEP_ky_nEPscan.f90
!
! PURPOSE:
!  Calculate the growth rate and frequency vs nEP
!------------------------------------------------------------

subroutine TGLFEP_ky_nEPscan

  use mpi
  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface

  implicit none
  integer :: id,np,ierr,STATUS(MPI_STATUS_SIZE)

  integer,parameter :: nmode_flag = 4
  integer,parameter :: nfactor = 14
  real :: factor(nfactor)

  integer :: i,ifactor,n
  real,dimension(nmode_flag,nfactor,nmodes) :: growthrate,growthrate_out &
                                              ,frequency,frequency_out
  real :: g(nmodes),f(nmodes)

  do i = 1,nfactor
    factor(i) = real(i)/10.0
  enddo

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)

  growthrate     = 0.0
  growthrate_out = 0.0
  frequency      = 0.0
  frequency_out  = 0.0

  do i = 1+id,nmode_flag*nfactor,np
    mode_in = (i-1)/nfactor+1
    ifactor = i-(mode_in-1)*nfactor
    factor_in = factor(ifactor)

    call TGLFEP_ky
    
    do n = 1,nmodes
      growthrate(mode_in, ifactor, n) = get_growthrate(n)
      frequency(mode_in, ifactor, n)  = get_frequency(n)
    enddo

  enddo

  call MPI_BARRIER(TGLFEP_COMM,ierr)

  call MPI_ALLREDUCE(growthrate                      &
                    ,growthrate_out                  &
                    ,nmode_flag*nfactor*nmodes       &
                    ,MPI_DOUBLE_PRECISION            &
                    ,MPI_SUM                         &
                    ,TGLFEP_COMM                     &
                    ,ierr)
  
  call MPI_ALLREDUCE(frequency                       &
                    ,frequency_out                   &
                    ,nmode_flag*nfactor*nmodes       &
                    ,MPI_DOUBLE_PRECISION            &
                    ,MPI_SUM                         &
                    ,TGLFEP_COMM                     &
                    ,ierr)

  if(id .eq. 0) then
    open(unit=20,file=trim('out.ky_nEPscan'//suffix),status='replace')

    do mode_in = 1,nmode_flag
      write(20,*) 'mode_flag ',mode_in,'ky ',ky_in,'width ',width_in
      write(20,*)"factor,(gamma(n),freq(n),n=1,nmodes_in)"
      do ifactor = 1,nfactor
        do n = 1,nmodes
          g(n) = growthrate_out(mode_in,ifactor,n)
          f(n) = frequency_out(mode_in,ifactor,n)
        enddo
        write(20,10)factor(ifactor),(g(n),f(n),n=1,nmodes)
      enddo
    enddo

    close(20)
  endif

10 format(F5.2,8F12.7)

end subroutine TGLFEP_ky_nEPscan
