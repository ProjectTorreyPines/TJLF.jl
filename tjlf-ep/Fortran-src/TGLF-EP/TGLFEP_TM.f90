!------------------------------------------------------------
! TGLFEP_TM.f90
!
! PURPOSE:
!  Calculate the growth rate and frequency spectra
!  Print out the fluxes if needed
!------------------------------------------------------------
module TGLFEP_ky_spectrum

  integer,parameter :: nky=30
  real,parameter :: ky=0.15

end module TGLFEP_ky_spectrum

subroutine TGLFEP_TM

  use mpi
  use tglf_interface
  use TGLFEP_interface
  use TGLFEP_ky_spectrum

  implicit none
  integer :: id,np,ierr,STATUS(MPI_STATUS_SIZE)
  logical :: write_flux_flag = .false.

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)

  call tglf_init('',TGLFEP_COMM)

  call TGLFEP_tglf_map

  tglf_use_transport_model_in = .true.

  tglf_kygrid_model_in = 0
  tglf_ky_in           = ky
  tglf_nky_in          = nky

  tglf_nmodes_in = nmodes

  tglf_nbasis_min_in = 32
  tglf_nbasis_max_in = 32
  tglf_nxgrid_in     = 32
  
  tglf_width_in      = width_in
  tglf_find_width_in = .false.
  
  call tglf_run_mpi

  if(id .eq. 0) then

    call write_eigenvalue_spectrum

    if(write_flux_flag) then
      print *,ir&
             ,tglf_elec_pflux_out,tglf_ion_pflux_out(1),tglf_ion_pflux_out(2)&
             ,tglf_elec_eflux_out,tglf_ion_eflux_out(1),tglf_ion_eflux_out(2)
      call write_flux_spectrum
    endif

    ! call write_potential_spectrum
    
  endif

end subroutine TGLFEP_TM

subroutine write_eigenvalue_spectrum

  use tglf_pkg
  use TGLFEP_interface
  use TGLFEP_ky_spectrum

  implicit none
  character(17) :: str_file
  integer :: i,n
  real :: ky_out,growthrate_out(nmodes),frequency_out(nmodes)
   
  write(str_file,'(A16,I1)')'out.eigenvalue_m',mode_in
  open(unit=15,file=trim(str_file//suffix),status='replace')

  write(15,*)"gyro-bohm normalized eigenvalue spectra for mode_flag ",&
              mode_in,"factor ",factor_in,"width ",width_in
  write(15,*)"ky,(gamma(n),freq(n),n=1,nmodes_in)"
  do i=1,nky
    ky_out = get_ky_spectrum_out(i)
    do n=1,nmodes
      growthrate_out(n) = get_eigenvalue_spectrum_out(1,i,n)
      frequency_out(n) = get_eigenvalue_spectrum_out(2,i,n)
    enddo
    write(15,10)ky_out,(growthrate_out(n),frequency_out(n),n=1,nmodes)
  enddo

  close(15)

10 format(F8.4,8F12.7)

end subroutine write_eigenvalue_spectrum

subroutine write_flux_spectrum

  use tglf_pkg
  use TGLFEP_interface
  use TGLFEP_ky_spectrum

  implicit none
  character(11) :: str_file
  integer :: i,j,is,imax,jmax
  real :: dky
  real :: dky0,dky1,ky0,ky1
  real :: pflux0,eflux0,pflux1,eflux1
  real :: pflux_out,eflux_out

  write(str_file,'(A10,I1)')'out.flux_m',mode_in
  open(unit=16,file=trim(str_file//suffix),status='replace')

  do is=1,3
    do j=1,2
      write(16,*)"species = ",is,"field =",j
      write(16,*)"ky,particle flux,energy flux"

      pflux0 = 0.0
      eflux0 = 0.0

      dky0=0.0
      ky0=0.0
      dky = get_ky_spectrum_out(1)
      do i=1,nky
        ky1 = get_ky_spectrum_out(i)
        if(i.eq.1)then
          dky1=ky1
        else
          dky = LOG(ky1/ky0)/(ky1-ky0)
          dky1 = ky1*(1.0 - ky0*dky)
          dky0 = ky0*(ky1*dky - 1.0)
        endif

        pflux1 = 0.0
        eflux1 = 0.0
        do imax = 1,nmodes
          pflux1 = pflux1 + get_flux_spectrum_out(1,is,j,i,imax)
          eflux1 = eflux1 + get_flux_spectrum_out(2,is,j,i,imax)
        enddo
        pflux_out = dky0*pflux0 + dky1*pflux1
        eflux_out = dky0*eflux0 + dky1*eflux1
        write(16,10)ky1,pflux_out,eflux_out
        pflux0 = pflux1
        eflux0 = eflux1
        ky0 = ky1
      enddo  ! i
    enddo  ! j
  enddo  ! is

  close(16)

10 format(F8.4,2F12.7)

end subroutine write_flux_spectrum

subroutine write_potential_spectrum

  use tglf_pkg
  use TGLFEP_interface
  use TGLFEP_ky_spectrum

  implicit none
  character(16) :: str_file
  integer :: i,n
  real :: phi

  write(str_file,'(A15,I1)')'out.potential_m',mode_in
  open(unit=17,file=trim(str_file//suffix),status='replace')

  write(17,*)"gyro-bohm normalized potential fluctuation amplitude spectra"
  write(17,*)"ky,potential"
  do i=1,nky
    phi = 0.0
    do n = 1,nmodes
      phi = phi + get_field_spectrum_out(2,i,n)
    enddo
    phi = sqrt(phi)
    write(17,10)get_ky_spectrum_out(i),phi
  enddo

  close(17)

 10 format(F8.4,F12.7)

end subroutine write_potential_spectrum
