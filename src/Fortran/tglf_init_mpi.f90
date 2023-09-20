!\
! module for MPI version of TGLF
!/

      MODULE tglf_mpi

        use mpi

        IMPLICIT NONE

        ! local communicator
        integer :: iCommTglf    = -1
        integer :: iProcTglf    = -1
        integer :: nProcTglf    = -1
        integer :: iProc0Tglf   = 0
        integer :: iGroupIDTglf = 1

      END MODULE tglf_mpi
!
!______________________________________________________
!
!--------------------------------------------------------------
! tglf_init.f90
!
! PURPOSE:
!  Initialize external TGLF interface.
!---------------------------------------------------------------

subroutine tglf_init(path_in,mpi_comm_in)

  use tglf_interface
  use tglf_mpi

  implicit none

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in
  integer,          intent(in) :: mpi_comm_in
  integer :: ierr

  ! Local variable
  logical :: inputdat_flag

  tglf_path_in = path_in
  iCommTglf    = mpi_comm_in

  call MPI_COMM_RANK(iCommTglf,iProcTglf,ierr)
  call MPI_COMM_SIZE(iCommTglf,nProcTglf,ierr)

  inquire(file=trim(path_in)//'input.tglf.gen',exist=inputdat_flag)

  if (inputdat_flag .eqv. .true.) then
     if (tglf_quiet_flag_in .eqv. .false.) then
        print '(a,a,a)', 'INFO: (TGLF) tglf_init reading ',trim(path_in),'input.tglf.gen'
     endif
     call tglf_read_input
  endif

end subroutine tglf_init

