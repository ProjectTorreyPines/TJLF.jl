subroutine TGLFEP_generate_input

  use TGLFEP_profile
  use TGLFEP_interface

  use mpi

  implicit none
  logical :: iexist
  integer :: i,j
  real, parameter :: pi=3.14159265359

  integer :: n_param

  integer :: i_proc
  integer :: n_proc
  integer :: ierr

  real, allocatable, dimension(:) :: as_min
  real, allocatable, dimension(:) :: as_max
  real, allocatable, dimension(:) :: taus_min
  real, allocatable, dimension(:) :: taus_max
  real, allocatable, dimension(:) :: rlts_min
  real, allocatable, dimension(:) :: rlts_max
  real, allocatable, dimension(:) :: rlns_min
  real, allocatable, dimension(:) :: rlns_max
  real, allocatable, dimension(:) :: ptot_scale
  real, allocatable, dimension(:) :: dptotdr_scale
  real :: dum_min
  real :: dum_max
  real :: rmin_min
  real :: rmin_max
  real :: rmaj_min
  real :: rmaj_max
  real :: q_min
  real :: q_max
  real :: shear_min
  real :: shear_max
  real :: shift_min
  real :: shift_max
  real :: kappa_min
  real :: kappa_max
  real :: s_kappa_min
  real :: s_kappa_max
  real :: delta_min
  real :: delta_max
  real :: s_delta_min
  real :: s_delta_max
  real :: zeta_min
  real :: zeta_max
  real :: s_zeta_min
  real :: s_zeta_max
  real :: betae_min
  real :: betae_max

  real :: x_dum
  integer :: n_dum
  integer :: seed_flag

  integer :: i_data_method
  integer, dimension(:), allocatable :: i_seed
  integer :: i_size

  real, allocatable, dimension(:) :: ran_arr
  integer :: i_dim

  n_param = 0

  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,ierr)

  inquire(file='input.TGLFEP.generate',exist=iexist)
  if(iexist) then
    if (i_proc .eq. 0) then
      open(unit=55,file='input.TGLFEP.generate',status='old')
      read(55,*) sign_bt
      read(55,*) sign_it
      read(55,*) nr
      read(55,*) seed_flag
      read(55,*) ns
    endif
    call MPI_BCAST(sign_bt,1,MPI_INTEGER,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(sign_it,1,MPI_INTEGER,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(nr,1,MPI_INTEGER,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(seed_flag,1,MPI_INTEGER,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(ns,1,MPI_INTEGER,0,TGLFEP_COMM_WORLD,ierr)

    geometry_flag = 1
    scan_n = nr
    irs = 1


    allocate(zs(ns))
    allocate(mass(ns))

    allocate(as_min(ns))
    allocate(as_max(ns))
    allocate(taus_min(ns))
    allocate(taus_max(ns))
    allocate(rlns_min(ns))
    allocate(rlns_max(ns))
    allocate(rlts_min(ns))
    allocate(rlts_max(ns))

    if (i_proc .eq. 0) then   ! File read from root process.
     do i = 1,ns
       read(55,*)
       read(55,*)
       read(55,*) x_dum
         zs(i) = x_dum
         print *, 'zs', i, 'at read', zs(i)
       read(55,*) x_dum
         mass(i) = x_dum

       if(i .ne. 1) then
         call readline(dum_min,dum_max,n_param)
         taus_min(i) = dum_min
         taus_max(i) = dum_max
         call readline(dum_min,dum_max,n_param)
         as_min(i) = dum_min
         as_max(i) = dum_max
       else
         as_min(i)   = 1.0  !electron density
         as_max(i)   = 1.0 
         taus_min(i) = 1.0  !electron temperature
         taus_max(i) = 1.0
       endif
 
       call readline(dum_min,dum_max,n_param)
       rlts_min(i) = dum_min
       rlts_max(i) = dum_max
       call readline(dum_min,dum_max,n_param)
       rlns_min(i) = dum_min
       rlns_max(i) = dum_max
     enddo
 
     read (55,*)
     read (55,*)

     call readline(rmin_min,rmin_max,n_param)
     call readline(rmaj_min,rmaj_max,n_param)
     call readline(q_min,q_max,n_param)
     call readline(shear_min,shear_max,n_param)
     call readline(shift_min,shift_max,n_param)
     call readline(kappa_min,kappa_max,n_param)
     call readline(s_kappa_min,s_kappa_max,n_param)
     call readline(delta_min,delta_max,n_param)
     call readline(s_delta_min,s_delta_max,n_param)
     call readline(zeta_min,zeta_max,n_param)
     call readline(s_zeta_min,s_zeta_max,n_param)
     call readline(betae_min,betae_max,n_param)

     close(55)

    endif ! Root process read.

    call MPI_BCAST(zs,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
!    print *, 'processor', i_proc, '  zs after bcast', zs
    call MPI_BCAST(mass,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(taus_min,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(taus_max,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(as_min,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(as_max,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rlts_min,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rlts_max,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rlns_min,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rlns_max,ns,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rmin_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rmin_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rmaj_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rmaj_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(q_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(q_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(shear_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(shear_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(shift_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(shift_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(kappa_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(kappa_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(s_kappa_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(s_kappa_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(delta_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(delta_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(s_delta_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(s_delta_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(zeta_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(zeta_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(s_zeta_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(s_zeta_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(betae_min,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(betae_max,1,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)

    if (nr .lt. 0) then
      nr = abs(nr)*n_param
      i_data_method = 1
    else
      i_data_method = 2
    endif

    allocate(as(nr,ns))
    allocate(taus(nr,ns))
    allocate(rlns(nr,ns))
    allocate(rlts(nr,ns))

    allocate(rmin(nr))
    allocate(rmaj(nr))
    allocate(q(nr))
    allocate(shear(nr))

    allocate(q_prime(nr))
    allocate(p_prime(nr))
    allocate(shift(nr))
    allocate(kappa(nr))
    allocate(s_kappa(nr))
    allocate(delta(nr))
    allocate(s_delta(nr))
    allocate(zeta(nr))
    allocate(s_zeta(nr))

    allocate(zeff(nr))
    allocate(betae(nr))

    allocate(rho_star(nr))
    allocate(omega_TAE(nr))
    allocate(omega_GAM(nr))
    allocate(gamma_E(nr))
    allocate(gamma_p(nr))

    allocate(ptot_scale(nr))
    allocate(dptotdr_scale(nr))

    i_dim = 12 + 4*ns
    allocate(ran_arr(i_dim))

! Generate the random "profile" on the root i_proc=0
    if (i_proc .eq. 0) then
! Initialize the random number sequence.
     call random_seed(size=i_size)
     allocate(i_seed(i_size))
     if (seed_flag .eq. 1) then
       call random_seed(put=i_seed)
     else
       call random_seed(get=i_seed)
!       print *, i_seed
     endif
!
     select case(i_data_method)
       case(1)
         print *, 'Grid data not yet operational.'
       case(2)
         do j = 1, nr
           call random_number(ran_arr)
           do i = 1, ns
             as(j,i) = as_min(i) + ran_arr(i)*(as_max(i)-as_min(i))
!             print *, j, i, 'as=', as(j,i)
             taus(j,i) = taus_min(i) + ran_arr(i+1)*(taus_max(i)-taus_min(i))
!             print *, j, i, 'taus=', taus(j,i), taus_min(i), taus_max(i), ran_arr(i+1)
             rlns(j,i) = rlns_min(i) + ran_arr(i+2)*(rlns_max(i)-rlns_min(i))
             rlts(j,i) = rlts_min(i) + ran_arr(i+3)*(rlts_max(i)-rlts_min(i))
           enddo ! Species index
           rmin(j) = rmin_min + ran_arr(4*ns+1)*(rmin_max-rmin_min)
           rmaj(j) = rmaj_min + ran_arr(4*ns+2)*(rmaj_max-rmaj_min)
           q(j) = q_min + ran_arr(4*ns+3)*(q_max-q_min)
           shear(j) = shear_min + ran_arr(4*ns+4)*(shear_max-shear_min)
           shift(j) = shift_min + ran_arr(4*ns+5)*(shift_max-shift_min)
           kappa(j) = kappa_min + ran_arr(4*ns+6)*(kappa_max-kappa_min)
           s_kappa(j) = s_kappa_min + ran_arr(4*ns+7)*(s_kappa_max-s_kappa_min)
           delta(j) = delta_min + ran_arr(4*ns+8)*(delta_max-delta_min)
           s_delta(j) = s_delta_min + ran_arr(4*ns+9)*(s_delta_max-s_delta_min)
           zeta(j) = zeta_min + ran_arr(4*ns+10)*(zeta_max-zeta_min)
           s_zeta(j) = s_zeta_min + ran_arr(4*ns+11)*(s_zeta_max-s_zeta_min)
           betae(j) = betae_min + ran_arr(4*ns+12)*(betae_max-betae_min)
           omega_TAE(j) = (2./betae(j))**0.5/2./q(j)/rmaj(j)
           omega_GAM(j) = (1./rmaj(j))*(1+taus(j,2))**0.5/(1.+1./(2.*q(j)))
           ptot_scale(j) = 0.
           dptotdr_scale(j) = 0.
           zeff(j) = 0.
           do i = 2, ns
             ptot_scale(j) = ptot_scale(j) + as(j,i)*taus(j,i)
             dptotdr_scale(j) = dptotdr_scale(j) + (rlns(j,i)+rlts(j,i))*as(j,i)*taus(j,i)  
             zeff(j) = zeff(j) + as(j,i)*zs(i)
           enddo ! Species sum
         enddo  ! Total number of data points
     end select
    endif ! i_proc=0
    
    rho_star(:) = 0.001

    call MPI_BCAST(taus,ns*nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(as,ns*nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rlts,ns*nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rlns,ns*nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rmin,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(rmaj,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(q,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(shear,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(shift,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(kappa,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(s_kappa,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(delta,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(s_delta,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(zeta,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(s_zeta,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
    call MPI_BCAST(betae,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)

    call MPI_BCAST(zeff,nr,MPI_DOUBLE_PRECISION,0,TGLFEP_COMM_WORLD,ierr)
!    if ((i_proc .eq. 0) .or. (i_proc .eq. 1)) print *, 'processor ', i_proc, ' :  zeff = ', zeff

    gamma_E(:) = 1.0E-7
    gamma_p(:) = 1.0E-7
    q_prime(:) = (q(:)/rmin(:))**2 * shear(:)
    p_prime(:) = -1.0*(abs(q(:))/rmin(:))*betae(:)/(8*pi)*dptotdr_scale(:)   ! (Ptot/Pe)*(1/Ptot) cancel, with Pe=1

  else

    print *, 'input.TGLFEP.generate file not found'
    stop

  endif

end subroutine TGLFEP_generate_input

subroutine readline(var_min,var_max,n_count)

  implicit none
  real, intent(out) :: var_min
  real, intent(out) :: var_max
  integer, intent(inout) :: n_count

  character(50) :: line
  character(10) :: entry_1
  character(10) :: entry_2
  character(1) :: lead_char

  integer :: comma_index
  integer :: close_index

  read (55,'(A50)') line
!  print *, line
!  read (line,'(A1)') lead_char
  lead_char = line(1:1)
  if (lead_char .eq. '(') then
!    print *, line
!    read (line,10) var_min, var_max
    comma_index = scan(line, ",", .true.)
    close_index = scan(line, ")", .true.)
    entry_1 = line(2:comma_index-1)
    entry_2 = line(comma_index+1:close_index)
    read (entry_1,*) var_min
    read (entry_2,*) var_max
    n_count = n_count + 1
  else
    read (line,*) var_min
    var_max = var_min
  endif
!  print *, line
!  print *, var_min, var_max

  return

!  10 format ('(',G,',',G)
  10 format('(',G,G)

end subroutine readline

