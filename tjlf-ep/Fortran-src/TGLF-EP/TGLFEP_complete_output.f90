!---------------------------------------------------
! TGLFEP_complete_output.f90
!
! PURPOSE: Write a complete output containing the
!  a value for all profile input radii. Those not
!  run or returning indeterminate values are 
!  omitted and replaced with inter/extrapolations.
!
! INPUT:  profile_in(scan_n)   ! Profile dimensions
! OUTPUT: profile_out(nr)      ! in TGLFEP_profile
!---------------------------------------------------

subroutine TGLFEP_complete_output(profile_in,profile_out,ir_min,ir_max,l_accept_profile)

  use TGLFEP_interface
  use TGLFEP_profile

  real, intent(in), dimension(scan_n) :: profile_in
  real, intent(out), dimension(nr) :: profile_out

  integer, intent(out) :: ir_min
  integer, intent(out) :: ir_max
  logical, intent(out), dimension(scan_n) :: l_accept_profile

  integer :: ir_1
  integer :: ir_2
  integer, dimension(nr) :: ir_mark

  integer :: i_r    ! Local iteration integer
  integer :: ir_counter
  integer :: ir_counter_max
  integer :: ir_in
  integer, dimension(scan_n) :: ir_out
  integer :: ir_max_0

  logical :: l_exclude
  real, dimension(nr) :: profile_mask

!! First populate profile_out with zeroes where profile_in is undefined.
!
!  do i_r = 1, irs-1     ! Done with do loop to avoid irs=1 error.
!    profile_out(i_r) = 0.0
!  enddo  ! i_r < irs

! Initialize all of profile_out to 0.

  profile_out(:) = 0.0

  l_accept_profile(:) = .true.
  do i_r = irs, irs+scan_n-1
    ir_in = i_r-irs+1
    if (input_profile_method .eq. 2) then
      ir_out(ir_in) = ir_exp(ir_in)
    else
      ir_out(ir_in) = i_r
    endif
    l_exclude = (profile_in(ir_in) .eq. profile_in(ir_in)+1.0)                ! Inf test
    l_exclude = (l_exclude .or. (profile_in(ir_in) .ne. profile_in(ir_in)) )  ! NaN test
!    l_exclude = ( l_exclude .or. (profile_in(ir_in) .gt. 9000.0) )        ! 10000. test (no result)
    l_exclude = (l_exclude .or. (profile_in(ir_in) .lt. 0.0) )            ! Negative value test
    if (.not. l_exclude) then
      profile_out(ir_out(ir_in)) = profile_in(ir_in)
    else
      profile_out(ir_out(ir_in)) = 0.0
      l_accept_profile(ir_in) = .false.
    endif
  enddo  ! irs <= i_r <= ir_max_0
    
  do i_r = ir_out(scan_n)+1, nr  ! Won't run if scan goes to max radius.
    profile_out(i_r) = 0.0
  enddo

  profile_mask(:) = 1.0
!  do i_r = 2, nr-1
!    if ( (profile_out(i_r-1) .eq. 0.0) .and. (profile_out(i_r+1) .eq. 0.0) ) then
!      profile_mask(i_r) = 0.0
!      ir_in = i_r-irs+1
!      if ((ir_in .ge. 1).and.(ir_in .le. scan_n)) l_accept_profile(i_r-irs+1) = .false.
!    endif
!  enddo
!  if (profile_out(2) .eq. 0.0) then
!    profile_mask(1) = 0.0
!    if (irs .eq. 1) l_accept_profile(1) = .false.
!  endif
!  if (profile_out(nr-1) .eq. 0.0) then
!    profile_mask(nr) = 0.0
!    if (irs+scan_n-1 .eq. nr) l_accept_profile(scan_n) = .false.
!  endif 
  do i_r = 2, scan_n-1
    if (.not.l_accept_profile(i_r-1) .and. .not.l_accept_profile(i_r+1)) then
      l_accept_profile(i_r) = .false.
      profile_mask(ir_out(i_r)) = 0.0
    endif
  enddo
  profile_out(:) = profile_mask(:)*profile_out(:)

!  print *, 'profile_out initial'
!  do i_r = 1, nr
!    print *, i_r, profile_out(i_r)
!  enddo
  
! Record indices for all non-zero entries in profile_out.

  ir_mark(:) = 0
  ir_counter = 0

  do i_r = 1, nr
    if (profile_out(i_r) .ne. 0.0) then
      ir_counter = ir_counter + 1
      ir_mark(ir_counter) = i_r
    endif
  enddo

  ir_counter_max = max(1,ir_counter)

  ir_min = ir_mark(1)
  ir_max = ir_mark(ir_counter_max)

!  print *, 'ir_counter_max', ir_counter_max
!  print *, 'ir_mark'
!  do i_r = 1, nr
!    print *, i_r, ir_mark(i_r)
!  enddo

! Repopulate zero entries with interpolation or input scale factor
! if blank space does not have a solution on both sides.

!! ----- Obsolete coding to smoothly fill in scale factor up to 10.0 -----
!!  if (ir_mark(1) .gt. 3) then
!!    profile_out(1:ir_mark(1)-3) =  10.0     !   factor_profile(1)
!!    ir_1 = ir_mark(1)-3
!!    ir_2 = ir_mark(1)
!!    profile_out(ir_1+1:ir_2-1) = profile_out(ir_1) + &
!!         (profile_out(ir_2)-profile_out(ir_1)) * &
!!         (rmin(ir_1+1:ir_2-1)-rmin(ir_1))/(rmin(ir_2)-rmin(ir_1))
!!  else
!!    ir_1 = ir_mark(1)
!!    profile_out(1:ir_mark(1)) = 10.0 -  & ! factor_profile(1) - &
!!!         (factor_profile(1)-profile_out(ir_1)) * &
!!         (10.0 - profile_out(ir_1)) * &
!!          (rmin(1:ir_1-1)-rmin(1))/(rmin(ir_1)-rmin(1))
!!  endif
!!-----------------------------------------------------------------------

  profile_out(1:ir_mark(1)) = profile_out(ir_mark(1))

!  ir_1 = ir_mark(1)
!  ir_2 = ir_mark(2)
!  profile_out(1:ir_2-1) = profile_out(ir_1) + &
!         (profile_out(ir_2)-profile_out(ir_1)) * &
!         (rmin(1:ir_2-1)-rmin(ir_1))/(rmin(ir_2)-rmin(ir_1))

!  print *, 'profile_out phase 1'
!  do i_r = 1, nr
!    print *, i_r, profile_out(i_r)
!  enddo

  ir_counter = 1

  do i_r = ir_mark(1)+1, ir_mark(ir_counter_max)
    ir_1 = ir_mark(ir_counter)
    ir_2 = ir_mark(ir_counter+1)
    if (profile_out(i_r) .eq. 0.0) then
      profile_out(i_r) = profile_out(ir_1) + &
         (profile_out(ir_2)-profile_out(ir_1)) * &
         (rmin(i_r)-rmin(ir_1))/(rmin(ir_2)-rmin(ir_1))
    else
       ir_counter = ir_counter + 1
    endif
  enddo

!  print *, 'profile_out phase 2'
!  do i_r = 1, nr
!    print *, i_r, profile_out(i_r)
!  enddo

!  ir_1 = ir_mark(ir_counter_max-1)
!  ir_2 = ir_mark(ir_counter_max)
!  profile_out(ir_1+1:nr) = profile_out(ir_1) + &
!         (profile_out(ir_2)-profile_out(ir_1)) * &
!         (rmin(ir_1:nr)-rmin(ir_1))/(rmin(ir_2)-rmin(ir_1))

!! ----- More obsolete coding to smoothly fill in scale factor up to 10.0 at edges. -----
!!  if (ir_mark(ir_counter_max)+3 .lt. nr) then
!!    profile_out(ir_mark(ir_counter_max)+3:nr) = 10.0   !  factor_profile(n_scan)
!!    ir_1 = ir_mark(ir_counter_max)
!!   ir_2 = ir_mark(ir_counter_max)+3
!!    profile_out(ir_1+1:ir_2-1) = profile_out(ir_1) + &
!!         (profile_out(ir_2)-profile_out(ir_1)) * &
!!         (rmin(ir_1+1:ir_2-1)-rmin(ir_1))/(rmin(ir_2)-rmin(ir_1))
!!  else if (ir_mark(ir_counter_max) .lt. nr) then
!!    ir_1 = ir_mark(ir_counter_max)
!!    profile_out(ir_1+1:nr) = profile_out(ir_1) - &
!!!         (factor_profile(1)-profile_out(ir_1)) * &`
!!         (10.0-profile_out(ir_1)) * &
!!         (rmin(1:ir_1-1)-rmin(1))/(rmin(ir_1)-rmin(1))
!!  endif
!! --------------------------------------------------------------------------------------

  profile_out(ir_mark(ir_counter_max):nr) = profile_out(ir_mark(ir_counter_max))

!  print *, 'profile_out phase 3'
!  do i_r = 1, nr
!    print *, i_r, profile_out(i_r)
!  enddo

end
