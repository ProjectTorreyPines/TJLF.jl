!---------------------------------------------------
! TGLFEP_interpolate.f90
!
! PURPOSE: Takes an array of length nr and repopulates
!   all values that are zero with interpolations or
!   extrapolations from adjacent points.
!
! INPUT:  profile_in(nr)       ! Profile dimensions
! OUTPUT: profile_out(nr)      ! in TGLFEP_profile
!---------------------------------------------------

subroutine TGLFEP_interpolate(profile_in,profile_out)

  use TGLFEP_interface
  use TGLFEP_profile

!  use EXPRO_interface
  use expro

  real, intent(in), dimension(nr) :: profile_in
  real, intent(out), dimension(nr) :: profile_out

  real, dimension(:), allocatable :: r_non0
  real, dimension(nr) :: r_list
  integer, dimension(:), allocatable :: ir_non0
  real, dimension(:), allocatable :: profile_in_non0

  integer :: i_r    ! Local iteration integer
  integer :: ir_counter
  integer :: ir_counter_max

! First, count the nonzero entries in profile_in

  ir_counter_max = 0
  do i_r = 1, nr
    if (profile_in(i_r) .ne. 0.0) ir_counter_max = ir_counter_max + 1
    r_list(i_r) = EXPRO_rmin(i_r)
  enddo

  allocate(profile_in_non0(ir_counter_max))
  allocate(r_non0(ir_counter_max))
  allocate(ir_non0(ir_counter_max))

! Now construct compressed arrays with zeros removed

  ir_counter = 0
  do i_r = 1, nr
    if (profile_in(i_r) .ne. 0.0) then
      ir_counter = ir_counter + 1
      profile_in_non0(ir_counter) = profile_in(i_r)
      r_non0(ir_counter) = EXPRO_rmin(i_r)
      ir_non0(ir_counter) = i_r
    endif
  enddo

! Extend first non-zero point into zero values at center.  
  do i_r = 1, ir_non0(1)
    profile_out(i_r) = profile_in(ir_non0(1))
  enddo

! Extend last non-zero point into zero values at edge.
  do i_r = ir_non0(ir_counter_max), nr
    profile_out(i_r) = profile_in(ir_non0(ir_counter_max))
  enddo

! Straight linear interpolation for remaining points.
  ir_counter = ir_non0(1)
  do i_r = ir_non0(1)+1, ir_non0(ir_counter_max)-1
    if (profile_in(i_r) .eq. 0.0) then
      profile_out(i_r) =    profile_in(ir_counter) +  (r_list(i_r) - r_non0(ir_counter)) *         &
                          (profile_in(ir_non0(ir_counter+1))-profile_in(ir_non0(ir_counter)) ) / &
                                   ( r_non0(ir_counter+1) - r_non0(ir_counter) )
    print *, 'Zero overwrite: ', i_r, r_list(i_r), ir_counter, ir_non0(ir_counter), r_non0(ir_counter), profile_in
    else
      profile_out(i_r) = profile_in(i_r)
      ir_counter = ir_counter + 1
    endif
  enddo

! ******  Cubic spline unsuitable, spurious min/max   ******
! Now just call the interpolation routine from $GACODE_ROOT/shared/math interpolating routine.
!  call cub_spline(r_non0,profile_in_non0,ir_counter_max,r_list,profile_out,nr)

  deallocate(profile_in_non0)
  deallocate(r_non0)
  deallocate(ir_non0)

end
