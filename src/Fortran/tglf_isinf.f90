
  logical function tglf_isinf(x)
  ! 
  ! test real variables for NAN 
  !
  implicit none
  real, intent(in) :: x

  tglf_isinf = .false.

  if (ABS(x) > HUGE(x)) tglf_isinf = .true.

  end function tglf_isinf
   
      
