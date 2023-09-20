  logical function tglf_isnan(x)
  ! 
  ! test function for platforms where isnan is not intrinsic.
  ! This numerical test may not work depending on the compiler
  !
  implicit none
  real, intent(in) :: x

  tglf_isnan = .false.

  if (x/=x) tglf_isnan = .true.

  end function tglf_isnan
   
      
