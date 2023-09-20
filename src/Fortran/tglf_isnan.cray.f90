  logical function tglf_isnan(x)
  ! 
  ! test real variables for NAN 
  !
  implicit none
  real, intent(in) :: x

  tglf_isnan = ISNAN(x)

  end function tglf_isnan
   
      
