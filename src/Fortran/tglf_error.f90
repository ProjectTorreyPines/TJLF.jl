  SUBROUTINE tglf_error(itype,error_message)
  
  ! error handeling routine for TGLF
  ! simply writes the error or informatin message
  ! if itype = 1 the error was fatal so TGLF is shutdown
  
  implicit none
  integer, intent(in):: itype
  character(len=*), intent(in) :: error_message
 
  print *,'ERROR: (tglf) ',error_message

  if (itype == 1) then
     CALL tglf_shutdown
     STOP
  endif
 
  END SUBROUTINE tglf_error   
      
