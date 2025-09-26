SUBROUTINE nm_sort2_dp(ARR1,ARR2)
  USE nm_type 
  USE nm_util, ONLY : check_length
  USE nm_nm, only: nm_get_index
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: ARR1,ARR2
  INTEGER(I4B) :: n
  INTEGER(I4B), DIMENSION(size(ARR1)) :: ind
  n=check_length(size(ARR1),size(ARR2),'nm_sort2')
  call nm_get_index(ARR1,ind)
  ARR1=ARR1(ind)
  ARR2=ARR2(ind)
END SUBROUTINE nm_sort2_dp

SUBROUTINE nm_sort2_sp(ARR1,ARR2)
  USE nm_type 
  USE nm_util, ONLY : check_length
  USE nm_nm, only: nm_get_index
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: ARR1,ARR2
  INTEGER(I4B) :: n
  INTEGER(I4B), DIMENSION(size(ARR1)) :: ind
  n=check_length(size(ARR1),size(ARR2),'nm_sort2')
  call nm_get_index(ARR1,ind)
  ARR1=ARR1(ind)
  ARR2=ARR2(ind)
END SUBROUTINE nm_sort2_sp
