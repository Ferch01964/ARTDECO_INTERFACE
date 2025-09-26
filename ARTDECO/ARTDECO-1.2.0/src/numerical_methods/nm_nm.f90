MODULE nm_nm
  INTERFACE
     SUBROUTINE nm_gaussNodes(x1,x2,x,A)
       USE nm_type
       REAL(DP), INTENT(IN) :: x1,x2
       REAL(DP), DIMENSION(:), INTENT(OUT) :: x,A
     END SUBROUTINE nm_gaussNodes
  END INTERFACE

  INTERFACE nm_get_index
     SUBROUTINE get_index_dp(arr,index)
       USE nm_type
       REAL(DP), DIMENSION(:), INTENT(IN) :: arr
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
     END SUBROUTINE get_index_dp
     SUBROUTINE get_index_sp(arr,index)
       USE nm_type
       REAL(SP), DIMENSION(:), INTENT(IN) :: arr
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
     END SUBROUTINE get_index_sp
     SUBROUTINE get_index_i4b(iarr,index)
       USE nm_type
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
     END SUBROUTINE get_index_i4b
  END INTERFACE nm_get_index

   INTERFACE nm_sort
     SUBROUTINE nm_sort_sp(arr)
       USE nm_type
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE nm_sort_sp
     SUBROUTINE nm_sort_dp(arr)
       USE nm_type
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE nm_sort_dp
     SUBROUTINE nm_sort_i(arr)
       USE nm_type
       INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE nm_sort_i
  END INTERFACE nm_sort

  INTERFACE nm_sort2
     SUBROUTINE nm_sort2_dp(ARR1,ARR2)
       USE nm_type
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: ARR1,ARR2
     END SUBROUTINE nm_sort2_dp
     SUBROUTINE nm_sort2_sp(ARR1,ARR2)
       USE nm_type
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: ARR1,ARR2
     END SUBROUTINE nm_sort2_sp
  END INTERFACE nm_sort2


END MODULE nm_nm
