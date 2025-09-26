SUBROUTINE nm_sort_dp(ARR)
    !Bubble sorting of integer array A
    USE nm_type
    USE nm_util, only :  Order
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: ARR
    INTEGER(I4B) :: n,i,j
    n = size(ARR)
    do i=1, n
        do j=n, i+1, -1
            call Order(ARR(j-1), ARR(j))
        end do
    end do
  
END SUBROUTINE nm_sort_dp

SUBROUTINE nm_sort_sp(ARR)
    !Bubble sorting of integer array A
    USE nm_type
    USE nm_util, only :  Order
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: ARR
    INTEGER(I4B) :: n,i,j
    n = size(ARR)
    do i=1, n
        do j=n, i+1, -1
            call Order(ARR(j-1), ARR(j))
        end do
    end do
  
END SUBROUTINE nm_sort_sp


SUBROUTINE nm_sort_i(ARR)
    !Bubble sorting of integer array A
    USE nm_type
    USE nm_util, only :  Order
    IMPLICIT NONE
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: ARR
    INTEGER(I4B) :: n,i,j
    n = size(ARR)
    do i=1, n
        do j=n, i+1, -1
            call Order(ARR(j-1), ARR(j))
        end do
    end do
  
END SUBROUTINE nm_sort_i
