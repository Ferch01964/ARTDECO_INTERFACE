      SUBROUTINE get_index_dp(ARR, ind)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  Returns the indices that would sort an array
      !!  It returns an array of indices of the same shape as ARR
      !!  (of type real and of kind DP) that index data along the 
      !!  given axis in sorted order
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE nm_type
        USE nm_util, only :  SWAP,check_length
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: ARR
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: ind
        INTEGER(I4B) :: N, i,j
        N = check_length(size(ARR),size(ind),'indexx_dp')
        ind = (/ (i, i= 1 , N) /)
        DO i = 1 , N
          DO j = i+1 , N
            IF (ARR(ind(i)) > ARR(ind(j))) then
                CALL   SWAP(ind(i), ind(j))                
            ENDIF
          ENDDO
        ENDDO

      END SUBROUTINE get_index_dp


      SUBROUTINE get_index_sp(ARR, ind)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  Returns the indices that would sort an array
      !!  It returns an array of indices of the same shape as ARR
      !!  (of type real and of kind SP) that index data along the 
      !!  given axis in sorted order
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE nm_type
        USE nm_util, only :  SWAP,check_length
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: ARR
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: ind
        INTEGER(I4B) :: N, i,j
        N = check_length(size(ARR),size(ind),'indexx_sp')
        ind = (/ (i, i= 1 , N) /)
        DO i = 1 , N
          DO j = i+1 , N
            IF (ARR(ind(i)) > ARR(ind(j))) then
                CALL   SWAP(ind(i), ind(j))                
            ENDIF
          ENDDO
        ENDDO

      END SUBROUTINE get_index_sp


      SUBROUTINE get_index_i4b(ARR, ind)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  Returns the indices that would sort an array
      !!  It returns an array of indices of the same shape as ARR
      !!  (of type integer and of kind I4B) that index data along the 
      !!  given axis in sorted order
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE nm_type
        USE nm_util, only :  SWAP,check_length
        IMPLICIT NONE
        INTEGER(I4B),DIMENSION(:), INTENT(IN) :: ARR
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: ind
        INTEGER(I4B) :: N, i,j
        N = check_length(size(ARR),size(ind),'indexx_i4b')
        ind = (/ (i, i= 1 , N) /)
        DO i = 1 , N
          DO j = i+1 , N
            IF (ARR(ind(i)) > ARR(ind(j))) then
                CALL   SWAP(ind(i), ind(j))                
            ENDIF
          ENDDO
        ENDDO

      END SUBROUTINE get_index_i4b





      
