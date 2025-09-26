MODULE nm_util
  USE nm_type
  INTERFACE swap
     MODULE PROCEDURE swap_i,swap_r, swap_r_dp ,swap_rv, &
          swap_iv, swap_rv_dp, swap_c,&
          swap_cv,swap_cm,swap_z,swap_zv,swap_zm,swap_im
          
  END INTERFACE swap

  INTERFACE order
     MODULE PROCEDURE order_i, order_r_dp ,order_r_sp          
  END INTERFACE order

 
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    SWAP Module         !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE swap_i(a,b)
    !!!! Permute between two integers a and b
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i
  !BL
  SUBROUTINE swap_r(a,b)
    !!!! Permute between two reels of kind sp a and b
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r
  !BL
  SUBROUTINE swap_r_dp(a,b)
    !!!! Permute between two reels of kind dp a and b
    REAL(DP), INTENT(INOUT) :: a,b
    REAL(DP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r_dp
  !BL
  SUBROUTINE swap_rv(a,b)
    !!!! Permute between two arrays a and b of type real and kind sp
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv
  !BL
  SUBROUTINE swap_iv(a,b)
    !!!! Permute between two arrays a and b of type integer and kind I4B
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: a,b
    INTEGER(I4B), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_iv
  SUBROUTINE swap_rv_dp(a,b)
    !!!! Permute between two arrays a and b of type real and kind dp
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(DP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv_dp
  !BL
  SUBROUTINE swap_c(a,b)
    !!!! Permute between two arrays a and b of type complex and kind SPC
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_c
  !BL
  SUBROUTINE swap_cv(a,b)
    !!! swap of array of type complex with a vector
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cv
  !BL
  SUBROUTINE swap_cm(a,b)
     !!!! Permute between two matrices a and b of type complex and kind SPC
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cm
  !BL
  SUBROUTINE swap_z(a,b)
    !!!! Permute between two complex of kind DPC  a and b 
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    COMPLEX(DPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z
  !BL
  SUBROUTINE swap_zv(a,b)
   !!!! Permute between two arrays a and b of type complex and kind DPC
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv
  !BL
  SUBROUTINE swap_zm(a,b)
    !!!! Permute between two matrices a and b of type complex and kind DPC
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm
  !BL

  SUBROUTINE swap_im(a,b)
    !!!! Permute between two matrices a and b of type integer and kind I4B
    INTEGER(I4B), DIMENSION(:,:), INTENT(INOUT) :: a,b
    INTEGER(I4B), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_im


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    ORDER Module        !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine Order_r_dp(p,q)
    !return p,q in ascending order
    REAL(DP), INTENT(INOUT) :: p,q
    if (p > q) then
        call Swap_r_dp(p, q)
    end if
  END SUBROUTINE Order_r_dp

  Subroutine Order_r_sp(p,q)
    !return p,q in ascending order
    REAL(SP), INTENT(INOUT) :: p,q
    if (p > q) then
        call Swap_r(p, q)
    end if
  END SUBROUTINE Order_r_sp

  Subroutine Order_i(p,q)
    !return p,q in ascending order
    INTEGER(I4B), INTENT(INOUT) :: p,q
    if (p > q) then
        call Swap_i(p, q)
    end if
  END SUBROUTINE Order_i





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    check length module   !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION check_length(n1,n2,string)
  !!! function used to check if the length n1
  !!! is eqaul  to n2 otherwise error 
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: check_length
    if (n1 == n2) then
       check_length=n1
    else
       write (*,*) 'Error: a check_length failed with this tag:', &
            string
       STOP 'program terminated by check_length'
    end if
  END FUNCTION check_length

  


END MODULE nm_util
