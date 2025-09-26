SUBROUTINE legendre(n,t, p, diffp)
    USE nm_type
    IMPLICIT NONE    
    REAL(DP), INTENT(IN) :: t
    INTEGER(I4B), INTENT(IN) :: n    
    REAL(DP),  INTENT(OUT) :: p, diffp
    REAL(DP)  :: p0, p1, eq1, eq2
    INTEGER(I4B) :: k
    

    p0 = 1.0_dp
    p1 = t
    DO k = 1, n -1 
       p = ((2.0_dp*k + 1.0_dp)*t*p1 - k*p0)/(1.0_dp + k )  
       p0 = p1
       p1 = p
    
    END DO       
    !print *, 'p, p0, p1 in legendre : ', p, p0, p1

    diffp = n*(p0 - t*p1)/(1.0_dp - t**2.0_dp)
   
END SUBROUTINE legendre


SUBROUTINE nm_gaussNodes(x1, x2, x, A)
    !!!!  Returns nodal abscissas {x(1...n)} and weights {A(1...n)}
    !!!!  of gauss legendre n-point quadrature, between the lower and upper
    !!!!  limits of integration x1 and x2
    !!!!  Please refer to "Numerical Methods in Engineering with Python 3"
    !!!!  p: 224-225 , for more information
    USE nm_type
    USE nm_util, only :  check_length

    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2
    REAL(DP), DIMENSION(:), INTENT(OUT) :: x,A
    REAL(DP), PARAMETER :: EPS=3.0D-14
    INTEGER(I4B) :: i, j, nRoots, n
    REAL(DP) :: xl,xm, t, dt, p, diffp

    n = check_length(size(x),size(A),'gaussNodes')

    xm = 0.5_dp*(x2+x1)
    xl = 0.5_dp*(x2-x1)
    nRoots = int((n+1)/2 )   ! number of non_neg roots

    DO  i = 1,nRoots
        t = cos(PI_D*(i - 0.25_dp)/(n+0.5_dp))       ! approx root
        !print *, 't : ',t
        DO j= 1,30
            CALL legendre(n,t, p, diffp)       ! Newton -raphson
         !   print *, 'p,diff ', p ,diffp 
            dt = -p/diffp
            t = t + dt
            IF ( abs(dt) < EPS  )  THEN
                x(i) = xm-xl*t
                x(n - i + 1) = xm+xl*t
                A(i) = 2.0_dp*xl/((1.0_dp - t**2.0_dp)*(diffp**2.0_dp) )       !  Eq(6.25)
                A(n- i +1) = A(i)
                EXIT
            END IF

        END DO
    END DO

END SUBROUTINE nm_gaussNodes
