
MODULE MMCRAD1D_INTERPOLATE

  ! MC merged mod_interp (mod_interp.f) 
  ! into the prsent file
  IMPLICIT NONE

  PUBLIC :: spline2, MCLININTPOL, MCBILININTPOL, MCTRILININTPOL

  PRIVATE :: bnsrch !, splint, spline, myspline, spline3, splint3

CONTAINS

  !=================================

  FUNCTION MCTRILININTPOL(fint, xint, yint, zint, nx, ny, nz, x, y, z)

    ! tri-linear interpolation 
    ! xint yint and zint assumed increasing
    ! NO EXTRAPOLATION 

    IMPLICIT NONE

    REAL (8) :: MCTRILININTPOL

    INTEGER, INTENT (IN) :: nx
    INTEGER, INTENT (IN) :: ny
    INTEGER, INTENT (IN) :: nz
    REAL(8), INTENT (IN) :: fint(nx,ny,nz)
    REAL(8), INTENT (IN) :: xint(nx)
    REAL(8), INTENT (IN) :: yint(ny)
    REAL(8), INTENT (IN) :: zint(nz)
    REAL(8), INTENT (IN) :: x
    REAL(8), INTENT (IN) :: y
    REAL(8), INTENT (IN) :: z
    INTEGER :: i, j, k
    REAL(8) :: t, u, v

    IF( (x < xint(1)) .OR. (x > xint(nx)) ) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. TRILININTPOL)  : ERROR '
       WRITE(*,*) '                          Extrapolation not allowed'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       write(*,*) 'x=',x
       write(*,*) xint(1)
       write(*,*) xint(nx)
       STOP 
    ENDIF
    IF( (y < yint(1)) .OR. (y > yint(ny)) ) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. TRILININTPOL)  : ERROR '
       WRITE(*,*) '                          Extrapolation not allowed'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       write(*,*) 'y=',y
       write(*,*) yint(1)
       write(*,*) yint(ny)
       STOP 
    ENDIF
    IF( (z < zint(1)) .OR. (z > zint(nz)) ) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. TRILININTPOL)  : ERROR '
       WRITE(*,*) '                          Extrapolation not allowed'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       write(*,*) 'z=',z
       write(*,*) zint(1)
       write(*,*) zint(nz)
       STOP 
    ENDIF

    i = 1
    DO WHILE (xint(i) <= x .AND. i /= nx)
       i = i + 1
    ENDDO
    i = i - 1

    j = 1
    DO WHILE (yint(j) <= y .AND. j /= ny)
       j = j + 1
    ENDDO
    j = j - 1

    k = 1
    DO WHILE (zint(k) <= z .AND. k /= nz)
       k = k + 1
    ENDDO
    k = k - 1

    t = (x-xint(i))/(xint(i+1)-xint(i))
    u = (y-yint(j))/(yint(j+1)-yint(j))
    v = (z-zint(k))/(zint(k+1)-zint(k))

    MCTRILININTPOL = (fint(i  ,j  ,k  ) * (1.0D0 - t) * (1.0D0 - u) * (1.0D0 - v)) &
         +         (fint(i+1,j  ,k  ) *    t        * (1.0D0 - u) * (1.0D0 - v))   &
         +         (fint(i  ,j+1,k  ) * (1.0D0 - t) *     u       * (1.0D0 - v))   &
         +         (fint(i  ,j  ,k+1) * (1.0D0 - t) * (1.0D0 - u) *       v    )   &
         +         (fint(i+1,j  ,k+1) *     t       * (1.0D0 - u) *       v    )   &
         +         (fint(i  ,j+1,k+1) * (1.0D0 - t) *     u       *       v    )   &
         +         (fint(i+1,j+1,k  ) *     t       *     u       * (1.0D0 - v))   &
         +         (fint(i+1,j+1,k+1) *     t       *     u       *       v    )

  END FUNCTION MCTRILININTPOL

  !=================================

  FUNCTION MCBILININTPOL(fint, xint, yint, nx, ny, x, y)

    ! bi-linear interoplation 
    ! xint and yint assumed increasing
    ! NO EXTAPOLATION 

    IMPLICIT NONE

    REAL (8) :: MCBILININTPOL

    INTEGER, INTENT (IN) :: nx
    INTEGER, INTENT (IN) :: ny
    REAL(8), INTENT (IN) :: fint(nx,ny)
    REAL(8), INTENT (IN) :: xint(nx)
    REAL(8), INTENT (IN) :: yint(ny)
    REAL(8), INTENT (IN) :: x
    REAL(8), INTENT (IN) :: y

    INTEGER :: i, j
    REAL(8) :: t, u

    IF( (x < xint(1)) .OR. (x > xint(nx)) ) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. BILININTPOL)  : ERROR '
       WRITE(*,*) '                          Extrapolation not allowed'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       write(*,*) 'x=',x
       write(*,*) xint(1)
       write(*,*) xint(nx)
       STOP 
    ENDIF
    IF( (y < yint(1)) .OR. (y > yint(ny)) ) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. BILININTPOL)  : ERROR '
       WRITE(*,*) '                          Extrapolation not allowed'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       write(*,*) 'y=',y
       write(*,*) yint(1)
       write(*,*) yint(ny)
       STOP 
    ENDIF

    i = 1
    DO WHILE (xint(i) <= x .AND. i /= nx)
       i = i + 1
    ENDDO
    i = i - 1

    j = 1
    DO WHILE (yint(j) <= y .AND. j /= ny)
       j = j + 1
    ENDDO
    j = j - 1

    t = (x-xint(i))/(xint(i+1)-xint(i))
    u = (y-yint(j))/(yint(j+1)-yint(j))

    MCBILININTPOL = ((1.0D0-t) * (1.0D0-u) * fint(i  , j)   ) &
         +        ( t        * (1.0D0-u) * fint(i+1, j)   ) &
         +        ( t        *    u      * fint(i+1, j+1) ) &
         +        ((1.0D0-t) *    u      * fint(i  , j+1) )

  END FUNCTION MCBILININTPOL

  !=================================

  FUNCTION MCLININTPOL(fint, xint, nx, xess)
    ! linear interoplation of fint @ xess
    ! xint assumed increasing, NO EXTAPOLATION 
    ! make sure xess belongs to [xint(1),xint(nx)]

    IMPLICIT NONE

    REAL (8) :: MCLININTPOL

    INTEGER, INTENT (IN)  :: nx
    REAL (8), INTENT (IN) :: fint(nx)
    REAL (8), INTENT (IN) :: xint(nx)
    REAL (8), INTENT (IN) :: xess
    INTEGER               :: i

    i = 1
    DO WHILE (xint(i) <= xess .AND. i /= nx)
       i = i + 1
    ENDDO
    i = i - 1

    ! NO extrapolation
    IF( (xess < xint(1)) .OR. (xess > xint(nx)) ) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. MCLININTPOL)  : ERROR '
       WRITE(*,*) '                          Extrapolation not allowed'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       write(*,*) 'x=',xess
       write(*,*) xint(1)
       write(*,*) xint(nx)
       STOP 
    ELSE
       MCLININTPOL = fint(i) * (xint(i+1) - xess) + fint(i+1) * (xess-xint(i))
       MCLININTPOL = MCLININTPOL / (xint(i+1) - xint(i))
    ENDIF

  END FUNCTION MCLININTPOL

  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE  SPLINE2(NTAU, TAU, Y, X, VAL, DERIV)

    !       NTAU= NUMBER OF DATA POINTS                   (SCALAR)
    !       X= WHERE THE INTERPOLATION IS TO BE PERFORMED (SCALAR)
    !       TAU= VECTOR OF X COORDINATES
    !       Y= DATA (NOT NECESSARILY EVENLY SPACED )      (VECTOR)
    !       VAL= INTERPOLATED VALUE OF FUNCTION           (SCALAR)
    IMPLICIT NONE

    INTEGER                                 :: M, NTAU
    REAL(8), dimension(NTAU)                :: TAU, Y
    REAL(8)                                 :: X
    REAL(8)                                 :: T1, T2, T3, T4, VAL, DERIV
    REAL(8)                                 :: F1, F2, F3, F4, A, B, C, D, DET

    call bnsrch(ntau ,tau, x, M  )
    !<TMP>
    !    WRITE(6,*)'M=',M
    !    WRITE(6,*)'Y(M-1), Y(M)=',Y(M-1),Y(M)
    !</TMP>

    ! check that x is always inside the permissible range of tau

    if (M .lt. (NTAU-1) .and. M .ne. 1) then
       T1=TAU(M-1)
       T2=TAU(M )
       T3=TAU(M+1)
       T4=TAU(M+2)
       F1=Y(M-1)
       F2=Y(M)
       F3=Y(M+1)
       F4=Y(M+2)
    endif

    if (M .EQ. NTAU-1) THEN
       T1=TAU(M-2)
       T2=TAU(M-1)
       T3=TAU(M)
       T4=TAU(M+1)
       F1=Y(M-2)
       F2=Y(M-1)
       F3=Y(M)
       F4=Y(M+1)
    endif

    if ((M .eq. NTAU) .OR. (M .EQ. NTAU-1)) THEN
       T1=TAU(ntau)
       T2=TAU(ntau-1)
       T3=TAU(ntau-2)
       T4=TAU(ntau-3)
       F1=Y(ntau)
       F2=Y(ntau-1)
       F3=Y(ntau-2)
       F4=Y(ntau-3)
    endif

    if (M .eq. 1)THEN
       T1=TAU(1)
       T2=TAU(2)
       T3=TAU(3)
       T4=TAU(4)
       F1=Y(1)
       F2=Y(2)
       F3=Y(3)
       F4=Y(4)
    endif

    DET =(T1-T4)*(T2-T3)*( (T1*T4)**2 +(T2*T3)**2) +                       &
         (T1-T2)*(T3-T4)*( (T1*T2)**2 +(T3*T4)**2) +                       &
         (T1-T3)*(T4-T2)*( (T1*T3)**2 +(T2*T4)**2)                    


    A= (F1-F2)*(T3-T4)*T3*T4 + (F3-F1)*(T2-T4)*T2*T4 +                     &
         (F1-F4)*(T2-T3)*T2*T3 + (F2-F3)*(T1-T4)*T1*T4 +                &
         (F4-F2)*(T1-T3)*T1*T3 + (F3-F4)*(T1-T2)*T1*T2

    A=A/DET

    B= (F2-F1)*(T3*T3-T4*T4)*T3*T4+(F1-F3)*(T2*T2-T4*T4)*T2*T4+            &
         (F4-F1)*(T2*T2-T3*T3)*T2*T3+(F3-F2)*(T1*T1-T4*T4)*T1*T4+            &
         (F2-F4)*(T1*T1-T3*T3)*T1*T3+(F4-F3)*(T1*T1-T2*T2)*T1*T2

    B=B/DET

    C= (F1-F2)*(T3-T4)*(T3*T4)**2 +(F3-F1)*(T2-T4)*(T2*T4)**2+             &
         (F1-F4)*(T2-T3)*(T2*T3)**2 +(F2-F3)*(T1-T4)*(T1*T4)**2+             &
         (F4-F2)*(T1-T3)*(T1*T3)**2 +(F3-F4)*(T1-T2)*(T1*T2)**2

    C=C/DET

    D=(T1*F2-T2*F1)*(T3-T4)*(T3*T4)**2 +                                   &
         (T3*F1-T1*F3)*(T2-T4)*(T2*T4)**2 +                                   &
         (T1*F4-T4*F1)*(T2-T3)*(T2*T3)**2 +                                   &
         (T2*F3-T3*F2)*(T1-T4)*(T1*T4)**2 +                                   &
         (T3*F4-T4*F3)*(T1-T2)*(T1*T2)**2 +                                   &
         (T4*F2-T2*F4)*(T1-T3)*(T1*T3)**2

    D=D/DET

    VAL=X*( X*(A*X+B)+C)+D
    DERIV=X*( X*(3.*A+2.*B))+C

    RETURN

  END SUBROUTINE SPLINE2


  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE bnsrch(ndat,xp,x,ip)

    ! BINARY search to locate index of coordinate
    ! This is the standard binary search routine and can be found in 
    ! many standard references.
    !
    ! Inputs:
    !      NDAT = dimension of XP
    !      XP = array containing vertex coordinates
    !      X  = Point to be interpolated
    !
    ! Outputs:
    !      IP = index of desired coordinates
    !
    !  written by Philip Gabriel 18 Jan 1996

    IMPLICIT NONE
    integer                                 :: hi, lo, ndat, ip
    real(8)  , dimension(ndat)              :: xp
    real(8)                                 :: x

10  hi=ndat
    lo=1   
    ip=(lo+hi)/2
20  if(hi-lo.le.1) go to 50
    if(xp(ip)-x)30,9000,40
30  lo=ip
    ip=(lo+hi)/2
    if(xp(lo).lt.xp(ip)) goto 20
40  hi=ip
    ip=(lo+hi)/2
    if(xp(ip).lt.xp(hi)) goto 20
50  if(xp(lo).eq.x) goto 100
    if(x-xp(lo+1))100,70,70
70  ip=lo+1
    goto 9000
100 ip=lo
9000 return

  END SUBROUTINE bnsrch

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

!!$  SUBROUTINE myspline(n,x,y,nn,xn,yn)
!!$    !---------spline fit to derive the yn value at point xn
!!$    !  Inputs:
!!$    !     n:    the length of x and y
!!$    !     x(n): the x values which x(1) < x(2) ... < x(n)
!!$    !     y(n): the y value which correspondent to x(n)
!!$    !     nn:  the length of vector xx and yy
!!$    !     xn:  the x value at which y value is wanted
!!$    !
!!$    !  Outputs:
!!$    !     yn: the wanted y value from the fitting
!!$    !
!!$    !  Internal variables:
!!$    !     yp1: the derivative of y over x at x(1), for natural bc, yp1=1.e31
!!$    !     ypn: the derivative of y over x at x(n), for natural bc, ypn=1.e31
!!$    !     y2(n): the second derivatives
!!$    !
!!$    IMPLICIT NONE
!!$    integer, parameter                      :: ny2 = 5000
!!$    integer                                 :: n, nn, i
!!$
!!$    real(8)                                 :: xx, yy, yp1, ypn
!!$    real(8), dimension(n)                   :: x, y
!!$    real(8), dimension(nn)                  :: xn, yn
!!$    real(8), dimension(ny2)                 :: y2
!!$
!!$    !--------the sorting which makes sure x(1)<x(2)<...<x(n)-------
!!$    !    CALL sort2(n,x,y)
!!$    !--------start spline------------
!!$    yp1 = 1.D31
!!$    ypn = 1.D31
!!$    CALL spline3(x,y,n,yp1,ypn,y2)
!!$
!!$    DO i = 1,nn
!!$       xx = xn(i)
!!$       CALL splint3(x,y,y2,n,xx,yy)
!!$       yn(i) = yy
!!$    ENDDO !i
!!$
!!$    RETURN
!!$
!!$  END SUBROUTINE myspline

  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

!!$  SUBROUTINE spline3(x,y,n,yp1,ypn,y2)
!!$
!!$    implicit none
!!$
!!$    INTEGER, PARAMETER                      :: NMAX = 5000, ny2 = 5000
!!$    INTEGER                                 :: n, i, k
!!$
!!$    REAL(8)                                 :: yp1, ypn, p, qn, sig, un
!!$    REAL(8), DIMENSION(n)                   :: x, y
!!$    REAL(8), DIMENSION(ny2)                 :: y2
!!$    REAL(8), DIMENSION(NMAX)                :: u
!!$
!!$    IF (yp1 .gt. .99D30) THEN
!!$       y2(1) = 0.D0
!!$       u(1) = 0.D0
!!$    ELSE
!!$       y2(1) = -0.5D0
!!$       u(1) = (3.D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
!!$    ENDIF
!!$
!!$    DO i = 2,n-1
!!$       sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
!!$       p = sig*y2(i-1)+2.D0
!!$       y2(i) = (sig-1.D0)/p
!!$       u(i) = (6.D0*((y(i+1)-y(i))/(x(i+1)                                   &
!!$            -x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))             &
!!$            -sig*u(i-1))/p
!!$    ENDDO
!!$
!!$    IF (ypn .gt. .99D30) THEN
!!$       qn = 0.D0
!!$       un = 0.D0
!!$    ELSE
!!$       qn = 0.5D0
!!$       un = (3.D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
!!$    ENDIF
!!$
!!$    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.D0)
!!$    DO k = n-1,1,-1
!!$       y2(k) = y2(k)*y2(k+1)+u(k)
!!$    ENDDO
!!$
!!$    RETURN
!!$
!!$  END SUBROUTINE spline3

  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

!!$  SUBROUTINE splint3(xa,ya,y2a,n,x,y)
!!$
!!$    implicit none
!!$
!!$    INTEGER, PARAMETER                      :: ny2 = 5000
!!$    INTEGER                                 :: n, k, khi, klo
!!$
!!$    REAL(8)                                 :: a, b, h
!!$    REAL(8)                                 :: x, y
!!$    REAL(8), DIMENSION(n)                   :: xa, ya 
!!$    REAL(8), DIMENSION(ny2)                 :: y2a
!!$
!!$    klo = 1
!!$    khi = n
!!$1   IF (khi-klo .gt. 1) THEN
!!$       k=(khi+klo)/2
!!$       IF (xa(k) .gt. x) THEN
!!$          khi=k
!!$       ELSE
!!$          klo=k
!!$       ENDIF
!!$       GOTO 1
!!$    ENDIF
!!$
!!$    h = xa(khi)-xa(klo)
!!$
!!$    IF (h .eq. 0.) pause 'bad xa input in splint3'
!!$
!!$    a = (xa(khi)-x)/h
!!$    b = (x-xa(klo))/h
!!$    y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo) +                            &
!!$         (b**3-b)*y2a(khi))*(h*h)/6.D0
!!$
!!$    RETURN
!!$  END SUBROUTINE splint3

  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

!!$  SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
!!$    IMPLICIT NONE
!!$    INTEGER :: N
!!$    REAL(8) :: XA(N),YA(N),Y2A(N)
!!$    REAL(8) :: X,Y,H, A, B
!!$    INTEGER :: KLO, KHI, K
!!$
!!$    KLO=1
!!$    KHI=N
!!$1   IF (KHI-KLO.GT.1) THEN
!!$       K=(KHI+KLO)/2
!!$       IF(XA(K).GT.X)THEN
!!$          KHI=K
!!$       ELSE
!!$          KLO=K
!!$       ENDIF
!!$       GOTO 1
!!$    ENDIF
!!$    H=XA(KHI)-XA(KLO)
!!$    IF (H.EQ.0.) PAUSE 'Bad XA input. in splint'
!!$    A=(XA(KHI)-X)/H
!!$    B=(X-XA(KLO))/H
!!$    Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
!!$    RETURN
!!$  END SUBROUTINE SPLINT

  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

!!$  SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
!!$
!!$    IMPLICIT NONE
!!$    INTEGER,PARAMETER :: NMAX=2000
!!$    INTEGER :: N
!!$    INTEGER :: K,I
!!$    REAL(8) :: X(N),Y(N),Y2(N),U(NMAX)
!!$    REAL(8) :: YP1,YPN, SIG
!!$    REAL(8) :: QN, UN, P
!!$
!!$    IF (YP1.GT..99E30) THEN
!!$       Y2(1) = 0.
!!$       U(1)  = 0.
!!$    ELSE
!!$       Y2(1) = -0.5
!!$       U(1)  = (3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
!!$    ENDIF
!!$    DO I=2,N-1
!!$       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
!!$       P=SIG*Y2(I-1)+2.
!!$       Y2(I)=(SIG-1.)/P
!!$       U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/&
!!$            (X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
!!$    ENDDO
!!$    IF (YPN.GT..99E30) THEN
!!$       QN=0.
!!$       UN=0.
!!$    ELSE
!!$       QN=0.5
!!$       UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
!!$    ENDIF
!!$    Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
!!$    DO K=N-1,1,-1
!!$       Y2(K)=Y2(K)*Y2(K+1)+U(K)
!!$    ENDDO
!!$    RETURN
!!$  END SUBROUTINE SPLINE


END MODULE MMCRAD1D_INTERPOLATE
