

MODULE MUTILITY

  USE MCONSTANTS, only : dp, max_len, xpi

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GAUSSIAN, REVERSE, LININTPOL, XINTEG2, PLINT, TRIMCAT, &
       TRILININTPOL, BILININTPOL, spline, splint

CONTAINS 


  FUNCTION TRILININTPOL(fint, xint, yint, zint, nx, ny, nz, x, y, z)

    ! tri-linear interpolation 
    ! xint yint and zint assumed increasing
    ! NO EXTRAPOLATION 

    IMPLICIT NONE

    REAL (kind=dp) :: TRILININTPOL

    INTEGER, INTENT (IN) :: nx
    INTEGER, INTENT (IN) :: ny
    INTEGER, INTENT (IN) :: nz
    REAL(kind=dp), INTENT (IN) :: fint(nx,ny,nz)
    REAL(kind=dp), INTENT (IN) :: xint(nx)
    REAL(kind=dp), INTENT (IN) :: yint(ny)
    REAL(kind=dp), INTENT (IN) :: zint(nz)
    REAL(kind=dp), INTENT (IN) :: x
    REAL(kind=dp), INTENT (IN) :: y
    REAL(kind=dp), INTENT (IN) :: z
    INTEGER :: i, j, k
    REAL(kind=dp) :: t, u, v

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

    TRILININTPOL = (fint(i  ,j  ,k  ) * (1.0D0 - t) * (1.0D0 - u) * (1.0D0 - v)) &
         +         (fint(i+1,j  ,k  ) *    t        * (1.0D0 - u) * (1.0D0 - v)) &
         +         (fint(i  ,j+1,k  ) * (1.0D0 - t) *     u       * (1.0D0 - v)) &
         +         (fint(i  ,j  ,k+1) * (1.0D0 - t) * (1.0D0 - u) *       v    ) &
         +         (fint(i+1,j  ,k+1) *     t       * (1.0D0 - u) *       v    ) &
         +         (fint(i  ,j+1,k+1) * (1.0D0 - t) *     u       *       v    ) &
         +         (fint(i+1,j+1,k  ) *     t       *     u       * (1.0D0 - v)) &
         +         (fint(i+1,j+1,k+1) *     t       *     u       *       v    )

  END FUNCTION TRILININTPOL

  !=================================

  FUNCTION BILININTPOL(fint, xint, yint, nx, ny, x, y)

    ! bi-linear interoplation 
    ! xint and yint assumed increasing
    ! NO EXTAPOLATION 

    IMPLICIT NONE

    REAL (kind=dp) :: BILININTPOL

    INTEGER, INTENT (IN) :: nx
    INTEGER, INTENT (IN) :: ny
    REAL(kind=dp), INTENT (IN) :: fint(nx,ny)
    REAL(kind=dp), INTENT (IN) :: xint(nx)
    REAL(kind=dp), INTENT (IN) :: yint(ny)
    REAL(kind=dp), INTENT (IN) :: x
    REAL(kind=dp), INTENT (IN) :: y

    INTEGER :: i, j
    REAL(kind=dp) :: t, u

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

    BILININTPOL = ((1.0D0-t) * (1.0D0-u) * fint(i  , j)   ) &
         +        ( t        * (1.0D0-u) * fint(i+1, j)   ) &
         +        ( t        *    u      * fint(i+1, j+1) ) &
         +        ((1.0D0-t) *    u      * fint(i  , j+1) )

  END FUNCTION BILININTPOL

  !--------------------------------------------------------------------------

  FUNCTION LININTPOL(fint, xint, ni, xess)
    ! linear interoplation of fint @ xess
    ! xint assumed increasing, NO EXTAPOLATION 
    ! make sure xess belongs to [xint(1),xint(ni)]

    IMPLICIT NONE

    REAL (KIND=dp) :: LININTPOL

    INTEGER, INTENT (IN)        :: ni
    REAL (KIND=dp), INTENT (IN) :: fint(ni)
    REAL (KIND=dp), INTENT (IN) :: xint(ni)
    REAL (KIND=dp), INTENT (IN) :: xess
    INTEGER                     :: i

    i = 1
    DO WHILE (xint(i) <= xess .AND. i /= ni)
       i = i + 1
    ENDDO
    i = i - 1

    ! NO extrapolation
    IF( (xess < xint(1)) .OR. (xess > xint(ni)) ) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. INTPOL)  : ERROR '
       WRITE(*,*) ' Extrapolation not allowed'
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP 
    ELSE
       LININTPOL = fint(i) * (xint(i+1) - xess) + fint(i+1) * (xess-xint(i))
       LININTPOL = LININTPOL / (xint(i+1) - xint(i))
    ENDIF

  END FUNCTION LININTPOL

  !--------------------------------------------------------------------------

  FUNCTION REVERSE(n, tab)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp), INTENT(IN) :: tab(n)
    REAL(KIND=dp) :: reverse(n)
    INTEGER :: i

    DO i = 1, n 
       reverse(n+1-i) = tab(i) 
    END DO

  END FUNCTION REVERSE


  !--------------------------------------------------------------------------

  FUNCTION GAUSSIAN(x, mu, sigma)

    ! This function gives as a result
    ! the value of the normalized gaussian function in x
    ! mu is the expected value (mode)
    ! sigma is the standard deviation

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(IN) :: x
    REAL(KIND=dp), INTENT(IN) :: mu
    REAL(KIND=dp), INTENT(IN) :: sigma
    REAL(KIND=dp) :: gaussian

    !---------------

    gaussian = ( 1.0_dp / ( sigma * SQRT(2.0_dp * xpi)) ) * EXP( - ((x - mu)**2.0_dp) / (2.0_dp * sigma**2.0_dp) )

  END FUNCTION GAUSSIAN

  !--------------------------------------------------------------------------

  FUNCTION XINTEG2(imin, imax, n, xin, yin)
    ! computes integral of yin, variable step 
    ! make sure you have at least two points of integration
    ! Trapeze method

    IMPLICIT NONE

    REAL (KIND=dp) :: XINTEG2

    INTEGER, INTENT (IN)        :: imin
    INTEGER, INTENT (IN)        :: imax
    INTEGER, INTENT (IN)        :: n
    REAL (KIND=dp), INTENT (IN) :: xin(n)
    REAL (KIND=dp), INTENT (IN) :: yin(n)
    REAL (KIND=dp)              :: xa, ya, xb, yb
    REAL (KIND=dp)              :: primitaux
    INTEGER                     :: i

    xa = xin(imin)
    ya = yin(imin)
    primitaux = 0.0_dp
    DO i=imin+1,imax
       xb = xin(i)
       yb = yin(i)
       primitaux = primitaux + (xb - xa) * (ya + yb)
       xa = xb
       ya = yb
    ENDDO
    XINTEG2 = primitaux * 0.5_dp

  END FUNCTION XINTEG2

  !--------------------------------------------------------------------------

  subroutine plint ( ftab, xtab, ntab, a_in, b_in, result )

    ! PLINT approximates the integral of unequally spaced data.
    !
    !  Discussion:
    !
    !    The method uses piecewise linear interpolation.
    !
    !  Modified:
    !
    !    10 February 2006
    !
    !  Reference:
    !
    !    Philip Davis, Philip Rabinowitz,
    !    Methods of Numerical Integration,
    !    Second Edition,
    !    Dover, 2007,
    !    ISBN: 0486453391,
    !    LC: QA299.3.D28.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NTAB, the number of entries in FTAB and
    !    XTAB.  NTAB must be at least 2.
    !
    !    Input, real ( kind = 8 ) XTAB(NTAB), the abscissas at which the
    !    function values are given.  The XTAB's must be distinct
    !    and in ascending order.
    !
    !    Input, real ( kind = 8 ) FTAB(NTAB), the function values, 
    !    FTAB(I) = F(XTAB(I)).
    !
    !    Input, real ( kind = 8 ) A, the lower limit of integration.  A should
    !    be, but need not be, near one endpoint of the interval
    !    (X(1), X(NTAB)).
    !
    !    Input, real ( kind = 8 ) B, the upper limit of integration.  B should
    !    be, but need not be, near one endpoint of the interval
    !    (X(1), X(NTAB)).
    !
    !    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.

    implicit none

    integer, intent(in) :: ntab
    real ( kind = dp ), intent(in) :: ftab(ntab)
    real ( kind = dp ), intent(in) :: xtab(ntab)
    real ( kind = dp ), intent(in) :: a_in
    real ( kind = dp ), intent(in) :: b_in
    real ( kind = dp ), intent(out) :: result

    real ( kind = dp ) :: fa
    real ( kind = dp ) :: fb
    real ( kind = dp ) :: a
    real ( kind = dp ) :: b   
    integer :: i
    integer :: ihi
    integer :: ilo
    integer :: ind
    real ( kind = dp ) :: slope
    real ( kind = dp ) :: syl


    a = a_in
    b = b_in
    result = 0.0_dp
    !
    !  Check the parameters:
    !
    if ( ntab < 2 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PLINT - Fatal error!'
       write ( *, '(a,i8)' ) '  NTAB < 2, NTAB = ', ntab
       stop
    end if

    do i = 2, ntab
       if ( xtab(i) <= xtab(i-1) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PLINT - Fatal error!'
          write ( *, '(a)' ) '  Nodes not in strict increasing order.'
          write ( *, '(a,i8)' ) '  XTAB(I) <= XTAB(I-1) for I = ', i
          write ( *, '(a,g14.6)' ) '  XTAB(I)   = ', xtab(i)
          write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
          stop
       end if
    end do

    if ( a == b ) then
       return
    end if
    !
    !  If B < A, temporarily switch A and B, and store sign.
    !
    if ( b < a ) then
       syl = b
       b = a
       a = syl
       ind = -1
    else
       syl = a
       ind = 1
    end if
    !
    !  Find ILO and IHI so that A <= XTAB(ILO) <= XTAB(IHI) <= B
    !  with the possible exception that A and B may be in the same
    !  interval, or completely to the right or left of the XTAB's.
    !
    ilo = ntab + 1

    do i = 1, ntab
       if ( a <= xtab(i) ) then
          ilo = i
          exit
       end if
    end do

    ihi = 0

    do i = ntab, 1, -1
       if ( xtab(i) <= b ) then
          ihi = i
          exit
       end if
    end do
    !
    !  Treat special cases where A, B lie both to left or both to right
    !  of XTAB interval, or in between same pair of XTAB's.
    !
    if ( ihi == 0 ) then

       slope = ( ftab(2) - ftab(1) ) / ( xtab(2) - xtab(1) )
       fa = ftab(1) + slope * ( a - xtab(1) )
       fb = ftab(1) + slope * ( b - xtab(1) )
       result = 0.5_dp * ( b - a ) * ( fa + fb )

    else if ( ilo == ntab + 1 ) then

       slope = ( ftab(ntab) - ftab(ntab-1) ) / ( xtab(ntab) - xtab(ntab-1) )
       fa = ftab(ntab-1) + slope * ( a - xtab(ntab-1) )
       fb = ftab(ntab-1) + slope * ( b - xtab(ntab-1) )
       result = 0.5_dp * ( b - a ) * ( fa + fb )

    else if ( ihi + 1 == ilo ) then

       slope = ( ftab(ilo) - ftab(ihi) ) / ( xtab(ilo) - xtab(ihi) )
       fa = ftab(ihi) + slope * ( a - xtab(ihi) )
       fb = ftab(ihi) + slope * ( b - xtab(ihi) )
       result = 0.5_dp * ( b - a ) * ( fa + fb )

    else
       !
       !  Carry out approximate integration.  We know that ILO is no greater
       !  than IHI-1, but equality is possible; A and B may be on either side
       !  of a single XTAB(I).  That's OK, then the loop below won't be executed
       !  at all.
       !
       result = 0.0_dp
       do i = ilo, ihi-1
          result = result + 0.5_dp * ( xtab(i+1) - xtab(i) ) &
               * ( ftab(i) + ftab(i+1) )
       end do
       !
       !  Add contribution from A-ILO and IHI-B.
       !  Still have to watch out if ILO = 1 or IHI=NTAB...
       !
       if ( ilo == 1 ) then
          slope = ( ftab(2) - ftab(1) ) / ( xtab(2) - xtab(1) )
          fa = ftab(1) + slope * ( a - xtab(1) )
          result = result + 0.5_dp * ( xtab(ilo) - a ) * ( fa + ftab(ilo) )
       else
          slope = ( ftab(ilo) - ftab(ilo-1) ) / ( xtab(ilo) - xtab(ilo-1) )
          fa = ftab(ilo-1) + slope * ( a - xtab(ilo-1) )
          result = result + 0.5_dp * ( xtab(ilo) - a ) * ( fa + ftab(ilo) )
       end if

       if ( ihi == ntab ) then
          slope = ( ftab(ntab) - ftab(ntab-1) ) / ( xtab(ntab) - xtab(ntab-1) )
          fb = ftab(ntab-1) + slope * ( b - xtab(ntab-1) )
          result = result + 0.5_dp * ( b - xtab(ntab) ) * ( fb + ftab(ntab) )
       else
          slope = ( ftab(ihi+1) - ftab(ihi) ) / ( xtab(ihi+1) - xtab(ihi) )
          fb = ftab(ihi) + slope * ( b - xtab(ihi) )
          result = result + 0.5_dp * ( b - xtab(ihi) ) * ( fb + ftab(ihi) )
       end if

    end if
    !
    !  Restore original values of A and B, reverse sign of integral
    !  because of earlier switch.
    ! 
    if ( ind /= 1 ) then
       ind = 1
       syl = b
       b = a
       a = syl
       result = -result
    end if

    return

  end subroutine plint

  !----------------------------------------------------------------

  FUNCTION TRIMCAT(ch1, ch2)
    ! concatenate character strings ch1 & ch2
    IMPLICIT NONE
    CHARACTER (LEN=max_len)        :: TRIMCAT
    CHARACTER (LEN=*), INTENT (IN) :: ch1
    CHARACTER (LEN=*), INTENT (IN) :: ch2
    INTEGER                        :: i1, i2
    CHARACTER (LEN=max_len)        :: lch1, lch2, ch3
    lch1 = ADJUSTL(ch1)
    lch2 = ADJUSTL(ch2)
    i1 = LEN_TRIM(lch1)
    i2 = LEN_TRIM(lch2)
    ch3(1:i1) = TRIM(lch1)
    ch3(i1+1:i1+i2) = TRIM(lch2)
    ch3(i1+i2+1:max_len) = REPEAT(' ',max_len-i1-i2)
    TRIMCAT = TRIM(ch3)
  END FUNCTION TRIMCAT

  !----------------------------------------------------------------

  SUBROUTINE splint(xa, ya, y2a, n, x, y)

    ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
    ! (with the xa(i) in order), and given the array y2a(1:n), which is the output
    ! from the subroutine spline, and given a value of x, this routine returns a
    ! cubic spline interpolated value y.
    ! (adopted from Numerical Recipes in FORTRAN 77)
    
    REAL(DP), INTENT(IN)  :: xa(n)
    REAL(DP), INTENT(IN)  :: ya(n)
    REAL(DP), INTENT(IN)  :: y2a(n)
    INTEGER,  INTENT(IN)  :: n
    REAL(DP), INTENT(IN)  :: x
    REAL(DP), INTENT(OUT) :: y

    INTEGER :: k, khi, klo
    REAL(DP) :: a, b, h

    klo=1
    khi=n
1   if (khi-klo.gt.1) then
       k=(khi+klo)/2
       if (xa(k).gt.x) then
          khi=k
       else
          klo=k
       endif
       goto 1
    endif

    h=xa(khi)-xa(klo)
    !if (h.eq.0.) pause 'bad xa input in splint'
    if (h.eq.0.0_dp) pause 'bad xa input in splint'

    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    !y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
     y=a*ya(klo)+b*ya(khi)+((a**3.0_dp-a)*y2a(klo)+(b**3.0_dp-b)*y2a(khi))*(h**2.0_dp)/6.0_dp

    return
  END SUBROUTINE splint

  !----------------------------------------------------------------
  
  SUBROUTINE spline(x, y, n, yp1, ypn, y2)

    ! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
    ! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
    ! the first derivative of the interpolating function at points 1 and n,
    ! respectively, this routine returns an array y2(1:n) of length n which
    ! contains the second derivatives of the interpolating function at the
    ! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
    ! the routine is signaled to set the corresponding boundary condition for a
    ! natural spline with zero second derivative on that boundary.
    ! Parameter: nmax is the largest anticipiated value of n
    ! (adopted from Numerical Recipes in FORTRAN 77)
    
    REAL(DP), INTENT(IN)  :: x(n)
    REAL(DP), INTENT(IN)  :: y(n)
    INTEGER, INTENT(IN)   :: n
    REAL(DP), INTENT(IN)  :: yp1
    REAL(DP), INTENT(IN)  :: ypn
    REAL(DP), INTENT(OUT) :: y2(n)
    
    !INTEGER, PARAMETER:: nmax=500
     INTEGER, PARAMETER:: nmax=5000
    INTEGER:: i, k
    REAL(DP):: p, qn, sig, un, u(nmax)

    !if (yp1.gt..99e30) then
    !   y2(1)=0.
    !   u(1)=0.
    !else
    !   y2(1)=-0.5
    !   u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    !endif

    if (yp1.gt..99d30) then
       y2(1) = 0.0_dp
       u(1)  = 0.0_dp
    else
       y2(1)=-0.5_dp
       u(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif

    do i=2, n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      ! p=sig*y2(i-1)+2.
      ! y2(i)=(sig-1.)/p
      ! u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
      !      & (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      p=sig*y2(i-1)+2.0_dp
       y2(i)=(sig-1.0_dp)/p
       u(i)=(6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
            & (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p

    enddo

!    if (ypn.gt..99e30) then
!       qn=0.
!       un=0.
!    else
!       qn=0.5
!       un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
!    endif

!    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    if (ypn.gt..99d30) then
       qn=0.0_dp
       un=0.0_dp
    else
       qn=0.5_dp
       un=(3.0_dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif

    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0_dp)

    do k=n-1, 1, -1
       y2(k)=y2(k)*y2(k+1)+u(k)
    enddo

    return
  END SUBROUTINE spline

END MODULE MUTILITY
