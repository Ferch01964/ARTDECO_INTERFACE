

SUBROUTINE trunc_potter(ngauss, mu1, mu2, u, wth, F, nbetal,&
     betal, coeftr)
 
  implicit none

  !--- Input Arguments
  integer, intent(in) :: ngauss
  integer, intent(in) :: nbetal
  real(8), intent(in) :: mu1
  real(8), intent(in) :: mu2
  real(8), dimension(ngauss), intent(in) :: u
  real(8), dimension(ngauss), intent(in) :: wth
  real(8), dimension(6,ngauss), intent(in) :: F

  !--- Local Arguments
  real(8), dimension(6,ngauss)           :: F_trunc

  !--- Output Arguments
  real(8), intent(out) :: coeftr
  real(8), dimension(6,0:nbetal), intent(out) :: betal

  !**************************************************************************
  !                         POTTER TRUNCATION
  !     Now truncate scattering matrix in F_trunc (renormalize it also) 
  !     and compute truncation coefficient (coeftr).
  !**************************************************************************
  CALL truncation(ngauss,u,wth,F,mu1,mu2,F_trunc(:,:),coeftr)

  !**************************************************************************
  !     Compute Legendre polynomiale coefficient
  !**************************************************************************
  CALL devel_bl(nbetal, ngauss, u, wth, F_trunc(:,:),betal(:,:))

  RETURN

END SUBROUTINE trunc_potter

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

SUBROUTINE truncation(ngauss,u,wth,F,mu1,mu2,F_trunc,coeftr)

  implicit none

  !--- Input Arguments
  integer                                 :: ngauss

  real(8)                                 :: mu1, mu2
  real(8), dimension(ngauss)              :: u, wth
  real(8), dimension(6,ngauss)            :: F

  !--- Local Arguments
  logical                                 :: test

  integer                                 :: j, k1, k2

  real(8)                                 :: a, x1, y1, y, air_trunc

  !--- Output Arguments
  real(8)                                 :: coeftr
  real(8), dimension(6,ngauss)           :: F_trunc

  F_trunc(:,1:ngauss) = F(:,:)
  !**************************************************************************
  !    We need mu1 > mu2, so test for this constraint
  !**************************************************************************
  IF (mu2 .GE. mu1) THEN
     WRITE(6,*)'bad value for mu1 and mu2, mu1 should be > than mu2!'
     WRITE(6,*)'mu1 =',mu1
     WRITE(6,*)'mu2 =',mu2
     STOP
  ENDIF

  !**************************************************************************
  !   Seek for the value of u(j) (finally seek for j!) which is just before
  !   mu1 and mu2
  !**************************************************************************
    test = .true.
    DO 301 j = ngauss/2,ngauss
      IF (.not. (test .and. (u(j) .gt. mu2))) GOTO 301
      k2 = j-1
      test=.false.
301 CONTINUE

    test = .true.
    DO 302 j = ngauss/2,ngauss
      IF (.not. (test .and. (u(j) .gt. mu1))) GOTO 302
      k1 = j-1
      test=.false.
302 CONTINUE

  !<TMP>
  !    WRITE(6,*)'k1=',k1,' u(k1)=',u(k1)
  !    WRITE(6,*)'k2=',k2,' u(k2)=',u(k2)
  !</TMP>

  !**************************************************************************
  !   Compute the slope a and the troncated phase function between 
  !   the first gauss point (forward scattering) and the one just before 
  !   mu1
  !**************************************************************************
  a = (DLOG10(F(1,k1)) - DLOG10(F(1,k2))) / (DACOS(u(k1)) - DACOS(u(k2)))
  y1 = DLOG10(F(1,k1))
  x1 = DACOS(u(k1))
  !<TMP>
  !    WRITE(6,*)'a=',a
  !    WRITE(6,*)'y1=',y1
  !    WRITE(6,*)'x1=',x1
  !    y = y1 + a*(0.0D0 - x1)
  !    WRITE(6,*)'Ftrunc(0)=',10.0D0**y
  !</TMP>
  DO j = k1+1,ngauss
     y = y1 + a*(DACOS(u(j)) - x1)
     F_trunc(1,j) = 10.0D0**y
  ENDDO !j

  !**************************************************************************
  !   Compute truncation coefficient
  !**************************************************************************
  air_trunc = 0.0D0
  DO j = 1,ngauss
     air_trunc = air_trunc + F_trunc(1,j)*wth(j)
  ENDDO !j
  coeftr = (2.0D0 - air_trunc) / 2.0D0 !c'est reellement f et non 2f!!!!
  !WRITE(6,*)'air_trunc=',air_trunc

  !**************************************************************************
  !    Normalize the truncated phase function, and apply the truncation to 
  !    the other coefficient of the phase matrix
  !**************************************************************************
  !CALL normal(ngauss,u,wth,F_trunc(1,:))
  F_trunc(1,:) = F_trunc(1,:) / (1.0D0 - coeftr)

  F_trunc(2,1:ngauss) = F(2,1:ngauss) * F_trunc(1,1:ngauss) / F(1,1:ngauss)
  F_trunc(3,1:ngauss) = F(3,1:ngauss) * F_trunc(1,1:ngauss) / F(1,1:ngauss)
  F_trunc(4,1:ngauss) = F(4,1:ngauss) * F_trunc(1,1:ngauss) / F(1,1:ngauss)
  !F_trunc(5,1:ngauss) = F(5,1:ngauss) * F_trunc(1,1:ngauss) / F(1,1:ngauss)
  !F_trunc(6,1:ngauss) = F(6,1:ngauss) * F_trunc(1,1:ngauss) / F(1,1:ngauss)
  F_trunc(5,1:ngauss) = F(5,1:ngauss) / (1.0D0 - coeftr)
  F_trunc(6,1:ngauss) = F(6,1:ngauss) / (1.0D0 - coeftr)

  RETURN

END SUBROUTINE truncation

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

SUBROUTINE normal(ngauss,u,wth,F)

  implicit none

  !--- Input Arguments
  integer                                 :: ngauss

  real(8), dimension(ngauss)	            :: u, wth

  !--- Local Arguments
  integer                                 :: j

  real(8)                                 :: som, fact

  !--- inoutput Arguments
  real(8), dimension(ngauss)	            :: F

  som = 0.0D0
  DO j = 1,ngauss
     som = som + F(j)*wth(j)
  ENDDO !j
  fact = 2.0D0 / som
  F(:) = fact * F(:)

  RETURN

END SUBROUTINE normal

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

subroutine devel_Bl(nbetal, nangle, u, w, F, betal)

  !  Calculate the expansion coefficients of the scattering matrix in    
  !  generalized spherical functions by numerical integration over the   
  !  scattering angle.                                                   
  !  From the adding paper 
  !  de Haan, Bosma, Hovenier, 1987, A&A, 183, 371     

  implicit none

  INTEGER, INTENT(IN) :: nbetal
  INTEGER, INTENT(IN) :: nangle
  REAL(KIND=8), INTENT(IN) :: u(nangle)
  REAL(KIND=8), INTENT(IN) :: w(nangle)
  REAL(KIND=8), INTENT(IN) :: F(6,nangle)
  !                            F(1) = F11
  !                            F(2) = F22
  !                            F(3) = F33
  !                            F(4) = F44
  !                            F(5) = F12
  !                            F(6) = F34

  REAL(KIND=8), INTENT(OUT) :: betal(6,0:nbetal)
  !                             betal(1) = alpha1
  !                             betal(2) = alpha2
  !                             betal(3) = alpha3
  !                             betal(4) = alpha4
  !                             betal(5) = beta1
  !                             betal(6) = beta2

  ! local variables
  REAL(KIND=8) :: P00(nangle,2) 
  REAL(KIND=8) :: P22(nangle,2)
  REAL(KIND=8) :: P2m2(nangle,2)
  REAL(KIND=8) :: P02(nangle,2)
  REAL(KIND=8) :: qroot6
  REAL(KIND=8) :: fac1
  REAL(KIND=8) :: fac2
  REAL(KIND=8) :: fac3
  REAL(KIND=8) :: sql41
  REAL(KIND=8) :: sql4
  REAL(KIND=8) :: twol1
  REAL(KIND=8) :: tmp2
  REAL(KIND=8) :: tmp1
  REAL(KIND=8) :: denom
  REAL(KIND=8) :: alfap
  REAL(KIND=8) :: alfam
  REAL(KIND=8) :: fl

  INTEGER :: i, l
  INTEGER :: lnew
  INTEGER :: lold
  INTEGER :: itmp

  ! ----------------
  !  Initialization 

  qroot6        = -0.25D0*dsqrt(6.D0)
  betal(:,:)    = 0.0D0

  !  Start loop over the coefficient index l                             *
  !  first update generalized spherical functions, then calculate betal. *
  !  lold and lnew are pointer-like indices used in recurrence           *
  lnew = 1
  lold = 2

  do l= 0, nbetal

     if (l .eq. 0) then
        ! Adding paper Eq. (77) with m=0                           
        do i=1, nangle
           P00(i,lold) = 1.D0
           P00(i,lnew) = 0.D0
           P02(i,lold) = 0.D0
           P22(i,lold) = 0.D0
           P2m2(i,lold)= 0.D0
           P02(i,lnew) = 0.D0
           P22(i,lnew) = 0.D0
           P2m2(i,lnew)= 0.D0
        enddo
     else
        fac1 = (2.D0*l-1.D0)/dble(l)
        fac2 = dble(l-1.D0)/dble(l)
        ! Adding paper Eq. (81) with m=0                           *
        do i=1, nangle
           P00(i,lold) = fac1*u(i)*P00(i,lnew) - fac2*P00(i,lold)
        enddo
     endif

     if (l .eq. 2) then
        ! Adding paper Eqs. (78) and (80)  with m=2                *
        ! sql4 contains the factor dsqrt(l*l-4) needed in          *
        ! the recurrence Eqs. (81) and (82)                        *
        do i=1, nangle
           P02(i,lold) = qroot6*(1.D0-u(i)*u(i))
           P22(i,lold) = 0.25D0*(1.D0+u(i))*(1.D0+u(i))
           P2m2(i,lold)= 0.25D0*(1.D0-u(i))*(1.D0-u(i))
           P02(i,lnew) = 0.D0
           P22(i,lnew) = 0.D0
           P2m2(i,lnew)= 0.D0
        enddo
        sql41 = 0.D0
     endif

     if (l .gt. 2) then
        ! Adding paper Eq. (82) with m=0 and m=2                   *
        sql4  = sql41
        sql41 = dsqrt(dble(l*l)-4.D0)
        twol1 = 2.D0*dble(l)-1.D0
        tmp1  = twol1/sql41
        tmp2  = sql4/sql41
        denom = (dble(l)-1.D0)*(dble(l*l)-4.D0)
        fac1  = twol1*(dble(l)-1.D0)*dble(l)/denom
        fac2  = 4.D0*twol1/denom
        fac3  = dble(l)*((dble(l)-1.D0)*(dble(l)-1.D0)-4.D0)/denom
        do i=1, nangle
           P02(i,lold) = tmp1*u(i)*P02(i,lnew) - tmp2*P02(i,lold)
           P22(i,lold) = (fac1*u(i)-fac2)*P22(i,lnew) - fac3*P22(i,lold)
           P2m2(i,lold)= (fac1*u(i)+fac2)*P2m2(i,lnew) - fac3*P2m2(i,lold)
        enddo
     endif

     !  Switch indices so that lnew indicates the function with      *
     !  the present index value l, this mechanism prevents swapping  *
     !  of entire arrays.                                            *
     itmp = lnew
     lnew = lold
     lold = itmp

     !  Now calculate the coefficients by integration over angle     *
     !  See de Haan et al. (1987) Eqs. (68)-(73).                    
     alfap = 0.0D0
     alfam = 0.0D0
     do i=1, nangle
        betal(1,l) = betal(1,l) + P00(i,lnew) * w(i) * F(1,i)
        alfap = alfap           + P22(i,lnew) * w(i) * (F(2,i)+F(3,i))
        alfam = alfam           + P2m2(i,lnew)* w(i) * (F(2,i)-F(3,i))
        betal(4,l) = betal(4,l) + P00(i,lnew) * w(i) * F(4,i)
        betal(5,l) = betal(5,l) + P02(i,lnew) * w(i) * F(5,i)
        betal(6,l) = betal(6,l) + P02(i,lnew) * w(i) * F(6,i)
     enddo

     ! Multiply with trivial factors like 0.5D0*(2*l+1)             
     fl = dble(l)+0.5D0
     betal(1,l)   =  fl*betal(1,l)
     alfap        =  fl*alfap
     alfam        =  fl*alfam
     betal(4,l)   =  fl*betal(4,l)
     betal(5,l)   =  fl*betal(5,l)
     betal(6,l)   =  fl*betal(6,l)
     betal(2,l)   =  0.5D0*(alfap+alfam)
     betal(3,l)   =  0.5D0*(alfap-alfam)

  enddo ! End of loop over index l                                        

end subroutine devel_Bl
