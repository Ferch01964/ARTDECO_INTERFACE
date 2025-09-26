
SUBROUTINE comp_Bl(flag, nangle, u, F, w, nbetal, apl, betal)

  !  Calculate the expansion coefficients of the scattering matrix in    
  !  generalized spherical functions by numerical integration over the   
  !  scattering angle.                                                   
  !  From the adding paper 
  !  de Haan, Bosma, Hovenier, 1987, A&A, 183, 371     

  implicit none

  INTEGER, INTENT(IN) :: flag
  INTEGER, INTENT(IN) :: nbetal
  INTEGER, INTENT(IN) :: nangle
  REAL(8), INTENT(IN) :: u(nangle)
  REAL(8), INTENT(IN) :: w(nangle)
  REAL(8), INTENT(IN) :: F(6,nangle)
  !                            F(1) = F11
  !                            F(2) = F22
  !                            F(3) = F33
  !                            F(4) = F44
  !                            F(5) = F12
  !                            F(6) = F34

  REAL(8), INTENT(INOUT) :: apl(4,0:nbetal,nangle) 
  REAL(8), INTENT(OUT)   :: betal(6,0:nbetal)
  !                             betal(1) = alpha1
  !                             betal(2) = alpha2
  !                             betal(3) = alpha3
  !                             betal(4) = alpha4
  !                             betal(5) = beta1
  !                             betal(6) = beta2

  ! local variables
  REAL(8) :: P00(nangle,2) 
  REAL(8) :: P22(nangle,2)
  REAL(8) :: P2m2(nangle,2)
  REAL(8) :: P02(nangle,2)
  REAL(8) :: qroot6
  REAL(8) :: fac1
  REAL(8) :: fac2
  REAL(8) :: fac3
  REAL(8) :: sql41
  REAL(8) :: sql4
  REAL(8) :: twol1
  REAL(8) :: tmp2
  REAL(8) :: tmp1
  REAL(8) :: denom
  REAL(8) :: alfap
  REAL(8) :: alfam

  INTEGER :: i, l
  INTEGER :: lnew
  INTEGER :: lold
  INTEGER :: itmp

  ! ----------------
  ! compute Legendre Polynomials

  if (flag .eq. 0) then

     apl(:,:,:) = -32768.0D0

     qroot6     = -0.25D0*dsqrt(6.D0)

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

        apl(1,l,:) = P00(:,lnew)
        apl(2,l,:) = P02(:,lnew)
        apl(3,l,:) = P22(:,lnew)
        apl(4,l,:) = P2m2(:,lnew)

     ENDDO

  endif

  !--------------
  ! compute Betal

  if (flag .eq. 0) then

     betal(:,:)    = 0.0D0

     do l= 0, nbetal
        do i=1, nangle
           betal(1,l) = betal(1,l) + (apl(1,l,i) * w(i) * F(1,i) * 0.5D0) 
        enddo
     enddo ! End of loop over index l                                        

  else

     do l= 0, nbetal
        !  Now calculate the coefficients by integration over angle     *
        !  See de Haan et al. (1987) Eqs. (68)-(73).                    
        alfap = 0.0D0
        alfam = 0.0D0
        do i=1, nangle
           alfap = alfap           + apl(3,l,i) * w(i) * (F(2,i)+F(3,i)) * 0.5D0
           alfam = alfam           + apl(4,l,i) * w(i) * (F(2,i)-F(3,i)) * 0.5D0
           betal(4,l) = betal(4,l) + apl(1,l,i) * w(i) * F(4,i) *0.5D0
           betal(5,l) = betal(5,l) + apl(2,l,i) * w(i) * F(5,i) *0.5D0
           betal(6,l) = betal(6,l) + apl(2,l,i) * w(i) * F(6,i) *0.5D0
        enddo
        betal(2,l) = ( alfap + alfam ) / 2.0D0
        betal(3,l) = ( alfap - alfam ) / 2.0D0
     enddo ! End of loop over index l                                        

  endif

end SUBROUTINE comp_Bl
