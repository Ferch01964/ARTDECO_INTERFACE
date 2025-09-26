
SUBROUTINE trunc_delta_fit(nmug, u_init, w_init, F_init, flag_DM, &
     thetac, fitall, nbetal, betal, coeftr)

  use mtrunc_constant, only: D2R

  implicit none

  !--- Input Argument
  integer, INTENT(IN)                     :: nmug
  integer, INTENT(IN)                     :: flag_DM
  integer, INTENT(IN)                     :: nbetal
  real(8), dimension(nmug), INTENT(IN)    :: u_init
  real(8), dimension(nmug), INTENT(IN)    :: w_init
  real(8), dimension(6,nmug), INTENT(IN)  :: F_init
  !                                          F_init(1) = F11
  !                                          F_init(2) = F22
  !                                          F_init(3) = F33
  !                                          F_init(4) = F44
  !                                          F_init(5) = F12
  !                                          F_init(6) = F34
  real(8), INTENT(IN)                     :: thetac
  logical, INTENT(IN)                     :: fitall

  !--- Local Argument
  integer                                 :: l
  real(8), dimension(6,0:nbetal)          :: beta_hu
  real(8), dimension(6,nmug)              :: F
  real(8), dimension(nmug)                :: u
  real(8), dimension(nmug)                :: w

  !--- Output Argument
  real(8), INTENT(OUT)                        :: coeftr
  real(8), dimension(6,0:nbetal), INTENT(OUT) :: betal
  !                                              betal(1) = alpha1
  !                                              betal(2) = alpha2
  !                                              betal(3) = alpha3
  !                                              betal(4) = alpha4
  !                                              betal(5) = beta1
  !                                              betal(6) = beta2

  ! *********************************************************************
  !                 Compute the delta-fit coefficients
  ! *********************************************************************
  F(:,:) = F_init(:,:)
  u(:)   = u_init(:)
  w(:)   = w_init(:)

  !--- all elements of the phase matrix
  CALL b_fit(nmug, u, w, F, flag_DM,  thetac, fitall, nbetal, beta_hu,    &
       coeftr)

  ! *********************************************************************
  !          Convert beta moment in legendre coefficient 
  !          (multiply by factor 2l+1)
  ! *********************************************************************
  DO l = 0,nbetal
     betal(:,l) = DBLE(2*l+1)*beta_hu(:,l)
  ENDDO !l

  RETURN

END SUBROUTINE trunc_delta_fit

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

SUBROUTINE b_fit(nmug, u0, w0, F0, flag_DM,  thetac, fitall,&
     nbetal, Bfit, ftrunc)
  
  use mtrunc_constant, only: Pi, D2R

  implicit none

  !--- Input Arguments
  integer                     :: nmug, nbetal, flag_DM
  real(8), dimension(nmug)    :: u0, w0
  real(8), dimension(6,nmug)  :: F0
  real(8), INTENT(IN)                     :: thetac
  logical, INTENT(IN)                     :: fitall

  !--- Local Arguments
  integer :: i, k, flag, l

  real(8), dimension(nmug)               :: ang0
  real(8)                                :: Betmp
  real(8), dimension(6,0:nbetal)         :: Bl
  real(8), dimension(6,0:nbetal)         :: Bfitdm

  integer, parameter          :: nang = 361
  real(8), dimension(nang)    :: ang, u
  real(8), dimension(6,nang)  :: F
  real(8) :: step
  real(8) :: resint
  real(8) :: sum_recomp
  real(8) :: sum

  real(8), dimension(nmug)             ::  F_recomp0
  real(8), dimension(nang)             ::  F_recomp

  real(8), dimension(4,0:nbetal,nang)  :: apl
  real(8), dimension(4,0:nbetal,nmug)  :: apl0

  !--- Output Arguments
  real(8)   :: ftrunc
  real(8), dimension(6,0:nbetal)          :: Bfit

  !---------

  F(:,:)      = -32768.0D0
  apl(:,:,:)  = -32768.0D0
  apl0(:,:,:) = -32768.0D0

  !--- Compute scattering angle in degrees
  ang0(:) = acos(u0(:))/D2R

  ! We create a common regular grid :
  ! for the SVD fitting 
  step = ( maxval(ang0) - minval(ang0)) / (nang-1)
  DO i= 1, nang
     ang(i) = minval(ang0) + (dble(i-1)*step)
     u(i)   = COS(ang(i)*D2R)
     CALL INTPOL(LOG10(F0(1,:)), u0, nmug, u(i), resint)
     F(1,i) = 10.0D0**resint
  ENDDO
  if (fitall) then
     DO i= 1, nang
        CALL INTPOL(F0(2,:)/F0(1,:), u0, nmug, u(i), resint)
        F(2,i) = resint * F(1,i)
        CALL INTPOL(F0(3,:)/F0(1,:), u0, nmug, u(i), resint)
        F(3,i) = resint * F(1,i)
        CALL INTPOL(F0(4,:)/F0(1,:), u0, nmug, u(i), resint)
        F(4,i) = resint * F(1,i)
        CALL INTPOL(F0(5,:)/F0(1,:), u0, nmug, u(i), resint)
        F(5,i) = resint * F(1,i)
        CALL INTPOL(F0(6,:)/F0(1,:), u0, nmug, u(i), resint)
        F(6,i) = resint * F(1,i)
     ENDDO
  endif

  !######################################################################
  !                                  F11
  !   RQ: flag = 0 compute the delta fitting for F11 only
  !       flag = 1 compute it for other elements of the phase matrix
  !######################################################################
  flag = 0

  ! *********************************************************************
  !	  Compute the moments Bl (legendre) for F11 only
  !     F11     in Bl(1,:)   
  !   RQ: flag = 0 compute the delta fitting for F11 only
  !       flag = 1 compute it for other elements of the phase matrix
  ! *********************************************************************
  CALL comp_Bl(flag, nmug, u0, F0, w0, nbetal, apl0, Bl)
  DO i= 1, nang
     DO l= 0, nbetal
        CALL INTPOL(apl0(1, l, :), u0, nmug, u(i), resint)
        apl(1,l,i)=resint
     enddo
  ENDDO

  ! *********************************************************************
  !	  Compute the delta fitting for the phase function (F11) only
  !     in order to renormalize the other elements (cf Potter)
  !     F'ij = Fij * F'11/F11
  !     F11     in Bfitdm(1,:)
  !    
  !   RQ: flag = 0 compute the delta fitting for F11 only
  !       flag = 1 compute it for other elements of the phase matrix
  ! *********************************************************************
  CALL comp_dfit(flag, nang, ang, F, thetac, nbetal, apl, Bl, Bfitdm)

  !! Normalize
  Bfitdm(1, :) = Bfitdm(1, :) / Bfitdm(1, 0)

  if (flag_DM .eq. 1) then
     !--- Apply the delta-M scaling
     ftrunc    = Bfitdm(1,nbetal)
     Bfit(1,:) = (Bfitdm(1,:) - ftrunc) / (1.D0 - ftrunc)
  else
     write(*,*) 'Pbl!?'
     stop
     !ftrunc    = 1.0D0 - (Bfitdm(1,0))
     !Bfit(1,:) = Bfitdm(1,:) / (1.0D0 - ftrunc)
  endif

  !--- Now recompute F11 from the betals to recompute Fij accordingly
  F_recomp0(:) = 0.0D0
  F_recomp(:) = 0.0D0
  DO k = 0,nbetal
     DO i = 1,nmug
        F_recomp0(i) = F_recomp0(i) + DBLE(2*k+1)*Bfit(1,k)*apl0(1,k,i) ! F11
     ENDDO
  ENDDO
!!$  sum_recomp  = 0.0
!!$  sum = 0.0
!!$  DO i = 1,nmug
!!$     sum_recomp = sum_recomp + w0(i) *F_recomp0(i) 
!!$     sum        = sum + w0(i)*F0(1,i) 
!!$  ENDDO
!!$  write(*,*) sum_recomp
!!$  write(*,*) sum
  F0(2,:) = F0(2,:) * F_recomp0(:) / F0(1,:)
  F0(3,:) = F0(3,:) * F_recomp0(:) / F0(1,:)
  F0(4,:) = F0(4,:) * F_recomp0(:) / F0(1,:)
  !F0(5,:) = F0(5,:) * F_recomp0(:) / F0(1,:)
  !F0(6,:) = F0(6,:) * F_recomp0(:) / F0(1,:)
  F0(5,:) = F0(5,:) / (1.0D0 - ftrunc)
  F0(6,:) = F0(6,:) / (1.0D0 - ftrunc)
  F0(1,:) = F_recomp0(:)

  !######################################################################
  !
  !                                  Fij
  !
  !######################################################################
  flag = 1

  ! *********************************************************************
  !     Now compute the moments Bl (legendre) for all the elements of the 
  !     phase matrix but F11
  !     F22+F33 in Bl(2,:)
  !     F22-F33 in Bl(3,:)
  !     F44     in Bl(4,:)
  !     F12     in Bl(5,:)
  !     F34     in Bl(6,:)
  ! ********************************************************************
  CALL comp_Bl(flag, nmug, u0, F0, w0, nbetal, apl0, Bl, betmp)

  if (fitall) then

     DO i= 1, nang
        CALL INTPOL(LOG10(F_recomp0(:)), u0, nmug, u(i), resint)
        F_recomp(i) = 10.0D0**resint
     ENDDO
!!$  CALL myspline(nmug, ang0, LOG10(F_recomp0(:)), nang, ang, F_recomp(:))
!!$  F_recomp(:) = 10.0D0**F_recomp(:)
     F(2,:) = F(2,:) * F_recomp(:) / F(1,:)
     F(3,:) = F(3,:) * F_recomp(:) / F(1,:)
     F(4,:) = F(4,:) * F_recomp(:) / F(1,:)
     !F(5,:) = F(5,:) * F_recomp(:) / F(1,:)
     !F(6,:) = F(6,:) * F_recomp(:) / F(1,:)
     F(5,:) = F(5,:) / (1.0D0 - ftrunc)
     F(6,:) = F(6,:) / (1.0D0 - ftrunc)
     F(1,:) = F_recomp(:)

     DO i= 1, nang
        DO l= 0, nbetal
           DO k= 2, 4
              CALL INTPOL(apl0(k, l, :), u0, nmug, u(i), resint)
              apl(k,l,i)=resint
           enddo
        enddo
     ENDDO

     ! *********************************************************************
     !     Compute the delta fitting for all element but F11 of the phase
     !     matrix
     !     F22+F33 in Bfitdm(2,:)
     !     F22-F33 in Bfitdm(3,:)
     !     F44     in Bfitdm(4,:)
     !     F12     in Bfitdm(5,:)
     !     F34     in Bfitdm(6,:)
     ! *********************************************************************
     CALL comp_dfit(flag, nang, ang, F, thetac, nbetal, apl, Bl, Bfitdm)
     Bfit(2:6,:) = Bfitdm(2:6,:)

  else

     ! just use the Bl for regular devel (no fitting)
     Bfit(2:6,:) = Bl(2:6,:)

  endif

  RETURN

END SUBROUTINE b_fit

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

SUBROUTINE comp_dfit(flag,ngauss,ang,y,thetac,nbetal,apl,Bl,Bfitdm)

  IMPLICIT NONE

  !--- Input  Arguments
  integer                                   :: ngauss, nbetal, flag
  real(8), dimension(6,0:nbetal)            :: Bl
  real(8), dimension(6,ngauss)              :: y
  real(8)                                   :: thetac
  real(8), dimension(ngauss)                :: ang
  real(8), dimension(4,0:nbetal,ngauss)     :: apl

  !--- Local Arguments
  integer                                  :: ndata, ma, i, k, l, nkp
  real(8), dimension(ngauss)               :: wght 
  real(8), dimension(6,ngauss,nbetal+1)    :: u
  real(8), dimension(6,ngauss)             :: b
  real(8), dimension(6,nbetal+1)           :: a
  real(8), dimension(ngauss)               :: y2my3
  real(8), dimension(0:nbetal)             :: alfap
  real(8), dimension(0:nbetal)             :: alfam

  !--- Output Arguments
  real(8), dimension(6,0:nbetal)          :: Bfitdm

  ! ***********************************************************************
  !  keep the first few moments Bl(0:nkp)
  ! ***********************************************************************
  nkp = 0

  ! ***********************************************************************
  !  compute b and u (see svdfit for their meaning) for the fits
  ! ***********************************************************************
  ndata = ngauss
  ma = nbetal - nkp

  ! set-up the weight for each angle
  wght(:) = 1.0D0
  if (flag .eq. 0) then
     DO i = 1,ngauss
        if (ang(i) .lt. thetac) wght(i) = 0.0d0
     ENDDO
  endif
  !  if (flag .eq. 1) then
  !     DO i = 1,ngauss
  !        if (ang(i) .gt. 178.0) wght(i) = 0.0d0
  !     ENDDO
  !  endif

  IF (flag .EQ. 0) THEN !F11 only

     DO k = 1,ngauss
        b(1,k) = 1.D0
        DO l = 0,nkp
           b(1,k) = b(1,k) - DBLE(2*l+1)*Bl(1,l)*apl(1,l,k)/y(1,k) !F11
        ENDDO !l
        u(1,k,ma) = 0.D0
        DO l = 0,nbetal-1
           u(1,k,ma) = u(1,k,ma) - DBLE(2*l+1)*apl(1,l,k)/y(1,k)   !F11
        ENDDO !l
        DO l = nkp+1,nbetal-1
           u(1,k,l-nkp) = DBLE(2*l+1)*apl(1,l,k)/y(1,k)            !F11
        ENDDO !l
        b(1,k) = b(1,k) * wght(k)
     ENDDO !k

  ELSE                 !Fij

     alfap(:) = Bl(2,:) + Bl(3,:)
     alfam(:) = Bl(2,:) - Bl(3,:)

     DO k = 1,ngauss
        y2my3(k) = y(2,k)-y(3,k)
        b(2:6,k) = 1.D0
        DO l = 0,nkp
           b(2,k) = b(2,k) - DBLE(2*l+1)*alfap(l)*apl(3,l,k)/(y(2,k)+y(3,k)) ! F22+F33
           b(3,k) = b(3,k) - DBLE(2*l+1)*alfam(l)*apl(4,l,k)/(y2my3(k))      ! F22-F33
           b(4,k) = b(4,k) - DBLE(2*l+1)*Bl(4,l)*apl(1,l,k)/y(4,k)    !F44
           b(5,k) = b(5,k) - DBLE(2*l+1)*Bl(5,l)*apl(2,l,k)/y(5,k)    !F12
           b(6,k) = b(6,k) - DBLE(2*l+1)*Bl(6,l)*apl(2,l,k)/y(6,k)    !F34
        ENDDO !l
        u(2:6,k,ma) = 0.D0
        DO l = 0,nbetal-1
           u(2,k,ma) = u(2,k,ma) - DBLE(2*l+1)*apl(3,l,k)/(y(2,k)+y(3,k))    ! F22+F33
           u(3,k,ma) = u(3,k,ma) - DBLE(2*l+1)*apl(4,l,k)/(y2my3(k))         ! F22-F33
           u(4,k,ma) = u(4,k,ma) - DBLE(2*l+1)*apl(1,l,k)/y(4,k)       !F44
           u(5,k,ma) = u(5,k,ma) - DBLE(2*l+1)*apl(2,l,k)/y(5,k)       !F12
           u(6,k,ma) = u(6,k,ma) - DBLE(2*l+1)*apl(2,l,k)/y(6,k)       !F34
        ENDDO !l
        DO l = nkp+1,nbetal-1
           u(2,k,l-nkp) = DBLE(2*l+1)*apl(3,l,k)/(y(2,k)+y(3,k))   !F22+F33
           u(3,k,l-nkp) = DBLE(2*l+1)*apl(4,l,k)/(y2my3(k))        !F22-F33
           u(4,k,l-nkp) = DBLE(2*l+1)*apl(1,l,k)/y(4,k)            !F44
           u(5,k,l-nkp) = DBLE(2*l+1)*apl(2,l,k)/y(5,k)            !F12
           u(6,k,l-nkp) = DBLE(2*l+1)*apl(2,l,k)/y(6,k)            !F34
        ENDDO !l

        b(:,k) = b(:,k) * wght(k)

        ! check wheteher y(2,:) - y(3,:)
        ! is not too low ( 1/(y(2,:) - y(3,:)) --> inf. )
        if (y2my3(k) .lt. 1.d-100) then
           u(3,k,:)  = 0.0
           b(3,k)    = 0.0
        endif

     ENDDO !k

  ENDIF

  ! *************************************************************************
  !  singular value decomposition fitting to derive b
  ! *************************************************************************
  IF (flag .EQ. 0) THEN
     a(1,:) = 0.0D0
     CALL svdfit(nbetal,ndata,a(1,:),ma,u(1,:,:),b(1,:))       !F11
  ELSE
     a(2:6,:) = 0.0D0
     CALL svdfit(nbetal,ndata,a(2,:),ma,u(2,:,:),b(2,:))       !F22+F33
     CALL svdfit(nbetal,ndata,a(3,:),ma,u(3,:,:),b(3,:))       !F22-F33
     CALL svdfit(nbetal,ndata,a(4,:),ma,u(4,:,:),b(4,:))       !F44
     CALL svdfit(nbetal,ndata,a(5,:),ma,u(5,:,:),b(5,:))       !F12
     CALL svdfit(nbetal,ndata,a(6,:),ma,u(6,:,:),b(6,:))       !F34
  ENDIF

  ! *************************************************************************
  !  the moments which we want to keep
  ! *************************************************************************
  IF (flag .EQ. 0) THEN
     DO i = 0,nkp
        Bfitdm(1,i) = Bl(1,i)
     ENDDO !i
  ELSE
     DO i = 0,nkp
        Bfitdm(2:6,i) = Bl(2:6,i)
     ENDDO !i
  ENDIF

  ! *************************************************************************
  !  other moments (the delta-truncation needs to be considered)
  !  the fitted values are the ones after truncation factor
  !  So, the true moments should add the truncation factor
  ! *************************************************************************
  IF (flag .EQ. 0) THEN
     DO i = nkp+1,nbetal
        Bfitdm(1,i) = a(1,i-nkp)
     ENDDO !i
  ELSE
     DO i = nkp+1,nbetal
        Bfitdm(4:6,i) = a(4:6,i-nkp)
        Bfitdm(2,i)   = (a(2,i-nkp) + a(3,i-nkp)) / 2.0D0
        Bfitdm(3,i)   = (a(2,i-nkp) - a(3,i-nkp)) / 2.0D0
     ENDDO !i
  ENDIF

  RETURN

END SUBROUTINE comp_dfit

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

SUBROUTINE svdfit(nbetal,ndata,a,ma,u,b)

  ! this subroutine performs least square fitting:
  !
  !        bb_i = a(1) * u_i(1) + a(2) * u_i(2) + ... + a(N) * u_i(N)
  !                 (i=1,ndata) (N: ma)
  !        so that   Sum ( b - bb )^2 is minimum
  !
  ! input:
  !
  !         u_i(k); (k=1,ma; i=1,ndata)   the x values,
  !                                 2-dimensional variables
  !
  !         b_i;          ( i=1,ndata )  the y values
  !
  ! output:
  !
  !         a(i); i=1,ma
  !
  !--- Input Arguments
  INTEGER                                 :: ma, ndata, nbetal
  REAL(8), DIMENSION(ndata)               :: b
  REAL(8), DIMENSION(ndata,nbetal+1)      :: u

  !--- Local Arguments
  INTEGER                                 :: i,j

  REAL(8), PARAMETER                      :: TOL=1.D-9
  REAL(8)                                 :: thresh, wmax
  REAL(8), dimension(ndata)               :: bb, b0
  REAL(8), dimension(nbetal+1)            :: w
  REAL(8), dimension(nbetal+1,nbetal+1)   :: v
  REAL(8), dimension(ndata,nbetal+1)      :: u0

  !--- Output Arguments
  REAL(8), DIMENSION(nbetal+1)            :: a

  !--- Initialized some parameters
  a(:) = 0.0D0
  b0(:) = b(:)
  DO j = 1,ma
     DO i = 1,ndata
        u0(i,j) = u(i,j)
     ENDDO !i
  ENDDO !j

  ! *************************************************************************
  !    evaluate the least square problem and sorting the eigenvalues
  ! *************************************************************************
  CALL svdcmp(u,ndata,ma,ndata,nbetal+1,w,v)

  ! *************************************************************************
  !    finding the largest eigenvalues and the cutoff eigenvalues
  ! *************************************************************************
  wmax = 0.D0
  DO j = 1,ma
     if (w(j) .gt. wmax) wmax = w(j)
  ENDDO !j
  thresh = TOL*wmax

  ! *************************************************************************
  !    zero the eigenvalues smaller than the cutoff values
  ! *************************************************************************
  DO j = 1,ma
     if (w(j) .lt. thresh) w(j) = 0.D0
  ENDDO !j

  ! *************************************************************************
  !    solve for a
  ! *************************************************************************
  CALL svbksb(u,w,v,ndata,ma,ndata,nbetal+1,b,a)
  DO i = 1,ndata
     bb(i) = 0.D0
     DO j = 1,ma
        bb(i) = bb(i)+a(j)*u0(i,j)
     ENDDO !j
  ENDDO !i

  RETURN

END SUBROUTINE svdfit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTPOL(fint, xint, ni, xess, fess)
  ! linear interoplation of fint @ xess
  ! xint assumed increasing, NO EXTAPOLATION 
  ! make sure xess belongs to [xint(1),xint(ni)]

  IMPLICIT NONE

  INTEGER, INTENT (IN)  :: ni
  REAL (8), INTENT (IN) :: fint(ni)
  REAL (8), INTENT (IN) :: xint(ni)
  REAL (8), INTENT (IN) :: xess
  REAL (8), INTENT (OUT) :: fess
  INTEGER                :: i

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
     fess = fint(i) * (xint(i+1) - xess) + fint(i+1) * (xess-xint(i))
     fess = fess / (xint(i+1) - xint(i))
  ENDIF

END SUBROUTINE INTPOL
