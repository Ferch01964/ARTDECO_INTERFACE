
SUBROUTINE trunc_DM(ngauss, u, wth, F, nbetal,&
     betal, coeftr)

  implicit none

  !--- Input Arguments
  integer, intent(in)                                 :: ngauss
  integer, intent(in)                                 :: nbetal
  real(8), dimension(ngauss), intent(in)	      :: u
  real(8), dimension(ngauss), intent(in)	      :: wth
  real(8), dimension(6,ngauss), intent(in)            :: F

  !--- Local Arguments
  integer :: k, i, flag
  real(8), dimension(4,0:nbetal,ngauss)   :: apl
  real(8), dimension(ngauss)              :: F_trunc_recomp
  real(8), dimension(6,ngauss)            :: F_trunc
  
  !--- Output Arguments
  real(8), intent(out)                            :: coeftr
  real(8), dimension(6,0:nbetal), intent(out)	  :: betal

  !**************************************************************************
  !                              DELTA-M
  !**************************************************************************
  !	  Compute the moments Bl (legendre) for F11 only
  ! *********************************************************************

  flag = 0
  CALL comp_Bl(flag, ngauss, u, F, wth, nbetal, apl, betal)

  ! We renormalize the Betal so that B0 of the phase function is 1
    betal(1,:) = betal(1,:) / betal(1,0)  

  ! apply D-M
  coeftr     = betal(1,nbetal)
  betal(1,:) = (betal(1,:) - coeftr) / (1.D0 - coeftr)

  !--- Now recompute F11 from the betals to recompute Fij accordingly
  F_trunc_recomp(:) = 0.0D0
  DO k = 0, nbetal
     DO i = 1, ngauss
        F_trunc_recomp(i) = F_trunc_recomp(i)  + DBLE(2*k+1)*betal(1,k)*apl(1,k,i) ! F11
     ENDDO
  ENDDO
 
  F_trunc(:,:) = 0.0D0
  !--- rescale Fij, Fii
  F_trunc(1,1:ngauss) = F_trunc_recomp(1:ngauss)
  F_trunc(2,1:ngauss) = F(2,:) * F_trunc_recomp(1:ngauss) / F(1,:)
  F_trunc(3,1:ngauss) = F(3,:) * F_trunc_recomp(1:ngauss) / F(1,:)
  F_trunc(4,1:ngauss) = F(4,:) * F_trunc_recomp(1:ngauss) / F(1,:) 
  !F_trunc(5,1:ngauss) = F(5,:) * F_trunc_recomp(1:ngauss) / F(1,:) 
  !F_trunc(6,1:ngauss) = F(6,:) * F_trunc_recomp(1:ngauss) / F(1,:) 
  F_trunc(5,1:ngauss) = F(5,:) / (1.0D0 - coeftr)
  F_trunc(6,1:ngauss) = F(6,:) / (1.0D0 - coeftr)

  ! rescaled phase matrix terms does not have a strong forward peak anymore
  ! so that the Betal expansion is working well (no renormalisation to do)

  flag = 1
  CALL comp_Bl(flag, ngauss, u, F_trunc, wth, nbetal, apl, betal)

  ! *********************************************************************
  !          Convert beta moment in legendre coefficient 
  !          (multiply by factor 2k+1)
  ! *********************************************************************
  DO k = 0,nbetal
     betal(:,k) = DBLE(2*k+1)*betal(:,k)
  ENDDO

  RETURN

END SUBROUTINE trunc_DM
