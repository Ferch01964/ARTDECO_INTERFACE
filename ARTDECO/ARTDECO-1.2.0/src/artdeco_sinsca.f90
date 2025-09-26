
MODULE MSINSCA

  IMPLICIT NONE

  PUBLIC  :: CALL_SINSCA, SINSCA, PRAY

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CALL_SINSCA

    USE MCONSTANTS, only : dp, undef_dp
    USE MCOMMON, only : verbose, &
         nmat, &
         mode, &
         nsza, &
         sza, &
         nvza, &
         vza, &
         nvaa, &
         vaa, &
         nlambda, &
         wlambda, &
         depol, &
         nlayers, &
         nalt_atm, &
         Svin, &
         F0, &
         nptcle, & 
         ptcle_nang_phasemat_layers, &
         ptcle_u_phasemat_layers, &
         ptcle_phasemat_layers, &
         ptcle_dtau_layers, &
         ptcle_ssa_layers, &
         surface_albedo, &
         surface_type, &
         surface_family, &
         n_aitaui_total, &
         n_aitaui_mono_layers, &
         nmax_aitaui_mono_layers, &
         ai_mono_layers, &
         gas_dtau_mono_layers, &
         gas_ssa_mono_layers, &
         SvR, &
         SvR_sig, &
         flux_out, &
         rt_cputime, &
         thermal, &
         varsol_fact

    USE MUTILITY, only : REVERSE
    USE MLAYERS, only : GET_MONO_LAYERS

    IMPLICIT NONE

    INTEGER :: sinsca_nmu
    REAL(kind=dp), ALLOCATABLE :: sinsca_mu(:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_p11_ptcle(:,:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_21_ptcle(:,:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_tau_gas(:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_ssa_gas(:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_tau_ptcle(:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_ssa_ptcle(:)
    REAL(kind=dp) :: sinsca_Stoke_in(nmat)

    REAL(kind=dp) :: sinsca_surface_albedo

    REAL(kind=dp), ALLOCATABLE :: sinsca_Stoke(:,:,:)

    INTEGER :: i, j, k, u, iai

    ! time
    REAL :: tcpu1
    REAL :: tcpu2

    ! -------------

    IF (thermal) THEN
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (sub. call_sinsca)  ERROR            '
       WRITE(*,*) '  Thermal emission is not implemented '
       WRITE(*,*) '  in call_sinsca                      '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       STOP
    ENDIF

    IF (verbose) THEN
       WRITE (*,*) ''
       WRITE (*, *)                   '  (sub. call_sinsca) set arguments'
       WRITE (*,*) ''
       WRITE (*,FMT='(A,1x,I6,1x,A)') '                    SINSCA will be called ', n_aitaui_total, ' times'
       if (mode.eq.'kdis') &
            WRITE (*,FMT='(A,1x,I6,A)') '                  (maximum number of call per wavelength is ', nmax_aitaui_mono_layers,')'
    ENDIF

    ALLOCATE(sinsca_tau_gas(nlayers), sinsca_ssa_gas(nlayers),&
         sinsca_tau_ptcle(nlayers), sinsca_ssa_ptcle(nlayers), &
         sinsca_Stoke(nmat, nvza, nvaa))

    ! init in case no lambertian surface
    sinsca_surface_albedo = 0.0_dp

    DO j = 1, nlambda

       CALL GET_MONO_LAYERS(j)

       sinsca_nmu = ptcle_nang_phasemat_layers(j)

       ALLOCATE(sinsca_mu(sinsca_nmu), &
            sinsca_p11_ptcle(nlayers,sinsca_nmu), &
            sinsca_21_ptcle(nlayers,sinsca_nmu))

       if (surface_family .eq. 'lambert')  sinsca_surface_albedo = surface_albedo(j)

       ! we revesre mu for it to be sorted in increasing order
       sinsca_mu(:) = REVERSE(sinsca_nmu, ptcle_u_phasemat_layers(1:sinsca_nmu,j))

       if (nptcle .gt. 0) then
          do i = 1, nlayers
             ! we revert P11 and 21 cause mu was reverted to be 
             ! sorted in increasing order
             sinsca_p11_ptcle(i,:)  = REVERSE(sinsca_nmu, ptcle_phasemat_layers(i,1,1:sinsca_nmu,j)) ! P11
             sinsca_21_ptcle(i,:)   = REVERSE(sinsca_nmu, ptcle_phasemat_layers(i,5,1:sinsca_nmu,j)) ! P12 = P21
             sinsca_tau_ptcle(i)    = ptcle_dtau_layers(i,j) 
             sinsca_ssa_ptcle(i)    = ptcle_ssa_layers(i,j) 
          enddo
       else
          sinsca_p11_ptcle(:,:) = 0.0_dp
          sinsca_21_ptcle(:,:) = 0.0_dp
          sinsca_tau_ptcle(:)   = 0.0_dp
          sinsca_ssa_ptcle(:)   = 0.0_dp
       endif

       SvR(:,:,:,:,j)    = 0.0_dp
       flux_out(:,:,:,j) = 0.0_dp

       IF (verbose) then
          if (mode.eq.'kdis') then
             WRITE (*, FMT='(A,2x,F12.5,1x,A, 1x, I5)')  &
                  &'  (sub. call_sinsca) run SINSCA for wlambda=', wlambda(j),' microns with nai=', n_aitaui_mono_layers
          else
             WRITE (*, FMT='(A45,2x,F12.5,1x,A8)')  &
                  &'  (sub. call_sinsca) run SINSCA for wlambda=', wlambda(j),' microns'  
          endif
       END IF

       CALL CPU_TIME(tcpu1)

       do iai = 1, n_aitaui_mono_layers

          ! init result arrays
          sinsca_Stoke(:,:,:) = undef_dp 
          sinsca_tau_gas(:)   = gas_dtau_mono_layers(:,iai) 
          sinsca_ssa_gas(:)   = gas_ssa_mono_layers(:,iai)

          DO i = 1, nsza

             sinsca_Stoke_in(1:nmat) = Svin(1:nmat) * F0(j) * varsol_fact(i)

             CALL SINSCA(verbose, &
                  wlambda(j), &
                  nmat, &
                  sza(i), nvza, vza, nvaa, vaa, &
                  nlayers, &
                  sinsca_nmu, sinsca_mu, sinsca_p11_ptcle, sinsca_21_ptcle, &
                  sinsca_tau_ptcle, &
                  sinsca_ssa_ptcle, &
                  sinsca_tau_gas,   &
                  sinsca_ssa_gas,   &
                  depol(j), &
                  sinsca_Stoke_in, &
                  surface_family, &
                  sinsca_surface_albedo, &
                  sinsca_Stoke)

             ! implement the common variable for stokes vector for reflected intensities for each geometry
             DO u = 1, nvza
                DO k = 1, nvaa
                   SvR(i,u,k,1,j)    = SvR(i,u,k,1,j) + sinsca_Stoke(1,u,k) * ai_mono_layers(iai)
                   if (nmat .gt. 1)then
                      SvR(i,u,k,2,j) = SvR(i,u,k,2,j) + sinsca_Stoke(2,u,k) * ai_mono_layers(iai)
                      SvR(i,u,k,3,j) = SvR(i,u,k,3,j) + sinsca_Stoke(3,u,k) * ai_mono_layers(iai)
                   endif
                   if (nmat .gt. 3) then 
                      SvR(i,u,k,4,j) = SvR(i,u,k,4,j) + sinsca_Stoke(4,u,k) * ai_mono_layers(iai)
                   endif
                ENDDO
             ENDDO

          ENDDO ! on sza

       end do

       CALL CPU_TIME(tcpu2)
       rt_cputime(j) = tcpu2-tcpu1
       ! IF (verbose) THEN
       !    WRITE (*,FMT='(1x,A35,1PE14.5, A9)') '  (sub. call_sinsca) CPU time     ', tcpu2-tcpu1, ' sec(s) -'
       ! ENDIF

       DEALLOCATE(sinsca_mu, &
            sinsca_p11_ptcle, &
            sinsca_21_ptcle)

    ENDDO ! on lambda

    SvR_sig(:,:,:,:,:)  = 0.0_dp

    DEALLOCATE(sinsca_tau_gas, sinsca_ssa_gas,&
         sinsca_tau_ptcle, sinsca_ssa_ptcle, &
         sinsca_Stoke)

  END SUBROUTINE CALL_SINSCA

  !===================================================================================================

  SUBROUTINE SINSCA(verbose, &
       lambda, &
       nmat, &
       tetas, ntetav, tetav, nphiv, phiv, &
       nlay, &
       nmu, mu, p11_ptcle, p21_ptcle, &
       tau_ptcle, &
       ssa_ptcle, &
       tau_gas,   &
       ssa_gas,   &
       delta_ray, &
       I0, &
       surface, &
       surface_albedo, &
       SvR_res)

    !  Compute RT under single scattering assumption

    USE MCONSTANTS, only : dp, xpi, undef_dp, deg2rad
    USE MUTILITY, only : LININTPOL
    USE MSURFACE_BRDF

    IMPLICIT NONE

    LOGICAL, INTENT(IN)       :: verbose
    REAL(KIND=dp), INTENT(IN) :: lambda    
    !    LOGICAL, INTENT(IN)       :: thermal
    INTEGER, INTENT(IN)       :: nmat
    REAL(KIND=dp), INTENT(IN) :: tetas
    INTEGER, INTENT(IN)       :: ntetav
    REAL(KIND=dp), INTENT(IN) :: tetav(ntetav)
    INTEGER, INTENT(IN)       :: nphiv
    REAL(KIND=dp), INTENT(IN) :: phiv(nphiv)
    INTEGER, INTENT(IN)       :: nlay
    INTEGER, INTENT(IN)       :: nmu
    REAL(KIND=dp), INTENT(IN) :: mu(nmu)
    REAL(KIND=dp), INTENT(IN) :: p11_ptcle(nlay,nmu)
    REAL(KIND=dp), INTENT(IN) :: p21_ptcle(nlay,nmu)
    REAL(KIND=dp), INTENT(IN) :: tau_ptcle(nlay)
    REAL(KIND=dp), INTENT(IN) :: ssa_ptcle(nlay)
    REAL(KIND=dp), INTENT(IN) :: tau_gas(nlay)
    REAL(KIND=dp), INTENT(IN) :: ssa_gas(nlay)
    REAL(KIND=dp), INTENT(IN) :: delta_ray
    REAL(KIND=dp), INTENT(IN) :: I0(nmat)
    CHARACTER(len=*), INTENT(IN) :: surface
    REAL(KIND=dp), INTENT(IN)  :: surface_albedo
    REAL(KIND=dp), INTENT(OUT) :: SvR_res(nmat, ntetav, nphiv)

    ! local arguments
    INTEGER :: itetav
    INTEGER :: iphiv
    INTEGER :: ilay
    !INTEGER :: i
    INTEGER :: i, j

    REAL(kind=dp) :: umu0
    REAL(kind=dp) :: umu
    REAL(kind=dp) :: ctetad
    REAL(kind=dp) :: stetad
    REAL(kind=dp) :: i2, ci2

    REAL(KIND=dp) :: TAU( 0:nlay )
    REAL(KIND=dp) :: tau_ly(nlay)
    REAL(KIND=dp) :: tausca_ly(nlay)
    REAL(KIND=dp) :: ptcle_tausca_ly(nlay)
    REAL(KIND=dp) :: gas_tausca_ly(nlay)
    REAL(KIND=dp) :: ssa_ly(nlay)
    REAL(KIND=dp) :: p11_ly(nlay)
    REAL(KIND=dp) :: p21_ly(nlay)

    REAL(KIND=DP) :: tau_tot

    REAL(kind=dp) :: p11ray
    REAL(kind=dp) :: p21ray

    REAL(KIND=dp) ::  EXP0, EXP1
    REAL(KIND=dp) ::  I_SINSCA, Q_SINSCA
    REAL(KIND=dp) ::  I_SURF, Q_SURF, U_SURF
    REAL(kind=dp) ::  brdf(nmat,nmat)

    !-----------

    if ( nmat .gt. 1) then
       if ( (I0(2) .ne. 0.0_dp) .or. (I0(3) .ne. 0.0_dp) .or. (I0(nmat) .ne. 0.0_dp) )  then
          write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write(*,*) ' (sub. sinsca) :   ERROR   '
          WRITE(*,*) '     The -sinsca- RT model only apply for unpolarized '
          WRITE(*,*) '     incoming radiation at the top of the atmopshere (SvR={1,0,0,0}) '
          write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          STOP
       endif
    endif

    !do j = 1, nmu
    !   print*, mu(j)
    !end do

    ! mu must be sorted in increasing order
    do i = 1, nmu-1
       if (mu(i) .ge. mu(i+1)) then
          write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write(*,*) ' (sub. sinsca) :   ERROR'
          write(*,*) '   mu must be sorted in increasing order'
          write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          STOP
       endif
    enddo

    SvR_res(:,:,:) = 0.0_dp

    ! get the total optical depth of layers (ptcle+gas)
    tau( 0 ) = 0.0_dp
    tau_tot  = 0.0_dp
    do ilay = 1, nlay
       tau_ly(ilay)    = tau_ptcle(ilay) + tau_gas(ilay) 
       tausca_ly(ilay) = (tau_ptcle(ilay) * ssa_ptcle(ilay)) + &
            &            (tau_gas(ilay)   * ssa_gas(ilay)) 
       ssa_ly(ilay)    = tausca_ly(ilay) / tau_ly(ilay)
       ptcle_tausca_ly(ilay) = tau_ptcle(ilay) * ssa_ptcle(ilay)
       gas_tausca_ly(ilay)   = tau_gas(ilay)   * ssa_gas(ilay)
       tau(ilay)             = tau(ilay-1) + tau_ly(ilay) 
       tau_tot               = tau_tot + tau_ly(ilay)  
    end do

    ! start to loop over geometries
    umu0 = COS(tetas* deg2rad)
    DO itetav = 1, ntetav

       UMU = COS(tetav(itetav) * deg2rad)

       DO iphiv = 1, nphiv

          ! Get the cosine of scattering angle
          ctetad = -UMU0 * UMU + SQRT( 1.0_dp - UMU0**2.0_dp ) * &
               SQRT( 1.0_dp - UMU**2.0_dp )  * COS( phiv(iphiv)*deg2rad )
          if (nmat .gt. 1) then
             STETAD  = SQRT(1.0_dp - CTETAD**2.0_dp)              
             if (STETAD.eq.0.0_dp) then
                CI2  = 1.0_dp
             else
                CI2  = ( - UMU0*SQRT( 1.0_dp - UMU**2.0_dp ) &
                     - UMU*SQRT( 1.0_dp - UMU0**2.0_dp ) * COS( phiv(iphiv)*deg2rad ) ) / &
                     STETAD
                IF (ci2 .GT.  1.0D0) ci2 =  1.0D0
                IF (ci2 .LT. -1.0D0) ci2 = -1.0D0                
             endif
             !--- The sign (-) is set to make a clockwise direction rotation to go from 
             !    the scattering plan to the local meridian plan
             if (phiv(iphiv).le.180.0_dp.and.phiv(iphiv).ge.0.0_dp)   i2 = ACOS(ci2) 
             if (phiv(iphiv).gt.180.0_dp.and.phiv(iphiv).le.360.0_dp) i2 = -ACOS(ci2)
          else
             i2 = undef_dp
          endif

          ! get the phase matrix terms P11 and P21 for 
          ! each layer of the atmosphere (gas + particles)
          p11_ly(:) = 0.0_dp
          p21_ly(:) = 0.0_dp
          CALL pray(ctetad, stetad, delta_ray, p11ray, p21ray)
          do ilay = 1, nlay
             ! ptcle contribution
             if (ptcle_tausca_ly(ilay).gt.0.0_dp) then
                p11_ly(ilay) = ptcle_tausca_ly(ilay) / tausca_ly(ilay) *    &
                     LININTPOL(p11_ptcle(ilay,:), mu(:), nmu, ctetad)
                if (nmat.gt.1) then
                   p21_ly(ilay) = ptcle_tausca_ly(ilay) / tausca_ly(ilay) * &
                        LININTPOL(p21_ptcle(ilay,:), mu(:), nmu, ctetad)
                endif
             endif
             ! add rayleigh contribution
             if (gas_tausca_ly(ilay).gt.0.0_dp) then
                p11_ly(ilay) =  p11_ly(ilay) + &
                     gas_tausca_ly(ilay) / tausca_ly(ilay) * p11ray
                if (nmat.gt.1) then
                   p21_ly(ilay) =  p21_ly(ilay) + &
                        gas_tausca_ly(ilay) / tausca_ly(ilay) * p21ray
                endif
             end if
          end do

          I_SINSCA = 0.0_dp
          Q_SINSCA = 0.0_dp
          EXP0     = 1.0_dp
          IF(UMU.GT.0.0_dp) THEN
             DO ilay = 1, nlay
                ! ** Upward intensity, Eq. STWL (65b)
                EXP1 = EXP( -( ( TAU( ilay ) )/UMU + TAU( ilay )/UMU0 ) )
                I_SINSCA = I_SINSCA + SSA_LY( ilay )*P11_LY( ilay )*( EXP0 - EXP1 )
                if (nmat .gt. 1) Q_SINSCA = Q_SINSCA + SSA_LY( ilay )*P21_ly( ilay )*( EXP0 - EXP1 )
                EXP0 = EXP1
             ENDDO
          ELSE
             write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             write(*,*) ' (sub. sinsca) :   ERROR   '
             WRITE(*,*) '         Downward intensity to be implemented'
             write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             write(*,*) ' '
             STOP
          END IF

          ! surface contribution
          IF (surface.eq.'lambert') then
             IF(surface_albedo.eq.0.0_dp) then
                I_SURF = 0.0_dp
             ELSE
                I_SURF = I0(1) * UMU0 / XPI * exp(-tau_tot/UMU0) * surface_albedo * exp(-tau_tot/UMU) 
             ENDIF
             IF (nmat.gt.1) then
                Q_SURF = 0.0_dp
                U_SURF = 0.0_dp
             END IF
          ELSE
             CALL surface_brdf(nmat, lambda, tetas, tetav(itetav), phiv(iphiv)+180.0, brdf(:,:))
             I_SURF = I0(1) * UMU0 / XPI * exp(-tau_tot/UMU0) * brdf(1,1) * exp(-tau_tot/UMU)
             IF (nmat.gt.1) then
                Q_SURF = I0(1) * UMU0 / XPI * exp(-tau_tot/UMU0) * brdf(2,1) * exp(-tau_tot/UMU)
                U_SURF = I0(1) * UMU0 / XPI * exp(-tau_tot/UMU0) * brdf(3,1) * exp(-tau_tot/UMU)
             ENDIF
          ENDIF

          SvR_res(1,itetav,iphiv) = I_SURF + ( I0(1) / ( 4.0_dp*XPI * ( 1.0_dp + UMU/UMU0 ) ) * I_SINSCA )
          if (nmat .gt. 1) then
             SvR_res(2,itetav,iphiv) = Q_SURF + &
                  (I0(1) / ( 4.0_dp*XPI * ( 1.0_dp + UMU/UMU0 ) ) * Q_SINSCA * COS(2.0_dp * i2 ))
             SvR_res(3,itetav,iphiv) = U_SURF - &
                  (I0(1) / ( 4.0_dp*XPI * ( 1.0_dp + UMU/UMU0 ) ) * Q_SINSCA * SIN(2.0_dp * i2 ))
          endif
          if (nmat .eq. 4) SvR_res(4,itetav,iphiv) = 0.0_dp

       END DO
    END DO

  END SUBROUTINE SINSCA

  !==================================================================

  SUBROUTINE pray(cteta, steta, delta, p11, p21)

    USE MCONSTANTS, only : dp

    implicit none

    !--- Input Arguments
    real(kind=dp), intent(in)  :: cteta
    real(kind=dp), intent(in)  :: steta
    real(kind=dp), intent(in)  :: delta
    real(kind=dp), intent(out) :: p11
    real(kind=dp), intent(out) :: p21

    !--- Local Arguments
    real(kind=dp) :: Gdelta
    real(kind=dp) :: Gdeltap

    !******************************************
    !     Compute Rayleigh phase matrix
    !******************************************

    if (delta .ne. 0D0) then

       GDelta  = (1.0D0 - delta) / (1.0D0 + delta / 2.0D0)
       GDeltap = (1.0D0 - 2.0D0*delta) / (1.0D0 - delta)

       p11 = GDelta * 3.0D0 * (1.0D0 + cteta**2.D0) / 4.0D0 &
            + (1.0D0 - GDelta)                                    ! P11
       P21 = - GDelta * 3.0D0 * steta**2.D0 /4.0D0  ! P21=P12

    else

       p11 = 3.0D0 / 4.0D0 *( 1.0D0 + cteta**2.D0 )   !P11
       p21 = 3.0D0 / 4.0D0 *( cteta**2.0D0 - 1.0D0 )  !P21=P12

    endif

    RETURN 

  END SUBROUTINE pray



END MODULE MSINSCA
