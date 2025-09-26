
MODULE MGET_BETAL

  IMPLICIT NONE

  ! suppress those attribute for f2py
  !  PRIVATE :: DEVEL_BL, CALL_DELTAFIT, CALL_POTTER, CALL_DELTAM, GET_F_RECOMP 
  !  PUBLIC  :: GET_BETAL

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE GET_BETAL

    USE MCONSTANTS, only : dp, &
         coeff_trunc_min, &
         eps_betal, &
         undef_dp, &
         nang_min_betal_exp


    !-- Line added on  19/02/14 for ptcle_use_betal, ...

    !--- Line added on 24/02/14 for ptcle_ilambda, ptcle_opt_exist
    USE MCOMMON, only : verbose, &
         nlambda, &
         wlambda, &
         nptcle,  &
         ptcle_hg, &
         ptcle_nbetal, &
         ptcle_betal, &
         
         ptcle_use_betal, &
         ptcle_nbetal_betal_store, &
         ptcle_opt_betal_store, &
         ptcle_betal_store, &
         ptcle_trunc_coeff_betal_store,&
         
         ptcle_opt,&
         
         ptcle_opt_exist,&
         ptcle_ilambda,&
         
         nang_max,&
         ptcle_nang_phasemat, &
         ptcle_u_phasemat, &
         ptcle_type, &
         ptcle_trunccoeff,   &
         ptcle_phasemat, &
         trunc_method, &
         ptcle_trunc_phasemat, &
         ptcle_trunc_normphasemat, &
         rt_model,&
         warning, &
         nmaxbetal, &
         nbetal_in

    USE MUTILITY, only : LININTPOL, spline, splint 
    !USE nr, ONLY : gauleg
    USE nm_nm, ONLY : nm_gaussNodes


    IMPLICIT NONE

    ! local variables
    INTEGER :: j, k, u, v
    INTEGER :: ipart

    INTEGER :: nang
    REAL(kind=dp), allocatable :: gauss_w(:)
    REAL(kind=dp), allocatable :: gauss_u(:)
    REAL(kind=dp), allocatable :: F(:,:)
    INTEGER :: nbetal
    REAL(kind=dp), allocatable :: betal(:,:)
    REAL(kind=dp), allocatable :: u_t2(:)
    REAL(kind=dp) :: trunccoeff

    INTEGER :: ngauss_betal         ! number of gauss points used for the mu grid
    !                                  to perform betal expansion


    REAL(kind=dp) :: norm
    REAL(kind=dp) :: yp1, ypn


    !------------

    IF (verbose) WRITE (*,*) ''

    ALLOCATE(ptcle_nbetal(nptcle,nlambda),&
         ptcle_betal(nptcle, 6, 0:nmaxbetal, nlambda))
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! It may be dependant on ptcle and lambda at some point
    ptcle_nbetal(:,:)    = nbetal_in
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ptcle_betal(:,:,:,:) = undef_dp

    ALLOCATE(ptcle_trunc_normphasemat(nptcle, nlambda),&
         ptcle_trunccoeff(nptcle,nlambda),&
         ptcle_trunc_phasemat(nptcle, 6, nang_max,nlambda))
    ptcle_trunc_normphasemat(:,:) = undef_dp
    ptcle_trunccoeff(:,:)         = 0.0_dp
    ptcle_trunc_phasemat(:,:,:,:) = undef_dp

    ! if rt_model is mcrad1d or sinsca or if  no trunc, no need to go further
    if (((rt_model .ne. 'sinsca') .and. (rt_model .ne. 'mcrad1d')) .or. &
         (trunc_method .ne. 'none') ) then
       
       DO ipart = 1, nptcle
          !--- Line added on 19/02/14
          !---    See if user has their betal , no need to go further
          IF(ptcle_use_betal(ipart).eqv..true.)   then
             DO j = 1, nlambda
                IF (verbose) THEN
                   WRITE (*, FMT='(1x,A,A20,A2,F12.5,1x,A8)') &
                        '  (sub. get_betal) Read from  Betal_l for ', &
                        TRIM(ADJUSTL(ptcle_type(ipart))), ' @',wlambda(j),' microns'             
                ENDIF
                
                if (ptcle_opt_exist(ipart,j).eqv..FALSE.) then
                   write(*,*) ''
                   WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   WRITE(*,*)               ' (in sub. get_betal)  : ERROR '
                   write(*,*)               '   For particle ', TRIM(ADJUSTL(ptcle_type(ipart)))
                   write(*,FMT='(1x,A,F12.5,2x,A)') '   the working wavelenght ', wlambda(j), 'microns'
                   write(*,*)               '   require optical properties such as Cext, albedo '
                   write(*,*)               '   simple difusion, ...'
                   WRITE(*,*)               '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   write(*,*) ''
                   STOP
                else 
                   !---Line added on 20/02/14 for setting trunc coeff
                   u = ptcle_ilambda(ipart, j)
                   ptcle_trunccoeff(ipart, j)           = ptcle_trunc_coeff_betal_store(ipart, u)

                   nbetal = ptcle_nbetal_betal_store(ipart, j) 
                   do k = 0, nbetal
                      ! WRITE(*,*) '!!!!!  betal  !!!!!!!',ptcle_betal_store(ipart,:,k +1,u)
                      !                                        ptcle_betal_store(ipart,2,k + 1,u),&
                      !                                        ptcle_betal_store(ipart,3,k + 1,u),&
                      !                                      ptcle_betal_store(ipart,4,k+1,u)

                      ptcle_betal(ipart, 1, k, j) = ptcle_betal_store(ipart,1,k,u)
                      ptcle_betal(ipart, 2, k, j) = ptcle_betal_store(ipart,2,k,u)
                      ptcle_betal(ipart, 3, k, j) = ptcle_betal_store(ipart,3,k,u)
                      ptcle_betal(ipart, 4, k, j) = ptcle_betal_store(ipart,4,k,u)
                      ptcle_betal(ipart, 5, k, j) = ptcle_betal_store(ipart,5,k,u)
                      ptcle_betal(ipart, 6, k, j) = ptcle_betal_store(ipart,6,k,u)
                   enddo

                endif

             ENDDO
!!!! End modification
          ELSE
             IF ((ptcle_hg(ipart)).and.(trunc_method.eq.'none')) then

                ! Betal corresponding to a H-G phase function
                DO j = 1, nlambda
                   IF (verbose) THEN
                      WRITE (*, FMT='(1x,A,A20,A2,F12.5,1x,A8)') &
                           '  (sub. get_betal) Get H-G Betal_l for ', &
                           TRIM(ADJUSTL(ptcle_type(ipart))), ' @',wlambda(j),' microns'             
                   ENDIF

                   nbetal = ptcle_nbetal(ipart, j)
                   nang   = ptcle_nang_phasemat(ipart, j)

                   do k = 0, nbetal 
                      ptcle_betal(ipart, 1, k, j) = (2.0_dp * k + 1.0_dp) * (ptcle_opt(ipart,3,j)**DBLE(k))
                   enddo

                   ! get the recomposed truncated phase matrix
                   CALL get_F_recomp(nbetal, nang, ptcle_u_phasemat(ipart,1:nang,j), &
                        ptcle_betal(ipart, :, 0:nbetal, j), ptcle_trunc_phasemat(ipart,:, 1:nang,j), &
                        ptcle_trunc_normphasemat(ipart, j))
                   ptcle_trunc_phasemat(ipart,2:6, :,j) = undef_dp
                   ptcle_trunccoeff(ipart, j)           = 0.0_dp

                END DO

             ELSE

                DO j = 1, nlambda

                   nbetal = ptcle_nbetal(ipart, j)
                   nang   = ptcle_nang_phasemat(ipart, j)

                   if (nang.lt.nang_min_betal_exp )then
                      ngauss_betal = nang_min_betal_exp 
                   else
                      ngauss_betal = nang
                   endif

                   ALLOCATE(u_t2(nang))

                   ! !!!!!!!!!!!!!!!!!!!!!!               
                   !ngauss_betal = nbetal * 2 
                   !if (ngauss_betal > nang) ngauss_betal = nang               
                   !if (ngauss_betal < 128) ngauss_betal = 128               
                   ! ngauss_betal = nang
                   ! !!!!!!!!!!!!!!!!!!!!!!!

                   ALLOCATE(gauss_w(ngauss_betal), gauss_u(ngauss_betal), F(6, ngauss_betal), betal(6,0:nbetal))

                   ! we create a gauss points grid needed to perform the Bl expansion
                   !CALL gauleg(minval(ptcle_u_phasemat(ipart,1:nang,j)),maxval(ptcle_u_phasemat(ipart,1:nang,j)), gauss_u, gauss_w)
                   CALL nm_gaussNodes(minval(ptcle_u_phasemat(ipart,1:nang,j)), &
                        maxval(ptcle_u_phasemat(ipart,1:nang,j)),gauss_u, gauss_w)
                   ! we interpolate the phase matrix on the gauss points grid
                   ! INTPOL(fint, xint, ni, xess)
                   ! linear interoplation of fint @ xess
                   ! xint assumed increasing, NO EXTAPOLATION 
                   ! make sure xess belongs to [xint(1),xint(ni)]
                   DO u = 1, 6
                      DO k = 1, ngauss_betal
                         ! F(u,k) = LININTPOL(ptcle_phasemat(ipart,u,1:nang,j), ptcle_u_phasemat(ipart,1:nang,j), nang, gauss_u(k))
                         ! SPLINE interpolation 
                         yp1 = (ptcle_phasemat(ipart,u,2,j) - ptcle_phasemat(ipart,u,1,j) ) / &
                              &  ( ptcle_u_phasemat(ipart,2,j)-ptcle_u_phasemat(ipart,1,j) )
                         ypn = (ptcle_phasemat(ipart,u,nang,j) - ptcle_phasemat(ipart,u,nang-1,j) ) / &
                              &( ptcle_u_phasemat(ipart,nang,j)-ptcle_u_phasemat(ipart,nang-1,j) ) 

                         CALL SPLINE( ptcle_u_phasemat(ipart,1:nang,j), ptcle_phasemat(ipart,u,1:nang,j), &
                              &nang, yp1, ypn,  u_t2)
                         CALL SPLINT( ptcle_u_phasemat(ipart,1:nang,j), ptcle_phasemat(ipart,u,1:nang,j), &
                              &u_t2, nang, gauss_u(k), F(u,k))

                      ENDDO
                   ENDDO

                   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !print*, ' get_betal'
                   !print*, TRIM(ADJUSTL(ptcle_type(ipart))),'  ',wlambda(j)
                   !DO k = 1, ngauss_betal
                   !   print*, gauss_u(k),gauss_w(k), F(1,k), F(4,k), F(5,k), F(6,k)
                   !ENDDO
                   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                   ! TRUNCATION (OR Bl EXPENSION ONLY)
                   select case (trunc_method)

                   case('none')
                      ! Bl expansion without truncation
                      IF (verbose) THEN
                         WRITE (*, FMT='(1x,A51,A20,A2,F12.5,1x,A8)') &
                              '  (sub. get_betal) Get Betal_l (no truncation) for ', &
                              TRIM(ADJUSTL(ptcle_type(ipart))), ' @',wlambda(j),' microns'             
                      ENDIF
                      CALL DEVEL_BL(nbetal, ngauss_betal, gauss_u, gauss_w, F, betal)
                      betal(:, :) = betal(:, :) / betal(1, 0)
                      trunccoeff = 0.0_dp

                   case('dfit')
                      ! get Bl from Delta-fit
                      ! only apply if we work with scalar radiative values and if 
                      ! the particle is spherical
                      IF (verbose) THEN
                         WRITE (*,  FMT='(1x,A47,A20,A2,F12.5,1x,A8)')&
                              '  (sub. get_betal) Get Betal_l (delta-fit) for ', &
                              TRIM(ADJUSTL(ptcle_type(ipart))), ' @',wlambda(j),' microns'             
                      ENDIF
                      CALL CALL_DELTAFIT(ngauss_betal, F, gauss_u, gauss_w, nbetal, betal, trunccoeff)

                   case('potter')
                      ! get Bl from Potter truncation
                      ! only apply if the particle is spherical or to the phase function only
                      IF (verbose) THEN
                         WRITE (*, FMT='(1x,A44,A20,A2,F12.5,1x,A8)')&
                              '  (sub. get_betal) Get Betal_l (Potter) for ', &
                              TRIM(ADJUSTL(ptcle_type(ipart))), ' @',wlambda(j),' microns'             
                      ENDIF
                      CALL CALL_POTTER(ngauss_betal, F, gauss_u, gauss_w, nbetal, betal, trunccoeff)

                   case('dm')
                      ! get Bl from Delta-M truncation
                      ! only apply to the phase function
                      IF (verbose) THEN
                         WRITE (*, FMT='(1x,A45,A20,A2,F12.5,1x,A8)') &
                              '  (sub. get_betal) Get Betal_l (Delta-M) for ', &
                              TRIM(ADJUSTL(ptcle_type(ipart))), ' @',wlambda(j),' microns'             
                      ENDIF
                      CALL CALL_DELTAM(ngauss_betal, F, gauss_u, gauss_w, nbetal, betal, trunccoeff)

                   end select

                   ! SOME SANITY CHECK :
                   if (trunccoeff .lt. coeff_trunc_min) then
                      if (warning) then
                         WRITE(*,*) ' '
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         WRITE(*,*) '  (in sub. get_betal)     WARNING'
                         WRITE (*, *)'       The truncation coefficient is too low.'
                         WRITE (*, *)'       The truncation may not be needed or is badly done'
                         WRITE (*, *)'       (e.g. bad extrapolation in Potter case)'
                         WRITE (*, FMT='(A22,2x,F12.6)') &
                              &      '        trunc_coeff = ',  trunccoeff
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      endif
                   endif
                   if (betal(1,nbetal) .eq. 0.0_dp) then
                      if (warning) then
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         write(*,*) '  (sub. get_betal)  WARNING:'
                         write(*,*) '         if betal(1,nbetal) = 0'
                         write(*,*) '         we force betal(2:6,nbetal) = 0 '
                         write(*,*) '         This may change at some point !! '
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      endif
                      betal(2:6,nbetal) = 0.0_dp
                   endif

                   IF( (betal(1,0) .gt. (1.0_dp+EPS_BETAL) ) .or. (betal(1,0) .lt. (1.0_dp-EPS_BETAL))) THEN
                      !write(*,*) 'betal, EPS_BETAL, 1+EPS_BETAL, 1-EPS_BETAL',betal(1,0),EPS_BETAL,&
                      !(1.0_dp+EPS_BETAL), (1.0_dp-EPS_BETAL)
                      !if (warning) then
                      !   WRITE(*,*) ' '
                      !   WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      !   WRITE(*,*) '  (in sub. get_betal)     WARNING'
                      !   WRITE(*,*) '    Beta_l should be normalized so that Beta_0^11 = 1'
                      !   WRITE (*, FMT='(A20,2x,F12.5,1x,A8)') &
                      !        &        '       wavelength  =', wlambda(j),' microns'  
                      !   WRITE(*,*) '       trunc_method ', TRIM(ADJUSTL(trunc_method))
                      !   WRITE (*, FMT='(A17,2x,F12.6)') &
                      !        &        '       Beta_0^11=',  betal(1,0)
                      !   WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      !endif
                      WRITE(*,*) ' '
                      WRITE(*,*) '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      WRITE(*,*) '  (in sub. get_betal)     ERROR'
                      WRITE(*,*) '    Beta_l should be normalized so that Beta_0^11 = 1'
                      WRITE (*, FMT='(A20,2x,F12.5,1x,A8)') &
                           &        '       wavelength  =', wlambda(j),' microns'  
                      WRITE(*,*) '       trunc_method ', TRIM(ADJUSTL(trunc_method))
                      WRITE (*, FMT='(A17,2x,F12.6)') &
                           &        '       Beta_0^11=',  betal(1,0)
                      WRITE(*,*) '  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      STOP
                   ENDIF

                   ! We always renormalize the Betal
                   !betal(:,:) = betal(:,:) / betal(1,0)

                   ! get the recomposed truncated phase matrix
                   CALL get_F_recomp(nbetal, nang, ptcle_u_phasemat(ipart,1:nang,j), &
                        betal, ptcle_trunc_phasemat(ipart,:, 1:nang,j), ptcle_trunc_normphasemat(ipart, j))


                   if (verbose) then
                      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      norm = 0.0_DP
                      DO k = 1, ngauss_betal
                         norm = norm +  F(1,k) * gauss_w(k)
                      END DO
                      print*, ""
                      print*, "    (get_betal) Phasemat int  (gauss)       =", norm
                      print*, '    (get_betal) Phasemat recomp int (trapz) =', ptcle_trunc_normphasemat(ipart, j)     
                      print*, '    (get_betal) trunc coeff=', trunccoeff
                      !  DO k = 0, nbetal
                      !     print*, k, betal(1,k), betal(2,k), betal(3,k), betal(4,k), betal(5,k), betal(6,k)
                      !  ENDDO
                      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   endif



                   ! implement the common variables
                   ptcle_betal(ipart, :, 0:nbetal, j) = betal(:,0:nbetal)
                   ptcle_trunccoeff(ipart, j)         = trunccoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ! print*, ' get_betal'
                   ! print*, TRIM(ADJUSTL(ptcle_type(ipart))),'  ',wlambda(j)
                   ! print*, trunccoeff
                   ! DO k = 0, nbetal
                   !    print*, betal(1,k), betal(2,k), betal(3,k), betal(4,k), betal(5,k), betal(6,k)
                   ! ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                   DEALLOCATE(gauss_w, gauss_u, F, betal)
                   DEALLOCATE(u_t2)


                ENDDO ! loop on wavelength

             ENDIF ! HG or not
          ENDIF ! not read betal

       ENDDO ! end loop on particles population

    endif

  END SUBROUTINE GET_BETAL

  !----------------------------------------------------------------

  SUBROUTINE CALL_DELTAFIT(ngauss, F, gauss_u, gauss_w, nbetal, betal, trunccoeff)

    USE MCONSTANTS, only : dp 
    USE MCOMMON, only : dfit_thetac, dfit_fitall
    USE MUTILITY, only : LININTPOL

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ngauss
    REAL(kind=dp), INTENT(IN) :: gauss_w(ngauss)
    REAL(kind=dp), INTENT(IN) :: gauss_u(ngauss)
    REAL(kind=dp), INTENT(IN) :: F(6,ngauss)
    INTEGER, INTENT(IN) :: nbetal
    REAL(kind=dp), INTENT(OUT) :: betal(6,0:nbetal)
    REAL(kind=dp), INTENT(OUT) :: trunccoeff

    !--- To call delta-fit trunc
    !- Input Argument
    real(kind=dp) :: trunc_u(ngauss)
    real(kind=dp) :: trunc_w(ngauss)
    real(kind=dp) :: trunc_F(6,ngauss)
    integer :: trunc_flag_DM 

    !---------------

    ! NOTE:  dfit_thetac and dfit_fitall are common variables

    ! init arrays
    trunc_u(:)   = 0.0_dp
    trunc_w(:)   = 0.0_dp
    trunc_F(:,:) = 0.0_dp

    !--- here flag_DM tell whether apply the delta-M scaling or no
    !    trunc_flag_DM = 1 apply the delta-M scaling in the delta fitting
    !    trunc_flag_DM = 0 no delta_M scaling
    trunc_flag_DM  = 1

    trunc_w(:) = gauss_w(:)
    trunc_u(:) = gauss_u(:)
    trunc_F(1,:)  = F(1, :)   ! F11 
    trunc_F(2,:)  = F(2, :)   ! F22
    trunc_F(3,:)  = F(3, :)   ! F33
    trunc_F(4,:)  = F(4, :)   ! F44 
    trunc_F(5,:)  = F(5, :)   ! F21 = F12
    trunc_F(6,:)  = F(6, :)   ! F34 = -F43

    CALL trunc_delta_fit(ngauss,& 
         trunc_u,  &
         trunc_w,  &
         trunc_F,  &
         trunc_flag_DM, &
         dfit_thetac, &
         dfit_fitall, &
         nbetal, &
         betal, &
         trunccoeff)

  END SUBROUTINE CALL_DELTAFIT

  !----------------------------------------------------------------

  SUBROUTINE CALL_POTTER(ngauss, F, gauss_u, gauss_w, nbetal, betal, trunccoeff)

    USE MCONSTANTS, only : dp, deg2rad
    USE MCOMMON, only : verbose, &
         potter_theta_min, &
         potter_theta_max

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ngauss
    REAL(kind=dp), INTENT(IN) :: gauss_w(ngauss)
    REAL(kind=dp), INTENT(IN) :: gauss_u(ngauss)
    REAL(kind=dp), INTENT(IN) :: F(6,ngauss)
    INTEGER, INTENT(IN) :: nbetal
    REAL(kind=dp), INTENT(OUT) :: betal(6,0:nbetal)
    REAL(kind=dp), INTENT(OUT) :: trunccoeff

    !--- To call delta-fit trunc
    !- Input Argument
    real(kind=dp) :: trunc_mu1
    real(kind=dp) :: trunc_mu2
    real(kind=dp) :: trunc_u(ngauss)
    real(kind=dp) :: trunc_w(ngauss)
    real(kind=dp) :: trunc_F(6,ngauss)

    !---------------

    ! init arrays
    trunc_u(:)   = 0.0_dp
    trunc_w(:)   = 0.0_dp
    trunc_F(:,:) = 0.0_dp

    trunc_w(:) = gauss_w(:)
    trunc_u(:) = gauss_u(:)
    trunc_F(1,:)  = F(1, :)   ! F11 
    trunc_F(2,:)  = F(2, :)   ! F22
    trunc_F(3,:)  = F(3, :)   ! F33
    trunc_F(4,:)  = F(4, :)   ! F44 
    trunc_F(5,:)  = F(5, :)   ! F21 = F12
    trunc_F(6,:)  = F(6, :)   ! F34 = -F43

    trunc_mu1      = COS(potter_theta_min*deg2rad)  
    trunc_mu2      = COS(potter_theta_max*deg2rad) 
    IF (verbose .eqv. .true.) then
       write(*,*) '---------------------'
       write(*,*) 'Potter truncation :' 
       write(*,FMT='(A10,F12.5)') ' theta1 = ', potter_theta_min
       write(*,FMT='(A10,F12.5)') ' theta2=  ', potter_theta_max
       write(*,*) '---------------------'
    ENDIF

    CALL trunc_potter(ngauss,    & 
         trunc_mu1, &
         trunc_mu2, &
         trunc_u,   &
         trunc_w,   &
         trunc_F,   &
         nbetal,    &
         betal,     &
         trunccoeff)

  END SUBROUTINE CALL_POTTER

  !----------------------------------------------------------------

  SUBROUTINE CALL_DELTAM(ngauss, F, gauss_u, gauss_w, nbetal, betal, trunccoeff)

    USE MCONSTANTS, only : dp 

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ngauss
    REAL(kind=dp), INTENT(IN) :: gauss_w(ngauss)
    REAL(kind=dp), INTENT(IN) :: gauss_u(ngauss)
    REAL(kind=dp), INTENT(IN) :: F(6,ngauss)
    INTEGER, INTENT(IN) :: nbetal
    REAL(kind=dp), INTENT(OUT) :: betal(6,0:nbetal)
    REAL(kind=dp), INTENT(OUT) :: trunccoeff

    !--- To call delta-fit trunc
    !- Input Argument
    real(kind=dp) :: trunc_u(ngauss)
    real(kind=dp) :: trunc_w(ngauss)
    real(kind=dp) :: trunc_F(6,ngauss)

    !---------------

    ! init arrays
    trunc_u(:)   = 0.0_dp
    trunc_w(:)   = 0.0_dp
    trunc_F(:,:) = 0.0_dp

    trunc_w(:) = gauss_w(:)
    trunc_u(:) = gauss_u(:)
    trunc_F(1,:)  = F(1, :)   ! F11 
    trunc_F(2,:)  = F(2, :)   ! F22
    trunc_F(3,:)  = F(3, :)   ! F33
    trunc_F(4,:)  = F(4, :)   ! F44 
    trunc_F(5,:)  = F(5, :)   ! F21 = F12
    trunc_F(6,:)  = F(6, :)   ! F34 = -F43

    CALL trunc_DM(ngauss,    & 
         trunc_u,   &
         trunc_w,   &
         trunc_F,   &
         nbetal,    &
         betal,     &
         trunccoeff)

  END SUBROUTINE CALL_DELTAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine devel_Bl(nbetal, nangle, u, w, F, betal)

    !  Calculate the expansion coefficients of the scattering matrix in    
    !  generalized spherical functions by numerical integration over the   
    !  scattering angle.                                                   
    !  From the adding paper 
    !  de Haan, Bosma, Hovenier, 1987, A&A, 183, 371     
    USE MCONSTANTS, only : dp
    USE MUTILITY, only : xinteg2

    implicit none

    INTEGER, INTENT(IN) :: nbetal
    INTEGER, INTENT(IN) :: nangle
    REAL(KIND=dp), INTENT(IN) :: u(nangle)
    REAL(KIND=dp), INTENT(IN) :: w(nangle)
    REAL(KIND=dp), INTENT(IN) :: F(6,nangle)
    !                            F(1) = F11
    !                            F(2) = F22
    !                            F(3) = F33
    !                            F(4) = F44
    !                            F(5) = F12
    !                            F(6) = F34

    REAL(KIND=dp), INTENT(OUT) :: betal(6,0:nbetal)
    !                             betal(1) = alpha1
    !                             betal(2) = alpha2
    !                             betal(3) = alpha3
    !                             betal(4) = alpha4
    !                             betal(5) = beta1
    !                             betal(6) = beta2

    ! local variables
    REAL(KIND=dp) :: P00(nangle,2) 
    REAL(KIND=dp) :: P22(nangle,2)
    REAL(KIND=dp) :: P2m2(nangle,2)
    REAL(KIND=dp) :: P02(nangle,2)
    REAL(KIND=dp) :: qroot6
    REAL(KIND=dp) :: fac1
    REAL(KIND=dp) :: fac2
    REAL(KIND=dp) :: fac3
    REAL(KIND=dp) :: sql41
    REAL(KIND=dp) :: sql4
    REAL(KIND=dp) :: twol1
    REAL(KIND=dp) :: tmp2
    REAL(KIND=dp) :: tmp1
    REAL(KIND=dp) :: denom
    REAL(KIND=dp) :: alfap
    REAL(KIND=dp) :: alfam
    REAL(KIND=dp) :: fl

    INTEGER :: i, l
    INTEGER :: lnew
    INTEGER :: lold
    INTEGER :: itmp

    ! ----------------
    !  Initialization 

    qroot6        = -0.25_DP*dsqrt(6._DP)
    betal(:,:)    =  0.0_DP

    !  Start loop over the coefficient index l                             *
    !  first update generalized spherical functions, then calculate betal. *
    !  lold and lnew are pointer-like indices used in recurrence           *
    lnew = 1
    lold = 2

    do l= 0, nbetal

       if (l .eq. 0) then
          ! Adding paper Eq. (77) with m=0                           
          do i=1, nangle
             P00(i,lold) = 1._DP
             P00(i,lnew) = 0._DP
             P02(i,lold) = 0._DP
             P22(i,lold) = 0._DP
             P2m2(i,lold)= 0._DP
             P02(i,lnew) = 0._DP
             P22(i,lnew) = 0._DP
             P2m2(i,lnew)= 0._DP
          enddo
       else
          fac1 = (2._DP*l-1._dp)/dble(l)
          fac2 = dble(l-1._dp)/dble(l)
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
             P02(i,lold) = qroot6*(1._DP-u(i)*u(i))
             P22(i,lold) = 0.25_DP*(1._DP+u(i))*(1._DP+u(i))
             P2m2(i,lold)= 0.25_DP*(1._DP-u(i))*(1._DP-u(i))
             P02(i,lnew) = 0._DP
             P22(i,lnew) = 0._DP
             P2m2(i,lnew)= 0._DP
          enddo
          sql41 = 0._DP
       endif

       if (l .gt. 2) then
          ! Adding paper Eq. (82) with m=0 and m=2                   *
          sql4  = sql41
          sql41 = dsqrt(dble(l*l)-4._dp)
          twol1 = 2._DP*dble(l)-1._dp
          tmp1  = twol1/sql41
          tmp2  = sql4/sql41
          denom = (dble(l)-1._dp)*(dble(l*l)-4._dp)
          fac1  = twol1*(dble(l)-1._dp)*dble(l)/denom
          fac2  = 4._DP*twol1/denom
          fac3  = dble(l)*((dble(l)-1._dp)*(dble(l)-1._dp)-4._dp)/denom
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
       alfap = 0.0_DP
       alfam = 0.0_DP
       do i=1, nangle
          betal(1,l) = betal(1,l) + P00(i,lnew) * w(i) * F(1,i)
          alfap = alfap           + P22(i,lnew) * w(i) * (F(2,i)+F(3,i))
          alfam = alfam           + P2m2(i,lnew)* w(i) * (F(2,i)-F(3,i))
          betal(4,l) = betal(4,l) + P00(i,lnew) * w(i) * F(4,i)
          betal(5,l) = betal(5,l) + P02(i,lnew) * w(i) * F(5,i)
          betal(6,l) = betal(6,l) + P02(i,lnew) * w(i) * F(6,i)
       enddo
       !betal(1,l) = XINTEG2(1, nangle, nangle, u , P00(:,lnew) * F(1,:)          )
       !alfap      = XINTEG2(1, nangle, nangle, u , P22(:,lnew) * (F(2,:)+F(3,:)) )
       !alfam      = XINTEG2(1, nangle, nangle, u , P2m2(:,lnew)* (F(2,:)-F(3,:)) )
       !betal(4,l) = XINTEG2(1, nangle, nangle, u , P00(:,lnew) * F(4,:) )
       !betal(5,l) = XINTEG2(1, nangle, nangle, u , P02(:,lnew) * F(5,:) )
       !betal(6,l) = XINTEG2(1, nangle, nangle, u , P02(:,lnew) * F(6,:) )

       ! Multiply with trivial factors like 0.5_DP*(2*l+1)             
       fl = dble(l)+0.5_DP
       betal(1,l)   =  fl*betal(1,l)
       alfap        =  fl*alfap
       alfam        =  fl*alfam
       betal(4,l)   =  fl*betal(4,l)
       betal(5,l)   =  fl*betal(5,l)
       betal(6,l)   =  fl*betal(6,l)
       betal(2,l)   =  0.5_DP*(alfap+alfam)
       betal(3,l)   =  0.5_DP*(alfap-alfam)

    enddo ! End of loop over index l                                        

  end subroutine devel_Bl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_F_recomp(nbetal, nangle, u, betal, F_recomp, norm)

    !  Calculate the recomposed phase matrix from 
    !  Legendre expension coefficients
    !  From the adding paper 
    !  de Haan, Bosma, Hovenier, 1987, A&A, 183, 371     
    USE MCONSTANTS, only : dp
    USE MUTILITY, only : xinteg2

    implicit none

    INTEGER, INTENT(IN) :: nbetal
    INTEGER, INTENT(IN) :: nangle
    REAL(KIND=dp), INTENT(IN) :: u(nangle)
    REAL(KIND=dp), INTENT(IN) :: betal(6,0:nbetal)
    !                             betal(1) = alpha1
    !                             betal(2) = alpha2
    !                             betal(3) = alpha3
    !                             betal(4) = alpha4
    !                             betal(5) = beta1
    !                             betal(6) = beta2
    REAL(KIND=dp), INTENT(OUT) :: F_recomp(6,nangle)
    REAL(kind=dp), INTENT(OUT) :: norm

    ! local variables
    REAL(KIND=dp) :: F22p33(nangle)
    REAL(KIND=dp) :: F22m33(nangle)
    REAL(KIND=dp) :: P00(nangle,2) 
    REAL(KIND=dp) :: P22(nangle,2)
    REAL(KIND=dp) :: P2m2(nangle,2)
    REAL(KIND=dp) :: P02(nangle,2)
    REAL(KIND=dp) :: qroot6
    REAL(KIND=dp) :: fac1
    REAL(KIND=dp) :: fac2
    REAL(KIND=dp) :: fac3
    REAL(KIND=dp) :: sql41
    REAL(KIND=dp) :: sql4
    REAL(KIND=dp) :: twol1
    REAL(KIND=dp) :: tmp2
    REAL(KIND=dp) :: tmp1
    REAL(KIND=dp) :: denom

    INTEGER :: i, l
    INTEGER :: lnew
    INTEGER :: lold
    INTEGER :: itmp

    ! ====================
    !  Initialization    

    qroot6        = -0.25_DP*dsqrt(6._DP)
    F_recomp(:,:) = 0.0_DP
    F22p33(:)     = 0.0_DP
    F22m33(:)     = 0.0_DP

    !  Start loop over the coefficient index l                             *
    !  first update generalized spherical functions, then calculate betal. *
    !  lold and lnew are pointer-like indices used in recurrence           *
    lnew = 1
    lold = 2

    do l= 0, nbetal

       if (l .eq. 0) then
          ! Adding paper Eq. (77) with m=0                           
          do i=1, nangle
             P00(i,lold) = 1._DP
             P00(i,lnew) = 0._DP
             P02(i,lold) = 0._DP
             P22(i,lold) = 0._DP
             P2m2(i,lold)= 0._DP
             P02(i,lnew) = 0._DP
             P22(i,lnew) = 0._DP
             P2m2(i,lnew)= 0._DP
          enddo
       else
          fac1 = (2._DP*l-1._dp)/dble(l)
          fac2 = dble(l-1._dp)/dble(l)
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
             P02(i,lold) = qroot6*(1._DP-u(i)*u(i))
             P22(i,lold) = 0.25_DP*(1._DP+u(i))*(1._DP+u(i))
             P2m2(i,lold)= 0.25_DP*(1._DP-u(i))*(1._DP-u(i))
             P02(i,lnew) = 0._DP
             P22(i,lnew) = 0._DP
             P2m2(i,lnew)= 0._DP
          enddo
          sql41 = 0._DP
       endif

       if (l .gt. 2) then
          ! Adding paper Eq. (82) with m=0 and m=2                   *
          sql4  = sql41
          sql41 = dsqrt(dble(l*l)-4._dp)
          twol1 = 2._DP*dble(l)-1._dp
          tmp1  = twol1/sql41
          tmp2  = sql4/sql41
          denom = (dble(l)-1._dp)*(dble(l*l)-4._dp)
          fac1  = twol1*(dble(l)-1._dp)*dble(l)/denom
          fac2  = 4._DP*twol1/denom
          fac3  = dble(l)*((dble(l)-1._dp)*(dble(l)-1._dp)-4._dp)/denom
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

       ! compute truncated phase matrix
       do i = 1, nangle
          F_recomp(1,i) = F_recomp(1,i) + betal(1,l)                * P00(i,lnew)
          F22p33(i)     = F22p33(i)     + (betal(2,l)+betal(3,l))   * P22(i,lnew)
          F22m33(i)     = F22m33(i)     + (betal(2,l)-betal(3,l))   * P2m2(i,lnew)
          F_recomp(4,i) = F_recomp(4,i) + betal(4,l)                * P00(i,lnew)
          F_recomp(5,i) = F_recomp(5,i) + betal(5,l)                * P02(i,lnew)
          F_recomp(6,i) = F_recomp(6,i) + betal(6,l)                * P02(i,lnew)
       enddo

    enddo ! End of loop over index l                                        

    F_recomp(2,:) = ( F22p33(:) + F22m33(:) ) / 2.0
    F_recomp(3,:) = ( F22p33(:) - F22m33(:) ) / 2.0

    ! compute the integral of F11
    norm = XINTEG2(1, nangle, nangle, u, F_recomp(1,:))

  end subroutine get_F_recomp

END MODULE MGET_BETAL
