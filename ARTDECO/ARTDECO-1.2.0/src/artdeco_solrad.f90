

MODULE MSOLRAD

  PRIVATE

  PUBLIC :: GET_SOLRAD

CONTAINS

  SUBROUTINE GET_SOLRAD

    USE MCONSTANTS , only : dp, &
         undef_dp
    USE MCOMMON, only : mode, &
         F0 ,&
         nlambda,&
         svin, &
         verbose,&
         wlambda, &
         solrad_nwvl, &
         solrad_wvl, &
         solrad_F0, &
         solrad_solcste, &
         kdis_solrad, &
         kdis_solrad_nwvl, &
         kdis_solrad_wvl, &
         kdis_solrad_F0, &
         kdis_solrad_solcste, &
         kdis_model, &
         kdis_nwvl, &
         kdis_wvlband, &
         lambda_ikdis, &
         thermal_only, &
         get_possol, &
         nsza, &
         sza, &
         lon, &
         lat, &
         day, &
         h_tu, &
         warning, &
         varsol_fact,&
         Fbeam_user_flag, &
         Fbeam_user_value

    !USE nr, ONLY : SORT
    USE nm_nm , ONLY : NM_SORT
    USE MUTILITY, only : LININTPOL, PLINT, XINTEG2
    USE MSOL, only : possol, varsol

    IMPLICIT NONE

    integer :: j
    real(kind=dp) :: resint
    real(kind=dp) :: phisun

    !--------------------

    ALLOCATE(varsol_fact(nsza))
    varsol_fact(:) = 1.0_dp

    ALLOCATE(F0(nlambda))
    F0(:)   = undef_dp 
    svin(1) = 1.0_dp
    svin(2) = 0.0_dp
    svin(3) = 0.0_dp
    svin(4) = 0.0_dp

    ! get solar zenith angle
    if (get_possol) then
       ! compute the sun geometry
       do j = 1, nsza
          call possol(day(j), h_tu(j), lon(j), lat(j), sza(j), phisun)
          if ((sza(j).gt.70.0).and.warning) then
             WRITE(*,*) ''
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ' (in sub. get_solrad)  : WARNING '
             WRITE(*,*) '  Solar zenithal angle must be <= 70 deg'
             WRITE(*,*) '  for the plan-parallel geometry approximation'
             WRITE(*,*) '  to remain accurate '
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ''
          endif
          ! compute scaling factor for the solar constant
          if (mode.ne.'mono') then
             call varsol(day(j), varsol_fact(j))
          endif
       end do
       !CALL SORT(sza)
       CALL NM_SORT(sza)
    endif

    ! set wavelength dependent solar properties
    if (thermal_only) then

       IF (verbose) THEN
          WRITE(*,*) ' '
          WRITE(*,*) '  (sub. get_solrad) No incoming radiation --> F0 = 0.0'
       ENDIF
       do j = 1, nlambda
          F0(j) = 0.0_dp
       end do

    else

       if (mode.eq.'kdis') then

          if (kdis_solrad) then
             ! the solrad integrated in Kdis bands already exist and we just need
             ! to select the used wavelengths (i.e. good index)
             IF (verbose) THEN
                WRITE(*,*) ' '
                WRITE(*,*) '  (sub. get_solrad) Select the solar radiation for the band to be used '
             ENDIF
             DO j = 1, nlambda
                F0(j) = kdis_solrad_F0(lambda_ikdis(j))
             ENDDO
          else
             IF (verbose) THEN
                WRITE(*,*) ' '
                WRITE(*,*) '  (sub. get_solrad) Integration of the solar radiation over each k-dis band '
             ENDIF
             kdis_solrad_nwvl = kdis_nwvl
             ALLOCATE(kdis_solrad_wvl(kdis_solrad_nwvl), kdis_solrad_F0(kdis_solrad_nwvl))
             kdis_solrad_wvl(:) = undef_dp
             kdis_solrad_F0(:)  = undef_dp
             kdis_solrad_solcste = solrad_solcste
             ! We start to integrate the solar radiation over the kdis bands
             ! we do the job for all bands even those not used later in the run 
             ! because we will store it in a file
             do j = 1, kdis_solrad_nwvl
                CALL PLINT(solrad_F0, solrad_wvl, solrad_nwvl, &
                     kdis_wvlband(2,j), kdis_wvlband(3,j), resint )
                kdis_solrad_F0(j)  = resint
                kdis_solrad_wvl(j) = kdis_wvlband(1,j)
             end do
             IF (verbose) THEN
                WRITE(*,*) ' '
                WRITE(*,*) '  (sub. get_solrad) Select the solar radiation for the band to be used '
             ENDIF
             DO j = 1, nlambda
                F0(j) = kdis_solrad_F0(lambda_ikdis(j))
             ENDDO
          endif

       elseif (mode.eq.'mono') then

          IF (verbose) THEN
             WRITE(*,*) ' '
             WRITE(*,*) '  (sub. get_solrad) Reflectivity are computed in mono mode --> F0 = 1.0'
          ENDIF
          do j = 1, nlambda
                !!!Line added on 15/09/14 for F0
                if ( Fbeam_user_flag ) then
                    F0(j)  = Fbeam_user_value
                else 
                    F0(j)  = 1.0_dp
               endif
             !F0(j) = 1.0_dp
             !write(*,*) ' in solrad --> F0 ' , F0(j) 
          end do

       endif

    endif

    IF (verbose) THEN
       WRITE(*,*) ' '
       WRITE(*,*) '  (sub. get_solrad) exit OK'
       WRITE(*,*) ' '
    ENDIF

  END SUBROUTINE GET_SOLRAD

END MODULE MSOLRAD
