

MODULE MSURFACE_BRDF

  IMPLICIT NONE

  PUBLIC  :: surface_brdf
  PRIVATE 

contains

  subroutine surface_brdf(nmat, wl, ts_in, tv_in, fi_in, brdf)

    ! This subroutine is the one to be called from RT codes.
    ! It call must be 'surface type' independant.
    ! All surface type dependant variables are then accessed
    ! from the inside of the subroutine
    ! through the MCOMMONS (invisible from RT codes).

    USE MBRDF_OCEAN, only : brdf_ocean
    USE MBRDF_PAVEL, only : BRM 
    !USE MCONSTANTS, only : dp, deg2rad
    USE MCONSTANTS, only : dp, deg2rad, eps_lamb
    USE MCOMMON, only : warning, &
         nlambda, &
         wlambda, &
         surface_n_par_brdf, &
         surface_par_brdf, & ! surface_par_brdf(surface_n_par_brdf, nlambda) 
         surface_n_par_bpdf, &
         surface_par_bpdf, & ! surface_par_bpdf(surface_n_par_bpdf, nlambda) 
         surface_brdf_model, &
         wind_spd, &
         ocean_distr, &
         ocean_xsal, &
         ocean_pcl, &
         ocean_shadow, &
         ocean_Rsw, &
         rt_model

    implicit none

    integer      , INTENT(IN) :: nmat   ! number of Stokes parameters
    real(kind=dp), INTENT(IN) :: wl     ! wavelength in microns
    real(kind=dp), INTENT(IN) :: ts_in  ! incident zenith angle (in deg. between 0 and 90)
    real(kind=dp), INTENT(IN) :: tv_in  ! reflected zenith angle (in deg. between 0 and 90)
    real(kind=dp), INTENT(IN) :: fi_in  ! relative azimuth in degre (incident-reflected)
    real(kind=dp), INTENT(OUT), dimension(nmat,nmat) :: brdf

    ! local variables

    real(kind=dp) :: azw ! azim. of the sun - azim. of the wind (in deg. between 0 and 360)

    ! for Pavel's BRM
    integer :: indwl
    Real(kind=dp)    :: par_brdf(3)
    Real(kind=dp)    :: Cv
    complex(kind=dp) :: m
    REAL(kind=dp)    :: Rip(4,4)

    ! --------------

    brdf(:,:) = 0.0_dp

    select case (surface_brdf_model)

    case('ocean')
       azw = 0.0_dp ! unused in case of isotropic surface       
       !call brdf_ocean(nmat, wl, wind_spd, ocean_distr, ocean_shadow, ocean_xsal, ocean_pcl, azw, ts_in, tv_in, -fi_in, brdf)
       call brdf_ocean(nmat, wl, wind_spd, ocean_distr, ocean_shadow, ocean_xsal, ocean_pcl, azw, ocean_Rsw, &
            ts_in, tv_in, -fi_in, brdf)


       ! The PHI rotation convention is taken to be the one of DISORT, DOAD, SINSCA
       ! (Note: the one of MCRAD1D is the opposite)
       if (nmat.gt.1) then
          brdf(1,3) = - brdf(1,3)
          brdf(2,3) = - brdf(2,3)
          brdf(3,1) = - brdf(3,1)
          brdf(3,2) = - brdf(3,2)
       end if
       if (nmat.gt.3) then
          brdf(1,4) = - brdf(1,4)
          brdf(2,4) = - brdf(2,4)
          brdf(4,1) = - brdf(4,1)
          brdf(4,2) = - brdf(4,2)
       endif

    case('lisp_ross')

       if (nmat.gt.1) then
          write(*,*) ' '
          write(*,*) ' ERROR :'
          write(*,*) ' The azimuthal angle signe convention must be tested before Maignan BPDF use' 
          write(*,*) ' '
          STOP
       endif

       indwl=1
       !DO WHILE(wlambda(indwl).ne.wl)
       DO WHILE(ABS(wlambda(indwl)-wl)>eps_lamb)
          indwl = indwl + 1
       END DO
       m  = CMPLX(surface_par_bpdf(2,indwl), surface_par_bpdf(3,indwl), kind=dp)
       Cv = surface_par_bpdf(1,indwl) 
       par_brdf(1) = surface_par_brdf(1,indwl) ! klambda
       par_brdf(2) = surface_par_brdf(2,indwl) ! k2
       par_brdf(3) = surface_par_brdf(3,indwl) ! k1
       CALL BRM(warning, nmat, 1, COS(ts_in*deg2rad), COS(tv_in*deg2rad), fi_in, m, 3, par_brdf(1:3), Cv, Rip)
       brdf(1:nmat, 1:nmat) = Rip(1:nmat, 1:nmat)

    case default 

       write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write(*,*) '  (sub. surface_brdf) ERROR '
       write(*,*) '    The BRDF ', TRIM(ADJUSTL(surface_brdf_model))
       write(*,*) '    is not implemented'
       write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       STOP

    end select


  end subroutine surface_brdf


END MODULE MSURFACE_BRDF
