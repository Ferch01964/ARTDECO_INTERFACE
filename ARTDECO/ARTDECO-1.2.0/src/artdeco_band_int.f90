

MODULE MBAND_INT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GET_BAND_INT

CONTAINS

  SUBROUTINE GET_BAND_INT

    USE MCONSTANTS, only : dp
    USE MCOMMON, only : verbose, &
         mode, &
         filter, &
         fluxint_out, &
         lambda_bound, &
         kdis_wvlband, &
         lambda_ikdis, &
         nsza,&
         nalt_atm, &
         flux_out, &
         nlambda, &
         wlambda, &
         nfilter, &
         bound_lamb_filter, &
         trans_filter, &
         lamb_filter, &
         nlamb_filter, &
         SvR_filter, &
         SvR_sig_filter, &
         flux_filter, &
         SvR, &
         SvR_sig, &
         rt_cputime_filter, &
         rt_cputime, &
         flux_only

    USE MUTILITY, only : LININTPOL

    IMPLICIT NONE

    integer :: ialt
    integer :: i
    integer :: j
    REAL(kind=dp) :: wght
    REAL(kind=dp) :: trans

    ! !!!!!!!!!!!!!

    if (mode.eq.'kdis') then

       if (filter) then

          IF (verbose) THEN
             WRITE (*,*) ''
             WRITE (*, *)'  (sub. get_band_int) compute radiative values integrated over filter'
          ENDIF
          rt_cputime_filter(:)      = 0.0_dp
          SvR_filter(:,:,:,:,:)     = 0.0_dp
          SvR_sig_filter(:,:,:,:,:) = 0.0_dp
          flux_filter(:,:,:,:)      = 0.0_dp

          DO i = 1 , nfilter
             DO j = 1, nlambda
                if ( (wlambda(j).ge.bound_lamb_filter(1,i)).and.&
                     (wlambda(j).le.bound_lamb_filter(2,i))) then

                   trans = LININTPOL(trans_filter(i,1:nlamb_filter(i)), &
                        lamb_filter(i,1:nlamb_filter(i)), &
                        nlamb_filter(i), wlambda(j))
                   if (flux_only.eqv..FALSE.) then
                     !WRITE (*,*) ''
                     !WRITE (*, *)'  (sub. get_band_int) trans', trans
                     !WRITE (*,*) ''
                     !WRITE (*, *)'  (sub. get_band_int) SvR(:,:,:,:,j)', SvR(:,:,:,:,j)
                     !WRITE (*,*) ''


                      SvR_filter(:,:,:,:,i)     = SvR_filter(:,:,:,:,i)     + &
                           trans * SvR(:,:,:,:,j)
                      SvR_sig_filter(:,:,:,:,i) = SvR_sig_filter(:,:,:,:,i) + &
                           (trans * SvR_sig(:,:,:,:,j))**2.0_dp
                   endif
                   flux_filter(:,:,:,i)      =  flux_filter(:,:,:,i)     +  &
                        trans * flux_out(:,:,:,j)  

                   rt_cputime_filter(i) = rt_cputime_filter(i) + rt_cputime(j)

                end if
             END DO
              if (flux_only.eqv..FALSE.) SvR_sig_filter(:,:,:,:,i) = SQRT(SvR_sig_filter(:,:,:,:,i))
          END DO

       else

          IF (verbose) THEN
             WRITE (*,*) ''
             WRITE (*, *)'  (sub. get_band_int) compute integrated radiative values '
          ENDIF

          fluxint_out(:,:,:) = 0.0_dp
          ! We compute the integrated (sum indeed) flux over the wavelength range

          ! special treatement for first band
          if (lambda_bound(1).ne.kdis_wvlband(2,lambda_ikdis(1))) then
             wght = (kdis_wvlband(3,lambda_ikdis(1)) - lambda_bound(1)) / &
                    (kdis_wvlband(3,lambda_ikdis(1)) - kdis_wvlband(2,lambda_ikdis(1)))
          else
             wght = 1.0_dp
          endif
          !print*, lambda_bound(1), kdis_wvlband(2,lambda_ikdis(1)), &
          !     kdis_wvlband(3,lambda_ikdis(1)), wght
          do i = 1, nsza
             DO ialt = 1, nalt_atm
                fluxint_out(i,1,ialt) = fluxint_out(i,1,ialt) + flux_out(i,1,ialt,1) * wght
                fluxint_out(i,2,ialt) = fluxint_out(i,2,ialt) + flux_out(i,2,ialt,1) * wght
                fluxint_out(i,3,ialt) = fluxint_out(i,3,ialt) + flux_out(i,3,ialt,1) * wght
             end DO
          end do
          ! mid band
          do j = 2, nlambda-1
             do i = 1, nsza
                DO ialt = 1, nalt_atm
                   fluxint_out(i,1,ialt) = fluxint_out(i,1,ialt) + flux_out(i,1,ialt,j)
                   fluxint_out(i,2,ialt) = fluxint_out(i,2,ialt) + flux_out(i,2,ialt,j) 
                   fluxint_out(i,3,ialt) = fluxint_out(i,3,ialt) + flux_out(i,3,ialt,j) 
                end DO
             end do
          end do
          ! last band
          if (lambda_bound(2).ne.kdis_wvlband(3,lambda_ikdis(nlambda))) then
             wght = ( lambda_bound(2) - kdis_wvlband(2,lambda_ikdis(nlambda)) ) / &
                     (kdis_wvlband(3,lambda_ikdis(nlambda)) - kdis_wvlband(2,lambda_ikdis(nlambda)))
          else
             wght = 1.0_dp
          endif
          !print*, lambda_bound(2), kdis_wvlband(2,lambda_ikdis(nlambda)), &
          !     kdis_wvlband(3,lambda_ikdis(nlambda)), wght
          !write(*,*) ' wght in artdeco_band_int ', wght
          do i = 1, nsza
             DO ialt = 1, nalt_atm
                fluxint_out(i,1,ialt) = fluxint_out(i,1,ialt) + flux_out(i,1,ialt,nlambda) * wght
                fluxint_out(i,2,ialt) = fluxint_out(i,2,ialt) + flux_out(i,2,ialt,nlambda) * wght
                fluxint_out(i,3,ialt) = fluxint_out(i,3,ialt) + flux_out(i,3,ialt,nlambda) * wght
             end DO
          end do

       endif

    endif

  END SUBROUTINE GET_BAND_INT

END MODULE MBAND_INT
