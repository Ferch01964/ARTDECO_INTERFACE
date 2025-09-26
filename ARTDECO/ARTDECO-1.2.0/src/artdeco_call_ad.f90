
MODULE MCALL_AD

  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: CALL_AD

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE CALL_AD

    USE MCONSTANTS, only : dp, xpi

    !!!! Line added on 25/03/14 for ad_flag_user_Fbeam
    USE MCOMMON, only : verbose, &
         ad_nmug, &
         ad_epsilon, &
         ad_nfoumax, &
         ad_mstop, &
         nmat, &
         nvaa, &
         nvza, &
         nsza, &
         vaa, &
         vza, &
         sza, &
         surface_type, &
         surface_family, &
         surface_albedo, &
         wlambda, &
         nlambda, &
         nlayers, &
         nmaxbetal, &
         nbetal_layers, &
         betal_layers, &
         tot_dtau_mono_layers, &
         ssa_mono_layers, &
         n_aitaui_mono_layers,&
         ai_mono_layers,&
         n_aitaui_total, &
         nmax_aitaui_mono_layers, &
         SvR, &
         F0, &
         Svin, &
         rt_cputime, &
         mode, &
         varsol_fact
         !,&
         !Fbeam_user_flag, &
         !Fbeam_user_value


    USE MLAYERS, only : GET_MONO_LAYERS

    IMPLICIT NONE
    
    !!! Line added on 25/03/14
    !REAL(kind=dp) :: ad_fbeam

    ! MUST BE EQUAL to the values in ad/addoub.f
    INTEGER, PARAMETER :: ad_nmu    = 2
    INTEGER, PARAMETER :: ad_nbelow = 2 
    INTEGER, PARAMETER :: ad_nall   = 3 

    ! adding inputs
    ! array dimensions
    INTEGER :: ad_ndmu
    INTEGER :: ad_ndsup
    INTEGER :: ad_ndlay
    INTEGER :: ad_ndcoef
    INTEGER :: ad_ndgeom
    !------
    REAL(kind=dp), allocatable :: ad_a(:)
    REAL(kind=dp), allocatable :: ad_b(:)
    REAL(kind=dp), allocatable :: ad_coefs(:,:,:,:)
    INTEGER, allocatable       :: ad_ncoefs(:)
    INTEGER                    :: ad_nlay
    REAL(kind=dp), allocatable :: ad_theta(:)
    REAL(kind=dp), allocatable :: ad_theta0(:)
    REAL(kind=dp), allocatable :: ad_phi(:)
    INTEGER :: ad_ngeom
    INTEGER :: ad_nmat        
    INTEGER :: ad_iface
    REAL(kind=dp) :: ad_xm

    ! adding outputs
    REAL(kind=dp), allocatable :: ad_R(:,:,:)
    REAL(kind=dp), allocatable :: ad_T(:,:,:)
    REAL(kind=dp), allocatable :: ad_R1(:,:,:)
    REAL(kind=dp), allocatable :: ad_T1(:,:,:)
    REAL(kind=dp), allocatable :: ad_Rst(:,:,:)
    REAL(kind=dp), allocatable :: ad_Tst(:,:,:)
    REAL(kind=dp), allocatable :: ad_R1st(:,:,:)
    REAL(kind=dp), allocatable :: ad_T1st(:,:,:)
    REAL(kind=dp), allocatable :: ad_Ri(:,:,:)
    REAL(kind=dp), allocatable :: ad_R1i(:,:,:)
    REAL(kind=dp), allocatable :: ad_Di(:,:,:)
    REAL(kind=dp), allocatable :: ad_D1i(:,:,:)
    REAL(kind=dp), allocatable :: ad_Rflux(:,:, :, :, :, :)
    REAL(kind=dp), allocatable :: ad_Tflux(:,:, :, :, :, :)
    REAL(kind=dp), allocatable :: ad_URU(  :,:, :, :)
    REAL(kind=dp), allocatable :: ad_UTU(  :,:,:, :)
    REAL(kind=dp), allocatable :: ad_Riflux(:,:, :, :, :, :)
    REAL(kind=dp), allocatable :: ad_Tiflux(:,:, :, :, :, :)
    REAL(kind=dp), allocatable :: ad_Diflux(:,:, :, :, :, :)
    REAL(kind=dp), allocatable :: ad_URUi(  :,:, :, :)
    REAL(kind=dp), allocatable :: ad_UDUi(  :,:, :, :)
    REAL(kind=dp), allocatable :: ad_UTUi(  :,:, :, :)

    ! addlam outputs
    REAL(kind=dp), allocatable :: ad_RL( :,:, :)
    REAL(kind=dp), allocatable :: ad_R1L(:,:, :)
    REAL(kind=dp), allocatable :: ad_DL( :,:, :)
    REAL(kind=dp), allocatable :: ad_D1L( :,:, :)
    REAL(kind=dp), allocatable :: ad_RLflux(:,:, :, :, :, :)
    REAL(kind=dp), allocatable :: ad_DLflux(:,:, :, :, :, :)
    REAL(kind=dp), allocatable :: ad_URUL(  :,:, :, :)

    ! stokout output
    REAL(kind=dp), allocatable :: ad_SvR(:,:)

    INTEGER :: i, u, k, j
    INTEGER :: igeom
    INTEGER :: iai

    ! time
    REAL :: tcpu1
    REAL :: tcpu2

    !-------------------
    ! NOTE : ad_epsilon, ad_nmug, ad_mstop and ad_nfoumax were read in a file

    if (ad_nmug.lt.2) then
       WRITE(*,*) ' '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. call_ad) : ERROR   '
       WRITE(*,*) ' Nstream must be > 2         '
       WRITE(*,*) ' Change it in ad_spec.dat    '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' '
       STOP
    endif

    IF (verbose) THEN
       WRITE (*,*) ''
       WRITE (*, *)'  (sub. call_ad) set ADDOUB arguments'
       WRITE (*,*) ''
       WRITE (*,FMT='(A,1x,I5,1x,A)') '                  ADDOUB will be called ',n_aitaui_total, ' times'
       WRITE (*,FMT='(A,1x,I6,A)') '                  (maximum number of call per wavelength is ', nmax_aitaui_mono_layers,')'
    ENDIF

    if (surface_family .eq. 'brdf') then
       WRITE (*,*) ''
       WRITE (*,*) '  (sub. call_ad) surface BRDF is not implemented in ADDOUB'
       WRITE (*,*) ''
       STOP
    endif

    !------- ALLOCATE ARRAYS WHOSE DIMENSION DOES NOT DEPEND ON WAVELENGTH
    ! array dimensions (that does not depend on wlambda)
    ! It was PARAMETERS in adding previously
    ! and no are arguments. We then need to allocate
    ! arrays with the right dimension
    ad_ndlay  = nlayers
    ad_ndgeom = nsza * nvza * nvaa
    ad_ndmu   = ad_nmug + ad_ndgeom + 1 ! as defined in ad_addoub.f
    ad_ndsup  = ad_ndmu * nmat      + 1    ! as defined in ad_addoub.f
    allocate(ad_a(ad_ndlay),&
         ad_b(ad_ndlay),&
         ad_ncoefs(ad_ndlay),&
         ad_theta(ad_ndgeom),&
         ad_theta0(ad_ndgeom),&
         ad_phi(ad_ndgeom),&
         ad_R(4,4,ad_ndgeom),&
         ad_T(4,4,ad_ndgeom),&
         ad_R1(4,4,ad_ndgeom),&
         ad_T1(4,4,ad_ndgeom),&
         ad_Rst(4,4,ad_ndgeom),&
         ad_Tst(4,4,ad_ndgeom),&
         ad_R1st(4,4,ad_ndgeom),&
         ad_T1st(4,4,ad_ndgeom),&
         ad_Ri(4,4,ad_ndgeom),&
         ad_R1i(4,4,ad_ndgeom),&
         ad_Di(4,4,ad_ndgeom),&
         ad_D1i(4,4,ad_ndgeom),&
         ad_Rflux(4,4, ad_ndgeom, ad_nmu, ad_nbelow, ad_nall),&
         ad_Tflux(4,4, ad_ndgeom, ad_nmu, ad_nbelow, ad_nall),&
         ad_URU(  4,4,              ad_nbelow, ad_nall),&
         ad_UTU(  4,4,              ad_nbelow, ad_nall),&
         ad_Riflux(4,4, ad_ndgeom, ad_nmu, ad_nbelow, ad_nall),&
         ad_Tiflux(4,4, ad_ndgeom, ad_nmu, ad_nbelow, ad_nall),&
         ad_Diflux(4,4, ad_ndgeom, ad_nmu, ad_nbelow, ad_nall),&
         ad_URUi(  4,4,              ad_nbelow, ad_nall),&
         ad_UDUi(  4,4,              ad_nbelow, ad_nall),&
         ad_UTUi(  4,4,              ad_nbelow, ad_nall),&
         ad_RL( 4,4, ad_ndgeom),&
         ad_R1L(4,4, ad_ndgeom),&
         ad_DL( 4,4, ad_ndgeom),&
         ad_D1L( 4,4, ad_ndgeom),&
         ad_RLflux(4,4, ad_ndgeom, ad_nmu, ad_nbelow, ad_nall),&
         ad_DLflux(4,4, ad_ndgeom, ad_nmu, ad_nbelow, ad_nall),&
         ad_URUL(  4,4,                    ad_nbelow, ad_nall),&
         ad_SvR(4, ad_ndgeom))
    !------------
    ! we use nmaxbetal to dimension ad_coefs array but we will send
    ! ad_coefs(:,:,0:ad_ndcoef,:) array to ADDOUB
    ! We must dimension this array once for all here 
    ! cause we can not allocate it multiple times
    ! with the right number of Nmom depending on the 
    ! wavelength when perfoming parallel computing.
    allocate(ad_coefs(4,4,0:nmaxbetal,ad_ndlay))
    !------------

    !-------- IMPLEMENT VARIABLES THAT DOES NOT DEPEND ON WAVELENGTH
    ad_nmat = nmat
    ad_nlay = nlayers
    ! unused variables
    ! interface specification
    ad_iface = 0 ! we assume no interface
    ad_xm = 1.33_dp
    !------------
    ! get the geometry variable to feed adding with
    ad_ngeom = nsza * nvza * nvaa
    igeom = 0
    DO i = 1, nsza
       DO u = 1, nvza
          DO k = 1, nvaa
             igeom = igeom + 1
             ad_theta0(igeom) = sza(i)
             ad_theta(igeom)  = vza(u)
             ad_phi(igeom)    = vaa(k)
          ENDDO
       ENDDO
    ENDDO

    IF (verbose) THEN
       WRITE (*,*) ''
    ENDIF

    DO j = 1, nlambda

       CALL GET_MONO_LAYERS(j)

       !----------
       ! array dimensions (that depends on wlambda)
       ! It was PARAMETERS in adding previously
       ! and no are arguments. We then need to allocate
       ! arrays with the right dimension
       ad_ndcoef = nbetal_layers(j)

       !------------
       ! get layer properties for adding
       ad_a(:)           = 0.0_dp
       ad_b(:)           = 0.0_dp
       ad_coefs(:,:,:,:) = 0.0_dp
       ad_ncoefs(:)      = 0

       DO i = 1, ad_nlay
          ! layers must be sorted in ascending order for adding
          ad_ncoefs(i) = nbetal_layers(j)
          DO u = 0, ad_ncoefs(i)
             ad_coefs(1,1,u,i) = betal_layers(nlayers+1-i, 1, u, j)    ! alpha1
             ad_coefs(1,2,u,i) = betal_layers(nlayers+1-i, 5, u, j)    ! beta1
             ad_coefs(2,1,u,i) = betal_layers(nlayers+1-i, 5, u, j)    ! beta1
             ad_coefs(2,2,u,i) = betal_layers(nlayers+1-i, 2, u, j)    ! alpha2 
             ad_coefs(3,3,u,i) = betal_layers(nlayers+1-i, 3, u, j)    ! alpha3
             ad_coefs(3,4,u,i) = betal_layers(nlayers+1-i, 6, u, j)    ! beta2 
             ad_coefs(4,3,u,i) = -betal_layers(nlayers+1-i, 6, u, j)   ! -beta2
             ad_coefs(4,4,u,i) = betal_layers(nlayers+1-i, 4, u, j)    ! alpha4 
          ENDDO
       ENDDO

       ! call the adding code
       IF (verbose) then
          if (mode.eq.'kdis') then
             WRITE (*, FMT='(1x,A,2x,F10.5,1x,A, 1x, I4)')  &
                  &'                 run ADDOUB for wlambda=',wlambda(j),' microns with nai=', n_aitaui_mono_layers
          else
             WRITE (*, FMT='(1x,A,2x,F10.5,1x,A)')  &
                  &'                 run ADDOUB for wlambda=',wlambda(j),' microns'  
          endif
       ENDIF

       CALL CPU_TIME(tcpu1)

       SvR(:,:,:,:,j) = 0.0_dp

       do iai = 1, n_aitaui_mono_layers

          DO i = 1, ad_nlay
             ad_a(i) = ssa_mono_layers(nlayers+1-i,iai)
             ad_b(i) = tot_dtau_mono_layers(nlayers+1-i,iai)
          END DO

          !write(*,*) 'ad_ndmu, ad_ndsup, ad_ndlay, ad_ndcoef, ad_ndgeom'
          !write(*,*) ad_ndmu, ad_ndsup, ad_ndlay, ad_ndcoef, ad_ndgeom
          CALL adding(ad_ndmu, ad_ndsup, ad_ndlay, ad_ndcoef, ad_ndgeom         &
               , ad_a, ad_b, ad_coefs(:,:,0:ad_ndcoef,:), ad_ncoefs, ad_nlay    &  
               , ad_theta, ad_theta0, ad_phi, ad_ngeom                          &
               , ad_nmug, ad_nfoumax, ad_mstop, ad_nmat, ad_epsilon             &
               , ad_iface, ad_xm                                                &
               , ad_R, ad_T, ad_R1, ad_T1, ad_Rst, ad_Tst, ad_R1st, ad_T1st     &
               , ad_Ri, ad_R1i, ad_Di, ad_D1i                                   &
               , ad_Rflux, ad_Tflux, ad_URU, ad_UTU                             &
               , ad_Riflux, ad_Tiflux, ad_Diflux, ad_URUi, ad_UTUi, ad_UDUi )
          !  Add a Lambertian surface, if desired     
          if ((surface_family.eq.'lambert').and.(surface_albedo(j).gt.0.0_dp)) then
             CALL addlam( ad_ndgeom, surface_albedo(j), ad_R, ad_R1, ad_T, ad_T1    & 
                  , ad_Rflux, ad_Tflux, ad_URU, ad_UTU                              &
                  , ad_ngeom, ad_nmat                                               &
                  , ad_RL, ad_R1L, ad_DL, ad_D1L, ad_RLflux, ad_DLflux, ad_URUL ) 
             CALL Stokout( ad_ndgeom, ad_theta0, ad_nmat, ad_ngeom, Svin, ad_RL, ad_SvR )
          else
             CALL Stokout( ad_ndgeom, ad_theta0, ad_nmat, ad_ngeom, Svin, ad_R, ad_SvR )
          endif
          igeom = 0
          ! implement the common variable for stokes vector for reflected intensities for each geometry
          DO i = 1, nsza
             DO u = 1, nvza
                DO k = 1, nvaa
                   igeom = igeom + 1
                   !!! Line added on  25/03/14 for Fbeam
                   !if( Fbeam_user_flag) then
                   !     ad_fbeam =  Fbeam_user_value
                   !else 
                   !     ad_fbeam = F0(j) * varsol_fact(i)/xpi
                   !endif
                   !WRITE(*,*) ' '
                   !WRITE(*,*) '  (sub. artdeco_call_ad) '
                   !WRITE(*,*) ' '
                   !WRITE(*,FMT='(1x,A35,1x,f10.5)')'>>>>>>                       ad_fbeam : ',&
                   !ad_fbeam                   
                   
                   !SvR(i,u,k,1,j)                  = SvR(i,u,k,1,j)  + ad_SvR(1,igeom) * ad_fbeam &
                   !     * ai_mono_layers(iai) 
                   !if (nmat .gt. 1) SvR(i,u,k,2,j) = SvR(i,u,k,2,j)  + ad_SvR(2,igeom) * ad_fbeam & 
                   !     * ai_mono_layers(iai) 
                   ! we multiply U by -1 to get the same convention as in other RT codes
                   !if (nmat .gt. 1) SvR(i,u,k,3,j) = SvR(i,u,k,3,j)  + ad_SvR(3,igeom) * ad_fbeam & 
                   !     * ai_mono_layers(iai) 
                   !if (nmat .gt. 3) SvR(i,u,k,4,j) = SvR(i,u,k,4,j)  + ad_SvR(4,igeom) * ad_fbeam & 
                   !     * ai_mono_layers(iai) 

                   SvR(i,u,k,1,j)                  = SvR(i,u,k,1,j)  + ad_SvR(1,igeom) * F0(j) * varsol_fact(i) &
                        * ai_mono_layers(iai) / xpi
                   if (nmat .gt. 1) SvR(i,u,k,2,j) = SvR(i,u,k,2,j)  + ad_SvR(2,igeom) * F0(j) * varsol_fact(i) & 
                        * ai_mono_layers(iai) / xpi
                   ! we multiply U by -1 to get the same convention as in other RT codes
                   if (nmat .gt. 1) SvR(i,u,k,3,j) = SvR(i,u,k,3,j)  + ad_SvR(3,igeom) * F0(j) * varsol_fact(i) & 
                        * ai_mono_layers(iai) / xpi
                   if (nmat .gt. 3) SvR(i,u,k,4,j) = SvR(i,u,k,4,j)  + ad_SvR(4,igeom) * F0(j) * varsol_fact(i) & 
                        * ai_mono_layers(iai) / xpi
                ENDDO
             ENDDO
          ENDDO

       end do! iai

       CALL CPU_TIME(tcpu2)
       rt_cputime(j) = tcpu2-tcpu1
       ! IF (verbose) THEN
       !    WRITE (*,FMT='(1x,A36,1PE14.5, A9)') '  (sub. call_ad) CPU time     ', tcpu2-tcpu1, ' sec(s) -'
       ! ENDIF

    ENDDO ! wlambda

    DEALLOCATE(ad_coefs)
    DEALLOCATE(ad_a,ad_b,ad_ncoefs,ad_theta,ad_theta0,&
         ad_phi,ad_R,ad_T,ad_R1,ad_T1,ad_Rst,ad_Tst,ad_R1st,&
         ad_T1st,ad_Ri,ad_R1i,ad_Di,ad_D1i,ad_Rflux,ad_Tflux,&
         ad_URU,ad_UTU,ad_Riflux,ad_Tiflux,ad_Diflux,ad_URUi,&
         ad_UDUi,ad_UTUi,ad_RL,ad_R1L,ad_DL,ad_D1L,ad_RLflux,&
         ad_DLflux,ad_URUL,ad_SvR)

    IF (verbose) THEN
       WRITE (*,*)                          '  '
       WRITE (*,* )                         '  (sub. call_ad) finished to run ADDOUB'
       WRITE (*,*)                          '  '
    ENDIF

  END SUBROUTINE CALL_AD

END MODULE MCALL_AD
