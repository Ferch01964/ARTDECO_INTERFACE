
MODULE MCALL_DOAD

  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: CALL_DOAD

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE CALL_DOAD

    USE MCONSTANTS, only : dp, deg2rad, undef_dp, xpi, plkavg_res
    !!! line added on 25/03/14 for fbeam
    USE MCOMMON, only : verbose, &
         doad_eps, &
         doad_nmug, &
         doad_nfoumx,  &
         doad_ncoef_min_bdrf, &
         surface_albedo, &
         surface_type, &
         surface_family, &
         nlambda, &
         nmat, &
         svin, &
         F0, &
         nalt_atm, &
         nlayers, &
         nbetal_layers, &
         betal_layers, &
         tot_dtau_mono_layers, &
         ssa_mono_layers, &
         n_aitaui_mono_layers,&
         ai_mono_layers,&
         n_aitaui_total, &
         nmax_aitaui_mono_layers, &
         nsza, &
         nvza, &
         nvaa, &
         sza, &
         vza, &
         vaa, &
         thermal, &
         t_atm, &
         surface_temp, &
         wlambda, &
         SvR, &
         SvR_DW, print_down_rad, &
         SvR_sig,&
         flux_out, &    ! flux_out(nsza,3,nalt_atm,nlambda) down/up/directional flux for given altitude (W m-2) and for each wavelength band
         rt_cputime, &
         mode, &
         ptcle_ssa_layers, &           ! ptcle_ssa_layers(nlayers,nlambda)   
         ptcle_phasemat_layers, &      ! ptcle_phasemat_layers(nlayers, 6, nang_max,nlambda)
         ptcle_ssa_layers_unscaled, &  ! ptcle_ssa_layers_unscaled(nlayers,nlambda)   
         ptcle_phasemat_tms_layers, &  ! ptcle_phasemat_tms_layers(nlayers, 6, nang_max,nlambda)
         nptcle,&
         ptcle_nang_phasemat_layers, &
         ptcle_u_phasemat_layers, & 
         ptcle_dtau_layers, &
         gas_dtau_mono_layers, &
         gas_ssa_mono_layers, &
         depol, &
         SvR_notms, &     ! SvR_notms(nsza, nvza, nvaa, nnmat, nlambda)
         SvR_sig_notms, & ! SvR_sig_notms(nsza, nvza, nvaa, nnmat, nlambda) 
         sinsca_corint, &
         lambda_ikdis, &
         kdis_wvlband, &
         varsol_fact
         !,&
         !Fbeam_user_flag,&
         !Fbeam_user_value

    USE MLAYERS, only : GET_MONO_LAYERS
    USE MUTILITY, only : REVERSE
    USE MSINSCA, only : SINSCA

    IMPLICIT NONE
    !!! Line added on 25/03
    !REAL(kind=dp) :: doad_fbeam
 
    ! doad arguments

    logical :: doad_verbose
    integer :: doad_ndmu
    integer :: doad_ndmuext
    integer :: doad_ndsup
    integer :: doad_ndcoef
    integer :: doad_ndmu0
    integer :: doad_ndlay
    integer :: doad_ndirf
    integer :: doad_ndphiext

    integer :: doad_nlayer
    integer :: doad_nmat
    integer :: doad_igrnd
    real(kind=dp) :: doad_svin(4)
    integer :: doad_nmuext
    real(kind=dp), allocatable :: doad_xmuext(:)
    integer :: doad_nphout
    real(kind=dp), allocatable :: doad_phout(:)
    integer :: doad_nmu0
    real(kind=dp), allocatable :: doad_xmu0(:)
    real(kind=dp), allocatable :: doad_b(:)
    real(kind=dp), allocatable :: doad_ssa(:)
    real(kind=dp), allocatable :: doad_alfbet(:,:,:,:)
    integer, allocatable :: doad_ncoef(:)
    real(kind=dp) :: doad_surf_alb
    integer :: doad_irf
    REAL(kind=dp) :: doad_wvnmlo
    REAL(kind=dp) :: doad_wvnmhi
    real(kind=dp), allocatable :: doad_tlev(:)
    real(kind=dp) :: doad_tsurf

    real(kind=dp), allocatable :: doad_SvR(:,:,:,:)
    real(kind=dp), allocatable :: doad_SvR_DW(:,:,:,:)

    real(kind=dp), allocatable :: doad_flux(:,:,:)
    real(kind=dp), allocatable :: doad_SvR_th(:,:)
    real(kind=dp), allocatable :: doad_flux_th(:,:)

    ! those variables are used for the TMS correction
    INTEGER :: sinsca_nmu
    REAL(kind=dp), ALLOCATABLE :: sinsca_mu(:)

    REAL(kind=dp), ALLOCATABLE :: sinsca_tau_ptcle(:)

    REAL(kind=dp), ALLOCATABLE :: sinsca_p11_ptcle_add(:,:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_p12_ptcle_add(:,:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_ssa_ptcle_add(:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_p11_ptcle_sub(:,:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_p12_ptcle_sub(:,:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_ssa_ptcle_sub(:)

    REAL(kind=dp), ALLOCATABLE :: sinsca_tau_gas(:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_ssa_gas(:)

    REAL(kind=dp), ALLOCATABLE :: sinsca_Stoke_add(:,:,:)
    REAL(kind=dp), ALLOCATABLE :: sinsca_Stoke_sub(:,:,:)

    REAL(kind=dp) :: sinsca_Stoke_in(nmat)

    INTEGER :: i, j, u, k
    integer :: ialt
    integer :: iai

    ! time
    REAL :: tcpu1
    REAL :: tcpu2

    EXTERNAL doad

    ! ====================
    !doad_nmug = 2*doad_nmug
    if (doad_nmug.lt.2) then
       WRITE(*,*) ' '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. call_doad) : ERROR   '
       WRITE(*,*) ' Nstream must be > 2         '
       WRITE(*,*) ' Change it in doad_spec.dat    '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' '
       STOP
    endif

    !doad_verbose = .false.
    doad_verbose = .true.

    ! NOTE : doad_eps, doad_nmug and doad_nfoumax were read in a file
    IF (verbose) THEN
       WRITE (*,*) ''
       WRITE (*, *)'  (sub. call_ad) set DOAD arguments'
       WRITE (*,*) ''
       WRITE (*,FMT='(A,1x,I5,1x,A)') '                  DOAD will be called ',n_aitaui_total, ' times'
       if (mode.eq.'kdis') &
            WRITE (*,FMT='(A,1x,I6,A)') '                  (maximum number of call per wavelength is ', nmax_aitaui_mono_layers,')'
    ENDIF

    ! ARRAY DIMENSIONS
    doad_ndmuext   = nvza + nsza    ! we add solar zenith angles cause it
    !                                 could be added to muext if no
    !                                 muext = mu0
    doad_ndmu      = doad_ndmuext + doad_nmug
    !doad_ndmu      = doad_ndmuext + 2* doad_nmug

    doad_ndsup     = nmat * doad_ndmu        
    doad_ndmu0     = nsza
    doad_ndlay     = nlayers
    doad_ndphiext  = nvaa

    doad_ndirf   = 1 ! must be > 0 cause it is an array size

    !! line added on 30/04/15 for doad_svr_dw
    ALLOCATE(doad_xmuext(doad_ndmuext), &
         doad_phout(doad_ndphiext), &
         doad_xmu0(doad_ndmu0), &
         doad_b(0:doad_ndlay), &
         doad_ssa(doad_ndlay), &
         doad_ncoef(doad_ndlay), &
         doad_tlev(0:doad_ndlay), &
         doad_SvR(nmat,doad_NDmu0,doad_Ndmuext-doad_NDmu0,doad_NDphiext), &
         doad_SvR_DW(nmat,doad_NDmu0,doad_Ndmuext-doad_NDmu0,doad_NDphiext), &
         doad_flux(doad_NDmu0, 3, doad_ndlay+1), &
         doad_SvR_th(nmat,doad_Ndmuext-doad_NDmu0), &
         doad_flux_th(2, doad_ndlay+1) )

    doad_SvR(:,:,:,:) = undef_dp
    doad_SvR_DW(:,:,:,:) = undef_dp
    doad_flux(:,:,:)  = undef_dp
    doad_SvR_th(:,:)  = undef_dp
    doad_flux_th(:,:) = undef_dp

    doad_nlayer = nlayers
    doad_nmat   = nmat

    ! Geometry
    doad_nmuext = nvza
    doad_xmuext(1:doad_nmuext) = cos(vza(:)*deg2rad)
    
    doad_nphout = nvaa
    ! in doad phi = phiv - phiext
    
    doad_phout(1:doad_nphout) = vaa(:)
    doad_nmu0   = nsza
    doad_xmu0(1:doad_nmu0) = cos(sza(:)*deg2rad) 

    ! type of surface
    if (surface_family .eq. 'brdf') then
       doad_igrnd    = 3
       doad_surf_alb = undef_dp
    endif

    if (thermal .eqv. .true.) then

       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. call_doad)  : ERROR      '
       WRITE(*,*) ' Thermal source emission           '
       WRITE(*,*) ' is not yet implemented in DOAD    '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ''
       STOP 

       doad_irf = 1
       do i = 0, doad_nlayer
          doad_tlev(i) = t_atm(doad_nlayer-i+1)
          write(*,*) i, doad_tlev(i)
       enddo
       ! surface temperature
       doad_tsurf = surface_temp
       ! if (mode.ne.'kdis') then
       !    WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       !    WRITE(*,*) ' (in sub. call_doad)  : ERROR '
       !    WRITE(*,*) ' doad_wvnmlo and doad_wvnmhi are not defined if '
       !    WRITE(*,*) ' mode.ne."kdis" '
       !    WRITE(*,*) ' It is needed if thermal emission is ON '
       !    WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       !    WRITE(*,*) ''
       !    STOP 
       ! endif
    else
       ! no thermal emission
       doad_irf      = 0
       doad_tlev(:)  = 0.0_dp
       doad_tsurf    = 0.0_dp
    endif

    ! for the TMS correction
    IF (sinsca_corint) then
       ALLOCATE(sinsca_tau_gas(nlayers), &
            sinsca_ssa_gas(nlayers),&
            sinsca_tau_ptcle(nlayers), &
            sinsca_ssa_ptcle_sub(nlayers),&
            sinsca_ssa_ptcle_add(nlayers), &
            sinsca_Stoke_add(nmat, nvza, nvaa), &
            sinsca_Stoke_sub(nmat, nvza, nvaa))
       sinsca_tau_ptcle(:)     = 0.0_dp
       sinsca_ssa_ptcle_sub(:) = 0.0_dp
       sinsca_ssa_ptcle_add(:) = 0.0_dp
    end IF

    IF (verbose) THEN
       WRITE (*,*) ''
    ENDIF

    DO j = 1, nlambda

       CALL GET_MONO_LAYERS(j)

       ! for BB thermal emission
       if (mode.eq.'kdis') then
          doad_wvnmlo = 1.0_dp / (kdis_wvlband(3, lambda_ikdis(j)) * 1.d-4 ) 
          doad_wvnmhi = 1.0_dp / (kdis_wvlband(2, lambda_ikdis(j)) * 1.d-4 )
       else
          doad_wvnmlo = 1.0_dp / (wlambda(j) * (1.0_dp + plkavg_res) * 1.d-4 )  
          doad_wvnmhi = 1.0_dp / (wlambda(j) * (1.0_dp - plkavg_res) * 1.d-4 ) 
       endif

       ! wavelength dependent surface characteristics
       if (surface_family .eq. 'lambert') then
          if ((surface_albedo(j).eq.0.0_dp).and.(thermal.eqv..false.)) then
             doad_igrnd    = 0
             doad_surf_alb = undef_dp
          else
             doad_igrnd = 1 
             doad_surf_alb = surface_albedo(j)     
          end if
       end if

       ! we set the array dimensioning value doad_ndcoef
       if (nbetal_layers(j).lt.doad_ncoef_min_bdrf) then
          doad_ndcoef = doad_ncoef_min_bdrf
       else
          doad_ndcoef = nbetal_layers(j)
       endif

       ALLOCATE(doad_alfbet(4,4,0:doad_ndcoef,doad_ndlay))

       doad_b(:)            = 0.0_dp
       doad_ssa(:)          = 0.0_dp
       doad_alfbet(:,:,:,:) = 0.0_dp
       doad_ncoef(:)        = 0       

       DO i = 1, doad_nlayer
          ! layers must be sorted in ascending order for adding
          doad_ncoef(i) = nbetal_layers(j)
          DO u = 0, doad_ncoef(i)
             doad_alfbet(1,1,u,i) = betal_layers(doad_nlayer+1-i, 1, u, j)    ! alpha1
             doad_alfbet(1,2,u,i) = betal_layers(doad_nlayer+1-i, 5, u, j)    ! beta1
             doad_alfbet(2,1,u,i) = betal_layers(doad_nlayer+1-i, 5, u, j)    ! beta1
             doad_alfbet(2,2,u,i) = betal_layers(doad_nlayer+1-i, 2, u, j)    ! alpha2 
             doad_alfbet(3,3,u,i) = betal_layers(doad_nlayer+1-i, 3, u, j)    ! alpha3
             doad_alfbet(3,4,u,i) = betal_layers(doad_nlayer+1-i, 6, u, j)    ! beta2 
             doad_alfbet(4,3,u,i) = -betal_layers(doad_nlayer+1-i, 6,u, j)    ! -beta2
             doad_alfbet(4,4,u,i) = betal_layers(doad_nlayer+1-i, 4, u, j)    ! alpha4 

           ENDDO
       ENDDO

       doad_svin(:) = svin(:) 

       ! for the TMS correction
       if (sinsca_corint) then
          sinsca_nmu = ptcle_nang_phasemat_layers(j)
          ALLOCATE(sinsca_mu(sinsca_nmu),      &
               sinsca_p11_ptcle_sub(nlayers,sinsca_nmu), &
               sinsca_p12_ptcle_sub(nlayers,sinsca_nmu), &
               sinsca_p11_ptcle_add(nlayers,sinsca_nmu), &
               sinsca_p12_ptcle_add(nlayers,sinsca_nmu))
          if(nptcle.eq.0) then
             sinsca_p11_ptcle_add(:,:) = 0.0_dp
             sinsca_p12_ptcle_add(:,:) = 0.0_dp
             sinsca_p11_ptcle_sub(:,:) = 0.0_dp
             sinsca_p12_ptcle_sub(:,:) = 0.0_dp
          endif
          ! we revesre mu for it to be sorted in increasing order
          sinsca_mu(:) = REVERSE(sinsca_nmu, ptcle_u_phasemat_layers(1:sinsca_nmu,j))

          if (nptcle .gt. 0) then
             do i = 1, nlayers
                ! we revert P11 and P12 cause mu was reverted to be 
                ! sorted in increasing order
                sinsca_p11_ptcle_sub(i,:) = REVERSE(sinsca_nmu, ptcle_phasemat_layers(i,1,1:sinsca_nmu,j)) ! P11
                sinsca_p12_ptcle_sub(i,:) = REVERSE(sinsca_nmu, ptcle_phasemat_layers(i,5,1:sinsca_nmu,j)) ! P12 = P21
                sinsca_ssa_ptcle_sub(i)   = ptcle_ssa_layers(i,j)
                sinsca_p11_ptcle_add(i,:) = REVERSE(sinsca_nmu, ptcle_phasemat_tms_layers(i,1,1:sinsca_nmu,j)) ! P11
                sinsca_p12_ptcle_add(i,:) = REVERSE(sinsca_nmu, ptcle_phasemat_tms_layers(i,5,1:sinsca_nmu,j)) ! P12 = P21
                sinsca_ssa_ptcle_add(i)   = ptcle_ssa_layers_unscaled(i,j)
                sinsca_tau_ptcle(i)       = ptcle_dtau_layers(i,j)
             enddo
          endif
       endif

       ! call the adding code
       IF (verbose) then
          if (mode.eq.'kdis') then
             WRITE (*, FMT='(1x,A,2x,F10.5,1x,A, 1x, I4)')  &
                  &'                 run DOAD for lambda=',wlambda(j),' microns with nai=', n_aitaui_mono_layers
          else
             WRITE (*, FMT='(1x,A,2x,F10.5,1x,A)')  &
                  &'                 run DOAD for lambda=',wlambda(j),' microns'  
          endif
       ENDIF

       CALL CPU_TIME(tcpu1)

       SvR(:,:,:,:,j)       = 0.0_dp
       !!! line added on 30/04/15
       IF ( print_down_rad) THEN
          SvR_DW(:,:,:,:,j)       = 0.0_dp
       END IF
       if (sinsca_corint) SvR_notms(:,:,:,:,j) = 0.0_dp
       flux_out(:,:,:,j)    = 0.0_dp

       do iai = 1, n_aitaui_mono_layers

          DO i = 1, doad_nlayer
             doad_ssa(i)   = ssa_mono_layers(nlayers+1-i,iai)
             doad_b(i)     = tot_dtau_mono_layers(nlayers+1-i,iai)            
          END DO


          





!          print *, 'doad_ndmu, doad_ndmuext, doad_ndsup,doad_ndcoef,doad_ndmu0,doad_ndlay,doad_ndirf,doad_ndphiext',&
!           doad_ndmu, doad_ndmuext, doad_ndsup,doad_ndcoef,doad_ndmu0,doad_ndlay,doad_ndirf,doad_ndphiext
!          print *, 'doad_nlayer,doad_nfoumx,doad_nmug,doad_nmat,',doad_nlayer,doad_nfoumx,doad_nmug,doad_nmat
!          print *, 'doad_igrnd',doad_igrnd
!          print *, 'doad_ncoef_min_bdrf',doad_ncoef_min_bdrf
!          print *, 'doad_svin',doad_svin
!          print *, 'doad_eps,doad_nmuext',doad_eps,doad_nmuext
!          print *, 'doad_xmuext',doad_xmuext
!          print *, 'doad_nphout,doad_nmu0',doad_nphout,doad_nmu0
!          print *, 'doad_phout',doad_phout
!          print *, 'doad_b',doad_b
!          print *, 'doad_ssa',doad_ssa
!          print *, 'doad_alfbet',doad_alfbet
!          print *, 'doad_ncoef',doad_ncoef
!          print *, 'doad_surf_alb',doad_surf_alb
!          print *, 'doad_irf',doad_irf
!          print *, 'doad_tlev',doad_tlev
!          print *, 'doad_tsurf',doad_tsurf
!          print *, 'doad_wvnmlo',doad_wvnmlo
!          print *, 'doad_wvnmhi',doad_wvnmhi
          
          !! line "SvR_DW, print_down_rad " added on 30/04/15 
          CALL doad( doad_verbose,&
               doad_ndmu,&
               doad_ndmuext,&
               doad_ndsup,&
               doad_ndcoef,&
               doad_ndmu0,&
               doad_ndlay,&
               doad_ndirf,&
               doad_ndphiext,&
               wlambda(j), &
               doad_nlayer,&
               doad_nfoumx,&
               doad_nmug,&
               doad_nmat,&
               doad_igrnd,&
               doad_ncoef_min_bdrf, &
               doad_svin,&
               doad_eps,&
               doad_nmuext,&
               doad_xmuext,&
               doad_nphout,&
               doad_phout,&
               doad_nmu0,&
               doad_xmu0,&
               doad_b, &
               doad_ssa,& 
               doad_alfbet,& 
               doad_ncoef,&
               doad_surf_alb, &
               doad_irf, &
               doad_tlev, &
               doad_tsurf,&
               doad_wvnmlo, &
               doad_wvnmhi, &
               doad_flux,&
               doad_SvR,&
               doad_SvR_DW, &
               doad_flux_th,&
               doad_SvR_th)

         


          ! for the TMS correction
          if (sinsca_corint) then
             ! init result arrays
             sinsca_Stoke_add(:,:,:) = undef_dp 
             sinsca_Stoke_sub(:,:,:) = undef_dp 
             sinsca_tau_gas(:)   = gas_dtau_mono_layers(:,iai)
             sinsca_ssa_gas(:)   = gas_ssa_mono_layers(:,iai)
          endif

          DO i = 1, nsza

            !if( Fbeam_user_flag) then
            !     doad_fbeam =  Fbeam_user_value
            !else 
            !      doad_fbeam = F0(j) * varsol_fact(i)
            !endif
         !    WRITE(*,*) ' '
         !    WRITE(*,*) '  (sub. artdeco_call_doad) '
         !    WRITE(*,*) ' '
             !WRITE(*,*) '>>>>>>   doad_fbeam : ',doad_fbeam
             
             !doad_flux(i,:,:)  = doad_flux(i,:,:)  * doad_fbeam       
             !doad_SvR(:,i,:,:) = doad_SvR(:,i,:,:) * doad_fbeam / xpi 

             doad_flux(i,:,:)  = doad_flux(i,:,:)  * F0(j) * varsol_fact(i)       
             doad_SvR(:,i,:,:) = doad_SvR(:,i,:,:) * F0(j) * varsol_fact(i) / xpi

             !!! line added on 30/04/15
             IF ( print_down_rad) THEN
                doad_SvR_DW(:,i,:,:) = doad_SvR_DW(:,i,:,:) * F0(j) * varsol_fact(i) / xpi 
             END IF
             
             
             ! for the TMS correction
             if (sinsca_corint) then

                sinsca_Stoke_in(1:nmat) = Svin(1:nmat) * F0(j) * varsol_fact(i)

                CALL SINSCA(verbose, &
                     wlambda(j), &
                     nmat, &
                     sza(i), nvza, vza, nvaa, vaa, &
                     nlayers, &
                     sinsca_nmu, sinsca_mu, &
                     sinsca_p11_ptcle_sub, sinsca_p12_ptcle_sub, &
                     sinsca_tau_ptcle, &
                     sinsca_ssa_ptcle_sub, &
                     sinsca_tau_gas,   &
                     sinsca_ssa_gas,   &
                     depol(j), &
                     sinsca_Stoke_in, &
                     'lambert', &
                     0.0_dp, &
                     sinsca_Stoke_sub)

                CALL SINSCA(verbose, &
                     wlambda(j), &
                     nmat, &
                     sza(i), nvza, vza, nvaa, vaa, &
                     nlayers, &
                     sinsca_nmu, sinsca_mu, &
                     sinsca_p11_ptcle_add, sinsca_p12_ptcle_add, &
                     sinsca_tau_ptcle, &
                     sinsca_ssa_ptcle_add, &
                     sinsca_tau_gas,   &
                     sinsca_ssa_gas,   &
                     depol(j), &
                     sinsca_Stoke_in, &
                     'lambert', &
                     0.0_dp, &
                     sinsca_Stoke_add)

                DO u = 1, nvza
                   DO k = 1, nvaa
                      SvR(i,u,k,1,j) = SvR(i,u,k,1,j) + &
                           (doad_SvR_th(1,u) + doad_SvR(1,i,u,k) - sinsca_Stoke_sub(1,u,k) + sinsca_Stoke_add(1,u,k)) *    &
                           ai_mono_layers(iai)
                      SvR_notms(i,u,k,1,j) = SvR_notms(i,u,k,1,j) + (doad_SvR_th(1,u) + doad_SvR(1,i,u,k) ) * ai_mono_layers(iai)

                                           
                      if (nmat .gt. 1) then
                         SvR(i,u,k,2,j) = SvR(i,u,k,2,j) + &
                              (doad_SvR_th(2,u) + doad_SvR(2,i,u,k) - sinsca_Stoke_sub(2,u,k) + sinsca_Stoke_add(2,u,k)) * &
                              ai_mono_layers(iai)
                         SvR(i,u,k,3,j) = SvR(i,u,k,3,j) + &
                              (doad_SvR_th(3,u) + doad_SvR(3,i,u,k) - sinsca_Stoke_sub(3,u,k) + sinsca_Stoke_add(3,u,k)) * &
                              ai_mono_layers(iai)
                         SvR_notms(i,u,k,2,j) = SvR_notms(i,u,k,2,j) + (doad_SvR_th(2,u) + doad_SvR(2,i,u,k) ) * ai_mono_layers(iai)
                         SvR_notms(i,u,k,3,j) = SvR_notms(i,u,k,3,j) + (doad_SvR_th(3,u) + doad_SvR(3,i,u,k) ) * ai_mono_layers(iai)
                      endif
                      if (nmat .gt. 3) then
                         SvR(i,u,k,4,j) = SvR(i,u,k,4,j) + &
                              (doad_SvR_th(4,u) + doad_SvR(4,i,u,k) - sinsca_Stoke_sub(4,u,k) + sinsca_Stoke_add(4,u,k)) * &
                              ai_mono_layers(iai)
                         SvR_notms(i,u,k,4,j) = SvR_notms(i,u,k,4,j) + (doad_SvR_th(4,u) + doad_SvR(4,i,u,k) ) * ai_mono_layers(iai)
                      endif
                   ENDDO
                ENDDO

             else

                DO u = 1, nvza
                   DO k = 1, nvaa
                      SvR(i,u,k,1,j) = SvR(i,u,k,1,j) + doad_SvR(1,i,u,k) * ai_mono_layers(iai)
                     
                       !!line added on 30/04/15
                      IF ( print_down_rad) THEN
                          !write(*,*) 'i,u,k, doad_SvR_DW',i,u,k,doad_SvR_DW(1,i,u,k)
                          SvR_DW(i,u,k,1,j) = SvR_DW(i,u,k,1,j) + doad_SvR_DW(1,i,u,k) * ai_mono_layers(iai)  
                      END IF

                      if (nmat .gt. 1) then
                         SvR(i,u,k,2,j) = SvR(i,u,k,2,j) + doad_SvR(2,i,u,k) * ai_mono_layers(iai)
                         SvR(i,u,k,3,j) = SvR(i,u,k,3,j) + doad_SvR(3,i,u,k) * ai_mono_layers(iai)
                          !!line added on 30/04/15
                          IF ( print_down_rad) THEN
                             SvR_DW(i,u,k,2,j) = SvR_DW(i,u,k,2,j) + doad_SvR_DW(2,i,u,k) * ai_mono_layers(iai)
                             SvR_DW(i,u,k,3,j) = SvR_DW(i,u,k,3,j) + doad_SvR_DW(3,i,u,k) * ai_mono_layers(iai)
                          END IF

                      endif
                      if (nmat .gt. 3) then
                        SvR(i,u,k,4,j) = SvR(i,u,k,4,j) + doad_SvR(4,i,u,k) * ai_mono_layers(iai)
                        IF ( print_down_rad) THEN
                           SvR_DW(i,u,k,4,j) = SvR_DW(i,u,k,4,j) + doad_SvR_DW(4,i,u,k) * ai_mono_layers(iai)
                        ENDIF
                      end if
                   ENDDO
                ENDDO

             endif

             ! COMPUTE FLUXES
             DO ialt = 1, nalt_atm 
                ! NB: layers ordering of DOAD is inverted regarding ARTDECO
                flux_out(i,1,nalt_atm+1-ialt,j) = flux_out(i,1,nalt_atm+1-ialt,j) + ai_mono_layers(iai) * doad_flux(i,1,ialt)
                flux_out(i,2,nalt_atm+1-ialt,j) = flux_out(i,2,nalt_atm+1-ialt,j) + ai_mono_layers(iai) * doad_flux(i,2,ialt)                  
                flux_out(i,3,nalt_atm+1-ialt,j) = flux_out(i,3,nalt_atm+1-ialt,j) + ai_mono_layers(iai) * doad_flux(i,3,ialt)                  
             ENDDO

          ENDDO ! loop on sza

       end do ! on ai 

       CALL CPU_TIME(tcpu2)
       rt_cputime(j) = tcpu2-tcpu1
       ! IF (verbose) THEN
       !    WRITE (*,FMT='(1x,A38,1PE14.5, A9)') '  (sub. call_doad) CPU time     ', tcpu2-tcpu1, ' sec(s) -'
       ! ENDIF

       DEALLOCATE(doad_alfbet)

       if (sinsca_corint) then
          DEALLOCATE(sinsca_mu, &
               sinsca_p11_ptcle_sub, &
               sinsca_p12_ptcle_sub,&
               sinsca_p11_ptcle_add, &
               sinsca_p12_ptcle_add)
       endif

    END DO ! loop on wavelength

    SvR_sig(:,:,:,:,:) = 0.0_dp
    if (sinsca_corint) SvR_sig_notms(:,:,:,:,:)  = 0.0_dp

    DEALLOCATE(doad_xmuext, doad_phout, doad_xmu0, doad_b, &
         doad_ssa, doad_ncoef, doad_tlev, &
         doad_SvR,doad_SvR_DW )

    if (sinsca_corint) then
       DEALLOCATE(sinsca_tau_gas, &
            sinsca_ssa_gas,&
            sinsca_tau_ptcle, &
            sinsca_ssa_ptcle_add, &
            sinsca_ssa_ptcle_sub,&
            sinsca_Stoke_add,&
            sinsca_Stoke_sub)
    endif

    IF (verbose) THEN
       WRITE (*,* )                         '  (sub. call_doad) finished to run DOAD'
       WRITE (*,*)                          '  '
    ENDIF

  END SUBROUTINE CALL_DOAD

END MODULE MCALL_DOAD
