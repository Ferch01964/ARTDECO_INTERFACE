
MODULE MCALL_OD

  IMPLICIT NONE

  PRIVATE 

  PUBLIC :: CALL_OD

CONTAINS

  SUBROUTINE CALL_OD

    USE MCONSTANTS, only : dp, &
         deg2rad, &
         tiniestdp, &
         undef_dp, &
         plkavg_res
    !!! Line added on 25/03/14 for Fbeam value
    USE MCOMMON, only : verbose, &
         F0, &
         Svin, &
         wlambda, &
         surface_type, &
         surface_family, &
         surface_albedo, &
         nmat, &
         nvza, &
         nvza_up, &
         nvza_dw, &
         ind_vza_up, &
         ind_vza_dw, &
         nvaa, &
         nsza, &
         vza, &
         vaa, &
         sza, &
         t_atm, &
         alt_atm, &
         nlambda, &
         thermal, &
         surface_temp, &
         SvR, &         
         SvR_DW, &
         print_down_rad,&
         SvR_sig, &
         od_accur, &
         od_nstr, &
         rt_cputime, &
         mode, &
         kdis_wvlband, &
         warning, &
         flux_out, &    ! flux_out(nsza,3,nalt_atm,nlambda) down/up/directional flux for given altitude (W m-2) and for each wavelength band
         nalt_atm, &
         lambda_bound, &
         nlayers, &
         n_aitaui_total, &
         nmax_aitaui_mono_layers, &
         n_aitaui_mono_layers, & 
         ai_mono_layers,&        ! ai_mono_layers(nmax_aitaui_mono_layers) 
         tot_dtau_mono_layers, & ! tot_dtau_mono_layers(nlayers, nmax_aitaui_mono_layers) 
         ssa_mono_layers, &      ! ssa_mono_layers(nlayers, nmax_aitaui_mono_layers) 
         nbetal_layers, &
         betal_layers,&       ! betal_layers(nlayers,6,0:nmaxbetal,nlambda) 
         nmaxbetal,&
         lambda_ikdis, &
         sinsca_corint, &
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
         SvR_notms, &    ! SvR_notms(nsza, nvza, nvaa, nnmat, nlambda)
         SvR_sig_notms,&   ! SvR_sig_notms(nsza, nvza, nvaa, nnmat, nlambda) 
         flux_only, &
         varsol_fact, &
         od_no_check,&
         print_down_rad
         !,&
         !Fbeam_user_flag, &
         !Fbeam_user_value


    USE MLAYERS,  only : GET_MONO_LAYERS
    USE MUTILITY, only : REVERSE
    USE MSINSCA,  only : SINSCA

    IMPLICIT NONE

    ! DISORT inputs
    !----------
    ! array dimension variables that was PARAMETERS in DISORT previously
    ! and are now arguments for optimization
    INTEGER :: od_mxcly
    INTEGER :: od_mxulv
    INTEGER :: od_mxcmu
    INTEGER :: od_mxumu
    INTEGER :: od_mxphi
    INTEGER :: od_mi
    INTEGER :: od_mi9m2
    INTEGER :: od_nnlyri
    !----------
    INTEGER :: od_nly
    REAL(kind=dp), allocatable :: od_dtauc(:)
    REAL(kind=dp), allocatable :: od_ssalb(:)
    INTEGER :: od_nmom
    REAL(kind=dp), allocatable :: od_pmom(:, :)
    REAL(kind=dp), allocatable :: od_temper(:)
    REAL(kind=dp) :: od_wvnmlo
    REAL(kind=dp) :: od_wvnmhi
    LOGICAL :: od_usrtau
    INTEGER :: od_ntau
    REAL(kind=dp), allocatable :: od_utau(:)
    LOGICAL :: od_usrang
    INTEGER :: od_numu
    REAL(kind=dp), allocatable :: od_umu(:)
    INTEGER :: od_nphi
    REAL(kind=dp), allocatable :: od_phi(:)
    INTEGER :: od_ibcnd
    REAL(kind=dp) :: od_fbeam
    REAL(kind=dp) :: od_umu0 
    REAL(kind=dp) :: od_phi0 
    REAL(kind=dp) :: od_fisot
    LOGICAL :: od_lamber
    REAL(kind=dp) :: od_surface_albedo
    REAL(kind=dp) :: od_btemp
    REAL(kind=dp) :: od_ttemp
    REAL(kind=dp) :: od_temis
    LOGICAL :: od_planck
    LOGICAL :: od_onlyfl
    LOGICAL :: od_prnt(5)
    CHARACTER(LEN=127) :: od_header
    INTEGER :: od_maxcly
    INTEGER :: od_maxulv
    INTEGER :: od_maxumu
    INTEGER :: od_maxphi
    INTEGER :: od_maxmom
    !DISORT outputs
    REAL(kind=dp), allocatable :: od_rfldir(:)
    REAL(kind=dp), allocatable :: od_rfldn(:)
    REAL(kind=dp), allocatable :: od_flup(:)
    REAL(kind=dp), allocatable :: od_dfdt(:)
    REAL(kind=dp), allocatable :: od_uavg(:)
    REAL(kind=dp), allocatable :: od_uu(:,:,:)
    REAL(kind=dp), allocatable :: od_albmed(:)
    REAL(kind=dp), allocatable :: od_trnmed(:)

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

    INTEGER :: i, u, k, j
    INTEGER :: iai
    INTEGER :: ialt

    !    INTEGER :: iband

    ! time
    REAL :: tcpu1
    REAL :: tcpu2

    
    !-------------------
    ! NOTE : od_nstr, od_accur and were read in a file

    if (od_nstr.lt.2) then
       WRITE(*,*) ' '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' (in sub. call_od) : ERROR   '
       WRITE(*,*) ' Nstream must be > 2         '
       WRITE(*,*) ' Change it in od_spec.dat    '
       WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,*) ' '
       STOP
    endif

    IF (verbose) THEN
       WRITE (*,*) ''
       WRITE (*,*) '  (sub. call_od) set DISORT arguments'
       WRITE (*,*) ''
       WRITE (*,FMT='(A,1x,I6,1x,A)') '                  DISORT will be called ', n_aitaui_total, ' times'
       if (mode.eq.'kdis') &
            WRITE (*,FMT='(A,1x,I6,A)') '                  (maximum number of call per wavelength is ', nmax_aitaui_mono_layers,')'
    ENDIF

    !------- ALLOCATE ARRAYS WHOSE DIMENSION DOES NOT DEPEND ON WAVELENGTH
    !----------
    ! array dimension variables that was PARAMETERS in DISORT previously
    ! and are now arguments for optimization
    od_mxcly  = nlayers   ! = Max no. of computational layers
    od_mxulv  = nalt_atm  ! = Max no. of output levels
    od_mxcmu  = od_nstr   ! = Max no. of computation polar angles
    od_mxumu  = nvza      ! = Max no. of output polar angles
    !od_mxumu  =  od_nstr

    od_mxphi  = nvaa      ! = Max no. of output azimuthal angles
    od_mi     = od_mxcmu / 2        ! as was defined in DISORT
    od_mi9m2  = 9 * od_mi - 2       ! as was defined in DISORT
    od_nnlyri = od_mxcmu * od_mxcly ! as was defined in DISORT
    !----------
    allocate(od_dtauc(od_mxcly), od_ssalb(od_mxcly),  &
         od_temper(0:od_mxcly),od_umu(od_mxumu),od_albmed(od_mxumu),od_trnmed(od_mxumu), &
         od_phi(od_mxphi),                                                               &
         od_utau(od_mxulv),od_rfldir(od_mxulv),od_rfldn(od_mxulv),od_flup(od_mxulv),     &
         od_dfdt(od_mxulv),od_uavg(od_mxulv),od_uu(od_mxumu, od_mxulv, od_mxphi))
    !----------
    ! we use nmaxbetal to dimension the array but we will send
    ! od_pmom(0:od_maxmom, :) array to DISORT
    ! We must dimension this array once for all here 
    ! cause we can not allocate it multiple times
    ! with the right number of Nmom depending on the 
    ! wavelength when perfoming parallel computing.
    allocate(od_pmom(0:nmaxbetal, od_mxcly))

    !-------- IMPLEMENT VARIABLES THAT DOES NOT DEPEND ON WAVELENGTH
    od_maxcly = nlayers
    od_nly    = nlayers
    !-----------
    ! we only want radiant quantities at the top of the atmopshere
    od_maxulv  = nalt_atm
    od_usrtau  = .true.
    od_ntau    = nalt_atm
    od_utau(1) = 0.0_dp 
    

    !-----------
    ! geometry
    od_usrang = .true.
    !od_usrang = .false.
    od_numu   = nvza
    od_maxumu = nvza
    !od_maxumu = od_nstr

    
    
    ! od_umu : cosines of user polar angles : must be sorted in increasing order
    ! (so angles in decreasing number between 90. and 0.0 degrees)
    ! It must also be different from 0
    
    DO i = 1, od_numu
       od_umu(i) = COS(vza(i)*deg2rad)     
      ! IF (od_umu(i) .lt. tiniestdp) THEN
     !     WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     !     WRITE(*,*) ' (in sub. call_od)  : ERROR '
     !     WRITE(*,*) ' Cosine of view polar angle can not be 0 for DISORT'
     !     WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     !     WRITE(*,*) ''
     !     STOP 
     !  ENDIF
       IF (i .gt. 1) THEN
          IF (od_umu(i) .le. od_umu(i-1)) THEN
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ' (in sub. call_od)  : ERROR '
             WRITE(*,*) ' Cosine of view polar angle must be sorted in increasing order for DISORT'
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ''
             STOP 
          ENDIF
       ENDIF      

    ENDDO        


    od_maxphi = nvaa
    od_nphi   = nvaa
    DO i= 1, od_nphi
       od_phi(i) = vaa(i)
    ENDDO
    ! we always consider relative azimuth angles --> Phi_solar = 0
    od_phi0 = 0.0_dp

    ! flag variable (see DISORT doc)
    od_ibcnd = 0 
    

    od_prnt(:) = .FALSE.

    if (thermal .eqv. .true.) then
       od_planck    = .true.
       ! atmosphere temperature profile
       do i = 0, od_nly
          od_temper(i) = t_atm(i+1)
       enddo
       ! surface temperature
       od_btemp = surface_temp
       ! if (mode.ne.'kdis') then
       !    WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       !    WRITE(*,*) ' (in sub. call_od)  : ERROR '
       !    WRITE(*,*) ' od_wvnmlo and od_wvnmhi are not defined if '
       !    WRITE(*,*) ' mode.ne."kdis" '
       !    WRITE(*,*) ' It is needed if thermal emission is ON '
       !    WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       !    WRITE(*,*) ''
       !    STOP 
       ! endif
    else
       od_planck = .false.
    endif

    !-----------
    ! unused inputs
    od_fisot     = 0.0_dp
    od_ttemp     = 0.0_dp
    od_temis     = 0.0_dp
    od_onlyfl    = flux_only
    od_header    = '----------------------------------------------------------&
         &-------------------------------------------------------------------'

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

    IF (verbose)  WRITE (*,*) ''

    !------ START LOOPING OVER LAMBDA
    DO j = 1, nlambda

       CALL GET_MONO_LAYERS(j)

       ! for BB thermal emission
       if (mode.eq.'kdis') then
          od_wvnmlo = 1.0_dp / (kdis_wvlband(3, lambda_ikdis(j)) * 1.d-4 ) 
          od_wvnmhi = 1.0_dp / (kdis_wvlband(2, lambda_ikdis(j)) * 1.d-4 )
       else
          od_wvnmlo = 1.0_dp / (wlambda(j) * (1.0_dp + plkavg_res) * 1.d-4 )  
          od_wvnmhi = 1.0_dp / (wlambda(j) * (1.0_dp - plkavg_res) * 1.d-4 ) 
       endif

       od_maxmom = nbetal_layers(j)
       od_nmom   = nbetal_layers(j)

       IF (od_nstr .gt. od_nmom) THEN      
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' (in sub. call_od)  : ERROR '
          WRITE(*,*) ' od_nmom must be .gt. od_nstr'
          WRITE(*,*) ' od_nmom = ', od_nmom
          WRITE(*,*) ' od_nstr = ', od_nstr
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ''
          STOP 
       ENDIF

      ! write(*,*) 'od_nmom,od_nly ',od_nmom,od_nly
       DO i = 1, od_nly
          DO u = 0, od_nmom
             ! PMOM = alpha1_u / ( 2u+1) 
       !      write(*,*) 'betal_layers(i, 1, u,j)',betal_layers(i, 1, u,j)
             od_pmom(u, i) = betal_layers(i, 1, u,j) / (2.0_dp * u + 1.0_dp) 
            ! od_pmom(u, i) = 0.86**u
          ENDDO
       ENDDO

       ! Surface 
       SELECT CASE(surface_family)
       CASE('lambert')
          od_lamber = .true.
          od_surface_albedo = surface_albedo(j)
       CASE('brdf')
          od_lamber = .false.
          od_surface_albedo = 0.0_dp
       END SELECT

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

       IF (verbose) then
          if (mode.eq.'kdis') then
             WRITE (*, FMT='(1x,A,2x,F10.5,1x,A, 1x, I5)')  &
                  &'                 run DISORT for lambda=',wlambda(j),' microns with nai=', n_aitaui_mono_layers
          else
             WRITE (*, FMT='(1x,A,2x,F10.5,1x,A)')  &
                  &'                 run DISORT for lambda=',wlambda(j),' microns'  
          endif
       ENDIF

       CALL CPU_TIME(tcpu1)

       SvR(:,:,:,:,j)    = 0.0_dp  
       IF( print_down_rad ) THEN
          SvR_DW(:,:,:,:,j)    = 0.0_dp
       END IF
       flux_out(:,:,:,j) = 0.0_dp

       !iband = 76
       !if (lambda_ikdis(j).eq.iband)then
       !print*, wlambda(j), kdis_wvlband(1, lambda_ikdis(j))
       !end if

       do iai = 1, n_aitaui_mono_layers
          !print*, iai
          od_utau(1) = 0.0_dp
          DO i = 1, od_nly
             od_dtauc(i)  = tot_dtau_mono_layers(i,iai)
             od_ssalb(i)  = ssa_mono_layers(i,iai)
             ! this condtion was taken from GAME 
             ! DISORT may return strange result for very low optical depth
             if (od_dtauc(i).le.5.0D-6) then
                od_dtauc(i) = 0.0d0
                od_ssalb(i) = 0.0d0
                
             end if
             od_utau(i+1) = od_utau(i) + od_dtauc(i)
            ! od_utau(i+1) = 0.5
            ! od_ssalb(i)  = 0.9
            ! od_dtauc(i) = 0.5
             !write(*,*)  'ssa_mono_layers,tot_dtau_mono_layers,od_dtauc,od_ssalb ',ssa_mono_layers(i,iai), &
             !& tot_dtau_mono_layers(i,iai),od_dtauc(i),od_ssalb(i)

             !print*, od_ssalb(i), od_dtauc(i)  
            !print*,od_utau(i)
          ENDDO

          ! for the TMS correction
          if (sinsca_corint) then
             ! init result arrays
             sinsca_Stoke_add(:,:,:) = undef_dp 
             sinsca_Stoke_sub(:,:,:) = undef_dp 
             sinsca_tau_gas(:)   = gas_dtau_mono_layers(:,iai)
             sinsca_ssa_gas(:)   = gas_ssa_mono_layers(:,iai)
          endif

          ! if (lambda_ikdis(j).eq.iband)then
          !    print*, ai_mono_layers(iai)
          !    DO i = od_nly, 1, -1
          !       print*, od_ssalb(i),od_dtauc(i)
          !    END DO
          ! end if

          DO i = 1, nsza ! loop on different solar polar angles

             od_umu0 = COS(sza(i)*deg2rad)

             ! Flux in the incident direction 
             od_fbeam = F0(j) * varsol_fact(i)

             !!!! Line added on 25/03/14 for Fbeam
             !if ( Fbeam_user_flag ) then
             !    od_fbeam = Fbeam_user_value
             !else 
             !   od_fbeam = F0(j) * varsol_fact(i)
             !endif

             !print *, 'artdeco_call_od.f90) line 470  ------>  od_fbeam :  ',  od_fbeam
             !print*,  od_utau(od_mxulv)
!             print*, 'od_ssalb', od_ssalb
!             print *, 'od_mxcly, od_mxulv,od_mxcmu,od_maxcly,od_maxulv,od_maxumu,od_maxphi,od_maxmom',&
!             od_mxcly,od_mxulv, od_mxcmu,od_maxcly,od_maxulv,od_maxumu,od_maxphi,od_maxmom
!             print *, 'od_mxumu, od_mxphi ', od_mxumu,od_mxphi
!             print *, 'od_mi',od_mi
!             print *,'od_mi9m2',od_mi9m2
!             print *,'od_nnlyri,od_nly', od_nnlyri,od_nly
!             print *,'od_dtauc',od_dtauc
!             print *, 'n_aitaui_mono_layers', n_aitaui_mono_layers
!             print *,'od_ssalb',od_ssalb
!             print *,'od_nmom',od_nmom
!              print *,'od_pmom', od_pmom
!             print *,'od_temper', od_temper
!             print *,'od_wvnmlo',od_wvnmlo
!             print *,'od_wvnmhi',od_wvnmhi
!             print *,'od_usrtau,od_usrang',od_usrtau,od_usrang
!             print *,'od_ntau,od_numu,od_nphi',od_ntau,od_numu,od_nphi
!              print *,'od_utau',od_utau
!             print *,'od_nstr',od_nstr
!             print *,'od_umu',od_umu
!             print *,'od_phi',od_phi
!             print *,'od_ibcnd',od_ibcnd
!             print *,'od_wvnmlo,od_fbeam',od_wvnmlo,od_fbeam
!             print *,'od_umu0',od_umu0
!             print *,'od_phi0',od_phi0
!             print *,'od_fisot',od_fisot
!             print *,'od_lamber',od_lamber
!             print *,'od_surface_albedo',od_surface_albedo
!             print *,'od_btemp',od_btemp
!             print *,'od_ttemp',od_ttemp
!             print *,'od_temis',od_temis
!             print *,'od_planck',od_planck
!             print *,'od_onlyfl',od_onlyfl
!             print *,'od_accur',od_accur
!             print *,'od_rfldir', od_rfldir

!             print *,'od_rfldn',od_rfldn
!             print *,'od_flup',od_flup
!             print *,'od_dfdt',od_dfdt
!             print *,'od_uavg',od_uavg


             CALL DISORT(od_no_check, warning, od_mxcly, od_mxulv, od_mxcmu, od_mxumu, od_mxphi, &
                  od_mi, od_mi9m2, od_nnlyri,                                             &
                  od_nly, od_dtauc, od_ssalb, od_nmom, od_pmom(0:od_maxmom,:), od_temper, &
                  od_wvnmlo, od_wvnmhi, od_usrtau, od_ntau, od_utau, od_nstr,             &
                  od_usrang, od_numu, od_umu, od_nphi, od_phi, od_ibcnd, od_fbeam,        &
                  od_umu0, od_phi0, od_fisot, od_lamber, od_surface_albedo, od_btemp,     &
                  od_ttemp, od_temis, od_planck, od_onlyfl, od_accur, od_prnt,            &
                  od_header, od_maxcly, od_maxulv, od_maxumu, od_maxphi,                  &
                  od_maxmom, od_rfldir, od_rfldn, od_flup, od_dfdt, od_uavg, od_uu,       &
                  od_albmed, od_trnmed ) 
             ! if (od_utau(1).ne.0.0_dp) then
             !    print*,  od_utau  
             !    call sleep(2)
             ! endif
              
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

                ! implement the result array 
                DO u = 1, nvza
                   DO k = 1, nvaa
                      ! od_uu(u,1,k) --> 1 means utau(1) = TOA
                      SvR(i,u,k,1,j) = SvR(i,u,k,1,j) + &
                           (od_uu(u,1,k) - sinsca_Stoke_sub(1, u, k) + sinsca_Stoke_add(1, u, k) ) * ai_mono_layers(iai)
                      SvR_notms(i,u,k,1,j) = SvR_notms(i,u,k,1,j) + od_uu(u,1,k)  * ai_mono_layers(iai)
                   ENDDO
                ENDDO

             else

              

               
                !DO u = 1, nvza
                !   DO k = 1, nvaa
                !      ! od_uu(u,1,k) --> 1 means utau(1) = TOA 
                !      SvR(i,u,k,1,j) = SvR(i,u,k,1,j) + od_uu(u,1,k) * ai_mono_layers(iai)                     
                !   ENDDO
                !ENDDO

                DO k = 1, nvaa
                  IF (print_down_rad) THEN
                    DO u = 1, nvza_dw
                     SvR_DW(i,u,k,1,j) = SvR_DW(i,u,k,1,j) + od_uu(ind_vza_dw(u),od_maxulv,k) * ai_mono_layers(iai)
                    ! write(*,*) 'UMU, SvR_DW',od_umu(ind_vza_dw(u)) ,  SvR_DW(i,u,k,1,j)
                     END DO !!! u
                  END IF
                  
                  IF (nvza_up .gt. 0) THEN
                    DO u = 1, nvza_up
                      SvR(i,u,k,1,j) = SvR(i,u,k,1,j) + od_uu(ind_vza_up(u),1,k) * ai_mono_layers(iai)
                    END DO  !!! u 
                  END IF

                END DO !!! k
             endif
                          ! COMPUTE FLUXES
             DO ialt = 1, nalt_atm 
                !write(*,*) 'od_rfldn,od_rfldir,od_flup',od_rfldn(ialt),od_rfldir(ialt),od_flup(ialt)
                flux_out(i,1,ialt,j) = flux_out(i,1,ialt,j) + ai_mono_layers(iai) * (od_rfldn(ialt)+od_rfldir(ialt))
                !write(*,*)  'ialt, flux_dw', ialt, flux_out(i,1,ialt,j)
                flux_out(i,2,ialt,j) = flux_out(i,2,ialt,j) + ai_mono_layers(iai) * od_flup(ialt)
                flux_out(i,3,ialt,j) = flux_out(i,3,ialt,j) + ai_mono_layers(iai) * od_rfldir(ialt)
             ENDDO
          ENDDO ! solar angles

       end do
       


       CALL CPU_TIME(tcpu2)
       rt_cputime(j) = tcpu2-tcpu1
       !IF (verbose) THEN
       !   WRITE (*,FMT='(1x,A30,1PE14.5, A9)') '  (sub. call_od) CPU time     ', (tcpu2-tcpu1), ' sec(s) -'
       !ENDIF

       if (sinsca_corint) then
          DEALLOCATE(sinsca_mu, &
               sinsca_p11_ptcle_sub, &
               sinsca_p12_ptcle_sub,&
               sinsca_p11_ptcle_add, &
               sinsca_p12_ptcle_add)
       endif

    ENDDO! lambda

    SvR_sig(:,:,:,:,:)  = 0.0_dp
    if (sinsca_corint) SvR_sig_notms(:,:,:,:,:)  = 0.0_dp

    deallocate(od_dtauc, od_ssalb,od_pmom,     &
         od_temper,od_umu,od_albmed,od_trnmed, &
         od_phi,                               &
         od_utau,od_rfldir,od_rfldn,od_flup,   &
         od_dfdt,od_uavg,od_uu)

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
       WRITE (*,*)                          '  '
       WRITE (*,* )                         '  (sub. call_od) finished to run DISORT'
       WRITE (*,*)                          '  '
    ENDIF

  END SUBROUTINE CALL_OD

END MODULE MCALL_OD
