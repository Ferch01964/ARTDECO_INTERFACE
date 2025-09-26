MODULE MLAYERS

  IMPLICIT NONE

  ! suppress those attribute for f2py
  !  PRIVATE :: TAU_RAYLEIGH, SET_GASABS, get_aitauieff_layers
  !  PUBLIC :: GET_LAYERS, get_mono_layers

CONTAINS

  SUBROUTINE GET_LAYERS

    USE MCONSTANTS , only : dp, &
         deg2rad, &
         eps_betal, &
         min_non0_alt, &
         nalt_distrib, &
         undef_dp, &
         undef_i

    USE MCOMMON , only : verbose, &
         nptcle, &
         nlambda, &
         wlambda, &
         nlayers, &
         mc_nmu, &
         depol, & 
         trunc_method, &
         od_nstr, &
         rt_model, &
         nmaxbetal, &
         ptcle_nbetal, &
         ptcle_betal, &
         nang_max, &
         ptcle_nang_phasemat, &
         ptcle_phasemat, &
         ptcle_trunc_phasemat, &
         ptcle_u_phasemat, &
         ptcle_Cext_reflamb, &
         ptcle_tau_reflamb, &
         ptcle_opt, &
         ptcle_trunccoeff, &
         ptcle_user_vdist_nalt, &
         ptcle_user_vdist_alt, &
         ptcle_user_vdist, &
         ptcle_vdist_type, &
         ptcle_h_vdist, &
         ptcle_sd_gauss_vdist, &
         ptcle_nang_phasemat_layers, &
         ptcle_phasemat_layers, &
         ptcle_u_phasemat_layers, &
         ptcle_dtau_layers, &
         ptcle_ssa_layers, &
         gas_dtauabs_layers, &
         gas_dtausca_layers, &
         alt_atm, &
         nalt_atm, &
         sinsca_corint, &
         p_atm, &
         mode, &
         kdis_nmaxai, &
         kdis_nai, &
         kdis_nai_c, &
         kdis_nsp, &
         kdis_nsp_c, &
         warning, &
         nesp_layers, &
         naimax_layers, &
         nai_layers, &
         ai_layers, &
         ptcle_dtauabs_layers, &
         tot_dtausca_layers, &
         nbetal_layers, &
         betal_layers, &
         nmax_aitaui_mono_layers, &
         n_aitaui_mono_layers, &
         n_aitaui_total, &
         ai_mono_layers, &
         tot_dtau_mono_layers, &
         ssa_mono_layers, &
         gas_dtau_mono_layers,&
         gas_ssa_mono_layers,&
         lambda_ikdis, &
         ptcle_dtau_layers_unscaled, &      
         ptcle_ssa_layers_unscaled, &     
         ptcle_phasemat_TMS_layers, &
         no_rayleigh, &
         print_aitaui

    USE MUTILITY, only : gaussian, linintpol, plint
    !USE nr, ONLY : gauleg
    USE nm_nm, ONLY : nm_gaussNodes
    
    IMPLICIT NONE

    REAL(KIND=dp) :: gl_ptcle_ssa(nptcle,nlambda) 
    REAL(KIND=dp) :: gl_ptcle_tau(nptcle,nlambda) 
    REAL(KIND=dp) :: gl_gas_cumul_tausca(nalt_atm,nlambda)                ! the cumulative gas scattering optical thickness from the top of the atmosphere  
    REAL(KIND=dp) :: gl_ptcle_cumul_tau(nptcle, nalt_atm,nlambda)         ! cumulative particle optical thickness from the top of the atmosphere 
    !                                                                       renormalized ssa and optical depth following the truncation
    ! Same but unscaled (NB : scaling is done in accordance with phase matrix truncation)
    REAL(KIND=dp) :: gl_ptcle_ssa_unscaled(nptcle,nlambda) 
    REAL(KIND=dp) :: gl_ptcle_tau_unscaled(nptcle,nlambda) 
    REAL(KIND=dp) :: gl_ptcle_cumul_tau_unscaled(nptcle, nalt_atm,nlambda)  

    REAL(KIND=dp) :: max_alt 
    REAL(KIND=dp) :: gl_distrib(nalt_distrib)
    REAL(KIND=dp) :: gl_alt_distrib(nalt_distrib)
    REAL(KIND=dp) :: gl_dalt
    REAL(KIND=dp) :: gl_res_int

    ! for Betal coeff. of the Rayleigh
    REAL(kind=dp) :: C
    REAL(kind=dp) :: D

    ! local layer properties
    REAL(KIND=dp) :: gl_fgas(nlayers,nlambda)                         !  scattering fraction of gas
    REAL(KIND=dp) :: gl_fptcle(nptcle, nlayers,nlambda)               !  scattering fraction of each ptcle  
    REAL(KIND=dp) :: gl_gas_dtausca_layers(nlayers,nlambda)           !  Gas scattering optical thickness for each layer
    REAL(KIND=dp), ALLOCATABLE :: gl_gas_dtauabs_layers(:,:,:,:)      !  gl_gas_dtauabs_layers(nlayers, nesp_layers, naimax_layers, nlambda) : Gas absorption optical thickness 
    INTEGER, ALLOCATABLE :: gl_nai_layers(:,:)                        !  gl_nai_layers(nesp_layers, nlambda)
    REAL(KIND=dp), ALLOCATABLE :: gl_ai_layers(:,:,:)                 !  gl_ai_layers(nesp_layers, nesp_layers, naimax_layers,nlambda) : kdis weight 
    REAL(KIND=dp) :: gl_ptcle_dtausca_layers(nptcle,nlayers,nlambda)  !  scattering optical thickness for each layer and each particle population
    REAL(KIND=dp) :: gl_ptcle_dtauabs_layers(nptcle,nlayers,nlambda)  !  absorption optical thickness for each layer and each particle population
    ! Same but unscaled (NB : scaling is done in accordance with phase matrix truncation)
    REAL(KIND=dp) :: gl_ptcle_dtausca_layers_unscaled(nptcle,nlayers,nlambda)  
    REAL(KIND=dp) :: gl_ptcle_dtauabs_layers_unscaled(nptcle,nlayers,nlambda)  


    REAL(KIND=dp) :: gl_tot_dtausca_layers(nlayers,nlambda)           !  Total scattering optical thickness for each layer 
    REAL(kind=dp) :: gl_fsca_ptcle
    REAL(kind=dp) :: gl_fsca_ptcle_unscaled

    REAL(kind=dp) :: gauss_w(nang_max)
    REAL(kind=dp) :: P_tmp

    ! loop index
    INTEGER :: i, k, u, j, l
    integer :: ilamb, isp, ikdislamb

    integer :: iesp
    integer :: iai

    integer :: i_test

    !----------- 



    IF (verbose) THEN
       WRITE (*,*) ''
       WRITE (*, *)'  (sub. get_layers) allocate arrays...'
    ENDIF

    ! dimension and allocate arrays
    if (mode.eq.'kdis') then
       nesp_layers   = kdis_nsp + kdis_nsp_c
       naimax_layers = kdis_nmaxai
       nmax_aitaui_mono_layers = 1
       n_aitaui_total = 0  
       DO ilamb = 1, nlambda
          n_aitaui_mono_layers = 1
          ikdislamb = lambda_ikdis(ilamb)
          DO isp = 1, kdis_nsp
             n_aitaui_mono_layers = n_aitaui_mono_layers * kdis_nai(isp,ikdislamb)
          ENDDO
          DO isp = 1, kdis_nsp_c
             n_aitaui_mono_layers = n_aitaui_mono_layers * kdis_nai_c(isp,ikdislamb)
          ENDDO
          n_aitaui_total = n_aitaui_total + n_aitaui_mono_layers
          if (n_aitaui_mono_layers.gt.nmax_aitaui_mono_layers) then
             nmax_aitaui_mono_layers = n_aitaui_mono_layers
          endif
       ENDDO
       n_aitaui_mono_layers = 1
    else
       ! the following are default values used in case no k-dist. is performed
       nesp_layers   = 1
       naimax_layers = 1
       nmax_aitaui_mono_layers = 1
       n_aitaui_total = nlambda
    endif

    ! start allocate COMMON arrays
    ALLOCATE(nai_layers(nesp_layers, nlambda), &
         ai_layers(nesp_layers, naimax_layers, nlambda))
    nai_layers(:, :)              = undef_i
    ai_layers(:, :, :)            = undef_dp

    ALLOCATE(gas_dtauabs_layers(nlayers, nesp_layers, naimax_layers, nlambda))
    gas_dtauabs_layers(:,:,:,:)    = 0.0_dp

    IF ( (rt_model .eq. 'mcrad1d') .or. &
         (rt_model .eq. 'sinsca') .or. &
         (sinsca_corint)) THEN
       ALLOCATE(ptcle_dtau_layers(nlayers,nlambda), &
            gas_dtausca_layers(nlayers,nlambda), &
            ptcle_ssa_layers(nlayers,nlambda), & 
            ptcle_nang_phasemat_layers(nlambda), & 
            ptcle_phasemat_layers(nlayers, 6, nang_max,nlambda), &
            ptcle_u_phasemat_layers(nang_max, nlambda),&
            gas_dtau_mono_layers(nlayers,nmax_aitaui_mono_layers) ,&
            gas_ssa_mono_layers(nlayers,nmax_aitaui_mono_layers))
       ptcle_dtau_layers(:,:)         = 0.0_dp
       gas_dtausca_layers(:,:)        = 0.0_dp
       ptcle_ssa_layers(:,:)          = undef_dp
       ptcle_nang_phasemat_layers(:)  = 0
       ptcle_phasemat_layers(:,:,:,:) = 0.0_dp
       ptcle_u_phasemat_layers(:,:)   = 0.0_dp
       gas_dtau_mono_layers(:,:)      = undef_dp
       gas_ssa_mono_layers(:,:)       = undef_dp

    ELSEIF(print_aitaui) then

       ALLOCATE(gas_dtausca_layers(nlayers,nlambda), &
            gas_dtau_mono_layers(nlayers,nmax_aitaui_mono_layers) ,&
            gas_ssa_mono_layers(nlayers,nmax_aitaui_mono_layers))
       gas_dtausca_layers(:,:)        = 0.0_dp
       gas_dtau_mono_layers(:,:)      = undef_dp
       gas_ssa_mono_layers(:,:)       = undef_dp

    ENDIF

    IF ( (rt_model.eq.'disort') .or. &
         (rt_model.eq.'addoub') .or. &
         (rt_model.eq.'doad') ) THEN
       ! print*, 8.0* DBLE(nesp_layers)*  DBLE(nlambda) / 1e9
       ! print*, 8.0* DBLE(nesp_layers)* DBLE(naimax_layers)*  DBLE(nlambda) / 1e9
       ! print*, 8.0* DBLE(nlayers)*  DBLE(nesp_layers)*  DBLE(naimax_layers)*  DBLE(nlambda) / 1e9
       ! print*, 8.0* DBLE(nlayers)*  DBLE(nlambda) / 1e9
       ! print*, 8.0* DBLE(nlambda) / 1e9
       ! print*, 8.0* DBLE(nlayers)*  DBLE(nlambda) *  DBLE(nmaxbetal+1) * 6 / 1e9
       ! print*, 8.0* DBLE(nmax_aitaui_mono_layers) / 1e9
       ! print*, 8.0* DBLE(nlayers) *  DBLE(nmax_aitaui_mono_layers) / 1e9
       ! print*, 8.0* DBLE(nlayers) *  DBLE(nmax_aitaui_mono_layers) / 1e9
       ! print*, 8.0* DBLE(nlayers) *  DBLE(nmaxbetal+1) * 6 / 1e9
       ALLOCATE(ptcle_dtauabs_layers(nlayers, nlambda), &
            tot_dtausca_layers(nlayers, nlambda), &
            nbetal_layers(nlambda), &
            betal_layers(nlayers,6,0:nmaxbetal,nlambda), &
            tot_dtau_mono_layers(nlayers, nmax_aitaui_mono_layers), &
            ssa_mono_layers(nlayers, nmax_aitaui_mono_layers))
       ptcle_dtauabs_layers(:, :)    = 0.0_dp
       tot_dtausca_layers(:,:)       = 0.0_dp
       nbetal_layers(:)              = undef_i
       betal_layers(:,:,:,:)         = undef_dp
       tot_dtau_mono_layers(:,:)     = undef_dp
       ssa_mono_layers(:,:)          = undef_dp
    ENDIF

    ALLOCATE(ai_mono_layers(nmax_aitaui_mono_layers))
    ai_mono_layers(:)             = undef_dp

    if (sinsca_corint) then
       ALLOCATE(ptcle_ssa_layers_unscaled(nlayers,nlambda), &
            ptcle_dtau_layers_unscaled(nlayers,nlambda), & 
            ptcle_phasemat_tms_layers(nlayers, 6, nang_max,nlambda))
       ptcle_ssa_layers_unscaled(:,:)      = 0.0_dp
       ptcle_dtau_layers_unscaled(:,:)      = 0.0_dp
       ptcle_phasemat_tms_layers(:,:,:,:)  = 0.0_dp
    endif
    ! end dimension and allocate COMMON layers definition arrays

    ALLOCATE(gl_gas_dtauabs_layers(nlayers,nesp_layers, naimax_layers, nlambda), &
         gl_ai_layers(nesp_layers, naimax_layers,nlambda),&
         gl_nai_layers(nesp_layers, nlambda))
    gl_gas_dtauabs_layers(:,:,:,:) = 0.0_dp
    gl_nai_layers(:,:)             = 0
    gl_ai_layers(:,:,:)            = 0.0_dp
    IF (nptcle .gt. 0) THEN
       gl_ptcle_ssa(:,:)              = 0.0_dp 
       gl_ptcle_tau(:,:)              = 0.0_dp 
       gl_ptcle_cumul_tau(:,:,:)      = 0.0_dp 
       gl_ptcle_dtausca_layers(:,:,:) = 0.0_dp 
       gl_ptcle_dtauabs_layers(:,:,:) = 0.0_dp 
       gl_ptcle_ssa_unscaled(:,:)              = 0.0_dp 
       gl_ptcle_tau_unscaled(:,:)              = 0.0_dp 
       gl_ptcle_cumul_tau_unscaled(:,:,:)      = 0.0_dp 
       gl_ptcle_dtausca_layers_unscaled(:,:,:) = 0.0_dp 
       gl_ptcle_dtauabs_layers_unscaled(:,:,:) = 0.0_dp 
       gl_fptcle(:,:,:)               = 0.0_dp 
    ENDIF
    gl_fgas(:,:)                 = 0.0_dp
    gl_gas_cumul_tausca(:,:)     = 0.0_dp 
    gl_gas_dtausca_layers(:,:)   = 0.0_dp 
    gl_tot_dtausca_layers(:,:)   = 0.0_dp 

    IF (verbose) THEN
       WRITE (*,*) ''
       WRITE (*, *)'  (sub. get_layers) set layer properties'
    ENDIF

    ! gl_distrib and gl_alt_distrib are intermediary 
    ! arrays for the distribution of particles 
    ! to be integrated to get the cumulative optical depth
    IF (nptcle .gt. 0) THEN
       max_alt = MAXVAL(alt_atm)
       gl_alt_distrib(:)       = undef_dp
       ! we create a well sampled logarithmic altitude grid 
       ! for the integration of particule distribution 
       gl_dalt = ( LOG10(max_alt) - LOG10(min_non0_alt) ) / (nalt_distrib - 3)
       DO i = 2, nalt_distrib-1
          gl_alt_distrib(i) = 10.0_dp**( LOG10(min_non0_alt) + (i-2) * gl_dalt)
       ENDDO
       gl_alt_distrib(nalt_distrib) = maxval(alt_atm)*1.01_dp
       gl_alt_distrib(1) = 0.0_dp
    ENDIF

    CALL SET_GASABS(gl_gas_dtauabs_layers, gl_ai_layers, gl_nai_layers)

    !----- START LOOPING OVER WAVELENGHTS
    ! in this loop, we compute the following variables :
    ! gl_fgas(nlayers,nlambda)                         !  scattering fraction of gas
    ! gl_fptcle(nptcle, nlayers,nlambda)               !  scattering fraction of each ptcle  
    ! gl_gas_dtausca_layers(nlayers,nlambda)           !  Gas scattering optical thickness for each layer
    ! gl_ptcle_dtausca_layers(nptcle,nlayers,nlambda)  !  scattering optical thickness for each layer and each particle population
    ! gl_ptcle_dtauabs_layers(nptcle,nlayers,nlambda)  !  absorption optical thickness for each layer and each particle population
    ! gl_tot_dtausca_layers(nlayers,nlambda)           !  Total scattering optical thickness for each layer 
    DO j = 1, nlambda

       ! set the particle properties (SSA & TAU) at a given wavelength
       ! regarding the reference wavelength and accounting for possible
       ! truncation of the phase matrix 
       IF (nptcle .gt. 0) THEN
          DO i = 1, nptcle

             if (trunc_method.ne.'none') then

                ! compute renormalized ssa and optical depth following the truncation
                gl_ptcle_ssa(i,j) = (1.0_dp - ptcle_trunccoeff(i,j)) * ptcle_opt(i,2,j) / &
                     (1.0_dp - (ptcle_opt(i,2,j) * ptcle_trunccoeff(i,j)))
                ! compute optical depth at the working wavelengths (based on the one at reference wavelengths)
                gl_ptcle_tau(i,j) = ptcle_tau_reflamb(i) * ptcle_opt(i,1,j) / ptcle_Cext_reflamb(i)
                ! computed renormalized optical depth following truncation
                gl_ptcle_tau(i,j) = gl_ptcle_tau(i,j) * (1.0_dp - (ptcle_trunccoeff(i,j) * ptcle_opt(i,2,j)))                

                if (sinsca_corint) then
                   gl_ptcle_ssa_unscaled(i,j) = ptcle_opt(i,2,j)
                   ! compute optical depth at the working wavelengths (based on the one at reference wavelengths)
                   gl_ptcle_tau_unscaled(i,j) = ptcle_tau_reflamb(i) * ptcle_opt(i,1,j) / ptcle_Cext_reflamb(i)
                endif

             else
                gl_ptcle_ssa(i,j) =  ptcle_opt(i,2,j)
                ! compute optical depth at the working wavelengths (based on the one at reference wavelengths)
                gl_ptcle_tau(i,j) = ptcle_tau_reflamb(i) * ptcle_opt(i,1,j) / ptcle_Cext_reflamb(i)
             endif

             ! print*, i, j, ptcle_trunccoeff(i,j), ptcle_opt(i,2,j),  gl_ptcle_ssa(i,j), &
             !      ptcle_tau_reflamb(i) * ptcle_opt(i,1,j) / ptcle_Cext_reflamb(i),  gl_ptcle_tau(i,j)


          ENDDO

       ENDIF

       ! compute the cumulative optical thicknesses for the gas
       if (no_rayleigh.eqv..FALSE.) then
          DO i = 1, nalt_atm
             gl_gas_cumul_tausca(i,j) = TAU_RAYLEIGH(p_atm(i), wlambda(j)) 
          ENDDO
       endif

       ! compute the cumulative optical thickness for the particle populations
       ! to be integrated and used to distribute among the different layers
       ! (not to be done if the distribution is 'layer')
       IF (nptcle .gt. 0) THEN
          DO i = 1, nptcle
             SELECT CASE(ptcle_vdist_type(i))

             CASE('homogeneous') 

                DO k = 1, nalt_atm                   
                   IF (alt_atm(k) .ge. ptcle_h_vdist(i)) THEN
                      gl_ptcle_cumul_tau(i,k,j) = 0.0_dp
                   ELSE IF( (alt_atm(k) .lt. ptcle_h_vdist(i) ) .and. &
                        (alt_atm(k) .gt. (ptcle_h_vdist(i) - ptcle_sd_gauss_vdist(i))) ) THEN
                      gl_ptcle_cumul_tau(i,k,j) = 1.0_dp - &
                           ((alt_atm(k) - (ptcle_h_vdist(i) - ptcle_sd_gauss_vdist(i))) / ptcle_sd_gauss_vdist(i))
                   ELSE
                      gl_ptcle_cumul_tau(i,k,j) = 1.0_dp
                   ENDIF

                ENDDO

             CASE('layer') 
                ! ----- just place the particle in a single layer
                DO k = 1, nalt_atm
                   IF (alt_atm(k) .gt. ptcle_h_vdist(i)) THEN
                      gl_ptcle_cumul_tau(i,k,j) = 0.0_dp
                   ELSE
                      gl_ptcle_cumul_tau(i,k,j) = 1.0_dp
                   ENDIF
                ENDDO

             CASE('gauss') 
                ! ----- get the particule altitude distribution
                gl_distrib(:) = 0.0_dp
                DO k = 1, nalt_distrib
                   gl_distrib(k) = gaussian(gl_alt_distrib(k), ptcle_h_vdist(i), ptcle_sd_gauss_vdist(i))
                   !if (gl_distrib(k).ne.0.0) print*, gl_alt_distrib(k),gl_distrib(k)
                ENDDO
                ! get the cumulative optical depth
                gl_ptcle_cumul_tau(i,1,j) = 0.0_dp
                DO k = 2, nalt_atm
                   gl_res_int = 0.0_dp             
                   call plint(gl_distrib, gl_alt_distrib, nalt_distrib, alt_atm(k), maxval(alt_atm(:)), gl_res_int )     
                   gl_ptcle_cumul_tau(i,k,j) = gl_res_int
                ENDDO

             CASE('sh')
                ! ----- get the particule altitude distribution
                gl_distrib(:) = 0.0_dp
                DO k = 1, nalt_distrib
                   gl_distrib(k) = EXP( - gl_alt_distrib(k) / ptcle_h_vdist(i) )
                   !if (gl_distrib(k).ne.0.0) print*, gl_alt_distrib(k),gl_distrib(k)
                ENDDO
                ! get the cumulative optical depth
                gl_ptcle_cumul_tau(i,1,j) = 0.0_dp
                DO k = 2, nalt_atm
                   gl_res_int = 0.0_dp             
                   call plint(gl_distrib, gl_alt_distrib, nalt_distrib, alt_atm(k), maxval(alt_atm(:)), gl_res_int )     
                   gl_ptcle_cumul_tau(i,k,j) = gl_res_int
                ENDDO

             CASE('user')
                ! DO k = 1, ptcle_user_vdist_nalt(i)
                !    print*, ptcle_user_vdist_alt(i,k),  ptcle_user_vdist(i,k)
                ! ENDDO
                ! get the cumulative optical depth
                gl_ptcle_cumul_tau(i,1,j) = 0.0_dp
                DO k = 2, nalt_atm
                   gl_res_int = 0.0_dp             
                   call plint(ptcle_user_vdist(i,1:ptcle_user_vdist_nalt(i)), &
                        ptcle_user_vdist_alt(i,1:ptcle_user_vdist_nalt(i)), &
                        ptcle_user_vdist_nalt(i), &
                        alt_atm(k), maxval(alt_atm(:)), gl_res_int ) 
                   gl_ptcle_cumul_tau(i,k,j) = gl_res_int
                ENDDO

             END SELECT
             ! normalization
             if (sinsca_corint) then
                gl_ptcle_cumul_tau_unscaled(i,:,j) = gl_ptcle_cumul_tau(i,:,j) * gl_ptcle_tau_unscaled(i,j) / &
                     maxval(gl_ptcle_cumul_tau(i,:,j))
             endif
             gl_ptcle_cumul_tau(i,:,j) = gl_ptcle_cumul_tau(i,:,j) * gl_ptcle_tau(i,j) / &
                  maxval(gl_ptcle_cumul_tau(i,:,j))
          ENDDO
       ENDIF

       ! - get layer properties
       DO i = 1, nlayers
          IF (nptcle .gt. 0) THEN
             DO k = 1, nptcle
                if (sinsca_corint) then
                   gl_ptcle_dtausca_layers_unscaled(k,i,j) = &
                        (gl_ptcle_cumul_tau_unscaled(k,i+1,j) - gl_ptcle_cumul_tau_unscaled(k,i,j)) * &
                        gl_ptcle_ssa_unscaled(k,j)
                   gl_ptcle_dtauabs_layers_unscaled(k,i,j) = &
                        (gl_ptcle_cumul_tau_unscaled(k,i+1,j) - gl_ptcle_cumul_tau_unscaled(k,i,j)) * &
                        (1.0_dp-gl_ptcle_ssa_unscaled(k,j))
                endif
                gl_ptcle_dtausca_layers(k,i,j) = (gl_ptcle_cumul_tau(k,i+1,j) - gl_ptcle_cumul_tau(k,i,j)) * &
                     gl_ptcle_ssa(k,j)
                gl_ptcle_dtauabs_layers(k,i,j) = (gl_ptcle_cumul_tau(k,i+1,j) - gl_ptcle_cumul_tau(k,i,j)) * &
                     (1.0_dp-gl_ptcle_ssa(k,j))
                ! add ptcle dtau to layers dtau
                gl_tot_dtausca_layers(i,j)     = gl_tot_dtausca_layers(i,j) + gl_ptcle_dtausca_layers(k,i,j)
             ENDDO
          ENDIF
          ! gas optical thicknesses
          if (no_rayleigh) then
             gl_gas_dtausca_layers(i,j) = 0.0_dp
          else
             gl_gas_dtausca_layers(i,j) = gl_gas_cumul_tausca(i+1,j) - gl_gas_cumul_tausca(i,j)
          endif
          ! add gas dtau to layer dtau
          gl_tot_dtausca_layers(i,j)     = gl_tot_dtausca_layers(i,j)     + gl_gas_dtausca_layers(i,j)
          ! get the scattering fraction of gas and each particles
          if (gl_tot_dtausca_layers(i,j).eq.0.0_dp) then
             ! to avoid problems in building Betal
             gl_fgas(i,j)      = 1.0_dp 
             gl_fptcle(:,i,j)  = 0.0_dp
          else
             gl_fgas(i,j) = gl_gas_dtausca_layers(i,j) / gl_tot_dtausca_layers(i,j)  
             DO k = 1, nptcle
                gl_fptcle(k,i,j) = gl_ptcle_dtausca_layers(k,i,j) / gl_tot_dtausca_layers(i,j) 
             END DO
          endif
       ENDDO! loop on layers
    END DO  ! lambda

    if (no_rayleigh.eqv..FALSE.) then
       IF (verbose) THEN
          WRITE (*,*) ''
          WRITE (*, FMT='(A)') '   (sub. get_layers)     lambda       Tot. Ray. Opacity'
          DO j = 1, nlambda
             WRITE (*, FMT='(21x,F12.7,5x,F12.7)') wlambda(j),  gl_gas_cumul_tausca(nalt_atm,j)
          ENDDO
       ENDIF
    endif

    !============================================
    !          implement common variables
    !============================================
    
    ! Sanity check on gas abs  
    DO j = 1, nlambda
       do iesp = 1, nesp_layers
          do iai = 1, gl_nai_layers(iesp,j)         
             do i = 1, nlayers      
                if (gl_gas_dtauabs_layers(i,iesp,iai,j) < 0.0_dp) then
                   if (warning) then
                      WRITE(*,*) ''
                      WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      WRITE(*,*) ' (sub. get_layers)    WARNING                       '
                      WRITE(*,*) '    Negative gas abs. opacity found                      '
                      WRITE(*,*) '       Value =  ', gl_gas_dtauabs_layers(i,iesp,iai,j)
                      WRITE(*,*) '     Reset it to 0.0                     '
                      WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                      WRITE(*,*) ''
                   endif
                   gl_gas_dtauabs_layers(i,iesp,iai,j) = 0.0_dp                
                endif
             end do
          end do
       end do
    end DO

    IF ( (rt_model .eq. 'mcrad1d') .or. &
         (rt_model .eq. 'sinsca') .or. &
         (sinsca_corint)) THEN

       DO j = 1, nlambda

          ! effective k-dis weight for each wavelength
          nai_layers(:,j)  = gl_nai_layers(:,j)
          ai_layers(:,:,j) = gl_ai_layers(:,:,j)
          do i=1, nlayers
             gas_dtauabs_layers(i,:,:,j) = gl_gas_dtauabs_layers(i,:,:,j)
             gas_dtausca_layers(i,j)     = gl_gas_dtausca_layers(i,j)
             if (nptcle.gt.0) then
                ptcle_dtau_layers(i,j) = SUM(gl_ptcle_dtausca_layers(:,i,j) + gl_ptcle_dtauabs_layers(:,i,j))
                ptcle_ssa_layers(i,j)  = SUM(gl_ptcle_dtausca_layers(:,i,j)) / ptcle_dtau_layers(i,j)
                if (ptcle_dtau_layers(i,j).eq.0.0_dp)  ptcle_ssa_layers(i,j) = 0.0_dp
                if (sinsca_corint) then
                   ptcle_dtau_layers_unscaled(i,j) = &
                        SUM(gl_ptcle_dtausca_layers_unscaled(:,i,j) + gl_ptcle_dtauabs_layers_unscaled(:,i,j))
                   ptcle_ssa_layers_unscaled(i,j)  = SUM(gl_ptcle_dtausca_layers_unscaled(:,i,j)) / ptcle_dtau_layers_unscaled(i,j)
                   if (ptcle_dtau_layers_unscaled(i,j).eq.0.0_dp)  ptcle_ssa_layers_unscaled(i,j) = 0.0_dp
                end if
             endif
          end do

          ! For Monte-Carlo we get the phase matrix for each layer.
          ! We only compute the layer phase matrix if their are ptcle
          ! because -ptcle_phasemat_layers- must account for ptcle ONLY
          ! NO Rayleigh contribution
          IF (nptcle .gt. 0) then

             !========================================
             ! we do the interpolation of
             ! particles phase matrix

             ! get the number of mu 
             mc_nmu = maxval(ptcle_nang_phasemat(:, j)) ! cause we did not read mcrad1d input
             ptcle_nang_phasemat_layers(j) = mc_nmu

             ! We create a common gauss grid :
             !CALL gauleg(1.0_dp, -1.0_dp, ptcle_u_phasemat_layers(2:mc_nmu-1,j), gauss_w(2:mc_nmu-1))
             CALL nm_gaussNodes(1.0_dp, -1.0_dp, ptcle_u_phasemat_layers(2:mc_nmu-1,j), gauss_w(2:mc_nmu-1))

             ! set the extremum values by hand cause gauss grid does not include
             ptcle_u_phasemat_layers(1,j)      = +1.0_dp
             ptcle_u_phasemat_layers(mc_nmu,j) = -1.0_dp

             IF (trunc_method .eq. 'none') then

                ! We use the native particle phase matrix 
                DO i = 1, nlayers
                   DO k = 1, nptcle
                      gl_fsca_ptcle = gl_ptcle_dtausca_layers(k,i,j) / &
                           (ptcle_dtau_layers(i,j)*ptcle_ssa_layers(i,j))
                      if ( gl_fsca_ptcle .gt. 0) then
                         DO u = 1, mc_nmu
                            DO l = 1, 6 
                               ! we interpolate the phase matrix on the gauss points grid
                               ! INTPOL(fint, xint, ni, xess)
                               ! linear interoplation of fint @ xess
                               ! xint assumed increasing, NO EXTAPOLATION 
                               ! make sure xess belongs to [xint(1),xint(ni)]
                               write(*,*) ptcle_u_phasemat_layers(u, j)
                               !i_test = ptcle_u_phasemat_layers(u, j)
                               !P_tmp = ptcle_phasemat(k, l, ptcle_u_phasemat_layers(u, j), j)
                               !write(*,*) "P_tmp_value = ", P_tmp
                               P_tmp = LININTPOL(ptcle_phasemat(k,l,1:ptcle_nang_phasemat(k,j),j), &
                                    &            ptcle_u_phasemat(k,1:ptcle_nang_phasemat(k,j),j), & 
                                    &            ptcle_nang_phasemat(k,j), ptcle_u_phasemat_layers(u,j))    
                               write(*,*) "P_tmp_interpolation = ", P_tmp
                               ! contribution of the ptcle to the scattering in the layer:
                               ptcle_phasemat_layers(i, l, u, j) = &
                                    ptcle_phasemat_layers(i, l, u, j) + gl_fsca_ptcle *  P_tmp
                            ENDDO
                         ENDDO
                      endif
                   ENDDO
                ENDDO

             ELSE

                ! we use the truncated particle phase matrix 
                DO i = 1, nlayers
                   DO k = 1, nptcle
                      gl_fsca_ptcle = gl_ptcle_dtausca_layers(k,i,j) / &
                           (ptcle_dtau_layers(i,j)*ptcle_ssa_layers(i,j))
                      if (sinsca_corint) then
                         gl_fsca_ptcle_unscaled = gl_ptcle_dtausca_layers_unscaled(k,i,j) / &
                              (ptcle_dtau_layers_unscaled(i,j)*ptcle_ssa_layers_unscaled(i,j))
                      endif
                      if ( gl_fsca_ptcle .gt. 0) then
                         DO u = 1, mc_nmu
                            DO l = 1, 6 
                               ! we interpolate the phase matrix on the gauss points grid
                               ! INTPOL(fint, xint, ni, xess)
                               ! linear interoplation of fint @ xess
                               ! xint assumed increasing, NO EXTAPOLATION 
                               ! make sure xess belongs to [xint(1),xint(ni)]
                               P_tmp = LININTPOL(ptcle_trunc_phasemat(k,l,1:ptcle_nang_phasemat(k,j),j), &
                                    ptcle_u_phasemat(k,1:ptcle_nang_phasemat(k,j),j), & 
                                    ptcle_nang_phasemat(k,j), ptcle_u_phasemat_layers(u,j))  
                               ! contribution of the ptcle to the scattering in the layer:
                               ptcle_phasemat_layers(i, l, u, j) = &
                                    ptcle_phasemat_layers(i, l, u, j) + gl_fsca_ptcle *  P_tmp                            
                               if (sinsca_corint) then                                                                
                                  ! phase function used in TMS correction; actual phase
                                  ! function divided by (1-FLYR*SSALB)
                                  P_tmp = LININTPOL(ptcle_phasemat(k,l,1:ptcle_nang_phasemat(k,j),j), &
                                       &            ptcle_u_phasemat(k,1:ptcle_nang_phasemat(k,j),j), & 
                                       &            ptcle_nang_phasemat(k,j), ptcle_u_phasemat_layers(u,j))
                                  P_tmp = P_tmp * gl_fsca_ptcle_unscaled
                                  ptcle_phasemat_tms_layers(i, l, u, j) = &
                                       ptcle_phasemat_tms_layers(i, l, u, j) + &
                                       P_tmp  / (1.0_dp - (ptcle_trunccoeff(k,j) * ptcle_opt(k,2,j)))
                               endif
                            ENDDO
                         ENDDO
                      endif
                   ENDDO
                ENDDO
             ENDIF

          ENDIF! (nptcle .gt. 0)
       ENDDO! lambda

    ELSEIF(print_aitaui) then

       DO j = 1, nlambda
          ! effective k-dis weight for each wavelength
          nai_layers(:,j)  = gl_nai_layers(:,j)
          ai_layers(:,:,j) = gl_ai_layers(:,:,j)
          do i=1, nlayers
             gas_dtauabs_layers(i,:,:,j) = gl_gas_dtauabs_layers(i,:,:,j)
             gas_dtausca_layers(i,j)     = gl_gas_dtausca_layers(i,j)
          end do
       END DO

    ENDIF

    IF ( (rt_model.eq.'disort') .or. &
         (rt_model.eq.'addoub') .or. &
         (rt_model.eq.'doad') ) THEN

       DO j=1, nlambda

          ! effective k-dis weight for each wavelength
          nai_layers(:,j)  = gl_nai_layers(:,j)
          ai_layers(:,:,j) = gl_ai_layers(:,:,j)
          do i=1, nlayers
             gas_dtauabs_layers(i,:,:,j) = gl_gas_dtauabs_layers(i,:,:,j)
             if (nptcle.gt.0) then
                ptcle_dtauabs_layers(i,j) = SUM(gl_ptcle_dtauabs_layers(:,i,j))
             end if
             tot_dtausca_layers(i,j)     = gl_tot_dtausca_layers(i,j)
          end do

          ! get Legendre coefficients for alpha1, alpha2 
          ! alpha3, alpha4, beta1, beta2 elements of the 
          ! single scattering phase matrix for each layers 
          ! see Section 4 of de Haan et al, ApJ, 1987, 183, 371
          ! get the highest order of Legendre expension coefficient
          ! among all atmospheric species (gas, particle)
          ! NOTE : this highest value will be the one used  in 
          !        setting the Legendre coeff. of each layers.
          IF (nptcle .gt. 0 ) THEN
             nbetal_layers(j) = maxval(ptcle_nbetal(:,j))
          ELSE
             IF (rt_model .eq. 'disort') THEN
                ! cause in case we only have Rayleigh, nbetal = 2
                ! but in DISORT, nbetal must be >= nstr
                nbetal_layers(j) = od_nstr
             ELSE
                nbetal_layers(j) = 2
             ENDIF

          ENDIF
          betal_layers(:,:,:, j) = 0.0_dp
          DO i = 1, nlayers
             IF (nptcle .gt. 0) THEN             
                DO k = 1, nptcle
                   DO u = 0, ptcle_nbetal(k,j)
                      betal_layers(i,:,u, j) = betal_layers(i,:,u, j) + gl_fptcle(k,i,j) * ptcle_betal(k,:,u,j)
                   ENDDO
                ENDDO
             ENDIF

             ! Rayleigh
             ! Formulae are given in de Rooij (1985) page 40    
             C = 2.0_dp*(1.0_dp - depol(j)) / ( 2.0_dp + depol(j))
             D = 2.0_dp*(1.0_dp - 2.0_dp * depol(j)) / ( 2.0_dp + depol(j))
             ! alpha1 n=0
             betal_layers(i, 1, 0, j) = betal_layers(i, 1, 0, j) + gl_fgas(i,j)
             ! alpha4 n=1
             betal_layers(i, 4, 1, j) = betal_layers(i, 4, 1, j) + gl_fgas(i,j) * 1.5_dp * D
             ! alpha1 n=2
             betal_layers(i, 1, 2, j) = betal_layers(i, 1, 2, j) + gl_fgas(i,j) * 0.5_dp * C
             ! alpha2 n=2
             betal_layers(i, 2, 2, j) = betal_layers(i, 2, 2, j) + gl_fgas(i,j) * 3.0_dp * C
             ! beta1 n=2  
             betal_layers(i, 5, 2, j) = betal_layers(i, 5, 2, j) + gl_fgas(i,j) * SQRT(1.5_dp) * C
             
             IF (betal_layers(i, 1, 0, j) .gt.  (1.0_dp+EPS_BETAL) .or. betal_layers(i, 1, 0, j) .lt.(1.0_dp-EPS_BETAL)) THEN
                WRITE(*,*) ' '
                WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                WRITE(*,*) ' (in sub. get_layers)     ERROR'
                WRITE(*,*) '    For each layer, Betal '
                WRITE(*,*) '    should be normalized so that Beta_0^11 = 1'
                WRITE(*,*) '    layer      ', i
                WRITE (*, FMT='(A17,2x,F12.5,1x,A8)') &
                     &     '    wavelength  =', wlambda(j),' microns'  
                WRITE (*, FMT='(A14,2x,F12.6)') &
                     &     '    Beta_0^11=',  betal_layers(i, 1, 0, j)
                WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                STOP
             ENDIF
             ! we always renormalize
             betal_layers(i, :, :, j) = betal_layers(i, :, :, j) / betal_layers(i, 1, 0, j)
          ENDDO ! en loop over layers
       END DO ! lambda

    END IF

    DEALLOCATE(gl_gas_dtauabs_layers, &
         gl_ai_layers,&
         gl_nai_layers)

    IF (verbose) THEN
       WRITE(*,*) ' '
       WRITE(*,*) '  (sub. get_layers) exit OK'
       WRITE(*,*) ' '
    ENDIF

  END SUBROUTINE GET_LAYERS

  !--------------------------------------------------------------------------

  FUNCTION TAU_RAYLEIGH(P, lamb) 

    USE MCONSTANTS, only : dp
    IMPLICIT NONE

    ! This routine gives the Rayleigh scattering optical depth given an atmospheric pressure 
    ! and a wavelength (see Hansen & Travis, 1974, Space Science Reviews for the use dnumbers and formula)

    REAL(KIND=dp), INTENT(IN) :: P           ! mb
    REAL(KIND=dp), INTENT(IN) :: lamb        ! microns
    REAL(KIND=dp) :: tau_rayleigh

    REAL(KIND=dp), PARAMETER :: P0 = 1013.25_dp ! mb   
    REAL(KIND=dp) :: tau0                       ! reference optical depth for pressure for P0 at "lambda"

    !-------

    tau0 = (0.008569_dp * lamb**(-4.0_dp) ) * & 
         (1.0_dp + (0.0113_dp*lamb**(-2.0_dp)) + (0.00013_dp*lamb**(-4.0_dp)))
    tau_rayleigh = tau0 * P / p0

    if (tau_rayleigh .lt. 0.0_dp) tau_rayleigh = 0.0_dp


  END FUNCTION TAU_RAYLEIGH

  !--------------------------------------------------------------------------

  SUBROUTINE SET_GASABS(gl_gas_dtauabs_layers, gl_ai_layers, gl_nai_layers)

    ! this routine will access COMMON variables but is supposed 
    ! to modify only gl_gas_dtauabs_layers, gl_ai_layers, gl_nai_layers
    ! we compute the (ai, taui) for each layer, each specy and each wavelength
    ! --> gl_nai_layers(nesp_layers, nlambda)
    !     gl_ai_layers(nesp_layrs,kdis_nmaxai, nlambda)
    !     gl_gas_dtauabs_layers(nlayers,nesp_layers,kdis_nmaxai,nlambda)

    USE MUTILITY, only : BILININTPOL, LININTPOL, TRILININTPOL, SPLINE, SPLINT, XINTEG2
    USE MCONSTANTS, only : dp, &
         undef_dp, &
         undef_i 
    USE MCOMMON, only : mode, &
         nlambda, &
         nlayers, &
         naimax_layers, &
         nalt_atm, &
         alt_atm, &      ! alt_atm(nalt_atm)   : altitude grid (km)
         t_atm, &        ! t_atm(nalt_atm)     : temperature (K)
         p_atm, &        ! p_atm(nalt_atm)     : pressure (mb)
         u_air_atm, &
         ngas_u_atm, &   
         gas_u_atm, &    
         u_atm,&          ! (ngas_u_atm, nalt_atm) concentration profil (cm-3)
         kdis_nmaxai, &
         kdis_nt, &      ! number of temperatures defined in the kdis model
         kdis_np, &      ! number of pressures defined in the kdis model
         kdis_nc, &      ! number of concentrations defined in the kdis model
         kdis_nwvl, &    ! number of wavelengths defined in the kdis model
         kdis_t, &       ! temperatures defined in the kdis model (K)
         kdis_p, &       ! pressures defined in the kdis model (mb)
         kdis_c, &       ! concentrations defined in the kdis model (cm-3)
         kdis_wvlband, & ! kdis wavelengths, kdis_wvlband(3,kdis_nwvl) 
         kdis_nsp, &     ! number of gas species accounted for in the kdis method 
         kdis_species, & ! gas species described in the kdis model
         kdis_nai, &     ! kdis_nai(kdis_nsp,kdis_nwvl) 
         kdis_ai, &      ! kdis_ai(kdis_nsp,kdis_nwvl,kdis_nmaxia) 
         kdis_ki, &      ! kdis_ki(kdis_nsp,kdis_nwvl,kdis_nmaxia,kdis_np,kdis_nt) 
         kdis_xsect,&      ! kdis_xsect(kdis_nsp,kdis_nwvl) 
         kdis_nsp_c, &     ! number of gas species with concentration dependency accounted for in the kdis method 
         kdis_species_c, & ! gas species with concentration dependency described in the kdis model
         kdis_nai_c, &     ! kdis_nai_c(kdis_nsp_c,kdis_nwvl) 
         kdis_ai_c, &      ! kdis_ai_c(kdis_nsp_c,kdis_nwvl,kdis_nmaxia) 
         kdis_ki_c, &      ! kdis_ki_c(kdis_nsp_c,kdis_nwvl,kdis_nmaxia,kdis_np,kdis_nt,kdis_c) 
         kdis_xsect_c,&    ! kdis_xsect_c(kdis_nsp_c,kdis_nwvl) 
         kdis_fcont, &
         kdis_fcont_c, &
         wlambda, &
         warning, &
         nesp_layers, &
         lambda_ikdis
 
    IMPLICIT NONE

    !    REAL(kind=dp), INTENT(OUT) :: gl_gas_dtauabs_layers(nlayers,nesp_layers,naimax_layers, nlambda)
    !    REAL(kind=dp), INTENT(OUT) :: gl_ai_layers(nesp_layers,naimax_layers,nlambda)
    !    INTEGER, INTENT(OUT) :: gl_nai_layers(nesp_layers,nlambda)

    REAL(kind=dp), INTENT(OUT) :: gl_gas_dtauabs_layers(:,:,:, :)
    REAL(kind=dp), INTENT(OUT) :: gl_ai_layers(:,:,:)
    INTEGER, INTENT(OUT) :: gl_nai_layers(:,:)

    integer :: i
    integer :: isp
    integer :: ispatm
    integer :: ispatm_h2o
    integer :: ilamb
    integer :: ikdislamb
    integer :: ilay
    integer :: ialt
    integer :: iai
    logical :: flag_sp
    logical :: flag_sp_h2o

    REAL(kind=dp) :: ki_tmp(nalt_atm)
    REAL(kind=dp) :: tau_tmp(nalt_atm)
    REAL(kind=dp) :: coldens_tmp(nalt_atm)
    logical :: coldens_computed

    integer       :: isp_tmp

    ! single precision variable required to call continuum
    ! subroutines
    REAL(kind=dp) :: anu
    REAL(kind=dp) :: p
    REAL(kind=dp) :: t
    REAL(kind=dp) :: delta_p
    REAL(kind=dp) :: delta_z
    REAL(kind=dp) :: rho
    REAL(kind=dp) :: rhoh2o
    REAL(kind=dp) :: tau_cont

    ! For T interpolation
    integer :: it
    REAL(kind=dp) :: ki_t(kdis_nt) 
    REAL(kind=dp) :: ki_t2(kdis_nt)
    REAL(kind=dp) :: yp1
    REAL(kind=dp) :: ypn

    ! Temporary values
    REAL(kind=dp) :: current_p
    REAL(kind=dp) :: current_t
    REAL(kind=dp) :: current_c

    ! ------
    ! NOTE : if mode=mono the gas abs will be 0
    !       --> no line by line implemented
    gl_gas_dtauabs_layers(:,:,:,:) = 0.0_dp 
    gl_ai_layers(:,:,:)            = 1.0_dp 
    gl_nai_layers(:,:)             = 1

    if (mode.eq.'kdis') then

       ispatm_h2o = 0 
       flag_sp_h2o = .FALSE.

       isp_tmp = 0

       ! first treat the case of gases with no concentration dependency
       do isp = 1, kdis_nsp  ! isp is the index of the specy in kdis_* arrays

          coldens_tmp(:) = 0.0_dp
          coldens_computed = .false.

          ! look for the specy concentration profil
          flag_sp = .FALSE.
          DO i = 1, ngas_u_atm
             if (kdis_species(isp).eq.gas_u_atm(i)) then
                ! the kdis specy has a defined concentration profil
                ! ispatm is the index of the specy in atm profil arrays
                ispatm  = i
                flag_sp = .TRUE.
                exit
             endif
          ENDDO
          IF (flag_sp.eqv..FALSE.) then
             WRITE(*,*) ''
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ' (sub. set_gasabs)    ERROR '
             WRITE(*,*) '   The gas concentration for k-dis gas specy ', &
                  TRIM(ADJUSTL( kdis_species(isp)))
             WRITE(*,*) '   is not defined.'
             WRITE(*,*) '   (see the atmospheric profil file or '
             WRITE(*,*) '   the uniformly distributed gas concentration file)'
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ''
             STOP
          ENDIF
          isp_tmp = isp_tmp + 1  ! index of the specy in gl_nai_layers, gl_ai_layers and gl_gas_dtauabs_layers

          ! for the current specy :
          ! get the (ai,ki) for each layer's (P,T)
          ! then get the taui (do that for each lambda)
          DO ilamb = 1, nlambda

             ikdislamb = lambda_ikdis(ilamb) ! ilamb corresponding index in kdis_ arrays
             gl_nai_layers(isp_tmp, ilamb) = kdis_nai(isp, ikdislamb)

             if (gl_nai_layers(isp_tmp, ilamb).eq.1) then

                gl_ai_layers(isp_tmp, 1, ilamb) = 1.0_dp

                if (kdis_xsect(isp, ikdislamb).gt.0.0_dp) then

                   if (coldens_computed .eqv. .false.) then
                      do ialt = 2, nalt_atm
                         coldens_tmp(ialt) = XINTEG2(1, ialt, nalt_atm, alt_atm*1.0D5, u_atm(ispatm,:) )
                         coldens_computed = .true.
                      end do
                   end if

                   do ilay = 1, nlayers
                      gl_gas_dtauabs_layers(ilay,isp_tmp,1,ilamb) = abs(coldens_tmp(ilay+1)-coldens_tmp(ilay)) &
                           & * kdis_xsect(isp, ikdislamb)
                   end do

                endif

             else

                DO iai = 1, gl_nai_layers(isp_tmp, ilamb)

                   gl_ai_layers(isp_tmp, iai, ilamb) = kdis_ai(isp, ikdislamb, iai)

                   ! compute ki for each atmosphere level
                   do ialt = 1, nalt_atm

                      if (p_atm(ialt).gt.maxval(kdis_p(:))) then
                         current_p = maxval(kdis_p(:))
                      else if (p_atm(ialt).lt.minval(kdis_p(:))) then
                         current_p = minval(kdis_p(:))
                      else
                         current_p = p_atm(ialt)
                      endif
                      
                      if (t_atm(ialt).gt.maxval(kdis_t(:))) then
                         current_t = maxval(kdis_t(:))
                      else if (t_atm(ialt).lt.minval(kdis_t(:))) then
                         current_t = minval(kdis_t(:))
                      else
                         current_t = t_atm(ialt)
                      endif

                      ! NOTE : A LINERA INTERPOLATION ON P AND SPLINE ON T
                      !        IS MADE FOR CONSISTENCY WITH GAME KDIS COEFFICIENTS
                      ! linear interpolation on P
                      DO it = 1, kdis_nt
                         ki_t(it) = LININTPOL(kdis_ki(isp, ikdislamb, iai, :, it),  kdis_p(:), kdis_np, current_p)
                      ENDDO
                      ! SPLINE interpolation on T
                      yp1 = (ki_t(2)-ki_t(1))               / (kdis_t(2) - kdis_t(1))                                           
                      ypn = (ki_t(kdis_nt)-ki_t(kdis_nt-1)) / (kdis_t(kdis_nt) - kdis_t(kdis_nt-1))

                      CALL SPLINE(kdis_t, ki_t, kdis_nt, yp1, ypn, ki_t2)
                      CALL SPLINT(kdis_t, ki_t, ki_t2, kdis_nt, current_t, ki_tmp(ialt))

                      ! bilinear interpolation on P and T 
                      ! ki_tmp(ialt) = BILININTPOL(kdis_ki(isp, ikdislamb, iai, :,:), kdis_p(:), kdis_t(:), &
                      !     kdis_np, kdis_nt, p_atm(ialt), t_atm(ialt))                                       

                   ENDDO

                   ! Compute integrated opacity
                   tau_tmp(:) = 0.0D0
                   DO ialt = 2, nalt_atm                      
                      !                                             z(cm)         
                      tau_tmp(ialt) = XINTEG2(1, ialt, nalt_atm, alt_atm*1.0D5,  ki_tmp * u_atm(ispatm,:) )
                   END DO

                   ! compute tau layers
                   DO ilay = 1, nlayers
                      gl_gas_dtauabs_layers(ilay, isp_tmp, iai, ilamb) = abs(tau_tmp(ilay+1) - tau_tmp(ilay))
                   END DO

                END DO

             endif
             ! ADD CONTINUUM CONTRIBUTION
             if (kdis_fcont(isp).gt.0.0_dp) then
                DO ilay = 1, nlayers
                   tau_cont = 0.0
                   anu      = 1.0d4/wlambda(ilamb)
                   !if ((anu.gt.50000.0).or.(anu.lt.200.0)) then
                   !   WRITE(*,*) 'The validity of continuum sub. outside GAME wavelength must be checked'
                   !   STOP
                   !endif

                   p       = (p_atm(ilay+1) + p_atm(ilay)) / 2.0_dp
                   t       = (t_atm(ilay+1) + t_atm(ilay)) / 2.0_dp
                   delta_p = p_atm(ilay+1) - p_atm(ilay) 
                   delta_z = alt_atm(ilay) - alt_atm(ilay+1) 

                   if (kdis_species(isp).eq.'h2o') then
                      rho     = ((u_atm(ispatm,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      CALL CONT_H2O(anu,rho,delta_p,p,t,delta_z,tau_cont)

                   elseif (kdis_species(isp).eq.'o2') then
                      if (ispatm_h2o.eq.0) then
                         ! we need h2o concentration to get o2 continuum
                         DO i = 1, ngas_u_atm
                            if (gas_u_atm(i).eq.'h2o') then
                               ! the kdis specy has a defined concentration profil
                               ! ispatm is the index of the specy in atm profil arrays
                               ispatm_h2o  = i
                               flag_sp_h2o = .TRUE.
                               exit
                            endif
                         ENDDO
                      endif
                      IF (flag_sp_h2o.eqv..FALSE.) then
                         WRITE(*,*) ''
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         WRITE(*,*) ' (sub. set_gasabs)    ERROR '
                         WRITE(*,*) '   The gas concentration for h2o '
                         WRITE(*,*) '   is not defined.'
                         WRITE(*,*) '   It is needed to compute o2 continuum'
                         WRITE(*,*) '   (see the atmospheric profil file or '
                         WRITE(*,*) '   the uniformly distributed gas concentration file)'
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         WRITE(*,*) ''
                         STOP
                      ENDIF
                      rho     = ((u_atm(ispatm,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      rhoh2o  = ((u_atm(ispatm_h2o,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm_h2o,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      CALL CONT_O2(anu, p, t, rho, delta_z, rhoh2o, tau_cont)

                   elseif (kdis_species(isp).eq.'n2') then
                      if (ispatm_h2o.eq.0) then
                         ! we need h2o concentration to get n2 continuum
                         DO i = 1, ngas_u_atm
                            if (gas_u_atm(i).eq.'h2o') then
                               ! the kdis specy has a defined concentration profil
                               ! ispatm is the index of the specy in atm profil arrays
                               ispatm_h2o  = i
                               flag_sp_h2o = .TRUE.
                               exit
                            endif
                         ENDDO
                      endif
                      IF (flag_sp_h2o.eqv..FALSE.) then
                         WRITE(*,*) ''
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         WRITE(*,*) ' (sub. set_gasabs)    ERROR '
                         WRITE(*,*) '   The gas concentration for h2o '
                         WRITE(*,*) '   is not defined.'
                         WRITE(*,*) '   It is needed to compute n2 continuum'
                         WRITE(*,*) '   (see the atmospheric profil file or '
                         WRITE(*,*) '   the uniformly distributed gas concentration file)'
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         WRITE(*,*) ''
                         STOP
                      ENDIF
                      rho     = ((u_atm(ispatm,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      rhoh2o  = ((u_atm(ispatm_h2o,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm_h2o,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      CALL CONT_N2(anu, p, t, rho, delta_z, rhoh2o, tau_cont)

                   elseif (kdis_species(isp).eq.'co2') then
                      rho     = (u_atm(ispatm,ilay+1)+u_atm(ispatm,ilay)) / 2.0_dp
                      rho     = rho * delta_z *1.0d5 
                      CALL CONT_CO2(anu,rho,p,t,tau_cont)

                   elseif (kdis_species(isp).eq.'o3') then
                      rho     = (u_atm(ispatm,ilay+1)+u_atm(ispatm,ilay)) / 2.0_dp
                      rho     = rho * delta_z *1.0d5 
                      CALL CONT_O3(anu,rho,t,tau_cont)

                   elseif (kdis_species(isp).eq.'no2') then
                      rho     = (u_atm(ispatm,ilay+1)+u_atm(ispatm,ilay)) / 2.0_dp
                      rho     = rho * delta_z *1.0d5 
                      CALL CONT_NO2(anu,rho,tau_cont)

                   elseif (kdis_species(isp).eq.'so2') then
                      rho     = (u_atm(ispatm,ilay+1)+u_atm(ispatm,ilay)) / 2.0_dp
                      rho     = rho * delta_z *1.0d5 
                      CALL CONT_SO2(anu,rho,tau_cont)

                   else
                      WRITE(*,*) ' ' 
                      WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
                      WRITE(*,*) '  (sub set_gasabs)   ERROR            '
                      WRITE(*,*) '   No continuum definition for  ', kdis_species(isp)
                      WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
                      WRITE(*,*) ' ' 
                      STOP
                   endif

                   ! if (tau_cont.gt.0.0_dp) then
                   !    WRITE(*,*) ''
                   !    WRITE(*,*)                   '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   !    WRITE(*,*)                   ' (sub. set_gasabs)    Continuum '
                   !    WRITE(*,*)                   '   The opacity for ', ADJUSTL(TRIM(kdis_species(isp)))
                   !    WRITE(*,FMT='(A,2x,F10.5,2x,A)') '   at wavelength ', wlambda(ilamb), 'microns'
                   !    WRITE(*,FMT='(A,2x,I3)')         '   for layer ', ilay
                   !    WRITE(*,FMT='(A,2x,E18.9)')      '   Tau = ', tau_cont
                   !    WRITE(*,*)                   '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   !    WRITE(*,*) ''                      
                   ! end if

                   !if (tau_cont.lt.0.0_dp) then
                   !   WRITE(*,*) ''
                   !   WRITE(*,*)                   '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   !   WRITE(*,*)                   ' (sub. set_gasabs)    ERROR '
                   !   WRITE(*,*)                   '   The opacity for ', ADJUSTL(TRIM(kdis_species(isp)))
                   !   WRITE(*,FMT='(A,2x,F10.5,2x,A)') '   is negative at wavelength ', wlambda(ilamb), 'microns'
                   !   WRITE(*,FMT='(A,2x,I3)')         '   for layer ', ilay
                   !   WRITE(*,FMT='(A,2x,E18.9)')      '   Tau = ', tau_cont
                   !   WRITE(*,*)                   '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   !   WRITE(*,*) ''
                   !   STOP
                   !endif
                   ! the continuum apply to all ai
                   gl_gas_dtauabs_layers(ilay, isp_tmp, :, ilamb) = gl_gas_dtauabs_layers(ilay, isp_tmp, :, ilamb) &
                        + (kdis_fcont(isp)*tau_cont)
                end DO
             endif ! continuum
          END DO ! lambda
       enddo! isp


       ! same on species that have concentration dependency
       do isp = 1, kdis_nsp_c

          coldens_tmp(:) = 0.0_dp
          coldens_computed = .false.

          flag_sp = .FALSE.
          DO i = 1, ngas_u_atm
             if (kdis_species_c(isp).eq.gas_u_atm(i)) then
                ispatm  = i
                flag_sp = .TRUE.
                exit
             endif
          ENDDO
          IF (flag_sp.eqv..FALSE.) then
             WRITE(*,*) ''
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ' (sub. set_gasabs)    ERROR '
             WRITE(*,*) '   The gas concentration for k-dis gas specy ', &
                  TRIM(ADJUSTL( kdis_species_c(isp)))
             WRITE(*,*) '   is not defined.'
             WRITE(*,*) '   (see the atmospheric profil file or '
             WRITE(*,*) '   the uniformly distributed gas concentration file)'
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ''
             STOP
          ENDIF
          isp_tmp = isp_tmp + 1

          DO ilamb = 1, nlambda

             ikdislamb = lambda_ikdis(ilamb)
             gl_nai_layers(isp_tmp, ilamb) = kdis_nai_c(isp, ikdislamb)

             if (gl_nai_layers(isp_tmp, ilamb).eq.1) then

                gl_ai_layers(isp_tmp, 1, ilamb) = 1.0_dp
                if (kdis_xsect_c(isp, ikdislamb).gt.0.0_dp) then
                   if (coldens_computed .eqv. .false.) then
                      do ialt = 2, nalt_atm
                         coldens_tmp(ialt) = XINTEG2(1, ialt, nalt_atm, alt_atm*1.0D5, u_atm(ispatm,:) )
                         coldens_computed = .true.
                      end do
                   end if
                   do ilay = 1, nlayers
                      gl_gas_dtauabs_layers(ilay,isp_tmp,1,ilamb) = abs(coldens_tmp(ilay+1)-coldens_tmp(ilay)) &
                           & * kdis_xsect_c(isp, ikdislamb)
                   end do
                end if

             else

                DO iai = 1, gl_nai_layers(isp_tmp, ilamb)

                   gl_ai_layers(isp_tmp, iai, ilamb) = kdis_ai_c(isp, ikdislamb, iai)


                   DO ialt = 1, nalt_atm

                      if (p_atm(ialt).gt.maxval(kdis_p(:))) then
                         current_p = maxval(kdis_p(:))
                      else if (p_atm(ialt).lt.minval(kdis_p(:))) then
                         current_p = minval(kdis_p(:))
                      else
                         current_p = p_atm(ialt)
                      endif
                      
                      if (t_atm(ialt).gt.maxval(kdis_t(:))) then
                         current_t = maxval(kdis_t(:))
                      else if (t_atm(ialt).lt.minval(kdis_t(:))) then
                         current_t = minval(kdis_t(:))
                      else
                         current_t = t_atm(ialt)
                      endif
                      
                      if (u_atm(ispatm, ialt).gt.maxval(kdis_c(:))) then
                         current_c = maxval(kdis_c(:))
                      else if (u_atm(ispatm, ialt).lt.minval(kdis_c(:))) then
                         current_c = minval(kdis_c(:))
                      else
                         current_c = u_atm(ispatm, ialt)
                      endif

                      ki_tmp(ialt) = TRILININTPOL(kdis_ki_c(isp,ikdislamb,iai,:,:,:), kdis_p(:), kdis_t(:), kdis_c(:), &
                           kdis_np, kdis_nt, kdis_nc, current_p, current_t, current_c)

                   ENDDO

                   ! Compute integrated opacity
                   tau_tmp(:) = 0.0D0
                   DO ialt = 2, nalt_atm                      
                      !                                             z(cm)         
                      tau_tmp(ialt) = XINTEG2(1, ialt, nalt_atm, alt_atm*1.0D5,  ki_tmp * u_atm(ispatm,:) )
                   END DO

                   ! compute tau layers
                   DO ilay = 1, nlayers
                      gl_gas_dtauabs_layers(ilay, isp_tmp, iai, ilamb) = abs(tau_tmp(ilay+1) - tau_tmp(ilay))
                   END DO

                END DO
             endif
             ! ADD CONTINUUM CONTRIBUTION
             if (kdis_fcont_c(isp).gt.0.0_dp) then
                DO ilay = 1, nlayers
                   tau_cont = 0.0
                   anu      = 1.0d4/wlambda(ilamb)
                   !if ((anu.gt.50000.0).or.(anu.lt.200.0)) then
                   !   WRITE(*,*) 'The validity of continuum sub. outside GAME wavelength must be checked'
                   !   STOP
                   !endif
                   p       = (p_atm(ilay+1) + p_atm(ilay)) / 2.0_dp
                   t       = (t_atm(ilay+1) + t_atm(ilay)) / 2.0_dp
                   delta_p = p_atm(ilay+1) - p_atm(ilay) 
                   delta_z = alt_atm(ilay) - alt_atm(ilay+1) 

                   if (kdis_species_c(isp).eq.'h2o') then
                      rho     = ((u_atm(ispatm,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      CALL CONT_H2O(anu,rho,delta_p,p,t,delta_z,tau_cont)

                   elseif (kdis_species_c(isp).eq.'o2') then
                      if (ispatm_h2o.eq.0) then
                         ! we need h2o concentration to get o2 continuum
                         DO i = 1, ngas_u_atm
                            if (gas_u_atm(i).eq.'h2o') then
                               ! the kdis specy has a defined concentration profil
                               ! ispatm is the index of the specy in atm profil arrays
                               ispatm_h2o  = i
                               flag_sp_h2o = .TRUE.
                               exit
                            endif
                         ENDDO
                      endif
                      IF (flag_sp_h2o.eqv..FALSE.) then
                         WRITE(*,*) ''
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         WRITE(*,*) ' (sub. set_gasabs)    ERROR '
                         WRITE(*,*) '   The gas concentration for h2o '
                         WRITE(*,*) '   is not defined.'
                         WRITE(*,*) '   It is needed to compute o2 continuum'
                         WRITE(*,*) '   (see the atmospheric profil file or '
                         WRITE(*,*) '   the uniformly distributed gas concentration file)'
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         WRITE(*,*) ''
                         STOP
                      ENDIF
                      rho     = ((u_atm(ispatm,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      rhoh2o  = ((u_atm(ispatm_h2o,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm_h2o,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      CALL CONT_O2(anu, p, t, rho, delta_z, rhoh2o, tau_cont)

                   elseif (kdis_species_c(isp).eq.'n2') then
                      if (ispatm_h2o.eq.0) then
                         ! we need h2o concentration to get n2 continuum
                         DO i = 1, ngas_u_atm
                            if (gas_u_atm(i).eq.'h2o') then
                               ! the kdis specy has a defined concentration profil
                               ! ispatm is the index of the specy in atm profil arrays
                               ispatm_h2o  = i
                               flag_sp_h2o = .TRUE.
                               exit
                            endif
                         ENDDO
                      endif
                      IF (flag_sp_h2o.eqv..FALSE.) then
                         WRITE(*,*) ''
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         WRITE(*,*) ' (sub. set_gasabs)    ERROR '
                         WRITE(*,*) '   The gas concentration for h2o '
                         WRITE(*,*) '   is not defined.'
                         WRITE(*,*) '   It is needed to compute n2 continuum'
                         WRITE(*,*) '   (see the atmospheric profil file or '
                         WRITE(*,*) '   the uniformly distributed gas concentration file)'
                         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                         WRITE(*,*) ''
                         STOP
                      ENDIF
                      rho     = ((u_atm(ispatm,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      rhoh2o  = ((u_atm(ispatm_h2o,ilay+1)/ u_air_atm(ilay+1)) + &
                           &     (u_atm(ispatm_h2o,ilay)  / u_air_atm(ilay))  ) / 2.0_dp
                      CALL CONT_N2(anu, p, t, rho, delta_z, rhoh2o, tau_cont)

                   elseif (kdis_species_c(isp).eq.'co2') then
                      rho     = (u_atm(ispatm,ilay+1)+u_atm(ispatm,ilay)) / 2.0_dp
                      rho     = rho * delta_z *1.0d5 
                      CALL CONT_CO2(anu,rho,p,t,tau_cont)

                   elseif (kdis_species_c(isp).eq.'o3') then
                      rho     = (u_atm(ispatm,ilay+1)+u_atm(ispatm,ilay)) / 2.0_dp
                      rho     = rho * delta_z *1.0d5 
                      CALL CONT_O3(anu,rho,t,tau_cont)

                   elseif (kdis_species_c(isp).eq.'no2') then
                      rho     = (u_atm(ispatm,ilay+1)+u_atm(ispatm,ilay)) / 2.0_dp
                      rho     = rho * delta_z *1.0d5 
                      CALL CONT_NO2(anu,rho,tau_cont)

                   elseif (kdis_species_c(isp).eq.'so2') then
                      rho     = (u_atm(ispatm,ilay+1)+u_atm(ispatm,ilay)) / 2.0_dp
                      rho     = rho * delta_z *1.0d5 
                      CALL CONT_SO2(anu,rho,tau_cont)
                   else
                      WRITE(*,*) ' ' 
                      WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
                      WRITE(*,*) '  (sub set_gasabs)   ERROR            '
                      WRITE(*,*) '   No continuum definition for  ', kdis_species_c(isp)
                      WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
                      WRITE(*,*) ' ' 
                      STOP
                   endif
                   !if (tau_cont.lt.0.0_dp) then
                   !   WRITE(*,*) ''
                   !   WRITE(*,*)                   '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   !   WRITE(*,*)                   ' (sub. set_gasabs)    ERROR '
                   !   WRITE(*,*)                   '   The opacity for ', ADJUSTL(TRIM(kdis_species_c(isp)))
                   !   WRITE(*,FMT='(A,2x,F10.5,2x,A)') '   is negative at wavelength ', wlambda(ilamb), 'microns'
                   !   WRITE(*,FMT='(A,2x,I3)')         '   for layer ', ilay
                   !   WRITE(*,FMT='(A,2x,E18.9)')      '   Tau = ', tau_cont
                   !   WRITE(*,*)                   '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                   !   WRITE(*,*) ''
                   !   STOP
                   !endif
                   ! the continuum apply to all ai
                   gl_gas_dtauabs_layers(ilay, isp_tmp, :, ilamb) = gl_gas_dtauabs_layers(ilay, isp_tmp, :, ilamb) + &
                        (kdis_fcont_c(isp)*tau_cont)
                end DO
             endif ! continuum

          END DO ! lambda
       enddo! isp

    endif ! kdis mode

  END SUBROUTINE SET_GASABS


  ! ===============================================================================

  SUBROUTINE GET_MONO_LAYERS(ilamb)

    USE MCONSTANTS, only : dp, undef_i, undef_dp

    USE MCOMMON, only : nptcle, &
         nlayers, &
         sinsca_corint, &
         rt_model, &
         warning, &
         mode,    &
         nesp_layers, &
         nai_layers,  &
         ai_layers,   &
         tot_dtausca_layers, &
         ptcle_dtauabs_layers,&
         nbetal_layers, &
         betal_layers,  &
         n_aitaui_mono_layers, &
         ai_mono_layers,  &
         tot_dtau_mono_layers, &
         ssa_mono_layers, &
         nmax_aitaui_mono_layers, &
         gas_dtauabs_layers, &    ! gas_dtauabs_layers(nlayers, nesp_layers, naimax_layers, nlambda) 
         gas_dtausca_layers, &    ! gas_dtausca_layers(nlayers, nlambda)
         gas_dtau_mono_layers, &  ! gas_dtau_mono_layers(nlayers,nmax_aitaui_mono_layers) 
         gas_ssa_mono_layers, &   ! gas_ssa_mono_layers(nlayers, nmax_aitaui_mono_layers) 
         lambda_ikdis, &
         print_aitaui

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ilamb

    integer :: ilay
    integer :: iai_tmp
    integer :: isp_tmp
    integer :: ikdislamb 
    real(kind=dp) :: ai_tmp
    real(kind=dp) :: taui_tmp
    real(kind=dp) :: gl_gas_dtauabs_mono_layers(nlayers,nmax_aitaui_mono_layers)

    ! -----------
    
    if (mode.eq.'kdis') then

       ! kdis
       
       ikdislamb = lambda_ikdis(ilamb)

       n_aitaui_mono_layers = 1
       DO isp_tmp = 1, nesp_layers
          n_aitaui_mono_layers = n_aitaui_mono_layers * nai_layers(isp_tmp,ilamb)          
       END DO

       ! Here we perform the combinasion of the (ai,taui)
       ! from the different species to obtain
       ! the effective (ai, taui) 
       DO ilay = 1, nlayers
          iai_tmp   = 0
          isp_tmp   = 0
          ai_tmp    = 1.0_dp
          taui_tmp  = 0.0_dp
          call get_aitauieff_layers(nesp_layers, nai_layers(:,ilamb), ai_layers(:,:,ilamb), gas_dtauabs_layers(ilay,:,:,ilamb), &
               isp_tmp, iai_tmp, ai_tmp, taui_tmp, ai_mono_layers(:), gl_gas_dtauabs_mono_layers(ilay,:))
       ENDDO

       IF (SUM(ai_mono_layers(1:n_aitaui_mono_layers)).ne.1.0_dp) THEN
          if (warning) then
             WRITE(*,*) ''
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ' (sub. set_gasabs)    WARNING '
             WRITE(*,FMT='(A,1x,I4)') '   For k-dis band', ikdislamb
             WRITE(*,*) '   The weight of the effective (ai,ki) is not normalized to 1.0'
             WRITE(*,FMT='(A,1x,F20.14)') '   SUM(ai_eff) = ', SUM(ai_mono_layers(1:n_aitaui_mono_layers))
             WRITE(*,*) ' We renormalize it '
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) ''
          endif
          ai_mono_layers(1:n_aitaui_mono_layers) = ai_mono_layers(1:n_aitaui_mono_layers) &
               / SUM(ai_mono_layers(1:n_aitaui_mono_layers))
       ENDIF

    else
       
       ! mono

       n_aitaui_mono_layers            = 1
       ai_mono_layers(:)               = 1.0_dp
       gl_gas_dtauabs_mono_layers(:,1) = 0.0_dp

    endif

    IF ( (rt_model.eq.'disort') .or. &
         (rt_model.eq.'addoub') .or. &
         (rt_model.eq.'doad') ) THEN
       DO ilay = 1, nlayers  
          DO iai_tmp = 1, n_aitaui_mono_layers

             if (nptcle.gt.0) then
                tot_dtau_mono_layers(ilay,iai_tmp) = tot_dtausca_layers(ilay,ilamb) + &
                     &                               gl_gas_dtauabs_mono_layers(ilay,iai_tmp) + &
                     &                               ptcle_dtauabs_layers(ilay, ilamb)
             else
                tot_dtau_mono_layers(ilay,iai_tmp) = tot_dtausca_layers(ilay,ilamb) + &
                     &                               gl_gas_dtauabs_mono_layers(ilay,iai_tmp)
             endif


             if (tot_dtau_mono_layers(ilay,iai_tmp).eq.0.0_dp) then
                ssa_mono_layers(ilay,iai_tmp)      = 0.0_dp
             else
                ssa_mono_layers(ilay,iai_tmp)      = tot_dtausca_layers(ilay,ilamb) / &
                     &                               tot_dtau_mono_layers(ilay,iai_tmp) 
             endif

          ENDDO
       END DO
    ENDIF
    IF ( (rt_model .eq. 'mcrad1d') .or. &
         (rt_model .eq. 'sinsca') .or.  &
         (sinsca_corint).or.            &
         (print_aitaui) ) THEN
       DO ilay = 1, nlayers  
          DO iai_tmp = 1, n_aitaui_mono_layers
             gas_dtau_mono_layers(ilay,iai_tmp) = gl_gas_dtauabs_mono_layers(ilay,iai_tmp) + &
                  gas_dtausca_layers(ilay, ilamb)
             if (gas_dtau_mono_layers(ilay,iai_tmp) .eq. 0.0_dp) then
                gas_ssa_mono_layers(ilay,iai_tmp) = 0.0_dp
             else
                gas_ssa_mono_layers(ilay,iai_tmp) = gas_dtausca_layers(ilay, ilamb) / &
                     gas_dtau_mono_layers(ilay,iai_tmp)
             endif
          ENDDO
       END DO
    ENDIF

  END SUBROUTINE GET_MONO_LAYERS

  ! ===============================================================================

  recursive subroutine get_aitauieff_layers(nesp, nai_in, ai_in, taui_in, &
       isp_tmp, iai_tmp, ai_tmp, taui_tmp, ai_out, taui_out)

    ! With this routine we compute the effective ai,taui
    ! resulting from the combination of each gas species ai, taui
    ! see e.g. eq.31 in Lacis and Oinas, 1991

    USE MCONSTANTS, only : dp
    USE MCOMMON, only : naimax_layers, &
         nmax_aitaui_mono_layers

    IMPLICIT NONE

    INTEGER, INTENT(IN)          :: nesp
    !    INTEGER, INTENT(IN)          :: nai_in(nesp)
    !    REAL(kind=dp), INTENT(IN)    :: ai_in(nesp, naimax_layers)
    !    REAL(kind=dp), INTENT(IN)    :: taui_in(nesp, naimax_layers)
    INTEGER, INTENT(IN)          :: nai_in(:)
    REAL(kind=dp), INTENT(IN)    :: ai_in(:,:)
    REAL(kind=dp), INTENT(IN)    :: taui_in(:,:)
    integer, INTENT(INOUT)       :: isp_tmp
    integer, INTENT(INOUT)       :: iai_tmp
    REAL(kind=dp), INTENT(INOUT) :: ai_tmp
    REAL(kind=dp), INTENT(INOUT) :: taui_tmp
    !    REAL(kind=dp), INTENT(INOUT) :: taui_out(nmax_aitaui_mono_layers)
    !    REAL(kind=dp), INTENT(INOUT) :: ai_out(nmax_aitaui_mono_layers)
    REAL(kind=dp), INTENT(INOUT) :: taui_out(:)
    REAL(kind=dp), INTENT(INOUT) :: ai_out(:)

    INTEGER :: iai

    ! ------ 

    isp_tmp = isp_tmp + 1

    DO iai = 1, nai_in(isp_tmp)

       if (isp_tmp.lt.nesp)  then

          ai_tmp     = ai_tmp   * ai_in(isp_tmp,iai)
          taui_tmp   = taui_tmp + taui_in(isp_tmp,iai)

          call get_aitauieff_layers(nesp, nai_in, ai_in, taui_in, &
               isp_tmp, iai_tmp, ai_tmp, taui_tmp, ai_out, taui_out )

          ai_tmp     = ai_tmp   / ai_in(isp_tmp,iai)
          taui_tmp   = taui_tmp - taui_in(isp_tmp,iai)

       else

          iai_tmp           = iai_tmp + 1
          ai_out(iai_tmp)   = ai_tmp    * ai_in(isp_tmp,iai)
          taui_out(iai_tmp) = taui_tmp  + taui_in(isp_tmp,iai)
          ! print*, iai_tmp, ai_out(iai_tmp), taui_out(iai_tmp)

       endif

    END DO
    isp_tmp = isp_tmp - 1

  END subroutine get_aitauieff_layers


END MODULE MLAYERS
