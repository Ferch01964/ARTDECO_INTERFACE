!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  See reference for :
!!!     C. Cornet, L. C-Labonnote, F. Szczap,
!!      Using a 3-D radiative transfer Monte–Carlo model to assess
!!      radiative effects on polarized reflectances above cloud scenes
!!      Three-dimensional polarized Monte Carlo atmospheric radiative
!!      transfer model(3DMCPOL):3D effects on polarized visible
!!      reflectances of a cirrus cloud
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



SUBROUTINE mcrad1d(verbose, lambda_in, ngeom, thetas, thetav, phir, &
     nlay, nmu, mu, p_ptcle_ly_in,                                  &
     tau_ptcle_ly, ssa_ptcle_ly,                    &                                            
     tau_gas_ly, ssa_gas_ly,                        &                                            
     z_ref,                                         & 
     delta_ray,                                     &
     thermal, temp,                                 &
     Stoke_in,                                      &
     nphot,                                         &
     surface_in,                                    &
     surface_albedo,                                &
     nmat,                                          &
     nmax_interact,                                 &
     vrm_in, vrm_epsilon_ddis, vrm_n_firstcp,       &
     vrm_n_sccp, vrm_n_lecp,                        &
     vrm_wspl_crit,                                 &
     vrm_nsplmax,                                   &
     vrm_wrr_crit,                                  &
     vrm_wrr_min,                                   &
     Stoke_out, Stoke_out_sig)

  ! Stoke_in is a stoke vector for the incoming flux
  !          in the direction parallel to the beam

  !$ USE OMP_LIB
  USE MMCRAD1D_CONSTANTS
  USE MMCRAD1D_SUBMAIN
  USE MMCRAD1D_ZIGGURAT
  USE MMCRAD1D_BRDF
  USE MMCRAD1D_INTERPOLATE, only : MCLININTPOL

  ! ARDECO MODULE
  USE MSURFACE_BRDF

  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: verbose
  REAL(8), INTENT(IN) :: lambda_in
  INTEGER, INTENT(IN) :: ngeom 
  REAL(8), INTENT(IN) :: thetas
  REAL(8), INTENT(IN) :: thetav(ngeom)
  REAL(8), INTENT(IN) :: phir(ngeom)
  INTEGER, INTENT(IN) :: nlay
  INTEGER, INTENT(IN) :: nmu
  REAL(8), INTENT(IN) :: mu(nmu)
  REAL(8), INTENT(IN) :: p_ptcle_ly_in(nlay,6,nmu)
  !                      p_ptcle_ly_in(:,1,:)  = P11
  !                      p_ptcle_ly_in(:,2,:)  = P22
  !                      p_ptcle_ly_in(:,3,:)  = P12 = P21
  !                      p_ptcle_ly_in(:,4,:)  = P33
  !                      p_ptcle_ly_in(:,5,:)  = P44
  !                      p_ptcle_ly_in(:,6,:)  = P34 = -P43
  REAL(8), INTENT(IN) :: tau_ptcle_ly(nlay)
  REAL(8), INTENT(IN) :: ssa_ptcle_ly(nlay)
  REAL(8), INTENT(IN) :: tau_gas_ly(nlay)
  REAL(8), INTENT(IN) :: ssa_gas_ly(nlay)
  REAL(8), INTENT(IN) :: z_ref(0:nlay)
  REAL(8), INTENT(IN) :: delta_ray
  LOGICAL, INTENT(IN) :: thermal
  REAL(8), INTENT(IN) :: temp(0:nlay)
  REAL(8), INTENT(IN) :: Stoke_in(4)
  INTEGER(8), INTENT(IN) :: nphot
  CHARACTER(len=*), INTENT(IN) :: surface_in
  REAL(8), INTENT(IN) :: surface_albedo
  INTEGER, INTENT(IN) :: nmat
  INTEGER, INTENT(IN) :: nmax_interact
  INTEGER, INTENT(IN) :: vrm_in
  REAL(8), INTENT(IN) :: vrm_epsilon_ddis
  INTEGER, INTENT(IN) :: vrm_n_firstcp
  INTEGER, INTENT(IN) :: vrm_n_sccp
  INTEGER, INTENT(IN) :: vrm_n_lecp
  REAL(8), INTENT(IN) :: vrm_wspl_crit
  INTEGER, INTENT(IN) :: vrm_nsplmax                                       
  REAL(8), INTENT(IN) :: vrm_wrr_crit                                    
  REAL(8), INTENT(IN) :: vrm_wrr_min                                       
  REAL(8), INTENT(OUT) :: Stoke_out(ngeom, nmat)
  REAL(8), INTENT(OUT) :: Stoke_out_sig(ngeom, nmat)

  !----- local variables

  CHARACTER(len=20) :: surface
  integer :: nelem

  logical :: ptcle_flag

  ! Atmosphere definition
  real(8), allocatable :: p_ptcle_ly(:,:,:) ! p_ptcle_ly(nlay,nelem,nmu) 
  real(8)              :: thetad(nmu)
  real(8)              :: prf_ptcle_ly(nlay,nmu)
  real(8) :: cext_ly(nlay)
  real(8) :: ssa_ly(nlay)
  real(8) :: cext_ptcle_ly(nlay)
  real(8) :: cext_gas_ly(nlay)
  real(8) :: fracext_gas_ly(nlay)
  real(8) :: tau_tot

  !  VRM 
  integer  :: vrm
  REAL(8)  :: rn_rr
  real(8)  :: p_ptcle_ddis(nmu) 
  real(8)  :: prf_ptcle_ddis(nmu)
  logical  :: ddis_mp
  integer  :: i_cp

  ! geometry
  real(8) :: cos_dir_sensor(ngeom, 3)
  real(8) :: costhetas, sinthetas

  ! for Ziggurate RNG in parallel
  integer, allocatable :: seed(:)
  integer :: grainsize = 32
  real(8) :: r
  integer :: nthreads
  integer :: nproc
  integer :: ithread

  ! single photon characteristics
  real(8) :: Stk(nmat)
  real(8) :: cos_dir(3)
  real(8) :: u(3)
  real(8) :: v(3)
  real(8) :: wght
  real(8) :: z_phot
  integer :: ilayer
  real(8) :: d
  ! NOTE   : count_scat is the total number of scattering including reflection by the surface
  integer :: count_scat
  ! NOTE   : count_reflect is the number of surface reflection
  integer :: count_reflect
  integer :: leave

  ! cloned photon characteristics
  REAL(8) :: Stk_cp(nmat)
  REAL(8) :: wght_cp
  REAL(8) :: cos_dir_cp(3)
  REAL(8) :: u_cp(3)
  REAL(8) :: v_cp(3)
  REAL(8) :: z_phot_cp
  INTEGER :: ilayer_cp
  INTEGER :: leave_cp
  INTEGER :: count_scat_cp
  REAL(8) :: wspl_cp

  logical :: flag_stop

  integer :: ind
  integer :: ind_cp

  REAL(8) :: le(nmat)

  integer :: split_level

  ! various counters
  integer(8) :: leave_bot, leave_top
  integer(8) :: interact, reflect, vanished
  integer(8) :: tot_interact
  integer(8) :: tot_reflect

  !  to compute statistics
  integer, parameter :: npack = 20
  integer(8) :: nphot_pack
  real(8), dimension(npack,ngeom,nmat) :: Stoke_pack
  integer(8), dimension(npack) :: count_phot_pack
  integer :: ipack
  real(8) :: stoke_out_mean(ngeom, nmat)
  !real(8) :: sig_one(ngeom, nmat)

  ! BRDF definition
  real(8), allocatable :: L_int_surf(:,:)
  real(8), allocatable :: L_int_teta(:)

  integer(8) :: iphot
  integer :: k, i, j
  integer :: igeom
  REAL(8) :: summ

  !-------------------------------
  !--- Initialization of some vectors and variables

  nphot_pack         = nphot / npack
  Stoke_pack(:,:,:)  = 0.0D0
  count_phot_pack(:) = 0
  leave_bot  = 0
  leave_top  = 0
  interact   = 0
  vanished   = 0
  reflect    = 0
  tot_interact = 0
  tot_reflect  = 0
  costhetas = cos(thetas*DTR)
  sinthetas = sin(thetas*DTR)
  ! get nelem
  IF (nmat .EQ. 1) THEN
     nelem = 1
  ELSE IF (nmat .EQ. 3) THEN
     nelem = 4
  ELSE IF (nmat .EQ. 4) THEN
     nelem = 6
  ELSE
     WRITE(*,*)'(mcrad1d) ERROR invalid value for nmat:',nmat
     STOP
  ENDIF

  ! test the presence of particles in the atmosphere
  ptcle_flag = .false.
  do k = 1, nlay
     if (tau_ptcle_ly(k) .gt. 0.0D0) ptcle_flag = .true.
  enddo

  ! mu must be sorted in decreasing order
  do k = 1, nmu-1
     if (mu(k) .le. mu(k+1)) then
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) '  MCRAD1D : ERROR'
        write(*,*) '      mu must be sorted in decreasing order'
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        STOP
     endif
  enddo

  ! compute Cext_ly for each layers (km-1)
  tau_tot = 0.0D0
  do k = 1, nlay
     cext_ly(k)       = (tau_ptcle_ly(k) + tau_gas_ly(k)) / (z_ref(k)-z_ref(k-1)) 
     cext_gas_ly(k)   = tau_gas_ly(k)    / (z_ref(k)-z_ref(k-1)) 
     cext_ptcle_ly(k) = tau_ptcle_ly(k)  / (z_ref(k)-z_ref(k-1)) 
     ssa_ly(k)        = (tau_ptcle_ly(k)*ssa_ptcle_ly(k) + tau_gas_ly(k)*ssa_gas_ly(k)) / &
          (tau_ptcle_ly(k) + tau_gas_ly(k))
     fracext_gas_ly(k) = cext_gas_ly(k) /  cext_ly(k)
     tau_tot           = tau_tot + (cext_ly(k)* (z_ref(k)-z_ref(k-1)) )
  enddo

  ALLOCATE(p_ptcle_ly(nlay,nelem,nmu))
  p_ptcle_ly(:,:,:) = 0.0D0
  if (ptcle_flag) then
     thetad(:) = ACOS(mu(:))/DTR
     ! get ptcle PRFs
     call get_prfs(nmu, nlay, mu, p_ptcle_ly_in(:,1,:), prf_ptcle_ly)
     if (nelem .ge.1) then
        p_ptcle_ly(:,1,:) = p_ptcle_ly_in(:,1,:)                          ! P11
     endif
     if (nelem .ge. 4) then
        p_ptcle_ly(:,2,:) = p_ptcle_ly_in(:,2,:) / p_ptcle_ly_in(:,1,:)   ! P22 / P11
        p_ptcle_ly(:,3,:) = p_ptcle_ly_in(:,3,:) / p_ptcle_ly_in(:,1,:)   ! P12 = P21 / P11
        p_ptcle_ly(:,4,:) = p_ptcle_ly_in(:,4,:) / p_ptcle_ly_in(:,1,:)   ! P33 / P11
     endif
     if (nelem .eq. 6) then
        p_ptcle_ly(:,5,:) = p_ptcle_ly_in(:,5,:) / p_ptcle_ly_in(:,1,:)   ! P44 / P11
        p_ptcle_ly(:,6,:) = p_ptcle_ly_in(:,6,:) / p_ptcle_ly_in(:,1,:)   ! P34 = -P43 / P11
     endif
     ! make sure we have no crazy value in layers where no particles are present
     DO  k = 1, nlay
        IF (fracext_gas_ly(k) .EQ. 1.0D0) THEN 
           p_ptcle_ly(k,:,:) = 0.0d0
           prf_ptcle_ly(k,:) = 0.0d0
        ENDIF
     ENDDO
  ENDIF

  ! ------------------- VRM ------------------ 
  ! no VRM needs to be applied if no particles present in the atm
  if (ptcle_flag) then
     vrm = vrm_in
  else
     vrm = 0
  endif

  if ( (vrm .eq. 1) .and. (vrm_n_firstcp .gt. 1) ) then
     ddis_mp = .true.
  else
     ddis_mp = .false.
  endif
  ! get a modified Phase function for DDISing      
  ! and set VRM variables
  if (vrm .ne. 0) then
     ! compute the p_ddis and prf_ddis only among particle phase functions (Rayleigh excluded)
     ! Define Pddis as specified in Buras and Mayer, 2011, eq. 10
     ! Set it to the P11 of the first layer
     p_ptcle_ddis(:) = P_ptcle_ly(1,1,:)
     ! then look to the other layers for the highest value
     IF (nlay .gt. 1) then
        DO  k = 2, nlay
           IF (fracext_gas_ly(k) .LT. 1.0D0) THEN 
              DO i = 1, nmu
                 if (p_ptcle_ddis(i) .lt. P_ptcle_ly(k,1,i)) p_ptcle_ddis(i) = P_ptcle_ly(k,1,i)
              ENDDO
           ENDIF
        ENDDO
     ENDIF
     ! renormalize Pddis
     summ = -MC_XINTEG2(1, nmu, nmu, mu, p_ptcle_ddis)
     p_ptcle_ddis = p_ptcle_ddis / summ * 2.0D0
     prf_ptcle_ddis(1) = 0.0D0
     DO k = 2, nmu
        prf_ptcle_ddis(k) =  -MC_XINTEG2(1, k, nmu, mu, p_ptcle_ddis) / 2.0D0
     ENDDO
     !DO i = 1, nmu
     ! write(*,*) thetad(i), p_ptcle_ddis(i), prf_ptcle_ddis(i)
     !ENDDO
     if ((ngeom .gt. 1) .and. (ddis_mp .eqv. .true.)) then
        WRITE(*,*) ''
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' MCRAD1D ERROR :'
        WRITE(*,*) ' VRM for multiple geom. require vrm_n_firstcp = 1'
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ''
        STOP
     endif
     if ( vrm_n_sccp .lt. vrm_n_lecp ) then
        WRITE(*,*) ''
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' MCRAD1D ERROR :'
        WRITE(*,*) ' vrm_n_sccp must be greater or equal to vrm_n_lecp '
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ''
        STOP
     endif
     if ( vrm_wrr_crit .le. vrm_wrr_min ) then
        WRITE(*,*) ''
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' MCRAD1D ERROR :'
        WRITE(*,*) ' in VRM = 2 (N-TULPE) : vrm_wrr_crit must be > vrm_wrr_min'
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ''
        STOP
     endif
     if ( vrm_wspl_crit .le. vrm_wrr_crit ) then
        WRITE(*,*) ''
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' MCRAD1D ERROR :'
        WRITE(*,*) ' in VRM = 2 (N-TULPE) : vrm_wspl_crit must be > vrm_wrr_crit'
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ''
        STOP
     endif
     if ( vrm_wrr_crit .gt. 1 ) then
        WRITE(*,*) ''
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' MCRAD1D ERROR :'
        WRITE(*,*) ' in VRM = 2 (N-TULPE) : vrm_wrr_crit must be < 1'
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ''
        STOP
     endif
     if ( vrm_wrr_min .gt. 1 ) then
        WRITE(*,*) ''
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' MCRAD1D ERROR :'
        WRITE(*,*) ' in VRM = 2 (N-TULPE) : vrm_wrr_min must be < 1'
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ''
        STOP
     endif
     if ( vrm_wspl_crit .lt. 1) then
        WRITE(*,*) ''
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' MCRAD1D ERROR :'
        WRITE(*,*) ' in VRM = 2 (N-TULPE) : vrm_wspl_crit must be > 1'
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ''
        STOP
     endif
     if ( ptcle_flag .eqv. .false.) then
        WRITE(*,*) ''
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' MCRAD1D ERROR :'
        WRITE(*,*) ' VRM does not work if no ptcle in the atmosphere'
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ''
        STOP
     endif
     if ( vrm_n_firstcp .gt. 1 .and. (ddis_mp .eqv. .false.)) then
        WRITE(*,*) ''
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' MCRAD1D WARNING :'
        WRITE(*,*) '  VRM is inefficient if vrm_n_firstcp > 1 and  '
        WRITE(*,*) '  ngeom > 1 (cause no DDISing is applied to MP in that case)'
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ''
        STOP
     endif

  endif

  if ((surface_in.eq.'lambert') .and. (surface_albedo.lt.1d-50)) then
     surface = 'none'
     if (verbose) then
        write(*,*) ' '
        write(*,*) '        (mcrad1d) no surface'
     endif
  else
     surface = TRIM(ADJUSTL(surface_in))
  endif
  !----------------------
  ! set surface BRDF variables
  if (surface .eq. 'brdf') then

     if (verbose) then
        write(*,*) ' '
        write(*,*) '        (mcrad1d) setting-up BRDF variables'
     endif
     ! we set the number of phi and theta grid points 
     nphisurf  = 181
     ntetasurf = 45
     ALLOCATE( phisurf(nphisurf), tetasurf(ntetasurf),brdf(nmat,nmat,ntetasurf,ntetasurf,nphisurf), &
          L_int_surf(ntetasurf,ntetasurf), prf_phi_surf(ntetasurf,ntetasurf, nphisurf), &
          prf_mu_surf(ntetasurf,ntetasurf),alb_surf(ntetasurf), L_int_teta(nphisurf) )
     ! compute angles to define the BRDF on
     dtetasurf = 90.0D0 / (ntetasurf-1)
     dphisurf  = 360.0D0 / (nphisurf-1)
     do i = 1, ntetasurf
        tetasurf(i) = (i-1) * dtetasurf
     enddo
     do i = 1, nphisurf
        phisurf(i)  = (i-1) * dphisurf
     enddo
     ! compute the BRDF 
     ! and the integrated reflectivity on phi
     if (verbose) then
        write(*,*) ' '
        write(*,*) '        (mcrad1d) --> compute BRDF'
     endif
     do i = 1, ntetasurf   
        do j = 1, ntetasurf 
           do k = 1, nphisurf   ! phi
              ! we put the result in brdf(:,:,i,j,nphisurf-k+1) soo that phi is inverted regarding the
              ! surface_brdf convention
              CALL surface_brdf(nmat, lambda_in, tetasurf(i), tetasurf(j), 180.0+phisurf(k), brdf(:,:,i,j,nphisurf-k+1))
           enddo
        enddo
     enddo

     if (verbose) then
        write(*,*) ' '
        write(*,*) '        (mcrad1d) --> compute - prf_phi_surf(theta_inc, theta_refl, phi)'
        write(*,*) '                              - prf_mu_surf(theta_inc, theta_refl)'
        write(*,*) '                              - alb_surf(theta_inc)'
     endif
     prf_phi_surf(:,:,:) = 0.0D0
     prf_mu_surf(:,:)    = 0.0D0
     L_int_surf(:,:)     = -32768.0D0
     alb_surf(:)         = -32768.0D0
     theta_incident : do i = 1, ntetasurf

        ! compute cumulative probability distributions
        do j = 1, ntetasurf   
           do k = 2, nphisurf
              prf_phi_surf(i,j,k) = MC_XINTEG2(1,k,nphisurf,phisurf*DTR,brdf(1,1,i,j,:))
           enddo
           L_int_surf(i,j)     = prf_phi_surf(i,j,nphisurf)
           prf_phi_surf(i,j,:) = prf_phi_surf(i,j,:) / L_int_surf(i,j) ! prf_phi_surf renormalization
        enddo
        do j = 2, ntetasurf 
           prf_mu_surf(i,j) = -MC_XINTEG2(1,j,ntetasurf,COS(tetasurf*DTR),COS(tetasurf*DTR)*L_int_surf(i,:))
        enddo
        prf_mu_surf(i,:) = prf_mu_surf(i,:) / prf_mu_surf(i,ntetasurf) ! prf_mu_surf renormalization

        ! compute the albedo
        L_int_teta(:) = 0.0D0
        do j = 1, nphisurf
           L_int_teta(j) = -MC_XINTEG2(1,ntetasurf,ntetasurf,COS(tetasurf*DTR),COS(tetasurf*DTR)*brdf(1,1,i,:,j))
        enddo
        alb_surf(i)      =  MC_XINTEG2(1,nphisurf,nphisurf,phisurf*DTR,L_int_teta(:)) / pi 
        !write(*,*) tetasurf(i), alb_surf(i) 

        ! BRDF renormalisation 
        brdf(:,:,i,:,:) =  brdf(:,:,i,:,:) / alb_surf(i)

     enddo theta_incident

     DEALLOCATE(L_int_surf,L_int_teta)
  endif

  !write(*,*) nphot

  if (verbose) then
     write(*,*) ' '
     write(*,*) '        (mcrad1d) start to through photons'
  endif

  ! define the number of threads that will be running in parallel
  !$ nproc = omp_get_num_procs() 
  nthreads = 1 ! so that nthreads = 1 if we do not compile the code with openmp option (for the ziggurat RNG initialisation)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$ nthreads = nproc 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$ call OMP_SET_NUM_THREADS(nthreads)
  !  Initialize random number function
  allocate(seed(nthreads))
  CALL INIT_RANDOM_SEED()
  do k = 1,nthreads
     call random_number(r)
     seed(k) = INT(123456789*r)
  enddo
  call par_zigset(nthreads, seed, grainsize)

  !   Compute the sensor direction for each geometry
  DO igeom = 1,ngeom
     cos_dir_sensor(igeom,1) = SIN(DTR*thetav(igeom))*COS(DTR*phir(igeom))
     cos_dir_sensor(igeom,2) = SIN(DTR*thetav(igeom))*SIN(DTR*phir(igeom))
     cos_dir_sensor(igeom,3) = COS(DTR*thetav(igeom))
  ENDDO !igeom

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !                       Loop on photons number
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !$OMP PARALLEL                                                   &                
  !$OMP DEFAULT(PRIVATE)                                           &
  !$OMP SHARED(ngeom, thetas, thetav, phir, nelem)                 &
  !$OMP SHARED(nlay, nmu, mu, prf_ptcle_ly,thetad, p_ptcle_ly)     &  
  !$OMP SHARED(tau_ptcle_ly, ssa_ptcle_ly,tau_gas_ly, ssa_gas_ly)   &      
  !$OMP SHARED(z_ref, delta_ray, thermal, temp, Stoke_in)           &                                   
  !$OMP SHARED(nphot,nmat)                                          &
  !$OMP SHARED(leave_bot, leave_top, cos_dir_sensor)                &
  !$OMP SHARED(cext_ly,cext_ptcle_ly)                               &
  !$OMP SHARED(cext_gas_ly,ssa_ly,fracext_gas_ly)                   &
  !$OMP SHARED(costhetas, sinthetas, interact)                      &
  !$OMP SHARED(surface, surface_albedo, tau_tot, vanished)          &
  !$OMP SHARED(nphot_pack,Stoke_pack,count_phot_pack)               &
  !$OMP SHARED(ddis_mp, vrm_epsilon_ddis,  p_ptcle_ddis, prf_ptcle_ddis)  &
  !$OMP SHARED(nphisurf,ntetasurf,phisurf,tetasurf)                       &
  !$OMP SHARED(prf_phi_surf, prf_mu_surf, alb_surf)                       &
  !$OMP SHARED(dphisurf, dtetasurf, brdf)                                 &
  !$OMP SHARED(vrm, vrm_wspl_crit, vrm_n_firstcp, vrm_n_sccp)             &
  !$OMP SHARED(tot_interact, tot_reflect, reflect, nmax_interact)         &  
  !$OMP SHARED(vrm_nsplmax, vrm_wrr_crit, vrm_wrr_min, vrm_n_lecp)      

  ithread    = 0
  !$ ithread = omp_get_thread_num()      

  !$OMP DO SCHEDULE(STATIC)       &
  !$OMP REDUCTION(+:vanished)     &
  !$OMP REDUCTION(+:leave_top)    &
  !$OMP REDUCTION(+:leave_bot)    &
  !$OMP REDUCTION(+:interact)     &
  !$OMP REDUCTION(+:Stoke_pack)   &
  !$OMP REDUCTION(+:tot_interact) &
  !$OMP REDUCTION(+:tot_reflect)  &
  !$OMP REDUCTION(+:reflect)      &
  !$OMP REDUCTION(+:count_phot_pack) 
  photon: DO iphot = 1, nphot

     ! choose in which packet to store the photon contribution (to get statistics)
     ipack = 1
     DO WHILE(((ipack*nphot_pack) .lt. iphot) .and. (ipack .lt. npack))
        ipack=ipack+1
     ENDDO
     count_phot_pack(ipack) = count_phot_pack(ipack) + 1 

     ! create a new photon (called "Mother Photon" in case of N-tulpe)
     call get_new_photon(thermal, sinthetas, costhetas, Stoke_in, nmat, nlay, z_ref, &
          Stk, wght, cos_dir, u, v, z_phot, ilayer, d, leave, count_scat,            &
          count_reflect)

     IF (vrm .eq. 0) THEN

        !***************************************
        !***************************************
        !*********    NO  VRM          *********
        !***************************************
        !***************************************

        flag_stop = .false.

        ! Start following the photon until it leaves the atmosphere
        ! or it vanished (due to absorption)          
        DO WHILE ((leave .EQ. 0) .AND. (wght .gt. min_wght) .and. (flag_stop .eqv. .false.))

           ! Propagate the photon to its new location
           ! determine the kind of target that stops it
           CALL MOVE_PHOTON(ithread, surface,          &
                nlay, z_ref, cext_ly, fracext_gas_ly,  &
                cos_dir, ilayer, z_phot, leave,        &
                ind)

           IF (LEAVE .EQ. 0) THEN ! photon interacts

              count_scat = count_scat + 1
              if (ind .eq. 2) count_reflect = count_reflect + 1
              if ( (nmax_interact .ne. -1) .and. & 
                   (count_scat .ge. nmax_interact) ) flag_stop = .true.

              !---  Modifiy the photon weight to account for absorption
              IF (ind .EQ. 0) THEN
                 wght = wght * ssa_gas_ly(ilayer)                  
              ELSE IF (ind .EQ. 1) THEN
                 wght = wght * ssa_ptcle_ly(ilayer)
              ELSE IF ((ind .EQ. 2) .AND. (TRIM(ADJUSTL(surface)) .eq. 'lambert')) then
                 wght = wght * surface_albedo                    
              ELSE IF ((ind .EQ. 2) .AND. (TRIM(ADJUSTL(surface)) .eq. 'brdf')) then
                 ! we need to interpolate to get the value of the albedo for the incident theta on the BRDF
                 wght = wght * MCLININTPOL(alb_surf, tetasurf, ntetasurf,  ACOS(-(cos_dir(3)))/DTR)
              ENDIF

              ! LOCAL ESTIMATE
              DO igeom = 1, ngeom
                 CALL local_estimate(nmat, nelem, surface, nlay, cext_ly, z_ref, tau_tot, &
                      ind, nmu ,thetad, p_ptcle_ly(ilayer,:,:), delta_ray,                &
                      wght, Stk, z_phot, cos_dir, v, ilayer,                              &  
                      thetav(igeom), cos_dir_sensor(igeom, :), le)
                 Stoke_pack(ipack, igeom, :) = Stoke_pack(ipack, igeom, :) + le(:)
              ENDDO

              ! Compute new stoke vector, new cos_dir, u and v after scattering
              if (flag_stop .eqv. .false.) then
                 CALL scatter(ithread, nmat, nelem, surface,                                           &  
                      .false., vrm_epsilon_ddis,  p_ptcle_ddis, prf_ptcle_ddis,                        &
                      cos_dir_sensor(1, :),                                                            &
                      ind, nmu, mu, thetad, p_ptcle_ly(ilayer,:,:), prf_ptcle_ly(ilayer,:), delta_ray, &
                      wght, Stk, v, u, cos_dir)
              endif

           ENDIF

        ENDDO ! stop following the photon

        tot_interact = tot_interact + count_scat
        if (count_scat .ne. 0) interact = interact + 1
        tot_reflect = tot_reflect + count_reflect
        if (count_reflect .ne. 0) reflect = reflect + 1
        IF (wght .le. min_wght) THEN
           vanished = vanished + 1
        ENDIF
        IF (leave .EQ. 2) THEN
           leave_top = leave_top + 1
        ELSE IF (leave .EQ. 1) THEN
           leave_bot = leave_bot + 1
        ENDIF

        !***************************************
        !***************************************
        !*******  END   NO  VRM          *******
        !***************************************
        !***************************************

     ELSE

        !***************************************
        !***************************************
        !*********        VRM          *********
        !***************************************
        !***************************************

        split_level = 0 ! test variable that must be initialized to 0
        i_cp        = 0 ! init CP counter per MP

        ! Start following the Mother Photon     
        DO WHILE ( leave .EQ. 0 )  

           ! Propagate the MP to its new location
           ! determine the kind of target that stops it
           CALL MOVE_PHOTON(ithread, surface,          &
                nlay, z_ref, cext_ly, fracext_gas_ly,  &
                cos_dir, ilayer, z_phot, leave,        &
                ind)

           IF (LEAVE .EQ. 0) THEN ! MP interacts

              ! count scat is the number of interaction the MP undergone (including surface reflection)
              count_scat = count_scat + 1

              !---  Modifiy the photon weight to account for absorption
              IF (ind .EQ. 0) THEN
                 wght = wght * ssa_gas_ly(ilayer)                  
              ELSE IF (ind .EQ. 1) THEN
                 wght = wght * ssa_ptcle_ly(ilayer)
              ELSE IF ((ind .EQ. 2) .AND. (TRIM(ADJUSTL(surface)) .eq. 'lambert')) then
                 wght = wght * surface_albedo   
              ELSE IF ((ind .EQ. 2) .AND. (TRIM(ADJUSTL(surface)) .eq. 'brdf')) then
                 ! we need to interpolate to get the value of the albedo for the incident theta on the BRDF
                 wght = wght * MCLININTPOL(alb_surf, tetasurf, ntetasurf,  ACOS(-(cos_dir(3)))/DTR)
              ENDIF

              ! LOCAL ESTIMATE on MP are performed for its vrm_n_firstcp first interactions
              if (count_scat .le. vrm_n_firstcp) then 
                 DO igeom = 1, ngeom
                    CALL local_estimate(nmat, nelem, surface, nlay, cext_ly, z_ref, tau_tot,  &
                         ind, nmu ,thetad, p_ptcle_ly(ilayer,:,:), delta_ray, &
                         wght, Stk, z_phot, cos_dir, v, ilayer,               & 
                         thetav(igeom), cos_dir_sensor(igeom, :), le)
                    Stoke_pack(ipack, igeom, :) = Stoke_pack(ipack, igeom, :) + le(:)
                    !write(*,*) 'MP     ', count_scat, le(1)
                 ENDDO
              endif

              ! create a CP photon 
              IF ( count_scat .eq. (vrm_n_firstcp + i_cp*vrm_n_lecp) ) THEN

                 i_cp = i_cp + 1
                 !write(*,*) 'create CP number', i_cp

                 ! we create a CP for each geometry because
                 ! DDIsing is (of course) geometry dependant
                 DO igeom = 1, ngeom

                    ! 1 - We create a cloned photon that will be DDISED
                    !     for the given geometry 
                    ! 2 - We perform the local estimate for that geometry
                    Stk_cp           = stk
                    wght_cp          = wght 
                    cos_dir_cp       = cos_dir
                    u_cp             = u
                    v_cp             = v
                    z_phot_cp        = z_phot
                    ilayer_cp        = ilayer
                    leave_cp         = leave
                    ind_cp           = ind
                    count_scat_cp    = 0

                    ! compute new stoke vector, new cos_dir, u and v after scattering (or reflection)
                    ! this happen from the same location than MP 
                    ! NB : CP is always DDISed 
                    CALL scatter(ithread, nmat, nelem, surface,                                                    &  
                         .true., vrm_epsilon_ddis, p_ptcle_ddis, prf_ptcle_ddis,                                   &
                         cos_dir_sensor(igeom, :),                                                                 &
                         ind_cp, nmu, mu, thetad, p_ptcle_ly(ilayer_cp,:,:), prf_ptcle_ly(ilayer_cp,:), delta_ray, &
                         wght_cp, Stk_cp, v_cp, u_cp, cos_dir_cp)
                    ! this first scattering from the MP position is accounted for
                    count_scat_cp = count_scat_cp + 1

                    ! then CP starts to live its own life...
                    DO WHILE ((leave_cp .EQ. 0) .AND. (count_scat_cp .le. vrm_n_sccp))   

                       ! compute Wspl
                       CALL get_wspl(nmu, thetad, p_ptcle_ddis, &
                            wght_cp, cos_dir_cp,                &      
                            cos_dir_sensor(igeom, :),           &
                            wspl_cp)

                       IF (wspl_cp .gt. vrm_wspl_crit) THEN

                          !write (*,*) ' split photon = ', wspl_cp
                          CALL split_photon(ithread, nmat, nelem, surface, surface_albedo, tau_tot,          & 
                               .true., vrm_epsilon_ddis, vrm_wspl_crit, vrm_nsplmax, vrm_n_sccp, vrm_n_lecp, &
                               vrm_wrr_crit, vrm_wrr_min, &    
                               nmu, mu, nlay, thetad, p_ptcle_ddis, prf_ptcle_ddis, p_ptcle_ly, prf_ptcle_ly, delta_ray,    &
                               z_ref, cext_ly, fracext_gas_ly, ssa_ptcle_ly, ssa_gas_ly,                                    &
                               cos_dir_sensor(igeom,:), thetav(igeom),                                                      &
                               i_cp, Stk_cp, cos_dir_cp, u_cp, v_cp, wght_cp, z_phot_cp, ilayer_cp, count_scat_cp, &
                               leave_cp, ind_cp, MIN(INT(wspl_cp), vrm_nsplmax),                   &
                               nmax_interact, count_scat,                                          &
                               split_level,                                                                        &
                               Stoke_pack(ipack, igeom, :))

                          leave_cp = 3

                       ELSE IF (  (wspl_cp .lt. vrm_wrr_crit)            .and. & 
                            ((wght_cp / MAX(wspl_cp, vrm_wrr_min)) .le. 1.0d0) ) THEN

                          !*********************
                          !  Russion roulette
                          !CALL RANDOM_NUMBER(rn_rr)
                          rn_rr = par_ziguni(ithread)
                          IF (rn_rr .le. (1.0D0 - MAX(wspl_cp, vrm_wrr_min)) ) THEN
                             ! kill the CP photon
                             !write(*,*) ' kill CP', &
                             !     wspl_cp, rn_rr, (1.0D0 - MAX(wspl_cp, vrm_wrr_min))
                             leave_cp = 3
                          ELSE
                             ! Increase the photon weight
                             wght_cp = wght_cp / MAX(wspl_cp, vrm_wrr_min)
                             ! Propagate the CP photon to its new location
                             ! determine the kind of target that stops the photon
                             CALL MOVE_PHOTON(ithread, surface,               &
                                  nlay, z_ref, cext_ly, fracext_gas_ly,       &
                                  cos_dir_cp, ilayer_cp, z_phot_cp, leave_cp, &
                                  ind_cp)
                          ENDIF
                          !*********************

                       ELSE

                          ! Propagate the CP photon to its new location
                          ! determine the kind of target that stops the photon
                          CALL MOVE_PHOTON(ithread, surface,               &
                               nlay, z_ref, cext_ly, fracext_gas_ly,       &
                               cos_dir_cp, ilayer_cp, z_phot_cp, leave_cp, &
                               ind_cp)

                       ENDIF

                       IF (leave_cp .eq. 0) THEN  ! CP interacts

                          count_scat_cp = count_scat_cp + 1                         

                          !---  Modifiy the CP weight to account for absorption
                          IF (ind_cp .EQ. 0) THEN
                             wght_cp = wght_cp * ssa_gas_ly(ilayer_cp)                  
                          ELSE IF (ind_cp .EQ. 1) THEN
                             wght_cp = wght_cp * ssa_ptcle_ly(ilayer_cp)
                          ELSE IF ((ind_cp .EQ. 2) .AND. (TRIM(ADJUSTL(surface)) .eq. 'lambert')) then
                             wght_cp = wght_cp * surface_albedo   
                          ELSE IF ((ind_cp .EQ. 2) .AND. (TRIM(ADJUSTL(surface)) .eq. 'brdf')) then
                             ! we need to interpolate to get the value of the albedo for the incident theta on the BRDF
                             wght_cp = wght_cp * MCLININTPOL(alb_surf, tetasurf, ntetasurf,  ACOS(-(cos_dir_cp(3)))/DTR)
                          ENDIF

                          ! LE on CP
                          IF ( (i_cp .eq. 1) .or. &
                               ((vrm_n_sccp+1-count_scat_cp) .lt. vrm_n_lecp) ) THEN
                             ! for the first CP of the MP, we must perform the LE for each scattering
                             ! to complet the Von-Neumann serie
                             ! for the following CP, we perform the LE for the last vrm_n_lecp scattering only
                             CALL local_estimate(nmat, nelem, surface, nlay, cext_ly, z_ref, tau_tot,  &
                                  ind_cp, nmu, thetad, p_ptcle_ly(ilayer_cp,:,:), delta_ray,           &
                                  wght_cp, Stk_cp, z_phot_cp, cos_dir_cp, v_cp, ilayer_cp,             &  
                                  thetav(igeom), cos_dir_sensor(igeom, :), le)
                             Stoke_pack(ipack, igeom, :) = Stoke_pack(ipack, igeom, :) + le(:)
                             !write(*,*) ' CP    ', count_scat, count_scat_cp, i_cp, le(1)
                          ENDIF

                          IF (count_scat_cp .le. vrm_n_sccp) then
                             ! compute new stoke vector, new cos_dir, u and v after scattering
                             ! NB : CP is always DDISed 
                             CALL scatter(ithread, nmat, nelem, surface,                                              &  
                                  .true., vrm_epsilon_ddis, p_ptcle_ddis, prf_ptcle_ddis,                             &
                                  cos_dir_sensor(igeom, :),                                                           &
                                  ind_cp, nmu, mu, thetad, p_ptcle_ly(ilayer_cp,:,:), prf_ptcle_ly(ilayer_cp,:), delta_ray, &
                                  wght_cp, Stk_cp, v_cp, u_cp, cos_dir_cp)
                          ENDIF

                       ENDIF !condition on CP interacting

                       !if (leave_cp .ne. 0) write(*,*) ' CP leave=', leave_cp

                    ENDDO ! stop following the CP photon

                 ENDDO ! loop on geometries

              ENDIF ! ==== CP

              ! FOR the MP : Compute new stoke vector, new cos_dir, u and v after scattering (or reflection)
              CALL scatter(ithread, nmat, nelem, surface,                                           &  
                   ddis_mp, vrm_epsilon_ddis,  p_ptcle_ddis, prf_ptcle_ddis,                        &
                   cos_dir_sensor(1, :),                                                            &
                   ind, nmu, mu, thetad, p_ptcle_ly(ilayer,:,:), prf_ptcle_ly(ilayer,:), delta_ray, &
                   wght, Stk, v, u, cos_dir)

           ENDIF ! MP interacted in the atmosphere

           !if (leave .ne. 0) write(*,*) 'MP leave=', leave

        ENDDO ! stop following the MP photon

        !***************************************
        !***************************************
        !*********    END      VRM     *********
        !***************************************
        !***************************************

     ENDIF !test on vrm (0 or 1)

  ENDDO photon
  !$OMP END DO
  !$OMP END PARALLEL 

  !--- Normalize the final Stoke vector for intensity
  stoke_pack(:, :, :) = stoke_pack(:, :, :) * costhetas * stoke_in(1)

  ! Compute statistics 
  stoke_out_mean(:,:) = 0.0D0
  do ipack = 1, npack
     do igeom = 1, ngeom
        stoke_pack(ipack, igeom, :) = stoke_pack(ipack, igeom, :) / count_phot_pack(ipack) 
        stoke_out_mean(igeom,:)     = stoke_out_mean(igeom,:) + stoke_pack(ipack, igeom, :)    
     enddo
  enddo
  stoke_out_mean(:,:) = stoke_out_mean(:,:) / npack
  stoke_out(:,:)      = stoke_out_mean(:,:)
  ! evaluate the uncertainty over all the packs
  stoke_out_sig(:,:)  = 0.0D0
  do ipack = 1, npack
     do igeom = 1, ngeom
        stoke_out_sig(igeom,:) = stoke_out_sig(igeom,:) + ( stoke_pack(ipack, igeom, :)**2.0D0 - stoke_out_mean(igeom,:)**2.0D0 )
     enddo
  enddo
  stoke_out_sig(:,:) = SQRT(stoke_out_sig(:,:) / npack ) 
!!$  ! evaluate the error on one pack
!!$  ! just assume the uncertainty goes like 1/SQRT(nphot)
!!$  sig_one(:,:) = stoke_out_sig(:,:) * SQRT(DBLE(npack))
!!$  ! then recompute the uncertainty estimate
!!$  stoke_out_sig(:,:)  = 0.0D0
!!$  do ipack = 1, npack
!!$     do igeom = 1, ngeom
!!$        stoke_out_sig(igeom,:) = stoke_out_sig(igeom,:) + ( (sig_one(igeom,:)**2.0D0+stoke_pack(ipack, igeom, :)**2.0D0) &
!!$             - stoke_out_mean(igeom,:)**2.0D0 )
!!$     enddo
!!$  enddo
!!$  stoke_out_sig(:,:) = SQRT(stoke_out_sig(:,:) / npack ) 


  !--- Write some interesting quantities
  if (verbose .eqv. .true.) then 
     write(*,*) '========== MCRAD1D ============='
     write(*,*) ' # limit of interaction          :', nmax_interact
     write(*,*) '================================'
     write(*,*) ' # of photon                     :',nphot
     if (vrm .ne. 2) then
        WRITE(*,*) ' # of photons leave from the top :',leave_top
        WRITE(*,*) '        -           from bottom  :',leave_bot
        WRITE(*,*) ' # of vanished photons           :',vanished
        write(*,*) ' '
        write(*,*) ' # of interacting photons        :',interact
        write(*,*) ' '
        WRITE(*,*) ' total # of interactions         :',tot_interact
        WRITE(*,*) ' total # of surface reflections  :',tot_reflect
     endif
     write(*,*) '================================'
  endif

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !          FOR TESTING
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! if (ngeom .eq. 1) then
  !    write(*,*) '(MCRAD1D) Write Pack radiances results'
  !    OPEN (UNIT = 10, FILE = '/Users/compiegne/work/artdeco/out/res_mcrad1d_pack.dat', STATUS = 'unknown')
  !    WRITE (UNIT=10, FMT=*) nmat
  !    WRITE (UNIT=10, FMT=*) npack
  !    DO i = 1, npack
  !       SELECT CASE (nmat)
  !       CASE(1) 
  !          WRITE (UNIT=10, FMT='(2x,I8,200(2x,E16.8))') count_phot_pack(i), Stoke_pack(i, 1, 1)
  !       CASE(3) 
  !          WRITE (UNIT=10, FMT='(2x,I8,200(2x,E16.8))') count_phot_pack(i), Stoke_pack(i, 1, 1), &
  !               Stoke_pack(i, 1, 2) , &
  !               Stoke_pack(i, 1, 3)
  !       CASE(4) 
  !          WRITE (UNIT=10, FMT='(2x,I8,200(2x,E16.8))') count_phot_pack(i), Stoke_pack(i, 1, 1), &
  !               Stoke_pack(i, 1, 2) , &
  !               Stoke_pack(i, 1, 3) , &
  !               Stoke_pack(i, 1, 4)
  !       END SELECT
  !    ENDDO
  !    CLOSE(UNIT = 10)
  ! else
  !    write(*,*) '(MCRAD1D) ngeom must be one if you want to get a file with individual pack results'
  ! endif
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !             END TO TEST
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (surface .eq. 'brdf') then
     DEALLOCATE( phisurf, tetasurf, prf_phi_surf, &
          prf_mu_surf,alb_surf, brdf)
  endif

  if (verbose)  write(*,*) ' '

  CALL par_zigunset()
  DEALLOCATE(p_ptcle_ly)
  !DEALLOCATE(p_ptcle_ly_reg)
  deallocate(seed)

END SUBROUTINE mcrad1d

!================================================================================================












