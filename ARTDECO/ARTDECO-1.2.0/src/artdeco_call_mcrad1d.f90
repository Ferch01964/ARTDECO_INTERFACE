

MODULE MCALL_MCRAD1D

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CALL_MCRAD1D

CONTAINS

!--------------------------

  SUBROUTINE CALL_MCRAD1D

    USE MCONSTANTS, only : dp, max_len, undef_dp
    USE MCOMMON, only : verbose, &
         mc_wrr_crit, &
         mc_vrm_n_sccp, &
         mc_vrm_n_lecp, &
         mc_vrm_n_firstcp, &
         mc_vrm, &
         mc_nsplmax, &
         mc_nphot, &
         mc_nmax_interact, &
         mc_epsilon_ddis, &
         mc_wspl_crit, &
         mc_wrr_min, &
         nvaa, &
         nvza, &
         nsza, &
         vaa, &
         vza, &
         sza, &
         nptcle, &
         nmat, &
         nlambda, &
         wlambda, &
         depol, &
         thermal, &
         surface_type, &
         surface_family, &
         surface_albedo, &
         nlayers, &
         ptcle_dtau_layers, &
         ptcle_ssa_layers, &
         ptcle_nang_phasemat_layers, &
         ptcle_phasemat_layers, &
         ptcle_u_phasemat_layers, &
         gas_dtauabs_layers, &
         gas_dtausca_layers, &
         t_atm, &
         alt_atm, &
         SvR, &
         SvR_sig, &
         F0, &
         Svin, &
         rt_cputime, &
         mode, &
         n_aitaui_total, &
         varsol_fact

    IMPLICIT NONE

    INTEGER  :: mc_ngeom 
    REAL(kind=dp) :: mc_thetas
    REAL(kind=dp), allocatable  :: mc_thetav(:)
    REAL(kind=dp), allocatable  :: mc_phir(:)
    INTEGER  :: mc_nlay
    INTEGER  :: mc_nmu_ly
    REAL(kind=dp), allocatable  :: mc_mu_ly(:)
    REAL(kind=dp), allocatable  :: mc_p_ly(:,:,:)
    REAL(kind=dp), allocatable  :: mc_tau_ptcle_ly(:)
    REAL(kind=dp), allocatable  :: mc_ssa_ptcle_ly(:)
    REAL(kind=dp), allocatable  :: mc_tau_gas_ly(:)
    REAL(kind=dp), allocatable  :: mc_ssa_gas_ly(:)
    REAL(kind=dp), allocatable  :: mc_z_ref(:)
    REAL(kind=dp)  :: mc_delta_ray
    LOGICAL  :: mc_thermal
    REAL(kind=dp), allocatable  :: mc_temp(:)
    REAL(kind=dp)  :: mc_Stoke_in(4)
    INTEGER  :: mc_nmat
    CHARACTER(len=max_len) :: mc_surface
    REAL(kind=dp) :: mc_surface_albedo
    REAL(kind=dp), allocatable  :: mc_Stoke_out(:,:)
    REAL(kind=dp), allocatable  :: mc_Stoke_sig(:,:)

    INTEGER :: igeom
    INTEGER :: i, u, k, j

    ! time
    REAL :: tcpu1
    REAL :: tcpu2

    !------------

    IF (verbose) THEN
       WRITE (*,*) ''
       WRITE (*, *)'  (sub. call_mcrad1d) set MCRAD1D arguments'
       WRITE (*,FMT='(A,1x,I6,1x,A)') '                  MCRAD1D will be called ', n_aitaui_total, ' times'
    ENDIF

    ! IF ((nmat.gt.1).and.(surface_family .eq. 'brdf')) then
    !    print*, 'ERROR : BDRF not working for nmat>1 in MCRAD1D'
    ! endif

    ! ARGUMENTS that does not depend on lambda
    mc_ngeom = nvza * nvaa
    mc_nmat  = nmat
    mc_nlay  = nlayers
    ALLOCATE(mc_thetav(mc_ngeom), &
         mc_phir(mc_ngeom), &
         mc_z_ref(0:mc_nlay), &
         mc_temp(0:mc_nlay), &
         mc_tau_ptcle_ly(mc_nlay), &
         mc_ssa_ptcle_ly(mc_nlay), &
         mc_tau_gas_ly(mc_nlay), &
         mc_ssa_gas_ly(mc_nlay), &
         mc_Stoke_out(mc_ngeom, mc_nmat), &
         mc_Stoke_sig(mc_ngeom, mc_nmat))
    !------------
    ! get the geometry variable
    igeom = 0
    DO u = 1, nvza
       DO k = 1, nvaa
          igeom = igeom + 1
          mc_thetav(igeom) = vza(u)
          mc_phir(igeom)   = vaa(k) 
       ENDDO
    ENDDO

    mc_thermal  = thermal
    mc_surface  = surface_family

    ! layer must be in ascending order for mcrad1d
    do i = 0, mc_nlay
       mc_z_ref(i) = alt_atm(mc_nlay+1-i)
       mc_temp(i)  = t_atm(mc_nlay+1-i) 
    enddo

    IF (verbose) THEN
       WRITE (*,*) ''
    ENDIF

    ! init in case no lambertian surface
    mc_surface_albedo = 0.0_dp

    DO j = 1, nlambda

       IF (verbose) WRITE (*, FMT='(1x,A45,2x,F12.5,1x,A8)')  &
            &'  (sub. call_mcrad1d) run MCRAD1D for wlambda=', wlambda(j),' microns'  

       if (surface_family .eq. 'lambert') mc_surface_albedo = surface_albedo(j)
       mc_nmu_ly = ptcle_nang_phasemat_layers(j)

       ALLOCATE(mc_mu_ly(mc_nmu_ly), &
            mc_p_ly(mc_nlay,6,mc_nmu_ly))

       mc_mu_ly(:)   = ptcle_u_phasemat_layers(1:mc_nmu_ly,j)

       ! layer must be in ascending order for mcrad1d
       if (mode.ne.'kdis') then
          do i = 1, mc_nlay
             ! no gas absorption in MC for now
             mc_tau_gas_ly(i)   = gas_dtausca_layers(mc_nlay+1-i,j) 
             mc_ssa_gas_ly(i)   = 1.0
          enddo
       else
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' (sub. call_mcrad1d)  ERROR        '
          WRITE(*,*) ' K-dis is not implemented for MCRAD1D '
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          STOP
       endif

       if (nptcle .gt. 0) then
          do i = 1, mc_nlay
             mc_p_ly(i,1,:)     = ptcle_phasemat_layers(mc_nlay+1-i,1,1:mc_nmu_ly,j) ! P11
             mc_p_ly(i,2,:)     = ptcle_phasemat_layers(mc_nlay+1-i,2,1:mc_nmu_ly,j) ! P22
             mc_p_ly(i,3,:)     = ptcle_phasemat_layers(mc_nlay+1-i,5,1:mc_nmu_ly,j) ! P12 = P21
             mc_p_ly(i,4,:)     = ptcle_phasemat_layers(mc_nlay+1-i,3,1:mc_nmu_ly,j) ! P33
             mc_p_ly(i,5,:)     = ptcle_phasemat_layers(mc_nlay+1-i,4,1:mc_nmu_ly,j) ! P44
             mc_p_ly(i,6,:)     = ptcle_phasemat_layers(mc_nlay+1-i,6,1:mc_nmu_ly,j) ! P34 = -P43
             mc_tau_ptcle_ly(i) = ptcle_dtau_layers(mc_nlay+1-i,j)
             mc_ssa_ptcle_ly(i) = ptcle_ssa_layers(mc_nlay+1-i,j)
          enddo
       else
          mc_p_ly(:,:,:)     = 0.0_dp
          mc_tau_ptcle_ly(:) = 0.0_dp
          mc_ssa_ptcle_ly(:) = 0.0_dp
       endif

       mc_delta_ray   = depol(j)

       mc_Stoke_out(:,:)   = undef_dp 
       mc_Stoke_sig(:,:)   = undef_dp 

       CALL CPU_TIME(tcpu1)

       DO i = 1, nsza

          mc_Stoke_in(:) = Svin(:) * F0(j) * varsol_fact(i) 
          mc_thetas = sza(i)
          call mcrad1d(verbose, wlambda(j), mc_ngeom, mc_thetas, mc_thetav, mc_phir,   &
               mc_nlay, mc_nmu_ly, mc_mu_ly, mc_p_ly,                              &
               mc_tau_ptcle_ly,                                                &
               mc_ssa_ptcle_ly,                                                &
               mc_tau_gas_ly,                                                  &
               mc_ssa_gas_ly,                                                  &
               mc_z_ref,                                                       & 
               mc_delta_ray,                                                   &
               mc_thermal, mc_temp,                                            &
               mc_Stoke_in,                                                    &
               mc_nphot,                                                       &
               mc_surface,                                                     &
               mc_surface_albedo,                                              &
               mc_nmat, mc_nmax_interact,                                      &
               mc_vrm, mc_epsilon_ddis,                                        &
               mc_vrm_n_firstcp, mc_vrm_n_sccp, mc_vrm_n_lecp,                 &
               mc_wspl_crit,  &
               mc_nsplmax,    &
               mc_wrr_crit,   &
               mc_wrr_min,    &
               mc_Stoke_out, mc_Stoke_sig)
          igeom = 0
          ! implement the common variable for stokes vector for reflected intensities for each geometry
          DO u = 1, nvza
             DO k = 1, nvaa
                igeom = igeom + 1
                SvR(i,u,k,1,j)     = mc_Stoke_out(igeom,1)
                SvR_sig(i,u,k,1,j) = mc_Stoke_sig(igeom,1)
                if (nmat .gt. 1)then
                   SvR(i,u,k,2,j)     = mc_Stoke_out(igeom,2) 
                   SvR_sig(i,u,k,2,j) = mc_Stoke_sig(igeom,2) 
                   SvR(i,u,k,3,j)     = -mc_Stoke_out(igeom,3) 
                   SvR_sig(i,u,k,3,j) = mc_Stoke_sig(igeom,3) 
                endif
                if (nmat .gt. 3) then 
                   SvR(i,u,k,4,j)     = -mc_Stoke_out(igeom,4) 
                   SvR_sig(i,u,k,4,j) = mc_Stoke_sig(igeom,4) 
                endif
             ENDDO
          ENDDO
       ENDDO ! on sza

       CALL CPU_TIME(tcpu2)
       rt_cputime(j) = tcpu2-tcpu1
       ! IF (verbose) THEN
       !    WRITE (*,FMT='(1x,A35,1PE14.5, A9)') '  (sub. call_mcrad1d) CPU time     ', tcpu2-tcpu1, ' sec(s) -'
       ! ENDIF

       DEALLOCATE(mc_mu_ly, mc_p_ly)

    ENDDO ! on lambda

    DEALLOCATE(mc_thetav, &
         mc_phir, &
         mc_z_ref, &
         mc_temp, &
         mc_tau_ptcle_ly, &
         mc_ssa_ptcle_ly, &
         mc_tau_gas_ly, &
         mc_ssa_gas_ly, &
         mc_Stoke_out, &
         mc_Stoke_sig)

    IF (verbose) THEN
       WRITE (*,*)                          '  '
       WRITE (*,*)                          '  (sub. call_mcrad1d) finished to run MCRAD1D'
       WRITE (*,*)                          '  '
    ENDIF

  END SUBROUTINE CALL_MCRAD1D

END MODULE MCALL_MCRAD1D
