MODULE mmcrad1d_submain

  IMPLICIT NONE

  PUBLIC :: target_choice, new_position, get_new_photon, get_prfs, INIT_RANDOM_SEED, &
       scatter, local_estimate, mc_xinteg2, move_photon,  get_wspl, &
       split_photon

  PRIVATE :: print_phot, newdirection, rot_euler, stoke_rot, stoke, stoke2, Compute_theta, &
       compute_psi, compute_trans, get_pray, getp11, &
       surface_le, surface_scatter, surface_stoke

CONTAINS

  FUNCTION MC_XINTEG2(imin, imax, n, xin, yin)
    ! computes integral of yin, variable step 
    ! make sure you have at least two points of integration
    ! Trapeze method

    IMPLICIT NONE

    REAL (8) :: MC_XINTEG2

    INTEGER, INTENT(IN)  :: imin
    INTEGER, INTENT(IN)  :: imax
    INTEGER, INTENT(IN)  :: n
    REAL (8), INTENT(IN) :: xin(n)
    REAL (8), INTENT(IN) :: yin(n)
    REAL (8)             :: xa, ya, xb, yb
    REAL (8)             :: primitaux
    INTEGER              :: i

    xa = xin(imin)
    ya = yin(imin)
    primitaux = 0.0D0
    DO i=imin+1,imax
       xb = xin(i)
       yb = yin(i)
       primitaux = primitaux + (xb - xa) * (ya + yb)
       xa = xb
       ya = yb
    ENDDO
    MC_XINTEG2 = primitaux * 0.5D0

  END FUNCTION MC_XINTEG2

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine print_phot(iphot, nmat, stoke_0, wght, cos_dir, u, v, z_phot, ilayer, &
       count_scat, leave, count_reflect)

    implicit none
    integer, intent(in) :: iphot
    integer, intent(in) :: nmat
    real(8), intent(in) :: stoke_0(nmat)
    real(8), intent(in) :: wght
    real(8), intent(in) :: cos_dir(3)
    real(8), intent(in) :: u(3)
    real(8), intent(in) :: v(3)
    real(8), intent(in) :: z_phot
    integer, intent(in) :: ilayer
    integer, intent(in) :: count_scat
    integer, intent(in) :: leave
    integer, intent(in) :: count_reflect

    integer :: j

    !-------

    write(*,*) '==================='
    write(*,*) 'phot #', iphot
    write(*,*) 'stoke_0       =', (stoke_0(j), j=1,nmat)
    write(*,*) 'wght          =', wght
    write(*,*) 'cos_dir       =', (cos_dir(j), j=1,3)
    write(*,*) 'u             =', (u(j), j=1,3)
    write(*,*) 'v             =', (v(j), j=1,3)
    write(*,*) 'z_phot        =', z_phot
    write(*,*) 'ilayer        =', ilayer
    write(*,*) 'count_scat    = ', count_scat
    write(*,*) 'count_reflect = ', count_reflect
    write(*,*) 'leave         = ', leave
    write(*,*) '==================='

  end subroutine print_phot

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE INIT_RANDOM_SEED()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE get_new_photon(thermal, sinthetas, costhetas, Stoke_in, nmat, nlay, z_ref,&
       stoke_0, wght, cos_dir, u, v, z_phot, ilayer, d, leave, count_scat, &
       count_reflect)

    USE MMCRAD1D_CONSTANTS
    implicit none

    logical, intent(in) :: thermal
    real(8), intent(in) :: costhetas 
    real(8), intent(in) :: sinthetas 
    real(8), intent(in) :: stoke_in(4)
    integer, intent(in) :: nmat
    integer, intent(in) :: nlay
    real(8), intent(in) :: z_ref(0:nlay)

    real(8), intent(out) :: stoke_0(nmat)
    real(8), intent(out) :: wght
    real(8), intent(out) :: cos_dir(3)
    real(8), intent(out) :: u(3)
    real(8), intent(out) :: v(3)
    real(8), intent(out) :: z_phot
    integer, intent(out) :: ilayer
    real(8), intent(out) :: d
    integer, intent(out) :: leave
    integer, intent(out) :: count_scat
    integer, intent(out) :: count_reflect

    !--------------------

    IF (thermal .eqv. .false.) then
       !***********************************
       !    Initialization of Stoke vectors
       !    we work with the normalized stoke vector
       IF (nmat .EQ. 1) THEN
          Stoke_0(1)      = Stoke_in(1) / Stoke_in(1)
       ELSE
          Stoke_0(1:nmat) = Stoke_in(1:nmat) / Stoke_in(1)
       ENDIF
       wght = 1.0D0

       ! Initialization of direction cosines
       ! We take the photon to be incident in the Oxz plane
       cos_dir(:) = (/ sinthetas, 0.0D0, -costhetas /) 
       !--- initialization of the photon polarization reference frame
       !    u is taken along the photon direction = cos_dir
       !    v is taken to be in the Oy direction (perp. to W since initial v is in the Oxz plane)
       !    w would result from the crossed produce of w and v
       u(:) = cos_dir(:)
       v(1) = 0.0D0
       v(2) = 1.0D0
       v(3) = 0.0D0

       !--- Initialization of the vector position
       z_phot = z_ref(nlay)

       ilayer = nlay    ! start at the first layer on top of the atmosphere
       d      = 0.0D0   ! distance of interaction in the cloud in Km

       leave          = 0
       count_scat     = 0 
       count_reflect  = 0 

    endif

  END SUBROUTINE get_new_photon

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE get_pray(theta, nelem, delta, pray)

    USE MMCRAD1D_CONSTANTS

    implicit none

    !--- Input Arguments
    real(8), intent(in)  :: theta
    integer, intent(in)  :: nelem
    real(8), intent(in)  :: delta
    real(8), intent(out) :: pray(nelem)

    !--- Local Arguments
    real(8) :: Gdelta, Gdeltap
    real(8) :: mu

    !******************************************
    !     Compute Rayleigh phase matrix
    !******************************************

    if (delta .ne. 0D0) then

       GDelta  = (1.0D0 - delta) / (1.0D0 + delta / 2.0D0)
       GDeltap = (1.0D0 - 2.0D0*delta) / (1.0D0 - delta)

       if (nelem .ge.1) then
          Pray(1) = GDelta * 3.0D0 * (1.0D0 + (COS(DTR*theta))**2.D0) / 4.0D0 &
               + (1.0D0 - GDelta)                                          !P11
       endif
       if (nelem .ge. 4) then    
          Pray(2) = Pray(1) - (1.0D0 - GDelta)                        !P22
          Pray(3) = - GDelta * 3.0D0 * (SIN(DTR*theta))**2.D0 /4.0D0  !P21=P12
          Pray(4) = GDelta * 3.0D0 * COS(DTR*theta) / 2.0D0           !P33
       endif
       if (nelem .eq. 6) then
          Pray(5) = GDeltap * Pray(4)                                  !P44
          Pray(6) = 0.0D0                                              !P34=-P43
       endif

    else

       mu = COS(DTR*theta)

       if (nelem .ge.1) then
          Pray(1) =  1.0D0 + mu**2.D0    !P11
       endif
       if (nelem .ge. 4) then    
          Pray(2) =  1.0D0 + mu**2.D0    !P22
          Pray(3) =  mu**2.0D0 - 1.0D0   !P21=P12
          Pray(4) =  2.0D0*mu            !P33
       endif
       if (nelem .eq. 6) then
          Pray(5) =  2.0D0*mu  !P44
          Pray(6) =  0.0       !P34=-P43
       endif

       Pray(:) = 3.0D0 / 4.0D0 * Pray(:)

    endif


    RETURN 

  END SUBROUTINE get_pray

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE target_choice(ithread,wgas,ind)

    USE MMCRAD1D_ZIGGURAT

    implicit none

    !--- Input Arguments
    real(8), intent(in)    :: wgas ! Cext(gas) / C_ext(gas+particles)
    integer, intent(in)    :: ithread
    integer, intent(out)   :: ind
    !--- Local Arguments
    real(8)           :: xr

    !--- get a random number
    !CALL RANDOM_NUMBER(xr)
    xr=par_ziguni(ithread)  

    IF (xr .LE. wgas) THEN
       ind = 0   ! target = gas
    ELSE
       ind = 1   ! target = ptcle 
    ENDIF

    RETURN

  END SUBROUTINE target_choice

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE New_position(old_z_phot, costheta, d, z_phot)

    implicit none

    !--- Input Arguments
    real(8), intent(in)               :: old_z_phot
    real(8), intent(in)               :: costheta
    real(8), intent(in)               :: d
    !--- Output Arguments
    real(8), intent(out)              :: z_phot

    z_phot = old_z_phot + d * costheta

    RETURN

  END SUBROUTINE New_position

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE newdirection(ithread, ind, nmu, mu, prf_ptcle, ctheta, theta, phi)

    use MMCRAD1D_CONSTANTS, only: DTR, Pi, EPS
    use MMCRAD1D_ZIGGURAT

    implicit none

    !--- Input Arguments
    integer, intent(in) :: ithread
    integer, intent(in) :: ind
    integer, intent(in) :: nmu
    real(8), intent(in) :: mu(nmu)
    real(8), intent(in) :: prf_ptcle(nmu)

    !--- Output Arguments
    real(8), intent(out)                     :: ctheta
    real(8), intent(out)                     :: theta
    real(8), intent(out)                     :: phi

    !--- Local Arguments
    integer                                 :: N1, N2, J0
    real(8)                                 :: q, expo, u ,v
    real(8)                                 :: xt, b, a, arg
    real(8)                                 :: rpd

    !-------------------

    IF (ind .EQ. 0) THEN  ! Rayleigh compute exact theta

       ! get random number
       !CALL RANDOM_NUMBER(xt)
       xt = par_ziguni(ithread)
       q = 4.0D0*(2.0D0*xt-1.0D0)
       expo = 1.0D0/3.0D0
       u = -q/2.0D0 + DSQRT(q*q/4.0D0+1.0D0)
       v = -q/2.0D0 - DSQRT(q*q/4.0D0+1.0D0)
       IF (v .LT. 0.0D0) THEN
          ctheta = u**expo - DABS(v)**expo
       ELSE
          ctheta = u**expo + v**expo
       ENDIF
       IF (ctheta .LT. -1.0D0) ctheta = -1.0D0
       IF (ctheta .GT.  1.0D0) ctheta =  1.0D0
       theta = ACOS(ctheta)/DTR

    ELSE   !Compute theta from a dichotomic method...

       !*************************************************************************
       !    Evaluate the scattering angle theta from PRF and a random number
       !*************************************************************************
       ! get random number
       !CALL RANDOM_NUMBER(xt)
       xt = par_ziguni(ithread)
       !--- Use a dichotomic method to find the 2 indice (J0+1 and J0) that 
       !    surround xt
       N1 = 1
       N2 = nmu
       DO WHILE ((N2-N1) .GT. 1)
          J0 = FLOOR(REAL((N1 + N2) / 2))
          IF (xt .GE. prf_ptcle(J0)) THEN
             N1 = J0
          ELSE
             N2 = J0
          ENDIF
       ENDDO
       IF ((xt .LT. prf_ptcle(J0)) .AND. (J0 .GT. 1)) J0 = J0 - 1
       !--- Now compute theta
       a = (prf_ptcle(J0+1) - prf_ptcle(J0)) / (mu(J0+1) - mu(J0)) !slope
       IF (DABS(a) .LE. EPS) THEN
          theta = ACOS(mu(J0)) / DTR
          ctheta = mu(J0)
       ELSE
          b = prf_ptcle(J0) - a * mu(J0)     ! ordonnee a l'origine
          arg = (xt - b) / a
          IF (arg .GT. 1.0D0) arg = 1.0D0
          IF (arg .LT. -1.0D0) arg = -1.0D0
          theta = ACOS(arg) / DTR
          ctheta = arg
       ENDIF

    ENDIF

    IF (theta .GT. 180.0D0) theta = 180.0D0
    IF (theta .LT.   0.0D0) theta =   0.0D0

    ! between 0 and 2pi
    ! get random number
    !CALL RANDOM_NUMBER(rpd)
    rpd = par_ziguni(ithread)
    phi = 2.0D0*Pi*rpd
    !--- degree
    phi = phi / DTR

    RETURN

  END SUBROUTINE newdirection

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE rot_euler(vecti,k,angle,vectd)

    use MMCRAD1D_CONSTANTS, only: DTR, Pi

    implicit none

    !--- Input Arguments
    real(8), intent(in)                 :: angle
    real(8), intent(in), dimension(3)   :: vecti, k

    !--- Local Arguments
    real(8)                                 :: cs, sn, nu
    real(8)                                 :: ui, vi, wi
    real(8)                                 :: kx, ky, kz
    real(8)                                 :: den

    !--- Output Arguments
    real(8), intent(out), dimension(3)     :: vectd

    !**************************************************************************
    !    cs=Cos(angle), sn=sin(angle) and nu = 1 - cos(angle)
    !**************************************************************************
    cs = COS(DTR*angle)
    sn = SIN(DTR*angle)
    nu = 1.0D0 - cs

    !**************************************************************************
    !   Expression of the composante of the initial vector and k
    !**************************************************************************
    ui = vecti(1)
    vi = vecti(2)
    wi = vecti(3)
    kx = k(1)
    ky = k(2)
    kz = k(3)

    !**************************************************************************
    !   From Ramella-Roman et al., 2005a, compute the new coordinate after
    !   rotation of angle around k
    !**************************************************************************
    vectd(1) = (kx*kx*nu+cs)    * ui + (ky*kx*nu-kz*sn) * vi + (kz*kx*nu+ky*sn) * wi
    vectd(2) = (kx*ky*nu+kz*sn) * ui + (ky*ky*nu+cs)    * vi + (ky*kz*nu-kx*sn) * wi
    vectd(3) = (kx*kz*nu-ky*sn) * ui + (ky*kz*nu+kx*sn) * vi + (kz*kz*nu+cs)    * wi
    den = DSQRT(vectd(1)**2+vectd(2)**2+vectd(3)**2)
    vectd(:) = vectd(:) / den

    RETURN

  END SUBROUTINE rot_euler

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE Stoke_rot(nmat,cs,sn,Stoke_old,Stoke)

    implicit none

    !--- Input Args.
    integer, intent(in)                  :: nmat
    real(8), intent(in)                  :: cs, sn
    real(8), intent(in), dimension(nmat) :: Stoke_old

    !--- Output Args.
    real(8), intent(out), dimension(nmat) :: Stoke

    !-------

    IF (nmat .EQ. 3) THEN
       Stoke(1) =  Stoke_old(1)
       Stoke(2) =  cs*Stoke_old(2) + sn*Stoke_old(3)
       Stoke(3) = -sn*Stoke_old(2) + cs*Stoke_old(3)
    ELSE IF (nmat .EQ. 4) THEN
       Stoke(1) =  Stoke_old(1)
       Stoke(2) =  cs*Stoke_old(2) + sn*Stoke_old(3)
       Stoke(3) = -sn*Stoke_old(2) + cs*Stoke_old(3)
       Stoke(4) =  Stoke_old(4)
    ELSE
       WRITE(6,*)'Pb, normalement  pas de rotation car nmat =',nmat
    ENDIF

    RETURN

  END SUBROUTINE Stoke_rot

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE Stoke(ind,nmat,nelem,nmu,thetad,P,Stoke_0,theta,delta,Stoke_1,p11)

    use MMCRAD1D_CONSTANTS, only: DTR, flag_interpol, Pi, EPS1
    use MMCRAD1D_INTERPOLATE

    implicit none

    !--- Input Arguments
    integer,intent(in)                           :: ind
    integer,intent(in)                           :: nmat
    integer,intent(in)                           :: nelem
    integer,intent(in)                           :: nmu
    real(8),intent(in), dimension(nmu)           :: thetad
    real(8),intent(in), dimension(nelem,nmu)     :: P
    real(8),intent(in), dimension(nmat)          :: Stoke_0
    real(8),intent(in)                           :: theta
    real(8),intent(in)                           :: delta
    real(8),intent(out), dimension(nmat)         :: Stoke_1
    real(8),intent(out)                          :: p11

    !--- Local Arguments
    integer                                 :: i
    integer                                 :: N1, N2, J0
    integer, dimension(1)                   :: ival
    integer, dimension(nmu)                 :: mask

    real(8)                                 :: P1
    real(8)                                 :: deriv

    real(8), dimension(nmu)                 :: Pii
    real(8), dimension(nmat,nmat)           :: Mat_P
    real(8), dimension(nelem)               :: coef_P

    !===================

    ! <MC 06/02/2012> 
    IF (nmat .eq. 1 ) THEN 
       Stoke_1 = Stoke_0
       p11     = -32768.0D0
       RETURN
    ENDIF
    ! and suppress any IF (nmat .eq. 1 )
    ! in the rest of the routine
    ! <MC 06/02/2012> 

    IF (ind .eq. 1) THEN ! PARTICLE

       !**************************************************************************
       !        Look for existing values (don't need to interpolate
       !**************************************************************************
       mask(:) = 0
       WHERE ((thetad(:) .GT. (theta - EPS1)) .AND.         &
            (thetad(:) .LT. (theta + EPS1)))
          mask(:) = 1
       END WHERE

       !**************************************************************************
       !    Interpolation of the coefficient of the phase matrix in theta
       !**************************************************************************
       IF (sum(mask(:)).GT. 0) THEN ! We don't need to interpolate

          ival = MAXLOC(mask(:))

          DO i = 1, nelem
             IF (i .EQ. 1) THEN
                coef_P(i) = LOG10(P(i,ival(1)))
             ELSE
                coef_P(i) = P(i,ival(1))
             ENDIF
          ENDDO !i

       ELSE  ! WE NEED TO INTERPOLATE

          DO i = 1,nelem
             IF (i .EQ. 1) THEN
                Pii(:) = LOG10(P(i,:))
             ELSE
                Pii(:) = P(i,:)
             ENDIF
             !**************************** ! Linear interpolation ! ************************
             IF (flag_interpol .EQ. 0) THEN  ! linear interpolation
                IF (i .EQ. 1) THEN ! Compute J0 only once!!!!
                   !--- Use a dichotomic method to find the 2 indice (J0+1 and J0) that 
                   !    surround theta
                   N1 = 1
                   N2 = nmu
                   DO WHILE ((N2-N1) .GT. 1)
                      J0 = FLOOR(REAL(N1 + N2) / 2.0D0)
                      IF (theta .GE. thetad(J0)) THEN
                         N1 = J0
                      ELSE
                         N2 = J0
                      ENDIF
                   ENDDO
                   IF ((theta .LT. thetad(J0)) .AND. (J0 .GT. 1)) J0 = J0 - 1
                ENDIF
                !--- Linear interpolation of the phase matrix on mu
                coef_P(i) = Pii(J0) + (Pii(J0+1) - Pii(J0)) * (theta - thetad(J0))&
                     / (thetad(J0+1) - thetad(J0))
                !************************* ! Polynomial interpolation !************************
             ELSE
                CALL SPLINE2(nmu,thetad,Pii,theta,P1,deriv)
                coef_P(i) = P1
             ENDIF
          ENDDO !i

       ENDIF
       p11 = 10.0D0**coef_P(1)

    ELSE !RAYLEIGH PARTICLE

       !*************************************************************************
       !             Compute exact rayleigth coefficients
       !*************************************************************************
       call get_pray(theta, nelem, delta, coef_P)
       p11 = coef_P(1)

    ENDIF

    !--- Build the phase matrix
    !    Remember that: elem = 1 --> P11
    !                          2 --> P22
    !                          3 --> P12 = P21
    !                          4 --> P33
    !                          5 --> P44
    !                          6 --> P34 = -P43
    Mat_P(:,:) = 0.0D0
    IF (ind .EQ. 1) THEN
       IF (nmat .EQ. 3) THEN
          Mat_P(1,1) =   1.0D0
          Mat_P(1,2) =   coef_P(3)
          Mat_P(2,1) =   coef_P(3)
          Mat_P(2,2) =   coef_P(2)
          Mat_P(3,3) =   coef_P(4)
       ELSE
          Mat_P(1,1) =   1.0D0
          Mat_P(1,2) =   coef_P(3)
          Mat_P(2,1) =   coef_P(3)
          Mat_P(2,2) =   coef_P(2)
          Mat_P(3,3) =   coef_P(4)
          Mat_P(3,4) =   coef_P(6)
          Mat_P(4,3) = - coef_P(6)
          Mat_P(4,4) =   coef_P(5)
       ENDIF
    ELSE
       IF (nmat .EQ. 3) THEN
          Mat_P(1,1) =   1.0D0
          Mat_P(1,2) =   coef_P(3) / coef_P(1)
          Mat_P(2,1) =   coef_P(3) / coef_P(1)
          Mat_P(2,2) =   coef_P(2) / coef_P(1)
          Mat_P(3,3) =   coef_P(4) / coef_P(1)
       ELSE
          Mat_P(1,1) =   1.0D0
          Mat_P(1,2) =   coef_P(3) / coef_P(1)
          Mat_P(2,1) =   coef_P(3) / coef_P(1)
          Mat_P(2,2) =   coef_P(2) / coef_P(1)
          Mat_P(3,3) =   coef_P(4) / coef_P(1)
          Mat_P(3,4) =   coef_P(6) / coef_P(1)
          Mat_P(4,3) = - coef_P(6) / coef_P(1)
          Mat_P(4,4) =   coef_P(5) / coef_P(1)
       ENDIF
    ENDIF

    !*************************************************************************
    !       Now compute the new Stoke vector
    !*************************************************************************
    IF (nmat .EQ. 3) THEN
       Stoke_1(1) = Mat_P(1,1)*Stoke_0(1) + Mat_P(1,2)*Stoke_0(2)
       Stoke_1(2) = Mat_P(2,1)*Stoke_0(1) + Mat_P(2,2)*Stoke_0(2)
       Stoke_1(3) = Mat_P(3,3)*Stoke_0(3)
    ELSE
       Stoke_1(1) = Mat_P(1,1)*Stoke_0(1) + Mat_P(1,2)*Stoke_0(2)
       Stoke_1(2) = Mat_P(2,1)*Stoke_0(1) + Mat_P(2,2)*Stoke_0(2)
       Stoke_1(3) = Mat_P(3,3)*Stoke_0(3) + Mat_P(3,4)*Stoke_0(4)
       Stoke_1(4) = Mat_P(4,3)*Stoke_0(3) + Mat_P(4,4)*Stoke_0(4)
    ENDIF

    ! Renormalization of the Stoke Vector so that I = 1.0
    ! Stoke_1(:) = Stoke_1(:) * Stoke_0(1) / Stoke_1(1) 

    RETURN

  END SUBROUTINE Stoke

  ! ========================================================================

  SUBROUTINE Stoke2(nmat,ind,nelem,nmu,thetad,P,Stoke_0,theta,delta,Stoke_1)

    use MMCRAD1D_CONSTANTS, only: DTR, flag_interpol, EPS1
    use MMCRAD1D_INTERPOLATE

    implicit none

    !--- Input Arguments
    integer,intent(in)                            :: nmat
    integer,intent(in)                            :: ind
    integer,intent(in)                            :: nelem
    integer,intent(in)                            :: nmu
    real(8),intent(in), dimension(nmu)            :: thetad
    real(8),intent(in), dimension(nelem,nmu)      :: P
    real(8),intent(in), dimension(nmat)           :: Stoke_0
    real(8),intent(in)                            :: theta
    real(8),intent(in)                            :: delta
    real(8),intent(out), dimension(nmat)          :: Stoke_1

    !--- Local Arguments
    integer                                 :: i 
    integer                                 :: N1, N2, J0
    integer, dimension(1)                   :: ival
    integer, dimension(nmu)                 :: mask

    real(8)                                 :: P1
    real(8)                                 :: deriv

    real(8), dimension(nmu)                :: Pii
    real(8), dimension(nmat,nmat)          :: Mat_P
    real(8), dimension(nelem)              :: coef_P
    real(8)                                :: coef

    !---------------

    IF (ind .EQ. 1) THEN   ! PARTICLE

       !**************************************************************************
       !        Look for existing values (don't need to interpolate
       !**************************************************************************
       mask(:) = 0
       WHERE ((thetad(:) .GT. (theta - EPS1)) .AND.       &
            (thetad(:) .LT. (theta + EPS1)))
          mask(:) = 1
       END WHERE

       !**************************************************************************
       !    Interpolation of the coefficient of the phase matrix in theta
       !**************************************************************************
       IF (sum(mask(:)).GT. 0) THEN ! We don't need to interpolate

          ival = MAXLOC(mask(:))
          DO i = 1,nelem
             IF (i .EQ. 1) THEN
                coef_P(i) = LOG10(P(i,ival(1)))
             ELSE
                coef_P(i) = P(i,ival(1))
             ENDIF
          ENDDO !i

       ELSE    !WE NEED TO INTERPOLATE

          DO i = 1,nelem
             IF (i .EQ. 1) THEN
                Pii(:) = LOG10(P(i,:))
             ELSE
                Pii(:) = P(i,:)
             ENDIF
             !**************************** ! Linear interpolation ! ************************
             IF (flag_interpol .EQ. 0) THEN  ! linear interpolation
                !--- Linear interpolation of the phase matrix on mu
                IF (i .EQ. 1) THEN !Compute J0 only once!!!!
                   !--- Use a dichotomic method to find the 2 indice (J0+1 and J0) that 
                   !    surround theta
                   N1 = 1
                   N2 = nmu
                   DO WHILE ((N2-N1) .GT. 1)
                      J0 = FLOOR(REAL(N1 + N2) / 2.0D0)
                      IF (theta .GE. thetad(J0)) THEN
                         N1 = J0
                      ELSE
                         N2 = J0
                      ENDIF
                   ENDDO
                   IF ((theta .LT. thetad(J0)) .AND. (J0 .GT. 1)) J0 = J0 - 1
                ENDIF
                coef_P(i) = Pii(J0) + (Pii(J0+1) - Pii(J0)) * (theta - thetad(J0))&
                     / (thetad(J0+1) - thetad(J0))
                !************************* ! Polynomial interpolation !************************
             ELSE
                CALL SPLINE2(nmu,thetad,Pii,theta,P1,deriv)
                coef_P(i) = P1
             ENDIF
          ENDDO !i

       ENDIF

    ELSE  !RAYLEIGH PARTICLE
       !*************************************************************************
       !             Compute exact rayleigth coefficients
       !*************************************************************************
       call get_pray(theta, nelem, delta, coef_P)

    ENDIF

    !--- Build the phase matrix
    !    Remember that: elem = 1 --> P11
    !                          2 --> P22
    !                          3 --> P12 = P21
    !                          4 --> P33
    !                          5 --> P44
    !                          6 --> P34 = -P43
    Mat_P(:,:) = 0.0D0
    IF (ind .EQ. 1) THEN
       coef = 10.D0**coef_P(1)
       IF (nmat .EQ. 1) THEN
          Mat_P(1,1) =   coef
       ELSE IF (nmat .EQ. 3) THEN
          Mat_P(1,1) =   coef
          Mat_P(1,2) =   coef_P(3) * Mat_P(1,1)
          Mat_P(2,1) =   coef_P(3) * Mat_P(1,1)
          Mat_P(2,2) =   coef_P(2) * Mat_P(1,1)
          Mat_P(3,3) =   coef_P(4) * Mat_P(1,1)
       ELSE
          Mat_P(1,1) =   coef
          Mat_P(1,2) =   coef_P(3) * Mat_P(1,1)
          Mat_P(2,1) =   coef_P(3) * Mat_P(1,1)
          Mat_P(2,2) =   coef_P(2) * Mat_P(1,1)
          Mat_P(3,3) =   coef_P(4) * Mat_P(1,1)
          Mat_P(3,4) =   coef_P(6) * Mat_P(1,1)
          Mat_P(4,3) = - coef_P(6) * Mat_P(1,1)
          Mat_P(4,4) =   coef_P(5) * Mat_P(1,1)
       ENDIF
    ELSE
       coef = coef_P(1)
       IF (nmat .EQ. 1) THEN
          Mat_P(1,1) =   coef
       ELSE IF (nmat .EQ. 3) THEN
          Mat_P(1,1) =   coef
          Mat_P(1,2) =   coef_P(3)
          Mat_P(2,1) =   coef_P(3)
          Mat_P(2,2) =   coef_P(2)
          Mat_P(3,3) =   coef_P(4)
       ELSE
          Mat_P(1,1) =   coef
          Mat_P(1,2) =   coef_P(3)
          Mat_P(2,1) =   coef_P(3)
          Mat_P(2,2) =   coef_P(2)
          Mat_P(3,3) =   coef_P(4)
          Mat_P(3,4) =   coef_P(6)
          Mat_P(4,3) = - coef_P(6)
          Mat_P(4,4) =   coef_P(5)
       ENDIF
    ENDIF

    !*************************************************************************
    !       Now compute the new Stoke vector
    !*************************************************************************
    IF (nmat .EQ. 1) THEN
       Stoke_1(1) = Mat_P(1,1)*Stoke_0(1)
    ELSE IF (nmat .EQ. 3) THEN
       Stoke_1(1) = Mat_P(1,1)*Stoke_0(1) + Mat_P(1,2)*Stoke_0(2)
       Stoke_1(2) = Mat_P(2,1)*Stoke_0(1) + Mat_P(2,2)*Stoke_0(2)
       Stoke_1(3) = Mat_P(3,3)*Stoke_0(3)
    ELSE
       Stoke_1(1) = Mat_P(1,1)*Stoke_0(1) + Mat_P(1,2)*Stoke_0(2)
       Stoke_1(2) = Mat_P(2,1)*Stoke_0(1) + Mat_P(2,2)*Stoke_0(2)
       Stoke_1(3) = Mat_P(3,3)*Stoke_0(3) + Mat_P(3,4)*Stoke_0(4)
       Stoke_1(4) = Mat_P(4,3)*Stoke_0(3) + Mat_P(4,4)*Stoke_0(4)
    ENDIF

    RETURN

  END SUBROUTINE Stoke2

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE getp11(ind,nmu,thetad,P,theta,delta,p11)

    use MMCRAD1D_CONSTANTS, only: DTR, flag_interpol, EPS1
    use MMCRAD1D_INTERPOLATE

    implicit none

    !--- Input Arguments
    integer,intent(in)                            :: ind
    integer,intent(in)                            :: nmu
    real(8),intent(in), dimension(nmu)            :: thetad
    real(8),intent(in), dimension(nmu)            :: P
    real(8),intent(in)                            :: theta
    real(8),intent(in)                            :: delta
    real(8),intent(out)                           :: p11

    !--- Local Arguments
    integer                                 :: N1, N2, J0
    integer, dimension(1)                   :: ival
    integer, dimension(nmu)                 :: mask

    real(8)                                 :: deriv
    real(8)                                 :: p11ray(1)

    !---------------

    IF (ind .EQ. 1) THEN   ! PARTICLE

       !**************************************************************************
       !        Look for existing values (don't need to interpolate
       !**************************************************************************
       mask(:) = 0
       WHERE ((thetad(:) .GT. (theta - EPS1)) .AND.       &
            (thetad(:) .LT. (theta + EPS1)))
          mask(:) = 1
       END WHERE

       !**************************************************************************
       !    Interpolation of the coefficient of the phase matrix in theta
       !**************************************************************************
       IF (sum(mask(:)).GT. 0) THEN ! We don't need to interpolate

          ival = MAXLOC(mask(:))
          p11 = P(ival(1))

       ELSE    !WE NEED TO INTERPOLATE

          !**************************** ! Linear interpolation ! ************************
          IF (flag_interpol .EQ. 0) THEN  ! linear interpolation
             !--- Linear interpolation of the phase matrix on mu
             !--- Use a dichotomic method to find the 2 indice (J0+1 and J0) that 
             !    surround theta
             N1 = 1
             N2 = nmu
             DO WHILE ((N2-N1) .GT. 1)
                J0 = FLOOR(REAL(N1 + N2) / 2.0D0)
                IF (theta .GE. thetad(J0)) THEN
                   N1 = J0
                ELSE
                   N2 = J0
                ENDIF
             ENDDO
             IF ((theta .LT. thetad(J0)) .AND. (J0 .GT. 1)) J0 = J0 - 1

             p11 = LOG10(P(J0)) + (LOG10(P(J0+1)) - LOG10(P(J0))) * (theta - thetad(J0))&
                  / (thetad(J0+1) - thetad(J0))

             !************************* ! Polynomial interpolation !************************
          ELSE
             CALL SPLINE2(nmu,thetad,LOG10(P(:)),theta,P11,deriv)
          ENDIF

          p11 = 10.0D0**p11

       ENDIF

    ELSE  !RAYLEIGH PARTICLE
       !*************************************************************************
       !             Compute exact rayleigth coefficients
       !*************************************************************************
       call get_pray(theta, 1, delta, p11ray)
       p11 = p11ray(1)

    ENDIF

    RETURN

  END SUBROUTINE getp11

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  SUBROUTINE Compute_theta(cos_dir,cos_dir_sensor,theta)

    use MMCRAD1D_CONSTANTS, only: DTR, Pi

    implicit none

    !--- Input Arguments
    real(8), intent(in), dimension(3)    :: cos_dir
    real(8), intent(in), dimension(3)    :: cos_dir_sensor

    !--- Local Arguments
    real(8)                              :: x0, x1, y0, y1, z0, z1, cs

    !--- Output Arguments
    real(8),intent(out)                  :: theta

    !-------

    x0 = cos_dir(1)
    y0 = cos_dir(2)
    z0 = cos_dir(3)
    x1 = cos_dir_sensor(1)
    y1 = cos_dir_sensor(2)
    z1 = cos_dir_sensor(3)

    !--- Find cos(theta) and sin(theta)
    cs =  x1*x0 + y1*y0 + z1*z0  ! produit scalaire

    !--- Now compute theta
    !IF (cs .GT.  1.0D0) cs =  1.0D0
    !IF (cs .LT. -1.0D0) cs = -1.0D0
    theta = ACOS(cs) / DTR
    !IF (theta .GT. 180.0D0) theta = 180.0D0
    !IF (theta .LT.   0.0D0) theta =   0.0D0

    RETURN

  END SUBROUTINE Compute_theta

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE Compute_psi(v,cos_dir,cos_dir_sensor,N,cspsi,snpsi)

    ! Compute the rotation angle betwen the last and the 
    ! current scattering plane

    use MMCRAD1D_CONSTANTS, only: DTR, Pi

    implicit none

    !--- Input Arguments
    real(8),intent(in), dimension(3)                  :: v
    real(8),intent(in), dimension(3)                  :: cos_dir
    real(8),intent(in), dimension(3)                  :: cos_dir_sensor
    !--- Output Arguments
    real(8),intent(out)                               :: snpsi
    real(8),intent(out)                               :: cspsi
    real(8),intent(out), dimension(3)                 :: N

    !--- Local Arguments
    real(8), parameter                      :: EPS6 = 1.0D-06
    real(8)                                 :: x0, y0, y1, z0, den, prod_scal
    real(8),dimension(3)                    :: u


    !**************************************************************************
    !   Compute u and w
    !**************************************************************************
    u(:) = cos_dir(:)

    !--- chg name
    x0 = cos_dir_sensor(1)
    y0 = cos_dir_sensor(2)
    z0 = cos_dir_sensor(3)

    !**************************************************************************
    !   Compute the normal to the new scattering plan (sensor)
    !**************************************************************************
    N(1) = u(2)*z0 - u(3)*y0
    N(2) = u(3)*x0 - u(1)*z0
    N(3) = u(1)*y0 - u(2)*x0
    den = SQRT(N(1)*N(1) + N(2)*N(2) + N(3)*N(3))
    IF (den .LT. EPS6) THEN
       prod_scal = u(1)*x0+u(2)*y0+u(3)*z0
       IF (prod_scal .LT. 0.0D0) THEN  !diffusion vers l'avant
          N(:) = v(:)
       ELSE                            !diffusion vers l'arriere
          N(:) = -v(:)
       ENDIF
    ELSE
       N(:) = N(:) / den
    ENDIF

    !**************************************************************************
    !   Compute the scalar product between N and v, should give cspsi
    !**************************************************************************
    cspsi = v(1)*N(1) + v(2)*N(2) + v(3)*N(3)
    IF (cspsi .GT.  1.0D0) cspsi =  1.0D0
    IF (cspsi .LT. -1.0D0) cspsi = -1.0D0
    snpsi = DSQRT(1.0D0 - cspsi*cspsi)
    IF (snpsi .GT.  1.0D0) snpsi =  1.0D0

    !**************************************************************************
    !   Change of reference frame, look for the expression of cos_dir_sensor
    !   in the base linked to the photon before the scattering event = old
    !   scattering plan. 
    !   Let say that x1, y1, z1 are the new coordinate of cos_dir_sensor
    !**************************************************************************
    y1 = v(1)*x0 + v(2)*y0 + v(3)*z0

    !--- Sign of snpsi
    IF (y1 .LT. 0.0D0) THEN
       snpsi = -snpsi
    ENDIF

    RETURN

  END SUBROUTINE Compute_psi

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE compute_Trans(nlay,ilayer,thetav,z_photon,z_ref,cext,Trans)

    use MMCRAD1D_CONSTANTS, only: DTR

    implicit none

    !--- Input Arguments
    integer,intent(in)                      :: nlay
    integer,intent(in)                      :: ilayer
    real(8),intent(in)                      :: thetav
    real(8),intent(in)                      :: z_photon
    real(8),intent(in) , dimension(0:nlay)  :: z_ref
    real(8),intent(in) , dimension(nlay)    :: cext
    real(8),intent(out)                     :: Trans

    !--- Local Arguments
    real(8) :: tau
    integer :: i

    !------------------------

    ! First compute the contribution from the layer 
    ! the photon is embedded in
    tau =  cext(ilayer) * (z_ref(ilayer) - z_photon) / cos(thetav*DTR)

    ! add contribution form layers on top of the one 
    ! the photon is located in
    DO i = ilayer+1, nlay
       tau = tau + ( cext(i) * (z_ref(i) - z_ref(i-1))  / cos(thetav*DTR) )
    ENDDO

    Trans = DEXP(-tau)

    RETURN

  END SUBROUTINE Compute_Trans

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE get_prfs(nmu, nlay, mu, p_ptcle, prf_ptcle)

    USE MMCRAD1D_CONSTANTS

    implicit none

    !--- Input Arguments
    integer, intent(in)  :: nmu 
    integer, intent(in)  :: nlay
    real(8), intent(in)  :: mu(nmu)
    real(8), intent(in)  :: p_ptcle(nlay,nmu)

    real(8), intent(out) :: prf_ptcle(nlay,nmu)


    !--- Local Arguments
    integer :: k, j
    real(8) :: summ

    !--------------------

    !  Compute the Probability repartition function for ptcle
    DO k = 1, nlay
       !--- Is the phase function normalized ?
       summ = -MC_XINTEG2(1, nmu, nmu, mu, p_ptcle(k,:))
       !WRITE(*,*)'int(P11_ptcle(',k,'))=', summ
       prf_ptcle(k,1) = 0.0D0
       DO j = 2, nmu
          prf_ptcle(k,j) = -MC_XINTEG2(1, j, nmu, mu, p_ptcle(k,:)) / summ
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE get_prfs

  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE scatter(ithread, nmat, nelem, surface, &  
       ddis, epsilon_ddis, p_ddis, prf_ddis,        &
       cos_dir_sensor,                              &
       ind, nmu, mu, thetad, p, prf, delta_ray,     &
       wght, Stk, v, u, cos_dir)

    !=============================================
    !  for DDISing see
    !   see Buras and Mayer, 2011, JQSRT, 112, 437
    !=============================================

    USE MMCRAD1D_CONSTANTS
    USE MMCRAD1D_ZIGGURAT

    IMPLICIT NONE 

    INTEGER, INTENT(IN)  :: ithread
    INTEGER, INTENT(IN)  :: nmat
    INTEGER, INTENT(IN)  :: nelem
    CHARACTER(len=*), INTENT(IN) :: surface
    LOGICAL, INTENT(IN)  :: ddis
    REAL(8), INTENT(IN)  :: epsilon_ddis
    REAL(8), INTENT(IN)  :: p_ddis(nmu)
    REAL(8), INTENT(IN)  :: prf_ddis(nmu)
    REAL(8), INTENT(IN)  :: cos_dir_sensor(3)
    INTEGER, INTENT(IN)  :: ind
    INTEGER, INTENT(IN)  :: nmu
    REAL(8), INTENT(IN)  :: mu(nmu)
    REAL(8), INTENT(IN)  :: thetad(nmu)
    REAL(8), INTENT(IN)  :: prf(nmu)
    REAL(8), INTENT(IN)  :: p(nelem,nmu)
    REAL(8), INTENT(IN)  :: delta_ray
    REAL(8), INTENT(INOUT)  :: wght
    REAL(8), INTENT(INOUT)  :: stk(nmat)
    REAL(8), INTENT(INOUT)  :: v(3)
    REAL(8), INTENT(INOUT)  :: u(3)
    REAL(8), INTENT(INOUT)  :: cos_dir(3)

    ! local variables
    REAL(8) :: ctheta, theta 
    REAL(8) :: ctheta_ddis, theta_ddis
    REAL(8) :: phi
    REAL(8) :: c1, s1, c2, s2
    REAL(8) :: P11, P11ddis
    REAL(8) :: rn

    LOGICAL :: natural_scat
    REAL(8) :: wght_coef

    REAL(8) :: stk_in(nmat), stk_tmp(nmat)
    REAL(8) :: cos_dir_in(3)
    REAL(8) :: u_in(3)
    REAL(8) :: v_in(3)
    REAL(8) :: wght_in

    REAL(8) :: v_sensor(3)

    !--------------   

    if (ind .lt. 2) then 

       ! gas or particle scattering
       cos_dir_in = cos_dir
       u_in       = u
       v_in       = v
       stk_in     = stk
       wght_in    = wght

       ! will the photon be ddised ?
       natural_scat = .true.
       IF (ddis .eqv. .true.) then
          !CALL RANDOM_NUMBER(rn)
          rn = par_ziguni(ithread)
          if (rn .le. epsilon_ddis) then
             natural_scat = .false.
          endif
       endif

       if (natural_scat .eqv. .true.) then 

          !********************************************************************
          !********************  Natural scattering ***************************
          !********************************************************************
          !    Compute the new direction of the photon from the scattering 
          !    matrix of the particle reached.
          !    (1) Compute the angle between incident and scattered photon
          !        theta (in degrees) 
          !    (2) Compute the new direction cosines, actually 
          !        this is done in (5).
          !    (3) Rotate the incident Stocke vector to bring it in the 
          !        scattering plan. Rotation of angle i1=phi (e.g. Euler's
          !        method in Ramella-Roman et al., 2005a).
          !    (4) Compute the new scattered Stoke vector from the incident one 
          !        and the phase matrix of the particle at scattering angle 
          !        theta (in the scattering plan).
          !    (5) Keep track of these rotations on the (u,v,w) vectors
          !        u belong parallel to the reference plan
          !        v        perpendicular
          !        w belong in the direction of photon propagation
          !********************************************************************
          !--- (1)
          CALL newdirection(ithread, ind, nmu, mu, prf, ctheta, theta, phi)

          !<CHG> 04/09/07
          !      cos_dir is equal to the unit vector w of the Euler's method
          !      I am using the rotational matrix define by Ramella-Roman because
          !      it is an expression more general than the one I am using from 
          !      Hovenier or other.... Therefore cos_dir is computed in (5)!
          !--- (2)
          !</CHG>

          IF (nmat .EQ. 1) THEN ! Don't need to take into account the rotation of the Stoke vector !!!

             !--- (4)
             stk(:) = stk_in(:)

             !--- (5)
             ! I am using the rotational matrix define in Ramella-Roman et al., 
             ! 2005a
             ! (a) rotation of v around u of phi
             ! (b) rotation of u around new v of Theta
             !-- (a)
             CALL rot_euler(v_in, u_in, phi, v)
             !-- (b)
             CALL rot_euler(u_in, v, theta, u)
             cos_dir(:) = u(:)

          ELSE ! (nmat>1) Need to take into account the rotation of the Stoke vector!!!

             !--- (3)
             c1 = COS(2.0D0*phi*DTR)
             s1 = SIN(2.0D0*phi*DTR)
             CALL Stoke_rot(nmat,c1,s1,Stk_in,Stk_tmp)
             !--- (4)
             call stoke(ind,nmat,nelem,nmu,thetad,p,Stk_tmp,theta,delta_ray,stk,p11)

             !--- (5)
             ! I am using the rotational matrix define in Ramella-Roman et al., 
             ! 2005a
             ! (a) rotation of v around u of phi
             ! (b) rotation of u around new v of Theta
             !-- (a)
             CALL rot_euler(v_in, u_in, phi, v)
             !-- (b)
             CALL rot_euler(u_in, v, theta, u)
             cos_dir(:) = u(:)

          ENDIF ! on nmat

       ELSE

          !********************************************************************
          !********************       DDISing      ***************************
          !********************************************************************
          ! we get new director cosines by getting 
          ! a (theta, phi) from P_ddis scattering and starting from the 
          ! detector direction.
          ! Note we need a V vector corresponding to the 
          !      cos_dir_sensor to perform the (theta, phi) rotation. 
          !      Since phi is random, we can just take the vectorial
          !      produce cos_dir_sensor by Oz
          v_sensor(1) =  cos_dir_sensor(2)
          v_sensor(2) = -cos_dir_sensor(1)
          v_sensor(3) =  0.0D0
          ! NOTE : ind = 1 cause we want to use the prf_ddis to get new direction
          !        i.e. a ptcle probability repartition function
          CALL newdirection(ithread, 1, nmu, mu, prf_ddis, ctheta_ddis, theta_ddis, phi)
          CALL rot_euler(v_sensor, cos_dir_sensor, phi, v)
          CALL rot_euler(cos_dir_sensor, v, theta_ddis, u)
          cos_dir = u

          if (nmat .eq. 1) then 

             ! no need to rotate the Stoke vector
             stk(:) = stk_in(:)

          else

             CALL Compute_psi(v_in, cos_dir_in, cos_dir, v, &
                  c1, s1)
             ! v_final must be v_n the normal to the 
             ! scattering angle
             ! We then do not need to rotate v_in by phi (psi)
             ! around cos_dir_in
             ! Rotate the stoke Vector
             c2 = 2.0D0*c1*c1 - 1.0D0
             s2 = 2.0D0*c1*s1
             CALL Stoke_rot(nmat,c2,s2,Stk_in,Stk_tmp)

             ! compute scattering angle between the old and new photon direction
             CALL Compute_theta(cos_dir,cos_dir_in,theta)
             ! update the stoke vector as if it had been scattered naturally
             call stoke(ind,nmat,nelem,nmu,thetad,p,Stk_tmp,theta,delta_ray,stk,p11)

          endif

       ENDIF ! TEST on natural scattering or DDISing

       ! update photon weight
       IF (ddis .eqv. .false.) THEN

          wght = wght_in 

       ELSE 

          IF (natural_scat .eqv. .true.) THEN

             IF (nmat .eq. 1) then
                ! get the phase function value p11 at theta
                ! for the weight coeff.
                call getp11(ind, nmu, thetad, P(1,:), theta, delta_ray, p11)
             ENDIF

             ! compute scattering angle between new phot direction and sensor dir.
             CALL Compute_theta(cos_dir,cos_dir_sensor,theta_ddis)
             ! get the phase function value p11ddis at theta_ddis
             ! for the weight coeff.
             call getp11(1, nmu, thetad, p_ddis, theta_ddis, delta_ray, p11ddis)

          ELSE

             ! get the phase function value p11ddis at theta_ddis
             ! for the weight coeff.
             call getp11(1,nmu,thetad,p_ddis,theta_ddis,delta_ray,p11ddis)

             IF (nmat .eq. 1) then
                ! compute scattering angle between the old and new phot direction
                CALL Compute_theta(cos_dir,cos_dir_in,theta)
                ! get the phase function value p11 at theta
                ! for the weight coeff.
                call getp11(ind, nmu, thetad, P(1,:), theta, delta_ray, p11)
             ENDIF

          ENDIF

          wght_coef = p11 / ( (1.0d0-epsilon_ddis)*p11  +  epsilon_ddis*p11ddis )
          wght      = wght_in * wght_coef

          !write(*,*) natural_scat, theta, p11, theta_ddis, p11ddis, wght_coef, wght_in      

       ENDIF

    ELSE

       ! surface reflection

       call surface_scatter(ithread, surface,  &
            nmat, cos_dir, v, u, Stk)

    ENDIF

  END SUBROUTINE scatter

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  subroutine surface_scatter(ithread, surface,    &
       nmat, cos_dir, v, u, stk)

    USE MMCRAD1D_CONSTANTS
    USE MMCRAD1D_ZIGGURAT
    USE MMCRAD1D_BRDF
    USE MMCRAD1D_INTERPOLATE, only : MCLININTPOL, MCBILININTPOL

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ithread
    CHARACTER(len=*), INTENT(IN) :: surface
    INTEGER, INTENT(IN) :: nmat
    REAL(8), INTENT(INOUT) :: cos_dir(3)
    REAL(8), INTENT(INOUT) :: v(3) 
    REAL(8), INTENT(INOUT) :: u(3)
    REAL(8), INTENT(INOUT) :: stk(nmat)

    ! local variables
    REAL(8) :: cos_dir_old(3)
    REAL(8) :: Stk_tmp(nmat)
    REAL(8) :: xr, rn
    REAL(8) :: teta_in
    REAL(8) :: phi_in
    REAL(8) :: teta_ref
    REAL(8) :: phi_ref
    REAL(8) :: dphi
    integer :: i
    real(8) :: v_tmp(3)
    real(8) :: oz(3)
    real(8) :: cspsi,snpsi
    real(8), allocatable :: prf_mu(:)
    real(8), allocatable :: prf_phi(:)

    !----------------

    IF (TRIM(ADJUSTL(surface)) .eq. 'lambert') then

       ! PHOTON REACHING LAMBERTIAN SURFACE
       ! Isotropic restitution of the intensity
       ! zenith angle
       !CALL RANDOM_NUMBER(xr)
       !xr = ACOS(SQRT(xr))
       xr = ACOS(SQRT(par_ziguni(ithread))) ! See for example Mayer (2009)
       ! azimuth angle
       !CALL RANDOM_NUMBER(rn)
       rn = par_ziguni(ithread) * PI * 2.0D0
       cos_dir(1) = SIN(xr)*COS(rn)
       cos_dir(2) = SIN(xr)*SIN(rn)
       cos_dir(3) = COS(xr)
       u(:) = cos_dir(:)
       ! !!!!!!!!!!!!!!!!
       ! this "v" computation is only good for a Lambertian surface
       ! we just need v to be orthogonal to cos_dir since we do not
       ! want to keep track of polarization...
       ! (However we still need v to properly rotate cos_dir: Euler method)
       ! Here we take the produce cos_dir by Oz
       v(1) = u(2)
       v(2) = -u(1)   
       v(3) =  0.0D0     
       ! renormalization
       v(:) = v(:) / SQRT(v(1)**2.0D0+v(2)**2.0D0)
       ! !!!!!!!!!!!!!!!!
       if (nmat .gt. 1) Stk(2:nmat) = 0.0D0  ! depolarisation

    elseif (TRIM(ADJUSTL(surface)) .eq. 'brdf') then

       allocate(prf_mu(ntetasurf), prf_phi(nphisurf))  

       cos_dir_old(:) = cos_dir(:)
       teta_in = ACOS(-cos_dir_old(3)) / DTR
       phi_in  = ATAN2(cos_dir_old(2), cos_dir_old(1)) / dtr

       ! we first compute the new direction
       if (flag_interpol.eq.0) then
          ! linear interpolation
          ! interpolate prf_mu_surf to the desired incident angle
          do i = 1, ntetasurf
             prf_mu(i) = MCLININTPOL(prf_mu_surf(:,i), tetasurf, ntetasurf, teta_in )
          enddo
          ! first get a theta direction
          xr       = par_ziguni(ithread) 
          teta_ref = MCLININTPOL(tetasurf(:), prf_mu(:), ntetasurf, xr )
          if ((teta_ref.gt.90.0).or.(teta_ref.lt.0.0)) then
             write(*,*) teta_ref, xr
             write(*,*) ' '
             do i = 1, ntetasurf
                write(*,*) tetasurf(i), prf_mu(i)
             enddo
          endif
          ! interpolate prf_phi_surf to the desired incident and reflected angles
          do i = 1, nphisurf
             prf_phi(i) = MCBILININTPOL(prf_phi_surf(:,:,i), tetasurf, tetasurf, ntetasurf, ntetasurf,&
                  teta_in, teta_ref)
          enddo
          ! get dphi 
          xr = par_ziguni(ithread) 
          ! get the corresponding teta_ref value
          dphi  =  MCLININTPOL(phisurf(:), prf_phi(:), nphisurf, xr )
       else
          write(*,*)' non-linear interp not implemented for surface BRDF reflect'
          STOP
       endif

       phi_ref = phi_in + dphi
       phi_ref = MOD(phi_ref, 360.0D0)
       !write(*,*) teta_in, phi_in, teta_ref, phi_ref, dphi
       ! Compute the new direction cosine
       cos_dir(1) = SIN(DTR*teta_ref)*COS(DTR*phi_ref)
       cos_dir(2) = SIN(DTR*teta_ref)*SIN(DTR*phi_ref)
       cos_dir(3) = COS(DTR*teta_ref)
       u(:) = cos_dir(:)

       ! ============================================================
       ! compute the v normal to the meridian
       ! and get the psi to rotate the Stoke vector with
       oz(1)=0.0
       oz(2)=0.0
       oz(3)=1.0
       CALL Compute_psi(v,cos_dir_old,oz,v_tmp,&
            cspsi,snpsi)

       ! rotate  the v around oz to get the new v
       ! (still in normal to the meridian plane but also normal to the new cos_dir)
       CALL rot_euler(v_tmp, oz, dphi, v)

       IF (nmat .gt. 1) THEN 

          !--- (3)
          CALL Stoke_rot(nmat,cspsi,snpsi,Stk,Stk_tmp)
          !--- (4)
          CALL surface_stoke(nmat, .true., Stk_tmp(:), teta_in, teta_ref, dphi, Stk(:))

       ENDIF ! on nmat
       ! ============================================================
       deallocate(prf_mu, prf_phi)  

    endif ! lambertian or BRDF

  END subroutine surface_scatter

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE local_estimate(nmat, nelem, surface, nlay, cext_ly, z_ref, tau_tot, &
       ind, nmu, thetad, p, delta_ray,                                           &
       wght, stk, z_phot, cos_dir, v, ilayer,                                    &    
       thetav, cos_dir_sensor, le)

    USE MMCRAD1D_CONSTANTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nmat
    INTEGER, INTENT(IN) :: nelem
    CHARACTER(len=*), INTENT(IN) :: surface
    INTEGER, INTENT(IN) :: nlay
    REAL(8), INTENT(IN) :: cext_ly(nlay)
    REAL(8), INTENT(IN) :: z_ref(nlay+1)
    REAL(8), INTENT(IN) :: tau_tot
    INTEGER, INTENT(IN) :: ind
    INTEGER, INTENT(IN) :: nmu
    REAL(8), INTENT(IN) :: thetad(nmu)
    REAL(8), INTENT(IN) :: p(nelem,nmu)
    REAL(8), INTENT(IN) :: delta_ray
    REAL(8), INTENT(IN) :: wght
    REAL(8), INTENT(IN) :: stk(nmat)
    REAL(8), INTENT(IN) :: z_phot
    REAL(8), INTENT(IN) :: cos_dir(3)
    REAL(8), INTENT(IN) :: v(3)
    INTEGER, INTENT(IN) :: ilayer
    REAL(8), INTENT(IN) :: thetav
    REAL(8), INTENT(IN) :: cos_dir_sensor(3)
    REAL(8), INTENT(OUT) :: le(nmat)

    !-- local variables
    real(8) :: stk_tmp(nmat), stk_old(nmat)
    real(8) :: theta_sensor
    real(8) :: v_tmp(3)
    real(8) :: w_tmp(3)
    real(8) :: c2, s2, c3, s3
    real(8) :: trans
    real(8) :: cspsi,snpsi
    real(8) :: epsi

    if (ind .lt. 2) then

       ! ===========================
       !    gas or particles LE
       ! ===========================

       !***********************************************
       !   As in Noel et al., 2002, each scattering event is contributing
       !   to the recorded Light (reach the sensor).
       !   (1) Compute the angles theta_sensor, scattering angle from the
       !       cloud particle to the sensor.
       !   (2) From cos_dir_old, cos_dir_sensor and theta_sensor, compute the 
       !       angle (psi) between the old scattering plan and the new one.
       !   (3) Rotate the Stoke vector of psi to bring it in the new 
       !       scattering plan
       !   (4) Compute the new stoke vector from the phase matrix at 
       !       theta_sensor.
       !   (5) Compute the angle epsi (equivalent to i2) between the last 
       !       scattering plan and the meridian plan and rotate the Stoke 
       !       vector. 
       !   (6) Compute the transmission from this event to the Sensor
       !*************************************************

       !--- (1)
       CALL Compute_theta(cos_dir,cos_dir_sensor,theta_sensor)

       IF (nmat .EQ. 1) THEN   ! Don't need to take into account the rotation of the Stoke vector

          !--- (4)
          call Stoke2(nmat,ind,nelem,nmu,thetad,p,Stk,theta_sensor,delta_ray,&
               Stk_tmp)

       ELSE ! Need to take into account the rotation of the Stoke vector

          !--- (2) here v_tmp is the normal to the new scattering plan that 
          !        contain cos_dir_sensor....
          CALL Compute_psi(v,cos_dir,cos_dir_sensor,v_tmp,&
               cspsi,snpsi)

          !--- (3)
          c2 = 2.0D0*cspsi*cspsi - 1.0D0
          s2 = 2.0D0*cspsi*snpsi
          CALL Stoke_rot(nmat,c2,s2,Stk,Stk_tmp)

          !--- (4)
          Stk_old(:) = Stk_tmp(:)
          call Stoke2(nmat,ind,nelem,nmu,thetad,p,Stk_old,theta_sensor,delta_ray,&
               Stk_tmp)

          !--- (5)
          !- Before compute w_tmp, vector parallel to the new scattering plan
          !  here cos_dir_sensor and u_tmp are identical
          w_tmp(1) = v_tmp(2)*cos_dir_sensor(3) -                    &
               v_tmp(3)*cos_dir_sensor(2) 
          w_tmp(2) = v_tmp(3)*cos_dir_sensor(1) -                    &
               v_tmp(1)*cos_dir_sensor(3) 
          w_tmp(3) = v_tmp(1)*cos_dir_sensor(2) -                    &
               v_tmp(2)*cos_dir_sensor(1)
          !- angle epsi = i2
          epsi=ATAN2(v_tmp(3), w_tmp(3))

          !- Want v, after rotation to the new meridian plan, in the same 
          !  direction that the Normal to this plan
          IF (w_tmp(3) .GT. 0.0D0) THEN
             IF (v_tmp(3) .LT. 0.0D0) THEN
                !  WRITE(6,*)'epsi cas 2=',epsi
                epsi = Pi + epsi
             ELSE
                !  WRITE(6,*)'epsi cas 4=',epsi
                epsi = epsi - Pi
             ENDIF
          ENDIF
          !- Rotation of the Stoke vector
          Stk_old(:) = Stk_tmp(:)
          c3 = COS(2.0D0*epsi)
          s3 = SIN(2.0D0*epsi)
          CALL Stoke_rot(nmat,c3,s3,Stk_old,Stk_tmp)

       ENDIF ! test on nmat

       !--- (6)
       CALL Compute_trans(nlay, ilayer, thetav, z_phot, z_ref, cext_ly, trans)
       !--- Now compute the transmission
       !    Rq: * we divided by cos(thetav) to account for the inclination between areas
       !          thetav is the angle between the direction of the virtual photon and the 
       !        * we divided by the factor 4 pi because we need to deal with 
       !          normalized phase function
       trans = trans / COS(DTR*thetav) / (4.0D0 * pi)

       !--- Add the Stoke vector reaching the Sensor
       le(:) = Stk_tmp(:) * wght * trans

    else

       ! LOCAL ESTIMATE ON SURFACE  
       CALL surface_le(surface,             &
            nmat, thetav, tau_tot,          &
            wght, stk,  cos_dir, v,         &
            cos_dir_sensor,                 &
            le)

    endif ! surface LE

  END SUBROUTINE local_estimate

  ! ========================================================================

  SUBROUTINE surface_le(surface,       &
       nmat, thetav, tau_tot,          &
       wght, stk, cos_dir, v,          &
       cos_dir_sensor,                 &
       le) 

    USE MMCRAD1D_CONSTANTS

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: surface
    INTEGER, INTENT(IN) :: nmat
    REAL(8), INTENT(IN) :: thetav
    REAL(8), INTENT(IN) :: tau_tot
    REAL(8), INTENT(IN) :: wght
    REAL(8), INTENT(IN) :: stk(nmat)
    REAL(8), INTENT(IN) :: cos_dir(3)
    REAL(8), INTENT(IN) :: v(3)
    REAL(8), INTENT(IN) :: cos_dir_sensor(3)
    REAL(8), INTENT(OUT) :: le(nmat)

    ! local variables
    real(8) :: stk_tmp(nmat), stk_old(nmat)
    real(8) :: v_tmp(3)
    real(8) :: w_tmp(3)
    real(8) :: oz(3)
    real(8) :: c2, s2
    real(8) :: trans
    real(8) :: cspsi,snpsi

    real(8) :: phi_sensor
    real(8) :: phi_phot
    real(8) :: dphi

    trans = EXP(-tau_tot/cos(thetav*DTR)) / pi

    IF (TRIM(ADJUSTL(surface)) .eq. 'lambert') then

       ! Lambertion surface
       ! In case of the lambertian surface, only the first stoke 
       ! parameter is modified because of the depolarization
       le(1) = Stk(1) * trans * wght  
       if (nmat .gt. 1) le(2:nmat) = 0.0D0  ! depolarisation

    elseif (TRIM(ADJUSTL(surface)) .eq. 'brdf') then

       ! BRDF surface
       ! get the dphi azimuthal angle between incident photon and the detector direction
       phi_sensor = ATAN2(cos_dir_sensor(2), cos_dir_sensor(1))
       phi_phot   = ATAN2(cos_dir(2), cos_dir(1))
       dphi       = (phi_sensor - phi_phot) / dtr
       if (dphi.lt.0.0D0) dphi = 360.0D0+dphi

       if (nmat.eq.1) then

          !--- (4)
          CALL surface_stoke(nmat, .false., Stk(1), ACOS(-cos_dir(3))/dtr, ACOS(cos_dir_sensor(3))/dtr, dphi, Stk_tmp(1))
          le(1) = Stk_tmp(1) * trans * wght 

       else

          ! ============================================================

          !--- (2) here v_tmp is normal to the meridian
          oz(1) = 0.0
          oz(2) = 0.0
          oz(3) = 1.0
          CALL Compute_psi(v,cos_dir,w_tmp,v_tmp,&
               cspsi,snpsi)

          !--- (3)
          c2 = 2.0D0*cspsi*cspsi - 1.0D0
          s2 = 2.0D0*cspsi*snpsi
          CALL Stoke_rot(nmat,c2,s2,Stk,Stk_tmp)

          !--- (4)
          Stk_old(:) = Stk_tmp(:) 
          CALL surface_stoke(nmat, .false., Stk_old(:), ACOS(-cos_dir(3))/dtr, ACOS(cos_dir_sensor(3))/dtr, dphi, Stk_tmp(:))

          le(:) = Stk_tmp(:) * trans * wght 
          ! ============================================================

       endif

    endif ! on lambertian surface or not

  END SUBROUTINE surface_le

  ! ========================================================================

  SUBROUTINE surface_stoke(nmat, norm, Stoke_0, teta_i, teta_r, phi, Stoke_1)

    USE MMCRAD1D_CONSTANTS, only : flag_interpol
    USE MMCRAD1D_BRDF, only : brdf, ntetasurf, tetasurf, nphisurf, phisurf
    USE MMCRAD1D_INTERPOLATE, only : MCTRILININTPOL

    IMPLICIT NONE

    integer, intent(in)  :: nmat
    logical, intent(in)  :: norm
    real(8), intent(in)  :: Stoke_0(nmat)
    real(8), intent(in)  :: teta_i ! in degree
    real(8), intent(in)  :: teta_r ! in degree
    real(8), intent(in)  :: phi    ! in degree
    real(8), intent(out) :: Stoke_1(nmat)

    ! local variables 
    real(8) :: brdf_local(nmat,nmat)

    ! compute BRDF_LOCAL
    IF (flag_interpol .EQ. 0) THEN  ! linear interpolation

       brdf_local(1,1) = MCTRILININTPOL(BRDF(1,1,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
            teta_i, teta_r, phi)
       if (nmat.gt.1) then
          brdf_local(1,2) = MCTRILININTPOL(BRDF(1,2,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(1,3) = MCTRILININTPOL(BRDF(1,3,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(2,1) = MCTRILININTPOL(BRDF(2,1,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(2,2) = MCTRILININTPOL(BRDF(2,2,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(2,3) = MCTRILININTPOL(BRDF(2,3,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(3,1) = MCTRILININTPOL(BRDF(3,1,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(3,2) = MCTRILININTPOL(BRDF(3,2,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(3,3) = MCTRILININTPOL(BRDF(3,3,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
       endif
       if (nmat.gt.3) then
          brdf_local(1,4) = MCTRILININTPOL(BRDF(1,4,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(2,4) = MCTRILININTPOL(BRDF(2,4,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(3,4) = MCTRILININTPOL(BRDF(3,4,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(4,4) = MCTRILININTPOL(BRDF(4,4,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(4,3) = MCTRILININTPOL(BRDF(4,3,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(4,2) = MCTRILININTPOL(BRDF(4,2,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
          brdf_local(4,1) = MCTRILININTPOL(BRDF(4,1,:,:,:), tetasurf, tetasurf, phisurf, ntetasurf, ntetasurf, nphisurf, &
               teta_i, teta_r, phi)
       endif

    ELSE

       write(*,*) ''
       write(*,*) '      (surface_stoke) : ERROR'
       write(*,*) '                        Only linear interpolation of the BRDF is allowed'
       write(*,*) ''
       STOP

    ENDIF

    if (norm)  brdf_local(:,:) = brdf_local(:,:) / brdf_local(1,1)

    ! apply the reflection matrix to the stokes vector
    if (nmat.eq.1) then
       Stoke_1(1) = Stoke_0(1) * brdf_local(1,1)
    elseif(nmat.eq.3) then
       Stoke_1(1) = Stoke_0(1) * brdf_local(1,1) + Stoke_0(2) * brdf_local(1,2) + Stoke_0(3) * brdf_local(1,3)
       Stoke_1(2) = Stoke_0(1) * brdf_local(2,1) + Stoke_0(2) * brdf_local(2,2) + Stoke_0(3) * brdf_local(2,3)        
       Stoke_1(3) = Stoke_0(1) * brdf_local(3,1) + Stoke_0(2) * brdf_local(3,2) + Stoke_0(3) * brdf_local(3,3)        
    elseif(nmat.eq.4) then
       Stoke_1(1) = Stoke_0(1) * brdf_local(1,1) + Stoke_0(2) * brdf_local(1,2) + Stoke_0(3) * brdf_local(1,3) &
            + Stoke_0(4) * brdf_local(1,4)
       Stoke_1(2) = Stoke_0(1) * brdf_local(2,1) + Stoke_0(2) * brdf_local(2,2) + Stoke_0(3) * brdf_local(2,3) &
            + Stoke_0(4) * brdf_local(2,4)
       Stoke_1(3) = Stoke_0(1) * brdf_local(3,1) + Stoke_0(2) * brdf_local(3,2) + Stoke_0(3) * brdf_local(3,3) &
            + Stoke_0(4) * brdf_local(3,4)
       Stoke_1(4) = Stoke_0(1) * brdf_local(4,1) + Stoke_0(2) * brdf_local(4,2) + Stoke_0(3) * brdf_local(4,3) &
            + Stoke_0(4) * brdf_local(4,4)
    endif

  END SUBROUTINE surface_stoke

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE get_wspl(nmu, thetad, pddis,  &
       wght, cos_dir,                      &      
       cos_dir_sensor,                     &
       wspl)

    ! see Buras and Mayer, 2011 (sect 2.4)

    use MMCRAD1D_CONSTANTS
    use MMCRAD1D_INTERPOLATE

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nmu
    REAL(8), INTENT(IN)  :: thetad(nmu)
    REAL(8), INTENT(IN)  :: pddis(nmu)
    REAL(8), INTENT(IN)  :: wght
    REAL(8), INTENT(IN)  :: cos_dir(3)
    REAL(8), INTENT(IN)  :: cos_dir_sensor(3)
    REAL(8), INTENT(OUT) :: wspl

    !-- local variables
    real(8) :: theta
    integer, dimension(1)                  :: ival
    integer, dimension(nmu)                :: mask
    integer                                :: N1, N2, J0
    real(8)                                :: deriv
    real(8), dimension(nmu)                :: Pii
    real(8)                                :: coef_P

    !***********************************************

    ! ------------
    ! compute the scattering angle (angle between photon direction and the sensor)
    CALL Compute_theta(cos_dir,cos_dir_sensor,theta)

    ! ------------
    ! Get the value for the scattering function corresponding to the 
    ! scattering angle

    ! look for existing values (don't need to interpolate)
    mask(:) = 0
    WHERE ((thetad(:) .GT. (theta - EPS1)) .AND.       &
         (thetad(:) .LT. (theta + EPS1)))
       mask(:) = 1
    END WHERE

    IF (sum(mask(:)).GT. 0) THEN ! We don't need to interpolate

       ival   = MAXLOC(mask(:))
       coef_P = LOG10(pddis(ival(1)))

    ELSE    !WE NEED TO INTERPOLATE

       Pii(:) = LOG10(pddis(:))

       IF (flag_interpol .EQ. 0) THEN  ! linear interpolation
          !--- Linear interpolation of the phase matrix on mu
          N1 = 1
          N2 = nmu
          DO WHILE ((N2-N1) .GT. 1)
             J0 = FLOOR(REAL(N1 + N2) / 2.0D0)
             IF (theta .GE. thetad(J0)) THEN
                N1 = J0
             ELSE
                N2 = J0
             ENDIF
          ENDDO
          IF ((theta .LT. thetad(J0)) .AND. (J0 .GT. 1)) J0 = J0 - 1
          coef_P = Pii(J0) + (Pii(J0+1) - Pii(J0)) * (theta - thetad(J0))&
               / (thetad(J0+1) - thetad(J0))
       ELSE ! Polynomial interpolation
          CALL SPLINE2(nmu,thetad,Pii,theta,coef_P,deriv)
       ENDIF

    ENDIF

    wspl = 10.0D0**coef_P * wght

    ! if (wspl .gt. 3.0) then
    !    write(*,*) 'wspl =', wspl
    !    write(*,*) 'cos_dir        =', cos_dir 
    !    write(*,*) 'cos_dir_sensor =', cos_dir_sensor
    !    write(*,*) 'theta        =', theta
    !    write(*,*) 'Pddis(theta) =', 10.0D0**coef_P
    !    write(*,*) 'wght =', wght
    !    STOP 
    ! endif

  END SUBROUTINE get_wspl

  !SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  SUBROUTINE move_photon(ithread,                &
       surface,                                  &
       nlay, z_ref, cext_ly, fracext_gas_ly,     &
       cos_dir, ilayer, z_phot, leave,           &
       ind)

    ! This subroutine is used to propagate the photon to its next location
    ! Also provided is the photon status:
    !     leave = 0 --> photon still in the atmosphere (it will then interact an other time)
    !     leave = 1 --> photon leaves atmosphere from the bottom
    !     leave = 2 --> photon leaves atmosphere from the top
    ! and the kind of target that stopped the photon (ind)
    !     ind = 0  gas 
    !     ind = 1  particle (cloud or aerosols)
    !     ind = 2  surface

    USE MMCRAD1D_ZIGGURAT

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ithread
    CHARACTER(len=*), INTENT(IN) :: surface
    INTEGER, INTENT(IN) :: nlay
    REAL(8), INTENT(IN) :: z_ref(0:nlay)
    REAL(8), INTENT(IN) :: cext_ly(nlay)
    REAL(8), INTENT(IN) :: fracext_gas_ly(nlay)
    REAL(8), INTENT(IN) :: cos_dir(3)
    INTEGER, INTENT(INOUT) :: ilayer
    REAL(8), INTENT(INOUT) :: z_phot
    INTEGER, INTENT(OUT) :: leave
    INTEGER, INTENT(OUT) :: ind

    ! local variables
    REAL(8) :: taulimit
    REAL(8) :: dlimit
    REAL(8) :: tau
    REAL(8) :: xr
    REAL(8) :: z_to_limit
    REAL(8) :: old_z_phot
    REAL(8) :: d

    !------------------

    taulimit = 0.0D0
    dlimit   = 0.0D0

    !***********************************************************************
    !      Compute where should the photon bit a particle or molecule (Monte Carlo)
    !      Here cext_ly is the total extinction coefficient (rayleight + particles)
    !***********************************************************************
    !--- Need a random number between 0 and 1
    !CALL RANDOM_NUMBER(xr)
    xr = par_ziguni(ithread)

    !--- distance d = (- ln(xr) / cext_ly) of interaction
    !    optical depth tau of interaction
    tau = -DLOG(xr)

500 CONTINUE

    !--- rescale tau (remove the optical depth covered in the previous layer)
    tau = tau - taulimit

    !***********************************************************************
    !      Compare it to the limit of the layer from the photon position
    !***********************************************************************
    !--- layer limit depend of the photon direction, ie ascending or descending
    !    z_to_limit = distance to the top or base of the layer on the z axis
    !    dlimit     = distance to the top or the base in the photon direction
    IF (cos_dir(3) .GT. 0.0D0) THEN !ascending
       z_to_limit = z_ref(ilayer) - z_phot    ! sign does not matter
    ELSE  ! descending
       z_to_limit = z_ref(ilayer-1) - z_phot  ! sign does not matter
    ENDIF
    dlimit   =  DABS(z_to_limit / cos_dir(3))
    taulimit = dlimit * cext_ly(ilayer)

    !--- Photon interact in the layer ilayer or no, leaving atmosphere or no
    IF (tau .GE. taulimit) THEN 

       !*************************************
       !   PHOTON IS LEAVING THE LAYER
       !*************************************
       !    leave = 0 --> photon still in the atmosphere goes in an other layer
       !    leave = 1 --> photon leaves atmosphere from the bottom
       !    leave = 2 --> photon leaves atmosphere from the top

       IF (cos_dir(3) .GT. 0.0D0) THEN ! photon in the z-positive direction
          !--- Need to increment the layer indice
          ilayer = ilayer + 1
          IF (ilayer .GT. nlay) THEN 
             ! photon is leaving the atmosphere from the top
             leave = 2
          ELSE ! photon still in the atmosphere at layer ilayer
             !--- Compute the new position of the photon at the border between the 2 layers
             old_z_phot = z_phot
             CALL New_position(old_z_phot,cos_dir(3),dlimit,z_phot)
             GOTO 500
          ENDIF
       ELSE ! photon in the z-negative direction
          !--- Need to de-increment the layer indice
          ilayer = ilayer - 1
          IF (ilayer .LT. 1) THEN 
             ! photon leaves the atmosphere from the bottom
             leave = 1
          ELSE ! photon still in the atmosphere at layer ilayer
             !--- Compute the new position of the photon at the border between the 2 layers
             old_z_phot = z_phot
             CALL New_position(old_z_phot,cos_dir(3),dlimit,z_phot)
             GOTO 500
          ENDIF
       ENDIF

    ELSE

       !********************************
       !   PHOTON STAY IN THE LAYER
       !********************************
       !--- Keep track of photon path inside the layer
       !    (distance that the photon crossed in its direction)
       d = tau / cext_ly(ilayer)

       !--- Compute the new position of the photon
       old_z_phot = z_phot
       CALL New_position(old_z_phot,cos_dir(3),d,z_phot)

       !--- Choose which kind of target the photon is beating
       !    If index ind = 0 --> rayleight
       !                 = 1 --> particles
       IF (fracext_gas_ly(ilayer) .EQ. 1) THEN 
          ! gas 
          ind = 0 
       ELSE
          ! choose between gas and particle
          CALL target_choice(ithread, fracext_gas_ly(ilayer),ind)
       ENDIF

    ENDIF

    !********************************
    ! the surface is the target
    IF ( ( leave .eq. 1 .and. TRIM(ADJUSTL(surface)) .ne. 'none') ) then 
       ind    = 2
       z_phot = 0.0D0
       ilayer = 1   
       leave  = 0
    ENDIF

  END SUBROUTINE move_photon

  !================================================================================================

  recursive subroutine split_photon(ithread, nmat, nelem, surface, surface_albedo, tau_tot,       & 
       ddis, vrm_epsilon_ddis, vrm_wspl_crit, vrm_nsplmax, vrm_n_sccp,  vrm_n_lecp,               & 
       vrm_wrr_crit, vrm_wrr_min,                                                                 &                                    
       nmu, mu, nlay, thetad, p_ptcle_ddis, prf_ptcle_ddis, p_ptcle_ly, prf_ptcle_ly, delta_ray,  &
       z_ref,cext_ly, fracext_gas_ly, ssa_ptcle_ly, ssa_gas_ly,                                   &
       cos_dir_sensor, thetav,                                                                    &
       i_cp, Stk_in, cos_dir_in, u_in, v_in, wght_in, z_phot_in, ilayer_in, count_scat_in,        &
       leave_in, ind_in,                                                                          &
       nspl, nmax_interact, count_scat, split_level, Stoke_out)

    USE MMCRAD1D_CONSTANTS
    USE MMCRAD1D_ZIGGURAT
    USE MMCRAD1D_BRDF, only : alb_surf, tetasurf, ntetasurf
    USE MMCRAD1D_INTERPOLATE, only : MCLININTPOL

    IMPLICIT NONE

    integer, INTENT(IN) :: ithread
    integer, INTENT(IN) :: nmat
    integer, INTENT(IN) :: nelem
    CHARACTER(len=*), INTENT(IN) :: surface
    REAL(8), INTENT(IN) :: surface_albedo
    REAL(8), INTENT(IN) :: tau_tot
    logical, INTENT(IN) :: ddis
    REAL(8), INTENT(IN) :: vrm_epsilon_ddis
    REAL(8), INTENT(IN) :: vrm_wspl_crit
    integer, INTENT(IN) :: vrm_nsplmax
    integer, INTENT(IN) :: vrm_n_sccp
    integer, INTENT(IN) :: vrm_n_lecp
    REAL(8), INTENT(IN) :: vrm_wrr_crit                                    
    REAL(8), INTENT(IN) :: vrm_wrr_min                                    
    integer, INTENT(IN) :: nmu 
    real(8), INTENT(IN) :: mu(nmu)
    integer, INTENT(IN) :: nlay
    real(8), INTENT(IN) :: thetad(nmu)
    real(8), INTENT(IN) :: p_ptcle_ddis(nmu)
    real(8), INTENT(IN) :: prf_ptcle_ddis(nmu)
    real(8), INTENT(IN) :: p_ptcle_ly(nlay,nelem,nmu) 
    real(8), INTENT(IN) :: prf_ptcle_ly(nlay,nmu)
    REAL(8), INTENT(IN) :: delta_ray
    REAL(8), INTENT(IN) :: z_ref(0:nlay)
    real(8), INTENT(IN) :: cext_ly(nlay)
    real(8), INTENT(IN) :: fracext_gas_ly(nlay)
    REAL(8), INTENT(IN) :: ssa_ptcle_ly(nlay)
    REAL(8), INTENT(IN) :: ssa_gas_ly(nlay)
    real(8), INTENT(IN) :: cos_dir_sensor(3)
    real(8), INTENT(IN) :: thetav
    INTEGER, INTENT(IN) :: i_cp
    real(8), INTENT(IN) :: Stk_in(nmat)
    real(8), INTENT(IN) :: cos_dir_in(3)
    real(8), INTENT(IN) :: u_in(3)
    real(8), INTENT(IN) :: v_in(3)
    real(8), INTENT(IN) :: wght_in
    real(8), INTENT(IN) :: z_phot_in
    integer, INTENT(IN) :: ilayer_in
    integer, INTENT(IN) :: count_scat_in
    integer, INTENT(IN) :: leave_in
    integer, INTENT(IN) :: ind_in
    integer, INTENT(IN) :: nspl
    integer, INTENT(IN) :: nmax_interact
    integer, INTENT(IN) :: count_scat

    integer, INTENT(INOUT) :: split_level
    real(8), INTENT(INOUT) :: Stoke_out(nmat)

    ! --- local variables
    REAL(8) :: rn_rr
    real(8) :: le(nmat)

    integer :: ispl

    ! splitted cloned photon characteristics
    REAL(8) :: Stk_spl(nmat)
    REAL(8) :: wght_spl
    REAL(8) :: cos_dir_spl(3)
    REAL(8) :: u_spl(3)
    REAL(8) :: v_spl(3)
    REAL(8) :: z_phot_spl
    INTEGER :: ilayer_spl
    INTEGER :: leave_spl
    INTEGER :: count_scat_spl
    REAL(8) :: wspl_spl

    integer :: ind_spl

    ! -----------------

    split_level = split_level + 1 

    do ispl = 1, nspl

       Stk_spl           = stk_in
       wght_spl          = wght_in / nspl
       cos_dir_spl       = cos_dir_in
       u_spl             = u_in
       v_spl             = v_in
       z_phot_spl        = z_phot_in
       ilayer_spl        = ilayer_in
       leave_spl         = leave_in
       ind_spl           = ind_in
       count_scat_spl    = count_scat_in

       DO WHILE ((leave_spl .EQ. 0) .AND. (count_scat_spl .le. vrm_n_sccp))   

          ! compute Wspl
          CALL get_wspl(nmu, thetad, p_ptcle_ddis, &
               wght_spl, cos_dir_spl,              &      
               cos_dir_sensor,                     &
               wspl_spl)

          IF (wspl_spl .gt. vrm_wspl_crit) THEN

             CALL split_photon(ithread, nmat, nelem, surface, surface_albedo,  tau_tot,                      & 
                  ddis, vrm_epsilon_ddis, vrm_wspl_crit, vrm_nsplmax, vrm_n_sccp, vrm_n_lecp,                &
                  vrm_wrr_crit, vrm_wrr_min, &    
                  nmu, mu, nlay, thetad, p_ptcle_ddis, prf_ptcle_ddis, p_ptcle_ly, prf_ptcle_ly, delta_ray,  &
                  z_ref, cext_ly, fracext_gas_ly, ssa_ptcle_ly, ssa_gas_ly,                                  &
                  cos_dir_sensor, thetav,                                                                    &
                  i_cp,Stk_spl, cos_dir_spl, u_spl, v_spl, wght_spl, z_phot_spl, ilayer_spl, count_scat_spl, &
                  leave_spl, ind_spl, MIN(INT(wspl_spl), vrm_nsplmax),                    &
                  nmax_interact, count_scat, split_level, Stoke_out)

             leave_spl = 3

          ELSEIF ( (wspl_spl .lt. vrm_wrr_crit)       .and. & 
               ((wght_spl / MAX(wspl_spl, vrm_wrr_min)) .le. 1.0d0) ) THEN

             !***********************************
             !      Russion roulette
             !CALL RANDOM_NUMBER(rn_rr)
             rn_rr = par_ziguni(ithread)
             IF (rn_rr .le. (1.0D0 - MAX(wspl_spl, vrm_wrr_min)) ) THEN
                ! kill the photon
                !write(*,*) ' kill SPL', &
                !     wspl_spl, rn_rr, (1.0D0 - MAX(wspl_spl, vrm_wrr_min))
                leave_spl = 3
             ELSE
                ! Increase the photon weight
                wght_spl = wght_spl / MAX(wspl_spl, vrm_wrr_min)
                ! Propagate the photon to its new location
                ! determine the kind of target that stops the photon
                CALL MOVE_PHOTON(ithread,   surface,                 &
                     nlay, z_ref, cext_ly, fracext_gas_ly,           &
                     cos_dir_spl, ilayer_spl, z_phot_spl, leave_spl, &
                     ind_spl)
             ENDIF
             !***********************************
          ELSE

             ! Propagate the photon to its new location
             ! determine the kind of target that stops it
             CALL MOVE_PHOTON(ithread,  surface,                  &
                  nlay, z_ref, cext_ly, fracext_gas_ly,           &
                  cos_dir_spl, ilayer_spl, z_phot_spl, leave_spl, &
                  ind_spl)

          ENDIF

          IF (leave_spl .eq. 0) THEN  ! photon  interacts

             count_scat_spl = count_scat_spl + 1

             !---  Modifiy the photon weight to account for absorption
             IF (ind_spl .EQ. 0) THEN
                wght_spl = wght_spl * ssa_gas_ly(ilayer_spl)                  
             ELSE IF (ind_spl .EQ. 1) THEN
                wght_spl = wght_spl * ssa_ptcle_ly(ilayer_spl)
             ELSE IF ((ind_spl .EQ. 2) .AND. (TRIM(ADJUSTL(surface)) .eq. 'lambert')) then
                wght_spl = wght_spl * surface_albedo   
             ELSE IF ((ind_spl .EQ. 2) .AND. (TRIM(ADJUSTL(surface)) .eq. 'brdf')) then
                ! we need to interpolate to get the value of the albedo for the incident theta on the BRDF
                wght_spl = wght_spl * MCLININTPOL(alb_surf, tetasurf, ntetasurf,  ACOS(-(cos_dir_spl(3)))/DTR)
             ENDIF

             ! LE on SPL
             IF ( i_cp .eq. 1 .or. &
                  ((vrm_n_sccp+1-count_scat_spl) .lt. vrm_n_lecp) ) THEN
                ! for the first CP of the MP, we must perform the LE for each scattering
                ! to complet the Von-Neumann serie
                ! for the following CP, we perform the LE for the last vrm_n_lecp scattering only
                CALL local_estimate(nmat, nelem, surface, nlay, cext_ly, z_ref, tau_tot,  &
                     ind_spl, nmu, thetad, p_ptcle_ly(ilayer_spl,:,:), delta_ray, &
                     wght_spl, Stk_spl, z_phot_spl, cos_dir_spl, v_spl, ilayer_spl,   & 
                     thetav, cos_dir_sensor(:), le)
                Stoke_out(:) = Stoke_out(:) + le(:)
                !write(*,*) '   SCP ', count_scat, count_scat_spl, i_cp, ispl, split_level, le(1)
             ENDIF

             IF (count_scat_spl .le. vrm_n_sccp) then
                ! compute new stoke vector, new cos_dir, u and v after scattering
                CALL scatter(ithread, nmat, nelem,  surface,                                                      &  
                     ddis, vrm_epsilon_ddis, p_ptcle_ddis, prf_ptcle_ddis,                                        &
                     cos_dir_sensor,                                                                              &
                     ind_spl, nmu, mu, thetad, p_ptcle_ly(ilayer_spl,:,:), prf_ptcle_ly(ilayer_spl,:), delta_ray, &
                     wght_spl, Stk_spl, v_spl, u_spl, cos_dir_spl)
             ENDIF

          ENDIF !condition on photon interacting

          !if (leave_spl .ne. 0) write(*,*) '   SP leave=', leave_spl

       ENDDO ! stop following the photon

    ENDDO ! loop on splitted cloned photons

    split_level = split_level - 1 

  end subroutine split_photon

END MODULE mmcrad1d_submain
