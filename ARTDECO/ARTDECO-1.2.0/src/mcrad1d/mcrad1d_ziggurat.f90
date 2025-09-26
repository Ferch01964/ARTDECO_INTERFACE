! Marsaglia & Tsang generator for random numbers
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001

! Parallel version - October 2006

! This version has been customised for parallel processing use,
! specifically with OpenMP.  Each thread uses its own pseudo-random 
! sequence. (Gib Bogle)
!--------------------------------------------------------------------------

MODULE mmcrad1d_ziggurat

  IMPLICIT NONE

  PRIVATE

  INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
  REAL(DP), PARAMETER  ::  m1=2147483648.0_DP 
  REAL(DP), PARAMETER  ::  m2=2147483648.0_DP
  REAL(DP), PARAMETER  ::  half=0.5_DP
  REAL(DP)             ::  dn0=3.442619855899_DP
  REAL(DP)             ::  tn0=3.442619855899_DP
  REAL(DP)             ::  vn=0.00991256303526217_DP
  REAL(DP)             ::  q
  REAL(DP)             ::  de0=7.697117470131487_DP
  REAL(DP)             ::  te0=7.697117470131487_DP
  REAL(DP)             ::  ve=0.003949659822581572_DP
  !   INTEGER,  SAVE       ::  iz, jz, jsr=123456789, kn(0:127),              &
  !                            ke(0:255), hz
  !   REAL(DP), SAVE       ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)
  !   LOGICAL,  SAVE       ::  initialized=.FALSE.

  integer, save :: par_n = 0
  integer, save :: par_step
  integer, allocatable, save ::  par_jsr(:)
  integer, allocatable, save ::  par_kn(:,:)
  integer, allocatable, save ::  par_ke(:,:)
  real(DP), allocatable, save :: par_wn(:,:)
  real(DP), allocatable, save :: par_fn(:,:)
  real(DP), allocatable, save :: par_we(:,:)
  real(DP), allocatable, save :: par_fe(:,:)

  PUBLIC  :: par_zigset, par_shr3, par_ziguni, par_zigunset

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE par_zigunset()

    IMPLICIT NONE

    deallocate(par_jsr)
    deallocate(par_kn)
    deallocate(par_ke)
    deallocate(par_wn)
    deallocate(par_fn)
    deallocate(par_we)
    deallocate(par_fe)

  END SUBROUTINE par_zigunset

!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE par_zigset( npar, par_jsrseed, grainsize)

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: npar, grainsize, par_jsrseed(0:npar)

    INTEGER  :: i, kpar
    REAL(DP) dn, tn, de, te

    par_n = npar
    par_step = grainsize
    ! First we need to allocate all the non-volatile arrays with the size npar
    allocate(par_jsr(0:npar*par_step))
    allocate(par_kn(0:127,0:npar-1))
    allocate(par_ke(0:255,0:npar-1))
    allocate(par_wn(0:127,0:npar-1))
    allocate(par_fn(0:127,0:npar-1))
    allocate(par_we(0:255,0:npar-1))
    allocate(par_fe(0:255,0:npar-1))

    ! Now treat each instance separately
    do kpar = 0,npar-1
       !  Set the seed
       par_jsr(kpar*par_step) = par_jsrseed(kpar)

       !  Tables for RNOR
       dn = dn0
       tn = tn0
       q = vn*EXP(half*dn*dn)
       par_kn(0,kpar) = (dn/q)*m1
       par_kn(1,kpar) = 0
       par_wn(0,kpar) = q/m1
       par_wn(127,kpar) = dn/m1
       par_fn(0,kpar) = 1.0_DP
       par_fn(127,kpar) = EXP( -half*dn*dn )
       DO  i = 126, 1, -1
          dn = SQRT( -2.0_DP * LOG( vn/dn + EXP( -half*dn*dn ) ) )
          par_kn(i+1,kpar) = (dn/tn)*m1
          tn = dn
          par_fn(i,kpar) = EXP(-half*dn*dn)
          par_wn(i,kpar) = dn/m1
       END DO

       !  Tables for REXP
       de = de0
       te = te0
       q = ve*EXP( de )
       par_ke(0,kpar) = (de/q)*m2
       par_ke(1,kpar) = 0
       par_we(0,kpar) = q/m2
       par_we(255,kpar) = de/m2
       par_fe(0,kpar) = 1.0_DP
       par_fe(255,kpar) = EXP( -de )
       DO  i = 254, 1, -1
          de = -LOG( ve/de + EXP( -de ) )
          par_ke(i+1,kpar) = m2 * (de/te)
          te = de
          par_fe(i,kpar) = EXP( -de )
          par_we(i,kpar) = de/m2
       END DO
    enddo
    RETURN
  END SUBROUTINE par_zigset

!!!!!!!!!!!!!!!!!!!!!!!!!

  !  Generate random 32-bit integers
  FUNCTION par_shr3(kpar) RESULT( ival )
    IMPLICIT NONE
    INTEGER  ::  ival, kpar
    integer :: jz, jsr

    jsr = par_jsr(kpar*par_step)
    jz = jsr
    jsr = IEOR( jsr, ISHFT( jsr,  13 ) )
    jsr = IEOR( jsr, ISHFT( jsr, -17 ) )
    jsr = IEOR( jsr, ISHFT( jsr,   5 ) )
    par_jsr(kpar*par_step) = jsr
    ival = jz + jsr
    RETURN
  END FUNCTION par_shr3

!!!!!!!!!!!!!!!!!!!!!!!!!

  !  Generate uniformly distributed random numbers, sequence kpar
  FUNCTION par_ziguni(kpar) RESULT( fn_val )
    IMPLICIT NONE
    integer :: kpar
    REAL(DP)  ::  fn_val

    if (kpar >= par_n) then
       write(*,*) 'thread number exceeds initialized max: ',kpar,par_n-1
       stop
    endif
    fn_val = half + 0.2328306e-9_DP * par_shr3(kpar)
    RETURN
  END FUNCTION par_ziguni

END MODULE mmcrad1d_ziggurat


