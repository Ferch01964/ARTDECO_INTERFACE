
REAL(8) FUNCTION DOAD_PLKAVG( WNUMLO, WNUMHI, T )

  ! ===================================================================
  !       This routine was taken from DISORT 2.
  ! ===================================================================
  !        Computes Planck function integrated between two wavenumbers

  !  INPUT :  WNUMLO : Lower wavenumber (inv cm) of spectral interval

  !           WNUMHI : Upper wavenumber

  !           T      : Temperature (K)

  !  OUTPUT : DOAD_PLKAVG : Integrated Planck function ( Watts/sq m )
  !                      = Integral (WNUMLO to WNUMHI) of
  !                        2h c**2  nu**3 / ( EXP(hc nu/kT) - 1)
  !                        (where h=Plancks constant, c=speed of
  !                         light, nu=wavenumber, T=temperature,
  !                         and k = Boltzmann constant)

  !  Reference : Specifications of the Physical World: New Value
  !                 of the Fundamental Constants, Dimensions/N.B.S.,
  !                 Jan. 1974

  !  Method :  For WNUMLO close to WNUMHI, a Simpson-rule quadrature
  !            is done to avoid ill-conditioning; otherwise

  !            (1)  For WNUMLO or WNUMHI small,
  !                 integral(0 to WNUMLO/HI) is calculated by expanding
  !                 the integrand in a power series and integrating
  !                 term by term;

  !            (2)  Otherwise, integral(WNUMLO/HI to INFINITY) is
  !                 calculated by expanding the denominator of the
  !                 integrand in powers of the exponential and
  !                 integrating term by term.

  !  Accuracy :  At least 6 significant digits, assuming the
  !              physical constants are infinitely accurate

  !  ERRORS WHICH ARE NOT TRAPPED:

  !      * power or exponential series may underflow, giving no
  !        significant digits.  This may or may not be of concern,
  !        depending on the application.

  !      * Simpson-rule special case is skipped when denominator of
  !        integrand will cause overflow.  In that case the normal
  !        procedure is used, which may be inaccurate if the
  !        wavenumber limits (WNUMLO, WNUMHI) are close together.

  !  LOCAL VARIABLES

  !        A1,2,... :  Power series coefficients
  !        C2       :  h * c / k, in units cm*K (h = Plancks constant,
  !                      c = speed of light, k = Boltzmann constant)
  !        D(I)     :  Exponential series expansion of integral of
  !                       Planck function from WNUMLO (i=1) or WNUMHI
  !                       (i=2) to infinity
  !        EPSIL    :  Smallest number such that 1+EPSIL .GT. 1 on
  !                       computer
  !        EX       :  EXP( - V(I) )
  !        EXM      :  EX**M
  !        MMAX     :  No. of terms to take in exponential series
  !        MV       :  Multiples of V(I)
  !        P(I)     :  Power series expansion of integral of
  !                       Planck function from zero to WNUMLO (I=1) or
  !                       WNUMHI (I=2)
  !        PI       :  3.14159...
  !        SIGMA    :  Stefan-Boltzmann constant (W/m**2/K**4)
  !        SIGDPI   :  SIGMA / PI
  !        SMALLV   :  Number of times the power series is used (0,1,2)
  !        V(I)     :  C2 * (WNUMLO(I=1) or WNUMHI(I=2)) / temperature
  !        VCUT     :  Power-series cutoff point
  !        VCP      :  Exponential series cutoff points
  !        VMAX     :  Largest allowable argument of EXP function

  ! ----------------------------------------------------------------------

  implicit none

  REAL(8), INTENT(IN) :: T
  REAL(8), INTENT(IN) :: WNUMHI
  REAL(8), INTENT(IN) :: WNUMLO

  !     .. Parameters ..
  REAL(8), PARAMETER :: A1 =  1. / 3.
  REAL(8), PARAMETER :: A2 = -1. / 8.
  REAL(8), PARAMETER :: A3 =  1. / 60.
  REAL(8), PARAMETER :: A4 = -1. / 5040.
  REAL(8), PARAMETER :: A5 =  1. / 272160.
  REAL(8), PARAMETER :: A6 = -1. / 13305600. 

  !     .. Local Scalars ..

  INTEGER :: I, K, M, MMAX, N, SMALLV
  REAL(8) :: C2, CONC, DEL, EPSIL, EX, EXM, HH, MV, OLDVAL, PI
  REAL(8) :: SIGDPI, SIGMA, VAL, VAL0, VCUT, VMAX, VSQ, X

  !     .. Local Arrays ..
  REAL(8)  :: D( 2 ), P( 2 ), V( 2 ), VCP( 7 )

  !     .. Statement Functions ..
  REAL(8)   ::    PLKF
  !     .. Statement Function definitions ..
  PLKF( X ) = X**3 / ( EXP( X ) - 1 )
  ! ---


  C2    = 1.438786  
  SIGMA = 5.67032D-8  
  VCUT  = 1.5 
  VCP   = (/ 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /)

  PI     = 2.*ASIN( 1.0 )
  VMAX   = LOG( HUGE(1.0D0) ) 
  EPSIL  = TINY(1.0D0)!RADIX(1.0D0) ** (1-DIGITS(1.0D0))
  SIGDPI = SIGMA / PI
  CONC   = 15. / PI**4

  IF( T.LT.0.0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LT.0. ) THEN
     WRITE(*,*) 'DOAD_PLKAVG--temperature or wavenums. wrong'
     WRITE(*,*) ' T=',T
     WRITE(*,*) ' WNUMLO=',WNUMLO
     WRITE(*,*) ' WNUMHI=',WNUMHI
     STOP
  ENDIF

  IF( T .LT. 1.E-4 ) THEN

     DOAD_PLKAVG = 0.0
     RETURN

  END IF

  V( 1 ) = C2*WNUMLO / T
  V( 2 ) = C2*WNUMHI / T

  IF( V( 1 ).GT.EPSIL .AND. V( 2 ).LT.VMAX .AND. &
       ( WNUMHI - WNUMLO ) / WNUMHI .LT. 1.E-2 ) THEN

     ! ** Wavenumbers are very close.  Get integral
     ! ** by iterating Simpson rule to convergence.

     HH     = V( 2 ) - V( 1 )
     OLDVAL = 0.0
     VAL0   = PLKF( V( 1 ) ) + PLKF( V( 2 ) )

     DO N = 1, 10

        DEL  = HH / ( 2*N )
        VAL  = VAL0

        DO K = 1, 2*N - 1
           VAL  = VAL + 2*( 1 + MOD( K,2 ) )*  &
                PLKF( V( 1 ) + K*DEL )
        ENDDO

        VAL  = DEL / 3.*VAL
        IF( ABS( ( VAL - OLDVAL ) / VAL ).LE.1.E-6 ) GO TO  30
        OLDVAL = VAL

     END DO

     WRITE(*,*) 'DOAD_PLKAVG--Simpson rule didnt converge'

30   CONTINUE

     DOAD_PLKAVG = SIGDPI * T**4 * CONC * VAL

     RETURN

  END IF

  ! *** General case ***
  SMALLV = 0

  DO I = 1, 2

     IF( V( I ).LT.VCUT ) THEN
        ! ** Use power series
        SMALLV = SMALLV + 1
        VSQ    = V( I )**2
        P( I ) = CONC*VSQ*V( I )*( A1 + &
             V( I )*( A2 + V( I )*( A3 + VSQ*( A4 + VSQ*( A5 + &
             VSQ*A6 ) ) ) ) )

     ELSE
        ! ** Use exponential series
        MMAX  = 0
        ! ** Find upper limit of series
40      CONTINUE
        MMAX  = MMAX + 1

        IF( V(I) .LT. VCP( MMAX ) ) GO TO  40

        EX     = EXP( - V(I) )
        EXM    = 1.0
        D( I ) = 0.0

        DO M = 1, MMAX
           MV     = M*V( I )
           EXM    = EX*EXM
           D( I ) = D( I ) + EXM*( 6.+ MV*( 6.+ MV*( 3.+ MV ) ) ) / M**4
        ENDDO
        D( I ) = CONC*D( I )

     END IF

  END DO

  DOAD_PLKAVG= -32768.0

  ! ** Handle ill-conditioning
  IF( SMALLV.EQ.2 ) THEN
     ! ** WNUMLO and WNUMHI both small
     DOAD_PLKAVG = P( 2 ) - P( 1 )
  ELSE IF( SMALLV.EQ.1 ) THEN
     ! ** WNUMLO small, WNUMHI large
     DOAD_PLKAVG = 1.- P( 1 ) - D( 2 )
  ELSE
     ! ** WNUMLO and WNUMHI both large
     DOAD_PLKAVG = D( 1 ) - D( 2 )
  END IF

  DOAD_PLKAVG = SIGDPI * T**4 * DOAD_PLKAVG

  IF( DOAD_PLKAVG.EQ.0.0 ) &
       WRITE(*,*) 'DOAD_PLKAVG--returns zero; possible underflow'

  RETURN

END FUNCTION DOAD_PLKAVG


