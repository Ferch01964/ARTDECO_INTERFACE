
MODULE MSOL

  IMPLICIT NONE

  PUBLIC  :: varsol, airmass, possol
  PRIVATE :: pos_fft

  REAL(KIND=8), PRIVATE, PARAMETER :: xpi     = 3.1415926535897932384D0
  REAL(KIND=8), PRIVATE, PARAMETER :: deg2rad = xpi/180.0D0

CONTAINS

  Subroutine airmass(zenith, mass)

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !    Get airmass from solar zenith angle
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN)  :: zenith
    REAL(kind=8), INTENT(OUT) :: mass
    real(kind=8) :: a
    real(kind=8) :: b
    real(kind=8) :: c
    real(kind=8) :: z
    real(kind=8) :: el
    !----
    a  = 0.50572D0
    b  = 6.07995D0
    c  = 1.6364D0
    el = 90.0D0-zenith
    z  = xpi*el/180D0
    mass = 1.0D0 /(sin(z) + a * (el+ b)**(-c) )

  end subroutine airmass

  ! ===========================================================

  subroutine varsol (jday, dsol)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: jday
    real(kind=8), INTENT(OUT) :: dsol

    real(kind=8) :: om

    ! ------

    om   = (0.9856D0*float(jday-4))*xpi/180.D0
    dsol = 1.0D0/((1.0D0-0.01673D0*cos(om))**2.0D0)

  end subroutine varsol

  ! ===========================================================

  subroutine possol (jday, tu, xlon, xlat, teta, phi)

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !  Subroutine from 6S package to get Sun position
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! solar position (zenithal angle teta, azimuthal angle phi
    ! in degrees)
    ! jday is the day of the year

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: jday
    REAL(kind=8), INTENT(IN) :: tu
    REAL(kind=8), INTENT(IN) :: xlon
    REAL(kind=8), INTENT(IN) :: xlat
    REAL(kind=8), INTENT(OUT) :: teta
    REAL(kind=8), INTENT(OUT) :: phi

    !-----

    call pos_fft (jday, tu, xlon, xlat, teta, phi)

  end subroutine possol

  ! ===========================================================

  subroutine pos_fft (j,tu,xlon,xlat,teta,phi)

    ! solar position (zenithal angle teta,azimuthal angle phi
    ! in degrees)
    ! j is the day number in the year
    !
    ! mean solar time (heure decimale)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: j
    REAL(kind=8), INTENT(IN) :: tu
    REAL(kind=8), INTENT(IN) :: xlon
    REAL(kind=8), INTENT(IN) :: xlat
    REAL(kind=8), INTENT(OUT) :: teta
    REAL(kind=8), INTENT(OUT) :: phi

    real(kind=8) :: tsm, xla, xj, tet
    real(kind=8) :: a1, a2, a3, a4, a5, et, tsv, ah, b1, b2, b3, b4
    real(kind=8) :: b5, b6, b7, delta, amuzero, elev, az, caz, azim

    ! -------

    tsm = tu+xlon/15.0D0
    xla = xlat*deg2rad
    xj  = float(j)
    tet = 2.0D0*xpi*xj/365.D0

    ! time equation (in mn.dec)
    a1 = 0.000075D0
    a2 = 0.001868D0
    a3 = 0.032077D0
    a4 = 0.014615D0
    a5 = 0.040849D0
    et = a1+a2*cos(tet)-a3*sin(tet)-a4*cos(2.*tet)-a5*sin(2.*tet)
    et = et*12.D0*60.D0/xpi

    ! true solar time
    tsv = tsm+et/60.D0
    tsv = (tsv-12.D0)

    ! hour angle
    ah = tsv*15.D0*deg2rad

    ! solar declination (in radian)
    b1 = 0.006918D0
    b2 = 0.399912D0
    b3 = 0.070257D0
    b4 = 0.006758D0
    b5 = 0.000907D0
    b6 = 0.002697D0
    b7 = 0.001480D0
    delta = b1-b2*cos(tet)+b3*sin(tet)-b4*cos(2.D0*tet)+b5*sin(2.D0*tet)-&
         b6*cos(3.D0*tet)+b7*sin(3.D0*tet)

    ! elevation,azimuth
    amuzero = sin(xla)*sin(delta)+cos(xla)*cos(delta)*cos(ah)
    elev    = asin(amuzero)
    az      = cos(delta)*sin(ah)/cos(elev)
    if ( (abs(az)-1.000D0).gt.0.00000D0) az = sign(1.D0,az)
    caz = (-cos(xla)*sin(delta)+sin(xla)*cos(delta)*cos(ah))/cos(elev)
    azim = asin(az)
    if(caz.le.0.D0) azim = xpi-azim
    if((caz.gt.0.D0).and.(az.le.0D0)) azim = 2D0*xpi+azim
    azim = azim+xpi
    if(azim.gt.(2D0*xpi)) azim = azim-(2D0*xpi)
    elev = elev*180.D0/xpi
    ! conversion in degrees
    teta = 90.D0 - elev
    phi  = azim/deg2rad

  end subroutine pos_fft


  ! subroutine day_number(jday,month,ia,j)

  !   integer jday, month, ia, j

  !   if (month.le.2) then
  !      j=31*(month-1)+jday
  !      return
  !   endif
  !   if (month.gt.8) then
  !      j=31*(month-1)-((month-2)/2)-2+jday
  !   else
  !      j=31*(month-1)-((month-1)/2)-2+jday
  !   endif
  !   if(ia.ne.0 .and. mod(ia,4).eq.0) j=j+1
  !   return

  ! end subroutine day_number

END MODULE MSOL
