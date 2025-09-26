
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! RCS version control information:
! $Header: BDREF.f,v 2.1 2000/03/27 21:40:51 laszlo Exp $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REAL FUNCTION  BDREF( WVNMLO, &
     WVNMHI, MU, MUP, DPHI )

  USE MSURFACE_BRDF

  !      Supplies surface bi-directional reflectivity.
  !
  !      NOTE 1: Bidirectional reflectivity in DISORT is defined
  !              by Eq. 39 in STWL.
  !      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
  !              angles) are positive.
  !
  !  INPUT:
  !
  !    WVNMLO : Lower wavenumber (inv cm) of spectral interval
  !
  !    WVNMHI : Upper wavenumber (inv cm) of spectral interval
  !
  !    MU     : Cosine of angle of reflection (positive)
  !
  !    MUP    : Cosine of angle of incidence (positive)
  !
  !    DPHI   : Difference of azimuth angles of incidence and reflection
  !                (radians)  
  !   Called by- DREF, SURFAC
  ! +-------------------------------------------------------------------+

  IMPLICIT NONE

  !     .. Scalar Arguments ..
  REAL      wl
  REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
  REAL      brdf(1,1)

  REAL dtr, pi 
  PARAMETER (pi  = 4.0*atan(1.0))
  PARAMETER (dtr = pi/180.0)

  ! get a mean wavelength
  !wl  = 2.0 / (WVNMHI+WVNMLO) * 1e4
  wl = ((1.0 / (WVNMHI))+ (1.0 / (WVNMLO))) / 2.0 * 1e4


  CALL surface_brdf(1, wl, acos(mup)/dtr, &
       acos(mu)/dtr, dphi/dtr+180.0, brdf)

  BDREF = brdf(1,1)

  RETURN

END FUNCTION BDREF
