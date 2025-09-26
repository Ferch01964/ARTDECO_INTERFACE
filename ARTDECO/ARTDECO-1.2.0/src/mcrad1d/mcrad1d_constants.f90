
MODULE MMCRAD1D_CONSTANTS

  implicit none

  real(8), public, parameter                      :: DTR   = 0.0174532925199432957D0
  real(8), public, parameter                      :: EPS   = 1.0D-06
  real(8), public, parameter                      :: EPS12 = 1.0D-12
  real(8), public, parameter                      :: EPS1  = 1.0D-1

  real(8), public, parameter                      :: Pi  = 3.1415926535897932385D0

  !--- interpolation flag               
  integer, public, parameter                      :: flag_interpol = 0
  !                                                  = 0 -> linear interpolation
  !                                                  = 1 -> polynomial

  real(8), public, parameter                      :: min_wght = 1.0d-100
  !                                               we stop following photons whose weight is less than min_weight

END MODULE MMCRAD1D_CONSTANTS
