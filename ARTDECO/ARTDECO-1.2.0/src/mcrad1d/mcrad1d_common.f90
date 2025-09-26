
MODULE MMCRAD1D_BRDF

  implicit none

  ! surface BRDF related variables
  integer, PUBLIC :: nphisurf 
  integer, PUBLIC :: ntetasurf
  real(8), PUBLIC :: dphisurf
  real(8), PUBLIC :: dtetasurf
  real(8), PUBLIC, allocatable :: phisurf(:)
  real(8), PUBLIC, allocatable :: tetasurf(:)
  real(8), PUBLIC, allocatable :: brdf(:,:,:,:,:)
  real(8), PUBLIC, allocatable :: prf_phi_surf(:,:,:)
  real(8), PUBLIC, allocatable :: prf_mu_surf(:,:)
  real(8), PUBLIC, allocatable :: alb_surf(:)

END MODULE MMCRAD1D_BRDF
