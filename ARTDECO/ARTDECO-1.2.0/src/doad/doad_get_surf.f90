

subroutine doad_get_surf(nmat, lambda, nmug, xmug, NDcoef, refl)

  ! Mathieu Compiegne (inspired from micrad.F from Michele Vesperini)  

  USE MSURFACE_BRDF

  IMPLICIT NONE

  INTEGER, INTENT(IN)       :: nmat
  REAL(kind=8), INTENT(IN)  :: lambda
  INTEGER, INTENT(IN)       :: nmug
  REAL(kind=8), INTENT(IN)  :: xmug(nmug)
  INTEGER, INTENT(IN)       :: NDcoef
  REAL(kind=8), INTENT(OUT) :: refl(nmat,nmat,nmug,nmug,0:NDcoef)

  ! local variables
  real(kind=8), parameter :: pi = 3.1415926535897932384D0
  real(kind=8), parameter :: deg2rad = pi/180.0D0

  INTEGER :: nphi
  real(kind=8), allocatable :: bdrfphi(:,:,:)
  real(kind=8), allocatable :: phi(:)
  real(kind=8), allocatable :: a(:), b(:)
  real(kind=8) :: dphi
  INTEGER :: iphi
  INTEGER :: imu
  INTEGER :: jmu
  INTEGER :: im
  INTEGER :: imat
  INTEGER :: jmat
  INTEGER :: imsurf

  INTEGER :: m
  INTEGER :: kmat

  REAL(kind=8) :: ts
  REAL(kind=8) :: tv

  !------------------

  refl(:,:,:,:,:) = 0.0D0

  dphi = 2.0d0 ! choix arbitraire
  nphi=nint(360.0D0/dphi)+1
  allocate(phi(nphi))

  ! order of the highest Fourier term
  m = min(nphi/2,ndcoef)! a verifier freq niquist
  ALLOCATE(a(0:m), b(0:m))

  ALLOCATE(bdrfphi(nmat,nmat,nphi))

  do imu = nmug, 1, -1
     ! reflected angle
     tv = acos(xmug(imu))/deg2rad

     do jmu = nmug, 1, -1
        ! incident angle
        ts = acos(xmug(jmu))/deg2rad

        ! compute the BDRF for all phi angles
        do iphi=1,nphi
           phi(iphi)=(iphi-1)*dphi
           CALL surface_brdf(nmat, lambda, ts, tv, 180.0D0+phi(iphi), bdrfphi(:,:,iphi))
        enddo

        !-----------------------------------------
        ! -- decomposition en serie de Fourier
        !-----------------------------------------
        ! --  the phase matrix 
        ! --  equivalence for kmat
        !        | z11 z12 z13 z14 |			| Z1  Z2  Z3  Z4 |
        !        | z21 z22 z23 z24 |			| Z5  Z6  Z7  Z8 |
        !        | z31 z32 z33 z34 |	equivalent to 	| Z9  Z10 Z11 Z12|
        !        | z41 z42 z43 z44 |			| Z13 Z14 Z15 Z16|

        do imat=1,nmat
           do jmat=1,nmat
              kmat=(imat-1)*4+jmat
              ! -- terms 4 and 13 are nul (f_to_matrix)
              if ((kmat.eq.4).or.(kmat.eq.13)) then
                 do im=0,m
                    refl(imat,jmat,imu,jmu,im) = 0.0D0
                 enddo
                 ! -- for doad: odd  terms of Stokes matrix are developped in sine (coeff b)
              elseif ((kmat.eq.3).or.(kmat.eq.7).or.(kmat.eq.8)) then
                 call doad_foursin(m, nphi, phi(:)*deg2rad, bdrfphi(imat,jmat,:), b)
                 do im=0,m
                    ! -- les termes en haut a droite sont negatifs
                    refl(imat,jmat,imu,jmu,im) = -b(im)
                 enddo
              elseif ((kmat.eq.9).or.(kmat.eq.10).or.(kmat.eq.14)) then
                 call doad_foursin(m, nphi, phi(:)*deg2rad, bdrfphi(imat,jmat,:), b)
                 do im=0,m
                    refl(imat,jmat,imu,jmu,im) = b(im)
                 enddo
              else 
                 ! -- for doad: even terms of Stokes matrix are developped in cosine (coeff a)
                 call doad_fourcos(m, nphi, phi(:)*deg2rad, bdrfphi(imat,jmat,:), a)
                 do im = 0, m
                    refl(imat,jmat,imu,jmu,im) = a(im)
                 enddo
              endif
           enddo
        enddo

     enddo ! loop on incident mug
  enddo ! looop on reflected mug

 
  ! -- test 
  imsurf=0
  do jmu=1,nmug
     do imu=1,nmug
        do im=0,ndcoef
           if (refl(1,1,imu,jmu,im).gt.1d-6) imsurf=max(im,imsurf)
        enddo
     enddo
  enddo
  !write(*,*) ' '
  !write(*,*) "        terme de fourier max pour la surface ",imsurf

  DEALLOCATE(bdrfphi)
  deallocate(phi)
  deallocate(a, b)

end subroutine doad_get_surf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine doad_fourcos(nfou, nmu, x, f, res)

  ! see Hansen, Journal of Atmospheric Sciences, 1971, 28, 120-125
  ! eq. 27

  IMPLICIT NONE

  integer, INTENT(IN) :: nfou
  integer, INTENT(IN) :: nmu
  real(kind=8), INTENT(IN)  :: x(nmu) ! the angle given in rad
  real(kind=8), INTENT(IN)  :: f(nmu)
  real(kind=8), INTENT(OUT) :: res(0:nfou)
  real(kind=8), parameter :: pi = 3.1415926535897932384D0

  ! local arguments
  integer :: m
  real(kind=8) :: sum

  do m = 0, nfou
     call doad_xinteg2(1, nmu, nmu, x, f*cos(m*x), sum)
     res(m)  =  1.0D0 / 2.0D0 / pi * sum
  enddo

end subroutine doad_fourcos


subroutine doad_foursin(nfou, nmu, x, f, res)

  ! see Hansen, Journal of Atmospheric Sciences, 1971, 28, 120-125
  ! eq. 28

  IMPLICIT NONE

  integer, INTENT(IN) :: nfou
  integer, INTENT(IN) :: nmu
  real(kind=8), INTENT(IN)  :: x(nmu) ! the angle given in rad
  real(kind=8), INTENT(IN)  :: f(nmu)
  real(kind=8), INTENT(OUT) :: res(0:nfou)
  real(kind=8), parameter :: pi = 3.1415926535897932384D0

  ! local arguments
  integer :: m
  real(kind=8) :: sum

  m       = 0
  res(m)  = 0.0D0

  ! m > 0
  do m = 1, nfou
     call doad_xinteg2(1, nmu, nmu, x, f*sin(m*x), sum)
     res(m)  =  1.0D0 / 2.0D0 / pi * sum
  enddo

end subroutine doad_foursin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE doad_xinteg2(imin, imax, n, xin, yin, res)
  ! computes integral of yin, variable step 
  ! make sure you have at least two points of integration
  ! Trapeze method

  IMPLICIT NONE

  INTEGER, INTENT (IN)        :: imin
  INTEGER, INTENT (IN)        :: imax
  INTEGER, INTENT (IN)        :: n
  REAL (KIND=8), INTENT (IN)  :: xin(n)
  REAL (KIND=8), INTENT (IN)  :: yin(n)
  REAL (KIND=8), INTENT (OUT) ::  res
  REAL (KIND=8)               :: xa, ya, xb, yb
  REAL (KIND=8)               :: primitaux
  INTEGER                     :: i

  xa = xin(imin)
  ya = yin(imin)
  primitaux = 0.0D0
  DO i=imin+1,imax
     if (ISNAN(yin(i)) .eqv. .FALSE.) then
        xb = xin(i)
        yb = yin(i)
        primitaux = primitaux + (xb - xa) * (ya + yb)
        xa = xb
        ya = yb
     ! else

     !    write(*,*) '(doad_xinteg2) : WARNING:'
     !    write(*,*) '                 A NaN values was encountered in the integrated function'
     !    write(*,*) '                 Compuation was done avoiding this value'

     endif
  ENDDO
  res = primitaux * 0.5D0

  !CONTAINS
  !  logical function isnan(x)
  !    real(KIND=8) :: x
  !    isnan = x .ne. x
  !    return
 !  end function isnan
END SUBROUTINE doad_xinteg2
