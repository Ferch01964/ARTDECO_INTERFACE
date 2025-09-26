
subroutine doad(verbose,                                   &
     NDmu,NDmuext,NDsup,NDcoef,NDmu0,NDlay,NDirf,NDphiext, &
     lambda,                                               &
     nlayer,nfoumx,nmug_in,nmat,igrnd_in, ncoef_min_bdrf,  &
     Svin, eps,                                            &
     nmuext_in,xmuext_in,nphout,phout,nmu0,xmu0,           &
     b_in, ssa, alfbet_in, ncoef,                          &
     surf_alb,                                             &
     irf, tlev, tsurf,                                     &
     wvnmlo, wvnmhi,                                       &
     flux,                                                 &
     SvR,                                                  &
     SvR_DW,                                               &
     flux_th,                                              &
     SvR_th )

  !-----------------------------------------------------------------------
  ! See : 
  !      de Haan, J. F.; Bosma, P. B.; Hovenier, J. W,
  !      The adding method for multiple scattering calculations 
  !      of polarized ligh
  !      Publication Date:	09/1987
  ! -----------------------------------------------------------------------
  !      
  !   Author:      Michele Vesperini
  !   Version:     96 04 12
  !
  !   Mathieu Compiegne, fall 2012 :
  !   - input/output modified
  !   - add some code from micrad.F 
  !   - no more implicit declaration 
  !   - all array dimensions are passing through arguments   
  !      (the static declaration using parameters was very problematic)
  !   - translation to f90
  !   - arguments are INTENT(IN or OUT) (more secure)
  !   - add NDmuext intead of a dimension of 11
  !
  !  Adding and doubling RT code including internal source and polarization
  !  MILADI change on 2015 :
  !  ADD of downward intensities
  ! -----------------------------------------------------------------------
  !
  !  ===== INPUT :
  !
  !   verbose : logical
  !
  !   Array dimensions
  !           NDmu     : integer
  !           NDmuext  : integer
  !           NDsup    : integer
  !           NDcoef   : integer
  !           NDmu0    : integer
  !           NDlay    : integer
  !           NDirf    : integer
  !           NDphiext : integer
  ! 
  !   lambda : real, wavelength (in microns)
  !
  !   nlayer : int, number of atmospheric layers
  !
  !   nfoumx : int, maximum number of Fourier term to use
  !
  !   nmug_in   : int, number of gauss points for the mu-integration
  !
  !   nmat   : int, number of elements of the Stokes vector taken into     
  !                account (4 = full polarization, 3 = 3x3 approximation, 
  !                2 = illegal, 1 = scalar)                               
  !
  !   igrnd_in   : int, index for ground model  
  !             = 0 : non surf (in fact lambert with alb=0.0)
  !             = 1 : Lambert surfaces          ) 0-th term only   
  !             = 2 : Lommel-Seeliger surfaces  )                          
  !                  see Van de Hulst, 1980 section 18.1.4 page 605                 
  !             = 3 : bidirectional reflectance )     
  !
  !  ncoef_min_bdrf : minimum number of Fourier term that must be used for 
  !                   surface BRDF/BPDF
  !  
  !   Svin     : real(4), normalized Stokes vector for incoming light
  !
  !   eps     : real, desired accuracy  
  !
  !   nmuext_in  : int, the number of extra (view) mu in XMUEXT_IN()   
  !
  !   xmuext_in  : real(NDmuext), contains extra (view) zenith angles cosines
  !             values between [0,1]
  !           
  !   nphout  : int, the number of (phi-phi_0) values of the emerging          
  !             radiation between [0,180] in XPHOUT() 
  !
  !   phout   : real(NDphiext), contains (phi-phi_0) values of the emerging          
  !             radiation between [0,180] 
  !
  !   nmu0    : int, the number of solar zenith angles mu_0 values 
  !
  !   xmu0    : real(nmu0), cosines of solar zenith angles
  !
  !   b_in       : real(0:NDlay), the optical thickness of the layer 
  !             b(1:nlayer)   bottom first
  !
  !   ssa     : real(NDlay) the single scattering albedo of layers
  !     
  !   alfbet_in  : real(4,4,0:NDcoef,NDlay) expansion coefficients of scattering
  !             matrix for each layer
  !             alfbet(4,4,0:ncoef,nlayer), bottom first
  !
  !   ncoef    : int(NDlay) number of Legendre coef for each layer
  !                   
  !   irf     : int, 0 or 1, whether the internal radiation field is computed or not
  !
  !   tlev    : temperature (K) of atmosphere levels 0 to NDlay 
  !             0 is the ground level
  !
  !   tsurf   : surface temperature (K)
  !
  !    wvnmlo, wvnmhi : Lower and upper wavenumber (inv cm) of spectral interval to feed DOAD_PLKAVG()
  ! 
  !  ===== OUTPUT :
  !
  !   flux    : real(NDmu0, 3, NDlay+1)  fluxes at each layer interface
  !                   related to incoming light
  !             flux(:,1,:) : downward flux (diffuse+direct)  
  !             flux(:,2,:) : upward flux  
  !             flux(:,3,:) : downward direct flux (transmitted)
  !
  !   SvR     : real(nmat,NDmu0,Ndmuext-NDmu0,NDphiext)
  !             Stoke vector of up radiance field ontop of the atmosphere
  !             related to incoming light 
  !   SvR_DW  : real(nmat,NDmu0,Ndmuext-NDmu0,NDphiext)
  !             Stoke vector of downward radiance field on surface
  !             related to incoming light 
  !
  !   flux_th : real(2, NDlay+1)  thermal in-situ fluxes at each layer interface
  !             flux(1,:) : downward flux (diffuse+direct)  
  !             flux(2,:) : upward flux  
  !
  !   SvR_th    : real(nmat, Ndmuext-NDmu0)
  !               Intensity of up radiance ontop of the atmosphere for thermal in-situ 
  !

  implicit none

  !     input
  logical, intent(in) :: verbose
  integer, intent(in) :: NDmu    
  integer, intent(in) :: NDmuext 
  integer, intent(in) :: NDsup   
  integer, intent(in) :: NDcoef  
  integer, intent(in) :: NDmu0   
  integer, intent(in) :: NDlay   
  integer, intent(in) :: NDirf   
  integer, intent(in) :: NDphiext  
  real(kind=8), intent(in) :: lambda
  integer, intent(in) :: nlayer
  integer, intent(in) :: nfoumx
  integer, intent(in) :: nmug_in
  integer, intent(in) :: nmat
  integer, intent(in) :: igrnd_in
  integer, intent(in) :: ncoef_min_bdrf
  real(kind=8), intent(in) :: Svin(4)
  real(kind=8), intent(in) :: eps
  integer, intent(in) :: nmuext_in
  real(kind=8), intent(in) :: xmuext_in(NDmuext)
  integer, intent(in) :: nphout
  real(kind=8), intent(in) :: phout(NDphiext)
  integer, intent(in) :: nmu0
  real(kind=8), intent(in) :: xmu0(NDmu0)
  real(kind=8), intent(in) :: b_in(0:NDlay)
  real(kind=8), intent(in) :: ssa(NDlay)
  real(kind=8), intent(in) :: alfbet_in(4,4,0:NDcoef,NDlay)
  integer, intent(in) :: ncoef(NDlay)
  real(kind=8), intent(in) :: surf_alb
  integer, intent(in) :: irf
  real(kind=8), intent(in) :: tlev(0:NDlay)
  real(kind=8), intent(in) :: tsurf
  real(kind=8), intent(in) :: wvnmlo
  real(kind=8), intent(in) :: wvnmhi

  !     output
  real(kind=8), intent(out) :: SvR(nmat,NDmu0,Ndmuext-NDmu0,NDphiext)
  real(kind=8), intent(out) :: SvR_DW(nmat,NDmu0,Ndmuext-NDmu0,NDphiext)
  real(kind=8), intent(out) :: flux(NDmu0, 3, NDlay+1)
  real(kind=8), intent(out) :: SvR_th(nmat,Ndmuext-NDmu0)
  real(kind=8), intent(out) :: flux_th(2, NDlay+1)

  ! ========================
  !     LOCAL VARIABLES      
  real(kind=8), parameter :: pi=3.1415926535897932384d0
  !  real(kind=8), parameter :: deg2rad = pi/180.0D0

  integer      :: imu0(NDmu0)
  real(kind=8) :: xmu(NDmu)
  real(kind=8) :: wmu(NDmu)
  ! -- surface reflectance and emissivity
  real(kind=8) :: refl(4,4,NDmu,NDmu,0:NDcoef)
  real(kind=8) :: emis(4,NDmu)
  !***********************************************************************
  ! internal field super matrices, dimension NDsup*NDsup                 *
  ! NDsup .ge. nmat*nmutot                                               *
  !***********************************************************************
  real(kind=8) :: Zmplus(NDsup,NDsup,NDlay)
  real(kind=8) :: Zmmin(NDsup,NDsup,NDlay)
  real(kind=8) :: en(NDmu,NDlay+1)
  real(kind=8) :: eb(NDmu,NDlay)
  real(kind=8) :: ek(NDmu,NDirf,NDlay)
  real(kind=8) :: taun(0:NDirf+1,0:NDlay)
  real(kind=8) :: etmu(NDmu,0:NDirf+1,0:NDlay)
  real(kind=8) :: ebtmu(NDmu,0:NDirf+1,NDlay)
  real(kind=8) :: ebmu(NDmu)
  real(kind=8) :: eblow(NDmu)
  real(kind=8) :: r(NDsup,NDsup)
  real(kind=8) :: t(NDsup,NDsup)
  real(kind=8) :: radd(NDsup,NDsup)
  real(kind=8) :: tadd(NDsup,NDsup)
  real(kind=8) :: ukk(NDsup,NDsup,NDirf,NDlay)
  real(kind=8) :: dkk(NDsup,NDsup,NDirf,NDlay)
  real(kind=8) :: unn(NDsup,NDsup,NDlay+1)
  real(kind=8) :: dnn(NDsup,NDsup,NDlay+1)

  !***********************************************************************
  ! the supervectors of the internal sources                             *
  !***********************************************************************
  real(kind=8) :: sr(NDsup)
  real(kind=8) :: st(NDsup)
  real(kind=8) :: sradd(NDsup)
  real(kind=8) :: stadd(NDsup)
  real(kind=8) :: sukk(NDsup,NDirf,NDlay)
  real(kind=8) :: sdkk(NDsup,NDirf,NDlay)
  real(kind=8) :: sunn(NDsup,NDlay+1)
  real(kind=8) :: sdnn(NDsup,NDlay+1)
  !***********************************************************************
  ! layer arrays, dimension NDlay                                        *
  ! number of homogeneous layers nlayer .le. NDlay                       *
  !***********************************************************************
  integer :: ndoubl(NDlay)
  integer :: M0(NDlay)
  integer :: M1(NDlay)
  integer :: M2(NDlay)
  !***********************************************************************
  !           former outputs
  !***********************************************************************
  real(kind=8) :: tau(0:NDirf+1,0:NDlay+1)
  real(kind=8) :: dsum(NDsup,NDmu0,NDphiext,0:NDirf+1,NDlay+1)
  real(kind=8) :: usum(NDsup,NDmu0,NDphiext,0:NDirf+1,NDlay+1)
  real(kind=8) :: fluxu(NDmu0,0:NDirf,NDlay+1)
  real(kind=8) :: fluxd(NDmu0,0:NDirf,NDlay+1)
  real(kind=8) :: etmu0(NDmu0,0:NDirf,NDlay+1)
  real(kind=8) :: sdsum(NDsup,0:NDirf+1,NDlay+1)
  real(kind=8) :: susum(NDsup,0:NDirf+1,NDlay+1)
  real(kind=8) :: sfluxu(0:NDirf,NDlay+1)
  real(kind=8) :: sfluxd(0:NDirf,NDlay+1)

  integer      :: nmuext
  integer      :: nmug
  real(kind=8) :: xmuext(NDmuext)
  integer      :: nlirf(0:NDlay+1)   !  real(0:NDlay+1), the number of internal field levels
  real(kind=8) :: b(0:NDlay)
  integer      :: igrnd
  real(kind=8) :: alfbet(4,4,0:NDcoef,NDlay)
  real(kind=8) :: ssor(0:NDlay)! real(0:NDlay), the strength of the internal sources 
  !                              ( 1 - a ) Planck(lambda,Temp) 
  !                              ssor(0) = ground level

  real(kind=8) :: alb, alb0, alb1, alb2, alb5, alb9
  real(kind=8) :: wi, wj

  integer :: Mmax

  logical :: flag_mu0
  integer :: i,j,k,m,il
  integer :: nsup
  integer :: nmutot
  integer :: imu
  integer :: iphi
  integer :: jmu0

  integer :: ncoef_brdf

  REAL(kind=8) :: DOAD_PLKAVG

  ! ======================

!!$  if (verbose) then
!!$     write(*,*) ''
!!$     write(*,*) ' (doad) '
!!$     write(*,*) ''
!!$  endif

  nmuext     = nmuext_in
  nmug       = nmug_in
  xmuext(:)  = xmuext_in(:)
  b(:)       = b_in(:)
  igrnd      = igrnd_in
  alfbet(:,:,:,:) = alfbet_in(:,:,:,:)

  if (igrnd.eq.3) then
     ! check array dimensioning
     if (ncoef_min_bdrf.gt.ndcoef) then
        write(*,*) ' '
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) ' (in sub doad)  :  ERROR '
        write(*,*) '  NDcoef must be >= ncoef_min_bdrf '
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) ' '        
        STOP
     end if
     ! set ncoef for BRDF/BPDF
     if (maxval(ncoef(:)).lt.ncoef_min_bdrf) then
        ncoef_brdf = ncoef_min_bdrf
     else
        ncoef_brdf = maxval(ncoef(:))
     end if
  end if

  if (maxval(ncoef).gt.ndcoef) then
     ! check array dimensioning
     write(*,*) ' '
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) ' (in sub doad)  :  ERROR '
     write(*,*) '  NDcoef must be >= MAXVAL(ncoef(:)) '
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) ' '        
     STOP
  end if

  ! init that to 0
  ssor(:)    = 0.0D0
  emis(:,:)  = 0.0D0
  nlirf(:)   = 0
  b(0)       = 0.d0

  ! the number of gaussian division points should be even,               
  ! then zero is not a gaussian division point                           
  if (((nmug/2)*2) .ne. nmug) nmug = nmug - 1

  ! set imu0 array form mu0 input     
  ! imu0    : int(NDmu0), index of solar zenith angles XMUEXT(IMU0())
  do i = 1, nmu0
     flag_mu0 = .false.
     do j = 1, nmuext
        if (xmu0(i) .eq. xmuext(j)) then
           imu0(i) = j               
           flag_mu0       = .true.
        endif
     enddo
     if (flag_mu0 .eqv. .false.) then
        !           we need to add an element to xmuext
        nmuext           = nmuext+1
        xmuext(nmuext)   = xmu0(i)
        imu0(i)          = nmuext
     endif
  enddo
  nmutot = nmug+nmuext
  nsup   = nmat*nmutot

  !-----------------------------------------------------------------------
  ! -- compute gauss xmu and weights for adding doubling computation
  call doad_gau_xmu(NDmu,NDmuext,nmug,nmuext,xmuext,xmu,wmu)
  ! -- xmu contains gauss xmu and external mu values
  ! -- wmu contains in fact sqrt(2.wi.xmu)

  ! =========  TEST ON ARGUMENTS =========
  !     taken from former sreadn.F       
  if (eps .lt. 1.d-10) then
     write(*,*) ''
     write(*,*) '(doad)'
     write(*,*) '(doad) eps should be > 1d-10'
     write(*,*) '(doad)'
     write(*,*) ''
     STOP
  endif

  if (nlayer .gt. NDlay) then
     print*,'(doad) The dimension of the layer matrices is too small'
     print*,'(doad)number of layers = ',nlayer,' but dim. is ', NDlay
     stop 
  endif

  if ((igrnd .lt. 0) .or. (igrnd .gt. 3)) then
     print*,'(doad) invalid value for igrnd (must be 1, 2 or 3)'
     stop 
  endif

  if ((irf .ne. 0) .and. (irf .ne. 1)) then
     print*,'(doad) invalid value for irf (must be 0 or 1)'
     stop 
  endif

  if (nmuext .gt. NDmuext) then
     print*,'(doad) The dimension of xmuext is too small'
     print*,'(doad) nmuext = ',nmuext,' while dimension is ',NDmuext
     stop 
  endif

  if (nmutot .gt. NDmu) then
     print*,'(doad) The dimension of the mu-matrices is too small'
     print*,'(doad) nmutot = ',nmutot,' while dim. is ',NDmu
     stop 
  endif

  if (nsup .gt. NDsup) then
     print*,'(doad) The dimension of the supermatrices is too small'
     print*,'(doad) nsup = ',nsup,' while dim. is ',NDsup
     stop
  endif

  if (nphout .gt. NDphiext) then
     print*,'(doad) The dimension of phout is too small'
     print*,'(doad) nphout = ',nphout,' while dimension is ',NDphiext
     stop 
  endif

  if (nmu0 .gt. NDmu0) then
     print*,'(doad) The dimension of imu0 is too small'
     print*,'(doad) nmu0 = ',nmu0,' while dimension is ',NDmu0
     stop 
  endif

  !     account for single scattering albedo in scattering matrix
  DO i = 1, nlayer         
     alfbet(1:nmat, 1:nmat, 0:ncoef(i), i) = alfbet(1:nmat, 1:nmat, 0:ncoef(i), i) * ssa(i)
  enddo

  !----------------------------------
  ! set ground reflectance properties
  refl(:,:,:,:,:) = 0.0D0    
  select case(igrnd)
  case(1)
     refl(1,1,1,1,0) = surf_alb

  case(2) 
     write(*,*) ' No Lommel-Seeliger surface implemented'

  case(3)
     !write(*,*) ' (doad) call doad_get_surf'
     call doad_get_surf(nmat, lambda, nmutot, xmu, ncoef_brdf, refl(1:nmat,1:nmat,1:nmutot,1:nmutot,0:NDcoef))

  end select

  !---------------------------
  !  Thermal emission settings
  ! if (irf .eq. 1 ) then

  ! print*, 'check in-situ thermal emission first :'
  ! print*, '     It seems to work only when there is no temperature gradiant within layers...'
  ! print*, '     and/or low opacities'
  ! print*, '        Bad interpolation ?'
  ! STOP 'doad_doad.f90' 

  !    !  Set thermal emission properties of layers
  !    do il = 1, nlayer
  !       ssor(il)  = (1.0D0 - ssa(il)) * &
  !            0.5D0 * &
  !            ( DOAD_PLKAVG( wvnmlo, wvnmhi, tlev(il-1)) + &
  !            DOAD_PLKAVG( wvnmlo, wvnmhi, tlev(il))       )
  !       if (ssor(il) .ne. 0.0D0) nlirf(il) = 1
  !    end do
  !    ! Surface
  !    if ((igrnd.eq.0).or.(tsurf.eq.0.0D0)) then
  !       ssor(0)  = 0.0D0
  !       nlirf(0) = 0
  !    ELSE

  !       print*, 'check surface emissivity first'
  !       STOP 'doad_doad.f90' 

  !       ssor(0)  = DOAD_PLKAVG( wvnmlo, wvnmhi, tsurf )        
  !       nlirf(0) = 1
  !       ! compute emissivities (angular integration)
  !       ! for each reflexion direction
  !       alb=0.
  !       alb0=0.
  !       do i=1,nmug
  !          alb1=0.
  !          alb5=0.
  !          alb9=0.
  !          alb2=0.
  !          ! integrate over incident directions
  !          do j=1,nmug
  !             wi=wmu(i)
  !             wj=wmu(j)
  !             alb=alb+refl(1,1,i,j,0)*wi*xmu(i)*xmu(j)*wj
  !             alb0=alb0+xmu(i)*xmu(j)*wi*wj
  !             alb1=alb1+refl(1,1,i,j,0)*wj
  !             alb5=alb5+refl(2,1,i,j,0)*wj
  !             alb9=alb5+refl(3,1,i,j,0)*wj
  !             alb2=alb2+refl(1,2,i,j,0)*wj ! to check
  !          enddo
  !          emis(1,i)=1.-alb1
  !          emis(2,i)=-alb5
  !          emis(3,i)=-alb9
  !          emis(4,i)=0
  !       enddo
  !       alb=alb*2./pi
  !       print*,"albedo surface =",alb,alb0
  !    ENDIF ! surface thermal emission
  ! endif ! thermal emission

  do il = 1, nlayer
     if (nlirf(il) .gt. NDirf) then
        print*,'(doad) The dimension of irf arrays is too small'
        print*,'(doad)  nlirf = ',nlirf(il),' while dimension is ',NDirf
        stop
     endif
     !     number of internal field levels must be uneven
     if ((2*int(nlirf(il)/2) .eq. nlirf(il)) .and. (nlirf(il).gt.0)) then
        write(*,*) '(doad) number of internal field levels must be uneven' 
        nlirf(il) = nlirf(il)-1
     endif
  enddo

  !***********************************************************************
  ! determine the number of fourier terms needed                         *
  !***********************************************************************
  call doad_fouMi(M0,M1,M2,nlayer,NDlay,alfbet,NDcoef,ncoef,eps)
  ! if BRDF/BPDF use, we may need to force M0(1)
  if ( (igrnd.eq.3)        .and. &
       (maxval(M0).lt.ncoef_brdf) ) then
     M0(1) =  ncoef_brdf
  endif
  call doad_maxM(M0,nlayer,NDlay,Mmax)

  !***********************************************************************
  !  fill some usefull arrays, e.g. optical depth, and exp(-tau/mu)       *
  !***********************************************************************
  call doad_fillar(b,NDlay,tau,taun,NDirf,nlayer,nlirf,xmu, & 
       NDmu,nmutot,                                         &
       en,eb,ek,etmu,ebtmu)
  if (Mmax .gt. nfoumx) Mmax = nfoumx

!   if (verbose) then
!      write(*,*) ' '
!      write(*,*) '  (doad) The number of necessary Fourier terms is ', Mmax
!      write(*,*) '  (doad) The max number of Fourier terms is ',nfoumx
!      write(*,*) ' '
!      do il=1,nlayer
!         write(*,85) il,M2(il),M1(il),M0(il)
! 85      format('   (doad) layer',i2,4x,'M2=',i4,4x,'M1=',i4,4x,'M0=',i5)
!      enddo
!      write(*,*) ' '
!   endif

  !***********************************************************************
  ! loop over the fourier terms m  
  !***********************************************************************
  !  Changed by MC
  !do m = 0,nfoumx
  do m = 0, Mmax

     call doad_initsm(radd,NDsup)
     call doad_initsm(tadd,NDsup)
     !    albedo changed to refl
     call doad_rgrndm(igrnd,m,refl,xmu,nmutot,nmat,wmu,radd, &
          NDcoef,NDmu,NDsup)
     if (m .eq. 0) then
        !   emis added
        call doad_ssurf(nmat,nmutot,emis,wmu,ssor,NDcoef,NDsup, &
             NDmu, NDlay, sradd,stadd)
     endif
     !***********************************************************************
     ! loop over the homogeneous layers il in loop 100                      *
     !***********************************************************************
     do  il=1,nlayer
        if (m .gt. M0(il)) then
           !***********************************************************************
           ! this fourier terms does not contribute to the radiation field in     *
           ! this layer except the direct contribution of incident radiation      *
           !***********************************************************************
           do k=1,nlirf(il)
              do j=1,nsup
                 do i=1,nsup
                    ukk(i,j,k,il) = 0.d0
                    dkk(i,j,k,il) = 0.d0
                 enddo
              enddo
           enddo
           do j=1,nsup
              do i=1,nsup
                 r(i,j) = 0.d0
                 t(i,j) = 0.d0
              
              enddo
           enddo
           

           goto 99
        endif
        !***********************************************************************
        ! calculate the m-th fouriercoefficient of the phasematrix             *
        !***********************************************************************
        call doad_expZm(m,il,NDlay,alfbet,NDcoef,ncoef,xmu,wmu,NDmu, &
             nmug,nmutot,nmat,Zmmin,Zmplus,NDsup,eps)
        call doad_bstart(m,il,NDlay,alfbet,NDcoef,ncoef,M0,    &
             b,eps,dabs(xmu(1)),                           &
             ndoubl,nlirf)
        !***********************************************************************
        ! reflection and transmission of an isolated homogeneous layer         *
        !***********************************************************************
        if (m .eq. 0) then
           call doad_shomly(nmat,nmutot,nmug,xmu,wmu,b,m,il,ebmu,r,t, & 
                nlirf,ndoubl,ukk,dkk,Zmmin,Zmplus, &
                NDsup,NDirf,NDmu,NDlay,ssor,sr,st,sukk,sdkk)
           !***********************************************************************
           ! internal radiation field of an isolated homogeneous layer            *
           !***********************************************************************
           call doad_shomin(nlirf,il,nmat,nmutot,nmug,ukk,dkk,ek, &
                sukk,sdkk,NDsup,NDirf,NDlay,NDmu)
        elseif (m .le. M1(il)) then
           call doad_homlay(nmat,nmutot,nmug,xmu,wmu,b,m,il,ebmu,r,t, &
                nlirf,ndoubl,ukk,dkk,Zmmin,Zmplus,NDsup,NDirf,NDmu,NDlay)
           call doad_homirf(nlirf,il,nmat,nmutot,nmug,ukk,dkk,ek, &
                NDsup,NDirf,NDlay,NDmu)
        elseif (m .le. M0(il)) then
           call doad_expbmu(b(il),xmu,NDmu,nmutot,ebmu)
           call doad_firstm(nmat,nmutot,xmu,wmu,b(il),ebmu,m,il,NDlay, &
                Zmmin,Zmplus,NDmu,NDsup,r,t)
           !write(*,*) 'After call to doad_firstm, r--> ', r
           call doad_ord1m(nmat,nmutot,xmu,wmu,taun,etmu,m,  &
                Zmmin,Zmplus, &
                ebtmu,nlirf,il,NDsup,NDmu,NDlay,NDirf,ukk,dkk)
        endif
        !***********************************************************************
        ! add the il-th layer, with r and t, to the one(s) below, with radd    *
        ! and tadd. On exit radd and tadd contain the reflection and transmis- *
        ! sion of the added layers, and r and t contain u and d of interface   *
        !***********************************************************************
99      call doad_expbmu(b(il),xmu,NDmu,nmutot,ebmu)
        if (il .eq. 1) then
           do i=1,nmutot
              eblow(i) = 1.d0
           enddo
        else
           call doad_expbmu(b(il-1),xmu,NDmu,nmutot,eblow)
        endif
        if (m .eq. 0) then
           call doad_sadd(nmat,nmutot,nmug,NDmu,NDsup,xmu,ebmu,eblow, &
                r,t,radd,tadd,sr,st,sradd,stadd)
           do j=1,nsup
              sunn(j,il) = sr(j)
              sdnn(j,il) = st(j)
           enddo
        else
           call doad_addtot(nmat,nmutot,nmug,NDmu,NDsup,xmu,ebmu,eblow, &
                r,t,radd,tadd)
        endif
        !***********************************************************************
        ! store u and d supermatrices for use in mullay                        *
        !***********************************************************************
        do j=1,nsup
           do i=1,nsup
              unn(i,j,il) = r(i,j)
              dnn(i,j,il) = t(i,j)              
           enddo
        enddo
     enddo ! loop over homogeneous layers
     if (m .eq. 0) then
        do j=1,nsup
           sunn(j,nlayer+1) = sradd(j)
           sdnn(j,nlayer+1) = 0.d0
        enddo
     endif
     do j=1,nsup
        do i=1,nsup
           unn(i,j,nlayer+1) = radd(i,j)
           dnn(i,j,nlayer+1) = 0.d0                      
        enddo

     enddo
     !***********************************************************************
     ! internal radiation field at interfaces of multi-layered atmosphere   *
     !***********************************************************************
     if (m .eq. 0) then
        call doad_smulay(nlayer,nmat,nmutot,nmug,unn,dnn,en,eb, &
             sunn,sdnn,NDsup,NDlay,NDmu)
     else
        call doad_mullay(nlayer,nmat,nmutot,nmug,unn,dnn,en,eb, &
             NDsup,NDlay,NDmu)
     endif
     !***********************************************************************
     ! internal radiation field of layered atmosphere                       *
     !***********************************************************************
     do il=1,nlayer
        if (m .eq. 0) then
           call doad_soremb(nlirf,il,nmat,nmutot,nmug,ukk,dkk, &
                unn,dnn,ek, &
                en,sukk,sdkk,sunn,sdnn,NDsup,NDirf,NDlay,NDmu)
        else
           call doad_irfemb(nlirf,il,nmat,nmutot,nmug, &
                ukk,dkk,unn,dnn,ek,en,NDsup,NDirf,NDlay,NDmu)
        endif
     enddo
     if (m .eq. 0) then
        call doad_fluxrt(2,nmat,nmutot,nmug,xmu,wmu,m,nlayer,nlirf, &
             ukk,dkk,unn,dnn,tau,en,NDmu,NDsup,NDirf,NDlay)
        call doad_fluxes(nmat,nmutot,nmug,xmu,wmu,m,nlayer,nlirf, &
             imu0,nmu0,fluxu,fluxd,etmu0, &
             ukk,dkk,unn,dnn,tau,en,NDmu,NDsup,NDirf,NDlay,NDmu0)
        call doad_sflux(nmat,nmutot,nmug,xmu,wmu,m,nlayer,nlirf, &
             sfluxu,sfluxd, &
             sukk,sdkk,sunn,sdnn,tau,NDmu,NDsup,NDirf,NDlay)
     endif
     !***********************************************************************
     ! add the m-th fourier term                                            *
     !***********************************************************************
     call doad_sumfou(nmat,nmug,nmuext,xmu,wmu,nphout, & 
          phout,nlirf,nlayer, &
          nmu0,imu0,m,ukk,unn,dkk,dnn,usum,dsum,Svin, &
          NDmu,NDmu0,NDsup,NDphiext,NDlay,NDirf)
     if (m .eq. 0) then
        call doad_sumsfm(nmat,nmug,nmuext,xmu,wmu,nlirf, &
             nlayer,m,sukk,sunn,sdkk,sdnn,susum,sdsum, &
             NDmu,NDsup,NDlay,NDirf)
     endif


  enddo ! loop on fourier term m

  !***********************************************************************
  ! write the result in the output variable
  do iphi= 1, nphout
     do imu = 1, nmuext_in
        do jmu0 = 1, nmu0
            
           
           !if ( xmuext_in(imu) .lt. 0) then
           !   write(*,*) 'UMU negatif, SVR', xmuext_in(imu) , &
           !   dsum(nmug+imu,jmu0,iphi,0,1)
            
           !end if
           !!! Line added on 30/04/15 for downward radiation
           SvR_DW(1,jmu0,imu,iphi) = dsum(nmug+imu,jmu0,iphi,0,1)
           !write(*,*) 'SvR_DW', SvR_DW(1,jmu0,imu,iphi)

           
           SvR(1,jmu0,imu,iphi)    = usum(nmug+imu,jmu0,iphi,0,nlayer+1)
                     
           !write(*,*) 'usum', usum(nmug+imu,jmu0,iphi,0,nlayer+1)
           if (nmat .gt.1) then           
              SvR(2,jmu0,imu,iphi) = usum(nmutot+nmug+imu,jmu0,iphi,0,nlayer+1)
              SvR(3,jmu0,imu,iphi) = usum(2*nmutot+nmug+imu,jmu0,iphi,0,nlayer+1)

              !! line added on 30/04/15
              SvR_DW(2,jmu0,imu,iphi) = dsum(nmutot+nmug+imu,jmu0,iphi,0,1)
              SvR_DW(3,jmu0,imu,iphi) = dsum(2*nmutot+nmug+imu,jmu0,iphi,0,1)
           endif
           if (nmat .gt.3) then
              SvR(4,jmu0,imu,iphi) = usum(3*nmutot+nmug+imu,jmu0,iphi,0,nlayer+1)
              SvR_DW(4,jmu0,imu,iphi) = dsum(3*nmutot+nmug+imu,jmu0,iphi,0,1)
           end if
        enddo
     enddo
  enddo
  if (irf.eq.1) then
     do imu = 1, nmuext_in
        SvR_th(1, imu)   = susum(nmug+imu,0,nlayer+1)
        if (nmat .gt.1) then           
           SvR_th(2,imu) = susum(nmutot+nmug+imu,0,nlayer+1)
           SvR_th(3,imu) = susum(2*nmutot+nmug+imu,0,nlayer+1)
        endif
        if (nmat .gt.3) SvR_th(4,imu) = susum(3*nmutot+nmug+imu,0,nlayer+1)
     enddo
  else
     SvR_th(:,:) = 0.0
  endif

  ! write the flux in output variable
  k = 0
  do jmu0 = 1, nmu0
     do il = 1, NDlay+1
        flux(jmu0,1,il) = (fluxd(jmu0,k,il) + etmu0(jmu0,k,il))
        flux(jmu0,2,il) = fluxu(jmu0,k,il)
        flux(jmu0,3,il) = etmu0(jmu0,k,il)
     end do
  end do
  if (irf.eq.1) then
     do il = 1, NDlay+1
        flux_th(1,il) = sfluxu(0,il)
        flux_th(2,il) = sfluxd(0,il) 
     end do
  else
     flux_th(:,:) = 0.0
  endif

  ! call doad_writ(usum,dsum,nmat,nphout,nlirf,nlayer, &
  !      nmug,nmuext,nmu0,xmu,imu0,phout,tau,b,fluxu,fluxd,etmu0, &
  !      NDmu,NDsup,NDmu0,NDphiext,NDirf,NDlay) 
  ! call doad_swrit(susum,sdsum,nmat,nlirf,nlayer, &
  !      nmug,nmuext,xmu,tau,b,sfluxu,sfluxd, &
  !      NDmu,NDsup,NDirf,NDlay)

  return

end subroutine doad
