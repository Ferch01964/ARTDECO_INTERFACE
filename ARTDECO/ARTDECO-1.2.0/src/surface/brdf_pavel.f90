
! This module was provided by Pavel Litvinov 
! at the Laboratoire d'Optique Atmospherique
! Two land BRDF are implemented (see Maignan, 2004, Remote Sensing of Environment, 2004, 90, 210-220):
!   The Ross-li with hot spot 
!   The Rahman – Pinty – Verstraete with hot spot
! The BPDF of Maignan, 2009, Remote Sensing of Environment, 2009, 113, 2642 - 2650
! is also accessible.


Module MBRDF_PAVEL

  Real(8), parameter, private :: dpi     = 3.141592653589793d0 
  Real(8), parameter, private :: d180_pi = 180.d0/dpi
  Real(8), parameter, private :: dpi_180 = dpi/180.d0

  PUBLIC  :: BRM
  PRIVATE :: BRDF_RPV, &
       LiSp_Ross, &
       Maign_Breon

Contains

  Subroutine BRM(warn, nmat, i_BRDF, mu0, muv, phi, m, n_par_brdf, par_brdf, Cv, Rip)

    Implicit none

    logical, intent(in) :: warn   ! whether the routine print warning or not
    Integer, intent(in) :: nmat   ! 1: BRDF only, 3:3x3 reflectivity matrix, 4: 4x4 reflectivity matrix
    Integer, intent(in) :: i_BRDF ! select the BRDF model to be used :
    !                               1 : The Rahman – Pinty – Verstraete with hot spot
    !                               2 : The Ross-li with hot spot 
    Real(8), intent(in) :: mu0    ! cosine of incident zenith angle
    Real(8), intent(in) :: muv    ! cosine of reflected zenith angle
    Real(8), intent(in) :: phi    ! relative azimuth angle
    Integer, intent(in) :: n_par_brdf ! number of parameters :
    !                                   must be 4 for   Rahman – Pinty – Verstraete
    !                                   and 3 for  Ross-li 
    Real(8), intent(in) :: par_brdf(n_par_brdf) ! Parameters for the BRDF :
    !                                    RPV model with HotSpot
    !                                      rho   = par_brdf(1)
    !                                      vk    = par_brdf(2)
    !                                      gteta = par_brdf(3)
    !                                      rhoH  = par_brdf(4)
    !                                    Ross_Li Sparse model with HotSpot
    !                                      par_brdf(1) = k(lambda)
    !                                      par_brdf(2) = k_2
    !                                      par_brdf(3) = k_1
    Real(8), intent(in)    :: Cv ! parameter of the Maignan 2009 BPDF model
    Complex(8), intent(in) :: m  ! complex refractive index for the land
    Real(8), intent (out)  :: Rip(1:nmat,1:nmat)

    Real(8) :: R11
    Real(8) :: R(4,4)

    select case(i_BRDF)
    case(0)
       if (n_par_brdf.ne.4) then
          write(*,*) ' ==============================================='
          write(*,*) '  (in sub. BRM)    ERROR              '
          write(*,*) '  Number of parameters for RPV model'
          write(*,*) '  must be 4 '
          write(*,*) ' ==============================================='       
          stop 'stop in MBRDF_PAVEL'       
       endif
       ! RPV model with HotSpot
       ! rho   = par_brdf(1)
       ! vk    = par_brdf(2)
       ! gteta = par_brdf(3)
       ! rhoH  = par_brdf(4)
       R11   = BRDF_RPV(mu0,muv,phi,par_brdf(1),par_brdf(2),par_brdf(3),par_brdf(4))

    case(1)

       if (n_par_brdf.ne.3) then
          write(*,*) ' ==============================================='
          write(*,*) '  (in sub. BRM)    ERROR              '
          write(*,*) '  Number of parameters for Li Sparse-Ross model'
          write(*,*) '  must be 3 '
          write(*,*) ' ==============================================='       
          stop 'stop in MBRDF_PAVEL'       
       endif
       ! Ross_Li Sparse model with HotSpot
       ! par_brdf(1) = k(lambda)
       ! par_brdf(2) = k_2
       ! par_brdf(3) = k_1
       call LiSp_Ross(mu0,0.d0,muv,phi,par_brdf(1),par_brdf(2),par_brdf(3), R11)

    case default 
       write(*,*) ' ====================================='
       write(*,*) '  (in sub. BRM)    ERROR              '
       write(*,*) ' i_BRDF=',i_BRDF,' - unknown value'
       write(*,*) ' ====================================='       
       stop 'stop in MBRDF_PAVEL'

    end select ! i_BRDF

    if (nmat.gt.1) then
       ! one parametric BPDF model of Maignan and Breon
       call Maign_Breon(warn,mu0,muv,phi,m,Cv,R)
       ! BRDF + Maignan_Fr_matr
       Rip(1:nmat,1:nmat) = R(1:nmat,1:nmat)
       Rip(1,1)           = Rip(1,1) + R11
    else
       Rip(1,1) = R11
    endif

  End Subroutine BRM

  !------------------------------------------------------

  !  Rahman Pinty Verstraete model
  Real(8) FUNCTION BRDF_RPV(cs,cv,phi,rho,vk,gt,rhoH)

    Implicit none

    Real(8), intent(in):: cs,cv,phi,rho,vk,gt,rhoH
    Real(8) xx,yy,zz,FF1,ww,aa,FF2,vv,G,FF3 

    xx=dabs(cs)**(vk-1.)
    yy=dabs(cv)**(vk-1.)
    zz=(dabs(cs)+dabs(cv))**(1.-vk)
    FF1=rho*xx*yy/zz
    xx=dsqrt(1.-cs*cs)
    yy=dsqrt(1.-cv*cv)
    ww=cs*cv+xx*yy*dcos(phi*dpi_180)
    aa=1+gt*gt+2*gt*ww
    FF2=(1.-gt*gt)/(aa**1.5)
    vv=xx/cs
    ww=yy/cv
    G=dsqrt(vv*vv+ww*ww-2*vv*ww*dcos(phi*dpi_180))
    !CD      FF3=1+(1-rho)/(1+G)
    FF3=1+(1.-rhoH)/(1.+G)
    BRDF_RPV=FF1*FF2*FF3
    return

  end FUNCTION BRDF_RPV

  !------------------------------------------------------

  ! Li Sparse- Ross model with azimuth function
  Subroutine LiSp_Ross(mu0, phi0, mu, Phir, klamb, k2, k1, LiD_R)

    Implicit none

    !  dl_par1=2.d0  !veg 
    Real(8), intent(in)  :: phi0
    Real(8), intent(in)  :: Phir
    Real(8), intent(in)  :: klamb
    Real(8), intent(in)  :: k2
    Real(8), intent(in)  :: k1
    Real(8), intent(out) :: LiD_R


    !  Parameter dl_par1 is diferent for soil and vegetation surfaces
    Real(8), parameter:: dl_par1=2., alpha0=1.5    !soil dl_par1=1.

    Real(8) cos_scat,sin_scat,tet_sc,Ross_K,LiD_K,mu0,mu 
    Real(8) tan_i,tan_v,D,sec_i_v,cos_t,sin_t,t,O
    Real(8) ni_x,ni_y,ni_z, nv_x,nv_y,nv_z, phi_inc, phi_v, sin_v,sin_inc,del_phi1
    Real(8) Hot_Sp,alpha

    tan_i=dsqrt(1.-mu0*mu0)/mu0
    tan_v=dsqrt(1.-mu*mu)/mu

    phi_v=Phir*dpi_180
    phi_inc=phi0*dpi_180

    sin_inc=dsqrt(1.-mu0*mu0)
    ni_x=sin_inc*dcos(phi_inc)
    ni_y=sin_inc*dsin(phi_inc)
    ni_z=mu0 !(cos of solar zenith angle)

    sin_v=dsqrt(1.-mu*mu)
    nv_x=sin_v*dcos(phi_v)
    nv_y=sin_v*dsin(phi_v)
    nv_z=mu !(cos of solar zenith angle)

    del_phi1=phi_v-phi_inc
    cos_scat=-(ni_x*nv_x+ni_y*nv_y+ni_z*nv_z)
    sin_scat=dsqrt(1-cos_scat*cos_scat)
    tet_sc=dacos(cos_scat)

    alpha=(dpi-tet_sc)*d180_pi
    Hot_Sp=(1.d0+1.d0/(1.d0+alpha/alpha0))

    Ross_K=Hot_Sp*((dpi*0.5d0-tet_sc)*cos_scat+sin_scat)/(mu+mu0)-dpi*0.25

    ! here positive angle corresponds to backscattering direction
    sec_i_v=(1.d0/mu+1.d0/mu0)
    !  D=dsqrt(tan_i*tan_i+tan_v*tan_v+2.d0*tan_i*tan_v*dcos(del_phi1))
    D=dsqrt(tan_i*tan_i+tan_v*tan_v-2.d0*tan_i*tan_v*dcos(del_phi1))

    cos_t=dsqrt(D*D+(tan_i*tan_v*dsin(del_phi1))**2)/sec_i_v*dl_par1
    If (dabs(cos_t) .gt. 1.) cos_t=1.

    ! cos_t=0.
    sin_t=dsqrt(1.-cos_t*cos_t)
    t=dacos(cos_t)
    O=1.d0/dpi*(t-sin_t*cos_t)*sec_i_v

    LiD_K=O-sec_i_v+(1.-cos_scat)/mu/mu0*0.5d0

    LiD_R=(LiD_K*k1+Ross_K*k2+1.)*klamb
    !LiD_R=(LiD_K*k1+Ross_K+k1)*klamb

    if (LiD_R .lt. 0.) LiD_R=0.

  End Subroutine LiSp_Ross

  !------------------------------------------------------

  ! one parametric BPDF model of Maignan and Breon
  Subroutine Maign_Breon(warn,mu01,muv1,phi,m,Cmgn,R_fr)

    Implicit none

    logical, intent(in) :: warn
    Complex(8), intent(in):: m
    Real(8), intent(in):: mu01,muv1,phi,Cmgn
    Real(8), intent(out):: R_fr(1:4,1:4)

    Real(8), parameter:: st_min=0.001,sin_eps=0.00000001
    Real(8) zz,alf,uu,koef,mu_alf !,pi
    Complex(8) r1,r2,r_par1,r_par2,r_perp1,r_perp2,r_par,r_perp
    Complex(8) r_par_per1,r_par_per2
    Real(8) r_par_sq,r_perp_sq
    Real(8) del_phi,sin_ti,sin_tv,sin_tet,cos_tet,cos_etv,cos_eti
    Real(8) sin_eti,sin_etv,cos_2etv,sin_2etv,cos_2eti,sin_2eti
    Real(8) L_v(1:4,1:4), L_inc(1:4,1:4), F_fr(1:4,1:4)
    Real(8) mu0,muv,tet0v,sin_phi
  
    mu0=mu01
    muv=muv1

    !  sin_eps=dsin(st_min*dpi_180)

    sin_tv=dsqrt(1.d0-muv*muv)
    sin_ti=dsqrt(1.d0-mu0*mu0)
    zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
    cos_tet=-zz
    sin_tet=dsqrt(1.d0-cos_tet*cos_tet)


    If (dabs(sin_tv) .lt. sin_eps) then
       tet0v=dacos(muv)  
       tet0v=dabs(tet0v-st_min*dpi_180)
       muv=dcos(tet0v)
       sin_tv=dsqrt(1.d0-muv*muv)
       zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
       cos_tet=-zz
       sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
    endif

    If (dabs(sin_ti) .lt. sin_eps) then
       tet0v=dacos(mu0)  
       tet0v=dabs(tet0v-st_min*dpi_180)
       mu0=dcos(tet0v)
       sin_ti=dsqrt(1.d0-mu0*mu0)
       zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
       cos_tet=-zz
       sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
    endif

    if (dabs(sin_tet) .lt. sin_eps) then
       tet0v=dacos(mu0)  
       tet0v=dabs(tet0v-st_min*dpi_180)
       muv=dcos(tet0v)
       sin_tv=dsqrt(1.d0-muv*muv)
       zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
       cos_tet=-zz
       sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
    endif

    if (warn) then
       if (dabs(sin_tet) .lt. sin_eps) write(*,*) 'sin_tet',sin_tet
       if (dabs(sin_tv)  .lt. sin_eps) write(*,*) 'sin_tv',sin_tv
       if (dabs(sin_ti)  .lt. sin_eps) write(*,*) 'sin_ti',sin_ti
    endif

    alf=dacos(zz)
    uu=dtan(alf*0.5d0)
    mu_alf=dcos(alf*0.5d0)
    koef=Cmgn*dexp(-uu)*0.25d0/(mu0+muv)

    r1=m*m*mu_alf
    r2=cdsqrt(m*m-1.d0+mu_alf*mu_alf)
    r_par1=r1-r2
    r_par2=r1+r2

    r_perp1=mu_alf-r2
    r_perp2=mu_alf+r2

    r_par=r_par1/r_par2
    r_perp=r_perp1/r_perp2

    r_par_sq=dble(r_par*dconjg(r_par))
    r_perp_sq=dble(r_perp*dconjg(r_perp))
    r_par_per1=r_par*dconjg(r_perp)
    r_par_per2=dconjg(r_par_per1)

    F_fr=0.d0
    F_fr(1,1)=(r_par_sq+r_perp_sq)*koef*0.5d0
    F_fr(2,2)= F_fr(1,1)
    F_fr(2,1)=(r_par_sq-r_perp_sq)*koef*0.5d0
    F_fr(1,2)= F_fr(2,1)
    F_fr(3,3)=dble(r_par_per1+r_par_per2)*koef*0.5d0
    F_fr(4,4)=F_fr(3,3)
    F_fr(3,4)=dble((0,-1.d0)*(r_par_per1-r_par_per2))*koef*0.5d0
    F_fr(4,3)=-F_fr(3,4)

    del_phi=phi*dpi_180-dpi
    cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
    cos_eti= (muv+mu0*cos_tet)/sin_tet/sin_ti !!! checked "-"

    sin_etv=dsin(del_phi)*sin_ti/sin_tet
    sin_eti=dsin(del_phi)*sin_tv/sin_tet

    cos_2etv=cos_etv*cos_etv-sin_etv*sin_etv
    sin_2etv=2.d0*cos_etv*sin_etv

    cos_2eti=cos_eti*cos_eti-sin_eti*sin_eti
    sin_2eti=2.d0*cos_eti*sin_eti

    L_v=0.d0
    L_inc=0.d0

    sin_phi=dsin(phi*dpi_180)
    If (dabs(sin_phi) .lt. sin_eps/10.d0) then
       L_v(1,1)=1.d0
       L_v(2,2)=1.d0
       L_v(2,3)=0
       L_v(3,3)=1.d0
       L_v(3,2)=0
       L_v(4,4)=1.d0

       L_inc(1,1)=1.d0
       L_inc(2,2)=1.d0
       L_inc(2,3)=0
       L_inc(3,3)=1.d0
       L_inc(3,2)=0
       L_inc(4,4)=1.d0
    else
       L_v(1,1)=1.d0
       L_v(2,2)=cos_2etv
       L_v(2,3)=sin_2etv
       L_v(3,3)=cos_2etv
       L_v(3,2)=-sin_2etv
       L_v(4,4)=1.d0

       L_inc(1,1)=1.d0
       L_inc(2,2)=cos_2eti
       L_inc(2,3)=sin_2eti
       L_inc(3,3)=cos_2eti
       L_inc(3,2)=-sin_2eti
       L_inc(4,4)=1.d0
    endif

    R_fr=Matmul(L_v,Matmul(F_fr,L_inc))

  end Subroutine Maign_Breon

End Module MBRDF_PAVEL
