      subroutine adding(NDmu, NDsup, NDlay, NDcoef, NDgeom
     +                 , a, b, coefs, ncoefs, nlayer
     +                 , theta, theta0, phi, ngeom
     +                 , nmug, nfoumax, mstop, nmat, epsilon
     +                 , iface, xm
     +                 , R, T, R1, T1, Rst, Tst, R1st, T1st
     +                 , Ri, R1i, Di, D1i
     +                 , Rflux, Tflux, URU, UTU
     +                 , Riflux, Tiflux, Diflux, URUi, UTUi, UDUi )
*----------------------------------------------------------------------*
*                                                                      *
*                     A D D I N G    M E T H O D                       *
*               F O R    P O L A R I Z E D   L I G H T                 *
*                                                                      *
*                   V.L. Dolman and J. Chowdhary                       *
*                          Free University                             *
*               Department of Physics and Astronomy                    *
*                         De Boelelaan 1081                            *
*                 1081 HV Amsterdam, The Netherlands                   *
*                                                                      *
*  Documentation can be found in :                                     *
*     [1] J.F. de Haan, P.B. Bosma and J.W. Hovenier (1987):           *
*         Astronomy and Astrophysics 183, pages 371-391.               *
*     [2] J.F. de Haan (1987):                                         *
*         'Effects of aerosols on the brightness and polarization of   *
*         cloudless planetary atmospheres', Thesis Free University     *
*         Amsterdam.                                                   *
*     [3] J.W. Hovenier (1970):                                        *
*         'Polarized light in planetary atmospheres',                  *
*         Thesis University of Leiden.                                 *
*     [4] J.W. Hovenier (1969):                                        *
*         Journal of Atmospheric Sciences 26, pages 488-499.           *
*     [5] W.A. de Rooij (1985):                                        *
*         'Reflection and transmission of polarized light by planetary *
*         atmospheres', Thesis Free University Amsterdam.              *
*     [6] J. Chowdhary (1990):                                         *
*         Incorporation of a smooth water-air interface in multiple    *
*         scattering calculations using the adding method for          *
*         polarized light, Graduation report, Free University Amsterdam*
*         Astronomy Department, pp. 60.                                *
*  In the comments, equation numbers generally refer to reference [1]  *
*  which is called the 'adding paper' unless stated otherwise.         *
*                                                                      *
*  Last modification:                                                  *
*      Reorganize datastructures, particularly for fluxes.             *
*      Cosmetic improvements.                                          *
*                                                                      *
*                              V.L. Dolman February 1993               *
*----------------------------------------------------------------------*
*  In the 'adding' subroutine, subroutines are called according to     *
*  the following scheme:                                               *
*                                                                      *
*          adding  ( setmu  ( member                                   *
*                  (        ( ad_gauleg                                *
*                  ( setfou                                            *
*                  ( prcoef                                            *
*                  ( initRT ( first  ( setgsf                          *
*                  (                 ( setexp                          *
*                  (                 ( fasemx                          *
*                  ( initi  ( firsti                                   *
*                  (        ( Fmatri                                   *
*                  (                                                   *
*                  ( setZm  ( transf                                   *
*                  (        ( scalZm                                   *
*                  (                                                   *
*                  ( renorm ( brack                                    *
*                  ( expbmu                                            *
*                  ( firstm                                            *
*                  ( layerm ( expbmu                                   *
*                  (        ( ord1m                                    *
*                  (        ( ord2m  ( fillup                          *
*                  (        ( bstart                                   *
*                  (        ( double                                   *
*                  ( add1m                                             *
*                  ( addlay ( notop                                    *
*                  (        ( nobot                                    *
*                  (        ( oneadd                                   *
*                  ( top2bot                                           *
*                  ( newfou                                            *
*                  ( endfou                                            *
*                  ( ad_fluxes ( intang                                *
*                  (                                                   *
*                  ( addi   ( rRprod                                   *
*                  (        ( lRprod                                   *
*                  ( addim1                                            *
*                  ( fluxi  ( intang                                   *
*                  (        ( lTprod                                   *
*                  (        ( rTprod                                   *
*                  (        ( fxdiri ( dirmat                          *
*                  (                 ( sstari                          *
*                  ( newfoui                                           *
*                                                                      *
*  Not indicated are calls to subroutines from the                     *
*  subroutines from the supermatrix package (in supermatrix.f):        *
*         addSM                                                        *
*         assign                                                       *
*         prod                                                         *
*         star                                                         *
*         trace                                                        *
*         rdiapr                                                       *
*         ldiapr                                                       *
*         Tstar                                                        *
*         zero                                                         *
*----------------------------------------------------------------------*
*  On entry:                                                           *
*          
*     NDmu    : maximum number of distinct (!) mu values, counting     *
*               both integration points and extra points               *
*               NDmu must be > (NDgeom + nstr)                         *
*     NDsup   : size of a supermatrix (number of rows or columns)      *
*                NDsup  must be > (NDmu * nmat=3)                      *
*     NDlay   : maximum number of layers in the atmosphere             *
*     NDcoef  : maximum order of terms in the expansion of the         *
*               scattering matrix in generalized spherical functions.  *
*     NDgeom  : maximum number of different combinations of viewing    *
*               directions and directions of incidence for which       *
*               results are calculated                                 *
*  Layer data (layers are numbered from bottom to top!):               *
*     a       : single scattering albedo                               *
*     b       : optical thickness                                      *
*     coefs   : expansion coefficients of the scattering matrix in     *
*               generalized spherical functions                        *
*               normalized so that coefs(1,1,0,layer) = 1              *
*     ncoefs  : order of the highest nonzero term in the expansion     *
*     nlayer  : number of layers in the atmosphere                     *
*                                                                      *
*  Geometry data :                                                     *
*     theta   : the angle between the viewing direction and the        *
*               normal, in degrees                                     *
*     theta0  : solar zenith angle in degrees                          *
*     phi     : azimuth of the emerging direction relative to the      *
*               direction of incidence in degrees                      *
*               measured clockwise when looking upward                 *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*                                                                      *
*  Technical data:                                                     *
*     nmug    : number of Gauss points for the mu integrations         *
*     nfoumax : maximum Fourier index                                  *
*     mstop   : indicator of Fourier series convergence should be      *
*               reached for all directions (0) or only the extra       *
*               geometries (1)                                         *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     epsilon : requested accuracy                                     *
*                                                                      *
*  Interface data :                                                    *
*     iface   : indictor that interface is present (iface = 1) or      *
*               absent (iface = 0)                                     *
*     xm      : index of refraction of water with respect to air       *
*----------------------------------------------------------------------*

      implicit double precision (a-h,o-z)
      dimension a(NDlay), b(NDlay)
     +        , coefs(4,4,0:NDcoef,NDlay), ncoefs(NDlay)
      dimension theta(NDgeom), theta0(NDgeom), phi(NDgeom)
*----------------------------------------------------------------------*
*  On exit:                                                            *
*                                                                      *
*  Radiative transfer matrices for each geometry for the atmosphere    *
*  alone (i.e. no interface, no surface):                              *
*     R    : reflection matrix                                         *
*     T    : transmission matrix                                       *
*     R1   : first order scattering contribution to reflection matrix  *
*     T1   : first order scattering contribution to transmission matrix*
*     Rst  : same as R but for illumination from below                 *
*     Tst  : same as T but for illumination from below                 *
*     R1st : same as R1 but for illumination from below                *
*     T1st : same as T1 but for illumination from below                *
*                                                                      *
*  Radiative transfer matrices for each geometry for the atmosphere    *
*  above a water-air interface (nothing below the interface):          *
*     Ri   : diffuse reflection matrix (no sunglint)                   *
*     Di   : downward radiation matrix just above the interface        *
*     R1i  : first order scattering contribution to Ri                 *
*     D1i  : first order scattering contribution to Di                 *
*                                                                      *
*  Flux data :                                                         *
*  Fluxes generally are of the form :                                  *
*                                                                      *
*                 1                                                    *
*     AU(mu)  = 2 int dmu0 mu0 Am(mu,mu0)                              *
*                 0                                                    *
*                                                                      *
*                 1                                                    *
*     UA(mu0) = 2 int dmu mu Am(mu,mu0)                                *
*                 0                                                    *
*                                                                      *
*                 1            1                                       *
*     UAU    = 2 int dmu mu 2 int dmu0 mu0 Am(mu,mu0)                  *
*                 0            0                                       *
*                                                                      *
*  where A is the reflection or transmission matrix and Am its         *
*  m=0 Fourier component in azimuth. We use mu=cos(theta) and          *
*  mu0=cos(theta0).                                                    *
*  The fluxes that have an argument mu or mu0 are calculated for all   *
*  geometries. We use the following datastructure to store them:       *
*                                                                      *
*      Aflux( k,l, igeom, muarg, illum, iorder ).                      *
*                                                                      *
*  The meaning of the indices is as follows:                           *
*                                                                      *
*     index    values           meaning                                *
*   -----------------------------------------------------------------  *
*     k,l      1,...,nmat       indices in a (nmat X nmat) matrix      *
*                                                                      *
*     igeom    1,...,ngeom      index of the geometry considered       *
*                                                                      *
*     muarg    Nmu0             integral over outgoing direction, so   *
*                               we have UA(mu0)                        *
*              Nmu              integral over incoming direction, so   *
*                               we have AU(mu)                         *
*                                                                      *
*    illum     Nabove           illumination from above                *
*              Nbelow           illumination from below                *
*                                                                      *
*    iorder    Nzero            only zero order scattering             *
*              Nfirst           only first order scattering            *
*              Nall             all orders of scattering (including 0) *
*   -----------------------------------------------------------------  *
*                                                                      *
*  The values Nmu0, Nmu, Nabove, Nbelow, Nzero, Nfirst and Nall  are   *
*  coded in a parameter statement.                                     *
*  For the spherical flux UAU we use the following datastructure:      *
*                                                                      *
*      UAU( k,l, illum, iorder ).                                      *
*                                                                      *
*  The meaning of the indices is the same as before.                   *
*  In the case that a water-air interface is present, there is a       *
*  subtlety regarding the arguments mu and mu0. Arguments of fluxes    *
*  that pertain to directions below the interface are always in the    *
*  form xmut(mu) which is defined in Eq. (7) of Chowdhary (1991).      *
*  For example we have                                                 *
*                                                                      *
*    Tiflux( .., .., .., muarg, Nabove, .. ) = TU(xmut(mu))            *
*                                                                      *
*  Physically this corresponds to isotropic illumination from above,   *
*  and a detector below the interface looking in a direction with      *
*  a cosine of zenith angle given by xmut(mu). (No water medium is     *
*  assumed present.)                                                   *
*                                                                      *
*  We have the following flux arrays:                                  *
*     Rflux and URU   : reflection by the atmosphere alone             *
*     Tflux and UTU   : transmission by the atmosphere alone           *
*     Riflux and URUi : reflection by the atmosphere-interface system  *
*     Tiflux and UTUi : transmission by the atmosphere-interface system*
*     Diflux and UDUi : downward radiation just above the interface    *
*----------------------------------------------------------------------*
      dimension T(4,4,NDgeom),   Tst(4,4,NDgeom)
     +        , R(4,4,NDgeom),   Rst(4,4,NDgeom)
     +        , T1(4,4,NDgeom),  T1st(4,4,NDgeom)
     +        , R1(4,4,NDgeom),  R1st(4,4,NDgeom)
      dimension Ri(4,4,NDgeom),  Di(4,4,NDgeom)
     +        , R1i(4,4,NDgeom), D1i(4,4,NDgeom)
      parameter( Nmu=2
     +         , Nall=3
     +         , Nbelow=2 )
      dimension Rflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , Tflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , URU(  4,4,              Nbelow, Nall)
     +        , UTU(  4,4,              Nbelow, Nall)
      dimension Riflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , Diflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , Tiflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , URUi(  4,4,              Nbelow, Nall)
     +        , UDUi(  4,4,              Nbelow, Nall)
     +        , UTUi(  4,4,              Nbelow, Nall)
*----------------------------------------------------------------------*
*  Constants important internally                                      *
*     NDsup   : size of a supermatrix (number of rows or columns)      *
*     NDmu    : maximum number of distinct (!) mu values, counting     *
*               both integration points and extra points               *
*     NDfou   : maximum number of Fourier terms, equal to NDcoef       *
*               (MC - obsolete : NDfou --> NDcoef)                     *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*  Layer data used internally:
*     e       : direct transmission through layer in viewing direction *
*     e0      : same as e but for solar direction                      *
*----------------------------------------------------------------------*
      dimension e(NDlay,NDgeom), e0(NDlay,NDgeom)
*----------------------------------------------------------------------*
*  Geometry data used internally:                                      *
*     xmu     : all different mu values (integration and extra points) *
*     smf     : supermatrix factors dsqrt(2*w*mu), or 1 for extra pts. *
*     nmug    : number of Gauss points for mu integrations             *
*     nmutot  : total number of distinct mu points                     *
*     imu     : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu value can be found  *
*     imu0    : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu0 value can be found *
*----------------------------------------------------------------------*
      dimension xmu(NDmu), smf(NDmu), imu(NDgeom), imu0(NDgeom)
*----------------------------------------------------------------------*
*  Technical data used internally:                                     *
*     M0      : highest Fourier term that must be treated with         *
*               doubling                                               *
*     M1      : highest Fourier term that must be treated with the     *
*               sum of first and second order scattering               *
*     M2      : highest Fourier term that must be treated with         *
*               single scattering (higher terms have no scattering)    *
*     iadd    : option indicating kind of adding procedure :           *
*               1 = normal adding                                      *
*               2 = top layer has no scattering                        *
*               3 = bottom layer has no scattering                     *
*----------------------------------------------------------------------*
      dimension M0(NDlay), M1(NDlay), M2(NDlay), iadd(NDlay,0:NDcoef)
*----------------------------------------------------------------------*
*  Bottom layer data used internally:                                  *
*     Rmbot   : Fourier component of reflection supermatrix            *
*     Rmsbot  : same as Rmbot but for illumination from below          *
*     Rm1bot  : first order contribution to Rmbot                      *
*     Rm1sbot : first order contribution to Rmsbot                     *
*                                                                      *
*     Tmbot   : Fourier component of transmission supermatrix          *
*     Tm1bot  : first order contribution to Tmbot                      *
*               (the transmission for illumination from below is       *
*               obtained from symmetry relations)                      *
*     ebbot   : dexp(-b/mu) where b is the optical thickness of        *
*               the bottom layer                                       *
*----------------------------------------------------------------------*
      dimension Rmbot(NDsup,NDsup),  Rmsbot(NDsup,NDsup) 
     +        , Rm1bot(NDsup,NDsup), Rm1sbot(NDsup,NDsup)
     +        , Tmbot(NDsup,NDsup),  Tm1bot(NDsup,NDsup), ebbot(NDmu)
*----------------------------------------------------------------------*
*  Top layer data used internally:                                     *
*     Rmtop  : Fourier component of reflection supermatrrix            *
*     Rm1top : first order contribution to Rmtop                       *
*                                                                      *
*     Tmtop  : Fourier component of transmission supermatrix           *
*     Tm1top : first order contribution to Tmtop                       *
*              (data for illumination from below is obtained from      *
*              symmetry relations for a homogeneous layer)             *
*     ebtop   : dexp(-b/mu) where b is the optical thickness of        *
*               the top layer                                          *
*----------------------------------------------------------------------*
      dimension Rmtop(NDsup,NDsup), Rm1top(NDsup,NDsup)
     +        , Tmtop(NDsup,NDsup), Tm1top(NDsup,NDsup), ebtop(NDmu)
*----------------------------------------------------------------------*
*     xRm1top ... xTm1sbot : scratch space for treatment of first      *
*            order m-th Fourier component                              *
*----------------------------------------------------------------------*
      dimension xRm1top(4,4,NDgeom), xRm1stop(4,4,NDgeom)
     +        , xTm1top(4,4,NDgeom), xTm1stop(4,4,NDgeom)
     +        , xRm1bot(4,4,NDgeom), xRm1sbot(4,4,NDgeom)
     +        , xTm1bot(4,4,NDgeom), xTm1sbot(4,4,NDgeom)
      dimension xRm1i(4,4,NDgeom), xDm1i(4,4,NDgeom)
*----------------------------------------------------------------------*
*  Interface data used internally:                                     *
*     xm    : index of refraction                                      *
*     xmut  : mut(mu) defined by Chowdhary (1991) Eq. (7) for each     *
*             value of mu                                              *
*     Rf    : Fresnel matrix describing refraction by interface        *
*     Tf    : Fresnel matrix describing transmission through interface *
*     Tfs   : same as Tf but for illumination from beneath             *
*     Rfs   : same as Rf but for illumination from beneath             *
*----------------------------------------------------------------------*
      dimension xmut(NDmu)
      dimension Rf(NDmu,3), Tf(NDmu,3), Tfs(NDmu,3), Rfs(NDmu,3)
*----------------------------------------------------------------------*
*  Combined interface-layers data used internally:                     *
*     Wms  : dummy-supermatrix, calculated when adding the bottomlayer *
*            to the interface and used to calculate fluxes             *
*     Wm1s : same as Wms, used to calculate 1st order fluxes.          *
*----------------------------------------------------------------------*
      dimension Wms(NDsup,NDsup), Wm1s(NDsup,NDsup)
*----------------------------------------------------------------------*
*  Phase matrix data used internally:                                  *
*     Zmplus : Fourier component of the phase matrix Zm(mu,mu')        *
*     Zmmin  : Fourier component of the phase matrix Zm(-mu,mu')       *
*----------------------------------------------------------------------*
      dimension Zmplus(NDsup,NDsup), Zmmin(NDsup,NDsup)
*----------------------------------------------------------------------*
*  Auxiliary data used internally:                                     *
*----------------------------------------------------------------------*
      logical nextm, verbo
*----------------------------------------------------------------------*
*     Test of constant values (array dimensions)
      call test_constant(NDmu, NDsup, NDgeom, nmat, nmug)
*----------------------------------------------------------------------*
      verbo = .false.
      call setmu( NDmu, NDgeom, theta, theta0, ngeom, nmug
     +          , imu, imu0, xmu, smf, nmutot )
      call setfou(NDmu, NDlay, NDcoef, coefs, ncoefs, nlayer, a, b
     +           , xmu, nmutot, epsilon, M0, M1, M2, iadd )
      if (verbo) call prcoef(NDlay, NDcoef, nlayer, nmat, coefs, ncoefs)
      call initRT(NDmu, NDlay, NDgeom, NDcoef
     +           , coefs, ncoefs, xmu, imu0, imu, phi
     +           , ngeom, nmat
     +           , a, b, nlayer
     +           , R1, T1, R1st, T1st, e, e0 
     +           , R , T , Rst , Tst )
*----------------------------------------------------------------------*
*  If needed, make preparations for treatment of interface             *
*----------------------------------------------------------------------*
      if (iface .eq. 1)
     +    call initi( NDmu, NDlay, NDgeom 
     +              , xmu, imu, imu0, nmutot, ngeom, nmat, nlayer, xm
     +              , R1, T1, R1st, T1st, e, e0
     +              , Rf, Tf, Tfs, Rfs, xmut
     +              , R1i, D1i, Ri, Di )
*----------------------------------------------------------------------*
*  Fourier loop :                                                      *
*----------------------------------------------------------------------*
      m = -1
 1000     m = m+1
          if (verbo) print *,' adding: treat Fourier term m=',m
          do 900 layer=1, nlayer
              if (verbo) print *,' adding: treat layer ', layer
              call setZm(NDmu, NDsup, NDlay, NDcoef
     +                  , m, layer, coefs, ncoefs, xmu, smf
     +                  , epsilon
     +                  , nmug, nmutot,nmat, Zmmin, Zmplus )
                if (m .eq. 0) call renorm(NDmu, NDsup 
     +                                 , Zmmin, Zmplus, nmug, nmutot
     +                                 , nmat, xmu, smf, epsilon )
              call expbmu(NDmu, b(layer), xmu, nmutot, ebtop )
              call firstm(NDmu, NDsup, NDgeom
     +                  , m, xmu, imu, imu0, smf, nmutot, ngeom, nmat
     +                  , Zmplus, Zmmin, a(layer), b(layer), ebtop
     +                  , xRm1top, xRm1stop, xTm1top, xTm1stop )
              call layerm(NDmu, NDsup, NDlay, NDcoef
     +               , m, M0, M1, M2, layer, xmu, smf
     +               , nmug, nmutot, nmat, epsilon
     +               , coefs, ncoefs, Zmplus, Zmmin, a(layer), b(layer)
     +               , ebtop, Rm1top, Rmtop, Tm1top, Tmtop )
              if (verbo) print *,' adding: layer properties calculated'
              if (layer .gt. 1) then
                  call addm1(NDmu, NDgeom, ngeom, imu, imu0, nmat
     +                      , xRm1top, xRm1stop, xTm1top, xTm1stop
     +                      , ebtop, ebbot
     +                      , xRm1bot, xRm1sbot, xTm1bot, xTm1sbot )
                  call addlay(NDmu, NDsup
     +                       , Rmtop, Tmtop, Rmbot, Tmbot, Rmsbot
     +                       , Rm1top, Tm1top, Rm1bot, Rm1sbot, Tm1bot
     +                       , ebtop, ebbot, iadd(layer,m)
     +                       , nmutot, nmug, nmat )
                  if (verbo) print *,' adding: layer is added to total'
              else
                  call top2bot(NDmu, NDsup, NDgeom, nmat, nmutot, ngeom
     +                        , Rmtop, Tmtop, Rm1top, Tm1top, ebtop
     +                        , Rmbot, Rmsbot, Tmbot
     +                        , Rm1bot, Rm1sbot, Tm1bot, ebbot
     +                        , xRm1top, xRm1stop, xTm1top, xTm1stop
     +                        , xRm1bot, xRm1sbot, xTm1bot, xTm1sbot )
                  if (verbo) print *,' adding: toplayer data'
     +                              ,' copied to bottomlayer'
              endif
  900     continue
*----------------------------------------------------------------------*
*         Add new Fourier term to the summation for reflection and     *
*         transmission for illumination from above and below           *
*----------------------------------------------------------------------*
          call newfou(NDmu, NDsup, NDgeom 
     +               , m, Rmbot, Tmbot, xRm1bot, xTm1bot, xTm1sbot
     +               , Rmsbot, xRm1sbot
     +               , xmu, imu, imu0, smf, phi
     +               , nmat, ngeom
     +               , R, T, Rst, Tst 
     +               , R1, T1, R1st, T1st )
*----------------------------------------------------------------------*
*         Check convergence of Fourier series                          *
*         This concerns the atmosphere only                            *
*----------------------------------------------------------------------*
          call endfou(NDmu, NDsup, NDlay, NDgeom
     +               , m, M1, M0, nfoumax, nlayer, epsilon, mstop
     +               , xmu, imu, imu0, nmutot, nmat, ngeom
     +               , Rm1bot, Rmbot, Tm1bot, Tmbot, nextm )
*----------------------------------------------------------------------*
*         Calculate fluxes                                             *
*----------------------------------------------------------------------*
          if (m .eq. 0) 
     +              call ad_fluxes(NDmu, NDsup,NDgeom
     +                  , xmu, smf, imu, imu0
     +                  , ngeom, nmug, nmutot, nmat, ebbot
     +                  , Rmbot, Rmsbot, Tmbot, Rm1bot, Rm1sbot, Tm1bot
     +                  , Rflux, Tflux, URU, UTU )
*----------------------------------------------------------------------*
*         If desired, add atmosphere and interface and calculate       *
*         fluxes pertaining to the atmosphere-interface system         *
*----------------------------------------------------------------------*
          if (iface .eq. 1) then
              call addi(NDmu, NDsup
     +                 , m, nmutot, nmug, nmat, smf, Rf, xmu
     +                 , Rmbot,  Tmbot,  Rmsbot, Wms, ebbot
     +                 , Rm1bot, Tm1bot, Rm1sbot, Wm1s )
              call addim1(NDmu, NDgeom, xmu, imu, imu0, ngeom, nmat
     +                 , Rf, ebbot
     +                 , xRm1bot, xRm1sbot, xTm1bot, xTm1sbot
     +                 , xRm1i, xDm1i )
              if (m .eq. 0) 
     +                  call fluxi(NDmu, NDsup, NDgeom
     +                   , xmu, smf, imu, imu0
     +                   , ngeom, nmug, nmat, nmutot, ebbot
     +                   , xmut, Tf, Tfs, Rf, Rfs, xm
     +                   , Rmbot,  Tmbot,  Rmsbot,  Wms
     +                   , Rm1bot, Tm1bot, Rm1sbot, Wm1s
     +                   , Riflux, Diflux, Tiflux, URUi, UDUi, UTUi )
              call newfoui(Ndmu, NDsup, NDgeom
     +                    , m, Rmbot, Tmbot, xRm1i, xDm1i
     +                    , xmu, imu, imu0, smf, phi
     +                    , nmat, ngeom
     +                    , Ri, Di, R1i, D1i )
          endif
          if (nextm) goto 1000
*----------------------------------------------------------------------*
*  End of Fourier loop                                                 *
*----------------------------------------------------------------------*
      return
      end subroutine adding

      subroutine test_constant(NDmu, NDsup, NDgeom, nmat, nmug)
*----------------------------------------------------------------------*
*     MC
*      This routine is here to ensure that the 
*      constants used to set-up array dimensions
*      are consistent between them.
*     Ex:
*      - NDgeom = 589 = (31 thetav * 19 phi)
*      - nmug   = 32 (nombre de points d'integration)
*      - NDmu   = 625 > (NDgeom + nmug)
*      - NDsup  = 1875 > (NDmu * nmat=3) 
*----------------------------------------------------------------------*
      if (NDmu.le.(NDgeom+nmug)) then
         write(*,*) ''
         write(*,*) '!!!           program STOP                !!!'
         write(*,*) '!!! NDmu MUST be greater than NDgeom+nmug !!!'
         write(*,*) ''
         STOP
      endif
      if (NDsup.le.(NDmu*nmat)) then
         write(*,*) ''
         write(*,*) '!!!           program STOP               !!!'
         write(*,*) '!!! NDsup MUST be greater than NDmu*nmat !!!'
         write(*,*) ''
         STOP
      endif
      end
      subroutine addi(NDmu, NDsup
     +               , m, nmutot, nmug, nmat, smf, Rf, xmu
     +               , Rm,  Tm,  Rms, Wms, ebmu
     +               , Rm1, Tm1, Rm1s, Wm1s )
*----------------------------------------------------------------------*
*  Add atmosphere and water-air interface, using the algorithm         *
*  described by Chowdhary (1991) section 3.3. Equation numbers in this *
*  subroutine refer to Chowdhary (1991).                               *
*  Only the m-th Fourier component is considered.                      *
*  On entry:                                                           *
*     m       : Fourier index                                          *
*     nmutot  : total number of mu values                              *
*     nmug    : number of Gauss points for integration over mu         *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     smf     : supermatrix factors dsqrt(2*w*mu), or 1 for extra pts. *
*     Rf      : Fresnel matrix given by Eq. (22), coded as follows     *
*                  1,1 element of RF(mu) in Rf(mu,1)                   *
*                  1,2 element of RF(mu) in Rf(mu,2)                   *
*                  3,3 element of RF(mu) in Rf(mu,3)                   *
*               other elements can be found via symmetries of Eq. (22) *
*     xmu     : all different mu values (integration and extra points) *
*     Rm      : Fourier component of reflection supermatrix for the    *
*               atmosphere alone                                       *
*     Tm      : Fourier component of transmission supermatrix for the  *
*               atmosphere alone                                       *
*     Rms     : same as Rm, but for illuminatiojn from below           *
*     ebmu    : dexp(-b/mu) for all mu values                          *
*     Rm1     : first order contribution to Rm                         *
*     Tm1     : first order contribution to Tm                         *
*     Rm1s    : first order contribution to Rms                        *
*     Tm1s    : first order contribution to Tms                        *
*  On exit:                                                            *
*     Rm      : Fourier component of reflection supermatrix for the    *
*               atmosphere-interface system                            *
*     Tm      : Fourier component of transmission supermatrix for the  *
*               downward radiation just above the interface            *
*     Rms     : supermatrix Vm* given by Eq. (90)                      *
*     Wms     : supermatrix Wm* given by Eq. (91)                      *
*     Rm1     : first order contribution to Rm given by Eq. (82)       *
*     Tm1     : first order contribution to Dm given by Eq. (79)       *
*     Rm1s    : unchanged! (first order reflection by atmosphere alone)*
*     Wm1s    : first order contribution to Wm* given by Eq. (95)      *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter(maxrep=15,trmin=1.D-10)
      dimension Rm(NDsup,NDsup),  Tm(NDsup,NDsup),  Rms(NDsup,NDsup)
     +        , Rm1(NDsup,NDsup), Tm1(NDsup,NDsup), Rm1s(NDsup,NDsup)
     +        , Wms(NDsup,NDsup), Wm1s(NDsup,NDsup)
     +        , ebmu(NDmu), Etop(NDsup), smf(NDmu), xmu(NDmu)
      dimension X(NDsup,NDsup), Rf(NDmu,3)
      logical verbo
      verbo = .false.
*----------------------------------------------------------------------*
*  calculate Qm by Eqs. (47)-(49) and put the result into X            *
*  use the product method described by De Haan et al. (1987)           *
*----------------------------------------------------------------------*
      call rRprod(NDmu, NDsup,  Wms, Rms, Rf, xmu, nmat, nmutot )
*                                                    Wms = R'*R" = C1
      call assign( NDsup,  Wm1s, Wms, nmat,nmutot)
*                                                    Wm1s = C1 = S1
      ir = 0
  100     ir = ir+1
          call prod( NDsup,  X, Wms, Wms, nmat,nmutot,nmug)
*                                                    X = Cr Cr = Cr+1
          call assign( NDsup,  Wms, X, nmat,nmutot)
*                                                    Wms = Cr+1
          call prod( NDsup,  X, Wm1s, Wms, nmat,nmutot,nmug)
*                                                    X = Sr Cr+1
          call addSM( NDsup, Wm1s, Wm1s, X, nmat,nmutot)
*                                                    Wm1s = Sr + Sr Cr+1
          call addSM( NDsup, Wm1s, Wm1s, Wms, nmat,nmutot)
*                                      Wm1s = Sr + Sr Cr+1 + Cr+1 = Sr+1
          call trace(NDsup, Wms, trC, nmat,nmug)
*                                                    trC = trace(Cr+1)
          if (verbo) print *,' addi: r = ',ir,' trace = ',trC
          if ((trC .gt. trmin) .and. (ir .le. maxrep)) goto 100
*----------------------------------------------------------------------*
      if (ir .gt. maxrep) then
          write(*,*) ' addi: WARNING repeated reflections did not'
     +              ,' converge after',maxrep,' steps'
          write(*,*) '         proceed anyway !'
      endif
      call assign( NDsup,  X, Wm1s, nmat,nmutot)
*                                                    X    = Q
*----------------------------------------------------------------------*
*  Create the diagonal matrix Tdirect, Eq. (61), so that it can be     *
*  used in the supermatrix formalism to deal with mu integrations      *
*  this means dropping the factor 2mu in denominator and the delta     *
*  function                                                            *
*----------------------------------------------------------------------*
      do 200 i=1, nmutot
           do 150 k=1, nmat
               Etop((i-1)*nmat+k) = ebmu(i)
  150      continue
  200 continue
*----------------------------------------------------------------------*
*  Calculate Dm from Eq. (71) and put it into Wm1s                     *
*----------------------------------------------------------------------*
      call prod( NDsup,  Wm1s, X, Tm, nmat,nmutot,nmug)
*                                                    Wm1s = QT' 
      call addSM( NDsup, Wm1s, Wm1s, Tm, nmat,nmutot)
*                                                    Wm1s = QT' + T' 
      call rdiapr(NDsup, Wms, X, Etop, nmat,nmutot)
*                                                    Wms  = QE
      call addSM( NDsup, Wm1s, Wm1s, Wms, nmat,nmutot)
*                                                    Wm1s = D 
*----------------------------------------------------------------------*
*  Find Tms by symmetry : T* = q3 T~ q3 (also inhomogeneous!)          *
*  and put the result into Tm                                          *
*----------------------------------------------------------------------*
      call Tstar(NDsup,  Wms, Tm, nmat, nmutot )
      call assign( NDsup,  Tm, Wms, nmat,nmutot)
*----------------------------------------------------------------------*
*  See Chowdhary (1991) section 3 for the exception of m=0             *
*  Calculate Q* from Eqs. (54)-(55) and put it into X                  *
*----------------------------------------------------------------------*
      if ( m .eq. 0 ) then
          call prod( NDsup,  Wms, X, Rms, nmat,nmutot,nmug)
*                                                     Wms  = QR'*
          call addSM( NDsup, Wms, Wms, Rms, nmat,nmutot) 
*                                                     Wms  = R'* + QR'*
          call lRprod(NDmu, NDsup,  X, Rf, Wms, xmu, nmat, nmutot )
*                                                     X    = Q* 
*----------------------------------------------------------------------*
*  Calculate Vms from Eq. (90) and put it into Rms                     *
*----------------------------------------------------------------------*
          call prod( NDsup,  Wms, Rms, X, nmat,nmutot,nmug)
*                                                     Wms  = R'*Q* 
          call addSM( NDsup, Rms, Rms, Wms, nmat,nmutot) 
*                                                     Rms  = V*
*----------------------------------------------------------------------*
*  Calculate Wms from Eq. (91) and put it into Wms                     *
*----------------------------------------------------------------------*
          call prod( NDsup,  Wms, Tm, X, nmat,nmutot,nmug)
*                                                     Wms  = T'*Q*
          call ldiapr(NDsup, X, Etop, X, nmat,nmutot)
*                                                     X    = EQ*
          call addSM( NDsup, Wms, Wms, X, nmat,nmutot) 
*                                                     Wms  = T'*Q* + EQ*
          call addSM( NDsup, Wms, Wms, Tm, nmat,nmutot) 
*                                                     Wms  = W*
      endif
*----------------------------------------------------------------------*
*  Calculate Rm from Eq. (73) and put it into Rm                       *
*----------------------------------------------------------------------*
      call rRprod(NDmu, NDsup,  X, Tm, Rf, xmu, nmat, nmutot )
*                                                     X    = T'*R"
      call assign( NDsup,  Tm, Wm1s, nmat,nmutot)
*                                                     Tm   = D
      call rdiapr(NDsup, Wm1s, X, Etop, nmat,nmutot)
*                                                     Wm1s = T'*R"E
      call addSM( NDsup, Rm, Rm, Wm1s, nmat,nmutot) 
*                                                     Rm   = R' + T'*R"E
      call prod( NDsup,  Wm1s, X, Tm, nmat,nmutot,nmug)
*                                                     Wm1s =T'*R"D =T'*U
      call addSM( NDsup, Rm, Rm, Wm1s, nmat,nmutot) 
*                                                Rm = R' + T'*U + T'*R"E
      call lRprod(NDmu, NDsup,  X, Rf, Tm, xmu, nmat, nmutot )
*                                                     X    = R"D = U
      call ldiapr(NDsup, X, Etop, X, nmat,nmutot)
*                                                     X    = EU
      call addSM( NDsup, Rm, Rm, X, nmat,nmutot) 
*                                                     Rm   = R
*---------------------------------------------------------------------*
*  Calculate Dm1 from Eq. (79) and put it into X                      *
*---------------------------------------------------------------------*
      call rRprod(NDmu, NDsup,  X, Rm1s, Rf, xmu, nmat, nmutot )
*                                                     X    = R1'*R"
      call rdiapr(NDsup, X, X, Etop, nmat,nmutot)
*                                                     X    = R1'*R"E
      call addSM( NDsup, X, Tm1, X, nmat,nmutot) 
*                                                     X    = D1
*---------------------------------------------------------------------*
*  Calculate Rm1 from Eq. (82) and put it into Rm1                    *
*  Note that in Eq. (82) an asterisk is missing on T'1 in the 3rd term*
*---------------------------------------------------------------------*
      call lRprod(NDmu, NDsup,  Wm1s, Rf, X, xmu, nmat, nmutot )
*                                                     Wm1s = R"D1
      call ldiapr(NDsup, Wm1s, Etop, Wm1s, nmat,nmutot)
*                                                     Wm1s = ER"D1
      call addSM( NDsup, Rm1, Rm1, Wm1s, nmat,nmutot) 
*                                                     Rm1  = R1 + ER"D1
*----------------------------------------------------------------------*
*  Find T1ms by symmetry : T1* = q3 T1~ q3 (also inhomogeneous!)       *
*  and put the result into Tm1                                         *
*----------------------------------------------------------------------*
      call Tstar(NDsup,  Wm1s, Tm1, nmat, nmutot )
      call assign( NDsup,  Tm1, Wm1s, nmat,nmutot)
*----------------------------------------------------------------------*
      call rRprod(NDmu, NDsup,  Wm1s, Tm1, Rf, xmu, nmat, nmutot )
*                                                      Wm1s = T1'*R"
      call rdiapr(NDsup, Wm1s, Wm1s, Etop, nmat,nmutot)
*                                                      Wm1s = T1'*R"E
      call addSM( NDsup, Rm1, Rm1, Wm1s, nmat,nmutot) 
*                                                      Rm1  = R1
*----------------------------------------------------------------------*
*  Calculate Wm1s from Eq. (95)                                        *
*----------------------------------------------------------------------*
      if (m .eq. 0) then
          call assign( NDsup,  Wm1s, Tm1, nmat,nmutot)
*                                                      Wm1s = T1'*
          call lRprod(NDmu, NDsup,  Tm1, Rf, Rm1s, xmu, nmat, nmutot )
*                                                      Tm1  = R"R1'*
          call ldiapr(NDsup, Tm1, Etop, Tm1, nmat,nmutot)
*                                                      Tm1  = ER"R1'*
          call addSM( NDsup, Wm1s, Wm1s, Tm1, nmat,nmutot) 
*                                                      Wm1s = W1*
      endif
*----------------------------------------------------------------------*
*  finally, D1 is put into Tm1                                         *
*----------------------------------------------------------------------*
      call assign( NDsup,  Tm1, X, nmat,nmutot)
*                                                      Tm1  = D1
      return
      end
      subroutine addim1( NDmu, NDgeom, xmu, imu, imu0, ngeom, nmat
     +                 , Rf, ebbot
     +                 , xRm1bot, xRm1sbot, xTm1bot, xTm1sbot
     +                 , xRm1i, xDm1i )
*----------------------------------------------------------------------*
*  Add atmosphere to interface for first order scattering, for each    *
*  geometry. Do this only for the m-th Fourier component.              *
*  Equations are given by Chowdhary (1991) Eqs. (79)-(82).             *
*  On entry:                                                           *
*     xmu     : all different mu values (integration and extra points) *
*     imu     : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu value can be found  *
*     imu0    : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu0 value can be found *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*      Rf     : Fresnel matrix given by Eq. (22), coded as follows     *
*                  1,1 element of RF(mu) in Rf(mu,1)                   *
*                  1,2 element of RF(mu) in Rf(mu,2)                   *
*                  3,3 element of RF(mu) in Rf(mu,3)                   *
*               other elements can be found via symmetries of Eq. (22) *
*     ebbot   : dexp(-b/mu) where b is the optical thickness of        *
*               the atmosphere                                         *
*     xRm1bot : first order scattering contribution to the reflection  *
*               matrix of the atmosphere for each geometry             *
*     xRm1sbot: idem, but for illumination from below                  *
*     xTm1bot : first order scattering contribution to the transission *
*               matrix of the atmosphere for each geometry             *
*     xTm1sbot: idem, but for illumination from below                  *
*  On exit:                                                            *
*     xRm1i   : first order contribution to the reflection matrix of   *
*               the atmosphere-interface combination for each geometry *
*     xDm1i   : first order contribution to the downward radiation     *
*               matrix just above the interface for each geometry      *
*----------------------------------------------------------------------*
      implicit double precision(a-h,o-z)
      dimension xmu(NDmu), imu(NDgeom), imu0(NDgeom)
     +        , Rf(NDmu,3), ebbot(NDmu)
     +        , xRm1bot(4,4,NDgeom), xRm1sbot(4,4,NDgeom)
     +        , xTm1bot(4,4,NDgeom), xTm1sbot(4,4,NDgeom)
      dimension xRm1i(4,4,NDgeom), xDm1i(4,4,NDgeom)
      dimension a(NDgeom), a0(NDgeom)
     +        , b(NDgeom), b0(NDgeom)
     +        , c(NDgeom), c0(NDgeom), X(4,4)
      do 200 i=1, ngeom
          mu = imu(i)
          mu0= imu0(i)
          ebmu = ebbot(mu)
          ebmu0= ebbot(mu0)
*----------------------------------------------------------------------*
*  We use Eqs. (79)-(82). Using Eqs. (38) and (61) we obtain from      *
*  Eq. (79):                                                           *
*                                                                      *
*  Dm1(mu,mu0) = Tm1'(mu,mu0) + Rm1'*(mu,mu0) 2 mu0 Rf(mu0) ebmu0      *
*                                                                      *
*                                                 ( a0 b0 0  0 )       *
*                  = Tm1'(mu,mu0) + Rm1'*(mu,mu0) ( b0 a0 0  0 )       *
*                                                 ( 0  0  c0 0 )       *
*                                                 ( 0  0  0 c0 )       *
*                                                                      *
*  There is an error in Eq. (82): T1' must be T1'* .                   *
*  Using Eqs. (38), (61), (80) and (82) we find:                       *
*                                                                      *
*                               ( a  b  0  0 )                         *
*  Rm1(mu,mu0) = Rm1'(mu,mu0) + ( b  a  0  0 ) Dm1(mu,mu0) + X(mu,mu0) *
*                               ( 0  0  c  0 )                         *
*                               ( 0  0  0  c )                         *
*                                                                      *
*  where a, b and c are elements of   2 mu Rf(mu) ebmu , and           *
*                                                                      *
*  X(mu,mu0) = Tm1'*(mu,mu0) 2 mu0 Rf(mu0) ebmu0                       *
*----------------------------------------------------------------------*
          a0(i) = 2.d0*xmu(imu0(i))*ebmu0*Rf(imu0(i),1)
          b0(i) = 2.d0*xmu(imu0(i))*ebmu0*Rf(imu0(i),2)
          c0(i) = 2.d0*xmu(imu0(i))*ebmu0*Rf(imu0(i),3)
          a(i)  = 2.d0*xmu(imu(i)) *ebmu *Rf(imu(i), 1)
          b(i)  = 2.d0*xmu(imu(i)) *ebmu *Rf(imu(i), 2)
          c(i)  = 2.d0*xmu(imu(i)) *ebmu *Rf(imu(i), 3)
  200 continue
*----------------------------------------------------------------------*
*  Handle case without polarization                                    *
*----------------------------------------------------------------------*
      if (nmat.eq.1) then
          do 300 i=1, ngeom
            xDm1i(1,1,i) = xRm1sbot(1,1,i)*a0(i) + xTm1bot(1,1,i)
            xRm1i(1,1,i) = xRm1bot(1,1,i) + a(i)*xDm1i(1,1,i) 
     +                                   + xTm1sbot(1,1,i)*a0(i)
  300     continue
*----------------------------------------------------------------------*
*  Handle 3X3 approximation                                            *
*----------------------------------------------------------------------*
      else if (nmat.eq.3) then
        do 800 i=1, ngeom
          do 400 k=1, nmat
            xDm1i(k,1,i) = xRm1sbot(k,1,i)*a0(i) + xRm1sbot(k,2,i)*b0(i)
            xDm1i(k,2,i) = xRm1sbot(k,1,i)*b0(i) + xRm1sbot(k,2,i)*a0(i)
            xDm1i(k,3,i) = xRm1sbot(k,3,i)*c0(i)
            X(k,1)     = xTm1sbot(k,1,i)*a0(i) + xTm1sbot(k,2,i)*b0(i)
            X(k,2)     = xTm1sbot(k,1,i)*b0(i) + xTm1sbot(k,2,i)*a0(i)
            X(k,3)     = xTm1sbot(k,3,i)*c0(i)
  400     continue
          do 600 l=1, nmat
            do 500 k=1, nmat
              xDm1i(k,l,i) = xDm1i(k,l,i) + xTm1bot(k,l,i)
              xRm1i(k,l,i) = xRm1bot(k,l,i) + X(k,l)
  500       continue
  600     continue
          do 700 l=1, nmat
            xRm1i(1,l,i) = xRm1i(1,l,i) + a(i)*xDm1i(1,l,i) 
     +                                  + b(i)*xDm1i(2,l,i)
            xRm1i(2,l,i) = xRm1i(2,l,i) + a(i)*xDm1i(2,l,i) 
     +                                  + b(i)*xDm1i(1,l,i)
            xRm1i(3,l,i) = xRm1i(3,l,i) + c(i)*xDm1i(3,l,i)
  700     continue
  800   continue
*----------------------------------------------------------------------*
*  Handle case with polarization fully included                        *
*----------------------------------------------------------------------*
      else if (nmat.eq.4) then
        do 1050 i=1, ngeom
          do 850 k=1, nmat
            xDm1i(k,1,i) = xRm1sbot(k,1,i)*a0(i) + xRm1sbot(k,2,i)*b0(i)
            xDm1i(k,2,i) = xRm1sbot(k,1,i)*b0(i) + xRm1sbot(k,2,i)*a0(i)
            xDm1i(k,3,i) = xRm1sbot(k,3,i)*c0(i)
            xDm1i(k,4,i) = xRm1sbot(k,4,i)*c0(i)
            X(k,1)    = xTm1sbot(k,1,i)*a0(i) + xTm1sbot(k,2,i)*b0(i)
            X(k,2)    = xTm1sbot(k,1,i)*b0(i) + xTm1sbot(k,2,i)*a0(i)
            X(k,3)    = xTm1sbot(k,3,i)*c0(i)
            X(k,4)    = xTm1sbot(k,4,i)*c0(i)
  850     continue
          do 950 l=1, nmat
            do 900 k=1, nmat
              xDm1i(k,l,i) = xDm1i(k,l,i) + xTm1bot(k,l,i)
              xRm1i(k,l,i) = xRm1bot(k,l,i) + X(k,l)
  900       continue
  950     continue
          do 1000 l=1, nmat
            xRm1i(1,l,i) = xRm1i(1,l,i) + a(i)*xDm1i(1,l,i) 
     +                                  + b(i)*xDm1i(2,l,i)
            xRm1i(2,l,i) = xRm1i(2,l,i) + a(i)*xDm1i(2,l,i) 
     +                                  + b(i)*xDm1i(1,l,i)
            xRm1i(3,l,i) = xRm1i(3,l,i) + c(i)*xDm1i(3,l,i)
            xRm1i(4,l,i) = xRm1i(4,l,i) + c(i)*xDm1i(4,l,i)
 1000     continue
 1050   continue
      else
          write(*,*) ' addim1: illegal value for nmat = ',nmat
          stop 'in addim1 illegal value for nmat'
      endif
      return
      end
      subroutine addlay( NDmu, NDsup, Rmtop, Tmtop, Rmbot, Tmbot, Rmsbot
     +                 , Rm1top, Tm1top, Rm1bot, Rm1sbot, Tm1bot
     +                 , ebtop, ebbot, iadd
     +                 , nmutot, nmug, nmat )
*----------------------------------------------------------------------*
*  Use the adding method to calculate reflection and transmission of   *
*  of the combination of top and bottom layer. The top layer is        *
*  assumed homogeneous. The resulting reflection and transmission      *
*  supermatrices are returned through arrays Rmbot, Tmbot and Rmsbot.  *
*  BEWARE:  Rmtop and Tmtop are used as scratch space !                *
*  The option variable iadd indicates the following cases :            *
*     iadd = 1 : normal adding                                         *
*     iadd = 2 : no scattering in top layer                            *
*     iadd = 3 : no scattering in bottom layer                         *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter( maxrep=15, trmin=1.D-10 )
      dimension W(NDsup,NDsup), X(NDsup,NDsup)
     +        , Y(NDsup,NDsup), Z(NDsup,NDsup)
     +        , Rmtop(NDsup,NDsup), Tmtop(NDsup,NDsup)
     +        , Rmbot(NDsup,NDsup), Tmbot(NDsup,NDsup)
     +        , Rmsbot(NDsup,NDsup)
     +        , Rm1top(NDsup,NDsup), Tm1top(NDsup,NDsup)
     +        , Rm1bot(NDsup,NDsup), Tm1bot(NDsup,NDsup)
     +        , Rm1sbot(NDsup,NDsup)
     +        , Etop(NDsup), Ebot(NDsup), ebtop(NDmu), ebbot(NDmu)
      logical verbo
      verbo = .false.
*----------------------------------------------------------------------*
*  Handle the cases of no scattering in top or bottom layer            *
*----------------------------------------------------------------------*
      if (iadd .eq. 2) then
          call notop(NDmu, NDsup, Rmbot, Tmbot, Rm1bot, Tm1bot
     +              , ebtop, ebbot, nmutot, nmat )
          goto 999
      else if (iadd .eq. 3) then
          call nobot(NDmu, NDsup, Rmtop, Tmtop, Rmbot, Tmbot, Rmsbot
     +               , Rm1top, Tm1top, Rm1bot, Rm1sbot, Tm1bot
     +               , ebtop, ebbot, nmutot, nmat )
          goto 999
      endif
*----------------------------------------------------------------------*
*  Use the product method to calculate the repeated reflections        *
*  between the layers, se de Haan et al. (1987) Eqs. (111)-(115)       *
*----------------------------------------------------------------------*
      call star(NDsup,  X, Rmtop, nmat,nmutot)
*                                                    X = R'*
      call prod( NDsup,  Y, X, Rmbot, nmat,nmutot,nmug)
*                                                    Y = R'*R" = C1
      call assign( NDsup,  Z, Y, nmat,nmutot)
*                                                    Z = C1 = S1
      ir = 0
  100     ir = ir+1
          call prod( NDsup,  X, Y, Y, nmat,nmutot,nmug)
*                                                    X = Cr Cr = Cr+1
          call assign( NDsup,  Y, X, nmat,nmutot)
*                                                    Y = Cr+1
          call prod( NDsup,  X, Z, Y, nmat,nmutot,nmug)
*                                                    X = Sr Cr+1
          call addSM( NDsup, Z, Z, X, nmat,nmutot)
*                                                    Z = Sr + Sr Cr+1
          call addSM( NDsup, Z, Z, Y, nmat,nmutot)
*                                         Z = Sr + Sr Cr+1 + Cr+1 = Sr+1
          call trace(NDsup, Y, trC, nmat,nmug)
*                                                    trC = trace(Cr+1)
          if (verbo) print *,' addlay: r = ',ir,' trace = ',trC
          if ((trC .gt. trmin) .and. (ir .le. maxrep)) goto 100
*----------------------------------------------------------------------*
      if (ir .gt. maxrep) then
          print *,' addlay: WARNING repeated reflections did not'
     +           ,' converge after',maxrep,' steps'
          print *,'         proceed anyway !'
      endif
*----------------------------------------------------------------------*
*  Now Z contains the matrix Q Eq. (115)                               *
*  Use the adding Eqs. (85)-(91)                                       *
*----------------------------------------------------------------------*
      do 200 i=1, nmutot
           do 150 k=1, nmat
               Etop((i-1)*nmat+k) = ebtop(i)
               Ebot((i-1)*nmat+k) = ebbot(i)
  150     continue
  200 continue
      call prod( NDsup,  X, Z, Tmtop, nmat,nmutot,nmug)
*                                              X = QT'
      call rdiapr(NDsup, W, Z, Etop, nmat,nmutot)
*                                              W = QE'
      call addSM( NDsup, X, X, W, nmat,nmutot)
*                                              X = QT' + QE'
      call addSM( NDsup, X, X, Tmtop, nmat,nmutot)
*                                              X = T' + QT' +QE' = D
*----------------------------------------------------------------------*
      call prod( NDsup,  W, Rmbot, X, nmat,nmutot,nmug)
*                                              W = R"D
      call rdiapr(NDsup, Y, Rmbot, Etop, nmat,nmutot)
*                                              Y = R"E'
      call addSM( NDsup, W, W, Y, nmat,nmutot)
*                                              W = R"E' + R"D' = U
*----------------------------------------------------------------------*
      call prod( NDsup,  Y, Tmbot, X, nmat,nmutot,nmug)
*                                              Y = T"D
      call ldiapr(NDsup, X, Ebot, X, nmat,nmutot)
*                                              X = E"D
      call addSM( NDsup, Y, Y, X, nmat,nmutot)
*                                              Y = T"D + E"D
      call rdiapr(NDsup, X, Tmbot, Etop, nmat,nmutot)
*                                              X = T"E'
      call addSM( NDsup, Y, Y, X, nmat,nmutot)
*                                             Y = T"D + E"D + T"E' =Ttot
*----------------------------------------------------------------------*
      call star(NDsup,  Tmtop, Tmtop, nmat,nmutot)
*                                              Tmtop = T'*
      call prod( NDsup,  X, Tmtop, W, nmat,nmutot,nmug)
*                                              X = T'*U
      call addSM( NDsup, Rmbot, Rmtop, X, nmat,nmutot)
*                                              Rmbot = R' + T'*U
      call ldiapr(NDsup, W, Etop, W, nmat,nmutot)
*                                              W = E'U
      call addSM( NDsup, Rmbot, Rmbot, W, nmat,nmutot)
*                                         Rmbot = R' + T'*U + E'U = Rtot
*----------------------------------------------------------------------*
      call star(NDsup,  Rmtop, Rmtop, nmat,nmutot)
*                                              Rmtop = R'*
      call prod( NDsup,  X, Z, Rmtop, nmat,nmutot,nmug)
*                                              X = QR'*
      call prod( NDsup,  Z, X, Tmbot, nmat,nmutot,nmug)
*                                              Z = QR'*T"*
      call rdiapr(NDsup, X, X, Ebot, nmat,nmutot)
*                                              X = QR'*E"
      call addSM( NDsup, X, X, Z, nmat,nmutot)
*                                              X = QR'*E" + QR'*T"*
      call Tstar(NDsup,  W, Tmbot, nmat,nmutot)
*                                              W = T"*
      call prod( NDsup,  Z, Rmtop, W, nmat,nmutot,nmug)
*                                              Z = R'*T"*
      call addSM( NDsup, X, X, Z, nmat,nmutot)
*                                          X = QR'*E" + QR'*T"* + R'*T"*
      call rdiapr(NDsup, Z, Rmtop, Ebot, nmat,nmutot)
*                                              Z = R'*E"
      call addSM( NDsup, X, X, Z, nmat,nmutot)
*                             X = QR'*E" + QR'*T"* + R'*T"* + R'*E" = U*
*----------------------------------------------------------------------*
      call prod( NDsup,  Z, Tmbot, X, nmat,nmutot,nmug)
*                                              Z = T"U*
      call addSM( NDsup, Rmsbot, Rmsbot, Z, nmat,nmutot)
*                                              Rmsbot = R"* + T"U*
      call ldiapr(NDsup, X, Ebot, X, nmat,nmutot)
*                                              X = E"U*
      call addSM( NDsup, Rmsbot, Rmsbot, X, nmat,nmutot)
*                                     Rmsbot = R"* + T"U* + E"U* = Rtot*
*----------------------------------------------------------------------*
      call assign( NDsup,  Tmbot, Y, nmat,nmutot)
*                                        Tmbot = T"D + E"D + T"E' =Ttot
*----------------------------------------------------------------------*
*  Use adding equations for first and zero order scattering            *
*----------------------------------------------------------------------*
      call oneadd( NDmu, NDsup, Rm1top, Tm1top, Rm1bot, Rm1sbot, Tm1bot
     +           , nmat, nmutot, ebtop, ebbot )
  999 return
      end
      subroutine addm1( NDmu, NDgeom, ngeom, imu, imu0, nmat
     +                , xRm1top, xRm1stop, xTm1top, xTm1stop
     +                , ebtop, ebbot
     +                , xRm1bot, xRm1sbot, xTm1bot, xTm1sbot )
*----------------------------------------------------------------------*
*  Use 'adding' equations for first order only for each geometry.      *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension imu(NDgeom), imu0(NDgeom)
      dimension xRm1top(4,4,NDgeom), xRm1stop(4,4,NDgeom)
     +        , xTm1top(4,4,NDgeom), xTm1stop(4,4,NDgeom)
     +        , ebtop(NDmu), ebbot(NDmu)
     +        , xRm1bot(4,4,NDgeom), xRm1sbot(4,4,NDgeom)
     +        , xTm1bot(4,4,NDgeom), xTm1sbot(4,4,NDgeom)
      do 1000 igeom=1, ngeom
          mu = imu(igeom)
          mu0= imu0(igeom)
          do 200 k=1, nmat
              do 100 l=1, nmat
                  xRm1bot(l,k,igeom) = xRm1top(l,k,igeom) 
     +                       + ebtop(mu)*xRm1bot(l,k,igeom)*ebtop(mu0)
                  xTm1bot(l,k,igeom) = xTm1bot(l,k,igeom)*ebtop(mu0) 
     +                       + ebbot(mu)*xTm1top(l,k,igeom)
                  xRm1sbot(l,k,igeom) = xRm1sbot(l,k,igeom) 
     +                       + ebbot(mu)*xRm1stop(l,k,igeom)*ebbot(mu0)
                  xTm1sbot(l,k,igeom) = ebtop(mu)*xTm1sbot(l,k,igeom)
     +                       + xTm1stop(l,k,igeom)*ebbot(mu0)
  100         continue
  200     continue
 1000 continue
      return
      end
      subroutine brack( el, array, n, ND, i1, i2 )
*----------------------------------------------------------------------*
*  Find a bracket around el in the first n elements of array which     *
*  must already be sorted in increasing order.                         *
*  On entry:                                                           *
*     el      : element to be bracketed                                *
*     array   : array in which bracket must be found                   *
*     n       : number of elements in array to be considered           *
*     ND      : dimension of array                                     *
*  On exit:                                                            *
*     i1, i2  : indices of bracket elements                            *
*               if el is outside the range of array, i1 and i2         *
*               are set to the appropriate extreme index value         *
*----------------------------------------------------------------------*
      implicit double precision(a-h,o-z)
      dimension array(ND)
      if (el .le. array(1)) then
          i1 = 1
          i2 = 1
      else if (el .gt. array(n)) then
          i1 = n
          i2 = n
      else
          do 100 i=2, n
              if ((el .le. array(i)) .and. (el .gt. array(i-1))) then
                  i2 = i
                  i1 = i - 1
              endif
  100     continue
      endif
      return
      end
      subroutine bstart( NDlay, NDcoef
     +                 , m, layer, coefs, ncoef, M0, xmumin
     +                 , a, b, epsilon, b0, ndoubl )
*----------------------------------------------------------------------*
*  Calculate the optical thickness b0 at which doubling should be      *
*  started to obtain an error less than epsilon in the total layer     *
*  with thickness b. This is done for the m-th Fourier component.      *
*  The values of b0 and ndoubl returned are such that b = b0*2**ndoubl *
*  The algorithm is described in de Haan et al. (1987).                *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension coefs(4,4,0:NDcoef,NDlay), ncoef(NDlay), M0(NDlay)
      logical verbo
      verbo = .false.
      if (m .gt. M0(layer)) then
          print *,' bstart: m = ',m,' larger than M0 = ',M0(layer)
          stop 'in bstart m larger than M0'
      endif
      if (ncoef(layer) .lt. M0(layer)) then
          print *,' bstart: ncoef =',ncoef(layer)
     +           ,' larger than M0 =',M0(layer)
          stop 'in bstart ncoef larger than M0'
      endif
*----------------------------------------------------------------------*
*  Calculate the effective albedo am, Eq. (143)                        *
*----------------------------------------------------------------------*
      am = 0.D0
      do 100 k=m, M0(layer)
          amtry = a*dabs(coefs(1,1,k,layer))/(2.D0*k+1.D0)
          if (am .lt. amtry) am = amtry
  100 continue
*----------------------------------------------------------------------*
*  Right hand side of Eq. (142)                                        *
*----------------------------------------------------------------------*
      rhs = 4.D0*epsilon/(9.D0*b*am**3*dble(2*m+1))
      ndoubl = -1
      b0     = 2.D0*b
  200     ndoubl = ndoubl+1
          b0     = 0.5D0*b0
          fb0mu  = b0/xmumin
          if (fb0mu .gt. 1.D0) fb0mu = 1.D0
          if (b0 .ge. (rhs/fb0mu)) goto 200
          if (b0 .ge. (1.D0/3.D0)) goto 200
      if (verbo) print *,' bstart: start optical thickness = ',b0
     +                                                ,' ndoubl=',ndoubl
      return
      end
      subroutine dirmat( NDmu,NDgeom
     +                 , Aflux, igeom, muarg, illum, mu, nmat, fac, Af )
*----------------------------------------------------------------------*
*  Put direct flux fac * Af into array Aflux.                          *
*  Use symmetries:      ( Af(1)  Af(2)   0      0    )                 *
*                       ( Af(2)  Af(1)   0      0    )                 *
*                       (  0      0     Af(3)   0    )                 *
*                       (  0      0      0     Af(3) )                 *
*  Assume that Aflux has been initialized to zero.                     *
*----------------------------------------------------------------------*
      implicit double precision(a-h,o-z)
      parameter( Nmu=2
     +         , Nzero=1,  Nall=3
     +         , Nbelow=2 )
      dimension Aflux(4,4, NDgeom, Nmu, Nbelow, Nall), Af(NDmu,3)
      if (nmat .eq. 1) then
          Aflux(1,1,igeom,muarg,illum,Nzero) = fac*Af(mu,1)
      else if (nmat .eq. 3) then
          Aflux(1,1,igeom,muarg,illum,Nzero) = fac*Af(mu,1)
          Aflux(2,2,igeom,muarg,illum,Nzero) = fac*Af(mu,1)
          Aflux(3,3,igeom,muarg,illum,Nzero) = fac*Af(mu,3)
          Aflux(1,2,igeom,muarg,illum,Nzero) = fac*Af(mu,2)
          Aflux(2,1,igeom,muarg,illum,Nzero) = fac*Af(mu,2)
      else if (nmat .eq. 4) then
          Aflux(1,1,igeom,muarg,illum,Nzero) = fac*Af(mu,1)
          Aflux(2,2,igeom,muarg,illum,Nzero) = fac*Af(mu,1)
          Aflux(3,3,igeom,muarg,illum,Nzero) = fac*Af(mu,3)
          Aflux(4,4,igeom,muarg,illum,Nzero) = fac*Af(mu,3)
          Aflux(1,2,igeom,muarg,illum,Nzero) = fac*Af(mu,2)
          Aflux(2,1,igeom,muarg,illum,Nzero) = fac*Af(mu,2)
      else
          write(*,*) ' dirmat: illlegal value for nmat=',nmat
          stop 'in dirmat illegal value for nmat'
      endif
      return
      end
      subroutine double(NDmu, NDsup, Rm, Tm, ebmu, nmutot, nmug, nmat )
*----------------------------------------------------------------------*
*  Calculate the m-th Fourier term of reflection ans transmission of   *
*  a homogeneous layer from the reflection and transmission of a layer *
*  with only half the optical thickness.                               *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter( maxrep=15, trmin=1.D-8 )
      dimension X(NDsup,NDsup), Y(NDsup,NDsup), Z(NDsup,NDsup)
     +        , Rm(NDsup,NDsup), Tm(NDsup,NDsup), ebmu(NDmu)
     +        , E(NDsup)
      logical verbo
      verbo = .false.
*----------------------------------------------------------------------*
*  Use the product method to calculate the repeated reflections        *
*  between the layers, se de Haan et al. (1987) Eqs. (111)-(115)       *
*----------------------------------------------------------------------*
      call star(NDsup,  X, Rm, nmat,nmutot)
*                                                    X = R*
      call prod( NDsup,  Y, X, Rm, nmat,nmutot,nmug)
*                                                    Y = R*R = C1
      call assign( NDsup,  Z, Y, nmat,nmutot)
*                                                    Z = C1 = S1
      ir = 0
  100     ir = ir+1
          call prod( NDsup,  X, Y, Y, nmat,nmutot,nmug)
*                                                    X = Cr Cr = Cr+1
          call assign( NDsup,  Y, X, nmat,nmutot)
*                                                    Y = Cr+1
          call prod( NDsup,  X, Z, Y, nmat,nmutot,nmug)
*                                                    X = Sr Cr+1
          call addSM( NDsup, Z, Z, X, nmat,nmutot)
*                                                    Z = Sr + Sr Cr+1
          call addSM( NDsup, Z, Z, Y, nmat,nmutot)
*                                         Z = Sr + Sr Cr+1 + Cr+1 = Sr+1
          call trace(NDsup, Y, trC, nmat,nmug)
*                                                    trC = trace(Cr+1)
          if (verbo) print *,' double: r = ',ir,' trace = ',trC
          if ((trC .gt. trmin) .and. (ir .le. maxrep)) goto 100
*----------------------------------------------------------------------*
      if (ir .gt. maxrep) then
          print *,' double: WARNING repeated reflections did not'
     +           ,' converge after',maxrep,' steps'
          print *,'         proceed anyway !'
      endif
*----------------------------------------------------------------------*
*  Now Z contains the matrix Q Eq. (115)                               *
*  Use the adding Eqs. (85)-(91) with identical layers                 *
*----------------------------------------------------------------------*
      do 200 i=1, nmutot
           do 150 k=1, nmat
               E((i-1)*nmat+k) = ebmu(i)
  150     continue
  200 continue
      call prod( NDsup,  X, Z, Tm, nmat,nmutot,nmug)
*                                              X = QT
      call rdiapr(NDsup, Z, Z, E, nmat,nmutot)
*                                              Z = QE
      call addSM( NDsup, X, X, Z, nmat,nmutot)
*                                              X = QT + QE
      call addSM( NDsup, X, X, Tm, nmat,nmutot)
*                                              X = T + QT +QE = D
*----------------------------------------------------------------------*
      call prod( NDsup,  Z, Rm, X, nmat,nmutot,nmug)
*                                              Z = RD
      call rdiapr(NDsup, Y, Rm, E, nmat,nmutot)
*                                              Y = RE
      call addSM( NDsup, Z, Z, Y, nmat,nmutot)
*                                              Z = RE + RD = U
*----------------------------------------------------------------------*
      call prod( NDsup,  Y, Tm, X, nmat,nmutot,nmug)
*                                              Y = TD
      call ldiapr(NDsup, X, E, X, nmat,nmutot)
*                                              X = ED
      call addSM( NDsup, Y, Y, X, nmat,nmutot)
*                                              Y = TD + ED
      call rdiapr(NDsup, X, Tm, E, nmat,nmutot)
*                                              X = TE
      call addSM( NDsup, Y, Y, X, nmat,nmutot)
*                                              Y = TD + ED + TE = Ttot
*----------------------------------------------------------------------*
      call star(NDsup,  Tm, Tm, nmat,nmutot)
*                                              Tm = T*
      call prod( NDsup,  X, Tm, Z, nmat,nmutot,nmug)
*                                              X = T*U
      call addSM( NDsup, Rm, Rm, X, nmat,nmutot)
*                                              Rm = R + T*U
      call ldiapr(NDsup, Z, E, Z, nmat,nmutot)
*                                              Z = EU
      call addSM( NDsup, Rm, Rm, Z, nmat,nmutot)
*                                              Rm = R + T*U + EU = Rtot
*----------------------------------------------------------------------*
      call assign( NDsup,  Tm, Y, nmat,nmutot)
*                                              Tm = Ttot
      return
      end
      subroutine endfou( NDmu, NDsup, NDlay, NDgeom
     +                 , m, M1, M0, nfoumax, nlayer, epsilon, mstop
     +                 , xmu, imu, imu0, nmutot, nmat, ngeom
     +                 , Rm1bot, Rmbot, Tm1bot, Tmbot, nextm )
*----------------------------------------------------------------------*
*  Decide if the next Fourier term (m+1) is needed, if so return       *
*  nextm = .true. else return nextm = .false.                          *
*  Because first order is treated separately without Fourier expansion *
*  it is enough to stop when only first order is needed.               *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension M1(NDlay), M0(NDlay)
     +        , xmu(NDmu), imu(NDgeom), imu0(NDgeom)
     +        , Rm1bot(NDsup,NDsup), Rmbot(NDsup,NDsup)
     +        , Tm1bot(NDsup,NDsup), Tmbot(NDsup,NDsup)
      logical nextm, verbo, vertic, almost, except
      save almost
      verbo = .false.
      except = .true.
      if (m .eq. 0) almost = .false.
      nextm = .true.
*----------------------------------------------------------------------*
*  First test whether the user supplied bound nfoumax has been reached *
*----------------------------------------------------------------------*
      if (m .ge. nfoumax) then
          nextm = .false.
          if (verbo) print *,' endfou: stop Fourier series after m = '
     +                      ,m,' (nfoumax reached)'
          goto 999
      endif
*----------------------------------------------------------------------*
*  Test if M1 has been reached (see de Haan et al. (1987) section 7.2) *
*  so that higher Fourier terms will only have single scattering.      *
*  Of course we should take the maximum M1 over all layers.            *
*  If first order scattering is not excepted (except = .false.) then   *
*  we must sum the Fourier series all the way to M0.                   *
*----------------------------------------------------------------------*
      maxM = 0
      if (except) then
          do 100 layer=1, nlayer
              if (maxM .lt. M1(layer)) maxM = M1(layer)
  100     continue
      else
          do 150 layer=1, nlayer
              if (maxM .lt. M0(layer)) maxM = M0(layer)
  150     continue
      endif
      if ( m .ge. maxM) then
          nextm = .false.
          if (verbo) print *,' endfou: stop Fourier series after m = '
     +                      ,m,' (M1 reached)'
          goto 999
      endif
*----------------------------------------------------------------------*
*  Skip other covergence tests when mstop = 0                          *
*----------------------------------------------------------------------*
      if (mstop .eq. 0) goto 999
*----------------------------------------------------------------------*
*  For polarization m = 2 is always needed due to rotations of plane   *
*  of reference                                                        *
*----------------------------------------------------------------------*
      if ((nmat .gt. 1) .and. (m .lt. 2)) then
          nextm = .true.
          goto 999
      endif
*----------------------------------------------------------------------*
*  Check the case of vertical directions which only need m = 0         *
*----------------------------------------------------------------------*
      vertic = .true.
      do 200 igeom=1, ngeom
          if ((dabs(xmu(imu(igeom)) -1.D0) .gt. 1.D-6) .and.
     +        (dabs(xmu(imu0(igeom))-1.D0) .gt. 1.D-6)) vertic = .false.
  200 continue
      if (vertic) then
          nextm = .false.
          if (verbo) print *,' endfou: stop Fourier series after m = '
     +                      ,m,' (vertical)'
          goto 999
      endif
*----------------------------------------------------------------------*
*  Check convergence for extra points                                  *
*----------------------------------------------------------------------*
      nextm = .false.
      if (except) then
        do 300 igeom=1, ngeom
          isup = (imu(igeom)-1)*nmat+1
          jsup = (imu0(igeom)-1)*nmat+1
          if (xmu(imu0(igeom))*dabs(Rmbot(isup,jsup)-Rm1bot(isup,jsup)) 
     +                                                     .gt. epsilon)
     +        nextm = .true.
          if (xmu(imu0(igeom))*dabs(Tmbot(isup,jsup)-Tm1bot(isup,jsup)) 
     +                                                     .gt. epsilon)
     +        nextm = .true.
  300   continue
      else
          isup = (imu(igeom)-1)*nmat+1
          jsup = (imu0(igeom)-1)*nmat+1
          if (xmu(imu0(igeom))*dabs(Rmbot(isup,jsup)) .gt. epsilon)
     +        nextm = .true.
          if (xmu(imu0(igeom))*dabs(Tmbot(isup,jsup)) .gt. epsilon)
     +        nextm = .true.
      endif
*----------------------------------------------------------------------*
*  Logical 'almost' indicates if convergence has already been suspected*
*  for the previous Fourier component. Only if two successive Fourier  *
*  components indicate convergence do we return nextm = .true.         *
*----------------------------------------------------------------------*
      if ((.not. nextm) .and. (.not. almost)) then
          nextm  = .true.
          almost = .true.
      else
          almost = .not. nextm
      endif
      if (verbo .and. .not. nextm) 
     +     print *,' endfou: stop Fourier series after m = '
     +                      ,m,' (for extra points)'
  999 return
      end
      subroutine expbmu(NDmu, b, xmu, nmutot, ebmu )
*----------------------------------------------------------------------*
*  Calculate dexp(-b/mu) for all mu in xmu and put result in array ebmu*
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), ebmu(NDmu)
      do 100 i=1, nmutot
          if (xmu(i) .gt. (b/100.D0)) then
              ebmu(i) = dexp(-b/xmu(i))
          else
              ebmu(i) = 0.D0
          endif
  100 continue
      return
      end
      subroutine fasemx(NDmu, NDlay, NDgeom, NDcoef
     +                 , layer, coefs, ncoefs
     +                 , xmu, imu, imu0, phi, ngeom
     +                 , gsfpl, gsfmi, nmat
     +                 , Zmin, Zplus )
*----------------------------------------------------------------------*
*  Calculate the phase matrix Z(+-mu,mu0,phi) for all geometries       *
*  for a specific layer.                                               *
*  This amounts to summing the expansion in generalized spherical      *
*  functions and rotating the reference planes of incident and         *
*  scattered Stokes vector. The scattering matrix uses the scattering  *
*  plane as reference, the phase matrix uses the local meridional      *
*  planes which are in general different from the scattering plane and *
*  also different for incident and scattered light.                    *
*  On entry, the layer number, the expansion coefficients, the mu and  *
*            phi values and the generalized spherical functions should *
*            be supplied.                                              *
*  On exit,  the phase matrix is returned through the arrays Zmin and  *
*            Zplus.                                                    *
*                                                                      *
*     Zmin   : phase matrix for extra mu points Z(-mu,mu',phi-phi')    *
*     Zplus  : phase matrix for extra mu points Z(+mu,mu',phi-phi')    *
*                                                                      *
*  Referring to the functions Plmn(u) in de Haan et al. (1987), we have*
*  for uplus = mu*mu0 + sqrt(1-mu**2) sqrt(1-mu0**2) cos(phi-phi0)     *
*  where mu, mu0 and phi-phi0 refer to to geometry number i :          *
*                                                                      *
*     gsfpl(l,1,i) = Pl00(uplus)                                       *
*     gsfpl(l,2,i) = Pl02(uplus)                                       *
*     gsfpl(l,3,i) = Pl22(uplus)                                       *
*     gsfpl(l,4,i) = Pl2-2(uplus)                                      *
*                                                                      *
*  and for umin = -mu*mu0 + sqrt(1-mu**2) sqrt(1-mu0**2) cos(phi-phi0) *
*  where mu, mu0 an phi-phi0 refer to to geometry number i :           *
*                                                                      *
*     gsfmi(l,1,i) = Pl00(umin)                                        *
*     gsfmi(l,2,i) = Pl02(umin)                                        *
*     gsfmi(l,3,i) = Pl22(umin)                                        *
*     gsfmi(l,4,i) = Pl2-2(umin)                                       *
*                                                                      *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter ( pi=3.1415926535897932384D0 )
      parameter ( twopi=2.D0*pi )
      parameter ( radfac=pi/180.D0 )
      dimension coefs(4,4,0:NDcoef,NDlay), ncoefs(NDlay)
     +        , imu(NDgeom), imu0(NDgeom), phi(NDgeom), xmu(NDmu)
     +        , Zmin(4,4,NDgeom), Zplus(4,4,NDgeom)
     +        , gsfmi(0:NDcoef,4,NDgeom), gsfpl(0:NDcoef,4,NDgeom)
*----------------------------------------------------------------------*
*  Local variables :                                                   *
*     Fmin   : scattering matrix for scattering angles > 90 degrees    *
*     Fplus  : scattering matrix for scattering angles < 90 degrees    *
*----------------------------------------------------------------------*
      dimension Fmin(4,4), Fplus(4,4)
      logical onemu, onemu0

      if (nmat .eq. 1) then
*----------------------------------------------------------------------*
*         Case without polarization                                    *
*----------------------------------------------------------------------*
          do 200 i=1, ngeom
              summi = 0.D0
              sumpl = 0.D0
              do 100 l=0, ncoefs(layer)
                  summi = summi + coefs(1,1,l,layer)*gsfmi(l,1,i)
                  sumpl = sumpl + coefs(1,1,l,layer)*gsfpl(l,1,i)
  100         continue
              Zmin(1,1,i)  = summi
              Zplus(1,1,i) = sumpl
  200     continue
      else
*----------------------------------------------------------------------*
*         Case with polarization                                       *
*----------------------------------------------------------------------*
          do 300 i=1, ngeom
            do 500 j=1, nmat
                do 600 k=1, nmat
                    Fmin(k,j)    = 0.D0
                    Fplus(k,j)   = 0.D0
                    Zmin(k,j,i)  = 0.D0
                    Zplus(k,j,i) = 0.D0
  600           continue
  500       continue
*----------------------------------------------------------------------*
*           Eqs. (68)-(73) with (69) combined with (70)                *
*----------------------------------------------------------------------*
            do 400 l=0, ncoefs(layer)
              Fmin(1,1) = Fmin(1,1) + coefs(1,1,l,layer)*gsfmi(l,1,i)
              Fplus(1,1)= Fplus(1,1)+ coefs(1,1,l,layer)*gsfpl(l,1,i)
              sum23 = coefs(2,2,l,layer) + coefs(3,3,l,layer)
              dif23 = coefs(2,2,l,layer) - coefs(3,3,l,layer)
              Fmin(2,2) = Fmin(2,2) + sum23*gsfmi(l,3,i)
              Fplus(2,2)= Fplus(2,2)+ sum23*gsfpl(l,3,i)
              Fmin(3,3) = Fmin(3,3) + dif23*gsfmi(l,4,i)
              Fplus(3,3)= Fplus(3,3)+ dif23*gsfpl(l,4,i)
              Fmin(1,2) = Fmin(1,2) + coefs(1,2,l,layer)*gsfmi(l,2,i)
              Fplus(1,2)= Fplus(1,2)+ coefs(1,2,l,layer)*gsfpl(l,2,i)
              if (nmat .eq. 4) then
                 Fmin(4,4) = Fmin(4,4) + coefs(4,4,l,layer)*gsfmi(l,1,i)
                 Fplus(4,4)= Fplus(4,4)+ coefs(4,4,l,layer)*gsfpl(l,1,i)
                 Fmin(3,4) = Fmin(3,4) + coefs(3,4,l,layer)*gsfmi(l,2,i)
                 Fplus(3,4)= Fplus(3,4)+ coefs(3,4,l,layer)*gsfpl(l,2,i)
              endif
  400       continue
*----------------------------------------------------------------------*
*           Solve F22 and F33 from combined results of Eqs. (69)-(70)  *
*----------------------------------------------------------------------*
            Fmin(2,2)  = 0.5D0*(Fmin(2,2)  + Fmin(3,3))
            Fplus(2,2) = 0.5D0*(Fplus(2,2) + Fplus(3,3))
            Fmin(3,3)  = Fmin(2,2)  - Fmin(3,3)
            Fplus(3,3) = Fplus(2,2) - Fplus(3,3)
*----------------------------------------------------------------------*
*           Use symmetry relations contained in Eq. (3)                *
*----------------------------------------------------------------------*
            Fmin(2,1)  = Fmin(1,2)
            Fplus(2,1) = Fplus(1,2)
            Fmin(4,3)  = -Fmin(3,4)
            Fplus(4,3) = -Fplus(3,4)
*----------------------------------------------------------------------*
*  Now the scattering matrix is contained in arrays Fmin and Fplus.    *
*----------------------------------------------------------------------*
*  Start the rotations of reference planes.                            *
*  We use Eqs. (6)-(10) from the thesis of Hovenier page I-2 and I-3.  *
*       Hovenier               'plus-case'         'min-case'          *
*     u, u', phi-phi'          mu, mu0, phi       -mu, -mu0, phi       *
*     +-cos(i1)                 cosi1p              cos1im             *
*     +-cos(i2)                 cosi2p              cosi2m             *
*     cos(theta)                costhp              costhm             *
*     cos(2*i1)                 cs2i1p              cs2i1m             *
*     sin(2*i1)                 s2i1p               s2i1m              *
*     cos(2*i2)                 cs2i2p              cs2i2m             *
*     sin(2*i2)                 s2i2p               s2i2m              *
*  If the azimuth change is a mutiple of pi, the reference planes      *
*  coincide with the scattering plane: no rotations needed.            *
*  Otherwise, special cases occur when mu=1 or mu0=1. When mu=mu0=1,   *
*  the rotations are undefined : we choose to do no rotations at all.  *
*  The signs in Hovenier's Eqs. (7)-(8) are treated separately later.  *
*----------------------------------------------------------------------*
           ph = phi(i)*radfac
           if (dabs(ph) .lt. 1.D-10) goto 700
           if (dabs(ph-pi) .lt. 1.D-10) goto 700
           if (dabs(ph-twopi) .lt. 1.D-10) goto 700
           emu0 = xmu(imu0(i))
           emu  = xmu(imu(i))
           onemu0 = .false.
           onemu  = .false.
           if (dabs(emu0-1.D0) .lt. 1.D-6) onemu0 = .true.
           if (dabs(emu -1.D0) .lt. 1.D-6) onemu  = .true.
           if (onemu .and. onemu0) goto 700
           cosph = dcos(ph)
           if (onemu0) then
               cosi1m =  cosph
               cosi1p =  cosph
               cosi2m = -1.D0
               cosi2p = -1.D0
           else if (onemu) then
               cosi1m =  1.D0
               cosi1p = -1.D0
               cosi2m = -cosph
               cosi2p =  cosph
           else
*----------------------------------------------------------------------*
*  Substitute Hovenier's Eq. (6) into Eqs. (7)-(8) and divide out the  *
*  factors sqrt(1-mu**2) and sqrt(1-mu0**2).                           *
*  Realize that cos(theta) is equal to the first order Legendre        *
*  polynomial contained in gsfmi(1,1,i) and gsfpl(1,1,i).              *
*----------------------------------------------------------------------*
               rtmu0  = dsqrt(1.D0-emu0*emu0)
               rtmu   = dsqrt(1.D0-emu*emu)
               rmu0   = rtmu0*emu
               rmu    = rtmu*emu0
               costhm = gsfmi(1,1,i)
               costhp = gsfpl(1,1,i)
               thfacm = 1.D0/dsqrt(1.D0-costhm*costhm)
               thfacp = 1.D0/dsqrt(1.D0-costhp*costhp)
               cosi1m = ( rmu0 +  rmu*cosph)*thfacm
               cosi1p = (-rmu0 +  rmu*cosph)*thfacp
               cosi2m = (-rmu  - rmu0*cosph)*thfacm
               cosi2p = (-rmu  + rmu0*cosph)*thfacp
          endif
*----------------------------------------------------------------------*
*  The signs in Hovenier's Eqs. (7)-(8) only affect his Eq. (9) and    *
*  are combined with the factor 2 in Eq. (9) into variable sgn2.       *
*  Use shorthand sin2 for 1-(cos(beta))**2 in Eqs. (9)-(10).           *
*----------------------------------------------------------------------*
            if (ph .ge. pi) then
                sgn2 = 2.D0
            else
                sgn2 = -2.D0
            endif
            sin2  = 1.D0-cosi1m*cosi1m
            s2i1m = sgn2*dsqrt(sin2)*cosi1m
            c2i1m = 1.D0-sin2-sin2

            sin2  = 1.D0-cosi1p*cosi1p
            s2i1p = sgn2*dsqrt(sin2)*cosi1p
            c2i1p = 1.D0-sin2-sin2

            sin2  = 1.D0-cosi2m*cosi2m
            s2i2m = sgn2*dsqrt(sin2)*cosi2m
            c2i2m = 1.D0-sin2-sin2

            sin2  = 1.D0-cosi2p*cosi2p
            s2i2p = sgn2*dsqrt(sin2)*cosi2p
            c2i2p = 1.D0-sin2-sin2
*----------------------------------------------------------------------*
*  Perform the rotations according to Hovenier Eqs. (1) and (2).       *
*----------------------------------------------------------------------*
            Zmin(1,1,i)  = Fmin(1,1)
            Zmin(2,1,i)  = c2i2m*Fmin(2,1)
            Zmin(3,1,i)  = s2i2m*Fmin(2,1)

            Zmin(1,2,i)  = c2i1m*Fmin(2,1)
            Zmin(2,2,i)  = c2i1m*c2i2m*Fmin(2,2)-s2i1m*s2i2m*Fmin(3,3)
            Zmin(3,2,i)  = c2i1m*s2i2m*Fmin(2,2)+s2i1m*c2i2m*Fmin(3,3)

            Zmin(1,3,i)  =-s2i1m*Fmin(2,1)
            Zmin(2,3,i)  =-s2i1m*c2i2m*Fmin(2,2)-c2i1m*s2i2m*Fmin(3,3)
            Zmin(3,3,i)  =-s2i1m*s2i2m*Fmin(2,2)+c2i1m*c2i2m*Fmin(3,3)

            Zplus(1,1,i) = Fplus(1,1)
            Zplus(2,1,i) = c2i2p*Fplus(2,1)
            Zplus(3,1,i) = s2i2p*Fplus(2,1)

            Zplus(1,2,i) = c2i1p*Fplus(2,1)
            Zplus(2,2,i) = c2i1p*c2i2p*Fplus(2,2)-s2i1p*s2i2p*Fplus(3,3)
            Zplus(3,2,i) = c2i1p*s2i2p*Fplus(2,2)+s2i1p*c2i2p*Fplus(3,3)

            Zplus(1,3,i) =-s2i1p*Fplus(2,1)
            Zplus(2,3,i) =-s2i1p*c2i2p*Fplus(2,2)-c2i1p*s2i2p*Fplus(3,3)
            Zplus(3,3,i) =-s2i1p*s2i2p*Fplus(2,2)+c2i1p*c2i2p*Fplus(3,3)

            if (nmat .eq. 4) then
                Zmin(4,2,i)  = -s2i1m*Fmin(3,4)
                Zmin(4,3,i)  = -c2i1m*Fmin(3,4)

                Zmin(2,4,i)  = -s2i2m*Fmin(3,4)
                Zmin(3,4,i)  = c2i2m*Fmin(3,4)
                Zmin(4,4,i)  = Fmin(4,4)

                Zplus(4,2,i)  = -s2i1p*Fplus(3,4)
                Zplus(4,3,i)  = -c2i1p*Fplus(3,4)

                Zplus(2,4,i)  = -s2i2p*Fplus(3,4)
                Zplus(3,4,i)  = c2i2p*Fplus(3,4)
                Zplus(4,4,i)  = Fplus(4,4)
            endif
*----------------------------------------------------------------------*
*  End of rotations, take next geometry (if there is one).             *
*----------------------------------------------------------------------*
            goto 300
*----------------------------------------------------------------------*
*  Skip the rotations, copy scattering matrix F into phase matrix Z.   *
*----------------------------------------------------------------------*
  700       do 710 k=1, nmat
                do 720 j=1, nmat
                    Zmin(k,j,i)  = Fmin(k,j)
                    Zplus(k,j,i) = Fplus(k,j)
  720           continue
  710       continue
  300     continue
      endif
      return
      end
      subroutine fillup( NDsup, Rm ,Tm, nmat, nmutot )
*----------------------------------------------------------------------*
*  Fill the upper triangle of a supermatrix using symmetry relations   *
*  for a vertically homogeneous layer.                                 *
*  By the upper triangle we mean mu < mu0 : nmat by nmat submatrices   *
*  are not split.                                                      *
*  See de Haan et al. (1987): Astron. Astrophys. 183, p. 371           *
*  Eqs. (96)-(97)                                                      *
*                 R = q3 R~ q3                                         *
*                 T = q4 T~ q4                                         *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension Rm(NDsup,NDsup), Tm(NDsup,NDsup)
     +        , q3(4), q4(4)
      do 10 k=1, 4
          q3(k) = 1.D0
          q4(k) = 1.D0
   10 continue
      q3(3) = -1.D0
      q4(4) = -1.D0
      do 400 mu0=1, nmutot
          jbase = (mu0-1)*nmat
          do 300 mu=1, mu0-1
              ibase = (mu-1)*nmat
              do 200 ki=1, nmat
                  i = ibase+ki
                  do 100 kj=1, nmat
                      j = jbase+kj
                      Rm(i,j) = q3(ki)*Rm(j,i)*q3(kj)
                      Tm(i,j) = q4(ki)*Tm(j,i)*q4(kj)
  100             continue
  200         continue
  300     continue
  400 continue
      return
      end
      subroutine first(NDmu, NDlay, NDgeom, NDcoef
     +                , coefs, ncoefs, xmu, imu0, imu
     +                , phi, ngeom, nmat
     +                , a, b, nlayer
     +                , R1, T1, R1st, T1st, e, e0 )
*----------------------------------------------------------------------*
*  Calculate the contribution of first order scattering to reflection  *
*  and transmission for all specified geometries. The formulae can be  *
*  found in de Haan et al. (1987) Eqs. (16)-(18), (5),(3) and (68)-(73)*
*  On entry:                                                           *
*     coefs   : expansion coefficients of the scattering matrix in     *
*               generalized spherical functions                        *
*     ncoefs  : order of the highest nonzero term in the expansion     *
*     xmu     : all different mu values (integration and extra points) *
*     imu0    : array indicating, for each geometry, the               *
*               position in array xmu where the mu0 value can be found *
*     phi     : azimuth of the emerging direction relative to the      *
*               direction of incidence in degrees                      *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     a       : single scattering albedo for each layer                *
*     b       : optical thickness for each layer                       *
*     nlayer  : number of layers in the atmosphere                     *
*  On exit:                                                            *
*     R1      : first order reflection matrix of the atmosphere        *
*     T1      : first order transmission matrix of the atmosphere      *
*     R1st    : same as R1 but for illumination from below             *
*     T1st    : same as T1 but for illumination from below             *
*     e       : factor exp(-b/mu) for each layer                       *
*     e0      : factor exp(-b/mu0) for each layer                      *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension coefs(4,4,0:NDcoef,NDlay), ncoefs(NDlay)
     +        , xmu(NDmu), imu(NDgeom), imu0(NDgeom), phi(NDgeom)
     +        , a(NDlay), b(NDlay)
      dimension e(NDlay,NDgeom), e0(NDlay,NDgeom)
     +        , R1(4,4,NDgeom), R1st(4,4,NDgeom)
     +        , T1(4,4,NDgeom), T1st(4,4,NDgeom)
     +        , Tdir(NDgeom),   Tdir0(NDgeom)
*----------------------------------------------------------------------*
*  Phase matrix data :                                                 *
*     Zmin   : phase matrix for extra points Z(-mu,mu',phi-phi')       *
*     Zplus  : phase matrix for extra points Z(+mu,mu',phi-phi')       *
*----------------------------------------------------------------------*
      dimension Zmin(4,4,NDgeom), Zplus(4,4,NDgeom)
     +        , gsfmi(0:NDcoef,4,NDgeom), gsfpl(0:NDcoef,4,NDgeom)
     +        , sgn(4,4)
*----------------------------------------------------------------------*
*  Precompute e-powers dexp(-b/mu) abs dexp(-b/mu0), and the           *
*  Generalized Spherical Functions                                     *
*----------------------------------------------------------------------*
      call setexp(NDmu, NDlay, NDgeom
     +           , xmu, imu, imu0, ngeom, b, nlayer, e, e0 )
      lmax = 0
      do 10 layer=1, nlayer
          if (ncoefs(layer) .gt. lmax) lmax = ncoefs(layer)
   10 continue
      call setgsf(NDmu, NDgeom, NDcoef, lmax, xmu, imu0, imu, phi, ngeom
     +            , nmat, gsfpl, gsfmi)
*----------------------------------------------------------------------*
*  Initialize reflection and transmission matrices for all geometries. *
*  The factors Tdir and Tdir0 represent the direct transmission        *
*  through the layers already treated (none at the start) for          *
*  emerging and incident directions respectively.                      *
*----------------------------------------------------------------------*
      do 300 i=1, ngeom
          do 200 k=1, nmat
              do 100 l=1, nmat
                  R1(k,l,i)  = 0.D0
                  R1st(k,l,i)= 0.D0
                  T1(k,l,i)  = 0.D0
                  T1st(k,l,i)= 0.D0
  100         continue
  200     continue
          Tdir(i)  = 1.D0
          Tdir0(i) = 1.D0
  300 continue
*----------------------------------------------------------------------*
*  Fill array sgn with the signs needed to obtain R* and T* from R and *
*  T. Eqs. (98)-(99) or Hovenier's thesis p. I-6 Eqs. (q)-(r) :        *
*                       1   1  -1  -1                                  *
*             sgn  =    1   1  -1  -1                                  *
*                      -1  -1   1   1                                  *
*                      -1  -1   1   1                                  *
*----------------------------------------------------------------------*
      do 320 k=1, 4
          do 310 l=1, 4
              sgn(k,l) = 1.D0
              if (k .gt. 2) sgn(k,l) = -sgn(k,l)
              if (l .gt. 2) sgn(k,l) = -sgn(k,l)
  310     continue
  320 continue
*----------------------------------------------------------------------*
*  Start loop over the layers. Fasemx calculates the phase matrix for  *
*  all geometries. Then start a loop over the geometries to update     *
*  reflection and transmission.                                        *
*----------------------------------------------------------------------*
      do 700 j=1, nlayer
          call fasemx( NDmu, NDlay, NDgeom, NDcoef
     +               , j, coefs, ncoefs, xmu, imu, imu0
     +                               , phi, ngeom
     +                               , gsfpl, gsfmi, nmat, Zmin, Zplus )
          do 600 i=1, ngeom
              emu  = xmu(imu(i))
              emu0 = xmu(imu0(i))
*----------------------------------------------------------------------*
*             Special case mu = mu0 = 0 for reflection                 *
*----------------------------------------------------------------------*
              if (dabs(emu+emu0) .lt. 1.D-10) then
                  hr   = 0.D0
                  hrTT = 0.D0
              else
                  hr  = 0.25D0*a(j)*(1.D0-e(j,i)*e0(j,i))/(emu+emu0)
                  hrTT= hr*Tdir(i)*Tdir0(i)
              endif
*----------------------------------------------------------------------*
*             Special case mu = mu0 for transmission                   *
*----------------------------------------------------------------------*
              if (dabs(emu-emu0) .gt. 1.D-10) then
                  ht   = 0.25D0*a(j)*(e(j,i)-e0(j,i))/(emu-emu0)
              else
                  ht   = 0.25D0*a(j)*b(j)*e(j,i)/(emu*emu)
              endif
              htT  = ht*Tdir(i)
              htT0 = ht*Tdir0(i)
*----------------------------------------------------------------------*
*             Update reflection and transmission                       *
*----------------------------------------------------------------------*
              eji  = e(j,i)
              e0ji = e0(j,i)
              do 500 k=1, nmat
                 do 400 l=1, nmat
                    R1(k,l,i)   = R1(k,l,i)*eji*e0ji + Zmin(k,l,i)*hr
                    R1st(k,l,i) = R1st(k,l,i)
     +                                      + sgn(k,l)*Zmin(k,l,i)*hrTT
                    T1(k,l,i)   = T1(k,l,i)*e0ji     + Zplus(k,l,i)*htT
                    T1st(k,l,i) = T1st(k,l,i)*eji
     +                                      + sgn(k,l)*Zplus(k,l,i)*htT0
  400            continue
  500         continue
*----------------------------------------------------------------------*
*             Update direct transmission of already treated layers     *
*----------------------------------------------------------------------*
              Tdir(i)  = Tdir(i)*e(j,i)
              Tdir0(i) = Tdir0(i)*e0(j,i)
  600     continue
  700 continue
*----------------------------------------------------------------------*
*  End of loop over layers                                             *
*----------------------------------------------------------------------*
      return
      end
      subroutine firsti( NDmu, NDlay, NDgeom
     +                 , xmu, imu, imu0, ngeom, nmat, nlayer
     +                 , e, e0, Rf
     +                 , R1, T1, R1st, T1st, R1i, D1i )
*----------------------------------------------------------------------*
*  Calculate first order reflection matrix and first order matrix      *
*  describing the downwards radiation just above the interface.        *
*  Equation numbers in this subroutine refer to Chowdhary (1991).      *
*  On entry:                                                           *
*     xmu     : all different mu values (integration and extra points) *
*     imu     : index where mu can be found in array xmu for each      *
*               geometry                                               *
*     imu0    : index where mu0 can be found in array xmu for each     *
*               geometry                                               *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     nlayer  : number of layers in the atmosphere                     *
*     e       : exp(-b/mu) for each layer and for each geometry        *
*     e0      : exp(-b/mu0) for each layer and for each geometry       *
*     Rf      : Fresnel matrix given by Eq. (22) coded as follows      *
*                  1,1 element of RF(mu) in Rf(mu,1)                   *
*                  1,2 element of RF(mu) in Rf(mu,2)                   *
*                  3,3 element of RF(mu) in Rf(mu,3)                   *
*               for all possible mu values                             *
*     R1      : first order reflection matrix of the atmosphere        *
*     T1      : first order transmission matrix of the atmosphere      *
*     R1st    : same as R1 but for illumination from below             *
*     T1st    : same as T1 but for illumination from below             *
*  On exit:                                                            *
*     R1i     : first order reflection matrix of atmosphere-interface  *
*               system                                                 *
*     D1i     : first order downward radiation field just above the    *
*               interface                                              *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), imu(NDgeom), imu0(NDgeom)
     +        , e(NDlay,NDgeom), e0(NDlay,NDgeom), Rf(NDmu,3)
     +        , R1(4,4,NDgeom),  T1(4,4,NDgeom)
     +        , R1st(4,4,NDgeom), T1st(4,4,NDgeom)
     +        , R1i(4,4,NDgeom), D1i(4,4,NDgeom)
      dimension a(NDgeom), a0(NDgeom)
     +        , b(NDgeom), b0(NDgeom)
     +        , c(NDgeom), c0(NDgeom), X(4,4)
      do 200 i=1, ngeom
*----------------------------------------------------------------------*
*  Calculate ebmu = exp(-b/mu) and ebmu0 = exp(-b/mu0)                 *
*----------------------------------------------------------------------*
          ebmu  = 1.d0
          ebmu0 = 1.d0 
          do 100 n=1, nlayer
              ebmu  = ebmu*e(n,i)
              ebmu0 = ebmu0*e0(n,i)
  100     continue
*----------------------------------------------------------------------*
*  We use Eqs. (79)-(82) but with the Fourier series summed according  *
*  to Eq. (37). Using Eqs. (38) and (61) we obtain from Eq. (79):      *
*                                                                      *
*  D1(mu,mu0,phi-phi0) =                                               *
*  = T1'(mu,mu0,phi-phi0) + R1'*(mu,mu0,phi-phi0) 2 mu0 Rf(mu0) ebmu0  *
*                                                                      *
*                                                 ( a0 b0 0  0 )       *
*  = T1'(mu,mu0,phi-phi0) + R1'*(mu,mu0,phi-phi0) ( b0 a0 0  0 )       *
*                                                 ( 0  0  c0 0 )       *
*                                                 ( 0  0  0 c0 )       *
*                                                                      *
*  There is an error in Eq. (82): T1' must be T1'* .                   *
*  Using Eqs. (37), (38), (61), (80) and (82) we find:                 *
*                                                                      *
*  R1(mu,mu0,phi-phi0) = R1'(mu,mu0,phi-phi0) +                        *
*                                                                      *
*          ( a  b  0  0 )                                              *
*        + ( b  a  0  0 ) D1(mu,mu0,phi-phi0) + X(mu,mu0,phi-phi0)     *
*          ( 0  0  c  0 )                                              *
*          ( 0  0  0  c )                                              *
*                                                                      *
*  where a, b and c are elements of   2 mu Rf(mu) ebmu , and           *
*                                                                      *
*  X(mu,mu0,phi-phi0) = T1'*(mu,mu0,phi-phi0) 2 mu0 Rf(mu0) ebmu0      *
*----------------------------------------------------------------------*
          a0(i) = 2.d0*xmu(imu0(i))*ebmu0*Rf(imu0(i),1)
          b0(i) = 2.d0*xmu(imu0(i))*ebmu0*Rf(imu0(i),2)
          c0(i) = 2.d0*xmu(imu0(i))*ebmu0*Rf(imu0(i),3)
          a(i)  = 2.d0*xmu(imu(i)) *ebmu *Rf(imu(i), 1)
          b(i)  = 2.d0*xmu(imu(i)) *ebmu *Rf(imu(i), 2)
          c(i)  = 2.d0*xmu(imu(i)) *ebmu *Rf(imu(i), 3)
  200 continue
*----------------------------------------------------------------------*
*  Treat case without polarization                                     *
*----------------------------------------------------------------------*
      if (nmat.eq.1) then
          do 300 i=1, ngeom
            D1i(1,1,i) = R1st(1,1,i)*a0(i) + T1(1,1,i)
            R1i(1,1,i) = R1(1,1,i) + a(i)*D1i(1,1,i) + T1st(1,1,i)*a0(i)
  300     continue
*----------------------------------------------------------------------*
*  Treat case of (3 X 3) approximation                                 *
*----------------------------------------------------------------------*
      else if (nmat.eq.3) then
        do 800 i=1, ngeom
          do 400 k=1, nmat
            D1i(k,1,i) = R1st(k,1,i)*a0(i) + R1st(k,2,i)*b0(i)
            D1i(k,2,i) = R1st(k,1,i)*b0(i) + R1st(k,2,i)*a0(i)
            D1i(k,3,i) = R1st(k,3,i)*c0(i)
            X(k,1)     = T1st(k,1,i)*a0(i) + T1st(k,2,i)*b0(i)
            X(k,2)     = T1st(k,1,i)*b0(i) + T1st(k,2,i)*a0(i)
            X(k,3)     = T1st(k,3,i)*c0(i)
  400     continue
          do 600 l=1, nmat
            do 500 k=1, nmat
              D1i(k,l,i) = D1i(k,l,i) + T1(k,l,i)
              R1i(k,l,i) = R1(k,l,i) + X(k,l)
  500       continue
  600     continue
          do 700 l=1, nmat
            R1i(1,l,i) = R1i(1,l,i) + a(i)*D1i(1,l,i) + b(i)*D1i(2,l,i) 
            R1i(2,l,i) = R1i(2,l,i) + a(i)*D1i(2,l,i) + b(i)*D1i(1,l,i) 
            R1i(3,l,i) = R1i(3,l,i) + c(i)*D1i(3,l,i)
  700     continue
  800   continue 
*----------------------------------------------------------------------*
*  Treat case with polarization fully included                         *
*----------------------------------------------------------------------*
      else if (nmat.eq.4) then
        do 1050 i=1, ngeom
          do 850 k=1, nmat
            D1i(k,1,i) = R1st(k,1,i)*a0(i) + R1st(k,2,i)*b0(i)
            D1i(k,2,i) = R1st(k,1,i)*b0(i) + R1st(k,2,i)*a0(i)
            D1i(k,3,i) = R1st(k,3,i)*c0(i)
            D1i(k,4,i) = R1st(k,4,i)*c0(i)
            X(k,1)    = T1st(k,1,i)*a0(i) + T1st(k,2,i)*b0(i)
            X(k,2)    = T1st(k,1,i)*b0(i) + T1st(k,2,i)*a0(i)
            X(k,3)    = T1st(k,3,i)*c0(i)
            X(k,4)    = T1st(k,4,i)*c0(i)
  850     continue
          do 950 l=1, nmat
            do 900 k=1, nmat
              D1i(k,l,i) = D1i(k,l,i) + T1(k,l,i)
              R1i(k,l,i) = R1(k,l,i) + X(k,l)
  900       continue
  950     continue
          do 1000 l=1, nmat
            R1i(1,l,i) = R1i(1,l,i) + a(i)*D1i(1,l,i) + b(i)*D1i(2,l,i) 
            R1i(2,l,i) = R1i(2,l,i) + a(i)*D1i(2,l,i) + b(i)*D1i(1,l,i) 
            R1i(3,l,i) = R1i(3,l,i) + c(i)*D1i(3,l,i)
            R1i(4,l,i) = R1i(4,l,i) + c(i)*D1i(4,l,i)
 1000     continue
 1050   continue 
      else
          write(*,*) ' firsti: illegal value for nmat = ',nmat
          stop 'in firsti illegal value for nmat'
      endif
      return
      end 
      subroutine firstm( NDmu, NDsup, NDgeom
     +                 , m, xmu, imu, imu0, smf, nmutot, ngeom, nmat
     +                 , Zmplus, Zmmin, a, b, ebmu
     +                 , xRm1, xRm1s, xTm1, xTm1s )
*----------------------------------------------------------------------*
*  Calculate the first order scattering contribution to the m-th       *
*  Fourier component of reflection and transmission of a homogeneous   *
*  layer. Formulae can be found in de Haan et al. (1987)               *
*  Eqs. (A-21)-(A-23)                                                  *
*  The calculation is performed for all geometries.
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), imu(NDgeom), imu0(NDgeom), smf(NDmu)
     +        , Zmmin(NDsup,NDsup), Zmplus(NDsup,NDsup), ebmu(NDmu)
     +        , xRm1(4,4,NDgeom), xTm1(4,4,NDgeom)
     +        , xRm1s(4,4,NDgeom), xTm1s(4,4,NDgeom)
      dimension sgn(4,4)
      quarta = 0.25D0*a
*----------------------------------------------------------------------*
*  Fill array sgn with the signs needed to obtain R* and T* from R and *
*  T. Eqs. (98)-(99) or Hovenier's thesis p. I-6 Eqs. (q)-(r) :        *
*                       1   1  -1  -1                                  *
*             sgn  =    1   1  -1  -1                                  *
*                      -1  -1   1   1                                  *
*                      -1  -1   1   1                                  *
*----------------------------------------------------------------------*
      do 200 l=1, nmat
          do 100 k=1, nmat
              sgn(k,l) = 1.0d0
              if (l .ge. 3) sgn(k,l) = -sgn(k,l)
              if (k .ge. 3) sgn(k,l) = -sgn(k,l)
  100     continue
  200 continue
*----------------------------------------------------------------------*
*  Start loop over all geometries                                      *
*----------------------------------------------------------------------*
      do 500 igeom=1, ngeom
          emu0 = xmu(imu0(igeom))
          emu  = xmu(imu(igeom))
          e0   = ebmu(imu0(igeom))
          e    = ebmu(imu(igeom))
          e0e  = 1.D0 - e*e0
          im   = (imu(igeom)-1)*nmat
          jm   = (imu0(igeom)-1)*nmat
          if (e0e .gt. 1.D-3) then
              h = emu0+emu
              if (h .gt. 1.D-10) h = 1.D0/h
              hR = quarta*h*e0e
              if (dabs(emu0-emu) .gt. 1.D-10) then
                  hT = quarta*(e-e0)/(emu-emu0)
              else
                  h = 0.D0
                  if (emu .gt. 1.D-10) h = b*e /(emu*emu)
                  hT = quarta*h
              endif
          else
*----------------------------------------------------------------------*
*  use Taylor series to avoid loss of accuarcy when b << 1             *
*----------------------------------------------------------------------*
              bmu0   = b/emu0
              bmu    = b/emu
              bmumu0 = bmu0+bmu
              h      = 1.D0-0.5D0*bmumu0*(1.D0-bmumu0/3.D0)
              y      = quarta*bmu/emu0
              hR     = y*h
              hT     = y*(h-bmu*bmu0/6.D0)
          endif
          do 400 k=1, nmat
              ik = im+k
              do 300 l=1, nmat
                  jl = jm+l
                  xRm1(k,l,igeom) = hR*Zmmin(ik,jl)
                  xTm1(k,l,igeom) = hT*Zmplus(ik,jl)
                  xRm1s(k,l,igeom)= sgn(l,k)*xRm1(k,l,igeom)
                  xTm1s(k,l,igeom)= sgn(l,k)*xTm1(k,l,igeom)
  300         continue
  400     continue
  500 continue
      return
      end
      subroutine ad_fluxes( NDmu, NDsup, NDgeom, xmu, smf, imu, imu0
     +                    , ngeom, nmug, nmutot, nmat, ebmu
     +                    , Rm, Rmst, Tm, Rm1, Rm1st, Tm1
     +                    , Rflux, Tflux, URU, UTU )
*----------------------------------------------------------------------*
*  Calculate reflected and transmitted fluxes of the atmosphere.       *
*  On entry:                                                           *
*     xmu     : all different mu values (integration and extra points) *
*     smf     : supermatrix factors dsqrt(2*w*mu), or 1 for extra pts. *
*     imu     : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu value can be found  *
*     imu0    : array filled by setmu with, for each geometry, the     *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*     nmug    : number of Gauss points for the mu integrations         *
*     nmutot  : total number of mu points (Gauss and extra points)     *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     ebmu    : array with factor dexp(-b/mu) for each mu value in xmu *
*     Rm      : supermatrix for reflection of the atmosphere           *
*     Rmst    : idem, but for illumination from below                  *
*     Tm      : supermatrix for diffuse transmission of the atmosphere *
*     Rm1     : first order scattering contribution to Rm              *
*     Rm1st   : first order scattering contribution to Rmst            *
*     Tm1     : first order scattering contribution to Tm              *
*  On exit:                                                            *
*     Rflux   : reflected flux (see below)                             *
*     Tflux   : transmitted flux (see below)                           *
*     URU     : spherical reflected flux (see below)                   *
*     UTU     : spherical transmitted flux (see below)                 *
*                                                                      *
*  Calculate reflected and transmitted fluxes by integrating the       *
*  reflection and transmission matrix. Hence fluxes are (nmat X nmat)  *
*  matrices. At least half of the elements is zero, but we do not want *
*  to use this fact to devise a more compact (and complicated) way to  *
*  store fluxes; the memory penalty is small.                          *
*                                                                      *
*  The m=0 Fourier component of the reflection matrix is Rm(mu,mu0).   *
*  We define the following fluxes as integrals over the incident       *
*  direction or the outgoing direction or both:                        *
*                 1                                                    *
*     RU(mu)  = 2 int dmu0 mu0 Rm(mu,mu0)                              *
*                 0                                                    *
*                                                                      *
*                 1                                                    *
*     UR(mu0) = 2 int dmu mu Rm(mu,mu0)                                *
*                 0                                                    *
*                                                                      *
*                 1            1                                       *
*     URU    = 2 int dmu mu 2 int dmu0 mu0 Rm(mu,mu0)                  *
*                 0            0                                       *
*                                                                      *
*  For the transmission we have similar definitions but with R         *
*  replaced by T. All definitions have an obvious counterpart          *
*  for illumination from below.                                        *
*  We also calculate for each flux the zero-th and first order         *
*  scattering contribution. For the reflection the zero-th order is    *
*  of course zero.                                                     *
*                                                                      *
*  The fluxes that have an argument mu or mu0 are calculated for all   *
*  geometries. We use the following datastructure to store them:       *
*                                                                      *
*      Rflux( k,l, igeom, muarg, illum, iorder ) and                   *
*      Tflux( k,l, igeom, muarg, illum, iorder ).                      *
*                                                                      *
*  Here Rflux is a reflected flux and Tflux a transmitted flux. The    *
*  meaning of the indices is as follows:                               *
*                                                                      *
*     index    values           meaning                                *
*   -----------------------------------------------------------------  *
*     k,l      1,...,nmat       indices in a (nmat X nmat) matrix      *
*                                                                      *
*     igeom    1,...,ngeom      index of the geometry considered       *
*                                                                      *
*     muarg    Nmu0             integral over outgoing direction, so   *
*                               we have UR(mu0) and UT(mu0)            *
*              Nmu              integral over incoming direction, so   *
*                               we have RU(mu) and TU(mu)              *
*                                                                      *
*    illum     Nabove           illumination from above                *
*              Nbelow           illumination from below                *
*                                                                      *
*    iorder    Nzero            only zero order scattering             *
*              Nfirst           only first order scattering            *
*              Nall             all orders of scattering (including 0) *
*   -----------------------------------------------------------------  *
*                                                                      *
*  The values Nmu0, Nmu, Nabove, Nbelow, Nzero, Nfirst and Nall  are   *
*  coded in a parameter statement.                                     *
*  For the spherical fluxes URU and UTU we use the following           *
*  datastructure:                                                      *
*                                                                      *
*      URU( k,l, illum, iorder ) and                                   *
*      UTU( k,l, illum, iorder ).                                      *
*                                                                      *
*  The meaning of the indices is the same as before.                   *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter( Nmu0=1,   Nmu=2
     +         , Nzero=1,  Nfirst=2 , Nall=3
     +         , Nabove=1, Nbelow=2 )
      dimension xmu(NDmu), ebmu(NDmu), smf(NDmu)
     +        , imu(NDgeom), imu0(NDgeom)
     +        , Rm(NDsup,NDsup),  Rmst(NDsup,NDsup),  Tm(NDsup,NDsup)
     +        , Rm1(NDsup,NDsup), Rm1st(NDsup,NDsup), Tm1(NDsup,NDsup)
     +        , Tmst(NDsup,NDsup)
      dimension Rflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , Tflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , URU(  4,4,              Nbelow, Nall)
     +        , UTU(  4,4,              Nbelow, Nall)
*----------------------------------------------------------------------*
* Initialize everything to zero                                        *
*----------------------------------------------------------------------*
      do 600 iorder=Nzero, Nall
          do 500 illum=Nabove, Nbelow
              do 400 l=1, 4
                  do 300 k=1, 4
                      do 200 igeom=1, NDgeom
                          do 100 muarg=Nmu0, Nmu
                             Rflux(k,l,igeom,muarg,illum,iorder) = 0.0d0
                             Tflux(k,l,igeom,muarg,illum,iorder) = 0.0d0
  100                     continue
  200                 continue
                      URU(k,l,illum,iorder) = 0.0d0
                      UTU(k,l,illum,iorder) = 0.0d0
  300             continue
  400         continue
  500     continue
  600 continue
*----------------------------------------------------------------------*
*  Reflection for illumination from above                              *
*----------------------------------------------------------------------*
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Rm, Nabove, Nall,  Rflux, URU )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Rm1, Nabove, Nfirst,  Rflux, URU )
*----------------------------------------------------------------------*
*  Transmission for illumination from above                            *
*----------------------------------------------------------------------*
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Tm, Nabove, Nall,  Tflux, UTU )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Tm1, Nabove, Nfirst,  Tflux, UTU )
*----------------------------------------------------------------------*
*  Reflection for illumination from below                              *
*----------------------------------------------------------------------*
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Rmst, Nbelow, Nall,  Rflux, URU )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Rm1st, Nbelow, Nfirst,  Rflux, URU )
*----------------------------------------------------------------------*
*  Transmission for illumination from below                            *
*----------------------------------------------------------------------*
      call Tstar(NDsup,  Tmst, Tm, nmat, nmutot ) 
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Tmst, Nbelow, Nall,  Tflux, UTU )
      call Tstar(NDsup,  Tmst, Tm1, nmat, nmutot ) 
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Tmst, Nbelow, Nfirst,  Tflux, UTU )
*----------------------------------------------------------------------*
*  Calculate the zero order (i.e. direct) transmitted fluxes           *
*----------------------------------------------------------------------*
      do 900 illum=Nabove, Nbelow
          do 800 igeom=1, ngeom
              mu0 = imu0(igeom)
              mu  = imu(igeom)
              do 700 k=1, nmat
                  Tflux(k,k,igeom,Nmu0,illum,Nzero) = ebmu(mu0)
                  Tflux(k,k,igeom,Nmu ,illum,Nzero) = ebmu(mu )
  700         continue
  800     continue
  900 continue
      do 2100 mu=1, nmug
          w = smf(mu)**2
          do 1700 k=1, nmat
              UTU(k,k,Nabove,Nzero) = UTU(k,k,Nabove,Nzero) + w*ebmu(mu)
 1700     continue
 2100 continue
*----------------------------------------------------------------------*
*  Use symmetry to obtain UTU for illumination from below              *
*----------------------------------------------------------------------*
      do 2400 iorder=Nzero, Nall
          do 2300 l=1, nmat
              do 2200 k=1, nmat
                  UTU(k,l,Nbelow,iorder) = UTU(l,k,Nabove,iorder) 
 2200         continue
 2300     continue
 2400 continue
*----------------------------------------------------------------------*
*  Add the zero order contributions to the higher order contributions  *
*  already obtained to get the sum over all orders.                    *
*----------------------------------------------------------------------*
      do 3100 illum=Nabove, Nbelow
          do 2800 muarg=Nmu0, Nmu
              do 2700 igeom=1, ngeom
                  do 2600 l=1, nmat
                      do 2500 k=1, nmat
                          Tflux(k,l,igeom,muarg,illum,Nall) =
     +                               Tflux(k,l,igeom,muarg,illum,Nall )
     +                             + Tflux(k,l,igeom,muarg,illum,Nzero)
 2500                 continue
 2600             continue
 2700         continue
 2800     continue
          do 3000 l=1, nmat
              do 2900 k=1, nmat
                  UTU(k,l,illum,Nall) = UTU(k,l,illum,Nall) 
     +                                + UTU(k,l,illum,Nzero)
 2900         continue
 3000     continue
 3100 continue
      return
      end
      subroutine fluxi( NDmu,NDsup, NDgeom, xmu, smf, imu, imu0
     +                   , ngeom, nmug, nmat, nmutot, ebmu
     +                   , xmut, Tf, Tfs, Rf, Rfs, xm
     +                   , Rm,  Dm,  Vm,  Wm 
     +                   , Rm1, Dm1, Vm1, Wm1 
     +                   , Riflux, Diflux, Tiflux, URUi, UDUi, UTUi )
*----------------------------------------------------------------------*
*  Calculate reflected and transmitted fluxes of the atmosphere-       *
*  water-air interface system, and downward fluxes just above the      *
*  interface.                                                          *
*  Equation numbers refer to Chowdhary (1991).                         *
*  On entry:                                                           *
*     xmu     : all different mu values (integration and extra points) *
*     smf     : supermatrix factors dsqrt(2*w*mu), or 1 for extra pts. *
*     imu     : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu value can be found  *
*     imu0    : array filled by setmu with, for each geometry, the     *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*     nmug    : number of Gauss points for the mu integrations         *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     ebmu    : array with factor dexp(-b/mu) for each mu value in xmu *
*     xmut   : factor defined in Eq. (7) for every mu value            *
*     Tf     : Fresnel matrix given by Eq. (23), coded as follows      *
*                  1,1 element of TF(mu) in Tf(mu,1)                   *
*                  1,2 element of TF(mu) in Tf(mu,2)                   *
*                  3,3 element of TF(mu) in Tf(mu,3)                   *
*     Tfs    : Fresnel matrix given by Eq. (27), coded as follows      *
*                  1,1 element of TF*(xmut(mu0)) in Tfs(mu0,1)         *
*                  1,2 element of TF*(xmut(mu0)) in Tfs(mu0,2)         *
*                  3,3 element of TF*(xmut(mu0)) in Tfs(mu0,3)         *
*     Rfs    : same as Tfs, but for reflection from below              *
*     xm     : refractive index of water with respect to air           *
*     Rm     : supermatrix for diffuse reflection of the atmosphere-   *
*              interface system                                        *
*     Dm     : supermatrix for diffuse downward radiation just above   *
*              the interface                                           *
*     Vm     : supermatrix defined in Eq. (90)                         *
*     Wm     : supermatrix defined in Eq. (91)                         *
*     Rm1    : first order scattering contribution to Rm               *
*     Dm1    : first order scattering contribution to Dm               *
*     Vm1    : first order scattering contribution to Vm               *
*     Wm1    : first order scattering contribution to Wm               *
*  On exit:                                                            *
*     Riflux  : reflected flux (see below)                             *
*     Diflux  : downward flux just above interface (see below)         *
*     Tiflux  : transmitted flux (see below)                           *
*     URUi    : spherical reflected flux (see below)                   *
*     UDUi    : spherical downward flux above interface (see below)    *
*     UTUi    : spherical transmitted flux (see below)                 *
*----------------------------------------------------------------------*
*  Using formulae in Chowdhary (1991) chapter 2 we can write the       *
*  fluxes in a computationally convenient form. First we define:       *
*                                                                      *
*  Xt(mu,mu0)  = 2 mut(mu) TF(mu) Dm(mu,mu0)                           *
*                                                                      *
*  Xu*(mu,mu0) = Vm*(mu,mu0) 2 mu0 TF*(mut(mu)) / xm**2                *
*                                                                      *
*  Xt*(mu,mu0) = Wm*(mu,mu0) 2 mu0 TF*(mut(mu0)) / xm**2               *
*                                                                      *
*  Xr*(mu,mu0) = 2 mut(mu) TF(mu) Vm*(mu,mu0) 2 mu0 TF*(mut(mu0))/xm**2*
*              = 2 mut(mu) TF(mu) Xu*(mu,mu0)                          *
*                                                                      *
*  where mut(.) is defined in Eq. (7),  TF(.) is defined in Eq. (23),  *
*        TF*(.) is defined in Eq. (27), Dm(.) is defined in Eq. (71),  *
*        Vm*(.) is defined in Eq. (90), Wm*(.)is defined in Eq. (91),  *
*  and xm is the refractive index of water with respect to air.        *
*  We readily obtain for the fluxes the following expressions.         *
*----------------------------------------------------------------------*
*                               1                                      *
*  UR(mu0)          =        2 int dmu mu Rm(mu,mu0)                   *
*                               0                                      *
*                                                                      *
*                               1                                      *
*  RU(mu)           =        2 int dmu0 mu0 Rm(mu,mu0)                 *
*                               0                                      *
*                                                                      *
*                               1            1                         *
*  URUi              =        2 int dmu mu 2 int dmu0 mu0 Rm(mu,mu0)   *
*                               0            0                         *
*----------------------------------------------------------------------*
*                               1                                      *
*  UT(mu0)          =        2 int dmu mu Xt(mu,mu0)                   *
*                               0                                      *
*                                                                      *
*                               1                                      *
*  TU(mut(mu))      =  xm**2 2 int dmu0 mu0 Xt(mu,mu0)                 *
*                               0                                      *
*                                                                      *
*                               1            1                         *
*  UTUi              =        2 int dmu mu 2 int dmu0 mu0 Xt(mu,mu0)   *
*                               0            0                         *
*----------------------------------------------------------------------*
*                               1                                      *
*  UD(mu0)          =        2 int dmu mu Dm(mu,mu0)                   *
*                               0                                      *
*                                                                      *
*                               1                                      *
*  DU(mu)           =        2 int dmu0 mu0 Dm(mu,mu0)                 *
*                               0                                      *
*                                                                      *
*                               1            1                         *
*  UDU              =        2 int dmu mu 2 int dmu0 mu0 Dm(mu,mu0)    *
*                               0            0                         *
*----------------------------------------------------------------------*
*                               1                                      *
*  U(U*)(mut(mu0))  =  xm**2 2 int dmu mu Xu*(mu,mu0)                  *
*                               0                                      *
*                                                                      *
*                               1                                      *
*  (U*)U(mu)        =        2 int dmu0 mu0 Xu*(mu,mu0)                *
*                               0                                      *
*                                                                      *
*                               1            1                         *
*  U(U*)U           =        2 int dmu mu 2 int dmu0 mu0 Xu*(mu,mu0)   *
*                               0            0                         *
*----------------------------------------------------------------------*
*                                                                      *
*                               1                                      *
*  U(T*)(mut(mu0))  =  xm**2 2 int dmu mu Xt*(mu,mu0)                  *
*                               0                                      *
*                                                                      *
*                               1                                      *
*  (T*)U(mu)        =        2 int dmu0 mu0 Xt*(mu,mu0)                *
*                               0                                      *
*                                                                      *
*                               1            1                         *
*  U(T*)U           =        2 int dmu mu 2 int dmu0 mu0 Xt*(mu,mu0)   *
*                               0            0                         *
*----------------------------------------------------------------------*
*                               1                                      *
*  U(R*)(mut(mu0))  =  xm**2 2 int dmu mu Xr*(mu,mu0)                  *
*                               0                                      *
*                                                                      *
*                               1                                      *
*  (R*)U(mut(mu))   =  xm**2 2 int dmu0 mu0 Xr*(mu,mu0)                *
*                               0                                      *
*                                                                      *
*                               1            1                         *
*  U(R*)U           =        2 int dmu mu 2 int dmu0 mu0 Xr*(mu,mu0)   *
*                               0            0                         *
*----------------------------------------------------------------------*
*  Fluxes are (nmat X nmat) matrices. At least half of the elements is *
*  zero, but we do not want to use this fact to devise a more compact  *
*  (and complicated) way to store fluxes; the memory penalty is small. *
*                                                                      *
*  We also calculate for each flux the zero-th and first order         *
*  scattering contributions.                                           *
*                                                                      *
*  The fluxes that have an argument are calculated for all geometries. *
*  We use the following datastructure to store them:                   *
*                                                                      *
*     Riflux( k,l, igeom, muarg, illum, iorder ),                      *
*     Diflux( k,l, igeom, muarg, illum, iorder ), and                  *
*     Tiflux( k,l, igeom, muarg, illum, iorder ).                      *
*                                                                      *
*  Here Riflux is a reflected flux, Tiflux a transmitted flux, and     *
*  Diflux a downward flux just above the interface. The meaning of the *
*  indices is as follows:                                              *
*                                                                      *
*     index    values           meaning                                *
*   -----------------------------------------------------------------  *
*     k,l      1,...,nmat       indices in a (nmat X nmat) matrix      *
*                                                                      *
*     igeom    1,...,ngeom      index of the geometry considered       *
*                                                                      *
*     muarg    Nmu0             integral over outgoing direction       *
*              Nmu              integral over incoming direction       *
*                                                                      *
*    illum     Nabove           illumination from above                *
*              Nbelow           illumination from below                *
*                                                                      *
*    iorder    Nzero            only zero order scattering             *
*              Nfirst           only first order scattering            *
*              Nall             all orders of scattering (including 0) *
*   -----------------------------------------------------------------  *
*                                                                      *
*  The values Nmu0, Nmu, Nabove, Nbelow, Nzero, Nfirst and Nall  are   *
*  coded in a parameter statement.                                     *
*  For the spherical fluxes URUi, UDUi and UTUi we use the following   *
*  datastructure:                                                      *
*                                                                      *
*      URUi( k,l, illum, iorder ),                                     *
*      UDUi( k,l, illum, iorder ), and                                 *
*      UTUi( k,l, illum, iorder ).                                     *
*                                                                      *
*  The meaning of the indices is the same as before.                   *
*                                                                      *
*  REMARK ABOUT ARGUMENTS                                              *
*  Arguments of ad_fluxes that pertain to directions below the         *
*  interface are always in the form xmut(mu) which is defined in       *
*  Eq. (7). For example we have                                        *
*                                                                      *
*    Tiflux( .., .., .., muarg, Nabove, .. ) = TU(xmut(mu))            *
*                                                                      *
*  Physically this corresponds to isotropic illumination from above,   *
*  and a detector below the interface looking in a direction with      *
*  a cosine of zenith angle given by xmut(mu). (No water medium is     *
*  assumed present.)                                                   *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter( Nmu0=1,   Nmu=2
     +         , Nzero=1,  Nfirst=2 , Nall=3
     +         , Nabove=1, Nbelow=2 )
      dimension xmu(NDmu), ebmu(NDmu), smf(NDmu) 
     +        , imu(NDgeom), imu0(NDgeom), xmut(NDmu)
     +        , Tf(NDmu,3), Tfs(NDmu,3), Rf(NDmu,3), Rfs(NDmu,3)
     +        , Rm(NDsup,NDsup),  Dm(NDsup,NDsup)
     +        , Vm(NDsup,NDsup),  Wm(NDsup,NDsup)
     +        , Rm1(NDsup,NDsup), Dm1(NDsup,NDsup)
     +        , Vm1(NDsup,NDsup), Wm1(NDsup,NDsup)
      dimension X(NDsup,NDsup)
      dimension Riflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , Diflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , Tiflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , URUi(  4,4,              Nbelow, Nall)
     +        , UDUi(  4,4,              Nbelow, Nall)
     +        , UTUi(  4,4,              Nbelow, Nall)
*----------------------------------------------------------------------*
* Initialize everything to zero                                        *
*----------------------------------------------------------------------*
      do 600 iorder=Nzero, Nall
          do 500 illum=Nabove, Nbelow
              do 400 l=1, 4
                  do 300 k=1, 4
                      do 200 igeom=1, NDgeom
                          do 100 muarg=Nmu0, Nmu
                            Riflux(k,l,igeom,muarg,illum,iorder) = 0.0d0
                            Diflux(k,l,igeom,muarg,illum,iorder) = 0.0d0
                            Tiflux(k,l,igeom,muarg,illum,iorder) = 0.0d0
  100                     continue
  200                 continue
                      URUi(k,l,illum,iorder) = 0.0d0
                      UDUi(k,l,illum,iorder) = 0.0d0
                      UTUi(k,l,illum,iorder) = 0.0d0
  300             continue
  400         continue
  500     continue
  600 continue
*----------------------------------------------------------------------*
*  Reflection for illumination from above                              *
*----------------------------------------------------------------------*
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Rm, Nabove, Nall,  Riflux, URUi )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Rm1, Nabove, Nfirst,  Riflux, URUi )
*----------------------------------------------------------------------*
*  Downward radiation above interface for illumination from above      *
*----------------------------------------------------------------------*
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Dm, Nabove, Nall,  Diflux, UDUi )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , Dm1, Nabove, Nfirst,  Diflux, UDUi )
*----------------------------------------------------------------------*
*  Transmission for illumination from above (use X = Xt)               *
*----------------------------------------------------------------------*
      call lTprod(NDmu, NDsup,  X, Tf, Dm, xmut, nmat, nmutot )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , X, Nabove, Nall,  Tiflux, UTUi )
      call lTprod(NDmu, NDsup,  X, Tf, Dm1, xmut, nmat, nmutot )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , X, Nabove, Nfirst,  Tiflux, UTUi )
*----------------------------------------------------------------------*
*  Transmission for illumination from below (use X = Xt*)              *
*----------------------------------------------------------------------*
      call rTprod(NDmu, NDsup,  X, Tfs, Wm, xmu, xm, nmat, nmutot )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , X, Nbelow, Nall,  Tiflux, UTUi )
      call rTprod(NDmu, NDsup,  X, Tfs, Wm1, xmu, xm, nmat, nmutot )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , X, Nbelow, Nfirst,  Tiflux, UTUi )
*----------------------------------------------------------------------*
*  Downward radiation above interface for illum from below (use X=Xu*) *
*----------------------------------------------------------------------*
      call rTprod(NDmu, NDsup,  X, Tfs, Vm, xmu, xm, nmat, nmutot )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , X, Nbelow, Nall,  Diflux, UDUi )
*----------------------------------------------------------------------*
*  Reflection for illumination from below (use X = Xr*)                *
*----------------------------------------------------------------------*
      call lTprod(NDmu, NDsup,  X, Tf, X, xmut, nmat, nmutot )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , X, Nbelow, Nall,  Riflux, URUi )
*----------------------------------------------------------------------*
*  First order Downward radiation above interface for illumination     *
*  from below (use X=Xu*)                                              *
*----------------------------------------------------------------------*
      call rTprod(NDmu, NDsup,  X, Tfs, Vm1, xmu, xm, nmat, nmutot )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , X, Nbelow, Nfirst,  Diflux, UDUi )
*----------------------------------------------------------------------*
*  First order reflection for illumination from below (use X=Xu*)      *
*----------------------------------------------------------------------*
      call lTprod(NDmu, NDsup,  X, Tf, X, xmut, nmat, nmutot )
      call intang(NDmu,NDsup, NDgeom
     +           , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +           , X, Nbelow, Nfirst,  Riflux, URUi )
*----------------------------------------------------------------------*
*  Calculate direct fluxes i.e. zero order scattering                  *
*----------------------------------------------------------------------*
      call flxdiri(NDmu, NDgeom,  xmu, smf, imu, imu0
     +            , ngeom, nmug, nmat, ebmu
     +            , xmut, Tf, Tfs, Rf, Rfs, xm
     +            , Riflux, Diflux, Tiflux, URUi, UDUi, UTUi )
*----------------------------------------------------------------------*
*  Add the zero order contributions to the higher order contributions  *
*  already obtained to get the sum over all orders.                    *
*----------------------------------------------------------------------*
      do 3100 illum=Nabove, Nbelow
          do 2600 l=1, nmat
              do 2500 k=1, nmat
                  do 2800 muarg=Nmu0, Nmu
                      do 2700 igeom=1, ngeom
                          Riflux(k,l,igeom,muarg,illum,Nall) =
     +                               Riflux(k,l,igeom,muarg,illum,Nall )
     +                             + Riflux(k,l,igeom,muarg,illum,Nzero)
                          Diflux(k,l,igeom,muarg,illum,Nall) =
     +                               Diflux(k,l,igeom,muarg,illum,Nall )
     +                             + Diflux(k,l,igeom,muarg,illum,Nzero)
                          Tiflux(k,l,igeom,muarg,illum,Nall) =
     +                               Tiflux(k,l,igeom,muarg,illum,Nall )
     +                             + Tiflux(k,l,igeom,muarg,illum,Nzero)
 2700                 continue
 2800             continue
                  URUi(k,l,illum,Nall) = URUi(k,l,illum,Nall) 
     +                                 + URUi(k,l,illum,Nzero)
                  UDUi(k,l,illum,Nall) = UDUi(k,l,illum,Nall) 
     +                                 + UDUi(k,l,illum,Nzero)
                  UTUi(k,l,illum,Nall) = UTUi(k,l,illum,Nall) 
     +                                 + UTUi(k,l,illum,Nzero)
 2500         continue
 2600     continue
 3100 continue
      return
      end
      subroutine flxdiri(NDmu, NDgeom, xmu, smf, imu, imu0
     +                  , ngeom, nmug, nmat, ebmu
     +                  , xmut, Tf, Tfs, Rf, Rfs, xm
     +                  , Riflux, Diflux, Tiflux, URUi, UDUi, UTUi )
*----------------------------------------------------------------------*
*  Calculate the zero order (i.e. direct) fluxes for the atmosphere    *
*  interface system, and store them in arrays Riflux, ..., UTUi.       *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter( Nmu0=1,   Nmu=2
     +         , Nzero=1,  Nall=3
     +         , Nabove=1, Nbelow=2 )
      dimension xmu(NDmu), ebmu(NDmu), smf(NDmu)
     +        , imu(NDgeom), imu0(NDgeom), xmut(NDmu)
     +        , Tf(NDmu,3), Tfs(NDmu,3), Rf(NDmu,3), Rfs(NDmu,3)
      dimension Riflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , Diflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , Tiflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , URUi(  4,4,              Nbelow, Nall)
     +        , UDUi(  4,4,              Nbelow, Nall)
     +        , UTUi(  4,4,              Nbelow, Nall)
*----------------------------------------------------------------------*
*  Calculate fluxes with one angular argument, for each geometry       *
*----------------------------------------------------------------------*
      do 800 igeom=1, ngeom
          mu0 = imu0(igeom)
          mu  = imu(igeom)
*----------------------------------------------------------------------*
*         Reflection for illumination from above                       *
*----------------------------------------------------------------------*
          fac = ebmu(mu0) * 2.0d0*xmu(mu0) * ebmu(mu0)
          call dirmat(NDmu,NDgeom
     +              , Riflux, igeom, Nmu0, Nabove, mu0, nmat, fac, Rf )
          fac = ebmu(mu)  * 2.0d0*xmu(mu)  * ebmu(mu)
          call dirmat(NDmu, NDgeom
     +              , Riflux, igeom, Nmu , Nabove, mu , nmat, fac, Rf )
*----------------------------------------------------------------------*
*         Transmission for illumination from above                     *
*----------------------------------------------------------------------*
          fac = 2.0d0*xmut(mu0) * ebmu(mu0)
          call dirmat(NDmu,NDgeom
     +              ,  Tiflux, igeom, Nmu0, Nabove, mu0, nmat, fac, Tf )
          fac = xm**2 * 2.0d0*xmut(mu) * ebmu(mu)
          call dirmat(NDmu,NDgeom
     +              ,  Tiflux, igeom, Nmu , Nabove, mu , nmat, fac, Tf )
*----------------------------------------------------------------------*
*         Downward field above interface for illumination from above   *
*         For illumination from below there is no direct down field.   *
*----------------------------------------------------------------------*
          do 700 k=1, nmat
              Diflux(k,k,igeom,Nmu0,Nabove,Nzero) = ebmu(mu0)
              Diflux(k,k,igeom,Nmu ,Nabove,Nzero) = ebmu(mu )
  700     continue
*----------------------------------------------------------------------*
*         Reflection for illumination from below                       *
*----------------------------------------------------------------------*
          fac = 2.0d0*xmut(mu0)
          call dirmat(NDmu,NDgeom
     +              ,  Riflux, igeom, Nmu0, Nbelow, mu0, nmat, fac, Rfs)
          fac = 2.0d0*xmut(mu)
          call dirmat(NDmu,NDgeom
     +              ,  Riflux, igeom, Nmu , Nbelow, mu , nmat, fac, Rfs)
*----------------------------------------------------------------------*
*         Transmission for illumination from below                     *
*----------------------------------------------------------------------*
          fac = ebmu(mu0) * 2.0d0*xmu(mu0) 
          call dirmat(NDmu,NDgeom
     +              ,  Tiflux, igeom, Nmu0, Nbelow, mu0, nmat, fac, Tfs)
          fac = ebmu(mu)  * 2.0d0*xmu(mu) / xm**2
          call dirmat(NDmu,NDgeom
     +              ,  Tiflux, igeom, Nmu , Nbelow, mu , nmat, fac, Tfs)
  800 continue
*----------------------------------------------------------------------*
*  Calculate the spherical fluxes URUi, UDUi and UTUi                  *
*----------------------------------------------------------------------*
      do 2100 mu=1, nmug
*----------------------------------------------------------------------*
*         Reflection for illumination from above                       *
*----------------------------------------------------------------------*
        fac = smf(mu)**2 * ebmu(mu)  * 2.0d0*xmu(mu)  * ebmu(mu)
        URUi(1,1,Nabove,Nzero) = URUi(1,1,Nabove,Nzero) + fac*Rf(mu,1)
        if (nmat .ge. 3) then
          URUi(1,2,Nabove,Nzero) = URUi(1,2,Nabove,Nzero) + fac*Rf(mu,2)
          URUi(3,3,Nabove,Nzero) = URUi(3,3,Nabove,Nzero) + fac*Rf(mu,3)
        endif
*----------------------------------------------------------------------*
*         Transmission for illumination from above                     *
*----------------------------------------------------------------------*
        fac = smf(mu)**2 * 2.0d0*xmut(mu)  * ebmu(mu)
        UTUi(1,1,Nabove,Nzero) = UTUi(1,1,Nabove,Nzero) + fac*Tf(mu,1)
        if (nmat .ge. 3) then
          UTUi(1,2,Nabove,Nzero) = UTUi(1,2,Nabove,Nzero) + fac*Tf(mu,2)
          UTUi(3,3,Nabove,Nzero) = UTUi(3,3,Nabove,Nzero) + fac*Tf(mu,3)
        endif
*----------------------------------------------------------------------*
*         Downward field for illumination from above                   *
*----------------------------------------------------------------------*
        fac = smf(mu)**2 * ebmu(mu)
        UDUi(1,1,Nabove,Nzero) = UDUi(1,1,Nabove,Nzero) + fac 
        if (nmat .ge. 3) UDUi(3,3,Nabove,Nzero) = UDUi(1,1,Nabove,Nzero)
*----------------------------------------------------------------------*
*         Transmission for illumination from below                     *
*----------------------------------------------------------------------*
        fac = smf(mu)**2 * ebmu(mu) * 2.0d0*xmu(mu) / xm**2
        UTUi(1,1,Nbelow,Nzero) = UTUi(1,1,Nbelow,Nzero) + fac*Tfs(mu,1)
        if (nmat .ge. 3) then
          UTUi(1,2,Nbelow,Nzero)= UTUi(1,2,Nbelow,Nzero) + fac*Tfs(mu,2)
          UTUi(3,3,Nbelow,Nzero)= UTUi(3,3,Nbelow,Nzero) + fac*Tfs(mu,3)
        endif
 2100 continue
*----------------------------------------------------------------------*
*  Reflection for illumination from below is difficult due to          *
*  integration beyond the critical angle; do it in separate routine.   *
*----------------------------------------------------------------------*
      call sstari(NDmu, xmu, smf, xmut, Rfs, nmug, nmat, xm, URUi )
*----------------------------------------------------------------------*
*  Use symmetry of Fresnel matrices to fill in rest of URUi UDUi UTUi  *
*----------------------------------------------------------------------*
      do 2200 illum=Nabove, Nbelow
          if (nmat .ge. 3) then
              URUi(2,2,illum,Nzero) = URUi(1,1,illum,Nzero)
              UDUi(2,2,illum,Nzero) = UDUi(1,1,illum,Nzero)
              UTUi(2,2,illum,Nzero) = UTUi(1,1,illum,Nzero)

              URUi(2,1,illum,Nzero) = URUi(1,2,illum,Nzero)
              UDUi(2,1,illum,Nzero) = UDUi(1,2,illum,Nzero)
              UTUi(2,1,illum,Nzero) = UTUi(1,2,illum,Nzero)
              if (nmat .eq. 4) then
                  URUi(4,4,illum,Nzero) = URUi(3,3,illum,Nzero)
                  UDUi(4,4,illum,Nzero) = UDUi(3,3,illum,Nzero)
                  UTUi(4,4,illum,Nzero) = UTUi(3,3,illum,Nzero)
              endif
          endif
 2200 continue
      return
      end
      subroutine Fmatri(NDmu, xmu, nmutot, xm, Rf, Tf, Tfs, Rfs, xmut )
*----------------------------------------------------------------------*
*  Calculate Fresnel matrices and xmut = mut(mu) for all mu values.    *
*  Equation numbers refer to Chowdhary (1991).                         *
*  The Fresnel matrices are given by Eqs. (22), (23), (26) and (27).   *
*  The Fresnel matrix for reflection from below is not calculated,     *
*  because it is not needed. The factor mut(mu) is given by Eq. (7).   *
*                                                                      *
*  On entry:                                                           *
*      xmu    : array containing mu values                             *
*      nmutot : number of mu values                                    *
*      xm     : index of refraction of the water medium w.r.t. air     *
*                                                                      *
*  On exit:                                                            *
*      Rf     : Fresnel matrix given by Eq. (22), coded as follows     *
*                  1,1 element of RF(mu) in Rf(mu,1)                   *
*                  1,2 element of RF(mu) in Rf(mu,2)                   *
*                  3,3 element of RF(mu) in Rf(mu,3)                   *
*               other elements can be found via symmetries of Eq. (22) *
*      Tf     : Fresnel matrix for transmission given by Eq. (23)      *
*               coded as Rf                                            *
*      Tfs    : Fresnel matrix for transmission from below, Eq. (27)   *
*               but with different argument: TF(mut(mu)) !!!           *
*               The reason is that we want Gauss points ABOVE the      *
*               interface to add atmosphere and interface.             *
*               Tfs is coded as Rf.                                    *
*      Rfs    : same as Tfs, but for reflection from below, Eq. (26)   *
*      xmut   : array with factors mut(mu), Eq. (7)                    *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), xmut(NDmu) 
     +        , Rf(NDmu,3), Tf(NDmu,3), Tfs(NDmu,3), Rfs(NDmu,3)
*----------------------------------------------------------------------*
*  Loop over all mu values                                             *
*----------------------------------------------------------------------*
      do 100 i=1, nmutot
*----------------------------------------------------------------------*
*         calculate mut(mu) by Eq. (7)                                 *
*----------------------------------------------------------------------*
          emut   = dsqrt(1.d0 + (xmu(i)**2-1.0d0)/xm**2)
          xmut(i) = emut
*----------------------------------------------------------------------*
*         calculate Fresnel coefficients using Eqs. (3)-(6) and        *
*         (8)-(11) Use Eq. (12) so that muts(mut(mu)) = mu !!          *
*----------------------------------------------------------------------*
          tl  =     2.0d0*xmu(i)  /(xm*xmu(i) + emut     )
          tr  =     2.0d0*xmu(i)  /(xmu(i)    + xm*emut  )
          rl  = (xm*xmu(i) - emut)/(xm*xmu(i) + emut     )
          rr  = (xmu(i) - xm*emut)/(xmu(i)    + xm*emut  )
          tls =    2.0d0*xm*emut  /(emut      + xm*xmu(i))
          trs =    2.0d0*xm*emut  /(xm*emut   + xmu(i)   )
          rls = (emut - xm*xmu(i))/(emut      + xm*xmu(i))
          rrs = (xm*emut - xmu(i))/(xm*emut   + xmu(i)   )
*----------------------------------------------------------------------*
*         Calculate the (1,1), (1,2) and (3,3) elements of Fresnel     *
*         matrices using Eqs. (22), (23), (26) and (27).               *
*----------------------------------------------------------------------*
          emu4   = 4.d0*xmu(i)
          Rf(i,1)  = (rl**2 + rr**2)/emu4
          Rf(i,2)  = (rl**2 - rr**2)/emu4
          Rf(i,3)  = 2.d0*rl*rr/emu4

          Tf(i,1)  = xm*(tl**2 + tr**2)/emu4
          Tf(i,2)  = xm*(tl**2 - tr**2)/emu4
          Tf(i,3)  = xm*2.d0*tl*tr/emu4
          
          emut4m = 4.d0*emut*xm
          Tfs(i,1)  = (tls**2 + trs**2)/emut4m
          Tfs(i,2)  = (tls**2 - trs**2)/emut4m
          Tfs(i,3)  = 2.d0*tls*trs/emut4m

          emut4  = 4.0d0*emut
          Rfs(i,1)  = (rls**2 + rrs**2)/emut4
          Rfs(i,2)  = (rls**2 - rrs**2)/emut4
          Rfs(i,3)  = 2.0d0*rls*rrs/emut4
  100 continue
*----------------------------------------------------------------------*
*  End of loop over mu values.                                         *
*----------------------------------------------------------------------*
      return
      end
      subroutine ad_gauleg( NDx, ngauss, a, b, x, w )
*----------------------------------------------------------------------*
*   Given the lower and upper limits of integration a and b, and given *
*   the number of Gauss-Legendre points ngauss, this routine returns   *
*   through array x the abscissas and through array w the weights of   *
*   the Gauss-Legendre quadrature formula. Eps is the desired accuracy *
*   of the abscissas. This routine is documented further in :          *
*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)   *
*   page 125. ISBN 0-521-30811-9.                                      *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter ( pi=3.1415926535897932384D0 )
      parameter( eps = 1.D-12 )
      dimension x(NDx), w(NDx)
      if (ngauss .gt. NDx) then
          print *,' ad_gauleg: arrays to small for ngauss =',ngauss
     +           ,' maximum NDx =',NDx
          stop 'in ad_gauleg too many Gauss points'
      endif
      m  = (ngauss+1)/2
      xm = 0.5D0*(a+b)
      xl = 0.5D0*(b-a)
      do 12 i=1, m
*         THIS IS A REALLY CLEVER ESTIMATE :
          z = dcos(pi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
    1     continue
              P1 = 1.D0
              P2 = 0.D0
              do 11 j=1, ngauss
                  P3 = P2
                  P2 = P1
                  P1 = (dble(2*j-1)*z*P2 - dble(j-1)*P3)/dble(j)
   11         continue
              Pa = dble(ngauss)*(z*P1-P2)/(z*z-1.D0)
              z1 = z
              z  = z1-P1/Pa
          if (dabs(z-z1) .gt. eps) goto 1
          x(i)          = xm-xl*z
          x(ngauss+1-i) = xm+xl*z
          w(i)          = 2.D0*xl/((1.D0-z*z)*Pa*Pa)
          w(ngauss+1-i) = w(i)
   12 continue
      return
      end
      subroutine initi( NDmu, NDlay, NDgeom
     +                , xmu, imu, imu0, nmutot, ngeom, nmat, nlayer, xm
     +                , R1, T1, R1st, T1st, e, e0
     +                , Rf, Tf, Tfs, Rfs, xmut
     +                , R1i, D1i, Ri, Di )
*----------------------------------------------------------------------*
*  Initiatilize the matrices for reflected radiation and downward      *
*  radiation just above the interface with the value for               *
*  first order scattering.                                             *
*  Equation numbers refer to Chowdhary (1991).                         *
*  On entry:                                                           *
*      xmu     : all different mu values (integration and extra points)*
*      imu     : index where mu can be found in array xmu for each     *
*                geometry                                              *
*      imu0    : index where mu0 can be found in array xmu for each    *
*                geometry                                              *
*      nmutot  : total number of mu-values                             *
*      ngeom   : number of combinations of viewing and incident        *
*                directions (number of geometries)                     *
*      nmat    : number of elements of the Stokes vector taken into    *
*                account (4 = full polarization, 3 = 3x3 approximation,*
*                2 = illegal, 1 = scalar)                              *
*      nlayer  : number of layers in the atmosphere                    *
*      xm      : index of refraction of water w.r.t. air               *
*      R1      : first order reflection matrix of the atmosphere       *
*      T1      : first order transmission matrix of the atmosphere     *
*      R1st    : same as R1 but for illumination from below            *
*      T1st    : same as T1 but for illumination from below            *
*      e       : exp(-b/mu) for each layer and for each geometry       *
*      e0      : exp(-b/mu0) for each layer and for each geometry      *
*  On exit:                                                            *
*      R1i     : first order reflection matrix of the atmosphere-      *
*                interface system for each geometry                    *
*      D1i     : first order matrix for the downard radiation just     *
*                above the interface                                   *
*      Ri      : reflection matrix of the atmosphere-interface system  *
*                for each geometry initialized to first order          *
*      Di      : matrix for the downard radiation just above interface *
*                for each geometry initialized to first order          *
*      Rf      : Fresnel matrix given by Eq. (22) coded as follows     *
*                   1,1 element of RF(mu) in Rf(mu,1)                  *
*                   1,2 element of RF(mu) in Rf(mu,2)                  *
*                   3,3 element of RF(mu) in Rf(mu,3)                  *
*                for all possible mu values                            *
*      Tf     : Fresnel matrix for transmission given by Eq. (23)      *
*               coded as Rf                                            *
*      Tfs    : Fresnel matrix for transmission from below, Eq. (27)   *
*               but with different argument: TF(mut(mu)) !!!           *
*               The reason is that we want Gauss points ABOVE the      *
*               interface to add atmosphere and interface.             *
*               Tfs is coded as Rf.                                    *
*      Rfs    : same as Tfs, but for reflection from below             *
*      xmut   : array with factors mut(mu), Eq. (7)                    *
*----------------------------------------------------------------------*
      implicit double precision(a-h,o-z)
      dimension xmu(NDmu), imu(NDgeom), imu0(NDgeom)
     +        , R1(4,4,NDgeom),   T1(4,4,NDgeom)
     +        , R1st(4,4,NDgeom), T1st(4,4,NDgeom)
     +        , e(NDlay,NDgeom),  e0(NDlay,NDgeom)
     +        , Rf(NDmu,3), Tf(NDmu,3), Tfs(NDmu,3), Rfs(NDmu,3)
     +        , xmut(NDmu)
      dimension R1i(4,4,NDgeom), D1i(4,4,NDgeom)
     +        , Ri(4,4,NDgeom),  Di(4,4,NDgeom)
*----------------------------------------------------------------------*
*  Except first order scattering as described by De Haan et al. (1987)?*
*----------------------------------------------------------------------*
      logical except
      except = .true.
*----------------------------------------------------------------------*
*  If we except first order scattering, initialize Ri and Di on first  *
*  order contribution, otherwise set Ri and Di to zero.                *
*----------------------------------------------------------------------*
      call Fmatri( NDmu, xmu, nmutot, xm, Rf, Tf, Tfs, Rfs, xmut )
      if (except) then
          call firsti( NDmu, NDlay, NDgeom 
     +           , xmu, imu, imu0, ngeom, nmat, nlayer, e, e0, Rf
     +           , R1, T1, R1st, T1st, R1i, D1i )
          do 400 igeom=1, ngeom
              do 300 k=1, nmat
                  do 200 l=1, nmat
                      Ri(l,k,igeom) = R1i(l,k,igeom)
                      Di(l,k,igeom) = D1i(l,k,igeom)
  200             continue
  300         continue
  400     continue
      else
          do 700 igeom=1, ngeom
              do 600 k=1, nmat
                  do 500 l=1, nmat
                      Ri(l,k,igeom) = 0.0d0
                      Di(l,k,igeom) = 0.0d0
  500             continue
  600         continue
  700     continue
      endif
      return
      end
      subroutine initRT(NDmu, NDlay, NDgeom, NDcoef
     +                 , coefs, ncoefs, xmu, imu0, imu
     +                 , phi, ngeom, nmat
     +                 , a, b, nlayer
     +                 , R1, T1, R1st, T1st, e, e0
     +                 , R , T , Rst , Tst )
*----------------------------------------------------------------------*
*  Initialize the reflection and transmission matrix for each          *
*  geometry with first order scattering.                               *
*  On entry:                                                           *
*     coefs   : expansion coefficients of the scattering matrix in     *
*               generalized spherical functions                        *
*     ncoefs  : order of the highest nonzero term in the expansion     *
*     xmu     : all different mu values (integration and extra points) *
*     imu0    : array indicating, for each geometry, the               *
*               position in array xmu where the mu0 value can be found *
*     phi     : azimuth of the emerging direction relative to the      *
*               direction of incidence in degrees                      *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     a       : single scattering albedo for each layer                *
*     b       : optical thickness for each layer                       *
*     nlayer  : number of layers in the atmosphere                     *
*  On exit:                                                            *
*     R1      : first order reflection matrix of the atmosphere        *
*     T1      : first order transmission matrix of the atmosphere      *
*     R1st    : same as R1 but for illumination from below             *
*     T1st    : same as T1 but for illumination from below             *
*     e       : factor exp(-b/mu) for each layer                       *
*     e0      : factor exp(-b/mu0) for each layer                      *
*     R       : equal to R1                                            *
*     T       : equal to T1                                            *
*     Rst     : equal to R1st                                          *
*     Tst     : equal to T1st                                          *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension coefs(4,4,0:NDcoef,NDlay), ncoefs(NDlay)
     +        , xmu(NDmu), imu(NDgeom), imu0(NDgeom), phi(NDgeom)
     +        , a(NDlay), b(NDlay)
      dimension e(NDlay,NDgeom), e0(NDlay,NDgeom)
     +        , R1(4,4,NDgeom), R1st(4,4,NDgeom)
     +        , T1(4,4,NDgeom), T1st(4,4,NDgeom)
     +        , T(4,4,NDgeom), Tst(4,4,NDgeom)
     +        , R(4,4,NDgeom), Rst(4,4,NDgeom)
*----------------------------------------------------------------------*
*  Except first order scattering as described by De Haan et al. (1987)?*
*----------------------------------------------------------------------*
      logical except 
      except = .true.
*----------------------------------------------------------------------*
*  If we except first order scattering, initialize R, T, Rst and Tst   *
*  on first order contribution, otherwise set them to zero.            *
*----------------------------------------------------------------------*
      if (except) then
          call first( NDmu, NDlay,NDgeom, NDcoef
     +        , coefs, ncoefs, xmu, imu0, imu
     +        , phi, ngeom 
     +        , nmat
     +        , a, b, nlayer 
     +        , R1, T1, R1st, T1st, e, e0 )
          do 300 igeom=1, ngeom
              do 200 k=1, nmat
                  do 100 l=1, nmat
                      R(l,k,igeom) = R1(l,k,igeom)
                      T(l,k,igeom) = T1(l,k,igeom)
                      Rst(l,k,igeom) = R1st(l,k,igeom)
                      Tst(l,k,igeom) = T1st(l,k,igeom)
  100             continue
  200         continue
  300     continue
      else
          do 600 igeom=1, ngeom
              do 500 k=1, nmat
                  do 400 l=1, nmat
                      R(l,k,igeom) = 0.0d0
                      T(l,k,igeom) = 0.0d0
                      Rst(l,k,igeom) = 0.0d0
                      Tst(l,k,igeom) = 0.0d0
  400             continue
  500         continue
  600     continue
      endif
      return
      end
      subroutine intang( NDmu, NDsup, NDgeom
     +                 , xmu, smf, imu, imu0, ngeom, nmug, nmat
     +                 , Am, illum, iorder,  Aflux, UAU )
*----------------------------------------------------------------------*
*  INTegrate over ANGles to obtain fluxes.                             *
*  On entry:                                                           *
*     xmu     : all different mu values (integration and extra points) *
*     smf     : supermatrix factors dsqrt(2*w*mu), or 1 for extra pts. *
*     imu     : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu value can be found  *
*     imu0    : array filled by setmu with, for each geometry, the     *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*     nmug    : number of Gauss points for the mu integrations         *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     Am      : supermatrix to be integrated                           *
*     illum,                                                           *
*     iorder  : indices that indicate where to store the calculated    *
*               fluxes in arrays Aflux and UAU (see below)             *
*               The precise meaning of illum and iorder is irrelevant. *
*  On exit:                                                            *
*     Aflux   : fluxes with only one angular integral (see below)      *
*     UAU     : spherical flux, i.e. with two angular integrals        *
*               (see below)                                            *
*                                                                      *
*  Calculate fluxes by integrating the matrix Am. Hence fluxes are     *
*  (nmat X nmat) matrices. At least half of the elements is zero, but  *
*  we do not want to use this fact to devise a more compact (and       *
*  complicated) way to store fluxes; the memory penalty is small.      *
*                                                                      *
*  The m=0 Fourier component of the matrix is Am(mu,mu0).              *
*  We define the following fluxes as integrals over the incident       *
*  direction or the outgoing direction or both:                        *
*                 1                                                    *
*     AU(mu)  = 2 int dmu0 mu0 Am(mu,mu0)                              *
*                 0                                                    *
*                                                                      *
*                 1                                                    *
*     UA(mu0) = 2 int dmu mu Am(mu,mu0)                                *
*                 0                                                    *
*                                                                      *
*                 1            1                                       *
*     UAU    = 2 int dmu mu 2 int dmu0 mu0 Am(mu,mu0)                  *
*                 0            0                                       *
*                                                                      *
*  The fluxes that have an argument mu or mu0 are calculated for all   *
*  geometries. We use the following datastructure to store them:       *
*                                                                      *
*      Aflux( k,l, igeom, muarg, illum, iorder ).                      *
*                                                                      *
*  The meaning of the indices is as follows:                           *
*                                                                      *
*     index    values           meaning                                *
*   -----------------------------------------------------------------  *
*     k,l      1,...,nmat       indices in a (nmat X nmat) matrix      *
*                                                                      *
*     igeom    1,...,ngeom      index of the geometry considered       *
*                                                                      *
*     muarg    Nmu0             integral over outgoing direction, so   *
*                               we have UA(mu0)                        *
*              Nmu              integral over incoming direction, so   *
*                               we have AU(mu)                         *
*                                                                      *
*    illum     irrelevant       irrelevant                             *
*                                                                      *
*    iorder    irrelevant       irrelevant                             *
*   -----------------------------------------------------------------  *
*                                                                      *
*  The values Nmu0 and Nmu  are coded in a parameter statement.        *
*  For the spherical flux UAU and UTU we use the following             *
*  datastructure:                                                      *
*                                                                      *
*      UAU( k,l, illum, iorder ).                                      *
*                                                                      *
*  The meaning of the indices is the same as before.                   *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter( Nmu0=1,   Nmu=2
     +         , Nbelow=2, Nall=3 )
      dimension xmu(NDmu), smf(NDmu)
     +        , imu(NDgeom), imu0(NDgeom)
     +        , Am(NDsup,NDsup)
      dimension Aflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , UAU(  4,4,              Nbelow, Nall)
*----------------------------------------------------------------------*
* Initialize everything to zero                                        *
*----------------------------------------------------------------------*
      do 400 l=1, 4
          do 300 k=1, 4
              do 200 igeom=1, NDgeom
                  do 100 muarg=Nmu0, Nmu
                     Aflux(k,l,igeom,muarg,illum,iorder) = 0.0d0
  100             continue
  200         continue
              UAU(k,l,illum,iorder) = 0.0d0
  300     continue
  400 continue
*----------------------------------------------------------------------*
* Start loop over the geometries                                       *
*----------------------------------------------------------------------*
      do 1600 igeom=1, ngeom
*----------------------------------------------------------------------*
*         calculate fluxes with argument mu0                           *
*         (i.e. integrate over outgoing directions specified by mu)    *
*----------------------------------------------------------------------*
          mu0  = imu0(igeom)
          jbase = (mu0-1)*nmat
          do 1200 mu=1, nmug
              w    = smf(mu)/smf(mu0)
              ibase = (mu-1)*nmat
              do 1100 l=1, nmat
                  do 1000 k=1, nmat
                      Aflux(k,l,igeom,Nmu0,illum,iorder) = 
     +                Aflux(k,l,igeom,Nmu0,illum,iorder) 
     +                                          + w*Am(ibase+k,jbase+l)
 1000             continue
 1100         continue
 1200     continue
*----------------------------------------------------------------------*
*         calculate fluxes with argument mu                            *
*         (i.e. integrate over incident directions specified by mu0)   *
*----------------------------------------------------------------------*
          mu   = imu(igeom)
          ibase = (mu-1)*nmat
          do 1500 mu0=1, nmug
              w    = smf(mu0)/smf(mu)
              jbase = (mu0-1)*nmat
              do 1400 l=1, nmat
                  do 1300 k=1, nmat
                      Aflux(k,l,igeom,Nmu,illum,iorder) = 
     +                Aflux(k,l,igeom,Nmu,illum,iorder) 
     +                                          + w*Am(ibase+k,jbase+l)
 1300             continue
 1400         continue
 1500     continue
 1600 continue
*----------------------------------------------------------------------*
*  End of loop over geometries                                         *
*  Calculate the spherical flux UAU                                    *
*----------------------------------------------------------------------*
      do 2100 mu=1, nmug
          ibase = (mu-1)*nmat
          do 2000 mu0=1, nmug
              jbase = (mu0-1)*nmat
              w = smf(mu)*smf(mu0)
              do 1900 l=1, nmat
                  do 1800 k=1, nmat
                      UAU(k,l,illum,iorder) = UAU(k,l,illum,iorder) 
     +                                         + w*Am(ibase+k,jbase+l)
 1800             continue
 1900         continue
 2000     continue
 2100 continue
      return
      end
      subroutine layerm(NDmu, NDsup, NDlay, NDcoef
     +                 , m, M0, M1, M2, layer, xmu, smf
     +                 , nmug, nmutot, nmat, epsilon
     +                 , coefs, ncoef, Zmplus, Zmmin, a, b
     +                 , ebtop, Rm1top, Rmtop, Tm1top, Tmtop )
*----------------------------------------------------------------------*
*  Calculate the m-th Fourier component of reflection and transmission *
*  of a homogeneous layer.                                             *
*     When  0 <= m <= M2      use doubling,                            *
*     when  M2 < m <= M1      use first plus second order scattering,  *
*     when  M1 < m <= M0      use first order scattering and           *
*     when  M0 < m            use no scattering at all.                *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension coefs(4,4,0:NDcoef,NDlay), ncoef(NDlay)
     +        , M0(NDlay), M1(NDlay), M2(NDlay)
     +        , xmu(NDmu), smf(NDmu), ebtop(NDmu)
     +        , Rmtop(NDsup,NDsup),  Rm1top(NDsup,NDsup)
     +        , Tmtop(NDsup,NDsup),  Tm1top(NDsup,NDsup)
     +        , Zmmin(NDsup,NDsup),  Zmplus(NDsup,NDsup)
      dimension ebmu(NDmu)
      logical verbo
      verbo = .false.
*----------------------------------------------------------------------*
      if (verbo)
     +    print *,' layerm: start calculating R and T for layer ', layer
      if (verbo)
     +    print *,'     M0=',M0(layer),' M1=',M1(layer),' M2=',M2(layer)
      nsup = nmutot*nmat
      if (m .gt. M0(layer)) then
          if (verbo) print *,' layerm: no scattering for m = ',m
          call zero(NDsup, Rmtop, nsup )
          call zero(NDsup, Tmtop, nsup )
          call zero(NDsup, Rm1top, nsup )
          call zero(NDsup, Tm1top, nsup )
          call expbmu(NDmu, b, xmu, nmutot, ebtop )
      else if (m .gt. M1(layer)) then
          if (verbo) print *,' layerm: one order suffices for m = ',m
          call ord1m(NDmu
     +              , NDsup, m, xmu, smf, nmutot, nmat, Zmplus, Zmmin
     +              , a, b, ebtop, Rm1top, Tm1top )
          call assign( NDsup,  Rmtop, Rm1top, nmat, nmutot )
          call assign( NDsup,  Tmtop, Tm1top, nmat, nmutot )
      else if (m .gt. M2(layer)) then
          if (verbo) print *,' layerm: two orders suffice for m = ',m
          call ord1m(NDmu
     +              , NDsup, m, xmu, smf, nmutot, nmat, Zmplus, Zmmin
     +              , a, b, ebtop, Rm1top, Tm1top )
          call assign( NDsup,  Rmtop, Rm1top, nmat, nmutot )
          call assign( NDsup,  Tmtop, Tm1top, nmat, nmutot )
          call ord2m( NDmu, NDsup
     +              , m, xmu, smf, nmug, nmutot, nmat, Zmplus, Zmmin
     +              , a, b, ebtop, Rmtop, Tmtop )
      else
          if (verbo) print *,' layerm: doubling needed for m = ',m
          call ord1m(NDmu
     +              , NDsup, m, xmu, smf, nmutot, nmat, Zmplus, Zmmin
     +              , a, b, ebtop, Rm1top, Tm1top )
          xmumin = xmu(1)
          call bstart(NDlay, NDcoef
     +               , m, layer, coefs, ncoef, M0, xmumin, a, b
     +               , epsilon
     +               , b0, ndoubl )
          call expbmu(NDmu, b0, xmu, nmutot, ebmu )
          call ord1m(NDmu
     +              , NDsup, m, xmu, smf, nmutot, nmat, Zmplus, Zmmin
     +              , a, b0, ebmu, Rmtop, Tmtop )
          call ord2m(NDmu, NDsup
     +              , m, xmu, smf, nmug, nmutot, nmat, Zmplus, Zmmin
     +              , a, b0, ebmu, Rmtop, Tmtop )
          if (ndoubl .gt. 0) then
              bb = b0
              do 100 i=1, ndoubl
                  call double( NDmu,NDsup
     +                       , Rmtop, Tmtop, ebmu, nmutot, nmug, nmat )
                  bb = 2.D0*bb
                  call expbmu(NDmu, bb, xmu, nmutot, ebmu )
  100         continue
          endif
      endif
      return
      end
      subroutine lRprod( NDmu, NDsup,  A, R, B, xmu, nmat, nmutot )
************************************************************************
*
*
*
*
*
*
*
************************************************************************
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup, NDsup), R(NDmu,3)
     +        , xmu(NDmu)

      if (nmat.eq.1) then
          do 200 mu=1, nmutot
              i = (mu-1)*nmat + 1
              x2mu = xmu(mu) + xmu(mu)
              do 100 mu0=1, nmutot
                  j = (mu0-1)*nmat + 1
                  A(i,j) = x2mu*( R(mu,1)*B(i,j) )
  100         continue
  200     continue

      else if (nmat.eq.3) then
        do 500 mu=1,nmutot
          i = (mu-1)*nmat
          x2mu = xmu(mu) + xmu(mu)
          do 400 mu0=1,nmutot
            j = (mu0-1)*nmat 
              do 300 k=1, nmat
                kj = k+j
                A(i+1,kj) = x2mu*(R(mu,1)*B(i+1,kj) + R(mu,2)*B(i+2,kj))
                A(i+2,kj) = x2mu*(R(mu,2)*B(i+1,kj) + R(mu,1)*B(i+2,kj))
                A(i+3,kj) = x2mu*(R(mu,3)*B(i+3,kj)) 
  300         continue
  400       continue
  500     continue

      else
        do 800 mu=1,nmutot
          i = (mu-1)*nmat
          x2mu = xmu(mu) + xmu(mu)
          do 700 mu0=1,nmutot
            j = (mu0-1)*nmat 
              do 600 k=1, nmat
                kj = k+j
                A(i+1,kj) = x2mu*(R(mu,1)*B(i+1,kj) + R(mu,2)*B(i+2,kj))
                A(i+2,kj) = x2mu*(R(mu,2)*B(i+1,kj) + R(mu,1)*B(i+2,kj))
                A(i+3,kj) = x2mu*(R(mu,3)*B(i+3,kj)) 
                A(i+4,kj) = x2mu*(R(mu,3)*B(i+4,kj))
  600         continue
  700       continue
  800     continue
      endif
      return
      end
      subroutine lTprod(NDmu, NDsup,  A, Tf, B, xmut, nmat, nmutot )
*----------------------------------------------------------------------*
*  Calculate the following 'left transmission product' :               *
*                                                                      *
*    A(mu,mu0) = 2 xmut(mu) TF(mu) B(mu,mu0)                           *
*                                                                      *
*  where xmut(mu) is defined in Chowdhary (1991) Eq. (7), and TF(mu)   *
*  is the Fresnel matrix defined in Chowdhary (1991) Eq. (23).         *
*                                                                      *
*  On entry:                                                           *
*     Tf     : Fresnel matrix given by Eq. (23), coded as follows      *
*                  1,1 element of TF(mu) in Tf(mu,1)                   *
*                  1,2 element of TF(mu) in Tf(mu,2)                   *
*                  3,3 element of TF(mu) in Tf(mu,3)                   *
*     B      : supermatrix B to be multiplied                          *
*     xmut   : factor defined in Eq. (7) for every mu value            *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     nmutot  : total number of distinct mu points                     *
*                                                                      *
*  On exit:                                                            *
*     A      : supermatrix given by above formula                      *
*                                                                      *
*  Note that A and B may be the same when this subroutine is called!   *
*  This is the reason for the awkward use of variable A1 (see below).  *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup, NDsup), Tf(NDmu,3), xmut(NDmu)
      if (nmat .eq. 1) then
*----------------------------------------------------------------------*
*  Treat case without polarization                                     *
*----------------------------------------------------------------------*
          do 200 mu=1, nmutot
              xmT1 = 2.0d0*xmut(mu)*Tf(mu,1)
              do 100 mu0=1, nmutot
                  A(mu,mu0) = xmT1*B(mu,mu0)
  100         continue
  200     continue
      else if (nmat .eq. 3) then
*----------------------------------------------------------------------*
*  Treat case with 3 X 3 approximation                                 *
*----------------------------------------------------------------------*
        do 500 mu=1, nmutot
          ibase = (mu-1)*nmat
          xmT1 = 2.0d0*xmut(mu)*Tf(mu,1)
          xmT2 = 2.0d0*xmut(mu)*Tf(mu,2)
          xmT3 = 2.0d0*xmut(mu)*Tf(mu,3)
          do 400 mu0=1, nmutot
            jbase = (mu0-1)*nmat
            do 300 k=1, nmat
              kj = jbase + k
              A1            = xmT1*B(ibase+1,kj) + xmT2*B(ibase+2,kj)
              A(ibase+2,kj) = xmT2*B(ibase+1,kj) + xmT1*B(ibase+2,kj)
              A(ibase+3,kj) = xmT3*B(ibase+3,kj)
              A(ibase+1,kj) = A1
  300       continue
  400     continue
  500   continue
      else
*----------------------------------------------------------------------*
*  Treat case with polarization fully included                         *
*----------------------------------------------------------------------*
        do 800 mu=1, nmutot
          ibase = (mu-1)*nmat
          xmT1 = 2.0d0*xmut(mu)*Tf(mu,1)
          xmT2 = 2.0d0*xmut(mu)*Tf(mu,2)
          xmT3 = 2.0d0*xmut(mu)*Tf(mu,3)
          do 700 mu0=1, nmutot
            jbase = (mu0-1)*nmat
            do 600 k=1, nmat
              kj = jbase + k
              A1            = xmT1*B(ibase+1,kj) + xmT2*B(ibase+2,kj)
              A(ibase+2,kj) = xmT2*B(ibase+1,kj) + xmT1*B(ibase+2,kj)
              A(ibase+3,kj) = xmT3*B(ibase+3,kj)
              A(ibase+4,kj) = xmT3*B(ibase+4,kj)
              A(ibase+1,kj) = A1
  600       continue
  700     continue
  800   continue
      endif
      return
      end
      subroutine member( e, a, NDa, n, found, index )
*----------------------------------------------------------------------*
*  Check if element e is contained in the first n positions of array a *
*  On entry :                                                          *
*      e        : element to be looked for                             *
*      a        : array to be searched                                 *
*      NDa      : dimension of array a                                 *
*      n        : number of elements of a that is checked              *
*  On exit :                                                           *
*      found    : .true. if e occurs in a, .false. otherwise           *
*      index    : position where e occurs in a if found,               *
*                 otherwise set to zero. If e occurs more than once,   *
*                 index refers to the first occurrence.                *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter( eps=1.D-8 )
      dimension a(NDa)
      logical found
      n = min0(n,NDa)
      found = .false.
      index = 0
      do 100 i=n, 1, -1
          if (dabs(e-a(i)) .lt. eps) index = i
  100 continue
      if (index .gt. 0) found = .true.
      return
      end
      subroutine newfou( NDmu, NDsup, NDgeom
     +                 , m, Rmbot, Tmbot, xRm1bot, xTm1bot, xTm1sbot
     +                 , Rmsbot, xRm1sbot
     +                 , xmu, imu, imu0, smf, phi
     +                 , nmat, ngeom
     +                 , R, T, Rst, Tst 
     +                 , R1, T1, R1st, T1st )
*----------------------------------------------------------------------*
*  Add the m-th Fourier term to the reflection and transmission        *
*  matrices for each geometry. Because first order scattering is       *
*  excepted we subtract the first order contribution to the m-th       *
*  Fourier term. We use Eqs. (127) and (128) of De Haan et al. (1987). *
*  On entry:                                                           *
*     m       : Fourier index                                          *
*     Rmbot   : Fourier component of reflection supermatrix            *
*     Tmbot   : Fourier component of transmission supermatrix          *
*    xRm1bot  : first order contribution to Rmbot                      *
*    xTm1bot  : first order contribution to Tmbot                      *
*               (the transmission for illumination from below is       *
*               obtained from symmetry relations)                      *
*    xTm1sbot : first order contribution to Fourier component of       *
*               transmission supermatrix for illumination from below   *
*     Rmsbot  : same as Rmbot but for illumination from below          *
*    xRm1sbot : first order contribution to Rmsbot                     *
*     xmu     : all different mu values (integration and extra points) *
*     imu     : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu value can be found  *
*     imu0    : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu0 value can be found *
*     smf     : supermatrix factors dsqrt(2*w*mu), or 1 for extra pts. *
*     phi     : azimuth of the emerging direction relative to the      *
*               direction of incidence in degrees                      *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*                                                                      *
*  On exit the m-th Fourier term has been added to :                   *
*     T       : transmission matrix of the atmosphere                  *
*     Tst     : same as T but for illumination from below              *
*     R       : reflection matrix of the atmosphere                    *
*     Rst     : same as R but for illumination from below              *
*  If first order is not excepted (except=.false.) the m-th Fourier    *
*  is also added to :                                                  *
*     T1      : 1st order transmission matrix of the atmosphere        *
*     T1st    : same as T1 but for illumination from below             *
*     R1      : 1st order reflection matrix of the atmosphere          *
*     R1st    : same as R1 but for illumination from below             *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter ( pi=3.1415926535897932384D0 )

      parameter ( radfac=pi/180.D0 )
      dimension Rmbot(NDsup,NDsup),  Rmsbot(NDsup,NDsup)
     +        , xRm1bot(4,4,NDgeom), xRm1sbot(4,4,NDgeom)
     +        , Tmbot(NDsup,NDsup),  Tmsbot(4,4)
     +        , xTm1bot(4,4,NDgeom), xTm1sbot(4,4,NDgeom)
     +        , xmu(NDmu), imu(NDgeom), imu0(NDgeom), phi(NDgeom)
     +        , smf(NDmu)
     +        , R(4,4,NDgeom),       Rst(4,4,NDgeom)
     +        , T(4,4,NDgeom),       Tst(4,4,NDgeom)
     +        , R1(4,4,NDgeom),      R1st(4,4,NDgeom)
     +        , T1(4,4,NDgeom),      T1st(4,4,NDgeom)
      dimension Bplus(4), Bmin(4), delmi(4), delpl(4)
*----------------------------------------------------------------------*
*  Except first order scattering as described by De Haan et al. (1987)?*
*----------------------------------------------------------------------*
      logical except
      except = .true.
*----------------------------------------------------------------------*
*  Start loop over geometries                                          *
*----------------------------------------------------------------------*
      do 1000 igeom=1, ngeom
*----------------------------------------------------------------------*
*         Calculate the diagonal B-matrices in Eqs. (48)-(49)          *
*----------------------------------------------------------------------*
          cosmph = dcos(m*phi(igeom)*radfac)
          sinmph = dsin(m*phi(igeom)*radfac)
          Bplus(1) = cosmph
          Bplus(2) = cosmph
          Bplus(3) = sinmph
          Bplus(4) = sinmph
          Bmin(1)  = -sinmph
          Bmin(2)  = -sinmph
          Bmin(3)  = cosmph
          Bmin(4)  = cosmph
*----------------------------------------------------------------------*
*         Calculate the diagonal matrices 1-DELTA and 1+DELTA          *
*----------------------------------------------------------------------*
          delmi(1) = 0.0d0
          delmi(2) = 0.0d0
          delmi(3) = 2.0d0
          delmi(4) = 2.0d0
          delpl(1) = 2.0d0
          delpl(2) = 2.0d0
          delpl(3) = 0.0d0
          delpl(4) = 0.0d0
*----------------------------------------------------------------------*
          ibase = (imu(igeom)-1)*nmat
          jbase = (imu0(igeom)-1)*nmat
*----------------------------------------------------------------------*
*         Find T* by symmetry : T* = q3 T~ q3 (also inhomogeneous!)    *
*----------------------------------------------------------------------*
          do 85 ki=1, nmat
              i = ibase+ki
              do 83 kj=1, nmat
                  j = jbase+kj
                  Tmsbot(ki,kj)  = Tmbot(j,i)
                  if (ki .eq. 3) then
                      Tmsbot(ki,kj)  = - Tmsbot(ki,kj)
                  endif
                  if (kj .eq. 3) then
                      Tmsbot(ki,kj)  = - Tmsbot(ki,kj)
                  endif
   83         continue
   85     continue
*----------------------------------------------------------------------*
*         Remove the supermatrix factors                               *
*         Insert the factor 0.5*(2-delta(m,0)) in Eqs. (127)-(128)     *
*----------------------------------------------------------------------*
          w = 1.d0/(smf(imu(igeom))*smf(imu0(igeom)))
          fac = 1.0d0
          if (m .eq. 0) fac = 0.5d0
*----------------------------------------------------------------------*
*         Except first order scattering                                *
*----------------------------------------------------------------------*
          if (except) then
            do 200 ki=1, nmat
              i = ibase+ki
              do 100 kj=1, nmat
                j = jbase+kj
                R(ki,kj,igeom) = R(ki,kj,igeom)
     +           + Bplus(ki)*fac*(w*Rmbot(i,j)
     +                                 -xRm1bot(ki,kj,igeom))*delpl(kj)
     +           + Bmin(ki) *fac*(w*Rmbot(i,j)
     +                                 -xRm1bot(ki,kj,igeom))*delmi(kj)
                T(ki,kj,igeom) = T(ki,kj,igeom)
     +           + Bplus(ki)*fac*(w*Tmbot(i,j)
     +                                 -xTm1bot(ki,kj,igeom))*delpl(kj)
     +           + Bmin(ki) *fac*(w*Tmbot(i,j)
     +                                 -xTm1bot(ki,kj,igeom))*delmi(kj) 
                Rst(ki,kj,igeom) = Rst(ki,kj,igeom)
     +           + Bplus(ki)*fac*(w*Rmsbot(i,j)
     +                                 -xRm1sbot(ki,kj,igeom))*delpl(kj)
     +           + Bmin(ki) *fac*(w*Rmsbot(i,j)
     +                                 -xRm1sbot(ki,kj,igeom))*delmi(kj)
                Tst(ki,kj,igeom) = Tst(ki,kj,igeom)
     +           + Bplus(ki)*fac*(w*Tmsbot(ki,kj)
     +                                 -xTm1sbot(ki,kj,igeom))*delpl(kj)
     +           + Bmin(ki) *fac*(w*Tmsbot(ki,kj)
     +                                 -xTm1sbot(ki,kj,igeom))*delmi(kj)
  100           continue
  200       continue
*----------------------------------------------------------------------*
*         Do not except first order scattering                         *
*----------------------------------------------------------------------*
          else
            do 400 ki=1, nmat
              i = ibase+ki
              do 300 kj=1, nmat
*----------------------------------------------------------------------*
*               sum of first and higher order scattering               *
*----------------------------------------------------------------------*
                j = jbase+kj
                R(ki,kj,igeom) = R(ki,kj,igeom)
     +                     + Bplus(ki)*fac*w*Rmbot(i,j)*delpl(kj)
     +                     + Bmin(ki) *fac*w*Rmbot(i,j)*delmi(kj)
                T(ki,kj,igeom) = T(ki,kj,igeom) 
     +                     + Bplus(ki)*fac*w*Tmbot(i,j)*delpl(kj)
     +                     + Bmin(ki) *fac*w*Tmbot(i,j)*delmi(kj)
                Rst(ki,kj,igeom) = Rst(ki,kj,igeom)
     +                     + Bplus(ki)*fac*w*Rmsbot(i,j)*delpl(kj)
     +                     + Bmin(ki) *fac*w*Rmsbot(i,j)*delmi(kj)
                Tst(ki,kj,igeom) = Tst(ki,kj,igeom)
     +                     + Bplus(ki)*fac*w*Tmsbot(ki,kj)*delpl(kj)
     +                     + Bmin(ki) *fac*w*Tmsbot(ki,kj)*delmi(kj)
*----------------------------------------------------------------------*
*               first order scattering                                 *
*----------------------------------------------------------------------*
                R1(ki,kj,igeom) = R1(ki,kj,igeom)
     +                 + Bplus(ki)*fac*xRm1bot(ki,kj,igeom)*delpl(kj)
     +                 + Bmin(ki) *fac*xRm1bot(ki,kj,igeom)*delmi(kj)
                T1(ki,kj,igeom) = T1(ki,kj,igeom) 
     +                 + Bplus(ki)*fac*xTm1bot(ki,kj,igeom)*delpl(kj)
     +                 + Bmin(ki) *fac*xTm1bot(ki,kj,igeom)*delmi(kj)
                R1st(ki,kj,igeom) = R1st(ki,kj,igeom)
     +                 + Bplus(ki)*fac*xRm1sbot(ki,kj,igeom)*delpl(kj)
     +                 + Bmin(ki) *fac*xRm1sbot(ki,kj,igeom)*delmi(kj)
                T1st(ki,kj,igeom) = T1st(ki,kj,igeom)
     +                 + Bplus(ki)*fac*xTm1sbot(ki,kj,igeom)*delpl(kj)
     +                 + Bmin(ki) *fac*xTm1sbot(ki,kj,igeom)*delmi(kj)
  300           continue
  400       continue
          endif
 1000 continue
      return
      end
      subroutine newfoui(NDmu, NDsup, NDgeom
     +                 , m, Rmbot, Tmbot, xRm1i, xDm1i 
     +                 , xmu, imu, imu0, smf, phi
     +                 , nmat, ngeom
     +                 , Ri, Di, R1i, D1i )
*----------------------------------------------------------------------*
*  Add the m-th Fourier term to the reflection matrix of the atmosphere*
*  interface system and the matrix for downward radiation just above   *
*  the interface for each geometry. Because first order scattering is  *
*  excepted we subtract the first order contribution to the m-th       *
*  Fourier term. We use Eqs. (127) and (128) of De Haan et al. (1987). *
*  On entry:                                                           *
*     m       : Fourier index                                          *
*     Rmbot   : Fourier component of reflection supermatrix            *
*     Tmbot   : Fourier component of supermatrix for downward          *
*               radiation just above the interface                     *
*    xRm1bot  : first order contribution to Rmbot                      *
*    xDm1i    : first order contribution to Tmbot                      *
*     xmu     : all different mu values (integration and extra points) *
*     imu     : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu value can be found  *
*     imu0    : array filled by setmu with, for each geometry, the     *
*               position in array xmu where the mu0 value can be found *
*     smf     : supermatrix factors dsqrt(2*w*mu), or 1 for extra pts. *
*     phi     : azimuth of the emerging direction relative to the      *
*               direction of incidence in degrees                      *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*                                                                      *
*  On exit the m-th Fourier term has been added to :                   *
*     Ri      : reflection matrix of the atmosphere interface system   *
*     Di      : matrix for downward radiation just above the           *
*               interface                                              *
*  If first order scattering is not excepted (except = .false.) then   *
*  the m-th Fourier term is also added to :                            *
*     R1i     : 1st order scattering contribution to Ri                *
*     D1i     : 1st order scattering contribution to Di                *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter ( pi=3.1415926535897932384D0 )

      parameter ( radfac=pi/180.D0 )
      dimension Rmbot(NDsup,NDsup), xRm1i(4,4,NDgeom)
     +        , Tmbot(NDsup,NDsup), xDm1i(4,4,NDgeom)
     +        , xmu(NDmu), imu(NDgeom), imu0(NDgeom), phi(NDgeom)
     +        , smf(NDmu)
     +        , Ri(4,4,NDgeom),  Di(4,4,NDgeom)
     +        , R1i(4,4,NDgeom), D1i(4,4,NDgeom)
      dimension Bplus(4), Bmin(4), delmi(4), delpl(4)
*----------------------------------------------------------------------*
*  Except first order scattering as described by De Haan et al. (1987)?*
*----------------------------------------------------------------------*
      logical except
      except = .true.
*----------------------------------------------------------------------*
*  Start loop over geometries                                          *
*----------------------------------------------------------------------*
      do 1000 igeom=1, ngeom
*----------------------------------------------------------------------*
*         Calculate the diagonal B-matrices in Eqs. (48)-(49)          *
*----------------------------------------------------------------------*
          cosmph = dcos(m*phi(igeom)*radfac)
          sinmph = dsin(m*phi(igeom)*radfac)
          Bplus(1) = cosmph
          Bplus(2) = cosmph
          Bplus(3) = sinmph
          Bplus(4) = sinmph
          Bmin(1)  = -sinmph
          Bmin(2)  = -sinmph
          Bmin(3)  = cosmph
          Bmin(4)  = cosmph
*----------------------------------------------------------------------*
*         Calculate the diagonal matrices 1-DELTA and 1+DELTA          *
*----------------------------------------------------------------------*
          delmi(1) = 0.0d0
          delmi(2) = 0.0d0
          delmi(3) = 2.0d0
          delmi(4) = 2.0d0
          delpl(1) = 2.0d0
          delpl(2) = 2.0d0
          delpl(3) = 0.0d0
          delpl(4) = 0.0d0
*----------------------------------------------------------------------*
          ibase = (imu(igeom)-1)*nmat
          jbase = (imu0(igeom)-1)*nmat
*----------------------------------------------------------------------*
*         Remove the supermatrix factors                               *
*         Insert the factor 0.5*(2-delta(m,0)) in Eqs. (127)-(128)     *
*----------------------------------------------------------------------*
          w = 1.d0/(smf(imu(igeom))*smf(imu0(igeom)))
          fac = 1.0d0
          if (m .eq. 0) fac = 0.5d0
*----------------------------------------------------------------------*
*         Except first order scattering                                *
*----------------------------------------------------------------------*
          if (except) then
            do 200 ki=1, nmat
              i = ibase+ki
              do 100 kj=1, nmat
                j = jbase+kj
                Ri(ki,kj,igeom) = Ri(ki,kj,igeom)
     +           + fac*( Bplus(ki)*(w*Rmbot(i,j)
     +                                  - xRm1i(ki,kj,igeom))*delpl(kj)
     +                 + Bmin(ki) *(w*Rmbot(i,j)
     +                                  - xRm1i(ki,kj,igeom))*delmi(kj))
                Di(ki,kj,igeom) = Di(ki,kj,igeom)
     +           + fac*( Bplus(ki)*(w*Tmbot(i,j)
     +                                  - xDm1i(ki,kj,igeom))*delpl(kj)
     +                 + Bmin(ki) *(w*Tmbot(i,j)
     +                                  - xDm1i(ki,kj,igeom))*delmi(kj))
  100         continue
  200       continue
*----------------------------------------------------------------------*
*         Do not except first order scattering                         *
*----------------------------------------------------------------------*
          else
            do 400 ki=1, nmat
              i = ibase+ki
              do 300 kj=1, nmat
                j = jbase+kj
*----------------------------------------------------------------------*
*               sum of first and higher order scattering               *
*----------------------------------------------------------------------*
                Ri(ki,kj,igeom) = Ri(ki,kj,igeom)
     +                       + fac*w*( Bplus(ki)*Rmbot(i,j)*delpl(kj)
     +                               + Bmin(ki) *Rmbot(i,j)*delmi(kj))
                Di(ki,kj,igeom) = Di(ki,kj,igeom)
     +                       + fac*w*( Bplus(ki)*Tmbot(i,j)*delpl(kj)
     +                               + Bmin(ki) *Tmbot(i,j)*delmi(kj))
*----------------------------------------------------------------------*
*               first order scattering                                 *
*----------------------------------------------------------------------*
                R1i(ki,kj,igeom) = R1i(ki,kj,igeom)
     +                  + fac*( Bplus(ki)*xRm1i(ki,kj,igeom)*delpl(kj)
     +                        + Bmin(ki) *xRm1i(ki,kj,igeom)*delmi(kj))
                D1i(ki,kj,igeom) = D1i(ki,kj,igeom)
     +                  + fac*( Bplus(ki)*xDm1i(ki,kj,igeom)*delpl(kj)
     +                        + Bmin(ki) *xDm1i(ki,kj,igeom)*delmi(kj))
  300         continue
  400       continue
          endif
 1000 continue
      return
      end
      subroutine nobot( NDmu,NDsup, Rmtop, Tmtop, Rmbot, Tmbot, Rmsbot
     +                , Rm1top, Tm1top, Rm1bot, Rm1sbot, Tm1bot
     +                , ebtop, ebbot, nmutot, nmat )
*----------------------------------------------------------------------*
*  Use adding equations when there is no scattering in the bottom layer*
*  Do this also for first and zero-th order scattering.                *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension Rmtop(NDsup,NDsup), Tmtop(NDsup,NDsup)
     +        , Rmbot(NDsup,NDsup), Tmbot(NDsup,NDsup)
     +        , Rmsbot(NDsup,NDsup)
     +        , Rm1top(NDsup,NDsup), Tm1top(NDsup,NDsup)
     +        , Rm1bot(NDsup,NDsup), Rm1sbot(NDsup,NDsup)
     +        , Tm1bot(NDsup,NDsup)
     +        , ebtop(NDmu), ebbot(NDmu), Etop(NDsup), Ebot(NDsup)
      do 200 mu=1, nmutot
          do 100 k=1, nmat
              Etop((mu-1)*nmat+k) = ebtop(mu)
              Ebot((mu-1)*nmat+k) = ebbot(mu)
  100     continue
  200 continue
*----------------------------------------------------------------------*
*  Reflection  R = R'                                                  *
*----------------------------------------------------------------------*
      call assign( NDsup,  Rmbot,  Rmtop,  nmat, nmutot )
      call assign( NDsup,  Rm1bot, Rm1top, nmat, nmutot )
*----------------------------------------------------------------------*
*  Transmission T = E"T'                                               *
*----------------------------------------------------------------------*
      call ldiapr(NDsup, Tmbot,  Ebot, Tmtop,  nmat, nmutot)
      call ldiapr(NDsup, Tm1bot, Ebot, Tm1top, nmat, nmutot)
*----------------------------------------------------------------------*
*  Reflection star R* = E"R'*E"                                        *
*----------------------------------------------------------------------*
      call star(NDsup,  Rmsbot,  Rmtop,  nmat, nmutot )
      call star(NDsup,  Rm1sbot, Rm1top, nmat, nmutot )
      call rdiapr(NDsup, Rmsbot,  Rmsbot,  Ebot, nmat, nmutot)
      call rdiapr(NDsup, Rm1sbot, Rm1sbot, Ebot, nmat, nmutot)
      call ldiapr(NDsup, Rmsbot,  Ebot, Rmsbot, nmat, nmutot)
      call ldiapr(NDsup, Rm1sbot, Ebot, Rm1sbot, nmat, nmutot)
*----------------------------------------------------------------------*
*  Direct transmission exp(-b/mu) = exp(-bbot/mu)*exp(-btop/mu)        *
*----------------------------------------------------------------------*
      do 300 mu=1, nmutot
          ebbot(mu) = ebbot(mu)*ebtop(mu)
  300 continue
      return
      end
      subroutine notop( NDmu, NDsup, Rmbot, Tmbot, Rm1bot, Tm1bot
     +                , ebtop, ebbot, nmutot, nmat )
*----------------------------------------------------------------------*
*  Use adding equations when there is no scattering in the top layer   *
*  Do this also for first and zero-th order scattering                 *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension Rmbot(NDsup,NDsup),  Tmbot(NDsup,NDsup)
     +        , Rm1bot(NDsup,NDsup), Tm1bot(NDsup,NDsup)
     +        , ebtop(NDmu), ebbot(NDmu), Etop(NDsup), Ebot(NDsup)
      do 200 mu=1, nmutot
          do 100 k=1, nmat
              Etop((mu-1)*nmat+k) = ebtop(mu)
              Ebot((mu-1)*nmat+k) = ebbot(mu)
  100     continue
  200 continue
*----------------------------------------------------------------------*
*  Reflection  R =  E'R"E'                                             *
*----------------------------------------------------------------------*
      call rdiapr(NDsup, Rmbot,  Rmbot,  Etop, nmat, nmutot)
      call rdiapr(NDsup, Rm1bot, Rm1bot, Etop, nmat, nmutot)
      call ldiapr(NDsup, Rmbot,  Etop, Rmbot,  nmat, nmutot)
      call ldiapr(NDsup, Rm1bot, Etop, Rm1bot, nmat, nmutot)
*----------------------------------------------------------------------*
*  Transmission T = T"E'                                               *
*----------------------------------------------------------------------*
      call rdiapr(NDsup, Tmbot,  Tmbot,  Etop, nmat, nmutot)
      call rdiapr(NDsup, Tm1bot, Tm1bot, Etop, nmat, nmutot)
*----------------------------------------------------------------------*
*  Reflection star R* = R"* is trivial !                               *
*----------------------------------------------------------------------*
*  Direct transmission exp(-b/mu) = exp(-bbot/mu)*exp(-btop/mu)        *
*----------------------------------------------------------------------*
      do 300 mu=1, nmutot
          ebbot(mu) = ebbot(mu)*ebtop(mu)
  300 continue
      return
      end
      subroutine oneadd( NDmu, NDsup, Rm1top, Tm1top, Rm1bot, Rm1sbot
     +                 , Tm1bot, nmat, nmutot 
     +                 , ebtop, ebbot )
*----------------------------------------------------------------------*
*  Use 'adding' equations for zero and first order only.               *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension Rm1top(NDsup,NDsup), Tm1top(NDsup,NDsup)
     +        , Rm1bot(NDsup,NDsup), Rm1sbot(NDsup,NDsup)
     +        , Tm1bot(NDsup,NDsup)
     +        , ebtop(NDmu), ebbot(NDmu), Etop(NDsup), Ebot(NDsup)
      do 200 mu=1, nmutot
          do 100 k=1, nmat
              Etop((mu-1)*nmat+k) = ebtop(mu)
              Ebot((mu-1)*nmat+k) = ebbot(mu)
  100     continue
  200 continue
*----------------------------------------------------------------------*
*  Reflection  R = R' + E'R"E'                                         *
*----------------------------------------------------------------------*
      call rdiapr(NDsup, Rm1bot, Rm1bot, Etop, nmat, nmutot)
      call ldiapr(NDsup, Rm1bot, Etop, Rm1bot, nmat, nmutot)
      call addSM( NDsup,  Rm1bot, Rm1bot, Rm1top, nmat, nmutot)
*----------------------------------------------------------------------*
*  Transmission T = E"T' + T"E'                                        *
*----------------------------------------------------------------------*
      call ldiapr(NDsup, Tm1top, Ebot, Tm1top, nmat, nmutot)
      call rdiapr(NDsup, Tm1bot, Tm1bot, Etop, nmat, nmutot)
      call addSM( NDsup,  Tm1bot, Tm1bot, Tm1top, nmat, nmutot)
*----------------------------------------------------------------------*
*  Reflection star R* = R"* + E"R'*E"                                  *
*----------------------------------------------------------------------*
      call star(NDsup,  Rm1top, Rm1top, nmat, nmutot )
      call rdiapr(NDsup, Rm1top, Rm1top, Ebot, nmat, nmutot)
      call ldiapr(NDsup, Rm1top, Ebot, Rm1top, nmat, nmutot)
      call addSM( NDsup,  Rm1sbot, Rm1sbot, Rm1top, nmat, nmutot)
*----------------------------------------------------------------------*
*  Direct transmission exp(-b/mu) = exp(-bbot/mu)*exp(-btop/mu)        *
*----------------------------------------------------------------------*
      do 300 mu=1, nmutot
          ebbot(mu) = ebbot(mu)*ebtop(mu)
  300 continue
      return
      end
      subroutine ord1m(NDmu, NDsup, m, xmu, smf, nmutot, nmat
     +                , Zmplus, Zmmin, a, b, ebmu, Rm1, Tm1 )
*----------------------------------------------------------------------*
*  Calculate the first order scattering contribution to the m-th       *
*  Fourier component of reflection and transmission of a homogeneous   *
*  layer. Formulae can be found in de Haan et al. (1987)               *
*  Eqs. (A-21)-(A-23)                                                  *
*  The resulting supermatrices are returned through the Rm1 and Tm1.   *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), smf(NDmu), ebmu(NDmu)
     +        , Zmmin(NDsup,NDsup), Zmplus(NDsup,NDsup)
     +        , Rm1(NDsup,NDsup),   Tm1(NDsup,NDsup)
      quarta = 0.25D0*a
      do 300 i=1, nmutot
          xmui = xmu(i)
          ei   = ebmu(i)
          awi  = quarta*smf(i)
          im   = (i-1)*nmat
          do 200 j=1, nmutot
              xmuj= xmu(j)
              ej  = ebmu(j)
              awij= awi*smf(j)
              jm  = (j-1)*nmat
              eiej= 1.D0-ei*ej
              if (eiej .gt. 1.D-3) then
                  h = xmui+xmuj
                  if (h .gt. 1.D-10) h = 1.D0/h
                  hR = awij*h*eiej
                  if (dabs(xmui-xmuj) .gt. 1.D-10) then
                      hT = awij*(ei-ej)/(xmui-xmuj)
                  else
                      h = 0.D0
                      if (xmui .gt. 1.D-10) h = b*ei/(xmui*xmui)
                      hT = awij*h
                  endif
              else
*----------------------------------------------------------------------*
*  use Taylor series to avoid loss of accuarcy when b << 1             *
*----------------------------------------------------------------------*
                  bperi  = b/xmui
                  bperj  = b/xmuj
                  bperij = bperi+bperj
                  h      = 1.D0-0.5D0*bperij*(1.D0-bperij/3.D0)
                  y      = awij*bperi/xmuj
                  hR     = y*h
                  hT     = y*(h-bperi*bperj/6.D0)
              endif
              do 100 k=1, nmat
                  ik = im+k
                  do 50 l=1, nmat
                      jl = jm+l
                      Rm1(ik,jl) = hR*Zmmin(ik,jl)
                      Tm1(ik,jl) = hT*Zmplus(ik,jl)
   50             continue
  100         continue
  200     continue
  300 continue
      return
      end
      subroutine ord2m( NDmu, NDsup, m, xmu, smf, nmug, nmutot, nmat
     +                , Zmplus, Zmmin, a, b, ebmu, Rm, Tm )
*----------------------------------------------------------------------*
*  Calculate second order scattering contribution to the m-th Fourier  *
*  component of reflection and transmission of a homogeneous layer.    *
*  Formulae can be found in                                            *
*      Hovenier, (1971) Astron. Astrophys. 13, p. 7-29                 *
*      Eqs. (A-26)-(A-38)                                              *
*  The second order scattering is added to whatever is in Rm and Tm.   *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), smf(NDmu), ebmu(NDmu)
     +        , Zmmin(NDsup,NDsup), Zmplus(NDsup,NDsup)
     +        , Rm(NDsup,NDsup),    Tm(NDsup,NDsup)
      dimension yik(4,4)
      aa16=a*a/16.D0
      do 200 kl=1,nmat
          do 100 il=1,nmat
              yik(il,kl) = 1.D0
              if ((kl.ge.3) .and. (il.le.2)) yik(il,kl) = -1.D0
              if ((il.ge.3) .and. (kl.le.2)) yik(il,kl) = -1.D0
  100     continue
  200 continue
      do 1100 i = 1,nmutot
*
      xi = xmu(i)
      ei = ebmu(i)
      asmfi = aa16*smf(i)
      im = (i-1)*nmat
      iz = 0
      if (xi .gt. 1.D-10) iz = 1
      bpi = 0.D0
      if (iz .eq .1D0) bpi = b/xi
      beigxi = bpi*ei
      do 1000 j = 1,i
*
      xj = xmu(j)
      ej = ebmu(j)
      asmfij = asmfi*smf(j)
      jm = (j-1)*nmat
      jz = 0
      if (xj .gt. 1.D-10) jz = 1
*----------------------------------------------------------------------*
*  When both mu(i) and mu(j) are almost zero : set reflection          *
*  and transmission to zero and go to next mu(j) if there is one.      *
*----------------------------------------------------------------------*
      if ((jz.eq.0) .and. (iz.eq.0)) then
          do 300 k=1, nmat
              ik = im+k
              do 85 l=1, nmat
                  jl = jm+l
                  Rm(ik,jl) = 0.D0
                  Tm(ik,jl) = 0.D0
   85         continue
  300     continue
      endif
*
      bpj = 0.D0
      if (jz .eq. 1) bpj = b/xj
      bij = bpi+bpj
      xipxj = xi+xj
      ximxj = xi-xj
      xjximxj = 0.D0
      if (i .ne. j) xjximxj = xj/ximxj
      eimej = ei-ej
      eiej = 1.D0-ei*ej
*----------------------------------------------------------------------*
*  Taylor series to avoid loss of accuracy when b << 1.                *
*----------------------------------------------------------------------*
      if (eiej .lt. 1.D-3) then
          eiej = bij*(1.D0-bij/2.D0*(1.D0-bij/3.D0))
          eimej = bpj*(1.D0-bpj/2.D0*(1.D0-bpj/3.D0))-
     +          bpi*(1.D0-bpi/2.D0*(1.D0-bpi/3.D0))
      endif
      e1 = xjximxj*eimej
      g1 = xj/xipxj*eiej
*----------------------------------------------------------------------*
*  Start integration over mu'                                          *
*----------------------------------------------------------------------*
      do 900 k=1, nmug
*
      xk = xmu(k)
      ek = ebmu(k)
      km = (k-1)*nmat
      bpk = b/xk
      bik = bpi+bpk
      bjk = bpj+bpk
*----------------------------------------------------------------------*
*  Calculate the functions e, f, g and h in Hovenier (1971)            *
*----------------------------------------------------------------------*
      if ((bpi.lt.1.D-3) .and. (bpj.lt.1.D-3) .and. (bpk.lt.1.D-3)) then
*----------------------------------------------------------------------*
*         Use Taylor series to avoid loss of accuracy when b << 1.     *
*----------------------------------------------------------------------*
          if (iz .eq. 0) then
              z = b/xk/xj
              e = 0.D0
              g = z*(1.D0-bjk/2.D0*(1.D0-bjk/3.D0))
              f = g-z*z*b/6.D0
              h = 0.D0
          else if (jz .eq. 0) then
              z = b/xk/xi
              e = 0.D0
              h = z*(1.D0-bik/2.D0*(1.D0-bik/3.D0))
              f = h-z*z*b/6.D0
              g = 0.D0
          else
              z = bpi*bpj*bpk/b/2.D0
              e = z*(1.D0-(bpk+2.D0*(bpi+bpj))/3.D0)
              f = z*(1.D0-(bpk+bpi+bpj)/3.D0)
              g = z*(1.D0-(bpk+bpi+2.D0*bpj)/3.D0)
              h = z*(1.D0-(bpk+bpj+2.D0*bpi)/3.D0)
          endif
      else
*----------------------------------------------------------------------*
*         No Taylor series is needed.                                  *
*----------------------------------------------------------------------*
          xipxk = xi+xk
          xjpxk = xj+xk
          ximxk = xi-xk
          xjmxk = xj-xk
          eimek = ei-ek
          ejmek = ej-ek
          eiek = 1.D0-ei*ek
          ejek = 1.D0-ej*ek
*----------------------------------------------------------------------*
*         Use Taylor series to avoid loss of accuracy when b << 1.     *
*----------------------------------------------------------------------*
          if (eiek .lt. 1.D-3) then
              eiek = bik*(1.D0-bik/2.D0*(1.D0-bik/3.D0))
              eimek = bpk*(1.D0-bpk/2.D0*(1.D0-bpk/3.D0))-
     +                       bpi*(1.D0-bpi/2.D0*(1.D0-bpi/3.D0))
          endif
          if (ejek .lt. 1.D-3) then
              ejek = bjk*(1.D0-bjk/2.D0*(1.D0-bjk/3.D0))
              ejmek = bpk*(1.D0-bpk/2.D0*(1.D0-bpk/3.D0))-
     +                       bpj*(1.D0-bpj/2.D0*(1.D0-bpj/3.D0))
          endif
          if (i .eq. j) then
              if (i .ne. k) then
                  e = (b/xj*ej-xk/xipxk*ej*ejek)/xjpxk
                  f = (b/xj*ej-xk/xjmxk*ejmek)/xjmxk
                  g = (g1-xk/ximxk*ej*eimek)/xjpxk
                  h = (g1-xk/xipxk*eiek)/xjmxk
              else
                  e = (b/xj*ej-xk/xipxk*ej*ejek)/xjpxk
                  f = b*b/2.D0/xj/xj/xj*ej
                  g = (g1-b/xi*ei*ej)/xipxj
                  h = g
              endif
          else if (i .eq .k) then
              e = (e1-xk/xipxk*ej*eiek)/xjpxk
              f = (beigxi-e1)/xj*xjximxj
              g = (g1-beigxi*ej)/xipxj
              h = (g1-xk/xipxk*eiek)/xjmxk
          else if (j .eq .k) then
              e = (e1-xk/xipxk*ej*eiek)/xjpxk
              f = (xi/xj*e1-b/xj*ej)/xj*xjximxj
              g = (g1-xk/ximxk*ej*eimek)/xjpxk
              h = (xi/xipxj*eiej-b/xj*ei*ej)/xipxj
          else
              e = (e1-xk/xipxk*ej*eiek)/xjpxk
              f = (e1-xk/ximxk*eimek)/xjmxk
              g = (g1-xk/ximxk*ej*eimek)/xjpxk
              h = (g1-xk/xipxk*eiek)/xjmxk
          endif
      endif
*----------------------------------------------------------------------*
*  End of calculation of functions e, f, g and h                       *
*----------------------------------------------------------------------*
      y = smf(k)*smf(k)*asmfij/xk
      do 800 il=1, nmat
         iil = im+il
         do 700 jl=1, nmat
            jjl = jm+jl
            do 600 kl=1, nmat
               kkl = km+kl
               Zpnik = Zmplus(iil,kkl)
               Zpnkj = Zmplus(kkl,jjl)
               Zmnik = Zmmin(iil,kkl)
               Zmnkj = Zmmin(kkl,jjl)
               Zmsik = yik(il,kl)*Zmnik
               Zpsik = yik(il,kl)*Zpnik
               Rm(iil,jjl) = Rm(iil,jjl)+y*(Zpsik*Zmnkj*g+Zmnik*Zpnkj*h)
               Tm(iil,jjl) = Tm(iil,jjl)+y*(Zmsik*Zmnkj*e+Zpnik*Zpnkj*f)
  600       continue
  700    continue
  800 continue
  900 continue
*----------------------------------------------------------------------*
*  End of integration over mu'                                         *
*----------------------------------------------------------------------*
 1000 continue
*----------------------------------------------------------------------*
*  End of loop over column index j                                     *
*----------------------------------------------------------------------*
 1100 continue
*----------------------------------------------------------------------*
*  End of loop over row index i                                        *
*  Fill upper triangle of the supermatrices using symmetry relations.  *
*----------------------------------------------------------------------*
      call fillup( NDsup, Rm, Tm, nmat, nmutot )
      return
      end
      subroutine prcoef(NDlay, NDcoef, nlayer, nmat, coefs, ncoefs )
*----------------------------------------------------------------------*
*  Print the coefficients for the expansion of the scattering matrix   *
*  in generalized spherical functions.                                 *
*  On entry, the coefficients should be supplied through array coefs.  *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension coefs(4,4,0:NDcoef,NDlay), ncoefs(NDlay)
*----------------------------------------------------------------------*
*  Loop over the layers :                                              *
*----------------------------------------------------------------------*
      print *,'1PRCOEF: print the coefficients of the scattering matrix'
      print *,'         for each layer in the atmosphere.'
      do 400 layer=1, nlayer
          write(*,401) layer
          if (nmat .eq. 1) then
              write(*,402)
              do 100 i=0, ncoefs(layer)
                  write(*,403) i, coefs(1,1,i,layer)
  100         continue
          else if (nmat .eq. 3) then
              write(*,404)
              do 200 i=0, ncoefs(layer)
                  write(*,405) i, coefs(1,1,i,layer), coefs(2,2,i,layer)
     +                          , coefs(3,3,i,layer), coefs(1,2,i,layer)
  200         continue
          else if (nmat .eq. 4) then
              write(*,406)
              do 300 i=0, ncoefs(layer)
                  write(*,407) i, coefs(1,1,i,layer), coefs(2,2,i,layer)
     +                          , coefs(3,3,i,layer), coefs(4,4,i,layer)
     +                          , coefs(1,2,i,layer), coefs(3,4,i,layer)
  300         continue
          else
              print *,' prcoef: ERROR illegal value nmat = ',nmat
              stop 'prcoef: illegal value of nmat'
          endif
  400 continue
*----------------------------------------------------------------------*
*  End of loop over layers.                                            *
*----------------------------------------------------------------------*
  401 format(' The expansion coefficients for layer ',i3,' are :')
  402 format(/,'   l',t14,'alpha1',/,' ---------------------------')
  403 format(i4,e18.11)
  404 format(/,'   l',t14,'alpha1',t32,'alpha2',t50,'alpha3'
     +      ,t86,'beta1',/
     +      ,' -------------------------------------------------------'
     +      ,'------------------')
  405 format(i4,4e18.11)
  406 format(/,'   l',t14,'alpha1',t32,'alpha2',t50,'alpha3'
     +      ,t68,'alpha4',t86,'beta1',t104,'beta2',/
     +      ,' -------------------------------------------------------'
     +      ,'--------------------------------------------------------')
  407 format(i4,6e18.11)
      return
      end
c$$$      subroutine prsup( NDsup, A, nmutot, nmat, iopt )
c$$$*----------------------------------------------------------------------*
c$$$*  Tool to print (part of) a supermatrix for debugging.                *
c$$$*      iopt = 1 : print only diagonal                                  *
c$$$*      iopt = 2 : print entire matrix (scalar only)                    *
c$$$*      iopt = 3 : print corner (nmat by nmat) submatrices              *
c$$$*      iopt = 4 : print (nmat by nmat) submatrices for last two mu     *
c$$$*----------------------------------------------------------------------*
c$$$      implicit double precision (a-h,o-z)
c$$$      parameter ( pi=3.1415926535897932384D0 )
c$$$      parameter ( twopi=2.D0*pi )
c$$$      parameter ( radfac=pi/180.D0 )
c$$$      dimension A(NDsup,NDsup)
c$$$      if (iopt .eq. 1) then
c$$$*----------------------------------------------------------------------*
c$$$*         option 1                                                     *
c$$$*----------------------------------------------------------------------*
c$$$          print *,' prsup: option 1 print diagonal'
c$$$          do 200 mu=1, nmutot
c$$$              i = (mu-1)*nmat
c$$$              print *,mu,' ',A(i,i)
c$$$  200     continue
c$$$      else if (iopt .eq. 2) then
c$$$*----------------------------------------------------------------------*
c$$$*         option 2                                                     *
c$$$*----------------------------------------------------------------------*
c$$$          print *,' prsup: option 2 print whole matrix (scalar)'
c$$$          if (nmat .ne. 1) then
c$$$              print *,' prsup: nmat = ',nmat,' INCOMPATIBLE with opt 2'
c$$$              return
c$$$          endif
c$$$          write(*,301) (i,i=1,nmutot)
c$$$          do 300 i=1, nmutot
c$$$              write(*,302) i,(A(i,j),j=1,nmutot)
c$$$  300     continue
c$$$  301     format(' ',100(i9))
c$$$  302     format(i3,100(1pe9.1))
c$$$      else if (iopt .eq. 3) then
c$$$*----------------------------------------------------------------------*
c$$$*         option 3                                                     *
c$$$*----------------------------------------------------------------------*
c$$$          print *,' prsup: option 3 print corner submatrices'
c$$$          write(*,401) (1,i=1,nmat),(nmutot,i=1,nmat)
c$$$          do 400 mu=1, nmutot,nmutot-1
c$$$              ibase = (mu-1)*nmat
c$$$              jbase = (nmutot-1)*nmat
c$$$              do 390 k=1, nmat
c$$$                  i = ibase+k
c$$$                  write(*,402) mu,(A(i,j),j=1,nmat)
c$$$     +                           ,(A(i,j),j=jbase+1,jbase+nmat)
c$$$  390         continue
c$$$              write(*,403)
c$$$  400     continue
c$$$  401     format(100(i9))
c$$$  402     format(i3,100(1pe9.1))
c$$$  403     format('  -------------------------------------'
c$$$     +          ,'  ----------------------------------')
c$$$      else if (iopt .eq. 4) then
c$$$*----------------------------------------------------------------------*
c$$$*         option 4                                                     *
c$$$*----------------------------------------------------------------------*
c$$$          print *,' prsup: option 4 print last two mu point submatrices'
c$$$          write(*,501) (nmutot-1,i=1,nmat),(nmutot,i=1,nmat)
c$$$          do 500 mu=nmutot-1, nmutot
c$$$              ibase = (mu-1)*nmat
c$$$              jbase = (nmutot-2)*nmat
c$$$              do 590 k=1, nmat
c$$$                  i = ibase+k
c$$$                  write(*,502) mu,(A(i,j),j=jbase+1,jbase+nmat)
c$$$     +                           ,(A(i,j),j=jbase+1+nmat,jbase+2*nmat)
c$$$  590         continue
c$$$              write(*,503)
c$$$  500     continue
c$$$  501     format(100(i9))
c$$$  502     format(i3,100(1pe9.1))
c$$$  503     format('  -------------------------------------'
c$$$     +          ,'  ----------------------------------')
c$$$      else
c$$$          print *,' prsup: illegal option ',iopt
c$$$          stop 'in prsup illegal option'
c$$$      endif
c$$$      return
c$$$      end
      subroutine rRprod( NDmu,NDsup,  A, B, R, xmu, nmat, nmutot )
************************************************************************
*
*
*
*
*
*
*
************************************************************************
      implicit double precision(a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup,NDsup), R(NDmu,3)
     +        , xmu(NDmu)

      if (nmat .eq. 1) then
          do 200 mu0=1, nmutot
              j = (mu0-1)*nmat + 1
              x2mu0 = xmu(mu0) + xmu(mu0)
              do 100 mu=1, nmutot
                  i = (mu-1)*nmat + 1
                  A(i,j) = x2mu0*( B(i,j)*R(mu0,1) )
  100         continue
  200     continue

      else if (nmat.eq.3) then
        do 500 mu0=1, nmutot
          j = (mu0-1)*nmat
          x2mu0 = xmu(mu0) + xmu(mu0)
          do 400 mu=1, nmutot
            i = (mu-1)*nmat
            do 300 k=1, nmat
              ki = k+i
              A(ki,j+1) = x2mu0*(B(ki,j+1)*R(mu0,1)+B(ki,j+2)*R(mu0,2))
              A(ki,j+2) = x2mu0*(B(ki,j+1)*R(mu0,2)+B(ki,j+2)*R(mu0,1))
              A(ki,j+3) = x2mu0*(B(ki,j+3)*R(mu0,3))
  300       continue
  400     continue
  500   continue

      else
        do 800 mu0=1, nmutot
          j = (mu0-1)*nmat
          x2mu0 = xmu(mu0) + xmu(mu0)
          do 700 mu=1, nmutot
            i = (mu-1)*nmat
            do 600 k=1, nmat
              ki = k+i
              A(ki,j+1) = x2mu0*(B(ki,j+1)*R(mu0,1)+B(ki,j+2)*R(mu0,2))
              A(ki,j+2) = x2mu0*(B(ki,j+1)*R(mu0,2)+B(ki,j+2)*R(mu0,1))
              A(ki,j+3) = x2mu0*(B(ki,j+3)*R(mu0,3))
              A(ki,j+4) = x2mu0*(B(ki,j+4)*R(mu0,3))
  600       continue
  700     continue
  800   continue
      endif
      return
      end
      subroutine rTprod(NDmu,NDsup,  A, Tfs, B, xmu, xm, nmat, nmutot )
*----------------------------------------------------------------------*
*  Calculate the following 'right transmission product' :              *
*                                                                      *
*    A(mu,mu0) = B(mu,mu0) 2 mu0 TF*(xmut(mu0)) / xm**2                *
*                                                                      *
*  where xmut(mu) is defined in Chowdhary (1991) Eq. (7), and TF*(mu)  *
*  is the Fresnel matrix defined in Chowdhary (1991) Eq. (27).         *
*  xm is the refractive index of water with respect to air.            *
*                                                                      *
*  On entry:                                                           *
*     Tfs    : Fresnel matrix given by Eq. (27), coded as follows      *
*                  1,1 element of TF*(xmut(mu0)) in Tfs(mu0,1)         *
*                  1,2 element of TF*(xmut(mu0)) in Tfs(mu0,2)         *
*                  3,3 element of TF*(xmut(mu0)) in Tfs(mu0,3)         *
*     B      : supermatrix B to be multiplied                          *
*     xmu    : all different mu values (integration and extra points)  *
*     xm     : refractive index of water with respect to air           *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     nmutot  : total number of distinct mu points                     *
*                                                                      *
*  On exit:                                                            *
*     A      : supermatrix given by above formula                      *
*                                                                      *
*  Note that A and B may be the same when this subroutine is called!   *
*  This is the reason for the awkward use of variable A1 (see below).  *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup, NDsup), Tfs(NDmu,3), xmu(NDmu)
      if (nmat .eq. 1) then
*----------------------------------------------------------------------*
*  Treat case without polarization                                     *
*----------------------------------------------------------------------*
          do 200 mu0=1, nmutot
              xmT1 = 2.0d0*xmu(mu0)*Tfs(mu0,1)/xm**2
              do 100 mu=1, nmutot
                  A(mu,mu0) = B(mu,mu0)*xmT1
  100         continue
  200     continue
      else if (nmat .eq. 3) then
*----------------------------------------------------------------------*
*  Treat case with 3 X 3 approximation                                 *
*----------------------------------------------------------------------*
        do 500 mu0=1, nmutot
          jbase = (mu0-1)*nmat
          xmT1 = 2.0d0*xmu(mu0)*Tfs(mu0,1)/xm**2
          xmT2 = 2.0d0*xmu(mu0)*Tfs(mu0,2)/xm**2
          xmT3 = 2.0d0*xmu(mu0)*Tfs(mu0,3)/xm**2
          do 400 mu=1, nmutot
            ibase = (mu-1)*nmat
            do 300 k=1, nmat
              ik = ibase + k
              A1            = B(ik,jbase+1)*xmT1 + B(ik,jbase+2)*xmT2
              A(ik,jbase+2) = B(ik,jbase+1)*xmT2 + B(ik,jbase+2)*xmT1
              A(ik,jbase+3) = B(ik,jbase+3)*xmT3
              A(ik,jbase+1) = A1
  300       continue
  400     continue
  500   continue
      else
*----------------------------------------------------------------------*
*  Treat case with polarization fully included                         *
*----------------------------------------------------------------------*
        do 800 mu0=1, nmutot
          jbase = (mu0-1)*nmat
          xmT1 = 2.0d0*xmu(mu0)*Tfs(mu0,1)/xm**2
          xmT2 = 2.0d0*xmu(mu0)*Tfs(mu0,2)/xm**2
          xmT3 = 2.0d0*xmu(mu0)*Tfs(mu0,3)/xm**2
          do 700 mu=1, nmutot
            ibase = (mu-1)*nmat
            do 600 k=1, nmat
              ik = ibase + k
              A1            = B(ik,jbase+1)*xmT1 + B(ik,jbase+2)*xmT2
              A(ik,jbase+2) = B(ik,jbase+1)*xmT2 + B(ik,jbase+2)*xmT1
              A(ik,jbase+3) = B(ik,jbase+3)*xmT3
              A(ik,jbase+4) = B(ik,jbase+4)*xmT3
              A(ik,jbase+1) = A1
  600       continue
  700     continue
  800   continue
      endif
      return
      end
      subroutine renorm( NDmu, NDsup
     +                 , Zmmin, Zmplus, nmug, nmutot, nmat, xmu, smf
     +                 , epsilon )
*----------------------------------------------------------------------*
*  Renormalize the phase matrix by updating the diagonal elements      *
*  of the phase matrix.                                                *
*  This routine only makes sense if called for m=0 Fourier component.  *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension Zmmin(NDsup,NDsup), Zmplus(NDsup,NDsup)
     +        , xmu(NDmu), smf(NDmu)
      dimension w(NDmu)
      logical verbo
      verbo = .false.
*----------------------------------------------------------------------*
*  Retrieve the weights from the supermatrix factors smf :             *
*----------------------------------------------------------------------*
      do 250 i=1, nmug
          w(i) = 0.5D0*smf(i)**2/xmu(i)
  250 continue
      fmax  =-1.0d0
      fmaxex=-1.0d0
      do 400 j=1, nmutot
*----------------------------------------------------------------------*
*         Calculate normalization integral r+t, which should be 2      *
*----------------------------------------------------------------------*
          r = 0.D0
          t = 0.D0
          jsup = (j-1)*nmat+1
          do 300 i=1, nmug
              isup = (i-1)*nmat+1
              r = r + Zmmin(isup,jsup)*w(i)
              t = t + Zmplus(isup,jsup)*w(i)
  300     continue
*----------------------------------------------------------------------*
*         Update digonal elements of phase matrix Fourier component    *
*         Only change the transmission part Zmplus, not Zmmin          *
*         For the extra mu-points the two points closest to the        *
*         forward direction are updated, each weighted by interpolation*
*----------------------------------------------------------------------*
          if (j .le. nmug) then
              fac = 1.0d0 + (2.0d0-r-t)/(Zmplus(jsup,jsup)*w(j))
              if (fac .gt. fmax) fmax = fac
              do 340 k=1, nmat
               Zmplus(jsup+k-1,jsup+k-1) = fac*Zmplus(jsup+k-1,jsup+k-1)
  340         continue
          else
            call brack( xmu(j), xmu, nmug, NDmu, i1, i2 )
            if (i1 .ne. i2) then
                relw1 = (xmu(i2)-xmu(j))/(xmu(i2)-xmu(i1))
                relw2 = (xmu(j)-xmu(i1))/(xmu(i2)-xmu(i1))
            else
                relw1 = 1.0d0
                relw2 = 0.0d0
            endif
            isup1 = (i1-1)*nmat+1
            isup2 = (i2-1)*nmat+1
            fac1 = 1.0d0 + relw1*(2.0d0-r-t)/(Zmplus(isup1,jsup)*w(i1))
            fac2 = 1.0d0 + relw2*(2.0d0-r-t)/(Zmplus(isup2,jsup)*w(i2))
            if (fac1 .gt. fmaxex) fmaxex = fac1
            do 440 k=1, nmat
              Zmplus(isup1+k-1,jsup+k-1)=fac1*Zmplus(isup1+k-1,jsup+k-1)
              Zmplus(jsup+k-1,isup1+k-1)=Zmplus(isup1+k-1,jsup+k-1)
              Zmplus(isup2+k-1,jsup+k-1)=fac2*Zmplus(isup2+k-1,jsup+k-1)
              Zmplus(jsup+k-1,isup2+k-1)=Zmplus(isup2+k-1,jsup+k-1)
  440       continue
          endif
  400 continue
*----------------------------------------------------------------------*
*  Print maximum renormalization factor (verbose only)                 *
*----------------------------------------------------------------------*
      if (verbo) then
        print *,' renorm: max renormalization factor ',fmax,' Gauss pts'
        print *,'         ',fmaxex,' (extra points)'
      endif
      return
      end
      subroutine scalZm(NDmu, NDsup, NDlay, NDcoef
     +                 , m, layer, coefs, ncoefs, xmu, smf, epsilon
     +                 , nmug, nmutot, nmat, Zmmin, Zmplus )
*----------------------------------------------------------------------*
*  Calculate the m-th Fourier component of the phase matrix Zm(mu0,mu) *
*  from the expansion coefficients of the scattering matrix into       *
*  generalized spherical functions. The formulae can be found in :     *
*                                                                      *
*    J.F. de Haan et al.: 1987, Astron. Astrophys. 183, pp. 371-391.   *
*                                                                      *
*  Essentially, Eqs. (66)-(82) are used.                               *
*                                                                      *
*  NO POLARIZATION IS INCLUDED !!                                      *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension Plm(NDmu,3,2), sqlm(0:NDcoef), smf(NDmu)
     +        , coefs(4,4,0:NDcoef,NDlay), ncoefs(NDlay), xmu(NDmu)
     +        , Zmmin(NDsup,NDsup), Zmplus(NDsup,NDsup)
      logical odd
*----------------------------------------------------------------------*
*  Precompute the factor dsqrt(l**2-m**2) needed in Eqs. (81)-(82)     *
*----------------------------------------------------------------------*
      m2 = m*m
      sqlm(m) = 0.D0
      do 100 l=m+1, ncoefs(layer)
          sqlm(l) = dsqrt(dble(l)**2-dble(m2))
  100 continue
*----------------------------------------------------------------------*
*  Compute the binomial factor 2m above m needed in Eq. (77)           *
*----------------------------------------------------------------------*
      bin2mm = 1.D0
      do 200 n=1, m
          bin2mm = bin2mm*dble(n+m)/dble(n)
  200 continue
      binfac = 2.D0**(-m)*dsqrt(bin2mm)
*----------------------------------------------------------------------*
*     Initialize Plm for l=m, Eq.(77) without factor i**-m             *
*     this factor drops out anyway in Eq. (66)                         *
*----------------------------------------------------------------------*
      lold = 1
      lnew = 2
      do 300 i=1, nmutot
          u  = xmu(i)
          rootu = dsqrt(dabs(1.D0-u*u))
          Plm(i,1,lold) = 0.D0
          if (m .ne. 0) then
              Plm(i,1,lnew) = binfac*rootu**m
          else
              Plm(i,1,lnew) = binfac
          endif
  300 continue
*----------------------------------------------------------------------*
*         Initialize phase matrix with term l=m, Eq.(66)               *
*----------------------------------------------------------------------*
      do 500 i=1,nmutot
          SP = coefs(1,1,m,layer)*Plm(i,1,lnew)
          do 400 j=1,nmutot
              Zmplus(i,j) = SP*Plm(j,1,lnew)
              Zmmin(i,j)  = Zmplus(i,j)
  400     continue
  500 continue
      if (ncoefs(layer) .gt. m) then
*----------------------------------------------------------------------*
*         Start loop over l (summation index in Eq. (66))              *
*----------------------------------------------------------------------*
          odd = .false.
          do 1200 l=m+1, ncoefs(layer)
              odd = .not. odd
*----------------------------------------------------------------------*
*             Do one step in recurrence for Plm, Eq.(81)               *
*----------------------------------------------------------------------*
              c1 = dble(l+l-1)/sqlm(l)
              c2 = sqlm(l-1)/sqlm(l)
              do 700 i=1, nmutot
                  u  = xmu(i)
                  Plm(i,1,lold) =c1*u*Plm(i,1,lnew)-c2*Plm(i,1,lold)
  700         continue
              itmp = lnew
              lnew = lold
              lold = itmp
*----------------------------------------------------------------------*
*             Add a new term to Zm, Eq.(66)                            *
*----------------------------------------------------------------------*
              if (odd) then
                  do 900 i=1,nmutot
                      SP = coefs(1,1,l,layer)*Plm(i,1,lnew)
                      do 800 j=1,nmutot
* OPTION: REWRITE IN CASE THIS LOOP DOES NOT VECTORIZE
*                         term = SP*Plm(j,1,lnew)
*                          Zmplus(i,j) = Zmplus(i,j)+ term
*                          Zmmin(i,j)  = Zmmin(i,j) - term
* END OPTION
                          Zmplus(i,j)= Zmplus(i,j)+ SP*Plm(j,1,lnew)
                          Zmmin(i,j) = Zmmin(i,j) - SP*Plm(j,1,lnew)
  800                 continue
  900             continue
              else
                  do 1100 i=1,nmutot
                      SP = coefs(1,1,l,layer)*Plm(i,1,lnew)
                      do 1000 j=1,nmutot
* OPTION: REWRITE IN CASE THIS LOOP DOES NOT VECTORIZE
*                          term = SP*Plm(j,1,lnew)
*                          Zmplus(i,j) = Zmplus(i,j)+ term
*                          Zmmin(i,j)  = Zmmin(i,j) + term
* END OPTION
                          Zmplus(i,j)= Zmplus(i,j)+ SP*Plm(j,1,lnew)
                          Zmmin(i,j) = Zmmin(i,j) + SP*Plm(j,1,lnew)
 1000                 continue
 1100             continue
              endif
 1200     continue
*----------------------------------------------------------------------*
*             End of summation loop over l                             *
*----------------------------------------------------------------------*
      endif
      return
      end
      subroutine setZm(NDmu, NDsup, NDlay, NDcoef
     +                , m, layer, coefs, ncoefs, xmu, smf, epsilon
     +                , nmug, nmutot, nmat, Zmmin, Zmplus )
*----------------------------------------------------------------------*
*  Calculate the m-th Fourier component of the phase matrix Zm(mu0,mu) *
*  from the expansion coefficients of the scattering matrix into       *
*  generalized spherical functions. The formulae can be found in :     *
*                                                                      *
*    J.F. de Haan et al.: 1987, Astron. Astrophys. 183, pp. 371-391.   *
*                                                                      *
*  Essentially, Eqs. (66)-(82) are used. The suggestion below Eq. (82) *
*  to diagonalize the matrix Plm is followed here.                     *
*  To this end we define the two matrices :                            *
*                                                                      *
*        ( 1     0     0     0 )             ( 1     0     0     0 )   *
*   D1 = ( 0     1     1     0 )        D2 = ( 0    0.5   0.5    0 )   *
*        ( 0     1    -1     0 )             ( 0    0.5  -0.5    0 )   *
*        ( 0     0     0     1 )             ( 0     0     0     1 )   *
*                                                                      *
*  It is clear that D1 and D2 are each other's inverse.                *
*  We diagonalize the Plm given in Eq. (74) :                          *
*                                                                      *
*           D1*Plm*D2  =  diag( Plm0,  Plm-2,  Plm+2,  Plm0 )          *
*                                                                      *
*  The sum in Eq. (66) is written :                                    *
*                                                                      *
*  sum( Plm*Sl*Plm ) = D1* sum( D2*Plm*D1 * D2*Sl*D2 * D1*Plm*D2 ) *D1 *
*                                                                      *
*  The matrix D2*Sl*D2 for a given l is given by :                     *
*                                                                      *
*     (  alpha1        beta1/2             beta1/ 2          0     )   *
*     (  beta1/2   (alpha2+alpha3)/4   (alpha2-alpha3)/4   beta2/2 )   *
*     (  beta1/2   (alpha2-alpha3)/4   (alpha2+alpha3)/4  -beta2/2 )   *
*     (    0          -beta2/2             beta2/2         alpha4  )   *
*                                                                      *
*  where the alpha's and beta's are given by Eqs. (68)-(73)            *
*                                                                      *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension coefs(4,4,0:NDcoef,NDlay), ncoefs(NDlay), xmu(NDmu)
     +        , smf(NDmu), Zmmin(NDsup,NDsup), Zmplus(NDsup,NDsup)
      dimension Plm(NDmu,3,2), DPDpl(NDsup), DPDmi(NDsup)
     +        , DSD(4,4), sqlm(0:NDcoef), sql4(NDcoef), rootu(NDmu)
      logical verbo 
      verbo = .false.

      qroot6 = -0.25D0*dsqrt(6.D0)
      if (nmat .eq. 1) then
          if (verbo) print *,' setZm: use scalar phase matrix'
          call scalZm(NDmu, NDsup, NDlay, NDcoef
     +               ,  m, layer, coefs, ncoefs, xmu, smf, epsilon
     +               , nmug, nmutot, nmat, Zmmin, Zmplus )
          goto 999
      endif
*----------------------------------------------------------------------*
*  Precompute the factor dsqrt(l**2-m**2) needed in Eqs. (81)-(82)     *
*  and also the factor dsqrt(l**2-4) needed in Eq. (82)                *
*----------------------------------------------------------------------*
      m2 = m*m
      do 100 l=m, ncoefs(layer)
          sqlm(l) = dsqrt(dabs(dble(l)**2-dble(m2)))
  100 continue
      do 110 l=2, ncoefs(layer)
          sql4(l) = dsqrt(dabs(dble(l)**2-4.D0))
  110 continue
*----------------------------------------------------------------------*
*  Compute the binomial factor 2m above m needed in Eq. (77)           *
*  and also the binomial factor 2m above m-2 needed in Eq. (80)        *
*----------------------------------------------------------------------*
      bin2mm = 1.D0
      do 200 n=1, m
          bin2mm = bin2mm*dble(n+m)/dble(n)
  200 continue
      twom   = 2.D0**(-m)
      binfac = twom*dsqrt(bin2mm)
      if (m .ge. 2) then
          bin2m2 = bin2mm*dble(m)*dble(m-1)/(dble(m+1)*dble(m+2))
          binf2  = -twom*dsqrt(bin2m2)
      endif
*----------------------------------------------------------------------*
*  Initialize phase matrix to zero                                     *
*----------------------------------------------------------------------*
      nsup = nmat*nmutot
      do 510 j=1,nsup
          do 500 i=1, nsup
              Zmplus(i,j) = 0.D0
              Zmmin(i,j)  = 0.D0
  500     continue
  510 continue
*----------------------------------------------------------------------*
*  Initialize Plm0 for l=m, Eq.(77), and for l=m-1 Eq. (76)            *
*----------------------------------------------------------------------*
      lold = 1
      lnew = 2
      do 300 i=1, nmutot
          rootu(i) = dsqrt(dabs(1.D0-xmu(i)**2))
          if (m .ne. 0) then
              Plm(i,1,lnew) = binfac*rootu(i)**m
          else
              Plm(i,1,lnew) = binfac
          endif
          Plm(i,1,lold) = 0.D0
  300 continue
*----------------------------------------------------------------------*
*  Set Plm2 and Plm-2 to zero: initialization will be done inside loop *
*----------------------------------------------------------------------*
      do 360 i=1, nmutot
          Plm(i,2,lnew) = 0.D0
          Plm(i,2,lold) = 0.D0
          Plm(i,3,lnew) = 0.D0
          Plm(i,3,lold) = 0.D0
  360 continue
*----------------------------------------------------------------------*
*     Start loop over l (summation index in Eq. (66))                  *
*     Parity of Plm is (-1)**(l-m)                                     *
*----------------------------------------------------------------------*
      parity = -1.D0
      do 1200 l=m, ncoefs(layer)
          parity = -parity
*----------------------------------------------------------------------*
*  Initialize Plm for l=max(m,2) Eqs.(78)-(80) without factor i**-m    *
*  This factor is cancelled in Eq. (66) by the factor -1**m            *
*  The exception for m=2 is needed to handle u=+-1 in Eq. (80)         *
*----------------------------------------------------------------------*
          if (l .eq. max0(m,2)) then
              if (m .eq. 0) then
                  do 310 i=1, nmutot
                      Plm(i,2,lnew) = qroot6*rootu(i)*rootu(i)
                      Plm(i,3,lnew) = Plm(i,2,lnew)
  310             continue
              else if (m .eq. 1) then
                  do 320 i=1, nmutot
                      u = xmu(i)
                      Plm(i,2,lnew) = -0.5D0*rootu(i)*(1.D0-u)
                      Plm(i,3,lnew) =  0.5D0*rootu(i)*(1.D0+u)
  320             continue
              else if (m .eq. 2) then
                  do 330 i=1, nmutot
                      u = xmu(i)
                      Plm(i,2,lnew) = -0.25D0*(1.D0-u)**2
                      Plm(i,3,lnew) = -0.25D0*(1.D0+u)**2
  330             continue
              else
                  do 340 i=1, nmutot
                      u = xmu(i)
                      urootm = rootu(i)**(m-2)
                      Plm(i,2,lnew) = binf2*urootm*(1.D0-u)*(1.D0-u)
                      Plm(i,3,lnew) = binf2*urootm*(1.D0+u)*(1.D0+u)
  340             continue
              endif
          endif
*----------------------------------------------------------------------*
*         Construct supervectors corresponding to the diagonal         *
*         elements of the matrix D1 * Plm * D2                         *
*----------------------------------------------------------------------*
          do 610 i=1, nmutot
              isup = nmat*(i-1)
              DPDpl(isup+1) = Plm(i,1,lnew)
              DPDpl(isup+2) = Plm(i,2,lnew)
              DPDpl(isup+3) = Plm(i,3,lnew)
              DPDmi(isup+1) = parity*Plm(i,1,lnew)
              DPDmi(isup+2) = parity*Plm(i,3,lnew)
              DPDmi(isup+3) = parity*Plm(i,2,lnew)
  610     continue
          if (nmat .eq. 4) then
              do 620 i=4, nsup, 4
                  DPDpl(i) = DPDpl(i-3)
                  DPDmi(i) = DPDmi(i-3)
  620         continue
          endif
*----------------------------------------------------------------------*
*         Construct the matrix D2 * S * D2                             *
*----------------------------------------------------------------------*
          DSD(1,1) = coefs(1,1,l,layer)
          DSD(2,1) = 0.5D0*coefs(1,2,l,layer)
          DSD(2,2) = 0.25D0*(coefs(2,2,l,layer)+coefs(3,3,l,layer))
          DSD(3,2) = 0.25D0*(coefs(2,2,l,layer)-coefs(3,3,l,layer))
          DSD(3,1) = DSD(2,1)
          DSD(1,2) = DSD(2,1)
          DSD(1,3) = DSD(2,1)
          DSD(2,3) = DSD(3,2)
          DSD(3,3) = DSD(2,2)
          if (nmat .eq. 4) then
              DSD(1,4) = 0.D0
              DSD(2,4) = 0.5D0*coefs(3,4,l,layer)
              DSD(3,4) = -DSD(2,4)
              DSD(4,4) = coefs(4,4,l,layer)
              DSD(4,1) = 0.D0
              DSD(4,2) = -DSD(2,4)
              DSD(4,3) = -DSD(3,4)
          endif
*----------------------------------------------------------------------*
*         Add a new term to the sum in Eq. (66)                        *
*         The factor (-1)**m is cancelled by i**m in the Plm           *
*----------------------------------------------------------------------*
          do 710 k2=1, nmat
              do 720 k1=1, nmat
                  do 730 j=k2, nsup, nmat
                      SPj = DSD(k1,k2)*DPDpl(j)
                      do 740 i=k1, nsup, nmat
                          Zmplus(i,j) = Zmplus(i,j) + DPDpl(i)*SPj
                          Zmmin(i,j)  = Zmmin(i,j)  + DPDmi(i)*SPj
  740                 continue
  730             continue
  720         continue
  710     continue
*----------------------------------------------------------------------*
*         When last coefficient has been treated : skip recurrence     *
*----------------------------------------------------------------------*
          if (l .eq. ncoefs(layer) ) goto 1200
*----------------------------------------------------------------------*
*         Do one step in recurrence for Plm0, Eq. (81)                 *
*----------------------------------------------------------------------*
          twol1 = 2.D0*l+1.D0
          f1new = twol1/sqlm(l+1)
          f1old = sqlm(l)/sqlm(l+1)
          do 700 i=1, nmutot
              u  = xmu(i)
              Plm(i,1,lold) = f1new*u*Plm(i,1,lnew)- f1old*Plm(i,1,lold)
  700     continue
*----------------------------------------------------------------------*
*         Do one step in recurrence for Plm2 and Plm-2, Eq. (82)       *
*         only when they have been initialized : l >= max(m,2) !!!     *
*----------------------------------------------------------------------*
          if (l .ge. max0(m,2)) then
              tmp   = 1.D0/(dble(l)*sql4(l+1)*sqlm(l+1))
              f2new = twol1*dble(l)*dble(l+1)*tmp
              f2newa= twol1*dble(2*m)*tmp
              f2old = dble(l+1)*sql4(l)*sqlm(l)*tmp
              do 800 i=1, nmutot
                  u = xmu(i)
                  Plm(i,2,lold) = (f2new*u+f2newa)*Plm(i,2,lnew)
     +                                             - f2old*Plm(i,2,lold)
                  Plm(i,3,lold) = (f2new*u-f2newa)*Plm(i,3,lnew)
     +                                             - f2old*Plm(i,3,lold)
  800         continue
          endif
          itmp = lnew
          lnew = lold
          lold = itmp
 1200 continue
*----------------------------------------------------------------------*
*     End of summation loop over l                                     *
*     Calculate D1 * sum * D1                                          *
*----------------------------------------------------------------------*
      call transf(NDsup, Zmmin,  nsup, nmat )
      call transf(NDsup, Zmplus, nsup, nmat )
  999 return
      end
      subroutine setexp(NDmu, NDlay, NDgeom, xmu, imu, imu0, ngeom
     +                 , b, nlayer, e, e0 )
*----------------------------------------------------------------------*
*  Compute the e-powers of the form exp(-b/mu) and exp(-b/mu0) for all *
*  layers and for all geometries.                                      *
*  On entry, the mu-values should be supplied through the array xmu,   *
*  their indexes for different geometries through imu and imu0, and    *
*  the optical thicknesses through bmsca, bmabs and baer and finally   *
*  the number of geometries ngeom and the number of layers nlayer      *
*  should be suplied.                                                  *
*  On exit, the e-powers are returned through arrays e and e0.         *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), imu(NDgeom), imu0(NDgeom), b(NDlay)
      dimension e(NDlay,NDgeom), e0(NDlay,NDgeom)
c      integer nlayer, ngeom
      do 200 layer=1, nlayer
          perb = 1.D0/b(layer)
          do 100 i=1, ngeom
              pmub = xmu(imu(i))*perb
              pmu0b= xmu(imu0(i))*perb
              if (pmub .gt. 1.D-3) then
                  e(layer,i) = dexp(-1.D0/pmub)
              else
                  e(layer,i) = 0.D0
              endif
              if (pmu0b .gt. 1.D-2) then
                  e0(layer,i) = dexp(-1.D0/pmu0b)
              else
                  e0(layer,i) = 0.D0
              endif
  100     continue
  200 continue
      return
      end
      subroutine setfou(NDmu, NDlay, NDcoef, coefs, ncoefs, nlayer, a, b
     +                 , xmu, nmutot, epsilon, M0, M1, M2, iadd )
*----------------------------------------------------------------------*
*  Calculate bounds M0, M1 and M2 on the Fourier index m such that     *
*  for  0  <= m <= M2 doubling must be used,                           *
*  for  M2 <  m <= M1 first plus second order scattering suffice and   *
*  for  M1 <  m <= M0 only first order scattering is needed.           *
*  For  M0 < m        there is no scattering at all !                  *
*  The criteria are described in de Haan et al. (1987) p. 386          *
*  Also determine the options indicating for each layer which type of  *
*  adding procedure should be used :                                   *
*     iadd = 1 : normal adding                                         *
*     iadd = 2 : top layer has no scattering                           *
*     iadd = 3 : bottom layer has no scattering                        *
*  This follows directly from the values of M0 for all layers.         *
*  The number of coefficients is truncated such that lmax <= M0, see   *
*  de Haan et al. (1987) first few lines of section 7.2.               *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension coefs(4,4,0:NDcoef,NDlay), ncoefs(NDlay), xmu(NDmu)
     +        , a(NDlay), b(NDlay)
     +        , M0(NDlay), M1(NDlay), M2(NDlay), iadd(NDlay,0:NDcoef)
      logical verbo
      verbo = .false.
*----------------------------------------------------------------------*
*  Determine minimum mu value to be used in criteria                   *
*----------------------------------------------------------------------*
      xmumin = 1.D0
      do 10 mu=1, nmutot
          if (xmu(mu) .lt. xmumin) xmumin = xmu(mu)
   10 continue
*----------------------------------------------------------------------*
*  Start loop over layers                                              *
*----------------------------------------------------------------------*
      do 1000 layer=1, nlayer
*----------------------------------------------------------------------*
*         Handle case of no scattering in this layer                   *
*----------------------------------------------------------------------*
          if (a(layer) .lt. 1.D-10) then
              M0(layer) = -1
              M1(layer) = -1
              M2(layer) = -1
          else
             fbmu = b(layer)/xmumin
              if (fbmu .gt. 1.D0) fbmu = 1.D0
              f3b = 3.D0*b(layer)
              if (f3b .gt. 1.D0) f3b = 1.D0
              qf   = 0.25D0*fbmu
              qff  = 0.25D0*fbmu*f3b
              qff2 = 0.25D0*fbmu*f3b**2
*----------------------------------------------------------------------*
*  Start loop over Fourier index m to find M0 using Eqs. below (140)   *
*----------------------------------------------------------------------*
              m = ncoefs(layer)
  100             am = a(layer)*coefs(1,1,m,layer)/dble(2*m+1)
                  h  = qf*am*dble(2*m+1)
                  if ((h .le. epsilon) .and. (m .gt. 0)) then
                      m = m-1
                      goto 100
                  endif
              M0(layer) = m
*----------------------------------------------------------------------*
*  Continue loop over Fourier index m to find M1 using Eqs. below (140)*
*----------------------------------------------------------------------*
  200             am = a(layer)*coefs(1,1,m,layer)/dble(2*m+1)
                  h  = qff*am**2*dble(2*m+1)
                  if ((h .le. epsilon) .and. (m .gt. 0)) then
                      m = m-1
                      goto 200
                  endif
              M1(layer) = m
*----------------------------------------------------------------------*
*  Continue loop over Fourier index m to find M2 using Eqs. below (140)*
*----------------------------------------------------------------------*
  300             am = a(layer)*coefs(1,1,m,layer)/dble(2*m+1)
                  h  = qff2*am**3*dble(2*m+1)
                  if ((h .le. epsilon) .and. (m .gt. 0)) then
                      m = m-1
                      goto 300
                  endif
              M2(layer) = m
          endif
          if (verbo) print *,' setfou: layer',layer,'  M0 =',M0(layer)
     +                      ,'  M1 =',M1(layer),'  M2 =',M2(layer)
          ncoefs(layer) = M0(layer)
 1000 continue
*----------------------------------------------------------------------*
*  End of loop over layers                                             *
*  Now determine the adding option iadd for each layer and each        *
*  Fourier index m.                                                    *
*----------------------------------------------------------------------*
      do 3000 m=0, NDcoef
          M0bot = -1
          do 2000 layer=1, nlayer
              if ((m .gt. M0(layer)) .or. (a(layer) .lt. 1.d-10)) then
                  iadd(layer,m) = 2
              else if (m .gt. M0bot) then
                  iadd(layer,m) = 3
              else
                  iadd(layer,m) = 1
              endif
              M0bot = max0(M0bot,M0(layer))
 2000     continue
 3000 continue
      return
      end
      subroutine setgsf(NDmu, NDgeom, NDcoef, lmax, xmu, imu0, imu, phi
     +                       , ngeom, nmat, gsfpl, gsfmi)
*----------------------------------------------------------------------*
*  Calculate the generalized spherical functions (gsf) used in         *
*  Eqs. (68) - (73) of                                                 *
*                                                                      *
*       de Haan et al. ; 1987, Astron. Astrophys. 183  pp. 371-391.    *
*                                                                      *
*  On entry, the different geometries for which the gsf will be        *
*  calculated should be supplied.                                      *
*  On exit, the gsf are returned through the arrays gsfpl and gsfmi.   *
*                                                                      *
*  Referring to the functions Plmn(u) in the above reference, we have  *
*  for uplus = mu*mu0 + sqrt(1-mu**2) sqrt(1-mu0**2) cos(phi-phi0)     *
*  where mu, mu0 and phi-phi0 refer to to geometry number i :          *
*                                                                      *
*     gsfpl(l,1,i) = Pl00(uplus)                                       *
*     gsfpl(l,2,i) = Pl02(uplus)                                       *
*     gsfpl(l,3,i) = Pl22(uplus)                                       *
*     gsfpl(l,4,i) = Pl2-2(uplus)                                      *
*                                                                      *
*  and for umin = -mu*mu0 + sqrt(1-mu**2) sqrt(1-mu0**2) cos(phi-phi0) *
*  where mu, mu0 and phi-phi0 refer to to geometry number i :          *
*                                                                      *
*     gsfmi(l,1,i) = Pl00(umin)                                        *
*     gsfmi(l,2,i) = Pl02(umin)                                        *
*     gsfmi(l,3,i) = Pl22(umin)                                        *
*     gsfmi(l,4,i) = Pl2-2(umin)                                       *
*                                                                      *
*  The method used is recurrence as described in the above reference,  *
*  Eqs. (76) to (82). All equation numbers refer to the above          *
*  reference.                                                          *
*  Notice that all these functions are real !!                         *
*                                            V.L. Dolman March 17 1989 *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter ( pi=3.1415926535897932384D0 )
      parameter ( radfac=pi/180.D0 )
      dimension xmu(NDmu), imu(NDgeom), imu0(NDgeom), phi(NDgeom)
     +        , gsfpl(0:NDcoef,4,NDgeom), gsfmi(0:NDcoef,4,NDgeom)
     +        , sql4(2:NDcoef)
      qroot6=0.25D0*dsqrt(6.D0)
      if (ngeom .le. 0) then
          print *,' setgsf: illegal number of geometries ngeom = ',ngeom
          stop 'setgsf: illegal number of geometries'
      endif
      if (lmax .gt. NDcoef) then
          print *,' setgsf: too many l values requested : ',lmax
          print *,'         maximum NDcoef = ',NDcoef
          stop 'setgsf: too l values requested'
      endif
*----------------------------------------------------------------------*
*   C A S E   W I T H O U T   P O L A R I Z A T I O N                  *
*----------------------------------------------------------------------*
      if (nmat .eq. 1) then
*----------------------------------------------------------------------*
*         Start loop over various geometries                           *
*----------------------------------------------------------------------*
          do 200 i=1,ngeom
              emu = xmu(imu(i))
              emu0= xmu(imu0(i))
              uplus=emu*emu0 + dsqrt((1.D0-emu*emu)*(1.D0-emu0*emu0))
     +                                            * dcos(radfac*phi(i))
              umin=uplus-2.D0*emu*emu0
*----------------------------------------------------------------------*
*             Eq. (77) with m=0                                        *
*----------------------------------------------------------------------*
              gsfpl(0,1,i)=1.D0
              gsfmi(0,1,i)=1.D0
              if (lmax .le. 0) goto 200
*----------------------------------------------------------------------*
*             Eq. (81) with l=0, m=0 and Eq. (76)                      *
*----------------------------------------------------------------------*
              gsfpl(1,1,i)=uplus*gsfpl(0,1,i)
              gsfmi(1,1,i)=umin *gsfmi(0,1,i)
              if (lmax .le. 1) goto 200
*----------------------------------------------------------------------*
*             Eq. (81) with l=1 and m=0                                *
*----------------------------------------------------------------------*
              gsfpl(2,1,i)=1.5D0*uplus*gsfpl(1,1,i)-0.5D0*gsfpl(0,1,i)
              gsfmi(2,1,i)=1.5D0*umin *gsfmi(1,1,i)-0.5D0*gsfmi(0,1,i)
              if (lmax .le. 2) goto 200
*----------------------------------------------------------------------*
*             Start recurrence over index l Eq. (81) with m=0          *
*----------------------------------------------------------------------*
              do 100 l=2, lmax-1
                  fl2 = dble(2*l+1)
                  perl1 =1.D0/dble(l+1)
                  gsfpl(l+1,1,i)=(fl2*uplus*gsfpl(l,1,i)
     +                                  - dble(l)*gsfpl(l-1,1,i) )*perl1
                  gsfmi(l+1,1,i)=(fl2*umin *gsfmi(l,1,i)
     +                                  - dble(l)*gsfmi(l-1,1,i) )*perl1
  100         continue
  200     continue
      else
*----------------------------------------------------------------------*
*     C A S E     W I T H    P O L A R I Z A T I O N                   *
*----------------------------------------------------------------------*
*  Precompute the factors sqrt(l**2-4) needed in Eq. (82) with m=0     *
*----------------------------------------------------------------------*
          sql4(2)=0.D0
          do 300 l=3,lmax
              sql4(l)=sqrt(dble(l)**2-4.D0)
  300     continue
*----------------------------------------------------------------------*
*         Start loop over various geometries                           *
*----------------------------------------------------------------------*
          do 500 i=1,ngeom
              emu = xmu(imu(i))
              emu0= xmu(imu0(i))
              uplus=emu*emu0 + dsqrt((1.D0-emu*emu)*(1.D0-emu0*emu0))
     +                                             * dcos(radfac*phi(i))
              umin=uplus-2.D0*emu*emu0
*----------------------------------------------------------------------*
*             Eq. (77) with m=0                                        *
*----------------------------------------------------------------------*
              gsfpl(0,1,i)=1.D0
              gsfmi(0,1,i)=1.D0
*----------------------------------------------------------------------*
*             Eq. (76)                                                 *
*----------------------------------------------------------------------*
              gsfpl(0,2,i)=0.D0  
              gsfmi(0,2,i)=0.D0  
              gsfpl(0,3,i)=0.D0  
              gsfmi(0,3,i)=0.D0  
              gsfpl(0,4,i)=0.D0  
              gsfmi(0,4,i)=0.D0  
              if (lmax .le. 0) goto 500
*----------------------------------------------------------------------*
*             Eq. (81) with l=0, m=0 and Eq. (76)                      *
*----------------------------------------------------------------------*
              gsfpl(1,1,i)=uplus*gsfpl(0,1,i)
              gsfmi(1,1,i)=umin *gsfmi(0,1,i)
*----------------------------------------------------------------------*
*             Eq. (76)                                                 *
*----------------------------------------------------------------------*
              gsfpl(1,2,i)=0.D0  
              gsfmi(1,2,i)=0.D0  
              gsfpl(1,3,i)=0.D0  
              gsfmi(1,3,i)=0.D0  
              gsfpl(1,4,i)=0.D0  
              gsfmi(1,4,i)=0.D0  
              if (lmax .le. 1) goto 500
*----------------------------------------------------------------------*
*             Eq. (81) with l=1 and m=0                                *
*----------------------------------------------------------------------*
              gsfpl(2,1,i)=1.5D0*uplus*gsfpl(1,1,i)-0.5D0*gsfpl(0,1,i)
              gsfmi(2,1,i)=1.5D0*umin *gsfmi(1,1,i)-0.5D0*gsfmi(0,1,i)
*----------------------------------------------------------------------*
*             Eq. (78) with the plus sign                              *
*----------------------------------------------------------------------*
              gsfpl(2,2,i)=-qroot6*(1.D0-uplus*uplus)
              gsfmi(2,2,i)=-qroot6*(1.D0-umin *umin )
*----------------------------------------------------------------------*
*             Eq. (80) with m=2                                        *
*----------------------------------------------------------------------*
              gsfpl(2,3,i)= 0.25D0*(1.D0+uplus)*(1.D0+uplus)
              gsfmi(2,3,i)= 0.25D0*(1.D0+umin )*(1.D0+umin )
              gsfpl(2,4,i)= 0.25D0*(1.D0-uplus)*(1.D0-uplus)
              gsfmi(2,4,i)= 0.25D0*(1.D0-umin )*(1.D0-umin )
              if (lmax .le. 2) goto 500
*----------------------------------------------------------------------*
*             Start recurrence over index l                            *
*----------------------------------------------------------------------*
              do 400 l=2, lmax-1
*----------------------------------------------------------------------*
*                 Eq. (81) with m=0                                    *
*----------------------------------------------------------------------*
                  fl2 = dble(l+l+1)
                  perl1 =1.D0/dble(l+1)
                  gsfpl(l+1,1,i)=(fl2*uplus*gsfpl(l,1,i)
     +                                  - dble(l)*gsfpl(l-1,1,i) )*perl1
                  gsfmi(l+1,1,i)=(fl2*umin *gsfmi(l,1,i)
     +                                  - dble(l)*gsfmi(l-1,1,i) )*perl1
*----------------------------------------------------------------------*
*                 Eq. (82) with m=0 and the upper sign                 *
*----------------------------------------------------------------------*
                  persq= 1.D0/sql4(l+1)
                  gsfpl(l+1,2,i)=(fl2*uplus*gsfpl(l,2,i)
     +                                  - sql4(l)*gsfpl(l-1,2,i) )*persq
                  gsfmi(l+1,2,i)=(fl2*umin *gsfmi(l,2,i)
     +                                  - sql4(l)*gsfmi(l-1,2,i) )*persq
*----------------------------------------------------------------------*
*                 Eq. (82) with m=2 and both signs                     *
*----------------------------------------------------------------------*
                  fll1=dble(l)*dble(l+1)
                  perll4=1.D0/(dble(l)*sql4(l+1)**2)
                  wplus=fl2*(fll1*uplus-4.D0)
                  wmin =fl2*(fll1*umin -4.D0)
                  q    =dble(l+1)*(dble(l)**2-4.D0)
                  gsfpl(l+1,3,i)=( wplus*gsfpl(l,3,i) 
     +                                       - q*gsfpl(l-1,3,i) )*perll4
                  gsfmi(l+1,3,i)=( wmin *gsfmi(l,3,i) 
     +                                       - q*gsfmi(l-1,3,i) )*perll4
                  wplus=fl2*(fll1*uplus+4.D0)
                  wmin =fl2*(fll1*umin +4.D0)
                  gsfpl(l+1,4,i)=( wplus*gsfpl(l,4,i) 
     +                                       - q*gsfpl(l-1,4,i) )*perll4
                  gsfmi(l+1,4,i)=( wmin *gsfmi(l,4,i) 
     +                                       - q*gsfmi(l-1,4,i) )*perll4
  400         continue
  500     continue
      endif
      return
      end
      subroutine setmu(NDmu,NDgeom, theta, theta0, ngeom, nmug
     +                 ,imu, imu0, xmu, smf, nmutot )
*----------------------------------------------------------------------*
*  Initialize the mu-values and the supermatrixfactors.                *
*  On entry :                                                          *
*      theta    : array containing the viewing angle in degrees        *
*      theta0   : array containing the solar zenith angle in degrees   *
*      ngeom    : number of requested geometries                       *
*      nmug     : number of desired Gauss points                       *
*  On exit :                                                           *
*      xmu      : array containing the mu-values, the Gauss points     *
*                 from 1 to nmug, the extra points from nmug+1 to      *
*                 nmutot. All mu values are different !!!              *
*      smf      : array containing the supermatrix factors,            *
*                 dsqrt(2*w*mu) for Gauss points, 1 for extra points.  *
*      imu      : array containing the index of the position in which  *
*                 the mu value of a certain geometry can be found in   *
*                 array xmu.                                           *
*      imu0     : array containing the index of the position in which  *
*                 the mu0 value of a certain geometry can be found in  *
*                 array xmu.                                           *
*      numtot   : total number of different mu values                  *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter ( pi=3.1415926535897932384D0 )
      parameter ( radfac=pi/180.D0 )
      dimension theta(NDgeom), theta0(NDgeom), xmu(NDmu), smf(NDmu)
     +        , imu(NDgeom), imu0(NDgeom)
      logical found, verbo
      verbo = .false.
      call ad_gauleg( NDmu, nmug, 0.D0, 1.D0, xmu, smf )
      if (verbo) then
          print *,' setmu: the integration mu values are :'
          print *,' '
          print *,'   i       mu(i)              w(i)'
          print *,' -----------------------------------------------'
          do 10 i=1,nmug
              print '(i4,2f20.14)',i, xmu(i), smf(i)
   10     continue
      endif
      nmutot=nmug
      do 200 i=1,ngeom
          emu0 = dabs(dcos(radfac*theta0(i)))
          call member( emu0, xmu, NDmu, nmutot, found, index )
          if (.not. found) then
              nmutot = nmutot+1
              if (nmutot .gt. NDmu) then
                  print *,' setmu: too many mu values; max NDmu=',NDmu
                  stop 'in setmu too many mu values'
              endif
              xmu(nmutot) = emu0
              smf(nmutot) = 1.D0
              imu0(i) = nmutot
          else
              imu0(i) = index
          endif
          emu = dabs(dcos(radfac*theta(i)))
          call member( emu, xmu, NDmu, nmutot, found, index )
          if (.not. found) then
              nmutot = nmutot+1
              if (nmutot .gt. NDmu) then
                  print *,' setmu: too many mu values; max NDmu=',NDmu
                  stop 'in setmu too many mu values'
              endif
              xmu(nmutot) = emu
              smf(nmutot) = 1.D0
              imu(i) = nmutot
          else
              imu(i) = index
          endif
  200 continue
      if (verbo) then
          print *,' '
          print *,' setmu: the extra mu values are :'
          print *,' '
          print *,'   i      mu(i)             w(i)'
          do 210 i=nmug+1,nmutot
              print '(i4,2f20.14)',i, xmu(i), smf(i)
  210     continue
      endif
      do 300 i=1,nmug
          smf(i) = dsqrt(2.D0*smf(i)*xmu(i))
  300 continue
      return
      end
      subroutine sstari(NDmu,xmu,smf,xmut,Rfs,nmug,nmat,xm,URUi )
*----------------------------------------------------------------------*
*  Calculate the spherical reflectance of the water-air interface, and *
*  put it in the right place in array URUi.                            *
*  The desired result can be written as follows [cf. Chowdhary (1991)  *
*  Eq. (44)].                                                          *
*            1                                                         *
*    s* = 2 int dmu mu 2 mu RF*(mu)                                    *
*            0                                                         *
*                                                                      *
*           muc                              1                         *
*       = 2 int dmu mu 2 mu RF*(mu)   +   2 int dmu mu 2 mu RF*(mu)    *
*            0                              muc                        *
*                                                                      *
*         ( muc**2  0      0    0  )         1                         *
*       = (   0    muc**2  0    0  )  +   2 int dmu mu 2 mu RF*(mu)    *
*         (   0     0      Hr   Hi )        muc                        *
*         (   0     0     -Hi   Hr )                                   *
*                                                                      *
*  where muc is defined in Chowdhary (1991) page 37 and Hr and Hi are  *
*  convenient short hands.                                             *
*  Hr and Hi are treated by transform mu" = sqrt(1 - mu/muc) and       *
*  the second term is treated by the transform mu' = mut*(mu) where    *
*  mut*(mu) is defined in Chowdhary (1991) Eq. (12). These transforms  *
*  remove branching points.                                            *
*----------------------------------------------------------------------*
      implicit double precision(a-h,o-z)
      parameter( Nzero=1,  Nall=3
     +         , Nbelow=2 )
      dimension xmu(NDmu), smf(NDmu), xmut(NDmu), Rfs(NDmu,3)
     +        , URUi(4,4,Nbelow,Nall)
      xmuc = dsqrt(1.0d0 - 1.0d0/xm**2)
      Hr = 0.0d0
      Hi = 0.0d0
      do 100 mu=1, nmug
*----------------------------------------------------------------------*
*         Calculate the factors Re(...) and Im(...) occurring in       *
*         Eq. (26) and store them in xRe and xIm. The argument mu is   *
*         replaced by x.                                               *
*----------------------------------------------------------------------*
        if (nmat .gt. 1) then
          x = xmuc * (1.0d0 - xmu(mu)**2)
          xts = xm*dsqrt(xmuc**2 - x**2)
          denom = (x**2 + (xm*xts)**2) * ((xm*x)**2 + xts**2)
          xRe = ( (xm*(x**2 + xts**2))**2 - ((1.0d0-xm**2)*x*xts)**2 )
     +          / denom
          xIm = 2.0d0*xm*(1.0d0-xm**2)*x*xts*(x**2 + xts**2) / denom
*----------------------------------------------------------------------*
*         Integrate to obtain Hr and Hi                                *
*----------------------------------------------------------------------*
          Hr = Hr + 2.0d0*xmuc**2 * smf(mu)**2 * (1.0d0-xmu(mu)**2) *xRe
          Hi = Hi + 2.0d0*xmuc**2 * smf(mu)**2 * (1.0d0-xmu(mu)**2) *xIm
*----------------------------------------------------------------------*
*         Integrate second term                                        *
*----------------------------------------------------------------------*
          URUi(1,2,Nbelow,Nzero) = URUi(1,2,Nbelow,Nzero)
     +         + smf(mu)**2 * 2.0d0*xmut(mu) * Rfs(mu,2) /xm**2
          URUi(3,3,Nbelow,Nzero) = URUi(3,3,Nbelow,Nzero)
     +         + smf(mu)**2 * 2.0d0*xmut(mu) * Rfs(mu,3) /xm**2
        endif
        URUi(1,1,Nbelow,Nzero) = URUi(1,1,Nbelow,Nzero)
     +       + smf(mu)**2 * 2.0d0*xmut(mu) * Rfs(mu,1) /xm**2
  100 continue
*----------------------------------------------------------------------*
*  Combine first and second term                                       *
*  Use symmetry to obtain other non-zero elements                      *
*----------------------------------------------------------------------*
      URUi(1,1,Nbelow,Nzero) = xmuc**2 + URUi(1,1,Nbelow,Nzero)
      URUi(3,3,Nbelow,Nzero) = Hr      + URUi(3,3,Nbelow,Nzero)
      URUi(2,1,Nbelow,Nzero) = URUi(1,2,Nbelow,Nzero)
      URUi(2,2,Nbelow,Nzero) = URUi(1,1,Nbelow,Nzero)
      if (nmat .ge. 4) then
          URUi(3,4,Nbelow,Nzero) =  Hi      
          URUi(4,3,Nbelow,Nzero) = -Hi
          URUi(4,4,Nbelow,Nzero) = URUi(3,3,Nbelow,Nzero)
      endif
      return
      end
c$$$      subroutine terms(NDgeom,  xmu, imu0, imu, nmat, ngeom, Svin
c$$$     +                , R, R1, Tst, T1st, Ri, R1i, Di, D1i
c$$$     +                , Rf, ebbot
c$$$     +                , SLpa, SLskyg, SLbsun, SLbsky 
c$$$     +                , SL1pa, SL1skyg, SL1bsun, SL1bsky )
c$$$*----------------------------------------------------------------------*
c$$$*  Calculate contributions to the Stokes vector reaching the satellite.*
c$$$*  On entry:                                                           *
c$$$*     xmu     : all different mu values (integration and extra points) *
c$$$*     imu0    : array filled by setmu with, for each geometry, the     *
c$$$*               position in array xmu where the mu0 value can be found *
c$$$*     imu     : array filled by setmu with, for each geometry, the     *
c$$$*               position in array xmu where the mu value can be found  *
c$$$*     nmat    : number of elements of the Stokes vector taken into     *
c$$$*               account (4 = full polarization, 3 = 3x3 approximation, *
c$$$*               2 = illegal, 1 = scalar)                               *
c$$$*     ngeom   : number of combinations of viewing and incident         *
c$$$*               directions (number of geometries)                      *
c$$$*     Svin    : incident Stokes vector normalized so that its first    *
c$$$*               component equals 1. The incident light can be written  *
c$$$*               I(0,mu,phi) = pi*F*delta(mu-mu0)*delta(phi-phi0)*Svin  *
c$$$*               pi*F is the flux in the direction of incidence !       *
c$$$*     R       : reflection matrix of the atmosphere                    *
c$$$*     R1      : first order scatterering part of R                     *
c$$$*     Tst     : same as T but for illumination from below              *
c$$$*     T1st    : first order scatterering part of T1st                  *
c$$$*     Ri      : reflection matrix of the atmosphere interface system   *
c$$$*     R1i     : first order scatterering part of Ri                    *
c$$$*     Di      : downward radiation matrix just above the interface     *
c$$$*     D1i     : first order scatterering part of Di                    *
c$$$*      Rf     : Fresnel matrix given by Eq. (22), coded as follows     *
c$$$*                  1,1 element of RF(mu) in Rf(mu,1)                   *
c$$$*                  1,2 element of RF(mu) in Rf(mu,2)                   *
c$$$*                  3,3 element of RF(mu) in Rf(mu,3)                   *
c$$$*               other elements can be found via symmetries of Eq. (22) *
c$$$*     ebbot   : dexp(-b/mu) where b is the optical thickness of        *
c$$$*               the atmosphere                                         *
c$$$*  On exit:                                                            *
c$$$*     SLpa     : atmospheric path radiance vector                      *
c$$$*     SLskyg   : skyglint radiance vector                              *
c$$$*     SLbsun   : background sun radiance vector                        *
c$$$*     SLbsky   : background sky radiance vector                        *
c$$$*     SL1pa    : first order contribution to SLpa                      *
c$$$*     SL1skyg  : first order contribution to SLskyg                    *
c$$$*     SL1bsun  : first order contribution to SLbsun                    *
c$$$*     SL1bsky  : first order contribution to SLbsky                    *
c$$$*----------------------------------------------------------------------*
c$$$      implicit double precision (a-h,o-z)
c$$$      parameter ( pi=3.1415926535897932384D0 )
c$$$      parameter ( twopi=2.D0*pi )
c$$$      parameter ( radfac=pi/180.D0 )
c$$$      dimension xmu(NDmu), imu0(NDgeom), imu(NDgeom)
c$$$      dimension R(4,4,NDgeom),   Tst(4,4,NDgeom)
c$$$     +        , R1(4,4,NDgeom),  T1st(4,4,NDgeom)
c$$$     +        , Ri(4,4,NDgeom),  Di(4,4,NDgeom)
c$$$     +        , R1i(4,4,NDgeom), D1i(4,4,NDgeom)
c$$$     +        , Rf(NDmu,3), ebbot(NDmu), Svin(4)
c$$$      dimension SLpa(4,NDgeom),   SLskyg(4,NDgeom)
c$$$     +        , SLbsun(4,NDgeom), SLbsky(4,NDgeom)
c$$$     +        , SL1pa(4,NDgeom),  SL1skyg(4,NDgeom)
c$$$     +        , SL1bsun(4,NDgeom),SL1bsky(4,NDgeom)
c$$$      dimension SvRi(4),  SvDi(4), SvUi(4), SvR1i(4), SvD1i(4)
c$$$     +        , RwTdir(4,4)
c$$$      do 20 i=1, nmat
c$$$          do 10 j=1, nmat
c$$$              RwTdir(i,j) = 0.0d0
c$$$   10     continue
c$$$   20 continue
c$$$      do 1000 igeom=1, ngeom
c$$$          mu0  = imu0(igeom)
c$$$          mu   = imu(igeom)
c$$$          emu0 = xmu(mu0)
c$$$          emu  = xmu(mu)
c$$$          ebmu = ebbot(mu)
c$$$          ebmu0= ebbot(mu0)
c$$$*----------------------------------------------------------------------*
c$$$*         Set RwTdir to (Rwai T'direct), use symmetries of Fresnel     *
c$$$*         matrix, Eq. (22)                                             *
c$$$*----------------------------------------------------------------------*
c$$$          RwTdir(1,1) = 2.0d0*Rf(mu0,1)*ebmu0
c$$$          RwTdir(2,1) = 2.0d0*Rf(mu0,2)*ebmu0
c$$$          RwTdir(3,3) = 2.0d0*Rf(mu0,3)*ebmu0
c$$$          RwTdir(1,2) = RwTdir(2,1)
c$$$          RwTdir(2,2) = RwTdir(1,1)
c$$$          RwTdir(4,4) = RwTdir(3,3)
c$$$*----------------------------------------------------------------------*
c$$$*         Calculate SLpa, and some auxiliary results :                 *
c$$$*         SvRi : reflected Stokes vector of atmosphere & interface     *
c$$$*         SvDi : downward radiation just above the interface           *
c$$$*         SvUi : unscattered upward radiation just above interface     *
c$$$*----------------------------------------------------------------------*
c$$$          do 200 i=1, nmat
c$$$              SLpa(i,igeom)  = 0.0d0
c$$$              SL1pa(i,igeom) = 0.0d0
c$$$              SvRi(i)  = 0.0d0
c$$$              SvR1i(i) = 0.0d0
c$$$              SvDi(i)  = 0.0d0
c$$$              SvD1i(i) = 0.0d0
c$$$              SvUi(i)  = 0.0d0
c$$$              do 300 j=1, nmat
c$$$                  SLpa(i,igeom)  = SLpa(i,igeom)  
c$$$     +                                      + emu0*R(i,j,igeom)*Svin(j)
c$$$                  SL1pa(i,igeom) = SL1pa(i,igeom) 
c$$$     +                                      + emu0*R1(i,j,igeom)*Svin(j)
c$$$                  SvRi(i)  = SvRi(i)  + emu0*Ri(i,j,igeom) *Svin(j)
c$$$                  SvR1i(i) = SvR1i(i) + emu0*R1i(i,j,igeom)*Svin(j)
c$$$                  SvDi(i)  = SvDi(i)  + emu0*Di(i,j,igeom) *Svin(j)
c$$$                  SvD1i(i) = SvD1i(i) + emu0*D1i(i,j,igeom)*Svin(j)
c$$$                  SvUi(i)  = SvUi(i)  + emu0*RwTdir(i,j)   *Svin(j)
c$$$  300         continue
c$$$  200     continue
c$$$*----------------------------------------------------------------------*
c$$$*         Set RwTdir to (T'direct* Rwai), use symmetries of Fresnel    *
c$$$*         matrix, Eq. (22)                                             *
c$$$*----------------------------------------------------------------------*
c$$$          RwTdir(1,1) = 2.0d0*Rf(mu,1)*ebmu
c$$$          RwTdir(2,1) = 2.0d0*Rf(mu,2)*ebmu
c$$$          RwTdir(3,3) = 2.0d0*Rf(mu,3)*ebmu
c$$$          RwTdir(1,2) = RwTdir(2,1)
c$$$          RwTdir(2,2) = RwTdir(1,1)
c$$$          RwTdir(4,4) = RwTdir(3,3)
c$$$*----------------------------------------------------------------------*
c$$$*         Calculate SLskyg and SLbsun. SLbsky follows from subtraction *
c$$$*         of all other terms from the total signal.                    *
c$$$*----------------------------------------------------------------------*
c$$$          do 500 i=1, nmat
c$$$              SLskyg(i,igeom)  = 0.0d0
c$$$              SLbsun(i,igeom)  = 0.0d0
c$$$              SL1skyg(i,igeom) = 0.0d0
c$$$              SL1bsun(i,igeom) = 0.0d0
c$$$              do 400 j=1, nmat
c$$$                  SLskyg(i,igeom)  = SLskyg(i,igeom) 
c$$$     +                                    + emu*RwTdir(i,j)*SvDi(j)
c$$$                  SL1skyg(i,igeom) = SL1skyg(i,igeom) 
c$$$     +                                    + emu*RwTdir(i,j)*SvD1i(j)
c$$$                  SLbsun(i,igeom)  = SLbsun(i,igeom) 
c$$$     +                                    + emu0*Tst(i,j,igeom)*SvUi(j)
c$$$                  SL1bsun(i,igeom) = SL1bsun(i,igeom) 
c$$$     +                                    + emu0*T1st(i,j,igeom)*SvUi(j)
c$$$  400         continue
c$$$              SLbsky(i,igeom) = SvRi(i) 
c$$$     +               - SLpa(i,igeom) - SLskyg(i,igeom) - SLbsun(i,igeom)
c$$$              SL1bsky(i,igeom) = SvR1i(i) 
c$$$     +            - SL1pa(i,igeom) - SL1skyg(i,igeom) - SL1bsun(i,igeom)
c$$$  500     continue
c$$$ 1000 continue
c$$$      return
c$$$      end
      subroutine top2bot(NDmu, NDsup, NDgeom, nmat, nmutot, ngeom
     +            , Rmtop, Tmtop, Rm1top, Tm1top, ebtop
     +            , Rmbot, Rmsbot, Tmbot, Rm1bot, Rm1sbot, Tm1bot, ebbot
     +            , xRm1top, xRm1stop, xTm1top, xTm1stop
     +            , xRm1bot, xRm1sbot, xTm1bot, xTm1sbot )
*----------------------------------------------------------------------*
*  Copy data pertaining to the top layer to the bottom layer.          *
*----------------------------------------------------------------------*
      implicit double precision(a-h,o-z)
      dimension Rmtop(NDsup,NDsup), Rm1top(NDsup,NDsup)
     +        , Tmtop(NDsup,NDsup), Tm1top(NDsup,NDsup), ebtop(NDmu)
      dimension Rmbot(NDsup,NDsup),  Rmsbot(NDsup,NDsup)
     +        , Rm1bot(NDsup,NDsup), Rm1sbot(NDsup,NDsup)
     +        , Tmbot(NDsup,NDsup),  Tm1bot(NDsup,NDsup), ebbot(NDmu)
      dimension xRm1top(4,4,NDgeom), xRm1stop(4,4,NDgeom)
     +        , xTm1top(4,4,NDgeom), xTm1stop(4,4,NDgeom)
     +        , xRm1bot(4,4,NDgeom), xRm1sbot(4,4,NDgeom)
     +        , xTm1bot(4,4,NDgeom), xTm1sbot(4,4,NDgeom)
      call assign( NDsup,  Rmbot,  Rmtop,  nmat, nmutot )
      call assign( NDsup,  Rm1bot, Rm1top, nmat, nmutot )
      call assign( NDsup,  Tmbot,  Tmtop,  nmat, nmutot )
      call assign( NDsup,  Tm1bot, Tm1top, nmat, nmutot )
      call star(NDsup,  Rmsbot,  Rmbot, nmat, nmutot )
      call star(NDsup,  Rm1sbot, Rm1bot, nmat, nmutot )
      do 300 igeom=1, ngeom
          do 200 k=1, nmat
              do 100 l=1, nmat
                  xRm1bot(l,k,igeom) = xRm1top(l,k,igeom)
                  xTm1bot(l,k,igeom) = xTm1top(l,k,igeom)
                  xRm1sbot(l,k,igeom)= xRm1stop(l,k,igeom)
                  xTm1sbot(l,k,igeom)= xTm1stop(l,k,igeom)
  100         continue
  200    continue
  300 continue
      do 800 mu=1, nmutot
          ebbot(mu) = ebtop(mu)
  800 continue
      return
      end
      subroutine transf(NDsup, S, nsup, nmat )
*----------------------------------------------------------------------*
*  Calculate the product D1 * SM * D1 where D1 is a supermatrix        *
*  containing in each submatrix :                                      *
*         ( 1   0   0   0 )          ( 1   0   0 )                     *
*         ( 0   1   1   0 )    or    ( 0   1   1 )                     *
*         ( 0   1  -1   0 )          ( 0   1  -1 )                     *
*         ( 0   0   0   1 )                                            *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension S(NDsup,NDsup)
      if ((nmat .ne. 3) .and. (nmat .ne. 4)) return
      do 300 i=2, nsup, nmat
          do 100 j=1, nsup
              s1 = S(i,j)
              s2 = S(i+1,j)
              S(i,j)   = s1+s2
              S(i+1,j) = s1-s2
  100     continue
          do 200 j=1, nsup
              s1 = S(j,i)
              s2 = S(j,i+1)
              S(j,i)   = s1+s2
              S(j,i+1) = s1-s2
  200     continue
  300 continue
      return
      end

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

      subroutine addlam( NDgeom
     +                 , A, R, R1, T, T1, Rflux, Tflux, URU, UTU
     +                 , ngeom, nmat
     +                 , RL, R1L, DL, D1L, RLflux, DLflux, URUL )
*----------------------------------------------------------------------*
*  Add an atmosphere to a homogeneous Lambertian surface with          *
*  reflectance A.                                                      *
*  On entry:                                                           *
*     A       : reflectance of Lambertina surface                      *
*     R       : reflection matrix of atmosphere                        *
*     R1      : first order scattering contribution to R               *
*     T       : transmission matrix of atmosphere                      *
*     T1      : first order scattering contribution to T               *
*     Rflux   : reflected fluxes for atmosphere                        *
*     Tflux   : transmitted fluxes for atmosphere                      *
*     URU     : spherical reflectance for atmosphere                   *
*     UTU     : spherical transmittance for the atmosphere             *
*     ngeom   : number of geometries                                   *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*  On exit:                                                            *
*     For the system of atmosphere and Lambertian surface:             *
*     RL      : reflection matrix                                      *
*     R1L     : first order scattering contribution to RL              *
*     DL      : skylight matrix                                        *
*     D1L     : first order scattering contribution to DL              *
*     RLflux  : reflected flux                                         *
*     DLflux  : downward flux above the Lambertian surface             *
*     URUL    : spherical reflectance                                  *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter( Nmu0=1,   Nmu=2
     +         , Nzero=1,  Nfirst=2 , Nall=3
     +         , Nabove=1, Nbelow=2 )
      dimension R(4,4,NDgeom),  T(4,4,NDgeom)
     +        , R1(4,4,NDgeom), T1(4,4,NDgeom)
      dimension Rflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , Tflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , URU(  4,4,              Nbelow, Nall)
     +        , UTU(  4,4,              Nbelow, Nall)
      dimension RL(4,4,NDgeom),  DL(4,4,NDgeom)
     +        , R1L(4,4,NDgeom), D1L(4,4,NDgeom)
      dimension RLflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , DLflux(4,4, NDgeom, Nmu, Nbelow, Nall)
     +        , URUL(  4,4,              Nbelow, Nall)
*----------------------------------------------------------------------*
* Initialize everything to value without Lambertian surface for        *
* illuminatin from above and to zero for illumination from below.      *
*----------------------------------------------------------------------*
      do 1000 igeom=1, NDgeom
          do 200 l=1, nmat
              do 100 k=1, nmat
                  RL(k,l,igeom)  = R(k,l,igeom)
                  DL(k,l,igeom)  = T(k,l,igeom)
                  R1L(k,l,igeom) = R1(k,l,igeom)
                  D1L(k,l,igeom) = T1(k,l,igeom)
  100         continue
  200     continue
          do 900 iorder=Nzero, Nall
              do 500 muarg=Nmu0, Nmu
                  do 400 l=1, nmat
                      do 300 k=1, nmat
                          RLflux(k,l,igeom,muarg,Nabove,iorder)=
     +                            Rflux(k,l,igeom,muarg,Nabove,iorder)
                          DLflux(k,l,igeom,muarg,Nabove,iorder)=
     +                            Tflux(k,l,igeom,muarg,Nabove,iorder)
                          RLflux(k,l,igeom,muarg,Nbelow,iorder)= 0.0d0
                          DLflux(k,l,igeom,muarg,Nbelow,iorder)= 0.0d0
  300                 continue
  400             continue
  500         continue
              do 700 l=1, nmat
                  do 600 k=1, nmat
                      URUL(k,l,Nabove,iorder) = URU(k,l,Nabove,iorder)
                      URUL(k,l,Nbelow,iorder) = 0.0d0
  600             continue
  700         continue
  900     continue
 1000 continue
*----------------------------------------------------------------------*
*  Preparations before loop over geometries. Adding a Lambertian       *
*  surface affects the (1,1), (1,2), (2,1) and (2,2) elements of a     *
*  radiative transfer matrix. For the scalar case only (1,1). This is  *
*  indicated by nStokes.                                               *
*----------------------------------------------------------------------*
      nStokes = 1
      if (nmat .gt. 1) nStokes = 2
      rep  = A/(1.0d0 - A*URU(1,1,Nbelow,Nall))
      rep1 = A/(1.0d0 - A*URU(1,1,Nbelow,Nfirst))
*----------------------------------------------------------------------*
*  Start loop over geometries.                                         *
*----------------------------------------------------------------------*
      do 1300 igeom=1, ngeom
*----------------------------------------------------------------------*
*  Loop over upper left (nStokes X nStokes) block with indices (k,l).  *
*----------------------------------------------------------------------*
          do 1200 l=1, nStokes
              ArepT  = rep  * Tflux(1,l,igeom,Nmu0,Nabove,Nall)
              Arep1T = rep1 * ( Tflux(1,l,igeom,Nmu0,Nabove,Nzero) +
     +                           Tflux(1,l,igeom,Nmu0,Nabove,Nfirst) )
              do 1100 k=1, nStokes
*----------------------------------------------------------------------*
*                 Reflection                                           *
*----------------------------------------------------------------------*
                  RL(k,l,igeom) = R(k,l,igeom) + 
     +                 Tflux(k,1,igeom,Nmu ,Nbelow,Nall) *ArepT
                  R1L(k,l,igeom) = R1(k,l,igeom) + 
     +               ( Tflux(k,1,igeom,Nmu ,Nbelow,Nzero) +
     +                 Tflux(k,1,igeom,Nmu ,Nbelow,Nfirst) )*Arep1T
*----------------------------------------------------------------------*
*                 Skylight                                             *
*----------------------------------------------------------------------*
                  DL(k,l,igeom) = T(k,l,igeom) +
     +                 Rflux(k,1,igeom,Nmu ,Nbelow,Nall) * ArepT
                  D1L(k,l,igeom) = T1(k,l,igeom) +
     +               ( Rflux(k,1,igeom,Nmu ,Nbelow,Nzero) +
     +                 Rflux(k,1,igeom,Nmu ,Nbelow,Nfirst) )*Arep1T
*----------------------------------------------------------------------*
*                 Reflected plane flux                                 *
*----------------------------------------------------------------------*
                  RLflux(k,l,igeom,Nmu0,Nabove,Nall) =  
     +                 Rflux(k,l,igeom,Nmu0 ,Nabove,Nall) +
     +                 UTU(k,1,Nbelow,Nall) * ArepT
                  RLflux(k,l,igeom,Nmu0,Nabove,Nfirst) =  
     +                 Rflux(k,l,igeom,Nmu0,Nabove,Nfirst) +
     +                 ( UTU(k,1,Nbelow,Nzero) 
     +                 + UTU(k,1,Nbelow,Nfirst) ) * Arep1T
*----------------------------------------------------------------------*
*                 Skylight plane flux                                  *
*----------------------------------------------------------------------*
                  DLflux(k,l,igeom,Nmu0,Nabove,Nall) =  
     +                 Tflux(k,l,igeom,Nmu0,Nabove,Nall) +
     +                 URU(k,1,Nbelow,Nall) * ArepT
                  DLflux(k,l,igeom,Nmu0,Nabove,Nfirst) =  
     +                 Tflux(k,l,igeom,Nmu0,Nabove,Nfirst) +
     +                 ( URU(k,1,Nbelow,Nzero) 
     +                 + URU(k,1,Nbelow,Nfirst) ) * Arep1T
*----------------------------------------------------------------------*
*                 Spherical fluxes                                     *
*----------------------------------------------------------------------*
                  if (igeom .eq. 1) then
                    URUL(k,l,Nabove,Nall) = URU(k,l,Nabove,Nall)+
     +                 UTU(k,1,Nbelow,Nall) * rep
     +               * UTU(1,l,Nabove,Nall)
                    URUL(k,l,Nabove,Nfirst) = URU(k,l,Nabove,Nfirst)+
     +                 ( UTU(k,1,Nbelow,Nzero) +
     +                   UTU(k,1,Nbelow,Nfirst) ) * rep1
     +               * ( UTU(1,l,Nabove,Nzero) + 
     +                   UTU(1,l,Nabove,Nfirst) )
                  endif
 1100         continue
 1200     continue   
 1300 continue   
      return 
      end

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

      subroutine Stokout( NDgeom, theta0, nmat, ngeom, Svin, A, SvA )
*----------------------------------------------------------------------*
*  Calculate the Stokes vector SvA obtained by letting the radiative   *
*  transfer matrix A operate on the incident Stokes vector Svin.       *
*  This is done for each geometry.                                     *
*  On entry:                                                           *
*     theta0  : solar zenith angle in degrees for each geometry        *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     ngeom   : number of combinations of viewing and incident         *
*               directions (number of geometries)                      *
*     Svin    : incident Stokes vector normalized so that its first    *
*               component equals 1. The incident light can be written  *
*               I(0,mu,phi) = pi*F*delta(mu-mu0)*delta(phi-phi0)*Svin  *
*               pi*F is the flux in the direction of incidence !       *
*     A       : radiative transfer matrix                              *
*                                                                      *
*  On exit :                                                           *
*     SvA     : Stokes vector mu0 * A * Svin                           *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      parameter ( pi=3.1415926535897932384D0 )
      parameter ( radfac=pi/180.D0 )
      dimension theta0(NDgeom)
      dimension A(4,4,NDgeom), Svin(4), SvA(4,NDgeom)
*----------------------------------------------------------------------*
*  Start loop over geometries                                          *
*----------------------------------------------------------------------*
      do 1000 igeom=1, ngeom
          emu0 = dcos(radfac*theta0(igeom))
          do 200 i=1, nmat
              SvA(i,igeom) = 0.0d0
              do 100 j=1, nmat
                 SvA(i,igeom) = SvA(i,igeom) + emu0*A(i,j,igeom)*Svin(j)
  100         continue
  200     continue
 1000 continue
      return
      end
