
MODULE MBRDF_OCEAN

  USE MCONSTANTS, only : xpi, deg2rad, dp

  IMPLICIT NONE

  PUBLIC  :: brdf_ocean
  PRIVATE :: indwat, &
       morcasiwat, &
       getbound, &
       RMATR

contains 

  !subroutine brdf_ocean(nmat, wl, wind_spd, distr, shadow, xsal, pcl, azw, ts_in, tv_in, fi_in, brdf)
  subroutine brdf_ocean(nmat, wl, wind_spd, distr, shadow, xsal, pcl, azw, Rsw_in, &
       ts_in, tv_in, fi_in, brdf)
    ! Modfied by Mathieu Compiegne sept 2012 : set a proper Fortan 90 syntax
    ! These routines are originally from 6SV code

    ! MC 07/01/2014 : Change of the routine that is used to compute the 
    !                 Stokes reflection matrix for the glitter
    !                 We now use the RMATR routine developped by Mishchenko
    !                 As a consequence we do not use the following routines : 
    !                      sunglint
    !                      Fresnel_mat
    !                      f_to_zmatrix
    !                      Shad_func
    !                      erfc_func
    ! Reference : 
    !               BPDF: Fabienne Maignan, François-Marie Bréon, Emilie Fédèle, Marc Bouvier,
    !               Polarized reflectances of natural surfaces: Spaceborne measurements and analytical modeling 
    !               BRDF avec HOT-SPOT:  F. Maignan, F.-M. Breon, R. Lacaze
    !               Bidirectional reflectance of Earth targets: Evaluation of analytical models using a large 
    !               set of spaceborne measurements with emphasis on the Hot Spot  Remote Sensing of Environment 
    !               90 (2004) 210–220
    !
    !       
    implicit none

    integer,       INTENT(IN) :: nmat
    real(kind=dp), INTENT(IN) :: wl        ! wavelength in microns
    real(kind=dp), INTENT(IN) :: wind_spd  ! wind speed in m/s 
    integer,       INTENT(IN) :: distr     ! wave facet statistical distribution
    !                                        1 : Cox & Munk anisotropic with Gram Charlier series correction term
    !                                        2 : Cox & Munk anisotropic
    !                                        3 : Cox & Munk isotropic
    logical, INTENT(IN) :: shadow
    real(kind=dp), INTENT(IN) :: xsal   ! ocean salinity (ppt)
    real(kind=dp), INTENT(IN) :: pcl    ! pigment concentration (in mg/m^3)
    real(kind=dp), INTENT(IN) :: azw    ! azim. of the sun - azim. of the wind (in deg. between 0 and 360)
    real(kind=dp), INTENT(IN) :: Rsw_in !  Reflectance emerging from sea water (above surface)
    !                                      if -1 --> computed in the routine
    real(kind=dp), INTENT(IN) :: ts_in  ! incident zenith angle (in deg. between 0 and 90)
    real(kind=dp), INTENT(IN) :: tv_in  ! reflected zenith angle (in deg. between 0 and 90)
    real(kind=dp), INTENT(IN) :: fi_in  ! relative azimuth (incident-reflected)
    real(kind=dp), INTENT(OUT), dimension(nmat,nmat) :: brdf ! reflectance of the ocean

    ! local variables 
    real(kind=dp), dimension(4,4) :: res_rmatr ! glint brdf from Mishchenko

    real(kind=dp) :: f_wc 
    real(kind=dp) :: wlp
    integer :: iwl
    real(kind=dp) :: Refwc(39)
    real(kind=dp) :: Rwc 
    real(kind=dp) :: ts, tv, wspd, fi
    real(kind=dp) :: teta_lim

    real(kind=dp) :: nr
    real(kind=dp) :: ni    

    real(kind=dp) :: angbnd(5)
    real(kind=dp) :: wsbnd(6)
    real(kind=dp) :: tdsbnd(5,6)
    real(kind=dp) :: tdvbnd(5,6)
    real(kind=dp) :: rw

    real(kind=8) :: tdv
    real(kind=8) :: tds
    integer :: iws1, iws2
    integer :: ivz1, ivz2
    integer :: isz1, isz2
    real(kind=dp) :: n12
    real(kind=dp) :: a
    real(kind=dp) :: Rsw

    ! FOR CALL TO RMATR
    complex(kind=dp) :: CN1
    complex(kind=dp) :: CN2
    REAL(kind=dp) :: SIGMA2
    REAL(kind=dp) :: DMUI
    REAL(kind=dp) :: DMUR
    REAL(kind=dp) :: PHIR
    REAL(kind=dp) :: PHII

    ! --------------
    ! effective reflectance of the whitecaps (Koepke, 1984)
    Refwc = (/0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190, &
         0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,      &
         0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/)

    angbnd        = (/0.0,45.0,60.0,75.0,85.0/)
    wsbnd         = (/1.0,3.0,5.0,7.0,9.0,20.0/)
    tdsbnd(1:5,1) = (/ 0.9787803,0.9787738,0.9787626,0.9787467,0.9787264/)
    tdsbnd(1:5,2) = (/ 0.9785573,0.9706900,0.9698871,0.9691746,0.9685547/)
    tdsbnd(1:5,3) = (/ 0.9680276,0.9666586,0.9479931,0.9404608,0.9385692/)
    tdsbnd(1:5,4) = (/ 0.9381815,0.9384519,0.9430056,0.9690591,0.9275920/)
    tdsbnd(1:5,5) = (/ 0.9058769,0.8951812,0.8899654,0.8892645,0.9980542/)
    tdsbnd(1:5,6) = (/ 0.9602273,0.9114283,0.8713799,0.8417820,0.7800314/)
    tdvbnd(1:5,1) = (/0.9787764,0.9787535,0.9787106,0.9786453,0.9785548/)
    tdvbnd(1:5,2) = (/0.9775019,0.9692680,0.9637051,0.9564344,0.9495727/)
    tdvbnd(1:5,3) = (/0.9438773,0.9288712,0.9225163,0.9069787,0.9044844/)
    tdvbnd(1:5,4) = (/0.9052351,0.9068328,0.9153687,0.8048478,0.8479503/)
    tdvbnd(1:5,5) = (/0.8678726,0.8797889,0.8878716,0.9091171,0.7294627/)
    tdvbnd(1:5,6) = (/0.8137348,0.8453338,0.8629867,0.8745421,0.9036854/)
    ! --------------

    ! SAME TESTS AS IN CELINE'S CODE
    ! Note: cos(88.3955084734844120 = 0.028 
    !teta_lim = 88.395508473484412
    teta_lim = 90.0 ! means no cut



    if (fi_in.lt.0.0) then
       fi = fi_in + 360.0D0
    else if (fi_in .gt. 360.0D0) then
       fi = fi_in - 360.0D0
    else
       fi = fi_in
    endif
    if (ts_in .gt. teta_lim) then
       ts = teta_lim
       !write(*,*) ' ' 
       !write(*,*) ' (brdf_ocean) Reset theta inc. to higher limit : ', teta_lim
       !write(*,*) ' ' 
    else
       ts = ts_in
    endif
    if (tv_in .gt. teta_lim) then
       tv = teta_lim
       !write(*,*) ' ' 
       !write(*,*) ' (brdf_ocean) Reset theta inc. to higher limit : ', teta_lim
       !write(*,*) ' ' 
    else
       tv = tv_in
    endif
    if (wind_spd .lt. 0.25) then 
       wspd = 0.25
       !write(*,*) ' ' 
       !write(*,*) ' (brdf_ocean) Reset wind speed to lower limit (i.e. 0.25)' 
       !write(*,*) ' ' 
    else
       wspd = wind_spd
    endif
    if ((wl .gt. 4.0) .or. (wl .lt. 0.4)) then
       write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write(*,*) '  (sub. brdf_ocean) ERROR '
       write(*,*) '    The ocean water BRDF is defined '
       write(*,*) '    for 0.4 to 4 microns only'
       write(*,FMT='(A, F10.5)') '    wl = ', wl
       write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       STOP
    endif

    if (distr.ne.3) then
       write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write(*,*) '  (sub. brdf_ocean) ERROR '
       write(*,*) '   Only isotropic Cox & Munk slope distribution '
       write(*,*) '   is available'
       write(*,*) '   --> distr = 3 '
       write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       STOP       
    endif

    ! ----------------------
    ! whitecaps contribution (lambertian)
    ! get the white caps surface fraction coverage
    ! following Monahan and Muircheartaigh (1980):
    f_wc = 2.951d-6 * wspd**3.52d0
    if (f_wc .gt. 1.0d0) f_wc = 1.0d0
    iwl = 1 + int( (wl-0.2) / 0.1 )
    wlp = 0.5_dp + (iwl-1)*0.1_dp
    Rwc = Refwc(iwl+1)+(wl-wlp)/0.1_dp*(Refwc(iwl)-Refwc(iwl+1))
    Rwc = Rwc * f_wc
    !write(*,*) 'Rwc, f_wc=', Rwc, f_wc

    !---------------------
    ! get refractive index
    call indwat(wl,xsal,nr,ni)
    !write(*,*) ''
    !write(*,*) ' in (brdf_ocean) nr = ', nr
    !write(*,*) ' in (brdf_ocean) ni = ', ni
    !write(*,*) ''
    n12 = sqrt(nr*nr+ni*ni)

    ! ----------------------
    ! COMPUTE BACKSCATTERED REFLECTANCE FROM THE SEA WATER (LAMBERTIAN)
    !  water reflectance below the sea surface
   
    if (Rsw_in .eq. -1.0_dp) then
       call morcasiwat(wl, pcl, rw)
       call getbound(wsbnd,1,6,wspd,iws1,iws2)
       call getbound(angbnd,1,5,ts,isz1,isz2)
       call getbound(angbnd,1,5,tv,ivz1,ivz2)
       tds = tdsbnd(isz1,iws1)
       tdv = tdvbnd(ivz1,iws1)
       ! water reflectance above the sea surface
       ! for explanation on value of a see 6SV doc
       a = 0.485_dp
       ! add change in solid angle from under to above to surface
       ! that account for 1/(n12*n12) decrease in sea water directional
       ! reflectance
       Rsw = (1.0_dp/(n12*n12)) * tds * tdv * Rw / (1.0_dp-a*Rw)
       !print*, 'Morel Rsw=', Rsw
    else
       Rsw = Rsw_in
       !print*, '      Rsw=', Rsw
    end if



    !call morcasiwat(wl, pcl, rw)
    !call getbound(wsbnd,1,6,wspd,iws1,iws2)
    !call getbound(angbnd,1,5,ts,isz1,isz2)
    !call getbound(angbnd,1,5,tv,ivz1,ivz2)
    !tds = tdsbnd(isz1,iws1)
    !tdv = tdvbnd(ivz1,iws1)
    ! water reflectance above the sea surface
    ! for explanation on value of a see 6SV doc
    !a = 0.485_dp
    ! add change in solid angle from under to above to surface
    ! that account for 1/(n12*n12) decrease in sea water directional
    ! reflectance
    !Rsw = (1.0_dp/(n12*n12)) * tds * tdv * Rw / (1.0_dp-a*Rw)

    ! ----------------------
    ! glint brdf contribution
    CN1 = CMPLX(1.0, 0.0, kind=dp)
    CN2 = CMPLX(nr, ni, kind=dp)
    SIGMA2 = (0.003 + 5.12d-3 * wspd) / 2.0D0
    DMUI = COS(TS*DEG2RAD)
    PHII = xpi
    DMUR = COS(TV*DEG2RAD)
    PHIR = fi*deg2rad  
    CALL RMATR ( CN1, CN2, SIGMA2, DMUI, PHII, DMUR, PHIR, shadow, res_rmatr) 
 
    ! brdf(1:nmat,1:nmat) = res_rmatr(1:nmat, 1:nmat) ! sun glint only 

    brdf(1:nmat,1:nmat) = 0.0_dp
    brdf(1:nmat,1:nmat) = (1.0_dp-f_wc) * res_rmatr(1:nmat, 1:nmat) ! sun glint 
    brdf(1,1)           = brdf(1,1) + Rwc + ((1-Rwc) * Rsw) ! add foam + undersurface reflectance
    
  end subroutine brdf_ocean

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine indwat(wl,xsal,nr,ni)

    ! Modfied by Mathieu Compiegne sept2012 : set a proper Fortan 90 syntax
    ! input parameters:  wl=
    !                    
    ! output parameters: nr=index of refraction of sea water
    !                    ni=extinction coefficient of sea water

    implicit none
    REAL(kind=dp), INTENT(IN) :: wl   ! wavelength (in micrometers)
    REAL(kind=dp), INTENT(IN) :: xsal ! ocean salinity (ppt)
    REAL(kind=dp), INTENT(OUT) :: nr  ! real part of refractive index
    REAL(kind=dp), INTENT(OUT) :: ni  ! imaginary part of refractive index

    ! local variables
    real(kind=dp) :: twl(62),tnr(62),tni(62)
    real(kind=dp) :: xwl,yr,yi,nrc,nic

    integer i

    ! Indices of refraction for pure water from Hale and Querry, 
    ! Applied Optique, March 1973, Vol. 12,  No. 3, pp. 555-563
    twl = (/0.250,0.275,0.300,0.325,0.345,0.375,0.400,0.425,0.445,0.475,&
         0.500,0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,&
         0.750,0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,&
         1.000,1.200,1.400,1.600,1.800,2.000,2.200,2.400,2.600,2.650,&
         2.700,2.750,2.800,2.850,2.900,2.950,3.000,3.050,3.100,3.150,&
         3.200,3.250,3.300,3.350,3.400,3.450,3.500,3.600,3.700,3.800,&
         3.900,4.000/)
    tnr = (/ 1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336,&
         1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330,&
         1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327,&
         1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219,&
         1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483,&
         1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364,&
         1.357,1.351/)
    tni = (/ 3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09,&
         3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10,&
         1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08,&
         1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08,&
         1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07,&
         3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06,&
         2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04,&
         1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03,&
         1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01,&
         2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01,&
         9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02,&
         1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03,&
         3.80E-03,4.60E-03/)
    i=2
10  if (wl.lt.twl(i)) goto 20
    if (i.lt.62) then
       i=i+1
       goto 10
    endif
20  xwl=twl(i)-twl(i-1)        
    yr=tnr(i)-tnr(i-1)        
    yi=tni(i)-tni(i-1)        
    nr=tnr(i-1)+(wl-twl(i-1))*yr/xwl
    ni=tni(i-1)+(wl-twl(i-1))*yi/xwl

    ! Correction to be applied to the index of refraction and to the extinction 
    ! coefficients of the pure water to obtain the ocean water one (see for 
    ! example Friedman). By default, a typical sea water is assumed 
    ! (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup. 
    ! In that case there is no correction for the extinction coefficient between 
    ! 0.25 and 4 microns. For the index of refraction, a correction of +0.006 
    ! has to be applied (McLellan). For a chlorinity of 19.0ppt the correction 
    ! is a linear function of the salt concentration. Then, in 6S users are able 
    ! to enter the salt concentration (in ppt).
    ! REFERENCES:
    ! Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
    ! McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
    !        New-York, 1965, p 129.
    ! Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
    !        N.J., 1942, p 173.

    nrc=0.006
    nic=0.000
    if (xsal .ge. 0.0) then
       nr=nr+nrc*(xsal/34.3)
       ni=ni+nic*(xsal/34.3)
    else
       nr=nr+nrc
       ni=ni+nic
    endif

  end subroutine indwat

  !------------------------------------------------------------------------------

  subroutine morcasiwat(wl,C,R2)

    ! Spectral diffuse attenuation coefficient of Case I Waters as Predicted 
    ! by MOREL within the spectral range 400-700nm (1988, Journal of Geophysical 
    ! Research, Vol.93, No C9, pp 10749-10768)
    !
    ! input parameters:	wl wavelength (IN MICROMETERS)
    !			C  pigment concentration
    ! output parameter:	R2  reflectance of water
    !
    ! According Morel,1988, we use:
    !
    ! Kd	spectral value of the attenuation coefficient for 
    !	 downwelling irradiance
    !	 with: Kd=Kw+Xc*C**e
    ! Kw	spectral value of the diffuse attenuation coefficient 
    !	 for pure oceanic water
    ! Xc, e	spectral coefficients to compute the diffuse attenuation 
    !	 coefficient for pigment
    ! bb	total backscattering coefficient
    !	 with: bb=0.5*bw+bbt*b
    ! bw	spectral value of the molecular scattering coefficient of water
    ! bbt,b	parameters to compute the scattering coefficients of pigments
    !
    ! R2	reflectance of water below the surface
    !	 with: R2=(0.33/u)*(bb/Kd)	where u is depending of R2

    implicit none

    real(kind=dp), intent(in)  :: wl
    real(kind=dp), intent(in)  :: C
    real(kind=dp), intent(out) :: R2

    real(kind=dp) :: Kw
    real(kind=dp) :: Kd
    real(kind=dp) :: tKw(61)
    real(kind=dp) :: tXc(61)
    real(kind=dp) :: te(61)
    real(kind=dp) :: tbw(61)
    real(kind=dp) :: xc
    real(kind=dp) :: e
    real(kind=dp) :: bw
    real(kind=dp) :: bb
    real(kind=dp) :: b
    real(kind=dp) :: bbt
    real(kind=dp) :: u1
    real(kind=dp) :: r1
    real(kind=dp) :: u2
    real(kind=dp) :: err
    integer       :: iwl

    tKw = (/0.0209,0.0200,0.0196,0.0189,0.0183,&
         0.0182,0.0171,0.0170,0.0168,0.0166,&
         0.0168,0.0170,0.0173,0.0174,0.0175,&
         0.0184,0.0194,0.0203,0.0217,0.0240,&
         0.0271,0.0320,0.0384,0.0445,0.0490,&
         0.0505,0.0518,0.0543,0.0568,0.0615,&
         0.0640,0.0640,0.0717,0.0762,0.0807,&
         0.0940,0.1070,0.1280,0.1570,0.2000,&
         0.2530,0.2790,0.2960,0.3030,0.3100,&
         0.3150,0.3200,0.3250,0.3300,0.3400,&
         0.3500,0.3700,0.4050,0.4180,0.4300,&
         0.4400,0.4500,0.4700,0.5000,0.5500,&
         0.6500/)
    tXc = (/0.1100,0.1110,0.1125,0.1135,0.1126,&
         0.1104,0.1078,0.1065,0.1041,0.0996,&
         0.0971,0.0939,0.0896,0.0859,0.0823,&
         0.0788,0.0746,0.0726,0.0690,0.0660,&
         0.0636,0.0600,0.0578,0.0540,0.0498,&
         0.0475,0.0467,0.0450,0.0440,0.0426,&
         0.0410,0.0400,0.0390,0.0375,0.0360,&
         0.0340,0.0330,0.0328,0.0325,0.0330,&
         0.0340,0.0350,0.0360,0.0375,0.0385,&
         0.0400,0.0420,0.0430,0.0440,0.0445,&
         0.0450,0.0460,0.0475,0.0490,0.0515,&
         0.0520,0.0505,0.0440,0.0390,0.0340,&
         0.0300/)
    te = (/0.668,0.672,0.680,0.687,0.693,&
         0.701,0.707,0.708,0.707,0.704,&
         0.701,0.699,0.700,0.703,0.703,&
         0.703,0.703,0.704,0.702,0.700,&
         0.700,0.695,0.690,0.685,0.680,&
         0.675,0.670,0.665,0.660,0.655,&
         0.650,0.645,0.640,0.630,0.623,&
         0.615,0.610,0.614,0.618,0.622,&
         0.626,0.630,0.634,0.638,0.642,&
         0.647,0.653,0.658,0.663,0.667,&
         0.672,0.677,0.682,0.687,0.695,&
         0.697,0.693,0.665,0.640,0.620,&
         0.600/)
    tbw = (/0.0076,0.0072,0.0068,0.0064,0.0061,&
         0.0058,0.0055,0.0052,0.0049,0.0047,&
         0.0045,0.0043,0.0041,0.0039,0.0037,&
         0.0036,0.0034,0.0033,0.0031,0.0030,&
         0.0029,0.0027,0.0026,0.0025,0.0024,&
         0.0023,0.0022,0.0022,0.0021,0.0020,&
         0.0019,0.0018,0.0018,0.0017,0.0017,&
         0.0016,0.0016,0.0015,0.0015,0.0014,&
         0.0014,0.0013,0.0013,0.0012,0.0012,&
         0.0011,0.0011,0.0010,0.0010,0.0010,&
         0.0010,0.0009,0.0008,0.0008,0.0008,&
         0.0007,0.0007,0.0007,0.0007,0.0007,&
         0.0007/)

    if (wl.lt.0.400_dp.or.wl.gt.0.700_dp) then

       R2 = 0.000_dp

    else

       iwl = 1 + nint((wl-0.400_dp)/0.005_dp)
       Kw  = tKw(iwl)
       Xc  = tXc(iwl)
       e   = te(iwl)
       bw  = tbw(iwl)

       if (abs(C).lt.0.0001_dp)then
          bb = 0.5_dp*bw
          Kd = Kw
       else
          b   = 0.30_dp*C**0.62_dp
          bbt = 0.002_dp+0.02_dp*(0.5_dp-0.25_dp*log10(C))*0.550_dp/wl
          bb  = 0.5_dp*bw+bbt*b
          Kd  = Kw+Xc*C**e
       endif

       u1 = 0.75_dp
       R1 = 0.33_dp*bb/u1/Kd

50     u2  = 0.90_dp*(1.0_dp-R1)/(1.0_dp+2.25_dp*R1)
       R2  = 0.33_dp*bb/u2/Kd
       err = abs((R2-R1)/R2)
       if (err.lt.0.0001_dp) goto 60
       R1 = R2
       goto 50

    endif

60  return

  end subroutine morcasiwat

  !=========================================================

  subroutine getbound(xvals,ifirst,ilast,x,ind1,ind2)

    !     DESCRIPTION: FINDS THE BOUNDS XVALS(IND1) AND XVALS(IND2)
    !                  FOR WHICH X IS BOUNDED.

    implicit none

    real(kind=dp), intent(in) :: xvals(*)
    integer, intent(in) :: ifirst
    integer, intent(in) :: ilast
    real(kind=dp), intent(in) :: x
    integer, intent(out) :: ind1
    integer, intent(out) :: ind2

    integer :: i
    integer :: imid

    !-------

    imid = ilast/2 + 1

    if(xvals(ilast).gt.xvals(ifirst))then
       if(x.gt.xvals(imid))then
          do i=imid,ilast-1
             if(xvals(i).le.x.and.xvals(i+1).ge.x)then
                ind1=i
                ind2=i+1
                goto 15
             endif
          enddo
       else
          do i=ifirst,imid
             if(xvals(i).le.x.and.xvals(i+1).ge.x)then
                ind1=i
                ind2=i+1
                goto 15
             endif
          enddo
       endif
       if(x.lt.xvals(ifirst))then
          ind1=ifirst
          ind2=ifirst+1
       else
          ind1=ilast-1
          ind2=ilast
       endif

    else

       if(x.lt.xvals(imid))then
          do i=imid,ilast-1
             if(xvals(i).ge.x.and.xvals(i+1).le.x)then
                ind1=i
                ind2=i+1
                goto 15
             endif
          enddo
       else
          do i=ifirst,imid
             if(xvals(i).ge.x.and.xvals(i+1).le.x)then
                ind1=i
                ind2=i+1
                goto 15
             endif
          enddo
       endif
       if(x.gt.xvals(ifirst))then
          ind1=ifirst
          ind2=ifirst+1
       else
          ind1=ilast-1
          ind2=ilast
       endif
    endif
15  continue
    return
  end subroutine getbound

  !**************************************************************
  !   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
  !   ILLUMINATION FROM ABOVE FOR
  !   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
  !   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
  !   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS INCLUDED.
  !
  !   FOR ALL FORMULAS AND DEFINITIONS, SEE THE PAPER
  !   M. I. Mishchenko and L. D. Travis, Satellite retrieval
  !   of aerosol properties over the ocean using polarization as well as
  !   intensity of reflected sunlight.  J. Geophys. Res. 102, 16989-
  !   17013 (1997).
  !
  !   INPUT INFORMATION:
  ! 
  !   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER)
  !   DMUI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
  !   PHII = INCIDENT AZIMUTH ANGLE (IN RADIANS)
  !   DMUR = ABS(COSINE OF THE REFLECTION ZENITH ANGLE)
  !   PHIR = REFLECTION AZIMUTH ANGLE IN RADIANS
  !
  !   OUTPUT PARAMETERS:
  !
  !   R - (4X4) REFLECTION MATRIX [R(1,1)=REFLECTION FUNCTION]

  SUBROUTINE RMATR (CN1,CN2,SIGMA2,DMUI_IN,PHII,DMUR_IN,PHIR, SHADOWFLAG, R)

    IMPLICIT NONE

    complex(kind=dp), intent(in) :: CN1
    complex(kind=dp), intent(in) :: CN2
    real(kind=dp), intent(in) :: SIGMA2
    real(kind=dp), intent(in) :: DMUI_IN
    real(kind=dp), intent(in) :: PHII
    real(kind=dp), intent(in) :: DMUR_IN
    real(kind=dp), intent(in) :: PHIR
    logical, intent(in) :: SHADOWFLAG
    real(kind=dp), intent(OUT) :: R(4,4)

    real(kind=dp) :: AF
    real(kind=dp) :: AF11
    real(kind=dp) :: AF12
    real(kind=dp) :: AF21
    real(kind=dp) :: AF22
    real(kind=dp) :: E1
    real(kind=dp) :: E2
    real(kind=dp) :: E3
    real(kind=dp) :: E4
    real(kind=dp) :: SHADOW
    real(kind=dp) :: S1
    real(kind=dp) :: S2
    real(kind=dp) :: S3
    real(kind=dp) :: RDZ4
    real(kind=dp) :: RDZ2
    real(kind=dp) :: PR2
    real(kind=dp) :: PR3
    real(kind=dp) :: PRKI
    real(kind=dp) :: SHADOWI
    real(kind=dp) :: SHADOWR
    real(kind=dp) :: T1
    real(kind=dp) :: T2
    real(kind=dp) :: TI1
    real(kind=dp) :: TI2
    real(kind=dp) :: TI3
    real(kind=dp) :: TIKR
    real(kind=dp) :: TR1
    real(kind=dp) :: TR2
    real(kind=dp) :: TR3
    real(kind=dp) :: TRKI
    real(kind=dp) :: UNIT1
    real(kind=dp) :: UNIT2
    real(kind=dp) :: UNIT3
    real(kind=dp) :: VI1
    real(kind=dp) :: VI2
    real(kind=dp) :: VI3
    real(kind=dp) :: VP1
    real(kind=dp) :: VP2
    real(kind=dp) :: VP3
    real(kind=dp) :: VR1
    real(kind=dp) :: VR2
    real(kind=dp) :: VR3
    real(kind=dp) :: XI
    real(kind=dp) :: XI1
    real(kind=dp) :: XXI
    real(kind=dp) :: DCOEFF
    real(kind=dp) :: DCOSR
    real(kind=dp) :: DCOT
    real(kind=dp) :: DEX
    real(kind=dp) :: DMOD
    real(kind=dp) :: DSI
    real(kind=dp) :: DSINI
    real(kind=dp) :: DSINR
    real(kind=dp) :: DCOSI
    real(kind=dp) :: DSR
    real(kind=dp) :: PI2
    real(kind=dp) :: FACTOR
    real(kind=dp) :: FACT1
    real(kind=dp) :: PI3
    real(kind=dp) :: PIKR
    real(kind=dp) :: PR1
    real(kind=dp) :: P
    real(kind=dp) :: PI1

    real(kind=dp) :: DMUR
    real(kind=dp) :: DMUI

    complex(kind=dp) :: C1
    complex(kind=dp) :: C2
    complex(kind=dp) :: C21
    complex(kind=dp) :: C22
    complex(kind=dp) :: CF11
    complex(kind=dp) :: CF21
    complex(kind=dp) :: CF22
    complex(kind=dp) :: CI
    complex(kind=dp) :: CPTPP
    complex(kind=dp) :: CRPAR
    complex(kind=dp) :: CRPER
    complex(kind=dp) :: CTPPP
    complex(kind=dp) :: CTPPT
    complex(kind=dp) :: CF12
    complex(kind=dp) :: CTTPP
    complex(kind=dp) :: CTTPT
    complex(kind=dp) :: CTTTP
    complex(kind=dp) :: CXI2

    integer :: I, J 

    !   CARTEZIAN COMPONENTS OF THE UNIT VECTORS OF THE INCIDENT AND
    !   SCATTERED BEAMS

    IF (ABS(DMUI_IN-1D0).LT.1D-10) THEN
       DMUI=0.999999999999D0
    ELSE
       DMUI=DMUI_IN
    END IF
    IF (ABS(DMUR_IN-1D0).LT.1D-10) THEN
       DMUR=0.999999999999D0
    ELSE
       DMUR = DMUR_IN
    END IF
    DCOSI=COS(PHII)
    DSINI=DSIN(PHII)
    DCOSR=COS(PHIR)
    DSINR=DSIN(PHIR)
    DSI=DSQRT(1D0-DMUI*DMUI)
    DSR=DSQRT(1D0-DMUR*DMUR)
    VI1=DSI*DCOSI
    VI2=DSI*DSINI
    VI3=-DMUI
    VR1=DSR*DCOSR
    VR2=DSR*DSINR
    VR3=DMUR

    !    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION

    UNIT1=VI1-VR1
    UNIT2=VI2-VR2
    UNIT3=VI3-VR3
    FACT1=UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
    FACTOR=DSQRT(1D0/FACT1)

    !    FRESNEL REFLECTION COEFFICIENTS

    XI1=FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
    CXI2=1D0 - (1D0-XI1*XI1)*CN1*CN1/(CN2*CN2)
    CXI2=CDSQRT(CXI2)
    C1=CN1*XI1
    C2=CN2*CXI2
    CRPER=(C1-C2)/(C1+C2)
    C1=CN2*XI1
    C2=CN1*CXI2
    CRPAR=(C1-C2)/(C1+C2)

    !    CALCULATION OF THE AMPLITUDE SCATTERING MATRIX

    TI1=-DMUI*DCOSI
    TI2=-DMUI*DSINI
    TI3=-DSI

    TR1=DMUR*DCOSR
    TR2=DMUR*DSINR
    TR3=-DSR

    PI1=-DSINI
    PI2=DCOSI
    PI3=0D0

    PR1=-DSINR
    PR2=DCOSR
    PR3=0D0

    PIKR=PI1*VR1+PI2*VR2+PI3*VR3
    PRKI=PR1*VI1+PR2*VI2+PR3*VI3
    TIKR=TI1*VR1+TI2*VR2+TI3*VR3
    TRKI=TR1*VI1+TR2*VI2+TR3*VI3

    E1=PIKR*PRKI
    E2=TIKR*TRKI
    E3=TIKR*PRKI
    E4=PIKR*TRKI

    CF11=E1*CRPER+E2*CRPAR
    CF12=-E3*CRPER+E4*CRPAR
    CF21=-E4*CRPER+E3*CRPAR
    CF22=E2*CRPER+E1*CRPAR

    !   CALCULATION OF THE STOKES REFLECTION MATRIX

    VP1=VI2*VR3-VI3*VR2
    VP2=VI3*VR1-VI1*VR3
    VP3=VI1*VR2-VI2*VR1
    DMOD=VP1*VP1+VP2*VP2+VP3*VP3
    DMOD=DMOD*DMOD

    RDZ2=UNIT3*UNIT3
    RDZ4=RDZ2*RDZ2

    DCOEFF=1D0/(4D0*DMUI*DMUR*DMOD*RDZ4*2D0*SIGMA2)
    DEX= -(UNIT1*UNIT1 + UNIT2*UNIT2)/(2D0*SIGMA2*RDZ2)
    DEX=DEXP(DEX)
    DCOEFF=DCOEFF*FACT1*FACT1*DEX

    AF=0.5D0*DCOEFF
    AF11=ABS(CF11)
    AF12=ABS(CF12)
    AF21=ABS(CF21)
    AF22=ABS(CF22)
    AF11=AF11*AF11
    AF12=AF12*AF12
    AF21=AF21*AF21
    AF22=AF22*AF22

    R(1,1)=(AF11+AF12+AF21+AF22)*AF
    R(1,2)=(AF11-AF12+AF21-AF22)*AF
    R(2,1)=(AF11-AF22+AF12-AF21)*AF
    R(2,2)=(AF11-AF12-AF21+AF22)*AF

    CI=(0D0, -1D0)

    C21=DCONJG(CF21)
    C22=DCONJG(CF22)
    CTTTP=CF11*DCONJG(CF12)
    CTTPT=CF11*C21
    CTTPP=CF11*C22
    CTPPT=CF12*C21
    CTPPP=CF12*C22
    CPTPP=CF21*C22

    R(1,3)=    (-CTTTP-CPTPP)*DCOEFF
    R(1,4)=-CI*( CTTTP+CPTPP)*DCOEFF
    R(2,3)=    (-CTTTP+CPTPP)*DCOEFF
    R(2,4)=-CI*( CTTTP-CPTPP)*DCOEFF
    R(3,1)=    (-CTTPT-CTPPP)*DCOEFF
    R(3,2)=    (-CTTPT+CTPPP)*DCOEFF
    R(3,3)=    ( CTTPP+CTPPT)*DCOEFF
    R(3,4)= CI*( CTTPP-CTPPT)*DCOEFF
    R(4,1)= CI*( CTTPT+CTPPP)*DCOEFF
    R(4,2)= CI*( CTTPT-CTPPP)*DCOEFF
    R(4,3)=-CI*( CTTPP+CTPPT)*DCOEFF
    R(4,4)=    ( CTTPP-CTPPT)*DCOEFF

    !  SHADOWING
    IF (SHADOWFLAG) THEN
       P=DACOS(-1D0)
       S1=DSQRT(2D0*SIGMA2/P)
       S3=1D0/(DSQRT(2D0*SIGMA2))
       S2=S3*S3
       XI=DMUI
       XXI=XI*XI
       DCOT=XI/DSQRT(1D0-XXI)
       T1=DEXP(-DCOT*DCOT*S2)
       T2=DERFC(DCOT*S3)
       SHADOWI=0.5D0*(S1*T1/DCOT-T2)
       XI=DMUR
       XXI=XI*XI
       DCOT=XI/DSQRT(1D0-XXI)
       T1=DEXP(-DCOT*DCOT*S2)
       T2=DERFC(DCOT*S3)
       SHADOWR=0.5D0*(S1*T1/DCOT-T2)
       SHADOW=1D0/(1D0+SHADOWI+SHADOWR)
       DO I=1,4
          DO J=1,4
             R(I,J)=R(I,J)*SHADOW
          END DO
       END DO
    END IF

    RETURN
  END SUBROUTINE RMATR

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!$  subroutine sunglint(nmat, wspd, distr, nr, ni, xsal, wl, azw, ts, tv, fi, rogmat)
!!$
!!$    ! Modfied by Mathieu Compiegne sept2012 : set a proper Fortan 90 syntax
!!$
!!$    implicit none
!!$
!!$    integer, INTENT(IN)       :: nmat
!!$    real(kind=dp), INTENT(IN) :: wspd   ! wind speed (m/s)
!!$    integer, INTENT(IN)       :: distr
!!$    real(kind=dp), INTENT(IN) :: nr     ! real part of sea water refactive index
!!$    real(kind=dp), INTENT(IN) :: ni     ! imaginary part of sea water refactive index
!!$    real(kind=dp), INTENT(IN) :: xsal   ! ocean salinity (ppt)
!!$    real(kind=dp), INTENT(IN) :: wl     ! wavelength in microns
!!$    real(kind=dp), INTENT(IN) :: azw    ! azim. of the sun - azim. of the wind (in deg.)
!!$    real(kind=dp), INTENT(IN) :: ts     ! incident zenith angle (in deg.)
!!$    real(kind=dp), INTENT(IN) :: tv     ! reflected zenith angle (in deg.)
!!$    real(kind=dp), INTENT(IN) :: fi     ! relative azimuth (incident-reflected)
!!$    real(kind=dp), INTENT(OUT), dimension(nmat,nmat) :: rogmat ! reflectance of the ocean
!!$
!!$    real(kind=dp) :: phw
!!$    real(kind=dp) :: cs,cv,ss,sv,phi,zx,zy,tantilt,tilt,proba,xe,xn,xe2,xn2
!!$    real(kind=dp) :: coef,cos2chi,coschi,sinchi
!!$    real(kind=dp) :: sigmaC,sigmaU,C21,C03,C40,C04,C22
!!$    real(kind=dp) :: r11,r21,r33,r43
!!$    real(kind=dp), dimension(nmat,nmat) :: zmat0
!!$    complex(kind=dp) :: epsilo
!!$    real(kind=dp) :: sigma2
!!$    real(kind=dp) :: fact
!!$
!!$    integer i, j
!!$
!!$    ! ==========
!!$
!!$    if ((wl .gt. 4.0) .or. (wl .lt. 0.25)) then
!!$       write(*,*) '(sunglint) The ocean water refractive index is defined '
!!$       write(*,*) '           for 0.25 to 4 microns only'
!!$       STOP
!!$    endif
!!$
!!$    phw = azw*deg2rad
!!$    cs  = cos(ts*deg2rad)
!!$    cv  = cos(tv*deg2rad)
!!$    ss  = sin(ts*deg2rad)
!!$    sv  = sin(tv*deg2rad)
!!$    phi = fi*deg2rad
!!$
!!$    ! ===============================
!!$    ! Cox and Munk slope distribution
!!$    Zx = -sv*sin(phi)/(cs+cv)
!!$    Zy = (ss+sv*cos(phi))/(cs+cv)
!!$    tantilt = sqrt(zx*zx+zy*zy)
!!$    tilt    = atan(tantilt)
!!$
!!$    if (distr .eq. 3) then
!!$       ! isotropic surface
!!$       sigma2 = 0.003 + 5.12d-3 * wspd
!!$       proba = exp(-tantilt**2.0d0 / sigma2) / (xpi * sigma2)
!!$    else
!!$       !  Anisotropic Gaussian distribution
!!$       !    phw = phi_sun - phi_wind
!!$       sigmaC = 0.003d0 + 0.00192d0*wspd
!!$       sigmaU = 0.00316d0*wspd
!!$       xe=(cos(phw)*Zx+sin(phw)*Zy)/sqrt(sigmaC)
!!$       xn=(-sin(phw)*Zx+cos(phw)*Zy)/sqrt(sigmaU)
!!$       xe2 = xe * xe
!!$       xn2 = xn * xn
!!$       if (distr .eq. 1) then
!!$          ! Gram Charlier series correction
!!$          C21 = 0.01d0 - 0.0086d0*wspd
!!$          C03 = 0.04d0 - 0.033d0*wspd
!!$          C40 = 0.40d0
!!$          C22 = 0.12d0
!!$          C04 = 0.23d0
!!$          coef = 1.0d0 - C21/2.0d0*(xe2-1.0d0)*xn-C03/6.0d0*(xn2-3.0d0)*xn
!!$          coef = coef+c40/24.0d0*(xe2*xe2-6.0d0*xe2+3.0d0)
!!$          coef = coef+C04/24.0d0*(xn2*xn2-6.0d0*xn2+3.0d0)
!!$          coef = coef+C22/4.0d0*(xe2-1.0d0)*(xn2-1.0d0)   
!!$       else
!!$          ! no Gram Charlier series correction
!!$          coef = 1.0d0
!!$       endif
!!$       ! slope distribution 
!!$       proba = 1.0d0 / (2.0d0 *xpi *sqrt(sigmaU) * sqrt(sigmaC)) * coef * exp(-(xe2+xn2)/2.d0)
!!$    endif
!!$
!!$    ! ===============================
!!$
!!$    ! Compute Fresnel's coefficient R11
!!$    cos2chi=cv*cs+sv*ss*cos(phi)
!!$    if (cos2chi.gt.1.0d0) cos2chi  =  0.99999999999d0
!!$    if (cos2chi.lt.-1.0d0) cos2chi = -0.99999999999d0
!!$
!!$    coschi=sqrt(0.5d0 * (1.0d0 + cos2chi))
!!$    sinchi=sqrt(0.5d0 * (1.0d0 - cos2chi))
!!$
!!$    ! ===============================
!!$    !  compute  sea Fresnel stokes matrix
!!$    epsilo=cmplx(nr,ni)**2.0d0
!!$    call Fresnel_mat(nmat,coschi,sinchi*sinchi,epsilo,r11,r21,r33,r43)
!!$    ! pour transposer les coefficients de fresnel dans le plan
!!$    ! de reflexion: rotation de (xpi-sig2) and (-sig1)
!!$    call f_to_zmatrix(nmat,r11,r21,r33,r43,cs,ss,cv,sv,phi,zmat0)                   
!!$    ! Compute Reflectance of the sun glint for a rough surface
!!$
!!$    fact = xpi / (cs*cv) * proba / (4.d0*(cos(tilt)**4d0)) 
!!$    do i=1,nmat
!!$       do j=1,nmat
!!$          rogmat(i,j) = fact * zmat0(i,j)   
!!$          !<MC>
!!$          if (abs(rogmat(i,j)) .lt. 1d-300) rogmat(i,j) = 0.0d0
!!$          !<END MC>
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine sunglint
!!$
!!$  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$  subroutine Fresnel_mat(nmat, cost, sin2, epsilo, r11, r21, r33, r43)
!!$
!!$    ! Modfied by Mathieu Compiegne sept2012 : set a proper Fortan 90 syntax
!!$    !  06/06/2012 Modified by C. Cornet to compute all the reflexion Stokes matrix components
!!$    !  to compute the Fresnel's coefficient of reflection (see for
!!$    !  example M. Born and E. Wolf, Principles of Optics, Pergamon Press, fifth
!!$    !  edition, 1975, pp 628
!!$
!!$    implicit none
!!$
!!$    integer, INTENT(IN)      :: nmat
!!$    REAL(kind=dp), INTENT(IN) :: cost      ! cosine of the incident radiation 
!!$    !                                       with respect of the wave facet normal.
!!$    REAL(kind=dp), INTENT(IN) :: sin2      ! sine^2 of the incident radiation 
!!$    !                                       with respect of the wave facet normal.
!!$    COMPLEX(kind=dp), INTENT(IN) :: epsilo ! dielectric constant 
!!$    ! Fresnel's coefficient for reflection :
!!$    REAL(kind=dp), INTENT(OUT) :: r11
!!$    REAL(kind=dp), INTENT(OUT) :: r21
!!$    REAL(kind=dp), INTENT(OUT) :: r33
!!$    REAL(kind=dp), INTENT(OUT) :: r43
!!$
!!$    ! local variables
!!$    complex(kind=dp) :: croot
!!$    complex(kind=dp) :: rr
!!$    complex(kind=dp) :: rl
!!$
!!$    croot=sqrt(epsilo-sin2*(1.0d0,0.d0))
!!$    ! -- fresnel reflexion coefficients (for amplitude)
!!$    rr=(cost-croot)/(cost+croot)
!!$    rl=(croot-epsilo*cost)/(croot+epsilo*cost)    
!!$    ! -- reflexion Stokes matrix components (for intensities)
!!$    r11 = 0.5d0*(rl*conjg(rl)+rr*conjg(rr))        ! is real
!!$    if (nmat .gt. 1) then
!!$       r21 = 0.5d0*(rl*conjg(rl)-rr*conjg(rr))        ! is real
!!$       r33 = 0.5d0*(rl*conjg(rr)+rr*conjg(rl))        ! is real
!!$       r43 = (0.d0,0.5d0)*(rl*conjg(rr)-rr*conjg(rl)) ! is real
!!$    endif
!!$
!!$  end subroutine Fresnel_mat
!!$
!!$  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$  subroutine f_to_zmatrix(nmat,r11,r21,r33,r43,cs,ss,cv,sv,phi,z)
!!$
!!$    ! Modfied by Mathieu Compiegne sept2012 : set a proper Fortan 90 syntax
!!$    !...................................................................
!!$    ! cs,ss,cv,sv cosinus and sinus of the incident and scattered angle relative to the meridian
!!$    !   Fresnel matrix to Scattering matrix
!!$    ! -- rotate polarization planes 
!!$    ! -- sig1: rotation angle from plane of incidence to scattering plane
!!$    ! -- sig2: rotation angle from scattering plane to plane of reflexion
!!$    !...................................................................
!!$
!!$    implicit none
!!$
!!$    ! Fresnel's coefficient for reflection :
!!$
!!$    integer, INTENT(IN)      :: nmat
!!$    REAL(kind=dp), INTENT(IN) :: r11 
!!$    REAL(kind=dp), INTENT(IN) :: r21
!!$    REAL(kind=dp), INTENT(IN) :: r33
!!$    REAL(kind=dp), INTENT(IN) :: r43
!!$    REAL(kind=dp), INTENT(IN) :: cs  ! cos incident zenith angle
!!$    REAL(kind=dp), INTENT(IN) :: ss  ! sin incident zenith angle
!!$    REAL(kind=dp), INTENT(IN) :: cv  ! cos reflected zenith angle
!!$    REAL(kind=dp), INTENT(IN) :: sv  ! sin reflected zenith angle
!!$    REAL(kind=dp), INTENT(IN) :: phi ! relative azimuth incident-reflected in rad
!!$    REAL(kind=dp), INTENT(OUT) :: z(nmat,nmat) ! Fresnel reflection matrix
!!$
!!$    ! local variables
!!$    REAL(kind=dp) :: costheta
!!$    REAL(kind=dp) :: sintheta
!!$    REAL(kind=dp) :: cosi1
!!$    REAL(kind=dp) :: cosi2
!!$    REAL(kind=dp) :: c1, s1
!!$    REAL(kind=dp) :: c2, s2 
!!$
!!$    z(1,1) = r11
!!$
!!$    if (nmat .gt. 1) then 
!!$
!!$       ! -- scattering angle   
!!$       costheta=cv*cs+sv*ss*cos(phi)
!!$       sintheta=sin(acos(costheta))
!!$       ! compute cosi2, hovenier    
!!$
!!$       cosi1 = -32768.0d0
!!$       cosi2 = -32768.0d0
!!$
!!$       if (ss .eq. 0.d0) then ! case thetav=0
!!$          cosi1=dcos(phi)
!!$          if (cs .eq. 0.d0) cosi1=-cosi1
!!$       elseif (sv .eq. 0.d0) then  ! case thetas=0
!!$          cosi2=dcos(phi)
!!$          if (cv .le. 0.0d0) cosi2=-cosi2
!!$       elseif (sintheta .ne. 0.d0) then   
!!$          if (phi .ge. 0.d0) then
!!$             cosi1=(cv-cs*costheta)/(ss*sintheta)
!!$             cosi2=(cs-cv*costheta)/(sv*sintheta)
!!$          else
!!$             cosi1=-(cv-cs*costheta)/(ss*sintheta)
!!$             cosi2=-(cs-cv*costheta)/(sv*sintheta)    
!!$          endif
!!$       else   ! case costheta=1
!!$          cosi1=1.d0 
!!$          cosi2=1.d0
!!$       endif
!!$
!!$       if (cosi2 .gt. 1.d0)  cosi2= 1.d0
!!$       if (cosi2 .lt. -1.d0) cosi2=-1.d0      
!!$       c1=(2.0d0*cosi1*cosi1-1.0d0)            ! cos 2a for the second rotation
!!$       s1=2.0d0*dsqrt(1.0d0-cosi1*cosi1)*cosi1 ! sin 2a for the second rotation
!!$       c2=(2.0d0*cosi2*cosi2-1.0d0)            ! cos 2a for the second rotation
!!$       s2=2.0d0*dsqrt(1.0d0-cosi2*cosi2)*cosi2 ! sin 2a for the second rotation
!!$
!!$       if  (abs(c1) .gt. 1.d0) then   
!!$          !if (abs(c1) .gt. 1.000001) write(*,*) 'c1 > 1',c1
!!$          c1=c1/abs(c1)
!!$          s1=0.d0
!!$       endif
!!$       if  (abs(c2) .gt. 1.d0) then   
!!$          !if (abs(c2) .gt. 1.000001) write(*,*) 'c2 > 1',c2
!!$          c2=c2/abs(c2)
!!$          s2=0.d0
!!$       endif
!!$
!!$       ! -- compute z-matrix = L(xpi=sig2).F.L(-sig1)
!!$
!!$       z(2,1)=r21*c1
!!$       z(1,2)=r21*c2
!!$       z(2,2)=r11*c2*c1-r33*s1*s2
!!$       z(3,1)=-r21*s1
!!$       z(3,2)=-r11*c2*s1-r33*s2*c1
!!$       z(1,3)=r21*s2
!!$       z(2,3)=r11*c1*s2+r33*c2*s1
!!$       z(3,3)=-r11*s1*s2+r33*c2*c1
!!$
!!$       if (nmat .eq. 4) then  
!!$          z(4,1)=0.0d0
!!$          z(4,2)=r43*s2     
!!$          z(4,3)=-r43*c2
!!$          z(1,4)=0.0d0
!!$          z(2,4)=r43*s1
!!$          z(3,4)=r43*c1
!!$          z(4,4)=r33   
!!$       endif
!!$    endif
!!$
!!$  end subroutine f_to_zmatrix
!!$ !=========================================================
!!$
!!$  Real(kind=dp) Function erfc_func(z)
!!$
!!$    Implicit none 
!!$
!!$    Real(kind=dp), intent(in) :: z
!!$
!!$    Real(kind=dp), parameter:: eps=0.000000000001d0
!!$    Real(kind=dp) :: erf
!!$    Real(kind=dp) :: erf_old
!!$    Real(kind=dp) :: erf_n
!!$    Real(kind=dp) :: rel_er
!!$    Integer :: n
!!$
!!$    If (dabs(z) .le. 1.d-50) then
!!$       erf=0.d0
!!$    else 
!!$       If (dabs(z) .gt. 6) then
!!$          if (z .lt. 0) erf=-1.d0
!!$          if (z .gt. 0) erf=1.d0
!!$       else
!!$          rel_er=1.d0
!!$          erf_n=z
!!$          erf=z
!!$          n=0
!!$          Do While (rel_er .gt. eps)
!!$             n=n+1
!!$             erf_old=erf
!!$             erf_n=erf_n*(1.d0-2.d0*n)*z*z/n/(2.d0*n+1.d0)
!!$             erf=erf+erf_n
!!$             rel_er=dabs((erf-erf_old)/erf)
!!$          End do
!!$          erf=2.d0*erf/sqrt(xpi)
!!$       end if
!!$    end if
!!$    erfc_func=1.d0-erf
!!$  End Function erfc_func
!!$
!!$  ! Shadowing function
!!$  ! Taken from Pavel Litvinov
!!$  Subroutine Shad_func(mu0,muv,sigma2,Shadow)
!!$
!!$    Real(kind=dp), intent(in) :: mu0
!!$    Real(kind=dp), intent(in) :: muv
!!$    Real(kind=dp), intent(in) :: sigma2
!!$    Real(kind=dp), intent(out):: Shadow
!!$    Real(kind=dp) :: S1
!!$    Real(kind=dp) :: s3
!!$    Real(kind=dp) :: s2
!!$    Real(kind=dp) :: XI
!!$    Real(kind=dp) :: XXI
!!$    Real(kind=dp) :: DCOT
!!$    Real(kind=dp) :: t1
!!$    Real(kind=dp) :: t2
!!$    Real(kind=dp) :: SHADOWI
!!$    Real(kind=dp) :: SHADOWR
!!$
!!$    ! another shadowing
!!$    S1=SQRT(2.D0*SIGMA2/xpi)
!!$    S3=1.D0/(SQRT(2D0*SIGMA2))
!!$    S2=S3*S3
!!$    XI=mu0
!!$    XXI=XI*XI
!!$    DCOT=XI/SQRT(1D0-XXI)
!!$    T1=EXP(-DCOT*DCOT*S2)
!!$    !CVF doesn't have DERFC. It was replaced by erfc_func
!!$    T2=erfc_func(DCOT*S3)
!!$    SHADOWI=0.5D0*(S1*T1/DCOT-T2)
!!$    XI=muv
!!$    XXI=XI*XI
!!$    DCOT=XI/SQRT(1D0-XXI)
!!$    T1=EXP(-DCOT*DCOT*S2)
!!$    !CVF doesn't have DERFC. It was replaced by erfc_func
!!$    T2=erfc_func(DCOT*S3)
!!$    SHADOWR=0.5D0*(S1*T1/DCOT-T2)
!!$    SHADOW=1D0/(1D0+SHADOWI+SHADOWR)
!!$
!!$  End Subroutine Shad_func

END MODULE MBRDF_OCEAN
