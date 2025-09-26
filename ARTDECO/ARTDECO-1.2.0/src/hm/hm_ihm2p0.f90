!  ##############################################################################
! 
!  CODE 'IHM' (Inhomogeneous Hexagonal Monocrystal):
!  CALCUL DES PROPRIETES OPTIQUES D'UNE PARTICULE DE GLACE
!
!  Base sur la version 2.4 (14/10/2006) du code PHM
!
!  Version IHMcode20 finalisee le 24 fevrier 2007
!  ##############################################################################

subroutine hm_ihm2p0(NPHINC,lambda,RN,RM,LL,RR,RNB,RMB,RNS,RMS,l,KSUIE,    &
     reff_b,veff_b,reff_s,veff_s, ndgt, libdir, &
     IHM_RES, THETAOUT, PGT)

  !<TMP>
  !USE IEEE_ARITHMETIC
  !</TMP>

  !        ************************************************************
  !        *      C.-LABONNOTE Laurent, BROGNIEZ Gerard               *
  !        *     Laboratoire d'Optique Atmospherique  UMR 8518        *
  !        *    Universite des Sciences et Technologies de Lille      *
  !        *         59655 Villeneuve d'Ascq Cedex - France           *
  !        *         E-mail: labon@loa.univ-lille1.fr                 *
  !        ************************************************************
  ! 
  !                           Fevrier 2007
  !---------------------------------------------------------------------
  !            *
  !             \                                /
  !              \                             /
  !               \                          /
  !               _\_______________________/
  !              /\ \                    /  \
  !             /  \                   *     \
  !            /    \_________________________\
  !  *____________   |                        |______________________*
  !           |      |________________________| 
  !            \    /                         /
  !             \  /         /    |     \    / 
  !              \/_________/_____|______\__/
  !                        /      |       \
  !                       /       |        \
  !                      /        |         \
  !                     /         |          \
  !                    *          *           *
  !---------------------------------------------------------------------


  IMPLICIT NONE

  !a Parametres: ################################################################

  REAL(KIND=8), PARAMETER                    :: ATOTi  = 1.00D+00              ! amplitude (totale) incidente;
  REAL(KIND=8), PARAMETER                    :: EPS12  = 1.00D-12              ! 'precision' de plusieurs calculs;
  REAL(KIND=8), PARAMETER                    :: EPS6   = 1.00D-06              ! 'precision' de plusieurs calculs;
  REAL(KIND=8), PARAMETER                    :: EPS8   = 1.00D-08              ! 'precision' de plusieurs calculs;

  !a Variables d'entree et de sortie: ###########################################

  !  MC : some variables are now argument of the subroutine :

  INTEGER, INTENT(IN)                        :: NPHINC                                             ! Nombre de PHotons INCidents;
  INTEGER, INTENT(IN)                        :: NDGT                            ! Nombre total des Domaines Grand-Theta;
  REAL(KIND=8), intent(in)                   :: lambda                          ! longueur d'onde (mic)
  REAL(KIND=8), intent(in)                   :: LL                              ! longueur de la colonne hexagonale [1E-6 m];
  REAL(KIND=8), intent(in)                   :: RR                              ! rayon de la colonne hexagonale [1E-6 m];
  REAL(KIND=8), intent(in)                   :: RM                              ! parties imaginaire et reelle de
  REAL(KIND=8), intent(in)                   :: RN                              ! l'indice de refraction de la glace;

  REAL(KIND=8), intent(in)                   :: reff_b, veff_b                  ! effective variance and radius for air;
  REAL(KIND=8), intent(in)                   :: reff_s, veff_s                  ! for soot;
  REAL(KIND=8), intent(in)                   :: RNB, RMB                        ! real and imag. part of the air bub. relative
                                                                                ! refractive index;
  REAL(KIND=8), intent(in)                   :: RNS, RMS                        ! real and imag. part of the soot bub. relative
                                                                                ! refractive index;
  REAL(KIND=8), intent(in)                   :: KSUIE                           ! Proportion de suie / air (0.0 = only bubbles);
  REAL(KIND=8), intent(in)                   :: l                               ! mean free path length;
  CHARACTER(200), intent(in)                 :: libdir

  !   REAL(KIND=8), DIMENSION(30)                :: readLL                          ! longueur de la colonne hexagonale [1E-6 m];
  !   REAL(KIND=8), DIMENSION(1643)              :: readNU                          ! nombre d'onde (=10000/LAMBDA) [cm-1];
  !   REAL(KIND=8), DIMENSION(1643)              :: readRM                          ! parties imaginaire et reelle de
  !   REAL(KIND=8), DIMENSION(1643)              :: readRN                          ! l'indice de refraction de la glace;
  !   REAL(KIND=8), DIMENSION(30)                :: readRR                          ! rayon de la colonne hexagonale [1E-6 m];

  REAL(KIND=8), intent(out),DIMENSION(4,4,NDGT)  :: PGT                             ! matrice de diffusion, comme fonction de 
  ! l'angle GRANDTHETA de diffusion; 

  REAL(KIND=8), INTENT(OUT), DIMENSION(8)    :: IHM_RES                         ! TABLE CONTAINING RESULTING SCALAR VALUES
  !  IHM_RES(1) =  NPHDIF : Nombre total de PHotons DIFfuses
  !  IHM_RES(2) =  FINC : Flux_d_energie_incidente
  !  IHM_RES(3) =  FDRT : Flux_d_energie_diffusee
  !  IHM_RES(4) =  FDFR : Flux_d_energie_diffractee
  !  IHM_RES(5) =  CEXTGO : Section_efficace_d_extinction
  !  IHM_RES(6) =  CSCAGO : Section_efficace_de diffusion
  !  IHM_RES(7) =  PIZERO : Albedo_pour_une_diffusion
  !  IHM_RES(8) =  GG     : Facteur_d_asymetrie

  REAL(KIND=8), intent(out)                  :: thetaout(NDGT)
  !   CHARACTER(LEN=28)                          :: SORTIE 

  !a Variables intermediaires, fonctions et constantes: #########################

  REAL(KIND=8)                               :: PIZERO                          ! albedo pour une diffusion [adim.];

  INTEGER                                    :: NPHDIF                          ! Nombre total de PHotons DIFfuses;
  REAL(KIND=8)                               :: CSCAGO                          ! sec. effic. totale de diffusion [1E-12 m2];
  REAL(KIND=8)                               :: CEXTGO                          ! sec. effic. totale d'extinction [1E-12 m2];
  REAL(KIND=8)                               :: FDFR                            ! flux d'energie diffractee par la particule;
  REAL(KIND=8)                               :: FINC                            ! flux d'energie incidente sur la particule;
  REAL(KIND=8)                               :: FDRT                            ! flux d'energie diffusee (RT) par la part.;
  REAL(KIND=8)                               :: GG                              ! facteur d'asymetrie, total [adim.];
  REAL(KIND=8)                               :: NU                              ! nombre d'onde (=10000/LAMBDA) [cm-1];

  REAL(KIND=8)                               :: RESOLGT                         ! resolution angulaire de la matrice de
  ! diffusion, en GRANDTHETA [degres];
  INTEGER                                     :: NDGTFR                          ! NDGT 'difFRaction' (resolution RESOLGT);

  REAL(KIND=8)                               :: IPARi,IPERi                     ! intensite incidente: composantes PAR et PER;
  REAL(KIND=8)                               :: IPARd,IPERd                     ! intensite diffusee: composantes PAR et PER;
  REAL(KIND=8)                               :: ReAPARi,ReAPERi                 ! amplitude incidente: parties reelles;
  REAL(KIND=8)                               :: ReAPARd,ReAPERd                 ! amplitude diffusee: parties reelles;
  REAL(KIND=8)                               :: ImAPARd,ImAPERd                 ! amplitude diffusee: parties imaginaires;
  REAL(KIND=8)                               :: ReCRPAR,ImCRPAR                 ! coeff. Fresnel, reflexion, PARallele;
  REAL(KIND=8)                               :: ReCRPER,ImCRPER                 ! coeff. Fresnel, reflexion, PERpendiculaire;
  REAL(KIND=8)                               :: ReCTPAR,ImCTPAR                 ! coeff. Fresnel, transmission, PAR;
  REAL(KIND=8)                               :: ReCTPER,ImCTPER                 ! coeff. Fresnel, transmission, PER;
  REAL(KIND=8), DIMENSION(2,2)               :: ReCSRD                          ! matrice d'amplitudes diffusees
  REAL(KIND=8), DIMENSION(2,2)               :: ImCSRD                          ! (entree de la routine STOKES);
  REAL(KIND=8), DIMENSION(2,2)               :: ReCR                            ! matrice d'amplitudes reflechies
  REAL(KIND=8), DIMENSION(2,2)               :: ImCR                            ! vers l'interieur de la particule;
  REAL(KIND=8)                               :: auxCS                           ! matrice d'amplitudes diffusees
  REAL(KIND=8), DIMENSION(2,2)               :: ReCS                            ! (transmises vers l'interieur puis
  REAL(KIND=8), DIMENSION(2,2)               :: ImCS                            ! refractees vers l'exterieur);
  REAL(KIND=8), DIMENSION(2,2)               :: ReCSA                           ! matrice d'amplitudes diffusees
  REAL(KIND=8), DIMENSION(2,2)               :: ImCSA                           ! (rotation A); 
  REAL(KIND=8), DIMENSION(2,2)               :: ReCSD                           ! matrice d'amplitudes diffusees
  REAL(KIND=8), DIMENSION(2,2)               :: ImCSD                           ! (rotation D); 
  REAL(KIND=8)                               :: auxCT                           ! matrice d'amplitudes transmises
  REAL(KIND=8), DIMENSION(2,2)               :: ReCT                            ! a l'interieur de la particule
  REAL(KIND=8), DIMENSION(2,2)               :: ImCT                            ! (considerant l'attenuation);

  REAL(KIND=8)                               :: ALPHAe                          ! 1er cosinus directeur, faisceau emergent;
  REAL(KIND=8), DIMENSION(0:7)               :: ALPHAf                          ! 1er cosinus directeur, faces de la colonne;
  REAL(KIND=8)                               :: ALPHAi                          ! 1er cosinus directeur, faisceau incident;
  REAL(KIND=8)                               :: ALPHAr                          ! 1er cosinus directeur, faisceau reflechi;
  REAL(KIND=8)                               :: ALPHAt                          ! 1er cosinus directeur, faisceau transmis;
  REAL(KIND=8), DIMENSION(NDGT)              :: ANGLESOLIDE                     ! angles solides elementaires
  ! (GRANDTHETA), integres sur GRANDPHI;
  REAL(KIND=8), DIMENSION(4,4)               :: aux44GTGP                       ! variable auxiliaire (sortie de la routine
  ! STOKES);
  REAL(KIND=8), DIMENSION(4,4,NDGT)          :: aux44GT                         ! variable auxiliaire (matrice de diffusion
  ! brute, avant normalisation);
  REAL(KIND=8)                               :: auxauxDELTA3                    ! variable auxiliaire au calcul de DELTA3;
  REAL(KIND=8), DIMENSION(0:7)               :: auxDELTA3                       ! distances internes a partir d'un impact;
  REAL(KIND=8)                               :: auxSURFAPP                      ! valeur effective de la surface apparente;
  REAL(KIND=8)                               :: auxXimp                         ! ordonnees d'un point d'impact autant
  REAL(KIND=8)                               :: auxYimp                         ! que ORIGINE d'un parcours du photon
  REAL(KIND=8)                               :: auxZimp                         ! a l'interieur de la particule;
  REAL(KIND=8)                               :: AX                              ! = RR * sqrt( 3 ) / 2;
  REAL(KIND=8)                               :: COSANGLEHETERO                  ! cosinus de l'angle d'heterogeneite;
  REAL(KIND=8)                               :: BBETAe                          ! 2eme cosinus directeur, faisceau emergent;
  REAL(KIND=8), DIMENSION(0:7)               :: BBETAf                          ! 2eme cosinus directeur, faces;
  REAL(KIND=8)                               :: BBETAi                          ! 2eme cosinus directeur, faisceau incident;
  REAL(KIND=8)                               :: BBETAr                          ! 2eme cosinus directeur, faisceau reflechi;
  REAL(KIND=8)                               :: BBETAt                          ! 2eme cosinus directeur, faisceau transmis;
  REAL(KIND=8), DIMENSION(NDGT,2)            :: COSDGT                          ! cosinus de DGT(IDGT,1 et 2);
  REAL(KIND=8)                               :: COSETA                          ! cosinus de l'angle ETA;
  REAL(KIND=8)                               :: COSPHI                          ! cosinus de l'angle PHI;
  REAL(KIND=8)                               :: COSTHETA                        ! cosinus de l'angle THETA;
  REAL(KIND=8)                               :: COSTHETAi                       ! cosinus de l'angle THETAi d'incidence;
  REAL(KIND=8)                               :: COSTHETAt                       ! cosinus de l'angle THETAt de refraction; 
  REAL(KIND=8)                               :: COSGRANDTHETA                   ! cosinus de l'angle GRANDTHETA;
  REAL(KIND=8)                               :: GRANDTHETA                      ! angle GRANDTHETA;
  REAL(KIND=8)                               :: DATT                            ! partie reelle de l'attenuation;
  REAL(KIND=8)                               :: DCOSA                           ! var. aux. (rotation A);
  REAL(KIND=8)                               :: DCOSB                           ! var. aux. (rotation B);
  REAL(KIND=8)                               :: DCOSD                           ! var. aux. (rotation D);
  REAL(KIND=8)                               :: DDDD                            ! distance interne maximale [1E-12 m2];
  REAL(KIND=8)                               :: DELTA3                          ! distance entre deux impacts internes;
  REAL(KIND=8)                               :: DENOMGG                         ! integrations au DENOMinateur et au
  REAL(KIND=8)                               :: NUMERGG                         ! NUMERateur dans le facteur d'asymetrie;
  REAL(KIND=8), DIMENSION(NDGT,2)            :: DGT                             ! Domaines de l'angle GRANDTHETA;
  REAL(KIND=8)                               :: DSINA                           ! var. aux. (rotation A);
  REAL(KIND=8)                               :: DSINB                           ! var. aux. (rotation B);
  REAL(KIND=8)                               :: DSIND                           ! var. aux. (rotation D);
  REAL(KIND=8), allocatable                  :: EDFR(:)                           ! distr. angul. de l'energie diffractee;
  REAL(KIND=8), allocatable                  :: maxEDFR(:)                         ! var. aux. (EDFR);
  REAL(KIND=8), allocatable                  :: minEDFR(:)                        ! var. aux. (EDFR);
  REAL(KIND=8)                               :: ETA                             ! angle entre le faisceau incident et l'axe
  ! OZ' [degres];
  REAL(KIND=8)                               :: maxETA                          ! var. aux. (EDFR);
  REAL(KIND=8)                               :: minETA                          ! var. aux. (EDFR);

  REAL(KIND=8)                               :: GAMMAe                          ! 3eme cosinus directeur, faisceau emergent;
  REAL(KIND=8), DIMENSION(0:7)               :: GAMMAf                          ! 3eme cosinus directeur, faces;
  REAL(KIND=8)                               :: GAMMAi                          ! 3eme cosinus directeur, faisceau incident;
  REAL(KIND=8)                               :: GAMMAr                          ! 3eme cosinus directeur, faisceau reflechi;
  REAL(KIND=8)                               :: GAMMAt                          ! 3eme cosinus directeur, faisceau transmis;
  REAL(KIND=8)                               :: GRANDK                          ! variable intermediaire (indice equiv.);
  REAL(KIND=8)                               :: GRANDN                          ! variable intermediaire (indice equiv.);
  REAL(KIND=8)                               :: KKKK                            ! = 2*PI/(LongueurD'Onde);
  REAL(KIND=8)                               :: LL2                             ! la moitie de la longueur de la colonne;
  REAL(KIND=8)                               :: PHI                             ! angle entre la projection OXp du faisceau
  ! incident 
  ! et l'axe OX [adim.];
  REAL(KIND=8)                               :: maxPHI                          ! var. aux. (EDFR);
  REAL(KIND=8)                               :: minPHI                          ! var. aux. (EDFR);
  REAL(KIND=8)                               :: rnETA                           ! "Random Number" de l'angle ETA;
  REAL(KIND=8)                               :: rnPHI                           ! "Random Number" de l'angle PHI;
  REAL(KIND=8)                               :: rnPOLAR                         ! polarisation aleatoire du champ incident;
  REAL(KIND=8)                               :: rnZETA                          ! "Random Number" de l'angle ZETA;
  REAL(KIND=8), DIMENSION(3,10)              :: SE                              ! coordonnees XYZ des dix sommets eclaires;
  REAL(KIND=8)                               :: SINPHI                          ! sinus de l'angle PHI;
  REAL(KIND=8)                               :: SINTHETA                        ! sinus de l'angle THETA;
  REAL(KIND=8)                               :: sumP11                          ! var. aux. (normalis. de P11 ray-tracing);
  REAL(KIND=8)                               :: XND,YND,ZND                     ! variables auxiliaires
  REAL(KIND=8)                               :: XNI,YNI,ZNI                     ! (cosinus directeurs divers);
  REAL(KIND=8)                               :: XNP,YNP,ZNP                     !     "
  REAL(KIND=8)                               :: XNR,YNR,ZNR                     !     "
  REAL(KIND=8)                               :: XNS,YNS,ZNS                     !     "
  REAL(KIND=8)                               :: XTD,YTD,ZTD                     !     "
  REAL(KIND=8)                               :: XTI,YTI,ZTI                     !     "
  REAL(KIND=8)                               :: XTR,YTR,ZTR                     !     "
  REAL(KIND=8)                               :: XTS,YTS,ZTS                     !     "
  REAL(KIND=8)                               :: Ximp,Yimp,Zimp                  ! (position) impact d'un photon sur une face;
  REAL(KIND=8)                               :: ZETA                            ! angle entre les axes OZ et OZ' [adim.];
  REAL(KIND=8)                               :: maxZETA                         ! var. aux. (EDFR);
  REAL(KIND=8)                               :: minZETA                         ! var. aux. (EDFR);

  INTEGER                                    :: auxFACEimp                      ! variable auxiliaire (identif. de FACEimp);
  INTEGER                                    :: auxNPNPT                        ! variable auxiliaire (routine DIFFRACT);
  INTEGER                                    :: ICH                             ! indicateur des colonnes hexagonales;
  INTEGER                                    :: IDGT                            ! indicateur du Domaine GrandTheta (NDGT);
  INTEGER                                    :: IFACE                           ! indicateur des face de la particule;
  INTEGER                                    :: INU                             ! indicateur des nombres d'onde d'interet;
  INTEGER                                    :: IPHINC                          ! Indicateur des PHotons INCidents;
  INTEGER                                    :: IREFL                           ! Indicateur de REFLexion totale;
  INTEGER                                    :: JP                              ! indicateurs dans les matrices
  INTEGER                                    :: KP                              ! AUX22GTGP, aux44GTGP, aux44GT et PGT;
  INTEGER                                    :: FACEimp                         ! indicateur de la face atteinte (impactee);

  REAL(KIND=8)                               :: surfapp                         ! fonction interne 'surface apparente de 
  ! la colonne hexagonale';
  REAL(KIND=8)                               :: sumSURFAPP                      ! somme des valeurs auxSURFAPP;
  REAL(KIND=8)                               :: maxSURFAPP                      ! SURFace APParente maximale; [1E-12 m2];
  REAL(KIND=8)                               :: intSURFAPP                      ! SURFace APParente intermediaire [1E-12 m2];
  REAL(KIND=8)                               :: minSURFAPP                      ! SURFace APParente minimale; [1E-12 m2];

  !<CHG,PHM version 2.4>
  REAL(KIND=8)                               :: rnmu
  REAL(KIND=8)                               :: mu
  REAL(KIND=8)                               :: THETA
  REAL(KIND=8)                               :: SINETA
  REAL(KIND=8)                               :: COSZETA
  !</CHG,PHM version 2.4>

  ! IHM variables ##############################################################
  CHARACTER(50)                              :: wavdir, slam
  CHARACTER(250)                             :: commande
  INTEGER, PARAMETER                         :: numr = 101, nmumax = 130,      &
       nmu = 125
  INTEGER                                    :: ISB                             ! inclusion type indicator;
  INTEGER                                    :: jr                              ! radius indicator;
  INTEGER                                    :: jt                              ! scattered angle indicator
  INTEGER                                    :: j
  INTEGER                                    :: Flag_face                       ! Flag that indicate weather or not the 
  ! photon is on a cristal face;
  INTEGER                                    :: cont                            ! Counter of reached inclusion
  REAL(KIND=8)                               :: rmin_b, rmax_b                  ! min and max effec. radius for air;
  REAL(KIND=8)                               :: rmin_s, rmax_s                  ! min and max effec. radius for soot;
  REAL(KIND=8)                               :: rmin, rmax                      ! min and max effec. radius for inclusion;
  REAL(KIND=8), DIMENSION(2,numr)            :: tab_r                           ! table of radius for air and soot;
  REAL(KIND=8)                               :: rnk, rnl                        ! Random number between [0,1];

  REAL(KIND=8)                               :: DRNI, DINI                      ! real and imag. part of the inclusion relative
  ! refractive index;
  REAL(KIND=8)                               :: HL, HLtot                       ! Distance parcourue avant de rencontrer une 
  ! inclusion, et distance total;
  REAL(KIND=8)                               :: r_inc                           ! Rayon de l'inclusion apres tirage aleatoire;
  REAL(KIND=8)                               :: ALPHAn                          ! 1er cos. directeur apres diff. par une inclu; 
  REAL(KIND=8)                               :: BBETAn                          ! 2ie cos. directeur apres diff. par une inclu; 
  REAL(KIND=8)                               :: GAMMAn                          ! 3ie cos. directeur apres diff. par une inclu;
  REAL(KIND=8)                               :: XNRn, YNRn, ZNRn                ! Cos. director of the normal to the new plan;
  REAL(KIND=8)                               :: XTRn, YTRn, ZTRn                ! Cos. director of the orthonormal  "       ";
  REAL(KIND=8)                               :: COSAEi, SINAEi                  ! Cos. and sin. of the angle needed to rotate 
  ! the incident electric field on the inclusion
  ! in the scattering plan;
  REAL(KIND=8)                               :: thetad, phid                    ! Angles de diffusion apres une inclusion;
  REAL(KIND=8)                               :: NRJ_Ei, NRJ_Ed                  ! Incident and scattered elec. field energy;
  REAL(KIND=8)                               :: NRJ_Et, NRJ_Etb, NRJ_Eta
  REAL(KIND=8), DIMENSION(numr)              :: NNr_b, NNr_s, NNr               ! Cumulative inclusion size distrib.;
  REAL(KIND=8), DIMENSION(2,numr)            :: PIZ                             ! Tableau d'w0 des inclusions;
  REAL(KIND=8), DIMENSION(-nmumax:nmumax)    :: RMU, PMU                        ! Gauss point and weigth;
  REAL(KIND=8), DIMENSION(2,numr,-nmumax:nmumax)                              &
       :: PTETA,                        & ! Tableau de Fonction de phase inclusion;
       Pcum,                         & ! Tableau de fct de phase cumulee;
       DRES1,                        & ! Tableau part. reel. de S1 (coef1 mat. ampl.);
       DRES2,                        & ! Tableau part. reel. de S2 (coef1 mat. ampl.);
       DIMS1,                        & ! Tableau part. imag. de S1 (coef1 mat. ampl.);
       DIMS2                           ! Tableau part. imag. de S2 (coef1 mat. ampl.);
  REAL(KIND=8), DIMENSION(2,2)               :: Matrot                          ! Matrice de rotation
  REAL(KIND=8), DIMENSION(2,2)               :: ReCS_rot, ImCS_rot              ! real and imaginary part of the phase matrix 
  ! after rotation

  REAL(KIND=8), PARAMETER                    ::                               &
       PI     = 3.1415926535897932384626433832795029D+00  ! = constante PI;
  REAL(KIND=8), PARAMETER                    :: PI2    = PI / 2.0D+00           ! = PI / 2;
  REAL(KIND=8), PARAMETER                    :: PI3    = PI / 3.0D+00           ! = PI / 3;
  REAL(KIND=8), PARAMETER                    :: PI6    = PI / 6.0D+00           ! = PI / 6;
  REAL(KIND=8), PARAMETER                    :: PI180  = PI / 180.0D+00         ! = PI / 180;
  REAL(KIND=8), PARAMETER                    ::                               &
       R32    = 0.8660254037844385965883020617184229D+00  ! = sqrt( 3 )/2;

  !<TMP>
  REAL(KIND=8)                               :: ALPHAe_prev,BBETAe_prev,GAMMAe_prev
  REAL(KIND=8), DIMENSION(3)                 :: v_old, v, w_old, w, u_old, u
  !</TMP>

  !  added by MC to get Gauss angle and weight
  REAL(KIND=8) :: rmu_tmp(2*nmu+1)
  REAL(KIND=8) :: pmu_tmp(2*nmu+1)

  !======================================================================

  !---------------------
  !  modif MC to allow for NDGT control as an input argument 
  RESOLGT = 180.D0 / (NDGT - 1)
  NDGTFR  = idnint( 1.0D+00 +  90.0D+00 / RESOLGT )
  allocate(EDFR(NDGTFR),maxEDFR(NDGTFR),minEDFR(NDGTFR))

  !a Lecture des variables d'entree (1): ########################################
  !  MC : to replace the Gauss weight and angles that were read in a file
  call ihm_gauleg(-1.0,1.0,RMU_tmp,PMU_tmp,2*nmu+1,2*nmu+1)
  do j = 1, 2*nmu+1
     !write(*,*) RMU_tmp(j), PMU_tmp(j)
     pmu(j-1-nmu) = PMU_tmp(j)
     rmu(j-1-nmu) = RMU_tmp(j)
  enddo

  NU   = 1.0D4 / lambda
  KKKK = 2.0D+00 * PI / ( lambda )
  DDDD = dsqrt( LL * LL + 4.0D+00 * RR * RR )

  LL2 = LL / 2.0D+00
  AX  = RR * R32

  !******************************************************************************
  !   Create the wave lib directory for the given wavelength
  !******************************************************************************

  !--- wav string and dir
  IF (lambda .LT. 10.0D0) THEN
     WRITE(slam,'(a3,i1,f5.3)')'wav',0,lambda
  ELSE
     WRITE(slam,'(a3,f6.3)')'wav',lambda
  ENDIF
  wavdir = libdir(:len_trim(libdir)) // slam(:len_trim(slam))

  !--- If this directory does not exist create it!
  commande = 'mkdir ' // wavdir(:len_trim(wavdir))
  CALL SYSTEM (commande)

  !a Ouverture du fichier de sortie: ############################################

!!$  write(*,*) ''
!!$  write(*,*) '############################################################'
!!$  write(*,*) '                          IHM '
!!$  write(*,*) ''
!!$  write(*,*)
!!$  write(*,*) ' CARACTERISTIQUES_GEOMETRIQUES_DU_PROBLEME:'
!!$  write(*,'('' Longueur_de_la_colonne_(1E-6_m)_= '',1e12.6)') &
!!$       LL
!!$  write(*,'('' Rayon_de_la_colonne_(1E-6_m)_= '',1e12.6)') &
!!$       RR
!!$  write(*,'('' Rayon_de_la_sphere_equiv._(1E-6_m)_= '',1e12.6)') &
!!$       ( 9.0D+00 * RR * RR * LL * dsqrt( 3.0D+00 ) /        &
!!$       ( 8.0D+00 * PI ) )**( 1.0D+00 / 3.0D+00 )
!!$  write(*,'('' Facteur_de_forme_=_Longueur_/_Largeur_= '',1e12.6)') &
!!$       LL / ( 2.0D+00 * RR )
!!$  IF (KSUIE .EQ. 0.0D0) THEN
!!$     WRITE(*,*)"type d'inclusion: air"
!!$     WRITE(*,'(a,f6.3)')" rayon effecif = ",reff_b
!!$     WRITE(*,'(a,f6.3)')" variance effective = ",veff_b
!!$  ELSE
!!$     WRITE(*,*)"type d'inclusion: suie"
!!$     WRITE(*,'(a,f6.3)')"rayon effectif = ",reff_s
!!$     WRITE(*,'(a,f6.3)')"varriance effective = ",veff_s
!!$  ENDIF
!!$  WRITE(*,'(a,f6.2)')" libre parcours moyen : ",l
!!$  write(*,*)
!!$  write(*,*) ' CARACTERISTIQUES_OPTIQUES_DU_PROBLEME:'
!!$  write(*,'('' Longueur d onde = '',f7.3)')&
!!$       lambda
!!$  write(*,'('' Nombre_d_onde_(cm-1)_= '',1e12.6)') &
!!$       NU
!!$  write(*,'('' Partie_reelle_de_l_indice_de_refr._= '',1e12.6)') &
!!$       RN
!!$  write(*,'('' Partie_imaginaire_de_l_indice_de_refr._= '',1e12.6)') &
!!$       RM
!!$
!!$  write(*,*)
!!$  write(*,*) ' CARACTERISTIQUES_DU_FAISCEAU_INCIDENT:'
!!$  write(*,'('' Nombre_de_photons_incidents_= '',i9)') &
!!$       NPHINC
!!$  write(*,'('' Amplitude_du_champ_electrique_incident_= '',1e12.6)') &
!!$       1.0D+00
!!$  write(*,*) 'Polarisation_du_champ_electrique_incident_= ', &
!!$       'aleatoire'
!!$  write(*,*)
!!$  write(*,*) ' ORIENTATIONS_DES_FAISCEAUX_EMERGENTS:'
!!$  write(*,'('' Resol._de_la_matrice_de_diffusion_(degres)_= '',1e12.6)') &
!!$       RESOLGT
  !</DOIT BOUGER DS LA SUBROUTINE Write_outputs>!!!!!!!!!!!!!!!!!!

  !a Initialisation de plusieurs variables: #####################################

  CALL Init_var(nmumax,numr,NDGT,NDGTFR,NPHDIF,sumSURFAPP,FDRT,FINC,FDFR,    &
       NNr_b,NNr_s,PIZ,PTETA,Pcum,DRES1,DRES2,DIMS1,DIMS2,aux44GT,  &
       EDFR)

  !a Angle de diffusion GRANDTHETA: #############################################
  !a
  !a    DGT(IDGT,1) = angle inferieur du domaine;
  !a    DGT(IDGT,2) = angle superieur du domaine;

  DGT(1,1) = 0.0D+00
  DGT(1,2) = RESOLGT / 2.0D+00
  DO IDGT = 2, NDGT-1
     DGT(IDGT,1) = DFLOAT( IDGT - 1 ) * RESOLGT - RESOLGT / 2.0D+00
     DGT(IDGT,2) = DFLOAT( IDGT     ) * RESOLGT - RESOLGT / 2.0D+00
  ENDDO !IDGT
  DGT(NDGT,1) = 180.0D+00 - RESOLGT / 2.0D+00
  DGT(NDGT,2) = 180.0D+00
  DO IDGT = 1, NDGT
     COSDGT(IDGT,1) = DCOS( DGT(IDGT,1) * PI180 )
     COSDGT(IDGT,2) = DCOS( DGT(IDGT,2) * PI180 )
  ENDDO !IDGT

  !a Les sommets de la particule qui sont eclaires: #############################
  !a
  !a    SE(1,I), SE(2,I) et SE(3,I): ordonnees X, Y et Z du I-eme sommet eclaire'

  CALL Sommets(RR,LL,AX,SE)    

  !a Cosinus directeurs des faces de la particule: ##############################

  DO IFACE = 0, 5
     ALPHAf(IFACE) = DCOS( DFLOAT( IFACE ) * PI3 )
     BBETAf(IFACE) = DSIN( DFLOAT( IFACE ) * PI3 )
     GAMMAf(IFACE) = 0.0D+00
  ENDDO !IFACE
  DO IFACE = 6, 7
     ALPHAf(IFACE) = 0.0D+00
     BBETAf(IFACE) = 0.0D+00
     GAMMAf(IFACE) = DCOS( DFLOAT( IFACE - 6 ) * PI )
  ENDDO !IFACE

  !a Surfaces apparentes extremes et leurs orientations: ########################

  CALL EXTREMES(RR,LL,PI180,R32,maxSURFAPP,minSURFAPP,maxETA,minETA,maxPHI,  &
       minPHI,maxZETA,minZETA)

  !a Contributions extremes associees a la diffraction: #########################

  auxNPNPT = 8 
  IF ( minSURFAPP .GT. (1.0D+06) ) auxNPNPT = 10

  CALL DIFFRACT(NU,KKKK,LL,RR,SE,maxSURFAPP,maxETA,maxZETA,maxPHI,NDGT,      &
       NDGTFR,RESOLGT,COSDGT,EPS12,PI,PI180,PI3,auxNPNPT,maxEDFR)

  CALL DIFFRACT(NU,KKKK,LL,RR,SE,minSURFAPP,minETA,minZETA,minPHI,NDGT,      &
       NDGTFR,RESOLGT,COSDGT,EPS12,PI,PI180,PI3,auxNPNPT,minEDFR)

  !  Initialisation de la fonction random
  CALL RANDOM_SEED

  !******************************************************************************
  !   Calcul des granulometries cumulees des bulles d'air (NNr_b) et/ou de la 
  !   suie (NNr_s)
  !******************************************************************************
  !   air bubbles
  CALL granulo_cumul(numr,reff_b,veff_b,lambda,NNr_b,rmin_b,rmax_b,          &
       tab_r(1,:))
  !   soot bubbles
  CALL granulo_cumul(numr,reff_s,veff_s,lambda,NNr_s,rmin_s,rmax_s,          &
       tab_r(2,:))

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !a                  Boucle sur les photons incidents ##########################
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  cont = 0          ! Counter of reached inclusions

  DO IPHINC = 1, NPHINC

     !<TMP>
     !      IF ( IPHINC .EQ. 100000) WRITE(6,*)'IPHINC=',IPHINC
     !      IF ( IPHINC .EQ. 300000) WRITE(6,*)'IPHINC=',IPHINC
     !      IF ( IPHINC .EQ. 700000) WRITE(6,*)'IPHINC=',IPHINC
     !</TMP>

     !--- Initialisation de la distance total parcourue ds le cristal par un photon
     HLtot = 0.0D0

     !a Initialisation de l'attenuation a l'interieur de la particule: #############
     DATT = 1.0D+00

     !a Choix aleatoire de l'orientation du faisceau incident: #####################
     !  Ce choix etait mal fait, maintenant utilisation de la densite de proba
     !  f(theta) = sin(theta) / 2

     !<CHG,28/09/07>
     !      call random_number( rnETA )
     !      ETA = ( rnETA * 90.0D+00 ) * PI180
     !</CHG>

     call random_number( rnPHI )
     PHI = ( rnPHI * 60.0D+00 - 30.0D+00 ) * PI180

     !<CHG,28/09/07>
     !      call random_number( rnZETA )
     !      ZETA = ( rnZETA * 90.0D+00 ) * PI180
     !</CHG>

     CALL RANDOM_NUMBER(rnmu)
     mu = rnmu                            ! parceque l on veut theta
     ! entre 0 et 90 degres

     !a Section geometrique totale: ################################################

     auxSURFAPP = surfapp(RR,LL,ETA,ZETA,PHI,R32)

     sumSURFAPP = sumSURFAPP + auxSURFAPP

     !a Cosinus directeurs du faisceau incident: ###################################

     COSPHI   = DCOS( PHI )
     SINPHI   = DSIN( PHI )

     !******************************************************************************
     !  Je ne calcul plus aleatoirement d'angle ETA n'y d'angle ZETA, mais je dois, 
     !  parceque le code a ete ecrit a partir de ces angles, les recalculer a partir
     !  de THETA et PHI...
     !
     !      COSETA   = dcos( ETA )
     !      SINTHETA = dcos( ZETA ) * COSETA
     !      COSTHETA = dsqrt( 1.0D+00 - SINTHETA * SINTHETA )

     THETA = PI / 2.0D0 - DACOS(mu)     !parceque theta est defini par rapport
     !a OXY et non par rapport a OZ
     COSTHETA = DCOS(THETA)
     IF (COSTHETA .GT.  1.0D0) COSTHETA =  1.0D0
     IF (COSTHETA .LT. -1.0D0) COSTHETA = -1.0D0
     SINTHETA = DSQRT(1.0D0 - COSTHETA * COSTHETA)

     !--- Calcul de ETA et ZETA a partir de theta et phi
     ! ETA
     SINETA = COSTHETA * COSPHI
     IF (SINETA .GT.  1.0D0) SINETA =  1.0D0
     IF (SINETA .LT. -1.0D0) SINETA = -1.0D0
     COSETA = DSQRT(1.0D0 - SINETA * SINETA)
     ETA = DACOS(COSETA)
     ! ZETA
     ! si COSETA = 0 => faisceau incident suivant OX
     ! et dans ce cas ZETA = 0
     IF (COSETA .GT. 1.0D-4) THEN
        COSZETA = SINTHETA / COSETA
        ZETA = DACOS(COSZETA)
     ELSE
        ZETA = 0.0D0
     ENDIF

     ALPHAi = - COSTHETA * COSPHI
     BBETAi = - COSTHETA * SINPHI
     GAMMAi = - SINTHETA

     !a Impact du photon incident sur une face externe: ############################

     FACEimp = -999

     CALL IMPACT(LL2,AX,COSPHI,COSTHETA,SINPHI,SINTHETA,EPS12,SE,R32,         &
          Ximp,Yimp,Zimp,FACEimp)

     !a Verification des resultats d'IMPACT, concernant la face atteinte: ##########

     CALL TEST_IMPACT(NPHDIF,FACEimp,ETA,ZETA,PHI,PI180,EPS12)

     !a APARi et APERi: composantes parallele et perpendiculaire
     !a                 de l'amplitude incidente; polarisation aleatoire... ########

     CALL random_number( rnPOLAR )

     ReAPARi = ATOTi * DCOS( 2.0D+00*PI*rnPOLAR )
     ReAPERi = ATOTi * DSIN( 2.0D+00*PI*rnPOLAR )

     !a IPARi et IPERi = composantes parallele et perpendiculaire
     !a                  de l'intensite incidente... ###############################

     IPARi = ReAPARi*ReAPARi
     IPERi = ReAPERi*ReAPERi

     !a FINC = flux d'energie incidente (intensite totale incidente)... ############

     FINC = FINC + IPARi + IPERi

     !a Cosinus directeurs des faisceaux reflechi et transmis: #####################

     !a Coefficients de reflexion et de transmission de Fresnel: ###################

     CALL CHANG_AG(RN,RM,EPS12,                                               &
          ALPHAf,BBETAf,GAMMAf,FACEimp,                              &
          ALPHAi,BBETAi,GAMMAi,COSTHETAi,                            &
          ALPHAr,BBETAr,GAMMAr,                                      &
          ALPHAt,BBETAt,GAMMAt,COSTHETAt,                            &
          GRANDN,GRANDK,COSANGLEHETERO,                              &
          ReCRPAR,ReCRPER,ReCTPAR,ReCTPER,                           &
          ImCRPAR,ImCRPER,ImCTPAR,ImCTPER)

     !a Matrice d'amplitude diffusee, par reflexion sur la face externe: ###########

     ReCSRD(1,1) = ReCRPAR
     ImCSRD(1,1) = ImCRPAR

     ReCSRD(1,2) = 0.0D+00
     ImCSRD(1,2) = 0.0D+00

     ReCSRD(2,1) = 0.0D+00
     ImCSRD(2,1) = 0.0D+00

     ReCSRD(2,2) = ReCRPER
     ImCSRD(2,2) = ImCRPER

     !a APARd et APERd: composantes parallele et perpendiculaire
     !a                 de l'amplitude diffusee... #################################

     ReAPARd = ReCSRD(1,1)*ReAPARi
     ReAPERd = ReCSRD(2,2)*ReAPERi
     ImAPARd = ImCSRD(1,1)*ReAPARi
     ImAPERd = ImCSRD(2,2)*ReAPERi

     !a IPARd et IPERd = composantes parallele et perpendiculaire
     !a                  de l'intensite diffusee... ################################

     IPARd = ReAPARd*ReAPARd + ImAPARd*ImAPARd
     IPERd = ReAPERd*ReAPERd + ImAPERd*ImAPERd

     !a FDRT = flux d'energie diffusee par Ray-Tracing... ##########################

     FDRT = FDRT + IPARd + IPERd

     !a Angles de diffusion entre les faisceaux incident et diffuse': ##############

     COSGRANDTHETA = ALPHAi * ALPHAr +                                        &
          BBETAi * BBETAr +                                        &
          GAMMAi * GAMMAr
     !<CHG,27/09/2007>
     IF (COSGRANDTHETA .GT. 1.0D0)  COSGRANDTHETA = 1.0D0
     IF (COSGRANDTHETA .LT. -1.0D0) COSGRANDTHETA = -1.0D0
     !</CHG,27/09/2007>
     GRANDTHETA = dacos( COSGRANDTHETA ) / PI180

     if ( GRANDTHETA .LT. ( RESOLGT / 2.0D+00 ) ) then
        IDGT = 1
     else
        if( GRANDTHETA .GT. ( 180.0D+00 - RESOLGT / 2.0D+00 ) )then
           IDGT = NDGT
        else
           IDGT = idint( ( GRANDTHETA +                                     &
                RESOLGT / 2.0D+00 ) / RESOLGT ) + 1
        endif
     endif

     !a Matrice de diffusion brute (= non-normalisee) aux44GTGP, calculee a partir #
     !a des amplitudes diffusees par ray-tracing, AUX22GTGP: #######################

     call STOKES(ReCSRD,ImCSRD,aux44GTGP)

     aux44GT(:,:,IDGT) = aux44GT(:,:,IDGT) + aux44GTGP(:,:)

     !a Le faisceau transmis (refracte') vers l'interieur de la particule devient ##
     !a celui qui emerge, pour un parcours interne, a partir de la face interne: ###

     ALPHAe  = ALPHAt
     BBETAe  = BBETAt
     GAMMAe  = GAMMAt

     !a Cosinus directeurs de la normale (XNI,YNI,ZNI) et de
     !a l'orthonormale (XTI,YTI,ZTI) au plan d'incidence du dioptre atteint: ######

     call NORMALEPI(ALPHAi,BBETAi,GAMMAi,                                     &
          EPS12,FACEimp,PI,ALPHAf,BBETAf,GAMMAf,                    & 
          XNI,XTI,YNI,YTI,ZNI,ZTI)

     !a Cosinus directeurs de la normale et
     !a de l'orthonormale au plan de refraction: ###################################
     !  Rq: Normalement (XNI,YNI,ZNI) et (XNR,YNR,ZNR) sont identique, car rayon 
     !      incident, reflechi et refracte st ds le meme plan (1ere loi de snell 
     !      descarte). Ca a ete verifie et c'est OK! (voir PHMcode24)

     call NORMALE(ALPHAe,BBETAe,GAMMAe,XNI,YNI,ZNI,                           &
          XNR,XTR,YNR,YTR,ZNR,ZTR)

     !a Matrice d'amplitudes transmises vers l'interieur: ##########################

     auxCS = dsqrt(GRANDN*COSTHETAt/COSTHETAi)

     ReCS(1,1) = auxCS*ReCTPAR
     ImCS(1,1) = auxCS*ImCTPAR

     ReCS(1,2) = 0.0D+00
     ImCS(1,2) = 0.0D+00

     ReCS(2,1) = 0.0D+00
     ImCS(2,1) = 0.0D+00

     ReCS(2,2) = auxCS*ReCTPER
     ImCS(2,2) = auxCS*ImCTPER


     !******************************************************************************
     !             IHM Part, prise en compte d'inclusion avant d'atteindre
     !             la face oppposee!!!
     !
     !     1) Je dois qd meme d'abord calculer la distance parcourue (DELTA3) avant 
     !        d'atteindre la face opposee (FACEimp), c'est le calcul ci-dessous!
     !
     !     2) Ensuite calcul du rayon de l'inclusion afin de determiner si DELTA3
     !        est suffisante pour que l'on puisse rester ds l'approximation d'onde
     !        plane (pour la theorie de Mie), voir article Macke et Michshenko, 
     !        JGR, 96
     !
     !     3) Calcul egalement de la distance avant de rencontrer une inclusion
     !
     !     4)
     !
     !     5)
     !
     !     6)
     !******************************************************************************

     !******************************************************************************
     !a 1) Point d'impact, et distance apres un parcours interne: ##################
     !******************************************************************************
     !--- Flag_face = flag qui indique si le photon est sur une face du cristal
     !    Flag_face = 1 -> photon sur la face FACEimp
     !    Flag_face = 0 -> photon ds le cristal
     Flag_face = 1

2020 continue

     auxXimp = Ximp
     auxYimp = Yimp
     auxZimp = Zimp

     do IFACE = 0, 5
        auxauxDELTA3 = ALPHAe * ALPHAf(IFACE) + BBETAe * BBETAf(IFACE)
        if( dabs( auxauxDELTA3 ).gt.EPS12 )then
           auxDELTA3(IFACE) = ( AX -                                     &
                auxXimp * ALPHAf(IFACE) -                &
                auxYimp * BBETAf(IFACE) ) / auxauxDELTA3
        else
           auxDELTA3(IFACE) = -DDDD
        endif
     enddo

     if( dabs( GAMMAe ).gt.EPS12 )then
        auxDELTA3(6) =  ( LL2 - auxZimp ) / GAMMAe
        auxDELTA3(7) = -( LL2 + auxZimp ) / GAMMAe
     else
        auxDELTA3(6) = -DDDD
        auxDELTA3(7) = -DDDD
     endif

     DELTA3 = +DDDD
     IF (Flag_face .EQ. 1) THEN
        auxDELTA3(FACEimp) = -DDDD
     ENDIF

     do IFACE = 0, 7
        if( auxDELTA3(IFACE).gt.EPS12 .and. auxDELTA3(IFACE).lt.DELTA3 )then
           DELTA3 = auxDELTA3(IFACE)
           auxFACEimp = IFACE
        endif
     enddo

     !--- Rq: les coordonnees du point d'impact ainsi que la face d'impact sont 
     !        determines plus tard, cad lorsque l'on sait que le photon ne 
     !        rencontre pas d'inclusion

     !******************************************************************************
     ! 2) Choix du type d'inclusion, et calcul de son rayon grace aux courbe de 
     !    distribution en taille cumulees 
     !******************************************************************************
     !--- tirage d'un nombre aleatoire que l'on compare a Ksuie (pourcentage suie 
     !    bulle d'air) afin de savoir si l'inclusion est de la suie ou de l'air.
     !    ISB = 1 -> air
     !    ISB = 2 -> suie
     CALL RANDOM_NUMBER(rnk)
     IF (rnk .GE. Ksuie) THEN
        ISB = 1
        DRNI = RNB
        DINI = RMB
        rmin = rmin_b
        rmax = rmax_b
        NNr(:)=NNr_b(:)
     ELSE
        ISB = 2
        DRNI = RNS
        DINI = RMS     
        rmin = rmin_s
        rmax = rmax_s
        NNr(:) = NNr_s(:)
     ENDIF

2100 CONTINUE

     !--- Calcul du rayon de l'inclusion (r_inc) grace a NNr et tab_r
     CALL incl_jradius(numr,NNr,jr)
     r_inc = tab_r(ISB,jr)

     !******************************************************************************
     ! 3) Calcul de la distance (HL) parcourue avant d'atteindre une inclusion (MC)
     !    a partir du libre parcour moyen l
     !******************************************************************************
     CALL RANDOM_NUMBER(rnl)
     HL = -l * DLOG(rnl)

     !******************************************************************************
     !   Does the photon hit an inclusion??
     !******************************************************************************
     IF (HL .GE. (DELTA3-4.0D0*r_inc)) GOTO 2200 ! photon tape la face opposee
     IF (DELTA3 .LE. 8.0D0*r_inc) GOTO 2200 ! distance a la face opposee
     ! trop courte pour que l'inclusion 
     ! tiree soit sur le trajet
     IF (HL .LE. 4.0D0*r_inc) GOTO 2100 ! distance a l'inclusion trop courte
     !      IF ((DELTA3-HL) .LE. 4.0D0*r_inc) GOTO 2100 !distance entre face opposee
     !                                                  !et inclusion trop courte


     !******************************************************************************
     !   YES the photon hit an inclusion!!!
     !******************************************************************************
     !--- Count it
     cont = cont + 1
     !--- Coordonnees de l'inclusion qui devient le nouveau point de depart
     Ximp = auxXimp + HL * ALPHAe
     Yimp = auxYimp + HL * BBETAe
     Zimp = auxZimp + HL * GAMMAe

     !--- Verification que l'on est tjs ds le cristal
     IF (DABS(Ximp) .GT. AX) THEN
        WRITE(6,*)'Pb inclusion endehors du cristal!'
        WRITE(6,*)'Ximp =',Ximp,' alors que Ax=',AX
        WRITE(6,*)'HL=',HL,' DELTA=',DELTA3
        WRITE(6,*)'auxX:',auxXimp,auxYimp,auxZimp
        WRITE(6,*)'Alphae:',ALPHAe,BBETAe,GAMMAe
        WRITE(6,*)'Alphae_p:',ALPHAe_prev,BBETAe_prev,GAMMAe_prev
        WRITE(6,*)'auxDELTA3=',auxDELTA3(:)
        WRITE(6,*)'flag_face=',Flag_face
        STOP
     ENDIF
     IF (DABS(Yimp) .GT. RR) THEN
        WRITE(6,*)'Pb inclusion endehors du cristal!'
        WRITE(6,*)'Yimp =',Yimp,' alors que RR=',RR
        WRITE(6,*)'HL=',HL,' DELTA=',DELTA3
        WRITE(6,*)'auxX:',auxXimp,auxYimp,auxZimp
        WRITE(6,*)'Alphae:',ALPHAe,BBETAe,GAMMAe
        WRITE(6,*)'Alphae_p:',ALPHAe_prev,BBETAe_prev,GAMMAe_prev
        WRITE(6,*)'auxDELTA3=',auxDELTA3(:)
        WRITE(6,*)'flag_face=',Flag_face
        STOP
     ENDIF
     IF (DABS(Zimp) .GT. LL2) THEN
        WRITE(6,*)'Pb inclusion endehors du cristal!'
        WRITE(6,*)'Zimp =',Zimp,' alors que LL2=',LL2
        WRITE(6,*)'HL=',HL,' DELTA=',DELTA3
        WRITE(6,*)'auxX:',auxXimp,auxYimp,auxZimp
        WRITE(6,*)'Alphae:',ALPHAe,BBETAe,GAMMAe
        WRITE(6,*)'Alphae_p:',ALPHAe_prev,BBETAe_prev,GAMMAe_prev
        WRITE(6,*)'auxDELTA3=',auxDELTA3(:)
        WRITE(6,*)'flag_face=',Flag_face
        STOP
     ENDIF


     !--- Distance parcourue ds le cristal
     HLtot = HLtot + HL

     !a Mise a jour de l'attenuation: ##############################################

     DATT = DATT * dexp( - HL * KKKK * GRANDK * COSANGLEHETERO )

     !******************************************************************************
     !   Compute the phase function, coefficient Sij of the phase matrix, and
     !   the single scattering albedo of the inclusion
     !******************************************************************************
     !--- Before look if those optical properties are not already in tables
     IF (PIZ(ISB,jr) .GT. 900.0D0) THEN ! Need to compute it or read it from
        ! library
        CALL opt_properties(libdir, ISB,jr,numr,nmu,nmumax,lambda,r_inc,RMU,PMU,DRNI,  &
             DINI,PIZ,PTETA,DRES1,DRES2,DIMS1,DIMS2)
        !      ENDIF

        !******************************************************************************
        !   Now compute the cumulative phase function (c'est une fonction de       
        !   partition!!!!)
        !******************************************************************************
        !--- Look before if it does not already exist
        !      IF (Pcum(ISB,jr,floor(nmu/2.)) .GT. 900.0D0) THEN !Compute it
        CALL P11_cumul(nmu,nmumax,RMU,PMU,PTETA(ISB,jr,:),Pcum(ISB,jr,:))
     ENDIF

     !******************************************************************************
     !   Absorption by the inclusion, DATT*(1-PIZ) is absorbed and the rest is 
     !   scattered...
     !******************************************************************************
     DATT = DATT * PIZ(ISB,jr)

     !******************************************************************************
     !   Evaluation des angles (thetad, phid) entre le photon diffusee par 
     !   l'inclusion et le photon incident, grace a la fct de phase cumulee (Pcum),
     !   ainsi que des champs incidents sur l'inclusion et de la matrice de phase
     !   de l'inclusion
     !   Rq: Ces angles sont calcules ds le repere lie au photon incident, cad 
     !       ds le repere du plan de reference avant diffusion (Oxyz).
     !******************************************************************************
     CALL Angle_diff(nmu,nmumax,RMU,Pcum(ISB,jr,:),DRES1(ISB,jr,:),           &
          DRES2(ISB,jr,:),DIMS1(ISB,jr,:),DIMS2(ISB,jr,:),ReCS,    &
          ImCS,ReAPARi,ReAPERi,jt,thetad,phid)
     !<TMP>
     IF (phid .NE. phid) THEN
        !      IF (IEEE_IS_NAN(phid)) THEN
        !         WRITE(6,*)IEEE_IS_NAN(phid)
        WRITE(6,*)'Apres Angle_diff!'
        WRITE(6,*)'phid=',phid
        WRITE(6,*)'thetad=',thetad
        WRITE(6,*)'jt=',jt
        WRITE(6,*)ReCS(:,:)
        WRITE(6,*)ImCS(:,:)
        WRITE(6,*)ReCS_rot(:,:)
        WRITE(6,*)DRES1(ISB,jr,:)
        WRITE(6,*)DRES2(ISB,jr,:)
        WRITE(6,*)DIMS1(ISB,jr,:)
        WRITE(6,*)DIMS2(ISB,jr,:)
        WRITE(6,*)COSAEi,SINAEi
        WRITE(6,*)'Alphaf:',ALPHAt,BBETAt,GAMMAt
        WRITE(6,*)'Alphae:',ALPHAe,BBETAe,GAMMAe
        WRITE(6,*)'Alphae_p:',ALPHAe_prev,BBETAe_prev,GAMMAe_prev
        pause
     ENDIF
     !</TMP>

     !******************************************************************************
     !         Compute the new cosine director of the scattered photon
     !         (ALPHAn,BBETAn,GAMMAn) in the base linked to the cristal (OXYZ)
     !
     !   Rq: les angles tetad et phid trouves ne sont pas les 
     !       angles d'euler du photon diffuse ds le repere (OXYZ), mais les
     !       angles du photon diffuse par rapport au repere (0xyz) lie 
     !       au photon incident, cad repere lie au plan de reference avant 
     !       diffusion.
     !******************************************************************************
     !<CHG,02/10/07>
     !      CALL New_direction(thetad,phid,XTR,YTR,ZTR,XNR,YNR,ZNR,ALPHAe,BBETAe,    &
     !                         GAMMAe,ALPHAn,BBETAn,GAMMAn)
     !</CHG,02/10/07>

     ! I am using the rotational matrix define in Ramella-Roman et al., 
     ! 2005a
     ! (a) rotation of v around w of phid
     !     Rq: v = Normal to the reference plan (old scattering plan)
     !         w = cosine director of the incident photon
     !         v' = Reuler(w,phid)*v
     !         w' = w
     ! (b) rotation of w around v' of Thetad
     !         w" = Reuler(v',Thetad)*w'
     !         v" = v'
     !-- (a)
     v_old(1) = XNR
     v_old(2) = YNR
     v_old(3) = ZNR
     w(1)     = ALPHAe
     w(2)     = BBETAe
     w(3)     = GAMMAe
     CALL rot_uvw(v_old,w,phid,v)
     !-- (b)
     w_old(:) = w(:)
     CALL rot_uvw(w_old,v,thetad,w)
     ALPHAn = w(1)
     BBETAn = w(2)
     GAMMAn = w(3)
     XNRn   = v(1)
     YNRn   = v(2)
     ZNRn   = v(3)
     XTRn   = v(2)*w(3) - v(3)*w(2)
     YTRn   = v(3)*w(1) - v(1)*w(3)
     ZTRn   = v(1)*w(2) - v(2)*w(1)

     !******************************************************************************
     !     (1) Compute the angle needed to rotate the incident electric field in 
     !         the new scattering plan.
     !         Actually I compute the sinus (SINAEi) and cosinus (COSAEi) of this 
     !         angle.
     !         Rq: Si je ne me suis pas trompe cet angle correspond uniquement 
     !             a phid!!!!!!!!!!!!!!!!!!! A VERIFIER !!!!!!!!!!!!!!!!!!!!!
     !             C'EST VERIFIE, ET L'ANGLE CORRESPOND BIEN A PHID!!!!!!!!!!
     !         
     !     (2) Compute also the new direction of the normal and orthonormal to the 
     !         new incident plan which is actually the scattering plan....
     !         (XNRn,YNRn,ZNRn) and (XTRn,YTRn,ZTRn)
     !
     !     CHG, 2/10/07
     !     (1) n'a plus de raison d'etre car car cet angle correspond a phid
     !     (2) est fait ci-avant grace a la subroutine de rotation rot_uvw 
     !         provenant de Ramella, 2005 
     !******************************************************************************
     !<CHG,02/10/07>
     !      CALL Angle_rotE(thetad,ALPHAe,BBETAe,GAMMAe,ALPHAn,BBETAn,GAMMAn,XNR,    &
     !                      YNR,ZNR,XTR,YTR,ZTR,SINAEi,COSAEi,XNRn,YNRn,ZNRn,XTRn,   &
     !                      YTRn,ZTRn)
     !</CHG,02/10/07>
     COSAEi = DCOS(phid)
     SINAEi = DSIN(phid)

     !--- Replace the old direction by the new one!
     ALPHAe = ALPHAn
     BBETAe = BBETAn
     GAMMAe = GAMMAn

     !--- Replace the old normal and orthonormal by the new one's
     XNR = XNRn
     YNR = YNRn
     ZNR = ZNRn
     XTR = XTRn
     YTR = YTRn
     ZTR = ZTRn

     !******************************************************************************
     !    Rotation of the electric field in the scattering plan and multiply it 
     !    by the phase matrix of the inclusion.
     !    (a) normalize the electric field
     !    (b) rotate it
     !    (c) multiply by the inclusion phase matrix
     !    (d) normalize it again and multiply by the incident norm
     !        and stoke the new phase matrix in the old one
     !******************************************************************************
     !--- (a)
     CALL Field_nrj(ReCS,ImCS,ReAPARi,ReAPERi,NRJ_Ei)
     IF (NRJ_Ei*DATT*DATT .LT. EPS8) THEN ! NRJ trop faible on ne doit plus 
        ! suivre le photon
        WRITE(6,*)'Attention NRJ faible!'
        WRITE(6,*)'NRJ_Ei=',NRJ_Ei*DATT*DATT
        !        GOTO 2900
     ENDIF
     ReCS(:,:) = ReCS(:,:) / DSQRT(NRJ_Ei)
     ImCS(:,:) = ImCS(:,:) / DSQRT(NRJ_Ei)
     !--- (b)
     Matrot(1,1) =  COSAEi
     Matrot(1,2) =  SINAEi
     Matrot(2,1) = -SINAEi
     Matrot(2,2) =  COSAEi
     CALL DMRMUL(Matrot,ReCS,ReCS_rot,2,2,2)
     CALL DMRMUL(Matrot,ImCS,ImCS_rot,2,2,2)

     !--- (c)
     ReCS(1,:) = ReCS_rot(1,:)*DRES2(ISB,jr,jt) -                             &
          ImCS_rot(1,:)*DIMS2(ISB,jr,jt)
     ReCS(2,:) = ReCS_rot(2,:)*DRES1(ISB,jr,jt) -                             &
          ImCS_rot(2,:)*DIMS1(ISB,jr,jt)
     ImCS(1,:) = ImCS_rot(1,:)*DRES2(ISB,jr,jt) +                             &
          ReCS_rot(1,:)*DIMS2(ISB,jr,jt)
     ImCS(2,:) = ImCS_rot(2,:)*DRES1(ISB,jr,jt) +                             &
          ReCS_rot(2,:)*DIMS1(ISB,jr,jt)

     !--- (d)
     CALL Field_nrj(ReCS,ImCS,ReAPARi,ReAPERi,NRJ_Ed)
     IF (NRJ_Ed*DATT*DATT .LT. EPS8) THEN ! NRJ trop faible on ne doit plus 
        ! suivre le photon
        WRITE(6,*)'Attention NRJ faible!'
        WRITE(6,*)'NRJ_Ed=',NRJ_Ed*DATT*DATT
        !        GOTO 2900
     ENDIF
     ReCS(:,:) = ReCS(:,:) / DSQRT(NRJ_Ed) * DSQRT(NRJ_Ei)
     ImCS(:,:) = ImCS(:,:) / DSQRT(NRJ_Ed) * DSQRT(NRJ_Ei)

     !******************************************************************************
     !   Continue to follow the photon, i.e. compute the distance to the cristal 
     !   border... etc, i.e. come back to 2020
     !******************************************************************************
     !--- Put the Flag_face to zero
     Flag_face = 0

     GOTO 2020

2200 CONTINUE

     Flag_face = 1

     !--- Determination du point d'impact ainsi que de la face d'impact
     Ximp = auxXimp + DELTA3 * ALPHAe
     Yimp = auxYimp + DELTA3 * BBETAe
     Zimp = auxZimp + DELTA3 * GAMMAe
     FACEimp = auxFACEimp

     !a Mise a jour de l'attenuation: ##############################################
     DATT = DATT * dexp( - DELTA3 * KKKK * GRANDK * COSANGLEHETERO )

     !a Cosinus directeurs de la normale et de l'orthonormale au plan d'incidence de
     !a ce nouveau dioptre atteint: ################################################

     call NORMALEM(ALPHAe,BBETAe,EPS12,GAMMAe,          &
          FACEimp,ALPHAf,BBETAf,GAMMAf,        &
          XNR,YNR,ZNR,XND,XTD,YND,YTD,ZND,ZTD)

     !a Angle de rotation A de la normale interne pour l'amener dans le repere du ##
     !a nouveau plan d'incidence atteint, Nr -> Nd: ################################

     DCOSA = XNR * XND + YNR * YND + ZNR * ZND
     DSINA = XNR * XTD + YNR * YTD + ZNR * ZTD

     ReCSA(1,1) = DCOSA*ReCS(1,1) + DSINA*ReCS(2,1)
     ImCSA(1,1) = DCOSA*ImCS(1,1) + DSINA*ImCS(2,1)

     ReCSA(1,2) = DCOSA*ReCS(1,2) + DSINA*ReCS(2,2)
     ImCSA(1,2) = DCOSA*ImCS(1,2) + DSINA*ImCS(2,2)

     ReCSA(2,1) = -DSINA*ReCS(1,1) + DCOSA*ReCS(2,1)
     ImCSA(2,1) = -DSINA*ImCS(1,1) + DCOSA*ImCS(2,1)

     ReCSA(2,2) = -DSINA*ReCS(1,2) + DCOSA*ReCS(2,2)
     ImCSA(2,2) = -DSINA*ImCS(1,2) + DCOSA*ImCS(2,2)

     !a Cosinus directeurs des faisceaux reflechi et transmis (refracte'), apres un
     !a (ou plusieurs) trajet(s) interne(s) a la particule: ########################

     !a Coefficients de reflexion et de transmission de Fresnel: ###################

     call CHANG_GA(RN,RM,                                &
          ALPHAf,BBETAf,GAMMAf,FACEimp,         &
          ALPHAe,BBETAe,GAMMAe,COSTHETAi,       &
          ALPHAr,BBETAr,GAMMAr,                 &
          ALPHAt,BBETAt,GAMMAt,COSTHETAt,IREFL, &
          GRANDN,GRANDK,COSANGLEHETERO,         &
          ReCRPAR,ReCRPER,ReCTPAR,ReCTPER,      &
          ImCRPAR,ImCRPER,ImCTPAR,ImCTPER)

     !<TMP>
     ALPHAe_prev = ALPHAe
     BBETAe_prev = BBETAe
     GAMMAe_prev = GAMMAe
     !</TMP>

     !a Cas particulier de reflexion totale: #######################################

     if( IREFL.eq.1 ) goto 12

     !a Matrice d'amplitudes diffusees apres deux refractions: #####################

     auxCT = DATT*dsqrt(COSTHETAt/(GRANDN*COSTHETAi))
     !<TMP>
     IF (COSTHETAi .LT. 0.0D0) THEN
        WRITE(6,*)'IPHINC=',IPHINC
        WRITE(6,*)'rnphi=',rnphi
        WRITE(6,*)'rnmu=',rnmu
        WRITE(6,*)'COSt=',COSTHETAt
        WRITE(6,*)'GRANDN=',GRANDN
        WRITE(6,*)'FACEimp=',FACEimp
        WRITE(6,*)'Flag_face=',Flag_face
        WRITE(6,*)'HL=',HL,' DELTA=',DELTA3
        WRITE(6,*)'Ximp:',Ximp,Yimp,Zimp
        WRITE(6,*)'ALPHAf:',ALPHAf(FACEimp),BBETAf(FACEimp),GAMMAf(FACEimp)
        WRITE(6,*)'Alphae_n:',ALPHAr,BBETAr,GAMMAr
        WRITE(6,*)'Alphae_p:',ALPHAe_prev,BBETAe_prev,GAMMAe_prev
        WRITE(6,*)'COSi=',COSTHETAi
        STOP
     ENDIF
     !</TMP>

     ReCT(1,1) = auxCT*ReCTPAR
     ImCT(1,1) = auxCT*ImCTPAR

     ReCT(1,2) = 0.0D+00
     ImCT(1,2) = 0.0D+00

     ReCT(2,1) = 0.0D+00
     ImCT(2,1) = 0.0D+00

     ReCT(2,2) = auxCT*ReCTPER
     ImCT(2,2) = auxCT*ImCTPER

     ReCS(1,1) = ReCT(1,1)*ReCSA(1,1) + ReCT(1,2)*ReCSA(2,1) - &
          ImCT(1,1)*ImCSA(1,1) - ImCT(1,2)*ImCSA(2,1)
     ImCS(1,1) = ReCT(1,1)*ImCSA(1,1) + ImCT(1,1)*ReCSA(1,1) + &
          ReCT(1,2)*ImCSA(2,1) + ImCT(1,2)*ReCSA(2,1)

     ReCS(1,2) = ReCT(1,1)*ReCSA(1,2) + ReCT(1,2)*ReCSA(2,2) - &
          ImCT(1,1)*ImCSA(1,2) - ImCT(1,2)*ImCSA(2,2)
     ImCS(1,2) = ReCT(1,1)*ImCSA(1,2) + ImCT(1,1)*ReCSA(1,2) + &
          ReCT(1,2)*ImCSA(2,2) + ImCT(1,2)*ReCSA(2,2)

     ReCS(2,1) = ReCT(2,1)*ReCSA(1,1) + ReCT(2,2)*ReCSA(2,1) - &
          ImCT(2,1)*ImCSA(1,1) - ImCT(2,2)*ImCSA(2,1)
     ImCS(2,1) = ReCT(2,1)*ImCSA(1,1) + ImCT(2,1)*ReCSA(1,1) + &
          ReCT(2,2)*ImCSA(2,1) + ImCT(2,2)*ReCSA(2,1)

     ReCS(2,2) = ReCT(2,1)*ReCSA(1,2) + ReCT(2,2)*ReCSA(2,2) - &
          ImCT(2,1)*ImCSA(1,2) - ImCT(2,2)*ImCSA(2,2)
     ImCS(2,2) = ReCT(2,1)*ImCSA(1,2) + ImCT(2,1)*ReCSA(1,2) + &
          ReCT(2,2)*ImCSA(2,2) + ImCT(2,2)*ReCSA(2,2)

     !a Directions de polarisation du champ emergent: ##############################

     call NORMALE(ALPHAt,BBETAt,GAMMAt,XND,YND,ZND, &
          XNS,XTS,YNS,YTS,ZNS,ZTS)

     !a Cosinus directeurs de la normale au plan de diffusion: #####################

     call RCDIRPD(ALPHAi,BBETAi,GAMMAi,                   & 
          ALPHAt,BBETAt,GAMMAt,EPS12,XNI,YNI,ZNI, &
          XNP,YNP,ZNP)

     !a Angle de rotation B de la normale au plan emergent pour l'amener dans le ###
     !a reperer du plan de diffusion, Ns -> Np: ####################################

     DCOSB = XNP * XNS + YNP * YNS + ZNP * ZNS
     DSINB = XNP * XTS + YNP * YTS + ZNP * ZTS

     ReCSD(1,1) = DCOSB*ReCS(1,1) - DSINB*ReCS(2,1)
     ImCSD(1,1) = DCOSB*ImCS(1,1) - DSINB*ImCS(2,1)

     ReCSD(1,2) = DCOSB*ReCS(1,2) - DSINB*ReCS(2,2)
     ImCSD(1,2) = DCOSB*ImCS(1,2) - DSINB*ImCS(2,2)

     ReCSD(2,1) = DSINB*ReCS(1,1) + DCOSB*ReCS(2,1)
     ImCSD(2,1) = DSINB*ImCS(1,1) + DCOSB*ImCS(2,1)

     ReCSD(2,2) = DSINB*ReCS(1,2) + DCOSB*ReCS(2,2)
     ImCSD(2,2) = DSINB*ImCS(1,2) + DCOSB*ImCS(2,2)

     !a Angle de rotation D de la normale Np au plan de diffusion et la normale Ni #
     !a au plan d'incidence: #######################################################

     DCOSD = XNP * XNI + YNP * YNI + ZNP * ZNI
     DSIND = XNP * XTI + YNP * YTI + ZNP * ZTI

     ReCSRD(1,1) = DCOSD*ReCSD(1,1) - DSIND*ReCSD(1,2)
     ImCSRD(1,1) = DCOSD*ImCSD(1,1) - DSIND*ImCSD(1,2)

     ReCSRD(1,2) = DSIND*ReCSD(1,1) + DCOSD*ReCSD(1,2)
     ImCSRD(1,2) = DSIND*ImCSD(1,1) + DCOSD*ImCSD(1,2)

     ReCSRD(2,1) = DCOSD*ReCSD(2,1) - DSIND*ReCSD(2,2)
     ImCSRD(2,1) = DCOSD*ImCSD(2,1) - DSIND*ImCSD(2,2)

     ReCSRD(2,2) = DSIND*ReCSD(2,1) + DCOSD*ReCSD(2,2)
     ImCSRD(2,2) = DSIND*ImCSD(2,1) + DCOSD*ImCSD(2,2)

     !a APARd et APERd: composantes parallele et perpendiculaire
     !a                 de l'amplitude diffusee... #################################

     ReAPARd = ReCSRD(1,1)*ReAPARi + ReCSRD(1,2)*ReAPERi
     ReAPERd = ReCSRD(2,1)*ReAPARi + ReCSRD(2,2)*ReAPERi
     ImAPARd = ImCSRD(1,1)*ReAPARi + ImCSRD(1,2)*ReAPERi
     ImAPERd = ImCSRD(2,1)*ReAPARi + ImCSRD(2,2)*ReAPERi

     !a IPARd et IPERd = composantes parallele et perpendiculaire
     !a                  de l'intensite diffusee... ################################

     IPARd = ReAPARd*ReAPARd + ImAPARd*ImAPARd
     IPERd = ReAPERd*ReAPERd + ImAPERd*ImAPERd

     !a FDRT = flux d'energie diffusee par Ray-Tracing... ##########################

     FDRT = FDRT + IPARd + IPERd

     !a Angle de diffusion entre les faisceaux incident et diffuse': ###############

     COSGRANDTHETA = ALPHAi * ALPHAt + &
          BBETAi * BBETAt + &
          GAMMAi * GAMMAt
     IF (COSGRANDTHETA .GT.  1.0D0) COSGRANDTHETA =  1.0D0
     IF (COSGRANDTHETA .LT. -1.0D0) COSGRANDTHETA = -1.0D0
     GRANDTHETA = dacos( COSGRANDTHETA ) / PI180

     if( GRANDTHETA.lt.( RESOLGT / 2.0D+00 ) )then
        IDGT = 1
     else
        if( GRANDTHETA.gt.( 180.0D+00 - RESOLGT / 2.0D+00 ) )then
           IDGT = NDGT
        else
           IDGT = idint( ( GRANDTHETA +                        &
                RESOLGT / 2.0D+00 ) / RESOLGT ) + 1
        endif
     endif

     !a Matrice de diffusion brute (= non-normalisee) aux44GTGP, calculee a partir #
     !a des amplitudes diffusees par ray-tracing, AUX22GTGP: #######################

     call STOKES(ReCSRD,ImCSRD,aux44GTGP)

     !<TMP>
     IF (IDGT .LT. 1) THEN
        WRITE(6,*)'IPHINC=',IPHINC
        WRITE(6,*)'rnphi=',rnphi
        WRITE(6,*)'rnmu=',rnmu
        WRITE(6,*)'Flag_face=',Flag_face
        WRITE(6,*)'IDGT2=',IDGT
        WRITE(6,*)'Grandtheta=',GRANDTHETA
        WRITE(6,*)'cosGrandtheta=',COSGRANDTHETA
        WRITE(6,*)ALPHAi,BBETAi,GAMMAi
        WRITE(6,*)ALPHAt,BBETAt,GAMMAt
        WRITE(6,*)ALPHAe,BBETAe,GAMMAe
        WRITE(6,*)'thetad, phid =',thetad,phid
        WRITE(6,*)XNR,YNR,ZNR
        WRITE(6,*)XTR,YTR,ZTR
        WRITE(6,*)ReCS(:,:)
        WRITE(6,*)ImCS(:,:)
        WRITE(6,*)ReCS_rot(:,:)
        WRITE(6,*)ISB,jr,jt
        WRITE(6,*)DRES2(ISB,jr,jt),DIMS2(ISB,jr,jt)
        WRITE(6,*)DRES1(ISB,jr,jt),DIMS1(ISB,jr,jt)
        WRITE(6,*)'cosai=',COSAEi,SINAEi
        WRITE(6,*)'datt=',DATT
        pause
     ENDIF
     !</TMP>

     do JP = 1, 4
        do KP = 1, 4
           aux44GT(JP,KP,IDGT) = aux44GT(JP,KP,IDGT) + aux44GTGP(JP,KP)
        enddo
     enddo

     !a Les etapes ci-dessus ne sont pas prises en compte lors du cas particulier de
     !a reflexion totale par une face interne: #####################################

12   continue

     !a Le faisceau reflechi a l'interieur de la particule devient le faisceau emer-
     !a gent a partir de cette face interne: #######################################

     ALPHAe  = ALPHAr
     BBETAe  = BBETAr
     GAMMAe  = GAMMAr

     !a Directions de la normale et de l'orthonormale au plan reflechi: ############

     call NORMALE(ALPHAe,BBETAe,GAMMAe,XND,YND,ZND, &
          XNR,XTR,YNR,YTR,ZNR,ZTR)

     !a Matrice d'amplitudes reflechies vers l'interieur de la particule: ##########

     ReCR(1,1) = ReCRPAR
     ImCR(1,1) = ImCRPAR

     ReCR(1,2) = 0.0D+00
     ImCR(1,2) = 0.0D+00

     ReCR(2,1) = 0.0D+00
     ImCR(2,1) = 0.0D+00

     ReCR(2,2) = ReCRPER
     ImCR(2,2) = ImCRPER

     ReCS(1,1) = ReCR(1,1)*ReCSA(1,1) + ReCR(1,2)*ReCSA(2,1) - &
          ImCR(1,1)*ImCSA(1,1) - ImCR(1,2)*ImCSA(2,1)
     ImCS(1,1) = ReCR(1,1)*ImCSA(1,1) + ImCR(1,1)*ReCSA(1,1) + &
          ReCR(1,2)*ImCSA(2,1) + ImCR(1,2)*ReCSA(2,1)

     ReCS(1,2) = ReCR(1,1)*ReCSA(1,2) + ReCR(1,2)*ReCSA(2,2) - &
          ImCR(1,1)*ImCSA(1,2) - ImCR(1,2)*ImCSA(2,2)
     ImCS(1,2) = ReCR(1,1)*ImCSA(1,2) + ImCR(1,1)*ReCSA(1,2) + &
          ReCR(1,2)*ImCSA(2,2) + ImCR(1,2)*ReCSA(2,2)

     ReCS(2,1) = ReCR(2,1)*ReCSA(1,1) + ReCR(2,2)*ReCSA(2,1) - &
          ImCR(2,1)*ImCSA(1,1) - ImCR(2,2)*ImCSA(2,1)
     ImCS(2,1) = ReCR(2,1)*ImCSA(1,1) + ImCR(2,1)*ReCSA(1,1) + &
          ReCR(2,2)*ImCSA(2,1) + ImCR(2,2)*ReCSA(2,1)

     ReCS(2,2) = ReCR(2,1)*ReCSA(1,2) + ReCR(2,2)*ReCSA(2,2) - &
          ImCR(2,1)*ImCSA(1,2) - ImCR(2,2)*ImCSA(2,2)
     ImCS(2,2) = ReCR(2,1)*ImCSA(1,2) + ImCR(2,1)*ReCSA(1,2) + &
          ReCR(2,2)*ImCSA(2,2) + ImCR(2,2)*ReCSA(2,2)

     !a APARd et APERd: composantes parallele et perpendiculaire
     !a                 de l'amplitude diffusee... #################################

     ReAPARd = ReCS(1,1)*ReAPARi + ReCS(1,2)*ReAPERi
     ReAPERd = ReCS(2,1)*ReAPARi + ReCS(2,2)*ReAPERi
     ImAPARd = ImCS(1,1)*ReAPARi + ImCS(1,2)*ReAPERi
     ImAPERd = ImCS(2,1)*ReAPARi + ImCS(2,2)*ReAPERi

     !a IPARd et IPERd = composantes parallele et perpendiculaire
     !a                  de l'intensite diffusee... ################################

     IPARd = DATT*DATT*(ReAPARd * ReAPARd + ImAPARd * ImAPARd)
     IPERd = DATT*DATT*(ReAPERd * ReAPERd + ImAPERd * ImAPERd)

     !a Critere d'arret pour les parcours internes: ################################

     if( ( ( IPARd+IPERd )/( IPARi+IPERi ) ).gt.(EPS6) ) goto 2020

2900 CONTINUE

     !a Nombre des photons diffuses: ###############################################

     NPHDIF = NPHDIF + 1

     !a Repartition angulaire de l'energie diffractee: #############################

     do IDGT = 1, NDGTFR
        EDFR(IDGT) = EDFR(IDGT) +                          &
             ( minEDFR(IDGT) +                     &
             ( auxSURFAPP - minSURFAPP ) *       &
             ( maxEDFR(IDGT) - minEDFR(IDGT) ) / &
             ( maxSURFAPP - minSURFAPP ) )
        !aaaa    EDFR(IDGT) = 0.0D+00
     enddo

     !a La fin de la troisieme boucle: les photon incidents ########################

  enddo ! Fin de la boucle sur les photons incidents

  !a Flux d'energie diffractee: #################################################

  FDFR = 0.0D+00
  do IDGT = 1, NDGTFR
     FDFR = FDFR + ATOTi * ATOTi * EDFR(IDGT)
  enddo

  !a Section efficaces totales (GO=RT+FR) de diffusion et d'extinction: #########

  CEXTGO = ( FDFR / FINC + 1.0D+00     ) * sumSURFAPP / NPHDIF
  CSCAGO = ( FDFR / FINC + FDRT / FINC ) * sumSURFAPP / NPHDIF

  !a Albedo pour une diffusion, ray-tracing seulement: ##########################

  PIZERO = FDRT / FINC

  !a Albedo pour une diffusion, total (eq.15, Mishchenko & Macke 1998): #########

  PIZERO = ( PIZERO + 1.0D+00 ) / 2.0D+00

  !a La contribution due a la diffraction n'intervient que dans les elements dia-
  !a gonaux de la matrice de diffusion: #########################################

  do IDGT = 1, NDGTFR
     aux44GT(1,1,IDGT) = aux44GT(1,1,IDGT) + EDFR(IDGT)
     aux44GT(2,2,IDGT) = aux44GT(2,2,IDGT) + EDFR(IDGT)
     aux44GT(3,3,IDGT) = aux44GT(3,3,IDGT) + EDFR(IDGT)
     aux44GT(4,4,IDGT) = aux44GT(4,4,IDGT) + EDFR(IDGT)
  enddo

  !a Normalisation des elements 'non-P11' de la matrice de diffusion: ###########

  do IDGT = 1, NDGT
     PGT(1,2,IDGT) = aux44GT(1,2,IDGT) / aux44GT(1,1,IDGT)
     PGT(1,3,IDGT) = aux44GT(1,3,IDGT) / aux44GT(1,1,IDGT)
     PGT(1,4,IDGT) = aux44GT(1,4,IDGT) / aux44GT(1,1,IDGT)
     PGT(2,1,IDGT) = aux44GT(2,1,IDGT) / aux44GT(1,1,IDGT)
     PGT(2,2,IDGT) = aux44GT(2,2,IDGT) / aux44GT(1,1,IDGT)
     PGT(2,3,IDGT) = aux44GT(2,3,IDGT) / aux44GT(1,1,IDGT)
     PGT(2,4,IDGT) = aux44GT(2,4,IDGT) / aux44GT(1,1,IDGT)
     PGT(3,1,IDGT) = aux44GT(3,1,IDGT) / aux44GT(1,1,IDGT)
     PGT(3,2,IDGT) = aux44GT(3,2,IDGT) / aux44GT(1,1,IDGT)
     PGT(3,3,IDGT) = aux44GT(3,3,IDGT) / aux44GT(1,1,IDGT)
     PGT(3,4,IDGT) = aux44GT(3,4,IDGT) / aux44GT(1,1,IDGT)
     PGT(4,1,IDGT) = aux44GT(4,1,IDGT) / aux44GT(1,1,IDGT)
     PGT(4,2,IDGT) = aux44GT(4,2,IDGT) / aux44GT(1,1,IDGT)
     PGT(4,3,IDGT) = aux44GT(4,3,IDGT) / aux44GT(1,1,IDGT)
     PGT(4,4,IDGT) = aux44GT(4,4,IDGT) / aux44GT(1,1,IDGT)
  enddo

  !a Facteur d'asymetrie et normalisation de la fonction de phase: ##############

  NUMERGG = 0.0D+00
  DENOMGG = 0.0D+00
  do IDGT = 1, NDGT
     ANGLESOLIDE(IDGT) = 2.0D+00 * PI *                      &
          ( COSDGT(IDGT,1) - COSDGT(IDGT,2) )
     aux44GT(1,1,IDGT) = aux44GT(1,1,IDGT) / ANGLESOLIDE(IDGT)
     NUMERGG = NUMERGG +                                     &
          aux44GT(1,1,IDGT) *                           &
          ( COSDGT(IDGT,1) * COSDGT(IDGT,1) -           &
          COSDGT(IDGT,2) * COSDGT(IDGT,2) ) / 2.0D+00
     DENOMGG = DENOMGG +                           &
          aux44GT(1,1,IDGT) *                 &
          ( COSDGT(IDGT,1) - COSDGT(IDGT,2) )
  enddo
  GG = NUMERGG / DENOMGG

  do IDGT = 1, NDGT
     PGT(1,1,IDGT) = 2.0D+00 * aux44GT(1,1,IDGT) / DENOMGG
  enddo

  !a Impression des parametres de sortie: #######################################

!!$  write(*,*) 
!!$  write(*,*) ' PARAMETRES_DE_SORTIE:'
!!$  write(*,'('' Nombre_total_des_photons_diffuses_= '',i9)') &
!!$       NPHDIF
!!$  write(*,'('' Flux_d_energie_incidente_= '',1e12.6)') &
!!$       FINC
!!$  write(*,'('' Flux_d_energie_diffusee_(RT)_= '',1e12.6)') &
!!$       FDRT
!!$  write(*,'('' Flux_d_energie_diffractee_= '',1e12.6)') &
!!$       FDFR
!!$  write(*,'('' Section_efficace_d_extinction_(1E-12_m2)_= '',1e12.6)') &
!!$       CEXTGO
!!$  write(*,'('' Section_efficace_de_diffusion_(1E-12_m2)_= '',1e12.6)') &
!!$       CSCAGO
!!$  write(*,'('' Albedo_pour_une_diffusion_= '',1e12.6)') &
!!$       PIZERO
!!$  write(*,'('' Facteur_d_asymetrie_= '',1e12.6)') &
!!$       GG
!!$  write(*,'('' Nombre d inclusion atteinte= '',i9)') &
!!$       cont
!!$      write(*,*)
!!$      write(*,*) ' Matrice_de_diffusion: '
!!$      write(*,'(a14,a38,4a39,a7)') '  Thetad      ',          &
!!$                     'P11         P12/P11      P13/P11      ', &
!!$                     'P14/P11      P21/P11      P22/P11      ',&
!!$                     'P23/P11      P24/P11      P31/P11      ',&
!!$                     'P32/P11      P33/P11      P34/P11      ',&
!!$                     'P41/P11      P42/P11      P43/P11      ',&
!!$                     'P44/P11'
!!$      do IDGT = 1, NDGT
!!$         write(*,'(1x,f8.4,16(1x,1e12.6))')                      &
!!$                   ( DGT(IDGT,1) + DGT(IDGT,2) )/2.0D+00,        &
!!$                   ( ( PGT(JP,KP,IDGT), KP = 1, 4 ), JP = 1, 4 )
!!$      enddo
!!$  write(*,*) ''
!!$  write(*,*) '                          IHM OUT '
!!$  write(*,*) '############################################################'
!!$  write(*,*) ''

  IHM_RES(1) =  NPHDIF ! Nombre total de PHotons DIFfuses
  IHM_RES(2) =  FINC   ! Flux_d_energie_incidente
  IHM_RES(3) =  FDRT   ! Flux_d_energie_diffusee
  IHM_RES(4) =  FDFR   ! Flux_d_energie_diffractee
  IHM_RES(5) =  CEXTGO ! Section_efficace_d_extinction
  IHM_RES(6) =  CSCAGO ! Section_efficace_de diffusion
  IHM_RES(7) =  PIZERO ! Albedo_pour_une_diffusion
  IHM_RES(8) =  GG     ! Facteur_d_asymetrie

  do IDGT = 1, NDGT
     THETAOUT(IDGT) = ( DGT(IDGT,1) + DGT(IDGT,2) )/2.0D+00
  enddo

  !a Fin de la deuxieme grande boucle: les particules d'interet #################
  ! 2000 continue
  !a Fin de la premiere grande boucle: les nombres d'onde d'interet #############
  ! 1000 continue

  DEALLOCATE(EDFR, maxEDFR, minEDFR)

end subroutine hm_ihm2p0

!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################
!############################################################################

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Initialisation de qqs tableaux et variables
!
SUBROUTINE Init_var(nmumax,numr,NDGT,NDGTFR,NPHDIF,sumSURFAPP,FDRT,FINC,FDFR,  &
     NNr_b,NNr_s,PIZ,PTETA,Pcum,DRES1,DRES2,DIMS1,DIMS2,aux44GT,&
     EDFR)

  implicit none

  !--- Input  Args.
  integer                                :: nmumax, numr, NDGT, NDGTFR

  !--- Output Args.
  integer                                :: NPHDIF

  real(kind=8)                           :: sumSURFAPP, FDRT, FINC, FDFR
  real(kind=8), dimension(numr)          :: NNr_b, NNr_s, NNr
  real(kind=8), dimension(NDGTFR)        :: EDFR
  real(kind=8), dimension(2,numr)        :: PIZ
  real(kind=8), dimension(4,4,NDGT)      :: aux44GT
  real(kind=8), dimension(2,numr,-nmumax:nmumax)                             &
       :: PTETA, Pcum, DRES1, DRES2,       &
       DIMS1, DIMS2

  NNr_b(:) = 0.0D0
  NNr_s(:) = 0.0D0
  NNr(:) = 0.0D0
  PIZ(:,:) = 999.0D0
  PTETA(:,:,:) = 999.0D0
  Pcum(:,:,:)  = 999.0D0
  DRES1(:,:,:) = 999.0D0
  DRES2(:,:,:) = 999.0D0
  DIMS1(:,:,:) = 999.0D0
  DIMS2(:,:,:) = 999.0D0
  FDFR = 0.0D+00
  FINC = 0.0D+00
  FDRT = 0.0D+00
  NPHDIF = 0
  aux44GT(:,:,:) = 0.0D+00
  EDFR(:) = 0.0D+00
  sumSURFAPP = 0.0D+00

  RETURN

END SUBROUTINE Init_var


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Retourne la position des sommets des faces du cristal
!
SUBROUTINE Sommets(RR,LL,AX,SE)

  implicit none

  !--- Input Args.
  real(kind=8)                            :: RR, LL, AX

  !--- Output Args.
  real(kind=8), dimension(3,10)           :: SE    

  SE(1,1) =  0.0D+00
  SE(2,1) =  RR
  SE(3,1) =  LL / 2.0D+00
  SE(1,2) = -AX
  SE(2,2) =  RR / 2.0D+00
  SE(3,2) =  SE(3,1)
  SE(1,3) =  SE(1,2)
  SE(2,3) = -SE(2,2)
  SE(3,3) =  SE(3,1)
  SE(1,4) =  SE(1,1)
  SE(2,4) = -SE(2,1)
  SE(3,4) =  SE(3,1)
  SE(1,5) =  SE(1,1)
  SE(2,5) = -SE(2,1)
  SE(3,5) = -SE(3,1)
  SE(1,6) = -SE(1,2)
  SE(2,6) = -SE(2,2)
  SE(3,6) = -SE(3,1)
  SE(1,7) = -SE(1,2)
  SE(2,7) =  SE(2,2)
  SE(3,7) = -SE(3,1)
  SE(1,8) =  SE(1,1)
  SE(2,8) =  SE(2,1)
  SE(3,8) = -SE(3,1)
  SE(1,9) = -SE(1,2)
  SE(2,9) =  SE(2,2)
  SE(3,9) =  SE(3,1)
  SE(1,10)= -SE(1,2)
  SE(2,10)= -SE(2,2)
  SE(3,10)=  SE(3,1)

  RETURN

END SUBROUTINE Sommets


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!  subroutine qui permet de calculer la distribution en taille cumulees
!  Cette subroutine provient du Meerhoff program
SUBROUTINE granulo_cumul(numr,par1,par2,rlam,NNr,rmin,rmax,tab_r)

  implicit none

  !--- Input Args.
  integer                                :: numr

  real(kind=8)                           :: par1, par2, rlam

  !--- Local Args.
  integer                                :: i

  real(kind=8), parameter                :: eps = 1.0D-10, cutoff = 1.0D-12
  real(kind=8)                           :: ref, sef, rref, r0, r1, N
  real(kind=8), dimension(1)             :: r, nwithr

  !--- Output Args.
  real(kind=8)                           :: rmin, rmax
  real(kind=8), dimension(numr)          :: NNr, tab_r

  !******************************************************************************
  !                                REMARQUE
  !   par1=reff et par2=veff
  !   rn=partie reelle de l'indice
  !   in=partie imaginaire de l'indice
  !******************************************************************************      
  !******************************************************************************
  !*  search for a value of r such that the size distribution
  !*  is less than the cutoff. Start the search at ref+sef which          
  !*  guarantees that such a value will be found on the TAIL of the       
  !*  distribution.                                                       
  !******************************************************************************
  ref = par1
  sef = dsqrt(par2)
  rref= ref
  r(1) = ref + sef
  r0   = ref
200 CALL ihm_sizedis( par1, par2, r, 1, nwithr )
  IF (nwithr(1) .GT. cutoff) THEN
     r0   = r(1)
     r(1) = 2.D0*r(1)
     GOTO 200
  ENDIF
  r1 = r(1)

  !******************************************************************************
  !*  Now the size distribution assumes the cutoff value somewhere        
  !*  between r0 and r1  Use bisection to find the corresponding r        
  !******************************************************************************
300 r(1) = 0.5D0*(r0+r1)
  CALL ihm_sizedis( par1, par2, r, 1, nwithr )
  IF (nwithr(1) .GT. cutoff) THEN
     r0 = r(1)
  ELSE
     r1 = r(1)
  ENDIF
  IF ((r1-r0) .GT. eps) GOTO 300
  rmax = 0.5D0*(r0+r1)

  !******************************************************************************
  !*  Search for a value of r on the low end of the size distribution     
  !*  such that the distribution falls below the cutoff. There is no      
  !*  guarantee that such a value exists, so use an extra test to see if  
  !*  the search comes very near to r = 0                                 
  !******************************************************************************
  r1 = rref
  r0 = 0.D0
400 r(1) = 0.5D0*r1
  CALL ihm_sizedis( par1, par2, r, 1, nwithr )
  IF (nwithr(1) .GT. cutoff) then
     r1 = r(1)
     IF (r1 .GT. eps) GOTO 400
  ELSE
     r0 = r(1)
  ENDIF

  !******************************************************************************
  !*  Possibly the size distribution goes through cutoff between r0       
  !*  and r1 try to find the exact value of r where this happens by       
  !*  bisection.                                                          
  !*  In case there is no solution, the algorithm will terminate soon.    
  !******************************************************************************
500 r(1) = 0.5D0*(r0+r1)
  CALL ihm_sizedis( par1, par2, r, 1, nwithr )
  IF (nwithr(1) .GT. cutoff) THEN
     r1 = r(1)
  ELSE
     r0 = r(1)
  ENDIF
  IF ((r1-r0) .GT. eps) GOTO 500
  IF (r1 .LE. eps) THEN
     rmin = 0.D0
  ELSE
     rmin = 0.5D0*(r0+r1)
  ENDIF

  !******************************************************************************
  !   Calcul de l'integral de nwithr (Nr) (doit etre egale a 1 car elle a ete 
  !   normee).
  !******************************************************************************      
  CALL integre(numr,rmin,rmax,par1,par2,NNr,N,tab_r)
  !--- normalize it
  NNr(:) = NNr(:) / N

  RETURN

END SUBROUTINE granulo_cumul


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
! Subroutine qui vient du MEERHOFF MIE PROGRAM VERSION 3.0, et qui calcul 
! la distribution en taille gamma standart a 2 parametres
SUBROUTINE ihm_sizedis( par1, par2, r, numr, nwithr )

  !******************************************************************************
  !*  Calculate the size distribution n(r) for the numr radius values     
  !*  contained in array r and return the results through the array nwithr
  !*  The size distributions are normalized such that the integral over   
  !*  all r is equal to one.                                              
  !******************************************************************************

  implicit none

  !--- Input Args.
  integer                                 :: numr

  real(kind=8)                            :: par1, par2
  real(kind=8), dimension(numr)           :: r

  !--- Local Args.
  integer                                 :: i

  real(kind=8)                            :: logC, alpha, b, alpha1

  !--- Other Args.
  real(kind=8)                            :: ihm_gammln

  !--- Output Args.
  real(kind=8), dimension(numr)           :: nwithr

  alpha  = 1.D0/par2 - 3.D0
  b      = 1.D0/(par1*par2)
  alpha1 = alpha+1.D0
  logC   = alpha1*dlog(b)-ihm_gammln(alpha1)
  DO i = 1, numr
     nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
  ENDDO !i

  RETURN

END SUBROUTINE ihm_sizedis


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Subroutine qui calcul l'integral cumulee de nwithr
!
SUBROUTINE integre(numr,rmin,rmax,par1,par2,NNr,N,r)

  implicit none

  !--- Input Args.
  integer                                 :: numr

  real(kind=8)                            :: rmin, rmax, par1, par2

  !--- Local Args.
  integer                                 :: i, j

  real(kind=8)                            :: NN, pas
  real(kind=8), dimension(numr)           :: Nr

  !--- Output Args.
  real(kind=8)                            :: N
  real(kind=8), dimension(numr)           :: NNr, r

  !******************************************************************************
  !                             Initialisation
  !******************************************************************************
  N = 0.0D0
  NN = 0.0D0
  r(:) = 0.0D0
  Nr(:) = 0.0D0
  NNr(:) = 0.0D0

  !******************************************************************************
  !                           Compute the table r
  !******************************************************************************
  pas = (rmax - rmin) / DBLE(numr - 1)
  DO i = 1, numr
     r(i) = rmin + DBLE(i-1) * pas
  ENDDO !i

  !******************************************************************************
  !                      Compute the size distribution (Nr)
  !******************************************************************************
  CALL ihm_sizedis( par1, par2, r, numr, Nr )

  !******************************************************************************
  !                     Compute the cumulative integral (NNr)
  !******************************************************************************
  N = Nr(1)
  NNr(1) = Nr(1)
  DO  j = 2, numr
     N = N+Nr(j)
     NN = 0.0D0
     DO i = 1,j
        NN = NN + Nr(i)
     ENDDO !i
     NNr(j) = NN
  ENDDO !j

  RETURN

END SUBROUTINE integre


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!  Return the value of the inclusion radius position in the table of r
SUBROUTINE incl_jradius(numr,NNr,jr)

  implicit none

  !--- Input Args.
  integer                                 :: numr

  real(kind=8), dimension(numr)           :: NNr

  !--- Local Args.
  integer                                 :: n1, n2, j
  real(kind=8)                            :: ran

  !--- Output Args.
  integer                                 :: jr

  n1 = 1
  n2 = numr
  CALL RANDOM_NUMBER(ran)
610 continue
  j = floor((n1+n2)/2.)
  IF ((n2-n1) .EQ. 1) goto 620
  IF (ran .LT. NNr(j)) goto 615
  n1 = j
  goto 610
615 continue
  n2 = j
  goto 610
620 continue

  jr =j

  IF (j .lt. 1) THEN
     WRITE(6,*)'jr=',j
     STOP 'bugg ds le calcul du r de l inclusion'
  ENDIF

  RETURN

END SUBROUTINE incl_jradius


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Return inclusion optical properties for a given wavelength and radius
!
SUBROUTINE opt_properties(libdir,isb,jr,numr,nmu,nmumax,lambda,r,rmu,pmu,drni,dini,   &
     piz,pteta,dres1,dres2,dims1,dims2)

  implicit none

  !--- Input Args.
  integer                                 :: nmumax, nmu, isb, numr, jr

  real(kind=8)                            :: lambda, r
  real(kind=8)                            :: drni, dini
  real(kind=8), dimension(-nmumax:nmumax) :: rmu, pmu

  !--- Local Args.
  character(200)                          :: libdir
  character(50)                           :: wavdir
  character(50)                           :: slam, sr, snr
  character(200)                          :: file_opt_prop

  integer                                 :: iunit

  logical                                 :: exts

  real(kind=8)                            :: w0, Qscat
  real(kind=8), dimension(-nmumax:nmumax) :: p11, res1, res2, ims1, ims2

  !--- Output Args.
  real(kind=8), dimension(2,numr)         :: piz
  real(kind=8), dimension(2,numr,-nmumax:nmumax)                             &
       :: pteta, dres1, dres2, dims1,     &
       dims2

  !******************************************************************************
  !   First of all verify that the optical properties file does not exist
  !******************************************************************************
  !--- wav string and dir
  IF (lambda .LT. 10.0D0) THEN
     WRITE(slam,'(a3,i1,f5.3)')'wav',0,lambda
  ELSE
     WRITE(slam,'(a3,f6.3)')'wav',lambda
  ENDIF
  wavdir = libdir(:len_trim(libdir)) // slam(:len_trim(slam))

  !--- build the opticals properties file name
  IF (r .LT. 10.0D0) THEN
     WRITE(sr,'(a1,i1,f5.3)')'r',0,r
  ELSE
     WRITE(sr,'(a1,f6.3)')'r',r
  ENDIF
  WRITE(snr,'(a1,f5.3)')'n',drni

  file_opt_prop = wavdir(:len_trim(wavdir)) // '/mie_' //                    &
       slam(:len_trim(slam)) // '_' // sr(:len_trim(sr)) // '_'   &
       // snr(:len_trim(snr)) // '.dat'

  !--- Does this file exist?
  INQUIRE(FILE = file_opt_prop(:len_trim(file_opt_prop)), EXIST = exts)

  IF (.NOT. exts) THEN
     !******************************************************************************
     !   If not, compute it and stoke it in the library file
     !******************************************************************************
     !--- First compute optical properties with a mie code (subroutine pmie)
     CALL pmie(nmu,nmumax,r,lambda,drni,dini,rmu,pmu,w0,p11,res1,res2,ims1,   &
          ims2,Qscat)
     piz(isb,jr) = w0
     pteta(isb,jr,:) = p11(:)
     dres1(isb,jr,:) = res1(:)
     dres2(isb,jr,:) = res2(:)
     dims1(isb,jr,:) = ims1(:)
     dims2(isb,jr,:) = ims2(:)

     !--- Open the output file
     iunit = 2
     OPEN(iunit,FILE=file_opt_prop,FORM='formatted',STATUS='new')

     !--- Now stoke results in file file_opt_prop
     CALL write_coefs(iunit,nmumax,nmu,lambda,r,drni,dini,rmu,w0,p11,res1,    &
          res2,ims1,ims2,Qscat)
     CLOSE(iunit)

  ELSE
     !******************************************************************************
     !   the optical properties file should already exist, so Read it!!!
     !******************************************************************************
     iunit = 2
     OPEN(iunit,FILE=file_opt_prop,FORM='formatted',STATUS='old')
     CALL read_coefs(iunit,nmumax,nmu,w0,p11,res1,res2,ims1,ims2)
     piz(isb,jr) = w0
     pteta(isb,jr,:) = p11(:)
     dres1(isb,jr,:) = res1(:)
     dres2(isb,jr,:) = res2(:)
     dims1(isb,jr,:) = ims1(:)
     dims2(isb,jr,:) = ims2(:)
     CLOSE(iunit)

  ENDIF

  RETURN

END SUBROUTINE opt_properties


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Mie program from BB
!
SUBROUTINE pmie(nmu,nmumax,r,lambda,rn,in,rmu,pmu,w0,p11,res1,res2,ims1,ims2,  &
     Qscat)

  implicit none

  !--- Input Atgs.
  integer                                 :: nmu, nmumax

  real(kind=8)                            :: r, lambda, rn, in
  real(kind=8), dimension(-nmumax:nmumax) :: rmu, pmu

  !--- Local Args.
  integer                                 :: n

  real(kind=8), parameter                 :: EPS = 1.0D-06,                  &
       pi = 3.1415926535897932384626433832795029D+00
  real(kind=8)                            :: alpha, Qext, asym
  real(kind=8), dimension(:), allocatable :: ra, rb, ia, ib


  !--- Output Args.
  real(kind=8)                            :: w0, Qscat
  real(kind=8), dimension(-nmumax:nmumax) :: p11, res1, res2, ims1, ims2

  !******************************************************************************
  !   Compute the Mie serie, scattering and extinction coefficient as well as
  !   asymetry parameter
  !******************************************************************************
  alpha = 2.0d+00 * pi * r / lambda
  n = floor(alpha + alpha + 5)
  ALLOCATE(ra(0:n+1),rb(0:n+1),ia(0:n+1),ib(0:n+1))
  CALL mie_serie(n,nmumax,nmu,rmu,pmu,alpha,lambda,rn,in,Qext,Qscat,asym,ra, &
       rb,ia,ib)

  !--- Compute singlescattering albedo
  w0 = Qscat / Qext
  IF (w0 .LT. 0.0D0) THEN
     WRITE(6,*)'Something is strange, w0 < 0 :',w0
     w0 = 0.0D0
  ELSE IF (w0 .GT. (1.0D0+EPS)) THEN
     WRITE(6,*)'Something is strange, w0 > 1 :',w0
     w0 = 1.0D0
  ENDIF

  !******************************************************************************
  !   Compute the coefficients of the phase matrix, and the phase function
  !******************************************************************************
  CALL phase_mat_coef(nmumax,nmu,n,ra,rb,ia,ib,rmu,pmu,alpha,Qscat,p11,res1, &
       res2,ims1,ims2)

  DEALLOCATE(ra,rb,ia,ib)

  RETURN

END SUBROUTINE pmie


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Calcul des series de Mie, facteurs de diffusion et d'extinction,
!   et facteur d'assymetrie pour un parametre de Mie alpha donne :
SUBROUTINE mie_serie(n2,nmumax,nmu,rmu,pmu,alpha,lambda,rn,in,Qext,Qscat,asym, &
     ra,rb,ia,ib)

  implicit none

  !--- Input Args.
  integer                                 :: n2
  integer                                 :: nmumax, nmu

  real(kind=8)                            :: lambda, rn, in, alpha
  real(kind=8), dimension(-nmumax:nmumax) :: rmu, pmu

  !--- Local Args.
  integer                                 :: i, kk, j, i1, n, n22, n11
  integer                                 :: r1, n1, n2p1, a2
  integer                                 :: irl, test, un, l1

  real(kind=8), parameter                 ::                                 &
       pi = 3.1415926535897932384626433832795029D+00
  real(kind=8)                            :: rbeta, ibeta
  real(kind=8)                            :: rsurl, frc
  real(kind=8)                            :: x, z, y, w, q, t
  real(kind=8)                            :: xp, yp, zp, tp
  real(kind=8)                            :: x1, x2, x3, x4, x5, x6, x7, x8, &
       x9
  real(kind=8)                            :: y1, y2, y3, y4, y5, y6, y7, y8
  real(kind=8), dimension(:), allocatable :: cna, rgna, igna, sna
  real(kind=8), dimension(:), allocatable :: rdna, rdnb, idnb   

  !--- Output Args.
  real(kind=8)                            :: Qext, Qscat, asym
  real(kind=8), dimension(0:n2+1)         :: ra, rb, ia, ib   

  r1 = 20
  rbeta = rn * alpha
  ibeta = in * alpha
  n1 = floor(alpha + alpha + r1)
  n2p1 = n2 + 1

  !--- Allocate table
  ALLOCATE(cna(-1:n1),rgna(-1:n1),igna(-1:n1),sna(-1:n1),rdna(0:n1),         &
       rdnb(0:n1),idnb(0:n1))

  !--- Initialize those tables
  cna(:) = 0.0D0
  rgna(:) = 0.0D0
  igna(:) = 0.0D0
  sna(:) = 0.0D0
  rdna(:) = 0.0D0
  rdnb(:) = 0.0D0
  idnb(:) = 0.0D0
  ra(:) = 0.0D0
  rb(:) = 0.0D0
  ia(:) = 0.0D0
  ib(:) = 0.0D0

  rsurl = alpha / 2.0D0 / pi + 1.0D-10
  irl = floor(rsurl)
  frc = rsurl - DBLE(irl)
  IF (frc .LT. 1.0D-08) alpha = alpha * (1.0D0 + 1.0d-08)
  cna(-1) = -DSIN(alpha)
  rgna(-1) = 0.0D0
  rgna(0) = 0.0D0
  cna(0) = DCOS(alpha)
  igna(-1) = 0.0D0
  igna(0) = -1.0D0

  !--- Initialisation de n22 et n11
  n22 = n2
  n11 = n1
  DO i = 1,n2
     cna(i) = DBLE(2*i-1) * cna(i-1) / alpha - cna(i-2)
     x = rgna(i-1)
     z = DBLE(i) / alpha
     y = igna(i-1)
     w = ((z-x) * (z-x) + (y*y))
     rgna(i) = (z-x) / w - z
     igna(i) = y / w
     IF (cna(i) .GE. 1.0d+100) THEN
        WRITE(6,*)'utilisation de n22!'
        n22 = i
        n2p1 = i + 1
        n11 = i + 15
        GOTO 26
     ENDIF
  ENDDO !i
26 CONTINUE

  rdna(n11) = 0.0D0
  rdnb(n11) = 0.0D0
  idnb(n11) = 0.0D0
  x1 = rbeta*rbeta + ibeta*ibeta
  x2 = rbeta/x1
  x3 = ibeta/x1
  sna(n11) = 0.0D0
  sna(n11-1) = 1.0D0

  DO kk = 1,n11
     i = n11 - kk
     x = rdnb(i+1)
     y = idnb(i+1)
     z = x + DBLE(i+1)*x2
     w = y - DBLE(i+1)*x3
     x4 = z*z + w*w
     rdnb(i) = DBLE(i+1)*x2 - z/x4
     idnb(i) = -DBLE(i+1)*x3 + w/x4
     z = DBLE(i+1)/alpha
     x = rdna(i+1)
     rdna(i) = z- 1.0D0/(x+z)
     sna(i-1) = DBLE(2*i + 1)*sna(i)/alpha - sna(i+1)
     IF (sna(i-1) .GT. 1.0d+60) THEN
        test = i-1
        x = sna(test)
        DO j = test,n22
           sna(j) = sna(j)/x
        ENDDO !j
     ENDIF
  ENDDO !kk

  q = -sna(0)/cna(-1)
  DO i1 = 1,n2p1
     sna(i1-1) = sna(i1-1)/q
  ENDDO !i1

  test = 0
  un = 1
  DO i = 1,n22
     x1 = sna(i)
     x2 = cna(i)
     x3 = rdnb(i)
     x4 = idnb(i)
     x5 = rdna(i)
     x6 = rgna(i)
     x7 = igna(i)
     y1 = x3 - rn*x5
     y2 = x4 - in*x5
     y3 = x3 - rn*x6 + in*x7
     y4 = x4 - rn*x7 - in*x6
     y5 = rn*x3 - in*x4 - x5
     y6 = in*x3 + rn*x4
     y7 = rn*x3 - in*x4 - x6
     y8 = in*x3 + rn*x4 - x7
     x4 = y2*y3 - y1*y4
     x3 = y1*y3 + y2*y4
     x5 = x1*x1 + x2*x2
     x6 = y3*y3 + y4*y4
     x7 = y5*y7 + y6*y8
     x8 = y6*y7 - y5*y8
     x9 = y7*y7 + y8*y8
     q = DBLE(i+i+1)/DBLE(i)/DBLE(i+1)*DBLE(un)
     y1 = x1*(x1*x3 + x2*x4)/x5/x6
     y2 = x1*(x1*x4 - x2*x3)/x5/x6
     y3 = x1*(x1*x7 + x2*x8)/x5/x9
     y4 = x1*(x1*x8 - x2*x7)/x5/x9
     ra(i) = y2*q
     ib(i) = y3*q
     q = -q
     rb(i) = y4*q
     ia(i) = y1*q
     un = -un
  ENDDO !i

  !--- Deallocate table
  DEALLOCATE(sna,cna,rdnb,idnb,rdna,rgna,igna)

  ra(0) = 0.0D0
  ia(0) = 0.0D0
  rb(0) = 0.0D0
  ib(0) = 0.0D0
  ra(n2p1) = 0.0D0
  ia(n2p1) = 0.0D0
  rb(n2p1) = 0.0D0
  ib(n2p1) = 0.0D0
  Qext = 0.0D0
  Qscat = 0.0D0
  asym = 0.0D0
  j = -1
  DO n = 1,n22
     l1 = n*n
     a2 = (n+1)*(n+1)
     x = ra(n)
     y = ia(n)
     z = rb(n)
     t = ib(n)
     xp = ra(n+1)
     yp = ia(n+1)
     zp = rb(n+1)
     tp = ib(n+1)
     Qext = Qext + DBLE(n*(n+1)*j)*(y-t)
     Qscat = Qscat + DBLE(l1*a2)/DBLE(n+n+1)*(x*x+y*y+z*z+t*t)
     y1 = DBLE(l1*(n+2)*(n+2)*(n+1))/DBLE(2*n+1)/DBLE(2*n+3)
     asym = asym - y1*(x*xp + y*yp + z*zp + t*tp)
     asym = asym - DBLE(n*(n+1))/DBLE(2*n+1)*(x*z + y*t)
     j = -j
  ENDDO !n

  Qext = Qext*2.0D0/alpha/alpha
  Qscat = Qscat*2.0D0/alpha/alpha
  asym = asym*4.0D0/alpha/alpha/Qscat

  RETURN

END SUBROUTINE mie_serie


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Return coefficients of the phase matrix as well as the phase function 
!
SUBROUTINE phase_mat_coef(nmumax,nmu,n,ra,rb,ia,ib,rmu,pmu,alpha,Qscat,p11,    &
     res1,res2,ims1,ims2)

  implicit none

  !--- Input Args.
  integer                                 :: nmumax, nmu, n

  real(kind=8)                            :: alpha, Qscat
  real(kind=8), dimension(0:n)            :: ra, rb, ia, ib
  real(kind=8), dimension(-nmumax:nmumax) :: rmu, pmu

  !--- Local Args.
  integer                                 :: j, i

  real(kind=8), parameter                 ::                                 &
       Pii = 3.1415926535897932384626433832795029D0
  real(kind=8)                            :: coef, x, tau, pip, pi, pim
  real(kind=8)                            :: rs1, rs2, is1, is2
  real(kind=8)                            :: y1, y2, y3, y4

  !--- Output Args.
  real(kind=8), dimension(-nmumax:nmumax) :: p11, res1, res2, ims1, ims2

  !******************************************************************************
  !   Initialisation de qqs tableaux et coef
  !******************************************************************************
  coef = 2.0D0 / Qscat / (alpha * alpha)
  p11(:) = 0.0D0
  res1(:) = 0.0D0
  res2(:) = 0.0D0
  ims1(:) = 0.0D0
  ims2(:) = 0.0D0

  !******************************************************************************
  !   Compute phase matrix coefficients
  !******************************************************************************
  DO j = -nmu,nmu
     x = -rmu(j)
     pim = 0.0D0
     pi = 1.0D0
     tau = x
     rs1 = 0.0D0
     rs2 = 0.0D0
     is1 = 0.0D0
     is2 = 0.0D0
     DO i = 1,n
        rs1 = rs1 - ia(i)*pi - ib(i)*tau
        rs2 = rs2 + ia(i)*tau + ib(i)*pi
        is1 = is1 + ra(i)*pi + rb(i)*tau
        is2 = is2 - ra(i)*tau - rb(i)*pi
        pip = (DBLE(2*i+1)*x*pi - DBLE(i+1)*pim)/DBLE(i)
        pim = pi
        pi = pip
        tau = DBLE(i+1)*x*pi - DBLE(i+2)*pim
     ENDDO !i
     y1 = rs1**2 + is1**2
     y2 = rs2**2 + is2**2
     y3 = 2.0D0*rs2*rs1
     y4 = 2.0D0*is2*is1
     res1(j) = rs1 / DSQRT(Pii)
     ims1(j) = is1 / DSQRT(Pii)
     res2(j) = rs2 / DSQRT(Pii)
     ims2(j) = is2 / DSQRT(Pii)
     p11(j) = coef*(y1 + y2)
  ENDDO !j

  RETURN

END SUBROUTINE phase_mat_coef


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Stoke the phase matrix coefficient of inclusion
!
SUBROUTINE write_coefs(iunit,nmumax,nmu,lambda,r,nr,ni,rmu,w0,p11,res1,res2,   &
     ims1,ims2,Qscat)

  implicit none

  integer                                 :: j
  integer                                 :: iunit, nmumax, nmu

  real(kind=8), parameter                 :: R2D = 57.295779513D0
  real(kind=8)                            :: w0, lambda, r, nr, ni, Qscat
  real(kind=8), dimension(-nmumax:nmumax) :: p11, res1, res2, ims1, ims2,    &
       rmu

  WRITE(iunit,*)' INCLUSION OPTICAL PROPERTIES FOR THE IHM PROGRAM'
  WRITE(iunit,*)' ------------------------------------------------'
  WRITE(iunit,*)' '
  WRITE(iunit,'(a,f6.3,a)')' Inclusion radius :',r,' mic'
  WRITE(iunit,'(a,f6.3,a)')' Wavelength :',lambda,' nm'
  WRITE(iunit,'(a,f6.4)')' Real part of nr :',nr
  WRITE(iunit,'(a,f6.4)')' Imaginary part of nr :',ni
  WRITE(iunit,'(a,f8.6)')' Qscat :',Qscat
  WRITE(iunit,'(a,f8.6)')' w0 :',w0
  WRITE(iunit,*)' '
  WRITE(iunit,'(a85)')                                                       &
       '  j    thetad         P11           ReS1           ImS1           ReS2           ImS2'
  DO j = nmu,-nmu,-1
     WRITE(iunit,100)j,DACOS(rmu(j))*R2D,p11(j),res1(j),ims1(j),res2(j),      &
          ims2(j)
  ENDDO !j

100 FORMAT(1x,i4,2x,f7.3,3x,e13.7,2x,e13.7,3(2x,e13.7))

  RETURN

END SUBROUTINE write_coefs


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Read the phase matrix coefficient of inclusion
!
SUBROUTINE read_coefs(iunit,nmumax,nmu,w0,p11,res1,res2,ims1,ims2)

  implicit none

  !--- Input Args.
  integer                                 :: iunit, nmumax, nmu

  !--- Local Args.
  character(200)                          :: commentaire
  character(5)                            :: strg

  integer, parameter                      :: nligne2pass = 8
  integer                                 :: i, j, l

  real(kind=8)                            :: thetad, r

  !--- Output Args.
  real(kind=8)                            :: w0
  real(kind=8), dimension(-nmumax:nmumax) :: p11, res1, res2, ims1, ims2

  !--- Read the single scattering albedo
  DO i = 1,nligne2pass
     READ(iunit,'(a200)')commentaire
  ENDDO !i
  READ(iunit,'(a5,f8.6)')strg,w0
  READ(iunit,'(a200)')commentaire
  READ(iunit,'(a200)')commentaire

  !--- Read phase function and coefficients of the phase matrix
  DO j = nmu,-nmu,-1
     READ(iunit,200)l,thetad,p11(j),res1(j),ims1(j),res2(j),      &
          ims2(j)
  ENDDO !j

200 FORMAT(1x,i4,2x,f7.3,3x,e13.7,2x,e13.7,3(2x,e13.7))

  RETURN

END SUBROUTINE read_coefs


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Return the cumulative phase function
!
SUBROUTINE P11_cumul(nmu,nmumax,rmu,pmu,P11,Pcum)

  implicit none

  !--- Input Args.
  integer                                 :: nmumax, nmu

  real(kind=8), dimension(-nmumax:nmumax) :: rmu, pmu, P11

  !--- Local Args.
  integer                                 :: j, i

  real(kind=8)                            :: som

  !--- Output Args.
  real(kind=8), dimension(-nmumax:nmumax) :: Pcum

  !******************************************************************************
  !   Compute the integrated P11
  !******************************************************************************
  som = 0.0D0
  DO j = -nmu,nmu
     som = som + P11(j) * pmu(j)
  ENDDO !j
  !<TMP>
  !    WRITE(6,*)'som=',som
  !    pause
  !</TMP>
  IF (som .LT. 1.9D0) THEN
     WRITE(6,*)'Strange the integrated phase function is less than 1.9 :',som
  ENDIF

  !******************************************************************************
  !   Now compute the cumulative phase function
  !******************************************************************************
  Pcum(:) = 0.0D0
  DO j = -nmu,nmu
     DO i = j,nmu
        Pcum(j) = Pcum(j) + P11(i) * pmu(i)
     ENDDO !i
  ENDDO !j
  Pcum(:) = Pcum(:) / som
  !<TMP>
  !    WRITE(6,*)'Pcum(-125:-115)=',Pcum(-125:-115)
  !    WRITE(6,*)'Pcum(115:125)=',Pcum(115:125)
  !</TMP>

  RETURN

END SUBROUTINE P11_cumul


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Return the angle thetad and phid between the incident photon and the 
!   scatterred one by an inclusion
!   Rq: thetad and phid are in radian!
!
SUBROUTINE Angle_diff(nmu,nmumax,rmu,Pcum,res1,res2,ims1,ims2,ReCS,ImCS,        &
     ReAPARi,ReAPERi,jt,thetad,phid)

  implicit none

  !--- Input Args.
  integer                                 :: nmu, nmumax

  real(kind=8)                            :: ReAPARi, ReAPERi
  real(kind=8), dimension(2,2)            :: ReCS, ImCS
  real(kind=8), dimension(-nmumax:nmumax) :: rmu, Pcum
  real(kind=8), dimension(-nmumax:nmumax) :: res1, res2, ims1, ims2

  !--- Local Args.
  integer                                 :: i
  integer                                 :: n1, n2, jp

  real(kind=8), parameter                 :: eps = 1.0D-08
  real(kind=8)                            :: pinmu
  real(kind=8)                            :: rtd, rpd
  real(kind=8)                            :: pente, mud, B
  real(kind=8)                            :: IS1, IS2, ISS
  real(kind=8)                            :: ReAPAR, ReAPER, ImAPAR, ImAPER
  real(kind=8)                            :: NRJ, Qxy, Uxy
  real(kind=8), dimension(2*nmumax + 1)   :: Fphi

  !--- Output Args.
  integer                                 :: jt
  real(kind=8)                            :: thetad, phid

  !******************************************************************************
  !   First find the indice (jt) that gives Pcum(jt) = rtd (random number)
  !******************************************************************************
  CALL RANDOM_NUMBER(rtd)
  n1 = -nmu
  n2 = nmu
83 CONTINUE
  jt = floor((n1+n2)/2.)
  IF (n2-n1 .EQ. 1) GOTO 85
  IF (rtd .LT. Pcum(jt)) GOTO 84
  n2 = jt
  GOTO 83
84 CONTINUE
  n1 = jt
  GOTO 83
85 CONTINUE
  IF ((jt .LT. 0) .AND. (jt .NE. -nmu)) jt = jt - 1

  !******************************************************************************
  !   Interpolate on the cosine of the scattering angle et calcul de ISS, 
  !   utilise par la suite pour le calcul de phid
  !******************************************************************************
  pente = (Pcum(jt+1) - Pcum(jt)) / (rmu(jt+1) - rmu(jt))
  IF (DABS(pente) .LT. eps) THEN
     mud = rmu(jt)
  ELSE
     B = Pcum(jt) - pente * rmu(jt)
     mud = (rtd - B) / pente
  ENDIF
  IF (mud .GT. 1.0D+00) mud = 1.0D+00
  IF (mud .LT. -1.0D+00) mud = -1.0D+00

  !--- Look for the forward direction
  IF (mud .GE. (1.0D0 - eps)) THEN 
     thetad = 0.0D0
     phid   = 0.0D0
     RETURN
  ENDIF
  thetad = DACOS(mud)
  IS1 = res1(jt)*res1(jt) + ims1(jt)*ims1(jt)
  IS2 = res2(jt)*res2(jt) + ims2(jt)*ims2(jt)
  IF ((IS1+IS2) .LT. eps) THEN
     WRITE(6,*)' Attention IS1+IS2 petit!! c est bizarre!'
     WRITE(6,*)' IS1 =',IS1
     WRITE(6,*)' IS2 =',IS2
     ISS = 0.0D0
     pause
  ELSE
     ISS = (IS1 - IS2) / (IS1 + IS2)
  ENDIF

  !******************************************************************************
  !   Normalisation des champs parallele et perpendiculaire et calcul des 
  !   quantites Qxy = Ex**2 - Ey**2 et Uxy = Ex*conj(Ey) + conj(Ex)*Ey
  !   Rq: 
  !    Ex correspond au champ parallele au plan de reference avant diffusion
  !    Ey correspond au champ perpendiculaire 
  !******************************************************************************
  !--- D'abord calcul des champs perpendiculaires et parallele incident sur 
  !    l'inclusion (a partir de la matrice d'amplitude lie a la transmission et 
  !    des champs incidents)
  ReAPAR = ReAPARi * ReCS(1,1) + ReAPERi * ReCS(1,2)
  ReAPER = ReAPARi * ReCS(2,1) + ReAPERi * ReCS(2,2)
  ImAPAR = ReAPARi * ImCS(1,1) + ReAPERi * ImCS(1,2)
  ImAPER = ReAPARi * ImCS(2,1) + ReAPERi * ImCS(2,2)

  !--- Normalisation
  NRJ = ReAPAR*ReAPAR + ReAPER*ReAPER + ImAPAR*ImAPAR + ImAPER*ImAPER
  IF (NRJ .LT. eps) THEN 
     WRITE(6,*)" Photon's NRJ to law =",NRJ
     WRITE(6,*)" We should not be in this subroutine Angle_diff"
     STOP
  ENDIF
  ReAPAR = ReAPAR / DSQRT(NRJ)
  ReAPER = ReAPER / DSQRT(NRJ)
  ImAPAR = ImAPAR / DSQRT(NRJ)
  ImAPER = ImAPER / DSQRT(NRJ)

  !--- Compute Qxy and Uxy
  Qxy = ReAPAR*ReAPAR + ImAPAR*ImAPAR - ReAPER*ReAPER - ImAPER*ImAPER
  Uxy = 2.0D0*(ReAPAR*ReAPER + ImAPAR*ImAPER)

  !******************************************************************************
  !   Compute phid, but before find Fphi the conditional repartition function
  !   that gives phid knowing thetad...
  !******************************************************************************
  !--- Compute Fphi
  CALL PHIALT(nmumax,nmu,Qxy,Uxy,ISS,pinmu,Fphi)

  !--- Find jp that gives Fphi(jp) = rpd (random number)
  CALL RANDOM_NUMBER(rpd)
  n1 = 1
  n2 = 2*nmu+1
93 CONTINUE
  jp = floor((n1+n2)/2.)
  IF ((n2-n1) .EQ. 1) GOTO 95
  IF (rpd .LT. Fphi(jp)) GOTO 94
  n1 = jp
  GOTO 93
94 CONTINUE
  n2 = jp
  GOTO 93
95 CONTINUE

  !--- Fine interpolation to get the scattered azimuthal angle phid...
  pente = Fphi(jp+1) - Fphi(jp)
  IF (DABS(pente) .LT. EPS) THEN
     phid = DBLE(jp) * pinmu
  ELSE
     B = Fphi(jp) - pente * DBLE(jp)
     phid = pinmu * (rpd - B) / pente
  ENDIF

  RETURN

END SUBROUTINE Angle_diff


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Return the cos. director of the scattered photon in the based linked to the 
!   cristal
!
SUBROUTINE New_direction(thetad,phid,XTR,YTR,ZTR,XNR,YNR,ZNR,ALPHAe,BBETAe,    &
     GAMMAe,ALPHAn,BBETAn,GAMMAn)

  implicit none

  !--- Input Args.
  real(kind=8)                            :: thetad, phid
  real(kind=8)                            :: XTR, YTR, ZTR
  real(kind=8)                            :: XNR, YNR, ZNR
  real(kind=8)                            :: ALPHAe, BBETAe, GAMMAe

  !--- Local Args.
  real(8), parameter                      :: EPS8 = 1.0D-08
  real(8)                                 :: ALPHAd, BBETAd, GAMMAd

  !--- Output Args.
  real(kind=8)                            :: ALPHAn, BBETAn, GAMMAn

  !******************************************************************************
  !   Compute the cos director of the scattered photon (ALPHAd,BBETAd,GAMMAd) in 
  !   the base linked to the reference plan (Oxyz) before the scattering event.
  !   Rq: In this base Oz = direction of incident photon (ALPHAe,BBETAe,GAMMAe)
  !                    Oy = Normal to the reference plan (XNR,YNR,ZNR)
  !                    Ox = Orthonormal to this plan (XTR,YTR,ZTR)
  !******************************************************************************
  ALPHAd = DSIN(thetad)*DCOS(phid)
  BBETAd = DSIN(thetad)*DSIN(phid)
  GAMMAd = DCOS(thetad)

  !<TMP> Verifions que la base lie au plan de reference avant diffusion est bien
  !      Orthonormee
  IF ((XTR*XNR+YTR*YNR+ZTR*ZNR .GT. EPS8) .OR.                               &
       (XNR*ALPHAe+YNR*BBETAe+ZNR*GAMMAe .GT. EPS8) .OR.                      &
       (ALPHAe*XTR+BBETAe*YTR+GAMMAe*ZTR .GT. EPS8)) THEN
     WRITE(6,*)'Pb ds New_direction, repere non orthonormee!!'
     WRITE(6,*)XTR*XNR+YTR*YNR+ZTR*ZNR
     WRITE(6,*)XNR*ALPHAe+YNR*BBETAe+ZNR*GAMMAe
     WRITE(6,*)ALPHAe*XTR+BBETAe*YTR+GAMMAe*ZTR
  ENDIF
  !</TMP>

  !******************************************************************************
  !   Compute the coordinate of this scattered photon (ALPHAn,BBETAn,GAMMAn) in 
  !   the base linked to the cristal (OXYZ).
  !
  !   Column of the matrix (M) that allow to go from the reference plan based 
  !   (Oxyz) to the one linked to the cristal (OXYZ) is just composed of the 
  !   cos director of the orthonormal, normal and incident photon in the base 
  !   (OXYZ), i.e.
  !
  !                         ( XTR  XNR  ALPHAe )
  !         M(Oxyz->OXYZ) = ( YTR  YNR  BBETAe )
  !                         ( ZTR  ZNR  GAMMAe )
  !   Therefore
  !              ( ALPHAn )       ( ALPHAd )
  !              ( BBETAn ) = M * ( BBETAd )
  !              ( GAMMAn )       ( GAMMAd )
  !******************************************************************************
  ALPHAn = XTR*ALPHAd + XNR*BBETAd + ALPHAe*GAMMAd
  BBETAn = YTR*ALPHAd + YNR*BBETAd + BBETAe*GAMMAd
  GAMMAn = ZTR*ALPHAd + ZNR*BBETAd + GAMMAe*GAMMAd

  RETURN

END SUBROUTINE New_direction


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Return the function Fphi, which is the conditional repartirion function
!
SUBROUTINE PHIALT(nmumax,nmu,Qxy,Uxy,ISS,pinmu,Fphi)

  implicit none

  !--- Input Args.
  integer                                 :: nmumax, nmu

  real(kind=8)                            :: Qxy, Uxy, ISS

  !--- Local Args.
  integer                                 :: j

  real(kind=8), parameter                 ::                                  &
       pi = 3.1415926535897932384626433832795029D+00
  real(kind=8)                            :: phi

  !--- Output Args.
  real(kind=8)                            :: pinmu
  real(kind=8), dimension(2*nmumax + 1)   :: Fphi

  pinmu = pi / DBLE(nmu)
  DO j = 1,2*nmu+1
     phi = DBLE(j-1)*pinmu
     Fphi(j) = phi / (2.0D0*pi) - (Qxy*DSIN(2.0D0*phi)/2.0D0 +                &
          Uxy*DSIN(phi)*DSIN(phi))*ISS/(2.0D0*pi)
  ENDDO !j 

  RETURN

END SUBROUTINE PHIALT


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Return the cosine and sine of the angle needed to write the electric field
!   in the scattered plan. 
!   Return also the normal and orthonormal to the new reference plan which is
!   the scattering plan
!
SUBROUTINE Angle_rotE(thetad,ALPi,BETi,GAMi,ALPd,BETd,GAMd,XNR,YNR,ZNR,XTR,    &
     YTR,ZTR,SINAEi,COSAEi,XNRn,YNRn,ZNRn,XTRn,YTRn,ZTRn)

  implicit none

  !--- Input Args.
  real(kind=8)                            :: thetad
  real(kind=8)                            :: ALPi, BETi, GAMi
  real(kind=8)                            :: ALPd, BETd, GAMd
  real(kind=8)                            :: XNR, YNR, ZNR,                  &
       XTR, YTR, ZTR

  !--- Local Args.
  real(kind=8), parameter                 :: EPP = 1.0D-12, EPS = 1.0D-08,   &
       pi = 3.1415926535897932384626433832795029D+00
  real(kind=8)                            :: XNd, YNd, ZNd
  real(kind=8)                            :: NNd, NNr, NNt, NNtn

  !--- Output Args.
  real(kind=8)                            :: SINAEi, COSAEi
  real(kind=8)                            :: XNRn, YNRn, ZNRn,               &
       XTRn, YTRn, ZTRn

  !******************************************************************************
  !   Compute the cosine director (XNd,YNd,ZNd) of the normal (Nd) to the
  !   scattering plan, it is just the vectorial product between the incident 
  !   direction and the scattered one!
  !******************************************************************************
  XNd = GAMd*BETi - BETd*GAMi
  YNd = ALPd*GAMi - GAMd*ALPi
  ZNd = BETd*ALPi - ALPd*BETi

  !--- Normalize it, if the norm is null, then incident and scattered photon
  !    are colinear.
  !    1) In the same direction if thetad = 0
  !    2) In the opposite direction if thetad = pi
  NNd = DSQRT(XNd*XNd + YNd*YNd + ZNd*ZNd)
  IF (NNd .LE. EPP) THEN
     IF (thetad .LT. pi/2.0D0) THEN
        COSAEi = 1.0D0            ! No rotation
        SINAEi = 0.0D0
        XNRn = XNR
        YNRn = YNR
        ZNRn = ZNR
        XTRn = XTR
        YTRn = YTR
        ZTRn = ZTR
     ELSE
        COSAEi = -1.0D0           ! Rotation of pi
        SINAEi = 0.0D0
        XNRn =  XNR
        YNRn =  YNR
        ZNRn =  ZNR
        XTRn = -XTR
        YTRn = -YTR
        ZTRn = -ZTR
     ENDIF
     RETURN
  ENDIF
  XNd = XNd / NNd
  YNd = YNd / NNd
  ZNd = ZNd / NNd

  !******************************************************************************
  !   Rotation angle given by his cosine (COSAEi) and his sine (SINAEi)
  !******************************************************************************
  !--- Verify if (XNR,YNR,ZNR) is normalize...
  NNr = XNR*XNR + YNR*YNR + ZNR*ZNR
  IF ((NNr .GT. 1.0D0) .AND. (NNr .LT. (1.0D0 - EPS))) THEN
     WRITE(6,*)'Bizarre la normale au plan de reference avant inclusion ne ',&
          'semble pas normee!!!!'
     WRITE(6,*)'NNr =',NNr
     STOP
  ENDIF

  !--- Same with (XTR,YTR,ZTR) is normalize...
  NNt = XTR*XTR + YTR*YTR + ZTR*ZTR
  IF ((NNt .GT. 1.0D0) .AND. (NNt .LT. (1.0D0 - EPS))) THEN
     WRITE(6,*)"Bizarre l'orthonormale au plan de reference avant ",         &
          "l'inclusion ne semble pas normee!!!!"
     WRITE(6,*)'NNt =',NNt
     STOP
  ENDIF

  !--- Cosine given by the scalar product of the normal to the scattering plan 
  !    and the normal to the plan of reference before the inclusion
  COSAEi = XNd*XNR + YNd*YNR + ZNd*ZNR
  IF (COSAEi .GT.  1.0D0) COSAEi =  1.0D0
  IF (COSAEi .LT. -1.0D0) COSAEi = -1.0D0

  !--- Sine given by the scalar product of the normal to the scattering plan
  !    and the orthonormal to the plan of reference before the inclusion
  !    Rq: le signe moins provient du fait que ce produit scalaire ns donne 
  !        le cosinus de pi-angle_recherche 
  SINAEi = -(XNd*XTR + YNd*YTR + ZNd*ZTR)
  IF (SINAEi .GT.  1.0D0) SINAEi =  1.0D0
  IF (SINAEi .LT. -1.0D0) SINAEi = -1.0D0

  !******************************************************************************
  !   Compute the normal and orthonormal to the new reference plan (scattering
  !   plan)
  !******************************************************************************
  !--- The normal of the scattering plan is actually the normal of the new 
  !    reference plan
  XNRn = XNd
  YNRn = YNd
  ZNRn = ZNd

  !--- The orthonormal is the vectorial product between the above normal and 
  !    the new direction of the photon ....
  XTRn = GAMd*YNRn - BETd*ZNRn
  YTRn = ALPd*ZNRn - GAMd*XNRn
  ZTRn = BETd*XNRn - ALPd*YNRn
  NNTn = XTRn*XTRn + YTRn*YTRn + ZTRn*ZTRn
  IF (NNTn .LT. 0.5D0) THEN
     WRITE(6,*)' Something strange happen!'
     WRITE(6,*)' The new cosine director or the new orthonormal are not'
     WRITE(6,*)' normalized!!!!!'
     WRITE(6,*)' XNRn=',XNRn,' YNRn=',YNRn,' ZNRn=',ZNRn
     WRITE(6,*)' ALPd=',ALPd,' BETd=',BETd,' GAMd=',GAMd
     STOP
  ENDIF
  IF (NNTn .LT. (1.0D0-EPS)) THEN
     XTRn = XTRn / DSQRT(NNTn)
     YTRn = YTRn / DSQRT(NNTn)
     ZTRn = ZTRn / DSQRT(NNTn)
  ENDIF

  RETURN

END SUBROUTINE Angle_rotE


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
!   Return the norm (or energy) of the electric field
!
SUBROUTINE Field_nrj(ReCS,ImCS,ReAPARi,ReAPERi,NRJ)

  implicit none

  !--- Input Args.
  real(kind=8)                            :: ReAPARi, ReAPERi
  real(kind=8), dimension(2,2)            :: ReCS, ImCS

  !--- Local Args.
  real(kind=8)                            :: ReAPAR, ReAPER, ImAPAR, ImAPER

  !--- Output Args.
  real(kind=8)                            :: NRJ

  !******************************************************************************
  !     Compute the inident electric field (ReAPAR,ReAPER) and (ImAPAR,ImAPER)
  !******************************************************************************
  ReAPAR = ReCS(1,1)*ReAPARi + ReCS(1,2)*ReAPERi
  ReAPER = ReCS(2,1)*ReAPARi + ReCS(2,2)*ReAPERi
  ImAPAR = ImCS(1,1)*ReAPARi + ImCS(1,2)*ReAPERi
  ImAPER = ImCS(2,1)*ReAPARi + ImCS(2,2)*ReAPERi

  !******************************************************************************
  !     Now compute the norm or NRJ
  !******************************************************************************
  NRJ = ReAPAR*ReAPAR + ImAPAR*ImAPAR + ReAPER*ReAPER + ImAPER*ImAPER

  RETURN

END SUBROUTINE Field_nrj


!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

function ihm_gammln(xarg)
  !************************************************************************
  !*  Return the value of the natural logarithm of the gamma function.    *
  !*  The argument xarg must be real and positive.                        *
  !*  This function is documented in :                                    *
  !*                                                                      *
  !*  W.H. Press et al. 1986, 'Numerical Recipes' Cambridge Univ. Pr.     *
  !*  page 157 (ISBN 0-521-30811)                                         *
  !*                                                                      *
  !*  When the argument xarg is between zero and one, the relation (6.1.4)*
  !*  on page 156 of the book by Press is used.                           *
  !*                                         V.L. Dolman April 18 1989    *
  !************************************************************************
  implicit none

  integer                                 :: j

  real(kind=8), parameter                 :: eps = 1.0D-7, one = 1.0D0,      &
       two = 2.0D0, half = 0.5D0,      &
       fpf = 5.5D0,                    &
       stp = 2.50662827465D0,          &
       pi = 3.1415926535897932384626433832795029D+00
  real(kind=8)                            :: xx, x, tmp, ser, gtmp, pix
  real(kind=8)                            :: xarg, ihm_gammln
  real(kind=8), dimension(6)              :: cof

  cof(:) = (/ 76.18009173D0,-86.50532033D0, 24.01409822D0,-1.231739516D0,    &
       0.120858003D-2, -0.536382D-5 /)

  IF (xarg .LE. 0.D0) THEN
     write(6,*) ' ihm_gammln: called with negative argument xarg = ',xarg
     stop 'function ihm_gammln called with negative value'
  ENDIF

  IF (dabs(xarg-one) .LT. eps) THEN
     write(6,*) ' ihm_gammln: argument too close to one for algorithm'
     stop ' in function ihm_gammln argument too close to one'
  ENDIF

  IF (xarg .GE. one) THEN
     xx = xarg
  ELSE
     xx = xarg+two
  ENDIF

  x = xx - one
  tmp = x+fpf
  tmp = (x+half)*dlog(tmp)-tmp
  ser = one
  DO j = 1, 6
     x = x+one
     ser = ser+cof(j)/x
  ENDDO !j
  gtmp = tmp+dlog(stp*ser)
  IF  (xarg .GT. one) THEN
     ihm_gammln = gtmp
  ELSE
     pix = pi*(one-xarg)
     ihm_gammln = dlog(pix/dsin(pix))-gtmp
  ENDIF

  RETURN

end function ihm_gammln

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
SUBROUTINE rot_uvw(vecti,k,angle,vectd)

  implicit none

  !--- Input Arguments
  real(8)                                 :: angle
  real(8), dimension(3)                   :: vecti, k

  !--- Local Arguments
  real(8), parameter                      :: EPS = 1.0D-012
  real(8)                                 :: cs, sn, nu
  real(8)                                 :: ui, vi, wi
  real(8)                                 :: kx, ky, kz
  real(8)                                 :: den

  !--- Output Arguments
  real(8), dimension(3)                   :: vectd


  !**************************************************************************
  !    cs=Cos(angle), sn=sin(angle) and nu = 1 - cos(angle)
  !**************************************************************************
  cs = DCOS(angle)
  sn = DSIN(angle)
  IF (cs .GT.  1.0D0) cs =  1.0D0
  IF (cs .LT. -1.0D0) cs = -1.0D0
  IF (sn .GT.  1.0D0) sn =  1.0D0
  IF (sn .LT. -1.0D0) sn = -1.0D0
  nu = 1.0D0 - cs

  !**************************************************************************
  !   Expression of the composante of the initial vector and k
  !**************************************************************************
  ui = vecti(1)
  vi = vecti(2)
  wi = vecti(3)
  kx = k(1)
  ky = k(2)
  kz = k(3)

  !**************************************************************************
  !   From Ramella-Roman et al., 2005a, compute the new coordinate after
  !   rotation of angle around k
  !**************************************************************************
  vectd(1) = (kx*kx*nu+cs)   *ui + (ky*kx*nu-kz*sn)*vi + (kz*kx*nu+ky*sn)*wi
  vectd(2) = (kx*ky*nu+kz*sn)*ui + (ky*ky*nu+cs)   *vi + (kz*ky*nu-kx*sn)*wi
  vectd(3) = (kx*kz*nu-ky*sn)*ui + (ky*kz*nu+kx*sn)*vi + (kz*kz*nu+cs)   *wi
  den = DSQRT(vectd(1)**2+vectd(2)**2+vectd(3)**2)
  vectd(:) = vectd(:) / den
  !<TMP>
  !    WRITE(6,*)'inside rot cs, sn=',cs,sn
  !    WRITE(6,*)'k=',k(:)
  !    WRITE(6,*)'vectd=',vectd(:)
  !    pause
  !</TMP>

  RETURN

END SUBROUTINE rot_uvw


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! added by MC (rather than read it in a file)

SUBROUTINE ihm_gauleg(x1,x2,x,w,n,ndim)

  implicit none

  !--- Input Arguments
  INTEGER, intent(in) ::  n
  INTEGER, intent(in) :: ndim
  DOUBLE PRECISION, intent(in) :: x1,x2

  !--- Local Arguments
  INTEGER i,j,m

  DOUBLE PRECISION EPS
  PARAMETER (EPS=3.d-14)
  DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1

  !--- Output Arguments
  DOUBLE PRECISION, intent(out) :: x(ndim),w(ndim)

  m = (n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)
  DO i = 1,m
     z = dcos(3.141592654d0*(DBLE(i)-.25d0)/(DBLE(n)+.5d0)) 
1    continue
     p1 = 1.d0
     p2 = 0.d0
     DO j = 1,n
        p3 = p2
        p2 = p1
        p1 = ((2.d0*DBLE(j)-1.d0)*z*p2-(DBLE(j)-1.d0)*p3)/DBLE(j)
     enddo
     pp = n*(z*p1-p2)/(z*z-1.d0)
     z1 = z
     z = z1-p1/pp
     IF(dabs(z-z1) .gt. EPS) goto 1
     x(i) = xm - xl*z
     x(n+1-i) = xm + xl*z
     w(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
     w(n+1-i) = w(i)
  enddo

  RETURN

END SUBROUTINE ihm_gauleg
