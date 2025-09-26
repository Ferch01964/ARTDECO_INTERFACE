! ########################################################################################################### 20061009
!
! CODE 'PHM' (Pristine Hexagonal Monocrystal):
! CALCUL DES PROPRIETES OPTIQUES D'UNE PARTICULE DE GLACE

!  ------ version 2.4: 14/10/2006
!
! * Calcul de theta associe au photon incident, pour une orientation aleatoire
!   a partir de la distribution de probabilite f(theta) = sin(theta) / 2 
!   (vient d'une repartition isotrope en angle solide)
!   D'ou la fct de repartition associee F(mu) = (1 - mu) / 2   mu = [-1,1]
!
! * Initialisation de la fct aleatoire sur l'heure machine grace a la sub
!   RANDOM_SEED
!
! * Declaration des variables en fortran 90 egalement....
!
! ------ version 2.3:
!
! Particule: 
! colonne / plaquette hexagonale de glace pure
!
! Orientation:
! aleatoire dans l'espace
!      
! Methode: 
! 'Geometric Optics' (GO) = 'Ray-Tracing' (RT) avec polarisation, 
!                           + diffraction de Fraunhofer (FR)
!
! ------ sous-version 'simple': lectrure d'un fichier d'entree contenant les
!                               valeurs du nombre d'onde (NU), des parties
!                               reelle (RN) et imaginaire (RM) de l'indice de 
!                               refraction de la glace pure, et de la longueur
!                               (LL) et du rayon (RR) de la particule
!
! ############################################################################################################

subroutine hm_phm2p4(NPHINC, lam, rn, rm, ll, rr, NDGT, & 
     PHM_RES, THETAOUT, PGT)

  !        ************************************************************
  !        *  BROGNIEZ Gerard, C.-LABONNOTE Laurent, DEROO Christine  *
  !        *               et PLANA-FATTORI Artemio                   *
  !        *     Laboratoire d'Optique Atmospherique  UMR 8518        *
  !        *    Universite des Sciences et Technologies de Lille      *
  !        *         59655 Villeneuve d'Ascq Cedex - France           *
  !        *         E-mail: gerard.brogniez@univ-lille1.fr           *
  !        ************************************************************
  !
  !                          Octobre 2006
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

  ! Parametres: ################################################################


  REAL(KIND=8), PARAMETER                    :: ATOTi   = 1.00D+00                                 ! amplitude (totale) incidente;
  REAL(KIND=8), PARAMETER                    :: EPS12   = 1.00D-12                                 ! 'precision' de plusieurs calculs;
  REAL(KIND=8), PARAMETER                    :: EPS6    = 1.00D-06                                 ! 'precision' de plusieurs calculs;


  ! Variables d'entree et de sortie: ###########################################

  !  MC (unsused variable)
  !   REAL(KIND=8), DIMENSION(30)                :: readLL                                             ! longueur de la colonne hexagonale [1E-6 m];
  !   REAL(KIND=8), DIMENSION(1643)              :: readNU                                             ! nombre d'onde (=10000/ALAMB) [cm-1];
  !   REAL(KIND=8), DIMENSION(1643)              :: readRM                                             ! parties imaginaire et reelle de
  !   REAL(KIND=8), DIMENSION(1643)              :: readRN                                             ! l'indice de refraction de la glace;
  !   REAL(KIND=8), DIMENSION(30)                :: readRR                                             ! rayon de la colonne hexagonale [1E-6 m];

  INTEGER, INTENT(IN)                        :: NPHINC                                             ! Nombre de PHotons INCidents;
  REAL(KIND=8), INTENT(IN)                   :: LAM                                                ! longueur d'onde (microns);
  REAL(KIND=8), INTENT(IN)                   :: RM                                                 ! parties imaginaire et reelle de
  REAL(KIND=8), INTENT(IN)                   :: RN                                                 ! l'indice de refraction de la glace;
  REAL(KIND=8), INTENT(IN)                   :: LL                                                 ! longueur de la colonne hexagonale [1E-6 m];
  REAL(KIND=8), INTENT(IN)                   :: RR                                                 ! rayon de la colonne hexagonale [1E-6 m];
  !  modif MC to allow for NDGT control as an input argument 
  INTEGER, INTENT(IN)                        :: NDGT                                               ! Nombre total des Domaines Grand-Theta;


  !  modif MC to allow for NDGT control as an input argument 
  REAL(KIND=8)                               :: RESOLGT                                            ! resolution angulaire de la matrice de
  ! diffusion, en GRANDTHETA [degres];
  INTEGER                                    :: NDGTFR                                             ! NDGT 'difFRaction' (resolution RESOLGT);

  REAL(KIND=8)                               :: NU                                                 ! nombre d'onde (=10000/ALAMB) [cm-1];
  REAL(KIND=8)                               :: Re                                                 ! Rayon de la sphere equivalente en volume;
  REAL(KIND=8)                               :: QQ                                                 ! Facteur de forme;
  REAL(KIND=8)                               :: CSCAGO                                             ! sec. effic. totale de diffusion [1E-12 m2];
  REAL(KIND=8)                               :: CEXTGO                                             ! sec. effic. totale d'extinction [1E-12 m2];
  REAL(KIND=8)                               :: FDFR                                               ! flux d'energie diffractee par la particule;
  REAL(KIND=8)                               :: FINC                                               ! flux d'energie incidente sur la particule;
  REAL(KIND=8)                               :: FDRT                                               ! flux d'energie diffusee (RT) par la part.;
  REAL(KIND=8)                               :: GG                                                 ! facteur d'asymetrie, total [adim.];

  REAL(kind=8), INTENT(OUT)                  :: THETAOUT(NDGT)
  REAL(KIND=8), INTENT(OUT)                  :: PGT(4,4,NDGT)                                      ! matrice de diffusion, comme fonction de 
  ! l'angle GRANDTHETA de diffusion; 
  REAL(KIND=8)                               :: PIZERO                                             ! albedo pour une diffusion [adim.];

  INTEGER                                    :: NPHDIF                                             ! Nombre total de PHotons DIFfuses;

  CHARACTER(LEN=28)                          :: SORTIE 

  REAL(KIND=8), INTENT(OUT), DIMENSION(8)    :: PHM_RES                                            ! TABLE CONTAINING RESULTING SCALAR VALUES
  !  PHM_RES(1) =  NPHDIF : Nombre total de PHotons DIFfuses
  !  PHM_RES(2) =  FINC : Flux_d_energie_incidente
  !  PHM_RES(3) =  FDRT : Flux_d_energie_diffusee
  !  PHM_RES(4) =  FDFR : Flux_d_energie_diffractee
  !  PHM_RES(5) =  CEXTGO : Section_efficace_d_extinction
  !  PHM_RES(6) =  CSCAGO : Section_efficace_de diffusion
  !  PHM_RES(7) =  PIZERO : Albedo_pour_une_diffusion
  !  PHM_RES(8) =  GG     : Facteur_d_asymetrie

  ! Variables intermediaires, fonctions et constantes: #########################

  REAL(KIND=8)                               :: IPARi,IPERi                                        ! intensite incidente: composantes PAR et PER;
  REAL(KIND=8)                               :: IPARd,IPERd                                        ! intensite diffusee: composantes PAR et PER;
  REAL(KIND=8)                               :: ReAPARi,ReAPERi                                    ! amplitude incidente: parties reelles;
  REAL(KIND=8)                               :: ReAPARd,ReAPERd                                    ! amplitude diffusee: parties reelles;
  REAL(KIND=8)                               :: ImAPARd,ImAPERd                                    ! amplitude diffusee: parties imaginaires;
  REAL(KIND=8)                               :: ReCRPAR,ImCRPAR                                    ! coeff. Fresnel, reflexion, PARallele;
  REAL(KIND=8)                               :: ReCRPER,ImCRPER                                    ! coeff. Fresnel, reflexion, PERpendiculaire;
  REAL(KIND=8)                               :: ReCTPAR,ImCTPAR                                    ! coeff. Fresnel, transmission, PAR;
  REAL(KIND=8)                               :: ReCTPER,ImCTPER                                    ! coeff. Fresnel, transmission, PER;
  REAL(KIND=8), DIMENSION(2,2)               :: ReCSRD                                             ! matrice d'amplitudes diffusees
  REAL(KIND=8), DIMENSION(2,2)               :: ImCSRD                                             ! (entree de la routine STOKES);
  REAL(KIND=8), DIMENSION(2,2)               :: ReCR                                               ! matrice d'amplitudes reflechies
  REAL(KIND=8), DIMENSION(2,2)               :: ImCR                                               ! vers l'interieur de la particule;
  REAL(KIND=8)                               :: auxCS                                              ! matrice d'amplitudes diffusees
  REAL(KIND=8), DIMENSION(2,2)               :: ReCS                                               ! (transmises vers l'interieur puis
  REAL(KIND=8), DIMENSION(2,2)               :: ImCS                                               ! refractees vers l'exterieur);
  REAL(KIND=8), DIMENSION(2,2)               :: ReCSA                                              ! matrice d'amplitudes diffusees
  REAL(KIND=8), DIMENSION(2,2)               :: ImCSA                                              ! (rotation A); 
  REAL(KIND=8), DIMENSION(2,2)               :: ReCSD                                              ! matrice d'amplitudes diffusees
  REAL(KIND=8), DIMENSION(2,2)               :: ImCSD                                              ! (rotation D); 
  REAL(KIND=8)                               :: auxCT                                              ! matrice d'amplitudes transmises
  REAL(KIND=8), DIMENSION(2,2)               :: ReCT                                               ! a l'interieur de la particule
  REAL(KIND=8), DIMENSION(2,2)               :: ImCT                                               ! (considerant l'attenuation);

  REAL(KIND=8)                               :: ALPHAe                                             ! 1er cosinus directeur, faisceau emergent;
  REAL(KIND=8), DIMENSION(0:7)               :: ALPHAf                                             ! 1er cosinus directeur, faces de la colonne;
  REAL(KIND=8)                               :: ALPHAi                                             ! 1er cosinus directeur, faisceau incident;
  REAL(KIND=8)                               :: ALPHAr                                             ! 1er cosinus directeur, faisceau reflechi;
  REAL(KIND=8)                               :: ALPHAt                                             ! 1er cosinus directeur, faisceau transmis;
  REAL(KIND=8)                               :: ANGLESOLIDE(NDGT)                                        ! angles solides elementaires
  ! (GRANDTHETA), integres sur GRANDPHI;
  REAL(KIND=8), DIMENSION(4,4)               :: aux44GTGP                                          ! variable auxiliaire (sortie de la routine STOKES);
  REAL(KIND=8)                               :: aux44GT(4,4,NDGT)                                  ! variable auxiliaire (matrice de diffusion
  ! brute, avant normalisation);
  REAL(KIND=8)                               :: auxauxDELTA3                                       ! variable auxiliaire au calcul de DELTA3;
  REAL(KIND=8), DIMENSION(0:7)               :: auxDELTA3                                          ! distances internes a partir d'un impact;
  REAL(KIND=8)                               :: auxSURFAPP                                         ! valeur effective de la surface apparente;
  REAL(KIND=8)                               :: auxXimp                                            ! ordonnees d'un point d'impact autant
  REAL(KIND=8)                               :: auxYimp                                            ! que ORIGINE d'un parcours du photon
  REAL(KIND=8)                               :: auxZimp                                            ! a l'interieur de la particule;
  REAL(KIND=8)                               :: AX                                                 ! = RR * sqrt( 3 ) / 2;
  REAL(KIND=8)                               :: COSANGLEHETERO                                     ! cosinus de l'angle d'heterogeneite;
  REAL(KIND=8)                               :: BBETAe                                             ! 2eme cosinus directeur, faisceau emergent;
  REAL(KIND=8), DIMENSION(0:7)               :: BBETAf                                             ! 2eme cosinus directeur, faces;
  REAL(KIND=8)                               :: BBETAi                                             ! 2eme cosinus directeur, faisceau incident;
  REAL(KIND=8)                               :: BBETAr                                             ! 2eme cosinus directeur, faisceau reflechi;
  REAL(KIND=8)                               :: BBETAt                                             ! 2eme cosinus directeur, faisceau transmis;
  REAL(KIND=8)                               :: COSDGT(NDGT,2)                                     ! cosinus de DGT(IDGT,1 et 2);
  REAL(KIND=8)                               :: COSETA                                             ! cosinus de l'angle ETA;
  REAL(KIND=8)                               :: COSPHI                                             ! cosinus de l'angle PHI;
  REAL(KIND=8)                               :: COSTHETA                                           ! cosinus de l'angle THETA;
  REAL(KIND=8)                               :: COSTHETAi                                          ! cosinus de l'angle THETAi d'incidence;
  REAL(KIND=8)                               :: COSTHETAt                                          ! cosinus de l'angle THETAt de refraction; 
  REAL(KIND=8)                               :: COSGRANDTHETA                                      ! cosinus de l'angle GRANDTHETA;
  REAL(KIND=8)                               :: GRANDTHETA                                         ! angle GRANDTHETA;
  REAL(KIND=8)                               :: DATT                                               ! partie reelle de l'attenuation;
  REAL(KIND=8)                               :: DCOSA                                              ! var. aux. (rotation A);
  REAL(KIND=8)                               :: DCOSB                                              ! var. aux. (rotation B);
  REAL(KIND=8)                               :: DCOSD                                              ! var. aux. (rotation D);
  REAL(KIND=8)                               :: DDDD                                               ! distance interne maximale [1E-12 m2];
  REAL(KIND=8)                               :: DELTA3                                             ! distance entre deux impacts internes;
  REAL(KIND=8)                               :: DENOMGG                                            ! integrations au DENOMinateur et au
  REAL(KIND=8)                               :: NUMERGG                                            ! NUMERateur dans le facteur d'asymetrie;
  REAL(KIND=8)                               :: DGT(NDGT,2)                                        ! Domaines de l'angle GRANDTHETA;
  REAL(KIND=8)                               :: DSINA                                              ! var. aux. (rotation A);
  REAL(KIND=8)                               :: DSINB                                              ! var. aux. (rotation B);
  REAL(KIND=8)                               :: DSIND                                              ! var. aux. (rotation D);
  REAL(KIND=8), allocatable                  :: EDFR(:)                                              ! distr. angul. de l'energie diffractee;
  REAL(KIND=8), allocatable                  :: maxEDFR(:)                                            ! var. aux. (EDFR);
  REAL(KIND=8), allocatable                  :: minEDFR(:)                                            ! var. aux. (EDFR);
  REAL(KIND=8)                               :: ETA                                                ! angle entre le faisceau incident et l'axe OZ' [degres];
  REAL(KIND=8)                               :: maxETA                                             ! var. aux. (EDFR);
  REAL(KIND=8)                               :: minETA                                             ! var. aux. (EDFR);

  REAL(KIND=8)                               :: GAMMAe                                             ! 3eme cosinus directeur, faisceau emergent;
  REAL(KIND=8), DIMENSION(0:7)               :: GAMMAf                                             ! 3eme cosinus directeur, faces;
  REAL(KIND=8)                               :: GAMMAi                                             ! 3eme cosinus directeur, faisceau incident;
  REAL(KIND=8)                               :: GAMMAr                                             ! 3eme cosinus directeur, faisceau reflechi;
  REAL(KIND=8)                               :: GAMMAt                                             ! 3eme cosinus directeur, faisceau transmis;
  REAL(KIND=8)                               :: GRANDK                                             ! variable intermediaire (indice equiv.);
  REAL(KIND=8)                               :: GRANDN                                             ! variable intermediaire (indice equiv.);
  REAL(KIND=8)                               :: KKKK                                               ! = 2*PI/(LongueurD'Onde);
  REAL(KIND=8)                               :: LL2                                                ! la moitie de la longueur de la colonne;
  REAL(KIND=8)                               :: PHI                                                ! angle entre la projection OXp du faisceau incident et l'axe OX [adim.];
  REAL(KIND=8)                               :: maxPHI                                             ! var. aux. (EDFR);
  REAL(KIND=8)                               :: minPHI                                             ! var. aux. (EDFR);
  REAL(KIND=8)                               :: rnETA                                              ! "Random Number" de l'angle ETA;
  REAL(KIND=8)                               :: rnPHI                                              ! "Random Number" de l'angle PHI;
  REAL(KIND=8)                               :: rnPOLAR                                            ! polarisation aleatoire du champ incident;
  REAL(KIND=8)                               :: rnZETA                                             ! "Random Number" de l'angle ZETA;
  REAL(KIND=8), DIMENSION(3,10)              :: SE                                                 ! coordonnees XYZ des dix sommets eclaires;
  REAL(KIND=8)                               :: SINPHI                                             ! sinus de l'angle PHI;
  REAL(KIND=8)                               :: SINTHETA                                           ! sinus de l'angle THETA;
  REAL(KIND=8)                               :: sumP11                                             ! var. aux. (normalis. de P11 ray-tracing);
  REAL(KIND=8)                               :: XND,YND,ZND                                        ! variables auxiliaires
  REAL(KIND=8)                               :: XNI,YNI,ZNI                                        ! (cosinus directeurs divers);
  REAL(KIND=8)                               :: XNP,YNP,ZNP                                        !     "
  REAL(KIND=8)                               :: XNR,YNR,ZNR                                        !     "
  REAL(KIND=8)                               :: XNS,YNS,ZNS                                        !     "
  REAL(KIND=8)                               :: XTD,YTD,ZTD                                        !     "
  REAL(KIND=8)                               :: XTI,YTI,ZTI                                        !     "
  REAL(KIND=8)                               :: XTR,YTR,ZTR                                        !     "
  REAL(KIND=8)                               :: XTS,YTS,ZTS                                        !     "
  REAL(KIND=8)                               :: Ximp,Yimp,Zimp                                     ! (position) impact d'un photon sur une face;
  REAL(KIND=8)                               :: ZETA                                               ! angle entre les axes OZ et OZ' [adim.];
  REAL(KIND=8)                               :: maxZETA                                            ! var. aux. (EDFR);
  REAL(KIND=8)                               :: minZETA                                            ! var. aux. (EDFR);

  INTEGER                                    :: auxFACEimp                                         ! variable auxiliaire (identif. de FACEimp);
  INTEGER                                    :: auxNPNPT                                           ! variable auxiliaire (routine DIFFRACT);
  INTEGER                                    :: ICH                                                ! indicateur des colonnes hexagonales;
  INTEGER                                    :: IDGT                                               ! indicateur du Domaine GrandTheta (NDGT);
  INTEGER                                    :: IFACE                                              ! indicateur des face de la particule;
  INTEGER                                    :: INU                                                ! indicateur des nombres d'onde d'interet;
  INTEGER                                    :: IPHINC                                             ! Indicateur des PHotons INCidents;
  INTEGER                                    :: IREFL                                              ! Indicateur de REFLexion totale;
  INTEGER                                    :: JP                                                 ! indicateurs dans les matrices
  INTEGER                                    :: KP                                                 ! AUX22GTGP, aux44GTGP, aux44GT et PGT;
  INTEGER                                    :: FACEimp                                            ! indicateur de la face atteinte (impactee);

  REAL(KIND=8)                               :: surfapp                                            ! fonction interne 'surface apparente de la colonne hexagonale';
  REAL(KIND=8)                               :: sumSURFAPP                                         ! somme des valeurs auxSURFAPP;
  REAL(KIND=8)                               :: maxSURFAPP                                         ! SURFace APParente maximale; [1E-12 m2];
  REAL(KIND=8)                               :: intSURFAPP                                         ! SURFace APParente intermediaire [1E-12 m2];
  REAL(KIND=8)                               :: minSURFAPP                                         ! SURFace APParente minimale; [1E-12 m2];

  !<CHG,version 2.4>
  REAL(KIND=8)                               :: rnmu
  REAL(KIND=8)                               :: mu
  REAL(KIND=8)                               :: THETA
  REAL(KIND=8)                               :: SINETA
  REAL(KIND=8)                               :: COSZETA
  !</CHG,version 2.4>

  REAL(KIND=8), PARAMETER                    :: PI     = 3.1415926535897932384626433832795029D+00  ! = constante PI;
  REAL(KIND=8), PARAMETER                    :: PI2    = PI / 2.0D+00                              ! = PI / 2;
  REAL(KIND=8), PARAMETER                    :: PI3    = PI / 3.0D+00                              ! = PI / 3;
  REAL(KIND=8), PARAMETER                    :: PI6    = PI / 6.0D+00                              ! = PI / 6;
  REAL(KIND=8), PARAMETER                    :: PI180  = PI / 180.0D+00                            ! = PI / 180;
  REAL(KIND=8), PARAMETER                    :: R32    = 0.8660254037844385965883020617184229D+00  ! = sqrt( 3 )/2;

  !---------------------
  !  modif MC to allow for NDGT control as an input argument 
  RESOLGT = 180.D0 / (NDGT - 1)
  NDGTFR  = idnint( 1.0D+00 +  90.0D+00 / RESOLGT )
  ALLOCATE(EDFR(NDGTFR), maxEDFR(NDGTFR), minEDFR(NDGTFR))

  !--- Compute wavenumber
  NU = 10000.0D0 / LAM

  !--- Calcul de la longueur LL et du rayon RR du cristal
  !      RR = Re*(8.0D0*PI/(9.0D0*DSQRT(3.0D0)*2.0D0*QQ))**(1.0D0/3.0D0)
  !      LL = 2.0D0*RR*QQ

  KKKK = 2.0D+00 * PI / ( 1.0D+04 / NU )
  DDDD = dsqrt( LL * LL + 4.0D+00 * RR * RR )

  LL2 = LL / 2.0D+00
  AX = RR * R32

  ! Ouverture du fichier de sortie: ############################################

!!$      write(*,*) ''
!!$      write(*,*) '############################################################'
!!$      write(*,*) '                          PHM '
!!$      write(*,*) ''
!!$      write(*,*) ' CARACTERISTIQUES_GEOMETRIQUES_DU_PROBLEME:'
!!$      write(*,'('' Longueur_de_la_colonne_(1E-6_m)_= '',1e12.6)') &
!!$                LL
!!$      write(*,'('' Rayon_de_la_colonne_(1E-6_m)_= '',1e12.6)') &
!!$                RR
!!$      write(*,'('' Rayon_de_la_sphere_equiv._(1E-6_m)_= '',1e12.6)') &
!!$                ( 9.0D+00 * RR * RR * LL * dsqrt( 3.0D+00 ) /        &
!!$                  ( 8.0D+00 * PI ) )**( 1.0D+00 / 3.0D+00 )
!!$      write(*,'('' Facteur_de_forme_=_Longueur_/_Largeur_= '',1e12.6)') &
!!$                LL / ( 2.0D+00 * RR )
!!$      write(*,*)
!!$      write(*,*) ' CARACTERISTIQUES_OPTIQUES_DU_PROBLEME:'
!!$      write(*,'('' Nombre_d_onde_(cm-1)_= '',1e12.6)') &
!!$                NU
!!$      write(*,'('' Partie_reelle_de_l_indice_de_refr._= '',1e12.6)') &
!!$                RN
!!$      write(*,'('' Partie_imaginaire_de_l_indice_de_refr._= '',1e12.6)') &
!!$                RM
!!$
!!$      write(*,*)
!!$      write(*,*) ' CARACTERISTIQUES_DU_FAISCEAU_INCIDENT:'
!!$      write(*,'('' Nombre_de_photons_incidents_= '',i9)') &
!!$                NPHINC
!!$      write(*,'('' Amplitude_du_champ_electrique_incident_= '',1e12.6)') &
!!$                1.0D+00
!!$      write(*,*) 'Polarisation_du_champ_electrique_incident_= ', &
!!$                 'aleatoire'
!!$      write(*,*)
!!$      write(*,*) ' ORIENTATIONS_DES_FAISCEAUX_EMERGENTS:'
!!$      write(*,'('' Resol._de_la_matrice_de_diffusion_(degres)_= '',1e12.6)') &
!!$                RESOLGT

  ! Initialisation de plusieurs variables: #####################################

  FDFR = 0.0D+00
  FINC = 0.0D+00
  FDRT = 0.0D+00
  NPHDIF = 0

  do JP = 1, 4
     do KP = 1, 4
        do IDGT = 1, NDGT
           aux44GT(JP,KP,IDGT) = 0.0D+00
        enddo
     enddo
  enddo

  do IDGT = 1, NDGTFR
     EDFR(IDGT) = 0.0D+00
  enddo

  sumSURFAPP = 0.0D+00

  ! Angle de diffusion GRANDTHETA: #############################################
  !
  !    DGT(IDGT,1) = angle inferieur du domaine;
  !    DGT(IDGT,2) = angle superieur du domaine;

  DGT(1,1) = 0.0D+00
  DGT(1,2) = RESOLGT / 2.0D+00
  do IDGT = 2, NDGT-1
     DGT(IDGT,1) = dfloat( IDGT - 1 ) * RESOLGT - RESOLGT / 2.0D+00
     DGT(IDGT,2) = dfloat( IDGT     ) * RESOLGT - RESOLGT / 2.0D+00
  enddo
  DGT(NDGT,1) = 180.0D+00 - RESOLGT / 2.0D+00
  DGT(NDGT,2) = 180.0D+00
  do IDGT = 1, NDGT
     COSDGT(IDGT,1) = dcos( DGT(IDGT,1) * PI180 )
     COSDGT(IDGT,2) = dcos( DGT(IDGT,2) * PI180 )
  enddo

  ! Les sommets de la particule qui sont eclaires: #############################
  !
  !    SE(1,I), SE(2,I) et SE(3,I): ordonnees X, Y et Z du I-eme sommet eclaire'

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

  ! Cosinus directeurs des faces de la particule: ##############################

  do IFACE = 0, 5
     ALPHAf(IFACE) = dcos( dfloat( IFACE ) * PI3 )
     BBETAf(IFACE) = dsin( dfloat( IFACE ) * PI3 )
     GAMMAf(IFACE) = 0.0D+00
  enddo
  do IFACE = 6, 7
     ALPHAf(IFACE) = 0.0D+00
     BBETAf(IFACE) = 0.0D+00
     GAMMAf(IFACE) = dcos( dfloat( IFACE - 6 ) * PI )
  enddo

  ! Surfaces apparentes extremes et leurs orientations: ########################

  call EXTREMES(RR,LL,PI180,R32,maxSURFAPP,minSURFAPP,       &
       maxETA,minETA,maxPHI,minPHI,maxZETA,minZETA)

  ! Contributions extremes associees a la diffraction: #########################

  auxNPNPT = 8 
  if( minSURFAPP.gt.(1.0D+06) ) auxNPNPT = 10
  !<TMP,21/01>
  call DIFFRACT(NU,KKKK,LL,RR,SE,maxSURFAPP,maxETA,maxZETA,maxPHI, &
       NDGT,NDGTFR,RESOLGT,COSDGT,                        &
       EPS12,PI,PI180,PI3,auxNPNPT,maxEDFR)

  call DIFFRACT(NU,KKKK,LL,RR,SE,minSURFAPP,minETA,minZETA,minPHI, &
       NDGT,NDGTFR,RESOLGT,COSDGT,                        &
       EPS12,PI,PI180,PI3,auxNPNPT,minEDFR)
  !</TMP,21/01>

  !<CHG,version 2.4>
  !  Initialisation de la fonction random
  CALL RANDOM_SEED
  !</CHG,version 2.4>

  ! La troisieme grande boucle: les photons incidents ##########################

  do IPHINC = 1, NPHINC

     ! Initialisation de l'attenuation a l'interieur de la particule: #############

     DATT = 1.0D+00

     ! Choix aleatoire de l'orientation du faisceau incident: #####################
     !<CHG,version 2.4>
     !  Ce choix etait mal fait, maintenant utilisation de la densite de proba
     !  f(theta) = sin(theta) / 2

     !      call random_number( rnETA )
     !      ETA = ( rnETA * 90.0D+00 ) * PI180

     call random_number( rnPHI )
     PHI = ( rnPHI * 60.0D+00 - 30.0D+00 ) * PI180

     !      call random_number( rnZETA )
     !      ZETA = ( rnZETA * 90.0D+00 ) * PI180

     CALL RANDOM_NUMBER(rnmu)
     mu = rnmu                        !equivalent a dabs(1-2*rnmu), parceque 
     !l on veut theta entre 0 et 90 degres
     !</CHG,version 2.4>

     ! Section geometrique totale: ################################################

     auxSURFAPP = surfapp(RR,LL,ETA,ZETA,PHI,R32)

     sumSURFAPP = sumSURFAPP + auxSURFAPP

     ! Cosinus directeurs du faisceau incident: ###################################

     COSPHI   = dcos( PHI )
     SINPHI   = dsin( PHI )
     !<CHG,version 2.4>
     !      COSETA   = dcos( ETA )
     !      SINTHETA = dcos( ZETA ) * COSETA
     !      COSTHETA = dsqrt( 1.0D+00 - SINTHETA * SINTHETA )
     THETA = PI / 2.0D0 - DACOS(mu) !parceque theta est defini par rapport
     !au OXY et non par rapport a OZ
     COSTHETA = DCOS(THETA)
     IF (COSTHETA .GT.  1.0D0) COSTHETA =  1.0D0
     IF (COSTHETA .LT. -1.0D0) COSTHETA = -1.0D0
     SINTHETA = DSQRT(1.0D0 - COSTHETA * COSTHETA)

     !--- Calcul de ETA et ZETA a partir de THETA et PHI
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
     !</CHG,version 2.4>

     ALPHAi = - COSTHETA * COSPHI
     BBETAi = - COSTHETA * SINPHI
     GAMMAi = - SINTHETA

     ! Impact du photon incident sur une face externe: ############################

     FACEimp = -999

     call IMPACT(LL2,AX,COSPHI,COSTHETA,SINPHI,SINTHETA,EPS12,SE,R32, &
          Ximp,Yimp,Zimp,FACEimp)

     ! Verification des resultats d'IMPACT, concernant la face atteinte: ##########

     call TEST_IMPACT(NPHDIF,FACEimp,ETA,ZETA,PHI,PI180,EPS12)

     ! APARi et APERi: composantes parallele et perpendiculaire
     !                 de l'amplitude incidente; polarisation aleatoire... ########

     call random_number( rnPOLAR )

     ReAPARi = ATOTi * dcos( 2.0D+00*PI*rnPOLAR )
     ReAPERi = ATOTi * dsin( 2.0D+00*PI*rnPOLAR )

     ! IPARi et IPERi = composantes parallele et perpendiculaire
     !                  de l'intensite incidente... ###############################

     IPARi = ReAPARi*ReAPARi
     IPERi = ReAPERi*ReAPERi

     ! FINC = flux d'energie incidente (intensite totale incidente)... ############

     FINC = FINC + IPARi + IPERi

     ! Cosinus directeurs des faisceaux reflechi et transmis: #####################

     ! Coefficients de reflexion et de transmission de Fresnel: ###################

     call CHANG_AG(RN,RM,EPS12,                     &
          ALPHAf,BBETAf,GAMMAf,FACEimp,    &
          ALPHAi,BBETAi,GAMMAi,COSTHETAi,  &
          ALPHAr,BBETAr,GAMMAr,            &
          ALPHAt,BBETAt,GAMMAt,COSTHETAt,  &
          GRANDN,GRANDK,COSANGLEHETERO,    &
          ReCRPAR,ReCRPER,ReCTPAR,ReCTPER, &
          ImCRPAR,ImCRPER,ImCTPAR,ImCTPER)

     ! Matrice d'amplitude diffusee, par reflexion a la face externe: #############

     ReCSRD(1,1) = ReCRPAR
     ImCSRD(1,1) = ImCRPAR

     ReCSRD(1,2) = 0.0D+00
     ImCSRD(1,2) = 0.0D+00

     ReCSRD(2,1) = 0.0D+00
     ImCSRD(2,1) = 0.0D+00

     ReCSRD(2,2) = ReCRPER
     ImCSRD(2,2) = ImCRPER

     ! APARd et APERd: composantes parallele et perpendiculaire
     !                 de l'amplitude diffusee... #################################

     ReAPARd = ReCSRD(1,1)*ReAPARi
     ReAPERd = ReCSRD(2,2)*ReAPERi
     ImAPARd = ImCSRD(1,1)*ReAPARi
     ImAPERd = ImCSRD(2,2)*ReAPERi

     ! IPARd et IPERd = composantes parallele et perpendiculaire
     !                  de l'intensite diffusee... ################################

     IPARd = ReAPARd*ReAPARd + ImAPARd*ImAPARd
     IPERd = ReAPERd*ReAPERd + ImAPERd*ImAPERd

     ! FDRT = flux d'energie diffusee par Ray-Tracing... ##########################

     FDRT = FDRT + IPARd + IPERd

     ! Angles de diffusion entre les faisceaux incident et diffuse': ##############

     COSGRANDTHETA = ALPHAi * ALPHAr + &
          BBETAi * BBETAr + &
          GAMMAi * GAMMAr

     !<CHG,4/10/2006>
     IF (COSGRANDTHETA .GT. 1.0D0)  COSGRANDTHETA = 1.0D0
     IF (COSGRANDTHETA .LT. -1.0D0) COSGRANDTHETA = -1.0D0
     !</CHG,4/10/2006>
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

     ! Matrice de diffusion brute (= non-normalisee) aux44GTGP, calculee a partir #
     ! des amplitudes diffusees par ray-tracing, AUX22GTGP: #######################

     call STOKES(ReCSRD,ImCSRD,aux44GTGP)

     do JP = 1, 4
        do KP = 1, 4
           aux44GT(JP,KP,IDGT) = aux44GT(JP,KP,IDGT) + aux44GTGP(JP,KP)
        enddo
     enddo

     ! Le faisceau transmis (refracte') vers l'interieur de la particule devient ##
     ! celui qui emerge, pour un parcours interne, a partir de la face interne: ###

     ALPHAe  = ALPHAt
     BBETAe  = BBETAt
     GAMMAe  = GAMMAt

     ! Cosinus directeurs de la normale et de
     ! l'orthonormale au plan d'incidence du dioptre atteint: #####################

     call NORMALEPI(ALPHAi,BBETAi,GAMMAi,                  &
          EPS12,FACEimp,PI,ALPHAf,BBETAf,GAMMAf, & 
          XNI,XTI,YNI,YTI,ZNI,ZTI)

     ! Cosinus directeurs de la normale et
     ! de l'orthonormale au plan de refraction: ###################################

     call NORMALE(ALPHAe,BBETAe,GAMMAe,XNI,YNI,ZNI, &
          XNR,XTR,YNR,YTR,ZNR,ZTR)

     ! Verification que (XNI,YNI,ZNI) et (XNR,YNR,ZNR) sont identiques.....
     !      WRITE(6,*)'NI =',XNI,YNI,ZNI
     !      WRITE(6,*)'NR =',XNR,YNR,ZNR
     !      STOP

     ! Matrice d'amplitudes transmises vers l'interieur: ##########################

     auxCS = dsqrt(GRANDN*COSTHETAt/COSTHETAi)

     ReCS(1,1) = auxCS*ReCTPAR
     ImCS(1,1) = auxCS*ImCTPAR

     ReCS(1,2) = 0.0D+00
     ImCS(1,2) = 0.0D+00

     ReCS(2,1) = 0.0D+00
     ImCS(2,1) = 0.0D+00

     ReCS(2,2) = auxCS*ReCTPER
     ImCS(2,2) = auxCS*ImCTPER

     ! Point d'impact apres un parcours interne: ##################################

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
     auxDELTA3(FACEimp) = -DDDD

     do IFACE = 0, 7
        if( auxDELTA3(IFACE).gt.EPS12 .and. auxDELTA3(IFACE).lt.DELTA3 )then
           DELTA3 = auxDELTA3(IFACE)
           auxFACEimp = IFACE
        endif
     enddo

     Ximp = auxXimp + DELTA3 * ALPHAe
     Yimp = auxYimp + DELTA3 * BBETAe
     Zimp = auxZimp + DELTA3 * GAMMAe
     FACEimp = auxFACEimp

     ! Mise a jour de l'attenuation: ##############################################

     DATT = DATT * dexp( - DELTA3 * KKKK * GRANDK * COSANGLEHETERO )

     ! Cosinus directeurs de la normale et de l'orthonormale au plan d'incidence de
     ! ce nouveau dioptre atteint: ################################################

     call NORMALEM(ALPHAe,BBETAe,EPS12,GAMMAe,          &
          FACEimp,ALPHAf,BBETAf,GAMMAf,        &
          XNR,YNR,ZNR,XND,XTD,YND,YTD,ZND,ZTD)

     ! Angle de rotation A de la normale interne pour l'amener dans le repere du ##
     ! nouveau plan d'incidence atteiny, Nr -> Nd: ################################

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

     ! Cosinus directeurs des faisceaux reflechi et transmis (refracte'), apres un
     ! (ou plusieurs) trajet(s) interne(s) a la particule: ########################

     ! Coefficients de reflexion et de transmission de Fresnel: ###################

     call CHANG_GA(RN,RM,                                &
          ALPHAf,BBETAf,GAMMAf,FACEimp,         &
          ALPHAe,BBETAe,GAMMAe,COSTHETAi,       &
          ALPHAr,BBETAr,GAMMAr,                 &
          ALPHAt,BBETAt,GAMMAt,COSTHETAt,IREFL, &
          GRANDN,GRANDK,COSANGLEHETERO,         &
          ReCRPAR,ReCRPER,ReCTPAR,ReCTPER,      &
          ImCRPAR,ImCRPER,ImCTPAR,ImCTPER)

     ! Cas particulier de reflexion totale: #######################################

     if( IREFL.eq.1 ) goto 12

     ! Matrice d'amplitudes diffusees apres deux refractions: #####################

     auxCT = DATT*dsqrt(COSTHETAt/(GRANDN*COSTHETAi))

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

     ! Directions de polarisation du champ emergent: ##############################

     call NORMALE(ALPHAt,BBETAt,GAMMAt,XND,YND,ZND, &
          XNS,XTS,YNS,YTS,ZNS,ZTS)

     ! Cosinus directeurs de la normale au plan de diffusion: #####################

     call RCDIRPD(ALPHAi,BBETAi,GAMMAi,                   & 
          ALPHAt,BBETAt,GAMMAt,EPS12,XNI,YNI,ZNI, &
          XNP,YNP,ZNP)

     ! Angle de rotation B de la normale au plan emergent pour l'amener dans le ###
     ! reperer du plan de diffusion, Ns -> Np: ####################################

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

     ! Angle de rotation D de la normale Np au plan de diffusion et la normale Ni #
     ! au plan d'incidence: #######################################################
     !<TMP>
     !      ReAPARd = ReCSD(1,1)*ReAPARi + ReCSD(1,2)*ReAPERi
     !      ReAPERd = ReCSD(2,1)*ReAPARi + ReCSD(2,2)*ReAPERi
     !      ImAPARd = ImCSD(1,1)*ReAPARi + ImCSD(1,2)*ReAPERi
     !      ImAPERd = ImCSD(2,1)*ReAPARi + ImCSD(2,2)*ReAPERi

     ! IPARd et IPERd = composantes parallele et perpendiculaire
     !                  de l'intensite diffusee... ################################

     !      IPARd = ReAPARd*ReAPARd + ImAPARd*ImAPARd
     !      IPERd = ReAPERd*ReAPERd + ImAPERd*ImAPERd
     !      WRITE(6,*)'avant D=',IPARd+IPERd
     !</TMP>

     DCOSD = XNP * XNI + YNP * YNI + ZNP * ZNI
     DSIND = XNP * XTI + YNP * YTI + ZNP * ZTI
     !<TMP>
     !      WRITE(6,*)'DCOSD=',DCOSD
     !      WRITE(6,*)'DSIND=',DSIND
     !</TMP>

     ReCSRD(1,1) = DCOSD*ReCSD(1,1) - DSIND*ReCSD(1,2)
     ImCSRD(1,1) = DCOSD*ImCSD(1,1) - DSIND*ImCSD(1,2)

     ReCSRD(1,2) = DSIND*ReCSD(1,1) + DCOSD*ReCSD(1,2)
     ImCSRD(1,2) = DSIND*ImCSD(1,1) + DCOSD*ImCSD(1,2)

     ReCSRD(2,1) = DCOSD*ReCSD(2,1) - DSIND*ReCSD(2,2)
     ImCSRD(2,1) = DCOSD*ImCSD(2,1) - DSIND*ImCSD(2,2)

     ReCSRD(2,2) = DSIND*ReCSD(2,1) + DCOSD*ReCSD(2,2)
     ImCSRD(2,2) = DSIND*ImCSD(2,1) + DCOSD*ImCSD(2,2)

     ! APARd et APERd: composantes parallele et perpendiculaire
     !                 de l'amplitude diffusee... #################################

     ReAPARd = ReCSRD(1,1)*ReAPARi + ReCSRD(1,2)*ReAPERi
     ReAPERd = ReCSRD(2,1)*ReAPARi + ReCSRD(2,2)*ReAPERi
     ImAPARd = ImCSRD(1,1)*ReAPARi + ImCSRD(1,2)*ReAPERi
     ImAPERd = ImCSRD(2,1)*ReAPARi + ImCSRD(2,2)*ReAPERi

     ! IPARd et IPERd = composantes parallele et perpendiculaire
     !                  de l'intensite diffusee... ################################

     IPARd = ReAPARd*ReAPARd + ImAPARd*ImAPARd
     IPERd = ReAPERd*ReAPERd + ImAPERd*ImAPERd

     ! FDRT = flux d'energie diffusee par Ray-Tracing... ##########################

     FDRT = FDRT + IPARd + IPERd
     !<TMP>
     !      WRITE(6,*)'IPHINC=',IPHINC
     !      WRITE(6,*)'FDRTi=',IPARd+IPERd
     !      WRITE(6,*)'FDRT=',FDRT
     !      pause
     !</TMP>

     ! Angle de diffusion entre les faisceaux incident et diffuse': ###############

     COSGRANDTHETA = ALPHAi * ALPHAt + &
          BBETAi * BBETAt + &
          GAMMAi * GAMMAt

     !<CHG,4/10/2006>
     IF (COSGRANDTHETA .GT. 1.0D0)  COSGRANDTHETA = 1.0D0
     IF (COSGRANDTHETA .LT. -1.0D0) COSGRANDTHETA = -1.0D0
     !</CHG,4/10/2006>
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

     ! Matrice de diffusion brute (= non-normalisee) aux44GTGP, calculee a partir #
     ! des amplitudes diffusees par ray-tracing, AUX22GTGP: #######################

     call STOKES(ReCSRD,ImCSRD,aux44GTGP)

     do JP = 1, 4
        do KP = 1, 4
           aux44GT(JP,KP,IDGT) = aux44GT(JP,KP,IDGT) + aux44GTGP(JP,KP)
        enddo
     enddo

     ! Les etapes ci-dessus ne sont pas prises en compte lors du cas particulier de
     ! reflexion totale par une face interne: #####################################

12   continue

     ! Le faisceau reflechi a l'interieur de la particule devient le faisceau emer-
     ! gent a partir de cette face interne: #######################################

     ALPHAe  = ALPHAr
     BBETAe  = BBETAr
     GAMMAe  = GAMMAr

     ! Directions de la normale et de l'orthonormale au plan reflechi: ############

     call NORMALE(ALPHAe,BBETAe,GAMMAe,XND,YND,ZND, &
          XNR,XTR,YNR,YTR,ZNR,ZTR)

     ! Matrice d'amplitudes reflechies vers l'interieur de la particule: ##########

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

     ! APARd et APERd: composantes parallele et perpendiculaire
     !                 de l'amplitude diffusee... #################################

     ReAPARd = ReCS(1,1)*ReAPARi + ReCS(1,2)*ReAPERi
     ReAPERd = ReCS(2,1)*ReAPARi + ReCS(2,2)*ReAPERi
     ImAPARd = ImCS(1,1)*ReAPARi + ImCS(1,2)*ReAPERi
     ImAPERd = ImCS(2,1)*ReAPARi + ImCS(2,2)*ReAPERi

     ! IPARd et IPERd = composantes parallele et perpendiculaire
     !                  de l'intensite diffusee... ################################

     IPARd = ReAPARd * ReAPARd + ImAPARd * ImAPARd
     IPERd = ReAPERd * ReAPERd + ImAPERd * ImAPERd

     ! Critere d'arret pour les parcours internes: ################################

     if( ( ( IPARd+IPERd )/( IPARi+IPERi ) ).gt.(1.0D-06) ) goto 2020

     ! Nombre des photons diffuses: ###############################################

     NPHDIF = NPHDIF + 1

     ! Repartition angulaire de l'energie diffractee: #############################

     do IDGT = 1, NDGTFR
        !<TMP,21/01>
        EDFR(IDGT) = EDFR(IDGT) +                          &
             ( minEDFR(IDGT) +                     &
             ( auxSURFAPP - minSURFAPP ) *       &
             ( maxEDFR(IDGT) - minEDFR(IDGT) ) / &
             ( maxSURFAPP - minSURFAPP ) )
        !         EDFR(IDGT) = 0.0D+00
        !</TMP,21/01>
     enddo

     ! La fin de la troisieme boucle: les photon incidents ########################
  enddo

  ! Flux d'energie diffractee: #################################################

  FDFR = 0.0D+00
  do IDGT = 1, NDGTFR
     FDFR = FDFR + ATOTi * ATOTi * EDFR(IDGT)
  enddo

  ! Section efficaces totales (GO=RT+FR) de diffusion et d'extinction: #########

  CEXTGO = ( FDFR / FINC + 1.0D+00     ) * sumSURFAPP / NPHDIF
  CSCAGO = ( FDFR / FINC + FDRT / FINC ) * sumSURFAPP / NPHDIF

  ! Albedo pour une diffusion, ray-tracing seulement: ##########################

  PIZERO = FDRT / FINC

  ! Albedo pour une diffusion, total (eq.15, Mishchenko & Macke 1998): #########

  PIZERO = ( PIZERO + 1.0D+00 ) / 2.0D+00

  ! La contribution due a la diffraction n'intervient que dans les elements dia-
  ! gonaux de la matrice de diffusion: #########################################

  do IDGT = 1, NDGTFR
     aux44GT(1,1,IDGT) = aux44GT(1,1,IDGT) + EDFR(IDGT)
     aux44GT(2,2,IDGT) = aux44GT(2,2,IDGT) + EDFR(IDGT)
     aux44GT(3,3,IDGT) = aux44GT(3,3,IDGT) + EDFR(IDGT)
     aux44GT(4,4,IDGT) = aux44GT(4,4,IDGT) + EDFR(IDGT)
  enddo

  ! Normalisation des elements 'non-P11' de la matrice de diffusion: ###########

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

  ! Facteur d'asymetrie et normalisation de la fonction de phase: ##############

  NUMERGG = 0.0D+00
  DENOMGG = 0.0D+00
  do IDGT = 1, NDGT
     ANGLESOLIDE(IDGT) = 2.0D+00 * PI *                      &
          ( COSDGT(IDGT,1) - COSDGT(IDGT,2) )
     !<CHG, 12/02/08>
     aux44GT(1,1,IDGT) = aux44GT(1,1,IDGT) / ANGLESOLIDE(IDGT)
     NUMERGG = NUMERGG +                                     &
          aux44GT(1,1,IDGT) *                           &
          ( COSDGT(IDGT,1) * COSDGT(IDGT,1) -           &
          !                     COSDGT(IDGT,2) * COSDGT(IDGT,2) ) / 2.0D+00
          COSDGT(IDGT,2) * COSDGT(IDGT,2) ) * PI
     DENOMGG = DENOMGG +                           &
          aux44GT(1,1,IDGT) *                 &
          ANGLESOLIDE(IDGT)
     !                   ( COSDGT(IDGT,1) - COSDGT(IDGT,2) )
     !</CHG, 12/02/08>
  enddo
  GG = NUMERGG / DENOMGG

  do IDGT = 1, NDGT
     PGT(1,1,IDGT) = 4.0D+00 * PI * aux44GT(1,1,IDGT) / DENOMGG
  enddo

  ! Impression des parametres de sortie: #######################################

!!$      write(*,*) 
!!$      write(*,*) ' PARAMETRES_DE_SORTIE:'
!!$      write(*,'('' Nombre_total_des_photons_diffuses_= '',i9)') &
!!$                NPHDIF
!!$      write(*,'('' Flux_d_energie_incidente_= '',1e12.6)') &
!!$                FINC
!!$      write(*,'('' Flux_d_energie_diffusee_(RT)_= '',1e12.6)') &
!!$                FDRT
!!$      write(*,'('' Flux_d_energie_diffractee_= '',1e12.6)') &
!!$                FDFR
!!$      write(*,'('' Section_efficace_d_extinction_(1E-12_m2)_= '',1e12.6)') &
!!$                CEXTGO
!!$      write(*,'('' Section_efficace_de_diffusion_(1E-12_m2)_= '',1e12.6)') &
!!$                CSCAGO
!!$      write(*,'('' Albedo_pour_une_diffusion_= '',1e12.6)') &
!!$                PIZERO
!!$      write(*,'('' Facteur_d_asymetrie_= '',1e12.6)') &
!!$                GG

!!$      write(*,*)
!!$      write(*,*) ' Matrice_de_diffusion: '
!!$      write(*,'(a39,4a39,a21)') '  Thetad      P11         P12/P11      ', &
!!$                 'P13/P11      P14/P11      P21/P11      ', &
!!$                 'P22/P11      P23/P11      P24/P11      ', &
!!$                 'P31/P11      P32/P11      P33/P11      ', &
!!$                 'P34/P11      P41/P11      P42/P11      ', &
!!$                 'P43/P11      P44/P11 '
!!$      do IDGT = 1, NDGT
!!$         write(*,'(1x,f8.4,16(1x,1e12.6))')                      &
!!$                   ( DGT(IDGT,1) + DGT(IDGT,2) )/2.0D+00,        &
!!$                   ( ( PGT(JP,KP,IDGT), KP = 1, 4 ), JP = 1, 4 )
!!$      enddo
!!$      write(*,*) ''
!!$      write(*,*) '                          PHM OUT '
!!$      write(*,*) '############################################################'
!!$      write(*,*) ''

  do IDGT = 1, NDGT
     THETAOUT(IDGT) = ( DGT(IDGT,1) + DGT(IDGT,2) )/2.0D+00
  enddo

  PHM_RES(1) =  NPHDIF ! Nombre total de PHotons DIFfuses
  PHM_RES(2) =  FINC ! Flux_d_energie_incidente
  PHM_RES(3) =  FDRT ! Flux_d_energie_diffusee
  PHM_RES(4) =  FDFR ! Flux_d_energie_diffractee
  PHM_RES(5) =  CEXTGO ! Section_efficace_d_extinction
  PHM_RES(6) =  CSCAGO ! Section_efficace_de diffusion
  PHM_RES(7) =  PIZERO ! Albedo_pour_une_diffusion
  PHM_RES(8) =  GG     ! Facteur_d_asymetrie

  ! Commented by MC
  ! Fin de la deuxieme grande boucle: les particules d'interet #################
  ! 2000 continue
  ! Fin de la premiere grande boucle: les nombres d'onde d'interet #############
  ! 1000 continue

  DEALLOCATE(EDFR, maxEDFR, minEDFR)

end subroutine hm_phm2p4
