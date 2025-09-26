!a ############################################################################
!a
!a CODE 'RHM' (Rought Hexagonal Monocrystal):
!a CALCUL DES PROPRIETES OPTIQUES D'UNE PARTICULE DE GLACE
!
!  Version 1.1, voir <CHG,RHM1.1>
!  J'ai juste rajoute des tests sur la direction du photon incident, 
!  reflechie et transmis par rapport a la normale perturbee
!                     7 May 2007 
!  J'ai compare au calcul de Hess et ca marche!!!!!!
!
!  Base sur la version 2.4 (14/10/2006) du code PHM 
!a ############################################################################

subroutine hm_rhm1p1(NPHINC,lam, rn, rm, ll, rr, tilt, ndgt, &
     RHM_RES, thetaout, PGT)

  !        ************************************************************
  !        *       C.-LABONNOTE Laurent and BROGNIER Gerard           *
  !        *     Laboratoire d'Optique Atmospherique  UMR 8518        *
  !        *    Universite des Sciences et Technologies de Lille      *
  !        *         59655 Villeneuve d'Ascq Cedex - France           *
  !        *         E-mail: gerard.brogniez@univ-lille1.fr           *
  !        ************************************************************
  !a
  !a                          Octobre 2006
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

  REAL(KIND=8), PARAMETER                    :: ATOTi   = 1.00D+00                                 ! amplitude (totale) incidente;
  REAL(KIND=8), PARAMETER                    :: EPS12   = 1.00D-12                                 ! 'precision' de plusieurs calculs;


  !a Variables d'entree et de sortie: ###########################################

  INTEGER, INTENT(IN)                        :: NPHINC                                             ! Nombre de PHotons INCidents;
  REAL(KIND=8), INTENT(IN)                   :: LAM                                                ! wavelength (microns)                 
  REAL(KIND=8), INTENT(IN)                   :: LL                                                 ! longueur de la colonne hexagonale [1E-6 m];
  REAL(KIND=8), INTENT(IN)                   :: RR                                                 ! rayon de la colonne hexagonale [1E-6 m];
  REAL(KIND=8), INTENT(IN)                   :: RM                                                 ! parties imaginaire et reelle de
  REAL(KIND=8), INTENT(IN)                   :: RN                                                 ! l'indice de refraction de la glace;
  REAL(KIND=8), INTENT(IN)                   :: TILT                                               !
  INTEGER, INTENT(IN)                        :: NDGT                                               ! Nombre total des Domaines Grand-Theta;

  ! MC comments unused variables
  !   REAL(KIND=8), DIMENSION(30)                :: readLL                                             ! longueur de la colonne hexagonale [1E-6 m];
  !   REAL(KIND=8), DIMENSION(1643)              :: readRM                                             ! parties imaginaire et reelle de
  !   REAL(KIND=8), DIMENSION(1643)              :: readRN                                             ! l'indice de refraction de la glace;
  !   REAL(KIND=8), DIMENSION(1643)              :: readNU                                             ! nombre d'onde (=10000/ALAMB) [cm-1];
  !   REAL(KIND=8), DIMENSION(30)                :: readRR                                             ! rayon de la colonne hexagonale [1E-6 m];

  REAL(KIND=8), INTENT(OUT), DIMENSION(8)    :: RHM_RES                                            ! TABLE CONTAINING RESULTING SCALAR VALUES
  !  RHM_RES(1) =  NPHDIF : Nombre total de PHotons DIFfuses
  !  RHM_RES(2) =  FINC : Flux_d_energie_incidente
  !  RHM_RES(3) =  FDRT : Flux_d_energie_diffusee
  !  RHM_RES(4) =  FDFR : Flux_d_energie_diffractee
  !  RHM_RES(5) =  CEXTGO : Section_efficace_d_extinction
  !  RHM_RES(6) =  CSCAGO : Section_efficace_de diffusion
  !  RHM_RES(7) =  PIZERO : Albedo_pour_une_diffusion
  !  RHM_RES(8) =  GG     : Facteur_d_asymetrie
  REAL(KIND=8), intent(out)                 :: PGT(4,4,NDGT)                                       ! matrice de diffusion, comme fonction de 
  ! l'angle GRANDTHETA de diffusion; 
  REAL(KIND=8), intent(out)                 :: thetaout(NDGT)

  !a Variables intermediaires, fonctions et constantes: #########################

  REAL(KIND=8)                               :: CSCAGO                                             ! sec. effic. totale de diffusion [1E-12 m2];
  REAL(KIND=8)                               :: CEXTGO                                             ! sec. effic. totale d'extinction [1E-12 m2];
  REAL(KIND=8)                               :: FDFR                                               ! flux d'energie diffractee par la particule;
  REAL(KIND=8)                               :: FINC                                               ! flux d'energie incidente sur la particule;
  REAL(KIND=8)                               :: FDRT                                               ! flux d'energie diffusee (RT) par la part.;
  REAL(KIND=8)                               :: GG                                                 ! facteur d'asymetrie, total [adim.];

  REAL(KIND=8)                               :: PIZERO                                             ! albedo pour une diffusion [adim.];

  INTEGER                                    :: NPHDIF                                             ! Nombre total de PHotons DIFfuses;

  !   CHARACTER(LEN=28)                          :: SORTIE 

  REAL(KIND=8)                               :: NU                                                 ! nombre d'onde (=10000/ALAMB) [cm-1];
  REAL(KIND=8)                               :: RESOLGT                                            ! resolution angulaire de la matrice de
  ! diffusion, en GRANDTHETA [degres];
  INTEGER                                    :: NDGTFR                                             ! NDGT 'difFRaction' (resolution RESOLGT);

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
  REAL(KIND=8)                               :: anglesolide(NDGT)                                        ! angles solides elementaires
  ! (GRANDTHETA), integres sur GRANDPHI;
  REAL(KIND=8), DIMENSION(4,4)               :: aux44GTGP                                          ! variable auxiliaire (sortie de la routine STOKES);
  REAL(KIND=8)                               :: aux44GT(4,4,NDGT)                                            ! variable auxiliaire (matrice de diffusion
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
  REAL(KIND=8)                               :: COSDGT(NDGT,2)                                             ! cosinus de DGT(IDGT,1 et 2);
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
  REAL(KIND=8), allocatable                  :: EDFR(:)                                               ! distr. angul. de l'energie diffractee;
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

  !<CHG,RHM>
  INTEGER                                    :: indic
  REAL(KIND=8), DIMENSION(0:7)               :: THETAf
  REAL(KIND=8), DIMENSION(0:7)               :: PHIf
  REAL(KIND=8), DIMENSION(0:7)               :: ALPHAf_nopert                                      ! 3eme cosinus directeur, faces;
  REAL(KIND=8), DIMENSION(0:7)               :: BBETAf_nopert                                      ! 3eme cosinus directeur, faces;
  REAL(KIND=8), DIMENSION(0:7)               :: GAMMAf_nopert                                      ! 3eme cosinus directeur, faces;
  REAL(KIND=8)                               :: PSCAL
  !</CHG,RHM>

  REAL(KIND=8), PARAMETER                    :: PI     = 3.1415926535897932384626433832795029D+00  ! = constante PI;
  REAL(KIND=8), PARAMETER                    :: PI2    = PI / 2.0D+00                              ! = PI / 2;
  REAL(KIND=8), PARAMETER                    :: PI3    = PI / 3.0D+00                              ! = PI / 3;
  REAL(KIND=8), PARAMETER                    :: PI6    = PI / 6.0D+00                              ! = PI / 6;
  REAL(KIND=8), PARAMETER                    :: PI180  = PI / 180.0D+00                            ! = PI / 180;
  REAL(KIND=8), PARAMETER                    :: R32    = 0.8660254037844385965883020617184229D+00  ! = sqrt( 3 )/2;

  !===================================================

  !---------------------
  !  modif MC to allow for NDGT control as an input argument 
  RESOLGT = 180.0D0 / (NDGT - 1)
  NDGTFR  = idnint( 1.0D+00 + 90.0D+00 / RESOLGT )

  !--- Compute wavenumber
  NU = 10000.0D0 / LAM ! lam in microns

  allocate(EDFR(NDGTFR),maxEDFR(NDGTFR),minEDFR(NDGTFR))

  KKKK = 2.0D+00 * PI / ( 1.0D+04 / NU )
  DDDD = dsqrt( LL * LL + 4.0D+00 * RR * RR )

  LL2 = LL / 2.0D+00
  AX = RR * R32

  !a Ouverture du fichier de sortie: ############################################

!!$      write(*,*) ''
!!$      write(*,*) '############################################################'
!!$      write(*,*) '                          RHM '
!!$      write(*,*) ''
!!$      write(*,*)
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
!!$      write(*,'('' Tilt angle (degres)                   = '',1e12.6)') &
!!$                TILT
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

  !a Initialisation de plusieurs variables: #####################################

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

  !a Angle de diffusion GRANDTHETA: #############################################
  !a
  !a    DGT(IDGT,1) = angle inferieur du domaine;
  !a    DGT(IDGT,2) = angle superieur du domaine;

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

  !a Les sommets de la particule qui sont eclaires: #############################
  !a
  !a    SE(1,I), SE(2,I) et SE(3,I): ordonnees X, Y et Z du I-eme sommet eclaire'

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

  !a Cosinus directeurs des faces de la particule: ##############################

  !<CHG,RHM>
  !--- Compute THETA and PHI attached to each face's cosine
  !    director in the crystal basis
  do IFACE = 0, 5
     THETAf(IFACE) = 0.0D0
     PHIf(IFACE)   = DBLE(IFACE) * PI3
     ALPHAf_nopert(IFACE) = dcos( dfloat( IFACE ) * PI3 )
     BBETAf_nopert(IFACE) = dsin( dfloat( IFACE ) * PI3 )
     GAMMAf_nopert(IFACE) = 0.0D+00
     !         ALPHAf(IFACE) = dcos( dfloat( IFACE ) * PI3 )
     !         BBETAf(IFACE) = dsin( dfloat( IFACE ) * PI3 )
     !         GAMMAf(IFACE) = 0.0D+00
  enddo
  do IFACE = 6, 7
     IF (IFACE .EQ. 6) THEN
        THETAf(IFACE) = PI2
        PHIf(IFACE)   = 0.0D0
     ELSE
        THETAf(IFACE) = -PI2
        PHIf(IFACE)   = 0.0D0
     ENDIF
     ALPHAf_nopert(IFACE) = 0.0D+00
     BBETAf_nopert(IFACE) = 0.0D+00
     GAMMAf_nopert(IFACE) = dcos( dfloat( IFACE - 6 ) * PI )
     !         ALPHAf(IFACE) = 0.0D+00
     !         BBETAf(IFACE) = 0.0D+00
     !         GAMMAf(IFACE) = dcos( dfloat( IFACE - 6 ) * PI )
  enddo
  !</CHG,RHM>

  !a Surfaces apparentes extremes et leurs orientations: ########################

  call EXTREMES(RR,LL,PI180,R32,maxSURFAPP,minSURFAPP,       &
       maxETA,minETA,maxPHI,minPHI,maxZETA,minZETA)

  !a Contributions extremes associees a la diffraction: #########################

  auxNPNPT = 8 
  if( minSURFAPP.gt.(1.0D+06) ) auxNPNPT = 10

  call DIFFRACT(NU,KKKK,LL,RR,SE,maxSURFAPP,maxETA,maxZETA,maxPHI, &
       NDGT,NDGTFR,RESOLGT,COSDGT,                        &
       EPS12,PI,PI180,PI3,auxNPNPT,maxEDFR)

  call DIFFRACT(NU,KKKK,LL,RR,SE,minSURFAPP,minETA,minZETA,minPHI, &
       NDGT,NDGTFR,RESOLGT,COSDGT,                        &
       EPS12,PI,PI180,PI3,auxNPNPT,minEDFR)

  !<CHG,version 2.4>
  !  Initialisation de la fonction random
  CALL RANDOM_SEED
  !</CHG,version 2.4>

  !a La troisieme grande boucle: les photons incidents ##########################

  do IPHINC = 1, NPHINC

     !<TMP>
     !        WRITE(6,*)'IPHINC=',IPHINC
     !</TMP>

     !a Initialisation de l'attenuation a l'interieur de la particule: #############

     DATT = 1.0D+00

     !<CHG,RHM1.1>
2400 CONTINUE
     !</CHG,RHM1.1>

     !a Choix aleatoire de l'orientation du faisceau incident: #####################
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
     mu = rnmu                        !parceque l on veut theta
     !entre 0 et 90 degres
     !</CHG,version 2.4>

     !a Section geometrique totale: ################################################

     auxSURFAPP = surfapp(RR,LL,ETA,ZETA,PHI,R32)

     sumSURFAPP = sumSURFAPP + auxSURFAPP

     !a Cosinus directeurs du faisceau incident: ###################################

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
     !</CHG,version 2.4>

     ALPHAi = - COSTHETA * COSPHI
     BBETAi = - COSTHETA * SINPHI
     GAMMAi = - SINTHETA

     !a Impact du photon incident sur une face externe: ############################

     FACEimp = -999

     call IMPACT(LL2,AX,COSPHI,COSTHETA,SINPHI,SINTHETA,EPS12,SE,R32, &
          Ximp,Yimp,Zimp,FACEimp)

     !a Verification des resultats d'IMPACT, concernant la face atteinte: ##########

     call TEST_IMPACT(NPHDIF,FACEimp,ETA,ZETA,PHI,PI180,EPS12)

     !a APARi et APERi: composantes parallele et perpendiculaire
     !a                 de l'amplitude incidente; polarisation aleatoire... ########

     call random_number( rnPOLAR )

     ReAPARi = ATOTi * dcos( 2.0D+00*PI*rnPOLAR )
     ReAPERi = ATOTi * dsin( 2.0D+00*PI*rnPOLAR )

     !a IPARi et IPERi = composantes parallele et perpendiculaire
     !a                  de l'intensite incidente... ###############################

     IPARi = ReAPARi*ReAPARi
     IPERi = ReAPERi*ReAPERi

     !a FINC = flux d'energie incidente (intensite totale incidente)... ############

     FINC = FINC + IPARi + IPERi

     !a Cosinus directeurs des faisceaux reflechi et transmis: #####################
     !a Coefficients de reflexion et de transmission de Fresnel: ###################

     !<CHG,RHM>
     !---- Dans un premier temps calcul des cosinus directeur des faces en 
     !     considerant une variation aleatoire de la normale dans un cone d angle 
     !     au sommet egal a 2*TILT.

     CALL COS_PERTUB(TILT,PI,PI180,FACEimp,THETAf,PHIf,ALPHAf,BBETAf,GAMMAf)
     !<TMP> OK ca marche !!!!
     !      WRITE(6,*)'produit scalaire entre 2 normales'
     !      DO iface = 0,7
     !        WRITE(6,*)iface,ALPHAf(iface)*ALPHAf_nopert(iface)+&
     !                  BBETAf(iface)*BBETAf_nopert(iface)+&
     !                  GAMMAf(iface)*GAMMAf_nopert(iface)
     !      ENDDO
     !      pause  
     !</CHG,RHM>

     !<CHG,RHM1.1>
     !---- Il faut absolument verifier que le rayon incident rentre ds le cristal
     !     par rapport a la normal perturbee.
     !     Pour cela je verifie que le produit scalaire du faisceau incident avec
     !     la normale perturbee est < 0 (car on est a l'exterieure du cristal)
     PSCAL = ALPHAi*ALPHAf(FACEimp) +                         &
          BBETAi*BBETAf(FACEimp) +                         &
          GAMMAi*GAMMAf(FACEimp)
     IF (PSCAL .GE. 0.0D0) THEN  !Je relance un autre rayon
        GOTO 2400
     ENDIF
     !</CHG,RHM1.1>

     call CHANG_AG(RN,RM,EPS12,                     &
          ALPHAf,BBETAf,GAMMAf,FACEimp,    &
          ALPHAi,BBETAi,GAMMAi,COSTHETAi,  &
          ALPHAr,BBETAr,GAMMAr,            &
          ALPHAt,BBETAt,GAMMAt,COSTHETAt,  &
          GRANDN,GRANDK,COSANGLEHETERO,    &
          ReCRPAR,ReCRPER,ReCTPAR,ReCTPER, &
          ImCRPAR,ImCRPER,ImCTPAR,ImCTPER)

     !<CHG,RHM1.1> 
     !---- Verifions que la direction du faisceau transmis ds le cristal va bien 
     !     ds le cristal!!
     !     Il suffit pour cela de verifier que le produit scalaire entre ce faisceau 
     !     transmis et la normale non perturbee soit < 0
     PSCAL = ALPHAt*ALPHAf_nopert(FACEimp) +                         &
          BBETAt*BBETAf_nopert(FACEimp) +                         &
          GAMMAt*GAMMAf_nopert(FACEimp)
     IF (PSCAL .GE. 0.0D0) THEN
        WRITE(6,*)'AG:',FACEimp,PSCAL
        pause
     ENDIF
     !</CHG,RHM1.1>

     !<TMP>
     !      IF (COSTHETAt*COSTHETAi .LT. 0.0D0) THEN
     !        WRITE(6,*)'FACEimp=',FACEimp
     !        WRITE(6,*)ALPHAf(FACEimp),BBETAf(FACEimp),GAMMAf(FACEimp)
     !        WRITE(6,*)ALPHAf_nopert(FACEimp),BBETAf_nopert(FACEimp),GAMMAf_nopert(FACEimp)
     !        WRITE(6,*)'COSTHETAt=',COSTHETAt
     !        WRITE(6,*)'COSTHETAi=',COSTHETAi
     !        WRITE(6,*)'PSCAL=',PSCAL
     !        WRITE(6,*)ALPHAi,BBETAi,GAMMAi
     !        WRITE(6,*)ALPHAr,BBETAr,GAMMAr
     !        WRITE(6,*)ALPHAt,BBETAt,GAMMAt
     !        pause
     !      ENDIF
     !</TMP>

     !a Matrice d'amplitude diffusee, par reflexion a la face externe: #############

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

     !a Matrice de diffusion brute (= non-normalisee) aux44GTGP, calculee a partir #
     !a des amplitudes diffusees par ray-tracing, AUX22GTGP: #######################

     call STOKES(ReCSRD,ImCSRD,aux44GTGP)

     do JP = 1, 4
        do KP = 1, 4
           aux44GT(JP,KP,IDGT) = aux44GT(JP,KP,IDGT) + aux44GTGP(JP,KP)
        enddo
     enddo

     !a Le faisceau transmis (refracte') vers l'interieur de la particule devient ##
     !a celui qui emerge, pour un parcours interne, a partir de la face interne: ###

     ALPHAe  = ALPHAt
     BBETAe  = BBETAt
     GAMMAe  = GAMMAt

     !a Cosinus directeurs de la normale et de
     !a l'orthonormale au plan d'incidence du dioptre atteint: #####################

     call NORMALEPI(ALPHAi,BBETAi,GAMMAi,                  &
          EPS12,FACEimp,PI,ALPHAf,BBETAf,GAMMAf, & 
          XNI,XTI,YNI,YTI,ZNI,ZTI)

     !a Cosinus directeurs de la normale et
     !a de l'orthonormale au plan de refraction: ###################################

     call NORMALE(ALPHAe,BBETAe,GAMMAe,XNI,YNI,ZNI, &
          XNR,XTR,YNR,YTR,ZNR,ZTR)

     !a Matrice d'amplitudes transmises vers l'interieur: ##########################

     auxCS = dsqrt(GRANDN*COSTHETAt/COSTHETAi)

     !<TMP>
     !      IF (auxCS .NE. auxCS) THEN
     !        WRITE(6,*)'COSTHETAi=',COSTHETAi
     !        WRITE(6,*)'GRANDN=',GRANDN
     !        WRITE(6,*)'COSTHETAt=',COSTHETAt
     !      ENDIF
     !</TMP>

     ReCS(1,1) = auxCS*ReCTPAR
     ImCS(1,1) = auxCS*ImCTPAR

     ReCS(1,2) = 0.0D+00
     ImCS(1,2) = 0.0D+00

     ReCS(2,1) = 0.0D+00
     ImCS(2,1) = 0.0D+00

     ReCS(2,2) = auxCS*ReCTPER
     ImCS(2,2) = auxCS*ImCTPER

     !a Point d'impact apres un parcours interne: ##################################

     !<TMP>
     !      indic = 0
     !</TMP>

2020 continue

     !<TMP>
     !      indic = indic + 1
     !</TMP>

     auxXimp = Ximp
     auxYimp = Yimp
     auxZimp = Zimp

     do IFACE = 0, 5
        !<CHG,RHM>
        !         auxauxDELTA3 = ALPHAe * ALPHAf(IFACE) + BBETAe * BBETAf(IFACE)
        auxauxDELTA3 = ALPHAe * ALPHAf_nopert(IFACE) +                    &
             BBETAe * BBETAf_nopert(IFACE)
        if( dabs( auxauxDELTA3 ).gt.EPS12 )then
           auxDELTA3(IFACE) = ( AX -                                     &
                !                                  auxXimp * ALPHAf(IFACE) -                &
                !                                  auxYimp * BBETAf(IFACE) ) / auxauxDELTA3
                auxXimp * ALPHAf_nopert(IFACE) -         &
                auxYimp * BBETAf_nopert(IFACE) )         &
                / auxauxDELTA3
           !</CHG,RHM>
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

     !a Mise a jour de l'attenuation: ##############################################

     DATT = DATT * dexp( - DELTA3 * KKKK * GRANDK * COSANGLEHETERO )

     !a Cosinus directeurs de la normale et de l'orthonormale au plan d'incidence de
     !a ce nouveau dioptre atteint: ################################################

     !<CHG,RHM>
2500 CONTINUE
     !---- Dans un premier temps calcul des cosinus directeur des faces en 
     !     considerant une variation aleatoire de la normale dans un cone d angle 
     !     au sommet egal a 2*TILT.

     CALL COS_PERTUB(TILT,PI,PI180,FACEimp,THETAf,PHIf,ALPHAf,BBETAf,GAMMAf)

     !<CHG,RHM1.1>
     !---- Il faut absolument verifier que le rayon incident prend une direction 
     !     qui resorte du cristal par rapport a la normal perturbee.
     !     Pour cela je verifie que le produit scalaire du faisceau incident avec
     !     la normale perturbee est > 0 (car on est ds le cristal)
     PSCAL = ALPHAe*ALPHAf(FACEimp) +                         &
          BBETAe*BBETAf(FACEimp) +                         &
          GAMMAe*GAMMAf(FACEimp)
     IF (PSCAL .LE. 0.0D0) THEN
        GOTO 2500
     ENDIF
     !</CHG,RHM1.1>
     !</CHG,RHM>

     call NORMALEM(ALPHAe,BBETAe,EPS12,GAMMAe,          &
          FACEimp,ALPHAf,BBETAf,GAMMAf,        &
          XNR,YNR,ZNR,XND,XTD,YND,YTD,ZND,ZTD)

     !a Angle de rotation A de la normale interne pour l'amener dans le repere du ##
     !a nouveau plan d'incidence atteiny, Nr -> Nd: ################################

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

     !<TMP>
     !      WRITE(6,*)FACEimp,'****costheta2=',ALPHAf(FACEimp)*ALPHAf_nopert(FACEimp) +      &
     !	                        BBETAf(FACEimp)*BBETAf_nopert(FACEimp) +           &
     !	                        GAMMAf(FACEimp)*GAMMAf_nopert(FACEimp)
     !      WRITE(6,*)'COSTHETAi2=',ALPHAe*ALPHAf_nopert(FACEimp) +      &
     !	                        BBETAe*BBETAf_nopert(FACEimp) +           &
     !	                        GAMMAe*GAMMAf_nopert(FACEimp)
     !</TMP>

     call CHANG_GA(RN,RM,                                &
          ALPHAf,BBETAf,GAMMAf,FACEimp,         &
          ALPHAe,BBETAe,GAMMAe,COSTHETAi,       &
          ALPHAr,BBETAr,GAMMAr,                 &
          ALPHAt,BBETAt,GAMMAt,COSTHETAt,IREFL, &
          GRANDN,GRANDK,COSANGLEHETERO,         &
          ReCRPAR,ReCRPER,ReCTPAR,ReCTPER,      &
          ImCRPAR,ImCRPER,ImCTPAR,ImCTPER)

     !<CHG,RHM1.1> 
     !---- Verifions que la direction du faisceau reflechi ds le cristal va bien 
     !     ds le cristal!!
     !     Il suffit pour cela de verifier que le produit scalaire entre ce faisceau 
     !     reflechi et la normale non perturbee soit < 0
     PSCAL = ALPHAr*ALPHAf_nopert(FACEimp) +                         &
          BBETAr*BBETAf_nopert(FACEimp) +                         &
          GAMMAr*GAMMAf_nopert(FACEimp)
     IF (PSCAL .GE. 0.0D0) THEN
        !        WRITE(6,*)'GA:',FACEimp,PSCAL
        GOTO 2500
     ENDIF
     !</CHG,RHM1.1>

     !a Cas particulier de reflexion totale: #######################################

     if( IREFL.eq.1 ) goto 12

     !a Matrice d'amplitudes diffusees apres deux refractions: #####################

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

     !a Matrice de diffusion brute (= non-normalisee) aux44GTGP, calculee a partir #
     !a des amplitudes diffusees par ray-tracing, AUX22GTGP: #######################

     call STOKES(ReCSRD,ImCSRD,aux44GTGP)

     !<TMP>
     !      IF (aux44GTGP(1,1) .NE. aux44GTGP(1,1)) THEN
     !      WRITE(6,*)'IPHINC=',IPHINC,indic
     !      WRITE(6,*)'ReCSRD=',ReCSRD(:,:)
     !      WRITE(6,*)'ImCSRD=',ImCSRD(:,:)
     !      WRITE(6,*)'ReCSD=',ReCSD(:,:)
     !      WRITE(6,*)'ReCS=',ReCS(:,:)
     !      WRITE(6,*)'ReCT=',ReCT(:,:)
     !      WRITE(6,*)'ReCSA=',ReCSA(:,:)
     !      WRITE(6,*)'DCOSA=',DCOSA
     !      WRITE(6,*)'auxCT=',auxCT
     !      WRITE(6,*)'auxCS=',auxCS
     !      WRITE(6,*)'ReCTPAR=',ReCTPAR      
     !      WRITE(6,*)'COSTHETAt=',COSTHETAt
     !      WRITE(6,*)'DATT=',DATT
     !      WRITE(6,*)'ReCTPAR=',ReCTPAR
     !      WRITE(6,*)'ImCTPAR=',ImCTPAR
     !      WRITE(6,*)'ReCTPER=',ReCTPER
     !      WRITE(6,*)'ImCTPER=',ImCTPER
     !      auxCT = DATT*dsqrt(COSTHETAt/(GRANDN*COSTHETAi))
     !      WRITE(6,*)'ReCS=',ReCS(:,:)
     !      WRITE(6,*)'ImCS=',ImCS(:,:)
     !      WRITE(6,*)'GRANDN=',GRANDN
     !      WRITE(6,*)'COSTHETAi=',COSTHETAi
     !      pause
     !      ENDIF
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

     !<TMP> 
     !      IF (ReCS(1,1) .NE. ReCS(1,1)) THEN
     !        WRITE(6,*)'IPHINC=',IPHINC
     !        WRITE(6,*)'ReCS=',ReCS(:,:)
     !        WRITE(6,*)'ReCR=',ReCR(:,:)
     !        WRITE(6,*)'****************'
     !      ENDIF
     !</TMP>

     !a APARd et APERd: composantes parallele et perpendiculaire
     !a                 de l'amplitude diffusee... #################################

     ReAPARd = ReCS(1,1)*ReAPARi + ReCS(1,2)*ReAPERi
     ReAPERd = ReCS(2,1)*ReAPARi + ReCS(2,2)*ReAPERi
     ImAPARd = ImCS(1,1)*ReAPARi + ImCS(1,2)*ReAPERi
     ImAPERd = ImCS(2,1)*ReAPARi + ImCS(2,2)*ReAPERi

     !a IPARd et IPERd = composantes parallele et perpendiculaire
     !a                  de l'intensite diffusee... ################################

     IPARd = ReAPARd * ReAPARd + ImAPARd * ImAPARd
     IPERd = ReAPERd * ReAPERd + ImAPERd * ImAPERd

     !<TMP>
     !      IF ((IPARd .EQ. NaNQ
     !      WRITE(6,*)'IPHINC=',IPHINC
     !	  WRITE(6,*)'IPARd=',IPARd
     !	  WRITE(6,*)'IPERd=',IPERd
     !	  WRITE(6,*)'aux44GT=',aux44GT(1,1,10:15)
     !	  pause
     !</TMP>

     !a Critere d'arret pour les parcours internes: ################################

     if( ( ( IPARd+IPERd )/( IPARi+IPERi ) ).gt.(1.0D-06) ) goto 2020

     !a Nombre des photons diffuses: ###############################################

     NPHDIF = NPHDIF + 1

     !a Repartition angulaire de l'energie diffractee: #############################

     do IDGT = 1, NDGTFR
        !<TMP>
        EDFR(IDGT) = EDFR(IDGT) +                          &
             ( minEDFR(IDGT) +                     &
             ( auxSURFAPP - minSURFAPP ) *       &
             ( maxEDFR(IDGT) - minEDFR(IDGT) ) / &
             ( maxSURFAPP - minSURFAPP ) )
        !          EDFR(IDGT) = 0.0D+00
        !</TMP>
     enddo

     !a La fin de la troisieme boucle: les photon incidents ########################
  ENDDO

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
!!$      write(*,*) ' thetad      ',              &
!!$                  'P11        P12/P11      P13/P11      P14/P11      ', &
!!$                  'P21/P11      P22/P11      P23/P11      P24/P11      ', &
!!$                  'P31/P11      P32/P11      P33/P11      P34/P11      ', &
!!$                  'P41/P11      P42/P11      P43/P11      P44/P11 '
!!$      do IDGT = 1, NDGT
!!$         write(*,'(1x,f8.4,16(1x,1e12.6))')                      &
!!$                   ( DGT(IDGT,1) + DGT(IDGT,2) )/2.0D+00,        &
!!$                   ( ( PGT(JP,KP,IDGT), KP = 1, 4 ), JP = 1, 4 )
!!$      enddo

!!$      write(*,*) ''
!!$      write(*,*) '                          RHM OUT '
!!$      write(*,*) '############################################################'
!!$      write(*,*) ''

  do IDGT = 1, NDGT
     thetaout(IDGT) = (DGT(IDGT,1) + DGT(IDGT,2) )/2.0D+00
  enddo

  RHM_RES(1) =  NPHDIF ! Nombre total de PHotons DIFfuses
  RHM_RES(2) =  FINC   ! Flux_d_energie_incidente
  RHM_RES(3) =  FDRT   ! Flux_d_energie_diffusee
  RHM_RES(4) =  FDFR   ! Flux_d_energie_diffractee
  RHM_RES(5) =  CEXTGO ! Section_efficace_d_extinction
  RHM_RES(6) =  CSCAGO ! Section_efficace_de diffusion
  RHM_RES(7) =  PIZERO ! Albedo_pour_une_diffusion
  RHM_RES(8) =  GG     ! Facteur_d_asymetrie

  ! Commented by MC
  !!a Fin de la deuxieme grande boucle: les particules d'interet #################
  ! 2000 continue
  !!a Fin de la premiere grande boucle: les nombres d'onde d'interet #############
  ! 1000 continue

  deallocate(EDFR, maxEDFR, minEDFR)

end subroutine hm_rhm1p1

!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################


!   Suboutine qui permet de calculer les cosinus directeur des faces, 
!   perturbes. C'est a dire que la normale est choisie aleatoirement ds un 
!   cone d'angle au sommet 2*TILT centre sur la normale non perturbee. 

SUBROUTINE COS_PERTUB(TILT,PI,PI180,FACEimp,THETAf,PHIf,ALPHAf,BBETAf,GAMMAf)

  IMPLICIT NONE

  !--- Input Var
  INTEGER                                    :: FACEimp
  REAL(KIND=8)                               :: TILT
  REAL(KIND=8)                               :: PI, PI180
  REAL(KIND=8), DIMENSION(0:7)               :: THETAf
  REAL(KIND=8), DIMENSION(0:7)               :: PHIf

  !--- Local Var
  INTEGER                                    :: IFACE

  REAL(KIND=8)                               :: tet, bet
  REAL(KIND=8)                               :: rtet, rbet
  REAL(KIND=8)                               :: csteta, snteta
  REAL(KIND=8)                               :: csphi, snphi
  REAL(KIND=8)                               :: ALPHAf_pert
  REAL(KIND=8)                               :: BBETAf_pert
  REAL(KIND=8)                               :: GAMMAf_pert

  !--- Output Var
  REAL(KIND=8), DIMENSION(0:7)               :: ALPHAf
  REAL(KIND=8), DIMENSION(0:7)               :: BBETAf
  REAL(KIND=8), DIMENSION(0:7)               :: GAMMAf

  DO IFACE = 0,7

     !--- choose bet randomly between 0 and 2*pi
     !           tet                  0 and TILT
     CALL RANDOM_NUMBER(rtet)
     tet = TILT * rtet * PI180         !to have it in radian

     CALL RANDOM_NUMBER(rbet)
     bet = 2.0D0 * PI * rbet

     !--- compute new cosine director in the base attached to old cosine
     ALPHAf_pert = DSIN(tet) * DCOS(bet)
     BBETAf_pert = DSIN(tet) * DSIN(bet)
     GAMMAf_pert = DCOS(tet)

     !--- Rotation of the pertub cosine director from the basis attached
     !    to each normal to the crystal basis
     !    There are 2 rotations:
     !      1) rotation de THETA-PI/2 autour de Oy"
     !      2) rotation de -PHI autour de oz'
     csteta = DCOS(-THETAf(IFACE)+PI/2.0D0)
     snteta = DSIN(-THETAf(IFACE)+PI/2.0D0)
     csphi  = DCOS(-PHIf(IFACE))
     snphi  = DSIN(-PHIf(IFACE))

     ALPHAf(IFACE) =  ALPHAf_pert*csteta*csphi +                              &
          GAMMAf_pert*snteta*csphi +                              &
          BBETAf_pert*snphi
     BBETAf(IFACE) = -ALPHAf_pert*csteta*snphi -                              &
          GAMMAf_pert*snteta*snphi +                              &
          BBETAf_pert*csphi
     GAMMAf(IFACE) = -ALPHAf_pert*snteta + GAMMAf_pert*csteta

  ENDDO !IFACE

  RETURN

END SUBROUTINE COS_PERTUB

!a ############################################################################
