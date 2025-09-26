
!a Fonction SURFAPP:
!a Calcul de la surface apparente totale de la colonne / plaquette hexagonale, 
!a associee a une orientation donnee (ETA,PHI,ZETA) du faisceau incident.

function SURFAPP(RR,LL,ETA,ZETA,PHI,R32)

  IMPLICIT NONE

  REAL(KIND=8)                               :: COSPHI                                         ! cosinue de l'angle PHI;
  REAL(KIND=8)                               :: COSTHETA                                       ! cosinus de l'angle THETA;
  REAL(KIND=8)                               :: ETA                                            ! angle entre le faisceau incident et l'axe OZ' [adim.];
  REAL(KIND=8)                               :: LL                                             ! longueur de la colonne hexagonale [1E-6 m];
  REAL(KIND=8)                               :: PHI                                            ! angle entre la projection OXp du faisceau incident et l'axe OX [adim.];
  REAL(KIND=8)                               :: R32                                            ! = sqrt( 3 )/2;
  REAL(KIND=8)                               :: RR                                             ! rayon de la colonne hexagonale [1E-6 m];
  REAL(KIND=8)                               :: SINTHETA                                       ! sinus de l'angle THETA;
  REAL(KIND=8)                               :: ZETA                                           ! angle entre les axes OZ et OZ' [adim.];

  REAL(KIND=8)                               :: SURFAPP                                        ! fonction interne 'surface apparente de la colonne hexagonale';

  COSPHI   = dcos( PHI )
  SINTHETA = dcos( ZETA ) * dcos( ETA )
  COSTHETA = dsqrt( 1.0D+00 - SINTHETA * SINTHETA )
  SURFAPP = 2.0D+00 * RR * LL * ( COSTHETA * COSPHI ) + &
       3.0D+00 * R32 * RR * RR * ( SINTHETA )

  return
end function SURFAPP

! ############################################################################

subroutine DMRMUL(A,B,R,N,M,L)                         

  IMPLICIT NONE 

  REAL(KIND=8), DIMENSION(*)                 :: A                                              ! matrice reelle, N lignes et M colonnes
  REAL(KIND=8), DIMENSION(*)                 :: B                                              ! matrice reelle, M lignes et L colonnes
  REAL(KIND=8), DIMENSION(*)                 :: R                                              ! matrice reelle, N lignes et L colonnes

  INTEGER                                    :: I,IB,IK,IR,J,JI,K,L,M,N

  IR = 0                                                                      
  IK = - M                                                                     
  do 11, K = 1, L                                                                
     IK = IK + M                                                                   
     do 12, J = 1, N                                                                
        IR = IR + 1                                                                   
        JI = J - N                                                                    
        IB = IK                                                                     
        R(IR) = 0.0D+00        
        do 13, I = 1, M                                                                
           JI = JI + N                                                                   
           IB = IB + 1                                                                   
           R(IR) = R(IR) + A(JI) * B(IB)                                              
13      enddo
12   enddo
11 enddo

  return                                                                    
end subroutine DMRMUL

! ############################################################################

!a Integration angulaire par la methode de Romberg

subroutine HM_INTANG(A,B,F,L,NN,RED)

  IMPLICIT NONE

  REAL(KIND=8)                               :: A
  REAL(KIND=8)                               :: B
  REAL(KIND=8), DIMENSION(0:2048)            :: F
  REAL(KIND=8)                               :: FF
  REAL(KIND=8)                               :: HI
  REAL(KIND=8)                               :: RED
  REAL(KIND=8), DIMENSION(0:10,0:10)         :: T
  REAL(KIND=8)                               :: TT

  INTEGER                                    :: I,IP,IP1,J,K,K1,L,NN

  IP = 1
  HI = ( B - A ) / dfloat( L )
  FF = ( F(0) + F(L) ) / 2.0D+00
  do I = NN, 0, -1
     TT = 0.0D+00
     do J = IP, L - IP, IP
        TT = TT + F(J)
     enddo
     T(I,0) = ( TT + FF ) * HI
     HI = HI * 2.0D+00
     IP = IP * 2
  enddo
  IP = 4
  IP1 = IP - 1
  do K = 1, NN
     K1 = K - 1
     do I = 0, NN - K
        T(I,K) = ( dfloat( IP ) * T(I+1,K1) - T(I,K1) ) / dfloat( IP1 )
     enddo
     IP = IP * 4
     IP1 = IP - 1
  enddo
  RED = T(0,NN)

  return
end subroutine HM_INTANG



! ############################################################################

! CALCUL DES COSINUS DIRECTEURS DE LA NORMALE ET DE L'ORTHONORMALE
! AU PLAN D'INCIDENCE D'IMPACT
! (XNI,YNI,ZNI): COSINUS DIRECTEURS DE LA NORMALE AU PLAN D'INCIDENCE
! (XTI,YTI,ZTI): COSINUS DIRECTEURS DE L'ORTHONORMALE AU PLAN D'INCIDENCE

      subroutine NORMALEPI(ALPHAi,BBETAi,GAMMAi,EPS,IFACE,PI,ALPHAf,BBETAf,GAMMAf,XNI,XTI,YNI,YTI,ZNI,ZTI) 

      IMPLICIT NONE

      REAL(KIND=8)                               :: PI                                             ! = constante PI;
      REAL(KIND=8), DIMENSION(0:7)               :: ALPHAf                                         ! 1er cosinus directeur, faces de la colonne;
      REAL(KIND=8)                               :: ALPHAi                                         ! 1er cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: BBETAi                                         ! 2eme cosinus directeur, faisceau incident;
      REAL(KIND=8), DIMENSION(0:7)               :: BBETAf                                         ! 2eme cosinus directeur, faces;
      REAL(KIND=8)                               :: DN                                             ! variable auxiliaire;
      REAL(KIND=8)                               :: EPS                                            ! 'precision' de plusieurs calculs;
      REAL(KIND=8), DIMENSION(0:7)               :: GAMMAf                                         ! 3eme cosinus directeur, faces;
      REAL(KIND=8)                               :: GAMMAi                                         ! 3eme cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: PHI                                            ! angle entre la projection OXp du faisceau
      REAL(KIND=8)                               :: SNPH                                           ! variable auxiliaire;
      REAL(KIND=8)                               :: XNI,YNI,ZNI                                    ! cosinus directeurs (normale); 
      REAL(KIND=8)                               :: XTI,YTI,ZTI                                    ! cosinus directeurs (orthonormale);

      INTEGER                                    :: IFACE                                          ! compteur des faces de la particule;

      XNI = BBETAi * GAMMAf(IFACE) - GAMMAi * BBETAf(IFACE)
      YNI = GAMMAi * ALPHAf(IFACE) - ALPHAi * GAMMAf(IFACE)
      ZNI = ALPHAi * BBETAf(IFACE) - BBETAi * ALPHAf(IFACE)
      if( XNI.lt.(0.0D+00) )then 
          XNI = - XNI
          YNI = - YNI
          ZNI = - ZNI
      endif
      DN = XNI * XNI + YNI * YNI + ZNI * ZNI       
      if( DN.le.(0.0D+00) ) DN = 0.0D+00
      DN = dsqrt(DN)
      if( DN.le.EPS )then
          call random_number( SNPH )
          PHI = 2.0D+00 * PI * SNPH
          if( IFACE.le.5 )then
              XNI = -BBETAf(IFACE) * dsin( PHI ) 
              YNI = ALPHAf(IFACE) * dsin( PHI )
              ZNI = dcos( PHI )
          else
              XNI = dcos( PHI )
              YNI = dsin( PHI )
              ZNI = 0.0D+00
          endif
      else
          XNI = XNI / DN
          YNI = YNI / DN
          ZNI = ZNI / DN
      endif
      XTI = YNI * GAMMAi - ZNI * BBETAi
      YTI = ZNI * ALPHAi - XNI * GAMMAi
      ZTI = XNI * BBETAi - YNI * ALPHAi

      return
      end

!###################################################################################

!a Routine NORMALE:

! CALCUL DES COSINUS DIRECTEURS DE LA NORMALE ET DE L'ORTHONORMALE
! AU PLAN D'INCIDENCE D'IMPACT POUR LE CHAMP EMERGENT OU REFRACTE
! (XN,YN,ZN) : COSINUS DIRECTEURS DE LA NORMALE AU PLAN D'D'ERGENCE
! (XT,YT,ZT) : COSINUS DIRECTEURS DE L'ORTHONORMALE AU PLAN D'EMERGENCE
! (XNI,YNI,ZNI) :COSINUS DIRECTEURS DE LA NORMALE AU PLAN D'INCIDENCE
!                   T = N^W

      subroutine NORMALE(ALP,BET,GAM,XNI,YNI,ZNI,XN,XT,YN,YT,ZN,ZT) 

      IMPLICIT NONE

      REAL(KIND=8)                               :: ALP,BET,GAM,XN,YN,ZN,XNI,YNI,ZNI,XT,YT,ZT

      XN = XNI
      YN = YNI
      ZN = ZNI
      XT = YN * GAM - ZN * BET
      YT = ZN * ALP - XN * GAM
      ZT = XN * BET - YN * ALP

      return
      end

! ############################################################################

! CALCUL DES COSINUS DIRECTEURS DE LA NORMALE ET DE L'ORTHONORMALE
! AU PLAN D'INCIDENCE D'IMPACT 
! (XNI,YNI,ZNI) : COSINUS DIRECTEURS DE LA NORMALE AU PLAN D'INCIDENCE
! (XN,YN,ZN) : COSINUS DIRECTEURS DE LA NORMALE AU PLAN D'EMERGENCE
! (XT,YT,ZT) : COSINUS DIRECTEURS DE L'ORTHONORMALE AU PLAN D'EMERGENCE

      subroutine NORMALEM(ALP,BET,EPS,GAM,IFACE,ALPHAf,BBETAf,GAMMAf,XNI,YNI,ZNI,XN,XT,YN,YT,ZN,ZT) 

      IMPLICIT NONE

      REAL(KIND=8)                               :: ALP,BET,GAM                                    ! 
      REAL(KIND=8), DIMENSION(0:7)               :: ALPHAf                                         ! 1er cosinus directeur, faces de la colonne;
      REAL(KIND=8), DIMENSION(0:7)               :: BBETAf                                         ! 2eme cosinus directeur, faces;
      REAL(KIND=8)                               :: DN                                             ! variable auxiliaire;
      REAL(KIND=8)                               :: EPS                                            ! 'precision' de plusieurs calculs
      REAL(KIND=8), DIMENSION(0:7)               :: GAMMAf                                         ! 3eme cosinus directeur, faces;
      REAL(KIND=8)                               :: PS                                             ! variable auxiliaire;
      REAL(KIND=8)                               :: XN,YN,ZN                                       ! 
      REAL(KIND=8)                               :: XNI,YNI,ZNI                                    ! 
      REAL(KIND=8)                               :: XT,YT,ZT                                       ! 

      INTEGER                                    :: IFACE                                          ! compteur des faces de la particule

      XN = BET * GAMMAf(IFACE) - GAM * BBETAf(IFACE)
      YN = GAM * ALPHAf(IFACE) - ALP * GAMMAf(IFACE)
      ZN = ALP * BBETAf(IFACE) - BET * ALPHAf(IFACE)
      PS = XNI * XN + YNI * YN + ZNI * ZN
      if( PS.lt.(0.0D+00) )then 
          XN = - XN
          YN = - YN
          ZN = - ZN
      endif
      DN = XN * XN + YN * YN + ZN * ZN
      if( DN.le.(0.0D+00) ) DN = 0.0D+00
      DN = dsqrt( DN )
      if( DN.le.EPS )then
          XN = XNI
          YN = YNI
          ZN = ZNI
      else
          XN = XN / DN
          YN = YN / DN
          ZN = ZN / DN
      endif
      XT = YN * GAM - ZN * BET
      YT = ZN * ALP - XN * GAM
      ZT = XN * BET - YN * ALP

      return
      end

! ############################################################################

! Cosinus directeurs de la normale au plan de diffusion.

      subroutine RCDIRPD(ALPHAi,BBETAi,GAMMAi,ALPD,BETD,GAMD,EPS,XN,YN,ZN,XP,YP,ZP)

      IMPLICIT NONE

      REAL(KIND=8)                               :: ALP,BET,GAM                                    !
      REAL(KIND=8)                               :: ALPD,BETD,GAMD                                 ! 
      REAL(KIND=8)                               :: ALPHAi                                         ! 1er cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: BBETAi                                         ! 2eme cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: GAMMAi                                         ! 3eme cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: DP                                             ! variable auxiliaire;
      REAL(KIND=8)                               :: EPS                                            ! 'precision' de plusieurs calculs;
      REAL(KIND=8)                               :: PS                                             ! variable auxliaire;
      REAL(KIND=8)                               :: XN,YN,ZN                                       !
      REAL(KIND=8)                               :: XP,YP,ZP                                       !

      XP = GAMMAi * BETD - BBETAi * GAMD
      YP = ALPHAi * GAMD - GAMMAi * ALPD
      ZP = BBETAi * ALPD - ALPHAi * BETD
      PS = XP * XN + YP * YN + ZP * ZN
      if( PS.lt.(0.0D+00) )then
          XP = - XP
          YP = - YP
          ZP = - ZP
      endif
      DP = XP * XP + YP * YP + ZP * ZP
      if( DP.lt.(0.0D+00) ) DP = 0.0D+00
      DP = dsqrt( DP )
      if( DP.lt.EPS )then
          XP = XN
          YP = YN
          ZP = ZN
      else
          XP = XP / DP
          YP = YP / DP
          ZP = ZP / DP
      endif
      XP = XP
      YP = YP
      ZP = ZP

      return
      end

! ############################################################################

!a Routine IMPACT:
!a
!a Identification de la position du point d'impact d'un photon incident sur la
!a surface de la colonne / plaquette hexagonale, et de la face en question.
!a
!a Deux systemes de coordonnees sont utilises: OXYZ, attache' a la particule,
!a et OX'Y'Z' (ici, OXincYincZinc) attache' au faisceau incident.

      subroutine IMPACT(LL2,AX,COSPHI,COSTHETA,SINPHI,SINTHETA,EPS,SE,R32,Ximp,Yimp,Zimp,FACEimp)

      IMPLICIT NONE
 

      REAL(KIND=8), DIMENSION(13)                :: AD                                             ! coefficients AD(K)=A(d);

      REAL(KIND=8)                               :: ALPHAi                                         ! 1er cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: AX                                             ! = rayon * sqrt( 3 ) / 2;
      REAL(KIND=8)                               :: Ayinc                                          ! valeur aleatoire sur l'axe OYinc; 
      REAL(KIND=8)                               :: Azinc                                          ! valeur aleatoire sur l'axe OZinc; 
      REAL(KIND=8)                               :: BBETAi                                         ! 2eme cosinus directeur, faisceau incident;

      REAL(KIND=8), DIMENSION(13)                :: BD                                             ! coefficients BD(K)=B(d);

      REAL(KIND=8)                               :: COSPHI                                         ! cosinus de l'angle PHI;
      REAL(KIND=8)                               :: COSTHETA                                       ! cosinus de l'angle THETA;

      REAL(KIND=8), DIMENSION(13)                :: D                                              ! coefficients D(K)=-C(d);

      REAL(KIND=8), DIMENSION(9:13)              :: E                                              ! ecarts (identification de la face);

      REAL(KIND=8)                               :: EE                                             ! ecart (abandon ou non d'un photon);
      REAL(KIND=8)                               :: EPS                                            ! 'precision' de plusieurs calculs;
      REAL(KIND=8)                               :: GAMMAi                                         ! 3eme cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: LL2                                            ! la moitie de la longueur de la colonne;

      REAL(KIND=8), DIMENSION(3,3)               :: M                                              ! matrice entre OXYZ et OXincYincZinc;
      REAL(KIND=8), DIMENSION(3,3)               :: MI                                             ! matrice Inverse de 'M';

      REAL(KIND=8)                               :: R32                                            ! = sqrt( 3 )/2;
      REAL(KIND=8)                               :: RHO                                            ! parametre 'rho' (conversion XYZarr>XYZimp);
      REAL(KIND=8)                               :: SINPHI                                         ! sinus de l'angle PHI;
      REAL(KIND=8)                               :: SINTHETA                                       ! sinus de l'angle THETA;

      REAL(KIND=8), DIMENSION(3,10)              :: SE                                             ! coordonnees XYZ des dix sommets eclaires;
      REAL(KIND=8), DIMENSION(3)                 :: auxSE                                          ! coordonnees XYZ d'un sommet eclaire';
      REAL(KIND=8), DIMENSION(3,10)              :: SEinc                                          ! coordonnees XincYincZinc des dix sommets;
      REAL(KIND=8), DIMENSION(3)                 :: auxSEinc                                       ! coordonnees XincYincZinc d'un sommet;

      REAL(KIND=8)                               :: MAXyinc                                        ! valeur maximale des Yinc=SEinc(2,I);
      REAL(KIND=8)                               :: MINyinc                                        ! valeur minimale des Yinc=SEinc(2,I);
      REAL(KIND=8)                               :: MAXzinc                                        ! valeur maximale des Zinc=SEinc(3,I);
      REAL(KIND=8)                               :: MINzinc                                        ! valeur minimale des Zinc=SEinc(3,I);

      REAL(KIND=8), DIMENSION(3)                 :: XYZincarr                                      ! positions Xincarr(=0), Yincarr, Zincarr;

      REAL(KIND=8)                               :: Yincarr                                        ! position Yincarr d'arrivee (OXincYincZinc);
      REAL(KIND=8)                               :: Zincarr                                        ! position Zincarr d'arrivee (OXincYincZinc);

      REAL(KIND=8), DIMENSION(3)                 :: XYZarr                                         ! positions Xarr, Yarr, Zarr;

      REAL(KIND=8)                               :: Xarr                                           ! position Xarr d'arrivee (OXYZ);
      REAL(KIND=8)                               :: Yarr                                           ! position Xarr d'arrivee (OXYZ);
      REAL(KIND=8)                               :: Zarr                                           ! position Xarr d'arrivee (OXYZ);

      REAL(KIND=8)                               :: Ximp                                           ! ordonnee Ximp de l'impact;
      REAL(KIND=8)                               :: Yimp                                           ! ordonnee Yimp de l'impact; 
      REAL(KIND=8)                               :: Zimp                                           ! ordonnee Zimp de l'impact;

      INTEGER                                    :: I, J, K                                        ! 
      INTEGER                                    :: FACEimp                                        ! face sur laquelle a lieu l'impact;


!a Matrice 'M' de changement entre les systemes de coordonnees OXYZ et OXincYincZinc: #########################

      M(1,1) = COSTHETA * COSPHI
      M(1,2) = COSTHETA * SINPHI
      M(1,3) = SINTHETA
      M(2,1) = - SINPHI
      M(2,2) = COSPHI
      M(2,3) = 0.0D+00
      M(3,1) = - SINTHETA * COSPHI
      M(3,2) = - SINTHETA * SINPHI
      M(3,3) = COSTHETA

!a Matrice Inverse de 'M', pour le retour
!a du systeme OXincYincZinc au systeme OXYZ de coordonnees: ###################

      MI(1,1) = M(1,1)
      MI(1,2) = M(2,1)
      MI(1,3) = M(3,1)
      MI(2,1) = M(1,2)
      MI(2,2) = M(2,2)
      MI(2,3) = M(3,2)
      MI(3,1) = M(1,3)
      MI(3,2) = M(2,3)
      MI(3,3) = M(3,3)

!a Cosinus directeurs du faisceau incident: ###################################

      ALPHAi = M(1,1)
      BBETAi = M(1,2)
      GAMMAi = M(1,3)

!a Projection de la colonne hexagonale
!a sur le plan perpendiculaire au faisceau incident: ##########################

      do I = 1, 10
         auxSE(1) = SE(1,I) ! coordonnee X du I-eme sommet
         auxSE(2) = SE(2,I) ! coordonnee Y du I-eme sommet
         auxSE(3) = SE(3,I) ! coordonnee Z du I-eme sommet
         call DMRMUL(M,auxSE,auxSEinc,3,3,1)
         SEinc(1,I) = 0.0D+00     ! coordonnee Xinc du I-eme sommet
         SEinc(2,I) = auxSEinc(2) ! coordonnee Yinc du I-eme sommet
         SEinc(3,I) = auxSEinc(3) ! coordonnee Zinc du I-eme sommet
      enddo

!a Les frontieres externes de la particule projetee: ##########################

      call ADBD(AD,BD,D,SEinc,1,2,1)
      call ADBD(AD,BD,D,SEinc,2,3,2)
      call ADBD(AD,BD,D,SEinc,3,4,3)
      call ADBD(AD,BD,D,SEinc,4,5,4)
      call ADBD(AD,BD,D,SEinc,5,6,5)
      call ADBD(AD,BD,D,SEinc,6,7,6)
      call ADBD(AD,BD,D,SEinc,7,8,7)
      call ADBD(AD,BD,D,SEinc,8,1,8)

!a Les frontieres internes de la particule projetee: ########################## 

      call ADBD(AD,BD,D,SEinc, 1, 9, 9)
      call ADBD(AD,BD,D,SEinc, 9,10,10)
      call ADBD(AD,BD,D,SEinc,10, 4,11)
      call ADBD(AD,BD,D,SEinc,10, 6,12)
      call ADBD(AD,BD,D,SEinc, 9, 7,13)

!a Limites extremes de la particule projetee: ################################# 

      MINyinc = SEinc(2,1)
      MAXyinc = MINyinc
      MINzinc = SEinc(3,1)
      MAXzinc = MINzinc
      do K = 2, 10
         if( SEinc(2,K).le.MINyinc )then
           MINyinc = SEinc(2,K)
         else
           if( SEinc(2,K).ge.MAXyinc ) MAXyinc = SEinc(2,K)
         endif
         if( SEinc(3,K).le.MINzinc )then
             MINzinc = SEinc(3,K)
         else
             if( SEinc(3,K).ge.MAXzinc ) MAXzinc = SEinc(3,K)
         endif
      enddo

!a Coordonnees d'arrivee du photon sur la particule projetee,
!a encore reperees dans le systeme de coordonnees OXincYincZinc: ##############

    7 continue
      call random_number( Ayinc )
      Yincarr = ( MINyinc + ( MAXyinc - MINyinc ) * Ayinc )
      call random_number( Azinc )
      Zincarr = ( MINzinc + ( MAXzinc - MINzinc ) * Azinc )

!a Abandon des photons qui ne touchent pas la particule: ######################

      do I = 1, 8
         EE = Zincarr * D(I) + Yincarr * AD(I) + BD(I)
         if( EE.gt.EPS ) goto 7
      enddo

!a Coordonnees d'arrivee du photon sur la particule projetee,
!a enfin dans le systeme de coordonnees de la particule, OXYZ: ################

      XYZincarr(1) = 0.0D+00
      XYZincarr(2) = Yincarr
      XYZincarr(3) = Zincarr
      call DMRMUL(MI,XYZincarr,XYZarr,3,3,1)
      Xarr = XYZarr(1)
      Yarr = XYZarr(2)
      Zarr = XYZarr(3)
      
!a Identification de la face atteinte: ########################################

      do K = 9, 13
         E(K) = Zincarr * D(K) + Yincarr * AD(K) + BD(K)
      enddo

!a Premiere possibilite: la face numerotee 6 ##################################

      if( E( 9).gt.EPS .and. &
          E(10).gt.EPS .and. &
          E(11).gt.EPS )then 
          FACEimp = 6
          RHO = ( Zarr - LL2 ) / GAMMAi
          Zimp = LL2
      else

!a Seconde possibilite: la face numerotee 0 ###################################

          if( E(10).lt.EPS .and. &
              E(12).lt.EPS .and. &
              E(13).gt.EPS )then 
              FACEimp = 0
              RHO = ( Xarr - AX ) / ALPHAi
          else

!a Troisieme possibilite: la face numerotee 1 #################################

              if( E( 9).lt.EPS .and. &
                  E(13).lt.EPS )then
                  FACEimp = 1
                  RHO = ( Xarr * 0.5D+00 + Yarr * R32 - AX ) / & 
                        ( ALPHAi * 0.5D+00 + BBETAi * R32 )
              else

!a Quatrieme possibilite: la face numerotee 5 #################################

                  if( E(11).lt.EPS .and. &
                      E(12).gt.EPS )then
                      FACEimp = 5
                      RHO = ( Xarr * 0.5D+00 - Yarr * R32 - AX ) / & 
                            ( ALPHAi * 0.5D+00 - BBETAi * R32 )
                  endif
              endif 
          endif
          Zimp = Zarr - RHO * GAMMAi
      endif
      Yimp = Yarr - RHO * BBETAi
      Ximp = Xarr - RHO * ALPHAi

      return
      end

! ############################################################################

! Routine ADBD:
!
! Calcul des coefficients de l'equation d'une droite, sous la forme:
!
!    Zf = Yf.A(d)/C(d) + B(d)/C(d) = - Yf.AD(k)/D(k) - BD(k)/D(k)
!
!    A(d) = Zinc(s=j) - Zinc(s=i)
!    B(d) = Zinc(s=i).Yinc(s=j) - Zinc(s=j).Yinc(s=i)
!    C(d) = Yinc(s=j) - Yinc(s=i)
!
!    Yinc(s) = coordonnee Y d'un sommet de la colonne hexagonale, dans 
!              le systeme OXincYincZinc attache au faisceau incident
!    Zinc(s) = coordonnee Z d'un sommet de la colonne hexagonale, dans 
!              le systeme OXincYincZinc attache au faisceau incident
!
! ou d=1...13 et s=1...10.

      subroutine ADBD(AD,BD,D,SEinc,I,J,K)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(13), INTENT(OUT)   :: AD                                             ! coefficients AD(K)=A(d);
      REAL(KIND=8), DIMENSION(13), INTENT(OUT)   :: BD                                             ! coefficients BD(K)=B(d);
      REAL(KIND=8), DIMENSION(13), INTENT(OUT)   :: D                                              ! coefficients D(K)=-C(d);

      REAL(KIND=8), DIMENSION(3,10), INTENT(IN)  :: SEinc                                          ! coordonnees XincYincZinc des dix sommets;

      INTEGER, INTENT(IN)                        :: I, J, K                                        ! 

      AD(K) = SEinc(3,J) - SEinc(3,I)
      BD(K) = SEinc(3,I) * SEinc(2,J) - SEinc(2,I) * SEinc(3,J)
      D(K)  = SEinc(2,I) - SEinc(2,J)

      end

! ############################################################################

! Routine EXTREMES:
!
! Identification des valeurs extremes maxSURFAPP et minSURFAPP de la surface
! apparente de la colonne / plaquette hexagonale, ainsi que leurs respectives
! orientations (maxETA,maxZETA) et (minETA,minZETA).

      subroutine EXTREMES(RR,LL,PI180,R32,maxSURFAPP,minSURFAPP,maxETA,minETA,maxPHI,minPHI,maxZETA,minZETA) 

      IMPLICIT NONE

      INTEGER                                    :: IETA                                           ! Indicateur de l'angle ETA;
      INTEGER                                    :: IPHI                                           ! Indicateur de l'angle PHI;
      INTEGER                                    :: IZETA                                          ! Indicateur de l'angle ZETA;

      REAL(KIND=8)                               :: ETA                                            ! angle entre le faisceau incident et l'axe OZ' [adim.];
      REAL(KIND=8)                               :: LL                                             ! longueur de la colonne hexagonale [1E-6 m];
      REAL(KIND=8)                               :: maxETA                                         ! valeur maximale de l'angle ETA;
      REAL(KIND=8)                               :: maxPHI                                         ! valeur maximale de l'angle PHI;
      REAL(KIND=8)                               :: maxSURFAPP                                     ! valeur max. de SURFAPP; [1E-12 m2];
      REAL(KIND=8)                               :: maxZETA                                        ! valeur maximale de l'angle ZETA;
      REAL(KIND=8)                               :: minETA                                         ! valeur minimale de l'angle ETA;
      REAL(KIND=8)                               :: minPHI                                         ! valeur minimale de l'angle PHI;
      REAL(KIND=8)                               :: minSURFAPP                                     ! valeur max. de SURFAPP; [1E-12 m2];
      REAL(KIND=8)                               :: minZETA                                        ! valeur minimale de l'angle ZETA;
      REAL(KIND=8)                               :: PHI                                            ! angle entre la projection OXp du faisceau incident et l'axe OX [adim.];
      REAL(KIND=8)                               :: PI180                                          ! = PI / 180;
      REAL(KIND=8)                               :: preSURFAPP                                     ! valeurs preliminaires de max/minSURFAPP;
      REAL(KIND=8)                               :: R32                                            ! = sqrt( 3 )/2;
      REAL(KIND=8)                               :: RPAS                                           ! resol. angulaire de l'orientation du faisceau incident [degres];
      REAL(KIND=8)                               :: RR                                             ! rayon de la colonne hexagonale [1E-6 m];
      REAL(KIND=8)                               :: ZETA                                           ! angle entre les axes OZ et OZ' [adim.];

      REAL(KIND=8)                               :: SURFAPP                                        ! fonction interne 'surface apparente de la colonne hexagonale';

      RPAS = 1.0D+00

! Valeur initiale de 'minSURFAPP': la surface TOTALE des huit faces ##########

      minSURFAPP = 3.0D+00 * RR * ( 2.0D+00 * LL + RR )
      maxSURFAPP = 0.0D+00

      do 11 IETA = 0, 90
         ETA = dfloat( IETA ) * RPAS * PI180
         do 12 IPHI = 0, 30
            PHI = dfloat( IPHI ) * RPAS * PI180
            do 13 IZETA = 0, 90
               ZETA = dfloat( IZETA ) * RPAS * PI180
               preSURFAPP = SURFAPP(RR,LL,ETA,ZETA,PHI,R32)
               if( preSURFAPP.gt.maxSURFAPP )then
!aaa              write(6,*) ETA/PI180,PHI/PI180,      &
!aaa                         ZETA/PI180,preSURFAPP,maxSURFAPP
                   maxETA  = ETA
                   maxPHI  = PHI
                   maxZETA = ZETA
                   maxSURFAPP = preSURFAPP
               endif
               if( preSURFAPP.lt.minSURFAPP )then
!aaa              write(6,*) ETA/PI180,PHI/PI180,      &
!aaa                         ZETA/PI180,preSURFAPP,minSURFAPP
                   minETA  = ETA
                   minPHI  = PHI
                   minZETA = ZETA
                   minSURFAPP = preSURFAPP
               endif
   13       enddo
   12    enddo
   11 enddo

      return
      end


! ############################################################################

! Routine TEST_IMPACT:

! Verification des regles de numerotation des faces, apres IMPACT.

      subroutine TEST_IMPACT(NPHDIF,NFACE,ETA,ZETA,PHI,PI180,EPS)

      IMPLICIT NONE

      INTEGER                                    :: NPHDIF                                         ! nombre total de photons diffuses [adim.];
      INTEGER                                    :: NFACE                                          ! Numero de la FACE concernee;

      REAL(KIND=8)                               :: EPS                                            ! 'precision' de certains calculs;
      REAL(KIND=8)                               :: ETA                                            ! angle entre le faisceau incident et l'axe OZ' [adim.];
      REAL(KIND=8)                               :: PHI                                            ! angle entre la projection OXp du faisceau incident et l'axe OX [adim.];
      REAL(KIND=8)                               :: PI180                                          ! = PI / 180;
      REAL(KIND=8)                               :: ZETA                                           ! angle entre les axes OZ et OZ' [adim.];

      if( NFACE.lt.0 )then
          write(6,'(''NPHDIF = '',i9.9,'' NFACE = '',i1)') &
                    NPHDIF,NFACE
          write(6,'(''routine IMPACT: face non identifiee'')')
          return
      endif

      if( ( dabs(  ETA/PI180-90.0D+00 ).lt.EPS ) .and. &
          ( dabs( ZETA/PI180-90.0D+00 ).lt.EPS ) .and. &
          ( dabs(  PHI/PI180- 0.0D+00 ).lt.EPS ) )then
          if( NFACE.eq.6 )then
              write(6,'(''NPHDIF = '',i9.9,'' NFACE = '',i1)') &
                        NPHDIF,NFACE
              write(6,'(''routine IMPACT: face interdite (a)'')')
              return
          endif
      endif

      if( ( dabs(  ETA/PI180-90.0D+00 ).lt.EPS ) .and. &
          ( dabs( ZETA/PI180-90.0D+00 ).lt.EPS ) .and. &
          ( dabs(  PHI/PI180-30.0D+00 ).lt.EPS ) )then
          if( NFACE.eq.5 .or. NFACE.eq.6 )then
              write(6,'(''NPHDIF = '',i9.9,'' NFACE = '',i1)') &
                        NPHDIF,NFACE
              write(6,'(''routine IMPACT: face interdite (b)'')')
              return
          endif
      endif

      if( ( (  ETA/PI180 ).lt.EPS ) .and. &
          ( ( ZETA/PI180 ).lt.EPS ) )then
          if( NFACE.ne.6 )then
              write(6,'(''NPHDIF = '',i9.9,'' NFACE = '',i1)') &
                        NPHDIF,NFACE
              write(6,'(''routine IMPACT: face interdite (c)'')')
              return
          endif
      endif

      return
      end



! ############################################################################

! Routine CHANG_AG:

! Cosinus directeurs des faisceaux reflechi (ALPHAr,BBETAr,GAMMAr)
!                   et transmis (refracte') (ALPHAt,BBETAt,GAMMAt), 
! a partir des cosinus directeurs des faces (ALPHAf,BBETAf,GAMMAf)
!          et des ceux du faisceau incident (ALPHAi,BBETAi,GAMMAi);
! et coefficients de reflexion (CRPAR,CRPER)
!           et de transmission (CTPAR,CTPER) de Fresnel.

! Cette version considere un faisceau incident sur une face externe de la
! particule, donc depuis le milieu '1' (air) vers le milieu '2' (glace).

! L'absorption est prise en compte selon l'approche de Chang et al. (2005);
 
      subroutine CHANG_AG(RN,RM,EPS,                       &
                          ALPHAf,BBETAf,GAMMAf,NFACE,      &
                          ALPHAi,BBETAi,GAMMAi,COSTHETAi,  &
                          ALPHAr,BBETAr,GAMMAr,            &
                          ALPHAt,BBETAt,GAMMAt,COSTHETAt,  &
                          GRANDN2,GRANDK2,COSANGLEHETERO,  &
                          ReCRPAR,ReCRPER,ReCTPAR,ReCTPER, &
                          ImCRPAR,ImCRPER,ImCTPAR,ImCTPER)

      IMPLICIT NONE

      REAL(KIND=8)                               :: ReNuF,ImNuF                                    ! coefficients de Fresnel: numerateur;
      REAL(KIND=8)                               :: ReDeF,ImDeF                                    ! coefficients de Fresnel: denominateur;
      REAL(KIND=8)                               :: ReCRPAR,ImCRPAR                                ! CRPAR: parties reelle et imaginaire;
      REAL(KIND=8)                               :: ReCRPER,ImCRPER                                ! CRPER: parties reelle et imaginaire;
      REAL(KIND=8)                               :: ReCTPAR,ImCTPAR                                ! CTPAR: parties reelle et imaginaire;
      REAL(KIND=8)                               :: ReCTPER,ImCTPER                                ! CTPER: parties reelle et imaginaire;

      REAL(KIND=8), DIMENSION(0:7)               :: ALPHAf                                         ! 1er cosinus directeur, faces de la colonne;
      REAL(KIND=8)                               :: ALPHAi                                         ! 1er cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: ALPHAr                                         ! 1er cosinus directeur, faisceau reflechi;
      REAL(KIND=8)                               :: ALPHAt                                         ! 1er cosinus directeur, faisceau transmis;
      REAL(KIND=8)                               :: BBBB                                           ! variable intermediaire (calcul de GRANDN2);

      REAL(KIND=8), DIMENSION(0:7)               :: BBETAf                                         ! 2eme cosinus directeur, faces;
      REAL(KIND=8)                               :: BBETAi                                         ! 2eme cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: BBETAr                                         ! 2eme cosinus directeur, faisceau reflechi;
      REAL(KIND=8)                               :: BBETAt                                         ! 2eme cosinus directeur, faisceau transmis;
      REAL(KIND=8)                               :: CCCC                                           ! variable intermediaire (calcul de GRANDN2);
      REAL(KIND=8)                               :: COSANGLEHETERO                                 ! cosinus de l'angle d'heterogeneite;
      REAL(KIND=8)                               :: COSTHETAi                                      ! cosinus de l'angle THETAi d'incidence;
      REAL(KIND=8)                               :: COSTHETAt                                      ! cosinus de l'angle THETAt de refraction; 
      REAL(KIND=8)                               :: EPS                                            ! 'precision' de plusieurs calculs;

      REAL(KIND=8), DIMENSION(0:7)               :: GAMMAf                                         ! 3eme cosinus directeur, faces;
      REAL(KIND=8)                               :: GAMMAi                                         ! 3eme cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: GAMMAr                                         ! 3eme cosinus directeur, faisceau reflechi;
      REAL(KIND=8)                               :: GAMMAt                                         ! 3eme cosinus directeur, faisceau transmis;
      REAL(KIND=8)                               :: GRANDN2                                        ! valeur apparente de 'petitn2';
      REAL(KIND=8)                               :: GRANDNs                                        ! variable intermediaire (calcul de GRANDN2);
      REAL(KIND=8)                               :: GRANDK2                                        ! valeur apparente de 'petitk2';
      REAL(KIND=8)                               :: petitk2                                        ! partie imaginaire de l'indice, milieu 2;
      REAL(KIND=8)                               :: petitn2                                        ! partie reelle de l'indice, milieu 2;
      REAL(KIND=8)                               :: RM                                             ! parties imaginaire et reelle de
      REAL(KIND=8)                               :: RN                                             ! l'indice de refraction de la glace;
      REAL(KIND=8)                               :: SINTHETAi                                      ! sinus de l'angle d'incidence;

      INTEGER                                    :: NFACE                                          ! indicateur de la face atteinte;

! Angle THETA d'incidence sur la face en question: ###########################

      COSTHETAi = ALPHAi * ALPHAf(NFACE) + &
                  BBETAi * BBETAf(NFACE) + &
                  GAMMAi * GAMMAf(NFACE) 
      IF (COSTHETAi .GT.  1.0D0) COSTHETAi =  1.0D0
      IF (COSTHETAi .LT. -1.0D0) COSTHETAi = -1.0D0
      SINTHETAi = dsqrt( 1.0D+00 - COSTHETAi * COSTHETAi )

! Cosinus directeurs du faisceau reflechi: ###################################

      ALPHAr = ALPHAi - 2.0D+00 * COSTHETAi * ALPHAf(NFACE)
      BBETAr = BBETAi - 2.0D+00 * COSTHETAi * BBETAf(NFACE)
      GAMMAr = GAMMAi - 2.0D+00 * COSTHETAi * GAMMAf(NFACE)

! Parties reelle et imaginaire 'apparentes' 
! de l'indice de refraction du deuxieme milieu: ##############################

      GRANDNs = SINTHETAi

      petitn2 = RN
      petitk2 = RM

      CCCC = GRANDNs * GRANDNs * ( petitn2 * petitn2 - petitk2 * petitk2 ) - &
             petitn2 * petitk2 * petitn2 * petitk2
      BBBB = - ( GRANDNs * GRANDNs + petitn2 * petitn2 - petitk2 * petitk2 )

      GRANDN2 = dsqrt( ( 1.0D+00 / 2.0D+00 ) *                              &
                       ( - BBBB + dsqrt( BBBB * BBBB - 4.0D+00 * CCCC ) ) )

      GRANDK2 = GRANDN2 * GRANDN2 - ( petitn2 * petitn2 - petitk2 * petitk2 )

      if( GRANDK2.lt.EPS )then
          GRANDK2 = petitk2
      else
          GRANDK2 = dsqrt( GRANDK2 )
      endif
      
! Angle THETA de refraction sur la face en question: #########################

      COSTHETAt = - dsqrt( 1.0D+00 -                           &
                           ( 1.0D+00 / GRANDN2 ) * SINTHETAi * &
                           ( 1.0D+00 / GRANDN2 ) * SINTHETAi ) 
      IF (COSTHETAt .GT.  1.0D0) COSTHETAt =  1.0D0
      IF (COSTHETAt .LT. -1.0D0) COSTHETAt = -1.0D0

! Cosinus directeurs du faisceau transmis (refracte'): #######################

      ALPHAt = ALPHAi / GRANDN2 +                                  & 
               ( COSTHETAt - COSTHETAi / GRANDN2 ) * ALPHAf(NFACE)
      BBETAt = BBETAi / GRANDN2 +                                  & 
               ( COSTHETAt - COSTHETAi / GRANDN2 ) * BBETAf(NFACE)
      GAMMAt = GAMMAi / GRANDN2 +                                  & 
               ( COSTHETAt - COSTHETAi / GRANDN2 ) * GAMMAf(NFACE)

! Angle HETERO d'heterogeneite associe' au faisceau transmis: ################

      COSANGLEHETERO = - COSTHETAt

! Coefficients de reflexion et de transmission de Fresnel: ###################

      ReNuF = (petitn2*petitn2-petitk2*petitk2)*COSTHETAi-GRANDN2*COSTHETAt
      ImNuF = 2.0D+00*petitn2*petitk2*COSTHETAi+GRANDK2
      ReDeF = (petitn2*petitn2-petitk2*petitk2)*COSTHETAi+GRANDN2*COSTHETAt
      ImDeF = 2.0D+00*petitn2*petitk2*COSTHETAi-GRANDK2

      ReCRPAR = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
      ImCRPAR = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

      ReNuF = 2.0D+00*petitn2*COSTHETAi
      ImNuF = 2.0D+00*petitk2*COSTHETAi

      ReCTPAR = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
      ImCTPAR = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

      ReNuF = COSTHETAi-GRANDN2*COSTHETAt
      ImNuF = GRANDK2
      ReDeF = COSTHETAi+GRANDN2*COSTHETAt
      ImDeF = -GRANDK2

      ReCRPER = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
      ImCRPER = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

      ReNuF = 2.0D+00*COSTHETAi
      ImNuF = 0.0D+00

      ReCTPER = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
      ImCTPER = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

      return
      end


! ############################################################################

! Routine CHANG_GA:

! Cosinus directeurs des faisceaux reflechi (ALPHAr,BBETAr,GAMMAr)
!                   et transmis (refracte') (ALPHAt,BBETAt,GAMMAt), 
! a partir des cosinus directeurs des faces (ALPHAf,BBETAf,GAMMAf)
!          et des ceux du faisceau incident (ALPHAi,BBETAi,GAMMAi);
! et coefficients de reflexion (CRPAR,CRPER)
!           et de transmission (CTPAR,CTPER) de Fresnel.

! Cette version considere un faisceau incident sur une face interne de la
! particule, donc depuis le milieu '1' (glace) vers le milieu '2' (air).

! Les parties reelles et imaginaires 'petitn2' et 'petitk2' correspondent a
! l'air tandis que 'petitn1' et 'petitk1' correspondent a la glace pure;

! L'absorption est prise en compte selon l'approche de Chang et al. (2005);
 
      subroutine CHANG_GA(petitn1,petitk1,                      &
                          ALPHAf,BBETAf,GAMMAf,NFACE,           &
                          ALPHAi,BBETAi,GAMMAi,COSTHETAi,       &
                          ALPHAr,BBETAr,GAMMAr,                 &
                          ALPHAt,BBETAt,GAMMAt,COSTHETAt,IREFL, &
                          GRANDN1,GRANDK1,COSANGLEHETERO,       &
                          ReCRPAR,ReCRPER,ReCTPAR,ReCTPER,      &
                          ImCRPAR,ImCRPER,ImCTPAR,ImCTPER)

      IMPLICIT NONE

      REAL(KIND=8)                               :: ReNuF,ImNuF            ! coefficients de Fresnel: numerateur;
      REAL(KIND=8)                               :: ReDeF,ImDeF            ! coefficients de Fresnel: denominateur;
      REAL(KIND=8)                               :: ReCRPAR,ImCRPAR        ! CRPAR: parties reelle et imaginaire;
      REAL(KIND=8)                               :: ReCRPER,ImCRPER        ! CRPER: parties reelle et imaginaire;
      REAL(KIND=8)                               :: ReCTPAR,ImCTPAR        ! CTPAR: parties reelle et imaginaire;
      REAL(KIND=8)                               :: ReCTPER,ImCTPER        ! CTPER: parties reelle et imaginaire;

      REAL(KIND=8), DIMENSION(0:7)               :: ALPHAf                 ! 1er cosinus directeur, faces de la colonne;
      REAL(KIND=8)                               :: ALPHAi                 ! 1er cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: ALPHAr                 ! 1er cosinus directeur, faisceau reflechi;
      REAL(KIND=8)                               :: ALPHAt                 ! 1er cosinus directeur, faisceau transmis;

      REAL(KIND=8), DIMENSION(0:7)               :: BBETAf                 ! 2eme cosinus directeur (faces);
      REAL(KIND=8)                               :: BBETAi                 ! 2eme cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: BBETAr                 ! 2eme cosinus directeur, faisceau reflechi;
      REAL(KIND=8)                               :: BBETAt                 ! 2eme cosinus directeur, faisceau transmis;

      REAL(KIND=8), DIMENSION(0:7)               :: GAMMAf                 ! 3eme cosinus directeur (faces);
      REAL(KIND=8)                               :: GAMMAi                 ! 3eme cosinus directeur, faisceau incident;
      REAL(KIND=8)                               :: GAMMAr                 ! 3eme cosinus directeur, faisceau reflechi;
      REAL(KIND=8)                               :: GAMMAt                 ! 3eme cosinus directeur, faisceau transmis;

      REAL(KIND=8)                               :: COSANGLEHETERO         ! cosinus de l'angle d'heterogeneite;
      REAL(KIND=8)                               :: COSPSIi                ! cosinus de l'angle PSIi d'incidence;
      REAL(KIND=8)                               :: COSTHETAi              ! cosinus de l'angle THETAi d'incidence;
      REAL(KIND=8)                               :: COSTHETAt              ! cosinus de l'angle THETAt de refraction; 
      REAL(KIND=8)                               :: delta             

      REAL(KIND=8)                               :: GRANDK1                ! valeur apparente de 'petitk1';
      REAL(KIND=8)                               :: GRANDN1                ! valeur apparente de 'petitn1';
      REAL(KIND=8)                               :: GRANDNi                ! var. aux. (reflexion totale);
      REAL(KIND=8)                               :: petitk1                ! partie imaginaire de l'indice, milieu 1;
      REAL(KIND=8)                               :: petitn1                ! partie reelle de l'indice, milieu 1;
      REAL(KIND=8)                               :: SINTHETAi              ! sinus de l'angle d'incidence;

      INTEGER                                    :: IREFL                  ! Indicateur de REFLexion totale ou non;
      INTEGER                                    :: NFACE                  ! indicateur de la face atteinte;


! Le cas particulier de reflexion totale est a priori absent: ################

      IREFL = 0

! Angle THETA d'incidence sur la face en question: ###########################

      COSTHETAi = ALPHAi * ALPHAf(NFACE) + &
                  BBETAi * BBETAf(NFACE) + &
                  GAMMAi * GAMMAf(NFACE) 
      IF (COSTHETAi .GT.  1.0D0) COSTHETAi =  1.0D0
      IF (COSTHETAi .LT. -1.0D0) COSTHETAi = -1.0D0

      SINTHETAi = dsqrt( 1.0D+00 - COSTHETAi * COSTHETAi )

      COSPSIi = dcos( dacos( COSANGLEHETERO ) - &
                      dacos( COSTHETAi ) )

! Cosinus directeurs du faisceau reflechi: ###################################

      ALPHAr = ALPHAi - 2.0D+00 * COSTHETAi * ALPHAf(NFACE)
      BBETAr = BBETAi - 2.0D+00 * COSTHETAi * BBETAf(NFACE)
      GAMMAr = GAMMAi - 2.0D+00 * COSTHETAi * GAMMAf(NFACE)

! Le cas particulier de reflexion totale peut avoir lieu: ####################

      if( SINTHETAi.gt.( 1.0D+00/GRANDN1 ) )then
          IREFL = 1
          GRANDNi = dsqrt( GRANDN1*GRANDN1*                &
                           SINTHETAi*SINTHETAi - 1.0D+00 )

          ReNuF = GRANDN1*COSTHETAi+2.0D+00*petitn1*petitk1*GRANDNi
          ImNuF = GRANDK1*COSPSIi-(petitn1*petitn1-petitk1*petitk1)*GRANDNi
          ReDeF = GRANDN1*COSTHETAi-2.0D+00*petitn1*petitk1*GRANDNi
          ImDeF = GRANDK1*COSPSIi+(petitn1*petitn1-petitk1*petitk1)*GRANDNi

          ReCRPAR = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
          ImCRPAR = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

          ReCTPAR = 0.0D+00
          ImCTPAR = 0.0D+00

          ReNuF = GRANDN1*COSTHETAi
          ImNuF = GRANDK1*COSPSIi-GRANDNi
          ReDeF = GRANDN1*COSTHETAi
          ImDeF = GRANDK1*COSPSIi+GRANDNi

          ReCRPER = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
          ImCRPER = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

          ReCTPER = 0.0D+00
          ImCTPER = 0.0D+00

! Les coefficients de reflexion de Fresnel EN INTENSITE, CRPAR*dconjg(CRPAR) #
! et CRPER*dconjg(CRPER), doivent etre unitaires: ############################

          delta = 1.0D+00 / sqrt( ReCRPAR*ReCRPAR + ImCRPAR*ImCRPAR )
          ReCRPAR = delta*ReCRPAR
          ImCRPAR = delta*ImCRPAR
          delta = 1.0D+00 / sqrt( ReCRPER*ReCRPER + ImCRPER*ImCRPER )
          ReCRPER = delta*ReCRPER
          ImCRPER = delta*ImCRPER
      else

! Cosinus directeurs du faisceau transmis (refracte'): #######################

          COSTHETAt = dsqrt( 1.0D+00 -               &
                             GRANDN1 * GRANDN1 *     &
                             SINTHETAi * SINTHETAi )
          IF (COSTHETAt .GT.  1.0D0) COSTHETAt =  1.0D0
          IF (COSTHETAt .LT. -1.0D0) COSTHETAt = -1.0D0
          ALPHAt = GRANDN1 * ALPHAi +                                  &
                   ( COSTHETAt - GRANDN1 * COSTHETAi ) * ALPHAf(NFACE)
          BBETAt = GRANDN1 * BBETAi +                                  &
                   ( COSTHETAt - GRANDN1 * COSTHETAi ) * BBETAf(NFACE)
          GAMMAt = GRANDN1 * GAMMAi +                                  &
                   ( COSTHETAt - GRANDN1 * COSTHETAi ) * GAMMAf(NFACE)


! Coefficients de reflexion et de transmission de Fresnel: ###################

          ReNuF = GRANDN1*COSTHETAi-(petitn1*petitn1-petitk1*petitk1)*COSTHETAt
          ImNuF = GRANDK1*COSPSIi-2.0D+00*petitn1*petitk1*COSTHETAt
          ReDeF = GRANDN1*COSTHETAi+(petitn1*petitn1-petitk1*petitk1)*COSTHETAt
          ImDeF = GRANDK1*COSPSIi+2.0D+00*petitn1*petitk1*COSTHETAt

          ReCRPAR = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
          ImCRPAR = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

          ReNuF = 2.0D+00*(petitn1*GRANDN1*COSTHETAi-petitk1*GRANDK1*COSPSIi) 
          ImNuF = 2.0D+00*(petitn1*GRANDK1*COSPSIi+petitk1*GRANDN1*COSTHETAi)

          ReCTPAR = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
          ImCTPAR = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

          ReNuF = GRANDN1*COSTHETAi-COSTHETAt
          ImNuF = GRANDK1*COSPSIi
          ReDeF = GRANDN1*COSTHETAi+COSTHETAt
          ImDeF = GRANDK1*COSPSIi

          ReCRPER = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
          ImCRPER = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

          ReNuF = 2.0D+00*GRANDN1*COSTHETAi
          ImNuF = 2.0D+00*GRANDK1*COSPSIi

          ReCTPER = ( ReNuF*ReDeF+ImNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )
          ImCTPER = ( ImNuF*ReDeF-ReNuF*ImDeF )/( ReDeF*ReDeF+ImDeF*ImDeF )

      endif

      return
      end


! ############################################################################

! Routine STOKES:
!
! Calcul des contributions a la matrice de diffusion, celle qui transforme les
! parametres de Stokes du faisceau incident dans ceux du faisceau diffuse'.
!
! Les calculs considerent une formulation encore plus explicite que celle pris
! en compte par Mishchenko et al. (2002), expressions (2.106-2.121).

      subroutine STOKES(ReAUX22GTGP,ImAUX22GTGP,aux44GTGP)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(2,2)               :: ReAUX22GTGP                                    ! variable auxiliaire (entree de STOKES);
      REAL(KIND=8), DIMENSION(2,2)               :: ImAUX22GTGP                                    ! variable auxiliaire (entree de STOKES);
      REAL(KIND=8), DIMENSION(4,4)               :: AUX44GTGP                                      ! variable auxiliaire (sortie de STOKES);
      REAL(KIND=8)                               :: RS11,RS12                                      ! partie reelle des amplitudes;
      REAL(KIND=8)                               :: RS21,RS22                                      ! "
      REAL(KIND=8)                               :: IS11,IS12                                      ! partie imaginaire des amplitudes;
      REAL(KIND=8)                               :: IS21,IS22                                      ! "


      RS22 = ReAUX22GTGP(2,2)
      RS11 = ReAUX22GTGP(1,1)
      RS12 = ReAUX22GTGP(1,2)
      RS21 = ReAUX22GTGP(2,1)

      IS22 = ImAUX22GTGP(2,2)
      IS11 = ImAUX22GTGP(1,1)
      IS12 = ImAUX22GTGP(1,2)
      IS21 = ImAUX22GTGP(2,1)

      aux44GTGP(1,1) = ( RS11*RS11 + IS11*IS11 +           &
                         RS12*RS12 + IS12*IS12 +           &
                         RS21*RS21 + IS21*IS21 +           &
                         RS22*RS22 + IS22*IS22 ) / 2.0D+00

      aux44GTGP(1,2) = ( RS11*RS11 + IS11*IS11 -           &
                         RS12*RS12 - IS12*IS12 +           &
                         RS21*RS21 + IS21*IS21 -           &
                         RS22*RS22 - IS22*IS22 ) / 2.0D+00

      aux44GTGP(1,3) = - ( RS11*RS12 + RS21*RS22 + IS11*IS12 + IS21*IS22 )

      aux44GTGP(1,4) = - ( RS12*IS11 - RS21*IS22 - RS11*IS12 + RS22*IS21 )

      aux44GTGP(2,1) = ( RS11*RS11 + IS11*IS11 +           &
                         RS12*RS12 + IS12*IS12 -           &
                         RS21*RS21 - IS21*IS21 -           &
                         RS22*RS22 - IS22*IS22 ) / 2.0D+00

      aux44GTGP(2,2) = ( RS11*RS11 + IS11*IS11 -           &
                         RS12*RS12 - IS12*IS12 -           &
                         RS21*RS21 - IS21*IS21 +           &
                         RS22*RS22 + IS22*IS22 ) / 2.0D+00

      aux44GTGP(2,3) = - ( RS11*RS12 - RS21*RS22 + IS11*IS12 - IS21*IS22 )

      aux44GTGP(2,4) = - ( RS12*IS11 + RS21*IS22 - RS11*IS12 - RS22*IS21 )

      aux44GTGP(3,1) = - ( RS11*RS21 + RS12*RS22 + IS11*IS21 + IS12*IS22 )

      aux44GTGP(3,2) = - ( RS11*RS21 - RS12*RS22 + IS11*IS21 - IS12*IS22 )

      aux44GTGP(3,3) = RS11*RS22 + RS12*RS21 + IS11*IS22 + IS12*IS21

      aux44GTGP(3,4) = RS22*IS11 + RS12*IS21 - RS11*IS22 - RS21*IS12

      aux44GTGP(4,1) = - ( RS11*IS21 + RS12*IS22 - RS21*IS11 - RS22*IS12 )

      aux44GTGP(4,2) = - ( RS11*IS21 - RS12*IS22 - RS21*IS11 + RS22*IS12 )

      aux44GTGP(4,3) = RS11*IS22 - RS21*IS12 - RS22*IS11 + RS12*IS21

      aux44GTGP(4,4) = RS11*RS22 - RS12*RS21 + IS11*IS22 - IS12*IS21

      return
      end

! ############################################################################

! Routine DIFFRACT:

! Calcul de la distribution angulaire de l'energie diffractee par une colonne
! hexagonale, en suivant l'approche presentee par Brogniez (1992, pp.101-106).

! La notation a ete fortement modifiee, en se raprochant de celle adoptee par
! G. Brogniez; p.ex., les cosinus directeurs du faisceau diffracte (ALPH,BETA)
! ont ete renommes MMU et ETA (voir page 102).

subroutine DIFFRACT(NU,KKKK,     & ! les photons..........
     LL,RR,SE,auxSURFAPP,        & ! la particule.........
     angleETA,ZETA,PHI,          & ! orientation..........
     NDGT,NDGTFR,RESOLGT,COSDGT, & ! grille GRANDTHETA....
     EPS,PI,PI180,PI3,           & ! constantes...........
     auxNPNPT,                   & ! "precision"..........
     auxEDFR)                      ! resultats............

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(2,3)               :: A                                              ! 
  REAL(KIND=8)                               :: LL                                             ! longueur de la colonne hexagonale [1E-6 m];
  REAL(KIND=8)                               :: ETA                                            ! 2eme cosinus directeur, faisc. diffracte';
  REAL(KIND=8)                               :: MMU                                            ! 1er cosinus directeur, faisc. diffracte';
  REAL(KIND=8)                               :: auxSURFAPP                                     ! surface apparente prise en compte; 

  REAL(KIND=8), DIMENSION(NDGT,2)            :: COSDGT                                         ! cosinus de DGT(IDGT,1 et 2);

  REAL(KIND=8)                               :: COSTETD                                        ! 
  REAL(KIND=8)                               :: COSTHETA                                       ! cosinus de l'angle THETA;

  REAL(KIND=8), DIMENSION(9)                 :: D                                              !  

  REAL(KIND=8)                               :: ED                                             !            
  REAL(KIND=8)                               :: EDT                                            !

  REAL(KIND=8), DIMENSION(NDGTFR)            :: auxEDFR                                        ! variable auxiliaire au calcul d'EDFR;

  REAL(KIND=8)                               :: sumEDFR                                        ! energie totale diffractee;
  REAL(KIND=8)                               :: EPS                                            ! 'precision' de plusieurs calculs;
  REAL(KIND=8)                               :: angleETA                                       ! angle entre le faisceau incident et l'axe OZ' [degres];
  REAL(KIND=8)                               :: FINCFR                                         ! (densite de) flux d'energie incidente;
  REAL(KIND=8)                               :: NU                                             ! nombre d'onde (=10000/ALAMB) [cm-1];
  REAL(KIND=8)                               :: PHI                                            ! angle entre la projection OXp du faisceau incident et l'axe OX [adim.];
  REAL(KIND=8)                               :: PHID                                           ! 
  REAL(KIND=8)                               :: PHIDI                                          ! 
  REAL(KIND=8)                               :: PHIDS                                          ! 
  REAL(KIND=8)                               :: PI                                             ! = constante PI;
  REAL(KIND=8)                               :: PI180                                          ! = PI / 180;
  REAL(KIND=8)                               :: PI3                                            !

  REAL(KIND=8), DIMENSION(2,8)               :: R                                              ! variable intermediaire (entre SE et X,Y);

  REAL(KIND=8)                               :: RESOLGT                                        ! resolution angulaire de la matrice de diffusion, en GRANDTHETA [degres];
  REAL(KIND=8)                               :: RR                                             ! rayon de la colonne hexagonale [1E-6 m];
  REAL(KIND=8)                               :: KKKK                                           ! = 2*PI/(LongueurD'Onde);

  REAL(KIND=8), DIMENSION(3,10)              :: SE                                             ! coordonnees XYZ des dix sommets eclaires;

  REAL(KIND=8)                               :: SINTETD                                        ! 
  REAL(KIND=8)                               :: SINTHETA                                       ! sinus de l'angle THETA;
  REAL(KIND=8)                               :: TETDF                                          ! 
  REAL(KIND=8)                               :: TETDI                                          ! 

  REAL(KIND=8), DIMENSION(5)                 :: X                                              ! coordonnee 'Xi' des cinq sommets;
  REAL(KIND=8), DIMENSION(5)                 :: Y                                              ! coordonnee 'Yi' des cinq sommets;

  REAL(KIND=8)                               :: ZETA                                           ! angle entre les axes OZ et OZ' [adim.];

  REAL(KIND=8), DIMENSION(0:2048)            :: AINTN                                          !
  REAL(KIND=8), DIMENSION(0:2048)            :: EDD                                            !
  REAL(KIND=8), DIMENSION(0:2048)            :: PHID1                                          ! 
  REAL(KIND=8), DIMENSION(0:2048)            :: PHID2                                          !

  REAL(KIND=8)                               :: ReCAMP                                         ! partie reelle de 'CAMP';
  REAL(KIND=8)                               :: ReCINCI                                        ! partie reelle de 'CINCI';

  INTEGER                                    :: auxNPNPT                                       ! variable auxiliaire ("precision" NP=NPT);
  INTEGER                                    :: I                                              ! 
  INTEGER                                    :: IDGT                                           ! indicateur du Domaine GrandTheta (NDGTFR);
  INTEGER                                    :: ILB                                            ! 
  INTEGER                                    :: ILT                                            ! 
  INTEGER                                    :: NDGT                                           ! Nombre total des Domaines GrandTheta;
  INTEGER                                    :: NDGTFR                                         ! NDGT 'difFRaction' (resolution RESOLGT);

  INTEGER                                    :: NP,NPT,L,LT                                    ! grilles en GRANDPHI et GRANDTHETA;


  ! Grille en GRANDPHI entre 0 et PI: ##########################################

  !write(*,*) 'start-DIFFRACT' 

  NP = auxNPNPT
  L  = 2 * ( 2**NP )

  ! Grille en GRANDTHETA entre COSDGT(IDGT,1) et COSDGT(IDGT,2): ###############

  NPT = auxNPNPT
  LT  = 2 * ( 2**NPT )

  ! Les sommets eclaires externes de la particule (numerotes de 1 a 8) sont ####
  ! exprimes dans le systeme OX'Y'Z' attache au faisceau incident: #############

  SINTHETA = dcos( ZETA ) * dcos( angleETA )
  COSTHETA = dsqrt( 1.0D+00 - SINTHETA * SINTHETA )

  A(1,1) = - dsin( PHI ) 
  A(2,1) = - dcos( PHI ) * SINTHETA
  A(1,2) =   dcos( PHI )
  A(2,2) = - dsin( PHI ) * SINTHETA
  A(1,3) =  0.0D+00
  A(2,3) =  COSTHETA

  call DMRMUL(A,SE,R,2,3,8)

  do I = 1, 5
     X(I) = R(1,I)
     Y(I) = R(2,I)
  enddo

  ! Flux d'energie incidente sur la particule projetee, par unite de surface: ##

  FINCFR = 1.0D+00 / auxSURFAPP

  FINCFR = FINCFR / ( ( 1.0D+04 / NU ) * ( 1.0D+04 / NU ) )

  ! Angle de diffusion GRANDPHI: ###############################################

  PHIDI = 0.0D+00
  PHIDS = PI
  do ILB = 0, L
     PHID       = PI * dfloat( ILB ) / dfloat( L )
     PHID1(ILB) = dcos( PHID )
     PHID2(ILB) = dsin( PHID )
  enddo

  ! Debut de la premiere boucle: 
  ! tous les domaines DGT de l'angle GRANDTHETA de diffusion ###################

  do IDGT = 1, NDGTFR

     TETDI = COSDGT(IDGT,1)
     TETDF = COSDGT(IDGT,2) 

     ! Debut de la deuxieme boucle: 
     ! angles GRANDTHETA de diffusion entre DGT(IDGT,1) et DGT(IDGT,2): ###########

     do ILT = 0, LT
        COSTETD = TETDI +                          &
             ( TETDF - TETDI ) *              & 
             ( dfloat( ILT ) / dfloat( LT ) )
        SINTETD = dsqrt( 1.0D+00 - COSTETD * COSTETD )

        ! Debut de la troisieme boucle: 
        ! angles GRANDPHI de diffusion entre 0 et PI #################################

        do ILB = 0, L

           ! Cosinus directeurs du faisceau diffracte': #################################

           MMU = SINTETD * PHID1(ILB)
           ETA = SINTETD * PHID2(ILB)

           do I = 1, 5
              D(I) = MMU * X(I) + ETA * Y(I)
           enddo

           ReCAMP = 0.0D+00

           do I = 1, 4
              if( dabs( X(I)-X(I+1) ).lt.EPS .and. &
                   dabs( Y(I)-Y(I+1) ).lt.EPS )then
                 exit
              endif

              ! Integration sur chaque triangle, D(I)=0 et D(I+1) non nulle: ###############

              if( dabs( D(I+1) ).gt.EPS .and. &
                   dabs( D(I) )  .lt.EPS )then

                 ReCINCI = ( ( X(I+1)*Y(I)-X(I)*Y(I+1) )/      &
                      ( KKKK*KKKK*D(I+1)*D(I+1) ) )*    &
                      ( - 1.0D+00 + dcos( KKKK*D(I+1) ) )
              endif

              ! Integration sur chaque triangle, D(I+1)=0 et D(I) non nulle: ###############

              if( dabs( D(I+1) ).lt.EPS .and. &
                   dabs( D(I) )  .gt.EPS )then

                 ReCINCI = ( ( X(I+1)*Y(I)-X(I)*Y(I+1) )/  &
                      ( -KKKK*KKKK*D(I)*D(I) ) )*   &
                      ( 1.0D+00 - dcos( KKKK*D(I) ) )
              endif

              ! Integration sur chaque triangle, D(I)=D(I+1), non nulles : #################

              if( dabs( D(I+1)-D(I) ).lt.EPS .and. &
                   dabs( D(I+1) )     .gt.EPS .and. &
                   dabs( D(I) )       .gt.EPS )then

                 ReCINCI = ( ( X(I+1)*Y(I)-X(I)*Y(I+1) )/      &
                      ( KKKK*KKKK*D(I)*D(I) ) )*        &
                      ( 1.0D+00 - dcos( KKKK*D(I) )     &
                      - KKKK*D(I)*dsin( KKKK*D(I) ) )
              endif

              ! Integration sur chaque triangle, D(I)=D(I+1), nulles : #####################

              if( dabs( D(I+1)-D(I) ).lt.EPS .and. &
                   dabs( D(I+1) )     .lt.EPS .and. &
                   dabs( D(I) )       .lt.EPS )then

                 ReCINCI = ( X(I)*Y(I+1)-X(I+1)*Y(I) )/2.0D+00

              endif

              ! Integration sur chaque triangle, cas general: ##############################

              if( dabs( D(I+1)-D(I) ).gt.EPS .and. &
                   dabs( D(I+1) )     .gt.EPS .and. &
                   dabs( D(I) )       .gt.EPS )then

                 ReCINCI = ( ( X(I+1)*Y(I)-X(I)*Y(I+1) )/                 &
                      ( KKKK*KKKK*D(I)*D(I+1)*( D(I+1)-D(I) ) ) )* &
                      ( D(I+1) - D(I) +                              &
                      D(I)*dcos( KKKK*D(I+1) ) -                   &
                      D(I+1)*dcos( KKKK*D(I) ) )
              endif

              ReCAMP = ReCAMP + ReCINCI

           enddo

           ! La prise en compte des quatre derniers triangles est remplacee par la double
           ! prise en compte des quatre derniers triangles: #############################

           ReCAMP = 2.0D+00 * ReCAMP

           ! Intensite diffractee dans la direction (MMU,ETA): ##########################

           AINTN(ILB) = ReCAMP * ReCAMP

           ! Fin de la troisieme boucle: ################################################

        enddo

        ! Integration sur GRANDPHI entre 0 et PI: ####################################

        call HM_INTANG(PHIDI,PHIDS,AINTN,L,NP,ED)

        EDD(ILT) = 2.0D+00 * ED * COSTETD

        ! Fin de la deuxieme boucle: #################################################

     enddo

     ! Integration sur le cosinus de l'angle GRANDTHETA de diffusion, pour chaque #
     ! domaine delimite' par COSDGT(IDGT,1) et COSDGT(IDGT,2): ####################

     call HM_INTANG(TETDF,TETDI,EDD,LT,NPT,EDT)

     auxEDFR(IDGT) = FINCFR * EDT

     ! Fin de la premiere boucle: #################################################

  enddo

  ! Energie totale diffractee par la particule: ################################
  sumEDFR = 0.0D+00
  do IDGT = 1, NDGTFR
     sumEDFR = sumEDFR + auxEDFR(IDGT)
     !write(*,*) auxEDFR(IDGT)
     !write(*,*) sumEDFR
  enddo
 
  !write(6,*)
  !write(6,'(''Routine_DIFFRACT: '')')
  !write(6,'(''  surfapp_= '',1e12.6,'' 1E-12_m2'')') auxSURFAPP
  !write(6,'(''  sumEDFR_= '',f6.2,'' %'')') 100.0D+00*sumEDFR
  !     write(6,'(''EDFR_(sans_normalisation):'')')
  !     do IDGT = 1, NDGTFR
  !        write(6,*) IDGT,dacos( COSDGT(IDGT,1) )/PI180,         &
  !                   dacos( COSDGT(IDGT,2) )/PI180,auxEDFR(IDGT)
  !aaa enddo

  return
end subroutine DIFFRACT
