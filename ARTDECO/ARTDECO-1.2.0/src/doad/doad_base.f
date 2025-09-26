c-----------------------------------------------------------------------
c      
c  Author:      Michele Vesperini
c  Version:     96 04 12
c
c Mathieu Compiegne (fall 2012) :
c      - rename all routine with a prefix doad_      
c      - rename NDcfi  to NDcoef
c               NDmui  to NDmu
c               NDsupi to NDsup
c               NDirfi to NDirf
c      - no more use of include dimension_sub.h since
c        all that dimensions (NDcoef, NDmu, NDsup, NDirf) are passed as argument
c
c
c Adding and doubling RT code including internal source and polarization
c-----------------------------------------------------------------------


      subroutine doad_addsm(w,x,y,nsup,NDsup)
************************************************************************
*** add two NDsup*NDsup supermatrices w=x+y                          ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension w(NDsup,NDsup), x(NDsup,NDsup), y(NDsup,NDsup)

      do 10 j=1,nsup
      do 10 i=1,nsup
          w(i,j) = x(i,j)+y(i,j)
   10 continue
      return
      end


      subroutine doad_addsv(w,x,y,nsup,NDsup)
************************************************************************
*** add two NDsup supervectors w=x+y                                 ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension w(NDsup), x(NDsup), y(NDsup)

      do 10 i=1,nsup
          w(i) = x(i)+y(i)
   10 continue
      return
      end


      subroutine doad_addtot(nmat,nmutot,nmug,NDmu,NDsup,xmu,ebmu,eblow,
     +    r,t,radd,tadd)
************************************************************************
*** adding-algorithm for 2 inhomogeneous layers                      ***
*** on entrance r and t contain the reflection and transmission of   ***
*** the upper layer with ebmu, and radd and tadd contain the         ***
*** reflection and transmission of the lower layer(s) with eblow     ***
*** on exit radd and tadd contain the reflection and transmission of ***
*** the added layer and r and t contain u and d at the interface     ***
*** see de Haan et al. (1987) Eqs. (85)-(91)                         ***
************************************************************************
      implicit double precision (a-h,o-z)

      dimension xmu(NDmu), ebmu(NDmu), eblow(NDmu)
      dimension r(NDsup,NDsup), t(NDsup,NDsup)
      dimension radd(NDsup,NDsup), tadd(NDsup,NDsup)
      dimension s(NDsup,NDsup), e(NDsup,NDsup), p(NDsup,NDsup)
cmv   common /mat/ s,e,p

      ng = nmug*nmat
      nsup = nmat*nmutot
      do 10 j=1,nsup
      do 10 i=1,nsup
          p(i,j) = r(i,j)
   10 continue
      call doad_transhm(p,nmutot,nmat,NDsup)
************************************************************************
* calculate Q_1, Eq. (85)                                              *
* and e and s are initialized to Q_1                                   *
************************************************************************
      call doad_prodsm(e,p,radd,nmutot,nmug,nmat,NDsup)
      do 20 j=1,nsup
      do 20 i=1,nsup
          s(i,j) = e(i,j)
   20 continue
      delta = 0.d0
      do 30 i=1,ng
          delta = delta+e(i,i)
   30 continue
      if ((delta*delta) .lt. 1.d-12) goto 100
      ik = 0
************************************************************************
* productmethod to calculate S, Eqs. (112)-(115)                       *
* e = C_r*C_r                                                          *
************************************************************************
   40 call doad_prodsm(p,e,e,nmutot,nmug,nmat,NDsup)
      do 50 j=1,nsup
      do 50 i=1,nsup
          e(i,j) = p(i,j)
   50 continue 
      call doad_prodsm(p,s,e,nmutot,nmug,nmat,NDsup)
************************************************************************
* P = S_r*C_r+1                                                         *
************************************************************************
      do 60 j=1,nsup
      do 60 i=1,nsup
          s(i,j) = p(i,j)+e(i,j)+s(i,j)
   60 continue
************************************************************************
* S = S_r+1                                                            *
************************************************************************
      delta = 0.d0
      do 70 i=1,ng
          delta = delta+e(i,i)
   70 continue
************************************************************************
* delta=trace(C_r+1), Eq. (124), is an estimate of the error           *
************************************************************************
      ik = ik+1
      if (ik .gt. 20) then
          write(*,90)
          write(*,80)
          write(*,90)
   80     format('*** stop : S-series converge to slow ***')
   90     format('****************************************')
          stop 'slow convergence in ADDTOT'
      endif
      if ((delta*delta) .gt. 1.d-12) goto 40
************************************************************************
* calculate D, Eq. (88)                                                *
* the repeated reflections matrix is contained in s                    *
************************************************************************
  100 call doad_prodsm(p,s,t,nmutot,nmug,nmat,NDsup)
      do 110 i=1,nmutot
          im = nmat*(i-1)
          do 110 j=1,nmutot
              jm = nmat*(j-1)
              ej = ebmu(j)
              do 110 k=1,nmat
                  ik = im+k
                  do 110 l=1,nmat
                      jl = jm+l
                      p(ik,jl) = p(ik,jl)+ej*s(ik,jl)+t(ik,jl)
  110 continue
************************************************************************
* calculate U, Eq. (89)                                                *
************************************************************************
      call doad_prodsm(s,radd,p,nmutot,nmug,nmat,NDsup)
      do 120 i=1,nmutot
          im = nmat*(i-1)
          do 120 j=1,nmutot
              jm = nmat*(j-1)
              ej = ebmu(j)
              do 120 k=1,nmat
                  ik = im+k
                  do 120 l=1,nmat
                      jl = jm+l
                      s(ik,jl) = s(ik,jl)+ej*radd(ik,jl)
  120 continue
************************************************************************
* calculate R, Eq. (90)                                                *
************************************************************************
      call doad_transhm(t,nmutot,nmat,NDsup)
      call doad_prodsm(e,t,s,nmutot,nmug,nmat,NDsup)
      do 130 i=1,nmutot
          im = nmat*(i-1)
          ei = ebmu(i)
          do 130 j=1,nmutot
              jm = nmat*(j-1)
              do 130 k=1,nmat
                  ik = im+k
                  do 130 l=1,nmat
                      jl = jm+l
                      radd(ik,jl) = e(ik,jl)+ei*s(ik,jl)+r(ik,jl)
  130 continue
************************************************************************
* calculate T, Eq. (91)                                                *
************************************************************************
      call doad_prodsm(e,tadd,p,nmutot,nmug,nmat,NDsup)
      do 140 i=1,nmutot
          im = nmat*(i-1)
          ei = eblow(i)
          do 140 j=1,nmutot
              jm = nmat*(j-1)
              ej = ebmu(j)
              do 140 k=1,nmat
                  ik = im+k
                  do 140 l=1,nmat
                      jl = jm+l
                      tadd(ik,jl) = ej*tadd(ik,jl)+ei*p(ik,jl)+e(ik,jl)
  140 continue
      do 150 j=1,nsup
      do 150 i=1,nsup
          t(i,j) = p(i,j)
          r(i,j) = s(i,j)
  150 continue
      return
      end


      subroutine doad_bstart(m,il,NDlay,alfbet,NDcoef,ncoef,M0,b,
     +    eps,xmumin,
     +    ndoubl,nlirf)
************************************************************************
*** Calculate the optical thickness bs at which doubling should be   ***
*** started to obtain an an error less than eps in each individual   ***
*** homogeneous il with thickness b(il)                              ***
*** This is done for the m-th Fourier component                      ***
*** The value of bs is given by: bs = b*2**(-ndoubl)                 ***
*** see de Haan et al. (1987) p. 387                                 ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension alfbet(4,4,0:NDcoef,NDlay), b(0:NDlay), M0(NDlay)
      dimension ncoef(NDlay), ndoubl(NDlay), nlirf(0:NDlay+1)

************************************************************************
*  Calculate the effective albedo am, Eq. (143)                        *
************************************************************************
      am = 0.d0
      do 10 k=m,M0(il)
          tmp = dabs(alfbet(1,1,k,il))/dble(2*k+1)
          if (am .lt. tmp) am = tmp
   10 continue
************************************************************************
*  Right hand side of Eq. (142)                                        *
************************************************************************
      rhs = 4.d0*eps/(9.d0*b(il)*am**3*dble(2*m+1))
      ndoubl(il) = -1
      b0 = 2.d0*b(il)
   20     ndoubl(il) = ndoubl(il)+1
          b0 = 0.5d0*b0
          fb0mu = b0/xmumin
          if (fb0mu .gt. 1.d0) fb0mu = 1.d0
          if ((b0 .ge. (rhs/fb0mu)) .or. (b0 .ge. (1.d0/3.d0))) goto 20
************************************************************************
* avoid a starting thickness greater then the smallest irf level       *
************************************************************************
      if (nlirf(il) .gt. (2*ndoubl(il))) ndoubl(il) = (nlirf(il)+1)/2
      return
      end


      subroutine doad_double(nmat,nmutot,nmug,NDmu,NDsup,ebmu,r,t,u,d)
************************************************************************
*** adding 2 identical homogeneous layers                            ***
*** on entrance r and t contain the reflection and transmission      ***
*** matrices of the identical layers                                 ***
*** on exit r and t contain the reflection and transmission matrices ***
*** of the combined layer and u and d contain the up- and downward   ***
*** radiation at the interface (middle) of the two layers            ***
*** see de Haan et al. (1987) Eqs. (85)-(91)                         ***
************************************************************************
      implicit double precision (a-h,o-z)

      dimension ebmu(NDmu)
      dimension r(NDsup,NDsup), t(NDsup,NDsup)
      dimension u(NDsup,NDsup), d(NDsup,NDsup)
      dimension s(NDsup,NDsup), e(NDsup,NDsup), p(NDsup,NDsup)
cmv   common /mat/ s,e,p

      nsup = nmat*nmutot
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine DOUBLE too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop ' dimension error in DOUBLE'
      endif
************************************************************************
* calculate Q_1, Eq. (85)                                              *
************************************************************************
      do 10 j=1,nsup
      do 10 i=1,nsup
          p(i,j) = r(i,j)
   10 continue
      call doad_transhm(p,nmutot,nmat,NDsup)
      call doad_prodsm(e,p,r,nmutot,nmug,nmat,NDsup)
      do 20 j=1,nsup
      do 20 i=1,nsup
          s(i,j) = e(i,j)
   20 continue
      delta = 0.d0
      do 30 i=1,nmug
          delta = delta+e(i,i)
   30 continue
      if ((delta*delta) .lt. 1.d-12) goto 100
      ik = 0
************************************************************************
* productmethod for calculating S, Eqs. (112)-(115)                    *
************************************************************************
   40 call doad_prodsm(p,e,e,nmutot,nmug,nmat,NDsup)
      do 50 j=1,nsup
      do 50 i=1,nsup
          e(i,j) = p(i,j)
   50 continue
      call doad_prodsm(p,s,e,nmutot,nmug,nmat,NDsup)
************************************************************************
* P = S_r*C_r+1                                                        *
************************************************************************
      do 60 j=1,nsup
      do 60 i=1,nsup
          s(i,j) = p(i,j)+e(i,j)+s(i,j)
   60 continue
************************************************************************
* S = S_r+1                                                            *
************************************************************************
      delta = 0.d0
      do 70 i=1,nmug
          delta = delta+e(i,i)
   70 continue
************************************************************************
* delta=trace(c_r+1), Eq. (124), is an estimate of the error           *
************************************************************************
      ik = ik+1
      if (ik .gt. 20) then
          write(*,90)
          write(*,80)
          write(*,90)
   80     format('*** stop : S-series converge to slow ***')
   90     format('****************************************')
          stop 'slow convergence in DOUBLE'
      endif
      if ((delta*delta) .gt. 1.d-12) goto 40
************************************************************************
* calculate D, Eq. (88)                                                *
************************************************************************
  100 call doad_prodsm(p,s,t,nmutot,nmug,nmat,NDsup)
      do 110 i=1,nmutot
          im=nmat*(i-1)
          do 110 j=1,nmutot
              jm=nmat*(j-1)
              ej=ebmu(j)
              do 110 k=1,nmat
                  ik=im+k
                  do 110 l=1,nmat
                      jl=jm+l
                      d(ik,jl)=p(ik,jl)+ej*s(ik,jl)+t(ik,jl)
  110 continue
************************************************************************
* calculate U, Eq. (89)                                                *
************************************************************************
      call doad_prodsm(s,r,d,nmutot,nmug,nmat,NDsup)
      do 120 j=1,nmutot
          jm=nmat*(j-1)
          ej=ebmu(j)
          do 120 i=1,nmutot
              im=nmat*(i-1)
              do 120 k=1,nmat
                  ik=im+k
                  do 120 l=1,nmat
                      jl=jm+l
                      u(ik,jl)=s(ik,jl)+ej*r(ik,jl)
  120 continue
************************************************************************
* calculate T, Eq. (91)                                                *
************************************************************************
      call doad_prodsm(e,t,d,nmutot,nmug,nmat,NDsup)
      do 130 i=1,nmutot
          im=nmat*(i-1)
          ei=ebmu(i)
          do 130 j=1,i
              jm=nmat*(j-1)
              ej=ebmu(j)
              do 130 k=1,nmat
                  ik=im+k
                  do 130 l=1,nmat
                      jl=jm+l
                      e(ik,jl)=e(ik,jl)+ei*d(ik,jl)+ej*t(ik,jl)
  130 continue
************************************************************************
* calculate R, Eq. (90)                                                *
************************************************************************
      call doad_transhm(t,nmutot,nmat,NDsup)
      call doad_prodsm(p,t,u,nmutot,nmug,nmat,NDsup)
      do 140 i=1,nmutot
          im=nmat*(i-1)
          ei=ebmu(i)
          do 140 j=1,i
              jm=nmat*(j-1)
              do 140 k=1,nmat
                  ik=im+k
                  do 140 l=1,nmat
                      jl=jm+l
                      r(ik,jl)=p(ik,jl)+ei*u(ik,jl)+r(ik,jl)
  140 continue
      call doad_symmut(r,nmat,nmutot,NDsup,3)
      call doad_symmut(e,nmat,nmutot,NDsup,4)
      call doad_xeqy(t,e,NDsup)
      return
      end


      subroutine doad_expbmu(b,xmu,NDmu,nmutot,ebmu)
************************************************************************
*** fill array ebmu with dimension NDmu with exp[-b/xmu(i)]          ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), ebmu(NDmu)

      do 10 i=1,nmutot
        if (xmu(i) .gt. 1.d-10) then
          if (b/xmu(i) .lt. 2.d+2) then
            ebmu(i) = dexp(-b/xmu(i))
          else
            ebmu(i) = 0.d0
          endif
        else
          ebmu(i) = 0.d0
        endif
   10 continue
      return
      end


      subroutine doad_expZm(m,il,NDlay,alfbet,NDcoef,ncoef,xmu,wmu,NDmu,
     +    nmug,nmutot,nmat,Zmmin,Zmplus,NDsup,eps)
************************************************************************
*  Calculate the m-th Fourier component of the phase matrix Zm(mu0,mu) *
*  from the expansion coefficients of the scattering matrix into       *
*  generalized spherical functions. The formulae can be found in :     *
*                                                                      *
*    J.F. de Haan et al.: 1987, Astron. Astrophys. 183, pp. 371-391.   *
*                                                                      *
*  Essentially, Eqs. (66)-(82) are used. The suggestion below Eq. (82) *
*  to diagonalize the matriz Plm is followed here.                     *
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
************************************************************************
      implicit double precision (a-h,o-z)

      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension alfbet(4,4,0:NDcoef,NDlay), ncoef(NDlay)
      dimension xmu(NDmu), wmu(NDmu)
      dimension Plm(NDmu,3,2), DPDpl(NDsup), DPDmi(NDsup), DSD(4,4)
      dimension rootu(NDmu), sqlm(0:NDcoef), sql4(NDcoef) 

      if (nmat .eq. 1) then
          call doad_scalZm(m,il,NDlay,alfbet,NDcoef,ncoef,xmu,wmu,NDmu,
     +        nmug,nmutot,nmat,Zmmin,Zmplus,NDsup,eps)
          return
      endif
      if (nmutot .gt. NDmu) then
          print*,'dimension in subroutine EXPZM is too small'
          print*,'dim. is NDmu=',NDmu,' but nmutot=',nmutot
          stop 'dimension error in EXPZM'
      endif
      nsup = nmutot*nmat
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine EXPZM is too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop 'dimension error in EXPZM'
      endif
      if (ncoef(il) .gt. NDcoef) then
          print*,'dimension in subroutine EXPZM is too small'
          print*,'dim. is NDcoef=',NDcoef,' but ncoef=',ncoef(il)
          stop 'dimension error in EXPZM'
      endif
************************************************************************
* Initialize phase matrix to zero                                      *
************************************************************************
      do 100 j=1,nsup
      do 100 i=1,nsup
          Zmplus(i,j,il) = 0.d0
          Zmmin(i,j,il)  = 0.d0
  100 continue
      qroot6 = -0.25d0*dsqrt(6.d0)
************************************************************************
* Precompute the factor sqrt(l**2-m**2) needed in Eqs. (81)-(82)       *
* and also the factor sqrt(l**2-4) needed in Eq. (82)                  *
************************************************************************
      m2 = m*m
      do 10 l=m,ncoef(il)
          sqlm(l) = dsqrt(dble(l*l-m2))
   10 continue
      do 20 l=2,ncoef(il)
          sql4(l) = dsqrt(dble(l*l-4))
   20 continue
************************************************************************
* Compute the binomial factor 2m above m needed in Eq. (77)            *
* and also the binomial factor 2m above m-2 needed in Eq. (80)         *
************************************************************************
      bin2mm = 1.d0
      do 30 n=1,m
          bin2mm = bin2mm*dble(n+m)/dble(n)
   30 continue
      twom   = 2.d0**(-m)
      binfac = twom*dsqrt(bin2mm)
      if (m .ge. 2) then
          bin2m2 = bin2mm*dble(m*(m-1))/dble((m+1)*(m+2))
          binf2  = -twom*dsqrt(bin2m2)
      endif
************************************************************************
* Initialize Plm for l=m-1 to zero, Eq. (76)                           *
************************************************************************
      lold = 1
      lnew = 2
      do 40 k=1,3
      do 40 i=1,nmutot
          Plm(i,k,lold) = 0.d0
   40 continue
      do 50 i=1,nmutot
          rootu(i) = dsqrt(1.d0-xmu(i)*xmu(i))
          if (m .ne. 0) then
              Plm(i,1,lnew) = binfac*rootu(i)**m
          else
              Plm(i,1,lnew) = binfac
          endif
   50 continue
************************************************************************
* Start loop over l in loop 200 (summation index in Eq. (66))          *
* Parity of Plm is (-1)**(l-m)                                         *
************************************************************************
      parity = -1.d0
      do 200 l=m,ncoef(il)
        parity = -parity
************************************************************************
* Initialize Plm for l=m, Eqs.(77)-(80) without factor i**-m           *
* This factor is cancelled in Eq. (66) by the factor -1**m             *
* The exception for m=2 is needed to handle u=+-1 in Eq. (80)          *
************************************************************************
          if (l .eq. max0(m,2)) then
            if (m .eq. 0) then
              do 60 i=1,nmutot
                Plm(i,2,lnew) = qroot6*rootu(i)*rootu(i)
                Plm(i,3,lnew) = Plm(i,2,lnew)
   60         continue
            elseif (m .eq. 1) then
              do 70 i=1,nmutot
                u = xmu(i)
                Plm(i,2,lnew) = -0.5d0*rootu(i)*(1.d0-u)
                Plm(i,3,lnew) =  0.5d0*rootu(i)*(1.d0+u)
   70         continue
            elseif (m .eq. 2) then
              do 80 i=1,nmutot
                u = xmu(i)
                Plm(i,2,lnew) = -0.25d0*(1.d0-u)**2
                Plm(i,3,lnew) = -0.25d0*(1.d0+u)**2
   80         continue
            else
              do 90 i=1,nmutot
                u = xmu(i)
                urootm = rootu(i)**(m-2)
                Plm(i,2,lnew) = binf2*urootm*(1.d0-u)*(1.d0-u)
                Plm(i,3,lnew) = binf2*urootm*(1.d0+u)*(1.d0+u)
   90         continue
            endif
          endif
************************************************************************
* Construct supervectors corresponding to the diagonal elements of the *
* matrix D1 * Plm * D2 = D2 * Plm * D1                                 *
************************************************************************
          do 110 i=1,nmutot
              isup = nmat*(i-1)
              DPDpl(isup+1) = Plm(i,1,lnew)
              DPDpl(isup+2) = Plm(i,2,lnew)
              DPDpl(isup+3) = Plm(i,3,lnew)
              DPDmi(isup+1) = parity*Plm(i,1,lnew)
              DPDmi(isup+2) = parity*Plm(i,3,lnew)
              DPDmi(isup+3) = parity*Plm(i,2,lnew)
  110     continue
          if (nmat .eq. 4) then
              do 120 i=4,nsup,4
                  DPDpl(i) = DPDpl(i-3)
                  DPDmi(i) = DPDmi(i-3)
  120         continue
          endif
************************************************************************
* Construct the matrix D2 * S * D2                                     *
************************************************************************
          DSD(1,1) = alfbet(1,1,l,il)
          DSD(2,1) = 0.5d0*alfbet(1,2,l,il)
          DSD(2,2) = 0.25d0*(alfbet(2,2,l,il)+alfbet(3,3,l,il))
          DSD(3,2) = 0.25d0*(alfbet(2,2,l,il)-alfbet(3,3,l,il))
          DSD(3,1) = DSD(2,1)
          DSD(1,2) = DSD(2,1)
          DSD(1,3) = DSD(2,1)
          DSD(2,3) = DSD(3,2)
          DSD(3,3) = DSD(2,2)
          if (nmat .eq. 4) then
              DSD(1,4) = 0.d0
              DSD(2,4) = 0.5d0*alfbet(3,4,l,il)
              DSD(3,4) = -DSD(2,4)
              DSD(4,4) = alfbet(4,4,l,il)
              DSD(4,1) = 0.d0
              DSD(4,2) = -DSD(2,4)
              DSD(4,3) = -DSD(3,4)
          endif
************************************************************************
* Add a new term to the sum in Eq. (66)                                *
* The factor (-1)**m is cancelled by i**m in the Plm                   *
************************************************************************
          do 140 k2=1,nmat
          do 140 k1=1,nmat
          do 140 j=k2,nsup,nmat
              SPj = DSD(k1,k2)*DPDpl(j)
              do 130 i=k1,nsup,nmat
                  Zmplus(i,j,il) = Zmplus(i,j,il)+DPDpl(i)*SPj
                  Zmmin(i,j,il) = Zmmin(i,j,il)+DPDmi(i)*SPj
  130         continue
  140     continue
************************************************************************
* Do one step in recurrence for Plm0, Eq.(81)                          *
************************************************************************
          if (l .eq. ncoef(il)) goto 200
          twol1 = dble(2*l+1)
          f1new = twol1/sqlm(l+1)
          f1old = sqlm(l)/sqlm(l+1)
          do 150 i=1,nmutot
              u = xmu(i)
              Plm(i,1,lold) = f1new*u*Plm(i,1,lnew)
     +                          -f1old*Plm(i,1,lold)
  150     continue
************************************************************************
* Do one step in recurrence for Plm2, Eq.(82)                          *
* only when they have been initialized : l >= max(m,2)                 *
************************************************************************
          if (l .ge. max0(m,2)) then
              tmp   = 1.d0/(dble(l)*sql4(l+1)*sqlm(l+1))
              f2new = twol1*dble(l*(l+1))*tmp
              f2newa= twol1*2.d0*dble(m)*tmp
              f2old = dble(l+1)*sql4(l)*sqlm(l)*tmp
              do 160 i=1,nmutot
                  u = xmu(i)
                  Plm(i,2,lold) = (f2new*u+f2newa)*Plm(i,2,lnew)
     +                          -f2old*Plm(i,2,lold)
                  Plm(i,3,lold) = (f2new*u-f2newa)*Plm(i,3,lnew)
     +                          -f2old*Plm(i,3,lold)
  160         continue
          endif
          ltmp = lnew
          lnew = lold
          lold = ltmp
  200 continue
************************************************************************
* End of summation loop over l                                         *
* Calculate D1 * sum * D1                                              *
************************************************************************
      call doad_transf(Zmmin,il,NDlay,NDsup,nsup,nmat)
      call doad_transf(Zmplus,il,NDlay,NDsup,nsup,nmat)
************************************************************************
* renormalization of the phase matrix for m=0                          *
************************************************************************
      if (m .eq. 0) then
        a = alfbet(1,1,0,il)
        call doad_renorm(a,Zmmin,Zmplus,nmug,nmutot,nmat,xmu,wmu,il,
     +    NDlay,NDmu,NDsup,eps)
      endif
      return
      end


      subroutine doad_fillar(b,NDlay,tau,taun,NDirf,nlayer,nlirf,
     +     xmu,NDmu,
     +     nmutot,en,eb,ek,etmu,ebtmu)
************************************************************************
*** fill arrays tau, taun, en, eb, ek, etmu and ebtmu                ***
*** b(n) contains the optical thickness of isolated layer n          ***
*** tau(0,n) contains the optical depth at interface n in atmosphere ***
*** tau(0,1)=b is the groundlevel and tau(0,nlayer+1)=0 the top      ***
*** tau(k,n) contains the optical depth at level k in layer n        ***
*** tau(nlirf(n)+1,n) contains the optical depth at interface n+1    ***
*** taun(k,n) the optical depth at level k in isolated layer n       ***
*** en(i,n) contains exp[-tau(0,n)/xmu(i)] = E(tau_n,N)              ***
*** eb(i,n) contains exp[-b(n)/xmu(i)] = E(tau_n,N-tau_n+1,N)        ***
*** ek(i,k,n) contains                                               ***
*** exp[-b(n)/xmu(i)*2**(k-M-1)] = E(tau_n^k,K-tau_n^k+1,K) if k<=M  ***
*** exp[-b(n)/xmu(i)*2**(-k+M-1)] = E(tau_n^k,K)            if k>=M  ***
*** etmu(i,k,n) contains exp(-taun(k,n)/xmu(i))                      ***
*** ebtmu(i,k,n) contains exp((taun(k,n)-taun(0,n))/xmu(i))          ***
*** see Stammes et al. (1989)                                        ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension en(NDmu,NDlay+1), eb(NDmu,NDlay), ek(NDmu,NDirf,NDlay)
      dimension nlirf(0:NDlay+1), b(0:NDlay), tau(0:NDirf+1,0:NDlay+1)
      dimension xmu(NDmu), taun(0:NDirf+1,0:NDlay)
      dimension etmu(NDmu,0:NDirf+1,0:NDlay)
      dimension ebtmu(NDmu,0:NDirf+1,NDlay)

      do 10 n=0,nlayer+1
        tau(0,n) = 0.d0
   10 continue
      do 20 n=1,nlayer
        tau(0,1) = tau(0,1)+b(n)
   20 continue
      tau(0,0) = tau(0,1)
      tau(1,0) = tau(0,1)
      do 30 n=2,nlayer
        tau(0,n) = tau(0,n-1)-b(n-1)
   30 continue
      do 40 n=1,nlayer
      do 40 i=1,nmutot
        if (xmu(i) .gt. 1.d-10) then
          if (tau(0,n)/xmu(i) .lt. 2.d+2) then
            en(i,n) = dexp(-tau(0,n)/xmu(i))
          else
            en(i,n) = 0.d0
          endif
        else
          en(i,n) = 0.d0
        endif
   40 continue
      do 45 i=1,nmutot
        en(i,nlayer+1) = 1.d0
   45 continue
      do 50 n=1,nlayer
      do 50 i=1,nmutot
        if (xmu(i) .gt. 1.d-10) then
          if (b(n)/xmu(i) .lt. 2.d+2) then
            eb(i,n) = dexp(-b(n)/xmu(i))
          else
            eb(i,n) = 0.d0
          endif
        else
          eb(i,n) = 0.d0
        endif
   50 continue
      do 80 n=1,nlayer
        M = (nlirf(n)+1)/2
        do 70 k=1,M
          kup = nlirf(n)-k+1
          do 60 i=1,nmutot
            if (xmu(i) .gt. 1.d-10) then
              if (b(n)/xmu(i)*2.d0**(k-M-1) .lt. 2.d+2) then
                ek(i,k,n) = dexp(-b(n)/xmu(i)*2.d0**(k-M-1))
                ek(i,kup,n) = dexp(-b(n)/xmu(i)*2.d0**(k-M-1))
              else
                ek(i,k,n) = 0.d0
                ek(i,kup,n) = 0.d0
              endif
            else
              ek(i,k,n) = 0.d0
              ek(i,kup,n) = 0.d0
            endif
   60     continue
   70   continue
   80 continue
      do 100 n=1,nlayer
        M = (nlirf(n)+1)/2
        do 90 k=1,M
          kup = nlirf(n)-k+1
          tau(k,n) = tau(0,n+1)+b(n)*(1.d0-2.d0**(k-M-1))
          tau(kup,n) = tau(0,n+1)+b(n)*2.d0**(-kup+M-1)
   90   continue
        tau(nlirf(n)+1,n) = tau(0,n+1)
  100 continue
      do 110 n=0,nlayer
        M =(nlirf(n)+1)/2
        taun(0,n) = b(n)
        taun(nlirf(n)+1,n) = 0.d0
        do 110 k=1,M
          kup = nlirf(n)-k+1
          taun(k,n) = b(n)*(1.d0-2.d0**(k-M-1))
          taun(kup,n) = b(n)*(2.d0**(M-kup-1))
  110 continue
      do 130 n=1,nlayer
        do 120 k=0,nlirf(n)
          do 120 i=1,nmutot
            if (xmu(i) .gt. 1.d-10) then
              if (taun(k,n)/xmu(i) .lt. 2.d+2) then
                etmu(i,k,n) = dexp(-taun(k,n)/xmu(i))
              else
                etmu(i,k,n) = 0.d0
              endif
            else
              etmu(i,k,n) = 0.d0
            endif
  120   continue
        do 130 i=1,nmutot
          etmu(i,nlirf(n)+1,n) = 1.d0
  130 continue
      do 140 i=1,nmutot
        etmu(i,0,0) = 1.d0
        etmu(i,1,0) = 1.d0
  140 continue
      do 160 n=1,nlayer
        do 150 k=1,nlirf(n)+1
          do 150 i=1,nmutot
            if (xmu(i) .gt. 1.d-10) then
              if ((taun(k,n)-taun(0,n))/xmu(i) .lt. 2.d+2) then
                ebtmu(i,k,n) = dexp((taun(k,n)-taun(0,n))/xmu(i))
              else
                ebtmu(i,k,n) = 0.d0
              endif
            else
              ebtmu(i,k,n) = 0.d0
            endif
  150   continue
        do 160 i=1,nmutot
          ebtmu(i,0,n) = 1.d0
  160 continue
      return
      end


      subroutine doad_firstm(nmat,nmutot,xmu,wmu,b,ebmu,m,il,NDlay,
     +    Zmmin,Zmplus,NDmu,NDsup,r,t)
************************************************************************
*** calculate the first order m-th fourier coefficient of the        ***
*** reflection and transmission matrices r and t in supermatrix form ***
*** see Hovenier (1971)                                              ***
*** the m-th fourier coefficient of the phasematrix of the layer is  ***
*** contained in the arrays Zmmin(i,j,il) and Zmplus(i,j,il)         ***
*** ebmu contains exp[-b/mu(i)] with b the optical thickness of the  ***
*** layer (starting, or total layer in case of single scat. only)    ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu), ebmu(NDmu)
      dimension r(NDsup,NDsup), t(NDsup,NDsup)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)

      ag4 = .25d0
************************************************************************
* loop 20 (i,j) is the loop over each matrix element contained in      *
* the lower triangle of the supermatrix                                *
************************************************************************
      do 20 i=1,nmutot
          xi = xmu(i)
          ei = ebmu(i)
          awi = ag4*wmu(i)
          im = nmat*(i-1)
          do 20 j=1,i
              xj = xmu(j)
              ej = ebmu(j)
              awij = awi*wmu(j)
              jm = nmat*(j-1)
              eiej = 1.d0-ei*ej
              if (eiej .gt. 1.d-3) then
                  z = xi+xj
                  if (z .gt. 1.d-10) z = 1.d0/z
************************************************************************
* Eq. (A-21)                                                           *
************************************************************************
                  yr = awij*z*eiej
                  if (dabs(xi-xj) .gt. 1.d-10) then
************************************************************************
* Eq. (A-22)                                                           *
************************************************************************
                      yt = awij*(ei-ej)/(xi-xj)
                  else
************************************************************************
* special case of the transmission matrix if mu(i)=mu(j)  Eq. (A-23)   *
* the reflection matrix is not defined for mu(i)=mu(j)=0               *
************************************************************************
                      z = 0.d0
                      if (xi .gt. 1.d-10) z = b*ei/xi/xi
                      yt = awij*z
                  endif
              else
************************************************************************
* expansion in taylor series to avoid loss of accuracy if b<<1         *
************************************************************************
                  bgi = b/xi
                  bgj = b/xj
                  bij = bgi+bgj
                  z = (1.d0-.5d0*bij*(1.d0-bij/3.d0))
                  y = awij*bgi/xj
                  yr = y*z
                  yt = y*(z-bgi*bgj/6.d0)
              endif
              do 10 l=1,nmat
                  jl = jm+l
                  do 10 k=1,nmat
                      ik = im+k
                      r(ik,jl) = yr*Zmmin(ik,jl,il)
                      t(ik,jl) = yt*Zmplus(ik,jl,il)
   10         continue
   20 continue
************************************************************************
* fill the upper triangle using symmetry relations                     *
************************************************************************
      call doad_symmut(r,nmat,nmutot,NDsup,3)
      call doad_symmut(t,nmat,nmutot,NDsup,4)
      return
      end

      subroutine doad_fouMi(M0,M1,M2,nlayer,NDlay,alfbet,NDcoef,
     +                      ncoef,eps)
************************************************************************
*** all orders of scatt. for  0 < m <= M2                            ***
*** two orders of scatt. for M2 < m <= M1                            ***
*** single scattering    for M1 < m <= M0                            ***
*** see de Haan et. al. (1987) page 386                              ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension alfbet(4,4,0:NDcoef,NDlay), ncoef(NDlay)
      dimension M0(NDlay), M1(NDlay), M2(NDlay)

      do 20 il=1,nlayer
          M0(il) = 0
          M1(il) = 0
          M2(il) = 0
          do 10 l=0,ncoef(il)
************************************************************************
* see Eqs. (138)-(140)                                                 *
************************************************************************
              x = dabs(alfbet(1,1,l,il))
              x0 = x*.25d0
              x1 = x0*x/(2.d0*l+1.d0)
              x2 = x1*x/(2.d0*l+1.d0)
              if (x0 .gt. eps) M0(il) = l
              if (x1 .gt. eps) M1(il) = l
              if (x2 .gt. eps) M2(il) = l
   10     continue
   20 continue
      return
      end


      subroutine doad_gauleg(ndim,ngauss,a,b,x,w,eps)
************************************************************************
*** given the lower and upper limits of integration a and b,         ***
*** and given the number of gauss-legendre points ,ngauss,           ***
*** this routine returns through array x the abscissas and           ***
*** through array w the weights of the Gauss-Legendre                ***
*** quadrature formula.                                              ***
*** eps is the desired accuracy of the abscissas.                    ***
*** this subroutine is documented further in :                       ***
*** W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr.        ***
*** (1987) page 125. ISBN 0-521-30811-9.                             ***
*** refenrences are given to formulae in A&S:                        ***
*** M. Abramowitz, I.A. Stegun 'Handbook of Mathematical             ***
*** Functions' Dover Pub. Inc. (1965)                                ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension x(ndim), w(ndim)

      pi = 4.d0*datan(1.d0)
      m = (ngauss+1)/2
      xmid = 0.5d0*(a+b)
      xlen = 0.5d0*(b-a)
      do 30 i=1,m
************************************************************************
* estimate of the root of Legendre polynomial of degree ngauss         *
************************************************************************
          z = dcos(pi*(i-0.25d0)/(ngauss+0.5d0))
   10     continue
          pn1 = 1.d0
          pn2 = 0.d0
          do 20 n=1,ngauss
************************************************************************
* A&S p.782 recurrence relation of Legendre polynomials                *
************************************************************************
              pn3 = pn2
              pn2 = pn1
              pn1 = ((2.d0*n-1.d0)*z*pn2-(n-1.d0)*pn3)/dble(n)
   20     continue
************************************************************************
* A&S p.783 deriviate of Legendre polynomials                          *
************************************************************************
          pnder = ngauss*(z*pn1-pn2)/(z*z-1.d0)
************************************************************************
* A&S p.18 newton-raphson iteration to root                            *
************************************************************************
          zold = z
          z = zold-pn1/pnder
          if (dabs(z-zold) .gt. eps) goto 10
************************************************************************
* A&S p.887 abscissa and root on interval a,b                          *
************************************************************************
          x(i) = xmid-xlen*z
          x(ngauss+1-i) = xmid+xlen*z
c	  print *,'x(',i,') = ',x(i),'  x(',ngauss+1-i,') = ',
c     &	  	x(ngauss+1-i)
          w(i) = 2.d0*xlen/((1.d0-z*z)*pnder*pnder)
          w(ngauss+1-i) = w(i)
   30 continue
      return
      end



      subroutine doad_homirf(nlirf,nl,nmat,nmutot,nmug,ukk,dkk,ek,
     +    NDsup,NDirf,NDlay,NDmu)
************************************************************************
*** internal radiation field of the isolated homogeneous layer nl    ***
*** the method is described in:                                      ***
*** Stammes, P., de Haan, J.F., Hovenier, J.W.:                      ***
*** 1989, Astron. Astrophys. xxx, p. xxx                             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension ukk(NDsup,NDsup,NDirf,NDlay)
      dimension dkk(NDsup,NDsup,NDirf,NDlay)
      dimension ek(NDmu,NDirf,NDlay), nlirf(0:NDlay+1)
************************************************************************
* on entrance the matrix ukk(i,j,k,n) contains the supermatrix U, at   *
* middle (level k) of a partial layer n build in k doubling steps      *
* on exit ukk(i,j,k,n) contains the supermatrix U, at level k in the   *
* entire isolated homogeneous layer n                                  *
* ek(i,k,n) contains the factor exp[-b(n)/u(i)*2**(k-M-1)], where      *
* b(n) is the optical thickness of layer n                             *
************************************************************************

      dimension uk(NDsup,NDsup), dk(NDsup,NDsup)
      dimension ukstar(NDsup,NDsup), dkstar(NDsup,NDsup)
      dimension ukp1(NDsup,NDsup), dkp1(NDsup,NDsup)
      dimension ukm1(NDsup,NDsup)
      dimension ekp1(NDmu), em(NDmu), e(NDmu)
      dimension x(NDsup,NDsup), y(NDsup,NDsup), z(NDsup,NDsup)
cmv   common /mat/ x,y,z,ukp1,dkp1,uk,dk,ukstar,dkstar,ukm1

      if (nlirf(nl) .eq. 0) return
      if (nmutot .gt. NDmu) then
          print*,'dimension in subroutine HOMIRF too small'
          print*,'dim. is NDmu=',NDmu,' but nmutot=',nmutot
          stop ' dimension error in HOMIRF'
      endif
      nsup = nmat*nmutot
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine HOMIRF too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop ' dimension error in HOMIRF'
      endif
      nl2 = (nlirf(nl)+1)/2
      do 80 k=nl2-1,1,-1
          kup = nlirf(nl)+1-k
************************************************************************
* initializing arrays                                                  *
************************************************************************
          do 10 i=1,nmutot
              em(i) = ek(i,k,nl)
              ekp1(i) = ek(i,k+1,nl)
              e(i) = ekp1(i)**(2.d0**(nl2-k)-1.d0)
   10     continue
          do 20 i=1,nsup
          do 20 j=1,nsup
              ukp1(j,i) = ukk(j,i,k+1,nl)
              dkp1(j,i) = dkk(j,i,k+1,nl)
              uk(j,i) = ukk(j,i,k,nl)
              dk(j,i) = dkk(j,i,k,nl)
   20     continue
************************************************************************
* internal radiation field of lower half of isolated layer             *
* D field, Eq. (64)                                                    *
************************************************************************
          call doad_prodex(x,dkp1,em,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_prodsm(y,dk,dkp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(z,x,y,nsup,NDsup)
          call doad_prodex(x,dk,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_addsm(y,x,z,nsup,NDsup)
          do 30 i=1,nsup
          do 30 j=1,nsup
              dkk(j,i,k,nl) = y(j,i)
   30     continue
************************************************************************
* U field, Eq. (65)                                                    *
************************************************************************
          call doad_prodex(x,uk,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_prodsm(y,uk,dkp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(z,x,y,nsup,NDsup)
          do 40 i=1,nsup
          do 40 j=1,nsup
              ukk(j,i,k,nl) = z(j,i)
   40     continue
************************************************************************
* internal radiation field of upper half of isolated layer             *
************************************************************************
          do 50 i=1,nsup
          do 50 j=1,nsup
              ukm1(j,i) = ukk(j,i,kup-1,nl)
   50     continue
************************************************************************
* calculate U*kk and D*kk according to Eqs. (68) and (69)              *
************************************************************************
          call doad_xeqy(ukstar,uk,NDsup)
          call doad_xeqy(dkstar,dk,NDsup)
          call doad_transhm(ukstar,nmutot,nmat,NDsup)
          call doad_transhm(dkstar,nmutot,nmat,NDsup)
************************************************************************
* D field, Eq. (70)                                                    *
************************************************************************
          call doad_prodsm(z,ukstar,ukm1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(x,z,dk,nsup,NDsup)
          do 60 i=1,nsup
          do 60 j=1,nsup
              dkk(j,i,kup,nl) = x(j,i)
   60     continue
************************************************************************
* U field, Eq. (71)                                                    *
************************************************************************
          call doad_prodex(z,ukm1,em,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_prodsm(x,dkstar,ukm1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(y,x,z,nsup,NDsup)
          call doad_addsm(z,y,uk,nsup,NDsup)
          do 70 i=1,nsup
          do 70 j=1,nsup
              ukk(j,i,kup,nl) = z(j,i)
   70     continue
   80 continue
      return
      end


      subroutine doad_homlay(nmat,nmutot,nmug,xmu,wmu,b,m,il,ebmu,
     +    r,t,nlirf,ndoubl,ukk,dkk,Zmmin,Zmplus,
     +    NDsup,NDirf,NDmu,NDlay)
************************************************************************
*** calculate reflection and transmission r and t of an isolated     ***
*** homogeneous layer including all orders of scatt.  0 < m <= M2    ***
*** see de Haan et al. (1987) page 387                               ***
************************************************************************
      implicit double precision (a-h,o-z)

      dimension xmu(NDmu), wmu(NDmu), b(0:NDlay), ebmu(NDmu)
      dimension ndoubl(NDlay), nlirf(0:NDlay+1)
      dimension ukk(NDsup,NDsup,NDirf,NDlay)
      dimension dkk(NDsup,NDsup,NDirf,NDlay)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension r(NDsup,NDsup), t(NDsup,NDsup)
      dimension u(NDsup,NDsup), d(NDsup,NDsup)
cmv   common /mat/ u,d

      nsup = nmat*nmutot
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine HOMLAY too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop 'dimension error in HOMLAY'
      endif
      bs = b(il)*2.d0**(-ndoubl(il))
      call doad_expbmu(bs,xmu,NDmu,nmutot,ebmu)
cmv   write(*,*) 'fourier term = ',m,' layer il = ',il
cmv   write(*,*) 'bstart = ',bs,' ndoubl(il) = ',ndoubl(il)
      call doad_firstm(nmat,nmutot,xmu,wmu,bs,ebmu,m,il,NDlay,
     +    Zmmin,Zmplus,NDmu,NDsup,r,t)
      call doad_seconm(nmat,nmutot,nmug,xmu,wmu,bs,ebmu,m,il,
     +    r,t,Zmmin,Zmplus,NDmu,NDsup,NDlay)
      nl2 = (nlirf(il)+1)/2
      do 40 n=1,ndoubl(il)
          call doad_double(nmat,nmutot,nmug,NDmu,NDsup,ebmu,r,t,u,d)
************************************************************************
* fill arrays ukk and dkk at all irf-levels in the isolated homogeneous*
* layer start when b=bs*2**(ndoubl(il)-nl2+1)                          *
* and calculate the reflection and transmission, r and t of that layer *
************************************************************************
          if ((n.gt.(ndoubl(il)-nl2)) .and. (nlirf(il).ne.0)) then
              nn = n-(ndoubl(il)-nl2)
              do 30 j=1,nsup
              do 30 i=1,nsup
                  ukk(i,j,nn,il) = u(i,j)
                  dkk(i,j,nn,il) = d(i,j)
   30         continue
          endif
          bs = 2.d0*bs
          call doad_expbmu(bs,xmu,NDmu,nmutot,ebmu)
   40 continue
      return
      end


      subroutine doad_initsm(s,NDsup)
************************************************************************
*** initialize the supermatrix s to zero                             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension s(NDsup,NDsup)

      do 10 j=1,NDsup
      do 10 i=1,NDsup
          s(i,j) = 0.d0
   10 continue
      return
      end


      subroutine doad_irfemb(nlirf,nl,nmat,nmutot,nmug,ukk,dkk,
     +    unn,dnn,ek,
     +    en,NDsup,NDirf,NDlay,NDmu)
************************************************************************
*** internal radiation field of an embedded homogeneous layer nl     ***
*** the method is described in:                                      ***
*** Stammes, P., de Haan, J.F., Hovenier, J.W.:                      ***
*** 1989, Astron. Astrophys. xxx, p. xxx                             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension ukk(NDsup,NDsup,NDirf,NDlay), unn(NDsup,NDsup,NDlay+1)
      dimension dkk(NDsup,NDsup,NDirf,NDlay), dnn(NDsup,NDsup,NDlay+1)
      dimension ek(NDmu,NDirf,NDlay), en(NDmu,NDlay+1), nlirf(0:NDlay+1)
************************************************************************
* on entrance the matrix ukk(i,j,k,n) contains the supermatrix U, at   *
* level k in the isolated homogeneous layer n                          *
* on exit ukk(i,j,k,n) contains the supermatrix U, at level k in the   *
* homogeneous layer n of the complete atmosphere                       *
* unn(i,j,n) contains the supermatrix U, at level n in the complete    *
* atmosphere                                                           *
* ek(i,k,n) contains the factor exp[-b(n)/u(i)*2**(k-M-1)] if k <= M   *
* and the factor exp[-b(n)/u(i)*2**(-k+M-1)] if k > M                  *
* en(i,n) contains the factor exp[-tau(0,n)/u(i)]                      *
* where b(n) is the optical thickness of layer n                       *
* and tau(0,n) is the optical depth at level n                         *
************************************************************************

      dimension uiso(NDsup,NDsup,NDirf), diso(NDsup,NDsup,NDirf)
      dimension uk(NDsup,NDsup), dk(NDsup,NDsup)
      dimension un(NDsup,NDsup), dnp1(NDsup,NDsup)
      dimension ukstar(NDsup,NDsup), dkstar(NDsup,NDsup)
      dimension eu(NDmu), ed(NDmu), e(NDmu)
      dimension x(NDsup,NDsup), y(NDsup,NDsup), z(NDsup,NDsup)
cmv   common /mat/ x,y,z,dnp1,un,uk,dk,ukstar,dkstar,uiso,diso

      if (nlirf(nl) .eq. 0) return
      if (nmutot .gt. NDmu) then
          print*,'dimension in subroutine IRFEMB too small'
          print*,'dim. is NDmu=',NDmu,' but nmutot=',nmutot
          stop ' dimension error in IRFEMB'
      endif
      nsup = nmat*nmutot
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine IRFEMB too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop ' dimension error in IRFEMB'
      endif
      if (nlirf(nl) .gt. NDirf) then
          print*,'dimension in subroutine IRFEMB too small'
          print*,'dim. is NDirf=',NDirf,' but nlirf(nl)=',nlirf(nl)
          stop ' dimension error in IRFEMB'
      endif
      do 5 k=1,nlirf(nl)
      do 5 i=1,nsup
      do 5 j=1,nsup
          uiso(j,i,k) = ukk(j,i,k,nl)
          diso(j,i,k) = dkk(j,i,k,nl)
    5 continue
      nl2 = (nlirf(nl)+1)/2
      do 50 k=1,nlirf(nl)
          kstar = nlirf(nl)+1-k
************************************************************************
* initializing arrays                                                  *
************************************************************************
          do 10 i=1,nmutot
              e(i) = en(i,nl+1)
              if (k .le. nl2) then
                  ed(i) = ek(i,k,nl)**(2.d0**(nl2-k+1)-1.d0)
                  eu(i) = ek(i,k,nl)
              else
                  ed(i) = ek(i,k,nl)
                  eu(i) = ek(i,k,nl)**(2.d0**(-nl2+k+1)-1.d0)
              endif
   10     continue
          do 20 i=1,nsup
          do 20 j=1,nsup
              dnp1(j,i) = dnn(j,i,nl+1)
              un(j,i) = unn(j,i,nl)
              uk(j,i) = uiso(j,i,k)
              dk(j,i) = diso(j,i,k)
              ukstar(j,i) = uiso(j,i,kstar)
              dkstar(j,i) = diso(j,i,kstar)
   20     continue
************************************************************************
* internal radiation field of embedded homogeneous layer nl            *
* D field, Eq. (72)                                                    *
************************************************************************
          call doad_transhm(ukstar,nmutot,nmat,NDsup)
          call doad_prodsm(x,ukstar,un,nmutot,nmug,nmat,NDsup)
          call doad_prodex(y,dnp1,ed,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_addsm(z,x,y,nsup,NDsup)
          call doad_prodsm(x,dk,dnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(y,z,x,nsup,NDsup)
          call doad_prodex(x,dk,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_addsm(z,y,x,nsup,NDsup)
          do 30 i=1,nsup
          do 30 j=1,nsup
              dkk(j,i,k,nl) = z(j,i)
   30     continue
************************************************************************
* U field, Eq. (73)                                                    *
************************************************************************
          call doad_transhm(dkstar,nmutot,nmat,NDsup)
          call doad_prodsm(x,dkstar,un,nmutot,nmug,nmat,NDsup)
          call doad_prodex(y,un,eu,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_addsm(z,x,y,nsup,NDsup)
          call doad_prodsm(x,uk,dnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(y,z,x,nsup,NDsup)
          call doad_prodex(x,uk,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_addsm(z,y,x,nsup,NDsup)
          do 40 i=1,nsup
          do 40 j=1,nsup
              ukk(j,i,k,nl) = z(j,i)
   40     continue
   50 continue
      return
      end


      subroutine doad_maxM(M0,nlayer,NDlay,Mmax)
************************************************************************
*** determine the maximal M0 value, Mmax                               *
************************************************************************
      dimension M0(NDlay)

      Mmax = 0
      do 10 il=1,nlayer
          if (M0(il) .gt. Mmax) Mmax = M0(il)
   10 continue
      return
      end


      subroutine doad_mullay(nlayer,nmat,nmutot,nmug,unn,dnn,en,eb,
     +    NDsup,NDlay,NDmu)
************************************************************************
*** radiation at the interfaces of a multilayered atmosphere         ***
*** the method is described in:                                      ***
*** Stammes, P., de Haan, J.F., Hovenier, J.W.:                      ***
*** 1989, Astron. Astrophys. xxx, p. xxx                             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension unn(NDsup,NDsup,NDlay+1)
      dimension dnn(NDsup,NDsup,NDlay+1)
      dimension en(NDmu,NDlay+1), eb(NDmu,NDlay)
************************************************************************
* on entrance the matrix unn(i,j,n) contains the supermatrix U, at     *
* level n of a partial atmosphere build in n adding steps              *
* on exit unn(i,j,n) contains the supermatrix U, at level n in the     *
* complete atmosphere                                                  *
* en(i,n) contains the factor exp[-tau(0,n)/u(i)]                      *
* eb(i,n) contains the factor exp[-b(n)/u(i)]                          *
* where, b(n) is the optical thickness of layer n                      *
* and tau(0,n) is the optical depth at level n                         *
************************************************************************

      dimension un(NDsup,NDsup), dn(NDsup,NDsup)
      dimension dnp1(NDsup,NDsup)
      dimension em(NDmu), e(NDmu)
      dimension x(NDsup,NDsup), y(NDsup,NDsup), z(NDsup,NDsup)
cmv   common /mat/ x,y,z,dnp1,un,dn

      if (nlayer .le. 1) return
      nsup = nmat*nmutot
      if (nmutot .gt. NDmu) then
          print*,'dimension in subroutine MULLAY too small'
          print*,'dim. is NDmu=',NDmu,' but nmutot=',nmutot
          stop ' dimension error in MULLAY'
      endif
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine MULLAY too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop ' dimension error in MULLAY'
      endif
      do 50 n=nlayer-1,1,-1
************************************************************************
* initializing arrays                                                  *
************************************************************************
          do 10 i=1,nmutot
              e(i) = en(i,n+1)
              em(i) = eb(i,n)
   10     continue
          do 20 i=1,nsup
          do 20 j=1,nsup
              dnp1(j,i) = dnn(j,i,n+1)
              un(j,i) = unn(j,i,n)
              dn(j,i) = dnn(j,i,n)
   20     continue
************************************************************************
* diffuse radiation field in the complete atmosphere                   *
* D field, Eq. (62)                                                    *
************************************************************************
          call doad_prodex(x,dnp1,em,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_prodsm(y,dn,dnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(z,x,y,nsup,NDsup)
          call doad_prodex(x,dn,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_addsm(y,x,z,nsup,NDsup)
          do 30 i=1,nsup
          do 30 j=1,nsup
              dnn(j,i,n) = y(j,i)
   30     continue
************************************************************************
* U field, Eq. (63)                                                    *
************************************************************************
          call doad_prodex(x,un,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_prodsm(y,un,dnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(z,x,y,nsup,NDsup)
          do 40 i=1,nsup
          do 40 j=1,nsup
              unn(j,i,n) = z(j,i)
   40     continue
   50 continue
      return
      end

      subroutine doad_oneord(nmat,nmutot,xmu,wmu,taun,etmu,m,
     +    Zmmin,Zmplus,
     +    ebtmu,nlirf,n,NDsup,NDmu,NDlay,NDirf,r1,t1)
************************************************************************
*** calculate the first order internal radiation field of the m-th   ***
*** fourier coefficient of an isolated homogeneous layer             ***
*** see Hovenier (1971) and Wauben (1989)                            ***
*** the m-th fourier coefficient of the phasematrix of layer n is    ***
*** contained in the arrays Zmmin(j,i,n) and Zmplus(j,i,n)           ***
*** etmu contains exp[-taun(kn,n)/xmu(i)] with taun(kn,n) the optical***
*** depth of level kn in isolated layer n                            ***
*** ebtmu contains exp[-(taun(0,n)-taun(kn,n))/xmu(i)]               ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu), nlirf(0:NDlay+1)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension etmu(NDmu,0:NDirf+1,0:NDlay)
      dimension ebtmu(NDmu,0:NDirf+1,NDlay)
      dimension taun(0:NDirf+1,0:NDlay)
      dimension r1(NDsup,NDsup,0:NDirf+1,NDlay)
      dimension t1(NDsup,NDsup,0:NDirf+1,NDlay)

      ag4 = .25d0
************************************************************************
* loop 10 (j,i,kn) is the loop over each matrix element (j,i) and      *
* for all the irf levels kn and the layer n                            *
* determine the reflection and transmission function yr and yt, see    *
* Wauben (1989) Eqs. (3)-(6) with mu0->i and mu->j                     *
************************************************************************
      do 10 kn=0,nlirf(n)+1
        do 10 j=1,nmutot
          jm = nmat*(j-1)
          xj = xmu(j)
          etj = etmu(j,kn,n)
          ebtj = ebtmu(j,kn,n)
          awj = ag4*wmu(j)
          do 10 i=1,nmutot
            im = nmat*(i-1)
            xi = xmu(i)
            eti = etmu(i,kn,n)
            ebti = ebtmu(i,kn,n)
            awji = awj*wmu(i)
            xsum = xi+xj
            if (dabs(xi-xj) .gt. 1.d-10) then
************************************************************************
* Eq. (3)                                                              *
************************************************************************
              yt = awji*(eti-etj)/(xi-xj)
            else
************************************************************************
* special case of the transmission matrix if mu(i)=mu(j)  Eq. (4)      *
************************************************************************
              yt = 0.d0
              if (xi .gt. 1.d-10) yt = awji/xi*taun(kn,n)*eti/xi
            endif
            if (xsum .gt. 1.d-10) then
************************************************************************
* Eq. (5)                                                              *
************************************************************************
              yr = awji*eti/xsum*(1.d0-ebti*ebtj)
            else
************************************************************************
* the reflection matrix is not defined for mu(i)=mu(j)=0, Eq.(6)       *
************************************************************************
              yr = 0.d0
            endif
            do 10 k=1,nmat
              ik = im+k
              do 10 l=1,nmat
                jl = jm+l
                r1(jl,ik,kn,n) = yr*Zmmin(jl,ik,n)
                t1(jl,ik,kn,n) = yt*Zmplus(jl,ik,n)
   10 continue
      return
      end


      subroutine doad_ord1m(nmat,nmutot,xmu,wmu,taun,etmu,m,
     +    Zmmin,Zmplus,
     +    ebtmu,nlirf,il,NDsup,NDmu,NDlay,NDirf,r1,t1)
************************************************************************
*** calculate the first order internal radiation field of the m-th   ***
*** fourier coefficient of an isolated homogeneous layer             ***
*** see Hovenier (1971) and Wauben (1989)                            ***
*** the m-th fourier coefficient of the phasematrix of layer il is   ***
*** contained in the arrays Zmmin(j,i,il) and Zmplus(j,i,il)         ***
*** etmu contains exp[-taun(kn,il)/xmu(i)] with taun(kn,il) the      ***
*** optical depth of level kn in isolated layer il                   ***
*** ebtmu contains exp[-(taun(0,il)-taun(kn,il))/xmu(i)]             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu), nlirf(0:NDlay+1)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension etmu(NDmu,0:NDirf+1,0:NDlay)
      dimension ebtmu(NDmu,0:NDirf+1,NDlay)
      dimension taun(0:NDirf+1,0:NDlay)
      dimension r1(NDsup,NDsup,NDirf,NDlay)
      dimension t1(NDsup,NDsup,NDirf,NDlay)

      ag4 = .25d0
************************************************************************
* loop 10 (j,i,kn) is the loop over each matrix element (j,i) and      *
* for all the irf levels kn                                            *
* determine the reflection and transmission function yr and yt, see    *
* Wauben (1989) Eqs. (3)-(6) with mu0->i and mu->j                     *
************************************************************************
      do 10 kn=1,nlirf(il)
        do 10 j=1,nmutot
          jm = nmat*(j-1)
          xj = xmu(j)
          etj = etmu(j,kn,il)
          ebtj = ebtmu(j,kn,il)
          awj = ag4*wmu(j)
          do 10 i=1,nmutot
            im = nmat*(i-1)
            xi = xmu(i)
            eti = etmu(i,kn,il)
            ebti = ebtmu(i,kn,il)
            awji = awj*wmu(i)
            xsum = xi+xj
            if (dabs(xi-xj) .gt. 1.d-10) then
************************************************************************
* Eq. (3)                                                              *
************************************************************************
              yt = awji*(eti-etj)/(xi-xj)
            else
************************************************************************
* special case of the transmission matrix if mu(i)=mu(j)  Eq. (4)      *
************************************************************************
              yt = 0.d0
              if (xi .gt. 1.d-10) yt = awji/xi*taun(kn,il)*eti/xi
            endif
            if (xsum .gt. 1.d-10) then
************************************************************************
* Eq. (5)                                                              *
************************************************************************
              yr = awji*eti/xsum*(1.d0-ebti*ebtj)
            else
************************************************************************
* the reflection matrix is not defined for mu(i)=mu(j)=0, Eq.(6)       *
************************************************************************
              yr = 0.d0
            endif
            do 10 k=1,nmat
              ik = im+k
              do 10 l=1,nmat
                jl = jm+l
                r1(jl,ik,kn,il) = yr*Zmmin(jl,ik,il)
                t1(jl,ik,kn,il) = yt*Zmplus(jl,ik,il)
   10 continue
      return
      end


      subroutine doad_ord2m(nmat,nmutot,nmug,xmu,wmu,taun,etmu,m,Zmmin,
     +    Zmplus,ebtmu,nlirf,n,NDsup,NDmu,NDlay,NDirf,r2,t2)
************************************************************************
*** calculate the second order internal radiation field of the m-th  ***
*** fourier coefficient of an isolated homogeneous layer             ***
*** see Hovenier (1971) and Wauben (1989)                            ***
*** the m-th fourier coefficient of the phasematrix of layer n is    ***
*** contained in the arrays Zmmin(j,i,n) and Zmplus(j,i,n)           ***
*** etmu contains exp[-taun(kn,n)/xmu(i)] with taun(kn,n) the optical***
*** depth of level kn in isolated layer n                            ***
*** ebtmu contains exp[-(taun(0,n)-taun(kn,n))/xmu(i)]               ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu), nlirf(0:NDlay+1)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension etmu(NDmu,0:NDirf+1,0:NDlay)
      dimension ebtmu(NDmu,0:NDirf+1,NDlay)
      dimension taun(0:NDirf+1,0:NDlay)
      dimension r2(NDsup,NDsup,NDirf,NDlay)
      dimension t2(NDsup,NDsup,NDirf,NDlay)
      dimension yik(4,4)

      do 20 kl=1,nmat
      do 20 il=1,nmat
        yik(il,kl) = 1.d0
        if ((kl .ge. 3) .and. (il .le. 2)) yik(il,kl) = -1.d0
        if ((il .ge. 3) .and. (kl .le. 2)) yik(il,kl) = -1.d0
   20 continue
      ag16 = 1.d0/16.d0
************************************************************************
* loop 60 (j,i,kn,n) is the loop over each matrix element (j,i) and    *
* for all the irf levels kn and the layers n                           *
************************************************************************
      do 60 kn=1,nlirf(n)
        do 60 j=1,nmutot
          jm = nmat*(j-1)
          xj = xmu(j)
          etj = etmu(j,kn,n)
          ebtj = ebtmu(j,kn,n)
          awj = ag16*wmu(j)
          jz = 0
          if (xj .gt. 1.d-10) jz = 1
          tdj = 0.d0
          if (jz .eq. 1) tdj = taun(kn,n)/xj
          btdj = 0.d0
          if (jz .eq. 1) btdj = (taun(0,n)-taun(kn,n))/xj
          do 60 i=1,nmutot
            im = nmat*(i-1)
            xi = xmu(i)
            eti = etmu(i,kn,n)
            ebi = etmu(i,0,n)
            ebti = ebtmu(i,kn,n)
            awji = awj*wmu(i)
            iz = 0
            if (xi .gt. 1.d-10) iz = 1
************************************************************************
* the reflection and transmission matrices are not defined for         *
* mu(i)=mu(j)=0, r2 and t2 are set to zero for that case               *
************************************************************************
            if ((jz .eq. 0) .and. (iz .eq. 0)) then
              do 30 k=1,nmat
                ik = im+k
                do 30 l=1,nmat
                  jl = jm+l
                  r2(jl,ik,kn,n) = 0.d0
                  t2(jl,ik,kn,n) = 0.d0
   30         continue
            else
              tdi = 0.d0
              if (iz .eq. 1) tdi = taun(kn,n)/xi
              btdi = 0.d0
              if (iz .eq. 1) btdi = (taun(0,n)-taun(kn,n))/xi
              xipxj = xi+xj
              ximxj = xi-xj
              xiximj = 0.d0
              if (i .ne. j) xiximj = xi/ximxj
              etimj = eti-etj
              ebtij = 1.d0-ebti*ebtj
              tx = xiximj*etimj
              rx = xi/xipxj*eti*ebtij
************************************************************************
* gaussian integration over mu(k) in loop 50                           *
************************************************************************
              do 50 k=1,nmug
                km = nmat*(k-1)
                xk = xmu(k)
                etk = etmu(k,kn,n)
                ebk = etmu(k,0,n)
                ebtk = ebtmu(k,kn,n)
                xipxk = xi+xk
                xjpxk = xj+xk
                ximxk = xi-xk
                xjmxk = xj-xk
                etimk = eti-etk
                etjmk = etj-etk
                ebtik = 1.d0-ebti*ebtk
                ebtjk = 1.d0-ebtj*ebtk
************************************************************************
* determine the transmission and reflection functions tm, tp, rm, rp   *
* see Wauben (1989) Eqs. (7)-(17) with mu->j, mu0->i and mu'->k        *
************************************************************************
                if (i .eq. j) then
                  if (i .ne. k) then
                    tm = (tdi*eti+xk/xipxk*ebi*(ebk*eti-ebtk))/xipxk
                    tp = (tdi*eti-xk/ximxk*etimk)/ximxk
                    rm = (rx+xk/xjmxk*ebi*(ebtk-ebtj))/xipxk
                    rp = (rx-xk/xjpxk*etk*ebtjk)/ximxk
                  else
                    tm = (tdi*eti+xk/xipxk*ebi*(ebk*eti-ebtk))/xipxk
                    tp = .5d0*tdi*tdi*eti/xi
                    rm = (rx-btdj*ebi*ebtj)/xipxj
                    rp=((xj/xipxj+tdi)*ebtij-btdi*ebtj*ebti)*eti/xipxj
                  endif
                elseif (i .eq. k) then
                  tm = (tx+xk/xjpxk*ebi*(ebk*etj-ebtk))/xipxk
                  tp = (tdi*eti-xj/ximxj*etimj)/ximxj
                  rm = (rx+xk/xjmxk*ebi*(ebtk-ebtj))/xipxk
                  rp = ((xj/xipxj+tdi)*ebtij-btdi*ebtj*ebti)*eti/xipxj
                elseif (j .eq. k) then
                  tm = (tx+xk/xjpxk*ebi*(ebk*etj-ebtk))/xipxk
                  tp = (tx-tdj*etj)/ximxj
                  rm = (rx-btdj*ebi*ebtj)/xipxj
                  rp = (rx-xk/xjpxk*etk*ebtjk)/ximxk
                else
                  tm = (tx+xk/xjpxk*ebi*(ebk*etj-ebtk))/xipxk
                  tp = (tx-xk/xjmxk*etjmk)/ximxk
                  rm = (rx+xk/xjmxk*ebi*(ebtk-ebtj))/xipxk
                  rp = (rx-xk/xjpxk*etk*ebtjk)/ximxk
                endif
                y = wmu(k)*wmu(k)*awji/xk
************************************************************************
* the matrixproduct z^m_{j,k}*z^m_{k,i}, where zmsik and zpsik take    *
* account of changing the signs of mu(i) and mu(k) in z^m_{i,k}        *
* see Hovenier (1969) relation C which leads to a minus sign for the   *
* matrix elements that are uneven functions of (phi-phi')              *
************************************************************************
                do 40 il=1,nmat
                  iil = im+il
                  do 40 jl=1,nmat
                    jjl = jm+jl
                    do 40 kl=1,nmat
                      kkl = km+kl
                      zpjk = Zmplus(jjl,kkl,n)
                      zpki = Zmplus(kkl,iil,n)
                      zmjk = Zmmin(jjl,kkl,n)
                      zmki = Zmmin(kkl,iil,n)
                      zmsjk = yik(jl,kl)*zmjk
                      zpsjk = yik(jl,kl)*zpjk
       t2(jjl,iil,kn,n)=t2(jjl,iil,kn,n)+y*(zmsjk*zmki*tm+zpjk*zpki*tp)
       r2(jjl,iil,kn,n)=r2(jjl,iil,kn,n)+y*(zpsjk*zmki*rm+zmjk*zpki*rp)
   40           continue
   50         continue
            endif
   60 continue
      return
      end



      subroutine doad_prodes(w,e,s,nmutot,nmug,nmat,NDsup,NDmu)
************************************************************************
*** multiply supervector S(NDsup) with diagonal matrix E(NDmu)       ***
*** W = E*S this is a full matrix multiplication AP. page 379        ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension w(NDsup), s(NDsup), e(NDmu)

      do 10 i=1,nmutot
	im = (i-1)*nmat
        do 10 j=1,nmat
          k = im+j
          w(k) = s(k)*e(i)
   10 continue
      return
      end


      subroutine doad_prodex(w,x,e,nmutot,nmug,nmat,NDsup,NDmu,ind)
************************************************************************
*** multiply supermatrix X(NDsup*NDsup) with diagonal matrix E(NDmu) ***
*** this is a full matrix multiplication AP. page 379                ***
***    ind.eq.1: W = X*E                                             ***
***    ind.ne.1: W = E*X                                             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension w(NDsup,NDsup), x(NDsup,NDsup), e(NDmu)

      nsup = nmat*nmutot
      if (ind .eq. 1) then
          do 20 i=1,nmutot
          do 20 j=1,nmat
              k = (i-1)*nmat+j
              do 10 l=1,nsup
                  w(l,k) = x(l,k)*e(i)
   10         continue
   20     continue
      else
          do 30 l=1,nsup
          do 30 i=1,nmutot
          do 30 j=1,nmat
              k = (i-1)*nmat+j
              w(k,l) = x(k,l)*e(i)
   30     continue
      endif
      return
      end


      subroutine doad_prodmv(w,x,s,nmutot,nmug,nmat,NDsup)
************************************************************************
*** calculate the supermatrix/vector product W = X*S                 ***
*** this is a truncated matrix multiplication AP. Eq. (95)           ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension w(NDsup), x(NDsup,NDsup), s(NDsup)

      ng = nmug*nmat
      nsup = nmutot*nmat
************************************************************************
* compute W = X*S ,note that ng=nmug*nmat should be even               *
************************************************************************
      do 10 i=1,nsup
          w(i) = 0.d0
   10 continue
      do 30 k=1,ng
          do 20 i=1,nsup
              w(i) = w(i)+x(i,k)*s(k)
   20     continue
   30 continue
      return
      end


      subroutine doad_prodsm(w,x,y,nmutot,nmug,nmat,NDsup)
************************************************************************
*** calculate the supermatrix product W = X*Y                        ***
*** this is a truncated matrix multiplication AP. Eq. (95)           ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension w(NDsup,NDsup), x(NDsup,NDsup), y(NDsup,NDsup)

      ng = nmug*nmat
      nsup = nmutot*nmat
************************************************************************
* compute W = X*Y ,note that ng=nmug*nmat should be even               *
************************************************************************
      do 40 j=1,nsup
          do 10 i=1,nsup
              w(i,j) = 0.d0
   10     continue
          do 30 k=1,ng
              do 20 i=1,nsup
                  w(i,j) = w(i,j)+x(i,k)*y(k,j)
   20         continue
   30     continue
   40 continue
      return
      end

      subroutine doad_renorm(a,Zmmin,Zmplus,nmug,nmutot,nmat,xmu,wmu,il,
     +  NDlay,NDmu,NDsup,eps)
************************************************************************
*  Renormalize the phase matrix as indicated in de Haan et al. (1987)  *
*  section 7.5 to ensure that Eq. (144) is valid.                      *
*  The procedure is described in :                                     *
*  J.E. Hansen, 1971, J. Atmos. Sci. 28, page 1400, especially the     *
*  discussion on page 1422 and appendix B.                             *
*  We rewrite Hansen's iteration for f(i,j) as follows :               *
*                                                                      *
*            1                                                         *
*           f(i,j) = 1                      (i,j=1,...,nmutot)         *
*                                                                      *
*           r(j) = sum Zmmin(i,j)*w(i)      (sum over j=1,...,nmug)    *
*                                                                      *
*           t(j) = sum Zmplus(i,j)*w(i)     (sum over j=1,...,nmug)    *
*                                                                      *
*           eps(j) = |2 - r(j) - t(j)|                                 *
*                                                                      *
*            k        k-1           k-1        k-1                     *
*           f(i,j) = f(i,j) * [ratio(j) + ratio(i)]                    *
*                                                                      *
*            k          k                                              *
*           t(j) = sum f(i,j)*Zmplus(i,j)*w(i)                         *
*                                                                      *
*              k                  k                                    *
*           eps(j) = |2 - r(j) - t(j)|                                 *
*                                                                      *
*                                                                      *
*  where                                                               *
*                 k-1              k-1                                 *
*           ratio(i) = 0.5*[2-r(i)] / t(i)   (i=1,...,nmutot; k=2,...) *
*                                                                      *
*  These differ slightly from Hansen's appendix B because Hansen is    *
*  simply WRONG (he forgot the factor 1/2 in eq. (144) of de Haan et   *
*  al. (1987)) !!!                                                     *
************************************************************************
      implicit double precision (a-h,o-z)

      parameter ( maxit = 100 )

      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension xmu(NDmu), wmu(NDmu)
      dimension f(NDmu,NDmu), r(NDmu), t(NDmu), ratio(NDmu)
      dimension w(NDmu)

      if (nmutot .gt. NDmu) then
          print*,'dimension in subroutine RENORM too small'
          print*,'dim. is NDmu=',NDmu,' but nmutot=',nmutot
          stop ' dimension error in RENORM'
      endif
************************************************************************
*  initialize the normalization factors f on one                       *
************************************************************************
      do 10 j=1,nmutot
      do 10 i=1,nmutot
          f(i,j) = 1.d0
   10 continue
************************************************************************
*  retrieve the weights from the supermatrix factors wmu               *
************************************************************************
      do 20 i=1,nmug
          w(i) = 0.5d0*wmu(i)**2/xmu(i)
   20 continue
      if (nmutot .gt. nmug) then
          do 30 i=nmug+1,nmutot
              w(i) = 0.d0
   30     continue
      endif
************************************************************************
*  Calculate r(j), t(j) and epsj according to Appendix B of Hansen 1971*
*  for epsj we take the maximum of Hansen's eps(j). NOTE FACTOR 2 !!!  *
************************************************************************
      epsj = 0.d0
      do 50 j=1,nmutot
          r(j) = 0.d0
          t(j) = 0.d0
          jsup = (j-1)*nmat+1
          do 40 i=1,nmug
              isup = (i-1)*nmat+1
              r(j) = r(j)+Zmmin(isup,jsup,il)*w(i)
              t(j) = t(j)+Zmplus(isup,jsup,il)*w(i)
   40     continue
          diffj = dabs(2.d0*a-r(j)-t(j))
          if (epsj .lt. diffj) epsj=diffj
          ratio(j) = 0.5d0*(2.d0*a-r(j))/t(j)
   50 continue
************************************************************************
*  If epsj is too large : start the iteration as described in          *
*  appendix B of Hansen (1971)                                         *
************************************************************************
      iter = 1
   55 if ((epsj .gt. eps) .and. (iter .lt. maxit)) then
          iter = iter+1
          do 60 j=1,nmutot
          do 60 i=1,nmutot
              f(i,j) = f(i,j)*(ratio(j)+ratio(i))
   60     continue
          epsj = 0.d0
          do 80 j=1,nmutot
              t(j) = 0.d0
              jsup = (j-1)*nmat+1
              do 70 i=1,nmug
                  isup = (i-1)*nmat+1
                  t(j) = t(j)+f(i,j)*Zmplus(isup,jsup,il)*w(i)
   70         continue
              diffj = dabs(2.d0*a-r(j)-t(j))
              if (epsj .lt. diffj) epsj = diffj
              ratio(j) = 0.5d0*(2.d0*a-r(j))/t(j)
   80     continue
          goto 55
      endif
************************************************************************
*  Handle the case that iteration failed to converge                   *
************************************************************************
      if (iter .ge. maxit) then
          stop 'renorm: iteration for renormalization not converged'
      endif
************************************************************************
*  Perform the renormalization on the 'transmission' part only         *
************************************************************************
      do 100 k1=1,nmat
      do 100 k2=1,nmat
      do 100 j=1,nmutot
          jsup = (j-1)*nmat+k2
          do 90 i=1,nmutot
              isup = (i-1)*nmat+k1
              Zmplus(isup,jsup,il) = f(i,j)*Zmplus(isup,jsup,il)
   90     continue
  100 continue
      return
      end


      subroutine doad_sadd(nmat,nmutot,nmug,NDmu,NDsup,xmu,ebmu,eblow,
     +    r,t,radd,tadd,sr,st,sradd,stadd)
************************************************************************
*** adding-algorithm for 2 inhomogeneous layers                      ***
*** on entrance r and t contain the reflection and transmission of   ***
*** the upper layer with ebmu, and radd and tadd contain the         ***
*** reflection and transmission of the lower layer(s) with eblow     ***
*** on exit radd and tadd contain the reflection and transmission of ***
*** the added layer and r and t contain u and d at the interface     ***
*** see de Haan et al. (1987) Eqs. (85)-(91)                         ***
*** similar for the supervectors sr, st, sradd and stadd             ***
************************************************************************
      implicit double precision (a-h,o-z)


      dimension xmu(NDmu), ebmu(NDmu), eblow(NDmu)
      dimension r(NDsup,NDsup), t(NDsup,NDsup)
      dimension radd(NDsup,NDsup), tadd(NDsup,NDsup)
      dimension sr(NDsup), st(NDsup), sradd(NDsup), stadd(NDsup)
      dimension s(NDsup,NDsup), e(NDsup,NDsup), p(NDsup,NDsup)
      dimension u(NDsup,NDsup), d(NDsup,NDsup)
      dimension su(NDsup), sd(NDsup), se(NDsup), sp(NDsup)
cmv   common /mat/ s,e,p,u,d

      ng = nmug*nmat
      nsup = nmat*nmutot
      if (NDsup .ne. NDsup) then
          print*,'dimension in subroutine SADD is wrong'
          print*,'dim. is NDsup=',NDsup,' but NDsup=',NDsup
          stop ' dimension error in SADD'
      endif
      do 10 j=1,nsup
      do 10 i=1,nsup
          p(i,j) = r(i,j)
   10 continue
      call doad_transhm(p,nmutot,nmat,NDsup)
************************************************************************
* calculate Q_1, Eq. (85)                                              *
* and e and s are initialized to Q_1                                   *
************************************************************************
      call doad_prodsm(e,p,radd,nmutot,nmug,nmat,NDsup)
      do 20 j=1,nsup
      do 20 i=1,nsup
          s(i,j) = e(i,j)
   20 continue
      delta = 0.d0
      do 30 i=1,ng
          delta = delta+e(i,i)
   30 continue
      if ((delta*delta) .lt. 1.d-12) goto 100
      ik = 0
************************************************************************
* productmethod to calculate S, Eqs. (112)-(115)                       *
* e = C_r*C_r                                                          *
************************************************************************
   40 call doad_prodsm(p,e,e,nmutot,nmug,nmat,NDsup)
      do 50 j=1,nsup
      do 50 i=1,nsup
          e(i,j) = p(i,j)
   50 continue 
      call doad_prodsm(p,s,e,nmutot,nmug,nmat,NDsup)
************************************************************************
* P = S_r*C_r+1                                                        *
************************************************************************
      do 60 j=1,nsup
      do 60 i=1,nsup
          s(i,j) = p(i,j)+e(i,j)+s(i,j)
   60 continue
************************************************************************
* S = S_r+1                                                            *
************************************************************************
      delta = 0.d0
      do 70 i=1,ng
          delta = delta+e(i,i)
   70 continue
************************************************************************
* delta=trace(C_r+1), Eq. (124), is an estimate of the error           *
************************************************************************
      ik = ik+1
      if (ik .gt. 20) then
          write(*,90)
          write(*,80)
          write(*,90)
   80     format('*** stop : S-series converge too slow ***')
   90     format('****************************************')
          stop 'slow convergence in ADDTOT'
      endif
      if ((delta*delta) .gt. 1.d-12) goto 40
************************************************************************
* calculate D, Eq. (88)                                                *
* the repeated reflections matrix is contained in s                    *
************************************************************************
  100 call doad_prodsm(d,s,t,nmutot,nmug,nmat,NDsup)
      do 110 i=1,nmutot
          im = nmat*(i-1)
          do 110 j=1,nmutot
              jm = nmat*(j-1)
              ej = ebmu(j)
              do 110 k=1,nmat
                  ik = im+k
                  do 110 l=1,nmat
                      jl = jm+l
                      d(ik,jl) = d(ik,jl)+ej*s(ik,jl)+t(ik,jl)
  110 continue
************************************************************************
* calculate U, Eq. (89)                                                *
************************************************************************
      call doad_prodsm(u,radd,d,nmutot,nmug,nmat,NDsup)
      do 120 i=1,nmutot
          im = nmat*(i-1)
          do 120 j=1,nmutot
              jm = nmat*(j-1)
              ej = ebmu(j)
              do 120 k=1,nmat
                  ik = im+k
                  do 120 l=1,nmat
                      jl = jm+l
                      u(ik,jl) = u(ik,jl)+ej*radd(ik,jl)
  120 continue
************************************************************************
* calculate R, Eq. (90)                                                *
************************************************************************
      call doad_transhm(t,nmutot,nmat,NDsup)
      call doad_prodsm(e,t,u,nmutot,nmug,nmat,NDsup)
      do 130 i=1,nmutot
          im = nmat*(i-1)
          ei = ebmu(i)
          do 130 j=1,nmutot
              jm = nmat*(j-1)
              do 130 k=1,nmat
                  ik = im+k
                  do 130 l=1,nmat
                      jl = jm+l
                      e(ik,jl) = e(ik,jl)+ei*u(ik,jl)+r(ik,jl)
  130 continue
************************************************************************
* calculate T, Eq. (91)                                                *
************************************************************************
      call doad_prodsm(p,tadd,d,nmutot,nmug,nmat,NDsup)
      do 140 i=1,nmutot
          im = nmat*(i-1)
          ei = eblow(i)
          do 140 j=1,nmutot
              jm = nmat*(j-1)
              ej = ebmu(j)
              do 140 k=1,nmat
                  ik = im+k
                  do 140 l=1,nmat
                      jl = jm+l
                      p(ik,jl) = ej*tadd(ik,jl)+ei*d(ik,jl)+e(ik,jl)
  140 continue
************************************************************************
* calculate the radiation due to the sources see Wauben (1990)         *
* calculate sd according to Eq.(36)                                    *
************************************************************************
      call doad_prodmv(se,r,sradd,nmutot,nmug,nmat,NDsup)
      call doad_addsv(su,st,se,nsup,NDsup)
      call doad_prodmv(se,s,su,nmutot,nmug,nmat,NDsup)
      call doad_addsv(sd,se,su,nsup,NDsup)
************************************************************************
* calculate su according to Eq.(37)                                    *
************************************************************************
      call doad_prodmv(se,radd,sd,nmutot,nmug,nmat,NDsup)
      call doad_addsv(su,sradd,se,nsup,NDsup)
************************************************************************
* calculate sr according to Eq.(38)                                    *
************************************************************************
      call doad_prodmv(se,t,su,nmutot,nmug,nmat,NDsup)
      call doad_addsv(sp,se,sr,nsup,NDsup)
      call doad_prodes(se,ebmu,su,nmutot,nmug,nmat,NDsup,NDmu)
      call doad_addsv(sradd,sp,se,nsup,NDsup)
************************************************************************
* calculate st according to Eq.(39)                                    *
************************************************************************
      call doad_prodmv(se,tadd,sd,nmutot,nmug,nmat,NDsup)
      call doad_addsv(sp,se,stadd,nsup,NDsup)
      call doad_prodes(se,ebmu,sd,nmutot,nmug,nmat,NDsup,NDmu)
      call doad_addsv(stadd,sp,se,nsup,NDsup)
************************************************************************
* calculate sr and su from symmetry relation                           *
************************************************************************
      call doad_xeqy(r,p,NDsup)
      call doad_xeqy(t,e,NDsup)
      do 150 j=1,nsup
        sr(j) = su(j)
        st(j) = sd(j)
        do 150 i=1,nsup
          t(i,j) = d(i,j)
          r(i,j) = u(i,j)
          tadd(i,j) = p(i,j)
          radd(i,j) = e(i,j)
  150 continue
      return
      end


      subroutine doad_scalZm(m,il,NDlay,alfbet,NDcoef,ncoef,
     +    xmu,wmu,NDmu,
     +    nmug,nmutot,nmat,Zmmin,Zmplus,NDsup,eps)
************************************************************************
*  Calculate the m-th Fourier component of the phase matrix Zm(mu0,mu) *
*  from the expansion coefficients of the scattering matrix into       *
*  generalized spherical functions. The formulae can be found in :     *
*                                                                      *
*    J.F. de Haan et al.: 1987, Astron. Astrophys. 183, pp. 371-391.   *
*                                                                      *
*  Essentially, Eqs. (66)-(82) are used.                               *
*                                                                      *
*  no polarization !                                                   *
************************************************************************
      implicit double precision (a-h,o-z)

      dimension alfbet(4,4,0:NDcoef,NDlay), ncoef(NDlay)
      dimension wmu(NDmu), xmu(NDmu)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension Plm(NDmu,3,2), sqlm(0:NDcoef)

      if (ncoef(il) .gt. NDcoef) then
          print*,'dimension in subroutine SCALZM is too small'
          print*,'dim. is NDcoef=',NDcoef,' but ncoef=',ncoef(il)
          stop 'dimension error in SCALZM'
      endif
      if (nmutot .gt. NDmu) then
          print*,'dimension in subroutine SCALZM is too small'
          print*,'dim. is NDmu=',NDmu,' but nmutot=',nmutot
          stop 'dimension error in SCALZM'
      endif
************************************************************************
* Precompute the factor sqrt(l**2-m**2) needed in Eqs. (81)            *
************************************************************************
      m2 = m*m
      do 10 l=m,ncoef(il)
          sqlm(l) = dsqrt(dble(l*l-m2))
   10 continue
************************************************************************
* Compute the binomial factor 2m above m needed in Eq. (77)            *
************************************************************************
      bin2mm = 1.d0
      do 20 n=1,m
          bin2mm = bin2mm*dble(n+m)/dble(n)
   20 continue
      binfac = 2.d0**(-m)*dsqrt(bin2mm)
************************************************************************
* Initialize Plm for l=m, Eq.(77) without factor i**-m                 *
* this factor is cancelled in Eq. (66) by the factor (-1)**m           *
************************************************************************
      lold = 1
      lnew = 2
      do 30 i=1,nmutot
          u  = xmu(i)
          rootu = dsqrt(1.d0-u*u)
          Plm(i,1,lold) = 0.d0
          if (m .ne. 0) then
              Plm(i,1,lnew) = binfac*rootu**m
          else
              Plm(i,1,lnew) = binfac
          endif
   30 continue
************************************************************************
* Initialize phase matrix with term l=m, Eq.(66)                       *
************************************************************************
      do 50 i=1,nmutot
          SP = alfbet(1,1,m,il)*Plm(i,1,lnew)
          do 40 j=1,nmutot
              Zmplus(i,j,il) = SP*Plm(j,1,lnew)
              Zmmin(i,j,il) = Zmplus(i,j,il)
   40     continue
   50 continue
************************************************************************
* Start loop over l loop 100 (summation index in Eq. (66))             *
* The parity of Plm is (-1)**(l-m)                                     *
************************************************************************
      parity = 1.d0
      do 100 l=m+1,ncoef(il)
          parity = -parity
************************************************************************
* Do one step in recurrence for Plm, Eq.(81)                           *
************************************************************************
          c1 = dble(l+l-1)/sqlm(l)
          c2 = sqlm(l-1)/sqlm(l)
          do 60 i=1,nmutot
              u = xmu(i)
              Plm(i,1,lold) = c1*u*Plm(i,1,lnew)-c2*Plm(i,1,lold)
   60     continue
          ltmp = lnew
          lnew = lold
          lold = ltmp
************************************************************************
* Add a new term to Zm, Eq.(66)                                        *
************************************************************************
          do 80 i=1,nmutot
            SP = alfbet(1,1,l,il)*Plm(i,1,lnew)
            do 70 j=1,nmutot
             Zmplus(i,j,il) = Zmplus(i,j,il)+SP*Plm(j,1,lnew)
             Zmmin(i,j,il) = Zmmin(i,j,il)+parity*SP*Plm(j,1,lnew)
   70       continue
   80     continue
  100 continue
************************************************************************
* End of summation loop over l                                         *
************************************************************************
      if (m .eq. 0) then
************************************************************************
* renomalization of phasematrix for m=0                                *
************************************************************************
        a = alfbet(1,1,0,il)
        call doad_renorm(a,Zmmin,Zmplus,nmug,nmutot,nmat,xmu,wmu,il,
     +    NDlay,NDmu,NDsup,eps)
      endif
      return
      end


      subroutine doad_sdoubl(nmat,nmutot,nmug,NDmu,NDsup,ebmu,r,t,u,d,
     +  sr,st,su,sd)
************************************************************************
*** adding 2 identical homogeneous layers                            ***
*** on entrance r and t contain the reflection and transmission      ***
*** matrices of the identical layers                                 ***
*** sr and st contain the emerging radiation due to the sources      ***
*** on exit r and t contain the reflection and transmission matrices ***
*** of the combined layer and u and d contain the up- and downward   ***
*** radiation at the interface (middle) of the two layers            ***
*** see de Haan et al. (1987) Eqs. (85)-(91)                         ***
*** similarly sr,st,su and sd contain the radiation due to sources   ***
************************************************************************
      implicit double precision (a-h,o-z)

      dimension ebmu(NDmu)
      dimension r(NDsup,NDsup), t(NDsup,NDsup)
      dimension u(NDsup,NDsup), d(NDsup,NDsup)
      dimension sr(NDsup), st(NDsup), su(NDsup), sd(NDsup)
      dimension s(NDsup,NDsup), e(NDsup,NDsup), p(NDsup,NDsup)
      dimension se(NDsup)
cmv   common /mat/ s,e,p

      nsup = nmat*nmutot
      if (NDsup .ne. NDsup) then
          print*,'dimension in subroutine SDOUBL too small'
          print*,'dim. is NDsup=',NDsup,' but NDsup=',NDsup
          stop ' dimension error in SDOUBL'
      endif
************************************************************************
* calculate Q_1, Eq. (85)                                              *
************************************************************************
      do 10 j=1,nsup
      do 10 i=1,nsup
          p(i,j) = r(i,j)
   10 continue
      call doad_transhm(p,nmutot,nmat,NDsup)
      call doad_prodsm(e,p,r,nmutot,nmug,nmat,NDsup)
      do 20 j=1,nsup
      do 20 i=1,nsup
          s(i,j) = e(i,j)
   20 continue
      delta = 0.d0
      do 30 i=1,nmug
          delta = delta+e(i,i)
   30 continue
      if ((delta*delta) .lt. 1.d-12) goto 100
      ik = 0
************************************************************************
* productmethod for calculating S, Eqs. (112)-(115)                    *
************************************************************************
   40 call doad_prodsm(p,e,e,nmutot,nmug,nmat,NDsup)
      do 50 j=1,nsup
      do 50 i=1,nsup
          e(i,j) = p(i,j)
   50 continue
      call doad_prodsm(p,s,e,nmutot,nmug,nmat,NDsup)
************************************************************************
* P = S_r*C_r+1                                                        *
************************************************************************
      do 60 j=1,nsup
      do 60 i=1,nsup
          s(i,j) = p(i,j)+e(i,j)+s(i,j)
   60 continue
************************************************************************
* S = S_r+1                                                            *
************************************************************************
      delta = 0.d0
      do 70 i=1,nmug
          delta = delta+e(i,i)
   70 continue
************************************************************************
* delta=trace(c_r+1), Eq. (124), is an estimate of the error           *
************************************************************************
      ik = ik+1
      if (ik .gt. 20) then
          write(*,90)
          write(*,80)
          write(*,90)
   80     format('*** stop : S-series converge to slow ***')
   90     format('****************************************')
          stop 'slow convergence in DOUBLE'
      endif
      if ((delta*delta) .gt. 1.d-12) goto 40
************************************************************************
* calculate D, Eq. (88)                                                *
************************************************************************
  100 call doad_prodsm(p,s,t,nmutot,nmug,nmat,NDsup)
      do 110 i=1,nmutot
          im=nmat*(i-1)
          do 110 j=1,nmutot
              jm=nmat*(j-1)
              ej=ebmu(j)
              do 110 k=1,nmat
                  ik=im+k
                  do 110 l=1,nmat
                      jl=jm+l
                      d(ik,jl)=p(ik,jl)+ej*s(ik,jl)+t(ik,jl)
  110 continue
************************************************************************
* calculate U, Eq. (89)                                                *
************************************************************************
      call doad_prodsm(p,r,d,nmutot,nmug,nmat,NDsup)
      do 120 j=1,nmutot
          jm=nmat*(j-1)
          ej=ebmu(j)
          do 120 i=1,nmutot
              im=nmat*(i-1)
              do 120 k=1,nmat
                  ik=im+k
                  do 120 l=1,nmat
                      jl=jm+l
                      u(ik,jl)=p(ik,jl)+ej*r(ik,jl)
  120 continue
************************************************************************
* calculate T, Eq. (91)                                                *
************************************************************************
      call doad_prodsm(e,t,d,nmutot,nmug,nmat,NDsup)
      do 130 i=1,nmutot
          im=nmat*(i-1)
          ei=ebmu(i)
          do 130 j=1,i
              jm=nmat*(j-1)
              ej=ebmu(j)
              do 130 k=1,nmat
                  ik=im+k
                  do 130 l=1,nmat
                      jl=jm+l
                      e(ik,jl)=e(ik,jl)+ei*d(ik,jl)+ej*t(ik,jl)
  130 continue
************************************************************************
* calculate R, Eq. (90)                                                *
************************************************************************
      call doad_transhm(t,nmutot,nmat,NDsup)
      call doad_prodsm(p,t,u,nmutot,nmug,nmat,NDsup)
      do 140 i=1,nmutot
          im=nmat*(i-1)
          ei=ebmu(i)
          do 140 j=1,i
              jm=nmat*(j-1)
              do 140 k=1,nmat
                  ik=im+k
                  do 140 l=1,nmat
                      jl=jm+l
                      p(ik,jl)=p(ik,jl)+ei*u(ik,jl)+r(ik,jl)
  140 continue
      call doad_symmut(p,nmat,nmutot,NDsup,3)
      call doad_symmut(e,nmat,nmutot,NDsup,4)
************************************************************************
* calculate the radiation due to the sources see Wauben (1990)         *
* calculate sd according to Eq.(36)                                    *
************************************************************************
      call doad_prodmv(se,r,sr,nmutot,nmug,nmat,NDsup)
      call doad_addsv(su,st,se,nsup,NDsup)
      call doad_prodmv(se,s,su,nmutot,nmug,nmat,NDsup)
      call doad_addsv(sd,se,su,nsup,NDsup)
************************************************************************
* calculate st according to Eq.(39)                                    *
************************************************************************
      call doad_prodmv(se,t,sd,nmutot,nmug,nmat,NDsup)
      call doad_addsv(su,se,st,nsup,NDsup)
      call doad_prodes(se,ebmu,sd,nmutot,nmug,nmat,NDsup,NDmu)
      call doad_addsv(st,su,se,nsup,NDsup)
************************************************************************
* calculate sr and su from symmetry relation                           *
************************************************************************
      do 160 i=1,nsup
        sr(i) = st(i)
        su(i) = sd(i)
  160 continue
      call doad_xeqy(r,p,NDsup)
      call doad_xeqy(t,e,NDsup)
      return
      end


      subroutine doad_seconm(nmat,nmutot,nmug,xmu,wmu,b,
     +    ebmu,m,n,r,t,Zmmin,Zmplus,NDmu,NDsup,NDlay)
************************************************************************
*** calculate the second order m-th fourier coefficient of the       ***
*** reflection and transmission matrices and add them to r and t in  ***
*** supermatrix form see Hovenier (1971)                             ***
*** the m-th fourier coefficient of the phasematrix of the layer is  ***
*** contained in the arrays Zmmin(i,j,n) and Zmplus(i,j,n)           ***
*** ebmu contains exp[-b/mu(i)] with b the optical thickness of the  ***
*** layer (starting, or total layer in case of second scat. only)    ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu), ebmu(NDmu)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension r(NDsup,NDsup), t(NDsup,NDsup)
      dimension yik(4,4)

      aag16 = 1.d0/16.d0
      do 10 kl=1,nmat
      do 10 il=1,nmat
        yik(il,kl) = 1.d0
        if ((kl .ge. 3) .and. (il .le. 2)) yik(il,kl) = -1.d0
        if ((il .ge. 3) .and. (kl .le. 2)) yik(il,kl) = -1.d0
   10 continue
************************************************************************
* loop 80 (i,j) is the loop over each matrix element contained in      *
* the lower triangle of the supermatrix                                *
************************************************************************
      do 80 i=1,nmutot
        xi = xmu(i)
        ei = ebmu(i)
        awi = aag16*wmu(i)
        im = (i-1)*nmat
        iz = 0
        if (xi .gt. 1.d-10) iz = 1
        bgi = 0.d0
        if (iz .eq. 1) bgi = b/xi
        beigxi = bgi*ei
        do 80 j=1,i
          xj = xmu(j)
          ej = ebmu(j)
          awij = awi*wmu(j)
          jm = (j-1)*nmat
          jz = 0
          if (xj .gt. 1.d-10) jz = 1
************************************************************************
* the reflection and transmission matrices are not defined for         *
* mu(i)=mu(j)=0, r and t are set to zero for that nmat*nmat submatrix  *
************************************************************************
          if ((jz .eq. 0) .and. (iz .eq. 0)) then
            do 20 l=1,nmat
              jl = jm+l
              do 20 k=1,nmat
                ik = im+k
                r(ik,jl) = 0.d0
                t(ik,jl) = 0.d0
   20       continue
          else
            bgj = 0.d0
            if (jz .eq. 1) bgj = b/xj
            bij = bgi+bgj
            xipxj = xi+xj
            ximxj = xi-xj
            xjximj = 0.d0
            if (i .ne. j) xjximj = xj/ximxj
            eimej = ei-ej
            eiej = 1.d0-ei*ej
************************************************************************
* taylor expansion of exponent to avoid loss of accuracy if b << 1     *
************************************************************************
            if (eiej .lt. 1.d-3) then
              eiej = bij*(1.d0-bij*.5d0*(1.d0-bij/3.d0))
              eimej = bgj*(1.d0-bgj*.5d0*(1.d0-bgj/3.d0))-
     +          bgi*(1.d0-bgi*.5d0*(1.d0-bgi/3.d0))
            endif
            e1 = xjximj*eimej
            g1 = xj/xipxj*eiej
************************************************************************
* gaussian integration over mu(k) in loop 60                           *
************************************************************************
            do 60 k=1,nmug
              xk = xmu(k)
              ek = ebmu(k)
              km = (k-1)*nmat
              bgk = b/xk
              bik = bgi+bgk
              bjk = bgj+bgk
************************************************************************
* b << 1, with special cases mu(i)=0 and mu(j)=0                       *
************************************************************************
              if ((bgi .lt. 1.d-3) .and. (bgj .lt. 1.d-3) .and.
     +          (bgk .lt. 1d-3)) then
                if (iz .eq. 0) then
                  z = b/xk/xj
                  e = 0.d0
                  g = z*(1.d0-bjk*.5d0*(1.d0-bjk/3.d0))
                  f = g-z*z*b/6.d0
                  h = 0.d0
                elseif (jz .eq. 0) then
                  z = b/xk/xi
                  e = 0.d0
                  h = z*(1.d0-bik*.5d0*(1.d0-bik/3.d0))
                  f = h-z*z*b/6.d0
                  g = 0.d0
                else
                  z = bgi*bgj*bgk/b*.5d0
                  e = z*(1.d0-(bgk+2.d0*(bgi+bgj))/3.d0)
                  f = z*(1.d0-(bgk+bgi+bgj)/3.d0)
                  g = z*(1.d0-(bgk+bgi+2.d0*bgj)/3.d0)
                  h = z*(1.d0-(bgk+bgj+2.d0*bgi)/3.d0)
                endif
              else   
                xipxk = xi+xk
                xjpxk = xj+xk
                ximxk = xi-xk
                xjmxk = xj-xk
                eimek = ei-ek
                ejmek = ej-ek
                eiek = 1.d0-ei*ek
                ejek = 1.d0-ej*ek
************************************************************************
* taylor expansion of exponent to avoid loss of accuracy if b << 1     *
************************************************************************
                if (eiek .lt. 1.d-3) then
                  eiek = bik*(1.d0-bik*.5d0*(1.d0-bik/3.d0))
                  eimek = bgk*(1.d0-bgk*.5d0*(1.d0-bgk/3.d0))-
     +              bgi*(1.d0-bgi*.5d0*(1.d0-bgi/3.d0))
                endif
                if (ejek .lt. 1.d-3) then
                  ejek = bjk*(1.d0-bjk*.5d0*(1.d0-bjk/3.d0))
                  ejmek = bgk*(1.d0-bgk*.5d0*(1.d0-bgk/3.d0))-
     +              bgj*(1.d0-bgj*.5d0*(1.d0-bgj/3.d0))
                endif
************************************************************************
* determine the functions e,f,g and h for different cases              *
* see Hovenier (1971) Eqs. (A-28)-(A-38) with mu->i, mu0->j and mu'->k *
************************************************************************
                if (i .eq. j) then
                  if (i .ne. k) then
                    e = (b/xj*ej-xk/xipxk*ej*ejek)/xjpxk
                    f = (b/xj*ej-xk/xjmxk*ejmek)/xjmxk
                    g = (g1-xk/ximxk*ej*eimek)/xjpxk
                    h = (g1-xk/xipxk*eiek)/xjmxk
                  else
                    e = (b/xj*ej-xk/xipxk*ej*ejek)/xjpxk
                    f = b*b*.5d0/xj/xj/xj*ej
                    g = (g1-beigxi*ej)/xipxj
                    h = g
                  endif
                elseif (i .eq. k) then
                  e = (e1-xk/xipxk*ej*eiek)/xjpxk
                  f = (e1-beigxi)/xjmxk
                  g = (g1-beigxi*ej)/xipxj
                  h = (g1-xk/xipxk*eiek)/xjmxk
                elseif (j .eq. k) then
                  e = (e1-xk/xipxk*ej*eiek)/xjpxk
                  f = (xi/xj*e1-bgj*ej)/ximxk
                  g = (g1-xk/ximxk*ej*eimek)/xjpxk
                  h = (xi/xipxj*eiej-bgj*ei*ej)/xipxj
                else
                  e = (e1-xk/xipxk*ej*eiek)/xjpxk
                  f = (e1-xk/ximxk*eimek)/xjmxk
                  g = (g1-xk/ximxk*ej*eimek)/xjpxk
                  h = (g1-xk/xipxk*eiek)/xjmxk
                endif
              endif
              y = wmu(k)*wmu(k)*awij/xk
              do 50 jl=1,nmat
                jjl = jm+jl
                do 50 il=1,nmat
                  iil = im+il
************************************************************************
* the matrixproduct z^m_{i,k}*z^m_{k,j}, where zmsik and zpsik take    *
* account of changing the signs of mu(i) and mu(k) in z^m_{i,k}        *
* see Hovenier (1969) relation C which leads to a minus sign for the   *
* matrix elements that are uneven functions of (phi-phi')              *
************************************************************************
                  do 40 kl=1,nmat
                    kkl = km+kl
                    zpnik = Zmplus(iil,kkl,n)
                    zpnkj = Zmplus(kkl,jjl,n)
                    zmnik = Zmmin(iil,kkl,n)
                    zmnkj = Zmmin(kkl,jjl,n)
                    zmsik = yik(il,kl)*zmnik
                    zpsik = yik(il,kl)*zpnik
                   r(iil,jjl)=r(iil,jjl)+y*(zpsik*zmnkj*g+zmnik*zpnkj*h)
                   t(iil,jjl)=t(iil,jjl)+y*(zmsik*zmnkj*e+zpnik*zpnkj*f)
   40             continue
   50         continue
   60       continue
          endif
   80 continue
      call doad_symmut(r,nmat,nmutot,NDsup,3)
      call doad_symmut(t,nmat,nmutot,NDsup,4)
      return
      end


      subroutine doad_sflux(nmat,nmutot,nmug,xmu,wmu,m,nlayer,nlirf,
     +    sfluxu,sfluxd,
     +    sukk,sdkk,sunn,sdnn,tau,NDmu,NDsup,NDirf,NDlay)
************************************************************************
*** calculates the upward (fluxu) and the downward (fluxd) fluxes    ***
*** at all the interface layers and irf levels                       ***
*** see Stammes et al. 1989, here F(tau) is calculated               ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu)
      dimension sukk(NDsup,NDirf,NDlay), sunn(NDsup,NDlay+1)
      dimension sdkk(NDsup,NDirf,NDlay), sdnn(NDsup,NDlay+1)
      dimension tau(0:NDirf+1,0:NDlay+1)
      dimension nlirf(0:NDlay+1)
      dimension sfluxu(0:NDirf,NDlay+1)
      dimension sfluxd(0:NDirf,NDlay+1)

      if (m .ne. 0) return
      do 5 il=1,nlayer+1
      do 5 k=0,nlirf(il)
        sfluxu(k,il) = 0.d0
        sfluxd(k,il) = 0.d0
    5 continue
      do 60 il=1,nlayer+1
        do 10 i=1,nmug
          i1 = nmat*(i-1)+1
          wi = wmu(i)
          sfluxu(0,il) = sfluxu(0,il)+wi*sunn(i1,il)
          sfluxd(0,il) = sfluxd(0,il)+wi*sdnn(i1,il)
   10   continue
        if (nlirf(il) .ne. 0) then
          do 50 k=1,nlirf(il)
            do 30 i=1,nmug
              i1 = nmat*(i-1)+1
              wi = wmu(i)
              sfluxu(k,il) = sfluxu(k,il)+wi*sukk(i1,k,il)
              sfluxd(k,il) = sfluxd(k,il)+wi*sdkk(i1,k,il)
   30       continue
   50     continue
        endif
   60 continue
      return
      end


      subroutine doad_shomin(nlirf,nl,nmat,nmutot,nmug,ukk,dkk,ek,
     +    sukk,sdkk,NDsup,NDirf,NDlay,NDmu)
************************************************************************
*** internal radiation field of the isolated homogeneous layer nl    ***
*** the method is described in:                                      ***
*** Stammes, P., de Haan, J.F., Hovenier, J.W.:                      ***
*** 1989, Astron. Astrophys. xxx, p. xxx                             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension ukk(NDsup,NDsup,NDirf,NDlay)
      dimension dkk(NDsup,NDsup,NDirf,NDlay)
      dimension sukk(NDsup,NDirf,NDlay), sdkk(NDsup,NDirf,NDlay)
      dimension ek(NDmu,NDirf,NDlay), nlirf(0:NDlay+1)
************************************************************************
* on entrance the matrix ukk(i,j,k,n) contains the supermatrix U, at   *
* middle (level k) of a partial layer n build in k doubling steps      *
* on exit ukk(i,j,k,n) contains the supermatrix U, at level k in the   *
* entire isolated homogeneous layer n                                  *
* ek(i,k,n) contains the factor exp[-b(n)/u(i)*2**(k-M-1)], where      *
* b(n) is the optical thickness of layer n                             *
************************************************************************

      dimension uk(NDsup,NDsup), dk(NDsup,NDsup)
      dimension ukstar(NDsup,NDsup), dkstar(NDsup,NDsup)
      dimension ukp1(NDsup,NDsup), dkp1(NDsup,NDsup)
      dimension ukm1(NDsup,NDsup)
      dimension ekp1(NDmu), em(NDmu), e(NDmu)
      dimension x(NDsup,NDsup), y(NDsup,NDsup), z(NDsup,NDsup)
      dimension suk(NDsup), sdk(NDsup), sdkp1(NDsup)
      dimension sx(NDsup), sy(NDsup), sz(NDsup)
cmv   common /mat/ x,y,z,ukp1,dkp1,uk,dk,ukstar,dkstar,ukm1

      if (nlirf(nl) .eq. 0) return
      if (nmutot .gt. NDmu) then
          print*,'dimension in subroutine SHOMIN too small'
          print*,'dim. is NDmu=',NDmu,' but nmutot=',nmutot
          stop ' dimension error in SHOMIN'
      endif
      nsup = nmat*nmutot
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine SHOMIN too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop ' dimension error in SHOMIN'
      endif
      nl2 = (nlirf(nl)+1)/2
      do 80 k=nl2-1,1,-1
          kup = nlirf(nl)+1-k
************************************************************************
* initializing arrays                                                  *
************************************************************************
          do 10 i=1,nmutot
              em(i) = ek(i,k,nl)
              ekp1(i) = ek(i,k+1,nl)
              e(i) = ekp1(i)**(2.d0**(nl2-k)-1.d0)
   10     continue
          do 20 i=1,nsup
            sdkp1(i) = sdkk(i,k+1,nl)
            suk(i) = sukk(i,k,nl)
            sdk(i) = sdkk(i,k,nl)
            do 20 j=1,nsup
              ukp1(j,i) = ukk(j,i,k+1,nl)
              dkp1(j,i) = dkk(j,i,k+1,nl)
              uk(j,i) = ukk(j,i,k,nl)
              dk(j,i) = dkk(j,i,k,nl)
   20     continue
************************************************************************
* internal radiation field of lower half of isolated layer             *
* D field, Eq. (64)                                                    *
************************************************************************
          call doad_prodex(x,dkp1,em,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_prodsm(y,dk,dkp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(z,x,y,nsup,NDsup)
          call doad_prodex(x,dk,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_addsm(y,x,z,nsup,NDsup)
          do 30 i=1,nsup
          do 30 j=1,nsup
              dkk(j,i,k,nl) = y(j,i)
   30     continue
************************************************************************
* internal radiation field of lower half of isolated layer             *
* D field for sources, Eq. (52)                                        *
************************************************************************
          call doad_prodes(sx,em,sdkp1,nmutot,nmug,nmat,NDsup,NDmu)
          call doad_prodmv(sy,dk,sdkp1,nmutot,nmug,nmat,NDsup)
          call doad_addsv(sz,sx,sy,nsup,NDsup)
          call doad_addsv(sy,sdk,sz,nsup,NDsup)
          do 35 j=1,nsup
              sdkk(j,k,nl) = sy(j)
              sukk(j,kup,nl) = sy(j)
   35     continue
************************************************************************
* U field, Eq. (65)                                                    *
************************************************************************
          call doad_prodex(x,uk,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_prodsm(y,uk,dkp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(z,x,y,nsup,NDsup)
          do 40 i=1,nsup
          do 40 j=1,nsup
              ukk(j,i,k,nl) = z(j,i)
   40     continue
************************************************************************
* U field for sources Eq. (53)                                         *
************************************************************************
          call doad_prodmv(sy,uk,sdkp1,nmutot,nmug,nmat,NDsup)
          call doad_addsv(sz,suk,sy,nsup,NDsup)
          do 45 j=1,nsup
              sukk(j,k,nl) = sz(j)
              sdkk(j,kup,nl) = sz(j)
   45     continue
************************************************************************
* internal radiation field of upper half of isolated layer             *
************************************************************************
          do 50 i=1,nsup
          do 50 j=1,nsup
              ukm1(j,i) = ukk(j,i,kup-1,nl)
   50     continue
************************************************************************
* calculate U*kk and D*kk according to Eqs. (68) and (69)              *
************************************************************************
          call doad_xeqy(ukstar,uk,NDsup)
          call doad_xeqy(dkstar,dk,NDsup)
          call doad_transhm(ukstar,nmutot,nmat,NDsup)
          call doad_transhm(dkstar,nmutot,nmat,NDsup)
************************************************************************
* D field, Eq. (70)                                                    *
************************************************************************
          call doad_prodsm(z,ukstar,ukm1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(x,z,dk,nsup,NDsup)
          do 60 i=1,nsup
          do 60 j=1,nsup
              dkk(j,i,kup,nl) = x(j,i)
   60     continue
************************************************************************
* U field, Eq. (71)                                                    *
************************************************************************
          call doad_prodex(z,ukm1,em,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_prodsm(x,dkstar,ukm1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(y,x,z,nsup,NDsup)
          call doad_addsm(z,y,uk,nsup,NDsup)
          do 70 i=1,nsup
          do 70 j=1,nsup
              ukk(j,i,kup,nl) = z(j,i)
   70     continue
   80 continue
      return
      end
 

      subroutine doad_shomly(nmat,nmutot,nmug,xmu,wmu,b,m,il,ebmu,
     +    r,t,nlirf,ndoubl,ukk,dkk,Zmmin,Zmplus,
     +    NDsup,NDirf,NDmu,NDlay,ssor,sr,st,sukk,sdkk)
************************************************************************
*** calculate reflection and transmission r and t of an isolated     ***
*** homogeneous layer with internal sources                          ***
*** including all orders of scatt.                                   ***
*** see de Haan et al. (1987) page 387                               ***
************************************************************************
      implicit double precision (a-h,o-z)

      dimension xmu(NDmu), wmu(NDmu), b(0:NDlay), ebmu(NDmu)
      dimension ssor(0:NDlay), ndoubl(NDlay), nlirf(0:NDlay+1)
      dimension ukk(NDsup,NDsup,NDirf,NDlay)
      dimension dkk(NDsup,NDsup,NDirf,NDlay)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension r(NDsup,NDsup), t(NDsup,NDsup)
      dimension sukk(NDsup,NDirf,NDlay), sdkk(NDsup,NDirf,NDlay)
      dimension sr(NDsup), st(NDsup)
      dimension u(NDsup,NDsup), d(NDsup,NDsup)
      dimension su(NDsup), sd(NDsup)
cmv   common /mat/ u,d

      nsup = nmat*nmutot
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine SHOMLY too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop 'dimension error in SHOMLY'
      endif
************************************************************************
* calculate also the radiation due to isotropically radiating sources  *
************************************************************************
      bs = b(il)*2.d0**(-ndoubl(il))
      call doad_expbmu(bs,xmu,NDmu,nmutot,ebmu)
cmv   write(*,*) 'fourier term = ',m,' layer il = ',il
cmv   write(*,*) 'bstart = ',bs,' ndoubl(il) = ',ndoubl(il)
      call doad_szero(nmat,nmutot,xmu,wmu,ebmu,bs,il,
     +    ssor,NDsup,NDmu,NDlay,sr,st)
      call doad_sone(nmat,nmutot,nmug,xmu,wmu,ebmu,bs,Zmmin,Zmplus,
     +    il,ssor,NDsup,NDmu,NDlay,sr,st)
      call doad_stwo(nmat,nmutot,nmug,xmu,wmu,ebmu,bs,Zmmin,Zmplus,
     +    il,ssor,NDsup,NDmu,NDlay,sr,st)
      call doad_firstm(nmat,nmutot,xmu,wmu,bs,ebmu,m,il,NDlay,
     +    Zmmin,Zmplus,NDmu,NDsup,r,t)
      call doad_seconm(nmat,nmutot,nmug,xmu,wmu,bs,ebmu,m,il,
     +    r,t,Zmmin,Zmplus,NDmu,NDsup,NDlay)
      nl2 = (nlirf(il)+1)/2
      do 20 n=1,ndoubl(il)
          call doad_sdoubl(nmat,nmutot,nmug,NDmu,NDsup,ebmu,r,t,u,d,
     +        sr,st,su,sd)
************************************************************************
* fill arrays ukk and dkk at all irf-levels in the isolated homogeneous*
* layer start when b=bs*2**(ndoubl(il)-nl2+1)                          *
* and calculate the reflection and transmission, r and t of that layer *
************************************************************************
          if ((n.gt.(ndoubl(il)-nl2)) .and. (nlirf(il).ne.0)) then
              nn = n-(ndoubl(il)-nl2)
              do 10 j=1,nsup
                sukk(j,nn,il) = su(j)
                sdkk(j,nn,il) = sd(j)
                do 10 i=1,nsup
                  ukk(i,j,nn,il) = u(i,j)
                  dkk(i,j,nn,il) = d(i,j)
   10         continue
          endif
          bs = 2.d0*bs
          call doad_expbmu(bs,xmu,NDmu,nmutot,ebmu)
   20 continue
      return
      end


      subroutine doad_smulay(nlayer,nmat,nmutot,nmug,unn,dnn,en,eb,
     +    sunn,sdnn,NDsup,NDlay,NDmu)
************************************************************************
*** radiation at the interfaces of a multilayered atmosphere         ***
*** the method is described in:                                      ***
*** Stammes, P., de Haan, J.F., Hovenier, J.W.:                      ***
*** 1989, Astron. Astrophys. xxx, p. xxx                             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension unn(NDsup,NDsup,NDlay+1)
      dimension dnn(NDsup,NDsup,NDlay+1)
      dimension sunn(NDsup,NDlay+1), sdnn(NDsup,NDlay+1)
      dimension en(NDmu,NDlay+1), eb(NDmu,NDlay)
************************************************************************
* on entrance the matrix unn(i,j,n) contains the supermatrix U, at     *
* level n of a partial atmosphere build in n adding steps              *
* on exit unn(i,j,n) contains the supermatrix U, at level n in the     *
* complete atmosphere                                                  *
* en(i,n) contains the factor exp[-tau(0,n)/u(i)]                      *
* eb(i,n) contains the factor exp[-b(n)/u(i)]                          *
* where, b(n) is the optical thickness of layer n                      *
* and tau(0,n) is the optical depth at level n                         *
************************************************************************

      dimension un(NDsup,NDsup), dn(NDsup,NDsup)
      dimension dnp1(NDsup,NDsup)
      dimension sdnp1(NDsup), sun(NDsup), sdn(NDsup)
      dimension em(NDmu), e(NDmu)
      dimension x(NDsup,NDsup), y(NDsup,NDsup), z(NDsup,NDsup)
      dimension sx(NDsup), sy(NDsup), sz(NDsup)
cmv   common /mat/ x,y,z,dnp1,un,dn

      if (nlayer .le. 1) return
      nsup = nmat*nmutot
      if (nmutot .gt. NDmu) then
          print*,'dimension in subroutine SMULAY too small'
          print*,'dim. is NDmu=',NDmu,' but nmutot=',nmutot
          stop ' dimension error in SMULAY'
      endif
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine SMULAY too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop ' dimension error in SMULAY'
      endif
      do 70 n=nlayer-1,1,-1
************************************************************************
* initializing arrays                                                  *
************************************************************************
          do 10 i=1,nmutot
              e(i) = en(i,n+1)
              em(i) = eb(i,n)
   10     continue
          do 20 i=1,nsup
            sdnp1(i) = sdnn(i,n+1)
            sun(i) = sunn(i,n)
            sdn(i) = sdnn(i,n)
            do 20 j=1,nsup
              dnp1(j,i) = dnn(j,i,n+1)
              un(j,i) = unn(j,i,n)
              dn(j,i) = dnn(j,i,n)
   20     continue
************************************************************************
* diffuse radiation field in the complete atmosphere                   *
* D field, Eq. (62)                                                    *
************************************************************************
          call doad_prodex(x,dnp1,em,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_prodsm(y,dn,dnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(z,x,y,nsup,NDsup)
          call doad_prodex(x,dn,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_addsm(y,x,z,nsup,NDsup)
          do 30 i=1,nsup
          do 30 j=1,nsup
              dnn(j,i,n) = y(j,i)
   30     continue
************************************************************************
* U field, Eq. (63)                                                    *
************************************************************************
          call doad_prodex(x,un,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_prodsm(y,un,dnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(z,x,y,nsup,NDsup)
          do 40 i=1,nsup
          do 40 j=1,nsup
              unn(j,i,n) = z(j,i)
   40     continue
************************************************************************
* diffuse radiation field from sources in the complete atmosphere      *
* D field, Eq. (50)                                                    *
************************************************************************
          call doad_prodes(sx,em,sdnp1,nmutot,nmug,nmat,NDsup,NDmu)
          call doad_prodmv(sy,dn,sdnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsv(sz,sx,sy,nsup,NDsup)
          call doad_addsv(sy,sdn,sz,nsup,NDsup)
          do 50 j=1,nsup
              sdnn(j,n) = sy(j)
   50     continue
************************************************************************
* U field, Eq. (51)                                                    *
************************************************************************
          call doad_prodmv(sy,un,sdnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsv(sz,sun,sy,nsup,NDsup)
          do 60 j=1,nsup
              sunn(j,n) = sz(j)
   60     continue
   70 continue
      return
      end


      subroutine doad_sone(nmat,nmutot,nmug,xmu,wmu,ebmu,b,Zmmin,Zmplus,
     +    il,ssor,NDsup,NDmu,NDlay,sr,st)
************************************************************************
*** calculate the first order radiation field of an isolated         ***
*** homogeneous layer containing internal sources                    ***
*** see Wauben (1990) Eqs. (23)-(31)                                 ***
*** the strength of the sources of layer il is contained in ssor(il) ***
*** the 0-th fourier coefficient of the phasematrix of layer il is   ***
*** contained in the arrays Zmmin(j,i,il) and Zmplus(j,i,il)         ***
*** ebmu contains exp[-b/xmu(i)] with b the optical thickness of     ***
*** isolated layer il                                                ***
*** on entrance sr and st contain the zero order emerging radiation  ***
*** on exit they contain the sum of the zero and first order         ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension ssor(0:NDlay), ebmu(NDmu)
      dimension sr(NDsup), st(NDsup)

      nsup = nmat*nmutot
      if (dabs(ssor(il)) .lt. 1.d-10) return
      if (b .lt. 1.d-12) return
      ag2 = ssor(il)/4.d0
************************************************************************
* loop 30 is the loop over each direction mu=mu_i                      *
************************************************************************
      do 30 i=1,nmutot
        im = nmat*(i-1)
        fact = ag2*wmu(i)
        xi = xmu(i)
        iz = 0
        if (xi .gt. 1.d-10) iz = 1
        bdi = 0.d0
        if (iz .eq. 1) bdi = b/xi
        ebi = ebmu(i)
        ei = 1.d0-ebi
        if (ei .lt. 1.d-3) ei=bdi*(1.d0-bdi/2.d0*(1.d0-bdi/3.d0))
************************************************************************
* gaussian integration over mu'=mu(k) in loop 20                       *
************************************************************************
        do 20 k=1,nmug
          xk = xmu(k)
          bdk = b/xk
          bik = bdi+bdk
          ebk = ebmu(k)
          eimek = ebi-ebk
          eiek = 1.d0-ebi*ebk
************************************************************************
* taylor expansion of exponent to avoid loss of accuracy if b << 1     *
************************************************************************
          if (eiek .lt. 1.d-3) then
            eiek = bik*(1.d0-bik*.5d0*(1.d0-bik/3.d0))
            eimek = bdk*(1.d0-bdk*.5d0*(1.d0-bdk/3.d0))-
     +        bdi*(1.d0-bdi*.5d0*(1.d0-bdi/3.d0))
          endif
          if (iz .eq. 0) then
            rp = 0.d0
            rm = 1.d0-ebk
          elseif ((bdi .lt. 1.d-3) .and.( bdk .lt. 1d-3)) then
************************************************************************
* taylor expansion of exponent to avoid loss of accuracy if b << 1     *
************************************************************************
            rp = 0.5d0*bdi*bdk*(1.d0-(bdk+2.d0*bdi)/3.d0)
            rm = 0.5d0*bdi*bdk*(1.d0-(bdk+bdi)/3.d0)
          else
************************************************************************
* determine the functions rp and rm see Wauben (1990) Eqs. (25)-(27)   *
* with mu->i and mu'->k                                                *
************************************************************************
            if (i .ne. k) then
              rm = ei+xk/(xk-xi)*eimek
            else
              rm = ei-bdi*ebi
            endif
            rp = ei-xk/(xk+xi)*eiek
          endif
          y = fact*wmu(k)*wmu(k)/xk
************************************************************************
* the 2*2 phase matrix obeys the symmetry relation                     *
* Z(-mu,-mu')=Z(mu,mu')  see Hovenier relation G                       *
************************************************************************
          k1 = nmat*(k-1)+1
          do 10 in=1,nmat
            ii = im+in
            zpik = Zmplus(ii,k1,il)
            zmik = Zmmin(ii,k1,il)
            sr(ii)=sr(ii)+y*(zmik*rp+zpik*rm)
   10     continue
   20   continue
   30 continue
************************************************************************
* loop 40 is the loop over each direction mu_i                         *
* using symmetry relation Eq.(22) the emerging radiation at the        *
* bottom equals the emerging radiation at the top                      *
************************************************************************
      do 40 i=1,nsup
        st(i) = sr(i)
   40 continue
      return
      end


      subroutine doad_soremb(nlirf,nl,nmat,nmutot,nmug,ukk,dkk,unn,
     +    dnn,ek,
     +    en,sukk,sdkk,sunn,sdnn,NDsup,NDirf,NDlay,NDmu)
************************************************************************
*** internal radiation field of an embedded homogeneous layer nl     ***
*** the method is described in:                                      ***
*** Stammes, P., de Haan, J.F., Hovenier, J.W.:                      ***
*** 1989, Astron. Astrophys. xxx, p. xxx                             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension ukk(NDsup,NDsup,NDirf,NDlay), unn(NDsup,NDsup,NDlay+1)
      dimension dkk(NDsup,NDsup,NDirf,NDlay), dnn(NDsup,NDsup,NDlay+1)
      dimension sukk(NDsup,NDirf,NDlay), sunn(NDsup,NDlay+1)
      dimension sdkk(NDsup,NDirf,NDlay), sdnn(NDsup,NDlay+1)
      dimension ek(NDmu,NDirf,NDlay), en(NDmu,NDlay+1), nlirf(0:NDlay+1)
************************************************************************
* on entrance the matrix ukk(i,j,k,n) contains the supermatrix U, at   *
* level k in the isolated homogeneous layer n                          *
* on exit ukk(i,j,k,n) contains the supermatrix U, at level k in the   *
* homogeneous layer n of the complete atmosphere                       *
* unn(i,j,n) contains the supermatrix U, at level n in the complete    *
* atmosphere                                                           *
* ek(i,k,n) contains the factor exp[-b(n)/u(i)*2**(k-M-1)] if k <= M   *
* and the factor exp[-b(n)/u(i)*2**(-k+M-1)] if k > M                  *
* en(i,n) contains the factor exp[-tau(0,n)/u(i)]                      *
* where b(n) is the optical thickness of layer n                       *
* and tau(0,n) is the optical depth at level n                         *
************************************************************************

      dimension uiso(NDsup,NDsup,NDirf), diso(NDsup,NDsup,NDirf)
      dimension uk(NDsup,NDsup), dk(NDsup,NDsup)
      dimension un(NDsup,NDsup), dnp1(NDsup,NDsup)
      dimension ukstar(NDsup,NDsup), dkstar(NDsup,NDsup)
      dimension eu(NDmu), ed(NDmu), e(NDmu)
      dimension x(NDsup,NDsup), y(NDsup,NDsup), z(NDsup,NDsup)
      dimension suk(NDsup), sdk(NDsup)
      dimension sun(NDsup), sdnp1(NDsup)
      dimension sx(NDsup), sy(NDsup), sz(NDsup)
cmv   common /mat/ x,y,z,dnp1,un,uk,dk,ukstar,dkstar,uiso,diso

      if (nlirf(nl) .eq. 0) return
      if (nmutot .gt. NDmu) then
          print*,'dimension in subroutine SOREMB too small'
          print*,'dim. is NDmu=',NDmu,' but nmutot=',nmutot
          stop ' dimension error in SOREMB'
      endif
      nsup = nmat*nmutot
      if (nsup .gt. NDsup) then
          print*,'dimension in subroutine SOREMB too small'
          print*,'dim. is NDsup=',NDsup,' but nsup=',nsup
          stop ' dimension error in SOREMB'
      endif
      if (nlirf(nl) .gt. NDirf) then
          print*,'dimension in subroutine SOREMB too small'
          print*,'dim. is NDirf=',NDirf,' but nlirf(nl)=',nlirf(nl)
          stop ' dimension error in SOREMB'
      endif
      do 5 k=1,nlirf(nl)
      do 5 i=1,nsup
      do 5 j=1,nsup
          uiso(j,i,k) = ukk(j,i,k,nl)
          diso(j,i,k) = dkk(j,i,k,nl)
    5 continue
      nl2 = (nlirf(nl)+1)/2
      do 50 k=1,nlirf(nl)
          kstar = nlirf(nl)+1-k
************************************************************************
* initializing arrays                                                  *
************************************************************************
          do 10 i=1,nmutot
              e(i) = en(i,nl+1)
              if (k .le. nl2) then
                  ed(i) = ek(i,k,nl)**(2.d0**(nl2-k+1)-1.d0)
                  eu(i) = ek(i,k,nl)
              else
                  ed(i) = ek(i,k,nl)
                  eu(i) = ek(i,k,nl)**(2.d0**(-nl2+k+1)-1.d0)
              endif
   10     continue
          do 20 i=1,nsup
            sdnp1(i) = sdnn(i,nl+1)
            sun(i) = sunn(i,nl)
            suk(i) = sukk(i,k,nl)
            sdk(i) = sdkk(i,k,nl)
            do 20 j=1,nsup
              dnp1(j,i) = dnn(j,i,nl+1)
              un(j,i) = unn(j,i,nl)
              uk(j,i) = uiso(j,i,k)
              dk(j,i) = diso(j,i,k)
              ukstar(j,i) = uiso(j,i,kstar)
              dkstar(j,i) = diso(j,i,kstar)
   20     continue
************************************************************************
* internal radiation field of embedded homogeneous layer nl            *
* D field, Eq. (72)                                                    *
************************************************************************
          call doad_transhm(ukstar,nmutot,nmat,NDsup)
          call doad_prodsm(x,ukstar,un,nmutot,nmug,nmat,NDsup)
          call doad_prodex(y,dnp1,ed,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_addsm(z,x,y,nsup,NDsup)
          call doad_prodsm(x,dk,dnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(y,z,x,nsup,NDsup)
          call doad_prodex(x,dk,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_addsm(z,y,x,nsup,NDsup)
          do 30 i=1,nsup
          do 30 j=1,nsup
              dkk(j,i,k,nl) = z(j,i)
   30     continue
************************************************************************
* internal radiation field of embedded homogeneous layer nl            *
* D field for sources Eq. (59)                                         *
************************************************************************
          call doad_prodmv(sx,ukstar,sun,nmutot,nmug,nmat,NDsup)
          call doad_prodes(sy,ed,sdnp1,nmutot,nmug,nmat,NDsup,NDmu)
          call doad_addsv(sz,sx,sy,nsup,NDsup)
          call doad_prodmv(sx,dk,sdnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsv(sy,sz,sx,nsup,NDsup)
          call doad_addsv(sz,sy,sdk,nsup,NDsup)
          do 35 j=1,nsup
              sdkk(j,k,nl) = sz(j)
   35     continue
************************************************************************
* U field, Eq. (73)                                                    *
************************************************************************
          call doad_transhm(dkstar,nmutot,nmat,NDsup)
          call doad_prodsm(x,dkstar,un,nmutot,nmug,nmat,NDsup)
          call doad_prodex(y,un,eu,nmutot,nmug,nmat,NDsup,NDmu,0)
          call doad_addsm(z,x,y,nsup,NDsup)
          call doad_prodsm(x,uk,dnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsm(y,z,x,nsup,NDsup)
          call doad_prodex(x,uk,e,nmutot,nmug,nmat,NDsup,NDmu,1)
          call doad_addsm(z,y,x,nsup,NDsup)
          do 40 i=1,nsup
          do 40 j=1,nsup
              ukk(j,i,k,nl) = z(j,i)
   40     continue
************************************************************************
* U field for sources Eq. (60)                                         *
************************************************************************
          call doad_prodmv(sx,dkstar,sun,nmutot,nmug,nmat,NDsup)
          call doad_prodes(sy,eu,sun,nmutot,nmug,nmat,NDsup,NDmu)
          call doad_addsv(sz,sx,sy,nsup,NDsup)
          call doad_prodmv(sx,uk,sdnp1,nmutot,nmug,nmat,NDsup)
          call doad_addsv(sy,sz,sx,nsup,NDsup)
          call doad_addsv(sz,sy,suk,nsup,NDsup)
          do 45 j=1,nsup
              sukk(j,k,nl) = sz(j)
   45     continue
   50 continue
      return
      end


      subroutine doad_stwo(nmat,nmutot,nmug,xmu,wmu,ebmu,b,Zmmin,Zmplus,
     +  il,ssor,NDsup,NDmu,NDlay,sr,st)
************************************************************************
*** calculate the second order internal radiation of an isolated     ***
*** homogeneous layer containing internal sources                    ***
*** the strength of the sources of layer il is contained in ssor(il) ***
*** the 0-th fourier coefficient of the phasematrix of layer il is   ***
*** contained in the arrays Zmmin(j,i,il) and Zmplus(j,i,il)         ***
*** ebmu contains exp[-b/xmu(i)] with b the optical thickness of     ***
*** isolated layer il                                                ***
*** on entrance sr and st contain the zero and first order radiation ***
*** on exit they contain the sum of the zero, first and second order ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu)
      dimension Zmmin(NDsup,NDsup,NDlay), Zmplus(NDsup,NDsup,NDlay)
      dimension ebmu(NDmu), ssor(0:NDlay)
      dimension sr(NDsup), st(NDsup)

      nsup = nmat*nmutot
      if (dabs(ssor(il)) .lt. 1.d-10) return
      if (b .lt. 1.d-12) return
************************************************************************
* loop 60 is also the loop over directions mu=mu_j                     *
* only the downward radiation is calculated due to symmetry            *
************************************************************************
      ag4 = ssor(il)/16.d0
      do 60 j=1,nmutot
        adj = ag4*wmu(j)
        jm = nmat*(j-1)
        xj = xmu(j)
        ebj = ebmu(j)
        jz = 0
        if (xj .gt. 1.d-10) jz = 1
        bdj = 0.d0
        if (jz .eq. 1) bdj = b/xj
************************************************************************
* gaussian integration over mu''=mu_k in loop 50                       *
************************************************************************
        do 50 k=1,nmug
          xk = xmu(k)
          adjk = adj*wmu(k)*wmu(k)/xk
          km = nmat*(k-1)
          ebk = ebmu(k)
          xjpxk = xj+xk
          xjmxk = xj-xk
          bdk = b/xk
          bdkj = bdj+bdk
          ebkmj = ebk-ebj
          ebkj = 1.d0-ebk*ebj
          if (ebkj .lt. 1.d-3) then
            ebkj = bdkj*(1.d0-bdkj*.5d0*(1.d0-bdkj/3.d0))
            ebkmj = bdj*(1.d0-bdj*.5d0*(1.d0-bdj/3.d0))-
     +          bdk*(1.d0-bdk*.5d0*(1.d0-bdk/3.d0))
          endif
************************************************************************
* gaussian integration over mu'=mu_i in loop 50                        *
************************************************************************
          do 50 i=1,nmug
            im = nmat*(i-1)
            xi = xmu(i)
            ebi = ebmu(i)
            bdi = b/xi
************************************************************************
* mu(j)=0 special case                                                 *
************************************************************************
            if (jz .eq. 0) then
              xipxk = xi+xk
              ximxk = xi-xk
              rpp = 0.d0
              rpm = 0.d0
              rmp = 1.d0-xi/xipxk-(1.d0-xi/xipxk*ebi)*ebk
              if (i .ne. k) then
                rmm = 1.d0-xi/ximxk*ebi+xk/ximxk*ebk
              else
                rmm = 1.d0-ebi-bdi*ebi
              endif
************************************************************************
* b<<1 special case                                                    *
************************************************************************
            elseif ((bdj .lt. 1.d-3) .and. (bdi .lt. 1.d-3) .and.
     +        (bdk .lt. 1d-3)) then
              rpp = bdj*bdi*bdk/6.d0
              rpm = rpp
              rmp = rpp
              rmm = rpp
            else
              xipxj = xi+xj
              ximxj = xi-xj
              xipxk = xi+xk
              ximxk = xi-xk
              ebkmi = ebk-ebi
              ebimj = ebi-ebj
              ebij = 1.d0-ebi*ebj
              ebki = 1.d0-ebk*ebi
	      bdij = bdi+bdj
	      bdki = bdi+bdk
              if (ebij .lt. 1.d-3) then
                ebij = bdij*(1.d0-bdij*.5d0*(1.d0-bdij/3.d0))
                ebimj = bdj*(1.d0-bdj*.5d0*(1.d0-bdj/3.d0))-
     +            bdi*(1.d0-bdi*.5d0*(1.d0-bdi/3.d0))
              endif
              if (ebki .lt. 1.d-3) then
                ebki = bdki*(1.d0-bdki*.5d0*(1.d0-bdki/3.d0))
                ebkmi = bdi*(1.d0-bdi*.5d0*(1.d0-bdi/3.d0))-
     +            bdk*(1.d0-bdk*.5d0*(1.d0-bdk/3.d0))
              endif
************************************************************************
* determine the functions r++, r+-, r-+, r--                           *
* see Wauben (????) Eqs. (?)-(??) with mu->j, mu'->i and mu''->k       *
************************************************************************
              if (i .eq. j) then
                if (i .ne. k) then
                  rpp = 1.d0-ebj-xi/ximxk*xi/xipxj*ebij
     +               +xk/ximxk*xk/xjpxk*ebkj
                  rpm = 1.d0-ebj-bdj*ebj*xj/xjpxk
     +               -(1.d0-xj/xjpxk*ebj)*xk/xjpxk*ebkj
                  rmp = 1.d0-ebj-xi/xipxj*xi/xipxk*ebij
     +               +(1.d0-xi/xipxk*ebi)*xk/xjmxk*ebkmj
                  rmm = 1.d0-ebj-b/xjmxk*ebj
     +               -xk/xjmxk*xk/xjmxk*ebkmj
                else
                  rpp = 1.d0-ebj+b/xipxj-((2.d0*xi+b)/xipxj
     +               -xi/xipxj*xi/xipxj)*ebij
                  rpm = 1.d0-ebj-bdj*ebj*xj/xjpxk
     +               -(1.d0-xj/xjpxk*ebj)*xk/xjpxk*ebkj
                  rmp = 1.d0-ebj-bdj*ebj*(1.d0-ebi*xi/xipxj)
     +               -xi/xipxj*xi/xipxj*ebij
                  rmm = 1.d0-ebj-bdj*ebj-0.5d0*bdj*bdj*ebj
                endif
              elseif (i .eq. k) then
                rpp = 1.d0-ebj+b/xipxj-((2.d0*xi+b)/xipxj
     +             -xi/xipxj*xi/xipxj)*ebij
                rpm = 1.d0-ebj-xi/xipxk*xi/ximxj*ebimj
     +             -(1.d0-xi/xipxk*ebi)*xk/xjpxk*ebkj
                rmp = 1.d0-ebj-xi/xipxj*xi/xipxk*ebij
     +             +(1.d0-xi/xipxk*ebi)*xk/xjmxk*ebkmj
                rmm = 1.d0-ebj-b/ximxj*ebi
     +             +(xi/ximxj*xi/ximxj-2.d0*xi/ximxj)*ebimj
              elseif (j .eq. k) then
                rpp = 1.d0-ebj-xi/ximxk*xi/xipxj*ebij
     +             +xk/ximxk*xk/xjpxk*ebkj
                rpm = 1.d0-ebj-xi/xipxk*xi/ximxj*ebimj
     +             -(1.d0-xi/xipxk*ebi)*xk/xjpxk*ebkj
                rmp = 1.d0-ebj-bdj*ebj*(1.d0-ebi*xi/xipxj)
     +             -xi/xipxj*xi/xipxj*ebij
                rmm = 1.d0-ebj+b/ximxj*ebj
     +             -xi/ximxj*xi/ximxj*ebimj
              else
                rpp = 1.d0-ebj-xi/ximxk*xi/xipxj*ebij
     +             +xk/ximxk*xk/xjpxk*ebkj
                rpm = 1.d0-ebj-xi/xipxk*xi/ximxj*ebimj
     +             -(1.d0-xi/xipxk*ebi)*xk/xjpxk*ebkj
                rmp = 1.d0-ebj-xi/xipxj*xi/xipxk*ebij
     +             +(1.d0-xi/xipxk*ebi)*xk/xjmxk*ebkmj
                rmm = 1.d0-ebj-xi/ximxj*xi/ximxk*ebimj
     +             -xk/ximxk*xk/xjmxk*ebkmj
              endif
            endif
            y = wmu(i)*wmu(i)*adjk/xi
************************************************************************
* the matrixproduct z^m_{j,k}*z^m_{k,i}*e_1                            *
* where Z^m(u,u') = Z^m(-u,-u')                                        *
************************************************************************
            do 40 jl=1,nmat
              jjl = jm+jl
              do 40 kl=1,nmat
                kkl = km+kl
                iil = im+1
                zpjk = Zmplus(jjl,kkl,il)
                zpki = Zmplus(kkl,iil,il)
                zmjk = Zmmin(jjl,kkl,il)
                zmki = Zmmin(kkl,iil,il)
                sr(jjl)=sr(jjl)+y*( zmjk*zpki*rpp+zmjk*zmki*rpm
     +            +zpjk*zmki*rmp+zpjk*zpki*rmm )
   40       continue
   50   continue
   60 continue
************************************************************************
* use symmetry relation I(b-tau,-mu) = I(tau,mu)                       *
************************************************************************
      do 70 j=1,nsup
        st(j) = sr(j)
   70 continue
      return
      end


      subroutine doad_sumfou(nmat,nmug,nmuext,xmu,wmu,nphout,
     +    phout,nlirf,
     +    nlayer,nmu0,imu0,m,ukk,unn,dkk,dnn,usum,dsum,fin,
     +    NDmu,NDmu0,NDsup,NDphi0,NDlay,NDirf)
************************************************************************
*** summation of the fourier series at phout and with mu0            ***
*** for mu0 some(nmu0) values are selected from the extra mu-points  ***
*** for mu all nmug+nmuext mu-points are used                        ***
*** the sums are stored in usum and dsum (mu,mu0,phout,irf,il)       ***
*** the first nmutot elements of the sums contain the intensity I    ***
*** the next nmutot elements the Stokes parameter Q etc.             ***
*** the incident flux vector is given by fin{1,2,3,4}                ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu), phout(NDphi0), imu0(NDmu0)
      dimension ukk(NDsup,NDsup,NDirf,NDlay), unn(NDsup,NDsup,NDlay+1)
      dimension dkk(NDsup,NDsup,NDirf,NDlay), dnn(NDsup,NDsup,NDlay+1)
      dimension dsum(NDsup,NDmu0,NDphi0,0:NDirf+1,NDlay+1)
      dimension usum(NDsup,NDmu0,NDphi0,0:NDirf+1,NDlay+1)
      dimension nlirf(0:NDlay+1), fin(4)
      nmutot = nmug+nmuext
      nsup = nmat*nmutot
      ns = nmug+1
      if (ns .gt. nmutot) return
      pi = 4.d0*datan(1.d0)
      if (m .eq. 0) then
        do 10 il=1,NDlay+1
        do 10 kirf=0,NDirf
        do 10 ip=1,nphout
        do 10 j=1,nmu0
        do 10 i=1,nsup
          usum(i,j,ip,kirf,il) = 0.d0
          dsum(i,j,ip,kirf,il) = 0.d0
   10   continue
      endif

      do 20 ip=1,nphout
        cmph = dcos(m*phout(ip)*pi/180.d0)
        smph = dsin(m*phout(ip)*pi/180.d0)
        if (m .ne. 0) then
          cmph = 2.d0*cmph
          smph = 2.d0*smph
        endif
        do 20 i=1,nmutot
          im = (i-1)*nmat
          do 20 j=1,nmu0
            ji = imu0(j)+nmug
            xm0 = xmu(ji)
            wij = wmu(i)*wmu(ji)
************************************************************************
* loop over elements of incident flux vector                           *
************************************************************************
            do 20 ic=1,nmat
              j1 = (ji-1)*nmat+ic
              do 20 k=1,nmat
                ik = im+k
                iik = (k-1)*nmutot+i
                y = cmph
                if (k .gt. 2) y = smph
                if (ic .gt. 2) then
                  y = -smph
                  if (k .gt. 2) y = cmph
                endif
                fac = y*fin(ic)/wij
                do 20 il=1,nlayer+1

!       if( ( il .eq. nlayer + 1 ).and.(iik .gt. nmug).and. 
!     &  (usum(iik,j,ip,0,il).lt.0)) then
!             write(*,*) 'usum before ', usum(iik,j,ip,0,il)
!        end if

                  usum(iik,j,ip,0,il) = usum(iik,j,ip,0,il)+
     +              xm0*unn(ik,j1,il)*fac

!      if( ( il .eq. nlayer + 1 ).and.(iik .gt. nmug).and. 
!     &  (usum(iik,j,ip,0,il).lt.0)) then
!       write(*,*) 'iik,j,ip, dsum,nmuext, xm0,dnn,fac, unn, usum ',
!     &   iik,j,ip,il, usum(iik,j,ip,0,il),nmuext, xm0,unn(ik,j1,il),fac
!      end if
!        if( ( il .eq.nlayer +  1 ).and.(iik .gt. nmug).and. 
!     & ( usum(iik,j,ip,0,il) .gt. 1000 ) ) then
!           write(*,*) 'iik,j,ip, dsum,dnn,fac',iik,j,ip,il,
!     & usum(iik,j,ip,0,il),unn(ik,j1,il),fac
!        end if


     
!      if( ( il .eq. 1 ).and.(iik .gt. nmug).and. 
!     &  (dsum(iik,j,ip,0,il).lt.0)) then
!             write(*,*) 'dsum before ', dsum(iik,j,ip,0,il)
!        end if

                  dsum(iik,j,ip,0,il) = dsum(iik,j,ip,0,il)+
     +              xm0*dnn(ik,j1,il)*fac

!        if( ( il .eq. 1 ).and.(iik .gt. nmug).and. 
!     &  (dsum(iik,j,ip,0,il).lt.0)) then
!       write(*,*) 'iik,j,ip, dsum,nmuext, xm0,dnn,fac, unn, usum ',
!     &   iik,j,ip,il, dsum(iik,j,ip,0,il),nmuext, xm0,dnn(ik,j1,il),fac
!        end if
!        if( ( il .eq. 1 ).and.(iik .gt. nmug).and. 
!     & ( dsum(iik,j,ip,0,il) .gt. 1000 ) ) then
!           write(*,*) 'iik,j,ip, dsum,dnn,fac',iik,j,ip,il,
!     & dsum(iik,j,ip,0,il),dnn(ik,j1,il),fac
!        end if
                  if (nlirf(il) .ne. 0) then
                    do 30 kirf=1,nlirf(il)
                      usum(iik,j,ip,kirf,il) = usum(iik,j,ip,kirf,il)+
     +                  xm0*ukk(ik,j1,kirf,il)*fac
                      dsum(iik,j,ip,kirf,il) = dsum(iik,j,ip,kirf,il)+
     +                  xm0*dkk(ik,j1,kirf,il)*fac
   30               continue
                  endif
   20   continue
      return
      end


      subroutine doad_sumsfm(nmat,nmug,nmuext,xmu,wmu,nlirf,
     +    nlayer,m,sukk,sunn,sdkk,sdnn,susum,sdsum,
     +    NDmu,NDsup,NDlay,NDirf)
************************************************************************
*** summation of the fourier series of source super vectors          ***
*** for mu all nmug+nmuext mu-points are used                        ***
*** the sums are stored in susum and sdsum (mu,irf,il)               ***
*** the first nmutot elements of the sums contain the intensity I    ***
*** the next nmutot elements the Stokes parameter Q etc.             ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu)
      dimension sukk(NDsup,NDirf,NDlay), sunn(NDsup,NDlay+1)
      dimension sdkk(NDsup,NDirf,NDlay), sdnn(NDsup,NDlay+1)
      dimension sdsum(NDsup,0:NDirf+1,NDlay+1)
      dimension susum(NDsup,0:NDirf+1,NDlay+1)
      dimension nlirf(0:NDlay+1)

      nmutot = nmug+nmuext
      nsup = nmat*nmutot
cmv -- modified to obtain correct output (sor.out) for I Q U V
cmv -- The variable ns is an artifact from the case for external illumination, 
cmv -- since only additional zenith angles are considered as solar zenith 
cmv -- angles. Therefore, the statements
cmv --ns = nmug+1
cmv --if (ns .gt. nmutot) return
cmv -- in sumsfm should be deleted, because otherwise a run without additional
cmv -- zenith angles (i.e. nmuext=0 and therefore nmutot=nmug) would result
cmv -- in zeros
cmv   ns = nmug+1
cmv   if (ns .gt. nmutot) return
      pi = 4.d0*datan(1.d0)
      if (m .eq. 0) then
        do 10 il=1,NDlay+1
        do 10 kirf=0,NDirf
        do 10 i=1,nsup
          susum(i,kirf,il) = 0.d0
          sdsum(i,kirf,il) = 0.d0
   10   continue
      endif
      cmph = 1.d0
      do 20 i=1,nmutot
        im = (i-1)*nmat
        do 20 k=1,nmat
          ik = im+k
          iik = (k-1)*nmutot+i
          do 20 il=1,nlayer+1
            susum(iik,0,il) = susum(iik,0,il)+
     +        sunn(ik,il)*cmph/wmu(i)
            sdsum(iik,0,il) = sdsum(iik,0,il)+
     +        sdnn(ik,il)*cmph/wmu(i)
            if (nlirf(il) .ne. 0) then
              do 30 kirf=1,nlirf(il)
                susum(iik,kirf,il) = susum(iik,kirf,il)+
     +            sukk(ik,kirf,il)*cmph/wmu(i)
                sdsum(iik,kirf,il) = sdsum(iik,kirf,il)+
     +            sdkk(ik,kirf,il)*cmph/wmu(i)
   30         continue
          endif
   20 continue
      return
      end



      subroutine doad_symmut(S,nmat,nmutot,NDsup,id)
************************************************************************
*** fill upper-triangle of supermatrix S, using symmetry relations   ***
*** see Hovenier, 1969                                               ***
*** id = 3: reflection      relation g                               ***
*** id = 4: transmission    relation h                               ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension S(NDsup,NDsup)

      nsup = nmat*nmutot
************************************************************************
* fill upper-triangle with its transposed values                       *
* and change signs of each id-th column                                *
************************************************************************
      do 30 i=1,nmutot
          ni = nmat*(i-1)
          do 20 k=1,nmat
              nk = ni+k
              y = 1.d0
              if (k .eq. id) y = -1.d0
              do 10 j=1,nk
                  S(j,nk) = y*S(nk,j)
   10         continue
   20     continue
   30 continue
      if (nmat .lt. id) return
************************************************************************
* change id-th row of sign                                             *
************************************************************************
      do 50 i=1,nmutot
          nid = (i-1)*nmat+id
          do 40 j=nid,nsup
              S(nid,j) = -S(nid,j)
   40     continue
   50 continue
      return
      end


      subroutine doad_szero(nmat,nmutot,xmu,wmu,ebmu,b,il,
     +    ssor,NDsup,NDmu,NDlay,sr,st)
************************************************************************
*** calculate the zero order radiation field of an isolated          ***
*** homogeneous layer due to isotropically radiating internal sources***
*** the strength of the sources of layer il is contained in ssor(il) ***
*** see Wauben (1990)                                                ***
*** results are stored in supervectors sr(i,il) and st(i,il)         ***
*** ebmu contains exp[-b/xmu(i)] with b the optical thickness of     ***
*** isolated layer il                                                ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu)
      dimension ssor(0:NDlay), ebmu(NDmu)
      dimension sr(NDsup), st(NDsup)

************************************************************************
* initializing the supervectors sr and st                              *
************************************************************************
      nsup = nmat*nmutot
      do 10 i=1,nsup
        sr(i) = 0.d0
        st(i) = 0.d0
   10 continue
      if (dabs(ssor(il)) .lt. 1.d-10) return
************************************************************************
* infinitesimally narrow source layer                                  *
************************************************************************
      if (b .lt. 1.d-12) then
        do 15 i=1,nmutot
          if (xmu(i) .gt. 1.d-10) then
            sr(nmat*(i-1)+1) = wmu(i)*ssor(il)/xmu(i)
            st(nmat*(i-1)+1) = wmu(i)*ssor(il)/xmu(i)
          else
            sr(nmat*(i-1)+1) = 5.d55
            st(nmat*(i-1)+1) = 5.d55
          endif
   15   continue
        return
      endif
************************************************************************
* loop 20 over the directions mu_i                                     *
************************************************************************
      do 20 i=1,nmutot
        im = nmat*(i-1)+1
        fact = ssor(il)*wmu(i)
        ei = 1.d0-ebmu(i)
        if (ei .gt. 1.d-3) then
************************************************************************
* unpolarized zero order emerging radiation see Wauben (1990), Eq.(10) *
* radiation emerging from top equals radiation emerging from bottom    *
************************************************************************
          sr(im) = fact*ei
          st(im) = sr(im)
        else
************************************************************************
* expansion in Taylor series to avoid loss of accuracy if b<<1         *
************************************************************************
          bdi = b/xmu(i)
          sr(im) = fact*bdi*(1.d0-bdi/2.d0*(1.d0-bdi/3.d0))
          st(im) = sr(im)
        endif
   20 continue
      return
      end


      subroutine doad_transf(Z,il,NDlay,NDsup,nsup,nmat)
************************************************************************
* Calculate the product D1 * Z * D1 where D1 is a supermatrix          *
* containing in each submatrix :                                       *
*         ( 1   0   0   0 )          ( 1   0   0 )                     *
*         ( 0   1   1   0 )    or    ( 0   1   1 )                     *
*         ( 0   1  -1   0 )          ( 0   1  -1 )                     *
*         ( 0   0   0   1 )                                            *
************************************************************************
      implicit double precision (a-h,o-z)
      dimension Z(NDsup,NDsup,NDlay)

      do 30 i=2,nsup,nmat
          do 10 j=1,nsup
              s1 = Z(i,j,il)
              s2 = Z(i+1,j,il)
              Z(i,j,il) = s1+s2
              Z(i+1,j,il) = s1-s2
   10     continue
          do 20 j=1,nsup
              s1 = Z(j,i,il)
              s2 = Z(j,i+1,il)
              Z(j,i,il) = s1+s2
              Z(j,i+1,il) = s1-s2
   20     continue
   30 continue
      return
      end


      subroutine doad_transhm(x,nmutot,nmat,NDsup)
************************************************************************
*** transformation of "light incident at the top" to "light incident ***
*** at the bottom", for reflection and transmission in case of       ***
*** homogeneous layers only                                          ***
*** see Hovenier (1969) relations q and r                            ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension x(NDsup,NDsup)

      if (nmat .lt. 3) return
************************************************************************
* change signs of each non-diagonal element of the 3-th and 4-th       *
* row and column                                                       *
************************************************************************
      nt = nmutot*nmat
      do 30 i=1,nmutot
          ni = nmat*(i-1)
          do 20 k=3,nmat
              nik = ni+k
              do 10 j=1,nt
                  x(nik,j) = -x(nik,j)
                  x(j,nik) = -x(j,nik)
   10         continue
   20     continue
   30 continue
      return
      end



      subroutine doad_xeqy(x,y,NDsup)
************************************************************************
*** set supermatrix x(NDsup,NDsup) equal to y                        ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension x(NDsup,NDsup), y(NDsup,NDsup)

      do 10 j=1,NDsup
      do 10 i=1,NDsup
          x(i,j) = y(i,j)
   10 continue
      return
      end


      subroutine doad_fluxes(nmat,nmutot,nmug,xmu,wmu,m,nlayer,nlirf,
     +    imu0,nmu0,fluxu,fluxd,etmu0,
     +    ukk,dkk,unn,dnn,tau,en,NDmu,NDsup,NDirf,NDlay,NDmu0)
************************************************************************
*** calculates the upward (fluxu) and the downward (fluxd) fluxes    ***
*** and the directly transmitted light (etmu0) for incident light    ***
*** with mu0 at all the interface layers and irf levels              ***
*** see Stammes et al. 1989, here F(tau) is calculated               ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu)
      dimension ukk(NDsup,NDsup,NDirf,NDlay), unn(NDsup,NDsup,NDlay+1)
      dimension dkk(NDsup,NDsup,NDirf,NDlay), dnn(NDsup,NDsup,NDlay+1)
      dimension tau(0:NDirf+1,0:NDlay+1), en(NDmu,NDlay+1)
      dimension nlirf(0:NDlay+1)
      dimension fluxu(NDmu0,0:NDirf,NDlay+1)
      dimension fluxd(NDmu0,0:NDirf,NDlay+1)
      dimension etmu0(NDmu0,0:NDirf,NDlay+1), imu0(NDmu0)

      if (m .ne. 0) return
      do 5 il=1,nlayer+1
      do 5 k=0,nlirf(il)
      do 5 j=1,nmu0
        fluxu(j,k,il) = 0.d0
        fluxd(j,k,il) = 0.d0
        etmu0(j,k,il) = 0.d0
    5 continue
      do 60 il=1,nlayer+1
        do 20 j=1,nmu0
          ij = imu0(j)+nmug
          ij1 = nmat*(ij-1)+1
          xj = xmu(ij)
          wj = wmu(ij)
          do 10 i=1,nmug
            i1 = nmat*(i-1)+1
            wi = wmu(i)
            fluxu(j,0,il) = fluxu(j,0,il)+wi*unn(i1,ij1,il)
            fluxd(j,0,il) = fluxd(j,0,il)+wi*dnn(i1,ij1,il)
   10     continue
          fluxu(j,0,il) = xj*fluxu(j,0,il)/wj
          fluxd(j,0,il) = xj*fluxd(j,0,il)/wj
          etmu0(j,0,il) = xj*en(ij,il)
   20   continue
        if (nlirf(il) .ne. 0) then
          do 50 k=1,nlirf(il)
            do 40 j=1,nmu0
              ij = imu0(j)+nmug
              ij1 = nmat*(ij-1)+1
              xj = xmu(ij)
              wj = wmu(ij)
              do 30 i=1,nmug
                i1 = nmat*(i-1)+1
                wi = wmu(i)
                fluxu(j,k,il) = fluxu(j,k,il)+wi*ukk(i1,ij1,k,il)
                fluxd(j,k,il) = fluxd(j,k,il)+wi*dkk(i1,ij1,k,il)
   30         continue
              fluxu(j,k,il) = xj*fluxu(j,k,il)/wj
              fluxd(j,k,il) = xj*fluxd(j,k,il)/wj
              if (xj .gt. 1.d-10) then
                if (tau(k,il)/xj .lt. 2.d+2) then
                  etmu0(j,k,il) = xj*dexp(-tau(k,il)/xj)
                else
                  etmu0(j,k,il) = 0.d0
                endif
              else
                etmu0(j,k,il) = 0.d0
              endif
   40       continue
   50     continue
        endif
   60 continue
      return
      end


      subroutine doad_fluxrt(iu,nmat,nmutot,nmug,xmu,wmu,m,nlayer,nlirf,
     +    ukk,dkk,unn,dnn,tau,en,NDmu,NDsup,NDirf,NDlay)
************************************************************************
*** calculates the upward (frj) and the downward (ftj) fluxes for    ***
*** the total atmosphere and the directly transmitted light (ej)     ***
*** for incident light with direction cosine mu_j                    ***
*** see Stammes et al. 1989 calculated is F/mu0                      ***
*** the bond albedo is calculated according to V.d. Hulst 1980 p.603 ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu)
      dimension ukk(NDsup,NDsup,NDirf,NDlay), unn(NDsup,NDsup,NDlay+1)
      dimension dkk(NDsup,NDsup,NDirf,NDlay), dnn(NDsup,NDsup,NDlay+1)
      dimension tau(0:NDirf+1,0:NDlay+1), en(NDmu,NDlay+1)
      dimension nlirf(0:NDlay+1)

      if (m .ne. 0) return
c$$$      write(iu,5) tau(0,1), nlayer
c$$$    5 format(1h1,'Fluxes for atmosphere with optical thickness ',
c$$$     +  f4.1,/,' consisting of ',i2,' homegeneous layer(s)')
c$$$      write(iu,15)
c$$$   15 format(1h0,11x,'mu0',15x,'fr(mu0)',13x,'ft(mu0)',13x,'sum')
      do 20 j=nmug+1,nmutot
        xj = xmu(j)
        j1 = nmat*(j-1)+1
        wj = wmu(j)
        frj = 0.d0
        ftj = 0.d0
        do 10 i=1,nmug
          i1 = nmat*(i-1)+1
          wi = wmu(i)
          frj = frj+wi*unn(i1,j1,nlayer+1)
          ftj = ftj+wi*dnn(i1,j1,1)
   10   continue
        frj = frj/wj
        ftj = ftj/wj
        ej = en(j,1)
        sum = frj+ftj+ej
c$$$        write(iu,25) xj,frj,ftj,sum
c$$$	write(43,25) xj,frj
   25   format(4f20.12)
   20 continue
************************************************************************
* the bond (spherical) albedo                                          *
************************************************************************
      ar = 0.d0
      at = 0.d0
      do 30 i=1,nmug
        i1 = nmat*(i-1)+1
        wi = wmu(i)
        do 30 j=1,nmug
          j1 = nmat*(j-1)+1
          wj = wmu(j)
          ar = ar+wi*unn(i1,j1,nlayer+1)*wj
          at = at+wi*dnn(i1,j1,1)*wj
   30 continue
c$$$      write(iu,35) ar,at
c$$$   35 format(1h0,' uru(bond albedo)=',f20.14,10x,'utu=',f20.14)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            BELOW ARE STORED UNUSED ROUTINES 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c$$$
c$$$
c$$$      subroutine doad_writsc(iunit,alfbet,NDcoef,NDlay,il,lmax,a,csth,
c$$$     +    Qext,Qsca,G,V,s1,s2,s3,s4,s5,s6)
c$$$************************************************************************
c$$$*  On entry :                                                          *
c$$$*      iunit     unit number of the coefficient file                   *
c$$$*      alfbet    array containing the coefficients, according to :     *
c$$$*                                                                      *
c$$$*                         ( alpha1(l) beta1(l)     0         0      )  *
c$$$*      alfbet(i,j,l)   =  ( beta1(l)  alpha2(l)    0         0      )  *
c$$$*                         (   0         0       alpha3(l) beta2(l)  )  *
c$$$*                         (   0         0      -beta2(l)  alpha4(l) )  *
c$$$*                                                                      *
c$$$*      NDcoef    third dimension of array alfbet (NDcoef >= lmax)      *
c$$$*      lmax      maximum index l for which coefficients are non-zero   *
c$$$*      a         albedo for single scattering                          *
c$$$*      csth      asymmetry parameter <cos th>                          *
c$$$*      Qext      extinction efficiency                                 *
c$$$*      Qsca      scattering efficiency                                 *
c$$$*      G         average geometrical cross section                     *
c$$$*      V         average volume                                        *
c$$$*      s1,...,s6 strings with description of scattering matrix         *
c$$$*  On exit :                                                           *
c$$$*      The coefficients and other information is written on the file   *
c$$$*      with unit number iunit in the standard format.                  *
c$$$*                                                 VLD January 9, 1989  *
c$$$************************************************************************
c$$$      implicit double precision (a-h,o-z)
c$$$      parameter ( eps = 1.d-6 )
c$$$
c$$$      dimension alfbet(4,4,0:NDcoef,NDlay)
c$$$      character*60 s1,s2,s3,s4,s5,s6
c$$$
c$$$      if (lmax .gt. NDcoef) then
c$$$          print *,' writsc: lmax = ',lmax,' maximum NDcoef = ',NDcoef
c$$$          stop 'in writsc too many coefficients requested'
c$$$      endif
c$$$      if (dabs(csth-alfbet(1,1,1,il)/3.d0) .gt. eps) then
c$$$          print *,' writsc: warning:    asymmetry parameter = ',csth
c$$$          print *,'         1/3 first order legendre coeff. = ',
c$$$     +          alfbet(1,1,1,il)/3.d0
c$$$      endif
c$$$      if ((a .lt. 0.d0) .or. (a .gt. 1.d0)) then
c$$$          print *,' writsc: warning: single scattering albedo = ',a,
c$$$     +          ' outside interval [0,1]'
c$$$      endif
c$$$      if (Qext .lt. 0.d0) then
c$$$          print *,' writsc: warning: extinction efficiency = ',Qext,
c$$$     +          ' negative !'
c$$$      endif
c$$$      if (Qsca .lt. 0.d0) then
c$$$          print *,' writsc: warning: scattering efficiency = ',Qsca,
c$$$     +          ' negative !'
c$$$      endif
c$$$      if (Qsca .gt. Qext) then
c$$$          print *,' writsc: warning:    scattering efficiency = ',Qsca
c$$$          print *,'         larger than extinction efficiency = ',Qext
c$$$      endif
c$$$      if (G .lt. 0.) then
c$$$          print *,' writsc: warning: geometrical cross section = ',G
c$$$     +                                                   ,' negative !'
c$$$      endif
c$$$      if (V .lt. 0.) then
c$$$          print *,' writsc: warning: average volume = ',V
c$$$     +                                                   ,' negative !'
c$$$      endif
c$$$      if (dabs(alfbet(1,1,0,il)-1.d0) .gt. eps) then
c$$$          print *,' writsc: warning: zero legendre coeff = ',
c$$$     +          alfbet(1,1,0,il)
c$$$          print *,'         not equal to one !'
c$$$      endif
c$$$      write(iunit,1)
c$$$    1 format(1h1,' EXPANSION COEFFICIENTS SCATTERING MATRIX')
c$$$      write(iunit,2) s1
c$$$    2 format(' ',a60)
c$$$      write(iunit,2) s2
c$$$      write(iunit,2) s3
c$$$      write(iunit,2) s4
c$$$      write(iunit,2) s5
c$$$      write(iunit,2) s6
c$$$      write(iunit,3) a,csth,Qext,Qsca,G,V
c$$$    3 format(' single scattering albedo  a  =',e22.14,/
c$$$     +      ,' asymmetry parameter <cos th> =',e22.14,/
c$$$     +      ,' extinction efficiency  Qext  =',e22.14,/
c$$$     +      ,' scattering efficiency  Qsca  =',e22.14,/
c$$$     +      ,' geometrical cross section G  =',e22.14,/
c$$$     +      ,' average volume            V  =',e22.14,/
c$$$     +      ,'   l'
c$$$     +      ,'       alpha1             alpha2      '
c$$$     +      ,'       alpha3             alpha4      '
c$$$     +      ,'       beta1              beta2       ')
c$$$      do 10 l=0,lmax
c$$$          write(iunit,4) l,alfbet(1,1,l,il),alfbet(2,2,l,il),
c$$$     +        alfbet(3,3,l,il),alfbet(4,4,l,il),
c$$$     +        alfbet(1,2,l,il),alfbet(3,4,l,il)
c$$$    4     format(i4,6f19.14)
c$$$          if ((dabs(alfbet(1,3,l,il)) .gt. eps) .or.
c$$$     +        (dabs(alfbet(1,4,l,il)) .gt. eps) .or.
c$$$     +        (dabs(alfbet(2,3,l,il)) .gt. eps) .or.
c$$$     +        (dabs(alfbet(2,4,l,il)) .gt. eps) .or.
c$$$     +        (dabs(alfbet(3,1,l,il)) .gt. eps) .or.
c$$$     +        (dabs(alfbet(3,2,l,il)) .gt. eps) .or.
c$$$     +        (dabs(alfbet(4,1,l,il)) .gt. eps) .or.
c$$$     +        (dabs(alfbet(4,2,l,il)) .gt. eps) )
c$$$     +            print *,' writsc: warning: nonzero coeff. in upper'
c$$$     +                   ,' right or lower left corner l = ',l
c$$$   10 continue
c$$$      write(iunit,4) -1,0.,0.,0.,0.,0.,0.
c$$$      return
c$$$      end
c$$$
c$$$
c$$$
c$$$      subroutine doad_clossc(iunit)
c$$$************************************************************************
c$$$*  On entry :                                                          *
c$$$*      iunit     number of the unit to be closed                       *
c$$$*  On exit :                                                           *
c$$$*      The file is closed.                                             *
c$$$*                                                 VLD January 9, 1989  *
c$$$************************************************************************
c$$$      close(unit=iunit,err=999)
c$$$      return
c$$$  999 print *,' clossc: error in closing file, unit=',iunit
c$$$      stop 'in clossc error in closing file'
c$$$      end
c$$$
c$$$
c$$$      subroutine doad_hgsc(a,g,NDcoef,NDlay,il,alfbet,lmax,eps)
c$$$************************************************************************
c$$$*** calculates the expansion coefficients in Legendre polynomials of ***
c$$$*** a Henyey-Greenstein scattering function (see Van de Hulst, 1980  ***
c$$$*** page 308)                                                        ***
c$$$*** with single scattering albedo (a) and asymmetry factor (g)       ***
c$$$************************************************************************
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension alfbet(4,4,0:NDcoef,NDlay)
c$$$      character*60 s1, s2, skip
c$$$
c$$$      gpowl = 1.d0
c$$$      do 10 l=0,NDcoef
c$$$          do 10 j=1,4
c$$$          do 10 i=1,4
c$$$          alfbet(i,j,l,il) = 0.d0
c$$$   10 continue
c$$$      do 20 l=0,NDcoef
c$$$          alfbet(1,1,l,il) = dble(l+l+1)*gpowl
c$$$          if (alfbet(1,1,l,il) .lt. eps*1.d-1) goto 25
c$$$          gpowl = gpowl*g
c$$$   20 continue
c$$$      print *,' warning, the desired accurary of the expansion ',
c$$$     +  'coefficients for '
c$$$      print *,' Henyey-Greenstein scattering is not reached '
c$$$   25 lmax = l
c$$$************************************************************************
c$$$* prepare for writing to output file                                   *
c$$$************************************************************************
c$$$      s1 = 'Henyey-Greenstein scattering'
c$$$      s2 = 'asymmetry factor g=<cos th>'
c$$$      skip = ' '
c$$$      call doad_writsc(7,alfbet,NDcoef,NDlay,il,lmax,a,g,
c$$$     +     0.d0,0.d0,0.d0,0.d0,skip,s1,s2,skip,skip,skip)
c$$$      do 30 l=0,lmax
c$$$          alfbet(1,1,l,il) = a*alfbet(1,1,l,il)
c$$$   30 continue
c$$$      return
c$$$      end
c$$$
c$$$      subroutine doad_tthgsc(a,g1,frac,g2,NDcoef,NDlay,il,alfbet,
c$$$     +                       lmax,eps)
c$$$************************************************************************
c$$$*** calculates the expansion coefficients in Legendre polynomials of ***
c$$$*** a two terms Henyey-Greenstein scattering function                ***
c$$$************************************************************************
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension alfbet(4,4,0:NDcoef,NDlay)
c$$$      character*60 s1, s2, s3, s4, skip
c$$$
c$$$      do 10 l=0,NDcoef
c$$$      do 10 j=1,4
c$$$      do 10 i=1,4
c$$$        alfbet(i,j,l,il) = 0.d0
c$$$   10 continue
c$$$      g1powl = 1.d0
c$$$      g2powl = 1.d0
c$$$      do 20 l=0,NDcoef
c$$$        alfbet(1,1,l,il) = dble(l+l+1)*(frac*g1powl+(1.d0-frac)*g2powl)
c$$$        if (alfbet(1,1,l,il) .lt. eps*1.d-1) goto 25
c$$$        g1powl = g1powl*g1
c$$$        g2powl = g2powl*g2
c$$$   20 continue
c$$$      print *,' warning, the desired accurary of the expansion ',
c$$$     +  'coefficients for '
c$$$      print *,' Henyey-Greenstein scattering is not reached '
c$$$   25 lmax = l
c$$$************************************************************************
c$$$* prepare for writing to output file                                   *
c$$$************************************************************************
c$$$      s1 = 'Two Terms Henyey-Greenstein scattering'
c$$$      write(s2,'(a,f6.4)') 'asymmetry factor g_1 = ',g1
c$$$      write(s3,'(a,f6.4)') 'asymmetry factor g_2 = ',g2
c$$$      write(s4,'(a,f6.4)') 'fraction f = ',frac
c$$$      skip = ' '
c$$$      call doad_writsc(7,alfbet,NDcoef,NDlay,il,lmax,a,0.d0,
c$$$     +     0.d0,0.d0,0.d0,0.d0,skip,s1,s2,s3,s4,skip)
c$$$      do 30 l=0,lmax
c$$$          alfbet(1,1,l,il) = a*alfbet(1,1,l,il)
c$$$   30 continue
c$$$      return
c$$$      end
c$$$
c$$$
c$$$
c$$$      subroutine doad_isotrp(a,NDcoef,NDlay,il,alfbet,lmax)
c$$$************************************************************************
c$$$*** calculate the expansion coefficients for isotropic scattering    ***
c$$$*** with single scattering albedo (a)                                ***
c$$$************************************************************************
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension alfbet(4,4,0:NDcoef,NDlay)
c$$$      character*60 s1, skip
c$$$
c$$$      do 10 l=0,NDcoef
c$$$      do 10 j=1,4
c$$$      do 10 i=1,4
c$$$          alfbet(i,j,l,il) = 0.d0
c$$$   10 continue
c$$$      alfbet(1,1,0,il) = 1.d0
c$$$      lmax = 0
c$$$************************************************************************
c$$$* prepare for writing to output file                                   *
c$$$************************************************************************
c$$$      s1 = 'isotropic scattering'
c$$$      skip = ' '
c$$$      call doad_writsc(7,alfbet,NDcoef,NDlay,il,lmax,a,0.d0,
c$$$     +     0.d0,0.d0,0.d0,0.d0,skip,s1,skip,skip,skip,skip)
c$$$      alfbet(1,1,0,il) = a*alfbet(1,1,0,il)
c$$$      return
c$$$      end
c$$$
c$$$
c$$$      subroutine doad_molsc(a,dpf,alfbet,NDcoef,NDlay,il,lmax)
c$$$************************************************************************
c$$$*** molsc computes the expansion coefficients in generalized spheri- ***
c$$$*** cal functions of the scattering matrix for molecular scattering  ***
c$$$*** including depolarization (dpf) and absorption (a)                ***
c$$$*** (see Van de Hulst, 1980 page 532-534)                            ***
c$$$************************************************************************
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension alfbet(4,4,0:NDcoef,NDlay)
c$$$      character*60 s1, s2, skip
c$$$
c$$$      c = 2.d0*(1.d0-dpf)/(2.d0+dpf)
c$$$      d = 2.d0*(1.d0-2.d0*dpf)/(2.d0+dpf)
c$$$      do 10 l=0,NDcoef
c$$$      do 10 j=1,4
c$$$      do 10 i=1,4
c$$$          alfbet(i,j,l,il) = 0.d0
c$$$   10 continue
c$$$************************************************************************
c$$$* fill alfbet according to de Rooij (1985) table 2.1                   *
c$$$************************************************************************
c$$$      alfbet(1,1,0,il) = 1.d0
c$$$      alfbet(4,4,1,il) = 1.5d0*d
c$$$      alfbet(1,1,2,il) = 0.5d0*c
c$$$      alfbet(2,2,2,il) = 3.d0*c
c$$$      alfbet(2,1,2,il) = dsqrt(1.5d0)*c
c$$$      alfbet(1,2,2,il) = alfbet(2,1,2,il)
c$$$      lmax = 2
c$$$      if (lmax .gt. NDcoef) lmax = NDcoef
c$$$************************************************************************
c$$$* prepare for writing to output file                                   *
c$$$************************************************************************
c$$$      s1 = 'Molecular scattering including depolarization'//
c$$$     +     ' and absorption'
c$$$      write(s2,'(a,f6.4)') 'depolarization factor (dpf) = ',dpf
c$$$      skip = ' '
c$$$      call doad_writsc(7,alfbet,NDcoef,NDlay,il,lmax,a,0.d0,0.d0,0.d0,
c$$$     +    0.d0,0.d0,skip,s1,s2,skip,skip,skip)
c$$$      do 20 l=0,lmax
c$$$      do 20 j=1,4
c$$$      do 20 i=1,4
c$$$          alfbet(i,j,l,il) = a*alfbet(i,j,l,il)
c$$$   20 continue
c$$$      return
c$$$      end
c$$$
c$$$
c$$$      subroutine doad_opensc(fname,iunit,newold)
c$$$************************************************************************
c$$$*  On entry :                                                          *
c$$$*      fname      name of the file to be opened                        *
c$$$*      iunit      number of the unit to be assigned to the file        *
c$$$*      newold     file status 'new' or 'old'                           *
c$$$*  On exit :                                                           *
c$$$*      The file fname is opened and assigned unit number iunit.        *
c$$$*                                                 VLD January 9, 1989  *
c$$$************************************************************************
c$$$      character*255 fname
c$$$      character*3 newold
c$$$
c$$$      open(unit=iunit,status=newold,err=999,file=trim(adjustl(fname)))
c$$$      return
c$$$  999 print *,' opensc: error in opening file ',trim(fname)
c$$$      print *,'         unit=',iunit,' requested status=',newold
c$$$      stop 'in opensc error in opening file'
c$$$      end
c$$$
c$$$
c$$$      subroutine doad_swrit(susum,sdsum,nmat,nlirf,nlayer,
c$$$     +  nmug,nmuext,xmu,tau,b,sfluxu,sfluxd,
c$$$     +  NDmu,NDsup,NDirf,NDlay)
c$$$************************************************************************
c$$$cmv                 for internal sources 
c$$$*** write arrays usum and dsum on unit 3                             ***
c$$$*** usum and dsum contain the upward and downward radiation resp.    ***
c$$$************************************************************************
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension susum(NDsup,0:NDirf+1,NDlay+1)
c$$$      dimension sdsum(NDsup,0:NDirf+1,NDlay+1)
c$$$      dimension nlirf(0:NDlay+1), xmu(NDmu)
c$$$      dimension tau(0:NDirf+1,0:NDlay+1), b(0:NDlay)
c$$$      dimension sfluxu(0:NDirf,NDlay+1)
c$$$      dimension sfluxd(0:NDirf,NDlay+1)
c$$$
c$$$      nmutot = nmug+nmuext
c$$$cfh      nsup = nmat*nmutot
c$$$      write(3,5) tau(0,1), nlayer
c$$$    5 format(1h1,'The Reflection and Transmission of the total ',
c$$$     +  'atmosphere with sources',/,' with optical thickness ',f4.1,
c$$$     +  ' and consisting of ',i2,' homogeneous layer(s)')
c$$$cfh       write(3,15)
c$$$   15 format(1h0,'The Reflection function')
c$$$cfh       write(3,25) 0.d0
c$$$   25 format(1h0,22x,'phi-phi0 = ',f10.3,/,/,t38,'I             ',
c$$$     +  'Q             U             V')
c$$$      do 10 i=1,nmutot
c$$$        x = xmu(i)
c$$$cfh        write(3,35) 0.d0,x,(susum((i+(k-1)*nmutot),0,nlayer+1)
c$$$cfh     +    ,k=1,nmat)
c$$$   35   format(1h ,'mu0 =',f8.5,3x,'mu =',f8.5,4f14.8)
c$$$   10 continue
c$$$cfh      write(3,45)
c$$$   45 format(1h0,'The Transmission function')
c$$$cfh      write(3,25) 0.d0
c$$$      do 20 i=1,nmutot
c$$$        x = xmu(i)
c$$$cfh        write(3,35) 0.d0,x,(sdsum((i+(k-1)*nmutot),0,1)
c$$$cfh     +    ,k=1,nmat)
c$$$   20 continue
c$$$      do 40 il=1,nlayer
c$$$cfh        write(3,55) il,b(il),nlirf(il)
c$$$   55   format(1h1,/,/,/,' Homogeneous layer number ',i2,
c$$$     +    ' with optical thickness ',f4.1,/,
c$$$     +    ' consisting of ',i2,' internal field level(s)')
c$$$        do 40 kirf=0,nlirf(il)
c$$$          if (kirf .eq. 0) then
c$$$cfh            write(3,64) il,tau(0,il)
c$$$   64       format(1h1,' U field at interface layer ',i2,
c$$$     +        ' at optical depth ',f8.5)
c$$$          else
c$$$cfh            write(3,65) kirf,il,tau(kirf,il)
c$$$   65       format(1h1,' U field at internal level ',i2,' of layer ',
c$$$     +        i2,' at optical depth ',f8.5)
c$$$          endif
c$$$cfh          write(3,25) 0.d0
c$$$          do 30 i=1,nmutot
c$$$            x = xmu(i)
c$$$cfh            write(3,35) 0.d0,x,(susum((i+(k-1)*nmutot),kirf,il)
c$$$cfh     +        ,k=1,nmat)
c$$$   30     continue
c$$$          if (kirf .eq. 0) then
c$$$cfh            write(3,74) il,tau(0,il)
c$$$   74       format(1h0,' D field at interface layer ',i2,
c$$$     +        ' at optical depth ',f8.5)
c$$$          else
c$$$cfh            write(3,75) kirf,il,tau(kirf,il)
c$$$   75       format(1h0,' D field at internal level ',i2,' of layer ',
c$$$     +        i2,' at optical depth ',f8.5)
c$$$          endif
c$$$cfh          write(3,25) 0.d0
c$$$          do 40 i=1,nmutot
c$$$            x = xmu(i)
c$$$cfh            write(3,35) 0.d0,x,(sdsum((i+(k-1)*nmutot),kirf,il)
c$$$cfh     +        ,k=1,nmat)
c$$$   40 continue
c$$$cfh      write(3,85) 
c$$$   85 format(1h1,' U field at the top')
c$$$cfh      write(3,25) 0.d0
c$$$      do 50 i=1,nmutot
c$$$        x = xmu(i)
c$$$cfh        write(3,35) 0.d0,x,(susum((i+(k-1)*nmutot),0,nlayer+1)
c$$$cfh     +    ,k=1,nmat)
c$$$   50 continue
c$$$      write(3,95) 
c$$$   95 format(1h0,' D field at the top')
c$$$cfh      write(3,25) 0.d0
c$$$      do 60 i=1,nmutot
c$$$        x = xmu(i)
c$$$cfh        write(3,35) 0.d0,x,(sdsum((i+(k-1)*nmutot),0,nlayer+1)
c$$$cfh     +    ,k=1,nmat)
c$$$   60 continue
c$$$************************************************************************
c$$$* writing the fluxes to unit 2                                         *
c$$$************************************************************************
c$$$cfh      write(3,105) 0.d0
c$$$  105 format(1h1,' The internal fluxes for sources',f8.6)
c$$$cfh      write(3,115)
c$$$  115 format(1h0,'     tau           flux u        flux d ',
c$$$     +  '     exp(-t/mu0)   nettod  flux')
c$$$      do 100 il=1,nlayer+1
c$$$      do 100 k=0,nlirf(il)
c$$$        xnetf = sfluxd(k,il)-sfluxu(k,il)
c$$$cfh        write(3,125) tau(k,il),sfluxu(k,il),sfluxd(k,il),
c$$$cfh     +    0.d0,xnetf
c$$$  125   format(5f14.8)
c$$$  100 continue
c$$$      return
c$$$      end
c$$$
c$$$
c$$$      subroutine doad_writ(usum,dsum,nmat,nphout,nlirf,nlayer,
c$$$     +  nmug,nmuext,nmu0,xmu,imu0,phout,tau,b,fluxu,fluxd,etmu0,
c$$$     +  NDmu,NDsup,NDmu0,NDphi0,NDirf,NDlay)
c$$$************************************************************************
c$$$*** write arrays usum and dsum on unit 2                             ***
c$$$*** usum and dsum contain the upward and downward radiation resp.    ***
c$$$************************************************************************
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension usum(NDsup,NDmu0,NDphi0,0:NDirf+1,NDlay+1)
c$$$      dimension dsum(NDsup,NDmu0,NDphi0,0:NDirf+1,NDlay+1)
c$$$      dimension nlirf(0:NDlay+1), xmu(NDmu), imu0(NDmu0), phout(NDphi0)
c$$$      dimension tau(0:NDirf+1,0:NDlay+1), b(0:NDlay)
c$$$      dimension fluxu(NDmu0,0:NDirf,NDlay+1)
c$$$      dimension fluxd(NDmu0,0:NDirf,NDlay+1)
c$$$      dimension etmu0(NDmu0,0:NDirf,NDlay+1)
c$$$
c$$$      nmutot = nmug+nmuext
c$$$      nsup = nmat*nmutot
c$$$      do 80 ip=1,nphout
c$$$        write(2,5) tau(0,1), nlayer
c$$$    5   format(1h1,'The Reflection and Transmission of the total ',
c$$$     +    'atmosphere',/,' with optical thickness ',f4.1,
c$$$     +    ' and consisting of ',i2,' homogeneous layer(s)')
c$$$        write(2,15)
c$$$   15   format(1h0,'The Reflection function')
c$$$        write(2,25) phout(ip)
c$$$   25   format(1h0,22x,'phi-phi0 = ',f10.3,/,/,t38,'I             ',
c$$$     +    'Q             U             V')
c$$$        do 10 j=1,nmu0
c$$$          ji = imu0(j)+nmug
c$$$          xj = xmu(ji)
c$$$          do 10 i=1,nmutot
c$$$            x = xmu(i)
c$$$            write(2,35) xj,x,(usum((i+(k-1)*nmutot),j,ip,0,nlayer+1)
c$$$     +        ,k=1,nmat)
c$$$   35       format(1h ,'mu0 =',f8.5,3x,'mu =',f8.5,4f14.8)
c$$$   10   continue
c$$$        write(2,45)
c$$$   45   format(1h0,'The Transmission function')
c$$$        write(2,25) phout(ip)
c$$$        do 20 j=1,nmu0
c$$$          ji = imu0(j)+nmug
c$$$          xj = xmu(ji)
c$$$          do 20 i=1,nmutot
c$$$            x = xmu(i)
c$$$            write(2,35) xj,x,(dsum((i+(k-1)*nmutot),j,ip,0,1)
c$$$     +        ,k=1,nmat)
c$$$   20   continue
c$$$        do 40 il=1,nlayer
c$$$          write(2,55) il,b(il),nlirf(il)
c$$$   55     format(1h1,/,/,/,' Homogeneous layer number ',i2,
c$$$     +      ' with optical thickness ',f4.1,/,
c$$$     +      ' consisting of ',i2,' internal field level(s)')
c$$$          do 40 kirf=0,nlirf(il)
c$$$            if (kirf .eq. 0) then
c$$$              write(2,64) il,tau(0,il)
c$$$   64         format(1h1,' U field at interface layer ',i2,
c$$$     +          ' at optical depth ',f8.5)
c$$$            else
c$$$              write(2,65) kirf,il,tau(kirf,il)
c$$$   65         format(1h1,' U field at internal level ',i2,' of layer ',
c$$$     +          i2,' at optical depth ',f8.5)
c$$$            endif
c$$$            write(2,25) phout(ip)
c$$$            do 30 j=1,nmu0
c$$$              ji = imu0(j)+nmug
c$$$              xj = xmu(ji)
c$$$              do 30 i=1,nmutot
c$$$                x = xmu(i)
c$$$                write(2,35) xj,x,(usum((i+(k-1)*nmutot),j,ip,kirf,il)
c$$$     +            ,k=1,nmat)
c$$$   30       continue
c$$$            if (kirf .eq. 0) then
c$$$              write(2,74) il,tau(0,il)
c$$$   74         format(1h0,' D field at interface layer ',i2,
c$$$     +          ' at optical depth ',f8.5)
c$$$            else
c$$$              write(2,75) kirf,il,tau(kirf,il)
c$$$   75         format(1h0,' D field at internal level ',i2,' of layer ',
c$$$     +          i2,' at optical depth ',f8.5)
c$$$            endif
c$$$            write(2,25) phout(ip)
c$$$            do 40 j=1,nmu0
c$$$              ji = imu0(j)+nmug
c$$$              xj = xmu(ji)
c$$$              do 40 i=1,nmutot
c$$$                x = xmu(i)
c$$$                write(2,35) xj,x,(dsum((i+(k-1)*nmutot),j,ip,kirf,il)
c$$$     +            ,k=1,nmat)
c$$$   40   continue
c$$$        write(2,85) 
c$$$   85   format(1h1,' U field at the top')
c$$$        write(2,25) phout(ip)
c$$$        do 50 j=1,nmu0
c$$$          ji = imu0(j)+nmug
c$$$          xj = xmu(ji)
c$$$          do 50 i=1,nmutot
c$$$            x = xmu(i)
c$$$            write(2,35) xj,x,(usum((i+(k-1)*nmutot),j,ip,0,nlayer+1)
c$$$     +        ,k=1,nmat)
c$$$   50   continue
c$$$        write(2,95) 
c$$$   95   format(1h0,' D field at the top')
c$$$        write(2,25) phout(ip)
c$$$        do 60 j=1,nmu0
c$$$          ji = imu0(j)+nmug
c$$$          xj = xmu(ji)
c$$$          do 60 i=1,nmutot
c$$$            x = xmu(i)
c$$$            write(2,35) xj,x,(dsum((i+(k-1)*nmutot),j,ip,0,nlayer+1)
c$$$     +        ,k=1,nmat)
c$$$   60   continue
c$$$   80 continue
c$$$************************************************************************
c$$$* writing the fluxes to unit 2                                         *
c$$$************************************************************************
c$$$        do 100 j=1,nmu0
c$$$          write(2,105) xmu(imu0(j)+nmug)
c$$$  105     format(1h1,' The internal fluxes for mu0 is ',f8.6)
c$$$          write(2,115)
c$$$  115     format(1h0,'      tau           flux u        flux d ',
c$$$     +    '     exp(-t/mu0)   nettod flux')
c$$$          do 100 il=1,nlayer+1
c$$$          do 100 k=0,nlirf(il)
c$$$            xnetf = fluxd(j,k,il)+etmu0(j,k,il)-fluxu(j,k,il)
c$$$            write(2,125) tau(k,il),fluxu(j,k,il),fluxd(j,k,il),
c$$$     +        etmu0(j,k,il),xnetf
c$$$  125       format(5f14.8)
c$$$  100 continue
c$$$      return
c$$$      end
