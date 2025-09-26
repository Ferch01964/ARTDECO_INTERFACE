************************************************************************
*                                                                      *
*             M E E R H O F F      M I E     P R O G R A M             *
*                                                                      *
*                             VERSION 3.0                              *
*                                                                      *
*         Astronomy Department of the Free University Amsterdam        *
*                                                                      *
*  This program calculates the scattering properties of a homogeneous  *
*  sphere or ensemble of spheres with different radii according to     *
*  the Mie theory. Documentation can be found in :                     *
*                                                                      *
*  - H.C. van de Hulst (1957) 'Light Scattering by Small Particles'    *
*         Wiley, New York; also Dover, New York 1981                   *
*  - W.A. de Rooij, C. van der Stap (1984) 'Expansion of Mie Scattering*
*         Matrices in Generalized Spherical Functions'                 *
*         Astron. Astrophys. 131, p. 237                               *
*  - W.J. Wiscombe (1980) Applied Optics 19, p.1505                    *
*  - J.F. de Haan, P.B. Bosma and J.W. Hovenier (1987):                *
*         Astronomy and Astrophysics 183, pages 371-391.               *
*  - V.L. Dolman (1989) 'Meerhoff Mie Program User Guide'              *
*         internal report Astronomy Dept. Free University              *
*                                                                      *
*   modification: bug fix initialisation of P02 array in expand        *
*                     (version 2.2 --> 2.3)                            *
*                     parameter nrunit instead of '88' in code         *
*                     all print * replaced by write(*,*)               *
*   modification: calculate the volume and print it on output,         *
*                     write the geometric cross section and the        *
*                     volume onto the coefficients file 'mie_sc_sil'   *
*                     (version 2.3 --> 3.0)  November, 26 1990         *
*                                                                      *
*                                                 February 19 1990     *
************************************************************************
      subroutine mie3(NDang, delta, cutoff
     +               , wavel, m, nsubr, ngaur, idis
     +               , par1, par2, par3, igaussu, nangle, F, u, w, miec)
      double precision delta, cutoff, wavel, 
     +                 par1, par2, par3, F, u, w, miec, 
     +                 rmin, rmax
      integer nsubr, ngaur, idis, nangle, igaussu
      double complex m
      dimension F(4,NDang)
     +        , miec(10),u(NDang),w(NDang)

************************************************************************
*  Determine the integration bounds for the radius integration         *
************************************************************************
      call rminmax( idis
     +   , par1, par2, par3, cutoff
     +   , rmin, rmax )

************************************************************************
*  Calculate the scattering matrix : this is most of the work !        *
************************************************************************
      call scatmat(NDang, u, w, m, wavel, idis
     +            , igaussu
     +            , nsubr, ngaur, rmin, rmax
     +            , par1, par2, par3
     +            , delta, F, miec, nangle )

      end subroutine mie3


      subroutine rminmax( idis, par1, par2, par3, cutoff
     +                  , rmin, rmax )
************************************************************************
*  Find the integration bounds rmin and rmax for the integration over  *
*  a size distribution. These bounds are chosen such that the size     *
*  distribution falls below the user specified cutoff. It is essential *
*  that the size distribution is normalized such that the integral     *
*  over all r is equal to one !                                        *
*  This is programmed rather clumsy and will in the future be changed  *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps=1.d-10 )
      double precision nwithr
      dimension r(1), nwithr(1)

      if (idis.eq.0) then
         rmin= par1
         rmax= par1
      else
          goto (10,20,30,40,50,60,70,80) idis
          write(*,*) ' rminmax: illegal size distribution index :',idis
          stop 'in rminmax illegal size distribution index'
************************************************************************
   10     sef = 1.D0/dsqrt(par2+3.D0)
          ref = 1.D0/(sef*sef*par2)
          rref= ref
          goto 100

   20     ref = par1
          sef = dsqrt(par2)
          rref= ref
          goto 100

   30     sef = dsqrt(par3)
          ref = dmax1(par1,par2)+sef
          rref= dmin1(par1,par2)
          goto 100

   40     sef = dsqrt(dexp((par2)**2)-1.d0)
          ref = par1*(1.D0+sef*sef)**0.4D0
          rref= ref
          goto 100

   50     ref = par1
!<BIZARRE!!!!!!!!!!!!!!!!!!!!!,23/09/04>
          sef = dsqrt(ref)
!</BIZ>
          rref= ref
          goto 100

   60     rmin= par2
          rmax= par3
          goto 999

   70     ref = par2
          sef = 2.D0*sef
          rref=0.5D0*ref
          goto 100

   80     ref = (par1/(par2*par3))**par3
          sef = 2.D0*ref
          rref= 0.5D0*ref
************************************************************************
*  search for a value of r such that the size distribution
*  is less than the cutoff. Start the search at ref+sef which          *
*  guarantees that such a value will be found on the TAIL of the       *
*  distribution.                                                       *
************************************************************************
  100     r(1) = ref+sef
          r0   = ref
  200          call sizedis( idis, par1, par2, par3, r, 1, nwithr )
          if (nwithr(1) .gt. cutoff) then
              r0   = r(1)
              r(1) = 2.D0*r(1)
              goto 200
          endif
          r1 = r(1)
************************************************************************
*  Now the size distribution assumes the cutoff value somewhere        *
*  between r0 and r1  Use bisection to find the corresponding r        *
************************************************************************
  300     r(1) = 0.5D0*(r0+r1)
          call sizedis( idis, par1, par2, par3, r, 1, nwithr )
          if (nwithr(1) .gt. cutoff) then
              r0 = r(1)
          else
              r1 = r(1)
          endif
          if ((r1-r0) .gt. eps) goto 300
          rmax = 0.5D0*(r0+r1)
************************************************************************
*  Search for a value of r on the low end of the size distribution     *
*  such that the distribution falls below the cutoff. There is no      *
*  guarantee that such a value exists, so use an extra test to see if  *
*  the search comes very near to r = 0                                 *
************************************************************************
          r1 = rref
          r0 = 0.D0
  400     r(1) = 0.5D0*r1
          call sizedis( idis, par1, par2, par3, r, 1, nwithr )
          if (nwithr(1) .gt. cutoff) then
              r1 = r(1)
              if (r1 .gt. eps) goto 400
          else
              r0 = r(1)
          endif
************************************************************************
*  Possibly the size distribution goes through cutoff between r0       *
*  and r1 try to find the exact value of r where this happens by       *
*  bisection.                                                          *
*  In case there is no solution, the algorithm will terminate soon.    *
************************************************************************
  500     r(1) = 0.5D0*(r0+r1)
          call sizedis( idis, par1, par2, par3, r, 1, nwithr )
          if (nwithr(1) .gt. cutoff) then
              r1 = r(1)
          else
              r0 = r(1)
          endif
          if ((r1-r0) .gt. eps) goto 500
          if(r1 .le. eps) then
              rmin = 0.D0
          else
              rmin = 0.5D0*(r0+r1)
          endif
      endif
 
c      write(*,*) ' rminmax: found rmin = ',rmin,' rmax = ',rmax
  999 return
      end


      subroutine sizedis( idis, par1, par2, par3, r, numr, nwithr )
************************************************************************
*  Calculate the size distribution n(r) for the numr radius values     *
*  contained in array r and return the results through the array nwithr*
*  The size distributions are normalized such that the integral over   *
*  all r is equal to one.                                              *
************************************************************************
      implicit double precision (a-h,o-z)
      double precision nwithr, logC, logC1, logC2
      dimension r(numr),nwithr(numr)
*
      pi     = dacos(-1.d0)
      root2p = dsqrt(pi+pi)
      if (idis .eq. 0) return
      goto (10, 20, 30, 40, 50, 60, 70, 80 ) idis
      write(*,*) ' sizedis: illegal index : ',idis
      stop 'illegal index in sizedis'
************************************************************************
*  1 TWO PARAMETER GAMMA with alpha and b given                        *
************************************************************************
   10 alpha = par1
      b     = par2
      alpha1= alpha+1.D0
      logC  = alpha1*dlog(b)-gammln(alpha1)
CL
CL je demande au code d'ecrire ici logC
CL
      write(6,*)'logC=',logC
      do 11 i=1, numr
          nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
   11 continue
      goto 999
************************************************************************
*  2 TWO PARAMETER GAMMA with par1= reff and par2= veff given          *
************************************************************************
   20 alpha = 1.D0/par2 - 3.D0
      b     = 1.D0/(par1*par2)
      alpha1= alpha+1.D0
      logC  = alpha1*dlog(b)-gammln(alpha1)
      do 21 i=1, numr
          nwithr(i) = dexp(logC+alpha*dlog(r(i))-b*r(i))
   21 continue
      goto 999
************************************************************************
*  3 BIMODAL GAMMA with equal mode weights                             *
************************************************************************
   30 alpha = 1.D0/par3 - 3.D0
      b1    = 1.D0/(par1*par3)
      b2    = 1.D0/(par2*par3)
      gamlna= gammln(alpha+1.D0)
      logC1 = (alpha+1.D0)*dlog(b1)-gamlna
      logC2 = (alpha+1.D0)*dlog(b2)-gamlna
      do 31 i=1, numr
          nwithr(i) = 0.5D0*( dexp(logC1+alpha*dlog(r(i))-b1*r(i))
     +                      + dexp(logC2+alpha*dlog(r(i))-b2*r(i)) )
   31 continue
      goto 999
************************************************************************
*  4 LOG NORMAL with rg and sigma given                                *
************************************************************************
   40 flogrg = dlog(par1)
      flogsi = dabs(dlog(par2))
      C      = 1.D0/(root2p*flogsi)
      fac    = -0.5D0/(flogsi*flogsi)
      do 41 i=1, numr
          nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
   41 continue
      goto 999
************************************************************************
*  5 LOG NORMAL with reff and veff given                               *
************************************************************************
   50 rg     = par1/(1.D0+par2)**2.5D0
      flogrg = dlog(rg)
      flogsi = dsqrt(dlog(1.D0+par2))
      C      = 1.D0/(root2p*flogsi)
      fac    = -0.5D0/(flogsi*flogsi)
      do 51 i=1, numr
          nwithr(i) = C * dexp( fac*(dlog(r(i))-flogrg)**2 ) / r(i)
   51 continue
      goto 999
************************************************************************
*  6 POWER LAW                                                         *
************************************************************************
   60 alpha = par1
      rmin  = par2
      rmax  = par3
* DEBUG
*     write(*,*) ' sizedis : power law with alpha = ',alpha
*     write(*,*) '                          rmin  = ',rmin
*     write(*,*) '                          rmax  = ',rmax
* END DEBUG
      if (dabs(alpha+1.D0) .lt. 1.d-10) then
          C = 1.D0/dlog(rmax/rmin)
      else
          alpha1 = alpha-1.d0
          C = alpha1 * rmax**alpha1 / ((rmax/rmin)**alpha1-1.d0)
      endif
      do 61 i=1, numr
          if ((r(i) .lt. rmax) .and. (r(i) .gt. rmin)) then
              nwithr(i) = C*r(i)**(-alpha)
          else
              nwithr(i) = 0.D0
          endif
   61 continue
      goto 999
************************************************************************
*  7 MODIFIED GAMMA with alpha, rc and gamma given                     *
************************************************************************
   70 alpha = par1
      rc    = par2
      gamma = par3
      b     = alpha / (gamma*(rc**gamma))
      aperg = (alpha+1.D0)/gamma
      logC  = dlog(gamma) + aperg*dlog(b) - gammln(aperg)
      do 71 i=1, numr
          nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
   71 continue
      goto 999
************************************************************************
*  8 MODIFIED GAMMA with alpha, b and gamma given                      *
************************************************************************
   80 alpha = par1
      b     = par2
      gamma = par3
      aperg = (alpha+1.D0)/gamma
      logC  = dlog(gamma) + aperg*dlog(b) - gammln(aperg)
      do 81 i=1, numr
          nwithr(i) = dexp( logC + alpha*dlog(r(i)) - b*r(i)**gamma )
   81 continue
************************************************************************
  999 if (numr .le. 1) return
* DEBUG
*     write(*,*) ' sizedis:'
*     write(*,*) '   i             r(i)               n(r(i))'
*     write(*,*) ' ----------------------------------------------------'
*     do 1000 i=1, numr
*         write(*,1001) i,r(i),nwithr(i)
*1000 continue
*1001 format(i4,1pe24.14,1pe22.14)
* END DEBUG
      return
      end
 
      function gammln(xarg)
************************************************************************
*  Return the value of the natural logarithm of the gamma function.    *
*  The argument xarg must be real and positive.                        *
*  This function is documented in :                                    *
*                                                                      *
*  W.H. Press et al. 1986, 'Numerical Recipes' Cambridge Univ. Pr.     *
*  page 157 (ISBN 0-521-30811)                                         *
*                                                                      *
*  When the argument xarg is between zero and one, the relation (6.1.4)*
*  on page 156 of the book by Press is used.                           *
*                                         V.L. Dolman April 18 1989    *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps=1.d-7, one=1.D0, two=2.D0, half=0.5D0
     +         , fpf=5.5D0 )
      dimension cof(6)
      data cof,stp/76.18009173D0,-86.50532033D0, 24.01409822D0
     +   ,-1.231739516D0, 0.120858003D-2, -0.536382D-5, 2.50662827465D0/
      pi = 4.D0*datan(1.D0)
      if (xarg .le. 0.D0) then
        write(*,*) ' gammln: called with negative argument xarg = ',xarg
        stop 'function gammln called with negative value'
      endif
      if (dabs(xarg-one) .lt. eps) then
          write(*,*) ' gammln: argument too close to one for algorithm'
          stop ' in function gammln argument too close to one'
      endif
      if (xarg .ge. one) then
          xx = xarg
      else
          xx = xarg+two
      endif
      x = xx-one
      tmp = x+fpf
      tmp = (x+half)*dlog(tmp)-tmp
      ser = one
      do 11 j=1, 6
          x = x+one
          ser = ser+cof(j)/x
   11 continue
      gtmp = tmp+dlog(stp*ser)
      if  (xarg .gt. one) then
          gammln = gtmp
      else
          pix = pi*(one-xarg)
          gammln = dlog(pix/dsin(pix))-gtmp
      endif
      return
      end

      subroutine scatmat( NDang, u, wth, m, lambda, ndis
     +                  , igaussu
     +                  , nsub, ngauss, rmin, rmax
     +                  , par1, par2, par3
     +                  , delta, F, miec, nangle )
************************************************************************
*  Calculate the scattering matrix of an ensemble of homogenous        *
*  spheres. On entry, the following must be supplied :                 *
*     m            : complex index of refraction                       *
*     lambda       : wavelength                                        *
*     ndis         : index of the size distribution                    *
*     nsub         : number of subintervals for integration over r     *
*     ngauss       : number of Gauss points used per subinterval       *
*     rmin         : lower bound for integration over r                *
*     rmax         : upper bound for integration over r                *
*     par1,2,3     : parameters of the size distribution               *
*     delta        : cutoff used in truncation of the Mie sum          *
*     igaussu      : When =1 u contains Gauss points,                  *
*                    when =0 u contains cosines of equidistant         *
*                     angles between thmin and thmax with step.        *
*  On exit, the following results are returned :                       *
*     u            : cosines of scattering angles                      *
*     wth          : Gaussian weights associated with u                *
*     F            : scattering matrix for all cosines in u            *
*     miec         : array containing cross sections etc.              *
*     nangle       : the number of scattering angles                   *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( ndn=4000, NDr=1000 )
      double complex   m, ci, d, Splusf, Sminf, cSplusf
      double complex   cSminf, Splusb, Sminb, cSplusb, cSminb
      double complex   an(NDn), bn(NDn)
      double precision lambda, nwithr, miec, numpar, thmin, thmax, step
      integer     igaussu
      dimension   u(NDang), wth(NDang), F(4,NDang)
      dimension   pi(NDn), tau(NDn), fi(0:NDn), chi(0:NDn)
      dimension   D(NDn), miec(10)
      dimension   r(NDr), w(NDr), nwithr(NDr)
      dimension   facf(NDn), facb(NDn)
      logical     symth
************************************************************************
*  Initialization                                                      *
************************************************************************
      thmin = 0.
      thmax = 180.
      step  = (thmax-thmin)/(nangle-1)
      do 10 j=1,NDang
          do 10 k=1,4
   10         F(k,j)=0.D0
      Csca  = 0.D0
      Cext  = 0.D0
      numpar= 0.D0
      G     = 0.D0
      reff  = 0.D0
      nfou  = 0
      fac90 = 1.D0
      ci    = dcmplx(0.D0,1.D0)
      call tstsym( igaussu, thmin, thmax, step, symth )
************************************************************************
*  Constants                                                           *
************************************************************************
      pie   = dacos(-1.d0)
      radfac= pie/180.D0
      rtox  = 2.D0*pie/lambda
      fakt  = lambda*lambda/(2.D0*pie)
* nfac is the number of precomputed factors (2n+1)/(n*(n+1))
      nfac  = 0
************************************************************************
*  distinguish between distribution or not                             *
************************************************************************
      if (ndis.eq.0) then
         w(1)     = 1.D0
         r(1)     = rmin
         nwithr(1)= 1.D0
         nsub     = 1
         ngauss   = 1
         dr       = 0.D0
      else
         dr = (rmax-rmin)/dble(nsub)
*        call gausspt( ngauss, ngauss, (rmax-dr) , rmax ,r ,w )
         call mie3_gauleg(  ngauss, ngauss, (rmax-dr) , rmax ,r ,w )
         call sizedis( ndis, par1, par2, par3, r, ngauss, nwithr )
      endif
************************************************************************
*  Start integration over radius r with largest radius !               *
************************************************************************
      do 60 l=nsub,1,-1
      do 50 i=ngauss,1,-1
*
      sw   = nwithr(i)*w(i)
      x    = rtox*r(i)
      nmax = INT(x + 4.05D0*x**(1.D0/3.D0) + 20)
      nfi  = nmax+60
      zabs = x*cdabs(m)
      nD   = INT(zabs + 4.05D0*zabs**(1.D0/3.D0) + 70)
      if ((nD.gt.NDn) .or. (nfi.gt.NDn)) then
          write(*,*) ' scatmat: estimated number of Mie-terms:',nD
          write(*,*) '          for particle sizeparameter   :',x
          write(*,*) '          maximum NDn is only          : ',NDn
          stop 'in scatmat no room for all Mie terms'
      endif
      call fichid( m, x, nfi, nmax, nD, fi, chi, D )
      call anbn( m, x, nmax, fi, chi, D, an, bn )
************************************************************************
*  Precompute the factor (2n+1)/(n*(n+1)) needed in Mie sum over n     *
************************************************************************
      if (nmax .gt. nfac) then
          do 26 n=nfac+1, nmax
              facf(n) = dble(2*n+1)/dble(n*(n+1))
              facb(n) = facf(n)
              if (mod(n,2) .eq. 1) facb(n) = -facb(n)
   26     continue
          nfac = nmax
      endif
************************************************************************
*  Calculate extinction and scattering cross section                   *
*  Use the convergence criterion to determine the number of terms that *
*  will later be used in the mie sum for the scattering matrix itself  *
************************************************************************
      Cextsum = 0.D0
      Cscasum = 0.D0
      nstop = nmax
      do 52 n=1, nmax
          aux = (2.D0*dble(n)+1.D0) *
     +             dabs(dble(an(n)*conjg(an(n)) + bn(n)*conjg(bn(n))))
          Cscasum = Cscasum + aux
          Cextsum = Cextsum + (2.D0*n+1.D0)*dble(an(n)+bn(n))
          if (aux .lt. delta) then
              nstop = n
              goto 53
          endif
   52 continue
   53 nfou = nstop
      if (nfou .ge. nmax) then
          write(*,*) ' WARNING from scatmat : Mie sum not converged for'
     +           ,' scattering cross section'
          write(*,*) '   radius r = ',r(i),' sizeparameter x = ',x
     +           ,' sizedistribution nr. ',ndis
          write(*,*) '   Re(m) = ',dble(m),' Im(m) = ',dimag(m)
          write(*,*) '   a priori estimate of number of Mie terms:',nmax
          write(*,*) '   term ',nmax,' for Csca was ',aux
          write(*,*) '   should have been less than ',delta
          write(*,*) '   the a priori estimate will be used'
      endif
************************************************************************
*  Only for the first run through the loop set points in u= dcos(th)   *
************************************************************************
      if ((l.eq.nsub) .and. (i.eq.ngauss)) then
************************************************************************
*  In case of expansion in GSF : set Gauss points for dcos(th)         *
************************************************************************
          if (igaussu .eq. 1) then
************************************************************************
*  Ensure accurate integrations: add two terms: nangle = 2*nfou+2      *
*  One should be sufficient, but total should be even !                *
************************************************************************
c              nangle = 2*nfou+2
*             if (nangle .lt. 2*nfou+2) then
c              write(*,*) 'nstop   =, nmax ', nstop, nmax
              if ((nangle-2) .lt. 2*nfou+2) then
                write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write(*,*) ' ERROR (in mie3)'
*                write(*,*) ' scatmat: nangle should be at least 2*nfou+2
*     +          for a good accuracy'
                write(*,*) ' scatmat: nangle-2 should be at least 2*nfou
     +          +2 for a good accuracy'
*                write(*,*) 'nangle   = ', nangle
                write(*,*) 'nangle   = ', nangle - 2
                write(*,*) '2*nfou+2 = ', 2*nfou+2
                write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                stop
             endif
              if (nangle .gt. NDang) then
                 write(*,*) ' scatmat: need too many integration angles'
     +                  ,' nangle=',nangle
                 write(*,*) '       maximum array size NDang = ',NDang
                 stop 'too many integration angles'
              endif
*             call gausspt(nangle,nangle,-1.d0,1.D0,u,wth)
*              call mie3_gauleg(nangle,nangle,-1.d0,1.D0,u,wth)
                call mie3_gauleg_boundin(nangle,nangle,-1.d0,1.D0,u,wth)
          else
************************************************************************
*  In case no expansion in GSF is desired : set u= dcos(th) for        *
*  for equidistant angles between thmin and thmax.                     *
************************************************************************
c              nangle = idnint((thmax-thmin)/step) + 1
              if (nangle .gt. NDang) then
                 write(*,*) ' scatmat : need too many scattering angles'
     +                  ,' nangle=',nangle
                 write(*,*) '       maximum array size NDang = ',NDang
                 stop 'too many scattering angles'
              endif
              wfac = -32768.0 ! <MC> 2.d0/dble(nangle)
              do 260 iang=1, nangle
                  th = thmin + dble(iang-1)*step
                  u(nangle+1-iang) = dcos( radfac*th )
                  wth(iang) = wfac
  260         continue
          endif
      endif
************************************************************************
*  Integration for normalization of size distibution, geometrical      *
*  cross section and effective radius                                  *
************************************************************************
      numpar = numpar+sw
      G      = G     +sw*r(i)*r(i)
      reff   = reff  +sw*r(i)*r(i)*r(i)
      if (symth) then
************************************************************************
*  Start loop over scattering angles, do only half and use symmetry    *
*  between forward and backward scattering angles                      *
*  The factor fac90 will later be used to correct for the fact that    *
*  for a symmetrical set of scattering angles with an odd number of    *
*  angles the scattering matrix is a factor 2 too big at 90 degrees    *
*  due to the way we programmed the symmetry relations                 *
************************************************************************
        if (mod(nangle,2) .eq. 1) then
            nhalf = (nangle+1)/2
            fac90 = 0.5D0
        else
            nhalf = nangle/2
        endif
*
        do 40 j=1, nhalf
            call pitau( u(j), nmax, pi, tau )
            Splusf = dcmplx(0.D0,0.D0)
            Sminf  = dcmplx(0.D0,0.D0)
            Splusb = dcmplx(0.D0,0.D0)
            Sminb  = dcmplx(0.D0,0.D0)
*  THIS IS THE INNERMOST LOOP !! (Mie sum)
*  can be programmed more efficiently by taking the facf multiplication
*  outside the angle loop over index j 
            do 20 n=1,nfou
                Splusf = Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
                Sminf  = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
                Splusb = Splusb + facb(n)*(an(n)+bn(n)) * (pi(n)-tau(n))
                Sminb  = Sminb  + facb(n)*(an(n)-bn(n)) * (pi(n)+tau(n))
   20       continue
            cSplusf = conjg(Splusf)
            cSminf  = conjg(Sminf )
            cSplusb = conjg(Splusb)
            cSminb  = conjg(Sminb )
            k = nangle-j+1 
*  the forward scattering elements
            F(1,j) = F(1,j) + 
     + DBLE(sw*(Splusf*cSplusf + Sminf *cSminf))
            F(2,j) = F(2,j) - 
     + DBLE(sw*(Sminf *cSplusf + Splusf*cSminf))
            F(3,j) = F(3,j) + 
     + DBLE(sw*(Splusf*cSplusf - Sminf *cSminf))
            F(4,j) = F(4,j) + 
     + DBLE(ci*sw*(Sminf *cSplusf - Splusf*cSminf))
*  the backward scattering elements
            F(1,k) = F(1,k) + 
     + DBLE(sw*(Splusb*cSplusb + Sminb *cSminb))
            F(2,k) = F(2,k) - 
     + DBLE(sw*(Sminb *cSplusb + Splusb*cSminb))
            F(3,k) = F(3,k) + 
     + DBLE(sw*(Splusb*cSplusb - Sminb *cSminb))
            F(4,k) = F(4,k) + 
     + DBLE(ci*sw*(Sminb *cSplusb - Splusb*cSminb))
   40   continue
      else
************************************************************************
*  Start loop over scattering angles, do all angles                    *
************************************************************************
        do 400 j=1, nangle
            call pitau( u(j), nmax, pi, tau )
            Splusf = dcmplx(0.D0,0.D0)
            Sminf  = dcmplx(0.D0,0.D0)
*  THIS IS THE INNERMOST LOOP !! (Mie sum)
            do 200 n=1,nfou
                Splusf = Splusf + facf(n)*(an(n)+bn(n)) * (pi(n)+tau(n))
                Sminf  = Sminf  + facf(n)*(an(n)-bn(n)) * (pi(n)-tau(n))
  200       continue
            cSplusf = conjg(Splusf)
            cSminf  = conjg(Sminf )
            k = nangle-j+1
*  the forward scattering elements
            F(1,j) = F(1,j) +      
     + DBLE(sw*(Splusf*cSplusf + Sminf *cSminf))
            F(2,j) = F(2,j) -      
     + DBLE(sw*(Sminf *cSplusf + Splusf*cSminf))
            F(3,j) = F(3,j) +      
     + DBLE(sw*(Splusf*cSplusf - Sminf *cSminf))
            F(4,j) = F(4,j) + 
     + DBLE(ci*sw*(Sminf *cSplusf - Splusf*cSminf))
  400   continue
      endif
************************************************************************
*  Integration for cross sections, shift radius to next subinterval    *
************************************************************************
      Csca = Csca + sw*Cscasum
      Cext = Cext + sw*Cextsum
      r(i) = r(i) - dr
   50 continue
      if (l .ne. 1)
     +         call sizedis( ndis, par1, par2, par3, r, ngauss, nwithr )
   60 continue
************************************************************************
*  End of integration over size distribution                           *
*  Some final corrections :                                            *
************************************************************************
      do 70 j=1,nangle
          do 70 k=1,4
   70         F(k,j)= F(k,j)/(2.D0*Csca)
      if (symth) then
          do 80 k=1,4
   80         F(k,nhalf) = fac90*F(k,nhalf)
      endif

      G     = pie*G            ! geometrical cross-section microns^2
      Csca  = Csca*fakt
      Cext  = Cext*fakt
      Qsca  = Csca/G
      Qext  = Cext/G
      albedo= Csca/Cext
      volume= (4.d0/3.d0)*pie*reff
      reff  = pie*reff/G
      xeff  = rtox*reff

      miec(1) = Csca
      miec(2) = Cext
      miec(3) = Qsca
      miec(4) = Qext
      miec(5) = albedo
      miec(6) = G 
      miec(7) = reff
      miec(8) = xeff
      miec(9) = numpar
      miec(10)= volume

      return
      end


      subroutine tstsym( igaussu, thmin, thmax, step, symth )
************************************************************************
*  Test if the set of theta points is symmetrical around 90 degrees    *
*  and return this information through logical symth                   *
*  In case of development in GSF we have a symmetrical Gauss set !!    *
************************************************************************
      implicit double precision(a-h,o-z)
      parameter( eps=1.d-6, heps=0.5d0*eps )
      integer igaussu
      logical symth
      symth = .false.
      if (igaussu .eq. 1) then
          symth = .true.
      else
          if ( (dabs( 180.d0 - thmin - thmax ) .lt. eps)  .and.
     +         (dmod( thmax-thmin+heps, step ) .lt. eps) )
     +              symth = .true.
      endif
* DEBUG
c      if (symth) then
c          write(*,*) ' tstsym: theta points symmetrical'
c      else
c          write(*,*) ' tstsym: theta points NOT symmetrical !!'
c      endif
* END DEBUG
      return
      end



      subroutine fichid( m, x, nchi, nmax, nD, psi, chi, D )
************************************************************************
*  Calculate functions psi(x)  chi(x) and D(z) where z = mx.           *
*  On entry, the following should be supplied :                        *
*      m      : complex index of refraction                            *
*      x      : sizeparameter                                          *
*      nchi   : starting index for backwards recurrence of chi         *
*      nmax   : number of chi, psi and D that must be available        *
*      nD     : starting index for backwards recurrence of D           *
*  On exit, the desired functions are returned through psi, chi and D  *
************************************************************************
      implicit double precision (a-h,o-z)
      double complex D,m,z, perz, zn1
      dimension psi(0:nchi), chi(0:nmax+1), D(nd)
*
      z = m*x
      perz = 1.D0/z
      perx = 1.D0/x
      sinx = dsin(x)
      cosx = dcos(x)
************************************************************************
*  (mis-) use the array psi to calculate the functions rn(x)
*  De Rooij and van der Stap Eq. (A6)
************************************************************************
      do 10 n=nchi-1, 0, -1
          psi(n) = 1.D0 / (dble(2*n+1)/x - psi(n+1))
   10 continue
************************************************************************
*  Calculate functions D(z) by backward recurrence
*  De Rooij and van der Stap Eq. (A11)
************************************************************************
      D(nD) = dcmplx(0.D0,0.D0)
      do 20 n=nD - 1, 1, -1
          zn1 = dble(n+1)*perz
          D(n) = zn1 - 1.D0/(D(n+1)+zn1)
   20 continue
************************************************************************
*  De Rooij and van der Stap Eqs. (A3) and (A1)
*  De Rooij and van der Stap Eq. (A5) where psi(n) was equal to r(n)
*  and Eq. (A2)
************************************************************************
      psi(0) = sinx
      psi1   = psi(0)*perx - cosx
      if (dabs(psi1) .gt. 1.d-4) then
          psi(1) = psi1
          do 30 n=2,nmax
              psi(n) = psi(n)*psi(n-1)
   30     continue
      else
          do 35 n=1,nmax
              psi(n) = psi(n)*psi(n-1)
   35     continue
      endif
************************************************************************
*  De Rooij and van der Stap Eqs. (A4) and (A2)
************************************************************************
      chi(0) = cosx
      chi(1) = chi(0)*perx + sinx
      do 40 n=1, nmax-1
          chi(n+1) = dble(2*n+1)*chi(n)*perx - chi(n-1)
   40 continue
* DEBUG
*     write(*,*) ' fichid: x = ',x
*     write(*,12)
*     do 26 n=0, nchi
*         write(*,11) n,psi(n),chi(n),D(n)
*  26 continue
*  11 format(i4,4e24.14)
*  12 format('   n',t20,'psi(n)',t44,'chi(n)',t68
*    +      ,'Re(D(n))',t92,'Im(D(n))',/
*    +      ,' ----------------------------------------------------'
*    +      ,'-----------------------------------------------------')
* END DEBUG
      return
      end


      subroutine anbn( m, x, nmax, psi, chi, d, an, bn )
************************************************************************
*  Calculate the Mie coefficients an and bn.                           *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( ndn=4000 )
      double complex m, zn, znm1, saved, perm
      double complex an(NDn), bn(NDn), D(NDn)
      dimension psi(0:NDn), chi(0:NDn)
      perm = 1.D0/m
      perx = 1.D0/x
      xn   = 0.D0
* DEBUG
*     write(*,*) ' anbn:'
*     write(*,*) '     Re(an)           Im(an)       Re(bn)      Im(bn)'
* END DEBUG
      do 100 n=1, nmax
          zn   = dcmplx( psi(n),   chi(n))
          znm1 = dcmplx( psi(n-1), chi(n-1))
          xn   = dble(n)*perx
          saved = D(n)*perm+xn
          an(n)= (saved*psi(n)-psi(n-1)) / (saved*zn-znm1)
          saved = m*D(n)+xn
          bn(n)= (saved*psi(n)-psi(n-1)) / (saved*zn-znm1)
* DEBUG
*          write(*,*) n,an(n),bn(n)
* END DEBUG
  100 continue
      return
      end

      subroutine pitau(u,nmax,pi,tau)
c     ***********************************************************
c     calculates pi,n(u) and tau,n(u) with upward recursion
c
c     ***********************************************************
      implicit double precision (a-h,o-z)
      dimension pi(nmax),tau(nmax)
c
c       starting values:
      pi(1) = 1.D0
      pi(2) = 3.D0*u
      delta = 3.D0*u*u-1.d0
      tau(1)= u
      tau(2)= 2.D0*delta-1.d0

c       upward recursion:
      do 10 n=2, nmax-1
          pi(n+1) = dble(n+1)/dble(n)*delta + u*pi(n)
          delta   = u*pi(n+1) - pi(n)
          tau(n+1)= dble(n+1)*delta - pi(n)
   10 continue
      return
      end

c      subroutine gausspt(ndim,ngauss,a,b,x,w)
c
c     ***********************************************************
c
c     put the gauss-legendre points of order ndim in array x,
c     the weights in array w. the points are transformed to the
c     interval [a,b]
c
c     ***********************************************************
c
c      implicit double precision (a-h,o-z)
c      dimension x(ndim),w(ndim)
c
c     ***********************************************************
c
c                     find starting values
c
c      gn=0.5D0/dble(ngauss)
c      extra=1.0D0/(.4D0*dble(ngauss)*dble(ngauss)+5.0D0)
c      xz=-gn
c      nt=0
c      nteken=0
c    5 pnm2=1.0D0
c      pnm1= xz
c      do 10 i=2,ngauss
c          pnm1xz= pnm1*xz
c          pn=2.0D0*pnm1xz-pnm2-(pnm1xz-pnm2)/dble(i)
c          pnm2= pnm1
c          pnm1= pn
c   10 continue
c      mteken=1
c      if (pn .le. 0.0D0) mteken=-1
c      if ((mteken+nteken) .eq. 0) then
c          nt=nt+1
c          x(nt)= xz
c      endif
c      nteken=mteken
c      if ((1.0D0-xz) .le. extra) go to 30
c      xz= xz+(1.D0-xz*xz)*gn+extra
c      go to 5
c   30 continue
c
c     ***********************************************************
c
c                determine zero's and weights
c
c      do 60 i=1,nt
c          xz= x(i)
c          delta2=1.D0
c   35         pnm2=1.0D0
c              pnm1= xz
c              pnm1af=1.0D0
c              z=.5D0+1.5D0*xz*xz
c              do 40 k=2,ngauss
c                  pnm1xz= pnm1*xz
c                  pn=2.0D0*pnm1xz-pnm2-(pnm1xz-pnm2)/dble(k)
c                  pnaf= xz*pnm1af+k*pnm1
c                  z= z+(dble(k)+0.5D0)*pn*pn
c                  pnm2= pnm1
c                  pnm1= pn
c                  pnm1af= pnaf
c   40         continue
c              delta1= pn/pnaf
c              xz= xz-delta1
c              if(delta1.lt.0.0D0) delta1=-delta1
c              if((delta1.ge.delta2).and.(delta2.lt.1.d-14)) go to 50
c              delta2= delta1
c          go to 35
c   50     x(i)= xz
c          w(i)=1.0D0/z
c   60 continue
c
c     ***********************************************************
c
c                 transform to the interval [a,b]
c
c      nghalf=ngauss/2
c      ngp1=ngauss+1
c      ntp1=nt+1
c      apb= a+b
c      bmag2=(b-a)/2.0D0
c      do 90 i=1,nghalf
c          x(ngp1-i)= b-bmag2*(1.0D0-x(ntp1-i))
c   90     w(ngp1-i)= bmag2*w(ntp1-i)
c      if (nghalf .ne. nt) then
c          x(nt)= apb/2.0D0
c          w(nt)= w(1)*bmag2
c      endif
c      do 120 i=1,nghalf
c          x(i)= apb-x(ngp1-i)
c          w(i)= w(ngp1-i)
c  120 continue
c      return
c      end


        subroutine mie3_gauleg_boundin(ndimin,ngaussin,a,b,xout,wout)
************************************************************************
*     Given the lower and upper limits of integration a and b, and given *
*     the number of Gauss-Legendre points ngauss, this routine returns   *
*     through array x the abscissas and through array w the weights of   *
*     the Gauss-Legendre quadrature formula. Eps is the desired accuracy *
*     of the abscissas. This routine is documented further in :          *
*     W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)   *
*     page 125 ISBN 0-521-30811-9                                        *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps=1.d-14 )
      double precision x(ndimin-2), w(ndimin-2)
      double precision xout(ndimin), wout(ndimin)
      double precision a,b,xm,xl,z,p1,p2,p3,pp,z1,pi
      ndim   = ndimin - 2
      ngauss = ngaussin - 2
      
      pi=4.D0*datan(1.D0)
      m=(ngauss+1)/2
      xm=0.5D0*(a+b)
      xl=0.5D0*(b-a)
      do 92 i=1,m
*     THIS IS A REALLY CLEVER ESTIMATE :
         z= dcos(pi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
   81     continue
              p1=1.D0
              p2=0.D0
              do 91 j=1,ngauss
                  p3= p2
                  p2= p1
                  p1=((dble(2*j)-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
   91         continue
              pp=ngauss*(z*p1-p2)/(z*z-1.d0)
              z1= z
              z= z1-p1/pp
          if (dabs(z-z1) .gt. eps) goto 81
          x(i)= xm-xl*z
          x(ngauss+1-i)= xm+xl*z
          w(i)=2.D0*xl/((1.D0-z*z)*pp*pp)
          w(ngauss+1-i)= w(i)
   92 continue
      
      xout(1)          = a
      xout(ndimin)     = b
      xout(2:ndimin-1) = x(:)
      wout(1)          = 0.0
      wout(ndimin)     = 0.0
      wout(2:ndimin-1) = w(:)
      
c      do 43 i = 1, ndimin
c         write(*,*) ' mie3_gauleg: ',i,'  x=',xout(i),' w=',wout(i)
c 43      continue

      return
      end



      subroutine mie3_gauleg(ndim,ngauss,a,b,x,w)
************************************************************************
*   Given the lower and upper limits of integration a and b, and given *
*   the number of Gauss-Legendre points ngauss, this routine returns   *
*   through array x the abscissas and through array w the weights of   *
*   the Gauss-Legendre quadrature formula. Eps is the desired accuracy *
*   of the abscissas. This routine is documented further in :          *
*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)   *
*   page 125 ISBN 0-521-30811-9                                        *
************************************************************************
      implicit double precision (a-h,o-z)
      parameter( eps=1.d-14 )
      double precision x(ndim), w(ndim)
      double precision a,b,xm,xl,z,p1,p2,p3,pp,z1,pi
      pi=4.D0*datan(1.D0)
      m=(ngauss+1)/2
      xm=0.5D0*(a+b)
      xl=0.5D0*(b-a)
      do 12 i=1,m
*         THIS IS A REALLY CLEVER ESTIMATE :
          z= dcos(pi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
    1     continue
              p1=1.D0
              p2=0.D0
              do 11 j=1,ngauss
                  p3= p2
                  p2= p1
                  p1=((dble(2*j)-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
   11         continue
              pp=ngauss*(z*p1-p2)/(z*z-1.d0)
              z1= z
              z= z1-p1/pp
          if (dabs(z-z1) .gt. eps) goto 1
          x(i)= xm-xl*z
          x(ngauss+1-i)= xm+xl*z
          w(i)=2.D0*xl/((1.D0-z*z)*pp*pp)
******          write(*,*) ' mie3_gauleg: ',i,'  x=',x(i),' w=',w(i)
          w(ngauss+1-i)= w(i)
   12 continue
      return
      end
