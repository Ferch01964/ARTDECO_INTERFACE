c-----------------------------------------------------------------------
c      
c  Author:      Michele Vesperini
c  Version:     96 04 12
c
c Adding and doubling RT code including internal source and polarization
c-----------------------------------------------------------------------


      subroutine doad_gau_xmu(NDmu,NDmuext,nmug,nmuext,xmuext,xmu,wmu)
c-----------------------------------------------------------------------
c -- compute xmu and corresponding gaussian weights
c -- for mu integration
c -- xmu array is filled with gauss xmu values and external xmu
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

      dimension xmuext(NDmuext)
      dimension xmu(NDmu), wmu(NDmu)

      nmutot = nmug+nmuext
************************************************************************
* test the Gauss-Legendre points and weights for zenith integration    *
************************************************************************
      call doad_gauleg(NDmu,nmug,0.d0,1.d0,xmu,wmu,1.d-11)
      do 10 i=1,nmuext
        xmu(nmug+i) = xmuext(i)
   10 continue
************************************************************************
* test the Gauss-Legendre points and weights for zenith integration    *
************************************************************************
      sum = 0.d0
c$$$      write(*,5)
c$$$      write(*,5)
c$$$    5 format(1h1,'Gaussian division points for zenith integration',/)
      do 20 i=1,nmug
        sum = sum+wmu(i)*xmu(i)
c$$$        write(*,15) i,xmu(i),wmu(i)
c$$$        write(*,15) i,xmu(i),wmu(i)
c$$$   15   format(1h ,'i =',i5,2x,'xmu =',f17.15,4x,'wmu =',f17.15)
   20 continue
c$$$      write(*,25) sum
c$$$      write(*,25) sum
c$$$   25 format(1h0,' int_0^1 x dx = ',f17.15)
      if (dabs(sum-.5d0) .gt. 1.d-11) write(*,35)
      if (dabs(sum-.5d0) .gt. 1.d-11) write(*,35)
   35 format(1h0,'Gauss integration over zenith is inaccurate !')
************************************************************************
* introduce wmu(nmutot), the weight factor of a supermatrix            *
************************************************************************
      do 30 i=1,nmutot
        if (i .le. nmug) then
          wmu(i) = dsqrt(2.d0*wmu(i)*xmu(i))
        else
          wmu(i) = 1.d0
        endif
   30 continue
      return
      end
