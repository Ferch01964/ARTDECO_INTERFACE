c-----------------------------------------------------------------------
c      
c  Author:      Michele Vesperini
c  Version:     96 04 12
c
c Adding and doubling RT code including internal source and polarization
c-----------------------------------------------------------------------


      subroutine doad_ssurf(nmat,nmutot,emis,wmu,ssor,NDcoef,
     &                                         NDsup,NDmu,NDlay,sr,st)
************************************************************************
*** calculate the zero order radiation field of a radiating surface  ***
*** the strength of the sources is contained in ssor(0)              ***
*** see Wauben (1990)                                                ***
*** results are stored in supervector sr(i)                          ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension wmu(NDmu)
      dimension ssor(0:NDlay)
      dimension sr(NDsup), st(NDsup)
cmv   added
      dimension emis(4,NDmu)

************************************************************************
* initializing the supervectors sr and st                              *
************************************************************************
      nsup = nmat*nmutot
      do 10 i=1,nsup
        sr(i) = 0.d0
        st(i) = 0.d0
   10 continue
************************************************************************
* loop 20 over the directions mu_i                                     *
************************************************************************

      do 20 i=1,nmutot
        im = nmat*(i-1)
cmv     sr(im) = ssor(0)*wmu(i)
cmv-added
        do imat=1,nmat
           sr(im+imat) = ssor(0)*wmu(i)*emis(imat,i)
        enddo 
cmv-added
   20 continue
      return
      end
