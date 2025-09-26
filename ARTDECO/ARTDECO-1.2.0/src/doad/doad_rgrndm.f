c-----------------------------------------------------------------------
c      
c  Author:      Michele Vesperini
c  Version:     96 04 12
c
c Adding and doubling RT code including internal source and polarization
c-----------------------------------------------------------------------

      subroutine doad_rgrndm(igrnd,m,refl,xmu,nmutot,nmat,wmu,r,
     &                                              NDcoef,NDmu,NDsup)
************************************************************************
*** fill reflection supermatrix, r, of groundsurface, m-th fourier   ***
*** coefficient                                                      ***
*** igrnd=1 : Lambert surfaces          )    0-th term only          ***
*** igrnd=2 : Lommel-Seeliger surfaces  )                            ***
*** see Van de Hulst, 1980 section 18.1.4 page 605                   ***
c -- modified to account for a directional reflectance M. Vesperini 960621
c - igrnd=3 : bidirectional reflectance )                            ***
************************************************************************
      implicit double precision (a-h,o-z)
      dimension xmu(NDmu), wmu(NDmu), r(NDsup,NDsup)
      dimension refl(4,4,NDmu,NDmu,0:NDcoef)

      nsup = nmutot*nmat

      do 10 j=1,nsup
      do 10 i=1,nsup
          r(i,j) = 0.d0
   10 continue

      if (igrnd.lt.3) then
      if (m .ne. 0) return

      albedo=refl(1,1,1,1,0)
      do 20 i=1,nmutot
          i1 = (i-1)*nmat+1
          do 20 j=1,nmutot
              j1 = (j-1)*nmat+1
              y = 1.d0
              if (igrnd .eq. 2) then
                if ((xmu(i).lt.1.d-10).and.(xmu(j).lt.1.d-10)) then
                  y = 0.d0
                else
                  y = 1.d0/(xmu(i)+xmu(j))
                endif
              endif
              r(i1,j1) = albedo*y*wmu(i)*wmu(j)
   20 continue
      else
c -- bidirectional reflectance
c -- attention
c -- j 	direction incidente
c -- i  direction reflechie
c -- il faut verifier que imat et jmat sont coherents avec ca !
      do imat=1,nmat
      do jmat=1,nmat
         do i=1,nmutot
            i1 = (i-1)*nmat+imat
            do j=1,nmutot
              j1 = (j-1)*nmat+jmat
              r(i1,j1)=refl(imat,jmat,i,j,m)*wmu(i)*wmu(j)
            enddo
         enddo
      enddo
      enddo
      endif

      return
      end
