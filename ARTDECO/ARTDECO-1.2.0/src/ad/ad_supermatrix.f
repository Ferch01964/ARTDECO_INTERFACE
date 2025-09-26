*----------------------------------------------------------------------*
*  Here are the subroutine from the supermatrix package to be called
*        by the addoub code  
*----------------------------------------------------------------------*
      subroutine addSM(NDsup, C, A, B, nmat, nmutot )
*----------------------------------------------------------------------*
*  Calculate the supermatrix sum C = A+B                               *
*  On entry:                                                           *
*     A, B    : supermatrices to be added                              *
*     nmat    : number of elements of the Stokes vector taken into     *
*               account (4 = full polarization, 3 = 3x3 approximation, *
*               2 = illegal, 1 = scalar)                               *
*     nmutot  : total number of distinct mu points                     *
*  On exit:                                                            *
*     C       : supermatrix sum of A and B                             *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup,NDsup), C(NDsup,NDsup)
      nsup = nmutot*nmat
      do 200 j=1, nsup
          do 100 i=1, nsup
              C(i,j) = A(i,j)+B(i,j)
  100     continue
  200 continue
      return
      end
      subroutine assign(NDsup, A, B, nmat, nmutot )
*----------------------------------------------------------------------*
*  Calculate the supermatrix assignment A = B                          *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup,NDsup)
      nsup = nmutot*nmat
      do 200 j=1, nsup
          do 100 i=1, nsup
              A(i,j) = B(i,j)
  100     continue
  200 continue
      return
      end
      subroutine prod(NDsup, A, B, C, nmat, nmutot, nmug )
*----------------------------------------------------------------------*
*  Calculate the supermatrix product A = B * C                         *
*  Usually a large fraction of the execution time is spent in this     *
*  subroutine, especially with polarization.                           *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup,NDsup), C(NDsup,NDsup)
      nsup = nmutot*nmat
      ng   = nmug*nmat
      do 200 j=1, nsup
          do 100 i=1, nsup
              A(i,j) = 0.D0
              do 300 k=1, ng
                  A(i,j) = A(i,j) + B(i,k)*C(k,j)
  300         continue
  100     continue
  200 continue
      return
      end
      subroutine star( NDsup, A, B, nmat, nmutot )
*----------------------------------------------------------------------*
*  Calculate A = B* where B stands for a reflection or transmission    *
*  supermatrix of a HOMOGENEOUS layer, and the star denotes            *
*  illumination from below :  A = q4q3 B q3q4.                         *
*  Eqs. (98)-(99) of de Haan et al. (1987)                             *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup,NDsup)
      call assign(NDsup, A, B, nmat, nmutot )
      if (nmat .ge. 3) then
          nsup = nmutot*nmat
          do 300 mu=1, nmutot
              ibase = (mu-1)*nmat
              do 200 k=3, nmat
                  i = ibase+k
                  do 100 j=1, nsup
                      A(i,j) = -A(i,j)
                      A(j,i) = -A(j,i)
  100             continue
  200         continue
  300     continue
      endif
      return
      end
      subroutine trace(NDsup, A, trA, nmat, nmug )
*----------------------------------------------------------------------*
*  Calculate the truncated supermatrix trace(A).                       *
*  The sum runs over the integration points only, see remark below     *
*  Eq. (124) of de Haan et al. (1987)                                  *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup)
      ng = nmug*nmat
      trA = 0.D0
      do 100 i=1, ng
          trA = trA + A(i,i)
  100 continue
      return
      end
      subroutine rdiapr(NDsup, A, B, E, nmat, nmutot )
*----------------------------------------------------------------------*
*  Calculate the full EXTENDED supermatrix product A = B E where       *
*  E is diagonal. (rdiapr = 'right diagonal product')                  *
*  The reason why the product is not limited to the integration        *
*  points is explained below Eq. (95) of de Haan et al. (1987).        *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup,NDsup), E(NDsup)
      nsup = nmutot*nmat
      do 200 j=1, nsup
          do 100 i=1, nsup
              A(i,j) = B(i,j)*E(j)
  100     continue
  200 continue
      return
      end
      subroutine ldiapr( NDsup, A, E, B, nmat, nmutot )
*----------------------------------------------------------------------*
*  Calculate the full EXTENDED supermatrix product A = E B where       *
*  E is diagonal. (ldiapr = 'left diagonal product')                   *
*  The reason why the product is not limited to the integration        *
*  points is explained below Eq. (95) of de Haan et al. (1987).        *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup), B(NDsup,NDsup), E(NDsup)
      nsup = nmutot*nmat
      do 200 j=1, nsup
          do 100 i=1, nsup
              A(i,j) = E(i)*B(i,j)
  100     continue
  200 continue
      return
      end
      subroutine Tstar(NDsup, Ts, T, nmat, nmutot )
*----------------------------------------------------------------------*
*  Calculate the transmission supermatrix Ts for illumination from     *
*  below from the normal transmission supermatrix T by symmetry :      *
*                         Ts = q3 T~ q3                                *
*  where T~ is the transpose of T, and q3 is defined above Eq. (96)    *
*  This symmetry is also valid for vertically inhomogeneous            *
*  atmospheres. It is described in Hovenier (1970) page I-6, however,  *
*  one should be aware of the difference between Hovenier's operator   *
*  T(mu,mu0,phi-phi0) and our supermatrix Tm(mu,mu0) !!                *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension T(NDsup,NDsup), Ts(NDsup,NDsup)
      nsup = nmutot*nmat
*----------------------------------------------------------------------*
*  Transpose T and put it in Ts                                        *
*----------------------------------------------------------------------*
      do 200 j=1, nsup
          do 100 i=1, nsup
              Ts(i,j) = T(j,i)
  100     continue
  200 continue
*----------------------------------------------------------------------*
*  Put a minus sign in every third row and in every third column       *
*----------------------------------------------------------------------*
      if (nmat .ge. 3) then
          do 400 i=3, nsup, nmat
              do 300 j=1, nsup
                  Ts(i,j) = -Ts(i,j)
                  Ts(j,i) = -Ts(j,i)
  300         continue
  400     continue
      endif
      return
      end
      subroutine zero(NDsup, A, nsup )
*----------------------------------------------------------------------*
*  Set supermatrix A to zero.                                          *
*----------------------------------------------------------------------*
      implicit double precision (a-h,o-z)
      dimension A(NDsup,NDsup)
      do 200 j=1, nsup
          do 100 i=1, nsup
              A(i,j) = 0.D0
  100     continue
  200 continue
      return
      end
