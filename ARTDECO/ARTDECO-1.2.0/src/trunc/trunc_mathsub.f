

cFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
C  (C) Copr. 1986-92 Numerical Recipes Software 71a.
      FUNCTION pythag(a,b)
      REAL(8) a,b,pythag
      REAL(8) absa,absb

      absa=dabs(a)
      absb=dabs(b)
      if(absa.gt.absb)then
         pythag=absa*dsqrt(1.D0+(absb/absa)**2)
      else
         if(absb.eq.0.D0)then
           pythag=0.D0
         else
           pythag=absb*dsqrt(1.D0+(absa/absb)**2)
         endif
      endif
      return
      END FUNCTION pythag

cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C  (C) Copr. 1986-92 Numerical Recipes Software 71a.
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)

      INTEGER m,mp,n,np,NMAX
      REAL(8) a(mp,np),v(np,np),w(np)
      PARAMETER(NMAX=1001)
C     U    USES pythag
      INTEGER i,its,j,jj,k,l,nm
      REAL(8) anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag

      if (np .gt. NMAX) then
         write(*,*) '' 
         write(*,*) '================================================' 
         write(*,*) ' (in trunc_mathsub.f -svdcmp-) : ERROR'
         write(*,*) ' nbetal+1 (=np) MUST BE smaller or equal to NMAX'
         write(*,*) '================================================' 
         write(*,*) '' 
         stop
      endif
      
      g=0.0D0
      scale=0.0D0
      anorm=0.0D0

      do 25 i=1,n
         l=i+1
         rv1(i)=scale*g
         g=0.0D0
         s=0.0D0
         scale=0.0D0
         if(i.le.m)then
            do 11 k=i,m
               scale=scale+dabs(a(k,i))
 11         continue
            if(scale.ne.0.0D0)then
               do 12 k=i,m
                  a(k,i)=a(k,i)/scale
                  s=s+a(k,i)*a(k,i)
 12            continue
               f=a(i,i)
               g=-sign(dsqrt(s),f)
               h=f*g-s
               a(i,i)=f-g
               do 15 j=l,n
                  s=0.0D0
                  do 13 k=i,m
                     s=s+a(k,i)*a(k,j)
 13               continue
                  f=s/h
                  do 14 k=i,m
                     a(k,j)=a(k,j)+f*a(k,i)
 14               continue
 15            continue
               do 16 k=i,m
                  a(k,i)=scale*a(k,i)
 16            continue
            endif
         endif
         w(i)=scale *g
         g=0.0D0
         s=0.0D0
         scale=0.0D0
         if((i.le.m).and.(i.ne.n))then
            do 17 k=l,n
               scale=scale+dabs(a(i,k))
 17         continue
            if(scale.ne.0.0D0)then
               do 18 k=l,n
                  a(i,k)=a(i,k)/scale
                  s=s+a(i,k)*a(i,k)
 18            continue
               f=a(i,l)
               g=-sign(dsqrt(s),f)
               h=f*g-s
               a(i,l)=f-g
               do 19 k=l,n
                  rv1(k)=a(i,k)/h
 19            continue
               do 23 j=l,m
                  s=0.0D0
                  do 21 k=l,n
                     s=s+a(j,k)*a(i,k)
 21               continue
                  do 22 k=l,n
                     a(j,k)=a(j,k)+s*rv1(k)
 22               continue
 23            continue
               do 24 k=l,n
                  a(i,k)=scale*a(i,k)
 24            continue
            endif
         endif
         anorm=max(anorm,(dabs(w(i))+dabs(rv1(i))))
 25   continue
      do 32 i=n,1,-1
         if(i.lt.n)then
            if(g.ne.0.0D0)then
               do 26 j=l,n
                  v(j,i)=(a(i,j)/a(i,l))/g
 26            continue
               do 29 j=l,n
                  s=0.0D0
                  do 27 k=l,n
                     s=s+a(i,k)*v(k,j)
 27               continue
                  do 28 k=l,n
                     v(k,j)=v(k,j)+s*v(k,i)
 28               continue
 29            continue
            endif
            do 31 j=l,n
               v(i,j)=0.0D0
               v(j,i)=0.0D0
 31         continue
         endif
         v(i,i)=1.0D0
         g=rv1(i)
         l=i
 32   continue
      do 39 i=min(m,n),1,-1
         l=i+1
         g=w(i)
         do 33 j=l,n
            a(i,j)=0.0
 33      continue
         if(g.ne.0.0D0)then
            g=1.0D0/g
            do 36 j=l,n
               s=0.0D0
               do 34 k=l,m
                  s=s+a(k,i)*a(k,j)
 34            continue
               f=(s/a(i,i))*g
               do 35 k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
 35            continue
 36         continue
            do 37 j=i,m
               a(j,i)=a(j,i)*g
 37         continue
         else
            do 38 j= i,m
               a(j,i)=0.0D0
 38         continue
         endif
         a(i,i)=a(i,i)+1.0D0
 39   continue
      do 49 k=n,1,-1
         do 48 its=1,30
            do 41 l=k,1,-1
               nm=l-1
               if((dabs(rv1(l))+anorm).eq.anorm)  goto 2
               if((dabs(w(nm))+anorm).eq.anorm)  goto 1
 41         continue 
 1          c=0.0D0
            s=1.0D0
            do 43 i=l,k
               f=s*rv1(i)
               rv1(i)=c*rv1(i)
               if((dabs(f)+anorm).eq.anorm) goto 2
               g=w(i)
               h=pythag(f,g)
               w(i)=h
               h=1.0D0/h
               c= (g*h)
               s=-(f*h)
               do 42 j=1,m
                  y=a(j,nm)
                  z=a(j,i)
                  a(j,nm)=(y*c)+(z*s)
                  a(j,i)=-(y*s)+(z*c)
 42            continue
 43         continue 
 2          z=w(k)
            if(l.eq.k)then
               if(z.lt.0.0D0)then
                  w(k)=-z
                  do 44 j=1,n
                     v(j,k)=-v(j,k)
 44               continue
               endif
               goto 3
            endif
            if(its.eq.100) pause 'no convergence in svdcmp'
            x=w(l)
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0D0*h*y)
            g=pythag(f,1.0D0)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c=1.0D0
            s=1.0D0
            do 47 j=l,nm
               i=j+1
               g=rv1(i)
               y=w(i)
               h=s*g
               g=c*g
               z=pythag(f,h)
               rv1(j)=z
               c=f/z
               s=h/z
               f= (x*c)+(g*s)
               g=-(x*s)+(g*c)
               h=y*s
               y=y*c
               do 45 jj=1,n
                  x=v(jj,j)
                  z=v(jj,i)
                  v(jj,j)= (x*c)+(z*s)
                  v(jj,i)=-(x*s)+(z*c)
 45            continue
               z=pythag(f,h)
               w(j)=z
               if(z.ne.0.0D0)then
                  z=1.0D0/z
                  c=f*z
                  s=h*z
               endif
               f= (c*g)+(s*y)
               x=-(s*g)+(c*y)
               do 46 jj=1,m
                  y=a(jj,j)
                  z=a(jj,i)
                  a(jj,j)= (y*c)+(z*s)
                  a(jj,i)=-(y*s)+(z*c)
 46            continue
 47         continue
            rv1(l)=0.0D0
            rv1(k)=f
            w(k)=x
 48      continue 
 3       continue
 49   continue

      return

      END SUBROUTINE svdcmp


cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C  (C) Copr. 1986-92 Numerical Recipes Software 71a.
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      INTEGER m,mp,n,np,NMAX
      REAL(8) b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=10000)
      INTEGER i,j,jj
      REAL(8) s,tmp(NMAX)

      do 12 j=1,n
         s=0.D0
         if(w(j).ne.0.D0)then
            do 11 i=1,m
               s=s+u(i,j)*b(i)
 11         continue
            s=s/w(j)
         endif
         tmp(j)=s
 12   continue
      do 14 j=1,n
         s=0.D0
         do 13 jj=1,n
            s=s+v(j,jj)*tmp(jj)
 13      continue
         x(j)=s
 14   continue
      return

      END SUBROUTINE svbksb

cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C  (C) Copr. 1986-92 Numerical Recipes Software 71a.
C
      SUBROUTINE fleg(x,pl,nl)

      INTEGER nl
      REAL(8) x,pl(nl)
      INTEGER j
      REAL(8) d,f1,f2,twox

      pl(1) = 1.D0
      pl(2) = x
      if (nl .gt. 2) then
         twox = 2.D0*x
         f2 = x
         d = 1.D0
         do 11 j = 3,nl
            f1 = d
            f2 = f2+twox
            d = d+1.D0
            pl(j) = (f2*pl(j-1)-f1*pl(j-2))/d
 11      continue
      endif

      return

      END SUBROUTINE fleg


cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      SUBROUTINE fleg2(x,pl,nl)

      INTEGER nl
      REAL(8) x,pl(0:nl)
      INTEGER j
      REAL(8) d,f1,f2,qroot6

      qroot6 = -dsqrt(6.0D0) / 4.0D0
      pl(0) = 0.0D0
      pl(1) = 0.0D0
      pl(2) = qroot6 * (1.D0-x**2)
      if (nl .gt. 2) then
         do 11 j = 2,nl-1
c     d = dsqrt(dble((j-2)*(j+2)))
c     f1 = dble(2*j-1) / d
c     f2 = dsqrt(dble((j-1)**2-4)) / d
            d = dsqrt(dble((j-1)*(j+3)))
            f1 = dble(2*j+1) / d
            f2 = dsqrt(dble(j**2-4)) / d
            pl(j+1) = f1*x*pl(j) - f2*pl(j-1)
 11      continue
      endif

      return

      END SUBROUTINE fleg2


cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      SUBROUTINE fleg3(x,pl,nl)

      INTEGER nl
      REAL(8) x,pl(0:nl)
      INTEGER j
      REAL(8) d,f1,f2,f3

      pl(0) = 0.D0
      pl(1) = 0.D0
      pl(2) = (1.D0 + x**2) / 4.D0
      if (nl .gt. 2) then
         do 11 j = 2,nl-1
            d = dble(j*(j-1)*(j+3))
            f1 = dble(j*(j+1)*(2*j+1)) / d
            f2 = -4.D0*dble(2*j+1) / d
            f3 = dble((j+1)*(j-2)*(j+2)) / d
c     d = dble((j-1)*(j-2)*(j+2))
c     f1 = dble((j-1)*(j)*(2*j-1)) / d
c     f2 = -4.D0*dble(2*j-1) / d
c     f3 = dble((j)*(j-3)*(j+1)) / d
            pl(j+1) = (f1*x+f2)*pl(j) - f3*pl(j-1)
 11      continue
      endif

      return

      END SUBROUTINE fleg3


cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      SUBROUTINE fleg4(x,pl,nl)

      INTEGER nl
      REAL(8) x,pl(0:nl)
      INTEGER j
      REAL(8) d,f1,f2,f3

      pl(0) = 0.D0
      pl(1) = 0.D0
      pl(2) = (1.D0 - x**2) / 4.D0
      if (nl .gt. 2) then
         do 11 j = 2,nl-1
            d = dble(j*(j-1)*(j+3))
            f1 = dble(j*(j+1)*(2*j+1)) / d
            f2 = 4.D0*dble(2*j+1) / d
            f3 = dble((j+1)*(j-2)*(j+2)) / d
c     d = dble((j-1)*(j-2)*(j+2))
c     f1 = dble((j-1)*(j)*(2*j-1)) / d
c     f2 = 4.D0*dble(2*j-1) / d
c     f3 = dble((j)*(j-3)*(j+1)) / d
            pl(j+1) = (f1*x+f2)*pl(j) - f3*pl(j-1)
 11      continue
      endif

      return

      END SUBROUTINE fleg4



cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C  (C) Copr. 1986-92 Numerical Recipes Software 71a.
c
      SUBROUTINE sort2(n,arr,brr)

      INTEGER n,M,NSTACK
      REAL(8) arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL(8) a,b,temp

      jstack = 0
      l = 1
      ir = n

 1    if (ir-l .lt. M) then
         do 12 j=l+1,ir
            a=arr(j)
            b=brr(j)
            do 11 i=j-1,1,-1
               if(arr(i).le.a)goto 2
               arr(i+1)=arr(i)
               brr(i+1)=brr(i)
 11         continue
            i=0 
 2          arr(i+1)=a
            brr(i+1)=b
 12      continue
         if (jstack .eq. 0) RETURN
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
         temp=brr(k)
         brr(k)=brr(l+1)
         brr(l+1)=temp
         if (arr(l+1) .gt. arr(ir)) then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
            temp=brr(l+1)
            brr(l+1)=brr(ir)
            brr(ir)=temp
         endif
         if (arr(l) .gt. arr(ir)) then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
            temp=brr(l)
            brr(l)=brr(ir)
            brr(ir)=temp
         endif
         if (arr(l+1) .gt. arr(l))then
            temp=arr(l+1)
            arr(l+1)=arr(l)
            arr(l)=temp
            temp=brr(l+1)
            brr(l+1)=brr(l)
            brr(l)=temp
         endif
         i=l+1
         j=ir
         a=arr(l)
         b=brr(l) 
 3       continue
         i=i+1
         if (arr(i) .lt. a) goto 3
 4       continue
         j=j-1
         if (arr(j) .gt. a) goto 4
         if (j .lt. i) goto 5
         temp=arr(i)
         arr(i)=arr(j)
         arr(j)=temp
         temp=brr(i)
         brr(i)=brr(j)
         brr(j)=temp
         goto 3
 5       arr(l)=arr(j)
         arr(j)=a
         brr(l)=brr(j)
         brr(j)=b
         jstack=jstack+2
         if (jstack .gt. NSTACK) pause 'NSTACK too small in sort2'
         if (ir-i+1 .ge. j-l) then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif

      goto 1

      END SUBROUTINE sort2


cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE myspline(n,x_in,y_in,nn,xn,yn)
c---------spline fit to derive the yn value at point xn
c     Inputs:
c     n:    the length of x and y
c     x(n): the x values which x(1) < x(2) ... < x(n)
c     y(n): the y value which correspondent to x(n)
c     nn:  the length of vector xx and yy
c     xn:  the x value at which y value is wanted
c     
c     Outputs:
c     yn: the wanted y value from the fitting
c     
c     Internal variables:
c     yp1: the derivative of y over x at x(1), for natural bc, yp1=1.e31
c     ypn: the derivative of y over x at x(n), for natural bc, ypn=1.e31
c     y2(n): the second derivatives
c     
      integer n,nn,ny2,i
      parameter (ny2=5000)
      real(8) xn(nn),yn(nn),y2(ny2),xx,yy,yp1,ypn
      real(8) x(n),y(n),x_in(n),y_in(n)

c--------the sorting which makes sure x(1)<x(2)<...<x(n)-------
      x = x_in
      y = y_in
      CALL sort2(n,x,y)

c--------start spline------------
      yp1 = 1.D31
      ypn = 1.D31
      CALL spline(x,y,n,yp1,ypn,y2)

      DO i = 1,nn
         xx = xn(i)
         CALL splint(x,y,y2,n,xx,yy)
         yn(i) = yy
      ENDDO                     !i

      RETURN

      END SUBROUTINE myspline
	
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

      INTEGER n,NMAX, ny2
      PARAMETER (NMAX=5000)
      parameter (ny2=5000)
      REAL(8) yp1,ypn,x(n),y(n),y2(ny2)
      INTEGER i,k
      REAL(8) p,qn,sig,un,u(NMAX)

      IF (yp1 .gt. .99D30) THEN
         y2(1) = 0.D0
         u(1) = 0.D0
      ELSE
         y2(1) = -0.5D0
         u(1) = (3.D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      ENDIF

      DO 11 i = 2,n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.D0
         y2(i) = (sig-1.D0)/p
         u(i) = (6.D0*((y(i+1)-y(i))/(x(i+1)
     +        -x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))
     +        -sig*u(i-1))/p
 11   CONTINUE

      IF (ypn .gt. .99D30) THEN
         qn = 0.D0
         un = 0.D0
      ELSE
         qn = 0.5D0
         un = (3.D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      ENDIF

      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.D0)
      DO 12 k = n-1,1,-1
         y2(k) = y2(k)*y2(k+1)+u(k)
 12   CONTINUE
      
      RETURN

      END SUBROUTINE spline


cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE splint(xa,ya,y2a,n,x,y)

      INTEGER n, ny2
      parameter (ny2=5000)
      REAL(8) x,y,xa(n),y2a(ny2),ya(n)
      INTEGER k,khi,klo
      REAL(8) a,b,h

      klo = 1
      khi = n
 1    IF (khi-klo .gt. 1) THEN
         k=(khi+klo)/2
         IF (xa(k) .gt. x) THEN
            khi=k
         ELSE
            klo=k
         ENDIF
         GOTO 1
      ENDIF

      h = xa(khi)-xa(klo)

      IF (h .eq. 0.) pause 'bad xa input in splint'
      
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+
     +     (b**3-b)*y2a(khi))*(h*h)/6.D0

      RETURN

      END SUBROUTINE splint


 
