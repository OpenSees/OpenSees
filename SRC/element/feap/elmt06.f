c.... *****************************************************************<-70
      subroutine ELMT05(d,ul,xl,ix,tl,s,r,ndf,ndm,nst,isw,dm,
     1                  nen,n,nh1,nh2,nh3,h,ctan,ior,iow)

c .....................................................................
c     elmt05 - Two dimensional truss element
c
c     written:  fmk
c     created:  03/99
c     revision: A
c .....................................................................
      implicit none

      integer ix(*),ndf,ndm,nst,isw,nh1,nh2,nh3,nen,n,ior,iow
      real*8  d(*), ul(ndf,nen,*), xl(ndm,*), tl(*), s(nst,*), r(*)
      
      real*8  p(*), h(*), ctan(*)
      
      integer  i,j
      real*8  cs, sn, L, dx, dy, k, m, A, E, rho, tran(4), eps, force

      if(isw.eq.0) then
c     output element type
         if(iow.lt.0) then
            write(*,*) '   Elmt05: 2d Truss Element.'
         else
            write(iow,*) '   Elmt05: 2d Truss Element.'
         endif

      elseif(isw.eq.1) then
c     input material properties & set nh1(=nh2), nh3
c     d(1) = A, d(2)=E, d(3) = rho, nh1 = 0, nh2 = 0

      elseif(isw.eq.2) then
c     check of mesh if desired (chec)


      elseif(isw.eq.3) then
c     compute element residual and tangent matrix

c ...... compute element length
         dx = xl(1,2) - xl(1,1)
         dy = xl(2,2) - xl(2,1)
         L = sqrt(dx*dx + dy*dy)
         cs = dx/L
         sn = dy/L
         tran(1) = cs
         tran(2) = sn
         tran(3) = -cs
         tran(4) = -sn

c ...... add stiffness portion to matrix
         A = d(1)
         E = d(2)
         if(ctan(1).ne.0.0d0) then
            k = ctan(1)*A*E/L
            do i = 1,4
               do j = 1,4
                  s(i,j) = s(i,j) + k*tran(i)*tran(j)
               enddo
            enddo
         endif

c ...... add mass portion to matrix
         rho = d(3)
         if (ctan(3).ne.0.0d0) then
            m = ctan(3)*rho*A*L/2
            do i = 1,4
               s(i,i) = s(i,i) + m
            enddo
         endif

c ...... compute the strain
         eps = cs*(ul(1,2,1)-ul(1,1,1))
         eps = eps + sn*ul(2,2,1)-ul(2,1,1)
         eps = eps/L

c .....  compute the residual 
         force = E*A*L*eps
         r(1) = cs*force
         r(2) = sn*force
         r(3) = -r(1)
         r(4) = -r(2)

c .....  add the inertia effects
         m = rho*A*L/2
         r(1) = r(1) - m*ul(1,1,5)
         r(2) = r(2) - m*ul(1,2,5)
         r(3) = r(3) - m*ul(2,1,5)
         r(4) = r(4) - m*ul(2,2,5)

      elseif(isw.eq.4) then
c     output element quantities
         write(iow,1001) n,d(1),d(2),d(3),ix(1),ix(2)

1001  format(5x,'Linear 2d Truss Element'//
     1  5x,'A:      ',e12.5/5x,'E: ',e12.5/
     2  5x,'rho: ',e12.5/5x,'end1: ',i5,'end2: ',i5)

      elseif(isw.eq.5) then
c     compute element mass matrix

c        compute element length
         dx = xl(1,2) - xl(1,1)
         dy = xl(2,2) - xl(2,1)
         L = sqrt(dx*dx + dy*dy)

c ...... compute mass matrix
         A = d(1)
         rho = d(3)
         m = ctan(3)*rho*A*L/2
         do i = 1,4
            s(i,i) = m
         enddo
         
         r(1) = m
         r(2) = m
         r(3) = m
         r(4) = m

      elseif(isw.eq.6) then
c     compute element residual

c ...... compute element length
         dx = xl(1,2) - xl(1,1)
         dy = xl(2,2) - xl(2,1)
         L = sqrt(dx*dx + dy*dy)
         cs = dx/L
         sn = dy/L

c ...... compute the strain
         eps = cs*(ul(1,2,1)-ul(1,1,1))
         eps = eps + sn*(ul(2,2,1)-ul(2,1,1))
         eps = eps/L

c .....  compute the residual 
         A = d(1)
         E = d(2)
         force = E*A*L*eps
         r(1) = cs*force
         r(2) = sn*force
         r(3) = -cs*force
         r(4) = -sn*force

c .....  add the inertia effects
         m = rho*A*L/2
         r(1) = r(1) - m*ul(1,1,5)
         r(2) = r(2) - m*ul(1,2,5)
         r(3) = r(3) - m*ul(2,1,5)
         r(4) = r(4) - m*ul(2,2,5)

      elseif(isw.eq.7) then
c     output surface loading

      elseif(isw.eq.8) then
c     compute stress projections at nodes

      endif


      end

