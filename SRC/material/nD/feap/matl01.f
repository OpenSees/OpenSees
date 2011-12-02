      subroutine MATL01(eps,trace,td,d,ud,hn,h1,nh,sig,dd,isw)

c      * * F E A P * * A Finite Element Analysis Program

c     Michael Scott
c     CE 233

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Linear Elastic Material Model

c     Input:
c          eps(*)  -  Current strains at point
c          trace   -  Trace of strain at point
c          td      -  Temperature change
c          d(*)    -  Program material parameters (ndd)
c          ud(*)   -  User material parameters (nud)
c          hn(nh)  -  History terms at point: t_n
c          h1(nh)  -  History terms at point: t_n+1
c          nh      -  Number of history terms
c          isw     -  Solution option from element

c     Output: (N.B. Arrays are set to zero before call)
c          sig(6)  -  Stresses at point.
c          dd(6,6) -  Current material tangent moduli
c-----[--.----+----.----+----.-----------------------------------------]

      implicit none
      
      integer  nh, isw
      real*8   trace, td
      real*8   eps(6),d(*),ud(*),hn(nh),h1(nh),sig(6),dd(6,6)

      real*8   E,nu
      real*8   mu,mu2,lam
      integer  i,j

c     // Read in material parameters
      E = ud(1)
      nu = ud(2)

      mu2 = E/(1.d0+nu)
      mu  = 0.5d0*mu2
      lam = mu2*nu/(1.d0-2.d0*nu)

      sig(1) = mu2*eps(1) + lam*trace
      sig(2) = mu2*eps(2) + lam*trace
      sig(3) = mu2*eps(3) + lam*trace
      sig(4) =  mu*eps(4)
      sig(5) =  mu*eps(5)
      sig(6) =  mu*eps(6)

      do i = 1,3
         dd(i  ,i  ) = mu2
         dd(i+3,i+3) = mu
         do j = 1,3
            dd(i,j) = dd(i,j) + lam
         end do
      end do
      
      end
