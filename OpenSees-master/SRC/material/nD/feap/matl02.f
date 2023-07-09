      subroutine MATL02(eps,trace,td,d,ud,hn,h1,nh,sig,dd,isw)

c      * * F E A P * * A Finite Element Analysis Program

c     Michael Scott
c     CE 233

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Linear Viscoelastic Material Model

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

      include  'iofile.h'
      include  'tdata.h'

      integer  nh, isw
      real*8   trace, td
      real*8   K,G,muK,muG,lamK,lamG,theta
      real*8   eps(6),d(*),ud(7),hn(nh),h1(nh),sig(6),dd(6,6)

      integer  i,j
      real*8   tmp1,tmp2
      real*8   qvn,evn,qvn1,evn1
      real*8   pn1
      real*8   qn(6),en(6),qn1(6),en1(6)
      real*8   sn1(6)
      real*8   Kve,Gve
      real*8   one3,two3

c     // Compute useful constants
      one3 = 1.d0/3.d0
      two3 = 2.d0*one3
      
c     // Read in material parameters
      K = ud(1)
      G = ud(2)
      muK = ud(3)
      muG = ud(4)
      lamK = ud(5)
      lamG = ud(6)
      theta = ud(7)

c     // Get the strain trace
      evn1 = trace

c     // If there is a bulk viscoelastic effect
      if (muK.gt.0.d0) then
	
c       // Put strain trace in trial history
         h1(7) = evn1

c       // Get volumetric strain and internal variable from committed history
         evn = hn(7)
         qvn = hn(14)

c       // Compute temporary quantities for integrating rate equation
         tmp1 = 1.d0 + theta*dt/lamK
         tmp2 = 1.d0 - (1.d0-theta)*dt/lamK

c       // Compute q_{v,n+1} and put in trial history
         qvn1 = (tmp2*qvn + evn1 - evn) / tmp1
         h1(14) = qvn1

c       // Compute pressure p_{n+1}
         pn1 = K*((1.d0-muK)*evn1 + muK*qvn1)

c       // Compute bulk viscoelastic tangent term
         Kve = K*((1.d0-muK) + muK/tmp1)
	
c     // Otherwise, no bulk viscoelastic effect
      else

         Kve = K
         pn1 = Kve*evn1
	
      end if

c     // Add pressure contribution to material stress
c     // \sigma += m p
      do i = 1,3
         sig(i) = pn1
      end do

c     // Add bulk term to material tangent
c     // D_t += Kve * mm^T
      do i = 1,3
         do j = 1,3
            dd(i,j) = Kve
         end do      
      end do


c     // Compute trial deviatoric strains
c     // e = I_{dev} \epsilon = (I - 1/3 mm^T) \epsilon
      evn1 = evn1*one3
      en1(1) = eps(1) - evn1
      en1(2) = eps(2) - evn1
      en1(3) = eps(3) - evn1
      en1(4) = eps(4)
      en1(5) = eps(5)
      en1(6) = eps(6)

c     // If there is a shear viscoelastic effect
      if (muG.gt.0.d0) then

c       // Put deviatoric strains into trial history
         do i = 1,6
            h1(i) = en1(i)
         end do
         
c       // Compute temporary quantities for integrating rate equation
         tmp1 = 1.d0 + theta*dt/lamG
         tmp2 = 1.d0 - (1.d0-theta)*dt/lamG

c       // Get deviatoric strains and internal variables from committed history
         do i = 1,6
            en(i) = hn(i)
            qn(i) = hn(i+7)
         end do

c       // Compute internal variables and put in trial history
c       // q_{n+1} = \frac{tmp2*q_n + e_{n+1} - e_n}{tmp1}
         do i = 1,6
            qn1(i) = (tmp2*qn(i) + en1(i) - en(i)) / tmp1
            h1(i+7) = qn1(i)
         end do

c       // Compute deviatoric stress
c       // s_{n+1} = 2 G I_o (\mu_0 e_{n+1} + \mu_1 q_{n+1})
c       // Do the bulk terms ...
         do i = 1,3
            sn1(i) = G*(1.d0-muG)*2.d0*en1(i)
            sn1(i) = sn1(i) + G*muG*2.d0*qn1(i)
         end do

c       // ... and now do the shear terms
         do i = 4,6
            sn1(i) = G*(1.d0-muG)*en1(i)
            sn1(i) = sn1(i) + G*muG*qn1(i)
         end do

c       // Compute shear viscoelastic tangent term
         Gve = G*((1.d0-muG) + muG/tmp1)

c     // Otherwise, no shear viscoelastic effect
      else

         Gve = G

c       // Do the bulk terms ...
         do i = 1,3
	    sn1(i) = G*2.d0*en1(i)
         end do

c       // ... and now do the shear terms
         do i = 4,6
	    sn1(i) = G*en1(i)
         end do

      end if

c     // Add deviatoric stress to total stress
      do i = 1,6
         sig(i) = sig(i) + sn1(i)
      end do

c     // Add deviatoric quantities to material tangent
c     // D_t += 2 G^{ve} I_o (I - 1/3 mm^T)
c     // Do the bulk terms ...
      do i = 1,3
         dd(i,i) = dd(i,i) + 2.d0*Gve*two3
      end do

c     // ... and do the deviatoric terms ...
      do i = 4,6
         dd(i,i) = dd(i,i) + Gve
      end do

c     // ... and finally do the off-diagonal terms
      Gve = Gve*two3            ! Gve = 2.d0*Gve*one3
      dd(1,2) = dd(1,2) - Gve
      dd(1,3) = dd(1,3) - Gve
      dd(2,3) = dd(2,3) - Gve
      dd(2,1) = dd(2,1) - Gve
      dd(3,1) = dd(3,1) - Gve
      dd(3,2) = dd(3,2) - Gve


      end
