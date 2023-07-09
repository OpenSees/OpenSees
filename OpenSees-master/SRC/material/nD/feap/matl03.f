      subroutine MATL03(eps,trace,td,d,ud,hn,h1,nh,sig,dd,isw)

c      * * F E A P * * A Finite Element Analysis Program

c     Michael Scott
c     CE 233

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Linear Isotropic Hardening Material Model

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
      include  'counts.h'

      integer  nh, isw
      real*8   trace, td
      real*8   eps(6),d(*),ud(4),hn(nh),h1(nh),sig(6),dd(6,6)

      real*8   K,G,sigY,H,twoG
      real*8   epn(9),en1(9),syn,kpn,kpn1
      real*8   pn1,sn1(9)
      real*8   one3,two3,root23

      real*8   evn1,nsn1,Fn1,nn1(9),dlam

      integer  i,j
      real*8   tmp1,tmp2,a(6)

c     // Compute useful constants
      one3 = 1.d0/3.d0
      two3 = 2.d0*one3
      root23 = sqrt(two3)
      
c     // Read in material parameters
      K = ud(1)
      G = ud(2)
      sigY = ud(3)
      H = ud(4)
      twoG = 2.d0*G

c     // Add pressure contribution to material stress
c     // \sigma += m p
      pn1 = K*trace
      do i = 1,3
         sig(i) = pn1
      end do
      
c     // Add bulk term to material tangent
c     // D_t += K * mm^T
      do i = 1,3
         do j = 1,i
            dd(i,j) = K
         end do      
      end do
      
c     // Compute trial deviatoric strains in 9-component form
c     // e_{n+1} = I_{dev} \epsilon = (I - 1/3 mm^T) \epsilon
      evn1 = trace*one3
      en1(1) = eps(1) - evn1
      en1(2) = eps(2) - evn1
      en1(3) = eps(3) - evn1
      en1(4) = 0.5d0*eps(4)
      en1(5) = en1(4)
      en1(6) = 0.5d0*eps(5)
      en1(7) = en1(6)
      en1(8) = 0.5d0*eps(6)
      en1(9) = en1(8)

c     // Get committed deviatoric plastic strains from history
c     // and map from 6- to 9-component form
      do i = 1,3
         epn(i) = hn(i)
      end do
      epn(4) = 0.5d0*hn(4)
      epn(5) = epn(4)
      epn(6) = 0.5d0*hn(5)
      epn(7) = epn(6)
      epn(8) = 0.5d0*hn(6)
      epn(9) = epn(8)

c     // Get committed effective plastic strain and compute yield stress
      kpn = hn(7)
      syn = sigY + H*kpn

c     // Compute trial deviatoric stress in 9-component form
c     // s_{n+1} = 2 G (e_{n+1} - e^p_n)
c     // Do the bulk terms ...
      do i = 1,9
         sn1(i) = twoG*(en1(i)-epn(i))
      end do

c     // Compute norm of trial deviatoric stress
      nsn1 = 0.d0
      do i = 1,9
         nsn1 = nsn1 + sn1(i)*sn1(i)
      end do
      nsn1 = sqrt(nsn1)

c     // Check yield condition
c     // F = \|s\| - sqrt(2/3)*\sigma_y
      Fn1 = nsn1 - root23*syn

c     // First iteration or an elastic step
      if (niter.eq.0 .or. Fn1.lt.0.d0) then

c       // No change in history variables n+1 <-- n
         do i = 1,7
	    h1(i) = hn(i)
         end do

c       // Add deviatoric part to material stress after mapping
c       // to 6-component form
         do i = 1,3
	    sig(i) = sig(i) + sn1(i)
         end do
         sig(4) = sig(4) + 0.5d0*(sn1(4)+sn1(5))
         sig(5) = sig(5) + 0.5d0*(sn1(6)+sn1(7))
         sig(6) = sig(6) + 0.5d0*(sn1(8)+sn1(9))
         
c       // Add deviatoric part to material tangent after mapping to
c       // 6-component form
c       // D_t += 2 G I_o I_{dev}
c       // Do the bulk terms ...
         tmp1 = twoG*two3
         do i = 1,3
            dd(i,i) = dd(i,i) + tmp1
         end do
c       // ... and do the deviatoric terms ...
         do i = 4,6
            dd(i,i) = dd(i,i) + G
         end do
c       // ... and finally do the off-diagonal terms
         tmp1 = G*two3          ! G*two3 == 2.d0*G*one3
         dd(2,1) = dd(2,1) - tmp1
         dd(3,1) = dd(3,1) - tmp1
         dd(3,2) = dd(3,2) - tmp1

c     // Plastic step
      else

c       // Compute normal to yield surface
c       // n_{n+1} = \frac{s_{n+1}}{\|s_{n+1}\|}
         tmp1 = 1.d0/nsn1       ! Save some flops
         do i = 1,9
	    nn1(i) = sn1(i)*tmp1
         end do
        
c       // Compute consistency parameter
         dlam = Fn1 / (twoG+two3*H)

c       // Update plastic deviatoric strains in 9-component form
         do i = 1,9
	    epn(i) = epn(i) + dlam*nn1(i)
         end do

c       // Update effective plastic strain
c       // \kappa{n+1} = \kappa_n + sqrt(2/3) (\Delta \lambda)
         kpn1 = kpn + root23*dlam

c       // Put trial history variables into vector after mapping back
c       // to 6-component form
         do i = 1,3
	    h1(i) = epn(i)
         end do
         h1(4) = 2.d0*epn(4)
         h1(5) = 2.d0*epn(6)
         h1(6) = 2.d0*epn(8)
         h1(7) = kpn1

c       // Return deviatoric stress, in 9-component form, to
c       // yield surface
         tmp1 = twoG*dlam
         do i = 1,9
            sn1(i) = sn1(i) - tmp1*nn1(i)
         end do

c       // Add deviatoric part to material stress after mapping to
c       // 6-component form
         do i = 1,3
	    sig(i) = sig(i) + sn1(i)
         end do
         sig(4) = sig(4) + 0.5d0*(sn1(4)+sn1(5))
         sig(5) = sig(5) + 0.5d0*(sn1(6)+sn1(7))
         sig(6) = sig(6) + 0.5d0*(sn1(8)+sn1(9))

c       // Add deviatoric part to material tangent
c       // D_t += tmp_1 I_o I_{dev}
         tmp1 = twoG*(1.d0-twoG*dlam/nsn1)
c       // Do the bulk terms ...
         tmp2 = two3*tmp1
         do i = 1,3
            dd(i,i) = dd(i,i) + tmp2
         end do
c       // ... and do the deviatoric terms ...
         tmp2 = 0.5d0*tmp1
         do i = 4,6
            dd(i,i) = dd(i,i) + tmp2
         end do
c       // ... and finally do the off-diagonal terms
         tmp2 = one3*tmp1
         dd(2,1) = dd(2,1) - tmp2
         dd(3,1) = dd(3,1) - tmp2
         dd(3,2) = dd(3,2) - tmp2

c       // Map normal vector to 6-component form
c       // n_6 = P^T n_9
         nn1(4) = 0.5d0*(nn1(4)+nn1(5))
         nn1(5) = 0.5d0*(nn1(6)+nn1(7))
         nn1(6) = 0.5d0*(nn1(8)+nn1(9))

c       // Compute temporary 6-component vector
c       // a = n^T I_{dev} = n - 1/3 (n^T m) m
         tmp1 = one3*(nn1(1)+nn1(2)+nn1(3))
         do i = 1,3
	    a(i) = nn1(i)-tmp1
         end do
         do i = 4,6
	    a(i) = nn1(i)
         end do

c       // Add deviatoric part to material tangent
c       // D_t += tmp_2 n n^T I_{dev} = tmp2 n a^T
         tmp2 = twoG*G*(1.d0/(G+one3*H)-2.d0*dlam/nsn1)
         do i = 1,6
	    do j = 1,i
               dd(i,j) = dd(i,j) - tmp2*nn1(i)*a(j)
	    end do
         end do
         
      end if

c     // Complete symmetric tangent matrix
      do i = 1,6
         do j = 1,i-1
	    dd(j,i) = dd(i,j)
         end do
      end do
      
      end
