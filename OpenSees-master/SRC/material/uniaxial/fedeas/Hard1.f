      subroutine HARD_1(matpar,hstvP,hstv,epsP,sigP,deps,sig,tang,ist)
c I  matpar contains fixed properties (4)
c     E    = Elastic modulus             --> matpar(1)
c     sigY = Yield stress                --> matpar(2)
c     Hiso = Isotropic hardening modulus --> matpar(3)
c     Hkin = Kinematic hardening modulus --> matpar(4)

c I  hstvP contains committed history variables:
c     hstvP(1) = ep    --> effective plastic strain
c     hstvP(2) = alpha --> internal hardening variable
c     hstvP(3) = kappa --> back stress for kinematic hardening
c	 
c O  hstv will be set to the corresponding trial values of hstvP

c I  epsP: strain at last committed state
c I  sigP: stress at last committed state
c I  deps: current strain increment
c O  sig : updated stress
c O  tang: updated tangent
c I  ist : tangent calculation switch 
c             1 = tangent, 2 = incremental secant, 3 = total secant

      implicit none
 
c     Arguments
      integer ist
      real*8  matpar(4),hstvP(3),hstv(3)
      real*8  epsP,sigP,deps
      real*8  sig,tang

c     Local variables
      real*8  E,sigY,Hiso,Hkin
      real*8  ep,alpha,kappa   
      real*8  eps,f,xsi,dGamma
      integer sgn

c     Material parameters
      E    = matpar(1)
      sigY = matpar(2)
      Hiso = matpar(3)
      Hkin = matpar(4)

c     History variables
      ep    = hstvP(1)
      alpha = hstvP(2)
      kappa = hstvP(3)

c     Current strain
      eps = epsP + deps    

c     Elastic predictor
      sig = E * (eps - ep)

c     Stress relative to back stress
      xsi = sig - kappa

c     Yield function
      f = dabs(xsi) - (sigY + Hiso*alpha)

c     Inside yield surface
      if (f.le.0.0) then
         tang = E
c     Outside yield surface ... do return mapping
      else
c     Consistency parameter
         dGamma = f / (E+Hiso+Hkin)

c     Normal to yield surface
         if (xsi.le.0.d0) then
            sgn = -1
         else
            sgn = 1
         endif

c     Updated stress
         sig = sig - dGamma*E*sgn
	
c     Updated plastic strain
         ep = ep + dGamma*sgn

c     Updated back stress
         kappa = kappa + dGamma*Hkin*sgn
	
c     Updated internal hardening variable
         alpha = alpha + dGamma

c     Elasto-plastic tangent
         tang = E*(Hkin+Hiso) / (E+Hkin+Hiso)
      endif

c     Update history variables
      hstv(1) = ep
      hstv(2) = alpha
      hstv(3) = kappa

c     Compute requested tangent
      if (ist.eq.2.and.deps.ne.0.d0) then
      	tang = (sig-sigP)/deps
      else if (ist.eq.3.and.eps.ne.0.d0)then
      	tang = sig/eps
      else
c     add additional cases, if needed
      endif	         

      return

      end 
