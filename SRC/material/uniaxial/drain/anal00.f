c **********************************************************************
      SUBROUTINE RESP00(kresis,ksave,kgem,kstep,ndof,kst,kenr,ener,
     1     ened,enso,beto,relas,rdamp,rinit,ddise,
     2     dise,vele)
c *********************************************************************
c  DRAIN-2DX ROTATIONAL SPRING ELEMENT - Bilinear return map
c ---------------------------------------------------------------------
c  PURPOSE
c     The purpose of the subroutine resp00() is to update the element 
c     state, form the static and damping resisting forces, perform
c     energy calculations, update response quantities, and put 
c     element results in /thelm/ for saving or printing.
c ---------------------------------------------------------------------
c  CALLED FROM:    respxx 
c ---------------------------------------------------------------------
c  ARGUMENTS
c   Input:
c     kresis      - indicator for calculating resisting forces
c                   (1=static, 2=static and damping)
c     ksave       - indicator for saving element results
c                   (1=yes, 0=no)
c     kgem        - second-order analysis code 
c                   (not used for this element)
c     ndof        - number of element DOFs
c     kenr        - energy calculation indicator
c                   (0=none, 1=static, 2=static and dynamic)
c     beto        - initial stiffness damping factor
c     ddise(ndof) - element nodal incremental displacement vector
c     dise(ndof)  - element total nodal displacement vector
c     vele(ndof)  - element nodal velocity vector
c   Output:
c     ener        - change in element elasto-plastic energy
c     ened        - change in element damping energy
c                   (not computed)
c     enso        - change in element second-order energy
c                   (always zero for this element)
c     relas(ndof) - element static resisting force vector
c     rdamp(ndof) - element damping resisting force vector
c                   (not computed)
c     rinit(ndof) - element initial resisting force vector
c                   (not computed)
c   Modify:
c     kstep       - step number in this segment
c     kst         - stiffness formation code
c                   (1=yes, 0=no)
c ---------------------------------------------------------------------
c  DOUBLE PRECISION
      implicit none
c ---------------------------------------------------------------------
c  ARGUMENT TYPES
      integer   kresis,ksave,kgem,kstep,ndof,kst,kenr
      real*8    ener,ened,enso,beto
      real*8    relas(ndof),rdamp(ndof),rinit(ndof)
      real*8    ddise(ndof),dise(ndof),vele(ndof)
      
c ---------------------------------------------------------------------
c     LABELLED COMMONS
      real*8		E,sigY,Hiso,Hkin
      real*8		ep,alpha,kappa
      real*8		epsP,sigP,tangP
      real*8		tang
      include		'infel00.h'
      
      integer sgn
      real*8	eps,sig,xsi,f,dGamma
      
c     Current strain
      eps = dise(2)-dise(1)    
      
c     Elastic predictor
      sig = E * (eps - ep)
      
c     Stress relative to back stress
      xsi = sig - kappa
      
c     Yield function
      f = dabs(xsi) - (sigY + Hiso*alpha)

c     Inside yield surface
c					 if (f <= 0.0) then
      if (f.le.0.0) then
         tang = E
c     Outside yield surface ... do return mapping
      else
c     Consistency parameter
         dGamma = f / (E+Hiso+Hkin)

c     Normal to yield surface
c         if (xsi <= 0.d0) then
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

c     Populate force vectors
      relas(2) =  sig
      relas(1) = -sig
      
c     Not computed
      rinit(2) = 0.d0
      rinit(1) = 0.d0
      
c     Not computed
      rdamp(2) = 0.d0
      rdamp(1) = 0.d0
      
      ener = 0.5d0*(sig+sigP)*(ddise(2)-ddise(1))
      
c     Not computed
      ened = 0.d0
      enso = 0.d0
      
      epsP = eps
      sigP = sig
      
      END
   

c **********************************************************************
      SUBROUTINE STIF00 (kstt,ktype,ndof,fk)
c *********************************************************************
c  DRAIN-2DX ROTATIONAL SPRING ELEMENT - Bilinear return map
c ---------------------------------------------------------------------
c  PURPOSE
c     The purpose of the subroutine stif00() is to form the stiffness
c     matrix and return it to the variable fk.
c ---------------------------------------------------------------------
c  CALLED FROM:  stifxx
c ---------------------------------------------------------------------
c  ARGUMENTS
c    Input:
c      kstt           - (1=total, 0=change in) stiffness
c      ktype          - 1 = elastic stiffness only,
c                       0 = elastic + geometric stiffness
c                      -1 = geometric stiffness only
c      ndof           - number of element DOFs
c    Output:
c      fk(ndof, ndof) - stiffness matrix
c ---------------------------------------------------------------------
c  DOUBLE PRECISION
      implicit none
c ---------------------------------------------------------------------
c  ARGUMENT TYPES
      integer   kstt,ktype,ndof
      real*8    fk(ndof,ndof)

c ---------------------------------------------------------------------
c  LABELLED COMMONS
      real*8		E,sigY,Hiso,Hkin
      real*8		ep,alpha,kappa
      real*8		epsP,sigP,tangP
      real*8		tang
      include		'infel00.h'

c
c TOTAL STIFFNESS MATRIX
c
      fk(1,1) =  tang
      fk(1,2) = -tang
      fk(2,1) = -tang
      fk(2,2) =  tang
	  
c
c CHANGE IN STIFFNESS MATRIX
c
      if(kstt.EQ.0) then
         fk(1,1) =  tang-tangP
         fk(1,2) = -fk(1,1)
         fk(2,1) =  fk(1,2)
         fk(2,2) =  fk(1,1)
      endif
      
      tang = tangP
c
c NO GEOMETRIC STIFFNESS
c
      if(ktype.EQ.-1) then
         fk(1,1) = 0.d0
         fk(1,2) = 0.d0
         fk(2,1) = 0.d0
         fk(2,2) = 0.d0
      end if
      
      RETURN
      END



