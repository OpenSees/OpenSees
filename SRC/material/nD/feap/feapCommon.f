      subroutine FEAPCOMMON (mydt, myniter)

c .....................................................................
c     feapCommon - a fortran subroutine to fill the named common block
c              used by Feap material models
c     ADD MORE ARGUMENTS AS NEEDED!!!
c
c     written:  MHS
c     created:  June 2001
c .....................................................................
	
      implicit none

      real*8         mydt
      integer        myniter
      
      integer         nstep,niter,naugm, titer,taugm, iaugm, iform
      common /counts/ nstep,niter,naugm, titer,taugm, iaugm, iform
      
      real*8          ttim,dt,c1,c2,c3,c4,c5, chi
      common /tdata/  ttim,dt,c1,c2,c3,c4,c5, chi
      
      dt    = mydt
      niter = myniter
      
c     ADD MORE ASSIGNMENTS AS NEEDED!!!
      
      end
