c.... *****************************************************************<-70
      subroutine GETCOMMON(mynh1, mynh3, sumnh, myh)
c .....................................................................
c     getCommon - a fortran subroutine to read info from the feap
c                 common blocks
c
c     written:  fmk
c     created:  03/99
c     revision: A
c .....................................................................

      implicit none

c.... common declarations
      character*4 o,head
      common /bdata/ o,head(20)

      integer         numnp,numel,nummat,nen,neq,ipr
      common /cdata/  numnp,numel,nummat,nen,neq,ipr

      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel

      real *8         bpr, ctan
      common /eltran/ bpr(3), ctan(3)

      integer         ior,iow
      common /iofile/ ior,iow

      integer         nh1,nh2,nh3
      common /hdata/  nh1,nh2,nh3

      real*8  h
      common  h(1000)

c ... subroutine arguments
      integer mynh1, mynh3, sumnh
      real *8 myh(*)

c ... local variables
      integer  i

c ... simply retrieve the variables from the common blocks     
      mynh1 = nh1
      mynh3 = nh3

c ... copy the stuff from the common block to the h array
      do i=1,sumnh
         myh(i) = h(i)
      enddo

      end

