c.... *****************************************************************<-70
      subroutine FILLCOMMON(mynen, mydm, myn, myior, myiow, mynh1, 
     1                      mynh2, mynh3, sumnh, myh, myctan, nrcount)
c .....................................................................
c     fillCommon - a fortran subroutine to fill the named common block
c                  used by the elements .. NOT YET TAKING ALL ARGS
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

      real*8          tol,rnmax,shift
      logical                         linear,shflg
      common /rdata/  tol,rnmax,shift,linear,shflg

      real*8  hr
      common  hr(10000)


c ... subroutine arguments
c.... *****************************************************************<-70
      integer mynen, myn, myior, myiow, mynh1, mynh2, mynh3, sumnh
      integer nrcount
      real *8 mydm, myh(*), myctan(*)

c ... local variables
      integer  i

c ... simply set the variables in the common blocks      
      nen = mynen
      n = myn
      ior = myior
      iow = myiow
      nh1 = mynh1
      nh2 = mynh2
      nh3 = mynh3
      rnmax = nrcount

c ... copy the stuff in the h array to the common block
      if (sumnh.gt.10000) then
         write(*,*)'fillCommon.f - allocated common block of'
         write(*,*)'needs to be of size: ',sumnh
         stop
      endif
      

      do i=1,sumnh
         hr(i) = myh(i)
      enddo

      ctan(1) = myctan(1)
      ctan(2) = myctan(2)
      ctan(3) = myctan(3)

      end

