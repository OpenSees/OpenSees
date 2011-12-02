      subroutine ELMT05(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

      implicit none

      real*8  d(*), ul(*), xl(*), tl(*), s(*), p(*)
      integer ix(*), ndf, ndm, nst, isw

      if(isw.gt.0) write(*,1000)
1000  format('WARNING: elmt05()-dummy subroutine, no elmt05() linked')
      end

