      subroutine elmt04(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw,
     1                  n,nh1,nh2,nh3,h,ctan,iow)

      implicit none

      real*8  d(*), ul(*), xl(*), tl(*), s(*), p(*), h(*), ctan(*)
      integer ix(*), ndf, ndm, nst, isw, nh1, nh2, nh3, iow, n

      if(isw.gt.0) write(*,1000)
1000  format('WARNING: elmt04()-dummy subroutine, no elmt04() linked')
      end

