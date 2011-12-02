      subroutine naccel(n,itr,mvec,tol,u,  f)
************************************************************************
*
*   NACCEL -- Newton iteration accelerator.
*
*   Argument  I/O/M/T  Description
*   --------  -------  -----------
*     N          I     Vector size.
*
*     ITR        I     Newton iteration count.  Initialization is done
*                      for ITR=1 and ITR is ignored thereafter.
*
*     MVEC       I     Maximum number of vectors to use in the
*                      acceleration algorithm.  May change from call
*                      to call but should be greater than 0 and no
*                      more than the internal parameter MAX (=10).
*
*     TOL        I     Tolerance for dropping vectors.  We drop the
*                      pair (z_k,w_k) if the sine of the angle between
*                      w_k and span{w_1, ..., w_(k-1)} is less than TOL.
*
*     U          T     Work space of length at least N*(2*MVEC+2).
*                      It should not be modified between calls.
*
*     F          M     Vector of length N.  On entry, it is the value
*                      of the function f at the current iterate.  It
*                      is overwritten with the accelerated correction.
*                      The unaccelerated correction would simply be f
*                      itself.
*
************************************************************************
      integer max
      parameter(max=10)

      integer n,itr,mvec
      double precision tol,u(n,2,mvec+1),f(n)

      double precision h(max+1,max+1)
      integer nvec,head,next,link(max+1)
      save nvec,head,next,link,h

      double precision c(max),t
      integer jptr,kptr,km1ptr,last,tmp,i,j,k
*     ==================================================================
*      First call: save f and the (unaccelerated) correction.
*     ==================================================================
      if(itr .eq. 1) then
        head = 1
        do 10 j=1,n
         u(j,1,head) = f(j)
         u(j,2,head) = f(j)
   10   continue
        link(1) = 0
        nvec = 1

*       Free storage linked list.
        next = 2
        do 20 k=2,max
   20   link(k) = k + 1
        link(max+1) = 0
        return
        endif
*     ==================================================================
*      Compute w_1.
*     ==================================================================
      do 100 j=1,n
  100 u(j,2,head) = u(j,2,head) - f(j)
      t = 0.0d0
      do 101 j=1,n
  101 t = t + u(j,2,head)**2
      t = 1.0d0/sqrt(t)

*     Normalize w_1 and apply same factor to z_1.
      do 110 j=1,n
       u(j,1,head) = t*u(j,1,head)
       u(j,2,head) = t*u(j,2,head)
  110 continue

*     Update H.
      kptr = link(head)
      do 120 k=2,nvec
       h(1,k) = 0.0d0
       do 121 j=1,n
  121  h(1,k) = h(1,k) + u(j,2,head)*u(j,2,kptr)
       kptr = link(kptr)
  120 continue
*     ==================================================================
*      Compute the Choleski factorization of H.
*     ==================================================================
      k = 2
      h(1,1) = 1.d0
  200 if(k .gt. min(nvec,mvec)) go to 250

      do 210 j=1,k-1
       h(k,j) = h(j,k)
       do 211 i=1,j-1
  211  h(k,j) = h(k,j) - h(k,i)*h(j,i)
       h(k,j) = h(k,j)/h(j,j)
  210 continue

      h(k,k) = 1.d0
      do 220 j=1,k-1
  220 h(k,k) = h(k,k) - h(k,j)**2

      if(h(k,k) .lt. tol**2) then
*       -----------------------------------------------
*        w_k is nearly in span{w_1, ..., w_(k-1)}
*       -----------------------------------------------
*       Remove w_k from linked list.
        km1ptr = head
        do 230 j=2,k-1
  230   km1ptr = link(km1ptr)      
        kptr = link(km1ptr)
        link(km1ptr) = link(kptr)
        nvec = nvec - 1

*       Update free storage list.
        link(kptr) = next
        next = kptr

*       Update H.
        do 240 j=k,nvec
         do 241 i=1,k-1
  241    h(i,j) = h(i,j+1)
         do 242 i=k,j-1
  242    h(i,j) = h(i+1,j+1)
  240   continue
        go to 200

      else
        h(k,k) = sqrt(h(k,k))
        k = k + 1
        go to 200
        endif
*     ------------------------------
*      Retain at most MVEC vectors.
*     ------------------------------
  250 if(nvec .gt. mvec) then
*       truncate the linked list.
        last = head
        do 260 j=2,mvec
  260   last = link(last)
        tmp = link(last)
        link(last) = 0
        last = tmp

*       Update free storage list.
        do 270 j=mvec+2,nvec
  270   last = link(last)
        link(last) = next
        next = tmp
        nvec = mvec
        endif
*     ==================================================================
*      Compute the projection of f onto {w_1, ... , w_nvec}.
*     ==================================================================
      jptr = head
      do 300 j=1,nvec
       c(j) = 0.0d0
       do 301 i=1,n
  301  c(j) = c(j) + f(i)*u(i,2,jptr)
       jptr = link(jptr)
       do 310 i=1,j-1
  310  c(j) = c(j) - h(j,i)*c(i)
       c(j) = c(j)/h(j,j)
  300 continue

      do 320 j=nvec,1,-1
       do 321 i=j+1,nvec
  321  c(j) = c(j) - h(i,j)*c(i)
       c(j) = c(j)/h(j,j)
  320 continue
*     ==================================================================
*      Compute the accelerated correction.
*     ==================================================================
*     Save f for the next call.
      do 410 j=1,n
  410 u(j,2,next) = f(j)

      kptr = head
      do 420 k=1,nvec
       do 421 j=1,n
  421  f(j) = f(j) - c(k)*u(j,2,kptr) + c(k)*u(j,1,kptr)
       kptr = link(kptr)
  420 continue

*     Save the correction for the next call.
      do 430 j=1,n
  430 u(j,1,next) = f(j)
*     ==================================================================
*      Shift the vectors to the right.
*     ==================================================================
      tmp = next
      next = link(tmp)
      link(tmp) = head
      head = tmp

*     Update H.
      do 500 j=nvec,1,-1
      do 500 i=1,j-1
  500 h(i+1,j+1) = h(i,j)

      nvec = nvec + 1

      return
      end
