c.... *****************************************************************<-70
c.... ELEMENT 11: Truss element 
c....             Plasticity with multiple hardening rules
c.... *****************************************************************
      subroutine ELMT11(d,ul,x,ix,t,s,p,ndf,ndm,nst,isw)
c.... *****************************************************************<-70
c.... USER NOTES.......................................................
c.... -----------------------------------------------------------------
c.... material,i
c....     user,11
c....     E,Area,Density,Yield,Q#,a,b,c
c.... -----------------------------------------------------------------
c.... Q# = 1 -> Linear Q Function:  f = |sig| - yield + q
c....           q = -a*alpha
c.... -----------------------------------------------------------------
c.... Q# = 2 -> Dual Exponential:   f = |sig| - yield + q
c....           q = [1 + exp(a*alpha) -2*exp(b*alpha)]*c
c.... -----------------------------------------------------------------
c.... Q# = 3 -> Hyperbolic Secant:  f = |sig| - yield + q
c....           q = [1 - sech(a*alpha)]*b
c.... -----------------------------------------------------------------
c.... Q# = 4 -> ?????????????????:  f = |sig| - yield + q
c....           q = ?
c.... -----------------------------------------------------------------
c.... *****************************************************************
c.... PROGRAM NOTES.................................................... 
c.... d(1)  = Young's Modulus
c.... d(2)  = Cross Sectional Area
c.... d(3)  = Material Density
c.... d(4)  = Initial yield value
c.... d(5)  = Switch for "q" model
c....       = 1 = Linear Hardening          - [uses d(6) only]
c....       = 2 = Dual exponential   "q"    - [uses d(6 to 8)]
c....       = 3 = Hyperbolic secant  "q"    - [uses d(6 to 7)]
c....       = * = **Write your own** "q"    - [uses d(? to ?)]
c.... d(6)  = q-parameter 1 (a)
c.... d(7)  = q-parameter 2 (b)
c.... d(8)  = q-parameter 3 (c)
c....
c.... Notes: To add another q-function:
c....       (1) Add the computation of your q in the trial state 
c....           calculation (two (isw) places)
c....       (2) Add the solver for your residual equation
c....           -Pass back updated consistency parameter - (alpahn)
c....           -Pass back updated stress state          - (sig)
c....           -Pass back the consistent tangent        - (calg)
c.... *****************************************************************
      implicit none
c.... *****************************************************************
c.... Make common declarations
c.... *****************************************************************
      character*4 o,head
      common /bdata/ o,head(20)
      integer         numnp,numel,nummat,nen,neq,ipr
      common /cdata/  numnp,numel,nummat,nen,neq,ipr
      real*8          dm
      integer            n,ma,mct,iel,nel
      common /eldata/ dm,n,ma,mct,iel,nel
      integer         ior,iow
      common /iofile/ ior,iow
      integer         nh1,nh2,nh3
      common /hdata/  nh1,nh2,nh3
      real*8  h
      common  h(1000)
c.... *****************************************************************
c.... Input/Output declarations
c.... *****************************************************************
      integer ix(*)
      integer ndf,ndm,nst,isw
      real*8  d(8),ul(ndf,1),x(ndm,1),t(*),s(nst,nst),p(1)   
c.... *****************************************************************
c.... Local declarations
c.... *****************************************************************
      integer i1,ii,j1,jj,i,j,qmodel
      real*8 xx(3),dx(3),db(3)
      real*8 xlen,deps,alphan,sign,sigtr,asigtr,sig
      real*8 ftr,calg,v,dgam,xll,qval,qgam
      save 
c.... *****************************************************************
c.... ISW Routing
c.... *****************************************************************
      go to (1,2,3,4,5,4,2,2,2,2),isw
c.... *****************************************************************
c.... ISW 0: Greet User!
c.... *****************************************************************
      if(isw.eq.0 .and. ior.lt.0) then
        write(*,*) ' Hi! I am Element 11, therefore I exist.'
        write(*,*) ' Elastic-Plastic Truss with "q" function'
      endif
      return
c.... *****************************************************************
c.... ISW 1: Get material information from input file
c.... *****************************************************************
1     call dinput(d,8)
      write(iow,2000) d(1),d(2),d(3),d(4),d(5),d(6),d(7),d(8)
      if(ior.lt.0) then
         write(*,2000) d(1),d(2),d(3),d(4),d(5),d(6),d(7),d(8)
      endif
      d(3) = d(3)*d(2)
      nh1 = 2         
      qmodel = d(5)
      call pzero (xx,3)
c.... *****************************************************************
c.... Plotting Order
c.... *****************************************************************
c      inord(iel)    = 3
c      ipord( 1,iel) = 1
c      ipord( 2,iel) = 2
c      ipord( 3,iel) = 1
c.... *****************************************************************
c.... ISW 2: Nothing..........
c.... *****************************************************************
2     return
c.... *****************************************************************
c.... ISW 3: Compute element stiffness and internal force arrays
c.... *****************************************************************
3     continue
c.... *****************************************************************
c.... Initialize and compute the strain increment
c.... *****************************************************************
      qmodel  = d(5)
      xlen    = 0.0
      deps    = 0.0
      alphan  = h(nh1)
      sign    = h(nh1+1)
      do 31 i = 1,ndm
        dx(i) = x(i,2) - x(i,1)
        xlen  = xlen + dx(i)**2
	deps  = deps + dx(i)*(ul(i,2+nen)-ul(i,1+nen))
31    continue
      deps    = deps/xlen
c.... *****************************************************************
c.... Compute the trial stress and check yield surface violation
c.... *****************************************************************
      sigtr = sign + d(1)*deps
      asigtr= dabs(sigtr)
      if(qmodel.eq.1) then
         qval = -d(6)*alphan
      elseif(qmodel.eq.2) then
         qval = (1 + dexp(d(6)*alphan) - 2.d0*dexp(d(7)*alphan))*d(8)
      elseif(qmodel.eq.3) then
         qval = (1 - 1/dcosh(d(6)*alphan))*d(7)
      else
         write(*,*) 'Hello!!! q function',qmodel,' not defined...yet.'
         stop
      endif
      ftr   = asigtr - d(4) + qval
      if(ftr .le. 0.d0) then
c....   ***************************************************************
c....   If elastic update stress and go
c....   ***************************************************************
         calg = d(1)
         sig = sigtr
      else
c....   ***************************************************************
c....   Else check hardening type, and solve the plasticity equations
c....   ***************************************************************
         if(qmodel.eq.1) then
            v = d(1)+d(6)
            calg = d(1)*d(6)/v
            dgam = ftr/v - 1.d-10
            alphan = alphan + dgam
            sig = sigtr - dgam*d(1)*sigtr/asigtr
         elseif(qmodel.eq.2) then
            call qdualexp(d,ftr,sigtr,alphan,sig,calg)
         elseif(qmodel.eq.3) then
            call qsech(d,ftr,sigtr,alphan,sig,calg)
         else
            write(*,*) ' Hardening model undefined- stopping! ' 
            stop
         endif
      endif
c.... *****************************************************************
c.... Update the history variables
c.... *****************************************************************
c      if(aratfl) c = d(1)
      h(nh2)   = alphan
      h(nh2+1) = sig
c.... *****************************************************************
c.... Form the residual
c.... *****************************************************************
      xll = dsqrt(xlen)
      sig = sig*d(2)/xll
      do 39 i = 1,ndf
       p(i) = sig*dx(i)
       p(i+ndf) = -p(i)
39    continue
      xlen  = xlen*xll
      do 32 i = 1,ndm
	db(i) = d(2)*calg*dx(i)
        dx(i) = dx(i)/xlen
32    continue
c.... *****************************************************************
c.... Compute the stiffness
c.... *****************************************************************
      i1 = 0
      do 37 ii = 1,2
        j1 = i1
        do 36 jj = ii,2
          do 34 i = 1,ndm
            do 33 j = 1,ndm
              s(i+i1,j+j1) = db(i)*dx(j)
33          continue
34        continue
          j1 = j1 + ndf
c.... *****************************************************************
c.... Change sign of coordinate differences for other end
c.... *****************************************************************
          do 35 j = 1,ndm
            dx(j) = -dx(j)
35        continue
36      continue
          i1 = i1 + ndf
37    continue
c.... *****************************************************************
c.... Fill in using symmetry
c.... *****************************************************************
      do 38 i = 1,ndm
      do 38 j = 1,ndm
        s(i+ndf,j) = s(j,i+ndf)
38    continue
      return
c.... *****************************************************************
c.... ISW 4: Compute current stress and strain
c.... *****************************************************************
4     continue
c.... *****************************************************************
c.... Initialize and compute the strain increment
c.... *****************************************************************
      qmodel  = d(5)
      xlen    = 0.0d0
      deps    = 0.0d0
      alphan  = h(nh1)
      sign    = h(nh1+1)
      do 41 i = 1,ndm
        dx(i) = x(i,2) - x(i,1)
        xlen  = xlen + dx(i)**2
	deps  = deps + dx(i)*(ul(i,2+nen)-ul(i,1+nen))
        xx(i) = (x(i,2) + x(i,1))/2.
41    continue
      deps  = deps/xlen
c.... *****************************************************************
c.... Compute the trial stress and check yield surface violation
c.... *****************************************************************
      sigtr = sign + d(1)*deps
      asigtr= dabs(sigtr)
      if(qmodel.eq.1) then
         qval = -d(6)*alphan
      elseif(qmodel.eq.2) then
         qval = (1 + dexp(d(6)*alphan) - 2.d0*dexp(d(7)*alphan))*d(8)
      elseif(qmodel.eq.3) then
         qval = (1 - 1/dcosh(d(6)*alphan))*d(7)
      else
         write(*,*) 'Hello!!! q function',qmodel,' not defined...yet.'
         stop
      endif
      ftr   = asigtr - d(4) + qval
      if(ftr .le. 0.d0) then
c....   ***************************************************************
c....   If elastic update stress and go
c....   ***************************************************************
         calg = d(1)
         sig = sigtr
      else
c....   ***************************************************************
c....   Else check hardening type, and solve the plasticity equations
c....   ***************************************************************
         if(qmodel.eq.1) then
            v = d(1)+d(6)
            calg = d(1)*d(6)/v
            dgam = ftr/v - 1.d-10
            alphan = alphan + dgam
            sig = sigtr - dgam*d(1)*sigtr/asigtr
         elseif(qmodel.eq.2) then
            call qdualexp(d,ftr,sigtr,alphan,sig,calg)
         elseif(qmodel.eq.3) then
            call qsech(d,ftr,sigtr,alphan,sig,calg)
         else
            write(*,*) ' Hardening model undefined- stopping! ' 
            stop
         endif
      endif
c.... *****************************************************************
c.... Update the history variables
c.... *****************************************************************
      h(nh2)   = alphan
      h(nh2+1) = sig
c.... *****************************************************************
c.... Write out stress, etc, to file
c.... *****************************************************************
      if(isw.eq.4) then
        mct = mct - 1
        if(mct.le.0) then
          write(iow,2001) o,head
          if(ior.lt.0) write(*,2001) o,head
          mct = 50
        endif
        write(iow,2002) n,ma,xx,sig,alphan
        if(ior.lt.0) write(*,2002) n,ma,xx,sig,alphan
      elseif(isw.eq.6) then
c....   ***************************************************************
c....   ISW 6: Form the residual
c....   ***************************************************************
        sig = sig*d(2)/dsqrt(xlen)
        do 42 i = 1,ndf
          p(i) = sig*dx(i)
	  p(i+ndf) = -p(i)
42      continue
      endif
      return
c.... *****************************************************************
c.... ISW 5: Compute mass distribution
c.... *****************************************************************
5     continue
      xlen = 0.0
      do 51 i = 1,ndm
        xlen = xlen + (x(i,2)-x(i,1))**2
51    continue
      s(1,1) = d(3)*sqrt(xlen)/3.0
      do 52 i = 1,ndm
c.... *****************************************************************
c.... Consistent mass
c.... *****************************************************************
        s(i,i) = s(1,1)
        s(i+ndf,i+ndf) = s(1,1)
        s(i+ndf,i) = s(1,1)/2.
        s(i,i+ndf) = s(1,1)/2.
c.... *****************************************************************
c.... Lumped mass
c.... *****************************************************************
        p(i) = s(1,1)*1.5
        p(i+ndf) = p(i)
52    continue
      return
c.... *****************************************************************
c.... END OF MAIN
c.... *****************************************************************
c.... *****************************************************************
c.... Format statements
c.... *****************************************************************
2000  format(9x,'truss element'//10x,
     1 'Modulus           =',e15.5/10x,
     1 'Area              =',e15.5/10x,
     1 'Density           =',e15.5/10x,
     1 'Yield Stress      =',e15.5/10x,
     1 'q model number    =',e15.5/10x,
     1 'q-parameter 1     =',e15.5/10x,
     1 'q-parameter 2     =',e15.5/10x,
     1 'q-parameter 3     =',e15.5)
2001  format(a1,20a4//9x,'truss element'//' elmt matl    ',
     1 '1-coord    2-coord    3-coord',5x,'stress',8x,'alpha')
2002  format(2i5,1p3e11.3,1p2e14.5)
      end


c.... *****************************************************************
c.... *****************************************************************
c.... QDUALEXP SUBROUTINE BELOW
c.... *****************************************************************
c.... *****************************************************************
      subroutine qdualexp(d,ftr,sigtr,alphan,sig,calg)
      implicit none
c.... *****************************************************************
c.... Make common declarations
c.... *****************************************************************
      logical debug
      common /debugs/ debug 
c.... *****************************************************************
c.... Input/Output declarations
c.... *****************************************************************
      real*8  d(1)
      real*8  ftr,alphan,sigtr,sig,calg
c.... *****************************************************************
c.... Local declarations
c.... *****************************************************************
      integer iter,itmax,file
      real*8  tol,dgam,gfun,qo,qn,dqn,dgfun 
c.... *****************************************************************
c.... Set the tolerance for convergence
c.... *****************************************************************
      tol = dmax1(1.d0,ftr)
      tol = dexp((dlog10(tol) - 12.d0)*dlog(10.d0)) 
c.... *****************************************************************
c.... Initialize
c.... *****************************************************************
      file    = 2
      iter    = 0
      itmax   = 25
      dgam    = 0
      qo      = (1 + dexp(d(6)*alphan) - 2*dexp(d(7)*alphan))*d(8)
c.... *****************************************************************
c.... Start iteration 
c.... *****************************************************************
 111  iter  = iter +1
      if(debug) then
         qn = alphan + dgam
         write(file,*) 'gfun(',iter,') = ',gfun,';'
         write(file,*) 'alpha(',iter,') = ',qn,';'
         if(iter.eq.itmax) write(*,*) 'plot(gfun)'
      endif
      if(iter.eq.itmax) then
         write(*,*) 'Plasticity not converged.'
         write(*,*) 'qdualexp subroutine (111)'
         write(*,*) 'STOPPING'
         stop
      endif
      qn    = 1 + dexp(d(6)*(alphan+dgam)) - 2*dexp(d(7)*(alphan+dgam))
      qn    = qn*d(8)
      dqn   = d(6)*dexp(d(6)*(alphan+dgam))
      dqn   = dqn - 2*d(7)*dexp(d(7)*(alphan+dgam))
      dqn   = dqn*d(8)
      gfun  = ftr - d(1)*dgam + qn - qo
      dgfun = -d(1) + dqn
      if(dabs(gfun).le.tol) goto 222
      if(dabs(dgfun).le.1.d-16) then
         write(*,*) 'Derivative of g function almost zero!'
         write(*,*) 'Check validity of input parameters. '
         write(*,*) 'STOPPING'
         stop
      endif
      dgam  = dgam - gfun/dgfun
      goto 111
 222  continue
      dgam   = dgam - 1.d-10
      alphan = alphan + dgam
      sig    = sigtr - d(1)*dgam*sigtr/dabs(sigtr)
      calg   = dqn*d(1)/(dqn - d(1))
      return
      end
c.... *****************************************************************
c.... *****************************************************************
c.... **END QDUALEXP***************************************************
c.... *****************************************************************
c.... *****************************************************************


c.... *****************************************************************
c.... *****************************************************************
c.... QSECH SUBROUTINE BELOW
c.... *****************************************************************
c.... *****************************************************************
      subroutine qsech(d,ftr,sigtr,alphan,sig,calg)
      implicit none
c.... *****************************************************************
c.... Make common declarations
c.... *****************************************************************
      logical debug
      common /debugs/ debug 
c.... *****************************************************************
c.... Input/Output declarations
c.... *****************************************************************
      real*8  d(1)
      real*8  ftr,alphan,sigtr,sig,calg
c.... *****************************************************************
c.... Local declarations
c.... *****************************************************************
      integer iter,itmax,file
      real*8  tol,dgam,gfun,qo,qn,dqn,dgfun 
c.... *****************************************************************
c.... Set the tolerance for convergence
c.... *****************************************************************
      tol = dmax1(1.d0,ftr)
      tol = dexp((dlog10(tol) - 12.d0)*dlog(10.d0)) 
c.... *****************************************************************
c.... Initialize
c.... *****************************************************************
      file    = 3
      iter    = 0
      itmax   = 25
      dgam    = 0
      qo      = (1 - 1/dcosh(alphan*d(6)))*d(7)
c.... *****************************************************************
c.... Start iteration 
c.... *****************************************************************
 111  iter  = iter + 1
      if(debug) then
         qn = alphan + dgam
         write(file,*) 'gfun(',iter,') = ',gfun,';'
         write(file,*) 'alpha(',iter,') = ',qn,';'
         if(iter.eq.itmax) write(*,*) 'plot(gfun)'
      endif
      if(iter.eq.itmax) then
         write(*,*) 'Plasticity not converged.'
         write(*,*) 'qdualexp subroutine (111)'
         write(*,*) 'STOPPING'
         stop
      endif
      qn    = (1 - 1/dcosh(d(6)*(alphan+dgam)))*d(7)
      dqn   = d(6)*d(7)*dtanh(d(6)*(alphan+dgam))
      dqn   = dqn/dcosh(d(6)*(alphan+dgam))
      gfun  = ftr - d(1)*dgam + qn - qo
      dgfun = -d(1) + dqn
      if(dabs(gfun).le.tol) goto 222
      if(dabs(dgfun).le.1.d-16) then
         write(*,*) 'Derivative of g function almost zero!'
         write(*,*) 'Check validity of input parameters. '
         write(*,*) 'STOPPING'
         stop
      endif
      dgam  = dgam - gfun/dgfun
      goto 111
 222  continue
      dgam   = dgam - 1.d-10
      alphan = alphan + dgam
      sig    = sigtr - d(1)*dgam*sigtr/dabs(sigtr)
      calg   = dqn*d(1)/(dqn - d(1))
      return
      end
c.... *****************************************************************
c.... *****************************************************************
c.... **END QSECH******************************************************
c.... *****************************************************************
c.... *****************************************************************















