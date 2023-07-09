      subroutine PD (d,hstvp,hstv,epsp,sigp,deps,str,dd,ist)

c
c       One-D version of (ElastoVisco) Plastic-Damage Model
c
c    % Modified for OpenSees: June 2006
c    % Original Coded: July 1998
c    % By Jeeho Lee (jeeholee@berkeley.edu)
c
c    # Input:  d[*],delt,eps,stif_typ
c    # Output: str,dd
c    # History Variables (size=11): hstvp, hstv
c     ->ieps,peps,kapa[2],cohn[2],deg,vdeg,dplas,oldeps,phibound
c
c    ! ist = 1: dd is tangent stiffness in any case
c


      integer ist
      real*8 d(10),hstvp(11),hstv(11),delt,xl(2,2),eps,str,dd
      real*8 oldeps,phibound
	  real*8 epsp,sigp,deps

c  Local Variables ----------------------------------------------

      integer index,iter,maxitr,crmode,flag
      real*8 cohn(2),kapa(2),peps,ieps,dplas,dplas1,matpara(4)
      real*8 trstr,estr,deleps,delps
      real*8 kp,sign,chleng
      real*8 lam,e,fenergy
      real*8 tol,tol2,cmax,kmax
      real*8 resf,ktcrit,fbody
      real*8 deg,vdeg,degstr,d2_eps,fstr,fkp,ck
      real*8 mu,viscom,viscot
      real*8 temp


c set numerical parameters ------------------------------------------

      tol = 10.d-12
      tol2 = 10.d-8
	  maxitr = 10

c temporary setting ++++++++++-----------------------


	  
	  eps = epsp + deps
	  
	  xl(1,1) = 0.
	  xl(1,2) = 1.
	  xl(2,1) = 0.
	  xl(2,2) = 0.
	  
	  delt = 1.d0
	  

c initialize the material properties & variables --------------------

      e = d(1)

      ktcrit = d(7)
      mu = d(8)*delt


c  retrieve history variables ---------------------------------------
 
      ieps = hstvp(1)
      peps = hstvp(2)
      kapa(1) = hstvp(3)
      kapa(2) = hstvp(4)
      cohn(1) = hstvp(5)
      cohn(2) = hstvp(6)
      deg = hstvp(7)
      vdeg = hstvp(8)
      dplas = hstvp(9)
      oldeps = hstvp(10)
      phibound = hstvp(11)
  
c  ------------------------------------------------------------------

      crmode = 0

c      if (delt .le. 0.d0) then
	     delt = 1.d0
c      endif
      viscom = mu/(mu+delt)
      viscot = delt/(mu+delt)
      
      if ((kapa(1)+kapa(2)) .le. 0.d0) then
        cohn(1) = d(4)
        cohn(2) = -d(5)
      endif


         
c  compute trial effective stresses ---------------------------------

      deleps = eps - peps
      trstr = e*deleps

c  check whether tensile or compressive state

      if (trstr .ge. 0.d0) then
        sign = 1.d0
        index = 1
        kp = kapa(1)
        dplas1 = 1.d0 -dplas    
      else
        sign = -1.d0
        index = 2
        kp = kapa(2)
        dplas1 = 1.d0
      endif

      flag = 0
	  
	  
	  call setpara(d,matpara)
  
      call yield1(index,cohn,trstr,resf,temp)

      cmax = dsqrt(d(4)*d(4) + d(5)*d(5))

      if (resf .lt. tol*cmax) then

        if (kp.gt.0.d0 .and. index.eq.2) then
        deps = eps - oldeps
        if (deps.lt.0.d0) then
C+++++++++++++++++++++++++ Reloading +++++++++++++++++++++++++++++++
          chleng = dsqrt((xl(1,1)-xl(1,2))**2+(xl(2,1)-xl(2,2))**2) 
       
	      call reloading(chleng,kp,d,matpara,eps,deps,peps,phibound
     & 	                  ,cohn,tol,maxitr)
          kapa(2) = kp
        else
   
          call unloading(d,eps,deps,kp,cohn,peps,tol,maxitr)
          phibound = resf/cohn(2) + 1.d0
        endif
        endif

        call degrad1(0,index,d,matpara,kapa(1),kapa(2),ck,deg,degstr)

        ieps = viscot*peps + viscom*ieps
        deleps = eps - ieps
 
        vdeg = viscom*vdeg + viscot*deg
        
        call elastan1(e,vdeg,dplas1,dd)
        
        str = (1.d0-vdeg)*dplas1*e*deleps
    
      else
   
c  plastic loading state  ---------------------------------------------
        
	chleng = dsqrt((xl(1,1)-xl(1,2))**2+(xl(2,1)-xl(2,2))**2)

        if ((kp.gt.ktcrit).and.(index.eq.1)) then

          delps = 0.d0
          call crstr1(e,cohn,trstr,dplas,dplas1,d2_eps)        
          crmode = 1
c          write(*,*) '!!! crmode ON'

        else
          call plasto1(d,matpara,index,sign,chleng,eps,trstr
     &             ,lam,kp,cohn,fenergy,fstr,fkp,ck,dplas1,tol,maxitr)

          delps = lam*sign

          call algotan1(e,sign,lam,fenergy,fstr,fkp,ck,dd)

        endif

        peps = peps + delps
        estr = e*(eps-peps)

        if (index .eq. 1) then
          kapa(1) = kp
        else
          kapa(2) = kp
        endif
      
        call degrad1(1,index,d,matpara,kapa(1),kapa(2),ck,deg,degstr)

        trstr = dplas1*e*(eps - ieps)

        ieps = viscot*peps + viscom*ieps
        deleps = eps - ieps
        
        vdeg = viscom*vdeg + viscot*deg

        str = (1.d0-vdeg)*dplas1*e*deleps
        
        call vdtan1(crmode,index,d,dplas1,mu,delt
     &             ,trstr,estr,vdeg,degstr,d2_eps,dd)
 
      endif



c  Compute secant tangent stiffness (only for elastic part)

c  Post-processing part

c  save history variables ------------------------------------
 
      hstv(1) = ieps
      hstv(2) = peps
      hstv(3) = kapa(1)
      hstv(4) = kapa(2)
      hstv(5) = cohn(1)
      hstv(6) = cohn(2)
      hstv(7) = deg
      hstv(8) = vdeg
      hstv(9) = dplas
      hstv(10) = eps
      hstv(11) = phibound
  
c  ------------------------------------------------------------
      return
      end


      
	  
 
	  
	  
      subroutine setpara(d,matpara)

      real*8 d(*), matpara(4)
      real*8 temp


      deg_para1 = 0.7d0
	  deg_para2 = 0.5d0
	  deg_para3 = 1.0d0
c ----- compute the parameters for the material model

      temp = d(6)/d(5)
      matpara(1) = 2.d0*temp-1.d0+2.d0*dsqrt(temp*temp-temp)

      temp = (1.d0-deg_para1)**(1.d0/deg_para3)
      matpara(2) = deg_para3
      matpara(3) = dlog(1.d0-deg_para2)/dlog(.5d0/matpara(1) + .5d0)
      matpara(4) = (temp - 0.5d0)/(temp*temp - temp)

      return
      end


	  
	  



      subroutine crstr1(e,cohn,trstr,dplas,dplas1,d2_eps)

      real*8 e,trstr,cohn(*),dplas,dplas1,d2_eps

      
      dplas1 = dplas1*cohn(1)/trstr
      dplas = 1.d0 - dplas1
      d2_eps = -e*dplas/trstr

      return
      end






      subroutine plasto1(d,matpara,index,sign,chleng,eps,trstr
     &           ,lam,kp,cohn,fenergy,fstr,fkp,ck,dplas1,toler,maxitr)



      integer index

      integer iter,maxitr,switch

      real*8 d(*),matpara(*),sign,kp,cohn(2),eps,trstr
      real*8 f,fkp,lam,chleng,dplas1

      real*8 kpn,toler,e
      real*8 resq,error,q_kapa,q_lam,lamkp,fstr,fenergy,ck

      real*8 temp
      

         
c plastic or viscoplastic loading case ---------------------------------


      e = dplas1*d(1)
      
c initial setting --------------------------------------------------

      iter = 0
      switch = 1
      lam = 0.d0

      if (index .eq. 1) then
        fenergy = d(2)/chleng
      else
        fenergy = d(3)/chleng
      endif

      kpn = kp

c iteration for damage evolution ------------------------------------
      
100   continue

      iter = iter + 1

      call damg1(1,index,d,matpara,kp,cohn,fstr,fkp,ck)

      call coml1(index,e,trstr,cohn,ck,lam,lamkp)

c   check the convergence of the damage evolution eqn

      resq = kpn - kp + sign*lam*fstr/fenergy
      error = dabs(resq)


c      write(*,*) '#################################'
c      write(*,*) 'toler =',toler
c      write(*,*) 'error =', error
c      write(*,*) 'kp =',kp

      if (error .gt. toler) then

        if (iter .gt. maxitr) then
	  
	  write(*,*) 'toler =',toler
	  write(*,*) 'error =', error
	  write(*,*) 'kp =',kp
	  error = error/0.d0
          stop 'VEPD_2D: exceed the maximum iteration (iter)!'
        endif

        q_lam = sign*fstr/fenergy
        q_kapa = -1.d0 + sign*lam*fkp/fenergy
        q_kapa = q_kapa + q_lam*lamkp 
        resq = resq/q_kapa

        kp = kp - resq

	temp = 1.d0 - toler
        if (kp.lt.kpn) then 
          kp = kpn
        elseif (kp.gt.temp) then
          kp = temp 
	  switch = -1
	endif

        goto 100                 ! ---> go back to 100 

      endif
        
c     end of iteration loop ----------------------------------------
      
      
      return
      end





      subroutine algotan1(e,sign,lam,fenergy,fstr,fkp,ck,dd)
        

      real*8 sign,e,lam,fenergy,fstr,fkp,ck,dd
      
      real*8 omega
     
c   construct the algorithmic tangent stiffness: dd

      omega = fstr*ck/(fenergy - sign*lam*fkp)
      dd = e*omega/(omega - e)

      return
      end








      subroutine vdtan1(crmode,index,d,dplas1,mu,delt
     &                 ,tstr,estr,vdeg,degstr,d2_eps,dd)


      integer index,crmode

      real*8 d(*),dplas1,tstr,estr,lam,vdeg,degstr,d2_eps,dd

      real*8 e,mu,delt
   
      real*8 temp,temp2,vdeg1


      e = d(1)
      vdeg1 = 1.d0 - vdeg
c --- Modify Inviscid Algorithmic Tangent for Large Crack Mode ---

      if (crmode .ge. 1) then
        dd = vdeg1*(e*dplas1 - estr*d2_eps)

      else

        temp = mu/(mu+delt)
        temp2 = delt/(mu+delt)

        tstr = temp*tstr + temp2**dplas1*estr
       
        dd = temp*vdeg1*dplas1*e 
     &      + temp2*(vdeg1 - tstr*degstr)*dd

      endif

      return
      end






      subroutine elastan1(e,vdeg,dplas1,dd)


      real*8 e,vdeg,dd,dplas1
      
      
c    compute the elastic tangent stiffness


      dd =(1.d0-vdeg)*dplas1*e
 

      return
      end
      





      subroutine coml1(index,e,trstr,cohn,ck,lam,lamkp)
      

      integer index

      real*8 e,trstr,cohn(2),ck
      real*8 lam,lamkp
      
      real*8 temp

    
     
      if (index .eq. 1) then
        lam = (trstr - cohn(1))/e
        lamkp = -ck/e
      else
        lam = -(trstr + cohn(2))/e
        lamkp = -ck/e
      endif
      
      return
      end






      subroutine yield1(index,cohn,str,resf,fb)

      integer index
      real*8 cohn(2),resf,fb,str
       

      if (index .eq. 1) then
        fb = cohn(2)*str/cohn(1)
      else
        fb = -str
      endif

      resf = fb - cohn(2)

      
      return
      end






      subroutine damg1(flag,index,d,matpara,kapa,cohn,fstr,fkp,ck)
      
      integer flag,index

      real*8 d(*),kapa,matpara(*)
      real*8 cohn(2),fkp,ck
      
      real*8 pht,phc,fstr
      real*8 ty,cy,at,ac,rt,rc,rpht,rphc

      real*8 temp


      ac = matpara(1)
      rt = matpara(2)
      rc = matpara(3)
      at = matpara(4)

      if (index .eq. 1) then

C--- Tensile Damage

        ty = d(4)
 
        pht = 1.d0+at*(2.d0+at)*kapa
        rpht = dsqrt(pht)
        fstr = ty*((1.d0+at)*dsqrt(pht)-pht)/at

  
c----------- compute cohesions -----------------------------

          cohn(1) = ty*rpht*((1.d0+at-rpht)/at)**(1.d0-rt)
 
c----------- compute derivatives ---------------------------

          fkp = ty*(2.d0+at)*((1.d0+at)/(2.d0*dsqrt(pht))-1.d0)

          temp = (1.d0 + at - rpht)/at
          ck = .5d0*ty*at*(2.d0+at)*(temp**(1.d0-rt)/rpht -
     &           (1.d0 - rt)*temp**(-rt)/at)
  

C--- Compressive Damage

      else

        cy = d(5)

        phc = 1.d0+ac*(2.d0+ac)*kapa
        rphc = dsqrt(phc)
        fstr = cy*((1.d0+ac)*dsqrt(phc)-phc)/ac

   
c----------- compute cohesions -----------------------------

          cohn(2) = -cy*rphc*((1.d0+ac-rphc)/ac)**(1.d0-rc)
  
c----------- compute derivatives ---------------------------

          fkp = cy*(2.d0+ac)*((1.d0+ac)/(2.d0*dsqrt(phc))-1.d0)

          temp = (1.d0 + ac - rphc)/ac
          ck = -.5d0*cy*ac*(2.d0+ac)*(temp**(1.d0-rc)/rphc -
     &             (1.d0 - rc)*temp**(-rc)/ac)
          
 
       endif

      return
      end



      


      subroutine degrad1(flag,index,d,matpara,kapa1,kapa2,ck,deg,degstr)


      integer flag,index
      
      real*8 d(*),ck,kapa1,kapa2,deg,degkp,degstr,matpara(*)

      real*8 s,dt,dc,rt,rc,ac,at,phi,phi2,temp,temp2,cmax,dfac


      ac = matpara(1)
      rt = matpara(2)
      rc = matpara(3)
      at = matpara(4)
     
      dfac = 1.d0

      phi = dsqrt(1.d0 + at*(2.d0 + at)*kapa1)
      temp = (1.d0 + at - phi)/at
      dt = 1.d0 - temp**rt

      phi2 =  dsqrt(1.d0 + ac*(2.d0 + ac)*kapa2)
      temp2 = (1.d0 + ac - phi2)/ac
      dc = 1.d0 - temp2**rc

      if (index .eq. 1) then
        s = 1.d0
      else
        s = 1.d0 - dfac
      endif

      deg = 1.d0 - (1.d0-s*dt)*(1.d0-dc)

      if (flag. ge. 1) then

       if (index .eq. 1) then
        degkp = .5d0*s*(2.d0+at)*(1.d0-dc)*rt*temp**(rt-1.d0)/phi
        degkp = degkp/ck
        degstr = (1.d0-dc)*degkp
       else
        degkp = .5d0*(2.d0+ac)*(1.d0-s*dt)*rc*temp2**(rc-1.d0)/phi2
        degkp = -degkp/ck
        degstr = (1.d0-s*dt)*degkp
       endif  
   
      endif
      
      return
      end







      subroutine reloading(chleng,kp,d,matpara,eps,deps,peps,phib,cohn
     &                     ,toler,maxitr)


      integer index

      integer iter,maxitr,switch

      real*8 d(*),kp,cohn(2),eps,deps,peps,chleng,phib,matpara(*)

      real*8 f,fkp,e,ar,br,cr
      real*8 pepsn,estr,kpn,phi,resf,cmax
      real*8 toler,resq,error,q_kapa,ek,fstr,fenergy,ck

      real*8 temp
      

         
c plastic or viscoplastic loading case ---------------------------------
 
      e = d(1)
      ar = 0.d0
      br = 1.d0
      cr = 1.d0 
      if (br .gt. 1.d0) then
        br = 1.d0
      endif

c initial setting --------------------------------------------------
    
      index = 2
      iter = 0
      switch = 1

      fenergy = d(3)/chleng
      fenergy = fenergy/cr

      pepsn = peps
      kpn = kp

 
c iteration for damage evolution ------------------------------------

100   continue

      iter = iter + 1
      
      estr = e*(eps - peps)

      call damg1(1,index,d,matpara,kp,cohn,fstr,fkp,ck)          
      call yield1(index,cohn,estr,resf,temp)

c   check the convergence of the damage evolution eqn

      phi = resf/cohn(2) + 1.d0

      if (phi .lt. 0.d0) then
        stop 'RELOADING: Negative phi!'
      endif

      resq = -peps + pepsn + ar*(phi-phib)*deps
      error = dabs(resq)

      if (error .gt. toler) then

        if (iter .gt. maxitr) then
	  
	  write(*,*) 'toler =',toler
	  write(*,*) 'error =', error
	  write(*,*) 'kp =',kp
	  error = error/0.d0
          stop 'RELOADING: exceed the maximum iteration (iter)!'
        endif

        ek = (fenergy/(1.d0-br) + (peps-pepsn)*fkp)/fstr
        q_kapa = -ek + ar*deps*(e*ek + estr*ck/cohn(2))/cohn(2)

        resq = resq/q_kapa

        kp = kp - resq

	temp = 1.d0 - toler
        if (kp.lt.kpn) then 
          kp = kpn
        elseif (kp.gt.temp) then
          kp = temp 
	  switch = -1
	endif

        peps = pepsn + (kp-kpn)*fenergy/(fstr*(1.d0-br))

        goto 100                 ! ---> go back to 100 

      endif
        
c     end of iteration loop ----------------------------------------
 

      return
      end






      subroutine unloading(d,eps,deps,kp,cohn,peps,toler,maxitr)

      integer iter,maxitr,switch,index

      real*8 d(*),cohn(2),eps,deps,peps,kp

      real*8 f,fkp,e,ar,br
      real*8 resf,q_epsr,estr,pepsn
      real*8 toler,resq,error,q_kapa,ek,fstr,ck,temp
      

         
c plastic or viscoplastic loading case ---------------------------------
 
      e = d(1)
      ar = 0.d0
      br = 1.d0

c initial setting ---------------------------------------------------
      iter = 0
      switch = 1
      index = 2
      pepsn = peps
 
c computation for damage evolution ------------------------------------

      br = br

100   continue

      iter = iter + 1

      estr = e*(eps - peps)
               
      call yield1(index,cohn,estr,resf,temp)

c   check the convergence of the damage evolution eqn

      resq = -peps + pepsn - ar*deps*resf/cohn(2)
      error = dabs(resq)

      if (error .gt. toler) then

        if (iter .gt. maxitr) then
	  
	  write(*,*) 'toler =',toler
	  write(*,*) 'error =', error
	  error = error/0.d0
          stop 'UNLOADING: exceed the maximum iteration (iter)!'
        endif

        q_epsr = -1.d0 - ar*e*deps/cohn(2)

        resq = resq/q_epsr

        peps = peps - resq

        goto 100                 ! ---> go back to 100 

      endif

      return
      end
