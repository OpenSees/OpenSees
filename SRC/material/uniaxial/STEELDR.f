	subroutine STEELDR(e_s,f_s,Etangm,d,hr)

	implicit none    
	integer isw,igxp1,il1,ik1,iet1,nvgp1,ix(2,1)
	real*8 e_s,e_so,f_so,fps_so,rejoin,val1,val2,val3,ep_o1,ttim
	real*8 region,Em_s,f_y,e_su,f_su,e_sh,e_sh1,f_sh1,P_major
	real*8 P_minor,f_s,fps_s,d(14),ep_so,fp_so,sum1,sum2,yield1
	real*8 ep_a,ep_M,epp_N,ep_o,ep_p,ep_r,ep_s,ep_sh,ep_sh1,ep_su
	real*8 ep_sushift,Ep_u,fp_a,fp_p,fp_r,fp_s,fp_sh,fp_sh1,repeat
	real*8 fp_su,fps_su,fp_t,P,W,s1,fps_target,point(6,3),hr(30)
	real*8 fps_rejoin,fp_rejoin,ep_rejoin,hist1(2),sim1
	real*8 D_dam,Dam_new,Dam_old,de_pl,fp_smid,Etangm
	dimension ep_o(2),ep_sushift(2)
	integer load,is1,k,m,icheat,i1,m1
c
c
c	Function to define stress-strain behavior of steel under cyclic loading,
c     using the DODD-RESTREPO material model
c
c	Last updated 10/6/13 - Implementation by I. Koutromanos and J. Bowers, Virginia Tech

c     INPUT:
c     d     =   material parameter array
c     e_s   =   uniaxial strain
c     ix    =   element nodes (ONLY APPLICABLE IN FEAP!)
c     
c     NOTE: ix is not really necessary, it is only passed for debugging purposes
c
c     
c     OUTPUT:
c     Etangm =  tangent modulus
c     f_s   =   steel stress

c	List of variables used in the model

c     Dam_old, Dam_new =  old and current value of damage variable (used when
c           low-cycle fatigue is included in the model).
c	ep_a = natural strain at initiation of Bauschinger effect
c	ep_M = maximum plastic strain defined as greater of |ep_o(1)| or ep_o(2)
c	epp_N = transformed-normalized natural strain ("epsilon double-prime)
c	ep_o(k) = "shifted" origin strain (k = 1 for compression, 2 for tension)
c	ep_p = strain magnitude between inclined envelope and point at initiaion
c	   of Bauschinger effect in natural coordinate system
c	ep_r(m) = natural strain at reversal point (m = 1 for max tension, m = 2
c	   for max compression)
c	e_s = steel strain in engineering coordinate system
c	ep_s = steel strain in natural coordinate system
c	ep_sh = natural strain at initiation of work hardening
c	ep_sh1 = arbitrary natural strain in work hardening region of tension
c	   monotonic curve
c	ep_su = "ultimate" strain in natural coordinate system
c	ep_sushift(k) = natural strain associated with "ultimate" true stress
c	Ep_u = natural coordinate unload modulus
c	Em_s = initial elastic modulus
c	fp_a = true stress at initiation of Bauschinger effect
c	fp_p = stress magnitude between inclined envelope and point at initiation 
c	   of Bauschinger effect in natural coordinate system
c	fp_r(m) = true stress at reversal point (m = 1 for max tension, m = 2
c	   for max compression)
c	f_s = engineering stress of steel
c	fp_s = true stress of steel
c	fp_sh = true stress at onset of strain hardening
c	fp_sh1 = true stress at arbitrary point ep_sh1
c	fp_su = "ultimate stress" in natural coordinates
c	fps_su = slope at "ultimate stress" in natural coordinates, numerically
c	   equal to fo_su
c	fp_t = stress magnitude between top and bottom inclined envelopes in 
c	   natural coordinate system
c	f_y = yield stress of steel in engineering coordinate system
c	load = loading direction (1 for compression, 2 for tension)
c	P = exponential factor to define work hardening region of skeleton curve 
c	   and shape of Bauschinger effect
c	point: array with dimension (5x3), it includes information about 
c               points in the (ep_s, fp_s) required to capture history effects 
c				point(i,1)= ep_s of stored point "i"
c               point(i,2)= fp_s of stored point "i"
c               point(i,3)= slope of curve fp_s(ep_s) at stored point "i"
c	region = at which region of the stress-strain curve we are located:

c	         region=0: We load along the skeleton curve
c	         region=1: Major reversal, (s1=-1)
c	         region=2: Major reversal, (s1=1)
c	         region=3: Minor reversal, "minor point" on positive side (s1=-1)
c	         region=4: Minor reversal, "minor point" on negative side (s1=1)
c	         region=5: Simple reversal, s1=-1,inside a "region3" 
c	         region=6: Simple reversal, s1=-1,inside a "region4"
c	         region=7: Simple reversal, s1=1,inside a "region3" 
c	         region=8: Simple reversal, s1=1,inside a "region4" 

c	s1 = strain direction factor (1 for tension, -1 for compression)
c	virgin = 1 for virgin loading, 0 for non-virgin
c	W = (omega) partial area under normalized Bauschinger curve

c
c
c     HISTORY VARIABLES STORED 
c     hr(1)     engineering strain of previous step
c     hr(2)     engineering stress of previous step
c     hr(3)     yield1: parameter, =0 if we have not yielded, =1 if we have yielded but not reached
c                           the strain hardening region yet, =2 if we have also entered the strain hardening region
c     hr(4)     region of the previous step
c     hr(5)     point(1,1)
c     hr(6)     point(1,2)
c     hr(7)     point(1,3)
c     hr(8)     point(2,1)
c     hr(9)     point(2,2)
c     hr(10)     point(2,3)
c     hr(11)     point(3,1)
c     hr(12)     point(3,2)
c     hr(13)     point(3,3)
c     hr(14)     point(4,1)
c     hr(15)     point(4,2)
c     hr(16)     point(4,3)
c     hr(17)     point(5,1)
c     hr(18)     point(5,2)
c     hr(19)     point(5,3)
c     hr(20)     ep_o(1)
c     hr(21)     ep_o(2)
c     hr(22)     ep_M ("maximum strain" which affects "unloading elastic modulus")
c     hr(23)     slope of stress strain diagram in f-prime, epsilon-prime space for previous step
c     hr(24)     hist1(1)   (it is initially 0, it becomes 1 when we have major unloading with s=1 and "loose" the original envelope in positive direction)
c     hr(25)     hist1(2)   (it is initially 0, it becomes 1 when we have major unloading with s=-1 and "loose" the original envelope in negative direction)
c     hr(26)     point(6,1)
c     hr(27)     point(6,2)
c     hr(28)     point(6,3)
c     hr(29)     sim1
c
c     NOTE: Parameters (26-29) are auxiliary, they are used when we have simple reversal, then reverse the direction of loading for one step,
c           then in the next step we restore the original loading direction
c
c     hr(30)     Damage parameter used for low-cycle fatigue


c
c

c	---------------------------------

c     Check whether low-cycle fatigue induced rupture has occurred:
c
      if(d(14).gt.1.d-8.and.
     *hr(30).lt.-99.d0) then
      f_s=0.d0
      return
      endif
c     --------------------------------------------------


	e_so=hr(1)
	f_so=hr(2)
	yield1=hr(3)
	region=hr(4)
	point(1,1)=hr(5)
	point(1,2)=hr(6)
	point(1,3)=hr(7)
	point(2,1)=hr(8)
	point(2,2)=hr(9)
	point(2,3)=hr(10)
	point(3,1)=hr(11)
	point(3,2)=hr(12)
	point(3,3)=hr(13)
	point(4,1)=hr(14)
	point(4,2)=hr(15)
	point(4,3)=hr(16)
	point(5,1)=hr(17)
	point(5,2)=hr(18)
	point(5,3)=hr(19)
	ep_o(1)=hr(20)
	ep_o(2)=hr(21)
	ep_M=hr(22)
	fps_so=hr(23)
	hist1(1)=hr(24)
	hist1(2)=hr(25)
	point(6,1)=hr(26)
	point(6,2)=hr(27)
	point(6,3)=hr(28)
	sim1=hr(29)
      Dam_old=hr(30)
      Dam_new=Dam_old

	Em_s=d(1)
	f_y=d(3)
	e_sh=d(4)
	e_sh1=d(5)
	f_sh1=d(6)
	e_su=d(7)
	f_su=d(8)
	P_major=d(9)
	P_minor=d(10)


	ep_sh=dlog(1.d0+e_sh) 
	fp_sh=f_y*dexp(ep_sh) 
	ep_sh1=dlog(1.d0+e_sh1) 
	fp_sh1=f_sh1*dexp(ep_sh1)
	ep_su=dlog(1.d0+e_su) 
	fp_su=f_su*dexp(ep_su) 
	fps_su=dexp(ep_su)*f_su 

	ep_s = dlog(1.d0+e_s) 
	ep_so = dlog(1.d0+e_so)                                     
c	converting from engineering strain to true strain

	fp_so = f_so*(1.d0+e_so)

c	Calculate Ep_u (unloading-Reloading modulus)

	IF (YIELD1.GT.0.5d0) then

	
	Ep_u=Em_s*(0.82d0+1.d0/(5.55d0+1000.d0*ep_M))

	ELSE

	Ep_u=Em_s

	ENDIF

	ep_sushift(1)=ep_su+ep_o(1)
	ep_sushift(2)=-ep_su+ep_o(2)

c	if you want to cheat and set P = P_major, 0 if you want to use Eq 36 to solve for P	
	if (P_major.gt.0) then
	   icheat = 1
	else
	   icheat = 0
	endif
	
	

	IF((EP_S.GE.POINT(1,1)).OR.(EP_S.LE.POINT(2,1))
     *	.or.ep_s.ge.ep_sushift(1).or.ep_s.le.ep_sushift(2)) THEN

 1101   format(4(1x,e15.5))
		if(ep_s.ge.point(1,1).or.(ep_s).ge.ep_sushift(1)) then

			ep_o1=ep_o(1)
				if(hist1(1).gt.0.5d0) then
			    s1=1.d0
				k=1
				ep_a=point(2,1)+s1*f_y/Ep_u
				fp_a=point(2,2)+s1*f_y
				ep_sushift(k)=s1*ep_su+ep_o(k)
				ep_rejoin=ep_sushift(k)
				fp_rejoin=s1*fp_su
				fps_rejoin=fps_su
					if(ep_s.ge.ep_sushift(k)) then
 
					fp_s=s1*fp_su
					fps_s=fps_su
					ep_sushift(k)=ep_s
					ep_o(k)=ep_sushift(k)-s1*ep_su
c	
					else
			        call Bauschinger(1,icheat,P_major,P_minor,region,
     *   ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *   ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	 ep_o,ix,ep_so,fp_so)
					endif
                  else
	            call virginLoading(ep_s,ep_o1,Em_s,f_y,ep_sh,ep_su,
     *   fp_sh1,ep_o,fps_su,fp_su,ep_sh1,fp_sh,fp_s,yield1,fps_s,ix)
                  
                  
                          if(ep_s.ge.ep_sushift(1)) then
                  	    fp_s=fp_su
					    fps_s=fps_su
					    ep_sushift(1)=ep_s
					    ep_o(1)=ep_sushift(1)-ep_su
                          endif
                  
                  endif
				 

		else
c		
			ep_o1=ep_o(2)
				if(hist1(2).gt.0.5d0) then
			    s1=-1.d0
				k=2
				ep_a=point(1,1)+s1*f_y/Ep_u
				fp_a=point(1,2)+s1*f_y
				ep_sushift(k)=s1*ep_su+ep_o(k)
				ep_rejoin=ep_sushift(k)
				fp_rejoin=s1*fp_su
				fps_rejoin=fps_su

					if(ep_s.le.ep_sushift(k)) then
					fp_s=s1*fp_su
					fps_s=fps_su
					else
	    call Bauschinger(1,icheat,P_major,P_minor,region,
     *  ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)
					endif
				else

	            call virginLoading(ep_s,ep_o1,Em_s,f_y,ep_sh,ep_su,
     *	fp_sh1, ep_o,fps_su,fp_su,ep_sh1,fp_sh,fp_s,yield1,fps_s,ix)
                  
                          if(ep_s.le.ep_sushift(2)) then
                  	    fp_s=-fp_su
					    fps_s=fps_su
					    ep_sushift(2)=ep_s
					    ep_o(2)=ep_sushift(2)+ep_su
                          endif
                  
                  
				endif

			
		endif
c	
	
	
	
	

	region=0.d0
	point(3,1)=0.d0
	point(3,2)=0.d0
	point(3,3)=0.d0
	point(4,1)=0.d0
	point(4,2)=0.d0
	point(4,3)=0.d0
	point(5,1)=0.d0
	point(5,2)=0.d0
	point(5,3)=0.d0
		If(ep_s.gt.point(1,1).or.ep_s.ge.ep_sushift(1)) then
		point(1,1)=ep_s
		point(1,2)=fp_s
		point(1,3)=fps_s
		endif


		If(ep_s.lt.point(2,1).or.ep_s.le.ep_sushift(2)) then
		point(2,1)=ep_s
		point(2,2)=fp_s
		point(2,3)=fps_s
		endif
	
	ELSEIF(EP_S.LT.POINT(1,1).AND.EP_S.GT.POINT(2,1)) THEN


c-------------------------------------------------------------------------------

		IF(REGION.LT.0.5d0) THEN

c			
			if(yield1.lt.0.5d0) then
			fp_s=Em_s*ep_s
			fps_s=Em_s


				If(ep_s.gt.point(1,1)) then
				point(1,1)=ep_s
				point(1,2)=fp_s
				point(1,3)=fps_s
				endif


				If(ep_s.lt.point(2,1)) then
				point(2,1)=ep_s
				point(2,2)=fp_s
				point(2,3)=fps_s
	endif
		
			elseif(yield1.lt.1.5d0) then
c	If We have reversal but are still in the yield plateau
				if(ep_so.lt.point(1,1).and.ep_so.gt.point(2,1)) then
c	
					if(point(5,2).gt.0.d0) then
					s1=-1.d0
					k=2
					m=1
					region=1.d0
					else
					s1=1.d0
					k=1
					m=2
					region=2.d0
					endif
				else
				point(5,1)=ep_so
				point(5,2)=fp_so
				point(5,3)=fps_so
					if(fp_so.gt.0.d0) then
					s1=-1.d0
					k=2
					m=1
					region=1.d0
					else
					s1=1.d0
					k=1
					m=2
					region=2.d0
					endif
				endif


			ep_r=point(5,1)
			fp_r=point(5,2)
c			fp_rejoin=s1*f_y
			ep_rejoin=point(m,1)+s1*(ep_o(2)-ep_o(1)+2.d0*f_y/Ep_u)
              fp_rejoin=s1*f_y*dexp(ep_rejoin)
			fps_rejoin=f_y
			if(s1.lt.0.d0) then
c		If we have yielded in tension and now unloading
c		define the "rejoin point" as the minimum attained strain!
				point(2,1)=ep_rejoin
				point(2,2)=fp_rejoin
				point(2,3)=f_y
			else
c		If we have yielded in compression and now unloading
c		define the "rejoin point" as the maximum attained strain!
				point(1,1)=ep_rejoin
				point(1,2)=fp_rejoin
				point(1,3)=f_y
		endif
			ep_sushift(k)=s1*ep_su+ep_o(k)
			ep_a=ep_r+s1*f_y/Ep_u
			fp_a=fp_r+s1*f_y
		
				if(dabs(ep_s-ep_r).le.dabs(ep_a-ep_r)) then
c	
				fp_s=fp_r+Ep_u*(ep_s-ep_r)
			
				else
				call Bauschinger(0,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)

				endif

			If(ep_s.gt.point(1,1)) then
			point(1,1)=ep_s
			point(1,2)=fp_s
			point(1,3)=fps_s
			endif


			If(ep_s.lt.point(2,1)) then
			point(2,1)=ep_s
			point(2,2)=fp_s
			point(2,3)=fps_s
		endif
				
			elseif(yield1.gt.1.5d0) then
c	******		IF WE HAVE ENTERED THE HARDENING REGION				
			
				point(5,1)=ep_so
				point(5,2)=fp_so
				point(5,3)=fps_so
c			
c     
				if(ep_s.lt.ep_so) then
					s1=-1.d0
					k=2
					m=1
					region=1.d0
					hist1(2)=1.d0
					point(1,1)=ep_so
					point(1,2)=fp_so
					point(1,3)=fps_so
					else
					s1=1.d0
					k=1
					m=2
					region=2.d0
					hist1(1)=1.d0
					point(2,1)=ep_so
					point(2,2)=fp_so
					point(2,3)=fps_so
					endif

c	


			ep_r=ep_so
			fp_r=fp_so
			ep_o(k)=ep_r+s1*fp_r/Ep_u
			ep_sushift(k)=s1*ep_su+ep_o(k)
			ep_a=ep_r+s1*f_y/Ep_u
			fp_a=fp_r+s1*f_y
			ep_rejoin=ep_sushift(k)
			fp_rejoin=s1*fp_su
			fps_rejoin=fps_su

c		



				if(dabs(ep_s-ep_r).le.dabs(ep_a-ep_r)) then
c	
				fp_s=fp_r+Ep_u*(ep_s-ep_r)
				fps_s=Ep_u
c			
				else
			call Bauschinger(1,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)


				endif



				
		endif


		ELSE


c----------------------------------------------------------------------		
			IF(REGION.GT.0.5d0.AND.REGION.LT.1.5d0) THEN 
     			CALL REVERSE1(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

c----------------------------------------------------------------------
			ELSEIF(REGION.GT.1.5d0.AND.REGION.LT.2.5d0) THEN 
     			CALL REVERSE2(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)


c--------------------------------------------------------------------
			ELSEIF(REGION.GT.2.5d0.AND.REGION.LT.3.5d0) THEN
		
			if(ep_s.lt.point(3,1)) then
			write(*,*) 'ton hpiame!'
			write(*,*) ix(1,1),ix(2,1),ep_s
			stop 
			endif

			if(ep_s.gt.point(5,1).and.ep_s.le.point(1,1)) then

				point(5,1)=point(2,1)
				point(5,2)=point(2,2)
				point(5,3)=point(2,3)
				point(4,1)=0.d0
				point(4,2)=0.d0
				point(4,3)=0.d0
				region=2.d0
	    		CALL REVERSE2(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			else
				 
     				CALL REVERSE3(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			endif
		


c----------------------------------------------------------------
			ELSEIF(REGION.GT.3.5d0.AND.REGION.LT.4.5d0) THEN


			if(ep_s.gt.point(3,1)) then
			write(*,*) 'ton hpiame2!'
			write(*,*) point(3,1),ep_s
				write(*,*) ix(1,1),ix(2,1),ep_s
			stop
			endif

			if(ep_s.lt.point(5,1).and.ep_s.ge.point(2,1)) then

				point(5,1)=point(1,1)
				point(5,2)=point(1,2)
				point(5,3)=point(1,3)
				point(4,1)=0.d0
				point(4,2)=0.d0
				point(4,3)=0.d0
				region=1.d0
	    		CALL REVERSE1(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			else
     			CALL REVERSE4(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			endif



c-------------------------------------------------------------
			ELSEIF(REGION.GT.4.5d0.AND.REGION.LT.5.5d0) THEN





			if(ep_s.gt.point(4,1).and.ep_s.le.point(1,1)) then
			point(5,1)=point(2,1)
			point(5,2)=point(2,2)
			point(5,3)=point(2,3)
			point(4,1)=0.d0
			point(4,2)=0.d0
			point(4,3)=0.d0
			region=2.d0
			sim1=0.d0
	    		CALL REVERSE2(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			
		else
     			CALL REVERSE5(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *  FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *  ,sim1,ix)
			endif


c-------------------------------------------------------------
			ELSEIF(REGION.GT.5.5d0.AND.REGION.LT.6.5d0) THEN
			if(ep_s.lt.point(4,1).and.ep_s.ge.point(2,1)) then

				
				point(5,1)=point(1,1)
				point(5,2)=point(1,2)
				point(5,3)=point(1,3)
				point(4,1)=0.d0
				point(4,2)=0.d0
				point(4,3)=0.d0
				region=1.d0
				sim1=0.d0
	    		CALL REVERSE1(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

			elseif((ep_s.gt.point(5,1).and.ep_s.lt.point(1,1)).and.
     *	sim1.lt.0.5d0) then

c				write(*,*) 'point(5,1)', point(5,1),point(5,2)
				point(5,1)=point(4,1)
				point(5,2)=point(4,2)
				point(5,3)=point(4,3)
				region=4.d0
				sim1=0.d0

	    		CALL REVERSE4(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

			else
     			CALL REVERSE6(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			endif


c-----------------------------------------------------------------
			ELSEIF(REGION.GT.6.5d0.AND.REGION.LT.7.5d0) THEN
			if(ep_s.gt.point(4,1).and.ep_s.le.point(1,1)) then


				point(5,1)=point(2,1)
				point(5,2)=point(2,2)
				point(5,3)=point(2,3)
				point(4,1)=0.d0
				point(4,2)=0.d0
				point(4,3)=0.d0
			region=2.d0
				sim1=0.d0
	    		CALL REVERSE2(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

			elseif(ep_s.lt.point(6,1).and.ep_s.gt.point(2,1).and.
     *	sim1.lt.0.5d0) then


				point(5,1)=point(4,1)
				point(5,2)=point(4,2)
				point(5,3)=point(4,3)
				region=3.d0
				sim1=0.d0
c
	    		CALL REVERSE3(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

			else
     			CALL REVERSE7(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			endif


c----------------------------------------------------------------------
		     ELSEIF(REGION.GT.7.5d0.AND.REGION.LT.8.5d0) THEN
			
			if(ep_s.lt.point(4,1).and.ep_s.ge.point(2,1))then


				point(5,1)=point(1,1)
				point(5,2)=point(1,2)
				point(5,3)=point(1,3)
				point(4,1)=0.d0
				point(4,2)=0.d0
				point(4,3)=0.d0
				region=1.d0
				sim1=0.d0
	    		CALL REVERSE1(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)


			else

     			CALL REVERSE8(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			endif
			ELSE
		  WRITE(*,*) 'COULD NOT FIND AN APPROPRIATE VALUE OF REGION!!'
				write(*,*) ix(1,1),ix(2,1),ep_s
			STOP
			
			ENDIF

		ENDIF
c---------------------------------------------------------------------------------





	ENDIF



c     --------------     
c     Accumulation of tensile stress work if d(14)>0
      if(d(14).gt.1.d-8) then

        fp_smid=0.5d0*(fp_s+fp_so)

        if(fp_smid.gt.0.d0) then
            de_pl=ep_s-ep_so-(fp_s-fp_so)/Ep_u
        D_dam=((fp_smid*fp_smid/(2.d0*d(1)*d(12)))**d(13))*dabs(de_pl)
            Dam_new=D_dam+Dam_old
        endif


        if(Dam_new.gt.d(14)) then
        fp_s=0.d0
        Dam_new=-100.d0
c     NOTE: The value Dam_new=-100 is simply used to "inform" the routine that the material has ruptured
        endif


      endif
c     ----------
 	




c------------------------------------

c	CALCULATE NEW STRESS AND UPDATE HISTORY VARIABLES


c
	f_s = fp_s/dexp(ep_s)
      Etangm=(fps_s-fp_s)*dexp(-2.d0*ep_s)


	if(yield1.gt.0.5d0.and.ep_s.lt.0.d0.and.ep_s.le.point(2,1)
     *    .and.((ep_s+0.002d0).lt.(-1.d0*ep_M))) 
     *	ep_M=dabs(ep_s+0.002d0)
	if(yield1.gt.0.5d0.and.ep_s.gt.0.d0.and.ep_s.ge.point(1,1)
     *	.and.((ep_s-0.002d0).gt.ep_M)) 
     *	ep_M=(ep_s-0.002d0)

	if(yield1.gt.0.5) then
	Ep_u=Em_s*(0.82+1./(5.55+1000.*ep_M))
	else
	Ep_u=Em_s
	endif

	if(yield1.gt.0.5d0.and.yield1.lt.1.5d0.and.region.lt.0.5d0) then
	ep_o(2)=point(1,1)-point(1,2)/Ep_u
	ep_o(1)=point(2,1)-point(2,2)/Ep_u
	endif

	


	hr(1)=e_s
	hr(2)=f_s
	hr(3)=yield1
	hr(4)=region
	hr(5)=point(1,1)
	hr(6)=point(1,2)
	hr(7)=point(1,3)
	hr(8)=point(2,1)
	hr(9)=point(2,2)
	hr(10)=point(2,3)
	hr(11)=point(3,1)
	hr(12)=point(3,2)
	hr(13)=point(3,3)
	hr(14)=point(4,1)
	hr(15)=point(4,2)
	hr(16)=point(4,3)
	hr(17)=point(5,1)
	hr(18)=point(5,2)
	hr(19)=point(5,3)
	hr(20)=ep_o(1)
	hr(21)=ep_o(2)
	hr(22)=ep_M
	hr(23)=fps_s
	hr(24)=hist1(1)
	hr(25)=hist1(2)
	hr(26)=point(6,1)
	hr(27)=point(6,2)
	hr(28)=point(6,3)
	hr(29)=sim1
      hr(30)=Dam_new
	
	




 1001 format(4(1x,e20.5))

	return
	end subroutine STEELDR






c	=============================================================================================reverse1


	subroutine reverse1(region,point,ep_s,ep_so,fp_so,Ep_u,fp_s,
     *	fps_s,ep_a,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

	    


	implicit none
	real*8 region,point(6,3),ep_s,ep_so,fp_so,Ep_u,fp_s,fps_s,s1
	real*8 ep_a,yield1,P_major,P_minor,f_y,d(14),ep_rejoin,fp_rejoin
	real*8 fps_rejoin,fp_a,ep_o(2),fps_so,hist1(2),sim1
	real*8 ep_su,fp_su,fps_su,e_su,f_su,ep_sushift(2)
	integer k,m,icheat,ix(2,1)



	f_y=d(3)

	e_su=d(7)
	f_su=d(8)
	P_major=d(9)
	P_minor=d(10)



	ep_su=dlog(1.d0+e_su) 
	fp_su=f_su*dexp(ep_su) 
	fps_su=dexp(ep_su)*f_su 


	s1=-1.d0
	k=2
	m=1


c**********************************************************************************

	IF(YIELD1.LT.1.5d0) THEN
	


	
			If(yield1.lt.0.5d0) then
			write(*,*) 'YIELD1 VALUE OUT OF BOUNDS IN REVERSE1!!!'
				write(*,*) ix(1,1),ix(2,1),ep_s
			stop
			endif



		If((ep_s-ep_so).le.0.0d0) then
c		we continue unloading


		ep_sushift(k)=s1*ep_su+ep_o(k)
		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y
		

			if((point(5,1)-ep_s).le.(point(5,1)-ep_a)) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u
			else
c			fp_rejoin=s1*f_y
			ep_rejoin=point(m,1)+s1*(ep_o(2)-ep_o(1)+2.d0*f_y/Ep_u)
              fp_rejoin=s1*f_y*dexp(ep_rejoin)
			fps_rejoin=f_y

			call Bauschinger(0,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)


			endif


		else


		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y

			if(((point(5,1)-ep_so).le.(point(5,1)-ep_a)).and.
     *	((point(5,1)-ep_s).le.(point(5,1)-ep_a))) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u

			else

c		We have a MINOR reversal (unloading inside a major reloading in the YIELD PLATEAU
c		is ALWAYS minor!)
		
		point(3,1)=point(5,1)
		point(3,2)=point(5,2)
		point(3,3)=point(5,3)

		point(4,1)=ep_so
		point(4,2)=fp_so
		point(4,3)=fps_so

		point(5,1)=ep_so
		point(5,2)=fp_so
		point(5,3)=fps_so

		region=4.d0


     			CALL REVERSE4(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	 FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			endif

		endif


c**************************************************************************************


	ELSE

	ep_a=point(5,1)+s1*f_y/Ep_u
	fp_a=point(5,2)+s1*f_y


		IF(ep_s.le.ep_so) then






			ep_sushift(k)=s1*ep_su+ep_o(k)
			ep_a=point(5,1)+s1*f_y/Ep_u
			fp_a=point(5,2)+s1*f_y
			ep_rejoin=ep_sushift(k)
			fp_rejoin=s1*fp_su
			fps_rejoin=fps_su

c	

c	
c------------------------------------------------------------------
				IF((point(5,1)-ep_s).le.(point(5,1)-ep_a)) THEN
				fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u
				ELSE
c	
			call Bauschinger(1,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)

				ENDIF
c--------------------------------------------------------------------

		ELSE


			if(((point(5,1)-ep_so).le.(point(5,1)-ep_a)).and.
     *	((point(5,1)-ep_s).le.(point(5,1)-ep_a))) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u

			else

				if((point(5,2)-fp_so).ge.2.d0*f_y) then

c				We have a MAJOR reversal 


				point(5,1)=ep_so
				point(5,2)=fp_so
				point(5,3)=fps_so
				
				point(2,1)=ep_so
				point(2,2)=fp_so
				point(2,3)=fps_so


				region=2.d0


			    s1=1.d0
				k=1
				ep_o(k)=ep_so-s1*fp_so/Ep_u
				ep_sushift(k)=s1*ep_su+ep_o(k)


				ep_o(1)=ep_so+fp_so/Ep_u
				ep_sushift(1)=ep_o(1)+ep_su
				if(hist1(1).lt.0.5d0) hist1(1)=1.d0
c				
c			
c				
				

     				CALL REVERSE2(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	 FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	 ,sim1,ix)

				else

c				We have a MINOR reversal 


				region=4.d0
c
     
			point(3,1)=point(5,1)
		point(3,2)=point(5,2)
		point(3,3)=point(5,3)

		point(4,1)=ep_so
		point(4,2)=fp_so
		point(4,3)=fps_so

		point(5,1)=ep_so
		point(5,2)=fp_so
		point(5,3)=fps_so
     			
     				CALL REVERSE4(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	   FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	   ,sim1,ix)

				endif
			endif

		ENDIF


	ENDIF
c********************************************************************************


	return
	end subroutine reverse1

c	=============================================================================================	











c	===============================================================================================reverse2


	subroutine reverse2(region,point,ep_s,ep_so,fp_so,Ep_u,fp_s,
     *	fps_s,ep_a,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

	    

	implicit none
	real*8 region,point(6,3),ep_s,ep_so,fp_so,Ep_u,fp_s,fps_s,s1
	real*8 ep_a,yield1,P_major,P_minor,f_y,d(14),ep_rejoin,fp_rejoin
	real*8 fps_rejoin,fp_a,ep_o(2),fps_so,hist1(2),sim1
	real*8 ep_su,fp_su,fps_su,e_su,f_su,ep_sushift(2)
	integer k,m,icheat,ix(2,1)



	f_y=d(3)

	e_su=d(7)
	f_su=d(8)
	P_major=d(9)
	P_minor=d(10)



	ep_su=dlog(1.d0+e_su) 
	fp_su=f_su*dexp(ep_su) 
	fps_su=dexp(ep_su)*f_su 


	s1=1.d0
	k=1
	m=2


c**********************************************************************************

	IF(YIELD1.LT.1.5d0) THEN
	
c	write(*,*) 'yield1<1.5'

	
			If(yield1.lt.0.5d0) then
			write(*,*) 'YIELD1 VALUE OUT OF BOUNDS IN REVERSE2!!!'
				write(*,*) ix(1,1),ix(2,1),ep_s
			stop
			endif



		If(ep_s.ge.ep_so) then
c		we continue unloading


		ep_sushift(k)=s1*ep_su+ep_o(k)
		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y
		

			if(((point(5,1)-ep_s).ge.(point(5,1)-ep_a))) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u
c	
			else
c			fp_rejoin=s1*f_y
			ep_rejoin=point(m,1)+s1*(ep_o(2)-ep_o(1)+2.d0*f_y/Ep_u)
              fp_rejoin=s1*f_y*dexp(ep_rejoin)
			fps_rejoin=f_y
c		

			call Bauschinger(0,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)


			endif


		else


		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y

			if(((point(5,1)-ep_so).ge.(point(5,1)-ep_a)).and.
     *	((point(5,1)-ep_s).ge.(point(5,1)-ep_a))) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u

			else

c		We have a MINOR reversal (reversal inside a major reloading in the YIELD PLATEAU
c		is ALWAYS minor!)
		
		point(3,1)=point(5,1)
		point(3,2)=point(5,2)
		point(3,3)=point(5,3)

		point(4,1)=ep_so
		point(4,2)=fp_so
		point(4,3)=fps_so

		point(5,1)=ep_so
		point(5,2)=fp_so
		point(5,3)=fps_so

		region=3.d0
c			
     			CALL REVERSE3(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)
			endif

		endif


c**************************************************************************************


	ELSE
	ep_a=point(5,1)+s1*f_y/Ep_u
	fp_a=point(5,2)+s1*f_y
c	
		IF(ep_s.ge.ep_so) then






			ep_sushift(k)=s1*ep_su+ep_o(k)
			ep_a=point(5,1)+s1*f_y/Ep_u
			fp_a=point(5,2)+s1*f_y
			ep_rejoin=ep_sushift(k)
			fp_rejoin=s1*fp_su
			fps_rejoin=fps_su

c		

c------------------------------------------------------------------
				IF(((point(5,1)-ep_s).ge.(point(5,1)-ep_a))) THEN
				fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u
				ELSE

			call Bauschinger(1,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)
				ENDIF
c--------------------------------------------------------------------

		ELSE


			if(((point(5,1)-ep_so).ge.(point(5,1)-ep_a)).and.
     *	((point(5,1)-ep_s).ge.(point(5,1)-ep_a))) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u

			else

				if((fp_so-point(5,2)).ge.2.d0*f_y) then

c				We have a MAJOR reversal! 


				point(5,1)=ep_so
				point(5,2)=fp_so
				point(5,3)=fps_so

				point(1,1)=ep_so
				point(1,2)=fp_so
				point(1,3)=fps_so


				if(hist1(2).lt.0.5d0) hist1(2)=1.d0

				region=1.d0

			    s1=-1.d0
				k=2
				ep_o(k)=ep_so+s1*fp_so/Ep_u
				ep_sushift(k)=s1*ep_su+ep_o(k)

				

     				CALL REVERSE1(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

				else

c				We have a MINOR reversal 

		point(3,1)=point(5,1)
		point(3,2)=point(5,2)
		point(3,3)=point(5,3)

		point(4,1)=ep_so
		point(4,2)=fp_so
		point(4,3)=fps_so

		point(5,1)=ep_so
		point(5,2)=fp_so
		point(5,3)=fps_so

				region=3.d0
c	
				
     				CALL REVERSE3(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

				endif
			endif

		ENDIF


	ENDIF
c********************************************************************************

	return
	end subroutine reverse2

c	=============================================================================================	






c	===============================================================================================reverse3


	subroutine reverse3(region,point,ep_s,ep_so,fp_so,Ep_u,fp_s,
     *	fps_s,ep_a,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

	    


	implicit none
	real*8 region,point(6,3),ep_s,ep_so,fp_so,Ep_u,fp_s,fps_s,s1
	real*8 ep_a,yield1,P_major,P_minor,f_y,d(14),ep_rejoin,fp_rejoin
	real*8 fps_rejoin,fp_a,ep_o(2),fps_so,hist1(2),sim1
	real*8 ep_su,fp_su,fps_su,e_su,f_su,ep_sushift(2)
	integer k,m,icheat,ix(2)


	f_y=d(3)

	e_su=d(7)
	f_su=d(8)
	P_major=d(9)
	P_minor=d(10)



	ep_su=dlog(1.d0+e_su) 
	fp_su=f_su*dexp(ep_su) 
	fps_su=dexp(ep_su)*f_su 


	s1=-1.d0
	k=2
	m=1


c**************************************************************************************

	IF(ep_s.le.ep_so) then

		ep_sushift(k)=s1*ep_su+ep_o(k)
		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y
		ep_rejoin=point(3,1)
		fp_rejoin=point(3,2)
		fps_rejoin=point(3,3)
c------------------------------------------------------------------
			IF(((point(5,1)-ep_s).le.(point(5,1)-ep_a))) THEN
c	
				fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u
			ELSE

				call Bauschinger(0,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)

			ENDIF
c--------------------------------------------------------------------

	ELSE


		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y

		if(((point(5,1)-ep_so).le.(point(5,1)-ep_a)).and.
     *	((point(5,1)-ep_s).le.(point(5,1)-ep_a))) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u

		else

c				We have a SIMPLE reversal 

			point(4,1)=point(5,1)
			point(4,2)=point(5,2)
			point(4,3)=point(5,3)

			point(6,1)=ep_so
			point(6,2)=fp_so
			point(6,3)=fps_so

			region=7.d0

c	
		
     			CALL REVERSE7(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

		endif

	ENDIF


c********************************************************************************
	return
	end subroutine reverse3

c	=============================================================================================	




c	===============================================================================================reverse4


	subroutine reverse4(region,point,ep_s,ep_so,fp_so,Ep_u,fp_s,
     *	fps_s,ep_a,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

	    


	implicit none
	real*8 region,point(6,3),ep_s,ep_so,fp_so,Ep_u,fp_s,fps_s,s1
	real*8 ep_a,yield1,P_major,P_minor,f_y,d(14),ep_rejoin,fp_rejoin
	real*8 fps_rejoin,fp_a,ep_o(2),fps_so,hist1(2),sim1
	real*8 ep_su,fp_su,fps_su,e_su,f_su,ep_sushift(2)
	integer k,m,icheat,ix(2,1)

	

	f_y=d(3)

	e_su=d(7)
	f_su=d(8)
	P_major=d(9)
	P_minor=d(10)

c

	ep_su=dlog(1.d0+e_su) 
	fp_su=f_su*dexp(ep_su) 
	fps_su=dexp(ep_su)*f_su 


	s1=1.d0
	k=1
	m=2


c**************************************************************************************

	IF(ep_s.ge.ep_so) then

		ep_sushift(k)=s1*ep_su+ep_o(k)
		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y
		ep_rejoin=point(3,1)
		fp_rejoin=point(3,2)
		fps_rejoin=point(3,3)
c	

c------------------------------------------------------------------
			IF(((point(5,1)-ep_s).ge.(point(5,1)-ep_a))) THEN
				fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u
			ELSE

				call Bauschinger(0,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)

			ENDIF
c--------------------------------------------------------------------

	ELSE

		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y

		if(((point(5,1)-ep_so).ge.(point(5,1)-ep_a)).and.
     *	((point(5,1)-ep_s).ge.(point(5,1)-ep_a))) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u

		else

c				We have a SIMPLE reversal 

			point(4,1)=point(5,1)
			point(4,2)=point(5,2)
			point(4,3)=point(5,3)

			point(5,1)=ep_so
			point(5,2)=fp_so
			point(5,3)=fps_so

			region=6.d0
		
     			CALL REVERSE6(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

		endif

	ENDIF


c********************************************************************************

	return
	end subroutine reverse4

c	=============================================================================================	








c	===============================================================================================reverse5


	subroutine reverse5(region,point,ep_s,ep_so,fp_so,Ep_u,fp_s,
     *	fps_s,ep_a,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

	    


	implicit none
	real*8 region,point(6,3),ep_s,ep_so,fp_so,Ep_u,fp_s,fps_s,s1
	real*8 ep_a,yield1,P_major,P_minor,f_y,d(14),ep_rejoin,fp_rejoin
	real*8 fps_rejoin,fp_a,ep_o(2),fps_so,hist1(2),sim1
	real*8 ep_su,fp_su,fps_su,e_su,f_su,ep_sushift(2)
	integer k,m,icheat,ix(2,1)



	f_y=d(3)

	e_su=d(7)
	f_su=d(8)
	P_major=d(9)
	P_minor=d(10)



	ep_su=dlog(1.d0+e_su) 
	fp_su=f_su*dexp(ep_su) 
	fps_su=dexp(ep_su)*f_su 


	s1=-1.d0
	k=2
	m=1


c**************************************************************************************

	IF(ep_s.le.ep_so) then

		ep_sushift(k)=s1*ep_su+ep_o(k)
		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y
		ep_rejoin=point(3,1)
		fp_rejoin=point(3,2)
		fps_rejoin=point(3,3)
c------------------------------------------------------------------
			IF((point(5,1)-ep_s).le.(point(5,1)-ep_a)) THEN
				fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u
			ELSE
				
				if(sim1.gt.0.5d0) then
				sim1=0.d0
				point(6,1)=0.d0
				point(6,2)=0.d0
				point(6,3)=0.d0
				endif

				call Bauschinger(0,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)

			ENDIF
c--------------------------------------------------------------------

	ELSE

		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y

		if(((point(5,1)-ep_so).le.(point(5,1)-ep_a)).and.
     *	((point(5,1)-ep_s).le.(point(5,1)-ep_a)).and.
     *	(point(5,1).ge.ep_s)) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u

		else

c				We have a new SIMPLE reversal 

			if(sim1.lt.0.5d0) then
			point(6,1)=ep_so
			point(6,2)=fp_so
			point(6,3)=fps_so
			sim1=1.d0
			else
			point(5,1)=0.d0
			point(5,2)=0.d0
			point(5,3)=0.d0
			sim1=0.d0
			endif
			
			region=7.d0
		
     			CALL REVERSE7(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

		endif

	ENDIF


c********************************************************************************
	return
	end subroutine reverse5

c	=============================================================================================	




c	===============================================================================================reverse6


	subroutine reverse6(region,point,ep_s,ep_so,fp_so,Ep_u,fp_s,
     *	fps_s,ep_a,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

	    


	implicit none
	real*8 region,point(6,3),ep_s,ep_so,fp_so,Ep_u,fp_s,fps_s,s1
	real*8 ep_a,yield1,P_major,P_minor,f_y,d(14),ep_rejoin,fp_rejoin
	real*8 fps_rejoin,fp_a,ep_o(2),fps_so,hist1(2),sim1
	real*8 ep_su,fp_su,fps_su,e_su,f_su,ep_sushift(2)
	integer k,m,icheat,ix(2,1)

	

	f_y=d(3)

	e_su=d(7)
	f_su=d(8)
	P_major=d(9)
	P_minor=d(10)



	ep_su=dlog(1.d0+e_su) 
	fp_su=f_su*dexp(ep_su) 
	fps_su=dexp(ep_su)*f_su 


	s1=-1.d0
	k=2
	m=1


c**************************************************************************************

	IF(ep_s.le.ep_so) then

		ep_sushift(k)=s1*ep_su+ep_o(k)
		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y
		ep_rejoin=point(4,1)
		fp_rejoin=point(4,2)
		fps_rejoin=point(4,3)
c------------------------------------------------------------------
			IF((point(5,1)-ep_s).le.(point(5,1)-ep_a)) THEN
c			
				fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u
			ELSE

				if(sim1.gt.0.5d0) then
				sim1=0.
				point(6,1)=0.d0
				point(6,2)=0.d0
				point(6,3)=0.d0
				endif

				call Bauschinger(0,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)

			ENDIF
c--------------------------------------------------------------------

	ELSE

		ep_a=point(5,1)+s1*f_y/Ep_u
		fp_a=point(5,2)+s1*f_y

c	
		if(((point(5,1)-ep_so).le.(point(5,1)-ep_a)).and.
     *	 ((point(5,1)-ep_s).le.(point(5,1)-ep_a)).and.
     *	 (point(5,1).ge.ep_s)) then
			fp_s=point(5,2)+Ep_u*(ep_s-point(5,1))
			fps_s=Ep_u

		else

c				We have a new SIMPLE reversal in the same region but with opposite direction 


			if(sim1.lt.0.5d0) then
			point(6,1)=ep_so
			point(6,2)=fp_so
			point(6,3)=fps_so
			sim1=1.d0
			else
			point(5,1)=0.d0
			point(5,2)=0.d0
			point(5,3)=0.d0
			sim1=0.
			endif

			region=8.d0

			
		
     			CALL REVERSE8(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

		endif

	ENDIF


c********************************************************************************

c	=============================================================================================	

	return
	end subroutine reverse6













c	===============================================================================================reverse7


	subroutine reverse7(region,point,ep_s,ep_so,fp_so,Ep_u,fp_s,
     *	fps_s,ep_a,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

	    


	implicit none
	real*8 region,point(6,3),ep_s,ep_so,fp_so,Ep_u,fp_s,fps_s,s1
	real*8 ep_a,yield1,P_major,P_minor,f_y,d(14),ep_rejoin,fp_rejoin
	real*8 fps_rejoin,fp_a,ep_o(2),fps_so,hist1(2),sim1
	real*8 ep_su,fp_su,fps_su,e_su,f_su,ep_sushift(2)
	integer k,m,icheat,ix(2,1)



	f_y=d(3)

	e_su=d(7)
	f_su=d(8)
	P_major=d(9)
	P_minor=d(10)



	ep_su=dlog(1.d0+e_su) 
	fp_su=f_su*dexp(ep_su) 
	fps_su=dexp(ep_su)*f_su 


	s1=1.d0
	k=1
	m=2


c**************************************************************************************

	IF(ep_s.ge.ep_so) then

		ep_sushift(k)=s1*ep_su+ep_o(k)
		ep_a=point(6,1)+s1*f_y/Ep_u
		fp_a=point(6,2)+s1*f_y
		ep_rejoin=point(4,1)
		fp_rejoin=point(4,2)
		fps_rejoin=point(4,3)
c------------------------------------------------------------------
			IF((point(6,1)-ep_s).ge.(point(6,1)-ep_a)) THEN
				fp_s=point(6,2)+Ep_u*(ep_s-point(6,1))
			fps_s=Ep_u
			ELSE

				if(sim1.gt.0.5d0) then
				sim1=0.d0
				point(5,1)=0.d0
				point(5,2)=0.d0
				point(5,3)=0.d0
				endif


				call Bauschinger(0,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)

			ENDIF
c--------------------------------------------------------------------

	ELSE


		ep_a=point(6,1)+s1*f_y/Ep_u
		fp_a=point(6,2)+s1*f_y

		if(((point(6,1)-ep_so).ge.(point(6,1)-ep_a)).and.
     *	((point(6,1)-ep_s).ge.(point(6,1)-ep_a)).and.
     *	(point(6,1).le.ep_s)) then
			fp_s=point(6,2)+Ep_u*(ep_s-point(6,1))
			fps_s=Ep_u

		else

c				We have a new SIMPLE reversal inside a region 3, but with opposite direction


			if(sim1.lt.0.5d0) then
			point(5,1)=ep_so
			point(5,2)=fp_so
			point(5,3)=fps_so
			sim1=1.d0
			else
			point(6,1)=0.d0
			point(6,2)=0.d0
			point(6,3)=0.d0
			sim1=0.
			endif

			region=5.d0
		
     			CALL REVERSE5(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

		endif

	ENDIF


c********************************************************************************
	return
	end subroutine reverse7

c	=============================================================================================	



c	===============================================================================================reverse8


	subroutine reverse8(region,point,ep_s,ep_so,fp_so,Ep_u,fp_s,
     *	fps_s,ep_a,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

	    


	implicit none
	real*8 region,point(6,3),ep_s,ep_so,fp_so,Ep_u,fp_s,fps_s,s1
	real*8 ep_a,yield1,P_major,P_minor,f_y,d(14),ep_rejoin,fp_rejoin
	real*8 fps_rejoin,fp_a,ep_o(2),fps_so,hist1(2),sim1
	real*8 ep_su,fp_su,fps_su,e_su,f_su,ep_sushift(2)
	integer k,m,icheat,ix(2,1)



	f_y=d(3)

	e_su=d(7)
	f_su=d(8)
	P_major=d(9)
	P_minor=d(10)



	ep_su=dlog(1.d0+e_su) 
	fp_su=f_su*dexp(ep_su) 
	fps_su=dexp(ep_su)*f_su 


	s1=1.d0
	k=1
	m=2


c**************************************************************************************

	IF(ep_s.ge.ep_so) then

		ep_sushift(k)=s1*ep_su+ep_o(k)
		ep_a=point(6,1)+s1*f_y/Ep_u
		fp_a=point(6,2)+s1*f_y
		ep_rejoin=point(3,1)
		fp_rejoin=point(3,2)
		fps_rejoin=point(3,3)
c------------------------------------------------------------------
			IF((point(6,1)-ep_s).ge.(point(6,1)-ep_a)) THEN
				fp_s=point(6,2)+Ep_u*(ep_s-point(6,1))
				fps_s=Ep_u
c							
			ELSE

				if(sim1.gt.0.5d0) then
				sim1=0.d0
				point(5,1)=0.d0
				point(5,2)=0.d0
				point(5,3)=0.d0
				endif

				call Bauschinger(0,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,
     *	ep_s,k,ep_rejoin,fp_rejoin,fps_rejoin,ep_su,fp_s,fps_s,
     *	ep_o,ix,ep_so,fp_so)

			ENDIF
c--------------------------------------------------------------------

	ELSE

		ep_a=point(6,1)+s1*f_y/Ep_u
		fp_a=point(6,2)+s1*f_y

		if(((point(6,1)-ep_so).ge.(point(6,1)-ep_a)).and.
     *	((point(6,1)-ep_s).ge.(point(6,1)-ep_a)).and.
     *	(point(6,1).le.ep_s)) then
			fp_s=point(6,2)+Ep_u*(ep_s-point(6,1))
			fps_s=Ep_u

		else

c				We have a new SIMPLE reversal inside a region 3, but with opposite direction


			if(sim1.lt.0.5d0) then
			point(5,1)=ep_so
			point(5,2)=fp_so
			point(5,3)=fps_so
			sim1=1.d0
			else
			point(6,1)=0.d0
			point(6,2)=0.d0
			point(6,3)=0.d0
			sim1=0.
			endif

			region=6.d0
			
     			CALL REVERSE6(REGION,POINT,EP_S,EP_SO,FP_SO,EP_U,FP_S,
     *	FPS_S,EP_A,yield1,d,ep_sushift,ep_o,fps_so,icheat,hist1
     *	,sim1,ix)

		endif

	ENDIF


c********************************************************************************
	return
	end subroutine reverse8

c	=============================================================================================	






c	=============================================================================================virginLoading


	subroutine virginLoading(ep_s,ep_o1,Em_s,f_y,ep_sh,ep_su,fp_sh1,
     *	ep_o,fps_su,fp_su,ep_sh1,fp_sh,val1,val2,val3,ix)    


	implicit none
	real*8 ep_s,ep_o(2),ep_o1,val1,val2,val3,f_y,Em_s,ep_sh,fp_sh1
	real*8 yield1
	real*8 ep_su,fp_su,P,fps_su,fp_sh,s1,d(14),ep_sh1,fp_s,aux1,aux2
	integer k,ix(2,1)

c	
	if (ep_s.ge.ep_o1)  then                                       
c	if in tensile loading
        s1 = 1.d0
	  k = 1
	else                                                           
c	if in compressive loading
        s1 = -1.d0
	  k =2
	endif


                                                                            
c	Model modified to account for possible shift due to reversal in yield plateau
      
	if ((dabs(ep_s-ep_o1).ge.0.d0).and.((dabs(ep_s-ep_o1)).lt.
     *	(f_y/Em_s)))  then                   
c	Elastic branch
        val1 = Em_s*ep_s                                              
c	  Eq 15
        val3 = Em_s
        val2=0.d0
      elseif (((dabs(ep_s-ep_o1)).ge.(f_y/Em_s))
     *	.and.((dabs(ep_s-ep_o1)).lt.ep_sh)) then              
c	Yield plateau
        val1 = s1*f_y*dexp(ep_s)
c	                                         
c	  Eq 16
        val3 = f_y
        if(val2.lt.0.5d0) val2=1.d0
      elseif ((dabs(ep_s-ep_o1).ge.ep_sh)
     *	.and.(dabs(ep_s-ep_o1).le.ep_su)) then               
c     Strain hardening branch
        P = dlog10((fp_sh1+fps_su*(ep_su-ep_sh1)-fp_su)/(fp_sh+fps_su*    
     *       (ep_su-ep_sh)-fp_su))/dlog10((ep_su-ep_sh1)/(ep_su-ep_sh))
	if(val2.lt.1.5d0) val2=2.d0
c	

c	  Eq 18.b - This term provides the initial power of the strain hardening curve       
	aux1=s1*(fp_sh+fps_su*(ep_su-ep_sh)-fp_su)
	aux2=(ep_su-s1*(ep_s-ep_o1))/(ep_su-ep_sh)
	  
	  val1 = aux1*(aux2**P)-fps_su*(s1*ep_su-(ep_s-ep_o1))+s1*
     *fp_su



c	  %Eq 18.a

        val3 = fps_su+aux1*P*(aux2**(P-1.d0))*(-1.d0*s1)/(ep_su-ep_sh)
c
        val2=2.d0
	else                                                           
c	Post-ultimate branch, NEEDS DEFINING
        val1 = fp_su
	  if(ep_s.lt.0.d0) val1=-fp_su
c	  
	  val3 = 0.d0
c        val2 = 3.
	endif
	return
	end subroutine virginLoading

c	=============================================================================================	



	
c	=================================================================================================Bauschinger



	subroutine Bauschinger(major,icheat,P_major,P_minor,region,
     *	ep_sushift,s1,fp_su,fps_su,Ep_u,f_y,ep_a,fp_a,ep_s,
     *  k,ep_target,fp_target,fps_target,ep_su,fp_s,fps_s,ep_o,ix,
     *	ep_so,fp_so)

	implicit none
	integer major,alert,k,Bisection,n,icheat,iter,ix(2,1)
	real*8 fp_p,ep_a,fp_a,fp_t,ep_p,ep_su,region,ep_o(2)
	real*8 epp_nl,epp_nu,epp_nm,rm,fp_s,fps_s,ep_rejoin,fp_rejoin
	real*8 P_major,P,ep_target,fp_target,fps_target,tol,tol12
	real*8 fp_su,fps_su,ep_sushift(2),W,f_y,ep_r(2),fp_r(2),Ep_u
	real*8 P_minor,depp_N,epp_N,Rl,Ru,R,J,s1,valn1,ep_s,fnorm1,ep_o1
	real*8 slope2,ep_so,fp_so

	
c	

	if((ep_target-ep_a).eq.0.d0) then
	slope2=fps_target
	else
	slope2=0.8d0*(fp_target-fp_a)/(ep_target-ep_a)
	endif
	if(fps_target.lt.slope2) slope2=fps_target



           
		if(major.eq.0) then
		P=P_minor
		else
			if (icheat.eq.1) then
                P = P_major
		  
			else
                 fp_p = fp_su*(s1-ep_sushift(k)+ep_a)-fp_a             
c			   %Eq 33.a
                 fp_t = fp_su*(2.d0-ep_sushift(1)+ep_sushift(2))       
c	 Eq 33.b
c                 ep_p = abs(ep_sushift(k)-ep_a) 
	ep_p      = dabs((0.2d0*s1+ep_o(k) - ep_a)/0.2d0)              
c	 Eq 33.c
c                 W = (0.001+0.00108/(1.043-ep_p/ep_su))*(fp_p/fp_t-0.69)  
c     *                /0.18+0.085

c	 Eq 34

	fnorm1 = dabs(fp_p/fp_t)
              W = ((0.001d0+0.00108d0/(1.043d0-ep_p))/0.18d0*(fnorm1-
     *	0.69d0)+0.085d0)
              If (W.GT.0.085d0) W = 0.085d0
              If (W.LT.0.05d0) W = 0.05d0

c	W=1.1d0*W
c              Omega = Omega*OmegFac
              P = 56.689d0*(W-0.077d0)**2.d0-4.921d0*(W-0.077d0)+0.1d0 

			endif
		endif

 
        depp_N = 0.d0 
	  R = 1.d0
	  tol = 1.d-9 
	  n = 0 
	  epp_N = 0.1d0 
	  Bisection = 0


   31 epp_N=epp_N+depp_N
	
	if(epp_N.lt.0.d0) goto 32
	


      call rcalc(epp_N,P,ep_target,fp_target,slope2,ep_a,fp_a,Ep_u,
     *	ep_s,R,ix)
      J = -(fp_target-fp_a+Ep_u*(ep_a-ep_target))/(fp_target-fp_a+ 
     *          slope2*(ep_a-ep_target))-P*(1.d0-(epp_N-1.d0)**2.d0
     *)**(P-1.d0)*(2.d0*epp_N-2.d0)

c	First derivative of r w/ respect to epp_N
            depp_N=-R/J
            n = n+1


	
	if((dabs(R).gt.tol).and.n.lt.10) goto 31
	goto 33

   32 continue
	Bisection=1
		
   33 continue
	if(n.eq.10) Bisection=1

        if (Bisection.eq.1) then                                       
c	  %Bisection algorithm
c	
            Ru = 1.d0 
		  Rl = 1.d0 
		  epp_N = 0.d0

	n=0

c	Bracket solution

   34		continue
		n=n+1
		valn1=n
		epp_Nl = 0.d0
		epp_Nu = epp_N+0.001d0*valn1
		call rcalc(epp_Nl,P,ep_target,fp_target,slope2,
     *	ep_a,fp_a,Ep_u,ep_s,Rl,ix)
		call rcalc(epp_Nu,P,ep_target,fp_target,slope2,
     *	ep_a,fp_a,Ep_u,ep_s,Ru,ix)

		if((Ru*Rl).gt.0.d0.and.n.lt.1000) goto 34


	if(n.eq.1000.and.(Ru*Rl).gt.0.d0) then
	write(*,*) 'could not bracket solution!'
	write(*,*) epp_Nl,Rl
	write(*,*) epp_Nu,Ru
	write(*,*) ep_target,fp_target,slope2,ep_a,fp_a,Ep_u,ep_s
	write(100,*) ep_target,fp_target,slope2,ep_a,fp_a,Ep_u,ep_s
				write(*,*) ix(1,1),ix(2,1),ep_s
c	write(*,*) ep_r(1),fp_r(1),ep_r(2),fp_r(2),P
	stop
	endif

            if (dabs(Rl).lt.tol) then                                  
c		  %if the lower bound is a root
                epp_N = epp_Nl
            elseif (dabs(Ru).lt.tol) then                              
c		  %or if the upper bound is a root
                epp_N = epp_Nu
            else                                                       
c		  %otherwise do the bisectional algorithm
                tol = 1.d-10
				tol12=1.d-4
				iter=0



   35				continue
				iter=iter+1
                   epp_Nm = (epp_Nl+epp_Nu)/2.d0
                 call rcalc(epp_Nm,P,ep_target,fp_target,slope2,
     *	ep_a,fp_a,Ep_u,ep_s,Rm,ix)
                    if ((Rl*Rm).lt.0.d0) then 
                        epp_Nu = epp_Nm
                        Ru = Rm
                    elseif ((Rl*Rm).gt.0.d0) then
                        epp_Nl = epp_Nm
                        Rl = Rm
				  else
                    endif

				if((epp_Nu-epp_Nl).gt.tol.and.iter.lt.10000
     *	.and.dabs(Rm).gt.tol12) goto 35

	if(iter.eq.10000) then
	write(*,*) 'maximum number of iterations in Bisection method
     *exceeded!'
	write(*,*) epp_Nu,epp_Nl,Rl,Ru,epp_Nm,Rm,iter
	write(*,*) ' ' 
	write(*,*) ep_target,fp_target,slope2,ep_a,fp_a,Ep_u,ep_s
	write(*,*) ix(1,1),ix(2,1),ep_s
	stop
	endif
                epp_N = epp_Nm
            endif
        endif


c	call rcalc(epp_N,P,ep_target,fp_target,slope2,
c     *	ep_a,fp_a,Ep_u,ep_s,Rm,ix)

c	

        fp_s = (epp_N*((fp_target-fp_a)-Ep_u*(ep_target-ep_a))+Ep_u*    
     *            (ep_s-ep_a)+fp_a)
c	



	if (epp_N.lt.0.0001d0.OR.(Ep_u-slope2)/Ep_u.LT.0.01d0) then
        fps_s = Ep_u
      else
        fps_s = 2.d0*P*(1.d0-(1.d0-epp_N)**2.d0)**(P-1.d0)*(1.d0-epp_N) 
        fps_s =fps_s*((fp_target-fp_a)-slope2*(ep_target-ep_a))/
     *	  (((ep_target-ep_a)*Ep_u-(fp_target-fp_a))/(Ep_u-slope2))
        fps_s = fps_s*Ep_u/(fps_s+Ep_u) + slope2
      endif

 

c	
c	write(*,*) 'fp_s=', fp_s,'fps_s=',fps_s,Ep_u
	if(fps_s.lt.0.) then
c	write(401,*) 'NEGATIVE SLOPE!',ix(1,1),ix(2,1)
c	write(*,*) ep_a,ep_target,fps_s,P
	fps_s=(fp_s-fp_so)/(ep_s-ep_so)
c	STOP
	endif


c	tangent slope of curve
	return
	end subroutine Bauschinger

c	====================================================================================================


c	===================================================================================================== rcalc
	subroutine rcalc(epp_N,P,ep_target,fp_target,fps_target,ep_a,
     *	fp_a,Ep_u,ep_s,R,ix)
	implicit none
	integer ix(2,1)
	real*8 epp_N,P,ep_target,fp_target,fps_target,R
	real*8 fp_a,Ep_u,ep_a,ep_s,val1,val2


	val2=(fp_target-fp_a)-fps_target*(ep_target-ep_a)
c	val2=fp_target-fp_a-fps_target*ep_target+fps_target*ep_a

	val1=epp_N*((fp_target-fp_a)-Ep_u*(ep_target-ep_a))
     *	+(Ep_u-fps_target)*(ep_s-ep_a)


c	val1=epp_N*fp_target-epp_N*fp_a-epp_N*Ep_u*ep_target+
c     *epp_N*Ep_u*ep_a+Ep_u*ep_s-Ep_u*ep_a-fps_target*ep_s+
c     *fps_target*ep_a
	
	
	


	if(val2.eq.0.d0) then
	write(*,*) 'ERROR!'
	write(*,*) ix(1,1),ix(2,1),ep_s
	stop
	endif

	R=1.d0-epp_N
	R=R**2.d0
	R=1.d0-R
	R=R**P
	R=R-val1/val2
c        R = (1.-(1.-epp_N)**2.)**P-val1/val2
c	  %Eq 29, moved all to the righthand side to solve for r=0
	return
	end subroutine rcalc
