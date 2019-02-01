
c  =====================================================================
c  sdmuc = interface subroutine between OpenSees and SD model
c  =====================================================================

      subroutine sdmuc (strhs, strsg, props, stran, nmats, nstrp,
     &                  istep, iiter, ielem, strhs0, etahs, hdp, oths)

      implicit none
      integer,parameter::dbl=selected_real_kind(15,307)
      integer,intent(in)::nmats
      integer,intent(in)::nstrp
      integer,intent(in)::istep
      integer,intent(in)::iiter
      integer,intent(in)::ielem
      real(kind=dbl),dimension(4),intent(inout)::strsg
      real(kind=dbl),dimension(4),intent(inout)::stran
      real(kind=dbl),dimension(nmats),intent(inout)::props
      real(kind=dbl),dimension(nstrp),intent(inout)::strhs
      real(kind=dbl),dimension(280),intent(inout)::strhs0
      real(kind=dbl),dimension(3,40),intent(inout)::etahs
      real(kind=dbl),dimension(3,80),intent(inout)::hdp
      real(kind=dbl),dimension(12),intent(inout)::oths

      call dsmod(strsg, stran, props, nmats, strhs, nstrp, etahs, hdp,
     &           oths, strhs0)

      end subroutine sdmuc

c=======================================================================
c
c     
c        SDM-sub contains the subroutines of S-D Model
c        =======
c
c
c       Includes: dsmod, clamda, hsf, calhp1, cnjpnt, psspar,
c       --------  angle, esspar, strinc, strcon, fxinc, tsinc, threed
c
c      -----------------------------------------------------------------
c
       subroutine  dsmod ( sig, deps, props, nmats, strhs, nstrp,
     *                     etahs, hdp, oths, strhs0 )
c                                                                       
c                          ... main subroutine of S-D Model
c
c      -----------------------------------------------------------------
c                                                                      
c
       implicit real*8  ( a-h, o-z )                            
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3,igaus
       common / axil  / svin, rf0in, pwpr, ainp
       common / elpar / ae, be, xgi, xg2, eyng, pora, xn, coef2
       common / strn  / eta(3), etacum(3), etarev(3), etar(3), etad(3)
       common / cal   / ical
       common / tstr  / fmuf, depsx, depsy
       common /hdpp/ x(3,20), y(3,20), cx(3,20), cy(3,20)
c
       dimension   props(nmats), strhs(nstrp),
     1             sig(4), dsig(4), eps(4), deps(4),
     1             depse(4), depsp(4), ddeps(4),
     2            hdp(3,80), etahs(3,40), oths(12),
     3             iflag(3), rna(3), coef6(3), strhs0(280)
c                                                                       
c
 1967  format(3i5,3(e12.5,2x))
       depsx  = deps(2) - deps(1)                                       
       depsy  = deps(3)*2.d0                                            
       if ( dabs(deps(1)).lt.1.d-10 .and. dabs(deps(2)).lt.1.d-10
     1      .and. dabs(deps(3)).lt.1.d-10 ) return
c
       pi     = 3.141592654d0                                           
       refeta = 0.01d0

       iitr = oths(11)
       istp = oths(12)
c                                                                      
       do 5 i = 1,4                                                     
       dsig(i)  = 0.d0                                                  
       depsp(i) = 0.d0
       depse(i) = 0.d0
       ddeps(i) = 0.d0
       if ( istp.eq.1 .and. iitr.eq.0 ) eps(i) = 0.d0
    5  continue                                                        
c                                                                      
       ec     = props(11) / ( 1.d0 - props(11))
       prefe  = 98.1d0
       fact   = props(38)                                               
       c1     = 0.d0
       c2     = 0.d0
       fmuc   = props(27)
       fmuf   = props(29)                                               
       fmyu0  = props(30)                                              
       sc     = props(31)                                               
c                                                                       
          eqzxp    = oths(1)
	  eqzx     = oths(2)
	  eqyzp    = oths(3)
	  eqyz     = oths(4)
          svin     = oths(5)                            
	  if ( oths(8).gt.1.d-05 ) 
     1                               props(6) = oths(8)
c
c      fmuf1  = oths(6,ielem)                               
c      rvol   = oths(7,ielem)                           
c
         do 95 i = 1,3
          iflag(i) = int(strhs(i)+.1)                     
 95       rna(i)   = strhs(i+3)
c
          eta(1)    = strhs(7)
          eta(2)    = strhs(8)
          eta(3)    = strhs(9)
          etacum(1) = strhs(10)
          etacum(2) = strhs(11)
          etacum(3) = strhs(12)
          etarev(1) = strhs(13)
          etarev(2) = strhs(14)
          etarev(3) = strhs(15)
          etad(1)   = strhs(16)
          etad(2)   = strhs(17)
          etad(3)   = strhs(18)
          ceta      = strhs(19)
          fis       = strhs(20)
          xn        = strhs(21)
          etar(1)   = strhs(22)
          etar(2)   = strhs(23)
          etar(3)   = strhs(24)
          rvol      = strhs(25)
c
c
       k = 1
	do 65 i = 1,3  
	do 65 j = 1,80
	  m = (i-1)*80+40+j
	  if (j .le. 20)             x(i,k)  = strhs0(m)
	  if (j.gt.20 .and. j.le.40) y(i,k)  = strhs0(m)
	  if (j.gt.40 .and. j.le.60) cx(i,k) = strhs0(m)
	  if (j.gt.60 .and. j.le.80) cy(i,k) = strhs0(m)
        k = k + 1
	if (k .eq. 21) k = 1
 65     continue
c
c
        if ( eta(1).le.0.00100d0 .and. iflag(1) .gt. 1 ) then
            coef2 = 2.d0
          elseif ( eta(1).gt.0.00100d0 .and. eta(1).lt.0.00300d0 .and. 
     1             iflag(1).gt.1 ) then
	    coef2 =  2.d0 - (( eta(1) - 0.00100d0 ) / 0.00200d0 ) 
          else
            coef2  = 1.d0
	end if
c
ccc21JAN16
       if  (iflag(1).eq.1 .and. etarev(1).lt.1.d-05) then
          fmu0   = fmyu0
	if ( eta(1).le.0.000100d0 ) fmu0 = 0.05d0
	if (eta(1) .gt. 0.000100d0 .and. eta(1) .lt. 0.0011d0) then
	 fmu0 =  0.05d0 + (( eta(1) - 0.000100d0 ) / 0.00100d0 ) * fmyu0
        end if
          coef6(1)  = 1.0d0
c      elseif ( (iflag.eq.1 .and. etarev.gt.1.d-05) .or.
c    1         (iflag.eq.1 .and. rvol.gt.0.25d0) ) then
c            fmu0   = 0.d0
c            coef6  = 1.0d0
c      elseif (iflag.gt.1) then
       else
          fmu0   = fmuc
c      if (rvol.le.0.15d0) fmu0 = 0.d0
c      if (rvol.le.0.15d0) fmu0 = fmuc/2.d0
          coef6(1)  = 2.d0
       end if
ccc21JAN16
c
        rvolc = dsqrt( x(1,1)**2 + y(1,1)**2 )
	if (rvol.lt.rvolc) rvol = rvolc
c
c
c     write(*,*)'istep',istep
c     write(*,*)'eta(1)',eta(1)
c     write(*,*)'ceta',ceta
c     write(*,*)'etar(1)',etar(1)
	  do 85 i = 1,3
	  do 85 j = 1,20
	     if ( iflag(i).eq.j.and.eta(i).lt.etahs(i,j)) then
c        if ( iflag(i).ge.1 ) eta(i) = etahs(i,j,ielem)
            if ( iflag(i).gt.1 ) eta(i) = etahs(i,j)
                etar(i) = etahs(i,j+20)
         end if
 85   continue
c
c
       call strinc ( depsx, depsy, incrmt, ddeps, deps )
c                                                                       
  100  incfai = incfai + 1                             
c                                                                      
       p = (sig(1) + sig(2))*0.5d0
c
       call esspar ( props, nmats, ec, prefe, p, fis, oths )
       call psspar ( props, nmats, sig, ec, prefe, si,  p,
     1               ggmin, ggmax, rf0, fis )
       call strcon ( sig, c1, c2, rf0, pk, 1 , ialarm ,iala )
c                                                                      
       bigx   = 0.5d0*( sig(2) - sig(1) )                       
       bigy   = sig(3)                                                  
       bigxp  = bigx / pk                                        
       bigyp  = bigy / pk                                  
       qk     = dsqrt( bigx**2 + bigy**2 )                              
       qkk    = dsqrt( (bigx - c1*pk)**2 + (bigy - c2*pk)**2 )   
c
c      szxnp  = 0.5d0*( sig(3) - sig(1) )
c      szxtp  = sig(6) / pk
c      syznp  = 0.5d0*( sig(2) - sig(3) )
c      syztp  = sig(5) / pk
c                                                                       
       call fxinc ( ddeps, pk, pi, bigxp, bigyp, c1, c2, rf0, rna,
     1              dbigx, dbigy, dbigxp, dbigyp, eangp, pangp )
c                                                                       
cccc
c      write(*,*)'bigx,bigy',bigx,bigy
c      write(*,*)'dbigx,dbigy',dbigx,dbigy
       if ( (bigx**2 + bigy**2) .gt. 1.d-10 ) go to 15         
       if ( (dbigx**2 + dbigy**2) .lt. 1.d-10 ) go to 35          
       bs = datan2(dbigy,dbigx)                                 
       go to 45                                                         
c                                                                       
   35  bs     = 0.d0                                                 
c      write(*,*)'went to 35'
	do 125 i = 1,3
 125    rna(i) = 0.d0                                                
       go to 55                                                         
c                                                                       
   15  bs = datan2(bigy,bigx)                                           
   45  continue                                                         
       if (bs.lt.0.d0) bs = bs + 2.d0*pi                    
c                                                                       
         iffl = iflag(1)
c      write(*,*)'eta(1)',eta(1)
c      write(*,*)'ceta',ceta
c      write(*,*)'etar(1)',etar(1)
c
       call hsf ( bigxp, bigyp, dbigxp, dbigyp, iflag, rna,
     1            cxx, cyy, rvol, rf0, etahs, ceta, 1 )
	 if (iflag(1).gt.(iffl+1)) write (6,3289)iffl,iflag,ielem,istep
 3289  format ('Warning in S-D Model #dsmod# iffl>iflag;
     1          iffl,iflag,ielem,istep =',4i5)
c                                                                       
c
   55   continue
	   do 135 i = 1,3
ccccc       if ( rna(i).lt.1.d-5 ) then
       if ( rna(1).lt.1.d-5 ) then
       if (ceta .gt. etarev(i)) then
          etarev(i) = ceta
          etacum(1) = ( etarev(1) - 0.020d0 ) / 0.06d0
c          write(*,*)'etarev,etacum',etarev(1),etacum(1)
       end if
          if ( etacum(1) .gt. 1.d0 ) etacum(1) = 1.d0
       eta(i)    = 0.d0                                    
       ceta      = 0.d0
       etar(i)   = 0.d0
ccc21JAN16
       if  (iflag(1).eq.1 .and. etarev(1).lt.1.d-05) then
          fmu0   = fmyu0
	if (eta(1).le.0.000100d0 ) fmu0 = 0.05d0
	if (eta(1) .gt. 0.000100d0 .and. eta(1) .lt. 0.0011d0) then
	     fmu0 =  0.05d0+((eta(1) - 0.000100d0) / 0.00100d0)*fmyu0
        end if
          coef6(1)  = 1.0d0
c      elseif ( (iflag.eq.1 .and. etarev.gt.1.d-05) .or.
c    1         (iflag.eq.1 .and. rvol.gt.0.25d0) ) then
c            fmu0   = 0.d0
c            coef6  = 1.0d0
c      elseif (iflag.gt.1) then
       else
          fmu0   = fmuc
c      if (rvol.le.0.15d0) fmu0 = 0.d0
c      if (rvol.le.0.15d0) fmu0 = fmuc/2.d0
          coef6(1)  = 2.d0
       end if
ccc21JAN16
       end if
 135   continue
c
c
       qh = rna(1)*pk                                      
c                                                                      
       go to (9011,9012,9012,9012,9012,9012,9012,9012,9012,9012,
     1        9012,9012,9012,9012,9012,9012,9012,9012,9012,9012,
     2        9012) iflag(1)
 9011  continue
       xh = bigx - c1*pk                                               
       yh = bigy - c2*pk                                                
       go to 9015                                                       
 9012  xh = bigx - cxx*pk                                              
       yh = bigy - cyy*pk                                             
 9015  if (qh.lt.1.0d-10) then                       
c
       fhp = dabs(y(1,1))                                             
       if (fmuf.gt.0.15d0. and. dabs(bigy).gt.1.0d-10. and.
     1   dabs(deps(3)).gt.1.0d-10) 
     1 fhp = - dabs(y(1,1))*((bigy/dabs(bigy))*(deps(3)/dabs(deps(3))) )
       go to 9016                                                     
       end if                                                           
c                          
       fhp = ( -bigx*xh - bigy*yh ) / (qh*pk)                           
 9016  continue                                                         
c
	 ft = fact
	 et = eta(1)
c                                                                     
c
c      Back-bone curve
c122   if (iflag(1) .eq. 1) then
       if (iflag(1) .eq. 1) then
ccc21JAN16       	  coef6(1) = 1.d0
c  if ( rvol .le. 0.50d0 ) then
ccc21JAN16	     fmu0 = fmyu0
	     ceta = eta(1)
c         end if
c  if ( rvol .gt. 0.50d0 ) fmu0 = fmuc
       else
	  coef6(1) = 2.d0
	  fmu0     = fmuc
       end if
c
          if ( ((etar(1)/refeta)*fact) .gt. 20.d0 ) then
             gpls = ggmin
             go to 8102
          end if
c
       gpls  = (ggmax - ggmin)*dexp(-etar(1)/(coef6(1)*refeta)*fact)
     1          + ggmin   
8102   gplss = gpls - fact*(etar(1)/(coef6(1)*refeta))*(gpls - ggmin) 
c      write(*,*)'iflag(1)',iflag(1),incfai
c      write(*,*)'coef6(1)',coef6(1)
c      write(*,*)'eta(1)',eta(1)
c      write(*,*)'etar(1)',etar(1)
c      write(*,*)'gpls',gpls
c
       if ( fmuf .gt. 0.15d0 ) then
         hp = gplss*pk*(1.d0 - rna(1)/rf0)**2           
         if (iflag(1).eq.21) hp=pk*2.d0*rna(1) /etahs(1,40)
       elseif ( fmuf.le.0.15d0 ) then
         hp = gplss*pk*(1.d0 - rna(1)/rf0)**2 
         if (iflag(1).eq.21) hp=pk*2.d0*rna(1) /etahs(1,20)
       end if
         if ( hp .gt. ggmax*pk ) hp = ggmax*pk
c
       rc  = fmu0 + 2.d0/pi*(fmuf - fmu0)*datan(eta(1)/(sc*coef2))    
c      write(*,*)'hp',hp
c      write(*,*)'coef2',coef2
c      write(*,*)'fmu0',fmu0
c      write(*,*)'rc',rc
c      write(*,*)'fmuf',fmuf
c      write(*,*)'sc',sc
c      write(*,*)'ggmax',ggmax
c      write(*,*)'ggmin',ggmin
c      write(*,*)'rna(1)',rna(1)
c      write(*,*)'rf0',rf0
c      write(*,*)'pk',pk
c      write(*,*)'etarev(1)',etarev(1)
c      write(*,*)'ceta',ceta
c      write(*,*)'etahs',etahs(1,1),etahs(1,2),etahs(1,3),etahs(1,4)
c      write(*,*)'etahs',etahs(1,21),etahs(1,22),etahs(1,23),etahs(1,24)
c      write(*,*)' '
c
       iala = 0
 222   call clamda ( ddeps, depse, depsp, dsig, bigx, bigy, pk, qk, qkk,
     1               c1*pk, c2*pk, rf0, rc, hp, bs, eangp, pangp, fhp,
     2               rna, coax, rvol, gplss, iflag, ireload,strhs,nstrp)
c       if (ireload.eq.1) then
c	 do 195 kjk = 1,3
c          iflag(kjk) = int(strhs(kjk)+.1)                     
c	  eta(kjk)   = strhs(kjk+6)
c195       rna(kjk)   = strhs(kjk+3)
c       go to 122
c       end if
c                                     
c
       if (sig(1) .gt. 1.d-11) sig1 = sig(1)
       if (sig(2) .gt. 1.d-11) sig2 = sig(2)
       if (sig(1) .gt. 1.d-11 .and. sig(2). gt. 1.d-11) sig3 = sig(3)
c
c
       do 25 i = 1,3                                                    
       sig(i)  = sig(i) + dsig(i)                                      
       eps(i)  = eps(i) + ddeps(i)                                      
   25  continue                                                         
c                                                                     
       call strcon (sig, c1, c2, rf0, ppk, 2 , ialarm , iala )
       
       if ( ialarm .eq. 1 ) then
         sig (1) = sig1
         sig (2) = sig2
         sig (3) = sig3
c          write (14,*) ielem, iiter, fhp
       end if
c                                                                      
       eta(1)  = eta(1) + dsqrt ( ( ddeps(2)-ddeps(1) )**2
     &           + 4.d0*ddeps(3)**2 ) 
       ceta    = ceta + dsqrt ( ( ddeps(2)-ddeps(1) )**2
     &           + 4.d0*ddeps(3)**2 ) 
c      eta(2)  = eta(2) + dsqrt ( ( ddeps(2)-ddeps(3) )**2
c    &           + 4.d0*ddeps(5)**2 ) 
c      eta(3)  = eta(3) + dsqrt ( ( ddeps(3)-ddeps(1) )**2
c    &           + 4.d0*ddeps(6)**2 ) 
c
       etar(1) = etar(1) + dsqrt( ( depsp(2) - depsp(1) )**2
     &        + 4.d0*depsp(3)**2 ) 
c      etar(2) = etar(2) + dsqrt( ( depsp(2) - dezp )**2 + dgyp**2 )
c      etar(3) = etar(3) + dsqrt( ( dezp - depsp(1) )**2 + dgzp**2 )
c
c      etacum(1) = etacum(1) + 1.d0/incrmt*dsqrt( (deps(2)-deps(1) )**2
c    &        + 4.d0*deps(3)**2 )                                      
c      etacum(2) = etacum(2) + 1.d0/incrmt*dsqrt( (deps(2)-deps(3) )**2
c    &        + 4.d0*deps(5)**2 )                                      
c      etacum(3) = etacum(3) + 1.d0/incrmt*dsqrt( (deps(3)-deps(1) )**2
c    &        + 4.d0*deps(6)**2 )                                      
c                                                                       

c       write(*,*)'sig',sig
c       write(*,*)'eps',eps
c       write(*,*)'ddeps',ddeps
c       write(*,*)'depsp',depsp
c       write(*,*)'depse',depse
c       write(*,*)'eta(1)',eta(1)
c       write(*,*)'ceta',ceta
c       write(*,*)'etar(1)',etar(1)
c       write(*,*)' '
 188   if (incfai.lt.incrmt) go to 100                                 
          if ( ialarm .eq. 1) then
            go to 199
          end if
c                                                                       
c
          oths(1)   = eqzxp
	  oths(2)   = eqzx
	  oths(3)   = eqyzp
	  oths(4)   = eqyz
	  oths(5)   = svin
c
c      oths(6)  = fmuf1
c      oths(7)  = rvol
	  oths(8)   = props(6)
c                                                                       
         do 96 i = 1,3
          strhs(i)   = dfloat(iflag(i))
 96       strhs(i+3) = rna(i)
c
          strhs(7)   = eta(1)
          strhs(8)   = eta(2)
          strhs(9)   = eta(3)
          strhs(10)  = etacum(1)
          strhs(11)  = etacum(2)
          strhs(12)  = etacum(3)
          strhs(13)  = etarev(1)
          strhs(14)  = etarev(2)
          strhs(15)  = etarev(3)
          strhs(16)  = etad(1)
          strhs(17)  = etad(2)
          strhs(18)  = etad(3)
          strhs(19)  = ceta
          strhs(20)  = fis
          strhs(21)  = xn
          strhs(22)  = etar(1)
          strhs(23)  = etar(2)
          strhs(24)  = etar(3)
          strhs(25)  = rvol
c
c                                                                       
       k = 1
	do 66 i = 1,3  
	do 66 j = 1,80
	  if (j .le. 20)             hdp(i,j) = x(i,k)
	  if (j.gt.20 .and. j.le.40) hdp(i,j) = y(i,k)
	  if (j.gt.40 .and. j.le.60) hdp(i,j) = cx(i,k)
	  if (j.gt.60 .and. j.le.80) hdp(i,j) = cy(i,k)
        k = k + 1
	if (k .eq. 21) k = 1
 66     continue
c
c                                                                       
 199   continue
c
       return                                                           
       end                                                              
c
c=================================================================== 
       subroutine clamda ( deps, depse, depsp, dsig, bigx, bigy, pk,
     1                     qk, qkk, c11, c22, rf0, rc, hp, bs, eang,
     2   pang, fhp, rna, coax, rvol, gplss, iflag, ireload,strhs,nstrp)
c
c                      ... plastic strain and stress increments
c                                                             
c      -------------------------------------------------------------
c                                                                      
c
       implicit real*8 ( a-h, o-z )     
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3, igaus
       common / elpar / ae, be, xgi, xg2, eyng, pora, xn, coef2
       common / tstr  / fmuf, depsx, depsy
       common / cal   / ical
       dimension deps(4), dsig(4), depse(4), depsp(4), rna(3), iflag(3),
     1     strhs(nstrp)
c                                                                     
       dsigx = 0.d0
       dsigy = 0.d0
       dsigxy = 0.d0
       pi    = 3.141592654d0                                        
       bmin  = 2.d0*pi                                               
       bmin1 = bmin                                                    
       isw   = 0                                          
       bigx0 = bigx/pk                                                
       bigy0 = bigy/pk                                                
       c1    = c11/pk
       c2    = c22/pk
	 if ( fmuf.lt.0.15d0 ) then
	 dsig(3) = hp * depsy
         dsig(2) = hp * depsx
	 dsig(1) =-hp * depsx
c
       pkk   = pk + 0.5d0*(dsig(2) + dsig(1))                         
       bigx1 = ( bigx + 0.5d0*(dsig(2)-dsig(1)) ) / pkk                 
       bigy1 = ( bigy + dsig(3) ) / pkk                                 
c       rvo = dsqrt( (bigx0)**2 + (bigy0)**2 )         
c       rvn = dsqrt( (bigx1)**2 + (bigy1)**2 )       
c	if (rvn.gt.rvo .and. iflag(1).gt.(int(strhs(1)+.1)) .and.
c     & ireload.eq.0) then
c	 ireload = 1
c	 return
c	end if
c	 ireload = 0
	 return
	 end if
c                                                                     
c                                                                     
c      number of steps for the loop                                    
c      ----------------------------                                    
c                                                                       
       angd  = dabs(eang - pang)                                       
	angl = angd
	 if (angd.gt.pi) then
	   if ( eang .gt. pang ) then
             angl = 2.d0*pi + pang - eang
	     astr = eang
           elseif ( pang .gt. eang ) then
	     angl = 2.d0*pi - pang + eang
	     astr = pang
           end if
	 end if
       angdd = (angl/pi)*180.0d0                                       
       iiii  = idint(angdd)                                            
       if (iiii.lt.1) iiii = 1                                          
       go to 202                                                       
c                                                                     
 101   isw  = 101                                                      
       iiii = 10
 202   if (isw .ne. 101) bbdd  = ( pang + eang ) / 2.d0
 102   continue
c                                                                       
c                                   -------   start - loop   --------
       do 50 ii=1,iiii                                                  
c                                                                     
       if (isw.eq.101) go to 303                                        
c                                                                     
c      assumed direction of plastic strain increment                    
c      ---------------------------------------------                  
c                                                                     
       if (pang.le.eang .and. angd.lt.pi) 
     1      bbb = pang + (1.d0/180.0d0)*pi/5.d0 +        
     2                              ((ii-1)/180.0d0)*pi               
       if (pang.gt.eang .and. angd.lt.pi ) 
     1      bbb = eang + (1.d0/180.0d0)*pi/5.d0 +        
     2                              ((ii-1)/180.0d0)*pi            
       if (iiii.eq.1)    bbb = (pang + eang) / 2.0d0                    
c
       if (angd.ge.pi)   
     1      bbb = astr + (1.d0/180.0d0)*pi/5.d0 +
     2                         ((ii-1)/180.0d0)*pi
       go to 505                                                
c                                                                       
 303   bbb = bbdd- (1.d0/180.0d0)*pi + ii*(1.d0/180.0d0)*pi/5.d0
       if  (bbb.lt.0.d0) bbb = 0.d0
 505   if ( bbb .gt. 2.d0*pi ) bbb = bbb - 2.d0*pi
c                                                                       
c      assumed conjugate point and 2betadeps                           
c      ------------------------------------                          
c                                                                     
       xc = rf0*pk*dcos(bbb) + c11                                      
       yc = rf0*pk*dsin(bbb) + c22                                    
c       write(*,*)'xc,yc',xc,yc
c                                                                       
       if (dabs(xc).lt.1.d-5) then                                    
       bde = (yc/dabs(yc))*pi/2.0d0                                    
       go to 40                                                       
       end if                                                           
c                                                                       
       bde = datan(yc/xc)                                             
       if (xc.lt.0.0d0)  bde = bde + pi                                 
   40  if (bde.lt.0.0d0) bde = bde + 2.d0*pi                            
c                                                                     
c      noncoaxiality                                                    
c      -------------                                                 
c                                                                       
       psi  = dabs(bbb-bs)
       if (psi .gt. pi)  psi = 2.d0*pi - psi                           
       c    = dcos(psi)                                                 
       coax = c
c       write(*,*)'angd',angd
c       write(*,*)'eang,pang',eang,pang
c       write(*,*)'bbb',bbb
c       write(*,*)'bs',bs
c       write(*,*)'psi',psi
c       write(*,*)'coax',coax
c       write(*,*)'fhp',fhp
c       write(*,*)'qk',qk
c       write(*,*)'pk',pk
c                                                                      
c      calculation of lamda                                             
c      --------------------                                            
c                                                                      
       rr   = qkk/pk                                                   
       qf   = dsqrt ( (xc - c11)**2 + (yc - c22)**2 )                   
       if ( fmuf.lt.0.15d0 .or. iflag(1).eq.21) then
	  rc = rr*c
	  end if
       ccx  = 0.500d0*(rc - rr*c) - 0.5d0*(xc - c11)/qf                 
       ccy  = 0.500d0*(rc - rr*c) + 0.5d0*(xc - c11)/qf        
       ccxy = 0.5d0*yc/qf                                             
c
c      if ( fmuf.lt.0.15d0 .or. iflag(1).eq.21 ) then
ccc     1      (rvol.lt.0.1d0 .and. rr.lt.0.1d0) ) then
c  ccx = 0.d0
c      ccy = 0.d0
c      end if
c                                                                       
       aa = 0.5d0*( fhp - (xc - c11)/qf )                             
       bb = 0.5d0*( fhp + (xc - c11)/qf )                              
       cc = yc /qf                                                    
c                                                            
       a1 = aa*ae + bb*be                                             
       a2 = aa*be + bb*ae                                            
       a3 = cc*xg2                           
c                                                                     
       xlamda = (a1*deps(1) + a2*deps(2) + a3*deps(3)) /              
     *          (a1*ccx     + a2*ccy     + a3*ccxy + hp)               
c
c
       if(xlamda.lt.0.d0) then                                        
c                                                                     
       do 15 k = 1,3                                                  
   15  dsig(k) = 0.d0                                                  
       go to 50                                                       
c                                                                      
       end if                                                          
c                                                                    
c      plastic strain increment                                       
c      ------------------------                                       
c                             
       depsp(1) = xlamda*ccx                                          
       depsp(2) = xlamda*ccy                                           
       depsp(3) = xlamda*ccxy                                          
c       write(*,*)'xlamda',xlamda
c       write(*,*)'deps',deps
c       write(*,*)'depsp',depsp
c       write(*,*)'ccx,ccy,ccxy',ccx,ccy,ccxy
c                                                                      
c      elastic strain increment                                       
c      ------------------------                                       
c                                                                     
       do 10 i = 1,3                                                  
       depse(i) = deps(i) - depsp(i)                                   
   10  continue                                    
c                                                                     
c      stress increment                                                
c      ----------------                                              
c                                                                     
       dsig(1) = ae*depse(1) + be*depse(2)                             
       dsig(2) = be*depse(1) + ae*depse(2)                             
       dsig(3) = xg2*depse(3)                                         
       if ( fmuf.lt.0.15d0 ) dsig(3) = hp*2.d0*deps(3)
c      if ( fmuf.lt.0.15d0 ) then
c      dsig(3) = hp*2.d0*deps(3)
c         eyngt = 2.d0*hp*(1.d0 + pora)
c      aet = (eyngt*(1.d0 - pora)) / ((1.d0+pora)*(1.d0-2.d0*pora))
c      bet = (eyngt*pora) / ((1.d0+pora)*(1.d0-2.d0*pora))
c      dsig(1) = aet*deps(1) + bet*deps(2)               
c      dsig(2) = bet*deps(1) + aet*deps(2)         
c      end if
c
c      stresses in x/p - y/p plane                                      
c      ---------------------------                                    
c                                                                      
       pkk   = pk + 0.5d0*(dsig(2) + dsig(1))                         
       bigx0 = bigx/pk                                                
       bigy0 = bigy/pk                                                
c the above to lines are inserted based on the 3-D version (April2005)
       bigx1 = ( bigx + 0.5d0*(dsig(2)-dsig(1)) ) / pkk                 
       bigy1 = ( bigy + dsig(3) ) / pkk                                 
       dxc   = bigx1 - bigx0                                         
       dyc   = bigy1 - bigy0                                          
c       rvo = dsqrt( (bigx0)**2 + (bigy0)**2 )         
c       rvn = dsqrt( (bigx1)**2 + (bigy1)**2 )       
c                                                                       
c      calculated conjugate point                                      
c      --------------------------                                       
c                                                                       
       call cnjpnt (bigx0,bigy0,dxc,dyc,c1,c2,rf0,xc1,yc1,bde1,bde1p,pi,
     &              3)                                
c                                                                     
c      calculated vs. assumed conj. point                 
c      ----------------------------------                  
c                                                                     
       bdel = dabs(bde - bde1)                                          
       if (bdel.gt.pi) bdel = 2.0d0*pi-bdel                           
c                                                                       
c      storred stress and strain increments                           
c      ------------------------------------                            
c                                                
       if (isw.eq.101) go to 404                                        
c                                                                     
       if (bdel.lt.bmin1)  then                                         
         bmin1 = bdel                                             
         bbdd  = bde                                            
         xlam  = xlamda                                       
       end  if                                                        
c                                                                    
       if ( bdel.gt.bmin1) then
         go to 101                       
       end if                                                        
       go to 50                                                        
c                                                                       
 404   if (bdel.le.bmin)  then                                         
c                                                                      
       bmin  = bdel                                                     
       xlam  = xlamda                                          
       ccxwr = ccx
       ccywr = ccy
c                                                                       
       dsigx  = dsig(1)                                                 
       dsigy  = dsig(2)                                                 
       dsigxy = dsig(3)                                               
c                                                                      
       depspx = depsp(1)                                              
       depspy = depsp(2)                                                
       depspt = depsp(3)                                              
c                                                                      
       depsex = depse(1)                                              
       depsey = depse(2)                                                
       depset = depse(3)                                                
c                                                                      
       end  if                                                         
   50  continue                                                        
c
c                            --------  end loop  --------
c
       if (isw.ne.101) go to 101                                    
c
c
c      calculated stress and strain increments - output                
c      ------------------------------------------------                 
c                                                                      
c
       dsig(1) = dsigx                                              
       dsig(2) = dsigy                                               
       dsig(3) = dsigxy                                              
c       bigx0 = bigx/pk                                                
c       bigy0 = bigy/pk                                                
c     comment included for the above two lines (April2005)
c      pkk   = pk + 0.5d0*(dsig(2) + dsig(1))                         
c      bigx1 = ( bigx + 0.5d0*(dsig(2)-dsig(1)) ) / pkk                 
c      bigy1 = ( bigy + dsig(3) ) / pkk                                 
c       rvo = dsqrt( (bigx0)**2 + (bigy0)**2 )         
c       rvn = dsqrt( (bigx1)**2 + (bigy1)**2 )       
c	if (rvn.gt.rvo .and. iflag(1).gt.(int(strhs(1)+.1)) .and.
c     & ireload.eq.0) then
c	 ireload = 1
c	 return
c	end if
c	ireload = 0
c                                                                     
       depsp(1) = depspx                                              
       depsp(2) = depspy                                             
       depsp(3) = depspt                                              
c                                                                     
       depse(1) = depsex                                                
       depse(2) = depsey                                                
       depse(3) = depset                                              
c                                                                       
       return                                                          
       end                                                              
c                                                                     
c==================================================================
       subroutine hsf ( bx, by, dbx, dby, iflag, rna,
     2                  cxx, cyy, rvol, rf0, etahs, ceta, i )
c                                                                      
c                  ... calculation of the hardening surfaces
c
c      ------------------------------------------------------------
c                                                                       
c
        implicit real*8 ( a-h, o-z )                     
        common / cal   / ical
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3, igaus
        common / strn  / eta(3), etacum(3), etarev(3), etar(3), etad(3)
        common / hdpp  / x(3,20), y(3,20), cx(3,20), cy(3,20)
c
        dimension etahs (3,40), rna(3), iflag(3)
c                                                                    
c
c       next stress state                                            
c       -----------------                                             
c                                                                       
c
	do 77 j = iflag(i),20
	x(i,j) = 0.d0
        y(i,j) = 0.d0
	cx(i,j) = 0.d0
  77    cy(i,j) = 0.d0
c
        bnx = bx + dbx                                             
        bny = by + dby                                              
c                                                                     
        if (iflag(i).eq.0) iflag(i) = 1
	k   = iflag(i)
c                                                                      
       if (iflag(i) .eq. 1 ) then
	  go to 10
       elseif (iflag(i) .gt. 1 .and. iflag(i) .lt. 21 ) then
	  go to 20
       else
          rna(i) = dsqrt((x(i,k-1)-cx(i,k-1))**2 + 
     1                              (y(i,k-1)-cy(i,k-1))**2)
          go to 210
       end if
c                                                                     
c       virgin loading                                               
c       --------------                                                
c                                                                       
  10    rvo = dsqrt( (bx -cx(i,1))**2 + (by -cy(i,1))**2 )         
        rvn = dsqrt( (bnx-cx(i,1))**2 + (bny-cy(i,1))**2 )       
c        write(*,*)'virgin'
c                                                                       
        rna(i) = rvo                            
c
ccc        rvol = rvo
ccccccccc    if (rvo.ge.rf0) write (6,*) 'RVO greater than RF0'
ccccccccc    if (rvn.ge.rf0) then
ccccccccc       write (6,*) 'RVN > RF0'
ccccccccc       write (6,*) ielem, istep, rvn, rf0
ccccccccc       dsg = dsqrt(dbx**2 + dby**2)
ccccccccc       write (6,*) dsg,rvo, rf0
ccccccccc    end if
            if(rvn.lt.rvo) go to 110                           
cccccccc        if( rvn.lt.rvo .and. rvn.lt.rf0 ) go to 110           
c                                                                      
c       back-bone curve
ccc21JAN16	if (iflag(1) .eq. 2) ceta = eta(i)
ccc21JAN16        if (iflag(i) .eq. 2) eta(i) = etad(i)
        iflag(i) = 1                                                 
	x(i,1) = 0.d0
	y(i,1) = 0.d0
c                                                iflag = 1            
        return                                                       
c                                                                     
c       reversal point
c       --------------
c                                                                       
ccc21JAN16  110   if (iflag(i) .eq. 1) etad(i) = eta(i)
 110       iflag(i) = 2                                                  
c           write(*,*)'reversal'
	do 12 j = 2,20
	x(i,j) = 0.d0
  12    y(i,j) = 0.d0
c                                                                      
	if ( (x(i,1).eq.0.d0 .and. y(i,1).eq.0.d0) .or. ical.eq.1 ) then
           x(i,1)  = bx                                             
           y(i,1)  = by                                        
           rna(i)  = 0.d0                                               
        end if
c        write(*,*)'bx,by',bx,by
c        write(*,*)'x(1,1),y(1,1)',x(1,1),y(1,1)
c
c
        if (eta(i) .gt. etahs(i,1)) then
          etahs(i,1)  = eta(i)
          etahs(i,21) = etar(i)
        end if
c
        do 15 j=2,20
	etahs(i,j) = 0.d0
 15     etahs(i,j+20) = 0.d0
c                                                iflag = 2            
c
        return                                                        
c                                                                      
c
c       unloading, reloading, general loading ...
c       -----------------------------------------
c                                                                       
 20     call calhp1 (cx(i,k-1), cy(i,k-1), x(i,k-1), y(i,k-1),
     1               cxx, cyy, bx, by, r, iflag(i), i)       
c        write(*,*)'un/re/gen loading'
c        write(*,*)'k,ical',k,ical
c        write(*,*)'iflag(1)',iflag(1)
c                                                                       
        if ( ical.eq.1 .and. k.eq.2 ) go to 10
        if ( ical.eq.1 .and. k.gt.2 ) then
	   k = k - 1
           iflag(i) = k
	   go to 20
        end if
c
        ro = r                                                         
        rn = dsqrt( (bnx-cxx)**2 + (bny-cyy)**2 )                      
	rna(i) = r
c                                                                       
        if(rn.ge.ro) go to 210                                       
c                                                                      
c       reversal point (cxu,cyu) , (xu,yu)                            
c       ----------------------------------                              
c                                                                       
        iflag(i) = iflag(i) + 1
c
	do 22 j = k+1,20
	x(i,j) = 0.d0
        y(i,j) = 0.d0
	cx(i,j) = 0.d0
  22    cy(i,j) = 0.d0
c
c                                                                      
	if ( (x(i,k).eq.0.d0 .and. y(i,k).eq.0.d0) .or. ical.eq.1 ) then
           x(i,k)  = bx                                              
           y(i,k)  = by                                           
           cx(i,k) = cxx                                          
           cy(i,k) = cyy                                       
	end if
c                                                                      
        rna(i)  = 0.d0                                                
        if (iflag(i) .eq. 21)  rna(i) = dsqrt((x(i,k)-cx(i,k))**2 + 
     1                              (y(i,k)-cy(i,k))**2)
c
ccc     if ( eta(i) .gt. etahs(i,k,ielem) .and. iflag(i) .ne. 20 ) then
        if ( eta(i) .gt. etahs(i,k) ) then
	      etahs(i,k)    = eta(i)
          etahs(i,k+20) = etar(i)
        end if
c
        do 25 j=k+1,20
	      etahs(i,j) = 0.d0
 25       etahs(i,j+20) = 0.d0
c                                      
c                                                iflag = iflag + 1
        return                                                          
c                                                                     
c       previous loading continues                    
c       --------------------------
c                                                                       
  210   rgo = dsqrt((x(i,k-1)-cx(i,k-1))**2 + (y(i,k-1)-cy(i,k-1))**2) 
              rgn = dsqrt( (bx-cx(i,k-1))**2 + (by-cy(i,k-1))**2 )      
c        write(*,*)'previous'
ccccc   rgn = dsqrt( (bnx-cx(i,k-1))**2 + (bny-cy(i,k-1))**2 )      
c                                                                      
        if (rgn.ge.rgo) then
                 rna(i) = rgo
ccccccc   rna(i) = dsqrt( (bx-cx(i,k-1))**2 + (by-cy(i,k-1))**2 )      
ccccccc   rna(i) = rgn
           iflag(i) = iflag(i) - 1
	   if ( iflag(i) .eq. 20 ) eta(i) = 0.d0
	   k = k - 1
c
	do 32 j = k,20
	x(i,j) = 0.d0
        y(i,j) = 0.d0
	cx(i,j) = 0.d0
  32    cy(i,j) = 0.d0
c
ccccc   if ( k .eq. 1 ) return
	   if ( k .eq. 1 ) go to 10
	   go to 210
        end if
ccccif ( iflag(i) .eq. 21 ) rna(i) = rgn
c                                                                       
c                                                iflag = iflag
        return                                                          
        end
c                                                                     
c===============================================================       
       subroutine calhp1 (cxo,cyo,xo,yo,cxi,cyi,xi,yi,r,iflag,i)   
c                                                                    
c      ---------------------------------------------------------
c                                                                     
       implicit real*8  ( a-h, o-z )                           
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3, igaus
       common / cal   / ical
       dimension iflag(3)
c                                                                      
       ical = 0
       rl  = dsqrt ((xo-cxo)**2 + (yo-cyo)**2 )
       rs  = dsqrt ((xi-cxo)**2 + (yi-cyo)**2 )
c       write(*,*)'xo,yo',xo,yo
c       write(*,*)'xi,yi',xi,yi
c       write(*,*)'cxo,cyo',cxo,cyo
c       write(*,*)'rl',rl
c       write(*,*)'rs',rs

         if (rs.gt.rl .and. iflag(i).lt.21) then
           ical = 1
           return
           end if
c
c      if (dabs(xo-cxo).gt.1.d-5.and.dabs(xo-xi).gt.1.d-3) go to 100    
       if (dabs(xo-cxo).gt.1.d-20) go to 100    
c                                                                       
       cxi = cxo                                                        
c                                                                     
       if (dabs(yo - yi) .gt. 1.d-5) go to 101                          
       cyi = 0.50d0*(yo + yi)                                          
       go to 200                                                      
c                                                                     
  101  cyi = 0.5d0*(yo + yi) + 0.5d0*(cxo - xi)**2/(yi - yo)            
       go to 200                                                        
c                                                                       
 100   a1  = xo - cxo                                                 
       a2  = yo - cyo                                                  
       a3  = xo - xi                                                  
       a4  = yo - yi  
c                                                                       
       b   = a1*a3 + a2*a4                                              
       c   = 0.5d0*a1*(xo*xo+yo*yo-xi*xi-yi*yi)+a4*(cxo*yo-cyo*xo)      
c                                                                      
       if (b.eq.0.d0) then
       cxi = cxo
       cyi = cyo
       go to 200
       end if
c
       cxi = c/b                                                      
       cyi = a2/a1*(cxi - cxo) + cyo                                  
c                                                                       
 200   r   = dsqrt( (cxi - xi)**2 + (cyi - yi)**2 )                     
c                                                                   
       return                                                           
       end                                                            
c                                                                     
c================================================================
        subroutine cnjpnt ( x0, y0, dx0, dy0, cx, cy, r, xc, yc,
     &                      bet, betp, pi, idn) 
c
c                      ... conjugate point at the Failure Surface
c                                                                     
c       ---------------------------------------------------------
c                                                                      
        implicit real*8 ( a-h, o-z )                       
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3, igaus
c                                                                      
        if (dabs(dx0) .gt. 1.d-20 .and. dabs(dy0/dx0) .le. 1.d+05)
     1      go to 200                          
c   gen-ver           if (dabs(dx0) .gt. 1.d-10 ) go to 200
        if (dabs(dy0) .gt. 1.d-20) go to 100                            
c                                                                      
        write(6,*) '### / Stop in S-D model / ### cnjpnt:
     &  dx=dy=0,idn,ielem,istep',dx0,dy0,idn,ielem,istep 
        stop                                                           
c                                                                       
c       point at the Y-axis
c       -------------------
  100   d = r**2 - (x0 - cx)**2                                       
        if (d.ge.0.0d0) go to 10                                        
c                                                                     
        write(6,*) '### / Stop in S-D model / ### cnjpnt:
     &  no conjugate point(1); d,r,x0,cx,idn,ielem,istep',
     &                         d,r,x0,cx,idn,ielem,istep
        stop                                                           
c                                                                       
   10   xc = x0                                                      
        yc = cy + d**0.5d0
        if (dy0.lt.0.0d0) yc = cy - d**0.5d0
        go to 30                                                        
c                                                                      
c       other than the Y-axis
c       ---------------------
  200   s = dy0/dx0                                                   
        a = 1.0d0 + s**2                                              
        b = s*(y0 - cy) - cx - s*s*x0                                  
        c = (s*x0)**2 - 2.0d0*s*x0*(y0-cy) + (y0-cy)**2 - r*r + cx*cx 
        d = b*b - a*c                                                 
c                                                                      
        if (d.ge.0.0d0) go to 20                                   
c                                                                      
        write(6,*) '### / S-D model: Error 3 / ### cnjpnt:
     &  no conjugate point(2); dx0,dy0,idn,ielem,istep',
     &                         dx0,dy0,idn,ielem,istep
        stop                                                          
c                                                                     
   20   xc = (-1.d0*b + d**0.5d0)/a                                        
        yc = s*(xc - x0) + y0                                           
c                                                                       
        if((dx0*(xc - x0) + dy0*(yc - y0)).gt.0.0) go to 30           
c                                                                     
        xc = (-1.d0*b - d**0.5d0)/a                                        
        yc = s*(xc - x0) + y0                                          
c                                                                     
c       angle
c       -----
   30   call angle ( xc, yc, pi, bet, 1.d-20 )
        xcp = xc - cx                                                   
        ycp = yc - cy                                                   
        call angle ( xcp, ycp, pi, betp, 1.d-20 )
c                                                                
        return                                                        
        end                                                            
c
c=====================================================================
       subroutine psspar ( props, nmats, sig, ec, prefe, si,  p,
     1                     ggmin, ggmax, rf0, fis )
c
c                      ... plastic stress-strain parameters
c
c      ---------------------------------------------------------------
c
c
       implicit real*8 ( a-h, o-z )
       common / strn   /  eta(3), etacum(3), etarev(3), etar(3), etad(3)
       common / axil   /  svin, rf0in, pwpr, ainp
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3, igaus
       dimension props(nmats), sig(4)
c
c      reference lines
c      ---------------
c
       p1     = props(39)                                              
       p2     = 10.d0
       p3     = 30.d0
       p4     = 50.d0
       p5     = 100.d0
       p6     = 200.d0
       p7     = 400.d0
       p8     = 5000.d0
c
       eo     = props(40)                                              
       elim   = eo - 0.001d0
c
       esr    = 0.0d0
       alqs   = 0.0d0
       e1     = props(41)
       e2     = props(42)                                             
       e3     = props(43)                                              
       e4     = props(44)                                              
       e5     = props(45)                                             
       e6     = props(46)                                               
       e7     = props(47)                                             
       e8     = e7 - 23.d0 * (e6-e7)
c
c      plastic stress-strain coefficients
c      ----------------------------------
c
       d1     = props(32)                                               
       b1     = props(33)                                              
       d2     = props(34)                                               
       b2     = props(35)                                             
       d3     = props(36)                                               
       b3     = props(37)                                             
c
c      initial vertical stress and mean stress
c      ---------------------------------------
c
       if ( etarev(1).eq.0.d0 .and. eta(1).eq.0.d0 ) then
       svin = sig(2)
       end if
c
c
c      mean stress
c      -----------                                                      
c                                                                       
       p = (sig(1) + sig(2))*0.5d0                                      
       if (p.lt.1.d-11) then                                           
       write(6,*) '### / Error in S-D Model / ### 
     1             p is extension ',
     2            'sig(1),sig(2),sig(3)=',sig(1),sig(2),sig(3)
       write (6,*) 'elem step iter subinc = ',ielem,istep,iiter,incfai
cccxxxccc
        sig(1) = 2.d0
        sig(2) = 2.d0
        sig(3) = sig(3)
        p = (sig(1) + sig(2))*0.5d0        
	sig(4) = p
c  do 11 i = 4,6
c11         if ((dabs(sig(i))/p) .gt. rf0) 
c    1        sig(i) = ( sig(i)/dabs(sig(i)) )*(rf0-0.05d0*rf0)*sig(2)
ccccc          stop                                
            end if                                 
c
c
c      state index - is                                                
c      ----------------                                                
c                                                                      
       if (esr.ne.0.d0) then        
       ellls = dlog10(p/prefe)                                         
       es    = esr - alqs*ellls                                         
       end if                                                           
c                                                                       
       if (esr.ne.0.d0.and.p.lt.10.0d0) then            
       de10 = dlog10(10.0d0/prefe)                                      
       e10  = esr - alqs*de10                                           
       des  = (elim - e10)*(1.0d0 - p/10.0d0)                           
       es   = e10 + des                                                 
       end if                                                          
c                                           
       if (esr.eq.0.d0) then                
       if (p.le.p1) es = e1 + (p1-p)*(elim - e1)/p1                     
       if (p.gt.p1.and.p.le.p2) es = e2 + (p2 - p)*(e1 - e2)/(p2 - p1) 
       if (p.gt.p2.and.p.le.p3) es = e3 + (p3 - p)*(e2 - e3)/(p3 - p2) 
       if (p.gt.p3.and.p.le.p4) es = e4 + (p4 - p)*(e3 - e4)/(p4 - p3)  
       if (p.gt.p4.and.p.le.p5) es = e5 + (p5 - p)*(e4 - e5)/(p5 - p4) 
       if (p.gt.p5.and.p.le.p6) es = e6 + (p6 - p)*(e5 - e6)/(p6 - p5) 
       if (p.gt.p6.and.p.le.p7) es = e7 + (p7 - p)*(e6 - e7)/(p7 - p6)  
       if (p.gt.p7.and.p.le.p8) es = e8 + (p8 - p)*(e7 - e8)/(p8 - p7) 
       if (p.gt.p8) write (6,*) '### / Stop in S-D model / ### 
     &    p is out of range;  p, p8, ielem, istep = ',p,p8,ielem,istep
       if (p.gt.p8) stop                                                
       end if                                                           
c                                                                      
       si  = (eo - ec)/(eo - es)                                        
c       write(*,*)'p',p
c       write(*,*)'eo,ec,es',eo,ec,es
c       write(*,*)'si',si
       if ( fis .eq. 0.d0 ) fis = si
c                                                                      
c      stress-strain parameters
c      ------------------------                               
c
       rf0   = d1*si + b1                                             
       ggmin = d2*si + b2                                               
       ggmax = d3*si + b3                                               
c                                                                      
       return
       end
c
c==================================================
       subroutine angle ( x, y, pi, angl, cr )
c
c                   ... anti-clockwise direction
c                       with respect to X-axis
c
c      --------------------------------------------
c
c
       implicit real*8 ( a-h, o-z )
c
       if ( y .eq. 0.d0) then
	  angl = 0.d0
	  return
       end if
c
       if ( dabs(x).lt.cr .or. ( dabs(x).ne.0.d0 .and. 
     1        dabs(y/x).gt.1000.d0 ) ) then
        if ( dabs(y) .lt. cr ) write (6,100) x, y, cr
          angl = ( y/dabs(y) ) * pi/2.0d0
          go to 99
        end if
c
       angl = datan (y/x)
       if ( x.lt.0.d0 )     angl = angl + pi
   99  if ( angl.lt.0.d0 )  angl = angl + 2.0d0*pi
c
 100   format ('x = ',e12.5,'y = ',e12.5,'cr = ',e12.5/
     1    'Warning in ## angle ##')
c
       return
       end
c
c============================================================
       subroutine esspar ( props, nmats, ec, prefe, p, fis, oths )
c
c                     ... elastic parameters and coefficients
c
c      ------------------------------------------------------
c
       implicit real*8 (a-h, o-z)
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3, igaus
       common / elpar /  ae, be, xgi, xg2, eyng, pora, xn, coef2
       common / axil  /  svin, rf0in, pwpr, ainp
        common / strn  / eta(3), etacum(3), etarev(3), etar(3), etad(3)
       common / tstr  / fmuf, depsx, depsy
       dimension props(nmats), oths(10)

c
       iitr = oths(11)
c       write(*,*)'iitr',iitr
       if ( fmuf.gt.0.15d0 .and. etarev(1).gt..02d0 .and. 
     1           etacum(1).gt.0.d0 ) then
c      write (14,9821) istep,ielem,iiter,
c    1                      oths(8,ielem),etarev(1),etacum(1)
 9821   format ('befor',3i5,5x,f8.2,2(e12.5,2x))
           afc  =  (1.d0 / fis)
           if ( afc .gt. 0.5d0 ) afc = 0.5d0
           if (etacum(1) .eq. 1.d0 ) afc = afc/fis
           if (iitr.ne.1) then
               props(6)=oths(8)-afc*etacum(1)*oths(8)
           end if 
           etacum(1) = 0.d0
       end if
c
       gelr  = props( 6) * 98.1d0
       pora  = props( 4)
       xn    = props(28)
c
ccccccccccc    etad(i)   is    redirrected in use  CAUTION !!!
          if ( fmuf.gt.0.15d0 .and. etarev(1) .gt. 0.02d0 ) then
            fct   = ( etarev(1) / 0.05d0 )
            if ( fct.gt.1.d0 ) fct = 1.d0
            xn    = props(28) + ( 0.85d0 - props(28) )*fct
c           etarev(1) = 0.d0
          end if
c
c      elastic moduli and coefficients
c      -------------------------------
c
       xgi  = gelr*((2.17d0-ec)**2.0d0 / (1.0d0+ec))*(p/prefe)**xn
       eyng = 2.d0*xgi*(1.d0 + pora)
c                                                                      
c      plane stress                                                     
c      ------------
c
c      ae  = eyng / (1.d0 - pora**2)
c      be  = (eyng*pora) / (1.d0 - pora**2)
c      xg2 = xgi*2.d0                                                  
c                                                                      
c      plane strain
c      ------------
       ae  = (eyng*(1.d0 - pora)) / ((1.d0+pora)*(1.d0-2.d0*pora))
       be  = (eyng*pora) / ((1.d0+pora)*(1.d0-2.d0*pora))
       xg2 = xgi*2.d0        

c       write(*,*)'gelr',gelr
c       write(*,*)'ec',ec
c       write(*,*)'p',p
c       write(*,*)'prefe',prefe
c       write(*,*)'xn',xn
c       write(*,*)'eyng',eyng
c       write(*,*)'pora',pora
c       write(*,*)'xgi',xgi
c       write(*,*)'ae,be,xg2',ae,be,xg2
c
       return
       end
c
c==============================================================
       subroutine strinc ( depsx, depsy, incrmt, ddeps, deps )
c
c                                ... division to subincrements
c
c      ---------------------------------------------------------
c
c
       implicit real*8 ( a-h, o-z )
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3, igaus
       dimension ddeps(4), deps(4)
c
c                                                                       
c      subincrements
c      -------------
c
       delinc = dsqrt ((depsx)**2 + (depsy)**2)
c       write(*,*)'input deviatoric strain = ',delinc
       dincr  = delinc / 0.0001d0
       incrmt = idint(dincr)
c
       if (incrmt.lt.1) incrmt = 1
       if (incrmt.gt.100) 
     1 write (6,*) '*** / Warning in S-D model / * strinc * ',
     2 '  incrmt=',incrmt,
     3 '  ielem=',ielem,
     4 '  istep=',istep
       if (incrmt.gt.1000) then
       write (6,*) '### / Stop in S-D model / # strinc # incrmt=',incrmt
       stop
       end if
       incfai = 0                                                      
c                                                                      
       do 11 i = 1,3
 11    ddeps(i) = deps(i) / incrmt  
c       write(*,*)'dincr',dincr
c       write(*,*)'incrmt',incrmt
c       write(*,*)'ddeps',ddeps
c                                                                       
       return
       end
c
c=====================================================
       subroutine strcon ( sig, c1, c2, rf0, pk, idn , ialarm , iala)
c
c                      ... stress state control
c
c      -----------------------------------------------
c
       implicit real*8 ( a - h , o - z )
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3, igaus
       dimension sig(4)
c
       ialarm = 0
 123   pk = 0.5d0 * ( sig(1) + sig(2) )    
       c11  = c1*pk                                                   
       c22  = c2*pk                                                   
       qk   = dsqrt ( ( (sig(2)-sig(1))/2.d0 - c11)**2
     1              + (sig(3)-c22)**2 )
c      qz   = dsqrt ( ( (sig(3)-sig(1))/2.d0 - c11)**2
c    1              + (sig(6)-c22)**2 )
c      qy   = dsqrt ( ( (sig(2)-sig(3))/2.d0 - c11)**2
c    1              + (sig(5)-c22)**2 )
c
       rk   = dabs( qk/pk )                                           
c      rz   = dabs( qz/pk )                                           
c      ry   = dabs( qy/pk )                                           
       rl   = rf0 * (1.d0 - 1.d-10) 
c                                                                     
       if ( sig(1).lt.1.d-11 .or. sig(2).lt.1.d-11 ) then
c         write (14,*)' negative sigma', istep, iiter, ielem, idn
	  ialarm = 1
c         write (14,*) 'before correction',sig(1), sig(2), sig(3)
c        if ( sig(1).gt.1.d-11 .and. sig(2).lt.1.d-11 ) then
c         write (14,*)' negative sigma 1', istep, iiter, ielem, idn
c        elseif ( sig(1).lt.1.d-11 .and. sig(2).gt.1.d-11 ) then
c         write (14,*)' negative sigma 2', istep, iiter, ielem, idn
c        else
c         write (14,*)' negative sigma 3', istep, iiter, ielem, idn
c        end if
       end if
 9699    format (4(5x,f8.5))
 9799    format (6(2x,e11.4))
 9788    format (3(2x,i5))
c
        qnk  = rf0*(1.d0-0.0005d0)*pk
c
       if ( rk .gt. rl ) then
c
c      ialarm = 1
cccc          qnk  = rf0*(1.d0-0.0005d0)*pk
cccc        qnk  = rf0*(1.d0-0.05d0)*pk
       qrat = qnk/qk                                                 
c                                                                     
       stdf   = (sig(2) - sig(1)) / 2.0d0 
       sig(1) = dabs ( pk - qrat*stdf )  
       sig(2) = dabs ( pk + qrat*stdf )                                
       sig(3) = qrat * sig(3)        
       end if                                                           
c
c      if ( rz.gt.rl .or. ry.gt.rl) then
c   sig(6) = ( qnk/qz ) * sig(6)
c        write (11,*) 'x-z &/or y-z need correction'
c        write (11,9788) ielem, istep, iiter
c        write (11,9799) (sig(k1),k1=1,6), rk
c        write (11,9799) qnk,qrat,stdf, rl, rf0
c      end if
c
c
       return                                                           
       end                                                              
c
c======================================================================
       subroutine fxinc ( ddeps, p, pi, bigxp, bigyp, c1, c2, rf0, rna,
     1             dbigx, dbigy, dbigxp, dbigyp, eangp, pangp )
c
c                     ... first assumption for the stress increment
c
c      ----------------------------------------------------------------
c
c
       implicit real*8 (a-h, o-z)
       common / elpar / ae, be, xgi, xg2, eyng, pora, xn, coef2
       common / elmnt / ielem, istep,iiter,incfai,sig1,sig2,sig3, igaus
       dimension  ddeps(4), rna(3)
c
c                                                                     
c      direction of total strain increment - tang                    
c      -----------------------------------                              
c                                                                       
       ddepsx = ddeps(2) - ddeps(1)                                     
       ddepsy = ddeps(3)*2.0d0                                          
       btst = dsqrt( ddepsx**2 + ddepsy**2 )   
       call angle ( ddepsx, ddepsy, pi, tang, 1.d-20 )
c                                                                       
c      boundary conjugate points
c      -------------------------
c                                                                       
       call cnjpnt ( bigxp, bigyp, ddepsx, ddepsy, c1, c2, rf0,
     1               xecp, yecp, eang, eangp, pi, 1 )
       call cnjpnt ( 0.d0, 0.d0, ddepsx, ddepsy, 0.d0, 0.d0, rf0,
     1               xpcp, ypcp, pang, pangp, pi, 2 )
c                                                                    
c                                                                       
c      calculation of stress increments
c      --------------------------------                                 
c                                                                      
       fact = 0.2d0
ccc       dbigx  =  fact*(be - ae)*ddepsx
ccc       dbigy  =  fact*xg2*ddepsy                             
          dbigx  =  fact*(ae - be)*ddepsx/2.d0
          dbigy  =  fact*xg2*ddeps(3)
c
       dbigxp =  dbigx / p               
       dbigyp =  dbigy / p                    
c                                                                     
       return
       end
c
