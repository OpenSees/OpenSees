
      INTEGER FUNCTION BISRCH (N,K,L) 
C       
C ... BISRCH IS AN INTEGER FUNCTION WHICH USES A BISECTION SEARCH   
C     TO FIND THE ENTRY J IN THE ARRAY K SUCH THAT THE VALUE L IS   
C     GREATER THAN OR EQUAL TO K(J) AND STRICTLY LESS THAN K(J+1).  
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTOR K    
C          K      INTEGER VECTOR      
C          L      INTEGER CONSTANT SUCH THAT  K(J) .GE. L .LT. K(J+1) 
C                 WITH J RETURNED AS VALUE OF INTEGER FUNCTION BISRCH 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,L,K(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER JLEFT,JMID,JRIGHT       
C       
      JLEFT = 1   
      JRIGHT = N  
      IF (N.EQ.2) GO TO 40  
      JMID = (N+1)/2
C       
   10 IF (L.GE.K(JMID)) GO TO 20      
C       
C ...... L .GE. K(LEFT)  AND  L .LT. K(JMID)    
C       
      JRIGHT = JMID 
      GO TO 30    
C       
C ...... L .GE. K(JMID)  AND  L .LT. K(JRIGHT)  
C       
   20 JLEFT = JMID
C       
C ...... TEST FOR CONVERGENCE 
C       
   30 IF (JRIGHT-JLEFT.EQ.1) GO TO 40 
      JMID = JLEFT+(JRIGHT-JLEFT+1)/2 
      GO TO 10    
C       
C ...... BISECTION SEARCH FINISHED    
C       
   40 BISRCH = JLEFT
C       
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION CHEBY (QA,QT,RRR,IP,CME,SME)
C       
C     COMPUTES THE SOLUTION TO THE CHEBYSHEV EQUATION     
C       
C ... PARAMETER LIST:       
C       
C          QA     RATIO OF PSEUDO-RESIDUALS     
C          QT     VIRTUAL SPECTRAL RADIUS       
C          RRR    ADAPTIVE PARAMETER  
C          IP     NUMBER OF ITERATIONS SINCE LAST CHANGE OF 
C                     PARAMETERS      
C          CME,   ESTIMATES FOR THE LARGEST AND SMALLEST EIGEN-     
C          SME      VALUES OF THE ITERATION MATRIX
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IP  
      DOUBLE PRECISION CME,QA,QT,RRR,SME
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      DOUBLE PRECISION X,Y,Z
C       
      Z = .5D0*(QA+DSQRT(DABS(QA**2-QT**2)))*(1.D0+RRR**IP) 
      X = Z**(1.D0/DBLE(FLOAT(IP)))   
      Y = (X+RRR/X)/(1.D0+RRR)
C       
      CHEBY = .5D0*(CME+SME+Y*(2.D0-CME-SME))   
C       
      RETURN      
      END 
      SUBROUTINE CHGCON (TRI,GAMOLD,RHOOLD,IBMTH) 
C       
C     COMPUTES THE NEW ESTIMATE FOR THE LARGEST EIGENVALUE FOR      
C     CONJUGATE GRADIENT ACCELERATION.
C       
C ... PARAMETER LIST:       
C       
C          TRI    TRIDIAGONAL MATRIX ASSOCIATED WITH THE EIGENVALUES
C                    OF THE CONJUGATE GRADIENT POLYNOMIAL 
C          GAMOLD 
C            AND  
C          RHOOLD PREVIOUS VALUES OF ACCELERATION PARAMETERS
C          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY CG 
C                      IBMTH = 1,  JACOBI       
C                            = 2,  REDUCED SYSTEM 
C                            = 3,  SSOR 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IBMTH 
      DOUBLE PRECISION TRI(2,1),GAMOLD,RHOOLD   
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IB2,IB3,IER,IP
      DOUBLE PRECISION CMOLD,END,START,EIGVSS,EIGVNS      
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      GO TO (10,20,30), IBMTH 
C       
C ... JACOBI CONJUGATE GRADIENT       
C       
   10 START = CME 
      IP = IN     
      GO TO 40    
C       
C ... REDUCED SYSTEM CG     
C       
   20 START = CME**2
      IP = IN     
      GO TO 40    
C       
C ... SSOR CG     
C       
   30 IF (ADAPT) START = SPR
      IF (.NOT.ADAPT) START = SPECR   
      IP = IN-IS  
C       
C ... DEFINE THE MATRIX     
C       
   40 IF (IP.GE.2) GO TO 60 
      IF (IP.EQ.1) GO TO 50 
C       
C ... IP = 0      
C       
      END = 0.D0  
      CMOLD = 0.D0
      GO TO 110   
C       
C ... IP = 1      
C       
   50 END = 1.D0-1.D0/GAMMA 
      TRI(1,1) = END
      TRI(2,1) = 0.D0       
      GO TO 110   
C       
C ... IP > 1      
C       
   60 IF ((IP.GT.2).AND.(DABS(START-CMOLD).LE.ZETA*START)) GO TO 120
      CMOLD = START 
C       
C ... COMPUTE THE LARGEST EIGENVALUE  
C       
      TRI(1,IP) = 1.D0-1.D0/GAMMA     
      TRI(2,IP) = (RHO-1.D0)/(RHO*RHOOLD*GAMMA*GAMOLD)    
      IF (ISYM.NE.0) GO TO 80 
      END = EIGVSS(IP,TRI,START,ZETA,ITMAX,IER) 
      IF (IER.EQ.0) GO TO 100 
      IF (LEVEL.GE.2) WRITE (NOUT,70) IER       
   70 FORMAT (/10X,'DIFFICULTY IN COMPUTATION OF MAXIMUM EIGENVALUE'/15X
     *   ,'OF ITERATION MATRIX'/10X,'SUBROUTINE ZBRENT RETURNED IER =', 
     *   I5)      
      GO TO 100   
   80 IB2 = 1+IP  
      IB3 = IB2+IP/2+1      
      END = EIGVNS(IP,TRI,TRI(1,IB2),TRI(1,IB3),IER)      
      IF (IER.EQ.0) GO TO 100 
      IF (LEVEL.GE.2) WRITE (NOUT,90) IER       
   90 FORMAT (/10X,'DIFFICULTY IN COMPUTATION OF MAXIMUM EIGENVALUE'/15X
     *   ,'OF ITERATION MATRIX'/10X,'SUBROUTINE EQRT1S RETURNED IER =', 
     *   I5)      
  100 CONTINUE    
      IF (IER.NE.0) GO TO 130 
C       
C ... SET SPECTRAL RADIUS FOR THE VARIOUS METHODS 
C       
  110 IF (IBMTH.EQ.1) CME = END       
      IF (IBMTH.EQ.2) CME = DSQRT(DABS(END))    
      IF (IBMTH.EQ.3.AND.ADAPT) SPR = END       
      IF (IBMTH.EQ.3.AND..NOT.ADAPT) SPECR = END
      RETURN      
C       
C ... RELATIVE CHANGE IN CME IS LESS THAN ZETA.  THEREFORE STOP     
C     CHANGING.   
C       
  120 ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      RETURN      
C       
C ... ESTIMATE FOR CME > 1.D0.  THEREFORE NEED TO STOP ADAPTIVE     
C     PROCEDURE AND KEEP OLD VALUE OF CME.      
C       
  130 ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      IF (LEVEL.GE.2) WRITE (NOUT,140) IN,START 
  140 FORMAT (/10X,'ESTIMATE OF MAXIMUM EIGENVALUE OF JACOBI   '/15X, 
     *   'MATRIX (CME) NOT ACCURATE'/10X,       
     *   'ADAPTIVE PROCEDURE TURNED OFF AT ITERATION ',I5/10X,      
     *   'FINAL ESTIMATE OF MAXIMUM EIGENVALUE =',D15.7/) 
C       
      RETURN      
      END 
      SUBROUTINE CHGSI (DTNRM,IBMTH)  
C       
C ... COMPUTES NEW CHEBYSHEV ACCELERATION PARAMETERS ADAPTIVELY.    
C       
C ... PARAMETER LIST:       
C       
C          DTNRM  NUMERATOR OF RAYLEIGH QUOTIENT
C          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY SI 
C                      IBMTH = 1,   JACOBI      
C                            = 2,   REDUCED SYSTEM
C                            = 3,   SYMMETRIC SOR 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IBMTH 
      DOUBLE PRECISION DTNRM
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      DOUBLE PRECISION CMOLD,ZM1,ZM2  
C       
C ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
C       
      DOUBLE PRECISION CHEBY
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      GO TO (10,30,50), IBMTH 
C       
C     --------------------- 
C ... JACOBI SEMI-ITERATIVE 
C     --------------------- 
C       
C ... CHEBYSHEV EQUATION    
C       
   10 CONTINUE    
      IF (IN.EQ.0) ZM1 = CME
      IF (IN.NE.0) ZM1 = CHEBY(QA,QT,RRR,IN-IS,CME,SME)   
C       
C ... RAYLEIGH QUOTIENT     
C       
      ZM2 = DTNRM/DELNNM    
C       
C ... COMPUTATION OF ITERATIVE PARAMETERS       
C       
      CMOLD = CME 
      CME = DMAX1(ZM1,ZM2,CMOLD)      
      IF (CME.GE.1.D0) GO TO 20       
      IF (CASEII) SME = -CME
      SIGE = (CME-SME)/(2.D0-CME-SME) 
      GAMMA = 2.D0/(2.D0-CME-SME)     
      RRR = (1.D0-DSQRT(DABS(1.D0-SIGE*SIGE)))/(1.D0+DSQRT(DABS(1.D0- 
     *   SIGE*SIGE)))       
      IS = IN     
      DELSNM = DELNNM       
      RHO = 1.D0  
      IF (LEVEL.GE.2) WRITE (NOUT,90) IN,ZM1,ZM2,CME,GAMMA,CME      
      RETURN      
C       
C ... ADAPTIVE PROCEDURE FAILED FOR JACOBI SI   
C       
   20 CME = CMOLD 
      ADAPT = .FALSE.       
      IF (LEVEL.GE.2) WRITE (NOUT,110) IN,CME   
      RETURN      
C       
C     -----------------------------   
C ... REDUCED SYSTEM SEMI-ITERATIVE   
C     -----------------------------   
C       
C ... CHEBYSHEV EQUATION    
C       
   30 CONTINUE    
      IF (IN.EQ.0) ZM1 = CME
      IF (IN.NE.0) ZM1 = CHEBY(QA,QT,RRR,2*(IN-IS),0.D0,0.D0)       
C       
C ... RAYLEIGH QUOTIENT     
C       
      ZM2 = DSQRT(DABS(DTNRM/DELNNM)) 
C       
C ... COMPUTATION OF NEW ITERATIVE PARAMETERS   
C       
      CMOLD = CME 
      CME = DMAX1(ZM1,ZM2,CMOLD)      
      IF (CME.GE.1.D0) GO TO 40       
      SIGE = CME*CME/(2.D0-CME*CME)   
      GAMMA = 2.D0/(2.D0-CME*CME)     
      RRR = (1.D0-DSQRT(DABS(1.D0-CME*CME)))/(1.D0+DSQRT(DABS(1.D0-CME* 
     *   CME)))   
      IS = IN     
      DELSNM = DELNNM       
      RHO = 1.D0  
      IF (LEVEL.GE.2) WRITE (NOUT,90) IN,ZM1,ZM2,CME,GAMMA,CME      
      RETURN      
C       
C ... ADAPTIVE PROCEDURE FAILED FOR REDUCED SYSTEM SI     
C       
   40 CME = CMOLD 
      ADAPT = .FALSE.       
      IF (LEVEL.GE.2) WRITE (NOUT,110) IN,CME   
      RETURN      
C       
C     -----------------------------   
C ... SYMMETRIC SOR SEMI-ITERATIVE    
C     ----------------------------    
C       
   50 CONTINUE    
      IF (SPECR.EQ.0.D0) SPECR = .171572875D0   
      IF (IN.EQ.0) GO TO 60 
      ZM1 = CHEBY(QA,QT,RRR,IN-IS,SPECR,0.D0)   
      GO TO 70    
   60 ZM1 = SPECR 
      SPR = SPECR 
C       
C ... RAYLEIGH QUOTIENT     
C       
   70 ZM2 = DTNRM/DELNNM    
C       
C ... COMPUTATION OF NEW ESTIMATE FOR SPECTRAL RADIUS     
C       
      IF (ADAPT) GO TO 80   
C       
C ... PARTIALLY ADAPTIVE SSOR SI      
C       
      SPECR = DMAX1(ZM1,ZM2,SPECR)    
      IS = IN+1   
      DELSNM = DELNNM       
      IF (LEVEL.GE.2) WRITE (NOUT,100) IN,ZM1,ZM2,CME,SPECR 
      RETURN      
C       
C ... FULLY ADAPTIVE SSOR SI
C       
   80 SPR = DMAX1(ZM1,ZM2,SPR)
      RETURN      
C       
C ... FORMAT STATEMENTS     
C       
   90 FORMAT (/30X,'PARAMETERS WERE CHANGED AT ITERATION NO.',I5/35X, 
     *   'SOLUTION TO CHEBYSHEV EQN.       =',D15.7/35X,  
     *   'SOLUTION TO RAYLEIGH QUOTIENT    =',D15.7/35X,  
     *   'NEW ESTIMATE FOR CME             =',D15.7/35X,  
     *   'NEW ESTIMATE FOR GAMMA           =',D15.7/35X,  
     *   'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)     
C       
  100 FORMAT (/30X,'PARAMETERS WERE CHANGED AT ITERATION NO.',I5/35X, 
     *   'SOLUTION TO CHEBYSHEV EQN.       =',D15.7/35X,  
     *   'SOLUTION TO RAYLEIGH QUOTIENT    =',D15.7/35X,  
     *   'NEW ESTIMATE FOR CME             =',D15.7/35X,  
     *   'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)     
C       
  110 FORMAT (/10X,'ESTIMATE OF MAXIMUM EIGENVALUE OF JACOBI   '/15X, 
     *   'MATRIX (CME) TOO LARGE'/10X,
     *   'ADAPTIVE PROCEDURE TURNED OFF AT ITERATION ',I5/10X,      
     *   'FINAL ESTIMATE OF MAXIMUM EIGENVALUE =',D15.7/) 
C       
      END 
      LOGICAL FUNCTION CHGSME (OLDNRM,ICNT)     
C       
C ... THIS FUNCTION TESTS FOR JACOBI SI WHETHER SME SHOULD BE CHANGED 
C ... WHEN CASEII = .FALSE..  IF THE TEST IS POSITIVE THE NEW VALUE 
C ... OF SME IS COMPUTED.   
C       
C ... PARAMETER LIST:       
C       
C          OLDNRM SQUARE OF THE NORM OF THE PSEUDO-RESIDUAL 
C                    AT THE LAST ITERATION      
C          ICNT   NUMBER OF ITERATIONS SINCE LAST CHANGE OF 
C                    PARAMETERS       
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER ICNT
      DOUBLE PRECISION OLDNRM 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IP  
      DOUBLE PRECISION Q,RN,SM1,SM2,WP,Z
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      CHGSME = .FALSE.      
      RN = DSQRT(DELNNM/OLDNRM)       
      IF (.NOT.(QA.GT.1.D0.AND.RN.GT.1.D0)) RETURN
      IF (IN.LE.IS+2) RETURN
C       
      ICNT = ICNT+1 
      IF (ICNT.LT.3) RETURN 
C       
C ... CHANGE SME IN J-SI ADAPTIVE PROCEDURE     
C       
      CHGSME = .TRUE.       
      SM1 = 0.D0  
      SM2 = 0.D0  
      IF (SME.GE.CME) GO TO 10
C       
C ... COMPUTE SM1 
C       
      IP = IN-IS  
      Q = QA*(1.D0+RRR**IP)/(2.D0*DSQRT(RRR**IP)) 
      Z = (Q+DSQRT(Q**2-1.D0))**(1.D0/DBLE(FLOAT(IP)))    
      WP = (Z**2+1.D0)/(2.D0*Z)       
      SM1 = .5D0*(CME+SME-WP*(CME-SME)) 
C       
C ... COMPUTE SM2 
C       
      Q = RN*(1.D0+RRR**IP)/((1.D0+RRR**(IP-1))*DSQRT(RRR)) 
      WP = (Q**2+1.D0)/(2.D0*Q)       
      SM2 = .5D0*(CME+SME-WP*(CME-SME)) 
C       
   10 SME = DMIN1(1.25D0*SM1,1.25D0*SM2,SME,-1.D0)
      SIGE = (CME-SME)/(2.D0-CME-SME) 
      GAMMA = 2.D0/(2.D0-CME-SME)     
      RRR = (1.D0-DSQRT(1.D0-SIGE**2))/(1.D0+DSQRT(1.D0-SIGE**2))   
      IS = IN     
      DELSNM = DELNNM       
      RHO = 1.D0  
C       
      IF (LEVEL.GE.2) WRITE (NOUT,20) IN,SM1,SM2,SME      
C       
   20 FORMAT (/30X,'ESTIMATE OF SMALLEST EIGENVALUE OF JACOBI'/37X, 
     *   'MATRIX (SME) CHANGED AT ITERATION ',I5/35X,     
     *   'FIRST ESTIMATE OF SME            =',D15.7/35X,  
     *   'SECOND ESTIMATE OF SME           =',D15.7/35X,  
     *   'NEW ESTIMATE OF SME              =',D15.7/)     
C       
      RETURN      
      END 
      SUBROUTINE ItDAXPY (N,DA,DX,INCX,DY,INCY)   
C       
C     OVERWRITE DOUBLE PRECISION DY WITH DOUBLE PRECISION DA*DX + DY. 
C       
      DOUBLE PRECISION DX(1),DY(1),DA 
      IF (N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF (INCX.EQ.INCY) IF (INCX-1) 10 , 30 , 70
   10 CONTINUE    
C       
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.     
C       
      IX = 1      
      IY = 1      
      IF (INCX.LT.0) IX = (-N+1)*INCX+1 
      IF (INCY.LT.0) IY = (-N+1)*INCY+1 
      DO 20 I = 1,N 
         DY(IY) = DY(IY)+DA*DX(IX)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 CONTINUE    
      RETURN      
C       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1    
C       
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4. 
C       
   30 M = N-(N/4)*4 
      IF (M.EQ.0) GO TO 50  
      DO 40 I = 1,M 
         DY(I) = DY(I)+DA*DX(I)       
   40 CONTINUE    
      IF (N.LT.4) RETURN    
   50 MP1 = M+1   
      DO 60 I = MP1,N,4     
         DY(I) = DY(I)+DA*DX(I)       
         DY(I+1) = DY(I+1)+DA*DX(I+1) 
         DY(I+2) = DY(I+2)+DA*DX(I+2) 
         DY(I+3) = DY(I+3)+DA*DX(I+3) 
   60 CONTINUE    
      RETURN      
C       
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.    
C       
   70 CONTINUE    
      NS = N*INCX 
      DO 80 I = 1,NS,INCX   
         DY(I) = DA*DX(I)+DY(I)       
   80 CONTINUE    
      RETURN      
      END 
      SUBROUTINE ItDCOPY (N,DX,INCX,DY,INCY)      
C       
C     COPY DOUBLE PRECISION DX TO DOUBLE PRECISION DY.    
C       
      DOUBLE PRECISION DX(1),DY(1)    
      IF (N.LE.0) RETURN    
      IF (INCX.EQ.INCY) IF (INCX-1) 10 , 30 , 70
   10 CONTINUE    
C       
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.      
C       
      IX = 1      
      IY = 1      
      IF (INCX.LT.0) IX = (-N+1)*INCX+1 
      IF (INCY.LT.0) IY = (-N+1)*INCY+1 
      DO 20 I = 1,N 
         DY(IY) = DX(IX)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 CONTINUE    
      RETURN      
C       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1    
C       
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7. 
C       
   30 M = N-(N/7)*7 
      IF (M.EQ.0) GO TO 50  
      DO 40 I = 1,M 
         DY(I) = DX(I)      
   40 CONTINUE    
      IF (N.LT.7) RETURN    
   50 MP1 = M+1   
      DO 60 I = MP1,N,7     
         DY(I) = DX(I)      
         DY(I+1) = DX(I+1)  
         DY(I+2) = DX(I+2)  
         DY(I+3) = DX(I+3)  
         DY(I+4) = DX(I+4)  
         DY(I+5) = DX(I+5)  
         DY(I+6) = DX(I+6)  
   60 CONTINUE    
      RETURN      
C       
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.    
C       
   70 CONTINUE    
      NS = N*INCX 
      DO 80 I = 1,NS,INCX   
         DY(I) = DX(I)      
   80 CONTINUE    
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION ItDDOT (N,DX,INCX,DY,INCY)  
C       
C     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY.
C       
      DOUBLE PRECISION DX(1),DY(1)    
      DDOT = 0.D0 
      IF (N.LE.0) RETURN    
      IF (INCX.EQ.INCY) IF (INCX-1) 10 , 30 , 70
   10 CONTINUE    
C       
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.     
C       
      IX = 1      
      IY = 1      
      IF (INCX.LT.0) IX = (-N+1)*INCX+1 
      IF (INCY.LT.0) IY = (-N+1)*INCY+1 
      DO 20 I = 1,N 
         DDOT = DDOT+DX(IX)*DY(IY)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 CONTINUE    
      RETURN      
C       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.   
C       
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5. 
C       
   30 M = N-(N/5)*5 
      IF (M.EQ.0) GO TO 50  
      DO 40 I = 1,M 
         DDOT = DDOT+DX(I)*DY(I)      
   40 CONTINUE    
      IF (N.LT.5) RETURN    
   50 MP1 = M+1   
      DO 60 I = MP1,N,5     
         DDOT = DDOT+DX(I)*DY(I)+DX(I+1)*DY(I+1)+DX(I+2)*DY(I+2)+DX(I+3)
     *      *DY(I+3)+DX(I+4)*DY(I+4)  
   60 CONTINUE    
      RETURN      
C       
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.       
C       
   70 CONTINUE    
      NS = N*INCX 
      DO 80 I = 1,NS,INCX   
         DDOT = DDOT+DX(I)*DY(I)      
   80 CONTINUE    
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION DETERM (N,TRI,XLMDA)      
C       
C     THIS SUBROUTINE COMPUTES THE DETERMINANT OF A SYMMETRIC       
C     TRIDIAGONAL MATRIX GIVEN BY TRI. DET(TRI - XLMDA*I) = 0       
C       
C ... PARAMETER LIST
C       
C          N      ORDER OF TRIDIAGONAL SYSTEM   
C          TRI    SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N 
C          XLMDA  ARGUMENT FOR CHARACTERISTIC EQUATION    
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N   
      DOUBLE PRECISION TRI(2,1),XLMDA 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER ICNT,L,NM1    
      DOUBLE PRECISION D1,D2,D3       
C       
      NM1 = N-1   
      D2 = TRI(1,N)-XLMDA   
      D1 = D2*(TRI(1,NM1)-XLMDA)-TRI(2,N)       
      IF (N.EQ.2) GO TO 20  
C       
C ... BEGINNING OF LOOP     
C       
      DO 10 ICNT = 2,NM1    
         L = NM1-ICNT+2     
         D3 = D2  
         D2 = D1  
         D1 = (TRI(1,L-1)-XLMDA)*D2-D3*TRI(2,L) 
   10 CONTINUE    
C       
C ... DETERMINANT COMPUTED  
C       
   20 DETERM = D1 
C       
      RETURN      
      END 
      SUBROUTINE DFAULT (IPARM,RPARM) 
C       
C ... THIS SUBROUTINE SETS THE DEFAULT VALUES OF IPARM AND RPARM.   
C       
C ... PARAMETER LIST:       
C       
C          IPARM  
C           AND   
C          RPARM  ARRAYS SPECIFYING OPTIONS AND TOLERANCES
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IPARM(12)     
      DOUBLE PRECISION RPARM(12)      
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
C     DRELPR  - COMPUTER PRECISION (APPROX.)    
C     IF INSTALLER OF PACKAGE DOES NOT KNOW DRELPR VALUE, 
C     AN APPROXIMATE VALUE CAN BE DETERMINED FROM A SIMPLE
C     FORTRAN PROGRAM SUCH AS 
C       
C     DOUBLE PRECISION DRELPR, TEMP   
C     DRELPR = 1.0D0
C   2 DRELPR = 0.5D0*DRELPR 
C     TEMP = DRELPR + 1.0D0 
C     IF(TEMP .GT. 1.0D0)  GO TO 2    
C     WRITE(6,3) DRELPR     
C   3 FORMAT(5X,D15.8)      
C     STOP
C     END 
C       
C     SOME VALUES ARE:      
C       
C     DRELPR = 1.26D-29  FOR CDC CYBER 170/750  (APPROX.) 2**-96    
C            = 2.22D-16  FOR DEC 10             (APPROX.) 2**-52    
C            = 7.11D-15  FOR VAX 11/780         (APPROX.) 2**-47    
C            = 1.14D-13  FOR IBM 370/158        (APPROX.) 2**-43    
C       
C             *** SHOULD BE CHANGED FOR OTHER MACHINES ***
C       
C     TO FACILITATE CONVERGENCE, RPARM(1) SHOULD BE SET TO
C          500.*DRELPR OR LARGER      
C       
      DRELPR = 7.11D-15
C       
      IPARM(1) = 100
      IPARM(2) = 0
      IPARM(3) = 0
      IPARM(4) = 6
      IPARM(5) = 0
      IPARM(6) = 1
      IPARM(7) = 1
      IPARM(8) = 0
      IPARM(9) = -1 
      IPARM(10) = 0 
      IPARM(11) = 0 
      IPARM(12) = 0 
C       
      RPARM(1) = 0.5D-5     
      RPARM(2) = 0.D0       
      RPARM(3) = 0.D0       
      RPARM(4) = .75D0      
      RPARM(5) = 1.D0       
      RPARM(6) = 0.D0       
      RPARM(7) = .25D0      
      RPARM(8) = 1.D2*DRELPR
      RPARM(9) = 0.D0       
      RPARM(10) = 0.D0      
      RPARM(11) = 0.D0      
      RPARM(12) = 0.D0      
C       
      RETURN      
      END 
      SUBROUTINE ECHALL (NN,IA,JA,A,RHS,IPARM,RPARM,ICALL)
C       
C ... THIS ROUTINE INITIALIZES THE ITPACK COMMON BLOCKS FROM THE    
C ... INFORMATION CONTAINED IN IPARM AND RPARM. ECHALL ALSO PRINTS THE
C ... VALUES OF ALL THE PARAMETERS IN IPARM AND RPARM.    
C       
C ... PARAMETER LIST:       
C       
C          IPARM  
C           AND   
C          RPARM  ARRAYS OF PARAMETERS SPECIFYING OPTIONS AND       
C                    TOLERANCES       
C          ICALL  INDICATOR OF WHICH PARAMETERS ARE BEING PRINTED   
C                    ICALL = 1,  INITIAL PARAMETERS       
C                    ICALL = 2,  FINAL PARAMETERS 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),IPARM(12),NN,ICALL    
      DOUBLE PRECISION A(1),RHS(NN),RPARM(12)   
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,N,NP1,NZRO  
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      IF (ICALL.NE.1) GO TO 100       
      N = NN      
      NP1 = N+1   
      NZRO = IA(NP1)-1      
C       
C ... INITIALIZE ITPACK COMMON
C       
      ZETA = RPARM(1)       
      CME = RPARM(2)
      SME = RPARM(3)
      FF = RPARM(4) 
      OMEGA = RPARM(5)      
      SPECR = RPARM(6)      
      BETAB = RPARM(7)      
      ITMAX = IPARM(1)      
      LEVEL = IPARM(2)      
      ISYM = IPARM(5)       
C       
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      BETADT = .FALSE.      
      IF (IPARM(6).EQ.1.OR.IPARM(6).EQ.3) ADAPT = .TRUE.  
      IF (IPARM(6).EQ.1) BETADT = .TRUE.
      IF (IPARM(6).EQ.2) PARTAD = .TRUE.
C       
      CASEII = .FALSE.      
      IF (IPARM(7).EQ.2) CASEII = .TRUE.
      IF (CASEII) SME = -CME
      IF (.NOT.CASEII.AND.SME.EQ.0.D0) SME = -1.D0
      SPR = SME   
C       
C ... SET REST OF COMMON VARIABLES TO ZERO      
C       
      IN = 0      
      IS = 0      
      HALT = .FALSE.
      BDELNM = 0.D0 
      DELNNM = 0.D0 
      DELSNM = 0.D0 
      GAMMA = 0.D0
      QA = 0.D0   
      QT = 0.D0   
      RHO = 0.D0  
      RRR = 0.D0  
      SIGE = 0.D0 
      STPTST = 0.D0 
      UDNM = 0.D0 
C       
      IF (LEVEL.LE.4) GO TO 80
C       
C     THIS SECTION OF ECHALL CAUSES PRINTING OF THE LINEAR SYSTEM AND 
C     THE ITERATIVE PARAMETERS
C       
      WRITE (NOUT,10)       
   10 FORMAT (///30X,'THE LINEAR SYSTEM IS AS FOLLOWS')   
      WRITE (NOUT,20)       
   20 FORMAT (/2X,'IA ARRAY') 
      WRITE (NOUT,30) (IA(I),I=1,NP1) 
   30 FORMAT (2X,10(2X,I8)) 
      WRITE (NOUT,40)       
   40 FORMAT (/2X,'JA ARRAY') 
      WRITE (NOUT,30) (JA(I),I=1,NZRO)
      WRITE (NOUT,50)       
   50 FORMAT (/2X,' A ARRAY') 
      WRITE (NOUT,60) (A(I),I=1,NZRO) 
   60 FORMAT (2X,5(2X,D20.13))
      WRITE (NOUT,70)       
   70 FORMAT (/2X,'RHS ARRAY')
      WRITE (NOUT,60) (RHS(I),I=1,N)  
   80 WRITE (NOUT,90)       
   90 FORMAT (///30X,'INITIAL ITERATIVE PARAMETERS')      
      GO TO 120   
  100 WRITE (NOUT,110)      
  110 FORMAT (///30X,'FINAL ITERATIVE PARAMETERS')
  120 WRITE (NOUT,130) IPARM(1),LEVEL,IPARM(3),NOUT,ISYM,IPARM(6)   
  130 FORMAT (35X,'IPARM(1)  =',I15,4X,'(ITMAX)'/35X,'IPARM(2)  =',I15, 
     *   4X,'(LEVEL) '/35X,'IPARM(3)  =',I15,4X,'(IRESET)'/35X,     
     *   'IPARM(4)  =',I15,4X,'(NOUT)  '/35X,'IPARM(5)  =',I15,4X,  
     *   '(ISYM)  '/35X,'IPARM(6)  =',I15,4X,'(IADAPT)')  
      WRITE (NOUT,140) IPARM(7),IPARM(8),IPARM(9),IPARM(10),IPARM(11),
     *   IPARM(12)
  140 FORMAT (35X,'IPARM(7)  =',I15,4X,'(ICASE)'/35X,'IPARM(8)  =',I15, 
     *   4X,'(NWKSP)'/35X,'IPARM(9)  =',I15,4X,'(NB)    '/35X,      
     *   'IPARM(10) =',I15,4X,'(IREMOVE)'/35X,'IPARM(11) =',I15,4X, 
     *   '(ITIME)'/35X,'IPARM(12) =',I15,4X,'(IDGTS)')    
      WRITE (NOUT,150) ZETA,CME,SME,FF,OMEGA,SPECR
  150 FORMAT (35X,'RPARM(1)  =',D15.8,4X,'(ZETA)  '/35X,'RPARM(2)  =',
     *   D15.8,4X,'(CME)   '/35X,'RPARM(3)  =',D15.8,4X,'(SME)   '/35X, 
     *   'RPARM(4)  =',D15.8,4X,'(FF)    '/35X,'RPARM(5)  =',D15.8,4X,
     *   '(OMEGA) '/35X,'RPARM(6)  =',D15.8,4X,'(SPECR) ')
      WRITE (NOUT,160) BETAB,RPARM(8),RPARM(9),RPARM(10),RPARM(11), 
     *   RPARM(12)
  160 FORMAT (35X,'RPARM(7)  =',D15.8,4X,'(BETAB) '/35X,'RPARM(8)  =',
     *   D15.8,4X,'(TOL)'/35X,'RPARM(9)  =',D15.8,4X,'(TIME1)'/35X, 
     *   'RPARM(10) =',D15.8,4X,'(TIME2)'/35X,'RPARM(11) =',D15.8,4X, 
     *   '(DIGIT1)'/35X,'RPARM(12) =',D15.8,4X,'(DIGIT2)')
C       
      RETURN      
      END 
      SUBROUTINE ECHOUT (IPARM,RPARM,IMTHD)     
C       
C     THIS ROUTINE INITIALIZES THE ITPACK COMMON BLOCKS FROM THE    
C     INFORMATION CONTAINED IN IPARM AND RPARM. 
C       
C ... PARAMETER LIST:       
C       
C          IPARM  
C           AND   
C          RPARM  ARRAYS OF PARAMETERS SPECIFYING OPTIONS AND       
C                    TOLERANCES       
C          IMTHD  INDICATOR OF METHOD 
C                    IMTHD = 1,  JCG  
C                    IMTHD = 2,  JSI  
C                    IMTHD = 3,  SOR  
C                    IMTHD = 4,  SSORCG 
C                    IMTHD = 5,  SSORSI 
C                    IMTHD = 6,  RSCG 
C                    IMTHD = 7,  RSSI 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IPARM(12),IMTHD 
      DOUBLE PRECISION RPARM(12)      
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
C ... INITIALIZE ITPACK COMMON
C       
      ZETA = RPARM(1)       
      CME = RPARM(2)
      SME = RPARM(3)
      FF = RPARM(4) 
      OMEGA = RPARM(5)      
      SPECR = RPARM(6)      
      BETAB = RPARM(7)      
      ITMAX = IPARM(1)      
      LEVEL = IPARM(2)      
      ISYM = IPARM(5)       
C       
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      BETADT = .FALSE.      
      IF (IPARM(6).EQ.1.OR.IPARM(6).EQ.3) ADAPT = .TRUE.  
      IF (IPARM(6).EQ.1) BETADT = .TRUE.
      IF (IPARM(6).EQ.2) PARTAD = .TRUE.
C       
      CASEII = .FALSE.      
      IF (IPARM(7).EQ.2) CASEII = .TRUE.
      IF (CASEII) SME = -CME
      IF (.NOT.CASEII.AND.SME.EQ.0.D0) SME = -1.D0
      SPR = SME   
C       
C ... SET REST OF COMMON VARIABLES TO ZERO      
C       
      IN = 0      
      IS = 0      
      HALT = .FALSE.
      BDELNM = 0.D0 
      DELNNM = 0.D0 
      DELSNM = 0.D0 
      GAMMA = 0.D0
      QA = 0.D0   
      QT = 0.D0   
      RHO = 0.D0  
      RRR = 0.D0  
      SIGE = 0.D0 
      STPTST = 0.D0 
      UDNM = 0.D0 
      IF (LEVEL.LE.2) RETURN
C       
C ... THIS SECTION OF ECHOUT ECHOES THE INPUT VALUES FOR THE INITIAL
C     ITERATIVE PARAMETERS  
C       
      WRITE (NOUT,10) ISYM,ITMAX,ZETA,ADAPT,CASEII
   10 FORMAT (///30X,'INITIAL ITERATIVE PARAMETERS',3X,   
     *   'RELEVANT SWITCHES'/35X,'ISYM   =',I15,8X,'IPARM(5)'/35X,  
     *   'ITMAX  =',I15,8X,'IPARM(1)'/35X,'ZETA   =',D15.8,8X,'RPARM(1)'
     *   /35X,'ADAPT  =',L15,8X,'IPARM(6)'/35X,'CASEII =',L15,8X,   
     *   'IPARM(7)')
      GO TO (80,20,100,60,40,80,20), IMTHD      
C       
C ... JSI, RSSI   
C       
   20 WRITE (NOUT,30) FF,CME,SME      
   30 FORMAT (35X,'FF     =',D15.8,8X,'RPARM(4)'/35X,'CME    =',D15.8,8X
     *   ,'RPARM(2)'/35X,'SME    =',D15.8,8X,'RPARM(3)'///) 
      RETURN      
C       
C ... SSORSI      
C       
   40 WRITE (NOUT,50) PARTAD,FF,CME,OMEGA,SPECR,BETAB,BETADT
   50 FORMAT (35X,'PARTAD =',L15,8X,'IPARM(6)'/35X,'FF     =',D15.8,8X, 
     *   'RPARM(4)'/35X,'CME    =',D15.8,8X,'RPARM(2)'/35X,'OMEGA  =',
     *   D15.8,8X,'RPARM(5)'/35X,'SPECR  =',D15.8,8X,'RPARM(6)'/35X,
     *   'BETAB  =',D15.8,8X,'RPARM(7)'/35X,'BETADT =',L15,8X,'IPARM(6)'
     *   ///)     
      RETURN      
C       
C ... SSORCG      
C       
   60 WRITE (NOUT,70) PARTAD,CME,OMEGA,SPECR,BETAB,BETADT 
   70 FORMAT (35X,'PARTAD =',L15,8X,'IPARM(6)'/35X,'CME    =',D15.8,8X, 
     *   'RPARM(2)'/35X,'OMEGA  =',D15.8,8X,'RPARM(5)'/35X,'SPECR  =',
     *   D15.8,8X,'RPARM(6)'/35X,'BETAB  =',D15.8,8X,'RPARM(7)'/35X,
     *   'BETADT =',L15,8X,'IPARM(6)'///)       
      RETURN      
C       
C ... JCG, RSCG   
C       
   80 IF (ADAPT) RETURN     
      WRITE (NOUT,90) CME   
   90 FORMAT (35X,'CME    =',D15.8,8X,'RPARM(2)'///)      
C       
  100 CONTINUE    
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION EIGVNS (N,TRI,D,E2,IER)   
C       
C     COMPUTES THE LARGEST EIGENVALUE OF A SYMMETRIC TRIDIAGONAL MATRIX 
C     FOR CONJUGATE GRADIENT ACCELERATION.      
C       
C ... PARAMETER LIST:       
C       
C          N      ORDER OF TRIDIAGONAL SYSTEM   
C          TRI    SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N 
C          D      ARRAY FOR EQRT1S (NEGATIVE DIAGONAL ELEMENTS)     
C          E2     ARRAY FOR EQRT1S (SUPER DIAGONAL ELEMENTS)
C          IER    ERROR FLAG: ON RETURN, IER=0 INDICATES THAT       
C                    THE LARGEST EIGENVALUE OF TRI WAS FOUND.       
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,IER 
      DOUBLE PRECISION TRI(2,1),D(N),E2(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I   
C       
      EIGVNS = 0.D0 
C       
      D(1) = -TRI(1,1)      
      DO 10 I = 2,N 
         D(I) = -TRI(1,I)   
         E2(I) = DABS(TRI(2,I))       
   10 CONTINUE    
C       
      CALL EQRT1S (D,E2,N,1,0,IER)    
      EIGVNS = -D(1)
C       
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION EIGVSS (N,TRI,START,ZETA,ITMAX,IER) 
C       
C     COMPUTES THE LARGEST EIGENVALUE OF A SYMMETRIC TRIDIAGONAL MATRIX 
C     FOR CONJUGATE GRADIENT ACCELERATION.      
C     MODIFIED IMSL ROUTINE ZBRENT USED.
C       
C ... PARAMETER LIST:       
C       
C          N      ORDER OF TRIDIAGONAL SYSTEM   
C          TRI    SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N 
C          START  INITIAL LOWER BOUND OF INTERVAL CONTAINING ROOT   
C          ZETA   STOPPING CRITERIA   
C          IER    ERROR FLAG: ON RETURN, IER=0 INDICATES THAT       
C                    THE LARGEST EIGENVALUE OF TRI WAS FOUND.       
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,ITMAX,IER   
      DOUBLE PRECISION TRI(2,1),START,ZETA,A,B,EPS
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER MAXFN,NSIG,ITMP 
C       
      EIGVSS = 0.D0 
      ITMP = IFIX(SNGL(-DLOG10(DABS(ZETA))))    
      NSIG = MAX0(ITMP,4)   
      MAXFN = MAX0(ITMAX,50)
C       
C     EPS = DMIN1(ZETA,0.5D-4)
C       
      EPS = 0.0D0 
      A = START   
      B = 1.0D0   
      CALL ZBRENT (N,TRI,EPS,NSIG,A,B,MAXFN,IER)
      EIGVSS = B  
C       
      RETURN      
      END 
      SUBROUTINE EQRT1S (D,E2,NN,M,ISW,IERR)    
C       
C   MODIFIED IMSL ROUTINE NAME   - EQRT1S       
C       
C-----------------------------------------------------------------------
C       
C   COMPUTER            - CDC/SINGLE  
C       
C   LATEST REVISION     - JUNE 1, 1980
C       
C   PURPOSE             - SMALLEST OR LARGEST M EIGENVALUES OF A    
C                           SYMMETRIC TRIDIAGONAL MATRIX  
C       
C   USAGE               - CALL EQRT1S (D,E2,N,M,ISW,IER)  
C       
C   ARGUMENTS    D      - INPUT VECTOR OF LENGTH N CONTAINING       
C                           THE DIAGONAL ELEMENTS OF THE MATRIX.  THE 
C                           COMPUTED EIGENVALUES REPLACE THE FIRST M
C                           COMPONENTS OF THE VECTOR D IN NON-      
C                           DECREASING SEQUENCE, WHILE THE REMAINING
C                           COMPONENTS ARE LOST.
C                E2     - INPUT VECTOR OF LENGTH N CONTAINING       
C                           THE SQUARES OF THE OFF-DIAGONAL ELEMENTS
C                           OF THE MATRIX.  INPUT E2 IS DESTROYED.  
C                N      - INPUT SCALAR CONTAINING THE ORDER OF THE  
C                           MATRIX. (= NN)      
C                M      - INPUT SCALAR CONTAINING THE NUMBER OF     
C                           SMALLEST EIGENVALUES DESIRED (M IS      
C                           LESS THAN OR EQUAL TO N).     
C                ISW    - INPUT SCALAR MEANING AS FOLLOWS - 
C                           ISW=1 MEANS THAT THE MATRIX IS KNOWN TO BE
C                             POSITIVE DEFINITE.
C                           ISW=0 MEANS THAT THE MATRIX IS NOT KNOWN
C                             TO BE POSITIVE DEFINITE.    
C                IER    - ERROR PARAMETER. (OUTPUT) (= IERR)
C                           WARNING ERROR       
C                             IER = 601 INDICATES THAT SUCCESSIVE   
C                               ITERATES TO THE K-TH EIGENVALUE WERE NOT
C                               MONOTONE INCREASING. THE VALUE K IS 
C                               STORED IN E2(1).
C                           TERMINAL ERROR      
C                             IER = 602 INDICATES THAT ISW=1 BUT MATRIX 
C                               IS NOT POSITIVE DEFINITE  
C       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32 
C                       - SINGLE/H36,H48,H60    
C       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND       
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL  
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C       
C   REMARKS      AS WRITTEN, THE ROUTINE COMPUTES THE M SMALLEST    
C                EIGENVALUES. TO COMPUTE THE M LARGEST EIGENVALUES, 
C                REVERSE THE SIGN OF EACH ELEMENT OF D BEFORE AND   
C                AFTER CALLING THE ROUTINE. IN THIS CASE, ISW MUST  
C                EQUAL ZERO.
C       
C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.   
C       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.    
C       
C-----------------------------------------------------------------------
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
C                                  SPECIFICATIONS FOR ARGUMENTS     
C       
      INTEGER NN,M,ISW,IERR 
      DOUBLE PRECISION D(NN),E2(NN)   
C       
C                                  SPECIFICATIONS FOR LOCAL VARIABLES 
C       
      INTEGER II,I,JJ,J,K1,K,N,IER    
      DOUBLE PRECISION DELTA,DLAM,EP,ERR,F,P,QP,Q,R,S,TOT 
C       
C                                  DRELPR = MACHINE PRECISION       
C                                  FIRST EXECUTABLE STATEMENT       
C       
      N = NN      
      IER = 0     
      DLAM = 0.0D0
      ERR = 0.0D0 
      S = 0.0D0   
C       
C                                  LOOK FOR SMALL SUB-DIAGONAL ENTRIES
C                                  DEFINE INITIAL SHIFT FROM LOWER  
C                                  GERSCHGORIN BOUND.     
C       
      TOT = D(1)  
      Q = 0.0D0   
      J = 0       
      DO 30 I = 1,N 
         P = Q    
         IF (I.EQ.1) GO TO 10 
         IF (P.GT.DRELPR*(DABS(D(I))+DABS(D(I-1)))) GO TO 20
   10    E2(I) = 0.0D0      
C       
C                                  COUNT IF E2(I) HAS UNDERFLOWED   
C       
   20    IF (E2(I).EQ.0.D0) J = J+1   
         Q = 0.0D0
         IF (I.NE.N) Q = DSQRT(DABS(E2(I+1)))   
         TOT = DMIN1(D(I)-P-Q,TOT)    
   30 CONTINUE    
      IF (ISW.EQ.1.AND.TOT.LT.0.0D0) GO TO 50   
      DO 40 I = 1,N 
         D(I) = D(I)-TOT    
   40 CONTINUE    
      GO TO 60    
   50 TOT = 0.0D0 
   60 DO 200 K = 1,M
C       
C                                  NEXT QR TRANSFORMATION 
C       
   70    TOT = TOT+S
         DELTA = D(N)-S     
         I = N    
         F = DABS(DRELPR*TOT) 
         IF (DLAM.LT.F) DLAM = F      
         IF (DELTA.GT.DLAM) GO TO 90  
         IF (DELTA.GE.(-DLAM)) GO TO 170
         IER = 602
         IF (LEVEL.GE.1) WRITE (NOUT,80)
   80    FORMAT ('0','*** W A R N I N G ************'/' ',
     *      '    IN ITPACK ROUTINE EQRT1S  '/' ', 
     *      '    PARAMETER ISW = 1 BUT MATRIX   '/' ',    
     *      '    NOT POSITIVE DEFINITE')
         GO TO 210
C       
C                                  REPLACE SMALL SUB-DIAGONAL SQUARES 
C                                  BY ZERO TO REDUCE THE INCIDENCE OF 
C                                  UNDERFLOWS   
C       
   90    IF (K.EQ.N) GO TO 110
         K1 = K+1 
         DO 100 J = K1,N    
            IF (E2(J).LE.(DRELPR*(D(J)+D(J-1)))**2) E2(J) = 0.0D0   
  100    CONTINUE 
  110    F = E2(N)/DELTA    
         QP = DELTA+F       
         P = 1.0D0
         IF (K.EQ.N) GO TO 140
         K1 = N-K 
         DO 130 II = 1,K1   
            I = N-II
            Q = D(I)-S-F    
            R = Q/QP
            P = P*R+1.0D0   
            EP = F*R
            D(I+1) = QP+EP  
            DELTA = Q-EP    
            IF (DELTA.GT.DLAM) GO TO 120
            IF (DELTA.GE.(-DLAM)) GO TO 170     
            IER = 602       
            IF (LEVEL.GE.0) WRITE (NOUT,80)     
            GO TO 210       
  120       F = E2(I)/Q     
            QP = DELTA+F    
            E2(I+1) = QP*EP 
  130    CONTINUE 
  140    D(K) = QP
         S = QP/P 
         IF (TOT+S.GT.TOT) GO TO 70   
         IER = 601
         E2(1) = K
         IF (LEVEL.GE.1) WRITE (NOUT,150) K     
  150    FORMAT ('0','*** W A R N I N G ************'/'0',
     *      '    IN ITPACK ROUTINE EQRT1S  '/' ', 
     *      '    SUCCESSIVE ITERATES TO THE',I10/' ',     
     *      '    EIGENVALUE WERE NOT MONOTONE INCREASING ') 
C       
C                                  SET ERROR -- IRREGULAR END       
C                                  DEFLATE MINIMUM DIAGONAL ELEMENT 
C       
         S = 0.0D0
         DELTA = QP 
         DO 160 J = K,N     
            IF (D(J).GT.DELTA) GO TO 160
            I = J 
            DELTA = D(J)    
  160    CONTINUE 
C       
C                                  CONVERGENCE  
C       
  170    IF (I.LT.N) E2(I+1) = E2(I)*F/QP       
         IF (I.EQ.K) GO TO 190
         K1 = I-K 
         DO 180 JJ = 1,K1   
            J = I-JJ
            D(J+1) = D(J)-S 
            E2(J+1) = E2(J) 
  180    CONTINUE 
  190    D(K) = TOT 
         ERR = ERR+DABS(DELTA)
         E2(K) = ERR
  200 CONTINUE    
      IF (IER.EQ.0) GO TO 220 
  210 CONTINUE    
  220 IERR = IER  
      RETURN      
      END 
      INTEGER FUNCTION IPSTR (OMEGA)  
C       
C     FINDS THE SMALLEST INTEGER, IPSTR, GREATER THAN 5 SUCH THAT   
C          IPSTR * (OMEGA-1)**(IPSTR-1) .LE. 0.50. IPSTR WILL BE SET
C          IN LOOP. 
C       
C ... PARAMETER LIST:       
C       
C          OMEGA  RELAXATION FACTOR FOR SOR METHOD
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      DOUBLE PRECISION OMEGA
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IP  
      DOUBLE PRECISION WM1  
C       
      WM1 = OMEGA-1.D0      
C       
      DO 10 IP = 6,940      
         IF (DBLE(FLOAT(IP))*(WM1**(IP-1)).GT.0.50D0) GO TO 10      
         IPSTR = IP 
         RETURN   
   10 CONTINUE    
      IPSTR = 940 
      RETURN      
C       
      END 
      SUBROUTINE ITERM (NN,A,U,WK,IMTHDD)       
C       
C     THIS ROUTINE PRODUCES THE ITERATION SUMMARY LINE AT THE END   
C     OF EACH ITERATION. IF LEVEL = 5, THE LATEST APPROXIMATION     
C     TO THE SOLUTION WILL BE PRINTED.
C       
C ... PARAMETER LIST:       
C       
C          NN     ORDER OF SYSTEM OR, FOR REDUCED SYSTEM  
C                    ROUTINES, ORDER OF BLACK SUBSYSTEM   
C          A      ITERATION MATRIX    
C          U      SOLUTION ESTIMATE   
C          WK     WORK ARRAY OF LENGTH NN       
C          IMTHD  INDICATOR OF METHOD (=IMTHDD) 
C                    IMTHD = 1,  JCG  
C                    IMTHD = 2,  JSI  
C                    IMTHD = 3,  SOR  
C                    IMTHD = 4,  SSORCG 
C                    IMTHD = 5,  SSORSI 
C                    IMTHD = 6,  RSCG 
C                    IMTHD = 7,  RSSI 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER NN,IMTHD      
      DOUBLE PRECISION A(1),U(NN),WK(NN)
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IMTHDD,IP,N 
      DOUBLE PRECISION QTFF 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
      N = NN      
      IMTHD = IMTHDD
C       
C ... PRINT VARIOUS PARAMETERS AFTER EACH ITERATION       
C       
      IF (LEVEL.LT.2) RETURN
      GO TO (10,110,170,210,50,10,110), IMTHD   
   10 IF (IN.GT.0) GO TO 30 
C       
C ... PRINT HEADER FOR JCG AND RSCG   
C       
      WRITE (NOUT,20)       
   20 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',5X,'CONVERGENCE',7X,'CME ',11X,'RHO',12X,'GAMMA'/ 
     *   ' ITERATIONS',4X,'TEST '//)  
C       
C ... PRINT SUMMARY LINE    
C       
   30 WRITE (NOUT,40) IN,STPTST,CME,RHO,GAMMA   
   40 FORMAT (4X,I5,3X,4D15.7)
      IF (LEVEL.GE.4) GO TO 250       
C       
      RETURN      
C       
   50 IF (IN.GT.0) GO TO 70 
C       
C ... PRINT HEADER FOR SSOR-SI
C       
      WRITE (NOUT,60)       
   60 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',4X,'CONVERGENCE',7X,'PARAMETER CHANGE TEST',10X,
     *   'RHO',12X,'GAMMA'/' ITERATIONS',3X,'TEST ',11X,'LHS(QA)',7X, 
     *   'RHS(QT**FF)'//)   
C       
C ... PRINT SUMMARY LINE    
C       
   70 IP = IN-IS  
      IF (IMTHD.EQ.7) IP = 2*IP       
      IF (IP.LT.3) GO TO 90 
      QTFF = QT**FF 
      WRITE (NOUT,80) IN,STPTST,QA,QTFF,RHO,GAMMA 
   80 FORMAT (4X,I5,3X,5D15.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
   90 WRITE (NOUT,100) IN,STPTST,RHO,GAMMA      
  100 FORMAT (4X,I5,3X,D15.7,30X,2D15.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
  110 IF (IN.GT.0) GO TO 130
C       
C ... PRINT HEADER FOR J-SI AND RS-SI 
C       
      WRITE (NOUT,120)      
  120 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',4X,'CONVERGENCE',7X,'PARAMETER CHANGE TEST',10X,
     *   'RHO'/' ITERATIONS',3X,'TEST ',11X,'LHS(QA)',7X,'RHS(QT**FF)'//
     *   )
C       
C ... PRINT SUMMARY LINE    
C       
  130 IP = IN-IS  
      IF (IMTHD.EQ.7) IP = 2*IP       
      IF (IP.LT.3) GO TO 150
      QTFF = QT**FF 
      WRITE (NOUT,140) IN,STPTST,QA,QTFF,RHO    
  140 FORMAT (4X,I5,3X,5D15.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
  150 WRITE (NOUT,160) IN,STPTST,RHO  
  160 FORMAT (4X,I5,3X,D15.7,30X,D15.7) 
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
C ... PRINT VARIOUS PARAMETERS AFTER EACH ITERATION FOR SOR.
C       
  170 IF (IN.GT.0) GO TO 190
C       
C ... PRINT HEADER FOR SOR  
C       
      WRITE (NOUT,180)      
  180 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',4X,'CONVERGENCE',6X,'CME ',9X,'OMEGA',7X,     
     *   'SPECTRAL'/' ITERATIONS',3X,'TEST',38X,'RADIUS'//) 
C       
C ... PRINT SUMMARY LINE FOR SOR      
C       
  190 CONTINUE    
      WRITE (NOUT,200) IN,STPTST,CME,OMEGA,SPECR
  200 FORMAT (4X,I5,3X,4D14.7)
      IF (LEVEL.GE.4) GO TO 250       
C       
      RETURN      
C       
C ... PRINT VARIOUS PARAMETERS AFTER EACH ITERATION FOR SSOR-CG.    
C       
  210 IF (IN.GT.0) GO TO 230
C       
C ... PRINT HEADER FOR SSOR-CG
C       
      WRITE (NOUT,220)      
  220 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  
     *   ' NUMBER OF',4X,'CONVERGENCE',3X,' SPECTRAL',6X,'S-PRIME',9X,
     *   'RHO',10X,'GAMMA'/' ITERATIONS',3X,'TEST ',10X,'RADIUS'//) 
C       
C ... PRINT SUMMARY LINE FOR SSOR-CG  
C       
  230 CONTINUE    
      WRITE (NOUT,240) IN,STPTST,SPECR,SPR,RHO,GAMMA      
  240 FORMAT (4X,I5,3X,5D14.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
C       
  250 IF (IMTHD.GT.5) GO TO 270       
      WRITE (NOUT,260) IN   
  260 FORMAT ('0',2X,'ESTIMATE OF SOLUTION AT ITERATION ',I5)       
      GO TO 290   
  270 WRITE (NOUT,280) IN   
  280 FORMAT ('0',2X,'ESTIMATE OF SOLUTION AT BLACK POINTS ',       
     *   'AT ITERATION ',I5)
  290 DO 300 I = 1,N
         WK(I) = U(I)/A(I)  
  300 CONTINUE    
      WRITE (NOUT,310) (WK(I),I=1,N)  
  310 FORMAT (2X,5(2X,D20.13))
      WRITE (NOUT,320)      
  320 FORMAT (//) 
C       
      RETURN      
      END 
      SUBROUTINE IVFILL (N,IV,IVAL)   
C       
C     FILLS AN INTEGER VECTOR, IV, WITH AN INTEGER VALUE, IVAL.     
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTOR IV   
C          IV     INTEGER VECTOR      
C          IVAL   INTEGER CONSTANT THAT FILLS FIRST N LOCATIONS OF IV 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,IVAL,IV(N)  
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,M,MP1       
C       
      IF (N.LE.0) RETURN    
C       
C     CLEAN UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 10  
C       
      M = MOD(N,10) 
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         IV(I) = IVAL       
   10 CONTINUE    
      IF (N.LT.10) RETURN   
C       
   20 MP1 = M+1   
      DO 30 I = MP1,N,10    
         IV(I) = IVAL       
         IV(I+1) = IVAL     
         IV(I+2) = IVAL     
         IV(I+3) = IVAL     
         IV(I+4) = IVAL     
         IV(I+5) = IVAL     
         IV(I+6) = IVAL     
         IV(I+7) = IVAL     
         IV(I+8) = IVAL     
         IV(I+9) = IVAL     
   30 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE OMEG (DNRM,IFLAG)    
C       
C     COMPUTES NEW VALUES FOR  CME, OMEGA, AND SPECR FOR  
C     FULLY ADAPTIVE SSOR METHODS.    
C       
C ... PARAMETER LIST:       
C       
C          DNRM   NUMERATOR OF RAYLEIGH QUOTIENT
C          IFLAG  INDICATOR OF APPROPRIATE ENTRY POINT    
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IFLAG 
      DOUBLE PRECISION DNRM 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      DOUBLE PRECISION TEMP,ZM1,ZM2   
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      ZM1 = 0.D0  
      ZM2 = 0.D0  
      IF (IFLAG.EQ.1) GO TO 10
C       
C ... IFLAG .NE. 1, COMPUTE NEW ESTIMATE FOR CME
C       
      ZM1 = ((1.D0-SPR)*(1.D0+BETAB*OMEGA**2)-OMEGA*(2.D0-OMEGA))/(OMEGA
     *   *(OMEGA-1.D0-SPR)) 
C       
      IF (.NOT.CASEII) ZM2 = DNRM/BDELNM
      IF (CASEII) ZM2 = DSQRT(DABS(DNRM/BDELNM))
      CME = DMAX1(CME,ZM1,ZM2)
C       
C ... IFLAG = 1, OR CONTINUATION OF IFLAG .NE. 1
C       
C        COMPUTE NEW VALUES OF OMEGA AND SPECR BASED ON CME AND BETAB 
C       
   10 IS = IN+1   
      DELSNM = DELNNM       
      IF (CME.GE.(4.D0*BETAB)) GO TO 30 
C       
C ... CME .LT. 4.D0*BETAB   
C       
      TEMP = DSQRT(DABS(1.D0-2.D0*CME+4.D0*BETAB))
      OMEGA = DMAX1((2.D0/(1.D0+TEMP)),1.D0)    
      TEMP = (1.D0-CME)/TEMP
      SPECR = (1.D0-TEMP)/(1.D0+TEMP) 
      IF (DABS(OMEGA-1.D0).LT.DRELPR) SPECR = 0.D0
      IF (LEVEL.GE.2) WRITE (NOUT,20) IN,BETAB,ZM1,ZM2,CME,OMEGA,SPECR
   20 FORMAT (/30X,'PARAMETERS WERE CHANGED AT ITERATION NO.',I5/35X, 
     *   'NEW ESTIMATE OF BETAB            =',D15.7/35X,  
     *   'SOLUTION TO CHEBYSHEV EQN.       =',D15.7/35X,  
     *   'SOLUTION TO RAYLEIGH QUOTIENT    =',D15.7/35X,  
     *   'NEW ESTIMATE FOR CME             =',D15.7/35X,  
     *   'NEW ESTIMATE FOR OMEGA           =',D15.7/35X,  
     *   'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)     
C       
      RETURN      
C       
C ... CME .GE. 4.D0*BETAB   
C       
C ... OMEGA-STAR WILL BE CHOSEN       
C       
   30 CME = 2.D0*DSQRT(DABS(BETAB))   
      OMEGA = 2.D0/(1.D0+DSQRT(DABS(1.D0-4.D0*BETAB)))    
      SPECR = OMEGA-1.D0    
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      IF (LEVEL.GE.2) WRITE (NOUT,20) IN,BETAB,ZM1,ZM2,CME,OMEGA,SPECR
C       
      RETURN      
      END 
      LOGICAL FUNCTION OMGCHG (NDUMMY)
C       
C ... THIS FUNCTION TESTS TO SEE WHETHER OMEGA SHOULD BE CHANGED    
C ... FOR SSOR CG METHOD.   
C       
C ... PARAMETER LIST:       
C       
C          NDUMMY ARBITRARY INTEGER PARAMETER   
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER NDUMMY
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      DOUBLE PRECISION DEL1,DEL2,X    
C       
      DOUBLE PRECISION PHI  
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
C ... STATEMENT FUNCTION PHI(X)       
C       
      PHI(X) = (1.D0-DSQRT(DABS(1.D0-X)))/(1.D0+DSQRT(DABS(1.D0-X)))
C       
      OMGCHG = .FALSE.      
      IF (IN-IS.LT.3) RETURN
      IF (SPECR.EQ.0.D0) GO TO 10     
      IF (SPECR.GE.SPR) RETURN
      DEL1 = -DLOG(DABS(PHI(SPECR)/PHI(SPECR/SPR)))       
      DEL2 = -DLOG(DABS(PHI(SPR)))    
      IF ((DEL1/DEL2).GE.FF) RETURN   
C       
   10 OMGCHG = .TRUE.       
C       
      RETURN      
      END 
      LOGICAL FUNCTION OMGSTR (NDUMMY)
C       
C     TESTS FOR FULLY ADAPTIVE SSOR METHODS WHETHER OMEGA-STAR      
C     SHOULD BE USED FOR OMEGA AND THE ADAPTIVE PROCESS TURNED      
C     OFF.
C       
C ... PARAMETER LIST:       
C       
C          NDUMMY ARBITRARY INTEGER PARAMETER   
C       
C ... SPECIFICATION FOR ARGUMENT      
C       
      INTEGER NDUMMY
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      DOUBLE PRECISION OMSTAR,TEMP,TEMP1,X      
C       
      DOUBLE PRECISION PHI  
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
C ... STATEMENT FUNCTION PHI(X)       
C       
      PHI(X) = (1.D0-DSQRT(DABS(1.D0-X)))/(1.D0+DSQRT(DABS(1.D0-X)))
C       
      OMGSTR = .FALSE.      
      IF (BETAB.GE..25D0.OR..NOT.ADAPT) RETURN  
      OMSTAR = 2.D0/(1.D0+DSQRT(DABS(1.D0-4.D0*BETAB)))   
C       
C ... TEST TO CHOSE OMEGA-STAR
C       
      IF ((OMSTAR.LE.1.D0).OR.(SPECR.LE.0.D0)) GO TO 10   
      TEMP = DLOG(DABS(PHI(OMSTAR-1.D0)))       
      TEMP1 = DLOG(DABS(PHI(SPECR)))  
      IF ((TEMP/TEMP1).LT.FF) RETURN  
C       
C ... OMEGA-STAR WAS CHOSEN 
C       
   10 OMEGA = OMSTAR
      SPECR = OMEGA-1.D0    
      OMGSTR = .TRUE.       
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      CME = 2.D0*DSQRT(DABS(BETAB))   
      RRR = PHI(1.D0-SPECR)**2
      GAMMA = 2.D0/(2.D0-SPECR)       
      SIGE = SPECR/(2.D0-SPECR)       
      RHO = 1.D0  
      IS = IN+1   
      DELSNM = DELNNM       
      IF (LEVEL.GE.2) WRITE (NOUT,20) IN,CME,OMEGA,SPECR  
   20 FORMAT (/30X,'OMEGA-STAR, AN ALTERNATE ESTIMATE OF',
     *   ' OMEGA, WAS CHOSEN AT ITERATION',I5/35X,
     *   'NEW ESTIMATE FOR CME             =',D15.7/35X,  
     *   'NEW ESTIMATE FOR OMEGA           =',D15.7/35X,  
     *   'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)     
C       
      RETURN      
      END 
      SUBROUTINE PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP,IBMTH)     
C       
C     COMPUTES ACCELERATION PARAMETERS FOR CONJUGATE GRADIENT       
C     ACCELERATED METHODS.  
C       
C ... PARAMETER LIST:       
C       
C          DTNRM  INNER PRODUCT OF RESIDUALS    
C          C1     OUTPUT: RHO*GAMMA   
C          C2     OUTPUT: RHO 
C          C3     OUTPUT: 1-RHO       
C          C4     OUTPUT: RHO*(1-GAMMA) 
C          GAMOLD OUTPUT: VALUE OF GAMMA AT PRECEDING ITERATION     
C          RHOTMP LAST ESTIMATE FOR VALUE OF RHO
C          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY CG 
C                      IBMTH = 1,   JACOBI      
C                            = 2,   REDUCED SYSTEM
C                            = 3,   SSOR
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IBMTH 
      DOUBLE PRECISION DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP    
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IP  
      DOUBLE PRECISION RHOOLD 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      IP = IN-IS  
C       
C ... SET RHOOLD AND GAMOLD 
C       
      RHOOLD = RHO
      GAMOLD = GAMMA
C       
C ... COMPUTE GAMMA (IN+1)  
C       
C ... FOR JACOBI OR REDUCED SYSTEM CG 
C       
      IF (IBMTH.LE.2) GAMMA = 1.D0/(1.D0-DTNRM/DELNNM)    
C       
C ... FOR SSOR CG 
C       
      IF (IBMTH.EQ.3) GAMMA = DELNNM/DTNRM      
C       
C ... COMPUTE RHO (IN+1)    
C       
      RHO = 1.D0  
      IF (IP.EQ.0) GO TO 20 
      IF (ISYM.EQ.0) GO TO 10 
      RHO = 1.D0/(1.D0-GAMMA*RHOTMP/DELSNM)     
      GO TO 20    
   10 RHO = 1.D0/(1.D0-GAMMA*DELNNM/(GAMOLD*DELSNM*RHOOLD)) 
C       
C ... COMPUTE CONSTANTS C1, C2, C3, AND C4      
C       
   20 DELSNM = DELNNM       
      RHOTMP = RHOOLD       
      C1 = RHO*GAMMA
      C2 = RHO    
      C3 = 1.D0-RHO 
      C4 = RHO*(1.D0-GAMMA) 
C       
      RETURN      
      END 
      SUBROUTINE PARSI (C1,C2,C3,IBMTH) 
C       
C     COMPUTES ACCELERATION PARAMETERS FOR SEMI-ITERATIVE 
C     ACCELERATED METHODS.  
C       
C ... PARAMETER LIST:       
C       
C          C1,C2  
C           AND   
C           C3    OUTPUT ACCELERATION PARAMETERS
C          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY SI 
C                      IBMTH = 1, JACOBI
C                            = 2, REDUCED SYSTEM
C                            = 3, SSOR
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IBMTH 
      DOUBLE PRECISION C1,C2,C3       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IP  
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      IP = IN-IS  
      IF (IP.EQ.0) GO TO 30 
      IF (IP.EQ.1) GO TO 10 
      RHO = 1.D0/(1.D0-SIGE*SIGE*RHO*.25D0)     
      GO TO 20    
   10 RHO = 1.D0/(1.D0-SIGE*SIGE*.5D0)
C       
   20 C1 = RHO*GAMMA
      C2 = RHO    
      C3 = 1.D0-RHO 
C       
      RETURN      
C       
C ... NONADAPTIVE INITIALIZATION FOR SEMI-ITERATIVE METHODS 
C       
   30 CONTINUE    
      GO TO (40,50,60), IBMTH 
C       
C ... JSI 
C       
   40 IF (CASEII) SME = -CME
      GAMMA = 2.D0/(2.D0-CME-SME)     
      SIGE = (CME-SME)/(2.D0-CME-SME) 
      GO TO 70    
C       
C ... REDUCED SYSTEM SI     
C       
   50 GAMMA = 2.D0/(2.D0-CME*CME)     
      SIGE = CME*CME/(2.D0-CME*CME)   
      RRR = (1.D0-DSQRT(DABS(1.D0-CME*CME)))/(1.D0+DSQRT(DABS(1.D0-CME* 
     *   CME)))   
      GO TO 70    
C       
C ... SSORSI      
C       
   60 GAMMA = 2.D0/(2.D0-SPECR)       
      SIGE = SPECR/(2.D0-SPECR)       
      RRR = (1.D0-DSQRT(DABS(1.D0-SIGE*SIGE)))/(1.D0+DSQRT(DABS(1.D0- 
     *   SIGE*SIGE)))       
C       
   70 RHO = 1.D0  
      C1 = GAMMA  
      C2 = 1.D0   
      C3 = 0.D0   
C       
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION PBETA (NN,IA,JA,A,V,W1,W2)
C       
C     ... COMPUTES THE NUMERATOR FOR THE COMPUTATION OF BETAB IN    
C     ...  SSOR METHODS.    
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          W1,W2  WORKSPACE VECTORS OF LENGTH N 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),V(NN),W1(NN),W2(NN) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IBGN,IEND,II,ITMP,JAI,JAJJ,JJ,K,N,NM1     
      DOUBLE PRECISION SUM,TEMP1,TEMP2
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      PBETA = 0.D0
      IF (ISYM.EQ.0) GO TO 110
C       
C     ************** NON - SYMMETRIC SECTION ********************   
C       
      DO 10 I = 1,N 
         W1(I) = V(I)       
   10 CONTINUE    
      TEMP1 = 0.D0
      TEMP2 = 0.D0
      ITMP = 2    
      IBGN = IA(1)
      IEND = IA(ITMP)-1     
      IF (IEND.LT.IBGN) GO TO 30      
      DO 20 I = IBGN,IEND   
         JAI = JA(I)
         TEMP1 = TEMP1-A(I)*W1(JAI)   
   20 CONTINUE    
   30 W1(1) = TEMP1 
      W2(1) = 0.D0
      NM1 = N-1   
      DO 70 K = 2,NM1       
         TEMP1 = 0.D0       
         TEMP2 = 0.D0       
         IBGN = IA(K)       
         IEND = IA(K+1)-1   
         IF (IEND.LT.IBGN) GO TO 60   
         DO 50 I = IBGN,IEND
            JAI = JA(I)     
            IF (JAI.GT.K) GO TO 40    
            TEMP2 = TEMP2-A(I)*W1(JAI)
            GO TO 50
   40       TEMP1 = TEMP1-A(I)*W1(JAI)
   50    CONTINUE 
   60    W1(K) = TEMP1      
         W2(K) = TEMP2      
   70 CONTINUE    
      TEMP2 = 0.D0
      IBGN = IA(N)
      IEND = IA(N+1)-1      
      IF (IEND.LT.IBGN) GO TO 90      
      DO 80 I = IBGN,IEND   
         JAI = JA(I)
         TEMP2 = TEMP2-A(I)*W1(JAI)   
   80 CONTINUE    
   90 W2(N) = TEMP2 
      DO 100 I = 1,N
         PBETA = PBETA+V(I)*W2(I)     
  100 CONTINUE    
      RETURN      
C       
C     **************** SYMMETRIC SECTION *************************  
C       
  110 DO 130 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = 0.D0 
         IF (IBGN.GT.IEND) GO TO 130  
         DO 120 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*V(JAJJ)   
  120    CONTINUE 
         PBETA = PBETA+SUM*SUM
  130 CONTINUE    
      RETURN      
C       
      END 
      SUBROUTINE PBSOR (NN,IA,JA,A,U,RHS)       
C       
C     ... THIS SUBROUTINE COMPUTES A BACKWARD SOR SWEEP.  
C       
C ... PARAMETER LIST:       
C       
C          N      ORDER OF SYSTEM (= NN)
C          OMEGA  RELAXATION FACTOR   
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          U      LATEST ESTIMATE OF SOLUTION   
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),RHS(NN)       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IBGN,IEND,II,JAJJ,JJ,N,NPL1     
      DOUBLE PRECISION OMM1,SUM,UI    
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      NPL1 = N+1  
      OMM1 = OMEGA-1.D0     
      IF (ISYM.EQ.0) GO TO 40 
C       
C     *************** NON - SYMMETRIC SECTION **********************
C       
      DO 30 I = 1,N 
         II = NPL1-I
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    U(II) = OMEGA*SUM-OMM1*U(II) 
   30 CONTINUE    
      RETURN      
C       
C     ***************** SYMMETRIC SECTION **************************
C       
   40 DO 60 II = 1,N
         UI = U(II) 
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UI      
   50    CONTINUE 
   60 CONTINUE    
C       
      DO 90 I = 1,N 
         II = NPL1-I
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 80   
         DO 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   70    CONTINUE 
   80    U(II) = OMEGA*SUM-OMM1*U(II) 
   90 CONTINUE    
      RETURN      
C       
      END 
      SUBROUTINE PERMAT (NN,IA,JA,A,P,NEWIA,ISYM,LEVEL,NOUT,IERR)   
C       
C*********************************************************************
C       
C ... SUBROUTINE PERMAT TAKES THE SPARSE MATRIX REPRESENTATION      
C     OF THE MATRIX STORED IN THE ARRAYS IA, JA, AND A AND
C     PERMUTES BOTH ROWS AND COLUMNS OVERWRITING THE PREVIOUS       
C     STRUCTURE.  
C       
C ... PARAMETER LIST:       
C       
C         N      ORDER OF SYSTEM (= NN) 
C         IA,JA  INTEGER ARRAYS OF THE SPARSE MATRIX REPRESENTATION 
C         A      D.P. ARRAY OF THE SPARSE MATRIX REPRESENTATION     
C         P      PERMUTATION VECTOR   
C         NEWIA  INTEGER WORK VECTOR OF LENGTH N
C         ISYM   SYMMETRIC/NONSYMMETRIC STORAGE SWITCH    
C         LEVEL  SWITCH CONTROLLING LEVEL OF OUTPUT       
C         NOUT OUTPUT UNIT NUMBER     
C         IER    OUTPUT ERROR FLAG (= IERR)     
C       
C                   IER =   0  NORMAL RETURN    
C                   IER = 301  NO ENTRY IN ITH ROW OF ORIGINAL      
C                              MATRIX. IF LEVEL IS GREATER THAN     
C                              0, I WILL BE PRINTED       
C                   IER = 302  THERE IS NO ENTRY IN THE ITH ROW     
C                              OF THE PERMUTED MATRIX     
C                   IER = 303  ERROR RETURN FROM QSORT IN 
C                              SORTING THE ITH ROW OF THE 
C                              PERMUTED MATRIX  
C ... IT IS ASSUMED THAT THE I-TH ENTRY OF THE PERMUTATION VECTOR   
C     P INDICATES THE ROW THE I-TH ROW GETS MAPPED INTO.  (I.E.     
C     IF ( P(I) = J ) ROW I GETS MAPPED INTO ROW J.)      
C       
C ... THE ARRAY NEWIA IS AN INTEGER WORK VECTOR OF LENGTH N WHICH   
C     KEEPS TRACK OF WHERE THE ROWS BEGIN IN THE PERMUTED STRUCTURE.
C       
C ... PERMAT IS CAPABLE OF PERMUTING BOTH THE SYMMETRIC AND NON-    
C     SYMMETRIC FORM OF IA, JA, AND A.  IF ( ISYM .EQ. 0 ) SYMMETRIC
C     FORM IS ASSUMED.      
C       
C ... TWO EXTERNAL MODULES ARE USED BY PERMAT.  THE FIRST IS INTEGER
C     FUNCTION BISRCH WHICH USES A BISECTION SEARCH ( ORDER LOG-BASE-2
C     OF N+1 ) THROUGH THE ARRAY IA TO FIND THE ROW INDEX OF AN ARBI- 
C     TRARY ENTRY EXTRACTED FROM THE ARRAY JA. THE SECOND IS SUBROUTINE 
C     QSORT WHICH PERFORMS A QUICK SORT TO PLACE THE ENTRIES IN     
C     THE PERMUTED ROWS IN COLUMN ORDER.
C       
C*********************************************************************
C       
      INTEGER NN,IA(1),JA(1),P(NN),NEWIA(NN),ISYM,IERR    
      DOUBLE PRECISION A(1) 
C       
C ... INTERNAL VARIABLES    
C       
      INTEGER BISRCH,I,IBGN,IEND,IP,IPP,J,JAJ,JP,IER,K,N,NELS,NEXT,NPL1 
C       
      DOUBLE PRECISION SAVE,TEMP      
C       
C*********************************************************************
C       
C ... PREPROCESSING PHASE   
C       
C ...... DETERMINE THE NUMBER OF NONZEROES IN THE ROWS OF THE PERMUTED
C        MATRIX AND STORE THAT IN NEWIA.  THEN SWEEP THRU NEWIA TO MAKE 
C        NEWIA(I) POINT TO THE BEGINNING OF EACH ROW IN THE PERMUTED
C        DATA STRUCTURE.  ALSO NEGATE ALL THE ENTRIES IN JA TO INDICATE 
C        THAT THOSE ENTRIES HAVE NOT BEEN MOVED YET.      
C       
      N = NN      
      IER = 0     
      NPL1 = N+1  
      NELS = IA(NPL1)-1     
      DO 10 I = 1,N 
         NEWIA(I) = 0       
   10 CONTINUE    
      DO 30 I = 1,N 
         IP = P(I)
         IBGN = IA(I)       
         IEND = IA(I+1)-1   
         IF (IBGN.GT.IEND) GO TO 90   
         DO 20 J = IBGN,IEND
            IPP = IP
            JAJ = JA(J)     
            JP = P(JAJ)     
            IF (ISYM.EQ.0.AND.IP.GT.JP) IPP = JP
            NEWIA(IPP) = NEWIA(IPP)+1 
            JA(J) = -JAJ    
   20    CONTINUE 
   30 CONTINUE    
      IBGN = 1    
      DO 40 I = 1,N 
         K = IBGN+NEWIA(I)  
         NEWIA(I) = IBGN    
         IBGN = K 
   40 CONTINUE    
C       
C ...... PREPROCESSING NOW FINISHED.  
C       
C ...... NOW PERMUTE JA AND A.  THIS PERMUTATION WILL PERFORM THE   
C        FOLLOWING STEPS    
C       
C           1.  FIND THE FIRST ENTRY IN JA NOT PERMUTED WHICH IS    
C               INDICATED BY AN NEGATIVE VALUE IN JA      
C           2.  COMPUTE WHICH ROW THE CURRENT ENTRY IS IN.  THIS    
C               IS COMPUTED BY A BISECTION SEARCH THRU THE ARRAY    
C               IA. 
C           3.  USING THE PERMUTATION ARRAY P AND THE ARRAY NEWIA   
C               COMPUTE WHERE THE CURRENT ENTRY IS TO BE PLACED.    
C           4.  THEN PICK UP THE ENTRY WHERE THE CURRENT ENTRY WILL 
C               GO.  PUT THE CURRENT ENTRY IN PLACE.  THEN MAKE THE 
C               DISPLACED ENTRY THE CURRENT ENTRY AND LOOP TO STEP 2. 
C           5.  THIS PROCESS WILL END WHEN THE NEXT ENTRY HAS ALREADY 
C               BEEN MOVED.  THEN LOOP TO STEP 1. 
C       
      DO 70 J = 1,NELS      
         IF (JA(J).GT.0) GO TO 70     
         JAJ = -JA(J)       
         SAVE = A(J)
         NEXT = J 
         JA(J) = JAJ
C       
   50    JP = P(JAJ)
         I = BISRCH(NPL1,IA,NEXT)     
         IP = P(I)
         IPP = IP 
         IF (ISYM.NE.0.OR.IP.LE.JP) GO TO 60    
         IPP = JP 
         JP = IP  
   60    NEXT = NEWIA(IPP)  
C       
         TEMP = SAVE
         SAVE = A(NEXT)     
         A(NEXT) = TEMP     
C       
         JAJ = -JA(NEXT)    
         JA(NEXT) = JP      
         NEWIA(IPP) = NEWIA(IPP)+1    
         IF (JAJ.GT.0) GO TO 50       
C       
   70 CONTINUE    
C       
C ...... THE MATRIX IS NOW PERMUTED BUT THE ROWS MAY NOT BE IN      
C        ORDER.  THE REMAINDER OF THIS SUBROUTINE PERFORMS
C        A QUICK SORT ON EACH ROW TO SORT THE ENTRIES IN  
C        COLUMN ORDER.  THE IA ARRAY IS ALSO CORRECTED FROM 
C        INFORMATION STORED IN THE NEWIA ARRAY.  NEWIA(I) NOW       
C        POINTS TO THE FIRST ENTRY OF ROW I+1.  
C       
      IA(1) = 1   
      DO 80 I = 1,N 
         IA(I+1) = NEWIA(I) 
         K = IA(I+1)-IA(I)  
         IF (K.EQ.1) GO TO 80 
         IF (K.LT.1) GO TO 110
C       
         IBGN = IA(I)       
         CALL QSORT (K,JA(IBGN),A(IBGN),IER)    
         IF (IER.NE.0) GO TO 130      
C       
   80 CONTINUE    
C       
C ...... END OF MATRIX PERMUTATION    
C       
      GO TO 150   
C       
C ... ERROR TRAPS 
C       
C ...... NO ENTRY IN ROW I IN THE ORIGINAL SYSTEM 
C       
   90 IER = 301   
      IF (LEVEL.GE.0) WRITE (NOUT,100) I
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE PERMAT  '/' ','    NO ENTRY IN ROW ',I10
     *   ,' OF ORIGINAL MATRIX ')     
      GO TO 150   
C       
C ...... NO ENTRY IN ROW I IN THE PERMUTED SYSTEM 
C       
  110 IER = 302   
      IF (LEVEL.GE.0) WRITE (NOUT,120) I
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE PRBNDX  '/' ','    NO ENTRY IN ROW ',I10
     *   ,' OF PERMUTED MATRIX ')     
      GO TO 150   
C       
C ...... ERROR RETURN FROM SUBROUTINE QSORT     
C       
  130 IER = 303   
      IF (LEVEL.GE.0) WRITE (NOUT,140) I
  140 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE QSORT   '/' ',  
     *   '    ERROR IN SORTING PERMUTED ROW ',I12/' ',    
     *   '    CALLED FROM ITPACK ROUTINE PRBNDX   ')      
C       
  150 CONTINUE    
      IERR = IER  
      RETURN      
      END 
      SUBROUTINE PERR (NN,IA,JA,A,RHS,U,W,DIGTT1,DIGTT2,IDGTTS)   
C       
C     PERR COMPUTES THE RESIDUAL, R = RHS - A*U.  THE USER
C     ALSO HAS THE OPTION OF PRINTING THE RESIDUAL AND/OR THE       
C     UNKNOWN VECTOR DEPENDING ON IDGTS.
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C          U      LATEST ESTIMATE OF SOLUTION   
C          W      WORKSPACE VECTOR    
C          DIGIT1 OUTPUT: MEASURE OF ACCURACY OF STOPPING TEST (= DIGTT1
C          DIGIT2 OUTPUT: MEASURE OF ACCURACY OF SOLUTION (= DIGTT2)
C          IDGTS   PARAMETER CONTROLING LEVEL OF OUTPUT (= IDGTTS)  
C                    IF IDGTS < 1 OR IDGTS > 4, THEN NO OUTPUT.     
C                            = 1, THEN NUMBER OF DIGITS IS PRINTED, PRO-
C                                 VIDED LEVEL .GE. 1      
C                            = 2, THEN SOLUTION VECTOR IS PRINTED, PRO- 
C                                 VIDED LEVEL .GE. 1      
C                            = 3, THEN RESIDUAL VECTOR IS PRINTED, PRO- 
C                                 VIDED LEVEL .GE. 1      
C                            = 4, THEN BOTH VECTORS ARE PRINTED, PRO- 
C                                 VIDED LEVEL .GE. 1      
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NN,IDGTTS   
      DOUBLE PRECISION A(1),RHS(NN),U(NN),W(NN),DIGTT1,DIGTT2       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IDGTS,N       
      DOUBLE PRECISION BNRM,DIGIT1,DIGIT2,RNRM,TEMP       
C       
C ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
C       
      DOUBLE PRECISION DDOT 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      IDGTS = IDGTTS
      DIGIT1 = 0.D0 
      DIGIT2 = 0.D0 
      IF (N.LE.0) GO TO 40  
C       
      DIGIT1 = -DLOG10(DABS(DRELPR))  
      IF (STPTST.GT.0.D0) DIGIT1 = -DLOG10(DABS(STPTST))  
      BNRM = DDOT(N,RHS,1,RHS,1)      
      IF (BNRM.EQ.0.D0) GO TO 10      
      CALL PMULT (N,IA,JA,A,U,W)      
      CALL WEVMW (N,RHS,W)  
      RNRM = DDOT(N,W,1,W,1)
      TEMP = RNRM/BNRM      
      IF (TEMP.EQ.0.D0) GO TO 10      
      DIGIT2 = -DLOG10(DABS(TEMP))/2.D0 
      GO TO 20    
C       
   10 DIGIT2 = -DLOG10(DABS(DRELPR))  
C       
   20 IF ((IDGTS.LT.1).OR.(LEVEL.LE.0)) GO TO 40
      WRITE (NOUT,30) DIGIT1,DIGIT2   
   30 FORMAT (/6X,'APPROX. NO. OF DIGITS (EST. REL. ERROR) =',F5.1,2X,
     *   '(DIGIT1)'/3X,'APPROX. NO. OF DIGITS (EST. REL. RESIDUAL) =',
     *   F5.1,2X,'(DIGIT2)')
C       
      IF (IDGTS.LE.1.OR.IDGTS.GT.4) GO TO 40    
      IF (IDGTS.NE.3) CALL VOUT (N,U,2,NOUT)    
      IF (IDGTS.GE.3) CALL VOUT (N,W,1,NOUT)    
C       
   40 CONTINUE    
      DIGTT1 = DIGIT1       
      DIGTT2 = DIGIT2       
      RETURN      
      END 
      SUBROUTINE PERVEC (N,V,P)       
C       
C     THIS SUBROUTINE PERMUTES A D.P. VECTOR AS DICTATED BY THE     
C     PERMUTATION VECTOR, P.  IF P(I) = J, THEN V(J) GETS V(I).     
C       
C ... PARAMETER LIST:       
C       
C          V      D.P. VECTOR OF LENGTH N       
C          P     INTEGER PERMUTATION VECTOR     
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,P(N)
      DOUBLE PRECISION V(N) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER II,NEXT,NOW   
      DOUBLE PRECISION SAVE,TEMP      
C       
      IF (N.LE.0) RETURN    
C       
      DO 20 II = 1,N
         IF (P(II).LT.0) GO TO 20     
C       
         NEXT = P(II)       
         SAVE = V(II)       
C       
   10    CONTINUE 
         IF (P(NEXT).LT.0) GO TO 20   
         TEMP = SAVE
         SAVE = V(NEXT)     
         V(NEXT) = TEMP     
C       
         NOW = NEXT 
         NEXT = P(NOW)      
         P(NOW) = -NEXT     
         GO TO 10 
C       
   20 CONTINUE    
C       
      DO 30 II = 1,N
         P(II) = -P(II)     
   30 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE PFSOR (NN,IA,JA,A,U,RHS)       
C       
C         THIS SUBROUTINE COMPUTES A FORWARD SOR SWEEP.   
C       
C ... PARAMETER LIST:       
C       
C         N       ORDER OF SYSTEM (= NN)
C          OMEGA  RELAXATION FACTOR   
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          U      LATEST ESTIMATE OF SOLUTION   
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),RHS(NN)       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,JAJJ,JJ,N  
      DOUBLE PRECISION OMM1,SUM,UI    
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      OMM1 = OMEGA-1.D0     
      IF (ISYM.EQ.0) GO TO 40 
C       
C     *********** NON - SYMMETRIC SECTION *********************     
C       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    UI = OMEGA*SUM-OMM1*U(II)    
         U(II) = UI 
   30 CONTINUE    
      RETURN      
C       
C     ************* SYMMETRIC SECTION *************************     
C       
   40 DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   50    CONTINUE 
   60    UI = OMEGA*SUM-OMM1*U(II)    
         U(II) = UI 
         IF (IBGN.GT.IEND) GO TO 80   
         DO 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UI      
   70    CONTINUE 
   80 CONTINUE    
      RETURN      
C       
      END 
      SUBROUTINE PFSOR1 (NN,IA,JA,A,U,RHS)      
C       
C         THIS SUBROUTINE COMPUTES A FORWARD SOR SWEEP ON U AND     
C         COMPUTES THE NORM OF THE PSEUDO-RESIDUAL VECTOR.
C       
C ... PARAMETER LIST:       
C       
C          N      ORDER OF SYSTEM (= NN)
C          OMEGA  RELAXATION FACTOR   
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          U      LATEST ESTIMATE OF SOLUTION   
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),RHS(NN)       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,JAJJ,JJ,N  
      DOUBLE PRECISION OMM1,SUM,SUMD,UI 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      OMM1 = OMEGA-1.D0     
      SUMD = 0.D0 
      IF (ISYM.EQ.0) GO TO 40 
C       
C     **************** NON - SYMMETRIC SECTION ******************   
C       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    CONTINUE 
         UI = OMEGA*SUM-OMM1*U(II)    
         SUMD = SUMD+(UI-U(II))**2    
         U(II) = UI 
   30 CONTINUE    
      GO TO 90    
C       
C     *************** SYMMETRIC SECTION ************************    
C       
   40 DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   50    CONTINUE 
   60    CONTINUE 
         UI = OMEGA*SUM-OMM1*U(II)    
         SUMD = SUMD+(UI-U(II))**2    
         U(II) = UI 
         IF (IBGN.GT.IEND) GO TO 80   
         DO 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UI      
   70    CONTINUE 
   80 CONTINUE    
C       
   90 DELNNM = DSQRT(SUMD)  
      RETURN      
C       
      END 
      SUBROUTINE PJAC (NN,IA,JA,A,U,RHS)
C       
C     ... THIS SUBROUTINE PERFORMS ONE JACOBI ITERATION.  
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          U      ESTIMATE OF SOLUTION OF A MATRIX PROBLEM
C          RHS    ON INPUT: CONTAINS THE RIGHT HAND SIDE OF 
C                    A MATRIX PROBLEM 
C                 ON OUTPUT: CONTAINS A*U + RHS 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),RHS(NN)       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,JAJJ,JJ,N  
      DOUBLE PRECISION RHSII,UII      
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      IF (ISYM.EQ.0) GO TO 30 
C       
C     *************** NON - SYMMETRIC SECTION ****************      
C       
      DO 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 20   
         RHSII = RHS(II)    
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHSII = RHSII-A(JJ)*U(JAJJ) 
   10    CONTINUE 
         RHS(II) = RHSII    
   20 CONTINUE    
      RETURN      
C       
C     ************** SYMMETRIC SECTION **********************       
C       
   30 DO 50 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 50   
         RHSII = RHS(II)    
         UII = U(II)
         DO 40 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHSII = RHSII-A(JJ)*U(JAJJ) 
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UII     
   40    CONTINUE 
         RHS(II) = RHSII    
   50 CONTINUE    
      RETURN      
C       
      END 
      SUBROUTINE PMULT (NN,IA,JA,A,U,W) 
C       
C     ... THIS SUBROUTINE PERFORMS ONE MATRIX-VECTOR MULTIPLICATION.
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          U      LATEST ESTIMATE OF SOLUTION   
C          W      ON RETURN W CONTAINS A*U      
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),W(NN) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,JJ,N       
      DOUBLE PRECISION SUM,UII,WII    
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      IF (N.LE.0) RETURN    
      IF (ISYM.EQ.0) GO TO 40 
C       
C     *************** NON - SYMMETRIC SECTION **********************
C       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = 0.0D0
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM+A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    W(II) = SUM
   30 CONTINUE    
      RETURN      
C       
C     ***************** SYMMETRIC SECTION **************************
C       
   40 CALL VFILL (N,W,0.D0) 
      DO 70 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         UII = U(II)
         WII = W(II)
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            WII = WII+A(JJ)*U(JAJJ)   
            W(JAJJ) = W(JAJJ)+A(JJ)*UII 
   50    CONTINUE 
   60    W(II) = WII
   70 CONTINUE    
      RETURN      
C       
      END 
      SUBROUTINE PRBNDX (NN,NBLACK,IA,JA,P,IP,LEVEL,NOUT,IER)       
C       
C**************************************************************     
C       
C     THIS SUBROUTINE COMPUTES THE RED-BLACK PERMUTATION  
C     VECTORS P ( AND ITS INVERSE IP ) IF POSSIBLE.       
C       
C     THE ALGORITHM IS TO MARK THE FIRST NODE AS RED (ARBITRARY).   
C     ALL OF ITS ADJACENT NODES ARE MARKED BLACK AND PLACED IN      
C     A STACK.  THE REMAINDER OF THE CODE PULLS THE FIRST NODE      
C     OFF THE TOP OF THE STACK AND TRIES TO TYPE ITS ADJACENT NODES.
C     THE TYPING OF THE ADJACENT POINT IS A FIVE WAY CASE STATEMENT 
C     WHICH IS WELL COMMENTED BELOW (SEE DO LOOP 100).    
C       
C     THE ARRAY P IS USED BOTH TO KEEP TRACK OF THE COLOR OF A NODE 
C     (RED NODE IS POSITIVE, BLACK IS NEGATIVE) BUT ALSO THE FATHER 
C     NODE THAT CAUSED THE COLOR MARKING OF THAT POINT.  SINCE      
C     COMPLETE INFORMATION ON THE ADJACENCY STRUCTURE IS HARD TO COME 
C     BY THIS FORMS A LINK TO ENABLE THE COLOR CHANGE OF A PARTIAL  
C     TREE WHEN A RECOVERABLE COLOR CONFLICT OCCURS.      
C       
C     THE ARRAY IP IS USED AS A STACK TO POINT TO THE SET OF NODES  
C     LEFT TO BE TYPED THAT ARE KNOWN TO BE ADJACENT TO THE CURRENT 
C     FATHER NODE.
C       
C*********************************************************************
C       
C     INPUT PARAMETERS      
C       
C        N      NUMBER OF NODES.  (INTEGER, SCALAR) (= NN)
C       
C        IA,JA  ADJACENCY STRUCTURE ARRAYS.  CAN BE EITHER THE      
C               SYMMETRIC OR NONSYMMETRIC FORM.  IT IS ASSUMED      
C               THAT FOR EVERY ROW WHERE ONLY ONE ELEMENT IS
C               STORED THAT ELEMENT CORRESPONDS TO THE DIAGONAL     
C               ENTRY.  THE DIAGONAL DOES NOT HAVE TO BE THE FIRST  
C               ENTRY STORED.  (INTEGER, ARRAYS)
C        LEVEL  SWITCH FOR PRINTING   
C        NOUT OUTPUT TAPE NUMBER      
C       
C     OUTPUT PARAMETERS     
C       
C        NBLACK NUMBER OF BLACK NODES.  NUMBER OF RED NODES IS      
C               N - NBLACK.  (INTEGER, SCALAR)  
C       
C        P, IP  PERMUTATION AND INVERSE PERMUTATION VECTORS.
C               (INTEGER, ARRAYS EACH OF LENGTH N)
C       
C        IER    ERROR FLAG. (INTEGER, SCALAR)   
C       
C               IER = 0, NORMAL RETURN.  INDEXING PERFORMED 
C                        SUCCESSFULLY 
C               IER =201, RED-BLACK INDEXING NOT POSSIBLE.
C       
C******************************************************************** 
C       
      INTEGER NN,NBLACK,IA(1),JA(1),P(NN),IP(NN),IER      
C       
      INTEGER FIRST,NEXT,LAST,I,OLD,YOUNG,IBGN,IEND,J,K,CURTYP,NXTTYP,
     *   TYPE,NRED,N
C       
C-----------------------------------------------------------------------
C       
      N = NN      
      IER = 0     
C       
C        IF ( N .LE. 0 ) GO TO 8000   
C       
      DO 10 I = 1,N 
         P(I) = 0 
         IP(I) = 0
   10 CONTINUE    
C       
C ... HANDLE THE FIRST SET OF POINTS UNTIL SOME ADJACENT POINTS     
C ... ARE FOUND   
C       
      FIRST = 1   
C       
   20 P(FIRST) = FIRST      
      IF (IA(FIRST+1)-IA(FIRST).GT.1) GO TO 40  
C       
C ... SEARCH FOR NEXT ENTRY THAT HAS NOT BEEN MARKED      
C       
      IF (FIRST.EQ.N) GO TO 130       
      IBGN = FIRST+1
      DO 30 I = IBGN,N      
         IF (P(I).NE.0) GO TO 30      
         FIRST = I
         GO TO 20 
   30 CONTINUE    
      GO TO 130   
C       
C ... FIRST SET OF ADJACENT POINTS FOUND
C       
   40 NEXT = 1    
      LAST = 1    
      IP(1) = FIRST 
C       
C ... LOOP OVER LABELED POINTS INDICATED IN THE STACK STORED IN     
C ... THE ARRAY IP
C       
   50 K = IP(NEXT)
      CURTYP = P(K) 
      NXTTYP = -CURTYP      
      IBGN = IA(K)
      IEND = IA(K+1)-1      
      IF (IBGN.GT.IEND) GO TO 110     
      DO 100 I = IBGN,IEND  
         J = JA(I)
         TYPE = P(J)
         IF (J.EQ.K) GO TO 100
C       
C================================================================== 
C       
C     THE FOLLOWING IS A FIVE WAY CASE STATEMENT DEALING WITH THE   
C     LABELING OF THE ADJACENT NODE.  
C       
C ... CASE I.  IF THE ADJACENT NODE HAS ALREADY BEEN LABELED WITH   
C              LABEL EQUAL TO NXTTYP, THEN SKIP TO THE NEXT ADJACENT
C              NODE.
C       
         IF (TYPE.EQ.NXTTYP) GO TO 100
C       
C ... CASE II.  IF THE ADJACENT NODE HAS NOT BEEN LABELED YET LABEL 
C               IT WITH NXTTYP AND ENTER IT IN THE STACK  
C       
         IF (TYPE.NE.0) GO TO 60      
         LAST = LAST+1      
         IP(LAST) = J       
         P(J) = NXTTYP      
         GO TO 100
C       
C ... CASE III.  IF THE ADJACENT NODE HAS ALREADY BEEN LABELED WITH 
C                OPPOSITE COLOR AND THE SAME FATHER SEED, THEN THERE
C                IS AN IRRECOVERABLE COLOR CONFLICT.      
C       
   60    IF (TYPE.EQ.CURTYP) GO TO 160
C       
C ... CASE IV.  IF THE ADJACENT NODE HAS THE RIGHT COLOR AND A DIFFERENT
C               FATHER NODE, THEN CHANGE ALL NODES OF THE YOUNGEST FATHE
C               NODE TO POINT TO THE OLDEST FATHER SEED AND RETAIN THE
C               SAME COLORS.
C       
         IF (TYPE*NXTTYP.LT.1) GO TO 80 
         OLD = MIN0(IABS(TYPE),IABS(NXTTYP))    
         YOUNG = MAX0(IABS(TYPE),IABS(NXTTYP))  
         DO 70 J = YOUNG,N  
            IF (IABS(P(J)).EQ.YOUNG) P(J) = ISIGN(OLD,P(J)) 
   70    CONTINUE 
         CURTYP = P(K)      
         NXTTYP = -CURTYP   
         GO TO 100
C       
C ... CASE V.  IF THE ADJACENT NODE HAS THE WRONG COLOR AND A DIFFERENT 
C              FATHER NODE, THEN CHANGE ALL NODES OF THE YOUNGEST FATHER
C              NODE TO POINT TO THE OLDEST FATHER NODE ALONG WITH   
C              CHANGING THEIR COLORS.  SINCE UNTIL THIS TIME THE    
C              YOUNGEST FATHER NODE TREE HAS BEEN INDEPENDENT NO OTHER
C              COLOR CONFLICTS WILL ARISE FROM THIS CHANGE. 
C       
   80    OLD = MIN0(IABS(TYPE),IABS(NXTTYP))    
         YOUNG = MAX0(IABS(TYPE),IABS(NXTTYP))  
         DO 90 J = YOUNG,N  
            IF (IABS(P(J)).EQ.YOUNG) P(J) = ISIGN(OLD,-P(J))
   90    CONTINUE 
         CURTYP = P(K)      
         NXTTYP = -CURTYP   
C       
C ... END OF CASE STATEMENT 
C       
C================================================================== 
C       
  100 CONTINUE    
C       
C ... ADVANCE TO NEXT NODE IN THE STACK 
C       
  110 NEXT = NEXT+1 
      IF (NEXT.LE.LAST) GO TO 50      
C       
C ... ALL NODES IN THE STACK HAVE BEEN REMOVED  
C       
C ... CHECK FOR NODES NOT LABELED.  IF ANY ARE FOUND      
C ... START THE LABELING PROCESS AGAIN AT THE FIRST       
C ... NODE FOUND THAT IS NOT LABELED. 
C       
      IBGN = FIRST+1
      DO 120 I = IBGN,N     
         IF (P(I).NE.0) GO TO 120     
         FIRST = I
         GO TO 20 
  120 CONTINUE    
C       
C===================================================================
C       
C ... ALL NODES ARE NOW TYPED EITHER RED OR BLACK 
C       
C ... GENERATE PERMUTATION VECTORS    
C       
  130 NRED = 0    
      NBLACK = 0  
      DO 150 I = 1,N
         IF (P(I).LT.0) GO TO 140     
C       
C       RED POINT 
C       
         NRED = NRED+1      
         IP(NRED) = I       
         P(I) = NRED
         GO TO 150
C       
C     BLACK POINT 
C       
  140    NBLACK = NBLACK+1  
         J = N-NBLACK+1     
         IP(J) = I
         P(I) = J 
C       
  150 CONTINUE    
C       
C ... SUCCESSFUL RED-BLACK ORDERING COMPLETED   
C       
      GO TO 180   
C       
C ........ ERROR TRAPS      
C       
C ...... N .LE. 0 
C       
C8000    IER = 200
C        GO TO 9000 
C       
C ...... TYPE CONFLICT      
C       
  160 IER = 201   
      IF (LEVEL.GE.0) WRITE (NOUT,170)
  170 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE PRBNDX  '/' ',  
     *   '    RED-BLACK INDEXING NOT POSSIBLE') 
C       
C ... RETURN      
C       
  180 CONTINUE    
      RETURN      
      END 
      SUBROUTINE PRSBLK (NNB,NNR,IA,JA,A,UR,VB) 
C       
C ... COMPUTE A BLACK-RS SWEEP ON A RED VECTOR INTO A BLACK VECTOR  
C       
C ... PARAMETER LIST:       
C       
C         NB      NUMBER OF BLACK POINTS (= NNB)
C         NR      NUMBER OF RED POINTS (= NNR)  
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          UR     ESTIMATE OF RED SOLUTION VECTOR 
C          VB     OUTPUT: PRESENT ESTIMATE OF BLACK SOLUTION
C                    VECTOR 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NNB,NNR     
      DOUBLE PRECISION A(1),UR(NNR),VB(NNB)     
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IBGN,IEND,INR,J,JAJ,NB,NR       
      DOUBLE PRECISION SUM,URI
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      NB = NNB    
      NR = NNR    
      IF (ISYM.EQ.0) GO TO 30 
C       
C     *************** NON - SYMMETRIC SECTION **********************
C       
      DO 20 I = 1,NB
         INR = I+NR 
         IBGN = IA(INR)     
         IEND = IA(INR+1)-1 
         SUM = VB(I)
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 J = IBGN,IEND
            JAJ = JA(J)     
            SUM = SUM-A(J)*UR(JAJ)    
   10    CONTINUE 
         VB(I) = SUM
   20 CONTINUE    
      RETURN      
C       
C     ***************** SYMMETRIC SECTION **************************
C       
   30 DO 50 I = 1,NR
         IBGN = IA(I)       
         IEND = IA(I+1)-1   
         IF (IBGN.GT.IEND) GO TO 50   
         URI = UR(I)
         DO 40 J = IBGN,IEND
            JAJ = JA(J)-NR  
            VB(JAJ) = VB(JAJ)-A(J)*URI
   40    CONTINUE 
   50 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE PRSRED (NNB,NNR,IA,JA,A,UB,VR) 
C       
C ... COMPUTES A RED-RS SWEEP ON A BLACK VECTOR INTO A RED VECTOR.  
C       
C ... PARAMETER LIST:       
C       
C         NB      NUMBER OF BLACK POINTS (= NNR)
C         NR      NUMBER OF RED POINTS (= NNB)  
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          UB     PRESENT ESTIMATE OF BLACK SOLUTION VECTOR 
C          VR     OUTPUT: PRESENT ESTIMATE OF RED SOLUTION VECTOR   
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NNB,NNR     
      DOUBLE PRECISION A(1),UB(NNB),VR(NNR)     
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,JAJJ,JJ,NB,NR
      DOUBLE PRECISION SUM  
C       
      NB = NNB    
      NR = NNR    
      DO 20 II = 1,NR       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 20   
         SUM = VR(II)       
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)-NR
            SUM = SUM-A(JJ)*UB(JAJJ)  
   10    CONTINUE 
         VR(II) = SUM       
   20 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE PSSOR1 (NN,IA,JA,A,U,RHS,FR,BR)
C       
C     ... COMPUTES COMPLETE SSOR SWEEP ON U.  U IS OVERWRITTEN      
C     ... WITH THE NEW ITERANT, FR AND BR WILL CONTAIN    
C     ... THE FORWARD AND BACKWARD RESIDUALS ON OUTPUT.   
C       
C ... PARAMETER LIST:       
C       
C         N       ORDER OF SYSTEM (= NN)
C          OMEGA  RELAXATION FACTOR   
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          U      ESTIMATE OF SOLUTION
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C          FR,BR  OUTPUT: FORWARD AND BACKWARD RESIDUALS RESPECTIVELY 
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),RHS(NN),FR(NN),BR(NN)   
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IBGN,IEND,II,JAJJ,JJ,N,NPL1     
      DOUBLE PRECISION OMM1,SUM,UII   
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      N = NN      
      NPL1 = N+1  
      OMM1 = OMEGA-1.D0     
      IF (ISYM.EQ.0) GO TO 40 
C       
C     *************** NON - SYMMETRIC SECTION **********************
C       
C     ... FORWARD SWEEP     
C       
      DO 30 II = 1,N
         BR(II) = U(II)     
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    UII = OMEGA*SUM-OMM1*U(II)   
         FR(II) = UII-U(II) 
         U(II) = UII
   30 CONTINUE    
      GO TO 90    
C       
C     ***************** SYMMETRIC SECTION **************************
C       
C     ... FORWARD SWEEP     
C       
   40 DO 80 II = 1,N
         BR(II) = U(II)     
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   50    CONTINUE 
   60    UII = OMEGA*SUM-OMM1*U(II)   
         FR(II) = UII-U(II) 
         U(II) = UII
         IF (IBGN.GT.IEND) GO TO 80   
         DO 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UII     
   70    CONTINUE 
   80 CONTINUE    
C       
C     ... BACKWARD SWEEP    
C       
   90 DO 120 I = 1,N
         II = NPL1-I
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         UII = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 110  
         DO 100 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            UII = UII-A(JJ)*U(JAJJ)   
  100    CONTINUE 
  110    U(II) = OMEGA*UII-OMM1*U(II) 
         BR(II) = U(II)-BR(II)
  120 CONTINUE    
C       
      RETURN      
C       
      END 
      SUBROUTINE PSTOP (N,U,DNRM,CCON,IFLAG,Q1) 
C       
C     THIS SUBROUTINE PERFORMS A TEST TO SEE IF THE ITERATIVE       
C     METHOD HAS CONVERGED TO A SOLUTION INSIDE THE ERROR 
C     TOLERANCE, ZETA.      
C       
C ... PARAMETER LIST:       
C       
C          N      ORDER OF SYSTEM     
C          U      PRESENT SOLUTION ESTIMATE     
C          DNRM   INNER PRODUCT OF PSEUDO-RESIDUALS AT PRECEDING    
C                    ITERATION
C          CON    STOPPING TEST PARAMETER (= CCON)
C          IFLAG  STOPPING TEST INTEGER FLAG    
C                    IFLAG = 0,  SOR ITERATION ZERO       
C                    IFLAG = 1,  NON-RS METHOD  
C                    IFLAG = 2,  RS METHOD      
C          Q1     STOPPING TEST LOGICAL FLAG    
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,IFLAG       
      DOUBLE PRECISION U(N),DNRM,CCON 
      LOGICAL Q1  
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      DOUBLE PRECISION CON,TL,TR,UOLD 
C       
C ... SPECIFICATIONS FOR ARGUMENT SUBROUTINES   
C       
      DOUBLE PRECISION DDOT 
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      CON = CCON  
      HALT = .FALSE.
C       
C     SPECIAL PROCEDURE FOR ZEROTH ITERATION    
C       
      IF (IN.GE.1) GO TO 10 
      Q1 = .FALSE.
      UDNM = 1.D0 
      STPTST = 1.D3 
      IF (IFLAG.LE.0) RETURN
C       
C ... TEST IF UDNM NEEDS TO BE RECOMPUTED       
C       
   10 CONTINUE    
      IF (Q1) GO TO 20      
      IF ((IN.GT.5).AND.(MOD(IN,5).NE.0)) GO TO 20
      UOLD = UDNM 
      UDNM = DDOT(N,U,1,U,1)
      IF (UDNM.EQ.0.D0) UDNM = 1.D0   
      IF ((IN.GT.5).AND.(DABS(UDNM-UOLD).LE.UDNM*ZETA)) Q1 = .TRUE. 
C       
C ... COMPUTE STOPPING TEST 
C       
   20 TR = DSQRT(UDNM)      
      TL = 1.D0   
      IF (CON.EQ.1.D0) GO TO 40       
      IF (IFLAG.EQ.2) GO TO 30
      TL = DSQRT(DNRM)      
      TR = TR*(1.D0-CON)    
      GO TO 40    
   30 TL = DSQRT(2.D0*DNRM) 
      TR = TR*(1.D0-CON*CON)
   40 STPTST = TL/TR
      IF (TL.GE.TR*ZETA) RETURN       
      HALT = .TRUE. 
C       
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION PVTBV (N,IA,JA,A,V)       
C       
C     THIS FUNCTION COMPUTES  (V**T)*A*V.       
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX 
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          V      D.P. VECTOR OF LENGTH N       
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),N 
      DOUBLE PRECISION A(1),V(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,JAJJ,JJ    
      DOUBLE PRECISION SUM,SUMR       
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
      PVTBV = 0.D0
      SUM = 0.D0  
      DO 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 20   
         SUMR = 0.D0
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUMR = SUMR-A(JJ)*V(JAJJ) 
   10    CONTINUE 
         SUM = SUM+V(II)*SUMR 
   20 CONTINUE    
C       
      IF (ISYM.EQ.0) SUM = 2.D0*SUM   
      PVTBV = SUM 
C       
      RETURN      
      END 
      SUBROUTINE QSORT (NN,KEY,DATA,ERROR)      
C       
C     ==================================================================
C       
C     Q U I C K S O R T     
C       
C         IN THE STYLE OF THE CACM PAPER BY BOB SEDGEWICK, OCTOBER 1978 
C       
C     INPUT:      
C         N    -- NUMBER OF ELEMENTS TO BE SORTED (= NN)  
C         KEY  -- AN ARRAY OF LENGTH  N  CONTAINING THE VALUES      
C                 WHICH ARE TO BE SORTED
C         DATA -- A SECOND ARRAY OF LENGTH  N  CONTAINING DATA      
C                 ASSOCIATED WITH THE INDIVIDUAL KEYS.    
C       
C     OUTPUT:     
C         KEY  -- WILL BE ARRANGED SO THAT VALUES ARE IN INCREASING 
C                 ORDER     
C         DATA -- REARRANGED TO CORRESPOND TO REARRANGED KEYS       
C         ERROR -- WILL BE ZERO UNLESS YOUR INPUT FILE WAS OF TRULY 
C                  ENORMOUS LENGTH, IN WHICH CASE IT WILL BE EQUAL TO 1.
C       
C     ==================================================================
C       
      INTEGER NN,ERROR,KEY(NN)
      DOUBLE PRECISION DATA(NN)       
C       
C     ------------------------
C       
      INTEGER TOP,LEFT,RIGHT,I,J,TINY,V,K,IP1,JM1,LLEN,RLEN,N       
      LOGICAL DONE
      DOUBLE PRECISION D    
      INTEGER STKLEN,STACK(30)
C       
      DATA TINY,STKLEN / 9,30 /       
C       
C     -----------------------------------       
C       
C     ... PROGRAM IS A DIRECT TRANSLATION INTO FORTRAN OF SEDGEWICK^S 
C         PROGRAM 2, WHICH IS NON-RECURSIVE, IGNORES FILES OF LENGTH
C         LESS THAN 'TINY' DURING PARTITIONING, AND USES MEDIAN OF THREE
C         PARTITIONING.     
C       
      N = NN      
      IF (N.EQ.1) RETURN    
      IF (N.LE.0) GO TO 240 
C       
      ERROR = 0   
      TOP = 1     
      LEFT = 1    
      RIGHT = N   
      DONE = (N.LE.TINY)    
C       
      IF (DONE) GO TO 150   
      CALL IVFILL (STKLEN,STACK,0)    
C       
C     ===========================================================   
C     QUICKSORT -- PARTITION THE FILE UNTIL NO SUBFILE REMAINS OF   
C     LENGTH GREATER THAN 'TINY'      
C     ===========================================================   
C       
C     ... WHILE NOT DONE DO ...       
C       
   10 IF (DONE) GO TO 150   
C       
C         ... FIND MEDIAN OF LEFT, RIGHT AND MIDDLE ELEMENTS OF CURRENT 
C             SUBFILE, WHICH IS  KEY(LEFT), ..., KEY(RIGHT) 
C       
      LFRH2 = (LEFT+RIGHT)/2
      K = KEY(LFRH2)
      D = DATA(LFRH2)       
      KEY(LFRH2) = KEY(LEFT)
      DATA(LFRH2) = DATA(LEFT)
      KEY(LEFT) = K 
      DATA(LEFT) = D
C       
      IF (KEY(LEFT+1).LE.KEY(RIGHT)) GO TO 20   
      K = KEY(LEFT+1)       
      D = DATA(LEFT+1)      
      KEY(LEFT+1) = KEY(RIGHT)
      DATA(LEFT+1) = DATA(RIGHT)      
      KEY(RIGHT) = K
      DATA(RIGHT) = D       
C       
   20 IF (KEY(LEFT).LE.KEY(RIGHT)) GO TO 30     
      K = KEY(LEFT) 
      D = DATA(LEFT)
      KEY(LEFT) = KEY(RIGHT)
      DATA(LEFT) = DATA(RIGHT)
      KEY(RIGHT) = K
      DATA(RIGHT) = D       
C       
   30 IF (KEY(LEFT+1).LE.KEY(LEFT)) GO TO 40    
      K = KEY(LEFT+1)       
      D = DATA(LEFT+1)      
      KEY(LEFT+1) = KEY(LEFT) 
      DATA(LEFT+1) = DATA(LEFT)       
      KEY(LEFT) = K 
      DATA(LEFT) = D
C       
   40 V = KEY(LEFT) 
C       
C         ... V IS NOW THE MEDIAN VALUE OF THE THREE KEYS.  NOW MOVE
C             FROM THE LEFT AND RIGHT ENDS SIMULTANEOUSLY, EXCHANGING 
C             KEYS AND DATA UNTIL ALL KEYS LESS THAN  V  ARE PACKED TO
C             THE LEFT, ALL KEYS LARGER THAN  V  ARE PACKED TO THE  
C             RIGHT.
C       
      I = LEFT+1  
      J = RIGHT   
C       
C         LOOP    
C             REPEAT I = I+1 UNTIL KEY(I) >= V; 
C             REPEAT J = J-1 UNTIL KEY(J) <= V; 
C         EXIT IF J < I;    
C             << EXCHANGE KEYS I AND J >>       
C         END     
C       
   50 CONTINUE    
   60 I = I+1     
      IF (KEY(I).LT.V) GO TO 60       
C       
   70 J = J-1     
      IF (KEY(J).GT.V) GO TO 70       
C       
      IF (J.LT.I) GO TO 80  
      K = KEY(I)  
      D = DATA(I) 
      KEY(I) = KEY(J)       
      DATA(I) = DATA(J)     
      KEY(J) = K  
      DATA(J) = D 
      GO TO 50    
C       
   80 K = KEY(LEFT) 
      D = DATA(LEFT)
      KEY(LEFT) = KEY(J)    
      DATA(LEFT) = DATA(J)  
      KEY(J) = K  
      DATA(J) = D 
C       
C         ... WE HAVE NOW PARTITIONED THE FILE INTO TWO SUBFILES,   
C             ONE IS (LEFT ... J-1)  AND THE OTHER IS (I...RIGHT).  
C             PROCESS THE SMALLER NEXT.  STACK THE LARGER ONE.      
C       
      LLEN = J-LEFT 
      RLEN = RIGHT-I+1      
      IF (MAX0(LLEN,RLEN).GT.TINY) GO TO 100    
C       
C             ... BOTH SUBFILES ARE TINY, SO UNSTACK NEXT LARGER FILE 
C       
      IF (TOP.EQ.1) GO TO 90
      TOP = TOP-2 
      LEFT = STACK(TOP)     
      RIGHT = STACK(TOP+1)  
      GO TO 10    
C       
   90 DONE = .TRUE. 
C       
      GO TO 10    
C       
C             ... ELSE ONE OR BOTH SUBFILES ARE LARGE     
C       
  100 IF (MIN0(LLEN,RLEN).GT.TINY) GO TO 120    
C       
C             ... ONE SUBFILE IS SMALL, ONE LARGE.  IGNORE THE SMALL ONE
C       
      IF (LLEN.GT.RLEN) GO TO 110     
      LEFT = I    
      GO TO 10    
C       
  110 RIGHT = J-1 
C       
      GO TO 10    
C       
C         ... ELSE BOTH ARE LARGER THAN TINY.  ONE MUST BE STACKED. 
C       
  120 IF (TOP.GE.STKLEN) GO TO 240    
      IF (LLEN.GT.RLEN) GO TO 130     
      STACK(TOP) = I
      STACK(TOP+1) = RIGHT  
      RIGHT = J-1 
      GO TO 140   
C       
  130 STACK(TOP) = LEFT     
      STACK(TOP+1) = J-1    
      LEFT = I    
C       
  140 TOP = TOP+2 
C       
      GO TO 10    
C       
C     ------------------------------------------------------------  
C     INSERTION SORT THE ENTIRE FILE, WHICH CONSISTS OF A LIST      
C     OF 'TINY' SUBFILES, LOCALLY OUT OF ORDER, GLOBALLY IN ORDER.  
C     ------------------------------------------------------------  
C       
C     ... FIRST, FIND LARGEST ELEMENT IN 'KEY'  
C       
  150 I = N-1     
      LEFT = MAX0(0,N-TINY) 
      K = KEY(N)  
      J = N       
C       
  160 IF (I.LE.LEFT) GO TO 180
      IF (KEY(I).LE.K) GO TO 170      
      K = KEY(I)  
      J = I       
C       
  170 I = I-1     
      GO TO 160   
C       
  180 IF (J.EQ.N) GO TO 190 
C       
C     ... LARGEST ELEMENT WILL BE IN  KEY(N)    
C       
      KEY(J) = KEY(N)       
      KEY(N) = K  
      D = DATA(N) 
      DATA(N) = DATA(J)     
      DATA(J) = D 
C       
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...       
C       
  190 I = N-1     
      IP1 = N     
C       
  200 IF (KEY(I).LE.KEY(IP1)) GO TO 220 
C       
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE 
C       
      K = KEY(I)  
      D = DATA(I) 
      J = IP1     
      JM1 = I     
C       
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR K FOUND'      
C       
  210 KEY(JM1) = KEY(J)     
      DATA(JM1) = DATA(J)   
      JM1 = J     
      J = J+1     
      IF (KEY(J).LT.K) GO TO 210      
C       
      KEY(JM1) = K
      DATA(JM1) = D 
C       
  220 IP1 = I     
      I = I-1     
      IF (I.GT.0) GO TO 200 
C       
  230 RETURN      
C       
  240 ERROR = 1   
      GO TO 230   
C       
      END 
      SUBROUTINE SBAGN (N,NZ,IA,JA,A,IWORK,LEVELL,NOUTT,IERR)       
C       
C ... THE ROUTINES SBINI, SBSIJ, AND SBEND CREATE A SPARSE
C     MATRIX STRUCTURE BY MEANS OF A LINKED LIST WHICH IS 
C     DESTROYED BY SBEND. SBAGN CREATES A NEW LINKED LIST 
C     SO THAT ELEMENTS MAY BE ADDED TO THE MATRIX AFTER SBEND       
C     HAS BEEN CALLED. SBAGN SHOULD BE CALLED WITH THE APPRO-       
C     PRIATE PARAMETERS, AND THEN SBSIJ AND SBEND CAN BE CALLED     
C     TO ADD THE ELEMENTS AND COMPLETE THE SPARSE MATRIX STRUC-     
C     TURE.       
C       
C ... PARAMETER LIST:       
C       
C           N       ORDER OF THE SYSTEM 
C           NZ      MAXIMUM NUMBER OF NON-ZERO ELEMENTS   
C                   IN THE SYSTEM     
C           IA, JA  INTEGER ARRAYS OF THE SPARSE
C                   MATRIX STRUCTURE  
C           A       D.P. ARRAY OF THE SPARSE MATRIX       
C                   STRUCTURE 
C           IWORK   WORK ARRAY OF DIMENSION NZ  
C           LEVEL   OUTPUT LEVEL CONTROL (= LEVELL)       
C           NOUT  OUTPUT FILE NUMBER (= NOUTT)  
C           IER     ERROR FLAG (= IERR). POSSIBLE RETURNS ARE       
C                      IER = 0, SUCCESSFUL COMPLETION     
C                          = 703, NZ TOO SMALL - NO MORE  
C                                 ELEMENTS CAN BE ADDED   
C       
C ... SPECIFICTIONS FOR ARGUMENTS     
C       
      INTEGER NZ,IA(1),JA(1),IWORK(NZ),N,LEVELL,NOUTT,IERR
      DOUBLE PRECISION A(NZ)
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IER,J,LEVEL,NOUT,NADD,NADDP1,NOW,NP1,NTO,NTN
C       
C ... INITIALIZE LOCAL VARIABLES AND MAKE ERROR CHECK     
C       
      NOW = IA(N+1)-1       
      NADD = NZ-NOW 
      IER = 0     
      LEVEL = LEVELL
      NOUT = NOUTT
      IF (NADD.LE.0) IER = 703
      IF (IER.EQ.0) GO TO 20
      IF (LEVEL.GE.0) WRITE (NOUT,10) IER       
   10 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE SBAGN   '/' ','    IER = ',I10/' ', 
     *   '    NZ TOO SMALL - NO ROOM FOR NEW ENTRY')      
      GO TO 90    
C       
C ... SHIFT ELEMENTS OF A AND JA DOWN AND ADD ZERO FILL   
C       
   20 NTO = NOW   
      NTN = NZ    
      DO 30 I = 1,NOW       
         JA(NTN) = JA(NTO)  
         A(NTN) = A(NTO)    
         NTO = NTO-1
         NTN = NTN-1
   30 CONTINUE    
      DO 40 I = 1,NADD      
         JA(I) = 0
         A(I) = 0.D0
   40 CONTINUE    
C       
C ... UPDATE IA TO REFLECT DOWNWARD SHIFT IN A AND JA     
C       
      NP1 = N+1   
      DO 50 I = 1,NP1       
         IA(I) = IA(I)+NADD 
   50 CONTINUE    
C       
C ... CREATE LINKED LIST    
C       
      NADDP1 = NADD+1       
      DO 60 I = NADDP1,NZ   
         IWORK(I) = I+1     
   60 CONTINUE    
      DO 70 I = 1,NADD      
         IWORK(I) = 0       
   70 CONTINUE    
      DO 80 I = 1,N 
         J = IA(I+1)-1      
         IWORK(J) = -I      
   80 CONTINUE    
C       
C ... INDICATE IN LAST POSITION OF IA HOW MANY SPACES     
C     ARE LEFT IN A AND JA FOR ADDITION OF ELEMENTS       
C       
      IA(N+1) = NADD
      RETURN      
C       
C ... ERROR RETURN
C       
   90 IERR = IER  
      RETURN      
      END 
      SUBROUTINE SBELM (NN,IA,JA,A,RHS,IW,RW,TOL,ISYM,LEVEL,NOUT,IER) 
C       
C ... SBELM IS DESIGNED TO REMOVE ROWS AND COLUMNS OF THE MATRIX    
C ... WHERE DABS(A(I,J))/A(I,I) .LE. TOL FOR J = 1 TO N AND A(I,I)  
C ... .GT. 0. THIS IS TO TAKE CARE OF MATRICES ARISING    
C ... FROM FINITE ELEMENT DISCRETIZATIONS OF PDE^S WITH DIRICHLET   
C ... BOUNDARY CONDITIONS.  ANY SUCH ROWS AND CORRESPONDING COLUMNS 
C ... ARE THEN SET TO THE IDENTITY AFTER CORRECTING RHS.  
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C          IW,RW  WORK ARRAYS OF LENGTH N       
C          TOL    TOLERANCE FACTOR    
C          ISYM   FLAG FOR TYPE OF STORAGE FOR SYSTEM     
C                 (0: SYMMETRIC, 1:NONSYMMETRIC)
C          LEVEL  PRINTING SWITCH FOR ERROR CONDITION     
C          NOUT OUTPUT TAPE NUMBER    
C          IER    ERROR FLAG: NONZERO VALUE ON RETURN MEANS 
C                    101 : DIAGONAL ENTRY NOT POSITIVE    
C                    102 : THERE IS NO DIAGONAL ENTRY IN ROW
C       
C********************************************************************** 
C       
C     UPDATE.  SBELM HAS BEEN REWRITTEN TO SPEED UP THE LOCATION OF 
C              OF ROWS WHICH ARE TO BE ELIMINATED.  THIS IS DONE BY 
C              FIRST STORING THE LARGEST ELEMENT OF EACH ROW IN     
C              THE ARRAY RW.  THE DIAGONAL ENTRY IS THEN COMPARED   
C              WITH THE CORRESPONDING ELEMENT IN RW.  IF IT IS      
C              DECIDED TO ELIMINATE THE ROW THEN IT IS MARKED FOR   
C              ELIMINATION. 
C       
C              WHEN A ROW IS TO BE ELIMINATED ITS DIAGONAL ENTRY    
C              IS STORED IN  RW  AND  IW IS MARKED BY A NONZERO     
C              (WHICH IS THIS ROW NUMBER)       
C       
C              ROWS WHICH HAVE ONLY DIAGONAL ENTRIES ARE NOT
C              ALTERED.     
C       
C*********************************************************************
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER NN,IA(1),JA(1),IW(NN),ISYM,LEVEL,NOUT,IER   
      DOUBLE PRECISION A(1),RHS(NN),RW(NN),TOL  
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,ICNT,IEND,JJ,JJDI,KK,N       
      DOUBLE PRECISION DI   
C       
      N = NN      
C       
C        IF (N .GE. 1) GO TO 10       
C           IER = 100       
C           RETURN
C 10     CONTINUE 
C       
C ... STORE THE LARGEST (DABSOLUTE VALUE) OFF DIAGONAL ENTRY FOR    
C ... ROW II IN RW(II).     
C       
      IER = 0     
      ICNT = 0    
      DO 10 II = 1,N
         RW(II) = 0.0D0     
         IW(II) = 0 
   10 CONTINUE    
      DO 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 140  
         DO 20 JJ = IBGN,IEND 
            KK = JA(JJ)     
            IF (KK.EQ.II) GO TO 20    
            RW(II) = DMAX1(RW(II),DABS(A(JJ)))  
            IF (ISYM.NE.0) GO TO 20   
            RW(KK) = DMAX1(RW(KK),DABS(A(JJ)))  
   20 CONTINUE    
C       
C ... FOR II = 1 TO N FIND THE DIAGONAL ENTRY IN ROW II   
C       
      DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         DO 40 JJ = IBGN,IEND 
            IF (JA(JJ).NE.II) GO TO 40
            DI = A(JJ)      
            JJDI = JJ       
            IF (DI.GT.0.D0) GO TO 50  
            IER = 101       
            IF (LEVEL.GE.0) WRITE (NOUT,30) II,DI 
   30       FORMAT ('0','*** F A T A L     E R R O R ************'/'0', 
     *         '    IN ITPACK ROUTINE SBELM   '/' ',      
     *         '    DIAGONAL ELEMENT',I10,' NOT POSITIVE  '/' ',    
     *         '    CURRENT VALUE = ',D15.8)    
            RETURN
   40    CONTINUE 
         GO TO 140
   50    CONTINUE 
C       
C ... CHECK THE SIZE OF THE LARGEST OFF DIAGONAL ELEMENT  
C ... ( STORED IN RW(II) ) AGAINST THE DIAGONAL ELEMENT DII.
C       
         IF (RW(II).NE.0.0D0) GO TO 60
         IF (1.0D0/DI.LE.TOL) GO TO 70
         GO TO 80 
   60    IF (RW(II)/DI.GT.TOL) GO TO 80 
C       
C ... THE OFF DIAGONAL ELEMENTS ARE SMALL COMPARED TO THE DIAGONAL  
C ... THEREFORE MARK IT FOR ELIMINATION AND PERFORM INITIAL 
C ... PROCESSING  
C       
   70    ICNT = ICNT+1      
         IW(II) = II
         RW(II) = DI
         A(JJDI) = 1.0D0    
         RHS(II) = RHS(II)/DI 
C       
   80 CONTINUE    
C       
C ... ELIMINATE THE ROWS AND COLUMNS INDICATED BY THE NONZERO       
C ... ENTRIES IN IW.  THERE ARE ICNT OF THEM    
C       
      IF (ICNT.EQ.0) GO TO 130
C       
C ... THE ELIMINATION IS AS FOLLOWS:  
C       
C     FOR II = 1 TO N DO    
C        IF ( IW(II) .NE. 0 ) THEN    
C           SET DIAGONAL VALUE TO 1.0  ( ALREADY DONE )   
C           SET RHS(II) = RHS(II) / RW(II)   ( ALREADY DONE )       
C           FIND NONZERO OFFDIAGONAL ENTRIES  KK
C           IF ( IW(KK) .EQ. 0 ) FIX UP RHS(KK)  WHEN USING SYMMETRIC ST
C           SET A(II,KK) = 0.0
C        ELSE ( I.E.  IW(II) .EQ. 0  )
C           FIND NONZERO OFFDIAGONAL ENTRIES   KK 
C           IF ( IW(KK) .NE. 0 ) FIX UP RHS(II) 
C                                AND SET A(II,KK) = 0.0   
C        END IF   
C     END DO      
C       
      DO 120 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IW(II).EQ.0) GO TO 100   
C       
C ... THE II-TH ROW IS TO BE ELIMINATED 
C       
         DO 90 JJ = IBGN,IEND 
            KK = JA(JJ)     
            IF (KK.EQ.II) GO TO 90    
            IF ((IW(KK).EQ.0).AND.(ISYM.EQ.0)) RHS(KK) = RHS(KK)-A(JJ)* 
     *         RHS(II)      
            A(JJ) = 0.0D0   
   90    CONTINUE 
         GO TO 120
C       
C ... THE II-TH ROW IS KEPT.  CHECK THE OFF-DIAGONAL ENTRIES
C       
  100    DO 110 JJ = IBGN,IEND
            KK = JA(JJ)     
            IF (KK.EQ.II.OR.IW(KK).EQ.0) GO TO 110
            RHS(II) = RHS(II)-A(JJ)*RHS(KK)     
            A(JJ) = 0.0D0   
  110    CONTINUE 
C       
  120 CONTINUE    
C       
  130 RETURN      
C       
C ... ERROR TRAPS -- NO DIAGONAL ENTRY IN ROW II (ROW MAY BE EMPTY).
C       
  140 CONTINUE    
      IER = 102   
      IF (LEVEL.GE.0) WRITE (NOUT,150) II       
  150 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE SBELM   '/' ',  
     *   '    NO DIAGONAL ENTRY IN ROW  ',I10)  
C       
      RETURN      
      END 
      SUBROUTINE SBEND (N,NZ,IA,JA,A,IWORK)     
C       
C***********************************************************************
C       
C     SBEND IS THE THIRD OF A SUITE OF SUBROUTINES TO AID THE       
C     USER TO CONSTRUCT THE  IA, JA, A DATA STRUCTURE USED IN       
C     ITPACK.     
C       
C     SBEND RESTRUCTURES THE LINKED LIST DATA STRUCTURE BUILT BY    
C     SBINI AND SBSIJ INTO THE FINAL DATA STRUCTURE REQUIRE BY      
C     ITPACK.  THE RESTRUCTURING CAN TAKE PLACE IN THE MINIMUM      
C     AMOUNT OF MEMORY REQUIRED TO HOLD THE NONZERO STRUCTURE OF    
C     THE SPARSE MATRIX BUT WILL RUN QUICKER IF MORE STORAGE
C     IS ALLOWED. 
C       
C     SBEND IS BASED ON SUBROUTINE BUILD OF THE SPARSE MATRIX       
C     PACKAGE SPARSPAK DEVELOPED BY ALAN GEORGE AND JOSEPH LUI      
C     OF THE UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO.   
C       
C ... PARAMETERS  
C       
C ...... INPUT    
C       
C     N       THE ORDER OF THE LINEAR SYSTEM    
C       
C     NZ      THE LENGTH OF THE ARRAYS JA, IWORK, AND A.  
C       
C ...... INPUT/OUTPUT       
C       
C     IA      INTEGER ARRAY OF LENGTH N+1.  THE FIRST N ENTRIES     
C             POINT TO THE BEGINNING OF THE LINKED LIST FOR EACH    
C             ROW.  IA(N+1)-1 IS THE TOP OF THE LINKED LISTS
C             CONTAINED IN JA, IWORK, AND A.  ON OUTPUT IA WILL     
C             POINT TO THE FIRST ENTRY OF EACH ROW IN THE FINAL     
C             DATA STRUCTURE. 
C       
C     JA      INTEGER ARRAY OF LENGTH NZ.  ON INPUT JA STORES THE   
C             COLUMN NUMBERS OF THE NONZERO ENTRIES AS INDICATED    
C             BY THE LINKED LISTS.  ON OUTPUT JA STORES THE 
C             COLUMN NUMBERS IN ROW ORDERED FORM. 
C       
C     A       D.P. ARRAY OF LENGTH NZ.  ON INPUT A STORES THE       
C             VALUE OF THE NOZERO ENTRIES AS INDICATED BY THE       
C             LINKED LISTS.  ON OUTPUT A STORES THE VALUES IN       
C             ROW ORDERED FORM.       
C       
C     IWORK    INTEGER ARRAY OF LENGTH NZ.  ON INPUT IWORK STORES THE 
C             THE LINKS OF THE LINKED LISTS.  ON OUTPUT IT IS       
C             DESTROYED.    
C       
C***********************************************************************
C       
      INTEGER N,NZ,IA(1),JA(NZ),IWORK(NZ)       
      DOUBLE PRECISION A(NZ)
C       
      INTEGER MAXTOP,NEXT,TOP,IDEG,NULINK,JAJ,HLINK,OHLINK,L,I,LINK,
     *   MHLINK   
      DOUBLE PRECISION VAL  
C       
C***********************************************************************
C       
C ... INITIALIZATION
C       
C ...... THE VARIABLES NEXT AND TOP RESPECTIVELY POINT TO THE       
C        NEXT AVAILABLE ENTRY FOR THE FINAL DATA STRUCTURE AND      
C        THE TOP OF THE REMAINDER OF THE LINKED LISTS.    
C       
      NEXT = 1    
      TOP = IA(N+1)+1       
      MAXTOP = NZ-IA(N+1)+1 
C       
C***********************************************************************
C       
C ... CONVERT EACH ROW INTO FINAL FORM
C       
      DO 90 I = 1,N 
         IDEG = 0 
         NULINK = IA(I)     
C       
C ... LOOP OVER EACH NODE IN THE LINKED LIST OF ROW I     
C       
   10    LINK = NULINK      
         IF (LINK.LE.0) GO TO 80      
         NULINK = IWORK(LINK) 
         JAJ = JA(LINK)     
         VAL = A(LINK)      
C       
C ... CHECK TO SEE IF A COLLISION BETWEEN THE LINKED LISTS
C     AND THE FINAL FORM HAS OCCURRED.
C       
         IF (NEXT.GE.TOP.AND.LINK.NE.TOP) GO TO 20
C       
C ... COLLISION HAS NOT OCCURRED.  FREE THE SPACE FOR THE TRIPLE    
C     (JA(LINK), A(LINK), IWORK(LINK))
C       
         JA(LINK) = 0       
         A(LINK) = 0.0D0    
         IWORK(LINK) = 0    
C       
C ... SPECIAL CASE TO MOVE  TOP  DOWN IF LINK .EQ. TOP    
C       
         IF (LINK.EQ.TOP) GO TO 60    
         GO TO 70 
C       
C***********************************************************************
C       
C ... COLLISION HAS OCCURRED.  CLEAR OFF SOME SPACE FOR THE CURRENT 
C     ENTRY BY MOVING THE TRIPLE ( JA(TOP),A(TOP),IWORK(TOP) )      
C     DOWNWARDS TO THE FREED TRIPLE ( JA(LINK),A(LINK),IWORK(LINK) ). 
C     THEN ADJUST THE LINK FIELDS.    
C       
C ...... PATCH UP THE LINKED LIST FOR THE CURRENT ROW I.  THEN      
C        TRAVERSE THE LINKED LIST CONTAINING TOP UNTIL THE POINTER  
C        POINTER BACK TO IA IS FOUND. 
C       
   20    IA(I) = LINK       
         HLINK = TOP
C       
   30    HLINK = IWORK(HLINK) 
         IF (HLINK.GT.0) GO TO 30     
C       
C ...... NOW FOLLOW THE LINKED LIST BACK TO TOP KEEPING TRACK       
C        OF THE OLD LINK.   
C       
C ......... SPECIAL CASE IF IA(-HLINK) = TOP    
C       
         MHLINK = -HLINK    
         IF (IA(MHLINK).NE.TOP) GO TO 40
C       
         IWORK(LINK) = IWORK(TOP)     
         JA(LINK) = JA(TOP) 
         A(LINK) = A(TOP)   
         IA(MHLINK) = LINK  
         IF (NULINK.EQ.TOP) NULINK = LINK       
         GO TO 60 
C       
C ......... USUAL CASE.     
C       
   40    HLINK = IA(MHLINK) 
   50    OHLINK = HLINK     
         HLINK = IWORK(OHLINK)
         IF (HLINK.NE.TOP) GO TO 50   
C       
         IWORK(LINK) = IWORK(TOP)     
         JA(LINK) = JA(TOP) 
         A(LINK) = A(TOP)   
         IF (OHLINK.NE.LINK) IWORK(OHLINK) = LINK 
         IF (NULINK.EQ.TOP) NULINK = LINK       
C       
C ... COLLAPSE TOP OF LINK LIST BY AS MUCH AS POSSIBLE    
C       
   60    TOP = TOP+1
         IF (TOP.GE.MAXTOP) GO TO 70  
         IF (IWORK(TOP).NE.0) GO TO 70
         GO TO 60 
C       
C***********************************************************************
C       
C ... PUT THE CURRENT TRIPLE INTO THE FINAL DATA STRUCTURE
C       
   70    JA(NEXT) = JAJ     
         A(NEXT) = VAL      
         NEXT = NEXT+1      
         IDEG = IDEG+1      
         GO TO 10 
C       
C ... FINAL STRUCTURE FOR ROW I IS COMPLETE.  LINKED LIST IS
C     DESTROYED AND WILL BE RECAPTURED AS NECESSARY BY THE
C     LOOP ON LABEL 60      
C       
   80    IA(I) = IDEG       
C       
   90 CONTINUE    
C       
C***********************************************************************
C       
C ... FINALIZE THE DATA STRUCTURE BY BUILDING THE FINAL VERSION OF  
C     IA. 
C       
      L = IA(1)+1 
      IA(1) = 1   
      DO 100 I = 1,N
         IDEG = IA(I+1)     
         IA(I+1) = L
         L = L+IDEG 
  100 CONTINUE    
C       
C ... FINAL IA, JA, A DATA STRUCTURE BUILT.     
C       
      RETURN      
      END 
      SUBROUTINE SBINI (N,NZ,IA,JA,A,IWORK)     
C       
C***********************************************************************
C       
C     SBINI IS THE FIRST OF A SUITE OF THREE SUBROUTINES TO AID     
C     THE USER TO CONSTRUCT THE IA, JA, A DATA STRUCTURE USED       
C     IN ITPACK.  
C       
C     SBINI INITIALIZES THE ARRAYS IA, JA, IWORK, AND A.  THE OTHER 
C     SUBROUTINES IN THE SUITE ARE SBSIJ ( WHICH BUILDS A LINKED    
C     LIST REPRESENTATION OF THE MATRIX STRUCTURE ) AND SBEND ( WHICH 
C     RESTRUCTURE THE LINKED LIST FORM INTO THE FINAL FORM ).       
C       
C ... PARAMETERS  
C       
C ...... INPUT    
C       
C     N          THE ORDER OF THE LINEAR SYSTEM 
C       
C     NZ         THE MAXIMUM NUMBER OF NONZEROES ALLOWED IN THE     
C                LINEAR SYSTEM.       
C       
C ...... OUTPUT   
C       
C     IA         INTEGER ARRAY OF LENGTH N+1.  SBINI SETS THIS ARRAY
C                TO -I FOR I = 1 THRU N.  IA(N+1) IS SET TO NZ.     
C       
C     JA         INTEGER ARRAY OF LENGTH NZ.  INITIALIZED TO ZERO HERE. 
C       
C     A          D.P. ARRAY OF LENGTH NZ.  INITIALIZED TO ZERO HERE.
C       
C     IWORK       INTEGER ARRAY OF LENGTH NZ.  INITIALIZED TO ZERO HERE.
C       
C***********************************************************************
C       
      INTEGER N,NZ,IA(1),JA(NZ),IWORK(NZ),I     
      DOUBLE PRECISION A(NZ)
C       
C***********************************************************************
C       
      DO 10 I = 1,N 
         IA(I) = -I 
   10 CONTINUE    
      IA(N+1) = NZ
C       
      CALL IVFILL (NZ,JA,0) 
      CALL IVFILL (NZ,IWORK,0)
      CALL VFILL (NZ,A,0.D0)
C       
      RETURN      
      END 
      SUBROUTINE SBSIJ (N,NZ,IA,JA,A,IWORK,II,JJ,VALL,MODE,LEVELL,NOUTT,
     *   IERR)    
C       
C***********************************************************************
C       
C     SBSIJ IS THE SECOND OF A SUITE OF THREE SUBROUTINES TO AID IN 
C     THE CONSTRUCTION OF THE IA, JA, A DATA STRUCTURE USED IN      
C     ITPACK.     
C       
C     SBSIJ TAKES THE INDIVIDUAL ENTRIES OF THE SPARSE MATRIX AS    
C     GIVEN TO IT AT EACH CALL VIA  (I,J,VAL) AND INSERTS IT INTO   
C     A LINKED LIST REPRESENTATION OF THE SPARSE MATRIX.  
C       
C     EACH ROW OF THE SPARSE MATRIX IS ASSOCIATED WITH A CIRCULAR   
C     LINKED LIST BEGINNING AT IA(I).  THE LAST ENTERED ELEMENT IN  
C     EACH LIST POINTS BACK TO IA(I) WITH THE VALUE -I.  THE LINKS  
C     ARE STORED IN THE ARRAY IWORK, WHILE JA AND A STORE THE COLUMN
C     NUMBER AND VALUE IN PARALLEL TO IWORK.  THE LINKED LISTED ARE 
C     STORED BEGINNING AT ENTRY NZ AND WORKING BACKWARDS TOWARDS 1. 
C       
C ... PARAMETERS  
C       
C ...... INPUT    
C       
C     N       THE ORDER OF THE LINEAR SYSTEM    
C       
C     NZ      THE LENGTH OF THE ARRAYS  JA, A, AND IWORK  
C       
C     I, J    THE ROW AND COLUMN NUMBERS OF THE ENTRY OF THE SPARSE 
C             LINEAR SYSTEM TO BE ENTERED IN THE DATA STRUCTURE(=II,JJ) 
C       
C     VAL     THE NONZERO VALUE ASSOCIATED WITH (I,J)  (= VALL)     
C       
C     MODE    IF THE (I,J) ENTRY HAS ALREADY BEEN SET, MODE SPECIFIES 
C             THE WAY IN WHICH THE ENTRY IS TO BE TREATED.
C             IF   MODE .LT. 0  LET THE VALUE REMAIN AS IS
C                       .EQ. 0  RESET IT TO THE NEW VALUE 
C                       .GT. 0  ADD THE NEW VALUE TO THE OLD VALUE  
C       
C     NOUT  OUTPUT FILE NUMBER (= NOUTT)
C       
C     LEVEL   OUTPUT FILE SWITCH (= LEVELL)     
C ... INPUT/OUTPUT
C       
C     IA      INTEGER ARRAY OF LENGTH N+1.  THE FIRST N ENTRIES     
C             POINT TO THE BEGINNING OF THE LINKED LIST FOR EACH    
C             ROW.  IA(N+1) POINTS TO THE NEXT ENTRY AVAILABLE FOR  
C             STORING THE CURRENT ENTRY INTO THE LINKED LIST.       
C       
C     JA      INTEGER ARRAY OF LENGTH NZ.  JA STORES THE COLUMN     
C             NUMBERS OF THE NONZERO ENTRIES.   
C       
C     A       D.P. ARRAY OF LENGTH NZ.  A STORES THE VALUE OF THE   
C             NONZERO ENTRIES.
C       
C     IWORK   INTEGER ARRAY OF LENGTH NZ. IWORK STORES THE LINKS.   
C       
C     IER     ERROR FLAG.(= IERR)  POSSIBLE RETURNS ARE   
C             IER =    0   SUCCESSFUL COMPLETION
C                 =  700   ENTRY WAS ALREADY SET,  VALUE HANDLED    
C                          AS SPECIFIED BY MODE.
C                 =  701   IMPROPER VALUE OF EITHER I OR J INDEX    
C                 =  702   NO ROOM REMAINING, NZ TOO SMALL. 
C       
C***********************************************************************
C       
      INTEGER N,NZ,IA(1),JA(NZ),IWORK(NZ),II,JJ,MODE,LEVELL,NOUTT,IERR
      DOUBLE PRECISION A(NZ),VALL     
C       
      INTEGER LINK,NEXT,NPL1,I,J,LEVEL,NOUT,IER 
      DOUBLE PRECISION VAL,TEMP       
C       
C***********************************************************************
C       
C ... CHECK THE VALIDITY OF THE (I,J) ENTRY     
C       
      I = II      
      J = JJ      
      VAL = VALL  
      LEVEL = LEVELL
      NOUT = NOUTT
      IER = 0     
      IF (I.LE.0.OR.I.GT.N) IER = 701 
      IF (J.LE.0.OR.J.GT.N) IER = 701 
      IF (IER.EQ.0) GO TO 20
      IF (LEVEL.GE.0) WRITE (NOUT,10) IER,I,J   
   10 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    ( ',I10,' , ',I10,' )'/' ',       
     *   '    IMPROPER VALUE FOR I OR J ')      
      GO TO 130   
C       
C ... TRAVERSE THE LINK LIST POINTED TO BY IA(I) UNTIL EITHER       
C ... THE J ENTRY OR THE END OF THE LIST HAS BEEN FOUND.  
C       
   20 NPL1 = N+1  
      LINK = IA(I)
C       
C ...... SPECIAL CASE FOR THE FIRST ENTRY IN THE ROW      
C       
      IF (LINK.GT.0) GO TO 30 
      NEXT = IA(NPL1)       
      IF (NEXT.LT.1) GO TO 110
C       
      IA(I) = NEXT
      JA(NEXT) = J
      A(NEXT) = VAL 
      IWORK(NEXT) = -I      
      IA(NPL1) = NEXT-1     
      GO TO 130   
C       
C ... FOLLOW THE LINK LIST UNTIL J OR THE END OF THE LIST IS FOUND  
C       
   30 IF (JA(LINK).EQ.J) GO TO 40     
      IF (IWORK(LINK).LE.0) GO TO 100 
      LINK = IWORK(LINK)    
      GO TO 30    
C       
C:      
C ... ENTRY (I,J) ALREADY HAS BEEN SET.  RESET VALUE DEPENDING ON MODE
C       
   40 IER = 700   
      IF (MODE.GE.0) GO TO 60 
      IF (LEVEL.GE.1) WRITE (NOUT,50) IER,I,J,A(LINK)     
   50 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    ( ',I10,' , ',I10,' )'/' ',       
     *   '    ENTRY ALREADY SET AND IS LEFT AS ',D15.8)   
      GO TO 130   
   60 IF (MODE.GE.1) GO TO 80 
      IF (LEVEL.GE.1) WRITE (NOUT,70) IER,I,J,A(LINK),VAL 
   70 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    ( ',I10,' , ',I10,' )'/' ',       
     *   '    ENTRY ALREADY SET - CURRENT VALUE OF',D15.8/' ',      
     *   '                                RESET TO',D15.8)
      A(LINK) = VAL 
      GO TO 130   
   80 TEMP = A(LINK)+VAL    
      IF (LEVEL.GE.1) WRITE (NOUT,90) IER,I,J,A(LINK),TEMP
   90 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    ( ',I10,' , ',I10,' )'/' ',       
     *   '    ENTRY ALREADY SET - CURRENT VALUE OF',D15.8/' ',      
     *   '                                RESET TO',D15.8)
      A(LINK) = TEMP
      GO TO 130   
C       
C ... ENTRY (I,J) HAS NOT BEEN SET.  ENTER IT INTO THE LINKED LIST  
C       
  100 NEXT = IA(NPL1)       
      IF (NEXT.LT.1) GO TO 110
C       
      IWORK(LINK) = NEXT    
      JA(NEXT) = J
      A(NEXT) = VAL 
      IWORK(NEXT) = -I      
      IA(NPL1) = NEXT-1     
      GO TO 130   
C       
C***********************************************************************
C       
C ... ERROR TRAP FOR NO ROOM REMAINING
C       
  110 IER = 702   
      IF (LEVEL.GE.0) WRITE (NOUT,120) IER      
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   
     *   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', 
     *   '    NZ TOO SMALL - NO ROOM FOR NEW ENTRY')      
C       
  130 CONTINUE    
      IERR = IER  
      RETURN      
      END 
      SUBROUTINE SCAL (NN,IA,JA,A,RHS,U,D,LEVEL,NOUT,IER) 
C       
C ... ORIGINAL MATRIX IS SCALED TO A UNIT DIAGONAL MATRIX.  RHS     
C ... AND U ARE SCALED ACCORDINGLY.  THE MATRIX IS THEN SPLIT AND   
C ... IA, JA, AND A RESHUFFLED.       
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX (= NN)    
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C          U      LATEST ESTIMATE OF SOLUTION   
C          D      OUTPUT VECTOR CONTAINING THE SQUARE ROOTS 
C                    OF THE DIAGONAL ENTRIES    
C          LEVEL  PRINTING SWITCH FOR ERROR CONDITION     
C          NOUT OUTPUT TAPE NUMBER    
C          IER    ERROR FLAG: ON RETURN NONZERO VALUES MEAN 
C                    401 : THE ITH DIAGONAL ELEMENT IS .LE. 0.      
C                    402 : NO DIAGONAL ELEMENT IN ROW I   
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),NN,LEVEL,NOUT,IER     
      DOUBLE PRECISION A(1),RHS(NN),U(NN),D(NN) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,IBGN,IEND,II,IM1,J,JADD,JAJJ,JJ,JJPI,N,NP1
      DOUBLE PRECISION DI   
C       
C ... EXTRACT SQUARE ROOT OF THE DIAGONAL OUT OF A AND SCALE U AND RHS
C       
      N = NN      
      IER = 0     
      DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 50   
         DO 40 JJ = IBGN,IEND 
            IF (JA(JJ).NE.II) GO TO 40
            DI = A(JJ)      
            IF (DI.GT.0.D0) GO TO 70  
            IF (DI.EQ.0.D0) GO TO 20  
            IER = 401       
            IF (LEVEL.GE.0) WRITE (NOUT,10) II  
   10       FORMAT ('0','*** F A T A L     E R R O R ************'/'0', 
     *         '    IN ITPACK ROUTINE SCAL    '/' ',      
     *         '    DIAGONAL ENTRY IN ROW ',I10,' NEGATIVE')
            RETURN
   20       IER = 401       
            IF (LEVEL.GE.0) WRITE (NOUT,30)     
   30       FORMAT ('0','*** F A T A L     E R R O R ************'/'0', 
     *         '    IN ITPACK ROUTINE SCAL    '/' ',      
     *         '    DIAGONAL ENTRY IN ROW ',I10,' IS ZERO') 
            RETURN
   40    CONTINUE 
   50    IER = 402
         IF (LEVEL.GE.0) WRITE (NOUT,60) II     
   60    FORMAT ('0','*** F A T A L     E R R O R ************'/'0',
     *      '    IN ITPACK ROUTINE SCAL    '/' ', 
     *      '    NO DIAGONAL ENTRY IN ROW',I10) 
         RETURN   
C       
   70    CONTINUE 
         DI = DSQRT(DABS(DI)) 
         RHS(II) = RHS(II)/DI 
         U(II) = U(II)*DI   
         D(II) = DI 
   80 CONTINUE    
C       
C ... SHIFT MATRIX TO ELIMINATE DIAGONAL ENTRIES
C       
      IF (N.EQ.1) GO TO 110 
      NP1 = N+1   
      DO 100 I = 1,N
         IM1 = I-1
         II = NP1-I 
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         JADD = IBGN+IEND   
         DO 90 J = IBGN,IEND
            JJ = JADD-J     
            JJPI = JJ+IM1   
            IF (JA(JJ).EQ.II) IM1 = I 
            A(JJPI) = A(JJ) 
            JA(JJPI) = JA(JJ) 
   90    CONTINUE 
         IA(II+1) = IA(II+1)+I-1      
  100 CONTINUE    
  110 IA(1) = IA(1)+N       
C       
C ... SCALE SHIFTED MATRIX AND STORE D ARRAY IN FIRST N ENTRIES OF A
C       
      DO 140 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         DI = D(II) 
         IF (IBGN.GT.IEND) GO TO 130  
         DO 120 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            A(JJ) = A(JJ)/(DI*D(JAJJ))
  120    CONTINUE 
  130    CONTINUE 
         A(II) = DI 
  140 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE SUM3 (N,C1,X1,C2,X2,C3,X3)     
C       
C ... COMPUTES X3 = C1*X1 + C2*X2 + C3*X3       
C       
C ... PARAMETER LIST:       
C       
C          N        INTEGER LENGTH OF VECTORS X1, X2, X3  
C          C1,C2,C3 D.P. CONSTANTS    
C          X1,X2,X3 D.P. VECTORS SUCH THAT      
C                   X3(I) = C1*X1(I) + C2*X2(I) + C3*X3(I)
C                   X3(I) = C1*X1(I) + C2*X2(I)  IF C3 = 0. 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I   
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
      DOUBLE PRECISION X1(N),X2(N),X3(N),C1,C2,C3 
C       
      IF (N.LE.0) RETURN    
      IF (DABS(C3).EQ.0.D0) GO TO 20  
C       
      DO 10 I = 1,N 
         X3(I) = C1*X1(I)+C2*X2(I)+C3*X3(I)     
   10 CONTINUE    
      RETURN      
C       
C ... COMPUTE X3 = C1*X1 + C2*X2      
C       
   20 DO 30 I = 1,N 
         X3(I) = C1*X1(I)+C2*X2(I)    
   30 CONTINUE    
C       
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION TAU (II)
C       
C ... THIS SUBROUTINE SETS TAU(II) FOR THE SOR METHOD.    
C       
C ... PARAMETER LIST:       
C       
C          II     NUMBER OF TIMES PARAMETERS HAVE BEEN CHANGED      
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER II  
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      DOUBLE PRECISION T(8) 
C       
      DATA T(1),T(2),T(3),T(4),T(5),T(6),T(7),T(8) / 1.5D0,1.8D0,1.85D0,
     *   1.9D0,1.94D0,1.96D0,1.975D0,1.985D0 /  
C       
      TAU = 1.992D0 
      IF (II.LE.8) TAU = T(II)
C       
      RETURN      
      END 
      FUNCTION TIMER (TIMDMY) 
C
C ... TIMER IS A ROUTINE TO RETURN THE EXECUTION TIME IN
C ... SECONDS.
C
C ... PARAMETERS -- 
C
C          TIMDMY   DUMMY ARGUMENT
C
C
C *********************************************
C **                                         **
C **   THIS ROUTINE IS NOT PORTABLE.         **
C **                                         **
C *********************************************
C
      REAL TIMDMY
C
C ... CRAY Y-MP.
C
C     TIMER = SECOND ()
C
C ... UNIX ETIME FACILITY.
C
C      EXTERNAL ETIME
      DIMENSION TARRAY(2)
C      REAL ETIME, TIMER
      REAL TIMER
C      TOTAL = ETIME (TARRAY)
      TOTAL = 0.0
      TIMER = TOTAL
C
C ... IBM RISC SYSTEM/6000.
C
C     TIMER = FLOAT(MCLOCK())/100.0
C
      RETURN
      END 
      LOGICAL FUNCTION TSTCHG (IBMTH) 
C       
C     THIS FUNCTION PERFORMS A TEST TO DETERMINE IF PARAMETERS      
C     SHOULD BE CHANGED FOR SEMI-ITERATION ACCELERATED METHODS.     
C       
C ... PARAMETER LIST:       
C       
C          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY SI 
C                      IBMTH = 1,   JACOBI      
C                            = 2,   REDUCED SYSTEM
C                            = 3,   SSOR
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IBMTH 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IP  
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
C       
      IP = IN-IS  
      IF (IBMTH.EQ.2) IP = 2*IP       
C       
      IF (IN.EQ.0) GO TO 10 
      IF (IP.LT.3) GO TO 20 
C       
      QA = DSQRT(DABS(DELNNM/DELSNM)) 
      QT = 2.D0*DSQRT(DABS(RRR**IP))/(1.D0+RRR**IP)       
      IF ((QA.GE.1.D0).OR.(QA.LT.QT**FF)) GO TO 20
C       
C ... TEST PASSES -- CHANGE PARAMETERS
C       
   10 TSTCHG = .TRUE.       
      RETURN      
C       
C ... TEST FAILS -- DO NOT CHANGE PARAMETERS    
C       
   20 TSTCHG = .FALSE.      
      RETURN      
C       
      END 
      SUBROUTINE UNSCAL (N,IA,JA,A,RHS,U,D)     
C       
C ... THIS SUBROUTINE REVERSES THE PROCESS OF SCAL.       
C       
C ... PARAMETER LIST:       
C       
C          N      DIMENSION OF MATRIX 
C          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
C          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
C          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
C          U      LATEST ESTIMATE OF SOLUTION   
C          D      VECTOR CONTAINING THE SQUARE ROOTS      
C                    OF THE DIAGONAL ENTRIES    
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER IA(1),JA(1),N 
      DOUBLE PRECISION A(1),RHS(N),U(N),D(N)    
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER IBGN,IEND,II,INEW,IS,JAJJ,JJ,JJPI 
      DOUBLE PRECISION DI   
C       
C ... EXTRACT DIAGONAL FROM SCALED A AND UNSCALE U AND RHS
C       
      DO 10 II = 1,N
         DI = A(II) 
         U(II) = U(II)/DI   
         RHS(II) = RHS(II)*DI 
         D(II) = DI 
   10 CONTINUE    
C       
C ... UNSCALE A   
C       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 30   
         DI = D(II) 
         DO 20 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            A(JJ) = A(JJ)*DI*D(JAJJ)  
   20    CONTINUE 
   30 CONTINUE    
C       
C ... INSERT DIAGONAL BACK INTO A     
C       
      DO 60 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IS = N-II
         INEW = IBGN-IS-1   
         A(INEW) = D(II)**2 
         JA(INEW) = II      
         IF (IS.EQ.0.OR.IBGN.GT.IEND) GO TO 50  
         DO 40 JJ = IBGN,IEND 
            JJPI = JJ-IS    
            A(JJPI) = A(JJ) 
            JA(JJPI) = JA(JJ) 
   40    CONTINUE 
   50    CONTINUE 
         IA(II) = INEW      
   60 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE VEVMW (N,V,W)
C       
C ... VEVMW COMPUTES V = V - W
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTORS V AND W       
C          V      D.P. VECTOR 
C          W      D.P. VECTOR SUCH THAT   V(I) = V(I) - W(I)
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N   
      DOUBLE PRECISION V(N),W(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,M,MP1       
C       
      IF (N.LE.0) RETURN    
      M = MOD(N,4)
C       
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         V(I) = V(I)-W(I)   
   10 CONTINUE    
      IF (N.LT.4) RETURN    
C       
   20 MP1 = M+1   
      DO 30 I = MP1,N,4     
         V(I) = V(I)-W(I)   
         V(I+1) = V(I+1)-W(I+1)       
         V(I+2) = V(I+2)-W(I+2)       
         V(I+3) = V(I+3)-W(I+3)       
   30 CONTINUE    
      RETURN      
C       
      END 
      SUBROUTINE VEVPW (N,V,W)
C       
C ... VPW COMPUTES    V = V + W       
C       
C ... PARAMETER LIST:       
C       
C          N      LENGTH OF VECTORS V AND W     
C          V      D.P. VECTOR 
C          W      D.P. VECTOR SUCH THAT   V(I) = V(I) + W(I)
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N   
      DOUBLE PRECISION V(N),W(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,M,MP1       
C       
      IF (N.LE.0) RETURN    
C       
      M = MOD(N,4)
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         V(I) = V(I)+W(I)   
   10 CONTINUE    
      IF (N.LT.4) RETURN    
C       
   20 MP1 = M+1   
      DO 30 I = MP1,N,4     
         V(I) = V(I)+W(I)   
         V(I+1) = V(I+1)+W(I+1)       
         V(I+2) = V(I+2)+W(I+2)       
         V(I+3) = V(I+3)+W(I+3)       
   30 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE VFILL (N,V,VAL)      
C       
C     FILLS A VECTOR, V, WITH A CONSTANT VALUE, VAL.      
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTOR V    
C          V      D.P. VECTOR 
C          VAL    D.P. CONSTANT THAT FILLS FIRST N LOCATIONS OF V   
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N   
      DOUBLE PRECISION V(N),VAL       
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,M,MP1       
C       
      IF (N.LE.0) RETURN    
C       
C     CLEAN UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 10  
C       
      M = MOD(N,10) 
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         V(I) = VAL 
   10 CONTINUE    
      IF (N.LT.10) RETURN   
C       
   20 MP1 = M+1   
      DO 30 I = MP1,N,10    
         V(I) = VAL 
         V(I+1) = VAL       
         V(I+2) = VAL       
         V(I+3) = VAL       
         V(I+4) = VAL       
         V(I+5) = VAL       
         V(I+6) = VAL       
         V(I+7) = VAL       
         V(I+8) = VAL       
         V(I+9) = VAL       
   30 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE VOUT (N,V,ISWT,NOUTT)
C       
C     THIS SUBROUTINE EFFECTS PRINTING OF RESIDUAL AND SOLUTION     
C     VECTORS - CALLED FROM PERR    
C       
C ... PARAMETER LIST:       
C       
C          V      VECTOR OF LENGTH N  
C          ISWT   LABELLING INFORMATION 
C          NOUT OUTPUT DEVICE NUMBER (= NOUTT)  
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N,ISWT,NOUTT  
      DOUBLE PRECISION V(N) 
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,J,JM1,K,KUPPER,NOUT   
C       
      NOUT = NOUTT
C       
C        IF (N .LE. 0) RETURN 
C       
      KUPPER = MIN0(N,8)    
      IF (ISWT.EQ.1) WRITE (NOUT,10)  
   10 FORMAT (//5X,'RESIDUAL VECTOR') 
      IF (ISWT.EQ.2) WRITE (NOUT,20)  
   20 FORMAT (//5X,'SOLUTION VECTOR') 
      WRITE (NOUT,30) (I,I=1,KUPPER)  
   30 FORMAT (10X,8I15)     
      WRITE (NOUT,40)       
   40 FORMAT (10X,120('-')/)
C       
      DO 60 J = 1,N,8       
         KUPPER = MIN0(J+7,N) 
         JM1 = J-1
         WRITE (NOUT,50) JM1,(V(K),K=J,KUPPER)  
   50    FORMAT (4X,I5,'+  ',8D15.5)  
   60 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE WEVMW (N,V,W)
C       
C ... WEVMW COMPUTES W = V - W
C       
C ... PARAMETER LIST:       
C       
C          N      INTEGER LENGTH OF VECTORS V AND W       
C          V      D.P. VECTOR 
C          W      D.P. VECTOR SUCH THAT   W(I) = V(I) - W(I)
C       
C ... SPECIFICATIONS FOR ARGUMENTS    
C       
      INTEGER N   
      DOUBLE PRECISION V(N),W(N)      
C       
C ... SPECIFICATIONS FOR LOCAL VARIABLES
C       
      INTEGER I,M,MP1       
C       
      IF (N.LE.0) RETURN    
      M = MOD(N,4)
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         W(I) = V(I)-W(I)   
   10 CONTINUE    
      IF (N.LT.4) RETURN    
C       
   20 MP1 = M+1   
      DO 30 I = MP1,N,4     
         W(I) = V(I)-W(I)   
         W(I+1) = V(I+1)-W(I+1)       
         W(I+2) = V(I+2)-W(I+2)       
         W(I+3) = V(I+3)-W(I+3)       
   30 CONTINUE    
C       
      RETURN      
      END 
      SUBROUTINE ZBRENT (N,TRI,EPS,NSIG,AA,BB,MAXFNN,IER) 
C       
C   MODIFIED IMSL ROUTINE NAME   - ZBRENT       
C       
C-----------------------------------------------------------------------
C       
C   COMPUTER            - CDC/SINGLE  
C       
C   LATEST REVISION     - JANUARY 1, 1978       
C       
C   PURPOSE             - ZERO OF A FUNCTION WHICH CHANGES SIGN IN A
C                           GIVEN INTERVAL (BRENT ALGORITHM)
C       
C   USAGE               - CALL ZBRENT (F,EPS,NSIG,A,B,MAXFN,IER)    
C       
C   ARGUMENTS    TRI    - A TRIDIAGONAL MATRIX OF ORDER N 
C                EPS    - FIRST CONVERGENCE CRITERION (INPUT).  A ROOT, 
C                           B, IS ACCEPTED IF DABS(F(B)) IS LESS THAN OR
C                           EQUAL TO EPS.  EPS MAY BE SET TO ZERO.  
C                NSIG   - SECOND CONVERGENCE CRITERION (INPUT).  A ROOT,
C                           B, IS ACCEPTED IF THE CURRENT APPROXIMATION 
C                           AGREES WITH THE TRUE SOLUTION TO NSIG   
C                           SIGNIFICANT DIGITS. 
C                A,B    - ON INPUT, THE USER MUST SUPPLY TWO POINTS, A
C                           AND B, SUCH THAT F(A) AND F(B) ARE OPPOSITE 
C                           IN SIGN. (= AA, BB) 
C                           ON OUTPUT, BOTH A AND B ARE ALTERED.  B 
C                           WILL CONTAIN THE BEST APPROXIMATION TO THE
C                           ROOT OF F. SEE REMARK 1.      
C                MAXFN  - ON INPUT, MAXFN SHOULD CONTAIN AN UPPER BOUND 
C                           ON THE NUMBER OF FUNCTION EVALUATIONS   
C                           REQUIRED FOR CONVERGENCE.  ON OUTPUT, MAXFN 
C                           WILL CONTAIN THE ACTUAL NUMBER OF FUNCTION
C                           EVALUATIONS USED. (= MAXFNN)  
C                IER    - ERROR PARAMETER. (OUTPUT)       
C                         TERMINAL ERROR
C                           IER = 501 INDICATES THE ALGORITHM FAILED TO 
C                             CONVERGE IN MAXFN EVALUATIONS.
C                           IER = 502 INDICATES F(A) AND F(B) HAVE THE
C                             SAME SIGN.
C       
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32 
C                       - SINGLE/H36,H48,H60    
C       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND       
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL  
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C       
C   REMARKS  1.  LET F(X) BE THE CHARACTERISTIC FUNCTION OF THE MATRIX
C                TRI EVALUATED AT X. FUNCTION DETERM EVALUATES F(X).
C                ON EXIT FROM ZBRENT, WHEN IER=0, A AND B SATISFY THE 
C                FOLLOWING, 
C                F(A)*F(B) .LE.0,     
C                DABS(F(B)) .LE. DABS(F(A)), AND
C                EITHER DABS(F(B)) .LE. EPS OR  
C                DABS(A-B) .LE. MAX(DABS(B),0.1)*10.0**(-NSIG).     
C                THE PRESENCE OF 0.1 IN THIS ERROR CRITERION CAUSES 
C                LEADING ZEROES TO THE RIGHT OF THE DECIMAL POINT TO BE 
C                COUNTED AS SIGNIFICANT DIGITS. SCALING MAY BE REQUIRED 
C                IN ORDER TO ACCURATELY DETERMINE A ZERO OF SMALL   
C                MAGNITUDE. 
C            2.  ZBRENT IS GUARANTEED TO REACH CONVERGENCE WITHIN   
C                K = (DLOG((B-A)/D)+1.0)**2 FUNCTION EVALUATIONS WHERE
C                  D=MIN(OVER X IN (A,B) OF     
C                    MAX(DABS(X),0.1)*10.0**(-NSIG)).     
C                THIS IS AN UPPER BOUND ON THE NUMBER OF EVALUATIONS. 
C                RARELY DOES THE ACTUAL NUMBER OF EVALUATIONS USED BY 
C                ZBRENT EXCEED DSQRT(K). D CAN BE COMPUTED AS FOLLOWS,
C                  P = DBLE(AMIN1(DABS(A),DABS(B)))       
C                  P = DMAX1(0.1,P)   
C                  IF ((A-0.1)*(B-0.1).LT.0.0) P = 0.1    
C                  D = P*10.0**(-NSIG)
C       
C   COPYRIGHT           - 1977 BY IMSL, INC. ALL RIGHTS RESERVED.   
C       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.    
C       
C-----------------------------------------------------------------------
C       
C *** BEGIN: ITPACK COMMON  
C       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
C       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
C       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, 
     *   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
C       
C *** END  : ITPACK COMMON  
C       
C     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
C       
C                                  SPECIFICATIONS FOR ARGUMENTS     
C       
      INTEGER NSIG,MAXFNN,IER 
      DOUBLE PRECISION TRI(2,1),EPS,AA,BB       
C       
C                                  SPECIFICATIONS FOR LOCAL VARIABLES 
C       
      INTEGER IC,MAXFN      
      DOUBLE PRECISION ZERO,HALF,ONE,THREE,TEN,A,B,T,FA,FB,C,FC,D,E,TOL,
     *   RM,S,P,Q,R,RONE,TEMP,DETERM  
      DATA ZERO / 0.D0 / ,HALF / 5.D-1 / ,ONE / 1.D0 / ,THREE / 3.D0 / ,
     *   TEN / 10.D0 /      
C       
C                                  FIRST EXECUTABLE STATEMENT       
C       
      A = AA      
      B = BB      
      MAXFN = MAXFNN
      IER = 0     
      T = TEN**(-NSIG)      
      IC = 2      
      FA = DETERM(N,TRI,A)  
      FB = DETERM(N,TRI,B)  
      S = B       
C       
C                                  TEST FOR SAME SIGN     
C       
      IF (FA*FB.GT.ZERO) GO TO 110    
   10 C = A       
      FC = FA     
      D = B-C     
      E = D       
   20 IF (DABS(FC).GE.DABS(FB)) GO TO 30
      A = B       
      B = C       
      C = A       
      FA = FB     
      FB = FC     
      FC = FA     
   30 CONTINUE    
      TOL = T*DMAX1(DABS(B),0.1D0)    
      RM = (C-B)*HALF       
C       
C                                  TEST FOR FIRST CONVERGENCE CRITERIA
C       
      IF (DABS(FB).LE.EPS) GO TO 80   
C       
C                                  TEST FOR SECOND CONVERGENCE CRITERIA 
C       
      IF (DABS(C-B).LE.TOL) GO TO 80  
C       
C                                  CHECK EVALUATION COUNTER 
C       
      IF (IC.GE.MAXFN) GO TO 90       
C       
C                                  IS BISECTION FORCED    
C       
      IF (DABS(E).LT.TOL) GO TO 60    
      IF (DABS(FA).LE.DABS(FB)) GO TO 60
      S = FB/FA   
      IF (A.NE.C) GO TO 40  
C       
C                                  LINEAR INTERPOLATION   
C       
      P = (C-B)*S 
      Q = ONE-S   
      GO TO 50    
C       
C                                  INVERSE QUADRATIC INTERPOLATION  
C       
   40 Q = FA/FC   
      R = FB/FC   
      RONE = R-ONE
      P = S*((C-B)*Q*(Q-R)-(B-A)*RONE)
      Q = (Q-ONE)*RONE*(S-ONE)
   50 IF (P.GT.ZERO) Q = -Q 
      IF (P.LT.ZERO) P = -P 
      S = E       
      E = D       
C       
C                                  IF DABS(P/Q).GE.75*DABS(C-B) THEN
C                                     FORCE BISECTION     
C       
      IF (P+P.GE.THREE*RM*Q) GO TO 60 
C       
C                                  IF DABS(P/Q).GE..5*DABS(S) THEN FORCE
C                                     BISECTION. S = THE VALUE OF P/Q 
C                                     ON THE STEP BEFORE THE LAST ONE 
C       
      IF (P+P.GE.DABS(S*Q)) GO TO 60  
      D = P/Q     
      GO TO 70    
C       
C                                  BISECTION    
C       
   60 E = RM      
      D = E       
C       
C                                  INCREMENT B  
C       
   70 A = B       
      FA = FB     
      TEMP = D    
      IF (DABS(TEMP).LE.HALF*TOL) TEMP = DSIGN(HALF*TOL,RM) 
      B = B+TEMP  
      S = B       
      FB = DETERM(N,TRI,S)  
      IC = IC+1   
      IF (FB*FC.LE.ZERO) GO TO 20     
      GO TO 10    
C       
C                                  CONVERGENCE OF B       
C       
   80 A = C       
      MAXFN = IC  
      GO TO 130   
C       
C                                  MAXFN EVALUATIONS      
C       
   90 IER = 501   
      A = C       
      MAXFN = IC  
      IF (LEVEL.GE.1) WRITE (NOUT,100) MAXFN    
  100 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE ZBRENT  '/' ',  
     *   '    ALGORITHM FAILED TO CONVERGE   '/' ','    IN',I6,     
     *   ' ITERATIONS ')    
      GO TO 130   
C       
C                                  TERMINAL ERROR - F(A) AND F(B) HAVE
C                                  THE SAME SIGN
C       
  110 IER = 502   
      MAXFN = IC  
      IF (LEVEL.GE.1) WRITE (NOUT,120)
  120 FORMAT ('0','*** W A R N I N G ************'/'0',   
     *   '    IN ITPACK ROUTINE ZBRENT  '/' ',  
     *   '    F(A) AND F(B) HAVE SAME SIGN   ') 
  130 CONTINUE    
      AA = A      
      BB = B      
      MAXFNN = MAXFN
      RETURN      
      END 
