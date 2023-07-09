C     Concrete Stress-Strain Model Including Unloading/Reloading Model
C     Proposed by Sakai and Kawashima (Envelope Curve: Mander Model)
C                                      Coded by J. Sakai    2002. 2.14
C ====================================================================
C ====================================================================
C ====================================================================
      SUBROUTINE NLU014( NTABL, MTB, NTB, D, MNLI, MNLR, NP, P,
     1                   NHST, HST, NIHST, IHST,
     1                   DEPS, DSIG, TK, DEPSV,
     1                   NEL, IPRT, INIT, LDBG, IT, ILOOP)
      IMPLICIT REAL*8(A-H,O-Z)
cDEC$ ATTRIBUTES DLLEXPORT:: NLU014, /MLSVAL/
      COMMON /MLSVAL/ VAL1, VAL2, VALS
      DIMENSION NTB(NTABL),D(MTB,NTABL),NP(MNLI),P(MNLR),
     1          HST(NHST),IHST(NIHST)
c  ---------------------------------------------------------------------
c     Material Type: Confined Concrete
c     Envelope Curve: Mander Model
c     Unloading and Reloading Curve: Sakai-Kawashima Model
c
c     Material Parameters from MATNL1 Data
c        P(1) : Initial Slope or Initial Modulus of Elasticity       YMc
c        P(2) : Peak Stress                                          Sigcc
c        P(3) : Strain at Peak Stress                                EPScc
c
c     State Variables
c        HST( 1) : Initial Stiffness
c        HST( 2) : Strain
c        HST( 3) : Stress
c        HST( 4) : Delta Epsilon of the Previous Step                  DE0
c        HST( 5) : Stress at the Unloading Point on Envelope           Sigule
c        HST( 6) : Strain at the Unloading Point on Envelope           EPSule
c        HST( 7) : Stress at the Unloading Point                       Sigul
c        HST( 8) : Strain at the Unloading Point                       EPSul
c        HST( 9) : Plastic Strain                                      EPSpl
c        HST(10) : Stress at the Unloading Point after Un/Reloading    Suln
c        HST(11) : Slope of Reloding Branch                            YMrl
c        HST(12) : Tangential Slope of the Previous Step               YMtan
c        HST(13) : Stress at the Reloading Point                       Sigrl
c        HST(14) : Strain at the Reloading Point                       EPSrl
c        HST(15) : HST( 9) EPSpl of Previous Unloading                 EPSpl0
c        HST(16) : HST(10) Suln of Previous Reloading                  Suln0
c        HST(17) : Partial Reloading Ratio GammaRL                     GamRL
c
c        IHST(1) : Index of Current Status                      Jcon
c        Jcon= 1 : On the Envelope Curve
c        Jcon= 2 : Zero Stress
c        Jcon= 3 : On the Unloading Curve (Parabolic)
c        Jcon= 4 : On the Reloading Line (Linear)   0.2< Normalized Strain <1.0
c        Jcon= 5 : On the Reloading Curve (Parabolic) 0< Normalized Strain <0.2
c        Jcon= 6 : On the Reloading Curve (Parabolic) 0.75< Normalized Stress <1
c        Jcon= 7 : On the Reloading Line (Linear)     Normalized Strain >1.0
c        Jcon= 8 : Zero Stress (No Compressive Loading)
c        IHST(2) : Number of Cyclcles                           Ncyc
c        IHST(3) : Index of Previous Status                     Jcon0
c        IHST(4) : Number of Cyclcles of Previous Status        Ncyc0
c
c        DEPS   (IN)   : Incremental Strain
c        DSIG   (OUT)  : Incremental Stress
c        TK     (OUT)  : Current Status
c  =====================================================================
c  ---------------------------------------------------------------------
c   Stress Tensor at the Start of This Increment
      Sig0=HST(3)
c  ---------------------------------------------------------------------
c   Reference Variables
      EPS0=HST(2)
      DE=DEPS
c  ---------------------------------------------------------------------
c   Strain at the End of This Increment
      EPS1=EPS0+DE
c  ---------------------------------------------------------------------
c   User Defined Material Parameters
      YMc=P(1)
      Sigcc=P(2)
      EPScc=P(3)
c
      r=YMc*EPScc/(YMc*EPScc-Sigcc)
c  ---------------------------------------------------------------------
c   Initial State Variables (INIT=1)
      IF(INIT.EQ.1) THEN
      NHST=17
      NIHST=4
      HST(1)=TK
c
      DO 10 I=4,NHST
      HST(I)=0.0
   10 CONTINUE
c
      IHST(1)=8
      IHST(2)=1
      IHST(3)=8
      IHST(4)=1
c
c  ---------------------------------------------------------------------
c   State Variables at the Start of This Increment
      ELSE
      DE0=HST(4)
      Sigule=HST(5)
      EPSule=HST(6)
      Sigul=HST(7)
      EPSul=HST(8)
      EPSpl=HST(9)
      Suln=HST(10)
      YMrl=HST(11)
      YMtan0=HST(12)
      Sigrl=HST(13)
      EPSrl=HST(14)
      EPSpl0=HST(15)
      Suln0=HST(16)
      GamRL=HST(17)
c
      Jcon=IHST(1)
      Ncyc=IHST(2)
      Jcon0=IHST(3)
      Ncyc0=IHST(4)
c
      r=YMc*EPScc/(YMc*EPScc-Sigcc)
C ********************************************************************
      IF(DE.EQ.0.0) THEN
      Sig1=Sig0
      YMtan=YMtan0
      ELSE IF(DE.LT.0.0) THEN
c   Compressive Loading
      CALL COMPR14 (EPS1,Sig1,YMtan,Jcon,EPS0,Sig0,YMc,EPScc,Sigcc,r,
     1              EPSule,Sigule,EPSpl,Suln,YMrl,EPSul,Sigul,
     1              EPSrl,Sigrl,EPSpl0,Suln0,GamRL,
     1              Ncyc,Jcon0,Ncyc0)
      ELSE
c   Unloading (Tensile) Loading
      CALL TENSI14 (EPS1,Sig1,YMtan,Jcon,EPS0,Sig0,
     1              EPSule,Sigule,EPSpl,Suln,EPSul,Sigul,
     1              EPSrl,Sigrl,EPSpl0,Suln0,GamRL,
     1              Ncyc,Jcon0,Ncyc0)
      END IF
C ********************************************************************
c   New State Variables
      HST(2)=EPS1
      HST(3)=Sig1
      DSIG=Sig1-Sig0
      TK=YMtan
      HST(4)=DE
      HST(5)=Sigule
      HST(6)=EPSule
      HST(7)=Sigul
      HST(8)=EPSul
      HST(9)=EPSpl
      HST(10)=Suln
      HST(11)=YMrl
      HST(12)=YMtan
      HST(13)=Sigrl
      HST(14)=EPSrl
      HST(15)=EPSpl0
      HST(16)=Suln0
      HST(17)=GamRL
c
      IHST(1)=Jcon
      IHST(2)=Ncyc
      IHST(3)=Jcon0
      IHST(4)=Ncyc0
c  ---------------------------------------------------------------------
      END IF
      RETURN
      END
c
C ====================================================================
C ------------------------ SUBROUTINE PROGRAM ------------------------
C ------------------- Steel Stress-Strain Relation   -----------------
C ========================== COMPRESSION =============================
C ====================================================================
C
      SUBROUTINE COMPR14 (EPS1,Sig1,YMtan,Jcon,EPS0,Sig0,YMc,EPScc,
     1                    Sigcc,r,EPSule,Sigule,EPSpl,Suln,YMrl,EPSul,
     1                    Sigul,EPSrl,Sigrl,EPSpl0,Suln0,GamRL,
     1                    Ncyc,Jcon0,Ncyc0)
      IMPLICIT REAL*8(A-H,O-Z)
c
C *** ENVELOPE CURVE (Jcon=1) ****************************************
      IF(Jcon.EQ.1) THEN
      GO TO 20
C
C *** ZERO STRESS (Jcon=2) *******************************************
      ELSE IF(Jcon.EQ.2) THEN
        IF(EPS1.LT.EPSpl) THEN
        Jcon0=2
C
        YMrl=(Suln-0.1*Suln0)/(EPSule-(0.2*EPSule+0.8*EPSpl))
        Sigrl=0.0
        EPSrl=EPSpl
          IF(EPS1.LT.EPSule) THEN
          Jcon=7
          GO TO 10
          ELSE IF(EPS1.LE.0.2*EPSule+0.8*EPSpl) THEN
          Jcon=4
          Sig1=YMrl*(EPS1-EPSule)+Suln
          YMtan=YMrl
          ELSE
          Jcon=5
          Sig1=2.5*Suln0*((EPS1-EPSpl)/(EPSule-EPSpl))**2.0
          YMtan=5.0*Suln0*(EPS1-EPSpl)/((EPSule-EPSpl))**2.0
          END IF
        ELSE
C
C  When hysteresis does not reach the plastic strain
        Sig1=0.0
        YMtan=0.0
        END IF
C
C *** RELOADING FROM UNLOADING PATH (Jcon=3) *************************
      ELSE IF(Jcon.EQ.3) THEN
C
C  When partial unloading ratio is smaller than 0.25
        IF(ABS(Sig0-Sigul).LT.ABS(Sigule*0.25)) THEN
          IF(EPS1.LE.EPSul) THEN
C
C  When strain has already exceeded the previous unloaded strain
            IF(Jcon0.EQ.7) THEN
            Jcon=7
            Jcon0=3
              IF(Sigrl.EQ.0.0) THEN
              IF(GamRL.GE.0.7) THEN
              EPSpl=EPSpl0
              Suln=Suln0
              Ncyc=Ncyc0
              END IF
              END IF
            GO TO 10
            ELSE IF(Jcon0.EQ.4) THEN
            Jcon=4
            Jcon0=3
              IF(Sigrl.EQ.0.0) THEN
              IF(GamRL.GE.0.7) THEN
              EPSpl=EPSpl0
              Suln=Suln0
              Ncyc=Ncyc0
              END IF
              END IF
            GO TO 5
            ELSE 
            Jcon=1
            Jcon0=3
            GO TO 20
            END IF
          ELSE
C
C  Otherwise, reloading along the previous unloading path
C  Jcon0 in not renewed here.
          Jcon=6
          Sig1=Sigul*((EPS1-EPSpl)/(EPSul-EPSpl))**2.0
          YMtan=2*Sigul*(EPS1-EPSpl)/((EPSul-EPSpl)**2.0)
          GO TO 30
          END IF
C
C  When partial unloading ratio is larger than 0.25
        ELSE
        Sigrl=Sig0
        EPSrl=EPS0
        YMrl=(Suln-Sigrl)/(EPSule-EPSrl)
        Jcon=4
        Jcon0=3
        GO TO 5
        END IF
C
C *** ON RELOADING PATH (Jcon=4) *************************************
      ELSE IF(Jcon.EQ.4) THEN
      GO TO 5
C
C *** ON RELOADING PATH, 0< Normalized Strain <0.2 (Jcon=5) **********
      ELSE IF(Jcon.EQ.5) THEN
        IF(EPS1.LT.EPSule) THEN
        Jcon=7
        GO TO 10
        ELSE IF(EPS1.LE.0.2*EPSule+0.8*EPSpl) THEN
        Jcon=4
        Sig1=YMrl*(EPS1-EPSule)+Suln
        YMtan=YMrl
        ELSE
        Sig1=2.5*Suln0*((EPS1-EPSpl)/(EPSule-EPSpl))**2.0
        YMtan=5.0*Suln0*(EPS1-EPSpl)/((EPSule-EPSpl)**2.0)
        END IF
C
C *** RELOADING ALONG THE PREVIOUS UNLOADING PATH (Jcon=6) ***********
      ELSE IF(Jcon.EQ.6) THEN
C
C  When strain has already exceeded the previous unloaded strain
        IF(EPS1.LE.EPSul) THEN
          IF(Jcon0.EQ.7) THEN
          Jcon=7
          Jcon0=3
            IF(Sigrl.EQ.0.0) THEN
            IF(GamRL.GE.0.7) THEN
            EPSpl=EPSpl0
            Suln=Suln0
            Ncyc=Ncyc0
            END IF
            END IF
          GO TO 10
          ELSE IF(Jcon0.EQ.4) THEN
          Jcon=4
          Jcon0=3
            IF(Sigrl.EQ.0.0) THEN
            IF(GamRL.GE.0.7) THEN
            EPSpl=EPSpl0
            Suln=Suln0
            Ncyc=Ncyc0
            END IF
            END IF
          GO TO 5
          ELSE 
          Jcon=1
          Jcon0=3
          GO TO 20
          END IF
C
C  Otherwise, still reloading along the previous unloading path
        ELSE
        Sig1=Sigul*((EPS1-EPSpl)/(EPSul-EPSpl))**2.0
        YMtan=2*Sigul*(EPS1-EPSpl)/((EPSul-EPSpl)**2.0)
        END IF
C
C *** RELOADING (NORMALIZED STRAIN IS LARGER THAN 1) (Jcon=7) *****
C *** Hysteresis has exceeded the unloaded strain from the envelope curve
C  Jcon0 in not renewed here. (Jcon0=3)
      ELSE IF(Jcon.EQ.7) THEN
      GO TO 10
C
C *** VIRGIN CONCRETE (Not Loaded in Compression) (Jcon=8) ***********
      ELSE
        IF(EPS1.GE.0.0) THEN
        Sig1=0.0
        YMtan=0.0
        ELSE
        Jcon=1
        Jcon0=8
        GO TO 20
        END IF
C --------------------------------------------------------------------
      END IF
      GO TO 30
C
C ********************************************************************
C *** Reloading Path (Linear Function, 0.2< Normalized Strain <1)
    5 CONTINUE
      IF(EPS1.LT.EPSule) THEN
      Jcon=7
      GO TO 10
      ELSE
      Sig1=YMrl*(EPS1-EPSule)+Suln
      YMtan=YMrl
      END IF
      GO TO 30
C
C ********************************************************************
C *** When Jcon=7
   10 CONTINUE
C
C *** Stress on the reloading path
      Sign0=YMrl*(EPS1-EPSule)+Suln
      YMtn0=YMrl
C
C *** Stress on the envelope curve (Mander Model)
      EPSxx=EPS1/EPScc
      AAA=EPSxx**r
      Senv=Sigcc*EPSxx*r/(r-1.0+AAA)
      YMtenv=Sigcc/EPScc*r*(r-1.0)*(1.0-AAA)/(r-1.0+AAA)**2.0
C
C *** If hysteresis has reached the envelope curve or not
      IF(Sign0.GT.Senv) THEN
      Sig1=Sign0
      YMtan=YMtn0
      ELSE
      Jcon=1
      Jcon0=7
      Sig1=Senv
      YMtan=YMtenv
      END IF
      GO TO 30
C
C *** ENVELOPE CURVE (Mander Model)
   20 CONTINUE
      EPSxx=EPS1/EPScc
      AAA=EPSxx**r
      Sig1=Sigcc*EPSxx*r/(r-1.0+AAA)
      YMtan=Sigcc/EPScc*r*(r-1.0)*(1.0-AAA)/(r-1.0+AAA)**2.0
      GO TO 30
C --------------------------------------------------------------------
C
   30 CONTINUE
      RETURN
      END
C ====================================================================
C
C ====================================================================
C ----------------------- SUBROUTINE PROGRAM -------------------------
C ======================= UNLOADING (TENSION) ========================
C ====================================================================
C
      SUBROUTINE TENSI14 (EPS1,Sig1,YMtan,Jcon,EPS0,Sig0,
     1                    EPSule,Sigule,EPSpl,Suln,EPSul,Sigul,
     1                    EPSrl,Sigrl,EPSpl0,Suln0,GamRL,
     1                    Ncyc,Jcon0,Ncyc0)
      IMPLICIT REAL*8(A-H,O-Z)
C
C *** UNLOADING FROM THE ENVELOPE CURVE (Jcon=1) *********************
      IF(Jcon.EQ.1) THEN
      Ncyc=1
      Jcon=3
      Jcon0=1
      GO TO 10
C
C *** ZERO STRESS (Jcon=2) *******************************************
      ELSE IF(Jcon.EQ.2) THEN
      Sig1=0.0
      YMtan=0.0
C
C *** ON UNLOADING PATH (Jcon=3) *************************************
      ELSE IF(Jcon.EQ.3) THEN
        IF(EPS1.GE.EPSpl) THEN
        Jcon=2
        Jcon0=3
        Sig1=0.0
        YMtan=0.0
        ELSE
        Sig1=Sigul*((EPS1-EPSpl)/(EPSul-EPSpl))**2.0
        YMtan=2.0*Sigul*(EPS1-EPSpl)/((EPSul-EPSpl)**2.0)
        END IF
C
C *** UNLOADING FROM RELOADING PATH (Jcon=4) *************************
      ELSE  IF(Jcon.EQ.4) THEN
      Jcon=3
      Jcon0=4
      GO TO 20
C
C *** ON UNLOADING PATH, 0< Normalized Strain <0.2 (Jcon=5) **********
      ELSE IF(Jcon.EQ.5) THEN
        IF(EPS1.GE.EPSpl) THEN
        Jcon=2
        Sig1=0.0
        YMtan=0.0
        ELSE
        Sig1=2.5*Suln0*((EPS1-EPSpl)/(EPSule-EPSpl))**2.0
        YMtan=5.0*Suln0*(EPS1-EPSpl)/((EPSule-EPSpl)**2.0)
        END IF
C
C *** UNLOADING FROM RELOADING ALONG THE PREVIOUS UNLOADING PATH (Jcon=6)
      ELSE IF(Jcon.EQ.6) THEN
C
        IF(EPS1.GE.EPSpl) THEN
        Jcon0=3
        Jcon=2
        Sig1=0.0
        YMtan=0.0
        ELSE
C  Jcon0 in not renewed here.
        Jcon=3
        Sig1=Sigul*((EPS1-EPSpl)/(EPSul-EPSpl))**2.0
        YMtan=2.0*Sigul*(EPS1-EPSpl)/((EPSul-EPSpl)**2.0)
        END IF
C
C *** UNLOADING FROM WHERE NORMALIZED STRAIN IS LARGER THAN 1 (Jcon=7)
      ELSE IF(Jcon.EQ.7) THEN
C *** (EPSule, Sigule) is renewed.
      EPSule=EPS0
      Sigule=Sig0
C
      Suln=Sig0
      Jcon=3
      Jcon0=7
      GO TO 20
C
C *** VIRGIN CONCRETE (Not Loaded in Compression) (Jcon=8) ***********
      ELSE
      Sig1=0.0
      YMtan=0.0
C --------------------------------------------------------------------
      END IF
      GO TO 30
C
C ********************************************************************
C *** UNLOADING FROM THE ENVELOPE CURVE
   10 CONTINUE
C *** Initialize the data about unloading point
      EPSule=EPS0
      Sigule=Sig0
C
      EPSul=EPSule
      Sigul=Sigule
C
      EPSrl=EPSule
      Sigrl=Sigule
C
C *** Prediction of the Plastic Strain from Eq.(9)
        IF(EPSule.GT.-1000.0e-6) THEN
        EPSpl=0.0
        ELSE IF(EPSule.GT.-3500.0e-6) THEN
        EPSpl=0.43*(EPSule+0.001)
        ELSE
        EPSpl=0.94*(EPSule+0.00235)
        END IF
C
        EPSpl0=EPSpl
C
C
C *** Prediction of the Unloaded Stress from Eqs.(18) and(20)
        IF(EPSule.GT.-1000.0e-6) THEN
        Suln=Sigule*1.0
        ELSE IF(EPSule.GT.-3500.0e-6) THEN
        Suln=Sigule*(1.0+32.0*(EPSule+0.001))
        ELSE
        Suln=Sigule*0.92
        END IF
C
        Suln0=Sigul
        Ncyc0=Ncyc
        GamRL=0.0
C
C *** Already reached the plastic strain
        IF(EPS1.GE.EPSpl) THEN
        Jcon=2
        Jcon0=3
        Sig1=0.0
        YMtan=0.0
C *** Otherwise, on unloading path
        ELSE
        Sig1=Sigul*((EPS1-EPSpl)/(EPSul-EPSpl))**2.0
        YMtan=2.0*Sigul*(EPS1-EPSpl)/((EPSul-EPSpl)**2.0)
        END IF
      GO TO 30
C
C ********************************************************************
C *** UNLOADING FROM RELOADING PATH
   20 CONTINUE
C *** Initialize the data about unloaded point
      EPSul=EPS0
      Sigul=Sig0
      GamRL=(EPS0-EPSpl)/(EPSule-EPSpl)
C
C ======= If it satisfies Eq. (6), the number of cycles is counted.
        IF(Sigrl.EQ.0.0) THEN
        IF(GamRL.GE.0.7) THEN
        Ncyc0=Ncyc
        Ncyc=Ncyc+1
C
          IF(GamRL.GT.1.0) THEN
          GamRL=1.0
          END IF
C
C *** Prediction of the Plastic Strain from Eqs.(11), (12) and (13)
        EPSpl0=EPSpl
          IF(Ncyc.LE.2) THEN
          Gamma=0.945+0.2*(1.0-GamRL)
          ELSE
          Gamma=0.965+0.005*(Ncyc-3)+0.2*(1.0-GamRL)
          END IF
C
          IF(Gamma.GT.1.0) THEN
          Gamma=1.0
          END IF
C
          IF(EPSule.GT.-1000.0e-6) THEN
          Gamma=1.0
          END IF
C
        EPSpl=(1.0-Gamma)*EPSule+Gamma*EPSpl
C
C
C *** Prediction of the Unloaded Stress from Eqs.(18) ~ (21)
          Suln0=Suln
C
          IF(Ncyc.LE.2) THEN
            IF(EPSule.GT.-3500.0e-6) THEN
            aBeta=1.0+(42.0-10.0*Ncyc)*(EPSule+0.001)
            ELSE
            aBeta=0.945
            END IF
          ELSE
            IF(EPSule.GT.-3500.0e-6) THEN
            aBeta=1.0+(20.0-2.0*Ncyc)*(EPSule+0.001)
            ELSE
            aBeta=0.965+0.005*(Ncyc-3)
            END IF
          END IF
          Beta=aBeta+0.2*(1.0-GamRL)
C
          IF(Beta.GT.1.0) THEN
          Beta=1.0
          END IF
C
          Suln=Suln*Beta
C
C ======= If it does not satisfy Eq. (6),
C         and if the partial reloading ratio is less than 0.7, then
        ELSE
        Suln=Suln
        EPSpl=EPSpl
        END IF
        GO TO 22
C
C         If the reloading is not from zero stress, then
        ELSE
        Suln=Suln
        EPSpl=EPSpl
        END IF
C
   22   IF(EPS1.GE.EPSpl) THEN
        Jcon=2
        Jcon0=3
        Sig1=0.0
        YMtan=0.0
        ELSE
        Sig1=Sigul*((EPS1-EPSpl)/(EPSul-EPSpl))**2.0
        YMtan=2.0*Sigul*(EPS1-EPSpl)/((EPSul-EPSpl)**2.0)
        END IF
      GO TO 30
C --------------------------------------------------------------------
   30 CONTINUE
      RETURN
      END
C ====================================================================
C ====================================================================
C ====================================================================
