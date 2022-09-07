C
C  This file is part of MUMPS 5.4.1, released
C  on Tue Aug  3 09:49:43 UTC 2021
C
C
C  Copyright 1991-2021 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
C
C History:
C -------
C This maximum transversal set of routines are
C based on the work done by Jacko Koster at CERFACS for 
C his PhD thesis from Institut National Polytechnique de Toulouse 
C at CERFACS (1995-1997) and includes modifications provided 
C by the author as well as work done by Stephane Pralet 
C first at CERFACS during his PhD thesis (2003-2004) then
C at INPT-IRIT (2004-2005) during his post-doctoral position.
C
C The main research publication references for this work are:
C  [1] I. S. Duff, (1981),
C      "Algorithm 575. Permutations for a zero-free diagonal",
C      ACM Trans. Math. Software 7(3), 387-390.
C  [2] I. S. Duff and J. Koster, (1998),
C      "The design and use of algorithms for permuting large
C      entries to the diagonal of sparse matrices",
C      SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
C  [3] I. S. Duff and J. Koster, (2001),
C      "On algorithms for permuting large entries to the diagonal 
C      of sparse matrices",
C      SIAM J. Matrix Anal. Appl., vol. 22, no. 4, pp. 973-996.
C
      SUBROUTINE DMUMPS_MTRANSI(ICNTL,CNTL)
      IMPLICIT NONE
      INTEGER NICNTL, NCNTL
      PARAMETER (NICNTL=10, NCNTL=10)
      INTEGER ICNTL(NICNTL)
      DOUBLE PRECISION CNTL(NCNTL)
      INTEGER I
      ICNTL(1) =  6
      ICNTL(2) =  6
      ICNTL(3) = -1
      ICNTL(4) = -1
      ICNTL(5) =  0
      DO 10 I = 6,NICNTL
        ICNTL(I) = 0
   10 CONTINUE
      CNTL(1) = 0.0D0
      CNTL(2) = 0.0D0
      DO 20 I = 3,NCNTL
        CNTL(I) = 0.0D0
   20 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_MTRANSI
      SUBROUTINE DMUMPS_MTRANSB
     &           (M,N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D,RINF) 
      IMPLICIT NONE
      INTEGER ::  M,N,NUM
      INTEGER(8), INTENT(IN)  :: NE
      INTEGER :: IRN(NE),IPERM(M),JPERM(N),Q(M),L(M)
      INTEGER(8), INTENT(IN)  :: IP(N+1)
      INTEGER(8), INTENT(OUT) :: PR(N)
      DOUBLE PRECISION  :: A(NE)
      DOUBLE PRECISION  :: D(M), RINF
      INTEGER    :: I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &              I0,UP,LOW, IK
      INTEGER(8) :: K,KK,KK1,KK2
      DOUBLE PRECISION    CSP,DI,DNEW,DQ0,AI,A0,BV,TBV,RLX
      DOUBLE PRECISION    ZERO,MINONE,ONE
      PARAMETER (ZERO=0.0D0,MINONE=-1.0D0,ONE=1.0D0)
      INTRINSIC abs,min
      EXTERNAL DMUMPS_MTRANSD, DMUMPS_MTRANSE,
     &         DMUMPS_MTRANSF, DMUMPS_MTRANSX
      RLX = D(1)
      NUM = 0
      BV = RINF
      DO 10 I = 1,N
        JPERM(I) = 0
        PR(I) = IP(I)
   10 CONTINUE
      DO 12 I = 1,M
        IPERM(I) = 0
        D(I) = ZERO
   12 CONTINUE
      DO 30 J = 1,N
        A0 = MINONE
        DO 20 K = IP(J),IP(J+1)-1_8
          I = IRN(K)
          AI = abs(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 20
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 20
            JPERM(J) = I 
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 20
            A0 = AI
            I0 = I
          ENDIF
   20   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 30
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   30 CONTINUE
      IF (M.EQ.N) THEN
        DO 35 I = 1,M
          BV = min(BV,D(I))
   35   CONTINUE
      ENDIF
      IF (NUM.EQ.N) GO TO 1000
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1_8
          I = IRN(K)
          AI = abs(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1_8
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (abs(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1_8
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1_8
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1_8
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,M
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE
      TBV = BV * (ONE-RLX)
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = M + 1
        UP = M + 1
        CSP = MINONE
        J = JORD
        PR(J) = -1_8
        DO 115 K = IP(J),IP(J+1)-1_8
          I = IRN(K)
          DNEW = abs(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.TBV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.TBV) THEN
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL DMUMPS_MTRANSD(I,M,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = int(J,8)
          ENDIF
  115   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            TBV = BV * (ONE-RLX)
            DO 152 IDUM = 1,M
              CALL DMUMPS_MTRANSE(QLEN,M,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).LT.TBV) GO TO 153
  152       CONTINUE
          ENDIF
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1_8
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
            DNEW = min(DQ0,abs(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.TBV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.TBV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.TBV) THEN
                IF (DI.NE.MINONE) THEN
                  CALL DMUMPS_MTRANSF(L(I),QLEN,M,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL DMUMPS_MTRANSD(I,M,Q,D,L,1)
              ENDIF
              JJ = IPERM(I)
              PR(JJ) = int(J,8)
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.MINONE) GO TO 190
        BV = min(BV,CSP)
        TBV = BV * (ONE-RLX)
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = int(PR(J))
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
  190   DO 191 IK = UP,M
          I = Q(IK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE 
        DO 192 IK = LOW,UP-1
          I = Q(IK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 IK = 1,QLEN
          I = Q(IK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
 1000 IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000
      CALL DMUMPS_MTRANSX(M,N,IPERM,L,JPERM)
 2000 RETURN
      END SUBROUTINE DMUMPS_MTRANSB
      SUBROUTINE DMUMPS_MTRANSD(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DI
      POS = L(I)
      IF (POS.LE.1) GO TO 20
      DI = D(I)
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END SUBROUTINE DMUMPS_MTRANSD
      SUBROUTINE DMUMPS_MTRANSE(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END SUBROUTINE DMUMPS_MTRANSE
      SUBROUTINE DMUMPS_MTRANSF(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        IF (POS.LE.1) GO TO 20
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
   20   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
      ELSE
        IF (POS.LE.1) GO TO 34
        DO 32 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 34
   32   CONTINUE
   34   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
      ENDIF
   40 Q(POS) = I
      L(I) = POS
      RETURN
      END SUBROUTINE DMUMPS_MTRANSF
      SUBROUTINE DMUMPS_MTRANSQ(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
      INTEGER ::WLEN,NVAL
      INTEGER :: LENL(*),LENH(*),W(*)
      INTEGER(8) :: IP(*)
      DOUBLE PRECISION :: A(*),VAL
      INTEGER XX,J,K,S,POS
      INTEGER(8) :: II
      PARAMETER (XX=10)
      DOUBLE PRECISION SPLIT(XX),HA
      NVAL = 0 
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+int(LENL(J),8),IP(J)+int(LENH(J)-1,8)
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)
      RETURN
      END SUBROUTINE DMUMPS_MTRANSQ
      SUBROUTINE DMUMPS_MTRANSR(N,NE,IP,IRN,A)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N
      INTEGER(8), INTENT(IN) :: NE
      INTEGER(8), INTENT(IN) :: IP(N+1)
      INTEGER, INTENT(INOUT) :: IRN(NE)
      DOUBLE PRECISION, INTENT(INOUT)    :: A(NE)
      INTEGER  :: THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
      INTEGER    :: J, LEN, HI
      INTEGER(8) :: K, IPJ, TD, FIRST, LAST, MID, R, S
      DOUBLE PRECISION       :: HA, KEY
      INTEGER(8) :: TODO(TDLEN)
      DO 100 J = 1,N
        LEN = int(IP(J+1) - IP(J))
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ +int(LEN,8)
        TD = 2_8
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
        TD = TD - 2_8
        GO TO 425
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2_8
  425   CONTINUE
        IF (TD.EQ.0_8) GO TO 400 
        IF (TODO(TD)-TODO(TD-1).GE.int(THRESH,8)) GO TO 500
        TD = TD - 2_8
        GO TO 425
  400   DO 200 R = IPJ+1_8,IPJ+int(LEN-1,8)
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1_8)
            IRN(R) = IRN(R-1_8)
            DO 300 S = R-1,IPJ+1_8,-1_8
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200 
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_MTRANSR
      SUBROUTINE DMUMPS_MTRANSS(M,N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4,RLX,RINF)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER(8), INTENT(IN) :: NE
      INTEGER, INTENT(OUT) :: NUMX
      INTEGER(8), INTENT(IN) :: IP(N+1)
      INTEGER  :: IRN(NE),IPERM(N), 
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(M),IW4(3*N+M)
      DOUBLE PRECISION A(NE),RLX,RINF
      INTEGER :: NUM,NVAL,WLEN,I,J,L,CNT,MOD, IDUM
      INTEGER(8) :: K, II, KDUM1, KDUM2
      DOUBLE PRECISION ::    BVAL,BMIN,BMAX
      EXTERNAL DMUMPS_MTRANSQ,DMUMPS_MTRANSU,DMUMPS_MTRANSX
      DO 20 J = 1,N
        FC(J) = J
        LEN(J) = int(IP(J+1) - IP(J))
   20 CONTINUE
      DO 21 I = 1,M
        IW(I) = 0
   21 CONTINUE
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL DMUMPS_MTRANSU(CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,
     &            NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(2*N+M+1))
      NUM = NUMX
      IF (NUM.NE.N) THEN
        BMAX = RINF
      ELSE
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0D0
          DO 25 K = IP(J),IP(J+1)-1_8
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001D0 * BMAX
      ENDIF
      BVAL = 0.0D0
      BMIN = 0.0D0
      WLEN = 0
      DO 48 J = 1,N
        L = int(IP(J+1) - IP(J))
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1_8
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
        K = IP(J+1)
   46   LENL(J) = int(K - IP(J))
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE
      DO 90 KDUM1 = 1_8,NE
        IF (NUM.EQ.NUMX) THEN
          DO 50 I = 1,M
            IPERM(I) = IW(I)
   50     CONTINUE
          DO 80 KDUM2 = 1_8,NE
            BMIN = BVAL
            IF (BMAX-BMIN .LE. RLX) GO TO 1000
            CALL DMUMPS_MTRANSQ(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 1000
            K = 1
            DO 70 IDUM = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+int(LEN(J)-1,8),
     &                   IP(J)+int(LENL(J),8),-1_8
                IF (A(II).GE.BVAL) GO TO 60 
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
              LEN(J) = int(II - IP(J) + 1)
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
   81     MOD = 1
        ELSE
          BMAX = BVAL
          IF (BMAX-BMIN .LE. RLX) GO TO 1000
          CALL DMUMPS_MTRANSQ(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 1000
          K = 1
          DO 87 IDUM = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+int(LEN(J),8),IP(J)+int(LENH(J)-1,8)
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = int(II - IP(J))
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL DMUMPS_MTRANSU(CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,
     &              NUM,NUMX,
     &              IW4(1),IW4(N+1),IW4(2*N+1),IW4(2*N+M+1))
   90 CONTINUE 
 1000 IF (M.EQ.N .and. NUMX.EQ.N) GO TO 2000
      CALL DMUMPS_MTRANSX(M,N,IPERM,IW,W)
 2000 RETURN
      END SUBROUTINE DMUMPS_MTRANSS
C
      SUBROUTINE DMUMPS_MTRANSU
     &           (ID,MOD,M,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER  :: ID,MOD,M,N,NUM,NUMX
      INTEGER(8), INTENT(IN) :: LIRN
      INTEGER  :: ARP(N),CV(M),IRN(LIRN),
     &        FC(N),IPERM(M),LENC(N),OUT(N),PR(N)
      INTEGER(8), INTENT(IN) :: IP(N)
      INTEGER I,J,J1,JORD,NFC,K,KK, 
     &        NUM0,NUM1,NUM2,ID0,ID1,LAST
      INTEGER(8) :: IN1, IN2, II
      IF (ID.EQ.1) THEN
        DO 5 I = 1,M
          CV(I) = 0
    5   CONTINUE
        DO 6 J = 1,N
          ARP(J) = 0
    6   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
        IF (MOD.EQ.1) THEN
          DO 8 J = 1,N
            ARP(J) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM
      NFC = 0
      ID0 = (ID-1)*N 
      DO 100 JORD = NUM0+1,N
        ID1 = ID0 + JORD
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + int(ARP(J),8)
          IN2 = IP(J) + int(LENC(J) - 1,8)
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = LENC(J)
   30     OUT(J) = LENC(J) - 1
          DO 60 KK = 1,JORD
            IN1 = int(OUT(J),8)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + int(LENC(J) - 1,8)
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = int(IN2 - II) - 1 
              GO TO 70
   40       CONTINUE
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
   70   CONTINUE
   80   IPERM(I) = J
        ARP(J) = int(II - IP(J)) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + int(LENC(J) - OUT(J) - 2,8)
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
   95   IF (NUM.EQ.NUM1) THEN
          LAST = JORD
          GO TO 101
        ENDIF
  100 CONTINUE
      LAST = N
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_MTRANSU
C
      SUBROUTINE DMUMPS_MTRANSW(M,N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,L32,OUT,PR,Q,L,U,D,RINF) 
      IMPLICIT NONE
      INTEGER :: M,N,NUM
      INTEGER(8), INTENT(IN) :: NE
      INTEGER  :: IRN(NE),IPERM(M),Q(M),L32(max(M,N))
      INTEGER(8)  :: IP(N+1), PR(N), L(M), JPERM(N), OUT(N)
      DOUBLE PRECISION A(NE),U(M),D(M),RINF,RINF3
      INTEGER  :: I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,JSP,
     &            UP,LOW,IK
      INTEGER(8) :: K, KK, KK1, KK2, K0, K1, K2, ISP
      DOUBLE PRECISION     :: CSP,DI,DMIN,DNEW,DQ0,VJ,RLX
      LOGICAL  :: LORD
      DOUBLE PRECISION    :: ZERO, ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      EXTERNAL DMUMPS_MTRANSD, DMUMPS_MTRANSE,
     &         DMUMPS_MTRANSF, DMUMPS_MTRANSX
      RLX = U(1)
      RINF3 = U(2)
      LORD = (JPERM(1).EQ.6)
      NUM = 0
      DO 10 I = 1,N
        JPERM(I) = 0_8
        PR(I) = IP(I)
        D(I) = RINF
   10 CONTINUE
      DO 15 I = 1,M
        U(I) = RINF3
        IPERM(I) = 0
        L(I) = 0_8
   15 CONTINUE
      DO 30 J = 1,N
         IF (int(IP(J+1)-IP(J)) .GT. N/10 .AND. N.GT.50) GO TO 30
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,M
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
        IF (JPERM(J).EQ.0_8) THEN
          JPERM(J) = L(I)
          D(J) = U(I)
          NUM = NUM + 1
        ELSEIF (D(J).GT.U(I)) THEN
          K = JPERM(J)
          II = IRN(K)
          IPERM(II) = 0
          JPERM(J) = L(I)
          D(J) = U(I)
        ELSE
          IPERM(I) = 0
        ENDIF
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 45 I = 1,M
        D(I) = ZERO
   45 CONTINUE
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1_8
        IF (K1.GT.K2) GO TO 95
        VJ = RINF
        DO 50 K = K1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60 
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1_8
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1_8
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1_8
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1_8
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,M
        D(I) = RINF
        Q(I) = 0
   99 CONTINUE
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        DMIN = RINF
        QLEN = 0
        LOW = M + 1
        UP = M + 1
        CSP = RINF
        J = JORD
        PR(J) = -1_8
        DO 115 K = IP(J),IP(J+1)-1_8
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            L(QLEN) = K
          ENDIF
  115   CONTINUE
        Q0 = QLEN
        QLEN = 0
        DO 120 IK = 1,Q0
          K = L(IK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            L32(LOW) = I
            Q(I) = LOW
          ELSE
            QLEN = QLEN + 1
            Q(I) = QLEN
            CALL DMUMPS_MTRANSD(I,M,L32,D,Q,2)
          ENDIF
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = int(J,8)
  120   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = L32(1)
            IF (D(I).LT.RINF) DMIN = D(I)*(ONE+RLX)
            IF (DMIN.GE.CSP) GO TO 160
  152       CALL DMUMPS_MTRANSE(QLEN,M,L32,D,Q,2)
            LOW = LOW - 1
            L32(LOW) = I
            Q(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = L32(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
  153     Q0 = L32(UP-1)
          DQ0 = D(Q0)
          IF (DQ0.GE.CSP) GO TO 160
          IF (DMIN.GE.CSP) GO TO 160
          UP = UP - 1
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          K1 = IP(J+1)-1_8
          IF (LORD) THEN
            IF (CSP.NE.RINF) THEN
              DI = CSP - VJ
              IF (A(K1).GE.DI) THEN
                K0 = JPERM(J)
                IF (K0.GE.K1-6) GO TO 178
  177           CONTINUE
                  K = (K0+K1)/2
                  IF (A(K).GE.DI) THEN 
                    K1 = K
                  ELSE 
                    K0 = K
                  ENDIF
                  IF (K0.GE.K1-6) GO TO 178
                GO TO 177
  178           DO 179 K = K0+1,K1          
                  IF (A(K).LT.DI) GO TO 179
                  K1 = K - 1
                  GO TO 181
  179           CONTINUE
              ENDIF
            ENDIF
  181       IF (K1.EQ.JPERM(J)) K1 = K1 - 1
          ENDIF
          K0 = IP(J)
          DI = CSP - VJ
          DO 155 K = K0,K1
            I = IRN(K)
            IF (Q(I).GE.LOW) GO TO 155
            DNEW = A(K) - U(I)
            IF (DNEW.GE.DI) GO TO 155
            DNEW = DNEW + VJ
            IF (DNEW.GT.D(I)) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = K
              JSP = J
              DI = CSP - VJ
            ELSE
              IF (DNEW.GE.D(I)) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                IF (Q(I).NE.0) THEN
                  CALL DMUMPS_MTRANSF(Q(I),QLEN,M,L32,D,Q,2)
                ENDIF
                LOW = LOW - 1
                L32(LOW) = I
                Q(I) = LOW
              ELSE   
                IF (Q(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  Q(I) = QLEN
                ENDIF
                CALL DMUMPS_MTRANSD(I,M,L32,D,Q,2)
              ENDIF
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = int(J,8)
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.RINF) GO TO 190
        NUM = NUM + 1
        I = IRN(ISP)
        J = JSP
        IPERM(I) = J
        JPERM(J) = ISP
        DO 170 JDUM = 1,NUM
          JJ = int(PR(J))
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
  180   DO 182 JJ = UP,M
          I = L32(JJ)
          U(I) = U(I) + D(I) - CSP
  182   CONTINUE 
  190   DO 191 JJ = UP,M
          I = L32(JJ)
          D(I) = RINF
          Q(I) = 0
  191   CONTINUE
        DO 192 JJ = LOW,UP-1
          I = L32(JJ)
          D(I) = RINF
          Q(I) = 0
  192   CONTINUE 
        DO 193 JJ = 1,QLEN
          I = L32(JJ)
          D(I) = RINF
          Q(I) = 0
  193   CONTINUE
  100 CONTINUE
 1000 CONTINUE
      DO 1200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
 1200 CONTINUE
      DO 1201 I = 1,M
        IF (IPERM(I).EQ.0) U(I) = ZERO
 1201 CONTINUE
      IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000
      CALL DMUMPS_MTRANSX(M,N,IPERM,Q,L32)
 2000 RETURN
      END SUBROUTINE DMUMPS_MTRANSW
      SUBROUTINE DMUMPS_MTRANSZ
     &           (M,N,IRN,LIRN,IP,LENC,IPERM,NUM,PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER :: M,N,NUM
      INTEGER(8), INTENT(IN) :: LIRN
      INTEGER :: ARP(N),CV(M),IRN(LIRN),IPERM(M),LENC(N),OUT(N),PR(N)
      INTEGER(8), INTENT(IN) :: IP(N)
C Local variables 
      INTEGER    :: I,J,J1,JORD,K,KK
      INTEGER(8) :: II, IN1, IN2
      EXTERNAL DMUMPS_MTRANSX
      DO 10 I = 1,M
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      DO 12 J = 1,N
        ARP(J) = LENC(J) - 1
   12 CONTINUE
      NUM = 0
      DO 1000 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
          IN1 = int(ARP(J),8)
          IF (IN1.LT.0_8) GO TO 30
          IN2 = IP(J) + int(LENC(J) - 1,8)
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = -1
   30     CONTINUE
          OUT(J) = LENC(J) - 1
          DO 60 KK = 1,JORD
            IN1 = int(OUT(J),8)
            IF (IN1.LT.0_8) GO TO 50
            IN2 = IP(J) + int(LENC(J) - 1,8)
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.JORD) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = int(IN2 - II - 1_8)
              GO TO 70
   40       CONTINUE
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 1000
   60     CONTINUE
   70   CONTINUE
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = int(IN2 - II - 1_8)
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 1000
          II = IP(J) + int(LENC(J) - OUT(J) - 2,8)
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
 1000 CONTINUE
      IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000
      CALL DMUMPS_MTRANSX(M,N,IPERM,CV,ARP)
 2000 RETURN
      END SUBROUTINE DMUMPS_MTRANSZ
      SUBROUTINE DMUMPS_MTRANSX(M,N,IPERM,RW,CW)
      IMPLICIT NONE
      INTEGER M,N
      INTEGER RW(M),CW(N),IPERM(M)
      INTEGER I,J,K
      DO 10 J = 1,N
        CW(J) = 0
   10 CONTINUE
      K = 0
      DO 20 I = 1,M
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          RW(K) = I
        ELSE
          J = IPERM(I)
          CW(J) = I
        ENDIF
   20 CONTINUE
      K = 0
      DO 30 J = 1,N
        IF (CW(J).NE.0) GO TO 30
        K = K + 1
        I = RW(K)
        IPERM(I) = -J
   30 CONTINUE
      DO 40 J = N+1,M
        K = K + 1
        I = RW(K)
        IPERM(I) = -J
   40 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_MTRANSX
