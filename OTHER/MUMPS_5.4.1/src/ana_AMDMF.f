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
      SUBROUTINE MUMPS_SYMQAMD_NEW
     &                ( JOB, THRESH, NDENSE, 
     &                 N, IWLEN, PE, PFREE, LEN, IW, NV, 
     &                 ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W, 
     &                 PERM, COMPLEM_LIST, SIZE_COMPLEM_LIST, 
     &                 AGG6 ) 
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, SIZE_COMPLEM_LIST
      INTEGER(8), INTENT(IN) :: IWLEN
      INTEGER, INTENT(IN)    :: THRESH
      LOGICAL, INTENT(IN)  ::  AGG6
      INTEGER, INTENT (IN) :: COMPLEM_LIST(max(1,SIZE_COMPLEM_LIST))
      INTEGER, INTENT(INOUT) :: JOB
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
      INTEGER, INTENT(INOUT)    :: PERM(N)
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: NV(N), LAST(N)
      INTEGER, INTENT(INOUT)   :: ELEN(N)
      INTEGER, INTENT(OUT) :: NDENSE(N), DEGREE(N), 
     &                        HEAD(N), NEXT(N), W(N)
      INTEGER THRESM, NDME, PERMeqN
      INTEGER NBD,NBED, NBDM, LASTD, NELME
      LOGICAL IDENSE
      INTEGER :: FDEG, ThresMin, ThresPrev, IBEGSchur, 
     &        ThresMinINIT
      LOGICAL :: AGG6_loc
      INTEGER :: THD_AGG
      LOGICAL :: SchurON
      INTEGER :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, ME, MINDEG, NEL, 
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER KNT1_UPDATED, KNT2_UPDATED
      INTEGER(8) MAXMEM, MEM, NEWMEM
      INTEGER :: MAXINT_N
      INTEGER(8) :: HASH, HMOD 
      INTEGER(8) :: P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, 
     &              PN, PSRC, PLN, PELN
      INTRINSIC max, min, mod
        IF (N.EQ.1) THEN
           ELEN(1) = 1
           LAST(1) = 1
           PE(1) = 0_8
           NV(1) = 1
           RETURN
        ENDIF
        AGG6_loc = AGG6
        THD_AGG = max(128, min(N/2048, 1024))
        IF ( SIZE_COMPLEM_LIST < 0 .OR. SIZE_COMPLEM_LIST > N ) THEN
          WRITE(*,*) "Internal MUMPS_SYMQAMD_NEW", SIZE_COMPLEM_LIST,N
          CALL MUMPS_ABORT()
        ENDIF
        IF (JOB.EQ.2) THEN
          SchurON = .FALSE.        
        ENDIF
        THRESM    = THRESH  
        IF (JOB.NE.2) THEN
          SchurON   = (SIZE_COMPLEM_LIST > 0)
          IF ((JOB.EQ.1) .AND. (.NOT.SchurON) .AND. (N .GT. 0)) THEN
          ENDIF
          IBEGSchur = N-SIZE_COMPLEM_LIST+1
          IF (THRESM.GT.N) THRESM = N
          IF (THRESM.LT.0) THRESM = 0
          IF ( SchurON )  THEN 
           DO I= 1, N
             IF ( PERM(I) .GE. IBEGSchur) THEN 
                 PERM(I) = N + 1
                IF (LEN(I) .EQ.0) THEN
                  PE(I) = 0_8
                ENDIF
             ENDIF
           ENDDO
          ENDIF
        ENDIF
        IF (SchurON) THEN
             THRESM    = N
             ThresMin  = N
             ThresPrev = N
        ELSE
             THRESM    = max(int(31*N/32),THRESM)
             THRESM    = max(THRESM,1)
             ThresMin  = max( 3*THRESM / 4, 1)
             ThresPrev = THRESM
        ENDIF
        ThresMinINIT = ThresMin/4
      IF (THRESM.GT.0) THEN
       IF ((THRESM.GT.N).OR.(THRESM.LT.2)) THEN 
          THRESM = N
       ENDIF
      ENDIF
      IF (JOB.EQ.2) THEN
      ENDIF
      PERMeqN = 0
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
      NEL   = 0
      WFLG   = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA  = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0  
      MEM  = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        NDENSE(I)= 0
        LAST (I) = 0
        HEAD (I) = 0
        NEXT (I) = 0
        NV (I) = 1
        W (I) = 1
   10 CONTINUE
      IF (JOB.EQ.2) THEN
        DO I = 1,SIZE_COMPLEM_LIST
             X       = COMPLEM_LIST(I)
             ELEN(X) = -I       
             NV(X)   = LEN(X)+1 
             DMAX = max(DMAX, LEN(X))
        ENDDO
        NEL = NEL + SIZE_COMPLEM_LIST  
        DO I=1,N
           DEGREE (I) = LEN (I) 
        ENDDO
      ELSE
        DO I=1, N
          ELEN (I) = 0
          DEGREE (I) = LEN (I)
        ENDDO
      ENDIF
      DO 20 I = 1, N
        IF (ELEN(I).LT.0) CYCLE   
        DEG = DEGREE (I)
        IF (PERM(I).EQ.N) THEN
           PERMeqN = I
           PERM(I) = N-1
        ENDIF
        FDEG = PERM(I)
        IF ( (DEG .GT. 0).OR.(PERM(I).EQ.N+1) ) THEN
          IF ( (THRESM.GT.0) .AND.
     &         (FDEG .GT.THRESM) ) THEN
            NBD = NBD+1
            IF (FDEG.NE.N+1) THEN
             DEGREE(I) = DEGREE(I)+N+2
             DEG = N
             INEXT = HEAD (DEG)
             IF (INEXT .NE. 0) LAST (INEXT) = I
             NEXT (I) = INEXT
             HEAD (DEG) = I 
             LAST(I)  = 0
             IF (LASTD.EQ.0) LASTD=I
            ELSE
             NBED = NBED+1
             DEGREE(I) = N+1
             DEG = N
             IF (LASTD.EQ.0) THEN
               LASTD     = I 
               HEAD(DEG) = I
               NEXT(I)   = 0 
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
            ENDIF
          ELSE
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            HEAD (FDEG) = I
          ENDIF
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0_8
          W (I) = 0
        ENDIF
   20 CONTINUE
          IF ((NBD.EQ.0).AND.(THRESM.GT.0)) THRESM = N
   30 IF (NEL .LT. N) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF ( (DEG.NE.N) .AND.
     &    (DEG.GT.THRESM+1) .AND. (NBD.GT.0) ) THEN
           MINDEG = N
           GOTO 30
        ENDIF
        IF (DEGREE(ME).LE.N)  THEN
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELSE
          MINDEG = 1
          NBDM = max(NBDM,NBD)
          IF (DEGREE(ME).GT.N+1) THEN
            IF (WFLG .GT. MAXINT_N) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN 
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + int(LEN(ME) -1,8)
            PLN       = P1
            PELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0_8) THEN
                X = E
  53            X = int(-PE(X))
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0 ) GOTO 53
                E = X
              ENDIF
              IF (ELEN(E).LT.0) THEN
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(PLN) = IW(PELN)
               IW(PELN) = E
               PLN  = PLN  + 1_8
               PELN = PELN + 1_8
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+LEN(E)-1
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(PLN)=E
               PLN = PLN+1_8
              ENDIF
  55        CONTINUE
            WFLG     = WFLG + 1
            LEN(ME)  = int(PLN-P1)
            ELEN(ME) = int(PELN-P1)
            NDME = NDENSE(ME)+NV(ME)
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
            DEGREE(ME) = NDENSE(ME)
            DEG = PERM(ME)
            MINDEG = min(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(N+1) ) GOTO 51
            ENDIF
            HEAD (N) = ME
            IF (THRESM.LT.N) THEN
             ThresMin  = max(THRESM+ThresMin,ThresPrev+ThresMin/2+1)
             ThresMin  = min(ThresMin, N)
             ThresPrev = ThresPrev+(N-ThresPrev)/2+ThresMinINIT
             THRESM    = max(
     &         THRESM + int(sqrt(dble(ThresMin)))+ ThresMinINIT ,
     &         ThresPrev)
             THRESM    = min(THRESM,N) 
             ThresMin  = min(THRESM, ThresMin)
             ThresPrev = THRESM
            ENDIF
            NBD    = NBED
            GOTO 30
          ENDIF
          IF (DEGREE(ME).EQ.N+1) THEN
             IF (NBD.NE.NBED) THEN
          write(6,*) ' ERROR in MUMPS_SYMQAMD_NEW ',
     &                ' quasi dense rows remains'
          CALL MUMPS_ABORT()
           ENDIF
           IF (JOB.EQ.1) THEN
            DO I = 1,SIZE_COMPLEM_LIST
             X       = COMPLEM_LIST(I)
             ELEN(X) = -(N-SIZE_COMPLEM_LIST+I)   
             NV(X)   = 1                
             PE(X)   = 0_8              
            ENDDO
            GOTO 265
           ENDIF
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0_8) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = int(-COMPLEM_LIST(1),8)
            ELSEIF (DEGREE(X).EQ.N+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = int(-ME,8)
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
   59      CONTINUE
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0_8
           IF (NEL.NE.N) THEN
            write(6,*) 'Internal ERROR 2 detected in QAMD'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            CALL MUMPS_ABORT()
           ENDIF
           IF (ME.NE. COMPLEM_LIST(1)) THEN
             DO I=1, SIZE_COMPLEM_LIST
               PE(COMPLEM_LIST(I)) = int(-COMPLEM_LIST(1),8)
             ENDDO
             PE(COMPLEM_LIST(1)) = 0_8
             NV( COMPLEM_LIST(1))= NV(ME)
             NV(ME)               = 0
             ELEN( COMPLEM_LIST(1)) = ELEN(ME)
             ELEN(ME)             = 0
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NDENSE(ME) = 0
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + int(LEN (ME) - 1,8)
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).LE.N) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (PERM(I)) = INEXT
              ENDIF
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
                  KNT1_UPDATED = 0
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0_8
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
                  KNT2_UPDATED = 0
                  IF (LEN (E) .EQ. 0) PE (E) = 0_8
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = int(IW (PN),8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1_8
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + int(LENJ - 1,8)
                      PSRC = PSRC + int(LENJ - 1,8)
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                IF (DEGREE(I).LE.N) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (PERM(I)) = INEXT
                ENDIF
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1_8)
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.N) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + int(ELN - 1,8)
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI - NDENSE(E)
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        AGG6_loc = (AGG6 .OR. (DEGREE(ME) .LT. THD_AGG))
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.N) GOTO 180
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (.NOT. AGG6_loc .AND. DEXT .EQ. 0) THEN
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (AGG6_loc .AND. (DEXT .EQ. 0) .AND.
     &            ((NDENSE(ME).EQ.NBD).OR.(NDENSE(E).EQ.0))) THEN
                PE (E) = int(-ME,8)
                W (E)  = 0
             ELSE IF (AGG6_loc .AND. DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + int(E,kind=8)
            ENDIF
  160     CONTINUE
          ELEN (I) = int(PN - P1 + 1_8)
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              IF (DEGREE(J).LE.N) DEG=DEG+NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
          IF (((ELEN(I).EQ.1).AND.(P3.EQ.PN))
     &     .OR.
     &         (AGG6_loc.AND.(DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD))
     &       )
     &    THEN
            PE (I) = int(-ME, 8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE(I) = min (DEG+NBD-NDENSE(ME), 
     &                       DEGREE(I))
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = int(PN - P1 + 1)
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.N) ) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
             X = I 
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + int(LN - 1,8)
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + int(LN - 1,8)
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                IF (PERM(J).GT.PERM(X)) THEN
                  PE (J) = int(-X,8)
                  NV (X) = NV (X) + NV (J)
                  NV (J) = 0
                  ELEN (J) = 0
                ELSE
                  PE (X) = int(-J,8)
                  NV (J) = NV (X) + NV (J)
                  NV (X) = 0
                  ELEN (X) = 0
                  X = J
                ENDIF
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        IF ( (THRESM .GT. 0).AND.(THRESM.LT.N) ) THEN 
          THRESM = max(ThresMin, THRESM-NVPIV)
        ENDIF
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).LE.N) THEN
            DEG = min (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
            IF (THRESM.GT.0) THEN
             IF (PERM(I) .GT. THRESM) THEN
               IDENSE = .TRUE.
               DEGREE(I) = DEGREE(I)+N+2
             ENDIF
             IF (IDENSE) THEN
               P1 = PE(I)
               P2 = P1 + int(ELEN(I) - 1, 8)
               IF (P2.GE.P1) THEN
               DO 264 PJ=P1,P2
                 E= IW(PJ)
                 NDENSE (E) = NDENSE(E) + NVI
 264           CONTINUE
               ENDIF
               NBD = NBD+NVI
               FDEG = N
               DEG = N
               INEXT = HEAD(DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               HEAD (DEG) = I
               LAST(I)    = 0
               IF (LASTD.EQ.0) LASTD=I
             ENDIF
            ENDIF
            IF (.NOT.IDENSE) THEN
            FDEG = PERM(I)
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (FDEG) = I
            ENDIF
            MINDEG = min (MINDEG, FDEG)
            ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = int(P - PME1)
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0_8
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + int(LEN (ME),8)
        ENDIF
      GO TO 30
      ENDIF
  265 CONTINUE
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = int(-PE (I))
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = int(-PE (J))
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = int(-PE (J))
              PE (J)= int(-E,8)
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = abs (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
      IF (.NOT.SchurON) THEN
        IF (PERMeqN.GT.0) PERM(PERMeqN) = N
      ENDIF
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_SYMQAMD_NEW
