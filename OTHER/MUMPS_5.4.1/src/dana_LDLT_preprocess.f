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
      SUBROUTINE DMUMPS_SET_CONSTRAINTS(
     &     N,PIV,FRERE,FILS,NFSIZ,IKEEP,
     &     NCST,KEEP,KEEP8, ROWSCA
     &     )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(OUT) :: NCST
      INTEGER    :: PIV(N),FRERE(N),FILS(N),NFSIZ(N),IKEEP(N)
      INTEGER    :: KEEP(500)
      INTEGER(8) :: KEEP8(150)
      DOUBLE PRECISION       :: ROWSCA(N)
      INTEGER I,P11,P1,P2,K1,K2,NLOCKED
      LOGICAL V1,V2
      NCST = 0
      NLOCKED = 0
      P11 = KEEP(93)
      DO I=KEEP(93)-1,1,-2
         P1 = PIV(I)
         P2 = PIV(I+1)
         K1 = IKEEP(P1)
         IF (K1 .NE. 0) THEN
           V1 = (K1+2*exponent(ROWSCA(P1)) .GE. -3)
         ELSE
           V1 = .FALSE.
         ENDIF
         K2 = IKEEP(P2)
         IF (K2 .NE. 0) THEN
           V2 =  (K2+exponent(ROWSCA(P2)**2) .GE. -3)
         ELSE
           V2 = .FALSE.
         ENDIF
         IF(V1 .AND. V2) THEN
            PIV(P11) = P1
            P11 = P11 - 1
            PIV(P11) = P2
            P11 = P11 - 1
         ELSE IF(V1) THEN
            NCST = NCST+1
            FRERE(NCST) = P1
            NCST = NCST+1
            FRERE(NCST) = P2
         ELSE IF(V2) THEN
            NCST = NCST+1
            FRERE(NCST) = P2                
            NCST = NCST+1
            FRERE(NCST) = P1
         ELSE
            NLOCKED = NLOCKED + 1
            FILS(NLOCKED) = P1
            NLOCKED = NLOCKED + 1
            FILS(NLOCKED) = P2                   
         ENDIF
      ENDDO
      DO I=1,NLOCKED
         PIV(I) = FILS(I)
      ENDDO
      KEEP(94) = KEEP(94) + KEEP(93) - NLOCKED
      KEEP(93) = NLOCKED
      DO I=1,NCST
         NLOCKED = NLOCKED + 1
         PIV(NLOCKED) = FRERE(I)
      ENDDO
      DO I=1,KEEP(93)/2
         NFSIZ(I) = 0
      ENDDO
      DO I=(KEEP(93)/2)+1,(KEEP(93)/2)+NCST,2
         NFSIZ(I) = I+1
         NFSIZ(I+1) = -1
      ENDDO
      DO I=(KEEP(93)/2)+NCST+1,(KEEP(93)/2)+KEEP(94)
         NFSIZ(I) = 0
      ENDDO
      END SUBROUTINE DMUMPS_SET_CONSTRAINTS
      SUBROUTINE DMUMPS_EXPAND_PERMUTATION(N,NCMP,N11,N22,PIV,
     &     INVPERM,PERM)
      IMPLICIT NONE
      INTEGER N11,N22,N,NCMP
      INTEGER, intent(in) :: PIV(N),PERM(N)
      INTEGER, intent (out):: INVPERM(N)
      INTEGER CMP_POS,EXP_POS,I,J,N2,K
      N2 = N22/2
      EXP_POS = 1
      DO CMP_POS=1,NCMP
         J = PERM(CMP_POS)
         IF(J .LE. N2) THEN
            K = 2*J-1
            I = PIV(K)
            INVPERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
            K = K+1
            I = PIV(K)
            INVPERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
         ELSE
            K = N2 + J
            I = PIV(K)
            INVPERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
         ENDIF
      ENDDO
      DO K=N22+N11+1,N
         I = PIV(K)
         INVPERM(I) = EXP_POS
         EXP_POS = EXP_POS+1
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_EXPAND_PERMUTATION
      SUBROUTINE DMUMPS_LDLT_COMPRESS(
     &     N,NZ, IRN, ICN, PIV,
     &     NCMP, IW, LW, IPE, LEN, IQ, 
     &     FLAG, ICMP, IWFR,
     &     IERROR, KEEP,KEEP8, ICNTL,INPLACE64_GRAPH_COPY)
      IMPLICIT NONE
      INTEGER, intent(in)     :: N
      INTEGER(8), intent(in)  :: NZ, LW
      INTEGER, intent(in)     :: IRN(NZ), ICN(NZ), PIV(N)
      INTEGER, intent(in)     :: ICNTL(60)
      INTEGER, intent(in)     :: KEEP(500)
      INTEGER(8), intent(in)  :: KEEP8(150)
      INTEGER, intent(out)    :: NCMP, IERROR
      INTEGER(8), intent(out) :: IWFR, IPE(N+1)
      INTEGER, intent(out)    :: IW(LW) 
      INTEGER, intent(out)    :: LEN(N)
      INTEGER(8), intent(out) :: IQ(N)
      INTEGER, intent(out)    :: FLAG(N), ICMP(N)
      LOGICAL, intent(inout) :: INPLACE64_GRAPH_COPY
      INTEGER    :: MP, N11, N22
      INTEGER    :: I, J, N1, K
      INTEGER(8) :: NDUP, L, K8, K1, K2, LAST
      INTRINSIC nint
      MP = ICNTL(2)
      IERROR = 0
      N22 = KEEP(93)
      N11 = KEEP(94)
      NCMP = N22/2 + N11
      DO I=1,NCMP
         IPE(I) = 0
      ENDDO
      K = 1
      DO I=1,N22/2
         J = PIV(K)
         ICMP(J) = I
         K = K + 1
         J = PIV(K)
         ICMP(J) = I
         K = K + 1
      ENDDO
      K = N22/2 + 1
      DO I=N22+1,N22+N11
         J = PIV(I)
         ICMP(J) = K
         K = K + 1
      ENDDO
      DO I=N11+N22+1,N
         J = PIV(I)
         ICMP(J) = 0
      ENDDO
      DO K8=1,NZ
         I = IRN(K8)
         J = ICN(K8)
         IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     &        .OR.(J.LT.1)) THEN
            IERROR = IERROR + 1
         ELSE
         I = ICMP(I)
         J = ICMP(J)
         IF ((I.NE.0).AND.(J.NE.0).AND.(I.NE.J)) THEN
               IPE(I) = IPE(I) + 1_8
               IPE(J) = IPE(J) + 1_8
            ENDIF
         ENDIF
      ENDDO
      IQ(1) = 1_8
      N1 = NCMP - 1
      IF (N1.GT.0) THEN
         DO I=1,N1
            IQ(I+1) = IPE(I) + IQ(I)
         ENDDO 
      ENDIF
      LAST = max(IPE(NCMP)+IQ(NCMP)-1_8,IQ(NCMP))
      DO I = 1,NCMP
         FLAG(I) = 0
         IPE(I)  = IQ(I)
      ENDDO
      IW(1:LAST) = 0
      IWFR = LAST + 1_8
      DO K8=1,NZ
         I = IRN(K8)
         J = ICN(K8)
         IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     &        .OR.(J.LT.1)) CYCLE
         I = ICMP(I)
         J = ICMP(J)
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
            IF ((I.GE.1).AND.(J.LE.N)) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1_8
            ENDIF
          ELSE
            IF ((J.GE.1).AND.(I.LE.N)) THEN
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1_8
            ENDIF
          ENDIF
         ENDIF
      ENDDO
      NDUP = 0_8
      DO I=1,NCMP
         K1 = IPE(I) 
         K2 = IQ(I) -1_8
         IF (K1.GT.K2) THEN
            LEN(I) = 0
         ELSE
            DO K8=K1,K2
               J     = -IW(K8)
               IF (J.LE.0) GO TO 250
               L     = IQ(J) 
               IQ(J) = L + 1_8
               IF (FLAG(J).EQ.I) THEN
                  NDUP = NDUP + 1_8
                  IW(L) = 0
                  IW(K8) = 0
               ELSE
                  IW(L)   = I
                  IW(K8)   = J
                  FLAG(J) = I
               ENDIF
            ENDDO
 250        LEN(I) = int(IQ(I) - IPE(I))
         ENDIF
      ENDDO
      IF (NDUP.NE.0_8) THEN
         IWFR = 1_8
         DO I=1,NCMP
            K1 = IPE(I) 
            IF (LEN(I).EQ.0) THEN
               IPE(I) = IWFR
               CYCLE
            ENDIF
            K2 = K1 + LEN(I) - 1
            L = IWFR
            IPE(I) = IWFR
            DO K8=K1,K2
               IF (IW(K8).NE.0) THEN
                  IW(IWFR) = IW(K8)
                  IWFR     = IWFR + 1_8
               ENDIF
            ENDDO
            LEN(I) = int(IWFR - L)
         ENDDO
      ENDIF
      IPE(NCMP+1) = IPE(NCMP) + int(LEN(NCMP),8)
      IWFR = IPE(NCMP+1)
      INPLACE64_GRAPH_COPY = (LW.GE.2*(IWFR-1_8))
      RETURN
      END SUBROUTINE DMUMPS_LDLT_COMPRESS
      SUBROUTINE DMUMPS_SYM_MWM(
     &     N, NE, IP, IRN, SCALING,LSC,CPERM, DIAG,
     &     ICNTL, WEIGHT,MARKED,FLAG,
     &     PIV_OUT, INFO)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N
      INTEGER(8), INTENT(IN) :: NE
      INTEGER :: ICNTL(10), INFO(10),LSC
      INTEGER :: CPERM(N),PIV_OUT(N), IRN(NE), DIAG(N)
      INTEGER(8), INTENT(IN) :: IP(N+1)
      DOUBLE PRECISION    :: SCALING(LSC),WEIGHT(N+2)
      INTEGER :: MARKED(N),FLAG(N)
      INTEGER :: NUM1,NUM2,NUMTOT,PATH_LENGTH,NLAST
      INTEGER :: I,BEST_BEG, CUR_EL,CUR_EL_PATH,CUR_EL_PATH_NEXT
      INTEGER :: L1,L2,TUP,T22
      INTEGER(8) :: PTR_SET1,PTR_SET2
      DOUBLE PRECISION    :: BEST_SCORE,CUR_VAL,TMP,VAL
      DOUBLE PRECISION INITSCORE, DMUMPS_UPDATESCORE, 
     &     DMUMPS_UPDATE_INVERSE, DMUMPS_METRIC2x2
      LOGICAL VRAI,FAUX,MAX_CARD_DIAG,USE_SCALING
      INTEGER SUM
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (SUM = 1, VRAI = .TRUE., FAUX = .FALSE.)
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      MAX_CARD_DIAG = .TRUE.
      NUM1 = 0
      NUM2 = 0
      NUMTOT = 0
      NLAST = N
      INFO = 0      
      MARKED = 1
      FLAG = 0
      VAL = ONE
      IF(LSC .GT. 1) THEN
         USE_SCALING = .TRUE.
      ELSE
         USE_SCALING = .FALSE.
      ENDIF
      TUP = ICNTL(2)
      IF(TUP .EQ. SUM) THEN
        INITSCORE = ZERO
      ELSE
        INITSCORE = ONE
      ENDIF
      IF(ICNTL(2) .GT. 2 .OR. ICNTL(2) .LE. 0) THEN
         WRITE(*,*)
     &        'ERROR: WRONG VALUE FOR ICNTL(2) = ',ICNTL(2)
         INFO(1) = -1
         RETURN
      ENDIF
      T22 = ICNTL(1)
      IF(ICNTL(1) .LT. 0 .OR. ICNTL(1) .GT. 2) THEN
         WRITE(*,*)
     &        'ERROR: WRONG VALUE FOR ICNTL(1) = ',ICNTL(1)
         INFO(1) = -1
         RETURN
      ENDIF
      DO CUR_EL=1,N
         IF(MARKED(CUR_EL) .LE. 0) THEN
            CYCLE
         ENDIF
         IF(CPERM(CUR_EL) .LT. 0) THEN
            MARKED(CUR_EL) = -1
            CYCLE
         ENDIF
         PATH_LENGTH = 2
         CUR_EL_PATH = CPERM(CUR_EL)
         IF(CUR_EL_PATH .EQ. CUR_EL) THEN
            MARKED(CUR_EL) = -1
            CYCLE
         ENDIF
         MARKED(CUR_EL) = 0
         WEIGHT(1) = INITSCORE
         WEIGHT(2) = INITSCORE
         L1 = int(IP(CUR_EL+1)-IP(CUR_EL))
         L2 = int(IP(CUR_EL_PATH+1)-IP(CUR_EL_PATH))
         PTR_SET1 = IP(CUR_EL)
         PTR_SET2 = IP(CUR_EL_PATH)
         IF(USE_SCALING) THEN
            VAL = -SCALING(CUR_EL_PATH) - SCALING(CUR_EL+N)
         ENDIF
         CUR_VAL = DMUMPS_METRIC2x2(
     &        CUR_EL,CUR_EL_PATH,
     &        IRN(PTR_SET1),IRN(PTR_SET2),
     &        L1,L2,
     &        VAL,DIAG,N,FLAG,FAUX,T22)
         WEIGHT(PATH_LENGTH+1) = 
     &        DMUMPS_UPDATESCORE(WEIGHT(1),CUR_VAL,TUP)
         DO
            IF(CUR_EL_PATH .EQ. CUR_EL) EXIT
            PATH_LENGTH = PATH_LENGTH+1
            MARKED(CUR_EL_PATH) = 0
            CUR_EL_PATH_NEXT = CPERM(CUR_EL_PATH)
            L1 = int(IP(CUR_EL_PATH+1)-IP(CUR_EL_PATH))
            L2 = int(IP(CUR_EL_PATH_NEXT+1)-IP(CUR_EL_PATH_NEXT))
            PTR_SET1 = IP(CUR_EL_PATH)
            PTR_SET2 = IP(CUR_EL_PATH_NEXT)
            IF(USE_SCALING) THEN
               VAL = -SCALING(CUR_EL_PATH_NEXT) 
     &              - SCALING(CUR_EL_PATH+N)
            ENDIF
            CUR_VAL = DMUMPS_METRIC2x2(
     &           CUR_EL_PATH,CUR_EL_PATH_NEXT,
     &           IRN(PTR_SET1),IRN(PTR_SET2),
     &           L1,L2,
     &           VAL,DIAG,N,FLAG,VRAI,T22)
            WEIGHT(PATH_LENGTH+1) = 
     &           DMUMPS_UPDATESCORE(WEIGHT(PATH_LENGTH-1),CUR_VAL,TUP)
            CUR_EL_PATH = CUR_EL_PATH_NEXT
         ENDDO
         IF(mod(PATH_LENGTH,2) .EQ. 1) THEN
            IF(WEIGHT(PATH_LENGTH+1) .GE. WEIGHT(PATH_LENGTH)) THEN
               CUR_EL_PATH = CPERM(CUR_EL)
            ELSE
               CUR_EL_PATH = CUR_EL
            ENDIF
            DO I=1,(PATH_LENGTH-1)/2
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
            NUMTOT = NUMTOT + PATH_LENGTH - 1
         ELSE
            IF(MAX_CARD_DIAG) THEN
               CUR_EL_PATH = CPERM(CUR_EL)
               IF(DIAG(CUR_EL) .NE. 0) THEN 
                  BEST_BEG = CUR_EL_PATH
                  GOTO 1000
               ENDIF
               DO I=1,(PATH_LENGTH/2)
                  CUR_EL_PATH_NEXT = CPERM(CUR_EL_PATH)
                  IF(DIAG(CUR_EL_PATH) .NE. 0) THEN 
                     BEST_BEG = CUR_EL_PATH_NEXT
                     GOTO 1000
                  ENDIF
               ENDDO
            ENDIF
            BEST_BEG = CUR_EL
            BEST_SCORE = WEIGHT(PATH_LENGTH-1)
            CUR_EL_PATH = CPERM(CUR_EL)
            DO I=1,(PATH_LENGTH/2)-1
               TMP = DMUMPS_UPDATESCORE(WEIGHT(PATH_LENGTH),
     &              WEIGHT(2*I-1),TUP)
               TMP = DMUMPS_UPDATE_INVERSE(TMP,WEIGHT(2*I),TUP)
               IF(TMP .GT. BEST_SCORE) THEN
                  BEST_SCORE = TMP
                  BEST_BEG = CUR_EL_PATH
               ENDIF
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               TMP = DMUMPS_UPDATESCORE(WEIGHT(PATH_LENGTH+1),
     &              WEIGHT(2*I),TUP)
               TMP = DMUMPS_UPDATE_INVERSE(TMP,WEIGHT(2*I+1),TUP)
               IF(TMP .GT. BEST_SCORE) THEN
                  BEST_SCORE = TMP
                  BEST_BEG = CUR_EL_PATH
               ENDIF
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
 1000       CUR_EL_PATH = BEST_BEG
            DO I=1,(PATH_LENGTH/2)-1
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
            NUMTOT = NUMTOT + PATH_LENGTH - 2
            MARKED(CUR_EL_PATH) = -1
         ENDIF
      ENDDO
      DO I=1,N
         IF(MARKED(I) .LT. 0) THEN
            IF(DIAG(I) .EQ. 0) THEN
               PIV_OUT(NLAST) = I
               NLAST = NLAST - 1
            ELSE
               NUM1 = NUM1 + 1
               PIV_OUT(NUM2+NUM1) = I
               NUMTOT = NUMTOT + 1
            ENDIF
         ENDIF
      ENDDO
      INFO(2) = NUMTOT
      INFO(3) = NUM1
      INFO(4) = NUM2
      RETURN
      END SUBROUTINE DMUMPS_SYM_MWM
      FUNCTION DMUMPS_UPDATESCORE(A,B,T)
      IMPLICIT NONE
      DOUBLE PRECISION DMUMPS_UPDATESCORE
      DOUBLE PRECISION A,B
      INTEGER T
      INTEGER SUM
      PARAMETER(SUM = 1)
      IF(T .EQ. SUM) THEN
         DMUMPS_UPDATESCORE = A+B
      ELSE
         DMUMPS_UPDATESCORE = A*B
      ENDIF
      END FUNCTION DMUMPS_UPDATESCORE
      FUNCTION DMUMPS_UPDATE_INVERSE(A,B,T)
      IMPLICIT NONE
      DOUBLE PRECISION DMUMPS_UPDATE_INVERSE
      DOUBLE PRECISION A,B
      INTEGER T
      INTEGER SUM
      PARAMETER(SUM = 1)
      IF(T .EQ. SUM) THEN
         DMUMPS_UPDATE_INVERSE = A-B
      ELSE
         DMUMPS_UPDATE_INVERSE = A/B
      ENDIF
      END FUNCTION DMUMPS_UPDATE_INVERSE
      FUNCTION DMUMPS_METRIC2x2(CUR_EL,CUR_EL_PATH,
     &     SET1,SET2,L1,L2,VAL,DIAG,N,FLAG,FLAGON,T)
      IMPLICIT NONE
      DOUBLE PRECISION DMUMPS_METRIC2x2
      INTEGER CUR_EL,CUR_EL_PATH,L1,L2,N
      INTEGER SET1(L1),SET2(L2),DIAG(N),FLAG(N)
      DOUBLE PRECISION VAL
      LOGICAL FLAGON
      INTEGER T
      INTEGER I,INTER,MERGE
      INTEGER STRUCT,MA47
      PARAMETER(STRUCT=0,MA47=1)
      IF(T .EQ. STRUCT) THEN
         IF(.NOT. FLAGON) THEN
            DO I=1,L1
               FLAG(SET1(I)) = CUR_EL
            ENDDO            
         ENDIF
         INTER = 0
         DO I=1,L2
            IF(FLAG(SET2(I)) .EQ. CUR_EL) THEN
               INTER = INTER + 1
               FLAG(SET2(I)) = CUR_EL_PATH
            ENDIF
         ENDDO
         MERGE = L1 + L2 - INTER
         DMUMPS_METRIC2x2 = dble(INTER) / dble(MERGE)
      ELSE IF (T .EQ. MA47) THEN
         MERGE = 3
         IF(DIAG(CUR_EL) .NE. 0) MERGE = 2
         IF(DIAG(CUR_EL_PATH) .NE. 0) MERGE = MERGE - 2
         IF(MERGE .EQ. 0) THEN
            DMUMPS_METRIC2x2 = dble(L1+L2-2)
            DMUMPS_METRIC2x2 = -(DMUMPS_METRIC2x2**2)/2.0D0
         ELSE IF(MERGE .EQ. 1) THEN
            DMUMPS_METRIC2x2 = - dble(L1+L2-4) * dble(L1-2)
         ELSE IF(MERGE .EQ. 2) THEN
            DMUMPS_METRIC2x2 = - dble(L1+L2-4) * dble(L2-2)
         ELSE
            DMUMPS_METRIC2x2 = - dble(L1-2) * dble(L2-2)
         ENDIF
      ELSE
         DMUMPS_METRIC2x2 = VAL
      ENDIF
      RETURN
      END FUNCTION 
      SUBROUTINE DMUMPS_EXPAND_PERM_SCHUR(NA, NCMP,
     &      INVPERM,PERM, 
     &      LISTVAR_SCHUR, SIZE_SCHUR, AOTOA)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: SIZE_SCHUR, LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN):: NA, NCMP
      INTEGER, INTENT(IN):: AOTOA(NCMP), PERM(NCMP)
      INTEGER, INTENT(OUT):: INVPERM(NA) 
      INTEGER CMP_POS, IO, I, K, IPOS
      DO CMP_POS=1, NCMP
        IO              = PERM(CMP_POS)
        INVPERM(AOTOA(IO)) = CMP_POS
      ENDDO
      IPOS = NCMP
      DO K =1,  SIZE_SCHUR
        I       = LISTVAR_SCHUR(K)
        IPOS    = IPOS+1
        INVPERM(I) = IPOS
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_EXPAND_PERM_SCHUR
      SUBROUTINE DMUMPS_GNEW_SCHUR
     & (NA, N, NZ, IRN, ICN, IW, LW, IPE, LEN,
     & IQ, FLAG, IWFR,
     & NRORM, NIORM, IFLAG,IERROR, ICNTL, 
     & symmetry, SYM, NBQD, AvgDens,
     & KEEP264, KEEP265, 
     & LISTVAR_SCHUR, SIZE_SCHUR, ATOAO, AOTOA,
     & INPLACE64_GRAPH_COPY
     & )
      IMPLICIT NONE
      INTEGER, intent(in)    :: NA
      INTEGER, intent(in)    :: N, SYM
      INTEGER(8), intent(in) :: NZ, LW
      INTEGER, intent(in)    :: ICNTL(60)
      INTEGER, intent(in)    :: IRN(NZ), ICN(NZ) 
      INTEGER, INTENT(IN)  :: SIZE_SCHUR, LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, intent(out)   :: IERROR, symmetry
      INTEGER, intent(out)   :: NBQD, AvgDens
      INTEGER, intent(out)   :: LEN(N), IW(LW)
      INTEGER(8), intent(out):: IWFR
      INTEGER(8), intent(out):: NRORM, NIORM
      INTEGER(8), intent(out):: IPE(N+1)
      INTEGER, INTENT(OUT) :: AOTOA(N) 
      INTEGER, INTENT(OUT) :: ATOAO(NA)
      INTEGER, intent(inout) :: IFLAG, KEEP264 
      INTEGER, intent(in)    :: KEEP265
      INTEGER(8), intent(out):: IQ(N)
      INTEGER, intent(out)   :: FLAG(N)
      LOGICAL, intent(inout) :: INPLACE64_GRAPH_COPY
      INTEGER    :: MP, MPG, I, J, N1
      INTEGER    :: NBERR, THRESH, IAO
      INTEGER(8) :: K8, K1, K2, LAST, NDUP
      INTEGER(8) :: NZOFFA, NDIAGA, L, N8
      DOUBLE PRECISION       :: RSYM
      INTRINSIC nint
      MP = ICNTL(2)
      MPG= ICNTL(3)
      ATOAO(1:NA) = 0
      DO I = 1, SIZE_SCHUR
        ATOAO(LISTVAR_SCHUR(I)) = -1
      ENDDO
      IAO = 0  
      DO I= 1, NA
        IF (ATOAO(I).LT.0) CYCLE
        IAO = IAO +1   
        ATOAO(I)   = IAO
        AOTOA(IAO) = I
      ENDDO
      NZOFFA = 0_8
      NDIAGA = 0
      IERROR = 0
      N8     = int(N,8)
      DO I=1,N+1
        IPE(I) = 0_8
      ENDDO
      IF (KEEP264.EQ.0) THEN
       IF ((SYM.EQ.0).AND.(KEEP265.EQ.-1)) THEN
        DO K8=1_8,NZ
         I = IRN(K8)
         J = ICN(K8)
         IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     &                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
         ELSE
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1_8
           NZOFFA  = NZOFFA + 1_8
          ELSE
           NDIAGA = NDIAGA + 1_8
          ENDIF
         ENDIF
        ENDDO
       ELSE
        DO K8=1_8,NZ
         I = IRN(K8)
         J = ICN(K8)
         IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     &                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
         ELSE
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1_8
           IPE(J) = IPE(J) + 1_8
           NZOFFA  = NZOFFA + 1_8
          ELSE
           NDIAGA = NDIAGA + 1_8
          ENDIF
         ENDIF
        ENDDO
       ENDIF
       IF (IERROR.GE.1) THEN
        KEEP264 = 0
       ELSE
        KEEP264 = 1
       ENDIF
      ELSE
       IF ((SYM.EQ.0).AND.(KEEP265.EQ.-1)) THEN
        DO K8=1_8,NZ
         I = IRN(K8)
         J = ICN(K8)
         I = ATOAO(I)
         J = ATOAO(J)
         IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
         IF (I.EQ.J) THEN
           NDIAGA = NDIAGA + 1_8
         ELSE 
           IPE(I) = IPE(I) + 1_8
           NZOFFA = NZOFFA + 1_8
         ENDIF
        ENDDO
       ELSE
        DO K8=1_8,NZ
         I = IRN(K8)
         J = ICN(K8)
         I = ATOAO(I)
         J = ATOAO(J)
         IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
         IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1_8
           IPE(J) = IPE(J) + 1_8
           NZOFFA  = NZOFFA + 1_8
         ELSE
           NDIAGA = NDIAGA + 1_8
         ENDIF
        ENDDO
       ENDIF
      ENDIF
      NIORM  = NZOFFA + 3_8*N8
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IF (mod(IFLAG,2) .EQ. 0) IFLAG = IFLAG+1
         IF ((MP.GT.0).AND.(ICNTL(4).GE.2))  THEN 
          WRITE (MP,99999) 
          DO 70 K8=1_8,NZ
           I = IRN(K8)
           J = ICN(K8)
           IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     &                            .OR.(J.LT.1)) THEN
            NBERR = NBERR + 1
            IF (NBERR.LE.10)  THEN
               IF (mod(K8,10_8).GT.3_8 .OR. mod(K8,10_8).EQ.0_8 .OR.
     &             (10_8.LE.K8 .AND. K8.LE.20_8)) THEN
                 WRITE (MP,'(I16,A,I10,A,I10,A)')
     &             K8,'th entry (in row',I,' and column',J,') ignored'
               ELSE
                 IF (mod(K8,10_8).EQ.1_8) 
     &             WRITE(MP,'(I16,A,I10,A,I10,A)')
     &             K8,'st entry (in row',I,' and column',J,') ignored'
                 IF (mod(K8,10_8).EQ.2_8) 
     &             WRITE(MP,'(I16,A,I10,A,I10,A)')
     &             K8,'nd entry (in row',I,' and column',J,') ignored'
                 IF (mod(K8,10_8).EQ.3_8) 
     &             WRITE(MP,'(I16,A,I10,A,I10,A)')
     &             K8,'rd entry (in row',I,' and column',J,') ignored'
               ENDIF
            ELSE
               GO TO 100
            ENDIF
           ENDIF
   70     CONTINUE
         ENDIF
      ENDIF
  100 NRORM = NIORM - 2_8*N8
      IQ(1) = 1_8
      N1 = N - 1
      IF (N1.GT.0) THEN
        DO I=1,N1
            IQ(I+1) = IPE(I) + IQ(I) 
        ENDDO
      ENDIF
      LAST = max(IPE(N)+IQ(N)-1,IQ(N))
      FLAG(1:N) = 0
      IPE(1:N)  = IQ(1:N)
      IW(1:LAST) = 0
      IWFR = LAST + 1_8
      IF (KEEP264 .EQ. 0) THEN
       IF ((SYM.EQ.0).AND.(KEEP265.EQ.-1)) THEN
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     &                           .OR.(J.LT.1)) CYCLE
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
          IF (I.NE.J) THEN
              IF ((J.GE.1).AND.(I.LE.N)) THEN
                IW(IQ(I)) = J
                IQ(I)     = IQ(I) + 1
              ENDIF
          ENDIF
        ENDDO
       ELSE IF (KEEP265.EQ.1) THEN
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     &                           .OR.(J.LT.1)) CYCLE
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
          IF (I.NE.J) THEN
              IF ((J.GE.1).AND.(I.LE.N)) THEN
                IW(IQ(J)) = I
                IQ(J)     = IQ(J) + 1
                IW(IQ(I)) = J
                IQ(I)     = IQ(I) + 1
              ENDIF
          ENDIF
        ENDDO
       ELSE 
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     &                           .OR.(J.LT.1)) CYCLE
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
          IF (I.NE.J) THEN
            IF (I.LT.J) THEN
              IF ((I.GE.1).AND.(J.LE.N)) THEN
                IW(IQ(I)) = -J
                IQ(I)     = IQ(I) + 1 
              ENDIF
            ELSE
              IF ((J.GE.1).AND.(I.LE.N)) THEN
                IW(IQ(J)) = -I
                IQ(J)     = IQ(J) + 1
              ENDIF
            ENDIF
          ENDIF
        ENDDO
       ENDIF 
      ELSE
       IF ((SYM.EQ.0).AND.(KEEP265.EQ.-1)) THEN
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
          IF (I.NE.J) THEN
               IW(IQ(I)) = J
               IQ(I)     = IQ(I) + 1
          ENDIF
        ENDDO
       ELSE IF (KEEP265.EQ.1) THEN
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
          IF (I.NE.J) THEN
               IW(IQ(J)) = I
               IQ(J)     = IQ(J) + 1
               IW(IQ(I)) = J
               IQ(I)     = IQ(I) + 1
          ENDIF
        ENDDO
       ELSE
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  
          IF (I.NE.J) THEN
            IF (I.LT.J) THEN
              IW(IQ(I)) = -J
              IQ(I)     = IQ(I) + 1 
            ELSE
              IW(IQ(J)) = -I
              IQ(J)     = IQ(J) + 1
            ENDIF
          ENDIF
        ENDDO
       ENDIF
      ENDIF
      IF (KEEP265.EQ.0) THEN
       NDUP = 0_8
       DO I=1,N
        K1 = IPE(I) 
        K2 = IQ(I) - 1_8
        IF (K1.GT.K2) THEN
         LEN(I) = 0
        ELSE
         DO K8=K1,K2
           J     = -IW(K8)
           IF (J.LE.0) EXIT
           IF (FLAG(J).EQ.I) THEN
            NDUP = NDUP + 1_8
            IW(K8) = 0
           ELSE
            L       = IQ(J) 
            IQ(J)   = L + 1
            IW(L)   = I
            IW(K8)  = J
            FLAG(J) = I
           ENDIF
         END DO  
         LEN(I) = int(IQ(I) - IPE(I))
        ENDIF
       ENDDO
       IF (NDUP.NE.0_8) THEN
        IWFR = 1_8
        DO I=1,N
         IF (LEN(I).EQ.0) THEN
            IPE(I) = IWFR
            CYCLE
         ENDIF
         K1 = IPE(I) 
         K2 = K1 + LEN(I) - 1
         L = IWFR
         IPE(I) = IWFR
         DO 270 K8=K1,K2
           IF (IW(K8).NE.0) THEN
            IW(IWFR) = IW(K8)
            IWFR     = IWFR + 1_8
           ENDIF
  270    CONTINUE
         LEN(I) = int(IWFR - L)
        ENDDO
       ENDIF
       IPE(N+1) = IPE(N) + int(LEN(N),8)
       IWFR = IPE(N+1)
      ELSE 
       IPE(1) = 1_8
       DO I = 1, N
        LEN(I) = int(IQ(I) - IPE(I))
       ENDDO
       DO I = 1, N
        IPE(I+1) = IPE(I) + int(LEN(I),8)
       ENDDO
       IWFR = IPE(N+1)
      ENDIF  
      symmetry = 100  
      IF (SYM.EQ.0) THEN
      RSYM =  dble(NDIAGA+2_8*NZOFFA - (IWFR-1_8))/
     &            dble(NZOFFA+NDIAGA) 
      IF ((KEEP265.EQ.0) .AND. (NZOFFA - (IWFR-1_8)).EQ.0_8) THEN
      ENDIF
      symmetry = nint (100.0D0*RSYM)
         IF ((MPG .GT. 0).AND.(ICNTL(4).GE.2) )
     &  write(MPG,'(A,A,I5)') 
     & ' Case of Schur:',
     & ' structural symmetry (in percent) of interior block=', 
     &   symmetry
        IF (MP.GT.0 .AND. MPG.NE.MP.AND. (ICNTL(4).GE.2) )
     &  write(MP,'(A,A,I5)') 
     & ' Case of Schur:',
     & ' structural symmetry (in percent) of interior block=', 
     &   symmetry
      ELSE
       symmetry = 100
      ENDIF
      INPLACE64_GRAPH_COPY = (LW.GE.2*(IWFR-1))
      AvgDens = nint(dble(IWFR-1_8)/dble(N))
      THRESH  = AvgDens*50 - AvgDens/10 + 1
      NBQD    = 0
      IF (N.GT.2) THEN
        DO I= 1, N
          J = max(LEN(I),1)
          IF (J.GT.THRESH) NBQD = NBQD+1
        ENDDO
      ENDIF
      IF (MPG .GT. 0.AND.(ICNTL(4).GE.2))
     &  write(MPG,'(A,1I5)') 
     &  ' Average density of rows/columns =', AvgDens
        IF (MP.GT.0 .AND. MPG.NE.MP.AND.(ICNTL(4).GE.2))
     &  write(MPG,'(A,1I5)') 
     &  ' Average density of rows/columns =', AvgDens
      RETURN
99999 FORMAT (/'*** Warning message from analysis routine ***')
      END SUBROUTINE DMUMPS_GNEW_SCHUR
      SUBROUTINE DMUMPS_GET_PERM_FROM_PE(N,PE,INVPERM,NFILS,WORK)
      IMPLICIT NONE
      INTEGER N
      INTEGER PE(N),INVPERM(N),NFILS(N),WORK(N)
      INTEGER I,FATHER,STKLEN,STKPOS,PERMPOS,CURVAR
      NFILS = 0
      DO I=1,N
         FATHER = -PE(I)
         IF(FATHER .NE. 0) NFILS(FATHER) = NFILS(FATHER) + 1
      ENDDO
      STKLEN = 0
      PERMPOS = 1
      DO I=1,N
         IF(NFILS(I) .EQ. 0) THEN
            STKLEN = STKLEN + 1
            WORK(STKLEN) = I
            INVPERM(I) = PERMPOS
            PERMPOS = PERMPOS + 1
         ENDIF
      ENDDO
      DO STKPOS = 1,STKLEN
         CURVAR = WORK(STKPOS)
         FATHER = -PE(CURVAR)
         DO
            IF(FATHER .EQ. 0) EXIT
            IF(NFILS(FATHER) .EQ. 1) THEN
               INVPERM(FATHER) = PERMPOS
               FATHER = -PE(FATHER)
               PERMPOS = PERMPOS + 1
            ELSE
               NFILS(FATHER) = NFILS(FATHER) - 1
               EXIT
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_GET_PERM_FROM_PE
      SUBROUTINE DMUMPS_GET_ELIM_TREE(N,PE,NV,WORK)
      IMPLICIT NONE
      INTEGER N
      INTEGER PE(N),NV(N),WORK(N)
      INTEGER I,FATHER,LEN,NEWSON,NEWFATHER
      DO I=1,N
         IF(NV(I) .GT. 0) CYCLE
         LEN = 1
         WORK(LEN) = I
         FATHER = -PE(I)
         DO
            IF(NV(FATHER) .GT. 0) THEN
               NEWSON = FATHER
               EXIT
            ENDIF
            LEN = LEN + 1
            WORK(LEN) = FATHER
            NV(FATHER) = 1
            FATHER = -PE(FATHER)
         ENDDO
         NEWFATHER = -PE(FATHER)
         PE(WORK(LEN)) = -NEWFATHER
         PE(NEWSON) = -WORK(1)
      ENDDO      
      END SUBROUTINE DMUMPS_GET_ELIM_TREE
