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
      SUBROUTINE DMUMPS_ELT_ASM_S_2_S_INIT(
     &    NELT, FRT_PTR, FRT_ELT,
     &    N, INODE, IW, LIW, A, LA, 
     &    NBROWS, NBCOLS,
     &    OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     &    RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR, 
     &    ICNTL, KEEP, KEEP8, MYID, LRGROUPS)
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR
      IMPLICIT NONE
      INTEGER NELT, N,LIW
      INTEGER(8) :: LA
      INTEGER KEEP(500), ICNTL(60)
      INTEGER(8) KEEP8(150)
      INTEGER INODE, MYID
      INTEGER NBROWS, NBCOLS 
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER IW(LIW), ITLOC(N+KEEP(253)), STEP(N),
     &        PTRIST(KEEP(28)), FILS(N)
      INTEGER(8), INTENT(IN) :: PTRARW(NELT+1), PTRAIW(NELT+1)
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(KEEP8(27))
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      DOUBLE PRECISION :: A(LA)
      DOUBLE PRECISION :: DBLARR(KEEP8(26))
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER, INTENT(IN)    :: LRGROUPS(N)
      INTEGER(8) :: POSELT
      DOUBLE PRECISION, DIMENSION(:), POINTER :: A_PTR
      INTEGER(8) :: LA_PTR
      INTEGER IOLDPS, NBCOLF, NBROWF, NSLAVES, HF,
     &        K1,K2,K,J,JPOS,NASS
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INCLUDE 'mumps_headers.h'
      IOLDPS  = PTRIST(STEP(INODE))
      CALL DMUMPS_DM_SET_DYNPTR( IW(IOLDPS+XXS), A, LA,
     &     PTRAST(STEP(INODE)), IW(IOLDPS+XXD), IW(IOLDPS+XXR),
     &     A_PTR, POSELT, LA_PTR )
      NBCOLF  = IW(IOLDPS+KEEP(IXSZ))
      NBROWF  = IW(IOLDPS+2+KEEP(IXSZ))
      NASS    = IW(IOLDPS+1+KEEP(IXSZ))
      NSLAVES = IW(IOLDPS+5+KEEP(IXSZ))
      HF      = 6 + NSLAVES + KEEP(IXSZ)
      IF (NASS.LT.0) THEN
          NASS         = -NASS
          IW(IOLDPS+1+KEEP(IXSZ)) = NASS
          CALL DMUMPS_ASM_SLAVE_ELEMENTS( INODE, N, NELT, IW, LIW,
     &    IOLDPS, A_PTR(POSELT), LA_PTR, 1_8, KEEP, KEEP8, ITLOC, FILS,
     &    PTRAIW, PTRARW,
     &    INTARR, DBLARR, KEEP8(27), KEEP8(26), FRT_PTR, FRT_ELT,
     &    RHS_MUMPS, LRGROUPS)
      ENDIF
      IF (NBROWS.GT.0) THEN
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = JPOS
           JPOS     = JPOS + 1
          ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_ELT_ASM_S_2_S_INIT
      SUBROUTINE DMUMPS_ASM_SLAVE_ELEMENTS( INODE, N, NELT, IW, LIW,
     &IOLDPS, A, LA, POSELT, KEEP, KEEP8, ITLOC, FILS, PTRAIW, PTRARW,
     &INTARR, DBLARR, LINTARR, LDBLARR, FRT_PTR, FRT_ELT, RHS_MUMPS,
     &LRGROUPS)
!$    USE OMP_LIB
      USE DMUMPS_ANA_LR,    ONLY : GET_CUT
      USE DMUMPS_LR_CORE,   ONLY : MAX_CLUSTER
      USE MUMPS_LR_COMMON,  ONLY : COMPUTE_BLR_VCS
      IMPLICIT NONE
      INTEGER, intent(in)    :: N, NELT, LIW, IOLDPS, INODE
      INTEGER(8), intent(in) :: LA, POSELT, LINTARR, LDBLARR
      INTEGER, intent(in)    :: IW(LIW)
      INTEGER, intent(in)    :: KEEP(500)
      INTEGER(8), intent(in) :: KEEP8(150)
      INTEGER, intent(inout) :: ITLOC(N+KEEP(253))
      DOUBLE PRECISION, intent(inout) :: A(LA)
      DOUBLE PRECISION, intent(in)    :: RHS_MUMPS(KEEP(255))
      INTEGER, intent(in)    :: INTARR(LINTARR)
      DOUBLE PRECISION, intent(in)    :: DBLARR(LDBLARR)
      INTEGER, intent(in)    :: FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER, intent(in)    :: FILS(N)
      INTEGER(8), intent(in) :: PTRAIW(NELT+1), PTRARW(NELT+1)
      INTEGER, INTENT(IN)    :: LRGROUPS(N)
!$    INTEGER :: NOMP
!$    INTEGER(8) :: CHUNK8  
!$    INTEGER    :: CHUNK
      INCLUDE 'mumps_headers.h'
      INTEGER    :: HF, NBROWF, NBCOLF, NASS, NSLAVES
      INTEGER    :: ILOC, IELL, ELTI, ELBEG, NUMELT
      INTEGER(8) :: SIZE_ELTI8 
      INTEGER    :: I, J, K, K1, K2
      INTEGER    :: IPOS, IPOS1, IPOS2, JPOS, IJROW
      INTEGER    :: IN
      INTEGER(8) :: II8, JJ8, J18, J28 
      INTEGER(8) :: AINPUT8
      INTEGER(8) :: AII8
      INTEGER(8) :: APOS, APOS2, ICT12
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_LS
      INTEGER :: NB_BLR_LS, NPARTSCB, NPARTSASS, MAXI_CLUSTER, 
     &           IBCKSZ2, MINSIZE, TOPDIAG
      INTEGER(8) :: JJ3
      INTEGER    :: K1RHS, K2RHS, JFirstRHS
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      NBCOLF  = IW(IOLDPS+KEEP(IXSZ))
      NBROWF  = IW(IOLDPS+2+KEEP(IXSZ))
      NASS    = IW(IOLDPS+1+KEEP(IXSZ))
      NSLAVES= IW(IOLDPS+5 + KEEP(IXSZ))
      HF      = 6 + NSLAVES + KEEP(IXSZ)
!$    NOMP = OMP_GET_MAX_THREADS()
      IF (KEEP(50) .EQ. 0 .OR. NBROWF .LT. KEEP(63)) THEN
!$      CHUNK8 = int(KEEP(361),8)
!$OMP   PARALLEL DO PRIVATE(JJ8) SCHEDULE(STATIC, CHUNK8)
!$OMP&  IF (int(NBROWF,8)*int(NBCOLF,8) > int(KEEP(361),8)
!$OMP&    .AND. NOMP .GT. 1)
        DO JJ8=POSELT, POSELT+int(NBROWF,8)*int(NBCOLF,8)-1_8
          A(JJ8) = ZERO
        ENDDO
!$OMP   END PARALLEL DO
      ELSE
        TOPDIAG = 0
        IF (IW(IOLDPS+XXLR).GE.1) THEN
          CALL GET_CUT(IW(IOLDPS+HF:IOLDPS+HF+NBROWF-1), 0,
     &                    NBROWF, LRGROUPS, NPARTSCB, 
     &                    NPARTSASS, BEGS_BLR_LS)
          NB_BLR_LS = NPARTSCB
          call MAX_CLUSTER(BEGS_BLR_LS,NB_BLR_LS+1,MAXI_CLUSTER)
          DEALLOCATE(BEGS_BLR_LS)
          CALL COMPUTE_BLR_VCS(KEEP(472), IBCKSZ2, KEEP(488), NASS)
          MINSIZE = int(IBCKSZ2 / 2)
          TOPDIAG = max(2*MINSIZE + MAXI_CLUSTER-1, TOPDIAG)
        ENDIF
!$      CHUNK = max( KEEP(360)/2,
!$   &               ((NBROWF+NOMP-1)/NOMP +2) / 3 )
!$OMP   PARALLEL DO PRIVATE(APOS,JJ3,JJ8) SCHEDULE(STATIC,CHUNK)
!$OMP&  IF (NBROWF .GT. KEEP(360) .AND. NOMP .GT. 1)
        DO JJ8 = 0_8, int(NBROWF-1,8)
          APOS = POSELT+ JJ8*int(NBCOLF,8)
          JJ3  = min( int(NBCOLF,8)  - 1_8, 
     &           JJ8 + int(NBCOLF-NBROWF,8) + TOPDIAG )
          A(APOS: APOS+JJ3) = ZERO
        ENDDO
!$OMP   END PARALLEL DO
      ENDIF
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = -JPOS
           JPOS     = JPOS + 1
          END DO
          K1 = IOLDPS + HF 
          K2 = K1 + NBROWF - 1
          JPOS = 1
          IF ((KEEP(253).GT.0).AND.(KEEP(50).NE.0)) THEN
           K1RHS = 0
           K2RHS = -1
           DO K = K1, K2
            J        = IW(K)
            ITLOC(J) = -ITLOC(J)*NBCOLF + JPOS
            IF ((K1RHS.EQ.0).AND.(J.GT.N)) THEN
             K1RHS = K
             JFirstRHS=J-N 
            ENDIF
            JPOS     = JPOS + 1
           ENDDO
           IF (K1RHS.GT.0) K2RHS=K2
           IF ( K2RHS.GE.K1RHS ) THEN
             IN = INODE
             DO WHILE (IN.GT.0) 
               IJROW = -ITLOC(IN)  
               DO K = K1RHS, K2RHS
                J    = IW(K)       
                I    = ITLOC(J)    
                ILOC = mod(I,NBCOLF) 
              APOS = POSELT+int(ILOC-1,8)*int(NBCOLF,8) + 
     &               int(IJROW-1,8) 
              A(APOS) = A(APOS) + RHS_MUMPS(
     &                 (JFirstRHS+(K-K1RHS)-1)*KEEP(254)+ IN)
             ENDDO
             IN = FILS(IN)
            ENDDO
           ENDIF
          ELSE  
           DO K = K1, K2
            J        = IW(K)
            ITLOC(J) = -ITLOC(J)*NBCOLF + JPOS
            JPOS     = JPOS + 1
           END DO
          ENDIF
          ELBEG  = FRT_PTR(INODE)
          NUMELT = FRT_PTR(INODE+1) - ELBEG
          DO IELL=ELBEG,ELBEG+NUMELT-1
           ELTI = FRT_ELT(IELL)
           J18= PTRAIW(ELTI)
           J28= PTRAIW(ELTI+1)-1_8
           AII8 = PTRARW(ELTI)
           SIZE_ELTI8 = J28 - J18 + 1_8
           DO II8=J18,J28
            I = ITLOC(INTARR(II8))
            IF (KEEP(50).EQ.0) THEN
             IF (I.LE.0) CYCLE
             AINPUT8    = AII8 + II8 - J18
             IPOS = mod(I,NBCOLF)
             ICT12 = POSELT + int(IPOS-1,8) * int(NBCOLF,8)
             DO JJ8 = J18, J28
              JPOS = ITLOC(INTARR(JJ8))
              IF (JPOS.LE.0) THEN 
                   JPOS = -JPOS
              ELSE
                   JPOS = JPOS/NBCOLF
              END IF
              APOS2    = ICT12 + int(JPOS - 1,8)
              A(APOS2) = A(APOS2) +  DBLARR(AINPUT8)
              AINPUT8   = AINPUT8 + SIZE_ELTI8
             END DO
            ELSE
              IF ( I .EQ. 0 ) THEN 
               AII8 = AII8 + J28 - II8 + 1_8
               CYCLE
              ENDIF
              IF ( I .LE. 0 ) THEN 
               IPOS1 = -I
               IPOS2 = 0
              ELSE 
               IPOS1 = I/NBCOLF
               IPOS2 = mod(I,NBCOLF)
              END IF
              ICT12 =  POSELT + int(IPOS2-1,8)*int(NBCOLF,8)
              DO JJ8=II8,J28
               AII8 = AII8 + 1_8
               J = ITLOC(INTARR(JJ8))
               IF ( J .EQ. 0 ) CYCLE
               IF ( IPOS2.EQ.0 .AND. J.LE.0) CYCLE
               IF ( J .LE. 0 ) THEN
                JPOS = -J
               ELSE
                JPOS = J/NBCOLF
               END IF
               IF ( (IPOS1.GE.JPOS) .AND. (IPOS2.GT.0) ) THEN
                 APOS2 = ICT12  + int(JPOS - 1,8)
                 A(APOS2) = A(APOS2) +  DBLARR(AII8-1_8)
               END IF
               IF ( (IPOS1.LT.JPOS) .AND. (J.GT.0) ) THEN
                 IPOS = mod(J,NBCOLF)
                 JPOS = IPOS1
                 APOS2 = POSELT + int(IPOS-1,8)*int(NBCOLF,8)
     &                          + int(JPOS - 1,8)
                 A(APOS2) = A(APOS2) +  DBLARR(AII8-1_8)
               END IF
              END DO
            END IF
           END DO
          END DO
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          DO K = K1, K2
           J = IW(K)
           ITLOC(J) = 0
          END DO
      END SUBROUTINE DMUMPS_ASM_SLAVE_ELEMENTS
