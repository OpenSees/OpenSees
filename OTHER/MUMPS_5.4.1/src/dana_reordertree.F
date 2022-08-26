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
      SUBROUTINE DMUMPS_REORDER_TREE(N,FRERE, STEP, FILS,
     &     NA,LNA,NE,ND, DAD, LDAD, USE_DAD,
     &     NSTEPS,PERM,SYM,INFO,LP,K215,K234,K55,K199,
     &     PROCNODE,SLAVEF, PEAK,SBTR_WHICH_M
     &     )
      IMPLICIT NONE
      INTEGER N,PERM,SYM, NSTEPS, LNA, LP,LDAD
      INTEGER FRERE(NSTEPS), FILS(N), STEP(N)
      INTEGER NA(LNA), NE(NSTEPS), ND(NSTEPS)
      INTEGER K215,K234,K55,K199
      INTEGER DAD(LDAD)
      LOGICAL USE_DAD
      INTEGER INFO(80)
      INTEGER SLAVEF,PROCNODE(NSTEPS)
      INTEGER :: SBTR_WHICH_M
      EXTERNAL MUMPS_PROCNODE
      INTEGER MUMPS_PROCNODE
      DOUBLE PRECISION PEAK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COST_TRAV
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH
      INTEGER IFATH,IN,NSTK,INODE,I,allocok,LOCAL_PERM
      INTEGER(8) NCB
      INTEGER(8) NELIM,NFR
      INTEGER NFR4,NELIM4
      INTEGER LEAF,NBLEAF,NBROOT, SIZE_TAB
      INTEGER, DIMENSION (:), ALLOCATABLE :: IPOOL,TNSTK
      INTEGER, DIMENSION (:), ALLOCATABLE,TARGET :: SON,TEMP
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: M,M_TOTAL, fact
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: TAB1,TAB2
      INTEGER, DIMENSION (:), POINTER :: TAB
      INTEGER dernier,fin
      INTEGER cour,II
      INTEGER CB_current,CB_MAX,ROOT_OF_CUR_SBTR
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: T1,T2
      INTEGER, DIMENSION (:), ALLOCATABLE :: RESULT
      INTEGER(8) MEM_SIZE,FACT_SIZE,SUM,MEM_SEC_PERM,FACT_SIZE_T,
     &     MEM_SIZE_T,TOTAL_MEM_SIZE,TMP_TOTAL_MEM_SIZE,TMP_SUM,
     &     SIZECB, SIZECB_LASTSON
      INTEGER(8) TMP8
      LOGICAL   SBTR_M
      INTEGER FIRST_LEAF,SIZE_SBTR
      EXTERNAL MUMPS_IN_OR_ROOT_SSARBR,MUMPS_INSSARBR
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR,MUMPS_INSSARBR
      DOUBLE PRECISION COST_NODE
      INCLUDE 'mumps_headers.h'
      TOTAL_MEM_SIZE=0_8
      ROOT_OF_CUR_SBTR=0
      IF((PERM.EQ.0).OR.(PERM.EQ.1).OR.
     &     (PERM.EQ.2).OR.(PERM.EQ.3).OR.(PERM.EQ.4).OR.
     &     (PERM.EQ.5).OR.(PERM.EQ.6))THEN
         LOCAL_PERM=0
      ENDIF
      SBTR_M=.FALSE.
      MEM_SIZE=0_8
      FACT_SIZE=0_8
      IF ((PERM.LT.0 .OR. PERM.GT.7)) THEN
         WRITE(*,*) "Internal Error in DMUMPS_REORDER_TREE",PERM
         CALL MUMPS_ABORT()
      END IF
      NBLEAF = NA(1)
      NBROOT = NA(2)
      IF((PERM.EQ.0).AND.(NBROOT.EQ.NBLEAF)) RETURN
      IF ((PERM.NE.7).AND.(SBTR_M.OR.(PERM.EQ.2))) THEN
         IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
            ALLOCATE(M_TOTAL(NSTEPS), stat=allocok )
            IF (allocok > 0) THEN
               IF ( LP .GT. 0 )
     &              WRITE(LP,*)'Memory allocation error in
     &              DMUMPS_REORDER_TREE'
               INFO(1)=-7
               INFO(2)=NSTEPS
               RETURN
            ENDIF
         ENDIF
      ENDIF
      IF(PERM.NE.7)THEN
         ALLOCATE(M(NSTEPS),stat=allocok )
         IF (allocok > 0) THEN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*)'Memory allocation error
     &in DMUMPS_REORDER_TREE'
            INFO(1)=-7
            INFO(2)=NSTEPS
            RETURN
         ENDIF
      ENDIF
      ALLOCATE( IPOOL(NBLEAF), fact(NSTEPS),TNSTK(NSTEPS),
     &     stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_REORDER_TREE'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      II=0
      DO I=1,NSTEPS
         TNSTK(I) = NE(I)
         IF(NE(I).GE.II) II=NE(I)
      ENDDO
      SIZE_TAB=max(II,NBROOT)
      ALLOCATE(SON(II), TEMP(II),
     &         TAB1(SIZE_TAB), TAB2(SIZE_TAB), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_REORDER_TREE'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      ALLOCATE(T1(SIZE_TAB),T2(SIZE_TAB),
     &         RESULT(SIZE_TAB),stat=allocok)
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_REORDER_TREE'
        INFO(1)=-7
        INFO(2)=SIZE_TAB
        RETURN
      ENDIF
      IF(PERM.EQ.7) THEN
         GOTO 001
      ENDIF
      IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
         ALLOCATE(COST_TRAV(NSTEPS), stat=allocok )
         IF (allocok > 0) THEN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*)'Memory allocation error
     &       in DMUMPS_REORDER_TREE'
            INFO(1)=-7
            INFO(2)=NSTEPS
            RETURN
         ENDIF
         COST_TRAV=0.0D0
         COST_NODE=0.0d0
      ENDIF
      IF(NBROOT.EQ.NBLEAF)THEN
         IF((PERM.NE.1).OR.(PERM.EQ.4).OR.(PERM.EQ.6))THEN
          WRITE(*,*)'Internal Error in reordertree:'
          WRITE(*,*)'  problem with perm parameter in reordertree'
          CALL MUMPS_ABORT()
        ENDIF
        DO I=1,NBROOT
          TAB1(I)=int(ND(STEP(NA(I+2+NBLEAF))),8)
          IPOOL(I)=NA(I+2+NBLEAF)
          M(STEP(IPOOL(I)))=TAB1(I)*TAB1(I)
        ENDDO
        CALL DMUMPS_FUSION_SORT(NA(2+NBLEAF+1),NBROOT,TAB1,TAB2,4,
     &    RESULT,T1,T2)
        GOTO 789
      ENDIF
      IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
         ALLOCATE(DEPTH(NSTEPS),stat=allocok)
         IF (allocok > 0) THEN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*)'Memory allocation error in
     &           DMUMPS_REORDER_TREE'
            INFO(1)=-7
            INFO(2)=NSTEPS
            RETURN
         ENDIF
         DEPTH=0
         NBROOT = NA(2)
         IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
         fin=NBROOT
         LEAF=NA(1)
 499     CONTINUE
         INODE=IPOOL(fin)
         IF(INODE.LT.0)THEN
            WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
            CALL MUMPS_ABORT()
         ENDIF
         IN=INODE
 4602    IN = FILS(IN)
         IF (IN .GT. 0 ) THEN
            GOTO 4602
         ENDIF
         IN=-IN
         DO I=1,NE(STEP(INODE))
            SON(I)=IN
            IN=FRERE(STEP(IN))
         ENDDO
         DO I=1,NE(STEP(INODE))
            IPOOL(fin)=SON(I)
            DEPTH(STEP(SON(I)))=DEPTH(STEP(INODE))+1
            SON(I)=0
            fin=fin+1
         ENDDO
         IF(NE(STEP(INODE)).EQ.0)THEN
            LEAF=LEAF-1
         ELSE
            fin=fin-1
            GOTO 499
         ENDIF
         fin=fin-1
         IF(fin.EQ.0) GOTO 489
         GOTO 499
 489     CONTINUE
      ENDIF
      DO I=1,NSTEPS
         M(I)=0_8
         IF (SBTR_M.OR.(PERM.EQ.2))  THEN
            IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
               M_TOTAL(I)=0_8
            ENDIF
         ENDIF
      ENDDO
      DO I=1,NSTEPS
         fact(I)=0_8
      ENDDO
      IPOOL(1:NBLEAF)=NA(3:2+NBLEAF)
      LEAF = NBLEAF + 1
 91   CONTINUE
        IF (LEAF.NE.1) THEN
           LEAF = LEAF -1
           INODE = IPOOL(LEAF)
        ENDIF
 96     CONTINUE
        NFR    = int(ND(STEP(INODE)),8)
        NSTK   = NE(STEP(INODE))
        NELIM4 = 0
        IN = INODE
 101    NELIM4 = NELIM4 + 1
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 101
        NELIM=int(NELIM4,8)
        IF(NE(STEP(INODE)).EQ.0) THEN
           M(STEP(INODE))=NFR*NFR
           IF (SBTR_M.OR.(PERM.EQ.2))  THEN
                 M_TOTAL(STEP(INODE))=NFR*NFR
           ENDIF
        ENDIF
        IF((PERM.EQ.4).OR.(PERM.EQ.3))THEN
           IF(MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &         K199))THEN
              DEPTH(STEP(INODE))=0
           ENDIF
        ENDIF
        IF ( SYM .eq. 0 ) THEN
          fact(STEP(INODE))=fact(STEP(INODE))+
     &      (2_8*NFR*NELIM)-(NELIM*NELIM)
        ELSE
          fact(STEP(INODE))=fact(STEP(INODE))+NFR*NELIM
        ENDIF
        IF (USE_DAD) THEN
          IFATH = DAD( STEP(INODE) )
        ELSE
          IN = INODE
 113      IN = FRERE(IN)
          IF (IN.GT.0) GO TO 113
          IFATH = -IN
        ENDIF
        IF (IFATH.EQ.0) THEN
           NBROOT = NBROOT - 1
           IF (NBROOT.EQ.0) GOTO 116
           GOTO 91
        ELSE
           fact(STEP(IFATH))=fact(STEP(IFATH))+fact(STEP(INODE))
           IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
              DEPTH(STEP(IFATH))=max(DEPTH(STEP(INODE)),
     &             DEPTH(STEP(IFATH)))
           ENDIF
        ENDIF
        TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
        IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
           INODE = IFATH        
           IN=INODE
           dernier=IN
           I=1
 5700      IN = FILS(IN)
           IF (IN .GT. 0 ) THEN
             dernier=IN
             I=I+1
             GOTO 5700
           ENDIF
           NCB=int(ND(STEP(INODE))-I,8)
           IN=-IN
           IF(PERM.NE.7)THEN
              DO I=1,NE(STEP(INODE))
                 SON(I)=IN
                 TEMP(I)=IN
                 IF(IN.GT.0) IN=FRERE(STEP(IN))
              ENDDO
           ELSE
              DO I=NE(STEP(INODE)),1,-1
                 SON(I)=IN
                 TEMP(I)=IN
                 IF(IN.GT.0) IN=FRERE(STEP(IN))
              ENDDO
           ENDIF
           NFR = int(ND(STEP(INODE)),8)
           DO II=1,NE(STEP(INODE))
             TAB1(II)=0_8
             TAB2(II)=0_8
             cour=SON(II)
             NELIM4=1
 151         cour=FILS(cour)
             IF(cour.GT.0) THEN
                NELIM4=NELIM4+1
                GOTO 151
             ENDIF
             NELIM=int(NELIM4,8)
             IF((SYM.EQ.0).OR.(K215.NE.0)) THEN
                SIZECB=(int(ND(STEP(SON(II))),8)-NELIM)
     &                *(int(ND(STEP(SON(II))),8)-NELIM)
             ELSE
                SIZECB=(int(ND(STEP(SON(II))),8)-NELIM)
     &                *(int(ND(STEP(SON(II))),8)-
     &               NELIM+1_8)/2_8
             ENDIF
             IF((PERM.EQ.0).OR.(PERM.EQ.5))THEN
                IF (K234 .NE. 0 .AND. K55.EQ.0 ) THEN
                   TMP8=NFR
                   TMP8=TMP8*TMP8
                   TAB1(II)=max(TMP8, M(STEP(SON(II)))) - SIZECB
                   TAB2(II)=SIZECB
                ELSE
                   TAB1(II)=M(STEP(SON(II)))- SIZECB
                   TAB2(II)=SIZECB
                ENDIF
             ENDIF
             IF((PERM.EQ.1).OR.(PERM.EQ.6)) THEN
                TAB1(II)=M(STEP(SON(II)))-SIZECB
                TAB1(II)=TAB1(II)-fact(STEP(SON(II)))
                TAB2(II)=SIZECB+fact(STEP(SON(II)))
             ENDIF
             IF(PERM.EQ.2)THEN
                IF (MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &               K199))THEN
                   TAB1(II)=M_TOTAL(STEP(SON(II)))-SIZECB
     &                  -fact(STEP(SON(II)))
                   TAB2(II)=SIZECB
                ELSE
                   TAB1(II)=M(STEP(SON(II)))-SIZECB
                   TAB2(II)=SIZECB             
                ENDIF
             ENDIF
             IF(PERM.EQ.3)THEN
                IF (MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &               K199))THEN
                   TAB1(II)=M(STEP(SON(II)))-SIZECB
                   TAB2(II)=SIZECB               
                ELSE
                   TAB1(II)=int(DEPTH(STEP(SON(II))),8)
                   TAB2(II)=M(STEP(SON(II)))
                ENDIF
             ENDIF
             IF(PERM.EQ.4)THEN
                IF (MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &               K199))THEN
                   TAB1(II)=M(STEP(SON(II)))-
     &                  SIZECB-fact(STEP(SON(II)))
                   TAB2(II)=SIZECB             
                ELSE
                   TAB1(II)=int(DEPTH(STEP(SON(II))),8)
                   TAB2(II)=M(STEP(SON(II)))
                ENDIF
             ENDIF
          ENDDO
          CALL DMUMPS_FUSION_SORT(SON,NE(STEP(INODE)),TAB1,TAB2,
     &         LOCAL_PERM
     &           ,RESULT,T1,T2)
          IF(PERM.EQ.0) THEN
             DO II=1,NE(STEP(INODE))
               cour=TEMP(II)
               NELIM4=1
 153           cour=FILS(cour)
               IF(cour.GT.0) THEN
                  NELIM4=NELIM4+1
                  GOTO 153
               ENDIF
               NELIM=int(NELIM4,8)
               IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                  SIZECB=(int(ND(STEP(TEMP(II))),8)-NELIM)*
     &                 (int(ND(STEP(TEMP(II))),8)-NELIM)
               ELSE
                  SIZECB=(int(ND(STEP(TEMP(II))),8)-NELIM)*
     &                    (int(ND(STEP(TEMP(II))),8)-NELIM+1_8)/2_8
               ENDIF
               TAB1(II)=SIZECB
             ENDDO
             CALL DMUMPS_FUSION_SORT(TEMP,NE(STEP(INODE)),TAB1,TAB2,3,
     &         RESULT,T1,T2)             
           ENDIF
           IF(PERM.EQ.1) THEN
              DO II=1,NE(STEP(INODE))
                cour=TEMP(II)
                NELIM4=1
 187            cour=FILS(cour)
                IF(cour.GT.0) THEN
                   NELIM4=NELIM4+1
                   GOTO 187
                ENDIF    
                NELIM=int(NELIM4,8)
                IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                   SIZECB=(int(ND(STEP(TEMP(II))),8)-NELIM)*
     &                    (int(ND(STEP(TEMP(II))),8)-NELIM)
                ELSE
                   SIZECB=(int(ND(STEP(TEMP(II))),8)-NELIM)*
     &                    (int(ND(STEP(TEMP(II))),8)-NELIM+1_8)/2_8
                ENDIF
                TAB1(II)=SIZECB+fact(STEP(TEMP(II)))
             ENDDO
             CALL DMUMPS_FUSION_SORT(TEMP,NE(STEP(INODE)),TAB1,TAB2,3,
     &         RESULT,T1,T2)
           ENDIF
           CONTINUE
           IFATH=INODE
           DO II=1,2
              SUM=0_8
              FACT_SIZE=0_8
              FACT_SIZE_T=0_8
              MEM_SIZE=0_8
              MEM_SIZE_T=0_8
              CB_MAX=0
              CB_current=0
              TMP_SUM=0_8
              IF(II.EQ.1) TAB=>SON 
              IF(II.EQ.2) TAB=>TEMP
              DO I=1,NE(STEP(INODE))
                 cour=TAB(I)
                 NELIM4=1
 149             cour=FILS(cour)
                 IF(cour.GT.0) THEN
                    NELIM4=NELIM4+1
                    GOTO 149
                 ENDIF    
                 NELIM=int(NELIM4, 8)
                 NFR=int(ND(STEP(TAB(I))),8)
                 IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                    SIZECB=(NFR-NELIM)*(NFR-NELIM)
                 ELSE
                    SIZECB=(NFR-NELIM)*(NFR-NELIM+1_8)/2_8
                 ENDIF
                 MEM_SIZE=max(MEM_SIZE,(M(STEP(TAB(I)))+SUM+FACT_SIZE))
                 IF (SBTR_M.OR.(PERM.EQ.2)) THEN
                       MEM_SIZE_T=max(MEM_SIZE_T,(M_TOTAL(STEP(TAB(I)))+
     &                      SUM+
     &                      FACT_SIZE_T))
                       FACT_SIZE_T=FACT_SIZE_T+fact(STEP(TAB(I)))
                 ENDIF
                 TOTAL_MEM_SIZE=max(TOTAL_MEM_SIZE,
     &                (M(STEP(TAB(I)))+SUM+FACT_SIZE))
                 TMP_SUM=TMP_SUM+fact(STEP(TAB(I)))
                 SUM=SUM+SIZECB
                 SIZECB_LASTSON = SIZECB
                 IF((PERM.EQ.1).OR.(PERM.EQ.4))THEN
                    FACT_SIZE=FACT_SIZE+fact(STEP(TAB(I)))
                 ENDIF
              ENDDO
              IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                 SIZECB=NCB*NCB
              ELSE
                 SIZECB=(NCB*(NCB+1_8))/2_8
              ENDIF
              IF (K234.NE.0 .AND. K55.EQ.0) THEN
                 TOTAL_MEM_SIZE=max(TOTAL_MEM_SIZE,
     &                ( (   int(ND(STEP(IFATH)),8)
     &                    * int(ND(STEP(IFATH)),8) )
     &                  + SUM-SIZECB_LASTSON+TMP_SUM )
     &           )
              ELSE IF (K234.NE.0 .AND. K55.NE.0) THEN
                 TOTAL_MEM_SIZE=max(TOTAL_MEM_SIZE,
     &                ( ( int(ND(STEP(IFATH)),8)
     &                  * int(ND(STEP(IFATH)),8) )
     &                  + SUM + TMP_SUM )
     &           )
              ELSE
                 TOTAL_MEM_SIZE=max(TOTAL_MEM_SIZE,
     &                ( ( int(ND(STEP(IFATH)),8)
     &                  * int(ND(STEP(IFATH)),8))
     &                  + max(SUM,SIZECB) + TMP_SUM )
     &                )
              ENDIF
              IF(II.EQ.1)THEN
                 TMP_TOTAL_MEM_SIZE=TOTAL_MEM_SIZE
              ENDIF
              IF(II.EQ.1)THEN
                 IF (K234.NE.0 .AND. K55.EQ.0) THEN
                   M(STEP(IFATH))=max(MEM_SIZE,((int(ND(STEP(IFATH)),8)
     &             *int(ND(STEP(IFATH)),8))+SUM-SIZECB_LASTSON+
     &             FACT_SIZE))
                 ELSE IF (K234.NE.0 .AND. K55.NE.0) THEN
                   M(STEP(IFATH))=max(MEM_SIZE,((int(ND(STEP(IFATH)),8)
     &             *int(ND(STEP(IFATH)),8))+SUM+FACT_SIZE))
                 ELSE
                   M(STEP(IFATH))=max(MEM_SIZE,((int(ND(STEP(IFATH)),8)
     &             *int(ND(STEP(IFATH)),8))+max(SUM,SIZECB)+FACT_SIZE))
                 ENDIF
                 IF (SBTR_M.OR.(PERM.EQ.2))  THEN
                       M_TOTAL(STEP(IFATH))=max(MEM_SIZE_T,
     &                      ((int(ND(STEP(IFATH)),8)
     &                      *int(ND(STEP(IFATH)),8))+max(SUM,SIZECB)+
     &                      FACT_SIZE_T))
                 ENDIF
              ENDIF
              IF((II.EQ.2).AND.(PERM.EQ.1).OR.(PERM.EQ.0).OR.
     &             (PERM.EQ.5).OR.(PERM.EQ.6).OR.
     &             (.NOT.SBTR_M.OR.(SBTR_WHICH_M.NE.1)))THEN
                 MEM_SEC_PERM=max(MEM_SIZE,((int(ND(STEP(IFATH)),8)
     &             *int(ND(STEP(IFATH)),8))+max(SUM,SIZECB)+FACT_SIZE))
              ENDIF
              IF((PERM.EQ.2).OR.(PERM.EQ.3).OR.(PERM.EQ.4))THEN
                 MEM_SEC_PERM=huge(MEM_SEC_PERM)
              ENDIF
           ENDDO
           IF(MEM_SEC_PERM.EQ.M(STEP(IFATH))) THEN
              TAB=>TEMP
           ELSE IF (MEM_SEC_PERM.LT.M(STEP(IFATH))) THEN
              WRITE(*,*)'Internal error 1 in DMUMPS_REORDER_TREE',
     &        MEM_SEC_PERM, M(STEP(IFATH))
              CALL MUMPS_ABORT()
           ELSE 
              TOTAL_MEM_SIZE=TMP_TOTAL_MEM_SIZE
              TAB=>SON
           ENDIF
           DO I=NE(STEP(INODE)),1,-1
              IF(I.EQ.NE(STEP(INODE))) THEN
                 FILS(dernier)=-TAB(I)
                 dernier=TAB(I)
                 GOTO 222
              ENDIF
              IF(I.EQ.1) THEN
                 FRERE(STEP(dernier))=TAB(I)
                 FRERE(STEP(TAB(I)))=-INODE
                 GOTO 222
              ENDIF
              IF(I.GT.1) THEN
                 FRERE(STEP(dernier))=TAB(I)
                 dernier=TAB(I)
                 GOTO 222
              ENDIF
 222          CONTINUE
           ENDDO
           GOTO 96
        ELSE
           GOTO 91
        ENDIF
 116    CONTINUE
        NBROOT = NA(2)
        IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
        IF (PERM.eq.1) THEN
          DO I=1,NBROOT
            TAB1(I)=M(STEP(NA(I+2+NBLEAF)))-fact(STEP(NA(I+2+NBLEAF)))
            TAB1(I)=-TAB1(I)
          ENDDO
          CALL DMUMPS_FUSION_SORT(NA(2+NBLEAF+1),NBROOT,TAB1,TAB2,4,
     &      RESULT,T1,T2)
          IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
        ENDIF
 001    CONTINUE
        fin=NBROOT
        LEAF=NA(1)
        FIRST_LEAF=-9999
        SIZE_SBTR=0
 999    CONTINUE
        INODE=IPOOL(fin)
        IF(INODE.LT.0)THEN
           WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
           CALL MUMPS_ABORT()
        ENDIF
        IN=INODE
 5602   IN = FILS(IN)
        IF (IN .GT. 0 ) THEN
           dernier=IN
           GOTO 5602
        ENDIF
        IN=-IN
        IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
           IF(SLAVEF.NE.1)THEN
              IF (USE_DAD) THEN
                 IFATH=DAD(INODE)
              ELSE
                IN = INODE
 395            IN = FRERE(IN)
                IF (IN.GT.0) GO TO 395
                IFATH = -IN
              ENDIF
              NFR4   = ND(STEP(INODE))
              NFR    = int(NFR4,8)
              NELIM4 = 0
              IN = INODE
 396          NELIM4 = NELIM4 + 1
              IN = FILS(IN)
              IF (IN .GT. 0 ) GOTO 396
              NELIM=int(NELIM4,8)
              IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                 SIZECB=(NFR-NELIM)*(NFR-NELIM)
              ELSE
                 SIZECB=(NFR-NELIM)*(NFR-NELIM+1_8)/2_8
              ENDIF
              CALL MUMPS_GET_FLOPS_COST(NFR4,NELIM4,NELIM4,
     &             SYM,1,COST_NODE)
              IF(IFATH.NE.0)THEN
                 IF(MUMPS_INSSARBR(PROCNODE(STEP(INODE)),K199))THEN
                    COST_TRAV(STEP(INODE))=COST_TRAV(STEP(
     &                   ROOT_OF_CUR_SBTR))
                 ELSE
                    COST_TRAV(STEP(INODE))=dble(COST_NODE)+
     &                   COST_TRAV(STEP(IFATH))+
     &                   dble(SIZECB*18_8)  
                 ENDIF
              ELSE
                 COST_TRAV(STEP(INODE))=dble(COST_NODE)
              ENDIF
           ENDIF
        ENDIF
        DO I=1,NE(STEP(INODE))
           TEMP(I)=IN
           IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
              IF((SLAVEF.NE.1).AND.(.NOT.MUMPS_IN_OR_ROOT_SSARBR(
     &             PROCNODE(STEP(INODE)),K199)))THEN
                 NFR4   = ND(STEP(INODE))
                 NFR    = int(NFR4,8)
                 NELIM4 = 0
                 II = TEMP(I)
 845             NELIM4 = NELIM4 + 1
                 II = FILS(II)
                 IF (II .GT. 0 ) GOTO 845
                 NELIM=int(NELIM4,8)
                 CALL MUMPS_GET_FLOPS_COST(NFR4,NELIM4,NELIM4,
     &                SYM,1,COST_NODE)
                 TAB1(I)=int(dble(COST_NODE)+
     &                COST_TRAV(STEP(INODE)),8)
                 TAB2(I)=0_8
              ELSE
                 SON(I)=IN
              ENDIF
           ELSE
              SON(I)=IN
           ENDIF
           IN=FRERE(STEP(IN))
        ENDDO
        IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
           IF((SLAVEF.NE.1).AND.(.NOT.MUMPS_IN_OR_ROOT_SSARBR(
     &          PROCNODE(STEP(INODE)),K199)))THEN
              CALL DMUMPS_FUSION_SORT(TEMP,NE(STEP(INODE)),TAB1,TAB2,
     &             LOCAL_PERM
     &             ,RESULT,T1,T2)
              TAB=>TEMP
              DO I=NE(STEP(INODE)),1,-1
                 IF(I.EQ.NE(STEP(INODE))) THEN
                    FILS(dernier)=-TAB(I)
                    dernier=TAB(I)
                    GOTO 221
                 ENDIF
                 IF(I.EQ.1) THEN
                    FRERE(STEP(dernier))=TAB(I)
                    FRERE(STEP(TAB(I)))=-INODE
                    GOTO 221
                 ENDIF
                 IF(I.GT.1) THEN
                    FRERE(STEP(dernier))=TAB(I)
                    dernier=TAB(I)
                    GOTO 221
                 ENDIF
 221             CONTINUE
                 SON(NE(STEP(INODE))-I+1)=TAB(I)
              ENDDO
           ENDIF
        ENDIF
        DO I=1,NE(STEP(INODE))
           IPOOL(fin)=SON(I)
           SON(I)=0
           fin=fin+1
        ENDDO
        IF(NE(STEP(INODE)).EQ.0)THEN
           IF(PERM.NE.7)THEN
              NA(LEAF+2)=INODE
           ENDIF
           LEAF=LEAF-1
        ELSE
           fin=fin-1
           GOTO 999
        ENDIF
        fin=fin-1
        IF(fin.EQ.0) THEN
           GOTO 789
        ENDIF
        GOTO 999
 789    CONTINUE
        IF(PERM.EQ.7) GOTO 5483
        NBROOT=NA(2)
        NBLEAF=NA(1)
        PEAK=0.0D0
        FACT_SIZE=0_8
        DO I=1,NBROOT
           PEAK=max(PEAK,dble(M(STEP(NA(2+NBLEAF+I)))))
           FACT_SIZE=FACT_SIZE+fact(STEP(NA(2+NBLEAF+I)))
        ENDDO
 5483   CONTINUE
        DEALLOCATE(IPOOL)
        DEALLOCATE(fact)
        DEALLOCATE(TNSTK)
        DEALLOCATE(SON)
        DEALLOCATE(TAB2)
        DEALLOCATE(TAB1)
        DEALLOCATE(T1)
        DEALLOCATE(T2)
        DEALLOCATE(RESULT)
        DEALLOCATE(TEMP)
        IF(PERM.NE.7)THEN
           DEALLOCATE(M)
        ENDIF
        IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
           DEALLOCATE(DEPTH)
        ENDIF
        IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
           DEALLOCATE(COST_TRAV)
        ENDIF
        IF ((PERM.NE.7).AND.(SBTR_M.OR.(PERM.EQ.2))) THEN
           IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1).OR.(PERM.EQ.2))THEN
              DEALLOCATE(M_TOTAL)
           ENDIF
        ENDIF
      RETURN
      END SUBROUTINE DMUMPS_REORDER_TREE
      SUBROUTINE DMUMPS_BUILD_LOAD_MEM_INFO(N,FRERE, STEP, FILS,
     &     NA,LNA,NE,ND, DAD, LDAD, USE_DAD,
     &     NSTEPS,PERM,SYM,INFO,LP,K47,K81,K76,K215,K234,K55,KEEP199,
     &     PROCNODE,MEM_SUBTREE,SLAVEF, SIZE_MEM_SBTR, PEAK
     &     ,SBTR_WHICH_M,SIZE_DEPTH_FIRST,SIZE_COST_TRAV,
     &     DEPTH_FIRST_TRAV,DEPTH_FIRST_SEQ,COST_TRAV,MY_FIRST_LEAF,
     &     MY_NB_LEAF,MY_ROOT_SBTR,SBTR_ID
     &     )
      IMPLICIT NONE
      INTEGER N,PERM,SYM, NSTEPS, LNA, LP, SIZE_MEM_SBTR,LDAD
      INTEGER FRERE(NSTEPS), FILS(N), STEP(N)
      INTEGER NA(LNA), NE(NSTEPS), ND(NSTEPS)
      INTEGER K47,K81,K76,K215,K234,K55,KEEP199
      INTEGER DAD(LDAD)
      LOGICAL USE_DAD
      INTEGER INFO(80)
      INTEGER SLAVEF,PROCNODE(NSTEPS)
      DOUBLE PRECISION, intent(out) :: MEM_SUBTREE(SIZE_MEM_SBTR,SLAVEF)
      INTEGER :: SBTR_WHICH_M
      INTEGER MY_FIRST_LEAF(SIZE_MEM_SBTR,SLAVEF),
     &     MY_ROOT_SBTR(SIZE_MEM_SBTR,SLAVEF),
     &     MY_NB_LEAF(SIZE_MEM_SBTR,SLAVEF)
      EXTERNAL MUMPS_ROOTSSARBR,MUMPS_PROCNODE
      LOGICAL MUMPS_ROOTSSARBR
      INTEGER MUMPS_PROCNODE
      DOUBLE PRECISION PEAK
      INTEGER SIZE_DEPTH_FIRST,DEPTH_FIRST_TRAV(SIZE_DEPTH_FIRST),
     &     DEPTH_FIRST_SEQ(SIZE_DEPTH_FIRST)
      INTEGER SIZE_COST_TRAV
      INTEGER SBTR_ID(SIZE_DEPTH_FIRST),OOC_CUR_SBTR
      DOUBLE PRECISION COST_TRAV(SIZE_COST_TRAV)
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH
      INTEGER IFATH,IN,INODE,I,allocok,LOCAL_PERM
      INTEGER(8) NELIM,NFR
      INTEGER NFR4,NELIM4
      INTEGER LEAF,NBLEAF,NBROOT, SIZE_TAB
      INTEGER, DIMENSION (:), ALLOCATABLE :: IPOOL,TNSTK
      INTEGER, DIMENSION (:), ALLOCATABLE,TARGET :: SON,TEMP
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: M,M_TOTAL, fact
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: TAB1,TAB2
      INTEGER x,dernier,fin,RANK_TRAV
      INTEGER II
      INTEGER ROOT_OF_CUR_SBTR
      INTEGER(8), DIMENSION (:), ALLOCATABLE :: T1,T2
      INTEGER, DIMENSION (:), ALLOCATABLE :: RESULT
      INTEGER(8) MEM_SIZE,FACT_SIZE,
     &     TOTAL_MEM_SIZE,
     &     SIZECB
      LOGICAL   SBTR_M
      INTEGER,DIMENSION(:),ALLOCATABLE :: INDICE
      INTEGER ID,FIRST_LEAF,SIZE_SBTR
      EXTERNAL MUMPS_IN_OR_ROOT_SSARBR,MUMPS_INSSARBR
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR,MUMPS_INSSARBR
      DOUBLE PRECISION COST_NODE
      INTEGER CUR_DEPTH_FIRST_RANK
      INCLUDE 'mumps_headers.h'
      TOTAL_MEM_SIZE=0_8
      ROOT_OF_CUR_SBTR=0
      ALLOCATE(INDICE( SLAVEF ), stat=allocok)
      IF (allocok > 0) THEN
         IF ( LP .GT. 0 )
     &        WRITE(LP,*)'Memory allocation error in
     &DMUMPS_REORDER_TREE'
         INFO(1)=-7
         INFO(2)=SLAVEF
         RETURN
      ENDIF
      IF((PERM.EQ.0).OR.(PERM.EQ.1).OR.
     &     (PERM.EQ.2).OR.(PERM.EQ.3).OR.(PERM.EQ.4).OR.
     &     (PERM.EQ.5).OR.(PERM.EQ.6))THEN
         LOCAL_PERM=0
      ENDIF
      IF (K47 == 4 .OR. ((K47.GE.2).AND.(K81.GE. 1))) THEN
        DO I=1,SLAVEF
          INDICE(I)=1
        ENDDO
        DO I=1,SLAVEF
          DO x=1,SIZE_MEM_SBTR
            MEM_SUBTREE(x,I)=-1.0D0
          ENDDO
        ENDDO
      ENDIF
      SBTR_M=((K47 == 4 .OR. ((K47.GE.2).AND.(K81 .GE. 1))))
      MEM_SIZE=0_8
      FACT_SIZE=0_8
      IF ((PERM.GT.7).AND.
     &     (.NOT.(K47 == 4 .OR. ((K47.GE.2).AND.(K81 .GE. 1))))) THEN
         WRITE(*,*) "Internal Error in DMUMPS_REORDER_TREE",PERM
         CALL MUMPS_ABORT()
      END IF
      NBLEAF = NA(1)
      NBROOT = NA(2)
      CUR_DEPTH_FIRST_RANK=1
      IF((PERM.EQ.0).AND.(NBROOT.EQ.NBLEAF)) THEN
        DEALLOCATE(INDICE)
        RETURN
      ENDIF
      IF (SBTR_M.OR.(PERM.EQ.2))  THEN
         IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
            ALLOCATE(M_TOTAL(NSTEPS), stat=allocok )
            IF (allocok > 0) THEN
               IF ( LP .GT. 0 )
     &              WRITE(LP,*)'Memory allocation error in
     &              DMUMPS_REORDER_TREE'
               INFO(1)=-7
               INFO(2)=NSTEPS
               RETURN
            ENDIF
         ENDIF
      ENDIF
      ALLOCATE( IPOOL(NBLEAF), M(NSTEPS), fact(NSTEPS),
     &          TNSTK(NSTEPS), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_REORDER_TREE'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      II=0
      DO I=1,NSTEPS
         TNSTK(I) = NE(I)
         IF(NE(I).GE.II) II=NE(I)
      ENDDO
      SIZE_TAB=max(II,NBROOT)
      ALLOCATE(SON(II), TEMP(II),
     &         TAB1(SIZE_TAB), TAB2(SIZE_TAB), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_REORDER_TREE'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      ALLOCATE(T1(SIZE_TAB),T2(SIZE_TAB),
     &         RESULT(SIZE_TAB),stat=allocok)
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in DMUMPS_REORDER_TREE'
        INFO(1)=-7
        INFO(2)=SIZE_TAB
        RETURN
      ENDIF
      IF(NBROOT.EQ.NBLEAF)THEN
         IF((PERM.NE.1).OR.(PERM.EQ.4).OR.(PERM.EQ.6))THEN
          WRITE(*,*)'Internal Error in reordertree:'
          WRITE(*,*)'  problem with perm parameter in reordertree'
          CALL MUMPS_ABORT()
        ENDIF
        DO I=1,NBROOT
          TAB1(I)=int(ND(STEP(NA(I+2+NBLEAF))),8)
          IPOOL(I)=NA(I+2+NBLEAF)
          M(STEP(IPOOL(I)))=TAB1(I)*TAB1(I)
        ENDDO
        CALL DMUMPS_FUSION_SORT(NA(2+NBLEAF+1),NBROOT,TAB1,TAB2,4,
     &    RESULT,T1,T2)
        GOTO 789
      ENDIF
      IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
         ALLOCATE(DEPTH(NSTEPS),stat=allocok)
         IF (allocok > 0) THEN
            IF ( LP .GT. 0 )
     &           WRITE(LP,*)'Memory allocation error in
     &           DMUMPS_REORDER_TREE'
            INFO(1)=-7
            INFO(2)=NSTEPS
            RETURN
         ENDIF
         DEPTH=0
         NBROOT = NA(2)
         IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
         fin=NBROOT
         LEAF=NA(1)
 499     CONTINUE
         INODE=IPOOL(fin)
         IF(INODE.LT.0)THEN
            WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
            CALL MUMPS_ABORT()
         ENDIF
         IN=INODE
 4602    IN = FILS(IN)
         IF (IN .GT. 0 ) THEN
            GOTO 4602
         ENDIF
         IN=-IN
         DO I=1,NE(STEP(INODE))
            SON(I)=IN
            IN=FRERE(STEP(IN))
         ENDDO
         DO I=1,NE(STEP(INODE))
            IPOOL(fin)=SON(I)
            DEPTH(STEP(SON(I)))=DEPTH(STEP(INODE))+1
            SON(I)=0
            fin=fin+1
         ENDDO
         IF(NE(STEP(INODE)).EQ.0)THEN
            LEAF=LEAF-1
         ELSE
            fin=fin-1
            GOTO 499
         ENDIF
         fin=fin-1
         IF(fin.EQ.0) GOTO 489
         GOTO 499
 489     CONTINUE
      ENDIF
      IF(K76.EQ.4.OR.(K76.EQ.6))THEN
         RANK_TRAV=NSTEPS
         DEPTH_FIRST_TRAV=0
         DEPTH_FIRST_SEQ=0
      ENDIF
      IF((K76.EQ.5).OR.(PERM.EQ.5).OR.(PERM.EQ.6))THEN
         COST_TRAV=0.0D0
         COST_NODE=0.0d0
      ENDIF        
      DO I=1,NSTEPS
         M(I)=0_8
         IF (SBTR_M.OR.(PERM.EQ.2))  THEN
            IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
               M_TOTAL(I)=0_8
            ENDIF
         ENDIF
      ENDDO
      DO I=1,NSTEPS
         fact(I)=0_8
      ENDDO
        NBROOT = NA(2)
        NBLEAF = NA(1)
        IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
        CONTINUE
        fin=NBROOT
        LEAF=NA(1)
        FIRST_LEAF=-9999
        SIZE_SBTR=0
 999    CONTINUE
        INODE=IPOOL(fin)
        IF(INODE.LT.0)THEN
           WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
           CALL MUMPS_ABORT()
        ENDIF
        IF(SIZE_SBTR.NE.0)THEN 
           IF(.NOT.MUMPS_INSSARBR(PROCNODE(STEP(INODE)),KEEP199))THEN
              IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
                 IF((SLAVEF.NE.1))THEN
                    MY_FIRST_LEAF(INDICE(ID+1)-1,ID+1)=FIRST_LEAF
                    MY_NB_LEAF(INDICE(ID+1)-1,ID+1)=SIZE_SBTR
                    FIRST_LEAF=-9999
                    SIZE_SBTR=0
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
        IF(MUMPS_ROOTSSARBR(PROCNODE(STEP(INODE)),KEEP199))THEN
           ROOT_OF_CUR_SBTR=INODE
        ENDIF
        IF (K76.EQ.4)THEN
           IF(SLAVEF.NE.1)THEN
              WRITE(*,*)'INODE=',INODE,'RANK',RANK_TRAV
              IF(MUMPS_INSSARBR(PROCNODE(STEP(INODE)),KEEP199))THEN
                 DEPTH_FIRST_TRAV(STEP(INODE))=DEPTH_FIRST_TRAV(STEP(
     &                ROOT_OF_CUR_SBTR))
              ELSE
                 DEPTH_FIRST_TRAV(STEP(INODE))=RANK_TRAV
              ENDIF
              RANK_TRAV=RANK_TRAV-1
           ENDIF
        ENDIF
        IF (K76.EQ.5)THEN
           IF(SLAVEF.NE.1)THEN
              IF (USE_DAD) THEN
                IFATH=DAD(INODE)
              ELSE
                IN = INODE
 395            IN = FRERE(IN)
                IF (IN.GT.0) GO TO 395
                IFATH = -IN
              ENDIF
              NFR4   = ND(STEP(INODE))
              NFR    = int(NFR4,8)
              NELIM4 = 0
              IN = INODE
 396          NELIM4 = NELIM4 + 1
              IN = FILS(IN)
              IF (IN .GT. 0 ) GOTO 396
              NELIM=int(NELIM4,8)
              IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                 SIZECB=(NFR-NELIM)*(NFR-NELIM)
              ELSE
                 SIZECB=(NFR-NELIM)*(NFR-NELIM+1_8)/2_8
              ENDIF
              CALL MUMPS_GET_FLOPS_COST(NFR4,NELIM4,NELIM4,
     &             SYM,1,COST_NODE)
              IF(IFATH.NE.0)THEN
                 IF(MUMPS_INSSARBR(PROCNODE(STEP(INODE)),KEEP199))THEN
                    COST_TRAV(STEP(INODE))=COST_TRAV(STEP(
     &                   ROOT_OF_CUR_SBTR))
                 ELSE
                    COST_TRAV(STEP(INODE))=dble(COST_NODE)+
     &                   COST_TRAV(STEP(IFATH))+
     &                   dble(SIZECB*18_8)  
                 ENDIF
              ELSE
                 COST_TRAV(STEP(INODE))=dble(COST_NODE)
              ENDIF
              IF(K76.EQ.5)THEN
                 WRITE(*,*)'INODE=',INODE,'COST=',COST_TRAV(STEP(INODE))
              ENDIF
           ENDIF
        ENDIF
        IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
              IF((SLAVEF.NE.1).AND.
     &          MUMPS_ROOTSSARBR(PROCNODE(STEP(INODE)),KEEP199))THEN
                IF (NE(STEP(INODE)).NE.0) THEN
                   ID=MUMPS_PROCNODE(PROCNODE(STEP(INODE)),KEEP199)
                   IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
                      MEM_SUBTREE(INDICE(ID+1),ID+1)=
     &                     dble(M_TOTAL(STEP(INODE)))
                   ELSE
                      MEM_SUBTREE(INDICE(ID+1),ID+1)=
     &                     dble(M(STEP(INODE)))
                   ENDIF
                   MY_ROOT_SBTR(INDICE(ID+1),ID+1)=INODE
                  INDICE(ID+1)=INDICE(ID+1)+1
                ENDIF
              ENDIF
              IF((SLAVEF.EQ.1).AND.FRERE(STEP(INODE)).EQ.0)THEN
                 ID=MUMPS_PROCNODE(PROCNODE(STEP(INODE)),KEEP199)
                 IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
                    MEM_SUBTREE(INDICE(ID+1),ID+1)=
     &                   dble(M_TOTAL(STEP(INODE)))
                 ELSE
                    MEM_SUBTREE(INDICE(ID+1),ID+1)=
     &                   dble(M(STEP(INODE)))
                 ENDIF
                 INDICE(ID+1)=INDICE(ID+1)+1
              ENDIF
        ENDIF
        IN=INODE
 5602   IN = FILS(IN)
        IF (IN .GT. 0 ) THEN
           dernier=IN
           GOTO 5602
        ENDIF
        IN=-IN
        DO I=1,NE(STEP(INODE))
           IPOOL(fin)=IN
           IF(IN.GT.0) IN=FRERE(STEP(IN))
           fin=fin+1
        ENDDO
        IF(NE(STEP(INODE)).EQ.0)THEN
           IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
              IF(SLAVEF.NE.1)THEN
                 IF(MUMPS_INSSARBR(PROCNODE(STEP(INODE)),KEEP199))THEN
                    IF(FIRST_LEAF.EQ.-9999)THEN
                       FIRST_LEAF=INODE
                    ENDIF
                    SIZE_SBTR=SIZE_SBTR+1
                 ENDIF
              ENDIF
           ENDIF
           IF(PERM.NE.7)THEN
              NA(LEAF+2)=INODE
           ENDIF
           LEAF=LEAF-1
        ELSE
           fin=fin-1
           GOTO 999
        ENDIF
        fin=fin-1
        IF(fin.EQ.0) THEN
           IF(SIZE_SBTR.NE.0)THEN 
              IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
                 IF((SLAVEF.NE.1))THEN
                    MY_FIRST_LEAF(INDICE(ID+1)-1,ID+1)=FIRST_LEAF
                    MY_NB_LEAF(INDICE(ID+1)-1,ID+1)=SIZE_SBTR
                    FIRST_LEAF=-9999
                    SIZE_SBTR=0
                 ENDIF
              ENDIF
           ENDIF
           GOTO 789
        ENDIF
        GOTO 999
 789    CONTINUE
        IF(K76.EQ.6)THEN
           OOC_CUR_SBTR=1
           DO I=1,NSTEPS
              TNSTK(I) = NE(I)
           ENDDO           
           NBROOT=NA(2)
           NBLEAF=NA(1)
           IPOOL(1:NBLEAF)=NA(3:2+NBLEAF)
           LEAF = NBLEAF + 1
 9100      CONTINUE
           IF (LEAF.NE.1) THEN
              LEAF = LEAF -1
              INODE = IPOOL(LEAF)
           ENDIF
 9600      CONTINUE
           IF(SLAVEF.NE.1)THEN
              ID=MUMPS_PROCNODE(PROCNODE(STEP(INODE)),KEEP199)
              DEPTH_FIRST_TRAV(STEP(INODE))=CUR_DEPTH_FIRST_RANK
              DEPTH_FIRST_SEQ(CUR_DEPTH_FIRST_RANK)=INODE
              CUR_DEPTH_FIRST_RANK=CUR_DEPTH_FIRST_RANK+1
              IF(MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &             KEEP199))THEN
                 SBTR_ID(STEP(INODE))=OOC_CUR_SBTR
              ELSE
                 SBTR_ID(STEP(INODE))=-9999
              ENDIF
              IF(MUMPS_ROOTSSARBR(PROCNODE(STEP(INODE)),
     &             KEEP199))THEN
                 OOC_CUR_SBTR=OOC_CUR_SBTR+1
              ENDIF
           ENDIF
           IF (USE_DAD) THEN
              IFATH = DAD( STEP(INODE) )
           ELSE
              IN = INODE
 1133         IN = FRERE(IN)
              IF (IN.GT.0) GO TO 1133
              IFATH = -IN
           ENDIF
           IF (IFATH.EQ.0) THEN
              NBROOT = NBROOT - 1
              IF (NBROOT.EQ.0) GOTO 1163
              GOTO 9100
           ENDIF
           TNSTK(STEP(IFATH))=TNSTK(STEP(IFATH))-1
           IF(TNSTK(STEP(IFATH)).EQ.0) THEN
              INODE=IFATH
              GOTO 9600
           ELSE
              GOTO 9100
           ENDIF        
 1163      CONTINUE
        ENDIF
        PEAK=0.0D0
        FACT_SIZE=0_8
        DO I=1,NBROOT
           PEAK=max(PEAK,dble(M(STEP(NA(2+NBLEAF+I)))))
           FACT_SIZE=FACT_SIZE+fact(STEP(NA(2+NBLEAF+I)))
        ENDDO
        CONTINUE
        DEALLOCATE(IPOOL)
        DEALLOCATE(M)
        DEALLOCATE(fact)
        DEALLOCATE(TNSTK)
        DEALLOCATE(SON)
        DEALLOCATE(TAB2)
        DEALLOCATE(TAB1)
        DEALLOCATE(T1)
        DEALLOCATE(T2)
        DEALLOCATE(RESULT)
        DEALLOCATE(TEMP)
        DEALLOCATE(INDICE)
        IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
           DEALLOCATE(DEPTH)
        ENDIF
        IF (SBTR_M.OR.(PERM.EQ.2))  THEN
           IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1).OR.(PERM.EQ.2))THEN
              DEALLOCATE(M_TOTAL)
           ENDIF
        ENDIF
      RETURN
      END SUBROUTINE DMUMPS_BUILD_LOAD_MEM_INFO
      RECURSIVE SUBROUTINE DMUMPS_FUSION_SORT(TAB,DIM,TAB1,TAB2,PERM,
     &  RESULT,TEMP1,TEMP2)
      IMPLICIT NONE
      INTEGER DIM
      INTEGER(8) TAB1(DIM),TAB2(DIM)
      INTEGER(8) TEMP1(DIM),TEMP2(DIM)
      INTEGER TAB(DIM), PERM,RESULT(DIM)
      INTEGER I,J,I1,I2
      IF(DIM.EQ.1) THEN
        RESULT(1)=TAB(1)
        TEMP1(1)=TAB1(1)
        TEMP2(1)=TAB2(1)
        RETURN
      ENDIF
      I=DIM/2
      CALL DMUMPS_FUSION_SORT(TAB(1),I,TAB1(1),TAB2(1),PERM,
     &  RESULT(1),TEMP1(1),TEMP2(1))
      CALL DMUMPS_FUSION_SORT(TAB(I+1),DIM-I,TAB1(I+1),TAB2(I+1),
     &  PERM,RESULT(I+1),TEMP1(I+1),TEMP2(I+1))
      I1=1
      I2=I+1
      J=1
      DO WHILE ((I1.LE.I).AND.(I2.LE.DIM))
        IF((PERM.EQ.3))THEN
          IF(TEMP1(I1).LE.TEMP1(I2))THEN
            TAB(J)=RESULT(I1)
            TAB1(J)=TEMP1(I1)
            J=J+1
            I1=I1+1
          ELSE
            TAB(J)=RESULT(I2)
            TAB1(J)=TEMP1(I2)
            J=J+1
            I2=I2+1
          ENDIF
          GOTO 3
        ENDIF
        IF((PERM.EQ.4).OR.(PERM.EQ.5))THEN
          IF (TEMP1(I1).GE.TEMP1(I2))THEN
            TAB(J)=RESULT(I1)
            TAB1(J)=TEMP1(I1)
            J=J+1
            I1=I1+1
          ELSE
            TAB(J)=RESULT(I2)
            TAB1(J)=TEMP1(I2)
            J=J+1
            I2=I2+1          
          ENDIF
          GOTO 3
        ENDIF
        IF((PERM.EQ.0).OR.(PERM.EQ.1).OR.(PERM.EQ.2)) THEN
          IF(TEMP1(I1).GT.TEMP1(I2))THEN
            TAB1(J)=TEMP1(I1)
            TAB2(J)=TEMP2(I1)
            TAB(J)=RESULT(I1)
            J=J+1
            I1=I1+1
            GOTO 3
          ENDIF
          IF(TEMP1(I1).LT.TEMP1(I2))THEN
            TAB1(J)=TEMP1(I2)
            TAB2(J)=TEMP2(I2)
            TAB(J)=RESULT(I2)
            J=J+1
            I2=I2+1
            GOTO 3
          ENDIF        
          IF((TEMP1(I1).EQ.TEMP1(I2)))THEN
            IF(TEMP2(I1).LE.TEMP2(I2))THEN
              TAB1(J)=TEMP1(I1)
              TAB2(J)=TEMP2(I1)
              TAB(J)=RESULT(I1)
              J=J+1
              I1=I1+1
            ELSE
              TAB1(J)=TEMP1(I2)
              TAB2(J)=TEMP2(I2)
              TAB(J)=RESULT(I2)
              J=J+1
              I2=I2+1
            ENDIF
          ENDIF
        ENDIF
  3   CONTINUE    
      ENDDO
      IF(I1.GT.I)THEN
        DO WHILE(I2.LE.DIM)
          TAB(J)=RESULT(I2)
          TAB1(J)=TEMP1(I2)
          TAB2(J)=TEMP2(I2)
          J=J+1
          I2=I2+1
        ENDDO
      ELSE
        IF(I2.GT.DIM)THEN
          DO WHILE(I1.LE.I)
            TAB1(J)=TEMP1(I1)
            TAB2(J)=TEMP2(I1)
            TAB(J)=RESULT(I1)
            J=J+1
            I1=I1+1
          ENDDO
        ENDIF
      ENDIF
      DO I=1,DIM
        TEMP1(I)=TAB1(I)
        TEMP2(I)=TAB2(I)
        RESULT(I)=TAB(I)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_FUSION_SORT
