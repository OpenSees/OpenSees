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
      SUBROUTINE DMUMPS_INIT_POOL_LAST3(IPOOL, LPOOL, LEAF)
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER LPOOL, LEAF
      INTEGER IPOOL(LPOOL)
      IPOOL(LPOOL-2) = 0
      IPOOL(LPOOL-1) = 0
      IPOOL(LPOOL)   = LEAF-1
      RETURN
      END SUBROUTINE DMUMPS_INIT_POOL_LAST3
      SUBROUTINE DMUMPS_INSERT_POOL_N
     &           (N, POOL, LPOOL, PROCNODE, SLAVEF, KEEP199,
     &           K28, K76, K80, K47, STEP, INODE)
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER N, INODE, LPOOL, K28, SLAVEF, K76, K80, K47, KEEP199
      INTEGER STEP(N), POOL(LPOOL), PROCNODE(K28)
      EXTERNAL MUMPS_IN_OR_ROOT_SSARBR
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR, ATM_CURRENT_NODE
      INTEGER NBINSUBTREE, NBTOP, INODE_EFF,POS_TO_INSERT
      INTEGER IPOS1, IPOS2, ISWAP
      INTEGER NODE,J,I
      ATM_CURRENT_NODE = ( K76 == 2 .OR. K76 ==3 .OR.
     &     K76==4 .OR. K76==5)
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      IF (INODE > N ) THEN
        INODE_EFF = INODE - N
      ELSE IF (INODE < 0) THEN
        INODE_EFF = - INODE
      ELSE
        INODE_EFF = INODE
      ENDIF
      IF(((INODE.GT.0).AND.(INODE.LE.N)).AND.(.NOT.
     &     MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE_EFF)),
     &               KEEP199))
     &  ) THEN
         IF ((K80 == 1 .AND. K47 .GE. 1) .OR.
     &     (( K80 == 2 .OR. K80==3 ) .AND.
     &          ( K47 == 4 ))) THEN
            CALL DMUMPS_REMOVE_NODE(INODE,1)
         ENDIF
      ENDIF
      IF ( MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE_EFF)),
     &                             KEEP199) ) THEN
        POOL(NBINSUBTREE + 1 ) = INODE
        NBINSUBTREE = NBINSUBTREE + 1
      ELSE
         POS_TO_INSERT=NBTOP+1
         IF((K76.EQ.4).OR.(K76.EQ.5).OR.(K76.EQ.6))THEN
            IF((INODE.GT.N).OR.(INODE.LE.0))THEN
               DO J=NBTOP,1,-1
                  IF((POOL(LPOOL-2-J).GT.0)
     &                 .AND.(POOL(LPOOL-2-J).LE.N))THEN
                     GOTO 333
                  ENDIF
                  IF ( POOL(LPOOL-2-J) < 0 ) THEN
                     NODE=-POOL(LPOOL-2-J)
                  ELSE IF ( POOL(LPOOL-2-J) > N ) THEN
                     NODE = POOL(LPOOL-2-J) - N
                  ELSE
                     NODE = POOL(LPOOL-2-J)
                  ENDIF
                  IF((K76.EQ.4).OR.(K76.EQ.6))THEN
                     IF(DEPTH_FIRST_LOAD(STEP(NODE)).GE.
     &                    DEPTH_FIRST_LOAD(STEP(INODE_EFF)))THEN
                        GOTO 333
                     ENDIF
                  ENDIF
                  IF(K76.EQ.5)THEN
                     IF(COST_TRAV(STEP(NODE)).LE.
     &                    COST_TRAV(STEP(INODE_EFF)))THEN
                        GOTO 333
                     ENDIF
                  ENDIF
                  POS_TO_INSERT=POS_TO_INSERT-1
               ENDDO
               IF(J.EQ.0) J=1
 333           CONTINUE
               DO I=NBTOP,POS_TO_INSERT,-1
                  POOL(LPOOL-2-I-1)=POOL(LPOOL-2-I)
               ENDDO
               POOL(LPOOL-2-POS_TO_INSERT)=INODE
               NBTOP = NBTOP + 1
               GOTO 20 
            ENDIF
            DO J=NBTOP,1,-1
               IF((POOL(LPOOL-2-J).GT.0).AND.(POOL(LPOOL-2-J).LE.N))THEN
                  GOTO 888
               ENDIF
               POS_TO_INSERT=POS_TO_INSERT-1
            ENDDO
 888        CONTINUE
            DO I=J,1,-1
               NODE=POOL(LPOOL-2-I)
               IF((K76.EQ.4).OR.(K76.EQ.6))THEN
                  IF(DEPTH_FIRST_LOAD(STEP(NODE)).GE.
     &                 DEPTH_FIRST_LOAD(STEP(INODE_EFF)))THEN
                     GOTO 999
                  ENDIF
               ENDIF
               IF(K76.EQ.5)THEN
                  IF(COST_TRAV(STEP(NODE)).LE.
     &                 COST_TRAV(STEP(INODE_EFF)))THEN
                     GOTO 999
                  ENDIF
               ENDIF
               POS_TO_INSERT=POS_TO_INSERT-1
            ENDDO
            IF(I.EQ.0) I=1
 999        CONTINUE
            DO J=NBTOP,POS_TO_INSERT,-1
               POOL(LPOOL-2-J-1)=POOL(LPOOL-2-J)
            ENDDO
            POOL(LPOOL-2-POS_TO_INSERT)=INODE
            NBTOP = NBTOP + 1
            GOTO 20
         ENDIF
         POOL( LPOOL - 2 - ( NBTOP + 1 ) ) = INODE
         NBTOP = NBTOP + 1
        IPOS1 = LPOOL - 2 - NBTOP
        IPOS2 = LPOOL - 2 - NBTOP + 1
 10     CONTINUE
        IF ( IPOS2 == LPOOL - 2 ) GOTO 20
        IF ( POOL(IPOS1) < 0 ) GOTO 20
        IF ( POOL(IPOS2) < 0 ) GOTO 30
        IF ( ATM_CURRENT_NODE ) THEN
          IF ( POOL(IPOS1) > N ) GOTO 20
          IF ( POOL(IPOS2) > N ) GOTO 30
        END IF
        GOTO 20
 30     CONTINUE
        ISWAP = POOL(IPOS1)
        POOL(IPOS1) = POOL(IPOS2)
        POOL(IPOS2) = ISWAP
        IPOS1 = IPOS1 + 1
        IPOS2 = IPOS2 + 1
        GOTO 10
 20     CONTINUE
      ENDIF
      POOL(LPOOL) = NBINSUBTREE 
      POOL(LPOOL - 1) = NBTOP
      RETURN
      END SUBROUTINE DMUMPS_INSERT_POOL_N
      LOGICAL FUNCTION DMUMPS_POOL_EMPTY(POOL, LPOOL)
      IMPLICIT NONE
      INTEGER LPOOL
      INTEGER POOL(LPOOL)
      INTEGER NBINSUBTREE, NBTOP
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      DMUMPS_POOL_EMPTY = (NBINSUBTREE + NBTOP == 0)
      RETURN
      END FUNCTION DMUMPS_POOL_EMPTY
      SUBROUTINE DMUMPS_EXTRACT_POOL( N, POOL, LPOOL, PROCNODE, SLAVEF,
     &           STEP, INODE, KEEP,KEEP8, MYID, ND,
     &           FORCE_EXTRACT_TOP_SBTR )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER INODE, LPOOL, SLAVEF, N
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N), POOL(LPOOL), PROCNODE(KEEP(28)),
     &        ND(KEEP(28))
      EXTERNAL MUMPS_INSSARBR, MUMPS_ROOTSSARBR, DMUMPS_POOL_EMPTY
      LOGICAL MUMPS_INSSARBR, MUMPS_ROOTSSARBR, DMUMPS_POOL_EMPTY
      EXTERNAL MUMPS_PROCNODE
      INTEGER MUMPS_PROCNODE
      INTEGER NBINSUBTREE, NBTOP, INSUBTREE, INODE_EFF, MYID
      LOGICAL LEFT, ATOMIC_SUBTREE,UPPER,FLAG_MEM,SBTR_FLAG,PROC_FLAG
      LOGICAL FORCE_EXTRACT_TOP_SBTR
      INTEGER NODE_TO_EXTRACT,I,J,MIN_PROC
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF ( KEEP(76) > 6 .OR. KEEP(76) < 0 ) THEN
         WRITE(*,*) "Error 2 in DMUMPS_EXTRACT_POOL: unknown strategy"
         CALL MUMPS_ABORT()
      ENDIF
      ATOMIC_SUBTREE =  ( KEEP(76) == 1 .OR. KEEP(76) == 3)
      IF ( DMUMPS_POOL_EMPTY(POOL, LPOOL) ) THEN
         WRITE(*,*) "Error 1 in DMUMPS_EXTRACT_POOL"
         CALL MUMPS_ABORT()
      ENDIF
      IF ( .NOT. ATOMIC_SUBTREE ) THEN
         LEFT = (NBTOP == 0)
         IF(.NOT.LEFT)THEN
            IF((KEEP(76).EQ.4).OR.(KEEP(76).EQ.5))THEN
               IF(NBINSUBTREE.EQ.0)THEN
                  LEFT=.FALSE.
               ELSE
                  IF ( POOL(NBINSUBTREE) < 0 ) THEN
                     I = -POOL(NBINSUBTREE)
                  ELSE IF ( POOL(NBINSUBTREE) > N ) THEN
                     I = POOL(NBINSUBTREE) - N
                  ELSE
                     I = POOL(NBINSUBTREE)
                  ENDIF
                  IF ( POOL(LPOOL-2-NBTOP) < 0 ) THEN
                     J = -POOL(LPOOL-2-NBTOP)
                  ELSE IF ( POOL(LPOOL-2-NBTOP) > N ) THEN
                     J = POOL(LPOOL-2-NBTOP) - N
                  ELSE
                     J = POOL(LPOOL-2-NBTOP)
                  ENDIF
                  IF(KEEP(76).EQ.4)THEN
                     IF(DEPTH_FIRST_LOAD(STEP(J)).GE.
     &                    DEPTH_FIRST_LOAD(STEP(I)))THEN
                        LEFT=.TRUE.
                     ELSE
                        LEFT=.FALSE.
                     ENDIF
                  ENDIF
                  IF(KEEP(76).EQ.5)THEN
                     IF(COST_TRAV(STEP(J)).LE.
     &                    COST_TRAV(STEP(I)))THEN
                        LEFT=.TRUE.
                     ELSE
                        LEFT=.FALSE.
                     ENDIF
                  ENDIF
               ENDIF           
            ENDIF
         ENDIF
      ELSE
         IF ( INSUBTREE == 1 ) THEN
            IF (NBINSUBTREE == 0) THEN
               WRITE(*,*) "Error 3 in DMUMPS_EXTRACT_POOL"
               CALL MUMPS_ABORT()
            ENDIF
            LEFT = .TRUE.
         ELSE
            LEFT = ( NBTOP == 0)
         ENDIF
      ENDIF
 222  CONTINUE
      IF ( LEFT ) THEN
         INODE = POOL( NBINSUBTREE )
         IF(KEEP(81).EQ.2)THEN
            IF((INODE.GE.0).AND.(INODE.LE.N))THEN
               CALL DMUMPS_MEM_NODE_SELECT(INODE,POOL,LPOOL,N,
     &              STEP,KEEP,KEEP8,PROCNODE,SLAVEF,MYID,SBTR_FLAG,
     &              PROC_FLAG,MIN_PROC)
               IF(.NOT.SBTR_FLAG)THEN
                  WRITE(*,*)MYID,': ca a change pour moi'
                  LEFT=.FALSE.
                  GOTO 222
               ENDIF
            ENDIF
         ELSEIF(KEEP(81).EQ.3)THEN
            IF((INODE.GE.0).AND.(INODE.LE.N))THEN
               NODE_TO_EXTRACT=INODE
               FLAG_MEM=.FALSE.
               CALL DMUMPS_LOAD_CHK_MEMCST_POOL(FLAG_MEM)
               IF(FLAG_MEM)THEN
                  CALL DMUMPS_MEM_NODE_SELECT(INODE,POOL,LPOOL,N,
     &                 STEP,KEEP,KEEP8,
     &                 PROCNODE,SLAVEF,MYID,SBTR_FLAG,
     &                 PROC_FLAG,MIN_PROC)
                  IF(.NOT.SBTR_FLAG)THEN
                     LEFT=.FALSE.
                     WRITE(*,*)MYID,': ca a change pour moi (2)'
                     GOTO 222
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         NBINSUBTREE = NBINSUBTREE - 1
         IF ( INODE < 0 ) THEN
            INODE_EFF = -INODE
         ELSE IF ( INODE > N ) THEN
            INODE_EFF = INODE - N
         ELSE
            INODE_EFF = INODE
         ENDIF
         IF ( MUMPS_INSSARBR( PROCNODE(STEP(INODE_EFF)),
     &                        KEEP(199)) ) THEN
            IF((KEEP(47).GE.2.AND.KEEP(81).EQ.1).AND.
     &           (INSUBTREE.EQ.0))THEN
               CALL DMUMPS_LOAD_SET_SBTR_MEM(.TRUE.)
            ENDIF
            INSUBTREE = 1
         ELSE IF ( MUMPS_ROOTSSARBR( PROCNODE(STEP(INODE_EFF)), 
     &           KEEP(199))) THEN
            IF((KEEP(47).GE.2.AND.KEEP(81).EQ.1).AND.
     &           (INSUBTREE.EQ.1))THEN
               CALL DMUMPS_LOAD_SET_SBTR_MEM(.FALSE.)
            ENDIF
            INSUBTREE = 0
         END IF
      ELSE
         IF (NBTOP < 1 ) THEN
            WRITE(*,*) "Error 5 in DMUMPS_EXTRACT_POOL", NBTOP
            CALL MUMPS_ABORT()
         ENDIF
         INODE = POOL( LPOOL - 2 - NBTOP )
         IF(KEEP(81).EQ.1)THEN
            CALL DMUMPS_LOAD_POOL_CHECK_MEM
     &           (INODE,UPPER,SLAVEF,KEEP,KEEP8,
     &            STEP,POOL,LPOOL,PROCNODE,N)
            IF(UPPER)THEN
               GOTO 666
            ELSE
               NBINSUBTREE=NBINSUBTREE-1
               IF ( MUMPS_INSSARBR( PROCNODE(STEP(INODE)), 
     &              KEEP(199)) ) THEN
                  INSUBTREE = 1
               ELSE IF ( MUMPS_ROOTSSARBR( PROCNODE(STEP(INODE)), 
     &                 KEEP(199))) THEN
                  INSUBTREE = 0
               ENDIF
               GOTO 777
            ENDIF
         ENDIF
         IF(KEEP(81).EQ.2)THEN
            CALL DMUMPS_MEM_NODE_SELECT(INODE,POOL,LPOOL,N,STEP,
     &           KEEP,KEEP8,
     &           PROCNODE,SLAVEF,MYID,SBTR_FLAG,PROC_FLAG,MIN_PROC)
            IF(SBTR_FLAG)THEN
               LEFT=.TRUE. 
               WRITE(*,*)MYID,': ca a change pour moi (3)'              
               GOTO 222
            ENDIF
         ELSE
            IF(KEEP(81).EQ.3)THEN
               IF((INODE.GE.0).AND.(INODE.LE.N))THEN
                  NODE_TO_EXTRACT=INODE
                  FLAG_MEM=.FALSE.
                  CALL DMUMPS_LOAD_CHK_MEMCST_POOL(FLAG_MEM)
                  IF(FLAG_MEM)THEN
                     CALL DMUMPS_MEM_NODE_SELECT(INODE,POOL,LPOOL,N,
     &                    STEP,KEEP,KEEP8,
     &                    PROCNODE,SLAVEF,MYID,SBTR_FLAG,
     &                    PROC_FLAG,MIN_PROC)
                     IF(SBTR_FLAG)THEN
                        LEFT=.TRUE.
                        WRITE(*,*)MYID,': ca a change pour moi (4)'
                        GOTO 222
                     ENDIF
                  ELSE
                     CALL DMUMPS_LOAD_CLEAN_MEMINFO_POOL(INODE)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
 666     CONTINUE
         NBTOP = NBTOP - 1
         IF((INODE.GT.0).AND.(INODE.LE.N))THEN
            IF ((( KEEP(80) == 2 .OR. KEEP(80)==3 ) .AND.
     &           ( KEEP(47) == 4 ))) THEN
               CALL DMUMPS_REMOVE_NODE(INODE,2)
            ENDIF
         ENDIF
         IF ( INODE < 0 ) THEN
            INODE_EFF = -INODE
         ELSE IF ( INODE > N ) THEN
            INODE_EFF = INODE - N
         ELSE
            INODE_EFF = INODE
         ENDIF
      END IF
 777  CONTINUE
      POOL(LPOOL)     = NBINSUBTREE 
      POOL(LPOOL - 1) = NBTOP
      POOL(LPOOL - 2) = INSUBTREE
      RETURN
      END SUBROUTINE DMUMPS_EXTRACT_POOL
      SUBROUTINE DMUMPS_MEM_CONS_MNG(INODE,POOL,LPOOL,N,STEP,
     &     KEEP,KEEP8,
     &     PROCNODE,SLAVEF,MYID,SBTR,FLAG_SAME_PROC,MIN_PROC)
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER INODE,LPOOL,N,MYID,SLAVEF,PROC,MIN_PROC
      INTEGER POOL(LPOOL),KEEP(500),STEP(N),PROCNODE(KEEP(28))
      INTEGER(8) KEEP8(150)
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      LOGICAL SBTR,FLAG_SAME_PROC
      INTEGER POS_TO_EXTRACT,NODE_TO_EXTRACT,NBTOP,I,INSUBTREE,
     &     NBINSUBTREE
      DOUBLE PRECISION MIN_COST, TMP_COST
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      MIN_COST=huge(MIN_COST) 
      TMP_COST=huge(TMP_COST)
      FLAG_SAME_PROC=.FALSE.
      SBTR=.FALSE.
      MIN_PROC=-9999
      IF((INODE.GT.0).AND.(INODE.LE.N))THEN
         POS_TO_EXTRACT=-1
         NODE_TO_EXTRACT=-1
         DO I=NBTOP,1,-1
            IF(NODE_TO_EXTRACT.LT.0)THEN
               POS_TO_EXTRACT=I
               NODE_TO_EXTRACT=POOL(LPOOL-2-I)
               CALL DMUMPS_LOAD_COMP_MAXMEM_POOL(NODE_TO_EXTRACT,
     &                                       TMP_COST,PROC)
               MIN_COST=TMP_COST
               MIN_PROC=PROC
            ELSE
               CALL DMUMPS_LOAD_COMP_MAXMEM_POOL(POOL(LPOOL-2-I),
     &                                       TMP_COST,PROC)
               IF((PROC.NE.MIN_PROC).OR.(TMP_COST.NE.MIN_COST))THEN
                  FLAG_SAME_PROC=.TRUE.
               ENDIF
               IF(TMP_COST.GT.MIN_COST)THEN
                  POS_TO_EXTRACT=I
                  NODE_TO_EXTRACT=POOL(LPOOL-2-I)
                  MIN_COST=TMP_COST
                  MIN_PROC=PROC
               ENDIF
            ENDIF
         ENDDO
         IF((KEEP(47).EQ.4).AND.(NBINSUBTREE.NE.0))THEN
            CALL DMUMPS_CHECK_SBTR_COST(NBINSUBTREE,INSUBTREE,NBTOP,
     &           MIN_COST,SBTR)
            IF(SBTR)THEN
               WRITE(*,*)MYID,': selecting from subtree'
               RETURN
            ENDIF
         ENDIF
         IF((.NOT.SBTR).AND.(.NOT.FLAG_SAME_PROC))THEN
            WRITE(*,*)MYID,': I must search for a task
     &           to save My friend'
            RETURN
         ENDIF
         INODE = NODE_TO_EXTRACT
         DO I=POS_TO_EXTRACT,NBTOP
            IF(I.NE.NBTOP)THEN
               POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
            ENDIF
         ENDDO
         POOL(LPOOL-2-NBTOP)=INODE
         CALL DMUMPS_LOAD_CLEAN_MEMINFO_POOL(INODE)
      ELSE
      ENDIF
      END SUBROUTINE DMUMPS_MEM_CONS_MNG
      SUBROUTINE DMUMPS_MEM_NODE_SELECT(INODE,POOL,LPOOL,N,STEP,
     &     KEEP,KEEP8,
     &     PROCNODE,SLAVEF,MYID,SBTR_FLAG,PROC_FLAG,MIN_PROC)
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER INODE,LPOOL,N,SLAVEF,MYID,MIN_PROC
      INTEGER POOL(LPOOL),KEEP(500),PROCNODE(KEEP(28)),STEP(N)
      INTEGER(8) KEEP8(150)
      LOGICAL SBTR_FLAG,PROC_FLAG
      EXTERNAL MUMPS_INSSARBR
      LOGICAL MUMPS_INSSARBR
      INTEGER NODE_TO_EXTRACT,I,POS_TO_EXTRACT,NBTOP,NBINSUBTREE
      NBTOP= POOL(LPOOL - 1)
      NBINSUBTREE = POOL(LPOOL)
      IF(NBTOP.GT.0)THEN
         WRITE(*,*)MYID,': NBTOP=',NBTOP
      ENDIF
      SBTR_FLAG=.FALSE.
      PROC_FLAG=.FALSE.
      CALL DMUMPS_MEM_CONS_MNG(INODE,POOL,LPOOL,N,STEP,KEEP,KEEP8,
     &     PROCNODE,SLAVEF,MYID,SBTR_FLAG,PROC_FLAG,MIN_PROC)
      IF(SBTR_FLAG)THEN
         RETURN
      ENDIF
      IF(MIN_PROC.EQ.-9999)THEN
         IF((INODE.GT.0).AND.(INODE.LT.N))THEN
            SBTR_FLAG=(NBINSUBTREE.NE.0)
         ENDIF
         RETURN
      ENDIF
      IF(.NOT.PROC_FLAG)THEN
         NODE_TO_EXTRACT=INODE
         IF((INODE.GE.0).AND.(INODE.LE.N))THEN
            CALL DMUMPS_FIND_BEST_NODE_FOR_MEM(MIN_PROC,POOL,
     &           LPOOL,INODE)
            IF(MUMPS_INSSARBR(PROCNODE(STEP(INODE)),
     &           KEEP(199)))THEN
               WRITE(*,*)MYID,': Extracting from a subtree
     &              for helping',MIN_PROC
               SBTR_FLAG=.TRUE.
               RETURN
            ELSE
               IF(NODE_TO_EXTRACT.NE.INODE)THEN
                  WRITE(*,*)MYID,': Extracting from top
     &                 inode=',INODE,'for helping',MIN_PROC
               ENDIF
               CALL DMUMPS_LOAD_CLEAN_MEMINFO_POOL(INODE)
            ENDIF
         ENDIF
         DO I=1,NBTOP
            IF (POOL(LPOOL-2-I).EQ.INODE)THEN
               GOTO 452
            ENDIF
         ENDDO
 452     CONTINUE
         POS_TO_EXTRACT=I
         DO I=POS_TO_EXTRACT,NBTOP-1
            POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
         ENDDO
         POOL(LPOOL-2-NBTOP)=INODE
      ENDIF
      END SUBROUTINE DMUMPS_MEM_NODE_SELECT
      SUBROUTINE DMUMPS_GET_INODE_FROM_POOL
     &           ( IPOOL, LPOOL, III, LEAF, 
     &             INODE, STRATEGIE )
            IMPLICIT NONE
      INTEGER, INTENT(IN) :: STRATEGIE, LPOOL
      INTEGER IPOOL (LPOOL)
      INTEGER III,LEAF
      INTEGER, INTENT(OUT) :: INODE
         LEAF  = LEAF - 1
         INODE = IPOOL( LEAF )
      RETURN
      END SUBROUTINE DMUMPS_GET_INODE_FROM_POOL
