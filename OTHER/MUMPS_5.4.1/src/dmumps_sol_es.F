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
      MODULE DMUMPS_SOL_ES
      PRIVATE
      PUBLIC:: PRUNED_SIZE_LOADED
      PUBLIC:: DMUMPS_CHAIN_PRUN_NODES
      PUBLIC:: DMUMPS_CHAIN_PRUN_NODES_STATS
      PUBLIC:: DMUMPS_INITIALIZE_RHS_BOUNDS
      PUBLIC:: DMUMPS_PROPAGATE_RHS_BOUNDS
      PUBLIC:: DMUMPS_TREE_PRUN_NODES
      PUBLIC:: DMUMPS_TREE_PRUN_NODES_STATS
      PUBLIC:: DMUMPS_SOL_ES_INIT
      INTEGER(8), POINTER, DIMENSION(:,:) :: SIZE_OF_BLOCK
      INTEGER(8) :: PRUNED_SIZE_LOADED
      INCLUDE 'mumps_headers.h'
      CONTAINS
      SUBROUTINE DMUMPS_SOL_ES_INIT(SIZE_OF_BLOCK_ARG, KEEP201)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KEEP201
      INTEGER(8), POINTER, DIMENSION(:,:) :: SIZE_OF_BLOCK_ARG
      IF (KEEP201 > 0) THEN
        SIZE_OF_BLOCK => SIZE_OF_BLOCK_ARG
      ELSE
        NULLIFY(SIZE_OF_BLOCK)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SOL_ES_INIT
      SUBROUTINE DMUMPS_TREE_PRUN_NODES( 
     &     fill,
     &     DAD, NE_STEPS, FRERE, KEEP28,
     &     FILS, STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     TO_PROCESS,
     &     nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &     Pruned_List, Pruned_Roots, Pruned_Leaves
     &     )
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: fill
      INTEGER, INTENT(IN) :: N, KEEP28
      INTEGER, INTENT(IN) :: DAD(KEEP28),NE_STEPS(KEEP28),FRERE(KEEP28)
      INTEGER, INTENT(IN) :: FILS(N), STEP(N)
      INTEGER, INTENT(IN) :: nodes_RHS(KEEP28),  nb_nodes_RHS
      INTEGER :: nb_prun_nodes
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_List(nb_prun_nodes)
      INTEGER :: nb_prun_roots
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Roots(nb_prun_roots)
      INTEGER :: nb_prun_leaves
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Leaves(nb_prun_leaves)
      LOGICAL :: TO_PROCESS(KEEP28) 
      INTEGER :: IN, I, ISTEP, TMP, TMPsave
      LOGICAL :: FILS_VISITED
      nb_prun_nodes = 0
      nb_prun_leaves = 0
      TO_PROCESS(:) = .FALSE.
      DO I = 1, nb_nodes_RHS
         TMP = nodes_RHS(I)
         TMPsave = TMP
         ISTEP = STEP(TMP)
         DO WHILE(.NOT.TO_PROCESS(ISTEP))
            TO_PROCESS(ISTEP) = .TRUE.
            nb_prun_nodes = nb_prun_nodes + 1
            IF(fill) THEN
               Pruned_List(nb_prun_nodes) = TMP
            END IF
            IN = FILS(TMP) 
            DO WHILE(IN.GT.0) 
               IN = FILS(IN)
            END DO
            FILS_VISITED = .FALSE.
            IF (IN.LT.0) THEN 
             FILS_VISITED = TO_PROCESS(STEP(-IN))
            ENDIF
            IF ( IN.LT.0.and..NOT.FILS_VISITED)
     &            THEN 
               TMP = -IN
               ISTEP = STEP(TMP)
            ELSE 
               IF (IN.EQ.0) THEN
                 nb_prun_leaves = nb_prun_leaves + 1
                 IF (fill) THEN
                    Pruned_Leaves(nb_prun_leaves) = TMP
                 END IF
               ELSE 
                 TMP = -IN
                 ISTEP = STEP(TMP)
               ENDIF
               DO WHILE (TMP.NE.TMPsave) 
                  TMP = abs(FRERE(ISTEP))
                  IF(TMP.NE.0) THEN 
                     ISTEP = STEP(TMP) 
                  ELSE 
                     exit
                  END IF
                  IF (.NOT.TO_PROCESS(ISTEP)) exit
               END DO
            END IF
         END DO
      END DO
      nb_prun_roots = 0
      DO I=1,nb_nodes_RHS
         TMP = nodes_RHS(I)
         ISTEP = STEP(TMP)
         IF(DAD(ISTEP).NE.0) THEN 
            IF(.NOT.TO_PROCESS(STEP(DAD(ISTEP)))) THEN
               nb_prun_roots = nb_prun_roots + 1
               IF(fill) THEN
                  Pruned_Roots(nb_prun_roots) = TMP
               END IF
            END IF
         ELSE 
            nb_prun_roots = nb_prun_roots + 1
            IF(fill) THEN
               Pruned_Roots(nb_prun_roots) = TMP
            END IF          
         END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_TREE_PRUN_NODES
      SUBROUTINE DMUMPS_CHAIN_PRUN_NODES(
     &     fill,
     &     DAD, KEEP28,
     &     STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     Pruned_SONS, TO_PROCESS,
     &     nb_prun_nodes,nb_prun_roots, nb_prun_leaves,
     &     Pruned_List, Pruned_Roots, Pruned_Leaves
     &     )
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: fill
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: STEP(N)
      INTEGER, INTENT(IN) :: KEEP28
      INTEGER, INTENT(IN) :: DAD(KEEP28)
      INTEGER, INTENT(IN) :: nb_nodes_RHS
      INTEGER, INTENT(IN) :: nodes_RHS(nb_nodes_RHS)
      INTEGER :: nb_prun_nodes
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_List(nb_prun_nodes)
      INTEGER :: nb_prun_roots
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Roots(nb_prun_roots)
      INTEGER :: nb_prun_leaves
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Leaves(nb_prun_leaves)
      INTEGER :: Pruned_SONS(KEEP28)
      LOGICAL :: TO_PROCESS(KEEP28)
      INTEGER :: IN, I, ISTEP, TMP
      nb_prun_nodes = 0
      nb_prun_roots = 0
      TO_PROCESS(:) = .FALSE.
      Pruned_SONS(:) = -1
      DO I = 1, nb_nodes_RHS
         TMP = nodes_RHS(I)
         ISTEP = STEP(TMP)
         TO_PROCESS(ISTEP) = .TRUE.
         IF (Pruned_SONS(ISTEP) .eq. -1) THEN
            Pruned_SONS(ISTEP) = 0
            nb_prun_nodes = nb_prun_nodes + 1
            IF(fill) THEN
               Pruned_List(nb_prun_nodes) = nodes_RHS(I)
            END IF
            IN = nodes_RHS(I)
            IN = DAD(STEP(IN))
            DO WHILE (IN.NE.0)
               TO_PROCESS(STEP(IN)) = .TRUE.
               IF (Pruned_SONS(STEP(IN)).eq.-1) THEN 
                  nb_prun_nodes = nb_prun_nodes + 1
                  IF(fill) THEN
                     Pruned_List(nb_prun_nodes) = IN
                  END IF
                  Pruned_SONS(STEP(IN)) = 1
                  TMP = IN
                  IN = DAD(STEP(IN))
               ELSE 
                  Pruned_SONS(STEP(IN)) = Pruned_SONS(STEP(IN)) + 1
                  GOTO 201
               ENDIF
            ENDDO
            nb_prun_roots = nb_prun_roots +1
            IF(fill) THEN
               Pruned_Roots(nb_prun_roots) = TMP
            END IF
         ENDIF
  201    CONTINUE
      ENDDO
      nb_prun_leaves = 0
      DO I = 1, nb_nodes_RHS
         TMP = nodes_RHS(I)
         ISTEP = STEP(TMP)
         IF (Pruned_SONS(ISTEP).EQ.0) THEN
            nb_prun_leaves = nb_prun_leaves +1
            IF(fill) THEN
              Pruned_Leaves(nb_prun_leaves) = TMP
            END IF
         END IF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_CHAIN_PRUN_NODES
      SUBROUTINE DMUMPS_INITIALIZE_RHS_BOUNDS(
     & STEP, N,
     & IRHS_PTR, NBCOL, IRHS_SPARSE, NZ_RHS,
     & JBEG_RHS, PERM_RHS, SIZE_PERM_RHS, K242, K243,
     & UNS_PERM_INV, SIZE_UNS_PERM_INV, K23,
     & RHS_BOUNDS, NSTEPS,
     & nb_sparse, MYID,
     & mode)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MYID, N, NSTEPS, K242, K243, K23
      INTEGER, INTENT(IN) :: JBEG_RHS, SIZE_PERM_RHS, nb_sparse
      INTEGER, INTENT(IN) :: NBCOL, NZ_RHS, SIZE_UNS_PERM_INV
      INTEGER, INTENT(IN) :: STEP(N), PERM_RHS(SIZE_PERM_RHS)
      INTEGER, INTENT(IN) :: IRHS_PTR(NBCOL+1),IRHS_SPARSE(NZ_RHS)
      INTEGER, INTENT(IN) :: UNS_PERM_INV(SIZE_UNS_PERM_INV)
      INTEGER, INTENT(INOUT):: RHS_BOUNDS(2*NSTEPS)
      INTEGER, INTENT(IN) :: mode 
      INTEGER :: I, ICOL, JPTR, J, JAM1, node, bound
      RHS_BOUNDS = 0
      ICOL = 0 
      DO I = 1, NBCOL
        IF ( (IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE
        ICOL = ICOL + 1
        bound = ICOL - mod(ICOL, nb_sparse) + 1
        IF(mod(ICOL, nb_sparse).EQ.0) bound = bound - nb_sparse
        IF(mode.EQ.0) THEN 
          IF ((K242.NE.0).OR.(K243.NE.0)) THEN
            JAM1 = PERM_RHS(JBEG_RHS+I-1)
          ELSE
            JAM1 = JBEG_RHS+I-1
          ENDIF
          node = abs(STEP(JAM1))
          IF(RHS_BOUNDS(2*node - 1).EQ.0) THEN 
            RHS_BOUNDS(2*node - 1) = bound                 
            RHS_BOUNDS(2*node)     = bound + nb_sparse - 1 
          ELSE
            RHS_BOUNDS(2*node) = bound + nb_sparse - 1
          END IF
        ELSE  
          DO JPTR = IRHS_PTR(I), IRHS_PTR(I+1)-1
            J = IRHS_SPARSE(JPTR)
            IF ( mode .EQ. 1 ) THEN
              IF (K23.NE.0) J = UNS_PERM_INV(J)
            ENDIF
            node = abs(STEP(J))
            IF(RHS_BOUNDS(2*node - 1).EQ.0) THEN
              RHS_BOUNDS(2*node - 1) = bound
              RHS_BOUNDS(2*node)     = bound + nb_sparse - 1
            ELSE
              RHS_BOUNDS(2*node) = bound + nb_sparse - 1
            END IF
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_INITIALIZE_RHS_BOUNDS
      SUBROUTINE DMUMPS_PROPAGATE_RHS_BOUNDS(
     & pruned_leaves, nb_pruned_leaves,
     & STEP, N, Pruned_SONS,
     & DAD, RHS_BOUNDS, NSTEPS,
     & MYID, COMM, KEEP485,
     & IW, LIW, PTRIST, KIXSZ,OOC_FCT_LOC, PHASE, LDLT, K38)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_headers.h'
      INTEGER, INTENT(IN) :: nb_pruned_leaves, N, NSTEPS
      INTEGER, INTENT(IN) :: STEP(N), DAD(NSTEPS), Pruned_SONS(NSTEPS)
      INTEGER, INTENT(IN) :: MYID, COMM, KEEP485
      INTEGER, INTENT(IN) :: pruned_leaves(nb_pruned_leaves)
      INTEGER, INTENT(IN) :: LIW, IW(LIW), PTRIST(NSTEPS)
      INTEGER, INTENT(IN) :: KIXSZ, OOC_FCT_LOC, PHASE, LDLT, K38
      INTEGER, INTENT(INOUT):: RHS_BOUNDS(2*NSTEPS)
      INTEGER :: I, node, father, size_pool, next_size_pool
      INTEGER :: IERR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: POOL, NBSONS
      ALLOCATE(POOL(nb_pruned_leaves),
     &         NBSONS(NSTEPS),
     &         STAT=IERR)
      IF (IERR.NE.0) THEN
         WRITE(6,*)'Allocation problem in DMUMPS_PROPAGATE_RHS_BOUNDS'
         CALL MUMPS_ABORT()
      END IF
      size_pool = nb_pruned_leaves
      POOL = pruned_leaves
      NBSONS = Pruned_SONS
      DO WHILE (size_pool.ne.0)
        next_size_pool =0 
        DO I=1, size_pool 
          node = STEP(POOL(I)) 
          IF (DAD(node).NE.0) THEN
            father = STEP(DAD(node))
            NBSONS(father) = NBSONS(father)-1
            IF (RHS_BOUNDS(2*father-1).EQ.0) THEN
              RHS_BOUNDS(2*father-1) = RHS_BOUNDS(2*node-1)
              RHS_BOUNDS(2*father)   = RHS_BOUNDS(2*node)
            ELSE
              RHS_BOUNDS(2*father-1) = min(RHS_BOUNDS(2*father-1),
     &                                     RHS_BOUNDS(2*node-1))
              RHS_BOUNDS(2*father) = max(RHS_BOUNDS(2*father),
     &                                     RHS_BOUNDS(2*node))
            END IF
            IF(NBSONS(father).EQ.0) THEN 
              next_size_pool = next_size_pool+1
              POOL(next_size_pool) = DAD(node)
            END IF
          END IF
        END DO
        size_pool = next_size_pool 
      END DO
      DEALLOCATE(POOL, NBSONS)
      RETURN
      END SUBROUTINE DMUMPS_PROPAGATE_RHS_BOUNDS
      INTEGER(8) FUNCTION DMUMPS_LOCAL_FACTOR_SIZE(IW,LIW,PTR,
     &                                 PHASE, LDLT, IS_ROOT)
        INTEGER, INTENT(IN) :: LIW, PTR, PHASE, LDLT
        INTEGER, INTENT(IN) :: IW(LIW)
        LOGICAL, INTENT(IN) :: IS_ROOT
        INTEGER(8) :: NCB, NELIM, LIELL, NPIV, NROW
        NCB   = int(IW(PTR),8)     
        NELIM = int(IW(PTR+1),8)   
        NROW  = int(IW(PTR+2),8)
        NPIV  = int(IW(PTR+3),8)
        LIELL = NPIV + NCB
        IF (IS_ROOT) THEN
          DMUMPS_LOCAL_FACTOR_SIZE = int(IW(PTR+1),8) *     
     &                               int(IW(PTR+2),8) / 2_8 
          RETURN
        ENDIF
        IF (NCB.GE.0_8) THEN 
          IF (PHASE.EQ.0   
     &      .OR. (PHASE.EQ.1.AND.LDLT.NE.0) 
     &  ) THEN
            DMUMPS_LOCAL_FACTOR_SIZE =
     &            NPIV*(NPIV-1_8)/2_8 + (NROW-NPIV)*NPIV
          ELSE
            DMUMPS_LOCAL_FACTOR_SIZE =
     &      NPIV*(NPIV+1_8)/2_8 + (LIELL-NPIV)*NPIV
          ENDIF
        ELSE
          DMUMPS_LOCAL_FACTOR_SIZE =
     &      -NCB*NELIM
        END IF
      RETURN
      END FUNCTION DMUMPS_LOCAL_FACTOR_SIZE
      INTEGER(8) FUNCTION DMUMPS_LOCAL_FACTOR_SIZE_BLR(IW,LIW,PTR,
     &                                LRSTATUS, IWHANDLER,
     &                                PHASE, LDLT, IS_ROOT)
      USE DMUMPS_LR_DATA_M
      USE DMUMPS_LR_TYPE
        INTEGER, INTENT(IN) :: LIW, PTR, PHASE, LDLT
        INTEGER, INTENT(IN) :: LRSTATUS, IWHANDLER
        INTEGER, INTENT(IN) :: IW(LIW)
        LOGICAL, INTENT(IN) :: IS_ROOT
        INTEGER(8) :: NCB, NELIM, LIELL, NPIV, NROW, FACTOR_SIZE
        INTEGER :: NB_PANELS, IPANEL, LorU, IBLOCK
        LOGICAL :: LR_ACTIVATED
        TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: LRB_PANEL
        NCB   = int(IW(PTR),8)     
        NELIM = int(IW(PTR+1),8)   
        NROW  = int(IW(PTR+2),8)
        NPIV  = int(IW(PTR+3),8)
        LIELL = NPIV + NCB
        LR_ACTIVATED=(LRSTATUS.GE.2)
        IF (LR_ACTIVATED) THEN
          FACTOR_SIZE = 0_8
          CALL DMUMPS_BLR_RETRIEVE_NB_PANELS(IWHANDLER, NB_PANELS)
          IF (LDLT.EQ.0) THEN 
            LorU = PHASE
          ELSE
            LorU = 0
          ENDIF
          DO IPANEL=1,NB_PANELS
            IF (IS_ROOT.AND.IPANEL.EQ.NB_PANELS) THEN
              CYCLE
            ENDIF
            IF (DMUMPS_BLR_EMPTY_PANEL_LORU(IWHANDLER, LorU, IPANEL))
     &        THEN
              CYCLE
            ENDIF
             CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(IWHANDLER, LorU,
     &                       IPANEL, LRB_PANEL)
            IF (size(LRB_PANEL).GT.0) THEN
              IF (PHASE.EQ.0) THEN
                FACTOR_SIZE = FACTOR_SIZE + 
     &            int(LRB_PANEL(1)%N,8)*(int(LRB_PANEL(1)%N,8)-1_8)/2_8
              ELSE
                FACTOR_SIZE = FACTOR_SIZE + 
     &           int(LRB_PANEL(1)%N,8)*(int(LRB_PANEL(1)%N,8)+1_8)/2_8
              ENDIF
            ENDIF
            DO IBLOCK=1,size(LRB_PANEL)
              IF (LRB_PANEL(IBLOCK)%ISLR) THEN
                FACTOR_SIZE = FACTOR_SIZE + int(LRB_PANEL(IBLOCK)%K,8)*
     &            int(LRB_PANEL(IBLOCK)%M+LRB_PANEL(IBLOCK)%M,8)
              ELSE
                FACTOR_SIZE = FACTOR_SIZE +
     &             int(LRB_PANEL(IBLOCK)%M*LRB_PANEL(IBLOCK)%N,8)
              ENDIF
            ENDDO
          ENDDO
          DMUMPS_LOCAL_FACTOR_SIZE_BLR = FACTOR_SIZE
        ELSE
          DMUMPS_LOCAL_FACTOR_SIZE_BLR = 
     &      DMUMPS_LOCAL_FACTOR_SIZE(IW, LIW, PTR, PHASE, LDLT, IS_ROOT)
        ENDIF
      RETURN
      END FUNCTION DMUMPS_LOCAL_FACTOR_SIZE_BLR
      SUBROUTINE DMUMPS_TREE_PRUN_NODES_STATS(MYID, N, KEEP28, KEEP201,
     &           FR_FACT,
     &           STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_LOC)
      INTEGER, intent(in) :: KEEP28, KEEP201, OOC_FCT_TYPE_LOC, MYID, N
      INTEGER(8), intent(in) :: FR_FACT
      INTEGER, intent(in) :: nb_prun_nodes
      INTEGER, intent(in) :: Pruned_List(nb_prun_nodes)
      INTEGER, intent(in) :: STEP(N)
      INTEGER I, ISTEP
      INTEGER(8) :: Pruned_Size
      IF (KEEP201 .GT. 0) THEN
        Pruned_Size = 0_8
        DO I = 1, nb_prun_nodes
          ISTEP = STEP(Pruned_List(I))
          Pruned_Size = Pruned_Size + SIZE_OF_BLOCK
     &                  (ISTEP, OOC_FCT_TYPE_LOC)
        ENDDO
        PRUNED_SIZE_LOADED = PRUNED_SIZE_LOADED +Pruned_Size
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_TREE_PRUN_NODES_STATS
      SUBROUTINE DMUMPS_CHAIN_PRUN_NODES_STATS
     &                (MYID, N, KEEP28, KEEP201, KEEP485, FR_FACT,
     &                STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_LOC
     & )
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP28, KEEP201, OOC_FCT_TYPE_LOC, N,
     &                       KEEP485
      INTEGER(8), intent(in) :: FR_FACT
      INTEGER, intent(in) :: nb_prun_nodes, MYID
      INTEGER, intent(in) :: Pruned_List(nb_prun_nodes)
      INTEGER, intent(in) :: STEP(N)
      INCLUDE 'mpif.h'
      INTEGER I, ISTEP
      INTEGER(8) :: Pruned_Size
      Pruned_Size = 0_8
      DO I = 1, nb_prun_nodes
        ISTEP = STEP(Pruned_List(I))
        IF (KEEP201 .GT. 0) THEN
            Pruned_Size = Pruned_Size + SIZE_OF_BLOCK
     &                    (ISTEP, OOC_FCT_TYPE_LOC)
        ENDIF
      ENDDO
      IF (KEEP201.GT.0) THEN
        IF (FR_FACT .NE. 0_8) THEN
          PRUNED_SIZE_LOADED = PRUNED_SIZE_LOADED +Pruned_Size
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_CHAIN_PRUN_NODES_STATS
      END MODULE DMUMPS_SOL_ES
      SUBROUTINE DMUMPS_PERMUTE_RHS_GS
     &          (LP, LPOK, PROKG, MPG, PERM_STRAT, 
     &           SYM_PERM, N, NRHS,
     &           IRHS_PTR, SIZE_IRHS_PTR, 
     &           IRHS_SPARSE, NZRHS, 
     &           PERM_RHS, IERR
     &         )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LP, MPG, PERM_STRAT, N, NRHS, 
     &                       SIZE_IRHS_PTR,
     &                       NZRHS
      LOGICAL, INTENT(IN) :: LPOK, PROKG
      INTEGER, INTENT(IN) :: SYM_PERM(N)
      INTEGER, INTENT(IN) :: IRHS_PTR(SIZE_IRHS_PTR)
      INTEGER, INTENT(IN) :: IRHS_SPARSE(NZRHS)
      INTEGER, INTENT(OUT) :: PERM_RHS(NRHS)
      INTEGER, INTENT(OUT) :: IERR
      INTEGER :: I,J,K, POSINPERMRHS, JJ,
     &           KPOS
      INTEGER, ALLOCATABLE :: ROW_REFINDEX(:)
      IERR = 0
      IF ((PERM_STRAT.NE.-1).AND.(PERM_STRAT.NE.1)) THEN
       IERR=-1
       IF (LPOK)
     & WRITE(LP,*) " INTERNAL ERROR -1 in ",
     &       " DMUMPS_PERMUTE_RHS_GS, PERM_STRAT =", PERM_STRAT, 
     &       " is out of range "
       RETURN
      ENDIF
      IF (PERM_STRAT.EQ.-1) THEN
       DO I=1,NRHS
        PERM_RHS(I) = I
       END DO
       GOTO 490
      ENDIF
      ALLOCATE(ROW_REFINDEX(NRHS), STAT=IERR)
      IF (IERR.GT.0) THEN
       IERR=-1
       IF (LPOK) THEN
          WRITE(LP,*) " ERROR -2 : ", 
     &         " ALLOCATE IN DMUMPS_PERMUTE_RHS_GS OF SIZE :",
     &         NRHS
       ENDIF
       RETURN
      ENDIF
      DO I=1,NRHS
        IF (IRHS_PTR(I+1)-IRHS_PTR(I).LE.0) THEN
          IERR =  1
          IF (I.EQ.1) THEN
            ROW_REFINDEX(I) = IRHS_SPARSE(IRHS_PTR(I))
          ELSE
            ROW_REFINDEX(I) = ROW_REFINDEX(I-1)
          ENDIF
        ELSE
          ROW_REFINDEX(I) = IRHS_SPARSE(IRHS_PTR(I))
        ENDIF
      END DO
      POSINPERMRHS = 0
      DO I=1,NRHS
       KPOS = N+1 
       JJ   = 0   
       DO J=1,NRHS
        K = ROW_REFINDEX(J)
        IF (K.LE.0) CYCLE 
        IF (SYM_PERM(K).LT.KPOS) THEN
         KPOS = SYM_PERM(K)
         JJ   = J
        ENDIF
       END DO
       IF (JJ.EQ.0) THEN
         IERR = -3 
         IF (LPOK)
     &   WRITE(LP,*) " INTERNAL ERROR -3 in ",
     &       " DMUMPS_PERMUTE_RHS_GS "
         GOTO 500
       ENDIF
       POSINPERMRHS           = POSINPERMRHS + 1
       PERM_RHS(POSINPERMRHS) = JJ
       ROW_REFINDEX(JJ)       = -ROW_REFINDEX(JJ)
      END DO
      IF (POSINPERMRHS.NE.NRHS) THEN
         IF (LPOK)
     &   WRITE(LP,*) " INTERNAL ERROR -4 in ",
     &       " DMUMPS_PERMUTE_RHS_GS ", maxval(ROW_REFINDEX)
         IERR = -4
         GOTO 500
      ENDIF
  490 CONTINUE
 500  CONTINUE
      IF (allocated(ROW_REFINDEX)) DEALLOCATE(ROW_REFINDEX)
      END SUBROUTINE DMUMPS_PERMUTE_RHS_GS
      SUBROUTINE DMUMPS_PERMUTE_RHS_AM1
     &          (PERM_STRAT, SYM_PERM,
     &           IRHS_PTR, NHRS,
     &           PERM_RHS, SIZEPERM, IERR
     &         )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: PERM_STRAT, NHRS, SIZEPERM
      INTEGER, INTENT(IN) :: SYM_PERM(SIZEPERM)
      INTEGER, INTENT(IN) :: IRHS_PTR(NHRS)
      INTEGER, INTENT(OUT):: IERR   
      INTEGER, INTENT(OUT):: PERM_RHS(SIZEPERM)
      DOUBLE PRECISION :: RAND_NUM
      INTEGER  I, J, STRAT
      IERR = 0
      STRAT = PERM_STRAT 
      IF( (STRAT.NE.-3).AND.
     &    (STRAT.NE.-2).AND.
     &    (STRAT.NE.-1).AND.
     &    (STRAT.NE. 1).AND.
     &    (STRAT.NE. 2).AND.
     &    (STRAT.NE. 6) ) THEN
        WRITE(*,*)"Warning: incorrect value for the RHS permutation; ",
     &            "defaulting to post-order"
        STRAT = 1
      END IF
      IF (STRAT .EQ. -3) THEN
         PERM_RHS(1:SIZEPERM)=0
         DO I=1, SIZEPERM 
           CALL random_number(RAND_NUM) 
           RAND_NUM = RAND_NUM*dble(SIZEPERM) 
           J = ceiling(RAND_NUM) 
           DO WHILE (PERM_RHS(J).NE.0) 
             CALL random_number(RAND_NUM)
             RAND_NUM = RAND_NUM*dble(SIZEPERM)
             J = ceiling(RAND_NUM)
           ENDDO
           PERM_RHS(J)=I
         ENDDO
      ELSEIF (STRAT .EQ. -2) THEN
         DO I=1, SIZEPERM
            PERM_RHS(SIZEPERM -I +1) = I
         ENDDO
      ELSEIF (STRAT .EQ. -1) THEN
         DO I=1, SIZEPERM
            PERM_RHS(I) = I
         ENDDO
      ELSEIF (STRAT .EQ.  1) THEN
         DO I=1, SIZEPERM
            PERM_RHS(SYM_PERM(I)) = I
         ENDDO
      ELSEIF (STRAT .EQ.  2) THEN
         DO I=1, SIZEPERM
            PERM_RHS(SIZEPERM-SYM_PERM(I)+1) = I
         ENDDO
      ENDIF
      END SUBROUTINE DMUMPS_PERMUTE_RHS_AM1
      SUBROUTINE DMUMPS_INTERLEAVE_RHS_AM1(
     &  PERM_RHS, SIZE_PERM,
     &  IPTR_WORKING, SIZE_IPTR_WORKING, WORKING, SIZE_WORKING,
     &  IRHS_PTR,
     &  STEP, SYM_PERM, N, NBRHS,
     &  PROCNODE, NSTEPS, SLAVEF, KEEP199,
     &  behaviour_L0, reorder, n_select, PROKG, MPG 
     &  )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  SIZE_PERM,
     &                        SIZE_IPTR_WORKING,
     &                        IPTR_WORKING(SIZE_IPTR_WORKING),
     &                        SIZE_WORKING,
     &                        WORKING(SIZE_WORKING),
     &                        N,
     &                        IRHS_PTR(N+1),
     &                        STEP(N),
     &                        SYM_PERM(N),
     &                        NBRHS,
     &                        NSTEPS,
     &                        PROCNODE(NSTEPS),
     &                        SLAVEF, KEEP199,
     &                        n_select, MPG
      LOGICAL, INTENT(IN) :: behaviour_L0,
     &                        reorder, PROKG
      INTEGER, INTENT(INOUT) :: PERM_RHS(SIZE_PERM)
      INTEGER :: I, J, K,
     &           entry,            
     &           node,             
     &           SIZE_PERM_WORKING,
     &           NB_NON_EMPTY,     
     &           to_be_found,      
     &           posintmprhs,      
     &           selected,         
     &           local_selected,   
     &           current_proc,     
     &           NPROCS,           
     &           n_pass,           
     &           pass,             
     &           nblocks,          
     &           n_select_loc,     
     &           IERR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_RHS,    
     &                                      PTR_PROCS,  
     &                                      LOAD_PROCS,  
     &                                      IPTR_PERM_WORKING,
     &                                      PERM_WORKING,
     &                                      MYTYPENODE,
     &                                      PERM_PO
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: USED
      LOGICAL :: allow_above_L0
      INTEGER, EXTERNAL :: MUMPS_TYPENODE_ROUGH
      NPROCS = SIZE_IPTR_WORKING - 1
      ALLOCATE(TMP_RHS(SIZE_PERM),
     &         PTR_PROCS(NPROCS),
     &         LOAD_PROCS(NPROCS),
     &         USED(SIZE_PERM),
     &         IPTR_PERM_WORKING(NPROCS+1),
     &         MYTYPENODE(NSTEPS),
     &         STAT=IERR)
      IF(IERR.GT.0) THEN
        WRITE(*,*)'Allocation error in DMUMPS_INTERLEAVE_RHS_AM1'
        CALL MUMPS_ABORT()
      END IF
      DO I=1, NSTEPS
        MYTYPENODE(I) = MUMPS_TYPENODE_ROUGH( PROCNODE(I), KEEP199 )
      ENDDO
      NB_NON_EMPTY = 0
      DO I=1,SIZE_PERM
        IF(IRHS_PTR(I+1)-IRHS_PTR(I).NE.0) THEN 
          NB_NON_EMPTY = NB_NON_EMPTY + 1
        END IF
      END DO
      K = 0
      IPTR_PERM_WORKING(1)=1
      DO I=1,NPROCS
        USED = .FALSE.
        DO J=IPTR_WORKING(I),IPTR_WORKING(I+1)-1
          USED(WORKING(J)) = .TRUE.        
        END DO
        DO J=1,N
          IF (USED(abs(STEP(PERM_RHS(J)))).AND.
     &      ((IRHS_PTR(PERM_RHS(J)+1)-IRHS_PTR(PERM_RHS(J))).NE.0))
     &    THEN
            K = K + 1
          END IF
        END DO
        IPTR_PERM_WORKING(I+1) = K+1
      END DO
      SIZE_PERM_WORKING = K
      ALLOCATE(PERM_WORKING(SIZE_PERM_WORKING),
     &         STAT=IERR)
      IF(IERR.GT.0) THEN
        WRITE(*,*)'Allocation error in DMUMPS_INTERLEAVE_RHS_AM1'
        CALL MUMPS_ABORT()
      END IF
      K = 0
      DO I=1,NPROCS
        USED = .FALSE.
        DO J=IPTR_WORKING(I),IPTR_WORKING(I+1)-1
          USED(WORKING(J)) = .TRUE.        
        END DO
        DO J=1,N
          IF (USED(abs(STEP(PERM_RHS(J)))).AND.
     &      ((IRHS_PTR(PERM_RHS(J)+1)-IRHS_PTR(PERM_RHS(J))).NE.0)) 
     &    THEN
            K = K + 1
            PERM_WORKING(K) = PERM_RHS(J)
          END IF
        END DO
      END DO      
      IF(behaviour_L0) THEN
        n_pass = 2
        allow_above_L0 = .false.
        to_be_found = 0
        DO I=1,SIZE_PERM
          IF((MYTYPENODE(abs(STEP(I))).LE.1).AND. 
     &    (IRHS_PTR(I+1)-IRHS_PTR(I).NE.0))       
     &    THEN
            to_be_found = to_be_found + 1
          END IF
        END DO
      ELSE
        n_pass = 1
        allow_above_L0 = .true.
        to_be_found = NB_NON_EMPTY
      END IF
      PTR_PROCS(1:NPROCS) = IPTR_PERM_WORKING(1:NPROCS)
      LOAD_PROCS = 0
      USED = .FALSE.
      current_proc = 1
      n_select_loc = n_select
      IF (n_select_loc.LE.0) THEN
       n_select_loc = 1
      ENDIF
      posintmprhs = 0
      DO pass=1,n_pass
        selected = 0
        DO WHILE(selected.LT.to_be_found)
          local_selected = 0
          DO WHILE(local_selected.LT.n_select_loc)
            IF(PTR_PROCS(current_proc).EQ.
     &        IPTR_PERM_WORKING(current_proc+1))
     &      THEN 
              EXIT
            ELSE 
              entry = PERM_WORKING(PTR_PROCS(current_proc))
              node  = abs(STEP(entry))
              IF(.NOT.USED(entry)) THEN
                IF(allow_above_L0.OR.(MYTYPENODE(node).LE.1)) THEN
                  USED(entry) = .TRUE.
                  selected = selected + 1
                  local_selected = local_selected + 1
                  posintmprhs = posintmprhs + 1
                  TMP_RHS(posintmprhs) = entry
                  IF(selected.EQ.to_be_found) EXIT
                END IF
              END IF
              PTR_PROCS(current_proc) = PTR_PROCS(current_proc) + 1
            END IF
          END DO
          current_proc = mod(current_proc,NPROCS)+1
        END DO
        to_be_found = NB_NON_EMPTY - to_be_found
        allow_above_L0 = .true.
        PTR_PROCS(1:NPROCS) = IPTR_PERM_WORKING(1:NPROCS)
      END DO
      DO I=1,SIZE_PERM
        IF(IRHS_PTR(PERM_RHS(I)+1)-IRHS_PTR(PERM_RHS(I)).EQ.0) THEN
          posintmprhs = posintmprhs+1
          TMP_RHS(posintmprhs) = PERM_RHS(I)
          IF(posintmprhs.EQ.SIZE_PERM) EXIT
        END IF
      END DO
      IF(reorder) THEN
        posintmprhs = 0 
        ALLOCATE(PERM_PO(N),STAT=IERR)
        IF(IERR.GT.0) THEN
          WRITE(*,*)'Allocation error in INTERLEAVE_RHS_AM1'
          CALL MUMPS_ABORT()
        END IF
        DO J=1,N
          PERM_PO(SYM_PERM(J))=J
        END DO
        nblocks = N/NBRHS
        DO I = 1, nblocks
          USED = .FALSE.
          DO J=1, NBRHS
            USED(TMP_RHS(NBRHS*(I-1)+J))=.TRUE.
          END DO
          DO J=1,N
            IF(USED(PERM_PO(J))) THEN
              posintmprhs = posintmprhs + 1
              PERM_RHS(posintmprhs) = PERM_PO(J)
            END IF
          END DO
        END DO
        IF(mod(N,NBRHS).NE.0) THEN
          USED = .FALSE.
          DO J=1, mod(N,NBRHS)
            USED(TMP_RHS(NBRHS*nblocks+J))=.TRUE.
          END DO
          DO J=1,N
            IF(USED(PERM_PO(J))) THEN
              posintmprhs = posintmprhs + 1
              PERM_RHS(posintmprhs) = PERM_PO(J)
            END IF
          END DO
        END IF
        DEALLOCATE(PERM_PO)
      ELSE
        PERM_RHS = TMP_RHS
      END IF
      DEALLOCATE(TMP_RHS,
     &           PTR_PROCS,
     &           LOAD_PROCS,
     &           USED,
     &           IPTR_PERM_WORKING,
     &           PERM_WORKING,
     &           MYTYPENODE)
      RETURN
      END SUBROUTINE DMUMPS_INTERLEAVE_RHS_AM1
