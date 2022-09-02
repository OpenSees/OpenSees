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
      SUBROUTINE DMUMPS_SOL_C(root, N, A, LA, IW, LIW, W, LWC, 
     & IWCB, LIWW, NRHS, NA, LNA, NE_STEPS, W2, MTYPE, ICNTL, FROM_PP,
     & STEP, FRERE, DAD, FILS, PTRIST, PTRFAC, IW1, LIW1, PTRACB,
     & LIWK_PTRACB, PROCNODE_STEPS, SLAVEF, INFO, KEEP,KEEP8, DKEEP,
     & COMM_NODES, MYID, MYID_NODES, BUFR, LBUFR, LBUFR_BYTES,
     & ISTEP_TO_INIV2, TAB_POS_IN_PERE, IBEG_ROOT_DEF, IEND_ROOT_DEF,
     & IROOT_DEF_RHS_COL1, RHS_ROOT, LRHS_ROOT, SIZE_ROOT, MASTER_ROOT,
     & RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD, POSINRHSCOMP_BWD,
     & NZ_RHS, NBCOL_INBLOC, NRHS_ORIG, JBEG_RHS, Step2node, LStep2node,
     & IRHS_SPARSE, IRHS_PTR, SIZE_PERM_RHS, PERM_RHS,
     & SIZE_UNS_PERM_INV, UNS_PERM_INV, NB_FS_IN_RHSCOMP_F,
     & NB_FS_IN_RHSCOMP_TOT, DO_NBSPARSE , RHS_BOUNDS, LRHS_BOUNDS
     & )
      USE DMUMPS_OOC
      USE DMUMPS_SOL_ES
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
#if defined(V_T)
      INCLUDE 'VT.inc'
#endif
      TYPE ( DMUMPS_ROOT_STRUC ) :: root
      INTEGER(8) :: LA
      INTEGER(8) :: LWC
      INTEGER    :: N,LIW,MTYPE,LIW1,LIWW,LNA
      INTEGER ICNTL(60),INFO(80), KEEP(500)
      DOUBLE PRECISION, intent(inout)     :: DKEEP(230)
      INTEGER(8) KEEP8(150)
      INTEGER IW(LIW),IW1(LIW1),NA(LNA),NE_STEPS(KEEP(28)),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)), FILS(N), PTRIST(KEEP(28)),
     &        DAD(KEEP(28))
      INTEGER(8) ::  PTRFAC(KEEP(28))
      INTEGER    ::  LIWK_PTRACB
      INTEGER(8) ::  PTRACB(LIWK_PTRACB)
      INTEGER NRHS, LRHSCOMP, NB_FS_IN_RHSCOMP_F, NB_FS_IN_RHSCOMP_TOT
      DOUBLE PRECISION    A(LA), W(LWC), 
     &        W2(KEEP(133))
      DOUBLE PRECISION ::  RHSCOMP(LRHSCOMP,NRHS)
      INTEGER SLAVEF, COMM_NODES, MYID, MYID_NODES
      INTEGER PROCNODE_STEPS(KEEP(28)), POSINRHSCOMP_FWD(N), 
     &        POSINRHSCOMP_BWD(N)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER IBEG_ROOT_DEF, IEND_ROOT_DEF, IROOT_DEF_RHS_COL1
      INTEGER SIZE_ROOT, MASTER_ROOT
      INTEGER(8) :: LRHS_ROOT
      DOUBLE PRECISION RHS_ROOT(LRHS_ROOT)
      LOGICAL, intent(in) :: FROM_PP 
      INTEGER, intent(in)      :: NZ_RHS, NBCOL_INBLOC, NRHS_ORIG
      INTEGER, intent(in)      :: SIZE_UNS_PERM_INV 
      INTEGER, intent(in)      :: SIZE_PERM_RHS 
      INTEGER, intent(in) :: JBEG_RHS
      INTEGER, intent(in) :: IRHS_SPARSE(NZ_RHS)
      INTEGER, intent(in) :: IRHS_PTR(NBCOL_INBLOC+1)
      INTEGER, intent(in) :: PERM_RHS(SIZE_PERM_RHS)
      INTEGER, intent(in) :: UNS_PERM_INV(SIZE_UNS_PERM_INV)
      INTEGER, intent(in) :: LStep2node
      INTEGER, intent(in) :: Step2node(LStep2node)
      LOGICAL, intent(in)  :: DO_NBSPARSE  
      INTEGER, intent(in)     :: LRHS_BOUNDS
      INTEGER, intent(inout)  :: RHS_BOUNDS (LRHS_BOUNDS) 
      INTEGER MP, LP, LDIAG
      INTEGER K,I,II
      INTEGER allocok
      INTEGER LPOOL,MYLEAF,MYROOT,NBROOT,LPANEL_POS
      INTEGER NSTK_S,IPOOL,IPANEL_POS,PTRICB
      INTEGER MTYPE_LOC
      INTEGER MODE_RHS_BOUNDS 
      INTEGER IERR
      INTEGER(8) :: IAPOS
      INTEGER       IOLDPS,
     &              LOCAL_M,
     &              LOCAL_N
#if defined(V_T)
      INTEGER soln_c_class, forw_soln, back_soln, root_soln
#endif
      LOGICAL DOFORWARD, DOROOT, DOBACKWARD
      LOGICAL I_WORKED_ON_ROOT, SPECIAL_ROOT_REACHED
      INTEGER IROOT
      LOGICAL DOROOT_FWD_OOC, DOROOT_BWD_PANEL
      LOGICAL SWITCH_OFF_ES
      LOGICAL DUMMY_BOOL
      INTEGER :: IDUMMY
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
      INCLUDE 'mumps_headers.h'
      INTEGER, DIMENSION(:), ALLOCATABLE ::  nodes_RHS
      INTEGER nb_nodes_RHS
      INTEGER nb_prun_leaves
      INTEGER, DIMENSION(:), ALLOCATABLE :: Pruned_Leaves
      INTEGER, DIMENSION(:), ALLOCATABLE ::  Pruned_List
      INTEGER  nb_prun_nodes
      INTEGER nb_prun_roots, JAM1
      INTEGER, DIMENSION(:), ALLOCATABLE ::  Pruned_SONS, Pruned_Roots
      INTEGER                            :: SIZE_TO_PROCESS
      LOGICAL, DIMENSION(:), ALLOCATABLE :: TO_PROCESS
      INTEGER ISTEP, INODE_PRINC
      LOGICAL AM1, DO_PRUN
      LOGICAL Exploit_Sparsity
      LOGICAL DO_NBSPARSE_BWD, PRUN_BELOW_BWD
      INTEGER :: OOC_FCT_TYPE_TMP
      INTEGER :: MUMPS_OOC_GET_FCT_TYPE
      EXTERNAL :: MUMPS_OOC_GET_FCT_TYPE
      DOUBLE PRECISION TIME_FWD,TIME_BWD,TIME_SpecialRoot
      INTEGER :: nb_sparse
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      LOGICAL, EXTERNAL :: MUMPS_IN_OR_ROOT_SSARBR
      MYLEAF = -1
      LP      = ICNTL(1)
      MP      = ICNTL(2)
      LDIAG   = ICNTL(4)
#if defined(V_T)
      CALL VTCLASSDEF( 'Soln_c',soln_c_class,ierr)
      CALL VTFUNCDEF( 'forw_soln',soln_c_class,forw_soln,ierr)
      CALL VTFUNCDEF( 'back_soln',soln_c_class,back_soln,ierr)
      CALL VTFUNCDEF( 'root_soln',soln_c_class,root_soln,ierr)
#endif
      IF (.NOT. FROM_PP) THEN
        CALL MUMPS_SECDEB(TIME_FWD)
      ENDIF
      NSTK_S   = 1
      PTRICB = NSTK_S + KEEP(28)
      IPOOL = PTRICB + KEEP(28)
      LPOOL = NA(1) + 1
      IPANEL_POS = IPOOL + LPOOL
      IF (KEEP(201).EQ.1) THEN
        LPANEL_POS = KEEP(228)+1
      ELSE
        LPANEL_POS = 1
      ENDIF
      IF (IPANEL_POS + LPANEL_POS -1 .ne. LIW1 )  THEN
         WRITE(*,*)  MYID, ": Internal Error 1 in DMUMPS_SOL_C",
     &   IPANEL_POS, LPANEL_POS, LIW1
         CALL MUMPS_ABORT()
      ENDIF
      DOFORWARD = .TRUE.
      DOBACKWARD= .TRUE.
      SPECIAL_ROOT_REACHED = .TRUE.
      SWITCH_OFF_ES    = .FALSE.
      IF ( KEEP(111).NE.0 .OR. KEEP(252).NE.0 ) THEN
        DOFORWARD = .FALSE.
      ENDIF
      IF (KEEP(221).eq.1) DOBACKWARD = .FALSE.
      IF (KEEP(221).eq.2) DOFORWARD  = .FALSE.
      IF ( KEEP(60).EQ.0 .AND.
     &    ( 
     &      (KEEP(38).NE.0 .AND.  root%yes) 
     &  .OR.
     &      (KEEP(20).NE.0 .AND. MYID_NODES.EQ.MASTER_ROOT)
     &    ) 
     &  .AND. KEEP(252).EQ.0
     &   )
     &THEN
        DOROOT = .TRUE.
      ELSE
        DOROOT = .FALSE.
      ENDIF
      DOROOT_BWD_PANEL = DOROOT .AND. MTYPE.NE.1 .AND. KEEP(50).EQ.0
     &                     .AND. KEEP(201).EQ.1
      DOROOT_FWD_OOC = DOROOT .AND. .NOT.DOROOT_BWD_PANEL
      AM1              = (KEEP(237) .NE. 0)
      Exploit_Sparsity = (KEEP(235) .NE. 0) .AND. (.NOT. AM1)
      DO_PRUN          = (Exploit_Sparsity.OR.AM1)
      IF (FROM_PP) THEN
        Exploit_Sparsity = .FALSE.
        DO_PRUN          = .FALSE.
           IF ( AM1 ) THEN
          WRITE(*,*) "Internal error 2 in DMUMPS_SOL_C"
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      IF ( DO_PRUN ) THEN
         ALLOCATE (Pruned_SONS(KEEP(28)), stat=I)
         IF(I.GT.0) THEN
               INFO(1)=-13
               INFO(2)=KEEP(28)
         END IF
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
      ENDIF
      IF ( DO_PRUN 
     &   ) THEN
         SIZE_TO_PROCESS = KEEP(28)
      ELSE
         SIZE_TO_PROCESS = 1
      ENDIF
      ALLOCATE (TO_PROCESS(SIZE_TO_PROCESS), stat=I)
      IF(I.GT.0) THEN
          INFO(1)=-13
          INFO(2)=KEEP(28)
      END IF
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
      IF(INFO(1).LT.0) GOTO 500
      IF ( DOFORWARD .AND. DO_PRUN ) THEN
         nb_prun_nodes = 0
         nb_prun_roots = 0
         Pruned_SONS(:) = -1
         IF ( Exploit_Sparsity ) THEN
            nb_nodes_RHS = 0
            DO I = 1, NZ_RHS
               ISTEP       = abs( STEP(IRHS_SPARSE(I)) )
               INODE_PRINC = Step2node( ISTEP )
               IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
                  nb_nodes_RHS = nb_nodes_RHS +1
                  Pruned_SONS(ISTEP) = 0 
               ENDIF
            ENDDO
            ALLOCATE(nodes_RHS(nb_nodes_RHS), STAT = allocok)
            IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_nodes_RHS
            END IF
            CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
            IF(INFO(1).LT.0) GOTO 500
            nb_nodes_RHS = 0
            Pruned_SONS = -1
            DO I = 1, NZ_RHS
               ISTEP       = abs( STEP(IRHS_SPARSE(I)) )
               INODE_PRINC = Step2node( ISTEP )
               IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
                  nb_nodes_RHS = nb_nodes_RHS +1
                  nodes_RHS(nb_nodes_RHS)  = INODE_PRINC
                  Pruned_SONS(ISTEP) = 0 
               ENDIF
            ENDDO
         ELSE IF ( AM1 ) THEN  
            nb_nodes_RHS = 0
            DO I = 1, NBCOL_INBLOC
              IF ( (IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE
              IF ( (KEEP(242) .NE. 0 ).OR. (KEEP(243).NE.0) ) THEN
                   JAM1 = PERM_RHS(JBEG_RHS+I-1)
              ELSE
                   JAM1 = JBEG_RHS+I-1
              ENDIF       
              ISTEP = abs(STEP(JAM1))
              INODE_PRINC = Step2node(ISTEP)
              IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
                 nb_nodes_RHS = nb_nodes_RHS +1
                 Pruned_SONS(ISTEP) = 0                 
              ENDIF
            ENDDO
            ALLOCATE(nodes_RHS(nb_nodes_RHS), STAT = allocok)
            IF(allocok.GT.0) THEN
               INFO(1)=-13
               INFO(2)=nb_nodes_RHS
            END IF
            CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
            IF(INFO(1).LT.0) GOTO 500
            nb_nodes_RHS = 0
            Pruned_SONS = -1
            DO I = 1, NBCOL_INBLOC
              IF ( (IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE
              IF ( (KEEP(242) .NE. 0 ).OR. (KEEP(243).NE.0) ) THEN
                   JAM1 = PERM_RHS(JBEG_RHS+I-1)
              ELSE
                   JAM1 = JBEG_RHS+I-1
              ENDIF
              ISTEP = abs(STEP(JAM1))
              INODE_PRINC = Step2node(ISTEP)
              IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
                 nb_nodes_RHS = nb_nodes_RHS +1
                 nodes_RHS(nb_nodes_RHS)  = INODE_PRINC
                 Pruned_SONS(ISTEP) = 0
              ENDIF
            ENDDO            
         ENDIF                  
         CALL DMUMPS_CHAIN_PRUN_NODES( 
     &        .FALSE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_RHS, nb_nodes_RHS,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots,
     &        nb_prun_leaves )  
         ALLOCATE(Pruned_List(nb_prun_nodes), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_nodes
         END IF
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         ALLOCATE(Pruned_Roots(nb_prun_roots), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_roots
         END IF
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         ALLOCATE(Pruned_Leaves(nb_prun_leaves), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_leaves
         END IF
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         CALL DMUMPS_CHAIN_PRUN_NODES( 
     &        .TRUE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_RHS, nb_nodes_RHS,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &        Pruned_List, Pruned_Roots, Pruned_Leaves )
         IF(allocated(nodes_RHS)) DEALLOCATE(nodes_RHS)
         CALL DMUMPS_OOC_SET_STATES_ES(N,
     &          KEEP(201), Pruned_List, nb_prun_nodes,
     &          STEP)
         IF ( KEEP(201) .GT. 0) THEN
           OOC_FCT_TYPE_TMP=MUMPS_OOC_GET_FCT_TYPE
     &                      ('F',MTYPE,KEEP(201),KEEP(50))
         ELSE
           OOC_FCT_TYPE_TMP = -5959 
         ENDIF
         CALL DMUMPS_CHAIN_PRUN_NODES_STATS(
     &        MYID_NODES, N, KEEP(28), KEEP(201), KEEP(485), 
     &        KEEP8(31), 
     &        STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_TMP
     &   )
         IF (DO_NBSPARSE) THEN
           nb_sparse = max(1,KEEP(497))
           MODE_RHS_BOUNDS = 0 
           IF (Exploit_Sparsity) MODE_RHS_BOUNDS = 2
           CALL DMUMPS_INITIALIZE_RHS_BOUNDS(
     &      STEP, N,
     &      IRHS_PTR, NBCOL_INBLOC, IRHS_SPARSE, NZ_RHS,
     &      JBEG_RHS, PERM_RHS, SIZE_PERM_RHS, KEEP(242), KEEP(243),
     &      UNS_PERM_INV, SIZE_UNS_PERM_INV, KEEP(23),
     &      RHS_BOUNDS, KEEP(28),
     &      nb_sparse, MYID_NODES,
     &      MODE_RHS_BOUNDS) 
           CALL DMUMPS_PROPAGATE_RHS_BOUNDS(
     &      Pruned_Leaves, nb_prun_leaves,
     &      STEP, N, Pruned_SONS,
     &      DAD, RHS_BOUNDS, KEEP(28),
     &      MYID_NODES, COMM_NODES, KEEP(485),
     &      IW, LIW, PTRIST,KEEP(IXSZ),OOC_FCT_TYPE_TMP,0,
     &      KEEP(50), KEEP(38))
         END IF
         SPECIAL_ROOT_REACHED = .FALSE.
         DO I= 1, nb_prun_roots
          IF ( (Pruned_Roots(I).EQ.KEEP(38)).OR.
     &         (Pruned_Roots(I).EQ.KEEP(20)) ) THEN
            SPECIAL_ROOT_REACHED = .TRUE.
            EXIT
          ENDIF
         ENDDO
         DEALLOCATE(Pruned_List)
      ENDIF  
      IF (KEEP(201).GT.0) THEN
        IF (DOFORWARD .OR. DOROOT_FWD_OOC) THEN
           CALL DMUMPS_SOLVE_INIT_OOC_FWD(PTRFAC,KEEP(28),MTYPE,
     &                                A,LA,DOFORWARD,IERR)
           IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            CALL MUMPS_ABORT()
          ENDIF
        ENDIF
      ENDIF
      IF (DOFORWARD) THEN
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          MTYPE_LOC = MTYPE
        ELSE
          MTYPE_LOC = 1
        ENDIF
#if defined(V_T)
        CALL VTBEGIN(forw_soln,ierr)
#endif
        IF ( .NOT. DO_PRUN ) THEN
          CALL MUMPS_INIT_NROOT_DIST(N, NBROOT, MYROOT, MYID_NODES,
     &           SLAVEF, NA, LNA, KEEP, STEP, PROCNODE_STEPS)
          DO ISTEP =1, KEEP(28)
            IW1(NSTK_S+ISTEP-1) = NE_STEPS(ISTEP)
          ENDDO
        ELSE
          CALL MUMPS_NBLOCAL_ROOTS_OR_LEAVES( N,
     &         nb_prun_roots, Pruned_Roots,
     &         MYROOT, MYID_NODES, SLAVEF, KEEP, STEP,
     &         PROCNODE_STEPS )
          IF (AM1) THEN 
            DEALLOCATE(Pruned_Roots)
          END IF
          IF ((Exploit_Sparsity).AND.(nb_prun_roots.EQ.NA(2))) THEN
            DEALLOCATE(Pruned_Roots)
            SWITCH_OFF_ES = .TRUE.   
          ENDIF
          DO ISTEP = 1, KEEP(28)
            IW1(NSTK_S+ISTEP-1) = Pruned_SONS(ISTEP)
          ENDDO
        ENDIF
          IF ( DO_PRUN ) THEN
            CALL MUMPS_INIT_POOL_DIST_NONA( N, MYLEAF, MYID_NODES,
     &            nb_prun_leaves, Pruned_Leaves, KEEP, KEEP8,
     &            STEP, PROCNODE_STEPS, IW1(IPOOL), LPOOL )
            MYLEAF = MYLEAF - 1
            DEALLOCATE(Pruned_Leaves)
          ELSE
            CALL MUMPS_INIT_POOL_DIST( N, MYLEAF, MYID_NODES,
     &             SLAVEF, NA, LNA, KEEP, KEEP8, STEP,
     &             PROCNODE_STEPS, IW1(IPOOL), LPOOL )
            MYLEAF = MYLEAF - 1 
          ENDIF
        CALL DMUMPS_SOL_R(N, A(1), LA, IW(1), LIW, W(1),
     &       LWC, NRHS,
     &       IW1(PTRICB), IWCB, LIWW,
     &       RHSCOMP,LRHSCOMP,POSINRHSCOMP_FWD,
     &       STEP, FRERE,DAD,FILS,
     &       IW1(NSTK_S),IW1(IPOOL),LPOOL,PTRIST,PTRFAC,
     &       MYLEAF, MYROOT, INFO,
     &       KEEP, KEEP8, DKEEP,
     &       PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     &       BUFR, LBUFR, LBUFR_BYTES,
     &       RHS_ROOT, LRHS_ROOT, MTYPE_LOC, 
     &       
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP
     &       )
        IF (DO_PRUN) THEN
          MYLEAF = -1
        ENDIF
#if defined(V_T)
        CALL VTEND(forw_soln,ierr)
#endif
      ENDIF
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
      IF ( INFO(1) .LT. 0 ) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP,*) MYID,
     &    ': ** ERROR RETURN FROM DMUMPS_SOL_R,INFO(1:2)=',
     &    INFO(1:2)
        END IF
        GOTO 500   
      END IF
      CALL MPI_BARRIER( COMM_NODES, IERR )
      IF (.NOT.FROM_PP) THEN
         CALL MUMPS_SECFIN(TIME_FWD)
         DKEEP(117)=TIME_FWD + DKEEP(117)
      ENDIF
      IF (DO_PRUN.AND.SWITCH_OFF_ES) THEN
         DO_PRUN          = .FALSE.
         Exploit_Sparsity = .FALSE.
           IF ( allocated(TO_PROCESS) .AND. SIZE_TO_PROCESS.NE.1 ) THEN
             DEALLOCATE (TO_PROCESS)
             SIZE_TO_PROCESS = 1
             ALLOCATE(TO_PROCESS(SIZE_TO_PROCESS),stat=I)
           ENDIF
      ENDIF
      IF ( DOBACKWARD .AND. DO_PRUN )  THEN
         nb_prun_leaves = 0
        IF ( Exploit_Sparsity .AND. (KEEP(111).EQ.0) ) THEN
          nb_nodes_RHS = nb_prun_roots
          ALLOCATE(nodes_RHS(nb_nodes_RHS), STAT = allocok)
          IF(allocok.GT.0) THEN
            WRITE(*,*)'Problem with allocation of nodes_RHS'
            INFO(1) = -13
            INFO(2) = nb_nodes_RHS
            CALL MUMPS_ABORT()
          END IF
          nodes_RHS(1:nb_prun_roots)=Pruned_Roots(1:nb_prun_roots)
          DEALLOCATE(Pruned_Roots)
        ELSE
         nb_nodes_RHS = 0
         Pruned_SONS(:) = -1 
         DO II = 1, NZ_RHS
            I = IRHS_SPARSE(II)
            IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)          
            ISTEP = abs(STEP(I))
            IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
               nb_nodes_RHS = nb_nodes_RHS +1
               Pruned_SONS(ISTEP) = 0
            ENDIF
         ENDDO
         ALLOCATE(nodes_RHS(nb_nodes_RHS), STAT = allocok)           
         IF(allocok.GT.0) THEN
           WRITE(*,*)'Problem with allocation of nodes_RHS'
           INFO(1) = -13
           INFO(2) = nb_nodes_RHS
           CALL MUMPS_ABORT()
         END IF
         nb_nodes_RHS = 0         
         Pruned_SONS(:) = -1  
         DO II = 1, NZ_RHS
            I = IRHS_SPARSE(II)
            IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
            ISTEP = abs(STEP(I))
            INODE_PRINC = Step2node(ISTEP)
            IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
               nb_nodes_RHS = nb_nodes_RHS +1
               nodes_RHS(nb_nodes_RHS)  = INODE_PRINC
               Pruned_SONS(ISTEP) = 0
            ENDIF
         ENDDO
        ENDIF
        IF ( Exploit_Sparsity ) THEN
           CALL DMUMPS_TREE_PRUN_NODES( 
     &     .FALSE.,
     &     DAD, NE_STEPS, FRERE, KEEP(28),
     &     FILS, STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     TO_PROCESS,
     &     nb_prun_nodes, nb_prun_roots, nb_prun_leaves
     &     )
           ALLOCATE(Pruned_List(nb_prun_nodes), STAT=allocok)
           IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_prun_nodes
           END IF
           CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
           IF(INFO(1).LT.0) GOTO 500
           ALLOCATE(Pruned_Roots(nb_prun_roots), STAT=allocok)
           IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_prun_roots
           END IF
           CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
           IF(INFO(1).LT.0) GOTO 500
           ALLOCATE(Pruned_Leaves(nb_prun_leaves), STAT=allocok)
           IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_prun_leaves
           END IF
           CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
           IF(INFO(1).LT.0) GOTO 500
           CALL DMUMPS_TREE_PRUN_NODES( 
     &     .TRUE.,
     &     DAD, NE_STEPS, FRERE, KEEP(28),
     &     FILS, STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     TO_PROCESS,
     &     nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &     Pruned_List, Pruned_Roots, Pruned_Leaves
     &     )
           CALL DMUMPS_OOC_SET_STATES_ES(N,
     &          KEEP(201), Pruned_List, nb_prun_nodes,
     &          STEP)
           IF(allocated(nodes_RHS)) DEALLOCATE(nodes_RHS)
           IF (KEEP(201).GT.0) THEN
             OOC_FCT_TYPE_TMP=MUMPS_OOC_GET_FCT_TYPE
     &                    ('B',MTYPE,KEEP(201),KEEP(50))
           ELSE
             OOC_FCT_TYPE_TMP = -5959 
           ENDIF
           CALL DMUMPS_TREE_PRUN_NODES_STATS(
     &          MYID_NODES, N, KEEP(28), KEEP(201),
     &        KEEP8(31), 
     &          STEP,
     &          Pruned_List,
     &          nb_prun_nodes, OOC_FCT_TYPE_TMP)
        ENDIF
      ENDIF
      IF(KEEP(201).EQ.1.AND.DOROOT_BWD_PANEL) THEN
         I_WORKED_ON_ROOT = .FALSE. 
         CALL DMUMPS_SOLVE_INIT_OOC_BWD(PTRFAC,KEEP(28),MTYPE,
     &   I_WORKED_ON_ROOT, IROOT, A, LA, IERR)
         IF (IERR .LT. 0) THEN
           INFO(1) = -90
           INFO(2) = IERR
         ENDIF 
      ENDIF
      IF (KEEP(201).EQ.1) THEN
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF ( INFO(1) .LT. 0 ) GOTO 500  
      ENDIF
      IF (KEEP(60).NE.0 .AND. KEEP(221).EQ.0
     &   .AND. MYID_NODES .EQ. MASTER_ROOT) THEN
        RHS_ROOT(1:NRHS*SIZE_ROOT) = ZERO
      ENDIF
      IF (.NOT. FROM_PP) THEN
        CALL MUMPS_SECDEB(TIME_SpecialRoot)
      ENDIF
      IF ( ( KEEP( 38 ) .NE. 0 ).AND. SPECIAL_ROOT_REACHED ) THEN
        IF ( KEEP(60) .EQ. 0 .AND. KEEP(252) .EQ. 0 ) THEN
          IF ( root%yes ) THEN
            IF (KEEP(201).GT.0) THEN
              IF ( (Exploit_Sparsity.AND.(KEEP(111).NE.0)) .and.
     &            (OOC_STATE_NODE(STEP(KEEP(38))).eq.-6) ) THEN
                  GOTO 1010
              ENDIF
            ENDIF
            IOLDPS = PTRIST(STEP(KEEP(38)))
            LOCAL_M = IW( IOLDPS + 2 + KEEP(IXSZ))
            LOCAL_N = IW( IOLDPS + 1 + KEEP(IXSZ))
            IF (KEEP(201).GT.0) THEN 
              CALL DMUMPS_SOLVE_GET_OOC_NODE(
     &           KEEP(38),PTRFAC,KEEP,A,LA,
     &           STEP,KEEP8,N,DUMMY_BOOL,IERR)
              IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                WRITE(*,*) '** ERROR after DMUMPS_SOLVE_GET_OOC_NODE',
     &          INFO(1)
                call MUMPS_ABORT()
              ENDIF
            ENDIF
            IAPOS   = PTRFAC(IW( IOLDPS + 4 + KEEP(IXSZ)))
            IF (LOCAL_M * LOCAL_N .EQ. 0) THEN
              IAPOS = min(IAPOS, LA)
            ENDIF
#if defined(V_T)
            CALL VTBEGIN(root_soln,ierr)
#endif
             CALL DMUMPS_ROOT_SOLVE( NRHS, root%DESCRIPTOR(1), 
     &       root%CNTXT_BLACS, LOCAL_M, LOCAL_N,
     &       root%MBLOCK, root%NBLOCK,
     &       root%IPIV(1), root%LPIV, MASTER_ROOT, MYID_NODES,
     &       COMM_NODES,
     &       RHS_ROOT(1),
     &       root%TOT_ROOT_SIZE, A( IAPOS ),
     &       INFO(1), MTYPE, KEEP(50), FROM_PP)
            IF(KEEP(201).GT.0)THEN 
              CALL DMUMPS_FREE_FACTORS_FOR_SOLVE(KEEP(38),
     &             PTRFAC,KEEP(28),A,LA,.FALSE.,IERR)
              IF(IERR.LT.0)THEN
                 INFO(1)=IERR
                 INFO(2)=0
                 WRITE(*,*)
     &           '** ERROR after DMUMPS_FREE_FACTORS_FOR_SOLVE ',
     &           INFO(1)
                 call MUMPS_ABORT()
              ENDIF
            ENDIF
          ENDIF  
        ENDIF
      ELSE IF ( ( KEEP(20) .NE. 0) .AND. SPECIAL_ROOT_REACHED ) THEN
        IF ( MYID_NODES .eq.  MASTER_ROOT ) THEN
        END IF 
      END IF 
      IF (.NOT.FROM_PP) THEN
         CALL MUMPS_SECFIN(TIME_SpecialRoot)
         DKEEP(119)=TIME_SpecialRoot + DKEEP(119)
      ENDIF
#if defined(V_T)
      CALL VTEND(root_soln,ierr)
#endif
 1010 CONTINUE
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      IF (DOBACKWARD) THEN
        IF ( KEEP(201).GT.0 .AND.  .NOT. DOROOT_BWD_PANEL )
     &    THEN
          I_WORKED_ON_ROOT = DOROOT
          IF (KEEP(38).gt.0 ) THEN
             IF ( ( Exploit_Sparsity.AND.(KEEP(111).EQ.0) )
     &            .OR. AM1 ) THEN
                IF (OOC_STATE_NODE(STEP(KEEP(38))).eq.-6) THEN
                   OOC_STATE_NODE(STEP(KEEP(38)))=-4
                ENDIF
             ENDIF
             IF (Exploit_Sparsity.AND.(KEEP(111).NE.0)) THEN
                IF (OOC_STATE_NODE(STEP(KEEP(38))).eq.-6) THEN
                   I_WORKED_ON_ROOT = .FALSE.
                ENDIF
             ENDIF
          ENDIF
        ENDIF                    
        IF (.NOT.AM1) THEN
          DO_NBSPARSE_BWD = .FALSE.
        ELSE
          DO_NBSPARSE_BWD = DO_NBSPARSE
        ENDIF
        PRUN_BELOW_BWD = AM1
        IF ( AM1 ) THEN
          CALL DMUMPS_CHAIN_PRUN_NODES( 
     &        .FALSE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_RHS, nb_nodes_RHS,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots,
     &        nb_prun_leaves)  
          ALLOCATE(Pruned_List(nb_prun_nodes), STAT=allocok)
          IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_nodes
          END IF
          CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
          IF(INFO(1).LT.0) GOTO 500
          ALLOCATE(Pruned_Roots(nb_prun_roots), STAT=allocok)
          IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_roots
          END IF
          CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
          IF(INFO(1).LT.0) GOTO 500
          ALLOCATE(Pruned_Leaves(nb_prun_leaves), STAT=allocok)
          IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_leaves
          END IF
          CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
          IF(INFO(1).LT.0) GOTO 500
          CALL DMUMPS_CHAIN_PRUN_NODES( 
     &        .TRUE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_RHS, nb_nodes_RHS,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &        Pruned_List, Pruned_Roots, Pruned_Leaves )
          CALL DMUMPS_OOC_SET_STATES_ES(N,
     &          KEEP(201), Pruned_List, nb_prun_nodes,
     &          STEP)
          IF (KEEP(201).GT.0) THEN
           OOC_FCT_TYPE_TMP=MUMPS_OOC_GET_FCT_TYPE
     &                    ('B',MTYPE,KEEP(201),KEEP(50))
          ELSE
           OOC_FCT_TYPE_TMP = -5959 
          ENDIF
          CALL DMUMPS_CHAIN_PRUN_NODES_STATS(
     &    MYID_NODES, N, KEEP(28), KEEP(201), KEEP(485), KEEP8(31),
     &    STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_TMP
     &    )
          IF (DO_NBSPARSE_BWD) THEN
            nb_sparse = max(1,KEEP(497))
            CALL DMUMPS_INITIALIZE_RHS_BOUNDS(
     &       STEP, N,
     &       IRHS_PTR, NBCOL_INBLOC, IRHS_SPARSE, NZ_RHS,
     &       JBEG_RHS, PERM_RHS, SIZE_PERM_RHS, KEEP(242), KEEP(243),
     &       UNS_PERM_INV, SIZE_UNS_PERM_INV, KEEP(23),
     &       RHS_BOUNDS, KEEP(28),
     &       nb_sparse, MYID_NODES,
     &       1) 
            CALL DMUMPS_PROPAGATE_RHS_BOUNDS(
     &       Pruned_Leaves, nb_prun_leaves,
     &       STEP, N, Pruned_SONS,
     &       DAD, RHS_BOUNDS, KEEP(28),
     &       MYID_NODES, COMM_NODES, KEEP(485),
     &       IW, LIW, PTRIST,KEEP(IXSZ),OOC_FCT_TYPE_TMP,1,
     &       KEEP(50), KEEP(38))
         END IF
        ENDIF
        IF ( KEEP(201).GT.0 ) THEN
          IROOT = max(KEEP(20),KEEP(38)) 
          CALL DMUMPS_SOLVE_INIT_OOC_BWD(PTRFAC,KEEP(28),MTYPE,
     &         I_WORKED_ON_ROOT, IROOT, A, LA, IERR)
        ENDIF
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          MTYPE_LOC = MTYPE
        ELSE
          MTYPE_LOC = 0
        ENDIF
#if defined(V_T)
        CALL VTBEGIN(back_soln,ierr)
#endif
        IF (.NOT.FROM_PP) THEN
          CALL MUMPS_SECDEB(TIME_BWD)
        ENDIF
        IF ( .NOT.SPECIAL_ROOT_REACHED ) THEN
          RHS_ROOT(1:NRHS*SIZE_ROOT) = ZERO
        ENDIF
        IF (AM1.AND.(NB_FS_IN_RHSCOMP_F.NE.NB_FS_IN_RHSCOMP_TOT)) THEN
         DO I =1, N
           II = POSINRHSCOMP_BWD(I)
           IF ((II.GT.0).AND.(II.GT.NB_FS_IN_RHSCOMP_F)) THEN
            DO K=1,NRHS
             RHSCOMP(II, K) = ZERO
            ENDDO
           ENDIF
         ENDDO
        ENDIF
        IF ( .NOT. DO_PRUN ) THEN
             CALL MUMPS_INIT_POOL_DIST_NA_BWD( N, MYROOT, MYID_NODES,
     &          NA, LNA, KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &          IW1(IPOOL), LPOOL )
             IF (MYLEAF .EQ. -1) THEN
               CALL MUMPS_NBLOCAL_ROOTS_OR_LEAVES( N,
     &         NA(1), 
     &         NA(3), 
     &         MYLEAF, MYID_NODES, SLAVEF, KEEP, STEP,
     &         PROCNODE_STEPS )
             ENDIF
        ELSE 
            CALL MUMPS_INIT_POOL_DIST_BWD(N, nb_prun_roots,
     &      Pruned_Roots,
     &      MYROOT, MYID_NODES, KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &      IW1(IPOOL), LPOOL)
            CALL MUMPS_NBLOCAL_ROOTS_OR_LEAVES( N,
     &       nb_prun_leaves, Pruned_Leaves,
     &       MYLEAF, MYID_NODES, SLAVEF, KEEP, STEP,
     &       PROCNODE_STEPS )
        ENDIF
        IF (KEEP(31) .EQ. 1) THEN
          DO I = 1, KEEP(28)
            IF (MUMPS_PROCNODE(PROCNODE_STEPS(I),KEEP(199)) .EQ.
     &      MYID_NODES) THEN
              IF ( .NOT. MUMPS_IN_OR_ROOT_SSARBR(PROCNODE_STEPS(I),
     &        KEEP(199)) ) THEN
                  IF ( DO_PRUN
     &               ) THEN
                    IF ( TO_PROCESS(I) ) THEN
                      KEEP(31) = KEEP(31) + 1
                    ENDIF
                  ELSE
                    KEEP(31) = KEEP(31) + 1
                  ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        CALL DMUMPS_SOL_S( N, A, LA, IW, LIW, W(1), LWC,
     &          NRHS,
     &          RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD,
     &          IW1(PTRICB),PTRACB,IWCB,LIWW, W2,
     &          NE_STEPS,
     &          STEP, FRERE,DAD,FILS,
     &          IW1(IPOOL),LPOOL,PTRIST,PTRFAC,MYLEAF,MYROOT,ICNTL,INFO,
     &          PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     &          BUFR, LBUFR, LBUFR_BYTES, KEEP, KEEP8, DKEEP,
     &          RHS_ROOT, LRHS_ROOT,
     &          MTYPE_LOC, 
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE, IW1(IPANEL_POS),
     &          LPANEL_POS, PRUN_BELOW_BWD, TO_PROCESS, SIZE_TO_PROCESS
     &        , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE_BWD
     &         , FROM_PP
     &          )
        CALL DMUMPS_CLEAN_PENDING( INFO(1), KEEP,
     &     BUFR, LBUFR,LBUFR_BYTES,
     &     COMM_NODES, IDUMMY,       
     &     SLAVEF, .TRUE., .FALSE. ) 
        CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
#if defined(V_T)
        CALL VTEND(back_soln,ierr)
#endif
        IF (.NOT.FROM_PP) THEN
          CALL MUMPS_SECFIN(TIME_BWD)
          DKEEP(118)=TIME_BWD+DKEEP(118)
        ENDIF
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
        IF (DOFORWARD) THEN
        K = min0(10,size(RHSCOMP,1))
        IF (LDIAG.EQ.4) K = size(RHSCOMP,1)
        IF ( .NOT. FROM_PP) THEN
          WRITE (MP,99992)
          IF (size(RHSCOMP,1).GT.0) 
     &    WRITE (MP,99993) (RHSCOMP(I,1),I=1,K)
          IF (size(RHSCOMP,1).GT.0.and.NRHS>1) 
     &              WRITE (MP,99994) (RHSCOMP(I,2),I=1,K)
          ENDIF
        ENDIF
      ENDIF
500   CONTINUE
      IF ( allocated(TO_PROCESS)) DEALLOCATE (TO_PROCESS)
      IF (Exploit_Sparsity.OR.AM1.OR.SWITCH_OFF_ES) THEN
         IF ( allocated(nodes_RHS)) DEALLOCATE (nodes_RHS)
         IF ( allocated(Pruned_SONS)) DEALLOCATE (Pruned_SONS)
         IF ( allocated(Pruned_Roots)) DEALLOCATE (Pruned_Roots)
         IF ( allocated(Pruned_List)) DEALLOCATE (Pruned_List)
         IF ( allocated(Pruned_Leaves)) DEALLOCATE (Pruned_Leaves)
      ENDIF
      RETURN 
99993 FORMAT (' RHS (internal, first column)'/(1X,1P,5D14.6))
99994 FORMAT (' RHS (internal, 2 nd  column)'/(1X,1P,5D14.6))
99992 FORMAT (//' LEAVING SOLVE (DMUMPS_SOL_C) WITH')
      END SUBROUTINE DMUMPS_SOL_C
      SUBROUTINE DMUMPS_GATHER_SOLUTION( NSLAVES, N, MYID, COMM,
     &           NRHS,
     &           MTYPE, RHS, LRHS, NCOL_RHS, JBEG_RHS, PTRIST,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, BUFFER,
     &           SIZE_BUF, SIZE_BUF_BYTES, CWORK, LCWORK,
     &           LSCAL, SCALING, LSCALING, 
     &           RHSCOMP, LRHSCOMP, NCOL_RHSCOMP, 
     &           POSINRHSCOMP, LPOS_N, PERM_RHS, SIZE_PERM_RHS )
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM, LIW, MTYPE, NCOL_RHS
      INTEGER NRHS, LRHS, LCWORK, LPOS_N, NCOL_RHSCOMP
      DOUBLE PRECISION RHS   (LRHS, NCOL_RHS)
      INTEGER, INTENT(in) :: JBEG_RHS
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION :: CWORK(LCWORK)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER SIZE_BUF, SIZE_BUF_BYTES
      INTEGER BUFFER(SIZE_BUF)
      INTEGER LRHSCOMP, POSINRHSCOMP(LPOS_N) 
      DOUBLE PRECISION, intent(in) :: RHSCOMP(LRHSCOMP, NCOL_RHSCOMP)
      LOGICAL, intent(in) :: LSCAL
      INTEGER, intent(in) :: LSCALING
      DOUBLE PRECISION, intent(in)    :: SCALING(LSCALING)
      INTEGER, INTENT(in) :: SIZE_PERM_RHS
      INTEGER, INTENT(in) :: PERM_RHS(SIZE_PERM_RHS)
      INTEGER I, II, J, J1, ISTEP, MASTER,
     &        MYID_NODES, TYPE_PARAL, N2RECV
      INTEGER LIELL, IPOS, NPIV, MAXNPIV_estim, MAXSurf
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR, allocok
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER RECORD_SIZE_P_1, SIZE1, SIZE2
      INTEGER POS_BUF, N2SEND, IPOSINRHSCOMP
      INTEGER :: JCOL_RHS
      INTEGER :: K242
!$    LOGICAL :: OMP_FLAG
!$    INTEGER :: CHUNK, NOMP
      INTEGER, PARAMETER :: FIN = -1
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IROWlist(:)
      INCLUDE 'mumps_headers.h'
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      TYPE_PARAL = KEEP(46)  
      I_AM_SLAVE = MYID .ne. MASTER .OR. TYPE_PARAL .eq. 1
      IF ( TYPE_PARAL == 1 ) THEN
        MYID_NODES = MYID
      ELSE
        MYID_NODES = MYID-1
      ENDIF
      IF (NSLAVES.EQ.1 .AND. TYPE_PARAL.EQ.1) THEN
           IF (LSCAL) THEN 
             IF (KEEP(350).EQ.2) THEN
             K242 = KEEP(242)
!$           NOMP = OMP_GET_MAX_THREADS()
!$           OMP_FLAG = .FALSE.
!$           CHUNK = max(N/2,1)
!$           IF (int(NRHS,8) * int(N,8) .GE. int(KEEP(363),8)) THEN
!$             OMP_FLAG = .TRUE.
!$         CHUNK=int((int(N,8)*int(NRHS,8)+int(NOMP-1,8))/int(NOMP,8))
!$             CHUNK = min(CHUNK,(N+KEEP(362)-1)/KEEP(362))
!$             CHUNK = max(KEEP(363)/2,CHUNK)
!$           ENDIF
!$OMP PARALLEL FIRSTPRIVATE(JBEG_RHS,N,K242)
!$OMP&  PRIVATE(IPOSINRHSCOMP,I,JCOL_RHS) IF (OMP_FLAG)
             DO J=1, NRHS
               IF (K242.EQ.0) THEN
                 JCOL_RHS = J+JBEG_RHS-1
               ELSE
                 JCOL_RHS = PERM_RHS(J+JBEG_RHS-1)
               ENDIF
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
               DO I=1, N
                 IPOSINRHSCOMP = POSINRHSCOMP(I)
                 IF (IPOSINRHSCOMP.GT.0) THEN 
                   RHS(I,JCOL_RHS) =
     &                 RHSCOMP(IPOSINRHSCOMP,J)*SCALING(I)
                 ELSE
                   RHS(I,JCOL_RHS) = ZERO
                 ENDIF
               ENDDO
!$OMP END DO NOWAIT
             ENDDO
!$OMP END PARALLEL
             ELSE
             DO J=1, NRHS
               IF (KEEP(242).EQ.0) THEN
                 JCOL_RHS = J+JBEG_RHS-1
               ELSE
                 JCOL_RHS = PERM_RHS(J+JBEG_RHS-1)
               ENDIF
               DO I=1, N
                 IPOSINRHSCOMP = POSINRHSCOMP(I)
                 IF (IPOSINRHSCOMP.GT.0) THEN 
                   RHS(I,JCOL_RHS) =
     &                 RHSCOMP(IPOSINRHSCOMP,J)*SCALING(I)
                 ELSE
                   RHS(I,JCOL_RHS) = ZERO
                 ENDIF
               ENDDO
             ENDDO
             ENDIF
           ELSE
             IF (KEEP(350).EQ.2) THEN
             K242 = KEEP(242)
!$           NOMP = OMP_GET_MAX_THREADS()
!$           OMP_FLAG = .FALSE.
!$           CHUNK = max(N/2,1)
!$           IF (NRHS * N .GE. KEEP(363)) THEN
!$             OMP_FLAG = .TRUE.
!$         CHUNK=int((int(N,8)*int(NRHS,8)+int(NOMP-1,8))/int(NOMP,8))
!$             CHUNK = min(CHUNK,(N+KEEP(362)-1)/KEEP(362))
!$             CHUNK = max(KEEP(363)/2,CHUNK)
!$           ENDIF
!$OMP PARALLEL FIRSTPRIVATE(JBEG_RHS,N,K242)
!$OMP&  PRIVATE(IPOSINRHSCOMP,I,JCOL_RHS) IF (OMP_FLAG)
             DO J=1, NRHS
               IF (K242.EQ.0) THEN
                 JCOL_RHS = J+JBEG_RHS-1
               ELSE
                 JCOL_RHS = PERM_RHS(J+JBEG_RHS-1)
               ENDIF
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
               DO I=1, N
                 IPOSINRHSCOMP = POSINRHSCOMP(I)
                 IF (IPOSINRHSCOMP.GT.0) THEN
                   RHS(I,JCOL_RHS) = RHSCOMP(IPOSINRHSCOMP,J)
                 ELSE
                   RHS(I,JCOL_RHS) = ZERO
                 ENDIF
               ENDDO
!$OMP END DO NOWAIT
             ENDDO
!$OMP END PARALLEL
             ELSE
             DO J=1, NRHS
               IF (KEEP(242).EQ.0) THEN
                 JCOL_RHS = J+JBEG_RHS-1
               ELSE
                 JCOL_RHS = PERM_RHS(J+JBEG_RHS-1)
               ENDIF
               DO I=1, N
                 IPOSINRHSCOMP = POSINRHSCOMP(I)
                 IF (IPOSINRHSCOMP.GT.0) THEN
                   RHS(I,JCOL_RHS) = RHSCOMP(IPOSINRHSCOMP,J)
                 ELSE
                   RHS(I,JCOL_RHS) = ZERO
                 ENDIF
               ENDDO
             ENDDO
             ENDIF
           ENDIF
        RETURN
      ENDIF
      MAXNPIV_estim = max(KEEP(246), KEEP(247))
      MAXSurf       = MAXNPIV_estim*NRHS
      IF (LCWORK .LT. MAXNPIV_estim) THEN
        WRITE(*,*) MYID, 
     &  ": Internal error 2 in DMUMPS_GATHER_SOLUTION:",
     &  TYPE_PARAL, LCWORK, KEEP(247), NRHS
        CALL MUMPS_ABORT()
      ENDIF
      IF (MYID.EQ.MASTER) THEN
         ALLOCATE(IROWlist(KEEP(247)),stat=allocok)
         IF(allocok.GT.0) THEN
            WRITE(*,*)'Problem with allocation of IROWlist'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      IF (NSLAVES .EQ. 1 .AND. TYPE_PARAL .EQ. 1) THEN
         CALL MUMPS_ABORT()
      ENDIF
      SIZE1=0
      CALL MPI_PACK_SIZE(MAXNPIV_estim+2,MPI_INTEGER, COMM, 
     &          SIZE1, IERR)
      SIZE2=0
      CALL MPI_PACK_SIZE(MAXSurf,MPI_DOUBLE_PRECISION, COMM,
     &                   SIZE2, IERR)
      RECORD_SIZE_P_1= SIZE1+SIZE2
      IF (RECORD_SIZE_P_1.GT.SIZE_BUF_BYTES) THEN
         write(6,*) MYID, 
     &    ' Internal error 3 in  DMUMPS_GATHER_SOLUTION '
         write(6,*) MYID, ' RECORD_SIZE_P_1, SIZE_BUF_BYTES=', 
     &                 RECORD_SIZE_P_1, SIZE_BUF_BYTES
         CALL MUMPS_ABORT()
      ENDIF
      N2SEND   =0
      N2RECV   =N
      POS_BUF  =0
      IF (I_AM_SLAVE) THEN
        POS_BUF = 0
        DO ISTEP = 1, KEEP(28)
          IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &          KEEP(199))) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &        NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N)
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              IF (MYID .EQ. MASTER) THEN
                   N2RECV=N2RECV-NPIV
                   IF (NPIV.GT.0) 
     &             CALL DMUMPS_NPIV_BLOCK_ADD ( .TRUE. )
              ELSE
                   IF (NPIV.GT.0) 
     &             CALL DMUMPS_NPIV_BLOCK_ADD ( .FALSE.)
              ENDIF
          ENDIF
        ENDDO
        CALL DMUMPS_NPIV_BLOCK_SEND()   
      ENDIF
      IF ( MYID .EQ. MASTER ) THEN
       DO WHILE (N2RECV .NE. 0)
        CALL MPI_RECV( BUFFER, SIZE_BUF_BYTES, MPI_PACKED,
     &                 MPI_ANY_SOURCE,
     &                 GatherSol, COMM, STATUS, IERR )
        POS_BUF = 0
        CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &                   NPIV, 1, MPI_INTEGER, COMM, IERR)
        DO WHILE (NPIV.NE.FIN)
          CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &             IROWlist, NPIV, MPI_INTEGER, COMM, IERR)
          DO J=1, NRHS
            IF (KEEP(242).EQ.0) THEN
              JCOL_RHS=J+JBEG_RHS-1
            ELSE
              JCOL_RHS=PERM_RHS(J+JBEG_RHS-1)
            ENDIF
            CALL MPI_UNPACK(BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &                 CWORK, NPIV, MPI_DOUBLE_PRECISION,
     &                 COMM, IERR)
            IF (LSCAL) THEN
              DO I=1,NPIV
                RHS(IROWlist(I),JCOL_RHS)=CWORK(I)*SCALING(IROWlist(I))
              ENDDO
            ELSE
              DO I=1,NPIV
                RHS(IROWlist(I),JCOL_RHS)=CWORK(I)
              ENDDO
            ENDIF
          ENDDO
          N2RECV=N2RECV-NPIV
          CALL MPI_UNPACK( BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &                   NPIV, 1, MPI_INTEGER, COMM, IERR)
        ENDDO
       ENDDO
       DEALLOCATE(IROWlist)
      ENDIF
      RETURN
      CONTAINS
        SUBROUTINE DMUMPS_NPIV_BLOCK_ADD ( ON_MASTER )
        LOGICAL, intent(in) ::  ON_MASTER     
        INTEGER :: JPOS, K242
        LOGICAL :: LOCAL_LSCAL
        IF (ON_MASTER) THEN
        IF (KEEP(350).EQ.2
     &  .AND. (NRHS.EQ.1.OR.((NPIV*NRHS*2*KEEP(16)).GE.KEEP(364)))) THEN
          LOCAL_LSCAL = LSCAL
          K242 = KEEP(242)
          DO J=1, NRHS
            IF (K242.EQ.0) THEN
              JPOS = J+JBEG_RHS-1
            ELSE
              JPOS = PERM_RHS(J+JBEG_RHS-1)
            ENDIF
            DO II=1,NPIV
              I=IW(J1+II-1)
              IPOSINRHSCOMP= POSINRHSCOMP(I) 
              IF (LOCAL_LSCAL) THEN
                RHS(I,JPOS) = RHSCOMP(IPOSINRHSCOMP,J)*SCALING(I)
              ELSE
                RHS(I,JPOS) = RHSCOMP(IPOSINRHSCOMP,J)
              ENDIF
            ENDDO
          ENDDO
        ELSE
         IF (KEEP(242).EQ.0) THEN
           IF (LSCAL) THEN
            DO II=1,NPIV
                I=IW(J1+II-1)
                IPOSINRHSCOMP= POSINRHSCOMP(I) 
                DO J=1, NRHS
                  RHS(I,J+JBEG_RHS-1) =
     &               RHSCOMP(IPOSINRHSCOMP,J)*SCALING(I)
                ENDDO
            ENDDO
           ELSE
            DO II=1,NPIV
                I=IW(J1+II-1)
                IPOSINRHSCOMP= POSINRHSCOMP(I)
                DO J=1, NRHS
                  RHS(I,J+JBEG_RHS-1) = RHSCOMP(IPOSINRHSCOMP,J)
                ENDDO
            ENDDO
           ENDIF
         ELSE
           IF (LSCAL) THEN
              DO II=1,NPIV
                I=IW(J1+II-1)
                IPOSINRHSCOMP= POSINRHSCOMP(I)
!DIR$ NOVECTOR
                DO J=1, NRHS
                  RHS(I,PERM_RHS(J+JBEG_RHS-1)) =
     &               RHSCOMP(IPOSINRHSCOMP,J)*SCALING(I)
                ENDDO
              ENDDO
           ELSE
              DO II=1,NPIV
                I=IW(J1+II-1)
                IPOSINRHSCOMP= POSINRHSCOMP(I)
!DIR$ NOVECTOR
                DO J=1, NRHS
                  RHS(I,PERM_RHS(J+JBEG_RHS-1)) =
     &               RHSCOMP(IPOSINRHSCOMP,J)
                ENDDO
              ENDDO
           ENDIF
         ENDIF
        ENDIF
         RETURN
        ENDIF
        CALL MPI_PACK(NPIV, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(IW(J1), NPIV, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        IPOSINRHSCOMP= POSINRHSCOMP(IW(J1)) 
        DO J=1,NRHS
          CALL MPI_PACK(RHSCOMP(IPOSINRHSCOMP,J), NPIV,
     &              MPI_DOUBLE_PRECISION,
     &              BUFFER, SIZE_BUF_BYTES, POS_BUF, COMM,
     &              IERR)
        ENDDO
        N2SEND=N2SEND+NPIV  
        IF ( POS_BUF + RECORD_SIZE_P_1 > SIZE_BUF_BYTES ) THEN
          CALL DMUMPS_NPIV_BLOCK_SEND()
        END IF
        RETURN
        END SUBROUTINE DMUMPS_NPIV_BLOCK_ADD
        SUBROUTINE DMUMPS_NPIV_BLOCK_SEND()
        IF (N2SEND .NE. 0) THEN
         CALL MPI_PACK(FIN, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
         CALL MPI_SEND(BUFFER, POS_BUF, MPI_PACKED, MASTER, 
     &                 GatherSol, COMM, IERR)
        ENDIF
        POS_BUF=0
        N2SEND=0
        RETURN
        END SUBROUTINE DMUMPS_NPIV_BLOCK_SEND
      END SUBROUTINE DMUMPS_GATHER_SOLUTION
      SUBROUTINE DMUMPS_GATHER_SOLUTION_AM1(NSLAVES, N, MYID, COMM,
     &           NRHS, RHSCOMP,  LRHSCOMP, NRHSCOMP_COL,
     &           KEEP, BUFFER,
     &           SIZE_BUF, SIZE_BUF_BYTES, 
     &           LSCAL, SCALING, LSCALING,
     &          IRHS_PTR_COPY, LIRHS_PTR_COPY, 
     &          IRHS_SPARSE_COPY, LIRHS_SPARSE_COPY,
     &          RHS_SPARSE_COPY, LRHS_SPARSE_COPY,
     &          UNS_PERM_INV, LUNS_PERM_INV,
     &          POSINRHSCOMP, LPOS_ROW, NB_FS_IN_RHSCOMP )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM
      INTEGER NRHS, LRHSCOMP, NRHSCOMP_COL 
      DOUBLE PRECISION, intent(in) :: RHSCOMP (LRHSCOMP, NRHSCOMP_COL)
      INTEGER KEEP(500)
      INTEGER SIZE_BUF, SIZE_BUF_BYTES, LPOS_ROW
      INTEGER BUFFER(SIZE_BUF)
      INTEGER, intent(in) :: LIRHS_PTR_COPY, LIRHS_SPARSE_COPY, 
     &                       LRHS_SPARSE_COPY, LUNS_PERM_INV, 
     &                       NB_FS_IN_RHSCOMP
      INTEGER :: IRHS_SPARSE_COPY(LIRHS_SPARSE_COPY), 
     &           IRHS_PTR_COPY(LIRHS_PTR_COPY), 
     &           UNS_PERM_INV(LUNS_PERM_INV), 
     &           POSINRHSCOMP(LPOS_ROW) 
      DOUBLE PRECISION :: RHS_SPARSE_COPY(LRHS_SPARSE_COPY)
      LOGICAL, intent(in) :: LSCAL
      INTEGER, intent(in) :: LSCALING
      DOUBLE PRECISION, intent(in)    :: SCALING(LSCALING)
      INTEGER COLSIZE, K, IZ, IPREV, NBCOL_INBLOC
      INTEGER I, II, J, MASTER,
     &         TYPE_PARAL, N2RECV, IPOSINRHSCOMP
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER RECORD_SIZE_P_1, SIZE1, SIZE2
      INTEGER POS_BUF, N2SEND
      INTEGER, PARAMETER :: FIN = -1
      INCLUDE 'mumps_headers.h'
      TYPE_PARAL = KEEP(46)  
      I_AM_SLAVE = MYID .ne. MASTER .OR. TYPE_PARAL .eq. 1
      NBCOL_INBLOC = size(IRHS_PTR_COPY)-1
      IF (NSLAVES.EQ.1 .AND. TYPE_PARAL.EQ.1) THEN 
        K=1              
        DO J = 1,  NBCOL_INBLOC
           COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
           IF (COLSIZE.EQ.0) CYCLE
           DO IZ=IRHS_PTR_COPY(J), IRHS_PTR_COPY(J+1)-1
             I = IRHS_SPARSE_COPY(IZ)
             IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
             IPOSINRHSCOMP = POSINRHSCOMP(I)
             IF (IPOSINRHSCOMP.GT.0) THEN 
                IF (LSCAL) THEN
                 RHS_SPARSE_COPY(IZ)=
     &                      RHSCOMP(IPOSINRHSCOMP,K)*SCALING(I)
                ELSE
                 RHS_SPARSE_COPY(IZ)=RHSCOMP(IPOSINRHSCOMP,K)
                ENDIF
             ENDIF          
           ENDDO
           K = K + 1
        ENDDO
        RETURN
      ENDIF
      IF (I_AM_SLAVE) THEN
        K=1              
        DO J = 1, NBCOL_INBLOC
           COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
           IF (COLSIZE.EQ.0) CYCLE
           DO IZ=IRHS_PTR_COPY(J), IRHS_PTR_COPY(J+1)-1
             I = IRHS_SPARSE_COPY(IZ)
             IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
             IPOSINRHSCOMP = POSINRHSCOMP(I)
             IF (IPOSINRHSCOMP.GT.0) THEN 
               RHS_SPARSE_COPY(IZ)=RHSCOMP(IPOSINRHSCOMP,K)
             ENDIF          
           ENDDO
           K = K + 1
        ENDDO
      ENDIF
      SIZE1=0
      CALL MPI_PACK_SIZE(3,MPI_INTEGER, COMM,  
     &          SIZE1, IERR)
      SIZE2=0
      CALL MPI_PACK_SIZE(1,MPI_DOUBLE_PRECISION, COMM,
     &                   SIZE2, IERR)
      RECORD_SIZE_P_1= SIZE1+SIZE2
      IF (RECORD_SIZE_P_1.GT.SIZE_BUF_BYTES) THEN
         write(6,*) MYID, 
     &    ' Internal error 3 in  DMUMPS_GATHER_SOLUTION_AM1 '
         write(6,*) MYID, ' RECORD_SIZE_P_1, SIZE_BUF_BYTES=', 
     &                 RECORD_SIZE_P_1, SIZE_BUF_BYTES
         CALL MUMPS_ABORT()
      ENDIF
      N2SEND   =0
      N2RECV   =size(IRHS_SPARSE_COPY)
      POS_BUF  =0
      IF (I_AM_SLAVE) THEN
        DO J = 1,  NBCOL_INBLOC
            COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
            IF (COLSIZE.LE.0) CYCLE
            K = 0 
            DO IZ=IRHS_PTR_COPY(J), IRHS_PTR_COPY(J+1)-1
             I = IRHS_SPARSE_COPY(IZ)
             II = I
             IF  (KEEP(23).NE.0) II = UNS_PERM_INV(I)
             IPOSINRHSCOMP = POSINRHSCOMP(II)
             IF (IPOSINRHSCOMP.GT.0) THEN
               IF (MYID .EQ. MASTER) THEN
                  N2RECV=N2RECV-1
                  IF (LSCAL) 
     &            CALL DMUMPS_AM1_BLOCK_ADD ( .TRUE. )
                  IRHS_SPARSE_COPY( IRHS_PTR_COPY(J) + K) =
     &               I
                  RHS_SPARSE_COPY( IRHS_PTR_COPY(J) + K) =
     &                RHS_SPARSE_COPY(IZ)
                  K = K+1 
               ELSE
                  CALL DMUMPS_AM1_BLOCK_ADD (  .FALSE. )
               ENDIF
              ENDIF          
            ENDDO
            IF (MYID.EQ.MASTER) 
     &             IRHS_PTR_COPY(J) = IRHS_PTR_COPY(J) + K
        ENDDO
        CALL DMUMPS_AM1_BLOCK_SEND()   
      ENDIF
      IF ( MYID .EQ. MASTER ) THEN
       DO WHILE (N2RECV .NE. 0)
        CALL MPI_RECV( BUFFER, SIZE_BUF_BYTES, MPI_PACKED,
     &                 MPI_ANY_SOURCE,
     &                 GatherSol, COMM, STATUS, IERR )
        POS_BUF = 0
        CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &                   J, 1, MPI_INTEGER, COMM, IERR)
        DO WHILE (J.NE.FIN)
          IZ = IRHS_PTR_COPY(J)
          CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &             I, 1, MPI_INTEGER, COMM, IERR)
          IRHS_SPARSE_COPY(IZ) = I
          CALL MPI_UNPACK(BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &             RHS_SPARSE_COPY(IZ), 1, MPI_DOUBLE_PRECISION,
     &             COMM, IERR)
          IF (LSCAL) THEN
              IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
              RHS_SPARSE_COPY(IZ) = RHS_SPARSE_COPY(IZ)*SCALING(I)    
          ENDIF
          N2RECV=N2RECV-1
          IRHS_PTR_COPY(J) = IRHS_PTR_COPY(J) + 1
          CALL MPI_UNPACK( BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &                   J, 1, MPI_INTEGER, COMM, IERR)
        ENDDO
       ENDDO
       IPREV = 1
       DO J=1, size(IRHS_PTR_COPY)-1
         I= IRHS_PTR_COPY(J) 
         IRHS_PTR_COPY(J) = IPREV
         IPREV = I
       ENDDO
      ENDIF
      RETURN
      CONTAINS
        SUBROUTINE DMUMPS_AM1_BLOCK_ADD ( SCALE_ONLY )
        LOGICAL, intent(in) ::  SCALE_ONLY    
        INTEGER III
        IF (SCALE_ONLY) THEN
         III = I
         IF (KEEP(23).NE.0) III = UNS_PERM_INV(I)
         IF (LSCAL) THEN
            RHS_SPARSE_COPY(IZ)=RHS_SPARSE_COPY(IZ)*SCALING(III)
         ENDIF
         RETURN
        ENDIF
        CALL MPI_PACK(J, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(I, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(RHS_SPARSE_COPY(IZ), 1, MPI_DOUBLE_PRECISION,
     &                BUFFER, SIZE_BUF_BYTES, POS_BUF, COMM,
     &                IERR)
        N2SEND=N2SEND+1  
        IF ( POS_BUF + RECORD_SIZE_P_1 > SIZE_BUF_BYTES ) THEN
          CALL DMUMPS_AM1_BLOCK_SEND()
        END IF
        RETURN
        END SUBROUTINE DMUMPS_AM1_BLOCK_ADD
        SUBROUTINE DMUMPS_AM1_BLOCK_SEND()
        IF (N2SEND .NE. 0) THEN
         CALL MPI_PACK(FIN, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
         CALL MPI_SEND(BUFFER, POS_BUF, MPI_PACKED, MASTER, 
     &                 GatherSol, COMM, IERR)
        ENDIF
        POS_BUF=0
        N2SEND=0
        RETURN
        END SUBROUTINE DMUMPS_AM1_BLOCK_SEND
      END SUBROUTINE DMUMPS_GATHER_SOLUTION_AM1
      SUBROUTINE DMUMPS_DISTSOL_INDICES(MTYPE, ISOL_LOC,
     &             PTRIST, KEEP,KEEP8,
     &             IW, LIW_PASSED, MYID_NODES, N, STEP,
     &             PROCNODE, NSLAVES, scaling_data, LSCAL
     &             , IRHS_loc_MEANINGFUL, IRHS_loc, Nloc_RHS
     &             )
      IMPLICIT NONE
      INTEGER MTYPE, MYID_NODES, N, NSLAVES
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE(KEEP(28))
      INTEGER ISOL_LOC(KEEP(89))
      INTEGER LIW_PASSED
      INTEGER IW(LIW_PASSED)
      INTEGER STEP(N)
      LOGICAL LSCAL
      LOGICAL :: IRHS_loc_MEANINGFUL
      INTEGER :: Nloc_RHS
      INTEGER :: IRHS_loc(Nloc_RHS)
      type scaling_data_t
        SEQUENCE
        DOUBLE PRECISION, dimension(:), pointer :: SCALING
        DOUBLE PRECISION, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      type (scaling_data_t) :: scaling_data
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INTEGER ISTEP, K
      INTEGER J1, IPOS, LIELL, NPIV, JJ
      LOGICAL :: CHECK_IRHS_loc
      INTEGER(8) :: DIFF_ADDR
      INCLUDE 'mumps_headers.h'
      CHECK_IRHS_loc=.FALSE.
      IF ( IRHS_loc_MEANINGFUL ) THEN
        IF (Nloc_RHS .GT. 0) THEN
          CALL MUMPS_SIZE_C( IRHS_loc(1), ISOL_loc(1),
     &                       DIFF_ADDR )
          IF (DIFF_ADDR .EQ. 0_8) THEN
            CHECK_IRHS_loc=.TRUE.
          ENDIF
        ENDIF
      ENDIF
      K=0
      DO ISTEP=1, KEEP(28)
          IF ( MYID_NODES == MUMPS_PROCNODE( PROCNODE(ISTEP),
     &                   KEEP(199))) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &        NPIV, LIELL, IPOS, IW, LIW_PASSED, PTRIST, STEP, N)
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              DO JJ=J1,J1+NPIV-1
                  K=K+1
                  IF (CHECK_IRHS_loc) THEN
                    IF (K.LE.Nloc_RHS) THEN
                      IF ( IW(JJ) .NE.IRHS_LOC(K) ) THEN
                      ENDIF
                    ENDIF
                  ENDIF
                  ISOL_LOC(K)=IW(JJ)
                  IF (LSCAL) THEN
                    scaling_data%SCALING_LOC(K)=
     &              scaling_data%SCALING(IW(JJ))
                  ENDIF
              ENDDO
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_DISTSOL_INDICES
      SUBROUTINE DMUMPS_DISTRIBUTED_SOLUTION(
     &           SLAVEF, N, MYID_NODES,
     &           MTYPE, RHSCOMP, LRHSCOMP, NBRHS_EFF,
     &           POSINRHSCOMP,
     &           ISOL_LOC, 
     &           SOL_LOC, NRHS, BEG_RHS, LSOL_LOC,
     &           PTRIST,
     &           PROCNODE_STEPS, KEEP,KEEP8, IW, LIW, STEP,
     &           scaling_data, LSCAL, NB_RHSSKIPPED,
     &           PERM_RHS, SIZE_PERM_RHS)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      type scaling_data_t
        SEQUENCE
        DOUBLE PRECISION, dimension(:), pointer :: SCALING
        DOUBLE PRECISION, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      TYPE (scaling_data_t) :: scaling_data
      LOGICAL LSCAL
      INTEGER SLAVEF, N, MYID_NODES, LIW, MTYPE, NBRHS_EFF, LRHSCOMP
      INTEGER POSINRHSCOMP(N), NB_RHSSKIPPED
      INTEGER LSOL_LOC, BEG_RHS
      INTEGER ISOL_LOC(LSOL_LOC)
      INTEGER, INTENT(in) :: NRHS 
      DOUBLE PRECISION SOL_LOC( LSOL_LOC, NRHS )
      DOUBLE PRECISION RHSCOMP( LRHSCOMP, NBRHS_EFF )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER, INTENT(in) :: SIZE_PERM_RHS
      INTEGER, INTENT(in) :: PERM_RHS( SIZE_PERM_RHS )
      INTEGER :: JJ, J1, ISTEP, K, KLOC, IPOSINRHSCOMP, JEMPTY
      INTEGER :: JCOL, JCOL_PERM
      INTEGER :: IPOS, LIELL, NPIV, JEND
      LOGICAL :: ROOT
!$    LOGICAL :: OMP_FLAG
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      K=0
      JEMPTY = BEG_RHS+NB_RHSSKIPPED-1
      JEND   = BEG_RHS+NB_RHSSKIPPED+NBRHS_EFF-1
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &      KEEP(199))) THEN
          ROOT=.false.
          IF (KEEP(38).ne.0) ROOT = STEP(KEEP(38))==ISTEP
          IF (KEEP(20).ne.0) ROOT = STEP(KEEP(20))==ISTEP
          IF ( ROOT ) THEN
                IPOS = PTRIST(ISTEP) + KEEP(IXSZ)
                LIELL = IW(IPOS+3)
                NPIV = LIELL
                IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
          ELSE
              IPOS = PTRIST(ISTEP) + 2 +KEEP(IXSZ)
              LIELL = IW(IPOS-2)+IW(IPOS+1)
              IPOS= IPOS+1
              NPIV = IW(IPOS)
              IPOS= IPOS+1
              IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +KEEP(IXSZ))
          END IF
          IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
               J1=IPOS+1+LIELL
          ELSE
               J1=IPOS+1
          END IF
            IF (NB_RHSSKIPPED.GT.0) THEN
              DO JCOL = BEG_RHS, JEMPTY
                IF (KEEP(242) .NE. 0) THEN
                JCOL_PERM = PERM_RHS(JCOL)
                ELSE
                 JCOL_PERM = JCOL
                ENDIF
                KLOC=K
                DO JJ=J1,J1+NPIV-1
                  KLOC=KLOC+1
                  SOL_LOC(KLOC, JCOL_PERM) = ZERO
                ENDDO
              ENDDO
            ENDIF
!$            OMP_FLAG = ( JEND-JEMPTY.GE.KEEP(362) .AND.
!$   &                     (NPIV*(JEND-JEMPTY) .GE. KEEP(363)/2 ) )
!$OMP PARALLEL DO PRIVATE(JCOL,JCOL_PERM,KLOC,JJ,IPOSINRHSCOMP) 
!$OMP&         IF(OMP_FLAG)
            DO JCOL = JEMPTY+1, JEND
              IF (KEEP(242) .NE. 0) THEN
                JCOL_PERM = PERM_RHS(JCOL)
              ELSE
                 JCOL_PERM = JCOL
              ENDIF
              DO JJ=J1,J1+NPIV-1
                 KLOC=K + JJ-J1 + 1
                IPOSINRHSCOMP = POSINRHSCOMP(IW(JJ))
                IF (LSCAL) THEN
                  SOL_LOC(KLOC,JCOL_PERM) =
     &            scaling_data%SCALING_LOC(KLOC)*
     &            RHSCOMP(IPOSINRHSCOMP,JCOL-JEMPTY)
                ELSE
                  SOL_LOC(KLOC,JCOL_PERM) =
     &            RHSCOMP(IPOSINRHSCOMP,JCOL-JEMPTY)
                ENDIF
              ENDDO
            ENDDO
!$OMP END PARALLEL DO
          K=K+NPIV
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_DISTRIBUTED_SOLUTION
      SUBROUTINE DMUMPS_SCATTER_RHS
     &           (NSLAVES, N, MYID, COMM,
     &           MTYPE, RHS, LRHS, NCOL_RHS, NRHS,
     &           RHSCOMP, LRHSCOMP, NCOL_RHSCOMP,
     &           POSINRHSCOMP_FWD, NB_FS_IN_RHSCOMP_F,
     &           PTRIST,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, 
     &           ICNTL, INFO)
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM, LIW, MTYPE
      INTEGER NRHS, LRHS, NCOL_RHS,  LRHSCOMP, NCOL_RHSCOMP
      INTEGER ICNTL(60), INFO(80)
      DOUBLE PRECISION, intent(in)  :: RHS (LRHS, NCOL_RHS)
      DOUBLE PRECISION, intent(out) :: RHSCOMP(LRHSCOMP, NCOL_RHSCOMP)
      INTEGER, intent(in)  :: POSINRHSCOMP_FWD(N), NB_FS_IN_RHSCOMP_F
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N) 
      INTEGER BUF_MAXSIZE, BUF_MAXREF
      PARAMETER (BUF_MAXREF=200000)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUF_INDX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BUF_RHS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BUF_RHS_2
      INTEGER ENTRIES_2_PROCESS, PROC_WHO_ASKS, BUF_EFFSIZE
      INTEGER INDX 
      INTEGER allocok
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER I, J, K, JJ, J1, ISTEP, MASTER,
     &        MYID_NODES, TYPE_PARAL
      INTEGER LIELL, IPOS, NPIV
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
!$    INTEGER :: CHUNK, NOMP
!$    LOGICAL :: OMP_FLAG
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      TYPE_PARAL = KEEP(46)
      I_AM_SLAVE = MYID .ne. 0 .OR. TYPE_PARAL .eq. 1
      IF ( TYPE_PARAL == 1 ) THEN
        MYID_NODES = MYID
      ELSE
        MYID_NODES = MYID-1
      ENDIF
      BUF_EFFSIZE = 0
      BUF_MAXSIZE = max(min(BUF_MAXREF,int(2000000/NRHS)),2000)
      IF ( KEEP(350).EQ.2 ) THEN
!$       NOMP = OMP_GET_MAX_THREADS()
         ALLOCATE (BUF_INDX(BUF_MAXSIZE),
     &        BUF_RHS_2(BUF_MAXSIZE*NRHS),
     &        stat=allocok)
      ELSE
         ALLOCATE (BUF_INDX(BUF_MAXSIZE),
     &        BUF_RHS(NRHS,BUF_MAXSIZE),
     &        stat=allocok)
      END IF
      IF (allocok .GT. 0) THEN
        INFO(1)=-13
        INFO(2)=BUF_MAXSIZE*(NRHS+1)
      ENDIF
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM, MYID )
      IF (INFO(1).LT.0) RETURN
      IF (MYID.EQ.MASTER) THEN
        ENTRIES_2_PROCESS = N - KEEP(89)
        IF (TYPE_PARAL.EQ.1.AND.ENTRIES_2_PROCESS.NE.0) THEN
         IF (NB_FS_IN_RHSCOMP_F.LT.LRHSCOMP) THEN
           DO K=1, NCOL_RHSCOMP
            DO I = NB_FS_IN_RHSCOMP_F +1, LRHSCOMP
             RHSCOMP (I, K) = ZERO
            ENDDO
           ENDDO
         ENDIF
        ENDIF
        IF ( KEEP(350).EQ.2 ) THEN
           DO WHILE ( ENTRIES_2_PROCESS .NE. 0)
              CALL MPI_RECV( BUF_INDX, BUF_MAXSIZE, MPI_INTEGER,
     &             MPI_ANY_SOURCE,
     &             ScatterRhsI, COMM, STATUS, IERR )
              CALL MPI_GET_COUNT(STATUS,MPI_INTEGER,BUF_EFFSIZE,IERR)
              PROC_WHO_ASKS = STATUS(MPI_SOURCE)
!$            OMP_FLAG = .FALSE.
!$            CHUNK = NRHS
!$            IF (BUF_EFFSIZE*NRHS .GE. KEEP(363)) THEN
!$              OMP_FLAG = .TRUE.
!$              CHUNK = max((BUF_EFFSIZE*NRHS+NOMP-1)/NOMP,KEEP(363)/2)
!$            ENDIF
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK) PRIVATE(INDX)
!$OMP&  IF (OMP_FLAG)
              DO K = 1, NRHS
                 DO I = 1, BUF_EFFSIZE
                    INDX = BUF_INDX( I )
                    BUF_RHS_2( I+(K-1)*BUF_EFFSIZE) = RHS( INDX, K )
                 ENDDO
              ENDDO
!$OMP END PARALLEL DO
              CALL MPI_SEND( BUF_RHS_2,
     &             NRHS*BUF_EFFSIZE,
     &             MPI_DOUBLE_PRECISION, PROC_WHO_ASKS,
     &             ScatterRhsR, COMM, IERR)
              ENTRIES_2_PROCESS = ENTRIES_2_PROCESS - BUF_EFFSIZE
           ENDDO
           BUF_EFFSIZE= 0
        ELSE
            DO WHILE ( ENTRIES_2_PROCESS .NE. 0)
              CALL MPI_RECV( BUF_INDX, BUF_MAXSIZE, MPI_INTEGER,
     &             MPI_ANY_SOURCE,
     &             ScatterRhsI, COMM, STATUS, IERR )
              CALL MPI_GET_COUNT( STATUS, MPI_INTEGER,BUF_EFFSIZE,IERR)
              PROC_WHO_ASKS = STATUS(MPI_SOURCE)
              DO I = 1, BUF_EFFSIZE
                 INDX = BUF_INDX( I )
                 DO K = 1, NRHS
                    BUF_RHS( K, I ) = RHS( INDX, K )
                 ENDDO
              ENDDO
              CALL MPI_SEND( BUF_RHS, NRHS*BUF_EFFSIZE,
     &             MPI_DOUBLE_PRECISION, PROC_WHO_ASKS,
     &             ScatterRhsR, COMM, IERR)
              ENTRIES_2_PROCESS = ENTRIES_2_PROCESS - BUF_EFFSIZE
           ENDDO
           BUF_EFFSIZE= 0
        ENDIF
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (MYID.NE.MASTER) THEN
         IF (NB_FS_IN_RHSCOMP_F.LT.LRHSCOMP) THEN
           DO K=1, NCOL_RHSCOMP
            DO I = NB_FS_IN_RHSCOMP_F +1, LRHSCOMP
             RHSCOMP (I, K) = ZERO
            ENDDO
           ENDDO
         ENDIF
        ENDIF
        DO ISTEP = 1, KEEP(28)
          IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &          KEEP(199))) THEN
             CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &       NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
              ELSE
                   J1=IPOS+1+LIELL
              END IF
              IF (MYID.EQ.MASTER) THEN
                INDX = POSINRHSCOMP_FWD(IW(J1))
                IF (KEEP(350).EQ.2 .AND.
     &   (NRHS.EQ.1.OR.((NPIV*NRHS*2*KEEP(16)).GE.KEEP(364)))) THEN
!$                 OMP_FLAG = .FALSE.
!$                 CHUNK = NRHS
!$                 IF (NPIV*NRHS .GE. KEEP(363)) THEN
!$                   OMP_FLAG = .TRUE.
!$                   CHUNK = max((NPIV*NRHS+NOMP-1)/NOMP,KEEP(363)/2)
!$                 ENDIF
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK) PRIVATE(J,JJ)
!$OMP&  FIRSTPRIVATE(INDX) IF (OMP_FLAG)
                   DO K = 1, NRHS
                      DO JJ=J1,J1+NPIV-1
                         J=IW(JJ)
                         RHSCOMP( INDX+JJ-J1, K ) = RHS( J, K )
                      ENDDO
                   ENDDO
!$OMP END PARALLEL DO
                ELSE 
                   DO JJ=J1,J1+NPIV-1
                      J=IW(JJ)
                      DO K = 1, NRHS
                         RHSCOMP( INDX+JJ-J1, K ) = RHS( J, K )
                      ENDDO
                   ENDDO
                END IF 
              ELSE
                DO JJ=J1,J1+NPIV-1
                  BUF_EFFSIZE = BUF_EFFSIZE + 1
                  BUF_INDX(BUF_EFFSIZE) = IW(JJ)
                  IF (BUF_EFFSIZE + 1 .GT. BUF_MAXSIZE) THEN
                   CALL DMUMPS_GET_BUF_INDX_RHS()
                  ENDIF
                ENDDO
              ENDIF
          ENDIF
        ENDDO
        IF ( BUF_EFFSIZE .NE. 0 .AND. MYID.NE.MASTER ) 
     &              CALL DMUMPS_GET_BUF_INDX_RHS()
      ENDIF
      IF (KEEP(350).EQ.2) THEN
        DEALLOCATE (BUF_INDX, BUF_RHS_2)
      ELSE
        DEALLOCATE (BUF_INDX, BUF_RHS)
      ENDIF
      RETURN
      CONTAINS
                  SUBROUTINE DMUMPS_GET_BUF_INDX_RHS()
                  CALL MPI_SEND(BUF_INDX, BUF_EFFSIZE, MPI_INTEGER,
     &            MASTER, ScatterRhsI, COMM, IERR )
                  IF (KEEP(350).EQ.2) THEN
                    CALL MPI_RECV(BUF_RHS_2, BUF_EFFSIZE*NRHS,
     &                 MPI_DOUBLE_PRECISION,
     &                 MASTER,
     &                 ScatterRhsR, COMM, STATUS, IERR )
!$                  OMP_FLAG = .FALSE.
!$                  CHUNK = NRHS
!$                  IF (BUF_EFFSIZE*NRHS .GE. KEEP(363)) THEN
!$                    OMP_FLAG = .TRUE.
!$              CHUNK = max((BUF_EFFSIZE*NRHS+NOMP-1)/NOMP,KEEP(363)/2)
!$                  ENDIF
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK) PRIVATE(INDX)
!$OMP&  IF (OMP_FLAG)
                    DO K = 1, NRHS
                        DO I = 1, BUF_EFFSIZE
                           INDX = POSINRHSCOMP_FWD(BUF_INDX(I))
                           RHSCOMP( INDX, K ) =
     &                              BUF_RHS_2( I+(K-1)*BUF_EFFSIZE)
                        ENDDO
                     ENDDO
!$OMP END PARALLEL DO
                   ELSE
                     CALL MPI_RECV(BUF_RHS, BUF_EFFSIZE*NRHS,
     &                 MPI_DOUBLE_PRECISION,
     &                 MASTER,
     &                 ScatterRhsR, COMM, STATUS, IERR )
                     DO I = 1, BUF_EFFSIZE
                        INDX = POSINRHSCOMP_FWD(BUF_INDX(I))
                        DO K = 1, NRHS
                           RHSCOMP( INDX, K ) = BUF_RHS( K, I )
                        ENDDO
                     ENDDO
                  END IF 
                  BUF_EFFSIZE = 0
                  RETURN
                  END SUBROUTINE DMUMPS_GET_BUF_INDX_RHS
      END SUBROUTINE DMUMPS_SCATTER_RHS
      SUBROUTINE DMUMPS_BUILD_POSINRHSCOMP
     &           (NSLAVES, N, MYID_NODES,
     &           PTRIST,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, 
     &           POSINRHSCOMP_ROW, POSINRHSCOMP_COL,
     &           POSINRHSCOMP_COL_ALLOC,
     &           MTYPE,
     &           NBENT_RHSCOMP, NB_FS_IN_RHSCOMP )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER, intent(in) :: NSLAVES, N, MYID_NODES, LIW
      INTEGER, intent(in) :: KEEP(500)
      INTEGER(8), intent(in) :: KEEP8(150)
      INTEGER, intent(in) :: PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER, intent(in) :: IW(LIW), STEP(N)
      INTEGER, intent(in) :: MTYPE 
      LOGICAL, intent(in) :: POSINRHSCOMP_COL_ALLOC
      INTEGER, intent(out):: POSINRHSCOMP_ROW(N),POSINRHSCOMP_COL(N) 
      INTEGER, intent(out):: NBENT_RHSCOMP, NB_FS_IN_RHSCOMP
      INTEGER ISTEP
      INTEGER NPIV
      INTEGER IPOS, LIELL
      INTEGER JJ, J1, JCOL
      INTEGER IPOSINRHSCOMP, IPOSINRHSCOMP_COL
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      POSINRHSCOMP_ROW = 0
      IF (POSINRHSCOMP_COL_ALLOC) POSINRHSCOMP_COL = 0
      IPOSINRHSCOMP   = 1     
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &     KEEP(199))) THEN
           CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP, NPIV, LIELL,
     &     IPOS, IW, LIW, PTRIST, STEP, N )
           IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
           ELSE
                   J1=IPOS+1+LIELL
           END IF
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
                   JCOL = IPOS+1+LIELL
           ELSE
                   JCOL = IPOS+1
           ENDIF
           DO JJ = J1, J1+NPIV-1
                POSINRHSCOMP_ROW(IW(JJ)) = IPOSINRHSCOMP+JJ-J1
           ENDDO
           IF (POSINRHSCOMP_COL_ALLOC) THEN
             DO JJ = JCOL, JCOL+NPIV-1
               POSINRHSCOMP_COL(IW(JJ)) = IPOSINRHSCOMP+JJ-JCOL
             ENDDO
           ENDIF
           IPOSINRHSCOMP       = IPOSINRHSCOMP + NPIV
        ENDIF
      ENDDO
      NB_FS_IN_RHSCOMP = IPOSINRHSCOMP -1
      IF (POSINRHSCOMP_COL_ALLOC) IPOSINRHSCOMP_COL=IPOSINRHSCOMP
      IF (IPOSINRHSCOMP.GT.N) GOTO 500 
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &     KEEP(199))) THEN
           CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &     NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
           IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
           ELSE
                   J1=IPOS+1+LIELL
           END IF
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
                   JCOL = IPOS+1+LIELL
           ELSE
                   JCOL = IPOS+1
           ENDIF
           IF (POSINRHSCOMP_COL_ALLOC) THEN
            DO JJ = NPIV, LIELL-1-KEEP(253)
              IF (POSINRHSCOMP_ROW(IW(J1+JJ)).EQ.0) THEN
               POSINRHSCOMP_ROW(IW(J1+JJ)) = - IPOSINRHSCOMP
               IPOSINRHSCOMP = IPOSINRHSCOMP + 1
              ENDIF
              IF (POSINRHSCOMP_COL(IW(JCOL+JJ)).EQ.0) THEN
               POSINRHSCOMP_COL(IW(JCOL+JJ)) = - IPOSINRHSCOMP_COL
               IPOSINRHSCOMP_COL = IPOSINRHSCOMP_COL + 1
              ENDIF
             ENDDO
           ELSE
             DO JJ = J1+NPIV, J1+LIELL-1-KEEP(253)
              IF (POSINRHSCOMP_ROW(IW(JJ)).EQ.0) THEN
               POSINRHSCOMP_ROW(IW(JJ)) = - IPOSINRHSCOMP
               IPOSINRHSCOMP = IPOSINRHSCOMP + 1
              ENDIF
             ENDDO
           ENDIF
        ENDIF
      ENDDO
 500  NBENT_RHSCOMP = IPOSINRHSCOMP - 1
      IF (POSINRHSCOMP_COL_ALLOC) 
     &     NBENT_RHSCOMP = max(NBENT_RHSCOMP, IPOSINRHSCOMP_COL-1)
      RETURN
      END SUBROUTINE DMUMPS_BUILD_POSINRHSCOMP
      SUBROUTINE DMUMPS_BUILD_POSINRHSCOMP_AM1
     &           (NSLAVES, N, MYID_NODES,
     &           PTRIST, DAD,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP,
     &           POSINRHSCOMP_ROW, POSINRHSCOMP_COL,
     &           POSINRHSCOMP_COL_ALLOC,
     &           MTYPE,
     &           IRHS_PTR, NBCOL_INBLOC, IRHS_SPARSE, NZ_RHS,
     &           PERM_RHS, SIZE_PERM_RHS, JBEG_RHS,
     &           NBENT_RHSCOMP,
     &           NB_FS_IN_RHSCOMP_FWD, NB_FS_IN_RHSCOMP_TOT,
     &           UNS_PERM_INV, SIZE_UNS_PERM_INV
     &            )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER, intent(in) :: NSLAVES, N, MYID_NODES, LIW, 
     &                       SIZE_UNS_PERM_INV
      INTEGER, intent(in) :: KEEP(500)
      INTEGER(8), intent(in) :: KEEP8(150)
      INTEGER, intent(in) :: PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER, intent(inout) :: DAD(KEEP(28))
      INTEGER, intent(in) :: IW(LIW), STEP(N)
      INTEGER, intent(in) :: NBCOL_INBLOC, IRHS_PTR(NBCOL_INBLOC+1)
      INTEGER, intent(in) :: NZ_RHS, IRHS_SPARSE(NZ_RHS)
      INTEGER, intent(in) :: SIZE_PERM_RHS, PERM_RHS(SIZE_PERM_RHS)
      INTEGER, intent(in) :: JBEG_RHS
      INTEGER, intent(in) :: MTYPE
      LOGICAL, intent(in) :: POSINRHSCOMP_COL_ALLOC
      INTEGER, intent(out):: POSINRHSCOMP_ROW(N),POSINRHSCOMP_COL(N)
      INTEGER, intent(out):: NBENT_RHSCOMP
      INTEGER, intent(out):: NB_FS_IN_RHSCOMP_FWD, NB_FS_IN_RHSCOMP_TOT
      INTEGER, intent(in) :: UNS_PERM_INV(SIZE_UNS_PERM_INV)
      INTEGER I, JAM1
      INTEGER ISTEP, OLDISTEP
      INTEGER NPIV
      INTEGER IPOS, LIELL
      INTEGER JJ, J1, JCOL, ABSJCOL
      INTEGER IPOSINRHSCOMP_ROW, IPOSINRHSCOMP_COL
      INTEGER NBENT_RHSCOMP_ROW, NBENT_RHSCOMP_COL
      LOGICAL GO_UP
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      IF(KEEP(237).EQ.0) THEN
        WRITE(*,*)'BUILD_POSINRHSCOMP_SPARSE available for A-1 only !'
        CALL MUMPS_ABORT()
      END IF
      POSINRHSCOMP_ROW = 0
      IF (POSINRHSCOMP_COL_ALLOC) POSINRHSCOMP_COL = 0
      IPOSINRHSCOMP_ROW = 0
      IPOSINRHSCOMP_COL = 0
      DO I = 1, NBCOL_INBLOC
        IF ((IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE 
        IF (KEEP(242).NE.0) THEN 
          JAM1 = PERM_RHS(JBEG_RHS+I-1)
        ELSE
          JAM1 = JBEG_RHS+I-1
        END IF
        ISTEP = abs(STEP(JAM1))
        GO_UP = .TRUE.
        DO WHILE(GO_UP) 
          IF(MYID_NODES.EQ.
     &      MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))) THEN
            CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &      NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
            IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
              J1=IPOS+1
            ELSE
              J1=IPOS+1+LIELL
            END IF
            IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
              JCOL = IPOS+1+LIELL
            ELSE
              JCOL = IPOS+1
            ENDIF
            IF(NPIV.GT.0) THEN 
              IF(POSINRHSCOMP_ROW(IW(J1)).EQ.0) THEN 
                DO JJ = J1, J1+NPIV-1
                  POSINRHSCOMP_ROW(IW(JJ))
     &            = IPOSINRHSCOMP_ROW + JJ - J1 + 1
                ENDDO
                IPOSINRHSCOMP_ROW = IPOSINRHSCOMP_ROW + NPIV
                IF (POSINRHSCOMP_COL_ALLOC) THEN
                  DO JJ = JCOL, JCOL+NPIV-1
                    POSINRHSCOMP_COL(IW(JJ))
     &              = - N - (IPOSINRHSCOMP_COL + JJ - JCOL + 1)
                  ENDDO
                  IPOSINRHSCOMP_COL = IPOSINRHSCOMP_COL + NPIV
                ENDIF
              ELSE
                GO_UP = .FALSE. 
              END IF
            END IF 
          END IF 
          IF(DAD(ISTEP).NE.0) THEN
            ISTEP = STEP(DAD(ISTEP))
          ELSE 
            GO_UP = .FALSE.
          END IF
        END DO 
      END DO 
      NB_FS_IN_RHSCOMP_FWD = IPOSINRHSCOMP_ROW
      IF(POSINRHSCOMP_COL_ALLOC) THEN 
        DO I =1, NZ_RHS
          JAM1 = IRHS_SPARSE(I)
          IF (KEEP(23).NE.0) JAM1 = UNS_PERM_INV(JAM1)
          ISTEP = abs(STEP(JAM1))
          GO_UP = .TRUE.
          DO WHILE(GO_UP) 
            IF(MYID_NODES.EQ.
     &        MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &        NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                J1=IPOS+1
              ELSE
                J1=IPOS+1+LIELL
              END IF
              IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
                JCOL = IPOS+1+LIELL
              ELSE
                JCOL = IPOS+1
              ENDIF
              ABSJCOL = abs(IW(JCOL))
              IF(NPIV.GT.0) THEN 
                IF(POSINRHSCOMP_COL(ABSJCOL).EQ.0)  THEN
                  DO JJ = JCOL, JCOL+NPIV-1
                    POSINRHSCOMP_COL(abs(IW(JJ))) = 
     &                     IPOSINRHSCOMP_COL+JJ-JCOL+1
                  END DO
                  IPOSINRHSCOMP_COL = IPOSINRHSCOMP_COL + NPIV
                ELSE IF (POSINRHSCOMP_COL(ABSJCOL).LT.-N) THEN
                  DO JJ = JCOL, JCOL+NPIV-1
                    POSINRHSCOMP_COL(abs(IW(JJ)))=
     &                -(N+POSINRHSCOMP_COL(abs(IW(JJ))))
                  END DO
                ELSE IF ((POSINRHSCOMP_COL(ABSJCOL).LT.0).AND.
     &                 (POSINRHSCOMP_COL(ABSJCOL).GE.-N))THEN
                  WRITE(*,*)'Internal error 7 in BUILD...SPARSE'
                  CALL MUMPS_ABORT()
                ELSE
                  GO_UP = .FALSE.
                END IF
              END IF 
            END IF 
            IF(DAD(ISTEP).NE.0) THEN
              ISTEP = STEP(DAD(ISTEP))
            ELSE 
              GO_UP = .FALSE.
            END IF
          END DO 
        END DO 
      END IF 
      NB_FS_IN_RHSCOMP_TOT = IPOSINRHSCOMP_COL
            IF (NSLAVES.NE.1) THEN
      DO I = 1, NBCOL_INBLOC
        IF ((IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE 
        IF (KEEP(242).NE.0) THEN 
          JAM1 = PERM_RHS(JBEG_RHS+I-1)
        ELSE
          JAM1 = JBEG_RHS+I-1
        END IF
        ISTEP = abs(STEP(JAM1))
        GO_UP = .TRUE.
        DO WHILE(GO_UP) 
          IF(MYID_NODES.EQ.
     &      MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))) THEN
            CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &      NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
            IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
              J1=IPOS+1
            ELSE
              J1=IPOS+1+LIELL
            END IF
            IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
              JCOL = IPOS+1+LIELL
            ELSE
              JCOL = IPOS+1
            ENDIF
            DO JJ = NPIV, LIELL-1-KEEP(253)
              IF(POSINRHSCOMP_ROW(IW(J1+JJ)).EQ.0) THEN 
                IPOSINRHSCOMP_ROW = IPOSINRHSCOMP_ROW + 1
                POSINRHSCOMP_ROW(IW(JJ+J1))
     &          = -IPOSINRHSCOMP_ROW
              END IF
            END DO
          END IF
          IF(DAD(ISTEP).GT.0) THEN
            OLDISTEP=ISTEP
            ISTEP = STEP(DAD(ISTEP))
            DAD(OLDISTEP)=-DAD(OLDISTEP)
          ELSE 
            GO_UP = .FALSE.
          END IF
        END DO 
      END DO 
      DAD=ABS(DAD)
      IF(POSINRHSCOMP_COL_ALLOC) THEN
        DO I =1, NZ_RHS
          JAM1 = IRHS_SPARSE(I)
          IF (KEEP(23).NE.0) JAM1 = UNS_PERM_INV(JAM1)
          ISTEP = abs(STEP(JAM1))
          GO_UP = .TRUE.
          DO WHILE(GO_UP) 
            IF(MYID_NODES.EQ.
     &        MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &        NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                J1=IPOS+1
              ELSE
                J1=IPOS+1+LIELL
              END IF
              IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
                JCOL = IPOS+1+LIELL
              ELSE
                JCOL = IPOS+1
              ENDIF
              IF (KEEP(23).NE.0) JAM1 = UNS_PERM_INV(JAM1)
              DO JJ = NPIV, LIELL-1-KEEP(253)
                IF(POSINRHSCOMP_COL(IW(JCOL+JJ)).EQ.0) THEN 
                  IPOSINRHSCOMP_COL = IPOSINRHSCOMP_COL + 1
                  POSINRHSCOMP_COL(IW(JCOL+JJ))
     &            = -IPOSINRHSCOMP_COL
                ELSE IF (POSINRHSCOMP_COL(IW(JCOL+JJ)).LT.-N) THEN
                  IPOSINRHSCOMP_COL = IPOSINRHSCOMP_COL + 1
                  POSINRHSCOMP_COL(IW(JCOL+JJ))
     &            = POSINRHSCOMP_COL(IW(JCOL+JJ)) + N              
                END IF
              END DO
            END IF
            IF(DAD(ISTEP).GT.0) THEN
              OLDISTEP=ISTEP
              ISTEP = STEP(DAD(ISTEP))
              DAD(OLDISTEP)=-DAD(OLDISTEP)
            ELSE 
              GO_UP = .FALSE.
            END IF
          END DO 
        END DO 
       DAD=ABS(DAD)
      END IF
      ENDIF
      NBENT_RHSCOMP_ROW = IPOSINRHSCOMP_ROW
      NBENT_RHSCOMP_COL = IPOSINRHSCOMP_COL
      NBENT_RHSCOMP = max(NBENT_RHSCOMP_ROW,NBENT_RHSCOMP_COL)
      RETURN
      END SUBROUTINE DMUMPS_BUILD_POSINRHSCOMP_AM1
