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
      SUBROUTINE DMUMPS_END_DRIVER( id )
      USE DMUMPS_OOC
      USE DMUMPS_STRUC_DEF
      USE DMUMPS_BUF
      IMPLICIT NONE
      include 'mpif.h'
      TYPE( DMUMPS_STRUC ) :: id
      LOGICAL I_AM_SLAVE
      INTEGER IERR
      INTEGER MASTER
      PARAMETER ( MASTER = 0 )
C     Explicit needed because of pointer arguments
      INTERFACE
      SUBROUTINE DMUMPS_FREE_ID_DATA_MODULES(id_FDM_F_ENCODING,
     &  id_BLRARRAY_ENCODING, KEEP8)
#     if defined(MUMPS_F2003)
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_FDM_F_ENCODING
#     else
      CHARACTER, DIMENSION(:), POINTER :: id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER :: id_FDM_F_ENCODING
#     endif
      INTEGER(8), intent(inout) :: KEEP8(150)
      END SUBROUTINE DMUMPS_FREE_ID_DATA_MODULES
      END INTERFACE
      I_AM_SLAVE = ( id%MYID .ne. MASTER .OR. id%KEEP(46) .NE. 0 )
C     ----------------------------------
C     Special stuff for implementations
C     where MPI_CANCEL does not exist or
C     is not correctly implemented.
C     At the moment, this is only
C     required for the slaves.
C     ----------------------------------
      IF (id%KEEP(201).GT.0 .AND. I_AM_SLAVE) THEN
        CALL DMUMPS_CLEAN_OOC_DATA(id,IERR)
        IF (IERR < 0) THEN
          id%INFO(1) = -90
          id%INFO(2) = 0
        ENDIF
      END IF
      CALL MUMPS_PROPINFO(id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID)
      IF (id%root%gridinit_done) THEN
        IF ( id%KEEP(38).NE.0 .and. id%root%yes ) THEN
          CALL blacs_gridexit( id%root%CNTXT_BLACS )
          id%root%gridinit_done = .FALSE.
        END IF
      END IF
      IF ( id%MYID .NE. MASTER .OR. id%KEEP(46) .ne. 0 ) THEN
C       Note that on some old platforms, COMM_NODES would have been
C       freed inside BLACS_GRIDEXIT, which may cause problems
C       in the call to MPI_COMM_FREE. (This was the case on the
C       old SP2 in Bonn.)
        CALL MPI_COMM_FREE( id%COMM_NODES, IERR )
C       Free communicator related to load messages.
        CALL MPI_COMM_FREE( id%COMM_LOAD, IERR )
      END IF
C     -----------------------------------
C     Right-hand-side is always user data
C     We do not free it.
C     -----------------------------------
      IF (associated(id%MEM_DIST))  THEN
         DEALLOCATE(id%MEM_DIST)
         NULLIFY(id%MEM_DIST)
      ENDIF
C
C
C
C  ---------------------------------
C  Allocated by DMUMPS, Used by user.
C  DMUMPS deallocates. User should
C  use them before DMUMPS_END_DRIVER or
C  copy.
C  ---------------------------------
      IF (associated(id%MAPPING)) THEN
        DEALLOCATE(id%MAPPING)
        NULLIFY(id%MAPPING)
      END IF
      NULLIFY(id%SCHUR_CINTERFACE)
C
C     -------------------------------------
C     Always deallocate scaling arrays
C     if they are associated, except
C     when provided by the user (on master)
C     -------------------------------------
      IF ( id%KEEP(52) .NE. -1 .or. id%MYID .ne. MASTER ) THEN
        IF (associated(id%COLSCA)) THEN
          DEALLOCATE(id%COLSCA)
          NULLIFY(id%COLSCA)
        ENDIF
        IF (associated(id%ROWSCA)) THEN
          DEALLOCATE(id%ROWSCA)
          NULLIFY(id%ROWSCA)
        ENDIF
      END IF
      IF (associated(id%PTLUST_S)) THEN
        DEALLOCATE(id%PTLUST_S)
        NULLIFY(id%PTLUST_S)
      END IF
      IF (associated(id%PTRFAC)) THEN
        DEALLOCATE(id%PTRFAC)
        NULLIFY(id%PTRFAC)
      END IF
      IF (associated(id%IS)) THEN
        DEALLOCATE(id%IS)
        NULLIFY(id%IS)
      ENDIF
      IF (associated(id%STEP))      THEN
        DEALLOCATE(id%STEP)
        NULLIFY(id%STEP)
      ENDIF
C     Begin PRUN_NODES
C     Info for pruning tree 
      IF (associated(id%Step2node))      THEN
        DEALLOCATE(id%Step2node)
        NULLIFY(id%Step2node)
      ENDIF
C     END PRUN_NODES
c     --------------------- 
      IF (associated(id%NE_STEPS))  THEN
        DEALLOCATE(id%NE_STEPS)
        NULLIFY(id%NE_STEPS)
      ENDIF
      IF (associated(id%ND_STEPS))  THEN
        DEALLOCATE(id%ND_STEPS)
        NULLIFY(id%ND_STEPS)
      ENDIF
      IF (associated(id%FRERE_STEPS))  THEN
        DEALLOCATE(id%FRERE_STEPS)
        NULLIFY(id%FRERE_STEPS)
      ENDIF
      IF (associated(id%DAD_STEPS))  THEN
        DEALLOCATE(id%DAD_STEPS)
        NULLIFY(id%DAD_STEPS)
      ENDIF
      IF (associated(id%SYM_PERM))  THEN
        DEALLOCATE(id%SYM_PERM)
        NULLIFY(id%SYM_PERM)
      ENDIF
      IF (associated(id%UNS_PERM))  THEN
        DEALLOCATE(id%UNS_PERM)
        NULLIFY(id%UNS_PERM)
      ENDIF
      IF (associated(id%PIVNUL_LIST))  THEN
        DEALLOCATE(id%PIVNUL_LIST)
        NULLIFY(id%PIVNUL_LIST)
      ENDIF
      IF (associated(id%FILS))      THEN
        DEALLOCATE(id%FILS)
        NULLIFY(id%FILS)
      ENDIF
      IF (associated(id%PTRAR))     THEN
        DEALLOCATE(id%PTRAR)
        NULLIFY(id%PTRAR)
      ENDIF
      IF (associated(id%FRTPTR))    THEN
        DEALLOCATE(id%FRTPTR)
        NULLIFY(id%FRTPTR)
      ENDIF
      IF (associated(id%FRTELT))    THEN
        DEALLOCATE(id%FRTELT)
        NULLIFY(id%FRTELT)
      ENDIF
      IF (associated(id%NA))        THEN
        DEALLOCATE(id%NA)
        NULLIFY(id%NA)
      ENDIF
      IF (associated(id%PROCNODE_STEPS)) THEN
        DEALLOCATE(id%PROCNODE_STEPS)
        NULLIFY(id%PROCNODE_STEPS)
      ENDIF
      IF (associated(id%RHSCOMP)) THEN
        DEALLOCATE(id%RHSCOMP)
        NULLIFY(id%RHSCOMP)
        id%KEEP8(25)=0_8
      ENDIF
      IF (associated(id%POSINRHSCOMP_ROW)) THEN
        DEALLOCATE(id%POSINRHSCOMP_ROW)
        NULLIFY(id%POSINRHSCOMP_ROW)
      ENDIF
      IF (id%POSINRHSCOMP_COL_ALLOC) THEN
        DEALLOCATE(id%POSINRHSCOMP_COL)
        NULLIFY(id%POSINRHSCOMP_COL)
        id%POSINRHSCOMP_COL_ALLOC = .FALSE.
      ENDIF
C     ------------------------------------------------
C     For hybrid host and element entry,
C     and DBLARR have not been allocated
C     on the master except if there was scaing.
C     ------------------------------------------------
      IF (id%KEEP(46).eq.1 .and.
     &    id%KEEP(55).ne.0 .and.
     &    id%MYID .eq. MASTER .and.
     &    id%KEEP(52) .eq. 0 ) THEN
        NULLIFY(id%DBLARR)
      ELSE
        IF (associated(id%DBLARR)) THEN
          DEALLOCATE(id%DBLARR)
          NULLIFY(id%DBLARR)
        ENDIF
      END IF
      IF (associated(id%INTARR))       THEN
        DEALLOCATE(id%INTARR)
        NULLIFY(id%INTARR)
      ENDIF
      IF (associated(id%root%RG2L_ROW))THEN
        DEALLOCATE(id%root%RG2L_ROW)
        NULLIFY(id%root%RG2L_ROW)
      ENDIF
      IF (associated(id%root%RG2L_COL))THEN
        DEALLOCATE(id%root%RG2L_COL)
        NULLIFY(id%root%RG2L_COL)
      ENDIF
C     IPIV is used both for ScaLAPACK and RR
C     Keep it outside DMUMPS_RR_FREE_POINTERS
      IF (associated(id%root%IPIV))    THEN
        DEALLOCATE(id%root%IPIV)
        NULLIFY(id%root%IPIV)
      ENDIF
      IF (associated(id%root%RHS_CNTR_MASTER_ROOT)) THEN
        DEALLOCATE(id%root%RHS_CNTR_MASTER_ROOT)
        NULLIFY(id%root%RHS_CNTR_MASTER_ROOT)
      ENDIF
      IF (associated(id%root%RHS_ROOT))THEN
        DEALLOCATE(id%root%RHS_ROOT)
        NULLIFY(id%root%RHS_ROOT)
      ENDIF
      CALL DMUMPS_RR_FREE_POINTERS(id)
      IF (associated(id%ELTPROC))     THEN
        DEALLOCATE(id%ELTPROC)
        NULLIFY(id%ELTPROC)
      ENDIF
C     id%CANDIDATES,id%I_AM_CAND and id%ISTEP_TO_INIV2
C     can be allocated on non-working master
C     in the case of arrowheads distribution
      IF (associated(id%CANDIDATES)) THEN
        DEALLOCATE(id%CANDIDATES)
        NULLIFY(id%CANDIDATES)
      ENDIF
      IF (associated(id%I_AM_CAND)) THEN
        DEALLOCATE(id%I_AM_CAND)
        NULLIFY(id%I_AM_CAND)
      ENDIF
      IF (associated(id%ISTEP_TO_INIV2)) THEN
        DEALLOCATE(id%ISTEP_TO_INIV2)
        NULLIFY(id%ISTEP_TO_INIV2)
      ENDIF
C     Node partitionning (only allocated on slaves)
      IF (I_AM_SLAVE) THEN
       IF (associated(id%TAB_POS_IN_PERE)) THEN
        DEALLOCATE(id%TAB_POS_IN_PERE)
        NULLIFY(id%TAB_POS_IN_PERE)
       ENDIF
       IF (associated(id%FUTURE_NIV2)) THEN
        DEALLOCATE(id%FUTURE_NIV2)
        NULLIFY(id%FUTURE_NIV2)
       ENDIF
      ENDIF
      IF(associated(id%DEPTH_FIRST))THEN
        DEALLOCATE(id%DEPTH_FIRST)
        NULLIFY(id%DEPTH_FIRST)
      ENDIF
      IF(associated(id%DEPTH_FIRST_SEQ))THEN
        DEALLOCATE(id%DEPTH_FIRST_SEQ)
        NULLIFY(id%DEPTH_FIRST_SEQ)
      ENDIF
      IF(associated(id%SBTR_ID))THEN
        DEALLOCATE(id%SBTR_ID)
        NULLIFY(id%SBTR_ID)
      ENDIF
      IF(associated(id%SCHED_DEP))THEN
        DEALLOCATE(id%SCHED_DEP)
        NULLIFY(id%SCHED_DEP)
      ENDIF
      IF(associated(id%SCHED_SBTR))THEN
        DEALLOCATE(id%SCHED_SBTR)
        NULLIFY(id%SCHED_SBTR)
      ENDIF
      IF(associated(id%SCHED_GRP))THEN
        DEALLOCATE(id%SCHED_GRP)
        NULLIFY(id%SCHED_GRP)
      ENDIF
      IF(associated(id%CROIX_MANU))THEN
        DEALLOCATE(id%CROIX_MANU)
        NULLIFY(id%CROIX_MANU)
      ENDIF
      IF (associated(id%MEM_SUBTREE)) THEN
        DEALLOCATE(id%MEM_SUBTREE)
        NULLIFY(id%MEM_SUBTREE)
      ENDIF
      IF (associated(id%MY_ROOT_SBTR)) THEN
        DEALLOCATE(id%MY_ROOT_SBTR)
        NULLIFY(id%MY_ROOT_SBTR)
      ENDIF
      IF (associated(id%MY_FIRST_LEAF)) THEN
        DEALLOCATE(id%MY_FIRST_LEAF)
        NULLIFY(id%MY_FIRST_LEAF)
      ENDIF
      IF (associated(id%MY_NB_LEAF)) THEN
        DEALLOCATE(id%MY_NB_LEAF)
        NULLIFY(id%MY_NB_LEAF)
      ENDIF
      IF (associated(id%COST_TRAV)) THEN
        DEALLOCATE(id%COST_TRAV)
        NULLIFY(id%COST_TRAV)
      ENDIF
      IF (associated(id%CB_SON_SIZE)) THEN
        DEALLOCATE(id%CB_SON_SIZE)
        NULLIFY(id%CB_SON_SIZE)
      ENDIF
      IF (associated(id%SUP_PROC)) THEN
         DEALLOCATE(id%SUP_PROC)
         NULLIFY(id%SUP_PROC)
      ENDIF
c     IF (id%KEEP(201).GT.0) THEN
      IF(associated (id%OOC_INODE_SEQUENCE))THEN
         DEALLOCATE(id%OOC_INODE_SEQUENCE)
         NULLIFY(id%OOC_INODE_SEQUENCE)
      ENDIF
      IF(associated (id%OOC_TOTAL_NB_NODES))THEN
         DEALLOCATE(id%OOC_TOTAL_NB_NODES)
         NULLIFY(id%OOC_TOTAL_NB_NODES)
      ENDIF
      IF(associated (id%OOC_SIZE_OF_BLOCK))THEN
         DEALLOCATE(id%OOC_SIZE_OF_BLOCK)
         NULLIFY(id%OOC_SIZE_OF_BLOCK)
      ENDIF
      IF(associated (id%OOC_VADDR))THEN
         DEALLOCATE(id%OOC_VADDR)
         NULLIFY(id%OOC_VADDR)
      ENDIF
      IF(associated (id%OOC_NB_FILES))THEN
         DEALLOCATE(id%OOC_NB_FILES)
         NULLIFY(id%OOC_NB_FILES)
      ENDIF
c     ENDIF
!     IF(id%KEEP(486).NE.0) THEN
        IF (associated(id%LRGROUPS)) THEN
           DEALLOCATE(id%LRGROUPS)
           NULLIFY(id%LRGROUPS)
        ENDIF
!     ENDIF
      CALL DMUMPS_FREE_ID_DATA_MODULES(id%FDM_F_ENCODING,
     &  id%BLRARRAY_ENCODING, id%KEEP8(1))
      IF (associated(id%MPITOOMP_PROCS_MAP)) THEN
        DEALLOCATE(id%MPITOOMP_PROCS_MAP)
        NULLIFY(id%MPITOOMP_PROCS_MAP)
      ENDIF
        IF (associated(id%SINGULAR_VALUES)) THEN
           DEALLOCATE(id%SINGULAR_VALUES)
           NULLIFY(id%SINGULAR_VALUES)
        ENDIF        
C     ----------------------------------------------
C     Deallocate S only after finishing the receives
C     (S is normally the largest memory available)
C     ----------------------------------------------
      IF (id%KEEP8(24).EQ.0_8) THEN
C       -- deallocate only when not provided/allocated by the user
        IF (associated(id%S))        DEALLOCATE(id%S)
      ENDIF
      NULLIFY(id%S)
      IF (I_AM_SLAVE) THEN
C       ------------------------
C       Deallocate buffer for
C       contrib-blocks (facto/
C       solve). Note that this
C       will cancel all possible
C       pending requests.
C       ------------------------
        CALL DMUMPS_BUF_DEALL_CB( IERR )
C       Deallocate buffer for integers (facto/solve)
        CALL DMUMPS_BUF_DEALL_SMALL_BUF( IERR )
      END IF
C Mapping information used during solve
      IF (associated(id%IPTR_WORKING)) THEN
        DEALLOCATE(id%IPTR_WORKING)
        NULLIFY(id%IPTR_WORKING)
      END IF
      IF (associated(id%WORKING)) THEN 
        DEALLOCATE(id%WORKING)
        NULLIFY(id%WORKING)
      END IF
      IF (associated(id%IPOOL_B_L0_OMP)) THEN
        DEALLOCATE(id%IPOOL_B_L0_OMP)
        NULLIFY(id%IPOOL_B_L0_OMP)
      END IF
      IF (associated(id%IPOOL_A_L0_OMP)) THEN
        DEALLOCATE(id%IPOOL_A_L0_OMP)
        NULLIFY(id%IPOOL_A_L0_OMP)
      END IF
      IF (associated(id%PHYS_L0_OMP)) THEN
        DEALLOCATE(id%PHYS_L0_OMP)
        NULLIFY(id%PHYS_L0_OMP)
      END IF
      IF (associated(id%VIRT_L0_OMP)) THEN
        DEALLOCATE(id%VIRT_L0_OMP)
        NULLIFY(id%VIRT_L0_OMP)
      END IF
      IF (associated(id%VIRT_L0_OMP_MAPPING)) THEN
        DEALLOCATE(id%VIRT_L0_OMP_MAPPING)
        NULLIFY(id%VIRT_L0_OMP_MAPPING)
      END IF
      IF (associated(id%PERM_L0_OMP)) THEN
        DEALLOCATE(id%PERM_L0_OMP)
        NULLIFY(id%PERM_L0_OMP)
      END IF
      IF (associated(id%PTR_LEAFS_L0_OMP)) THEN
        DEALLOCATE(id%PTR_LEAFS_L0_OMP)
        NULLIFY(id%PTR_LEAFS_L0_OMP)
      END IF
      IF (associated(id%L0_OMP_MAPPING)) THEN
        DEALLOCATE(id%L0_OMP_MAPPING)
        NULLIFY(id%L0_OMP_MAPPING)
      END IF
      IF (associated(id%I4_L0_OMP)) THEN
        DEALLOCATE(id%I4_L0_OMP)
        NULLIFY(id%I4_L0_OMP)
      END IF
      IF (associated(id%I8_L0_OMP)) THEN
        DEALLOCATE(id%I8_L0_OMP)
        NULLIFY(id%I8_L0_OMP)
      END IF
      RETURN
      END SUBROUTINE DMUMPS_END_DRIVER
      SUBROUTINE DMUMPS_FREE_ID_DATA_MODULES(id_FDM_F_ENCODING,
     &  id_BLRARRAY_ENCODING, KEEP8)
      USE MUMPS_FRONT_DATA_MGT_M, only : MUMPS_FDM_STRUC_TO_MOD,
     &                                   MUMPS_FDM_END
      USE DMUMPS_LR_DATA_M, only : DMUMPS_BLR_STRUC_TO_MOD,
     &                             DMUMPS_BLR_END_MODULE
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Free data from modules kept from one phase to the other
C     and referenced through the main MUMPS structure, id.
C
C     Both id%FDM_F_ENCODING and id%BLRARRAY_ENCODING
C     are concerned.
C
C
C
C     Arguments:
C     =========
C
#     if defined(MUMPS_F2003)
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_FDM_F_ENCODING
#     else
      CHARACTER, DIMENSION(:), POINTER :: id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER :: id_FDM_F_ENCODING
#     endif
      INTEGER(8), intent(inout) :: KEEP8(150)
C
      IF (associated(id_FDM_F_ENCODING)) THEN
C       Allow access to FDM_F data for BLR_END_MODULE
        CALL MUMPS_FDM_STRUC_TO_MOD('F', id_FDM_F_ENCODING)
        IF (associated(id_BLRARRAY_ENCODING)) THEN
C         Pass id_BLRARRAY_ENCODING control to module
C         and terminate BLR module of current instance
          CALL DMUMPS_BLR_STRUC_TO_MOD(id_BLRARRAY_ENCODING)
          CALL DMUMPS_BLR_END_MODULE(0, KEEP8, 
     &                               LRSOLVE_ACT_OPT=.TRUE.)
        ENDIF
C       ---------------------------------------
C       FDM data structures are still allocated
C       in the module and should be freed
C       ---------------------------------------
        CALL MUMPS_FDM_END('F')
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FREE_ID_DATA_MODULES
