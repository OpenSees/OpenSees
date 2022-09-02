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
      RECURSIVE SUBROUTINE DMUMPS_PROCESS_BLFAC_SLAVE(
     &   COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, PERM,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, OPASSW, OPELIW,
     &    ITLOC, RHS_MUMPS, FILS, DAD,
     &    PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL,KEEP,KEEP8,DKEEP,IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &               , LRGROUPS
     &     )
      USE DMUMPS_BUF
      USE DMUMPS_LOAD
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_TYPE
      USE DMUMPS_FAC_LR
      USE DMUMPS_LR_DATA_M
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR
      USE DMUMPS_FAC_FRONT_AUX_M, 
     &                ONLY : DMUMPS_COMPUTE_SIZE_SCHUR_IN_FRONT
#if defined(BLR_MT)
!$    USE OMP_LIB
#endif
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 60 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION       DKEEP(230)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS, LA
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER PERM(N), STEP(N), PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(N)
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER(8), INTENT(IN) :: PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28))
      INTEGER ITLOC( N + KEEP(253)), FILS( N ), DAD( KEEP(28) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER ND( KEEP(28) ), FRERE_STEPS( KEEP(28) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      DOUBLE PRECISION DBLARR( KEEP8(26) )
      INTEGER INTARR( KEEP8(27) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INTEGER INODE, IPOSK, JPOSK, NCOLU, NPIV, POSITION, IERR
      INTEGER(8) POSELT, POSBLOCFACTO
      INTEGER(8) LAELL
      INTEGER(8) :: LA_PTR 
      DOUBLE PRECISION, DIMENSION(:), POINTER :: A_PTR
      INTEGER IOLDPS, LCONT1, NROW1, NCOL1, NPIV1
      INTEGER NSLAVES_TOT, HS, DEST, NSLAVES_FOLLOW
      INTEGER FPERE
      INTEGER(8) CPOS, LPOS
      LOGICAL DYNAMIC
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER allocok
      INTEGER LR_ACTIVATED_INT
      LOGICAL LR_ACTIVATED, COMPRESS_CB
      INTEGER NB_BLR_U, CURRENT_BLR_U
      TYPE (LRB_TYPE), DIMENSION(:), ALLOCATABLE :: BLR_U
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_U
      TYPE (LRB_TYPE), DIMENSION(:), POINTER :: BLR_LS
      TYPE (LRB_TYPE), POINTER, DIMENSION(:,:) :: CB_LRB
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_LS, BEGS_BLR_COL
      INTEGER    :: NB_BLR_LS, IPANEL,
     &           MAXI_CLUSTER_LS, MAXI_CLUSTER, 
     &           NB_BLR_COL, MAXI_CLUSTER_COL, NPARTSASS_MASTER
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK, TAU
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JPVT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BLOCKLR
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RWORK
      INTEGER :: OMP_NUM, LWORK
      INTEGER :: II,JJ
       INTEGER :: NFS4FATHER, NASS1, NELIM, INFO_TMP(2)
       INTEGER :: NVSCHUR_K253, NSLAVES_L, IROW_L
       INTEGER :: NBROWSinF
       DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: M_ARRAY
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: UDYNAMIC
      DOUBLE PRECISION ONE,ALPHA
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      DYNAMIC = .FALSE.
      POSITION  = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPOSK, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, JPOSK, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF ( NPIV .LE. 0 ) THEN
        NPIV = - NPIV
        WRITE(*,*) MYID,':error, received negative NPIV in BLFAC'
        CALL MUMPS_ABORT()
      END IF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOLU, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 LR_ACTIVATED_INT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LR_ACTIVATED   = (LR_ACTIVATED_INT.EQ.1)
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 IPANEL, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF (LR_ACTIVATED) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NB_BLR_U, 1, MPI_INTEGER,
     &                 COMM, IERR )
          CURRENT_BLR_U = 1 
          ALLOCATE(BLR_U(max(NB_BLR_U,1)),
     &             BEGS_BLR_U(NB_BLR_U+2), stat=allocok)
          if (allocok .GT. 0) THEN
             IFLAG = -13 
             IERROR = max(NB_BLR_U,1) + NB_BLR_U+2
             GOTO 700
          endif
          CALL DMUMPS_MPI_UNPACK_LR(BUFR, LBUFR, LBUFR_BYTES,
     &                        POSITION, JPOSK-1, 0, 'V',
     &                        BLR_U, NB_BLR_U, 
     &                        BEGS_BLR_U(1),
     &                        KEEP8, COMM, IERR, IFLAG, IERROR)
          IF (IFLAG.LT.0) GOTO 700
      ELSE
      LAELL = int(NPIV,8) * int(NCOLU,8)
      CALL DMUMPS_GET_SIZE_NEEDED(
     &        0, LAELL, .FALSE.,
     &        KEEP(1), KEEP8(1),
     &        N, KEEP(28), IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &        KEEP(IXSZ),COMP,DKEEP(97),MYID, SLAVEF,
     &        PROCNODE_STEPS, DAD, 
     &        IFLAG, IERROR)
      IF (IFLAG.LT.0) GOTO 700
      LRLU  = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      KEEP8(67) = min(LRLUS, KEEP8(67))
      KEEP8(69) = KEEP8(69) + LAELL
      KEEP8(68) = max(KEEP8(69), KEEP8(68))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                           LA-LRLUS,0_8, LAELL,KEEP,KEEP8,LRLUS)
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*NCOLU,
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
      ENDIF
      IF (PTRIST(STEP( INODE )) .EQ. 0) DYNAMIC = .TRUE.
      IF ( (PTRIST(STEP( INODE )).NE.0) .AND.
     &  (IPOSK + NPIV -1 .GT.
     &   IW(PTRIST(STEP(INODE))+3+KEEP(IXSZ))) )THEN
        DYNAMIC = .TRUE.
      ENDIF
      IF (LR_ACTIVATED) THEN
       DYNAMIC = .FALSE.
      ENDIF
      IF (DYNAMIC)  THEN
        ALLOCATE(UDYNAMIC(LAELL), stat=allocok)
        if (allocok .GT. 0) THEN
          IFLAG = -13 
          CALL MUMPS_SET_IERROR(LAELL,IERROR)
          GOTO 700
        endif
        UDYNAMIC(1_8:LAELL) = A(POSBLOCFACTO:POSBLOCFACTO+LAELL-1_8)
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        KEEP8(69) = KEEP8(69) - LAELL
        POSFAC = POSFAC - LAELL
      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &          LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLUS)
      ENDIF
      IF ( PTRIST(STEP(INODE)) .EQ. 0 ) THEN
        CALL DMUMPS_TREAT_DESCBAND(INODE, COMM_LOAD,
     &    ASS_IRECV, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    PERM,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &     )
        IF ( IFLAG .LT. 0 ) GOTO 600
      ENDIF
      DO WHILE ( IPOSK + NPIV -1 .GT.
     &            IW( PTRIST(STEP( INODE )) + 3 +KEEP(IXSZ)) )
        MSGSOU = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE)),
     &                           KEEP(199) )
        SET_IRECV = .FALSE.
        BLOCKING  = .TRUE.
        MESSAGE_RECEIVED = .FALSE.
        CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MSGSOU, BLOC_FACTO_SYM, STATUS, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &      )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &    STATUS, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &     )
      IOLDPS  = PTRIST(STEP( INODE ))
      CALL DMUMPS_DM_SET_DYNPTR( IW(IOLDPS+XXS), A, LA,
     &     PTRAST(STEP(INODE)), IW(IOLDPS+XXD), IW(IOLDPS+XXR),
     &     A_PTR, POSELT, LA_PTR )
      LCONT1  = IW( IOLDPS + KEEP(IXSZ) )
      NROW1   = IW( IOLDPS + 2  + KEEP(IXSZ))
      NPIV1   = IW( IOLDPS + 3  + KEEP(IXSZ))
      NSLAVES_TOT = IW( IOLDPS + 5  + KEEP(IXSZ))
      HS      = 6 + NSLAVES_TOT + KEEP(IXSZ)
      NCOL1   = LCONT1 + NPIV1
      IF (LR_ACTIVATED) THEN
         CALL DMUMPS_BLR_DEC_AND_RETRIEVE_L (IW(IOLDPS+XXF), IPANEL, 
     &        BEGS_BLR_LS, BLR_LS)
         NB_BLR_LS = size(BEGS_BLR_LS)-2   
#if defined(BLR_MT)          
!$OMP PARALLEL
#endif
         CALL DMUMPS_BLR_UPDATE_TRAILING_I (
     &        A_PTR(POSELT), LA_PTR, 1_8,
     &        IFLAG, IERROR, NCOL1,
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS),
     &        BEGS_BLR_U(1), size(BEGS_BLR_U),
     &        CURRENT_BLR_U, 
     &        BLR_LS(1),  NB_BLR_LS+1,
     &        BLR_U(1), NB_BLR_U+1,
     &        0,       
     &        .TRUE.,  
     &        0,       
     &        2,       
     &        1,
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &         )
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          CALL DEALLOC_BLR_PANEL (BLR_U, NB_BLR_U, KEEP8)
          IF (allocated(BLR_U)) DEALLOCATE(BLR_U)
          IF (associated(BEGS_BLR_U)) DEALLOCATE(BEGS_BLR_U)
          IF (IFLAG.LT.0) GOTO 700
         IF (KEEP(486).EQ.3) THEN
         CALL DMUMPS_BLR_TRY_FREE_PANEL(IW(IOLDPS+XXF), IPANEL, 
     &                            KEEP8) 
         ENDIF
      ELSE
      CPOS = POSELT + int(JPOSK - 1,8)
      LPOS = POSELT + int(IPOSK - 1,8)
      IF ( NPIV .GT. 0 ) THEN
          IF (DYNAMIC) THEN
            CALL dgemm( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     &            UDYNAMIC(1), NPIV,
     &            A_PTR( LPOS ), NCOL1, ONE,
     &            A_PTR( CPOS ), NCOL1 )
          ELSE
            CALL dgemm( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     &            A( POSBLOCFACTO ), NPIV,
     &            A_PTR( LPOS ), NCOL1, ONE,
     &            A_PTR( CPOS ), NCOL1 )
          ENDIF
      ENDIF 
      ENDIF
      IF (NPIV .GT. 0) THEN
        FLOP1 = dble(NCOLU*NPIV)*dble(2*NROW1)
        FLOP1 = -FLOP1
        CALL DMUMPS_LOAD_UPDATE(1, .FALSE., FLOP1, KEEP,KEEP8 )
      ENDIF
      IF ( IW(IOLDPS+6+KEEP(IXSZ)).EQ.
     &    huge(IW(IOLDPS+6+KEEP(IXSZ))) ) THEN
          IW(IOLDPS+6+KEEP(IXSZ)) = 1
      ENDIF
      IW(IOLDPS+6+KEEP(IXSZ)) =
     &         IW(IOLDPS+6+KEEP(IXSZ)) + 1
      IF (.NOT.LR_ACTIVATED) THEN
      IF (DYNAMIC) THEN
       DEALLOCATE(UDYNAMIC)
      ELSE
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        KEEP8(69) = KEEP8(69) - LAELL
        POSFAC = POSFAC - LAELL
      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                      LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLUS)
      ENDIF
      ENDIF
      NSLAVES_FOLLOW = IW( IOLDPS + 5 +KEEP(IXSZ) ) - XTRA_SLAVES_SYM
      IF ( IW( IOLDPS + 6  +KEEP(IXSZ)) .eq. 0 .and.
     &     KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 0 )
     &     THEN
         DEST = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE)), KEEP(199) )
         CALL DMUMPS_BUF_SEND_1INT( INODE, DEST, END_NIV2_LDLT,
     &                             COMM, KEEP, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_BLFAC_SLAVE.'
           IFLAG = -99
           GOTO 700
         END IF
      END IF
      IF (IW(PTRIST(STEP(INODE)) + 6+KEEP(IXSZ) ) .eq. 0) THEN
           NPIV1 = IW( IOLDPS + 3  + KEEP(IXSZ))
           NASS1 = IW( IOLDPS + 4 + KEEP(IXSZ))  
           NELIM = NASS1 - NPIV1
          COMPRESS_CB= .FALSE.
          IF (LR_ACTIVATED) THEN
            COMPRESS_CB = ((IW(PTRIST(STEP(INODE))+XXLR).EQ.1).OR.
     &                     (IW(PTRIST(STEP(INODE))+XXLR).EQ.3))
             IF (COMPRESS_CB.AND.NPIV.EQ.0) THEN
              COMPRESS_CB = .FALSE.
              IW(IOLDPS+XXLR) = IW(IOLDPS+XXLR) -1
             ENDIF
            IF (COMPRESS_CB) THEN
              CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_C (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_COL,  NPARTSASS_MASTER)
              NB_BLR_COL   = size(BEGS_BLR_COL) - 1
              allocate(CB_LRB(NB_BLR_LS,NB_BLR_COL-NPARTSASS_MASTER),
     &                 stat=allocok)
              IF (allocok > 0) THEN
                IFLAG  = -13
                IERROR = NB_BLR_LS*(NB_BLR_COL-NPARTSASS_MASTER)
                GOTO 700
              ENDIF
              DO II=1,NB_BLR_LS
              DO JJ=1,NB_BLR_COL-NPARTSASS_MASTER
                NULLIFY(CB_LRB(II,JJ)%Q)
                NULLIFY(CB_LRB(II,JJ)%R)
                CB_LRB(II,JJ)%ISLR = .FALSE.
              ENDDO
              ENDDO
              CALL DMUMPS_BLR_SAVE_CB_LRB(IW(IOLDPS+XXF),CB_LRB)
              call MAX_CLUSTER(BEGS_BLR_LS,NB_BLR_LS+1,MAXI_CLUSTER_LS)
              CALL MAX_CLUSTER(BEGS_BLR_COL,NB_BLR_COL,MAXI_CLUSTER_COL)
              MAXI_CLUSTER = max(MAXI_CLUSTER_LS,
     &         MAXI_CLUSTER_COL+NELIM,NPIV)
              LWORK = MAXI_CLUSTER*MAXI_CLUSTER
              OMP_NUM = 1
#if defined(BLR_MT)
!$            OMP_NUM = OMP_GET_MAX_THREADS()
#endif
              ALLOCATE(BLOCKLR(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &            RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &            TAU(MAXI_CLUSTER*OMP_NUM),
     &            JPVT(MAXI_CLUSTER*OMP_NUM), 
     &            WORK(LWORK*OMP_NUM),
     &            stat=allocok)
              IF (allocok > 0 ) THEN
                IFLAG  = -13
                IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
                GOTO 700
              ENDIF
              NFS4FATHER = -9999
              IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2) ) THEN
               CALL DMUMPS_BLR_RETRIEVE_NFS4FATHER ( IW(IOLDPS+XXF),
     &             NFS4FATHER )
               NFS4FATHER = max(NFS4FATHER,0) + NELIM
              ENDIF
              ALLOCATE(M_ARRAY(max(NFS4FATHER,1)), stat=allocok)
              IF (allocok.gt.0) THEN
                IFLAG = -13
                IERROR = max(NFS4FATHER,1)
                GOTO 700
              ENDIF
              BEGS_BLR_COL(1+NPARTSASS_MASTER) = 
     &               BEGS_BLR_COL(1+NPARTSASS_MASTER) - NELIM
              NBROWSinF = 0
              IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2).AND.
     &             NFS4FATHER.GT.0  ) THEN
                CALL DMUMPS_COMPUTE_NBROWSinF (
     &                N, INODE, FPERE, KEEP, 
     &                IOLDPS, HS, 
     &                IW, LIW, 
     &                NROW1, NCOL1, NPIV1,
     &                NELIM, NFS4FATHER,
     &                NBROWSinF
     &                )
              ENDIF
              IF ((KEEP(114).EQ.1) .AND. (KEEP(116).GT.0)
     &            .AND. (KEEP(50).EQ.2)
     &           ) THEN
                 NSLAVES_L = IW(PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ))
                 IROW_L    = PTRIST(STEP(INODE)) + 6 + NSLAVES_L + 
     &                       KEEP(IXSZ)
                 CALL DMUMPS_COMPUTE_SIZE_SCHUR_IN_FRONT ( 
     &                N, 
     &                NROW1,
     &                KEEP(116), 
     &                IW(IROW_L),
     &                PERM, NVSCHUR_K253 )
              ELSE
                 NVSCHUR_K253 = 0
              ENDIF
#if defined(BLR_MT)          
!$OMP PARALLEL
#endif
              CALL DMUMPS_COMPRESS_CB_I(
     &          A_PTR(POSELT), LA_PTR, 1_8, NCOL1,
     &          BEGS_BLR_LS(1), size(BEGS_BLR_LS),
     &          BEGS_BLR_COL(1), size(BEGS_BLR_COL),
     &          NB_BLR_LS, NB_BLR_COL-NPARTSASS_MASTER,
     &          NPARTSASS_MASTER, 
     &          NROW1, NCOL1-NPIV1, INODE,
     &          IW(IOLDPS+XXF), 1, 2, IFLAG, IERROR,
     &          DKEEP(12), KEEP(466), KEEP(484), KEEP(489),
     &          CB_LRB(1,1),
     &          WORK, TAU, JPVT, LWORK, RWORK, BLOCKLR,
     &          MAXI_CLUSTER, KEEP8, OMP_NUM,
     &          NFS4FATHER, NPIV1, NVSCHUR_K253, KEEP(1), 
     &          M_ARRAY,
     &          NELIM, NBROWSinF )
#if defined(BLR_MT)
!$OMP END PARALLEL
#endif
              IF (IFLAG.LT.0) GOTO 650
              IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2).AND.
     &             NFS4FATHER.GT.0  ) THEN
                 INFO_TMP(1) = IFLAG
                 INFO_TMP(2) = IERROR
                 CALL DMUMPS_BLR_SAVE_M_ARRAY( IW(IOLDPS+XXF),
     &            M_ARRAY, INFO_TMP)
                 IFLAG  = INFO_TMP(1) 
                 IERROR = INFO_TMP(2) 
              ENDIF
 650          CONTINUE         
              IF (allocated(M_ARRAY)) DEALLOCATE(M_ARRAY)
              IF (allocated(BLOCKLR)) DEALLOCATE(BLOCKLR)
              IF (allocated(RWORK)) DEALLOCATE(RWORK)
              IF (allocated(TAU)) DEALLOCATE(TAU)
              IF (allocated(JPVT)) DEALLOCATE(JPVT)
              IF (allocated(WORK)) DEALLOCATE(WORK)           
              IF (IFLAG.LT.0) GOTO 700
            ENDIF 
          ENDIF
         CALL DMUMPS_END_FACTO_SLAVE( COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    root,
     &    MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, PERM,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD, PTRARW, PTRAIW,
     &    INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &               , LRGROUPS
     &     )
       ENDIF 
      RETURN
 700  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
 600  CONTINUE
      IF (allocated(BLR_U)) DEALLOCATE(BLR_U)
      IF (COMPRESS_CB) THEN
        IF (allocated(BLOCKLR)) DEALLOCATE(BLOCKLR)
        IF (allocated(RWORK)) DEALLOCATE(RWORK)
        IF (allocated(TAU)) DEALLOCATE(TAU)
        IF (allocated(JPVT)) DEALLOCATE(JPVT)
        IF (allocated(WORK)) DEALLOCATE(WORK)
      ENDIF
      IF (allocated(M_ARRAY)) DEALLOCATE(M_ARRAY)
      IF (DYNAMIC) THEN
       IF (allocated(UDYNAMIC)) DEALLOCATE(UDYNAMIC)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_BLFAC_SLAVE
