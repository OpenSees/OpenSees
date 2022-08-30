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
      RECURSIVE SUBROUTINE DMUMPS_PROCESS_BLOCFACTO(
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
     &    ICNTL, KEEP,KEEP8, DKEEP, 
     &    IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &    )
      USE DMUMPS_OOC, ONLY : IO_BLOCK
      USE MUMPS_OOC_COMMON, ONLY : TYPEF_L,
     &                       STRAT_WRITE_MAX,
     &                       STRAT_TRY_WRITE
      USE DMUMPS_LOAD
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_FAC_LR
      USE DMUMPS_ANA_LR, ONLY : GET_CUT
      USE DMUMPS_LR_DATA_M
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 60 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION DKEEP(230)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA
      INTEGER(8) :: POSFAC
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PERM(N), STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(N)
      INTEGER COMM, MYID
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PTLUST_S(KEEP(28)),
     &        ITLOC(N+KEEP(253)), FILS(N), DAD(KEEP(28)), ND(KEEP(28))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8), INTENT(IN) :: PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      INTEGER INTARR( KEEP8(27) )
      DOUBLE PRECISION DBLARR( KEEP8(26) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      LOGICAL :: I_HAVE_SET_K117
      INTEGER INODE, POSITION, NPIV, IERR, LP
      INTEGER NCOL
      INTEGER(8) :: POSBLOCFACTO
      INTEGER :: LD_BLOCFACTO 
      INTEGER(8) :: LA_BLOCFACTO 
      INTEGER(8) :: LA_PTR 
      INTEGER(8) :: POSELT
      DOUBLE PRECISION, DIMENSION(:), POINTER :: A_PTR
      INTEGER IOLDPS, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW
      INTEGER (8) :: LPOS, UPOS, LPOS2, IPOS, KPOS
      INTEGER ICT11
      INTEGER I, IPIV, FPERE
      LOGICAL LASTBL, KEEP_BEGS_BLR_L
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      DOUBLE PRECISION ONE,ALPHA
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INTEGER LIWFAC, STRAT, NextPivDummy
      TYPE(IO_BLOCK) :: MonBloc
      LOGICAL LAST_CALL
      INTEGER LRELAY_INFO
      INTEGER :: INFO_TMP(2)
      INTEGER :: NELIM, NPARTSASS_MASTER, NPARTSASS_MASTER_AUX,
     &           IPANEL, 
     &           CURRENT_BLR, 
     &           NB_BLR_L, NB_BLR_U, NB_BLR_COL
      TYPE (LRB_TYPE), POINTER, DIMENSION(:,:) :: CB_LRB
      TYPE (LRB_TYPE), DIMENSION(:), POINTER :: BLR_U, BLR_L
      LOGICAL :: LR_ACTIVATED, COMPRESS_CB, COMPRESS_PANEL
      LOGICAL OOCWRITE_COMPATIBLE_WITH_BLR
      INTEGER :: LR_ACTIVATED_INT
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_L, BEGS_BLR_U,
     & BEGS_BLR_COL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK, TAU
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JPVT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RWORK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BLOCK
      INTEGER :: OMP_NUM
      INTEGER NPARTSASS, NPARTSCB, MAXI_CLUSTER, LWORK,
     &        MAXI_CLUSTER_L, MAXI_CLUSTER_U, MAXI_CLUSTER_COL
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DYN_BLOCFACTO
      INTEGER, DIMENSION(:), ALLOCATABLE :: DYN_PIVINFO
      LOGICAL :: DYNAMIC_ALLOC
      INTEGER  :: allocok
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      KEEP_BEGS_BLR_L = .FALSE.
      nullify(BEGS_BLR_L)
      NB_BLR_U = -7654321
      NULLIFY(BEGS_BLR_U)
      I_HAVE_SET_K117 = .FALSE.
      DYNAMIC_ALLOC = .FALSE.
      FPERE    = -1
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NELIM, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 NPARTSASS_MASTER , 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPANEL,
     &                 1, MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, LR_ACTIVATED_INT,
     &                 1, MPI_INTEGER, COMM, IERR )
      LR_ACTIVATED    = (LR_ACTIVATED_INT.EQ.1)
      IF ( LR_ACTIVATED ) THEN
        LA_BLOCFACTO = int(NPIV,8) * int(NPIV+NELIM,8)
      ELSE
        LA_BLOCFACTO = int(NPIV,8) * int(NCOL,8)
      ENDIF
      CALL DMUMPS_GET_SIZE_NEEDED(
     &      NPIV, LA_BLOCFACTO, .FALSE.,
     &      KEEP(1), KEEP8(1),
     &      N, KEEP(28), IW, LIW, A, LA,
     &      LRLU, IPTRLU,
     &      IWPOS, IWPOSCB, PTRIST, PTRAST,
     &      STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &      KEEP(IXSZ),COMP,DKEEP(97),MYID,SLAVEF, PROCNODE_STEPS, 
     &      DAD, IFLAG, IERROR)
      IF (IFLAG.LT.0) GOTO 700
      LRLU  = LRLU - LA_BLOCFACTO
      LRLUS = LRLUS - LA_BLOCFACTO
      KEEP8(67) = min(LRLUS, KEEP8(67))
      KEEP8(69) = KEEP8(69) + LA_BLOCFACTO
      KEEP8(68) = max(KEEP8(69), KEEP8(68))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LA_BLOCFACTO
      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE., .FALSE.,
     &               LA-LRLUS,0_8,LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
      IF ((NPIV .EQ. 0) 
     &     ) THEN
        IPIV=1
      ELSE
        IPIV = IWPOS
        IWPOS = IWPOS + NPIV
        IF (NPIV .GT. 0) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 IW( IPIV ), NPIV,
     &                 MPI_INTEGER, COMM, IERR )
        ENDIF
        IF ( LR_ACTIVATED ) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*(NPIV+NELIM),
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
            LD_BLOCFACTO = NPIV+NELIM
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NB_BLR_U, 1, MPI_INTEGER,
     &                 COMM, IERR )
            ALLOCATE(BLR_U(max(NB_BLR_U,1)), stat=allocok)
            IF (allocok > 0 ) THEN
               IFLAG = -13
               IERROR = max(NB_BLR_U,1)
               LP = ICNTL(1)
               IF (ICNTL(4) .LE. 0) LP=-1
               IF (LP > 0) WRITE(LP,*) MYID,
     &              ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
               GOTO 700
            ENDIF
            ALLOCATE(BEGS_BLR_U(NB_BLR_U+2), stat=allocok)
            IF (allocok > 0 ) THEN
               IFLAG = -13
               IERROR = NB_BLR_U+2
               LP = ICNTL(1)
               IF (ICNTL(4) .LE. 0) LP=-1
               IF (LP > 0) WRITE(LP,*) MYID,
     &              ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
               GOTO 700
            ENDIF
            CALL DMUMPS_MPI_UNPACK_LR(BUFR, LBUFR, LBUFR_BYTES, 
     &                             POSITION, NPIV, NELIM, 'H',
     &                             BLR_U(1), NB_BLR_U, 
     &                             BEGS_BLR_U(1),
     &                             KEEP8, COMM, IERR, IFLAG, IERROR)
            IF (IFLAG.LT.0) GOTO 700
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), NPIV*NCOL,
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
          LD_BLOCFACTO = NCOL
        ENDIF
      ENDIF 
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 LRELAY_INFO, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
          CALL DMUMPS_TREAT_DESCBAND( INODE, COMM_LOAD,
     &    ASS_IRECV,
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
     &    INTARR, DBLARR, ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS, 
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &    )
          IF ( IFLAG .LT. 0 ) GOTO 600
      ENDIF
      IF ( IW( PTRIST(STEP(INODE)) + 3 +KEEP(IXSZ)) .EQ. 0 ) THEN
       DO WHILE ( IW(PTRIST(STEP(INODE)) + XXNBPR) .NE. 0)
        BLOCKING = .TRUE.
        SET_IRECV = .FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, CONTRIB_TYPE2,
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
     &    INTARR, DBLARR, ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &      )
        IF ( IFLAG .LT. 0 ) GOTO 600
       END  DO
      ENDIF
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
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
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &    )
      IOLDPS = PTRIST(STEP(INODE))
      CALL DMUMPS_DM_SET_DYNPTR( IW(IOLDPS+XXS), A, LA,
     &     PTRAST(STEP(INODE)), IW(IOLDPS+XXD), IW(IOLDPS+XXR),
     &     A_PTR, POSELT, LA_PTR )
      LCONT1 = IW( IOLDPS + KEEP(IXSZ))
      NASS1  = IW( IOLDPS + 1 + KEEP(IXSZ))
      COMPRESS_PANEL = (IW(IOLDPS+XXLR).GE.2)
      OOCWRITE_COMPATIBLE_WITH_BLR = 
     &          ( .NOT.LR_ACTIVATED.OR.  (.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(486).NE.2) 
     &          )
      IF ( NASS1 < 0 ) THEN
        NASS1 = -NASS1
        IW( IOLDPS + 1 + KEEP(IXSZ)) = NASS1
        IF (KEEP(55) .EQ. 0) THEN 
          CALL DMUMPS_ASM_SLAVE_ARROWHEADS(INODE, N, IW, LIW,
     &       IOLDPS, A_PTR(POSELT), LA_PTR, 1_8, KEEP, KEEP8, ITLOC,
     &       FILS, PTRAIW,
     &       PTRARW, INTARR, DBLARR, KEEP8(27), KEEP8(26), RHS_MUMPS, 
     &       LRGROUPS)
        ELSE
          CALL DMUMPS_ASM_SLAVE_ELEMENTS(INODE, N, NELT, IW, LIW,
     &       IOLDPS, A_PTR(POSELT), LA_PTR, 1_8, KEEP, KEEP8, ITLOC,
     &       FILS, PTRAIW,
     &       PTRARW, INTARR, DBLARR, KEEP8(27), KEEP8(26),
     &       FRTPTR, FRTELT, RHS_MUMPS, LRGROUPS)
        ENDIF
      ENDIF
      NROW1  = IW( IOLDPS + 2 +KEEP(IXSZ))
      NPIV1  = IW( IOLDPS + 3 +KEEP(IXSZ))
      NSLAV1 = IW( IOLDPS + 5 + KEEP(IXSZ))
      HS     = 6 + NSLAV1 + KEEP(IXSZ)
      NCOL1  = LCONT1 + NPIV1
      IF (NPIV.GT.0) THEN
        ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
        IF (DYNAMIC_ALLOC) THEN
          DO I = 1, NPIV
            IF (DYN_PIVINFO(I).EQ.I) CYCLE
            ISW = IW(ICT11+I)
            IW(ICT11+I) = IW(ICT11+DYN_PIVINFO(I))
            IW(ICT11+DYN_PIVINFO(I)) = ISW
            IPOS = POSELT + int(NPIV1 + I - 1,8)
            KPOS = POSELT + int(NPIV1 + DYN_PIVINFO(I) - 1,8)
            CALL dswap(NROW1, A_PTR(IPOS), NCOL1, A_PTR(KPOS), NCOL1)
          ENDDO
        ELSE
          DO I = 1, NPIV
            IF (IW(IPIV+I-1).EQ.I) CYCLE
            ISW = IW(ICT11+I)
            IW(ICT11+I) = IW(ICT11+IW(IPIV+I-1))
            IW(ICT11+IW(IPIV+I-1)) = ISW
            IPOS = POSELT + int(NPIV1 + I - 1,8)
            KPOS = POSELT + int(NPIV1 + IW(IPIV+I-1) - 1,8)
            CALL dswap(NROW1, A_PTR(IPOS), NCOL1, A_PTR(KPOS), NCOL1)
          ENDDO
        ENDIF
        LPOS2 = POSELT + int(NPIV1,8)
        LPOS  = LPOS2 + int(NPIV,8)
        IF ((.NOT. LR_ACTIVATED).OR.KEEP(475).EQ.0) THEN
          IF (DYNAMIC_ALLOC) THEN
            CALL dtrsm('L','L','N','N',NPIV, NROW1, ONE, 
     &           DYN_BLOCFACTO, LD_BLOCFACTO, A_PTR(LPOS2), NCOL1)
          ELSE
            CALL dtrsm('L','L','N','N',NPIV, NROW1, ONE, 
     &           A(POSBLOCFACTO), LD_BLOCFACTO,
     &           A_PTR(LPOS2), NCOL1)
          ENDIF
      ENDIF
      ENDIF 
      COMPRESS_CB = .FALSE.
      IF ( LR_ACTIVATED) THEN 
        COMPRESS_CB    = ((IW(IOLDPS+XXLR).EQ.1).OR.
     &                    (IW(IOLDPS+XXLR).EQ.3))
        IF (COMPRESS_CB.AND.NPIV.EQ.0) THEN
           COMPRESS_CB = .FALSE.
           IW(IOLDPS+XXLR) = IW(IOLDPS+XXLR) -1
        ENDIF
        IF (NPIV.NE.0) THEN
        IF ( (NPIV1.EQ.0) 
     &   ) THEN
          IOLDPS = PTRIST(STEP(INODE))
          CALL GET_CUT(IW(IOLDPS+HS:IOLDPS+HS+NROW1-1), 0,
     &                    NROW1, LRGROUPS, NPARTSCB, 
     &                    NPARTSASS, BEGS_BLR_L)
          CALL REGROUPING2(BEGS_BLR_L, NPARTSASS, 0, NPARTSCB,
     &                        NROW1-0, KEEP(488), .TRUE., KEEP(472))
          NB_BLR_L =  NPARTSCB
          IF (IPANEL.EQ.1) THEN
           BEGS_BLR_COL=>BEGS_BLR_U
          ELSE
           ALLOCATE(BEGS_BLR_COL(size(BEGS_BLR_U)+IPANEL-1),
     &               stat=allocok)
           IF (allocok > 0 ) THEN
               IFLAG = -13
               IERROR = size(BEGS_BLR_U)+IPANEL-1
               LP = ICNTL(1)
               IF (ICNTL(4) .LE. 0) LP=-1
               IF (LP > 0) WRITE(LP,*) MYID,
     &           ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
               GOTO 700
            ENDIF
            BEGS_BLR_COL(1:IPANEL-1) = 1
            DO I=1,size(BEGS_BLR_U)
               BEGS_BLR_COL(IPANEL+I-1) =  BEGS_BLR_U(I)
            ENDDO
          ENDIF
          INFO_TMP(1) = IFLAG
          INFO_TMP(2) = IERROR
          IF (IFLAG.LT.0) GOTO 700
          CALL DMUMPS_BLR_SAVE_INIT(IW(IOLDPS+XXF), 
     &           .FALSE.,       
     &           .TRUE.,        
     &           .TRUE.,        
     &           NPARTSASS_MASTER, 
     &           BEGS_BLR_L, 
     &           BEGS_BLR_COL, 
     &           huge(NPARTSASS_MASTER),
     &           INFO_TMP)
          IFLAG  = INFO_TMP(1) 
          IERROR = INFO_TMP(2) 
          IF (IPANEL.NE.1) THEN
            DEALLOCATE(BEGS_BLR_COL)
          ENDIF
          IF (IFLAG.LT.0) GOTO 700
        ELSE           
          CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_L (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_L)
          KEEP_BEGS_BLR_L = .TRUE.  
          NB_BLR_L  = size(BEGS_BLR_L) - 2
          NPARTSASS = 1
          NPARTSCB  = NB_BLR_L
        ENDIF 
      ENDIF
      ENDIF
      IF ( (NPIV .GT. 0)
     &   ) THEN
        IF (LR_ACTIVATED) THEN
        call MAX_CLUSTER(BEGS_BLR_L,NB_BLR_L+1,MAXI_CLUSTER_L)
        call MAX_CLUSTER(BEGS_BLR_U,NB_BLR_U+1,MAXI_CLUSTER_U)
        IF (LASTBL.AND.COMPRESS_CB) THEN
          MAXI_CLUSTER=max(MAXI_CLUSTER_U+NELIM,MAXI_CLUSTER_L)
        ELSE
          MAXI_CLUSTER=max(MAXI_CLUSTER_U,MAXI_CLUSTER_L)
        ENDIF
        LWORK = MAXI_CLUSTER*MAXI_CLUSTER
        OMP_NUM = 1
#if defined(BLR_MT)
!$      OMP_NUM = OMP_GET_MAX_THREADS()
#endif
        ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &       RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &       TAU(MAXI_CLUSTER*OMP_NUM),
     &       JPVT(MAXI_CLUSTER*OMP_NUM), 
     &       WORK(LWORK*OMP_NUM), stat=allocok)
        IF (allocok > 0 ) THEN
           IFLAG = -13
           IERROR = MAXI_CLUSTER*OMP_NUM*MAXI_CLUSTER
     &          + 2*MAXI_CLUSTER*OMP_NUM
     &          + MAXI_CLUSTER*OMP_NUM
     &          + MAXI_CLUSTER*OMP_NUM
     &          + LWORK*OMP_NUM
           LP = ICNTL(1)
           IF (ICNTL(4) .LE. 0) LP=-1
           IF (LP > 0) WRITE(LP,*) MYID,
     &          ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
           GOTO 700
        ENDIF
        CURRENT_BLR=1 
        ALLOCATE(BLR_L(NB_BLR_L), stat=allocok)
        IF (allocok > 0 ) THEN
              IFLAG = -13
              IERROR = NB_BLR_L
              LP = ICNTL(1)
              IF (ICNTL(4) .LE. 0) LP=-1
              IF (LP > 0) WRITE(LP,*) MYID,
     &             ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
              GOTO 700
           ENDIF
#if defined(BLR_MT)          
!$OMP PARALLEL
#endif
        CALL DMUMPS_COMPRESS_PANEL_I_NOOPT
     &        (A_PTR(POSELT), LA_PTR, 1_8,
     &        IFLAG, IERROR, NCOL1,
     &        BEGS_BLR_L(1), size(BEGS_BLR_L), NB_BLR_L+1,
     &        DKEEP(8), KEEP(466), KEEP(473),
     &        BLR_L(1), 
     &        CURRENT_BLR, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCK, MAXI_CLUSTER, NELIM, 
     &        .TRUE.,  
     &        NPIV, NPIV1,
     &        2, KEEP(483), KEEP8, OMP_NUM
     &        )
#if defined(BLR_MT)
!$OMP MASTER
#endif
        IF ( (KEEP(486).EQ.2) 
     &     ) THEN
          CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &         IW(IOLDPS+XXF),
     &         0, 
     &         IPANEL, BLR_L)
        ENDIF
#if defined(BLR_MT)          
!$OMP END MASTER
!$OMP BARRIER
#endif
          IF (IFLAG.LT.0) GOTO 300
          IF (KEEP(475).GE.1) THEN
            IF (DYNAMIC_ALLOC) THEN
              CALL DMUMPS_BLR_PANEL_LRTRSM(
     &              DYN_BLOCFACTO, LA_BLOCFACTO, 1_8, 
     &              LD_BLOCFACTO, -6666, 
     &              NB_BLR_L+1,  
     &              BLR_L, CURRENT_BLR, CURRENT_BLR+1, NB_BLR_L+1, 
     &              2, 0, 0,  
     &              .TRUE.) 
            ELSE
              CALL DMUMPS_BLR_PANEL_LRTRSM(A, LA, POSBLOCFACTO, 
     &              LD_BLOCFACTO, -6666, 
     &              NB_BLR_L+1,  
     &              BLR_L, CURRENT_BLR, CURRENT_BLR+1, NB_BLR_L+1, 
     &              2, 0, 0,  
     &              .TRUE.) 
            ENDIF
#if defined(BLR_MT)          
!$OMP BARRIER
#endif          
            IF (KEEP(486).NE.2) THEN
              CALL DMUMPS_DECOMPRESS_PANEL_I_NOOPT(
     &        A_PTR(POSELT), LA_PTR, 1_8,
     &        NCOL1, NCOL1,
     &        .TRUE.,    
     &        NPIV1+1,   
     &        1,         
     &        NB_BLR_L+1, BLR_L(1), CURRENT_BLR, 'V', 1)
            ENDIF
          ENDIF
 300      CONTINUE
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
        ENDIF
      ENDIF
      IF ( (KEEP(201).eq.1) .AND. 
     &    (OOCWRITE_COMPATIBLE_WITH_BLR .OR. NPIV.EQ.0) ) THEN
        MonBloc%INODE = INODE
        MonBloc%MASTER = .FALSE.
        MonBloc%Typenode = 2
        MonBloc%NROW = NROW1
        MonBloc%NCOL = NCOL1  
        MonBloc%NFS  = NASS1
        MonBloc%LastPiv = NPIV1 + NPIV 
        MonBloc%LastPanelWritten_L = -9999 
        MonBloc%LastPanelWritten_U = -9999 
        NULLIFY(MonBloc%INDICES)
        MonBloc%Last = LASTBL
        STRAT = STRAT_TRY_WRITE 
        NextPivDummy      = -8888 
        LIWFAC = IW(IOLDPS+XXI)
        LAST_CALL = .FALSE.
        CALL DMUMPS_OOC_IO_LU_PANEL_I( STRAT, TYPEF_L,
     &       A_PTR(POSELT),
     &       LA_PTR, MonBloc, NextPivDummy, NextPivDummy,
     &       IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG,LAST_CALL)
      ENDIF
      IF ( (NPIV .GT. 0)
     &   ) THEN
        IF (LR_ACTIVATED) THEN
          IF (NELIM.GT.0) THEN
            UPOS = 1_8+int(NPIV,8)
            IF (DYNAMIC_ALLOC) THEN
              CALL DMUMPS_BLR_UPD_NELIM_VAR_L_I(
     &        DYN_BLOCFACTO, LA_BLOCFACTO, UPOS,
     &        A_PTR(POSELT), LA_PTR, LPOS-POSELT+1_8,
     &        IFLAG, IERROR, LD_BLOCFACTO, NCOL1,
     &        BEGS_BLR_L(1), size(BEGS_BLR_L),
     &        CURRENT_BLR, BLR_L(1), NB_BLR_L+1, 
     &        CURRENT_BLR+1, NELIM, 'N')
            ELSE
              CALL DMUMPS_BLR_UPD_NELIM_VAR_L_I(
     &        A(POSBLOCFACTO), LA_BLOCFACTO, UPOS,
     &        A_PTR(POSELT), LA_PTR, LPOS-POSELT+1_8,
     &        IFLAG, IERROR, LD_BLOCFACTO, NCOL1,
     &        BEGS_BLR_L(1), size(BEGS_BLR_L),
     &        CURRENT_BLR, BLR_L(1), NB_BLR_L+1, 
     &        CURRENT_BLR+1, NELIM, 'N')
            ENDIF
          ENDIF
#if defined(BLR_MT)          
!$OMP PARALLEL
#endif
         CALL DMUMPS_BLR_UPDATE_TRAILING_I(
     &        A_PTR(POSELT), LA_PTR, 1_8, 
     &        IFLAG, IERROR, NCOL1,
     &        BEGS_BLR_L(1), size(BEGS_BLR_L),
     &        BEGS_BLR_U(1), size(BEGS_BLR_U), CURRENT_BLR, 
     &        BLR_L(1), NB_BLR_L+1, 
     &        BLR_U(1), NB_BLR_U+1, 
     &        0,      
     &        .TRUE., 
     &        NPIV1,  
     &        2, 0, 
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &        )
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
        ELSE 
            IF (DYNAMIC_ALLOC) THEN
              UPOS = int(NPIV+1,8)
              CALL dgemm('N','N', NCOL-NPIV, NROW1, NPIV,
     &             ALPHA,DYN_BLOCFACTO(UPOS), NCOL,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS), NCOL1)
            ELSE
              UPOS = POSBLOCFACTO+int(NPIV,8)
              CALL dgemm('N','N', NCOL-NPIV, NROW1, NPIV,
     &             ALPHA,A(UPOS), NCOL,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS), NCOL1)
            ENDIF
        ENDIF
      ENDIF
      IW(IOLDPS+KEEP(IXSZ) ) = IW(IOLDPS+KEEP(IXSZ) ) - NPIV
      IW(IOLDPS + 3+KEEP(IXSZ) ) = IW(IOLDPS+3+KEEP(IXSZ) ) + NPIV
      IF (LASTBL) THEN
        IW(IOLDPS+1+KEEP(IXSZ) ) = IW(IOLDPS + 3+KEEP(IXSZ) )
      ENDIF
      IF ( .not. LASTBL .AND. 
     &   (IW(IOLDPS+1+KEEP(IXSZ)) .EQ. IW(IOLDPS + 3+KEEP(IXSZ))) ) THEN
        write(*,*) 'Internal ERROR 1 **** IN BLACFACTO '
        CALL MUMPS_ABORT()
      ENDIF
      IF (LR_ACTIVATED) THEN
        IF ((NPIV.GT.0)
     &     ) THEN
          CALL DEALLOC_BLR_PANEL( BLR_U, NB_BLR_U, KEEP8)
          DEALLOCATE(BLR_U)
          IF (KEEP(486).EQ.3) THEN
            CALL DEALLOC_BLR_PANEL( BLR_L, NB_BLR_L, KEEP8)
            DEALLOCATE(BLR_L)
          ELSE
            CALL UPD_MRY_LU_LRGAIN(BLR_L, 0, NPARTSCB, 'V')
          ENDIF
      ENDIF 
      ENDIF 
      IF (DYNAMIC_ALLOC) THEN
        DEALLOCATE(DYN_BLOCFACTO)
        DEALLOCATE(DYN_PIVINFO)
      ELSE
        LRLU  = LRLU + LA_BLOCFACTO
        LRLUS = LRLUS + LA_BLOCFACTO
        KEEP8(69) = KEEP8(69) - LA_BLOCFACTO
        POSFAC = POSFAC - LA_BLOCFACTO
        CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &             LA-LRLUS,0_8,-LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
        IWPOS = IWPOS - NPIV
      ENDIF
      FLOP1 = dble( NPIV1*NROW1 ) +
     &        dble(NROW1*NPIV1)*dble(2*NCOL1-NPIV1-1)
     &   -
     &        dble((NPIV1+NPIV)*NROW1 ) -
     &        dble(NROW1*(NPIV1+NPIV))*dble(2*NCOL1-NPIV1-NPIV-1)
      CALL DMUMPS_LOAD_UPDATE( 1, .FALSE., FLOP1, KEEP,KEEP8 )
      IF (LASTBL) THEN
        IF (KEEP(486).NE.0) THEN
          IF (LR_ACTIVATED) THEN
            CALL STATS_COMPUTE_FLOP_SLAVE_TYPE2(NROW1, NCOL1, NASS1,
     &              KEEP(50), INODE)
          ELSE
            CALL UPD_FLOP_FRFRONT_SLAVE(NROW1, NCOL1, NASS1,
     &              KEEP(50), INODE)
          ENDIF
        ENDIF
       IF (LR_ACTIVATED) THEN
         IF (COMPRESS_CB) THEN
           CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_C (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_COL, NPARTSASS_MASTER_AUX)
           BEGS_BLR_COL(1+NPARTSASS_MASTER) = 
     &               BEGS_BLR_COL(1+NPARTSASS_MASTER) - NELIM
           NB_BLR_COL = size(BEGS_BLR_COL) - 1
           IF (NPIV.EQ.0) THEN
             call MAX_CLUSTER(BEGS_BLR_L,NB_BLR_L+1,MAXI_CLUSTER_L)
             call MAX_CLUSTER(BEGS_BLR_COL,NB_BLR_COL,MAXI_CLUSTER_COL)
             IF (COMPRESS_CB) THEN
              MAXI_CLUSTER=max(MAXI_CLUSTER_COL+NELIM,MAXI_CLUSTER_L)
             ELSE
              MAXI_CLUSTER=max(MAXI_CLUSTER_COL,MAXI_CLUSTER_L)
             ENDIF
             LWORK = MAXI_CLUSTER*MAXI_CLUSTER
             OMP_NUM = 1
#if defined(BLR_MT)
!$           OMP_NUM = OMP_GET_MAX_THREADS()
#endif
             ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &         RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &         TAU(MAXI_CLUSTER*OMP_NUM),
     &         JPVT(MAXI_CLUSTER*OMP_NUM), 
     &         WORK(LWORK*OMP_NUM), stat=allocok)
             IF (allocok > 0 ) THEN
               IFLAG = -13
               IERROR = MAXI_CLUSTER*OMP_NUM*MAXI_CLUSTER
     &           + 2*MAXI_CLUSTER*OMP_NUM
     &           + MAXI_CLUSTER*OMP_NUM
     &           + MAXI_CLUSTER*OMP_NUM
     &           + LWORK*OMP_NUM
               LP = ICNTL(1)
               IF (ICNTL(4) .LE. 0) LP=-1
               IF (LP > 0) WRITE(LP,*) MYID,
     &          ': ERROR allocation during DMUMPS_PROCESS_BLOCFACTO'
               GOTO 700
             ENDIF
           ENDIF
           allocate(CB_LRB(NB_BLR_L,NB_BLR_COL-NPARTSASS_MASTER),
     &                 stat=allocok)
           IF (allocok > 0) THEN
             IFLAG  = -13
             IERROR = NB_BLR_L*(NB_BLR_COL-NPARTSASS_MASTER)
             GOTO 700
           ENDIF
           CALL DMUMPS_BLR_SAVE_CB_LRB(IW(IOLDPS+XXF),CB_LRB)
         ENDIF
#if defined(BLR_MT)          
!$OMP PARALLEL
#endif
         IF (COMPRESS_CB) THEN
           CALL DMUMPS_COMPRESS_CB_I(
     &      A_PTR(POSELT), LA_PTR, 1_8, NCOL1,
     &      BEGS_BLR_L(1), size(BEGS_BLR_L),
     &      BEGS_BLR_COL(1), size(BEGS_BLR_COL),
     &      NB_BLR_L, NB_BLR_COL-NPARTSASS_MASTER,
     &      NPARTSASS_MASTER, 
     &      NROW1, NCOL1-NPIV1-NPIV, INODE,
     &      IW(IOLDPS+XXF), 0, 2, IFLAG, IERROR,
     &      DKEEP(12), KEEP(466), KEEP(484), KEEP(489),
     &      CB_LRB(1,1),
     &      WORK, TAU, JPVT, LWORK, RWORK, BLOCK,
     &      MAXI_CLUSTER, KEEP8, OMP_NUM,
     &      -9999, -9999, -9999, KEEP(1)
     &       )
#if defined(BLR_MT)
!$OMP BARRIER
#endif
         ENDIF
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif
         IF (IFLAG.LT.0) GOTO 700
       ENDIF
         CALL DMUMPS_END_FACTO_SLAVE(
     &    COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    root, 
     &    MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     &    PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, PERM,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &     )
      ENDIF 
      IF (LR_ACTIVATED) THEN
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(WORK)) DEALLOCATE(WORK)
        IF (allocated(TAU)) DEALLOCATE(TAU)
        IF (allocated(JPVT)) DEALLOCATE(JPVT)
        IF (allocated(BLOCK)) DEALLOCATE(BLOCK)
        IF (associated(BEGS_BLR_L)) THEN
            IF (.NOT. KEEP_BEGS_BLR_L) DEALLOCATE(BEGS_BLR_L)
        ENDIF
        IF ((NPIV.GT.0)
     &     ) THEN
          IF (associated(BEGS_BLR_U)) DEALLOCATE(BEGS_BLR_U)
        ENDIF
      ENDIF
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_BLOCFACTO
      SUBROUTINE DMUMPS_MPI_UNPACK_LR(
     &           BUFR, LBUFR, LBUFR_BYTES, POSITION,
     &                             NPIV, NELIM, DIR, 
     &                             BLR_U, NB_BLOCK_U,
     &                             BEGS_BLR_U, KEEP8,
     &                             COMM, IERR, IFLAG, IERROR)
      USE DMUMPS_LR_CORE, ONLY : LRB_TYPE, ALLOC_LRB
      USE DMUMPS_LR_TYPE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LBUFR
      INTEGER, INTENT(IN) :: LBUFR_BYTES
      INTEGER, INTENT(IN) :: BUFR(LBUFR)
      INTEGER, INTENT(INOUT) :: POSITION
      INTEGER, INTENT(IN)    :: NB_BLOCK_U, NELIM, NPIV
      CHARACTER(len=1) :: DIR
      INTEGER, INTENT(IN) :: COMM
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, INTENT(OUT) :: IERR
      TYPE (LRB_TYPE), INTENT(OUT), 
     &          DIMENSION(max(NB_BLOCK_U,1)):: BLR_U
      INTEGER, INTENT(OUT), DIMENSION(NB_BLOCK_U+2)  :: BEGS_BLR_U 
      INTEGER(8) :: KEEP8(150)
      LOGICAL :: ISLR
      INTEGER :: ISLR_INT, I
      INTEGER :: K, M, N
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      IERR = 0
      IF (size(BLR_U) .NE. 
     &     MAX(NB_BLOCK_U,1) ) THEN
        WRITE(*,*) "Internal error 1 in DMUMPS_MPI_UNPACK",
     &             NB_BLOCK_U,size(BLR_U)
        CALL MUMPS_ABORT()
      ENDIF
      BEGS_BLR_U(1) = 1
      BEGS_BLR_U(2) = NPIV+NELIM+1 
      DO I = 1, NB_BLOCK_U
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ISLR_INT, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 K, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 M, 1,
     &                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 N, 1,
     &                 MPI_INTEGER, COMM, IERR )
        BEGS_BLR_U(I+2) = BEGS_BLR_U(I+1) + M
        IF (ISLR_INT .eq. 1) THEN
          ISLR = .TRUE.
        ELSE
          ISLR = .FALSE.
        ENDIF
        CALL ALLOC_LRB( BLR_U(I), K, M, N, ISLR, 
     &             IFLAG, IERROR, KEEP8 )
        IF (IFLAG.LT.0) RETURN
        IF (ISLR) THEN
          IF (K .GT. 0) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%Q(1,1), M*K, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR )
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%R(1,1), N*K, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR)
          ENDIF
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     BLR_U(I)%Q(1,1), M*N, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR)
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_MPI_UNPACK_LR
