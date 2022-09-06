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
      SUBROUTINE DMUMPS_ASS_ROOT( root, KEEP50,
     &                     NROW_SON, NCOL_SON, INDROW_SON, 
     &                     INDCOL_SON, NSUPCOL, VAL_SON, VAL_ROOT,
     &                     LOCAL_M, LOCAL_N,
     &                     RHS_ROOT, NLOC_ROOT, CBP )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER, INTENT(IN) :: KEEP50
      INTEGER NCOL_SON, NROW_SON, NSUPCOL
      INTEGER, intent(in) :: CBP   
      INTEGER INDROW_SON( NROW_SON ), INDCOL_SON( NCOL_SON )
      INTEGER LOCAL_M, LOCAL_N
      DOUBLE PRECISION VAL_SON( NCOL_SON, NROW_SON )
      DOUBLE PRECISION VAL_ROOT( LOCAL_M, LOCAL_N )
      INTEGER NLOC_ROOT
      DOUBLE PRECISION RHS_ROOT( LOCAL_M, NLOC_ROOT )
      INTEGER I, J, INDROW, INDCOL, IPOSROOT, JPOSROOT
      IF (CBP .EQ. 0) THEN
        DO I = 1, NROW_SON
          INDROW = INDROW_SON(I)
          IPOSROOT = (root%NPROW*((INDROW-1)/root%MBLOCK)+root%MYROW)
     &             * root%MBLOCK + mod(INDROW-1,root%MBLOCK) + 1
          DO J = 1, NCOL_SON-NSUPCOL
          INDCOL = INDCOL_SON(J)
          IF (KEEP50.NE.0) THEN
            JPOSROOT = (root%NPCOL*((INDCOL-1)/root%NBLOCK)+root%MYCOL)
     &               * root%NBLOCK + mod(INDCOL-1,root%NBLOCK) + 1
            IF (IPOSROOT < JPOSROOT) THEN
              CYCLE
            ENDIF
          ENDIF
          VAL_ROOT( INDROW, INDCOL ) =
     &    VAL_ROOT( INDROW, INDCOL ) + VAL_SON(J,I)
          END DO
          DO J = NCOL_SON-NSUPCOL+1, NCOL_SON
            INDCOL = INDCOL_SON(J)
            RHS_ROOT( INDROW, INDCOL ) =
     &      RHS_ROOT( INDROW, INDCOL ) + VAL_SON(J,I)
          ENDDO
        END DO
      ELSE
        DO I=1, NROW_SON  
          DO J = 1, NCOL_SON  
           RHS_ROOT( INDROW_SON( I ), INDCOL_SON(J)) =
     &     RHS_ROOT(INDROW_SON(I),INDCOL_SON(J)) + VAL_SON(J,I)
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_ASS_ROOT
      RECURSIVE SUBROUTINE DMUMPS_BUILD_AND_SEND_CB_ROOT
     &  ( COMM_LOAD, ASS_IRECV, N, ISON, IROOT,
     &    PTRI, PTRR,
     &    root, NBROW, NBCOL, SHIFT_LIST_ROW_SON,
     &    SHIFT_LIST_COL_SON, SHIFT_VAL_SON_ARG, LDA_ARG, TAG,
     &    MYID, COMM, BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER,
     &    NSTK, COMP, IFLAG, IERROR, PERM,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,TRANSPOSE_ASM,
     &    ND, FRERE,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  
     &               , LRGROUPS
     &     )
      USE DMUMPS_OOC        
      USE DMUMPS_BUF
      USE DMUMPS_LOAD
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR
      IMPLICIT NONE
      INTEGER KEEP(500), ICNTL(60)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION  DKEEP(230)
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, ISON, IROOT, TAG
      INTEGER PTRI( KEEP(28) )
      INTEGER(8) :: PTRR( KEEP(28) )
      INTEGER NBROW, NBCOL
      INTEGER, INTENT(IN):: LDA_ARG
      INTEGER(8), INTENT(IN) :: SHIFT_VAL_SON_ARG
      INTEGER SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON
      INTEGER MYID, COMM
      LOGICAL TRANSPOSE_ASM
      INCLUDE 'mpif.h'
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS, LA
      INTEGER IWPOS, IWPOSCB
      INTEGER LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(N)
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST( KEEP(28) ), PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28)), NSTK( N )
      INTEGER COMP, IFLAG, IERROR
      INTEGER PERM(N)
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER NBFIN, SLAVEF
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER ITLOC( N + KEEP(253) ), FILS( N ), DAD(KEEP(28))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER(8), INTENT(IN) :: PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER INTARR( KEEP8(27) )
      DOUBLE PRECISION DBLARR( KEEP8(26) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION, DIMENSION(:), POINTER :: SONA_PTR
      INTEGER(8) :: LSONA_PTR, POSSONA_PTR
      INTEGER allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTRROW,  PTRCOL
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NSUPROW, NSUPCOL
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ROW_INDEX_LIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: COL_INDEX_LIST
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER I, POS_IN_ROOT, IROW, JCOL, IGLOB, JGLOB
      INTEGER PDEST, IERR
      INTEGER LOCAL_M, LOCAL_N
      INTEGER(8) :: POSROOT
      INTEGER NSUBSET_ROW, NSUBSET_COL
      INTEGER NRLOCAL, NCLOCAL
      INTEGER :: LDA
      INTEGER(8) :: SHIFT_VAL_SON
      LOGICAL SET_IRECV, BLOCKING, MESSAGE_RECEIVED
      INTEGER NBROWS_ALREADY_SENT
      INTEGER SIZE_MSG
      INTEGER LP
      INCLUDE 'mumps_headers.h'
      LOGICAL SKIPLAST_RHS_ROWS, BCP_SYM_NONEMPTY
      INTEGER BBPCBP
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      BBPCBP  = 0   
      LP = ICNTL(1)
      IF ( ICNTL(4) .LE. 0 ) LP = -1
      IF (LDA_ARG < 0) THEN 
        CALL DMUMPS_SET_LDA_SHIFT_VAL_SON(IW, LIW, PTRI(STEP(ISON)),
     &                             LDA, SHIFT_VAL_SON)
      ELSE
        LDA = LDA_ARG
        SHIFT_VAL_SON = SHIFT_VAL_SON_ARG
      ENDIF
      ALLOCATE(PTRROW(root%NPROW + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPROW + 1
      endif
      ALLOCATE(PTRCOL(root%NPCOL + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPCOL + 1
      endif
      ALLOCATE(NSUPROW(root%NPROW + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPROW + 1
      endif
      ALLOCATE(NSUPCOL(root%NPCOL + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPCOL + 1
      endif
      IF (IFLAG.LT.0) THEN
         IF (LP > 0) write(6,*) MYID, ' : MEMORY ALLOCATION ',
     &     'FAILURE in DMUMPS_BUILD_AND_SEND_CB_ROOT'
         CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
         RETURN
      ENDIF
      SKIPLAST_RHS_ROWS = ((KEEP(253).GT.0).AND.(KEEP(50).EQ.0))
      BCP_SYM_NONEMPTY = .FALSE.
      PTRROW = 0
      PTRCOL = 0
      NSUPROW = 0
      NSUPCOL = 0
      DO I = 1, NBROW                   
        IGLOB  =  IW( PTRI(STEP(ISON)) +
     &                          SHIFT_LIST_ROW_SON + I - 1 )
        IF (SKIPLAST_RHS_ROWS.AND.(IGLOB.GT.N)) CYCLE
        IF ( .NOT. TRANSPOSE_ASM ) THEN 
          IF (IGLOB.GT.N) THEN
            BCP_SYM_NONEMPTY = .TRUE.
            POS_IN_ROOT = IGLOB - N
            JCOL =  mod((POS_IN_ROOT-1)/root%NBLOCK,root%NPCOL) 
            NSUPCOL(JCOL+1) = NSUPCOL(JCOL+1) + 1
            PTRCOL( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
          ELSE  
            POS_IN_ROOT = root%RG2L_ROW( IGLOB ) 
            IROW  = mod((POS_IN_ROOT-1)/root%MBLOCK,root%NPROW)
            PTRROW ( IROW + 2 ) = PTRROW( IROW + 2 ) + 1
          ENDIF
        ELSE        
          IF (IGLOB .GT. N) THEN 
            POS_IN_ROOT = IGLOB - N
          ELSE  
            POS_IN_ROOT = root%RG2L_COL( IGLOB ) 
          ENDIF
          JCOL =  mod( ( POS_IN_ROOT - 1 ) / root%NBLOCK, root%NPCOL ) 
          IF (IGLOB.GT.N)  
     &               NSUPCOL(JCOL+1) = NSUPCOL(JCOL+1) + 1
          PTRCOL( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
        END IF
      END DO
      IF (KEEP(50).NE.0 .AND.(.NOT.TRANSPOSE_ASM).AND.BCP_SYM_NONEMPTY)
     &             BBPCBP = 1
      DO I = 1, NBCOL                   
        JGLOB   =  IW( PTRI(STEP(ISON)) +
     &                SHIFT_LIST_COL_SON + I - 1 ) 
        IF ((KEEP(50).GT.0) .AND. (JGLOB.GT.N)) CYCLE  
        IF ( .NOT. TRANSPOSE_ASM ) THEN
          IF (KEEP(50).EQ.0) THEN
            IF (JGLOB.LE.N) THEN
              POS_IN_ROOT = root%RG2L_COL(JGLOB)  
            ELSE
              POS_IN_ROOT = JGLOB-N
            ENDIF
            JCOL =  mod((POS_IN_ROOT-1) / root%NBLOCK, root%NPCOL ) 
            IF (JGLOB.GT.N) THEN
             NSUPCOL(JCOL+1) = NSUPCOL(JCOL+1)  + 1  
            ENDIF
            PTRCOL ( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
          ELSE 
            POS_IN_ROOT = root%RG2L_COL(JGLOB) 
            JCOL =  mod((POS_IN_ROOT-1) / root%NBLOCK, root%NPCOL )
            PTRCOL ( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
            IF (BCP_SYM_NONEMPTY) THEN
             POS_IN_ROOT = root%RG2L_ROW(JGLOB) 
             IROW  = mod((POS_IN_ROOT-1)/root%MBLOCK,root%NPROW)
             NSUPROW(IROW+1) = NSUPROW(IROW+1)+1
             PTRROW( IROW + 2 ) = PTRROW( IROW + 2 ) + 1
            ENDIF
          ENDIF
        ELSE  
          IF (JGLOB.LE.N) THEN
           POS_IN_ROOT = root%RG2L_ROW( JGLOB ) 
          ELSE
           POS_IN_ROOT = JGLOB-N
          ENDIF
          IROW        = mod( ( POS_IN_ROOT - 1 ) /
     &                  root%MBLOCK, root%NPROW )
          PTRROW ( IROW + 2 ) = PTRROW( IROW + 2 ) + 1
        END IF
      END DO
      PTRROW( 1 ) = 1
      DO IROW = 2, root%NPROW + 1
        PTRROW( IROW ) = PTRROW( IROW ) + PTRROW( IROW - 1 )
      END DO
      PTRCOL( 1 ) = 1
      DO JCOL = 2, root%NPCOL + 1
        PTRCOL( JCOL ) = PTRCOL( JCOL ) + PTRCOL( JCOL - 1 )
      END DO
      ALLOCATE(ROW_INDEX_LIST(PTRROW(root%NPROW+1)-1+1),
     &         stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = PTRROW(root%NPROW+1)-1+1
      endif
      ALLOCATE(COL_INDEX_LIST(PTRCOL(root%NPCOL+1)-1+1),
     &         stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = PTRCOL(root%NPCOL+1)-1+1
      endif
      DO I = 1, NBROW
        IGLOB  =  IW( PTRI(STEP(ISON)) +
     &                          SHIFT_LIST_ROW_SON + I - 1 )
        IF (SKIPLAST_RHS_ROWS.AND.(IGLOB.GT.N)) CYCLE
        IF ( .NOT. TRANSPOSE_ASM ) THEN
          IF (IGLOB.GT.N) CYCLE   
          POS_IN_ROOT = root%RG2L_ROW( IGLOB )
          IROW        = mod( ( POS_IN_ROOT - 1 ) / root%MBLOCK,
     &                       root%NPROW )
          ROW_INDEX_LIST( PTRROW( IROW + 1 ) ) = I 
          PTRROW ( IROW + 1 ) = PTRROW( IROW + 1 ) + 1
        ELSE
          IF (IGLOB.LE.N) THEN
           POS_IN_ROOT = root%RG2L_COL( IGLOB )
          ELSE
           POS_IN_ROOT = IGLOB - N  
          ENDIF
          JCOL        = mod( ( POS_IN_ROOT - 1 ) / root%NBLOCK,
     &                       root%NPCOL )
          COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ) = I 
          PTRCOL ( JCOL + 1 ) = PTRCOL( JCOL + 1 ) + 1 
        END IF
      END DO
      DO I = 1, NBCOL 
        JGLOB =  IW( PTRI(STEP(ISON))+SHIFT_LIST_COL_SON+I - 1 ) 
        IF ((KEEP(50).GT.0) .AND. (JGLOB.GT.N)) CYCLE  
        IF ( .NOT. TRANSPOSE_ASM ) THEN
          IF ( JGLOB.LE.N ) THEN
           POS_IN_ROOT = root%RG2L_COL( JGLOB )
          ELSE
           POS_IN_ROOT = JGLOB - N
          ENDIF
          JCOL        = mod( ( POS_IN_ROOT - 1 ) /
     &               root%NBLOCK, root%NPCOL )
          COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ) = I 
          PTRCOL ( JCOL + 1 ) = PTRCOL( JCOL + 1 ) + 1
        ELSE
          IF ( JGLOB.LE.N ) THEN
           POS_IN_ROOT = root%RG2L_ROW( JGLOB )
          ELSE
           POS_IN_ROOT = JGLOB - N
          ENDIF
          IROW        = mod( ( POS_IN_ROOT - 1 ) /
     &                root%MBLOCK, root%NPROW )
          ROW_INDEX_LIST( PTRROW( IROW + 1 ) ) = I    
          PTRROW( IROW + 1 ) = PTRROW( IROW + 1 ) + 1 
        END IF
      END DO
      IF (BCP_SYM_NONEMPTY) THEN
        DO I = 1, NBROW
          IGLOB  =  IW( PTRI(STEP(ISON)) +
     &                         SHIFT_LIST_ROW_SON + I - 1 )
          IF (IGLOB.LE.N) CYCLE  
          POS_IN_ROOT = IGLOB - N
          JCOL =  mod((POS_IN_ROOT-1)/root%NBLOCK,root%NPCOL)
          COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ) = I 
          PTRCOL ( JCOL + 1 ) = PTRCOL( JCOL + 1 ) + 1
        ENDDO
        DO I=1, NBCOL
         JGLOB =  IW( PTRI(STEP(ISON))+SHIFT_LIST_COL_SON+I - 1 ) 
         IF (JGLOB.GT.N) THEN
           EXIT
         ELSE
           POS_IN_ROOT = root%RG2L_ROW(JGLOB) 
         ENDIF
         IROW  = mod((POS_IN_ROOT-1)/root%MBLOCK,root%NPROW)
         ROW_INDEX_LIST( PTRROW( IROW + 1 ) ) = I    
         PTRROW( IROW + 1 ) = PTRROW( IROW + 1 ) + 1 
        ENDDO
      ENDIF
      DO IROW = root%NPROW, 2, -1
        PTRROW( IROW ) = PTRROW( IROW - 1 )
      END DO
      PTRROW( 1 ) = 1
      DO JCOL = root%NPCOL, 2, -1
        PTRCOL( JCOL ) = PTRCOL( JCOL - 1 )
      END DO
      PTRCOL( 1 ) = 1
      JCOL  = root%MYCOL
      IROW  = root%MYROW
      IF ( root%yes ) THEN
         if (IROW .ne. root%MYROW .or. JCOL.ne.root%MYCOL) then
        write(*,*) ' error in grid position buildandsendcbroot'
        CALL MUMPS_ABORT()
        end if
        IF ( PTRIST(STEP(IROOT)).EQ.0.AND.
     &       PTLUST_S(STEP(IROOT)).EQ.0) THEN
           CALL DMUMPS_ROOT_ALLOC_STATIC(root, IROOT, N, IW, LIW,
     &               A, LA,
     &               FILS, DAD, MYID, SLAVEF, PROCNODE_STEPS,
     &               LPTRAR, NELT, FRTPTR, FRTELT,
     &               PTRAIW, PTRARW, INTARR, DBLARR,
     &               LRLU, IPTRLU,
     &               IWPOS, IWPOSCB, PTRIST, PTRAST,
     &               STEP, PIMASTER, PAMASTER, ITLOC, RHS_MUMPS,
     &               COMP, LRLUS, IFLAG, KEEP,KEEP8,DKEEP, IERROR )
           KEEP(121) = -1 
           IF (IFLAG.LT.0) THEN
                CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
                RETURN
           ENDIF
        ELSE
           KEEP(121) = KEEP(121) - 1
           IF ( KEEP(121) .eq. 0 ) THEN 
              IF (KEEP(201).EQ.1) THEN
                 CALL DMUMPS_OOC_FORCE_WRT_BUF_PANEL(IERR)
              ELSE IF (KEEP(201).EQ.2) THEN
                 CALL DMUMPS_FORCE_WRITE_BUF(IERR)              
              ENDIF
              CALL DMUMPS_INSERT_POOL_N(N, IPOOL, LPOOL, PROCNODE_STEPS,
     &        SLAVEF, KEEP(199), KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     &        STEP, IROOT+N )
              IF (KEEP(47) .GE. 3) THEN
                 CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &                IPOOL, LPOOL, 
     &                PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &                MYID, STEP, N, ND, FILS )
              ENDIF
          END IF
        END IF
        CALL DMUMPS_DM_SET_DYNPTR( IW(PTRI(STEP(ISON))+XXS), A, LA,
     &      PTRR(STEP(ISON)), IW(PTRI(STEP(ISON))+XXD),
     &      IW(PTRI(STEP(ISON))+XXR),
     &      SONA_PTR, POSSONA_PTR, LSONA_PTR )
       IF (KEEP(60) .NE. 0 ) THEN
         LOCAL_M = root%SCHUR_LLD
         LOCAL_N = root%SCHUR_NLOC
            NRLOCAL = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
            NCLOCAL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
            CALL DMUMPS_ROOT_LOCAL_ASSEMBLY( N,
     &        root%SCHUR_POINTER(1),
     &        LOCAL_M, LOCAL_N,
     &        root%NPCOL, root%NPROW, root%MBLOCK, root%NBLOCK,
     &        NBCOL, NBROW,
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     &        LDA, SONA_PTR( POSSONA_PTR + SHIFT_VAL_SON ),
     &        ROW_INDEX_LIST( PTRROW( IROW + 1 ) ),
     &        COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     &        NRLOCAL,
     &        NCLOCAL,
     &        NSUPROW(IROW+1), NSUPCOL(JCOL+1),
     &        root%RG2L_ROW(1), root%RG2L_COL(1), TRANSPOSE_ASM,
     &        KEEP,
     &        root%RHS_ROOT(1,1), root%RHS_NLOC )
       ELSE
        IF ( PTRIST(STEP( IROOT )) .GE. 0 ) THEN
          IF ( PTRIST(STEP( IROOT )) .EQ. 0 ) THEN
            LOCAL_N = IW( PTLUST_S(STEP(IROOT)) + 1 + KEEP(IXSZ))
            LOCAL_M = IW( PTLUST_S(STEP(IROOT)) + 2 + KEEP(IXSZ))
            POSROOT = PTRFAC(IW( PTLUST_S(STEP(IROOT)) +4+KEEP(IXSZ) ))
          ELSE
            LOCAL_N = - IW( PTRIST(STEP(IROOT)) +KEEP(IXSZ))
            LOCAL_M = IW( PTRIST(STEP(IROOT)) + 1 +KEEP(IXSZ))
            POSROOT = PAMASTER(STEP( IROOT ))
          ENDIF
          NCLOCAL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
          NRLOCAL = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
          CALL DMUMPS_ROOT_LOCAL_ASSEMBLY( N, A( POSROOT ),
     &        LOCAL_M, LOCAL_N,
     &        root%NPCOL, root%NPROW, root%MBLOCK, root%NBLOCK,
     &        NBCOL, NBROW,
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     &        LDA, SONA_PTR( POSSONA_PTR + SHIFT_VAL_SON ),
     &        ROW_INDEX_LIST( PTRROW( IROW + 1 ) ),
     &        COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     &        NRLOCAL,
     &        NCLOCAL,
     &        NSUPROW(IROW+1), NSUPCOL(JCOL+1),
     &        root%RG2L_ROW(1), root%RG2L_COL(1), TRANSPOSE_ASM,
     &        KEEP,
     &        root%RHS_ROOT(1,1), root%RHS_NLOC )
        END IF
       ENDIF
      END IF
      DO IROW = 0, root%NPROW - 1
        DO JCOL = 0, root%NPCOL - 1
          PDEST = IROW * root%NPCOL + JCOL
          IF ( (root%MYROW.eq.IROW.and.root%MYCOL.eq.JCOL) .and.
     &         MYID.ne.PDEST) THEN
            write(*,*) 'error: myrow,mycol=',root%MYROW,root%MYCOL
            write(*,*) ' MYID,PDEST=',MYID,PDEST
            CALL MUMPS_ABORT()
          END IF
          IF ( root%MYROW .NE. IROW .OR. root%MYCOL .NE. JCOL) THEN
            NBROWS_ALREADY_SENT = 0
            IERR = -1
            DO WHILE ( IERR .EQ. -1 )
              NSUBSET_ROW = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
              NSUBSET_COL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
              IF ( LRLU .LT. int(NSUBSET_ROW,8) * int(NSUBSET_COL,8)
     &        .AND. LRLUS .GT. int(NSUBSET_ROW,8) * int(NSUBSET_COL,8) )
     &        THEN
                CALL DMUMPS_COMPRE_NEW(N, KEEP(28),
     &          IW, LIW, A, LA,
     &          LRLU, IPTRLU,
     &          IWPOS, IWPOSCB, PTRIST, PTRAST,
     &          STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &          KEEP(IXSZ), COMP, DKEEP(97),
     &          MYID, SLAVEF, KEEP(199), PROCNODE_STEPS, DAD)
                IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) MYID,": pb compress in",
     &                            "DMUMPS_BUILD_AND_SEND_CB_ROOT"
                  WRITE(*,*) MYID,': LRLU, LRLUS=',LRLU,LRLUS
                  CALL MUMPS_ABORT()
                END IF
              END IF
              CALL DMUMPS_DM_SET_DYNPTR(
     &        IW(PTRI(STEP(ISON))+XXS), A, LA,
     &        PTRR(STEP(ISON)), IW(PTRI(STEP(ISON))+XXD),
     &        IW(PTRI(STEP(ISON))+XXR),
     &        SONA_PTR, POSSONA_PTR, LSONA_PTR )
              CALL DMUMPS_BUF_SEND_CONTRIB_TYPE3_I( N, ISON,
     &        NBCOL, NBROW,
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     &        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     &        LDA, SONA_PTR( POSSONA_PTR + SHIFT_VAL_SON ),
     &        TAG,
     &        ROW_INDEX_LIST( PTRROW( IROW + 1 ) ),
     &        COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     &        NSUBSET_ROW, NSUBSET_COL,
     &        NSUPROW(IROW+1), NSUPCOL(JCOL+1),
     &        root%NPROW, root%NPCOL, root%MBLOCK,
     &        root%RG2L_ROW(1), root%RG2L_COL(1),
     &        root%NBLOCK, PDEST,
     &        COMM, IERR, A( POSFAC ), LRLU, TRANSPOSE_ASM,
     &        SIZE_MSG, NBROWS_ALREADY_SENT, KEEP, BBPCBP )
              IF ( IERR .EQ. -1 ) THEN
                  BLOCKING  = .FALSE.
                  SET_IRECV = .TRUE.
                  MESSAGE_RECEIVED = .FALSE.
                  CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV, 
     &            BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &            MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &            STATUS, BUFR, LBUFR,
     &            LBUFR_BYTES, PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB,
     &            IPTRLU, LRLU, LRLUS, N, IW, LIW, A, LA,
     &            PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     &            PIMASTER, PAMASTER, NSTK,
     &            COMP, IFLAG, IERROR, COMM, PERM, IPOOL, LPOOL,
     &            LEAF, NBFIN, MYID, SLAVEF, root,
     &            OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD,
     &            PTRARW,PTRAIW,INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,
     &            ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT, 
     &            ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &             )
                  IF ( IFLAG .LT. 0 ) GOTO 500
                  IF (LDA_ARG < 0) THEN 
                    CALL DMUMPS_SET_LDA_SHIFT_VAL_SON(
     &                    IW, LIW, PTRI(STEP(ISON)),
     &                    LDA, SHIFT_VAL_SON)
                  ENDIF
              END IF
            END DO
            IF ( IERR == -2 ) THEN
              IFLAG  = -17
              IERROR = SIZE_MSG
              IF (LP > 0) WRITE(LP, *) "FAILURE, SEND BUFFER TOO
     & SMALL DURING DMUMPS_BUILD_AND_SEND_CB_ROOT"
              CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
              GOTO 500
            ENDIF
            IF ( IERR == -3 ) THEN
              IF (LP > 0) WRITE(LP, *) "FAILURE, RECV BUFFER TOO
     & SMALL DURING DMUMPS_BUILD_AND_SEND_CB_ROOT"
              IFLAG  = -20
              IERROR = SIZE_MSG
              CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
              GOTO 500
            ENDIF
          END IF
        END DO
      END DO
 500  CONTINUE
      DEALLOCATE(PTRROW)
      DEALLOCATE(PTRCOL)
      DEALLOCATE(ROW_INDEX_LIST)
      DEALLOCATE(COL_INDEX_LIST)
      RETURN
      CONTAINS
        SUBROUTINE DMUMPS_SET_LDA_SHIFT_VAL_SON(IW, LIW, IOLDPS,
     &                             LDA, SHIFT_VAL_SON)
        INTEGER,    INTENT(IN)  :: LIW, IOLDPS
        INTEGER,    INTENT(IN)  :: IW(LIW)
        INTEGER,    INTENT(OUT) :: LDA
        INTEGER(8), INTENT(OUT) :: SHIFT_VAL_SON
        INCLUDE 'mumps_headers.h'
        INTEGER :: LCONT, NROW, NPIV, NASS, NELIM
        LCONT  = IW(IOLDPS+KEEP(IXSZ))
        NROW   = IW(IOLDPS+2+KEEP(IXSZ))
        NPIV   = IW(IOLDPS+3+KEEP(IXSZ))
        NASS   = IW(IOLDPS+4+KEEP(IXSZ))
        NELIM  = NASS-NPIV
        IF (IW(IOLDPS+XXS).EQ.S_NOLCBNOCONTIG38.OR.
     &      IW(IOLDPS+XXS).EQ.S_ALL) THEN
          SHIFT_VAL_SON      = int(NPIV,8)
          LDA                = LCONT + NPIV
        ELSE IF (IW(IOLDPS+XXS).EQ.S_NOLCBCONTIG38) THEN
          SHIFT_VAL_SON = int(NROW,8)*int(LCONT+NPIV-NELIM,8)
          LDA           = NELIM
        ELSE IF (IW(IOLDPS+XXS).EQ.S_NOLCLEANED38) THEN
          SHIFT_VAL_SON=0_8
          LDA = NELIM
       ELSE
          WRITE(*,*) MYID,
     &    ": internal error in DMUMPS_SET_LDA_SHIFT_VAL_SON",
     &    IW(IOLDPS+XXS), "ISON=",ISON
          CALL MUMPS_ABORT()
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_SET_LDA_SHIFT_VAL_SON
      END SUBROUTINE DMUMPS_BUILD_AND_SEND_CB_ROOT
      SUBROUTINE DMUMPS_ROOT_LOCAL_ASSEMBLY( N, VAL_ROOT,
     &   LOCAL_M, LOCAL_N,
     &   NPCOL, NPROW, MBLOCK, NBLOCK, NBCOL_SON, NBROW_SON, INDCOL_SON,
     &   INDROW_SON, LD_SON, VAL_SON, SUBSET_ROW, SUBSET_COL,
     &   NSUBSET_ROW, NSUBSET_COL, NSUPROW, NSUPCOL,
     &   RG2L_ROW, RG2L_COL, TRANSPOSE_ASM,
     &   KEEP, RHS_ROOT, NLOC  )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INTEGER N, LOCAL_M, LOCAL_N
      DOUBLE PRECISION VAL_ROOT( LOCAL_M, LOCAL_N )
      INTEGER NPCOL, NPROW, MBLOCK, NBLOCK
      INTEGER NBCOL_SON, NBROW_SON
      INTEGER INDCOL_SON( NBCOL_SON ), INDROW_SON( NBROW_SON )
      INTEGER LD_SON 
      INTEGER NSUPROW, NSUPCOL
      DOUBLE PRECISION VAL_SON( LD_SON, NBROW_SON )
      INTEGER KEEP(500)
      INTEGER NSUBSET_ROW, NSUBSET_COL
      INTEGER SUBSET_ROW( NSUBSET_ROW ), SUBSET_COL( NSUBSET_COL )
      INTEGER RG2L_ROW( N ), RG2L_COL( N )
      LOGICAL TRANSPOSE_ASM
      INTEGER NLOC
      DOUBLE PRECISION RHS_ROOT( LOCAL_M, NLOC)
      INTEGER ISUB, JSUB, I, J, IPOS_ROOT, JPOS_ROOT
      INTEGER ILOC_ROOT, JLOC_ROOT, IGLOB, JGLOB
      IF (KEEP(50).EQ.0) THEN
        DO ISUB = 1, NSUBSET_ROW
          I         = SUBSET_ROW( ISUB )
          IGLOB     = INDROW_SON( I )
          IPOS_ROOT = RG2L_ROW( IGLOB )   
          ILOC_ROOT = MBLOCK
     &            * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &            + mod( IPOS_ROOT - 1, MBLOCK ) + 1
          DO JSUB = 1, NSUBSET_COL-NSUPCOL
            J         = SUBSET_COL( JSUB )
            JGLOB     = INDCOL_SON( J )
            JPOS_ROOT = RG2L_COL( JGLOB )
            JLOC_ROOT = NBLOCK
     &              * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &              + mod( JPOS_ROOT - 1, NBLOCK ) + 1
            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     &           VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
          END DO
          DO JSUB = NSUBSET_COL-NSUPCOL+1, NSUBSET_COL
            J         = SUBSET_COL( JSUB )
            JGLOB     = INDCOL_SON( J )
             JPOS_ROOT = JGLOB - N  
             JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
             RHS_ROOT(ILOC_ROOT, JLOC_ROOT) =  
     &            RHS_ROOT(ILOC_ROOT, JLOC_ROOT) + VAL_SON( J, I )
          ENDDO
        END DO
      ELSE
        IF ( .NOT. TRANSPOSE_ASM ) THEN
          DO ISUB = 1, NSUBSET_ROW - NSUPROW 
            I         = SUBSET_ROW( ISUB )
            IGLOB     = INDROW_SON( I )
            IPOS_ROOT = RG2L_ROW( IGLOB )
            ILOC_ROOT = MBLOCK
     &            * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &            + mod( IPOS_ROOT - 1, MBLOCK ) + 1
            DO JSUB = 1, NSUBSET_COL -NSUPCOL
              J         = SUBSET_COL( JSUB )
              JGLOB     = INDCOL_SON( J )
              JPOS_ROOT = RG2L_COL( JGLOB )
              IF (KEEP(50).NE.0. AND. JPOS_ROOT .GT. IPOS_ROOT) CYCLE
              JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     &            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
            END DO
          END DO
          DO JSUB = NSUBSET_COL -NSUPCOL+1, NSUBSET_COL
            J         = SUBSET_COL( JSUB )
            JGLOB     = INDROW_SON( J )  
            JPOS_ROOT = JGLOB - N  
            JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
            DO ISUB = NSUBSET_ROW - NSUPROW +1, NSUBSET_ROW
              I         = SUBSET_ROW( ISUB )
              IGLOB     = INDCOL_SON( I )  
              IPOS_ROOT = RG2L_ROW(IGLOB)
              ILOC_ROOT = MBLOCK
     &            * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &            + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              RHS_ROOT(ILOC_ROOT, JLOC_ROOT) =  
     &            RHS_ROOT(ILOC_ROOT, JLOC_ROOT) + VAL_SON( I, J )
            END DO
          END DO
        ELSE
          DO ISUB = 1, NSUBSET_COL-NSUPCOL 
            I         = SUBSET_COL( ISUB )
            IGLOB     = INDROW_SON( I )
            JPOS_ROOT = RG2L_COL( IGLOB )
            JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
            DO JSUB = 1, NSUBSET_ROW
              J         = SUBSET_ROW( JSUB )
              JGLOB     = INDCOL_SON( J )
              IPOS_ROOT = RG2L_ROW( JGLOB )  
              ILOC_ROOT = MBLOCK
     &                * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     &            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
            END DO
           ENDDO
           DO ISUB = NSUBSET_COL-NSUPCOL+1, NSUBSET_COL
            I         = SUBSET_COL( ISUB )
            IGLOB     = INDROW_SON( I )
            JPOS_ROOT = IGLOB - N 
            JLOC_ROOT = NBLOCK
     &                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                + mod( JPOS_ROOT - 1, NBLOCK ) + 1
            DO JSUB = 1, NSUBSET_ROW
              J         = SUBSET_ROW( JSUB )
              JGLOB     = INDCOL_SON( J )
              IPOS_ROOT = RG2L_ROW( JGLOB )  
              ILOC_ROOT = MBLOCK
     &                * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              RHS_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     &            RHS_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
            END DO
           ENDDO
        END IF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_ROOT_LOCAL_ASSEMBLY
      SUBROUTINE DMUMPS_INIT_ROOT_ANA
     &( MYID, NPROCS, N, root, COMM_ROOT, IROOT, FILS,
     &  K50, K46, K51
     &     , K60, IDNPROW, IDNPCOL, IDMBLOCK, IDNBLOCK
     & )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INTEGER MYID, MYID_ROOT
      TYPE (DMUMPS_ROOT_STRUC)::root
      INTEGER COMM_ROOT
      INTEGER N, IROOT, NPROCS, K50, K46, K51
      INTEGER FILS( N )
      INTEGER K60, IDNPROW, IDNPCOL, IDMBLOCK, IDNBLOCK
      INTEGER INODE, NPROWtemp, NPCOLtemp
      LOGICAL SLAVE
      root%ROOT_SIZE     = 0
      root%TOT_ROOT_SIZE = 0
      SLAVE = ( MYID .ne. 0 .or.
     &        ( MYID .eq. 0 .and. K46 .eq. 1 ) )
      INODE = IROOT
      DO WHILE ( INODE .GT. 0 )
        INODE = FILS( INODE )
        root%ROOT_SIZE = root%ROOT_SIZE + 1
      END DO
      IF ( ( K60 .NE. 2 .AND. K60 .NE. 3 ) .OR.
     &       IDNPROW .LE. 0 .OR. IDNPCOL .LE. 0
     &      .OR. IDMBLOCK .LE.0 .OR. IDNBLOCK.LE.0
     &      .OR. IDNPROW * IDNPCOL .GT. NPROCS ) THEN
        root%MBLOCK = K51
        root%NBLOCK = K51
        CALL DMUMPS_DEF_GRID( NPROCS, root%NPROW, root%NPCOL,
     &                         root%ROOT_SIZE, K50 )
        IF  ( K60 .EQ. 2 .OR. K60 .EQ. 3 ) THEN
          IDNPROW = root%NPROW
          IDNPCOL = root%NPCOL
          IDMBLOCK = root%MBLOCK
          IDNBLOCK = root%NBLOCK
        ENDIF
      ELSE IF  ( K60 .EQ. 2 .OR. K60 .EQ. 3 ) THEN
        root%NPROW = IDNPROW
        root%NPCOL = IDNPCOL
        root%MBLOCK = IDMBLOCK
        root%NBLOCK = IDNBLOCK
      ENDIF
      IF  ( K60 .EQ. 2 .OR. K60 .EQ. 3 ) THEN
        IF (SLAVE) THEN
          root%LPIV = 0
          IF (K46.EQ.0) THEN
            MYID_ROOT=MYID-1
          ELSE
            MYID_ROOT=MYID
          ENDIF
          IF (MYID_ROOT < root%NPROW*root%NPCOL) THEN
            root%MYROW = MYID_ROOT / root%NPCOL
            root%MYCOL = mod(MYID_ROOT, root%NPCOL)
            root%yes  = .true.
          ELSE
            root%MYROW = -1
            root%MYCOL = -1
            root%yes  = .FALSE.
          ENDIF
        ELSE
          root%yes  = .FALSE.
        ENDIF
      ELSE IF ( SLAVE ) THEN
        IF ( root%gridinit_done) THEN
           IF (root%yes) THEN
             CALL blacs_gridexit( root%CNTXT_BLACS )
             root%gridinit_done = .FALSE.
           ENDIF
        END IF
        root%CNTXT_BLACS = COMM_ROOT
        CALL blacs_gridinit( root%CNTXT_BLACS, 'R',
     &                       root%NPROW, root%NPCOL )
        root%gridinit_done = .TRUE.
        CALL blacs_gridinfo( root%CNTXT_BLACS,
     &                       NPROWtemp, NPCOLtemp,
     &                       root%MYROW, root%MYCOL )
        IF ( root%MYROW .NE. -1 ) THEN
          root%yes = .true.
        ELSE
          root%yes = .false.
        END IF
        root%LPIV = 0
      ELSE
        root%yes = .FALSE.
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_INIT_ROOT_ANA
      SUBROUTINE DMUMPS_INIT_ROOT_FAC( N, root, FILS, IROOT,
     &                                 KEEP, INFO )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE ( DMUMPS_ROOT_STRUC ):: root
      INTEGER N, IROOT, INFO(80), KEEP(500)
      INTEGER FILS( N )
      INTEGER INODE, I, allocok
      IF ( associated( root%RG2L_ROW ) ) THEN
        DEALLOCATE( root%RG2L_ROW )
        NULLIFY( root%RG2L_ROW )
      ENDIF
      IF ( associated( root%RG2L_COL ) ) THEN
        DEALLOCATE( root%RG2L_COL )
        NULLIFY( root%RG2L_COL )
      ENDIF
      ALLOCATE( root%RG2L_ROW( N ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        RETURN
      ENDIF
      ALLOCATE( root%RG2L_COL( N ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        DEALLOCATE( root%RG2L_ROW ); NULLIFY( root%RG2L_ROW )
        INFO(1)=-13
        INFO(2)=N
        RETURN
      ENDIF
      INODE = IROOT
      I = 1
      DO WHILE ( INODE .GT. 0 )
        root%RG2L_ROW( INODE ) = I
        root%RG2L_COL( INODE ) = I
        I = I + 1
        INODE = FILS( INODE )
      END DO
      root%TOT_ROOT_SIZE=0
      RETURN
      END SUBROUTINE DMUMPS_INIT_ROOT_FAC
      SUBROUTINE DMUMPS_DEF_GRID( NPROCS, NPROW, NPCOL, SIZE, K50 )
      IMPLICIT NONE
      INTEGER NPROCS, NPROW, NPCOL, SIZE, K50
      INTEGER NPROWtemp, NPCOLtemp, NPROCSused, FLATNESS
      LOGICAL KEEPIT
      IF ( K50 .EQ. 1 ) THEN
        FLATNESS = 2
      ELSE
        FLATNESS = 3
      ENDIF
      NPROW  = int(sqrt(dble(NPROCS)))
      NPROWtemp = NPROW
      NPCOL  = int(NPROCS / NPROW)
      NPCOLtemp = NPCOL
      NPROCSused = NPROWtemp * NPCOLtemp
 10   CONTINUE
      IF ( NPROWtemp >= NPCOLtemp/FLATNESS .AND. NPROWtemp > 1) THEN
        NPROWtemp = NPROWtemp - 1
        NPCOLtemp = int(NPROCS / NPROWtemp)
        KEEPIT=.FALSE.
        IF ( NPROWtemp * NPCOLtemp .GE. NPROCSused ) THEN
          IF ( ( K50 .NE. 1 .AND. NPROWtemp >= NPCOLtemp/FLATNESS)
     &         .OR. NPROWtemp * NPCOLtemp .GT. NPROCSused )
     &         KEEPIT=.TRUE.
        END IF
        IF ( KEEPIT ) THEN
          NPROW = NPROWtemp
          NPCOL = NPCOLtemp
          NPROCSused = NPROW * NPCOL
        END IF
        GO TO 10
      END IF
      RETURN
      END SUBROUTINE DMUMPS_DEF_GRID
      SUBROUTINE DMUMPS_SCATTER_ROOT(MYID, M, N, ASEQ,
     &                    LOCAL_M, LOCAL_N,
     &                    MBLOCK, NBLOCK,
     &                    APAR,
     &                    MASTER_ROOT,
     &                    NPROW, NPCOL,
     &                    COMM)
      IMPLICIT NONE
      INTEGER MYID, MASTER_ROOT, COMM
      INTEGER M, N
      INTEGER NPROW, NPCOL
      INTEGER LOCAL_M, LOCAL_N
      INTEGER MBLOCK, NBLOCK
      DOUBLE PRECISION APAR( LOCAL_M, LOCAL_N )
      DOUBLE PRECISION ASEQ( M, N )
      INCLUDE 'mpif.h'
      INTEGER I, J, SIZE_IBLOCK, SIZE_JBLOCK, IDEST, IROW, ICOL
      INTEGER IBLOCK, JBLOCK, II, JJ, KK
      INTEGER IAPAR, JAPAR, IERR, allocok
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WK
      LOGICAL JUPDATE
      ALLOCATE(WK( MBLOCK * NBLOCK ), stat=allocok)
      if(allocok.ne.0) then
         WRITE(6,*) ' Allocation error of WK in '
     &       // 'routine DMUMPS_SCATTER_ROOT '
         CALL MUMPS_ABORT()
      endif
        IAPAR = 1
        JAPAR = 1
        DO J = 1, N, NBLOCK
          SIZE_JBLOCK = NBLOCK
          IF ( J + NBLOCK > N ) THEN
            SIZE_JBLOCK = N - J + 1
          END IF
          JUPDATE = .FALSE.
          DO I = 1, M, MBLOCK
            SIZE_IBLOCK = MBLOCK
            IF ( I + MBLOCK > M ) THEN
              SIZE_IBLOCK = M - I + 1
            END IF
            IBLOCK = I / MBLOCK
            JBLOCK = J / NBLOCK
            IROW = mod ( IBLOCK, NPROW )
            ICOL = mod ( JBLOCK, NPCOL )
            IDEST = IROW * NPCOL + ICOL
            IF ( IDEST .NE. MASTER_ROOT ) THEN
              IF ( MYID .EQ. MASTER_ROOT ) THEN
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  WK(KK)=ASEQ(II,JJ)
                  KK=KK+1
                END DO
                END DO
                CALL MPI_SSEND( WK, SIZE_IBLOCK*SIZE_JBLOCK,
     &                         MPI_DOUBLE_PRECISION,
     &                         IDEST, 128, COMM, IERR )
              ELSE IF ( MYID .EQ. IDEST ) THEN
                CALL MPI_RECV( WK(1),
     &                         SIZE_IBLOCK*SIZE_JBLOCK,
     &                         MPI_DOUBLE_PRECISION,
     &                         MASTER_ROOT,128,COMM,STATUS,IERR)
                KK=1
                DO JJ=JAPAR,JAPAR+SIZE_JBLOCK-1
                DO II=IAPAR,IAPAR+SIZE_IBLOCK-1
                  APAR(II,JJ)=WK(KK)
                  KK=KK+1
                END DO
                END DO
                JUPDATE = .TRUE.
                IAPAR = IAPAR + SIZE_IBLOCK
              END IF
            ELSE IF ( MYID.EQ. MASTER_ROOT ) THEN
              APAR( IAPAR:IAPAR+SIZE_IBLOCK-1,
     &              JAPAR:JAPAR+SIZE_JBLOCK-1 )
     &        = ASEQ(I:I+SIZE_IBLOCK-1,J:J+SIZE_JBLOCK-1)
              JUPDATE = .TRUE.
              IAPAR = IAPAR + SIZE_IBLOCK
            END IF
          END DO
          IF ( JUPDATE ) THEN
            IAPAR = 1
            JAPAR = JAPAR + SIZE_JBLOCK
          END IF
        END DO
        DEALLOCATE(WK)
      RETURN
      END SUBROUTINE DMUMPS_SCATTER_ROOT
      SUBROUTINE DMUMPS_GATHER_ROOT(MYID, M, N, ASEQ,
     &                    LOCAL_M, LOCAL_N,
     &                    MBLOCK, NBLOCK,
     &                    APAR,
     &                    MASTER_ROOT,
     &                    NPROW, NPCOL,
     &                    COMM)
      IMPLICIT NONE
      INTEGER MYID, MASTER_ROOT, COMM
      INTEGER M, N
      INTEGER NPROW, NPCOL
      INTEGER LOCAL_M, LOCAL_N
      INTEGER MBLOCK, NBLOCK
      DOUBLE PRECISION APAR( LOCAL_M, LOCAL_N )
      DOUBLE PRECISION ASEQ( M, N )
      INCLUDE 'mpif.h'
      INTEGER I, J, SIZE_IBLOCK, SIZE_JBLOCK, ISOUR, IROW, ICOL
      INTEGER IBLOCK, JBLOCK, II, JJ, KK
      INTEGER IAPAR, JAPAR, IERR, allocok
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: WK
      LOGICAL JUPDATE
      ALLOCATE(WK( MBLOCK * NBLOCK ), stat=allocok)
      if(allocok.ne.0) then
         WRITE(6,*) ' Allocation error of WK in '
     &       // 'routine DMUMPS_GATHER_ROOT '
         CALL MUMPS_ABORT()
      endif
        IAPAR = 1
        JAPAR = 1
        DO J = 1, N, NBLOCK
          SIZE_JBLOCK = NBLOCK
          IF ( J + NBLOCK > N ) THEN
            SIZE_JBLOCK = N - J + 1
          END IF
          JUPDATE = .FALSE.
          DO I = 1, M, MBLOCK
            SIZE_IBLOCK = MBLOCK
            IF ( I + MBLOCK > M ) THEN
              SIZE_IBLOCK = M - I + 1
            END IF
            IBLOCK = I / MBLOCK
            JBLOCK = J / NBLOCK
            IROW = mod ( IBLOCK, NPROW )
            ICOL = mod ( JBLOCK, NPCOL )
            ISOUR = IROW * NPCOL + ICOL
            IF ( ISOUR .NE. MASTER_ROOT ) THEN
              IF ( MYID .EQ. MASTER_ROOT ) THEN
                CALL MPI_RECV( WK(1), SIZE_IBLOCK*SIZE_JBLOCK,
     &                         MPI_DOUBLE_PRECISION,
     &                         ISOUR, 128, COMM, STATUS, IERR )
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  ASEQ(II,JJ)=WK(KK)
                  KK=KK+1
                END DO
                END DO
              ELSE IF ( MYID .EQ. ISOUR ) THEN
                KK=1
                DO JJ=JAPAR,JAPAR+SIZE_JBLOCK-1
                DO II=IAPAR,IAPAR+SIZE_IBLOCK-1
                  WK(KK)=APAR(II,JJ)
                  KK=KK+1
                END DO
                END DO
                CALL MPI_SSEND( WK( 1 ),
     &                         SIZE_IBLOCK*SIZE_JBLOCK,
     &                         MPI_DOUBLE_PRECISION,
     &                         MASTER_ROOT,128,COMM,IERR)
                JUPDATE = .TRUE.
                IAPAR = IAPAR + SIZE_IBLOCK
              END IF
            ELSE IF ( MYID.EQ. MASTER_ROOT ) THEN
              ASEQ(I:I+SIZE_IBLOCK-1,J:J+SIZE_JBLOCK-1)
     &        = APAR( IAPAR:IAPAR+SIZE_IBLOCK-1,
     &                JAPAR:JAPAR+SIZE_JBLOCK-1 )
              JUPDATE = .TRUE.
              IAPAR = IAPAR + SIZE_IBLOCK
            END IF
          END DO
          IF ( JUPDATE ) THEN
            IAPAR = 1
            JAPAR = JAPAR + SIZE_JBLOCK
          END IF
        END DO
        DEALLOCATE(WK)
      RETURN
      END SUBROUTINE DMUMPS_GATHER_ROOT
      SUBROUTINE DMUMPS_ROOT_ALLOC_STATIC(root, IROOT, N,
     &                  IW, LIW, A, LA,
     &                  FILS, DAD, MYID, SLAVEF, PROCNODE_STEPS,
     &                  LPTRAR, NELT, FRTPTR, FRTELT,
     &                  PTRAIW, PTRARW, INTARR, DBLARR,
     &                  LRLU, IPTRLU,
     &                  IWPOS, IWPOSCB, PTRIST, PTRAST,
     &                  STEP, PIMASTER, PAMASTER, ITLOC, RHS_MUMPS,
     &                  COMP, LRLUS, IFLAG, KEEP,KEEP8,DKEEP,IERROR )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INTEGER MYID
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION DKEEP(230)
      TYPE (DMUMPS_ROOT_STRUC ) :: root
      INTEGER(8) :: LA, LRLU, IPTRLU, LRLUS
      INTEGER IROOT, LIW, N, IWPOS, IWPOSCB
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, INTENT(IN) :: SLAVEF
      INTEGER, INTENT(IN) :: PROCNODE_STEPS(KEEP(28))
      INTEGER PTRIST(KEEP(28)), STEP(N)
      INTEGER(8) :: PTRAST(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PIMASTER(KEEP(28))
      INTEGER ITLOC( N + KEEP(253) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER COMP, IFLAG, IERROR
      INCLUDE 'mumps_headers.h'
      INTEGER FILS( N ), DAD(KEEP(28))
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1), FRTELT( NELT )
      INTEGER(8), INTENT(IN) :: PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER INTARR(KEEP8(27))
      DOUBLE PRECISION DBLARR(KEEP8(26))
      INTEGER numroc
      EXTERNAL numroc
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER(8) :: LREQA_ROOT
      INTEGER LREQI_ROOT, LOCAL_M, LOCAL_N, allocok
      LOGICAL :: EARLYT3ROOTINS
          LOCAL_M = numroc( root%ROOT_SIZE, root%MBLOCK,
     &              root%MYROW, 0, root%NPROW )
          LOCAL_M = max( 1, LOCAL_M )
          LOCAL_N = numroc( root%ROOT_SIZE, root%NBLOCK,
     &              root%MYCOL, 0, root%NPCOL )
          IF (KEEP(253).GT.0) THEN
            root%RHS_NLOC = numroc( KEEP(253), root%NBLOCK,
     &              root%MYCOL, 0, root%NPCOL )
            root%RHS_NLOC = max(1, root%RHS_NLOC)
          ELSE
            root%RHS_NLOC = 1
          ENDIF
          IF (associated( root%RHS_ROOT) ) 
     &             DEALLOCATE (root%RHS_ROOT)
          ALLOCATE(root%RHS_ROOT(LOCAL_M,root%RHS_NLOC), 
     &              stat=allocok)
          IF ( allocok.GT.0) THEN
            IFLAG=-13
            IERROR = LOCAL_M*root%RHS_NLOC
            RETURN
          ENDIF
          IF (KEEP(253).NE.0) THEN
            root%RHS_ROOT = ZERO
            CALL DMUMPS_ASM_RHS_ROOT (  N, FILS, 
     &           root, KEEP, RHS_MUMPS, 
     &           IFLAG, IERROR )
            IF ( IFLAG .LT. 0 ) RETURN
          ENDIF
          IF (KEEP(60) .NE. 0) THEN
            PTRIST(STEP(IROOT)) = -6666666
          ELSE
            LREQI_ROOT = 2 + KEEP(IXSZ)
            LREQA_ROOT = int(LOCAL_M,8) * int(LOCAL_N,8)
            IF (LREQA_ROOT.EQ.0_8) THEN
              PTRIST(STEP(IROOT)) = -9999999
              RETURN
            ENDIF
            CALL DMUMPS_ALLOC_CB(.FALSE.,0_8,.FALSE.,.FALSE.,
     &                     MYID,N,KEEP,KEEP8,DKEEP,IW,LIW,A,LA,
     &                     LRLU, IPTRLU,
     &                     IWPOS, IWPOSCB, SLAVEF, PROCNODE_STEPS, DAD,
     &                     PTRIST, PTRAST,
     &                     STEP, PIMASTER, PAMASTER, LREQI_ROOT,
     &                     LREQA_ROOT, IROOT, S_NOTFREE, .TRUE., COMP,
     &                     LRLUS, KEEP8(67), IFLAG, IERROR
     &           )
            IF ( IFLAG .LT. 0 ) RETURN
            PTRIST  ( STEP(IROOT) ) = IWPOSCB + 1
            PAMASTER( STEP(IROOT) ) = IPTRLU  + 1_8
            IW( IWPOSCB + 1 + KEEP(IXSZ)) = - LOCAL_N
            IW( IWPOSCB + 2 + KEEP(IXSZ)) =   LOCAL_M
          ENDIF
          EARLYT3ROOTINS = KEEP(200) .EQ.0
          IF (LOCAL_N > 0 .AND. .NOT. EARLYT3ROOTINS ) THEN
            IF (KEEP(60) .EQ. 0) THEN
              CALL DMUMPS_SET_TO_ZERO(A(IPTRLU+1_8), LOCAL_M,
     &        LOCAL_M, LOCAL_N, KEEP)
            ELSE
              CALL DMUMPS_SET_TO_ZERO(root%SCHUR_POINTER(1),
     &        root%SCHUR_LLD, LOCAL_M, LOCAL_N, KEEP)
            ENDIF
            IF (KEEP(55) .eq. 0) THEN
              IF (KEEP(60) .EQ. 0) THEN
                CALL DMUMPS_ASM_ARR_ROOT( N, root, IROOT,
     &          A(IPTRLU+1_8), LOCAL_M, LOCAL_M, LOCAL_N,
     &          FILS, PTRAIW, PTRARW, INTARR, DBLARR,
     &          KEEP8(27), KEEP8(26), MYID )
              ELSE
                CALL DMUMPS_ASM_ARR_ROOT( N, root, IROOT,
     &          root%SCHUR_POINTER(1), root%SCHUR_LLD, LOCAL_M, LOCAL_N,
     &          FILS, PTRAIW, PTRARW, INTARR, DBLARR,
     &          KEEP8(27), KEEP8(26), MYID )
              ENDIF
            ELSE
              IF (KEEP(60) .EQ. 0) THEN
                CALL DMUMPS_ASM_ELT_ROOT( N, root,
     &          A(IPTRLU+1_8), LOCAL_M, LOCAL_M, LOCAL_N,
     &          LPTRAR, NELT, FRTPTR, FRTELT,
     &          PTRAIW, PTRARW, INTARR, DBLARR,
     &          KEEP8(27), KEEP8(26), KEEP, KEEP8, MYID )
              ELSE
                CALL DMUMPS_ASM_ELT_ROOT( N, root,
     &           root%SCHUR_POINTER(1), root%SCHUR_LLD,
     &           root%SCHUR_MLOC, root%SCHUR_NLOC,
     &           LPTRAR, NELT, FRTPTR, FRTELT,
     &           PTRAIW, PTRARW, INTARR, DBLARR,
     &           KEEP8(27), KEEP8(26), KEEP, KEEP8, MYID )
              ENDIF
            ENDIF
          ENDIF
      RETURN
      END SUBROUTINE DMUMPS_ROOT_ALLOC_STATIC
      SUBROUTINE DMUMPS_ASM_ELT_ROOT( N, root, 
     &       VALROOT, LOCAL_M_LLD, LOCAL_M, LOCAL_N,
     &       LPTRAR, NELT, FRTPTR, FRTELT,
     &       PTRAIW, PTRARW,
     &       INTARR, DBLARR, LINTARR, LDBLARR,
     &       KEEP, KEEP8,
     &       MYID)
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER :: N, MYID, LOCAL_M, LOCAL_N, KEEP(500)
      INTEGER :: LOCAL_M_LLD
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION VALROOT(LOCAL_M_LLD,LOCAL_N)
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1), FRTELT( NELT )
      INTEGER(8), INTENT(IN) :: PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER(8), INTENT(IN) :: LINTARR, LDBLARR
      INTEGER, INTENT(INOUT) :: INTARR(LINTARR)
      DOUBLE PRECISION DBLARR(LDBLARR)
      INTEGER(8) :: J1, J2, K8, IPTR
      INTEGER    :: IELT, I, J, IGLOB, JGLOB, SIZEI, IBEG
      INTEGER    :: ARROW_ROOT
      INTEGER    :: IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID
      INTEGER    :: ILOCROOT, JLOCROOT
      ARROW_ROOT = 0
      DO IPTR = FRTPTR(KEEP(38)), FRTPTR(KEEP(38)+1) - 1
        IELT = FRTELT( IPTR )
        J1 = PTRAIW(IELT)
        J2 = PTRAIW(IELT+1)-1
        K8 = PTRARW(IELT)
        SIZEI=int(J2-J1)+1
        DO J=1, SIZEI
          JGLOB          = INTARR(J1+J-1)
          INTARR(J1+J-1) = root%RG2L_ROW(JGLOB)
        ENDDO
        DO J = 1, SIZEI
          JGLOB         = INTARR(J1+J-1)
          IF ( KEEP(50).eq. 0 ) THEN
            IBEG = 1
          ELSE
            IBEG = J
          END IF
          DO I = IBEG, SIZEI
            IGLOB = INTARR(J1+I-1)
            IF ( KEEP(50).eq.0 ) THEN
              IPOSROOT = INTARR(J1+I-1)
              JPOSROOT = INTARR(J1+J-1)
            ELSE
              IF ( INTARR(J1+I-1).GT. INTARR(J1+J-1) ) THEN
                IPOSROOT = INTARR(J1+I-1)
                JPOSROOT = INTARR(J1+J-1)
              ELSE
                IPOSROOT = INTARR(J1+J-1)
                JPOSROOT = INTARR(J1+I-1)
              END IF
            END IF
            IROW_GRID = mod( ( IPOSROOT - 1 )/root%MBLOCK,
     &                         root%NPROW )
            JCOL_GRID = mod( ( JPOSROOT - 1 )/root%NBLOCK,
     &                         root%NPCOL )
            IF ( IROW_GRID.EQ.root%MYROW .AND.
     &           JCOL_GRID.EQ.root%MYCOL ) THEN
              ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &               ( root%MBLOCK * root%NPROW ) )
     &             + mod( IPOSROOT - 1, root%MBLOCK ) + 1
              JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &               ( root%NBLOCK * root%NPCOL ) )
     &             + mod( JPOSROOT - 1, root%NBLOCK ) + 1
              VALROOT( ILOCROOT, JLOCROOT ) = 
     &        VALROOT( ILOCROOT, JLOCROOT ) + DBLARR(K8)
            ENDIF
            K8 = K8 + 1_8
          END DO
        END DO
        ARROW_ROOT = ARROW_ROOT + int(PTRARW(IELT+1_8)-PTRARW(IELT))
      END DO
      KEEP(49) = ARROW_ROOT
      RETURN
      END SUBROUTINE DMUMPS_ASM_ELT_ROOT
      SUBROUTINE DMUMPS_ASM_RHS_ROOT
     &           ( N, FILS, root, KEEP, RHS_MUMPS, 
     &             IFLAG, IERROR )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INTEGER N, KEEP(500), IFLAG, IERROR
      INTEGER FILS(N)
      TYPE (DMUMPS_ROOT_STRUC ) :: root
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER JCOL, IPOS_ROOT, JPOS_ROOT,
     &        IROW_GRID, JCOL_GRID, ILOCRHS, JLOCRHS,
     &        INODE
      INODE = KEEP(38)
      DO WHILE (INODE.GT.0)
        IPOS_ROOT = root%RG2L_ROW( INODE ) 
        IROW_GRID  = mod( ( IPOS_ROOT - 1 ) / root%MBLOCK, root%NPROW )
        IF (  IROW_GRID .NE. root%MYROW ) GOTO 100 
        ILOCRHS = root%MBLOCK * ( ( IPOS_ROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOS_ROOT - 1, root%MBLOCK ) + 1
        DO JCOL = 1, KEEP(253) 
          JPOS_ROOT = JCOL
          JCOL_GRID  = mod((JPOS_ROOT-1)/root%NBLOCK, root%NPCOL)
          IF (JCOL_GRID.NE.root%MYCOL ) CYCLE
           JLOCRHS = root%NBLOCK * ( ( JPOS_ROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOS_ROOT - 1, root%NBLOCK ) + 1
          root%RHS_ROOT(ILOCRHS, JLOCRHS) =
     &                 RHS_MUMPS(INODE+(JCOL-1)*KEEP(254))
        ENDDO
 100    CONTINUE
        INODE=FILS(INODE)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ASM_RHS_ROOT
      SUBROUTINE DMUMPS_ASM_ARR_ROOT( N, root, IROOT,
     &   VALROOT, LOCAL_M_LLD, LOCAL_M, LOCAL_N, FILS,
     &       PTRAIW, PTRARW,
     &       INTARR, DBLARR, LINTARR, LDBLARR,
     &       MYID)
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER :: N, MYID, IROOT, LOCAL_M, LOCAL_N
      INTEGER :: LOCAL_M_LLD
      INTEGER FILS( N )
      INTEGER(8), INTENT(IN) :: PTRARW( N ), PTRAIW( N )
      DOUBLE PRECISION VALROOT(LOCAL_M_LLD,LOCAL_N)
      INTEGER(8), INTENT(IN) :: LINTARR, LDBLARR
      INTEGER INTARR(LINTARR)
      DOUBLE PRECISION DBLARR(LDBLARR)
      DOUBLE PRECISION VAL
      INTEGER(8) :: JJ, J1,JK, J2,J3, J4, AINPUT
      INTEGER IORG, IBROT, NUMORG,
     &        IROW, JCOL
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID
      INTEGER ILOCROOT, JLOCROOT
      NUMORG = root%ROOT_SIZE
      IBROT  = IROOT
      DO IORG = 1, NUMORG
        JK = PTRAIW(IBROT)
        AINPUT = PTRARW(IBROT)
        IBROT = FILS(IBROT)
        JJ = JK + 1
        J1 = JJ + 1
        J2 = J1 + INTARR(JK)
        J3 = J2 + 1
        J4 = J2 - INTARR(JJ)
        JCOL = INTARR(J1)
        DO JJ = J1, J2
         IROW = INTARR(JJ)
         VAL  = DBLARR(AINPUT)
         AINPUT = AINPUT + 1
         IPOSROOT = root%RG2L_ROW( IROW )
         JPOSROOT = root%RG2L_COL( JCOL )
         IROW_GRID  = mod( ( IPOSROOT - 1 ) / root%MBLOCK, root%NPROW )
         JCOL_GRID  = mod( ( JPOSROOT - 1 ) / root%NBLOCK, root%NPCOL )
         IF ( IROW_GRID .EQ. root%MYROW .AND.
     &        JCOL_GRID .EQ. root%MYCOL ) THEN
            ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
            JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
            VALROOT( ILOCROOT, JLOCROOT ) =
     &      VALROOT( ILOCROOT, JLOCROOT ) + VAL
         END IF
        END DO
        IF (J3 .LE. J4) THEN
         IROW =  INTARR(J1)
         DO JJ= J3,J4
          JCOL = INTARR(JJ)
          VAL  = DBLARR(AINPUT)
          AINPUT = AINPUT + 1
          IPOSROOT = root%RG2L_ROW( IROW )
          JPOSROOT = root%RG2L_COL( JCOL )
          IROW_GRID= mod( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW)
          JCOL_GRID= mod( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL)
          IF ( IROW_GRID .EQ. root%MYROW .AND.
     &        JCOL_GRID .EQ. root%MYCOL ) THEN
            ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
            JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
            VALROOT( ILOCROOT, JLOCROOT ) = 
     &      VALROOT( ILOCROOT, JLOCROOT ) + VAL
          END IF
         END DO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ASM_ARR_ROOT
