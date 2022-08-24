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
      RECURSIVE SUBROUTINE DMUMPS_MAPLIG( COMM_LOAD, ASS_IRECV,
     &  BUFR, LBUFR, LBUFR_BYTES,
     &
     &  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     &  NFRONT_PERE, NASS_PERE, NFS4FATHER, LMAP, TROW,
     &  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &  LRLUS, N, IW,
     &  LIW, A, LA,
     &  PTRIST, PTLUST, PTRFAC,
     &  PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &  IFLAG, IERROR, MYID, COMM, PERM, IPOOL, LPOOL, LEAF,
     &  NBFIN, ICNTL, KEEP,KEEP8,DKEEP,
     &  root, OPASSW, OPELIW,
     &  ITLOC, RHS_MUMPS,
     &  FILS, DAD, PTRARW, PTRAIW, INTARR, DBLARR, ND, FRERE,
     &  LPTRAR, NELT, FRTPTR, FRTELT, 
     &
     &  ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &               , LRGROUPS
     &  )
      USE DMUMPS_BUF
      USE DMUMPS_LOAD
      USE DMUMPS_LR_DATA_M
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR
      USE DMUMPS_FAC_FRONT_AUX_M, 
     &                ONLY : DMUMPS_COMPUTE_SIZE_SCHUR_IN_FRONT
#if ! defined(NO_FDM_MAPROW)
      USE MUMPS_FAC_MAPROW_DATA_M
#endif
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
#if ! defined(NO_FDM_MAPROW)
#endif
      TYPE (DMUMPS_ROOT_STRUC ) :: root
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER ICNTL( 60 ), KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(230)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER SLAVEF, NBFIN
      INTEGER(8) :: LA, IPTRLU, LRLU, LRLUS, POSFAC
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(N)
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER COMP
      INTEGER NSTK( KEEP(28) )
      INTEGER PERM(N)
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER INODE_PERE, ISON
      INTEGER :: NFS4FATHER
      INTEGER NBROWS_ALREADY_SENT
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE( * )
      INTEGER LMAP 
      INTEGER TROW( LMAP )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION DBLARR(KEEP8(26))
      INTEGER INTARR(KEEP8(27))
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N ), DAD( KEEP(28) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8), INTENT(IN) ::  PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER NOSLA, I
      INTEGER I_POSMYIDIN_PERE
      INTEGER INDICE_PERE
      INTEGER PDEST, PDEST_MASTER
      LOGICAL :: LOCAL_ASSEMBLY_TO_BE_DONE
      INTEGER NROWS_TO_SEND
      INTEGER PDEST_MASTER_ISON, IPOS_IN_SLAVE
      LOGICAL DESCLU, SLAVE_ISON
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED 
      INTEGER MSGSOU, MSGTAG
      INTEGER LP
      LOGICAL PACKED_CB
      LOGICAL IS_ERROR_BROADCASTED, IS_ofType5or6
      INTEGER ITYPE_SON, TYPESPLIT
      INTEGER :: KEEP253_LOC 
      INTEGER :: NVSCHUR, NSLAVES_L, NROW_L, IROW_L, NASS_L, NELIM_L
      LOGICAL :: CB_IS_LR
      INTEGER :: IWXXF_HANDLER
      DOUBLE PRECISION :: ADummy(1)
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: SON_A
      INTEGER(8) :: IACHK, RECSIZE 
#if ! defined(NO_FDM_MAPROW)
      INTEGER :: INFO_TMP(2)
#endif
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP, PERM_LOC
      IS_ERROR_BROADCASTED = .FALSE.
      TYPESPLIT = MUMPS_TYPESPLIT( PROCNODE_STEPS(STEP(INODE_PERE)),
     &                             KEEP(199) )
      ITYPE_SON = MUMPS_TYPENODE( PROCNODE_STEPS(STEP(ISON)),
     &                            KEEP(199) )
      IS_ofType5or6 = ((TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6))
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      CB_IS_LR = (IW(PTRIST(STEP(ISON))+XXLR).EQ.1 .OR.
     &                   IW(PTRIST(STEP(ISON))+XXLR).EQ.3)
      IWXXF_HANDLER = IW(PTRIST(STEP(ISON))+XXF)
#if ! defined(NO_FDM_MAPROW)
#endif
      ALLOCATE(SLAVES_PERE(0:max(1,NSLAVES_PERE)), stat=allocok)
      if (allocok .GT. 0) THEN
         IF (LP > 0) THEN
            write(LP,*) MYID,
     &           ' : PB allocation SLAVES_PERE in DMUMPS_MAPLIG'
         ENDIF
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      endif
      IF (NSLAVES_PERE.GT.0) 
     &SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE_PERE)),
     &                 KEEP(199) )
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      if (allocok .GT. 0) THEN
         IF (LP>0) THEN
            write(LP,*) MYID,
     &           ' : PB allocation NBROW in DMUMPS_MAPLIG'
         ENDIF
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 670
      endif
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP>0) THEN
        write(LP,*) MYID, ' : PB allocation LMAP in DMUMPS_MAPLIG'
        ENDIF
        IFLAG  =-13
        IERROR = LMAP
        GOTO 680
      endif
      MAP( 1 : LMAP ) = TROW( 1 : LMAP )
      PDEST_MASTER_ISON = MUMPS_PROCNODE(PROCNODE_STEPS(STEP(ISON)),
     &                    KEEP(199))
      SLAVE_ISON = PDEST_MASTER_ISON .NE. MYID
      IF (SLAVE_ISON) THEN
        IF ( PTRIST(STEP( ISON )) .EQ. 0 ) THEN
          CALL DMUMPS_TREAT_DESCBAND( ISON, COMM_LOAD,
     &    ASS_IRECV, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &    IFLAG, IERROR, COMM,
     &    PERM,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE, LPTRAR,
     &    NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &    )
          IF ( IFLAG .LT. 0 ) THEN
            IS_ERROR_BROADCASTED = .TRUE.  
            GOTO 670                       
          ENDIF
        END IF
#if ! defined(NO_FDM_MAPROW)
        IF ( 
     &     ( IW( PTRIST(STEP(ISON)) + 1 + KEEP(IXSZ) ) .NE.
     &       IW( PTRIST(STEP(ISON)) + 3 + KEEP(IXSZ) ) ) .OR.
     &     ( KEEP(50) .NE. 0 .AND.
     &       IW( PTRIST(STEP(ISON)) + 6 + KEEP(IXSZ) ) .NE. 0 ) )
     &  THEN
          INFO_TMP=0 
          CALL MUMPS_FMRD_SAVE_MAPROW(
     &         IW(PTRIST(STEP(ISON))+XXA),
     &         INODE_PERE, ISON, NSLAVES_PERE, NFRONT_PERE,
     &         NASS_PERE, LMAP, NFS4FATHER,
     &         SLAVES_PERE(1:NSLAVES_PERE),
     &         MAP,        
     &         INFO_TMP)
               IF (INFO_TMP(1) < 0) THEN
                 IFLAG = INFO_TMP(1)
                 IERROR = INFO_TMP(2)
               ENDIF
          GOTO 670 
        ELSE
          GOTO 10
        ENDIF
#endif
        DO WHILE (
     &     ( IW( PTRIST(STEP(ISON)) + 1 + KEEP(IXSZ) ) .NE.
     &       IW( PTRIST(STEP(ISON)) + 3 + KEEP(IXSZ) ) ) .OR.
     &     ( KEEP(50) .NE. 0 .AND.
     &       IW( PTRIST(STEP(ISON)) + 6 + KEEP(IXSZ) ) .NE. 0 ) )
          IF ( KEEP(50).eq.0) THEN
            MSGSOU = PDEST_MASTER_ISON
            MSGTAG = BLOC_FACTO
          ELSE
            IF ( IW( PTRIST(STEP(ISON)) + 1 + KEEP(IXSZ) ) .NE.
     &           IW( PTRIST(STEP(ISON)) + 3 + KEEP(IXSZ) ) ) THEN
              MSGSOU = PDEST_MASTER_ISON
              MSGTAG = BLOC_FACTO_SYM
            ELSE
              MSGSOU = MPI_ANY_SOURCE
              MSGTAG = BLOC_FACTO_SYM_SLAVE
            END IF
          END IF
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MSGSOU, MSGTAG,
     &    STATUS, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &    IFLAG, IERROR, COMM,
     &    PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE, LPTRAR,
     &    NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &    )
          IF ( IFLAG .LT. 0 ) THEN
            IS_ERROR_BROADCASTED = .TRUE.  
            GOTO 670                       
          ENDIF
        END DO
      ENDIF
#if ! defined(NO_FDM_MAPROW)
 10   CONTINUE
#endif
      IF ( NSLAVES_PERE .EQ. 0 ) THEN
        NBROW( 0 ) = LMAP_LOC
      ELSE
        DO I = 0, NSLAVES_PERE
          NBROW( I ) = 0
        END DO
        DO I = 1, LMAP_LOC
          INDICE_PERE = MAP( I )
          CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           NASS_PERE,
     &           NFRONT_PERE - NASS_PERE,
     &           NSLAVES_PERE,
     &           INDICE_PERE,
     &           NOSLA,
     &           IPOS_IN_SLAVE )
          NBROW( NOSLA ) = NBROW( NOSLA ) + 1
        END DO
        DO I = 1, NSLAVES_PERE
          NBROW(I)=NBROW(I)+NBROW(I-1)
        ENDDO
      ENDIF
      ALLOCATE(PERM_LOC(LMAP_LOC), stat=allocok)
      IF (allocok .GT. 0) THEN
          IF (LP.GT.0) THEN
          write(LP,*) MYID,': PB allocation PERM_LOC in DMUMPS_MAPLIG'
          ENDIF
          IFLAG  =-13
          IERROR = LMAP_LOC
          GOTO 670
      ENDIF
      KEEP253_LOC   = 0
      DO I = LMAP_LOC, 1, -1
          INDICE_PERE = MAP( I )
          IF (INDICE_PERE > NFRONT_PERE - KEEP(253)) THEN
             KEEP253_LOC = KEEP253_LOC + 1
          ENDIF
          CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           NASS_PERE,
     &           NFRONT_PERE - NASS_PERE,
     &           NSLAVES_PERE,
     &           INDICE_PERE,
     &           NOSLA,
     &           IPOS_IN_SLAVE )
          PERM_LOC( NBROW( NOSLA ) ) = I
          NBROW( NOSLA ) = NBROW( NOSLA ) - 1
      ENDDO
      DO I = 0, NSLAVES_PERE
          NBROW(I)=NBROW(I)+1
      END DO
      IF ((KEEP(114).EQ.1) .AND. (KEEP(50).EQ.2) .AND.
     &    (KEEP(116).GT.0) .AND. ((LMAP_LOC-KEEP253_LOC).GT.0) 
     &   ) THEN
        IF (ITYPE_SON.EQ.1) THEN
          NELIM_L   =  IW(PTLUST(STEP(ISON))+1+KEEP(IXSZ)) 
          NASS_L    =  NELIM_L +
     &                 IW(PTLUST(STEP(ISON))+3+KEEP(IXSZ))    
          IROW_L    =  PTLUST(STEP(ISON))+6+KEEP(IXSZ)+NASS_L
          NROW_L    = LMAP_LOC
        ELSE
         NROW_L    = LMAP_LOC 
         NSLAVES_L = IW( PTRIST(STEP( ISON )) + 5 + KEEP(IXSZ) )
         IROW_L    = PTRIST(STEP(ISON)) + 6 + NSLAVES_L + KEEP(IXSZ)
        ENDIF
        CALL DMUMPS_COMPUTE_SIZE_SCHUR_IN_FRONT ( 
     &     N, 
     &     NROW_L-KEEP253_LOC, 
     &     KEEP(116), 
     &     IW(IROW_L),
     &     PERM, NVSCHUR )
      ELSE
         NVSCHUR = 0
      ENDIF
      PDEST_MASTER = SLAVES_PERE(0)
      I_POSMYIDIN_PERE = -99999
      LOCAL_ASSEMBLY_TO_BE_DONE = .FALSE.
      DO I = 0, NSLAVES_PERE
        IF (SLAVES_PERE(I) .EQ. MYID) THEN
          I_POSMYIDIN_PERE = I
          LOCAL_ASSEMBLY_TO_BE_DONE = .TRUE.
#if ! defined(NO_FDM_DESCBAND)
          IF (PTRIST(STEP(INODE_PERE)) .EQ. 0
     &      .AND. MYID .NE. PDEST_MASTER) THEN
            CALL DMUMPS_TREAT_DESCBAND( INODE_PERE, COMM_LOAD,
     &      ASS_IRECV, 
     &      BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &      IFLAG, IERROR, COMM,
     &      PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &    
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, DAD, PTRARW, PTRAIW,
     &      INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE, LPTRAR,
     &      NELT, FRTPTR, FRTELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &      )
            IF ( IFLAG .LT. 0 ) THEN
              IS_ERROR_BROADCASTED = .TRUE. 
              GOTO 600                      
            ENDIF
          ENDIF
#endif
        ENDIF
      END DO
      IF (KEEP(120).NE.0 .AND. LOCAL_ASSEMBLY_TO_BE_DONE) THEN
        CALL DMUMPS_LOCAL_ASSEMBLY_TYPE2(I_POSMYIDIN_PERE,
     &     SLAVES_PERE(I_POSMYIDIN_PERE),  
     &     MYID, PDEST_MASTER, ISON, INODE_PERE, 
     &     NSLAVES_PERE, NASS_PERE, NFRONT_PERE, NFS4FATHER,
     &     LMAP_LOC, MAP, NBROW, PERM_LOC,
     &     IS_ofType5or6, IFLAG, IERROR,
     &     N, SLAVEF, KEEP, IPOOL, LPOOL, STEP,
     &     PROCNODE_STEPS, COMM_LOAD, ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &     KEEP8, IW, LIW, A, LA, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &     PTRIST, PTLUST, PTRAST, PAMASTER, PIMASTER, ND,
     &     NELT, FRTPTR, FRTELT,
     &     OPASSW, OPELIW,
     &     ITLOC, RHS_MUMPS, KEEP253_LOC, NVSCHUR,
     &     FILS, DAD, LPTRAR, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &     ITYPE_SON, LRGROUPS)
        LOCAL_ASSEMBLY_TO_BE_DONE = .FALSE.
        IF (IFLAG < 0) THEN
          GOTO 600
        ENDIF
      ENDIF
      DO I = NSLAVES_PERE, 0, -1   
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
           DESCLU = .FALSE.
           NBROWS_ALREADY_SENT = 0
           IF (I == NSLAVES_PERE) THEN
             NROWS_TO_SEND=LMAP_LOC-NBROW(I)+1
           ELSE
             NROWS_TO_SEND=NBROW(I+1)-NBROW(I)
           ENDIF
           PACKED_CB=(IW(PTRIST(STEP(ISON))+XXS).EQ.S_CB1COMP)
           IERR = -1
           DO WHILE (IERR .EQ. -1)
             IF ( IW ( PTRIST(STEP(ISON) )+KEEP(IXSZ) )
     &             .GT. N + KEEP(253) ) THEN
               WRITE(*,*) MYID,': Internal error in Maplig'
               WRITE(*,*) MYID,': PTRIST(STEP(ISON))/N=',
     &                            PTRIST(STEP(ISON)), N
               WRITE(*,*) MYID,': I, NBROW(I)=',I, NBROW(I)
               WRITE(*,*) MYID,': NSLAVES_PERE=',NSLAVES_PERE
               WRITE(*,*) MYID,': ISON, INODE_PERE=',ISON,INODE_PERE
               WRITE(*,*) MYID,': Son header=',
     &         IW(PTRIST(STEP(ISON)): PTRIST(STEP(ISON))+5+KEEP(IXSZ))
               CALL MUMPS_ABORT()
             END IF
             IF (NROWS_TO_SEND .EQ. 0 .AND. PDEST.NE.PDEST_MASTER) THEN
                IERR = 0
                CYCLE
             ENDIF
             IF (CB_IS_LR) THEN
              CALL DMUMPS_BUF_SEND_CONTRIB_TYPE2( 
     &                 NBROWS_ALREADY_SENT,
     &                 DESCLU, INODE_PERE,
     &                 NFRONT_PERE, NASS_PERE, NFS4FATHER,
     &                 NSLAVES_PERE, ISON,
     &                 NROWS_TO_SEND, LMAP_LOC, MAP,
     &                 PERM_LOC(min(LMAP_LOC,NBROW(I))),
     &                 IW( PTRIST(STEP(ISON))),
     &                 ADummy, 1_8,
     &                 I, PDEST, PDEST_MASTER, 
     &                 COMM, IERR, 
     &                 KEEP,KEEP8, STEP, N, SLAVEF,
     &                 ISTEP_TO_INIV2, TAB_POS_IN_PERE, PACKED_CB,
     &                 KEEP253_LOC, NVSCHUR, 
     &                 ITYPE_SON, MYID,
     &                 NPIV_CHECK = IW(PTLUST(STEP(ISON))+3+KEEP(IXSZ)))
             ELSE
               CALL DMUMPS_DM_SET_DYNPTR(
     &                   IW(PTRIST(STEP(ISON))+XXS), 
     &                   A, LA,
     &                   PTRAST(STEP(ISON)),  
     &                   IW(PTRIST(STEP(ISON))+XXD),
     &                   IW(PTRIST(STEP(ISON))+XXR), 
     &                   SON_A, IACHK, RECSIZE ) 
               CALL DMUMPS_BUF_SEND_CONTRIB_TYPE2( NBROWS_ALREADY_SENT,
     &                 DESCLU, INODE_PERE,
     &                 NFRONT_PERE, NASS_PERE, NFS4FATHER,
     &                 NSLAVES_PERE, ISON,
     &                 NROWS_TO_SEND, LMAP_LOC, MAP,
     &                 PERM_LOC(min(LMAP_LOC,NBROW(I))),
     &                 IW( PTRIST(STEP(ISON))),
     &                 SON_A(IACHK:IACHK+RECSIZE-1_8),
     &                 RECSIZE,
     &                 I, PDEST, PDEST_MASTER, 
     &                 COMM, IERR, 
     &                KEEP,KEEP8, STEP, N, SLAVEF,
     &                 ISTEP_TO_INIV2, TAB_POS_IN_PERE, PACKED_CB,
     &                 KEEP253_LOC, NVSCHUR, 
     &                 ITYPE_SON, MYID)
             ENDIF
             IF ( IERR .EQ. -2 ) THEN
               IFLAG  = -17
               IF (LP .GT. 0) THEN
                 WRITE(LP,*)
     &           "FAILURE: SEND BUFFER TOO SMALL IN DMUMPS_MAPLIG"
               ENDIF
               IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &         NROWS_TO_SEND * IW(PTRIST(STEP(ISON))+KEEP(IXSZ))
     &        * KEEP( 35 )
               GO TO 600
             END IF
             IF ( IERR .EQ. -3 ) THEN
               IF (LP .GT. 0) THEN
                 WRITE(LP,*)
     &           "FAILURE: RECV BUFFER TOO SMALL IN DMUMPS_MAPLIG"
               ENDIF
               IFLAG  = -20
               IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &         NROWS_TO_SEND * IW(PTRIST(STEP(ISON))+KEEP(IXSZ))
     &         * KEEP( 35 )
               GOTO 600
             ENDIF
             IF (KEEP(219).NE.0) THEN
              IF ( IERR .EQ. -4 ) THEN
                IFLAG  = -13
               IERROR = NFS4FATHER
               IF (LP .GT. 0) THEN
                 WRITE(LP, *)
     & "FAILURE: MAX_ARRAY allocation failed IN DMUMPS_MAPLIG"
               ENDIF
               GO TO 600
              END IF
             END IF
             IF ( IERR .EQ. -1 ) THEN
               IF (LOCAL_ASSEMBLY_TO_BE_DONE) THEN
                 CALL DMUMPS_LOCAL_ASSEMBLY_TYPE2(I_POSMYIDIN_PERE,
     &           SLAVES_PERE(I_POSMYIDIN_PERE),  
     &           MYID, PDEST_MASTER, ISON, INODE_PERE, 
     &           NSLAVES_PERE, NASS_PERE, NFRONT_PERE, NFS4FATHER,
     &           LMAP_LOC, MAP, NBROW, PERM_LOC,
     &           IS_ofType5or6, IFLAG, IERROR,
     &           N, SLAVEF, KEEP, IPOOL, LPOOL, STEP,
     &           PROCNODE_STEPS, COMM_LOAD, ISTEP_TO_INIV2,
     &           TAB_POS_IN_PERE,
     &           KEEP8, IW, LIW, A, LA, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &           PTRIST, PTLUST, PTRAST, PAMASTER, PIMASTER, ND,
     &           NELT, FRTPTR, FRTELT,
     &           OPASSW, OPELIW,
     &           ITLOC, RHS_MUMPS, KEEP253_LOC, NVSCHUR,
     &           FILS, DAD,
     &           LPTRAR, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &           ITYPE_SON, LRGROUPS)
                 LOCAL_ASSEMBLY_TO_BE_DONE = .FALSE.
                 IF (IFLAG < 0) THEN
                   GOTO 600
                 ENDIF
               ELSE
                 BLOCKING = .FALSE.
                 SET_IRECV = .TRUE.
                 MESSAGE_RECEIVED = .FALSE.
                 CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &           ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &           MPI_ANY_SOURCE, MPI_ANY_TAG,
     &           STATUS,
     &           BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &           IWPOS, IWPOSCB, IPTRLU,
     &           LRLU, LRLUS, N, IW, LIW, A, LA,
     &           PTRIST, PTLUST, PTRFAC,
     &           PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &           IFLAG, IERROR, COMM,
     &           PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &           root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD,
     &           PTRARW, PTRAIW,
     &           INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,LPTRAR,
     &           NELT, FRTPTR, FRTELT, 
     &           ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &           )
                 IF ( IFLAG .LT. 0 ) THEN
                   IS_ERROR_BROADCASTED=.TRUE.
                   GOTO 600
                 ENDIF
               END IF
             END IF 
           ENDDO
        ENDIF 
      END DO
      IF (LOCAL_ASSEMBLY_TO_BE_DONE) THEN
        CALL DMUMPS_LOCAL_ASSEMBLY_TYPE2(I_POSMYIDIN_PERE,
     &     SLAVES_PERE(I_POSMYIDIN_PERE),  
     &     MYID, PDEST_MASTER, ISON, INODE_PERE, 
     &     NSLAVES_PERE, NASS_PERE, NFRONT_PERE, NFS4FATHER,
     &     LMAP_LOC, MAP, NBROW, PERM_LOC,
     &     IS_ofType5or6, IFLAG, IERROR,
     &     N, SLAVEF, KEEP, IPOOL, LPOOL, STEP,
     &     PROCNODE_STEPS, COMM_LOAD, ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &     KEEP8, IW, LIW, A, LA, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &     PTRIST, PTLUST, PTRAST, PAMASTER, PIMASTER, ND,
     &     NELT, FRTPTR, FRTELT,
     &     OPASSW, OPELIW,
     &     ITLOC, RHS_MUMPS, KEEP253_LOC, NVSCHUR,
     &     FILS, DAD, LPTRAR, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &     ITYPE_SON, LRGROUPS)
        LOCAL_ASSEMBLY_TO_BE_DONE = .FALSE.
        IF (IFLAG < 0) THEN
          GOTO 600
        ENDIF
      ENDIF
      IF (CB_IS_LR) THEN
        CALL DMUMPS_BLR_FREE_CB_LRB(IWXXF_HANDLER,
     &             .FALSE., 
     &             KEEP8)
        IF ((KEEP(486).EQ.3).OR.KEEP(486).EQ.0) THEN
         CALL DMUMPS_BLR_END_FRONT(IWXXF_HANDLER, IFLAG, KEEP8)
        ENDIF
      ENDIF
      IF (KEEP(214) .EQ. 2) THEN
        CALL DMUMPS_STACK_BAND( N, ISON,
     &    PTRIST, PTRAST, PTLUST, PTRFAC, IW, LIW, A, LA,
     &    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     &    IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER,
     &    IFLAG, IERROR, SLAVEF, PROCNODE_STEPS, DAD, MYID,
     &    COMM, KEEP,KEEP8, DKEEP, ITYPE_SON )
        IF (IFLAG .LT. 0) THEN
          IS_ERROR_BROADCASTED = .TRUE.
          GOTO 600
        ENDIF
      ENDIF
      CALL DMUMPS_FREE_BAND( N, ISON, PTRIST, PTRAST, IW, LIW,
     &             A, LA, LRLU, LRLUS, IWPOSCB, IPTRLU,
     &             STEP, MYID, KEEP, KEEP8, ITYPE_SON
     &)
 600  CONTINUE
      DEALLOCATE(PERM_LOC)
 670  CONTINUE
      DEALLOCATE(MAP)
 680  CONTINUE
      DEALLOCATE(NBROW)
      DEALLOCATE(SLAVES_PERE)
 700  CONTINUE
      IF (IFLAG .LT. 0 .AND. .NOT. IS_ERROR_BROADCASTED) THEN
        CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_MAPLIG
      SUBROUTINE DMUMPS_MAPLIG_FILS_NIV1( COMM_LOAD, ASS_IRECV, 
     &  BUFR, LBUFR, LBUFR_BYTES,
     &
     &  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     &  NFRONT_PERE, NASS_PERE, NFS4FATHER, LMAP, TROW,
     &  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &  LRLUS, N, IW,
     &  LIW, A, LA,
     &  PTRIST, PTLUST, PTRFAC,
     &  PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &  IFLAG, IERROR, MYID, COMM, PERM, IPOOL, LPOOL, LEAF,
     &  NBFIN, ICNTL, KEEP,KEEP8,DKEEP, root,
     &  OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &  FILS, DAD, PTRARW, PTRAIW, INTARR, DBLARR,
     &  ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT, 
     &
     &  ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &  )
      USE DMUMPS_BUF
      USE DMUMPS_LOAD
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_FAC_LR, ONLY: DMUMPS_DECOMPRESS_PANEL
      USE DMUMPS_FAC_FRONT_AUX_M, 
     &                ONLY : DMUMPS_COMPUTE_SIZE_SCHUR_IN_FRONT
      USE DMUMPS_LR_DATA_M
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR
     &                                  , DMUMPS_DM_FREE_BLOCK
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER ICNTL( 60 ), KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(230)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER SLAVEF, NBFIN
      INTEGER(8) :: LA, IPTRLU, LRLU, LRLUS, POSFAC
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(N)
      INTEGER COMP
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER INODE_PERE, ISON
      INTEGER NFS4FATHER
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: M_ARRAY
      LOGICAL :: M_ARRAY_RETRIEVED
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE(NSLAVES_PERE)
      INTEGER NELIM, LMAP, TROW( LMAP ), NASS
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION DBLARR(KEEP8(26))
      INTEGER INTARR(KEEP8(27))
      INTEGER LPTRAR, NELT
      INTEGER IW( LIW )
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL )
      INTEGER NSTK( KEEP(28) ), ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER PERM(N)
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N ), DAD( KEEP(28) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8), INTENT(IN) :: PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER LP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER NOSLA, I, ISTCHK, ISTCHK_LOC
      INTEGER NBROWS_ALREADY_SENT 
      INTEGER INDICE_PERE
      INTEGER INDICE_PERE_ARRAY_ARG(1)
      INTEGER PDEST, PDEST_MASTER, NFRONT
      LOGICAL SAME_PROC, DESCLU
      INTEGER(8) :: IACHK, POSROW, ASIZE, RECSIZE
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: SON_A
      INTEGER(8) :: DYNSIZE
      INTEGER NSLSON, NBCOLS, NROW, NROWS_TO_SEND,
     &        NPIV, NROWS_TO_STACK, II, IROW_SON,
     &        IPOS_IN_SLAVE, DECR, ITYPE_SON
      INTEGER NBCOLS_EFF
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      LOGICAL PACKED_CB
      LOGICAL :: CB_IS_LR
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR
      INTEGER :: NB_BLR_COLS, NB_BLR_ROWS,
     &           NB_BLR_SHIFT, PANEL2DECOMPRESS, 
     &           CURRENT_PANEL_SIZE, PANEL_BEG_OFFSET,
     &           NROWS_ALREADY_STACKED, NROWS_TO_STACK_LOC
      INTEGER :: NVSCHUR, IROW_L
      INTEGER(8) :: LA_TEMP
      DOUBLE PRECISION :: ADummy(1)
      DOUBLE PRECISION, ALLOCATABLE :: A_TEMP(:)
      TYPE (LRB_TYPE), POINTER :: CB_LRB(:,:)
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP, PERM_LOC
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      if (NSLAVES_PERE.le.0) then
       write(6,*) ' error 2 in maplig_fils_niv1 ', NSLAVES_PERE
       CALL MUMPS_ABORT()
      endif
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      IF (allocok .GT. 0) THEN
        IF (LP > 0)
     &  write(LP,*) MYID,
     &  ' : PB allocation NBROW in DMUMPS_MAPLIG_FILS_NIV1'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      ENDIF
      ALLOCATE(SLAVES_PERE(0:NSLAVES_PERE), stat =allocok)
      IF ( allocok .GT. 0 ) THEN
        IF (LP > 0) write(LP,*) MYID,
     &  ' : PB allocation SLAVES_PERE in DMUMPS_MAPLIG_FILS_NIV1'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      ENDIF
      SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = MUMPS_PROCNODE( 
     &                       PROCNODE_STEPS(STEP(INODE_PERE)),
     &                       KEEP(199) )
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) MYID,
     &   ' : PB allocation LMAP in DMUMPS_MAPLIG_FILS_NIV1'
        IFLAG  =-13
        IERROR = LMAP_LOC
        GOTO 700
      endif
      MAP( 1 : LMAP_LOC ) = TROW( 1 : LMAP_LOC )
      DO I = 0, NSLAVES_PERE
        NBROW( I ) = 0
      END DO
      IF (NSLAVES_PERE == 0) THEN
        NBROW(0) = LMAP_LOC
      ELSE
       DO I = 1, LMAP_LOC
        INDICE_PERE = MAP( I )
        CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &         NASS_PERE,
     &         NFRONT_PERE - NASS_PERE,
     &         NSLAVES_PERE,
     &         INDICE_PERE,
     &         NOSLA,
     &         IPOS_IN_SLAVE )
        NBROW( NOSLA ) = NBROW( NOSLA ) + 1
       END DO
        DO I = 1, NSLAVES_PERE
          NBROW(I)=NBROW(I)+NBROW(I-1)
        ENDDO
      ENDIF
      ALLOCATE(PERM_LOC(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
         IF (LP > 0) THEN
            write(LP,*) MYID,
     &     ': PB allocation PERM_LOC in DMUMPS_MAPLIG_FILS_NIV1'
         ENDIF
        IFLAG  =-13
        IERROR = LMAP_LOC
        GOTO 700
      endif
        ISTCHK     = PIMASTER(STEP(ISON))
        NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
      DO I = LMAP_LOC, 1, -1
          INDICE_PERE = MAP( I )
          CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           NASS_PERE,
     &           NFRONT_PERE - NASS_PERE,
     &           NSLAVES_PERE,
     &           INDICE_PERE,
     &           NOSLA,
     &           IPOS_IN_SLAVE )
          PERM_LOC( NBROW( NOSLA ) ) = I
          NBROW( NOSLA ) = NBROW( NOSLA ) - 1
      ENDDO
      DO I = 0, NSLAVES_PERE
          NBROW(I)=NBROW(I)+1
      END DO
      PDEST_MASTER = MYID
      IF ( SLAVES_PERE(0) .NE. MYID ) THEN
        WRITE(*,*) 'Error 1 in MAPLIG_FILS_NIV1:',MYID, SLAVES_PERE
        CALL MUMPS_ABORT()
      END IF
      PDEST        = PDEST_MASTER
        I = 0
        ISTCHK     = PIMASTER(STEP(ISON))
        NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
        NELIM      = IW(ISTCHK+1+KEEP(IXSZ))
        NROW       = IW(ISTCHK+2+KEEP(IXSZ))
        NPIV       = IW(ISTCHK+3+KEEP(IXSZ))
        NASS       = NPIV+NELIM
        IF (NPIV.LT.0) THEN
         write(6,*) ' Error 2 in DMUMPS_MAPLIG_FILS_NIV1 ', NPIV
         CALL MUMPS_ABORT()
        ENDIF
        NSLSON     = IW(ISTCHK+5+KEEP(IXSZ))
        NFRONT     = NPIV + NBCOLS
        PACKED_CB=(IW(PTRIST(STEP(ISON))+XXS) .eq. S_CB1COMP)
        IF (I == NSLAVES_PERE) THEN
          NROWS_TO_STACK=LMAP_LOC-NBROW(I)+1
        ELSE
          NROWS_TO_STACK=NBROW(I+1)-NBROW(I)
        ENDIF
      IF ((KEEP(114).EQ.1) .AND.  (KEEP(50).EQ.2) .AND.
     &    (KEEP(116).GT.0) .AND. ((NFRONT-NASS-KEEP(253)).GT.0) 
     &   ) THEN
         IROW_L    = PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ) + NASS
         CALL DMUMPS_COMPUTE_SIZE_SCHUR_IN_FRONT ( 
     &     N, 
     &     NFRONT-NASS-KEEP(253), 
     &     KEEP(116), 
     &     IW(IROW_L),
     &     PERM, NVSCHUR )
      ELSE
         NVSCHUR = 0
      ENDIF
        DECR=1
        IW(PTLUST(STEP(INODE_PERE))+XXNBPR) =
     &  IW(PTLUST(STEP(INODE_PERE))+XXNBPR) - DECR
        IW(PTRIST(STEP(ISON))+XXNBPR) =
     &  IW(PTRIST(STEP(ISON))+XXNBPR) - DECR
        CB_IS_LR = (IW(ISTCHK+XXLR).EQ.1 .OR.
     &              IW(ISTCHK+XXLR).EQ.3)
        NROWS_ALREADY_STACKED = 0
 100    CONTINUE   
        NROWS_TO_STACK_LOC = NROWS_TO_STACK
        PANEL_BEG_OFFSET = 0
        IF (CB_IS_LR.AND.NROWS_TO_STACK.GT.0) THEN
          CALL DMUMPS_BLR_RETRIEVE_CB_LRB(
     &                  IW(ISTCHK+XXF), CB_LRB)
          CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(
     &                  IW(ISTCHK+XXF), BEGS_BLR)
          NB_BLR_ROWS = size(BEGS_BLR) - 1
          CALL DMUMPS_BLR_RETRIEVE_NB_PANELS(IW(ISTCHK+XXF), 
     &                    NB_BLR_SHIFT)
          PANEL2DECOMPRESS = -1
          DO II=NB_BLR_SHIFT+1,NB_BLR_ROWS
            IF (BEGS_BLR(II+1)-1-NASS.GT.
     &                NROWS_ALREADY_STACKED+NBROW(I)-1) THEN
              PANEL2DECOMPRESS = II
              EXIT
            ENDIF
          ENDDO
          IF (PANEL2DECOMPRESS.EQ.-1) THEN
            write(*,*) 'Internal error: PANEL2DECOMPRESS not found'
            CALL MUMPS_ABORT()
          ENDIF
          IF (KEEP(50).EQ.0) THEN
            NB_BLR_COLS = size(BEGS_BLR)  - 1
          ELSE
            NB_BLR_COLS = PANEL2DECOMPRESS
          ENDIF
          CURRENT_PANEL_SIZE = BEGS_BLR(PANEL2DECOMPRESS+1)
     &              - BEGS_BLR(PANEL2DECOMPRESS) 
          PANEL_BEG_OFFSET = NBROW(I) + NROWS_ALREADY_STACKED
     &              - BEGS_BLR(PANEL2DECOMPRESS) + NASS
          NROWS_TO_STACK_LOC = 
     &               min(NROWS_TO_STACK-NROWS_ALREADY_STACKED,
     &                   CURRENT_PANEL_SIZE-PANEL_BEG_OFFSET)
          LA_TEMP = CURRENT_PANEL_SIZE*NBCOLS
          allocate(A_TEMP(LA_TEMP),stat=allocok)
          IF (allocok.GT.0) THEN
             CALL MUMPS_SETI8TOI4(LA_TEMP,IERROR)
             IFLAG  = -13
            GOTO 700
          ENDIF
#if defined(BLR_MT)             
!$OMP PARALLEL             
#endif
          CALL DMUMPS_DECOMPRESS_PANEL(A_TEMP, LA_TEMP, 1_8,
     &          NBCOLS, NBCOLS, .TRUE., 1, 1, 
     &          NB_BLR_COLS-NB_BLR_SHIFT, 
     &          CB_LRB(PANEL2DECOMPRESS-NB_BLR_SHIFT,
     &                 1:NB_BLR_COLS-NB_BLR_SHIFT),
     &          0, 'V', 5,
     &          CBASM_TOFIX_IN=.TRUE.,
     &          ONLY_NELIM_IN=CURRENT_PANEL_SIZE-PANEL_BEG_OFFSET)
#if defined(BLR_MT)             
!$OMP END PARALLEL             
#endif
        ENDIF
        CALL DMUMPS_DM_SET_DYNPTR(
     &                   IW(PTRIST(STEP(ISON))+XXS), 
     &                   A, LA,
     &                   PAMASTER(STEP(ISON)),
     &                   IW(PTRIST(STEP(ISON))+XXD),
     &                   IW(PTRIST(STEP(ISON))+XXR), 
     &                   SON_A, IACHK, RECSIZE )
        DO II = NROWS_ALREADY_STACKED+1,
     &              NROWS_ALREADY_STACKED+NROWS_TO_STACK_LOC
          IROW_SON=PERM_LOC(NBROW(I)+II-1)
          INDICE_PERE = MAP(IROW_SON)
          CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &         NASS_PERE,
     &         NFRONT_PERE - NASS_PERE,
     &         NSLAVES_PERE,
     &         INDICE_PERE,
     &         NOSLA,
     &         IPOS_IN_SLAVE )
          INDICE_PERE = IPOS_IN_SLAVE
          IF (PACKED_CB) THEN
            IF (NELIM.EQ.0) THEN
            POSROW = IACHK +
     &         int(IROW_SON,8)*int(IROW_SON-1,8)/2_8
            ELSE
            POSROW = IACHK +
     &         int(NELIM+IROW_SON,8)*int(NELIM+IROW_SON-1,8)/2_8
            ENDIF
          ELSE
            POSROW = IACHK +
     &             int(NELIM+IROW_SON-1,8)*int(NBCOLS,8)
          ENDIF
          IF (KEEP(50).NE.0) THEN
            NBCOLS_EFF = NELIM + IROW_SON
          ELSE
            NBCOLS_EFF = NBCOLS
          ENDIF
          INDICE_PERE_ARRAY_ARG(1) = INDICE_PERE
          IF (CB_IS_LR) THEN
            CALL DMUMPS_ASM_SLAVE_MASTER(N, INODE_PERE, IW, LIW, 
     &         A, LA, ISON, 1, NBCOLS_EFF, 
     &         INDICE_PERE_ARRAY_ARG,
     &         A_TEMP(1+(II+PANEL_BEG_OFFSET
     &                -NROWS_ALREADY_STACKED-1)*NBCOLS), 
     &         PTLUST, PTRAST,
     &         STEP, PIMASTER, OPASSW, IWPOSCB, 
     &         MYID, KEEP,KEEP8,.FALSE.,NBCOLS)
          ELSE
            CALL DMUMPS_ASM_SLAVE_MASTER(N, INODE_PERE, IW, LIW, 
     &         A, LA, ISON, 1, NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &         SON_A(POSROW), PTLUST, PTRAST,
     &         STEP, PIMASTER, OPASSW, IWPOSCB, 
     &         MYID, KEEP,KEEP8,.FALSE.,NBCOLS_EFF)
          ENDIF
        ENDDO
        IF (CB_IS_LR.AND.NROWS_TO_STACK.GT.0) THEN
          deallocate(A_TEMP)
          NROWS_ALREADY_STACKED = NROWS_ALREADY_STACKED 
     &                  + NROWS_TO_STACK_LOC
          IF (NROWS_ALREADY_STACKED.LT.NROWS_TO_STACK) THEN
            GOTO 100
          ENDIF
        ENDIF
        IF (KEEP(219).NE.0) THEN
         IF(NSLAVES_PERE.GT.0 .AND. KEEP(50).EQ.2) THEN
          IF (CB_IS_LR) THEN
             CALL DMUMPS_BLR_RETRIEVE_M_ARRAY (
     &            IW(ISTCHK+XXF), M_ARRAY)
              M_ARRAY_RETRIEVED = .TRUE.
          ELSE
           IF (PACKED_CB) THEN
             POSROW = IACHK
     &          + int(NELIM+NBROW(1),8)*int(NELIM+NBROW(1)-1,8)/2_8
             ASIZE  = int(LMAP_LOC+NELIM,8)*int(NELIM+LMAP_LOC+1,8)/2_8
     &          - int(NELIM+NBROW(1),8)*int(NELIM+NBROW(1)-1,8)/2_8
           ELSE
             POSROW = IACHK +
     &                 int(NELIM+NBROW(1)-1,8)*int(NBCOLS,8)
             ASIZE  = int(LMAP_LOC-NBROW(1)+1,8) * int(NBCOLS,8)
           ENDIF
           CALL DMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
           IF (IERR .NE.0) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     &    ": PB allocation MAX_ARRAY during DMUMPS_MAPLIG_FILS_NIV1"
              IFLAG=-13
              IERROR=NFS4FATHER
              GOTO 700
           ENDIF
           IF  ( LMAP_LOC-NBROW(1)+1-KEEP(253)-NVSCHUR.GT. 0 ) THEN
           CALL DMUMPS_COMPUTE_MAXPERCOL(
     &          SON_A(POSROW),ASIZE,NBCOLS,
     &          LMAP_LOC-NBROW(1)+1-KEEP(253)-NVSCHUR,
     &          BUF_MAX_ARRAY,NFS4FATHER,PACKED_CB,
     &          NELIM+NBROW(1))
           ELSE
                CALL DMUMPS_SETMAXTOZERO(BUF_MAX_ARRAY,
     &          NFS4FATHER)
           ENDIF
           M_ARRAY => BUF_MAX_ARRAY
           M_ARRAY_RETRIEVED = .FALSE.
          ENDIF 
          CALL DMUMPS_ASM_MAX(N, INODE_PERE, IW, LIW, 
     &          A, LA, ISON, NFS4FATHER,
     &          M_ARRAY(1), PTLUST, PTRAST,
     &          STEP, PIMASTER, OPASSW,
     &          IWPOSCB,MYID, KEEP,KEEP8)
           IF ( M_ARRAY_RETRIEVED ) 
     &     CALL DMUMPS_BLR_FREE_M_ARRAY ( IW(ISTCHK+XXF) )
         ENDIF
        ENDIF 
          IF (IW(PTRIST(STEP(ISON))+XXNBPR) .EQ. 0
     &       ) THEN
               ISTCHK_LOC = PIMASTER(STEP(ISON))
               SAME_PROC= ISTCHK_LOC .LT. IWPOSCB
               IF (SAME_PROC) THEN
                 CALL DMUMPS_RESTORE_INDICES(N, ISON, INODE_PERE,
     &            IWPOSCB, PIMASTER, PTLUST, IW, LIW, STEP,
     &            KEEP,KEEP8)
               ENDIF
          ENDIF
          IF ( IW(PTLUST(STEP(INODE_PERE))+XXNBPR) .EQ. 0
     &       ) THEN
            CALL DMUMPS_INSERT_POOL_N( N, IPOOL, LPOOL,
     &        PROCNODE_STEPS,
     &        SLAVEF, KEEP(199), KEEP(28), KEEP(76), KEEP(80),
     &        KEEP(47), STEP, INODE_PERE+N )
            IF (KEEP(47) .GE. 3) THEN
              CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &       IPOOL, LPOOL, 
     &       PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &       MYID, STEP, N, ND, FILS )
            ENDIF
          END IF
      DO I = 0, NSLAVES_PERE
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
           NBROWS_ALREADY_SENT = 0
 95        CONTINUE
           NFRONT = IW(PIMASTER(STEP(ISON))+KEEP(IXSZ))
           NELIM  = IW(PIMASTER(STEP(ISON))+1+KEEP(IXSZ))
           DESCLU = .TRUE.
           IF (I == NSLAVES_PERE) THEN
             NROWS_TO_SEND=LMAP_LOC-NBROW(I)+1
           ELSE
             NROWS_TO_SEND=NBROW(I+1)-NBROW(I)
           ENDIF
           IF ( NROWS_TO_SEND .EQ. 0) CYCLE
           ITYPE_SON = MUMPS_TYPENODE( PROCNODE_STEPS(STEP(ISON)),
     &                                 KEEP(199) )
           IF (CB_IS_LR) THEN
             CALL DMUMPS_BUF_SEND_CONTRIB_TYPE2(NBROWS_ALREADY_SENT,
     &               DESCLU, INODE_PERE,
     &               NFRONT_PERE, NASS_PERE, NFS4FATHER, 
     &                    NSLAVES_PERE,
     &               ISON, NROWS_TO_SEND, LMAP_LOC,
     &               MAP, PERM_LOC(min(LMAP_LOC,NBROW(I))),
     &               IW(PIMASTER(STEP(ISON))),
     &               ADummy, 1_8,
     &               I, PDEST, PDEST_MASTER, COMM, IERR,
     &               KEEP,KEEP8, STEP, N, SLAVEF,
     &               ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &               PACKED_CB, KEEP(253), NVSCHUR,
     &               ITYPE_SON, MYID,
     &               NPIV_CHECK = IW(PTLUST(STEP(ISON))+3+KEEP(IXSZ)))
           ELSE
             CALL DMUMPS_DM_SET_DYNPTR(
     &                   IW(PTRIST(STEP(ISON))+XXS), 
     &                   A, LA,
     &                   PAMASTER(STEP(ISON)),
     &                   IW(PTRIST(STEP(ISON))+XXD),
     &                   IW(PTRIST(STEP(ISON))+XXR), 
     &                   SON_A, IACHK, RECSIZE )
             CALL DMUMPS_BUF_SEND_CONTRIB_TYPE2(NBROWS_ALREADY_SENT,
     &               DESCLU, INODE_PERE,
     &               NFRONT_PERE, NASS_PERE, NFS4FATHER, 
     &                    NSLAVES_PERE,
     &               ISON, NROWS_TO_SEND, LMAP_LOC,
     &               MAP, PERM_LOC(min(LMAP_LOC,NBROW(I))),
     &               IW(PIMASTER(STEP(ISON))),
     &               SON_A(IACHK:IACHK+RECSIZE-1_8), 
     &               RECSIZE,
     &               I, PDEST, PDEST_MASTER, COMM, IERR,
     &
     &               KEEP,KEEP8, STEP, N, SLAVEF,
     &               ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &               PACKED_CB, KEEP(253), NVSCHUR,
     &               ITYPE_SON, MYID)
           ENDIF
            IF ( IERR .EQ. -2 ) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     &": FAILURE, SEND BUFFER TOO SMALL DURING DMUMPS_MAPLIG_FILS_NIV1"
              IFLAG  = -17
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &        NROWS_TO_SEND *  KEEP( 35 )
              GO TO 700
            END IF
            IF ( IERR .EQ. -3 ) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     &": FAILURE, RECV BUFFER TOO SMALL DURING DMUMPS_MAPLIG_FILS_NIV1"
              IFLAG  = -20
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &        NROWS_TO_SEND *  KEEP( 35 )
              GO TO 700
            ENDIF
            IF (KEEP(219).NE.0) THEN
             IF ( IERR .EQ. -4 ) THEN
               IFLAG  = -13
               IERROR = BUF_LMAX_ARRAY
              IF (LP > 0) WRITE(LP,*) MYID,
     &": FAILURE, MAX_ARRAY ALLOC FAILED DURING DMUMPS_MAPLIG_FILS_NIV1"
               GO TO 700
             ENDIF
            ENDIF
            IF ( IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV = .TRUE.
              MESSAGE_RECEIVED = .FALSE.
              CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &          ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &          MPI_ANY_SOURCE, MPI_ANY_TAG,
     &          STATUS, 
     &          BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &          IWPOS, IWPOSCB, IPTRLU,
     &          LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &          PTLUST, PTRFAC,
     &          PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &          IFLAG, IERROR, COMM,
     &          PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &          root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &          FILS, DAD, PTRARW, PTRAIW,
     &          INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,
     &          LPTRAR, NELT, FRTPTR, FRTELT, 
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &          )
              IF ( IFLAG .LT. 0 ) GOTO 600
              GO TO 95
            END IF
        END IF
      END DO
      ISTCHK = PTRIST(STEP(ISON))
      PTRIST(STEP( ISON )) = -77777777
      IF ( IW(ISTCHK+KEEP(IXSZ)) .GE. 0 ) THEN
         WRITE(*,*) 'error 3 in DMUMPS_MAPLIG_FILS_NIV1'
         CALL MUMPS_ABORT()
      ENDIF
      CALL MUMPS_GETI8(DYNSIZE,IW(ISTCHK+XXD))
      CALL DMUMPS_FREE_BLOCK_CB_STATIC(.FALSE., MYID, N, ISTCHK,
     &     IW, LIW, LRLU, LRLUS, IPTRLU,
     &     IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &     )
      IF (DYNSIZE .GT. 0_8) THEN
        CALL DMUMPS_DM_FREE_BLOCK( SON_A, DYNSIZE,
     &                             KEEP(405).EQ.1, KEEP8 )
      ENDIF
      GOTO 600 
 700  CONTINUE
      CALL DMUMPS_BDC_ERROR(MYID, SLAVEF, COMM, KEEP )
 600  CONTINUE
      IF (CB_IS_LR) THEN
        CALL DMUMPS_BLR_FREE_CB_LRB(IW(ISTCHK+XXF),
     &             .FALSE., 
     &             KEEP8)
        IF ((KEEP(486).EQ.3).OR.KEEP(486).EQ.0) THEN
         CALL DMUMPS_BLR_END_FRONT(IW(ISTCHK+XXF), IFLAG, KEEP8)
        ENDIF
      ENDIF
      IF (allocated(NBROW))       DEALLOCATE(NBROW)
      IF (allocated(MAP))         DEALLOCATE(MAP)
      IF (allocated(PERM_LOC))        DEALLOCATE(PERM_LOC)
      IF (allocated(SLAVES_PERE)) DEALLOCATE(SLAVES_PERE)
      RETURN
      END SUBROUTINE DMUMPS_MAPLIG_FILS_NIV1
      SUBROUTINE DMUMPS_LOCAL_ASSEMBLY_TYPE2(I, PDEST, MYID,
     &           PDEST_MASTER, ISON, IFATH, NSLAVES_PERE, NASS_PERE,
     &           NFRONT_PERE, NFS4FATHER, LMAP_LOC, MAP,
     &           NBROW, PERM, IS_ofType5or6, IFLAG, IERROR,
     &           N, SLAVEF, KEEP, 
     &           IPOOL, LPOOL, STEP,
     &           PROCNODE_STEPS, COMM_LOAD, ISTEP_TO_INIV2,
     &           TAB_POS_IN_PERE,
     &           KEEP8, IW, LIW, A, LA, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &           PTRIST, PTLUST, PTRAST, PAMASTER, PIMASTER, ND,
     &           NELT, FRTPTR, FRTELT,
     &           OPASSW, OPELIW,
     &           ITLOC, RHS_MUMPS, KEEP253_LOC, NVSCHUR,
     &           FILS, DAD,
     &           LPTRAR, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &           SON_NIV, LRGROUPS)
      USE DMUMPS_BUF, ONLY: DMUMPS_BUF_MAX_ARRAY_MINSIZE,
     &                              BUF_MAX_ARRAY
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_LR_DATA_M
      USE DMUMPS_FAC_LR, ONLY: DMUMPS_DECOMPRESS_PANEL
      USE DMUMPS_LOAD, ONLY : DMUMPS_LOAD_POOL_UPD_NEW_POOL
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR
     &                 , DMUMPS_DM_SET_PTR, DMUMPS_DM_FREE_BLOCK
      IMPLICIT NONE
      INTEGER ICNTL(60)
      INTEGER, intent(in) :: I, PDEST, MYID, PDEST_MASTER, IFATH, ISON
      INTEGER, intent(in) :: N, SLAVEF
      INTEGER, intent(in) :: NSLAVES_PERE, NASS_PERE, NFRONT_PERE
      INTEGER, intent(in) :: NFS4FATHER
      INTEGER, intent(in) :: KEEP(500), STEP(N)
      INTEGER, intent(in) :: LMAP_LOC
      INTEGER, intent(in) :: NBROW(0:NSLAVES_PERE)
      INTEGER, intent(in) :: MAP(LMAP_LOC), PERM(LMAP_LOC)
      INTEGER, intent(inout) :: IFLAG, IERROR
      INTEGER(8), intent(in) :: KEEP8(150)
      INTEGER, intent(in) :: LIW, NELT, LPTRAR
      INTEGER(8), intent(in) :: LA
      INTEGER(8), intent(inout) :: IPTRLU, LRLU, LRLUS
      INTEGER, intent(inout) :: IWPOSCB
      INTEGER, intent(inout) :: IW(LIW)
      DOUBLE PRECISION, intent(inout) :: A( LA )
      INTEGER(8) :: PTRAST(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER    :: PTRIST(KEEP(28)), PIMASTER(KEEP(28)), ND(KEEP(28))
      INTEGER    :: PTLUST(KEEP(28))
      INTEGER, intent(inout) :: ITLOC(N)
      INTEGER, intent(in) :: FRTPTR( N+1 ), FRTELT( NELT )
      DOUBLE PRECISION, intent(inout) :: OPASSW, OPELIW
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER, intent(in) :: KEEP253_LOC, NVSCHUR
      INTEGER, intent(in) :: FILS(N), DAD( KEEP(28) )
      INTEGER(8), intent(in) :: PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER, intent(in) :: PROCNODE_STEPS( KEEP(28) ), COMM_LOAD
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION DBLARR(KEEP8(26))
      INTEGER INTARR(KEEP8(27))
      INTEGER LPOOL
      INTEGER IPOOL( LPOOL )
      LOGICAL, intent(in) :: IS_ofType5or6
      INTEGER, intent(in) :: SON_NIV
      INTEGER, intent(in) :: LRGROUPS(N)
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mpif.h'
      INTEGER    :: ISTCHK, ISTCHK_LOC, NBCOLS,
     &              NROW, NPIV, NSLSON,
     &              NFRONT, LDA_SON, NROWS_TO_STACK, II, INDICE_PERE,
     &              NOSLA, COLLIST, IPOS_IN_SLAVE, IROW_SON, ITMP,
     &              NBCOLS_EFF, DECR, NELIM
      LOGICAL    :: PACKED_CB, SAME_PROC
      INTEGER(8) :: SIZFR, POSROW, SHIFTCB_SON
      INTEGER(8) :: IACHK
      INTEGER :: SON_XXS
      DOUBLE PRECISION, DIMENSION(:), POINTER :: SON_A
      DOUBLE PRECISION, DIMENSION(:), POINTER :: SON_A_MASTER
      INTEGER(8) :: DYN_SIZE
      INTEGER    :: IERR, LP
      INTEGER INDICE_PERE_ARRAY_ARG(1)
      INTEGER :: INBPROCFILS_SON
      LOGICAL :: CB_IS_LR
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: M_ARRAY
      LOGICAL :: M_ARRAY_RETRIEVED
      INTEGER(8) :: POSELT
      INTEGER    :: IOLDPS, PARPIV_T1
      LOGICAL    :: LR_ACTIVATED
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_ROW, BEGS_BLR_COL,
     &                    BEGS_BLR_STA
      INTEGER :: NB_BLR_COLS, NB_BLR_ROWS,
     &           NB_COL_SHIFT, PANEL2DECOMPRESS, 
     &           CURRENT_PANEL_SIZE, PANEL_BEG_OFFSET,
     &           allocok, NROWS_ALREADY_STACKED, NROWS_TO_STACK_LOC,
     &           NB_ROW_SHIFT, NASS_SHIFT, NCOL_SHIFT, NROW_SHIFT
      INTEGER(8) :: LA_TEMP
      DOUBLE PRECISION, ALLOCATABLE :: A_TEMP(:)
      TYPE (LRB_TYPE), POINTER :: CB_LRB(:,:)
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
            IF (I == NSLAVES_PERE) THEN
              NROWS_TO_STACK = LMAP_LOC - NBROW(I) + 1
            ELSE
              NROWS_TO_STACK = NBROW(I+1) - NBROW(I)
            ENDIF
            DECR = 1
            IF ( MYID .EQ. PDEST_MASTER ) THEN
              IW(PTLUST(STEP(IFATH))+XXNBPR) =
     &            IW(PTLUST(STEP(IFATH))+XXNBPR) - DECR
              IF ( PDEST .EQ. PDEST_MASTER .AND. DECR .NE. 0) THEN
                IW(PIMASTER(STEP(ISON))+XXNBPR) =
     &             IW(PIMASTER(STEP(ISON))+XXNBPR) - DECR
              ENDIF
            ENDIF
            ISTCHK     = PTRIST(STEP(ISON))
            NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
            NROW       = IW(ISTCHK+2+KEEP(IXSZ))
            NPIV       = IW(ISTCHK+3+KEEP(IXSZ))
            NSLSON     = IW(ISTCHK+5+KEEP(IXSZ))
            NFRONT     = NPIV + NBCOLS
            SON_XXS    = IW(ISTCHK+XXS)
            PACKED_CB = ( SON_XXS .EQ. S_CB1COMP )
            CALL DMUMPS_DM_SET_DYNPTR(
     &          SON_XXS, 
     &          A, LA,
     &          PTRAST(STEP(ISON)),
     &          IW(PTRIST(STEP(ISON))+XXD),
     &          IW(PTRIST(STEP(ISON))+XXR), 
     &          SON_A, IACHK, SIZFR)
            CB_IS_LR = (IW(ISTCHK+XXLR).EQ.1 .OR.
     &              IW(ISTCHK+XXLR).EQ.3)
            NELIM = -9999
            IF (CB_IS_LR.AND.(SON_NIV.EQ.1).AND.
     &          KEEP(50).NE.0) THEN
             ISTCHK_LOC = PTLUST(STEP(ISON))
             NELIM  = IW(ISTCHK_LOC+1+KEEP(IXSZ))
             NPIV = IW(ISTCHK_LOC+3+KEEP(IXSZ))
             NFRONT = IW(ISTCHK_LOC+2+KEEP(IXSZ))
             NROW = NFRONT - NPIV
             NFRONT = NBCOLS
             NPIV = 0
            ENDIF
            IF (CB_IS_LR) THEN
             LDA_SON = NBCOLS
             SHIFTCB_SON = -9999
            ELSE
             IF (SON_XXS.EQ.S_NOLCBCONTIG ) THEN
               LDA_SON     = NBCOLS
               SHIFTCB_SON = int(NPIV,8)*int(NROW,8)
             ELSE IF (IW(ISTCHK+XXS).EQ.S_NOLCLEANED) THEN
               LDA_SON     = NBCOLS
               SHIFTCB_SON = 0_8
             ELSE
               LDA_SON     = NFRONT
               SHIFTCB_SON = int(NPIV,8)
             ENDIF
            ENDIF
            IF (PDEST .NE. PDEST_MASTER) THEN
                IF ( KEEP(55) .eq. 0 ) THEN
                  CALL DMUMPS_ASM_SLAVE_TO_SLAVE_INIT
     &            (N, IFATH, IW, LIW,
     &            A, LA, NROWS_TO_STACK, NBCOLS,
     &            OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &            ITLOC, RHS_MUMPS,
     &            FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &            KEEP,KEEP8, MYID, LRGROUPS )
                ELSE
                  CALL DMUMPS_ELT_ASM_S_2_S_INIT(NELT, FRTPTR, FRTELT,
     &            N, IFATH, IW, LIW,
     &            A, LA, NROWS_TO_STACK, NBCOLS, 
     &            OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &            ITLOC, RHS_MUMPS,
     &            FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &            KEEP, KEEP8, MYID, LRGROUPS )
                ENDIF
            ENDIF
            NROWS_ALREADY_STACKED = 0
 100        CONTINUE   
            NROWS_TO_STACK_LOC = NROWS_TO_STACK
            PANEL_BEG_OFFSET = 0
            IF (CB_IS_LR.AND.NROWS_TO_STACK.GT.0) THEN
              CALL DMUMPS_BLR_RETRIEVE_CB_LRB(
     &                  IW(ISTCHK+XXF), CB_LRB)
              IF (SON_NIV.EQ.1) THEN
                CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(
     &                  IW(ISTCHK+XXF), BEGS_BLR_ROW)
                CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_DYN(
     &                  IW(ISTCHK+XXF), BEGS_BLR_COL)
                NB_BLR_ROWS = size(BEGS_BLR_ROW) - 1
                CALL DMUMPS_BLR_RETRIEVE_NB_PANELS(IW(ISTCHK+XXF), 
     &                    NB_COL_SHIFT)
                NB_ROW_SHIFT = NB_COL_SHIFT
                NASS_SHIFT   = BEGS_BLR_ROW(NB_ROW_SHIFT+1)-1
              ELSE 
                CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(
     &                      IW(ISTCHK+XXF), BEGS_BLR_STA)
                NB_BLR_ROWS = size(BEGS_BLR_STA) - 2
                BEGS_BLR_ROW => BEGS_BLR_STA(2:NB_BLR_ROWS+2)
                CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_C(
     &                     IW(ISTCHK+XXF), BEGS_BLR_COL,
     &                     NB_COL_SHIFT)
                NB_ROW_SHIFT = 0
                NASS_SHIFT = 0
              ENDIF
              PANEL2DECOMPRESS = -1
              DO II=NB_ROW_SHIFT+1,NB_BLR_ROWS
                IF (BEGS_BLR_ROW(II+1)-1-NASS_SHIFT.GT.
     &                    NROWS_ALREADY_STACKED+NBROW(I)-1) THEN
                  PANEL2DECOMPRESS = II
                  EXIT
                ENDIF
              ENDDO
              IF (PANEL2DECOMPRESS.EQ.-1) THEN
                write(*,*) 'Internal error: PANEL2DECOMPRESS not found'
                CALL MUMPS_ABORT()
              ENDIF
              IF (KEEP(50).EQ.0) THEN
                NB_BLR_COLS = size(BEGS_BLR_COL)  - 1
              ELSEIF (SON_NIV.EQ.1) THEN
                NB_BLR_COLS = PANEL2DECOMPRESS
              ELSE
                NB_BLR_COLS = -1
                NCOL_SHIFT = NPIV
                NROW_SHIFT = NBCOLS-NROW
                DO II=NB_COL_SHIFT+1,size(BEGS_BLR_COL)-1
                  IF (BEGS_BLR_COL(II+1)-NCOL_SHIFT.GT.
     &              BEGS_BLR_ROW(PANEL2DECOMPRESS+1)-1+NROW_SHIFT) THEN
                    NB_BLR_COLS = II
                    EXIT
                  ENDIF
                ENDDO
                IF (NB_BLR_COLS.EQ.-1) THEN
                  write(*,*) 'Internal error: NB_BLR_COLS not found'
                  CALL MUMPS_ABORT()
                ENDIF
              ENDIF
              CURRENT_PANEL_SIZE = BEGS_BLR_ROW(PANEL2DECOMPRESS+1)
     &                     - BEGS_BLR_ROW(PANEL2DECOMPRESS) 
              PANEL_BEG_OFFSET = NBROW(I) + NROWS_ALREADY_STACKED
     &                 - BEGS_BLR_ROW(PANEL2DECOMPRESS) + NASS_SHIFT
              NROWS_TO_STACK_LOC =
     &               min(NROWS_TO_STACK-NROWS_ALREADY_STACKED,
     &                   CURRENT_PANEL_SIZE-PANEL_BEG_OFFSET) 
              LA_TEMP = CURRENT_PANEL_SIZE*NBCOLS
              allocate(A_TEMP(LA_TEMP),stat=allocok)
              IF (allocok.GT.0) THEN
                CALL MUMPS_SETI8TOI4(LA_TEMP,IERROR)
                IFLAG  = -13
                RETURN
              ENDIF
#if defined(BLR_MT)             
!$OMP PARALLEL             
#endif
              CALL DMUMPS_DECOMPRESS_PANEL(A_TEMP, LA_TEMP, 1_8,
     &              NBCOLS, NBCOLS, .TRUE., 1, 1, 
     &              NB_BLR_COLS-NB_COL_SHIFT, 
     &              CB_LRB(PANEL2DECOMPRESS-NB_ROW_SHIFT,
     &              1:NB_BLR_COLS-NB_COL_SHIFT),
     &              0, 'V', 6,
     &              CBASM_TOFIX_IN=.TRUE.,
     &              ONLY_NELIM_IN=CURRENT_PANEL_SIZE-PANEL_BEG_OFFSET)
#if defined(BLR_MT)             
!$OMP END PARALLEL             
#endif
            ENDIF
            DO II = NROWS_ALREADY_STACKED+1,
     &              NROWS_ALREADY_STACKED+NROWS_TO_STACK_LOC
              IROW_SON = PERM(NBROW(I)+II-1)
              INDICE_PERE=MAP(IROW_SON)
              CALL MUMPS_BLOC2_GET_ISLAVE(
     &        KEEP,KEEP8, IFATH, STEP, N, SLAVEF,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &        NASS_PERE,
     &        NFRONT_PERE - NASS_PERE,
     &        NSLAVES_PERE,
     &        INDICE_PERE,
     &        NOSLA,
     &        IPOS_IN_SLAVE )
              INDICE_PERE = IPOS_IN_SLAVE
              IF ( PACKED_CB ) THEN
                IF (NBCOLS - NROW .EQ. 0 ) THEN
                  ITMP = IROW_SON 
                  POSROW = IACHK+
     &                     int(ITMP,8) * int(ITMP-1,8) / 2_8
                ELSE
                  ITMP = IROW_SON + NBCOLS - NROW
                  POSROW = IACHK
     &               + int(ITMP,8) * int(ITMP-1,8) / 2_8
     &               - int(NBCOLS-NROW,8) * int(NBCOLS-NROW+1,8)/2_8
                ENDIF
              ELSE 
                POSROW = IACHK + SHIFTCB_SON
     &               +int(IROW_SON-1,8)*int(LDA_SON,8)
              ENDIF
              IF (PDEST == PDEST_MASTER) THEN
                 IF (KEEP(50).NE.0) THEN
                   NBCOLS_EFF = IROW_SON + NBCOLS - NROW
                 ELSE
                   NBCOLS_EFF = NBCOLS
                 ENDIF
                 INDICE_PERE_ARRAY_ARG(1) = INDICE_PERE
                 IF ((IS_ofType5or6).AND.(KEEP(50).EQ.0)) THEN
                    IF (CB_IS_LR) THEN
                     write(*,*) 'Compress CB + Type5or6 fronts not',
     &                       'coded yet!!!'
                     CALL MUMPS_ABORT()
                   ENDIF
                   CALL DMUMPS_ASM_SLAVE_MASTER(N, IFATH, IW, LIW, 
     &             A, LA, ISON, NROWS_TO_STACK, NBCOLS_EFF, 
     &             INDICE_PERE_ARRAY_ARG,
     &             SON_A(POSROW), PTLUST, PTRAST,
     &             STEP, PIMASTER, OPASSW,
     &             IWPOSCB, MYID, KEEP,KEEP8,
     &             IS_ofType5or6, LDA_SON
     &             )
                   EXIT  
                 ELSE IF ( (KEEP(50).NE.0) .AND. 
     &              (.NOT.PACKED_CB).AND.(IS_ofType5or6) ) THEN
                    IF (CB_IS_LR) THEN
                     write(*,*) 'Compress CB + Type5or6 fronts not',
     &                       'coded yet!!!'
                     CALL MUMPS_ABORT()
                   ENDIF
                   CALL DMUMPS_ASM_SLAVE_MASTER(N, IFATH, IW, LIW, 
     &             A, LA, ISON, NROWS_TO_STACK,
     &             NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &             SON_A(POSROW), PTLUST, PTRAST,
     &             STEP, PIMASTER, OPASSW,
     &             IWPOSCB, MYID, KEEP,KEEP8,
     &             IS_ofType5or6, LDA_SON
     &)
                   EXIT
                 ELSE
                   IF (CB_IS_LR) THEN
                     CALL DMUMPS_ASM_SLAVE_MASTER(N, IFATH, IW, LIW, 
     &                 A, LA, ISON, 1, NBCOLS_EFF, 
     &                 INDICE_PERE_ARRAY_ARG,
     &                 A_TEMP(1+(II+PANEL_BEG_OFFSET
     &                        -NROWS_ALREADY_STACKED-1)*NBCOLS), 
     &                 PTLUST, PTRAST,
     &                 STEP, PIMASTER, OPASSW,
     &                 IWPOSCB, MYID, KEEP,KEEP8,
     &                 IS_ofType5or6, NBCOLS )
                   ELSE
                     CALL DMUMPS_ASM_SLAVE_MASTER(N, IFATH, IW, LIW, 
     &                 A, LA, ISON, 1, NBCOLS_EFF, 
     &                 INDICE_PERE_ARRAY_ARG,
     &                 SON_A(POSROW), PTLUST, PTRAST,
     &                 STEP, PIMASTER, OPASSW,
     &                 IWPOSCB, MYID, KEEP,KEEP8,
     &                 IS_ofType5or6, LDA_SON )
                   ENDIF
                 ENDIF
              ELSE
                ISTCHK  = PTRIST(STEP(ISON))
                COLLIST = ISTCHK + 6 + KEEP(IXSZ) 
     &                   + IW( ISTCHK + 5 +KEEP(IXSZ)) + NROW + NPIV
                IF (CB_IS_LR.AND.(SON_NIV.EQ.1).AND.
     &              KEEP(50).NE.0) THEN
                  ISTCHK_LOC = PTLUST(STEP(ISON))
                  COLLIST = ISTCHK_LOC + 6 + KEEP(IXSZ) 
     &                   + IW( ISTCHK + 5 +KEEP(IXSZ)) 
     &                   + IW(ISTCHK_LOC+2+KEEP(IXSZ)) 
     &                   + IW(ISTCHK_LOC+3+KEEP(IXSZ)) 
                ENDIF
                 IF (KEEP(50).NE.0) THEN
                   NBCOLS_EFF = IROW_SON + NBCOLS - NROW
                   IF (CB_IS_LR.AND.SON_NIV.EQ.1) 
     &             NBCOLS_EFF = IROW_SON + NBCOLS - (NROW-NELIM)
                 ELSE
                   NBCOLS_EFF = NBCOLS
                 ENDIF
                 INDICE_PERE_ARRAY_ARG(1) = INDICE_PERE
                 IF ( (IS_ofType5or6) .AND.
     &                 (
     &                  ( KEEP(50).EQ.0) 
     &                    .OR. 
     &                  ( (KEEP(50).NE.0).and. (.NOT.PACKED_CB) )
     &                 )
     &               ) THEN
                    IF (CB_IS_LR) THEN
                     write(*,*) 'Compress CB + Type5or6 fronts not',
     &                       'coded yet!!!'
                     CALL MUMPS_ABORT()
                   ENDIF
                   CALL DMUMPS_ASM_SLAVE_TO_SLAVE(N, IFATH,
     &             IW, LIW,
     &             A, LA, NROWS_TO_STACK, NBCOLS, 
     &             INDICE_PERE_ARRAY_ARG,
     &             IW( COLLIST ), SON_A(POSROW),
     &             OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &             ITLOC, RHS_MUMPS,
     &             FILS, ICNTL, KEEP,KEEP8,
     &             MYID, IS_ofType5or6, LDA_SON)
                   IW( PTRIST(STEP(IFATH))+XXNBPR) =
     &               IW( PTRIST(STEP(IFATH))+XXNBPR) - NROWS_TO_STACK
                   EXIT
                 ELSE
                   IF (CB_IS_LR) THEN
                     CALL DMUMPS_ASM_SLAVE_TO_SLAVE(N, IFATH,
     &                 IW, LIW,
     &                 A, LA, 1, NBCOLS_EFF, 
     &                 INDICE_PERE_ARRAY_ARG,
     &                 IW( COLLIST ), 
     &                 A_TEMP(1+(II+PANEL_BEG_OFFSET
     &                        -NROWS_ALREADY_STACKED-1)*NBCOLS), 
     &                 OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &                 ITLOC, RHS_MUMPS,
     &                 FILS, ICNTL, KEEP,KEEP8,
     &                 MYID, IS_ofType5or6, NBCOLS)
                   ELSE
                     CALL DMUMPS_ASM_SLAVE_TO_SLAVE(N, IFATH,
     &                 IW, LIW,
     &                 A, LA, 1, NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &                 IW( COLLIST ), SON_A(POSROW),
     &                 OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &                 ITLOC, RHS_MUMPS,
     &                 FILS, ICNTL, KEEP,KEEP8,
     &                 MYID, IS_ofType5or6, LDA_SON)
                   ENDIF
                   IW( PTRIST(STEP(IFATH))+XXNBPR) =
     &                 IW( PTRIST(STEP(IFATH))+XXNBPR) - 1
                 ENDIF
              ENDIF
            ENDDO
            IF (CB_IS_LR.AND.NROWS_TO_STACK.GT.0) THEN
              deallocate(A_TEMP)
              NROWS_ALREADY_STACKED = NROWS_ALREADY_STACKED 
     &                              + NROWS_TO_STACK_LOC
              IF (NROWS_ALREADY_STACKED.LT.NROWS_TO_STACK) THEN
                GOTO 100
              ENDIF
            ENDIF
            IF (PDEST.EQ.PDEST_MASTER) THEN 
             IF (KEEP(219).NE.0) THEN
               IF(NSLAVES_PERE.GT.0 .AND. KEEP(50).EQ.2) THEN
                IF (CB_IS_LR) THEN
                 CALL DMUMPS_BLR_RETRIEVE_M_ARRAY (
     &            IW(ISTCHK+XXF), M_ARRAY)
                 M_ARRAY_RETRIEVED = .TRUE.
                ELSE
                  IF (PACKED_CB) THEN
                    WRITE(*,*) "Error 1 in PARPIV/DMUMPS_MAPLIG"
                    CALL MUMPS_ABORT()
                  ELSE
                    POSROW = IACHK + SHIFTCB_SON+
     &                       int(NBROW(1)-1,8)*int(LDA_SON,8)
                  ENDIF
                  CALL DMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
                  IF (IERR .NE.0) THEN
                    IF (LP .GT. 0) THEN
                      WRITE(LP, *) "MAX_ARRAY allocation failed"
                    ENDIF
                    IFLAG=-13
                    IERROR=NFS4FATHER
                    RETURN
                  ENDIF
                  ITMP=-9999
                  IF (LMAP_LOC-NBROW(1)+1-KEEP253_LOC-NVSCHUR.NE.0)
     &            THEN
                  CALL DMUMPS_COMPUTE_MAXPERCOL(
     &                 SON_A(POSROW),
     &       SIZFR-SHIFTCB_SON-int(NBROW(1)-1,8)*int(LDA_SON,8),
     &                 LDA_SON, 
     &                 LMAP_LOC-NBROW(1)+1-KEEP253_LOC-NVSCHUR,
     &                 BUF_MAX_ARRAY,NFS4FATHER,PACKED_CB,ITMP)
                  ELSE
                       CALL DMUMPS_SETMAXTOZERO(
     &                 BUF_MAX_ARRAY, NFS4FATHER)
                  ENDIF
                  M_ARRAY => BUF_MAX_ARRAY(1:size(BUF_MAX_ARRAY))
                  M_ARRAY_RETRIEVED = .FALSE.
                ENDIF
                CALL DMUMPS_ASM_MAX(N, IFATH, IW, LIW, 
     &                 A, LA, ISON, NFS4FATHER,
     &                 M_ARRAY(1), PTLUST, PTRAST,
     &                 STEP, PIMASTER,
     &                 OPASSW,IWPOSCB,MYID, KEEP,KEEP8)
                IF ( M_ARRAY_RETRIEVED ) 
     &          CALL DMUMPS_BLR_FREE_M_ARRAY ( IW(ISTCHK+XXF) )
               ENDIF
             ENDIF 
             ISTCHK_LOC = PIMASTER(STEP(ISON))
               SAME_PROC= ISTCHK_LOC .LT. IWPOSCB
               IF ( SAME_PROC ) THEN
                 INBPROCFILS_SON = PTRIST(STEP(ISON))+XXNBPR
                 WRITE(*,*)
     &           "Internal error 0 in DMUMPS_LOCAL_ASSEMBLY_TYPE2",
     &           INBPROCFILS_SON, PIMASTER(STEP(ISON))
                 CALL MUMPS_ABORT()
               ELSE
                 INBPROCFILS_SON = PIMASTER(STEP(ISON))+XXNBPR
               ENDIF
               IF ( IW(INBPROCFILS_SON) .EQ. 0 ) THEN
                 IF (SAME_PROC) THEN
                   CALL DMUMPS_RESTORE_INDICES(N, ISON, IFATH,
     &               IWPOSCB, PIMASTER, PTLUST, IW, LIW, STEP,
     &               KEEP,KEEP8)
                 ENDIF
                 IF (SAME_PROC) THEN
                   ISTCHK_LOC = PTRIST(STEP(ISON))
                   PTRIST(STEP( ISON) ) = -99999999
                 ELSE
                   PIMASTER(STEP( ISON )) = -99999999
                 ENDIF
                 CALL MUMPS_GETI8(DYN_SIZE, IW(ISTCHK_LOC+XXD))
                 IF (DYN_SIZE .GT. 0_8) THEN
                   CALL DMUMPS_DM_SET_PTR( PAMASTER(STEP(ISON)),
     &             DYN_SIZE, SON_A_MASTER )
                 ENDIF
                 CALL DMUMPS_FREE_BLOCK_CB_STATIC(.FALSE., MYID, N,
     &            ISTCHK_LOC,
     &            IW, LIW, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &            LA, KEEP,KEEP8, .FALSE.
     &            )
                 IF (DYN_SIZE .GT. 0_8) THEN
                   CALL DMUMPS_DM_FREE_BLOCK( SON_A_MASTER, DYN_SIZE,
     &                                        KEEP(405).EQ.1, KEEP8 )
                 ENDIF
               ENDIF
             IF ( IW(PTLUST(STEP(IFATH))+XXNBPR) .EQ. 0
     &       ) THEN
             IOLDPS = PTLUST(STEP(IFATH))  
             IF (NSLAVES_PERE.EQ.0) THEN 
               POSELT = PTRAST(STEP(IFATH))
               PARPIV_T1 = -999 
               LR_ACTIVATED = (IW(IOLDPS+XXLR).GT.0)
               CALL DMUMPS_PARPIVT1_SET_NVSCHUR_and_MAX (
     &           N, IFATH, IW, LIW, A, LA, KEEP, PERM,
     &            IOLDPS, POSELT, 
     &            NFRONT_PERE, NASS_PERE, LR_ACTIVATED, PARPIV_T1)
             ENDIF
               CALL DMUMPS_INSERT_POOL_N( N, IPOOL, LPOOL,
     &           PROCNODE_STEPS,
     &           SLAVEF, KEEP(199), KEEP(28), KEEP(76), KEEP(80),
     &           KEEP(47), STEP, IFATH+N )
               IF (KEEP(47) .GE. 3) THEN
                 CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &          IPOOL, LPOOL, 
     &          PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &          MYID, STEP, N, ND, FILS )
               ENDIF
             END IF
            ELSE
               CALL DMUMPS_ASM_SLAVE_TO_SLAVE_END
     &         (N, IFATH, IW, LIW,
     &         NBROW(I), STEP, PTRIST, ITLOC, RHS_MUMPS,
     &         KEEP,KEEP8)
            END IF
      RETURN
      END SUBROUTINE DMUMPS_LOCAL_ASSEMBLY_TYPE2
