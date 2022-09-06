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
      SUBROUTINE DMUMPS_PROCESS_CONTRIB_TYPE2( COMM_LOAD, ASS_IRECV, 
     &   MSGLEN, BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, POSFAC,
     &   N, IW, LIW, A, LA, PTRIST, PTLUST, PTRFAC, PTRAST,
     &   STEP, PIMASTER, PAMASTER, PERM,
     &   COMP, root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, NSTK_S,
     &   FILS, DAD, PTRARW, PTRAIW, INTARR, DBLARR, NBFIN,
     &   MYID, COMM, ICNTL, KEEP,KEEP8,DKEEP, IFLAG, IERROR,
     &   IPOOL, LPOOL, LEAF, ND, FRERE_STEPS, LPTRAR, NELT,
     &   FRTPTR, FRTELT, 
     &   ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &               , LRGROUPS
     &   )
      USE DMUMPS_LOAD
      USE DMUMPS_BUF
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_FAC_LR, ONLY: DMUMPS_DECOMPRESS_PANEL
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR,
     &                                    DMUMPS_DM_FREE_BLOCK
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 60 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(230)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV, MSGLEN
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: LRLU, IPTRLU, LRLUS, LA, POSFAC
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER NBFIN
      INTEGER COMP
      INTEGER NELT, LPTRAR
      INTEGER PROCNODE_STEPS( KEEP(28) ), PTRIST(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER PTLUST( KEEP(28) )
      INTEGER PERM(N)
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(N)
      INTEGER ITLOC( N + KEEP(253)), NSTK_S( KEEP(28) )
      INTEGER :: FILS( N ), DAD(KEEP(28))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER ND(KEEP(28)), FRERE_STEPS( KEEP(28) )
      INTEGER(8), INTENT(IN) :: PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER INTARR( KEEP8(27) )
      DOUBLE PRECISION DBLARR( KEEP8(26) )
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER COMM, MYID, IFLAG, IERROR
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER FRTPTR(N+1), FRTELT( NELT )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER NFS4FATHER
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MUMPS_PROCNODE, MUMPS_TYPESPLIT
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPESPLIT
      INTEGER IERR
      INTEGER NBROWS_ALREADY_SENT, NBROWS_PACKET
      INTEGER I, INODE, ISON, POSITION, NBROW, LROW, IROW, INDCOL
      INTEGER LREQI
      INTEGER(8) :: LREQA, POSCONTRIB
      INTEGER ROW_LENGTH
      INTEGER MASTER
      INTEGER ISTCHK
      LOGICAL SAME_PROC
      LOGICAL SLAVE_NODE
      LOGICAL IS_ofType5or6
      INTEGER ISHIFT_BUFR, LBUFR_LOC, LBUFR_BYTES_LOC
      INTEGER TYPESPLIT
      INTEGER DECR
      INTEGER :: INBPROCFILS_SON 
      LOGICAL :: CB_IS_LR
      INTEGER :: CB_IS_LR_INT, NB_BLR_COLS, allocok,
     &           NBROWS_PACKET_2PACK, PANEL_BEG_OFFSET
      INTEGER(8) :: LA_TEMP
      DOUBLE PRECISION, ALLOCATABLE :: A_TEMP(:)
      TYPE (LRB_TYPE), POINTER :: LRB
      TYPE (LRB_TYPE), ALLOCATABLE, TARGET :: BLR_CB(:)
      INTEGER(8) :: IACHK, SIZFR8, DYN_SIZE
      DOUBLE PRECISION, DIMENSION(:), POINTER :: DYNPTR
      INTEGER    :: NSLAVES, NFRONT, NASS1, IOLDPS, PARPIV_T1
      LOGICAL    :: LR_ACTIVATED
      INTEGER(8) :: POSELT
      INCLUDE 'mumps_headers.h'
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, ISON, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NBROW, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, LROW, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NBROWS_ALREADY_SENT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NBROWS_PACKET, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 CB_IS_LR_INT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CB_IS_LR = (CB_IS_LR_INT.EQ.1)
      MASTER     = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE)),
     &                             KEEP(199) )
      SLAVE_NODE = MASTER .NE. MYID
      TYPESPLIT = MUMPS_TYPESPLIT( PROCNODE_STEPS(STEP(INODE)),
     &                             KEEP(199) )
      IS_ofType5or6 = ((TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6))
      IF (SLAVE_NODE .AND. PTRIST(STEP(INODE)) ==0) THEN
        ISHIFT_BUFR     = ( MSGLEN + KEEP(34) ) / KEEP(34)
        LBUFR_LOC       = LBUFR - ISHIFT_BUFR + 1
        LBUFR_BYTES_LOC = LBUFR_LOC * KEEP(34)
          CALL DMUMPS_TREAT_DESCBAND( INODE, COMM_LOAD, ASS_IRECV,
     &     BUFR(ISHIFT_BUFR), LBUFR_LOC, LBUFR_BYTES_LOC,
     &     PROCNODE_STEPS, POSFAC,
     &     IWPOS, IWPOSCB, IPTRLU,
     &     LRLU, LRLUS, N, IW, LIW, A, LA,
     &     PTRIST, PTLUST, PTRFAC,
     &     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &     IFLAG, IERROR, COMM,
     &     PERM, IPOOL, LPOOL, LEAF,
     &     NBFIN, MYID, SLAVEF,
     &
     &     root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD,
     &     PTRARW, PTRAIW,
     &     INTARR, DBLARR, ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &     LPTRAR, NELT, FRTPTR, FRTELT, 
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &     )
          IF (IFLAG.LT.0) RETURN
      ENDIF
      IF ( SLAVE_NODE ) THEN
         LREQI = LROW + NBROWS_PACKET
      ELSE
         LREQI = NBROWS_PACKET
      END IF
         LREQA = int(LROW,8)
         CALL DMUMPS_GET_SIZE_NEEDED(
     &           LREQI, LREQA, .FALSE.,
     &           KEEP(1), KEEP8(1), 
     &           N, KEEP(28), IW, LIW, A, LA,
     &           LRLU, IPTRLU,
     &           IWPOS, IWPOSCB, PTRIST, PTRAST,
     &           STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &           KEEP(IXSZ), COMP, DKEEP(97),
     &           MYID, SLAVEF, PROCNODE_STEPS, DAD, 
     &           IFLAG, IERROR )
         IF (IFLAG.LT.0) THEN
            CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
            RETURN
         ENDIF
         LRLU  = LRLU - LREQA
         LRLUS = LRLUS - LREQA
         POSCONTRIB = POSFAC
         POSFAC = POSFAC + LREQA
         KEEP8(67) = min(LRLUS, KEEP8(67))
         KEEP8(69) = KEEP8(69) + LREQA
         KEEP8(68) = max(KEEP8(69), KEEP8(68))
         CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &        LA-LRLUS,0_8,LREQA,KEEP,KEEP8,LRLUS)
         IF  ( SLAVE_NODE ) THEN
            IROW   = IWPOS
            INDCOL = IWPOS + NBROWS_PACKET
         ELSE
            IROW   = IWPOS
            INDCOL = -1
         END IF
         IWPOS = IWPOS + LREQI
         IF ( SLAVE_NODE ) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           IW( INDCOL ), LROW, MPI_INTEGER,
     &           COMM, IERR )
         END IF
         DO I = 1, NBROWS_PACKET
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           IW( IROW + I - 1 ), 1, MPI_INTEGER,
     &           COMM, IERR )
         END DO
         IF ( SLAVE_NODE ) THEN
            IF ( NBROWS_ALREADY_SENT + NBROWS_PACKET == NBROW ) THEN
              IW(PTRIST(STEP(INODE))+XXNBPR) =
     &        IW(PTRIST(STEP(INODE))+XXNBPR) - NBROW
            ENDIF
            IF ( KEEP(55) .eq. 0 ) THEN               
               CALL DMUMPS_ASM_SLAVE_TO_SLAVE_INIT
     &              (N, INODE, IW, LIW, A, LA,
     &              NBROW, LROW,
     &              OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &              ITLOC, RHS_MUMPS,
     &              FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &              KEEP,KEEP8, MYID, LRGROUPS )
            ELSE
               CALL DMUMPS_ELT_ASM_S_2_S_INIT(
     &              NELT, FRTPTR, FRTELT,
     &              N, INODE, IW, LIW, A, LA,
     &              NBROW, LROW,
     &              OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &              ITLOC, RHS_MUMPS,
     &              FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &              KEEP,KEEP8, MYID, LRGROUPS )
            ENDIF
            IF (CB_IS_LR) THEN
              CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &             NB_BLR_COLS, 1,
     &             MPI_INTEGER, COMM, IERR )
              CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &             PANEL_BEG_OFFSET, 1,
     &             MPI_INTEGER, COMM, IERR )
              allocate(BLR_CB(NB_BLR_COLS),stat=allocok)
              IF (allocok.GT.0) THEN
                IERROR = NB_BLR_COLS
                IFLAG  = -13
                CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
                RETURN
              ENDIF
              DO I=1,NB_BLR_COLS
                LRB => BLR_CB(I)
                CALL DMUMPS_MPI_UNPACK_LRB(BUFR, LBUFR, 
     &             LBUFR_BYTES, POSITION, LRB, KEEP8, 
     &             COMM, IFLAG, IERROR)
              ENDDO
              NBROWS_PACKET_2PACK = max(NBROWS_PACKET,BLR_CB(1)%M)
              LA_TEMP = NBROWS_PACKET_2PACK*LROW
              allocate(A_TEMP(LA_TEMP),stat=allocok)
              IF (allocok.GT.0) THEN
                CALL MUMPS_SETI8TOI4(LA_TEMP,IERROR)
                IFLAG  = -13
                CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
                RETURN
              ENDIF
#if defined(BLR_MT)             
!$OMP PARALLEL             
#endif
              CALL DMUMPS_DECOMPRESS_PANEL(A_TEMP, LA_TEMP, 1_8,
     &              LROW, LROW, .TRUE., 1, 1, 
     &              NB_BLR_COLS, BLR_CB, 0, 'V', 3,
     &              CBASM_TOFIX_IN=.TRUE.,
     &              ONLY_NELIM_IN=NBROWS_PACKET_2PACK-PANEL_BEG_OFFSET)
#if defined(BLR_MT)             
!$OMP END PARALLEL             
#endif
              DO I=1,NBROWS_PACKET
                IF (KEEP(50).EQ.0) THEN
                  ROW_LENGTH = LROW
                ELSE
                  CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ROW_LENGTH,
     &                 1,
     &                 MPI_INTEGER,
     &                 COMM, IERR )
                ENDIF
                CALL DMUMPS_ASM_SLAVE_TO_SLAVE(N, INODE, IW, LIW, A, LA,
     &              1, ROW_LENGTH, IW( IROW+I-1 ),IW(INDCOL),
     &              A_TEMP(1+(I-1+PANEL_BEG_OFFSET)*LROW),
     &              OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &              ITLOC, RHS_MUMPS,
     &              FILS, ICNTL, KEEP,KEEP8, MYID, IS_ofType5or6, 
     &              LROW)
              ENDDO
              CALL DEALLOC_BLR_PANEL(BLR_CB, NB_BLR_COLS, KEEP8)
              deallocate(A_TEMP, BLR_CB)
              GOTO 200
            ENDIF
            DO I=1,NBROWS_PACKET
               IF(KEEP(50).NE.0)THEN
                  CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ROW_LENGTH,
     &                 1,
     &                 MPI_INTEGER,
     &                 COMM, IERR )
               ELSE
                 ROW_LENGTH=LROW
               ENDIF
               CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              A(POSCONTRIB),
     &              ROW_LENGTH,
     &              MPI_DOUBLE_PRECISION,
     &              COMM, IERR )
               CALL DMUMPS_ASM_SLAVE_TO_SLAVE(N, INODE, IW, LIW, A, LA,
     &              1, ROW_LENGTH, IW( IROW+I-1 ),IW(INDCOL),
     &              A(POSCONTRIB),
     &              OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &              ITLOC, RHS_MUMPS,
     &              FILS, ICNTL, KEEP,KEEP8, MYID, IS_ofType5or6, 
     &              ROW_LENGTH )
            ENDDO
 200        CONTINUE            
            CALL DMUMPS_ASM_SLAVE_TO_SLAVE_END
     &           (N, INODE, IW, LIW,
     &           NBROWS_PACKET, STEP, PTRIST,
     &           ITLOC, RHS_MUMPS,KEEP,KEEP8)
         ELSE
            IF (CB_IS_LR) THEN
              CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &             NB_BLR_COLS, 1,
     &             MPI_INTEGER, COMM, IERR )
              CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &             PANEL_BEG_OFFSET, 1,
     &             MPI_INTEGER, COMM, IERR )
              allocate(BLR_CB(NB_BLR_COLS),stat=allocok)
              IF (allocok.GT.0) THEN
                IERROR = NB_BLR_COLS
                IFLAG  = -13
                CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
                RETURN
              ENDIF
              DO I=1,NB_BLR_COLS
                LRB => BLR_CB(I)
                CALL DMUMPS_MPI_UNPACK_LRB(BUFR, LBUFR, 
     &             LBUFR_BYTES, POSITION, LRB, KEEP8, 
     &             COMM, IFLAG, IERROR)
              ENDDO
              NBROWS_PACKET_2PACK = max(NBROWS_PACKET,BLR_CB(1)%M)
              LA_TEMP = NBROWS_PACKET_2PACK*LROW
              allocate(A_TEMP(LA_TEMP),stat=allocok)
              IF (allocok.GT.0) THEN
                CALL MUMPS_SETI8TOI4(LA_TEMP,IERROR)
                IFLAG  = -13
                CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
                RETURN
              ENDIF
#if defined(BLR_MT)             
!$OMP PARALLEL             
#endif
              CALL DMUMPS_DECOMPRESS_PANEL(A_TEMP, LA_TEMP, 1_8,
     &              LROW, LROW, .TRUE., 1, 1, 
     &              NB_BLR_COLS, BLR_CB, 0, 'V', 4,
     &              CBASM_TOFIX_IN=.TRUE., 
     &              ONLY_NELIM_IN=NBROWS_PACKET_2PACK-PANEL_BEG_OFFSET)
#if defined(BLR_MT)             
!$OMP END PARALLEL             
#endif
              DO I=1,NBROWS_PACKET
                IF(KEEP(50).NE.0)THEN
                  CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ROW_LENGTH,
     &                 1,
     &                 MPI_INTEGER,
     &                 COMM, IERR )
                ELSE
                  ROW_LENGTH=LROW
                ENDIF
                CALL DMUMPS_ASM_SLAVE_MASTER(N, INODE, IW, LIW, A, LA,
     &              ISON, 1, ROW_LENGTH, IW( IROW+I-1 ),
     &              A_TEMP(1+(I-1+PANEL_BEG_OFFSET)*LROW),
     &              PTLUST, PTRAST,
     &              STEP, PIMASTER, OPASSW,
     &              IWPOSCB, MYID, KEEP,KEEP8,
     &              IS_ofType5or6, LROW
     &          )
              ENDDO
              CALL DEALLOC_BLR_PANEL(BLR_CB, NB_BLR_COLS, KEEP8)
              deallocate(A_TEMP, BLR_CB)
              GOTO 300
            ENDIF
            DO I=1,NBROWS_PACKET
               IF(KEEP(50).NE.0)THEN
                  CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 ROW_LENGTH,
     &                 1,
     &                 MPI_INTEGER,
     &                 COMM, IERR )
               ELSE
                 ROW_LENGTH=LROW
               ENDIF
               CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              A(POSCONTRIB),
     &              ROW_LENGTH,
     &              MPI_DOUBLE_PRECISION,
     &              COMM, IERR )
               CALL DMUMPS_ASM_SLAVE_MASTER(N, INODE, IW, LIW, A, LA,
     &              ISON, 1, ROW_LENGTH, IW( IROW +I-1 ),
     &              A(POSCONTRIB), PTLUST, PTRAST,
     &              STEP, PIMASTER, OPASSW,
     &              IWPOSCB, MYID, KEEP,KEEP8,
     &              IS_ofType5or6, ROW_LENGTH
     &)
            ENDDO
 300      CONTINUE
          IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          IF (KEEP(219).NE.0) THEN
            IF(KEEP(50) .EQ. 2) THEN
               CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              NFS4FATHER,
     &              1,
     &              MPI_INTEGER,
     &              COMM, IERR )
               IF(NFS4FATHER .GT. 0) THEN
                  CALL DMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
                  IF (IERR .NE. 0) THEN
                      IERROR         = BUF_LMAX_ARRAY
                      IFLAG          = -13
                      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
                      RETURN
                  ENDIF
                  CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 BUF_MAX_ARRAY,
     &                 NFS4FATHER,
     &                 MPI_DOUBLE_PRECISION,
     &                 COMM, IERR )
                  CALL DMUMPS_ASM_MAX(N, INODE, IW, LIW, A, LA,
     &                 ISON, NFS4FATHER,
     &                 BUF_MAX_ARRAY, PTLUST, PTRAST,
     &                 STEP, PIMASTER, OPASSW,
     &                 IWPOSCB, MYID, KEEP,KEEP8)
               ENDIF
            ENDIF
          ENDIF
          ENDIF
          IF (NBROWS_ALREADY_SENT + NBROWS_PACKET == NBROW ) THEN
            DECR = 1
            ISTCHK = PIMASTER(STEP(ISON))
            SAME_PROC = ISTCHK .LT. IWPOSCB
            IW(PTLUST(STEP(INODE))+XXNBPR) =
     &      IW(PTLUST(STEP(INODE))+XXNBPR) - DECR
            IF (SAME_PROC) THEN
              INBPROCFILS_SON = PTRIST(STEP(ISON))+XXNBPR
            ELSE
              INBPROCFILS_SON = PIMASTER(STEP(ISON))+XXNBPR
            ENDIF
            IW(INBPROCFILS_SON) = IW(INBPROCFILS_SON) - DECR
            IF ( IW(INBPROCFILS_SON) .EQ. 0 ) THEN
               IF (SAME_PROC) THEN
                  CALL DMUMPS_RESTORE_INDICES(N, ISON, INODE, IWPOSCB,
     &                 PIMASTER, PTLUST, IW, LIW, STEP, KEEP,KEEP8)
               ENDIF
               IF (SAME_PROC) THEN
                  ISTCHK = PTRIST(STEP(ISON))
                  PTRIST(STEP( ISON) ) = -99999999
               ELSE
                  PIMASTER(STEP( ISON )) = -99999999
               ENDIF
               CALL DMUMPS_DM_SET_DYNPTR( IW(ISTCHK+XXS), A, LA,
     &              PAMASTER(STEP(ISON)), IW(ISTCHK+XXD),
     &              IW(ISTCHK+XXR), DYNPTR, IACHK, SIZFR8)
               CALL MUMPS_GETI8(DYN_SIZE, IW(ISTCHK+XXD))
               CALL DMUMPS_FREE_BLOCK_CB_STATIC(
     &              .FALSE., MYID, N, ISTCHK,
     &              IW, LIW, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &              LA, KEEP,KEEP8, .FALSE.
     &              )
               IF ( DYN_SIZE .GT. 0_8 ) THEN
                 CALL DMUMPS_DM_FREE_BLOCK( DYNPTR, DYN_SIZE,
     &                                      KEEP(405).EQ.1, KEEP8 )
               ENDIF
            ENDIF
            IF (IW(PTLUST(STEP(INODE))+XXNBPR) .EQ. 0) THEN
             IOLDPS = PTLUST(STEP(INODE))
             NSLAVES= IW(IOLDPS+5+KEEP(IXSZ))
             IF (NSLAVES.EQ.0) THEN 
               NFRONT = IW(IOLDPS+KEEP(IXSZ))
               NASS1  = iabs(IW(IOLDPS + 2+KEEP(IXSZ)))
               POSELT = PTRAST(STEP(INODE))
               PARPIV_T1 = -999 
               LR_ACTIVATED = (IW(IOLDPS+XXLR).GT.0)
               CALL DMUMPS_PARPIVT1_SET_NVSCHUR_and_MAX (
     &           N, INODE, IW, LIW, A, LA, KEEP, PERM,
     &            IOLDPS, POSELT, 
     &            NFRONT, NASS1, LR_ACTIVATED, PARPIV_T1)
             ENDIF
               CALL DMUMPS_INSERT_POOL_N( N, IPOOL, LPOOL,
     &              PROCNODE_STEPS,
     &              SLAVEF, KEEP(199), KEEP(28), KEEP(76), KEEP(80),
     &              KEEP(47), STEP, INODE+N )
               IF (KEEP(47) .GE. 3) THEN
                  CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &          IPOOL, LPOOL, 
     &                 PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &                 MYID, STEP, N, ND, FILS )
               ENDIF
            ENDIF
          ENDIF 
      END IF 
         IWPOS = IWPOS - LREQI
         LRLU = LRLU + LREQA
         LRLUS = LRLUS + LREQA
         KEEP8(69) = KEEP8(69) - LREQA
         POSFAC = POSFAC - LREQA
         CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &        LA-LRLUS,0_8,-LREQA,KEEP,KEEP8,LRLUS)
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_CONTRIB_TYPE2
