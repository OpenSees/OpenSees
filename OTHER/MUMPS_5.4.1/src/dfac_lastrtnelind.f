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
      SUBROUTINE DMUMPS_LAST_RTNELIND( COMM_LOAD, ASS_IRECV, 
     &    root, FRERE, IROOT, 
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
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  
     &               , LRGROUPS
     &     )
      USE DMUMPS_BUF
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER IROOT
      INTEGER ICNTL( 60 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(230)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: POSFAC,IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB
      INTEGER(8) :: LA
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER, intent(in) :: LRGROUPS(N)
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER PERM(N)
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC( N+KEEP(253) ), FILS( N ), DAD( KEEP(28) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8), INTENT(IN) :: PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND(KEEP(28)), FRERE(KEEP(28))
      DOUBLE PRECISION DBLARR( KEEP8(26) )
      INTEGER INTARR( KEEP8(27) )
      INTEGER I, NELIM, NB_CONTRI_GLOBAL, NUMORG, 
     &        NFRONT, IROW, JCOL, PDEST, HF, IOLDPS,
     &        IN, DEB_ROW, ILOC_ROW, IFSON, ILOC_COL,
     &        IPOS_SON, NELIM_SON, NSLAVES_SON, HS,
     &        IROW_SON, ICOL_SON, ISLAVE, IERR, 
     &        NELIM_SENT, IPOS_STATREC, TYPE_SON
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mumps_tags.h'
      NB_CONTRI_GLOBAL = KEEP(41)
      NUMORG    = root%ROOT_SIZE
      NELIM     = KEEP(42)
      NFRONT    = NUMORG + KEEP(42)
      DO IROW = 0, root%NPROW - 1
        DO JCOL = 0, root%NPCOL - 1
            PDEST = IROW * root%NPCOL + JCOL
          IF ( PDEST .NE. MYID ) THEN
           CALL DMUMPS_BUF_SEND_ROOT2SLAVE(NFRONT, 
     &     NB_CONTRI_GLOBAL, PDEST, COMM, KEEP, IERR)
           if (IERR.lt.0) then
                write(6,*) ' error detected by ',
     &          'DMUMPS_BUF_SEND_ROOT2SLAVE'
                CALL MUMPS_ABORT()
               endif
           ENDIF
        END DO
      END DO
      CALL  DMUMPS_PROCESS_ROOT2SLAVE( NFRONT,
     &    NB_CONTRI_GLOBAL, root,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, COMM_LOAD,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    PTRARW, PTRAIW,
     &    INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND )
       IF (IFLAG < 0 ) RETURN
      HF = 6 + KEEP(IXSZ)
      IOLDPS = PTLUST_S(STEP(IROOT))
      IN = IROOT
      DEB_ROW = IOLDPS + HF
      ILOC_ROW    = DEB_ROW
      DO WHILE (IN.GT.0) 
       IW(ILOC_ROW)           = IN
       IW(ILOC_ROW+NFRONT)    = IN
       ILOC_ROW = ILOC_ROW + 1
       IN = FILS(IN)
      END DO
      IFSON = -IN
      ILOC_ROW    = IOLDPS + HF + NUMORG
      ILOC_COL    = ILOC_ROW + NFRONT
      IF ( NELIM.GT.0 ) THEN
        IN = IFSON
        DO WHILE (IN.GT.0)
          IPOS_SON  = PIMASTER(STEP(IN))
          IF (IPOS_SON .EQ. 0) GOTO 100
          NELIM_SON   = IW(IPOS_SON+1+KEEP(IXSZ))
          if (NELIM_SON.eq.0) then
                write(6,*) ' error 1 in process_last_rtnelind'
                CALL MUMPS_ABORT()
              endif
          NSLAVES_SON = IW(IPOS_SON+5+KEEP(IXSZ))
          HS          = 6 + NSLAVES_SON + KEEP(IXSZ)
          IROW_SON    = IPOS_SON + HS
          ICOL_SON    = IROW_SON + NELIM_SON
          DO I = 1, NELIM_SON
            IW( ILOC_ROW+I-1 ) = IW( IROW_SON+I-1 )
          ENDDO
          DO I = 1, NELIM_SON
            IW( ILOC_COL+I-1 ) = IW( ICOL_SON+I-1 )
          ENDDO
          NELIM_SENT = ILOC_ROW - IOLDPS - HF + 1
          DO ISLAVE = 0,NSLAVES_SON
            IF (ISLAVE.EQ.0) THEN
             PDEST= MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IN)),KEEP(199))
            ELSE
             PDEST = IW(IPOS_SON + 5 + ISLAVE+KEEP(IXSZ))
            ENDIF
            IF (PDEST.NE.MYID) THEN
             CALL DMUMPS_BUF_SEND_ROOT2SON(IN, NELIM_SENT,
     &        PDEST, COMM, KEEP, IERR )
             if (IERR.lt.0) then
                write(6,*) ' error detected by ',
     &          'DMUMPS_BUF_SEND_ROOT2SLAVE'
                CALL MUMPS_ABORT()
               endif
            ELSE
             CALL DMUMPS_PROCESS_ROOT2SON( COMM_LOAD, ASS_IRECV,
     &       IN, NELIM_SENT, root,
     &
     &       BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &       IWPOS, IWPOSCB, IPTRLU,
     &       LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST, 
     &       PTLUST_S, PTRFAC,
     &       PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &       IFLAG, IERROR, COMM,
     &       PERM,
     &       IPOOL, LPOOL, LEAF,
     &       NBFIN, MYID, SLAVEF,
     &
     &       OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &       FILS, DAD, PTRARW, PTRAIW,
     &       INTARR, DBLARR, ICNTL, KEEP, KEEP8, DKEEP, ND, FRERE,
     &       LPTRAR, NELT, FRTPTR, FRTELT, 
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &       )
             IF ( ISLAVE .NE. 0 ) THEN
               IF (KEEP(50) .EQ. 0) THEN
                IPOS_STATREC = PTRIST(STEP(IN))+6+KEEP(IXSZ)
               ELSE
                IPOS_STATREC = PTRIST(STEP(IN))+8+KEEP(IXSZ)
               ENDIF
               IF (IW(IPOS_STATREC).EQ. S_REC_CONTSTATIC) THEN
                  IW(IPOS_STATREC) = S_ROOT2SON_CALLED
               ELSE
                IF (NSLAVES_SON .EQ. 0) THEN
                  TYPE_SON = 1
                ELSE
                  TYPE_SON = 2
                ENDIF
                CALL DMUMPS_FREE_BAND( N, IN, PTRIST, PTRAST,
     &          IW, LIW, A, LA, LRLU, LRLUS, IWPOSCB,
     &          IPTRLU, STEP, MYID, KEEP, KEEP8, TYPE_SON
     &        )
               ENDIF
             ENDIF
             IPOS_SON  = PIMASTER(STEP(IN))
            ENDIF
          END DO
          CALL  DMUMPS_FREE_BLOCK_CB_STATIC(
     &       .FALSE., MYID, N, IPOS_SON,
     &       IW, LIW,
     &       LRLU, LRLUS, IPTRLU,
     &       IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &         )
          ILOC_ROW = ILOC_ROW + NELIM_SON
          ILOC_COL = ILOC_COL + NELIM_SON
 100      CONTINUE
          IN = FRERE(STEP(IN))
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_LAST_RTNELIND
