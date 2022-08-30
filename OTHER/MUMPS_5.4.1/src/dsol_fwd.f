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
      SUBROUTINE DMUMPS_SOL_R(N, A, LA, IW, LIW, WCB, LWCB,
     &    NRHS,
     &    PTRICB, IWCB, LIWCB, 
     &    RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD, 
     &    STEP,
     &    FRERE, DAD, FILS,
     &    NSTK, IPOOL, LPOOL, PTRIST, PTRFAC, MYLEAF, MYROOT,
     &    INFO,
     &    KEEP, KEEP8, DKEEP,
     &    PROCNODE_STEPS,
     &    SLAVEF, COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,
     &    RHS_ROOT, LRHS_ROOT, MTYPE, 
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP
     &    )
      USE DMUMPS_STATIC_PTR_M, ONLY : DMUMPS_SET_STATIC_PTR,
     &                                DMUMPS_GET_TMP_PTR
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER(8), INTENT(IN) :: LA, LWCB
      INTEGER, INTENT(IN) :: N, LIW, LPOOL, LIWCB
      INTEGER, INTENT(IN) :: SLAVEF, MYLEAF, MYROOT, COMM, MYID
      INTEGER INFO( 80 ), KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER NRHS
      DOUBLE PRECISION A( LA ), WCB( LWCB )
      INTEGER(8), intent(in) :: LRHS_ROOT
      DOUBLE PRECISION RHS_ROOT( LRHS_ROOT )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER STEP( N ), FRERE( KEEP(28) ), FILS( N ),
     &        DAD( KEEP(28) )
      INTEGER NSTK(KEEP(28)), IPOOL( LPOOL )
      INTEGER PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PTRICB( KEEP(28) ) 
      LOGICAL, intent(in) :: DO_NBSPARSE
      INTEGER, intent(in) :: LRHS_BOUNDS
      INTEGER, intent(in) :: RHS_BOUNDS(LRHS_BOUNDS)
      INTEGER IW( LIW ), IWCB( LIWCB )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER, intent(in) ::  POSINRHSCOMP_FWD(N), LRHSCOMP 
      DOUBLE PRECISION, intent(inout) :: RHSCOMP(LRHSCOMP,NRHS)
      LOGICAL, intent(in) :: FROM_PP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY(1)
      LOGICAL FLAG
      DOUBLE PRECISION, DIMENSION(:), POINTER :: A_PTR
      INTEGER(8) :: LA_PTR
      INTEGER NBFIN, MYROOT_LEFT
      INTEGER POSIWCB
      INTEGER(8) :: POSWCB, PLEFTWCB
      INTEGER INODE, IFATH
      INTEGER III, LEAF
      LOGICAL BLOQ
      EXTERNAL MUMPS_PROCNODE
      INTEGER MUMPS_PROCNODE
      LOGICAL ERROR_WAS_BROADCASTED
      DUMMY(1) = 1
      KEEP(266)=0
      POSIWCB = LIWCB
      POSWCB  = LWCB
      PLEFTWCB= 1_8
      PTRICB = 0
      LEAF = MYLEAF + 1
      III    = 1
      NBFIN = SLAVEF
      MYROOT_LEFT = MYROOT
      IF ( MYROOT_LEFT .EQ. 0 ) THEN
        NBFIN = NBFIN - 1
        CALL DMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &       RACINE_SOLVE, SLAVEF, KEEP)
        IF (NBFIN.EQ.0) GOTO 260
      END IF
   50 CONTINUE
      IF (SLAVEF .EQ. 1) THEN
         CALL DMUMPS_GET_INODE_FROM_POOL
     &        ( IPOOL(1), LPOOL, III, LEAF, INODE,
     &          KEEP(208) )
        GOTO 60
      ENDIF
      BLOQ = ( ( III .EQ. LEAF )
     &     )
      CALL DMUMPS_SOLVE_RECV_AND_TREAT( BLOQ, FLAG,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MYID, SLAVEF, COMM,
     &     N, NRHS, IPOOL, LPOOL, LEAF,
     &     NBFIN, NSTK, IW, LIW, A, LA, PTRIST, PTRFAC,
     &     IWCB, LIWCB,
     &     WCB, LWCB, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP,
     &     PROCNODE_STEPS,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &     , FROM_PP
     &    )
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      IF (.not. FLAG) THEN
         IF (III .NE. LEAF) THEN
            CALL DMUMPS_GET_INODE_FROM_POOL
     &           (IPOOL(1), LPOOL, III, LEAF, INODE,
     &           KEEP(208) )
            GOTO 60
         ENDIF                  
      ENDIF                     
      GOTO 50
 60   CONTINUE
        CALL DMUMPS_SET_STATIC_PTR(A)
        CALL DMUMPS_GET_TMP_PTR(A_PTR)
        LA_PTR = LA
      CALL DMUMPS_SOLVE_NODE_FWD( INODE,
     &        huge(INODE), huge(INODE), 
     &        BUFR, LBUFR, LBUFR_BYTES,
     &        MYID, SLAVEF, COMM,  N,
     &        IPOOL, LPOOL, LEAF, NBFIN, NSTK,
     &        IWCB, LIWCB, WCB, LWCB, A_PTR(1), LA_PTR,
     &        IW, LIW, NRHS, 
     &        POSWCB, PLEFTWCB, POSIWCB,
     &        PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     &        FILS, STEP, FRERE, DAD,
     &        INFO, KEEP,KEEP8, DKEEP, RHS_ROOT, LRHS_ROOT, MTYPE,
     &        RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &        , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE 
     &        , FROM_PP
     &        , ERROR_WAS_BROADCASTED
     & )
      IF ( INFO(1) .LT. 0 ) THEN
        IF (.NOT. ERROR_WAS_BROADCASTED) THEN
          CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
        ENDIF
        GOTO 260
      ENDIF
      IFATH = DAD(STEP(INODE))
      IF ( IFATH .EQ. 0 ) THEN
        MYROOT_LEFT = MYROOT_LEFT - 1
        IF (MYROOT_LEFT .EQ. 0) THEN
          NBFIN = NBFIN - 1
          IF (SLAVEF .GT. 1) THEN
            CALL DMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID,
     &           COMM, RACINE_SOLVE, SLAVEF, KEEP)
          ENDIF
        END IF
      ELSE
        IF ( MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IFATH)), KEEP(199))
     &       .EQ. MYID ) THEN
           IF ( PTRICB(STEP(INODE)) .EQ. 1 .OR.
     &          PTRICB(STEP(INODE)) .EQ. -1 ) THEN
             NSTK(STEP(IFATH)) = NSTK(STEP(IFATH)) - 1
             IF (NSTK(STEP(IFATH)) .EQ. 0) THEN
               IPOOL(LEAF) = IFATH
               LEAF = LEAF + 1
               IF (LEAF .GT. LPOOL) THEN
                  WRITE(*,*)
     &            'Internal error DMUMPS_TRAITER_MESSAGE_SOLVE',
     &            LEAF, LPOOL
                  CALL MUMPS_ABORT()
               ENDIF
             ENDIF
             PTRICB(STEP(INODE)) = 0
           ENDIF
        ENDIF
      ENDIF
      IF ( NBFIN .EQ. 0 ) GOTO 260
      GOTO 50
  260 CONTINUE
      CALL DMUMPS_CLEAN_PENDING(INFO(1), KEEP, BUFR, LBUFR,LBUFR_BYTES,
     &     COMM, DUMMY(1),  
     &     SLAVEF, .TRUE., .FALSE.) 
      RETURN
      END SUBROUTINE DMUMPS_SOL_R
