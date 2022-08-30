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
      SUBROUTINE DMUMPS_SOL_S(N, A, LA, IW, LIW, W, LWC,
     &    NRHS, 
     &    RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD,
     &    PTRICB, PTRACB, IWCB, LIWW, W2, 
     &    NE_STEPS, STEP,
     &    FRERE, DAD, FILS, IPOOL, LPOOL, PTRIST, PTRFAC, 
     &    MYLEAF, MYROOT, ICNTL, INFO, 
     &    PROCNODE_STEPS,
     &    SLAVEF, COMM,MYID, BUFR, LBUFR, LBUFR_BYTES,
     &    KEEP,KEEP8, DKEEP, RHS_ROOT, LRHS_ROOT, MTYPE, 
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, PANEL_POS, LPANEL_POS
     &    , PRUN_BELOW, TO_PROCESS, SIZE_TO_PROCESS
     &    , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP
     &    )
      USE DMUMPS_STATIC_PTR_M, ONLY : DMUMPS_SET_STATIC_PTR,
     &                                DMUMPS_GET_TMP_PTR
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER(8), intent(in) :: LA
      INTEGER(8), intent(in) :: LWC
      INTEGER, intent(in) :: N,LIW,LIWW,LPOOL
      INTEGER, intent(in) :: SLAVEF,MYLEAF,MYROOT,COMM,MYID
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER NE_STEPS(KEEP(28))
      INTEGER IPOOL(LPOOL)
      INTEGER LPANEL_POS
      INTEGER PANEL_POS(LPANEL_POS)
      INTEGER ICNTL(60), INFO(80)
      INTEGER PTRIST(KEEP(28)),
     &        PTRICB(KEEP(28))
      INTEGER(8) :: PTRACB(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER NRHS
      DOUBLE PRECISION A(LA), W(LWC)
      DOUBLE PRECISION W2(KEEP(133))
      INTEGER IW(LIW),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)),DAD(KEEP(28)),FILS(N)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER LRHSCOMP, POSINRHSCOMP_BWD(N)
      DOUBLE PRECISION RHSCOMP(LRHSCOMP,NRHS)
      INTEGER(8), intent(in) :: LRHS_ROOT
      DOUBLE PRECISION RHS_ROOT( LRHS_ROOT )
      LOGICAL, INTENT(in) :: PRUN_BELOW
      INTEGER, intent(in)           :: SIZE_TO_PROCESS
      LOGICAL, intent(in)           :: TO_PROCESS(SIZE_TO_PROCESS)
      LOGICAL, intent(in) :: DO_NBSPARSE
      INTEGER, intent(in) :: LRHS_BOUNDS
      INTEGER, intent(in) :: RHS_BOUNDS(LRHS_BOUNDS)
      LOGICAL, intent(in) :: FROM_PP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      LOGICAL FLAG
      DOUBLE PRECISION, DIMENSION(:), POINTER :: A_PTR
      INTEGER(8) :: LA_PTR
      INTEGER(8) :: POSWCB, PLEFTW
      INTEGER POSIWCB
      INTEGER NBFINF
      INTEGER INODE
      INTEGER III,IIPOOL,MYLEAF_LEFT
      LOGICAL BLOQ
      INTEGER DUMMY(1)
      LOGICAL :: ERROR_WAS_BROADCASTED, DO_MCAST2_TERMBWD
      LOGICAL :: ALLOW_OTHERS_TO_LEAVE
      LOGICAL, DIMENSION(:), ALLOCATABLE :: DEJA_SEND
      INTEGER :: allocok
      DUMMY(1)=0
      KEEP(266)=0
      ALLOCATE(DEJA_SEND( 0:SLAVEF-1 ), stat=allocok)
      if(allocok.ne.0) then
         WRITE(6,*) ' Allocation error of DEJA_SEND in '
     &        //'routine DMUMPS_SOL_S '
         INFO(1)=-13
         INFO(2)=SLAVEF
      endif
      CALL MUMPS_PROPINFO( ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT.0 ) GOTO 340
      PLEFTW = 1_8
      POSIWCB = LIWW
      POSWCB = LWC
      III = 1
      IIPOOL = MYROOT + 1
      MYLEAF_LEFT = MYLEAF
      NBFINF = SLAVEF
      ALLOW_OTHERS_TO_LEAVE = ( MYLEAF_LEFT .EQ. 0 .AND.
     &                          KEEP(31) .EQ. 0 )
      ALLOW_OTHERS_TO_LEAVE = ALLOW_OTHERS_TO_LEAVE .OR.
     &                        KEEP(31) .EQ. 1
      IF (ALLOW_OTHERS_TO_LEAVE) THEN
        CALL DMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID, COMM, TERMBWD,
     &                  SLAVEF, KEEP)
        NBFINF = NBFINF - 1
        IF (NBFINF .EQ. 0 .AND. MYLEAF_LEFT .EQ. 0) THEN
          GOTO 340
        ENDIF
      ENDIF
      ERROR_WAS_BROADCASTED = .FALSE.
      DO_MCAST2_TERMBWD = .FALSE.
      DO WHILE ( NBFINF .NE. 0 .OR. MYLEAF_LEFT .NE. 0 )
        BLOQ = (  III .EQ. IIPOOL  )
        CALL DMUMPS_BACKSLV_RECV_AND_TREAT( BLOQ, FLAG, BUFR, LBUFR,
     &     LBUFR_BYTES, MYID, SLAVEF, COMM,
     &     N, IWCB, LIWW, POSIWCB,
     &     W, LWC, POSWCB,
     &     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &     IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &     STEP,  FRERE, FILS, PROCNODE_STEPS,
     &     PLEFTW, KEEP,KEEP8, DKEEP,
     &     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAF_LEFT, 
     &     NRHS, MTYPE, 
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &     , PRUN_BELOW, TO_PROCESS, SIZE_TO_PROCESS
     &     , FROM_PP
     &     )
        IF ( INFO(1) .LT. 0 ) GOTO 340 
        IF ( .NOT. FLAG ) THEN
          IF (III .NE. IIPOOL) THEN
            INODE = IPOOL(IIPOOL-1)
            IIPOOL = IIPOOL - 1
              CALL DMUMPS_SET_STATIC_PTR(A)
              CALL DMUMPS_GET_TMP_PTR(A_PTR)
              LA_PTR = LA
            CALL DMUMPS_SOLVE_NODE_BWD( INODE, 
     &        N, IPOOL, LPOOL, IIPOOL, NBFINF,
     &        A_PTR(1), LA_PTR, IW, LIW, W, LWC, NRHS,
     &        POSWCB, PLEFTW, POSIWCB,
     &        RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD,
     &        PTRICB, PTRACB, IWCB, LIWW, W2, 
     &        NE_STEPS, STEP,
     &        FRERE, FILS, PTRIST, PTRFAC,
     &        MYLEAF_LEFT, INFO,
     &        PROCNODE_STEPS, DEJA_SEND,
     &        SLAVEF, COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,
     &        KEEP,KEEP8, DKEEP, RHS_ROOT, LRHS_ROOT, MTYPE,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, PANEL_POS, LPANEL_POS,
     &        PRUN_BELOW, TO_PROCESS, SIZE_TO_PROCESS
     &        , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP
     &        , ERROR_WAS_BROADCASTED
     &        , DO_MCAST2_TERMBWD
     &        )
            IF ( INFO(1) .LT. 0 ) THEN
              IF (.NOT. ERROR_WAS_BROADCASTED) THEN
                IF (NBFINF .EQ. 0 ) THEN
                  CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
                ENDIF
              ENDIF
            ENDIF
            IF (DO_MCAST2_TERMBWD) THEN
              CALL DMUMPS_MCAST2( DUMMY, 1,  MPI_INTEGER, MYID, COMM,
     &                            TERMBWD, SLAVEF, KEEP )
            ENDIF
          ENDIF                   
        END IF                    
      ENDDO
  340 CONTINUE
      IF (ALLOCATED(DEJA_SEND)) DEALLOCATE(DEJA_SEND)
      RETURN
      END SUBROUTINE DMUMPS_SOL_S
