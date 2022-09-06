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
      SUBROUTINE DMUMPS_SCATTER_DIST_RHS(
     &           NSLAVES, N,
     &           MYID_NODES, COMM_NODES,
     &           NRHS_COL, NRHS_loc, LRHS_loc,
     &           MAP_RHS_loc,
     &           IRHS_loc, RHS_loc, RHS_loc_size,
     &           RHSCOMP, LD_RHSCOMP,
     &           POSINRHSCOMP_FWD, NB_FS_IN_RHSCOMP,
     &           LSCAL, scaling_data_dr,
     &           LP, LPOK, KEEP, NB_BYTES_LOC, INFO )
      USE DMUMPS_STRUC_DEF
!$    USE OMP_LIB
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NSLAVES, N, MYID_NODES
      INTEGER, INTENT(IN)    :: NRHS_loc, LRHS_loc
      INTEGER, INTENT(IN)    :: NRHS_COL
      INTEGER, INTENT(IN)    :: COMM_NODES
      INTEGER, INTENT(IN)    :: MAP_RHS_loc(max(1,NRHS_loc))
      INTEGER, INTENT(IN)    :: IRHS_loc(NRHS_loc)
      INTEGER(8), INTENT(IN) :: RHS_loc_size
      DOUBLE PRECISION, INTENT(IN)    :: RHS_loc(RHS_loc_size)
      INTEGER, INTENT(IN)    :: NB_FS_IN_RHSCOMP, LD_RHSCOMP
      INTEGER, INTENT(IN)    :: POSINRHSCOMP_FWD(N)
      DOUBLE PRECISION, INTENT(OUT)   :: RHSCOMP(LD_RHSCOMP, NRHS_COL)
      INTEGER                :: KEEP(500)
      LOGICAL, INTENT(IN)    :: LSCAL
      type scaling_data_t
        SEQUENCE
        DOUBLE PRECISION, dimension(:), pointer :: SCALING
        DOUBLE PRECISION, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      type(scaling_data_t), INTENT(IN) :: scaling_data_dr
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, INTENT(IN)    :: LP      
      INTEGER, INTENT(INOUT) :: INFO(2) 
      INTEGER(8), INTENT(OUT):: NB_BYTES_LOC
      INCLUDE 'mpif.h'
      INTEGER :: IERR_MPI
!$    LOGICAL :: OMP_FLAG
!$    INTEGER :: CHUNK, NOMP
!$    INTEGER(8) :: CHUNK8
      INTEGER :: allocok
      INTEGER :: MAXRECORDS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROWSTOSEND 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NEXTROWTOSEND 
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: BUFR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFRECI 
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BUFRECR
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: IS_SEND_ACTIVE, TOUCHED
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MPI_REQI, MPI_REQR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IRHS_loc_sorted
      INTEGER :: Iloc  
      INTEGER :: Iloc_sorted 
      INTEGER :: IREQ  
      INTEGER :: IMAP, IPROC_MAX 
      INTEGER :: IFS   
      INTEGER :: MAX_ACTIVE_SENDS
      INTEGER :: NB_ACTIVE_SENDS
      INTEGER :: NB_FS_TOUCHED
      INTEGER :: NBROWSTORECV                
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
!$    NOMP = OMP_GET_MAX_THREADS()
      NB_BYTES_LOC = 0_8
      ALLOCATE( NBROWSTOSEND    (NSLAVES),
     &          NEXTROWTOSEND   (NSLAVES),
     &          IRHS_loc_sorted (NRHS_loc),
     &          stat=allocok )
      IF (allocok > 0) THEN
        INFO(1) = -13
        INFO(2) = NSLAVES+NSLAVES+NRHS_loc
      ENDIF
      NB_BYTES_LOC = int(2*NSLAVES+NRHS_loc,8)*KEEP(34)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, allocok, 1,
     &                    MPI_INTEGER, MPI_SUM,
     &                    COMM_NODES, IERR_MPI )
      IF (allocok .GT. 0) RETURN
      NBROWSTOSEND(1:NSLAVES) = 0
      DO Iloc = 1, NRHS_loc
        IF (IRHS_loc(Iloc) .GE. 1 .AND.
     &      IRHS_loc(Iloc) .LE. N) THEN
          IMAP = MAP_RHS_loc(Iloc)
          NBROWSTOSEND(IMAP+1) = NBROWSTOSEND(IMAP+1)+1
        ENDIF
      ENDDO
      NEXTROWTOSEND(1)=1 
      DO IMAP=1, NSLAVES-1
        NEXTROWTOSEND(IMAP+1)=NEXTROWTOSEND(IMAP)+NBROWSTOSEND(IMAP)
      ENDDO
      NBROWSTOSEND=0
      DO Iloc = 1, NRHS_loc
        IF (IRHS_loc(Iloc) .GE. 1 .AND.
     &      IRHS_loc(Iloc) .LE. N) THEN
          IMAP = MAP_RHS_loc(Iloc)
          Iloc_sorted = NEXTROWTOSEND(IMAP+1)+NBROWSTOSEND(IMAP+1)
          IRHS_loc_sorted(Iloc_sorted) = Iloc
          NBROWSTOSEND(IMAP+1)=NBROWSTOSEND(IMAP+1)+1
        ENDIF
      ENDDO
      CALL DMUMPS_DR_BUILD_NBROWSTORECV()
      MAX_ACTIVE_SENDS = min(10, NSLAVES)
      IF (KEEP(72) .EQ.1 ) THEN
        MAXRECORDS = 15
      ELSE
        MAXRECORDS = min(200000,2000000/NRHS_COL)
        MAXRECORDS = min(MAXRECORDS,
     &                50000000 / MAX_ACTIVE_SENDS / NRHS_COL)
        MAXRECORDS = max(MAXRECORDS, 50)
      ENDIF
      ALLOCATE(BUFR(MAXRECORDS*NRHS_COL,
     &                        MAX_ACTIVE_SENDS),   
     &         MPI_REQI(MAX_ACTIVE_SENDS),         
     &         MPI_REQR(MAX_ACTIVE_SENDS),         
     &         IS_SEND_ACTIVE(MAX_ACTIVE_SENDS),
     &         BUFRECI(MAXRECORDS),                
     &         BUFRECR(MAXRECORDS*NRHS_COL),       
     &         TOUCHED(NB_FS_IN_RHSCOMP),          
     &         stat=allocok)
      IF (allocok .GT. 0) THEN
        IF (LP .GT. 0) WRITE(LP, '(A)')
     &    'Error: Allocation problem in DMUMPS_SCATTER_DIST_RHS'
        INFO(1)=-13
        INFO(2)=NRHS_COL*MAXRECORDS*MAX_ACTIVE_SENDS+
     &          3*MAX_ACTIVE_SENDS+MAXRECORDS*(1+NRHS_COL)
     &          + NB_FS_IN_RHSCOMP
      ENDIF
      NB_BYTES_LOC=NB_BYTES_LOC +
     &  KEEP(34) * ( int(2*MAX_ACTIVE_SENDS,8) + int(MAXRECORDS,8) ) +
     &  KEEP(34) * (int(MAX_ACTIVE_SENDS,8) + int(NB_FS_IN_RHSCOMP,8)) +
     &  KEEP(35) * (
     &      int( MAXRECORDS,8)*int(NRHS_COL,8)*int(MAX_ACTIVE_SENDS,8)
     &      + int(MAXRECORDS,8) * int(NRHS_COL,8) )
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, allocok, 1,
     &                    MPI_INTEGER, MPI_SUM,
     &                    COMM_NODES, IERR_MPI )
      IF (allocok .NE. 0) RETURN
      NB_ACTIVE_SENDS = 0
      DO IREQ = 1, MAX_ACTIVE_SENDS
        IS_SEND_ACTIVE(IREQ) = .FALSE.
      ENDDO
      NB_FS_TOUCHED = 0
      DO IFS = 1, NB_FS_IN_RHSCOMP
        TOUCHED(IFS) = .FALSE.
      ENDDO
      IPROC_MAX=maxloc(NBROWSTOSEND,DIM=1)-1
      DO WHILE (NBROWSTOSEND(IPROC_MAX+1) .NE. 0)
        IF (IPROC_MAX .EQ. MYID_NODES) THEN
          CALL DMUMPS_DR_ASSEMBLE_LOCAL()
        ELSE
          CALL DMUMPS_DR_TRY_SEND(IPROC_MAX)
        ENDIF
        CALL DMUMPS_DR_TRY_RECV()
        CALL DMUMPS_DR_TRY_FREE_SEND()
        IPROC_MAX=maxloc(NBROWSTOSEND,DIM=1)-1
      ENDDO
      DO WHILE ( NBROWSTORECV .NE. 0)
        CALL DMUMPS_DR_TRY_RECV()
        CALL DMUMPS_DR_TRY_FREE_SEND()
      ENDDO
      DO WHILE (NB_ACTIVE_SENDS .NE. 0)
        CALL DMUMPS_DR_TRY_FREE_SEND()
      ENDDO
      CALL DMUMPS_DR_EMPTY_ROWS()
      RETURN
      CONTAINS
        SUBROUTINE DMUMPS_DR_BUILD_NBROWSTORECV()
        INTEGER :: IPROC
        DO IPROC = 0, NSLAVES-1
          CALL MPI_REDUCE( NBROWSTOSEND(IPROC+1), NBROWSTORECV,
     &                   1, MPI_INTEGER,
     &                   MPI_SUM, IPROC, COMM_NODES, IERR_MPI )
        ENDDO
        END SUBROUTINE DMUMPS_DR_BUILD_NBROWSTORECV
        SUBROUTINE DMUMPS_DR_TRY_RECV()
        IMPLICIT NONE
        INCLUDE 'mumps_tags.h'
        INTEGER :: MPI_STATUS(MPI_STATUS_SIZE), MSGSOU
        INTEGER :: NBRECORDS
        LOGICAL :: FLAG
        CALL MPI_IPROBE( MPI_ANY_SOURCE, DistRhsI, COMM_NODES,
     &                   FLAG, MPI_STATUS, IERR_MPI )
        IF (FLAG) THEN
          MSGSOU = MPI_STATUS( MPI_SOURCE )
          CALL MPI_GET_COUNT(MPI_STATUS, MPI_INTEGER,
     &                       NBRECORDS, IERR_MPI)
          CALL MPI_RECV(BUFRECI(1), NBRECORDS, MPI_INTEGER,
     &                  MSGSOU, DistRhsI,
     &                  COMM_NODES, MPI_STATUS, IERR_MPI)
          CALL MPI_RECV(BUFRECR(1), NBRECORDS*NRHS_COL,
     &                  MPI_DOUBLE_PRECISION,
     &                  MSGSOU, DistRhsR,
     &                  COMM_NODES, MPI_STATUS, IERR_MPI)
          CALL DMUMPS_DR_ASSEMBLE_FROM_BUFREC(NBRECORDS,
     &                                        BUFRECI, BUFRECR)
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_DR_TRY_RECV
        SUBROUTINE DMUMPS_DR_ASSEMBLE_FROM_BUFREC
     &             (NBRECORDS, BUFRECI_ARG, BUFRECR_ARG)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NBRECORDS
        INTEGER, INTENT(INOUT) :: BUFRECI_ARG(NBRECORDS) 
        DOUBLE PRECISION, INTENT(IN) :: BUFRECR_ARG(NBRECORDS,
     &                                     NRHS_COL)
        INTEGER :: I, K, IRHSCOMP, IFIRSTNOTTOUCHED, ILASTNOTTOUCHED
        IFIRSTNOTTOUCHED = NBRECORDS+1
        ILASTNOTTOUCHED  = 0
        DO I = 1, NBRECORDS
          IF (BUFRECI(I) .LE. 0) THEN
            WRITE(*,*) "Internal error 1 in DMUMPS_DR_TRY_RECV",
     &      I, BUFRECI(I), BUFRECI(1)
            CALL MUMPS_ABORT()
          ENDIF
          IRHSCOMP=POSINRHSCOMP_FWD(BUFRECI(I))
          BUFRECI_ARG(I)=IRHSCOMP
          IF ( .NOT. TOUCHED(IRHSCOMP) ) THEN
            IFIRSTNOTTOUCHED=min(IFIRSTNOTTOUCHED,I)
            ILASTNOTTOUCHED=max(ILASTNOTTOUCHED,I)
          ENDIF
        ENDDO
!$      OMP_FLAG = ( NRHS_COL.GE.KEEP(362) .AND.
!$   &               NRHS_COL*NBRECORDS .GE. KEEP(363)/2)
!$OMP PARALLEL DO PRIVATE(I,IRHSCOMP) IF (OMP_FLAG)
        DO K = 1, NRHS_COL
          DO I = IFIRSTNOTTOUCHED, ILASTNOTTOUCHED
            IRHSCOMP=BUFRECI_ARG(I)
            IF (.NOT. TOUCHED(IRHSCOMP)) THEN
              RHSCOMP(IRHSCOMP,K)=ZERO
            ENDIF
          ENDDO
          DO I = 1, NBRECORDS
            IRHSCOMP=BUFRECI_ARG(I)
            RHSCOMP(IRHSCOMP,K) = RHSCOMP(IRHSCOMP,K) +
     &                            BUFRECR_ARG(I,K)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        DO I = 1, NBRECORDS
          IRHSCOMP = BUFRECI_ARG(I)
          IF (.NOT. TOUCHED(IRHSCOMP)) THEN
            NB_FS_TOUCHED = NB_FS_TOUCHED + 1
            TOUCHED(IRHSCOMP) = .TRUE.
          ENDIF
        ENDDO
        NBROWSTORECV = NBROWSTORECV - NBRECORDS
        RETURN
        END SUBROUTINE DMUMPS_DR_ASSEMBLE_FROM_BUFREC
        SUBROUTINE DMUMPS_DR_ASSEMBLE_LOCAL()
        INTEGER :: NBRECORDS, I, K, IFIRSTNOTTOUCHED
        INTEGER :: Iloc       
        INTEGER :: Iglob      
        INTEGER :: IRHSCOMP   
        INTEGER(8) :: ISHIFT
        IF ( NBROWSTOSEND(MYID_NODES+1) .EQ. 0) THEN
          WRITE(*,*) "Internal error in DMUMPS_DR_ASSEMBLE_LOCAL"
          CALL MUMPS_ABORT()
        ENDIF
        NBRECORDS=min(MAXRECORDS, NBROWSTOSEND(MYID_NODES+1))
        IFIRSTNOTTOUCHED=NBRECORDS+1
        DO I = 1, NBRECORDS
          IRHSCOMP = POSINRHSCOMP_FWD(IRHS_loc(
     &               IRHS_loc_sorted(NEXTROWTOSEND(MYID_NODES+1)+I-1)))
          IF (.NOT. TOUCHED(IRHSCOMP)) THEN
            IFIRSTNOTTOUCHED=I
            EXIT
          ENDIF
        ENDDO
        IF (LSCAL) THEN
!$        OMP_FLAG = (NRHS_COL.GE.KEEP(362) .AND.
!$   &                NRHS_COL*NBRECORDS .GE. KEEP(363)/2)
!$OMP PARALLEL DO PRIVATE(K, ISHIFT, I, IRHSCOMP, Iloc, Iglob)
!$OMP&  FIRSTPRIVATE(NBRECORDS) IF (OMP_FLAG)
          DO K = 1, NRHS_COL
            ISHIFT = (K-1) * LRHS_loc
            DO I = IFIRSTNOTTOUCHED, NBRECORDS
              IRHSCOMP = POSINRHSCOMP_FWD(IRHS_loc(
     &               IRHS_loc_sorted(NEXTROWTOSEND(MYID_NODES+1)+I-1)))
               IF (.NOT. TOUCHED(IRHSCOMP)) THEN
                 RHSCOMP(IRHSCOMP,K)=ZERO
               ENDIF
            ENDDO
            DO I = 1, NBRECORDS
              Iloc = IRHS_loc_sorted(NEXTROWTOSEND(MYID_NODES+1)+I-1)
              Iglob = IRHS_loc(Iloc)
              IRHSCOMP = POSINRHSCOMP_FWD(Iglob)
              RHSCOMP(IRHSCOMP,K) = RHSCOMP(IRHSCOMP,K)+
     &                              RHS_loc(Iloc+ISHIFT)*
     &                        scaling_data_dr%SCALING_LOC(Iloc)
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ELSE
!$        OMP_FLAG = (NRHS_COL.GE.KEEP(362) .AND.
!$   &                NRHS_COL*NBRECORDS .GE. KEEP(363)/2)
!$OMP PARALLEL DO PRIVATE(K, ISHIFT, I, IRHSCOMP, Iloc, Iglob)
!$OMP&  FIRSTPRIVATE(NBRECORDS) IF (OMP_FLAG)
          DO K = 1, NRHS_COL
            ISHIFT = (K-1) * LRHS_loc
            DO I = IFIRSTNOTTOUCHED, NBRECORDS
              IRHSCOMP = POSINRHSCOMP_FWD(IRHS_loc(
     &               IRHS_loc_sorted(NEXTROWTOSEND(MYID_NODES+1)+I-1)))
               IF (.NOT. TOUCHED(IRHSCOMP)) THEN
                 RHSCOMP(IRHSCOMP,K)=ZERO
               ENDIF
            ENDDO
            DO I = 1, NBRECORDS
              Iloc = IRHS_loc_sorted(NEXTROWTOSEND(MYID_NODES+1)+I-1)
              Iglob = IRHS_loc(Iloc)
              IRHSCOMP = POSINRHSCOMP_FWD(Iglob)
              RHSCOMP(IRHSCOMP,K) = RHSCOMP(IRHSCOMP,K)+
     &                              RHS_loc(Iloc+ISHIFT)
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDIF
        DO I = 1, NBRECORDS
          IRHSCOMP = POSINRHSCOMP_FWD(IRHS_loc(
     &               IRHS_loc_sorted(NEXTROWTOSEND(MYID_NODES+1)+I-1)))
          IF (.NOT. TOUCHED(IRHSCOMP)) THEN
            NB_FS_TOUCHED = NB_FS_TOUCHED + 1
            TOUCHED(IRHSCOMP) = .TRUE.
          ENDIF
        ENDDO
        NEXTROWTOSEND(MYID_NODES+1)=NEXTROWTOSEND(MYID_NODES+1)+
     &                              NBRECORDS
        NBROWSTOSEND(MYID_NODES+1)=NBROWSTOSEND(MYID_NODES+1)-
     &                             NBRECORDS
        NBROWSTORECV = NBROWSTORECV - NBRECORDS
        RETURN
        END SUBROUTINE DMUMPS_DR_ASSEMBLE_LOCAL
        SUBROUTINE DMUMPS_DR_GET_NEW_BUF( IBUF )
        INTEGER, INTENT(OUT) :: IBUF
        INTEGER :: I
        IBUF = -1
        IF (NB_ACTIVE_SENDS .NE. MAX_ACTIVE_SENDS) THEN
          DO I=1, MAX_ACTIVE_SENDS
            IF (.NOT. IS_SEND_ACTIVE(I)) THEN
              IBUF = I
              EXIT
            ENDIF
          ENDDO
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_DR_GET_NEW_BUF
        SUBROUTINE DMUMPS_DR_TRY_FREE_SEND()
        INTEGER :: MPI_STATUS(MPI_STATUS_SIZE)
        INTEGER :: I
        LOGICAL :: FLAG
        IF (NB_ACTIVE_SENDS .GT. 0) THEN
          DO I=1, MAX_ACTIVE_SENDS
            IF (IS_SEND_ACTIVE(I)) THEN
              CALL MPI_TEST( MPI_REQR(I), FLAG, MPI_STATUS, IERR_MPI )
              IF (FLAG) THEN
                 CALL MPI_WAIT(MPI_REQI(I), MPI_STATUS, IERR_MPI)
                 NB_ACTIVE_SENDS = NB_ACTIVE_SENDS - 1
                 IS_SEND_ACTIVE(I)=.FALSE.
                 IF (NB_ACTIVE_SENDS .EQ. 0) THEN
                   RETURN
                 ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_DR_TRY_FREE_SEND
        SUBROUTINE DMUMPS_DR_TRY_SEND(IPROC_ARG)
        IMPLICIT NONE
        INTEGER, INTENT(IN)    :: IPROC_ARG
        INCLUDE 'mumps_tags.h'
        INTEGER :: NBRECORDS, IBUF, I, K
        INTEGER(8) :: IPOSRHS 
        INTEGER    :: IPOSBUF 
        IF (IPROC_ARG .EQ. MYID_NODES) THEN
          WRITE(*,*) "Internal error 1 in DMUMPS_DR_TRY_SEND"
          CALL MUMPS_ABORT()
        ENDIF
        IF (NBROWSTOSEND(IPROC_ARG+1) .EQ. 0) THEN
          WRITE(*,*) "Internal error 2 in DMUMPS_DR_TRY_SEND"
          CALL MUMPS_ABORT()
        ENDIF
        CALL DMUMPS_DR_GET_NEW_BUF(IBUF)
        IF (IBUF .GT. 0) THEN
          NBRECORDS = min(MAXRECORDS,NBROWSTOSEND(IPROC_ARG+1))
!$        OMP_FLAG = .FALSE.
!$        CHUNK = NRHS_COL*NBRECORDS
!$        IF (CHUNK .GE. KEEP(363)) THEN
!$          OMP_FLAG = .TRUE.
!$          CHUNK = max((CHUNK+NOMP-1)/NOMP,KEEP(363)/2)
!$        ENDIF
          IF (LSCAL) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK)
!$OMP& PRIVATE(K, I, IPOSBUF, IPOSRHS, Iloc) IF (OMP_FLAG)
            DO K=1, NRHS_COL
              DO I = 1, NBRECORDS
                IPOSBUF = (K-1)*NBRECORDS
                IPOSRHS = int(K-1,8)*int(LRHS_loc,8)
                Iloc = IRHS_loc_sorted(NEXTROWTOSEND(IPROC_ARG+1)+I-1)
                BUFR( IPOSBUF + I, IBUF )
     &                = RHS_loc( IPOSRHS + Iloc ) *
     &                  scaling_data_dr%SCALING_LOC(Iloc)
              ENDDO
            ENDDO
!$OMP END PARALLEL DO
          ELSE
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK)
!$OMP& PRIVATE(K, I, IPOSBUF, IPOSRHS, Iloc) IF (OMP_FLAG)
            DO K=1, NRHS_COL
              DO I = 1, NBRECORDS
                IPOSBUF = (K-1)*NBRECORDS
                IPOSRHS = int(K-1,8)*int(LRHS_loc,8)
                Iloc = IRHS_loc_sorted(NEXTROWTOSEND(IPROC_ARG+1)+I-1)
                BUFR( IPOSBUF + I, IBUF )
     &              = RHS_loc( IPOSRHS + Iloc )
              ENDDO
            ENDDO
!$OMP END PARALLEL DO
          ENDIF
          DO I = 1, NBRECORDS
            Iloc = IRHS_loc_sorted(NEXTROWTOSEND(IPROC_ARG+1)+I-1)
            IRHS_loc_sorted(NEXTROWTOSEND(IPROC_ARG+1)+I-1)
     &           = IRHS_loc(Iloc)
          ENDDO
          CALL MPI_ISEND( IRHS_loc_sorted(NEXTROWTOSEND(IPROC_ARG+1)),
     &                    NBRECORDS, MPI_INTEGER, IPROC_ARG, DistRhsI,
     &                    COMM_NODES, MPI_REQI(IBUF), IERR_MPI )
          CALL MPI_ISEND( BUFR(1,IBUF), NBRECORDS*NRHS_COL,
     &                    MPI_DOUBLE_PRECISION,
     &                    IPROC_ARG, DistRhsR,
     &                    COMM_NODES, MPI_REQR(IBUF), IERR_MPI )
          NEXTROWTOSEND(IPROC_ARG+1)=NEXTROWTOSEND(IPROC_ARG+1)+
     &                               NBRECORDS
          NBROWSTOSEND(IPROC_ARG+1)=NBROWSTOSEND(IPROC_ARG+1)-NBRECORDS
          NB_ACTIVE_SENDS = NB_ACTIVE_SENDS + 1
          IS_SEND_ACTIVE(IBUF)=.TRUE.
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_DR_TRY_SEND
        SUBROUTINE DMUMPS_DR_EMPTY_ROWS()
        INTEGER :: K, IFS
        IF ( NB_FS_TOUCHED .NE. NB_FS_IN_RHSCOMP ) THEN
!$        OMP_FLAG = (NRHS_COL .GE. KEEP(362)) .AND.
!$   &    (NRHS_COL*NB_FS_IN_RHSCOMP >  KEEP(363)/2)
!$OMP PARALLEL DO FIRSTPRIVATE(NB_FS_IN_RHSCOMP) IF (OMP_FLAG)
          DO K = 1, NRHS_COL
            DO IFS = 1, NB_FS_IN_RHSCOMP
              IF ( .NOT. TOUCHED(IFS) ) THEN
                RHSCOMP( IFS, K) = ZERO
              ENDIF
            ENDDO
            DO IFS = NB_FS_IN_RHSCOMP +1, LD_RHSCOMP
              RHSCOMP (IFS, K) = ZERO
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ELSE
!$        OMP_FLAG = .FALSE.
!$        CHUNK8 = int(NRHS_COL,8)*int(LD_RHSCOMP-NB_FS_IN_RHSCOMP,8)
!$        CHUNK8 = max(CHUNK8,1_8)
!$        IF (CHUNK8 .GE. int(KEEP(363),8)) THEN
!$          OMP_FLAG = .TRUE.
!$          CHUNK8 = max((CHUNK8+NOMP-1)/NOMP,int(KEEP(363)/2,8))
!$        ENDIF
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK8)
!$OMP&  IF (OMP_FLAG)
          DO K = 1, NRHS_COL
            DO IFS = NB_FS_IN_RHSCOMP +1, LD_RHSCOMP
             RHSCOMP (IFS, K) = ZERO
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_DR_EMPTY_ROWS
      END SUBROUTINE DMUMPS_SCATTER_DIST_RHS
      SUBROUTINE DMUMPS_SOL_INIT_IRHS_loc(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
        TYPE (DMUMPS_STRUC) :: id
        INTEGER, PARAMETER :: MASTER = 0
        INTEGER            :: ROW_OR_COL_INDICES 
        INTEGER            :: IERR_MPI
        LOGICAL            :: I_AM_SLAVE
        INTEGER, POINTER   :: idIRHS_loc(:)
        INTEGER, POINTER   :: UNS_PERM(:)
        INTEGER            :: UNS_PERM_TO_BE_DONE, I, allocok
        INTEGER, TARGET    :: IDUMMY(1)
        INCLUDE 'mpif.h'
        NULLIFY(UNS_PERM)
        IF (id%JOB .NE. 9) THEN
          WRITE(*,*) "Internal error 1 in DMUMPS_SOL_INIT_IRHS_loc"
          CALL MUMPS_ABORT()
        ENDIF
        I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &               ( id%MYID .eq. MASTER .AND.
     &                 id%KEEP(46) .eq. 1 ) )
        IF (id%MYID .EQ. MASTER) THEN
          IF (id%ICNTL(20).EQ.10) THEN
            ROW_OR_COL_INDICES = 0 
          ELSE IF (id%ICNTL(20).EQ.11) THEN
            ROW_OR_COL_INDICES = 1 
          ELSE 
            ROW_OR_COL_INDICES = 0 
          ENDIF
          IF (id%ICNTL(9) .NE. 1) THEN
            ROW_OR_COL_INDICES = 1 - ROW_OR_COL_INDICES
          ENDIF
          IF (id%KEEP(23).NE.0 .AND. id%ICNTL(9) .NE.1) THEN
            UNS_PERM_TO_BE_DONE = 1
          ELSE
            UNS_PERM_TO_BE_DONE = 0
          ENDIF
        ENDIF
        CALL MPI_BCAST(ROW_OR_COL_INDICES,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR_MPI)
        CALL MPI_BCAST(UNS_PERM_TO_BE_DONE,1,MPI_INTEGER,MASTER,
     &               id%COMM,IERR_MPI)
        IF ( I_AM_SLAVE ) THEN
          IF (id%KEEP(89) .GT. 0) THEN
            IF (.NOT. associated(id%IRHS_loc)) THEN
              id%INFO(1)=-22
              id%INFO(2)=17
            ELSE IF (size(id%IRHS_loc) < id%KEEP(89) ) THEN
              id%INFO(1)=-22
              id%INFO(2)=17
            ENDIF
          ENDIF
        ENDIF
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                       id%INFO(1),
     &                       id%COMM, id%MYID )
        IF (id%INFO(1).LT.0) goto 500
        IF (I_AM_SLAVE) THEN
          IF (associated(id%IRHS_loc)) THEN
            IF (size(id%IRHS_loc) .GT. 0) THEN
              idIRHS_loc => id%IRHS_loc
            ELSE
              idIRHS_loc => IDUMMY
            ENDIF
          ELSE
            idIRHS_loc => IDUMMY
          ENDIF
          CALL MUMPS_BUILD_IRHS_loc(id%MYID_NODES, id%NSLAVES, id%N,
     &    id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1), id%IS(1),
     &    max(1, id%KEEP(32)),
     &    id%STEP(1), id%PROCNODE_STEPS(1), idIRHS_loc(1),
     &    ROW_OR_COL_INDICES)
        ENDIF
        IF (UNS_PERM_TO_BE_DONE .EQ. 1) THEN
          IF (id%MYID.NE.MASTER) THEN
            ALLOCATE(UNS_PERM(id%N),stat=allocok)
            IF (allocok > 0) THEN
                id%INFO(1)=-13
                id%INFO(2)=id%N
                GOTO 100
            ENDIF
          ENDIF
 100      CONTINUE
          CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                         id%INFO(1),
     &                         id%COMM, id%MYID )
          IF (id%INFO(1) .LT. 0) GOTO 500
          IF ( id%MYID .EQ. MASTER ) THEN
            UNS_PERM => id%UNS_PERM
          ENDIF
          CALL MPI_BCAST(UNS_PERM(1),id%N,MPI_INTEGER,MASTER,
     &               id%COMM,IERR_MPI)
          IF (I_AM_SLAVE .AND. id%KEEP(89) .NE.0) THEN
            DO I=1, id%KEEP(89)
              id%IRHS_loc(I)=UNS_PERM(id%IRHS_loc(I))
            ENDDO
          ENDIF
        ENDIF
 500    CONTINUE
        IF (id%MYID.NE.MASTER) THEN
           IF (associated(UNS_PERM)) DEALLOCATE(UNS_PERM)
        ENDIF
        NULLIFY(UNS_PERM)
        RETURN
      END SUBROUTINE DMUMPS_SOL_INIT_IRHS_loc
