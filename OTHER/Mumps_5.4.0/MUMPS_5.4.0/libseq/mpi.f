C
C  This file is part of MUMPS 5.4.0, released
C  on Tue Apr 13 15:26:30 UTC 2021
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
C*******************************************************************
C
C  This file contains stub MPI/BLACS/ScaLAPACK library functions for
C  linking/running MUMPS on a platform where MPI is not installed.
C
C*******************************************************************
C
C MPI
C
C******************************************************************
      SUBROUTINE MPI_BSEND( BUF, CNT, DATATYPE, DEST, TAG, COMM,
     &            IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, DATATYPE, DEST, TAG, COMM, IERR
      INTEGER BUF(*)
      WRITE(*,*) 'Error. MPI_BSEND should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_BSEND
C***********************************************************************
      SUBROUTINE MPI_BUFFER_ATTACH(BUF, CNT,  IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, IERR
      INTEGER BUF(*)
      IERR = 0
      RETURN
      END SUBROUTINE MPI_BUFFER_ATTACH
C***********************************************************************
      SUBROUTINE MPI_BUFFER_DETACH(BUF, CNT,  IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, IERR
      INTEGER BUF(*)
           IERR = 0
      RETURN
      END SUBROUTINE MPI_BUFFER_DETACH
      SUBROUTINE MPI_GATHER( SENDBUF, CNT, 
     &         DATATYPE, RECVBUF, RECCNT, RECTYPE,
     &         ROOT, COMM, IERR )
      IMPLICIT NONE
      INTEGER CNT, DATATYPE, RECCNT, RECTYPE, ROOT, COMM, IERR
      INTEGER SENDBUF(*), RECVBUF(*)
      IF ( RECCNT .NE. CNT ) THEN
        WRITE(*,*) 'ERROR in MPI_GATHER, RECCNT != CNT'
        STOP
      ELSE
        CALL MUMPS_COPY( CNT, SENDBUF, RECVBUF, DATATYPE, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) 'ERROR in MPI_GATHER, DATATYPE=',DATATYPE
          STOP
        END IF
      END IF
      IERR = 0
      RETURN
      END SUBROUTINE MPI_GATHER
C***********************************************************************
      SUBROUTINE MPI_GATHERV( SENDBUF, CNT, 
     &         DATATYPE, RECVBUF, RECCNT, DISPLS, RECTYPE,
     &         ROOT, COMM, IERR )
      IMPLICIT NONE
      INTEGER CNT, DATATYPE, RECTYPE, ROOT, COMM, IERR
      INTEGER RECCNT(1)
      INTEGER SENDBUF(*), RECVBUF(*)
      INTEGER DISPLS(*)
C
C     Note that DISPLS is ignored in this version. One may
C     want to copy in reception buffer with a shift DISPLS(1).
C     This requires passing the offset DISPLS(1) to
C     "MUMPS_COPY_DATATYPE" routines.
C
      IF ( RECCNT(1) .NE. CNT ) THEN
        WRITE(*,*) 'ERROR in MPI_GATHERV, RECCNT(1) != CNT'
        STOP
      ELSE
        CALL MUMPS_COPY( CNT, SENDBUF, RECVBUF, DATATYPE, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) 'ERROR in MPI_GATHERV, DATATYPE=',DATATYPE
          STOP
        END IF
      END IF
      IERR = 0
      RETURN
      END SUBROUTINE MPI_GATHERV
C***********************************************************************
      SUBROUTINE MPI_ALLREDUCE( SENDBUF, RECVBUF, CNT, DATATYPE,
     &                          OPERATION, COMM, IERR )
      IMPLICIT NONE
      INTEGER CNT, DATATYPE, OPERATION, COMM, IERR
      INTEGER SENDBUF(*), RECVBUF(*)
      LOGICAL, EXTERNAL :: MUMPS_IS_IN_PLACE
      IF (.NOT. MUMPS_IS_IN_PLACE(SENDBUF, CNT)) THEN
        CALL MUMPS_COPY( CNT, SENDBUF, RECVBUF, DATATYPE, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) 'ERROR in MPI_ALLREDUCE, DATATYPE=',DATATYPE
          STOP
        END IF
      ENDIF
      IERR = 0
      RETURN
      END SUBROUTINE MPI_ALLREDUCE
C***********************************************************************
      SUBROUTINE MPI_REDUCE( SENDBUF, RECVBUF, CNT, DATATYPE, OP,
     &           ROOT, COMM, IERR )
      IMPLICIT NONE
      INTEGER CNT, DATATYPE, OP, ROOT, COMM, IERR
      INTEGER SENDBUF(*), RECVBUF(*)
      LOGICAL, EXTERNAL :: MUMPS_IS_IN_PLACE
      IF (.NOT. MUMPS_IS_IN_PLACE(SENDBUF, CNT)) THEN
        CALL MUMPS_COPY( CNT, SENDBUF, RECVBUF, DATATYPE, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) 'ERROR in MPI_REDUCE, DATATYPE=',DATATYPE
          STOP
        END IF
      ENDIF
      IERR = 0
      RETURN
      END SUBROUTINE MPI_REDUCE
C***********************************************************************
      SUBROUTINE MPI_REDUCE_SCATTER( SENDBUF, RECVBUF, RCVCNT, 
     &           DATATYPE, OP, COMM, IERR )
      IMPLICIT NONE
      INTEGER RCVCNT, DATATYPE, OP, COMM, IERR
      INTEGER SENDBUF(*), RECVBUF(*)
      LOGICAL, EXTERNAL :: MUMPS_IS_IN_PLACE
      IF (.NOT. MUMPS_IS_IN_PLACE(SENDBUF, RCVCNT)) THEN
        CALL MUMPS_COPY( RCVCNT, SENDBUF, RECVBUF, DATATYPE, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) 'ERROR in MPI_REDUCE_SCATTER, DATATYPE=',DATATYPE
          STOP
        END IF
      ENDIF
      IERR = 0
      RETURN
      END SUBROUTINE MPI_REDUCE_SCATTER
C***********************************************************************
      SUBROUTINE MPI_ABORT( COMM, IERRCODE, IERR )
      IMPLICIT NONE
      INTEGER COMM, IERRCODE, IERR
      WRITE(*,*) "** MPI_ABORT called"
      STOP
      END SUBROUTINE MPI_ABORT
C***********************************************************************
      SUBROUTINE MPI_ALLTOALL( SENDBUF, SENDCNT, SENDTYPE,
     &                         RECVBUF, RECVCNT, RECVTYPE, COMM, IERR )
      IMPLICIT NONE
      INTEGER SENDCNT, SENDTYPE, RECVCNT, RECVTYPE, COMM, IERR
      INTEGER SENDBUF(*), RECVBUF(*)
      IF ( RECVCNT .NE. SENDCNT ) THEN
        WRITE(*,*) 'ERROR in MPI_ALLTOALL, RECVCNT != SENDCNT'
        STOP
      ELSE IF ( RECVTYPE .NE. SENDTYPE ) THEN
        WRITE(*,*) 'ERROR in MPI_ALLTOALL, RECVTYPE != SENDTYPE'
        STOP
      ELSE
        CALL MUMPS_COPY( SENDCNT, SENDBUF, RECVBUF, SENDTYPE, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) 'ERROR in MPI_ALLTOALL, SENDTYPE=',SENDTYPE
          STOP
        END IF
      END IF
      IERR = 0
      RETURN
      END SUBROUTINE MPI_ALLTOALL
C***********************************************************************
      SUBROUTINE MPI_ATTR_PUT( COMM, KEY, VAL, IERR )
      IMPLICIT NONE
      INTEGER COMM, KEY, VAL, IERR
      RETURN
      END SUBROUTINE MPI_ATTR_PUT
C***********************************************************************
      SUBROUTINE MPI_BARRIER( COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_BARRIER
C***********************************************************************
      SUBROUTINE MPI_GET_PROCESSOR_NAME( NAME, RESULTLEN, IERROR)
      CHARACTER (LEN=*) NAME
      INTEGER RESULTLEN,IERROR
      RESULTLEN = 1
      IERROR = 0
      NAME = 'X'
      RETURN
      END SUBROUTINE MPI_GET_PROCESSOR_NAME
C***********************************************************************
      SUBROUTINE MPI_BCAST( BUFFER, CNT, DATATYPE, ROOT, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, DATATYPE, ROOT, COMM, IERR
      INTEGER BUFFER( * )
      IERR = 0
      RETURN
      END SUBROUTINE MPI_BCAST
C***********************************************************************
      SUBROUTINE MPI_CANCEL( IREQ, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IREQ, IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_CANCEL
C***********************************************************************
      SUBROUTINE MPI_COMM_CREATE( COMM, GROUP, COMM2, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, GROUP, COMM2, IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_COMM_CREATE
C***********************************************************************
      SUBROUTINE MPI_COMM_DUP( COMM, COMM2, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, COMM2, IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_COMM_DUP
C***********************************************************************
      SUBROUTINE MPI_COMM_FREE( COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_COMM_FREE
C***********************************************************************
      SUBROUTINE MPI_COMM_GROUP( COMM, GROUP, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, GROUP, IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_COMM_GROUP
C***********************************************************************
      SUBROUTINE MPI_COMM_RANK( COMM, RANK, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, RANK, IERR
      RANK = 0
      IERR = 0
      RETURN
      END SUBROUTINE MPI_COMM_RANK
C***********************************************************************
      SUBROUTINE MPI_COMM_SIZE( COMM, SIZE, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, SIZE, IERR
      SIZE = 1
      IERR = 0
      RETURN
      END SUBROUTINE MPI_COMM_SIZE
C***********************************************************************
      SUBROUTINE MPI_COMM_SPLIT( COMM, COLOR, KEY, COMM2, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, COLOR, KEY, COMM2, IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_COMM_SPLIT
C***********************************************************************
c     SUBROUTINE MPI_ERRHANDLER_SET( COMM, ERRHANDLER, IERR )
c     IMPLICIT NONE
c     INCLUDE 'mpif.h'
c     INTEGER COMM, ERRHANDLER, IERR
c     IERR = 0
c     RETURN
c     END SUBROUTINE MPI_ERRHANDLER_SET
C***********************************************************************
      SUBROUTINE MPI_FINALIZE( IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_FINALIZE
C***********************************************************************
      SUBROUTINE MPI_GET_COUNT( STATUS, DATATYPE, CNT, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER DATATYPE, CNT, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      WRITE(*,*) 'Error. MPI_GET_CNT should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_GET_COUNT
C***********************************************************************
      SUBROUTINE MPI_GROUP_FREE( GROUP, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER GROUP, IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_GROUP_FREE
C***********************************************************************
      SUBROUTINE MPI_GROUP_RANGE_EXCL( GROUP, N, RANGES, GROUP2, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER GROUP, N, GROUP2, IERR
      INTEGER RANGES(*)
      IERR = 0
      RETURN
      END SUBROUTINE MPI_GROUP_RANGE_EXCL
C***********************************************************************
      SUBROUTINE MPI_GROUP_SIZE( GROUP, SIZE, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER GROUP, SIZE, IERR
      SIZE = 1 ! Or should it be zero ?
      IERR = 0
      RETURN
      END SUBROUTINE MPI_GROUP_SIZE
C***********************************************************************
      SUBROUTINE MPI_INIT(IERR)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_INIT
C***********************************************************************
      SUBROUTINE MPI_INITIALIZED( FLAG, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      LOGICAL FLAG
      INTEGER IERR
      FLAG = .TRUE.
      IERR = 0
      RETURN
      END SUBROUTINE MPI_INITIALIZED
C***********************************************************************
      SUBROUTINE MPI_IPROBE( SOURCE, TAG, COMM, FLAG, STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER SOURCE, TAG, COMM, IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      LOGICAL FLAG
      FLAG = .FALSE.
      IERR = 0
      RETURN
      END SUBROUTINE MPI_IPROBE
C***********************************************************************
      SUBROUTINE MPI_IRECV( BUF, CNT, DATATYPE, SOURCE, TAG, COMM,
     &           IREQ, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, DATATYPE, SOURCE, TAG, COMM, IREQ, IERR
      INTEGER BUF(*)
      IERR = 0
      RETURN
      END SUBROUTINE MPI_IRECV
C***********************************************************************
      SUBROUTINE MPI_ISEND( BUF, CNT, DATATYPE, DEST, TAG, COMM,
     &           IREQ, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, DATATYPE, DEST, TAG, COMM, IERR, IREQ
      INTEGER BUF(*)
      WRITE(*,*) 'Error. MPI_ISEND should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_ISEND
C***********************************************************************
      SUBROUTINE MPI_TYPE_COMMIT( NEWTYP, IERR_MPI )
      IMPLICIT NONE
      INTEGER NEWTYP, IERR_MPI
      RETURN
      END SUBROUTINE MPI_TYPE_COMMIT
C***********************************************************************
      SUBROUTINE MPI_TYPE_FREE( NEWTYP, IERR_MPI )
      IMPLICIT NONE
      INTEGER NEWTYP, IERR_MPI
      RETURN
      END SUBROUTINE MPI_TYPE_FREE
C***********************************************************************
      SUBROUTINE MPI_TYPE_CONTIGUOUS( LENGTH, DATATYPE, NEWTYPE,
     &                                IERR_MPI )
      IMPLICIT NONE
      INTEGER LENGTH, DATATYPE, NEWTYPE, IERR_MPI
      RETURN
      END SUBROUTINE MPI_TYPE_CONTIGUOUS
C***********************************************************************
      SUBROUTINE MPI_OP_CREATE( FUNC, COMMUTE, OP, IERR )
      IMPLICIT NONE
      EXTERNAL FUNC
      LOGICAL COMMUTE
      INTEGER OP, IERR
      OP = 0
      RETURN
      END SUBROUTINE MPI_OP_CREATE
C***********************************************************************
      SUBROUTINE MPI_OP_FREE( OP, IERR )
      IMPLICIT NONE
      INTEGER OP, IERR
      RETURN
      END SUBROUTINE MPI_OP_FREE
C***********************************************************************
      SUBROUTINE MPI_PACK( INBUF, INCNT, DATATYPE, OUTBUF, OUTCNT,
     &           POSITION, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER INCNT, DATATYPE, OUTCNT, POSITION, COMM, IERR
      INTEGER INBUF(*), OUTBUF(*)
      WRITE(*,*) 'Error. MPI_PACKED should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_PACK
C***********************************************************************
      SUBROUTINE MPI_PACK_SIZE( INCNT, DATATYPE, COMM, SIZE, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER INCNT, DATATYPE, COMM, SIZE, IERR
      WRITE(*,*) 'Error. MPI_PACK_SIZE should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_PACK_SIZE
C***********************************************************************
      SUBROUTINE MPI_PROBE( SOURCE, TAG, COMM, STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER SOURCE, TAG, COMM, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      WRITE(*,*) 'Error. MPI_PROBE should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_PROBE
C***********************************************************************
      SUBROUTINE MPI_RECV( BUF, CNT, DATATYPE, SOURCE, TAG, COMM,
     &           STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, DATATYPE, SOURCE, TAG, COMM, IERR
      INTEGER BUF(*), STATUS(MPI_STATUS_SIZE)
      WRITE(*,*) 'Error. MPI_RECV should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_RECV
C***********************************************************************
      SUBROUTINE MPI_REQUEST_FREE( IREQ, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IREQ, IERR
      IERR = 0
      RETURN
      END SUBROUTINE MPI_REQUEST_FREE
C***********************************************************************
      SUBROUTINE MPI_SEND( BUF, CNT, DATATYPE, DEST, TAG, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, DATATYPE, DEST, TAG, COMM, IERR
      INTEGER BUF(*)
      WRITE(*,*) 'Error. MPI_SEND should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_SEND
C***********************************************************************
      SUBROUTINE MPI_SSEND( BUF, CNT, DATATYPE, DEST, TAG, COMM, IERR)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, DATATYPE, DEST, TAG, COMM, IERR
      INTEGER BUF(*)
      WRITE(*,*) 'Error. MPI_SSEND should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_SSEND
C***********************************************************************
      SUBROUTINE MPI_TEST( IREQ, FLAG, STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IREQ, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      LOGICAL FLAG
      FLAG = .FALSE.
      IERR = 0
      RETURN
      END SUBROUTINE MPI_TEST
C***********************************************************************
      SUBROUTINE MPI_UNPACK( INBUF, INSIZE, POSITION, OUTBUF, OUTCNT,
     &           DATATYPE, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER INSIZE, POSITION, OUTCNT, DATATYPE, COMM, IERR
      INTEGER INBUF(*), OUTBUF(*)
      WRITE(*,*) 'Error. MPI_UNPACK should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_UNPACK
C***********************************************************************
      SUBROUTINE MPI_WAIT( IREQ, STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IREQ, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      WRITE(*,*) 'Error. MPI_WAIT should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_WAIT
C***********************************************************************
      SUBROUTINE MPI_WAITALL( CNT, ARRAY_OF_REQUESTS, STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER ARRAY_OF_REQUESTS( CNT )
      WRITE(*,*) 'Error. MPI_WAITALL should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_WAITALL
C***********************************************************************
      SUBROUTINE MPI_WAITANY( CNT, ARRAY_OF_REQUESTS, INDEX, STATUS,
     &           IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, INDEX, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER ARRAY_OF_REQUESTS( CNT )
      WRITE(*,*) 'Error. MPI_WAITANY should not be called.'
      STOP
      IERR = 0
      RETURN
      END SUBROUTINE MPI_WAITANY
C***********************************************************************
      DOUBLE PRECISION FUNCTION MPI_WTIME( )
C     elapsed time
      DOUBLE PRECISION VAL
C     write(*,*) 'Entering MPI_WTIME'
      CALL MUMPS_ELAPSE( VAL )
      MPI_WTIME = VAL
C     write(*,*) 'Exiting MPI_WTIME'
      RETURN
      END FUNCTION MPI_WTIME


C***********************************************************************
C
C  Utilities to copy data
C
C***********************************************************************

      SUBROUTINE MUMPS_COPY( CNT, SENDBUF, RECVBUF, DATATYPE, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER CNT, DATATYPE, IERR
      INTEGER SENDBUF(*), RECVBUF(*)
      IF ( DATATYPE .EQ. MPI_INTEGER ) THEN
         CALL MUMPS_COPY_INTEGER( SENDBUF, RECVBUF, CNT )
      ELSEIF ( DATATYPE .EQ. MPI_LOGICAL ) THEN
         CALL MUMPS_COPY_LOGICAL( SENDBUF, RECVBUF, CNT )
      ELSE IF ( DATATYPE .EQ. MPI_REAL ) THEN
         CALL MUMPS_COPY_REAL( SENDBUF, RECVBUF, CNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_PRECISION .OR.
     &        DATATYPE .EQ. MPI_REAL8 ) THEN
         CALL MUMPS_COPY_DOUBLE_PRECISION( SENDBUF, RECVBUF, CNT )
      ELSE IF ( DATATYPE .EQ. MPI_COMPLEX ) THEN
         CALL MUMPS_COPY_COMPLEX( SENDBUF, RECVBUF, CNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_COMPLEX ) THEN
         CALL MUMPS_COPY_DOUBLE_COMPLEX( SENDBUF, RECVBUF, CNT )
      ELSE IF ( DATATYPE .EQ. MPI_2DOUBLE_PRECISION) THEN
         CALL MUMPS_COPY_2DOUBLE_PRECISION( SENDBUF, RECVBUF, CNT )
      ELSE IF ( DATATYPE .EQ. MPI_2INTEGER) THEN
         CALL MUMPS_COPY_2INTEGER( SENDBUF, RECVBUF, CNT )
      ELSE IF ( DATATYPE .EQ. MPI_INTEGER8) THEN
        CALL MUMPS_COPY_INTEGER8( SENDBUF, RECVBUF, CNT )
      ELSE
        IERR=1
        RETURN
      END IF
      IERR=0
      RETURN
      END SUBROUTINE MUMPS_COPY

      SUBROUTINE MUMPS_COPY_INTEGER( S, R, N )
      IMPLICIT NONE
      INTEGER N
      INTEGER S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_INTEGER
      SUBROUTINE MUMPS_COPY_INTEGER8( S, R, N )
      IMPLICIT NONE
      INTEGER N
      INTEGER(8) S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_INTEGER8
      SUBROUTINE MUMPS_COPY_LOGICAL( S, R, N )
      IMPLICIT NONE
      INTEGER N
      LOGICAL S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END
      SUBROUTINE MUMPS_COPY_2INTEGER( S, R, N )
      IMPLICIT NONE
      INTEGER N
      INTEGER S(N+N),R(N+N)
      INTEGER I
      DO I = 1, N+N
        R(I) = S(I)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_2INTEGER
      SUBROUTINE MUMPS_COPY_REAL( S, R, N )
      IMPLICIT NONE
      INTEGER N
      REAL S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END
      SUBROUTINE MUMPS_COPY_2DOUBLE_PRECISION( S, R, N )
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION S(N+N),R(N+N)
      INTEGER I
      DO I = 1, N+N
        R(I) = S(I)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_2DOUBLE_PRECISION
      SUBROUTINE MUMPS_COPY_DOUBLE_PRECISION( S, R, N )
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END
      SUBROUTINE MUMPS_COPY_COMPLEX( S, R, N )
      IMPLICIT NONE
      INTEGER N
      COMPLEX S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_COMPLEX
      SUBROUTINE MUMPS_COPY_DOUBLE_COMPLEX( S, R, N )
      IMPLICIT NONE
      INTEGER N
C     DOUBLE COMPLEX S(N),R(N)
      COMPLEX(kind=kind(0.0D0)) :: S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END
      LOGICAL FUNCTION MUMPS_IS_IN_PLACE( SENDBUF, CNT )
      INTEGER SENDBUF(*), CNT
      INCLUDE 'mpif.h'
      INTEGER :: I
C     Check address using C code
      MUMPS_IS_IN_PLACE = .FALSE.
      IF ( CNT .GT. 0 ) THEN
        CALL MUMPS_CHECKADDREQUAL(SENDBUF(1), MPI_IN_PLACE, I)
        IF (I .EQ. 1) THEN
          MUMPS_IS_IN_PLACE = .TRUE.
        ENDIF
      ENDIF
C Begin old code which requires the MPI_IN_PLACE
C variable to have the F2003 attribute VOLATILE
C     IF ( CNT .GT. 0 ) THEN
C       MPI_IN_PLACE = -1
C       IF (SENDBUF(1) .EQ. MPI_IN_PLACE) THEN
C         MPI_IN_PLACE = -9876543
C         IF (MUMPS_CHECK_EQUAL(SENDBUF(1), MPI_IN_PLACE)) THEN
C           MUMPS_IS_IN_PLACE = .TRUE.
C         ENDIF
C       ENDIF
C     ENDIF
C End old code
      RETURN
      END FUNCTION MUMPS_IS_IN_PLACE
C Begin old code
C     LOGICAL FUNCTION MUMPS_CHECK_EQUAL(I,J)
C     INTEGER :: I,J
C     IF (I.EQ.J) THEN
C       MUMPS_CHECK_EQUAL = .TRUE. 
C     ELSE
C       MUMPS_CHECK_EQUAL = .FALSE. 
C     ENDIF
C     END FUNCTION MUMPS_CHECK_EQUAL
C End old code



C***********************************************************************
C
C     BLACS
C
C***********************************************************************
      SUBROUTINE blacs_gridinit( CNTXT, C, NPROW, NPCOL )
      IMPLICIT NONE
      INTEGER CNTXT, NPROW, NPCOL
      CHARACTER C
        WRITE(*,*) 'Error. BLACS_GRIDINIT should not be called.'
        STOP
      RETURN
      END SUBROUTINE blacs_gridinit
C***********************************************************************
      SUBROUTINE blacs_gridinfo( CNTXT, NPROW, NPCOL, MYROW, MYCOL )
      IMPLICIT NONE
      INTEGER CNTXT, NPROW, NPCOL, MYROW, MYCOL
        WRITE(*,*) 'Error. BLACS_GRIDINFO should not be called.'
        STOP
      RETURN
      END SUBROUTINE blacs_gridinfo
C***********************************************************************
      SUBROUTINE blacs_gridexit( CNTXT )
      IMPLICIT NONE
      INTEGER CNTXT
        WRITE(*,*) 'Error. BLACS_GRIDEXIT should not be called.'
        STOP
      RETURN
      END SUBROUTINE blacs_gridexit


C***********************************************************************
C
C     ScaLAPACK
C
C***********************************************************************
      SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC,
     &           ICTXT, LLD, INFO )
      IMPLICIT NONE
      INTEGER ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER DESC( * )
        WRITE(*,*) 'Error. DESCINIT should not be called.'
        STOP
      RETURN
      END SUBROUTINE DESCINIT
C***********************************************************************
      INTEGER FUNCTION numroc( N, NB, IPROC, ISRCPROC, NPROCS ) 
      INTEGER N, NB, IPROC, ISRCPROC, NPROCS
C     Can be called
      IF ( NPROCS .ne. 1 ) THEN
        WRITE(*,*) 'Error. Last parameter from NUMROC should be 1'
        STOP
      ENDIF
      IF ( IPROC .ne. 0 ) THEN
        WRITE(*,*) 'Error. IPROC should be 0 in NUMROC.'
        STOP
      ENDIF
      NUMROC = N
      RETURN
      END FUNCTION numroc
C***********************************************************************
      SUBROUTINE pcpotrf( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      COMPLEX            A( * )
        WRITE(*,*) 'Error. PCPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcpotrf
C***********************************************************************
      SUBROUTINE pcgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX            A( * )
        WRITE(*,*) 'Error. PCGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcgetrf
C***********************************************************************
      SUBROUTINE pctrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * )
        WRITE(*,*) 'Error. PCTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pctrtrs
C***********************************************************************
      SUBROUTINE pzpotrf( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
C     DOUBLE COMPLEX     A( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * )
        WRITE(*,*) 'Error. PZPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzpotrf
C***********************************************************************
      SUBROUTINE pzgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
C     DOUBLE COMPLEX     A( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * )
        WRITE(*,*) 'Error. PZGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzgetrf
C***********************************************************************
      SUBROUTINE pztrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
C     DOUBLE COMPLEX     A( * ), B( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
        WRITE(*,*) 'Error. PZTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pztrtrs
C***********************************************************************
      SUBROUTINE pspotrf( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      REAL               A( * )
        WRITE(*,*) 'Error. PSPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pspotrf
C***********************************************************************
      SUBROUTINE psgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      REAL               A( * )
        WRITE(*,*) 'Error. PSGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE psgetrf
C***********************************************************************
      SUBROUTINE pstrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
        WRITE(*,*) 'Error. PSTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pstrtrs
C***********************************************************************
      SUBROUTINE pdpotrf( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
        WRITE(*,*) 'Error. PDPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdpotrf
C***********************************************************************
      SUBROUTINE pdgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      DOUBLE PRECISION   A( * )
        WRITE(*,*) 'Error. PDGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdgetrf
C***********************************************************************
      SUBROUTINE pdtrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * )
        WRITE(*,*) 'Error. PDTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdtrtrs
C***********************************************************************
      SUBROUTINE INFOG2L( GRINDX, GCINDX, DESC, NPROW, NPCOL, MYROW,
     &                    MYCOL, LRINDX, LCINDX, RSRC, CSRC )
      IMPLICIT NONE
      INTEGER            CSRC, GCINDX, GRINDX, LRINDX, LCINDX, MYCOL,
     &                   MYROW, NPCOL, NPROW, RSRC
      INTEGER            DESC( * )
        WRITE(*,*) 'Error. INFOG2L should not be called.'
        STOP
      RETURN
      END SUBROUTINE INFOG2L
C***********************************************************************
      INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
      INTEGER            INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
        INDXG2P = 0
        WRITE(*,*) 'Error. INFOG2L should not be called.'
        STOP
      RETURN
      END FUNCTION INDXG2P
C***********************************************************************
      SUBROUTINE pcscal(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      COMPLEX            ALPHA
      COMPLEX            X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PCSCAL should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcscal
C***********************************************************************
      SUBROUTINE pzscal(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
C     DOUBLE COMPLEX     ALPHA
C     DOUBLE COMPLEX     X( * )
      COMPLEX(kind=kind(0.0D0)) :: ALPHA, X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PZSCAL should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzscal
C***********************************************************************
      SUBROUTINE pdscal(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PDSCAL should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdscal
C***********************************************************************
      SUBROUTINE psscal(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      REAL               ALPHA
      REAL               X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PSSCAL should not be called.'
        STOP
      RETURN
      END SUBROUTINE psscal
C***********************************************************************
      SUBROUTINE pzdot
     &    ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
      IMPLICIT NONE
      INTEGER N, IX, JX, IY, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
C     DOUBLE COMPLEX X(*), Y(*)
      COMPLEX(kind=kind(0.0D0)) :: X(*), Y(*)
      DOUBLE PRECISION DOT
        DOT = 0.0d0
        WRITE(*,*) 'Error. PZDOT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzdot
C***********************************************************************
      SUBROUTINE pcdot
     &    ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
      IMPLICIT NONE
      INTEGER N, IX, JX, IY, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      COMPLEX X(*), Y(*)
      REAL DOT
        DOT = 0.0e0
        WRITE(*,*) 'Error. PCDOT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcdot
C***********************************************************************
      SUBROUTINE pddot
     &    ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
      IMPLICIT NONE
      INTEGER N, IX, JX, IY, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      DOUBLE PRECISION X(*), Y(*), DOT
        DOT = 0.0d0
        WRITE(*,*) 'Error. PDDOT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pddot
C***********************************************************************
      SUBROUTINE psdot
     &    ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
      IMPLICIT NONE
      INTEGER N, IX, JX, IY, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      REAL X(*), Y(*), DOT
        DOT = 0.0e0
        WRITE(*,*) 'Error. PSDOT should not be called.'
        STOP
      RETURN
      END SUBROUTINE psdot
C***********************************************************************
      SUBROUTINE zgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
C     DOUBLE COMPLEX A(*)
      COMPLEX(kind=kind(0.0D0)) :: A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. ZGEBS2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE zgebs2d
C***********************************************************************
      SUBROUTINE cgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      COMPLEX A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. CGEBS2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE cgebs2d
C***********************************************************************
      SUBROUTINE sgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      REAL A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. SGEBS2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE sgebs2d
C***********************************************************************
      SUBROUTINE dgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      DOUBLE PRECISION A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. DGEBS2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE dgebs2d
C***********************************************************************
      SUBROUTINE zgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
C     DOUBLE COMPLEX A(*)
      COMPLEX(kind=kind(0.0D0)) :: A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. ZGEBR2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE zgebr2d
C***********************************************************************
      SUBROUTINE cgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      COMPLEX A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. CGEBR2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE cgebr2d
C***********************************************************************
      SUBROUTINE sgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      REAL A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. SGEBR2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE sgebr2d
C***********************************************************************
      SUBROUTINE dgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      DOUBLE PRECISION A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. DGEBR2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE dgebr2d
C***********************************************************************
      SUBROUTINE pcgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     &                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX            A( * ), B( * )
        WRITE(*,*) 'Error. PCGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcgetrs
C***********************************************************************
      SUBROUTINE pzgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     &                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
c     DOUBLE COMPLEX     A( * ), B( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
        WRITE(*,*) 'Error. PZGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzgetrs
C***********************************************************************
      SUBROUTINE psgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     &                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      REAL               A( * ), B( * )
        WRITE(*,*) 'Error. PSGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE psgetrs
C***********************************************************************
      SUBROUTINE pdgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     &                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      DOUBLE PRECISION   A( * ), B( * )
        WRITE(*,*) 'Error. PDGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdgetrs
C***********************************************************************
      SUBROUTINE pcpotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     &           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      COMPLEX         A( * ), B( * )
        WRITE(*,*) 'Error. PCPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcpotrs
C***********************************************************************
      SUBROUTINE pzpotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     &           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
c     DOUBLE COMPLEX     A( * ), B( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
        WRITE(*,*) 'Error. PZPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzpotrs
C***********************************************************************
      SUBROUTINE pspotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     &           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      REAL            A( * ), B( * )
        WRITE(*,*) 'Error. PSPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pspotrs
C***********************************************************************
      SUBROUTINE pdpotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     &           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      DOUBLE          PRECISION A( * ), B( * )
        WRITE(*,*) 'Error. PDPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdpotrs
C***********************************************************************
      SUBROUTINE pscnrm2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      REAL NORM2
      COMPLEX X( * )
        WRITE(*,*) 'Error. PCNRM2 should not be called.'
        STOP
      RETURN
      END SUBROUTINE pscnrm2
C***********************************************************************
      SUBROUTINE pdznrm2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      DOUBLE PRECISION NORM2
C     DOUBLE COMPLEX X( * )
      COMPLEX(kind=kind(0.0D0)) :: X( * )
        WRITE(*,*) 'Error. PZNRM2 should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdznrm2
C***********************************************************************
      SUBROUTINE psnrm2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      REAL    NORM2, X( * )
        WRITE(*,*) 'Error. PSNRM2 should not be called.'
        STOP
      RETURN
      END SUBROUTINE psnrm2
C***********************************************************************
      SUBROUTINE pdnrm2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      DOUBLE PRECISION NORM2, X( * )
        WRITE(*,*) 'Error. PDNRM2 should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdnrm2
C***********************************************************************
      REAL FUNCTION pclange( NORM, M, N, A, IA,  JA,
     &                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      COMPLEX      A( * ), WORK( * )
      PCLANGE = 0.0e0
        WRITE(*,*) 'Error. PCLANGE should not be called.'
        STOP
      RETURN
      END FUNCTION pclange
C***********************************************************************
      DOUBLE PRECISION FUNCTION pzlange( NORM, M, N, A, IA,  JA,
     &                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      REAL         A( * ), WORK( * )
      PZLANGE = 0.0d0
        WRITE(*,*) 'Error. PZLANGE should not be called.'
        STOP
      RETURN
      END FUNCTION pzlange
C***********************************************************************
      REAL FUNCTION pslange( NORM, M, N, A, IA,  JA,
     &                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      REAL         A( * ), WORK( * )
      PSLANGE = 0.0e0
        WRITE(*,*) 'Error. PSLANGE should not be called.'
        STOP
      RETURN
      END FUNCTION pslange
C***********************************************************************
      DOUBLE PRECISION FUNCTION pdlange( NORM, M, N, A, IA,  JA,
     &                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      DOUBLE       PRECISION A( * ), WORK( * )
      PDLANGE = 0.0d0
        WRITE(*,*) 'Error. PDLANGE should not be called.'
        STOP
      RETURN
      END FUNCTION pdlange
C***********************************************************************
      SUBROUTINE pcgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     &           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      REAL            ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      COMPLEX         A( * ), WORK( * )
        WRITE(*,*) 'Error. PCGECON should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcgecon
C***********************************************************************
      SUBROUTINE pzgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     &           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      DOUBLE PRECISION ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
C     DOUBLE COMPLEX  A( * ), WORK( * )
      COMPLEX(kind=kind(0.0D0)) :: A( * ), WORK( * )
        WRITE(*,*) 'Error. PZGECON should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzgecon
C***********************************************************************
      SUBROUTINE psgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     &           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      REAL            ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      REAL            A( * ), WORK( * )
        WRITE(*,*) 'Error. PSGECON should not be called.'
        STOP
      RETURN
      END SUBROUTINE psgecon
C***********************************************************************
      SUBROUTINE pdgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     &           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      DOUBLE          PRECISION ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      DOUBLE          PRECISION A( * ), WORK( * )
        WRITE(*,*) 'Error. PDGECON should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdgecon
C***********************************************************************
      SUBROUTINE pcgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     &           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER    IA, JA, INFO, LWORK, M, N
      INTEGER    DESCA( * ), IPIV( * )
      COMPLEX    A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PCGEQPF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcgeqpf
C***********************************************************************
      SUBROUTINE pzgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     &           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER    IA, JA, INFO, LWORK, M, N
      INTEGER    DESCA( * ), IPIV( * )
C     DOUBLE COMPLEX A( * ), TAU( * ), WORK( * )
      COMPLEX(kind=kind(0.0D0)) :: A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PZGEQPF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzgeqpf
C***********************************************************************
      SUBROUTINE psgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     &           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER         IA, JA, INFO, LWORK, M, N
      INTEGER         DESCA( * ), IPIV( * )
      REAL       A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PSGEQPF should not be called.'
        STOP
      RETURN
      END SUBROUTINE psgeqpf
C***********************************************************************
      SUBROUTINE pdgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     &           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER         IA, JA, INFO, LWORK, M, N
      INTEGER         DESCA( * ), IPIV( * )
      DOUBLE PRECISION A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PDGEQPF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdgeqpf
C***********************************************************************
      SUBROUTINE pcaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     &           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      COMPLEX A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PCAXPY should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcaxpy
C***********************************************************************
      SUBROUTINE pzaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     &           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
C     DOUBLE COMPLEX A(*),X(*),Y(*)
      COMPLEX(kind=kind(0.0D0)) :: A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PZAXPY should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzaxpy
C***********************************************************************
      SUBROUTINE psaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     &           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      REAL A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PSAXPY should not be called.'
        STOP
      RETURN
      END SUBROUTINE psaxpy
C***********************************************************************
      SUBROUTINE pdaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     &           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      DOUBLE PRECISION A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PDAXPY should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdaxpy
C***********************************************************************
      SUBROUTINE pctrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      COMPLEX            ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * )
        WRITE(*,*) 'Error. PCTRSM should not be called.'
        STOP
      RETURN
      END SUBROUTINE pctrsm 
C***********************************************************************
      SUBROUTINE pztrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
C     DOUBLE COMPLEX     ALPHA
      COMPLEX(kind=kind(0.0D0)) ::     ALPHA
      INTEGER            DESCA( * ), DESCB( * )
C     DOUBLE COMPLEX     A( * ), B( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
        WRITE(*,*) 'Error. PZTRSM should not be called.'
        STOP
      RETURN
      END SUBROUTINE pztrsm 
C***********************************************************************
      SUBROUTINE pstrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      REAL               ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
        WRITE(*,*) 'Error. PSTRSM should not be called.'
        STOP
      RETURN
      END SUBROUTINE pstrsm 
C***********************************************************************
      SUBROUTINE pdtrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      DOUBLE PRECISION   ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * )
        WRITE(*,*) 'Error. PDTRSM should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdtrsm 
C***********************************************************************
      SUBROUTINE pcunmqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     &                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     &                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER   DESCA( * ), DESCC( * )
      COMPLEX   A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PCUNMQR should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcunmqr
C***********************************************************************
      SUBROUTINE pzunmqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     &                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     &                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER   DESCA( * ), DESCC( * )
C     DOUBLE COMPLEX A(  *  ), C( * ), TAU( * ), WORK( * )
      COMPLEX(kind=kind(0.0D0)) :: A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PZUNMQR should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzunmqr
C***********************************************************************
      SUBROUTINE psormqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     &                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     &                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER   DESCA( * ), DESCC( * )
      REAL      A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PSORMQR should not be called.'
        STOP
      RETURN
      END SUBROUTINE psormqr
C***********************************************************************
      SUBROUTINE pdormqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     &                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     &                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER         IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER         DESCA( * ), DESCC( * )
      DOUBLE PRECISION  A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PDORMQR should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdormqr
C***********************************************************************
      SUBROUTINE chk1mat( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     &                    DESCAPOS0, INFO )
      IMPLICIT NONE
      INTEGER            DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA, NAPOS0
      INTEGER            DESCA( * )
        WRITE(*,*) 'Error. CHK1MAT should not be called.'
        STOP
      RETURN
      END SUBROUTINE chk1mat
C***********************************************************************
      SUBROUTINE pchk2mat( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     &                     DESCAPOS0, MB, MBPOS0, NB, NBPOS0, IB, JB,
     &                     DESCB, DESCBPOS0, NEXTRA, EX, EXPOS, INFO )
      IMPLICIT NONE
      INTEGER            DESCAPOS0, DESCBPOS0, IA, IB, INFO, JA, JB, MA,
     &                   MAPOS0, MB, MBPOS0, NA, NAPOS0, NB, NBPOS0,
     &                   NEXTRA
      INTEGER            DESCA( * ), DESCB( * ), EX( NEXTRA ),
     &                   EXPOS( NEXTRA )
        WRITE(*,*) 'Error. PCHK2MAT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pchk2mat
C***********************************************************************
      SUBROUTINE pxerbla( CONTXT, SRNAME, INFO )
      IMPLICIT NONE
      INTEGER CONTXT, INFO
      CHARACTER SRNAME
        WRITE(*,*) 'Error. PXERBLA should not be called.'
        STOP
      RETURN
      END SUBROUTINE pxerbla
C***********************************************************************
      SUBROUTINE descset( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
     &                    LLD )
      IMPLICIT NONE
      INTEGER            ICSRC, ICTXT, IRSRC, LLD, M, MB, N, NB
      INTEGER            DESC( * )
        WRITE(*,*) 'Error. DESCSET should not be called.'
        STOP
      RETURN
      END SUBROUTINE descset

