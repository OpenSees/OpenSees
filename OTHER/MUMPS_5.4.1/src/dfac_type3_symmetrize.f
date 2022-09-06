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
      SUBROUTINE DMUMPS_SYMMETRIZE( BUF, BLOCK_SIZE,
     &                           MYROW, MYCOL, NPROW, NPCOL,
     &                           A, LOCAL_M, LOCAL_N, N, MYID, COMM )
      IMPLICIT NONE
      INTEGER BLOCK_SIZE, NPROW, NPCOL, LOCAL_M, LOCAL_N, N, COMM
      INTEGER MYROW, MYCOL, MYID
      DOUBLE PRECISION BUF( BLOCK_SIZE * BLOCK_SIZE )
      DOUBLE PRECISION A( LOCAL_M, LOCAL_N )
      INTEGER NBLOCK, IBLOCK, JBLOCK, IBLOCK_SIZE, JBLOCK_SIZE
      INTEGER ROW_SOURCE, ROW_DEST, COL_SOURCE, COL_DEST
      INTEGER IGLOB, JGLOB
      INTEGER IROW_LOC_SOURCE, JCOL_LOC_SOURCE
      INTEGER IROW_LOC_DEST, JCOL_LOC_DEST
      INTEGER PROC_SOURCE, PROC_DEST
      NBLOCK = ( N - 1 ) / BLOCK_SIZE + 1
      DO IBLOCK = 1, NBLOCK
        IF ( IBLOCK .NE. NBLOCK
     &    ) THEN
          IBLOCK_SIZE = BLOCK_SIZE
        ELSE
          IBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
        END IF
        ROW_SOURCE = mod( IBLOCK - 1, NPROW ) 
        COL_DEST   = mod( IBLOCK - 1, NPCOL )
        IGLOB = ( IBLOCK - 1 ) * BLOCK_SIZE + 1
        IROW_LOC_SOURCE = BLOCK_SIZE *
     &                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     &                  + mod( IGLOB - 1, BLOCK_SIZE ) + 1
        JCOL_LOC_DEST   = BLOCK_SIZE *
     &                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     &                  + mod( IGLOB - 1, BLOCK_SIZE ) + 1
        DO JBLOCK = 1, IBLOCK
          IF ( JBLOCK .NE. NBLOCK
     &      ) THEN
            JBLOCK_SIZE = BLOCK_SIZE
          ELSE
            JBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
          END IF
          COL_SOURCE = mod( JBLOCK - 1, NPCOL )
          ROW_DEST   = mod( JBLOCK - 1, NPROW )
          PROC_SOURCE = ROW_SOURCE * NPCOL + COL_SOURCE
          PROC_DEST   = ROW_DEST   * NPCOL + COL_DEST
          IF ( PROC_SOURCE .eq. PROC_DEST ) THEN
           IF ( MYID .eq. PROC_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     &                  ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     &                  + mod( JGLOB - 1, BLOCK_SIZE ) + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     &                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     &                  + mod( JGLOB - 1, BLOCK_SIZE ) + 1
            IF ( IBLOCK .eq. JBLOCK ) THEN
               IF ( IBLOCK_SIZE .ne. JBLOCK_SIZE ) THEN
                WRITE(*,*) MYID,': Error in calling transdiag:unsym'
                CALL MUMPS_ABORT()
              END IF
              CALL DMUMPS_TRANS_DIAG( A( IROW_LOC_SOURCE,
     &                 JCOL_LOC_SOURCE),
     &                 IBLOCK_SIZE, LOCAL_M )
            ELSE
              CALL DMUMPS_TRANSPO(
     &           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ),
     &           A( IROW_LOC_DEST, JCOL_LOC_DEST ),
     &           IBLOCK_SIZE, JBLOCK_SIZE, LOCAL_M )
            END IF
           END IF
          ELSE IF (  MYROW .eq. ROW_SOURCE 
     &    .AND. MYCOL .eq. COL_SOURCE ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     &                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     &                  + mod( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL DMUMPS_SEND_BLOCK( BUF,
     &           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ), LOCAL_M,
     &           IBLOCK_SIZE, JBLOCK_SIZE, COMM, PROC_DEST )
          ELSE IF ( MYROW .eq. ROW_DEST 
     &    .AND.     MYCOL .eq. COL_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     &                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     &                  + mod( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL DMUMPS_RECV_BLOCK( BUF,
     &           A( IROW_LOC_DEST, JCOL_LOC_DEST ), LOCAL_M,
     &           JBLOCK_SIZE, IBLOCK_SIZE, COMM, PROC_SOURCE )
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE DMUMPS_SYMMETRIZE
      SUBROUTINE DMUMPS_SEND_BLOCK( BUF, A, LDA, M, N, COMM, DEST )
      IMPLICIT NONE
      INTEGER M, N, LDA, DEST, COMM
      DOUBLE PRECISION BUF(*), A(LDA,*)
      INTEGER I, IBUF, IERR
      INTEGER J
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      IBUF = 1
      DO J = 1, N
        BUF( IBUF: IBUF + M - 1 ) = A( 1 : M, J )
        DO I = 1, M
        END DO
        IBUF = IBUF + M
      END DO
      CALL MPI_SEND( BUF, M * N, MPI_DOUBLE_PRECISION,
     &     DEST, SYMMETRIZE, COMM, IERR )
      RETURN
      END SUBROUTINE DMUMPS_SEND_BLOCK
      SUBROUTINE DMUMPS_RECV_BLOCK( BUF, A, LDA, M, N, COMM, SOURCE )
      IMPLICIT NONE
      INTEGER LDA, M, N, COMM, SOURCE
      DOUBLE PRECISION BUF(*), A( LDA, *)
      INTEGER I, IBUF, IERR
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      CALL MPI_RECV( BUF(1), M * N, MPI_DOUBLE_PRECISION, SOURCE,
     &               SYMMETRIZE, COMM, STATUS, IERR )
      IBUF = 1
      DO I = 1, M
        CALL dcopy( N, BUF(IBUF), 1, A(I,1), LDA )
        IBUF = IBUF + N
      END DO
      RETURN
      END SUBROUTINE DMUMPS_RECV_BLOCK
      SUBROUTINE DMUMPS_TRANS_DIAG( A, N, LDA )
      IMPLICIT NONE
      INTEGER N,LDA
      DOUBLE PRECISION A( LDA, * )
      INTEGER I, J
      DO I = 2, N
        DO J = 1, I - 1
          A( J, I ) = A( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE DMUMPS_TRANS_DIAG
      SUBROUTINE DMUMPS_TRANSPO( A1, A2, M, N, LD )
      IMPLICIT NONE
      INTEGER M,N,LD
      DOUBLE PRECISION A1( LD,* ), A2( LD, * )
      INTEGER I, J
      DO J = 1, N
        DO I = 1, M
          A2( J, I ) = A1( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE DMUMPS_TRANSPO
