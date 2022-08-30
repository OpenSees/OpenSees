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
      SUBROUTINE DMUMPS_ELT_DISTRIB(
     &            N, NELT, NA_ELT8,
     &            COMM, MYID, SLAVEF,
     &            IELPTR_LOC8, RELPTR_LOC8,
     &            ELTVAR_LOC, ELTVAL_LOC,
     &            LINTARR, LDBLARR,
     &            KEEP,KEEP8, MAXELT_SIZE,
     &            FRTPTR, FRTELT, A, LA, FILS,
     &            id, root )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N, NELT
      INTEGER(8) :: NA_ELT8
      INTEGER COMM, MYID, SLAVEF, MAXELT_SIZE, MSGLEN
      INTEGER(8), intent(IN) :: LA
      INTEGER FRTPTR( N+1 )
      INTEGER FRTELT( NELT ), FILS ( N )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8), INTENT(IN)    :: IELPTR_LOC8( NELT + 1 )
      INTEGER(8), INTENT(INOUT) :: RELPTR_LOC8( NELT + 1 )
      INTEGER(8), INTENT(IN) :: LINTARR, LDBLARR
      INTEGER ELTVAR_LOC( LINTARR )
      DOUBLE PRECISION ELTVAL_LOC( LDBLARR )
      DOUBLE PRECISION A( LA )
      TYPE(DMUMPS_STRUC)     :: id
      TYPE(DMUMPS_ROOT_STRUC) :: root
      INTEGER numroc
      EXTERNAL numroc
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR_MPI
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER MSGTAG
      INTEGER allocok
      INTEGER I, DEST, MAXELT_REAL_SIZE, MPG, IEL, SIZEI, SIZER
      INTEGER NBRECORDS, NBUF
      INTEGER(8) :: RECV_IELTPTR8
      INTEGER(8) :: RECV_RELTPTR8
      INTEGER INODE
      INTEGER(8) :: IELTPTR8, RELTPTR8
      LOGICAL FINI, PROKG, I_AM_SLAVE, EARLYT3ROOTINS
      INTEGER(8) :: PTR_ROOT
      INTEGER LOCAL_M, LOCAL_N, LP, IBEG, IGLOB, JGLOB
      INTEGER ARROW_ROOT
      INTEGER IELT, J, NB_REC, IREC
      INTEGER(8) :: K8, IVALPTR8
      INTEGER ILOCROOT, JLOCROOT, IPOSROOT, JPOSROOT, IPTR
      INTEGER JCOL_GRID, IROW_GRID
      INTEGER NBELROOT
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      DOUBLE PRECISION  VAL
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER, DIMENSION( :, : ), ALLOCATABLE :: BUFI
      DOUBLE PRECISION, DIMENSION( :, : ), ALLOCATABLE :: BUFR
      DOUBLE PRECISION, DIMENSION( : ), ALLOCATABLE :: TEMP_ELT_R
      INTEGER, DIMENSION( : ), ALLOCATABLE :: TEMP_ELT_I
      INTEGER(8), DIMENSION( : ), ALLOCATABLE :: ELROOTPOS8
      INTEGER, DIMENSION( : ), ALLOCATABLE, TARGET :: RG2LALLOC
      INTEGER, DIMENSION( : ), POINTER     :: RG2L
      MPG = id%ICNTL(3)
      LP  = id%ICNTL(1)
      I_AM_SLAVE = ( KEEP(46) .eq. 1 .or. MYID .ne.MASTER )
      PROKG = ( MPG > 0 .and. MYID .eq. MASTER )
      PROKG   = (PROKG.AND.(id%ICNTL(4).GE.2))
      KEEP(49) = 0
      ARROW_ROOT = 0
      EARLYT3ROOTINS = KEEP(200) .EQ.0
      IF ( MYID .eq. MASTER ) THEN
        IF ( KEEP(46) .eq. 0 ) THEN
          NBUF = SLAVEF
        ELSE
          NBUF = SLAVEF - 1
        END IF
        NBRECORDS = KEEP(39)
        IF (NA_ELT8 < int(NBRECORDS,8)) THEN
          NBRECORDS = int(NA_ELT8)
        ENDIF
        IF ( KEEP(50) .eq. 0 ) THEN
          MAXELT_REAL_SIZE = MAXELT_SIZE * MAXELT_SIZE
        ELSE
          MAXELT_REAL_SIZE = MAXELT_SIZE * (MAXELT_SIZE+1)/2
        END IF
        IF ( MAXELT_REAL_SIZE .GT. KEEP(39) ) THEN
          NBRECORDS = MAXELT_REAL_SIZE
          IF ( MPG .GT. 0 ) THEN
            WRITE(MPG,*)
     & ' ** Warning : For element distrib NBRECORDS set to ',
     & MAXELT_REAL_SIZE,' because one element is large'
          END IF
        END IF
        ALLOCATE( BUFI( 2*NBRECORDS+1, NBUF ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = 2*NBRECORDS + 1
          GOTO 100
        END IF
        ALLOCATE( BUFR( NBRECORDS+1, NBUF ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = NBRECORDS + 1
          GOTO 100
        END IF
        IF ( KEEP(52) .ne. 0 ) THEN
          ALLOCATE( TEMP_ELT_R( MAXELT_REAL_SIZE ), stat =allocok )
          IF ( allocok .gt. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = MAXELT_REAL_SIZE
            GOTO 100
          END IF
        END IF
        ALLOCATE( TEMP_ELT_I( MAXELT_SIZE ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = MAXELT_SIZE
            GOTO 100
        END IF
        IF ( KEEP(38) .ne. 0 ) THEN
          NBELROOT = FRTPTR(KEEP(38)+1)-FRTPTR(KEEP(38))
          IF ( EARLYT3ROOTINS ) THEN
            ALLOCATE( ELROOTPOS8( max(NBELROOT,1) ),
     &              stat = allocok )
            IF ( allocok .gt. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = NBELROOT
              GOTO 100
            END IF
          ENDIF
          IF (KEEP(46) .eq. 0 ) THEN
           ALLOCATE( RG2LALLOC( N ), stat = allocok )
           IF ( allocok .gt. 0 ) THEN
               id%INFO(1) = -13
               id%INFO(2) = N
               GOTO 100
           END IF
           INODE = KEEP(38)
           I     = 1
           DO WHILE ( INODE .GT. 0 )
             RG2LALLOC( INODE ) = I
             INODE = FILS( INODE )
             I = I + 1
           END DO
           RG2L => RG2LALLOC
          ELSE 
           RG2L => root%RG2L_ROW
          END IF
        END IF
        DO I = 1, NBUF
          BUFI( 1, I ) = 0
          BUFR( 1, I ) = ZERO
        END DO
      END IF
 100  CONTINUE
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1), COMM, MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST( NBRECORDS, 1, MPI_INTEGER, MASTER,
     &                COMM, IERR_MPI )
      RECV_IELTPTR8 = 1_8
      RECV_RELTPTR8 = 1_8
      IF ( MYID .eq. MASTER ) THEN
        NBELROOT = 0
        RELTPTR8 = 1_8
        RELPTR_LOC8(1) = 1
        DO IEL = 1, NELT
          IELTPTR8 = int(id%ELTPTR( IEL ),8)
          SIZEI   = int(int(id%ELTPTR( IEL + 1 ),8) - IELTPTR8)
          IF ( KEEP( 50 ) .eq. 0 ) THEN
            SIZER = SIZEI * SIZEI
          ELSE
            SIZER = SIZEI * ( SIZEI + 1 ) / 2
          END IF
          DEST = id%ELTPROC( IEL )
          IF ( DEST .eq. -2 ) THEN
            NBELROOT   = NBELROOT   + 1
            FRTELT( FRTPTR(KEEP(38)) + NBELROOT - 1 ) = IEL
            ELROOTPOS8( NBELROOT ) = RELTPTR8
            GOTO 200
          END IF
          IF ( DEST .ge. 0 .and. KEEP(46) .eq. 0 ) DEST = DEST + 1
          IF ( KEEP(52) .ne. 0 ) THEN
            CALL DMUMPS_SCALE_ELEMENT( N, SIZEI, SIZER,
     &               id%ELTVAR( IELTPTR8 ), id%A_ELT( RELTPTR8 ),
     &               TEMP_ELT_R(1), MAXELT_REAL_SIZE,
     &               id%ROWSCA(1), id%COLSCA(1), KEEP(50) )
          END IF
          IF ( DEST .eq. 0 .or. ( DEST .eq. -1 .and. KEEP(46) .ne. 0 ) )
     &      THEN
            ELTVAR_LOC( RECV_IELTPTR8: RECV_IELTPTR8 + SIZEI - 1 )
     &      = id%ELTVAR( IELTPTR8: IELTPTR8 + SIZEI - 1 )
            RECV_IELTPTR8 = RECV_IELTPTR8 + SIZEI
            IF ( KEEP(52) .ne. 0 ) THEN
              ELTVAL_LOC( RECV_RELTPTR8: RECV_RELTPTR8 + SIZER - 1)
     &        = TEMP_ELT_R( 1: SIZER )
              RECV_RELTPTR8 = RECV_RELTPTR8 + SIZER
            END IF
          END IF
          IF ( DEST .NE. 0 .AND. DEST. NE. -3 ) THEN
            IF ( KEEP(52) .eq. 0 ) THEN
              CALL DMUMPS_ELT_FILL_BUF(
     &           id%ELTVAR(IELTPTR8),
     &           id%A_ELT (RELTPTR8),
     &           SIZEI, SIZER,
     &
     &           DEST, NBUF, NBRECORDS,
     &           BUFI, BUFR, COMM )
            ELSE
              CALL DMUMPS_ELT_FILL_BUF(
     &           id%ELTVAR(IELTPTR8),
     &           TEMP_ELT_R( 1 ),
     &           SIZEI, SIZER,
     &
     &           DEST, NBUF, NBRECORDS,
     &           BUFI, BUFR, COMM )
            END IF
          END IF
 200      CONTINUE
          RELTPTR8 = RELTPTR8 + SIZER
          IF ( KEEP(46) .eq. 0 .OR. KEEP(52) .eq. 0 ) THEN
            RELPTR_LOC8( IEL + 1 ) = RELTPTR8
          ELSE
            RELPTR_LOC8( IEL + 1 ) = RECV_RELTPTR8
          ENDIF
        END DO
        IF ( KEEP(46) .eq. 0 .OR. KEEP(52) .eq. 0 ) THEN
          KEEP8(26) = RELTPTR8 - 1_8
        ELSE
          KEEP8(26) = RECV_RELTPTR8 - 1_8
        ENDIF
        IF ( RELTPTR8 - 1_8 .NE. NA_ELT8 ) THEN
          WRITE(*,*) " ** Internal error in DMUMPS_ELT_DISTRIB",
     &               RELTPTR8 - 1_8, NA_ELT8
          CALL MUMPS_ABORT()
        END IF
        DEST = -2
        IELTPTR8 = 1_8
        RELTPTR8 = 1_8
        SIZEI   = 1
        SIZER   = 1
        CALL DMUMPS_ELT_FILL_BUF(
     &           id%ELTVAR(IELTPTR8),
     &           id%A_ELT (RELTPTR8),
     &           SIZEI, SIZER,
     &
     &           DEST, NBUF, NBRECORDS,
     &           BUFI, BUFR, COMM )
        IF ( KEEP(52) .NE. 0 ) DEALLOCATE( TEMP_ELT_R )
      ELSE
        FINI = ( RECV_IELTPTR8 .eq. IELPTR_LOC8( NELT+1 )
     &     .and. RECV_RELTPTR8 .eq. RELPTR_LOC8( NELT+1 ) )
        DO WHILE ( .not. FINI )
          CALL MPI_PROBE( MASTER, MPI_ANY_TAG,
     &                    COMM, STATUS, IERR_MPI )
          MSGTAG = STATUS( MPI_TAG    )
          SELECT CASE ( MSGTAG )
             CASE( ELT_INT )
               CALL MPI_GET_COUNT( STATUS, MPI_INTEGER,
     &                             MSGLEN, IERR_MPI )
               CALL MPI_RECV( ELTVAR_LOC( RECV_IELTPTR8 ), MSGLEN,
     &            MPI_INTEGER, MASTER, ELT_INT,
     &            COMM, STATUS, IERR_MPI )
               RECV_IELTPTR8 = RECV_IELTPTR8 + MSGLEN
             CASE( ELT_REAL )
                CALL MPI_GET_COUNT( STATUS, MPI_DOUBLE_PRECISION,
     &                              MSGLEN, IERR_MPI )
                CALL MPI_RECV( ELTVAL_LOC( RECV_RELTPTR8 ), MSGLEN,
     &            MPI_DOUBLE_PRECISION, MASTER, ELT_REAL,
     &            COMM, STATUS, IERR_MPI )
                RECV_RELTPTR8 = RECV_RELTPTR8 + MSGLEN
          END SELECT
          FINI = ( RECV_IELTPTR8 .eq. IELPTR_LOC8( NELT+1 )
     &       .and. RECV_RELTPTR8 .eq. RELPTR_LOC8( NELT+1 ) )
        END DO
      END IF
      IF ( KEEP(38) .NE. 0 .AND. EARLYT3ROOTINS ) THEN
        IF ( I_AM_SLAVE .and. root%yes ) THEN
          CALL DMUMPS_GET_ROOT_INFO(root,
     &                              LOCAL_M, LOCAL_N, PTR_ROOT, LA)
          CALL DMUMPS_SET_ROOT_TO_ZERO(root, KEEP, A, LA)
        END IF
        IF ( MYID .NE. MASTER ) THEN
          ALLOCATE( BUFI( NBRECORDS * 2 + 1, 1 ), stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = NBRECORDS * 2 + 1
            GOTO 250
          END IF
          ALLOCATE( BUFR( NBRECORDS, 1 )        , stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = NBRECORDS
          END IF
        END IF
 250    CONTINUE
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1), COMM, MYID )
        IF ( id%INFO(1) .LT. 0 ) RETURN
        IF ( MYID .eq. MASTER ) THEN
        DO IPTR = FRTPTR(KEEP(38)), FRTPTR(KEEP(38)+1) - 1
          IELT = FRTELT( IPTR )
          SIZEI = id%ELTPTR( IELT + 1 ) - id%ELTPTR( IELT )
          DO I = 1, SIZEI
            TEMP_ELT_I( I ) = RG2L
     &              ( id%ELTVAR( id%ELTPTR(IELT) + I - 1 ) )
          END DO
          IVALPTR8 = ELROOTPOS8( IPTR - FRTPTR(KEEP(38)) + 1 ) - 1
          K8 = 1_8
          DO J = 1, SIZEI
            JGLOB = id%ELTVAR( id%ELTPTR( IELT ) + J - 1 )
            IF ( KEEP(50).eq. 0 ) THEN
              IBEG = 1
            ELSE
              IBEG = J
            END IF
            DO I = IBEG, SIZEI
              IGLOB = id%ELTVAR( id%ELTPTR( IELT ) + I - 1 )
              IF ( KEEP(52) .eq. 0 ) THEN
                VAL = id%A_ELT( IVALPTR8 + K8 )
              ELSE
                VAL = id%A_ELT( IVALPTR8 + K8 ) *
     &                id%ROWSCA( IGLOB ) * id%COLSCA( JGLOB )
              END IF
              IF ( KEEP(50).eq.0 ) THEN
                IPOSROOT = TEMP_ELT_I( I )
                JPOSROOT = TEMP_ELT_I( J )
              ELSE
                IF ( TEMP_ELT_I(I) .GT. TEMP_ELT_I(J) ) THEN
                  IPOSROOT = TEMP_ELT_I(I)
                  JPOSROOT = TEMP_ELT_I(J)
                ELSE
                  IPOSROOT = TEMP_ELT_I(J)
                  JPOSROOT = TEMP_ELT_I(I)
                END IF
              END IF
              IROW_GRID = mod( ( IPOSROOT - 1 )/root%MBLOCK,
     &                           root%NPROW )
              JCOL_GRID = mod( ( JPOSROOT - 1 )/root%NBLOCK,
     &                           root%NPCOL )
              IF ( KEEP(46) .eq. 0 ) THEN
                DEST = IROW_GRID * root%NPCOL + JCOL_GRID + 1
              ELSE
                DEST = IROW_GRID * root%NPCOL + JCOL_GRID
              END IF
              IF ( DEST .eq. MASTER ) THEN
                ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
                JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
              ARROW_ROOT = ARROW_ROOT + 1
              IF (KEEP(60)==0) THEN
                A( PTR_ROOT
     &             + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &             + int(ILOCROOT - 1,8) )
     &          =  A( PTR_ROOT
     &             + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &             + int(ILOCROOT - 1,8) )
     &          + VAL
              ELSE
                root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                            * int(root%SCHUR_LLD,8)
     &                            + int(ILOCROOT,8) )
     &          = root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                            * int(root%SCHUR_LLD,8)
     &                            + int(ILOCROOT,8) )
     &          + VAL
              ENDIF
              ELSE
                CALL DMUMPS_ARROW_FILL_SEND_BUF(
     &          IPOSROOT, JPOSROOT, VAL, DEST, BUFI, BUFR, NBRECORDS,
     &          NBUF, LP, COMM, KEEP(46) )
              END IF
              K8 = K8 + 1_8
            END DO
          END DO
        END DO
        CALL DMUMPS_ARROW_FINISH_SEND_BUF(
     &          BUFI, BUFR, NBRECORDS,
     &          NBUF, LP, COMM, KEEP(46) )
        ELSE
          FINI = .FALSE.
          DO WHILE ( .not. FINI )
            CALL MPI_RECV( BUFI(1,1), 2*NBRECORDS+1,
     &                MPI_INTEGER, MASTER,
     &                ARROWHEAD,
     &                COMM, STATUS, IERR_MPI )
            NB_REC = BUFI(1,1)
            ARROW_ROOT = ARROW_ROOT + NB_REC
            IF (NB_REC.LE.0) THEN
              FINI = .TRUE.
              NB_REC = -NB_REC
            ENDIF
            IF (NB_REC.EQ.0) EXIT
            CALL MPI_RECV( BUFR(1,1), NBRECORDS, MPI_DOUBLE_PRECISION,
     &                     MASTER, ARROWHEAD,
     &                     COMM, STATUS, IERR_MPI )
            DO IREC = 1, NB_REC
              IPOSROOT = BUFI( IREC * 2, 1 )
              JPOSROOT = BUFI( IREC * 2 + 1, 1 )
              VAL      = BUFR( IREC, 1 )
              ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
              JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
              IF (KEEP(60).eq.0) THEN
                 A( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &                       + int(ILOCROOT-1,8))
     &        =  A( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &                       + int(ILOCROOT-1,8))
     &           + VAL
              ELSE
                root%SCHUR_POINTER(int(JLOCROOT-1,8)
     &                         * int(root%SCHUR_LLD,8)
     &                         + int(ILOCROOT,8) )
     &        = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                         * int(root%SCHUR_LLD,8)
     &                         + int(ILOCROOT,8))
     &          + VAL
              ENDIF
            END DO
          END DO
          DEALLOCATE( BUFI )
          DEALLOCATE( BUFR )
        END IF
      END IF
      IF ( MYID .eq. MASTER ) THEN
        DEALLOCATE( BUFI )
        DEALLOCATE( BUFR )
        IF (allocated(ELROOTPOS8)) DEALLOCATE(ELROOTPOS8)
        IF (KEEP(38).ne.0) THEN 
          IF (KEEP(46) .eq. 0 ) THEN
             DEALLOCATE(RG2LALLOC)
          ENDIF
        ENDIF
        DEALLOCATE( TEMP_ELT_I )
      END IF
      KEEP(49) = ARROW_ROOT
      RETURN
      END SUBROUTINE DMUMPS_ELT_DISTRIB
      SUBROUTINE DMUMPS_ELT_FILL_BUF(
     &       ELNODES, ELVAL, SIZEI, SIZER,
     &       DEST, NBUF, NBRECORDS, BUFI, BUFR, COMM )
      IMPLICIT NONE
      INTEGER SIZEI, SIZER, DEST, NBUF, NBRECORDS, COMM
      INTEGER ELNODES( SIZEI ), BUFI( 2*NBRECORDS + 1, NBUF )
      DOUBLE PRECISION ELVAL( SIZER ), BUFR( NBRECORDS + 1, NBUF )
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER I, IBEG, IEND, IERR_MPI, NBRECR
      INTEGER NBRECI
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      IF ( DEST .lt. 0 ) THEN
        IBEG = 1
        IEND = NBUF
      ELSE
        IBEG = DEST
        IEND = DEST
      END IF
      DO I = IBEG, IEND
        NBRECI = BUFI(1,I)
        IF ( NBRECI .ne.0  .and.
     &       ( DEST.eq.-2 .or.
     &         NBRECI + SIZEI .GT. 2*NBRECORDS ) ) THEN
           CALL MPI_SEND( BUFI(2, I), NBRECI, MPI_INTEGER,
     &                    I, ELT_INT, COMM, IERR_MPI )
           BUFI(1,I) = 0
           NBRECI    = 0
        END IF
        NBRECR = int(dble(BUFR(1,I))+0.5D0)
        IF ( NBRECR .ne.0  .and.
     &       ( DEST.eq.-2 .or.
     &         NBRECR + SIZER .GT. NBRECORDS ) ) THEN
           CALL MPI_SEND( BUFR(2, I), NBRECR, MPI_DOUBLE_PRECISION,
     &                    I, ELT_REAL, COMM, IERR_MPI )
           BUFR(1,I) = ZERO
           NBRECR    = 0
        END IF
        IF ( DEST .ne. -2 ) THEN
          BUFI( 2 + NBRECI : 2 + NBRECI + SIZEI - 1, I ) =
     &    ELNODES( 1: SIZEI )
          BUFR( 2 + NBRECR : 2 + NBRECR + SIZER - 1, I ) =
     &    ELVAL( 1: SIZER )
          BUFI(1,I) = NBRECI + SIZEI
          BUFR(1,I) = dble( NBRECR + SIZER )
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_ELT_FILL_BUF
      SUBROUTINE DMUMPS_MAXELT_SIZE( ELTPTR, NELT, MAXELT_SIZE )
      INTEGER NELT, MAXELT_SIZE
      INTEGER ELTPTR( NELT + 1 )
      INTEGER I, S
      MAXELT_SIZE = 0
      DO I = 1, NELT
        S = ELTPTR( I + 1 ) - ELTPTR( I )
        MAXELT_SIZE = max( S, MAXELT_SIZE )
      END DO
      RETURN
      END SUBROUTINE DMUMPS_MAXELT_SIZE
      SUBROUTINE DMUMPS_SCALE_ELEMENT( N, SIZEI, SIZER,
     &               ELTVAR, ELTVAL,
     &               SELTVAL, LSELTVAL,
     &               ROWSCA, COLSCA, K50 )
      INTEGER N, SIZEI, SIZER, LSELTVAL, K50
      INTEGER ELTVAR( SIZEI )
      DOUBLE PRECISION ELTVAL( SIZER )
      DOUBLE PRECISION SELTVAL( LSELTVAL )
      DOUBLE PRECISION ROWSCA( N ), COLSCA( N )
      INTEGER I, J, K
      K = 1
      IF ( K50 .eq. 0 ) THEN
        DO J = 1, SIZEI
          DO I = 1, SIZEI
            SELTVAL(K) = ELTVAL(K) *
     &                   ROWSCA(ELTVAR(I)) *
     &                   COLSCA(ELTVAR(J))
            K = K + 1
          END DO
        END DO
      ELSE
        DO J = 1, SIZEI
          DO I = J, SIZEI
            SELTVAL(K) = ELTVAL(K) *
     &                   ROWSCA(ELTVAR(I)) *
     &                   COLSCA(ELTVAR(J))
            K = K + 1
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_SCALE_ELEMENT
