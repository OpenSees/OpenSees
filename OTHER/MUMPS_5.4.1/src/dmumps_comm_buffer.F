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
        MODULE DMUMPS_BUF
        PRIVATE
        PUBLIC :: DMUMPS_BUF_TRY_FREE_CB, DMUMPS_BUF_INIT,
     &   DMUMPS_BUF_INI_MYID,
     &   DMUMPS_BUF_ALLOC_CB ,       DMUMPS_BUF_DEALL_CB ,
     &   DMUMPS_BUF_ALLOC_SMALL_BUF, DMUMPS_BUF_DEALL_SMALL_BUF,
     &   DMUMPS_BUF_ALLOC_LOAD_BUFFER,DMUMPS_BUF_DEALL_LOAD_BUFFER,
     &   DMUMPS_BUF_SEND_CB,     DMUMPS_BUF_SEND_VCB,
     &   DMUMPS_BUF_SEND_1INT,       DMUMPS_BUF_SEND_DESC_BANDE,
     &   DMUMPS_BUF_SEND_MAPLIG, DMUMPS_BUF_SEND_MAITRE2,
     &   DMUMPS_BUF_SEND_CONTRIB_TYPE2,
     &   DMUMPS_BUF_SEND_BLOCFACTO, DMUMPS_BUF_SEND_BLFAC_SLAVE,
     &   DMUMPS_BUF_SEND_MASTER2SLAVE,
     &   DMUMPS_BUF_SEND_CONTRIB_TYPE3, DMUMPS_BUF_SEND_RTNELIND,
     &   DMUMPS_BUF_SEND_ROOT2SLAVE, DMUMPS_BUF_SEND_ROOT2SON,
     &   DMUMPS_BUF_SEND_BACKVEC,DMUMPS_BUF_SEND_UPDATE_LOAD, 
     &   DMUMPS_BUF_DIST_IRECV_SIZE,
     &   DMUMPS_BUF_BCAST_ARRAY, DMUMPS_BUF_ALL_EMPTY,
     &   DMUMPS_BUF_BROADCAST, DMUMPS_BUF_SEND_NOT_MSTR,
     &   DMUMPS_BUF_SEND_FILS ,DMUMPS_BUF_DEALL_MAX_ARRAY
     &   ,DMUMPS_BUF_MAX_ARRAY_MINSIZE
     &   ,DMUMPS_BUF_TEST
         PUBLIC :: DMUMPS_BLR_PACK_CB_LRB
     &   ,DMUMPS_MPI_PACK_LRB
     &   ,DMUMPS_MPI_UNPACK_LRB
        INTEGER NEXT, REQ, CONTENT, OVHSIZE
        PARAMETER( NEXT = 0, REQ = 1, CONTENT = 2, OVHSIZE = 2 )
        INTEGER, SAVE :: SIZEofINT, SIZEofREAL, BUF_MYID
        TYPE DMUMPS_COMM_BUFFER_TYPE
          INTEGER LBUF, HEAD, TAIL,LBUF_INT, ILASTMSG
          INTEGER, DIMENSION(:),POINTER :: CONTENT
        END TYPE DMUMPS_COMM_BUFFER_TYPE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ), SAVE :: BUF_CB
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ), SAVE :: BUF_SMALL
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ), SAVE :: BUF_LOAD
        INTEGER, SAVE :: SIZE_RBUF_BYTES
        INTEGER, SAVE ::  BUF_LMAX_ARRAY
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE
     &       , SAVE, TARGET :: BUF_MAX_ARRAY
        PUBLIC :: BUF_LMAX_ARRAY, BUF_MAX_ARRAY
      CONTAINS
        SUBROUTINE DMUMPS_BUF_TRY_FREE_CB()
        CALL DMUMPS_BUF_TRY_FREE(BUF_CB)
        RETURN
        END SUBROUTINE DMUMPS_BUF_TRY_FREE_CB
        SUBROUTINE DMUMPS_BUF_TRY_FREE(B)
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: B
        INCLUDE 'mpif.h'
        LOGICAL :: FLAG
        INTEGER :: IERR_MPI
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG,
     &                   STATUS, IERR_MPI )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END iF
        RETURN
        END SUBROUTINE DMUMPS_BUF_TRY_FREE
        SUBROUTINE DMUMPS_BUF_INI_MYID( MYID )
        IMPLICIT NONE
        INTEGER MYID
        BUF_MYID  = MYID
        RETURN
        END SUBROUTINE DMUMPS_BUF_INI_MYID
        SUBROUTINE DMUMPS_BUF_INIT( IntSize, RealSize )
        IMPLICIT NONE
        INTEGER IntSize, RealSize
        SIZEofINT = IntSize
        SIZEofREAL = RealSize
        NULLIFY(BUF_CB  %CONTENT)
        NULLIFY(BUF_SMALL%CONTENT)
        NULLIFY(BUF_LOAD%CONTENT)
        BUF_CB%LBUF     = 0
        BUF_CB%LBUF_INT = 0
        BUF_CB%HEAD     = 1
        BUF_CB%TAIL     = 1
        BUF_CB%ILASTMSG = 1
        BUF_SMALL%LBUF     = 0
        BUF_SMALL%LBUF_INT = 0
        BUF_SMALL%HEAD     = 1
        BUF_SMALL%TAIL     = 1
        BUF_SMALL%ILASTMSG = 1
        BUF_LOAD%LBUF     = 0
        BUF_LOAD%LBUF_INT = 0
        BUF_LOAD%HEAD     = 1
        BUF_LOAD%TAIL     = 1
        BUF_LOAD%ILASTMSG = 1
        RETURN
        END SUBROUTINE DMUMPS_BUF_INIT
        SUBROUTINE DMUMPS_BUF_ALLOC_CB( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL BUF_ALLOC( BUF_CB, SIZE, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_ALLOC_CB
        SUBROUTINE DMUMPS_BUF_ALLOC_SMALL_BUF( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL BUF_ALLOC( BUF_SMALL, SIZE, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_ALLOC_SMALL_BUF
        SUBROUTINE DMUMPS_BUF_ALLOC_LOAD_BUFFER( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL BUF_ALLOC( BUF_LOAD, SIZE, IERR )        
        RETURN
        END SUBROUTINE DMUMPS_BUF_ALLOC_LOAD_BUFFER
        SUBROUTINE DMUMPS_BUF_DEALL_LOAD_BUFFER( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL BUF_DEALL( BUF_LOAD, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_DEALL_LOAD_BUFFER
        SUBROUTINE DMUMPS_BUF_DEALL_MAX_ARRAY()
        IMPLICIT NONE
        IF (allocated( BUF_MAX_ARRAY)) DEALLOCATE( BUF_MAX_ARRAY )
        RETURN
        END SUBROUTINE DMUMPS_BUF_DEALL_MAX_ARRAY
        SUBROUTINE DMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
        IMPLICIT NONE
        INTEGER IERR, NFS4FATHER
        IERR = 0
        IF (allocated( BUF_MAX_ARRAY)) THEN
          IF (BUF_LMAX_ARRAY .GE. NFS4FATHER) RETURN
          DEALLOCATE( BUF_MAX_ARRAY )
        ENDIF
        ALLOCATE(BUF_MAX_ARRAY(NFS4FATHER),stat=IERR)
        IF ( IERR .GT. 0 ) THEN
           IERR = -1
           RETURN
        END IF
        BUF_LMAX_ARRAY=NFS4FATHER
        RETURN
        END SUBROUTINE DMUMPS_BUF_MAX_ARRAY_MINSIZE
        SUBROUTINE DMUMPS_BUF_DEALL_CB( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL BUF_DEALL( BUF_CB, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_DEALL_CB
        SUBROUTINE DMUMPS_BUF_DEALL_SMALL_BUF( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL BUF_DEALL( BUF_SMALL, IERR )
        RETURN
        END SUBROUTINE DMUMPS_BUF_DEALL_SMALL_BUF
        SUBROUTINE BUF_ALLOC( BUF, SIZE, IERR )
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: BUF
        INTEGER SIZE, IERR
        IERR         = 0
        BUF%LBUF     = SIZE
        BUF%LBUF_INT = ( SIZE + SIZEofINT - 1 ) / SIZEofINT
        IF ( associated ( BUF%CONTENT ) ) DEALLOCATE( BUF%CONTENT )
        ALLOCATE( BUF%CONTENT( BUF%LBUF_INT ), stat = IERR )
       IF (IERR .NE. 0) THEN
          NULLIFY( BUF%CONTENT )
          IERR         = -1
          BUF%LBUF     =  0
          BUF%LBUF_INT =  0
        END IF
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE BUF_ALLOC
        SUBROUTINE BUF_DEALL( BUF, IERR )
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: BUF
        INTEGER :: IERR
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        LOGICAL :: FLAG
        IF ( .NOT. associated ( BUF%CONTENT ) ) THEN
          BUF%HEAD     = 1
          BUF%LBUF     = 0
          BUF%LBUF_INT = 0
          BUF%TAIL     = 1
          BUF%ILASTMSG = 1
          RETURN
        END IF
        DO WHILE ( BUF%HEAD.NE.0 .AND. BUF%HEAD .NE. BUF%TAIL )
          CALL MPI_TEST(BUF%CONTENT( BUF%HEAD + REQ ), FLAG,
     &                  STATUS, IERR_MPI)
          IF ( .not. FLAG ) THEN
            WRITE(*,*) '** Warning: trying to cancel a request.'
            WRITE(*,*) '** This might be problematic'
            CALL MPI_CANCEL( BUF%CONTENT( BUF%HEAD + REQ ), IERR_MPI )
            CALL MPI_REQUEST_FREE( BUF%CONTENT( BUF%HEAD + REQ ),
     &                             IERR_MPI )
          END IF
          BUF%HEAD = BUF%CONTENT( BUF%HEAD + NEXT )
        END DO
        DEALLOCATE( BUF%CONTENT )
        NULLIFY( BUF%CONTENT )
        BUF%LBUF     = 0
        BUF%LBUF_INT = 0
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE BUF_DEALL
        SUBROUTINE DMUMPS_BUF_SEND_CB( NBROWS_ALREADY_SENT,
     &                                INODE, FPERE, NFRONT, LCONT,
     &                                NASS, NPIV,
     &                                IWROW, IWCOL, A, PACKED_CB,
     &                                DEST, TAG, COMM, KEEP, IERR )
        IMPLICIT NONE
        INTEGER DEST, TAG, COMM, IERR
        INTEGER NBROWS_ALREADY_SENT
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER INODE, FPERE, NFRONT, LCONT, NASS, NPIV 
        INTEGER IWROW( LCONT ), IWCOL( LCONT )
        DOUBLE PRECISION A( * )
        LOGICAL PACKED_CB
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER NBROWS_PACKET
        INTEGER POSITION, IREQ, IPOS, I, J1
        INTEGER SIZE1, SIZE2, SIZE_PACK, SIZE_AV, SIZE_AV_REALS
        INTEGER IZERO, IONE
        INTEGER SIZECB
        INTEGER LCONT_SENT
        INTEGER DEST2(1)
        PARAMETER( IZERO = 0, IONE = 1 )
        LOGICAL RECV_BUF_SMALLER_THAN_SEND
        DOUBLE PRECISION TMP
        DEST2(1) = DEST
        IERR = 0
        IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          CALL MPI_PACK_SIZE( 11 + LCONT + LCONT, MPI_INTEGER,
     &                        COMM, SIZE1,  IERR_MPI)
        ELSE
          CALL MPI_PACK_SIZE( 5, MPI_INTEGER, COMM, SIZE1, IERR_MPI)
        ENDIF
        CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
        IF ( SIZE_AV .LT. SIZE_RBUF_BYTES ) THEN
          RECV_BUF_SMALLER_THAN_SEND = .FALSE.
        ELSE
          SIZE_AV = SIZE_RBUF_BYTES
          RECV_BUF_SMALLER_THAN_SEND = .TRUE.
        ENDIF
        SIZE_AV_REALS = ( SIZE_AV - SIZE1 ) / SIZEofREAL
        IF (SIZE_AV_REALS < 0 ) THEN
          NBROWS_PACKET = 0
        ELSE
          IF (PACKED_CB) THEN
            TMP=2.0D0*dble(NBROWS_ALREADY_SENT)+1.0D0
            NBROWS_PACKET = int(
     &                      ( sqrt( TMP * TMP
     &                        + 8.0D0 * dble(SIZE_AV_REALS)) - TMP )
     &                        / 2.0D0 )
          ELSE
            IF (LCONT.EQ.0) THEN
              NBROWS_PACKET = 0
            ELSE
              NBROWS_PACKET = SIZE_AV_REALS / LCONT
            ENDIF
          ENDIF
        ENDIF
 10     CONTINUE
        NBROWS_PACKET = max(0,
     &            min(NBROWS_PACKET, LCONT - NBROWS_ALREADY_SENT))
        IF (NBROWS_PACKET .EQ. 0 .AND. LCONT .NE. 0) THEN
           IF (RECV_BUF_SMALLER_THAN_SEND) THEN
            IERR = -3
            GOTO 100
         ELSE
            IERR = -1
            GOTO 100
          ENDIF
        ENDIF
        IF (PACKED_CB) THEN
          SIZECB = (NBROWS_ALREADY_SENT*NBROWS_PACKET)+(NBROWS_PACKET
     &             *(NBROWS_PACKET+1))/2
        ELSE
          SIZECB = NBROWS_PACKET * LCONT
        ENDIF
        CALL MPI_PACK_SIZE( SIZECB, MPI_DOUBLE_PRECISION,
     &                    COMM, SIZE2,  IERR_MPI )
        SIZE_PACK = SIZE1 + SIZE2
        IF (SIZE_PACK .GT. SIZE_AV ) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF (NBROWS_PACKET > 0) THEN
             GOTO 10
          ELSE
             IF (RECV_BUF_SMALLER_THAN_SEND) THEN
               IERR=-3
               GOTO 100
            ELSE
               IERR = -1
               GOTO 100
             ENDIF
          ENDIF
        ENDIF
        IF (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.LCONT .AND.
     &     SIZE_PACK  .LT. SIZE_RBUF_BYTES / 4
     &    .AND. 
     &    .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &       THEN
            IERR = -1
            GOTO 100
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE , DEST2)
        IF (IERR .EQ. -1 .OR. IERR .EQ. -2) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF ( NBROWS_PACKET > 0 )  GOTO 10
        ENDIF
        IF ( IERR .LT. 0 ) GOTO 100
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        IF (PACKED_CB) THEN
          LCONT_SENT=-LCONT
        ELSE
          LCONT_SENT=LCONT
        ENDIF
        CALL MPI_PACK( LCONT_SENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NBROWS_PACKET, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        IF (NBROWS_ALREADY_SENT == 0) THEN
          CALL MPI_PACK( LCONT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( NASS-NPIV, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( LCONT , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IONE,  1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IWROW, LCONT, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( IWCOL, LCONT, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
        ENDIF
        IF ( LCONT .NE. 0 ) THEN
          J1 = 1 + NBROWS_ALREADY_SENT * NFRONT
          IF (PACKED_CB) THEN
           DO I = NBROWS_ALREADY_SENT+1,
     &            NBROWS_ALREADY_SENT+NBROWS_PACKET
            CALL MPI_PACK( A( J1 ), I, MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
             J1 = J1 + NFRONT
           END DO
          ELSE
           DO I = NBROWS_ALREADY_SENT+1,
     &            NBROWS_ALREADY_SENT+NBROWS_PACKET
            CALL MPI_PACK( A( J1 ), LCONT, MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
             J1 = J1 + NFRONT
           END DO
          ENDIF
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ),
     &                  IERR_MPI )
        IF ( SIZE_PACK .LT. POSITION ) THEN
          WRITE(*,*) 'Error Try_send_cb: SIZE, POSITION=',SIZE_PACK,
     &               POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &    CALL BUF_ADJUST( BUF_CB, POSITION )
        NBROWS_ALREADY_SENT = NBROWS_ALREADY_SENT + NBROWS_PACKET
        IF (NBROWS_ALREADY_SENT .NE. LCONT ) THEN
          IERR = -1
          RETURN
        ENDIF
 100    CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_CB
        SUBROUTINE DMUMPS_BUF_SEND_MASTER2SLAVE( NRHS, INODE, IFATH,
     &             EFF_CB_SIZE, LD_CB, LD_PIV, NPIV, 
     &             JBDEB, JBFIN,
     &             CB, SOL,
     &             DEST, COMM, KEEP, IERR )
        IMPLICIT NONE
        INTEGER NRHS, INODE, IFATH, EFF_CB_SIZE, LD_CB, LD_PIV, NPIV 
        INTEGER DEST, COMM, IERR, JBDEB, JBFIN
        DOUBLE PRECISION CB( LD_CB*(NRHS-1)+EFF_CB_SIZE )
        DOUBLE PRECISION SOL( max(1, LD_PIV*(NRHS-1)+NPIV) )
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE, SIZE1, SIZE2, K
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 6, MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( NRHS * (EFF_CB_SIZE + NPIV),
     &                      MPI_DOUBLE_PRECISION, COMM,
     &                      SIZE2, IERR_MPI )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( IFATH, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( EFF_CB_SIZE  , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NPIV , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBDEB , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBFIN , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        DO K = 1, NRHS
               CALL MPI_PACK( CB ( 1 + LD_CB * (K-1) ),
     &                        EFF_CB_SIZE, MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        END DO
        IF ( NPIV .GT. 0 ) THEN
          DO K=1, NRHS
          CALL MPI_PACK( SOL(1+LD_PIV*(K-1)),
     &                         NPIV, MPI_DOUBLE_PRECISION,
     &                         BUF_CB%CONTENT( IPOS ), SIZE,
     &                         POSITION, COMM, IERR_MPI )
          ENDDO
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, Master2Slave, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_send_master2slave: SIZE, POSITION = ',
     &               SIZE, POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_MASTER2SLAVE
        SUBROUTINE DMUMPS_BUF_SEND_VCB( NRHS_B, NODE1, NODE2, NCB, LDW,
     &             LONG,
     &             IW, W, JBDEB, JBFIN,
     &             RHSCOMP, NRHS, LRHSCOMP, IPOSINRHSCOMP, NPIV,
     &             KEEP,
     &             DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER LDW, DEST, TAG, COMM, IERR
        INTEGER NRHS_B, NODE1, NODE2, NCB, LONG, JBDEB, JBFIN
        INTEGER IW( max( 1, LONG ) )
        INTEGER, INTENT(IN) :: LRHSCOMP, NRHS, IPOSINRHSCOMP, NPIV
        DOUBLE PRECISION W( max( 1, LDW * NRHS_B ) )
        DOUBLE PRECISION RHSCOMP(LRHSCOMP,NRHS)
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS
        INTEGER SIZE1, SIZE2, SIZE, K
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        IF ( NODE2 .EQ. 0 ) THEN
         CALL MPI_PACK_SIZE( 4+LONG, MPI_INTEGER, COMM, SIZE1,
     &                       IERR_MPI )
        ELSE
         CALL MPI_PACK_SIZE( 6+LONG, MPI_INTEGER, COMM, SIZE1,
     &                       IERR_MPI )
        END IF
        SIZE2 = 0
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK_SIZE( NRHS_B*LONG, MPI_DOUBLE_PRECISION,
     &                        COMM, SIZE2, IERR_MPI )
        END IF
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( NODE1, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        IF ( NODE2 .NE. 0 ) THEN
          CALL MPI_PACK( NODE2, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( NCB, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        ENDIF
        CALL MPI_PACK( JBDEB, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBFIN, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LONG,  1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK( IW, LONG, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
          IF (NODE2.EQ.0) THEN
            DO K=1, NRHS_B
              IF (NPIV.GT.0) THEN
              CALL MPI_PACK( RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1), NPIV,
     &                          MPI_DOUBLE_PRECISION,
     &                          BUF_CB%CONTENT( IPOS ), SIZE,
     &                          POSITION, COMM, IERR_MPI )
              ENDIF
              IF (LONG-NPIV .NE.0) THEN
                CALL MPI_PACK( W(NPIV+1+(K-1)*LDW), LONG-NPIV,
     &                          MPI_DOUBLE_PRECISION,
     &                          BUF_CB%CONTENT( IPOS ), SIZE,
     &                          POSITION, COMM, IERR_MPI )
              ENDIF
            END DO
          ELSE
            DO K=1, NRHS_B
              CALL MPI_PACK( W(1+(K-1)*LDW), LONG, MPI_DOUBLE_PRECISION,
     &                          BUF_CB%CONTENT( IPOS ), SIZE,
     &                          POSITION, COMM, IERR_MPI )
            END DO
          ENDIF
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ),
     &                  IERR_MPI )
        IF ( SIZE .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_VCB
        SUBROUTINE DMUMPS_BUF_SEND_1INT( I, DEST, TAG, COMM,
     &                                   KEEP, IERR )
        IMPLICIT NONE
        INTEGER I
        INTEGER DEST, TAG, COMM, IERR
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER IPOS, IREQ, MSG_SIZE, POSITION
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 1, MPI_INTEGER,
     &                      COMM, MSG_SIZE, IERR_MPI )
        CALL BUF_LOOK( BUF_SMALL, IPOS, IREQ, MSG_SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
         write(6,*) ' Internal error in DMUMPS_BUF_SEND_1INT',
     &       ' Buf size (bytes)= ',BUF_SMALL%LBUF
      RETURN
        ENDIF
        POSITION=0
        CALL MPI_PACK( I, 1,
     &                 MPI_INTEGER, BUF_SMALL%CONTENT( IPOS ),
     &                 MSG_SIZE,
     &                 POSITION, COMM, IERR_MPI )
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_SMALL%CONTENT(IPOS), MSG_SIZE,
     &                  MPI_PACKED, DEST, TAG, COMM,
     &                  BUF_SMALL%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_1INT
        SUBROUTINE DMUMPS_BUF_ALL_EMPTY(CHECK_COMM_NODES,
     &             CHECK_COMM_LOAD,FLAG)
        LOGICAL, INTENT(IN)  :: CHECK_COMM_NODES, CHECK_COMM_LOAD
        LOGICAL, INTENT(OUT) :: FLAG
        LOGICAL FLAG1, FLAG2, FLAG3
        FLAG = .TRUE.
        IF (CHECK_COMM_NODES) THEN
          CALL DMUMPS_BUF_EMPTY( BUF_SMALL, FLAG1 )
          CALL DMUMPS_BUF_EMPTY( BUF_CB, FLAG2 )
          FLAG = FLAG .AND. FLAG1 .AND. FLAG2
        ENDIF
        IF ( CHECK_COMM_LOAD ) THEN
          CALL DMUMPS_BUF_EMPTY( BUF_LOAD, FLAG3 )
          FLAG = FLAG .AND. FLAG3
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_BUF_ALL_EMPTY
        SUBROUTINE DMUMPS_BUF_EMPTY( B, FLAG )
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: B
        LOGICAL :: FLAG
        INTEGER SIZE_AVAIL
        CALL DMUMPS_BUF_SIZE_AVAILABLE(B, SIZE_AVAIL)
        FLAG = ( B%HEAD == B%TAIL )
        RETURN
        END SUBROUTINE DMUMPS_BUF_EMPTY
        SUBROUTINE DMUMPS_BUF_SIZE_AVAILABLE( B, SIZE_AV )
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: B
        INTEGER SIZE_AV
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        LOGICAL :: FLAG
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS,
     &                   IERR_MPI )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END IF
        IF ( B%HEAD .LE. B%TAIL ) THEN
           SIZE_AV = max( B%LBUF_INT - B%TAIL, B%HEAD - 2 )
        ELSE
           SIZE_AV = B%HEAD - B%TAIL - 1
        END IF
        SIZE_AV = min(SIZE_AV - OVHSIZE, SIZE_AV)
        SIZE_AV = SIZE_AV * SIZEofINT
        RETURN
        END SUBROUTINE DMUMPS_BUF_SIZE_AVAILABLE
        SUBROUTINE DMUMPS_BUF_TEST()
        INTEGER :: IPOS, IREQ, IERR
        INTEGER, PARAMETER :: IONE=1
        INTEGER :: MSG_SIZE
        INTEGER :: DEST2(1)
        DEST2=-10
        MSG_SIZE=1
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, MSG_SIZE, IERR, 
     &                 IONE , DEST2,.TRUE.)
        RETURN
        END SUBROUTINE DMUMPS_BUF_TEST
        SUBROUTINE BUF_LOOK( B, IPOS, IREQ, MSG_SIZE, IERR, 
     &    NDEST , PDEST, TEST_ONLY)
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: B
        INTEGER, INTENT(IN)        :: MSG_SIZE
        INTEGER, INTENT(OUT)       :: IPOS, IREQ, IERR
        LOGICAL, INTENT(IN), OPTIONAL :: TEST_ONLY
        INTEGER NDEST
        INTEGER, INTENT(IN)        :: PDEST(max(1,NDEST))
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER :: MSG_SIZE_INT
        INTEGER :: IBUF
        LOGICAL :: FLAG
        INTEGER :: STATUS(MPI_STATUS_SIZE)
        IERR = 0
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS,
     &                   IERR_MPI )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END iF
        MSG_SIZE_INT = ( MSG_SIZE + ( SIZEofINT - 1 ) ) / SIZEofINT
        MSG_SIZE_INT = MSG_SIZE_INT + OVHSIZE
        IF (present(TEST_ONLY)) RETURN
        FLAG = (     ( B%HEAD .LE. B%TAIL )
     &               .AND. (
     &                 ( MSG_SIZE_INT .LE. B%LBUF_INT - B%TAIL )
     &                 .OR. ( MSG_SIZE_INT .LE. B%HEAD - 2 ) ) )
     &         .OR.
     &               ( ( B%HEAD .GT. B%TAIL )
     &               .AND. ( MSG_SIZE_INT .LE. B%HEAD - B%TAIL - 1 ) )
        IF ( .NOT. FLAG
     &       ) THEN
        IERR = -1
        IF ( MSG_SIZE_INT .GT. B%LBUF_INT - 1 ) THEN
           IERR = -2
        ENDIF
        IPOS = -1
        IREQ = -1
        RETURN
      END IF
        IF ( B%HEAD .LE. B%TAIL ) THEN
          IF ( MSG_SIZE_INT .LE. B%LBUF_INT - B%TAIL + 1 ) THEN
            IBUF = B%TAIL
          ELSE IF ( MSG_SIZE_INT .LE. B%HEAD - 1 ) THEN
            IBUF = 1
          END IF
        ELSE
          IBUF = B%TAIL
        END IF
        B%CONTENT( B%ILASTMSG + NEXT ) = IBUF
        B%ILASTMSG = IBUF
        B%TAIL = IBUF + MSG_SIZE_INT
        B%CONTENT( IBUF + NEXT ) = 0
        IPOS = IBUF + CONTENT
        IREQ = IBUF + REQ
        RETURN
        END SUBROUTINE BUF_LOOK
        SUBROUTINE BUF_ADJUST( BUF, SIZE )
        IMPLICIT NONE
        TYPE ( DMUMPS_COMM_BUFFER_TYPE ) :: BUF
        INTEGER SIZE
        INTEGER SIZE_INT
        SIZE_INT = ( SIZE + SIZEofINT - 1 ) / SIZEofINT
        SIZE_INT = SIZE_INT + OVHSIZE
        BUF%TAIL = BUF%ILASTMSG + SIZE_INT
        RETURN
        END SUBROUTINE BUF_ADJUST
      SUBROUTINE DMUMPS_BUF_SEND_DESC_BANDE(
     &             INODE, NBPROCFILS, NLIG, ILIG, NCOL, ICOL,
     &             NASS, NSLAVES, LIST_SLAVES, 
     &             ESTIM_NFS4FATHER_ATSON,
     &             DEST, IBC_SOURCE, NFRONT, COMM, KEEP, IERR
     &             , LRSTATUS
     &)
      IMPLICIT NONE
        INTEGER COMM, IERR, NFRONT
        INTEGER, intent(in) :: INODE
        INTEGER, intent(in) :: NLIG, NCOL, NASS, NSLAVES
        INTEGER, intent(in) :: ESTIM_NFS4FATHER_ATSON
        INTEGER NBPROCFILS, DEST
        INTEGER ILIG( NLIG )
        INTEGER ICOL( NCOL )
        INTEGER, INTENT(IN) :: IBC_SOURCE
        INTEGER LIST_SLAVES( NSLAVES )
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER, INTENT(IN) :: LRSTATUS
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE_INT, SIZE_BYTES, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE_INT = ( 9 + NLIG + NCOL + NSLAVES + 1 )
        SIZE_BYTES = SIZE_INT * SIZEofINT
        IF (SIZE_INT.GT.SIZE_RBUF_BYTES ) THEN
         IERR = -3
      RETURN
        END IF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_BYTES, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = SIZE_INT
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = INODE
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NBPROCFILS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NLIG
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NCOL
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NASS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NFRONT
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = LRSTATUS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = ESTIM_NFS4FATHER_ATSON
        POSITION = POSITION + 1
        IF (NSLAVES.GT.0) THEN
         BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 ) = 
     &   LIST_SLAVES( 1: NSLAVES )
         POSITION = POSITION + NSLAVES
        ENDIF
        BUF_CB%CONTENT( POSITION:POSITION + NLIG - 1 ) = ILIG
        POSITION = POSITION + NLIG
        BUF_CB%CONTENT( POSITION:POSITION + NCOL - 1 ) = ICOL
        POSITION = POSITION + NCOL
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE_BYTES ) THEN
          WRITE(*,*) 'Error in DMUMPS_BUF_SEND_DESC_BANDE :',
     &               ' wrong estimated size'
          CALL MUMPS_ABORT()
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE_BYTES,
     &                  MPI_PACKED,
     &                  DEST, MAITRE_DESC_BANDE, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_DESC_BANDE
        SUBROUTINE DMUMPS_BUF_SEND_MAITRE2( NBROWS_ALREADY_SENT,
     &  IPERE, ISON, NROW,
     &  IROW, NCOL, ICOL, VAL, LDA, NELIM, TYPE_SON,
     &  NSLAVES, SLAVES, DEST, COMM, IERR, 
     & 
     &  SLAVEF, KEEP,KEEP8, INIV2, TAB_POS_IN_PERE )
        IMPLICIT NONE
        INTEGER NBROWS_ALREADY_SENT
        INTEGER LDA, NELIM, TYPE_SON
        INTEGER IPERE, ISON, NROW, NCOL, NSLAVES
        INTEGER IROW( NROW )
        INTEGER ICOL( NCOL )
        INTEGER SLAVES( NSLAVES )
        DOUBLE PRECISION VAL(LDA, *)
        INTEGER IPOS, IREQ, DEST, COMM, IERR
        INTEGER SLAVEF, KEEP(500), INIV2
        INTEGER(8) KEEP8(150)
        INTEGER TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE1, SIZE2, SIZE3, SIZE_PACK, POSITION, I
        INTEGER NBROWS_PACKET, NCOL_SEND
        INTEGER SIZE_AV
        LOGICAL RECV_BUF_SMALLER_THAN_SEND
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        IF ( NELIM .NE. NROW ) THEN
          WRITE(*,*) 'Error in TRY_SEND_MAITRE2:',NELIM, NROW
          CALL MUMPS_ABORT()
        END IF
        IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          CALL MPI_PACK_SIZE( NROW+NCOL+7+NSLAVES, MPI_INTEGER,
     &                      COMM, SIZE1, IERR_MPI )
          IF ( TYPE_SON .eq. 2 ) THEN
          CALL MPI_PACK_SIZE( NSLAVES+1, MPI_INTEGER,
     &                          COMM, SIZE3, IERR_MPI )
          ELSE
            SIZE3 = 0
          ENDIF
          SIZE1=SIZE1+SIZE3
        ELSE
          CALL MPI_PACK_SIZE(7, MPI_INTEGER,COMM,SIZE1,IERR_MPI)
        ENDIF
        IF ( KEEP(50).ne.0  .AND. TYPE_SON .eq. 2 ) THEN
          NCOL_SEND = NROW
        ELSE
          NCOL_SEND = NCOL
        ENDIF
        CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
        IF (SIZE_AV .LT. SIZE_RBUF_BYTES) THEN
          RECV_BUF_SMALLER_THAN_SEND = .FALSE.
        ELSE
          RECV_BUF_SMALLER_THAN_SEND = .TRUE.
          SIZE_AV = SIZE_RBUF_BYTES
        ENDIF
        IF (NROW .GT. 0 ) THEN 
         NBROWS_PACKET = (SIZE_AV - SIZE1) / NCOL_SEND / SIZEofREAL
         NBROWS_PACKET = min(NBROWS_PACKET, NROW - NBROWS_ALREADY_SENT)
         NBROWS_PACKET = max(NBROWS_PACKET, 0)
        ELSE
          NBROWS_PACKET =0
        ENDIF
        IF (NBROWS_PACKET .EQ. 0 .AND. NROW .NE. 0) THEN
           IF (RECV_BUF_SMALLER_THAN_SEND) THEN
              IERR=-3
              GOTO 100
           ELSE
              IERR=-1
              GOTO 100
          ENDIF
        ENDIF
 10     CONTINUE
        CALL MPI_PACK_SIZE( NBROWS_PACKET * NCOL_SEND,
     &           MPI_DOUBLE_PRECISION,
     &           COMM, SIZE2, IERR_MPI )
        SIZE_PACK = SIZE1 + SIZE2
        IF (SIZE_PACK .GT. SIZE_AV) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF ( NBROWS_PACKET .GT. 0 ) THEN
            GOTO 10
          ELSE
            IF (RECV_BUF_SMALLER_THAN_SEND) THEN
                IERR = -3
                GOTO 100
             ELSE
                IERR = -1
                GOTO 100
            ENDIF
          ENDIF
        ENDIF
       IF (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.NROW .AND.
     &   SIZE_PACK - SIZE1  .LT. ( SIZE_RBUF_BYTES - SIZE1 ) / 2
     &  .AND. 
     &   .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &       THEN
           IERR = -1
           GOTO 100
       ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          GOTO 100
        ENDIF
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( ISON,  1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NCOL, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NBROWS_PACKET, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
        IF (NBROWS_ALREADY_SENT .EQ. 0) THEN
          IF (NSLAVES.GT.0) THEN
            CALL MPI_PACK( SLAVES, NSLAVES, MPI_INTEGER,
     &                BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                POSITION, COMM, IERR_MPI )
          ENDIF
          CALL MPI_PACK( IROW, NROW, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( ICOL, NCOL, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
          IF ( TYPE_SON .eq. 2 ) THEN
            CALL MPI_PACK( TAB_POS_IN_PERE(1,INIV2), NSLAVES+1, 
     &                 MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI )
          ENDIF
        ENDIF
        IF (NBROWS_PACKET.GE.1) THEN
          DO I=NBROWS_ALREADY_SENT+1,
     &                   NBROWS_ALREADY_SENT+NBROWS_PACKET
            CALL MPI_PACK( VAL(1,I), NCOL_SEND, 
     &               MPI_DOUBLE_PRECISION,
     &               BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &               POSITION, COMM, IERR_MPI )
          ENDDO
        ENDIF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, MAITRE2, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        IF ( SIZE_PACK .LT. POSITION ) THEN
          write(*,*) 'Try_send_maitre2, SIZE,POSITION=',
     &                SIZE_PACK,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &    CALL BUF_ADJUST( BUF_CB, POSITION )
        NBROWS_ALREADY_SENT = NBROWS_ALREADY_SENT + NBROWS_PACKET
        IF ( NBROWS_ALREADY_SENT .NE. NROW ) THEN
          IERR = -1
        ENDIF
 100    CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_MAITRE2
        SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE2(NBROWS_ALREADY_SENT,
     &  DESC_IN_LU,
     &  IPERE, NFRONT_PERE, NASS_PERE, NFS4FATHER,
     &  NSLAVES_PERE,
     &  ISON, NBROW, LMAP, MAPROW, PERM, IW_CBSON, A_CBSON, LA_CBSON,
     &  ISLAVE, PDEST, PDEST_MASTER, COMM, IERR, 
     &  
     & KEEP,KEEP8, STEP, N, SLAVEF,
     & ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     & PACKED_CB, KEEP253_LOC, NVSCHUR,
     & SON_NIV, MYID, NPIV_CHECK )
        USE DMUMPS_LR_TYPE
        USE DMUMPS_LR_DATA_M
        IMPLICIT NONE
        INTEGER NBROWS_ALREADY_SENT
        INTEGER, INTENT (in) :: KEEP253_LOC, NVSCHUR
        INTEGER, INTENT (in) :: SON_NIV
        INTEGER, INTENT (in), OPTIONAL :: NPIV_CHECK
        INTEGER IPERE, ISON, NBROW, MYID
        INTEGER PDEST, ISLAVE, COMM, IERR
        INTEGER PDEST_MASTER, NASS_PERE, NSLAVES_PERE,
     &       NFRONT_PERE, LMAP
        INTEGER MAPROW( LMAP ), PERM( max(1, NBROW ))
        INTEGER IW_CBSON( * )
        DOUBLE PRECISION A_CBSON( : )
        INTEGER(8) :: LA_CBSON
        LOGICAL DESC_IN_LU, PACKED_CB
       INTEGER   KEEP(500), N , SLAVEF
       INTEGER(8) KEEP8(150)
       INTEGER   STEP(N), 
     &          ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR_MPI
      INTEGER NFS4FATHER,SIZE3,PS1,NCA,LROW1
      INTEGER(8) :: ASIZE
      LOGICAL COMPUTE_MAX
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: M_ARRAY
      INTEGER NBROWS_PACKET 
      INTEGER MAX_ROW_LENGTH
      INTEGER LROW, NELIM
      INTEGER(8) :: ITMP8
      INTEGER NPIV, NFRONT, HS
      INTEGER SIZE_PACK, SIZE0, SIZE1, SIZE2, POSITION,I
      INTEGER SIZE_INTEGERS, B, SIZE_REALS, TMPSIZE, ONEorTWO, SIZE_AV
      INTEGER NBINT, L
      INTEGER(8) :: APOS, SHIFTCB_SON, LDA_SON8
      INTEGER IPOS_IN_SLAVE
      INTEGER STATE_SON
      INTEGER INDICE_PERE, NROW, IPOS, IREQ, NOSLA
      INTEGER IONE, J, THIS_ROW_LENGTH
      INTEGER SIZE_DESC_BANDE, DESC_BANDE_BYTES
      LOGICAL RECV_BUF_SMALLER_THAN_SEND
      LOGICAL NOT_ENOUGH_SPACE
      INTEGER PDEST2(1)
      LOGICAL CB_IS_LR
      TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:)
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_ROW, BEGS_BLR_COL,
     &                    BEGS_BLR_STA
      INTEGER :: NB_ROW_SHIFT, NB_COL_SHIFT, NASS_SHIFT, PANEL2SEND,
     &           CURRENT_PANEL_SIZE, NB_BLR_ROWS, NB_BLR_COLS,
     &           CB_IS_LR_INT, NCOL_SHIFT, NROW_SHIFT,
     &           NBROWS_PACKET_2PACK,
     &           PANEL_BEG_OFFSET
      INTEGER :: NPIV_LR
      PARAMETER ( IONE=1 )
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
      CB_IS_LR = (IW_CBSON(1+XXLR).EQ.1
     &       .OR. IW_CBSON(1+XXLR).EQ.3)
      IF (CB_IS_LR) THEN
        CB_IS_LR_INT = 1
      ELSE
        CB_IS_LR_INT = 0
      ENDIF
      COMPUTE_MAX = (KEEP(219) .NE. 0) .AND.
     &              (KEEP(50) .EQ. 2) .AND.
     &              (PDEST.EQ.PDEST_MASTER)
      IF (NBROWS_ALREADY_SENT == 0) THEN 
        IF (COMPUTE_MAX) THEN
          CALL DMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
          IF (IERR .NE. 0) THEN
            IERR         = -4
            RETURN
          ENDIF
        ENDIF
      ENDIF
      PDEST2(1) = PDEST
      IERR   = 0
      LROW   = IW_CBSON( 1 + KEEP(IXSZ))
      NELIM  = IW_CBSON( 2 + KEEP(IXSZ))
      NPIV   = IW_CBSON( 4 + KEEP(IXSZ))
      IF ( NPIV .LT. 0 ) THEN
          NPIV = 0
      END IF
      NROW   = IW_CBSON( 3 + KEEP(IXSZ))
      NFRONT = LROW + NPIV
      HS     = 6 + IW_CBSON( 6 + KEEP(IXSZ)) + KEEP(IXSZ)
      IF (CB_IS_LR) THEN
        CALL DMUMPS_BLR_RETRIEVE_CB_LRB(IW_CBSON(1+XXF), CB_LRB)
        IF (SON_NIV.EQ.1) THEN
          CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(IW_CBSON(1+XXF),
     &                     BEGS_BLR_ROW)
          CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_DYN(IW_CBSON(1+XXF),
     &                     BEGS_BLR_COL)
          NB_BLR_ROWS = size(BEGS_BLR_ROW) - 1
          CALL DMUMPS_BLR_RETRIEVE_NB_PANELS(IW_CBSON(1+XXF), 
     &                    NB_COL_SHIFT)
          NB_ROW_SHIFT = NB_COL_SHIFT
          NASS_SHIFT = BEGS_BLR_ROW(NB_ROW_SHIFT+1)-1
          NPIV_LR  = BEGS_BLR_COL(NB_COL_SHIFT+1)-1
        ELSE
          NPIV_LR=NPIV
          CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(IW_CBSON(1+XXF),
     &                     BEGS_BLR_STA)
          NB_BLR_ROWS = size(BEGS_BLR_STA) - 2
          BEGS_BLR_ROW => BEGS_BLR_STA(2:NB_BLR_ROWS+2)
          CALL DMUMPS_BLR_RETRIEVE_BEGS_BLR_C(IW_CBSON(1+XXF),
     &                     BEGS_BLR_COL, NB_COL_SHIFT)
          NASS_SHIFT = 0
          NB_ROW_SHIFT = 0
        ENDIF
        PANEL2SEND = -1
        DO I=NB_ROW_SHIFT+1,NB_BLR_ROWS
          IF (BEGS_BLR_ROW(I+1)-1-NASS_SHIFT
     &         .GT.NBROWS_ALREADY_SENT+PERM(1)-1) THEN
            PANEL2SEND = I
            EXIT
          ENDIF
        ENDDO
        IF (PANEL2SEND.EQ.-1) THEN
          write(*,*) 'Internal error: PANEL2SEND not found'
          CALL MUMPS_ABORT()
        ENDIF
        IF (KEEP(50).EQ.0) THEN
          NB_BLR_COLS = size(BEGS_BLR_COL)  - 1
        ELSEIF (SON_NIV.EQ.1) THEN
          NB_BLR_COLS = PANEL2SEND
        ELSE
          NB_BLR_COLS = -1
          NCOL_SHIFT = NPIV_LR
          NROW_SHIFT = LROW - NROW
          DO I=NB_COL_SHIFT+1,size(BEGS_BLR_COL)-1
            IF (BEGS_BLR_COL(I+1)-NCOL_SHIFT.GT.
     &          BEGS_BLR_ROW(PANEL2SEND+1)-1+NROW_SHIFT) THEN
              NB_BLR_COLS = I
              EXIT
            ENDIF
          ENDDO
          IF (NB_BLR_COLS.EQ.-1) THEN
            write(*,*) 'Internal error: NB_BLR_COLS not found'
            CALL MUMPS_ABORT()
          ENDIF
          MAX_ROW_LENGTH = BEGS_BLR_ROW(PANEL2SEND+1)-1+NROW_SHIFT
        ENDIF
        CURRENT_PANEL_SIZE = BEGS_BLR_ROW(PANEL2SEND+1)
     &                     - BEGS_BLR_ROW(PANEL2SEND) 
        PANEL_BEG_OFFSET = PERM(1) + NBROWS_ALREADY_SENT -
     &                     BEGS_BLR_ROW(PANEL2SEND) + NASS_SHIFT
      ENDIF
      STATE_SON = IW_CBSON(1+XXS)
      IF (STATE_SON .EQ. S_NOLCBCONTIG) THEN
               LDA_SON8    = int(LROW,8)
               SHIFTCB_SON = int(NPIV,8)*int(NROW,8)
      ELSE IF (STATE_SON .EQ. S_NOLCLEANED) THEN
               LDA_SON8    = int(LROW,8)
               SHIFTCB_SON = 0_8
      ELSE
               LDA_SON8     = int(NFRONT,8)
               SHIFTCB_SON = int(NPIV,8)
      ENDIF
      CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
      IF (PDEST .EQ. PDEST_MASTER) THEN
        SIZE_DESC_BANDE=0 
      ELSE
        SIZE_DESC_BANDE=(7+SLAVEF+KEEP(127)*2)
        SIZE_DESC_BANDE=SIZE_DESC_BANDE+int(dble(KEEP(12))*
     &                  dble(SIZE_DESC_BANDE)/100.0D0)
        SIZE_DESC_BANDE=max(SIZE_DESC_BANDE,
     &     7+NSLAVES_PERE+NFRONT_PERE+NFRONT_PERE-NASS_PERE)
      ENDIF
      DESC_BANDE_BYTES=SIZE_DESC_BANDE*SIZEofINT
      IF ( SIZE_AV .LT. SIZE_RBUF_BYTES-DESC_BANDE_BYTES ) THEN
        RECV_BUF_SMALLER_THAN_SEND = .FALSE.
      ELSE
        RECV_BUF_SMALLER_THAN_SEND = .TRUE.
        SIZE_AV = SIZE_RBUF_BYTES-DESC_BANDE_BYTES
      ENDIF
      SIZE1=0
      IF (NBROWS_ALREADY_SENT==0) THEN
          IF(COMPUTE_MAX) THEN
               CALL MPI_PACK_SIZE(1, MPI_INTEGER,
     &            COMM, SIZE0, IERR_MPI )
               IF(NFS4FATHER .GT. 0) THEN
                CALL MPI_PACK_SIZE( NFS4FATHER, MPI_DOUBLE_PRECISION,
     &             COMM, SIZE1, IERR_MPI )
               ENDIF
               SIZE1 = SIZE1+SIZE0
          ENDIF
      ENDIF
      IF (KEEP(50) .EQ. 0) THEN
        ONEorTWO = 1
      ELSE
        ONEorTWO = 2
      ENDIF
      IF (PDEST .EQ.PDEST_MASTER) THEN
        L = 0
      ELSE IF (KEEP(50) .EQ. 0) THEN
        L = LROW
      ELSE
        L = LROW + PERM(1) - LMAP + NBROWS_ALREADY_SENT - 1
        ONEorTWO=ONEorTWO+1
      ENDIF
      NBINT = 6 + L + 1 
      IF (CB_IS_LR) THEN
        NBINT = NBINT + 4*(NB_BLR_COLS-NB_COL_SHIFT) + 2
      ENDIF
      CALL MPI_PACK_SIZE( NBINT, MPI_INTEGER,
     &                    COMM, TMPSIZE, IERR_MPI )
      SIZE1 = SIZE1 + TMPSIZE
      SIZE_AV = SIZE_AV - SIZE1
      NOT_ENOUGH_SPACE=.FALSE.
      IF (SIZE_AV .LT.0 ) THEN
        NBROWS_PACKET = 0
        NOT_ENOUGH_SPACE=.TRUE.
      ELSE
        IF ( KEEP(50) .EQ. 0 ) THEN
          NBROWS_PACKET =
     &       SIZE_AV / ( ONEorTWO*SIZEofINT+LROW*SIZEofREAL)
        ELSE
          B = 2 * ONEorTWO + 
     &      ( 1 + 2 *  LROW + 2 * PERM(1) + 2 * NBROWS_ALREADY_SENT )
     &      * SIZEofREAL / SIZEofINT
          NBROWS_PACKET=int((dble(-B)+sqrt((dble(B)*dble(B))+
     &        dble(4)*dble(2)*dble(SIZE_AV)/dble(SIZEofINT) *
     &        dble(SIZEofREAL/SIZEofINT)))*
     &        dble(SIZEofINT) / dble(2) / dble(SIZEofREAL))
        ENDIF
      ENDIF
 10   CONTINUE
      NBROWS_PACKET = max( 0, NBROWS_PACKET)
      NBROWS_PACKET = min(NBROW-NBROWS_ALREADY_SENT, NBROWS_PACKET)
      NOT_ENOUGH_SPACE = NOT_ENOUGH_SPACE .OR.
     &                   (NBROWS_PACKET .EQ.0.AND. NBROW.NE.0)
      NBROWS_PACKET_2PACK = NBROWS_PACKET
      IF (CB_IS_LR) THEN
        NBROWS_PACKET_2PACK = CURRENT_PANEL_SIZE
        CALL MUMPS_BLR_GET_SIZEREALS_CB_LRB(SIZE_REALS, CB_LRB, 
     &             NB_ROW_SHIFT,
     &             NB_COL_SHIFT, NB_BLR_COLS, PANEL2SEND)
        NOT_ENOUGH_SPACE = (SIZE_AV.LT.SIZE_REALS)
        IF (.NOT.NOT_ENOUGH_SPACE) THEN
          NBROWS_PACKET = min(NBROWS_PACKET,
     &              CURRENT_PANEL_SIZE-PANEL_BEG_OFFSET)
        ENDIF
      ENDIF
      IF (NOT_ENOUGH_SPACE) THEN
         IF (RECV_BUF_SMALLER_THAN_SEND) THEN
          IERR = -3
          GOTO 100
       ELSE
          IERR = -1
          GOTO 100
        ENDIF
      ENDIF
      IF (CB_IS_LR) THEN
        IF (KEEP(50).EQ.0) THEN
          MAX_ROW_LENGTH = -99999
        ELSEIF (SON_NIV.EQ.1) THEN
          MAX_ROW_LENGTH = LROW+PERM(1)-LMAP+NBROWS_ALREADY_SENT
     &                 + NBROWS_PACKET_2PACK-1
        ENDIF
      ELSE
        IF (KEEP(50).EQ.0) THEN
          MAX_ROW_LENGTH = -99999
          SIZE_REALS = NBROWS_PACKET_2PACK * LROW
        ELSE
          SIZE_REALS = (  LROW + PERM(1) + NBROWS_ALREADY_SENT ) *
     &      NBROWS_PACKET_2PACK + ( NBROWS_PACKET_2PACK *
     &     ( NBROWS_PACKET_2PACK + 1) ) / 2
          MAX_ROW_LENGTH = LROW+PERM(1)-LMAP+NBROWS_ALREADY_SENT
     &                 + NBROWS_PACKET_2PACK-1
        ENDIF
      ENDIF
      SIZE_INTEGERS = ONEorTWO* NBROWS_PACKET_2PACK
      CALL MPI_PACK_SIZE( SIZE_REALS, MPI_DOUBLE_PRECISION,
     &                    COMM, SIZE2,  IERR_MPI )
      CALL MPI_PACK_SIZE( SIZE_INTEGERS, MPI_INTEGER,
     &                    COMM, SIZE3,  IERR_MPI )
      IF (SIZE2 + SIZE3 .GT. SIZE_AV ) THEN
         NBROWS_PACKET = NBROWS_PACKET -1
         IF (NBROWS_PACKET .GT. 0 .AND..NOT.CB_IS_LR) THEN
           GOTO 10
         ELSE
           IF (RECV_BUF_SMALLER_THAN_SEND) THEN
             IERR = -3
             GOTO 100
          ELSE
             IERR = -1
             GOTO 100
           ENDIF
         ENDIF
      ENDIF
        SIZE_PACK = SIZE1 + SIZE2 + SIZE3
        IF (NBROWS_PACKET + NBROWS_ALREADY_SENT.NE.NBROW .AND.
     &       SIZE_PACK  .LT. SIZE_RBUF_BYTES / 4 .AND.
     &    .NOT. RECV_BUF_SMALLER_THAN_SEND .AND.
     &    .NOT. CB_IS_LR)
     &       THEN
            IERR = -1
            GOTO 100
        ENDIF
        IF (SIZE_PACK.GT.SIZE_RBUF_BYTES ) THEN
          IERR = -3
          GOTO 100
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE , PDEST2)
        IF (IERR .EQ. -1 .OR. IERR.EQ. -2) THEN
          NBROWS_PACKET = NBROWS_PACKET - 1
          IF (NBROWS_PACKET > 0 ) GOTO 10
        ENDIF
        IF ( IERR .LT. 0 ) GOTO 100
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        CALL MPI_PACK( NBROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        IF (KEEP(50)==0) THEN
          CALL MPI_PACK( LROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        ELSE
          CALL MPI_PACK( MAX_ROW_LENGTH, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        ENDIF
        CALL MPI_PACK( NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        CALL MPI_PACK( NBROWS_PACKET, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
        CALL MPI_PACK( CB_IS_LR_INT, 1, MPI_INTEGER,
     &          BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &          POSITION, COMM, IERR_MPI  )
        IF ( PDEST .NE. PDEST_MASTER ) THEN
          IF (KEEP(50)==0) THEN
          CALL MPI_PACK( IW_CBSON( HS + NROW +  NPIV + 1 ), LROW,
     &                 MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
          ELSE
           IF (MAX_ROW_LENGTH > 0) THEN
           CALL MPI_PACK( IW_CBSON( HS + NROW +  NPIV + 1 ),
     &                 MAX_ROW_LENGTH,
     &                 MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                 POSITION, COMM, IERR_MPI  )
           ENDIF
          ENDIF
        END IF
        DO J=NBROWS_ALREADY_SENT+1,NBROWS_ALREADY_SENT+NBROWS_PACKET
           I = PERM(J)
           INDICE_PERE=MAPROW(I)
           CALL MUMPS_BLOC2_GET_ISLAVE(
     &          KEEP,KEEP8, IPERE, STEP, N, SLAVEF,
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &          NASS_PERE,
     &          NFRONT_PERE - NASS_PERE,
     &          NSLAVES_PERE,
     &          INDICE_PERE,
     &          NOSLA,
     &          IPOS_IN_SLAVE )
           INDICE_PERE = IPOS_IN_SLAVE
           CALL MPI_PACK( INDICE_PERE, 1, MPI_INTEGER,
     &          BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &          POSITION, COMM, IERR_MPI )
        ENDDO
        IF (CB_IS_LR) THEN
          CALL DMUMPS_BLR_PACK_CB_LRB(CB_LRB, NB_ROW_SHIFT,
     &             NB_COL_SHIFT, NB_BLR_COLS, PANEL2SEND,
     &             PANEL_BEG_OFFSET,
     &             BUF_CB%CONTENT(IPOS:),
     &             SIZE_PACK, POSITION, COMM, IERR)
          IF (KEEP(50).ne.0) THEN
            DO J=NBROWS_ALREADY_SENT+1,NBROWS_ALREADY_SENT+NBROWS_PACKET
              I = PERM(J)
              THIS_ROW_LENGTH = LROW + I - LMAP
              CALL MPI_PACK( THIS_ROW_LENGTH, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                        POSITION, COMM, IERR_MPI )
            ENDDO
          ENDIF
          GOTO 200
        ENDIF
        DO J=NBROWS_ALREADY_SENT+1,NBROWS_ALREADY_SENT+NBROWS_PACKET
           I = PERM(J)
           INDICE_PERE=MAPROW(I)
           CALL MUMPS_BLOC2_GET_ISLAVE(
     &          KEEP,KEEP8, IPERE, STEP, N, SLAVEF,
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &          
     &          NASS_PERE,
     &          NFRONT_PERE - NASS_PERE,
     &          NSLAVES_PERE,
     &          INDICE_PERE,
     &          NOSLA,
     &          IPOS_IN_SLAVE )
          IF (KEEP(50).ne.0) THEN
            THIS_ROW_LENGTH = LROW + I - LMAP
            CALL MPI_PACK( THIS_ROW_LENGTH, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &                      POSITION, COMM, IERR_MPI )
         ELSE
            THIS_ROW_LENGTH = LROW
         ENDIF
         IF (DESC_IN_LU) THEN 
            IF ( PACKED_CB ) THEN
             IF (NELIM.EQ.0) THEN
               ITMP8 = int(I,8)
             ELSE
               ITMP8 = int(NELIM+I,8)
             ENDIF
             APOS = ITMP8 * (ITMP8-1_8) / 2_8 + 1_8
            ELSE
             APOS = int(I+NELIM-1, 8) * int(LROW,8) + 1_8
            ENDIF
         ELSE
            IF ( PACKED_CB ) THEN
             IF ( LROW .EQ. NROW )  THEN
               ITMP8 = int(I,8)
               APOS  = ITMP8 * (ITMP8-1_8)/2_8 + 1_8
             ELSE
               ITMP8 = int(I + LROW - NROW,8)
               APOS  = ITMP8 * (ITMP8-1_8)/2_8 + 1_8 -
     &                 int(LROW - NROW, 8) * int(LROW-NROW+1,8) / 2_8
             ENDIF
            ELSE
             APOS = int( I - 1, 8 ) * LDA_SON8 + SHIFTCB_SON + 1_8
            ENDIF
         ENDIF
         CALL MPI_PACK( A_CBSON( APOS ), THIS_ROW_LENGTH,
     &        MPI_DOUBLE_PRECISION,
     &        BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &        POSITION, COMM, IERR_MPI )
        ENDDO
 200    CONTINUE
      IF (NBROWS_ALREADY_SENT == 0) THEN
        IF (COMPUTE_MAX) THEN
           CALL MPI_PACK(NFS4FATHER,1,
     &          MPI_INTEGER,
     &          BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &          POSITION, COMM, IERR_MPI )
           IF (NFS4FATHER .GT. 0) THEN
            IF (CB_IS_LR) THEN
              CALL DMUMPS_BLR_RETRIEVE_M_ARRAY (
     &            IW_CBSON(1+XXF), M_ARRAY)
              CALL MPI_PACK(M_ARRAY(1), NFS4FATHER,
     &             MPI_DOUBLE_PRECISION,
     &             BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &             POSITION, COMM, IERR_MPI )
              CALL DMUMPS_BLR_FREE_M_ARRAY ( IW_CBSON(1+XXF) )
            ELSE
              BUF_MAX_ARRAY(1:NFS4FATHER) = ZERO
              IF(MAPROW(NROW) .GT. NASS_PERE) THEN
                 DO PS1=1,NROW
                    IF(MAPROW(PS1).GT.NASS_PERE) EXIT
                 ENDDO
                 IF (DESC_IN_LU) THEN
                   IF (PACKED_CB) THEN
                    APOS = int(NELIM+PS1,8) * int(NELIM+PS1-1,8) /
     &                     2_8 + 1_8
                    NCA  = -44444
                    ASIZE  = int(NROW,8) * int(NROW+1,8)/2_8 -
     &                       int(NELIM+PS1,8) * int(NELIM+PS1-1,8)/2_8
                    LROW1  = PS1 + NELIM
                   ELSE
                    APOS = int(PS1+NELIM-1,8) * int(LROW,8) + 1_8
                    NCA = LROW
                    ASIZE = int(NCA,8) * int(NROW-PS1+1,8)
                    LROW1 = LROW
                   ENDIF
                 ELSE
                    IF (PACKED_CB) THEN
                      IF (NPIV.NE.0) THEN
         WRITE(*,*) "Error in PARPIV/DMUMPS_BUF_SEND_CONTRIB_TYPE2"
         CALL MUMPS_ABORT()
                      ENDIF
                      LROW1=LROW-NROW+PS1
                      ITMP8 = int(PS1 + LROW - NROW,8)
                      APOS = ITMP8 * (ITMP8 - 1_8) / 2_8 + 1_8 -
     &                       int(LROW-NROW,8)*int(LROW-NROW+1,8)/2_8
                      ASIZE = int(LROW,8)*int(LROW+1,8)/2_8 -
     &                       ITMP8*(ITMP8-1_8)/2_8
                      NCA   = -555555
                    ELSE
                      APOS = int(PS1-1,8) * LDA_SON8 + 1_8 + SHIFTCB_SON
                      NCA = int(LDA_SON8)
                      ASIZE = LA_CBSON - APOS + 1_8
                      LROW1=-666666
                    ENDIF
                 ENDIF
                 IF ( NROW-PS1+1-KEEP253_LOC-NVSCHUR .GT. 0 ) THEN
                   CALL DMUMPS_COMPUTE_MAXPERCOL(
     &                A_CBSON(APOS),ASIZE,NCA,
     &                NROW-PS1+1-KEEP253_LOC-NVSCHUR,
     &                BUF_MAX_ARRAY,NFS4FATHER,PACKED_CB,LROW1)
                 ENDIF
              ENDIF
              CALL MPI_PACK(BUF_MAX_ARRAY, NFS4FATHER,
     &             MPI_DOUBLE_PRECISION,
     &             BUF_CB%CONTENT( IPOS ), SIZE_PACK,
     &             POSITION, COMM, IERR_MPI )
            ENDIF 
           ENDIF
        ENDIF 
      ENDIF  
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  PDEST, CONTRIB_TYPE2, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        IF ( SIZE_PACK.LT. POSITION ) THEN
          WRITE(*,*) ' contniv2: SIZE, POSITION =',SIZE_PACK, POSITION
          WRITE(*,*) ' NBROW, LROW = ', NBROW, LROW
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_CB, POSITION )
        NBROWS_ALREADY_SENT=NBROWS_ALREADY_SENT + NBROWS_PACKET
        IF (NBROWS_ALREADY_SENT .NE. NBROW ) THEN
           IERR = -1
        ENDIF
 100    CONTINUE
        RETURN
      END SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE2
      SUBROUTINE MUMPS_BLR_GET_SIZEREALS_CB_LRB(SIZE_OUT,
     &        CB_LRB, NB_ROW_SHIFT, NB_COL_SHIFT, NB_BLR_COLS, 
     &        PANEL2SEND)
        USE DMUMPS_LR_TYPE
        IMPLICIT NONE
      TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:)
      INTEGER, INTENT(IN) :: NB_ROW_SHIFT, NB_COL_SHIFT, NB_BLR_COLS, 
     &                       PANEL2SEND
      INTEGER, intent(out) :: SIZE_OUT
      INTEGER :: J
      TYPE(LRB_TYPE), POINTER :: LRB
          SIZE_OUT = 0
          DO J=1,NB_BLR_COLS-NB_COL_SHIFT
            LRB => CB_LRB(PANEL2SEND-NB_ROW_SHIFT,J)
            IF (LRB%ISLR) THEN
              SIZE_OUT = SIZE_OUT + LRB%K*(LRB%M+LRB%N)
            ELSE
              SIZE_OUT = SIZE_OUT + LRB%M*LRB%N
            ENDIF
          ENDDO
          RETURN
      END SUBROUTINE MUMPS_BLR_GET_SIZEREALS_CB_LRB
      SUBROUTINE DMUMPS_BLR_PACK_CB_LRB(
     &        CB_LRB, NB_ROW_SHIFT, NB_COL_SHIFT, NB_BLR_COLS, 
     &        PANEL2SEND, PANEL_BEG_OFFSET,
     &        BUF, LBUF, POSITION, COMM, IERR)
        USE DMUMPS_LR_TYPE
        IMPLICIT NONE
      TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:)
      INTEGER, INTENT(IN) :: NB_ROW_SHIFT, NB_COL_SHIFT, NB_BLR_COLS, 
     &                       PANEL2SEND, PANEL_BEG_OFFSET
      INTEGER, intent(out) :: IERR
      INTEGER, intent(in)  :: COMM, LBUF  
      INTEGER, intent(inout) :: POSITION
      INTEGER, intent(inout) :: BUF(:) 
      INTEGER :: J, IERR_MPI
      INCLUDE 'mpif.h'
          IERR = 0
          CALL MPI_PACK( NB_BLR_COLS-NB_COL_SHIFT, 1, MPI_INTEGER,
     &              BUF(1), LBUF, POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( PANEL_BEG_OFFSET, 1, MPI_INTEGER,
     &              BUF(1), LBUF, POSITION, COMM, IERR_MPI )
          DO J=1,NB_BLR_COLS-NB_COL_SHIFT
            CALL DMUMPS_MPI_PACK_LRB(
     &               CB_LRB(PANEL2SEND-NB_ROW_SHIFT,J),
     &               BUF, LBUF, POSITION, COMM, IERR )
          ENDDO
        END SUBROUTINE DMUMPS_BLR_PACK_CB_LRB
        SUBROUTINE DMUMPS_BUF_SEND_MAPLIG( 
     &                INODE, NFRONT, NASS1, NFS4FATHER,
     &                ISON, MYID, NSLAVES, SLAVES_PERE,
     &                TROW, NCBSON,
     &                COMM, IERR,
     &                DEST, NDEST, SLAVEF, 
     & 
     &                KEEP,KEEP8, STEP, N, 
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &
     &                                  )
        IMPLICIT NONE
      INTEGER INODE, NFRONT, NASS1, NCBSON, NSLAVES, 
     &          NDEST
      INTEGER SLAVEF, MYID, ISON
      INTEGER TROW( NCBSON )
      INTEGER DEST( NDEST )
      INTEGER SLAVES_PERE( NSLAVES )
      INTEGER COMM, IERR
      INTEGER KEEP(500), N
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N), 
     &        ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE_AV, IDEST, NSEND, SIZE, NFS4FATHER
        INTEGER TROW_SIZE, POSITION, INDX, INIV2
        INTEGER IPOS, IREQ
        INTEGER IONE
        PARAMETER ( IONE=1 )
        INTEGER NASS_SON
        NASS_SON = -99998
        IERR = 0
        IF ( NDEST .eq. 1 ) THEN
          IF ( DEST(1).EQ.MYID )  GOTO 500
          SIZE = SIZEofINT * ( 7 + NSLAVES + NCBSON )
          IF ( NSLAVES.GT.0 ) THEN
             SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
          ENDIF
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
            IERR = -3
            RETURN
          END IF
          CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST
     &                 )
          IF (IERR .LT. 0 ) THEN
            RETURN
          ENDIF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NCBSON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFS4FATHER
              POSITION = POSITION + 1
              IF ( NSLAVES.GT.0 ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     &          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     &          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+NCBSON-1 ) =
     &        TROW( 1: NCBSON )
              POSITION = POSITION + NCBSON
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
                WRITE(*,*) 'Error in DMUMPS_BUF_SEND_MAPLIG :',
     &                     ' wrong estimated size'
                CALL MUMPS_ABORT()
              END IF
              KEEP(266)=KEEP(266)+1
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     &                        MPI_PACKED,
     &                        DEST( NDEST ), MAPLIG, COMM,
     &                        BUF_CB%CONTENT( IREQ ),
     &                        IERR_MPI )
        ELSE
          NSEND = 0
          DO IDEST = 1, NDEST
            IF ( DEST( IDEST ) .ne. MYID ) NSEND = NSEND + 1
          END DO
          SIZE = SIZEofINT * 
     &         ( ( OVHSIZE + 7 + NSLAVES )* NSEND + NCBSON )
          IF ( NSLAVES.GT.0 ) THEN
           SIZE = SIZE + SIZEofINT * NSEND*( NSLAVES + 1 )
          ENDIF
          CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
          IF ( SIZE_AV .LT. SIZE ) THEN
            IERR = -1
            RETURN
          END IF
          DO IDEST= 1, NDEST
            CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                IDEST, NCBSON, 
     &                NDEST, 
     &                TROW_SIZE, INDX  )
            SIZE = SIZEofINT * ( NSLAVES + TROW_SIZE + 7 )
            IF ( NSLAVES.GT.0 ) THEN
             SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
            ENDIF
            IF ( MYID .NE. DEST( IDEST ) ) THEN
               IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
                IERR = -3
                RETURN
              ENDIF
              CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR,
     &                       IONE, DEST(IDEST) )
              IF ( IERR .LT. 0 )  THEN
                WRITE(*,*) 'Internal error DMUMPS_BUF_SEND_MAPLIG',
     &                     'IERR after BUF_LOOK=',IERR
                CALL MUMPS_ABORT()
              END IF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = TROW_SIZE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFS4FATHER
              POSITION = POSITION + 1
              IF ( NSLAVES.GT.0 ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     &          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     &          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+TROW_SIZE-1 ) =
     &        TROW( INDX: INDX + TROW_SIZE - 1 )
              POSITION = POSITION + TROW_SIZE
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
               WRITE(*,*) ' ERROR 1 in TRY_SEND_MAPLIG:',
     &          'Wrong estimated size'
               CALL MUMPS_ABORT()
              END IF
              KEEP(266)=KEEP(266)+1
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     &                        MPI_PACKED,
     &                        DEST( IDEST ), MAPLIG, COMM,
     &                        BUF_CB%CONTENT( IREQ ),
     &                        IERR_MPI )
            END IF
          END DO
        END IF
 500    CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_MAPLIG
        SUBROUTINE DMUMPS_BUF_SEND_BLOCFACTO( INODE, NFRONT,
     &             NCOL, NPIV, FPERE, LASTBL, IPIV, VAL,
     &             PDEST, NDEST, KEEP, NB_BLOC_FAC,
     &             NSLAVES_TOT,
     &             WIDTH, COMM,
     &             NELIM, NPARTSASS, CURRENT_BLR_PANEL, 
     &             LR_ACTIVATED, BLR_LorU,
     &
     &             IERR )
      USE DMUMPS_LR_TYPE
      IMPLICIT NONE
        INTEGER, intent(in) :: INODE, NCOL, NPIV, 
     &                         FPERE, NFRONT, NDEST
        INTEGER, intent(in) :: IPIV( NPIV )
        DOUBLE PRECISION, intent(in) :: VAL( NFRONT, * )
        INTEGER, intent(in) :: PDEST( NDEST ) 
        INTEGER, intent(inout) :: KEEP(500)
        INTEGER, intent(in) :: NB_BLOC_FAC,
     &                         NSLAVES_TOT, COMM, WIDTH
        LOGICAL, intent(in) :: LASTBL
        LOGICAL, intent(in) :: LR_ACTIVATED
        INTEGER, intent(in) :: NELIM, NPARTSASS, CURRENT_BLR_PANEL
        TYPE (LRB_TYPE), DIMENSION(:), intent(in) :: BLR_LorU
        INTEGER, intent(inout) :: IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE3, SIZET,
     &          IDEST, IPOSMSG, I
        INTEGER NPIVSENT
        INTEGER SSS
        INTEGER  :: NBMSGS
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  RELAY_INFO
        INTEGER :: LRELAY_INFO, DEST_BLOCFACTO, TAG_BLOCFACTO
        INTEGER :: LR_ACTIVATED_INT
        IERR = 0
        LRELAY_INFO = 0
        NBMSGS = NDEST
        IF ( LASTBL ) THEN
          IF ( KEEP(50) .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 4 + NPIV + ( NBMSGS - 1 ) * OVHSIZE +
     &                          1+LRELAY_INFO,
     &                          MPI_INTEGER, COMM, SIZE1, IERR_MPI )
          ELSE
            CALL MPI_PACK_SIZE( 6 + NPIV + ( NBMSGS - 1 ) * OVHSIZE + 
     &                          1+LRELAY_INFO,
     &                          MPI_INTEGER, COMM, SIZE1, IERR_MPI )
          END IF
        ELSE
          IF ( KEEP(50) .eq. 0 ) THEN
          CALL MPI_PACK_SIZE( 3 + NPIV + ( NBMSGS - 1 ) * OVHSIZE + 
     &                        1+LRELAY_INFO,
     &                        MPI_INTEGER, COMM, SIZE1, IERR_MPI )
          ELSE
            CALL MPI_PACK_SIZE( 4 + NPIV + ( NBMSGS - 1 ) * OVHSIZE + 
     &                          1+LRELAY_INFO,
     &                          MPI_INTEGER, COMM, SIZE1, IERR_MPI )
          END IF
        END IF
        SIZE2 = 0
        CALL MPI_PACK_SIZE( 4, MPI_INTEGER, COMM, SIZE3, IERR_MPI )
        SIZE2=SIZE2+SIZE3
        IF ( KEEP(50).NE.0 ) THEN
          CALL MPI_PACK_SIZE( 1, MPI_INTEGER, COMM, SIZE3, IERR_MPI )
          SIZE2=SIZE2+SIZE3
        ENDIF
        IF ((NPIV.GT.0)
     &     ) THEN
          IF (.NOT. LR_ACTIVATED) THEN
            CALL MPI_PACK_SIZE( NPIV*NCOL, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE3, IERR_MPI )
            SIZE2 = SIZE2+SIZE3
          ELSE
            CALL MPI_PACK_SIZE( NPIV*(NPIV+NELIM), MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE3, IERR_MPI )
            SIZE2 = SIZE2+SIZE3
              CALL MUMPS_MPI_PACK_SIZE_LR( BLR_LorU, SIZE3, COMM, IERR )
            SIZE2 = SIZE2+SIZE3
          ENDIF
        ENDIF
        SIZET = SIZE1 + SIZE2 
        IF (SIZET.GT.SIZE_RBUF_BYTES) THEN
          SSS = 0 
          IF ( LASTBL ) THEN
           IF ( KEEP(50) .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 4 + NPIV + 1+LRELAY_INFO,
     &                        MPI_INTEGER, COMM, SSS, IERR_MPI )
           ELSE
            CALL MPI_PACK_SIZE( 6 + NPIV + 1+LRELAY_INFO, 
     &                           MPI_INTEGER, COMM, SSS, IERR_MPI )
           END IF
          ELSE
           IF ( KEEP(50) .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 3 + NPIV + 1+LRELAY_INFO,
     &                        MPI_INTEGER, COMM, SSS, IERR_MPI )
           ELSE
            CALL MPI_PACK_SIZE( 4 + NPIV + 1+LRELAY_INFO,
     &                        MPI_INTEGER, COMM, SSS, IERR_MPI )
           END IF
          END IF
          SSS = SSS + SIZE2  
          IF (SSS.GT.SIZE_RBUF_BYTES) THEN
           IERR = -3 
            RETURN
          ENDIF
        ENDIF
        IF (LRELAY_INFO.GT.0) THEN
         CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZET, IERR, 
     &                 NBMSGS , RELAY_INFO(2))
        ELSE
         CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZET, IERR, 
     &                 NBMSGS , PDEST)
        ENDIF
        IF ( IERR .LT. 0 ) THEN
          RETURN
        ENDIF
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NBMSGS - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NBMSGS - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NBMSGS - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NBMSGS
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        NPIVSENT = NPIV
        IF (LASTBL) NPIVSENT = -NPIV
        CALL MPI_PACK( NPIVSENT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        IF ( LASTBL .or. KEEP(50).ne.0 ) THEN
          CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     &                   BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                   POSITION, COMM, IERR_MPI )
        END IF
        IF ( LASTBL .AND. KEEP(50) .NE. 0 ) THEN
            CALL MPI_PACK( NSLAVES_TOT, 1, MPI_INTEGER,
     &                   BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                   POSITION, COMM, IERR_MPI )
            CALL MPI_PACK( NB_BLOC_FAC, 1, MPI_INTEGER,
     &                   BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                   POSITION, COMM, IERR_MPI )
        END IF
        CALL MPI_PACK( NCOL, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NELIM, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NPARTSASS, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( CURRENT_BLR_PANEL, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        IF (LR_ACTIVATED) THEN
          LR_ACTIVATED_INT = 1
        ELSE
          LR_ACTIVATED_INT = 0
        ENDIF
        CALL MPI_PACK( LR_ACTIVATED_INT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        IF ( KEEP(50) .ne. 0 ) THEN
          CALL MPI_PACK( NSLAVES_TOT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        ENDIF
        IF ( (NPIV.GT.0)
     &     ) THEN
          IF (NPIV.GT.0) THEN
            CALL MPI_PACK( IPIV, NPIV, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
          ENDIF
          IF (LR_ACTIVATED) THEN
              DO I = 1, NPIV
              CALL MPI_PACK( VAL(1,I), NPIV+NELIM,
     &                        MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
              END DO
              CALL DMUMPS_MPI_PACK_LR( BLR_LorU,
     &         BUF_CB%CONTENT(IPOSMSG:
     &              IPOSMSG+(SIZET+KEEP(34)-1)/KEEP(34)-1),
     &         SIZET, POSITION, COMM, IERR) 
          ELSE
            DO I = 1, NPIV
              CALL MPI_PACK( VAL(1,I), NCOL,
     &                        MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
            END DO
          ENDIF
        ENDIF
        CALL MPI_PACK( LRELAY_INFO, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        IF ( LRELAY_INFO.GT.0) 
     &    CALL MPI_PACK( RELAY_INFO, LRELAY_INFO, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        DO IDEST = 1, NBMSGS
          IF (LRELAY_INFO .GT. 0) THEN
            DEST_BLOCFACTO = RELAY_INFO(IDEST+1)
          ELSE
            DEST_BLOCFACTO = PDEST(IDEST)
          ENDIF
          IF ( KEEP(50) .EQ. 0) THEN
            TAG_BLOCFACTO = BLOC_FACTO
            KEEP(266)=KEEP(266)+1
            CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, 
     &                MPI_PACKED,
     &                DEST_BLOCFACTO, TAG_BLOCFACTO, COMM,
     &                BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     &                IERR_MPI )
          ELSE
            KEEP(266)=KEEP(266)+1
            CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, 
     &                MPI_PACKED,
     &                DEST_BLOCFACTO, BLOC_FACTO_SYM, COMM,
     &                BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     &                IERR_MPI )
          END IF
        END DO
        SIZET = SIZET - ( NBMSGS - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZET .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blocfacto : size < position'
          WRITE(*,*) ' Size,position=',SIZET,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZET .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_BLOCFACTO
        SUBROUTINE DMUMPS_BUF_SEND_BLFAC_SLAVE( INODE,
     &             NPIV, FPERE, IPOSK, JPOSK, UIP21K, NCOLU,
     &             NDEST, PDEST, COMM, KEEP,
     &             LR_ACTIVATED, BLR_LS, IPANEL,
     &             A , LA, POSBLOCFACTO, LD_BLOCFACTO,
     &             IPIV, MAXI_CLUSTER, IERR )
      USE DMUMPS_LR_TYPE
        IMPLICIT NONE
        INTEGER INODE, NCOLU, IPOSK, JPOSK, NPIV, NDEST, FPERE
        DOUBLE PRECISION UIP21K( NPIV, * )
        INTEGER PDEST( NDEST ) 
        INTEGER   COMM, IERR
        INTEGER, INTENT(INOUT) :: KEEP(500)
        LOGICAL, intent(in) :: LR_ACTIVATED
        TYPE (LRB_TYPE), DIMENSION(:), POINTER :: BLR_LS
        INTEGER(8), intent(in)  :: LA, POSBLOCFACTO
        INTEGER, intent(in)     :: LD_BLOCFACTO, IPIV(NPIV), 
     &                             MAXI_CLUSTER, IPANEL
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER LR_ACTIVATED_INT
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZET,
     &          IDEST, IPOSMSG, SSS, SSLR
        IERR = 0
        CALL MPI_PACK_SIZE( 6 + ( NDEST - 1 ) * OVHSIZE,
     &                      MPI_INTEGER, COMM, SIZE1, IERR_MPI )
        SIZE2  = 0
        CALL MPI_PACK_SIZE(2, MPI_INTEGER, COMM, SSLR, IERR_MPI )
        SIZE2=SIZE2+SSLR
        IF (.NOT. LR_ACTIVATED) THEN
        CALL MPI_PACK_SIZE( abs(NPIV)*NCOLU, MPI_DOUBLE_PRECISION,
     &                      COMM, SSLR, IERR_MPI )
         SIZE2=SIZE2+SSLR
        ELSE
          CALL MUMPS_MPI_PACK_SIZE_LR( BLR_LS, SSLR, COMM, IERR )
          SIZE2=SIZE2+SSLR
        ENDIF
        SIZET = SIZE1 + SIZE2
        IF (SIZET.GT.SIZE_RBUF_BYTES) THEN
         CALL MPI_PACK_SIZE( 6 ,
     &                      MPI_INTEGER, COMM, SSS, IERR_MPI )
         SSS = SSS+SIZE2
         IF (SSS.GT.SIZE_RBUF_BYTES) THEN
           IERR = -2
      RETURN
         ENDIF
        END IF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZET, IERR, 
     &                 NDEST, PDEST)
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( IPOSK, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JPOSK, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NPIV, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NCOLU, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        IF (LR_ACTIVATED) THEN
          LR_ACTIVATED_INT = 1
        ELSE
          LR_ACTIVATED_INT = 0
        ENDIF
        CALL MPI_PACK( LR_ACTIVATED_INT, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( IPANEL, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        IF (LR_ACTIVATED) THEN
                 CALL MUMPS_MPI_PACK_SCALE_LR( BLR_LS,
     &           BUF_CB%CONTENT( IPOSMSG:
     &                   IPOSMSG+(SIZET+KEEP(34)-1)/KEEP(34)-1 ),
     &           SIZET, POSITION, COMM,
     &           A, LA, POSBLOCFACTO, LD_BLOCFACTO, 
     &           IPIV, NPIV, MAXI_CLUSTER, IERR )
        ELSE
        CALL MPI_PACK( UIP21K, abs(NPIV) * NCOLU,
     &                        MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOSMSG ), SIZET,
     &                        POSITION, COMM, IERR_MPI )
        ENDIF
        DO IDEST = 1, NDEST
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, MPI_PACKED,
     &                  PDEST(IDEST), BLOC_FACTO_SYM_SLAVE, COMM,
     &                  BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     &                  IERR_MPI )
        END DO
        SIZET = SIZET - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZET .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blfac slave : size < position'
          WRITE(*,*) ' Size,position=',SIZET,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZET .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_BLFAC_SLAVE
        SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE3( N, ISON,
     &             NBCOL_SON, NBROW_SON, INDCOL_SON, INDROW_SON,
     &             LD_SON, VAL_SON, TAG, SUBSET_ROW, SUBSET_COL,
     &             NSUBSET_ROW, NSUBSET_COL,
     &             NSUPROW, NSUPCOL,
     &             NPROW, NPCOL, MBLOCK, RG2L_ROW, RG2L_COL,
     &             NBLOCK, PDEST, COMM, IERR , 
     &             TAB, TABSIZE, TRANSP, SIZE_PACK,
     &             N_ALREADY_SENT, KEEP, BBPCBP ) 
        IMPLICIT NONE
        INTEGER N, ISON, NBCOL_SON, NBROW_SON, NSUBSET_ROW, NSUBSET_COL
        INTEGER NPROW, NPCOL, MBLOCK, NBLOCK, LD_SON
        INTEGER BBPCBP
        INTEGER PDEST, TAG, COMM, IERR
        INTEGER INDCOL_SON( NBCOL_SON ), INDROW_SON( NBROW_SON )
        INTEGER SUBSET_ROW( NSUBSET_ROW ), SUBSET_COL( NSUBSET_COL )
        INTEGER :: RG2L_ROW(N)
        INTEGER :: RG2L_COL(N)
        INTEGER NSUPROW, NSUPCOL
        INTEGER(8), INTENT(IN) :: TABSIZE
        INTEGER SIZE_PACK
        INTEGER KEEP(500)
        DOUBLE PRECISION VAL_SON( LD_SON, * ), TAB(*)
        LOGICAL TRANSP
        INTEGER N_ALREADY_SENT
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE1, SIZE2, SIZE_AV, POSITION
        INTEGER SIZE_CBP, SIZE_TMP
        INTEGER IREQ, IPOS, ITAB
        INTEGER ISUB, JSUB, I, J 
        INTEGER ILOC_ROOT, JLOC_ROOT
        INTEGER IPOS_ROOT, JPOS_ROOT
        INTEGER IONE
        LOGICAL RECV_BUF_SMALLER_THAN_SEND
        INTEGER PDEST2(1)
        PARAMETER ( IONE=1 )
        INTEGER N_PACKET
        INTEGER NSUBSET_ROW_EFF, NSUBSET_COL_EFF, NSUPCOL_EFF
        PDEST2(1) = PDEST
        IERR = 0
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          CALL DMUMPS_BUF_SIZE_AVAILABLE( BUF_CB, SIZE_AV )
          IF (SIZE_AV .LT. SIZE_RBUF_BYTES) THEN
            RECV_BUF_SMALLER_THAN_SEND = .FALSE.
          ELSE
            RECV_BUF_SMALLER_THAN_SEND = .TRUE.
            SIZE_AV = SIZE_RBUF_BYTES
          ENDIF
          SIZE_AV = min(SIZE_AV, SIZE_RBUF_BYTES)
          CALL MPI_PACK_SIZE(8 + NSUBSET_COL,
     &                      MPI_INTEGER, COMM, SIZE1, IERR_MPI )
          SIZE_CBP = 0
          IF (N_ALREADY_SENT .EQ. 0 .AND.
     &        min(NSUPROW,NSUPCOL) .GT.0) THEN
            CALL MPI_PACK_SIZE(NSUPROW, MPI_INTEGER, COMM,
     &           SIZE_CBP, IERR_MPI )
            CALL MPI_PACK_SIZE(NSUPCOL, MPI_INTEGER, COMM,
     &           SIZE_TMP, IERR_MPI )
            SIZE_CBP = SIZE_CBP + SIZE_TMP
            CALL MPI_PACK_SIZE(NSUPROW*NSUPCOL,
     &           MPI_DOUBLE_PRECISION, COMM,
     &           SIZE_TMP, IERR_MPI )
            SIZE_CBP = SIZE_CBP + SIZE_TMP
            SIZE1 = SIZE1 + SIZE_CBP
          ENDIF
          IF (BBPCBP.EQ.1) THEN
            NSUBSET_COL_EFF = NSUBSET_COL - NSUPCOL
            NSUPCOL_EFF = 0
          ELSE
            NSUBSET_COL_EFF = NSUBSET_COL
            NSUPCOL_EFF = NSUPCOL
          ENDIF
          NSUBSET_ROW_EFF = NSUBSET_ROW - NSUPROW
          N_PACKET =
     &    (SIZE_AV - SIZE1) / (SIZEofINT + NSUBSET_COL_EFF * SIZEofREAL)
 10       CONTINUE
          N_PACKET = min( N_PACKET,
     &                    NSUBSET_ROW_EFF-N_ALREADY_SENT )
          IF (N_PACKET .LE. 0 .AND.
     &        NSUBSET_ROW_EFF-N_ALREADY_SENT.GT.0) THEN
             IF (RECV_BUF_SMALLER_THAN_SEND) THEN
              IERR=-3
              GOTO 100
           ELSE
              IERR = -1
              GOTO 100
            ENDIF
          ENDIF
          CALL MPI_PACK_SIZE( 8 + NSUBSET_COL_EFF + N_PACKET,
     &                      MPI_INTEGER, COMM, SIZE1, IERR_MPI )
          SIZE1 = SIZE1 + SIZE_CBP
          CALL MPI_PACK_SIZE( N_PACKET * NSUBSET_COL_EFF,
     &                      MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR_MPI )
          SIZE_PACK = SIZE1 + SIZE2
          IF (SIZE_PACK .GT. SIZE_AV) THEN
            N_PACKET = N_PACKET - 1
            IF ( N_PACKET > 0 ) THEN
              GOTO 10
            ELSE
               IF (RECV_BUF_SMALLER_THAN_SEND) THEN
                IERR = -3
                GOTO 100
             ELSE
                IERR = -1
                GOTO 100
              ENDIF
            ENDIF
          ENDIF
          IF (N_PACKET + N_ALREADY_SENT .NE. NSUBSET_ROW - NSUPROW
     &         .AND.
     &         SIZE_PACK .LT. SIZE_RBUF_BYTES / 4
     &         .AND. .NOT. RECV_BUF_SMALLER_THAN_SEND)
     &         THEN
             IERR = -1
             GOTO 100
          ENDIF
        ELSE 
          N_PACKET = 0
          CALL MPI_PACK_SIZE(8,MPI_INTEGER, COMM, SIZE_PACK, IERR_MPI )
        END IF
        IF ( SIZE_PACK.GT.SIZE_RBUF_BYTES ) THEN
           IERR = -3
           GOTO 100
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE_PACK, IERR, 
     &                 IONE, PDEST2
     &               )
        IF ( IERR .LT. 0 ) GOTO 100
        POSITION = 0
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSUBSET_ROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSUPROW, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSUBSET_COL, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSUPCOL, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( N_ALREADY_SENT, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( N_PACKET, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( BBPCBP, 1, MPI_INTEGER,
     &                 BUF_CB%CONTENT( IPOS ),
     &                 SIZE_PACK, POSITION, COMM, IERR_MPI )
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          IF (N_ALREADY_SENT .EQ. 0 .AND.
     &          min(NSUPROW, NSUPCOL) .GT. 0) THEN
            DO ISUB = NSUBSET_ROW-NSUPROW+1, NSUBSET_ROW
              I =  SUBSET_ROW( ISUB )
              IPOS_ROOT = RG2L_ROW(INDCOL_SON( I ))
              ILOC_ROOT = MBLOCK
     &                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                 + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR_MPI )
            ENDDO
            DO ISUB = NSUBSET_COL-NSUPCOL+1, NSUBSET_COL
               J = SUBSET_COL( ISUB )
               JPOS_ROOT = INDROW_SON( J ) - N
               JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR_MPI )
            ENDDO
            IF ( TABSIZE.GE.int(NSUPROW,8)*int(NSUPCOL,8) ) THEN
              ITAB = 1
              DO JSUB = NSUBSET_ROW - NSUPROW+1, NSUBSET_ROW
                J = SUBSET_ROW(JSUB)
                DO ISUB = NSUBSET_COL - NSUPCOL+1, NSUBSET_COL
                  I = SUBSET_COL(ISUB)
                  TAB(ITAB) = VAL_SON(J, I)
                  ITAB = ITAB + 1
                ENDDO
              ENDDO
              CALL MPI_PACK(TAB(1), NSUPROW*NSUPCOL,
     &         MPI_DOUBLE_PRECISION, 
     &         BUF_CB%CONTENT( IPOS ),
     &         SIZE_PACK, POSITION, COMM, IERR_MPI )
            ELSE
              DO JSUB = NSUBSET_ROW - NSUPROW+1, NSUBSET_ROW
                J = SUBSET_ROW(JSUB)
                DO ISUB = NSUBSET_COL - NSUPCOL+1, NSUBSET_COL
                  I = SUBSET_COL(ISUB)
                  CALL MPI_PACK(VAL_SON(J,I), 1,
     &            MPI_DOUBLE_PRECISION, 
     &            BUF_CB%CONTENT( IPOS ),
     &            SIZE_PACK, POSITION, COMM, IERR_MPI )
                ENDDO
              ENDDO
            ENDIF
          ENDIF
          IF ( .NOT. TRANSP ) THEN
            DO ISUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
              I         = SUBSET_ROW( ISUB )
              IPOS_ROOT = RG2L_ROW( INDROW_SON( I ) )
              ILOC_ROOT = MBLOCK
     &                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                 + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR_MPI )
            END DO
            DO JSUB = 1, NSUBSET_COL_EFF - NSUPCOL_EFF
              J         = SUBSET_COL( JSUB )
              JPOS_ROOT = RG2L_COL( INDCOL_SON( J ) )
              JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR_MPI )
            END DO
            DO JSUB = NSUBSET_COL_EFF-NSUPCOL_EFF+1, NSUBSET_COL_EFF
               J = SUBSET_COL( JSUB )
               JPOS_ROOT = INDCOL_SON( J ) - N
               JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR_MPI )
            ENDDO
          ELSE
            DO JSUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
              J         = SUBSET_ROW( JSUB )
              IPOS_ROOT = RG2L_ROW( INDCOL_SON( J ) )
              ILOC_ROOT = MBLOCK
     &                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     &                 + mod( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     &                       BUF_CB%CONTENT( IPOS ),
     &                       SIZE_PACK, POSITION, COMM, IERR_MPI )
            END DO
            DO ISUB = 1, NSUBSET_COL_EFF - NSUPCOL_EFF
              I         = SUBSET_COL( ISUB )
              JPOS_ROOT = RG2L_COL( INDROW_SON( I ) )
              JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR_MPI )
            END DO
            DO ISUB = NSUBSET_COL_EFF - NSUPCOL_EFF + 1, NSUBSET_COL_EFF
              I         = SUBSET_COL( ISUB )
              JPOS_ROOT = INDROW_SON(I) - N
              JLOC_ROOT = NBLOCK
     &                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     &                  + mod( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     &                      BUF_CB%CONTENT( IPOS ),
     &                      SIZE_PACK, POSITION, COMM, IERR_MPI )
            ENDDO
          END IF
          IF ( TABSIZE.GE.int(N_PACKET,8)*int(NSUBSET_COL_EFF,8) ) THEN
            IF ( .NOT. TRANSP ) THEN
              ITAB = 1
              DO ISUB = N_ALREADY_SENT+1,
     &                  N_ALREADY_SENT+N_PACKET
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL_EFF
                  J              = SUBSET_COL( JSUB )
                  TAB( ITAB )    = VAL_SON(J,I)
                  ITAB           = ITAB + 1
                END DO
              END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL_EFF*N_PACKET,
     &         MPI_DOUBLE_PRECISION, 
     &         BUF_CB%CONTENT( IPOS ),
     &         SIZE_PACK, POSITION, COMM, IERR_MPI )
            ELSE
              ITAB = 1
              DO JSUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
                J = SUBSET_ROW( JSUB )
                DO ISUB = 1, NSUBSET_COL_EFF
                  I         = SUBSET_COL( ISUB )
                  TAB( ITAB ) = VAL_SON( J, I )
                  ITAB = ITAB + 1
                END DO
              END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL_EFF*N_PACKET,
     &         MPI_DOUBLE_PRECISION, 
     &         BUF_CB%CONTENT( IPOS ),
     &         SIZE_PACK, POSITION, COMM, IERR_MPI )
            END IF
          ELSE
            IF ( .NOT. TRANSP ) THEN
              DO ISUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL_EFF
                  J         = SUBSET_COL( JSUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     &            MPI_DOUBLE_PRECISION,
     &            BUF_CB%CONTENT( IPOS ),
     &            SIZE_PACK, POSITION, COMM, IERR_MPI )
                END DO
              END DO
            ELSE
              DO JSUB = N_ALREADY_SENT+1, N_ALREADY_SENT+N_PACKET
                J = SUBSET_ROW( JSUB )
                DO ISUB = 1, NSUBSET_COL_EFF
                  I         = SUBSET_COL( ISUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     &            MPI_DOUBLE_PRECISION,
     &            BUF_CB%CONTENT( IPOS ),
     &            SIZE_PACK, POSITION, COMM, IERR_MPI )
                END DO
              END DO
            END IF
          ENDIF
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                PDEST, TAG, COMM, BUF_CB%CONTENT( IREQ ),
     &                IERR_MPI )
        IF ( SIZE_PACK .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending contribution to root:Size<positn'
          WRITE(*,*) ' Size,position=',SIZE_PACK,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE_PACK .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_CB, POSITION )
        N_ALREADY_SENT = N_ALREADY_SENT + N_PACKET
        IF (NSUBSET_ROW * NSUBSET_COL .NE. 0) THEN
          IF ( N_ALREADY_SENT.NE.NSUBSET_ROW_EFF ) IERR = -1
        ENDIF
  100   CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE3
        SUBROUTINE DMUMPS_BUF_SEND_RTNELIND( ISON, NELIM,
     &             NELIM_ROW, NELIM_COL, NSLAVES, SLAVES,
     &             DEST, COMM, KEEP, IERR )
        INTEGER ISON, NELIM
        INTEGER NSLAVES, DEST, COMM, IERR
        INTEGER NELIM_ROW( NELIM ), NELIM_COL( NELIM )
        INTEGER SLAVES( NSLAVES )
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE = ( 3 + NSLAVES + 2 * NELIM ) * SIZEofINT
        IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
             IERR = -3
      RETURN
        ENDIF
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = ISON
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NELIM
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_ROW
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_COL
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 ) = SLAVES
        POSITION = POSITION + NSLAVES
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
          WRITE(*,*) 'Error in DMUMPS_BUF_SEND_ROOT_NELIM_INDICES:',
     &               'wrong estimated size'
          CALL MUMPS_ABORT()
        END IF
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE, 
     &                  MPI_PACKED,
     &                  DEST, ROOT_NELIM_INDICES, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_RTNELIND
        SUBROUTINE DMUMPS_BUF_SEND_ROOT2SON( ISON, NELIM_ROOT,
     &             DEST, COMM, KEEP, IERR )
        IMPLICIT NONE
        INTEGER ISON, NELIM_ROOT, DEST, COMM, IERR
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER IPOS, IREQ, SIZE
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        SIZE = 2 * SIZEofINT
        CALL BUF_LOOK( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          WRITE(*,*) 'Internal error 1 with small buffers '
          CALL MUMPS_ABORT()
        END IF
        IF ( IERR .LT. 0 ) THEN
          RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS )     = ISON
        BUF_SMALL%CONTENT( IPOS + 1 ) = NELIM_ROOT
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     &                  MPI_PACKED,
     &                  DEST, ROOT_2SON, COMM,
     &                  BUF_SMALL%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_ROOT2SON
        SUBROUTINE DMUMPS_BUF_SEND_ROOT2SLAVE
     &  ( TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, KEEP, IERR )
        IMPLICIT NONE
        INTEGER TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, IERR
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        SIZE = 2 * SIZEofINT
        CALL BUF_LOOK( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          WRITE(*,*) 'Internal error 2 with small buffers '
          CALL MUMPS_ABORT()
        END IF
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS     ) = TOT_ROOT_SIZE
        BUF_SMALL%CONTENT( IPOS + 1 ) = TOT_CONT2RECV
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     &                  MPI_PACKED,
     &                  DEST, ROOT_2SLAVE, COMM,
     &                  BUF_SMALL%CONTENT( IREQ ), IERR_MPI )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_ROOT2SLAVE
        SUBROUTINE DMUMPS_BUF_SEND_BACKVEC
     &             ( NRHS, INODE, W, LW, LD_W, DEST, MSGTAG,
     &               JBDEB, JBFIN, KEEP, COMM, IERR )
        IMPLICIT NONE
        INTEGER NRHS, INODE,LW,COMM,IERR,DEST,MSGTAG, LD_W
        INTEGER, intent(in) :: JBDEB, JBFIN
        DOUBLE PRECISION :: W(LD_W, *)
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INTEGER :: IERR_MPI
        INTEGER SIZE, SIZE1, SIZE2
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE, K
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        CALL MPI_PACK_SIZE( 4 , MPI_INTEGER, COMM, SIZE1, IERR_MPI )
        CALL MPI_PACK_SIZE( LW*NRHS, MPI_DOUBLE_PRECISION, COMM,
     &                      SIZE2, IERR_MPI )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LW   , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBDEB   , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( JBFIN   , 1, MPI_INTEGER,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        DO K=1, NRHS
        CALL MPI_PACK( W(1,K), LW, MPI_DOUBLE_PRECISION,
     &                        BUF_CB%CONTENT( IPOS ), SIZE,
     &                        POSITION, COMM, IERR_MPI )
        END DO
        KEEP(266)=KEEP(266)+1
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     &                  DEST, MSGTAG, COMM,
     &                  BUF_CB%CONTENT( IREQ ), IERR_MPI )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_update: SIZE, POSITION = ',
     &               SIZE, POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL BUF_ADJUST( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_BACKVEC
        SUBROUTINE DMUMPS_BUF_SEND_UPDATE_LOAD
     &             ( BDC_SBTR,BDC_MEM,BDC_MD, COMM, NPROCS, LOAD,
     &               MEM,SBTR_CUR,
     &               LU_USAGE,
     &               FUTURE_NIV2,
     &               MYID, KEEP, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER FUTURE_NIV2(NPROCS)
        DOUBLE PRECISION LU_USAGE
        DOUBLE PRECISION LOAD
        DOUBLE PRECISION MEM,SBTR_CUR
        LOGICAL BDC_MEM,BDC_SBTR,BDC_MD
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, WHAT, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        NDEST = 0
        DO I = 1, NPROCS
           IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
              NDEST = NDEST + 1
           ENDIF
        ENDDO
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR_MPI )
        NREALS = 1
        IF (BDC_MEM) THEN
          NREALS = 2
        ENDIf
        IF (BDC_SBTR)THEN
          NREALS = 3
        ENDIF
        IF(BDC_MD)THEN
           NREALS=NREALS+1
        ENDIF
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR_MPI )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        WHAT = 0  
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM, 1, MPI_DOUBLE_PRECISION,
     &                   BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                   POSITION, COMM, IERR_MPI )
        END IF
        IF (BDC_SBTR) THEN
          CALL MPI_PACK( SBTR_CUR, 1, MPI_DOUBLE_PRECISION,
     &                   BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                   POSITION, COMM, IERR_MPI )
        END IF
        IF(BDC_MD)THEN
           CALL MPI_PACK( LU_USAGE, 1, MPI_DOUBLE_PRECISION,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
            IDEST = IDEST + 1
            KEEP(267)=KEEP(267)+1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                      POSITION, MPI_PACKED, I,
     &                      UPDATE_LOAD, COMM,
     &                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                      IERR_MPI )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_SEND_UPDATE_LOAD'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_UPDATE_LOAD
        SUBROUTINE DMUMPS_BUF_BROADCAST
     &             ( WHAT, COMM, NPROCS, 
     &               FUTURE_NIV2,
     &               LOAD, UPD_LOAD,
     &               MYID, KEEP, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR, WHAT
        DOUBLE PRECISION LOAD,UPD_LOAD
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
        INTEGER FUTURE_NIV2(NPROCS)
        PARAMETER ( IZERO=0 )
        IERR = 0
        IF (WHAT .NE. 2 .AND. WHAT .NE. 3 .AND.
     &       WHAT.NE.6.AND. WHAT.NE.8 .AND.WHAT.NE.9.AND.
     &       WHAT.NE.17) THEN
          WRITE(*,*)
     &          "Internal error 1 in DMUMPS_BUF_BROADCAST",WHAT
        END IF
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR_MPI )
        IF((WHAT.NE.17).AND.(WHAT.NE.10))THEN
           NREALS = 1
        ELSE
           NREALS = 2
        ENDIF
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR_MPI )
        SIZE = SIZE1 + SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        IF((WHAT.EQ.17).OR.(WHAT.EQ.10))THEN
           CALL MPI_PACK( UPD_LOAD, 1, MPI_DOUBLE_PRECISION,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
          IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
            IDEST = IDEST + 1
            KEEP(267)=KEEP(267)+1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                      POSITION, MPI_PACKED, I,
     &                      UPDATE_LOAD, COMM,
     &                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                      IERR_MPI )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_BROADCAST'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_BROADCAST
        SUBROUTINE DMUMPS_BUF_SEND_FILS
     &             ( WHAT, COMM, NPROCS,
     &               FATHER_NODE,INODE,NCB,KEEP,
     &               MYID,REMOTE, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR, WHAT,REMOTE
        INTEGER FATHER_NODE,INODE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER :: IERR_MPI
        INTEGER POSITION, IREQ, IPOS, SIZE
        INTEGER NDEST, IDEST, IPOSMSG
        INTEGER IZERO,NCB,KEEP(500)
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        MYID2(1) = MYID
        NDEST = 1
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
           CALL MPI_PACK_SIZE( 4 + OVHSIZE, 
     &          MPI_INTEGER, COMM,
     &          SIZE, IERR_MPI )
        ELSE
           CALL MPI_PACK_SIZE( 2 + OVHSIZE, 
     &          MPI_INTEGER, COMM,
     &          SIZE, IERR_MPI )
        ENDIF
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( FATHER_NODE, 1, MPI_INTEGER,
     &                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &                 POSITION, COMM, IERR_MPI )
        IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
           CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
           CALL MPI_PACK( NCB, 1, MPI_INTEGER,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
        ENDIF
        IDEST = 1
        KEEP(267)=KEEP(267)+1
        CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                 POSITION, MPI_PACKED, REMOTE,
     &                 UPDATE_LOAD, COMM,
     &                 BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                 IERR_MPI )
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_SEND_FILS'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_FILS
        SUBROUTINE DMUMPS_BUF_SEND_NOT_MSTR( COMM, MYID, NPROCS,
     &  MAX_SURF_MASTER, KEEP, IERR)
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER COMM, MYID, IERR, NPROCS
        DOUBLE PRECISION MAX_SURF_MASTER
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER :: IERR_MPI
        INTEGER IPOS, IREQ, IDEST, IPOSMSG, POSITION, I
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        INTEGER NDEST, NINTS, NREALS, SIZE, SIZE1, SIZE2
        INTEGER WHAT
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        NINTS = 1 + ( NDEST-1 ) * OVHSIZE
        NREALS = 1
        CALL MPI_PACK_SIZE( NINTS, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR_MPI )
        CALL MPI_PACK_SIZE( NREALS,
     &                       MPI_DOUBLE_PRECISION, COMM,
     &                       SIZE2, IERR_MPI )
        SIZE=SIZE1+SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     &       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        WHAT = 4
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( MAX_SURF_MASTER, 1, MPI_DOUBLE_PRECISION,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        IDEST = 0
        DO I = 0, NPROCS - 1
           IF ( I .ne. MYID ) THEN
              IDEST = IDEST + 1
              KEEP(267)=KEEP(267)+1
              CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &             POSITION, MPI_PACKED, I,
     &             UPDATE_LOAD, COMM,
     &             BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &             IERR_MPI )
           END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_BCAST_ARRAY'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_SEND_NOT_MSTR
        SUBROUTINE DMUMPS_BUF_BCAST_ARRAY( BDC_MEM,
     &      COMM, MYID, NPROCS,
     &      FUTURE_NIV2,
     &      NSLAVES,
     &      LIST_SLAVES,INODE,
     &      MEM_INCREMENT, FLOPS_INCREMENT,CB_BAND, WHAT,
     &      KEEP,
     &      IERR )
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        LOGICAL BDC_MEM
        INTEGER COMM, MYID, NPROCS, NSLAVES, IERR
        INTEGER FUTURE_NIV2(NPROCS)
        INTEGER LIST_SLAVES(NSLAVES),INODE
        DOUBLE PRECISION MEM_INCREMENT(NSLAVES)
        DOUBLE PRECISION FLOPS_INCREMENT(NSLAVES)
        DOUBLE PRECISION CB_BAND(NSLAVES)
        INTEGER, INTENT(INOUT) :: KEEP(500)
        INTEGER :: IERR_MPI
        INTEGER NDEST, NINTS, NREALS, SIZE1, SIZE2, SIZE
        INTEGER IPOS, IPOSMSG, IREQ, POSITION
        INTEGER I, IDEST, WHAT
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        MYID2(1)=MYID
        IERR = 0
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
        IF ( NDEST == 0 ) THEN
           RETURN
        ENDIF
        NINTS = 2 +  NSLAVES + ( NDEST - 1 ) * OVHSIZE + 1
        NREALS = NSLAVES
        IF (BDC_MEM) NREALS = NREALS + NSLAVES
        IF(WHAT.EQ.19) THEN 
           NREALS = NREALS + NSLAVES
        ENDIF
        CALL MPI_PACK_SIZE( NINTS, 
     &                       MPI_INTEGER, COMM,
     &                       SIZE1, IERR_MPI )
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &       COMM, SIZE2, IERR_MPI )
        SIZE = SIZE1+SIZE2
        CALL BUF_LOOK( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     &       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) THEN
      RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     &    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( LIST_SLAVES, NSLAVES, MPI_INTEGER,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( FLOPS_INCREMENT, NSLAVES,
     &      MPI_DOUBLE_PRECISION,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM_INCREMENT, NSLAVES,
     &      MPI_DOUBLE_PRECISION,
     &      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &      POSITION, COMM, IERR_MPI )
        END IF
        IF(WHAT.EQ.19)THEN
           CALL MPI_PACK( CB_BAND, NSLAVES,
     &          MPI_DOUBLE_PRECISION,
     &          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     &          POSITION, COMM, IERR_MPI )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
            IDEST = IDEST + 1
            KEEP(267)=KEEP(267)+1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     &                      POSITION, MPI_PACKED, I,
     &                      UPDATE_LOAD, COMM,
     &                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     &                      IERR_MPI )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_BCAST_ARRAY'
          WRITE(*,*) ' Size,position=',SIZE,POSITION
          CALL MUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     &  CALL BUF_ADJUST( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_BUF_BCAST_ARRAY
        SUBROUTINE DMUMPS_BUF_DIST_IRECV_SIZE
     &             ( DMUMPS_LBUFR_BYTES)
        IMPLICIT NONE
        INTEGER DMUMPS_LBUFR_BYTES 
        SIZE_RBUF_BYTES = DMUMPS_LBUFR_BYTES
        RETURN
      END SUBROUTINE DMUMPS_BUF_DIST_IRECV_SIZE
      SUBROUTINE MUMPS_MPI_PACK_SIZE_LR( BLR_LorU, SIZE_OUT, COMM,
     &                                   IERR )
      USE DMUMPS_LR_TYPE
      INTEGER, intent(out) :: SIZE_OUT, IERR
      INTEGER, intent(in)  :: COMM
      TYPE (LRB_TYPE), DIMENSION(:), intent(in) :: BLR_LorU
      INTEGER :: I, SIZE_LOC, IERR_MPI
      INCLUDE 'mpif.h'
        IERR = 0
        SIZE_OUT = 0
        CALL MPI_PACK_SIZE( 1, MPI_INTEGER, COMM, SIZE_LOC,  IERR_MPI )
        SIZE_OUT = SIZE_OUT + SIZE_LOC
        DO I = 1, size(BLR_LorU)
          CALL MUMPS_MPI_PACK_SIZE_LRB(BLR_LorU(I), SIZE_LOC, COMM, 
     &                                 IERR )
          SIZE_OUT = SIZE_OUT + SIZE_LOC
        ENDDO
        RETURN
      END SUBROUTINE MUMPS_MPI_PACK_SIZE_LR
      SUBROUTINE MUMPS_MPI_PACK_SIZE_LRB(LRB, SIZE_OUT, COMM, IERR )
      USE DMUMPS_LR_TYPE
      INTEGER, intent(out) :: SIZE_OUT, IERR
      INTEGER, intent(in)  :: COMM
      TYPE (LRB_TYPE), intent(in) :: LRB
      INTEGER :: SIZE_LOC, IERR_MPI
      INCLUDE 'mpif.h'
      IERR = 0
      SIZE_OUT = 0
      CALL MPI_PACK_SIZE( 4,      
     &          MPI_INTEGER, COMM, SIZE_LOC,  IERR_MPI )
      SIZE_OUT = SIZE_OUT + SIZE_LOC
      IF ( LRB%ISLR ) THEN
        IF (LRB%K .GT. 0) THEN
          CALL MPI_PACK_SIZE( LRB%M * LRB%K,
     &       MPI_DOUBLE_PRECISION, COMM, SIZE_LOC,  IERR_MPI )
          SIZE_OUT = SIZE_OUT + SIZE_LOC
          CALL MPI_PACK_SIZE( LRB%K * LRB%N,
     &       MPI_DOUBLE_PRECISION, COMM, SIZE_LOC,  IERR_MPI )
          SIZE_OUT = SIZE_OUT + SIZE_LOC
        ENDIF
      ELSE
        CALL MPI_PACK_SIZE( LRB%M * LRB%N,
     &       MPI_DOUBLE_PRECISION, COMM, SIZE_LOC,  IERR_MPI )
        SIZE_OUT = SIZE_OUT + SIZE_LOC
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_MPI_PACK_SIZE_LRB
      SUBROUTINE DMUMPS_MPI_PACK_LR( BLR_LorU, BUF, LBUF, POSITION,
     &                              COMM, IERR )
      USE DMUMPS_LR_TYPE
      INTEGER, intent(out) :: IERR
      INTEGER, intent(in)  :: COMM, LBUF  
      INTEGER, intent(inout) :: POSITION
      INTEGER, intent(inout) :: BUF(:) 
      TYPE (LRB_TYPE), DIMENSION(:), intent(in) :: BLR_LorU
      INTEGER I
      INTEGER :: IERR_MPI
      INCLUDE 'mpif.h'
      IERR = 0
      CALL MPI_PACK( size(BLR_LorU), 1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      DO I = 1, size(BLR_LorU)
        CALL DMUMPS_MPI_PACK_LRB(BLR_LorU(I), BUF, LBUF, POSITION, 
     &                          COMM, IERR )
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_MPI_PACK_LR
      SUBROUTINE DMUMPS_MPI_PACK_LRB( LRB, BUF, LBUF, POSITION,
     &                              COMM, IERR )
      USE DMUMPS_LR_TYPE
      INTEGER, intent(out) :: IERR
      INTEGER, intent(in)  :: COMM, LBUF  
      INTEGER, intent(inout) :: POSITION
      INTEGER, intent(inout) :: BUF(:) 
      TYPE (LRB_TYPE), intent(in) :: LRB
      INTEGER ISLR_INT
      INTEGER :: IERR_MPI
      INCLUDE 'mpif.h'
      IERR = 0
      IF (LRB%ISLR) THEN
        ISLR_INT = 1
      ELSE
        ISLR_INT = 0
      ENDIF
      CALL MPI_PACK( ISLR_INT, 1, MPI_INTEGER,
     &     BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      CALL MPI_PACK( LRB%K,
     &     1, MPI_INTEGER,
     &     BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      CALL MPI_PACK( LRB%M, 
     &     1, MPI_INTEGER,
     &     BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      CALL MPI_PACK( LRB%N,
     &     1, MPI_INTEGER,
     &     BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      IF (LRB%ISLR) THEN
        IF (LRB%K .GT. 0) THEN
          CALL MPI_PACK( LRB%Q(1,1), 
     &      LRB%M*LRB%K, MPI_DOUBLE_PRECISION,
     &      BUF(1), LBUF, POSITION, COMM, IERR_MPI )
          CALL MPI_PACK( LRB%R(1,1),
     &      LRB%N*LRB%K, MPI_DOUBLE_PRECISION,
     &      BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        ENDIF
      ELSE
        CALL MPI_PACK( LRB%Q(1,1), LRB%M*LRB%N
     &     ,MPI_DOUBLE_PRECISION,
     &     BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_MPI_PACK_LRB
      SUBROUTINE DMUMPS_MPI_UNPACK_LRB(
     &           BUFR, LBUFR, LBUFR_BYTES, POSITION,
     &                             LRB, KEEP8,
     &                             COMM, IFLAG, IERROR)
      USE DMUMPS_LR_CORE, ONLY : LRB_TYPE, ALLOC_LRB
      USE DMUMPS_LR_TYPE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LBUFR
      INTEGER, INTENT(IN) :: LBUFR_BYTES
      INTEGER, INTENT(IN) :: BUFR(LBUFR)
      INTEGER, INTENT(INOUT) :: POSITION
      INTEGER, INTENT(IN) :: COMM
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      TYPE (LRB_TYPE), INTENT(OUT) :: LRB
      INTEGER(8) :: KEEP8(150)
      LOGICAL :: ISLR
      INTEGER :: ISLR_INT
      INTEGER :: K, M, N
      INTEGER :: IERR_MPI
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               ISLR_INT, 1, MPI_INTEGER, COMM, IERR_MPI )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               K, 1,
     &               MPI_INTEGER, COMM, IERR_MPI )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               M, 1,
     &               MPI_INTEGER, COMM, IERR_MPI )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               N, 1,
     &               MPI_INTEGER, COMM, IERR_MPI )
      IF (ISLR_INT .eq. 1) THEN
        ISLR = .TRUE.
      ELSE
        ISLR = .FALSE.
      ENDIF
      CALL ALLOC_LRB( LRB, K, M, N, ISLR, 
     &           IFLAG, IERROR, KEEP8 )
      IF (IFLAG.LT.0) RETURN
      IF (ISLR) THEN
        IF (K .GT. 0) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   LRB%Q(1,1), M*K, MPI_DOUBLE_PRECISION,
     &                   COMM, IERR_MPI )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   LRB%R(1,1), N*K, MPI_DOUBLE_PRECISION,
     &                   COMM, IERR_MPI )
        ENDIF
      ELSE
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   LRB%Q(1,1), M*N, MPI_DOUBLE_PRECISION,
     &                   COMM, IERR_MPI )
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_MPI_UNPACK_LRB
      SUBROUTINE MUMPS_MPI_PACK_SCALE_LR
     &                  ( BLR, BUF, LBUF, POSITION,
     &                    COMM, 
     &                    A , LA, POSELTD, LD_DIAG,
     &                    IPIV, NPIV, MAXI_CLUSTER, 
     &                    IERR )
      USE DMUMPS_LR_TYPE
      INTEGER, intent(out) :: IERR
      INTEGER, intent(in)  :: COMM, LBUF 
      INTEGER, intent(inout) :: POSITION
      INTEGER, intent(inout) :: BUF(:)  
      TYPE  (LRB_TYPE), DIMENSION(:), intent(in) :: BLR
      INTEGER(8), intent(in)  :: LA, POSELTD
      INTEGER, intent(in)     :: LD_DIAG, NPIV 
      INTEGER, intent(in)     :: IPIV(NPIV), MAXI_CLUSTER
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      INTEGER :: IERR_MPI
      INTEGER I, ISLR_INT, J, ALLOCOK
      DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) ::  SCALED
      DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:) ::  BLOCK
      DOUBLE PRECISION :: PIV1, PIV2, OFFDIAG
      INCLUDE 'mpif.h'
      IERR = 0
      CALL MPI_PACK( size(BLR), 1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
      allocate(BLOCK(MAXI_CLUSTER), STAT=ALLOCOK )
      IF ( ALLOCOK .GT. 0 ) THEN
             WRITE(*,*) 'pb allocation in mumps_mpi_pack_scale_lr'
             IERR = -1
             GOTO 500
      END IF
      allocate(SCALED(MAXI_CLUSTER,2), STAT=ALLOCOK )
      IF ( ALLOCOK .GT. 0 ) THEN
             WRITE(*,*) 'pb allocation in mumps_mpi_pack_scale_lr'
             IERR = -1
             GOTO 500
      END IF
      DO I = 1, size(BLR)
        IF (BLR(I)%ISLR) THEN
          ISLR_INT = 1
        ELSE
          ISLR_INT = 0
        ENDIF
        CALL MPI_PACK( ISLR_INT, 1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( BLR(I)%K,
     &       1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( BLR(I)%M, 
     &       1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        CALL MPI_PACK( BLR(I)%N,
     &       1, MPI_INTEGER,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
        IF (BLR(I)%ISLR) THEN
          IF (BLR(I)%K .GT. 0) THEN
            CALL MPI_PACK( BLR(I)%Q(1,1), BLR(I)%M*BLR(I)%K,
     &       MPI_DOUBLE_PRECISION,
     &       BUF(1), LBUF, POSITION, COMM, IERR_MPI )
            J =1
          DO WHILE (J <= BLR(I)%N)
              IF (IPIV(J) > 0) THEN
                SCALED(1:BLR(I)%K,1) = A(POSELTD+LD_DIAG*(J-1)+J-1) 
     &            * BLR(I)%R(1:BLR(I)%K,J)
                J = J+1
              CALL MPI_PACK( SCALED(1,1), BLR(I)%K,
     &           MPI_DOUBLE_PRECISION,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
              ELSE 
                PIV1    = A(POSELTD+LD_DIAG*(J-1)+J-1)
                PIV2    = A(POSELTD+LD_DIAG*J+J)
                OFFDIAG = A(POSELTD+LD_DIAG*(J-1)+J)
                BLOCK(1:BLR(I)%K)    = BLR(I)%R(1:BLR(I)%K,J)
                SCALED(1:BLR(I)%K,1) = PIV1 * BLR(I)%R(1:BLR(I)%K,J)
     &            + OFFDIAG * BLR(I)%R(1:BLR(I)%K,J+1)
                CALL MPI_PACK( SCALED(1,1), BLR(I)%K,
     &           MPI_DOUBLE_PRECISION,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                SCALED(1:BLR(I)%K,2) = OFFDIAG * BLOCK(1:BLR(I)%K)
     &            + PIV2 * BLR(I)%R(1:BLR(I)%K,J+1)
                 J =J+2
                CALL MPI_PACK( SCALED(1,2), BLR(I)%K,
     &           MPI_DOUBLE_PRECISION,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
              ENDIF
          END DO
        ENDIF
        ELSE
          J = 1
          DO WHILE (J <= BLR(I)%N)
              IF (IPIV(J) > 0) THEN
                SCALED(1:BLR(I)%M,1) = A(POSELTD+LD_DIAG*(J-1)+J-1) 
     &           * BLR(I)%Q(1:BLR(I)%M,J)
                CALL MPI_PACK( SCALED(1,1), BLR(I)%M,
     &           MPI_DOUBLE_PRECISION,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                J = J+1
              ELSE 
                PIV1    = A(POSELTD+LD_DIAG*(J-1)+J-1)
                PIV2    = A(POSELTD+LD_DIAG*J+J)
                OFFDIAG = A(POSELTD+LD_DIAG*(J-1)+J)
                BLOCK(1:BLR(I)%M)    = BLR(I)%Q(1:BLR(I)%M,J)
                SCALED(1:BLR(I)%M,1) = PIV1 * BLR(I)%Q(1:BLR(I)%M,J)
     &            + OFFDIAG * BLR(I)%Q(1:BLR(I)%M,J+1)
                CALL MPI_PACK( SCALED(1,1), BLR(I)%M,
     &           MPI_DOUBLE_PRECISION,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                SCALED(1:BLR(I)%M,2) = OFFDIAG * BLOCK(1:BLR(I)%M)
     &            + PIV2 * BLR(I)%Q(1:BLR(I)%M,J+1)
                CALL MPI_PACK( SCALED(1,2), BLR(I)%M,
     &           MPI_DOUBLE_PRECISION,
     &           BUF(1), LBUF, POSITION, COMM, IERR_MPI )
                 J=J+2
              ENDIF
          END DO
        ENDIF
      ENDDO
 500  CONTINUE
      IF (allocated(BLOCK)) deallocate(BLOCK)
      IF (allocated(SCALED)) deallocate(SCALED)
      RETURN
      END SUBROUTINE MUMPS_MPI_PACK_SCALE_LR
      END MODULE DMUMPS_BUF
