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
      SUBROUTINE DMUMPS_BUILD_MAPPING
     & ( N, MAPPING, NNZ, IRN, JCN, PROCNODE, STEP,
     &   SLAVEF, PERM, FILS,
     &   RG2L, KEEP,KEEP8, MBLOCK, NBLOCK, NPROW, NPCOL )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N, SLAVEF, MBLOCK, NBLOCK, NPROW, NPCOL
      iNTEGER(8) :: NNZ
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER IRN( NNZ ), JCN( NNZ ) 
      INTEGER MAPPING( NNZ ), STEP( N )
      INTEGER PROCNODE( KEEP(28) ), PERM( N ), FILS( N ), RG2L( N )
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE
      INTEGER K4, IOLD, JOLD, INEW, JNEW, ISEND, JSEND, IARR, INODE
      INTEGER(8) :: K8
      INTEGER TYPE_NODE, DEST
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID
      INODE = KEEP(38)
      K4 = 1
      DO WHILE ( INODE .GT. 0 )
        RG2L( INODE ) = K4
        INODE = FILS( INODE )
        K4 = K4 + 1
      END DO
      DO K8 = 1_8, NNZ
        IOLD = IRN( K8 )
        JOLD = JCN( K8 )
        IF ( IOLD .GT. N .OR. IOLD .LT. 1 .OR.
     &       JOLD .GT. N .OR. JOLD .LT. 1 ) THEN
           MAPPING( K8 ) = -1
           CYCLE
        END IF
        IF ( IOLD .eq. JOLD ) THEN
          ISEND = IOLD
          JSEND = JOLD
        ELSE
          INEW = PERM( IOLD )
          JNEW = PERM( JOLD )
          IF ( INEW .LT. JNEW ) THEN
            ISEND = IOLD
            IF ( KEEP(50) .ne. 0 ) ISEND = -IOLD
            JSEND = JOLD
          ELSE
            ISEND = -JOLD
            JSEND = IOLD
          END IF
        END IF
        IARR = abs( ISEND )
        TYPE_NODE = MUMPS_TYPENODE( PROCNODE(abs(STEP(IARR))),
     &                              KEEP(199) )
        IF ( TYPE_NODE .eq. 1 .or. TYPE_NODE .eq. 2 ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            DEST = MUMPS_PROCNODE( PROCNODE(abs(STEP(IARR))),
     &                             KEEP(199) ) + 1
          ELSE
            DEST = MUMPS_PROCNODE( PROCNODE(abs(STEP(IARR))),
     &                             KEEP(199) )
          END IF
        ELSE
          IF ( ISEND .LT. 0 ) THEN
            IPOSROOT = RG2L( JSEND )
            JPOSROOT = RG2L( IARR  )
          ELSE
            IPOSROOT = RG2L( IARR  )
            JPOSROOT = RG2L( JSEND )
          END IF
          IROW_GRID = mod( ( IPOSROOT - 1 )/MBLOCK, NPROW )
          JCOL_GRID = mod( ( JPOSROOT - 1 )/NBLOCK, NPCOL )
          IF ( KEEP( 46 ) .eq. 0 ) THEN
            DEST = IROW_GRID * NPCOL + JCOL_GRID + 1
          ELSE
            DEST = IROW_GRID * NPCOL + JCOL_GRID
          END IF
        END IF
        MAPPING( K8 ) = DEST
      END DO
      RETURN
      END SUBROUTINE DMUMPS_BUILD_MAPPING
      SUBROUTINE DMUMPS_REDISTRIBUTION(
     & N, NZ_loc8, id,
     & DBLARR, LDBLARR, INTARR, LINTARR,
     & PTRAIW, PTRARW, KEEP,KEEP8, MYID, COMM, NBRECORDS,
     &
     & A, LA, root, PROCNODE_STEPS, SLAVEF, PERM, STEP,
     & ICNTL, INFO, NSEND8, NLOCAL8,
     & ISTEP_TO_INIV2, CANDIDATES
     & )
!$    USE OMP_LIB
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N
      INTEGER(8) :: NZ_loc8
      TYPE (DMUMPS_STRUC) :: id
      INTEGER(8) :: LDBLARR, LINTARR
      DOUBLE PRECISION DBLARR( LDBLARR )
      INTEGER INTARR( LINTARR )
      INTEGER(8), INTENT(IN) :: PTRAIW( N ), PTRARW( N )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER MYID, COMM, NBRECORDS
      INTEGER(8) :: LA
      INTEGER SLAVEF
      INTEGER ISTEP_TO_INIV2(KEEP(71))
      INTEGER CANDIDATES(SLAVEF+1, max(1,KEEP(56)))
      DOUBLE PRECISION A( LA )
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER PROCNODE_STEPS(KEEP(28)), PERM( N ), STEP( N )
      INTEGER INFO( 80 ), ICNTL(60)
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE, numroc, 
     &        MUMPS_TYPESPLIT
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, numroc, 
     &        MUMPS_TYPESPLIT
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER :: IERR, MSGSOU
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IW4
      INTEGER END_MSG_2_RECV
      INTEGER I
      INTEGER(8) :: I18, IA8
      INTEGER(8) :: K8
      INTEGER TYPE_NODE, DEST
      INTEGER IOLD, JOLD, IARR, ISEND, JSEND
      INTEGER allocok,  TYPESPLIT, T4MASTER, INIV2, NCAND
      LOGICAL T4_MASTER_CONCERNED, EARLYT3ROOTINS
      DOUBLE PRECISION VAL
      INTEGER(8) :: PTR_ROOT
      INTEGER LOCAL_M, LOCAL_N, ARROW_ROOT
      INTEGER IROW_GRID, JCOL_GRID, IPOSROOT, JPOSROOT,
     &        ILOCROOT, JLOCROOT
      INTEGER MP,LP
      INTEGER KPROBE, FREQPROBE
      INTEGER(8) :: IS18, IIW8, IS8, IAS8
      INTEGER ISHIFT
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: BUFI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: BUFR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFRECI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BUFRECR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IACT, IREQI, IREQR
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: SEND_ACTIVE
      LOGICAL :: FLAG
      INTEGER(8), INTENT(OUT) :: NSEND8, NLOCAL8
      INTEGER MASTER_NODE, ISTEP
      LOGICAL :: DOIT, OMP_FLAG, OMP_FLAG_P
      INTEGER NOMP, NOMP_P, IOMP, P2
      NSEND8  = 0_8
      NLOCAL8 = 0_8
      LP = ICNTL(1)
      MP = ICNTL(2)
      END_MSG_2_RECV = SLAVEF
      ALLOCATE( IACT(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IACT in matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( IREQI(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IREQI in matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( IREQR(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IREQR in matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( SEND_ACTIVE(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating SEND_ACTIVE in matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( BUFI( NBRECORDS * 2 + 1, 2, SLAVEF ), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating int buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = ( NBRECORDS * 2 + 1 ) * SLAVEF * 2
        GOTO 20
      END IF
      ALLOCATE( BUFR( NBRECORDS, 2, SLAVEF), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating real buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * SLAVEF * 2
        GOTO 20
      END IF
      ALLOCATE( BUFRECI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * 2 + 1
        GOTO 20
      END IF
      ALLOCATE( BUFRECR( NBRECORDS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS
        GOTO 20
      END IF
      ALLOCATE( IW4( N, 2 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(LP,*) '** Error allocating IW4 for matrix distribution'
        INFO(1) = -13
        INFO(2) = N * 2
      END IF
 20   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 100
      ARROW_ROOT = 0
      DO I = 1, N
          I18 = PTRAIW( I )
          IA8 = PTRARW( I )
          IF ( IA8 .GT. 0_8 ) THEN
            DBLARR( IA8 ) = ZERO
            IW4( I, 1 ) = INTARR( I18 )
            IW4( I, 2 ) = -INTARR( I18 + 1_8 )
            INTARR( I18 + 2_8 ) = I
          END IF
      END DO
      EARLYT3ROOTINS = KEEP(200) .EQ.0
      IF ( KEEP(38) .NE. 0 .AND. EARLYT3ROOTINS ) THEN
        CALL DMUMPS_GET_ROOT_INFO(root,LOCAL_M, LOCAL_N, PTR_ROOT, LA)
        CALL DMUMPS_SET_ROOT_TO_ZERO(root, KEEP, A, LA)
      ELSE
        LOCAL_M = -19999; LOCAL_N = -29999; PTR_ROOT = -99999_8
      END IF
      DO I = 1, SLAVEF
        BUFI( 1, 1, I ) = 0
      END DO
      DO I = 1, SLAVEF
        BUFI( 1, 2, I ) = 0
      END DO
      DO I = 1, SLAVEF
        SEND_ACTIVE( I ) = .FALSE.
        IACT( I ) = 1
      END DO
      KPROBE = 0
      FREQPROBE = max(1,NBRECORDS/10)
      IF (SLAVEF .EQ. 1) FREQPROBE = huge(FREQPROBE)
      NOMP = 1
!$    NOMP=omp_get_max_threads()
      OMP_FLAG = KEEP(399).EQ.1 .AND. NOMP .GE.2 .AND. SLAVEF.EQ.1
!$OMP PARALLEL PRIVATE( K8, I, DEST,
!$OMP&          T4MASTER, T4_MASTER_CONCERNED,
!$OMP&          INIV2, NCAND, IROW_GRID, JCOL_GRID, IPOSROOT, JPOSROOT,
!$OMP&          ILOCROOT, JLOCROOT,
!$OMP&          TYPE_NODE, TYPESPLIT, MASTER_NODE,
!$OMP&          IA8, ISHIFT, IIW8, IS18, IS8, IAS8, VAL,
!$OMP&          IARR, ISTEP, ISEND, JSEND,
!$OMP&          IOLD, JOLD, IOMP, DOIT, P2, NOMP_P, OMP_FLAG_P )
!$OMP& REDUCTION(+:NSEND8, NLOCAL8, ARROW_ROOT) IF (OMP_FLAG)
      IOMP=0
!$    IOMP=omp_get_thread_num()
      NOMP_P=1
!$    NOMP_P=omp_get_num_threads()
      OMP_FLAG_P = .FALSE.
!$    OMP_FLAG_P = OMP_FLAG .AND. NOMP_P .GT. 1
      IF (OMP_FLAG_P) THEN
        IF ( NOMP_P .GE. 16 ) THEN
          NOMP_P=16
          P2 = 4
        ELSE IF (NOMP_P.GE.8) THEN
          NOMP_P=8
          P2 = 3
        ELSE IF (NOMP_P.GE.4) THEN
          NOMP_P=4
          P2 = 2
        ELSE IF (NOMP_P.GE.2) THEN
          NOMP_P=2
          P2 = 1
        ENDIF
      ELSE
        NOMP_P = 1
        P2 = 0
      ENDIF
      IF ( IOMP .LT. NOMP_P ) THEN
       DO K8 = 1_8, NZ_loc8
        IF ( SLAVEF .GT. 1 ) THEN
!$OMP     MASTER
          KPROBE = KPROBE + 1
          IF ( KPROBE .eq. FREQPROBE ) THEN
            KPROBE = 0
            CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     &                       FLAG, STATUS, IERR )
            IF ( FLAG ) THEN
              MSGSOU = STATUS( MPI_SOURCE )
              CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, 
     &                   MPI_INTEGER,
     &                   MSGSOU, ARR_INT, COMM, STATUS, IERR )
              CALL MPI_RECV( BUFRECR(1), NBRECORDS,
     &                   MPI_DOUBLE_PRECISION,
     &                   MSGSOU, ARR_REAL, COMM, STATUS, IERR )
              CALL DMUMPS_DIST_TREAT_RECV_BUF(
     &               BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     &               KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &               A, LA,
     &               END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &               PTRAIW, PTRARW, PERM, STEP,
     &               INTARR, LINTARR, DBLARR, LDBLARR
     &               )
            END IF
          END IF
!$OMP     END MASTER
        ENDIF
        IOLD = id%IRN_loc(K8)
        JOLD = id%JCN_loc(K8)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     &                 .OR.(JOLD.LT.1) ) THEN
          CYCLE
        ENDIF
        IF (OMP_FLAG_P) THEN
          IF (IOLD.EQ.JOLD) THEN
            IARR  = IOLD
          ELSE IF (PERM(IOLD).LT.PERM(JOLD)) THEN
            IARR  = IOLD
          ELSE
            IARR  = JOLD
          ENDIF
          DOIT = ( IOMP .EQ. ibits(IARR, P2-1, P2))
        ELSE
          DOIT = .TRUE.
        ENDIF
        IF (DOIT) THEN
          IF (IOLD.EQ.JOLD) THEN
            ISEND = IOLD
            JSEND = IOLD
            IARR  = IOLD
          ELSE IF (PERM(IOLD).LT.PERM(JOLD)) THEN
            IARR  = IOLD
            IF ( KEEP(50) .NE. 0 ) THEN
              ISEND = -IOLD
            ELSE
              ISEND = IOLD
            ENDIF
            JSEND = JOLD
          ELSE
            IARR  = JOLD
            ISEND = -JOLD
            JSEND = IOLD
          ENDIF
          ISTEP = abs(STEP(IARR))
          CALL MUMPS_TYPEANDPROCNODE( TYPE_NODE, MASTER_NODE,
     &    PROCNODE_STEPS(ISTEP), KEEP(199) )
          T4_MASTER_CONCERNED = .FALSE.
          T4MASTER            = -9999
          VAL = id%A_loc(K8)
          IF ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) THEN
            VAL = VAL * id%ROWSCA(IOLD)*id%COLSCA(JOLD)
          ENDIF
          IF ( TYPE_NODE .eq. 1 ) THEN
            DEST = MASTER_NODE
            IF (DEST.EQ.MYID) THEN
              NLOCAL8 = NLOCAL8 + 1_8
              IF (ISEND.EQ.JSEND) THEN
                IA8 = PTRARW(ISEND)
                DBLARR(IA8) = DBLARR(IA8) + VAL
              ELSE IF (ISEND.GE.0) THEN 
                IS18         = PTRAIW(IARR)
                ISHIFT       = INTARR(IS18) + IW4(IARR,2)
                INTARR(IS18+ISHIFT+2)       = JSEND
                DBLARR(PTRARW(IARR)+ISHIFT) = VAL
                IW4(IARR,2)  = IW4(IARR,2) - 1
              ELSE 
                ISHIFT = IW4(IARR,1)
                INTARR(PTRAIW(IARR)+ISHIFT+2)  = JSEND
                DBLARR(PTRARW(IARR)+ISHIFT) = VAL
                IW4(IARR,1)  = IW4(IARR,1) - 1
                IF ( IW4(IARR,1) .EQ. 0
     &             .AND. STEP(IARR) > 0 ) THEN
                  CALL DMUMPS_QUICK_SORT_ARROWHEADS( N, PERM,
     &               INTARR( PTRAIW(IARR) + 3 ),
     &               DBLARR( PTRARW(IARR) + 1 ),
     &               INTARR( PTRAIW(IARR) ), 1,
     &               INTARR( PTRAIW(IARR) ) )
                END IF
              ENDIF
              CYCLE
            ENDIF
          ELSE IF ( TYPE_NODE .eq. 2 ) THEN
            IF ( ISEND .LT. 0 ) THEN
              DEST = -1
            ELSE
              DEST = MASTER_NODE
            END IF
            INIV2         = ISTEP_TO_INIV2(ISTEP)
            IF ( KEEP(79) .GT. 0) THEN
              TYPESPLIT  = MUMPS_TYPESPLIT( PROCNODE_STEPS(ISTEP),
     &                                      KEEP(199) )
              IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
                T4_MASTER_CONCERNED = .TRUE.
                T4MASTER=CANDIDATES(CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
              ENDIF
            ENDIF
          ELSE 
            ARROW_ROOT = ARROW_ROOT + 1
            IF (EARLYT3ROOTINS) THEN
              IF ( ISEND < 0 ) THEN
                IPOSROOT = root%RG2L_ROW(JSEND)
                JPOSROOT = root%RG2L_ROW(IARR )
              ELSE
                IPOSROOT = root%RG2L_ROW(IARR )
                JPOSROOT = root%RG2L_ROW(JSEND)
              END IF
              IROW_GRID = mod( ( IPOSROOT-1 )/root%MBLOCK, root%NPROW )
              JCOL_GRID = mod( ( JPOSROOT-1 )/root%NBLOCK, root%NPCOL )
              DEST = IROW_GRID * root%NPCOL + JCOL_GRID
            ELSE
              DEST = -2
            ENDIF
            IF ( OMP_FLAG_P ) THEN
              IF ( EARLYT3ROOTINS ) THEN
                     ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                            ( root%MBLOCK * root%NPROW ) )
     &                          + mod( IPOSROOT - 1, root%MBLOCK ) + 1
                     JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                            ( root%NBLOCK * root%NPCOL ) )
     &                          + mod( JPOSROOT - 1, root%NBLOCK ) + 1
                     IF (KEEP(60)==0) THEN
                       A( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &                   + int(ILOCROOT-1,8)) =  A( PTR_ROOT
     &                   + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &                   + int(ILOCROOT - 1,8) )
     &                 + VAL
                     ELSE
                       root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                                 * int(root%SCHUR_LLD,8)
     &                                 + int(ILOCROOT,8) )
     &                 = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                                 * int(root%SCHUR_LLD,8)
     &                                 + int(ILOCROOT,8))
     &                 + VAL
                     ENDIF
              ELSE
                    IF (ISEND.EQ.JSEND) THEN
                      IA8 = PTRARW(ISEND)
                      DBLARR(IA8) = DBLARR(IA8) + VAL
                    ELSE IF (ISEND.GE.0) THEN 
                      IS18         = PTRAIW(IARR)
                      ISHIFT       = INTARR(IS18) + IW4(IARR,2)
                      IW4(IARR,2)  = IW4(IARR,2) - 1
                      IIW8         = IS18 + ISHIFT + 2
                      INTARR(IIW8) = JSEND
                      IS8          = PTRARW(IARR)
                      IAS8         = IS8 + ISHIFT
                      DBLARR(IAS8) = VAL
                    ELSE 
                      IS8          = PTRAIW(IARR)+IW4(IARR,1)+2
                      INTARR(IS8)  = JSEND
                      IAS8         = PTRARW(IARR)+IW4(IARR,1)
                      IW4(IARR,1)  = IW4(IARR,1) - 1
                      DBLARR(IAS8) = VAL
                      IF ( IW4(IARR,1) .EQ. 0
     &                   .AND. STEP(IARR) > 0 ) THEN
                        CALL DMUMPS_QUICK_SORT_ARROWHEADS( N, PERM,
     &                     INTARR( PTRAIW(IARR) + 3 ),
     &                     DBLARR( PTRARW(IARR) + 1 ),
     &                     INTARR( PTRAIW(IARR) ), 1,
     &                     INTARR( PTRAIW(IARR) ) )
                      END IF
                    ENDIF
              ENDIF
              CYCLE
            ENDIF
          END IF
          IF (DEST .eq. -1) THEN
            NLOCAL8 = NLOCAL8 + 1_8
            NSEND8  = NSEND8 + int(SLAVEF -1,8)
          ELSE IF (DEST .EQ. -2) THEN
            NLOCAL8 = NLOCAL8 + 1_8
            NSEND8  = NSEND8 + int(SLAVEF -1,8)
          ELSE
            IF (DEST .eq.MYID ) THEN
              NLOCAL8 = NLOCAL8 + 1_8
            ELSE
              NSEND8 = NSEND8 + 1_8
            ENDIF
          ENDIF
          IF ( DEST.EQ.-1) THEN
            INIV2 = ISTEP_TO_INIV2(ISTEP)
            NCAND = CANDIDATES(SLAVEF+1,INIV2)
            IF (KEEP(79) .GT. 0) THEN
              DO I=1, SLAVEF
                DEST=CANDIDATES(I,INIV2)
                IF (DEST.LT.0) EXIT 
                IF (I.EQ.NCAND+1) CYCLE
                CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &       BUFI, BUFR, BUFRECI, BUFRECR,
     &       NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &       SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &       N, PTRAIW, PTRARW, PERM, STEP, END_MSG_2_RECV,
     &       PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N,
     &       IW4(1,1), root, KEEP,KEEP8 )
              ENDDO
            ELSE
              DO I=1, NCAND
                DEST=CANDIDATES(I,INIV2)
                CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &          BUFI, BUFR, BUFRECI, BUFRECR,
     &          NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &          SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &          N, PTRAIW, PTRARW, PERM, STEP, END_MSG_2_RECV,
     &          PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N,
     &          IW4(1,1), root, KEEP,KEEP8 )
              ENDDO
            ENDIF
           DEST=MASTER_NODE
           CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &     BUFI, BUFR, BUFRECI, BUFRECR,
     &     NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &     SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &     N, PTRAIW, PTRARW, PERM, STEP, END_MSG_2_RECV,
     &     PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N,
     &     IW4(1,1), root, KEEP,KEEP8 )
           IF (T4_MASTER_CONCERNED) THEN
            DEST = T4MASTER
            CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &      BUFI, BUFR, BUFRECI, BUFRECR,
     &      NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &      SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &      N, PTRAIW, PTRARW, PERM, STEP, END_MSG_2_RECV,
     &      PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N,
     &      IW4(1,1), root, KEEP,KEEP8 )
           ENDIF
          ELSE IF (DEST .GE. 0) THEN
           CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &     BUFI, BUFR, BUFRECI, BUFRECR,
     &     NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &     SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &     N, PTRAIW, PTRARW, PERM, STEP, END_MSG_2_RECV,
     &     PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N,
     &     IW4(1,1), root, KEEP,KEEP8 )
           IF (T4_MASTER_CONCERNED) THEN
            DEST = T4MASTER
            CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &      BUFI, BUFR, BUFRECI, BUFRECR,
     &      NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &      SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &      N, PTRAIW, PTRARW, PERM, STEP, END_MSG_2_RECV,
     &      PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N,
     &      IW4(1,1), root, KEEP,KEEP8 )
           ENDIF
          ELSE IF (DEST .EQ. -2) THEN
            DO I = 0, SLAVEF-1
              DEST=I
              CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &        BUFI, BUFR, BUFRECI, BUFRECR,
     &        NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &        SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &        N, PTRAIW, PTRARW, PERM, STEP, END_MSG_2_RECV,
     &        PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N,
     &        IW4(1,1), root, KEEP, KEEP8 )
            ENDDO
          ENDIF
        ENDIF 
       END DO
      ENDIF
!$OMP END PARALLEL
      DEST = -3
      CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &BUFI, BUFR, BUFRECI, BUFRECR,
     &NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &N, PTRAIW, PTRARW, PERM, STEP, END_MSG_2_RECV,
     &PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, 
     &IW4(1,1), root, KEEP,KEEP8 )
      DO WHILE ( END_MSG_2_RECV .NE. 0 )
        CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, MPI_INTEGER,
     &                 MPI_ANY_SOURCE, ARR_INT, COMM, STATUS, IERR )
        MSGSOU = STATUS( MPI_SOURCE )
        CALL MPI_RECV( BUFRECR(1), NBRECORDS, MPI_DOUBLE_PRECISION,
     &                 MSGSOU, ARR_REAL, COMM, STATUS, IERR )
        CALL DMUMPS_DIST_TREAT_RECV_BUF(
     &           BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     &           KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &           A, LA,
     &           END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &           PTRAIW, PTRARW, PERM, STEP,
     &           INTARR, LINTARR, DBLARR, LDBLARR
     &           )
      END DO
      DO I = 1, SLAVEF
        IF ( SEND_ACTIVE( I ) ) THEN
          CALL MPI_WAIT( IREQI( I ), STATUS, IERR )
          CALL MPI_WAIT( IREQR( I ), STATUS, IERR )
        END IF
      END DO
      KEEP(49) = ARROW_ROOT
 100  CONTINUE
      IF (ALLOCATED(IW4))     DEALLOCATE( IW4 )
      IF (ALLOCATED(BUFI))    DEALLOCATE( BUFI )
      IF (ALLOCATED(BUFR))    DEALLOCATE( BUFR )
      IF (ALLOCATED(BUFRECI)) DEALLOCATE( BUFRECI )
      IF (ALLOCATED(BUFRECR)) DEALLOCATE( BUFRECR )
      IF (ALLOCATED(IACT))    DEALLOCATE( IACT )
      IF (ALLOCATED(IREQI))   DEALLOCATE( IREQI )
      IF (ALLOCATED(IREQR))   DEALLOCATE( IREQR )
      IF (ALLOCATED(SEND_ACTIVE)) DEALLOCATE( SEND_ACTIVE )
      RETURN
      END SUBROUTINE DMUMPS_REDISTRIBUTION
      SUBROUTINE DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &  BUFI, BUFR, BUFRECI, BUFRECR,
     &  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR, N,
     &  PTRAIW, PTRARW, PERM, STEP, END_MSG_2_RECV,
     &  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4, root,
     &  KEEP,KEEP8 )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ISEND, JSEND, DEST, NBRECORDS, SLAVEF, COMM, MYID, N
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER END_MSG_2_RECV, LOCAL_M, LOCAL_N
      INTEGER(8) :: LINTARR, LDBLARR
      INTEGER(8) :: LA, PTR_ROOT
      INTEGER BUFI( NBRECORDS * 2 + 1, 2, SLAVEF )
      INTEGER BUFRECI( NBRECORDS * 2 + 1 )
      INTEGER IREQI(SLAVEF), IREQR(SLAVEF), IACT(SLAVEF)
      INTEGER IW4( N, 2 )
      INTEGER(8) PTRAIW( N ), PTRARW( N )
      INTEGER PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER INTARR( LINTARR )
      DOUBLE PRECISION DBLARR( LDBLARR ), A( LA )
      LOGICAL SEND_ACTIVE(SLAVEF)
      DOUBLE PRECISION BUFR( NBRECORDS, 2, SLAVEF )
      DOUBLE PRECISION BUFRECR( NBRECORDS )
      DOUBLE PRECISION VAL
      INTEGER ISLAVE, IBEG, IEND, NBREC, IREQ
      INTEGER TAILLE_SEND_I, TAILLE_SEND_R, MSGSOU
      LOGICAL FLAG, SEND_LOCAL
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      IF ( DEST .eq. -3 ) THEN
        IBEG = 1
        IEND = SLAVEF
      ELSE
        IBEG = DEST + 1
        IEND = DEST + 1
      END IF
      SEND_LOCAL = .FALSE.
      DO ISLAVE = IBEG, IEND
        NBREC = BUFI(1,IACT(ISLAVE),ISLAVE)
        IF ( DEST .eq. -3 ) THEN
          BUFI(1,IACT(ISLAVE),ISLAVE) = - NBREC
        END IF
        IF ( DEST .eq. -3 .or. NBREC + 1 > NBRECORDS ) THEN
          DO WHILE ( SEND_ACTIVE( ISLAVE ) )
            CALL MPI_TEST( IREQR( ISLAVE ), FLAG, STATUS, IERR )
            IF ( .NOT. FLAG ) THEN
                CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     &                           FLAG, STATUS, IERR )
                IF ( FLAG ) THEN
                  MSGSOU = STATUS(MPI_SOURCE)
                  CALL MPI_RECV( BUFRECI(1), 2*NBRECORDS+1,
     &                  MPI_INTEGER, MSGSOU, ARR_INT, COMM,
     &                  STATUS, IERR )
                  CALL MPI_RECV( BUFRECR(1), NBRECORDS,
     &                  MPI_DOUBLE_PRECISION, MSGSOU,
     &                  ARR_REAL, COMM, STATUS, IERR )
                  CALL DMUMPS_DIST_TREAT_RECV_BUF(
     &              BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     &              KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &              A, LA,
     &              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &              PTRAIW, PTRARW, PERM, STEP,
     &              INTARR, LINTARR, DBLARR, LDBLARR
     &              )
                END IF
            ELSE
                CALL MPI_WAIT( IREQI( ISLAVE ), STATUS, IERR )
                SEND_ACTIVE( ISLAVE ) = .FALSE.
            END IF
          END DO
          IF ( ISLAVE - 1 .ne. MYID ) THEN
            TAILLE_SEND_I = NBREC * 2 + 1
            TAILLE_SEND_R = NBREC
            CALL MPI_ISEND( BUFI(1, IACT(ISLAVE), ISLAVE ),
     &          TAILLE_SEND_I,
     &          MPI_INTEGER, ISLAVE - 1, ARR_INT, COMM,
     &          IREQI( ISLAVE ), IERR )
            CALL MPI_ISEND( BUFR(1, IACT(ISLAVE), ISLAVE ),
     &          TAILLE_SEND_R,
     &          MPI_DOUBLE_PRECISION, ISLAVE - 1, ARR_REAL, COMM,
     &          IREQR( ISLAVE ), IERR )
            SEND_ACTIVE( ISLAVE ) = .TRUE.
          ELSE
            SEND_LOCAL = .TRUE.
          END IF
          IACT( ISLAVE ) = 3 - IACT( ISLAVE )
          BUFI( 1, IACT( ISLAVE ), ISLAVE ) = 0
        END IF
        IF ( DEST .ne. -3 ) THEN
          IREQ = BUFI(1,IACT(ISLAVE),ISLAVE) + 1
          BUFI(1,IACT(ISLAVE),ISLAVE) = IREQ
          BUFI(IREQ*2,IACT(ISLAVE),ISLAVE)  = ISEND
          BUFI(IREQ*2+1,IACT(ISLAVE),ISLAVE) = JSEND
          BUFR(IREQ,IACT(ISLAVE),ISLAVE )    = VAL
        END IF
      END DO
      IF ( SEND_LOCAL ) THEN
            ISLAVE = MYID + 1
            CALL DMUMPS_DIST_TREAT_RECV_BUF(
     &              BUFI(1,3-IACT(ISLAVE),ISLAVE),
     &              BUFR(1,3-IACT(ISLAVE),ISLAVE),
     &              NBRECORDS, N, IW4(1,1),
     &              KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &              A, LA,
     &              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &              PTRAIW, PTRARW, PERM, STEP,
     &              INTARR, LINTARR, DBLARR, LDBLARR
     &              )
      END IF
      RETURN
      END SUBROUTINE DMUMPS_DIST_FILL_BUFFER
      SUBROUTINE DMUMPS_DIST_TREAT_RECV_BUF
     &           ( BUFI, BUFR, NBRECORDS, N, IW4,
     &             KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT, A, LA,
     &             END_MSG_2_RECV, MYID, PROCNODE_STEPS,
     &             SLAVEF,
     &             PTRAIW, PTRARW, PERM, STEP,
     &             INTARR, LINTARR, DBLARR, LDBLARR )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER NBRECORDS, N, MYID, SLAVEF
      INTEGER BUFI( NBRECORDS * 2 + 1 )
      DOUBLE PRECISION BUFR( NBRECORDS )
      INTEGER IW4( N, 2 )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER END_MSG_2_RECV
      INTEGER(8) :: PTRAIW( N ), PTRARW( N )
      INTEGER :: PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER(8), INTENT(IN) :: LINTARR, LDBLARR
      INTEGER INTARR( LINTARR )
      INTEGER LOCAL_M, LOCAL_N
      INTEGER(8) :: PTR_ROOT, LA
      DOUBLE PRECISION A( LA ), DBLARR( LDBLARR )
      INTEGER MUMPS_TYPENODE, MUMPS_PROCNODE
      EXTERNAL MUMPS_TYPENODE, MUMPS_PROCNODE
      INTEGER IREC, NB_REC, NODE_TYPE, IPROC
      INTEGER IPOSROOT, JPOSROOT, ILOCROOT, JLOCROOT
      INTEGER(8) :: IA8, IS18, IIW8, IS8, IAS8
      INTEGER ISHIFT, IARR, JARR
      INTEGER TAILLE
      LOGICAL :: EARLYT3ROOTINS
      DOUBLE PRECISION VAL
      EARLYT3ROOTINS = KEEP(200) .EQ.0
      NB_REC = BUFI( 1 )
      IF ( NB_REC .LE. 0 ) THEN
        END_MSG_2_RECV = END_MSG_2_RECV - 1
        NB_REC = - NB_REC
      END IF
      IF ( NB_REC .eq. 0 ) GOTO 100
      DO IREC = 1, NB_REC
        IARR = BUFI( IREC * 2 )
        JARR = BUFI( IREC * 2 + 1 )
        VAL  = BUFR( IREC )
        NODE_TYPE = MUMPS_TYPENODE( 
     &              PROCNODE_STEPS(abs(STEP(abs( IARR )))),
     &              KEEP(199) )
        IF ( NODE_TYPE .eq. 3 .AND. EARLYT3ROOTINS ) THEN
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L_ROW( IARR )
            JPOSROOT = root%RG2L_COL( JARR )
          ELSE
            IPOSROOT = root%RG2L_ROW( JARR )
            JPOSROOT = root%RG2L_COL( -IARR )
          END IF
          ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
          JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
          IF (KEEP(60)==0) THEN
            A( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &        + int(ILOCROOT-1,8)) =  A( PTR_ROOT
     &        + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &        + int(ILOCROOT - 1,8) )
     &      + VAL
          ELSE
            root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                      * int(root%SCHUR_LLD,8)
     &                      + int(ILOCROOT,8) )
     &      = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                      * int(root%SCHUR_LLD,8)
     &                      + int(ILOCROOT,8))
     &      + VAL
          ENDIF
        ELSE IF (IARR.GE.0) THEN
          IF (IARR.EQ.JARR) THEN
            IA8 = PTRARW(IARR)
            DBLARR(IA8) = DBLARR(IA8) + VAL
          ELSE
            IS18         = PTRAIW(IARR)
            ISHIFT       = INTARR(IS18) + IW4(IARR,2)
            IW4(IARR,2)  = IW4(IARR,2) - 1
            IIW8         = IS18 + ISHIFT + 2
            INTARR(IIW8) = JARR
            IS8          = PTRARW(IARR)
            IAS8         = IS8 + ISHIFT
            DBLARR(IAS8) = VAL
          ENDIF
        ELSE
          IARR = -IARR
          IS8          = PTRAIW(IARR)+IW4(IARR,1)+2
          INTARR(IS8)  = JARR
          IAS8         = PTRARW(IARR)+IW4(IARR,1)
          IW4(IARR,1)  = IW4(IARR,1) - 1
          DBLARR(IAS8) = VAL
          IF ( IW4(IARR,1) .EQ. 0
     &         .AND. STEP(IARR) > 0 ) THEN
            IPROC = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(IARR)),
     &                            KEEP(199) )
            IF ( IPROC .EQ. MYID ) THEN
              TAILLE = INTARR( PTRAIW(IARR) )
              CALL DMUMPS_QUICK_SORT_ARROWHEADS( N, PERM,
     &           INTARR( PTRAIW(IARR) + 3 ),
     &           DBLARR( PTRARW(IARR) + 1 ),
     &           TAILLE, 1, TAILLE )
            ENDIF
          END IF
        ENDIF
      ENDDO
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_DIST_TREAT_RECV_BUF
