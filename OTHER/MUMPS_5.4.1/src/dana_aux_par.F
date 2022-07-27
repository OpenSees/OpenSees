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
      MODULE DMUMPS_PARALLEL_ANALYSIS
      USE DMUMPS_STRUC_DEF
      USE MUMPS_MEMORY_MOD
      USE MUMPS_ANA_ORD_WRAPPERS
      INCLUDE 'mpif.h'
      PUBLIC DMUMPS_ANA_F_PAR
      INTERFACE DMUMPS_ANA_F_PAR
      MODULE PROCEDURE DMUMPS_ANA_F_PAR
      END INTERFACE
      PRIVATE
      TYPE ORD_TYPE
      INTEGER           :: CBLKNBR, N
      INTEGER, POINTER  :: PERMTAB(:) => null()
      INTEGER, POINTER  :: PERITAB(:) => null()
      INTEGER, POINTER  :: RANGTAB(:) => null()
      INTEGER, POINTER  :: TREETAB(:) => null()
      INTEGER, POINTER  :: BROTHER(:) => null()
      INTEGER, POINTER  :: SON(:) => null()
      INTEGER, POINTER  :: NW(:) => null()
      INTEGER, POINTER  :: FIRST(:) => null()
      INTEGER, POINTER  :: LAST(:) => null()
      INTEGER, POINTER  :: TOPNODES(:) => null()
      INTEGER           :: COMM, COMM_NODES, NPROCS, NSLAVES, MYID
      INTEGER           :: TOPSTRAT, SUBSTRAT, ORDTOOL, TOPVARS
      LOGICAL           :: IDO
      END TYPE ORD_TYPE
      TYPE GRAPH_TYPE
      INTEGER(8)        :: NZ_LOC
      INTEGER           :: N, COMM
      INTEGER, POINTER  :: IRN_LOC(:) => null()
      INTEGER, POINTER  :: JCN_LOC(:) => null()
      END TYPE GRAPH_TYPE
      TYPE ARRPNT
      INTEGER, POINTER :: BUF(:) => null()
      END TYPE ARRPNT
      INTEGER    :: MP, MPG, LP, NRL, TOPROWS
      INTEGER(8) :: MEMCNT, MAXMEM
      LOGICAL    :: PROK, PROKG, LPOK
      CONTAINS
      SUBROUTINE DMUMPS_ANA_F_PAR(id, WORK1, WORK2, NFSIZ, FILS,
     &     FRERE)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC)   :: id
      INTEGER, TARGET  :: WORK1(:), WORK2(:)
      INTEGER          :: NFSIZ(:), FILS(:), FRERE(:)
      TYPE(ORD_TYPE)       :: ord
      INTEGER, POINTER     :: IPE(:), NV(:),
     &     NE(:), NA(:), NODE(:),
     &     ND(:), SUBORD(:), NAMALG(:),
     &     IPS(:), CUMUL(:),
     &     SAVEIRN(:), SAVEJCN(:)
      INTEGER              :: MYID, NPROCS, IERR, NEMIN, LDIAG
      LOGICAL              :: SPLITROOT
      INTEGER(8), PARAMETER :: K79REF=12000000_8
      INTEGER, PARAMETER :: LIDUMMY = 1
      INTEGER :: IDUMMY(1) 
      DOUBLE PRECISION      :: TIMEB
      nullify(IPE, NV, NE, NA, NODE, ND, SUBORD, NAMALG, IPS,
     &     CUMUL, SAVEIRN, SAVEJCN)
      CALL MPI_COMM_RANK (id%COMM, MYID, IERR)
      CALL MPI_COMM_SIZE (id%COMM, NPROCS, IERR)
      LP    = id%ICNTL(1)
      MP    = id%ICNTL(2)
      MPG   = id%ICNTL(3)
      PROK  = (MP.GT.0)
      PROKG = (MPG.GT.0) .AND. (MYID .EQ. 0)
      LPOK  = (LP.GT.0) .AND. (id%ICNTL(4).GE.1)
      LDIAG = id%ICNTL(4)
      ord%PERMTAB => WORK1(1        : id%N)
      ord%PERITAB => WORK1(id%N+1   : 2*id%N)
      ord%TREETAB => WORK1(2*id%N+1 : 3*id%N)
      IF(id%KEEP(54) .NE. 3) THEN
         IF(MYID.EQ.0) THEN
            SAVEIRN    => id%IRN_loc
            SAVEJCN    => id%JCN_loc
            id%IRN_loc => id%IRN
            id%JCN_loc => id%JCN
            id%KEEP8(29) = id%KEEP8(28)
         ELSE
            id%KEEP8(29)=0_8
         END IF
      END IF
      MAXMEM=0
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      CALL DMUMPS_SET_PAR_ORD(id, ord)
      id%INFOG(7) = id%KEEP(245)
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      IF (PROKG) CALL MUMPS_SECDEB( TIMEB )
      CALL DMUMPS_DO_PAR_ORD(id, ord, WORK2)
      IF (PROKG) THEN
       CALL MUMPS_SECFIN( TIMEB )
       WRITE(*,'(" ELAPSED time in parallel ordering       =",F12.4)') 
     &   TIMEB
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      IF(id%MYID .EQ. 0) THEN
         CALL MUMPS_REALLOC(IPE, id%N, id%INFO, LP, FORCE=.FALSE.,
     &        COPY=.FALSE., STRING='', 
     &        MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(NV, id%N, id%INFO, LP,
     &        MEMCNT=MEMCNT, ERRCODE=-7)
         IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      END IF
      ord%SUBSTRAT = 0
      ord%TOPSTRAT = 0
      CALL DMUMPS_PARSYMFACT(id, ord, IPE, NV, WORK2)
      IF(id%KEEP(54) .NE. 3) THEN
         IF(MYID.EQ.0) THEN
            id%IRN_loc => SAVEIRN
            id%JCN_loc => SAVEJCN
         END IF
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      NULLIFY(ord%PERMTAB)
      NULLIFY(ord%PERITAB)
      NULLIFY(ord%TREETAB)
      CALL MUMPS_IDEALLOC(ord%FIRST, ord%LAST, MEMCNT=MEMCNT)
      IF (MYID .EQ. 0) THEN
         IPS => WORK1(1:id%N)
         NE     => WORK1(id%N+1   : 2*id%N)
         NA     => WORK1(2*id%N+1 : 3*id%N)
         NODE   => WORK2(1        : id%N  )
         ND     => WORK2(id%N+1   : 2*id%N)
         SUBORD => WORK2(2*id%N+1 : 3*id%N)
         NAMALG => WORK2(3*id%N+1 : 4*id%N)
         CALL MUMPS_REALLOC(CUMUL, id%N, id%INFO, LP,
     &        STRING='CUMUL', MEMCNT=MEMCNT, ERRCODE=-7)
         IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
         NEMIN  = id%KEEP(1)
         CALL DMUMPS_ANA_LNEW(id%N, IPE(1), NV(1), IPS(1), NE(1),
     &        NA(1), NFSIZ(1), NODE(1), id%INFOG(6), FILS(1), FRERE(1),
     &        ND(1), NEMIN, SUBORD(1), id%KEEP(60), id%KEEP(20),
     &        id%KEEP(38), NAMALG(1), id%KEEP(104), CUMUL(1),
     &        id%KEEP(50), id%ICNTL(13), id%KEEP(37), id%NSLAVES,
     &        id%KEEP(250).EQ.1, .FALSE., IDUMMY, LIDUMMY)
         CALL MUMPS_DEALLOC(CUMUL, NV, IPE, MEMCNT=MEMCNT)
         CALL DMUMPS_ANA_M(NE(1), ND(1), id%INFOG(6), id%INFOG(5),
     &        id%KEEP(2), id%KEEP(50), id%KEEP8(101), id%KEEP(108),
     &        id%KEEP(5), id%KEEP(6), id%KEEP(226), id%KEEP(253))
         IF ( id%KEEP(53) .NE. 0 ) THEN
            CALL MUMPS_MAKE1ROOT(id%N, FRERE(1), FILS(1), NFSIZ(1),
     &           id%KEEP(20))
         END IF
         IF (  (id%KEEP(48) == 4 .AND. id%KEEP8(21).GT.0_8)
     &        .OR.
     &        (id%KEEP (48)==5 .AND. id%KEEP8(21) .GT. 0_8 )
     &        .OR.
     &        (id%KEEP(24).NE.0.AND.id%KEEP8(21).GT.0_8) ) THEN 
            CALL DMUMPS_SET_K821_SURFACE(id%KEEP8(21), id%KEEP(2),
     &           id%KEEP(48), id%KEEP(50), id%NSLAVES)
         END IF
         IF ((id%KEEP(210).LT.0) .OR. (id%KEEP(210).GT.2))
     &        id%KEEP(210)=0
         IF ((id%KEEP(210).EQ.0) .AND. (id%KEEP(201).GT.0))
     &        id%KEEP(210)=1    
         IF ((id%KEEP(210).EQ.0) .AND. (id%KEEP(201).EQ.0))
     &        id%KEEP(210)=2    
         IF (id%KEEP(210).EQ.2) id%KEEP8(79)=huge(id%KEEP8(79))
         IF ((id%KEEP(210).EQ.1) .AND. (id%KEEP8(79).LE.0_8)) THEN
           id%KEEP8(79)=K79REF * int(id%NSLAVES,8)
         ENDIF
         IF ( (id%KEEP(79).EQ.0).OR.(id%KEEP(79).EQ.2).OR.
     &        (id%KEEP(79).EQ.3).OR.(id%KEEP(79).EQ.5).OR.
     &        (id%KEEP(79).EQ.6)
     &   )  THEN
          IF (id%KEEP(210).EQ.1) THEN
            SPLITROOT = .FALSE. 
            IF ( id%KEEP(62).GE.1) THEN
               IDUMMY(1) = -1
               CALL DMUMPS_CUTNODES(id%N, FRERE(1), FILS(1),
     &              NFSIZ(1), IDUMMY, LIDUMMY, id%INFOG(6),
     &              id%NSLAVES, id%KEEP(1), id%KEEP8(1), SPLITROOT,
     &              MP, LDIAG, id%INFOG(1), id%INFOG(2))
               IF (id%INFOG(1).LT.0) RETURN
            ENDIF
          ENDIF
         ENDIF
         SPLITROOT = (((id%ICNTL(13).GT.0) .AND.
     &        (id%NSLAVES.GT.id%ICNTL(13))) .OR.
     &        (id%ICNTL(13).EQ.-1)) .AND. (id%KEEP(60).EQ.0)
         IF (SPLITROOT) THEN
            IDUMMY(1) = -1
            CALL DMUMPS_CUTNODES(id%N, FRERE(1), FILS(1), NFSIZ(1),
     &           IDUMMY, LIDUMMY, id%INFOG(6), 
     &           id%NSLAVES, id%KEEP(1), id%KEEP8(1),
     &           SPLITROOT, MP, LDIAG, id%INFOG(1), id%INFOG(2))
            IF (id%INFOG(1).LT.0) RETURN
         ENDIF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_ANA_F_PAR
      SUBROUTINE DMUMPS_SET_PAR_ORD(id, ord)
      TYPE(DMUMPS_STRUC)  :: id
      TYPE(ORD_TYPE)      :: ord
      INTEGER  :: IERR
#if defined(parmetis) || defined(parmetis3)
      INTEGER  :: I, COLOR, BASE, WORKERS
      LOGICAL  :: IDO
#endif
      IF(id%MYID .EQ. 0) id%KEEP(245) = id%ICNTL(29)
      CALL MPI_BCAST( id%KEEP(245), 1,
     &     MPI_INTEGER, 0, id%COMM, IERR )
      IF ((id%KEEP(245) .LT. 0) .OR. (id%KEEP(245) .GT. 2)) THEN
         id%KEEP(245) = 0
      END IF
      IF (id%KEEP(245) .EQ. 0) THEN
#if defined(ptscotch)
         IF(id%NSLAVES .LT. 2) THEN
            IF(PROKG) WRITE(MPG,'("Warning: older versions
     &of PT-SCOTCH require at least 2 processors.")')
         END IF
         ord%ORDTOOL    = 1
         ord%TOPSTRAT   = 0
         ord%SUBSTRAT   = 0
         ord%COMM       = id%COMM
         ord%COMM_NODES = id%COMM_NODES
         ord%NPROCS     = id%NPROCS
         ord%NSLAVES    = id%NSLAVES
         ord%MYID       = id%MYID
         ord%IDO        = (id%MYID .GE. 1) .OR. (id%KEEP(46) .EQ. 1)
         id%KEEP(245)   = 1
         IF(PROKG) WRITE(MPG,
     &           '("Parallel ordering tool set to PT-SCOTCH.")')
         RETURN
#endif
#if defined(parmetis) || defined(parmetis3)
         IF(id%N.LE.100) THEN
            WORKERS = 2
         ELSE
            WORKERS = min(id%NSLAVES,id%N/16)
         END IF
         I=1
         DO
            IF (I .GT. WORKERS) EXIT
            ord%NSLAVES = I
            I = I*2
         END DO
         BASE = id%NPROCS-id%NSLAVES
         ord%NPROCS  = ord%NSLAVES + BASE
         IDO = (id%MYID .GE. BASE) .AND.
     &        (id%MYID .LE. BASE+ord%NSLAVES-1)
         ord%IDO = IDO
         IF ( IDO ) THEN
            COLOR = 1
         ELSE
            COLOR = MPI_UNDEFINED
         END IF
         CALL MPI_COMM_SPLIT( id%COMM, COLOR, 0, 
     &        ord%COMM_NODES, IERR )
         ord%ORDTOOL  = 2
         ord%TOPSTRAT = 0
         ord%SUBSTRAT = 0
         ord%MYID     = id%MYID
         IF(PROKG) WRITE(MPG,
     &        '("Parallel ordering tool set to ParMETIS.")')
         id%KEEP(245) = 2
         RETURN
#endif         
         id%INFO(1)  = -38
         id%INFOG(1) = -38
         IF(id%MYID .EQ.0 ) THEN
            WRITE(LP,
     &           '("No parallel ordering tools available.")')
            WRITE(LP,
     &           '("Please install PT-SCOTCH or ParMETIS.")')
         END IF
         RETURN
      ELSE IF (id%KEEP(245) .EQ. 1) THEN
#if defined(ptscotch)
         IF(id%NSLAVES .LT. 2) THEN
            IF(PROKG) WRITE(MPG,'("Warning: older versions
     &of PT-SCOTCH require at least 2 processors.")')
         END IF
         ord%ORDTOOL    = 1
         ord%TOPSTRAT   = 0
         ord%SUBSTRAT   = 0
         ord%COMM       = id%COMM
         ord%COMM_NODES = id%COMM_NODES
         ord%NPROCS     = id%NPROCS
         ord%NSLAVES    = id%NSLAVES
         ord%MYID       = id%MYID
         ord%IDO        = (id%MYID .GE. 1) .OR. (id%KEEP(46) .EQ. 1)
         IF(PROKG) WRITE(MPG,
     &        '(" Using PT-SCOTCH for parallel ordering")')
         RETURN
#else
         id%INFOG(1) = -38
         id%INFO(1)  = -38
         IF(id%MYID .EQ.0 ) WRITE(LP,
     &        '(" PT-SCOTCH not available")')
         RETURN
#endif
      ELSE IF (id%KEEP(245) .EQ. 2) THEN
#if defined(parmetis) || defined(parmetis3)
         IF(id%N.LE.100) THEN
            WORKERS = 2
         ELSE
            WORKERS = min(id%NSLAVES,id%N/16)
         END IF
         I=1
         DO
            IF (I .GT. WORKERS) EXIT
            ord%NSLAVES = I
            I = I*2
         END DO
         BASE = id%NPROCS-id%NSLAVES
         ord%NPROCS  = ord%NSLAVES + BASE
         IDO = (id%MYID .GE. BASE) .AND.
     &        (id%MYID .LE. BASE+ord%NSLAVES-1)
         ord%IDO = IDO
         IF ( IDO ) THEN
            COLOR   = 1
         ELSE
            COLOR = MPI_UNDEFINED
         END IF
         CALL MPI_COMM_SPLIT( id%COMM, COLOR, 0, ord%COMM_NODES,
     &        IERR )
         ord%ORDTOOL  = 2
         ord%TOPSTRAT = 0
         ord%SUBSTRAT = 0
         ord%MYID     = id%MYID
         IF(PROKG) WRITE(MPG,
     &        '(" Using ParMETIS for parallel ordering")')
         RETURN
#else
         id%INFOG(1) = -38
         id%INFO(1)  = -38
         IF(id%MYID .EQ.0 ) WRITE(LP,
     &        '(" ParMETIS not available.")')
         RETURN
#endif
      END IF
      END SUBROUTINE DMUMPS_SET_PAR_ORD
      SUBROUTINE DMUMPS_DO_PAR_ORD(id, ord, WORK)
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC)            :: id
      TYPE(ORD_TYPE)                :: ord
      INTEGER                       :: WORK(:)
#if defined(parmetis) || defined(parmetis3)
      INTEGER                       :: IERR
#endif
      IF (ord%ORDTOOL .EQ. 1) THEN
#if defined(ptscotch)
         CALL DMUMPS_PTSCOTCH_ORD(id, ord, WORK)
#else
         id%INFOG(1) = -38
         id%INFO(1)  = -38
         WRITE(LP,*)'PT-SCOTCH not available. Aborting...'
         CALL MUMPS_ABORT()
#endif
      ELSE IF (ord%ORDTOOL .EQ. 2) THEN
#if defined(parmetis) || defined(parmetis3)
         CALL DMUMPS_PARMETIS_ORD(id, ord, WORK)
         if(ord%IDO) CALL MPI_COMM_FREE(ord%COMM_NODES, IERR)
#else
         id%INFOG(1) = -38
         id%INFO(1)  = -38
         WRITE(LP,*)'ParMETIS not available. Aborting...'
         CALL MUMPS_ABORT()
#endif
      END IF
      RETURN
      END SUBROUTINE DMUMPS_DO_PAR_ORD
#if defined(parmetis) || defined(parmetis3)
      SUBROUTINE DMUMPS_PARMETIS_ORD(id, ord, WORK)
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC)  :: id
      TYPE(ORD_TYPE)      :: ord
      INTEGER, TARGET    :: WORK(:)
      INTEGER             :: I, MYID, NPROCS, IERR, BASE, METIS_IDX_SIZE
      INTEGER, POINTER    :: FIRST(:),
     &     LAST(:), SWORK(:)
      INTEGER             :: BASEVAL, VERTLOCNBR,
     &     OPTIONS(10)
      INTEGER(8), POINTER :: VERTLOCTAB(:)
      INTEGER, POINTER    :: EDGELOCTAB(:), RCVCNTS(:)
      INTEGER(8)          :: EDGELOCNBR
      INTEGER, POINTER    :: SIZES(:), ORDER(:)
      nullify(FIRST, LAST, SWORK, VERTLOCTAB, EDGELOCTAB, RCVCNTS,
     &      SIZES, ORDER)
      CALL MPI_COMM_RANK (id%COMM, MYID, IERR)
      CALL MPI_COMM_SIZE (id%COMM, NPROCS, IERR)
      IERR=0
      IF(size(WORK) .LT. id%N*3) THEN
         WRITE(LP,
     &        '("Insufficient workspace inside DMUMPS_PARMETIS_ORD")')
         CALL MUMPS_ABORT()
      END IF
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      BASEVAL = 1
      BASE    = id%NPROCS-id%NSLAVES
      CALL MUMPS_REALLOC(FIRST, NPROCS+1, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(LAST, NPROCS+1, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      CALL DMUMPS_GRAPH_DIST(id, ord, FIRST,
     &     LAST, BASE, NPROCS, WORK(1: 2*id%N), TYPE=2)
      VERTLOCNBR = LAST(MYID+1)-FIRST(MYID+1) + 1
      CALL MUMPS_I8REALLOC(VERTLOCTAB, VERTLOCNBR+1, id%INFO,
     &        LP, STRING='VERTLOCTAB', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      SWORK => WORK(id%N+1:3*id%N)
      CALL DMUMPS_BUILD_DIST_GRAPH(id, FIRST, LAST, VERTLOCTAB,
     &     EDGELOCTAB, SWORK)
      IF(id%INFO(1).LT.0) RETURN
      EDGELOCNBR = VERTLOCTAB(VERTLOCNBR+1)-1_8
      OPTIONS(:) = 0
      ORDER => WORK(1:id%N)
      CALL MUMPS_REALLOC(SIZES, 2*ord%NSLAVES, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      IF(ord%IDO) THEN
         CALL MUMPS_METIS_IDXSIZE(METIS_IDX_SIZE) 
         IF (METIS_IDX_SIZE.EQ.32) THEN
           IF (id%KEEP(10).EQ.1) THEN
            id%INFO(1) = -52
            id%INFO(2) = 1
           ELSE
            CALL MUMPS_PARMETIS_MIXEDto32(id, BASE, VERTLOCNBR, FIRST,
     &           VERTLOCTAB, EDGELOCTAB, BASEVAL, OPTIONS, ORDER,
     &           SIZES, ord%COMM_NODES, IERR)
           ENDIF
         ELSE IF (METIS_IDX_SIZE.EQ.64) THEN
            CALL MUMPS_PARMETIS_MIXEDto64
     &           (id, BASE, VERTLOCNBR, FIRST,
     &           VERTLOCTAB, EDGELOCTAB, BASEVAL, OPTIONS, ORDER,
     &           SIZES, ord%COMM_NODES, IERR)
         ELSE
              WRITE(*,*)
     &        "Internal error in PARMETIS wrappers, METIS_IDX_SIZE=",
     &        METIS_IDX_SIZE
              CALL MUMPS_ABORT()
         END IF
      END IF
      CALL MUMPS_IDEALLOC(EDGELOCTAB, MEMCNT=MEMCNT)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      CALL MUMPS_I8DEALLOC(VERTLOCTAB)
      IF(IERR.GT.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 20
      CALL MPI_BCAST(SIZES(1), 2*ord%NSLAVES, MPI_INTEGER,
     &     BASE, id%COMM, IERR)
      ord%CBLKNBR = 2*ord%NSLAVES-1
      CALL MUMPS_REALLOC(RCVCNTS, id%NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      DO I=1, id%NPROCS
         RCVCNTS(I) = max(LAST(I)-FIRST(I)+1,0)
      END DO
      FIRST = FIRST-1
      IF(FIRST(1) .LT. 0) THEN
         FIRST(1)   = 0
      END IF
      CALL MPI_ALLGATHERV ( ORDER(1), VERTLOCNBR, MPI_INTEGER,
     &     ord%PERMTAB(1),
     &     RCVCNTS(1), FIRST(1), MPI_INTEGER, id%COMM, IERR )
      DO I=1, id%N
         ord%PERITAB(ord%PERMTAB(I)) = I
      END DO
      CALL MUMPS_REALLOC(ord%RANGTAB, 2*ord%NSLAVES, id%INFO,
     &     LP, STRING='RANGTAB', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      CALL DMUMPS_BUILD_TREETAB(ord%TREETAB, ord%RANGTAB,
     &     SIZES, ord%CBLKNBR)
      CALL MUMPS_DEALLOC(SIZES, FIRST, LAST,
     &     RCVCNTS, MEMCNT=MEMCNT)
      CALL MUMPS_REALLOC(ord%SON, ord%CBLKNBR, id%INFO,
     &     LP, STRING='SON', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ord%BROTHER, ord%CBLKNBR, id%INFO,
     &     LP, STRING='BROTHER', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ord%NW, ord%CBLKNBR, id%INFO,
     &     LP, STRING='NW', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      CALL DMUMPS_BUILD_TREE(ord)
      ord%N = id%N
      ord%COMM = id%COMM
      RETURN
 20   CONTINUE
      CALL MUMPS_DEALLOC(FIRST      , MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(LAST       , MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(SIZES      , MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(ord%RANGTAB, MEMCNT=MEMCNT)
      RETURN
      END SUBROUTINE DMUMPS_PARMETIS_ORD
#endif
#if defined(ptscotch)
      SUBROUTINE DMUMPS_PTSCOTCH_ORD(id, ord, WORK)
      IMPLICIT NONE
      INCLUDE 'ptscotchf.h'
      TYPE(DMUMPS_STRUC)            :: id
      TYPE(ORD_TYPE)                :: ord
      INTEGER, TARGET               :: WORK(:)
      INTEGER                       :: MYID, NPROCS, IERR
      INTEGER, POINTER              :: FIRST(:),
     &     LAST(:), SWORK(:)
      INTEGER                       :: BASEVAL, VERTLOCNBR,
     &     BASE, SCOTCH_INT_SIZE
      INTEGER(8)                    :: EDGELOCNBR
      INTEGER(8), POINTER           :: VERTLOCTAB(:)
      INTEGER, POINTER              :: EDGELOCTAB(:)
      nullify(FIRST, LAST, SWORK, VERTLOCTAB, EDGELOCTAB)
      IF (size(WORK) .LT. id%N*3) THEN
         WRITE(LP,
     &      '("Insufficient workspace inside DMUMPS_PTSCOTCH_ORD")')
         CALL MUMPS_ABORT()
      END IF
      CALL MPI_BARRIER(id%COMM, IERR)
      CALL MPI_COMM_RANK (id%COMM, MYID, IERR)
      CALL MPI_COMM_SIZE (id%COMM, NPROCS, IERR)
      BASE     = id%NPROCS-id%NSLAVES
      BASEVAL  = 1
      CALL MUMPS_REALLOC(FIRST, NPROCS+1, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(LAST, NPROCS+1, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      CALL DMUMPS_GRAPH_DIST(id, ord, FIRST,
     &     LAST, BASE, NPROCS, WORK(1: 2*id%N), TYPE=2)
      VERTLOCNBR = LAST(MYID+1)-FIRST(MYID+1) + 1
      CALL MUMPS_I8REALLOC(VERTLOCTAB, VERTLOCNBR+1, id%INFO,
     &        LP, STRING='VERTLOCTAB', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      SWORK => WORK(id%N+1:3*id%N)
      CALL DMUMPS_BUILD_DIST_GRAPH(id, FIRST, LAST, VERTLOCTAB,
     &     EDGELOCTAB, SWORK)
      IF(id%INFO(1).LT.0) RETURN
      EDGELOCNBR = VERTLOCTAB(VERTLOCNBR+1)-1_8
      CALL MUMPS_REALLOC(ord%PERMTAB, id%N, id%INFO,
     &     LP, STRING='PERMTAB', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ord%PERITAB, id%N, id%INFO,
     &     LP, STRING='PERITAB', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ord%RANGTAB, id%N+1, id%INFO,
     &     LP, STRING='RANGTAB', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ord%TREETAB, id%N, id%INFO,
     &     LP, STRING='TREETAB', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      IF(ord%IDO) THEN
         CALL MUMPS_SCOTCH_INTSIZE(SCOTCH_INT_SIZE)
         IF(SCOTCH_INT_SIZE.EQ.32) THEN
           IF (id%KEEP(10).EQ.1) THEN
            id%INFO(1) = -52
            id%INFO(2) = 2
           ELSE
            CALL MUMPS_PTSCOTCH_MIXEDto32(id, ord,
     &           BASEVAL, 
     &           VERTLOCNBR, VERTLOCTAB,
     &           EDGELOCNBR, EDGELOCTAB,
     &           IERR)
          ENDIF
         ELSE
            CALL MUMPS_PTSCOTCH_MIXEDto64(id, ord,
     &           BASEVAL, 
     &           VERTLOCNBR, VERTLOCTAB,
     &           EDGELOCNBR, EDGELOCTAB,
     &           IERR)
         END IF
      END IF
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 11
      CALL  MPI_BCAST (ord%CBLKNBR, 1,      MPI_INTEGER,
     &     BASE, id%COMM, IERR)
      CALL  MPI_BCAST (ord%PERMTAB(1), id%N,   MPI_INTEGER,
     &     BASE, id%COMM, IERR)
      CALL  MPI_BCAST (ord%PERITAB(1), id%N,   MPI_INTEGER,
     &     BASE, id%COMM, IERR)
      CALL  MPI_BCAST (ord%RANGTAB(1), id%N+1, MPI_INTEGER,
     &     BASE, id%COMM, IERR)
      CALL  MPI_BCAST (ord%TREETAB(1), id%N,   MPI_INTEGER,
     &     BASE, id%COMM, IERR)
      CALL MUMPS_REALLOC(ord%SON, ord%CBLKNBR, id%INFO,
     &     LP, STRING='SON', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ord%BROTHER, ord%CBLKNBR, id%INFO,
     &     LP, STRING='BROTHER', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ord%NW, ord%CBLKNBR, id%INFO,
     &     LP, STRING='NW', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL DMUMPS_BUILD_TREE(ord)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      ord%N = id%N
      ord%COMM = id%COMM
      CALL MUMPS_DEALLOC(EDGELOCTAB, MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(FIRST     , MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(LAST      , MEMCNT=MEMCNT)
      CALL MUMPS_I8DEALLOC(VERTLOCTAB, MEMCNT=MEMCNT)
      RETURN
 11   CONTINUE
      CALL MUMPS_DEALLOC(FIRST      , MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(LAST       , MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(ord%RANGTAB, MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(EDGELOCTAB, MEMCNT=MEMCNT)
      RETURN
      END SUBROUTINE DMUMPS_PTSCOTCH_ORD
#endif
      FUNCTION DMUMPS_STOP_DESCENT(id, ord, NACTIVE, ANODE, RPROC,
     &     ALIST, LIST, PEAKMEM, NNODES, CHECKMEM)
      IMPLICIT NONE
      LOGICAL              :: DMUMPS_STOP_DESCENT
      INTEGER              :: NACTIVE, RPROC, ANODE, PEAKMEM, NNODES
      INTEGER              :: ALIST(NNODES), LIST(NNODES)
      TYPE(ORD_TYPE)       :: ord
      TYPE(DMUMPS_STRUC)   :: id
      LOGICAL, OPTIONAL    :: CHECKMEM
      INTEGER              :: IPEAKMEM, BIG, MAX_NROWS, MIN_NROWS
      INTEGER              :: TOPROWS, NRL, HOSTMEM, SUBMEM
      INTEGER              :: I, NZ_ROW, WEIGHT
      LOGICAL              :: ICHECKMEM
      INTEGER              :: NZ4
      IF(present(CHECKMEM)) THEN
         ICHECKMEM = CHECKMEM
      ELSE
         ICHECKMEM = .FALSE.
      END IF
      DMUMPS_STOP_DESCENT = .FALSE.
      IF(NACTIVE .GE. RPROC) THEN
         DMUMPS_STOP_DESCENT = .TRUE.
         RETURN
      END IF
      IF(NACTIVE .EQ. 0) THEN
         DMUMPS_STOP_DESCENT = .TRUE.
         RETURN
      END IF
      IF(.NOT. ICHECKMEM) RETURN
      BIG = ALIST(NACTIVE)
      IF(NACTIVE .GT. 1) THEN
         MAX_NROWS = ord%NW(ALIST(NACTIVE-1))
         MIN_NROWS = ord%NW(ALIST(1))
      ELSE
         MAX_NROWS = 0
         MIN_NROWS = id%N
      END IF
      DO I=1, ANODE
         WEIGHT = ord%NW(LIST(I))
         IF(WEIGHT .GT. MAX_NROWS) MAX_NROWS = WEIGHT
         IF(WEIGHT .LT. MIN_NROWS) MIN_NROWS = WEIGHT
      END DO
      I = ord%SON(BIG)
      DO
         WEIGHT = ord%NW(I)
         IF(WEIGHT .GT. MAX_NROWS) MAX_NROWS = WEIGHT
         IF(WEIGHT .LT. MIN_NROWS) MIN_NROWS = WEIGHT
         IF(ord%BROTHER(I) .EQ. -1) EXIT
         I = ord%BROTHER(I)
      END DO
      TOPROWS = ord%TOPNODES(2)+ord%RANGTAB(BIG+1)-ord%RANGTAB(BIG)
      SUBMEM  = 7 *id%N 
      HOSTMEM = 12*id%N 
      NZ4=int(id%KEEP8(28))
      NZ_ROW = 2*(NZ4/id%N) 
      IF(id%KEEP(46) .EQ. 0) THEN
         NRL = 0
      ELSE
         NRL = MIN_NROWS
      END IF
      HOSTMEM = HOSTMEM + 2*TOPROWS*NZ_ROW
      HOSTMEM = HOSTMEM +NRL
      HOSTMEM = HOSTMEM + max(NRL,TOPROWS)*(NZ_ROW+2)
      HOSTMEM = HOSTMEM + 6*max(NRL,TOPROWS)
      HOSTMEM = HOSTMEM + 3*TOPROWS
      NRL = MAX_NROWS
      SUBMEM = SUBMEM +NRL
      SUBMEM = SUBMEM + NRL*(NZ_ROW+2)
      SUBMEM = SUBMEM + 6*NRL
      IPEAKMEM = max(HOSTMEM, SUBMEM)
      IF((IPEAKMEM .GT. PEAKMEM) .AND.
     &     (PEAKMEM .NE. 0)) THEN
         DMUMPS_STOP_DESCENT = .TRUE.
         RETURN
      ELSE
         DMUMPS_STOP_DESCENT = .FALSE.
         PEAKMEM = IPEAKMEM
         RETURN
      END IF
      END FUNCTION DMUMPS_STOP_DESCENT
      FUNCTION DMUMPS_CNT_KIDS(NODE, ord)
      IMPLICIT NONE
      INTEGER :: DMUMPS_CNT_KIDS
      INTEGER :: NODE
      TYPE(ORD_TYPE) :: ord
      INTEGER :: CURR
      DMUMPS_CNT_KIDS = 0
      IF(ord%SON(NODE) .EQ. -1) THEN
         RETURN
      ELSE
         DMUMPS_CNT_KIDS = 1
         CURR = ord%SON(NODE)
         DO
            IF(ord%BROTHER(CURR) .NE. -1) THEN
               DMUMPS_CNT_KIDS = DMUMPS_CNT_KIDS+1
               CURR = ord%BROTHER(CURR)
            ELSE
               EXIT
            END IF
         END DO
      END IF
      RETURN
      END FUNCTION DMUMPS_CNT_KIDS
      SUBROUTINE DMUMPS_GET_SUBTREES(ord, id)
      IMPLICIT NONE
      TYPE(ORD_TYPE)     :: ord
      TYPE(DMUMPS_STRUC) :: id
      INTEGER, ALLOCATABLE :: ALIST(:), AWEIGHTS(:), LIST(:), WORK(:)
      INTEGER  :: NNODES, BIG, CURR, ND, NACTIVE, RPROC, ANODE, BASE, I,
     &     NK, PEAKMEM, allocok
      LOGICAL  :: SD
      NNODES = ord%NSLAVES
      CALL MUMPS_REALLOC(ord%TOPNODES, 2*max(NNODES,2), id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ord%FIRST, id%NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ord%LAST, id%NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      ALLOCATE(ALIST(NNODES), AWEIGHTS(NNODES), LIST(NNODES),
     &     WORK(0:NNODES+1), stat=allocok)
      IF(allocok.GT.0) THEN
         id%INFO(1)=-13
         id%INFO(2)=4*NNODES+2
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1), id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GO TO 90
      NACTIVE = 0
      DO I=1, ord%CBLKNBR
         IF (ord%TREETAB(I).EQ.-1) THEN
            NACTIVE = NACTIVE+1
            IF(NACTIVE.LE.NNODES) THEN
               ALIST(NACTIVE) = I
               AWEIGHTS(NACTIVE) = ord%NW(I)
            END IF
         END IF
      END DO
      IF((ord%CBLKNBR .EQ. 1) .OR.
     &   (NACTIVE.GT.NNODES) .OR.
     &   ( NNODES .LT. DMUMPS_CNT_KIDS(ord%CBLKNBR, ord) )) THEN
         ord%TOPNODES(1) = 1
         ord%TOPNODES(2) = ord%RANGTAB(ord%CBLKNBR+1) - ord%RANGTAB(1)
         ord%TOPNODES(3) = ord%RANGTAB(1)
         ord%TOPNODES(4) = ord%RANGTAB(ord%CBLKNBR+1)-1
         ord%FIRST = 0
         ord%LAST  = -1
         RETURN
      END IF
      CALL DMUMPS_MERGESORT(NACTIVE, AWEIGHTS(1:NACTIVE),
     &     WORK(0:NACTIVE+1))
      CALL DMUMPS_MERGESWAP(NACTIVE, WORK(0:NACTIVE+1),
     &     AWEIGHTS(1:NACTIVE), 
     &     ALIST(1:NACTIVE))
      RPROC       = NNODES   
      ANODE       = 0
      PEAKMEM     = 0
      ord%TOPNODES = 0
      DO
         IF(NACTIVE .EQ. 0) EXIT
         BIG = ALIST(NACTIVE)
         NK  = DMUMPS_CNT_KIDS(BIG, ord)
         IF((NK .GT. (RPROC-NACTIVE+1)) .OR. (NK .EQ. 0)) THEN 
            ANODE       = ANODE+1
            LIST(ANODE) = BIG
            NACTIVE     = NACTIVE-1
            RPROC       = RPROC-1
            CYCLE
         END IF
         SD = DMUMPS_STOP_DESCENT(id, ord, NACTIVE, ANODE,
     &        RPROC, ALIST, LIST, PEAKMEM, NNODES, CHECKMEM=.TRUE.)
         IF ( SD ) 
     &        THEN
            IF(NACTIVE.GT.0) THEN
               LIST(ANODE+1:ANODE+NACTIVE) = ALIST(1:NACTIVE)
               ANODE = ANODE+NACTIVE
            END IF
            EXIT
         END IF
         ord%TOPNODES(1) = ord%TOPNODES(1)+1
         ord%TOPNODES(2) = ord%TOPNODES(2) +
     &        ord%RANGTAB(BIG+1) - ord%RANGTAB(BIG)
         ord%TOPNODES(2+2*(ord%TOPNODES(1)-1)+1) = ord%RANGTAB(BIG)
         ord%TOPNODES(2+2*(ord%TOPNODES(1)-1)+2) = 
     &        ord%RANGTAB(BIG+1)-1
         CURR              = ord%SON(BIG)
         ALIST(NACTIVE)    = CURR
         AWEIGHTS(NACTIVE) = ord%NW(CURR)
         DO
            IF(ord%BROTHER(CURR) .EQ. -1) EXIT
            NACTIVE           = NACTIVE+1
            CURR              = ord%BROTHER(CURR)
            ALIST(NACTIVE)    = CURR
            AWEIGHTS(NACTIVE) = ord%NW(CURR)
         END DO
         CALL DMUMPS_MERGESORT(NACTIVE, AWEIGHTS(1:NACTIVE),
     &        WORK(0:NACTIVE+1))
         CALL DMUMPS_MERGESWAP(NACTIVE, WORK(0:NACTIVE+1),
     &        AWEIGHTS(1:NACTIVE), 
     &        ALIST(1:NACTIVE))
      END DO
      DO I=1, ANODE
         AWEIGHTS(I) = ord%NW(LIST(I))
      END DO
      CALL DMUMPS_MERGESORT(ANODE, AWEIGHTS(1:ANODE), WORK(0:ANODE+1))
      CALL DMUMPS_MERGESWAP(ANODE, WORK(0:ANODE+1), AWEIGHTS(1:ANODE), 
     &     ALIST(1:ANODE))
      IF (id%KEEP(46) .EQ. 1) THEN
         BASE = 0
      ELSE
         ord%FIRST(1) = 0
         ord%LAST(1)  = -1
         BASE = 1
      END IF
      DO I=1, ANODE
         CURR = LIST(I)
         ND = CURR
         IF(ord%SON(ND) .NE. -1) THEN
            ND = ord%SON(ND)
            DO
               IF((ord%SON(ND) .EQ. -1) .AND. 
     &              (ord%BROTHER(ND).EQ.-1)) THEN
                  EXIT
               ELSE IF(ord%BROTHER(ND) .EQ. -1) THEN
                  ND = ord%SON(ND)
               ELSE 
                  ND = ord%BROTHER(ND)
               END IF
            END DO
         END IF
         ord%FIRST(BASE+I) = ord%RANGTAB(ND)
         ord%LAST(BASE+I)  = ord%RANGTAB(CURR+1)-1
      END DO
      DO I=ANODE+1, id%NSLAVES
         ord%FIRST(BASE+I) = id%N+1
         ord%LAST(BASE+I) = id%N
      END DO            
      DEALLOCATE(LIST, ALIST, AWEIGHTS, WORK)
 90   continue
      RETURN
      END SUBROUTINE DMUMPS_GET_SUBTREES
      SUBROUTINE DMUMPS_PARSYMFACT(id, ord, GPE, GNV, WORK)  
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC)   :: id
      TYPE(ORD_TYPE)       :: ord
      INTEGER, POINTER     :: GPE(:), GNV(:)
      INTEGER, TARGET      :: WORK(:)
      TYPE(GRAPH_TYPE)     :: top_graph
      INTEGER(8), POINTER  :: IPE(:), IPET(:),
     &     BUF_PE1(:), BUF_PE2(:), TMP1(:)
      INTEGER, POINTER     :: PE(:),
     &     LENG(:), I_HALO_MAP(:)
      INTEGER, POINTER     :: NDENSE(:), LAST(:),
     &     DEGREE(:), W(:), PERM(:),
     &     LISTVAR_SCHUR(:), NEXT(:),
     &     HEAD(:), NV(:), ELEN(:),
     &     RCVCNT(:), LSTVAR(:)
      INTEGER, POINTER     :: MYLIST(:),
     &     LPERM(:),
     &     LIPERM(:),
     &     NVT(:), BUF_NV1(:),
     &     BUF_NV2(:), ROOTPERM(:),
     &     TMP2(:), BWORK(:), NCLIQUES(:)
      INTEGER              :: MYNCLIQUES, MYMAXVARS, ICLIQUES,
     &     TOTNCLIQUES
      INTEGER(8)           :: MYNVARS, TOTNVARS
      INTEGER(8), POINTER  :: LVARPT(:)
      INTEGER              :: HIDX, NCMPA, I, J, SIZE_SCHUR, MYID,
     &     NPROCS, IERR, NROWS_LOC, GLOB_IDX, TMP,
     &     NTVAR, TGSIZE, MAXS, RHANDPE,
     &     RHANDNV, RIDX, PROC, JOB, K
      INTEGER(8)           :: PFREES, PFS_SAVE, PELEN, PFREET, PFT_SAVE
      INTEGER              :: STATUSPE(MPI_STATUS_SIZE)
      INTEGER              :: STATUSNV(MPI_STATUS_SIZE)
      INTEGER              :: STATUSCLIQUES(MPI_STATUS_SIZE)
      INTEGER, PARAMETER   :: ITAG=30
      LOGICAL              :: AGG6
      INTEGER              :: THRESH
      nullify(PE, IPE, LENG, I_HALO_MAP, NCLIQUES)
      nullify(NDENSE, LAST, DEGREE, W, PERM, LISTVAR_SCHUR,
     &     NEXT, HEAD, NV, ELEN, RCVCNT, LSTVAR)
      nullify(MYLIST, LVARPT,
     &     LPERM, LIPERM, IPET, NVT, BUF_PE1, BUF_PE2,
     &     BUF_NV1, BUF_NV2, ROOTPERM, TMP1, TMP2, BWORK)
      CALL MPI_COMM_RANK (id%COMM, MYID, IERR)
      CALL MPI_COMM_SIZE (id%COMM, NPROCS, IERR)
      IF(size(WORK) .LT. 4*id%N) THEN
         WRITE(LP,*)'Insufficient workspace in DMUMPS_PARSYMFACT'
         CALL MUMPS_ABORT()
      ELSE
         HEAD => WORK(       1 :   id%N)
         ELEN => WORK(  id%N+1 : 2*id%N)
         LENG => WORK(2*id%N+1 : 3*id%N)
         PERM => WORK(3*id%N+1 : 4*id%N)
      END IF
      CALL DMUMPS_GET_SUBTREES(ord, id)
      CALL MUMPS_IDEALLOC(ord%SON, ord%BROTHER, ord%NW,
     &     ord%RANGTAB, MEMCNT=MEMCNT)
      NROWS_LOC = ord%LAST(MYID+1)-ord%FIRST(MYID+1)+1
      NRL = NROWS_LOC
      TOPROWS = ord%TOPNODES(2)
      BWORK => WORK(1 : 2*id%N)
      CALL DMUMPS_BUILD_LOC_GRAPH(id, ord, HIDX, IPE, PE, LENG,
     &     I_HALO_MAP, top_graph, BWORK)
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF(id%INFO(1).lt.0) RETURN
      TMP = id%N
      DO I=1, NPROCS
         TMP = TMP-(ord%LAST(I)-ord%FIRST(I)+1)
      END DO
      TMP = ceiling(dble(TMP)*1.10D0)
      IF(MYID .EQ. 0) THEN
         TMP = max(max(TMP, HIDX),1)
      ELSE
         TMP = max(HIDX,1)
      END IF
      SIZE_SCHUR = HIDX - NROWS_LOC
      CALL MUMPS_REALLOC(NDENSE, TMP, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(LAST, TMP, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(NEXT, TMP, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(DEGREE, TMP, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(W, TMP, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(NV, TMP, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(LISTVAR_SCHUR, max(SIZE_SCHUR,1), id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      DO I=1, SIZE_SCHUR
         LISTVAR_SCHUR(I) = NROWS_LOC+I
      END DO
      THRESH = -1
      AGG6   = .FALSE.
      PFREES = IPE(NROWS_LOC+1)
      PFS_SAVE = PFREES
      PELEN = PFREES-1 + 2_8*int(NROWS_LOC+ord%TOPNODES(2),8)
         DO I=1, HIDX
            PERM(I) = I
         END DO
         IF(SIZE_SCHUR.EQ.0) THEN
            JOB = 0
         ELSE
            JOB = 1
         END IF
         IF(HIDX .GT.0) CALL MUMPS_SYMQAMD_NEW(JOB, THRESH, NDENSE(1),
     &        HIDX, PELEN, IPE(1), PFREES, LENG(1), PE(1), NV(1), 
     &        ELEN(1), LAST(1), NCMPA, DEGREE(1), HEAD(1), NEXT(1), 
     &        W(1), PERM(1), LISTVAR_SCHUR(1), SIZE_SCHUR, AGG6)
      MYNCLIQUES = 0
      MYNVARS    = 0
      MYMAXVARS  = 0
      DO I=1, HIDX
         IF(IPE(I) .GT. 0) THEN
            MYMAXVARS  = MAX(MYMAXVARS,LENG(I))
            MYNVARS    = MYNVARS+LENG(I)
            MYNCLIQUES = MYNCLIQUES+1
         END IF
      END DO
      CALL MPI_REDUCE(MYNVARS, TOTNVARS, 1, MPI_INTEGER8,
     &     MPI_SUM, 0, id%COMM, IERR)
      CALL MUMPS_REALLOC(NCLIQUES, NPROCS, id%INFO,
     &     LP, STRING='NCLIQUES', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MPI_GATHER(MYNCLIQUES, 1, MPI_INTEGER, NCLIQUES(1), 1, 
     &     MPI_INTEGER, 0, id%COMM, IERR)
      IF(id%MYID.EQ.0) THEN
         TOTNCLIQUES = sum(NCLIQUES)
         CALL MUMPS_I8REALLOC(LVARPT, TOTNCLIQUES+1, id%INFO,
     &        LP, STRING='LVARPT', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_IREALLOC8(LSTVAR, TOTNVARS, id%INFO,
     &        LP, STRING='LSTVAR', MEMCNT=MEMCNT, ERRCODE=-7)
         LVARPT(1) = 1_8
         ICLIQUES  = 0
         DO I=1, HIDX
            IF(IPE(I) .GT. 0) THEN
               ICLIQUES = ICLIQUES+1
               LVARPT(ICLIQUES+1) = LVARPT(ICLIQUES)+LENG(I)
               DO J=0, LENG(I)-1
                  LSTVAR(LVARPT(ICLIQUES)+J) =
     &                 I_HALO_MAP(PE(IPE(I)+J)-NROWS_LOC)
               END DO
            END IF
         END DO
         DO PROC=1, NPROCS-1
            DO I=1, NCLIQUES(PROC+1)
               ICLIQUES = ICLIQUES+1
               CALL MPI_RECV(K, 1, MPI_INTEGER, PROC, ITAG, id%COMM,
     &              STATUSCLIQUES, IERR)
               LVARPT(ICLIQUES+1) = LVARPT(ICLIQUES)+K
               CALL MPI_RECV(LSTVAR(LVARPT(ICLIQUES)), K, MPI_INTEGER,
     &              PROC, ITAG, id%COMM, STATUSCLIQUES, IERR)
            END DO
         END DO
         LPERM => WORK(3*id%N+1 : 4*id%N)
         NTVAR   = ord%TOPNODES(2)
         CALL DMUMPS_MAKE_LOC_IDX(id, ord%TOPNODES, LPERM, LIPERM, ord)
         CALL DMUMPS_ASSEMBLE_TOP_GRAPH(id, ord%TOPNODES(2), LPERM,
     &        top_graph, TOTNCLIQUES, LSTVAR, LVARPT, IPET, PE,
     &        LENG, ELEN)
         TGSIZE = ord%TOPNODES(2)+TOTNCLIQUES
         PFREET = IPET(TGSIZE+1)
         PFT_SAVE = PFREET
         nullify(LPERM)
      ELSE
         CALL MUMPS_REALLOC(MYLIST, MYMAXVARS, id%INFO,
     &        LP, STRING='MYLIST', MEMCNT=MEMCNT, ERRCODE=-7)
         IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
         DO I=1, HIDX
            IF(IPE(I) .GT. 0) THEN
               DO J=1, LENG(I)
                  MYLIST(J) = I_HALO_MAP(PE(IPE(I)+J-1)-NROWS_LOC)
               END DO
               CALL MPI_SEND(LENG(I), 1, MPI_INTEGER, 0, ITAG,
     &              id%COMM, IERR)
               CALL MPI_SEND(MYLIST(1), LENG(I), MPI_INTEGER, 0, ITAG,
     &              id%COMM, IERR)
            END IF
         END DO
      END IF
      CALL MUMPS_IDEALLOC(top_graph%IRN_LOC,
     &     top_graph%JCN_LOC, ord%TOPNODES, MEMCNT=MEMCNT)
      IF(MYID .EQ. 0) THEN
         CALL MUMPS_IREALLOC8(PE, max(PFREET+int(TGSIZE,8),1_8),id%INFO,
     &        LP, COPY=.TRUE., STRING='J2:PE', MEMCNT=MEMCNT,
     &        ERRCODE=-7)
         CALL MUMPS_REALLOC(NDENSE, max(TGSIZE,1), id%INFO, LP,
     &        STRING='J2:NDENSE', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(NVT, max(TGSIZE,1), id%INFO, LP,
     &        STRING='J2:NVT', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(LAST, max(TGSIZE,1), id%INFO, LP,
     &        STRING='J2:LAST', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(DEGREE, max(TGSIZE,1), id%INFO, LP,
     &        STRING='J2:DEGREE', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(NEXT, max(TGSIZE,1), id%INFO, LP,
     &        STRING='J2:NEXT', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(W, max(TGSIZE,1), id%INFO, LP,
     &        STRING='J2:W', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(LISTVAR_SCHUR, max(TOTNCLIQUES,1), id%INFO,
     &        LP, STRING='J2:LVSCH', MEMCNT=MEMCNT, ERRCODE=-7)
         IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
         DO I=1, TOTNCLIQUES
            LISTVAR_SCHUR(I) = NTVAR+I
         END DO
         THRESH = -1
            CALL MUMPS_REALLOC(HEAD, max(TGSIZE,1), id%INFO,
     &        LP, STRING='J2:HEAD', MEMCNT=MEMCNT, ERRCODE=-7)
            CALL MUMPS_REALLOC(PERM, max(TGSIZE,1), id%INFO,
     &           LP, COPY=.TRUE., STRING='J2:PERM',
     &           MEMCNT=MEMCNT, ERRCODE=-7)
            IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
            DO I=1, TGSIZE
               PERM(I) = I
            END DO
            PELEN = max(PFREET+int(TGSIZE,8),1_8)
            IF(TGSIZE.GT.0) CALL MUMPS_SYMQAMD_NEW(2, -1, NDENSE(1),
     &           TGSIZE, PELEN, IPET(1), PFREET, LENG(1), PE(1), 
     &           NVT(1), ELEN(1), LAST(1), NCMPA, DEGREE(1), HEAD(1), 
     &           NEXT(1), W(1), PERM(1), LISTVAR_SCHUR(1), TOTNCLIQUES, 
     &           AGG6)
      END IF
      CALL MPI_BARRIER(id%COMM, IERR)
      CALL MPI_BARRIER(id%COMM, IERR)
      CALL MUMPS_DEALLOC(LISTVAR_SCHUR, MEMCNT=MEMCNT) 
      CALL MUMPS_DEALLOC(PE, MEMCNT=MEMCNT) 
      IF(MYID .EQ. 0) THEN
         MAXS = NROWS_LOC
         DO I=2, NPROCS
            IF((ord%LAST(I)-ord%FIRST(I)+1) .GT. MAXS)
     &           MAXS = (ord%LAST(I)-ord%FIRST(I)+1)
         END DO
         CALL MUMPS_I8REALLOC(BUF_PE1, max(MAXS,1), id%INFO,
     &        LP, STRING='BUF_PE1', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_I8REALLOC(BUF_PE2, max(MAXS,1), id%INFO,
     &        LP, STRING='BUF_PE2', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(BUF_NV1, max(MAXS,1), id%INFO,
     &        LP, STRING='BUF_NV1', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(BUF_NV2, max(MAXS,1), id%INFO,
     &        LP, STRING='BUF_NV2', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(GPE, id%N, id%INFO,
     &        LP, STRING='GPE', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(GNV, id%N, id%INFO,
     &        LP, STRING='GNV', MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_REALLOC(ROOTPERM, TOTNCLIQUES, id%INFO,
     &        LP, STRING='ROOTPERM', MEMCNT=MEMCNT, ERRCODE=-7)
         IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
         RIDX = 0
         TMP1    => BUF_PE1
         TMP2    => BUF_NV1
         NULLIFY(BUF_PE1, BUF_NV1)
         BUF_PE1 => IPE
         BUF_NV1 => NV
         DO PROC=0, NPROCS-2
            CALL MPI_IRECV(BUF_PE2(1), ord%LAST(PROC+2)-
     &           ord%FIRST(PROC+2)+1, MPI_INTEGER8, PROC+1, PROC+1,
     &           id%COMM, RHANDPE, IERR)
            CALL MPI_IRECV(BUF_NV2(1), ord%LAST(PROC+2)-
     &           ord%FIRST(PROC+2)+1, MPI_INTEGER, PROC+1, PROC+1,
     &           id%COMM, RHANDNV, IERR)
            DO I=1, ord%LAST(PROC+1)-ord%FIRST(PROC+1)+1
               GLOB_IDX = ord%PERITAB(I+ord%FIRST(PROC+1)-1)
               IF(BUF_PE1(I) .GT. 0) THEN
                  RIDX=RIDX+1
                  ROOTPERM(RIDX) = GLOB_IDX
                  GNV(GLOB_IDX) = BUF_NV1(I)
               ELSE IF (BUF_PE1(I) .EQ. 0) THEN
                  GPE(GLOB_IDX) = 0
                  GNV(GLOB_IDX) = BUF_NV1(I)
               ELSE
                  GPE(GLOB_IDX) = -ord%PERITAB(-BUF_PE1(I)+
     &                 ord%FIRST(PROC+1)-1)
                  GNV(GLOB_IDX) = BUF_NV1(I)
               END IF
            END DO
            CALL MPI_WAIT(RHANDPE, STATUSPE, IERR)
            CALL MPI_WAIT(RHANDNV, STATUSNV, IERR)
            IF(PROC .NE. 0) THEN
               TMP1    => BUF_PE1
               TMP2    => BUF_NV1
            END IF
            BUF_PE1 => BUF_PE2
            BUF_NV1 => BUF_NV2
            NULLIFY(BUF_PE2, BUF_NV2)
            BUF_PE2 => TMP1
            BUF_NV2 => TMP2
            NULLIFY(TMP1, TMP2)
         END DO
         DO I=1, ord%LAST(PROC+1)-ord%FIRST(PROC+1)+1
            GLOB_IDX = ord%PERITAB(I+ord%FIRST(PROC+1)-1)
            IF(BUF_PE1(I) .GT. 0) THEN
               RIDX=RIDX+1
               ROOTPERM(RIDX) = GLOB_IDX
               GNV(GLOB_IDX) = BUF_NV1(I)
            ELSE IF (BUF_PE1(I) .EQ. 0) THEN
               GPE(GLOB_IDX) = 0
               GNV(GLOB_IDX) = BUF_NV1(I)
            ELSE
               GPE(GLOB_IDX) = -ord%PERITAB(-BUF_PE1(I)+
     &              ord%FIRST(PROC+1)-1)
               GNV(GLOB_IDX) = BUF_NV1(I)
            END IF
         END DO
         DO I=1, NTVAR
            GLOB_IDX = LIPERM(I)
            IF(IPET(I) .EQ. 0) THEN
               GPE(GLOB_IDX) = 0
               GNV(GLOB_IDX) = NVT(I)
            ELSE
               GPE(GLOB_IDX) = -LIPERM(-IPET(I))
               GNV(GLOB_IDX) = NVT(I)
            END IF
         END DO
         DO I=1, TOTNCLIQUES
            GLOB_IDX      = ROOTPERM(I)
            GPE(GLOB_IDX) = -LIPERM(-IPET(NTVAR+I))
         END DO
      ELSE
         CALL MPI_SEND(IPE(1), ord%LAST(MYID+1)-ord%FIRST(MYID+1)+1,
     &        MPI_INTEGER8, 0, MYID, id%COMM, IERR)
         CALL MPI_SEND(NV(1), ord%LAST(MYID+1)-ord%FIRST(MYID+1)+1,
     &        MPI_INTEGER, 0, MYID, id%COMM, IERR)
      END IF
      CALL MUMPS_DEALLOC(BUF_NV1, BUF_NV2, MEMCNT=MEMCNT)
      CALL MUMPS_I8DEALLOC(BUF_PE1, BUF_PE2, IPE, IPET,
     &     TMP1, LVARPT, MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(PE, I_HALO_MAP, NDENSE,
     &     LAST, DEGREE, MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(W, LISTVAR_SCHUR, NEXT,
     &     NV, MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(LSTVAR, NCLIQUES, MYLIST, 
     &     MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(LPERM, LIPERM, NVT, MEMCNT=MEMCNT)
      CALL MUMPS_DEALLOC(ROOTPERM, TMP2, MEMCNT=MEMCNT)
      NULLIFY(HEAD, ELEN, LENG, PERM, RCVCNT)
      RETURN
      END SUBROUTINE DMUMPS_PARSYMFACT
      SUBROUTINE DMUMPS_MAKE_LOC_IDX(id, TOPNODES, LPERM, LIPERM, ord)
      IMPLICIT NONE 
      TYPE(DMUMPS_STRUC)   :: id
      INTEGER, POINTER  :: TOPNODES(:), LPERM(:), LIPERM(:)
      TYPE(ORD_TYPE)    :: ord
      INTEGER           :: I, J, K, GIDX
      CALL MUMPS_REALLOC(LPERM , ord%N, id%INFO,
     &        LP, STRING='LIDX:LPERM', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(LIPERM, TOPNODES(2), id%INFO,
     &        LP, STRING='LIDX:LIPERM', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      LPERM = 0
      K = 1 
      DO I=TOPNODES(1), 1, -1
         DO J=TOPNODES(2*I+1), TOPNODES(2*I+2)
            GIDX        = ord%PERITAB(J) 
            LPERM(GIDX) = K
            LIPERM(K)   = GIDX
            K           = K+1
         END DO
      END DO
      RETURN
      END SUBROUTINE DMUMPS_MAKE_LOC_IDX
      SUBROUTINE DMUMPS_ASSEMBLE_TOP_GRAPH(id, NLOCVARS, LPERM,
     &     top_graph, NCLIQUES, LSTVAR, LVARPT, IPE, PE, LENG, ELEN)
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC)  :: id
      TYPE(GRAPH_TYPE)    :: top_graph
      INTEGER, POINTER    :: LPERM(:), LSTVAR(:), 
     &     PE(:), LENG(:), ELEN(:)
      INTEGER(8)          :: LVARPT(:)
      INTEGER             :: NCLIQUES
      INTEGER(8), POINTER :: IPE(:)
      INTEGER             :: I, IDX, NLOCVARS
      INTEGER(8)          :: INNZ, PNT, SAVEPNT
      CALL MUMPS_REALLOC(LENG, max(NLOCVARS+NCLIQUES,1)  , id%INFO,
     &        LP, STRING='ATG:LENG', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(ELEN, max(NLOCVARS+NCLIQUES,1)  , id%INFO,
     &        LP, STRING='ATG:ELEN', MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_I8REALLOC(IPE , NLOCVARS+NCLIQUES+1, id%INFO,
     &        LP, STRING='ATG:IPE', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      LENG = 0
      ELEN = 0
      DO INNZ=1, top_graph%NZ_LOC
         IF((LPERM(top_graph%JCN_LOC(INNZ)) .NE. 0) .AND.
     &        (top_graph%JCN_LOC(INNZ) .NE. top_graph%IRN_LOC(INNZ)))
     &        THEN
            LENG(LPERM(top_graph%IRN_LOC(INNZ))) =
     &           LENG(LPERM(top_graph%IRN_LOC(INNZ))) + 1
         END IF
      END DO
      DO I=1, NCLIQUES
         DO INNZ=LVARPT(I), LVARPT(I+1)-1
            ELEN(LPERM(LSTVAR(INNZ))) = ELEN(LPERM(LSTVAR(INNZ)))+1
            LENG(NLOCVARS+I) = LENG(NLOCVARS+I)+1
         END DO
      END DO
      IPE(1) = 1
      DO I=1, NLOCVARS+NCLIQUES
         IPE(I+1) = IPE(I)+int(LENG(I),8)+int(ELEN(I),8)
      END DO
      CALL MUMPS_IREALLOC8(PE, IPE(NLOCVARS+NCLIQUES+1)+
     &     int(NLOCVARS,8)+int(NCLIQUES,8),
     &     id%INFO, LP, STRING='ATG:PE', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      LENG = 0
      ELEN = 0
      DO I=1, NCLIQUES
         DO INNZ=LVARPT(I), LVARPT(I+1)-1
            IDX = LPERM(LSTVAR(INNZ))
            PE(IPE(IDX)+int(ELEN(IDX),8)) = NLOCVARS+I
            PE(IPE(NLOCVARS+I)+int(LENG(NLOCVARS+I),8)) = IDX
            ELEN(LPERM(LSTVAR(INNZ))) = ELEN(LPERM(LSTVAR(INNZ)))+1
            LENG(NLOCVARS+I) = LENG(NLOCVARS+I)+1
         end do
      end do
      DO INNZ=1, top_graph%NZ_LOC
         IF((LPERM(top_graph%JCN_LOC(INNZ)) .NE. 0) .AND.
     &        (top_graph%JCN_LOC(INNZ) .NE. top_graph%IRN_LOC(INNZ)))
     &        THEN
            PE(IPE(LPERM(top_graph%IRN_LOC(INNZ)))+
     &           ELEN(LPERM(top_graph%IRN_LOC(INNZ))) +
     &           LENG(LPERM(top_graph%IRN_LOC(INNZ)))) =
     &           LPERM(top_graph%JCN_LOC(INNZ))
            LENG(LPERM(top_graph%IRN_LOC(INNZ))) =
     &           LENG(LPERM(top_graph%IRN_LOC(INNZ))) + 1
         END IF
      END DO
      DO I=1, NLOCVARS+NCLIQUES
         LENG(I) = LENG(I)+ELEN(I)
      END DO
      SAVEPNT = 1
      PNT = 0
      LPERM(1:NLOCVARS+NCLIQUES) = 0
      DO I=1, NLOCVARS+NCLIQUES
         DO INNZ=IPE(I), IPE(I+1)-1
            IF(LPERM(PE(INNZ)) .EQ. I) THEN
               LENG(I) = LENG(I)-1
            ELSE
               LPERM(PE(INNZ)) = I 
               PNT = PNT+1
               PE(PNT) = PE(INNZ)
            END IF
         END DO
         IPE(I) = SAVEPNT
         SAVEPNT = PNT+1
      END DO
      IPE(NLOCVARS+NCLIQUES+1) = SAVEPNT
      RETURN
      END SUBROUTINE DMUMPS_ASSEMBLE_TOP_GRAPH
#if defined(parmetis) || defined(parmetis3)
      SUBROUTINE DMUMPS_BUILD_TREETAB(TREETAB, RANGTAB, SIZES, CBLKNBR)
      INTEGER, POINTER  :: TREETAB(:), RANGTAB(:), SIZES(:)
      INTEGER           :: CBLKNBR,allocok
      INTEGER           :: LCHILD, RCHILD, K, I
      INTEGER, POINTER  :: PERM(:)
      ALLOCATE(PERM(CBLKNBR),stat=allocok)
      if(allocok.GT.0) then
         write(*,*) "Allocation error of PERM in DMUMPS_BUILD_TREETAB"
         return
      endif
      TREETAB(CBLKNBR) = -1
      IF(CBLKNBR .EQ. 1) THEN
         DEALLOCATE(PERM)
         TREETAB(1) = -1
         RANGTAB(1) = 1
         RANGTAB(2)= SIZES(1)+1
         RETURN
      END IF
      LCHILD = CBLKNBR - (CBLKNBR+1)/2
      RCHILD = CBLKNBR-1
      K = 1
      PERM(CBLKNBR) = CBLKNBR
      PERM(LCHILD) = CBLKNBR+1 - (2*K+1)
      PERM(RCHILD) = CBLKNBR+1 - (2*K)
      TREETAB(RCHILD) = CBLKNBR
      TREETAB(LCHILD) = CBLKNBR
      IF(CBLKNBR .GT. 3) THEN
         CALL REC_TREETAB(TREETAB, PERM, (CBLKNBR-1)/2,
     &        LCHILD, CBLKNBR, 2*K+1)
         CALL REC_TREETAB(TREETAB, PERM, (CBLKNBR-1)/2,
     &        RCHILD, CBLKNBR, 2*K)
      END IF
      RANGTAB(1)=1
      DO I=1, CBLKNBR
         RANGTAB(I+1) = RANGTAB(I)+SIZES(PERM(I))
      END DO
      DEALLOCATE(PERM)
      RETURN
      CONTAINS
      RECURSIVE SUBROUTINE REC_TREETAB(TREETAB, PERM, SUBNODES,
     &     ROOTN, CBLKNBR, K)
      INTEGER, POINTER  :: TREETAB(:), PERM(:)
      INTEGER           :: SUBNODES, ROOTN, K, CBLKNBR
      INTEGER           :: LCHILD, RCHILD
      LCHILD = ROOTN - (SUBNODES+1)/2
      RCHILD = ROOTN-1
      PERM(LCHILD) = CBLKNBR+1 - (2*K+1)
      PERM(RCHILD) = CBLKNBR+1 - (2*K)
      TREETAB(RCHILD) = ROOTN
      TREETAB(LCHILD) = ROOTN
      IF(SUBNODES .GT. 3) THEN
         CALL REC_TREETAB(TREETAB, PERM, (SUBNODES-1)/2, LCHILD,
     &        CBLKNBR, 2*K+1)
         CALL REC_TREETAB(TREETAB, PERM, (SUBNODES-1)/2, RCHILD,
     &        CBLKNBR, 2*K)
      END IF
      END SUBROUTINE REC_TREETAB
      END SUBROUTINE DMUMPS_BUILD_TREETAB
#endif
#if defined(ptscotch) || defined(parmetis) || defined(parmetis3)
      SUBROUTINE DMUMPS_BUILD_DIST_GRAPH(id, FIRST, LAST, IPE,
     &     PE, WORK)
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC)      :: id
      INTEGER(8), POINTER     :: IPE(:)
      INTEGER, POINTER        :: FIRST(:), LAST(:), PE(:),
     &     WORK(:)
      INTEGER                 :: IERR, MYID, NPROCS
      INTEGER                 :: I, PROC, J, LOC_ROW
      INTEGER(8)              :: LOCNNZ, INNZ, NEW_LOCNNZ, OFFDIAG,
     &                           RCVPNT, PNT, SAVEPNT, DUPS, TOTDUPS
      INTEGER                 :: NROWS_LOC
      INTEGER                 :: STATUS(MPI_STATUS_SIZE)
      INTEGER, POINTER        :: MAPTAB(:), SDISPL(:)
      INTEGER(8), POINTER     :: MSGCNT(:), SNDCNT(:), RCVCNT(:)
      INTEGER, POINTER        :: RDISPL(:), BUFLEVEL(:),
     &     SIPES(:,:), LENG(:)
      INTEGER, POINTER        :: PCNT(:), TSENDI(:),
     &     TSENDJ(:), RCVBUF(:)
      TYPE(ARRPNT), POINTER   :: APNT(:)
      INTEGER                 :: BUFSIZE, SOURCE, MAXS, allocok
      INTEGER, PARAMETER      :: ITAG=30
      LOGICAL                 :: FLAG
      DOUBLE PRECISION        :: SYMMETRY
      INTEGER(KIND=8)         :: TLEN
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      INTEGER                 :: L
#endif
      nullify(MAPTAB, SNDCNT, RCVCNT, SDISPL)
      nullify(RDISPL, MSGCNT, SIPES, LENG, BUFLEVEL)
      nullify(PCNT, TSENDI, TSENDJ, RCVBUF, APNT)
      CALL MPI_COMM_RANK (id%COMM, MYID, IERR)
      CALL MPI_COMM_SIZE (id%COMM, NPROCS, IERR)
      IF(MUMPS_GETSIZE(WORK) .LT. id%N*2) THEN
         WRITE(LP,
     &        '("Insufficient workspace inside BUILD_SCOTCH_GRAPH")')
         CALL MUMPS_ABORT()
      END IF
      CALL MUMPS_I8REALLOC(SNDCNT, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_I8REALLOC(RCVCNT, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_I8REALLOC(MSGCNT, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(RDISPL, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      ALLOCATE(APNT(NPROCS), stat=allocok)
      IF(allocok.GT.0) THEN
         id%INFO(1)=-13
         id%INFO(2)=NPROCS
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1), id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GO TO 90
      SNDCNT = 0
      BUFSIZE = 1000
      BUFSIZE = id%KEEP(39)
      LOCNNZ = id%KEEP8(29)
      NROWS_LOC = LAST(MYID+1)-FIRST(MYID+1)+1
      MAPTAB => WORK(     1 :   id%N)
      LENG   => WORK(id%N+1 : 2*id%N)
      MAXS = 0
      DO I=1, NPROCS
         IF((LAST(I)-FIRST(I)+1) .GT. MAXS) THEN
            MAXS = LAST(I)-FIRST(I)+1
         END IF
         DO J=FIRST(I), LAST(I)
            MAPTAB(J) = I
         END DO
      END DO
      ALLOCATE(SIPES(max(1,MAXS), NPROCS), stat=allocok)
      IF(allocok.GT.0) THEN
         id%INFO(1)=-13
         id%INFO(2)=max(1,MAXS)*NPROCS
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1), id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GO TO 90
      OFFDIAG=0
      SIPES=0
      DO INNZ=1, LOCNNZ
         IF(id%IRN_loc(INNZ) .NE. id%JCN_loc(INNZ)) THEN
            OFFDIAG = OFFDIAG+1
            PROC = MAPTAB(id%IRN_loc(INNZ))
            LOC_ROW = id%IRN_loc(INNZ)-FIRST(PROC)+1
            SIPES(LOC_ROW, PROC) = SIPES(LOC_ROW, PROC)+1
            SNDCNT(PROC) = SNDCNT(PROC)+1
            PROC = MAPTAB(id%JCN_loc(INNZ))
            LOC_ROW = id%JCN_loc(INNZ)-FIRST(PROC)+1
            SIPES(LOC_ROW, PROC) = SIPES(LOC_ROW, PROC)+1
            SNDCNT(PROC) = SNDCNT(PROC)+1
         END IF
      END DO
      CALL MPI_ALLREDUCE (OFFDIAG, id%KEEP8(127), 1, MPI_INTEGER8,
     &     MPI_SUM, id%COMM, IERR)
      id%KEEP8(127) = id%KEEP8(127)+3*id%N
      id%KEEP8(126) = id%KEEP8(127)-2*id%N
      CALL MPI_ALLTOALL(SNDCNT(1), 1, MPI_INTEGER8, RCVCNT(1), 1,
     &     MPI_INTEGER8, id%COMM, IERR)
      CALL MUMPS_I8DEALLOC(SNDCNT, MEMCNT=MEMCNT)
      RDISPL(:) = MAXS
      CALL MPI_REDUCE_SCATTER ( SIPES(1,1), LENG(1), RDISPL(1), 
     &     MPI_INTEGER, MPI_SUM, id%COMM, IERR )
      DEALLOCATE(SIPES)
      TLEN = 0_8
      IPE(1) = 1_8
      DO I=1, NROWS_LOC
        IPE(I+1) = IPE(I) + int(LENG(I),8)
         TLEN = TLEN+int(LENG(I),8)
      END DO
      CALL MUMPS_IREALLOC8(PE, max(IPE(NROWS_LOC+1)-1_8,1_8), id%INFO,
     &        LP, STRING='PE', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      LENG(:) = 0
      CALL MUMPS_REALLOC(BUFLEVEL, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL DMUMPS_SEND_BUF(APNT, PROC, NPROCS, BUFSIZE, IPE, PE, LENG,
     &     RCVBUF, MSGCNT, BUFLEVEL, id%COMM)
      NEW_LOCNNZ = 0
      DO I=1, NPROCS
         NEW_LOCNNZ = NEW_LOCNNZ+RCVCNT(I)
         MSGCNT(I) = RCVCNT(I)/int(BUFSIZE,8)
      END DO
      RCVPNT = 1
      BUFLEVEL = 0 
      DO INNZ=1, LOCNNZ
         IF(mod(INNZ,int(BUFSIZE,8)/10_8) .EQ. 0) THEN
            CALL MPI_IPROBE( MPI_ANY_SOURCE, ITAG, id%COMM,
     &           FLAG, STATUS, IERR )
            IF(FLAG) THEN
               SOURCE = STATUS(MPI_SOURCE)
               CALL MPI_RECV(RCVBUF(1), 2*BUFSIZE, MPI_INTEGER, SOURCE,
     &              ITAG, id%COMM, STATUS, IERR)
               CALL DMUMPS_ASSEMBLE_MSG(BUFSIZE, RCVBUF, IPE, PE, LENG)
               MSGCNT(SOURCE+1)=MSGCNT(SOURCE+1)-1
               RCVPNT = RCVPNT + BUFSIZE
            END IF
         END IF
         IF(id%IRN_loc(INNZ) .NE. id%JCN_loc(INNZ)) THEN
            PROC = MAPTAB(id%IRN_loc(INNZ))
            APNT(PROC)%BUF(2*BUFLEVEL(PROC)+1) = id%IRN_loc(INNZ)-
     &           FIRST(PROC)+1
            APNT(PROC)%BUF(2*BUFLEVEL(PROC)+2) = id%JCN_loc(INNZ)
            BUFLEVEL(PROC) = BUFLEVEL(PROC)+1
            IF(BUFLEVEL(PROC) .EQ. BUFSIZE) THEN
               CALL DMUMPS_SEND_BUF(APNT, PROC, NPROCS, BUFSIZE, IPE,
     &              PE, LENG, RCVBUF, MSGCNT, BUFLEVEL, id%COMM)
            END IF
            PROC = MAPTAB(id%JCN_loc(INNZ))
            APNT(PROC)%BUF(2*BUFLEVEL(PROC)+1) = id%JCN_loc(INNZ)-
     &           FIRST(PROC)+1
            APNT(PROC)%BUF(2*BUFLEVEL(PROC)+2) = id%IRN_loc(INNZ)
            BUFLEVEL(PROC) = BUFLEVEL(PROC)+1
            IF(BUFLEVEL(PROC) .EQ. BUFSIZE) THEN
               CALL DMUMPS_SEND_BUF(APNT, PROC, NPROCS, BUFSIZE, IPE,
     &              PE, LENG, RCVBUF, MSGCNT, BUFLEVEL, id%COMM)
            END IF
         END IF
      END DO
      CALL DMUMPS_SEND_BUF(APNT, -1, NPROCS, BUFSIZE, IPE, PE, LENG,
     &     RCVBUF, MSGCNT, BUFLEVEL, id%COMM)
      DUPS = 0
      PNT = 0
      SAVEPNT = 1
      MAPTAB = 0
      DO I=1, NROWS_LOC
         DO INNZ=IPE(I),IPE(I+1)-1
            IF(MAPTAB(PE(INNZ)) .EQ. I) THEN
               DUPS = DUPS+1
            ELSE
               MAPTAB(PE(INNZ)) = I 
               PNT = PNT+1
               PE(PNT) = PE(INNZ)
            END IF
         END DO
         IPE(I) = SAVEPNT
         SAVEPNT = PNT+1
      END DO
      CALL MPI_REDUCE( DUPS, TOTDUPS, 1, MPI_INTEGER8, MPI_SUM,
     &     0,  id%COMM, IERR )
      IF(MYID .EQ. 0) THEN
         SYMMETRY = dble(TOTDUPS)/(dble(id%KEEP8(28))-dble(id%N))
         SYMMETRY = min(SYMMETRY,1.0d0)
         IF(id%KEEP(50) .GE. 1) SYMMETRY = 1.d0
         IF(PROKG) WRITE(MPG,'(" Structural symmetry is:",i3,"%")')
     &        ceiling(SYMMETRY*100.d0)
         id%INFOG(8) = ceiling(SYMMETRY*100.0d0)
      END IF
      IPE(NROWS_LOC+1) = SAVEPNT
      CALL MUMPS_DEALLOC(BUFLEVEL, RDISPL, MEMCNT=MEMCNT)
      CALL MUMPS_I8DEALLOC(MSGCNT, SNDCNT, RCVCNT, MEMCNT=MEMCNT)
      DEALLOCATE(APNT)
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      DO I=1, LAST(MYID+1)-FIRST(MYID+1)+1
         L = int(IPE(I+1)-IPE(I))
         CALL DMUMPS_MERGESORT(L,
     &        PE(IPE(I):IPE(I+1)-1),
     &        WORK(:))
         CALL DMUMPS_MERGESWAP1(L, WORK(:),
     &        PE(IPE(I):IPE(I+1)-1))
      END DO
#endif
 90   continue
      RETURN
      END SUBROUTINE DMUMPS_BUILD_DIST_GRAPH
#endif      
      SUBROUTINE DMUMPS_BUILD_LOC_GRAPH(id, ord, GSIZE, IPE, PE, LENG,
     &     I_HALO_MAP, top_graph, WORK)
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC)   :: id
      TYPE(ORD_TYPE)       :: ord
      TYPE(GRAPH_TYPE)     :: top_graph
      INTEGER(8), POINTER  :: IPE(:)
      INTEGER, POINTER     :: PE(:), LENG(:),
     &     I_HALO_MAP(:), WORK(:)
      INTEGER              :: GSIZE
      INTEGER                :: IERR, MYID, NPROCS
      INTEGER                :: I, PROC, J, LOC_ROW
      INTEGER(8)             :: LOCNNZ, INNZ, NEW_LOCNNZ, TOP_CNT, TIDX,
     &                          RCVPNT
      INTEGER                :: IIDX,JJDX
      INTEGER                :: HALO_SIZE, NROWS_LOC, DUPS
      INTEGER                :: STATUS(MPI_STATUS_SIZE)
      INTEGER(8), POINTER    :: MSGCNT(:), SNDCNT(:), RCVCNT(:)
      INTEGER, POINTER       :: MAPTAB(:),
     &     SDISPL(:), HALO_MAP(:), BUFLEVEL(:)
      INTEGER, POINTER       :: RDISPL(:),
     &     SIPES(:,:)
      INTEGER, POINTER       :: PCNT(:), TSENDI(:),
     &     TSENDJ(:), RCVBUF(:)
      TYPE(ARRPNT), POINTER  :: APNT(:)
      INTEGER                :: BUFSIZE, SOURCE, MAXS, allocok
      INTEGER(8)             :: PNT, SAVEPNT
      INTEGER, PARAMETER     :: ITAG=30
      INTEGER(KIND=8)        :: TLEN
      LOGICAL                :: FLAG
      nullify(MAPTAB, SNDCNT, RCVCNT, SDISPL, HALO_MAP)
      nullify(RDISPL, MSGCNT, SIPES, BUFLEVEL)
      nullify(PCNT, TSENDI, TSENDJ, RCVBUF, APNT)
      CALL MPI_COMM_RANK (id%COMM, MYID, IERR)
      CALL MPI_COMM_SIZE (id%COMM, NPROCS, IERR)
      IF(MUMPS_GETSIZE(WORK) .LT. id%N*2) THEN
         WRITE(LP,
     &        '("Insufficient workspace inside BUILD_LOC_GRAPH")')
         CALL MUMPS_ABORT()
      END IF
      MAPTAB   => WORK(     1 :   id%N)
      HALO_MAP => WORK(id%N+1 : 2*id%N)
      CALL MUMPS_I8REALLOC(SNDCNT, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_I8REALLOC(RCVCNT, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_I8REALLOC(MSGCNT, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_REALLOC(RDISPL, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      ALLOCATE(APNT(NPROCS), stat=allocok)
      IF(allocok.GT.0) THEN
         id%INFO(1)=-13
         id%INFO(2)=NPROCS
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1), id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GO TO 90
      SNDCNT = 0
      TOP_CNT = 0
      BUFSIZE = 10000
      LOCNNZ = id%KEEP8(29)
      NROWS_LOC = ord%LAST(MYID+1)-ord%FIRST(MYID+1)+1
      MAPTAB = 0
      MAXS = 0
      DO I=1, NPROCS
         IF((ord%LAST(I)-ord%FIRST(I)+1) .GT. MAXS) THEN
            MAXS = ord%LAST(I)-ord%FIRST(I)+1
         END IF
         DO J=ord%FIRST(I), ord%LAST(I)
            MAPTAB(ord%PERITAB(J)) = I
         END DO
      END DO
      ALLOCATE(SIPES(max(1,MAXS), NPROCS), stat=allocok)
      IF(allocok.GT.0) THEN
         id%INFO(1)=-13
         id%INFO(2)=max(1,MAXS)*NPROCS
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1), id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GO TO 90
      SIPES(:,:)  = 0
      TOP_CNT     = 0
      DO INNZ=1, LOCNNZ
         IF(id%IRN_loc(INNZ) .NE. id%JCN_loc(INNZ)) THEN
            PROC = MAPTAB(id%IRN_loc(INNZ))
            IF(PROC .EQ. 0) THEN
               TOP_CNT = TOP_CNT+1
            ELSE
               IIDX = ord%PERMTAB(id%IRN_loc(INNZ))
               LOC_ROW = IIDX-ord%FIRST(PROC)+1
               SIPES(LOC_ROW, PROC) = SIPES(LOC_ROW, PROC)+1
               SNDCNT(PROC) = SNDCNT(PROC)+1
            END IF
            PROC = MAPTAB(id%JCN_loc(INNZ))
            IF(PROC .EQ. 0) THEN
               TOP_CNT = TOP_CNT+1
            ELSE
               IIDX = ord%PERMTAB(id%JCN_loc(INNZ))
               LOC_ROW = IIDX-ord%FIRST(PROC)+1
               SIPES(LOC_ROW, PROC) = SIPES(LOC_ROW, PROC)+1
               SNDCNT(PROC) = SNDCNT(PROC)+1
            END IF
         END IF
      END DO
      CALL MPI_ALLTOALL(SNDCNT(1), 1, MPI_INTEGER8, RCVCNT(1), 1,
     &     MPI_INTEGER8, id%COMM, IERR)
      I = ceiling(dble(MAXS)*1.20D0)
      CALL MUMPS_REALLOC(LENG, max(I,1), id%INFO,
     &        LP, STRING='B_L_G:LENG', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      CALL MUMPS_I8DEALLOC(SNDCNT, MEMCNT=MEMCNT)
      RDISPL(:) = MAXS
      CALL MPI_REDUCE_SCATTER ( SIPES(1,1), LENG(1), RDISPL(1), 
     &     MPI_INTEGER, MPI_SUM, id%COMM, IERR )
      DEALLOCATE(SIPES)
      I = ceiling(dble(NROWS_LOC+1)*1.20D0)
      CALL MUMPS_I8REALLOC(IPE, max(I,1), id%INFO,
     &        LP, STRING='B_L_G:IPE', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      TLEN = 0_8
      IPE(1) = 1_8
      DO I=1, NROWS_LOC
         IPE(I+1) = IPE(I) + int(LENG(I),8)
         TLEN = TLEN+int(LENG(I),8)
      END DO
      CALL MUMPS_IREALLOC8(TSENDI, max(TOP_CNT,1_8), id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_IREALLOC8(TSENDJ, max(TOP_CNT,1_8), id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      LENG(:) = 0
      CALL MUMPS_REALLOC(BUFLEVEL, NPROCS, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL DMUMPS_SEND_BUF(APNT, PROC, NPROCS, BUFSIZE, IPE, PE,
     &     LENG, RCVBUF, MSGCNT, BUFLEVEL, id%COMM)
      NEW_LOCNNZ = 0
      DO I=1, NPROCS
         NEW_LOCNNZ = NEW_LOCNNZ + RCVCNT(I)
         MSGCNT(I) = RCVCNT(I)/int(BUFSIZE,8)
      END DO
      CALL MUMPS_IREALLOC8(PE, max(NEW_LOCNNZ+
     &     2_8*int(NROWS_LOC+ord%TOPNODES(2),8),1_8),
     &     id%INFO, LP, STRING='B_L_G:PE', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      RCVPNT   = 1
      BUFLEVEL = 0 
      TIDX     = 0
      DO INNZ=1, LOCNNZ
         IF(mod(INNZ,int(BUFSIZE/10,8)) .EQ. 0) THEN
            CALL MPI_IPROBE( MPI_ANY_SOURCE, ITAG, id%COMM,
     &           FLAG, STATUS, IERR )
            IF(FLAG) THEN
               SOURCE = STATUS(MPI_SOURCE)
               CALL MPI_RECV(RCVBUF(1), 2*BUFSIZE, MPI_INTEGER, SOURCE,
     &              ITAG, id%COMM, STATUS, IERR)
               CALL DMUMPS_ASSEMBLE_MSG(BUFSIZE, RCVBUF, IPE, PE, LENG)
               MSGCNT(SOURCE+1)=MSGCNT(SOURCE+1)-1
               RCVPNT = RCVPNT + BUFSIZE
            END IF
         END IF
         IF(id%IRN_loc(INNZ) .NE. id%JCN_loc(INNZ)) THEN
            PROC = MAPTAB(id%IRN_loc(INNZ))
            IF((MAPTAB(id%JCN_loc(INNZ)).NE.PROC) .AND.
     &           (MAPTAB(id%JCN_loc(INNZ)).NE.0)  .AND.
     &           (PROC.NE.0)) THEN
               IERR = -50
               id%INFO(1) = IERR
            END IF
            IF(PROC .EQ. 0) THEN
               TIDX = TIDX+1
               TSENDI(TIDX) = id%IRN_loc(INNZ)
               TSENDJ(TIDX) = id%JCN_loc(INNZ)
            ELSE
               IIDX = ord%PERMTAB(id%IRN_loc(INNZ))
               JJDX = ord%PERMTAB(id%JCN_loc(INNZ))
               APNT(PROC)%BUF(2*BUFLEVEL(PROC)+1)=IIDX-ord%FIRST(PROC)+1
               IF( (JJDX .GE. ord%FIRST(PROC)) .AND.
     &              (JJDX .LE. ord%LAST(PROC)) ) THEN
                  APNT(PROC)%BUF(2*BUFLEVEL(PROC)+2) =
     &            JJDX-ord%FIRST(PROC)+1
               ELSE
                  APNT(PROC)%BUF(2*BUFLEVEL(PROC)+2) = -id%JCN_loc(INNZ)
               END IF
               BUFLEVEL(PROC) = BUFLEVEL(PROC)+1
               IF(BUFLEVEL(PROC) .EQ. BUFSIZE) THEN
                  CALL DMUMPS_SEND_BUF(APNT, PROC, NPROCS, BUFSIZE, IPE,
     &                 PE, LENG, RCVBUF, MSGCNT, BUFLEVEL, id%COMM)
               END IF
            END IF
            PROC = MAPTAB(id%JCN_loc(INNZ))
            IF(PROC .EQ. 0) THEN
               TIDX = TIDX+1
               TSENDI(TIDX) = id%JCN_loc(INNZ)
               TSENDJ(TIDX) = id%IRN_loc(INNZ)
            ELSE
               IIDX = ord%PERMTAB(id%JCN_loc(INNZ))
               JJDX = ord%PERMTAB(id%IRN_loc(INNZ))
               APNT(PROC)%BUF(2*BUFLEVEL(PROC)+1) =
     &         IIDX-ord%FIRST(PROC)+1
               IF( (JJDX .GE. ord%FIRST(PROC)) .AND.
     &              (JJDX .LE. ord%LAST(PROC)) ) THEN
                  APNT(PROC)%BUF(2*BUFLEVEL(PROC)+2) =
     &            JJDX-ord%FIRST(PROC)+1
               ELSE
                  APNT(PROC)%BUF(2*BUFLEVEL(PROC)+2) = -id%IRN_loc(INNZ)
               END IF
               BUFLEVEL(PROC) = BUFLEVEL(PROC)+1
               IF(BUFLEVEL(PROC) .EQ. BUFSIZE) THEN
                  CALL DMUMPS_SEND_BUF(APNT, PROC, NPROCS, BUFSIZE, IPE,
     &                 PE, LENG, RCVBUF, MSGCNT, BUFLEVEL, id%COMM)
               END IF
            END IF
         END IF
      END DO
      CALL DMUMPS_SEND_BUF(APNT, -1, NPROCS, BUFSIZE, IPE, PE, LENG,
     &     RCVBUF, MSGCNT, BUFLEVEL, id%COMM)
      DUPS = 0
      PNT = 0
      SAVEPNT = 1
      MAPTAB(:) = 0
      HALO_MAP(:) = 0
      HALO_SIZE = 0
      DO I=1, NROWS_LOC
         DO INNZ=IPE(I),IPE(I+1)-1
            IF(PE(INNZ) .LT. 0) THEN
               IF(HALO_MAP(-PE(INNZ)) .EQ. 0) THEN
                  HALO_SIZE = HALO_SIZE+1
                  HALO_MAP(-PE(INNZ)) = NROWS_LOC+HALO_SIZE
               END IF
               PE(INNZ) = HALO_MAP(-PE(INNZ))
            END IF
            IF(MAPTAB(PE(INNZ)) .EQ. I) THEN
               DUPS = DUPS+1
               LENG(I) = LENG(I)-1
            ELSE
               MAPTAB(PE(INNZ)) = I 
               PNT = PNT+1
               PE(PNT) = PE(INNZ)
            END IF
         END DO
         IPE(I) = SAVEPNT
         SAVEPNT = PNT+1
      END DO
      IPE(NROWS_LOC+1) = SAVEPNT
      CALL MUMPS_REALLOC(I_HALO_MAP, HALO_SIZE, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      J=0
      DO I=1, id%N
         IF(HALO_MAP(I) .GT. 0) THEN
            J = J+1
            I_HALO_MAP(HALO_MAP(I)-NROWS_LOC) = I
         END IF
         IF(J .EQ. HALO_SIZE) EXIT 
      END DO
      CALL MUMPS_REALLOC(LENG, max(NROWS_LOC+HALO_SIZE,1), id%INFO,
     &     LP, COPY=.TRUE.,
     &     STRING='lcgrph:leng', MEMCNT=MEMCNT, ERRCODE=-7)
      LENG(NROWS_LOC+1:NROWS_LOC+HALO_SIZE) = 0
      CALL MUMPS_I8REALLOC(IPE, NROWS_LOC+HALO_SIZE+1, id%INFO,
     &     LP, COPY=.TRUE.,
     &     STRING='lcgrph:ipe', MEMCNT=MEMCNT, ERRCODE=-7)
      IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
      IPE(NROWS_LOC+2:NROWS_LOC+HALO_SIZE+1) = IPE(NROWS_LOC+1)
      GSIZE = NROWS_LOC + HALO_SIZE
      CALL MPI_GATHER(TOP_CNT, 1, MPI_INTEGER8, RCVCNT(1), 1, 
     & MPI_INTEGER8, 0, id%COMM, IERR)
      IF(MYID.EQ.0) THEN
         NEW_LOCNNZ = sum(RCVCNT)
         top_graph%NZ_LOC = NEW_LOCNNZ
         top_graph%COMM = id%COMM
         CALL MUMPS_IREALLOC8(top_graph%IRN_LOC, max(1_8,NEW_LOCNNZ), 
     &        id%INFO, LP, MEMCNT=MEMCNT, ERRCODE=-7)
         CALL MUMPS_IREALLOC8(top_graph%JCN_LOC, max(1_8,NEW_LOCNNZ), 
     &        id%INFO, LP, MEMCNT=MEMCNT, ERRCODE=-7)
         IF(MEMCNT .GT. MAXMEM) MAXMEM=MEMCNT
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1), id%COMM, id%MYID)
         IF ( id%INFO(1) .LT. 0 ) GO TO 90
      ELSE
         ALLOCATE(top_graph%IRN_LOC(1), top_graph%JCN_LOC(1),
     &        stat=allocok)
         IF(allocok.GT.0) THEN
            id%INFO(1)=-13
            id%INFO(2)=2
         ENDIF
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1), id%COMM, id%MYID)
         IF ( id%INFO(1) .LT. 0 ) GO TO 90
      END IF
      IF(MYID.EQ.0) THEN
         top_graph%IRN_LOC(1:TOP_CNT) = TSENDI(1:TOP_CNT)
         top_graph%JCN_LOC(1:TOP_CNT) = TSENDJ(1:TOP_CNT)
         DO PROC=2, NPROCS
            DO WHILE (RCVCNT(PROC) .GT. 0)
               I = int(min(int(BUFSIZE,8), RCVCNT(PROC)))
               CALL MPI_RECV(top_graph%IRN_LOC(TOP_CNT+1), I,
     &              MPI_INTEGER, PROC-1, ITAG, id%COMM, STATUS, IERR)
               CALL MPI_RECV(top_graph%JCN_LOC(TOP_CNT+1), I,
     &              MPI_INTEGER, PROC-1, ITAG, id%COMM, STATUS, IERR)
               RCVCNT(PROC) = RCVCNT(PROC)-I
               TOP_CNT = TOP_CNT+I
            END DO
         END DO
      ELSE
         DO WHILE (TOP_CNT .GT. 0)
            I = int(MIN(int(BUFSIZE,8), TOP_CNT))
            CALL MPI_SEND(TSENDI(TOP_CNT-I+1), I,
     &           MPI_INTEGER, 0, ITAG, id%COMM, IERR)
            CALL MPI_SEND(TSENDJ(TOP_CNT-I+1), I,
     &           MPI_INTEGER, 0, ITAG, id%COMM, IERR)
            TOP_CNT = TOP_CNT-I
         END DO
      END IF
      CALL MUMPS_DEALLOC(BUFLEVEL, RDISPL, TSENDI,
     &     TSENDJ, MEMCNT=MEMCNT)
      CALL MUMPS_I8DEALLOC(MSGCNT, SNDCNT, RCVCNT, MEMCNT=MEMCNT)
      DEALLOCATE(APNT)
 90   continue
      RETURN
      END SUBROUTINE DMUMPS_BUILD_LOC_GRAPH
      SUBROUTINE DMUMPS_SEND_BUF(APNT, PROC, NPROCS, BUFSIZE, IPE, PE,
     &     LENG, RCVBUF, MSGCNT, SNDCNT, COMM)
      IMPLICIT NONE
      INTEGER                 :: NPROCS, PROC, COMM, allocok
      TYPE(ARRPNT)            :: APNT(:)
      INTEGER                 :: BUFSIZE
      INTEGER, POINTER        :: RCVBUF(:), LENG(:), PE(:)
      INTEGER                 :: SNDCNT(:)
      INTEGER(8)              :: MSGCNT(:), IPE(:)
      LOGICAL, SAVE           :: INIT = .TRUE.
      INTEGER, POINTER, SAVE  :: SPACE(:,:,:)
      LOGICAL, POINTER, SAVE  :: PENDING(:)
      INTEGER, POINTER, SAVE  :: REQ(:), CPNT(:)
      INTEGER                 :: IERR, MYID, I, SOURCE
      INTEGER(8)              :: TOTMSG
      LOGICAL                 :: FLAG, TFLAG
      INTEGER                 :: STATUS(MPI_STATUS_SIZE)
      INTEGER                 :: TSTATUS(MPI_STATUS_SIZE)
      INTEGER, PARAMETER      :: ITAG=30, FTAG=31
      INTEGER, POINTER        :: TMPI(:), RCVCNT(:)
      CALL MPI_COMM_RANK (COMM, MYID, IERR)
      CALL MPI_COMM_SIZE (COMM, NPROCS, IERR)
      IF(INIT) THEN
         ALLOCATE(SPACE(2*BUFSIZE, 2, NPROCS), stat=allocok)
         IF(allocok.GT.0) THEN
            write(*,*) "Allocation error of SPACE in DMUMPS_SEND_BUF"
            return
         ENDIF
         ALLOCATE(RCVBUF(2*BUFSIZE), stat=allocok)
         IF(allocok.GT.0) THEN
            write(*,*) "Allocation error of RCVBUF in DMUMPS_SEND_BUF"
            return
         ENDIF
         ALLOCATE(PENDING(NPROCS), CPNT(NPROCS), stat=allocok)
         IF(allocok.GT.0) THEN
            write(*,*) "Allocation error of PENDING/CPNT"
     &           ," in DMUMPS_SEND_BUF"
            return
         ENDIF
         ALLOCATE(REQ(NPROCS), stat=allocok)
         IF(allocok.GT.0) THEN
            write(*,*) "Allocation error of REQ in DMUMPS_SEND_BUF"
            return
         ENDIF
         PENDING = .FALSE.
         DO I=1, NPROCS
            APNT(I)%BUF => SPACE(:,1,I)
            CPNT(I)   = 1
         END DO
         INIT = .FALSE.
         RETURN
      END IF
      IF(PROC .EQ. -1) THEN
         TOTMSG = sum(MSGCNT)
         DO
            IF(TOTMSG .EQ. 0) EXIT
            CALL MPI_RECV(RCVBUF(1), 2*BUFSIZE, MPI_INTEGER,
     &           MPI_ANY_SOURCE, ITAG, COMM, STATUS, IERR)
            CALL DMUMPS_ASSEMBLE_MSG(BUFSIZE, RCVBUF, IPE, PE, LENG)
            SOURCE = STATUS(MPI_SOURCE)
            TOTMSG = TOTMSG-1
            MSGCNT(SOURCE+1)=MSGCNT(SOURCE+1)-1
         END DO
         DO I=1, NPROCS
            IF(PENDING(I)) THEN
               CALL MPI_WAIT(REQ(I), TSTATUS, IERR)
            END IF
         END DO
         ALLOCATE(RCVCNT(NPROCS), stat=allocok)
         IF(allocok.GT.0) THEN
            write(*,*) "Allocation error of RCVCNT in DMUMPS_SEND_BUF"
            return
         ENDIF
         CALL MPI_ALLTOALL(SNDCNT(1), 1, MPI_INTEGER, RCVCNT(1), 1,
     &        MPI_INTEGER, COMM, IERR)
         DO I=1, NPROCS
            IF(SNDCNT(I) .GT. 0) THEN
               TMPI => APNT(I)%BUF(:)
               CALL MPI_ISEND(TMPI(1), 2*SNDCNT(I), MPI_INTEGER, I-1,
     &              FTAG, COMM, REQ(I), IERR)
            END IF
         END DO
         DO I=1, NPROCS
            IF(RCVCNT(I) .GT. 0) THEN
               CALL MPI_RECV(RCVBUF(1), 2*RCVCNT(I), MPI_INTEGER, I-1,
     &              FTAG, COMM, STATUS, IERR)
               CALL DMUMPS_ASSEMBLE_MSG(RCVCNT(I), RCVBUF,
     &              IPE, PE, LENG)
            END IF
         END DO
         DO I=1, NPROCS
            IF(SNDCNT(I) .GT. 0) THEN
               CALL MPI_WAIT(REQ(I), TSTATUS, IERR)
            END IF
         END DO
         DEALLOCATE(SPACE)
         DEALLOCATE(PENDING, CPNT)
         DEALLOCATE(REQ)
         DEALLOCATE(RCVBUF, RCVCNT)
         nullify(SPACE, PENDING, CPNT, REQ, RCVBUF, RCVCNT)
         INIT = .TRUE.
         RETURN
      END IF
      IF(PENDING(PROC)) THEN
         DO
            CALL MPI_TEST(REQ(PROC), TFLAG, TSTATUS, IERR)
            IF(TFLAG) THEN
               PENDING(PROC) = .FALSE.
               EXIT
            ELSE
               CALL MPI_IPROBE( MPI_ANY_SOURCE, ITAG, COMM,
     &              FLAG, STATUS, IERR )
               IF(FLAG) THEN
                  SOURCE = STATUS(MPI_SOURCE)
                  CALL MPI_RECV(RCVBUF(1), 2*BUFSIZE, MPI_INTEGER,
     &                 SOURCE, ITAG, COMM, STATUS, IERR)
                  CALL DMUMPS_ASSEMBLE_MSG(BUFSIZE, RCVBUF, IPE,
     &                 PE, LENG)
                  MSGCNT(SOURCE+1)=MSGCNT(SOURCE+1)-1
               END IF
            END IF
         END DO
      END IF
      TMPI => APNT(PROC)%BUF(:)
      CALL MPI_ISEND(TMPI(1), 2*BUFSIZE, MPI_INTEGER, PROC-1,
     &     ITAG, COMM, REQ(PROC), IERR)
      PENDING(PROC) = .TRUE.
      CPNT(PROC) = mod(CPNT(PROC),2)+1
      APNT(PROC)%BUF => SPACE(:,CPNT(PROC),PROC)
      SNDCNT(PROC)  = 0
      RETURN
      END SUBROUTINE DMUMPS_SEND_BUF
      SUBROUTINE DMUMPS_ASSEMBLE_MSG(BUFSIZE, RCVBUF, IPE, PE, LENG)
      IMPLICIT NONE
      INTEGER             :: BUFSIZE
      INTEGER             :: RCVBUF(:), PE(:), LENG(:)
      INTEGER(8)          :: IPE(:)
      INTEGER             :: I, ROW, COL
      DO I=1, 2*BUFSIZE, 2
         ROW = RCVBUF(I)
         COL = RCVBUF(I+1)
         PE(IPE(ROW)+LENG(ROW)) = COL
         LENG(ROW) = LENG(ROW) + 1
      END DO
      RETURN
      END SUBROUTINE DMUMPS_ASSEMBLE_MSG 
#if defined(ptscotch) || defined(parmetis) || defined(parmetis3)
      SUBROUTINE DMUMPS_BUILD_TREE(ord)
      TYPE(ORD_TYPE)  :: ord
      INTEGER :: I
      ord%SON     = -1
      ord%BROTHER = -1
      ord%NW      = 0
      DO I=1, ord%CBLKNBR
         ord%NW(I) = ord%NW(I)+ord%RANGTAB(I+1) - ord%RANGTAB(I)  
         IF (ord%TREETAB(I) .NE. -1) THEN
            IF (ord%SON(ord%TREETAB(I)) .EQ. -1) THEN
               ord%SON(ord%TREETAB(I)) = I
            ELSE
               ord%BROTHER(I) = ord%SON(ord%TREETAB(I))
               ord%SON(ord%TREETAB(I)) = I
            END IF
            ord%NW(ord%TREETAB(I)) = ord%NW(ord%TREETAB(I))+ ord%NW(I)
         END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_BUILD_TREE
      SUBROUTINE DMUMPS_GRAPH_DIST(id, ord, FIRST,
     &     LAST, BASE, NPROCS, WORK, TYPE)
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC)   :: id
      TYPE(ORD_TYPE)       :: ord
      INTEGER              :: FIRST(:), LAST(:), BASE, NPROCS, TYPE
      INTEGER, TARGET      :: WORK(:)
      INTEGER, POINTER     :: TMP(:), NZ_ROW(:)
      INTEGER              :: I, IERR, P, F, J
      INTEGER(8)           :: LOCNNZ, INNZ, LOCOFFDIAG,
     &     OFFDIAG, T, SHARE
      DO I=0, BASE-1
         FIRST(I+1) = 0
         LAST(I+1)  = -1
      END DO
      IF(TYPE.EQ.1) THEN
         SHARE = int(id%N/ord%NSLAVES,8)
         DO I=1, ord%NSLAVES
            FIRST(BASE+I) = (I-1)*int(SHARE)+1
            LAST (BASE+I) = (I)*int(SHARE)
         END DO
         LAST(BASE+ord%NSLAVES) = MAX(LAST(BASE+ord%NSLAVES), id%N)
         DO I = ord%NSLAVES+1, id%NSLAVES+1
            FIRST(BASE+I) = id%N+1
            LAST (BASE+I) = id%N
         END DO
      ELSE IF (TYPE.EQ.2) THEN
         TMP    => WORK(1:id%N)
         NZ_ROW => WORK(id%N+1:2*id%N)
         TMP = 0
         LOCOFFDIAG = 0_8
         LOCNNZ = id%KEEP8(29)
         DO INNZ=1, LOCNNZ
            IF(id%IRN_loc(INNZ) .NE. id%JCN_loc(INNZ)) THEN
               TMP(id%IRN_loc(INNZ)) = TMP(id%IRN_loc(INNZ))+1
               LOCOFFDIAG = LOCOFFDIAG+1
               IF(id%SYM.GT.0) THEN
                  TMP(id%JCN_loc(INNZ)) = TMP(id%JCN_loc(INNZ))+1
                  LOCOFFDIAG = LOCOFFDIAG+1
               END IF
            END IF
         END DO
         CALL MPI_ALLREDUCE(TMP(1), NZ_ROW(1), id%N,
     &        MPI_INTEGER, MPI_SUM, id%COMM, IERR)
         CALL MPI_ALLREDUCE(LOCOFFDIAG, OFFDIAG, 1,
     &        MPI_INTEGER8, MPI_SUM, id%COMM, IERR)
         nullify(TMP)
         SHARE = (OFFDIAG-1_8)/int(ord%NSLAVES,8) + 1_8
         P = 0
         T = 0_8
         F = 1
         DO I=1, id%N
            T = T+int(NZ_ROW(I),8)
            IF (
     &           (T .GE. SHARE) .OR.
     &           ((id%N-I).EQ.(ord%NSLAVES-P-1)) .OR.
     &           (I.EQ.id%N)
     &           ) THEN
               P             = P+1
               IF(P.EQ.ord%NSLAVES) THEN
                  FIRST(BASE+P) = F
                  LAST(BASE+P)  = id%N
                  EXIT
               ELSE
                  FIRST(BASE+P) = F
                  LAST(BASE+P)  = I
                  F             = I+1
                  T             = 0_8
               END IF
            END IF
         END DO
         DO J=P+1, NPROCS+1-BASE
            FIRST(BASE+J) = id%N+1
            LAST(BASE+J)  = id%N
         END DO
      END IF         
      RETURN
      END SUBROUTINE DMUMPS_GRAPH_DIST
#endif     
      SUBROUTINE DMUMPS_MERGESWAP(N, L, A1, A2)
      INTEGER   :: I, LP, ISWAP, N
      INTEGER   :: L(0:), A1(:), A2(:)
      LP = L(0)
      I  = 1
      DO 
         IF ((LP==0).OR.(I>N)) EXIT
         DO 
            IF (LP >= I) EXIT
            LP = L(LP)
         END DO
         ISWAP    = A1(LP)
         A1(LP)   = A1(I)
         A1(I)    = ISWAP
         ISWAP    = A2(LP)
         A2(LP)   = A2(I)
         A2(I)    = ISWAP
         ISWAP    = L(LP)
         L(LP) = L(I)
         L(I)  = LP
         LP = ISWAP 
         I  = I + 1
      ENDDO
      END SUBROUTINE DMUMPS_MERGESWAP
#if defined(DETERMINISTIC_PARALLEL_GRAPH)  
      SUBROUTINE DMUMPS_MERGESWAP1(N, L, A)
      INTEGER   :: I, LP, ISWAP, N
      INTEGER   :: L(0:), A(:)
      LP = L(0)
      I  = 1
      DO 
         IF ((LP==0).OR.(I>N)) EXIT
         DO 
            IF (LP >= I) EXIT
            LP = L(LP)
         END DO
         ISWAP    = A(LP)
         A(LP)   = A(I)
         A(I)    = ISWAP
         ISWAP    = L(LP)
         L(LP) = L(I)
         L(I)  = LP
         LP = ISWAP 
         I  = I + 1
      ENDDO
      END SUBROUTINE DMUMPS_MERGESWAP1
#endif
      SUBROUTINE DMUMPS_MERGESORT(N, K, L)
      INTEGER    :: N
      INTEGER    :: K(:), L(0:)
      INTEGER    :: P, Q, S, T
      CONTINUE
      L(0) = 1
      T = N + 1
      DO  P = 1,N - 1
         IF (K(P) <= K(P+1)) THEN
            L(P) = P + 1
         ELSE
            L(T) = - (P+1)
            T = P
       END IF
      END DO
      L(T) = 0
      L(N) = 0
      IF (L(N+1) == 0) THEN
         RETURN 
      ELSE
         L(N+1) = iabs(L(N+1))
      END IF
 200  CONTINUE
      S = 0
      T = N+1
      P = L(S)
      Q = L(T)
      IF(Q .EQ. 0) RETURN
 300  CONTINUE
      IF(K(P) .GT. K(Q)) GOTO 600 
      CONTINUE
      L(S) = sign(P,L(S))
      S = P
      P = L(P)
      IF (P .GT. 0) GOTO 300
      CONTINUE
      L(S) = Q
      S = T
      DO
         T = Q
         Q = L(Q)
         IF (Q .LE. 0) EXIT
      END DO
      GOTO 800
 600  CONTINUE
      L(S) = sign(Q, L(S))
      S = Q
      Q = L(Q)
      IF (Q .GT. 0) GOTO 300
      CONTINUE
      L(S) = P
      S = T
      DO
         T = P
         P = L(P)
         IF (P .LE. 0) EXIT
      END DO
 800  CONTINUE
      P = -P
      Q = -Q
      IF(Q.EQ.0) THEN
         L(S) = sign(P, L(S))
         L(T) = 0
         GOTO 200
      END IF
      GOTO 300
      END SUBROUTINE DMUMPS_MERGESORT
      FUNCTION MUMPS_GETSIZE(A)
      INTEGER, POINTER :: A(:)
      INTEGER          :: MUMPS_GETSIZE
      IF(associated(A)) THEN
         MUMPS_GETSIZE = size(A)
      ELSE
         MUMPS_GETSIZE = 0_8
      END IF
      RETURN
      END FUNCTION MUMPS_GETSIZE
#if defined(parmetis) || defined(parmetis3)
      SUBROUTINE MUMPS_PARMETIS_MIXEDto32(id, BASE, VERTLOCNBR, FIRST,
     &     VERTLOCTAB, EDGELOCTAB, BASEVAL, OPTIONS, ORDER,
     &     SIZES, COMM, IERR)
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC) :: id
      INTEGER            :: FIRST(:), EDGELOCTAB(:), OPTIONS(:)
      INTEGER            :: SIZES(:), ORDER(:)
      INTEGER(8)         :: VERTLOCTAB(:)
      INTEGER            :: VERTLOCNBR, BASEVAL, IERR, COMM, BASE
      INTEGER, POINTER   :: VERTLOCTAB_I4(:)
      IF( VERTLOCTAB(VERTLOCNBR+1).GT.huge(VERTLOCNBR)) THEN
         id%INFO(1) = -51
         CALL MUMPS_SET_IERROR(
     &    VERTLOCTAB(VERTLOCNBR+1), id%INFO(2))
         RETURN
      END IF
      nullify(VERTLOCTAB_I4)
      CALL MUMPS_REALLOC(VERTLOCTAB_I4, VERTLOCNBR+1, id%INFO,
     &     id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MUMPS_ICOPY_64TO32(VERTLOCTAB(1),
     &     VERTLOCNBR+1, VERTLOCTAB_I4(1))
      CALL MUMPS_PARMETIS(FIRST(1+BASE), VERTLOCTAB_I4(1),
     &     EDGELOCTAB(1), BASEVAL, OPTIONS(1), ORDER(1),
     &     SIZES(1), COMM, IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_DEALLOC(VERTLOCTAB_I4, MEMCNT=MEMCNT)
      RETURN
      END SUBROUTINE MUMPS_PARMETIS_MIXEDto32
      SUBROUTINE MUMPS_PARMETIS_MIXEDto64
     &     (id, BASE, VERTLOCNBR, FIRST,
     &     VERTLOCTAB, EDGELOCTAB, BASEVAL, OPTIONS, ORDER,
     &     SIZES, COMM, IERR)
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC) :: id
      INTEGER            :: FIRST(:), EDGELOCTAB(:), OPTIONS(:)
      INTEGER            :: SIZES(:), ORDER(:)
      INTEGER(8)         :: VERTLOCTAB(:)
      INTEGER            :: VERTLOCNBR, BASEVAL, IERR, COMM, BASE
      INTEGER(8), POINTER :: FIRST_I8(:), EDGELOCTAB_I8(:),
     &     SIZES_I8(:), ORDER_I8(:)
#if defined(parmetis)
      INTEGER(8), POINTER :: OPTIONS_I8(:)
      INTEGER(8)          :: BASEVAL_I8
      nullify(OPTIONS_I8)
      IF (id%KEEP(10).NE.1) THEN
       CALL MUMPS_I8REALLOC(OPTIONS_I8, size(OPTIONS), id%INFO,
     &     id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
       IF ( id%INFO(1) .LT. 0 ) RETURN
       CALL MUMPS_ICOPY_32TO64(OPTIONS(1), size(OPTIONS)
     &      , OPTIONS_I8(1))
       BASEVAL_I8 = int(BASEVAL,8)
      END IF
#endif
      nullify(FIRST_I8, EDGELOCTAB_I8, SIZES_I8, ORDER_I8)
      IF (id%KEEP(10).EQ.1) THEN
       CALL MUMPS_PARMETIS_64(FIRST(1+BASE), VERTLOCTAB(1),
     &     EDGELOCTAB(1),
     &     BASEVAL, OPTIONS(1),
     &     ORDER(1),
     &     SIZES(1), COMM, IERR)
      ELSE
       CALL MUMPS_I8REALLOC(FIRST_I8, size(FIRST), id%INFO,
     &     id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
       IF ( id%INFO(1) .LT. 0 ) GOTO 5
       CALL MUMPS_I8REALLOC8(EDGELOCTAB_I8, 
     &     VERTLOCTAB(VERTLOCNBR+1)-1_8,
     &     id%INFO, id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
       IF ( id%INFO(1) .LT. 0 ) GOTO 5
       CALL MUMPS_I8REALLOC(SIZES_I8, size(SIZES), id%INFO,
     &     id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
       IF ( id%INFO(1) .LT. 0 ) GOTO 5
       CALL MUMPS_I8REALLOC(ORDER_I8, size(ORDER), id%INFO,
     &     id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
       IF ( id%INFO(1) .LT. 0 ) GOTO 5
 5     CONTINUE
       CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     COMM, id%MYID )
       IF ( id%INFO(1) .LT. 0 ) RETURN
       CALL MUMPS_ICOPY_32TO64(FIRST(1), size(FIRST), FIRST_I8(1))
       CALL MUMPS_ICOPY_32TO64_64C(EDGELOCTAB(1),
     &     VERTLOCTAB(VERTLOCNBR+1)-1_8, EDGELOCTAB_I8(1))
       CALL MUMPS_PARMETIS_64(FIRST_I8(1+BASE), VERTLOCTAB(1),
     &     EDGELOCTAB_I8(1),
#if defined(parmetis3)
     &     BASEVAL, OPTIONS(1),
#else
     &     BASEVAL_I8, OPTIONS_I8(1),
#endif
     &     ORDER_I8(1),
     &     SIZES_I8(1), COMM, IERR)
      END IF
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      CALL MUMPS_ICOPY_64TO32(ORDER_I8(1),
     &     size(ORDER), ORDER(1))
      CALL MUMPS_ICOPY_64TO32(SIZES_I8(1),
     &     size(SIZES), SIZES(1))
 10   CONTINUE 
      CALL MUMPS_I8DEALLOC(FIRST_I8, MEMCNT=MEMCNT)
      CALL MUMPS_I8DEALLOC(SIZES_I8, MEMCNT=MEMCNT)
      CALL MUMPS_I8DEALLOC(EDGELOCTAB_I8, MEMCNT=MEMCNT)
      CALL MUMPS_I8DEALLOC(ORDER_I8, MEMCNT=MEMCNT)
#if defined(parmetis)
      CALL MUMPS_I8DEALLOC(OPTIONS_I8, MEMCNT=MEMCNT)
#endif
      RETURN
      END SUBROUTINE MUMPS_PARMETIS_MIXEDto64
#endif
#if defined(ptscotch)
      SUBROUTINE MUMPS_PTSCOTCH_MIXEDto32(id, ord,
     &     BASEVAL, 
     &     VERTLOCNBR, VERTLOCTAB,
     &     EDGELOCNBR, EDGELOCTAB,
     &     IERR)
      IMPLICIT NONE
      INCLUDE 'ptscotchf.h'
      TYPE(DMUMPS_STRUC)  :: id
      TYPE(ORD_TYPE)      :: ord
      INTEGER             :: BASEVAL, VERTLOCNBR
      INTEGER(8)          :: EDGELOCNBR
      INTEGER(8)          :: VERTLOCTAB(:)
      INTEGER             :: EDGELOCTAB(:)
      INTEGER             :: IERR
      INTEGER, POINTER    :: VERTLOCTAB_I4(:)
      INTEGER             :: EDGELOCNBR_I4, MYWORKID
      DOUBLE PRECISION    :: GRAPHDAT(SCOTCH_DGRAPHDIM),
     &     ORDEDAT(SCOTCH_DORDERDIM), STRADAT(SCOTCH_STRATDIM),
     &     CORDEDAT(SCOTCH_ORDERDIM)
      CHARACTER  STRSTRING*1024
      nullify(VERTLOCTAB_I4)
      CALL MUMPS_REALLOC(VERTLOCTAB_I4, VERTLOCNBR+1, id%INFO, LP,
     &     MEMCNT=MEMCNT, ERRCODE=-7)
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN      
      CALL MUMPS_ICOPY_64TO32(VERTLOCTAB(1),
     &     VERTLOCNBR+1, VERTLOCTAB_I4(1))
      EDGELOCNBR_I4 = int(EDGELOCNBR)
      IF(ord%SUBSTRAT .NE. 0) THEN
         STRSTRING='n{sep=m{asc=b{width=3,strat=q{strat=f}},'//
     &        'low=q{strat=h},vert=1000,dvert=100,dlevl=0,'//
     &        'proc=1,seq=q{strat=m{type=h,vert=100,'//
     &        'low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},'//
     &        'org=h{pass=10}f{bal=0.2}}}}},ole=s,ose=s,osq=s}'
      END IF
      IF(ord%IDO) THEN
         CALL MPI_COMM_RANK (ord%COMM_NODES, MYWORKID, IERR)
      ELSE
         MYWORKID = -1
      END IF
      CALL MUMPS_DGRAPHINIT(GRAPHDAT, ord%COMM_NODES, IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      CALL SCOTCHFDGRAPHBUILD(GRAPHDAT, BASEVAL, VERTLOCNBR,
     &     VERTLOCNBR, VERTLOCTAB_I4(1), VERTLOCTAB_I4(2),
     &     VERTLOCTAB_I4(1), VERTLOCTAB_I4(1), EDGELOCNBR_I4,
     &     EDGELOCNBR_I4, EDGELOCTAB(1), EDGELOCTAB(1),
     &     EDGELOCTAB(1), IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      CALL SCOTCHFSTRATINIT(STRADAT, IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      IF(ord%SUBSTRAT .NE. 0) THEN
         CALL SCOTCHFSTRATDGRAPHORDER(STRADAT, STRSTRING, IERR)
      END IF
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      CALL SCOTCHFDGRAPHORDERINIT(GRAPHDAT, ORDEDAT, IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      CALL SCOTCHFDGRAPHORDERCOMPUTE(GRAPHDAT, ORDEDAT, STRADAT,
     &     IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      IF(MYWORKID .EQ. 0) THEN
         CALL SCOTCHFDGRAPHCORDERINIT(GRAPHDAT, CORDEDAT,
     &        ord%PERMTAB(1), ord%PERITAB(1), ord%CBLKNBR,
     &        ord%RANGTAB(1), ord%TREETAB(1), IERR)
         IF(IERR.NE.0) THEN
            id%INFO(1:2) = -50
         END IF
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      IF(MYWORKID .EQ. 0) THEN
         CALL SCOTCHFDGRAPHORDERGATHER(GRAPHDAT, ORDEDAT,
     &        CORDEDAT, IERR)
         IF(IERR.NE.0) THEN
            id%INFO(1:2) = -50
         END IF
      ELSE
         CALL SCOTCHFDGRAPHORDERGATHER(GRAPHDAT, ORDEDAT,
     &        ORDEDAT, IERR)
         IF(IERR.NE.0) THEN
            id%INFO(1:2) = -50
         END IF
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      IF(MYWORKID .EQ. 0) 
     &     CALL SCOTCHFDGRAPHCORDEREXIT(GRAPHDAT, CORDEDAT)
      CALL SCOTCHFDGRAPHORDEREXIT(GRAPHDAT, ORDEDAT)
      CALL SCOTCHFSTRATEXIT(STRADAT) 
      CALL SCOTCHFDGRAPHEXIT(GRAPHDAT)
 10   CONTINUE
      CALL MUMPS_DEALLOC(VERTLOCTAB_I4, MEMCNT=MEMCNT)
      RETURN
      END SUBROUTINE MUMPS_PTSCOTCH_MIXEDto32
      SUBROUTINE MUMPS_PTSCOTCH_MIXEDto64(id, ord,
     &     BASEVAL, 
     &     VERTLOCNBR, VERTLOCTAB,
     &     EDGELOCNBR, EDGELOCTAB,
     &     IERR)
      IMPLICIT NONE
      INCLUDE 'ptscotchf.h'
      TYPE(DMUMPS_STRUC)  :: id
      TYPE(ORD_TYPE)      :: ord
      INTEGER             :: BASEVAL, VERTLOCNBR
      INTEGER(8)          :: EDGELOCNBR
      INTEGER(8)          :: VERTLOCTAB(:)
      INTEGER             :: EDGELOCTAB(:)
      INTEGER             :: IERR
      INTEGER             :: MYWORKID
      DOUBLE PRECISION    :: GRAPHDAT(SCOTCH_DGRAPHDIM),
     &     ORDEDAT(SCOTCH_DORDERDIM), STRADAT(SCOTCH_STRATDIM),
     &     CORDEDAT(SCOTCH_ORDERDIM)
      CHARACTER  STRSTRING*1024
      INTEGER(8), POINTER :: EDGELOCTAB_I8(:), PERMTAB_I8(:),
     &     PERITAB_I8(:), RANGTAB_I8(:), TREETAB_I8(:)
      INTEGER(8)          :: CBLKNBR_I8, VERTLOCNBR_I8, BASEVAL_I8
      IF(ord%SUBSTRAT .NE. 0) THEN
         STRSTRING='n{sep=m{asc=b{width=3,strat=q{strat=f}},'//
     &        'low=q{strat=h},vert=1000,dvert=100,dlevl=0,'//
     &        'proc=1,seq=q{strat=m{type=h,vert=100,'//
     &        'low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},'//
     &        'org=h{pass=10}f{bal=0.2}}}}},ole=s,ose=s,osq=s}'
      END IF
      IF(ord%IDO) THEN
         CALL MPI_COMM_RANK (ord%COMM_NODES, MYWORKID, IERR)
      ELSE
         MYWORKID = -1
      END IF
      nullify(EDGELOCTAB_I8, PERMTAB_I8, PERITAB_I8,
     &     RANGTAB_I8, TREETAB_I8)
      IF (id%KEEP(10).NE.1) THEN 
       CALL MUMPS_I8REALLOC8(EDGELOCTAB_I8, 
     &     VERTLOCTAB(VERTLOCNBR+1)-1_8,
     &     id%INFO, id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
       IF ( id%INFO(1) .LT. 0 ) GOTO 5
       IF (MYWORKID .EQ. 0) THEN
         CALL MUMPS_I8REALLOC(PERMTAB_I8, size(ord%PERMTAB),
     &        id%INFO, id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
         IF ( id%INFO(1) .LT. 0 ) GOTO 5
         CALL MUMPS_I8REALLOC(PERITAB_I8, size(ord%PERITAB),
     &        id%INFO, id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
         IF ( id%INFO(1) .LT. 0 ) GOTO 5
         CALL MUMPS_I8REALLOC(TREETAB_I8, size(ord%TREETAB),
     &        id%INFO, id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
         IF ( id%INFO(1) .LT. 0 ) GOTO 5
         CALL MUMPS_I8REALLOC(RANGTAB_I8, size(ord%RANGTAB),
     &        id%INFO, id%ICNTL(1), MEMCNT=MEMCNT, ERRCODE=-7)
       END IF
 5     CONTINUE
       CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
       IF ( id%INFO(1) .LT. 0 ) RETURN
       CALL MUMPS_ICOPY_32TO64_64C(EDGELOCTAB(1),
     &     VERTLOCTAB(VERTLOCNBR+1)-1_8, EDGELOCTAB_I8(1))
       BASEVAL_I8    = int(BASEVAL,8)
       VERTLOCNBR_I8 = int(VERTLOCNBR,8)
      ENDIF
      CALL MUMPS_DGRAPHINIT(GRAPHDAT, ord%COMM_NODES, IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      IF (id%KEEP(10).NE.1) THEN
       CALL SCOTCHFDGRAPHBUILD(GRAPHDAT, BASEVAL_I8, VERTLOCNBR_I8,
     &     VERTLOCNBR_I8, VERTLOCTAB(1), VERTLOCTAB(2),
     &     VERTLOCTAB(1), VERTLOCTAB(1), EDGELOCNBR,
     &     EDGELOCNBR, EDGELOCTAB_I8(1), EDGELOCTAB_I8(1),
     &     EDGELOCTAB_I8(1), IERR)
      ELSE
       CALL SCOTCHFDGRAPHBUILD(GRAPHDAT, BASEVAL, VERTLOCNBR,
     &     VERTLOCNBR, VERTLOCTAB(1), VERTLOCTAB(2),
     &     VERTLOCTAB(1), VERTLOCTAB(1), EDGELOCNBR,
     &     EDGELOCNBR, EDGELOCTAB(1), EDGELOCTAB(1),
     &     EDGELOCTAB(1), IERR)
      ENDIF
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      CALL SCOTCHFSTRATINIT(STRADAT, IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      IF(ord%SUBSTRAT .NE. 0) THEN
         CALL SCOTCHFSTRATDGRAPHORDER(STRADAT, STRSTRING, IERR)
      END IF
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      CALL SCOTCHFDGRAPHORDERINIT(GRAPHDAT, ORDEDAT, IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      CALL SCOTCHFDGRAPHORDERCOMPUTE(GRAPHDAT, ORDEDAT, STRADAT,
     &     IERR)
      IF(IERR.NE.0) THEN
         id%INFO(1:2) = -50
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      IF(MYWORKID .EQ. 0) THEN
        IF (id%KEEP(10).NE.1) THEN 
         CALL SCOTCHFDGRAPHCORDERINIT(GRAPHDAT, CORDEDAT,
     &        PERMTAB_I8(1), PERITAB_I8(1), CBLKNBR_I8, RANGTAB_I8(1),
     &        TREETAB_I8(1), IERR)
        ELSE
         CALL SCOTCHFDGRAPHCORDERINIT(GRAPHDAT, CORDEDAT,
     &        ord%PERMTAB(1), ord%PERITAB(1), ord%CBLKNBR,
     &         ord%RANGTAB(1),ord%TREETAB(1), IERR)
        ENDIF
         IF(IERR.NE.0) THEN
            id%INFO(1:2) = -50
         END IF
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      IF(MYWORKID .EQ. 0) THEN
         CALL SCOTCHFDGRAPHORDERGATHER(GRAPHDAT, ORDEDAT,
     &        CORDEDAT, IERR)
         IF(IERR.NE.0) THEN
            id%INFO(1:2) = -50
         END IF
      ELSE
         CALL SCOTCHFDGRAPHORDERGATHER(GRAPHDAT, ORDEDAT,
     &        ORDEDAT, IERR)
         IF(IERR.NE.0) THEN
            id%INFO(1:2) = -50
         END IF
      END IF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     ord%COMM_NODES, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 10
      CALL SCOTCHFDGRAPHORDEREXIT(GRAPHDAT, ORDEDAT)
      CALL SCOTCHFSTRATEXIT(STRADAT) 
      CALL SCOTCHFDGRAPHEXIT(GRAPHDAT)
 10   CONTINUE
      IF (id%KEEP(10).NE.1) THEN
       CALL MUMPS_I8DEALLOC(EDGELOCTAB_I8, MEMCNT=MEMCNT)
       IF(MYWORKID .EQ. 0) THEN
         CALL SCOTCHFDGRAPHCORDEREXIT(GRAPHDAT, CORDEDAT)
         CALL MUMPS_ICOPY_64TO32(PERMTAB_I8(1),
     &        size(ord%PERMTAB), ord%PERMTAB(1))
         CALL MUMPS_ICOPY_64TO32(PERITAB_I8(1),
     &        size(ord%PERITAB), ord%PERITAB(1))
         CALL MUMPS_ICOPY_64TO32(TREETAB_I8(1),
     &        size(ord%TREETAB), ord%TREETAB(1))
         CALL MUMPS_ICOPY_64TO32(RANGTAB_I8(1),
     &        size(ord%RANGTAB), ord%RANGTAB(1))
         ord%CBLKNBR = int(CBLKNBR_I8)
         CALL MUMPS_I8DEALLOC(PERMTAB_I8, MEMCNT=MEMCNT)
         CALL MUMPS_I8DEALLOC(PERITAB_I8, MEMCNT=MEMCNT)
         CALL MUMPS_I8DEALLOC(RANGTAB_I8, MEMCNT=MEMCNT)
         CALL MUMPS_I8DEALLOC(TREETAB_I8, MEMCNT=MEMCNT)
       END IF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_PTSCOTCH_MIXEDto64
#endif      
      END MODULE
