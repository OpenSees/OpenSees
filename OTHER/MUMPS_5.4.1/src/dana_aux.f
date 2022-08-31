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
      MODULE DMUMPS_ANA_AUX_M
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE DMUMPS_ANA_F(N, NZ8, IRN, ICN, LIWALLOC,
     &     IKEEP1, IKEEP2, IKEEP3,
     &     IORD, NFSIZ, FILS, FRERE, LISTVAR_SCHUR, SIZE_SCHUR,
     &     ICNTL, INFO, KEEP,KEEP8, NSLAVES, PIV, 
     &     CNTL4, COLSCA, ROWSCA
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)         
     &          , METIS_OPTIONS
#endif      
     &          , NORIG_ARG, SIZEOFBLOCKS, GCOMP_PROVIDED_IN, GCOMP
     & )
      USE MUMPS_ANA_ORD_WRAPPERS
      USE MUMPS_ANA_BLK_M, ONLY : COMPACT_GRAPH_T
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N,  SIZE_SCHUR, NSLAVES
      INTEGER(8), INTENT(IN) :: NZ8
      INTEGER(8), INTENT(IN) :: LIWALLOC
      INTEGER, INTENT(in) :: LISTVAR_SCHUR(:)
      INTEGER, POINTER :: IRN(:), ICN(:)
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, INTENT(INOUT) :: IORD
      INTEGER, INTENT(INOUT) :: INFO(80), KEEP(500)
      INTEGER(8), INTENT(INOUT) :: KEEP8(150)
      INTEGER, INTENT(OUT)   :: NFSIZ(:), FILS(:), FRERE(:)
      INTEGER, INTENT(INOUT) :: PIV(:)
      INTEGER, INTENT(INOUT) :: IKEEP1(:), IKEEP2(:), IKEEP3(:)
      DOUBLE PRECISION :: CNTL4
      DOUBLE PRECISION, POINTER :: COLSCA(:), ROWSCA(:)
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)      
      INTEGER, INTENT(IN)    :: METIS_OPTIONS(40)
#endif      
      INTEGER, INTENT(IN), OPTIONAL     :: NORIG_ARG
      INTEGER, INTENT(IN), OPTIONAL     :: SIZEOFBLOCKS(N)
      LOGICAL, INTENT(IN), OPTIONAL     :: GCOMP_PROVIDED_IN
      TYPE(COMPACT_GRAPH_T), OPTIONAL   :: GCOMP
      INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: IWALLOC
      INTEGER, DIMENSION(:), POINTER             :: IW
      INTEGER(8), DIMENSION(:), ALLOCATABLE, TARGET :: IPEALLOC
      INTEGER(8), DIMENSION(:), POINTER             :: IPE 
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IPQ8
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PTRAR
      INTEGER, DIMENSION(:), ALLOCATABLE   :: PARENT
      INTEGER, DIMENSION(:), ALLOCATABLE   :: IWL1
      INTEGER NBBUCK
      INTEGER, DIMENSION(:), ALLOCATABLE   :: WTEMP
      INTEGER IERR
      INTEGER I, K, NCMPA, IN, IFSON
      INTEGER(8) :: J8, I8
      INTEGER    :: NORIG
      INTEGER(8) :: IFIRST, ILAST
      INTEGER(8) IWFR8
      INTEGER NEMIN, LP, MP, LDIAG, ITEMP, symmetry
      INTEGER NBQD, AvgDens
      LOGICAL PROK, COMPRESS_SCHUR, LPOK, COMPUTE_PERM
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
#if defined(metis4) || defined(parmetis3)
      INTEGER NUMFLAG
#endif      
      INTEGER METIS_IDX_SIZE
      INTEGER OPT_METIS_SIZE
#endif
#if defined(scotch) || defined(ptscotch)
      INTEGER :: SCOTCH_INT_SIZE
#endif
#if defined(pord)
      INTEGER :: PORD_INT_SIZE
#endif
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COLSCA_TEMP
      INTEGER THRESH, IVersion
      LOGICAL AGG6
      INTEGER MINSYM
      PARAMETER (MINSYM=50)
      INTEGER(8) :: K79REF
      PARAMETER(K79REF=12000000_8)
      INTEGER, PARAMETER :: LIDUMMY = 1 
      INTEGER            :: IDUMMY(1)   
      INTEGER MTRANS, COMPRESS,NCMP,IERROR,J,JPERM,NCST
      INTEGER TOTEL 
#if defined(pord)
      INTEGER TOTW
#endif
      INTEGER WEIGHTUSED, WEIGHTREQUESTED
      LOGICAL IDENT,SPLITROOT
      LOGICAL FREE_CENTRALIZED_MATRIX
      LOGICAL GCOMP_PROVIDED
      LOGICAL INPLACE64_GRAPH_COPY, INPLACE64_RESTORE_GRAPH
      INTEGER(8) :: LIW8, NZG8
      DOUBLE PRECISION TIMEB
      EXTERNAL MUMPS_ANA_H, DMUMPS_ANA_J, 
     &     DMUMPS_ANA_K, DMUMPS_ANA_GNEW,
     &     DMUMPS_ANA_LNEW, DMUMPS_ANA_M
#if defined(OLDDFS)
      EXTERNAL DMUMPS_ANA_L
#endif
      EXTERNAL DMUMPS_GNEW_SCHUR
      EXTERNAL DMUMPS_LDLT_COMPRESS, DMUMPS_EXPAND_PERMUTATION,
     &     DMUMPS_SET_CONSTRAINTS
      IF (LIWALLOC.GT.0_8) THEN
        ALLOCATE( IWALLOC (LIWALLOC), stat = IERR )
        IF ( IERR .GT. 0 ) THEN
         INFO( 1 ) = -7
         CALL MUMPS_SET_IERROR(LIWALLOC,INFO(2))
         GOTO 90                
        ENDIF
      ENDIF
      ALLOCATE( IWL1 (N), stat = IERR )
      IF ( IERR .GT. 0 ) THEN
         INFO( 1 ) = -7
         INFO( 2 ) = N
         GOTO 90                
      ENDIF
      ALLOCATE( IPEALLOC(N+1), stat = IERR )
      IF ( IERR .GT. 0 ) THEN
         INFO( 1 ) = -7
         INFO( 2 ) = (N+1)*KEEP(10)
         GOTO 90                
      ENDIF
      ALLOCATE( PTRAR (N,3), stat = IERR )
      IF ( IERR .GT. 0 ) THEN
         INFO( 1 ) = -7
         INFO( 2 ) = 3*N
         GOTO 90                
      ENDIF
      symmetry      = INFO(8)
      NBQD          = 0
      GCOMP_PROVIDED=.FALSE.
      WEIGHTUSED    = 0
      NORIG         = N
      IF (present(NORIG_ARG)) THEN
        NORIG=NORIG_ARG
      ENDIF
      IF (present(GCOMP_PROVIDED_IN))
     &       GCOMP_PROVIDED = GCOMP_PROVIDED_IN
      IF (GCOMP_PROVIDED.AND.(.NOT. present(GCOMP))) THEN
        INFO(1) = -900
        WRITE(6,*) " INTERNAL ERROR in MUMPS(ANA_F) ", 
     &  GCOMP_PROVIDED_IN, present(GCOMP)
        INFO(2) = 1
        RETURN
      ENDIF
      IF ( (LIWALLOC.EQ.0_8).AND.(.not.GCOMP_PROVIDED)) THEN
        INFO(1) = -900
        WRITE(6,*) " INTERNAL ERROR in MUMPS(ANA_F) ", 
     &  "LIWALLOC, GCOMP_PROVIDED=", LIWALLOC, GCOMP_PROVIDED
        INFO(2) = 2
        RETURN
      ENDIF
      IF (GCOMP_PROVIDED) THEN
        NZG8 =  GCOMP%NZG
        LIW8 =  NZG8 + int(GCOMP%NG,8)+1_8
        IW   => GCOMP%ADJ(1:LIW8)
        IPE  => GCOMP%IPE(1:GCOMP%NG+1)
         DO I=1,GCOMP%NG
          PTRAR(I,2) = int(IPE(I+1)-IPE(I))
         ENDDO
      ELSE
        LIW8 =  LIWALLOC
        NZG8 = NZ8
        IW   => IWALLOC(1:LIW8)
        IPE  => IPEALLOC(1:N+1)
      ENDIF
      LP    = ICNTL(1)
      MP    = ICNTL(3)
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      PROK  = ((MP.GT.0).AND.(ICNTL(4).GE.2))
      LDIAG = ICNTL(4)
      COMPRESS_SCHUR = .FALSE.
      IF (PROK) THEN
        IF (present(GCOMP)) THEN
         WRITE(MP,'(A,I10,A,I13,A)') " Processing a graph of size:", N
     &                          ," with ", GCOMP%NZG, " edges"
        ELSE
         WRITE(MP,'(A,I10)') " Processing a graph of size:", N
        ENDIF
      ENDIF
      IF (GCOMP_PROVIDED) THEN 
       FREE_CENTRALIZED_MATRIX = .FALSE.
      ELSE
       FREE_CENTRALIZED_MATRIX = (
     &                    (KEEP(54).EQ.3).AND.
     &                    (KEEP(494).EQ.0).AND.
     &                    (KEEP(106).NE.2)
     &                          )
      ENDIF
      INPLACE64_GRAPH_COPY    = .FALSE.
      INPLACE64_RESTORE_GRAPH = .TRUE.
      IF (KEEP(1).LT.0) KEEP(1) = 0
      NEMIN = KEEP(1)
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
        IF (present(SIZEOFBLOCKS)) THEN
          K = min(10,GCOMP%NG)
          IF (LDIAG.EQ.4) K = GCOMP%NG
          WRITE (MP,99909) N, NZG8, INFO(1)
          I8= 0_8
          WRITE(MP,'(A)') "  Graph adjacency "
          DO J=1, K
           IFIRST = GCOMP%IPE(J)
           ILAST= min(GCOMP%IPE(J+1)-1,GCOMP%IPE(J)+K-1)
           write(MP,'(A,I10)') "  .... node/column:", J 
           write(MP,'(8X,10I9)') 
     &       (GCOMP%ADJ(I8),I8=IFIRST,ILAST)
          ENDDO
        ELSE
          J8 = min(NZG8, 10_8)
          IF (LDIAG .EQ.4) J8 = NZG8
          WRITE (MP,99999) N, NZG8, LIW8, INFO(1)
          IF (J8.GT.0_8) WRITE (MP,99998) (IRN(I8),ICN(I8),I8=1_8,J8)
        ENDIF
         K = min0(10,N)
         IF (LDIAG.EQ.4) K = N
         IF (IORD.EQ.1 .AND. K.GT.0) THEN
            WRITE (MP,99997) (IKEEP1(I),I=1,K)
         ENDIF
      ENDIF
      NCMP    = N   
      IF (KEEP(60).NE.0) THEN
         IF ((SIZE_SCHUR.LE.0 ).OR.
     &        (SIZE_SCHUR.GE.N) ) GOTO 90 
      ENDIF
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
      IF  ( ( KEEP(60).NE.0).AND.(SIZE_SCHUR.GT.0)
     &     .AND.
     &     ((IORD.EQ.7).OR.(IORD.EQ.5))
     &     )THEN
         COMPRESS_SCHUR=.TRUE.
         NCMP          = N-SIZE_SCHUR
         ALLOCATE(IPQ8(N),stat=IERR)
         IF ( IERR .GT. 0 ) THEN
           INFO( 1 ) = -7
           INFO( 2 ) = N*KEEP(10)
         ENDIF
         CALL DMUMPS_GNEW_SCHUR(N,NCMP,NZ8,IRN(1), ICN(1), IW(1), LIW8,
     &        IPE(1), PTRAR(1,2),
     &        IPQ8, IWL1, IWFR8, KEEP8(126), KEEP8(127),
     &        INFO(1), INFO(2), ICNTL, symmetry, 
     &        KEEP(50), NBQD, AvgDens, 
     &        KEEP(264), KEEP(265),
     &        LISTVAR_SCHUR(1), SIZE_SCHUR, FRERE(1), FILS(1), 
     &        INPLACE64_GRAPH_COPY)
         INFO(8) = symmetry
         INPLACE64_GRAPH_COPY = INPLACE64_GRAPH_COPY.AND.
     &                          (.NOT.FREE_CENTRALIZED_MATRIX)
         DEALLOCATE(IPQ8)
         IORD = 5
         KEEP(95) = 1
         NBQD     = 0           
      ELSE
#endif
        IF (GCOMP_PROVIDED) THEN
          IWFR8 = GCOMP%NZG+1_8
        ELSE
          ALLOCATE(IPQ8(N),stat=IERR)
          IF ( IERR .GT. 0 ) THEN
           INFO( 1 ) = -7
           INFO( 2 ) = N*KEEP(10)
          ENDIF
          IF ( PROK ) THEN
            CALL MUMPS_SECDEB( TIMEB )
          ENDIF
          CALL DMUMPS_ANA_GNEW(N,NZ8,IRN(1), ICN(1), IW(1), LIW8,
     &        IPE(1), PTRAR(1,2),
     &        IPQ8, IWL1, IWFR8, KEEP8(126), KEEP8(127),
     &        INFO(1), INFO(2), ICNTL, symmetry, 
     &        KEEP(50), NBQD, AvgDens, KEEP(264), KEEP(265), 
     &        .TRUE., INPLACE64_GRAPH_COPY)
          INFO(8) = symmetry
          INPLACE64_GRAPH_COPY = INPLACE64_GRAPH_COPY.AND.
     &                           (.NOT.FREE_CENTRALIZED_MATRIX)
          DEALLOCATE(IPQ8)
        ENDIF
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
      ENDIF
#endif
      IF(NBQD .GT. 0) THEN
         IF( KEEP(50) .EQ. 2 .AND. ICNTL(12) .EQ. 0 ) THEN
            IF(KEEP(95) .NE. 1) THEN
               IF ( PROK ) 
     &              WRITE( MP,*) 
     &              'Compressed/constrained ordering set OFF'
               KEEP(95) = 1   
            ENDIF
         ENDIF
      ENDIF
      IF ( (KEEP(60).NE.0) .AND. (IORD.GT.1) .AND.
     &     .NOT. COMPRESS_SCHUR ) THEN
         IORD = 0               
      ENDIF 
      IF ( (KEEP(50).EQ.2)
     & .AND. (KEEP(95) .EQ. 3)
     & .AND. (IORD .EQ. 7) ) THEN  
        IORD = 2 
      ENDIF
      CALL MUMPS_SET_ORDERING( NORIG, KEEP, 
     &     KEEP(50), NSLAVES, IORD, 
     &     NBQD, AvgDens,
     &     PROK, MP )
      IF(KEEP(50) .EQ. 2) THEN
         IF(KEEP(95) .EQ. 3 .AND. IORD .NE. 2) THEN
            IF (PROK) WRITE(MP,*)
     &      'WARNING: DMUMPS_ANA_F constrained ordering not '//
     &      ' available with selected ordering. Move to' //
     &      ' compressed ordering.'
            KEEP(95) = 2
         ENDIF
      ELSE
         KEEP(95) = 1
      ENDIF
      MTRANS = KEEP(23)
      COMPRESS = KEEP(95) - 1
      IF(COMPRESS .GT. 0 .AND. KEEP(52) .EQ. -2) THEN
         IF(CNTL4 .GE. 0.0D0) THEN
            IF (KEEP(1).LE.8) THEN
               NEMIN = 16   
            ELSE
               NEMIN = 2*KEEP(1)
            ENDIF
            IF (PROK) 
     &           WRITE(MP,*) 'Setting static pivoting ON, COMPRESS =',
     &           COMPRESS
         ENDIF
      ENDIF
      IF(MTRANS .GT. 0 .AND. KEEP(50) .EQ. 2) THEN
         KEEP(23) = 0
      ENDIF
      IF (COMPRESS .EQ. 2) THEN
         IF (IORD.NE.2) THEN
            WRITE(*,*) "IORD not compatible with COMPRESS:",
     &           IORD, COMPRESS
            CALL MUMPS_ABORT()
         ENDIF
         CALL  DMUMPS_SET_CONSTRAINTS(
     &        N,PIV(1),FRERE(1),FILS(1),NFSIZ(1),IKEEP1(1),
     &        NCST,KEEP,KEEP8, ROWSCA(1)
     &        )
      ENDIF
      IF ( IORD .NE. 1 ) THEN
         IF (COMPRESS .GE. 1) THEN
            ALLOCATE(IPQ8(N),stat=IERR)
            IF ( IERR .GT. 0 ) THEN
               INFO( 1 ) = -7
               INFO( 2 ) = N*KEEP(10)
            ENDIF
            CALL DMUMPS_LDLT_COMPRESS(
     &           N, NZ8, IRN(1), ICN(1), PIV(1),
     &           NCMP, IW(1), LIW8, IPE(1), PTRAR(1,2), IPQ8,
     &           IWL1, FILS(1), IWFR8,
     &           IERROR, KEEP, KEEP8, ICNTL, INPLACE64_GRAPH_COPY)
            DEALLOCATE(IPQ8)
            symmetry = 100
         ENDIF
         IF ( (symmetry.LT.MINSYM).AND.(KEEP(50).EQ.0) ) THEN
            IF(KEEP(23) .EQ. 7 ) THEN
               KEEP(23) = -5
              GOTO 90 
            ELSE IF(KEEP(23) .EQ. -9876543) THEN
               IDENT = .TRUE.
               KEEP(23) = 5
               IF (PROK) WRITE(MP,'(A)')
     &              ' ... Apply column permutation (already computed)'
               DO J=1,N
                  JPERM = PIV(J)
                  FILS(JPERM) = J 
                  IF (JPERM.NE.J) IDENT = .FALSE.
               ENDDO
               IF (.NOT.IDENT) THEN
                  DO J8=1_8,NZ8
                     J = ICN(J8)
                     IF ((J.LE.0).OR.(J.GT.N)) CYCLE
                     ICN(J8) = FILS(J)
                  ENDDO
                  ALLOCATE(COLSCA_TEMP(N), stat=IERR)
                  IF ( IERR > 0 ) THEN
                     INFO( 1 ) = -7
                     INFO( 2 ) = N
                     GOTO 90    
                  ENDIF
                  DO J = 1, N
                     COLSCA_TEMP(J)=COLSCA(J)
                  ENDDO
                  DO J=1, N
                     COLSCA(FILS(J))=COLSCA_TEMP(J)
                  ENDDO
                  DEALLOCATE(COLSCA_TEMP)
                  IF (PROK)
     &                 WRITE(MP,'(/A)')
     &                 ' WARNING input matrix data modified'
                  ALLOCATE(IPQ8(N),stat=IERR)
                  IF ( IERR .GT. 0 ) THEN
                     INFO( 1 ) = -7
                     INFO( 2 ) = N*KEEP(10)
                  ENDIF
                  CALL DMUMPS_ANA_GNEW
     &                 (N,NZ8,IRN(1), ICN(1), IW(1), LIW8, 
     &                  IPE(1), PTRAR(1,2), 
     &                 IPQ8, IWL1, IWFR8, KEEP8(126), KEEP8(127), 
     &                 INFO(1), INFO(2), ICNTL, symmetry, KEEP(50),
     &                 NBQD, AvgDens, KEEP(264), KEEP(265), 
     &                 .TRUE.,INPLACE64_GRAPH_COPY)
                  INFO(8) = symmetry
                  DEALLOCATE(IPQ8)
                  NCMP = N
               ELSE
                  KEEP(23) = 0
               ENDIF
            ENDIF
         ELSE IF (KEEP(23) .EQ. 7 .OR. KEEP(23) .EQ. -9876543 ) THEN
            IF (PROK) WRITE(MP,'(A)')
     &           ' ... No column permutation'
            KEEP(23) = 0
         ENDIF
      ENDIF                     
      IF (FREE_CENTRALIZED_MATRIX
     &     .AND.COMPRESS.EQ.0.AND.(.NOT.COMPRESS_SCHUR)) THEN
         deallocate(IRN)
         NULLIFY(IRN)
         deallocate(ICN)
         NULLIFY(ICN)
      ENDIF
      INPLACE64_RESTORE_GRAPH = 
     &               INPLACE64_RESTORE_GRAPH.AND.(COMPRESS.NE.1)
      ALLOCATE( PARENT ( N ), stat = IERR )
      IF ( IERR .GT. 0 ) THEN
           INFO( 1 ) = -7
           INFO( 2 ) = N
           GOTO 90              
      ENDIF
      IF (IORD.NE.1 .AND. IORD.NE.5) THEN
         IF ( KEEP(60) .NE. 0 ) THEN
          IORD =  0
         ENDIF
         IF (PROK) THEN
            IF (IORD.EQ.2) THEN
               WRITE(MP,'(A)') ' Ordering based on AMF '
#if defined(scotch) || defined(ptscotch)
            ELSE IF (IORD.EQ.3) THEN
               WRITE(MP,'(A)') ' Ordering based on SCOTCH '
#endif
#if defined(pord)
            ELSE IF (IORD.EQ.4) THEN
               WRITE(MP,'(A)') ' Ordering based on PORD '
#endif
            ELSE IF (IORD.EQ.6) THEN
               WRITE(MP,'(A)') ' Ordering based on QAMD '
            ELSE
               WRITE(MP,'(A)') ' Ordering based on AMD '
            ENDIF
         ENDIF
         IF ( PROK ) THEN
            CALL MUMPS_SECDEB( TIMEB )
         ENDIF
         IF ( KEEP(60) .NE. 0 ) THEN
            CALL MUMPS_HAMD(N, LIW8, IPE(1), IWFR8, PTRAR(1,2), IW(1), 
     &           IWL1, IKEEP1(1), 
     &           IKEEP2(1), NCMPA, FILS(1), IKEEP3(1),
     &           PTRAR, PTRAR(1,3),
     &           PARENT,
     &           LISTVAR_SCHUR(1), SIZE_SCHUR)
            IF (KEEP(60)==1) THEN
               KEEP(20) = LISTVAR_SCHUR(1)
            ELSE
               KEEP(38) = LISTVAR_SCHUR(1)
            ENDIF
         ELSE
            IF ( .FALSE. ) THEN
#if defined(pord)
            ELSEIF (IORD .EQ. 4) THEN
               CALL MUMPS_PORD_INTSIZE(PORD_INT_SIZE)
               TOTW = N
               IF ( (COMPRESS .EQ. 1) 
     &               .OR.
     &            ( (NORIG.NE.N).AND.present(SIZEOFBLOCKS) )
     &                 ) THEN
                 IF (COMPRESS .EQ. 1) THEN
                   DO I=1,KEEP(93)/2
                      IWL1(I) = 2
                   ENDDO
                   DO I=1+KEEP(93)/2,NCMP
                      IWL1(I) = 1
                   ENDDO
                 ELSE IF 
     &             ( (NORIG.NE.N).AND.present(SIZEOFBLOCKS) ) THEN
                   TOTW = NORIG
                   DO I= 1, N
                     IWL1(I) = SIZEOFBLOCKS(I)
                   ENDDO
                 ENDIF
                  IF (PORD_INT_SIZE .EQ. 64) THEN
                    CALL MUMPS_PORDF_WND_MIXEDto64(NCMP, IWFR8-1_8,
     &                   IPE, IW,
     &                   IWL1, NCMPA, TOTW, PARENT, 
     &                   INFO(1), LP, LPOK, KEEP(10),
     &                   INPLACE64_GRAPH_COPY
     &                   )
                  ELSE IF (PORD_INT_SIZE .EQ. 32) THEN
                    CALL MUMPS_PORDF_WND_MIXEDto32(NCMP, IWFR8-1_8,
     &                   IPE, IW,
     &                   IWL1, NCMPA, TOTW, PARENT, 
     &                   INFO(1), LP, LPOK, KEEP(10))
                  ELSE
                     WRITE(*,*)
     &              "Internal error in PORD wrappers, PORD_INT_SIZE=",
     &              PORD_INT_SIZE
                     CALL MUMPS_ABORT()
                  ENDIF
                  IF ( NCMPA .NE. 0 ) THEN
                     write(6,*) ' Out PORD, NCMPA=', NCMPA
                     INFO( 1 ) = -9999
                     INFO( 2 ) = 4
                     GOTO 90 
                  ENDIF
                  IF (INFO(1) .LT.0) GOTO 90
                  IF (COMPRESS.EQ.1) THEN
                   CALL DMUMPS_GET_ELIM_TREE(NCMP,PARENT,IWL1,FILS(1))
                   CALL DMUMPS_GET_PERM_FROM_PE(NCMP,PARENT,IKEEP1(1),
     &                 FRERE(1),PTRAR(1,1))
                   DO I=1,NCMP
                     IKEEP2(IKEEP1(I))=I
                   ENDDO
                  ENDIF
               ELSE  
                  IF (PORD_INT_SIZE.EQ.64) THEN
                  CALL MUMPS_PORDF_MIXEDto64(NCMP, IWFR8-1_8, IPE,
     &                 IW, 
     &                 IWL1, NCMPA, PARENT, 
     &                 INFO(1), LP, LPOK, KEEP(10),
     &                 INPLACE64_GRAPH_COPY
     &                 )
                  ELSE IF (PORD_INT_SIZE.EQ.32) THEN
                  CALL MUMPS_PORDF_MIXEDto32(NCMP, IWFR8-1_8, IPE,
     &                 IW, 
     &                 IWL1, NCMPA, PARENT, 
     &                 INFO(1), LP, LPOK, KEEP(10))
                  ELSE
                    WRITE(*,*)
     &              "Internal error in PORD wrappers, PORD_INT_SIZE=",
     &              PORD_INT_SIZE
                    CALL MUMPS_ABORT()
                  ENDIF
               ENDIF
               IF ( NCMPA .NE. 0 ) THEN
                  write(6,*) ' Out PORD, NCMPA=', NCMPA
                  INFO( 1 ) = -9999
                  INFO( 2 ) = 4
                  GOTO 90 
               ENDIF
               IF (INFO(1) .LT. 0) GOTO 90
#endif    
#if defined(scotch) || defined(ptscotch)
            ELSEIF (IORD .EQ. 3) THEN
               CALL MUMPS_SCOTCH_INTSIZE(SCOTCH_INT_SIZE)
              IF ( (COMPRESS .EQ. 1) 
     &               .OR.
     &            ( (NORIG.NE.N).AND.present(SIZEOFBLOCKS) )
     &                 ) THEN
                 WEIGHTREQUESTED=1
                 IF (COMPRESS .EQ. 1) THEN
                   DO I=1,KEEP(93)/2
                      IWL1(I) = 2
                   ENDDO
                   DO I=1+KEEP(93)/2,NCMP
                      IWL1(I) = 1
                   ENDDO
                 ELSE IF 
     &             ( (NORIG.NE.N).AND.present(SIZEOFBLOCKS) ) THEN
                   DO I= 1, N
                     IWL1(I) = SIZEOFBLOCKS(I)
                   ENDDO
                 ENDIF
               ELSE
                 WEIGHTREQUESTED = 0
                 DO I= 1, N
                     IWL1(I) = 1
                 ENDDO
               ENDIF
               IF (SCOTCH_INT_SIZE.EQ.32) THEN
                IF (KEEP(10).EQ.1) THEN
                 INFO(1)  = -52
                 INFO(2) = 2
                ELSE
                 CALL MUMPS_SCOTCH_MIXEDto32(NCMP, 
     &              IWFR8-1_8, IPE,
     &              PARENT, IWFR8,
     &              PTRAR(1,2), IW, IWL1, IKEEP1,
     &              IKEEP2, NCMPA, INFO, LP, LPOK,
     &              WEIGHTUSED, WEIGHTREQUESTED)
                ENDIF
               ELSE IF (SCOTCH_INT_SIZE.EQ.64) THEN
                 CALL MUMPS_SCOTCH_MIXEDto64(NCMP, 
     &              IWFR8-1_8, IPE,
     &              PARENT, IWFR8,
     &              PTRAR(1,2), IW, IWL1, IKEEP1,
     &              IKEEP2, NCMPA, INFO, LP, LPOK, KEEP(10),
     &              INPLACE64_GRAPH_COPY,
     &              WEIGHTUSED, WEIGHTREQUESTED)
               ELSE
                 WRITE(*,*)
     &           "Internal error in SCOTCH wrappers, SCOTCH_INT_SIZE=",
     &           SCOTCH_INT_SIZE
                 CALL MUMPS_ABORT()
               ENDIF
               IF (INFO(1) .LT. 0) GOTO 90
               IF ( (COMPRESS .EQ. 1)
     &             .OR. 
     &            ( (NORIG.NE.N).AND.present(SIZEOFBLOCKS).AND.
     &              (WEIGHTUSED.EQ.0) )
     &            ) THEN
                 CALL DMUMPS_GET_ELIM_TREE(NCMP,PARENT,IWL1,FILS(1))
                 CALL DMUMPS_GET_PERM_FROM_PE(NCMP,PARENT,IKEEP1(1),
     &                FRERE(1),PTRAR(1,1))
                 DO I=1,NCMP
                   IKEEP2(IKEEP1(I))=I
                 ENDDO
               ENDIF
#endif
            ELSEIF (IORD .EQ. 2) THEN
               NBBUCK = 2*N
               COMPUTE_PERM=.FALSE.
               IF(COMPRESS .GE. 1) THEN
                  COMPUTE_PERM=.TRUE.
                  DO I=1,KEEP(93)/2
                     IWL1(I) = 2
                  ENDDO
                  DO I=1+KEEP(93)/2,NCMP
                     IWL1(I) = 1
                  ENDDO
                  TOTEL = KEEP(93)+KEEP(94)
               ELSE
                  IWL1(1) = -1
                  TOTEL   = N
               ENDIF
               IF (present(SIZEOFBLOCKS)) THEN
                IF (COMPRESS.GE.1) THEN
                 CALL MUMPS_ABORT()
                ENDIF
                NBBUCK = max(NBBUCK, NORIG-N)
                NBBUCK = max(NBBUCK, 2*NORIG)
                NCMP  = N 
                TOTEL = NORIG
                DO I= 1, N
                  IWL1(I) = SIZEOFBLOCKS(I)
                ENDDO
               ENDIF
               ALLOCATE( WTEMP ( 0: NBBUCK + 1), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = NBBUCK+2
                  GOTO 90       
               ENDIF
               IF(COMPRESS .LE. 1) THEN
                  CALL MUMPS_HAMF4
     &               (TOTEL, NCMP, COMPUTE_PERM, NBBUCK, LIW8, IPE(1),
     &               IWFR8, PTRAR(1,2),
     &               IW(1), IWL1, IKEEP1(1), IKEEP2(1), NCMPA, FILS(1),
     &               IKEEP3(1), PTRAR, PTRAR(1,3), WTEMP, PARENT(1))
               ELSE
                  IF(PROK) WRITE(MP,'(A)') 
     &                 ' Constrained Ordering based on AMF'
                  CALL MUMPS_CST_AMF(NCMP, NBBUCK, LIW8, IPE(1),
     &             IWFR8, PTRAR(1,2), 
     &             IW(1), IWL1, IKEEP1(1), IKEEP2(1), NCMPA, FILS(1), 
     &             IKEEP3(1), PTRAR, PTRAR(1,3), WTEMP,
     &             NFSIZ(1), FRERE(1), PARENT(1))
               ENDIF
               DEALLOCATE(WTEMP)
            ELSEIF (IORD .EQ. 6) THEN
               ALLOCATE( WTEMP ( N ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = N
                  GOTO 90       
               ENDIF
               THRESH = 1
               IVersion = 2
               COMPUTE_PERM=.FALSE.
               IF(COMPRESS .EQ. 1) THEN
                  COMPUTE_PERM=.TRUE.
                  DO I=1,KEEP(93)/2
                     IWL1(I) = 2
                  ENDDO
                  DO I=1+KEEP(93)/2,NCMP
                     IWL1(I) = 1
                  ENDDO
                  TOTEL = KEEP(93)+KEEP(94)
               ELSE
                  IWL1(1) = -1
                  TOTEL = N
               ENDIF
               IF (present(SIZEOFBLOCKS)) THEN
                IF (COMPRESS.EQ.1) THEN
                 CALL MUMPS_ABORT()
                ENDIF
                NCMP  = N 
                TOTEL = NORIG
                DO I= 1, N
                  IWL1(I) = SIZEOFBLOCKS(I)
                ENDDO
               ENDIF
               CALL MUMPS_QAMD
     &              (TOTEL,COMPUTE_PERM,IVersion, THRESH, WTEMP,
     &              NCMP, LIW8, IPE(1), IWFR8, PTRAR(1,2), IW(1),
     &              IWL1, IKEEP1(1), IKEEP2(1), NCMPA, FILS(1),
     &              IKEEP3(1), PTRAR, PTRAR(1,3), PARENT(1))
               DEALLOCATE(WTEMP)
            ELSE
               COMPUTE_PERM=.FALSE.
               IF(COMPRESS .EQ. 1) THEN
                  COMPUTE_PERM=.TRUE.
                  DO I=1,KEEP(93)/2
                     IWL1(I) = 2
                  ENDDO
                  DO I=1+KEEP(93)/2,NCMP
                     IWL1(I) = 1
                  ENDDO
                  TOTEL = KEEP(93)+KEEP(94)
               ELSE
                  IWL1(1) = -1
                  TOTEL = N
               ENDIF
               IF (present(SIZEOFBLOCKS)) THEN
                IF (COMPRESS.EQ.1) THEN
                 CALL MUMPS_ABORT()
                ENDIF
                NCMP  = N 
                TOTEL = NORIG
                DO I= 1, N
                  IWL1(I) = SIZEOFBLOCKS(I)
                ENDDO
               ENDIF
               CALL MUMPS_ANA_H(TOTEL, COMPUTE_PERM,
     &              NCMP, LIW8, IPE(1), IWFR8, PTRAR(1,2),
     &              IW(1), IWL1, IKEEP1(1), IKEEP2(1), NCMPA, FILS(1),
     &              IKEEP3(1), PTRAR, PTRAR(1,3), PARENT(1))
            ENDIF
         ENDIF
         IF(COMPRESS .GE. 1) THEN
            CALL DMUMPS_EXPAND_PERMUTATION(N,NCMP,KEEP(94),KEEP(93),
     &           PIV(1),IKEEP1(1),IKEEP2(1))
            COMPRESS = -1
         ENDIF
         IF ( PROK ) THEN
          CALL MUMPS_SECFIN( TIMEB )
#if  defined(scotch) || defined(ptscotch)
          IF (IORD.EQ.3) THEN
            WRITE( MP, '(A,F12.4)' )
     &        ' ELAPSED TIME SPENT IN SCOTCH reordering =', TIMEB
          ENDIF
#endif
         ENDIF
      ENDIF  
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
      IF (IORD.EQ.5) THEN
         IF (PROK) THEN
            WRITE(MP,'(A)') ' Ordering based on METIS'
         ENDIF
         IF ( PROK ) THEN
            CALL MUMPS_SECDEB( TIMEB )
         ENDIF
         CALL MUMPS_METIS_IDXSIZE(METIS_IDX_SIZE)
         IF (KEEP(10).EQ.1.AND.METIS_IDX_SIZE.NE.64) THEN
           INFO(1) = -52
           INFO(2) = 1
           GOTO 90 
         ENDIF
#if defined(metis4) || defined(parmetis3)
         NUMFLAG = 1
         OPT_METIS_SIZE = 8
#else
         OPT_METIS_SIZE = 40
#endif
         IF (COMPRESS .EQ. 1) THEN
            DO I=1,KEEP(93)/2
               FRERE(I) = 2
            ENDDO
            DO I=KEEP(93)/2+1,NCMP
               FRERE(I) = 1
            ENDDO
#if defined(metis4) || defined(parmetis3)
            IF (METIS_IDX_SIZE .EQ.32) THEN
              CALL MUMPS_METIS_NODEWND_MIXEDto32(
     &           NCMP, IPE, IW, FRERE, 
     &           NUMFLAG, METIS_OPTIONS(1), OPT_METIS_SIZE,
     &           IKEEP2, IKEEP1, INFO(1), LP, LPOK )
            ELSE IF (METIS_IDX_SIZE .EQ.64) THEN
              CALL MUMPS_METIS_NODEWND_MIXEDto64(
     &           NCMP, IPE, IW, FRERE, 
     &           NUMFLAG, METIS_OPTIONS(1), OPT_METIS_SIZE,
     &           IKEEP2, IKEEP1, INFO(1), LP, LPOK, KEEP(10),
     &           INPLACE64_GRAPH_COPY )
            ELSE
              WRITE(*,*)
     &        "Internal error in METIS wrappers, METIS_IDX_SIZE=",
     &        METIS_IDX_SIZE
              CALL MUMPS_ABORT()
            ENDIF
         ELSE
           IF ((NORIG.NE.N).AND.present(SIZEOFBLOCKS)) THEN
            DO I=1, N
             FRERE(I) = SIZEOFBLOCKS(I)
            ENDDO
            IF (METIS_IDX_SIZE .EQ.32) THEN
              CALL MUMPS_METIS_NODEWND_MIXEDto32(
     &           NCMP, IPE, IW, FRERE, 
     &           NUMFLAG, METIS_OPTIONS(1), OPT_METIS_SIZE,
     &           IKEEP2, IKEEP1, INFO(1), LP, LPOK )
            ELSE IF (METIS_IDX_SIZE .EQ.64) THEN
              CALL MUMPS_METIS_NODEWND_MIXEDto64(
     &           NCMP, IPE, IW, FRERE, 
     &           NUMFLAG, METIS_OPTIONS(1), OPT_METIS_SIZE,
     &           IKEEP2, IKEEP1, INFO(1), LP, LPOK, KEEP(10),
     &           INPLACE64_GRAPH_COPY )
            ELSE
              WRITE(*,*)
     &        "Internal error in METIS wrappers, METIS_IDX_SIZE=",
     &        METIS_IDX_SIZE
              CALL MUMPS_ABORT()
            ENDIF
           ELSE
             IF (METIS_IDX_SIZE .EQ.32) THEN
              CALL MUMPS_METIS_NODEND_MIXEDto32(
     &           NCMP, IPE, IW, NUMFLAG, 
     &           METIS_OPTIONS(1), OPT_METIS_SIZE,
     &           IKEEP2, IKEEP1, INFO(1), LP, LPOK )
             ELSE IF (METIS_IDX_SIZE .EQ.64) THEN
              CALL MUMPS_METIS_NODEND_MIXEDto64(
     &           NCMP, IPE, IW, NUMFLAG, 
     &           METIS_OPTIONS(1), OPT_METIS_SIZE,
     &           IKEEP2, IKEEP1, INFO(1), LP,LPOK,KEEP(10), 
     &           LIW8, INPLACE64_GRAPH_COPY,
     &           INPLACE64_RESTORE_GRAPH)
             ELSE
              WRITE(*,*)
     &        "Internal error in METIS wrappers, METIS_IDX_SIZE=",
     &        METIS_IDX_SIZE
              CALL MUMPS_ABORT()
           ENDIF
           ENDIF
         ENDIF
#else
         ELSE
            IF (present(SIZEOFBLOCKS)) THEN
             DO I=1,N  
              FRERE(I) = SIZEOFBLOCKS(I)
             ENDDO
            ELSE
             DO I=1,NCMP
               FRERE(I) = 1
             ENDDO 
            ENDIF
         ENDIF
         IF (METIS_IDX_SIZE .EQ. 32) THEN
           CALL MUMPS_METIS_NODEND_MIXEDto32(
     &        NCMP, IPE, IW, FRERE,
     &        METIS_OPTIONS(1), OPT_METIS_SIZE,
     &        IKEEP2, IKEEP1, INFO(1), LP, LPOK )
         ELSE IF (METIS_IDX_SIZE .EQ. 64) THEN
           CALL MUMPS_METIS_NODEND_MIXEDto64(
     &        NCMP, IPE, IW, FRERE,
     &        METIS_OPTIONS(1), OPT_METIS_SIZE,
     &        IKEEP2, IKEEP1, INFO(1), LP,LPOK,KEEP(10),
     &        LIW8, INPLACE64_GRAPH_COPY,
     &        INPLACE64_RESTORE_GRAPH)
         ELSE
              IF (LPOK) WRITE(LP,*)
     &        "Internal error in METIS wrappers, METIS_IDX_SIZE=",
     &        METIS_IDX_SIZE
              CALL MUMPS_ABORT()
         ENDIF
#endif
         IF (INFO(1) .LT.0) GOTO 90
         IF ( PROK ) THEN
            CALL MUMPS_SECFIN( TIMEB )
            WRITE( MP, '(A,F12.4)' )
     &        ' ELAPSED TIME SPENT IN METIS reordering  =', TIMEB
         ENDIF
         IF ( COMPRESS_SCHUR ) THEN
            CALL DMUMPS_EXPAND_PERM_SCHUR(
     &           N, NCMP, IKEEP1(1),IKEEP2(1),
     &           LISTVAR_SCHUR(1), SIZE_SCHUR, FILS(1))
            COMPRESS = -1       
         ENDIF
         IF (COMPRESS .EQ. 1) THEN
            CALL DMUMPS_EXPAND_PERMUTATION(N,NCMP,KEEP(94),
     &           KEEP(93),PIV(1),IKEEP1(1),IKEEP2(1))
            COMPRESS = -1       
         ENDIF
      ENDIF                     
#endif
      IF (PROK) THEN
         IF (IORD.EQ.1) THEN
            WRITE(MP,'(A)') ' Ordering given is used'
         ENDIF
      ENDIF
      IF (IORD.EQ.1 .OR. IORD.EQ.5 .OR. COMPRESS.EQ.-1 
     &    .OR.  
     &   ( (NORIG.NE.N).AND.present(SIZEOFBLOCKS) .AND.(IORD.EQ.3)
     &     .AND. (WEIGHTUSED.EQ.0)
     &   )
     &   ) THEN
         IF ((KEEP(106).EQ.1).OR.(KEEP(106).EQ.3)
     &          .OR.(KEEP(60).NE.0)) THEN
            IF ( COMPRESS .EQ. -1 ) THEN
               ALLOCATE(IPQ8(N),stat=IERR)
               IF ( IERR .GT. 0 ) THEN
                     INFO( 1 ) = -7
                     INFO( 2 ) = N*KEEP(10)
               ENDIF
               CALL DMUMPS_ANA_GNEW(N,NZ8,IRN(1),ICN(1),IW(1),LIW8,
     &              IPE(1), PTRAR(1,2),
     &              IPQ8, IWL1, IWFR8, KEEP8(126), KEEP8(127),
     &              INFO(1), INFO(2), ICNTL, symmetry, KEEP(50),
     &              NBQD, AvgDens, KEEP(264),KEEP(265), .TRUE.,
     &              INPLACE64_GRAPH_COPY)
               DEALLOCATE(IPQ8)
            ENDIF
            COMPRESS = 0
           IF (KEEP(106).EQ.3.AND.KEEP(60).EQ.0) THEN
           ELSE
            ALLOCATE( WTEMP ( 2*N ), stat = IERR )
            IF ( IERR .GT. 0 ) THEN
               INFO( 1 ) = -7
               INFO( 2 ) = 2*N
               GOTO 90          
            ENDIF
            THRESH = -1
            IF (KEEP(60) == 0) THEN
               ITEMP = 0 
            ELSE
               ITEMP = SIZE_SCHUR
               IF (KEEP(60)==1) THEN
                  KEEP(20) = LISTVAR_SCHUR(1)
               ELSE
                  KEEP(38) = LISTVAR_SCHUR(1)
               ENDIF
            ENDIF
            AGG6 =.FALSE.
            IF (present(SIZEOFBLOCKS)) THEN
             DO I=1, N
              IWL1(I) = SIZEOFBLOCKS(I)
             ENDDO
             TOTEL = NORIG
            ELSE
             IWL1(1) = -1
             TOTEL = N
            ENDIF
            CALL MUMPS_SYMQAMD(THRESH, WTEMP,
     &           N, TOTEL, LIW8, IPE(1), IWFR8, PTRAR(1,2), IW(1),
     &           IWL1(1), WTEMP(N+1),
     &           IKEEP2(1), NCMPA, FILS(1), IKEEP3(1), PTRAR, 
     &           PTRAR(1,3),IKEEP1(1), LISTVAR_SCHUR(1), ITEMP, 
     &           AGG6, PARENT)
            DEALLOCATE(WTEMP)
           ENDIF  
         ELSE
            CALL DMUMPS_ANA_J(N, NZ8, IRN(1), ICN(1), IKEEP1(1), IW(1),
     &           LIW8, IPE(1),
     &           PTRAR(1,2), IWL1, IWFR8,
     &           INFO(1),INFO(2), MP)
            IF (KEEP(60) .EQ. 0) THEN
               ITEMP = 0 
            ELSE 
               ITEMP = SIZE_SCHUR
            ENDIF
            CALL DMUMPS_ANA_K(N, IPE(1), IW(1), LIW8, IWFR8, IKEEP1(1),
     &              IKEEP2(1), IWL1,
     &              PTRAR, NCMPA, ITEMP, PARENT)
            IF (KEEP(60) .EQ. 0) THEN
               IF (KEEP(60) .EQ. 1) THEN
                  KEEP(20) = LISTVAR_SCHUR(1)
               ELSE
                  KEEP(38) = LISTVAR_SCHUR(1)
               ENDIF
            ENDIF
         ENDIF                  
      ENDIF                     
#if defined(OLDDFS)
      CALL DMUMPS_ANA_L
     &     (N, PARENT, IWL1, IKEEP1(1), IKEEP2(1), IKEEP3(1),
     &     NFSIZ, INFO(6), FILS(1), FRERE(1), PTRAR(1,3),
     &     NEMIN, KEEP(60))
#else
      IF (allocated(IPEALLOC)) DEALLOCATE(IPEALLOC)
      ALLOCATE(WTEMP(N), stat=IERR)
      IF ( IERR .GT. 0 ) THEN
         INFO( 1 ) = -7
         INFO( 2 ) = N
         GOTO 90                
      ENDIF
      IF (present(SIZEOFBLOCKS)) THEN
       CALL DMUMPS_ANA_LNEW
     &     (N, PARENT, IWL1, IKEEP1(1), IKEEP2(1), IKEEP3(1),
     &     NFSIZ(1), PTRAR, INFO(6), FILS(1), FRERE(1), 
     &     PTRAR(1,3), NEMIN, WTEMP, KEEP(60),
     &     KEEP(20),KEEP(38),PTRAR(1,2),KEEP(104),IW(1),KEEP(50), 
     &     ICNTL(13), KEEP(37), NSLAVES, KEEP(250).EQ.1 
     &     , .TRUE. , SIZEOFBLOCKS, N
     &      )
      ELSE
       CALL DMUMPS_ANA_LNEW
     &     (N, PARENT, IWL1, IKEEP1(1), IKEEP2(1), IKEEP3(1),
     &     NFSIZ(1), PTRAR, INFO(6), FILS(1), FRERE(1), 
     &     PTRAR(1,3), NEMIN, WTEMP, KEEP(60),
     &     KEEP(20),KEEP(38),PTRAR(1,2),KEEP(104),IW(1),KEEP(50), 
     &     ICNTL(13), KEEP(37), NSLAVES, KEEP(250).EQ.1 
     &     , .FALSE., IDUMMY, LIDUMMY )
      ENDIF
      DEALLOCATE(WTEMP)
#endif
      IF (KEEP(60).NE.0)  THEN
         IF (KEEP(60)==1) THEN
            IN = KEEP(20)
         ELSE
            IN = KEEP(38)
         ENDIF
         DO WHILE (IN.GT.0) 
            IN = FILS (IN)
         END DO
         IFSON = -IN
         IF (KEEP(60)==1) THEN
            IN = KEEP(20)
         ELSE
            IN = KEEP(38)
         ENDIF
         DO I=2,SIZE_SCHUR
            FILS(IN) = LISTVAR_SCHUR (I)
            IN       = FILS(IN)
            FRERE (IN) = N+1
         ENDDO
         FILS(IN) = -IFSON
      ENDIF
      CALL DMUMPS_ANA_M(IKEEP2(1),
     &     PTRAR(1,3), INFO(6),
     &     INFO(5), KEEP(2), KEEP(50),
     &     KEEP8(101), KEEP(108), KEEP(5),
     &     KEEP(6), KEEP(226), KEEP(253))
      KEEP(59) = INFO(5)
      IF ( KEEP(53) .NE. 0 ) THEN
         CALL MUMPS_MAKE1ROOT( N, FRERE(1), FILS(1), NFSIZ(1),
     &                         KEEP(20) )
      END IF
      IF (  (KEEP(48) == 4 .AND. KEEP8(21).GT.0_8)
     &     .OR.
     &     (KEEP (48)==5 .AND. KEEP8(21) .GT. 0_8 )
     &     .OR.
     &     (KEEP(24).NE.0.AND.KEEP8(21).GT.0_8) ) THEN 
         CALL DMUMPS_SET_K821_SURFACE(KEEP8(21), KEEP(2),
     &        KEEP(48), KEEP(50), NSLAVES)
      END IF
      IF (KEEP(210).LT.0.OR.KEEP(210).GT.2) THEN
         KEEP(210)=0
      ENDIF
      IF (KEEP(210).EQ.0.AND.KEEP(201).GT.0) THEN
         KEEP(210)=1            
      ENDIF
      IF (KEEP(210).EQ.0.AND.KEEP(201).EQ.0) THEN
         KEEP(210)=2            
      ENDIF
      IF (KEEP(210).EQ.2) THEN
         KEEP8(79)=huge(KEEP8(79))
      ENDIF
      IF (KEEP(210).EQ.1.AND.KEEP8(79).LE.0_8) THEN
        KEEP8(79)=K79REF * int(NSLAVES,8)
      ENDIF
      IF ( (KEEP(79).EQ.0).OR.(KEEP(79).EQ.2).OR.
     &     (KEEP(79).EQ.3).OR.(KEEP(79).EQ.5).OR.
     &     (KEEP(79).EQ.6)
     &   )  THEN
       IF (KEEP(210).EQ.1) THEN
        SPLITROOT = .FALSE. 
        IF ( KEEP(62).GE.1) THEN
          IWL1(1) = -1
          IF (present(SIZEOFBLOCKS)) THEN
            DO I= 1, N
               IWL1(I) = SIZEOFBLOCKS(I)
            ENDDO
          ENDIF
          CALL DMUMPS_CUTNODES(N, FRERE(1), FILS(1), NFSIZ(1), 
     &                       IWL1(1), N, INFO(6),
     &                       NSLAVES, KEEP,KEEP8, SPLITROOT,
     &                       MP, LDIAG, INFO(1), INFO(2))
          IF (INFO(1).LT.0) GOTO 90 
          IF (PROK) THEN
               WRITE(MP,*) " Number of split nodes in pre-splitting=",
     &         KEEP(61)
          ENDIF
        ENDIF
       ENDIF
      ENDIF
      SPLITROOT = ((ICNTL(13).GT.0 .AND. NSLAVES.GT.ICNTL(13)) .OR.
     &     ICNTL(13).EQ.-1 )
      IF (KEEP(53) .NE. 0) THEN
         SPLITROOT = .TRUE.
      ENDIF
      SPLITROOT = (SPLITROOT.AND.( (KEEP(60).EQ.0) ))
      IF (SPLITROOT) THEN
         IWL1(1) = -1
         IF (present(SIZEOFBLOCKS)) THEN
            DO I= 1, N
               IWL1(I) = SIZEOFBLOCKS(I)
            ENDDO
         ENDIF
         CALL DMUMPS_CUTNODES(N, FRERE(1), FILS(1), NFSIZ(1), 
     &                    IWL1(1), N, INFO(6),
     &                    NSLAVES, KEEP,KEEP8, SPLITROOT,
     &                    MP, LDIAG, INFO(1), INFO(2))
         IF (INFO(1).LT.0) GOTO 90 
         IF ( KEEP(53) .NE. 0 ) THEN
          CALL MUMPS_MAKE1ROOT( N, FRERE(1), FILS(1), NFSIZ(1),
     &                          KEEP(20) )
         ENDIF
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
         K = min0(10,N)
         IF (LDIAG.EQ.4) K = N
         IF (K.GT.0) WRITE (MP,99987) (NFSIZ(I),I=1,K)
         IF (K.GT.0) WRITE (MP,99989) (FILS(I),I=1,K)
         IF (K.GT.0) WRITE (MP,99988) (FRERE(I),I=1,K)
      ENDIF
      GO TO 90
      INFO(1) = -4
      INFO(2) = K
      GOTO 90 
 90   CONTINUE 
      IF (INFO(1) .NE. 0) THEN
        IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) 
     &           WRITE (LP,99996) INFO(1), INFO(2)
      ENDIF
      IF (allocated(IWALLOC))  DEALLOCATE(IWALLOC)
      IF (allocated(IWL1))     DEALLOCATE(IWL1)
      IF (allocated(IPEALLOC)) DEALLOCATE(IPEALLOC)
      IF (allocated(PTRAR))    DEALLOCATE(PTRAR)
      IF (allocated(PARENT))   DEALLOCATE(PARENT)
      RETURN
99999 FORMAT (/'Entering ordering phase with ...'/
     &  '                N        NNZ       LIW       INFO(1)'/,
     &     6X, I10, I11, I12, I10)
99998 FORMAT ('Matrix entries:    IRN()   ICN()'/
     &     (I12, I9, I12, I9, I12, I9))
99909 FORMAT (/'Entering ordering phase with graph dimensions ...'/
     &  '              |V|        |E|      INFO(1)'/,
     &     10X, I10, I13, I10)
99997 FORMAT ('IKEEP1(.)=', 10I8/(12X, 10I8))
99996 FORMAT 
     &  (/'** Error/warning return ** from Analysis *  INFO(1:2)= ', 
     &       (I3, I16))
99989 FORMAT ('FILS (.)  =', 10I9/(11X, 10I9))
99988 FORMAT ('FRERE(.)  =', 10I9/(11X, 10I9))
99987 FORMAT ('NFSIZ(.)  =', 10I9/(11X, 10I9))
      END SUBROUTINE DMUMPS_ANA_F
      SUBROUTINE DMUMPS_ANA_N_DIST( id, PTRAR )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_STRUC
      IMPLICIT NONE
      include 'mpif.h'
      TYPE(DMUMPS_STRUC), INTENT(INOUT), TARGET :: id
      INTEGER(8), INTENT(OUT), TARGET        :: PTRAR(:)
      INTEGER          :: IERR, allocok
      INTEGER          :: IOLD, JOLD, INEW, JNEW
      INTEGER(8)       :: K, INZ
      INTEGER, POINTER :: IIRN(:), IJCN(:)
      INTEGER(8), POINTER :: IWORK1(:), IWORK2(:)
      LOGICAL          :: IDO
      IF(id%KEEP(54) .EQ. 3) THEN
            IIRN => id%IRN_loc
            IJCN => id%JCN_loc
            INZ = id%KEEP8(29)
            IWORK1 => PTRAR(id%N+1:id%N+id%N)
            allocate(IWORK2(id%N),stat=allocok)
            IF (allocok > 0 ) THEN
               id%INFO(1) = -7
               id%INFO(2) = id%N
               RETURN
            ENDIF
            IDO = .TRUE.
      ELSE
            IIRN => id%IRN
            IJCN => id%JCN
            INZ  =  id%KEEP8(28)
            IWORK1 => PTRAR(1:id%N)
            IWORK2 => PTRAR(id%N+1:id%N+id%N)
            IDO = id%MYID .EQ. 0
      END IF
      DO 50 IOLD=1,id%N
         IWORK1(IOLD) = 0_8
         IWORK2(IOLD) = 0_8
 50   CONTINUE
      IF(IDO) THEN
         DO 70 K=1_8,INZ
            IOLD = IIRN(K)
            JOLD = IJCN(K)
            IF ( (IOLD.GT.id%N).OR.(JOLD.GT.id%N).OR.(IOLD.LT.1)
     &           .OR.(JOLD.LT.1) ) GOTO 70
            IF (IOLD.NE.JOLD) THEN
               INEW = id%SYM_PERM(IOLD)
               JNEW = id%SYM_PERM(JOLD)
               IF ( id%KEEP( 50 ) .EQ. 0 ) THEN
                  IF (INEW.LT.JNEW) THEN
                     IWORK2(IOLD) = IWORK2(IOLD) + 1_8
                  ELSE
                     IWORK1(JOLD) = IWORK1(JOLD) + 1_8
                  ENDIF
               ELSE
                  IF ( INEW .LT. JNEW ) THEN
                     IWORK1( IOLD ) = IWORK1( IOLD ) + 1_8
                  ELSE 
                     IWORK1( JOLD ) = IWORK1( JOLD ) + 1_8
                  END IF
               ENDIF
            ENDIF
 70      CONTINUE
      END IF
      IF (id%KEEP(54) .EQ. 3) THEN
         CALL MPI_ALLREDUCE(IWORK1(1), PTRAR(1), id%N,
     &        MPI_INTEGER8, MPI_SUM, id%COMM, IERR )
         CALL MPI_ALLREDUCE(IWORK2(1), PTRAR(id%N+1), id%N,
     &        MPI_INTEGER8, MPI_SUM, id%COMM, IERR )
         deallocate(IWORK2)
      ELSE
         CALL MPI_BCAST( PTRAR(1), 2*id%N, MPI_INTEGER8,
     &        0, id%COMM, IERR )
      END IF
      RETURN
      END SUBROUTINE DMUMPS_ANA_N_DIST
      SUBROUTINE DMUMPS_ANA_O( N, NZ, MTRANS, PERM, IKEEPALLOC,
     &     idIRN, idJCN, idA, idROWSCA, idCOLSCA, WORK2, KEEP,
     &     ICNTL, INFO, INFOG )
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N
      INTEGER(8), INTENT(IN) :: NZ 
      INTEGER, INTENT(OUT)   :: PERM(:) 
      INTEGER, POINTER, DIMENSION(:) :: idIRN, idJCN
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: idA
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: idROWSCA, idCOLSCA
      INTEGER, TARGET        :: IKEEPALLOC(3*N)
      INTEGER, INTENT(INOUT) :: MTRANS
      INTEGER                :: KEEP(500)
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER, INTENT(INOUT) :: INFOG(80) 
      INTEGER, TARGET        :: WORK2(N)
      INTEGER    :: allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IW
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: S2
      TARGET :: S2
      INTEGER ICNTL64(10), INFO64(10)
      INTEGER ICNTL_SYM_MWM(10),INFO_SYM_MWM(10)
      DOUBLE PRECISION CNTL64(10)
      INTEGER MPRINT,LP, MP       
      INTEGER JPERM
      INTEGER NUMNZ, I, J, JPOS
      LOGICAL PROK, IDENT, DUPPLI
      INTEGER K50, KER_SIZE, NZER_DIAG, MTRANSLOC,RZ_DIAG
      INTEGER(8) :: LIWG
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IPE
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IPQ8
      INTEGER    :: LSC
      INTEGER(8) :: NZTOT, NZREAL, IPIW, LIW, LIWMIN, NZsave, 
     &              K, KPOS, LDW, LDWMIN, IRNW, RSPOS, CSPOS,
     &              LS2,J8, N8
      LOGICAL SCALINGLOC
      INTEGER,POINTER,DIMENSION(:) :: ZERODIAG
      INTEGER,POINTER,DIMENSION(:) :: STR_KER
      INTEGER,POINTER,DIMENSION(:) :: MARKED
      INTEGER,POINTER,DIMENSION(:) :: FLAG
      INTEGER,POINTER,DIMENSION(:) :: PIV_OUT
      DOUBLE PRECISION THEMIN, THEMAX, COLNORM,MAXDBL, ABSAK
      DOUBLE PRECISION ZERO,TWO,ONE
      PARAMETER(ZERO = 0.0D0,TWO = 2.0D0,ONE = 1.0D0)
      N8     = int(N,8)
      MPRINT = ICNTL(3)
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      PROK   = ((MPRINT.GT.0).AND.(ICNTL(4).GE.2))
      K50 = KEEP(50)
      SCALINGLOC = .FALSE.
      IF(KEEP(52) .EQ. -2) THEN
         IF(.not.associated(idA)) THEN
         ELSE
            SCALINGLOC = .TRUE.
         ENDIF
      ELSE IF(KEEP(52) .EQ. 77) THEN
         SCALINGLOC = .TRUE.
         IF( MTRANS .NE. 5 .AND. MTRANS .NE. 6 
     &           .AND. MTRANS .NE. 7) THEN
               SCALINGLOC = .FALSE.
         ENDIF
         IF(.not.associated(idA)) THEN
            SCALINGLOC = .FALSE.
            IF (PROK) 
     &           WRITE(MPRINT,*) 'Analysis: auto scaling OFF because ',
     &                           'A not provided at analysis '
         ENDIF
      ENDIF
      IF ( (KEEP(50).EQ.2).AND.(ICNTL(8).NE.-2).AND.
     &    (MTRANS .EQ. 7 .OR. KEEP(95) .EQ. 0) ) THEN
        ZERODIAG => IKEEPALLOC(1:N)
        ZERODIAG = 0
        NZER_DIAG = N
        RZ_DIAG = 0
        DO K=1,NZ
            I = idIRN(K)
            J = idJCN(K)
            IF (I.NE.J) CYCLE
            IF ( (J.LE.N).AND.(J.GE.1) ) THEN
              IF(ZERODIAG(I) .EQ. 0) THEN
                ZERODIAG(I) = 1
                IF(associated(idA)) THEN
                  ABSAK= abs(idA(K))
                  IF(ABSAK .EQ. dble(0.0D0)) THEN
                     RZ_DIAG = RZ_DIAG + 1
                  ENDIF
                ENDIF
                NZER_DIAG = NZER_DIAG - 1                     
              ENDIF
            ENDIF
        ENDDO
        IF( (NZER_DIAG+RZ_DIAG) .LT. (N/10) ) THEN
          MTRANS = 0
          KEEP(95) =1
          GOTO 500
        ENDIF
      ENDIF
      IF(SCALINGLOC) THEN
         IF (PROK) WRITE(MPRINT,*) 
     &        'Scaling will be computed during analysis'
      ENDIF
      IF( MTRANS.NE.0 .AND. (.NOT.associated(idA)) ) MTRANS=1
      MTRANSLOC = MTRANS
      IF (MTRANS.LT.0 .OR. MTRANS.GT.7) GO TO 500
      IF (K50 .EQ. 0) THEN
         IF(.NOT. SCALINGLOC .AND. MTRANS .EQ. 7) THEN 
            GO TO 500
         ENDIF
         IF(SCALINGLOC) THEN
            IF (MTRANSLOC.NE.6) THEN
              MTRANSLOC = 5
            ENDIF
         ENDIF
      ELSE
         IF (MTRANS .EQ. 7) MTRANSLOC = 5
      ENDIF
      IF(SCALINGLOC .AND. MTRANSLOC .NE. 5 .AND.
     &     MTRANSLOC .NE. 6 ) THEN
         IF (PROK) WRITE(MPRINT,*)
     &        'WARNING scaling required: set MTRANS option to 5'
         MTRANSLOC = 5
      ENDIF
      IF (N.EQ.1) THEN
        MTRANS=0
        GO TO 500
      ENDIF
      IF(K50 .NE. 0) THEN
         NZTOT = 2_8*NZ+N8
      ELSE
         NZTOT = NZ
      ENDIF
      ZERODIAG => IKEEPALLOC(1:N) 
      STR_KER => IKEEPALLOC(N+1:2*N)
      CALL DMUMPS_MTRANSI(ICNTL64,CNTL64)
      ICNTL64(1) = ICNTL(1)
      ICNTL64(2) = ICNTL(2)
      ICNTL64(3) = ICNTL(3)
      ICNTL64(4) = -1
      IF (ICNTL(4).EQ.3) ICNTL64(4) = 0
      IF (ICNTL(4).EQ.4) ICNTL64(4) = 1
      ICNTL64(5) = -1
      IF (PROK) THEN
         WRITE(MPRINT,'(A,I3)')
     &     'Compute maximum matching (Maximum Transversal):',
     &        MTRANSLOC
         IF (MTRANSLOC.EQ.1)
     &   WRITE(MPRINT,'(A,I3)')' ... JOB =',MTRANSLOC
         IF (MTRANSLOC.EQ.2)
     &   WRITE(MPRINT,'(A,I3,A)')
     &     ' ... JOB =',MTRANSLOC,': BOTTLENECK THESIS'
         IF (MTRANSLOC.EQ.3)
     &   WRITE(MPRINT,'(A,I3,A)')
     &     ' ... JOB =',MTRANSLOC,': BOTTLENECK SIMAX'
         IF (MTRANSLOC.EQ.4)
     &   WRITE(MPRINT,'(A,I3,A)')
     &     ' ... JOB =',MTRANSLOC,': MAXIMIZE SUM DIAGONAL'
         IF (MTRANSLOC.EQ.5 .OR. MTRANSLOC.EQ.6)
     &   WRITE(MPRINT,'(A,I3,A)')
     &     ' ... JOB =',MTRANSLOC,
     &     ': MAXIMIZE PRODUCT DIAGONAL AND SCALE'
      ENDIF
      INFOG(23) = MTRANSLOC
      CNTL64(2) = huge(CNTL64(2))
      IRNW = 1
      IPIW = IRNW + NZTOT
      IF (MTRANSLOC.EQ.1) LIWMIN = 5_8*N8
      IF (MTRANSLOC.EQ.2) LIWMIN = 3_8*N8
      IF (MTRANSLOC.EQ.3) LIWMIN = 10_8*N8 + NZTOT
      IF (MTRANSLOC.EQ.4) LIWMIN = 2_8*N8
      IF (MTRANSLOC.EQ.5) LIWMIN = 5_8*N8
      IF (MTRANSLOC.EQ.6) LIWMIN = 5_8*N8 + NZTOT
      LIW = LIWMIN
      LIWG  = LIW + NZTOT
      ALLOCATE(IW(LIWG), stat=allocok)
      IF (allocok .GT. 0 ) THEN
         GOTO 410
      ENDIF
      ALLOCATE( IPQ8(N), IPE(N+1), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
         INFO( 1 ) = -7
         INFO( 2 ) = (2*N+1)*KEEP(10)
         GOTO 500               
      ENDIF
      IF (MTRANSLOC.EQ.1) THEN
       LDWMIN = N8+3_8
      ENDIF
      IF (MTRANSLOC.EQ.2) LDWMIN = max( N8+NZTOT , N8+3_8 )
      IF (MTRANSLOC.EQ.3) LDWMIN = max( NZTOT+1_8 , N8+3_8 )
      IF (MTRANSLOC.EQ.4) LDWMIN = 2_8 * N8 + 
     &                             max( NZTOT , N8+3_8 )
      IF (MTRANSLOC.EQ.5) LDWMIN = 3_8*N8 + NZTOT
      IF (MTRANSLOC.EQ.6) LDWMIN = 4_8*N8 + NZTOT
      LDW   = LDWMIN
      ALLOCATE(S2(LDW), stat=allocok)
      IF (allocok .GT. 0 ) THEN
         GOTO 430
      ENDIF
      IF(MTRANSLOC .NE. 1) LDW = LDW-NZTOT
      RSPOS = NZTOT
      CSPOS = RSPOS+N8
      NZREAL = 0_8
      DO 5 J=1,N
        IPQ8(J) = 0_8
  5   CONTINUE
      IF(K50 .EQ. 0) THEN
         DO 10 K=1,NZ
            I = idIRN(K)
            J = idJCN(K)
            IF ( (J.LE.N).AND.(J.GE.1).AND.
     &           (I.LE.N).AND.(I.GE.1) ) THEN
               IPQ8(J) = IPQ8(J) + 1_8
               NZREAL = NZREAL + 1_8
            ENDIF
 10      CONTINUE
      ELSE
         ZERODIAG = 0
         NZER_DIAG = N
         RZ_DIAG = 0
         DO K=1,NZ
            I = idIRN(K)
            J = idJCN(K)
            IF ( (J.LE.N).AND.(J.GE.1).AND.
     &           (I.LE.N).AND.(I.GE.1) ) THEN
               IPQ8(J) = IPQ8(J) + 1_8
               NZREAL = NZREAL + 1_8
               IF(I .NE. J) THEN
                  IPQ8(I) = IPQ8(I) + 1_8
                  NZREAL = NZREAL + 1_8
               ELSE
                  IF (ZERODIAG(I) .EQ. 0) THEN
                     ZERODIAG(I) = 1
                     IF(associated(idA)) THEN
                        ABSAK= abs(idA(K))
                        IF(ABSAK .EQ. dble(0.0D0)) THEN
                           RZ_DIAG = RZ_DIAG + 1
                        ENDIF
                         ZERODIAG(I) = exponent(ABSAK)
                         if ( ZERODIAG(I).EQ.0)  ZERODIAG(I)=1
                     ENDIF
                     NZER_DIAG = NZER_DIAG - 1                     
                  ELSE
                     IF(associated(idA)) THEN
                        ABSAK= abs(idA(K))
                        ZERODIAG(I) = ZERODIAG(I)+ exponent(ABSAK)
                        if ( ZERODIAG(I).EQ.0)  ZERODIAG(I)=1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         IF(MTRANSLOC .GE. 4) THEN
            DO I =1, N
               IF(ZERODIAG(I) .EQ. 0) THEN
                  IPQ8(I) = IPQ8(I) + 1_8
                  NZREAL = NZREAL + 1_8
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      IPE(1)   = 1
      DO 20 J=1,N
        IPE(J+1)   = IPE(J)+IPQ8(J)
  20  CONTINUE
      DO 25 J=1, N
        IPQ8(J ) = IPE(J)
  25  CONTINUE
      IF(K50 .EQ. 0) THEN
         IF (MTRANSLOC.EQ.1) THEN
            DO K=1,NZ
               I = idIRN(K)
               J = idJCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     &              (I.LE.N).AND.(I.GE.1)) THEN
                  KPOS              = IPQ8(J)
                  IW(IRNW+KPOS-1_8) = I
                  IPQ8(J)   = IPQ8(J) + 1_8
               ENDIF
            END DO
         ELSE
            IF ( .not.associated(idA)) THEN
               INFO(1) = -22
               INFO(2) = 4
               GOTO 500
            ENDIF
            DO K=1,NZ
               I = idIRN(K)
               J = idJCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     &              (I.LE.N).AND.(I.GE.1)) THEN
                  KPOS            = IPQ8(J)
                  IW(IRNW+KPOS-1) = I
                  S2(KPOS)         = abs(idA(K))
                  IPQ8(J)   = IPQ8(J) + 1_8
               ENDIF
            END DO
         ENDIF
      ELSE
         IF (MTRANSLOC.EQ.1) THEN
            DO K=1,NZ
               I = idIRN(K)
               J = idJCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     &              (I.LE.N).AND.(I.GE.1)) THEN
                  KPOS            = IPQ8(J)
                  IW(IRNW+KPOS-1) = I
                  IPQ8(J)   = IPQ8(J) + 1_8
                  IF(I.NE.J) THEN
                     KPOS            = IPQ8(I)
                     IW(IRNW+KPOS-1) = J
                     IPQ8(I)   = IPQ8(I) + 1_8
                  ENDIF
               ENDIF
            ENDDO
         ELSE
            IF ( .not.associated(idA) ) THEN
               INFO(1) = -22
               INFO(2) = 4
               GOTO 500
            ENDIF
            THEMAX = ZERO
            THEMIN = huge(THEMIN)
            DO K=1,NZ
               I = idIRN(K)
               J = idJCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     &              (I.LE.N).AND.(I.GE.1)) THEN
                  KPOS            = IPQ8(J)
                  IW(IRNW+KPOS-1_8) = I
                  S2(KPOS)          = abs(idA(K))
                  IPQ8(J)   = IPQ8(J) + 1_8
                  IF(abs(idA(K)) .GT. THEMAX) THEN
                     THEMAX = abs(idA(K))
                  ELSE IF(abs(idA(K)) .LT. THEMIN 
     &                    .AND. abs(idA(K)).GT. ZERO) THEN
                     THEMIN = abs(idA(K))
                  ENDIF
                  IF(I.NE.J) THEN
                     KPOS            = IPQ8(I)
                     IW(IRNW+KPOS-1) = J
                     S2(KPOS)         = abs(idA(K))
                     IPQ8(I)   = IPQ8(I) + 1_8
                  ENDIF
               ENDIF
            ENDDO
            DO I =1, N
               IF(ZERODIAG(I) .EQ. 0) THEN
                  KPOS            = IPQ8(I)
                  IW(IRNW+KPOS-1) = I
                  S2(KPOS)         = ZERO
                  IPQ8(I)   = IPQ8(I) + 1_8
               ENDIF
            ENDDO
            IF ( THEMAX .NE. ZERO ) THEN
              CNTL64(2) = (log(THEMAX/THEMIN))*(dble(N))
     &             - log(THEMIN) + ONE
            ENDIF
         ENDIF
      ENDIF
      DUPPLI = .FALSE.
      NZsave = NZREAL
      FLAG => IKEEPALLOC(2*N+1:3*N)
      IF(MTRANSLOC.NE.1) THEN
         CALL DMUMPS_SUPPRESS_DUPPLI_VAL(N,NZREAL,IPE(1),IW(IRNW),S2,
     &        PERM(1),IPQ8(1))
      ELSE
         CALL DMUMPS_SUPPRESS_DUPPLI_STR(N,NZREAL,IPE(1),IW(IRNW),
     &        PERM(1))
      ENDIF
      IF(NZREAL .NE. NZsave) DUPPLI = .TRUE.
      LS2 = NZTOT
      IF ( MTRANSLOC .EQ. 1 ) THEN
         LS2 = 1_8
         LDW = 1_8
      ENDIF
      CALL DMUMPS_MTRANS_DRIVER(MTRANSLOC ,N, N, NZREAL, 
     &     IPE, IW(IRNW), S2(1), LS2,
     &     NUMNZ, PERM(1), LIW, IW(IPIW), LDW, S2(LS2+1),
     &     IPQ8,
     &     ICNTL64, CNTL64, INFO64, INFO)
      IF (INFO(1).LT.0) THEN 
         IF (LP.GT.0 .AND. ICNTL(4).GE.1)
     &        WRITE(LP,'(A,I5)')
     &   ' Not enough memory in MAXTRANS INFO(1)=',INFO(1)
         GOTO 500 
      ENDIF
      IF (INFO64(1).LT.0) THEN
         IF (LP.GT.0 .AND. ICNTL(4).GE.1)
     &        WRITE(LP,'(A,I5)')
     &   ' INTERNAL ERROR in MAXTRANS INFO(1)=',INFO64(1)
         INFO(1) = -9964
         INFO(2) = INFO64(1)
         GO TO 500
      ENDIF
      IF (INFO64(1).GT.0) THEN
         IF (MP.GT.0 .AND. ICNTL(4).GE.2)
     &        WRITE(MP,'(A,I5)')
     &        ' WARNING in MAXTRANS INFO(1)=',INFO64(1)
      ENDIF
      KER_SIZE = 0
      IF(K50 .EQ. 2) THEN
         DO I=1,N
            IF(ZERODIAG(I) .EQ. 0) THEN
               IF(PERM(I) .EQ. I) THEN
                  KER_SIZE = KER_SIZE + 1
                  PERM(I) = -I
                  STR_KER(KER_SIZE) = I
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      IF (NUMNZ.LT.N) GO TO 400
      IF(K50 .EQ. 0) THEN
         IDENT = .TRUE.
         IF (MTRANS .EQ. 0 ) GOTO 102
         DO 80 J=1,N
            JPERM = PERM(J)
            IW(IRNW+int(JPERM-1,8)) = J
            IF (JPERM.NE.J) IDENT = .FALSE.
 80      CONTINUE
         IF(IDENT) THEN 
            MTRANS = 0
         ELSE
            IF(MTRANS .EQ. 7) THEN
               MTRANS = -9876543
               GOTO 102
            ENDIF
            IF (PROK) WRITE(MPRINT,'(A)')
     &           ' ... Apply column permutation'
            DO 100 K=1,NZ
               J = idJCN(K)
               IF ((J.LE.0).OR.(J.GT.N)) GO TO 100
               idJCN(K) = IW(IRNW+int(J-1,8))
 100        CONTINUE
            IF (MP.GT.0 .AND. ICNTL(4).GE.2)
     &           WRITE(MP,'(/A)')
     &           ' WARNING input matrix data modified'
         ENDIF
 102     CONTINUE
         IF (SCALINGLOC) THEN
            IF ( associated(idCOLSCA))
     &           DEALLOCATE( idCOLSCA )
            IF ( associated(idROWSCA))
     &           DEALLOCATE( idROWSCA )
            ALLOCATE( idCOLSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               INFO(1)=-5
               INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in DMUMPS_ANA_O'
                  WRITE (LP,'(A)') 
     &                 '** Failure during allocation of COLSCA'
                  GOTO 500
               ENDIF
            ENDIF
            ALLOCATE( idROWSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               INFO(1)=-5
               INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in DMUMPS_ANA_O'
                  WRITE (LP,'(A)') 
     &                 '** Failure during allocation of ROWSCA'
                  GOTO 500
               ENDIF
            ENDIF
            KEEP(52) = -2
            KEEP(74) = 1
            MAXDBL = log(huge(MAXDBL))
            DO J=1,N
               IF(S2(RSPOS+J) .GT. MAXDBL) THEN
                  S2(RSPOS+J) = ZERO
               ENDIF
               IF(S2(CSPOS+J) .GT. MAXDBL) THEN
                  S2(CSPOS+J)= ZERO
               ENDIF
            ENDDO
            DO 105 J=1,N
               J8 = int(J,8)
               idROWSCA(J) = exp(S2(RSPOS+J8))
               IF(idROWSCA(J) .EQ. ZERO) THEN
                  idROWSCA(J) = ONE
               ENDIF
               IF ( MTRANS .EQ.  -9876543 .OR. MTRANS.EQ. 0 ) THEN
                 idCOLSCA(J)= exp(S2(CSPOS+J8))
                 IF(idCOLSCA(J) .EQ. ZERO) THEN
                   idCOLSCA(J) = ONE
                 ENDIF
               ELSE
                 idCOLSCA(IW(IRNW+J8-1_8))= exp(S2(CSPOS+J8))
                 IF(idCOLSCA(IW(IRNW+J8-1_8)) .EQ. ZERO) THEN
                   idCOLSCA(IW(IRNW+J8-1_8)) = ONE
                 ENDIF
               ENDIF
 105        CONTINUE
         ENDIF
      ELSE
         IDENT = .FALSE.         
         IF(SCALINGLOC) THEN
            IF ( associated(idCOLSCA)) DEALLOCATE( idCOLSCA )
            IF ( associated(idROWSCA)) DEALLOCATE( idROWSCA )
            ALLOCATE( idCOLSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               INFO(1)=-5
               INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in DMUMPS_ANA_O'
                  WRITE (LP,'(A)') 
     &                 '** Failure during allocation of COLSCA'
                  GOTO 500
               ENDIF
            ENDIF
            ALLOCATE( idROWSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               INFO(1)=-5
               INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in DMUMPS_ANA_O'
                  WRITE (LP,'(A)') 
     &                 '** Failure during allocation of ROWSCA'
                  GOTO 500
               ENDIF
            ENDIF
            KEEP(52) = -2
            KEEP(74) = 1
            MAXDBL = log(huge(MAXDBL))
            DO J=1,N
               J8 = int(J,8)
               IF(S2(RSPOS+J8)+S2(CSPOS+J8) .GT. MAXDBL) THEN
                  S2(RSPOS+J8) = ZERO
                  S2(CSPOS+J8)= ZERO
               ENDIF
            ENDDO
            DO J=1,N
               J8 = int(J,8)
               IF(PERM(J) .GT. 0) THEN
                  idROWSCA(J) = 
     &                 exp((S2(RSPOS+J8)+S2(CSPOS+J8))/TWO)
                  IF(idROWSCA(J) .EQ. ZERO) THEN
                     idROWSCA(J) = ONE
                  ENDIF
                  idCOLSCA(J)= idROWSCA(J)
               ENDIF
            ENDDO
            DO JPOS=1,KER_SIZE
               I = STR_KER(JPOS)
               COLNORM = ZERO
               DO K = IPE(I),IPE(I+1) - 1
                  IF ( PERM( IW( IRNW+K-1_8) ) > 0 ) THEN
                    COLNORM = max(COLNORM,S2(J))
                  ENDIF
               ENDDO
               COLNORM = exp(COLNORM) 
               idROWSCA(I) = ONE / COLNORM
               idCOLSCA(I) = idROWSCA(I)
            ENDDO
         ENDIF
         IF(MTRANS .EQ. 7 .OR. KEEP(95) .EQ. 0) THEN
            IF( (NZER_DIAG+RZ_DIAG) .LT. (N/10) 
     &           .AND. KEEP(95) .EQ. 0) THEN
               MTRANS = 0
               KEEP(95) = 1
               GOTO 390
            ELSE
               IF(KEEP(95) .EQ. 0) THEN
                 IF(SCALINGLOC) THEN
                  KEEP(95) = 3
                 ELSE
                  KEEP(95) = 2   
                 ENDIF
               ENDIF
               IF(MTRANS .EQ. 7) MTRANS = 5
            ENDIF
         ENDIF
         IF(MTRANS .EQ. 0) GOTO 390
         ICNTL_SYM_MWM = 0
         INFO_SYM_MWM = 0
         IF(MTRANS .EQ. 5 .OR. MTRANS .EQ. 6 .OR.
     &        MTRANS .EQ. 7) THEN
            ICNTL_SYM_MWM(1) = 0
            ICNTL_SYM_MWM(2) = 1
         ELSE IF(MTRANS .EQ. 4) THEN
            ICNTL_SYM_MWM(1) = 2
            ICNTL_SYM_MWM(2) = 1
         ELSE
            ICNTL_SYM_MWM(1) = 0
            ICNTL_SYM_MWM(2) = 1
         ENDIF
         MARKED => IKEEPALLOC(N+1:2*N)
         FLAG => IKEEPALLOC(2*N+1:3*N)
         PIV_OUT => WORK2(1:N)
         IF(MTRANSLOC .LT. 4) THEN
            LSC = 1
         ELSE
            LSC = 2*N
         ENDIF
         CALL DMUMPS_SYM_MWM(
     &        N, NZREAL, IPE, IW(IRNW), S2(1),LSC, PERM(1),
     &        ZERODIAG(1),
     &        ICNTL_SYM_MWM, S2(LSC+1),MARKED(1),FLAG(1),
     &        PIV_OUT(1), INFO_SYM_MWM)
         IF(INFO_SYM_MWM(1) .NE. 0) THEN
            WRITE(*,*) '** Error in DMUMPS_ANA_O'
            RETURN
         ENDIF
         IF(INFO_SYM_MWM(3) .EQ. N) THEN
            IDENT = .TRUE.
         ELSEIF  ( (ICNTL(12).EQ.0).AND.
     &           ( (N-INFO_SYM_MWM(4)-INFO_SYM_MWM(3)) .GT. N/10 )
     &           ) THEN
            IDENT = .TRUE.
            KEEP(95) = 1
         ELSE
            DO I=1,N
               PERM(I) = PIV_OUT(I)
            ENDDO
         ENDIF
         KEEP(93) = INFO_SYM_MWM(4)
         KEEP(94) = INFO_SYM_MWM(3)
         IF (IDENT) MTRANS=0
      ENDIF
 390  IF(MTRANS .EQ. 0) THEN
         KEEP(95) = 1 
         IF (PROK) THEN
           WRITE (MPRINT,'(A)')
     &  ' ... Column permutation not used'
         ENDIF
      ENDIF
      GO TO 500
 400  IF ((LP.GE.0).AND.(ICNTL(4).GE.1))
     &   WRITE (LP,'(/A)') '** Error: Matrix is structurally singular'
      INFO(1) = -6
      INFO(2) = NUMNZ
      GOTO 500
 410  IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
       WRITE (LP,'(/A)') '** Error in DMUMPS_ANA_O'
       WRITE (LP,'(A,I14)')
     & '** Failure during allocation of INTEGER array of size ',
     & LIWG
      ENDIF
      INFO(1) = -7
      CALL MUMPS_SET_IERROR(LIWG,INFO(2))
      GOTO 500
 430  IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
       WRITE (LP,'(/A)') '** Error in DMUMPS_ANA_O'
       WRITE (LP,'(A)') '** Failure during allocation of S2'
      ENDIF
      INFO(1) = -5
      CALL MUMPS_SET_IERROR(LDW,INFO(2))
 500  CONTINUE
      IF (allocated(IW)) DEALLOCATE(IW)
      IF (allocated(S2)) DEALLOCATE(S2)
      IF (allocated(IPE)) DEALLOCATE(IPE)
      IF (allocated(IPQ8)) DEALLOCATE(IPQ8)
      RETURN
      END SUBROUTINE DMUMPS_ANA_O
      END MODULE DMUMPS_ANA_AUX_M
      SUBROUTINE DMUMPS_ANA_K(N,IPE, IW, LW, IWFR, IPS, IPV, 
     &                           NV, FLAG,
     &                           NCMPA, SIZE_SCHUR, PARENT)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, SIZE_SCHUR
      INTEGER, INTENT(IN)    :: IPS(N)
      INTEGER(8), INTENT(IN) :: LW
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: IPV(N), NV(N), PARENT(N)
      INTEGER(8), INTENT(INOUT) :: IWFR
      INTEGER(8), INTENT(INOUT) :: IPE(N)
      INTEGER, INTENT(INOUT)    :: IW(LW)
      INTEGER, INTENT(OUT)      ::  FLAG(N)
      INTEGER I,J,ML,MS,ME,MINJS,IE,KDUMMY
      INTEGER LN,JS,JE
      INTEGER(8) :: JP, JP1, JP2, LWFR, IP
      DO 10 I=1,N
        FLAG(I) = 0
        NV(I) = 0
        J = IPS(I)
        IPV(J) = I
   10 CONTINUE
      NCMPA = 0
      DO 100 ML=1,N-SIZE_SCHUR 
        MS = IPV(ML)
        ME = MS
        FLAG(MS) = ME
        IP = IWFR
        MINJS = N
        IE = ME
        DO 70 KDUMMY=1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0_8) GO TO 60
          LN = IW(JP)
          DO 50 JP1=1_8,int(LN,8)
            JP = JP + 1_8
            JS = IW(JP)
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
            IPE(IE) = JP
            IW(JP) = LN - int(JP1)
            CALL DMUMPS_ANA_D(N, IPE, IW, IP-1_8, LWFR, NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP=IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1_8
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
   40       IW(IWFR) = JS
            MINJS = min0(MINJS,IPS(JS)+0)
            IWFR = IWFR + 1_8
   50     CONTINUE
   60     IPE(IE) = int(-ME,8)
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
        IPE(ME) = 0_8
        NV(ME) = 1
        GO TO 100
   90   MINJS = IPV(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
        IW(IWFR) = IW(IP)
        IW(IP) = int(IWFR - IP)
        IPE(ME) = IP
        IWFR = IWFR + 1_8
  100 CONTINUE
      IF (SIZE_SCHUR == 0) GOTO 500
      DO ML = N-SIZE_SCHUR+1,N
        ME = IPV(ML)
        IE = ME
        DO KDUMMY=1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0_8) GO TO 160
          LN = IW(JP)
  160     IPE(IE) = int(-IPV(N-SIZE_SCHUR+1),8)
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 190
        ENDDO
  190   NV(ME) = 0
        IPE(ME) = int(-IPV(N-SIZE_SCHUR+1),8)
      ENDDO
      ME = IPV(N-SIZE_SCHUR+1)
      IPE(ME) = 0_8
      NV(ME) = SIZE_SCHUR
  500 DO I=1,N
       PARENT(I) = int(IPE(I))
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ANA_K
      SUBROUTINE DMUMPS_ANA_J(N, NZ, IRN, ICN, PERM,
     & IW, LW, IPE, IQ, FLAG,
     & IWFR, IFLAG, IERROR, MP)
      INTEGER, INTENT(IN)    :: N
      INTEGER(8), INTENT(IN) :: NZ, LW
      INTEGER, INTENT(IN) :: IRN(NZ), ICN(NZ) 
      INTEGER, INTENT(IN) :: PERM(N)
      INTEGER, INTENT(IN)  ::  MP
      INTEGER(8), INTENT(OUT):: IWFR
      INTEGER, INTENT(OUT)   :: IERROR
      INTEGER, INTENT(OUT)   :: IQ(N)
      INTEGER(8), INTENT(OUT) :: IPE(N) 
      INTEGER, INTENT(OUT) :: IW(LW) 
      INTEGER, INTENT(OUT) :: FLAG(N)
      INTEGER, INTENT(INOUT) :: IFLAG
      INTEGER  :: I,J,LBIG,IN,LEN,JDUMMY,L1
      INTEGER(8) :: K, K1, K2, KL, KID
      IERROR = 0
      DO 10 I=1,N
        IQ(I) = 0
   10 CONTINUE
      DO 80 K=1_8,NZ
        I = IRN(K)
        J = ICN(K)
        IW(K) = -I
        IF (I.EQ.J) GOTO 40
        IF (I.GT.J) GOTO 30
        IF (I.GE.1 .AND. J.LE.N) GO TO 60
        GO TO 50
   30   IF (J.GE.1 .AND. I.LE.N) GO TO 60
        GO TO 50
   40   IW(K) = 0
        IF (I.GE.1 .AND. I.LE.N) GO TO 80
   50   IERROR = IERROR + 1
        IW(K) = 0
        IF (IERROR.LE.1 .AND. MP.GT.0) WRITE (MP,99999) 
        IF (IERROR.LE.10 .AND. MP.GT.0) WRITE (MP,99998) K, I, J
        GO TO 80
   60   IF (PERM(J).GT.PERM(I)) GO TO 70
        IQ(J) = IQ(J) + 1
        GO TO 80
   70   IQ(I) = IQ(I) + 1
   80 CONTINUE
      IF (IERROR.GE.1) THEN
        IF (mod(IFLAG,2) .EQ. 0) IFLAG = IFLAG+1
      ENDIF
      IWFR = 1_8
      LBIG = 0
      DO 100 I=1,N
        L1 = IQ(I)
        LBIG = max0(L1,LBIG)
        IWFR = IWFR + int(L1,8)
        IPE(I) = IWFR - 1_8
  100 CONTINUE
      DO 140 K=1_8,NZ
        I = -IW(K)
        IF (I.LE.0) GO TO 140
        KL = K
        IW(K) = 0
        DO 130 KID=1,NZ
          J = ICN(KL)
          IF (PERM(I).LT.PERM(J)) GO TO 110
          KL = IPE(J)
          IPE(J) = KL - 1_8
          IN = IW(KL)
          IW(KL) = I
          GO TO 120
  110     KL = IPE(I)
          IPE(I) = KL - 1_8
          IN = IW(KL)
          IW(KL) = J
  120     I = -IN
          IF (I.LE.0) GO TO 140
  130   CONTINUE
  140 CONTINUE
      K = IWFR - 1_8
      KL = K + int(N,8)
      IWFR = KL + 1_8
      DO 170 I=1,N
        FLAG(I) = 0
        J = N + 1 - I
        LEN = IQ(J)
        IF (LEN.LE.0) GO TO 160
        DO 150 JDUMMY=1,LEN
          IW(KL) = IW(K)
          K = K - 1_8
          KL = KL - 1_8
  150   CONTINUE
  160   IPE(J) = KL
        KL = KL - 1_8
  170 CONTINUE
      IF (LBIG.GE.huge(N)) GO TO 190
      DO 180 I=1,N
        K = IPE(I)
        IW(K) = IQ(I)
        IF (IQ(I).EQ.0) IPE(I) = 0_8
  180 CONTINUE
      GO TO 230
  190 IWFR = 1_8
      DO 220 I=1,N
        K1 = IPE(I) + 1_8
        K2 = IPE(I) + int(IQ(I),8)
        IF (K1.LE.K2) GO TO 200
        IPE(I) = 0_8
        GO TO 220
  200   IPE(I) = IWFR
        IWFR = IWFR + 1_8
        DO 210 K=K1,K2
          J = IW(K)
          IF (FLAG(J).EQ.I) GO TO 210
          IW(IWFR) = J
          IWFR = IWFR + 1_8
          FLAG(J) = I
  210   CONTINUE
        K = IPE(I)
        IW(K) = int(IWFR - K - 1_8)
  220 CONTINUE
  230 RETURN
99999 FORMAT (' *** WARNING MESSAGE FROM DMUMPS_ANA_J ***' )
99998 FORMAT (I6, ' NON-ZERO (IN ROW, I6, 11H AND COLUMN ', I6,
     & ') IGNORED')
      END SUBROUTINE DMUMPS_ANA_J
      SUBROUTINE DMUMPS_ANA_D(N, IPE, IW, LW, IWFR,NCMPA)
      INTEGER, INTENT(IN)      :: N
      INTEGER(8), INTENT(IN)   :: LW
      INTEGER(8), INTENT(OUT)  :: IWFR
      INTEGER(8), INTENT(INOUT):: IPE(N)
      INTEGER, INTENT(INOUT)   :: NCMPA
      INTEGER, INTENT(INOUT)   :: IW(LW)
      INTEGER    :: I, IR
      INTEGER(8) :: K1, K, K2, LWFR
      NCMPA = NCMPA + 1
      DO 10 I=1,N
        K1 = IPE(I)
        IF (K1.LE.0_8) GO TO 10
        IPE(I) = int(IW(K1), 8)
        IW(K1) = -I
   10 CONTINUE
      IWFR = 1_8
      LWFR = IWFR
      DO 60 IR=1,N
        IF (LWFR.GT.LW) GO TO 70
        DO 20 K=LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
   30   I = -IW(K)
        IW(IWFR) = int(IPE(I))
        IPE(I) = int(IWFR,8)
        K1 = K + 1_8
        K2 = K + int(IW(IWFR),8)
        IWFR = IWFR + 1_8
        IF (K1.GT.K2) GO TO 50
        DO 40 K=K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1_8
   40   CONTINUE
   50   LWFR = K2 + 1_8
   60 CONTINUE
   70 RETURN
      END SUBROUTINE DMUMPS_ANA_D
#if defined(OLDDFS)
      SUBROUTINE DMUMPS_ANA_L(N, IPE, NV, IPS, NE, NA, NFSIZ, 
     &                  NSTEPS,
     &                  FILS, FRERE,NDD,NEMIN, KEEP60)
      INTEGER N,NSTEPS
      INTEGER NDD(N)
      INTEGER FILS(N), FRERE(N)
      INTEGER IPS(N), NE(N), NA(N), NFSIZ(N)
      INTEGER IPE(N), NV(N)
      INTEGER NEMIN, KEEP60
      INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
      INTEGER K,L,ISON,IN,INP,IFSON,INC,INO
      INTEGER INOS,IB,IL
      DO 10 I=1,N
        IPS(I) = 0
        NE(I) = 0
   10 CONTINUE
      DO 20 I=1,N
        IF (NV(I).GT.0) GO TO 20
        IF = -IPE(I)
        IS = -IPS(IF)
        IF (IS.GT.0) IPE(I) = IS
        IPS(IF) = -I
   20 CONTINUE
      NR = N + 1
      DO 50 I=1,N
        IF (NV(I).LE.0) GO TO 50
        IF = -IPE(I)
        IF (IF.NE.0) THEN
         IS = -IPS(IF)
         IF (IS.GT.0) IPE(I) = IS
         IPS(IF) = -I
        ELSE
         NR = NR - 1
         NE(NR) = I
        ENDIF
   50 CONTINUE
      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE
      NR1 = NR
      INS = 0
 1000 IF (NR1.GT.N) GO TO 1151
      INS = NE(NR1)
      NR1 = NR1 + 1
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
 1080 IF (IPE(INS).LT.0) THEN
       INS       = -IPE(INS)
       FILS(INS) = 0
       GO TO 1080
      ENDIF
      IF (IPE(INS).EQ.0) THEN
       INS = 0
       GO TO 1000
      ENDIF
      INB = IPE(INS)
      IF (NV(INB).EQ.0) THEN
       INS = INB
       GO TO 1070
      ENDIF
      IF (NV(INB).GE.NV(INS)) THEN
       INS = INB
       GO TO 1070
      ENDIF
      INF = INB
 1090 INF = IPE(INF)
      IF (INF.GT.0) GO TO 1090
      INF  = -INF
      INFS = -FILS(INF)
      IF (INFS.EQ.INS) THEN
       FILS(INF) = -INB
       IPS(INF)  = -INB
       IPE(INS)  = IPE(INB)
       IPE(INB)  = INS
       INS       = INB
       GO TO 1070
      ENDIF
      INSW = INFS
 1100 INFS = IPE(INSW)
      IF (INFS.NE.INS) THEN
       INSW = INFS
       GO TO 1100
      ENDIF
      IPE(INS) = IPE(INB)
      IPE(INB) = INS
      IPE(INSW)= INB
      INS      =INB
      GO TO 1070
 1151 CONTINUE
      DO 51 I=1,N
       FRERE(I) = IPE(I)
       FILS(I) = IPS(I)
 51   CONTINUE
      IS = 1
      I  = 0
      IL = 0
      DO 160 K=1,N
        IF (I.GT.0) GO TO 60
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
   60   DO 70 L=1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
   80   IPS(I) = K
        NE(IS) = NE(IS) + 1
        IF (NV(I).GT.0) GO TO 89
      IN = I
 81   IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 81
      IF = -IN
      IN = IF
 82   INL = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 82
      IFSON = -IN
      FILS(INL) = I
      IN = I
 83   INP = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 83
      IF (IFSON .EQ. I) GO TO 86
      FILS(INP) = -IFSON
      IN = IFSON
 84   INC =IN
      IN = FRERE(IN)
      IF (IN.NE.I) GO TO 84
      FRERE(INC) = FRERE(I)
      GO TO 120
 86   IF (FRERE(I).LT.0) FILS(INP) = 0
      IF (FRERE(I).GT.0) FILS(INP) = -FRERE(I)
      GO TO 120
   89   IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        NDD(IS) = NV(I)
        NFSIZ(I) = NV(I)
        IF (NA(IS).LT.1) GO TO 110
        IF (   (KEEP60.NE.0).AND.
     &         (NE(IS).EQ.NDD(IS)) ) GOTO 110
        IF (NDD(IS-1)-NE(IS-1).EQ.NDD(IS)) GO TO 100
        IF ((NE(IS-1).GE.NEMIN).AND.
     &         (NE(IS).GE.NEMIN) ) GO TO 110
        IF (2*NE(IS-1)*(NDD(IS)-NDD(IS-1)+NE(IS-1)).GE.
     &    ((NDD(IS)+NE(IS-1))*
     &    (NDD(IS)+NE(IS-1))*NEMIN/100)) GO TO 110
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        NDD(IS-1) = NDD(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
      IN=I
 101  INL = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 101
      IFSON = -IN
      IN = IFSON
 102  INO = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 102
      FILS(INL) = INO
      NFSIZ(I) = NDD(IS-1)
      IN = INO
 103  INP = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 103
      INOS = -IN
      IF (IFSON.EQ.INO) GO TO 107
      IN = IFSON
      FILS(INP) = -IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
      IF (INOS.EQ.0) FRERE(INS) = -I
      IF (INOS.NE.0) FRERE(INS) =  INOS
      IF (INOS.EQ.0) GO TO 109
 107  IN = INOS
      IF (IN.EQ.0) GO TO 109
 108  INT = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 108
      FRERE(INT) = -I
 109  CONTINUE
        GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.LT.0) GOTO 150
        IF (IB.EQ.0) GOTO 140
        NA(IL) = 0
  140   I = IB
        GO TO 160
  150   I = -IB
        IL = IL + 1
  160 CONTINUE
      NSTEPS = IS - 1
      DO 170 I=1,N
        K = FILS(I)
        IF (K.GT.0) THEN
          FRERE(K)  = N + 1
          NFSIZ(K)  = 0
        ENDIF
 170  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_ANA_L
#else
      SUBROUTINE DMUMPS_ANA_LNEW(N, IPE, NV, IPS, NE, NA, NFSIZ,
     &               NODE, NSTEPS,
     &               FILS, FRERE, ND, NEMIN, SUBORD, KEEP60, 
     &               KEEP20, KEEP38, NAMALG,NAMALGMAX,
     &               CUMUL,KEEP50, ICNTL13, KEEP37, NSLAVES,
     &               ALLOW_AMALG_TINY_NODES 
     &               , BLKON, SIZEOFBLOCKS, LSIZEOFBLOCKS
     &               )
      IMPLICIT NONE
      INTEGER  N, NSTEPS, KEEP60, KEEP20, KEEP38, KEEP50
      INTEGER ND(N), NFSIZ(N)
      INTEGER IPE(N), FILS(N), FRERE(N), SUBORD(N)
      INTEGER NV(N), IPS(N), NE(N), NA(N), NODE(N)
      INTEGER NEMIN,AMALG_COUNT
      INTEGER NAMALG(N),NAMALGMAX, CUMUL(N)
      DOUBLE PRECISION SIZE_DADI_AMALGAMATED, PERCENT_FILL
      DOUBLE PRECISION ACCU, FLOPS_FATHER, FLOPS_SON,
     &                  FLOPS_AVANT, FLOPS_APRES
      INTEGER ICNTL13, KEEP37, NSLAVES
      LOGICAL ALLOW_AMALG_TINY_NODES
      LOGICAL, INTENT(IN) :: BLKON
      INTEGER, INTENT(IN) :: LSIZEOFBLOCKS
      INTEGER, INTENT(IN) :: SIZEOFBLOCKS(LSIZEOFBLOCKS) 
#if  defined(NOAMALGTOFATHER)
#else
#endif
      INTEGER I,IF,IS,NR,INS
      INTEGER K,L,ISON,IN,IFSON,INO
      INTEGER INOS,IB,IL
      INTEGER IPERM
      INTEGER MAXNODE
#if defined(NOAMALGTOFATHER)
      INTEGER INB,INF,INFS,INL,INSW,INT1,NR1
#else
      INTEGER DADI
#endif
      LOGICAL AMALG_TO_father_OK
      AMALG_COUNT = 0
      DO 10 I=1,N
        CUMUL(I)= 0
        IPS(I)  = 0
        NE(I)   = 0
        SUBORD(I) = 0
        NAMALG(I) = 0
   10 CONTINUE
      DO I=1,N
        IF (BLKON) THEN
         NODE(I) = SIZEOFBLOCKS(I)
        ELSE
         NODE(I) = 1
        ENDIF
      ENDDO
      FRERE(1:N) = IPE(1:N)
      NR = N + 1
      MAXNODE = 1   
      DO 50 I=1,N
        IF = -FRERE(I)
        IF (NV(I).EQ.0) THEN
          IF (SUBORD(IF).NE.0) SUBORD(I) = SUBORD(IF)
          SUBORD(IF) = I
          IF (BLKON) THEN
            NODE(IF) = NODE(IF)+SIZEOFBLOCKS(I)
          ELSE
            NODE(IF) = NODE(IF)+1
          ENDIF
          MAXNODE = max(NODE(IF),MAXNODE)
        ELSE
          IF (IF.NE.0) THEN
            IS = -IPS(IF)
            IF (IS.GT.0) FRERE(I) = IS
            IPS(IF) = -I
          ELSE
            NR = NR - 1
            NE(NR) = I
          ENDIF
        ENDIF
   50 CONTINUE
        MAXNODE = int(dble(MAXNODE)*dble(NEMIN) / dble(100))
        MAXNODE = max(MAXNODE,2000)
#if defined(NOAMALGTOFATHER)
      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE
      NR1 = NR
      INS = 0
 1000 IF (NR1.GT.N) GO TO 1151
      INS = NE(NR1)
      NR1 = NR1 + 1
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
 1080 IF (FRERE(INS).LT.0) THEN
       INS       = -FRERE(INS)
       FILS(INS) = 0
       GO TO 1080
      ENDIF
      IF (FRERE(INS).EQ.0) THEN
       INS = 0
       GO TO 1000
      ENDIF
      INB = FRERE(INS)
      IF (NV(INB).GE.NV(INS)) THEN
       INS = INB
       GO TO 1070
      ENDIF
      INF = INB
 1090 INF = FRERE(INF)
      IF (INF.GT.0) GO TO 1090
      INF  = -INF
      INFS = -FILS(INF)
      IF (INFS.EQ.INS) THEN
        FILS(INF) = -INB
        IPS(INF)  = -INB
        FRERE(INS)  = FRERE(INB)
        FRERE(INB)  = INS
      ELSE
        INSW = INFS
 1100   INFS = FRERE(INSW)
        IF (INFS.NE.INS) THEN
          INSW = INFS
          GO TO 1100
        ENDIF
        FRERE(INS) = FRERE(INB)
        FRERE(INB) = INS
        FRERE(INSW)= INB
      ENDIF
        INS      = INB
        GO TO 1070
 1151   CONTINUE
#endif
      DO 51 I=1,N
       FILS(I) = IPS(I)
 51   CONTINUE
      IS = 1
      I = 0
      IPERM = 1
      DO 160 K=1,N
        AMALG_TO_father_OK=.FALSE.
        IF (I.LE.0) THEN
         IF (NR.GT.N) EXIT
         I = NE(NR)
         NE(NR) = 0
         NR = NR + 1
         IL = N
         NA(N) = 0
        ENDIF
        DO 70 L=1,N
          IF (IPS(I).GE.0) EXIT
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
#if ! defined(NOAMALGTOFATHER)
        DADI = -IPE(I)  
        IF ( (DADI.NE.0) .AND.
     &      (
     &       (KEEP60.EQ.0).OR.
     &       ( (KEEP20.NE.DADI).AND.(KEEP38.NE.DADI) )
     &      )
     &     ) THEN
           ACCU = dble(2)*dble(NODE(I))*dble(NV(DADI)-NV(I)+NODE(I))
           SIZE_DADI_AMALGAMATED = 
     &           dble(NV(DADI)+NODE(I)) *
     &           dble(NV(DADI)+NODE(I)) 
           PERCENT_FILL = dble(100) * ACCU / SIZE_DADI_AMALGAMATED
           ACCU = ACCU + dble(CUMUL(I))
           AMALG_TO_father_OK =  ( 
     &           ( (NODE(I).LE.MAXNODE).AND.(NODE(DADI).LE.MAXNODE) ) 
     &         .OR. 
     &           ( (NODE(I).LE.NEMIN.and. NODE(DADI).GT. MAXNODE)
     &     .OR.(NODE(DADI).LE.NEMIN .and. NODE(I).GT.MAXNODE)))
           AMALG_TO_father_OK = ( AMALG_TO_father_OK .AND.
     &       ( PERCENT_FILL < dble(NEMIN) ) )
           AMALG_TO_father_OK = ( AMALG_TO_father_OK .AND.
     &     ( ACCU / SIZE_DADI_AMALGAMATED .LE. dble(NEMIN)) )
           IF (AMALG_TO_father_OK) THEN
              CALL MUMPS_GET_FLOPS_COST(NV(I),NODE(I),NODE(I),
     &                                  KEEP50,1,FLOPS_SON)
              CALL MUMPS_GET_FLOPS_COST(NV(DADI),NODE(DADI),
     &                             NODE(DADI),
     &                             KEEP50,1,FLOPS_FATHER)
              FLOPS_AVANT = FLOPS_FATHER+FLOPS_SON
     &                      + max(dble(200.0) * dble(NV(I)-NODE(I))
     &                            * dble(NV(I)-NODE(I)),
     &                            dble(10000.0))
              CALL MUMPS_GET_FLOPS_COST(NV(DADI)+NODE(I),
     &                             NODE(DADI)+NODE(I),
     &                             NODE(DADI)+NODE(I),
     &                             KEEP50,1,FLOPS_APRES)
              IF (FLOPS_APRES.GT.FLOPS_AVANT*
     &         (dble(1)+dble(max(8,NEMIN)-8)/dble(100))) THEN
                 AMALG_TO_father_OK = .FALSE.
              ENDIF
           ENDIF
           IF ( (NV(I).GT. 50*NV(DADI)).AND. (NSLAVES.GT.1) 
     &          .AND. (ICNTL13.LE.0)
     &          .AND. (NV(I).GT. KEEP37) )  THEN
             IF ( ( ACCU / SIZE_DADI_AMALGAMATED ) .LT. 0.2 ) THEN
               AMALG_TO_father_OK = .TRUE.
             ENDIF
           ENDIF
           IF ( ALLOW_AMALG_TINY_NODES .AND.
     &     NODE(I) * 900 .LE. NV(DADI) - NAMALG(DADI)) THEN
             IF ( NAMALG(DADI) < (NV(DADI)-NAMALG(DADI))/50 ) THEN
                AMALG_TO_father_OK = .TRUE.
                NAMALG(DADI) = NAMALG(DADI) + NODE(I)
             ENDIF
           ENDIF
           IF ( DADI .EQ. -FRERE(I) 
     &       .AND. -FILS(DADI).EQ.I  
     &       ) THEN
             AMALG_TO_father_OK = ( AMALG_TO_father_OK .OR.
     &                          ( NV(I)-NODE(I).EQ.NV(DADI)) )
           ENDIF
           IF (AMALG_TO_father_OK) THEN
             CUMUL(DADI)=CUMUL(DADI)+nint(ACCU)
             NAMALG(DADI) = NAMALG(DADI) + NAMALG(I)
             AMALG_COUNT = AMALG_COUNT+1
             IN = DADI
 75          IF (SUBORD(IN).EQ.0) GOTO 76
               IN = SUBORD(IN)
               GOTO 75
 76          CONTINUE
             SUBORD(IN) = I
             NV(I)      = 0
             IFSON = -FILS(DADI)
             IF (IFSON.EQ.I) THEN
              IF (FILS(I).LT.0) THEN
                FILS(DADI) =  FILS(I)
                GOTO 78
              ELSE
                IF (FRERE(I).GT.0) THEN
                  FILS(DADI) = -FRERE(I)  
                ELSE
                  FILS(DADI) = 0
                ENDIF
                GOTO 90
              ENDIF
             ENDIF
             IN = IFSON
  77         INS = IN
             IN = FRERE(IN)
             IF (IN.NE.I) GOTO 77
             IF (FILS(I) .LT.0) THEN
               FRERE(INS) = -FILS(I)
             ELSE
               FRERE(INS) = FRERE(I)  
               GOTO 90
             ENDIF
  78         CONTINUE
             IN = -FILS(I)
  79         INO = IN
             IN = FRERE(IN)
             IF (IN.GT.0) GOTO 79
             FRERE(INO) = FRERE(I)
  90         CONTINUE
             NODE(DADI) = NODE(DADI)+ NODE(I) 
             NV(DADI)   = NV(DADI) +  NODE(I) 
             NA(IL+1)   = NA(IL+1) + NA(IL)
             GOTO 120
           ENDIF
        ENDIF
#endif
        NE(IS) = NE(IS) + NODE(I) 
        IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        ND(IS) = NV(I)
        NODE(I) = IS
        IPS(I) = IPERM
        IPERM = IPERM + 1
        IN = I
  777   IF (SUBORD(IN).EQ.0) GO TO 778
          IN = SUBORD(IN)
          NODE(IN) = IS
          IPS(IN) = IPERM
          IPERM = IPERM + 1
          GO TO 777
  778   IF (NA(IS).LE.0) GO TO 110
#if defined(NOAMALGTOFATHER)
        IF (   (KEEP60.NE.0).AND.
     &         (NE(IS).EQ.ND(IS)) ) GOTO 110
        IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) THEN
           GO TO 100
        ENDIF
        IF(NAMALG(IS-1) .GE. NAMALGMAX) THEN
           GOTO 110
        ENDIF
        IF ((NE(IS-1).GE.NEMIN).AND.
     &         (NE(IS).GE.NEMIN) ) GO TO 110
        IF (2*NE(IS-1)*(ND(IS)-ND(IS-1)+NE(IS-1)).GE.
     &    ((ND(IS)+NE(IS-1))*
     &    (ND(IS)+NE(IS-1))*NEMIN/100)) GO TO 110
        NAMALG(IS-1) = NAMALG(IS-1)+1
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        ND(IS-1) = ND(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
        NODE(I) = IS-1
        IFSON = -FILS(I)
        IN = IFSON
 102    INO = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 102
        NV(INO) = 0
        IN = I
  888   IF (SUBORD(IN).EQ.0) GO TO 889
        IN = SUBORD(IN)
        GO TO 888
  889   SUBORD(IN) = INO
      INOS = -FILS(INO)
      IF (IFSON.EQ.INO) THEN 
         FILS(I) = -INOS
         GO TO 107
      ENDIF
      IN = IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
        IF (INOS.EQ.0) THEN
          FRERE(INS) = -I
          GO TO 120
        ELSE
          FRERE(INS) =  INOS
        ENDIF
 107    IN = INOS
        IF (IN.EQ.0) GO TO 120
 108    INT1 = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 108
        FRERE(INT1) = -I
        GO TO 120
#endif
  110   IS = IS + 1
  120   IB = FRERE(I)
        IF (IB.GE.0) THEN
          IF (IB.GT.0) NA(IL) = 0
          I = IB
        ELSE
          I = -IB
          IL = IL + 1
        ENDIF
  160 CONTINUE
      NSTEPS = IS - 1
      DO I=1, N
        IF (NV(I).EQ.0) THEN
          FRERE(I) = N+1
          NFSIZ(I) = 0
        ELSE
          NFSIZ(I) = ND(NODE(I))
          IF (SUBORD(I) .NE.0) THEN
           INOS = -FILS(I)  
           INO = I
           DO WHILE (SUBORD(INO).NE.0) 
             IS = SUBORD(INO)
             FILS(INO) = IS
             INO = IS
           END DO
           FILS(INO) = -INOS
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ANA_LNEW
#endif
      SUBROUTINE DMUMPS_ANA_M(NE, ND, NSTEPS,
     & MAXFR, MAXELIM, K50, SIZEFAC_TOT, MAXNPIV,
     & K5,K6,PANEL_SIZE,K253)
      IMPLICIT NONE
      INTEGER, INTENT(in)    :: NSTEPS, K50, K253, K5, K6
      INTEGER, INTENT(in)    :: NE(NSTEPS), ND(NSTEPS)
      INTEGER, INTENT(out)   :: MAXNPIV, PANEL_SIZE
      INTEGER, INTENT(out)   :: MAXFR, MAXELIM
      INTEGER(8), INTENT(out):: SIZEFAC_TOT
      INTEGER ITREE, NFR, NELIM
      INTEGER LKJIB
      INTEGER(8) :: SIZEFAC
      LKJIB        = max(K5,K6)
      MAXFR        = 0
      MAXELIM      = 0
      MAXNPIV      = 0
      PANEL_SIZE   = 0
      SIZEFAC_TOT  = 0_8
      DO ITREE=1,NSTEPS
        NELIM = NE(ITREE)
        NFR = ND(ITREE) + K253
        IF (NFR.GT.MAXFR)         MAXFR   = NFR
        IF (NFR-NELIM.GT.MAXELIM) MAXELIM = NFR - NELIM
        IF (NELIM .GT. MAXNPIV) THEN
           MAXNPIV = NELIM
        ENDIF
        IF (K50.EQ.0) THEN
          SIZEFAC = (2_8*int(NFR,8) - int(NELIM,8))*int(NELIM,8)
          PANEL_SIZE = max(PANEL_SIZE, NFR*(LKJIB+1))
        ELSE
         SIZEFAC = int(NFR,8) * int(NELIM,8)
         PANEL_SIZE = max(PANEL_SIZE, NELIM*(LKJIB+1))
         PANEL_SIZE = max(PANEL_SIZE, (NFR-NELIM)*(LKJIB+1))
        ENDIF
        SIZEFAC_TOT = SIZEFAC_TOT + SIZEFAC
      END DO
      RETURN
      END SUBROUTINE DMUMPS_ANA_M
      SUBROUTINE DMUMPS_ANA_R( N, FILS, FRERE,
     & NSTK, NA )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: FILS(N), FRERE(N)
      INTEGER, INTENT(OUT) ::  NSTK(N), NA(N) 
      INTEGER NBROOT, NBLEAF, ILEAF, I, IN, ISON
      NA   = 0
      NSTK = 0
      NBROOT  = 0
      ILEAF   = 1
      DO 11 I=1,N
         IF (FRERE(I).EQ. N+1) CYCLE
         IF (FRERE(I).EQ.0) NBROOT = NBROOT + 1
         IN = I
 12      IN = FILS(IN)
         IF (IN.GT.0) GO TO 12
         IF (IN.EQ.0) THEN 
            NA(ILEAF) = I
            ILEAF     = ILEAF + 1
            CYCLE
         ENDIF
         ISON = -IN
 13      NSTK(I) = NSTK(I) + 1
         ISON = FRERE(ISON)
         IF (ISON.GT.0) GO TO 13
 11   CONTINUE
      NBLEAF = ILEAF-1
      IF (N.GT.1) THEN
         IF (NBLEAF.GT.N-2) THEN
            IF (NBLEAF.EQ.N-1) THEN
               NA(N-1) = -NA(N-1)-1
               NA(N)   = NBROOT
            ELSE
               NA(N) = -NA(N)-1
            ENDIF
         ELSE
            NA(N-1) = NBLEAF
            NA(N)   = NBROOT
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_ANA_R
      SUBROUTINE DMUMPS_DIAG_ANA
     &( MYID, COMM, KEEP,KEEP8, INFO, INFOG, RINFO, RINFOG, ICNTL )
      IMPLICIT NONE
      INTEGER COMM, MYID, KEEP(500), INFO(80), ICNTL(60), INFOG(80)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION RINFO(40), RINFOG(40)
      INCLUDE 'mpif.h'
      INTEGER MASTER, MPG
      PARAMETER( MASTER = 0 )
      MPG = ICNTL(3)
      IF ( MYID.eq.MASTER.and.MPG.GT.0.AND.ICNTL(4).GE.2) THEN
       WRITE(MPG, 99992) INFO(1), INFO(2),
     &  KEEP8(109), KEEP8(111), INFOG(4),
     &  INFOG(5), KEEP(28), INFOG(32), INFOG(7), KEEP(23), 
     &  ICNTL(7), KEEP(12), 
     &  ICNTL(18),
     &  KEEP(56), KEEP(61), RINFOG(1)
       IF (KEEP(95).GT.1)             
     &      WRITE(MPG, 99993) KEEP(95) 
       IF (KEEP(54).GT.0) WRITE(MPG, 99994) KEEP(54)
       IF (KEEP(60).GT.0) WRITE(MPG, 99995) KEEP(60)
       IF (KEEP(253).GT.0)  WRITE(MPG, 99996) KEEP(253)
      ENDIF
      RETURN
99992 FORMAT(/'Leaving analysis phase with  ...'/
     &       ' INFOG(1)                                       =',I16/
     &       ' INFOG(2)                                       =',I16/
     &       ' -- (20) Number of entries in factors (estim.)  =',I16/
     &       ' --  (3) Real space for factors    (estimated)  =',I16/
     &       ' --  (4) Integer space for factors (estimated)  =',I16/
     &       ' --  (5) Maximum frontal size      (estimated)  =',I16/
     &       ' --  (6) Number of nodes in the tree            =',I16/
     &       ' -- (32) Type of analysis effectively used      =',I16/
     &       ' --  (7) Ordering option effectively used       =',I16/
     &       ' ICNTL (6) Maximum transversal option           =',I16/
     &       ' ICNTL (7) Pivot order option                   =',I16/
     &       ' ICNTL(14) Percentage of memory relaxation      =',I16/
     &       ' ICNTL(18) Distributed input matrix (on if >0)  =',I16/
     &       ' Number of level 2 nodes                        =',I16/
     &       ' Number of split nodes                          =',I16/
     &       ' RINFOG(1) Operations during elimination (estim)=',
     &                                                        1PD10.3)
99993 FORMAT(' Ordering compressed/constrained (ICNTL(12))    =',I16)
99994 FORMAT(' Distributed matrix entry format (ICNTL(18))    =',I16)
99995 FORMAT(' Effective Schur option (ICNTL(19))             =',I16)
99996 FORMAT(' Forward solution during factorization, NRHS    =',I16)
      END SUBROUTINE DMUMPS_DIAG_ANA
      SUBROUTINE DMUMPS_CUTNODES
     &           ( N, FRERE, FILS, NFSIZ, SIZEOFBLOCKS, LSIZEOFBLOCKS,
     &             NSTEPS, NSLAVES, 
     &             KEEP, KEEP8, SPLITROOT, MP, LDIAG, INFO1, INFO2 )
      IMPLICIT NONE
      INTEGER N, NSTEPS, NSLAVES, KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER FRERE( N ), FILS( N ), NFSIZ( N )
      INTEGER LSIZEOFBLOCKS
      INTEGER SIZEOFBLOCKS(LSIZEOFBLOCKS)
      LOGICAL SPLITROOT
      INTEGER MP, LDIAG  
      INTEGER INFO1, INFO2
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPOOL 
      INTEGER INODE, DEPTH, I, IBEG, IEND, IIPOOL, NROOT
      INTEGER MAX_DEPTH, ISON, TOT_CUT, MAX_CUT, STRAT
      INTEGER(8) :: K79
      INTEGER NFRONT, K82, allocok
      LOGICAL BLKON
      BLKON = .NOT.(SIZEOFBLOCKS(1).EQ.-1) 
      K79  = KEEP8(79)
      K82  = abs(KEEP(82))
      STRAT= KEEP(62)
      IF (KEEP(210).EQ.1) THEN
        MAX_DEPTH = 2*NSLAVES*K82
        STRAT     = STRAT/4
      ELSE
        IF (( NSLAVES .eq. 1 ).AND. (.NOT. SPLITROOT) ) RETURN
        IF (NSLAVES.EQ.1) THEN
          MAX_DEPTH=1
        ELSE
          MAX_DEPTH = int( log( dble( NSLAVES - 1 ) ) 
     &                 / log(2.0D0) )
        ENDIF
      ENDIF
      ALLOCATE(IPOOL(NSTEPS+1), stat=allocok)
      IF (allocok.GT.0) THEN
        INFO1= -7
        INFO2= NSTEPS+1
        RETURN
      ENDIF
      NROOT = 0
      DO INODE = 1, N
        IF ( FRERE(INODE) .eq. 0 ) THEN
          NROOT = NROOT + 1
          IPOOL( NROOT ) = INODE
        END IF
      END DO
      IBEG = 1
      IEND = NROOT
      IIPOOL   = NROOT + 1
      IF (SPLITROOT) THEN
         MAX_DEPTH=0 
      ENDIF
      DO DEPTH = 1, MAX_DEPTH
        DO I = IBEG, IEND
          INODE = IPOOL( I )
          ISON = INODE
          DO WHILE ( ISON .GT. 0 )
            ISON = FILS( ISON )
          END DO
          ISON = - ISON
          DO WHILE ( ISON .GT. 0 )
            IPOOL( IIPOOL ) = ISON
            IIPOOL = IIPOOL + 1
            ISON = FRERE( ISON )
          END DO
        END DO
        IPOOL( IBEG ) = -IPOOL( IBEG )
        IBEG = IEND + 1
        IEND = IIPOOL - 1
      END DO
      IPOOL( IBEG ) = -IPOOL( IBEG )
      TOT_CUT = 0
      IF (SPLITROOT) THEN
        MAX_CUT = NROOT*max(K82,2)
        INODE = abs(IPOOL(1))
        NFRONT = NFSIZ( INODE )
        K79 = max(
     &         int(NFRONT,8)*int(NFRONT,8)/(int(K82+1,8)*int(K82+1,8)),
     &         9_8)
        IF (KEEP(53).NE.0) THEN
          MAX_CUT =  NFRONT
          K79 = 121_8*121_8
        ELSE
          K79 = min(2000_8*2000_8,K79)
          IF (KEEP(376) .EQ. 1) THEN
            K79 = min(int(KEEP(9)+1,8)*int(KEEP(9)+1,8),K79)
          ENDIF
        ENDIF
      ELSE
         MAX_CUT = 2 * NSLAVES
         IF (KEEP(210).EQ.1) THEN
            MAX_CUT = 4 * (MAX_CUT + 4)
         ENDIF
      ENDIF
      DEPTH   = -1
      DO I = 1, IIPOOL - 1
        INODE = IPOOL( I )
        IF ( INODE .LT. 0 ) THEN
          INODE = -INODE
          DEPTH = DEPTH + 1
        END IF
        CALL DMUMPS_SPLIT_1NODE
     &           ( INODE, N, FRERE, FILS, NFSIZ, NSTEPS, NSLAVES,
     &             KEEP,KEEP8, TOT_CUT, STRAT, DEPTH, 
     &             K79, SPLITROOT, MP, LDIAG, 
     &             BLKON, SIZEOFBLOCKS, LSIZEOFBLOCKS )
        IF ( TOT_CUT > MAX_CUT )  EXIT
      END DO
      KEEP(61) = TOT_CUT
      DEALLOCATE(IPOOL)
      RETURN
      END SUBROUTINE DMUMPS_CUTNODES
      RECURSIVE SUBROUTINE DMUMPS_SPLIT_1NODE
     & ( INODE, N, FRERE, FILS, NFSIZ, NSTEPS, NSLAVES, KEEP,KEEP8,
     &   TOT_CUT, STRAT, DEPTH, K79, SPLITROOT, MP, LDIAG,
     &   BLKON, SIZEOFBLOCKS, LSIZEOFBLOCKS )
      IMPLICIT NONE
      INTEGER(8) :: K79
      INTEGER INODE, N, NSTEPS, NSLAVES, KEEP(500), STRAT, 
     &        DEPTH, TOT_CUT, MP, LDIAG
      INTEGER(8) KEEP8(150)
      INTEGER FRERE( N ), FILS( N ), NFSIZ( N )
      LOGICAL SPLITROOT
      LOGICAL BLKON
      INTEGER LSIZEOFBLOCKS
      INTEGER SIZEOFBLOCKS(LSIZEOFBLOCKS)
      INTEGER I, IN, NPIV, NFRONT, NSLAVES_ESTIM
      DOUBLE PRECISION WK_SLAVE, WK_MASTER
      INTEGER INODE_SON, INODE_FATH, IN_SON, IN_FATH, IN_GRANDFATH
      INTEGER NPIV_COMPG, NPIV_SON_COMPG, NPIV_FATH_COMPG
      INTEGER NPIV_SON, NPIV_FATH, NPIV_TEMP
      INTEGER NCB, NSLAVESMIN, NSLAVESMAX
      INTEGER  MUMPS_BLOC2_GET_NSLAVESMIN,
     &         MUMPS_BLOC2_GET_NSLAVESMAX
      EXTERNAL  MUMPS_BLOC2_GET_NSLAVESMIN,
     &         MUMPS_BLOC2_GET_NSLAVESMAX
      IF  ( (KEEP(210).EQ.1.AND.KEEP(60).EQ.0) .OR.
     &       (SPLITROOT) ) THEN
        IF ( FRERE ( INODE ) .eq. 0 ) THEN 
          NFRONT = NFSIZ( INODE )
          NPIV = NFRONT
          IF (BLKON) THEN
           IN = INODE
           NPIV_COMPG = 0 
           DO WHILE( IN > 0 )
             NPIV_COMPG = NPIV_COMPG + 1
             IN = FILS( IN )
           ENDDO
          ELSE 
            NPIV_COMPG = NPIV
          ENDIF
          NCB = 0
          IF ( int(NFRONT,8)*int(NFRONT,8).GT.K79
     &    ) THEN 
           GOTO 333
          ENDIF
        ENDIF
      ENDIF
      IF ( FRERE ( INODE ) .eq. 0 ) RETURN
      NFRONT = NFSIZ( INODE )
      IN = INODE
      NPIV       = 0
      NPIV_COMPG = 0 
      DO WHILE( IN > 0 )
        IF (BLKON) THEN
          NPIV = NPIV + SIZEOFBLOCKS(IN)
        ENDIF
        NPIV_COMPG = NPIV_COMPG + 1
        IN = FILS( IN )
      END DO
      IF (.NOT.BLKON) NPIV = NPIV_COMPG
      NCB = NFRONT - NPIV
      IF ( (NFRONT - (NPIV/2)) .LE. KEEP(9)) RETURN
      IF ((KEEP(50) == 0.and.int(NFRONT,8) * int(NPIV,8) > K79 ) .OR.
     &(KEEP(50) .NE.0.and.int(NPIV,8) * int(NPIV,8) > K79 )) GOTO 333
      IF (KEEP(210).EQ.1) THEN
        NSLAVESMIN    = 1   
        NSLAVESMAX    = 64  
        NSLAVES_ESTIM = 32+NSLAVES
      ELSE
        NSLAVESMIN = MUMPS_BLOC2_GET_NSLAVESMIN 
     &         ( NSLAVES, KEEP(48), KEEP8(21), KEEP(50),
     &         NFRONT, NCB, KEEP(375), KEEP(119))
        NSLAVESMAX = MUMPS_BLOC2_GET_NSLAVESMAX 
     &        ( NSLAVES, KEEP(48), KEEP8(21), KEEP(50),
     &          NFRONT, NCB, KEEP(375), KEEP(119))
        NSLAVES_ESTIM = max (1, 
     &   nint( dble(NSLAVESMAX-NSLAVESMIN)/dble(3) )
     &                    )
        NSLAVES_ESTIM = min (NSLAVES_ESTIM, NSLAVES-1)
      ENDIF
      IF ( KEEP(50) .eq. 0 ) THEN
       WK_MASTER = 0.6667D0 * 
     &                dble(NPIV)*dble(NPIV)*dble(NPIV) +
     &                dble(NPIV)*dble(NPIV)*dble(NCB)
       WK_SLAVE  = dble( NPIV ) * dble( NCB ) *
     &         ( 2.0D0 * dble(NFRONT) - dble(NPIV) )
     &         / dble(NSLAVES_ESTIM)
      ELSE
       WK_MASTER = dble(NPIV)*dble(NPIV)*dble(NPIV) / dble(3)
       WK_SLAVE  = 
     &           (dble(NPIV)*dble(NCB)*dble(NFRONT)) 
     &           / dble(NSLAVES_ESTIM)
      ENDIF
      IF (KEEP(210).EQ.1) THEN
        IF ( dble( 100 + STRAT )
     &        * WK_SLAVE / dble(100) .GE. WK_MASTER ) RETURN
      ELSE
        IF ( dble( 100 + STRAT * max( DEPTH-1, 1 ) )
     &        * WK_SLAVE / dble(100) .GE. WK_MASTER ) RETURN
      ENDIF
 333  CONTINUE
      IF (NPIV .LE. 1 ) RETURN
       NPIV_SON  = max(NPIV/2,1)
       NPIV_FATH = NPIV - NPIV_SON
       IF (SPLITROOT) THEN
          IF (NCB .NE .0) THEN
           WRITE(*,*) "Error splitting"
           CALL MUMPS_ABORT()
         ENDIF
         NPIV_FATH = min(int(sqrt(dble(K79))), int(NPIV/2))
         NPIV_SON  = NPIV - NPIV_FATH
       ENDIF
       INODE_SON = INODE
       IF (BLKON) THEN
         NPIV_TEMP       = 0
         NPIV_SON_COMPG  = 0
         IN_SON = INODE
         DO WHILE (IN_SON > 0) 
           NPIV_TEMP = NPIV_TEMP + SIZEOFBLOCKS(IN_SON)
           NPIV_SON_COMPG  = NPIV_SON_COMPG +1
           IF (NPIV_TEMP.GE.NPIV_SON) EXIT
           IN_SON = FILS( IN_SON )
         END DO
         NPIV_FATH_COMPG = NPIV_COMPG - NPIV_SON_COMPG
         NPIV_SON        = NPIV_TEMP
         NPIV_FATH       = NPIV - NPIV_SON 
       ELSE
         NPIV_SON_COMPG  = NPIV_SON
         NPIV_FATH_COMPG = NPIV_FATH
         IN_SON = INODE
         DO I = 1, NPIV_SON_COMPG - 1
           IN_SON = FILS( IN_SON )
         END DO
       ENDIF
       IF (NPIV_FATH_COMPG.EQ.0)  RETURN
       NSTEPS  = NSTEPS + 1
       TOT_CUT = TOT_CUT + 1
       INODE_FATH = FILS( IN_SON )
       IF ( INODE_FATH .LT. 0 ) THEN
       write(*,*) 'Error: INODE_FATH < 0 ', INODE_FATH
       END IF
       IN_FATH = INODE_FATH
       DO WHILE ( FILS( IN_FATH ) > 0 )
         IN_FATH = FILS( IN_FATH )
       END DO
       FRERE( INODE_FATH ) = FRERE( INODE_SON )
       FRERE( INODE_SON  ) = - INODE_FATH
       FILS ( IN_SON     ) = FILS( IN_FATH )
       FILS ( IN_FATH    ) = - INODE_SON
       IN = FRERE( INODE_FATH )
       DO WHILE ( IN > 0 )
           IN = FRERE( IN )
       END DO
       IF ( IN .eq. 0 )  GO TO 10
       IN = -IN
       DO WHILE ( FILS( IN ) > 0 )
           IN = FILS( IN )
       END DO
       IN_GRANDFATH = IN
       IF ( FILS( IN_GRANDFATH ) .eq. - INODE_SON ) THEN
           FILS( IN_GRANDFATH ) = -INODE_FATH
       ELSE
           IN = IN_GRANDFATH
           IN = - FILS ( IN )
           DO WHILE ( FRERE( IN ) > 0 )
             IF ( FRERE( IN ) .eq. INODE_SON ) THEN
               FRERE( IN ) = INODE_FATH
               GOTO 10
             END IF
             IN = FRERE( IN )
           END DO
           WRITE(*,*) 'ERROR 2 in SPLIT NODE',
     &          IN_GRANDFATH, IN, FRERE(IN)
       END IF
 10    CONTINUE
       NFSIZ(INODE_SON) = NFRONT
       NFSIZ(INODE_FATH) = NFRONT - NPIV_SON
       KEEP(2) = max( KEEP(2), NFRONT - NPIV_SON )
       IF (SPLITROOT) THEN
         RETURN
       ENDIF
        CALL DMUMPS_SPLIT_1NODE
     &  ( INODE_FATH, N, FRERE, FILS, NFSIZ, NSTEPS,
     &   NSLAVES, KEEP,KEEP8, TOT_CUT, STRAT, DEPTH, 
     &   K79, SPLITROOT, MP, LDIAG, 
     &   BLKON, SIZEOFBLOCKS, LSIZEOFBLOCKS )
      IF (.NOT. SPLITROOT) THEN
        CALL DMUMPS_SPLIT_1NODE
     &   ( INODE_SON, N, FRERE, FILS, NFSIZ, NSTEPS,
     &   NSLAVES, KEEP,KEEP8, TOT_CUT, STRAT, DEPTH, 
     &   K79, SPLITROOT, MP, LDIAG, 
     &   BLKON, SIZEOFBLOCKS, LSIZEOFBLOCKS )
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SPLIT_1NODE
      SUBROUTINE DMUMPS_ANA_GNEW
     & (N, NZ, IRN, ICN, IW, LW, IPE, LEN,
     & IQ, FLAG, IWFR,
     & NRORM, NIORM, IFLAG,IERROR, ICNTL, 
     & symmetry, SYM, NBQD, AvgDens,
     & KEEP264, KEEP265, PRINTSTAT, 
     & INPLACE64_GRAPH_COPY 
     & )
      IMPLICIT NONE
      INTEGER, intent(in)    :: N, SYM
      INTEGER(8), intent(in) :: LW
      INTEGER(8), intent(in) :: NZ
      INTEGER, intent(in)    :: ICNTL(60)
      INTEGER, intent(in)    :: IRN(NZ), ICN(NZ) 
      INTEGER, intent(out)   :: IERROR, symmetry
      INTEGER, intent(out)   :: NBQD, AvgDens
      INTEGER, intent(out)   :: LEN(N), IW(LW)
      INTEGER(8), intent(out):: IWFR
      INTEGER(8), intent(out):: NRORM, NIORM
      INTEGER(8), intent(out):: IPE(N+1)
      INTEGER, intent(inout) :: IFLAG, KEEP264, KEEP265
      INTEGER(8), intent(out):: IQ(N)
      INTEGER, intent(out)   :: FLAG(N)
      LOGICAL, intent(in)    :: PRINTSTAT
      LOGICAL, intent(inout) :: INPLACE64_GRAPH_COPY
      INTEGER    :: MP, MPG, I, J, N1
      INTEGER    :: NBERR, THRESH
      INTEGER(8) :: K8, K1, K2, LAST, NDUP
      INTEGER(8) :: NZOFFA, NDIAGA, L, N8
      DOUBLE PRECISION       :: RSYM
      INTRINSIC nint
      MP = ICNTL(2)
      MPG= ICNTL(3)
      NZOFFA = 0_8
      NDIAGA = 0
      IERROR = 0
      N8     = int(N,8)
      DO I=1,N+1
        IPE(I) = 0_8
      ENDDO
      IF (KEEP264.EQ.0) THEN
       IF ((SYM.EQ.0).AND.(KEEP265.EQ.-1)) THEN
        DO K8=1_8,NZ
         I = IRN(K8)
         J = ICN(K8)
         IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     &                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
         ELSE
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1_8
           NZOFFA  = NZOFFA + 1_8
          ELSE
           NDIAGA = NDIAGA + 1_8
          ENDIF
         ENDIF
        ENDDO
       ELSE
        DO K8=1_8,NZ
         I = IRN(K8)
         J = ICN(K8)
         IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     &                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
         ELSE
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1_8
           IPE(J) = IPE(J) + 1_8
           NZOFFA  = NZOFFA + 1_8
          ELSE
           NDIAGA = NDIAGA + 1_8
          ENDIF
         ENDIF
        ENDDO
       ENDIF
      ELSE
       IF ((SYM.EQ.0).AND.(KEEP265.EQ.-1)) THEN
        DO K8=1_8,NZ
         I = IRN(K8)
         J = ICN(K8)
         IF (I.EQ.J) THEN
           NDIAGA = NDIAGA + 1_8
         ELSE 
           IPE(I) = IPE(I) + 1_8
           NZOFFA = NZOFFA + 1_8
         ENDIF
        ENDDO
       ELSE
        DO K8=1_8,NZ
         I = IRN(K8)
         J = ICN(K8)
         IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1_8
           IPE(J) = IPE(J) + 1_8
           NZOFFA  = NZOFFA + 1_8
         ELSE
           NDIAGA = NDIAGA + 1_8
         ENDIF
        ENDDO
       ENDIF
      ENDIF
      NIORM  = NZOFFA + 3_8*N8
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IF (mod(IFLAG,2) .EQ. 0) IFLAG = IFLAG+1
         IF ((MP.GT.0).AND.(ICNTL(4).GE.2))  THEN 
          WRITE (MP,99999) 
          DO 70 K8=1_8,NZ
           I = IRN(K8)
           J = ICN(K8)
           IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     &                            .OR.(J.LT.1)) THEN
            NBERR = NBERR + 1
            IF (NBERR.LE.10)  THEN
               IF (mod(K8,10_8).GT.3_8 .OR. mod(K8,10_8).EQ.0_8 .OR.
     &             (10_8.LE.K8 .AND. K8.LE.20_8)) THEN
                 WRITE (MP,'(I16,A,I10,A,I10,A)')
     &             K8,'th entry (in row',I,' and column',J,') ignored'
               ELSE
                 IF (mod(K8,10_8).EQ.1_8) 
     &             WRITE(MP,'(I16,A,I10,A,I10,A)')
     &             K8,'st entry (in row',I,' and column',J,') ignored'
                 IF (mod(K8,10_8).EQ.2_8) 
     &             WRITE(MP,'(I16,A,I10,A,I10,A)')
     &             K8,'nd entry (in row',I,' and column',J,') ignored'
                 IF (mod(K8,10_8).EQ.3_8) 
     &             WRITE(MP,'(I16,A,I10,A,I10,A)')
     &             K8,'rd entry (in row',I,' and column',J,') ignored'
               ENDIF
            ELSE
               GO TO 100
            ENDIF
           ENDIF
   70     CONTINUE
         ENDIF
      ENDIF
  100 NRORM = NIORM - 2_8*N8
      IQ(1) = 1_8
      N1 = N - 1
      IF (N1.GT.0) THEN
        DO I=1,N1
            IQ(I+1) = IPE(I) + IQ(I) 
        ENDDO
      ENDIF
      LAST = max(IPE(N)+IQ(N)-1,IQ(N))
      FLAG(1:N) = 0
      IPE(1:N)  = IQ(1:N)
      IW(1:LAST) = 0
      IWFR = LAST + 1_8
      IF (KEEP264 .EQ. 0) THEN
       IF ((SYM.EQ.0).AND.(KEEP265.EQ.-1)) THEN
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          IF (I.NE.J) THEN
              IF ((J.GE.1).AND.(I.LE.N)) THEN
                IW(IQ(I)) = J
                IQ(I)     = IQ(I) + 1
              ENDIF
          ENDIF
        ENDDO
       ELSE IF (KEEP265.EQ.1) THEN
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          IF (I.NE.J) THEN
              IF ((J.GE.1).AND.(I.LE.N)) THEN
                IW(IQ(J)) = I
                IQ(J)     = IQ(J) + 1
                IW(IQ(I)) = J
                IQ(I)     = IQ(I) + 1
              ENDIF
          ENDIF
        ENDDO
       ELSE 
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          IF (I.NE.J) THEN
            IF (I.LT.J) THEN
              IF ((I.GE.1).AND.(J.LE.N)) THEN
                IW(IQ(I)) = -J
                IQ(I)     = IQ(I) + 1 
              ENDIF
            ELSE
              IF ((J.GE.1).AND.(I.LE.N)) THEN
                IW(IQ(J)) = -I
                IQ(J)     = IQ(J) + 1
              ENDIF
            ENDIF
          ENDIF
        ENDDO
       ENDIF 
      ELSE
       IF ((SYM.EQ.0).AND.(KEEP265.EQ.-1)) THEN
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          IF (I.NE.J) THEN
               IW(IQ(I)) = J
               IQ(I)     = IQ(I) + 1
          ENDIF
        ENDDO
       ELSE IF (KEEP265.EQ.1) THEN
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          IF (I.NE.J) THEN
               IW(IQ(J)) = I
               IQ(J)     = IQ(J) + 1
               IW(IQ(I)) = J
               IQ(I)     = IQ(I) + 1
          ENDIF
        ENDDO
       ELSE
        DO K8=1_8,NZ
          I = IRN(K8)
          J = ICN(K8)
          IF (I.NE.J) THEN
            IF (I.LT.J) THEN
              IW(IQ(I)) = -J
              IQ(I)     = IQ(I) + 1 
            ELSE
              IW(IQ(J)) = -I
              IQ(J)     = IQ(J) + 1
            ENDIF
          ENDIF
        ENDDO
       ENDIF
      ENDIF
      IF (KEEP265.EQ.0) THEN
       NDUP = 0_8
       DO I=1,N
        K1 = IPE(I) 
        K2 = IQ(I) - 1_8
        IF (K1.GT.K2) THEN
         LEN(I) = 0
        ELSE
         DO K8=K1,K2
           J     = -IW(K8)
           IF (J.LE.0) EXIT
           IF (FLAG(J).EQ.I) THEN
            NDUP = NDUP + 1_8
            IW(K8) = 0
           ELSE
            L     = IQ(J) 
            IW(L) = I
            IQ(J) = L + 1
            IW(K8)  = J
            FLAG(J) = I
           ENDIF
         END DO  
         LEN(I) = int((IQ(I) - IPE(I)))
        ENDIF
       ENDDO
       IF (NDUP.NE.0_8) THEN
        IWFR = 1_8
        DO I=1,N
         IF (LEN(I).EQ.0) THEN
            IPE(I) = IWFR
            CYCLE
         ENDIF
         K1 = IPE(I) 
         K2 = K1 + LEN(I) - 1
         L = IWFR
         IPE(I) = IWFR
         DO 270 K8=K1,K2
           IF (IW(K8).NE.0) THEN
            IW(IWFR) = IW(K8)
            IWFR     = IWFR + 1_8
           ENDIF
  270    CONTINUE
         LEN(I) = int(IWFR - L)
        ENDDO
       ELSE 
         KEEP265   = 1
       ENDIF
       IPE(N+1) = IPE(N) + int(LEN(N),8)
       IWFR = IPE(N+1)
      ELSE 
       IPE(1) = 1_8
       DO I = 1, N
        LEN(I) = int(IQ(I) - IPE(I))
       ENDDO
       DO I = 1, N
        IPE(I+1) = IPE(I) + int(LEN(I),8)
       ENDDO
       IWFR = IPE(N+1)
      ENDIF  
      symmetry = 100  
      IF (SYM.EQ.0) THEN
       RSYM =  dble(NDIAGA+2_8*NZOFFA - (IWFR-1_8))/
     &            dble(NZOFFA+NDIAGA) 
       IF ((KEEP265.EQ.0) .AND. (NZOFFA - (IWFR-1_8)).EQ.0_8) 
     &    THEN
         KEEP265 = -1
       ENDIF
       symmetry = min(nint (100.0D0*RSYM), 100)
       IF (PRINTSTAT) THEN
         IF ((MPG .GT. 0).AND.(ICNTL(4).GE.2) )
     &   write(MPG,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
         IF (MP.GT.0 .AND. MPG.NE.MP.AND. (ICNTL(4).GE.2) )
     &   write(MP,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
       ENDIF
      ELSE
      ENDIF
      AvgDens = nint(dble(IWFR-1_8)/dble(N))
      THRESH  = AvgDens*50 - AvgDens/10 + 1
      NBQD    = 0
      IF (N.GT.2) THEN
        DO I= 1, N
          J = max(LEN(I),1)
          IF (J.GT.THRESH) NBQD = NBQD+1
        ENDDO
      ENDIF
      INPLACE64_GRAPH_COPY = (LW.GE.2*(IWFR-1_8))
      IF (PRINTSTAT) THEN
       IF (MPG .GT. 0.AND.(ICNTL(4).GE.2))
     &  write(MPG,'(A,1I5)') 
     &  ' Average density of rows/columns =', AvgDens
        IF (MP.GT.0 .AND. MPG.NE.MP.AND.(ICNTL(4).GE.2))
     &  write(MP,'(A,1I5)') 
     &  ' Average density of rows/columns =', AvgDens
      ENDIF
      RETURN
99999 FORMAT (/'*** Warning message from analysis routine ***')
      END SUBROUTINE DMUMPS_ANA_GNEW
      SUBROUTINE DMUMPS_SET_K821_SURFACE
     &     (KEEP821, KEEP2, KEEP48 ,KEEP50, NSLAVES)
      IMPLICIT NONE
      INTEGER NSLAVES, KEEP2, KEEP48, KEEP50
      INTEGER (8) :: KEEP821
      INTEGER(8) KEEP2_SQUARE, NSLAVES8
      NSLAVES8= int(NSLAVES,8)
      KEEP2_SQUARE = int(KEEP2,8) * int(KEEP2,8)
      KEEP821 = max(KEEP821*int(KEEP2,8),1_8)
#if defined(t3e) 
      KEEP821 = min(1500000_8, KEEP821)
#elif defined(SP_)
      KEEP821 = min(3000000_8, KEEP821)
#else
      KEEP821 = min(2000000_8, KEEP821)
#endif
#if defined(t3e) 
      IF (NSLAVES .GT. 64) THEN
         KEEP821 = 
     &        min(8_8*KEEP2_SQUARE/NSLAVES8+1_8, KEEP821)
      ELSE
         KEEP821 = 
     &        min(4_8*KEEP2_SQUARE/NSLAVES8+1_8, KEEP821)
      ENDIF 
#else
      IF (NSLAVES.GT.64) THEN
         KEEP821 = 
     &        min(6_8*KEEP2_SQUARE/NSLAVES8+1_8, KEEP821)
      ELSE
         KEEP821 = 
     &        min(4_8*KEEP2_SQUARE/NSLAVES8+1_8, KEEP821)
      ENDIF
#endif
         IF (KEEP50 .EQ. 0 ) THEN
            KEEP821 = max(KEEP821,(7_8*KEEP2_SQUARE /
     &          4_8 / int(max(NSLAVES-1,1),8)) + int(KEEP2,8))
         ELSE
            KEEP821 = max(KEEP821,(7_8*KEEP2_SQUARE /
     &          4_8 / int(max(NSLAVES-1,1),8)) + int(KEEP2,8))
         ENDIF
      IF (KEEP50 .EQ. 0 ) THEN
#if defined(t3e)
         KEEP821 = max(KEEP821,200000_8)
#else 
         KEEP821 = max(KEEP821,300000_8)
#endif
      ELSE
#if defined(t3e)
         KEEP821 = max(KEEP821,40000_8)
#else 
         KEEP821 = max(KEEP821,80000_8)
#endif
      ENDIF
      KEEP821 = -KEEP821 
      RETURN
      END SUBROUTINE DMUMPS_SET_K821_SURFACE
      SUBROUTINE DMUMPS_MTRANS_DRIVER(JOB,M,N,NE,
     &     IP,IRN,A,LA,NUM,PERM,LIW,IW,LDW,DW,
     &     IPQ8,
     &     ICNTL,CNTL,INFO, INFOMUMPS)
      IMPLICIT NONE
      INTEGER :: NICNTL, NCNTL, NINFO, INFOMUMPS(80)
      PARAMETER (NICNTL=10, NCNTL=10, NINFO=10)
      INTEGER :: JOB,M,N,NUM
      INTEGER(8), INTENT(IN) :: NE, LIW,LDW, LA
      INTEGER(8)             :: IP(N+1), IPQ8(N)
      INTEGER :: IRN(NE),PERM(M),IW(LIW)
      INTEGER :: ICNTL(NICNTL),INFO(NINFO)
      DOUBLE PRECISION :: A(LA)
      DOUBLE PRECISION :: DW(LDW),CNTL(NCNTL)
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IWtemp8
      INTEGER    :: allocok
      INTEGER    :: I,J,WARN1,WARN2,WARN4
      INTEGER(8) :: K
      DOUBLE PRECISION       :: FACT,ZERO,ONE,RINF,RINF2,RINF3
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+0)
      EXTERNAL DMUMPS_MTRANSZ,DMUMPS_MTRANSB,DMUMPS_MTRANSR,
     &         DMUMPS_MTRANSS,DMUMPS_MTRANSW
      INTRINSIC abs,log
      RINF = CNTL(2)
      RINF2 = huge(RINF2)/dble(2*N)
      RINF3 = 0.0D0
      WARN1 = 0
      WARN2 = 0
      WARN4 = 0
      IF (JOB.LT.1 .OR. JOB.GT.6) THEN
         INFO(1) = -1
         INFO(2) = JOB
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
         GO TO 99
      ENDIF
      IF (M.LT.1 .OR. M.LT.N) THEN
         INFO(1) = -2
         INFO(2) = M
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
         GO TO 99
      ENDIF
      IF (N.LT.1) THEN
         INFO(1) = -2
         INFO(2) = N
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
         GO TO 99
      ENDIF
      IF (NE.LT.1) THEN
         INFO(1) = -3
         CALL MUMPS_SET_IERROR(NE,INFO(2))
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
         GO TO 99
      ENDIF
      IF (JOB.EQ.1) K = int(4*N +   M,8)
      IF (JOB.EQ.2) K = int(N + 2*M,8)
      IF (JOB.EQ.3) K = int(8*N + 2*M + NE,8)
      IF (JOB.EQ.4) K = int(N + M,8)
      IF (JOB.EQ.5) K = int(3*N + 2*M,8)
      IF (JOB.EQ.6) K = int(3*N + 2*M + NE,8)
      IF (LIW.LT.K) THEN
         INFO(1) = -4
         CALL MUMPS_SET_IERROR(K,INFO(2))
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
         GO TO 99
      ENDIF
      IF (JOB.GT.1) THEN
         IF (JOB.EQ.2) K = int(      M,8)
         IF (JOB.EQ.3) K = int(1,8)
         IF (JOB.EQ.4) K = int(    2*M,8)
         IF (JOB.EQ.5) K = int(N + 2*M,8)
         IF (JOB.EQ.6) K = int(N + 3*M,8)
         IF (LDW .LT. K) THEN
            INFO(1) = -5
            CALL MUMPS_SET_IERROR(K,INFO(2))
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
            GO TO 99
         ENDIF
      ENDIF
      IF (ICNTL(5).EQ.0) THEN
         DO 3 I = 1,M
            IW(I) = 0
 3       CONTINUE
         DO 6 J = 1,N
            DO 4 K = IP(J),IP(J+1)-1_8
               I = IRN(K)
               IF (I.LT.1 .OR. I.GT.M) THEN
                  INFO(1) = -6
                  INFO(2) = J
                  IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
                  GO TO 99
               ENDIF
               IF (IW(I).EQ.J) THEN
                  INFO(1) = -7
                  INFO(2) = J
                  IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I 
                  GO TO 99
               ELSE
                  IW(I) = J
               ENDIF
 4          CONTINUE
 6       CONTINUE
      ENDIF
      IF (ICNTL(3).GT.0) THEN
         IF (ICNTL(4).EQ.0 .OR. ICNTL(4).EQ.1) THEN
            WRITE(ICNTL(3),9020) JOB,M,N,NE
            IF (ICNTL(4).EQ.0) THEN
               WRITE(ICNTL(3),9021) (IP(J),J=1,min(10,N+1))
               WRITE(ICNTL(3),9022) (IRN(K),K=1_8,min(10_8,NE))
               IF (JOB.GT.1) WRITE(ICNTL(3),9023) 
     &                              (A(K),K=1_8,min(10_8,NE))
            ELSEIF (ICNTL(4).EQ.1) THEN
               WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
               WRITE(ICNTL(3),9022) (IRN(K),K=1_8,NE)
               IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(K),K=1_8,NE)
            ENDIF
            WRITE(ICNTL(3),9024) (ICNTL(J),J=1,NICNTL)
            WRITE(ICNTL(3),9025) (CNTL(J),J=1,NCNTL)
         ENDIF
      ENDIF
      DO 8 I=1,NINFO
         INFO(I) = 0
    8 CONTINUE
      IF (JOB.EQ.1) THEN
         DO 10 J = 1,N
            IW(J) = int(IP(J+1) - IP(J))
 10      CONTINUE
         CALL DMUMPS_MTRANSZ(M,N,IRN,NE,IP,IW(1),PERM,NUM,
     &        IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1))
         GO TO 90
      ENDIF
      IF (JOB.EQ.2) THEN
         DW(1) = max(ZERO,CNTL(1))
         CALL DMUMPS_MTRANSB(M,N,NE,IP,IRN,A,PERM,NUM,
     &        IW(1),IPQ8,IW(N+1),IW(N+M+1),DW,RINF2)
         GO TO 90
      ENDIF
      IF (JOB.EQ.3) THEN
         DO 20 K = 1,NE
            IW(K) = IRN(K)
 20      CONTINUE
         CALL DMUMPS_MTRANSR(N,NE,IP,IW,A)
         FACT = max(ZERO,CNTL(1))
         CALL DMUMPS_MTRANSS(M,N,NE,IP,IW(1),A,PERM,NUM,IW(NE+1),
     &        IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &        IW(NE+5*N+1),IW(NE+5*N+M+1),FACT,RINF2)
         GO TO 90
      ENDIF
      IF ((JOB.EQ.4).OR.(JOB.EQ.5).or.(JOB.EQ.6)) THEN
           ALLOCATE(IWtemp8(M+N+N), stat=allocok)
           IF (allocok.GT.0) THEN
             INFOMUMPS(1) = -7
             INFOMUMPS(2) = M+N+N
             GOTO 90
           ENDIF
      ENDIF
      IF (JOB.EQ.4) THEN
         DO 50 J = 1,N
            FACT = ZERO
            DO 30 K = IP(J),IP(J+1)-1_8
               IF (abs(A(K)).GT.FACT) FACT = abs(A(K))
 30         CONTINUE
            IF(FACT .GT. RINF3) RINF3 = FACT
            DO 40 K = IP(J),IP(J+1)-1_8
               A(K) = FACT - abs(A(K))
 40         CONTINUE
 50      CONTINUE
         DW(1)      = max(ZERO,CNTL(1))
         DW(2)      = RINF3
         IWtemp8(1) = int(JOB,8)
         CALL DMUMPS_MTRANSW(M,N,NE,IP,IRN,A,PERM,NUM,
     &        IWtemp8(1),IW(1),IWtemp8(N+1),IPQ8,IW(N+1),
     &        IWtemp8(2*N+1),
     &        DW(1),DW(M+1),RINF2)
         DEALLOCATE(IWtemp8)  
         GO TO 90
      ENDIF
      IF (JOB.EQ.5 .or. JOB.EQ.6) THEN
         RINF3=ONE
         IF (JOB.EQ.5) THEN
            DO 75 J = 1,N
               FACT = ZERO
               DO 60 K = IP(J),IP(J+1)-1_8
                  IF (A(K).GT.FACT) FACT = A(K)
 60            CONTINUE
               DW(2*M+J) = FACT
               IF (FACT.NE.ZERO) THEN
                  FACT = log(FACT)
                  IF(FACT .GT. RINF3) RINF3=FACT
                  DO 70 K = IP(J),IP(J+1)-1_8
                     IF (A(K).NE.ZERO) THEN
                        A(K) = FACT - log(A(K))
                        IF(A(K) .GT. RINF3) RINF3=A(K)
                     ELSE
                        A(K) = FACT + RINF
                     ENDIF
 70               CONTINUE
               ELSE
                  DO 71 K = IP(J),IP(J+1)-1_8
                     A(K) = ONE
 71               CONTINUE
               ENDIF
 75         CONTINUE
         ENDIF
         IF (JOB.EQ.6) THEN
            DO 175 K = 1,NE
               IW(3*N+2*M+K) = IRN(K)
 175        CONTINUE
            DO 61 I = 1,M
               DW(2*M+N+I) = ZERO
 61         CONTINUE
            DO 63 J = 1,N
               DO 62 K = IP(J),IP(J+1)-1_8
                  I = IRN(K)
                  IF (A(K).GT.DW(2*M+N+I)) THEN
                     DW(2*M+N+I) = A(K)
                  ENDIF
 62            CONTINUE
 63         CONTINUE
            DO 64 I = 1,M
               IF (DW(2*M+N+I).NE.ZERO) THEN
                  DW(2*M+N+I) = 1.0D0/DW(2*M+N+I)
               ENDIF
 64         CONTINUE
            DO 66 J = 1,N
               DO 65 K = IP(J),IP(J+1)-1
                  I = IRN(K)
                  A(K) = DW(2*M+N+I) * A(K)
 65            CONTINUE
 66         CONTINUE
            CALL DMUMPS_MTRANSR(N,NE,IP,IW(3*N+2*M+1),A)
            DO 176 J = 1,N
               IF (IP(J).NE.IP(J+1)) THEN
                  FACT = A(IP(J))
               ELSE
                  FACT = ZERO
               ENDIF
               DW(2*M+J) = FACT
               IF (FACT.NE.ZERO) THEN
                  FACT = log(FACT)
                  DO 170 K = IP(J),IP(J+1)-1_8
                     IF (A(K).NE.ZERO) THEN
                        A(K) = FACT - log(A(K))
                        IF(A(K) .GT. RINF3) RINF3=A(K)
                     ELSE
                        A(K) = FACT + RINF
                     ENDIF
 170              CONTINUE
               ELSE
                  DO 171 K = IP(J),IP(J+1)-1_8
                     A(K) = ONE
 171              CONTINUE
               ENDIF
 176        CONTINUE
         ENDIF
         DW(1) = max(ZERO,CNTL(1))
         RINF3 = RINF3+ONE
         DW(2) = RINF3
         IWtemp8(1) = int(JOB,8)
          IF (JOB.EQ.5) THEN
         CALL DMUMPS_MTRANSW(M,N,NE,IP,IRN,A,PERM,NUM,
     &        IWtemp8(1),IW(1),IWtemp8(N+1),IPQ8,IW(N+1),
     &        IWtemp8(2*N+1),
     &        DW(1),DW(M+1),RINF2)
         ENDIF
         IF (JOB.EQ.6) THEN
         CALL DMUMPS_MTRANSW(M,N,NE,IP,IW(3*N+2*M+1),A,PERM,NUM,
     &        IWtemp8(1),IW(1),IWtemp8(N+1),IPQ8,IW(N+1),
     &        IWtemp8(2*N+1),
     &        DW(1),DW(M+1),RINF2)
         ENDIF
         IF ((JOB.EQ.5).or.(JOB.EQ.6)) THEN
          DEALLOCATE(IWtemp8)
         ENDIF
         IF (JOB.EQ.6) THEN
            DO 79 I = 1,M
               IF (DW(2*M+N+I).NE.0.0D0) THEN
                  DW(I) = DW(I) + log(DW(2*M+N+I))
               ENDIF
 79         CONTINUE
         ENDIF
         IF (NUM.EQ.N) THEN
            DO 80 J = 1,N
               IF (DW(2*M+J).NE.ZERO) THEN
                  DW(M+J) = DW(M+J) - log(DW(2*M+J))
               ELSE
                  DW(M+J) = ZERO
               ENDIF
 80         CONTINUE
         ENDIF
         FACT = 0.5D0*log(RINF2)
         DO 86 I = 1,M
            IF (DW(I).LT.FACT) GO TO 86
            WARN2 = 2
            GO TO 90
 86      CONTINUE 
         DO 87 J = 1,N
            IF (DW(M+J).LT.FACT) GO TO 87
            WARN2 = 2
            GO TO 90
 87      CONTINUE 
      ENDIF
 90   IF (INFOMUMPS(1).LT.0) RETURN
      IF (NUM.LT.N) WARN1 = 1
      IF (JOB.EQ.4 .OR. JOB.EQ.5 .OR. JOB.EQ.6) THEN 
         IF (CNTL(1).LT.ZERO) WARN4 = 4
      ENDIF
      IF (INFO(1).EQ.0) THEN
         INFO(1) = WARN1 + WARN2 + WARN4
         IF (INFO(1).GT.0 .AND. ICNTL(2).GT.0) THEN
            WRITE(ICNTL(2),9010) INFO(1)
            IF (WARN1.EQ.1) WRITE(ICNTL(2),9011)
            IF (WARN2.EQ.2) WRITE(ICNTL(2),9012)
            IF (WARN4.EQ.4) WRITE(ICNTL(2),9014)
         ENDIF
      ENDIF
      IF (ICNTL(3).GE.0) THEN
         IF (ICNTL(4).EQ.0 .OR. ICNTL(4).EQ.1) THEN
            WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
            WRITE(ICNTL(3),9031) NUM
            IF (ICNTL(4).EQ.0) THEN
               WRITE(ICNTL(3),9032) (PERM(J),J=1,min(10,M))
               IF (JOB.EQ.5 .OR. JOB.EQ.6) THEN
                  WRITE(ICNTL(3),9033) (DW(J),J=1,min(10,M))
                  WRITE(ICNTL(3),9034) (DW(M+J),J=1,min(10,N))
               ENDIF
            ELSEIF (ICNTL(4).EQ.1) THEN
               WRITE(ICNTL(3),9032) (PERM(J),J=1,M)
               IF (JOB.EQ.5 .OR. JOB.EQ.6) THEN
                  WRITE(ICNTL(3),9033) (DW(J),J=1,M)
                  WRITE(ICNTL(3),9034) (DW(M+J),J=1,N)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
 99   RETURN
 9001 FORMAT (' ****** Error in DMUMPS_MTRANSA. INFO(1) = ',I2,
     &     ' because ',(A),' = ',I14)
 9004 FORMAT (' ****** Error in DMUMPS_MTRANSA. INFO(1) = ',I2/
     &     '        LIW too small, must be at least ',I14)
 9005 FORMAT (' ****** Error in DMUMPS_MTRANSA. INFO(1) = ',I2/
     &     '        LDW too small, must be at least ',I14)
 9006 FORMAT (' ****** Error in DMUMPS_MTRANSA. INFO(1) = ',I2/
     &     '        Column ',I8,
     &     ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in DMUMPS_MTRANSA. INFO(1) = ',I2/
     &     '        Column ',I8,
     &     ' contains two or more entries with row index ',I8)
 9010 FORMAT (' ****** Warning from DMUMPS_MTRANSA. INFO(1) = ',I2)
 9011 FORMAT ('        - The matrix is structurally singular.')
 9012 FORMAT ('        - Some scaling factors may be too large.')
 9014 FORMAT ('        - CNTL(1) is negative and was treated as zero.')
 9020 FORMAT (' ****** Input parameters for DMUMPS_MTRANSA:'/
     &     ' JOB =',I10/' M   =',I10/' N   =',I10/' NE  =',I14)
 9021 FORMAT (' IP(1:N+1)   = ',8I8/(15X,8I8))
 9022 FORMAT (' IRN(1:NE)   = ',8I8/(15X,8I8))
 9023 FORMAT (' A(1:NE)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9024 FORMAT (' ICNTL(1:10) = ',8I8/(15X,2I8))
 9025 FORMAT (' CNTL(1:10)  = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for DMUMPS_MTRANSA:'/
     &     ' INFO(1:2)   = ',2I8)
 9031 FORMAT (' NUM         = ',I8)
 9032 FORMAT (' PERM(1:M)   = ',8I8/(15X,8I8))
 9033 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
 9034 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END SUBROUTINE DMUMPS_MTRANS_DRIVER
      SUBROUTINE DMUMPS_SUPPRESS_DUPPLI_VAL(N,NZ,IP,IRN,A,FLAG,POSI)
      IMPLICIT NONE
      INTEGER, INTENT(IN)       :: N
      INTEGER(8), INTENT(INOUT) :: NZ
      INTEGER(8), INTENT(INOUT) :: IP(N+1)
      INTEGER, INTENT(INOUT)    :: IRN(NZ)
      DOUBLE PRECISION, INTENT(INOUT)       :: A(NZ)
      INTEGER, INTENT(OUT)         :: FLAG(N)
      INTEGER(8), INTENT(OUT)      :: POSI(N)
      INTEGER    :: ROW, COL
      INTEGER(8) :: K, WR_POS, BEG_COL, SV_POS
      FLAG = 0
      WR_POS = 1_8
      DO COL=1,N
         BEG_COL = WR_POS
         DO K=IP(COL),IP(COL+1)-1_8
            ROW = IRN(K)
            IF(FLAG(ROW) .NE. COL) THEN
               IRN(WR_POS) = ROW
               A(WR_POS) = A(K)
               FLAG(ROW) = COL
               POSI(ROW) = WR_POS
               WR_POS = WR_POS+1
            ELSE
               SV_POS = POSI(ROW)
               A(SV_POS) = A(SV_POS) + A(K)
            ENDIF
         ENDDO
         IP(COL) = BEG_COL
      ENDDO
      IP(N+1) = WR_POS
      NZ = WR_POS-1_8
      RETURN
      END SUBROUTINE DMUMPS_SUPPRESS_DUPPLI_VAL
      SUBROUTINE DMUMPS_SUPPRESS_DUPPLI_STR(N,NZ,IP,IRN,FLAG)
      IMPLICIT NONE
      INTEGER, INTENT(IN)       :: N
      INTEGER(8), INTENT(INOUT) :: NZ
      INTEGER(8), INTENT(INOUT) :: IP(N+1)
      INTEGER, INTENT(INOUT)    :: IRN(NZ)
      INTEGER, INTENT(OUT)         :: FLAG(N)
      INTEGER    :: ROW, COL
      INTEGER(8) :: K, WR_POS, BEG_COL
      FLAG = 0
      WR_POS = 1_8
      DO COL=1,N
         BEG_COL = WR_POS
         DO K=IP(COL),IP(COL+1)-1_8
            ROW = IRN(K)
            IF(FLAG(ROW) .NE. COL) THEN
               IRN(WR_POS) = ROW
               FLAG(ROW) = COL
               WR_POS = WR_POS+1_8
            ENDIF
         ENDDO
         IP(COL) = BEG_COL
      ENDDO
      IP(N+1) = WR_POS
      NZ = WR_POS-1_8
      RETURN
      END SUBROUTINE DMUMPS_SUPPRESS_DUPPLI_STR
      SUBROUTINE DMUMPS_SORT_PERM( N, NA, LNA, NE_STEPS,
     &          PERM, FILS, 
     &          DAD_STEPS, STEP, NSTEPS, 
     &          KEEP60, KEEP20, KEEP38,
     &          INFO)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  ::  N, NSTEPS, LNA
      INTEGER, INTENT(IN)  ::  FILS( N ), STEP(N), NA(LNA)
      INTEGER, INTENT(IN)  ::  DAD_STEPS ( NSTEPS ), NE_STEPS (NSTEPS)
      INTEGER, INTENT(IN)  ::  KEEP60, KEEP20, KEEP38
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER, INTENT(OUT) ::  PERM( N )
      INTEGER  :: IPERM, INODE, IN, ISCHUR
      INTEGER  :: INBLEAF, INBROOT, allocok
      INTEGER, ALLOCATABLE, DIMENSION (:) :: POOL, NSTK
      INBLEAF = NA(1) 
      INBROOT = NA(2) 
      ALLOCATE(POOL(INBLEAF), NSTK(NSTEPS), stat=allocok)
      IF (allocok > 0 ) THEN
        INFO(1) = -7
        INFO(2) = INBLEAF + NSTEPS
        RETURN
      ENDIF
      POOL(1:INBLEAF) = NA(3:2+INBLEAF)
      NSTK(1:NSTEPS) = NE_STEPS(1:NSTEPS)
      ISCHUR = 0
      IF ( KEEP60.GT.0 ) THEN
        ISCHUR = max (KEEP20, KEEP38)
      ENDIF
      IPERM = 1
      DO WHILE ( INBLEAF .NE. 0 )
        INODE = POOL( INBLEAF )
        INBLEAF = INBLEAF - 1
        IN = INODE
        IF (INODE.NE.ISCHUR) THEN
           DO WHILE ( IN .GT. 0 )
             PERM ( IN ) = IPERM
             IPERM = IPERM + 1
             IN = FILS( IN )
           END DO
        ENDIF
        IN = DAD_STEPS(STEP( INODE ))
        IF ( IN .eq. 0 ) THEN
          INBROOT = INBROOT - 1
        ELSE
          NSTK( STEP(IN) ) = NSTK( STEP(IN) ) - 1
          IF ( NSTK( STEP(IN) ) .eq. 0 ) THEN
            INBLEAF = INBLEAF + 1
            POOL( INBLEAF ) = IN
          END IF
        END IF
      END DO
      IF (IPERM.LE.N) THEN
       IF (ISCHUR.GT.0) THEN
         IN = ISCHUR
         DO WHILE ( IN .GT. 0 )
             PERM ( IN ) = IPERM
             IPERM = IPERM + 1
             IN = FILS( IN )
         END DO
       ENDIF
      ENDIF
      DEALLOCATE(POOL, NSTK)
      RETURN
      END SUBROUTINE DMUMPS_SORT_PERM
      SUBROUTINE DMUMPS_EXPAND_TREE_STEPS(  ICNTL, 
     &          N, NBLK, BLKPTR, BLKVAR,
     &          FILS_OLD, FILS_NEW, NSTEPS,
     &          STEP_OLD, STEP_NEW, PAR2_NODES, NB_NIV2,
     &          DAD_STEPS, FRERE_STEPS, 
     &          NA, LNA, LRGROUPS_OLD, LRGROUPS_NEW, 
     &          K20, K38
     &           ) 
      IMPLICIT NONE
      INTEGER, INTENT(IN)    ::  N, NBLK, ICNTL(60), NSTEPS, LNA, 
     &                           NB_NIV2
      INTEGER, INTENT(IN)    ::  BLKPTR(NBLK+1), BLKVAR(N)
      INTEGER, INTENT(IN)    ::  FILS_OLD(NBLK), STEP_OLD(NBLK), 
     &                           LRGROUPS_OLD(NBLK)
      INTEGER, INTENT(OUT)   ::  FILS_NEW(N), STEP_NEW(N), 
     &                           LRGROUPS_NEW(N)
      INTEGER, INTENT(INOUT) ::  DAD_STEPS(NSTEPS), FRERE_STEPS(NSTEPS)
      INTEGER, INTENT(INOUT) ::  NA(LNA), PAR2_NODES(NB_NIV2), K20, K38
      INTEGER :: IB, I, IBFS, IBNB, IFS, INB
      INTEGER NBLEAF, NBROOT, ISTEP, IGROUP
      INTEGER :: II
      IF (K20.GT.0) K20 =  BLKVAR(BLKPTR(K20))
      IF (K38.GT.0) K38 =  BLKVAR(BLKPTR(K38))
      NBLEAF = NA(1)
      NBROOT = NA(2)
      IF (NBLK.GT.1) THEN
       DO I= 3, 3+NBLEAF+NBROOT-1
         IBNB  = NA(I)
         INB   = BLKVAR(BLKPTR(IBNB))
         NA(I) = INB
       ENDDO
      ENDIF
      IF (PAR2_NODES(1).GT.0) THEN
       DO I=1, NB_NIV2
         IBNB = PAR2_NODES(I)
         INB          = BLKVAR(BLKPTR(IBNB))
         PAR2_NODES(I) = INB
       ENDDO
      ENDIF
      DO I= 1, NSTEPS
        IBNB = DAD_STEPS(I)
        IF (IBNB.EQ.0) THEN
         INB = 0
        ELSE
         INB   = BLKVAR(BLKPTR(IBNB))
        ENDIF
        DAD_STEPS(I) = INB
      ENDDO
      DO I= 1, NSTEPS
       IBNB = FRERE_STEPS(I)
       IF (IBNB.EQ.0) THEN
         INB = 0
       ELSE
        INB  = BLKVAR(BLKPTR(abs(IBNB)))
        IF (IBNB.LT.0) INB=-INB
       ENDIF
       FRERE_STEPS(I) = INB
      ENDDO
      DO IB=1, NBLK
        IBFS = FILS_OLD(IB)
        IF (IBFS.EQ.0) THEN
         IFS = 0
        ELSE
         IFS  = BLKVAR(BLKPTR(abs(IBFS)))
         IF (IBFS.LT.0) IFS=-IFS
        ENDIF
        IF (BLKPTR(IB+1)-BLKPTR(IB).EQ.0) CYCLE
        DO II=BLKPTR(IB), BLKPTR(IB+1)-1
          IF (II.LT. BLKPTR(IB+1)-1) THEN
            FILS_NEW(BLKVAR(II))= BLKVAR(II+1)
          ELSE
            FILS_NEW(BLKVAR(II))= IFS
          ENDIF
         ENDDO
      ENDDO
      DO IB=1, NBLK
        ISTEP = STEP_OLD(IB)
        IF (BLKPTR(IB+1)-BLKPTR(IB).EQ.0) CYCLE
        IF (ISTEP.LT.0) THEN 
         DO II=BLKPTR(IB), BLKPTR(IB+1)-1
           STEP_NEW(BLKVAR(II)) =  ISTEP
         ENDDO
        ELSE
          I    = BLKVAR(BLKPTR(IB))
          STEP_NEW(I) = ISTEP
          DO II=BLKPTR(IB)+1, BLKPTR(IB+1)-1
           STEP_NEW(BLKVAR(II)) = -ISTEP
          ENDDO
        ENDIF
      ENDDO
      DO IB=1, NBLK
        IGROUP = LRGROUPS_OLD(IB)
        IF (BLKPTR(IB+1)-BLKPTR(IB).EQ.0) CYCLE
        DO II=BLKPTR(IB), BLKPTR(IB+1)-1
           LRGROUPS_NEW(BLKVAR(II)) =  IGROUP
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_EXPAND_TREE_STEPS
      SUBROUTINE DMUMPS_DIST_AVOID_COPIES(N,NSLAVES,
     &     ICNTL,INFOG, NE, NFSIZ,
     &     FRERE, FILS,
     &     KEEP,KEEP8,PROCNODE,
     &     SSARBR,NBSA,PEAK,IERR
     &           , SIZEOFBLOCKS, LSIZEOFBLOCKS 
     &     )
      USE MUMPS_STATIC_MAPPING
      IMPLICIT NONE
      INTEGER N, NSLAVES, NBSA, IERR
      INTEGER ICNTL(60),INFOG(80),KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NE(N),NFSIZ(N),FRERE(N),FILS(N),PROCNODE(N)
      INTEGER SSARBR(N)
      DOUBLE PRECISION PEAK
      INTEGER, intent(IN) :: LSIZEOFBLOCKS
      INTEGER, intent(IN) :: SIZEOFBLOCKS(LSIZEOFBLOCKS)
      CALL MUMPS_DISTRIBUTE(N,NSLAVES,
     &     ICNTL,INFOG, NE, NFSIZ,
     &     FRERE, FILS,
     &     KEEP,KEEP8,PROCNODE,
     &     SSARBR,NBSA,dble(PEAK),IERR
     &           , SIZEOFBLOCKS, LSIZEOFBLOCKS 
     &     )
      RETURN
      END SUBROUTINE DMUMPS_DIST_AVOID_COPIES
      SUBROUTINE DMUMPS_SET_PROCNODE(INODE, PROCNODE, VALUE, FILS, N)
      INTEGER, intent(in) :: INODE, N, VALUE
      INTEGER, intent(in) :: FILS(N)
      INTEGER, intent(inout) :: PROCNODE(N)
      INTEGER IN
      IN=INODE
      DO WHILE ( IN > 0 )
         PROCNODE( IN ) = VALUE
         IN=FILS( IN )
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_SET_PROCNODE
