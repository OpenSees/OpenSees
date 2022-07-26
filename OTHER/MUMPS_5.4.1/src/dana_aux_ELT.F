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
      SUBROUTINE DMUMPS_ANA_F_ELT(N, NELT, ELTPTR, ELTVAR, LIW,
     &            IKEEP, 
     &            IORD, NFSIZ, FILS, FRERE, 
     &            LISTVAR_SCHUR, SIZE_SCHUR,
     &            ICNTL, INFO, KEEP,KEEP8,
     &            NSLAVES, 
     &            XNODEL, NODEL
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)         
     &            , METIS_OPTIONS
#endif      
     &            )
      USE MUMPS_ANA_ORD_WRAPPERS
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N,  NELT, SIZE_SCHUR, NSLAVES, LIW
      INTEGER, INTENT(IN)    :: ELTPTR(NELT+1) 
      INTEGER, INTENT(IN)    :: ELTVAR(ELTPTR(NELT+1)-1)
      INTEGER, INTENT(IN)    :: LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, INTENT(INOUT) :: IORD
      INTEGER, INTENT(INOUT) :: IKEEP(N,3)
      INTEGER, INTENT(INOUT) :: INFO(80), KEEP(500)
      INTEGER(8), INTENT(INOUT) :: KEEP8(150)
      INTEGER, INTENT(OUT)   :: NFSIZ(N), FILS(N), FRERE(N)
      INTEGER, INTENT(OUT)   :: XNODEL(N+1), NODEL(ELTPTR(NELT+1)-1)
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)      
      INTEGER, INTENT(IN)    :: METIS_OPTIONS(40)
#endif
      INTEGER K,I,L1,L2,NCMPA,IFSON,IN
      INTEGER NEMIN, MPRINT, LP, MP, LDIAG
      INTEGER(8)  :: NZ8, LLIW8, IWFR8
      INTEGER allocok, ITEMP
      LOGICAL PROK, NOSUPERVAR, LPOK
      INTEGER(8) :: K79REF
      PARAMETER(K79REF=12000000_8)
      LOGICAL SPLITROOT
      INTEGER, PARAMETER :: LIDUMMY = 1
      INTEGER :: IDUMMY(1) 
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW2
      INTEGER, DIMENSION(:), ALLOCATABLE :: PARENT
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IWtemp
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IPE8
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3) 
#if defined(metis4) || defined(parmetis3)
      INTEGER :: NUMFLAG
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: NUMFLAG
#endif
      INTEGER :: OPT_METIS_SIZE, METIS_IDX_SIZE
      INTEGER :: IERR
#endif
      INTEGER IDUM 
      EXTERNAL DMUMPS_ANA_G11_ELT, DMUMPS_ANA_G12_ELT,
     &         DMUMPS_ANA_G1_ELT, DMUMPS_ANA_G2_ELT, 
     &         DMUMPS_ANA_G2_ELTNEW,
     &         DMUMPS_ANA_J1_ELT, DMUMPS_ANA_J2_ELT,
     &         DMUMPS_ANA_K,
     &         DMUMPS_ANA_LNEW, DMUMPS_ANA_M,
     &         MUMPS_AMD_ELT
#if defined(OLDDFS)
      EXTERNAL DMUMPS_ANA_L
#endif
        ALLOCATE( IW ( LIW ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO( 1 ) = -7
          INFO( 2 ) = LIW
          GOTO 90 
        ENDIF
        ALLOCATE( IPE8 ( N + 1 ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO( 1 ) = -7
          INFO( 2 ) = (N+1)*KEEP(10)
          GOTO 90 
        ENDIF
        ALLOCATE( PARENT(N), IWtemp ( N, 3 ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO( 1 ) = -7
          INFO( 2 ) = 4*N
          GOTO 90 
        ENDIF
      MPRINT= ICNTL(3)
      LP    = ICNTL(1)
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MP    = ICNTL(3)
      PROK  = ((MP.GT.0).AND.(ICNTL(4).GE.2))
      LDIAG = ICNTL(4)
      IF (KEEP(60).NE.0) THEN
       NOSUPERVAR=.TRUE.
       IF (IORD.GT.1) IORD = 0
      ELSE
       NOSUPERVAR=.FALSE.
      ENDIF
      IF (IORD == 7) THEN
         IF ( N < 10000 ) THEN
           IORD = 0
         ELSE
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
           IORD = 5
#else
           IORD = 0
#endif
         ENDIF
      END IF
#if ! defined(metis) && ! defined(parmetis) && ! defined(metis4) && ! defined(parmetis3)
      IF (IORD == 5) IORD = 0
#endif
      IF (KEEP(1).LT.1) KEEP(1) = 1
      NEMIN = KEEP(1)
      IF (LDIAG.LE.2 .OR. MP.LE.0) GO TO 10
      WRITE (MP,99999) N, NELT, LIW, INFO(1)
      K = min0(10,NELT+1)
      IF (LDIAG.EQ.4) K = NELT+1
      IF (K.GT.0) WRITE (MP,99998) (ELTPTR(I),I=1,K)
      K = min0(10,ELTPTR(NELT+1)-1)
      IF (LDIAG.EQ.4) K = ELTPTR(NELT+1)-1
      IF (K.GT.0) WRITE (MP,99995) (ELTVAR(I),I=1,K)
      K = min0(10,N)
      IF (LDIAG.EQ.4) K = N
      IF (IORD.EQ.1 .AND. K.GT.0) THEN
        WRITE (MP,99997) (IKEEP(I,1),I=1,K)
      ENDIF
   10 L1 = 1
      L2 = L1 + N
      IF (LIW .LT. 3*N) THEN
          INFO(1) = -2002
          INFO(2) = LIW
      ENDIF
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
      IF ( IORD == 5 ) THEN
        IF (LIW .LT. N+N+1) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          GOTO 90 
        ENDIF
      ELSE
#endif
      IF (NOSUPERVAR) THEN
        IF ( LIW .LT. 2*N ) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          GOTO 90 
        END IF
      ELSE
        IF ( LIW .LT.  4*N+4 ) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          GOTO 90 
        END IF
      ENDIF
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
      ENDIF
#endif
      IDUM=0
      CALL DMUMPS_NODEL(NELT, N, ELTPTR(NELT+1)-1, ELTPTR, ELTVAR,
     &           XNODEL, NODEL, IW(L1), IDUM, ICNTL)
      IF (IORD.NE.1 .AND. IORD .NE. 5) THEN
        IORD = 0
        IF (NOSUPERVAR) THEN
          CALL DMUMPS_ANA_G1_ELT(N, NZ8, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              IWtemp(1,2), IW(L1))
        ELSE
         CALL DMUMPS_ANA_G11_ELT(N, NZ8, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              IWtemp(1,2), 4*N+4, IW(L1))
        ENDIF
        LLIW8 = max(NZ8,int(N,8))
        ALLOCATE( IW2(LLIW8), stat = allocok )
        IF (allocok.GT.0) THEN
          INFO(1) = -7
          CALL MUMPS_SET_IERROR(LLIW8, INFO(2))
          GOTO 90 
        ENDIF
        IF (NOSUPERVAR) THEN
         CALL DMUMPS_ANA_G2_ELT(N, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              IW2, LLIW8, IPE8, IWtemp(1,2),
     &              IW(L1), IWFR8)
        ELSE
         CALL DMUMPS_ANA_G12_ELT(N, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              IW2, LLIW8, IPE8, IWtemp(1,2),
     &              IW(L1), IWFR8)
        ENDIF
        IF (NOSUPERVAR) THEN
         CALL MUMPS_HAMD(N, LLIW8, IPE8, IWFR8, IWtemp(1,2), IW2,
     &   IW(L1), IKEEP,
     &   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), IW(L2), IWtemp(1,3),
     &   IWtemp,
     &   LISTVAR_SCHUR, SIZE_SCHUR)
         IF (KEEP(60) == 1) THEN
           KEEP(20) = LISTVAR_SCHUR(1)
         ELSEIF (KEEP(60) == 2 .OR. KEEP(60) == 3 ) THEN
           KEEP(38) = LISTVAR_SCHUR(1)
        ELSE
           WRITE(*,*) "Internal error in DMUMPS_ANA_F_ELT",KEEP(60)
           CALL MUMPS_ABORT()
         ENDIF
        ELSE
         CALL MUMPS_AMD_ELT(N, LLIW8, IPE8, IWFR8, IWtemp(1,2), IW2, 
     &   IW(L1), IKEEP, 
     &   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), IW(L2), IWtemp(1,3),
     &   IWtemp)
        ENDIF
      ELSE
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
        IF (IORD.EQ.5) THEN
         IF (PROK) THEN
          WRITE(MPRINT,'(A)') ' Ordering based on METIS'
         ENDIF
         CALL DMUMPS_ANA_G1_ELT(N, NZ8, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              IWtemp(1,2), IW(L1))
         LLIW8 = max(NZ8,int(N,8))
         ALLOCATE( IW2(LLIW8), stat = allocok )
         IF (allocok.GT.0) THEN
           INFO(1) = -7
           CALL MUMPS_SET_IERROR(LLIW8, INFO(2))
           GOTO 90 
         ENDIF
         CALL DMUMPS_ANA_G2_ELTNEW(N, NELT, ELTPTR(NELT+1)-1, 
     &              ELTPTR, ELTVAR, XNODEL, NODEL,
     &              IW2, LLIW8, IPE8, IWtemp(1,2),
     &              IW(L1), IWFR8)
#if defined(metis4) || defined(parmetis3)
         NUMFLAG = 1
         OPT_METIS_SIZE = 8
#else
         ALLOCATE( NUMFLAG ( N ), stat = IERR )
         IF ( IERR .GT. 0 ) THEN
            INFO( 1 ) = -7
            INFO( 2 ) = N
            GOTO 90             
         ENDIF
         DO I=1,N
            NUMFLAG(I) = 1
         ENDDO
         OPT_METIS_SIZE = 40
#endif
         CALL MUMPS_METIS_IDXSIZE(METIS_IDX_SIZE)
         IF (KEEP(10).EQ.1.AND.METIS_IDX_SIZE.NE.64) THEN
           INFO(1) = -52
           INFO(2) = 1
           GOTO 90
         ENDIF
         IF (METIS_IDX_SIZE .EQ. 32) THEN
           CALL MUMPS_METIS_NODEND_MIXEDto32(N, IPE8, IW2,
#if defined(metis4) || defined(parmetis3)
     &          NUMFLAG,  
#else
     &          NUMFLAG,  
#endif
     &          METIS_OPTIONS(1), OPT_METIS_SIZE,
     &          IKEEP(1:N,2), IKEEP(1:N,1), INFO(1), LP, LPOK)
         ELSE IF (METIS_IDX_SIZE .EQ. 64) THEN
           CALL MUMPS_METIS_NODEND_MIXEDto64(N, IPE8, IW2,
#if defined(metis4) || defined(parmetis3)
     &          NUMFLAG, 
#else
     &          NUMFLAG, 
#endif
     &          METIS_OPTIONS(1), OPT_METIS_SIZE,
     &          IKEEP(1:N,2), IKEEP(1:N,1), INFO(1), 
     &          LP, LPOK, KEEP(10), 
     &          LLIW8, .FALSE., .TRUE. )
         ELSE
            WRITE(*,*)
     &        "Internal error in METIS wrappers, METIS_IDX_SIZE=",
     &        METIS_IDX_SIZE
            CALL MUMPS_ABORT()
         ENDIF
         IF (INFO(1) .LT. 0) GOTO 90
         DEALLOCATE(IW2)
      ELSE IF (IORD.NE.1) THEN
          WRITE(*,*) IORD
          WRITE(*,*)  'bad option for ordering'
          CALL MUMPS_ABORT()
        ENDIF
#endif
       DO K=1,N
         IW(L1+K) = 0
       ENDDO
       DO K=1,N
         IF ((IKEEP(K,1).LE.0).OR.(IKEEP(K,1).GT.N)) 
     &    GO TO 40
         IF (IW(L1+IKEEP(K,1)).EQ.1) THEN
          GOTO 40
         ELSE
          IW(L1+IKEEP(K,1)) = 1
         ENDIF
       ENDDO
       CALL DMUMPS_ANA_J1_ELT(N, NZ8, NELT, ELTPTR(NELT+1)-1,
     &             ELTPTR, ELTVAR, XNODEL, NODEL, 
     &             IKEEP, IWtemp(1,2), IW(L1))
       LLIW8 = NZ8+int(N,8)
       ALLOCATE( IW2(LLIW8), stat = allocok )
       IF (allocok.GT.0) THEN
         INFO(1) = -7
         CALL MUMPS_SET_IERROR(LLIW8,INFO(2))
         GOTO 90 
       ENDIF
       CALL DMUMPS_ANA_J2_ELT(N, NELT, ELTPTR(NELT+1)-1,
     &             ELTPTR, ELTVAR, XNODEL, NODEL, 
     &             IKEEP, IW2, LLIW8, IPE8, IWtemp(1,2),
     &             IW(L1), IWFR8)
       IF (KEEP(60) == 0) THEN
         ITEMP = 0 
       ELSE
         ITEMP = SIZE_SCHUR
         IF (KEEP(60) == 1) THEN
           KEEP(20) = LISTVAR_SCHUR(1)
         ELSEIF (KEEP(60) == 2 .OR. KEEP(60) == 3 ) THEN
           KEEP(38) = LISTVAR_SCHUR(1)
        ELSE
           WRITE(*,*) "Internal error in DMUMPS_ANA_F_ELT",KEEP(60)
           CALL MUMPS_ABORT()
         ENDIF
       ENDIF
       CALL DMUMPS_ANA_K(N, IPE8, IW2, LLIW8, IWFR8, IKEEP, 
     &    IKEEP(1,2), IW(L1),
     &    IW(L2), NCMPA, ITEMP, IWtemp)
      ENDIF
#if defined(OLDDFS)
      CALL DMUMPS_ANA_L(N, IWtemp, IW(L1), IKEEP, IKEEP(1,2),
     &     IKEEP(1,3),
     &     NFSIZ, INFO(6), FILS, FRERE, IWtemp(1,3), NEMIN, KEEP(60))
#else
      CALL DMUMPS_ANA_LNEW(N, IWtemp, IW(L1), IKEEP, IKEEP(1,2),
     &     IKEEP(1,3),
     &     NFSIZ, IWtemp(1,2), 
     &     INFO(6), FILS, FRERE, IWtemp(1,3), NEMIN, 
     &     IW(L2),  KEEP(60), KEEP(20), KEEP(38),
     &     IW2,KEEP(104),IW(L2+N),KEEP(50),
     &     ICNTL(13), KEEP(37), NSLAVES, KEEP(250).EQ.1,
     &     .FALSE., IDUMMY, LIDUMMY)
#endif
      DEALLOCATE(IW2)
      IF (KEEP(60).NE.0) THEN
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
      CALL DMUMPS_ANA_M(IKEEP(1,2),
     &  IWtemp(1,3), INFO(6),
     &  INFO(5), KEEP(2),KEEP(50),
     &  KEEP8(101), KEEP(108),KEEP(5),
     &  KEEP(6), KEEP(226), KEEP(253))
      IF ( KEEP(53) .NE. 0 ) THEN
        CALL MUMPS_MAKE1ROOT( N, FRERE, FILS, NFSIZ, KEEP(20) )
      END IF
      IF ( KEEP(48) == 4 .OR.
     &   ( (KEEP(24).NE.0).AND.(KEEP8(21).GT.0_8) ) ) THEN
          CALL DMUMPS_SET_K821_SURFACE(KEEP8(21), KEEP(2),
     &    KEEP(48), KEEP(50), NSLAVES)
      END IF
      IF (KEEP(210).LT.0.OR.KEEP(210).GT.2) KEEP(210)=0
      IF (KEEP(210).EQ.0.AND.KEEP(201).GT.0) KEEP(210)=1 
      IF (KEEP(210).EQ.0.AND.KEEP(201).EQ.0) KEEP(210)=2 
      IF (KEEP(210).EQ.2) KEEP8(79)=huge(KEEP8(79))
      IF (KEEP(210).EQ.1.AND.KEEP8(79).LE.0_8) THEN
        KEEP8(79)=K79REF * int(NSLAVES,8)
      ENDIF
      IF (KEEP(79).EQ.0) THEN
       IF (KEEP(210).EQ.1) THEN
        SPLITROOT = .FALSE. 
        IF ( KEEP(62).GE.1) THEN
          IDUMMY(1)= -1
          CALL DMUMPS_CUTNODES(N, FRERE, FILS, NFSIZ,
     &                       IDUMMY, LIDUMMY, INFO(6),
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
         IDUMMY(1) = -1
         CALL DMUMPS_CUTNODES(N, FRERE, FILS, NFSIZ,
     &                    IDUMMY, LIDUMMY, INFO(6),
     &                    NSLAVES, KEEP,KEEP8, SPLITROOT,
     &                    MP, LDIAG, INFO(1), INFO(2))
         IF (INFO(1).LT.0) GOTO 90 
         IF ( KEEP(53) .NE. 0 ) THEN
          CALL MUMPS_MAKE1ROOT( N, FRERE, FILS, NFSIZ, KEEP(20) )
         ENDIF
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
         K = min0(10,N)
         IF (LDIAG.EQ.4) K = N
         IF (K.GT.0) WRITE (MP,99997) (IKEEP(I,1),I=1,K)
         IF (K.GT.0) WRITE (MP,99991) (IKEEP(I,2),I=1,K)
         IF (K.GT.0) WRITE (MP,99990) (IKEEP(I,3),I=1,K)
         IF (K.GT.0) WRITE (MP,99987) (NFSIZ(I),I=1,K)
         IF (K.GT.0) WRITE (MP,99989) (FILS(I),I=1,K)
         IF (K.GT.0) WRITE (MP,99988) (FRERE(I),I=1,K)
      ENDIF
      GO TO 90
   40 INFO(1) = -4
      INFO(2) = K
   90 CONTINUE
      IF (INFO(1) .LT.0) THEN
        IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99996) INFO(1)
        IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99982) INFO(2)
      ENDIF
      IF (allocated(IW)) DEALLOCATE(IW)
      IF (allocated(IPE8)) DEALLOCATE(IPE8)
      IF (allocated(IW2)) DEALLOCATE(IW2)
      IF (allocated(IWtemp)) DEALLOCATE(IWtemp)
      RETURN
99999 FORMAT (/'Entering analysis phase with ...'/
     & '                N         NELT       LIW       INFO(1)'/,
     & 9X, I10, I11, I12, I14)
99998 FORMAT ('Element pointers:  ELTPTR()   '/(9X, 7I10))
99995 FORMAT ('Element variables: ELTVAR()   '/(9X, 7I10))
99997 FORMAT ('IKEEP(.,1)=', 10I6/(12X, 10I6))
99996 FORMAT (/'** Error return ** from Analysis   *  INFO(1)=', I3)
99991 FORMAT ('IKEEP(.,2)=', 10I6/(12X, 10I6))
99990 FORMAT ('IKEEP(.,3)=', 10I6/(12X, 10I6))
99989 FORMAT ('FILS (.)  =', 10I6/(12X, 10I6))
99988 FORMAT ('FRERE(.)  =', 10I6/(12X, 10I6))
99987 FORMAT ('NFSIZ(.)  =', 10I6/(12X, 10I6))
99982 FORMAT ('Error in permutation array KEEP   INFO(2)=', I3)
      END SUBROUTINE DMUMPS_ANA_F_ELT
      SUBROUTINE DMUMPS_NODEL( NELT, N, NELNOD, XELNOD, ELNOD,
     &                        XNODEL, NODEL, FLAG, IERROR, ICNTL ) 
      IMPLICIT NONE
      INTEGER NELT, N, NELNOD, IERROR, ICNTL(60)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     &        FLAG(N)
      INTEGER I, J, K, MP, NBERR
      MP = ICNTL(2)
      FLAG(1:N) = 0
      XNODEL(1:N) = 0
      IERROR = 0
      DO I = 1, NELT
        DO K = XELNOD(I), XELNOD(I+1)-1
          J = ELNOD(K)
          IF ( J.LT.1 .OR. J.GT.N ) THEN
            IERROR = IERROR + 1
          ELSE
            IF ( FLAG(J).NE.I ) THEN
              XNODEL(J) = XNODEL(J) + 1
              FLAG(J) = I
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF ( IERROR.GT.0 .AND. MP.GT.0 .AND. ICNTL(4).GE.2 ) THEN
        NBERR = 0
        WRITE(MP,99999)
        DO I = 1, NELT
          DO K = XELNOD(I), XELNOD(I+1)-1
            J = ELNOD(K)
            IF ( J.LT.1 .OR. J.GT.N ) THEN
              NBERR = NBERR + 1
              IF (NBERR.LE.10) THEN
                WRITE(MP,'(A,I8,A,I8,A)')
     &          'Element ',I,' variable ',J,' ignored.' 
              ELSE
                GO TO 100
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
  100 CONTINUE
      K = 1
      DO I = 1, N
         K = K + XNODEL(I)
         XNODEL(I) = K
      ENDDO
      XNODEL(N+1) = XNODEL(N)
      FLAG(1:N) = 0
      DO I = 1, NELT
         DO K = XELNOD(I), XELNOD(I+1)-1
            J = ELNOD(K)
            IF (FLAG(J).NE.I) THEN   
              XNODEL(J) = XNODEL(J) - 1
              NODEL(XNODEL(J)) = I
              FLAG(J) = I
            ENDIF
         ENDDO
      ENDDO
      RETURN
99999 FORMAT (/'*** Warning message from subroutine DMUMPS_NODEL ***')
      END SUBROUTINE DMUMPS_NODEL
      SUBROUTINE DMUMPS_ANA_G1_ELT(N, NZ, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  LEN, FLAG)
      IMPLICIT NONE
      INTEGER N, NELT, NELNOD
      INTEGER(8), INTENT(OUT) :: NZ
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER  XNODEL(N+1), NODEL(NELNOD),
     &        FLAG(N)
      INTEGER I,J,K1,K2,K3
      FLAG(1:N) = 0
      LEN(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN 
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                LEN(I) = LEN(I) + 1
                LEN(J) = LEN(J) + 1
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZ = 0_8
      DO I = 1,N
        NZ = NZ + int(LEN(I),8)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ANA_G1_ELT
      SUBROUTINE DMUMPS_ANA_G2_ELTNEW(N, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NELT,NELNOD
      INTEGER(8), INTENT(IN) :: LW
      INTEGER(8), INTENT(OUT) :: IWFR 
      INTEGER(8), INTENT(OUT) :: IPE(N+1)
      INTEGER LEN(N)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER  XNODEL(N+1), NODEL(NELNOD), 
     &          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1_8
      DO I = 1,N
        IWFR = IWFR + int(LEN(I),8)
          IPE(I) = IWFR
      ENDDO 
      IPE(N+1)=IPE(N)
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                IPE(I) = IPE(I) - 1
                IW(IPE(I)) = J
                IPE(J) = IPE(J) - 1
                IW(IPE(J)) = I
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ANA_G2_ELTNEW
      SUBROUTINE DMUMPS_ANA_G2_ELT(N, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NELT,NELNOD
      INTEGER(8), INTENT(IN) ::  LW
      INTEGER(8), INTENT(OUT) :: IWFR 
      INTEGER(8), INTENT(OUT) :: IPE(N)
      INTEGER LEN(N)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER  XNODEL(N+1), NODEL(NELNOD), 
     &          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1_8
      DO I = 1,N
        IWFR = IWFR + int(LEN(I),8)
        IF (LEN(I).GT.0) THEN
          IPE(I) = IWFR
        ELSE
          IPE(I) = 0_8
        ENDIF
      ENDDO 
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                IPE(I) = IPE(I) - 1_8
                IW(IPE(I)) = J
                IPE(J) = IPE(J) - 1_8
                IW(IPE(J)) = I
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ANA_G2_ELT
      SUBROUTINE DMUMPS_ANA_J1_ELT(N, NZ, NELT, NELNOD,
     & XELNOD, ELNOD, XNODEL, NODEL, 
     & PERM, LEN, FLAG)
      IMPLICIT NONE
      INTEGER N,NELT,NELNOD
      INTEGER(8), INTENT(OUT) :: NZ
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER PERM(N)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD), FLAG(N)
      INTEGER I,J,K1,K2,K3
      FLAG(1:N) = 0
      LEN(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                IF (PERM(J).GT.PERM(I)) THEN
                  LEN(I) = LEN(I) + 1
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZ = 0_8
      DO I = 1,N 
        NZ = NZ + int(LEN(I),8)
      ENDDO   
      RETURN  
      END SUBROUTINE DMUMPS_ANA_J1_ELT
      SUBROUTINE DMUMPS_ANA_J2_ELT(N, NELT, NELNOD,
     & XELNOD, ELNOD, XNODEL, NODEL, 
     & PERM, IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NELT,NELNOD
      INTEGER(8), INTENT(IN) :: LW
      INTEGER(8), INTENT(OUT) :: IWFR 
      INTEGER(8), INTENT(OUT) :: IPE(N)
      INTEGER  XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER PERM(N)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD), IW(LW), 
     &          FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 0_8
      DO I = 1,N
        IWFR = IWFR + int(LEN(I) + 1,8)
        IPE(I) = IWFR 
      ENDDO
      IWFR = IWFR + 1_8
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                IF (PERM(J).GT.PERM(I)) THEN
                  IW(IPE(I)) = J
                  IPE(I) = IPE(I) - 1_8
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DO I = 1,N
        J = int(IPE(I))
        IW(J) = LEN(I)
        IF (LEN(I).EQ.0) IPE(I) = 0_8
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ANA_J2_ELT
      SUBROUTINE DMUMPS_ANA_DIST_ELEMENTS( MYID, SLAVEF, N, 
     &           PROCNODE, STEP, PTRAIW, PTRARW,
     &           NELT, FRTPTR, FRTELT, 
     &           KEEP,KEEP8, ICNTL, SYM )
      IMPLICIT NONE
      INTEGER MYID, SLAVEF, N, NELT, SYM
      INTEGER KEEP( 500 ), ICNTL( 60 )
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: PTRAIW( NELT+1 ), PTRARW( NELT+1 )
      INTEGER STEP( N )
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PROCNODE( KEEP(28) )
      INTEGER MUMPS_TYPENODE, MUMPS_PROCNODE
      EXTERNAL MUMPS_TYPENODE, MUMPS_PROCNODE
      INTEGER(8) :: IPTRI8, IPTRR8, NVAR8
      INTEGER ELT, I, K
      INTEGER TYPE_PARALL, ITYPE, IRANK
      LOGICAL :: EARLYT3ROOTINS
      TYPE_PARALL = KEEP(46)
      PTRAIW( 1:NELT ) = 0_8
      EARLYT3ROOTINS = KEEP(200) .EQ.0
      DO I = 1, N
        IF (STEP(I).LT.0) CYCLE
        ITYPE = MUMPS_TYPENODE( PROCNODE(abs(STEP(I))), KEEP(199) )
        IRANK = MUMPS_PROCNODE( PROCNODE(abs(STEP(I))), KEEP(199) )
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK = IRANK + 1
        END IF
        IF ( (ITYPE .EQ. 2) .OR.
     &       (ITYPE .EQ. 3 .AND. .NOT. EARLYT3ROOTINS ) .OR.
     &       (ITYPE .EQ. 1 .AND. IRANK .EQ. MYID) ) THEN
          DO K = FRTPTR(I),FRTPTR(I+1)-1
            ELT = FRTELT(K)
            PTRAIW( ELT ) = PTRARW(ELT+1) - PTRARW(ELT)
          ENDDO
        ELSE  
        END IF
      END DO
      IPTRI8 = 1_8
      DO ELT = 1,NELT
        NVAR8 = PTRAIW( ELT )
        PTRAIW( ELT ) = IPTRI8
        IPTRI8 = IPTRI8 + NVAR8
      ENDDO
      PTRAIW( NELT+1 ) = IPTRI8
      KEEP8(27) = IPTRI8 - 1
      IF ( .TRUE. ) THEN  
        IF (SYM .EQ. 0) THEN
          IPTRR8 = 1_8
          DO ELT = 1,NELT
            NVAR8 = PTRAIW( ELT+1 ) - PTRAIW( ELT )
            PTRARW( ELT ) = IPTRR8
            IPTRR8 = IPTRR8 + NVAR8*NVAR8
          ENDDO
          PTRARW( NELT+1 ) = IPTRR8
        ELSE
          IPTRR8 = 1_8
          DO ELT = 1,NELT
            NVAR8 = PTRAIW( ELT+1 ) - PTRAIW( ELT )
            PTRARW( ELT ) = IPTRR8
            IPTRR8 = IPTRR8 + (NVAR8*(NVAR8+1))/2
          ENDDO
          PTRARW( NELT+1 ) = IPTRR8
        ENDIF
      ELSE
        IF (SYM .EQ. 0) THEN
          IPTRR8 = 1_8
          DO ELT = 1,NELT
            NVAR8 = PTRARW( ELT+1 ) - PTRARW( ELT )
            PTRARW( ELT ) = IPTRR8
            IPTRR8 = IPTRR8 + NVAR8*NVAR8
          ENDDO
          PTRARW( NELT+1 ) = IPTRR8
        ELSE
          IPTRR8 = 1_8
          DO ELT = 1,NELT
            NVAR8 = PTRARW( ELT+1 ) - PTRARW( ELT )
            PTRARW( ELT ) = IPTRR8
            IPTRR8 = IPTRR8 + (NVAR8*(NVAR8+1))/2  
          ENDDO 
          PTRARW( NELT+1 ) = IPTRR8
        ENDIF
      ENDIF 
      KEEP8(26) = IPTRR8 - 1_8
      RETURN
      END SUBROUTINE DMUMPS_ANA_DIST_ELEMENTS
      SUBROUTINE DMUMPS_ELTPROC( N, NELT, ELTPROC, SLAVEF, PROCNODE,
     &                           KEEP)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, NELT, SLAVEF
      INTEGER, INTENT(IN)    :: PROCNODE( N )
      INTEGER, INTENT(INOUT) :: ELTPROC( NELT )
      INTEGER                :: KEEP(500)
      INTEGER ELT, I, ITYPE
      LOGICAL :: EARLYT3ROOTINS 
      INTEGER, EXTERNAL :: MUMPS_TYPENODE, MUMPS_PROCNODE
      EARLYT3ROOTINS = KEEP(200) .EQ.0
      DO ELT = 1, NELT
          I = ELTPROC(ELT)
          IF ( I .NE. 0) THEN
           ITYPE = MUMPS_TYPENODE(PROCNODE(I),KEEP(199))
           IF (ITYPE.EQ.1) THEN
             ELTPROC(ELT) = MUMPS_PROCNODE(PROCNODE(I),KEEP(199))
           ELSE IF ( ITYPE.EQ.2 .OR. .NOT. EARLYT3ROOTINS ) THEN
             ELTPROC(ELT) = -1
           ELSE
            ELTPROC(ELT) = -2
           ENDIF
          ELSE
           ELTPROC(ELT) = -3
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ELTPROC
      SUBROUTINE DMUMPS_FRTELT(N, NELT, NELNOD, FRERE, FILS, NA, NE,
     &           XNODEL, NODEL, FRTPTR, FRTELT, ELTNOD) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, NELT, NELNOD
      INTEGER, INTENT(IN) :: FRERE(N), FILS(N), NA(N), NE(N)
      INTEGER, INTENT(OUT):: FRTPTR(N+1), FRTELT(NELT), ELTNOD(NELT)
      INTEGER, INTENT(IN) :: XNODEL(N+1), NODEL(NELNOD) 
      INTEGER, DIMENSION(:), ALLOCATABLE :: TNSTK, IPOOL
      INTEGER I, K, IFATH, allocok
      INTEGER INODE, LEAF, NBLEAF, NBROOT, III, IN
      ALLOCATE(TNSTK( N ), stat=allocok)
      IF (allocok.ne.0) THEN
         WRITE(6,*) ' Allocation error of TNSTK in '
     &       // 'routine DMUMPS_FRTELT '
         CALL MUMPS_ABORT()
      ENDIF
      ALLOCATE(IPOOL( N ), stat=allocok)
      IF (allocok.ne.0) THEN
         WRITE(6,*) ' Allocation error of IPOOL in '
     &      // 'routine DMUMPS_FRTELT '
         CALL MUMPS_ABORT()
      ENDIF
      TNSTK = NE
      LEAF = 1
      IF (N.EQ.1) THEN
        NBROOT = 1
        NBLEAF = 1
        IPOOL(1) = 1
        LEAF = LEAF + 1
      ELSEIF (NA(N).LT.0) THEN
        NBLEAF = N
        NBROOT = N
        DO 20 I=1,NBLEAF-1
           INODE = NA(I)
           IPOOL(LEAF) = INODE
           LEAF        = LEAF + 1
 20     CONTINUE
        INODE = -NA(N)-1
        IPOOL(LEAF) = INODE
        LEAF        = LEAF + 1
      ELSEIF (NA(N-1).LT.0) THEN
        NBLEAF = N-1
        NBROOT = NA(N)
        IF (NBLEAF-1.GT.0) THEN
         DO 30 I=1,NBLEAF-1
          INODE = NA(I)
          IPOOL(LEAF) = INODE
          LEAF        = LEAF + 1
 30      CONTINUE
        ENDIF
        INODE = -NA(N-1)-1
        IPOOL(LEAF) = INODE
        LEAF        = LEAF + 1
      ELSE
        NBLEAF = NA(N-1)
        NBROOT = NA(N)
        DO 40 I = 1,NBLEAF
          INODE = NA(I)
          IPOOL(LEAF) = INODE
          LEAF        = LEAF + 1
 40     CONTINUE
      ENDIF
      ELTNOD(1:NELT) = 0
      III = 1
 90   CONTINUE
        IF (III.NE.LEAF) THEN
           INODE=IPOOL(III)
           III = III + 1
        ELSE
           WRITE(6,*) ' ERROR 1 in subroutine DMUMPS_FRTELT '
           CALL MUMPS_ABORT()
        ENDIF
 95     CONTINUE 
        IN = INODE
 100    CONTINUE
        DO K = XNODEL(IN),XNODEL(IN+1)-1
          I = NODEL(K)
          IF (ELTNOD(I).EQ.0) ELTNOD(I) = INODE
        ENDDO
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
        IN = INODE
 110    IN = FRERE(IN)
        IF (IN.GT.0) GO TO 110
        IF (IN.EQ.0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 115
         GOTO 90
        ELSE
         IFATH = -IN
        ENDIF
        TNSTK(IFATH) = TNSTK(IFATH) - 1
        IF ( TNSTK(IFATH) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
        ELSE
            GOTO 90
        ENDIF
  115 CONTINUE
      FRTPTR(1:N) = 0
      DO I = 1,NELT
        IF (ELTNOD(I) .NE. 0) THEN
         FRTPTR(ELTNOD(I)) = FRTPTR(ELTNOD(I)) + 1
        ENDIF
      ENDDO
      K = 1
      DO I = 1,N
        K = K + FRTPTR(I)
        FRTPTR(I) = K
      ENDDO
      FRTPTR(N+1) = FRTPTR(N)
      DO K = 1,NELT
        INODE = ELTNOD(K)
        IF (INODE .NE. 0) THEN
         FRTPTR(INODE) = FRTPTR(INODE) - 1
         FRTELT(FRTPTR(INODE)) = K
        ENDIF
      ENDDO
      DEALLOCATE(TNSTK, IPOOL)
      RETURN
      END SUBROUTINE DMUMPS_FRTELT
      SUBROUTINE DMUMPS_ANA_G11_ELT(N, NZ, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  LEN, LW, IW)
      IMPLICIT NONE
      INTEGER N,NELT,NELNOD,LW
      INTEGER(8), INTENT(OUT) :: NZ
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     &        IW(LW)
      INTEGER I,J,K1,K2,K3,LP,NSUP,SUPVAR
      INTEGER INFO44(6)
      EXTERNAL DMUMPS_SUPVAR
      LP = 6
      CALL DMUMPS_SUPVAR(N,NELT,XELNOD(NELT+1)-1,ELNOD,XELNOD,
     &           NSUP,IW(3*N+3+1),3*N+3,IW,LP,INFO44)
      IF (INFO44(1) .LT. 0) THEN
        IF (LP.GE.0) WRITE(LP,*) 
     &     'Error return from DMUMPS_SUPVAR. INFO(1) = ',INFO44(1) 
      ENDIF
      IW(1:NSUP) = 0
      LEN(1:N) = 0
      DO I = 1,N
        SUPVAR = IW(3*N+3+1+I)
        IF (SUPVAR .EQ. 0) CYCLE
        IF (IW(SUPVAR).NE.0) THEN
          LEN(I) = -IW(SUPVAR)
        ELSE
          IW(SUPVAR) = I
        ENDIF
      ENDDO
      IW(N+1:2*N) = 0
      NZ = 0_8
      DO SUPVAR = 1,NSUP
        I = IW(SUPVAR)
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF (LEN(J).GE.0) THEN
                IF ((I.NE.J) .AND. (IW(N+J).NE.I)) THEN
                  IW(N+J) = I
                  LEN(I) = LEN(I) + 1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        NZ = NZ + int(LEN(I),8)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ANA_G11_ELT
      SUBROUTINE DMUMPS_ANA_G12_ELT(N, NELT, NELNOD,
     &  XELNOD, ELNOD, XNODEL, NODEL, 
     &  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NELT,NELNOD
      INTEGER(8), INTENT(IN) :: LW
      INTEGER(8), INTENT(OUT) :: IWFR 
      INTEGER(8), INTENT(OUT) :: IPE(N)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     &          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1_8
      DO I = 1,N
        IF (LEN(I).GT.0) THEN
          IWFR = IWFR + int(LEN(I),8)
          IPE(I) = IWFR
        ELSE
          IPE(I) = 0_8
        ENDIF
      ENDDO 
      FLAG(1:N) = 0
      DO I = 1,N
        IF (LEN(I).LE.0) CYCLE
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF (LEN(J) .GT. 0) THEN
                IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                  IPE(I) = IPE(I) - 1
                  IW(IPE(I)) = J
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ANA_G12_ELT
      SUBROUTINE DMUMPS_SUPVAR(N,NELT,NZ,ELTVAR,ELTPTR,NSUP,SVAR,
     &                 LIW,IW,LP,INFO)
      INTEGER LIW,LP,N,NELT,NSUP,NZ
      INTEGER INFO(6)
      INTEGER ELTPTR(NELT+1),ELTVAR(NZ)
      INTEGER IW(LIW),SVAR(0:N)
      INTEGER FLAG,NEW,VARS
      EXTERNAL DMUMPS_SUPVARB
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      INFO(4) = 0
      IF (N.LT.1) GO TO 10
      IF (NELT.LT.1) GO TO 20
      IF (NZ.LT.ELTPTR(NELT+1)-1) GO TO 30
      IF (LIW.LT.6) THEN
         INFO(4) = 3*N + 3
         GO TO 40
      END IF
      NEW = 1
      VARS = NEW + LIW/3
      FLAG = VARS + LIW/3
      CALL DMUMPS_SUPVARB(N,NELT,ELTPTR,NZ,ELTVAR,SVAR,NSUP,LIW/3-1,
     &           IW(NEW),IW(VARS),IW(FLAG),INFO)
      IF (INFO(1).EQ.-4) THEN
         INFO(4) = 3*N + 3
         GO TO 40
      ELSE
         INFO(4) = 3*NSUP + 3
      END IF
      GO TO 50
   10 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   20 INFO(1) = -2
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   30 INFO(1) = -3
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   40 INFO(1) = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9010) INFO(4)
      END IF
   50 RETURN
 9000 FORMAT (/3X,'Error message from DMUMPS_SUPVAR: INFO(1) = ',I2)
 9010 FORMAT (3X,'LIW is insufficient. Upper bound on required work',
     &       'space is ',I8)
      END SUBROUTINE DMUMPS_SUPVAR
      SUBROUTINE DMUMPS_SUPVARB( N, NELT, ELTPTR, NZ, ELTVAR,
     &           SVAR, NSUP, MAXSUP, NEW, VARS, FLAG, INFO )
      INTEGER MAXSUP,N,NELT,NSUP,NZ
      INTEGER ELTPTR(NELT+1),ELTVAR(NZ)
      INTEGER INFO(6)
      INTEGER FLAG(0:MAXSUP), NEW(0:MAXSUP),SVAR(0:N),
     &          VARS(0:MAXSUP)
      INTEGER I,IS,J,JS,K,K1,K2
      DO 10 I = 0,N
         SVAR(I) = 0
   10 CONTINUE
      VARS(0) = N + 1
      NEW(0) = -1
      FLAG(0) = 0
      NSUP = 0
      DO 40 J = 1,NELT
         K1 = ELTPTR(J)
         K2 = ELTPTR(J+1) - 1
         DO 20 K = K1,K2
            I = ELTVAR(K)
            IF (I.LT.1 .OR. I.GT.N) THEN
               INFO(2) = INFO(2) + 1
               GO TO 20
            END IF
            IS = SVAR(I)
            IF (IS.LT.0) THEN
               ELTVAR(K) = 0
               INFO(3) = INFO(3) + 1
               GO TO 20
            END IF
            SVAR(I) = SVAR(I) - N - 2
            VARS(IS) = VARS(IS) - 1
   20    CONTINUE
         DO 30 K = K1,K2
            I = ELTVAR(K)
            IF (I.LT.1 .OR. I.GT.N) GO TO 30
            IS = SVAR(I) + N + 2
            IF (FLAG(IS).LT.J) THEN
               FLAG(IS) = J
               IF (VARS(IS).GT.0) THEN
                  NSUP = NSUP + 1
                  IF (NSUP.GT.MAXSUP) THEN
                     INFO(1) = -4
                     RETURN
                  END IF
                  VARS(NSUP) = 1
                  FLAG(NSUP) = J
                  NEW(IS) = NSUP
                  SVAR(I) = NSUP
               ELSE
                  VARS(IS) = 1
                  NEW(IS) = IS
                  SVAR(I) = IS
               END IF
            ELSE
               JS = NEW(IS)
               VARS(JS) = VARS(JS) + 1
               SVAR(I) = JS
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_SUPVARB
