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
      MODULE MUMPS_ANA_ORD_WRAPPERS
      IMPLICIT NONE
       CONTAINS
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
#if defined(metis4) || defined(parmetis3)
      SUBROUTINE MUMPS_METIS_NODEWND_MIXEDto32( NCMP, IPE8, IW, FRERE,
     & NUMFLAG,
     & OPTIONS_METIS, LOPTIONS_METIS, IKEEP2, IKEEP1, INFO, 
     & LP, LPOK )
      IMPLICIT NONE
      INTEGER    :: INFO(2), LOPTIONS_METIS
      INTEGER    :: NCMP, NUMFLAG, IKEEP1(:), IKEEP2(:), FRERE(:)
      INTEGER    :: OPTIONS_METIS(LOPTIONS_METIS), IW(:)
      INTEGER, INTENT(IN) :: LP
      LOGICAL, INTENT(IN) :: LPOK
      INTEGER(8) :: IPE8(:)
      INTEGER,ALLOCATABLE, DIMENSION(:) :: IPE
      INTEGER :: allocok
      IF (IPE8(NCMP+1) .GE. int(huge(IW),8)) THEN
        INFO(1) = -51
        CALL MUMPS_SET_IERROR(
     &   IPE8(NCMP+1), INFO(2))
        RETURN
      ENDIF
      ALLOCATE(IPE(NCMP+1), stat=allocok)
      IF (allocok > 0) THEN
        INFO(1)=-7
        INFO(2)=NCMP+1
        IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR memory allocation in METIS_NODEWND_MIXEDto32"
        RETURN
      ENDIF
      CALL MUMPS_ICOPY_64TO32(IPE8(1), NCMP+1, IPE)
      CALL METIS_NODEWND(NCMP, IPE, IW(1),FRERE(1),
     &           NUMFLAG, OPTIONS_METIS,
     &           IKEEP2(1), IKEEP1(1) )
      RETURN
      END SUBROUTINE MUMPS_METIS_NODEWND_MIXEDto32
      SUBROUTINE MUMPS_METIS_NODEND_MIXEDto32( NCMP, IPE8, IW, NUMFLAG,
     & OPTIONS_METIS, LOPTIONS_METIS, IKEEP2, IKEEP1, INFO, 
     & LP, LPOK)
      IMPLICIT NONE
      INTEGER    :: INFO(2), LOPTIONS_METIS
      INTEGER    :: NCMP, NUMFLAG, IKEEP1(:), IKEEP2(:), IW(:)
      INTEGER    :: OPTIONS_METIS(LOPTIONS_METIS)
      INTEGER(8) :: IPE8(:)
      INTEGER, INTENT(IN) :: LP
      LOGICAL, INTENT(IN) :: LPOK
      INTEGER,ALLOCATABLE, DIMENSION(:) :: IPE
      INTEGER :: allocok
      IF (IPE8(NCMP+1) .GE. int(huge(IW),8)) THEN
        INFO(1) = -51
        CALL MUMPS_SET_IERROR(
     &   IPE8(NCMP+1), INFO(2))
        RETURN
      ENDIF
      ALLOCATE(IPE(NCMP+1), stat=allocok)
      IF (allocok > 0) THEN
        INFO(1)=-7
        INFO(2)=NCMP+1
        IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR memory allocation in METIS_NODEND_MIXEDto32"
        RETURN
      ENDIF
      CALL MUMPS_ICOPY_64TO32(IPE8(1), NCMP+1, IPE)
      CALL METIS_NODEND(NCMP, IPE, IW(1),
     &           NUMFLAG, OPTIONS_METIS,
     &           IKEEP2(1), IKEEP1(1) )
      DEALLOCATE(IPE)
      RETURN
      END SUBROUTINE MUMPS_METIS_NODEND_MIXEDto32
#else
      SUBROUTINE MUMPS_METIS_NODEND_MIXEDto32( NCMP, IPE8, IW, FRERE,
     & OPTIONS_METIS, LOPTIONS_METIS, IKEEP2, IKEEP1, INFO, 
     & LP, LPOK )
      IMPLICIT NONE
      INTEGER    :: INFO(2), LOPTIONS_METIS
      INTEGER    :: NCMP, IKEEP1(:), IKEEP2(:), FRERE(:), IW(:)
      INTEGER    :: OPTIONS_METIS(LOPTIONS_METIS)
      INTEGER(8) :: IPE8(:)
      INTEGER, INTENT(IN) :: LP
      LOGICAL, INTENT(IN) :: LPOK
      INTEGER,ALLOCATABLE, DIMENSION(:) :: IPE
      INTEGER :: allocok
      IF (IPE8(NCMP+1) .GE. int(huge(IW),8)) THEN
        INFO(1) = -51
        CALL MUMPS_SET_IERROR(
     &   IPE8(NCMP+1), INFO(2))
        RETURN
      ENDIF
      ALLOCATE(IPE(NCMP+1), stat=allocok)
      IF (allocok > 0) THEN
        INFO(1)=-7
        INFO(2)=NCMP+1
        IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR memory allocation in METIS_NODEND_MIXEDto32"
        RETURN
      ENDIF
      CALL MUMPS_ICOPY_64TO32(IPE8(1), NCMP+1, IPE)
      CALL METIS_NODEND( NCMP, IPE, IW(1), FRERE(1),
     & OPTIONS_METIS, IKEEP2(1), IKEEP1(1))
      DEALLOCATE(IPE)
      RETURN
      END SUBROUTINE MUMPS_METIS_NODEND_MIXEDto32
#endif
#endif
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
#if defined(metis4) || defined(parmetis3)
      SUBROUTINE MUMPS_METIS_NODEWND_MIXEDto64( NCMP, IPE8, IW, FRERE,
     & NUMFLAG,
     & OPTIONS_METIS, LOPTIONS_METIS, IKEEP2, IKEEP1, INFO, 
     & LP, LPOK, KEEP10, INPLACE64_GRAPH_COPY )
      IMPLICIT NONE
      INTEGER    :: INFO(2), LOPTIONS_METIS
      INTEGER    :: NCMP, NUMFLAG, IKEEP1(:), IKEEP2(:), FRERE(:)
      INTEGER    :: OPTIONS_METIS(LOPTIONS_METIS), IW(:)
      INTEGER(8) :: IPE8(:)
      INTEGER, INTENT(IN) :: LP, KEEP10
      LOGICAL, INTENT(IN) :: LPOK
      LOGICAL, INTENT(IN) :: INPLACE64_GRAPH_COPY
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IW8, FRERE8,
     &                                         IKEEP18, IKEEP28
      INTEGER :: allocok
      IF (KEEP10.EQ.1) THEN
       CALL METIS_NODEWND(NCMP, IPE8(1), IW(1),FRERE,
     &           NUMFLAG, OPTIONS_METIS,
     &           IKEEP2(1), IKEEP1(1) )
      ELSE
       IF (INPLACE64_GRAPH_COPY)  THEN
          CALL MUMPS_ICOPY_32TO64_64C_IP(IW(1), IPE8(NCMP+1)-1_8)
       ELSE
         ALLOCATE(IW8(IPE8(NCMP+1)-1_8),
     &             stat=allocok)
         IF (allocok > 0) THEN
          INFO(1)=-7
          CALL MUMPS_SET_IERROR(
     &     int(KEEP10,8)* ( IPE8(NCMP+1)-1_8 )
     &             , INFO(2)
     &                      )
          IF (LPOK) WRITE(LP,'(A)') 
     &      "ERROR memory allocation in METIS_NODEWND_MIXEDto64"
          RETURN
         ENDIF
         CALL MUMPS_ICOPY_32TO64_64C(IW(1), IPE8(NCMP+1)-1_8, IW8   )
       ENDIF
       ALLOCATE(FRERE8(NCMP),
     &         IKEEP18(NCMP), IKEEP28(NCMP), stat=allocok)
       IF (allocok > 0) THEN
          INFO(1)=-7
          CALL MUMPS_SET_IERROR(
     &     int(KEEP10,8)* ( 3_8*int(NCMP,8) )
     &             , INFO(2)
     &                      )
          IF (LPOK) WRITE(LP,'(A)') 
     &      "ERROR memory allocation in METIS_NODEWND_MIXEDto64"
          RETURN
       ENDIF
       CALL MUMPS_ICOPY_32TO64    (FRERE, NCMP        , FRERE8)
       IF (INPLACE64_GRAPH_COPY)  THEN
         CALL METIS_NODEWND(NCMP, IPE8(1), IW(1),FRERE8,
     &           NUMFLAG, OPTIONS_METIS,
     &           IKEEP2(1), IKEEP1(1) )
       ELSE
         CALL METIS_NODEWND(NCMP, IPE8(1), IW8,FRERE8,
     &           NUMFLAG, OPTIONS_METIS,
     &           IKEEP2(1), IKEEP1(1) )
       ENDIF
       CALL MUMPS_ICOPY_64TO32(IKEEP18, NCMP, IKEEP1(1))
       CALL MUMPS_ICOPY_64TO32(IKEEP28, NCMP, IKEEP2(1))
       IF (INPLACE64_GRAPH_COPY)  THEN
        DEALLOCATE(FRERE8, IKEEP18, IKEEP28)
       ELSE
        DEALLOCATE(IW8, FRERE8, IKEEP18, IKEEP28)
       ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_METIS_NODEWND_MIXEDto64
      SUBROUTINE MUMPS_METIS_NODEND_MIXEDto64( NCMP, IPE8, IW, NUMFLAG,
     & OPTIONS_METIS, LOPTIONS_METIS, IKEEP2, IKEEP1, INFO,
     & LP, LPOK, KEEP10,
     & LIW8, INPLACE64_GRAPH_COPY, INPLACE64_RESTORE_GRAPH
     &   )
      IMPLICIT NONE
      INTEGER    :: INFO(2), LOPTIONS_METIS
      INTEGER    :: NCMP, NUMFLAG, IKEEP1(:), IKEEP2(:), IW(:)
      INTEGER    :: OPTIONS_METIS(LOPTIONS_METIS)
      INTEGER(8) :: IPE8(:)
      INTEGER, INTENT(IN) :: LP, KEEP10
      LOGICAL, INTENT(IN) :: LPOK
      INTEGER(8) :: LIW8
      LOGICAL, INTENT(IN) :: INPLACE64_GRAPH_COPY, 
     &                       INPLACE64_RESTORE_GRAPH
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IW8,
     &                                         IKEEP18, IKEEP28
      INTEGER :: allocok
      IF (KEEP10.EQ.1) THEN
       CALL METIS_NODEND(NCMP, IPE8(1), IW(1),
     &           NUMFLAG, OPTIONS_METIS,
     &           IKEEP2(1), IKEEP1(1) )
      ELSE
        IF (INPLACE64_GRAPH_COPY)  THEN
          CALL MUMPS_ICOPY_32TO64_64C_IP(IW(1), IPE8(NCMP+1)-1_8)
        ELSE
          ALLOCATE(IW8(IPE8(NCMP+1)-1_8), stat=allocok)
          IF (allocok > 0) THEN
            INFO(1)=-7
            CALL MUMPS_SET_IERROR(  int(KEEP10,8)*
     &         ( IPE8(NCMP+1)-1_8+2_8*int(NCMP,8) )
     &         , INFO(2) )
            IF (LPOK) WRITE(LP,'(A)') 
     &      "ERROR 1 memory allocation in METIS_METIS_NODEND_MIXEDto64"
            RETURN
          ENDIF
          CALL MUMPS_ICOPY_32TO64_64C(IW(1), IPE8(NCMP+1)-1_8, IW8 )
        ENDIF
        ALLOCATE(IKEEP18(NCMP), IKEEP28(NCMP), stat=allocok)
        IF (allocok > 0) THEN
            INFO(1)=-7
            CALL MUMPS_SET_IERROR(  int(KEEP10,8)*
     &         2_8*int(NCMP,8), INFO(2) )
            IF (LPOK) WRITE(LP,'(A)') 
     &      "ERROR 2 memory allocation in METIS_METIS_NODEND_MIXEDto64"
            RETURN
        ENDIF
        IF (INPLACE64_GRAPH_COPY)  THEN
          CALL METIS_NODEND(NCMP, IPE8(1), IW(1),
     &              NUMFLAG, OPTIONS_METIS,
     &              IKEEP28, IKEEP18 )
        ELSE
          CALL METIS_NODEND(NCMP, IPE8(1), IW8,
     &             NUMFLAG, OPTIONS_METIS,
     &             IKEEP28, IKEEP18 )
        ENDIF
        CALL MUMPS_ICOPY_64TO32(IKEEP18, NCMP, IKEEP1(1))
        CALL MUMPS_ICOPY_64TO32(IKEEP28, NCMP, IKEEP2(1))
        IF (INPLACE64_GRAPH_COPY) THEN
         IF (INPLACE64_RESTORE_GRAPH) THEN
          CALL MUMPS_ICOPY_64TO32_64C_IP(IW(1), IPE8(NCMP+1)-1_8)
         ENDIF
         DEALLOCATE(IKEEP18, IKEEP28)
        ELSE
          DEALLOCATE(IW8, IKEEP18, IKEEP28)
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_METIS_NODEND_MIXEDto64
#else
      SUBROUTINE MUMPS_METIS_NODEND_MIXEDto64( NCMP, IPE8, IW, FRERE,
     & OPTIONS_METIS, LOPTIONS_METIS, IKEEP2, IKEEP1, INFO,
     & LP, LPOK, KEEP10,
     & LIW8, INPLACE64_GRAPH_COPY, INPLACE64_RESTORE_GRAPH
     &  )
      IMPLICIT NONE
      INTEGER    :: INFO(2)
      INTEGER    :: LOPTIONS_METIS
      INTEGER    :: NCMP, IKEEP1(:), IKEEP2(:), FRERE(:), IW(:)
      INTEGER    :: OPTIONS_METIS(LOPTIONS_METIS)
      INTEGER(8) :: IPE8(:)
      INTEGER, INTENT(IN) :: LP, KEEP10
      LOGICAL, INTENT(IN) :: LPOK
      INTEGER(8) :: LIW8
      LOGICAL, INTENT(IN) :: INPLACE64_GRAPH_COPY,
     &                       INPLACE64_RESTORE_GRAPH
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IW8, FRERE8,
     &                                         IKEEP18, IKEEP28,
     &                                         OPTIONS_METIS8
      INTEGER :: allocok
      IF (KEEP10.EQ.1) THEN
       CALL METIS_NODEND( NCMP, IPE8(1), IW(1), FRERE(1),
     &      OPTIONS_METIS, IKEEP2(1), IKEEP1(1) )
      ELSE
        IF (INPLACE64_GRAPH_COPY)  THEN
          CALL MUMPS_ICOPY_32TO64_64C_IP(IW(1), IPE8(NCMP+1)-1_8)
        ELSE
          ALLOCATE(IW8(IPE8(NCMP+1)-1_8), stat=allocok)
          IF (allocok > 0) THEN
            INFO(1)=-7
            CALL MUMPS_SET_IERROR( int(KEEP10,8) * (IPE8(NCMP+1)-1_8)
     &         , INFO(2) )
            IF (LPOK) WRITE(LP,'(A)') 
     &      "ERROR 1 memory allocation in METIS_METIS_NODEND_MIXEDto64"
            RETURN
          ENDIF
          CALL MUMPS_ICOPY_32TO64_64C(IW(1), IPE8(NCMP+1)-1_8, IW8 )
        ENDIF
        ALLOCATE(FRERE8(NCMP),
     &         IKEEP18(NCMP), IKEEP28(NCMP),
     &         OPTIONS_METIS8(LOPTIONS_METIS), stat=allocok)
        IF (allocok > 0) THEN
          INFO(1)=-7
          CALL MUMPS_SET_IERROR( 
     &       int(KEEP10,8)*
     &           (3_8*int(NCMP,8)+int(LOPTIONS_METIS,8))
     &            , INFO(2))
          IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR 2 memory allocation in METIS_NODEND_MIXEDto64"
          RETURN
        ENDIF
        CALL MUMPS_ICOPY_32TO64    (FRERE(1), NCMP, FRERE8)
        CALL MUMPS_ICOPY_32TO64    (OPTIONS_METIS, LOPTIONS_METIS,
     &                                 OPTIONS_METIS8)
        IF (INPLACE64_GRAPH_COPY) THEN
          CALL METIS_NODEND( int(NCMP,8), IPE8(1), IW(1), FRERE8,
     &                       OPTIONS_METIS8, IKEEP28, IKEEP18 )
        ELSE
          CALL METIS_NODEND( int(NCMP,8), IPE8(1), IW8, FRERE8,
     &                       OPTIONS_METIS8, IKEEP28, IKEEP18 )
        ENDIF
        CALL MUMPS_ICOPY_64TO32(IKEEP18, NCMP, IKEEP1(1))
        CALL MUMPS_ICOPY_64TO32(IKEEP28, NCMP, IKEEP2(1))
        IF (INPLACE64_GRAPH_COPY) THEN
         IF (INPLACE64_RESTORE_GRAPH) THEN
            CALL MUMPS_ICOPY_64TO32_64C_IP(IW(1), IPE8(NCMP+1)-1_8)
         ENDIF
         DEALLOCATE(FRERE8, IKEEP18, IKEEP28, OPTIONS_METIS8)
        ELSE
          DEALLOCATE(IW8, FRERE8, IKEEP18, IKEEP28, OPTIONS_METIS8)
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_METIS_NODEND_MIXEDto64
#endif
#endif
#if defined(scotch) || defined(ptscotch)
      SUBROUTINE MUMPS_SCOTCH_MIXEDto32(NCMP, LIW8, IPE8, PARENT, IWFR8,
     &           PTRAR, IW, IWL1, IKEEP1,
     &           IKEEP2, NCMPA, INFO, LP, LPOK,
     &           WEIGHTUSED, WEIGHTREQUESTED)
      IMPLICIT NONE
      INTEGER,    INTENT(IN)    :: NCMP
      INTEGER(8), INTENT(IN)    :: LIW8
      INTEGER,    INTENT(OUT)   :: NCMPA
      INTEGER(8), INTENT(INOUT) :: IPE8(:) 
      INTEGER,    INTENT(OUT)   :: PARENT(NCMP)
      INTEGER(8), INTENT(IN)    :: IWFR8
      INTEGER                   :: PTRAR(NCMP)
      INTEGER                   :: IW(:) 
      INTEGER                   :: IWL1(NCMP)
      INTEGER,    INTENT(OUT)   :: IKEEP1(:)  
      INTEGER,    INTENT(OUT)   :: IKEEP2(:)  
      INTEGER,    INTENT(INOUT) :: INFO(2)
      INTEGER,    INTENT(IN)    :: LP
      LOGICAL,    INTENT(IN)    :: LPOK
      INTEGER,    INTENT(OUT)   :: WEIGHTUSED
      INTEGER,    INTENT(IN)    :: WEIGHTREQUESTED
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPE
      INTEGER :: allocok
      IF (IWFR8 .GE. int(huge(IW),8)) THEN
        INFO(1) = -51
        CALL MUMPS_SET_IERROR(IPE8(NCMP+1), INFO(2))
        RETURN
      ENDIF
      ALLOCATE(IPE(NCMP+1), stat=allocok)
      IF (allocok > 0) THEN
        IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR memory allocation in MUMPS_SCOTCH_MIXEDto32"
        INFO(1) = -7
        INFO(2) = NCMP+1
        RETURN
      ENDIF
      CALL MUMPS_ICOPY_64TO32(IPE8(1),NCMP+1,IPE)
      CALL MUMPS_SCOTCH( NCMP, int(LIW8), IPE, int(IWFR8),
     &                   PTRAR, IW(1), IWL1, IKEEP1(1),
     &                   IKEEP2(1), NCMPA,
     &                   WEIGHTUSED, WEIGHTREQUESTED )
      PARENT(1:NCMP)=IPE(1:NCMP)
      DEALLOCATE(IPE)
      RETURN
      END SUBROUTINE MUMPS_SCOTCH_MIXEDto32
      SUBROUTINE MUMPS_SCOTCH_MIXEDto64(
     &           NCMP, LIW8, IPE8, PARENT, IWFR8,
     &           PTRAR, IW, IWL1, IKEEP1,
     &           IKEEP2, NCMPA, INFO, LP, LPOK, KEEP10,
     &           INPLACE64_GRAPH_COPY,
     &           WEIGHTUSED, WEIGHTREQUESTED )
      IMPLICIT NONE
      INTEGER,    INTENT(IN)    :: NCMP
      INTEGER(8), INTENT(IN)    :: LIW8
      INTEGER,    INTENT(OUT)   :: NCMPA
      INTEGER(8), INTENT(INOUT) :: IPE8(:) 
      INTEGER,    INTENT(OUT)   :: PARENT(NCMP)
      INTEGER(8), INTENT(IN)    :: IWFR8
      INTEGER                   :: PTRAR(NCMP)
      INTEGER                   :: IW(:)
      INTEGER                   :: IWL1(NCMP)
      INTEGER,    INTENT(OUT)   :: IKEEP1(:) 
      INTEGER,    INTENT(OUT)   :: IKEEP2(:) 
      INTEGER,    INTENT(INOUT) :: INFO(2)
      INTEGER,    INTENT(IN)    :: LP
      LOGICAL,    INTENT(IN)    :: LPOK
      INTEGER,    INTENT(IN)    :: KEEP10
      LOGICAL, INTENT(IN)       :: INPLACE64_GRAPH_COPY
      INTEGER,    INTENT(OUT)   :: WEIGHTUSED
      INTEGER,    INTENT(IN)    :: WEIGHTREQUESTED
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: 
     &                     PTRAR8, IW8, IWL18, IKEEP18,
     &                     IKEEP28, IPE8_TEMP
      INTEGER :: allocok
      ALLOCATE( IPE8_TEMP(NCMP+1), stat=allocok )
      IF (allocok > 0) THEN
        IF (LPOK) WRITE(LP,*)
     &  "ERROR memory allocation in MUMPS_SCOTCH_MIXEDto64"
        INFO(1) = -7
        INFO(2) = NCMP+1
        RETURN
      ENDIF
      IPE8_TEMP(1:NCMP+1) = IPE8(1:NCMP+1)
      IF (KEEP10.EQ.1) THEN
        CALL MUMPS_SCOTCH_64( NCMP, LIW8,
     &                        IPE8_TEMP(1), 
     &                        IWFR8,
     &                        PTRAR, IW(1), IWL1, IKEEP1(1),
     &                        IKEEP2(1), NCMPA,
     &                        WEIGHTUSED, WEIGHTREQUESTED)
        PARENT(1:NCMP) = int(IPE8_TEMP(1:NCMP))
      ELSE
        IF (INPLACE64_GRAPH_COPY)  THEN
          CALL MUMPS_ICOPY_32TO64_64C_IP(IW(1), IPE8_TEMP(NCMP+1)-1_8)
        ELSE
         ALLOCATE( IW8(LIW8), stat=allocok )
         IF (allocok > 0) THEN
          IF (LPOK) WRITE(LP,*)
     &    "ERROR memory allocation in MUMPS_SCOTCH_MIXEDto64"
          INFO(1) = -7
          CALL MUMPS_SET_IERROR( int(KEEP10,8) * LIW8
     &                        , INFO(2) )
          GOTO 500
         ENDIF
         CALL MUMPS_ICOPY_32TO64_64C(IW(1),LIW8,IW8)
        ENDIF
        ALLOCATE( 
     &     PTRAR8(NCMP), IWL18(NCMP), IKEEP18(NCMP), IKEEP28(NCMP),
     &     stat=allocok )
        IF (allocok > 0) THEN
          IF (LPOK) WRITE(LP,*)
     &    "ERROR memory allocation in MUMPS_SCOTCH_MIXEDto64"
          INFO(1) = -7
          CALL MUMPS_SET_IERROR( int(KEEP10,8) *
     &                            ( int(NCMP,8)*4_8 )
     &                        , INFO(2) )
          GOTO 500
         ENDIF
        CALL MUMPS_ICOPY_32TO64(PTRAR,NCMP,PTRAR8)
        IF (WEIGHTREQUESTED.EQ.1) THEN
         CALL MUMPS_ICOPY_32TO64(IWL1,NCMP,IWL18)
        ENDIF
        IF (INPLACE64_GRAPH_COPY) THEN
         CALL MUMPS_SCOTCH_64( int(NCMP,8), LIW8,   
     &                        IPE8_TEMP(1),         
     &                        IWFR8,                
     &                        PTRAR8, IW(1), IWL18, 
     &                        IKEEP1(1), IKEEP2(1), NCMPA, 
     &                        WEIGHTUSED,           
     &                        WEIGHTREQUESTED )     
        ELSE
         CALL MUMPS_SCOTCH_64( int(NCMP,8), LIW8,   
     &                        IPE8_TEMP(1),         
     &                        IWFR8,                
     &                        PTRAR8, IW8, IWL18,   
     &                        IKEEP1(1), IKEEP2(1), NCMPA, 
     &                        WEIGHTUSED,           
     &                        WEIGHTREQUESTED )     
        ENDIF
        IF (NCMPA .LT. 0) THEN
           IF (LPOK) WRITE(LP,*) 
     &               ' Error on output from SCOTCH, NCMPA=', NCMPA
           INFO( 1 ) = -9999
           INFO( 2 ) = 3
           GOTO 500 
        ENDIF
        CALL MUMPS_ICOPY_64TO32(IWL18,NCMP,IWL1)
        CALL MUMPS_ICOPY_64TO32(IKEEP18,NCMP,IKEEP1(1))
        CALL MUMPS_ICOPY_64TO32(IKEEP28,NCMP,IKEEP2(1))
        CALL MUMPS_ICOPY_64TO32(IPE8_TEMP(1),NCMP,PARENT)
 500    CONTINUE
        IF (.NOT.INPLACE64_GRAPH_COPY) THEN
         IF (ALLOCATED(IW8)) DEALLOCATE(IW8)
        ENDIF
        IF (ALLOCATED(PTRAR8)) DEALLOCATE(PTRAR8)
        IF (ALLOCATED(IWL18)) DEALLOCATE(IWL18)
        IF (ALLOCATED(IKEEP18)) DEALLOCATE(IKEEP18)
        IF (ALLOCATED(IKEEP28)) DEALLOCATE(IKEEP28)
      ENDIF
      IF (ALLOCATED(IPE8_TEMP)) DEALLOCATE(IPE8_TEMP)
      RETURN
      END SUBROUTINE MUMPS_SCOTCH_MIXEDto64
#endif
#if defined (scotch) || defined (ptscotch)
      SUBROUTINE MUMPS_SCOTCH_KWAY_MIXEDto32(NHALO, HALOEDGENBR,
     &     IPTRHALO, JCNHALO,
     &     NBGROUPS, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
      IMPLICIT NONE
      include 'scotchf.h'
      INTEGER(8) :: HALOEDGENBR
      INTEGER    :: NHALO, NBGROUPS
      INTEGER    :: JCNHALO(HALOEDGENBR), PARTS(NHALO)
      INTEGER(8) :: IPTRHALO(NHALO+1)
      INTEGER, INTENT(IN)    :: LP, KEEP10
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      DOUBLE PRECISION :: GRAFDAT(SCOTCH_GRAPHDIM)
      DOUBLE PRECISION :: STRADAT(SCOTCH_STRATDIM)
      INTEGER :: BASEVAL, IERR, EDGENBR
      INTEGER, ALLOCATABLE    :: IPTRHALO_I4(:)
      INTEGER :: allocok
      IF (IPTRHALO(size(IPTRHALO)) .GE. int(huge(LP),8)) THEN
        IFLAG   = -51
        CALL MUMPS_SET_IERROR( IPTRHALO (size(IPTRHALO)), 
     &                         IERROR )
        RETURN
      ENDIF
      ALLOCATE(IPTRHALO_I4(size(IPTRHALO)), stat=allocok)
      IF (allocok > 0) THEN
        IFLAG   = -7
        IERROR  = size(IPTRHALO)
        IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR memory allocation in MUMPS_SCOTCH_KWAY_MIXEDto32"
        RETURN
      END IF
      CALL MUMPS_ICOPY_64TO32(IPTRHALO,
     &     size(IPTRHALO), IPTRHALO_I4)
      BASEVAL = 1
      EDGENBR = IPTRHALO_I4(NHALO+1)
      CALL SCOTCHFGRAPHBUILD(GRAFDAT(1), BASEVAL, NHALO,
     &     IPTRHALO_I4(1), IPTRHALO_I4(2), IPTRHALO_I4(1),
     &     IPTRHALO_I4(1), EDGENBR, JCNHALO(1), JCNHALO(1), IERR)
      CALL SCOTCHFSTRATINIT(STRADAT, IERR)
      CALL SCOTCHFGRAPHPART(GRAFDAT(1), NBGROUPS, STRADAT(1),
     &     PARTS(1), IERR)
      CALL SCOTCHFSTRATEXIT(STRADAT)
      CALL SCOTCHFGRAPHEXIT(GRAFDAT)
      PARTS(1:NHALO) = PARTS(1:NHALO)+1
      DEALLOCATE(IPTRHALO_I4)
      RETURN
      END SUBROUTINE MUMPS_SCOTCH_KWAY_MIXEDto32
      SUBROUTINE MUMPS_SCOTCH_KWAY_MIXEDto64(NHALO, HALOEDGENBR,
     &     IPTRHALO, JCNHALO,
     &     NBGROUPS, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
      IMPLICIT NONE
      include 'scotchf.h'
      INTEGER(8) :: HALOEDGENBR
      INTEGER    :: NHALO, NBGROUPS
      INTEGER    :: JCNHALO(HALOEDGENBR), PARTS(NHALO)
      INTEGER(8) :: IPTRHALO(NHALO+1)
      INTEGER, INTENT(IN)    :: LP, KEEP10
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      DOUBLE PRECISION :: GRAFDAT(SCOTCH_GRAPHDIM)
      DOUBLE PRECISION :: STRADAT(SCOTCH_STRATDIM)
      INTEGER :: IERR
      INTEGER(8), ALLOCATABLE :: JCNHALO_I8(:), PARTS_I8(:)
      INTEGER(8) :: NHALO_I8, NBGROUPS_I8, EDGENBR_I8,
     &     BASEVAL_I8
      INTEGER :: allocok
      ALLOCATE(JCNHALO_I8(IPTRHALO(NHALO+1)-1_8),
     &   PARTS_I8(size(PARTS)),  stat=allocok)
      IF (allocok > 0) THEN
        IFLAG  =-7
        CALL MUMPS_SET_IERROR(
     &   int(KEEP10,8)* (IPTRHALO(NHALO+1)-1_8
     &  +int(size(PARTS),8)),
     &   IERROR)
        IF (LPOK) WRITE(LP,'(A)') 
     &       "ERROR memory allocation in MUMPS_SCOTCH_KWAY_MIXEDto64 "
       ENDIF
      CALL MUMPS_ICOPY_32TO64_64C(JCNHALO,
     &     IPTRHALO(NHALO+1)-1, JCNHALO_I8)
      NHALO_I8    = int(NHALO,8)
      NBGROUPS_I8 = int(NBGROUPS,8)
      BASEVAL_I8  = 1_8
      EDGENBR_I8  = IPTRHALO(NHALO+1)
      CALL SCOTCHFGRAPHBUILD(GRAFDAT(1), BASEVAL_I8, NHALO_I8,
     &     IPTRHALO(1), IPTRHALO(2), IPTRHALO(1),
     &     IPTRHALO(1), EDGENBR_I8, JCNHALO_I8(1), JCNHALO_I8(1), IERR)
      CALL SCOTCHFSTRATINIT(STRADAT, IERR)
      CALL SCOTCHFGRAPHPART(GRAFDAT(1), NBGROUPS_I8, STRADAT(1),
     &     PARTS_I8(1), IERR)
      CALL SCOTCHFSTRATEXIT(STRADAT)
      CALL SCOTCHFGRAPHEXIT(GRAFDAT)
      CALL MUMPS_ICOPY_64TO32(PARTS_I8,
     &     size(PARTS), PARTS)
      DEALLOCATE(JCNHALO_I8, PARTS_I8)      
      PARTS(1:NHALO) = PARTS(1:NHALO)+1
      RETURN
      END SUBROUTINE MUMPS_SCOTCH_KWAY_MIXEDto64
#endif
#if defined (metis) || defined (parmetis) || defined (metis4) || defined (parmetis3)
      SUBROUTINE MUMPS_METIS_KWAY_MIXEDto32(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO, NBGROUPS, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
      IMPLICIT NONE
      INTEGER(8) :: HALOEDGENBR
      INTEGER    :: NHALO, NBGROUPS
      INTEGER    :: JCNHALO(HALOEDGENBR), PARTS(NHALO)
      INTEGER(8) :: IPTRHALO(NHALO+1)
      INTEGER, INTENT(IN)    :: LP, KEEP10
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, ALLOCATABLE    :: IPTRHALO_I4(:)
      INTEGER :: allocok
      IF (IPTRHALO(size(IPTRHALO)) .GE. int(huge(LP),8)) THEN
        IFLAG   = -51
        CALL MUMPS_SET_IERROR( IPTRHALO (size(IPTRHALO)), 
     &   IERROR)
        RETURN
      ENDIF
      ALLOCATE(IPTRHALO_I4(size(IPTRHALO)), stat=allocok)
      IF (allocok > 0) THEN
        IFLAG   = -7
        IERROR  = size(IPTRHALO)
        IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR memory allocation in MUMPS_METIS_KWAY_MIXEDto32"
        RETURN
      END IF
      CALL MUMPS_ICOPY_64TO32(IPTRHALO,
     &     size(IPTRHALO), IPTRHALO_I4)
      CALL MUMPS_METIS_KWAY(NHALO, IPTRHALO_I4(1),
     &               JCNHALO(1), NBGROUPS, PARTS(1))
      DEALLOCATE(IPTRHALO_I4)
      RETURN
      END SUBROUTINE MUMPS_METIS_KWAY_MIXEDto32
      SUBROUTINE MUMPS_METIS_KWAY_MIXEDto64(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO, NBGROUPS, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
      IMPLICIT NONE
      INTEGER(8) :: HALOEDGENBR
      INTEGER    :: NHALO, NBGROUPS
      INTEGER    :: JCNHALO(HALOEDGENBR), PARTS(NHALO)
      INTEGER(8) :: IPTRHALO(NHALO+1)
      INTEGER, INTENT(IN)    :: LP, KEEP10
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: JCNHALO_I8, PARTS_I8
      INTEGER(8) :: NHALO_I8, NBGROUPS_I8
      INTEGER :: allocok
      ALLOCATE(JCNHALO_I8(IPTRHALO(NHALO+1)-1_8),
     &         PARTS_I8(size(PARTS)), stat=allocok)
      IF (allocok > 0) THEN
        IFLAG  = -7
        CALL MUMPS_SET_IERROR(
     &   int(KEEP10,8)* (IPTRHALO(NHALO+1)-1_8+int(size(PARTS),8)),
     &   IERROR)
        IF (LPOK) WRITE(LP,'(A)') 
     &       "ERROR memory allocation in MUMPS_METIS_KWAY_MIXEDto64 "
       ENDIF
      NHALO_I8    = int(NHALO,8)
      NBGROUPS_I8 = int(NBGROUPS,8)
      CALL MUMPS_ICOPY_32TO64_64C(JCNHALO,
     &     IPTRHALO(NHALO+1)-1, JCNHALO_I8)
      CALL MUMPS_METIS_KWAY_64(NHALO_I8, IPTRHALO(1),
     &               JCNHALO_I8(1), NBGROUPS_I8, PARTS_I8(1))
      CALL MUMPS_ICOPY_64TO32(PARTS_I8,
     &     size(PARTS), PARTS)
      DEALLOCATE(JCNHALO_I8, PARTS_I8)
      RETURN
      END SUBROUTINE MUMPS_METIS_KWAY_MIXEDto64
      SUBROUTINE MUMPS_METIS_KWAY_AB_MIXEDto32(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO, NBGROUPS, PARTS, VWGT, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
      IMPLICIT NONE
      INTEGER(8) :: HALOEDGENBR
      INTEGER    :: NHALO, NBGROUPS
      INTEGER    :: JCNHALO(HALOEDGENBR), PARTS(NHALO), VWGT(NHALO)
      INTEGER(8) :: IPTRHALO(NHALO+1)
      INTEGER, INTENT(IN)    :: LP, KEEP10
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, ALLOCATABLE    :: IPTRHALO_I4(:)
      INTEGER :: allocok
      IF (IPTRHALO(size(IPTRHALO)) .GE. int(huge(LP),8)) THEN
        IFLAG   = -51
        CALL MUMPS_SET_IERROR( IPTRHALO (size(IPTRHALO)), 
     &   IERROR)
        RETURN
      ENDIF
      ALLOCATE(IPTRHALO_I4(size(IPTRHALO)), stat=allocok)
      IF (allocok > 0) THEN
        IFLAG   = -7
        IERROR  = size(IPTRHALO)
        IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR memory allocation in MUMPS_METIS_KWAY_AB_MIXEDto32"
        RETURN
      END IF
      CALL MUMPS_ICOPY_64TO32(IPTRHALO,
     &     size(IPTRHALO), IPTRHALO_I4)
      CALL MUMPS_METIS_KWAY_AB(NHALO, IPTRHALO_I4(1),
     &               JCNHALO(1), NBGROUPS, PARTS(1), VWGT(1))
      DEALLOCATE(IPTRHALO_I4)
      RETURN
      END SUBROUTINE MUMPS_METIS_KWAY_AB_MIXEDto32
      SUBROUTINE MUMPS_METIS_KWAY_AB_MIXEDto64(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO, NBGROUPS, PARTS, VWGT, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
      IMPLICIT NONE
      INTEGER(8) :: HALOEDGENBR
      INTEGER    :: NHALO, NBGROUPS
      INTEGER    :: JCNHALO(HALOEDGENBR), PARTS(NHALO), VWGT(NHALO)
      INTEGER(8) :: IPTRHALO(NHALO+1)
      INTEGER, INTENT(IN)    :: LP, KEEP10
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: JCNHALO_I8, PARTS_I8
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: VWGT_I8
      INTEGER(8) :: NHALO_I8, NBGROUPS_I8
      INTEGER :: allocok
      ALLOCATE(JCNHALO_I8(IPTRHALO(NHALO+1)-1_8),
     &         PARTS_I8(size(PARTS)), VWGT_I8(NHALO), stat=allocok)
      IF (allocok > 0) THEN
        IFLAG  = -7
        CALL MUMPS_SET_IERROR(
     &   int(KEEP10,8)* (IPTRHALO(NHALO+1)-1_8+int(size(PARTS),8))
     &   +int(NHALO,8), IERROR)
        IF (LPOK) WRITE(LP,'(A)') 
     &       "ERROR memory allocation in MUMPS_METIS_KWAY_AB_MIXEDto64 "
       ENDIF
      NHALO_I8    = int(NHALO,8)
      NBGROUPS_I8 = int(NBGROUPS,8)
      CALL MUMPS_ICOPY_32TO64_64C(JCNHALO,
     &     IPTRHALO(NHALO+1)-1, JCNHALO_I8)
      CALL MUMPS_ICOPY_32TO64_64C(VWGT,
     &     NHALO_I8, VWGT_I8)
      CALL MUMPS_METIS_KWAY_AB_64(NHALO_I8, IPTRHALO(1),
     &               JCNHALO_I8(1), NBGROUPS_I8, PARTS_I8(1), 
     &               VWGT_I8(1))
      CALL MUMPS_ICOPY_64TO32(PARTS_I8,
     &     size(PARTS), PARTS)
      DEALLOCATE(JCNHALO_I8, PARTS_I8, VWGT_I8)
      RETURN
      END SUBROUTINE MUMPS_METIS_KWAY_AB_MIXEDto64
#endif      
#if defined(pord)
      SUBROUTINE MUMPS_PORDF_MIXEDto32( NVTX, NEDGES8, XADJ8, IW,
     &                                  NV, NCMPA, PARENT,
     &                                  INFO, LP, LPOK, KEEP10 )
      IMPLICIT NONE
      INTEGER, INTENT(IN)     :: LP
      LOGICAL, INTENT(IN)     :: LPOK
      INTEGER, INTENT(INOUT)  :: INFO(2)
      INTEGER, INTENT(IN)     :: NVTX
      INTEGER, INTENT(OUT)    :: NCMPA
      INTEGER(8), INTENT(IN)  :: NEDGES8
      INTEGER(8)              :: XADJ8(:) 
      INTEGER, INTENT(OUT)    :: NV(NVTX)
      INTEGER                 :: IW(:) 
      INTEGER, INTENT(OUT)    :: PARENT(NVTX)
      INTEGER, INTENT(IN)     :: KEEP10
      INTEGER, DIMENSION(:), ALLOCATABLE :: XADJ
      INTEGER :: I, allocok
      IF (NEDGES8.GT. int(huge(IW),8)) THEN
        INFO(1) = -51
        CALL MUMPS_SET_IERROR(NEDGES8,INFO(2))
        RETURN
      ENDIF
      ALLOCATE(XADJ(NVTX+1), stat=allocok)
      IF (allocok > 0) THEN
        INFO(1)=-7
        INFO(2)=NVTX+1
        IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR memory allocation in MUMPS_PORD_MIXEDto32"
        RETURN
      ENDIF
      CALL MUMPS_ICOPY_64TO32(XADJ8(1), NVTX+1, XADJ)
      CALL MUMPS_PORDF( NVTX, int(NEDGES8), XADJ, IW(1),
     &                                    NV, NCMPA )
      DO I= 1, NVTX
        PARENT(I) = XADJ(I)
      ENDDO
      DEALLOCATE(XADJ)
      RETURN
      END SUBROUTINE MUMPS_PORDF_MIXEDto32
      SUBROUTINE MUMPS_PORDF_MIXEDto64( NVTX, NEDGES8, XADJ8, IW,
     &                                  NV, NCMPA, PARENT,
     &                                  INFO, LP, LPOK, KEEP10,
     &                                  INPLACE64_GRAPH_COPY )
      IMPLICIT NONE
      INTEGER, INTENT(IN)     :: LP
      LOGICAL, INTENT(IN)     :: LPOK
      INTEGER, INTENT(INOUT)  :: INFO(2)
      INTEGER, INTENT(IN)     :: NVTX
      INTEGER, INTENT(OUT)    :: NCMPA
      INTEGER(8), INTENT(IN)  :: NEDGES8
      INTEGER(8)              :: XADJ8(:) 
      INTEGER, INTENT(OUT)    :: NV(NVTX)
      INTEGER, INTENT(IN)     :: IW(:)
      INTEGER, INTENT(OUT)    :: PARENT(NVTX)
      INTEGER, INTENT(IN)     :: KEEP10
      LOGICAL, INTENT(IN)     :: INPLACE64_GRAPH_COPY
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IW8, NV8
      INTEGER :: I, allocok
      IF (KEEP10.EQ.1) THEN
        CALL MUMPS_PORDF( int(NVTX,8), NEDGES8, XADJ8(1), IW(1),
     &                                        NV, NCMPA )
        DO I=1, NVTX
          PARENT(I)=int(XADJ8(I))
        ENDDO
      ELSE
        IF (INPLACE64_GRAPH_COPY)  THEN
          CALL MUMPS_ICOPY_32TO64_64C_IP(IW(1), NEDGES8)
        ELSE
         ALLOCATE(IW8(NEDGES8), stat=allocok)
         IF (allocok > 0) THEN
          INFO(1)=-7
          CALL MUMPS_SET_IERROR(NEDGES8,INFO(2))
          IF (LPOK) WRITE(LP,'(A)') 
     &      "ERROR memory allocation in MUMPS_PORD_MIXEDto64"
          RETURN
         ENDIF
         CALL MUMPS_ICOPY_32TO64_64C(IW(1), NEDGES8, IW8)
        ENDIF
        ALLOCATE(NV8(NVTX), stat=allocok)
        IF (allocok > 0) THEN
          INFO(1)=-7
          CALL MUMPS_SET_IERROR(int(NVTX,8),INFO(2))
          IF (LPOK) WRITE(LP,'(A)') 
     &      "ERROR memory allocation in MUMPS_PORD_MIXEDto64"
          RETURN
        ENDIF
        IF (INPLACE64_GRAPH_COPY)  THEN
         CALL MUMPS_PORDF( int(NVTX,8), NEDGES8, XADJ8(1), IW(1),
     &                                        NV8, NCMPA )
        ELSE
         CALL MUMPS_PORDF( int(NVTX,8), NEDGES8, XADJ8(1), IW8,
     &                                        NV8, NCMPA )
         DEALLOCATE(IW8)
        ENDIF
        CALL MUMPS_ICOPY_64TO32(XADJ8(1), NVTX, PARENT)
        CALL MUMPS_ICOPY_64TO32(NV8, NVTX, NV)
        DEALLOCATE(NV8)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_PORDF_MIXEDto64
      SUBROUTINE MUMPS_PORDF_WND_MIXEDto32( NVTX, NEDGES8,
     &                                      XADJ8, IW,
     &                                      NV, NCMPA, N, PARENT, 
     &                                      INFO, LP, LPOK, KEEP10 )
      IMPLICIT NONE
      INTEGER, INTENT(IN)     :: LP
      LOGICAL, INTENT(IN)     :: LPOK
      INTEGER, INTENT(INOUT)  :: INFO(2)
      INTEGER, INTENT(IN)     :: NVTX, N         
      INTEGER, INTENT(OUT)    :: NCMPA           
      INTEGER, INTENT(INOUT)  :: NV(NVTX)           
      INTEGER(8)              :: XADJ8(:) 
      INTEGER(8), INTENT(IN)  :: NEDGES8         
      INTEGER                 :: IW(:) 
      INTEGER, INTENT(OUT)    :: PARENT(NVTX)
      INTEGER, INTENT(IN)     :: KEEP10
      INTEGER, DIMENSION(:), ALLOCATABLE :: XADJ
      INTEGER :: I, allocok
      IF (NEDGES8.GT. int(huge(IW),8)) THEN
        INFO(1) = -51
        CALL MUMPS_SET_IERROR(NEDGES8,INFO(2))
        RETURN
      ENDIF
      ALLOCATE(XADJ(NVTX+1), stat=allocok)
      IF (allocok > 0) THEN
        INFO(1)=-7
        INFO(2)=NVTX+1
        IF (LPOK) WRITE(LP,'(A)') 
     &    "ERROR memory allocation in MUMPS_PORDF_WND_MIXEDto32"
        RETURN
      ENDIF
      CALL MUMPS_ICOPY_64TO32(XADJ8(1),NVTX+1,XADJ)
      CALL MUMPS_PORDF_WND( NVTX, int(NEDGES8),
     &                             XADJ, IW(1),
     &                             NV, NCMPA, N )
      DO I= 1, NVTX
        PARENT(I) = XADJ(I)
      ENDDO
      DEALLOCATE(XADJ)
      RETURN
      END SUBROUTINE MUMPS_PORDF_WND_MIXEDto32
      SUBROUTINE MUMPS_PORDF_WND_MIXEDto64( NVTX, NEDGES8,
     &                                      XADJ8, IW,
     &                                      NV, NCMPA, N, PARENT, 
     &                                      INFO, LP, LPOK, KEEP10, 
     &                                      INPLACE64_GRAPH_COPY )
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: LP
      LOGICAL, INTENT(IN)      :: LPOK
      INTEGER, INTENT(INOUT)   :: INFO(2)
      INTEGER, INTENT(IN)      :: NVTX, N         
      INTEGER, INTENT(OUT)     :: NCMPA           
      INTEGER, INTENT(INOUT)   :: NV(NVTX)        
      INTEGER(8)               :: XADJ8(:) 
      INTEGER(8), INTENT(IN)   :: NEDGES8         
      INTEGER                  :: IW(:)
      INTEGER, INTENT(OUT)     :: PARENT(NVTX)
      INTEGER, INTENT(IN)      :: KEEP10
      LOGICAL, INTENT(IN)     :: INPLACE64_GRAPH_COPY
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IW8, NV8
      INTEGER :: allocok
      IF (KEEP10.EQ.1) THEN
        CALL MUMPS_PORDF_WND( int(NVTX,8), NEDGES8,
     &                             XADJ8(1), IW(1),
     &                             NV, NCMPA, int(N,8) )
        CALL MUMPS_ICOPY_64TO32(XADJ8(1), NVTX, PARENT)
      ELSE
        IF (INPLACE64_GRAPH_COPY)  THEN
          CALL MUMPS_ICOPY_32TO64_64C_IP(IW(1), NEDGES8)
        ELSE
          ALLOCATE(IW8(NEDGES8), stat=allocok)
          IF (allocok > 0) THEN
            INFO(1)=-7
            CALL MUMPS_SET_IERROR(NEDGES8,INFO(2))
            IF (LPOK) WRITE(LP,'(A)') 
     &      "ERROR memory allocation in MUMPS_PORDF_WND_MIXEDto64"
            RETURN
          ENDIF
          CALL MUMPS_ICOPY_32TO64_64C(IW(1), NEDGES8, IW8)
        ENDIF
        ALLOCATE(NV8(NVTX), stat=allocok)
        IF (allocok > 0) THEN
            INFO(1)=-7
            CALL MUMPS_SET_IERROR(int(NVTX,8),INFO(2))
            IF (LPOK) WRITE(LP,'(A)') 
     &      "ERROR memory allocation in MUMPS_PORDF_WND_MIXEDto64"
            RETURN
        ENDIF
        CALL MUMPS_ICOPY_32TO64(NV, NVTX, NV8)
        IF (INPLACE64_GRAPH_COPY)  THEN
          CALL MUMPS_PORDF_WND( int(NVTX,8), NEDGES8,
     &                             XADJ8(1), IW(1),
     &                             NV8, NCMPA, int(N,8) )
        ELSE
          CALL MUMPS_PORDF_WND( int(NVTX,8), NEDGES8,
     &                             XADJ8(1), IW8,
     &                             NV8, NCMPA, int(N,8) )
          DEALLOCATE(IW8)
        ENDIF
        CALL MUMPS_ICOPY_64TO32(XADJ8(1), NVTX, PARENT)
        CALL MUMPS_ICOPY_64TO32(NV8, NVTX, NV)
        DEALLOCATE(NV8)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_PORDF_WND_MIXEDto64
#endif
      SUBROUTINE MUMPS_ANA_WRAP_RETURN()
      RETURN
      END SUBROUTINE MUMPS_ANA_WRAP_RETURN
      END MODULE MUMPS_ANA_ORD_WRAPPERS
