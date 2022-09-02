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
      SUBROUTINE DMUMPS_GET_NS_OPTIONS_FACTO(N,KEEP,ICNTL,MPG)
      IMPLICIT NONE
          INTEGER N, KEEP(500), ICNTL(60), MPG
          KEEP(19)=0
          RETURN
      END SUBROUTINE DMUMPS_GET_NS_OPTIONS_FACTO
      SUBROUTINE DMUMPS_GET_NS_OPTIONS_SOLVE(ICNTL, KEEP, 
     &                                       NRHS, MPG, INFO)
      IMPLICIT NONE
      INTEGER, intent(in)   :: KEEP(500), NRHS, MPG, ICNTL(60)
      INTEGER, intent(inout):: INFO(80)
      IF (KEEP(19).EQ.0.AND.KEEP(110).EQ.0) THEN
        IF (KEEP(111).NE.0) THEN
         INFO(1) = -37
         INFO(2) = 56  
         IF (KEEP(110).EQ.0) INFO(2) = 24 
          IF(MPG.GT.0) THEN
           WRITE( MPG,'(A)')
     &'** ERROR  : Null space computation requirement'
          WRITE( MPG,'(A)')
     &'** not consistent with factorization options'
         ENDIF
         GOTO 333
        ENDIF
      ENDIF
       IF (ICNTL(9).NE.1) THEN
         IF (KEEP(111).NE.0) THEN
          INFO(1) = -37
          INFO(2) = 9  
          IF (MPG.GT.0) THEN
           WRITE(MPG,'(A)')
     &'** ERROR  ICNTL(25) incompatible with '
           WRITE( MPG,'(A)')
     &'** option transposed system (ICNTL(9).ne.1) '
          ENDIF
         ENDIF
         GOTO 333
       ENDIF
      IF (KEEP(19).EQ.2) THEN
       IF ((KEEP(111).NE.0).AND.(KEEP(50).EQ.0)) THEN
         INFO(1) = -37
         INFO(2) = 0  
         IF (MPG.GT.0) THEN
          WRITE(MPG,'(A)')
     &'** ERROR  ICNTL(25) incompatible with '
          WRITE( MPG,'(A)')
     &'** option RRQR (ICNLT(56)=2) and unsym. matrices '
         ENDIF
        ENDIF
        GOTO 333
      ENDIF
      IF (KEEP(111).eq.-1.AND.NRHS.NE.KEEP(112)+KEEP(17))THEN
            INFO(1)=-32
            INFO(2)=NRHS
            GOTO 333
      ENDIF
      IF (KEEP(111).gt.0 .AND. NRHS .NE. 1) THEN
            INFO(1)=-32
            INFO(2)=NRHS
            GOTO 333
      ENDIF
      IF (KEEP(248) .NE.0.AND.KEEP(111).NE.0) THEN
         IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     &           ' ERROR: ICNTL(20) and ICNTL(30) functionalities ',
     &           ' incompatible with null space'
         ENDIF
         INFO(1) = -37
         IF (KEEP(237).NE.0) THEN
            INFO(2) = 30    
            IF (MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     &           ' ERROR: ICNTL(30) functionality ',
     &              ' incompatible with null space'
            ENDIF
         ELSE
            IF (MPG.GT.0) THEN
               WRITE(MPG,'(A)')
     &              ' ERROR: ICNTL(20) functionality ',
     &              ' incompatible with null space'
            ENDIF
            INFO(2) = 20    
         ENDIF
         GOTO 333
      ENDIF
      IF (( KEEP(111) .LT. -1 ) .OR.
     &     (KEEP(111).GT.KEEP(112)+KEEP(17)) .OR.
     &     (KEEP(111) .EQ.-1 .AND. KEEP(112)+KEEP(17).EQ.0))
     &     THEN
         INFO(1)=-36
         INFO(2)=KEEP(111)
         GOTO 333
      ENDIF
      IF (KEEP(221).NE.0.AND.KEEP(111).NE.0) THEN
         INFO(1)=-37
         INFO(2)=26
         GOTO 333
      ENDIF
 333  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_GET_NS_OPTIONS_SOLVE
      SUBROUTINE DMUMPS_RR_INIT_POINTERS(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (DMUMPS_STRUC) id
      NULLIFY(id%root%QR_TAU)
      NULLIFY(id%root%SVD_U)
      NULLIFY(id%root%SVD_VT)
      NULLIFY(id%root%SINGULAR_VALUES)
      RETURN
      END SUBROUTINE DMUMPS_RR_INIT_POINTERS
      SUBROUTINE DMUMPS_RR_FREE_POINTERS(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (DMUMPS_STRUC) id
      IF (associated(id%root%QR_TAU))  THEN
        DEALLOCATE(id%root%QR_TAU)
        NULLIFY(id%root%QR_TAU)
      ENDIF
      IF (associated(id%root%SVD_U))  THEN
        DEALLOCATE(id%root%SVD_U)
        NULLIFY(id%root%SVD_U)
      ENDIF
      IF (associated(id%root%SVD_VT))  THEN
        DEALLOCATE(id%root%SVD_VT)
        NULLIFY(id%root%SVD_VT)
      ENDIF
      IF (associated(id%root%SINGULAR_VALUES))  THEN
        DEALLOCATE(id%root%SINGULAR_VALUES)
        NULLIFY(id%root%SINGULAR_VALUES)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_RR_FREE_POINTERS
