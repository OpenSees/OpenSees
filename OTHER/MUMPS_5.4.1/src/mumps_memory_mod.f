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
      MODULE MUMPS_MEMORY_MOD
      INTERFACE MUMPS_DEALLOC
      MODULE PROCEDURE MUMPS_IDEALLOC
      END INTERFACE
      INTERFACE MUMPS_REALLOC
      MODULE PROCEDURE MUMPS_IREALLOC
      MODULE PROCEDURE MUMPS_DREALLOC, MUMPS_SREALLOC, MUMPS_ZREALLOC
      MODULE PROCEDURE MUMPS_CREALLOC
      END INTERFACE
      INTEGER(8), PRIVATE :: ISIZE, I8SIZE, SSIZE, DSIZE, CSIZE, ZSIZE
      CONTAINS
      SUBROUTINE MUMPS_MEMORY_SET_DATA_SIZES()
      INTEGER             :: I(2)
      INTEGER(8)          :: I8(2)
      REAL(kind(1.e0))    :: S(2)
      REAL(kind(1.d0))    :: D(2)
      COMPLEX(kind(1.e0)) :: C(2)
      COMPLEX(kind(1.d0)) :: Z(2)
      CALL MUMPS_SIZE_C(I (1), I (2), ISIZE) 
      CALL MUMPS_SIZE_C(S (1), S (2), SSIZE) 
      CALL MUMPS_SIZE_C(D (1), D (2), DSIZE) 
      CALL MUMPS_SIZE_C(C (1), C (2), CSIZE) 
      CALL MUMPS_SIZE_C(Z (1), Z (2), ZSIZE) 
      CALL MUMPS_SIZE_C(I8(1), I8(2), I8SIZE) 
      RETURN
      END SUBROUTINE MUMPS_MEMORY_SET_DATA_SIZES
      SUBROUTINE MUMPS_IREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      INTEGER, POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      INTEGER, POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*ISIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ISIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ISIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+
     &           MINSIZE*ISIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_IREALLOC
      SUBROUTINE MUMPS_I8REALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      INTEGER(8), POINTER          :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      INTEGER(8), POINTER          :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*I8SIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*I8SIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*I8SIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+
     &           int(MINSIZE,8)*I8SIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_I8REALLOC
      SUBROUTINE MUMPS_IREALLOC8(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      INTEGER, POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: LP
      INTEGER(8)                   :: MINSIZE
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      INTEGER, POINTER             :: TEMP(:)
      INTEGER(8)                   :: I
      INTEGER                      :: IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = int(min(MINSIZE,huge(I)))
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = int(min(MINSIZE,huge(I)))
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((int(size(ARRAY),8) .LT. MINSIZE) .OR.
     &           ((int(size(ARRAY),8).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE*ISIZE
               END IF
               DO I=1, min(int(size(ARRAY),8), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ISIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((int(size(ARRAY),8) .LT. MINSIZE) .OR.
     &           ((int(size(ARRAY),8).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ISIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE*ISIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_IREALLOC8
      SUBROUTINE MUMPS_I8REALLOC8(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      INTEGER(8), POINTER          :: ARRAY(:)
      INTEGER                      :: INFO(:), LP
      INTEGER(8)                   :: MINSIZE
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      INTEGER(8), POINTER          :: TEMP(:)
      INTEGER                      :: IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      INTEGER(8)                   :: ASIZE, I
      ASIZE = int(size(ARRAY),8)
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = int(MINSIZE)
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = int(MINSIZE)
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((ASIZE .LT. MINSIZE) .OR.
     &           ((ASIZE.NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*I8SIZE
               END IF
               DO I=1, min(ASIZE, MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              ASIZE*I8SIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((ASIZE .LT. MINSIZE) .OR.
     &           ((ASIZE.NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              ASIZE*I8SIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+
     &           int(MINSIZE,8)*I8SIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_I8REALLOC8
      SUBROUTINE MUMPS_SREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      REAL(kind(1.E0)), POINTER    :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      REAL(kind(1.E0)), POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*SSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*SSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*SSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE*SSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_SREALLOC
      SUBROUTINE MUMPS_DREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      REAL(kind(1.D0)), POINTER    :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      REAL(kind(1.D0)), POINTER    :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*DSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*DSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*DSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+
     &           int(MINSIZE,8)*DSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_DREALLOC
      SUBROUTINE MUMPS_CREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      COMPLEX(kind((1.E0,1.E0))), POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      COMPLEX(kind((1.E0,1.E0))), POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*CSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*CSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*CSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+
     &           int(MINSIZE,8)*CSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_CREALLOC
      SUBROUTINE MUMPS_ZREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      COMPLEX(kind((1.D0,1.D0))), POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      COMPLEX(kind((1.D0,1.D0))), POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+int(MINSIZE,8)*16_8
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT =MEMCNT-
     &              int(size(ARRAY),8)*ZSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT =MEMCNT-
     &              int(size(ARRAY),8)*ZSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+
     &           int(MINSIZE,8)*ZSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_ZREALLOC
      SUBROUTINE MUMPS_IDEALLOC(A1, A2, A3, A4, A5, A6, A7, MEMCNT)
      INTEGER, POINTER :: A1(:)
      INTEGER, POINTER, OPTIONAL :: A2(:), A3(:), A4(:), A5(:),
     &     A6(:), A7(:)
      INTEGER(8), OPTIONAL :: MEMCNT
      INTEGER(8) :: IMEMCNT
      IMEMCNT = 0
      IF(associated(A1)) THEN
         IMEMCNT = IMEMCNT+int(size(A1),8)*ISIZE
         DEALLOCATE(A1)
         NULLIFY(A1)
      END IF
      IF(present(A2)) THEN
         IF(associated(A2)) THEN
            IMEMCNT = IMEMCNT+int(size(A2),8)*ISIZE
            DEALLOCATE(A2)
            NULLIFY(A2)
         END IF
      END IF
      IF(present(A3)) THEN
         IF(associated(A3)) THEN
            IMEMCNT = IMEMCNT+int(size(A3),8)*ISIZE
            DEALLOCATE(A3)
            NULLIFY(A3)
         END IF
      END IF
      IF(present(A4)) THEN
         IF(associated(A4)) THEN
            IMEMCNT = IMEMCNT+int(size(A4),8)*ISIZE
            DEALLOCATE(A4)
            NULLIFY(A4)
         END IF
      END IF
      IF(present(A5)) THEN
         IF(associated(A5)) THEN
            IMEMCNT = IMEMCNT+int(size(A5),8)*ISIZE
            DEALLOCATE(A5)
            NULLIFY(A5)
         END IF
      END IF
      IF(present(A6)) THEN
         IF(associated(A6)) THEN
            IMEMCNT = IMEMCNT+int(size(A6),8)*ISIZE
            DEALLOCATE(A6)
            NULLIFY(A6)
         END IF
      END IF
      IF(present(A7)) THEN
         IF(associated(A7)) THEN
            IMEMCNT = IMEMCNT+int(size(A7),8)*ISIZE
            DEALLOCATE(A7)
            NULLIFY(A7)
         END IF
      END IF
      IF(present(MEMCNT)) MEMCNT = MEMCNT-IMEMCNT
      RETURN
      END SUBROUTINE MUMPS_IDEALLOC
      SUBROUTINE MUMPS_I8DEALLOC(A1, A2, A3, A4, A5, A6, A7, MEMCNT)
      INTEGER(8), POINTER :: A1(:)
      INTEGER(8), POINTER, OPTIONAL :: A2(:), A3(:), A4(:), A5(:),
     &     A6(:), A7(:)
      INTEGER(8), OPTIONAL :: MEMCNT
      INTEGER(8) :: IMEMCNT
      IMEMCNT = 0
      IF(associated(A1)) THEN
         IMEMCNT = IMEMCNT+int(size(A1),8)*I8SIZE
         DEALLOCATE(A1)
         NULLIFY(A1)
      END IF
      IF(present(A2)) THEN
         IF(associated(A2)) THEN
            IMEMCNT = IMEMCNT+int(size(A2),8)*I8SIZE
            DEALLOCATE(A2)
            NULLIFY(A2)
         END IF
      END IF
      IF(present(A3)) THEN
         IF(associated(A3)) THEN
            IMEMCNT = IMEMCNT+int(size(A3),8)*I8SIZE
            DEALLOCATE(A3)
            NULLIFY(A3)
         END IF
      END IF
      IF(present(A4)) THEN
         IF(associated(A4)) THEN
            IMEMCNT = IMEMCNT+int(size(A4),8)*I8SIZE
            DEALLOCATE(A4)
            NULLIFY(A4)
         END IF
      END IF
      IF(present(A5)) THEN
         IF(associated(A5)) THEN
            IMEMCNT = IMEMCNT+int(size(A5),8)*I8SIZE
            DEALLOCATE(A5)
            NULLIFY(A5)
         END IF
      END IF
      IF(present(A6)) THEN
         IF(associated(A6)) THEN
            IMEMCNT = IMEMCNT+int(size(A6),8)*I8SIZE
            DEALLOCATE(A6)
            NULLIFY(A6)
         END IF
      END IF
      IF(present(A7)) THEN
         IF(associated(A7)) THEN
            IMEMCNT = IMEMCNT+int(size(A7),8)*I8SIZE
            DEALLOCATE(A7)
            NULLIFY(A7)
         END IF
      END IF
      IF(present(MEMCNT)) MEMCNT = MEMCNT-IMEMCNT
      RETURN
      END SUBROUTINE MUMPS_I8DEALLOC
      END MODULE
