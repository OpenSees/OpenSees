C
C  This file is part of MUMPS 5.4.0, released
C  on Tue Apr 13 15:26:30 UTC 2021
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
      INTEGER FUNCTION CMUMPS_IXAMAX(N,X,INCX,GRAIN)
!$    USE OMP_LIB
      IMPLICIT NONE
      COMPLEX, intent(in) :: X(*)
      INTEGER, intent(in) :: INCX,N
      INTEGER, intent(in) :: GRAIN
      REAL ABSMAX
      INTEGER :: I
      INTEGER(8) :: IX
!$    INTEGER :: NOMP, CHUNK
!$    INTEGER :: IMAX
!$    REAL    :: XMAX, VALABS
!$    REAL, PARAMETER :: RZERO = 0.0E0
!$    NOMP = OMP_GET_MAX_THREADS()
      CMUMPS_IXAMAX = 0
      IF ( N.LT.1 ) RETURN
      CMUMPS_IXAMAX = 1
      IF ( N.EQ.1 .OR. INCX.LE.0 ) RETURN
!$    IF (NOMP.GT.1 .AND. N.GE.GRAIN*2) THEN
!$    IF ( INCX.EQ.1 ) THEN
!$      CHUNK = max(GRAIN,(N+NOMP-1)/NOMP)
!$      ABSMAX = RZERO
!$OMP PARALLEL PRIVATE(I, VALABS, XMAX, IMAX)
!$OMP& FIRSTPRIVATE(N, CHUNK)
!$      XMAX = RZERO
!$OMP DO SCHEDULE(static, CHUNK)
!$      DO I = 1, N
!$        VALABS = abs(X(I))
!$        IF ( VALABS .GT. XMAX ) THEN
!$          XMAX = VALABS
!$          IMAX = I
!$        ENDIF
!$      ENDDO
!$OMP END DO
!$      IF (XMAX .GT. RZERO) THEN
!$OMP CRITICAL
!$        IF (XMAX .GT. ABSMAX) THEN
!$          CMUMPS_IXAMAX = IMAX
!$          ABSMAX = XMAX
!$        ENDIF
!$OMP END CRITICAL
!$      ENDIF
!$OMP END PARALLEL
!$    ELSE
!$      CHUNK = max(GRAIN,(N+NOMP-1)/NOMP)
!$      ABSMAX = RZERO
!$OMP PARALLEL PRIVATE(I, VALABS, XMAX, IMAX, IX)
!$OMP& FIRSTPRIVATE(N, CHUNK, INCX)
!$      XMAX = RZERO
!$OMP DO SCHEDULE(static, CHUNK)
!$      DO I = 1, N
!$        IX = 1 + int((I-1),8)*int(INCX,8)
!$        VALABS = abs(X(IX))
!$        IF ( VALABS .GT. XMAX ) THEN
!$          XMAX = VALABS
!$          IMAX = I
!$        ENDIF
!$      ENDDO
!$OMP END DO
!$      IF (XMAX .GT. RZERO) THEN
!$OMP CRITICAL
!$        IF (XMAX .GT. ABSMAX) THEN
!$          CMUMPS_IXAMAX = IMAX
!$          ABSMAX = XMAX
!$        ENDIF
!$OMP END CRITICAL
!$      ENDIF
!$OMP END PARALLEL
!$    ENDIF
!$    ELSE
      IF ( INCX.EQ.1 ) THEN
        ABSMAX = abs(X(1))
        DO I = 2, N
          IF ( abs(X(I)) .LE. ABSMAX ) CYCLE
          CMUMPS_IXAMAX = I
          ABSMAX = abs(X(I))
        ENDDO
      ELSE
        IX = 1
        ABSMAX = abs(X(1))
        IX = IX + INCX
        DO I = 2, N
           IF ( abs(X(IX)).LE.ABSMAX ) GOTO 5
           CMUMPS_IXAMAX = I
           ABSMAX = abs(X(IX))
   5       IX = IX + INCX
        ENDDO
      ENDIF
!$    ENDIF
      RETURN
      END FUNCTION CMUMPS_IXAMAX
