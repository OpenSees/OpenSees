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
      MODULE MUMPS_LR_COMMON
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE COMPUTE_BLR_VCS(K472, IBCKSZ, MAXSIZE, NASS)
        INTEGER, INTENT(IN) :: MAXSIZE, NASS, K472
        INTEGER, INTENT(OUT) :: IBCKSZ
        IF (K472.EQ.1) THEN
          IF (NASS.LE.1000) THEN
            IBCKSZ = 128
          ELSEIF (NASS.GT.1000.AND.NASS.LE.5000) THEN
            IBCKSZ = 256
          ELSEIF (NASS.GT.5000.AND.NASS.LE.10000) THEN
            IBCKSZ = 384
          ELSE
            IBCKSZ = 512
          ENDIF
          IBCKSZ = min(IBCKSZ,MAXSIZE)
        ELSE
          IBCKSZ = MAXSIZE
        ENDIF
      END SUBROUTINE COMPUTE_BLR_VCS
      SUBROUTINE MUMPS_UPD_TREE(NV, NSTEPS, N, FIRST, LPTR, RPTR, F,
     &     VLIST, FILS, FRERE_STEPS, STEP, DAD_STEPS, NE_STEPS, NA, LNA,
     &     PVS, K38, STEP_SCALAPACK_ROOT)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, NV, NSTEPS, LNA, F, VLIST(NV)
      INTEGER                :: FILS(:), FRERE_STEPS(:),
     &     DAD_STEPS(:), STEP(:), NE_STEPS(:), NA(:)
      INTEGER, INTENT(INOUT) :: PVS(NSTEPS), LPTR, RPTR
      INTEGER, INTENT(INOUT) :: K38
      INTEGER, INTENT(IN)    :: STEP_SCALAPACK_ROOT
      LOGICAL :: FIRST
      INTEGER :: PV, NODE, I
      PV        = VLIST(1)
      NODE      = ABS(STEP(PV))
      PVS(NODE) = PV
      IF(FIRST) THEN
         I = DAD_STEPS(NODE)
         DO WHILE(FILS(I).GT.0)
            I = FILS(I)
         END DO
         FILS(I) = -PV
      END IF
      IF(FRERE_STEPS(NODE) .GT. 0) THEN
         FRERE_STEPS(NODE) = PVS(ABS(STEP(FRERE_STEPS(NODE))))
      ELSE IF(FRERE_STEPS(NODE) .LT. 0) THEN
         FRERE_STEPS(NODE) = -PVS(ABS(STEP(DAD_STEPS(NODE))))
      END IF
      IF(DAD_STEPS(NODE) .EQ. 0) THEN
         NA(RPTR) = PV
         RPTR     = RPTR -1
      ELSE
         DAD_STEPS(NODE) = PVS(ABS(STEP(DAD_STEPS(NODE))))
      END IF
      IF(NE_STEPS(NODE) .EQ. 0) THEN
         NA(LPTR) = PV
         LPTR     = LPTR -1
      END IF
      STEP(VLIST(1)) = ABS(STEP(VLIST(1)))
      IF (STEP(VLIST(1)).EQ.STEP_SCALAPACK_ROOT) THEN
       K38 = VLIST(1)
      ENDIF
      DO I=1, NV-1
         IF(STEP(VLIST(I+1)).GT.0) STEP(VLIST(I+1)) = -STEP(VLIST(I+1))
         FILS(VLIST(I)) = VLIST(I+1)
      END DO
      FILS(VLIST(NV)) = F
      RETURN
      END SUBROUTINE MUMPS_UPD_TREE
      END MODULE MUMPS_LR_COMMON
