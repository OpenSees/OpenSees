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
      SUBROUTINE DMUMPS_MCAST2(DATA, LDATA, MPITYPE, ROOT, COMMW, TAG,
     &SLAVEF, KEEP)
      USE DMUMPS_BUF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER LDATA, ROOT, COMMW, TAG, MPITYPE, SLAVEF
      INTEGER DEST
      INTEGER, INTENT(INOUT) :: KEEP(500)
      INTEGER DATA(LDATA)
      DO 10 DEST = 0, SLAVEF - 1
        IF (DEST .NE. ROOT) THEN
          IF ( LDATA .EQ. 1 .and. MPITYPE .EQ. MPI_INTEGER ) THEN
            CALL DMUMPS_BUF_SEND_1INT( DATA(1), DEST, TAG, 
     &                                COMMW, KEEP, IERR )
          ELSE
            WRITE(*,*) 'Error : bad argument to DMUMPS_MCAST2'
            CALL MUMPS_ABORT()
          END IF
        ENDIF
   10 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_MCAST2
      SUBROUTINE DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      INTEGER MYID, SLAVEF, COMM
      INTEGER, INTENT(INOUT) :: KEEP(500)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY (1)
      DUMMY(1) = -98765
      CALL DMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID,
     &                 COMM, TERREUR, SLAVEF, KEEP )
      RETURN
      END SUBROUTINE DMUMPS_BDC_ERROR
