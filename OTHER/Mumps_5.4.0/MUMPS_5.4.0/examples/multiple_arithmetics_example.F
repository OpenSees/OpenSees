C
C  This file is part of MUMPS 5.4.0, released
C  on Tue Apr 13 15:26:30 UTC 2021
C
      PROGRAM MUMPS_MULTIPLE_ARITHMETICS_TEST
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'smumps_struc.h'
      INCLUDE 'dmumps_struc.h'
      INCLUDE 'cmumps_struc.h'
      INCLUDE 'zmumps_struc.h'
      TYPE (SMUMPS_STRUC) smumps_par
      TYPE (DMUMPS_STRUC) dmumps_par
      TYPE (CMUMPS_STRUC) cmumps_par
      TYPE (ZMUMPS_STRUC) zmumps_par
      INTEGER IERR
      CALL MPI_INIT(IERR)
C Define a communicator for the packages.
      smumps_par%COMM = MPI_COMM_WORLD
      dmumps_par%COMM = smumps_par%COMM
      cmumps_par%COMM = smumps_par%COMM
      zmumps_par%COMM = smumps_par%COMM
C  Initialize all instances of the package
C  for L U factorization (sym = 0, with working host)
      smumps_par%JOB = -1
      smumps_par%SYM = 0
      smumps_par%PAR = 1
      CALL SMUMPS(smumps_par)
      IF (smumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  smumps_par%INFOG(1)= ", smumps_par%INFOG(1), 
     &            "  smumps_par%INFOG(2)= ", smumps_par%INFOG(2) 
       GOTO 500
      END IF

      dmumps_par%JOB = smumps_par%JOB
      dmumps_par%SYM = smumps_par%SYM
      dmumps_par%PAR = smumps_par%PAR
      cmumps_par%JOB = smumps_par%JOB
      cmumps_par%SYM = smumps_par%SYM
      cmumps_par%PAR = smumps_par%PAR
      zmumps_par%JOB = smumps_par%JOB
      zmumps_par%SYM = smumps_par%SYM
      zmumps_par%PAR = smumps_par%PAR
      
      CALL DMUMPS(dmumps_par)
      IF (dmumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  dmumps_par%INFOG(1)= ", dmumps_par%INFOG(1), 
     &            "  dmumps_par%INFOG(2)= ", dmumps_par%INFOG(2) 
       GOTO 500
      END IF

      CALL CMUMPS(cmumps_par)
      IF (cmumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  cmumps_par%INFOG(1)= ", cmumps_par%INFOG(1), 
     &            "  cmumps_par%INFOG(2)= ", cmumps_par%INFOG(2) 
       GOTO 500
      END IF

      CALL ZMUMPS(zmumps_par)
      IF (zmumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  zmumps_par%INFOG(1)= ", zmumps_par%INFOG(1), 
     &            "  zmumps_par%INFOG(2)= ", zmumps_par%INFOG(2) 
       GOTO 500
      END IF

      IF ( smumps_par%MYID .eq. 0 )THEN
         write(6,'(A)') "Creation of all instaces went well"
      ENDIF
           
C     Destroy the instances (deallocate internal data structures)
      smumps_par%JOB = -2
      CALL SMUMPS(smumps_par)
      IF (smumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  smumps_par%INFOG(1)= ", smumps_par%INFOG(1), 
     &            "  smumps_par%INFOG(2)= ", smumps_par%INFOG(2) 
       GOTO 500
      END IF

      dmumps_par%JOB = smumps_par%JOB
      cmumps_par%JOB = smumps_par%JOB
      zmumps_par%JOB = smumps_par%JOB
      
      CALL DMUMPS(dmumps_par)
      IF (dmumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  dmumps_par%INFOG(1)= ", dmumps_par%INFOG(1), 
     &            "  dmumps_par%INFOG(2)= ", dmumps_par%INFOG(2) 
       GOTO 500
      END IF
      CALL CMUMPS(cmumps_par)
      IF (cmumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  cmumps_par%INFOG(1)= ", cmumps_par%INFOG(1), 
     &            "  cmumps_par%INFOG(2)= ", cmumps_par%INFOG(2) 
       GOTO 500
      END IF
      CALL ZMUMPS(zmumps_par)
      IF (zmumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  zmumps_par%INFOG(1)= ", zmumps_par%INFOG(1), 
     &            "  zmumps_par%INFOG(2)= ", zmumps_par%INFOG(2) 
       GOTO 500
      END IF
      
 500  CALL MPI_FINALIZE(IERR)
      STOP
      END PROGRAM MUMPS_MULTIPLE_ARITHMETICS_TEST

