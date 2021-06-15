C
C  This file is part of MUMPS 5.4.0, released
C  on Tue Apr 13 15:26:30 UTC 2021
C
      PROGRAM MUMPS_TEST_SAVE_RESTORE
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'smumps_struc.h'
      TYPE (SMUMPS_STRUC) mumps_par_save, mumps_par_restore
      INTEGER IERR, I
      CALL MPI_INIT(IERR)
C Define a communicator for the package.
      mumps_par_save%COMM = MPI_COMM_WORLD
C  Initialize an instance of the package
C  for L U factorization (sym = 0, with working host)
      mumps_par_save%JOB = -1
      mumps_par_save%SYM = 0
      mumps_par_save%PAR = 1
      CALL SMUMPS(mumps_par_save)
      IF (mumps_par_save%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &        "  mumps_par_save%INFOG(1)= ", mumps_par_save%INFOG(1), 
     &        "  mumps_par_save%INFOG(2)= ", mumps_par_save%INFOG(2) 
       GOTO 500
      END IF
C  Define problem on the host (processor 0)
      IF ( mumps_par_save%MYID .eq. 0 ) THEN
        READ(5,*) mumps_par_save%N
        READ(5,*) mumps_par_save%NZ
        ALLOCATE( mumps_par_save%IRN ( mumps_par_save%NZ ) )
        ALLOCATE( mumps_par_save%JCN ( mumps_par_save%NZ ) )
        ALLOCATE( mumps_par_save%A( mumps_par_save%NZ ) )
        DO I = 1, mumps_par_save%NZ
          READ(5,*) mumps_par_save%IRN(I),mumps_par_save%JCN(I)
     &          ,mumps_par_save%A(I)
        END DO
      END IF
C  Activate OOC
      mumps_par_save%ICNTL(22)=1
C  Call package for factorization
      mumps_par_save%JOB = 4
      CALL SMUMPS(mumps_par_save)
      IF (mumps_par_save%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &        "  mumps_par_save%INFOG(1)= ", mumps_par_save%INFOG(1), 
     &        "  mumps_par_save%INFOG(2)= ", mumps_par_save%INFOG(2) 
       GOTO 500
      END IF
C  Call package for save
      mumps_par_save%JOB = 7
      mumps_par_save%SAVE_DIR="/tmp"
      mumps_par_save%SAVE_PREFIX="mumps_simpletest_save"
      CALL SMUMPS(mumps_par_save)
      IF (mumps_par_save%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &        "  mumps_par_save%INFOG(1)= ", mumps_par_save%INFOG(1), 
     &        "  mumps_par_save%INFOG(2)= ", mumps_par_save%INFOG(2) 
       GOTO 500
      END IF
C  Deallocate user data
      IF ( mumps_par_save%MYID .eq. 0 )THEN
        DEALLOCATE( mumps_par_save%IRN )
        DEALLOCATE( mumps_par_save%JCN )
        DEALLOCATE( mumps_par_save%A   )
      END IF
C  Destroy the instance (deallocate internal data structures)
      mumps_par_save%JOB = -2
      CALL SMUMPS(mumps_par_save)
C  Now mumps_par_save has be destroyed
C  We use a new instance mumps_par_restore to finish the computation
C  Define a communicator for the package on the new instace.
      mumps_par_restore%COMM = MPI_COMM_WORLD
C  Initialize a new instance of the package
C  for L U factorization (sym = 0, with working host)
      mumps_par_restore%JOB = -1
      mumps_par_restore%SYM = 0
      mumps_par_restore%PAR = 1
      CALL SMUMPS(mumps_par_restore)
      IF (mumps_par_restore%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &        "  mumps_par_restore%INFOG(1)= ", 
     &        mumps_par_restore%INFOG(1), 
     &        "  mumps_par_restore%INFOG(2)= ", 
     &        mumps_par_restore%INFOG(2) 
       GOTO 500
      END IF
C  Call package for restore with OOC feature
      mumps_par_restore%JOB = 8
      mumps_par_restore%SAVE_DIR="/tmp"
      mumps_par_restore%SAVE_PREFIX="mumps_simpletest_save"
      CALL SMUMPS(mumps_par_restore)
      IF (mumps_par_restore%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &        "  mumps_par_restore%INFOG(1)= ",
     &        mumps_par_restore%INFOG(1), 
     &        "  mumps_par_restore%INFOG(2)= ", 
     &        mumps_par_restore%INFOG(2) 
       GOTO 500
      END IF
C  Define rhs on the host (processor 0)
      IF ( mumps_par_restore%MYID .eq. 0 ) THEN
        ALLOCATE( mumps_par_restore%RHS ( mumps_par_restore%N  ) )
        DO I = 1, mumps_par_restore%N
           READ(5,*) mumps_par_restore%RHS(I)
        END DO
      END IF
C  Call package for solution
      mumps_par_restore%JOB = 3
      CALL SMUMPS(mumps_par_restore)
      IF (mumps_par_restore%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  mumps_par_restore%INFOG(1)= ", 
     &        mumps_par_restore%INFOG(1), 
     &            "  mumps_par_restore%INFOG(2)= ", 
     &        mumps_par_restore%INFOG(2) 
       GOTO 500
      END IF
C  Solution has been assembled on the host
      IF ( mumps_par_restore%MYID .eq. 0 ) THEN
        WRITE( 6, * ) ' Solution is ',
     &        (mumps_par_restore%RHS(I),I=1,mumps_par_restore%N)
      END IF
C  Deallocate user data
      IF ( mumps_par_restore%MYID .eq. 0 )THEN
        DEALLOCATE( mumps_par_restore%RHS )
      END IF
C  Delete the saved files
C  Note mumps_par_restore%ICNTL(34) is kept to default (0) to suppress
C  also the OOC files.
      mumps_par_restore%JOB = -3
      CALL SMUMPS(mumps_par_restore)
C  Destroy the instance (deallocate internal data structures)
      mumps_par_restore%JOB = -2
      CALL SMUMPS(mumps_par_restore)
 500  CALL MPI_FINALIZE(IERR)
      STOP
      END
