C
C  This file is part of MUMPS 5.4.0, released
C  on Tue Apr 13 15:26:30 UTC 2021
C
      PROGRAM MUMPS_TEST
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_struc.h'
      TYPE (CMUMPS_STRUC) mumps_par
      INTEGER IERR, I
      INTEGER(8) I8
      CALL MPI_INIT(IERR)
C Define a communicator for the package.
      mumps_par%COMM = MPI_COMM_WORLD
C  Initialize an instance of the package
C  for L U factorization (sym = 0, with working host)
      mumps_par%JOB = -1
      mumps_par%SYM = 0
      mumps_par%PAR = 1
      CALL CMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), 
     &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
C  Define problem on the host (processor 0)
      IF ( mumps_par%MYID .eq. 0 ) THEN
        READ(5,*) mumps_par%N
        READ(5,*) mumps_par%NNZ
        ALLOCATE( mumps_par%IRN ( mumps_par%NNZ ) )
        ALLOCATE( mumps_par%JCN ( mumps_par%NNZ ) )
        ALLOCATE( mumps_par%A( mumps_par%NNZ ) )
        ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
        DO I8 = 1, mumps_par%NNZ
          READ(5,*) mumps_par%IRN(I8),mumps_par%JCN(I8),mumps_par%A(I8)
        END DO
        DO I = 1, mumps_par%N
          READ(5,*) mumps_par%RHS(I)
        END DO
      END IF
C  Call package for solution
      mumps_par%JOB = 6
      CALL CMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), 
     &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
C  Solution has been assembled on the host
      IF ( mumps_par%MYID .eq. 0 ) THEN
        WRITE( 6, * ) ' Solution is ',(mumps_par%RHS(I),I=1,mumps_par%N)
      END IF
C  Deallocate user data
      IF ( mumps_par%MYID .eq. 0 )THEN
        DEALLOCATE( mumps_par%IRN )
        DEALLOCATE( mumps_par%JCN )
        DEALLOCATE( mumps_par%A   )
        DEALLOCATE( mumps_par%RHS )
      END IF
C  Destroy the instance (deallocate internal data structures)
      mumps_par%JOB = -2
      CALL CMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), 
     &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
 500  CALL MPI_FINALIZE(IERR)
      STOP
      END
