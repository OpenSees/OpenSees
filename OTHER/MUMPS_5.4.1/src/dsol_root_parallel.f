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
      SUBROUTINE DMUMPS_ROOT_SOLVE( NRHS, DESCA_PAR, 
     &  CNTXT_PAR,LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     &  IPIV,LPIV,MASTER_ROOT,MYID,COMM,
     &  RHS_SEQ,SIZE_ROOT,A,INFO,MTYPE,LDLT )
      IMPLICIT NONE
      INTEGER NRHS, MTYPE
      INTEGER DESCA_PAR( 9 )
      INTEGER LOCAL_M, LOCAL_N, MBLOCK, NBLOCK
      INTEGER CNTXT_PAR, MASTER_ROOT, SIZE_ROOT
      INTEGER MYID, COMM
      INTEGER LPIV, IPIV( LPIV )
      INTEGER INFO(80), LDLT
      DOUBLE PRECISION RHS_SEQ( SIZE_ROOT *NRHS)
      DOUBLE PRECISION A( LOCAL_M, LOCAL_N )
      INTEGER IERR, NPROW, NPCOL, MYROW, MYCOL
      INTEGER LOCAL_N_RHS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,: ) ::RHS_PAR
      EXTERNAL numroc
      INTEGER  numroc
      INTEGER allocok
      CALL blacs_gridinfo( CNTXT_PAR, NPROW, NPCOL, MYROW, MYCOL )
      LOCAL_N_RHS = numroc(NRHS, NBLOCK, MYCOL, 0, NPCOL)
      LOCAL_N_RHS = max(1,LOCAL_N_RHS)
      ALLOCATE(RHS_PAR(LOCAL_M, LOCAL_N_RHS),stat=allocok)
      IF (allocok > 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root.'
        WRITE(*,*) ' Reduce number of right hand sides.'
        CALL MUMPS_ABORT()
      ENDIF
      CALL DMUMPS_SCATTER_ROOT( MYID, SIZE_ROOT, NRHS, RHS_SEQ,
     &      LOCAL_M, LOCAL_N_RHS,
     &      MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     &      NPROW, NPCOL, COMM )
      CALL DMUMPS_SOLVE_2D_BCYCLIC (SIZE_ROOT, NRHS, MTYPE,
     &     A, DESCA_PAR, LOCAL_M, LOCAL_N, LOCAL_N_RHS,
     &     IPIV, LPIV, RHS_PAR, LDLT, 
     &     MBLOCK, NBLOCK, CNTXT_PAR,
     &     IERR)
      CALL DMUMPS_GATHER_ROOT( MYID, SIZE_ROOT, NRHS,
     &    RHS_SEQ, LOCAL_M, LOCAL_N_RHS,
     &    MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     &    NPROW, NPCOL, COMM )
      DEALLOCATE(RHS_PAR)
      RETURN
      END SUBROUTINE DMUMPS_ROOT_SOLVE
      SUBROUTINE DMUMPS_SOLVE_2D_BCYCLIC (SIZE_ROOT, NRHS, MTYPE,
     &     A, DESCA_PAR, LOCAL_M, LOCAL_N, LOCAL_N_RHS,
     &     IPIV, LPIV, RHS_PAR, LDLT, 
     &     MBLOCK, NBLOCK, CNTXT_PAR,
     &     IERR)
      IMPLICIT NONE
      INTEGER, intent (in) :: SIZE_ROOT, NRHS, LDLT, LOCAL_M, 
     &                        LOCAL_N, LOCAL_N_RHS, 
     &                        MBLOCK, NBLOCK, CNTXT_PAR, MTYPE
      INTEGER, intent (in) :: DESCA_PAR( 9 ) 
      INTEGER, intent (in) :: LPIV, IPIV( LPIV )
      DOUBLE PRECISION, intent (in) :: A( LOCAL_M, LOCAL_N )
      DOUBLE PRECISION, intent (inout) :: RHS_PAR(LOCAL_M, LOCAL_N_RHS)
      INTEGER, intent (out) :: IERR
      INTEGER              :: DESCB_PAR( 9 )
      IERR = 0
      CALL DESCINIT( DESCB_PAR, SIZE_ROOT, 
     &      NRHS, MBLOCK, NBLOCK, 0, 0,
     &      CNTXT_PAR, LOCAL_M, IERR )
            IF (IERR.NE.0) THEN
              WRITE(*,*) 'After DESCINIT, IERR = ', IERR
              CALL MUMPS_ABORT()
            END IF
      IF ( LDLT .eq. 0 .OR. LDLT .eq. 2 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          CALL pdgetrs('N',SIZE_ROOT,NRHS,A,1,1,DESCA_PAR,IPIV,
     &      RHS_PAR,1,1,DESCB_PAR,IERR)
        ELSE
          CALL pdgetrs('T',SIZE_ROOT,NRHS,A,1,1,DESCA_PAR,IPIV,
     &      RHS_PAR, 1, 1, DESCB_PAR,IERR)
        END IF
      ELSE
        CALL pdpotrs( 'L', SIZE_ROOT, NRHS, A, 1, 1, DESCA_PAR,
     &    RHS_PAR, 1, 1, DESCB_PAR, IERR )
      END IF
      IF ( IERR .LT. 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root'
        CALL MUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_2D_BCYCLIC
