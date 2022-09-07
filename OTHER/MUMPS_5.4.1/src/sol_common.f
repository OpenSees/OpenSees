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
      SUBROUTINE MUMPS_SOL_GET_NPIV_LIELL_IPOS ( ISTEP, KEEP,
     &           NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: ISTEP, LIW, KEEP(500), N
      INTEGER, INTENT(IN)  :: IW( LIW )
      INTEGER, INTENT(IN)  :: STEP( N ), PTRIST( KEEP(28) )
      INTEGER, INTENT(OUT) :: NPIV, LIELL, IPOS
      INCLUDE 'mumps_headers.h'
      INTEGER :: SROOT
      IF (KEEP(38) .NE. 0) THEN
        SROOT = STEP(KEEP(38))
      ELSE IF (KEEP(20) .NE. 0) THEN
        SROOT = STEP(KEEP(20))
      ELSE
        SROOT = 0
      ENDIF
      IPOS = PTRIST(ISTEP)
      IF (IPOS .LE. 0) THEN
        WRITE(*,*) "Internal error 1 in MUMPS_SOL_GET_NPIV_LIELL_IPOS",
     &  ISTEP 
        CALL MUMPS_ABORT()
      ENDIF
      NPIV = IW(IPOS+3+KEEP(IXSZ))
      IF ( ISTEP.EQ.SROOT ) THEN
             IPOS = PTRIST(ISTEP)
             LIELL = IW(IPOS+3+KEEP(IXSZ))
             NPIV = LIELL
             IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
      ELSE
             IPOS = PTRIST(ISTEP) + 2+ KEEP(IXSZ)
             LIELL = IW(IPOS-2)+IW(IPOS+1)
             IPOS= IPOS+1
             NPIV = IW(IPOS)
             IPOS= IPOS+1
             IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +KEEP(IXSZ))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_SOL_GET_NPIV_LIELL_IPOS
      SUBROUTINE MUMPS_BUILD_IRHS_loc(MYID_NODES, NSLAVES, N,
     &           PTRIST, KEEP,KEEP8, IW, LIW, STEP, PROCNODE_STEPS,
     &           IRHS_loc, ROW_OR_COL_INDICES)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: KEEP(500)
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      INTEGER, INTENT(IN)    :: MYID_NODES, NSLAVES, N, LIW
      INTEGER, INTENT(IN)    :: PTRIST(KEEP(28))
      INTEGER, INTENT(IN)    :: IW(LIW), STEP(N)
      INTEGER, INTENT(IN)    :: PROCNODE_STEPS(KEEP(28))
      INTEGER, INTENT(OUT)   :: IRHS_loc(KEEP(89))
      INTEGER, INTENT(IN)    :: ROW_OR_COL_INDICES
      INTEGER :: ISTEP
      INTEGER :: NPIV, LIELL, IPOS
      INTEGER :: IIRHS_loc 
      INTEGER :: J1 
      INCLUDE 'mumps_headers.h'
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      IIRHS_loc = 0
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &                                   KEEP(199))) THEN
          CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS ( ISTEP, KEEP,
     &           NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
          IF ( ROW_OR_COL_INDICES .EQ. 0 .OR. KEEP(50).NE.0 ) THEN
            J1 = IPOS + 1
          ELSE IF (ROW_OR_COL_INDICES .EQ. 1 ) THEN
            J1 = IPOS + LIELL + 1
          ELSE
            WRITE(*,*) "Internal error 1 in MUMPS_BUILD_IRHS_loc",
     &      ROW_OR_COL_INDICES
            CALL MUMPS_ABORT()
          ENDIF
          IF (IIRHS_loc+NPIV .GT. KEEP(89)) THEN
            WRITE(*,*) "Internal error 2 in MUMPS_BUILD_IRHS_loc",
     &      IIRHS_loc, KEEP(89)
            CALL MUMPS_ABORT()
          ENDIF
          IRHS_loc(IIRHS_loc+1:IIRHS_loc+NPIV)=IW(J1:J1+NPIV-1)
          IIRHS_loc=IIRHS_loc+NPIV
        ENDIF
      ENDDO
      IF (IIRHS_loc .NE. KEEP(89)) THEN
        WRITE(*,*) "Internal error 3 in MUMPS_BUILD_IRHS_loc",
     &  IIRHS_loc, KEEP(89)
        CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_BUILD_IRHS_loc
      SUBROUTINE MUMPS_SOL_RHSMAPINFO( N, Nloc_RHS, INFO23,
     &                           IRHS_loc, MAP_RHS_loc,
     &                           POSINRHSCOMP_FWD,
     &                           NSLAVES, MYID_NODES, COMM_NODES,
     &                           ICNTL, INFO )
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, Nloc_RHS
      INTEGER, INTENT(IN)    :: INFO23
      INTEGER, INTENT(IN)    :: IRHS_loc   (max(1,Nloc_RHS))     
      INTEGER, INTENT(OUT)   :: MAP_RHS_loc(max(1,Nloc_RHS)) 
      INTEGER, INTENT(IN)    :: POSINRHSCOMP_FWD (N)
      INTEGER, INTENT(IN)    :: NSLAVES, MYID_NODES, COMM_NODES
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INCLUDE 'mpif.h'
      INTEGER :: I, NFS_LOC, NFS_TOT, IERR_MPI, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: GLOBAL_MAPPING 
      ALLOCATE(GLOBAL_MAPPING(N), stat=allocok)
      IF (allocok .GT. 0) THEN
        INFO(1)=-13
        INFO(2)= N
      ENDIF
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, allocok, 1,
     &                   MPI_INTEGER, MPI_SUM,
     &                   COMM_NODES, IERR_MPI)
      IF (allocok .NE. 0) RETURN
      NFS_LOC = 0
      NFS_TOT = 0
      DO I = 1, N
        IF (POSINRHSCOMP_FWD(I) .LE. 0) THEN
          GLOBAL_MAPPING(I) = 0
        ELSE
          GLOBAL_MAPPING(I) = MYID_NODES
          NFS_LOC = NFS_LOC + 1
        ENDIF
      ENDDO
      IF (NFS_LOC .NE. INFO23) THEN
        WRITE(*,*) "Internal error 1 in MUMPS_SOL_RHSMAPINFO",
     &  NFS_LOC, INFO23
        CALL MUMPS_ABORT()
      ENDIF
      CALL MPI_ALLREDUCE(NFS_LOC, NFS_TOT, 1, MPI_INTEGER,
     &                   MPI_SUM, COMM_NODES, IERR_MPI)
      IF (NFS_tot .NE. N) THEN
        WRITE(*,*) "Internal error 1 in MUMPS_SOL_RHSMAPINFO",
     &  NFS_LOC, NFS_TOT, N
        CALL MUMPS_ABORT()
      ENDIF
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, GLOBAL_MAPPING, N, MPI_INTEGER,
     &                   MPI_SUM, COMM_NODES, IERR_MPI)
      DO I = 1, Nloc_RHS
        IF (IRHS_loc(I) .GE.1 .AND. IRHS_loc(I) .LE. N) THEN
          MAP_RHS_loc(I) = GLOBAL_MAPPING(IRHS_loc(I))
        ELSE
          MAP_RHS_loc(I) = -87878787
        ENDIF
      ENDDO
      DEALLOCATE(GLOBAL_MAPPING)
      RETURN
      END SUBROUTINE MUMPS_SOL_RHSMAPINFO
