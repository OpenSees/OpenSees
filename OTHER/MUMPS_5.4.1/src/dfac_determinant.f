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
      SUBROUTINE DMUMPS_UPDATEDETER(PIV, DETER, NEXP)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: PIV
      DOUBLE PRECISION, intent(inout) :: DETER
      INTEGER, intent(inout) :: NEXP
      DETER=DETER*fraction(PIV)
      NEXP=NEXP+exponent(PIV)+exponent(DETER)
      DETER=fraction(DETER)
      RETURN
      END SUBROUTINE DMUMPS_UPDATEDETER
      SUBROUTINE DMUMPS_UPDATEDETER_SCALING(PIV, DETER, NEXP)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: PIV
      DOUBLE PRECISION, intent(inout) :: DETER
      INTEGER, intent(inout) :: NEXP
      DETER=DETER*fraction(PIV)
      NEXP=NEXP+exponent(PIV)+exponent(DETER)
      DETER=fraction(DETER)
      RETURN
      END SUBROUTINE DMUMPS_UPDATEDETER_SCALING
      SUBROUTINE DMUMPS_GETDETER2D(BLOCK_SIZE,IPIV,
     &                      MYROW, MYCOL, NPROW, NPCOL,
     &                      A, LOCAL_M, LOCAL_N, N, MYID,
     &                      DETER,NEXP,SYM)
      IMPLICIT NONE
      INTEGER, intent (in)    :: SYM
      INTEGER, intent (inout) :: NEXP
      DOUBLE PRECISION, intent (inout) :: DETER
      INTEGER, intent (in)    :: BLOCK_SIZE, NPROW, NPCOL,
     &                           LOCAL_M, LOCAL_N, N
      INTEGER, intent (in)    :: MYROW, MYCOL, MYID, IPIV(LOCAL_M)
      DOUBLE PRECISION, intent(in) :: A(*)
      INTEGER  I,IMX,DI,NBLOCK,IBLOCK,ILOC,JLOC,
     &         ROW_PROC,COL_PROC, K
      DI = LOCAL_M + 1 
      NBLOCK = ( N - 1 ) / BLOCK_SIZE
      DO IBLOCK = 0, NBLOCK
        ROW_PROC = mod( IBLOCK, NPROW ) 
        IF ( MYROW.EQ.ROW_PROC ) THEN
          COL_PROC = mod( IBLOCK, NPCOL )
          IF ( MYCOL.EQ.COL_PROC ) THEN
            ILOC = ( IBLOCK / NPROW ) * BLOCK_SIZE
            JLOC = ( IBLOCK / NPCOL ) * BLOCK_SIZE
            I   =   ILOC + JLOC *  LOCAL_M + 1
            IMX = min(ILOC+BLOCK_SIZE,LOCAL_M)
     &            + (min(JLOC+BLOCK_SIZE,LOCAL_N)-1)*LOCAL_M
     &            + 1
            K=1
            DO WHILE ( I .LT. IMX )
              CALL DMUMPS_UPDATEDETER(A(I),DETER,NEXP)
              IF (SYM.EQ.1) THEN
               CALL DMUMPS_UPDATEDETER(A(I),DETER,NEXP)
              ENDIF
              IF (SYM.NE.1) THEN 
                IF (IPIV(ILOC+K) .NE. IBLOCK*BLOCK_SIZE+K) THEN
                  DETER = -DETER
                ENDIF
              ENDIF
              K = K + 1
              I = I + DI
            END DO
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_GETDETER2D
      SUBROUTINE DMUMPS_DETER_REDUCTION(
     &           COMM, DETER_IN, NEXP_IN,
     &           DETER_OUT, NEXP_OUT, NPROCS)
      IMPLICIT NONE
      INTEGER, intent(in) :: COMM, NPROCS
      DOUBLE PRECISION, intent(in) :: DETER_IN
      INTEGER,intent(in) :: NEXP_IN
      DOUBLE PRECISION,intent(out):: DETER_OUT
      INTEGER,intent(out):: NEXP_OUT
      INTEGER            :: IERR_MPI
      EXTERNAL DMUMPS_DETERREDUCE_FUNC
      INTEGER TWO_SCALARS_TYPE, DETERREDUCE_OP
      DOUBLE PRECISION :: INV(2)
      DOUBLE PRECISION :: OUTV(2)
      INCLUDE 'mpif.h'
      IF (NPROCS .EQ. 1) THEN
        DETER_OUT = DETER_IN
        NEXP_OUT  = NEXP_IN
        RETURN
      ENDIF
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_DOUBLE_PRECISION,
     &                         TWO_SCALARS_TYPE,
     &                         IERR_MPI)
      CALL MPI_TYPE_COMMIT(TWO_SCALARS_TYPE, IERR_MPI)
      CALL MPI_OP_CREATE(DMUMPS_DETERREDUCE_FUNC,  
     &                   .TRUE.,             
     &                   DETERREDUCE_OP,     
     &                   IERR_MPI)
      INV(1)=DETER_IN
      INV(2)=dble(NEXP_IN)
      CALL MPI_ALLREDUCE( INV, OUTV, 1, TWO_SCALARS_TYPE,
     &                    DETERREDUCE_OP, COMM, IERR_MPI)
      CALL MPI_OP_FREE(DETERREDUCE_OP, IERR_MPI)
      CALL MPI_TYPE_FREE(TWO_SCALARS_TYPE, IERR_MPI)
      DETER_OUT = OUTV(1)
      NEXP_OUT  = int(OUTV(2))
      RETURN
      END SUBROUTINE DMUMPS_DETER_REDUCTION
      SUBROUTINE DMUMPS_DETERREDUCE_FUNC(INV, INOUTV, NEL, DATATYPE)
      IMPLICIT NONE
#if defined(WORKAROUNDINTELILP64MPI2INTEGER)
      INTEGER(4), INTENT(IN)    :: NEL, DATATYPE
#else
      INTEGER, INTENT(IN)    :: NEL, DATATYPE
#endif
      DOUBLE PRECISION, INTENT(IN)    :: INV    ( 2 * NEL )
      DOUBLE PRECISION, INTENT(INOUT) :: INOUTV ( 2 * NEL )
      INTEGER I, TMPEXPIN, TMPEXPINOUT
      DO I = 1, NEL
        TMPEXPIN    = int(INV   (I*2))
        TMPEXPINOUT = int(INOUTV(I*2))
        CALL DMUMPS_UPDATEDETER(INV(I*2-1), 
     &                          INOUTV(I*2-1), 
     &                          TMPEXPINOUT)   
        TMPEXPINOUT = TMPEXPINOUT + TMPEXPIN
        INOUTV(I*2) = dble(TMPEXPINOUT)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_DETERREDUCE_FUNC
      SUBROUTINE DMUMPS_DETER_SQUARE(DETER, NEXP)
      IMPLICIT NONE
      INTEGER, intent (inout) :: NEXP
      DOUBLE PRECISION, intent (inout) :: DETER
      DETER=DETER*DETER
      NEXP=NEXP+NEXP
      RETURN
      END SUBROUTINE DMUMPS_DETER_SQUARE
      SUBROUTINE DMUMPS_DETER_SCALING_INVERSE(DETER, NEXP)
      IMPLICIT NONE
      INTEGER, intent (inout) :: NEXP
      DOUBLE PRECISION, intent (inout) :: DETER
      DETER=1.0D0/DETER
      NEXP=-NEXP
      RETURN
      END SUBROUTINE DMUMPS_DETER_SCALING_INVERSE
      SUBROUTINE DMUMPS_DETER_SIGN_PERM(DETER, N, VISITED, PERM)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(inout) :: DETER
      INTEGER, intent(in)    :: N
      INTEGER, intent(inout) :: VISITED(N)
      INTEGER, intent(in)    :: PERM(N)
      INTEGER I, J, K
      K = 0
      DO I = 1, N
        IF (VISITED(I) .GT. N) THEN
          VISITED(I)=VISITED(I)-N-N-1
          CYCLE
        ENDIF
        J = PERM(I)
        DO WHILE (J.NE.I)
          VISITED(J) = VISITED(J) + N + N + 1
          K = K + 1
          J = PERM(J)
        ENDDO
      ENDDO
      IF (mod(K,2).EQ.1) THEN
        DETER = -DETER
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_DETER_SIGN_PERM
      SUBROUTINE DMUMPS_PAR_ROOT_MINMAX_PIV_UPD (
     &                      BLOCK_SIZE,IPIV,
     &                      MYROW, MYCOL, NPROW, NPCOL,
     &                      A, LOCAL_M, LOCAL_N, N, MYID,
     &                      DKEEP, KEEP, SYM)
       USE DMUMPS_FAC_FRONT_AUX_M, 
     &      ONLY : DMUMPS_UPDATE_MINMAX_PIVOT 
      IMPLICIT NONE
      INTEGER, intent (in)    :: BLOCK_SIZE, NPROW, NPCOL,
     &                           LOCAL_M, LOCAL_N, N, SYM
      INTEGER, intent (in)    :: MYROW, MYCOL, MYID, IPIV(LOCAL_M)
      DOUBLE PRECISION, intent(in) :: A(*)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER, INTENT(IN) :: KEEP(500)
      INTEGER  I,IMX,DI,NBLOCK,IBLOCK,ILOC,JLOC,
     &         ROW_PROC,COL_PROC, K
      DOUBLE PRECISION :: ABSPIVOT
      DI = LOCAL_M + 1 
      NBLOCK = ( N - 1 ) / BLOCK_SIZE
      DO IBLOCK = 0, NBLOCK
        ROW_PROC = mod( IBLOCK, NPROW ) 
        IF ( MYROW.EQ.ROW_PROC ) THEN
          COL_PROC = mod( IBLOCK, NPCOL )
          IF ( MYCOL.EQ.COL_PROC ) THEN
            ILOC = ( IBLOCK / NPROW ) * BLOCK_SIZE
            JLOC = ( IBLOCK / NPCOL ) * BLOCK_SIZE
            I   =   ILOC + JLOC *  LOCAL_M + 1
            IMX = min(ILOC+BLOCK_SIZE,LOCAL_M)
     &            + (min(JLOC+BLOCK_SIZE,LOCAL_N)-1)*LOCAL_M
     &            + 1
            K=1
            DO WHILE ( I .LT. IMX )
             IF (SYM.NE.1) THEN
              ABSPIVOT =  abs(A(I))
             ELSE
              ABSPIVOT =  abs(A(I)*A(I))
             ENDIF
             CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( ABSPIVOT,
     &             DKEEP, KEEP, .FALSE.)
              K = K + 1
              I = I + DI
            END DO
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_PAR_ROOT_MINMAX_PIV_UPD
