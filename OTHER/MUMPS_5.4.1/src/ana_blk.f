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
      SUBROUTINE MUMPS_AB_FREE_LMAT ( LMAT )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      TYPE(LMATRIX_T) :: LMAT
      INTEGER :: J
      IF (associated(LMAT%COL)) THEN
        DO J = 1,LMAT%NBCOL
         IF (associated(LMAT%COL(J)%IRN)) THEN
           DEALLOCATE(LMAT%COL(J)%IRN)
           NULLIFY(LMAT%COL(J)%IRN)
         ENDIF
        ENDDO
        DEALLOCATE(LMAT%COL)
        NULLIFY(LMAT%COL)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_AB_FREE_LMAT
      SUBROUTINE MUMPS_AB_FREE_GCOMP ( GCOMP )
      USE MUMPS_ANA_BLK_M, ONLY : COMPACT_GRAPH_T
      IMPLICIT NONE
      TYPE(COMPACT_GRAPH_T) :: GCOMP
      IF (associated(GCOMP%IPE)) THEN
          DEALLOCATE(GCOMP%IPE)
          NULLIFY(GCOMP%IPE)
      ENDIF
      IF (associated(GCOMP%ADJ)) THEN
          DEALLOCATE(GCOMP%ADJ)
          NULLIFY(GCOMP%ADJ)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_AB_FREE_GCOMP
      SUBROUTINE MUMPS_AB_COMPUTE_SIZEOFBLOCK (
     &          NBLK, NDOF, BLKPTR, BLKVAR,
     &          SIZEOFBLOCKS, DOF2BLOCK )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NBLK, NDOF
      INTEGER, INTENT(IN) :: BLKPTR(NBLK+1), BLKVAR(NDOF)
      INTEGER, INTENT(OUT):: SIZEOFBLOCKS(NBLK), DOF2BLOCK(NDOF)
      INTEGER :: IB, I, IDOF
      DO IB=1, NBLK
        SIZEOFBLOCKS(IB)= BLKPTR(IB+1)-BLKPTR(IB)
        DO I=BLKPTR(IB), BLKPTR(IB+1)-1
          IDOF = BLKVAR(I)
          DOF2BLOCK(IDOF) = IB
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_COMPUTE_SIZEOFBLOCK 
      SUBROUTINE MUMPS_AB_COORD_TO_LMAT (  MYID, 
     &     NBLK, NDOF, NNZ, IRN, JCN,
     &     DOF2BLOCK, 
     &     IFLAG, IERROR, LP, LPOK,
     &     LMAT)
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: MYID, NBLK, NDOF
      INTEGER(8), INTENT(IN) :: NNZ
      INTEGER, INTENT(IN) :: IRN(max(1_8,NNZ)), JCN(max(1_8,NNZ))
      INTEGER, INTENT(IN) :: DOF2BLOCK(NDOF)
      INTEGER             :: LP, IFLAG, IERROR
      LOGICAL, INTENT(IN) :: LPOK
      TYPE(LMATRIX_T)     :: LMAT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: FLAG
      INTEGER :: allocok
      INTEGER :: I, J, JJB, IIB, IB, JB, NB, PT
      INTEGER(8) :: I8
      LMAT%NBCOL = NBLK
      LMAT%NZL   = 0_8
      ALLOCATE(LMAT%COL(NBLK),FLAG(NBLK), STAT=allocok)
      IF (allocok.NE.0) THEN
           IFLAG = -7
           IERROR = 2*NBLK
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LMAT%COL"
           END IF
           RETURN
      ENDIF
      DO IB=1,NBLK
       LMAT%COL(IB)%NBINCOL = 0
       FLAG(IB)             = 0
      ENDDO
      IERROR = 0 
      DO I8=1, NNZ
         I = IRN(I8)
         J = JCN(I8)
         IF ( (I.GT.NDOF).OR.(J.GT.NDOF).OR.(I.LT.1)
     &                     .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
         ELSE
          IB  = DOF2BLOCK(I)
          JB  = DOF2BLOCK(J)
          JJB = min(IB,JB)
          IF (IB.NE.JB) THEN
           LMAT%NZL = LMAT%NZL+1_8
           LMAT%COL(JJB)%NBINCOL =  LMAT%COL(JJB)%NBINCOL + 1
          ENDIF
         ENDIF
      ENDDO
      IF (IERROR.GE.1) THEN
         IF (mod(IFLAG,2) .EQ. 0) IFLAG = IFLAG+1
      ENDIF
      DO JB=1,NBLK
       NB =  LMAT%COL(JB)%NBINCOL
       IF (NB.GT.0) THEN
        ALLOCATE(LMAT%COL(JB)%IRN(NB), STAT=allocok)
        IF (allocok.NE.0) THEN
           IFLAG  = -7
           IERROR = NBLK
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LMAT%COL"
           END IF
           RETURN
        ENDIF
       ENDIF
      ENDDO
      DO I8=1, NNZ
         I = IRN(I8)
         J = JCN(I8)
         IF ( (I.LE.NDOF).AND.(J.LE.NDOF).AND.(I.GE.1)
     &                     .AND.(J.GE.1)) THEN
          IB  = DOF2BLOCK(I)
          JB  = DOF2BLOCK(J)
          JJB = min(IB,JB)
          IIB = max(IB,JB)
          IF (IIB.NE.JJB) THEN
             PT        = FLAG(JJB)+1
             FLAG(JJB) = PT
             LMAT%COL(JJB)%IRN(PT) = IIB
          ENDIF
         ENDIF
      ENDDO
      CALL MUMPS_AB_LOCALCLEAN_LMAT ( MYID,
     &     NBLK, LMAT, FLAG(1), IFLAG, IERROR, LP, LPOK
     & )
      DEALLOCATE(FLAG)
      RETURN
      END SUBROUTINE MUMPS_AB_COORD_TO_LMAT
      SUBROUTINE MUMPS_AB_LOCALCLEAN_LMAT (  MYID,
     &     NBLK, LMAT, FLAG, IFLAG, IERROR, LP, LPOK
     &     )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MYID, NBLK, LP
      LOGICAL, INTENT(IN) :: LPOK
      INTEGER, INTENT(OUT) :: FLAG(NBLK)
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      TYPE(LMATRIX_T), INTENT(INOUT)  :: LMAT
      INTEGER, POINTER, DIMENSION(:)     :: PTCLEAN
      INTEGER :: allocok, IB, JB, NB
      DO JB=1, NBLK
        FLAG(JB) = 0
      ENDDO
      LMAT%NZL = 0_8
      DO JB=1, NBLK
       IF ( LMAT%COL(JB)%NBINCOL.EQ.0) CYCLE
       NB = 0
       DO IB=1,  LMAT%COL(JB)%NBINCOL
        IF (FLAG(LMAT%COL(JB)%IRN(IB)).EQ.JB) THEN
         LMAT%COL(JB)%IRN(IB)=0
        ELSE
         NB = NB+1
         LMAT%NZL = LMAT%NZL+1_8
         FLAG(LMAT%COL(JB)%IRN(IB)) = JB
        ENDIF
       ENDDO
       IF (NB.GT.0) THEN
         ALLOCATE(PTCLEAN(NB), STAT=allocok)
         IF (allocok.NE.0) THEN
           IFLAG     = -7
           IERROR    = NB
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate PTCLEAN of size", 
     &                     IERROR
           END IF
           RETURN
         ENDIF
         NB=0
         DO IB=1,  LMAT%COL(JB)%NBINCOL
          IF (LMAT%COL(JB)%IRN(IB).NE.0) THEN
           NB = NB+1
           PTCLEAN(NB)=LMAT%COL(JB)%IRN(IB)
          ENDIF
         ENDDO
         LMAT%COL(JB)%NBINCOL  = NB
         deallocate(LMAT%COL(JB)%IRN)
         LMAT%COL(JB)%IRN => PTCLEAN
         NULLIFY(PTCLEAN)
       ELSE
          deallocate(LMAT%COL(JB)%IRN)
          NULLIFY(LMAT%COL(JB)%IRN)
       ENDIF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_LOCALCLEAN_LMAT
      SUBROUTINE MUMPS_AB_LMAT_TO_LUMAT( 
     &     LMAT, LUMAT, INFO, ICNTL )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T, COMPACT_GRAPH_T
      IMPLICIT NONE
      TYPE(LMATRIX_T)        :: LMAT, LUMAT
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER    :: IB, IIB, JB, allocok, LP, MPG, NB, IERR
      LOGICAL LPOK, PROKG
      LP    = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MPG   = ICNTL( 3 )
      PROKG = ( MPG .GT. 0 .and. (ICNTL(4).GE.2) ) 
      LUMAT%NBCOL  = LMAT%NBCOL
      LUMAT%NZL    = 2_8*LMAT%NZL
      ALLOCATE( LUMAT%COL(LMAT%NBCOL),STAT=allocok)
      IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           INFO( 2 ) = LMAT%NBCOL
           IF ( LPOK ) THEN
             WRITE(LP, *) " ERROR allocating LUMAT%COL "
           END IF
           RETURN
      ENDIF
      DO JB=1,  LMAT%NBCOL
         LUMAT%COL(JB)%NBINCOL = LMAT%COL(JB)%NBINCOL
      ENDDO
      DO JB=1,  LMAT%NBCOL
        DO IB=1, LMAT%COL(JB)%NBINCOL
          IIB=LMAT%COL(JB)%IRN(IB)
          LUMAT%COL(IIB)%NBINCOL = LUMAT%COL(IIB)%NBINCOL + 1
        ENDDO
      ENDDO
      DO JB=1,  LMAT%NBCOL
        NB = LUMAT%COL(JB)%NBINCOL
        ALLOCATE(LUMAT%COL(JB)%IRN(NB), STAT=IERR)
        IF (IERR.NE.0) THEN
           INFO(1)  = -7
           INFO(2)  = NB
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocating columns of LUMAT"
           END IF
           RETURN
        ENDIF
      ENDDO
      DO JB=1, LMAT%NBCOL
        LUMAT%COL(JB)%NBINCOL = 0
      ENDDO
      DO JB=1, LMAT%NBCOL
        DO IB=1, LMAT%COL(JB)%NBINCOL
          IIB=LMAT%COL(JB)%IRN(IB)
          NB = LUMAT%COL(JB)%NBINCOL+1
          LUMAT%COL(JB)%NBINCOL = NB
          LUMAT%COL(JB)%IRN(NB)= IIB
          NB = LUMAT%COL(IIB)%NBINCOL+1
          LUMAT%COL(IIB)%NBINCOL = NB
          LUMAT%COL(IIB)%IRN(NB)= JB
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_LMAT_TO_LUMAT
      SUBROUTINE MUMPS_AB_PRINT_LMATRIX (LMAT, MYID, LP)
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      TYPE(LMATRIX_T), INTENT(IN) :: LMAT
      INTEGER, INTENT(IN)         :: MYID, LP
      INTEGER :: JB
      write(LP,*) MYID, " ... LMATRIX  %NBCOL, %NZL= ", 
     &                  LMAT%NBCOL, LMAT%NZL
      IF (LMAT%NBCOL.GE.0.AND.associated(LMAT%COL)) THEN
       DO JB=1,  LMAT%NBCOL
        IF (LMAT%COL(JB)%NBINCOL.GT.0) THEN
         WRITE(LP,*) MYID, " ... Column=", JB , " nb entries =", 
     &    LMAT%COL(JB)%NBINCOL, " List of entries:",
     &    LMAT%COL(JB)%IRN(1:LMAT%COL(JB)%NBINCOL)
        ENDIF
       ENDDO
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_AB_PRINT_LMATRIX
      SUBROUTINE MUMPS_AB_LMAT_TO_CLEAN_G( MYID, UNFOLD, 
     &     READY_FOR_ANA_F,
     &     LMAT, GCOMP, INFO, ICNTL )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T, COMPACT_GRAPH_T
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: MYID
      LOGICAL, INTENT(IN)    :: UNFOLD, READY_FOR_ANA_F
      TYPE(LMATRIX_T)        :: LMAT
      TYPE(COMPACT_GRAPH_T)  :: GCOMP
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER    :: IB, IIB, JJB, allocok, LP, MPG
      INTEGER(8) :: JPOS, SIZEGCOMPALLOCATED
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IQ
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: WORK
      INTEGER(8) :: IFIRST, ILAST
      INTEGER :: L
#endif
      LOGICAL LPOK, PROKG
      LP    = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MPG   = ICNTL( 3 )
      PROKG = ( MPG .GT. 0 .and. (ICNTL(4).GE.2) ) 
      GCOMP%NG  = LMAT%NBCOL
      IF (UNFOLD) THEN 
       GCOMP%NZG = 2_8*LMAT%NZL
       SIZEGCOMPALLOCATED = GCOMP%NZG + int(GCOMP%NG,8)+1_8
      ELSE IF (READY_FOR_ANA_F) THEN
       GCOMP%NZG = LMAT%NZL
       SIZEGCOMPALLOCATED = GCOMP%NZG + int(GCOMP%NG,8)+1_8
      ELSE
       GCOMP%NZG = LMAT%NZL
       SIZEGCOMPALLOCATED = GCOMP%NZG
      ENDIF
      GCOMP%SIZEADJALLOCATED= SIZEGCOMPALLOCATED
      ALLOCATE( GCOMP%ADJ(SIZEGCOMPALLOCATED), 
     &          GCOMP%IPE(GCOMP%NG+1), 
     &          IQ(GCOMP%NG),STAT=allocok)
      IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           CALL MUMPS_SET_IERROR(
     &        GCOMP%NZG + 3_8*int(GCOMP%NG,8)+1_8, INFO(2))
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocating graph in",
     &                     " MUMPS_AB_LMAT_TO_CLEAN_G"
           END IF
           RETURN
      ENDIF
      DO JJB=1, GCOMP%NG
         IQ(JJB)=0_8
      ENDDO
      IF (UNFOLD) THEN 
       DO JJB=1, GCOMP%NG
        DO IB=1, LMAT%COL(JJB)%NBINCOL
          IIB=LMAT%COL(JJB)%IRN(IB)
          IQ(JJB)=IQ(JJB)+1
          IQ(IIB)=IQ(IIB)+1
        ENDDO
       ENDDO
      ELSE
       DO JJB=1, GCOMP%NG
         IQ(JJB) = LMAT%COL(JJB)%NBINCOL
       ENDDO
      ENDIF
      GCOMP%IPE(1) = 1_8
      DO JJB=1, GCOMP%NG
        GCOMP%IPE(JJB+1) = GCOMP%IPE(JJB)+IQ(JJB)
      ENDDO
      IF (UNFOLD) THEN 
       DO JJB=1, GCOMP%NG
        IQ(JJB)= GCOMP%IPE(JJB)
       ENDDO
       DO JJB=1, GCOMP%NG
        DO IB=1, LMAT%COL(JJB)%NBINCOL
          IIB=LMAT%COL(JJB)%IRN(IB)
          GCOMP%ADJ(IQ(IIB))= JJB
          IQ(IIB)           = IQ(IIB)+1_8
          GCOMP%ADJ(IQ(JJB))= IIB
          IQ(JJB)           = IQ(JJB)+1_8
        ENDDO
       ENDDO
      ELSE
       DO JJB=1, GCOMP%NG
        JPOS =  GCOMP%IPE(JJB)
        DO IB=1, LMAT%COL(JJB)%NBINCOL
          IIB=LMAT%COL(JJB)%IRN(IB)
          GCOMP%ADJ(JPOS)= IIB
          JPOS           = JPOS+1_8
       ENDDO
      ENDDO
      ENDIF
      DEALLOCATE(IQ)
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      IF (.NOT.READY_FOR_ANA_F) THEN
        ALLOCATE(WORK(0:GCOMP%NG),stat=allocok)
        IF (allocok.NE.0) THEN
          INFO( 1 ) = -7
          INFO( 2 ) = GCOMP%NG
          IF ( LPOK ) THEN
             WRITE(LP, *) " ERROR allocating WORK in",
     &                    " MUMPS_AB_LMAT_TO_CLEAN_G"
          END IF
          RETURN
        ENDIF
        DO JJB=1, GCOMP%NG
          IFIRST = GCOMP%IPE(JJB)
          ILAST= GCOMP%IPE(JJB+1)-1
          L = int(ILAST-IFIRST+1)
          IF ( L .GE. 2 ) THEN 
            IF (L .GE. GCOMP%NG ) THEN
              WRITE(*,*) " Internal error in MUMPS_AB_LMAT_TO_CLEAN_G",
     &        L, GCOMP%NG
              CALL MUMPS_ABORT()
            ENDIF
            CALL MUMPS_MERGESORT( L,
     &      GCOMP%ADJ(IFIRST:ILAST), WORK(0:L+1) )
            CALL MUMPS_MERGESWAP1( L,
     &      WORK(0:L+1), GCOMP%ADJ(IFIRST:ILAST) )
          ENDIF
        ENDDO
        DEALLOCATE(WORK)
      ENDIF
#endif
      RETURN
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      CONTAINS
      SUBROUTINE MUMPS_MERGESORT(N, K, L)
      INTEGER    :: N
      INTEGER    :: K(:), L(0:)
      INTEGER    :: P, Q, S, T
      CONTINUE
      L(0) = 1
      T = N + 1
      DO  P = 1,N - 1
         IF (K(P) <= K(P+1)) THEN
            L(P) = P + 1
         ELSE
            L(T) = - (P+1)
            T = P
       END IF
      END DO
      L(T) = 0
      L(N) = 0
      IF (L(N+1) == 0) THEN
         RETURN 
      ELSE
         L(N+1) = iabs(L(N+1))
      END IF
 200  CONTINUE
      S = 0
      T = N+1
      P = L(S)
      Q = L(T)
      IF(Q .EQ. 0) RETURN
 300  CONTINUE
      IF(K(P) .GT. K(Q)) GOTO 600 
      CONTINUE
      L(S) = sign(P,L(S))
      S = P
      P = L(P)
      IF (P .GT. 0) GOTO 300
      CONTINUE
      L(S) = Q
      S = T
      DO
         T = Q
         Q = L(Q)
         IF (Q .LE. 0) EXIT
      END DO
      GOTO 800
 600  CONTINUE
      L(S) = sign(Q, L(S))
      S = Q
      Q = L(Q)
      IF (Q .GT. 0) GOTO 300
      CONTINUE
      L(S) = P
      S = T
      DO
         T = P
         P = L(P)
         IF (P .LE. 0) EXIT
      END DO
 800  CONTINUE
      P = -P
      Q = -Q
      IF(Q.EQ.0) THEN
         L(S) = sign(P, L(S))
         L(T) = 0
         GOTO 200
      END IF
      GOTO 300
      END SUBROUTINE MUMPS_MERGESORT
      SUBROUTINE MUMPS_MERGESWAP1(N, L, A)
      INTEGER   :: I, LP, ISWAP, N
      INTEGER   :: L(0:), A(:)
      LP = L(0)
      I  = 1
      DO 
         IF ((LP==0).OR.(I>N)) EXIT
         DO 
            IF (LP >= I) EXIT
            LP = L(LP)
         END DO
         ISWAP    = A(LP)
         A(LP)   = A(I)
         A(I)    = ISWAP
         ISWAP    = L(LP)
         L(LP) = L(I)
         L(I)  = LP
         LP = ISWAP 
         I  = I + 1
      ENDDO
      END SUBROUTINE MUMPS_MERGESWAP1
#endif
      END SUBROUTINE MUMPS_AB_LMAT_TO_CLEAN_G
      SUBROUTINE MUMPS_AB_COL_DISTRIBUTION ( OPTION,
     &     INFO, ICNTL, COMM, NBLK, MYID, NPROCS,
     &     LMAT, MAPCOL )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER, INTENT(IN) :: OPTION, NBLK
      INTEGER, INTENT(IN) :: ICNTL(60), COMM, MYID, NPROCS
      INTEGER             :: INFO(80)
      TYPE(LMATRIX_T)     :: LMAT
      INTEGER, INTENT(OUT):: MAPCOL(NBLK)
      INTEGER    :: LP, SIZE_NZROW, I
      LOGICAL    :: LPOK
      INTEGER(8) :: NZL, NNZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: NZ_ROW
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      IF (OPTION.EQ.1) THEN
       NNZ        = -9999 
       SIZE_NZROW = 1
      ELSE
       NZL  = LMAT%NZL
       SIZE_NZROW = NBLK
      ENDIF
      ALLOCATE(NZ_ROW(NBLK), STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1)  = -7
           INFO(2)  = SIZE_NZROW
           IF ( LPOK ) THEN
              WRITE(LP, *) 
     &     " ERROR allocate in MUMPS_AB_COL_DISTRIBUTION ", INFO(2)
           END IF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF (INFO(1).LT.0) GOTO 500
      IF (OPTION.NE.1) THEN
        DO I = 1, NBLK
         MAPCOL(I) = LMAT%COL(I)%NBINCOL
        ENDDO
        CALL MPI_ALLREDUCE(MAPCOL(1), NZ_ROW(1), NBLK, 
     &        MPI_INTEGER, MPI_SUM, COMM, IERR)
        CALL MPI_ALLREDUCE(NZL, NNZ, 1, 
     &        MPI_INTEGER8, MPI_SUM, COMM, IERR)
      ENDIF
      CALL MUMPS_AB_COMPUTE_MAPCOL (OPTION, INFO, ICNTL, MYID,
     &   NNZ, NZ_ROW(1), SIZE_NZROW, NBLK, NPROCS, MAPCOL(1))
 500  CONTINUE
      IF (allocated(NZ_ROW)) DEALLOCATE(NZ_ROW)
      RETURN
      END SUBROUTINE MUMPS_AB_COL_DISTRIBUTION
      SUBROUTINE MUMPS_AB_COMPUTE_MAPCOL (OPTION, INFO, ICNTL, 
     &    MYID, NNZ, NZ_ROW, SIZE_NZROW, NBLK,  NPROCS, MAPCOL )
      INTEGER, INTENT(IN) :: OPTION, MYID, SIZE_NZROW, NBLK
      INTEGER, INTENT(IN) :: ICNTL(60), NPROCS
      INTEGER             :: INFO(80)
      INTEGER(8)          :: NNZ
      INTEGER, INTENT(IN) :: NZ_ROW(SIZE_NZROW)
      INTEGER, INTENT(OUT):: MAPCOL(NBLK)
      INTEGER    :: I, J, P, F, LP, IERR
      LOGICAL    :: LPOK
      INTEGER(8) :: SHARE, T
      INTEGER, DIMENSION(:), ALLOCATABLE :: FIRST
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      ALLOCATE(FIRST(NPROCS+1), STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1)  = -7
           INFO(2)  = NPROCS+1
           IF ( LPOK ) THEN
              WRITE(LP, *) 
     &     " ERROR allocate in MUMPS_AB_COL_DISTRIBUTION ", INFO(2)
           END IF
           GOTO 500
      ENDIF
      DO I=1,NPROCS+1
       FIRST(I) = 0
      ENDDO
      IF (OPTION.EQ.1) THEN
       SHARE = int(NBLK/NPROCS,8)
       DO I=1, NPROCS
          FIRST(I) = (I-1)*int(SHARE)+1
       END DO
       FIRST(NPROCS+1)=NBLK+1
      ELSE
       SHARE = (NNZ-1_8)/int(NPROCS,8) + 1_8
         P = 0
         T = 0_8
         F = 1
         DO I=1, NBLK
            T = T+int(NZ_ROW(I),8)
            IF (
     &           (T .GE. SHARE) .OR.
     &           ((NBLK-I).EQ.(NPROCS-P-1)) .OR.
     &           (I.EQ.NBLK)
     &           ) THEN
               P             = P+1
               IF(P.EQ.NPROCS) THEN
                  FIRST(P) = F
                  EXIT
               ELSE
                  FIRST(P) = F
                  F             = I+1
                  T             = 0_8
               END IF
            END IF
            IF ((I.EQ.NBLK).AND.(P.NE.NPROCS)) THEN 
             DO J=P,NPROCS
              FIRST(J) = FIRST(P)
             ENDDO
            ENDIF
         END DO
         FIRST(NPROCS+1) = NBLK+1
      ENDIF
      DO I=1,NPROCS
        DO J=FIRST(I), FIRST(I+1)-1
          MAPCOL(J) = I-1
        ENDDO
      ENDDO
      IF (allocated(FIRST))  DEALLOCATE(FIRST)
 500  CONTINUE
      RETURN
      END SUBROUTINE MUMPS_AB_COMPUTE_MAPCOL
      SUBROUTINE MUMPS_AB_BUILD_DCLEAN_LUMATRIX ( 
     &     MAPCOLonLUMAT, MAPCOL_IN_NSTEPS,
     &     INFO, ICNTL, KEEP, COMM, MYID, NBLK, NPROCS,
     &     LMAT, MAPCOL, SIZEMAPCOL,
     &     STEP, SIZESTEP, 
     &     LUMAT)
      USE MUMPS_ANA_BLK_M 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      LOGICAL, INTENT(IN) :: MAPCOLonLUMAT, MAPCOL_IN_NSTEPS
      INTEGER, INTENT(IN) :: MYID, NPROCS, NBLK, SIZEMAPCOL
      INTEGER, INTENT(IN) :: ICNTL(60), COMM, KEEP(500)
      INTEGER, INTENT(IN) :: SIZESTEP
      INTEGER, INTENT(IN) :: STEP(SIZESTEP)
      INTEGER, INTENT(INOUT) ::  INFO(80)
      TYPE(LMATRIX_T), INTENT(IN)  :: LMAT
      INTEGER, INTENT(INOUT)       :: MAPCOL(SIZEMAPCOL)
      TYPE(LMATRIX_T), INTENT(OUT) :: LUMAT
      INTEGER ::      NBLKloc, IERR, JB, IB, LP, NB, I,
     &                NBRECORDS
      INTEGER(8)  ::  NNZ, NZ_locMAX8, NSEND8, NLOCAL8
      LOGICAL :: LPOK
      INTEGER, ALLOCATABLE, DIMENSION(:) :: WT, WNBINCOL
      INTEGER OPTION
      PARAMETER (OPTION=2)
      NBLKloc = LMAT%NBCOL
      IF (NBLKloc.NE.NBLK) THEN
       write(6,*) "Internal error in MUMPS_AB_BUILD_DCLEAN_LUMATRIX ",
     &       "NBLKloc, NBLK=", NBLKloc, NBLK
      ENDIF
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      ALLOCATE(WT(NBLK), WNBINCOL(NBLK), STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1) = -7
           INFO(2) = 2*NBLK
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT%COL; WT"
           END IF
           GOTO 500
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
       IF ( INFO(1) .LT. 0 ) GOTO 500
      DO JB=1, NBLK
       WT(JB) = LMAT%COL(JB)%NBINCOL
      ENDDO
      DO JB=1,NBLK
       IF ( LMAT%COL(JB)%NBINCOL.EQ.0) CYCLE
       DO IB=1,  LMAT%COL(JB)%NBINCOL
        I = LMAT%COL(JB)%IRN(IB)
        WT(I)= WT(I)+1
       ENDDO
      ENDDO
      CALL MPI_ALLREDUCE(WT(1), WNBINCOL(1), NBLK, 
     &        MPI_INTEGER, MPI_SUM, COMM, IERR)
      IF (allocated(WT)) DEALLOCATE(WT)
      IF (MAPCOLonLUMAT) THEN
       NNZ = 0_8
       DO I=1, NBLK
        NNZ=NNZ+WNBINCOL(I)
       ENDDO
       CALL  MUMPS_AB_COMPUTE_MAPCOL (OPTION, INFO, ICNTL,
     &     MYID, NNZ, WNBINCOL(1), NBLK, 
     &     NBLK, NPROCS, MAPCOL(1))
       CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
       IF ( INFO(1) .LT. 0 ) GOTO 500
      ENDIF
      LUMAT%NBCOL = NBLK
      LUMAT%NZL   = 0_8
      ALLOCATE(LUMAT%COL(NBLK), STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1) = -7
           INFO(2) = NBLK
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT%COL; WT"
           END IF
      ENDIF
      IF ( INFO(1) .GE. 0 ) THEN
       DO JB=1,NBLK
        NB =  WNBINCOL(JB)
        IF (MAPCOL_IN_NSTEPS) THEN
         IF (MAPCOL(abs(STEP(JB))).EQ.MYID) THEN
          LUMAT%NZL             = LUMAT%NZL + int(NB,8)
         ELSE
          NB = 0
         ENDIF
        ELSE
         IF (MAPCOL(JB).EQ.MYID) THEN
          LUMAT%NZL             = LUMAT%NZL + int(NB,8)
         ELSE
          NB = 0
         ENDIF
        ENDIF
        LUMAT%COL(JB)%NBINCOL = NB
        IF (NB.GT.0) THEN
         ALLOCATE(LUMAT%COL(JB)%IRN(NB), STAT=IERR)
         IF (IERR.NE.0) THEN
           INFO(1)  = -7
           INFO(2)  = NB
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LMAT%COL"
           END IF
           EXIT
         ENDIF
        ENDIF
       ENDDO
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      IF (allocated(WNBINCOL)) DEALLOCATE(WNBINCOL)
      CALL MPI_ALLREDUCE(LUMAT%NZL, NZ_locMAX8, 1, MPI_INTEGER8,
     &                   MPI_MAX, COMM, IERR)
      NBRECORDS = KEEP(39)
      IF (NZ_locMAX8 .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(NZ_locMAX8)
      ENDIF
      CALL MUMPS_AB_DIST_LMAT_TO_LUMAT ( 
     &  .TRUE.,   
     &  MAPCOL_IN_NSTEPS,  
     &  INFO, ICNTL, COMM, MYID, NBLK, NPROCS,
     &  LMAT, MAPCOL, SIZEMAPCOL, STEP, SIZESTEP, 
     &  LUMAT, NBRECORDS, NSEND8, NLOCAL8
     &  )
      CALL MUMPS_AB_FREE_LMAT(LMAT)
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      ALLOCATE(WT(NBLK), STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1) = -7
           INFO(2) = 2*NBLK
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT%COL; WT"
           END IF
           GOTO 500
      ENDIF
      CALL MUMPS_AB_LOCALCLEAN_LMAT ( MYID,
     &     NBLK, LUMAT, WT(1), INFO(1), INFO(2), LP, LPOK
     & )
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      DEALLOCATE(WT)
      GOTO 600
  500 CONTINUE
      IF (allocated(WT)) DEALLOCATE(WT)
      IF (allocated(WNBINCOL)) DEALLOCATE(WNBINCOL)
  600 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_AB_BUILD_DCLEAN_LUMATRIX
      SUBROUTINE MUMPS_INIALIZE_REDIST_LUMAT (
     &  INFO, ICNTL, KEEP, COMM, MYID, NBLK, 
     &  LUMAT, PROCNODE_STEPS, NSTEPS, MAPCOL,
     &  LUMAT_REMAP, NBRECORDS, STEP 
     &  )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER  :: IERR, MASTER
      PARAMETER (MASTER=0)
      INTEGER, INTENT(IN) :: MYID, NBLK, NSTEPS, KEEP(500)
      INTEGER, INTENT(IN)   :: ICNTL(60), COMM
      INTEGER               :: INFO(80)
      INTEGER, INTENT(IN) :: PROCNODE_STEPS(NSTEPS)
      TYPE(LMATRIX_T), INTENT(IN)  :: LUMAT
      INTEGER, INTENT(IN)          :: STEP(NBLK)
      TYPE(LMATRIX_T), INTENT(INOUT) :: LUMAT_REMAP
      INTEGER, INTENT(OUT) :: NBRECORDS
      INTEGER, INTENT(OUT) :: MAPCOL(NSTEPS)
      INTEGER :: LP, MP, ISTEP, JB, NB
      LOGICAL :: LPOK
      INTEGER, ALLOCATABLE, DIMENSION(:) :: WT, WNBINCOL
      INTEGER MUMPS_PROCNODE
      INTEGER(8) :: NZ_locMAX8
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      ALLOCATE(WT(NBLK), WNBINCOL(NBLK), STAT=IERR)
      IF (IERR.NE.0) THEN
        INFO(1) = -7
        INFO(2) = 2*NBLK
        IF ( LPOK ) THEN
           WRITE(LP, *) " ERROR allocate WT"
        END IF
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      DO JB=1, NBLK
       WT(JB) = LUMAT%COL(JB)%NBINCOL
      ENDDO
      CALL MPI_ALLREDUCE(WT(1), WNBINCOL(1), NBLK, 
     &        MPI_INTEGER, MPI_SUM, COMM, IERR)
      IF (allocated(WT)) DEALLOCATE(WT)
      IF (MYID.EQ.MASTER) THEN
       DO ISTEP=1, NSTEPS
        MAPCOL(ISTEP) = 
     &              MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))
       ENDDO
      ENDIF
      CALL MPI_BCAST( MAPCOL(1), NSTEPS, MPI_INTEGER,
     &     MASTER, COMM, IERR )
      CALL MPI_BCAST( STEP(1), NBLK, MPI_INTEGER,
     &     MASTER, COMM, IERR )
      LUMAT_REMAP%NBCOL = NBLK
      ALLOCATE(LUMAT_REMAP%COL(NBLK), STAT=IERR)
      IF (IERR.NE.0) THEN
           INFO(1) = -7
           INFO(2) = NBLK
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT_REMAP%COL"
           END IF
      ENDIF
      IF ( INFO(1) .GE. 0 ) THEN
        LUMAT_REMAP%NZL = 0_8
       DO JB=1,NBLK
        NB =  WNBINCOL(JB)
        IF (MAPCOL(abs(STEP(JB))).EQ.MYID) THEN
         LUMAT_REMAP%NZL             = LUMAT_REMAP%NZL + int(NB,8)
        ELSE
         NB = 0
        ENDIF
        LUMAT_REMAP%COL(JB)%NBINCOL = NB
        IF (NB.GT.0) THEN
         ALLOCATE(LUMAT_REMAP%COL(JB)%IRN(NB), STAT=IERR)
         IF (IERR.NE.0) THEN
           INFO(1)  = -7
           INFO(2)  = NB
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate of LUMAT_REMAP%COL"
           END IF
           EXIT
         ENDIF
        ENDIF
       ENDDO
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      IF (allocated(WNBINCOL)) DEALLOCATE(WNBINCOL)
       CALL MPI_ALLREDUCE(LUMAT_REMAP%NZL, NZ_locMAX8, 1, MPI_INTEGER8,
     &                   MPI_MAX, COMM, IERR)
        NBRECORDS = KEEP(39)
      IF (NZ_locMAX8 .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(NZ_locMAX8)
      ENDIF
      GOTO 600
  500 CONTINUE
      IF (allocated(WT)) DEALLOCATE(WT)
      IF (allocated(WNBINCOL)) DEALLOCATE(WNBINCOL)
  600 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_INIALIZE_REDIST_LUMAT
      SUBROUTINE MUMPS_AB_DCOORD_TO_DCOMPG ( 
     &     MYID, NPROCS, COMM, 
     &     NBLK, NDOF, NNZ, 
     &     IRN, JCN, DOF2BLOCK, 
     &     ICNTL, INFO, KEEP, 
     &     LUMAT, GCOMP, READY_FOR_ANA_F)
      USE MUMPS_ANA_BLK_M, ONLY: LMATRIX_T, COMPACT_GRAPH_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER, INTENT(IN)    :: MYID, NPROCS, NBLK, NDOF
      INTEGER(8), INTENT(IN) :: NNZ
      INTEGER, INTENT(IN)    :: IRN(max(1_8,NNZ)), JCN(max(1_8,NNZ))
      LOGICAL, INTENT(IN)    :: READY_FOR_ANA_F
      INTEGER, INTENT(INOUT) :: DOF2BLOCK(NDOF)
      INTEGER, INTENT(IN)    :: ICNTL(60), COMM
      INTEGER, INTENT(INOUT) :: KEEP(500), INFO(80)
      TYPE(COMPACT_GRAPH_T)  :: GCOMP
      TYPE(LMATRIX_T)        :: LUMAT  
      TYPE(LMATRIX_T)        :: LMAT
      INTEGER :: IDUMMY_ARRAY(1)
      INTEGER :: allocok, LP, MPG
      LOGICAL :: LPOK, PROKG
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAPCOL
      LOGICAL :: MAPCOLonLUMAT, MAPCOL_IN_NSTEPS
      INTEGER OPTION
      PARAMETER (OPTION=2)
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MPG = ICNTL( 3 )
      PROKG = ( MPG .GT. 0 .and. MYID .eq. MASTER )
      MAPCOLonLUMAT    = .FALSE.
      MAPCOL_IN_NSTEPS = .FALSE.
      IF (KEEP(14).EQ.1) THEN
         CALL MUMPS_ABORT()
      ENDIF
      IF (KEEP(14).EQ.0) THEN
         CALL MPI_BCAST( DOF2BLOCK, NDOF, MPI_INTEGER, MASTER, 
     &     COMM, IERR )
      ENDIF
      CALL MUMPS_AB_COORD_TO_LMAT (  MYID, 
     &          NBLK, NDOF, NNZ, IRN, JCN,
     &          DOF2BLOCK, 
     &          INFO(1), INFO(2), LP, LPOK, 
     &          LMAT)
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      ALLOCATE(MAPCOL(NBLK), STAT=allocok) 
      IF (allocok.NE.0) THEN
           INFO(1)  = -7
           INFO(2)  = NBLK
           IF ( LPOK ) THEN
              WRITE(LP, *) " ERROR allocate MAPCOL of size", 
     &           INFO(2)
           END IF
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      IF (.NOT.MAPCOLonLUMAT) THEN 
        CALL MUMPS_AB_COL_DISTRIBUTION (OPTION, 
     &     INFO, ICNTL, COMM, NBLK, MYID, NPROCS,
     &     LMAT, MAPCOL)
        CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
      ENDIF
      CALL MUMPS_AB_BUILD_DCLEAN_LUMATRIX (
     &    MAPCOLonLUMAT, MAPCOL_IN_NSTEPS,
     &    INFO, ICNTL, KEEP, COMM,  MYID, NBLK, NPROCS,
     &    LMAT, MAPCOL, NBLK, 
     &    IDUMMY_ARRAY, 1, 
     &    LUMAT)
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
       IF (allocated(MAPCOL))  DEALLOCATE(MAPCOL)
       CALL MUMPS_AB_LMAT_TO_CLEAN_G ( MYID, .FALSE., 
     &         READY_FOR_ANA_F,
     &         LUMAT, GCOMP, INFO, ICNTL
     &      )
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      IF (KEEP(494).EQ.0) THEN
       CALL MUMPS_AB_FREE_LMAT(LUMAT)
      ENDIF
      GOTO 600
  500 CONTINUE
      IF (allocated(MAPCOL))  DEALLOCATE(MAPCOL)
      CALL MUMPS_AB_FREE_LMAT(LMAT)
      CALL MUMPS_AB_FREE_LMAT(LUMAT)
  600 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_AB_DCOORD_TO_DCOMPG
      SUBROUTINE MUMPS_AB_DCOORD_TO_DTREE_LUMAT ( 
     &     MYID, NPROCS, COMM, 
     &     NBLK, NDOF, NNZ, 
     &     IRN, JCN, 
     &     PROCNODE_STEPS, NSTEPS, STEP,
     &     ICNTL, INFO, KEEP, 
     &     MAPCOL, LUMAT)
      USE MUMPS_ANA_BLK_M, ONLY: LMATRIX_T, COMPACT_GRAPH_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER, INTENT(IN)    :: MYID, NPROCS, NBLK, NDOF, NSTEPS
      INTEGER(8), INTENT(IN) :: NNZ
      INTEGER, INTENT(IN)    :: IRN(max(1_8,NNZ)), JCN(max(1_8,NNZ))
      INTEGER, INTENT(IN)    :: ICNTL(60), COMM
      INTEGER, INTENT(IN)    :: PROCNODE_STEPS(NSTEPS)
      INTEGER, INTENT(IN)    :: STEP(NBLK)
      INTEGER, INTENT(INOUT) :: KEEP(500), INFO(80)
      INTEGER, INTENT(OUT)   :: MAPCOL(NSTEPS)
      TYPE(LMATRIX_T)        :: LUMAT  
      INTEGER, DIMENSION(:), allocatable:: DOF2BLOCK
      TYPE(LMATRIX_T)        :: LMAT
      INTEGER :: allocok, LP
      LOGICAL :: LPOK
      INTEGER :: IDOF, ISTEP
      LOGICAL :: MAPCOL_IN_NSTEPS, MAPCOLonLUMAT
      INTEGER OPTION
      PARAMETER (OPTION=2)
      INTEGER MUMPS_PROCNODE
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MAPCOLonLUMAT    = .FALSE.
      MAPCOL_IN_NSTEPS = .TRUE.
      IF (KEEP(14).EQ.1) THEN
         CALL MUMPS_ABORT()
      ENDIF
      allocate(DOF2BLOCK(NDOF), STAT=allocok)
      IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           INFO( 2 ) = NDOF
           IF ( LPOK ) WRITE(LP, 150) ' DOF2BLOCK'
      ENDIF
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      DO IDOF=1, NDOF
       DOF2BLOCK(IDOF) = IDOF
      ENDDO
      CALL MUMPS_AB_COORD_TO_LMAT (  MYID, 
     &          NBLK, NDOF, NNZ, IRN, JCN,
     &          DOF2BLOCK, 
     &          INFO(1), INFO(2), LP, LPOK, 
     &          LMAT)
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      IF (allocated(DOF2BLOCK))  DEALLOCATE(DOF2BLOCK)
      IF (MYID.EQ.MASTER) THEN
       DO ISTEP=1, NSTEPS
        MAPCOL(ISTEP) = 
     &              MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))
       ENDDO
      ENDIF
      CALL MPI_BCAST( MAPCOL(1), NSTEPS, MPI_INTEGER,
     &     MASTER, COMM, IERR )
      CALL MPI_BCAST( STEP(1), NBLK, MPI_INTEGER,
     &     MASTER, COMM, IERR )
      CALL MUMPS_AB_BUILD_DCLEAN_LUMATRIX(
     &    MAPCOLonLUMAT, MAPCOL_IN_NSTEPS,
     &    INFO, ICNTL, KEEP, COMM,  MYID, NBLK, NPROCS,
     &    LMAT, MAPCOL, NSTEPS,
     &    STEP, NBLK, LUMAT)
      CALL MUMPS_PROPINFO(  ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      GOTO 600
  500 CONTINUE
      IF (allocated(DOF2BLOCK))  DEALLOCATE(DOF2BLOCK)
      CALL MUMPS_AB_FREE_LMAT(LMAT)
      CALL MUMPS_AB_FREE_LMAT(LUMAT)
  600 CONTINUE
      RETURN
 150  FORMAT(
     & /' ** FAILURE IN MUMPS_AB_DCOORD_TO_DTREE_LUMAT, ', 
     &  ' DYNAMIC ALLOCATION OF ',
     &     A30)
      END SUBROUTINE MUMPS_AB_DCOORD_TO_DTREE_LUMAT
      SUBROUTINE MUMPS_AB_DIST_LMAT_TO_LUMAT ( 
     &  UNFOLD,
     &  MAPCOL_IN_NSTEPS,
     &  INFO, ICNTL, COMM, MYID, NBLK, SLAVEF,
     &  LMAT, MAPCOL, SIZEMAPCOL, STEP, SIZESTEP, 
     &  LUMAT, NBRECORDS, NSEND8, NLOCAL8
     &  )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER  :: IERR, MASTER, MSGSOU
      PARAMETER (MASTER=0)
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      LOGICAL, INTENT(IN) :: UNFOLD, MAPCOL_IN_NSTEPS
      INTEGER, INTENT(IN) :: MYID, SLAVEF, NBLK 
      INTEGER, INTENT(IN) :: SIZEMAPCOL, SIZESTEP 
      INTEGER, INTENT(IN) :: ICNTL(60), COMM, NBRECORDS
      INTEGER             :: INFO(80)
      TYPE(LMATRIX_T), INTENT(IN)  :: LMAT
      INTEGER, INTENT(IN)          :: MAPCOL(SIZEMAPCOL)
      INTEGER, INTENT(IN)          :: STEP(SIZESTEP)
      TYPE(LMATRIX_T), INTENT(INOUT) :: LUMAT
      INTEGER(8), INTENT(OUT) :: NSEND8, NLOCAL8
      INTEGER :: LP, MP, allocok
      INTEGER :: IB, JB, I, II, ISEND, JSEND, ITOSEND
      LOGICAL :: LPOK 
      INTEGER :: NBTOSEND
      INTEGER END_MSG_2_RECV
      INTEGER KPROBE, FREQPROBE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTLOC
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: BUFI
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFRECI
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IACT, IREQI
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: SEND_ACTIVE
      INTEGER  :: DEST
      LOGICAL  :: FLAG
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      IF (UNFOLD) THEN
        NBTOSEND = 2
      ELSE
        NBTOSEND = 1
      ENDIF
      NSEND8  = 0_8
      NLOCAL8 = 0_8
      END_MSG_2_RECV = SLAVEF-1
      ALLOCATE( IACT(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IACT in matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( IREQI(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IREQI in matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( SEND_ACTIVE(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating SEND_ACTIVE in matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( BUFI( NBRECORDS * 2 + 1, 2, SLAVEF ), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating int buffer for matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = ( NBRECORDS * 2 + 1 ) * SLAVEF * 2
        GOTO 20
      END IF
      ALLOCATE( BUFRECI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = NBRECORDS * 2 + 1
        GOTO 20
      END IF
      ALLOCATE( PTLOC( NBLK ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -7
        INFO(2) = NBLK
        GOTO 20
      END IF
 20   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 100
      DO I = 1, SLAVEF
        BUFI( 1, 1, I ) = 0
      END DO
      DO I = 1, SLAVEF
        BUFI( 1, 2, I ) = 0
      END DO
      DO I = 1, SLAVEF
        SEND_ACTIVE( I ) = .FALSE.
        IACT( I ) = 1
      END DO
      DO I = 1, NBLK
        PTLOC(I) = 0
      END DO
      KPROBE = 0
      FREQPROBE = max(1,NBRECORDS/10)
      IF (SLAVEF .EQ. 1) FREQPROBE = huge(FREQPROBE)
      DO JB=1,NBLK
       IF ( LMAT%COL(JB)%NBINCOL.EQ.0) CYCLE
       DO II=1,  LMAT%COL(JB)%NBINCOL
        KPROBE = KPROBE + 1
        IF ( KPROBE .eq. FREQPROBE ) THEN
          KPROBE = 0
          CALL MPI_IPROBE( MPI_ANY_SOURCE, LMATDIST, COMM,
     &                     FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            MSGSOU = STATUS( MPI_SOURCE )
            CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, 
     &                 MPI_INTEGER,
     &                 MSGSOU, LMATDIST, COMM, STATUS, IERR )
            CALL MUMPS_AB_LMAT_TREAT_RECV_BUF(
     &             MYID, BUFRECI(1), NBRECORDS, LUMAT, 
     &             NBLK, PTLOC(1), END_MSG_2_RECV
     &             )
          END IF
        END IF
        IB = LMAT%COL(JB)%IRN(II)
        DO ITOSEND=1,NBTOSEND
         IF (ITOSEND.EQ.1) THEN
          IF (MAPCOL_IN_NSTEPS) THEN
           DEST  = MAPCOL(abs(STEP(JB)))
          ELSE
           DEST  = MAPCOL(JB)
          ENDIF
          ISEND = IB
          JSEND = JB
         ELSE
          IF (MAPCOL_IN_NSTEPS) THEN
           DEST  = MAPCOL(abs(STEP(IB)))
          ELSE
           DEST  = MAPCOL(IB)
          ENDIF
          ISEND = JB
          JSEND = IB
         ENDIF
         IF (DEST.EQ.MYID) THEN
          LUMAT%COL(JSEND)%IRN(1+PTLOC(JSEND))= ISEND
          PTLOC(JSEND)    = PTLOC(JSEND) + 1
          NLOCAL8 = NLOCAL8 + 1_8
         ELSE
          NSEND8 = NSEND8 + 1_8
          CALL MUMPS_AB_LMAT_FILL_BUFFER( 
     &     DEST, ISEND, JSEND, NBLK,
     &     BUFI, BUFRECI, PTLOC,
     &     NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, 
     &     SEND_ACTIVE, LMAT, LUMAT, END_MSG_2_RECV
     &      )
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      DEST = -3
      CALL MUMPS_AB_LMAT_FILL_BUFFER(DEST, ISEND, JSEND, 
     &     NBLK, BUFI, BUFRECI, PTLOC,
     &     NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, 
     &     SEND_ACTIVE, LMAT, LUMAT, END_MSG_2_RECV
     &      )
      DO WHILE ( END_MSG_2_RECV .NE. 0 )
        CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, MPI_INTEGER,
     &                 MPI_ANY_SOURCE, LMATDIST, COMM, STATUS, IERR )
        CALL MUMPS_AB_LMAT_TREAT_RECV_BUF(
     &             MYID, BUFRECI(1), NBRECORDS, LUMAT, 
     &             NBLK, PTLOC(1), END_MSG_2_RECV
     &             )
      END DO
      DO I = 1, SLAVEF
        IF ( SEND_ACTIVE( I ) ) THEN
          CALL MPI_WAIT( IREQI( I ), STATUS, IERR )
        END IF
      END DO
 100  CONTINUE
      IF (ALLOCATED(PTLOC))   DEALLOCATE( PTLOC )
      IF (ALLOCATED(BUFI))    DEALLOCATE( BUFI )
      IF (ALLOCATED(BUFRECI)) DEALLOCATE( BUFRECI )
      IF (ALLOCATED(IACT))    DEALLOCATE( IACT )
      IF (ALLOCATED(IREQI))   DEALLOCATE( IREQI )
      IF (ALLOCATED(SEND_ACTIVE)) DEALLOCATE( SEND_ACTIVE )
      RETURN
      END SUBROUTINE MUMPS_AB_DIST_LMAT_TO_LUMAT
      SUBROUTINE MUMPS_AB_LMAT_TREAT_RECV_BUF (
     &   MYID, BUFI, NBRECORDS, LUMAT,  
     &   NBLK, PTLOC, END_MSG_2_RECV
     & )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER, INTENT(IN)   :: NBLK, MYID, NBRECORDS
      INTEGER, INTENT(IN)   :: BUFI( NBRECORDS * 2 + 1 ) 
      INTEGER, INTENT(INOUT):: END_MSG_2_RECV, PTLOC(NBLK)
      TYPE(LMATRIX_T), INTENT(INOUT) :: LUMAT
      INTEGER :: IREC, NB_REC, IB, JB
      NB_REC = BUFI( 1 )
      IF ( NB_REC .LE. 0 ) THEN
        END_MSG_2_RECV = END_MSG_2_RECV - 1
        NB_REC = - NB_REC
      END IF
      IF ( NB_REC .eq. 0 ) RETURN
      DO IREC = 1, NB_REC
       IB = BUFI( IREC * 2 )
       JB = BUFI( IREC * 2 + 1 )
       LUMAT%COL(JB)%IRN(1+PTLOC(JB))= IB
       PTLOC(JB) = PTLOC(JB) + 1
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_LMAT_TREAT_RECV_BUF
      SUBROUTINE MUMPS_AB_LMAT_FILL_BUFFER (
     &     DEST, ISEND, JSEND, NBLK,
     &     BUFI, BUFRECI, PTLOC,
     &     NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, 
     &     SEND_ACTIVE, LMAT, LUMAT, END_MSG_2_RECV
     &      )
      USE MUMPS_ANA_BLK_M, ONLY : LMATRIX_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER, INTENT(IN)  :: DEST, ISEND, JSEND, SLAVEF, COMM, MYID, 
     &                        NBRECORDS, NBLK
      INTEGER, INTENT(INOUT) :: END_MSG_2_RECV, PTLOC(NBLK)
      TYPE(LMATRIX_T), INTENT(IN)     :: LMAT
      TYPE(LMATRIX_T), INTENT(INOUT)  :: LUMAT
      LOGICAL, INTENT(INOUT) ::  SEND_ACTIVE(SLAVEF)
      INTEGER, INTENT(INOUT) ::  IREQI(SLAVEF), IACT(SLAVEF)
      INTEGER, INTENT(INOUT) ::  BUFI( NBRECORDS * 2 + 1, 2, SLAVEF )
      INTEGER, INTENT(INOUT) ::  BUFRECI( NBRECORDS * 2 + 1)
      INTEGER :: IBEG, IEND, ISLAVE, TAILLE_SEND_I, IREQ, MSGSOU,
     &           NBREC, IERR
      LOGICAL :: FLAG
      IF ( DEST .eq. -3 ) THEN
        IBEG = 1
        IEND = SLAVEF
      ELSE
        IBEG = DEST + 1
        IEND = DEST + 1
      END IF
      DO ISLAVE = IBEG, IEND
        NBREC = BUFI(1,IACT(ISLAVE),ISLAVE)
        IF ( DEST .eq. -3 ) THEN
          BUFI(1,IACT(ISLAVE),ISLAVE) = - NBREC
        END IF
        IF ( DEST .eq. -3 .or. NBREC + 1 > NBRECORDS ) THEN
          DO WHILE ( SEND_ACTIVE( ISLAVE ) )
            CALL MPI_TEST( IREQI( ISLAVE ), FLAG, STATUS, IERR )
            IF ( .NOT. FLAG ) THEN
                CALL MPI_IPROBE( MPI_ANY_SOURCE, LMATDIST, COMM,
     &                           FLAG, STATUS, IERR )
                IF ( FLAG ) THEN
                  MSGSOU = STATUS(MPI_SOURCE)
                  CALL MPI_RECV( BUFRECI(1), 2*NBRECORDS+1,
     &                  MPI_INTEGER, MSGSOU, LMATDIST, COMM,
     &                  STATUS, IERR )
                  CALL MUMPS_AB_LMAT_TREAT_RECV_BUF(
     &             MYID, BUFRECI, NBRECORDS, LUMAT, 
     &             NBLK, PTLOC(1), END_MSG_2_RECV
     &             )
                END IF
            ELSE
              SEND_ACTIVE( ISLAVE ) = .FALSE.
            END IF
          END DO
          IF ( ISLAVE - 1 .ne. MYID ) THEN
            TAILLE_SEND_I = NBREC * 2 + 1
            CALL MPI_ISEND( BUFI(1, IACT(ISLAVE), ISLAVE ),
     &          TAILLE_SEND_I,
     &          MPI_INTEGER, ISLAVE - 1, LMATDIST, COMM,
     &          IREQI( ISLAVE ), IERR )
            SEND_ACTIVE( ISLAVE ) = .TRUE.
          ELSE
            IF (NBREC.NE.0) THEN 
              write(*,*) " Internal error in ", 
     &                    " MUMPS_AB_LMAT_FILL_BUFFER "
              CALL MUMPS_ABORT()
            ENDIF
          END IF
          IACT( ISLAVE ) = 3 - IACT( ISLAVE )
          BUFI( 1, IACT( ISLAVE ), ISLAVE ) = 0
        END IF
        IF ( DEST .ne. -3 ) THEN
          IREQ = BUFI(1,IACT(ISLAVE),ISLAVE) + 1
          BUFI(1,IACT(ISLAVE),ISLAVE) = IREQ
          BUFI(IREQ*2,IACT(ISLAVE),ISLAVE)  = ISEND
          BUFI(IREQ*2+1,IACT(ISLAVE),ISLAVE) = JSEND
        END IF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_AB_LMAT_FILL_BUFFER
      SUBROUTINE MUMPS_AB_GATHER_GRAPH (
     &    ICNTL, KEEP, COMM, MYID, NPROCS, INFO, 
     &    GCOMP_DIST, GCOMP)
      USE MUMPS_ANA_BLK_M, ONLY : COMPACT_GRAPH_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      TYPE(COMPACT_GRAPH_T), INTENT(IN)  :: GCOMP_DIST
      INTEGER, INTENT(IN)    :: MYID, NPROCS,  ICNTL(60), COMM,
     &                          KEEP(500)
      INTEGER, INTENT(INOUT) ::  INFO(80)
      TYPE(COMPACT_GRAPH_T) :: GCOMP
       INTEGER    :: NG, allocok, LP, MPG, I, J, K
       INTEGER    :: INDX, NB_BLOCK_SENT, MAX_NBBLOCK_loc, NRECV,
     &               BLOCKSIZE, SIZE_SENT, NB_BLOCKS, NBNONEMPTY,
     &               FIRSTNONEMPTY, LASTNONEMPTY, NBBLOCK_loc
       INTEGER(4) :: IOVFLO
       INTEGER(8) :: NZG, NZG_CENT, I8, IBEG8, IEND8, 
     &               SIZEGCOMPALLOCATED
       LOGICAL    :: LPOK, PROKG
       INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IQ
      INTEGER, ALLOCATABLE :: REQPTR(:)
      INTEGER(8), ALLOCATABLE :: GPTR(:), GPTR_cp(:)
      LP  = ICNTL( 1 )
      LPOK  = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      MPG = ICNTL( 3 )
      PROKG = ( MPG .GT. 0 .and. MYID .eq. MASTER )
      PROKG = (PROKG.AND.(ICNTL(4).GE.2))
      IOVFLO = huge(IOVFLO)  
      BLOCKSIZE = int(max(100000_8,int(IOVFLO,8)/200_8))
      NZG =  GCOMP_DIST%NZG
      NG  =  GCOMP_DIST%NG
      CALL MPI_REDUCE( NZG, NZG_CENT, 1, MPI_INTEGER8, 
     &     MPI_SUM, MASTER, COMM, IERR )
      IF (MYID.EQ.MASTER) THEN
        GCOMP%NZG = NZG_CENT
        GCOMP%NG  = NG
        SIZEGCOMPALLOCATED     = NZG_CENT+int(NG,8)+1_8
        GCOMP%SIZEADJALLOCATED = SIZEGCOMPALLOCATED
        ALLOCATE( GCOMP%ADJ(SIZEGCOMPALLOCATED),
     &          GCOMP%IPE(NG+1), 
     &          GPTR( NPROCS ),
     &          GPTR_cp( NPROCS ),
     &          REQPTR( NPROCS-1 ),
     &          IQ(NG+1),STAT=allocok)
        IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           CALL MUMPS_SET_IERROR(
     &        NZG_CENT + 3_8*int(NG,8)+3_8+3_8*int(NPROCS,8)-1_8, 
     &        INFO(2))
           IF ( LPOK )
     &     WRITE(LP, *) " ERROR allocating graph in",
     &                  " MUMPS_AB_GATHER_GRAPH"
        ENDIF
      ELSE
        ALLOCATE( IQ(NG+1), STAT=allocok)
        IF (allocok.NE.0) THEN
           INFO( 1 ) = -7
           INFO( 2 ) = NG+1
           IF ( LPOK )
     &     WRITE(LP, *) " ERROR allocating pointers",
     &                  " MUMPS_AB_GATHER_GRAPH"
        END IF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     COMM, MYID )
      IF (INFO(1).LT.0) GOTO 500
      FIRSTNONEMPTY = 0
      LASTNONEMPTY  = -1
      DO I=1,NG
       IQ(I) = int(GCOMP_DIST%IPE(I+1)-GCOMP_DIST%IPE(I))
       IF (IQ(I).NE.0) THEN
         IF (FIRSTNONEMPTY.EQ.0) FIRSTNONEMPTY=I
         LASTNONEMPTY = I
       ENDIF
      ENDDO
      NBNONEMPTY = LASTNONEMPTY-FIRSTNONEMPTY+1
      IF (MYID.EQ.MASTER) THEN
       DO J=1, NG
        GCOMP%IPE(J) = 0
       ENDDO
       J=FIRSTNONEMPTY   
       IF (NBNONEMPTY.GT.0) THEN
        DO I=FIRSTNONEMPTY, LASTNONEMPTY
          GCOMP%IPE(J) = IQ(I)
          J = J+1
        ENDDO
       ENDIF
       DO I = 1, NPROCS - 1
         CALL MPI_RECV( NBNONEMPTY, 1, 
     &           MPI_INTEGER, I,
     &           GATHERG_NB, COMM, STATUS, IERR )
         IF (NBNONEMPTY.GT.0) THEN
           CALL MPI_RECV( J, 1, 
     &           MPI_INTEGER, I,
     &           GATHERG_FIRST, COMM, STATUS, IERR )
           CALL MPI_RECV( GCOMP%IPE(J), NBNONEMPTY, 
     &           MPI_INTEGER8, I,
     &           GATHERG_IPE, COMM, STATUS, IERR )
         ENDIF
       ENDDO
      ELSE
        CALL MPI_SEND( NBNONEMPTY, 1, MPI_INTEGER, MASTER,
     &       GATHERG_NB, COMM, IERR )
        IF (NBNONEMPTY.GT.0) THEN
          CALL MPI_SEND( FIRSTNONEMPTY, 1, MPI_INTEGER, MASTER,
     &       GATHERG_FIRST, COMM, IERR )
          CALL MPI_SEND( IQ(FIRSTNONEMPTY), NBNONEMPTY, 
     &       MPI_INTEGER8, MASTER,
     &       GATHERG_IPE, COMM, IERR )
        ENDIF
      ENDIF
      IF (MYID.EQ.MASTER) THEN
       IQ(1) = 1_8
       DO I=1,NG
         IQ(I+1) = IQ(I) + GCOMP%IPE(I)
         GCOMP%IPE(I) = IQ(I)
       ENDDO
       GCOMP%IPE(NG+1) = IQ(NG+1)
       DEALLOCATE(IQ)
      ELSE
       DEALLOCATE(IQ)
      ENDIF
      IF (MYID.EQ.MASTER) THEN
        NB_BLOCK_SENT = 0
        MAX_NBBLOCK_loc  = 0
        DO I = 1, NPROCS - 1
            CALL MPI_RECV( GPTR( I+1 ), 1, 
     &           MPI_INTEGER8, I,
     &           GATHERG_NZG, COMM, STATUS, IERR )
         NBBLOCK_loc = ceiling(dble(GPTR(I+1))/dble(BLOCKSIZE))
         MAX_NBBLOCK_loc = max(MAX_NBBLOCK_loc, NBBLOCK_loc)
         NB_BLOCK_SENT = NB_BLOCK_SENT + NBBLOCK_loc
        ENDDO
        GPTR( 1 ) = NZG + 1_8
         DO I = 2, NPROCS
            GPTR( I ) = GPTR( I ) + GPTR( I-1 )
         END DO
      ELSE
        CALL MPI_SEND( NZG, 1, MPI_INTEGER8, MASTER,
     &        GATHERG_NZG, COMM, IERR )
      ENDIF
      IF (MYID.EQ.MASTER) THEN
        DO I=1, NPROCS
         GPTR_cp(I) = GPTR(I) 
        ENDDO
        IF (NZG.GT.0_8) THEN
         DO I8=1, NZG
          GCOMP%ADJ(I8) = GCOMP_DIST%ADJ(I8)
         ENDDO
        ENDIF
        NB_BLOCKS = 0
        DO K = 1, MAX_NBBLOCK_loc
         NRECV = 0
         DO I = 1, NPROCS - 1
            IBEG8     = GPTR_cp( I )
            IF ( IBEG8 .LT. GPTR(I+1))  THEN
              NRECV = NRECV + 1
              IEND8 = min(IBEG8+int(BLOCKSIZE,8)-1_8, 
     &                    GPTR(I+1)-1_8)
              GPTR_cp( I ) = IEND8 + 1_8
              SIZE_SENT   = int(IEND8 -  IBEG8 + 1_8)
              NB_BLOCKS   = NB_BLOCKS + 1
              CALL MPI_IRECV( GCOMP%ADJ(IBEG8), SIZE_SENT, 
     &           MPI_INTEGER,
     &           I, GATHERG_ADJ, COMM, REQPTR(I), IERR )
            ELSE
             REQPTR( I ) = MPI_REQUEST_NULL
            ENDIF
         END DO
         DO I = 1, NRECV
             CALL MPI_WAITANY
     &           ( NPROCS-1, REQPTR, INDX, 
     &           STATUS, IERR )
         ENDDO
      END DO
        DEALLOCATE( REQPTR )
        DEALLOCATE( GPTR )
        DEALLOCATE( GPTR_cp )
      ELSE
        IF (NZG.EQ.0) GOTO 600
        DO I8=1_8, NZG, int(BLOCKSIZE,8)
         SIZE_SENT = BLOCKSIZE
         IF (NZG-I8+1_8.LT.int(BLOCKSIZE,8)) THEN
              SIZE_SENT = int(NZG-I8+1_8)
         ENDIF
         CALL MPI_SEND( 
     &           GCOMP_DIST%ADJ(I8), SIZE_SENT,
     &           MPI_INTEGER, MASTER,
     &           GATHERG_ADJ, COMM, IERR )
        ENDDO
      ENDIF
       GOTO 600
 500  CONTINUE
      IF (MYID.EQ.MASTER) THEN
       IF (associated(GCOMP%ADJ)) THEN
          DEALLOCATE(GCOMP%ADJ)
          nullify(GCOMP%ADJ)
       ENDIF
       IF (associated(GCOMP%IPE)) THEN
        DEALLOCATE(GCOMP%IPE)
        nullify(GCOMP%IPE)
       ENDIF
      ENDIF
 600  CONTINUE
      IF (allocated(IQ)) DEALLOCATE(IQ)
      RETURN
      END SUBROUTINE MUMPS_AB_GATHER_GRAPH
