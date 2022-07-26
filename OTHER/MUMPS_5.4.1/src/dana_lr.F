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
      MODULE DMUMPS_ANA_LR
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_STATS
      USE MUMPS_LR_COMMON
      USE MUMPS_ANA_ORD_WRAPPERS
      USE MUMPS_ANA_BLK_M, ONLY: LMATRIX_T
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE GET_CUT(IWR, NASS, NCB, LRGROUPS, NPARTSCB,
     &   NPARTSASS, CUT)
      INTEGER, INTENT(IN) :: NASS, NCB
      INTEGER, INTENT(IN) :: IWR(*)
      INTEGER, INTENT(IN), DIMENSION(:) :: LRGROUPS
      INTEGER, INTENT(OUT) :: NPARTSCB, NPARTSASS
      INTEGER, POINTER, DIMENSION(:) :: CUT
      INTEGER :: I, CURRENT_PART, CUTBUILDER,allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BIG_CUT
      ALLOCATE(BIG_CUT(max(NASS,1)+NCB+1),stat=allocok)
      IF(allocok.GT.0) THEN
         write(*,*) "Allocation error of BIG_CUT in GET_CUT"
         CALL MUMPS_ABORT()
      ENDIF
      CURRENT_PART = LRGROUPS(IWR(1))
      BIG_CUT(1) = 1
      BIG_CUT(2) = 2
      CUTBUILDER = 2
      NPARTSASS = 0
      NPARTSCB  = 0 
      DO I = 2,NASS + NCB
        IF (LRGROUPS(IWR(I)) == CURRENT_PART) THEN
          BIG_CUT(CUTBUILDER) = BIG_CUT(CUTBUILDER) + 1
        ELSE
          CUTBUILDER = CUTBUILDER + 1
          BIG_CUT(CUTBUILDER) = BIG_CUT(CUTBUILDER-1) + 1
          CURRENT_PART = LRGROUPS(IWR(I))
        END IF
        IF (I == NASS) NPARTSASS = CUTBUILDER - 1
      END DO
      IF (NASS.EQ.1) NPARTSASS= 1
      NPARTSCB = CUTBUILDER - 1 - NPARTSASS
      ALLOCATE(CUT(max(NPARTSASS,1)+NPARTSCB+1),stat=allocok)
      IF(allocok.GT.0) THEN
         write(*,*) "Allocation error of CUT in GET_CUT"
         CALL MUMPS_ABORT()
      ENDIF
      IF (NPARTSASS.EQ.0) THEN
        CUT(1) = 1
        CUT(2:2+NPARTSCB) = BIG_CUT(1:1+NPARTSCB)
      ELSE
        CUT = BIG_CUT(1:NPARTSASS+NPARTSCB+1)
      ENDIF 
      if(allocated(BIG_CUT)) DEALLOCATE(BIG_CUT)
      END SUBROUTINE GET_CUT
      SUBROUTINE SEP_GROUPING(NV, VLIST, N, NZ, LRGROUPS, NBGROUPS, IW, 
     &   LW, IPE, LEN, GROUP_SIZE, HALO_DEPTH, TRACE, WORKH, NODE,
     &   GEN2HALO, K482, K472, K469, SEP_SIZE, 
     &   KEEP10, LP, LPOK, IFLAG, IERROR)
        INTEGER(8), INTENT(IN) :: NZ, LW
        INTEGER, INTENT(IN)    :: NV, N, GROUP_SIZE, HALO_DEPTH
        INTEGER, INTENT(IN)    :: IW(LW), LEN(N), NODE, K482
        INTEGER(8), INTENT(IN) :: IPE(N+1)
        INTEGER, INTENT(IN)    :: K472, K469, SEP_SIZE, KEEP10, LP
        LOGICAL                :: LPOK
        INTEGER, INTENT(INOUT) :: NBGROUPS, WORKH(N)
        INTEGER, INTENT(INOUT) :: VLIST(NV), TRACE(N)
        INTEGER                :: LRGROUPS(:) 
        INTEGER, INTENT(INOUT) :: GEN2HALO(N)
        INTEGER, INTENT(INOUT) :: IFLAG, IERROR
        INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IPTRHALO
        INTEGER,    ALLOCATABLE, DIMENSION(:) :: PARTS, JCNHALO
        INTEGER(8) :: HALOEDGENBR
        INTEGER    :: NHALO,
     &       NBGROUPS_KWAY, I, GROUP_SIZE2, LRGROUPS_SIGN, IERR
#if defined (metis) || defined (parmetis) || defined (metis4) || defined (parmetis3)
        INTEGER :: METIS_IDX_SIZE
#endif
#if defined (scotch) || defined (ptscotch)
        INTEGER :: SCOTCH_IDX_SIZE
#endif
        CALL COMPUTE_BLR_VCS(K472, GROUP_SIZE2, GROUP_SIZE, NV)
        NBGROUPS_KWAY = MAX(NINT(dble(NV)/dble(GROUP_SIZE2)),1)
        IF (NV .GE. SEP_SIZE) THEN
          LRGROUPS_SIGN = 1
        ELSE
          LRGROUPS_SIGN = -1
        ENDIF
        IF (NBGROUPS_KWAY > 1) THEN 
          IF (K469.EQ.3) THEN
!$OMP CRITICAL(gethalo_cri)
            CALL GETHALONODES(N, IW, LW, IPE, VLIST, NV, HALO_DEPTH,
     &                 NHALO, TRACE, WORKH, NODE, LEN, HALOEDGENBR,
     &                 GEN2HALO)
            ALLOCATE(PARTS(NHALO), IPTRHALO(NHALO+1),
     &      JCNHALO(HALOEDGENBR), STAT=IERR)
            IF (IERR.GT.0) THEN
              IF (LPOK) WRITE(LP,*) 
     &        " Error allocate integer array of size: ", 
     &         int(NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR
              IFLAG  = -7
              CALL MUMPS_SET_IERROR 
     &         (int(NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR,
     &          IERROR)
            ENDIF
            CALL GETHALOGRAPH(WORKH, NHALO, N, IW, LW, IPE, IPTRHALO,
     &     JCNHALO, HALOEDGENBR,TRACE,NODE, GEN2HALO)
!$OMP END CRITICAL(gethalo_cri)
            IF (IFLAG.LT.0) RETURN
          ELSE
            CALL GETHALONODES(N, IW, LW, IPE, VLIST, NV, HALO_DEPTH,
     &                 NHALO, TRACE, WORKH, NODE, LEN, HALOEDGENBR,
     &                 GEN2HALO)
            ALLOCATE(PARTS(NHALO), IPTRHALO(NHALO+1),
     &      JCNHALO(HALOEDGENBR), STAT=IERR)
            IF (IERR.GT.0) THEN
              IF (LPOK) WRITE(LP,*) 
     &        " Error allocate integer array of size: ", 
     &         int(NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR
              IFLAG  = -7
              CALL MUMPS_SET_IERROR 
     &         (int(NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR,
     &          IERROR)
              RETURN
            ENDIF
            CALL GETHALOGRAPH(WORKH, NHALO, N, IW, LW, IPE, IPTRHALO, 
     &     JCNHALO, HALOEDGENBR,TRACE,NODE, GEN2HALO)
          ENDIF
          IF (K482.EQ.1) THEN
#if defined (metis) || defined (parmetis) || defined (metis4) || defined (parmetis3)
             CALL MUMPS_METIS_IDXSIZE(METIS_IDX_SIZE)
             IF (METIS_IDX_SIZE .EQ. 64) THEN
                CALL MUMPS_METIS_KWAY_MIXEDto64(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
             ELSE
               IF (KEEP10.EQ.1) THEN
                IFLAG  = -52
                IERROR = 1
               ELSE
                CALL MUMPS_METIS_KWAY_MIXEDto32(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
               ENDIF
             ENDIF
#endif
          ELSE IF (K482.EQ.2) THEN
#if defined (scotch) || defined (ptscotch)
             CALL MUMPS_SCOTCH_INTSIZE(SCOTCH_IDX_SIZE)
             IF (SCOTCH_IDX_SIZE .EQ. 32) THEN
               IF (KEEP10.EQ.1) THEN
                IFLAG  = -52
                IERROR = 2
               ELSE
                CALL MUMPS_SCOTCH_KWAY_MIXEDto32(
     &               NHALO, HALOEDGENBR, IPTRHALO, JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
               ENDIF
             ELSE
                CALL MUMPS_SCOTCH_KWAY_MIXEDto64(
     &               NHALO, HALOEDGENBR, IPTRHALO, JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
             END IF
#endif
          ELSE
            WRITE(6,*) " Internal ERROR K482=", K482
            CALL MUMPS_ABORT()
          END IF
          IF (IFLAG.LT.0) GOTO 500
          CALL GET_GLOBAL_GROUPS(PARTS,VLIST,NV,
     &     NBGROUPS_KWAY, LRGROUPS, N, NBGROUPS, LRGROUPS_SIGN)
        ELSE 
!$OMP CRITICAL(lrgrouping_cri)
          DO I=1,NV
            LRGROUPS(VLIST(I)) = LRGROUPS_SIGN*(NBGROUPS + 1)
          END DO
            NBGROUPS = NBGROUPS + 1
!$OMP END CRITICAL(lrgrouping_cri)
        END IF
  500   IF (allocated(IPTRHALO)) then
           DEALLOCATE(IPTRHALO)
        ENDIF
        IF (allocated(PARTS)) then
           DEALLOCATE(PARTS) 
        ENDIF
        IF (allocated(JCNHALO)) then
           DEALLOCATE(JCNHALO   )
        ENDIF
        RETURN
      END SUBROUTINE SEP_GROUPING
      SUBROUTINE SEP_GROUPING_AB (NV, NVEXPANDED, 
     &   VLIST, N, LRGROUPS, NBGROUPS, LUMAT, SIZEOFBLOCKS,
     &   GROUP_SIZE, HALO_DEPTH, TRACE, WORKH, NODE,
     &   GEN2HALO, K482, K472, K469, SEP_SIZE, 
     &   KEEP10, LP, LPOK, IFLAG, IERROR)
        TYPE(LMATRIX_T)        :: LUMAT
        INTEGER, INTENT(IN)    :: NV, NVEXPANDED, 
     &                            N, GROUP_SIZE, HALO_DEPTH
        INTEGER, INTENT(IN)    :: SIZEOFBLOCKS(N)
        INTEGER, INTENT(IN)    :: NODE, K482
        INTEGER, INTENT(IN)    :: K472, K469, SEP_SIZE, KEEP10, LP
        LOGICAL                :: LPOK
        INTEGER, INTENT(INOUT) :: NBGROUPS, WORKH(N)
        INTEGER, INTENT(INOUT) :: VLIST(NV), TRACE(N)
        INTEGER                :: LRGROUPS(:) 
        INTEGER, INTENT(INOUT) :: GEN2HALO(N)
        INTEGER, INTENT(INOUT) :: IFLAG, IERROR
        INTEGER(8), ALLOCATABLE, DIMENSION(:) :: IPTRHALO
        INTEGER,    ALLOCATABLE, DIMENSION(:) :: PARTS, JCNHALO
        INTEGER,    ALLOCATABLE, DIMENSION(:) :: VWGT
        INTEGER(8) :: HALOEDGENBR
        INTEGER    :: NHALO,
     &       NBGROUPS_KWAY, I, GROUP_SIZE2, LRGROUPS_SIGN, IERR
        DOUBLE PRECISION :: COMPRESS_RATIO
#if defined (metis) || defined (parmetis) || defined (metis4) || defined (parmetis3)
        INTEGER :: METIS_IDX_SIZE
#endif
#if defined (scotch) || defined (ptscotch)
        INTEGER :: SCOTCH_IDX_SIZE
#endif
        CALL COMPUTE_BLR_VCS(K472, GROUP_SIZE2, GROUP_SIZE, NVEXPANDED)
        COMPRESS_RATIO= dble(NVEXPANDED)/dble(NV)
        NBGROUPS_KWAY = MAX(NINT(dble(NVEXPANDED)/dble(GROUP_SIZE2)),1)
        NBGROUPS_KWAY = min(NBGROUPS_KWAY, NV)
        IF (NVEXPANDED .GE. SEP_SIZE) THEN
          LRGROUPS_SIGN = 1
        ELSE
          LRGROUPS_SIGN = -1
        ENDIF
        IF (NBGROUPS_KWAY > 1) THEN 
          IF (K469.EQ.3) THEN
!$OMP CRITICAL(gethalo_cri)
            CALL GETHALONODES_AB(N, LUMAT, VLIST, NV, HALO_DEPTH,
     &                 NHALO, TRACE, WORKH, NODE, HALOEDGENBR,
     &                 GEN2HALO)
            ALLOCATE(PARTS(NHALO), IPTRHALO(NHALO+1),
     &      JCNHALO(HALOEDGENBR), VWGT(NHALO), STAT=IERR)
            IF (IERR.GT.0) THEN
              IF (LPOK) WRITE(LP,*) 
     &        " Error allocate integer array of size: ", 
     &         int(2*NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR
              IFLAG  = -7
              CALL MUMPS_SET_IERROR 
     &         (int(2*NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR,
     &          IERROR)
            ENDIF
            DO I=1, NHALO
             VWGT(I) = SIZEOFBLOCKS(WORKH(I))
            ENDDO
            CALL GETHALOGRAPH_AB(WORKH, NV,
     &      NHALO, N, LUMAT, IPTRHALO,
     &      JCNHALO, HALOEDGENBR,TRACE,NODE, GEN2HALO, PARTS)
!$OMP END CRITICAL(gethalo_cri)
            IF (IFLAG.LT.0) RETURN
          ELSE
            CALL GETHALONODES_AB(N, LUMAT, VLIST, NV, HALO_DEPTH,
     &                 NHALO, TRACE, WORKH, NODE, HALOEDGENBR,
     &                 GEN2HALO)
            ALLOCATE(PARTS(NHALO), IPTRHALO(NHALO+1),
     &      JCNHALO(HALOEDGENBR), VWGT(NHALO), STAT=IERR)
            IF (IERR.GT.0) THEN
              IF (LPOK) WRITE(LP,*) 
     &        " Error allocate integer array of size: ", 
     &         int(2*NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR
              IFLAG  = -7
              CALL MUMPS_SET_IERROR 
     &         (int(2*NHALO+(KEEP10*(NHALO+1)),8) + HALOEDGENBR,
     &          IERROR)
              RETURN
            ENDIF
            DO I=1, NHALO
             VWGT(I) = SIZEOFBLOCKS(WORKH(I))
            ENDDO
            CALL GETHALOGRAPH_AB(WORKH, NV, 
     &      NHALO, N, LUMAT, IPTRHALO, 
     &     JCNHALO, HALOEDGENBR,TRACE,NODE, GEN2HALO, PARTS)
          ENDIF
          IF (K482.EQ.1) THEN
#if defined (metis) || defined (parmetis) || defined (metis4) || defined (parmetis3)
             CALL MUMPS_METIS_IDXSIZE(METIS_IDX_SIZE)
             IF (METIS_IDX_SIZE .EQ. 64) THEN
                CALL MUMPS_METIS_KWAY_AB_MIXEDto64(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO,
     &               NBGROUPS_KWAY, PARTS, VWGT, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
             ELSE
               IF (KEEP10.EQ.1) THEN
                IFLAG  = -52
                IERROR = 1
               ELSE
                CALL MUMPS_METIS_KWAY_AB_MIXEDto32(NHALO, HALOEDGENBR,
     &               IPTRHALO,
     &               JCNHALO,
     &               NBGROUPS_KWAY, PARTS, VWGT, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
               ENDIF
             ENDIF
#endif
          ELSE IF (K482.EQ.2) THEN
#if defined (scotch) || defined (ptscotch)
             CALL MUMPS_SCOTCH_INTSIZE(SCOTCH_IDX_SIZE)
             IF (SCOTCH_IDX_SIZE .EQ. 32) THEN
               IF (KEEP10.EQ.1) THEN
                IFLAG  = -52
                IERROR = 2
               ELSE
                CALL MUMPS_SCOTCH_KWAY_MIXEDto32(
     &               NHALO, HALOEDGENBR, IPTRHALO, JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
               ENDIF
             ELSE
                CALL MUMPS_SCOTCH_KWAY_MIXEDto64(
     &               NHALO, HALOEDGENBR, IPTRHALO, JCNHALO,
     &               NBGROUPS_KWAY, PARTS, LP, LPOK, KEEP10,
     &               IFLAG, IERROR)
             END IF
#endif
          ELSE
            WRITE(6,*) " Internal ERROR K482=", K482
            CALL MUMPS_ABORT()
          END IF
          IF (IFLAG.LT.0) GOTO 500
          CALL GET_GLOBAL_GROUPS(PARTS,VLIST,NV,
     &     NBGROUPS_KWAY, LRGROUPS, N, NBGROUPS, LRGROUPS_SIGN)
        ELSE 
!$OMP CRITICAL(lrgrouping_cri)
          DO I=1,NV
            LRGROUPS(VLIST(I)) = LRGROUPS_SIGN*(NBGROUPS + 1)
          END DO
            NBGROUPS = NBGROUPS + 1
!$OMP END CRITICAL(lrgrouping_cri)
        END IF
  500   IF (allocated(IPTRHALO)) then
           DEALLOCATE(IPTRHALO)
        ENDIF
        IF (allocated(PARTS)) then
           DEALLOCATE(PARTS) 
        ENDIF
        IF (allocated(JCNHALO)) then
           DEALLOCATE(JCNHALO   )
        ENDIF
        IF (allocated(VWGT)) then
           DEALLOCATE(VWGT) 
        ENDIF
        RETURN
      END SUBROUTINE SEP_GROUPING_AB
      SUBROUTINE GETHALONODES_AB(N, LUMAT, IND, NIND, PMAX, 
     &                        NHALO, TRACE, WORKH, NODE, HALOEDGENBR,
     &                        GEN2HALO)
        TYPE(LMATRIX_T)        :: LUMAT
        INTEGER,DIMENSION(:),INTENT(IN) :: IND
        INTEGER, INTENT(IN)             :: N, NODE
        INTEGER, INTENT(IN)             :: PMAX,NIND
        INTEGER, INTENT(OUT)            :: NHALO
        INTEGER, INTENT(INOUT)          :: TRACE(N), WORKH(N)
        INTEGER                         :: GEN2HALO(N)
        INTEGER(8), INTENT(OUT)         :: HALOEDGENBR
        INTEGER                         :: I, J, II
        INTEGER                         :: HALOI, NB, NEWNHALO
        INTEGER(8)                      :: SEPEDGES_TOTAL,
     &                                     SEPEDGES_INTERNAL
        WORKH(1:NIND)  = IND
        NHALO          = NIND
        NEWNHALO       = 0
        HALOEDGENBR    = 0_8
        SEPEDGES_TOTAL = 0_8
        SEPEDGES_INTERNAL = 0_8
        DO I=1,NIND
          HALOI = WORKH(I)
          GEN2HALO(HALOI) = I
          IF (TRACE(HALOI) .NE. NODE) THEN
             TRACE(HALOI) = NODE
          END IF
        ENDDO
        DO I=1,NIND
          HALOI = WORKH(I)
          NB = LUMAT%COL(HALOI)%NBINCOL
          SEPEDGES_TOTAL = SEPEDGES_TOTAL + int(NB,8)
          DO J=1, NB
           II = LUMAT%COL(HALOI)%IRN(J)
           IF (TRACE(II).NE.NODE) THEN
            NEWNHALO               = NEWNHALO + 1
            WORKH(NHALO+NEWNHALO)  = II
            GEN2HALO(II)           = NHALO+NEWNHALO
            TRACE(II)              = NODE
           ELSE
            IF (GEN2HALO(II).LE.NHALO) THEN
             SEPEDGES_INTERNAL = SEPEDGES_INTERNAL + 1_8
            ENDIF
           ENDIF
          ENDDO
        END DO
        HALOEDGENBR = SEPEDGES_TOTAL + 
     &                (SEPEDGES_TOTAL - SEPEDGES_INTERNAL)
        NHALO  = NHALO + NEWNHALO
      END SUBROUTINE GETHALONODES_AB
      SUBROUTINE GETHALOGRAPH_AB(HALO,NSEP,NHALO,
     &     N,LUMAT,IPTRHALO,JCNHALO,
     &     HALOEDGENBR,TRACE,NODE, GEN2HALO, IQ)
      INTEGER, INTENT(IN) :: N
      TYPE(LMATRIX_T)     :: LUMAT
      INTEGER,INTENT(IN):: NSEP, NHALO, NODE
      INTEGER,INTENT(IN):: GEN2HALO(N)
      INTEGER,DIMENSION(NHALO),INTENT(IN) :: HALO
      INTEGER, INTENT(IN)     :: TRACE(N)
      INTEGER(8),INTENT(IN)   :: HALOEDGENBR
      INTEGER(8), INTENT(OUT) :: IPTRHALO(NHALO+1)
      INTEGER, INTENT(OUT)    :: JCNHALO(HALOEDGENBR)
      INTEGER                 :: IQ(NHALO)
      INTEGER::I,J,NB,II,JJ,HALOI,HALOJ
      DO I=NSEP+1, NHALO
       IQ(I) = 0
      ENDDO
      DO I=1,NSEP
        HALOI = HALO(I)
        NB    = LUMAT%COL(HALOI)%NBINCOL
        IQ(I) = NB
        DO JJ=1, NB
         II = LUMAT%COL(HALOI)%IRN(JJ)
         J  = GEN2HALO(II)
         IF (J.GT.NSEP) THEN
           IQ(J) = IQ(J) + 1
         ENDIF
        ENDDO
      ENDDO
      IPTRHALO(1) = 1_8
      DO I=1,NHALO
       IPTRHALO(I+1) = IPTRHALO(I)+int(IQ(I),8)
      ENDDO
      DO I=1,NSEP
        HALOI = HALO(I)
        NB    = LUMAT%COL(HALOI)%NBINCOL
        DO JJ=1, NB
         HALOJ = LUMAT%COL(HALOI)%IRN(JJ)
         J     = GEN2HALO(HALOJ)
         JCNHALO(IPTRHALO(I)) = J
         IPTRHALO(I)          = IPTRHALO(I) + 1
         IF (J.GT.NSEP) THEN
           JCNHALO(IPTRHALO(J)) = I
           IPTRHALO(J)          = IPTRHALO(J) + 1
         ENDIF
        ENDDO
      ENDDO
      IPTRHALO(1) = 1_8
      DO I=1,NHALO
       IPTRHALO(I+1) = IPTRHALO(I)+int(IQ(I),8)
      ENDDO
      END SUBROUTINE GETHALOGRAPH_AB
      SUBROUTINE GET_GLOBAL_GROUPS(PARTS, SEP, NSEP, NPARTS, 
     &                    LRGROUPS, N, NBGROUPS, LRGROUPS_SIGN)
      INTEGER,INTENT(IN) :: NSEP, N, LRGROUPS_SIGN
      INTEGER :: PARTS(:)
      INTEGER,DIMENSION(:),INTENT(INOUT) :: SEP 
      INTEGER, INTENT(INOUT) :: NPARTS 
      INTEGER, INTENT(INOUT) :: NBGROUPS
      INTEGER                :: LRGROUPS(:)  
      INTEGER::I,CNT,NB_PARTS_WITHOUT_SEP_NODE,allocok
      INTEGER,DIMENSION(:),ALLOCATABLE::SIZES, RIGHTPART
      INTEGER,DIMENSION(:),ALLOCATABLE::PARTPTR
      INTEGER,DIMENSION(:),ALLOCATABLE :: NEWSEP
      ALLOCATE( NEWSEP(NSEP),
     &          SIZES(NPARTS),
     &          RIGHTPART(NPARTS),
     &          PARTPTR(NPARTS+1),stat=allocok)
      IF(allocok.GT.0) THEN
         write(*,*) "Allocation error in GET_GLOBAL_GROUPS"
         CALL MUMPS_ABORT()
      ENDIF
      NB_PARTS_WITHOUT_SEP_NODE = 0
      RIGHTPART = 0
      SIZES = 0
      DO I=1,NSEP
        SIZES(PARTS(I)) = SIZES(PARTS(I)) + 1
      END DO
      CNT = 0
      PARTPTR(1)=1
      DO I=2,NPARTS+1
        PARTPTR(I) = PARTPTR(I-1) + SIZES(I-1)
        IF (SIZES(I-1)==0) THEN
          NB_PARTS_WITHOUT_SEP_NODE = NB_PARTS_WITHOUT_SEP_NODE + 1
        ELSE
          CNT = CNT + 1
          RIGHTPART(I-1) = CNT
        END IF
      END DO
      NPARTS = NPARTS - NB_PARTS_WITHOUT_SEP_NODE 
!$OMP CRITICAL(lrgrouping_cri)
      DO I=1,NSEP
        NEWSEP(PARTPTR(PARTS(I))) = SEP(I)  
        LRGROUPS(SEP(I)) = LRGROUPS_SIGN*(RIGHTPART(PARTS(I)) 
     &    + NBGROUPS)
        PARTPTR(PARTS(I)) = 
     &    PARTPTR(PARTS(I)) + 1
      END DO
      NBGROUPS = NBGROUPS + NPARTS
!$OMP END CRITICAL(lrgrouping_cri)
      SEP = NEWSEP
      DEALLOCATE(NEWSEP,SIZES,RIGHTPART,PARTPTR)
      END SUBROUTINE GET_GLOBAL_GROUPS
      SUBROUTINE GETHALONODES(N, IW, LW, IPE, IND, NIND, PMAX, 
     &                        NHALO, TRACE, WORKH, NODE, LEN, CNT,
     &                        GEN2HALO)
        INTEGER,DIMENSION(:),INTENT(IN) :: IND
        INTEGER(8), INTENT(IN)          :: LW
        INTEGER, INTENT(IN)             :: N, NODE
        INTEGER, INTENT(IN)             :: IW(LW), LEN(N)
        INTEGER(8), INTENT(IN)          :: IPE(N+1)
        INTEGER, INTENT(IN)             :: PMAX,NIND
        INTEGER, INTENT(OUT)            :: NHALO
        INTEGER, INTENT(INOUT)          :: TRACE(N), WORKH(N)
        INTEGER                         :: GEN2HALO(N)
        INTEGER(8), INTENT(OUT)         :: CNT
        INTEGER                         :: DEPTH, I, LAST_LVL_START
        INTEGER                         :: HALOI
        INTEGER(8)                      :: J
        WORKH(1:NIND) = IND
        LAST_LVL_START = 1
        NHALO = NIND
        CNT = 0
        DO I=1,NIND
          HALOI = WORKH(I)
          GEN2HALO(HALOI) = I
          IF (TRACE(HALOI) .NE. NODE) THEN
             TRACE(HALOI) = NODE
          END IF
          DO J=IPE(HALOI),IPE(HALOI+1)-1
            IF (TRACE(IW(J)).EQ.NODE) THEN
              CNT = CNT + 2
            END IF
          END DO
        END DO
        DO DEPTH=1,PMAX
          CALL NEIGHBORHOOD(WORKH, NHALO, N, IW, LW, IPE,
     &                      TRACE, NODE, LEN, CNT, LAST_LVL_START,
     &                      DEPTH, PMAX, GEN2HALO)
        END DO
      END SUBROUTINE GETHALONODES
      SUBROUTINE NEIGHBORHOOD(HALO, NHALO, N, IW, LW, IPE,
     &                        TRACE, NODE, LEN, CNT, LAST_LVL_START,
     &                        DEPTH, PMAX, GEN2HALO)
        INTEGER, INTENT(IN)                 :: N, NODE, DEPTH, PMAX
        INTEGER,INTENT(INOUT)               :: NHALO, GEN2HALO(N)
        INTEGER, INTENT(INOUT)              :: LAST_LVL_START
        INTEGER(8), INTENT(INOUT)           :: CNT
        INTEGER,DIMENSION(:),INTENT(INOUT)  :: HALO
        INTEGER(8), INTENT(IN)              :: LW
        INTEGER(8), INTENT(IN)              :: IPE(N+1)
        INTEGER, TARGET, INTENT(IN)         :: IW(LW)
        INTEGER, INTENT(IN)                 :: LEN(N)
        INTEGER,DIMENSION(:)                :: TRACE
        INTEGER :: AvgDens, THRESH
        INTEGER :: I,INEI,NADJI,NEWNHALO, NEIGH
        INTEGER, DIMENSION(:), POINTER :: ADJI
        INTEGER(8) :: J
        NEWNHALO = 0
        AvgDens = nint(dble(IPE(N+1)-1_8)/dble(N))
        THRESH  = AvgDens*10
        DO I=LAST_LVL_START,NHALO
          NADJI = LEN(HALO(I))
          IF (NADJI.GT.THRESH) CYCLE
          ADJI => IW(IPE(HALO(I)):IPE(HALO(I)+1)-1)
          DO INEI=1,NADJI
            IF (TRACE(ADJI(INEI)) .NE. NODE) THEN
              NEIGH = ADJI(INEI)
              IF (LEN(NEIGH).GT.THRESH) CYCLE
              TRACE(NEIGH) = NODE
              NEWNHALO = NEWNHALO + 1
              HALO(NHALO+NEWNHALO) = NEIGH
              GEN2HALO(NEIGH) = NHALO + NEWNHALO
              DO J=IPE(NEIGH),IPE(NEIGH+1)-1
                IF (TRACE(IW(J)).EQ.NODE) THEN
                  CNT = CNT + 2
                END IF
              END DO
            END IF
          END DO
        END DO
        LAST_LVL_START = NHALO + 1
        NHALO = NHALO + NEWNHALO
      END SUBROUTINE NEIGHBORHOOD
      SUBROUTINE GETHALOGRAPH(HALO,NHALO,N,IW,LW,IPE,IPTRHALO,JCNHALO,
     &     HALOEDGENBR,TRACE,NODE, GEN2HALO)
      INTEGER, INTENT(IN) :: N
      INTEGER,INTENT(IN):: NHALO, NODE
      INTEGER,INTENT(IN):: GEN2HALO(N)
      INTEGER,DIMENSION(NHALO),INTENT(IN) :: HALO
      INTEGER(8), INTENT(IN)              :: LW
      INTEGER(8), INTENT(IN)              :: IPE(N+1)
      INTEGER, INTENT(IN)     :: IW(LW), TRACE(N)
      INTEGER(8),INTENT(IN)   :: HALOEDGENBR
      INTEGER(8), INTENT(OUT) :: IPTRHALO(NHALO+1)
      INTEGER, INTENT(OUT)    :: JCNHALO(HALOEDGENBR)
      INTEGER::I,IPTR_CNT,JCN_CNT,HALOI
      INTEGER(8) :: J, CNT
      CNT = 0
      IPTR_CNT = 2
      JCN_CNT = 1
      IPTRHALO(1) = 1
      DO I=1,NHALO
         HALOI = HALO(I)
         DO J=IPE(HALOI),IPE(HALOI+1)-1
            IF (TRACE(IW(J))==NODE) THEN
               CNT = CNT + 1
               JCNHALO(JCN_CNT) = GEN2HALO(IW(J))
               JCN_CNT = JCN_CNT + 1
            END IF
         END DO
         IPTRHALO(IPTR_CNT) = CNT + 1
         IPTR_CNT = IPTR_CNT + 1
      END DO
      END SUBROUTINE GETHALOGRAPH
      SUBROUTINE GET_GROUPS(NHALO,PARTS,SEP,NSEP,NPARTS,
     &     CUT,NEWSEP,PERM,IPERM)
      INTEGER,INTENT(IN) :: NHALO,NSEP
      INTEGER,DIMENSION(:),INTENT(IN) :: SEP
      INTEGER,POINTER,DIMENSION(:)::PARTS
      INTEGER,POINTER,DIMENSION(:)::CUT,NEWSEP,PERM,
     &   IPERM
      INTEGER,INTENT(INOUT) :: NPARTS 
      INTEGER::I,CNT,NB_PARTS_WITHOUT_SEP_NODE,allocok
      INTEGER,DIMENSION(:),ALLOCATABLE::SIZES
      INTEGER,DIMENSION(:),ALLOCATABLE::PARTPTR
      ALLOCATE(NEWSEP(NSEP),stat=allocok)
      IF(allocok.GT.0) THEN
         write(*,*) "Allocation error in GET_GROUPS"
         CALL MUMPS_ABORT()
      ENDIF
      ALLOCATE(PERM(NSEP),stat=allocok)
      IF(allocok.GT.0) THEN
         write(*,*) "Allocation error in GET_GROUPS"
         CALL MUMPS_ABORT()
      ENDIF
      ALLOCATE(IPERM(NSEP),stat=allocok)
      IF(allocok.GT.0) THEN
         write(*,*) "Allocation error in GET_GROUPS"
         CALL MUMPS_ABORT()
      ENDIF
      ALLOCATE(SIZES(NPARTS),stat=allocok)
      IF(allocok.GT.0) THEN
         write(*,*) "Allocation error in GET_GROUPS"
         CALL MUMPS_ABORT()
      ENDIF
      ALLOCATE(PARTPTR(NPARTS+1),stat=allocok)
      IF(allocok.GT.0) THEN
         write(*,*) "Allocation error in GET_GROUPS"
         CALL MUMPS_ABORT()
      ENDIF
      NB_PARTS_WITHOUT_SEP_NODE = 0
      SIZES = 0
      DO I=1,NSEP
        SIZES(PARTS(I)) = 
     &     SIZES(PARTS(I))+1
      END DO
      PARTPTR(1)=1
      DO I=2,NPARTS+1
        PARTPTR(I) = PARTPTR(I-1) + SIZES(I-1)
        IF (SIZES(I-1)==0) THEN
          NB_PARTS_WITHOUT_SEP_NODE = NB_PARTS_WITHOUT_SEP_NODE + 1
        END IF
      END DO
      ALLOCATE(CUT(NPARTS-NB_PARTS_WITHOUT_SEP_NODE+1),stat=allocok)
      IF(allocok.GT.0) THEN
         write(*,*) "Allocation error in GET_GROUPS"
         CALL MUMPS_ABORT()
      ENDIF
      CUT(1) = 1
      CNT = 2
      DO I=2,NPARTS+1
        IF (SIZES(I-1).NE.0) THEN
          CUT(CNT) = PARTPTR(I) 
          CNT = CNT + 1
        END IF
      END DO
      NPARTS = NPARTS - NB_PARTS_WITHOUT_SEP_NODE 
      CUT(NPARTS+1) = NSEP+1
      DO I=1,NSEP
        NEWSEP(PARTPTR(PARTS(I))) = SEP(I)  
        PERM(PARTPTR(PARTS(I))) = I
        IPERM(I) = PARTPTR(PARTS(I))
        PARTPTR(PARTS(I)) = 
     &     PARTPTR(PARTS(I)) + 1
      END DO
      DEALLOCATE(SIZES,PARTPTR)
      END SUBROUTINE GET_GROUPS
      SUBROUTINE DMUMPS_LR_GROUPING(N, NZ8, NSTEPS, IRN, JCN, FILS,
     &     FRERE_STEPS, DAD_STEPS, NE_STEPS, STEP, NA, LNA,
     &     LRGROUPS, SYM, ICNTL, HALO_DEPTH, GROUP_SIZE, SEP_SIZE,
     &     K38, K20, K60,
     &     IFLAG, IERROR, K264, K265, K482, K472, MAXFRONT, K10, 
     &     K54, LPOK, LP)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, NSTEPS, LNA, SYM,
     &     HALO_DEPTH, SEP_SIZE, GROUP_SIZE
      INTEGER(8), INTENT(IN) :: NZ8
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, INTENT(INOUT) :: K38, K20, K264, K265 
      INTEGER, INTENT(IN)    :: K482, K10, K60, K54
      INTEGER, INTENT(IN)    :: LP
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, POINTER, DIMENSION(:) :: IRN, JCN
      INTEGER, INTENT(IN)    ::  NE_STEPS(:), ICNTL(60)
      INTEGER :: FILS(:), FRERE_STEPS(:), STEP(:),
     &     NA(:), DAD_STEPS(:), LRGROUPS(:)
      INTEGER, INTENT(IN) :: K472, MAXFRONT
      INTEGER :: K482_LOC, K38ou20
      INTEGER :: I, F, PV, NV, NLEAVES, NROOTS, PP, C, NF, NODE,
     &     SYMTRY, NBQD, AD
      INTEGER(8) :: LW, IWFR, NRORM, NIORM
      INTEGER :: LPTR, RPTR, NBGROUPS
      LOGICAL :: FIRST
      INTEGER, ALLOCATABLE, DIMENSION (:) :: POOL, PVS, WORK
      INTEGER, ALLOCATABLE, DIMENSION (:) :: LEN, IW
      INTEGER(8), ALLOCATABLE, DIMENSION (:) :: IPE, IQ
      INTEGER, ALLOCATABLE, DIMENSION (:) :: TRACE, WORKH, GEN2HALO
      INTEGER :: STEP_SCALAPACK_ROOT
      INTEGER :: GROUP_SIZE2, IERR
      LOGICAL :: INPLACE64_GRAPH_COPY
      K38ou20=max(K38,K20)
      IF (K38ou20.GT.0) THEN
       STEP_SCALAPACK_ROOT = STEP(K38ou20)
      ELSE
       STEP_SCALAPACK_ROOT = 0
      ENDIF
      IF((K482.LE.0) .OR. (K482.GT.3)) THEN
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
         K482_LOC = 1
#elif defined(ptscotch) || defined(scotch)
         K482_LOC = 2
#else
         K482_LOC = 3
#endif
      ELSE IF (K482.EQ.1) THEN
#if !defined(parmetis) && !defined(metis) && !defined(parmetis3) && !defined(metis4)
#if defined(ptscotch) || defined(scotch)     
         K482_LOC = 2
#else
         K482_LOC = 3
#endif
#else
         K482_LOC = 1
#endif
      ELSE IF (K482.EQ.2) THEN
#if !defined(ptscotch) && !defined(scotch)
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
         K482_LOC = 1
#else
         K482_LOC = 3
#endif
#else
         K482_LOC = 2
#endif
      ELSE IF (K482.EQ.3) THEN
         K482_LOC = 3
      END IF
      NBGROUPS = 0
      IF (K265.EQ.-1) THEN
       LW = NZ8
      ELSE
        LW = 2_8 * NZ8
      ENDIF
      ALLOCATE(IW(LW), IPE(N+1), LEN(N), IQ(N), 
     &         POOL(NA(1)), PVS(NSTEPS),
     &         STAT=IERR)
      IF (IERR.GT.0) THEN
       IF (LPOK) WRITE(LP,*) " Error allocate integer array of size: ", 
     *    LW+int(N,8)+int(K10*(2*N+1),8)
        IFLAG = -7
        CALL MUMPS_SET_IERROR(LW+int(N,8)+int(K10*(2*N+1),8),IERROR)
        GOTO 500
      ENDIF
      CALL DMUMPS_ANA_GNEW(N, NZ8, IRN(1), JCN(1), IW(1), LW, IPE(1),
     &     LEN(1), IQ(1), LRGROUPS(1), IWFR, NRORM, NIORM,
     &     IFLAG, IERROR,
     &     ICNTL(1) , SYMTRY, SYM, NBQD, AD, K264, K265,.FALSE.,
     &     INPLACE64_GRAPH_COPY)
      IF (K54.EQ.3) THEN
        deallocate(IRN)
        deallocate(JCN)
        NULLIFY(IRN)
        NULLIFY(JCN)
      ENDIF
      IF (allocated(IQ)) DEALLOCATE(IQ)
      LRGROUPS = -1
      NLEAVES = NA(1)
      NROOTS  = NA(2)
      LPTR = 2+NLEAVES
      RPTR = 2+NLEAVES+NROOTS
      DO I = 1, NROOTS
         POOL(I) = NA(2+NLEAVES+I)
      END DO
      PP = NROOTS
      ALLOCATE(WORK(MAXFRONT), TRACE(N), WORKH(N), GEN2HALO(N), 
     &         STAT=IERR)
      IF (IERR.GT.0) THEN
       IF (LPOK) WRITE(LP,*) " Error allocate integer array of size: ", 
     *    3*N+MAXFRONT
        IFLAG = -7
        IERROR = 3*N+MAXFRONT
        RETURN
      ENDIF
      TRACE = 0
      DO WHILE(PP .GT. 0)
         PV = ABS(POOL(PP))
         NODE = STEP(PV)
         FIRST = POOL(PP) .LT. 0
         NV = 0
         F = PV
         DO WHILE(F .GT. 0)
            NV = NV+1
            WORK(NV) = F
            F = FILS(F)
         END DO
         CALL COMPUTE_BLR_VCS(K472, GROUP_SIZE2, GROUP_SIZE, NV)
         IF (NV .GE. GROUP_SIZE2)  THEN
           IF ( (K482_LOC.EQ.3)
     &           .OR.
     &         ( (K60.NE.0).AND.(WORK(1).EQ.K38ou20) )
     &        ) 
     &     THEN
               DO I=1,NV
                  LRGROUPS(WORK(I))=NBGROUPS+1+(I-1)/GROUP_SIZE2
               END DO
               NBGROUPS = NBGROUPS + (NV-1)/GROUP_SIZE2 + 1
            ELSE
              CALL SEP_GROUPING(NV, WORK(1), N, NZ8,
     &              LRGROUPS, NBGROUPS, IW(1), LW, IPE(1), LEN(1),
     &              GROUP_SIZE, HALO_DEPTH, TRACE(1), WORKH(1), NODE,
     &              GEN2HALO(1), K482_LOC, K472, 0, SEP_SIZE, 
     &              K10, LP, LPOK, IFLAG, IERROR)
              IF (IFLAG.LT.0) GOTO 500
            END IF
         ELSE
           IF (NV .GE. SEP_SIZE) THEN
             DO I = 1, NV
             LRGROUPS( WORK(I) ) = (NBGROUPS + 1)
             ENDDO
           ELSE
             DO I = 1, NV
             LRGROUPS( WORK(I) ) = -(NBGROUPS + 1)
             ENDDO
           ENDIF
            NBGROUPS = NBGROUPS + 1
         ENDIF
         CALL MUMPS_UPD_TREE(NV, NSTEPS, N, FIRST, LPTR, RPTR, F,
     &        WORK(1),
     &        FILS, FRERE_STEPS, STEP, DAD_STEPS,
     &        NE_STEPS, NA, LNA, PVS(1), K38ou20, 
     &        STEP_SCALAPACK_ROOT)
          IF (STEP_SCALAPACK_ROOT.GT.0) THEN
           IF (K38.GT.0) THEN
             K38 = K38ou20
           ELSE
             K20 = K38ou20
           ENDIF
          ENDIF
         PP = PP-1
         NF = NE_STEPS(NODE)
         IF(NF .GT. 0) THEN
            PP = PP+1
            POOL(PP) = F
            C = STEP(-F)
            F = FRERE_STEPS(C)
            DO WHILE(F .GT. 0)
               PP = PP+1
               POOL(PP) = F
               C = STEP(F)
               F = FRERE_STEPS(C)
            END DO
         END IF
      END DO
 500  IF (allocated(POOL)) DEALLOCATE(POOL)
      IF (allocated(PVS)) DEALLOCATE(PVS)
      IF (allocated(WORK)) DEALLOCATE(WORK)
      IF (allocated(IPE)) DEALLOCATE(IPE)
      IF (allocated(LEN)) DEALLOCATE(LEN)
      IF (allocated(TRACE)) DEALLOCATE(TRACE)
      IF (allocated(WORKH)) DEALLOCATE(WORKH)
      IF (allocated(GEN2HALO)) DEALLOCATE(GEN2HALO)
      RETURN
      END SUBROUTINE DMUMPS_LR_GROUPING
      SUBROUTINE DMUMPS_LR_GROUPING_NEW(N, NZ8, NSTEPS, IRN, JCN, FILS,
     &     FRERE_STEPS, DAD_STEPS, STEP, NA, LNA, LRGROUPS, 
     &     SYM, ICNTL, HALO_DEPTH, GROUP_SIZE, SEP_SIZE, K38, K20,
     &     K60, IFLAG, IERROR, K264, K265, K482, K472, MAXFRONT, K469, 
     &     K10, K54, LPOK, LP)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, NSTEPS, LNA, SYM,
     &     HALO_DEPTH, SEP_SIZE, GROUP_SIZE
      INTEGER(8), INTENT(IN) :: NZ8
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, INTENT(INOUT) :: K38, K20, K264, K265 
      INTEGER, INTENT(IN)    :: K482, K10, MAXFRONT, K60, K54
      INTEGER, INTENT(IN)    :: LP
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, POINTER, DIMENSION(:)  :: IRN, JCN
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, POINTER :: FILS(:), FRERE_STEPS(:), STEP(:),
     &      NA(:), DAD_STEPS(:), LRGROUPS(:)
      INTEGER, INTENT(IN) :: K472, K469
      INTEGER :: K482_LOC,  K469_LOC, K38ou20
      INTEGER :: I, F, PV, NV, NODE,
     &     SYMTRY, NBQD, AD
      LOGICAL :: PVSCHANGED
      INTEGER(8) :: LW, IWFR, NRORM, NIORM
      INTEGER :: NBGROUPS, NBGROUPS_local
      INTEGER, ALLOCATABLE, DIMENSION (:) :: PVS, WORK
      INTEGER, ALLOCATABLE, DIMENSION (:) :: LEN, IW
      INTEGER(8), ALLOCATABLE, DIMENSION (:) :: IPE, IQ
      INTEGER, ALLOCATABLE, DIMENSION (:) :: TRACE, WORKH, 
     &                                       GEN2HALO
      INTEGER, ALLOCATABLE, DIMENSION (:) :: TRACE_PRV, WORKH_PRV, 
     &                                   GEN2HALO_PRV
      INTEGER :: STEP_SCALAPACK_ROOT
      INTEGER :: GROUP_SIZE2, IERR, OMP_NUM
      INTEGER :: IERR_PRIV 
      LOGICAL :: INPLACE64_GRAPH_COPY
      K38ou20=max(K38,K20)
      IF (K38ou20.GT.0) THEN
       STEP_SCALAPACK_ROOT = STEP(K38ou20)
      ELSE
       STEP_SCALAPACK_ROOT = 0
      ENDIF
      IF((K482.LE.0) .OR. (K482.GT.3)) THEN
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
         K482_LOC = 1
#elif defined(ptscotch) || defined(scotch)
         K482_LOC = 2
#else
         K482_LOC = 3
#endif
      ELSE IF (K482.EQ.1) THEN
#if !defined(parmetis) && !defined(metis) && !defined(parmetis3) && !defined(metis4)
#if defined(ptscotch) || defined(scotch)     
         K482_LOC = 2
#else
         K482_LOC = 3
#endif
#else
         K482_LOC = 1
#endif
      ELSE IF (K482.EQ.2) THEN
#if !defined(ptscotch) && !defined(scotch)
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
         K482_LOC = 1
#else
         K482_LOC = 3
#endif
#else
         K482_LOC = 2
#endif
      ELSE IF (K482.EQ.3) THEN
         K482_LOC = 3
      END IF
      IF (K482_LOC.EQ.2) THEN
        K469_LOC = 1
      ELSE
        K469_LOC = K469
      ENDIF
      NBGROUPS = 0
      LW = 2_8 * NZ8
      ALLOCATE(IW(LW), IPE(N+1), LEN(N), IQ(N), 
     &         PVS(NSTEPS),
     &         STAT=IERR)
      IF (IERR.GT.0) THEN
       IF (LPOK) WRITE(LP,*) " Error allocate integer array of size: ", 
     *    LW+int(N,8)+int(K10*(2*N+1),8)
        IFLAG = -7
        CALL MUMPS_SET_IERROR(LW+int(N,8)+int(K10*(2*N+1),8),IERROR)
        GOTO 501
      ENDIF
      CALL DMUMPS_ANA_GNEW(N, NZ8, IRN(1), JCN(1), IW(1), LW, IPE(1),
     &     LEN(1), IQ(1), LRGROUPS(1), IWFR, NRORM, NIORM,
     &     IFLAG, IERROR,
     &     ICNTL(1) , SYMTRY, SYM, NBQD, AD, K264, K265,.FALSE.,
     &     INPLACE64_GRAPH_COPY)
      IF (K54.EQ.3) THEN
        deallocate(IRN)
        deallocate(JCN)
        NULLIFY(IRN)
        NULLIFY(JCN)
      ENDIF
      IF (allocated(IQ)) DEALLOCATE(IQ)
      LRGROUPS = -1
      IF (K469_LOC.NE.2) THEN
        ALLOCATE(TRACE(N), WORKH(N), GEN2HALO(N), 
     &          STAT=IERR)
        IF (IERR.GT.0) THEN
          IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", 3*N
          IFLAG = -7
          IERROR = 3*N
          GOTO 501
        ENDIF
      ENDIF
      PVSCHANGED = .FALSE.
      OMP_NUM = 1
!$    OMP_NUM = OMP_GET_MAX_THREADS()
      OMP_NUM = min(OMP_NUM,8)
!$OMP PARALLEL PRIVATE(I, NODE, PV, NV, F, GROUP_SIZE2, WORK, IERR_PRIV,
!$OMP&         WORKH_PRV, TRACE_PRV, GEN2HALO_PRV, NBGROUPS_local
!$OMP&         )
!$OMP&         IF (K469_LOC.GT.1) NUM_THREADS(OMP_NUM)
      ALLOCATE(WORK(MAXFRONT), STAT=IERR_PRIV)
      IF (IERR_PRIV.GT.0) THEN
        IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", MAXFRONT
!$OMP ATOMIC WRITE
        IFLAG = -7
!$OMP END ATOMIC
!$OMP ATOMIC WRITE
        IERROR = MAXFRONT
!$OMP END ATOMIC
      ENDIF
      IF (IERR_PRIV .EQ. 0 .AND. K469_LOC.EQ.2) THEN
        ALLOCATE(TRACE_PRV(N), WORKH_PRV(N), GEN2HALO_PRV(N),
     &           STAT=IERR_PRIV)
        IF (IERR_PRIV.GT.0) THEN
          IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", 3*N
!$OMP ATOMIC WRITE
          IFLAG = -7
!$OMP END ATOMIC
!$OMP ATOMIC WRITE
          IERROR = 3*N
!$OMP END ATOMIC
        ENDIF
      ENDIF
!$OMP BARRIER
      IF (IFLAG .LT. 0 ) THEN
        GOTO 500
      ENDIF
      IF (K469_LOC.EQ.2) THEN
        TRACE_PRV = 0
      ELSE
!$OMP SINGLE
        TRACE = 0
!$OMP END SINGLE
      ENDIF
!$OMP DO
      DO I = 1,N
        IF (STEP(I).GT.0) PVS(STEP(I)) = I
      END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC,1)
      DO NODE=NSTEPS,1,-1
        IF (IFLAG.LT.0) CYCLE
        PV = PVS(NODE)
        NV = 0
        F = PV
        DO WHILE(F .GT. 0)
          NV = NV+1
          WORK(NV) = F
          F = FILS(F)
        END DO
        CALL COMPUTE_BLR_VCS(K472, GROUP_SIZE2, GROUP_SIZE, NV)
        IF (NV .GE. GROUP_SIZE2)  THEN
          IF ( (K482_LOC.EQ.3)
     &           .OR.
     &         ( (K60.NE.0).AND.(WORK(1).EQ.K38ou20) )
     &        ) 
     &    THEN
!$OMP CRITICAL(lrgrouping_cri)
            DO I=1,NV
              LRGROUPS(WORK(I))=NBGROUPS+1+(I-1)/GROUP_SIZE2
            END DO
            NBGROUPS = NBGROUPS + (NV-1)/GROUP_SIZE2 + 1
!$OMP END CRITICAL(lrgrouping_cri)
          ELSE
            IF (K469_LOC .EQ. 2) THEN
              CALL SEP_GROUPING(NV, WORK(1), N, NZ8,
     &              LRGROUPS, NBGROUPS, IW(1), LW, IPE(1), LEN(1),
     &              GROUP_SIZE, HALO_DEPTH, TRACE_PRV, WORKH_PRV, 
     &              NODE, GEN2HALO_PRV, K482_LOC, K472, K469_LOC, 
     &              SEP_SIZE, K10, LP, LPOK, IFLAG, IERROR)
            ELSE
              CALL SEP_GROUPING(NV, WORK(1), N, NZ8,
     &              LRGROUPS, NBGROUPS, IW(1), LW, IPE(1), LEN(1),
     &              GROUP_SIZE, HALO_DEPTH, TRACE, WORKH, 
     &              NODE, GEN2HALO, K482_LOC, K472, K469_LOC, 
     &              SEP_SIZE, K10, LP, LPOK, IFLAG, IERROR)
            ENDIF
            IF (IFLAG.LT.0) CYCLE
            PVS(NODE) = WORK(1)
!$OMP ATOMIC WRITE
            PVSCHANGED = .TRUE.
!$OMP END ATOMIC
            STEP(WORK(1)) = ABS(STEP(WORK(1)))
            IF (STEP(WORK(1)).EQ.STEP_SCALAPACK_ROOT) THEN
              IF (K38.GT.0) THEN
                K38 = WORK(1)
              ELSE
                K20 = WORK(1)
              ENDIF
            ENDIF
            DO I=1, NV-1
              STEP(WORK(I+1)) = -STEP(WORK(1))
              IF (FILS(WORK(I)).LE.0) THEN
                FILS(WORK(NV)) = FILS(WORK(I))
              ENDIF
              FILS(WORK(I)) = WORK(I+1)
            ENDDO
          ENDIF
        ELSE
!$OMP ATOMIC CAPTURE
          NBGROUPS = NBGROUPS + 1
          NBGROUPS_local = NBGROUPS
!$OMP END ATOMIC 
          IF (NV .GE. SEP_SIZE) THEN
            DO I = 1, NV
              LRGROUPS( WORK(I) ) = NBGROUPS_local
            ENDDO
          ELSE
            DO I = 1, NV
              LRGROUPS( WORK(I) ) = -NBGROUPS_local
            ENDDO
          ENDIF
        ENDIF
      ENDDO       
!$OMP END DO
      IF (IFLAG.LT.0) GOTO 500
      IF (.NOT.PVSCHANGED) GOTO 500
!$OMP DO
      DO NODE = 1,NSTEPS
        IF(FRERE_STEPS(NODE) .GT. 0) THEN
          FRERE_STEPS(NODE) = PVS(ABS(STEP(FRERE_STEPS(NODE))))
        ELSE IF(FRERE_STEPS(NODE) .LT. 0) THEN
          FRERE_STEPS(NODE) = -PVS(ABS(STEP(DAD_STEPS(NODE))))
        ENDIF
        IF(DAD_STEPS(NODE) .NE. 0) THEN
          DAD_STEPS(NODE) = PVS(ABS(STEP(DAD_STEPS(NODE))))
        END IF
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO
      DO I=3,LNA
        NA(I) = PVS(ABS(STEP(NA(I))))
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO
      DO I=1,N
        IF (FILS(I).LT.0) THEN
          FILS(I) = -PVS(ABS(STEP(-FILS(I))))
        ENDIF
      ENDDO
!$OMP END DO 
 500  CONTINUE
      IF (allocated(WORK)) DEALLOCATE(WORK)
      IF (K469_LOC.EQ.2) THEN
        IF (allocated(TRACE_PRV))    DEALLOCATE(TRACE_PRV)
        IF (allocated(WORKH_PRV))    DEALLOCATE(WORKH_PRV)
        IF (allocated(GEN2HALO_PRV)) DEALLOCATE(GEN2HALO_PRV)
      ENDIF
!$OMP END PARALLEL
 501  CONTINUE
      IF (K469_LOC.NE.2) THEN
        IF (allocated(TRACE))    DEALLOCATE(TRACE)
        IF (allocated(WORKH))    DEALLOCATE(WORKH)
        IF (allocated(GEN2HALO)) DEALLOCATE(GEN2HALO)
      ENDIF
      IF (allocated(PVS)) DEALLOCATE(PVS)
      IF (allocated(IPE)) DEALLOCATE(IPE)
      IF (allocated(LEN)) DEALLOCATE(LEN)
      RETURN
      END SUBROUTINE DMUMPS_LR_GROUPING_NEW
      SUBROUTINE DMUMPS_AB_LR_GROUPING(N, MAPCOL, SIZEMAPCOL,
     &     NSTEPS, LUMAT, FILS,
     &     FRERE_STEPS, DAD_STEPS, STEP, NA, LNA, LRGROUPS, 
     &     SIZEOFBLOCKS, SYM, ICNTL, HALO_DEPTH, GROUP_SIZE, 
     &     SEP_SIZE, K38, K20,
     &     K60, IFLAG, IERROR, K264, K265, K482, K472, MAXFRONT, K469, 
     &     K10, K54, LPOK, LP, MYID, COMM)
      IMPLICIT NONE
      INTEGER, INTENT(IN)     :: MYID, COMM
      TYPE(LMATRIX_T)         :: LUMAT
      INTEGER, INTENT(IN)    :: N, NSTEPS, LNA, SYM,
     &     HALO_DEPTH, SEP_SIZE, GROUP_SIZE
      INTEGER, INTENT(IN)    :: SIZEMAPCOL
      INTEGER, INTENT(IN)    :: MAPCOL(SIZEMAPCOL)
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, INTENT(INOUT) :: K38, K20, K264, K265 
      INTEGER, INTENT(IN)    :: K482, K10, MAXFRONT, K60, K54
      INTEGER, INTENT(IN)    :: LP
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, POINTER :: FILS(:), FRERE_STEPS(:), STEP(:),
     &      NA(:), DAD_STEPS(:), LRGROUPS(:)
      INTEGER, INTENT(IN)  :: SIZEOFBLOCKS(N)
      INTEGER, INTENT(IN) :: K472, K469
      INTEGER :: K482_LOC,  K469_LOC, K38ou20
      INTEGER :: I, F, PV, NV, NVEXPANDED, NODE
      DOUBLE PRECISION    :: COMPRESS_RATIO
      LOGICAL :: PVSCHANGED
      INTEGER :: NBGROUPS, NBGROUPS_local
      INTEGER, ALLOCATABLE, DIMENSION (:) :: PVS, WORK
      INTEGER, ALLOCATABLE, DIMENSION (:) :: TRACE, WORKH, 
     &                                       GEN2HALO
      INTEGER, ALLOCATABLE, DIMENSION (:) :: TRACE_PRV, WORKH_PRV, 
     &                                   GEN2HALO_PRV
      INTEGER :: STEP_SCALAPACK_ROOT
      INTEGER :: GROUP_SIZE2, IERR, OMP_NUM
      INTEGER :: IERR_PRIV 
      LOGICAL :: MAPCOL_PROVIDED
      MAPCOL_PROVIDED = (MAPCOL(1).GE.0)
      K38ou20=max(K38,K20)
      IF (K38ou20.GT.0) THEN
       STEP_SCALAPACK_ROOT = STEP(K38ou20)
      ELSE
       STEP_SCALAPACK_ROOT = 0
      ENDIF
      IF((K482.LE.0) .OR. (K482.GT.3)) THEN
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
         K482_LOC = 1
#elif defined(ptscotch) || defined(scotch)
         K482_LOC = 2
#else
         K482_LOC = 3
#endif
      ELSE IF (K482.EQ.1) THEN
#if !defined(parmetis) && !defined(metis) && !defined(parmetis3) && !defined(metis4)
#if defined(ptscotch) || defined(scotch)     
         K482_LOC = 2
#else
         K482_LOC = 3
#endif
#else
         K482_LOC = 1
#endif
      ELSE IF (K482.EQ.2) THEN
#if !defined(ptscotch) && !defined(scotch)
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
         K482_LOC = 1
#else
         K482_LOC = 3
#endif
#else
         K482_LOC = 2
#endif
      ELSE IF (K482.EQ.3) THEN
         K482_LOC = 3
      END IF
      IF (K482_LOC.EQ.2) THEN
        K469_LOC = 1
      ELSE
        K469_LOC = K469
      ENDIF
      NBGROUPS = 0
      ALLOCATE( PVS(NSTEPS), STAT=IERR)
      IF (IERR.GT.0) THEN
        IFLAG = -7
        IERROR = NSTEPS
        IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", IERROR
        GOTO 501
      ENDIF
      LRGROUPS = -1
      IF (K469_LOC.NE.2) THEN
        ALLOCATE(TRACE(N), WORKH(N), GEN2HALO(N), 
     &          STAT=IERR)
        IF (IERR.GT.0) THEN
          IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", 3*N
          IFLAG = -7
          IERROR = 3*N
          GOTO 501
        ENDIF
      ENDIF
      PVSCHANGED = .FALSE.
      OMP_NUM = 1
!$    OMP_NUM = OMP_GET_MAX_THREADS()
      OMP_NUM = min(OMP_NUM,8)
!$OMP PARALLEL PRIVATE(I, NODE, PV, NV, F, GROUP_SIZE2, WORK, IERR_PRIV,
!$OMP&         WORKH_PRV, TRACE_PRV, GEN2HALO_PRV, NBGROUPS_local,
!$OMP&         NVEXPANDED, COMPRESS_RATIO
!$OMP&         )
!$OMP&         IF (K469_LOC.GT.1) NUM_THREADS(OMP_NUM)
      ALLOCATE(WORK(MAXFRONT), STAT=IERR_PRIV)
      IF (IERR_PRIV.GT.0) THEN
        IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", MAXFRONT
!$OMP ATOMIC WRITE
        IFLAG = -7
!$OMP END ATOMIC
!$OMP ATOMIC WRITE
        IERROR = MAXFRONT
!$OMP END ATOMIC
      ENDIF
      IF (IERR_PRIV .EQ. 0 .AND. K469_LOC.EQ.2) THEN
        ALLOCATE(TRACE_PRV(N), WORKH_PRV(N), GEN2HALO_PRV(N),
     &           STAT=IERR_PRIV)
        IF (IERR_PRIV.GT.0) THEN
          IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", 3*N
!$OMP ATOMIC WRITE
          IFLAG = -7
!$OMP ATOMIC WRITE
          IERROR = 3*N
        ENDIF
      ENDIF
!$OMP BARRIER
      IF (IFLAG .LT. 0 ) THEN
        GOTO 500
      ENDIF
      IF (K469_LOC.EQ.2) THEN
        TRACE_PRV = 0
      ELSE
!$OMP SINGLE
        TRACE = 0
!$OMP END SINGLE
      ENDIF
!$OMP DO
      DO I = 1,N
        IF (STEP(I).GT.0) PVS(STEP(I)) = I
      END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC,1)
      DO NODE=NSTEPS,1,-1
        IF (IFLAG.LT.0) CYCLE
          IF (MAPCOL_PROVIDED) THEN
            IF (MAPCOL(NODE).NE.MYID) THEN
             PVS(NODE) = -999  
             CYCLE
            ENDIF
          ENDIF
        PV = PVS(NODE)
        NV         = 0
        NVEXPANDED = 0
        F = PV
        DO WHILE(F .GT. 0)
          NV = NV+1
          NVEXPANDED = NVEXPANDED+SIZEOFBLOCKS(F)
          WORK(NV) = F
          F = FILS(F)
        END DO
        COMPRESS_RATIO = dble(NVEXPANDED)/dble(NV)
        CALL COMPUTE_BLR_VCS(K472, GROUP_SIZE2, GROUP_SIZE, NVEXPANDED)
        IF (NVEXPANDED .GE. GROUP_SIZE2)  THEN
          IF ( (K482_LOC.EQ.3)
     &           .OR.
     &         ( (K60.NE.0).AND.(WORK(1).EQ.K38ou20) )
     &        ) 
     &    THEN
            GROUP_SIZE2 = max(int(dble(GROUP_SIZE2)/COMPRESS_RATIO), 1)
!$OMP CRITICAL(lrgrouping_cri)
            DO I=1,NV
              LRGROUPS(WORK(I))=NBGROUPS+1+(I-1)/GROUP_SIZE2
            END DO
            NBGROUPS = NBGROUPS + (NV-1)/GROUP_SIZE2 + 1
!$OMP END CRITICAL(lrgrouping_cri)
          ELSE
            IF (K469_LOC .EQ. 2) THEN
              CALL SEP_GROUPING_AB(NV, NVEXPANDED, WORK(1), N, 
     &              LRGROUPS, NBGROUPS, LUMAT, SIZEOFBLOCKS,
     &              GROUP_SIZE, HALO_DEPTH, TRACE_PRV, WORKH_PRV, 
     &              NODE, GEN2HALO_PRV, K482_LOC, K472, K469_LOC, 
     &              SEP_SIZE, K10, LP, LPOK, IFLAG, IERROR)
            ELSE
              CALL SEP_GROUPING_AB(NV, NVEXPANDED, WORK(1), N, 
     &              LRGROUPS, NBGROUPS, LUMAT,  SIZEOFBLOCKS,
     &              GROUP_SIZE, HALO_DEPTH, TRACE, WORKH, 
     &              NODE, GEN2HALO, K482_LOC, K472, K469_LOC, 
     &              SEP_SIZE, K10, LP, LPOK, IFLAG, IERROR)
            ENDIF
            IF (IFLAG.LT.0) CYCLE
            PVS(NODE) = WORK(1)
!$OMP ATOMIC WRITE
            PVSCHANGED = .TRUE.
!$OMP END ATOMIC
            STEP(WORK(1)) = ABS(STEP(WORK(1)))
            IF (STEP(WORK(1)).EQ.STEP_SCALAPACK_ROOT) THEN
              IF (K38.GT.0) THEN
                K38 = WORK(1)
              ELSE
                K20 = WORK(1)
              ENDIF
            ENDIF
            DO I=1, NV-1
              STEP(WORK(I+1)) = -STEP(WORK(1))
              IF (FILS(WORK(I)).LE.0) THEN
                FILS(WORK(NV)) = FILS(WORK(I))
              ENDIF
              FILS(WORK(I)) = WORK(I+1)
            ENDDO
          ENDIF
        ELSE
!$OMP ATOMIC CAPTURE
          NBGROUPS = NBGROUPS + 1
          NBGROUPS_local = NBGROUPS
!$OMP END ATOMIC 
          IF (NVEXPANDED .GE. SEP_SIZE) THEN
            DO I = 1, NV
              LRGROUPS( WORK(I) ) = NBGROUPS_local
            ENDDO
          ELSE
            DO I = 1, NV
              LRGROUPS( WORK(I) ) = -NBGROUPS_local
            ENDDO
          ENDIF
        ENDIF
      ENDDO       
!$OMP END DO
      IF (IFLAG.LT.0) GOTO 500
      IF (.NOT.PVSCHANGED) GOTO 500
!$OMP DO
      DO NODE = 1,NSTEPS
        IF(FRERE_STEPS(NODE) .GT. 0) THEN
          FRERE_STEPS(NODE) = PVS(ABS(STEP(FRERE_STEPS(NODE))))
        ELSE IF(FRERE_STEPS(NODE) .LT. 0) THEN
          FRERE_STEPS(NODE) = -PVS(ABS(STEP(DAD_STEPS(NODE))))
        ENDIF
        IF(DAD_STEPS(NODE) .NE. 0) THEN
          DAD_STEPS(NODE) = PVS(ABS(STEP(DAD_STEPS(NODE))))
        END IF
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO
      DO I=3,LNA
        NA(I) = PVS(ABS(STEP(NA(I))))
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO
      DO I=1,N
        IF (FILS(I).LT.0) THEN
          FILS(I) = -PVS(ABS(STEP(-FILS(I))))
        ENDIF
      ENDDO
!$OMP END DO 
 500  CONTINUE
      IF (allocated(WORK)) DEALLOCATE(WORK)
      IF (K469_LOC.EQ.2) THEN
        IF (allocated(TRACE_PRV))    DEALLOCATE(TRACE_PRV)
        IF (allocated(WORKH_PRV))    DEALLOCATE(WORKH_PRV)
        IF (allocated(GEN2HALO_PRV)) DEALLOCATE(GEN2HALO_PRV)
      ENDIF
!$OMP END PARALLEL
 501  CONTINUE
      IF (K469_LOC.NE.2) THEN
        IF (allocated(TRACE))    DEALLOCATE(TRACE)
        IF (allocated(WORKH))    DEALLOCATE(WORKH)
        IF (allocated(GEN2HALO)) DEALLOCATE(GEN2HALO)
      ENDIF
      IF (allocated(PVS)) DEALLOCATE(PVS)
      RETURN
      END SUBROUTINE DMUMPS_AB_LR_GROUPING
      SUBROUTINE DMUMPS_AB_LR_MPI_GROUPING(
     &     N, MAPCOL, SIZEMAPCOL,
     &     NSTEPS, LUMAT, FILS,
     &     FRERE_STEPS, DAD_STEPS, STEP, NA, LNA, LRGROUPS, 
     &     SIZEOFBLOCKS, SYM, ICNTL, HALO_DEPTH, GROUP_SIZE, 
     &     SEP_SIZE, K38, K20,
     &     K60, IFLAG, IERROR, K264, K265, K482, K472, MAXFRONT, K469, 
     &     K10, K54, LPOK, LP, 
     &     COMM, MYID, NPROCS 
     &     )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER  :: IERR_MPI, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER, INTENT(IN)     :: MYID, COMM, NPROCS
      TYPE(LMATRIX_T)         :: LUMAT
      INTEGER, INTENT(IN)    :: N, NSTEPS, LNA, SYM,
     &     HALO_DEPTH, SEP_SIZE, GROUP_SIZE
      INTEGER, INTENT(IN)    :: SIZEMAPCOL
      INTEGER, INTENT(IN)    :: MAPCOL(SIZEMAPCOL)
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER, INTENT(INOUT) :: K38, K20, K264, K265 
      INTEGER, INTENT(IN)    :: K482, K10, MAXFRONT, K60, K54
      INTEGER, INTENT(IN)    :: LP
      LOGICAL, INTENT(IN)    :: LPOK
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INTEGER, POINTER :: FILS(:), FRERE_STEPS(:), STEP(:),
     &      NA(:), DAD_STEPS(:), LRGROUPS(:)
      INTEGER, INTENT(IN)  :: SIZEOFBLOCKS(N)
      INTEGER, INTENT(IN) :: K472, K469
      INTEGER :: K482_LOC,  K469_LOC, K38ou20
      INTEGER :: I, F, PV, NV, NVEXPANDED, NODE
      DOUBLE PRECISION    :: COMPRESS_RATIO
      LOGICAL :: PVSCHANGED
      INTEGER :: PVSCHANGED_INT, PVSCHANGED_INT_GLOB, IPROC
      INTEGER :: NBGROUPS, NBGROUPS_local 
      INTEGER, ALLOCATABLE, DIMENSION (:) :: PVS, WORK
      INTEGER :: NBGROUPS_sent
      INTEGER :: NBNODES_LOC, SIZE_SENT, ISHIFT, 
     &           MSGSOU, ILOOP
      INTEGER, ALLOCATABLE, DIMENSION (:) :: TRACE, WORKH, 
     &                                       GEN2HALO
      INTEGER, ALLOCATABLE, DIMENSION (:) :: TRACE_PRV, WORKH_PRV, 
     &                                   GEN2HALO_PRV
      INTEGER :: STEP_SCALAPACK_ROOT
      INTEGER :: GROUP_SIZE2, IERR, OMP_NUM
      INTEGER :: IERR_PRIV 
      LOGICAL :: MAPCOL_PROVIDED
      MAPCOL_PROVIDED = (MAPCOL(1).GE.0)
      K38ou20=max(K38,K20)
      IF (K38ou20.GT.0) THEN
       STEP_SCALAPACK_ROOT = STEP(K38ou20)
      ELSE
       STEP_SCALAPACK_ROOT = 0
      ENDIF
      IF (MAPCOL_PROVIDED) THEN
       CALL MPI_BCAST( FILS(1), N, MPI_INTEGER,
     &     MASTER, COMM, IERR )
      ENDIF
      IF((K482.LE.0) .OR. (K482.GT.3)) THEN
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
         K482_LOC = 1
#elif defined(ptscotch) || defined(scotch)
         K482_LOC = 2
#else
         K482_LOC = 3
#endif
      ELSE IF (K482.EQ.1) THEN
#if !defined(parmetis) && !defined(metis) && !defined(parmetis3) && !defined(metis4)
#if defined(ptscotch) || defined(scotch)     
         K482_LOC = 2
#else
         K482_LOC = 3
#endif
#else
         K482_LOC = 1
#endif
      ELSE IF (K482.EQ.2) THEN
#if !defined(ptscotch) && !defined(scotch)
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
         K482_LOC = 1
#else
         K482_LOC = 3
#endif
#else
         K482_LOC = 2
#endif
      ELSE IF (K482.EQ.3) THEN
         K482_LOC = 3
      END IF
      IF (K482_LOC.EQ.2) THEN
        K469_LOC = 1
      ELSE
        K469_LOC = K469
      ENDIF
      NBGROUPS = 0
      ALLOCATE( PVS(NSTEPS), STAT=IERR)
      IF (IERR.GT.0) THEN
        IFLAG = -7
        IERROR = NSTEPS
        IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", IERROR
        GOTO 491
      ENDIF
      LRGROUPS = -1
      IF (K469_LOC.NE.2) THEN
        ALLOCATE(TRACE(N), WORKH(N), GEN2HALO(N), 
     &          STAT=IERR)
        IF (IERR.GT.0) THEN
          IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", 3*N
          IFLAG = -7
          IERROR = 3*N
          GOTO 491
        ENDIF
      ENDIF
491   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), IFLAG,
     &     COMM, MYID )
      IF (IFLAG.LT.0) GOTO 501
      PVSCHANGED = .FALSE.
      OMP_NUM = 1
!$    OMP_NUM = OMP_GET_MAX_THREADS()
      OMP_NUM = min(OMP_NUM,8)
!$OMP PARALLEL PRIVATE(I, NODE, PV, NV, F, GROUP_SIZE2, WORK, IERR_PRIV,
!$OMP&         WORKH_PRV, TRACE_PRV, GEN2HALO_PRV, NBGROUPS_local,
!$OMP&         NVEXPANDED, COMPRESS_RATIO, IPROC
!$OMP&         )
!$OMP&         IF (K469_LOC.GT.1) NUM_THREADS(OMP_NUM)
      ALLOCATE(WORK(2*MAXFRONT+1), STAT=IERR_PRIV)
      IF (IERR_PRIV.GT.0) THEN
        IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", 2*MAXFRONT+1
!$OMP ATOMIC WRITE
        IFLAG = -7
!$OMP END ATOMIC
!$OMP ATOMIC WRITE
        IERROR = 2*MAXFRONT+1
!$OMP END ATOMIC
      ENDIF
      IF (IERR_PRIV .EQ. 0 .AND. K469_LOC.EQ.2) THEN
        ALLOCATE(TRACE_PRV(N), WORKH_PRV(N), GEN2HALO_PRV(N),
     &           STAT=IERR_PRIV)
        IF (IERR_PRIV.GT.0) THEN
          IF (LPOK) WRITE(LP,*) " Error allocate integer array of ", 
     *    "size: ", 3*N
!$OMP ATOMIC WRITE
          IFLAG = -7
!$OMP END ATOMIC
!$OMP ATOMIC WRITE
          IERROR = 3*N
!$OMP END ATOMIC
        ENDIF
      ENDIF
!$OMP BARRIER
      IF (IFLAG .LT. 0 ) THEN
        GOTO 498
      ENDIF
      IF (K469_LOC.EQ.2) THEN
        TRACE_PRV = 0
      ELSE
!$OMP SINGLE
        TRACE = 0
!$OMP END SINGLE
      ENDIF
!$OMP DO
      DO I = 1,N
        IF (STEP(I).GT.0) PVS(STEP(I)) = I
      END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC,1)
      DO NODE=NSTEPS,1,-1
        IF (IFLAG.LT.0) CYCLE
        IF (MAPCOL_PROVIDED) THEN
            IPROC = MAPCOL(NODE)
            IF (IPROC.NE.MYID) THEN
             PVS(NODE) = -999  
             CYCLE
            ENDIF
        ENDIF
        PV = PVS(NODE)
        NV         = 0
        NVEXPANDED = 0
        F = PV
        DO WHILE(F .GT. 0)
          NV = NV+1
          NVEXPANDED = NVEXPANDED+SIZEOFBLOCKS(F)
          WORK(NV) = F
          F = FILS(F)
        END DO
        COMPRESS_RATIO = dble(NVEXPANDED)/dble(NV)
        CALL COMPUTE_BLR_VCS(K472, GROUP_SIZE2, GROUP_SIZE, NVEXPANDED)
        IF (NVEXPANDED .GE. GROUP_SIZE2)  THEN
          IF ( (K482_LOC.EQ.3)
     &           .OR.
     &         ( (K60.NE.0).AND.(WORK(1).EQ.K38ou20) )
     &        ) 
     &    THEN
            GROUP_SIZE2 = max(int(dble(GROUP_SIZE2)/COMPRESS_RATIO), 1)
!$OMP CRITICAL(lrgrouping_cri)
            DO I=1,NV
              LRGROUPS(WORK(I))=NBGROUPS+1+(I-1)/GROUP_SIZE2
            END DO
            NBGROUPS = NBGROUPS + (NV-1)/GROUP_SIZE2 + 1
!$OMP END CRITICAL(lrgrouping_cri)
          ELSE
            IF (K469_LOC .EQ. 2) THEN
              CALL SEP_GROUPING_AB(NV, NVEXPANDED, WORK(1), N, 
     &              LRGROUPS, NBGROUPS, LUMAT, SIZEOFBLOCKS,
     &              GROUP_SIZE, HALO_DEPTH, TRACE_PRV, WORKH_PRV, 
     &              NODE, GEN2HALO_PRV, K482_LOC, K472, K469_LOC, 
     &              SEP_SIZE, K10, LP, LPOK, IFLAG, IERROR)
            ELSE
              CALL SEP_GROUPING_AB(NV, NVEXPANDED, WORK(1), N, 
     &              LRGROUPS, NBGROUPS, LUMAT, SIZEOFBLOCKS,
     &              GROUP_SIZE, HALO_DEPTH, TRACE, WORKH, 
     &              NODE, GEN2HALO, K482_LOC, K472, K469_LOC, 
     &              SEP_SIZE, K10, LP, LPOK, IFLAG, IERROR)
            ENDIF
            IF (IFLAG.LT.0) CYCLE
            PVS(NODE) = WORK(1)
!$OMP ATOMIC WRITE
            PVSCHANGED = .TRUE.
!$OMP END ATOMIC
            STEP(WORK(1)) = ABS(STEP(WORK(1)))
            IF (STEP(WORK(1)).EQ.STEP_SCALAPACK_ROOT) THEN
              IF (K38.GT.0) THEN
                K38 = WORK(1)
              ELSE
                K20 = WORK(1)
              ENDIF
            ENDIF
            DO I=1, NV-1
              STEP(WORK(I+1)) = -STEP(WORK(1))
              IF (FILS(WORK(I)).LE.0) THEN
                FILS(WORK(NV)) = FILS(WORK(I))
              ENDIF
              FILS(WORK(I)) = WORK(I+1)
            ENDDO
          ENDIF
        ELSE
!$OMP ATOMIC CAPTURE
          NBGROUPS = NBGROUPS + 1
          NBGROUPS_local = NBGROUPS
!$OMP END ATOMIC 
          IF (NVEXPANDED .GE. SEP_SIZE) THEN
            DO I = 1, NV
              LRGROUPS( WORK(I) ) = NBGROUPS_local
            ENDDO
          ELSE
            DO I = 1, NV
              LRGROUPS( WORK(I) ) = -NBGROUPS_local
            ENDDO
          ENDIF
        ENDIF
      ENDDO       
!$OMP END DO
 498  CONTINUE
!$OMP MASTER
      CALL MUMPS_PROPINFO( ICNTL(1), IFLAG,
     &     COMM, MYID )
!$OMP END MASTER
!$OMP BARRIER
      IF (IFLAG.LT.0) GOTO 500
      IF (K469_LOC.EQ.2) THEN
        IF (allocated(TRACE_PRV))    DEALLOCATE(TRACE_PRV)
        IF (allocated(WORKH_PRV))    DEALLOCATE(WORKH_PRV)
        IF (allocated(GEN2HALO_PRV)) DEALLOCATE(GEN2HALO_PRV)
      ENDIF
!$OMP MASTER
      IF (K469_LOC.NE.2) THEN
          IF (allocated(WORKH)) DEALLOCATE(WORKH)
          IF (allocated(TRACE)) DEALLOCATE(TRACE)
          IF (allocated(GEN2HALO)) DEALLOCATE(GEN2HALO)
      ENDIF
!$OMP END MASTER
      IF (.NOT.MAPCOL_PROVIDED) THEN 
!$OMP MASTER
        IF (PVSCHANGED) THEN
         PVSCHANGED_INT_GLOB = 1
        ELSE
         PVSCHANGED_INT_GLOB = 0
        ENDIF
!$OMP END MASTER
      ELSE
!$OMP MASTER
        IF (PVSCHANGED) THEN
           PVSCHANGED_INT = 1
        ELSE
           PVSCHANGED_INT = 0
        ENDIF
        CALL MPI_ALLREDUCE( PVSCHANGED_INT, PVSCHANGED_INT_GLOB, 1, 
     &                 MPI_INTEGER,
     &                 MPI_MAX, COMM, IERR_MPI )
        PVSCHANGED_INT_GLOB = 1
        IF (PVSCHANGED_INT_GLOB.NE.0) THEN
          IF (NPROCS.GT.1) THEN 
            ALLOCATE(WORKH(2*N+3*NSTEPS+1), STAT=IERR_PRIV)
            IF (IERR_PRIV.GT.0) THEN
              IF (LPOK) WRITE(LP,*) 
     &        " Error allocate integer array of ", 
     &        "size: ", 2*MAXFRONT+1
              IFLAG = -7
              IERROR = 2*N+3*NSTEPS+1
            ENDIF
            CALL MUMPS_PROPINFO( ICNTL(1), IFLAG,
     &           COMM, MYID )
            IF (IFLAG.LT.0) GOTO 499
            IF (MYID.EQ.MASTER) THEN
              IPROC = 0
              DO WHILE (IPROC.NE.NPROCS-1) 
                IPROC = IPROC + 1
                CALL MPI_RECV( NBNODES_LOC, 1, MPI_INTEGER, 
     &                MPI_ANY_SOURCE,
     &                GROUPING, COMM, STATUS, IERR )
                MSGSOU = STATUS( MPI_SOURCE )
                IF (NBNODES_LOC.EQ.0) THEN
                  CYCLE
                ENDIF
                CALL MPI_RECV( NBGROUPS_sent, 1, MPI_INTEGER, 
     &                 MSGSOU, GROUPING, COMM, STATUS, IERR )
                CALL MPI_RECV( SIZE_SENT, 1, MPI_INTEGER, 
     &                 MSGSOU, GROUPING, COMM, STATUS, IERR )
                CALL MPI_RECV( WORKH, SIZE_SENT, MPI_INTEGER, 
     &                 MSGSOU, GROUPING, COMM, STATUS, IERR )
                ISHIFT = 0
                DO ILOOP=1, NBNODES_LOC
                  ISHIFT = ISHIFT+1
                  NODE   = WORKH (ISHIFT)
                  ISHIFT = ISHIFT+1
                  NV     = WORKH(ISHIFT)
                  PVS(NODE)           = WORKH(ISHIFT+1)
                  STEP(WORKH(ISHIFT+1)) = NODE
                  IF (STEP(WORKH(ISHIFT+1)).EQ.STEP_SCALAPACK_ROOT) THEN
                    IF (K38.GT.0) THEN
                      K38 = WORKH(ISHIFT+1)
                    ELSE
                      K20 = WORKH(ISHIFT+1)
                    END IF
                  END IF
                  DO I=2, NV
                    STEP(WORKH(I+ISHIFT)) = -NODE
                  END DO
                  DO I=1, NV
                   FILS(WORKH(I+ISHIFT))     = WORKH(I+1+ISHIFT)
                   IF (WORKH(NV+1+I+ISHIFT).LT.0) THEN
                     LRGROUPS(WORKH(I+ISHIFT)) = 
     &                            - NBGROUPS + WORKH(NV+1+I+ISHIFT)
                   ELSE
                     LRGROUPS(WORKH(I+ISHIFT)) = 
     &                            NBGROUPS + WORKH(NV+1+I+ISHIFT)
                   END IF
                  END DO
                  ISHIFT = ISHIFT + 2*NV +1
                END DO
                NBGROUPS = NBGROUPS + NBGROUPS_sent
              ENDDO
            ELSE
              NBNODES_LOC = 0
              SIZE_SENT   = 0
              ISHIFT      = 0
              DO NODE = 1,NSTEPS
               IPROC = MAPCOL(NODE)
               IF (IPROC.EQ.MYID) THEN
                 NBNODES_LOC   = NBNODES_LOC + 1
                 ISHIFT        = ISHIFT +1
                 WORKH(ISHIFT) = NODE
                 ISHIFT        = ISHIFT +1
                 NV            = 0
                 F             =  PVS(NODE)
                 DO WHILE (F.GT.0) 
                  NV               = NV + 1
                  WORKH(NV+ISHIFT) = F
                  F                = FILS(F)
                 ENDDO
                 WORKH(ISHIFT) = NV
                 WORKH(NV+1+ISHIFT) = F
                 DO I=1, NV
                  WORKH(NV+1+I+ISHIFT) = LRGROUPS(WORKH(I+ISHIFT))
                 ENDDO
                 ISHIFT = ISHIFT + 2*NV+1
               ENDIF
              ENDDO
              SIZE_SENT = ISHIFT
              CALL MPI_SEND( NBNODES_LOC, 1, MPI_INTEGER, MASTER,
     &               GROUPING, COMM, IERR )
              IF (NBNODES_LOC.GT.0) THEN
                CALL MPI_SEND( NBGROUPS, 1, MPI_INTEGER, MASTER,
     &               GROUPING, COMM, IERR )
                CALL MPI_SEND( SIZE_SENT, 1, MPI_INTEGER, MASTER,
     &                 GROUPING, COMM, IERR )
                CALL MPI_SEND( WORKH, SIZE_SENT, MPI_INTEGER, MASTER,
     &               GROUPING, COMM, IERR )
              ENDIF
            ENDIF
          ENDIF
        ENDIF 
 499  CONTINUE
!$OMP END MASTER
      ENDIF
!$OMP BARRIER
      IF (IFLAG.LT.0) GOTO 500
      IF (MYID.EQ.MASTER) THEN
        IF (PVSCHANGED_INT_GLOB.EQ.0) GOTO 500
!$OMP   DO
        DO NODE = 1,NSTEPS
          IF(FRERE_STEPS(NODE) .GT. 0) THEN
            FRERE_STEPS(NODE) = PVS(ABS(STEP(FRERE_STEPS(NODE))))
          ELSE IF(FRERE_STEPS(NODE) .LT. 0) THEN
            FRERE_STEPS(NODE) = -PVS(ABS(STEP(DAD_STEPS(NODE))))
          ENDIF
          IF(DAD_STEPS(NODE) .NE. 0) THEN
            DAD_STEPS(NODE) = PVS(ABS(STEP(DAD_STEPS(NODE))))
          END IF
        ENDDO
!$OMP   END DO NOWAIT
!$OMP   DO
        DO I=3,LNA
          NA(I) = PVS(ABS(STEP(NA(I))))
        ENDDO
!$OMP   END DO NOWAIT
!$OMP   DO
        DO I=1,N
          IF (FILS(I).LT.0) THEN
            FILS(I) = -PVS(ABS(STEP(-FILS(I))))
          ENDIF
        ENDDO
!$OMP   END DO 
      ENDIF 
 500  CONTINUE
      IF (allocated(WORK)) DEALLOCATE(WORK)
      IF (K469_LOC.EQ.2) THEN
        IF (allocated(TRACE_PRV))    DEALLOCATE(TRACE_PRV)
        IF (allocated(WORKH_PRV))    DEALLOCATE(WORKH_PRV)
        IF (allocated(GEN2HALO_PRV)) DEALLOCATE(GEN2HALO_PRV)
      ENDIF
!$OMP END PARALLEL
 501  CONTINUE
      IF (K469_LOC.NE.2) THEN
        IF (allocated(TRACE))    DEALLOCATE(TRACE)
        IF (allocated(WORKH))    DEALLOCATE(WORKH)
        IF (allocated(GEN2HALO)) DEALLOCATE(GEN2HALO)
      ENDIF
      IF (allocated(PVS)) DEALLOCATE(PVS)
      RETURN
      END SUBROUTINE DMUMPS_AB_LR_MPI_GROUPING
      END MODULE DMUMPS_ANA_LR
