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
      MODULE DMUMPS_FAC_FRONT_AUX_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC_H(NFRONT,NASS,IW,LIW,A,LA,
     &   INOPV,NOFFW,
     &     DET_EXPW, DET_MANTW, DET_SIGNW,
     &     IOLDPS,POSELT,UU,SEUIL,KEEP, KEEP8, DKEEP,
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &     PP_LastPIVRPTRFilled_U,MAXFROMN,IS_MAXFROMN_AVAIL, 
     &     Inextpiv, OOC_EFFECTIVE_ON_FRONT, NVSCHUR
     &)
!$    USE OMP_LIB
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER NFRONT,NASS,LIW,INOPV
      INTEGER(8) :: LA
      INTEGER    :: KEEP(500)
      INTEGER(8) :: KEEP8(150)
      DOUBLE PRECISION       :: DKEEP(230)
      DOUBLE PRECISION UU, SEUIL
      DOUBLE PRECISION A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION, intent(in) :: MAXFROMN
      LOGICAL, intent(inout) :: IS_MAXFROMN_AVAIL
      INTEGER, intent(inout) :: Inextpiv
      LOGICAL, intent(in)    :: OOC_EFFECTIVE_ON_FRONT
      INTEGER, intent(in)    :: NVSCHUR
      DOUBLE PRECISION AMROW
      DOUBLE PRECISION RMAX
      DOUBLE PRECISION  SWOP
      INTEGER(8) :: APOS, POSELT
      INTEGER(8) :: J1, J2, J3_8, JJ, IDIAG
      INTEGER(8) :: J1_ini
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS
      INTEGER NPIV,IPIV,IPIV_SHIFT
      INTEGER, intent(inout) :: NOFFW
      INTEGER, intent(inout) :: DET_EXPW, DET_SIGNW
      DOUBLE PRECISION, intent(inout) :: DET_MANTW
      INTEGER J, J3
      INTEGER NPIVP1,JMAX,ISW,ISWPS1
      INTEGER ISWPS2,KSW,XSIZE
      INTEGER I_PIVRPTR_L, I_PIVR_L, NBPANELS_L
      INTEGER I_PIVRPTR_U, I_PIVR_U, NBPANELS_U
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &        PP_LastPIVRPTRFilled_U
      INTEGER ISHIFT, K206
      INTEGER DMUMPS_IXAMAX
      INCLUDE 'mumps_headers.h'
      INTRINSIC max
      DOUBLE PRECISION, PARAMETER :: RZERO = 0.0D0
#if defined(_OPENMP)
      INTEGER :: NOMP, CHUNK, K360
      K360 = KEEP(360)
      NOMP    = OMP_GET_MAX_THREADS()
#endif
        NFRONT8 = int(NFRONT,8)
        INOPV   = 0
        XSIZE   = KEEP(IXSZ)
        NPIV    = IW(IOLDPS+1+XSIZE)
        NPIVP1  = NPIV + 1
        K206    = KEEP(206)
        IF ((KEEP(50).NE.1).AND.OOC_EFFECTIVE_ON_FRONT) THEN
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR_L, I_PIVR_L, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)
     &              +KEEP(IXSZ),
     &       IW, LIW)
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_U, NBPANELS_U, 
     &       I_PIVRPTR_U, I_PIVR_U, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
        ENDIF
        ISHIFT = 0   
        IF (K206.GE.1) THEN
          IF (Inextpiv.GT.NPIVP1.AND.Inextpiv.LE.NASS) THEN
            ISHIFT = Inextpiv - NPIVP1
          ENDIF
          IF (ISHIFT.GT.0.AND.IS_MAXFROMN_AVAIL) THEN
            IPIV = NPIVP1
            APOS = POSELT + NFRONT8*int(NPIV,8) + int(IPIV-1,8)
            IDIAG = APOS + int(IPIV - NPIVP1,8)*NFRONT8
            IF (abs(A(IDIAG)) .GE. UU*MAXFROMN .AND.
     &          abs(A(IDIAG)) .GT. max(SEUIL,tiny(RMAX)) ) THEN
              ISHIFT = 0
            ENDIF
          ENDIF
          IF ( ISHIFT .GT. 0) THEN
            IS_MAXFROMN_AVAIL = .FALSE.
          ENDIF
        ENDIF  
          DO 460 IPIV_SHIFT=NPIVP1+ISHIFT,NASS+ISHIFT
            IF (IPIV_SHIFT .LE. NASS) THEN
              IPIV=IPIV_SHIFT
            ELSE
              IPIV=IPIV_SHIFT-NASS-1+NPIVP1
            ENDIF
            APOS = POSELT + NFRONT8*int(NPIV,8) + int(IPIV-1,8)
            JMAX = 1
            AMROW = RZERO
            J1    = APOS
            J3    = NASS -NPIV
            JMAX  = DMUMPS_IXAMAX(J3,A(J1),NFRONT,KEEP(360))
            JJ    = J1 + int(JMAX-1,8)*NFRONT8
            AMROW = abs(A(JJ))
            RMAX  = AMROW
            J1    = APOS +  int(NASS-NPIV,8) * NFRONT8
            J3 = NFRONT - NASS - KEEP(253)-NVSCHUR
            IF (IS_MAXFROMN_AVAIL) THEN
              RMAX = max(MAXFROMN,RMAX)
              IS_MAXFROMN_AVAIL = .FALSE.
            ELSE
              IF (J3.EQ.0) GOTO 370
              IF (KEEP(351).EQ.1) THEN
                J1_ini = J1
!$              CHUNK = max(K360/2,(J3+NOMP-1)/NOMP)
!$OMP  PARALLEL DO schedule(static, CHUNK)
!$OMP& FIRSTPRIVATE(J1_ini,NFRONT8,J3)
!$OMP& REDUCTION(max:RMAX) IF (J3.GE.K360)
                DO J=1,J3
                  RMAX = max(abs(A(J1_ini + int(J-1,8) * NFRONT8)),
     &                       RMAX)
                END DO
!$OMP  END PARALLEL DO
              ELSE
                DO J=1,J3
                  RMAX = max(abs(A(J1)), RMAX)
                  J1 = J1 + NFRONT8
                END DO
              ENDIF
            END IF
  370       IF (RMAX.LE.tiny(RMAX)) GO TO 460
            IDIAG = APOS + int(IPIV - NPIVP1,8)*NFRONT8
            IF (abs(A(IDIAG)) .GE. UU*RMAX .AND.
     &          abs(A(IDIAG)) .GT. max(SEUIL,tiny(RMAX))) THEN
               JMAX = IPIV - NPIV
               GO TO 380
            ENDIF
            IF ( .NOT. (AMROW .GE. UU*RMAX .AND.
     &                  AMROW .GT. max(SEUIL,tiny(RMAX))) ) GO TO 460
            NOFFW = NOFFW + 1
  380       CONTINUE
            IF (K206.GE.1) THEN
              Inextpiv = IPIV + 1 
            ENDIF
           CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( abs(A(APOS + int(JMAX - 1,8) * NFRONT8 )), 
     &             DKEEP, KEEP, .FALSE.)
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER(
     &             A(APOS + int(JMAX - 1,8) * NFRONT8 ),
     &             DET_MANTW, DET_EXPW )
            ENDIF
            IF (IPIV.EQ.NPIVP1) GO TO 400
            IF (KEEP(405) .EQ.0) THEN
              KEEP8(80) = KEEP8(80)+1
            ELSE
!$OMP         ATOMIC UPDATE
              KEEP8(80) = KEEP8(80)+1
!$OMP         END ATOMIC
            ENDIF
            DET_SIGNW = - DET_SIGNW
            J1   = POSELT + int(NPIV,8)
            J3_8 = POSELT + int(IPIV-1,8)
            DO J= 1,NFRONT
              SWOP  = A(J1)
              A(J1) = A(J3_8)
              A(J3_8) = SWOP
              J1 = J1 + NFRONT8
              J3_8 = J3_8 + NFRONT8
            END DO
            ISWPS1 = IOLDPS + 5 + NPIVP1 + NFRONT + XSIZE
            ISWPS2 = IOLDPS + 5 + IPIV + NFRONT + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            DET_SIGNW = -DET_SIGNW
            J1 = POSELT + int(NPIV,8) * NFRONT8
            J2 = POSELT + int(NPIV + JMAX - 1,8) * NFRONT8
            DO KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + 1_8
              J2 = J2 + 1_8
            END DO
            ISWPS1 = IOLDPS + 5 + NPIV + 1 + XSIZE
            ISWPS2 = IOLDPS + 5 + NPIV + JMAX + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460    CONTINUE
       INOPV = 1
       GOTO 430
  420 CONTINUE
              IF (OOC_EFFECTIVE_ON_FRONT) THEN
                IF (KEEP(251).EQ.0) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, NPIV+JMAX,
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
                ENDIF
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, IPIV,
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
              ENDIF
 430  CONTINUE
      IS_MAXFROMN_AVAIL = .FALSE.
      RETURN
      END SUBROUTINE DMUMPS_FAC_H
      SUBROUTINE DMUMPS_FAC_M(IBEG_BLOCK,
     &     NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &     IOLDPS,POSELT,IFINB,LKJIB,LKJIT,XSIZE)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LIW,INODE,IFINB,LKJIB,IBEG_BLOCK
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    VALPIV
      INTEGER(8) :: APOS, POSELT, UUPOS, LPOS
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS
      INTEGER LKJIT, XSIZE
      DOUBLE PRECISION ONE, ALPHA
      INTEGER NPIV,JROW2
      INTEGER NEL2,NPIVP1,KROW,NEL
      INCLUDE 'mumps_headers.h'
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NFRONT8= int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        IFINB  = 0
        IF (IW(IOLDPS+3+XSIZE).LE.0) THEN
          IF (NASS.LT.LKJIT) THEN
           IW(IOLDPS+3+XSIZE) = NASS
          ELSE
           IW(IOLDPS+3+XSIZE) = min0(NASS,LKJIB)
          ENDIF
        ENDIF
        JROW2 = IW(IOLDPS+3+XSIZE)
        NEL2   = JROW2 - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
          IW(IOLDPS+3+XSIZE) = min0(JROW2+LKJIB,NASS)
          IBEG_BLOCK = NPIVP1+1
         ENDIF
        ELSE
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT8
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT8
 541     CONTINUE
         LPOS   = APOS + NFRONT8
         UUPOS  = APOS + 1_8
         CALL dger(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     &              A(LPOS+1_8),NFRONT)
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_FAC_M
      SUBROUTINE DMUMPS_FAC_N(NFRONT,NASS,IW,LIW,A,LA,
     &       IOLDPS,POSELT,IFINB,XSIZE,
     &       KEEP,MAXFROMN,IS_MAXFROMN_AVAIL,NVSCHUR)
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER NFRONT,NASS,LIW,IFINB
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    ALPHA,VALPIV
      INTEGER(8) :: APOS, POSELT, UUPOS, LPOS, IRWPOS
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS,NPIV,XSIZE
      INTEGER, intent(in) :: KEEP(500)
      DOUBLE PRECISION, intent(inout) :: MAXFROMN
      LOGICAL, intent(inout) :: IS_MAXFROMN_AVAIL
      INTEGER, intent(in)    :: NVSCHUR
      INTEGER NEL,IROW,NEL2,JCOL,NELMAXM
      INTEGER NPIVP1
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
#if defined(_OPENMP)
      LOGICAL:: OMP_FLAG
      INTEGER:: NOMP, K360, CHUNK
      NOMP = OMP_GET_MAX_THREADS()
      K360 = KEEP(360)
#endif
        NFRONT8=int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1 
        NELMAXM= NEL -KEEP(253)-NVSCHUR
        NEL2   = NASS - NPIVP1
        IFINB  = 0
        IF (NPIVP1.EQ.NASS) IFINB = 1
        APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
        VALPIV = ONE/A(APOS)
#if defined(_OPENMP)
        OMP_FLAG = .FALSE.
        CHUNK = max(NEL,1)
        IF (NOMP.GT.1) THEN
          IF (NEL.LT.K360) THEN
              IF (NEL*NEL2.GE.KEEP(361)) THEN
                OMP_FLAG = .TRUE.
                CHUNK = max(20, (NEL+NOMP-1)/NOMP)
              ENDIF
          ELSE
            OMP_FLAG = .TRUE.
            CHUNK = max(K360/2, (NEL+NOMP-1)/NOMP)
          ENDIF
        ENDIF
#endif
        IF (KEEP(351).EQ.2) THEN
          MAXFROMN = 0.0D0
          IF (NEL2 > 0) THEN
            IS_MAXFROMN_AVAIL = .TRUE.
          ENDIF
!$OMP PARALLEL DO schedule(static, CHUNK)
!$OMP& PRIVATE(LPOS, UUPOS, IRWPOS, ALPHA, JCOL)
!$OMP& FIRSTPRIVATE(APOS,NFRONT8,VALPIV,NEL,NEL2)
!$OMP& REDUCTION(max:MAXFROMN) IF(OMP_FLAG)
          DO IROW = 1, NEL
            LPOS = APOS + NFRONT8*int(IROW,8)
            A(LPOS) = A(LPOS)*VALPIV
            ALPHA   = -A(LPOS)
            IRWPOS  = LPOS + 1_8
            UUPOS  = APOS + 1_8
            IF (NEL2 > 0) THEN
              A(IRWPOS) = A(IRWPOS) + ALPHA*A(UUPOS)
              IF (IROW.LE.NELMAXM)
     &          MAXFROMN=max(MAXFROMN, abs(A(IRWPOS)))
              IRWPOS = IRWPOS+1_8
              UUPOS  = UUPOS+1_8
              DO JCOL = 2, NEL2
                A(IRWPOS) = A(IRWPOS) + ALPHA*A(UUPOS)
                IRWPOS = IRWPOS+1_8
                UUPOS  = UUPOS+1_8
              ENDDO
            ENDIF
          END DO
!$OMP END PARALLEL DO
        ELSE
!$OMP PARALLEL DO schedule(static, CHUNK)
!$OMP& FIRSTPRIVATE(APOS,NFRONT8,VALPIV,NEL,NEL2)
!$OMP& PRIVATE(LPOS, UUPOS, IRWPOS, ALPHA, JCOL) IF(OMP_FLAG)
          DO IROW = 1, NEL
            LPOS = APOS + NFRONT8*int(IROW,8)
            A(LPOS) = A(LPOS)*VALPIV
            ALPHA   = -A(LPOS)
            IRWPOS  = LPOS + 1_8
            UUPOS  = APOS + 1_8
            DO JCOL = 1, NEL2
              A(IRWPOS) = A(IRWPOS) + ALPHA*A(UUPOS)
              IRWPOS = IRWPOS+1_8
              UUPOS  = UUPOS+1_8
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FAC_N
      SUBROUTINE DMUMPS_FAC_PT_SETLOCK427( K427_OUT, K427,
     &           K405, K222, NEL1, NASS )
      INTEGER, INTENT(IN) :: K427, K405, K222, NEL1, NASS
      INTEGER, INTENT(OUT) :: K427_OUT
      K427_OUT = K427
      IF ( K405 .EQ. 1 ) THEN
        IF ( K427_OUT .GT. 0 ) K427_OUT = 0
        IF ( K427_OUT .LT. 0 ) K427_OUT = -1
      ENDIF
      IF ( K427_OUT .GT. 99 ) K427_OUT = 0
      IF ( K427_OUT .LT. -100 ) K427_OUT = -1
      RETURN
      END SUBROUTINE DMUMPS_FAC_PT_SETLOCK427
      SUBROUTINE DMUMPS_FAC_P(A,LA,NFRONT,
     &      NPIV,NASS,POSELT,CALL_UTRSM, KEEP, INODE,
     &      CALL_OOC, IWFAC, LIWFAC, LAFAC, MonBloc, MYID, KEEP8,
     &      LNextPiv2beWritten, UNextPiv2beWritten,
     &      IFLAG )
      USE DMUMPS_OOC, ONLY : IO_BLOCK, TYPEF_BOTH_LU,
     &                       DMUMPS_OOC_IO_LU_PANEL
      USE MUMPS_OOC_COMMON, ONLY : STRAT_TRY_WRITE
      IMPLICIT NONE
      INTEGER(8) :: LA,POSELT,LAFAC
      DOUBLE PRECISION    A(LA)
      INTEGER NFRONT, NPIV, NASS
      LOGICAL, INTENT(IN) :: CALL_UTRSM
      INTEGER, INTENT(INOUT) :: IFLAG
      LOGICAL, INTENT(IN) :: CALL_OOC
      INTEGER  LIWFAC, MYID,
     &      LNextPiv2beWritten, UNextPiv2beWritten
      INTEGER  IWFAC(LIWFAC)
      TYPE(IO_BLOCK) :: MonBloc
      INTEGER :: KEEP(500)
      INTEGER(8) :: KEEP8(150)
      INTEGER(8) :: LPOS, LPOS1, LPOS2, UPOS
      INTEGER NEL1, NEL11, IFLAG_OOC
      INTEGER :: INODE
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INCLUDE 'mumps_headers.h'
      NEL1   = NFRONT - NASS
      NEL11  = NFRONT - NPIV
      LPOS2  = POSELT + int(NASS,8)*int(NFRONT,8) 
      LPOS   = LPOS2 + int(NPIV,8)  
      LPOS1  = POSELT + int(NPIV,8) 
      UPOS   = POSELT + int(NASS,8) 
        IF ( CALL_UTRSM ) THEN
          CALL dtrsm('R', 'U', 'N', 'U', NEL1, NPIV, ONE,
     &                A(POSELT), NFRONT, A(UPOS), NFRONT)
        ENDIF
      CALL dtrsm('L','L','N','N',NPIV,NEL1,ONE,A(POSELT),NFRONT,
     &            A(LPOS2),NFRONT)
      IF (CALL_OOC) THEN
           CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT_TRY_WRITE, TYPEF_BOTH_LU, 
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IWFAC, LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,
     &           .FALSE. ) 
           IF (IFLAG_OOC .LT. 0) THEN
             IFLAG = IFLAG_OOC
             GOTO 500
           ENDIF
      ENDIF
      CALL dgemm('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     &            NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
      IF ((CALL_UTRSM).AND.(NASS-NPIV.GT.0)) THEN
        LPOS2 = POSELT + int(NPIV,8)*int(NFRONT,8) 
        LPOS  = LPOS2  + int(NASS,8) 
        CALL dgemm('N','N',NEL1,NASS-NPIV,NPIV,ALPHA,A(UPOS),
     &            NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
      ENDIF
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_P
      SUBROUTINE DMUMPS_FAC_T(A,LA,NPIVB,NFRONT,
     &                             NPIV,NASS,POSELT)
      IMPLICIT NONE
      INTEGER NPIVB,NASS
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER(8) :: APOS, POSELT
      INTEGER NFRONT, NPIV, NASSL
      INTEGER(8) :: LPOS, LPOS1, LPOS2
      INTEGER NEL1, NEL11, NPIVE
      DOUBLE PRECISION    ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        NPIVE  = NPIV - NPIVB
        NASSL  = NASS - NPIVB
        APOS   = POSELT + int(NPIVB,8)*int(NFRONT,8)
     &                  + int(NPIVB,8)
        LPOS2  = APOS + int(NASSL,8)
        CALL dtrsm('R','U','N','U',NEL1,NPIVE,ONE,A(APOS),NFRONT,
     &              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + int(NFRONT,8)*int(NPIVE,8)
        LPOS1  = APOS  + int(NFRONT,8)*int(NPIVE,8)
        CALL dgemm('N','N',NEL1,NEL11,NPIVE,ALPHA,A(LPOS2),
     &          NFRONT,A(LPOS1),NFRONT,ONE,A(LPOS),NFRONT)
        RETURN
        END SUBROUTINE DMUMPS_FAC_T
      SUBROUTINE DMUMPS_FAC_SQ(IBEG_BLOCK, IEND_BLOCK, NPIV,
     &    NFRONT, LAST_ROW, LAST_COL, A, LA, POSELT, 
     &    FIRST_COL, CALL_LTRSM, CALL_UTRSM, CALL_GEMM, 
     &    WITH_COMM_THREAD, LR_ACTIVATED
     &    )
!$    USE OMP_LIB
#if defined(_OPENMP)
      USE DMUMPS_BUF
#endif
      IMPLICIT NONE
      INTEGER, intent(in)     :: IBEG_BLOCK, IEND_BLOCK
      INTEGER, intent(in)     :: NPIV, NFRONT, LAST_ROW, LAST_COL
      INTEGER, intent(in)     :: FIRST_COL
      INTEGER(8), intent(in)  :: LA
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      INTEGER(8), intent(in)  :: POSELT 
      LOGICAL, intent(in)     :: CALL_LTRSM, CALL_UTRSM, CALL_GEMM
      LOGICAL, intent(in) :: WITH_COMM_THREAD, LR_ACTIVATED
      INTEGER(8) :: NFRONT8, LPOSN, LPOS2N
      INTEGER(8) :: LPOS, LPOS1, LPOS2, UPOS, POSELT_LOCAL
      INTEGER :: NELIM, LKJIW, NEL1, NEL11, UTRSM_NCOLS
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
!$    INTEGER :: NOMP
!$    LOGICAL :: TRSM_GEMM_FINISHED
!$    LOGICAL :: SAVE_NESTED, SAVE_DYNAMIC
      NFRONT8= int(NFRONT,8)
      NELIM  = IEND_BLOCK - NPIV
      NEL1   = LAST_ROW - IEND_BLOCK
      IF ( NEL1 < 0 ) THEN
        WRITE(*,*)
     &  "Internal error 1 in DMUMPS_FAC_SQ,IEND_BLOCK>LAST_ROW",
     &  IEND_BLOCK, LAST_ROW
        CALL MUMPS_ABORT()
      ENDIF
      LKJIW  = NPIV - IBEG_BLOCK + 1
      NEL11  = LAST_COL - NPIV
      LPOS2  = POSELT + int(IEND_BLOCK,8)*NFRONT8 + int(IBEG_BLOCK-1,8)
      UTRSM_NCOLS = LAST_COL - FIRST_COL
      UPOS   = POSELT + int(IBEG_BLOCK-1,8)*NFRONT8 + int(FIRST_COL,8)
      POSELT_LOCAL = POSELT + int(IBEG_BLOCK-1,8)*NFRONT8 
     &                      + int(IBEG_BLOCK-1,8)
      IF ((NEL1.NE.0).AND.(LKJIW.NE.0)) THEN
        IF (WITH_COMM_THREAD .EQV. .FALSE.) THEN
           IF (CALL_LTRSM) THEN
             CALL dtrsm('L','L','N','N',LKJIW,NEL1,ONE,
     &                A(POSELT_LOCAL),NFRONT,
     &                A(LPOS2),NFRONT)
           ENDIF
           IF (CALL_UTRSM) THEN
             CALL dtrsm('R','U','N','U',UTRSM_NCOLS,LKJIW,ONE,
     &                A(POSELT_LOCAL),NFRONT,
     &                A(UPOS),NFRONT)
             LPOS2N = POSELT + int(NPIV,8)*NFRONT8 + int(IBEG_BLOCK-1,8)
             LPOSN  = POSELT + int(NPIV,8)*NFRONT8 + int(FIRST_COL,8)
             CALL dgemm('N','N',UTRSM_NCOLS,NELIM,
     &               LKJIW,ALPHA,A(UPOS),NFRONT,A(LPOS2N),
     &               NFRONT,ONE,A(LPOSN),NFRONT)
           ENDIF
           IF (CALL_GEMM) THEN
            LPOS   = LPOS2 + int(LKJIW,8)
            LPOS1  = POSELT_LOCAL + int(LKJIW,8)
            CALL dgemm('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     &           NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
           ENDIF
        ELSE
!$        NOMP = OMP_GET_MAX_THREADS()
!$        CALL OMP_SET_NUM_THREADS(2)
!$        SAVE_NESTED = OMP_GET_NESTED()
!$        SAVE_DYNAMIC = OMP_GET_DYNAMIC()
!$        CALL OMP_SET_NESTED(.TRUE.)
!$        CALL OMP_SET_DYNAMIC(.FALSE.)
!$        TRSM_GEMM_FINISHED = .FALSE.
!$OMP     PARALLEL SHARED(TRSM_GEMM_FINISHED)
!$        IF (OMP_GET_THREAD_NUM() .EQ. 1) THEN
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$          CALL OMP_SET_NUM_THREADS(int(NOMP,4))
#else
!$          CALL OMP_SET_NUM_THREADS(NOMP)
#endif
          IF (CALL_LTRSM) THEN
            CALL dtrsm('L','L','N','N',LKJIW,NEL1,ONE,
     &                 A(POSELT_LOCAL),NFRONT,
     &                 A(LPOS2),NFRONT)
          ENDIF
          IF (CALL_UTRSM) THEN
            CALL dtrsm('R','U','N','U',UTRSM_NCOLS,LKJIW,ONE,
     &                 A(POSELT_LOCAL),NFRONT,
     &                 A(UPOS),NFRONT)
            LPOS2N = POSELT + int(NPIV,8)*NFRONT8 + int(IBEG_BLOCK-1,8)
            LPOSN  = POSELT + int(NPIV,8)*NFRONT8 + int(FIRST_COL,8)
            CALL dgemm('N','N',UTRSM_NCOLS,NELIM,
     &                LKJIW,ALPHA,A(UPOS),NFRONT,A(LPOS2N),
     &                NFRONT,ONE,A(LPOSN),NFRONT)
          ENDIF
          IF (CALL_GEMM) THEN
            LPOS   = LPOS2 + int(LKJIW,8)
            LPOS1  = POSELT_LOCAL + int(LKJIW,8)
            CALL dgemm('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     &             NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
          END IF
!$           TRSM_GEMM_FINISHED = .TRUE.
!$        ELSE
!$           DO WHILE (.NOT. TRSM_GEMM_FINISHED)
!$             CALL DMUMPS_BUF_TEST()
!$             CALL MUMPS_USLEEP(10000)
!$           END DO
!$        END IF
!$OMP     END PARALLEL
!$        CALL OMP_SET_NESTED(SAVE_NESTED)
!$        CALL OMP_SET_DYNAMIC(SAVE_DYNAMIC)
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$        CALL OMP_SET_NUM_THREADS(int(NOMP,4))
#else
!$        CALL OMP_SET_NUM_THREADS(NOMP)
#endif
        ENDIF
      ELSE
        IF (CALL_UTRSM.AND.UTRSM_NCOLS.NE.0) THEN
          CALL dtrsm('R','U','N','U',UTRSM_NCOLS,LKJIW,ONE,
     &                 A(POSELT_LOCAL),NFRONT,
     &                 A(UPOS),NFRONT)
          LPOS2N = POSELT + int(NPIV,8)*NFRONT8 + int(IBEG_BLOCK-1,8)
          LPOSN  = POSELT + int(NPIV,8)*NFRONT8 + int(FIRST_COL,8)
          CALL dgemm('N','N',UTRSM_NCOLS,NELIM,
     &              LKJIW,ALPHA,A(UPOS),NFRONT,A(LPOS2N),
     &              NFRONT,ONE,A(LPOSN),NFRONT)
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FAC_SQ
      SUBROUTINE DMUMPS_FAC_MQ(IBEG_BLOCK,IEND_BLOCK,
     &     NFRONT, NASS, NPIV, LAST_COL, A, LA, POSELT, IFINB,
     &     LR_ACTIVATED
     &     )
      IMPLICIT NONE
      INTEGER, intent(in)    :: IBEG_BLOCK, IEND_BLOCK, NFRONT, 
     &                          NASS, NPIV, LAST_COL
      INTEGER, intent(out)   ::  IFINB
      INTEGER(8), intent(in) :: LA, POSELT
      DOUBLE PRECISION, intent(inout) :: A(LA)
      LOGICAL, intent(in)    :: LR_ACTIVATED
      DOUBLE PRECISION    :: VALPIV
      INTEGER(8) :: APOS,  UUPOS, LPOS
      INTEGER(8) :: NFRONT8
      DOUBLE PRECISION    :: ONE, ALPHA
      INTEGER    :: NEL2,NPIVP1,KROW,NEL
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NFRONT8= int(NFRONT,8)
        NPIVP1 = NPIV + 1
        NEL  = LAST_COL - NPIVP1
        IFINB  = 0
        NEL2   = IEND_BLOCK - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (IEND_BLOCK.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
         ENDIF
        ELSE
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT8
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT8
 541     CONTINUE
         LPOS   = APOS + NFRONT8
         UUPOS  = APOS + 1_8
#if defined(MUMPS_USE_BLAS2)
         CALL dger(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     &              A(LPOS+1_8),NFRONT)
#else
         CALL dgemm('N','N',NEL,NEL2,1,ALPHA,A(UUPOS),NEL,
     &               A(LPOS),NFRONT,ONE,A(LPOS+1_8),NFRONT)
#endif
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_FAC_MQ
      SUBROUTINE DMUMPS_FAC_FR_UPDATE_CBROWS( INODE, NFRONT, NASS, 
     &     CALL_UTRSM, A, LA, LAFAC, POSELT, IW, LIW, IOLDPS,
     &     MonBloc, MYID, NOFFW,
     &     DET_EXPW, DET_MANTW, DET_SIGNW,
     &     LIWFAC,
     &     PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &     LNextPiv2beWritten, UNextPiv2beWritten, 
     &     PP_LastPIVRPTRFilled_L, PP_LastPIVRPTRFilled_U,
     &     
     &     XSIZE, SEUIL, UU, DKEEP, KEEP8, KEEP, IFLAG, 
     &     OOC_EFFECTIVE_ON_FRONT, NVSCHUR)
      USE DMUMPS_OOC, ONLY: IO_BLOCK
      IMPLICIT NONE
      INTEGER, intent(in)    :: INODE, NFRONT, NASS,
     &                          LIW, MYID, XSIZE, IOLDPS, LIWFAC
      INTEGER(8), intent(in) :: LA, POSELT
      INTEGER, intent(inout) :: NOFFW
      INTEGER, intent(inout) :: DET_EXPW, DET_SIGNW
      DOUBLE PRECISION, intent(inout) :: DET_MANTW
      INTEGER, intent(inout) :: PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &                   LNextPiv2beWritten, UNextPiv2beWritten,
     &                   PP_LastPIVRPTRFilled_L, PP_LastPIVRPTRFilled_U,
     &                   IFLAG
      LOGICAL, intent(in)    :: CALL_UTRSM
      INTEGER, intent(inout) :: IW(LIW)
      DOUBLE PRECISION, intent(inout) :: A(LA)
      DOUBLE PRECISION, intent(in)       :: SEUIL, UU, DKEEP(230)
      INTEGER, intent(in)    :: KEEP( 500 ) 
      INTEGER(8), intent(inout) :: LAFAC
      INTEGER(8)             :: KEEP8(150)
      INTEGER, intent(in)    :: NVSCHUR
      TYPE(IO_BLOCK), intent(inout)  :: MonBloc
      LOGICAL, intent(in)    :: OOC_EFFECTIVE_ON_FRONT
      INTEGER  :: NPIV, NEL1, IBEG_BLOCK, IFINB, INOPV
      INTEGER Inextpiv
      DOUBLE PRECISION :: MAXFROMN
      LOGICAL :: IS_MAXFROMN_AVAIL
      NPIV   = IW(IOLDPS+1+XSIZE)
      NEL1   = NFRONT - NASS
      IF (KEEP(206).GE.1) THEN
        Inextpiv = 1   
      ELSE 
        Inextpiv = 0   
      ENDIF
      IF ((NPIV.GT.0).AND.(NEL1.GT.0)) THEN
        IF (OOC_EFFECTIVE_ON_FRONT) THEN
          MonBloc%LastPiv = NPIV
        ENDIF
        CALL DMUMPS_FAC_P(A,LA,NFRONT, NPIV, NASS, POSELT, 
     &                    CALL_UTRSM, KEEP, INODE,
     &                    OOC_EFFECTIVE_ON_FRONT, IW(IOLDPS), 
     &                    LIWFAC, LAFAC,
     &                    MonBloc, MYID, KEEP8,
     &                    LNextPiv2beWritten, UNextPiv2beWritten,
     &                    IFLAG)
      ENDIF
        NPIV   = IW(IOLDPS+1+XSIZE)
        IBEG_BLOCK = NPIV
        IF (NASS.EQ.NPIV) GOTO 500
        IS_MAXFROMN_AVAIL = .FALSE.
 120    CALL DMUMPS_FAC_H(NFRONT,NASS,IW,LIW,A,LA,
     &     INOPV, NOFFW,
     &     DET_EXPW, DET_MANTW, DET_SIGNW,
     &     IOLDPS,POSELT,UU,SEUIL,
     &     KEEP, KEEP8, DKEEP,
     &     PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &     PP_LastPIVRPTRFilled_U, MAXFROMN, IS_MAXFROMN_AVAIL,
     &     Inextpiv, OOC_EFFECTIVE_ON_FRONT, NVSCHUR
     &     )
        IF (INOPV.NE.1) THEN
         CALL DMUMPS_FAC_N(NFRONT,NASS,IW,LIW,A,LA,
     &                 IOLDPS,POSELT,IFINB,XSIZE,
     &                 KEEP, MAXFROMN, IS_MAXFROMN_AVAIL,
     &                 NVSCHUR)
         IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
         IF (IFINB.EQ.0) GOTO 120
        ENDIF
        NPIV   = IW(IOLDPS+1+XSIZE)
        NEL1   = NFRONT - NASS
        IF ((NPIV.LE.IBEG_BLOCK).OR.(NEL1.EQ.0)) GO TO 500
        CALL DMUMPS_FAC_T(A,LA,IBEG_BLOCK,
     &                NFRONT,NPIV,NASS,POSELT)
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_FR_UPDATE_CBROWS
      SUBROUTINE DMUMPS_FAC_I(NFRONT,NASS,LAST_ROW,
     &    IBEG_BLOCK, IEND_BLOCK,
     &    N,INODE,IW,LIW,A,LA,
     &    INOPV,NOFFW,NBTINYW,
     &    DET_EXPW, DET_MANTW, DET_SIGNW,
     &    IFLAG,IOLDPS,POSELT,UU,SEUIL,KEEP,KEEP8,
     &    DKEEP,PIVNUL_LIST,LPN_LIST,
     &
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &     PP_LastPIVRPTRFilled_U,
     &     PIVOT_OPTION, LR_ACTIVATED, IEND_BLR, Inextpiv,
     &     OOC_EFFECTIVE_ON_FRONT, NVSCHUR, PARPIV_T1,
     &     TIPIV    
     &     )
!$    USE OMP_LIB
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER, intent(in)    :: IBEG_BLOCK, IEND_BLOCK
      INTEGER, intent(inout), OPTIONAL :: TIPIV(:)
      INTEGER(8), intent(in) :: LA
      DOUBLE PRECISION, intent(inout) :: A(LA)
      INTEGER, intent(in)    :: NFRONT,NASS,N,LIW,INODE,LAST_ROW
      INTEGER, intent(inout) :: IFLAG,INOPV,NOFFW, NBTINYW
      INTEGER, intent(inout) :: DET_EXPW, DET_SIGNW
      DOUBLE PRECISION, intent(inout) :: DET_MANTW
      DOUBLE PRECISION, intent(in)       :: UU, SEUIL
      INTEGER, intent(inout) :: IW(LIW)
      INTEGER, intent(in)    :: IOLDPS
      INTEGER(8), intent(in) :: POSELT
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in)    :: LPN_LIST
      INTEGER, intent(inout) :: PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION    DKEEP(230)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &        PP_LastPIVRPTRFilled_U
      INTEGER, intent(in)    :: PIVOT_OPTION, IEND_BLR
      LOGICAL, intent(in)    :: LR_ACTIVATED
      INTEGER, intent(inout) :: Inextpiv
      LOGICAL, intent(in)    :: OOC_EFFECTIVE_ON_FRONT
      INTEGER, intent(in)    :: NVSCHUR
      INTEGER, intent(in) :: PARPIV_T1
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION SWOP
      INTEGER XSIZE
      INTEGER(8) :: APOS, IDIAG
      INTEGER(8) :: J1, J2, JJ, J3
      INTEGER(8) :: NFRONT8
      INTEGER ILOC
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      DOUBLE PRECISION RZERO, RMAX, AMROW, MAX_PREV_in_PARPIV
      INTEGER(8) :: APOSMAX, APOSROW
      DOUBLE PRECISION       :: RMAX_NORELAX
      DOUBLE PRECISION PIVNUL
      DOUBLE PRECISION FIXA, CSEUIL
      INTEGER NPIV,IPIV, LRLOC
      INTEGER NPIVP1,JMAX,J,ISW,ISWPS1
      INTEGER ISWPS2,KSW, HF, IPIVNUL
      INTEGER DMUMPS_IXAMAX
      INTEGER :: ISHIFT, K206
      INTEGER :: IPIV_SHIFT,IPIV_END
      INTRINSIC max
      DATA RZERO /0.0D0/
#if defined(_OPENMP)
      INTEGER :: NOMP,CHUNK,K361
#endif
      INTEGER I_PIVRPTR_L, I_PIVR_L, NBPANELS_L
      INTEGER I_PIVRPTR_U, I_PIVR_U, NBPANELS_U
#if defined(_OPENMP)
        NOMP    = OMP_GET_MAX_THREADS()
        K361    = KEEP(361)
#endif
        PIVNUL  = DKEEP(1)
        FIXA    = DKEEP(2)
        CSEUIL  = SEUIL
        NFRONT8 = int(NFRONT,8)
        K206    = KEEP(206)
        XSIZE   = KEEP(IXSZ)
        NPIV    = IW(IOLDPS+1+XSIZE)
        HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
        NPIVP1  = NPIV + 1
        APOSMAX = POSELT+NFRONT8*NFRONT8-1_8
        IF (OOC_EFFECTIVE_ON_FRONT) THEN
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR_L, I_PIVR_L, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_U, NBPANELS_U, 
     &       I_PIVRPTR_U, I_PIVR_U, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
        ENDIF
        IF ( present(TIPIV) ) THEN
          ILOC    = NPIVP1 - IBEG_BLOCK + 1
          TIPIV(ILOC) = ILOC
        ENDIF
        IF (INOPV .EQ. -1) THEN
           APOS = POSELT + NFRONT8*int(NPIVP1-1,8) + int(NPIV,8)
           IDIAG = APOS
           CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( abs(A(APOS)), DKEEP, KEEP, .TRUE.)
           IF(abs(A(APOS)).LT.SEUIL) THEN
              IF (dble(A(APOS)) .GE. RZERO) THEN
                 A(APOS) = CSEUIL
              ELSE
                 A(APOS) = -CSEUIL
              ENDIF
              NBTINYW = NBTINYW + 1
           ELSE IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER(A(APOS), DET_MANTW, DET_EXPW )
           ENDIF
           IF (OOC_EFFECTIVE_ON_FRONT) THEN
             IF (KEEP(251).EQ.0) THEN 
               CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
             ENDIF
             CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
      ISHIFT   = 0            
      IPIV_END = IEND_BLOCK   
      IF (K206.GE.1) THEN
        IF (Inextpiv.GT.NPIVP1.AND.Inextpiv.LE.IEND_BLOCK) THEN
          ISHIFT = Inextpiv - NPIVP1
        ENDIF
        IF ( K206.EQ.1
     &      .OR.  (K206 .GT.1 .AND. IEND_BLOCK.EQ.IEND_BLR) ) THEN
          IPIV_END = IEND_BLOCK + ISHIFT
        ENDIF
      ENDIF
          DO 460 IPIV_SHIFT = NPIVP1+ISHIFT, IPIV_END
            IF (IPIV_SHIFT .LE. IEND_BLOCK) THEN
              IPIV=IPIV_SHIFT
            ELSE
              IPIV = IPIV_SHIFT-IEND_BLOCK-1+NPIVP1
              IF (IBEG_BLOCK.EQ.NPIVP1) THEN
                EXIT
              ENDIF
            ENDIF
            APOS = POSELT + NFRONT8*int(IPIV-1,8) + int(NPIV,8)
            JMAX = 1
            IF ((PIVOT_OPTION.EQ.0).OR.(UU.EQ.RZERO)) THEN
              IF (A(APOS).EQ.ZERO) GO TO 630 
              GO TO 380 
            ENDIF
            AMROW = RZERO
            J1 = APOS
            IF (PIVOT_OPTION.EQ.1 .OR. (LR_ACTIVATED .AND.
     &         (KEEP(480).GE.2 
     &         ))) THEN
              J = IEND_BLR - NPIV
            ELSE
              J = NASS - NPIV
            ENDIF
            J2 = J1 + J - 1_8
            JMAX  = DMUMPS_IXAMAX(J,A(J1),1,KEEP(361))
            JJ    = J1 + int(JMAX - 1,8)
            AMROW = abs(A(JJ))
            RMAX = AMROW
            IF (PIVOT_OPTION.GE.2) THEN
              J1 = J2 + 1_8
              IF (PIVOT_OPTION.GE.3
     &           ) THEN
                J2 = APOS +
     &          int(- NPIV + NFRONT - 1 - KEEP(253) - NVSCHUR,8)
              ELSE
                J2 = APOS +int(- NPIV + NASS - 1 ,8)
              ENDIF
              IF (J2.LT.J1) GO TO 370
              IF (KEEP(351).EQ.1) THEN
!$              CHUNK = max(K361/2,(int(J2-J1)+NOMP-1)/NOMP)
!$OMP  PARALLEL DO schedule(static, CHUNK) PRIVATE(JJ)
!$OMP& FIRSTPRIVATE(J1,J2)
!$OMP& REDUCTION(max:RMAX) IF ((J2-J1).GE.K361)
                DO JJ=J1,J2
                  RMAX = max(abs(A(JJ)),RMAX)
                ENDDO
!$OMP  END PARALLEL DO
              ELSE
                DO 360 JJ=J1,J2
                  RMAX = max(abs(A(JJ)),RMAX)
  360           CONTINUE
              ENDIF
  370         CONTINUE
            ENDIF
            IDIAG = APOS + int(IPIV - NPIVP1,8)
            IF (PARPIV_T1.NE.0) THEN
             RMAX_NORELAX = dble(A(APOSMAX+int(IPIV,8)))
            ELSE
             RMAX_NORELAX = RZERO
            ENDIF
            RMAX         = max(RMAX,RMAX_NORELAX)
            IF ( RMAX .LE. PIVNUL ) THEN
               IF (LAST_ROW.EQ.NFRONT) THEN
                 LRLOC = LAST_ROW -KEEP(253)-NVSCHUR
               ELSE
                 LRLOC = LAST_ROW
               ENDIF
               IF (NFRONT - KEEP(253) .EQ. NASS) THEN 
                 IF (IEND_BLOCK.NE.NASS ) THEN 
                   GOTO 460
                 ENDIF
                 J1=POSELT+int(IPIV-1,8)+int(NPIV,8)*NFRONT8
                 J2=POSELT+int(IPIV-1,8)+int(LRLOC-1,8)*NFRONT8
               ELSE
                 J1=POSELT+int(IPIV-1,8)
                 J2=POSELT+int(IPIV-1,8)+int(LRLOC-1,8)*NFRONT8
               ENDIF
               DO JJ=J1, J2, NFRONT8
                 IF ( abs(A(JJ)) .GT. PIVNUL ) THEN
                   GOTO 460
                 END IF
               ENDDO
               IF ((PARPIV_T1.NE.0)
     &             .AND.(PARPIV_T1.NE.-1)
     &             .AND.(RMAX_NORELAX.LT.0)
     &             .AND.(IPIV.GT.1)) THEN
                 MAX_PREV_in_PARPIV = RZERO
                 DO JJ=1,IPIV-1
                  MAX_PREV_in_PARPIV= max ( MAX_PREV_in_PARPIV, 
     &                dble(A(APOSMAX+int(JJ,8))) )
                 ENDDO
                 IF (MAX_PREV_in_PARPIV.GT.PIVNUL) THEN
                   APOSROW = POSELT + NFRONT8*int(IPIV-1,8)
                   DO JJ=1,IPIV-1
                     IF (abs(A(APOSROW+JJ-1)).GT.PIVNUL) GOTO 460
                   ENDDO
                 ENDIF
               ENDIF
               CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( abs(A(IDIAG)), DKEEP, KEEP, .TRUE.)
!$OMP ATOMIC CAPTURE
               KEEP(109) = KEEP(109)+1
               IPIVNUL = KEEP(109)
!$OMP END ATOMIC
               PIVNUL_LIST(IPIVNUL) = IW( IOLDPS+HF+NPIV+IPIV-NPIVP1 )
               IF(dble(FIXA).GT.RZERO) THEN
                  IF(dble(A(IDIAG)) .GE. RZERO) THEN
                     A(IDIAG) = FIXA
                  ELSE
                     A(IDIAG) = -FIXA
                  ENDIF
               ELSE
                 J1 = APOS
                  J2 = APOS +
     &                int(- NPIV + NFRONT - 1 - KEEP(253),8)
                 DO JJ=J1,J2
                   A(JJ) = ZERO
                 ENDDO
                 A(IDIAG) = -FIXA
               ENDIF
               JMAX = IPIV - NPIV
               GOTO 385   
            ENDIF
            RMAX         = max(RMAX,abs(RMAX_NORELAX))
            IF (abs(A(IDIAG)) .GE. UU*RMAX .AND.
     &          abs(A(IDIAG)) .GT. max(SEUIL,tiny(RMAX))) THEN
               JMAX = IPIV - NPIV
               GO TO 380
            ENDIF
            IF ( .NOT. (AMROW .GE. UU*RMAX .AND.
     &                  AMROW .GT. max(SEUIL,tiny(RMAX))) ) GO TO 460
            NOFFW = NOFFW + 1
  380       CONTINUE
          IF (K206.GE.1) THEN
             Inextpiv = IPIV + 1
          ENDIF
            CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( abs(A(APOS+int(JMAX-1,8))),
     &             DKEEP, KEEP, .FALSE.)
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER( A(APOS+int(JMAX-1,8)),
     &                                 DET_MANTW,
     &                                 DET_EXPW )
            ENDIF
  385       CONTINUE
            IF (IPIV.EQ.NPIVP1) GO TO 400
            IF (KEEP(405) .EQ. 0) THEN
              KEEP8(80) = KEEP8(80)+1
            ELSE
!$OMP         ATOMIC UPDATE
              KEEP8(80) = KEEP8(80)+1
!$OMP         END ATOMIC
            ENDIF
            IF (PARPIV_T1.NE.0) THEN
             SWOP = A(APOSMAX+int(NPIVP1,8))
             A(APOSMAX+int(NPIVP1,8))   = A(APOSMAX+int(IPIV,8))
             A(APOSMAX+int(IPIV,8)) = SWOP
            ENDIF
            DET_SIGNW = - DET_SIGNW
            J1 = POSELT + int(NPIV,8)*NFRONT8
            J2 = J1 + NFRONT8 - 1_8
            J3 = POSELT + int(IPIV-1,8)*NFRONT8
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3)
              A(J3) = SWOP
              J3 = J3 + 1_8
  390       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NPIVP1
            ISWPS2 = IOLDPS + HF - 1 + IPIV
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            DET_SIGNW = - DET_SIGNW
            IF ( present(TIPIV) ) THEN
              TIPIV(ILOC) = ILOC + JMAX - 1
            ENDIF
            J1 = POSELT + int(NPIV,8)
            J2 = POSELT + int(NPIV + JMAX - 1,8)
            DO 410 KSW=1,LAST_ROW
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT8
              J2 = J2 + NFRONT8
  410       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NFRONT + NPIV + 1
            ISWPS2 = IOLDPS + HF - 1 + NFRONT + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420  
  460     CONTINUE
          IF (K206 .GE. 1) THEN
            Inextpiv=IEND_BLOCK+1
          ENDIF
      IF (IEND_BLOCK.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 430
  630 CONTINUE
      IFLAG = -10
      GOTO 430
  420 CONTINUE
              IF (OOC_EFFECTIVE_ON_FRONT) THEN
                IF (KEEP(251).EQ.0) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, IPIV, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
                ENDIF
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIV+JMAX, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
              ENDIF
  430 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_I
      SUBROUTINE DMUMPS_FAC_I_LDLT 
     &   ( NFRONT,NASS,INODE,IBEG_BLOCK,IEND_BLOCK,
     &     IW,LIW, A,LA, INOPV,
     &     NNEGW, NB22T1W, NBTINYW,
     &     DET_EXPW, DET_MANTW, DET_SIGNW,
     &     IFLAG,IOLDPS,POSELT,UU, SEUIL,KEEP,KEEP8,PIVSIZ,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, XSIZE,
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk,
     &     PP_LastPIVRPTRIndexFilled,MAXFROMM,IS_MAXFROMM_AVAIL, 
     &     PIVOT_OPTION, IEND_BLR, Inextpiv, 
     &     OOC_EFFECTIVE_ON_FRONT,
     &     NVSCHUR, PARPIV_T1, LR_ACTIVATED
     &     )
!$    USE OMP_LIB
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER(8) :: POSELT, LA
      INTEGER NFRONT,NASS,LIW,INODE,IFLAG,INOPV,
     &        IOLDPS
      INTEGER, intent(inout) :: NNEGW, NB22T1W, NBTINYW
      INTEGER, intent(inout) :: DET_EXPW, DET_SIGNW
      DOUBLE PRECISION, intent(inout) :: DET_MANTW
      INTEGER, intent(in) :: IBEG_BLOCK, IEND_BLOCK
      INTEGER, intent(in)    :: PIVOT_OPTION,IEND_BLR
      INTEGER, intent(inout) :: Inextpiv
      LOGICAL, intent(in)    :: OOC_EFFECTIVE_ON_FRONT
      INTEGER PIVSIZ,LPIV, XSIZE
      DOUBLE PRECISION A(LA) 
      DOUBLE PRECISION UU, UULOC, SEUIL
      INTEGER IW(LIW)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(230)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk
      INTEGER PP_LastPIVRPTRIndexFilled
      DOUBLE PRECISION, intent(in) :: MAXFROMM
      LOGICAL, intent(inout) :: IS_MAXFROMM_AVAIL
      INTEGER, intent(in)    :: NVSCHUR
      INTEGER, intent(in) :: PARPIV_T1
      LOGICAL, intent(in) :: LR_ACTIVATED
      include 'mpif.h'
      INTEGER (8) :: POSPV1,POSPV2,OFFDAG,APOSJ
      INTEGER JMAX, LIM, LIM_SWAP
      DOUBLE PRECISION RMAX,AMAX,TMAX, MAX_PREV_in_PARPIV
      DOUBLE PRECISION RMAX_NORELAX, TMAX_NORELAX, UULOCM1
      INTEGER(8) :: APOSMAX, APOSROW
      DOUBLE PRECISION MAXPIV
      DOUBLE PRECISION PIVNUL
      DOUBLE PRECISION FIXA, CSEUIL
      DOUBLE PRECISION PIVOT,DETPIV
      INCLUDE 'mumps_headers.h'
      INTEGER :: HF, IPIVNUL
      INTEGER :: J
      INTEGER(8) :: APOS, J1, J2, JJ, NFRONT8, KK, J1_ini, JJ_ini
      INTEGER    :: LDA
      INTEGER(8) :: LDA8
      INTEGER NPIV,IPIV
      INTEGER NPIVP1,K 
      INTEGER :: ISHIFT, K206, IPIV_SHIFT, IPIV_END
      INTRINSIC max
      DOUBLE PRECISION ZERO, ONE
      PARAMETER( ZERO = 0.0D0 )
      PARAMETER( ONE = 1.0D0 )
      DOUBLE PRECISION RZERO,RONE
      PARAMETER(RZERO=0.0D0, RONE=1.0D0)
#if defined(_OPENMP)
      LOGICAL :: OMP_FLAG
      INTEGER :: NOMP, CHUNK, J1_end
#endif
      INTEGER I_PIVRPTR, I_PIVR, NBPANELS_L
!$    NOMP = OMP_GET_MAX_THREADS()
      PIVNUL = DKEEP(1)
      FIXA   = DKEEP(2)
      CSEUIL = SEUIL
      LDA     = NFRONT
      LDA8    = int(LDA,8)
      NFRONT8   = int(NFRONT,8)
      K206      = KEEP(206)
      UULOC = UU
      IF (UULOC.GT.RZERO) THEN 
        UULOCM1 = RONE/UULOC
      ELSE
        UULOCM1  = RONE
      ENDIF
      HF = 6 + XSIZE
      IF (KEEP(50).NE.1 .AND. OOC_EFFECTIVE_ON_FRONT) THEN
             CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR, I_PIVR, IOLDPS+2*NFRONT+6+KEEP(IXSZ),
     &       IW, LIW)
      ENDIF
      PIVSIZ = 1
      NPIV    = IW(IOLDPS+1+XSIZE)
      NPIVP1  = NPIV + 1
      APOSMAX = POSELT+LDA8*LDA8-1_8
      IF(INOPV .EQ. -1) THEN
         APOS = POSELT + (LDA8+1_8) * int(NPIV,8)
         POSPV1 = APOS
         CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( abs(A(APOS)), DKEEP, KEEP, .TRUE.)
         IF(abs(A(APOS)).LT.SEUIL) THEN
            IF(dble(A(APOS)) .GE. RZERO) THEN
               A(APOS) = CSEUIL
            ELSE
               A(APOS) = -CSEUIL
               NNEGW = NNEGW+1
            ENDIF
            NBTINYW = NBTINYW + 1
         ELSE IF (KEEP(258) .NE. 0) THEN
            CALL DMUMPS_UPDATEDETER( A(APOS), DET_MANTW, DET_EXPW )
         ENDIF
              IF (KEEP(50).NE.1 .AND. OOC_EFFECTIVE_ON_FRONT) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, NPIVP1,
     &               PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
              ENDIF
         GO TO 420
      ENDIF
      INOPV   = 0
      ISHIFT = 0              
      IPIV_END = IEND_BLOCK   
      IF (K206.GE.1) THEN
        IF (Inextpiv.GT.NPIVP1.AND.Inextpiv.LE.IEND_BLOCK) THEN
          ISHIFT = Inextpiv - NPIVP1
        ENDIF
        IF ( K206.EQ.1
     &      .OR.  (K206 .GT.1 .AND. IEND_BLOCK.EQ.IEND_BLR) ) THEN
          IPIV_END = IEND_BLOCK + ISHIFT
        ENDIF
        IF (ISHIFT.GT.0.AND.IS_MAXFROMM_AVAIL) THEN
           IPIV = NPIVP1
           APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIV,8)
           POSPV1 = APOS + int(IPIV - NPIVP1,8)
           PIVOT = A(POSPV1)
           IF ( MAXFROMM .GT. PIVNUL ) THEN
               IF ( abs(PIVOT) .GE. UULOC*MAXFROMM
     &         .AND. abs(PIVOT) .GT. max(SEUIL,tiny(MAXFROMM)) ) THEN
                 ISHIFT = 0
               ENDIF
           ENDIF
        ENDIF  
        IF ( ISHIFT .GT. 0) THEN
          IS_MAXFROMM_AVAIL = .FALSE.
        ENDIF
      ENDIF  
       DO 460 IPIV_SHIFT = NPIVP1+ISHIFT, IPIV_END
            IF (IPIV_SHIFT .LE. IEND_BLOCK) THEN
              IPIV=IPIV_SHIFT
            ELSE
              IPIV = IPIV_SHIFT-IEND_BLOCK-1+NPIVP1
              IF (IBEG_BLOCK.EQ.NPIVP1) THEN
                EXIT
              ENDIF
            ENDIF
         APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIV,8)
         POSPV1 = APOS + int(IPIV - NPIVP1,8)
         PIVOT = A(POSPV1)
         IF (UULOC.EQ.RZERO.OR.PIVOT_OPTION.EQ.0) THEN 
            IF (abs(A(APOS)).EQ.RZERO) GO TO 630
            IF (A(APOS).LT.RZERO) NNEGW = NNEGW+1
            CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( abs(A(APOS)), 
     &             DKEEP, KEEP, .FALSE.)
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER(A(APOS), DET_MANTW, DET_EXPW )
            ENDIF
            GO TO 420
         ENDIF
         IF ( IS_MAXFROMM_AVAIL ) THEN
            IF ( MAXFROMM .GT. PIVNUL ) THEN
               IF ( abs(PIVOT) .GE. UULOC*MAXFROMM
     &         .AND. abs(PIVOT) .GT. max(SEUIL,tiny(MAXFROMM)) ) THEN
                 IF (PIVOT .LT. RZERO) NNEGW = NNEGW+1
                 CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &                ( abs(PIVOT), 
     &                  DKEEP, KEEP, .FALSE.)
                 IF (KEEP(258) .NE. 0) THEN
                   CALL DMUMPS_UPDATEDETER(PIVOT, DET_MANTW, DET_EXPW )
                 ENDIF
                 GOTO 415
               ENDIF
            ENDIF
            IS_MAXFROMM_AVAIL = .FALSE.
         ENDIF
         AMAX = -RONE
         JMAX = 0
         IF (PIVOT_OPTION.EQ.3
     &      ) THEN
           LIM = NFRONT - KEEP(253)-NVSCHUR
         ELSEIF (PIVOT_OPTION.GE.2
     &          ) THEN
           LIM = NASS
         ELSEIF (PIVOT_OPTION.GE.1) THEN
           LIM = IEND_BLR
         ELSE
           write(*,*) 'Internal error in FAC_I_LDLT 1x1:',
     &                PIVOT_OPTION
           CALL MUMPS_ABORT()
         ENDIF
         J1 = APOS
         J2 = POSPV1 - 1_8
         DO JJ=J1,J2
            IF(abs(A(JJ)) .GT. AMAX) THEN
               AMAX = abs(A(JJ))
               JMAX = IPIV - int(POSPV1-JJ)
            ENDIF
         ENDDO
         J1 = POSPV1 + LDA8
         DO J=1, IEND_BLOCK - IPIV
            IF(abs(A(J1)) .GT. AMAX) THEN
               AMAX = abs(A(J1))
               JMAX = IPIV + J
            ENDIF
            J1 = J1 + LDA8
         ENDDO
           RMAX = RZERO
           J1_ini = J1
#if defined(_OPENMP)
           J1_end = LIM - IEND_BLOCK
           CHUNK = max(J1_end,1)
           IF ( J1_end.GE.KEEP(360)) THEN
             OMP_FLAG = .TRUE.
             CHUNK = max(KEEP(360)/2,(J1_end+NOMP-1)/NOMP)
           ELSE
             OMP_FLAG = .FALSE.
           ENDIF
#endif
!$OMP PARALLEL DO SCHEDULE(static, CHUNK) PRIVATE(J1)
!$OMP& REDUCTION(max:RMAX) IF(OMP_FLAG)
           DO J=1, LIM - IEND_BLOCK
              J1 = J1_ini + int(J-1,8) * LDA8
              RMAX = max(abs(A(J1)),RMAX)
           ENDDO
!$OMP END PARALLEL DO
            IF (PARPIV_T1.NE.0) THEN
             RMAX_NORELAX = dble(A(APOSMAX+int(IPIV,8)))
            ELSE
             RMAX_NORELAX = RZERO
            ENDIF
            RMAX         = max(RMAX,RMAX_NORELAX)
         IF (max(AMAX,RMAX,abs(PIVOT)).LE.PIVNUL) THEN
               IF ((PARPIV_T1.NE.0)
     &             .AND.(PARPIV_T1.NE.-1)
     &             .AND.(RMAX_NORELAX.LT.0)
     &             .AND.(IPIV.GT.1)) THEN
                 MAX_PREV_in_PARPIV = RZERO
                 DO JJ=1,IPIV-1
                  MAX_PREV_in_PARPIV= max ( MAX_PREV_in_PARPIV, 
     &                dble(A(APOSMAX+int(JJ,8))) )
                 ENDDO
                 IF (MAX_PREV_in_PARPIV.GT.PIVNUL) THEN
                   APOSROW = POSELT + NFRONT8*int(IPIV-1,8)
                   DO JJ=1,IPIV-1
                     IF (abs(A(APOSROW+JJ-1)).GT.PIVNUL) THEN
                       GOTO 460
                     ENDIF
                   ENDDO
                 ENDIF
               ENDIF
            CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( abs(A(POSPV1)), DKEEP, KEEP, .TRUE.)
!$OMP ATOMIC CAPTURE
            KEEP(109) = KEEP(109) + 1
            IPIVNUL = KEEP(109)
!$OMP END ATOMIC
            PIVNUL_LIST(IPIVNUL) = IW( IOLDPS+HF+NPIV+IPIV-NPIVP1 )
            IF(dble(FIXA).GT.RZERO) THEN
               IF(dble(PIVOT) .GE. RZERO) THEN
                  A(POSPV1) = FIXA
               ELSE
                  A(POSPV1) = -FIXA
               ENDIF
            ELSE
               J1 = APOS
               J2 = POSPV1 - 1_8
               DO JJ=J1,J2
                  A(JJ) = ZERO
               ENDDO
               J1 = POSPV1 + LDA8
               DO J=1, IEND_BLOCK - IPIV
                  A(J1) = ZERO
                  J1 = J1 + LDA8
               ENDDO
               DO J=1,LIM - IEND_BLOCK
                  A(J1) = ZERO
                  J1 = J1 + LDA8
               ENDDO
               A(POSPV1) = ONE
            ENDIF
            PIVOT = A(POSPV1)
            GO TO 415
         ENDIF
            RMAX         = max(RMAX,abs(RMAX_NORELAX))
         IF ( abs(PIVOT).GE.UULOC*max(RMAX,AMAX)
     &        .AND. abs(PIVOT) .GT. max(SEUIL,tiny(RMAX)) ) THEN
               IF (PIVOT .LT. ZERO) NNEGW = NNEGW+1
               CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &                ( abs(PIVOT), 
     &                  DKEEP, KEEP, .FALSE.)
               IF (KEEP(258) .NE.0 ) THEN
                 CALL DMUMPS_UPDATEDETER(PIVOT, DET_MANTW, DET_EXPW )
               ENDIF
               GO TO 415
         END IF
         IF (NPIVP1.EQ.IEND_BLOCK) THEN
           GOTO 460
         ELSE IF (JMAX.EQ.0) THEN
           GOTO 460
         ENDIF
         IF (max(abs(PIVOT),RMAX,AMAX).LE.tiny(RMAX)) THEN
           GOTO 460
         ENDIF
         IF (
     &   (KEEP(19).NE.0).AND.(max(AMAX,RMAX,abs(PIVOT)).LE.SEUIL)
     &      )
     &       THEN
           GO TO 460
         ENDIF
           IF (RMAX.LT.AMAX) THEN
             J1 = APOS
             J2 = POSPV1 - 1_8
             DO JJ=J1,J2
               IF(int(POSPV1-JJ) .NE. IPIV-JMAX) THEN
                 RMAX = max(RMAX,abs(A(JJ)))
               ENDIF
             ENDDO
             J1 = POSPV1 + LDA8
             DO J=1,IEND_BLOCK-IPIV
               IF(IPIV+J .NE. JMAX) THEN
                 RMAX = max(abs(A(J1)),RMAX)
               ENDIF
               J1 = J1 + LDA8
             ENDDO
           ENDIF
           APOSJ = POSELT + int(JMAX-1,8)*LDA8 + int(NPIV,8)
           POSPV2 = APOSJ + int(JMAX - NPIVP1,8)
           IF (IPIV.LT.JMAX) THEN
              OFFDAG = APOSJ + int(IPIV - NPIVP1,8)
           ELSE
              OFFDAG = APOS + int(JMAX - NPIVP1,8)
           END IF
           TMAX = RZERO
#if defined(_OPENMP)
           J1_end = LIM-JMAX
           CHUNK = max(J1_end,1)
           IF (J1_end.GE.KEEP(360)) THEN
             OMP_FLAG = .TRUE.
             CHUNK = max(KEEP(360)/2,(J1_end+NOMP-1)/NOMP)
           ELSE
             OMP_FLAG = .FALSE.
           ENDIF
#endif
           IF (JMAX .LT. IPIV) THEN
              JJ_ini = POSPV2
!$OMP PARALLEL DO SCHEDULE(static, CHUNK) IF (OMP_FLAG)
!$OMP& PRIVATE(JJ) REDUCTION(max:TMAX)
              DO K = 1, LIM - JMAX
                 JJ = JJ_ini+ int(K,8)*NFRONT8
                 IF (JMAX+K.NE.IPIV) THEN
                    TMAX=max(TMAX,abs(A(JJ)))
                 ENDIF
              ENDDO
!$OMP END PARALLEL DO
              DO KK =  APOSJ, POSPV2-1_8
                 TMAX = max(TMAX,abs(A(KK)))
              ENDDO
           ELSE
              JJ_ini = POSPV2
!$OMP PARALLEL DO SCHEDULE(static, CHUNK) PRIVATE(JJ) 
!$OMP& REDUCTION(max:TMAX) IF(OMP_FLAG)
              DO K = 1, LIM-JMAX
                 JJ = JJ_ini + int(K,8)*NFRONT8
                 TMAX=max(TMAX,abs(A(JJ)))
              ENDDO
!$OMP END PARALLEL DO
              DO KK =  APOSJ, POSPV2 - 1_8
                 IF (KK.NE.OFFDAG) THEN
                    TMAX = max(TMAX,abs(A(KK)))
                 ENDIF
              ENDDO
           ENDIF
           IF (PARPIV_T1.NE.0) THEN
             TMAX_NORELAX = max(SEUIL*UULOCM1, 
     &                          abs(dble(A(APOSMAX+int(JMAX,8))))
     &                          )
           ELSE
             TMAX_NORELAX = SEUIL*UULOCM1
           ENDIF
           TMAX = max (TMAX,TMAX_NORELAX)
           DETPIV = A(POSPV1)*A(POSPV2) - A(OFFDAG)**2
           IF (SEUIL.GT.RZERO) THEN
                IF (sqrt(abs(DETPIV)) .LE. SEUIL ) THEN
                  GOTO 460
                ENDIF
           ENDIF
           MAXPIV = max(abs(A(POSPV1)),abs(A(POSPV2)))
           IF (MAXPIV.EQ.RZERO) MAXPIV = RONE
           IF ((abs(A(POSPV2))*RMAX+AMAX*TMAX)*UULOC.GT.
     &          abs(DETPIV) .OR. abs(DETPIV) .EQ. RZERO) THEN
             GO TO 460
           ENDIF
           IF ((abs(A(POSPV1))*TMAX+AMAX*RMAX)*UULOC.GT.
     &          abs(DETPIV) .OR. abs(DETPIV) .EQ. RZERO) THEN
             GO TO 460
           ENDIF
           CALL DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( sqrt(abs(DETPIV)),
     &             DKEEP, KEEP, .FALSE.)
           IF (KEEP(258) .NE.0 ) THEN
             CALL DMUMPS_UPDATEDETER(DETPIV, DET_MANTW, DET_EXPW )
           ENDIF
           PIVSIZ = 2
           NB22T1W = NB22T1W + 1
           IF(DETPIV .LT. RZERO) THEN
             NNEGW = NNEGW+1
           ELSE IF(A(POSPV2) .LT. RZERO) THEN
             NNEGW = NNEGW+2
           ENDIF
 415       CONTINUE
          IF (K206.GE.1) THEN
             Inextpiv = max(NPIVP1+PIVSIZ, IPIV+1)
          ENDIF
           DO K=1,PIVSIZ
              IF (PIVSIZ .EQ. 2) THEN
                IF (K==1) THEN
                  LPIV = min(IPIV,JMAX)
                ELSE
                  LPIV   = max(IPIV,JMAX)
                ENDIF
              ELSE
                LPIV = IPIV
              ENDIF
              IF (LPIV.EQ.NPIVP1) GOTO 416 
              IF (KEEP(405) .EQ. 0) THEN
                KEEP8(80) = KEEP8(80)+1
              ELSE
!$OMP           ATOMIC UPDATE
                KEEP8(80) = KEEP8(80)+1
!$OMP           END ATOMIC
              ENDIF
              LIM_SWAP = NFRONT
              CALL DMUMPS_SWAP_LDLT( A, LA, IW, LIW,
     &             IOLDPS, NPIVP1, LPIV, POSELT, LIM_SWAP,
     &             LDA, NFRONT, 1, PARPIV_T1, KEEP(50),
     &             KEEP(IXSZ), -9999)
 416          CONTINUE
              IF (KEEP(50).NE.1 .AND. OOC_EFFECTIVE_ON_FRONT) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, LPIV, PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
              ENDIF
              NPIVP1 = NPIVP1 + 1
           ENDDO
           IF(PIVSIZ .EQ. 2) THEN
              A(POSELT+(LDA8+1_8)*int(NPIV,8)+1_8) = DETPIV
           ENDIF
           GOTO 420
  460   CONTINUE
          IF (K206 .GE. 1) THEN
            Inextpiv=IEND_BLOCK+1
          ENDIF
      IF (IEND_BLOCK.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      PIVSIZ = 0
      IFLAG = -10
  420 CONTINUE
      IS_MAXFROMM_AVAIL = .FALSE.
      RETURN
      END SUBROUTINE DMUMPS_FAC_I_LDLT
      SUBROUTINE DMUMPS_FAC_MQ_LDLT(IEND_BLOCK,
     &     NFRONT,NASS,NPIV,INODE,
     &     A,LA,LDA, 
     &     POSELT,IFINB,PIVSIZ,
     &     MAXFROMM, IS_MAXFROMM_AVAIL, IS_MAX_USEFUL,
     &     PARPIV_T1, LAST_ROW, IEND_BLR, NVSCHUR_K253,
     &     LR_ACTIVATED
     &     )
      IMPLICIT NONE
      INTEGER, intent(out):: IFINB
      INTEGER, intent(in) :: INODE, NFRONT, NASS, NPIV
      INTEGER, intent(in) :: IEND_BLOCK
      INTEGER, intent(in) :: LDA
      INTEGER(8), intent(in) :: LA
      DOUBLE PRECISION, intent(inout) :: A(LA)
      INTEGER, intent(in)    :: LAST_ROW
      INTEGER, intent(in)    :: IEND_BLR
      INTEGER(8) :: POSELT
      DOUBLE PRECISION, intent(out) :: MAXFROMM
      LOGICAL, intent(out) :: IS_MAXFROMM_AVAIL
      LOGICAL, intent(in) :: IS_MAX_USEFUL
      INTEGER, intent(in) :: PARPIV_T1
      INTEGER, INTENT(in) :: NVSCHUR_K253
      LOGICAL, intent(in) :: LR_ACTIVATED
      DOUBLE PRECISION    VALPIV
      DOUBLE PRECISION :: MAXFROMMTMP
      INTEGER  NCB1
      INTEGER(8) :: NFRONT8
      INTEGER(8) :: LDA8
      INTEGER(8) :: K1POS
      INTEGER NEL2, NEL
      DOUBLE PRECISION ONE, ZERO
      DOUBLE PRECISION A11,A22,A12
      INTEGER(8) :: APOS, LPOS, LPOS1, LPOS2
      INTEGER(8) :: POSPV1, POSPV2
      INTEGER PIVSIZ,NPIV_NEW,J2,I
      INTEGER(8) :: OFFDAG, OFFDAG_OLD, IBEG, IEND
      INTEGER(8) :: JJ, K1, K2, IROW
      INTEGER(8) :: ROW_SHIFT, JJ_LOC, IBEG_LOC, IEND_LOC
      DOUBLE PRECISION SWOP,DETPIV,MULT1,MULT2
      INTEGER(8) :: APOSMAX
      INCLUDE 'mumps_headers.h'
      PARAMETER(ONE  = 1.0D0,
     &          ZERO = 0.0D0)
      LDA8     = int(LDA,8)
      NFRONT8  = int(NFRONT,8)
      NPIV_NEW = NPIV + PIVSIZ
      NEL      = NFRONT - NPIV_NEW
      IFINB    = 0
      IS_MAXFROMM_AVAIL = .FALSE.
      NEL2   = IEND_BLOCK - NPIV_NEW
      IF (NEL2.EQ.0) THEN
        IF (IEND_BLOCK.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      MAXFROMM = 0.0D0
      IF(PIVSIZ .EQ. 1) THEN
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + LDA8
         IF (NEL2 > 0) THEN
           IF (.NOT. IS_MAX_USEFUL) THEN
             DO I=1, NEL2
               K1POS = LPOS + int(I-1,8)*LDA8
               A(APOS+int(I,8))=A(K1POS)
               A(K1POS) = A(K1POS) * VALPIV
               DO JJ=1_8, int(I,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDDO
           ELSE
             IS_MAXFROMM_AVAIL = .TRUE.
             DO I=1, NEL2
               K1POS = LPOS + int(I-1,8)*LDA8
               A(APOS+int(I,8))=A(K1POS)
               A(K1POS) = A(K1POS) * VALPIV
               A(K1POS+1_8)=A(K1POS+1_8) - A(K1POS) * A(APOS+1_8)
               MAXFROMM=max( MAXFROMM,abs(A(K1POS+1_8)) )
               DO JJ = 2_8, int(I,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDDO
           ENDIF
         ENDIF
         NCB1 = LAST_ROW - IEND_BLOCK
         IF (NCB1.GT.0) THEN
          IF (.NOT. IS_MAX_USEFUL) THEN
!$OMP       PARALLEL DO PRIVATE(JJ,K1POS) IF (NCB1 > 300)
            DO I=NEL2+1, NEL2 + NCB1
              K1POS = LPOS+ int(I-1,8)*LDA8
              A(APOS+int(I,8))=A(K1POS)
              A(K1POS) = A(K1POS) * VALPIV
              DO JJ = 1_8, int(NEL2,8)
                A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
          ELSE
            MAXFROMMTMP=0.0D0
!$OMP       PARALLEL DO PRIVATE(JJ,K1POS)
!$OMP&      REDUCTION(max:MAXFROMMTMP) IF (NCB1-NVSCHUR_K253>300)
            DO I=NEL2+1, NEL2 + NCB1 - NVSCHUR_K253
              K1POS = LPOS+ int(I-1,8)*LDA8
              A(APOS+int(I,8))=A(K1POS)
              A(K1POS) = A(K1POS) * VALPIV
              IF (NEL2 > 0) THEN
                A(K1POS+1_8) = A(K1POS+1_8) - A(K1POS) * A(APOS+1_8)
                MAXFROMMTMP=max(MAXFROMMTMP, abs(A(K1POS+1_8)))
                DO JJ = 2_8, int(NEL2,8)
                  A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
                ENDDO
              ENDIF
            ENDDO
!$OMP       END PARALLEL DO
            DO I = NEL2 + NCB1 - NVSCHUR_K253 + 1, NEL2 + NCB1
              K1POS = LPOS+ int(I-1,8)*LDA8
              A(APOS+int(I,8))=A(K1POS)
              A(K1POS) = A(K1POS) * VALPIV
              DO JJ = 1_8, int(NEL2,8)
                A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
              ENDDO
            ENDDO
            MAXFROMM=max(MAXFROMM, MAXFROMMTMP)
          ENDIF
         ENDIF
      ELSE
         POSPV1 = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         POSPV2 = POSPV1 + NFRONT8 + 1_8
         OFFDAG_OLD = POSPV2 - 1_8
         OFFDAG = POSPV1 + 1_8
         SWOP = A(POSPV2)
         DETPIV = A(OFFDAG)
          A22 = A(POSPV1)/DETPIV   
          A11 =  SWOP/DETPIV       
          A12 = -A(OFFDAG_OLD)/DETPIV   
          A(OFFDAG)     = A(OFFDAG_OLD)  
          A(OFFDAG_OLD) = ZERO
         LPOS1   = POSPV2 + LDA8 - 1_8
         LPOS2   = LPOS1 + 1_8
         CALL dcopy(LAST_ROW-NPIV_NEW, A(LPOS1), LDA, A(POSPV1+2_8), 1)
         CALL dcopy(LAST_ROW-NPIV_NEW, A(LPOS2), LDA, A(POSPV2+1_8), 1)
         JJ = POSPV2 + NFRONT8-1_8  
         IBEG = JJ + 2_8
         IEND = IBEG
         DO J2 = 1,NEL2
           K1 = JJ
           K2 = JJ+1_8
           MULT1 = - (A11*A(K1)+A12*A(K2))
           MULT2 = - (A12*A(K1)+A22*A(K2))
           K1 = POSPV1 + 2_8
           K2 = POSPV2 + 1_8
           DO IROW = IBEG, IEND
              A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
              K1 = K1 + 1_8
              K2 = K2 + 1_8
           ENDDO
           A( JJ       ) = -MULT1
           A( JJ + 1_8 ) = -MULT2
           IBEG = IBEG + NFRONT8
           IEND = IEND + NFRONT8 + 1_8
           JJ = JJ+NFRONT8
         ENDDO
         IEND = IEND-1_8
!$OMP    PARALLEL DO PRIVATE(J2, K1, K2, MULT1, MULT2, IROW, JJ_LOC,
!$OMP&   ROW_SHIFT, IBEG_LOC, IEND_LOC) IF (LAST_ROW-IEND_BLOCK>300)
         DO J2 = 1,LAST_ROW-IEND_BLOCK
           ROW_SHIFT = (J2-1_8)*NFRONT8
           JJ_LOC = JJ + ROW_SHIFT
           IBEG_LOC = IBEG + ROW_SHIFT
           IEND_LOC = IEND + ROW_SHIFT
           K1 = JJ_LOC
           K2 = JJ_LOC+1_8
           MULT1 = - (A11*A(K1)+A12*A(K2))
           MULT2 = - (A12*A(K1)+A22*A(K2))
           K1 = POSPV1 + 2_8
           K2 = POSPV2 + 1_8
           DO IROW = IBEG_LOC, IEND_LOC
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1_8
               K2 = K2 + 1_8
           ENDDO
           A( JJ_LOC       ) = -MULT1
           A( JJ_LOC + 1_8 ) = -MULT2
         ENDDO
!$OMP    END PARALLEL DO
      ENDIF
      IF ((IS_MAXFROMM_AVAIL).AND.(NEL2.GT.0)) THEN
          IF (PARPIV_T1.NE.0) THEN
             APOSMAX = POSELT+LDA8*LDA8-1_8 + int(NPIV_NEW+1,8)
             MAXFROMM = max(MAXFROMM,
     &                      dble(A(APOSMAX))
     &                      )
          ENDIF
       ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FAC_MQ_LDLT
      SUBROUTINE DMUMPS_FAC_SQ_LDLT(IBEG_BLOCK,IEND_BLOCK,NPIV,
     &    NFRONT,NASS,INODE,A,LA,
     &    LDA,
     &    POSELT,
     &    KEEP,KEEP8, 
     &    FIRST_ROW_TRSM, LAST_ROW_TRSM, 
     &    LAST_COL_GEMM, LAST_ROW_GEMM, 
     &    CALL_TRSM, CALL_GEMM, LR_ACTIVATED,
     &    IW, LIW, OFFSET_IW
     &    )
      IMPLICIT NONE
      INTEGER, intent(in) :: NPIV
      INTEGER, intent(in) :: NFRONT, NASS, IBEG_BLOCK, IEND_BLOCK
      INTEGER(8), intent(in) :: LA
      DOUBLE PRECISION, intent(inout) :: A(LA)
      INTEGER, intent(in) :: INODE
      INTEGER    :: KEEP(500)
      INTEGER(8) :: KEEP8(150)
      INTEGER(8), intent(in) :: POSELT
      INTEGER, intent(in) :: LDA
      INTEGER, intent(in) :: LAST_COL_GEMM
      INTEGER, intent(in) :: LAST_ROW_GEMM, LAST_ROW_TRSM,
     &                       FIRST_ROW_TRSM
      LOGICAL, intent(in) :: CALL_TRSM, CALL_GEMM, LR_ACTIVATED
      INTEGER :: OFFSET_IW, LIW
      INTEGER :: IW(LIW)
      INTEGER(8) :: LDA8
      INTEGER NPIV_BLOCK, NEL1
      INTEGER NRHS_TRSM
      INTEGER(8) :: LPOS, UPOS, APOS
      INTEGER IROW
      INTEGER Block
      INTEGER BLSIZE
      DOUBLE PRECISION ONE, ALPHA
      INCLUDE 'mumps_headers.h'
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0)
      LDA8 = int(LDA,8)
      NEL1 = LAST_COL_GEMM - IEND_BLOCK
      NRHS_TRSM = LAST_ROW_TRSM-FIRST_ROW_TRSM
      NPIV_BLOCK  = NPIV - IBEG_BLOCK + 1
      IF (NPIV_BLOCK.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF (CALL_TRSM) THEN
          APOS = POSELT + LDA8*int(IBEG_BLOCK-1,8) + int(IBEG_BLOCK-1,8)
          LPOS = POSELT + LDA8*int(FIRST_ROW_TRSM,8)+int(IBEG_BLOCK-1,8)
          UPOS = POSELT + LDA8*int(IBEG_BLOCK-1,8)+int(FIRST_ROW_TRSM,8)
          CALL dtrsm('L', 'U', 'T', 'U', NPIV_BLOCK, NRHS_TRSM,
     &              ONE, A(APOS), LDA, A(LPOS), LDA)
          CALL DMUMPS_FAC_LDLT_COPY2U_SCALEL(NRHS_TRSM, 1, KEEP(424),
     &             NFRONT, NPIV_BLOCK, LIW, IW, OFFSET_IW, LA, A,
     &             POSELT, LPOS, UPOS, APOS, .NOT.LR_ACTIVATED)
        ENDIF
        IF (CALL_GEMM) THEN
#if defined(GEMMT_AVAILABLE)
         IF ( KEEP(421).EQ. -1) THEN
           LPOS = POSELT + LDA8*int(IEND_BLOCK,8) + int(IBEG_BLOCK-1,8)
           UPOS = POSELT + LDA8*int(IBEG_BLOCK-1,8) + int(IEND_BLOCK,8)
           APOS = POSELT + LDA8*int(IEND_BLOCK,8) + int(IEND_BLOCK,8)
           CALL dgemmt( 'U','N','N', NEL1,
     &                NPIV_BLOCK,
     &                ALPHA, A( UPOS ), LDA,
     &                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         ELSE
#endif
          IF ( LAST_COL_GEMM - IEND_BLOCK > KEEP(7) ) THEN
           BLSIZE = KEEP(8)
          ELSE
           BLSIZE = LAST_COL_GEMM - IEND_BLOCK
          END IF
          IF ( LAST_COL_GEMM - IEND_BLOCK .GT. 0 ) THEN
#if defined(SAK_BYROW)
           DO IROW = IEND_BLOCK+1, LAST_COL_GEMM, BLSIZE
            Block = min( BLSIZE, LAST_COL_GEMM - IROW + 1 )
            LPOS = POSELT + int(IROW  - 1,8) * LDA8 +
     &                      int(IBEG_BLOCK - 1,8)
            UPOS = POSELT + int(IBEG_BLOCK - 1,8) * LDA8 +
     &                      int(IROW - 1,8)
            APOS = POSELT + int(IROW  - 1,8) * LDA8 + 
     &             int(IEND_BLOCK,8)
            CALL dgemm( 'N','N', IROW + Block - IEND_BLOCK - 1, 
     &             Block, NPIV_BLOCK,
     &                 ALPHA, A( UPOS ), LDA,
     &                 A( LPOS ), LDA, ONE, A( APOS ), LDA )
           ENDDO
#else
           DO IROW = IEND_BLOCK+1, LAST_COL_GEMM, BLSIZE
            Block = min( BLSIZE, LAST_COL_GEMM - IROW + 1 )
            LPOS = POSELT + int( IROW - 1,8) * LDA8 +
     &                      int(IBEG_BLOCK - 1,8)
            UPOS = POSELT + int(IBEG_BLOCK - 1,8) * LDA8 +
     &                      int( IROW - 1,8)
            APOS = POSELT + int( IROW - 1,8) * LDA8 + int( IROW - 1,8)
            CALL dgemm( 'N','N', Block, LAST_COL_GEMM - IROW + 1, 
     &                 NPIV_BLOCK, ALPHA, A( UPOS ), LDA,
     &                 A( LPOS ), LDA, ONE, A( APOS ), LDA )
           END DO
#endif
          END IF
#if defined(GEMMT_AVAILABLE)
         END IF 
#endif
         LPOS = POSELT + int(LAST_COL_GEMM,8)*LDA8 + int(IBEG_BLOCK-1,8)
         UPOS = POSELT + int(IBEG_BLOCK-1,8) * LDA8 +
     &                 int(IEND_BLOCK,8)
         APOS = POSELT + int(LAST_COL_GEMM,8)*LDA8 + int(IEND_BLOCK,8)
         IF (LAST_ROW_GEMM .GT. LAST_COL_GEMM) THEN
           CALL dgemm('N', 'N', NEL1, LAST_ROW_GEMM-LAST_COL_GEMM, 
     &              NPIV_BLOCK, ALPHA, A(UPOS), LDA, A(LPOS), LDA, 
     &              ONE, A(APOS), LDA)
         ENDIF
        ENDIF
      ENDIF
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_SQ_LDLT
        SUBROUTINE DMUMPS_SWAP_LDLT( A, LA, IW, LIW,
     &                       IOLDPS, NPIVP1, IPIV, POSELT, LASTROW2SWAP,
     &                       LDA, NFRONT, LEVEL, PARPIV, K50, XSIZE,
     &                       IBEG_BLOCK_TO_SEND )
        IMPLICIT NONE
      INTEGER(8) :: POSELT, LA
      INTEGER LIW, IOLDPS, NPIVP1, IPIV
      INTEGER LDA, NFRONT, LEVEL, PARPIV, K50, XSIZE
      INTEGER LASTROW2SWAP
      DOUBLE PRECISION A( LA )
      INTEGER IW( LIW )
      INTEGER, INTENT(IN) :: IBEG_BLOCK_TO_SEND
      INCLUDE 'mumps_headers.h'
      INTEGER :: IBEG
      INTEGER ISW, ISWPS1, ISWPS2, HF
      INTEGER(8) :: IDIAG, APOS
      INTEGER(8) :: LDA8
      DOUBLE PRECISION SWOP
            LDA8 = int(LDA,8)
            APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIVP1-1,8)
            IDIAG = APOS + int(IPIV - NPIVP1,8)
            HF = 6 + IW( IOLDPS + 5 + XSIZE) + XSIZE
            ISWPS1 = IOLDPS + HF + NPIVP1 - 1
            ISWPS2 = IOLDPS + HF + IPIV - 1
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            ISW = IW(ISWPS1+NFRONT)
            IW(ISWPS1+NFRONT) = IW(ISWPS2+NFRONT)
            IW(ISWPS2+NFRONT) = ISW
            IF ( LEVEL .eq. 2 ) THEN
              IBEG = IBEG_BLOCK_TO_SEND
              CALL dswap( NPIVP1 - 1 - IBEG + 1,
     &            A( POSELT + int(NPIVP1-1,8) +
     &                     int(IBEG-1,8) * LDA8), LDA,
     &            A( POSELT + int(IPIV-1,8)   +
     &                     int(IBEG-1,8) * LDA8), LDA )
            END IF
            CALL dswap( NPIVP1-1,
     &           A( POSELT+int(NPIVP1-1,8) * LDA8 ), 1,
     &           A( POSELT + int(IPIV-1,8) * LDA8 ), 1 )
             CALL dswap( IPIV - NPIVP1 - 1,
     &           A( POSELT+int(NPIVP1,8) * LDA8 + int(NPIVP1-1,8) ),
     &           LDA, A( APOS + 1_8 ), 1 )
            SWOP = A(IDIAG)
            A(IDIAG) = A( POSELT+int(NPIVP1-1,8)*LDA8+int(NPIVP1-1,8) )
            A( POSELT + int(NPIVP1-1,8)*LDA8 + int(NPIVP1-1,8) ) = SWOP
            CALL dswap( LASTROW2SWAP - IPIV,
     &      A( APOS  + LDA8 ), LDA,
     &      A( IDIAG + LDA8 ), LDA )
            IF (PARPIV.NE.0 .AND.K50.EQ.2) THEN
             IF ( LEVEL .eq. 2 .OR. LEVEL.eq.1) THEN
              APOS                 = POSELT+LDA8*LDA8-1_8
              SWOP                 = A(APOS+int(NPIVP1,8))
              A(APOS+int(NPIVP1,8))= A(APOS+int(IPIV,8))
              A(APOS+int(IPIV,8))  = SWOP
             ENDIF
            ENDIF
        RETURN
        END SUBROUTINE DMUMPS_SWAP_LDLT
      SUBROUTINE DMUMPS_FAC_LDLT_COPY2U_SCALEL( IROWMAX, IROWMIN,
     &                  SIZECOPY, LDA, NCOLS, LIW, IW, OFFSET_IW,
     &                  LA, A, POSELT, A_LPOS, A_UPOS, A_DPOS,
     &                  COPY_NEEDED )
!$    USE OMP_LIB
      INTEGER, INTENT(IN) :: IROWMAX, IROWMIN
      INTEGER, INTENT(IN) :: SIZECOPY
      INTEGER, INTENT(IN) :: LDA, NCOLS
      INTEGER, INTENT(IN) :: LIW
      INTEGER, INTENT(IN) :: IW(LIW)
      INTEGER, INTENT(IN) :: OFFSET_IW
      INTEGER(8), INTENT(IN) :: LA
      DOUBLE PRECISION, INTENT(INOUT) :: A(LA)
      INTEGER(8), INTENT(IN) :: POSELT, A_LPOS, A_UPOS, A_DPOS
      LOGICAL, INTENT(IN) :: COPY_NEEDED
      INTEGER(8) :: LPOS, UPOS
      INTEGER(8) :: DPOS, POSPV1, POSPV2, OFFDAG
      INTEGER(8) :: LDA8
      INTEGER :: IROWEND, IROW, Block2
      INTEGER :: I, J
      DOUBLE PRECISION :: MULT1, MULT2, A11, DETPIV, A22, A12
      INTEGER :: BLSIZECOPY
      DOUBLE PRECISION :: ONE
      PARAMETER (ONE = 1.0D0)
      INTEGER(8) :: LPOSI, UPOSI
      LOGICAL :: PIVOT_2X2
!$    LOGICAL :: OMP_FLAG
!$    INTEGER :: NOMP, CHUNK
      LDA8 = int(LDA,8)
      IF (SIZECOPY.NE.0) THEN
        BLSIZECOPY = SIZECOPY
      ELSE
        BLSIZECOPY = 250
      ENDIF
!$    NOMP = OMP_GET_MAX_THREADS()
!$    OMP_FLAG = .FALSE.
!$    CHUNK = (64/4)
!$    IF (NOMP .GT. 1 .AND. NCOLS .GE. 4*CHUNK) THEN
!$      OMP_FLAG = .TRUE.
!$      CHUNK = max(2*CHUNK, NCOLS/NOMP)
!$    ENDIF
      DO IROWEND = IROWMAX, IROWMIN, -BLSIZECOPY
        Block2 = min(BLSIZECOPY, IROWEND)
        IROW = IROWEND - Block2 + 1
        LPOS = A_LPOS + int(IROW-1,8)*LDA8
        UPOS = A_UPOS + int(IROW-1,8)
!$OMP  PARALLEL DO PRIVATE(PIVOT_2X2, A11, DPOS,
!$OMP&   POSPV1, POSPV2, OFFDAG, A22, A12, DETPIV, J, MULT1, MULT2
!$OMP&   , LPOSI, UPOSI
!$OMP&   ) FIRSTPRIVATE(Block2, LDA, LDA8, LPOS, UPOS, A_DPOS)
!$OMP& SCHEDULE(STATIC,CHUNK) IF(OMP_FLAG)
        DO I=1, NCOLS
          PIVOT_2X2 = .FALSE.
          IF(IW(OFFSET_IW+I-1) .LE. 0) THEN
            PIVOT_2X2 = .TRUE.
          ELSE
            IF (I .GT. 1) THEN
              IF (IW(OFFSET_IW+I-2) .LE. 0) THEN
                cycle
              ENDIF
            ENDIF
          ENDIF
          DPOS = A_DPOS + LDA8*int(I-1,8) + int(I-1,8)
          IF(.not. PIVOT_2X2) THEN
            A11 = ONE/A(DPOS)
            LPOSI = LPOS+int(I-1,8)
            IF (COPY_NEEDED) THEN
              UPOSI = UPOS+int(I-1,8)*LDA8
              DO J = 1, Block2
                A(UPOSI+int(J-1,8)) = A(LPOSI+int(J-1,8)*LDA8)
              END DO
            ENDIF
            DO J = 1, Block2
              A(LPOSI+int(J-1,8)*LDA8) = A(LPOSI+int(J-1,8)*LDA8)*A11
            END DO
          ELSE
            IF (COPY_NEEDED) THEN
              CALL dcopy(Block2, A(LPOS+int(I-1,8)),
     &                 LDA, A(UPOS+int(I-1,8)*LDA8), 1)
              CALL dcopy(Block2, A(LPOS+int(I,8)),
     &                 LDA, A(UPOS+int(I,8)*LDA8), 1)
            ENDIF
            POSPV1 = DPOS
            POSPV2 = DPOS + int(LDA+1,8)
            OFFDAG = POSPV1+1_8
            A11 = A(POSPV1)
            A22 = A(POSPV2)
            A12 = A(OFFDAG)
            DETPIV = A11*A22 - A12**2
            A22 = A11/DETPIV
            A11 = A(POSPV2)/DETPIV
            A12 = -A12/DETPIV
            DO J = 1,Block2
              MULT1 = A11*A(LPOS+int(J-1,8)*LDA8+int(I-1,8))
     &              + A12*A(LPOS+int(J-1,8)*LDA8+int(I,8))
              MULT2 = A12*A(LPOS+int(J-1,8)*LDA8+int(I-1,8))
     &              + A22*A(LPOS+int(J-1,8)*LDA8+int(I,8))
              A(LPOS+int(J-1,8)*LDA8+int(I-1,8)) = MULT1
              A(LPOS+int(J-1,8)*LDA8+int(I,8))   = MULT2
            ENDDO
          ENDIF
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
      END SUBROUTINE DMUMPS_FAC_LDLT_COPY2U_SCALEL
      SUBROUTINE DMUMPS_FAC_LDLT_COPYSCALE_U( IROWMAX, IROWMIN,
     &                  SIZECOPY, LDA, NCOLS, LIW, IW, OFFSET_IW,
     &                  LA, A, POSELT, A_LPOS, A_UPOS, A_DPOS )
!$    USE OMP_LIB
      INTEGER, INTENT(IN) :: IROWMAX, IROWMIN
      INTEGER, INTENT(IN) :: SIZECOPY
      INTEGER, INTENT(IN) :: LDA, NCOLS
      INTEGER, INTENT(IN) :: LIW
      INTEGER, INTENT(IN) :: IW(LIW)
      INTEGER, INTENT(IN) :: OFFSET_IW
      INTEGER(8), INTENT(IN) :: LA
      DOUBLE PRECISION, INTENT(INOUT) :: A(LA)
      INTEGER(8), INTENT(IN) :: POSELT, A_LPOS, A_UPOS, A_DPOS
      INTEGER(8) :: LPOS, UPOS
      INTEGER(8) :: DPOS, POSPV1, POSPV2, OFFDAG
      INTEGER(8) :: LDA8
      INTEGER :: IROWEND, IROW, Block2
      INTEGER :: I, J
      DOUBLE PRECISION :: MULT1, MULT2, A11, DETPIV, A22, A12
      INTEGER :: BLSIZECOPY
      DOUBLE PRECISION :: ONE
      PARAMETER (ONE = 1.0D0)
      INTEGER(8) :: LPOSI, UPOSI
      LOGICAL :: PIVOT_2X2
!$    LOGICAL :: OMP_FLAG
!$    INTEGER :: NOMP, CHUNK
      LDA8 = int(LDA,8)
      IF (SIZECOPY.NE.0) THEN
        BLSIZECOPY = SIZECOPY
      ELSE
        BLSIZECOPY = 250
      ENDIF
!$    NOMP = OMP_GET_MAX_THREADS()
!$    OMP_FLAG = .FALSE.
!$    CHUNK = (64/4)
!$    IF (NOMP .GT. 1 .AND. NCOLS .GE. 4*CHUNK) THEN
!$      OMP_FLAG = .TRUE.
!$      CHUNK = max(2*CHUNK, NCOLS/NOMP)
!$    ENDIF
      DO IROWEND = IROWMAX, IROWMIN, -BLSIZECOPY
        Block2 = min(BLSIZECOPY, IROWEND)
        IROW = IROWEND - Block2 + 1
        LPOS = A_LPOS + int(IROW-1,8)*LDA8
        UPOS = A_UPOS + int(IROW-1,8)
!$OMP  PARALLEL DO PRIVATE(PIVOT_2X2, A11, DPOS,
!$OMP&   POSPV1, POSPV2, OFFDAG, A22, A12, DETPIV, J, MULT1, MULT2
!$OMP&   , LPOSI, UPOSI
!$OMP&   ) FIRSTPRIVATE(Block2, LDA, LDA8, LPOS, UPOS, POSELT)
!$OMP& SCHEDULE(STATIC,CHUNK) IF(OMP_FLAG)
        DO I=1, NCOLS
          PIVOT_2X2 = .FALSE.
          IF(IW(OFFSET_IW+I-1) .LE. 0) THEN
            PIVOT_2X2 = .TRUE.
          ELSE
            IF (I .GT. 1) THEN
              IF (IW(OFFSET_IW+I-2) .LE. 0) THEN
                cycle
              ENDIF
            ENDIF
          ENDIF
          DPOS = A_DPOS + LDA8*int(I-1,8) + int(I-1,8)
          IF(.not. PIVOT_2X2) THEN
            A11 = A(DPOS)
            LPOSI = LPOS+int(I-1,8)
            UPOSI = UPOS+int(I-1,8)*LDA8
            DO J = 1, Block2
              A(UPOSI+int(J-1,8)) = A(LPOSI+int(J-1,8)*LDA8)*A11
            END DO
          ELSE
            POSPV1 = DPOS
            POSPV2 = DPOS + int(LDA+1,8)
            OFFDAG = POSPV1+1_8
            A11 = A(POSPV1)
            A22 = A(POSPV2)
            A12 = A(OFFDAG)
            DO J = 1,Block2
              MULT1 = A11*A(LPOS+int(J-1,8)*LDA8+int(I-1,8))
     &              + A12*A(LPOS+int(J-1,8)*LDA8+int(I,8))
              MULT2 = A12*A(LPOS+int(J-1,8)*LDA8+int(I-1,8))
     &              + A22*A(LPOS+int(J-1,8)*LDA8+int(I,8))
              A(UPOS+int(I-1,8)*LDA8+int(J-1,8)) = MULT1
              A(UPOS+int(I,8)*LDA8+int(J-1,8))   = MULT2
            ENDDO
          ENDIF
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_FAC_LDLT_COPYSCALE_U
      SUBROUTINE DMUMPS_FAC_T_LDLT(NFRONT,NASS,
     &    IW,LIW,A,LA,
     &    LDA,
     &    IOLDPS,POSELT,KEEP,KEEP8,
     &    POSTPONE_COL_UPDATE, ETATASS,
     &    TYPEFile, LAFAC, MonBloc, NextPiv2beWritten,
     &    LIWFAC, MYID, IFLAG, OFFSET_IW, INODE )
      USE DMUMPS_OOC
      IMPLICIT NONE
      INTEGER NFRONT, NASS,LIW
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW) 
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: POSELT
      INTEGER LDA
      INTEGER IOLDPS, ETATASS
      LOGICAL POSTPONE_COL_UPDATE
      INTEGER(8) :: LAFAC
      INTEGER TYPEFile, NextPiv2beWritten
      INTEGER LIWFAC, MYID, IFLAG
      TYPE(IO_BLOCK):: MonBloc
      INTEGER IDUMMY
      LOGICAL LAST_CALL
      INTEGER :: OFFSET_IW
      INTEGER, intent(in):: INODE
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: UPOS, APOS, LPOS
      INTEGER(8) :: LDA8
      INTEGER BLSIZE, BLSIZE2, Block, IROW, NPIV, IROWEND
      INTEGER I2, I2END, Block2
      DOUBLE PRECISION  ONE, ALPHA, BETA, ZERO
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      LDA8 = int(LDA,8)
      IF (ETATASS.EQ.1) THEN
        BETA = ZERO
      ELSE
        BETA = ONE
      ENDIF
      IF ( NFRONT - NASS > KEEP(58) ) THEN
        IF ( NFRONT - NASS > KEEP(57) ) THEN
          BLSIZE = KEEP(58)
        ELSE
          BLSIZE = (NFRONT - NASS)/2
        END IF
      ELSE
        BLSIZE = NFRONT - NASS
      END IF
      BLSIZE2 = KEEP(218)
      NPIV = IW( IOLDPS + 1 + KEEP(IXSZ))
      IF ( NFRONT - NASS .GT. 0 ) THEN
       IF ( POSTPONE_COL_UPDATE ) THEN
         LPOS = POSELT + LDA8 * int(NASS,8)
         CALL dtrsm( 'L', 'U', 'T', 'U',
     &               NPIV, NFRONT-NASS, ONE,
     &               A( POSELT ), LDA,
     &               A( LPOS ), LDA )
       ENDIF
#if defined(GEMMT_AVAILABLE)
       IF ( KEEP(421).EQ. -1) THEN
         LPOS = POSELT + int(NASS,8)*LDA8
         UPOS = POSELT + int(NASS,8)
         APOS = POSELT + int(NASS,8)*LDA8 + int(NASS,8)
         IF (POSTPONE_COL_UPDATE) THEN
           CALL DMUMPS_FAC_LDLT_COPY2U_SCALEL( NFRONT - NASS, 1,
     &                 KEEP(424), NFRONT, NPIV,
     &                 LIW, IW, OFFSET_IW, LA, A, POSELT, LPOS, UPOS,
     &                 POSELT, .TRUE. )
         ENDIF
         CALL dgemmt('U', 'N', 'N', NFRONT-NASS, NPIV,
     &                ALPHA, A( UPOS ), LDA,
     &                A( LPOS ), LDA,
     &                BETA,
     &                A( APOS ), LDA )
       ELSE
#endif
       DO IROWEND = NFRONT - NASS, 1, -BLSIZE
        Block = min( BLSIZE, IROWEND )
        IROW  = IROWEND - Block + 1
        LPOS = POSELT + int(NASS,8)*LDA8 + int(IROW-1,8) * LDA8
        APOS = POSELT + int(NASS,8)*LDA8 + int(IROW-1,8) * LDA8 +
     &                  int(NASS + IROW - 1,8)
        UPOS = POSELT + int(NASS,8)
        IF (.NOT. POSTPONE_COL_UPDATE) THEN
          UPOS = POSELT + int(NASS + IROW - 1,8)
        ENDIF
        IF (POSTPONE_COL_UPDATE) THEN
          CALL DMUMPS_FAC_LDLT_COPY2U_SCALEL( Block, 1,
     &                KEEP(424), NFRONT, NPIV,
     &                LIW, IW, OFFSET_IW, LA, A, POSELT, LPOS, UPOS,
     &                POSELT, .TRUE. )
        ENDIF
        DO I2END = Block, 1, -BLSIZE2
          Block2 = min(BLSIZE2, I2END)
          I2 = I2END - Block2+1
          CALL dgemm('N', 'N', Block2, Block-I2+1, NPIV, ALPHA,
     &               A(UPOS+int(I2-1,8)), LDA,
     &               A(LPOS+int(I2-1,8)*LDA8), LDA,
     &               BETA,
     &               A(APOS + int(I2-1,8) + int(I2-1,8)*LDA8), LDA)
          IF (KEEP(201).EQ.1) THEN
            IF (NextPiv2beWritten.LE.NPIV) THEN
              LAST_CALL=.FALSE.
              CALL DMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE, TYPEFile,
     &        A(POSELT), LAFAC, MonBloc,
     &        NextPiv2beWritten, IDUMMY,
     &        IW(IOLDPS), LIWFAC, MYID,
     &        KEEP8(31),
     &        IFLAG,LAST_CALL )
              IF (IFLAG .LT. 0 ) RETURN
            ENDIF
          ENDIF
        ENDDO
        IF ( NFRONT - NASS - IROW + 1 - Block > 0 ) THEN
        CALL dgemm( 'N', 'N', Block, NFRONT-NASS-Block-IROW+1, NPIV,
     &              ALPHA,  A( UPOS ), LDA,
     &              A( LPOS + LDA8 * int(Block,8) ), LDA,
     &              BETA,
     &              A( APOS + LDA8 * int(Block,8) ), LDA )
        ENDIF
       END DO
#if defined(GEMMT_AVAILABLE)
      END IF
#endif
       IF ( (POSTPONE_COL_UPDATE).AND.(NASS-NPIV.GT.0) ) THEN
         LPOS = POSELT + int(NPIV,8)*LDA8
         UPOS = POSELT + int(NPIV,8)
         CALL DMUMPS_FAC_LDLT_COPYSCALE_U( NASS-NPIV, 1, 
     &       KEEP(424), NFRONT, NPIV, 
     &       LIW, IW, OFFSET_IW, LA, A, POSELT, LPOS, UPOS, POSELT)
         LPOS = POSELT + LDA8 * int(NASS,8)
         CALL dgemm('N', 'N', NASS-NPIV, NFRONT-NASS, NPIV, ALPHA,
     &               A( POSELT + int(NPIV,8)), LDA,
     &               A( LPOS ), LDA,
     &               BETA,
     &               A( LPOS + int(NPIV,8) ), LDA)
       ENDIF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_FAC_T_LDLT
      SUBROUTINE DMUMPS_STORE_PERMINFO( PIVRPTR, NBPANELS, PIVR, NASS,
     &                                  K, P, LastPanelonDisk,
     &                                  LastPIVRPTRIndexFilled )
      IMPLICIT NONE
      INTEGER, intent(in) :: NBPANELS, NASS, K, P
      INTEGER, intent(inout) :: PIVRPTR(NBPANELS), PIVR(NASS)
      INTEGER LastPanelonDisk, LastPIVRPTRIndexFilled
      INTEGER I
      IF ( LastPanelonDisk+1 > NBPANELS ) THEN
           WRITE(*,*) "INTERNAL ERROR IN DMUMPS_STORE_PERMINFO!"
           WRITE(*,*) "NASS=",NASS,"PIVRPTR=",PIVRPTR(1:NBPANELS)
           WRITE(*,*) "K=",K, "P=",P, "LastPanelonDisk=",LastPanelonDisk
           WRITE(*,*) "LastPIVRPTRIndexFilled=", LastPIVRPTRIndexFilled
           CALL MUMPS_ABORT()
      ENDIF
      PIVRPTR(LastPanelonDisk+1) = K + 1
      IF (LastPanelonDisk.NE.0) THEN
        PIVR(K - PIVRPTR(1) + 1) = P
        DO I = LastPIVRPTRIndexFilled + 1, LastPanelonDisk
          PIVRPTR(I)=PIVRPTR(LastPIVRPTRIndexFilled)
        ENDDO
      ENDIF
      LastPIVRPTRIndexFilled = LastPanelonDisk + 1
      RETURN
      END SUBROUTINE DMUMPS_STORE_PERMINFO
      SUBROUTINE DMUMPS_UPDATE_MINMAX_PIVOT 
     &           ( DIAG, DKEEP, KEEP, NULLPIVOT)
!$    USE OMP_LIB
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)    :: DIAG
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      LOGICAL, INTENT(IN) :: NULLPIVOT
      INTEGER, INTENT(IN) :: KEEP(500)
      IF (KEEP(405).EQ.0) THEN
        DKEEP(21) = max(DKEEP(21), DIAG)
        DKEEP(19) = min(DKEEP(19), DIAG)
        IF (.NOT.NULLPIVOT) THEN
          DKEEP(20) = min(DKEEP(20), DIAG)
        ENDIF
      ELSE
!$OMP   ATOMIC UPDATE
        DKEEP(21) = max(DKEEP(21), DIAG)
!$OMP   END ATOMIC
!$OMP   ATOMIC UPDATE
        DKEEP(19) = min(DKEEP(19), DIAG)
!$OMP   END ATOMIC
        IF (.NOT.NULLPIVOT) THEN
!$OMP     ATOMIC UPDATE
          DKEEP(20) = min(DKEEP(20), DIAG)
!$OMP     END ATOMIC
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_UPDATE_MINMAX_PIVOT 
      SUBROUTINE DMUMPS_COMPUTE_SIZE_SCHUR_IN_FRONT (
     &           N, NCB, SIZE_SCHUR, ROW_INDICES, PERM, 
     &           NVSCHUR
     &           )
      IMPLICIT NONE
      INTEGER, intent(in) :: N, NCB, SIZE_SCHUR
      INTEGER, intent(in) :: ROW_INDICES(NCB), PERM(N)
      INTEGER, intent(out):: NVSCHUR
      INTEGER :: I, IPOS, IBEG_SCHUR
      IBEG_SCHUR = N - SIZE_SCHUR +1
      NVSCHUR    = 0
      IPOS       = NCB
      DO I= NCB,1,-1
       IF (abs(ROW_INDICES(I)).LE.N) THEN
        IF (PERM(ROW_INDICES(I)).LT.IBEG_SCHUR) EXIT 
       ENDIF
       IPOS = IPOS -1
      ENDDO
      NVSCHUR = NCB-IPOS
      RETURN
      END SUBROUTINE DMUMPS_COMPUTE_SIZE_SCHUR_IN_FRONT
      END MODULE DMUMPS_FAC_FRONT_AUX_M
