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
      RECURSIVE SUBROUTINE DMUMPS_TRAITER_MESSAGE_SOLVE
     &     ( BUFR, LBUFR, LBUFR_BYTES,
     &     MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &     N, NRHS, IPOOL, LPOOL, LEAF,
     &     NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     &     PTRFAC, IWCB, LIWCB,
     &     WCB, LWCB, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB,
     &     INFO, KEEP, KEEP8, DKEEP, STEP, PROCNODE_STEPS, 
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &     , FROM_PP
     &    )
      USE DMUMPS_OOC 
      USE DMUMPS_SOL_LR, ONLY: DMUMPS_SOL_SLAVE_LR_U
      USE DMUMPS_BUF 
      IMPLICIT NONE
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIW
      INTEGER(8), INTENT(IN) :: LA, LWCB
      INTEGER N, NRHS, LPOOL, LEAF, NBFIN, LRHSCOMP
      INTEGER LIWCB, POSIWCB
      INTEGER(8) :: POSWCB, PLEFTWCB
      INTEGER INFO( 80 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ),  NSTK_S( N )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      INTEGER PTRICB(KEEP(28)),PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      DOUBLE PRECISION WCB( LWCB ), A( LA )
      DOUBLE PRECISION RHSCOMP( LRHSCOMP, NRHS )
      INTEGER, intent(in) :: POSINRHSCOMP_FWD(N)
      LOGICAL, intent(in) :: FROM_PP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER(8) :: PTRX, PTRY, IFR8
      INTEGER IERR, K, JJ, JBDEB, JBFIN, NRHS_B
      INTEGER :: IWHDLR, LDA_SLAVE
      INTEGER :: MTYPE_SLAVE 
      INTEGER FINODE, FPERE, LONG, NCB, POSITION, NCV, NPIV
      INTEGER PDEST, I, IPOSINRHSCOMP
      INTEGER J1
      INTEGER(8) :: APOS
      LOGICAL DUMMY
      LOGICAL FLAG
!$    LOGICAL :: OMP_FLAG
      EXTERNAL MUMPS_PROCNODE
      INTEGER  MUMPS_PROCNODE
      LOGICAL COMPRESS_PANEL, LR_ACTIVATED
      LOGICAL OOCWRITE_COMPATIBLE_WITH_BLR 
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      INCLUDE 'mumps_headers.h'
      IF ( MSGTAG .EQ. RACINE_SOLVE ) THEN
         NBFIN = NBFIN - 1
         IF ( NBFIN .eq. 0 ) GOTO 270
      ELSE  IF (MSGTAG .EQ. ContVec ) THEN
         POSITION = 0
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        FINODE, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        FPERE, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        NCB, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        JBDEB, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        JBFIN, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        LONG, 1, MPI_INTEGER, COMM, IERR )
         NRHS_B = JBFIN-JBDEB+1
          IF ( NCB .eq. 0 ) THEN
             PTRICB(STEP(FINODE)) = -1
          ELSE
            IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
              PTRICB(STEP(FINODE)) = NCB + 1
            END IF
            IF ( ( POSIWCB - LONG ) .LT. 0 ) THEN
               INFO( 1 ) = -14
               INFO( 2 ) = LONG
               GOTO 260
            END IF
            IF ( POSWCB - PLEFTWCB + 1_8 .LT.
     &           int(LONG,8) * int(NRHS_B,8)) THEN
               INFO( 1 ) = -11
               CALL MUMPS_SET_IERROR(PLEFTWCB-POSWCB-1_8+
     &             int(LONG,8) * int(NRHS_B,8),
     &             INFO(2))
               GOTO 260
            END IF
            IF (LONG .GT. 0) THEN
              CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &             IWCB( 1 ),
     &             LONG, MPI_INTEGER, COMM, IERR )
              DO K = 1, NRHS_B
                 CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                WCB( PLEFTWCB ),
     &                LONG, MPI_DOUBLE_PRECISION, COMM, IERR )
                 DO I = 1, LONG
                 IPOSINRHSCOMP= abs(POSINRHSCOMP_FWD(IWCB(I)))
                 RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1) = 
     &           RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1) + 
     &               WCB(PLEFTWCB+I-1)
                 ENDDO
              END DO
              PTRICB(STEP(FINODE)) = PTRICB(STEP(FINODE)) - LONG
            ENDIF
          END IF  
          IF ( PTRICB(STEP(FINODE)) == 1 .OR.
     &         PTRICB(STEP(FINODE)) == -1 ) THEN
             NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             PTRICB(STEP(FINODE)) = 0
          END IF
          IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
             IPOOL( LEAF ) = FPERE
             LEAF = LEAF + 1
             IF ( LEAF > LPOOL ) THEN
                WRITE(*,*)
     &          'Internal error 1 DMUMPS_TRAITER_MESSAGE_SOLVE',
     &          LEAF, LPOOL
                CALL MUMPS_ABORT()
             END IF
          ENDIF
      ELSEIF ( MSGTAG .EQ. Master2Slave ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         NCV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         NPIV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         JBDEB, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         JBFIN, 1, MPI_INTEGER, COMM, IERR )
          NRHS_B = JBFIN-JBDEB+1
          PTRY = PLEFTWCB
          PTRX = PLEFTWCB + int(NCV,8) * int(NRHS_B,8)
          PLEFTWCB = PLEFTWCB + int(NPIV + NCV,8) * int(NRHS_B,8)
          IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
             INFO(1) = -11
             CALL MUMPS_SET_IERROR(-POSWCB+PLEFTWCB-1_8,INFO(2))
             GO TO 260
          END IF
          DO K=1, NRHS_B
             CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &            WCB( PTRY + (K-1) * NCV ), NCV,
     &            MPI_DOUBLE_PRECISION, COMM, IERR )
          ENDDO
          IF ( NPIV .GT. 0 ) THEN
             DO K=1, NRHS_B
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               WCB( PTRX + (K-1)*NPIV ), NPIV,
     &               MPI_DOUBLE_PRECISION, COMM, IERR )
             END DO
          END IF
          LR_ACTIVATED = (IW(PTRIST(STEP(FINODE))+XXLR).GT.0)
          COMPRESS_PANEL = (IW(PTRIST(STEP(FINODE))+XXLR).GE.2)
          OOCWRITE_COMPATIBLE_WITH_BLR =
     &     (.NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(485).EQ.0) 
     &     )
          IF (KEEP(201).GT.0.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
             CALL DMUMPS_SOLVE_GET_OOC_NODE(
     &            FINODE,PTRFAC,KEEP,A,LA,STEP,
     &            KEEP8,N,DUMMY,IERR)
              IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          IF ( IW(PTRIST(STEP(FINODE))+XXLR) .GE. 2 .AND.
     &         KEEP(485) .EQ. 1 ) THEN
            IWHDLR = IW(PTRIST(STEP(FINODE))+XXF)
            MTYPE_SLAVE = 1     
            CALL DMUMPS_SOL_SLAVE_LR_U( FINODE, IWHDLR,
     &              -9999,  
     &              WCB, LWCB,
     &              NPIV, NCV, 
     &              PTRX, PTRY,
     &              JBDEB, JBFIN,
     &              MTYPE_SLAVE, KEEP, 
     &              INFO(1), INFO(2) )
          ELSE
            APOS = PTRFAC(STEP(FINODE))
            IF (KEEP(201) .EQ. 1) THEN
              MTYPE_SLAVE = 0
              LDA_SLAVE = NCV
            ELSE
              MTYPE_SLAVE = 1     
              LDA_SLAVE = NPIV
            ENDIF
            CALL DMUMPS_SOLVE_GEMM_UPDATE
     &        ( A, LA, APOS, NPIV,
     &          LDA_SLAVE, 
     &          NCV,
     &          NRHS_B, WCB, LWCB,
     &          PTRX, NPIV,  
     &          PTRY, NCV,   
     &          MTYPE_SLAVE, KEEP, ONE )
          ENDIF
          IF ((KEEP(201).GT.0).AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
             CALL DMUMPS_FREE_FACTORS_FOR_SOLVE(FINODE,PTRFAC,
     &            KEEP(28),A,LA,.TRUE.,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          PLEFTWCB = PLEFTWCB - int(NPIV,8) * int(NRHS_B,8)
          PDEST = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(FPERE)),
     &                            KEEP(199) )
          IF ( PDEST .EQ. MYID ) THEN
             IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
                NCB = IW( PTRIST(STEP(FINODE)) + 2 + KEEP(IXSZ) )
                PTRICB(STEP(FINODE)) = NCB + 1
             END IF
             J1 = PTRIST(STEP(FINODE))+3+KEEP(IXSZ) 
!$           OMP_FLAG = ( JBFIN-JBDEB+1.GE.KEEP(362) .AND.
!$   &                     (NCV*(JBFIN-JBDEB+1) .GE. KEEP(363) ) )
!$OMP PARALLEL DO PRIVATE(I,JJ,IFR8,IPOSINRHSCOMP) IF(OMP_FLAG)
             DO K=1, NRHS_B
                  IFR8 = PTRY+int(K-1,8)*int(NCV,8)
                  DO I = 1,NCV
                  JJ = IW(J1+I)
                  IPOSINRHSCOMP= abs(POSINRHSCOMP_FWD(JJ))
                  RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1)= 
     &              RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1)
     &              + WCB(IFR8+int(I-1,8))
                  ENDDO
             ENDDO
!$OMP END PARALLEL DO
            PTRICB(STEP(FINODE)) = PTRICB(STEP(FINODE)) - NCV
             IF ( PTRICB( STEP( FINODE ) ) == 1 ) THEN
                NSTK_S(STEP(FPERE))  = NSTK_S(STEP(FPERE)) - 1
                PTRICB(STEP(FINODE)) = 0
             END IF
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
               IPOOL( LEAF ) = FPERE
               LEAF = LEAF + 1
               IF ( LEAF > LPOOL ) THEN
                 WRITE(*,*)
     &           'INTERNAL Error in DMUMPS_TRAITER_MESSAGE_SOLVE',
     &           LEAF, LPOOL
                 CALL MUMPS_ABORT()
               END IF
            ENDIF
          ELSE
 210        CONTINUE
            CALL DMUMPS_BUF_SEND_VCB( NRHS_B, FINODE, FPERE,
     &            IW(PTRIST(STEP( FINODE )) + 2 + KEEP(IXSZ) ), NCV,NCV,
     &            IW(PTRIST(STEP(FINODE))+4+ KEEP(IXSZ) ),
     &            WCB( PTRY ), JBDEB, JBFIN,
     &            RHSCOMP, 1, 1, -9999, -9999,
     &            KEEP, PDEST, ContVec, COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
              CALL DMUMPS_SOLVE_RECV_AND_TREAT( .FALSE., FLAG,
     &              BUFR, LBUFR, LBUFR_BYTES,
     &              MYID, SLAVEF, COMM,
     &              N, NRHS, IPOOL, LPOOL, LEAF,
     &              NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &              IWCB, LIWCB,
     &              WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &              PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP,
     &              PROCNODE_STEPS, 
     &              RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &              , FROM_PP
     &              )
              IF ( INFO( 1 )  .LT. 0 )  GOTO 270
              GOTO 210
            ELSE IF ( IERR .EQ. -2 ) THEN
              INFO( 1 ) = -17
              INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     &                    NCV * KEEP( 35 )
              GOTO 260
            ELSE IF ( IERR .EQ. -3 ) THEN
              INFO( 1 ) = -20
              INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     &               NCV * KEEP( 35 )
            END IF
         END IF
         PLEFTWCB = PLEFTWCB - int(NCV,8) * int(NRHS_B,8)
      ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GOTO 270
      ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     &         (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
      ELSE
          INFO(1)=-100
          INFO(2)=MSGTAG
          GO TO 260
      ENDIF
      GO TO 270
 260  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
 270  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_TRAITER_MESSAGE_SOLVE
      SUBROUTINE DMUMPS_SOLVE_NODE_FWD( INODE,
     &     LASTFSL0STA, LASTFSL0DYN, 
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MYID, SLAVEF, COMM,
     &     N, IPOOL, LPOOL, LEAF,
     &     NBFIN, NSTK_S,
     &     IWCB, LIWCB,
     &     WCB, LWCB, A, LA, IW, LIW,
     &     NRHS, POSWCB, PLEFTWCB, POSIWCB,
     &     PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     &     FILS, STEP, FRERE, DAD,
     &     INFO, KEEP,KEEP8, DKEEP, RHS_ROOT, LRHS_ROOT, MTYPE,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD,
     &     
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP
     &     , ERROR_WAS_BROADCASTED
     &    )
      USE DMUMPS_SOL_LR 
      USE DMUMPS_OOC
      USE DMUMPS_BUF
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER, INTENT( IN ) :: INODE, LASTFSL0STA, LASTFSL0DYN
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MYID, SLAVEF, COMM
      INTEGER LIWCB, LIW, POSIWCB
      INTEGER(8) :: POSWCB, PLEFTWCB, LWCB
      INTEGER(8) :: LA
      INTEGER N, LPOOL, LEAF, NBFIN
      INTEGER INFO( 80 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ), NSTK_S(KEEP(28))
      INTEGER IWCB( LIWCB ), IW( LIW )
      INTEGER NRHS
      DOUBLE PRECISION WCB( LWCB ), A( LA )
      INTEGER(8) :: LRHS_ROOT
      DOUBLE PRECISION RHS_ROOT( LRHS_ROOT )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER FILS( N ), STEP( N ), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &     TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER POSINRHSCOMP_FWD(N), LRHSCOMP
      DOUBLE PRECISION RHSCOMP(LRHSCOMP, NRHS)
      LOGICAL, intent(in) :: DO_NBSPARSE
      INTEGER, intent(in) :: LRHS_BOUNDS
      INTEGER, intent(in) :: RHS_BOUNDS(LRHS_BOUNDS)
      LOGICAL, intent(in) :: FROM_PP
      LOGICAL, intent(out) :: ERROR_WAS_BROADCASTED
      EXTERNAL dgemv, dtrsv, dgemm, dtrsm, MUMPS_PROCNODE
      INTEGER MUMPS_PROCNODE
      DOUBLE PRECISION ALPHA,ONE,ZERO
      PARAMETER (ZERO=0.0D0, ONE = 1.0D0, ALPHA=-1.0D0)
      INTEGER :: IWHDLR
      INTEGER JBDEB, JBFIN, NRHS_B 
      INTEGER LDADIAG
      INTEGER(8) :: APOS, APOS1, IFR8, IFR_ini8
      INTEGER I, J, K, IPOS, NSLAVES, J1, J2, J3, FPERE, FPERE_MAPPING,
     &     NPIV, NCB, LIELL, JJ, NELIM, IERR
      INTEGER(8) :: PCB_COURANT, PPIV_COURANT, PPIV_PANEL, PCB_PANEL
      INTEGER IPOSINRHSCOMP_TMP
      INTEGER Effective_CB_Size, NUPDATE, ISLAVE, PDEST, FirstIndex
      LOGICAL FLAG
      INTEGER :: NUPDATE_NONCRITICAL, IPOSINRHSCOMPLASTFSDYN 
!$    LOGICAL :: OMP_FLAG
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: APOSDEB
      INTEGER TempNROW, TempNCOL, PANEL_SIZE, 
     &     JFIN, NBJ, NUPDATE_PANEL,
     &     TYPEF
      INTEGER LD_WCBPIV         
      INTEGER LD_WCBCB          
      LOGICAL :: LDEQLIELLPANEL 
      LOGICAL :: CBINITZERO  
      INTEGER LDAJ, LDAJ_FIRST_PANEL
      INTEGER LDAtemp
      LOGICAL COMPRESS_PANEL, LR_ACTIVATED
      LOGICAL OOCWRITE_COMPATIBLE_WITH_BLR 
      INTEGER TMP_NBPANELS,
     &     I_PIVRPTR, I_PIVR, IPANEL
      LOGICAL MUST_BE_PERMUTED
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY( 1 )
      ERROR_WAS_BROADCASTED = .FALSE.
      DUMMY(1)=1
      LR_ACTIVATED   = (IW(PTRIST(STEP(INODE))+XXLR).GT.0)
      COMPRESS_PANEL = (IW(PTRIST(STEP(INODE))+XXLR).GE.2)
      OOCWRITE_COMPATIBLE_WITH_BLR =
     &     (.NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(485).EQ.0) 
     &     )
      IF (DO_NBSPARSE) THEN
       JBDEB= RHS_BOUNDS(2*STEP(INODE)-1)
       JBFIN= RHS_BOUNDS(2*STEP(INODE))
      ELSE
       JBDEB = 1
       JBFIN = NRHS
      ENDIF
      NRHS_B = JBFIN-JBDEB+1
      IF (DO_NBSPARSE) THEN
         if (JBDEB.GT.JBFIN) then
         write(6,*) " Internal error 1 in nbsparse :", 
     &    JBDEB, JBFIN
         CALL MUMPS_ABORT()
       endif
       IF (JBDEB.LT.1 .OR. JBDEB.GT.NRHS .or. 
     &      JBFIN.LT.1 .OR. JBFIN.GT.NRHS ) THEN
         write(6,*) " Internal error 2 in nbsparse :", 
     &    JBDEB, JBFIN
         CALL MUMPS_ABORT()
       endif
      ENDIF
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq.KEEP( 20 ) ) THEN
         LIELL = IW( PTRIST( STEP(INODE)) + 3 + KEEP(IXSZ))
         NPIV  = LIELL
         NELIM = 0
         NSLAVES = 0
         IPOS = PTRIST( STEP(INODE)) + 5 + KEEP(IXSZ)
      ELSE
        IPOS = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        NELIM = IW(IPOS-1)
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ) )
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        IF ((KEEP(201).GT.0).AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
           CALL DMUMPS_SOLVE_GET_OOC_NODE(
     &          INODE,PTRFAC,KEEP,A,LA,STEP,
     &          KEEP8,N,MUST_BE_PERMUTED,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              ERROR_WAS_BROADCASTED = .FALSE.
              GOTO 270
           ENDIF
           IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
           CALL DMUMPS_OOC_PP_CHECK_PERM_FREED(
     &                 IW(IPOS+1+2*LIELL+1+NSLAVES),
     &                 MUST_BE_PERMUTED )
           ENDIF
        ENDIF                     
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ))
        IPOS = IPOS + 1 + NSLAVES
      END IF
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
         J1 = IPOS + 1
         J2 = IPOS + LIELL
         J3 = IPOS + NPIV
      ELSE
         J1 = IPOS + LIELL + 1
         J2 = IPOS + 2 * LIELL
         J3 = IPOS + LIELL + NPIV
      END IF
      NCB = LIELL-NPIV
      IF (KEEP(50).NE.0) THEN
       LDADIAG = NPIV
      ELSE
       LDADIAG = LIELL
      ENDIF
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq. KEEP(20) ) THEN
        IFR8 = 0_8
        IPOSINRHSCOMP_TMP = POSINRHSCOMP_FWD(IW(J1)) 
        IFR_ini8 = IFR8
!$      OMP_FLAG = ( JBFIN-JBDEB+1.GE.KEEP(362) .AND.
!$   &                     (J3-J1+1)*(JBFIN-JBDEB+1) .GE. KEEP(363) )
!$OMP PARALLEL DO PRIVATE(IFR8,JJ) IF(OMP_FLAG)
        DO K=JBDEB,JBFIN
           IFR8 = IFR_ini8 + int(K-1,8)*int(NPIV,8)
           DO JJ = J1, J3
             RHS_ROOT(IFR8+int(JJ-J1+1,8)) = 
     &               RHSCOMP(IPOSINRHSCOMP_TMP+JJ-J1,K)
           ENDDO
        ENDDO
!$OMP END PARALLEL DO
        IF ( NPIV .LT. LIELL ) THEN
            WRITE(*,*) ' Internal error 1 in DMUMPS_SOLVE_NODE_FWD',
     &      NPIV, LIELL
            CALL MUMPS_ABORT()
        END IF
        GO TO 270
      END IF
      APOS = PTRFAC(STEP(INODE))
      IF ( (KEEP(201).EQ.1).AND.OOCWRITE_COMPATIBLE_WITH_BLR ) THEN  
        IF (MTYPE.EQ.1) THEN
            IF ((MTYPE.EQ.1).AND.NSLAVES.NE.0) THEN
              TempNROW= NPIV+NELIM
              TempNCOL= NPIV
              LDAJ_FIRST_PANEL=TempNROW
            ELSE
              TempNROW= LIELL
              TempNCOL= NPIV
              LDAJ_FIRST_PANEL=TempNROW
            ENDIF
            TYPEF=TYPEF_L
        ELSE 
            TempNCOL= LIELL
            TempNROW= NPIV
            LDAJ_FIRST_PANEL=TempNCOL
            TYPEF= TYPEF_U
        ENDIF
        PANEL_SIZE = DMUMPS_OOC_PANEL_SIZE( LDAJ_FIRST_PANEL )
      ENDIF                     
      PPIV_COURANT = PLEFTWCB
      PLEFTWCB = PLEFTWCB + int(LIELL,8) * int(NRHS_B,8)
      IF ( POSWCB - PLEFTWCB + 1_8 .LT. 0 ) THEN
         INFO(1) = -11
         CALL MUMPS_SET_IERROR(PLEFTWCB-POSWCB-1_8, INFO(2))
         ERROR_WAS_BROADCASTED = .FALSE.
         GOTO 270
      END IF
      IF (KEEP(201) .EQ. 1 .AND. OOCWRITE_COMPATIBLE_WITH_BLR) THEN
        LDEQLIELLPANEL = .TRUE.
        LD_WCBPIV      = LIELL
        LD_WCBCB       = LIELL
        PCB_COURANT    = PPIV_COURANT + NPIV
      ELSE
        LDEQLIELLPANEL = .FALSE.
        LD_WCBPIV   = NPIV
        LD_WCBCB    = NCB
        PCB_COURANT = PPIV_COURANT + int(NPIV,8)*int(NRHS_B,8)
      ENDIF
      FPERE = DAD(STEP(INODE))
      IF ( FPERE .NE. 0 ) THEN
        FPERE_MAPPING = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(FPERE)),
     &  KEEP(199) )
      ELSE
        FPERE_MAPPING = -1
      ENDIF
      IF ( LASTFSL0DYN .LE. N ) THEN
        CBINITZERO = .TRUE.  
      ELSE IF ( FPERE_MAPPING .EQ. MYID ) THEN
        CBINITZERO = .TRUE.  
      ELSE
        CBINITZERO = .FALSE. 
      ENDIF
      CALL DMUMPS_RHSCOMP_TO_WCB(
     &     NPIV, NCB, LIELL, CBINITZERO, LDEQLIELLPANEL,
     &     RHSCOMP(1, JBDEB), LRHSCOMP, NRHS_B,
     &     POSINRHSCOMP_FWD, N,
     &     WCB(PPIV_COURANT),
     &     IW, LIW, J1, J3, J2, KEEP, DKEEP)
      IF ( NPIV .NE. 0 ) THEN
         IF ((KEEP(201).EQ.1).AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN 
        APOSDEB = APOS
        J = 1
        IPANEL = 0
  10    CONTINUE
          IPANEL = IPANEL + 1
          JFIN    = min(J+PANEL_SIZE-1, NPIV)
          IF (IW(IPOS+ LIELL + JFIN) < 0) THEN
            JFIN=JFIN+1
          ENDIF
          NBJ     = JFIN-J+1
          LDAJ    = LDAJ_FIRST_PANEL-J+1 
          IF ( (KEEP(50).NE.1).AND. MUST_BE_PERMUTED ) THEN
           CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF, TMP_NBPANELS,
     &            I_PIVRPTR, I_PIVR, IPOS+1+2*LIELL, IW, LIW)
               IF (NPIV.EQ.(IW(I_PIVRPTR+IPANEL-1)-1)) THEN
                  MUST_BE_PERMUTED=.FALSE. 
               ELSE
                  CALL DMUMPS_PERMUTE_PANEL(
     &                 IW( I_PIVR+ IW(I_PIVRPTR+IPANEL-1)-
     &                 IW(I_PIVRPTR)),
     &                 NPIV-IW(I_PIVRPTR+IPANEL-1)+1, 
     &                 IW(I_PIVRPTR+IPANEL-1)-1,
     &                 A(APOSDEB),
     &                 LDAJ, NBJ, J-1 ) 
               ENDIF
            ENDIF 
            NUPDATE_PANEL = LDAJ - NBJ
            PPIV_PANEL = PPIV_COURANT+int(J-1,8)
            PCB_PANEL  = PPIV_PANEL+int(NBJ,8)
            APOS1 = APOSDEB+int(NBJ,8)
            IF  (MTYPE.EQ.1) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dtrsv( 'L', 'N', 'U', NBJ, A(APOSDEB), LDAJ, 
     &                 WCB(PPIV_PANEL), 1 )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL dgemv('N', NUPDATE_PANEL,NBJ,ALPHA, A(APOS1),
     &                    LDAJ,  WCB(PPIV_PANEL), 1, ONE,
     &                    WCB(PCB_PANEL), 1)
                  ENDIF
               ELSE
#endif
                  CALL dtrsm( 'L','L','N','U', NBJ, NRHS_B, ONE,
     &                 A(APOSDEB), LDAJ, WCB(PPIV_PANEL),
     &                 LIELL )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL dgemm('N', 'N', NUPDATE_PANEL, NRHS_B, NBJ, 
     &                    ALPHA,
     &                    A(APOS1), LDAJ, WCB(PPIV_PANEL), LIELL, ONE,
     &                    WCB(PCB_PANEL), LIELL)
                  ENDIF
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
            ELSE
#if defined(MUMPS_USE_BLAS2)
               IF (NRHS_B == 1) THEN
                  CALL dtrsv( 'L', 'N', 'N', NBJ, A(APOSDEB), LDAJ,
     &                 WCB(PPIV_PANEL), 1 )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL dgemv('N',NUPDATE_PANEL, NBJ, ALPHA, A(APOS1),
     &                    LDAJ, WCB(PPIV_PANEL), 1,
     &                    ONE, WCB(PCB_PANEL), 1 )
                  ENDIF
               ELSE
#endif
                  CALL dtrsm('L','L','N','N',NBJ, NRHS_B, ONE,
     &                 A(APOSDEB), LDAJ, WCB(PPIV_PANEL),
     &                 LIELL)
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL dgemm('N', 'N', NUPDATE_PANEL, NRHS_B, NBJ, 
     &                    ALPHA,
     &                    A(APOS1), LDAJ, WCB(PPIV_PANEL), LIELL, ONE,
     &             WCB(PCB_PANEL), LIELL)
                  ENDIF
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
            ENDIF
            APOSDEB = APOSDEB+int(LDAJ,8)*int(NBJ,8)
            J=JFIN+1
            IF ( J .LE. NPIV ) GOTO 10
         ELSE                   
          IF ( IW(PTRIST(STEP(INODE))+XXLR) .GE. 2 .AND.
     &         KEEP(485) .EQ. 1 ) THEN
            IWHDLR = IW(PTRIST(STEP(INODE))+XXF)
            CALL DMUMPS_SOL_FWD_LR_SU (
     &           INODE, N, IWHDLR, NPIV, NSLAVES, 
     &           IW, IPOS, LIW,
     &           LIELL, WCB, LWCB, 
     &           LD_WCBPIV, LD_WCBCB,
     &           PPIV_COURANT, PCB_COURANT, 
     &           RHSCOMP, LRHSCOMP, NRHS, 
     &           POSINRHSCOMP_FWD, JBDEB, JBFIN,
     &           MTYPE, KEEP, OOCWRITE_COMPATIBLE_WITH_BLR,
     &           INFO(1), INFO(2) )
            IF (INFO(1).LT.0) THEN
              ERROR_WAS_BROADCASTED = .FALSE.
              GOTO 270
            ENDIF
          ELSE
           CALL DMUMPS_SOLVE_FWD_TRSOLVE (
     &         A, LA, APOS, NPIV, LDADIAG, 
     &         NRHS_B, WCB, LWCB, LD_WCBPIV,
     &         PPIV_COURANT, MTYPE, KEEP)
          ENDIF
         END IF                 
      END IF                    
      NCB   = LIELL - NPIV
      IF ( MTYPE .EQ. 1 ) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            APOS1 = APOS  + int(NPIV,8) * int(LIELL,8)
         ELSE
            APOS1 = APOS + int(NPIV,8) * int(NPIV,8)
         END IF
         IF ( NSLAVES .EQ. 0 .OR. NPIV .eq. 0 ) THEN
            NUPDATE = NCB
         ELSE
            NUPDATE = NELIM
         END IF
      ELSE
         APOS1 = APOS + int(NPIV,8)
         NUPDATE = NCB
      END IF
      IF (KEEP(201).NE.1) THEN  
      IF ( IW(PTRIST(STEP(INODE))+XXLR) .LT. 2 .OR.
     &         KEEP(485).EQ.0) THEN
            IF (MTYPE .EQ. 1) THEN
              LDAtemp = NPIV
            ELSE
              LDAtemp = LIELL
            ENDIF
            CALL DMUMPS_SOLVE_GEMM_UPDATE 
     &           (A, LA, APOS1, NPIV, LDAtemp, NUPDATE, 
     &           NRHS_B, WCB, LWCB, PPIV_COURANT, LD_WCBPIV,
     &           PCB_COURANT, LD_WCBCB,
     &           MTYPE, KEEP, ONE)
      ENDIF
      END IF                    
      IF ( IW(PTRIST(STEP(INODE))+XXLR) .LT. 2 .OR.
     &         KEEP(485).EQ.0) THEN
        CALL DMUMPS_SOLVE_LD_AND_RELOAD (
     &    INODE, N, NPIV, LIELL, NELIM, NSLAVES,
     &    PPIV_COURANT, 
     &    IW, IPOS, LIW, 
     &    A, LA, APOS,
     &    WCB, LWCB, LD_WCBPIV, 
     &    RHSCOMP, LRHSCOMP, NRHS, 
     &    POSINRHSCOMP_FWD, JBDEB, JBFIN, 
     &    MTYPE, KEEP, OOCWRITE_COMPATIBLE_WITH_BLR
     &  )
      ENDIF
      IF ((KEEP(201).EQ.1).AND.OOCWRITE_COMPATIBLE_WITH_BLR) 
     &THEN 
         CALL DMUMPS_FREE_FACTORS_FOR_SOLVE(INODE,PTRFAC,KEEP(28),
     &        A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            ERROR_WAS_BROADCASTED = .FALSE.
            GOTO 270
         ENDIF
      END IF
      IF ( FPERE .EQ. 0 ) THEN
        PLEFTWCB = PLEFTWCB - int(LIELL,8) *int(NRHS_B,8)
        GOTO 270
      ENDIF
      IF ( NUPDATE .NE. 0 .OR. NCB.EQ.0 ) THEN
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(FPERE)),
     &        KEEP(199)) .EQ. MYID) THEN
          IF ( NCB .ne. 0 ) THEN
            PTRICB(STEP(INODE)) = NCB + 1
            NUPDATE_NONCRITICAL = NUPDATE
            IF (LASTFSL0DYN .LE. N) THEN
              IF ( LASTFSL0DYN .EQ. 0 ) THEN
                IPOSINRHSCOMPLASTFSDYN = 0
              ELSE
                IPOSINRHSCOMPLASTFSDYN =
     &                          abs(POSINRHSCOMP_FWD(LASTFSL0DYN))
              ENDIF
              DO I = 1, NUPDATE
                IF ( abs(POSINRHSCOMP_FWD( IW(J3+I) )) .GT. 
     &               IPOSINRHSCOMPLASTFSDYN ) THEN
                  IF (abs(STEP(IW(J3+I))) .GT.
     &                abs(STEP( LASTFSL0STA))
     &                .OR. KEEP(261) .NE. 1) THEN
                    NUPDATE_NONCRITICAL = I - 1
                    EXIT
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
!$          OMP_FLAG = ( NRHS_B.GE.KEEP(362) .AND.
!$   &                  (NUPDATE*NRHS_B .GE. KEEP(363)) )
!$OMP PARALLEL DO PRIVATE(I,IFR8,IPOSINRHSCOMP_TMP) IF(OMP_FLAG)
            DO K = JBDEB, JBFIN
              IFR8 = PCB_COURANT + int(K-JBDEB,8)*int(LD_WCBCB,8)
              DO I = 1, NUPDATE_NONCRITICAL
                IPOSINRHSCOMP_TMP = 
     &          abs(POSINRHSCOMP_FWD(IW(J3 + I)))
                RHSCOMP( IPOSINRHSCOMP_TMP, K ) = 
     &          RHSCOMP( IPOSINRHSCOMP_TMP, K )
     &          + WCB(IFR8 + int(I-1,8))
              ENDDO
            ENDDO
!$OMP END PARALLEL DO
            IF ( CBINITZERO ) THEN
              IF ( NUPDATE .NE. NUPDATE_NONCRITICAL) THEN
                IF (.NOT. CBINITZERO) THEN
                WRITE(*,*) ' Internal error 3 in DMUMPS_SOLVE_NODE_FWD',
     &          CBINITZERO, INODE, NUPDATE, NUPDATE_NONCRITICAL
                CALL MUMPS_ABORT()
                ENDIF
                DO K = JBDEB, JBFIN
                  IFR8 = PCB_COURANT + int(K-JBDEB,8)*int(LD_WCBCB,8)
!$OMP CRITICAL(DMUMPS_RHSCOMP_CRI)
                  DO I = NUPDATE_NONCRITICAL+1, NUPDATE
                    IPOSINRHSCOMP_TMP = 
     &                abs(POSINRHSCOMP_FWD(IW(J3 + I)))
                    RHSCOMP( IPOSINRHSCOMP_TMP, K ) = 
     &                   RHSCOMP( IPOSINRHSCOMP_TMP, K )
     &                   + WCB(IFR8 + int(I-1,8))
                  ENDDO
!$OMP END CRITICAL(DMUMPS_RHSCOMP_CRI)
                ENDDO
              ENDIF
            ENDIF
            PTRICB(STEP( INODE )) = PTRICB(STEP( INODE )) - NUPDATE
          ELSE
             PTRICB(STEP( INODE )) = -1
          ENDIF
        ELSE
 210      CONTINUE
          CALL DMUMPS_BUF_SEND_VCB( NRHS_B, INODE, FPERE, 
     &          NCB, LD_WCBCB,
     &          NUPDATE,
     &          IW( J3 + 1 ), WCB( PCB_COURANT ), JBDEB, JBFIN,
     &          RHSCOMP, 1, 1, -9999, -9999,
     &          KEEP,
     &          MUMPS_PROCNODE(PROCNODE_STEPS(STEP(FPERE)), KEEP(199)),
     &          ContVec,
     &          COMM, IERR )
          IF ( IERR .EQ. -1 ) THEN
            CALL DMUMPS_SOLVE_RECV_AND_TREAT( .FALSE., FLAG,
     &            BUFR, LBUFR, LBUFR_BYTES,
     &            MYID, SLAVEF, COMM,
     &            N, NRHS, IPOOL, LPOOL, LEAF,
     &            NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &            IWCB, LIWCB,
     &            WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &            PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP,
     &            PROCNODE_STEPS,
     &            RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &            , FROM_PP
     &            )
            IF ( INFO( 1 )  .LT. 0 )  THEN
              ERROR_WAS_BROADCASTED = .TRUE.
              GOTO 270
            ENDIF
            GOTO 210
          ELSE IF ( IERR .EQ. -2 ) THEN
            INFO( 1 ) = -17
            INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     &           ( NUPDATE + 3 ) * KEEP( 34 )
            ERROR_WAS_BROADCASTED = .FALSE.
            GOTO 270
          ELSE IF ( IERR .EQ. -3 ) THEN
            INFO( 1 ) = -20
            INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     &           ( NUPDATE + 3 ) * KEEP( 34 )
            ERROR_WAS_BROADCASTED = .FALSE.
            GOTO 270
          END IF
        ENDIF
      END IF
      IF ( NSLAVES .NE. 0 .AND. MTYPE .eq. 1
     &     .and. NPIV .NE. 0 ) THEN
         DO ISLAVE = 1, NSLAVES
            PDEST = IW( PTRIST(STEP(INODE)) + 5 + ISLAVE +KEEP(IXSZ))
            CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &           KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &           ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &           ISLAVE, NCB - NELIM, 
     &           NSLAVES, 
     &           Effective_CB_Size, FirstIndex )
 222        CONTINUE
                 CALL DMUMPS_BUF_SEND_MASTER2SLAVE( NRHS_B,
     &           INODE, FPERE,
     &           Effective_CB_Size, LD_WCBCB, LD_WCBPIV, NPIV,
     &           JBDEB, JBFIN, 
     &           WCB( PCB_COURANT + NELIM + FirstIndex - 1 ),
     &           WCB( PPIV_COURANT ),
     &           PDEST, COMM, KEEP, IERR )
            IF ( IERR .EQ. -1 ) THEN
               CALL DMUMPS_SOLVE_RECV_AND_TREAT( .FALSE., FLAG,
     &              BUFR, LBUFR, LBUFR_BYTES,
     &              MYID, SLAVEF, COMM,
     &              N, NRHS, IPOOL, LPOOL, LEAF,
     &              NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     &              IWCB, LIWCB,
     &              WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &              PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP,
     &              PROCNODE_STEPS,
     &              RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &              , FROM_PP
     &              )
               IF ( INFO( 1 )  .LT. 0 ) THEN
                 ERROR_WAS_BROADCASTED = .TRUE.
                 GOTO 270
               ENDIF
               GOTO 222
            ELSE IF ( IERR .EQ. -2 ) THEN
               INFO( 1 ) = -17
               INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS_B*KEEP(35) +
     &               6 * KEEP( 34 )
               ERROR_WAS_BROADCASTED = .FALSE.
               GOTO 270
            ELSE IF ( IERR .EQ. -3 ) THEN
               INFO( 1 ) = -20
               INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS_B*KEEP(35) +
     &              6 * KEEP( 34 )
               ERROR_WAS_BROADCASTED = .FALSE.
               GOTO 270
            END IF
         END DO
      END IF
      PLEFTWCB = PLEFTWCB - int(LIELL,8)*int(NRHS_B,8)
 270  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_NODE_FWD
      RECURSIVE SUBROUTINE DMUMPS_SOLVE_RECV_AND_TREAT( BLOQ, FLAG,
     &           BUFR, LBUFR, LBUFR_BYTES,
     &           MYID, SLAVEF, COMM,
     &           N, NRHS, IPOOL, LPOOL, LEAF,
     &           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     &           IWCB, LIWCB,
     &           WCB, LWCB, POSWCB,
     &           PLEFTWCB, POSIWCB,
     &           PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP, PROCNODE_STEPS,
     &           RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &            , FROM_PP
     &            )
      IMPLICIT NONE
      LOGICAL BLOQ
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, NRHS, LPOOL, LEAF, NBFIN
      INTEGER LIWCB, POSIWCB
      INTEGER(8) :: POSWCB, PLEFTWCB
      INTEGER LIW
      INTEGER(8), INTENT(IN) :: LA, LWCB
      INTEGER INFO( 80 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL)
      INTEGER NSTK_S( KEEP(28) )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      DOUBLE PRECISION WCB( LWCB ), A( LA )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      LOGICAL FLAG
      INTEGER LRHSCOMP, POSINRHSCOMP_FWD(N)
      DOUBLE PRECISION RHSCOMP(LRHSCOMP,NRHS)
      LOGICAL, intent(in) :: FROM_PP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER MSGSOU, MSGTAG, MSGLEN
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
         FLAG = .FALSE.
           CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     &                     COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     &                   FLAG, STATUS, IERR )
      END IF
      IF ( FLAG ) THEN
         KEEP(266) = KEEP(266) -1
         MSGSOU = STATUS( MPI_SOURCE )
         MSGTAG = STATUS( MPI_TAG )
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
         ELSE
           CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED,
     &                  MSGSOU, MSGTAG, COMM, STATUS, IERR )
           CALL DMUMPS_TRAITER_MESSAGE_SOLVE( BUFR, LBUFR, LBUFR_BYTES,
     &          MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &          N, NRHS, IPOOL, LPOOL, LEAF,
     &          NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &          IWCB, LIWCB,
     &          WCB, LWCB, POSWCB,
     &          PLEFTWCB, POSIWCB,
     &          PTRICB, INFO, KEEP,KEEP8, DKEEP, STEP,
     &          PROCNODE_STEPS, 
     &          RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &          , FROM_PP
     &          )
         END IF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_RECV_AND_TREAT
      SUBROUTINE DMUMPS_RHSCOMP_TO_WCB(
     &      NPIV, NCB, LIELL, CBINITZERO, LDEQLIELLPANEL,
     &      RHSCOMP, LRHSCOMP, NRHS_B,
     &      POSINRHSCOMP_FWD, N,
     &      WCB,
     &      IW, LIW, J1, J3, J2, KEEP, DKEEP)
      IMPLICIT NONE
      INTEGER, INTENT( IN )    :: NPIV, NCB, LIELL, N,
     &                            LRHSCOMP, NRHS_B,
     &                            LIW, J1, J2, J3
      LOGICAL, INTENT( IN )    :: LDEQLIELLPANEL
      LOGICAL, INTENT( IN )    :: CBINITZERO
      INTEGER, INTENT( IN )    :: POSINRHSCOMP_FWD( N ), IW( LIW )
      DOUBLE PRECISION, INTENT( INOUT ) :: RHSCOMP( LRHSCOMP, NRHS_B )
      DOUBLE PRECISION, INTENT( OUT )   :: WCB( int(LIELL,8)*
     &                                 int(NRHS_B,8) )
      INTEGER                  :: KEEP(500)
      DOUBLE PRECISION                     :: DKEEP(150)
      INTEGER, PARAMETER :: ZERO = 0.0D0
      INTEGER(8), PARAMETER :: PPIV_COURANT = 1_8
      INTEGER(8) :: PCB_COURANT
      INTEGER :: LD_WCBCB, LD_WCBPIV, J, JJ, K, IPOSINRHSCOMP
      INTEGER(8) :: IFR8, IFR_ini8 
      INCLUDE 'mpif.h' 
!$    LOGICAL :: OMP_FLAG
      IF ( LDEQLIELLPANEL ) THEN
        LD_WCBPIV = LIELL
        LD_WCBCB  = LIELL
        PCB_COURANT = PPIV_COURANT + NPIV
      ELSE
        LD_WCBPIV = NPIV
        LD_WCBCB  = NCB
        PCB_COURANT = PPIV_COURANT + NPIV * NRHS_B
      ENDIF
      IF ( LDEQLIELLPANEL ) THEN  
        DO K=1, NRHS_B
           IFR8 = PPIV_COURANT+int(K-1,8)*int(LD_WCBPIV,8)-1_8
           IPOSINRHSCOMP = POSINRHSCOMP_FWD(IW(J1)) 
           DO JJ = J1, J3
              IFR8 = IFR8 + 1_8
              WCB(IFR8) = RHSCOMP(IPOSINRHSCOMP,K) 
              IPOSINRHSCOMP = IPOSINRHSCOMP + 1
           ENDDO
           IF (NCB.GT.0 .AND. .NOT. CBINITZERO) THEN
              DO JJ = J3+1, J2
                 J = IW(JJ)
                 IFR8 = IFR8 + 1_8
                 IPOSINRHSCOMP = abs(POSINRHSCOMP_FWD(J))
                 WCB(IFR8) = RHSCOMP(IPOSINRHSCOMP,K) 
                 RHSCOMP (IPOSINRHSCOMP,K) = ZERO
              ENDDO
           ENDIF
        ENDDO
      ELSE                          
        PCB_COURANT = PPIV_COURANT + LD_WCBPIV*NRHS_B
        IFR8 = PPIV_COURANT - 1_8
        IFR_ini8 = IFR8
        IPOSINRHSCOMP = POSINRHSCOMP_FWD(IW(J1)) 
!$      OMP_FLAG = ( NRHS_B .GE. KEEP(362) .AND.
!$   &               int(NCB,8)*int(NRHS_B,8) .GE. KEEP(363) )
!$OMP PARALLEL DO PRIVATE(JJ,IFR8) IF(OMP_FLAG)
        DO K=1, NRHS_B
          IFR8 = IFR_ini8 + int(K-1,8)*int(NPIV,8)
          DO  JJ = J1, J3
           WCB(IFR8+int(JJ-J1+1,8)) = 
     &          RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        IFR8 = PCB_COURANT - 1_8
        IF (NCB.GT.0 .AND. .NOT. CBINITZERO) THEN
          IFR_ini8 = IFR8
!$        OMP_FLAG = ( NRHS_B.GE.KEEP(362) .AND.
!$   &                 NCB*NRHS_B .GE. KEEP(363) )
!$OMP PARALLEL DO PRIVATE (IFR8, JJ, J, IPOSINRHSCOMP) IF (OMP_FLAG)
          DO K=1, NRHS_B
            IFR8 = IFR_ini8+(K-1)*NCB
            DO JJ = J3 + 1, J2
              J = IW(JJ)
              IPOSINRHSCOMP = abs(POSINRHSCOMP_FWD(J))
              WCB(IFR8+int(JJ-J3,8)) = RHSCOMP(IPOSINRHSCOMP,K)
              RHSCOMP(IPOSINRHSCOMP,K)=ZERO
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDIF
      ENDIF                     
      IF ( CBINITZERO ) THEN
!$      OMP_FLAG = int(NCB,8)*int(NRHS_B,8) .GE. KEEP(363)
!$OMP PARALLEL DO COLLAPSE(2) IF ( OMP_FLAG )
        DO K = 1, NRHS_B
          DO JJ = 1, NCB
            WCB(PCB_COURANT+int(K-1,8)*int(LD_WCBCB,8)+JJ-1_8) = ZERO
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_RHSCOMP_TO_WCB
