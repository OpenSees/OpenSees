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
      SUBROUTINE DMUMPS_COMPRESS_LU(SIZE_INPLACE,
     &MYID,N,IOLDPS,TYPE,IW, LIW, A, LA,
     &POSFAC, LRLU, LRLUS, IWPOS, PTRAST, PTRFAC, STEP, KEEP,KEEP8,
     &SSARBR,INODE,IERR
     & , LRGROUPS, NASS
     &)
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
!$    USE OMP_LIB
      USE DMUMPS_LR_CORE
      IMPLICIT NONE
      INTEGER MYID
      INTEGER IOLDPS, TYPE, LIW, N, KEEP(500)
      INTEGER(8) :: SIZE_INPLACE, LA, POSFAC, LRLU, LRLUS
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) KEEP8(150)
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER IWPOS, LDLT
      INTEGER STEP( N )
      INTEGER (8) :: PTRFAC(KEEP(28))
      LOGICAL SSARBR
      INTEGER IOLDSHIFT, IPSSHIFT
      INTEGER LRGROUPS(N), NASS
      INCLUDE 'mumps_headers.h'
      INTEGER LCONT, NELIM, NROW, NPIV, INTSIZ
      INTEGER NFRONT, NSLAVES
      INTEGER IPS, IPSIZE
      INTEGER(8) :: SIZELU, SIZECB, IAPOS, I, SIZESHIFT, ITMP8
      LOGICAL MOVEPTRAST
      LOGICAL LRCOMPRESS_PANEL
      INTEGER INODE
      INTEGER IERR
      INTEGER PARPIV_T1
      LOGICAL LR_ACTIVATED
      IERR=0
      LDLT = KEEP(50)
      IOLDSHIFT = IOLDPS + KEEP(IXSZ)
      IF ( IW( IOLDSHIFT ) < 0 ) THEN
        write(*,*) ' ERROR 1 compressLU:Should not point to a band.'
        CALL MUMPS_ABORT()
      ELSE IF ( IW( IOLDSHIFT + 2 ) < 0 ) THEN
        write(*,*) ' ERROR 2 compressLU:Stack not performed yet',
     &  IW(IOLDSHIFT + 2)
        CALL MUMPS_ABORT()
      ENDIF
      LCONT  = IW( IOLDSHIFT )
      NELIM  = IW( IOLDSHIFT + 1 )
      NROW   = IW( IOLDSHIFT + 2 )
      NPIV   = IW( IOLDSHIFT + 3 )
      IAPOS  = PTRFAC(IW( IOLDSHIFT + 4 ))
      NSLAVES= IW( IOLDSHIFT + 5 )
      NFRONT = LCONT + NPIV
      INTSIZ = IW(IOLDPS+XXI)
      LRCOMPRESS_PANEL = (IW(IOLDPS+XXLR).GE.2)
      IF ( (NSLAVES > 0  .AND. TYPE .NE. 2) .OR. 
     &     (NSLAVES .eq. 0 .AND. TYPE .EQ. 2 ) ) THEN
          WRITE(*,*) ' ERROR 3 compressLU: problem with level of inode'
          CALL MUMPS_ABORT()
      END IF
      IF (LDLT.EQ.0) THEN
        SIZELU = int(LCONT + NROW, 8) * int(NPIV,8)
      ELSE
        SIZELU =   int(NROW,8) * int(NPIV,8)
      ENDIF
      IF ( TYPE .EQ. 2 ) THEN
        IF (LDLT.EQ.0) THEN
          SIZECB = int(NELIM,8) * int(LCONT,8)
        ELSE
          IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
            SIZECB = int(NELIM+1,8) * int(NELIM + NPIV,8)
          ELSE
            SIZECB = int(NELIM,8) * int(NELIM + NPIV,8)
          ENDIF
        ENDIF
      ELSE
        LR_ACTIVATED =  (IW(IOLDPS+XXLR).GT.0)
        IF (LDLT.EQ.0) THEN
         CALL DMUMPS_SET_PARPIVT1 ( INODE, NFRONT, NELIM+NPIV, 
     &     KEEP, LR_ACTIVATED, PARPIV_T1) 
         IF (PARPIV_T1.EQ.0) THEN
          SIZECB = int(LCONT,8) * int(LCONT,8)
         ELSE
          SIZECB = int(LCONT,8) * int(LCONT,8) + int(NELIM + NPIV,8)
         ENDIF
        ELSE
         CALL DMUMPS_SET_PARPIVT1 ( INODE, NFRONT, NELIM+NPIV, 
     &     KEEP, LR_ACTIVATED, PARPIV_T1) 
         IF (PARPIV_T1.EQ.0) THEN
           SIZECB = int(NROW,8) * int(LCONT,8)
         ELSE
           SIZECB = int(NROW,8) * int(LCONT,8) + int(NELIM + NPIV,8)
         ENDIF
        ENDIF
      END IF
      CALL MUMPS_SUBTRI8TOARRAY( IW(IOLDPS+XXR), SIZECB )
      IF ((KEEP(201).NE.0)
     &     .OR.(LRCOMPRESS_PANEL.AND.KEEP(486).EQ.2)
     &      ) THEN
        SIZESHIFT = SIZELU
      ELSE
        SIZESHIFT = 0_8
        IF (SIZECB.EQ.0_8) THEN
          GOTO 500
        ENDIF
      ENDIF
      IF (KEEP(201).EQ.2) THEN
         IF (KEEP(405) .EQ. 0) THEN
           KEEP8(31)=KEEP8(31)+SIZELU
           CALL DMUMPS_NEW_FACTOR(INODE,PTRFAC,KEEP,KEEP8,
     &          A,LA,SIZELU, IERR)
         ELSE
!$OMP CRITICAL(critical_old_ooc)
           KEEP8(31)=KEEP8(31)+SIZELU
           CALL DMUMPS_NEW_FACTOR(INODE,PTRFAC,KEEP,KEEP8,
     &          A,LA,SIZELU, IERR)
!$OMP END CRITICAL(critical_old_ooc)
         ENDIF
         IF(IERR.LT.0)THEN
            WRITE(*,*)MYID,': Internal error in DMUMPS_NEW_FACTOR'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      IF ( IOLDPS + INTSIZ .NE. IWPOS ) THEN
         IPS = IOLDPS + INTSIZ
         MOVEPTRAST = .FALSE.
         DO WHILE ( IPS .NE. IWPOS )
           IPSIZE = IW(IPS+XXI)
           IPSSHIFT = IPS + KEEP(IXSZ)
           IF ( IPSIZE .LE. 0 .OR. IPS .GT. IWPOS ) THEN
             WRITE(*,*) " Internal error 1 DMUMPS_COMPRESS_LU"
             WRITE(*,*) " IOLDPS, INTSIZ, IWPOS, LIW=",
     &                    IOLDPS, INTSIZ, IWPOS, LIW
             WRITE(*,*) " IWPOS, IPS, IPSIZE =", IWPOS, IPS, IPSIZE
             WRITE(*,*) " Header at IOLDPS =",
     &                  IW(IOLDPS:IOLDPS+KEEP(IXSZ)+5)
             WRITE(*,*) " Header at IPS =",
     &                  IW(IPS:IPS+KEEP(IXSZ)+5)
             CALL MUMPS_ABORT()
           ENDIF
           IF (IPS+IPSIZE .GT. IWPOS) THEN
             WRITE(*,*) " Internal error 2 DMUMPS_COMPRESS_LU"
             WRITE(*,*) " IOLDPS, INTSIZ, IWPOS, LIW=",
     &                    IOLDPS, INTSIZ, IWPOS, LIW
             WRITE(*,*) " IWPOS, IPS, IPSIZE =", IWPOS, IPS, IPSIZE
             WRITE(*,*) " Header at IOLDPS =",
     &                  IW(IOLDPS:IOLDPS+KEEP(IXSZ)+5)
             WRITE(*,*) " Header at IOLDPS+INTSIZ =",
     &                  IW(IOLDPS+INTSIZ:IOLDPS+INTSIZ+KEEP(IXSZ)+5)
             WRITE(*,*) " Header at IPS =",
     &                  IW(IPS:IPS+KEEP(IXSZ)+5)
             WRITE(*,*) " ========================== "
             WRITE(*,*) " Headers starting at IOLDPS:"
             IPS = IOLDPS
             DO WHILE (IPS .LE. IWPOS)
               WRITE(*,*) " -> new IW header at position" , IPS, ":",
     &                  IW(IPS:IPS+KEEP(IXSZ)+5)
               IPS = IPS + IW(IPS+XXI)
             ENDDO
             CALL MUMPS_ABORT()
           ENDIF
           IF ( IW( IPSSHIFT + 2 ) < 0 ) THEN
             NFRONT = IW( IPSSHIFT )
             IF (IW(IPSSHIFT+4) .LT. 0) THEN
               WRITE(*,*) " Internal error 3 DMUMPS_COMPRESS_LU"
               WRITE(*,*) " IPS,IPSSHIFT,IWPOS=" ,IPS,IPSSHIFT,IWPOS
               WRITE(*,*) " Header at IPS =", IW(IPS:IPS+KEEP(IXSZ)+5)
             ENDIF
             PTRFAC(IW(IPSSHIFT+4))=PTRFAC(IW(IPSSHIFT+4)) -
     &               SIZECB - SIZESHIFT
             MOVEPTRAST = .TRUE.
             PTRAST(IW(IPSSHIFT+4))=PTRAST(IW(IPSSHIFT+4))-SIZECB
     &               - SIZESHIFT
           ELSE IF ( IW( IPSSHIFT ) < 0 ) THEN
             IF (IW(IPSSHIFT+3) .LT. 0) THEN
               WRITE(*,*) " Internal error 4 DMUMPS_COMPRESS_LU"
               WRITE(*,*) " IPS,IPSSHIFT,IWPOS=" ,IPS,IPSSHIFT,IWPOS
               WRITE(*,*) " Header at IPS =", IW(IPS:IPS+KEEP(IXSZ)+5)
             ENDIF
             PTRFAC(IW(IPSSHIFT+3)) = PTRFAC(IW(IPSSHIFT+3))
     &                                  -SIZECB-SIZESHIFT
           ELSE
             NFRONT = IW( IPSSHIFT ) + IW( IPSSHIFT + 3 )
             IF (IW(IPSSHIFT+4) .LT. 0) THEN
               WRITE(*,*) " Internal error 4 DMUMPS_COMPRESS_LU"
               WRITE(*,*) " IPS,IPSSHIFT,IWPOS=" ,IPS,IPSSHIFT,IWPOS
               WRITE(*,*) " Header at IPS =", IW(IPS:IPS+KEEP(IXSZ)+5)
             ENDIF
             PTRFAC(IW( IPSSHIFT + 4 )) = 
     &           PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB - SIZESHIFT
           END IF
           IPS = IPS + IPSIZE
         END DO
         IF (SIZECB+SIZESHIFT .NE. 0_8) THEN
           DO I=IAPOS+SIZELU-SIZESHIFT, POSFAC-SIZECB-SIZESHIFT-1_8
             A( I ) = A( I + SIZECB + SIZESHIFT)
           END DO
         END IF
      ENDIF
      POSFAC = POSFAC  - (SIZECB+SIZESHIFT)
      LRLU   = LRLU    + (SIZECB+SIZESHIFT)
      ITMP8  = (SIZECB+SIZESHIFT) - SIZE_INPLACE
      LRLUS  = LRLUS   + ITMP8
      IF (KEEP(405) .EQ. 0) THEN
        KEEP8(69) = KEEP8(69) - ITMP8
      ELSE
!$OMP   ATOMIC UPDATE
        KEEP8(69) = KEEP8(69) - ITMP8
!$OMP   END ATOMIC
      ENDIF
 500  CONTINUE
      IF (LRCOMPRESS_PANEL.AND.KEEP(486).EQ.2) THEN
        CALL DMUMPS_LOAD_MEM_UPDATE(SSARBR,.FALSE.,
     &     LA-LRLUS,SIZELU-SIZESHIFT,-(SIZECB+SIZESHIFT)+SIZE_INPLACE,
     &     KEEP,KEEP8,LRLUS)
      ELSE
        CALL DMUMPS_LOAD_MEM_UPDATE(SSARBR,.FALSE.,
     &     LA-LRLUS,SIZELU,-SIZECB+SIZE_INPLACE,
     &     KEEP,KEEP8,LRLUS)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_COMPRESS_LU
      SUBROUTINE DMUMPS_STACK_BAND( N, ISON, 
     &    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA, 
     &    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP, 
     &    IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER,
     &    IFLAG, IERROR, SLAVEF, PROCNODE_STEPS, DAD, MYID, COMM,
     &    KEEP, KEEP8, DKEEP, TYPE_SON
     &     )
!$    USE OMP_LIB
      USE DMUMPS_OOC
      USE DMUMPS_LOAD
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_DYNPTR
      IMPLICIT NONE
      INTEGER(8) :: LA, LRLU, LRLUS, POSFAC, IPTRLU
      INTEGER N, ISON, LIW, IWPOS, IWPOSCB,
     &        COMP, IFLAG, IERROR, SLAVEF, MYID, COMM,
     &        TYPE_SON
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION DKEEP(230)
      INTEGER, INTENT(IN) :: PROCNODE_STEPS(KEEP(28)), DAD(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), STEP(N), 
     & PIMASTER(KEEP(28)), IW(LIW)
      INTEGER PTLUST_S(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      DOUBLE PRECISION OPELIW
      DOUBLE PRECISION FLOP1, FLOP1_EFFECTIVE
      DOUBLE PRECISION A( LA )
      INTEGER(8) :: LREQA, POSA, POSALOC, OLDPOS, JJ
      INTEGER  NFRONT, NCOL_L, NROW_L, LREQI, NSLAVES_L,
     &         POSI, I, IROW_L, ICOL_L, LDA_BAND, NASS
      LOGICAL NONEED_TO_COPY_FACTORS
      INTEGER(8) :: LREQA_HEADER
      INTEGER LIWFAC, STRAT, TYPEFile, NextPivDummy,
     &        IOLDPS_CB
      LOGICAL LAST_CALL
      TYPE(IO_BLOCK) :: MonBloc 
      INTEGER LRSTATUS
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0d0)
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: SON_A
      INTEGER(8) :: IACHK, SIZFR_SON_A, ITMP8
      FLOP1 = ZERO
      NCOL_L = IW( PTRIST(STEP( ISON )) + 3 + KEEP(IXSZ) )
      NROW_L = IW( PTRIST(STEP( ISON )) + 2 + KEEP(IXSZ) )
      NSLAVES_L = IW( PTRIST(STEP( ISON )) + 5 + KEEP(IXSZ) )
      LRSTATUS = IW( PTRIST(STEP( ISON )) + XXLR)
      LDA_BAND = NCOL_L + IW( PTRIST(STEP( ISON )) + KEEP(IXSZ) )
      IF  ( KEEP(50) .eq. 0 ) THEN
        NFRONT = LDA_BAND
      ELSE
        NFRONT = IW( PTRIST(STEP( ISON )) + 7 + KEEP(IXSZ) )
      END IF
      IF (KEEP(201).EQ.1) THEN 
          IOLDPS_CB = PTRIST(STEP( ISON ))
          LIWFAC    = IW(IOLDPS_CB+XXI)
          TYPEFile  = TYPEF_L
          NextPivDummy      = -8888 
          MonBloc%INODE    = ISON
          MonBloc%MASTER   = .FALSE.   
          MonBloc%Typenode =  2        
          MonBloc%NROW     = NROW_L
          MonBloc%NCOL     = LDA_BAND
          MonBloc%NFS      = IW(IOLDPS_CB+1+KEEP(IXSZ))
          MonBloc%LastPiv  = NCOL_L    
          MonBloc%LastPanelWritten_L=-9999 
          MonBloc%LastPanelWritten_U=-9999 
          NULLIFY(MonBloc%INDICES)
          STRAT        = STRAT_WRITE_MAX
          LAST_CALL    = .TRUE.
          MonBloc%Last = .TRUE.
          CALL DMUMPS_DM_SET_DYNPTR(IW(PTRIST(STEP(ISON))+XXS),
     &    A, LA, PTRAST(STEP(ISON)),
     &    IW(PTRIST(STEP(ISON))+XXD),
     &    IW(PTRIST(STEP(ISON))+XXR),
     &    SON_A, IACHK, SIZFR_SON_A)
          CALL DMUMPS_OOC_IO_LU_PANEL_I
     &          ( STRAT, TYPEFile, 
     &           SON_A(IACHK), SIZFR_SON_A, MonBloc,
     &           NextPivDummy, NextPivDummy,
     &           IW(IOLDPS_CB), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG,LAST_CALL )
          IF ((NCOL_L.EQ.0).OR.(NROW_L.EQ.0)) THEN 
          ENDIF
      ENDIF  
      NONEED_TO_COPY_FACTORS = ((KEEP(201).EQ.1) .OR. (KEEP(201).EQ.-1)
     &               .OR. (LRSTATUS.GE.2.AND.KEEP(486).EQ.2)
     &                         )  
      IF ((NCOL_L.EQ.0).OR.(NROW_L.EQ.0)) THEN 
        GOTO 80
      ENDIF
      LREQI   = 4 + NCOL_L + NROW_L + KEEP(IXSZ)
      LREQA_HEADER =  int(NCOL_L,8) * int(NROW_L,8)
      IF (NONEED_TO_COPY_FACTORS) THEN 
        LREQA = 0_8
      ELSE
        LREQA   = LREQA_HEADER
      ENDIF
      IF ( LRLU .LT. LREQA .OR.
     &  IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LREQA ) THEN
          IFLAG  = -9
          CALL MUMPS_SET_IERROR(LREQA - LRLUS, IERROR)
          GO TO 700
        END IF
        CALL DMUMPS_COMPRE_NEW( N,KEEP(28), IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &        KEEP(IXSZ), COMP, DKEEP(97),
     &        MYID, SLAVEF, KEEP(199), PROCNODE_STEPS, DAD )
        IF ( LRLU .NE. LRLUS ) THEN
               WRITE(*,*) 'PB compress DMUMPS_STACK_BAND:LRLU,LRLUS=',
     &         LRLU, LRLUS
               IFLAG = -9
               CALL MUMPS_SET_IERROR(LREQA - LRLUS, IERROR)
               GOTO 700
        END IF
        IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
          IFLAG  = -8
          IERROR = IWPOS + LREQI - 1 - IWPOSCB
          GOTO 700
        END IF
      END IF
      IF (.NOT. NONEED_TO_COPY_FACTORS) THEN
        POSA = POSFAC
        POSFAC = POSFAC + LREQA
        LRLU = LRLU - LREQA
        LRLUS = LRLUS - LREQA
        KEEP8(67) = min(LRLUS, KEEP8(67))
        KEEP8(69) = KEEP8(69) + LREQA
        KEEP8(68) = max(KEEP8(69), KEEP8(68))
        IF(KEEP(201).NE.2)THEN
           CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &          LA-LRLUS,LREQA,LREQA,KEEP,KEEP8,LRLUS)
        ELSE
           CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &          LA-LRLUS,0_8,LREQA,KEEP,KEEP8,LRLUS)
        ENDIF
      ENDIF
      POSI = IWPOS
      IWPOS = IWPOS + LREQI
      PTLUST_S(STEP( ISON )) = POSI
      IW(POSI:POSI+KEEP(IXSZ)-1)=-99999
      IW(POSI+XXS)=-9999
      IW(POSI+XXI)=LREQI
      CALL MUMPS_STOREI8(0_8, IW(POSI+XXD))
      CALL MUMPS_STOREI8(LREQA, IW(POSI+XXR))
      CALL MUMPS_STOREI8(LREQA_HEADER, IW(POSI+XXR))
      IW(POSI+XXLR) = LRSTATUS
      IW(POSI+XXF) = IW(PTRIST(STEP(ISON))+XXF)
      POSI=POSI+KEEP(IXSZ)
      IW( POSI     ) = - NCOL_L
      IW( POSI + 1 ) =   NROW_L
      IW( POSI + 2 ) =   NFRONT - NCOL_L
      IW( POSI + 3 ) =   STEP(ISON)
      IF (.NOT. NONEED_TO_COPY_FACTORS) THEN
        PTRFAC(STEP(ISON)) = POSA
      ELSE
        PTRFAC(STEP(ISON)) = -77777_8
      ENDIF
      IROW_L = PTRIST(STEP(ISON)) + 6 + NSLAVES_L + KEEP(IXSZ)
      ICOL_L = PTRIST(STEP(ISON)) + 6 + NROW_L + NSLAVES_L + KEEP(IXSZ)
      DO I = 1, NROW_L
        IW( POSI+3+I ) = IW( IROW_L+I-1 )
      ENDDO
      DO I = 1, NCOL_L
        IW( POSI+NROW_L+3+I) = IW( ICOL_L+I-1 )
      ENDDO
      IF (.NOT.NONEED_TO_COPY_FACTORS) THEN
        CALL DMUMPS_DM_SET_DYNPTR(IW(PTRIST(STEP(ISON))+XXS),
     &  A, LA, PTRAST(STEP(ISON)),
     &  IW(PTRIST(STEP(ISON))+XXD),
     &  IW(PTRIST(STEP(ISON))+XXR),
     &  SON_A, IACHK, SIZFR_SON_A)
        POSALOC = POSA
        DO I = 1, NROW_L
          OLDPOS =  IACHK + int(I-1,8)*int(LDA_BAND,8)
          DO JJ = 0_8, int(NCOL_L-1,8)
            A( POSALOC+JJ ) = SON_A( OLDPOS+JJ )
          ENDDO
          POSALOC = POSALOC + int(NCOL_L,8)
        END DO
      ENDIF
      IF (KEEP(201).EQ.2) THEN
       KEEP8(31)=KEEP8(31)+LREQA
      ENDIF
      ITMP8 = int(NCOL_L,8) * int(NROW_L,8)
      IF (KEEP(405) .EQ.1) THEN
!$OMP   ATOMIC UPDATE
        KEEP8(10) = KEEP8(10) + ITMP8
!$OMP   END ATOMIC
      ELSE
        KEEP8(10) = KEEP8(10) + ITMP8
      ENDIF
      IF (KEEP(201).EQ.2) THEN 
        CALL DMUMPS_NEW_FACTOR(ISON,PTRFAC,KEEP,KEEP8,A,LA,LREQA,IFLAG)
        IF(IFLAG.LT.0)THEN
          WRITE(*,*)MYID,': Internal error in DMUMPS_NEW_FACTOR'
          IERROR=0
          GOTO 700
        ENDIF
        POSFAC = POSFAC - LREQA
        LRLU = LRLU + LREQA
        LRLUS = LRLUS + LREQA
!$OMP ATOMIC UPDATE
        KEEP8(69) = KEEP8(69) - LREQA
!$OMP END ATOMIC
        CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &            LA-LRLUS,LREQA,0_8,KEEP,KEEP8,LRLUS)
      ENDIF
  80  CONTINUE
      IF (TYPE_SON == 1) THEN
         GOTO 90
      ENDIF
      IF ( KEEP(50) .eq. 0 ) THEN
         FLOP1 = dble( NCOL_L * NROW_L) +
     &     dble(NROW_L*NCOL_L)*dble(2*NFRONT-NCOL_L-1)
      ELSE
         FLOP1 = dble( NCOL_L ) * dble( NROW_L )
     &         * dble( 2 * LDA_BAND - NROW_L - NCOL_L + 1)
      END IF
      OPELIW = OPELIW + FLOP1
      FLOP1_EFFECTIVE = FLOP1
      NASS = IW( PTRIST(STEP( ISON )) + 4 + KEEP(IXSZ) )
      IF ( NCOL_L .NE. NASS ) THEN
        IF ( KEEP(50).eq.0 ) THEN
           FLOP1 = dble( NASS * NROW_L) +
     &     dble(NROW_L*NASS)*dble(2*NFRONT-NASS-1)
        ELSE
           FLOP1 = dble( NASS ) * dble( NROW_L ) *
     &     dble( 2 * LDA_BAND - NROW_L - NASS + 1)
        END IF
      END IF
      CALL DMUMPS_LOAD_UPDATE(1,.FALSE.,FLOP1_EFFECTIVE-FLOP1,
     &                        KEEP,KEEP8)
      CALL DMUMPS_LOAD_UPDATE(2,.FALSE.,-FLOP1,KEEP,KEEP8)
 90   CONTINUE
      RETURN
 700  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      RETURN
      END SUBROUTINE DMUMPS_STACK_BAND
      SUBROUTINE DMUMPS_FREE_BAND( N, ISON, 
     &    PTRIST, PTRAST, IW, LIW, A, LA, 
     &    LRLU, LRLUS, IWPOSCB,
     &    IPTRLU, STEP, MYID, KEEP, KEEP8, TYPE_SON
     &     )
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_PTR,
     &                                    DMUMPS_DM_FREE_BLOCK
      IMPLICIT NONE
      include 'mumps_headers.h'
      INTEGER(8) :: LRLU, LRLUS, IPTRLU, LA
      INTEGER ISON, MYID, N, IWPOSCB, TYPE_SON
      INTEGER KEEP(500), STEP(N)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER PTRIST(KEEP(28))
      INTEGER LIW
      INTEGER IW(LIW)
      DOUBLE PRECISION A(LA)
      INTEGER ISTCHK
      INTEGER(8) :: DYN_SIZE
      DOUBLE PRECISION, DIMENSION(:), POINTER :: FORTRAN_POINTER
      ISTCHK = PTRIST(STEP(ISON))
      CALL MUMPS_GETI8( DYN_SIZE, IW(ISTCHK+XXD) )
      IF (DYN_SIZE .GT. 0_8) THEN
        CALL DMUMPS_DM_SET_PTR( PTRAST(STEP(ISON)), 
     &       DYN_SIZE, FORTRAN_POINTER )
      ENDIF
      CALL DMUMPS_FREE_BLOCK_CB_STATIC(.FALSE.,MYID, N, ISTCHK,
     &     IW, LIW, LRLU, LRLUS, IPTRLU,
     &     IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &     )
      PTRIST(STEP( ISON )) = -9999888
      PTRAST(STEP( ISON )) = -9999888_8
      RETURN
      END SUBROUTINE DMUMPS_FREE_BAND
      SUBROUTINE DMUMPS_MEM_ESTIM_BLR_ALL( SUM_OF_PEAKS, KEEP, KEEP8,
     &           MYID, COMM, N, NELT, NA, LNA, NNZ8, NA_ELT8, NSLAVES,
     &           INFO, INFOG, PROK, MP, PROKG, MPG
     &            )
      IMPLICIT NONE
      LOGICAL, INTENT(IN)   :: PROK, PROKG, SUM_OF_PEAKS
      INTEGER , INTENT(IN)  :: MYID, COMM, N, NELT, NSLAVES, 
     &                         LNA, MP, MPG
      INTEGER(8), INTENT(IN):: NA_ELT8, NNZ8
      INTEGER, INTENT(IN):: NA(LNA)
      INTEGER    :: KEEP(500), INFO(80), INFOG(80)
      INTEGER(8) :: KEEP8(150)
      INTEGER, PARAMETER :: MASTER = 0
      INTEGER    :: OOC_STAT, BLR_STRAT, BLR_CASE
      INTEGER    :: IRANK
      LOGICAL    :: EFF, PERLU_ON, COMPUTE_MAXAVG
      INTEGER(8) :: TOTAL_BYTES
      INTEGER    :: TOTAL_MBYTES
      INTEGER, DIMENSION(3) :: LRLU_UD, OOC_LRLU_UD
      PERLU_ON  = .TRUE.     
      EFF       = .FALSE.    
      COMPUTE_MAXAVG = .NOT.(NSLAVES.EQ.1 .AND. KEEP(46).EQ.1)
      IF ( PROKG.AND.SUM_OF_PEAKS) THEN
       WRITE( MPG,'(A)') 
     & ' Estimations with BLR compression of LU factors:'
       WRITE( MPG,'(A,I6,A) ')
     & ' ICNTL(38) Estimated compression rate of LU factors =', 
     &   KEEP(464), '/1000'
      ENDIF
      OOC_STAT  =  0  
      BLR_STRAT =  1  
      BLR_CASE  =  1  
      CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        MYID, N, NELT, NA(1), LNA, KEEP8(28),
     &        KEEP8(30), NSLAVES, TOTAL_MBYTES, EFF, 
     &        OOC_STAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES, 
     &        BLR_CASE, SUM_OF_PEAKS, .FALSE. , 
     &        .FALSE. 
     &         )
      CALL MUMPS_MEM_CENTRALIZE( MYID, COMM,
     &                           TOTAL_MBYTES, LRLU_UD, IRANK )
      IF (SUM_OF_PEAKS) THEN
        INFO(30) = TOTAL_MBYTES 
        IF (MYID.EQ.MASTER) THEN 
          INFOG(36) = LRLU_UD(1)
          INFOG(37) = LRLU_UD(2)
        ENDIF
      ENDIF
      IF (MYID.EQ.MASTER) THEN
        IF ( KEEP(46) .eq. 0 ) THEN
          LRLU_UD(3) = (LRLU_UD(2)-TOTAL_MBYTES)/NSLAVES
        ELSE
          LRLU_UD(3) = LRLU_UD(2)/NSLAVES
        ENDIF
      ENDIF
      IF ( PROKG.AND.SUM_OF_PEAKS ) THEN
       IF (COMPUTE_MAXAVG) THEN
         WRITE( MPG,'(A,I12) ')
     & '    Maximum estim. space in Mbytes, IC facto.    (INFOG(36)):',
     &        INFOG(36)
       ENDIF
       WRITE(MPG,'(A,I12) ')
     & '    Total space in MBytes, IC factorization      (INFOG(37)):'
     &        ,INFOG(37)
      END IF
      OOC_STAT  =  1  
      BLR_STRAT =  1  
      BLR_CASE  =  1  
      CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        MYID, N, NELT, NA(1), LNA, KEEP8(28),
     &        KEEP8(30), NSLAVES, TOTAL_MBYTES, EFF, 
     &        OOC_STAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES, 
     &        BLR_CASE, SUM_OF_PEAKS, .FALSE. , 
     &        .FALSE. 
     &         )
      CALL MUMPS_MEM_CENTRALIZE( MYID, COMM,
     &                             TOTAL_MBYTES, OOC_LRLU_UD, IRANK )
      IF (SUM_OF_PEAKS) THEN
        INFO(31) = TOTAL_MBYTES 
        IF (MYID.EQ.MASTER) THEN 
          INFOG(38)= OOC_LRLU_UD(1) 
          INFOG(39)= OOC_LRLU_UD(2)
        ENDIF
      ENDIF
      IF (MYID.EQ.MASTER) THEN
        IF ( KEEP(46) .eq. 0 ) THEN
          OOC_LRLU_UD(3) = (OOC_LRLU_UD(2)-TOTAL_MBYTES)/NSLAVES
        ELSE
          OOC_LRLU_UD(3) = OOC_LRLU_UD(2)/NSLAVES
        ENDIF
      ENDIF
      IF ( PROKG.AND.SUM_OF_PEAKS  ) THEN
       IF (COMPUTE_MAXAVG) THEN
         WRITE( MPG,'(A,I12) ')
     & '    Maximum estim. space in Mbytes, OOC facto.   (INFOG(38)):',
     &        INFOG(38)
       ENDIF
       WRITE(MPG,'(A,I12) ')
     & '    Total space in MBytes, OOC factorization     (INFOG(39)):'
     &        ,INFOG(39)
      END IF
      END SUBROUTINE DMUMPS_MEM_ESTIM_BLR_ALL
      SUBROUTINE DMUMPS_MAX_MEM( KEEP, KEEP8,
     &           MYID, N, NELT, NA, LNA, NNZ8, NA_ELT8, NSLAVES,
     &           MEMORY_MBYTES, EFF, OOC_STRAT, BLR_STRAT, PERLU_ON, 
     &           MEMORY_BYTES, 
     &           BLR_CASE, SUM_OF_PEAKS, MEM_EFF_ALLOCATED,
     &           UNDER_L0_OMP
     &           )
      IMPLICIT NONE
      LOGICAL, INTENT(IN)    :: EFF, PERLU_ON, UNDER_L0_OMP
      INTEGER, INTENT(IN)    :: OOC_STRAT, BLR_STRAT
      INTEGER, INTENT(IN)    :: KEEP(500)
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      INTEGER, INTENT(IN)    :: MYID, N, NELT, NSLAVES, LNA
      INTEGER(8), INTENT(IN) :: NA_ELT8, NNZ8
      INTEGER, INTENT(IN)    :: NA(LNA)
      INTEGER(8), INTENT(OUT):: MEMORY_BYTES
      INTEGER, INTENT(OUT)   :: MEMORY_MBYTES
      INTEGER, INTENT(IN)    :: BLR_CASE
      LOGICAL, INTENT(IN)    :: SUM_OF_PEAKS
      LOGICAL, INTENT(IN)    :: MEM_EFF_ALLOCATED 
      INTEGER  :: MUMPS_GET_POOL_LENGTH
      EXTERNAL :: MUMPS_GET_POOL_LENGTH
      INTEGER(8) :: MemEstimGlobal
      LOGICAL    :: I_AM_SLAVE, I_AM_MASTER
      INTEGER    :: PERLU, NBRECORDS
      INTEGER(8) :: NB_REAL 
      INTEGER(8) :: TEMP, NB_BYTES, NB_INT
      INTEGER    :: DMUMPS_LBUF_INT
      INTEGER(8) :: DMUMPS_LBUFR_BYTES8, DMUMPS_LBUF8
      INTEGER    :: NBUFS
      INTEGER(8) :: TEMPI
      INTEGER(8) :: TEMPR
      INTEGER    :: MIN_PERLU
      INTEGER(8) :: BUF_OOC, BUF_OOC_PANEL, BUF_OOC_NOPANEL
      INTEGER(8) :: OOC_NB_FILE_TYPE
      INTEGER(8) :: NSTEPS8, N8, NELT8
      INTEGER(8) :: I8OVERI
      INTEGER(8) :: SUM_NRLADU_underL0,
     &              SUM_NRLADU_if_LR_LU_underL0, 
     &              SUM_NRLADULR_UD_underL0, 
     &              SUM_NRLADULR_WC_underL0
      I8OVERI   = int(KEEP(10),8)
      PERLU     = KEEP(12)
      NSTEPS8   = int(KEEP(28),8)
      N8        = int(N,8)
      NELT8     = int(NELT,8)
      IF (.NOT.PERLU_ON) PERLU = 0
      I_AM_MASTER = ( MYID .eq. 0 )
      I_AM_SLAVE  = ( KEEP(46).eq. 1 .or. MYID .ne. 0 )
      TEMP    = 0_8
      NB_REAL = 0_8
      NB_BYTES = 0_8
      NB_INT  = 0_8
      IF (KEEP(235) .NE. 0 .OR. KEEP(237) .NE. 0) THEN
         NB_INT  = NB_INT + NSTEPS8
      ENDIF
      NB_INT = NB_INT + 5_8 * NSTEPS8
      NB_INT = NB_INT + NSTEPS8 + int(KEEP(56),8)*int(NSLAVES+2,8)
      NB_INT = NB_INT + 3_8 * N8
      IF (KEEP(23).ne.0 .and. I_AM_MASTER) NB_INT=NB_INT + N8
      IF (KEEP(55).eq.0) THEN
        NB_INT = NB_INT + 2_8 * N8
      ELSE
        NB_INT = NB_INT + 2_8 * ( NELT8 + 1_8 )
      ENDIF
      IF (KEEP(55) .ne. 0 ) THEN
        NB_INT = NB_INT + N8 + 1_8 + NELT8
      END IF
      NB_INT = NB_INT + int(LNA,8)
      IF ( .NOT. EFF ) THEN
       IF (I_AM_SLAVE) THEN
         IF ( KEEP8(24).EQ.0_8 ) THEN
            SUM_NRLADU_underL0          = 0_8
            SUM_NRLADU_if_LR_LU_underL0 = 0_8
            SUM_NRLADULR_UD_underL0     = 0_8
            SUM_NRLADULR_WC_underL0     = 0_8
           CALL DMUMPS_SET_MEMESTIMGLOBAL (
     &          OOC_STRAT, BLR_STRAT, BLR_CASE, SUM_OF_PEAKS,
     &          KEEP8, 
     &           SUM_NRLADU_underL0, SUM_NRLADU_if_LR_LU_underL0,
     &           SUM_NRLADULR_UD_underL0, SUM_NRLADULR_WC_underL0,
     &           KEEP8(53),  
     &           KEEP8(54),  
     &           KEEP8(33),  
     &           KEEP8(34),  
     &           KEEP8(35),  
     &           KEEP8(50),  
     &           KEEP8(36),  
     &           KEEP8(47),  
     &           KEEP8(37),  
     &           KEEP8(38),  
     &           KEEP8(39),  
     &          MemEstimGlobal 
     &          )
           IF (KEEP(400).LE.0) THEN
              NB_REAL = NB_REAL + MemEstimGlobal +
     &             int(PERLU,8)*(MemEstimGlobal / 100_8 + 1_8 )
           ELSE IF (BLR_STRAT.EQ.0) THEN
               IF ( OOC_STRAT .GT. 0 .OR. OOC_STRAT .EQ. -1 ) THEN
                NB_REAL = NB_REAL + MemEstimGlobal +
     &             int(PERLU,8)*(KEEP8(14) / 100_8 + 1_8 )
               ELSE
                NB_REAL = NB_REAL + MemEstimGlobal +
     &             int(PERLU,8)*(KEEP8(12) / 100_8 + 1_8 )
               ENDIF
           ELSE
              NB_REAL = NB_REAL + MemEstimGlobal +
     &             int(PERLU,8)*(MemEstimGlobal / 100_8 + 1_8 )
           ENDIF
         ENDIF
       ELSE
           NB_REAL = NB_REAL + 1_8
       ENDIF
      ELSE IF (I_AM_SLAVE) THEN
        IF (UNDER_L0_OMP) THEN
          IF (MEM_EFF_ALLOCATED) THEN
           NB_REAL = NB_REAL  + KEEP8(63)
          ELSE
             NB_REAL = NB_REAL + KEEP8(62)
          ENDIF
        ELSE 
          IF (MEM_EFF_ALLOCATED) THEN
            NB_REAL = NB_REAL + KEEP8(23) + KEEP8(74)
          ELSE
           NB_REAL = NB_REAL + KEEP8(67) + KEEP8(74) 
          ENDIF
        ENDIF
      ENDIF
      IF ( OOC_STRAT .GT. 0 .AND. I_AM_SLAVE ) THEN
        BUF_OOC_NOPANEL = 2_8 * KEEP8(119)
        IF (KEEP(50).EQ.0)THEN
          BUF_OOC_PANEL = 8_8 * int(KEEP(226),8)
        ELSE
          BUF_OOC_PANEL = 4_8 * int(KEEP(226),8)
        ENDIF
        IF (OOC_STRAT .EQ. 2) THEN
          BUF_OOC = BUF_OOC_NOPANEL
        ELSE
          BUF_OOC = BUF_OOC_PANEL
        ENDIF
        NB_REAL = NB_REAL + min(BUF_OOC + int(max(PERLU,0),8) *
     &          (BUF_OOC/100_8+1_8),12000000_8)
        IF (OOC_STRAT .EQ. 2) THEN
          OOC_NB_FILE_TYPE = 1_8
        ELSE
          IF (KEEP(50).EQ.0) THEN
            OOC_NB_FILE_TYPE = 2_8
          ELSE
            OOC_NB_FILE_TYPE = 1_8
          ENDIF
        ENDIF
        NB_INT = NB_INT + OOC_NB_FILE_TYPE * NSTEPS8 * I8OVERI
        NB_INT = NB_INT + OOC_NB_FILE_TYPE * NSTEPS8 * I8OVERI
        NB_INT = NB_INT + OOC_NB_FILE_TYPE * NSTEPS8
      ENDIF
      NB_REAL = NB_REAL + KEEP8(26)
      IF (KEEP(252).EQ.1 .AND. .NOT. I_AM_MASTER) THEN
        NB_REAL = NB_REAL + N8
      ENDIF
      IF ( .not. ( I_AM_SLAVE .and. I_AM_MASTER .and. KEEP(52) .eq. 0
     &         .and. KEEP(55) .ne. 0 ) ) THEN
        NB_INT  = NB_INT  + KEEP8(27)
      END IF
      IF ( I_AM_SLAVE .and. KEEP(38) .ne. 0 ) THEN
        NB_INT = NB_INT + 2_8 * N8
      END IF
      TEMPI= 0_8
      TEMPR = 0_8
      NBRECORDS = KEEP(39)
      IF (KEEP(55).eq.0) THEN
        IF (NNZ8 < int(NBRECORDS,8)) THEN
          NBRECORDS=int(NNZ8)
        ENDIF
      ELSE
        IF (NA_ELT8 < int(NBRECORDS,8)) THEN
          NBRECORDS=int(NA_ELT8)
        ENDIF
      ENDIF
      IF ( KEEP(54) .eq. 0 ) THEN
        IF ( I_AM_MASTER ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            NBUFS = NSLAVES 
          ELSE
            NBUFS = NSLAVES - 1
            IF (KEEP(55) .eq. 0 )
     &      TEMPI = TEMPI + 2_8 * N8
          END IF
          TEMPI = TEMPI + 2_8 * int(NBRECORDS,8) * int(NBUFS,8)
          TEMPR = TEMPR + int(NBRECORDS,8) * int(NBUFS,8)
        ELSE
          IF ( KEEP(55) .eq. 0 )THEN
            TEMPI = TEMPI + 2_8 * int(NBRECORDS,8)
            TEMPR = TEMPR + int(NBRECORDS,8)
          END IF
        END IF
      ELSE
        IF ( I_AM_SLAVE ) THEN
          TEMPI = TEMPI + int(1+4*NSLAVES,8) * int(NBRECORDS,8)
          TEMPR = TEMPR + int(1+2*NSLAVES,8) * int(NBRECORDS,8)
        END IF
      END IF
      TEMP = max( NB_BYTES + (NB_INT + TEMPI) * int(KEEP(34),8)
     &           + (NB_REAL+TEMPR) * int(KEEP(35),8)
     &            , TEMP )
      IF ( I_AM_SLAVE ) THEN
        IF (BLR_STRAT.NE.0) THEN
         DMUMPS_LBUFR_BYTES8 = int(KEEP(380),8) * int(KEEP(35),8)
        ELSE
         DMUMPS_LBUFR_BYTES8 = int(KEEP(44),8) * int(KEEP(35),8)
        ENDIF
        DMUMPS_LBUFR_BYTES8 = max( DMUMPS_LBUFR_BYTES8,
     &                      100000_8 )
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU=2
        ELSE
          MIN_PERLU=0
        ENDIF
        IF (KEEP(489).GT.0) THEN
          DMUMPS_LBUFR_BYTES8 = DMUMPS_LBUFR_BYTES8
     &        + int( 0.5D0 * dble(max(PERLU,MIN_PERLU))*
     &        dble(DMUMPS_LBUFR_BYTES8)/100D0,8)
        ELSE
          DMUMPS_LBUFR_BYTES8 = DMUMPS_LBUFR_BYTES8
     &        + int( 2.0D0 * dble(max(PERLU,MIN_PERLU))*
     &        dble(DMUMPS_LBUFR_BYTES8)/100D0,8)
        ENDIF
        DMUMPS_LBUFR_BYTES8 = min(DMUMPS_LBUFR_BYTES8,
     &                            int(huge (KEEP(43))-100,8))
        NB_BYTES = NB_BYTES + DMUMPS_LBUFR_BYTES8
        IF (.NOT.UNDER_L0_OMP) THEN
          IF (BLR_STRAT.NE.0) THEN
           DMUMPS_LBUF8 = int( dble(KEEP(213)) / 100.0D0
     &                       * dble(KEEP( 379 ) * KEEP( 35 )), 8 )
          ELSE
           DMUMPS_LBUF8 = int( dble(KEEP(213)) / 100.0D0
     &                       * dble(KEEP( 43 ) * KEEP( 35 )), 8 )
          ENDIF
          DMUMPS_LBUF8 = max( DMUMPS_LBUF8, 100000_8 )
          DMUMPS_LBUF8 = DMUMPS_LBUF8
     &                   + int( 2.0D0 * dble(max(PERLU,MIN_PERLU))*
     &                     dble(DMUMPS_LBUF8)/100D0, 8)
          DMUMPS_LBUF8 = min(DMUMPS_LBUF8, int(huge (KEEP(43)-100),8))
          DMUMPS_LBUF8 = max(DMUMPS_LBUF8, DMUMPS_LBUFR_BYTES8+
     &                   3_8*int(KEEP(34),8))
          NB_BYTES = NB_BYTES + DMUMPS_LBUF8
        ENDIF
        DMUMPS_LBUF_INT = ( KEEP(56) + 
     &         NSLAVES * NSLAVES ) * 5
     &               * KEEP(34)
        NB_BYTES = NB_BYTES + int(DMUMPS_LBUF_INT,8)
        IF (.NOT.EFF) THEN
          IF (UNDER_L0_OMP) THEN
            NB_INT = NB_INT + N8*KEEP(400)
          ENDIF
          IF (OOC_STRAT .GT. 0) THEN
             NB_INT = NB_INT +  int(
     &            KEEP(138) + 2 * max(PERLU,10) *
     &            ( KEEP(138) / 100 + 1 )
     &                               ,8)
          ELSE
             NB_INT = NB_INT +  int(
     &            KEEP(137) + 2 * max(PERLU,10) *
     &            ( KEEP(137) / 100 + 1 )
     &                               ,8)
          ENDIF
        ENDIF
        IF (.NOT.UNDER_L0_OMP) THEN
          IF (OOC_STRAT .GT. 0) THEN
             NB_INT = NB_INT +  int(
     &            KEEP(225) + 2 * max(PERLU,10) *
     &            ( KEEP(225) / 100 + 1 )
     &                               ,8)
          ELSE
             NB_INT = NB_INT +  int(
     &            KEEP(15) + 2 * max(PERLU,10) *
     &            ( KEEP(15) / 100 + 1 )
     &                               ,8)
          ENDIF
        ENDIF
        NB_INT = NB_INT + NSTEPS8
        NB_INT = NB_INT + NSTEPS8 * I8OVERI
        NB_INT = NB_INT + N8 + 4_8 * NSTEPS8 +
     &           int(MUMPS_GET_POOL_LENGTH(NA(1), KEEP, KEEP8),8)
        NB_INT = NB_INT + 2_8 * NSTEPS8 * I8OVERI
        IF (KEEP(494).NE.0) THEN
         NB_INT = NB_INT + N8
        ENDIF
      ENDIF
      MEMORY_BYTES = NB_BYTES + NB_INT * int(KEEP(34),8) +
     &               NB_REAL * int(KEEP(35),8)
      MEMORY_BYTES = max( MEMORY_BYTES, TEMP )
      MEMORY_MBYTES = nint( dble(MEMORY_BYTES) / dble(1000000) )
      RETURN
      END SUBROUTINE DMUMPS_MAX_MEM
      SUBROUTINE DMUMPS_SET_MEMESTIMGLOBAL (
     &           OOC_STRAT, BLR_STRAT, BLR_CASE, SUM_OF_PEAKS,
     &           KEEP8, 
     &           SUM_NRLADU_underL0, SUM_NRLADU_if_LR_LU_underL0,
     &           SUM_NRLADULR_UD_underL0, SUM_NRLADULR_WC_underL0,
     &           PEAK_FR,              
     &           PEAK_FR_OOC,          
     &           NRLNEC_if_LR_LU,      
     &           NRLNEC_if_LR_LUCB,    
     &           NRLNECOOC_if_LR_LUCB, 
     &           NRLNEC_if_LR_CB,      
     &           NRLADULR_UD,          
     &           NRLADULR_WC,          
     &           NRLNECLR_CB_UD,       
     &           NRLNECLR_LUCB_UD,     
     &           NRLNECLR_LUCB_WC,     
     &           MemEstimGlobal  
     &           )
      INTEGER, INTENT(IN)    :: OOC_STRAT, BLR_STRAT
      INTEGER, INTENT(IN)    :: BLR_CASE
      LOGICAL, INTENT(IN)    :: SUM_OF_PEAKS
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      INTEGER(8), INTENT(IN) :: SUM_NRLADU_underL0, 
     &              SUM_NRLADU_if_LR_LU_underL0, 
     &              SUM_NRLADULR_UD_underL0, 
     &              SUM_NRLADULR_WC_underL0
      INTEGER(8), INTENT(IN) :: 
     &           PEAK_FR,              
     &           PEAK_FR_OOC,          
     &           NRLNEC_if_LR_LU,      
     &           NRLNEC_if_LR_LUCB,    
     &           NRLNECOOC_if_LR_LUCB, 
     &           NRLNEC_if_LR_CB,      
     &           NRLADULR_UD,          
     &           NRLADULR_WC,          
     &           NRLNECLR_CB_UD,       
     &           NRLNECLR_LUCB_UD,     
     &           NRLNECLR_LUCB_WC
      INTEGER(8), INTENT(OUT) :: MemEstimGlobal
      IF ( OOC_STRAT .GT. 0 .OR. OOC_STRAT .EQ. -1 ) THEN
        MemEstimGlobal = PEAK_FR_OOC   
      ELSE
        MemEstimGlobal = PEAK_FR 
      ENDIF
      IF (BLR_STRAT.GT.0) THEN
        IF (.NOT.SUM_OF_PEAKS) THEN
         IF (BLR_STRAT.EQ.1) THEN
          IF (BLR_CASE.LE.1) THEN
           IF (OOC_STRAT.EQ.0) THEN
            MemEstimGlobal = KEEP8(40)  
           ELSE
            MemEstimGlobal = KEEP8(41) 
           ENDIF
          ELSE IF (BLR_CASE.EQ.2) THEN
           IF (OOC_STRAT.EQ.0) THEN
            MemEstimGlobal = KEEP8(33)  
           ELSE
            MemEstimGlobal = KEEP8(54)  
           ENDIF
          ELSE
           IF (OOC_STRAT.EQ.0) THEN
            MemEstimGlobal = KEEP8(53)  
           ELSE
            MemEstimGlobal = KEEP8(42)  
           ENDIF
          ENDIF
         ELSE 
          IF (BLR_CASE.LE.1) THEN
           IF (OOC_STRAT.EQ.0) THEN
            MemEstimGlobal = KEEP8(43)
           ELSE
            MemEstimGlobal = KEEP8(45)
           ENDIF
          ELSE IF (BLR_CASE.EQ.2) THEN
           IF (OOC_STRAT.EQ.0) THEN
            MemEstimGlobal = KEEP8(34)
           ELSE
            MemEstimGlobal = KEEP8(35)
           ENDIF
          ELSE
           IF (OOC_STRAT.EQ.0) THEN
            MemEstimGlobal = KEEP8(44)
           ELSE
            MemEstimGlobal = KEEP8(46)
           ENDIF
          ENDIF
         ENDIF
        ELSE 
         IF (BLR_STRAT.EQ.1) THEN
          IF (OOC_STRAT.EQ.0) THEN
            MemEstimGlobal =  NRLNEC_if_LR_LU
     &            + SUM_NRLADU_if_LR_LU_underL0
          ELSE
            MemEstimGlobal =  PEAK_FR_OOC
          ENDIF
          IF (BLR_CASE.EQ.1) THEN
           MemEstimGlobal = MemEstimGlobal + NRLADULR_UD
     &                     + SUM_NRLADULR_UD_underL0
          ELSE IF (BLR_CASE.EQ.3) THEN
           MemEstimGlobal = MemEstimGlobal + NRLADULR_WC
     &                     + SUM_NRLADULR_WC_underL0
          ENDIF
         ELSE IF (BLR_STRAT.EQ.2) THEN
          IF (OOC_STRAT.EQ.0) THEN
            MemEstimGlobal = NRLNEC_if_LR_LUCB 
     &                    + SUM_NRLADU_if_LR_LU_underL0
          ELSE
            MemEstimGlobal = NRLNECOOC_if_LR_LUCB
          ENDIF
          IF (BLR_CASE.EQ.1) THEN
           MemEstimGlobal = MemEstimGlobal + NRLNECLR_LUCB_UD
     &         + SUM_NRLADULR_UD_underL0
          ELSE IF (BLR_CASE.EQ.3) THEN
           MemEstimGlobal = MemEstimGlobal + NRLNECLR_LUCB_WC
     &                     + SUM_NRLADULR_WC_underL0
          ENDIF
         ELSE
          IF (OOC_STRAT.EQ.0) THEN
            MemEstimGlobal = NRLNEC_if_LR_CB
     &         + SUM_NRLADU_underL0
          ELSE
            MemEstimGlobal =  NRLNECOOC_if_LR_LUCB
          ENDIF
          MemEstimGlobal = MemEstimGlobal + NRLNECLR_CB_UD
         ENDIF
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SET_MEMESTIMGLOBAL
      SUBROUTINE DMUMPS_SET_BLRSTRAT_AND_MAXS_K8 ( 
     &           MAXS_BASE8, MAXS_BASE_RELAXED8, 
     &           BLR_STRAT,
     &           KEEP, KEEP8)
      IMPLICIT NONE
      INTEGER(8), INTENT(OUT) :: MAXS_BASE8, MAXS_BASE_RELAXED8
      INTEGER, INTENT(OUT) :: BLR_STRAT
      INTEGER, INTENT(IN)    :: KEEP(500)
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      CALL DMUMPS_SET_BLRSTRAT_AND_MAXS (
     &           MAXS_BASE8, MAXS_BASE_RELAXED8,
     &           BLR_STRAT,
     &           KEEP(1), 
     &           KEEP8(12),
     &           KEEP8(14),
     &           KEEP8(33),
     &           KEEP8(34),
     &           KEEP8(35),
     &           KEEP8(50) )
      RETURN
      END SUBROUTINE DMUMPS_SET_BLRSTRAT_AND_MAXS_K8
      SUBROUTINE DMUMPS_SET_BLRSTRAT_AND_MAXS(
     &           MAXS_BASE8, MAXS_BASE_RELAXED8,
     &           BLR_STRAT, KEEP, 
     &           NRLNEC,
     &           NRLNEC_ACTIVE,
     &           NRLNEC_if_LR_LU,
     &           NRLNEC_if_LR_LUCB,
     &           NRLNECOOC_if_LR_LUCB,
     &           NRLNEC_if_LR_CB )
      IMPLICIT NONE
      INTEGER(8), INTENT(OUT) :: MAXS_BASE8, MAXS_BASE_RELAXED8
      INTEGER,    INTENT(OUT) :: BLR_STRAT
      INTEGER,    INTENT(IN)  :: KEEP(500)
      INTEGER(8), INTENT(IN)  :: NRLNEC,
     &                           NRLNEC_ACTIVE,
     &                           NRLNEC_if_LR_LU,
     &                           NRLNEC_if_LR_LUCB,
     &                           NRLNECOOC_if_LR_LUCB,
     &                           NRLNEC_if_LR_CB
      INTEGER :: PERLU
      PERLU = KEEP(12)
      IF (KEEP(201) .EQ. 0) THEN
        MAXS_BASE8 = NRLNEC
       ELSE
        MAXS_BASE8 = NRLNEC_ACTIVE
      ENDIF
      BLR_STRAT = 0
      IF (KEEP(486).EQ.2) THEN
        IF (KEEP(489).EQ.1) THEN
          BLR_STRAT  = 2
          IF (KEEP(201).NE.0) THEN
           MAXS_BASE8 = NRLNECOOC_if_LR_LUCB
          ELSE
           MAXS_BASE8 = NRLNEC_if_LR_LUCB
          ENDIF
        ELSE
          BLR_STRAT  = 1
          IF (KEEP(201).NE.0) THEN
            MAXS_BASE8 = NRLNEC_ACTIVE
          ELSE
            MAXS_BASE8 = NRLNEC_if_LR_LU
          ENDIF
        ENDIF
      ELSE IF (KEEP(486).EQ.3) THEN
        IF (KEEP(489).EQ.1) THEN
          BLR_STRAT  = 3
          IF (KEEP(201).NE.0) THEN
           MAXS_BASE8 = NRLNECOOC_if_LR_LUCB
          ELSE
           MAXS_BASE8 = NRLNEC_if_LR_CB
          ENDIF
        ENDIF
      ENDIF
      IF ( MAXS_BASE8 .GT. 0_8 ) THEN
          MAXS_BASE_RELAXED8 =
     &         MAXS_BASE8 + int(PERLU,8) * ( MAXS_BASE8 / 100_8 + 1_8)
          MAXS_BASE_RELAXED8 = max(MAXS_BASE_RELAXED8, 1_8)
      ELSE
        MAXS_BASE_RELAXED8 = 1_8
      END IF
      RETURN
      END SUBROUTINE DMUMPS_SET_BLRSTRAT_AND_MAXS
      SUBROUTINE DMUMPS_MEM_ALLOWED_SET_MAXS ( MAXS,
     &           BLR_STRAT, OOC_STRAT, MAXS_ESTIM_RELAXED8,
     &           KEEP, KEEP8, MYID, N, NELT, NA, LNA, 
     &           NSLAVES, ICNTL38, ICNTL39, IFLAG, IERROR
     &           )
      IMPLICIT NONE
      INTEGER(8), INTENT(OUT)    :: MAXS
      INTEGER, INTENT(INOUT)     :: IFLAG, IERROR
      INTEGER, INTENT(IN)    :: BLR_STRAT
      INTEGER, INTENT(IN)    :: OOC_STRAT
      INTEGER(8), INTENT(IN) :: MAXS_ESTIM_RELAXED8
      INTEGER, INTENT(IN)       :: KEEP(500)
      INTEGER(8), INTENT(INOUT) :: KEEP8(150)
      INTEGER, INTENT(IN)    :: MYID, N, NELT, NSLAVES, LNA
      INTEGER, INTENT(IN)    :: NA(LNA), ICNTL38, ICNTL39
      INTEGER(8) :: SMALLER_MAXS, UPDATED_DIFF
      LOGICAL    :: EFF, PERLU_ON, SUM_OF_PEAKS
      INTEGER    :: BLR_CASE
      INTEGER(8) :: TOTAL_BYTES, MEM_ALLOWED_BYTES, 
     &              MEM_DISPO_BYTES, MEM_DISPO
      INTEGER    :: TOTAL_MBYTES, PERLU
      INTEGER(8) :: MEM_DISPO_BYTES_NR, MEM_DISPO_NR, 
     &              TOTAL_BYTES_NR
      INTEGER    :: TOTAL_MBYTES_NR
      INTEGER, PARAMETER :: IDUMMY = -9999
      LOGICAL, PARAMETER :: BDUMMY =.FALSE.
      PERLU_ON     = .TRUE.
      PERLU        = KEEP(12)
      EFF          = .FALSE.
      SUM_OF_PEAKS = .TRUE.
      BLR_CASE     = 1          
      MEM_ALLOWED_BYTES  = KEEP8(4)  
      CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        MYID, N, NELT, NA(1), LNA,
     &        KEEP8(28), KEEP8(30),
     &        NSLAVES, TOTAL_MBYTES, EFF, OOC_STRAT,
     &        BLR_STRAT, PERLU_ON, TOTAL_BYTES, 
     &        BLR_CASE, SUM_OF_PEAKS, .FALSE. , 
     &        .FALSE. 
     &         )
      MEM_DISPO_BYTES = MEM_ALLOWED_BYTES-TOTAL_BYTES
      MEM_DISPO       = MEM_DISPO_BYTES/int(KEEP(35),8)
      IF (BLR_STRAT.EQ.0) THEN
           UPDATED_DIFF = 0_8
      ELSE IF (BLR_STRAT.EQ.1) THEN
           IF (KEEP(464).NE.0) THEN
            UPDATED_DIFF = int(
     &                          dble(KEEP8(36)) * ( 1.0D0 - 
     &                          dble(ICNTL38)/dble(KEEP(464)) )
     &                     , 8)
           ELSE
            UPDATED_DIFF = int (
     &                          -dble(KEEP8(11)-KEEP8(32))  * 
     &                           dble(ICNTL38) / 1000.0D0
     &                     , 8)
           ENDIF
      ELSE IF (BLR_STRAT.EQ.2) THEN
           IF (KEEP(464)+KEEP(465).NE.0) THEN
            UPDATED_DIFF = int(
     &                          dble(KEEP8(38)) * ( 1.0D0 -
     &                          dble(ICNTL38+ICNTL39)/
     &                          dble(KEEP(464)+KEEP(465)) ) 
     &                     , 8)
           ELSE
            UPDATED_DIFF = int(
     &                         -dble(KEEP8(39))*
     &                          dble(ICNTL38+ICNTL39)/1000.0D0
     &                     , 8)
           ENDIF
      ELSE IF (BLR_STRAT.EQ.3) THEN
           IF (KEEP(465).NE.0) THEN
            UPDATED_DIFF =  int(
     &                          dble(KEEP8(37)) * ( 1.0D0 -
     &                          dble(ICNTL39)/dble(KEEP(465)) )
     &                      , 8)
           ELSE
            UPDATED_DIFF = int(
     &                         -dble(KEEP8(39))*
     &                          dble(ICNTL39)/1000.0D0
     &                      , 8)
           ENDIF
      ELSE
        UPDATED_DIFF = 0_8
      ENDIF
      MEM_DISPO   = MEM_DISPO + UPDATED_DIFF
      MAXS        = MAXS_ESTIM_RELAXED8
      MEM_DISPO_NR = 0_8
      IF ( (MEM_DISPO.LT.0) .AND. MAXS_ESTIM_RELAXED8.GT. 
     &    (KEEP8(4)/int(KEEP(35),8)) ) THEN
         PERLU_ON     = .FALSE.
             CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        MYID, N, NELT, NA(1), LNA,
     &        KEEP8(28), KEEP8(30),
     &        NSLAVES, TOTAL_MBYTES_NR, EFF, OOC_STRAT,
     &        BLR_STRAT, PERLU_ON, TOTAL_BYTES_NR, 
     &        BLR_CASE, SUM_OF_PEAKS, .FALSE., 
     &        .FALSE. 
     &         )
         MEM_DISPO_BYTES_NR = MEM_ALLOWED_BYTES-TOTAL_BYTES_NR
         MEM_DISPO_NR       = 
     &            MEM_DISPO_BYTES_NR/int(KEEP(35),8)
     &             + UPDATED_DIFF
        IF ( MEM_DISPO_NR.LT.0 ) THEN
          IFLAG=-19
          CALL MUMPS_SET_IERROR(-MEM_DISPO_NR,IERROR)
             GOTO 500
        ELSE
          IF (BLR_STRAT.GE.2) THEN
            IFLAG=-19
            CALL MUMPS_SET_IERROR(-MEM_DISPO_NR,IERROR)
            GOTO 500
          ELSE
            MEM_DISPO_NR = MEM_DISPO_NR - 
     &                    (int(KEEP(12),8)/120_8)*
     &                    (KEEP8(11)/4_8)
            IF ( MEM_DISPO_NR.LT.0 ) THEN
              IFLAG=-19
              CALL MUMPS_SET_IERROR(-MEM_DISPO_NR,IERROR)
              GOTO 500
            ELSE
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      MAXS = MAXS_ESTIM_RELAXED8 
      IF (BLR_STRAT.EQ.0) THEN
           IF (MEM_DISPO.GT.0) THEN
            IF (OOC_STRAT.EQ.0) THEN
              MAXS = MAXS_ESTIM_RELAXED8+(MEM_DISPO/2_8)
            ELSE
              MAXS = MAXS_ESTIM_RELAXED8+(MEM_DISPO/2_8)
            ENDIF
           ELSE
              MAXS = MAXS_ESTIM_RELAXED8 + MEM_DISPO
           ENDIF
      ELSE IF (BLR_STRAT.EQ.1)  THEN
          IF ( MEM_DISPO .LT. 0) THEN
           IF (OOC_STRAT.EQ.0) THEN
             SMALLER_MAXS = KEEP8(34) + 
     &          int(PERLU,8) * ( KEEP8(34) / 100_8 + 1_8)
           ELSE
             SMALLER_MAXS = KEEP8(35) + 
     &          int(PERLU,8) * ( KEEP8(35) / 100_8 + 1_8)
           ENDIF
           MAXS = max(MAXS_ESTIM_RELAXED8+MEM_DISPO,
     &            SMALLER_MAXS)
          ENDIF
      ELSE IF (BLR_STRAT.EQ.2)  THEN
          IF ( MEM_DISPO.LT.0) THEN
           MAXS = max(
     &            MAXS_ESTIM_RELAXED8+MEM_DISPO,
     &            MAXS_ESTIM_RELAXED8)
          ENDIF
      ELSE IF (BLR_STRAT.EQ.3)  THEN
          IF ( MEM_DISPO.LT.0) THEN
           MAXS = max(
     &            MAXS_ESTIM_RELAXED8+MEM_DISPO,
     &            MAXS_ESTIM_RELAXED8)
          ENDIF
      ENDIF
      IF (MAXS .LE. 0_8) THEN
         IFLAG=-19
         IF (MEM_DISPO.LT.0) THEN
          CALL MUMPS_SET_IERROR(MEM_DISPO,IERROR)
         ELSE
          CALL MUMPS_SET_IERROR(MAXS_ESTIM_RELAXED8-MAXS,IERROR)
         ENDIF
      ENDIF
       CALL  DMUMPS_MEM_ALLOWED_SET_K75 ( 
     &           MAXS, MYID, 
     &           .FALSE., 
     &           N, NELT, NA, LNA, NSLAVES,
     &           BLR_STRAT, OOC_STRAT, 
     &           KEEP, KEEP8, IFLAG, IERROR
     &           )
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_MEM_ALLOWED_SET_MAXS
      SUBROUTINE DMUMPS_MEM_ALLOWED_SET_K75 ( 
     &           MAXS, MYID, UNDER_L0_OMP, 
     &           N, NELT, NA, LNA, NSLAVES,
     &           BLR_STRAT, OOC_STRAT, 
     &           KEEP, KEEP8, IFLAG, IERROR
     &           )
      IMPLICIT NONE
      INTEGER(8), INTENT(IN)    :: MAXS
      INTEGER, INTENT(IN)       :: MYID, N, NELT, LNA, NSLAVES, 
     &                             BLR_STRAT, OOC_STRAT
      LOGICAL, INTENT(IN)       :: UNDER_L0_OMP
      INTEGER, INTENT(IN)       :: NA(LNA), KEEP(500)
      INTEGER(8), INTENT(INOUT) :: KEEP8(150)
      INTEGER, INTENT(INOUT)    :: IFLAG, IERROR
      INTEGER(8) :: KEEP8_23_SAVETMP, TOTAL_BYTES
      INTEGER    :: TOTAL_MBYTES
      LOGICAL    :: PERLU_ON, MEM_EFF_ALLOCATED, EFF
      INTEGER, PARAMETER :: IDUMMY = -9999
      LOGICAL, PARAMETER :: BDUMMY =.FALSE.
      KEEP8_23_SAVETMP =  KEEP8(23)
      KEEP8(23) = MAXS
      PERLU_ON          =.TRUE.
      MEM_EFF_ALLOCATED = .TRUE.
      EFF               = .TRUE.
      KEEP8(74) = 0_8 
      KEEP8(63) = 0_8
      CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &    MYID, N, NELT, NA(1), LNA, KEEP8(28),
     &    KEEP8(30),
     &    NSLAVES, TOTAL_MBYTES, EFF , KEEP(201),
     &    BLR_STRAT, PERLU_ON, TOTAL_BYTES,
     &    IDUMMY, BDUMMY , MEM_EFF_ALLOCATED,
     &    UNDER_L0_OMP
     &    )
      KEEP8(23) = KEEP8_23_SAVETMP
      KEEP8(75) = KEEP8(4) - TOTAL_BYTES
      KEEP8(75) =  KEEP8(75)/int(KEEP(35),8)
      IF (KEEP8(75).LT.0_8) THEN
          IFLAG=-19
          CALL MUMPS_SET_IERROR(-KEEP8(75),IERROR)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_MEM_ALLOWED_SET_K75
      SUBROUTINE DMUMPS_SETMAXTOZERO(M_ARRAY, M_SIZE)
      IMPLICIT NONE
      INTEGER M_SIZE
      DOUBLE PRECISION M_ARRAY(M_SIZE)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      M_ARRAY=ZERO
      RETURN
      END SUBROUTINE DMUMPS_SETMAXTOZERO
      SUBROUTINE DMUMPS_COMPUTE_NBROWSinF (
     &      N, INODE, IFATH, KEEP, 
     &      IOLDPS, HF, IW, LIW, 
     &      NROWS, NCOLS, NPIV,
     &      NELIM, NFS4FATHER,
     &      NBROWSinF
     &  )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, INODE, IFATH
      INTEGER, INTENT(IN) :: KEEP(500)
      INTEGER, INTENT(IN) :: IOLDPS, HF, LIW, NROWS, NCOLS
      INTEGER, INTENT(IN) :: NPIV, NELIM, NFS4FATHER
      INTEGER, INTENT(IN) :: IW(LIW)
      INTEGER, INTENT(OUT)::  NBROWSinF
      INTEGER ::   ShiftFirstRowinFront
      NBROWSinF = 0
      IF ( (KEEP(219).EQ.0).OR.(KEEP(50).NE.2).OR.
     &     (NFS4FATHER.LE.0) ) THEN
       RETURN
      ENDIF
      ShiftFirstRowinFront = NCOLS-NPIV-NELIM-NROWS
      IF (ShiftFirstRowinFront.EQ.0) THEN
       NBROWSinF = min(NROWS, NFS4FATHER-NELIM)
      ELSE IF (ShiftFirstRowinFront.LT.NFS4FATHER-NELIM) THEN
       NBROWSinF = min(NROWS,NFS4FATHER-NELIM-ShiftFirstRowinFront)
      ELSE
       NBROWSinF=0
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_COMPUTE_NBROWSinF 
      SUBROUTINE DMUMPS_COMPUTE_ESTIM_NFS4FATHER (
     &      N, INODE, IFATH, FILS, PERM, KEEP, 
     &      IOLDPS, HF, IW, LIW, NFRONT, NASS1,
     &      ESTIM_NFS4FATHER_ATSON
     &  )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, INODE, IFATH
      INTEGER, INTENT(IN) :: FILS(N), PERM(N), KEEP(500)
      INTEGER, INTENT(IN) :: IOLDPS, HF, LIW, NFRONT, NASS1
      INTEGER, INTENT(IN) :: IW(LIW)
      INTEGER, INTENT(OUT)::  ESTIM_NFS4FATHER_ATSON
      INTEGER :: J, J_LASTFS, IN, NCB, I, IPOS
      ESTIM_NFS4FATHER_ATSON = 0
      IN = IFATH
      J_LASTFS = IN
      DO WHILE (IN.GT.0)
          J_LASTFS = IN
          IN = FILS(IN)
      ENDDO
      NCB   = NFRONT-NASS1
      IPOS  = IOLDPS + HF + NASS1 
      ESTIM_NFS4FATHER_ATSON = 0
      DO I=1, NCB
        J = IW(IPOS+ESTIM_NFS4FATHER_ATSON)
        IF (PERM(J).LE.PERM(J_LASTFS)) THEN
         ESTIM_NFS4FATHER_ATSON = 
     &        ESTIM_NFS4FATHER_ATSON+1
        ELSE
         EXIT
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_COMPUTE_ESTIM_NFS4FATHER
      SUBROUTINE DMUMPS_COMPUTE_MAXPERCOL(
     &     A,ASIZE,NCOL,NROW,
     &     M_ARRAY,NMAX,PACKED_CB,LROW1)
      IMPLICIT NONE
      INTEGER(8) :: ASIZE
      INTEGER NROW,NCOL,NMAX,LROW1
      LOGICAL PACKED_CB
      DOUBLE PRECISION A(ASIZE)
      DOUBLE PRECISION M_ARRAY(NMAX)
      INTEGER I
      INTEGER(8):: APOS, J, LROW
      DOUBLE PRECISION ZERO,TMP
      PARAMETER (ZERO=0.0D0)
      DO I=1, NMAX
          M_ARRAY(I) = ZERO
      ENDDO
      APOS = 0_8
      IF (PACKED_CB) THEN
        LROW=int(LROW1,8)
      ELSE
        LROW=int(NCOL,8)
      ENDIF
      DO I=1,NROW
         DO J=1_8,int(NMAX,8)
            TMP = abs(A(APOS+J))
            IF(TMP.GT.M_ARRAY(J)) M_ARRAY(J) = TMP
         ENDDO
         APOS = APOS + LROW
         IF (PACKED_CB) LROW=LROW+1_8
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_COMPUTE_MAXPERCOL
      SUBROUTINE DMUMPS_SIZE_IN_STRUCT( id, NB_INT, NB_CMPLX, NB_CHAR )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC) :: id
      INTEGER(8) NB_INT, NB_CMPLX 
      INTEGER(8) NB_REAL,NB_CHAR
      NB_INT        = 0_8
      NB_CMPLX      = 0_8
      NB_REAL       = 0_8
      NB_CHAR       = 0_8
      IF (associated(id%IS))          NB_INT=NB_INT+size(id%IS)
      NB_INT=NB_INT+size(id%KEEP)
      NB_INT=NB_INT+size(id%ICNTL)
      NB_INT=NB_INT+size(id%INFO)
      NB_INT=NB_INT+size(id%INFOG)
      IF (associated(id%MAPPING))     NB_INT=NB_INT+size(id%MAPPING)
      IF (associated(id%STEP))        NB_INT=NB_INT+size(id%STEP)
      IF (associated(id%NE_STEPS  ))  NB_INT=NB_INT+size(id%NE_STEPS  )
      IF (associated(id%ND_STEPS))    NB_INT=NB_INT+size(id%ND_STEPS)
      IF (associated(id%Step2node))   NB_INT=NB_INT+size(id%Step2node)
      IF (associated(id%FRERE_STEPS)) NB_INT=NB_INT+size(id%FRERE_STEPS)
      IF (associated(id%DAD_STEPS))   NB_INT=NB_INT+size(id%DAD_STEPS)
      IF (associated(id%FILS))        NB_INT=NB_INT+size(id%FILS)
      IF (associated(id%PTRAR))
     &     NB_INT=NB_INT+size(id%PTRAR)* id%KEEP(10)
      IF (associated(id%FRTPTR))      NB_INT=NB_INT+size(id%FRTPTR)
      NB_INT=NB_INT+size(id%KEEP8) * id%KEEP(10)
      IF (associated(id%PTRFAC)) NB_INT=NB_INT+size(id%PTRFAC) *
     &                                         id%KEEP(10)
      IF (associated(id%FRTELT))      NB_INT=NB_INT+size(id%FRTELT)
      IF (associated(id%NA))          NB_INT=NB_INT+size(id%NA)
      IF       (associated(id%PROCNODE_STEPS))
     &  NB_INT=NB_INT+size(id%PROCNODE_STEPS)
      IF (associated(id%PTLUST_S)) NB_INT=NB_INT+size(id%PTLUST_S)
      IF (associated(id%INTARR)) NB_INT=NB_INT+id%KEEP8(27)
      IF (associated(id%ELTPROC))  NB_INT=NB_INT+size(id%ELTPROC)
      IF (associated(id%CANDIDATES))
     &     NB_INT=NB_INT+size(id%CANDIDATES)
      IF (associated(id%SYM_PERM))
     &     NB_INT=NB_INT+size(id%SYM_PERM)
      IF (associated(id%UNS_PERM))
     &  NB_INT=NB_INT+size(id%UNS_PERM)
      IF       (associated(id%ISTEP_TO_INIV2))
     &  NB_INT=NB_INT+size(id%ISTEP_TO_INIV2)
      IF       (associated(id%FUTURE_NIV2))
     &  NB_INT=NB_INT+size(id%FUTURE_NIV2)
      IF (associated(id%TAB_POS_IN_PERE))
     &  NB_INT=NB_INT+size(id%TAB_POS_IN_PERE)
      IF (associated(id%I_AM_CAND))
     &  NB_INT=NB_INT+size(id%I_AM_CAND)
      IF (associated(id%MEM_DIST)) 
     &  NB_INT=NB_INT+size(id%MEM_DIST)
      IF (associated(id%POSINRHSCOMP_ROW))
     &     NB_INT=NB_INT+size(id%POSINRHSCOMP_ROW)
      IF(id%POSINRHSCOMP_COL_ALLOC.AND.associated(id%POSINRHSCOMP_COL))
     &     NB_INT=NB_INT+size(id%POSINRHSCOMP_COL)
      IF       (associated(id%MEM_SUBTREE))
     &  NB_REAL=NB_REAL+size(id%MEM_SUBTREE)*(id%KEEP(35)/id%KEEP(16))
      IF       (associated(id%MY_ROOT_SBTR))
     &  NB_INT=NB_INT+size(id%MY_ROOT_SBTR)
      IF       (associated(id%MY_FIRST_LEAF))
     &  NB_INT=NB_INT+size(id%MY_FIRST_LEAF)
      IF (associated(id%MY_NB_LEAF)) NB_INT=NB_INT+size(id%MY_NB_LEAF)
      IF (associated(id%DEPTH_FIRST)) NB_INT=NB_INT+size(id%DEPTH_FIRST)
      IF (associated(id%DEPTH_FIRST_SEQ))
     &     NB_INT=NB_INT+size(id%DEPTH_FIRST_SEQ)
      IF (associated(id%SBTR_ID)) NB_INT=NB_INT+size(id%SBTR_ID)
      IF (associated(id%SCHED_DEP)) NB_INT=NB_INT+size(id%SCHED_DEP)
      IF (associated(id%SCHED_GRP)) NB_INT=NB_INT+size(id%SCHED_GRP)
      IF (associated(id%SCHED_SBTR)) NB_INT=NB_INT+size(id%SCHED_SBTR)
      IF (associated(id%CROIX_MANU)) NB_INT=NB_INT+size(id%CROIX_MANU)
      IF (associated(id%COST_TRAV))
     &     NB_REAL=NB_REAL+size(id%COST_TRAV)*(id%KEEP(35)/id%KEEP(16))
      IF (associated(id%CB_SON_SIZE)) NB_INT=NB_INT+size(id%CB_SON_SIZE)
      IF       (associated(id%OOC_INODE_SEQUENCE))
     &  NB_INT=NB_INT+size(id%OOC_INODE_SEQUENCE)
      IF       (associated(id%OOC_SIZE_OF_BLOCK))
     &  NB_INT=NB_INT+size(id%OOC_SIZE_OF_BLOCK)*id%KEEP(10)
      IF       (associated(id%OOC_VADDR)) 
     &  NB_INT=NB_INT+size(id%OOC_VADDR)*id%KEEP(10)
      IF       (associated(id%OOC_TOTAL_NB_NODES))
     &  NB_INT=NB_INT+size(id%OOC_TOTAL_NB_NODES)
      IF       (associated(id%OOC_NB_FILES))
     &  NB_INT=NB_INT+size(id%OOC_NB_FILES)
      IF       (associated(id%OOC_FILE_NAME_LENGTH))
     &  NB_INT=NB_INT+size(id%OOC_FILE_NAME_LENGTH)
      IF (associated(id%PIVNUL_LIST)) NB_INT=NB_INT+size(id%PIVNUL_LIST)
      IF (associated(id%SUP_PROC))    NB_INT=NB_INT+size(id%SUP_PROC)
      IF (associated(id%IPTR_WORKING))
     &     NB_INT=NB_INT+size(id%IPTR_WORKING)
      IF (associated(id%WORKING)) NB_INT=NB_INT+size(id%WORKING)
      IF (associated(id%LRGROUPS))
     &     NB_INT=NB_INT+size(id%LRGROUPS)
      IF (associated(id%IPOOL_B_L0_OMP))
     &     NB_INT=NB_INT+size(id%IPOOL_B_L0_OMP)
      IF (associated(id%IPOOL_A_L0_OMP))
     &     NB_INT=NB_INT+size(id%IPOOL_A_L0_OMP)
      IF (associated(id%PHYS_L0_OMP))
     &     NB_INT=NB_INT+size(id%PHYS_L0_OMP)
      IF (associated(id%VIRT_L0_OMP))
     &     NB_INT=NB_INT+size(id%VIRT_L0_OMP)
      IF (associated(id%PERM_L0_OMP))
     &     NB_INT=NB_INT+size(id%PERM_L0_OMP)
      IF (associated(id%PTR_LEAFS_L0_OMP))
     &     NB_INT=NB_INT+size(id%PTR_LEAFS_L0_OMP)
      IF (associated(id%L0_OMP_MAPPING))
     &     NB_INT=NB_INT+size(id%L0_OMP_MAPPING)
      IF (associated(id%SINGULAR_VALUES))
     &     NB_REAL=NB_REAL+size(id%SINGULAR_VALUES)
      IF (associated(id%root%RG2L_COL))
     &     NB_INT=NB_INT+size(id%root%RG2L_COL)
      IF (associated(id%root%RG2L_ROW))
     &     NB_INT=NB_INT+size(id%root%RG2L_ROW)
      IF (associated(id%root%IPIV))
     &     NB_INT=NB_INT+size(id%root%IPIV)
      IF (associated(id%root%RHS_CNTR_MASTER_ROOT))
     &     NB_CMPLX=NB_CMPLX+size(id%root%RHS_CNTR_MASTER_ROOT)
      IF (associated(id%root%SCHUR_POINTER))
     &     NB_CMPLX=NB_CMPLX+size(id%root%SCHUR_POINTER)
      IF (associated(id%root%QR_TAU))
     &     NB_CMPLX=NB_CMPLX+size(id%root%QR_TAU)
      IF (associated(id%root%RHS_ROOT))
     &     NB_CMPLX=NB_CMPLX+size(id%root%RHS_ROOT)
      IF (associated(id%root%SVD_U))
     &     NB_CMPLX=NB_CMPLX+size(id%root%SVD_U)
      IF (associated(id%root%SVD_VT))
     &     NB_CMPLX=NB_CMPLX+size(id%root%SVD_VT)
      IF (associated(id%root%SINGULAR_VALUES))
     &     NB_REAL=NB_REAL+size(id%root%SINGULAR_VALUES)
      IF (associated(id%DBLARR))  NB_CMPLX=NB_CMPLX+id%KEEP8(26)
      IF (associated(id%RHSCOMP)) NB_CMPLX = NB_CMPLX + id%KEEP8(25)
      IF (associated(id%S))       NB_CMPLX=NB_CMPLX+id%KEEP8(23)
      IF (associated(id%COLSCA).AND.(id%KEEP(52).NE.-1))
     &     NB_REAL=NB_REAL+size(id%COLSCA)
      IF (associated(id%ROWSCA).AND.(id%KEEP(52).NE.-1))
     &     NB_REAL=NB_REAL+size(id%ROWSCA)
      NB_REAL=NB_REAL+size(id%CNTL)
      NB_REAL=NB_REAL+size(id%RINFO)
      NB_REAL=NB_REAL+size(id%RINFOG)
      NB_REAL=NB_REAL+size(id%DKEEP)
      NB_CHAR=NB_CHAR+len(id%VERSION_NUMBER)
      NB_CHAR=NB_CHAR+len(id%OOC_TMPDIR)
      NB_CHAR=NB_CHAR+len(id%OOC_PREFIX)
      NB_CHAR=NB_CHAR+len(id%WRITE_PROBLEM)
      NB_CHAR=NB_CHAR+len(id%SAVE_DIR)
      NB_CHAR=NB_CHAR+len(id%SAVE_PREFIX)
      NB_CMPLX = NB_CMPLX + NB_REAL
      NB_CMPLX = NB_CMPLX + id%KEEP8(71) + id%KEEP8(64)
      RETURN
      END SUBROUTINE DMUMPS_SIZE_IN_STRUCT 
      SUBROUTINE DMUMPS_COPYI8SIZE(N8,SRC,DEST)
      IMPLICIT NONE
      INTEGER(8) :: N8
      DOUBLE PRECISION, intent(in)  :: SRC(N8)
      DOUBLE PRECISION, intent(out) :: DEST(N8)
      INTEGER(8) :: SHIFT8, HUG8
      INTEGER    :: I, I4SIZE
      IF(int(huge(I4SIZE),8) .EQ. int(huge(HUG8),8)) THEN
         CALL dcopy(N8, SRC(1), 1, DEST(1), 1)
      ELSE
         HUG8=int(huge(I4SIZE),8)
         DO I = 1, int(( N8 + HUG8 - 1_8 ) / HUG8)
            SHIFT8 = 1_8 + int(I-1,8) * HUG8
            I4SIZE = int(min(HUG8, N8-SHIFT8+1_8))
            CALL dcopy(I4SIZE, SRC(SHIFT8), 1, DEST(SHIFT8), 1)
         ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_COPYI8SIZE
      SUBROUTINE DMUMPS_SET_TMP_PTR( THE_ADDRESS, THE_SIZE8 )
      USE DMUMPS_STATIC_PTR_M
      INTEGER(8), INTENT(IN) :: THE_SIZE8
      DOUBLE PRECISION,    INTENT(IN) :: THE_ADDRESS(THE_SIZE8)
      CALL DMUMPS_SET_STATIC_PTR(THE_ADDRESS(1:THE_SIZE8)) 
      RETURN
      END SUBROUTINE DMUMPS_SET_TMP_PTR
      SUBROUTINE DMUMPS_OOC_IO_LU_PANEL_I
     &     ( STRAT, TYPEFile, 
     &     AFAC, LAFAC, MonBloc,
     &     LNextPiv2beWritten, UNextPiv2beWritten,
     &     IW, LIWFAC,
     &     MYID, FILESIZE, IERR , LAST_CALL)
      USE DMUMPS_OOC, ONLY : IO_BLOCK,
     &                       DMUMPS_OOC_IO_LU_PANEL
      IMPLICIT NONE
      TYPE(IO_BLOCK), INTENT(INOUT):: MonBloc
      INTEGER(8) :: LAFAC
      INTEGER, INTENT(IN)   :: STRAT, LIWFAC, MYID, TYPEFile
      INTEGER, INTENT(INOUT)        :: IW(0:LIWFAC-1) 
      DOUBLE PRECISION, INTENT(IN) :: AFAC(LAFAC)
      INTEGER,   INTENT(INOUT) :: LNextPiv2beWritten, UNextPiv2beWritten
      INTEGER(8), INTENT(INOUT) :: FILESIZE
      INTEGER,   INTENT(OUT) :: IERR
      LOGICAL,   INTENT(IN)  :: LAST_CALL
      CALL DMUMPS_OOC_IO_LU_PANEL 
     &     ( STRAT, TYPEFile, 
     &     AFAC, LAFAC, MonBloc,
     &     LNextPiv2beWritten, UNextPiv2beWritten,
     &     IW, LIWFAC,
     &     MYID, FILESIZE, IERR , LAST_CALL)
      RETURN
      END SUBROUTINE DMUMPS_OOC_IO_LU_PANEL_I
      SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE3_I ( N, ISON,
     &             NBCOL_SON, NBROW_SON, INDCOL_SON, INDROW_SON,
     &             LD_SON, VAL_SON, TAG, SUBSET_ROW, SUBSET_COL,
     &             NSUBSET_ROW, NSUBSET_COL,
     &             NSUPROW, NSUPCOL,
     &             NPROW, NPCOL, MBLOCK, RG2L_ROW, RG2L_COL,
     &             NBLOCK, PDEST, COMM, IERR , 
     &             TAB, TABSIZE, TRANSP, SIZE_PACK,
     &             N_ALREADY_SENT, KEEP, BBPCBP ) 
      USE DMUMPS_BUF, ONLY : DMUMPS_BUF_SEND_CONTRIB_TYPE3
      IMPLICIT NONE
      INTEGER N, ISON, NBCOL_SON, NBROW_SON, NSUBSET_ROW, NSUBSET_COL
      INTEGER NPROW, NPCOL, MBLOCK, NBLOCK, LD_SON
      INTEGER BBPCBP
      INTEGER PDEST, TAG, COMM, IERR
      INTEGER INDCOL_SON( NBCOL_SON ), INDROW_SON( NBROW_SON )
      INTEGER SUBSET_ROW( NSUBSET_ROW ), SUBSET_COL( NSUBSET_COL )
      INTEGER :: RG2L_ROW(N)
      INTEGER :: RG2L_COL(N)
      INTEGER NSUPROW, NSUPCOL
      INTEGER(8), INTENT(IN) :: TABSIZE
      INTEGER SIZE_PACK
      INTEGER KEEP(500)
      DOUBLE PRECISION VAL_SON( LD_SON, * ), TAB(*)
      LOGICAL TRANSP
      INTEGER N_ALREADY_SENT
      CALL DMUMPS_BUF_SEND_CONTRIB_TYPE3( N, ISON,
     &             NBCOL_SON, NBROW_SON, INDCOL_SON, INDROW_SON,
     &             LD_SON, VAL_SON, TAG, SUBSET_ROW, SUBSET_COL,
     &             NSUBSET_ROW, NSUBSET_COL,
     &             NSUPROW, NSUPCOL,
     &             NPROW, NPCOL, MBLOCK, RG2L_ROW, RG2L_COL,
     &             NBLOCK, PDEST, COMM, IERR , 
     &             TAB, TABSIZE, TRANSP, SIZE_PACK,
     &             N_ALREADY_SENT, KEEP, BBPCBP ) 
      RETURN
      END SUBROUTINE DMUMPS_BUF_SEND_CONTRIB_TYPE3_I
      SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING_I(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR_L, sizeBEGS_BLR_L,
     &        BEGS_BLR_U, sizeBEGS_BLR_U, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        BLR_U,
     &        NB_BLR_U, NELIM, LBANDSLAVE, ISHIFT, NIV, SYM,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT)
      USE DMUMPS_LR_TYPE, ONLY : LRB_TYPE
      USE DMUMPS_FAC_LR, ONLY : DMUMPS_BLR_UPDATE_TRAILING
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR_L, NB_BLR_U, 
     &                                CURRENT_BLR,
     &                                NELIM, NIV, SYM, TOL_OPT
      INTEGER, intent(inout)         :: IFLAG, IERROR
      LOGICAL, intent(in)          :: LBANDSLAVE
      INTEGER, intent(in)          :: ISHIFT
      DOUBLE PRECISION, intent(inout) :: A(LA)
      TYPE(LRB_TYPE),intent(in) :: BLR_U(NB_BLR_U-CURRENT_BLR)
      TYPE(LRB_TYPE),intent(in) :: BLR_L(NB_BLR_L-CURRENT_BLR)
      INTEGER :: sizeBEGS_BLR_L, sizeBEGS_BLR_U
      INTEGER :: BEGS_BLR_L(sizeBEGS_BLR_L)
      INTEGER ::  BEGS_BLR_U(sizeBEGS_BLR_U)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      CALL DMUMPS_BLR_UPDATE_TRAILING(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR_L, BEGS_BLR_U, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        BLR_U,
     &        NB_BLR_U, NELIM, LBANDSLAVE, ISHIFT, NIV, SYM,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT)
      RETURN
      END SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING_I
      SUBROUTINE DMUMPS_COMPRESS_CB_I(A_PTR, LA_PTR, POSELT, LDA,
     &        BEGS_BLR, sizeBEGS_BLR, BEGS_BLR_U, sizeBEGS_BLR_U,
     &        NB_ROWS, NB_COLS, NB_INASM,
     &        NROWS, NCOLS, INODE,   
     &        IWHANDLER, SYM, NIV, IFLAG, IERROR,
     &        TOLEPS, TOL_OPT, KPERCENT, K489, CB_LRB,
     &        WORK, TAU, JPVT, LWORK, RWORK, BLOCK,
     &        MAXI_CLUSTER, KEEP8, OMP_NUM,
     &        NFS4FATHER, NPIV, NVSCHUR_K253, KEEP,  
     &        M_ARRAY,
     &        NELIM, 
     &        NBROWSinF
     &        )
      USE DMUMPS_LR_TYPE, ONLY : LRB_TYPE
      USE DMUMPS_FAC_LR, ONLY : DMUMPS_COMPRESS_CB
      IMPLICIT NONE
      INTEGER(8), intent(in)       :: LA_PTR
      DOUBLE PRECISION, intent(inout)       :: A_PTR(LA_PTR)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER :: sizeBEGS_BLR, sizeBEGS_BLR_U
      INTEGER, intent(in)          :: LDA, NB_ROWS, NB_COLS, NB_INASM
      INTEGER, INTENT(IN)          :: NIV, IWHANDLER, MAXI_CLUSTER, 
     &                                KPERCENT, TOL_OPT, LWORK, OMP_NUM
      INTEGER, INTENT(IN)          :: K489, NROWS, NCOLS, INODE, SYM
      INTEGER, intent(inout)         :: IFLAG, IERROR
      TYPE(LRB_TYPE), intent(inout) :: CB_LRB(NB_ROWS,NB_COLS)
      INTEGER :: BEGS_BLR(sizeBEGS_BLR), BEGS_BLR_U(sizeBEGS_BLR_U)
      DOUBLE PRECISION    :: RWORK(2*MAXI_CLUSTER*OMP_NUM)
      DOUBLE PRECISION :: BLOCK(MAXI_CLUSTER, MAXI_CLUSTER*OMP_NUM)
      DOUBLE PRECISION :: WORK(LWORK*OMP_NUM), TAU(MAXI_CLUSTER*OMP_NUM)
      INTEGER :: JPVT(MAXI_CLUSTER*OMP_NUM)
      INTEGER(8) :: KEEP8(150)
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER, INTENT(in) :: NFS4FATHER, NPIV, NVSCHUR_K253, KEEP(500)
      DOUBLE PRECISION :: M_ARRAY(max(NFS4FATHER,1))
      INTEGER, intent(in)         :: NELIM
      INTEGER, intent(in)         :: NBROWSinF
      CALL DMUMPS_COMPRESS_CB(A_PTR, LA_PTR, POSELT, LDA,
     &        BEGS_BLR, BEGS_BLR_U, 
     &        NB_ROWS, NB_COLS, NB_INASM,
     &        NROWS, NCOLS, INODE,   
     &        IWHANDLER, SYM, NIV, IFLAG, IERROR,
     &        TOLEPS, TOL_OPT, KPERCENT, K489, CB_LRB,
     &        WORK, TAU, JPVT, LWORK, RWORK, BLOCK,
     &        MAXI_CLUSTER, KEEP8, 
     &        NFS4FATHER, NPIV, NVSCHUR_K253, KEEP,  
     &        M_ARRAY=M_ARRAY,
     &        NELIM=NELIM, 
     &        NBROWSinF=NBROWSinF
     &        )
      RETURN
      END SUBROUTINE DMUMPS_COMPRESS_CB_I
      SUBROUTINE DMUMPS_COMPRESS_PANEL_I_NOOPT(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, sizeBEGS_BLR,
     &        NB_BLR, TOLEPS, TOL_OPT, K473, BLR_PANEL, 
     &        CURRENT_BLR,
     &        DIR, WORK, TAU, JPVT, 
     &        LWORK, RWORK, BLOCK,
     &        MAXI_CLUSTER, NELIM, 
     &        LBANDSLAVE, NPIV, ISHIFT, NIV, KPERCENT, 
     &        KEEP8, OMP_NUM
     &        )
      USE DMUMPS_LR_TYPE, ONLY : LRB_TYPE
      USE DMUMPS_FAC_LR, ONLY : DMUMPS_COMPRESS_PANEL
      IMPLICIT NONE
      INTEGER(8), intent(in)       :: LA, POSELT
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR, NIV
      INTEGER, intent(in)          :: OMP_NUM
      INTEGER, intent(inout)          :: IFLAG, IERROR
      TYPE(LRB_TYPE), intent(inout) :: BLR_PANEL(NB_BLR-CURRENT_BLR)
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      INTEGER :: MAXI_CLUSTER
      DOUBLE PRECISION    :: RWORK(2*MAXI_CLUSTER*OMP_NUM)
      DOUBLE PRECISION :: BLOCK(MAXI_CLUSTER,MAXI_CLUSTER*OMP_NUM)
      DOUBLE PRECISION :: WORK(MAXI_CLUSTER*MAXI_CLUSTER*OMP_NUM)
      DOUBLE PRECISION :: TAU(MAXI_CLUSTER*OMP_NUM)
      INTEGER :: JPVT(MAXI_CLUSTER*OMP_NUM)
      INTEGER :: sizeBEGS_BLR
      INTEGER :: BEGS_BLR(sizeBEGS_BLR)
      INTEGER(8) :: KEEP8(150)
      INTEGER, intent(in)          :: NPIV, ISHIFT, KPERCENT, K473,
     &                                TOL_OPT
      LOGICAL, intent(in)          :: LBANDSLAVE
      INTEGER                      :: LWORK, NELIM
      DOUBLE PRECISION,intent(in)              :: TOLEPS
      CHARACTER(len=1) :: DIR
      CALL DMUMPS_COMPRESS_PANEL(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, NB_BLR, TOLEPS, TOL_OPT, K473, BLR_PANEL, 
     &        CURRENT_BLR,
     &        DIR, WORK, TAU, JPVT, 
     &        LWORK, RWORK, BLOCK,
     &        MAXI_CLUSTER, NELIM, 
     &        LBANDSLAVE, NPIV, ISHIFT, NIV, KPERCENT, 
     &        KEEP8
     &        )
      RETURN
      END SUBROUTINE DMUMPS_COMPRESS_PANEL_I_NOOPT
      SUBROUTINE DMUMPS_DECOMPRESS_PANEL_I_NOOPT(
     &        A, LA, POSELT, LDA11,
     &        LDA21, COPY_DENSE_BLOCKS,
     &        BEGS_BLR_DIAG, BEGS_BLR_FIRST_OFFDIAG,
     &        NB_BLR, BLR_PANEL, CURRENT_BLR, DIR, 
     &        DECOMP_TIMER)
      USE DMUMPS_LR_TYPE, ONLY : LRB_TYPE
      USE DMUMPS_FAC_LR, ONLY : DMUMPS_DECOMPRESS_PANEL
      IMPLICIT NONE
      INTEGER(8), intent(in)       :: LA
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      LOGICAL, intent(in)          :: COPY_DENSE_BLOCKS  
      INTEGER, intent(in)          :: NB_BLR, CURRENT_BLR
      INTEGER, intent(in)          :: BEGS_BLR_DIAG, 
     &                                BEGS_BLR_FIRST_OFFDIAG
      TYPE(LRB_TYPE), intent(inout) :: BLR_PANEL(NB_BLR-CURRENT_BLR)
      CHARACTER(len=1) :: DIR
      INTEGER, intent(in) :: DECOMP_TIMER
      INTEGER, intent(in) :: LDA11, LDA21
      CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, LDA11,
     &        LDA21, COPY_DENSE_BLOCKS,
     &        BEGS_BLR_DIAG, BEGS_BLR_FIRST_OFFDIAG,
     &        NB_BLR, BLR_PANEL, CURRENT_BLR, DIR,
     &        DECOMP_TIMER)
      RETURN
      END SUBROUTINE DMUMPS_DECOMPRESS_PANEL_I_NOOPT
      SUBROUTINE DMUMPS_BLR_UPD_NELIM_VAR_L_I(
     &        A_U, LA_U, UPOS, A_L, LA_L, LPOS, IFLAG, IERROR, LDU, LDL,
     &        BEGS_BLR_L, sizeBEGS_BLR_L, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        FIRST_BLOCK, NELIM, UTRANS)
      USE DMUMPS_LR_TYPE, ONLY : LRB_TYPE
      USE DMUMPS_FAC_LR,  ONLY : DMUMPS_BLR_UPD_NELIM_VAR_L
      IMPLICIT NONE
      INTEGER(8), intent(in)       :: LA_U, LA_L
      INTEGER(8), intent(in)       :: UPOS, LPOS
      INTEGER, intent(in)          :: LDU, LDL, NB_BLR_L, CURRENT_BLR,
     &                                NELIM,  FIRST_BLOCK
      CHARACTER(len=1),INTENT(IN)  :: UTRANS
      INTEGER, intent(inout)         :: IFLAG, IERROR
      DOUBLE PRECISION, TARGET, intent(inout) :: A_L(LA_L), A_U(LA_U)
      TYPE(LRB_TYPE),intent(in)      :: BLR_L(NB_BLR_L-CURRENT_BLR)
      INTEGER, INTENT(in)            :: sizeBEGS_BLR_L
      INTEGER                        :: BEGS_BLR_L(sizeBEGS_BLR_L)
      CALL DMUMPS_BLR_UPD_NELIM_VAR_L(
     &        A_U, LA_U, UPOS, A_L, LA_L, LPOS, IFLAG, IERROR, LDU, LDL,
     &        BEGS_BLR_L, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        FIRST_BLOCK, NELIM, UTRANS)
      RETURN
      END SUBROUTINE DMUMPS_BLR_UPD_NELIM_VAR_L_I
      SUBROUTINE DMUMPS_BLR_SLV_UPD_TRAIL_LDLT_I(A, LA, POSELT, 
     &        IFLAG, IERROR, NCOL, NROW,
     &        A_BLOCFACTO, LA_BLOCFACTO, LD_BLOCFACTO, 
     &        BEGS_BLR_LM, sizeBEGS_BLR_LM,
     &        NB_BLR_LM, BLR_LM, ISHIFT_LM,
     &        BEGS_BLR_LS, sizeBEGS_BLR_LS,
     &        NB_BLR_LS, BLR_LS, ISHIFT_LS,
     &        CURRENT_BLR_LM, CURRENT_BLR_LS,
     &        IW2, BLOCK,
     &        MAXI_CLUSTER, OMP_NUM,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT
     &        )
      USE DMUMPS_LR_TYPE, ONLY : LRB_TYPE
      USE DMUMPS_FAC_LR,  ONLY : DMUMPS_BLR_SLV_UPD_TRAIL_LDLT
      IMPLICIT NONE
      INTEGER(8), intent(in)  :: LA, LA_BLOCFACTO
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      DOUBLE PRECISION, intent(in)     :: A_BLOCFACTO(LA_BLOCFACTO)
      INTEGER(8), intent(in)  :: POSELT 
      INTEGER, intent(inout)  :: IFLAG, IERROR
      INTEGER, intent(in)     :: NCOL, NROW, IW2(*), TOL_OPT,
     &                           MAXI_CLUSTER, OMP_NUM, LD_BLOCFACTO
      INTEGER, intent(in)     :: NB_BLR_LM, NB_BLR_LS, 
     &                           ISHIFT_LM, ISHIFT_LS, 
     &                           CURRENT_BLR_LM, CURRENT_BLR_LS
      DOUBLE PRECISION, INTENT(INOUT) ::
     &                      BLOCK(MAXI_CLUSTER,MAXI_CLUSTER*OMP_NUM)
      INTEGER :: sizeBEGS_BLR_LM, sizeBEGS_BLR_LS
      INTEGER :: BEGS_BLR_LM(sizeBEGS_BLR_LM)
      INTEGER :: BEGS_BLR_LS(sizeBEGS_BLR_LS)
      TYPE(LRB_TYPE),intent(in) :: BLR_LM(NB_BLR_LM-CURRENT_BLR_LM),
     &                             BLR_LS(NB_BLR_LS-CURRENT_BLR_LS)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      CALL  DMUMPS_BLR_SLV_UPD_TRAIL_LDLT(A, LA, POSELT, 
     &        IFLAG, IERROR, NCOL, NROW,
     &        A_BLOCFACTO, LA_BLOCFACTO, LD_BLOCFACTO, 
     &        BEGS_BLR_LM, NB_BLR_LM, BLR_LM, ISHIFT_LM,
     &        BEGS_BLR_LS, NB_BLR_LS, BLR_LS, ISHIFT_LS,
     &        CURRENT_BLR_LM, CURRENT_BLR_LS,
     &        IW2, BLOCK,
     &        MAXI_CLUSTER,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT
     &        )
      RETURN
      END SUBROUTINE DMUMPS_BLR_SLV_UPD_TRAIL_LDLT_I
