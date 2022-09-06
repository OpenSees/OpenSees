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
      MODULE MUMPS_BUILD_SORT_INDEX_ELT_M
      CONTAINS
      SUBROUTINE MUMPS_ELT_BUILD_SORT(
     &           NUMELT, LIST_ELT, 
     &           MYID, INODE, N, IOLDPS,
     &           HF, NFRONT, NFRONT_EFF, PERM, 
     &           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, IWPOS,
     &           IFSON, STEP, PIMASTER, PTRIST, PTRAIW, NELT, 
     &           IW, LIW, 
     &           INTARR, LINTARR, ITLOC,
     &           FILS, FRERE_STEPS, 
     &           KEEP, SON_LEVEL2, NIV1, IFLAG, 
     &           DAD, PROCNODE_STEPS, SLAVEF, 
     &           FRT_PTR, FRT_ELT, Pos_First_NUMORG,
     &           SONROWS_PER_ROW, LSONROWS_PER_ROW
     & )
      IMPLICIT NONE
      INTEGER NELT, INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     &        NUMSTK, NUMORG, IFSON, MYID, IFLAG,
     &        NUMELT
      INTEGER KEEP(500)
      INTEGER LIST_ELT(*)
      INTEGER(8), INTENT(IN) :: PTRAIW(NELT+1)
      INTEGER STEP(N), PIMASTER(KEEP(28)), PTRIST(KEEP(28)),
     &        ITLOC(N+KEEP(253)), FILS(N), FRERE_STEPS(KEEP(28)),
     &        PERM(N)
      INTEGER, TARGET :: IW(LIW)
      INTEGER, INTENT(IN), TARGET :: IWPOSCB
      INTEGER, INTENT(IN)         :: IWPOS
      INTEGER(8), INTENT(IN) :: LINTARR
      INTEGER :: INTARR(LINTARR)
      LOGICAL, intent(in)    :: NIV1
      LOGICAL, intent(out)   :: SON_LEVEL2
      INTEGER, intent(out)   :: NFRONT_EFF
      INTEGER, intent(in)    :: DAD (KEEP(28))
      INTEGER, intent(in) :: PROCNODE_STEPS(KEEP(28)), SLAVEF
      INTEGER, intent(in) :: FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER, intent(out) :: Pos_First_NUMORG
      INTEGER, intent(in)    :: LSONROWS_PER_ROW
      INTEGER, intent(out)   :: SONROWS_PER_ROW(LSONROWS_PER_ROW)
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     &        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     &        ITRANS, J, JT1, ISON, IELL, LSTK, 
     &        NROWS, HS, IP1, IP2, IBROT,
     &        I, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     &        LAST_J_ASS, JMIN, MIN_PERM
      INTEGER :: K, K1, K2, K3, KK
      INTEGER(8) :: JJ8, J18, J28
      LOGICAL LEVEL1_SON
      INTEGER INBPROCFILS_SON
      INTEGER TYPESPLIT
      INTEGER ELTI, NUMELT_IBROT
      INCLUDE 'mumps_headers.h'
      INTEGER, POINTER :: SON_IWPOSCB
      INTEGER, POINTER, DIMENSION(:) :: SON_IW
      INTEGER, POINTER, DIMENSION(:) :: PTTRI, PTLAST
      INTEGER                        :: LREQ, allocok
      INTEGER, ALLOCATABLE, TARGET   :: TMP_ALLOC_ARRAY(:)
      INTEGER  MUMPS_TYPESPLIT, MUMPS_TYPENODE
      EXTERNAL MUMPS_TYPESPLIT, MUMPS_TYPENODE 
      IW(IOLDPS+XXNBPR) = 0
      Pos_First_NUMORG = 1
      TYPESPLIT  = MUMPS_TYPESPLIT (PROCNODE_STEPS(STEP(INODE)), 
     &              KEEP(199))
      SON_LEVEL2 = .FALSE.
      IOLDP2     = IOLDPS + HF - 1
      ICT11      = IOLDP2 + NFRONT
      NFRONT_EFF = NASS1
      NTOTFS     = 0
      IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) ) THEN
        K2    = PIMASTER(STEP(IFSON))
        LSTK  = IW(K2    +KEEP(IXSZ))
        NELIM = IW(K2 + 1+KEEP(IXSZ))
        NPIVS  = IW(K2 + 3+KEEP(IXSZ))
        IF (NPIVS.LT.0) NPIVS = 0
        NSLSON = IW(K2 + 5+KEEP(IXSZ))
        IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
        LEVEL1_SON    = NSLSON.EQ.0
        NCOLS  = NPIVS + LSTK
        NROWS  = NCOLS
        ITRANS = NROWS
        IF (NIV1) THEN
          write(6,*) MYID, ':',
     &    ' Internal error 2 in MUMPS_ELT_BUILD_SORT ',
     &    ' interior split node of type 1 '
          CALL MUMPS_ABORT()
        ENDIF
        I= MUMPS_TYPENODE(PROCNODE_STEPS(STEP(IFSON)),KEEP(199))
        J= MUMPS_TYPESPLIT(PROCNODE_STEPS(STEP(IFSON)), 
     &              KEEP(199))
        IF (LEVEL1_SON.or.J.LT.4) THEN
           write(6,*) MYID, ':',
     &     ' Internal error 3 in MUMPS_ELT_BUILD_SORT ',
     &     ' son', IFSON, 
     &     ' of interior split node', INODE, ' of type 1 ', 
     &     ' NSLSON =', NSLSON, ' TYPE_SON=', I, 'TYPESPLIT_SON=', J
           CALL MUMPS_ABORT()
        ENDIF
        SON_IW      => IW
        SON_IWPOSCB => IWPOSCB
        IF (K2 .GT. SON_IWPOSCB) THEN
          INBPROCFILS_SON = K2 + XXNBPR
        ELSE
          INBPROCFILS_SON = PTRIST(STEP(IFSON))+XXNBPR
        ENDIF
        IW(IOLDPS+XXNBPR)=NSLSON
        SON_IW(INBPROCFILS_SON) = NSLSON
        SONROWS_PER_ROW(1:NFRONT-NASS1) = 1
        IF ( K2.GT. IWPOSCB ) THEN
          NROWS = IW(K2 + 2+KEEP(IXSZ))
          ITRANS = NPIVS + NROWS
        ENDIF
        HS = NSLSON + 6 + KEEP(IXSZ)
        K1 = K2 + HS + NROWS + NPIVS
        K2 = K1 + LSTK - 1
        K3 = K1 + NELIM - 1
        IF (NELIM.GT.0) THEN
         DO KK=K1,K3
          NTOTFS = NTOTFS + 1
          JT1 = IW(KK)
          IW(ICT11 + NTOTFS) = JT1
          IW(KK) = NTOTFS
          IW(IOLDP2 + NTOTFS) = IW(KK - ITRANS)
         ENDDO
        ENDIF
        DO KK =K3+1, K2
         NTOTFS = NTOTFS + 1
         JT1 = IW(KK)
         ITLOC(JT1) = NTOTFS 
         IW(KK) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        NFRONT_EFF = NTOTFS
        DO IELL=1,NUMELT
          ELTI = LIST_ELT(IELL)
          J18= PTRAIW(ELTI)
          J28= PTRAIW(ELTI+1)-1
          DO JJ8=J18,J28
           J = INTARR(JJ8)
            INTARR(JJ8) = ITLOC(J)
          ENDDO
        ENDDO
        Pos_First_NUMORG = ITLOC(INODE)
        K1 = IOLDPS+HF
        DO KK=K1+NELIM,K1+NFRONT_EFF-1
          ITLOC(IW(KK)) = 0
        ENDDO
        RETURN   
      ENDIF
      LREQ= 2*NUMSTK
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
       ALLOCATE(TMP_ALLOC_ARRAY(LREQ), stat=allocok)
       IF (allocok .GT. 0) THEN
        IFLAG = -13
        GOTO 800
       ENDIF
       PTTRI  => TMP_ALLOC_ARRAY(1:NUMSTK)
       PTLAST => TMP_ALLOC_ARRAY(NUMSTK+1:LREQ)
      ELSE
       PTTRI  => IW(IWPOS:IWPOS+NUMSTK-1)
       PTLAST => IW(IWPOS+NUMSTK:IWPOS+LREQ)
      ENDIF
      IF (.NOT. NIV1) SONROWS_PER_ROW(1:NFRONT-NASS1) = 0
      IN = INODE
      INEW = IOLDPS + HF
      INEW1 = 1
      DO WHILE (IN.GT.0)
       ITLOC(IN)        = INEW1
       IW(INEW)         = IN
       IW(INEW+NFRONT)  = IN
       INEW1     = INEW1 + 1
       INEW      = INEW + 1
       IN = FILS(IN)
      END DO
      NTOTFS = NUMORG
      IF (NUMSTK .NE. 0) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          K2 = PIMASTER(STEP(ISON))
          SON_IW => IW
          SON_IWPOSCB => IWPOSCB
          LSTK   = SON_IW(K2    +KEEP(IXSZ))
          NELIM  = SON_IW(K2 + 1+KEEP(IXSZ))
          NPIVS  = SON_IW(K2 + 3+KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = SON_IW(K2 + 5+KEEP(IXSZ))
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1_SON    = NSLSON.EQ.0
          NCOLS  = NPIVS + LSTK
          NROWS  = NCOLS
          ITRANS = NROWS
          IF (K2 .GT. SON_IWPOSCB) THEN
            INBPROCFILS_SON = K2+XXNBPR
          ELSE
            INBPROCFILS_SON = PTRIST(STEP(ISON))+XXNBPR
          ENDIF
          IF (NIV1) THEN
           SON_IW(INBPROCFILS_SON) = NSLSON
           IW(IOLDPS+XXNBPR) = IW(IOLDPS+XXNBPR) + NSLSON
          ELSE
           IF (LEVEL1_SON) THEN
            SON_IW(INBPROCFILS_SON) = 1
           ELSE
            SON_IW(INBPROCFILS_SON) = NSLSON
           ENDIF
           IW(IOLDPS+XXNBPR) = IW(IOLDPS+XXNBPR) +
     &                         SON_IW(INBPROCFILS_SON)
          ENDIF
          IF (K2.GT.SON_IWPOSCB) THEN
           NROWS = SON_IW(K2 + 2+KEEP(IXSZ))
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 + KEEP(IXSZ)
          K1 = K2 + HS + NROWS + NPIVS
          K2 = K1 + LSTK - 1 - KEEP(253)
          K3 = K1 + NELIM - 1
          IF (NELIM .NE. 0) THEN
            DO KK = K1, K3
              NTOTFS = NTOTFS + 1
              JT1 = SON_IW(KK)
              IW(ICT11 + NTOTFS) = JT1
              ITLOC(JT1) = NTOTFS
              SON_IW(KK) = NTOTFS
              IW(IOLDP2 + NTOTFS) = SON_IW(KK - ITRANS)
            ENDDO
          ENDIF
          PTTRI(IELL)  = K2+1
          PTLAST(IELL) = K2
          K1 = K3 + 1
          IF (NASS1 .NE. NFRONT - KEEP(253)) THEN
            DO KK = K1, K2
              J = SON_IW(KK)
              IF (ITLOC(J) .EQ. 0) THEN 
                PTTRI(IELL) = KK
                EXIT
              ENDIF
            ENDDO
          ELSE
            DO KK = K1, K2
              SON_IW(KK) = ITLOC(SON_IW(KK))
            ENDDO
            DO KK=K2+1, K2+KEEP(253)
              SON_IW(KK)=NFRONT-KEEP(253)+KK-K2
           ENDDO
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IF (NFRONT-KEEP(253).EQ.NASS1) GOTO 500
      MIN_PERM = N + 1
      JMIN     = -1
      DO IELL = 1, NUMSTK 
        SON_IW => IW
        ILOC = PTTRI( IELL )
        IF ( ILOC .LE. PTLAST( IELL ) ) THEN 
         IF ( PERM( SON_IW( ILOC ) ) .LT. MIN_PERM ) THEN
           JMIN     = SON_IW( ILOC )
           MIN_PERM = PERM( JMIN )
         END IF
        END IF
      END DO
      NEWEL = IOLDP2 + NASS1 + NFRONT
      DO WHILE ( MIN_PERM .NE. N + 1 )
          NEWEL  = NEWEL + 1
          NFRONT_EFF = NFRONT_EFF + 1
          IW( NEWEL ) = JMIN
          ITLOC( JMIN ) = NFRONT_EFF
          LAST_J_ASS = JMIN
          MIN_PERM = N + 1
          DO IELL = 1,  NUMSTK
            SON_IW => IW
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
              IF ( SON_IW( PTTRI( IELL ) ) .eq. LAST_J_ASS )
     &        PTTRI( IELL ) = PTTRI( IELL ) + 1
            ENDIF
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
             IF ( PERM(SON_IW( PTTRI( IELL )) ) .LT. MIN_PERM ) THEN
                JMIN        = SON_IW( PTTRI( IELL ) )
                MIN_PERM = PERM( JMIN )
             END IF
            END IF
          END DO
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT - KEEP(253)) THEN 
      DO IELL = 1,NUMELT
        ELTI = LIST_ELT(IELL)
         J18= PTRAIW(ELTI)
         J28= PTRAIW(ELTI+1)-1_8
         DO JJ8=J18,J28
           J     = INTARR( JJ8 )
           IF ( ITLOC( J ) .eq. 0 ) THEN
            NEWEL  = NEWEL + 1
            NFRONT_EFF = NFRONT_EFF + 1
            IW( NEWEL ) = J
            ITLOC( J ) = NFRONT_EFF
           END IF
         ENDDO
      ENDDO
       IF ( (TYPESPLIT.EQ.4).AND.
     &      (NFRONT_EFF.LT.NFRONT-KEEP(253)) ) THEN
         IBROT = INODE
         DO WHILE
     &      (
     &        ( MUMPS_TYPESPLIT 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),KEEP(199))
     &           .EQ.5 
     &        )
     &        .OR.
     &        ( MUMPS_TYPESPLIT 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),KEEP(199))
     &           .EQ.6  
     &        )
     &      )
          IBROT = DAD(STEP(IBROT))
          NUMELT_IBROT = FRT_PTR(IBROT+1) - FRT_PTR(IBROT)
          IF (NUMELT_IBROT.EQ.0) CYCLE 
          DO IELL = FRT_PTR(IBROT), FRT_PTR(IBROT+1)
            ELTI = FRT_ELT(IELL)
            J18= PTRAIW(ELTI)
            J28= PTRAIW(ELTI+1)-1
            DO JJ8 = J18, J28
              J     = INTARR( JJ8 )
              IF ( ITLOC( J ) .eq. 0 ) THEN
                NEWEL  = NEWEL + 1
                NFRONT_EFF = NFRONT_EFF + 1
                IW( NEWEL ) = J
                ITLOC( J ) = NFRONT_EFF
              END IF
            ENDDO
          ENDDO
          IF (NFRONT_EFF.EQ.NFRONT-KEEP(253)) EXIT
        ENDDO
        IF (NFRONT_EFF.NE.NFRONT-KEEP(253) .AND.
     &       .NOT. (KEEP(376).EQ.1 .AND. KEEP(79) .GE.1)) THEN
          write(6,*) MYID, ': INODE', INODE, ' of type 4 ', 
     &             ' not yet fully assembled ', 
     &             ' NFRONT_EFF, NFRONT =',  NFRONT_EFF, NFRONT
          CALL MUMPS_ABORT()
        ENDIF
       ENDIF
      ENDIF
      IF ( NEWEL1_SAVE .eq. NFRONT_EFF ) THEN
         DO KK=NASS1+1, NFRONT_EFF
           IW( IOLDP2+KK ) = IW( ICT11+KK )
         ENDDO
      ELSE
        CALL MUMPS_SORT( N, PERM, 
     &           IW( NEWEL_SAVE + 1 ), NFRONT_EFF - NEWEL1_SAVE )
        CALL MUMPS_SORTED_MERGE( N, NASS1, PERM, ITLOC,
     &    IW( NEWEL_SAVE + 1), NFRONT_EFF - NEWEL1_SAVE,
     &    IW( ICT11  + NASS1 + 1 ), NEWEL1_SAVE - NASS1,
     &    IW( IOLDP2 + NASS1 + 1 ), NFRONT_EFF - NASS1 )
        DO KK = NASS1+1, NFRONT_EFF
          IW(ICT11 + KK) = IW(IOLDP2+KK)
        ENDDO
      END IF
  500 CONTINUE
      IF ( KEEP(253).GT.0) THEN
        IP1 = IOLDPS +  HF + NFRONT_EFF  
        IP2 = IOLDPS + HF + NFRONT + NFRONT_EFF 
        DO I= 1, KEEP(253)
          IW(IP1+I-1) = N+I
          IW(IP2+I-1) = N+I
          ITLOC(N+I)  = NFRONT_EFF + I
        ENDDO
        NFRONT_EFF = NFRONT_EFF + KEEP(253)
      ENDIF
      IF (NFRONT.GT.NFRONT_EFF) THEN
        IP1 = IOLDPS + NFRONT + HF
        IP2 = IOLDPS + NFRONT_EFF + HF
        DO I=1, NFRONT_EFF
          IW(IP2+I-1)=IW(IP1+I-1)
        ENDDO
      ELSE IF (NFRONT .LT. NFRONT_EFF) THEN
        WRITE(*,*) "Internal error in MUMPS_ELT_BUILD_SORT",
     &             NFRONT, NFRONT_EFF
        IFLAG = -53
        GOTO 800
      ENDIF
      IF ( (NUMSTK .NE.0) 
     & .AND. (NFRONT-KEEP(253).GT.NASS1 )   
     &  ) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          K2 = PIMASTER(STEP(ISON))
          SON_IW => IW
          SON_IWPOSCB => IWPOSCB
          LSTK = SON_IW(K2+KEEP(IXSZ))
          NELIM = SON_IW(K2 + 1 +KEEP(IXSZ))
          NPIVS = SON_IW(K2 + 3 +KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = SON_IW(K2 + 5 +KEEP(IXSZ))
          LEVEL1_SON = (NSLSON .EQ. 0)
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          IF (K2.GT.SON_IWPOSCB) THEN
           NROWS = SON_IW(K2 + 2+KEEP(IXSZ))
          ENDIF
          HS = NSLSON + 6 +KEEP(IXSZ)
          K1 = K2 + HS + NROWS + NPIVS
          K2 = K1 + LSTK - 1
          K3 = K1 + NELIM - 1
          K1 = K3 + 1
          IF (NFRONT-KEEP(253).GT.NASS1) THEN
            DO KK = K1, K2
              J = SON_IW(KK)
              SON_IW(KK) = ITLOC(J)
              IF (NIV1 .AND. NSLSON.EQ.0) THEN
              ELSE
                IF (SON_IW(KK) .LE. NASS1 .OR. NIV1) THEN
                ELSE
                  SONROWS_PER_ROW(SON_IW(KK)-NASS1) =
     &                         SONROWS_PER_ROW(SON_IW(KK)-NASS1) + 1
                ENDIF
              ENDIF
            ENDDO
          ELSE
             IF (.not. NIV1) THEN
                WRITE(*,*) "Internal error 1 in MUMPS_ELT_BUILD_SORT"
                CALL MUMPS_ABORT() 
              ENDIF
              IF (.not.LEVEL1_SON) THEN
              ENDIF
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      DO IELL=1,NUMELT
        ELTI = LIST_ELT(IELL)
        J18 = PTRAIW(ELTI)
        J28 = PTRAIW(ELTI+1)-1
        DO JJ8=J18,J28
          J = INTARR(JJ8)
            INTARR(JJ8) = ITLOC(J)
        ENDDO
      ENDDO
        K1 = IOLDPS + HF + NUMORG
        K2 = K1 + NFRONT_EFF - 1 + NASS
        DO K = K1, K2
          I = IW(K)
          ITLOC(I) = 0
        ENDDO
  800 CONTINUE
      IF (allocated(TMP_ALLOC_ARRAY)) DEALLOCATE(TMP_ALLOC_ARRAY)
      RETURN
      END SUBROUTINE MUMPS_ELT_BUILD_SORT
      END MODULE MUMPS_BUILD_SORT_INDEX_ELT_M
