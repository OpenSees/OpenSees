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
      MODULE MUMPS_BUILD_SORT_INDEX_M
      CONTAINS
      SUBROUTINE MUMPS_BUILD_SORT_INDEX(
     &           MYID, INODE, N, IOLDPS, HF, LP, LPOK, 
     &           NFRONT, NFRONT_EFF, PERM, DAD,
     &           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, IWPOS,
     &           IFSON, STEP, PIMASTER, PTRIST, PTRAIW, IW, LIW, 
     &           INTARR, LINTARR, ITLOC, FILS, FRERE_STEPS, 
     &           SON_LEVEL2, NIV1, KEEP,KEEP8, IFLAG,
     &           ISON_IN_PLACE, PROCNODE_STEPS, SLAVEF,
     &           SONROWS_PER_ROW, LSONROWS_PER_ROW
     & )
      IMPLICIT NONE
      INTEGER INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     &        NUMSTK, NUMORG, IFSON, MYID, LP
      LOGICAL, intent(in) :: LPOK
      INTEGER, intent(in) :: ISON_IN_PLACE
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8), INTENT(IN) :: PTRAIW(N)
      INTEGER STEP(N), PIMASTER(KEEP(28)), PTRIST(KEEP(28)),
     &        ITLOC(N+KEEP(253)), FILS(N), FRERE_STEPS(KEEP(28)),
     &        PERM(N)
      INTEGER, TARGET :: IW(LIW)
      INTEGER, INTENT(IN), TARGET :: IWPOSCB
      INTEGER, INTENT(IN)         :: IWPOS
      INTEGER(8), INTENT(IN) :: LINTARR
      INTEGER :: INTARR(LINTARR)
      LOGICAL, intent(in)    :: NIV1
      INTEGER, intent(inout) :: IFLAG
      LOGICAL, intent(out)   :: SON_LEVEL2
      INTEGER, intent(out)   :: NFRONT_EFF
      INTEGER, intent(in)    :: DAD (KEEP(28))
      INTEGER, intent(in) :: PROCNODE_STEPS(KEEP(28)), SLAVEF
      INTEGER, intent(in)    :: LSONROWS_PER_ROW
      INTEGER, intent(out)   :: SONROWS_PER_ROW(LSONROWS_PER_ROW)
      INTEGER NELIM_SON_IN_PLACE 
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     &        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     &        ITRANS, J, JT1, ISON, IELL, LSTK, 
     &        NROWS, HS, IP1, IP2, IBROT, IORG, 
     &        I, K, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     &        LAST_J_ASS, JMIN, MIN_PERM
      LOGICAL LEVEL1_SON
      INTEGER :: K1, K2, K3, KK
      INTEGER(8) :: J18, J28, JJ8, JDEBROW8
      INTEGER INBPROCFILS_SON
      INTEGER TYPESPLIT
      INCLUDE 'mumps_headers.h'
      INTEGER, POINTER :: SON_IWPOSCB
      INTEGER, POINTER, DIMENSION(:) :: SON_IW
      INTEGER, POINTER, DIMENSION(:) :: PTTRI, PTLAST
      INTEGER                        :: LREQ, allocok
      INTEGER, ALLOCATABLE, TARGET   :: TMP_ALLOC_ARRAY(:)
      INTEGER  MUMPS_TYPESPLIT, MUMPS_TYPENODE
      EXTERNAL MUMPS_TYPESPLIT, MUMPS_TYPENODE 
      IW(IOLDPS+XXNBPR) = 0
      TYPESPLIT  = MUMPS_TYPESPLIT (PROCNODE_STEPS(STEP(INODE)), 
     &              KEEP(199))
      SON_LEVEL2 = .FALSE.
      IOLDP2     = IOLDPS + HF - 1
      ICT11      = IOLDP2 + NFRONT
      NTOTFS = 0
      NELIM_SON_IN_PLACE = 0
      IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) ) THEN
        K2    = PIMASTER(STEP(IFSON))
        LSTK  = IW(K2    +KEEP(IXSZ))
        NELIM = IW(K2 + 1+KEEP(IXSZ))
        IF ( ISON_IN_PLACE > 0 ) THEN
           IF (ISON_IN_PLACE.NE.IFSON) THEN
         write(6,*) MYID, ':',
     &   ' Internal error 1 in MUMPS_BUILD_SORT_INDEX ',
     &   ' in place node is not the first son a interior split node '
         CALL MUMPS_ABORT()
          ENDIF
          NELIM_SON_IN_PLACE = NELIM
        ENDIF
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
     &    ' Internal error 2 in MUMPS_BUILD_SORT_INDEX ',
     &    ' interior split node of type 1 '
          CALL MUMPS_ABORT()
        ENDIF
        I= MUMPS_TYPENODE(PROCNODE_STEPS(STEP(IFSON)),KEEP(199))
        J= MUMPS_TYPESPLIT(PROCNODE_STEPS(STEP(IFSON)), 
     &              KEEP(199))
        IF (LEVEL1_SON.or.J.LT.4) THEN
           write(6,*) MYID, ':',
     &     ' Internal error 3 in MUMPS_BUILD_SORT_INDEX ',
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
        DO KK =K3+1, K3+NUMORG  
         NTOTFS = NTOTFS + 1
         JT1 = IW(KK)
         ITLOC(JT1) = NTOTFS   
         IW(KK) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        DO KK =K3+NUMORG+1, K2
         NTOTFS = NTOTFS + 1
         JT1 = IW(KK)
         ITLOC(JT1) = NTOTFS 
         IW(KK) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        NFRONT_EFF = NTOTFS
        IBROT = INODE
        DO IORG = 1, NUMORG
          J18 = PTRAIW(IBROT) + 2
          JT1 = INTARR(J18)
          INTARR(J18) = ITLOC(JT1)
          IBROT = FILS(IBROT)
         J28 = J18 + INTARR(J18 - 2) - INTARR(J18 - 1)
         J18 = J18 + 1
         IF (J18 .LE. J28) THEN
          DO JJ8 = J18, J28
            J = INTARR(JJ8)
            INTARR(JJ8) = ITLOC(J)
          ENDDO
         ENDIF
        ENDDO
        K1 = IOLDPS+HF
        DO KK=K1+NELIM,K1+NFRONT_EFF-1
          ITLOC(IW(KK)) = 0
        ENDDO
        RETURN   
      ENDIF
      LREQ= 2*NUMSTK+2
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
       ALLOCATE(TMP_ALLOC_ARRAY(LREQ), stat=allocok)
       IF (allocok .GT. 0) THEN
        IFLAG = -13
        GOTO 800
       ENDIF
       PTTRI  => TMP_ALLOC_ARRAY(1:NUMSTK+1)
       PTLAST => TMP_ALLOC_ARRAY(NUMSTK+2:LREQ)
      ELSE
       PTTRI  => IW(IWPOS:IWPOS+NUMSTK)
       PTLAST => IW(IWPOS+NUMSTK+1:IWPOS+LREQ-1)
      ENDIF
      NFRONT_EFF = NASS1
      IF ( ISON_IN_PLACE > 0 ) THEN
        ISON   = ISON_IN_PLACE
        K2     = PIMASTER(STEP(ISON))
        LSTK   = IW(K2    +KEEP(IXSZ))
        NELIM  = IW(K2 + 1+KEEP(IXSZ))
        NPIVS  = IW(K2 + 3+KEEP(IXSZ))
        IF (NPIVS.LT.0) NPIVS = 0
        NSLSON = IW(K2 + 5+KEEP(IXSZ))
        NCOLS  = NPIVS + LSTK
        NROWS  = NCOLS
        ITRANS = NROWS
        IF ( K2 .GT. IWPOSCB ) THEN
          NROWS = IW(K2 + 2+KEEP(IXSZ))
          ITRANS = NPIVS + NROWS
        ENDIF
        HS = NSLSON + 6 + KEEP(IXSZ)
        K1 = K2 + HS + NROWS + NPIVS
        K2 = K1 + LSTK - 1
        K3 = K1 + NELIM - 1
        DO KK = K1, K3
          NTOTFS = NTOTFS + 1
          JT1 = IW(KK)
          IW(ICT11 + NTOTFS) = JT1
          ITLOC(JT1) = NTOTFS
          IW(KK) = NTOTFS
          IW(IOLDP2 + NTOTFS) = IW(KK - ITRANS)
        ENDDO
        NELIM_SON_IN_PLACE = NTOTFS
      ENDIF
      IF (.NOT. NIV1) SONROWS_PER_ROW(1:NFRONT-NASS1) = 0
      IN = INODE
      INEW = IOLDPS + HF +  NTOTFS
      INEW1 = NTOTFS + 1
      JDEBROW8 = PTRAIW(INODE)+3
      PTTRI(NUMSTK+1)  = 0
      PTLAST(NUMSTK+1) = 0 + INTARR(JDEBROW8-3) - 1
   50 CONTINUE
      J18 = PTRAIW(IN) + 2
      JT1 = INTARR(J18)
      INTARR(J18) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW)         = JT1
      IW(INEW+NFRONT)  = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      NTOTFS = NTOTFS + NUMORG
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
          IF ( K2 .GT. SON_IWPOSCB ) THEN
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
          IF (NELIM .NE. 0 .AND. ISON.NE.ISON_IN_PLACE) THEN
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
 199  CONTINUE
      IF ( PTTRI( NUMSTK + 1 ) .LE. PTLAST( NUMSTK + 1 ) ) THEN
      IF ( ITLOC( INTARR( JDEBROW8+PTTRI( NUMSTK + 1 ) ) ) .NE. 0 ) THEN
       PTTRI( NUMSTK + 1 ) = PTTRI( NUMSTK + 1 ) + 1
       GOTO 199
      END IF
      END IF
      MIN_PERM = N + 1
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
      IELL = NUMSTK + 1
      ILOC =  PTTRI( IELL )
      IF ( ILOC .LE. PTLAST( IELL ) ) THEN
        IF ( PERM( INTARR( JDEBROW8+ILOC ) ) .LT. MIN_PERM ) THEN
         JMIN        = INTARR( JDEBROW8+ILOC )
         MIN_PERM = PERM( JMIN )
        END IF
      END IF
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
          IELL = NUMSTK + 1
 145      CONTINUE
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
            IF ( INTARR( JDEBROW8+PTTRI( IELL ) ) .eq. LAST_J_ASS ) THEN
              PTTRI( IELL ) = PTTRI( IELL ) + 1 
              GOTO 145
            END IF
          END IF
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
            IF (PERM(INTARR( JDEBROW8+PTTRI(IELL) )) .LT. MIN_PERM) THEN
              JMIN        = INTARR( JDEBROW8+PTTRI(IELL) )
              MIN_PERM = PERM( JMIN )
            END IF
          END IF
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT - KEEP(253)) THEN 
      IBROT = INODE
      DO IORG = 1, NUMORG
         J18    = PTRAIW(IBROT) + 2
         J28    = J18 + INTARR(J18 - 2) - INTARR(J18-1)
         IBROT = FILS( IBROT )
         IF ( IORG.EQ. 1) THEN
           IF ( KEEP(50).NE.0 ) CYCLE
           J18 = J18 + 1 + INTARR(J18-2)
         ELSE
           J18 = J18 + 1
         ENDIF
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
          IN = IBROT
          DO WHILE (IN.GT.0.AND.NFRONT_EFF.LT.NFRONT-KEEP(253))
            J18    = PTRAIW(IN) + 2
            J28    = J18 + INTARR(J18 - 2) - INTARR(J18-1)
            IN = FILS( IN )
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
        IF (LPOK) THEN
         WRITE(LP,*) " Error in MUMPS_BUILD_SORT_INDEX:",
     &     " matrix structure might have changed,", 
     &     " analysis (JOB=1) should be performed again ", 
     &     " NFRONTexpected, NFRONTeffective=", NFRONT, NFRONT_EFF
        ENDIF
       IFLAG = -53
       GOTO 800
      ENDIF
      IF ( NUMSTK .NE. 0  
     &    .AND. (NFRONT-KEEP(253).GT.NASS1)
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
                WRITE(*,*) "Internal error 1 in MUMPS_BUILD_SORT_INDEX"
                CALL MUMPS_ABORT() 
              ENDIF
              IF (.not.LEVEL1_SON) THEN
              ENDIF
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IBROT = INODE
      DO IORG = 1, NUMORG
        J18 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J28 = J18 + INTARR(J18 - 2) - INTARR(J18 - 1)
        J18 = J18 + 1
        DO JJ8 = J18, J28
          J = INTARR(JJ8)
            INTARR(JJ8) = ITLOC(J)
        ENDDO
      ENDDO
        K1 = IOLDPS + HF
        K2 = K1 + NFRONT_EFF -1
        IF (KEEP(50).EQ.0) K2 = K2 + NELIM_SON_IN_PLACE
        DO K = K1, K2
          I = IW(K)
          ITLOC(I) = 0
        ENDDO
        IF (KEEP(50).EQ.0) THEN
          K1 = IOLDPS+HF+NFRONT_EFF+NELIM_SON_IN_PLACE+NUMORG
          K2 = K1 + NASS -NELIM_SON_IN_PLACE - 1
          DO K = K1, K2
            I = IW(K)
            ITLOC(I) = 0
          ENDDO
        ENDIF
  800 CONTINUE
      IF (allocated(TMP_ALLOC_ARRAY)) DEALLOCATE(TMP_ALLOC_ARRAY)
      RETURN
      END SUBROUTINE MUMPS_BUILD_SORT_INDEX
      END MODULE MUMPS_BUILD_SORT_INDEX_M
      SUBROUTINE MUMPS_SORT( N, PERM, IW, LIW )
      IMPLICIT NONE
      INTEGER N, LIW
      INTEGER PERM( N ), IW( LIW )
      INTEGER I, SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, LIW - 1
          IF ( PERM( IW( I ) ) .GT. PERM( IW( I + 1 ) ) ) THEN
            DONE = .FALSE.
            SWAP  = IW( I + 1 )
            IW( I + 1 ) = IW( I )
            IW( I ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT
      SUBROUTINE MUMPS_SORTED_MERGE( N, NASS1, PERM, ITLOC,
     &                             SMALL, LSMALL,
     &                             LARGE, LLARGE,
     &                             MERGE, LMERGE )
      IMPLICIT NONE
      INTEGER N, NASS1, LSMALL, LLARGE, LMERGE
      INTEGER PERM( N ), ITLOC( N ) 
      INTEGER SMALL(LSMALL), LARGE(LLARGE), MERGE(LMERGE)
      INTEGER PSMALL, PLARGE, PMERGE, VSMALL, VLARGE, VMERGE
      PSMALL = 1
      PLARGE = 1
      PMERGE = 1
      DO WHILE ( PSMALL .LE. LSMALL .or. PLARGE.LE. LLARGE )
        IF ( PSMALL .GT. LSMALL ) THEN
          VMERGE = LARGE( PLARGE )
          PLARGE = PLARGE + 1
        ELSE IF ( PLARGE .GT. LLARGE ) THEN
          VMERGE = SMALL( PSMALL )
          PSMALL = PSMALL + 1
        ELSE
          VSMALL = SMALL( PSMALL )
          VLARGE = LARGE( PLARGE )
          IF ( PERM( VSMALL ) .LT. PERM( VLARGE ) ) THEN
            VMERGE = VSMALL
            PSMALL   = PSMALL + 1
          ELSE
            VMERGE = VLARGE
            PLARGE   = PLARGE + 1
          END IF
        END IF
        MERGE( PMERGE ) = VMERGE
        ITLOC( VMERGE ) = PMERGE + NASS1
        PMERGE = PMERGE + 1
      END DO
      PMERGE = PMERGE - 1
      RETURN
      END SUBROUTINE MUMPS_SORTED_MERGE
