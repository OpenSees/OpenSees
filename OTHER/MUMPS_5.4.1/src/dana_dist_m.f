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
      SUBROUTINE DMUMPS_ANA_DISTM(MYID, N, STEP, FRERE, FILS, IPOOL,
     & LIPOOL, NE, DAD, ND, PROCNODE, SLAVEF, ABOVE_L0, SIZECB_UNDER_L0,
     & SIZECB_UNDER_L0_IF_LRCB, MAXFR_UNDER_L0,
     & MAX_FRONT_SURFACE_LOCAL_L0, MAX_SIZE_FACTOR_L0, 
     & ENTRIES_IN_FACTORS_UNDER_L0, ENTRIES_IN_FACTORS_MASTERS_LO,
     & COST_SUBTREES_UNDER_LO, OPSA_UNDER_L0, PEAK_FR, PEAK_FR_OOC, 
     & NRLADU, NIRADU, NIRNEC, NRLNEC, NRLNEC_ACTIVE, 
     & NRLADU_if_LR_LU, NRLNEC_if_LR_LU, NRLNEC_if_LR_LUCB, 
     & NRLNECOOC_if_LR_LUCB, NRLNEC_if_LR_CB, NRLADULR_UD, NRLADULR_WC,
     & NRLNECLR_CB_UD, NRLNECLR_LUCB_UD, NRLNECLR_LUCB_WC,
     & PEAK_LRLU_UD,PEAK_OOC_LRLU_UD,PEAK_OOC_LRLU_WC, PEAK_LRLUCB_UD,
     & PEAK_LRLUCB_WC,PEAK_OOC_LRLUCB_UD, PEAK_OOC_LRLUCB_WC,
     & PEAK_LRCB_UD, PEAK_OOC_LRCB_UD, NIRADU_OOC, NIRNEC_OOC, MAXFR,
     & OPSA, UU, KEEP,KEEP8, LOCAL_M, LOCAL_N, SBUF_RECOLD,
     & SBUF_SEND_FR, SBUF_REC_FR, SBUF_SEND_LR, SBUF_REC_LR, 
     & OPS_SUBTREE, NSTEPS, I_AM_CAND,NMB_PAR2, ISTEP_TO_INIV2,
     & CANDIDATES, IFLAG, IERROR, MAX_FRONT_SURFACE_LOCAL,
     & MAX_SIZE_FACTOR, ENTRIES_IN_FACTORS_LOC,
     & ENTRIES_IN_FACTORS_LOC_MASTERS, ROOT_yes, ROOT_NPROW, ROOT_NPCOL
     &     )
      USE DMUMPS_LR_CORE, ONLY : IS_FRONT_BLR_CANDIDATE
      USE DMUMPS_ANA_LR, ONLY : COMPUTE_BLR_VCS
      IMPLICIT NONE
      LOGICAL, intent(in) :: ROOT_yes
      INTEGER, intent(in) :: ROOT_NPROW, ROOT_NPCOL
      INTEGER, intent(in) :: MYID, N, LIPOOL
      LOGICAL, intent(in) :: ABOVE_L0
      INTEGER, intent(in) :: MAXFR_UNDER_L0
      INTEGER(8), intent(in) :: MAX_FRONT_SURFACE_LOCAL_L0, 
     &                          MAX_SIZE_FACTOR_L0, 
     &                          ENTRIES_IN_FACTORS_UNDER_L0,
     &                          ENTRIES_IN_FACTORS_MASTERS_LO
      DOUBLE PRECISION, intent(in) :: COST_SUBTREES_UNDER_LO, 
     &                                OPSA_UNDER_L0
      INTEGER(8), intent(inout) :: SIZECB_UNDER_L0, 
     &                             SIZECB_UNDER_L0_IF_LRCB
      INTEGER, intent(inout) ::  IFLAG, IERROR
      INTEGER  NIRADU, NIRNEC
      INTEGER(8) NRLADU, NRLNEC, NRLNEC_ACTIVE
      INTEGER(8), intent(out) :: NRLADU_if_LR_LU, 
     &  NRLADULR_UD, NRLADULR_WC,
     &  NRLNEC_if_LR_LU, NRLNEC_if_LR_LUCB, NRLNEC_if_LR_CB,
     &  NRLNECOOC_if_LR_LUCB, NRLNECLR_CB_UD,
     &  NRLNECLR_LUCB_UD, NRLNECLR_LUCB_WC
      INTEGER(8), intent(out)::
     &  PEAK_FR, PEAK_FR_OOC,
     &  PEAK_LRLU_UD, 
     &  PEAK_OOC_LRLU_UD, PEAK_OOC_LRLU_WC,
     &  PEAK_LRLUCB_UD, PEAK_LRLUCB_WC,
     &  PEAK_OOC_LRLUCB_UD, PEAK_OOC_LRLUCB_WC,
     &  PEAK_LRCB_UD, PEAK_OOC_LRCB_UD
      INTEGER(8) NRLADU_CURRENT, NRLADU_ROOT_3
      INTEGER NIRADU_OOC, NIRNEC_OOC
      INTEGER MAXFR, NSTEPS
      INTEGER(8) MAX_FRONT_SURFACE_LOCAL
      INTEGER STEP(N)
      INTEGER FRERE(NSTEPS), FILS(N), IPOOL(max(LIPOOL,1)), NE(NSTEPS),
     &        ND(NSTEPS), PROCNODE(NSTEPS), DAD(NSTEPS)
      DOUBLE PRECISION UU
      INTEGER  SLAVEF, KEEP(500), LOCAL_M, LOCAL_N
      INTEGER(8) KEEP8(150)
      INTEGER(8) ENTRIES_IN_FACTORS_LOC,
     &          ENTRIES_IN_FACTORS_LOC_MASTERS
      INTEGER  SBUF_SEND_FR, SBUF_REC_FR
      INTEGER  SBUF_SEND_LR, SBUF_REC_LR
      INTEGER(8) SBUF_RECOLD
      INTEGER  NMB_PAR2
      INTEGER  ISTEP_TO_INIV2( KEEP(71) )
      LOGICAL  I_AM_CAND(NMB_PAR2)
      INTEGER  CANDIDATES( SLAVEF+1, NMB_PAR2 )
      DOUBLE PRECISION OPSA
      DOUBLE PRECISION OPSA_LOC 
      INTEGER(8) MAX_SIZE_FACTOR
      DOUBLE PRECISION OPS_SUBTREE
      DOUBLE PRECISION OPS_SBTR_LOC 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TNSTK, LSTKI 
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: LSTKR
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: LSTKR_if_LRCB,
     &                                         LSTKRLR_CB_UD,
     &                                         LSTKRLR_CB_WC
      LOGICAL  OUTER_SENDS_FR  
      INTEGER(8) :: SAVE_SIZECB_UNDER_L0, 
     &              SAVE_SIZECB_UNDER_L0_IF_LRCB
      INTEGER SBUFR_FR, SBUFS_FR
      INTEGER SBUFR_LR, SBUFS_LR
      INTEGER(8) SBUFS_CB, SBUFR_CB
      INTEGER ITOP,NELIM,NFR
      INTEGER(8) ISTKR, LSTK
      INTEGER(8) :: NRLADU_CURRENT_MISSING
      INTEGER(8) :: ISTKR_if_LRCB, ISTKRLR_CB_UD, ISTKRLR_CB_WC,
     &              K464_8, K465_8
      INTEGER    :: LRSTATUS, IDUMMY
      INTEGER    :: NBNODES_BLR
      LOGICAL    :: COMPRESS_PANEL, COMPRESS_CB
      INTEGER ISTKI,  STKI, ISTKI_OOC
      INTEGER K,NSTK, IFATH
      INTEGER INODE, LEAF, IN
      INTEGER LEVEL, MAXITEMPCB
      INTEGER(8) CURRENT_ACTIVE_MEM, MAXTEMPCB
      INTEGER(8):: MAXTEMPCB_LR
      INTEGER   :: NB_BLR
      LOGICAL UPDATE, UPDATEF, MASTER, MASTERF, INSSARBR
      INTEGER LEVELF, NCB, SIZECBI
      INTEGER(8) NCB8
      INTEGER(8) NFR8, NELIM8
      INTEGER(8) SIZECB, SIZECBINFR, SIZECB_SLAVE
      INTEGER(8) SIZECB_if_LRCB, SIZECB_SLAVE_if_LRCB
      INTEGER(8) SIZECBLR_SLAVE_UD, SIZECBLR_SLAVE_WC
      INTEGER(8) SIZECBLR_UD, SIZECBLR_WC
      INTEGER(8) :: PEAK_DYN_LRLU_UD, PEAK_DYN_LRCB_UD,
     &              PEAK_DYN_LRLUCB_UD, PEAK_DYN_LRLU_WC,
     &              PEAK_DYN_LRLUCB_WC
      INTEGER SIZEHEADER, SIZEHEADER_OOC, XSIZE_OOC
      INTEGER EXTRA_PERM_INFO_OOC
      INTEGER NBROWMAX, NSLAVES_LOC, NSLAVES_PASSED,
     &         NELIMF, NFRF, NCBF,
     &         NBROWMAXF, LKJIB_FR, LKJIB_LR,
     &         NBR, NBCOLFAC
      INTEGER(8) LEV3MAXREC, CBMAXR, CBMAXS
      INTEGER ALLOCOK
      INTEGER PANEL_SIZE
      LOGICAL PACKED_CB
      DOUBLE PRECISION OPS_NODE, OPS_NODE_MASTER, OPS_NODE_SLAVE
      INTEGER(8) ENTRIES_NODE_UPPER_PART, ENTRIES_NODE_LOWER_PART
      INTEGER NBouter_MIN
      INCLUDE 'mumps_headers.h'
      INTEGER WHAT
      INTEGER(8) IDUMMY8
      INTRINSIC min, int, real
      INTEGER DMUMPS_OOC_GET_PANEL_SIZE
      EXTERNAL DMUMPS_OOC_GET_PANEL_SIZE
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      EXTERNAL MUMPS_MAX_SURFCB_NBROWS
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, 
     &         MUMPS_IN_OR_ROOT_SSARBR
      logical :: FORCE_CAND, CONCERNED, UPDATES, STACKCB, MASTERSON
      integer :: IFSON, LEVELSON
      IF (KEEP(50).eq.2) THEN
        EXTRA_PERM_INFO_OOC = 1
      ELSE IF (KEEP(50).eq.0) THEN
        EXTRA_PERM_INFO_OOC = 2
      ELSE
        EXTRA_PERM_INFO_OOC = 0
      ENDIF
      PACKED_CB=( KEEP(215).EQ.0 .AND. KEEP(50).NE.0 )
      MAX_FRONT_SURFACE_LOCAL=0_8
      MAX_SIZE_FACTOR=0_8
      ALLOCATE( LSTKR(NSTEPS), TNSTK(NSTEPS), 
     &          LSTKI(NSTEPS) , 
     &          LSTKR_if_LRCB(NSTEPS), LSTKRLR_CB_UD(NSTEPS),
     &          LSTKRLR_CB_WC(NSTEPS),
     &          stat=ALLOCOK)
      if (ALLOCOK .GT. 0) THEN
        IFLAG  =-7
        IERROR = 6*NSTEPS
        RETURN
      endif
      LKJIB_FR = max(KEEP(5),KEEP(6))
      OUTER_SENDS_FR = (KEEP(263).NE.0 .OR. 
     &      KEEP(50).EQ.0. AND. (KEEP(468).LT.3 .OR. UU.EQ.0.0D0))
      IF ( OUTER_SENDS_FR ) THEN
         LKJIB_FR = max(LKJIB_FR, KEEP(420))
      ENDIF
      LKJIB_LR = max(LKJIB_FR,KEEP(488))
      IF (KEEP(66).NE.0.AND.SLAVEF.GT.1) THEN
         IF ( KEEP(50).EQ.0 ) THEN
           NBouter_MIN = ceiling
     &    (
     &       (dble(KEEP(59))*dble(KEEP(108))*dble(KEEP(35)))
     &       /
     &       (dble(huge(KEEP(108))-10000000)) 
     &    )
         ELSE
           NBouter_MIN = ceiling
     &    (
     &       ( max (dble(KEEP(108))*dble(KEEP(108)), 
     &              dble(KEEP(59))*dble(KEEP(108)/2) 
     &             )
     &        *dble(KEEP(35)))
     &       /
     &       (dble(huge(KEEP(108))-10000000)) 
     &    )
         ENDIF
           NBouter_MIN = max (NBouter_MIN, 4)
           LKJIB_FR = min(KEEP(108)/NBouter_MIN, 4321)
      ENDIF
      TNSTK = NE
      LEAF = LIPOOL+1 
#if defined(OLD_OOC_NOPANEL)
      XSIZE_OOC=XSIZE_OOC_NOPANEL
#else
      IF (KEEP(50).EQ.0) THEN
              XSIZE_OOC=XSIZE_OOC_UNSYM
      ELSE
              XSIZE_OOC=XSIZE_OOC_SYM
      ENDIF
#endif
      SIZEHEADER_OOC = XSIZE_OOC+6 
      SIZEHEADER = XSIZE_IC + 6  
      ISTKR_if_LRCB   = 0_8
      ISTKRLR_CB_UD   = 0_8
      ISTKRLR_CB_WC   = 0_8
      ISTKR           = 0_8
      ISTKI           = 0
      ISTKI_OOC       = 0
      NBNODES_BLR     = 0
      OPSA_LOC   = 0.0D0
      ENTRIES_IN_FACTORS_LOC = 0_8
      ENTRIES_IN_FACTORS_LOC_MASTERS = 0_8
      OPS_SBTR_LOC = 0.0D0
      NRLADU     = 0_8
      NIRADU     = 0
      NIRADU_OOC = 0
      NRLADU_CURRENT = 0_8
      NRLADULR_UD = 0_8
      NRLADULR_WC = 0_8
      NRLADU_ROOT_3 = 0_8
      NRLNEC_ACTIVE = 0_8
      IF (ABOVE_L0) THEN
        SAVE_SIZECB_UNDER_L0 = SIZECB_UNDER_L0
        SAVE_SIZECB_UNDER_L0_IF_LRCB = SIZECB_UNDER_L0_IF_LRCB
      ELSE
        SAVE_SIZECB_UNDER_L0 = 0_8
        SAVE_SIZECB_UNDER_L0_IF_LRCB = 0_8
      ENDIF
      PEAK_DYN_LRLU_UD   = SAVE_SIZECB_UNDER_L0
      PEAK_DYN_LRCB_UD   = SAVE_SIZECB_UNDER_L0_IF_LRCB
      PEAK_DYN_LRLUCB_UD = SAVE_SIZECB_UNDER_L0_IF_LRCB
      PEAK_DYN_LRLU_WC   = SAVE_SIZECB_UNDER_L0
      PEAK_DYN_LRLUCB_WC = SAVE_SIZECB_UNDER_L0
      NRLNEC     = 0_8
      NRLADU_if_LR_LU      = 0_8
      NRLNEC_if_LR_LU      = 0_8
      NRLNEC_if_LR_CB      = 0_8
      NRLNEC_if_LR_LUCB    = 0_8
      NRLNECOOC_if_LR_LUCB = 0_8
      NRLNECLR_CB_UD       = 0_8
      NRLNECLR_LUCB_UD     = 0_8
      NRLNECLR_LUCB_WC     = 0_8
      NIRNEC     = 0
      NIRNEC_OOC = 0
      MAXFR      = 0
      PEAK_FR           = 0_8
      PEAK_FR_OOC       = 0_8
      PEAK_LRLU_UD      = 0_8
      PEAK_OOC_LRLU_UD  = 0_8
      PEAK_OOC_LRLU_WC  = 0_8
      PEAK_LRLUCB_UD    = 0_8
      PEAK_LRLUCB_WC    = 0_8
      PEAK_OOC_LRLUCB_UD= 0_8
      PEAK_OOC_LRLUCB_WC= 0_8
      PEAK_LRCB_UD      = 0_8
      PEAK_OOC_LRCB_UD  = 0_8
      ITOP       = 0
      MAXTEMPCB  = 0_8
      MAXTEMPCB_LR  = 0_8
      MAXITEMPCB = 0
      SBUFS_CB   = 1_8
      SBUFS_FR   = 1
      SBUFS_LR   = 1
      SBUFR_CB   = 1_8
      SBUFR_FR   = 1
      SBUFR_LR   = 1
      IF (KEEP(38) .NE. 0 .AND. KEEP(60).EQ.0) THEN
        INODE  = KEEP(38)
        NRLADU_ROOT_3 = int(LOCAL_M,8)*int(LOCAL_N,8)
        NRLADU = NRLADU_ROOT_3
        NRLNEC_ACTIVE = NRLADU_CURRENT
        MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_ROOT_3)
        NRLNEC = NRLADU
        NRLADU_if_LR_LU      = NRLADU_ROOT_3
        NRLNECOOC_if_LR_LUCB = NRLNEC_ACTIVE
        NRLNEC_if_LR_LU   = NRLADU
        NRLNEC_if_LR_CB   = NRLADU
        NRLNEC_if_LR_LUCB = NRLADU
        PEAK_LRLU_UD      = max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD + SIZECB_UNDER_L0)
        PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
        PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
        PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
        PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
        PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
        PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
        PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
        PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
        IF (MUMPS_PROCNODE(PROCNODE(STEP(INODE)),KEEP(199))
     &                                       .EQ. MYID) THEN
          NIRADU     = SIZEHEADER+2*(ND(STEP(INODE))+KEEP(253))
          NIRADU_OOC = SIZEHEADER_OOC+2*(ND(STEP(INODE))+KEEP(253))
        ELSE
          NIRADU     = SIZEHEADER
          NIRADU_OOC = SIZEHEADER_OOC
        ENDIF
        NIRNEC     = NIRADU
        NIRNEC_OOC = NIRADU_OOC
      ENDIF
      IF((KEEP(24).eq.0).OR.(KEEP(24).eq.1)) THEN
         FORCE_CAND=.FALSE.           
      ELSE
         FORCE_CAND=(mod(KEEP(24),2).eq.0)
      END IF
 90   CONTINUE
      IF (LEAF.NE.1) THEN
         LEAF = LEAF - 1
         INODE = IPOOL(LEAF)
      ELSE
         IF (LIPOOL.NE.0) THEN
           WRITE(MYID+6,*) ' ERROR 1 in DMUMPS_ANA_DISTM '
           CALL MUMPS_ABORT()
         ELSE
           GOTO 115
         ENDIF
      ENDIF
 95   CONTINUE 
      NFR    = ND(STEP(INODE))+KEEP(253)
      NFR8   = int(NFR,8)
      NSTK   = NE(STEP(INODE))
      NELIM = 0 
        IN = INODE
 100    NELIM = NELIM + 1 
      NELIM8=int(NELIM,8)
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
      IFSON = -IN
      IFATH = DAD(STEP(INODE))
      MASTER = MUMPS_PROCNODE(PROCNODE(STEP(INODE)),KEEP(199))
     &           .EQ. MYID
      LEVEL  = MUMPS_TYPENODE(PROCNODE(STEP(INODE)),KEEP(199))
      INSSARBR = MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &        KEEP(199))
      UPDATE=.FALSE.
       if(.NOT.FORCE_CAND) then
         UPDATE = ( (MASTER.AND.(LEVEL.NE.3) ).OR. LEVEL.EQ.2 )
       else
         if(MASTER.and.(LEVEL.ne.3)) then
            UPDATE = .TRUE.
         else if(LEVEL.eq.2) then
            if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(INODE)))) THEN
              UPDATE = .TRUE.
            end if
         end if
       end if
      NCB      = NFR-NELIM
      NCB8     = int(NCB,8)
      SIZECBINFR = NCB8*NCB8
      IF (KEEP(50).EQ.0) THEN
        SIZECB = SIZECBINFR
      ELSE
        IFATH = DAD(STEP(INODE))
        IF ( IFATH.NE.KEEP(38) .AND. PACKED_CB ) THEN
          SIZECB    = (NCB8*(NCB8+1_8))/2_8
        ELSE
          SIZECB    = SIZECBINFR
        ENDIF
      ENDIF
      CALL COMPUTE_BLR_VCS(KEEP(472), NB_BLR, KEEP(488), NELIM)
      IDUMMY = -99999
      CALL IS_FRONT_BLR_CANDIDATE (INODE, LEVEL, NFR, NELIM, 
     &                KEEP(494), 1, KEEP(490),
     &                KEEP(491), KEEP(492),
     &                KEEP(20), KEEP(60), DAD(STEP(INODE)), KEEP(38),
     &                LRSTATUS, IDUMMY)
      COMPRESS_PANEL = (LRSTATUS.GE.2)
      COMPRESS_CB    = ((LRSTATUS.EQ.1).OR.(LRSTATUS.EQ.3))
      IF (COMPRESS_PANEL.OR.COMPRESS_CB) NBNODES_BLR = NBNODES_BLR+1
      IF (COMPRESS_PANEL) THEN
        K464_8 = int(KEEP(464),8)
      ELSE
        K464_8 = 1000_8
      ENDIF
      IF (COMPRESS_CB) THEN
        K465_8 = int(KEEP(465),8)
        SIZECB_if_LRCB = 0_8
        SIZECBLR_UD    = SIZECB*K465_8/1000_8
        SIZECBLR_WC    = SIZECB
      ELSE
        K465_8 = 1000_8
        SIZECBLR_UD    = 0_8
        SIZECBLR_WC    = 0_8
        SIZECB_if_LRCB = SIZECB 
      ENDIF
      SIZECBI      = 2* NCB  + SIZEHEADER
      IF (LEVEL.NE.2) THEN
        NSLAVES_LOC     = -99999999
        SIZECB_SLAVE = -99999997_8
        SIZECB_SLAVE_if_LRCB = SIZECB_SLAVE
        NBROWMAX        = NCB
      ELSE
        IF (KEEP(48) .EQ. 5) THEN
          WHAT = 5 
          IF (FORCE_CAND) THEN
            NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     &                    ISTEP_TO_INIV2(STEP(INODE)))
          ELSE
            NSLAVES_LOC=SLAVEF-1
          ENDIF
          NSLAVES_PASSED=NSLAVES_LOC
        ELSE
          WHAT = 2 
          NSLAVES_PASSED=SLAVEF
          NSLAVES_LOC   =SLAVEF-1
        ENDIF
         CALL MUMPS_MAX_SURFCB_NBROWS(WHAT, KEEP,KEEP8,
     &     NCB, NFR, NSLAVES_PASSED, NBROWMAX, SIZECB_SLAVE
     &    )
       IF (COMPRESS_CB) THEN
        SIZECB_SLAVE_if_LRCB =  0_8
        SIZECBLR_SLAVE_UD = SIZECB_SLAVE*K465_8/1000_8
        SIZECBLR_SLAVE_WC = SIZECB_SLAVE
       ELSE
        SIZECB_SLAVE_if_LRCB = SIZECB_SLAVE
        SIZECBLR_SLAVE_UD    = 0_8
        SIZECBLR_SLAVE_WC    = 0_8
       ENDIF
      ENDIF
      IF (KEEP(60).GT.1) THEN
         IF (MASTER .AND. INODE.EQ.KEEP(38)) THEN
          NIRADU     = NIRADU+SIZEHEADER+2*(ND(STEP(INODE))+KEEP(253))
          NIRADU_OOC = NIRADU_OOC+SIZEHEADER_OOC+
     &                 2*(ND(STEP(INODE))+KEEP(253))
         ENDIF
      ENDIF
      IF (LEVEL.EQ.3) THEN
         IF ( 
     &     KEEP(60).LE.1
     &      ) THEN
           NRLADU_CURRENT = int(LOCAL_M,8)*int(LOCAL_N,8)
           NRLNEC = max(NRLNEC,NRLADU+ISTKR+
     &                 NRLADU_CURRENT)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_ROOT_3 + 
     &                        NRLADU_CURRENT+ISTKR)
           NRLNEC_if_LR_LU   = max(NRLNEC_if_LR_LU,
     &                         NRLADU_if_LR_LU+ISTKR+
     &                         NRLADU_CURRENT)
           NRLNEC_if_LR_CB   = max(NRLNEC_if_LR_CB  ,
     &                         NRLADU+ISTKR_if_LRCB+
     &                         NRLADU_CURRENT)
           NRLNEC_if_LR_LUCB = max(NRLNEC_if_LR_LUCB,
     &                         NRLADU_if_LR_LU+ISTKR_if_LRCB+
     &                         NRLADU_CURRENT)
           NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &                            NRLADU_ROOT_3 + 
     &                            NRLADU_CURRENT+ISTKR_if_LRCB)
           PEAK_LRLU_UD     = max(PEAK_LRLU_UD,
     &                        NRLNEC_if_LR_LU + NRLADULR_UD)
           PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
           PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
           PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
           PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
           PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
           PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
           PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
           PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
         ENDIF
         IF (MASTER) THEN 
            IF (NFR.GT.MAXFR) MAXFR = NFR
         ENDIF
      ENDIF
      IF(KEEP(86).EQ.1)THEN
         IF(MASTER.AND.(.NOT.MUMPS_IN_OR_ROOT_SSARBR(
     &        PROCNODE(STEP(INODE)), KEEP(199)))
     &     )THEN
            IF(LEVEL.EQ.1)THEN
               MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &              NFR8*NFR8)
               IF (KEEP(268).NE.0) THEN
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &              NFR8*NFR8+NELIM8)
               ENDIF
            ELSEIF(LEVEL.EQ.2)THEN
               IF(KEEP(50).EQ.0)THEN
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                 NFR8*NELIM8)
               ELSE
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                 NELIM8*NELIM8)
                 IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
                  MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                  NELIM8*(NELIM8+1_8))
                 ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF (LEVEL.EQ.2) THEN
        IF (MASTER) THEN
          IF (KEEP(50).EQ.0) THEN
             SBUFS_FR = max(SBUFS_FR, NFR*LKJIB_FR+LKJIB_FR+4)
             SBUFS_LR = max(SBUFS_LR, NFR*LKJIB_LR+LKJIB_LR+4)
          ELSE
             SBUFS_FR = max(SBUFS_FR, NELIM*LKJIB_FR+NELIM+6)
             SBUFS_LR = max(SBUFS_LR, NELIM*LKJIB_LR+NELIM+6)
          ENDIF
        ELSEIF (UPDATE) THEN
            if (KEEP(50).EQ.0) THEN
              SBUFR_FR   = max(SBUFR_FR, NFR*LKJIB_FR+LKJIB_FR+4)
              SBUFR_LR   = max(SBUFR_LR, NFR*LKJIB_LR+LKJIB_LR+4)
            else
              SBUFR_FR = max( SBUFR_FR, NELIM*LKJIB_FR+NELIM+6 )
              SBUFR_LR = max( SBUFR_LR, NELIM*LKJIB_LR+NELIM+6 )
              SBUFS_FR = max( SBUFS_FR, NBROWMAX*LKJIB_FR+6 )
              SBUFS_LR = max( SBUFS_LR, NBROWMAX*LKJIB_LR+6 )
              SBUFR_FR = max( SBUFR_FR, NBROWMAX*LKJIB_FR+6 )
              SBUFR_LR = max( SBUFR_LR, NBROWMAX*LKJIB_LR+6 )
            endif
        ENDIF
      ENDIF
      IF ( UPDATE ) THEN
          IF ( (MASTER) .AND. (LEVEL.EQ.1) ) THEN
            NIRADU     = NIRADU + 2*NFR + SIZEHEADER
            NIRADU_OOC = NIRADU_OOC + 2*NFR + SIZEHEADER_OOC
            PANEL_SIZE = DMUMPS_OOC_GET_PANEL_SIZE(
     &      2_8*int(KEEP(226),8), NFR, KEEP(227), KEEP(50))
            NIRADU_OOC = NIRADU_OOC +
     &      EXTRA_PERM_INFO_OOC*(2+NELIM + NELIM/PANEL_SIZE+1)
            IF (KEEP(50).EQ.0) THEN
             NRLADU_CURRENT = int(NELIM,8)*int(2*NFR-NELIM,8)
             NRLADU = NRLADU + NRLADU_CURRENT
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
            ELSE
             NRLADU_CURRENT = int(NELIM,8)*int(NFR,8)
             NRLADU = NRLADU + NRLADU_CURRENT
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
            ENDIF
             IF (COMPRESS_PANEL) THEN
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + 0_8
               NRLADU_CURRENT_MISSING = NRLADU_CURRENT
               NRLADULR_UD     = NRLADULR_UD + 
     &                           NRLADU_CURRENT*K464_8/1000_8
               NRLADULR_WC     = NRLADULR_WC + 
     &                           NRLADU_CURRENT
             ELSE
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + NRLADU_CURRENT
               NRLADU_CURRENT_MISSING = 0_8
             ENDIF
            SIZECBI        = 2* NCB  + SIZEHEADER
          ELSEIF (LEVEL.EQ.2) THEN
            IF (MASTER) THEN
              NIRADU     = NIRADU+SIZEHEADER +SLAVEF-1+2*NFR 
              NIRADU_OOC = NIRADU_OOC+SIZEHEADER_OOC +SLAVEF-1+2*NFR 
              IF (KEEP(50).EQ.0) THEN
                NBCOLFAC=NFR
              ELSE
                NBCOLFAC=NELIM
              ENDIF
              PANEL_SIZE = DMUMPS_OOC_GET_PANEL_SIZE(
     &        2_8*int(KEEP(226),8), NBCOLFAC, KEEP(227), KEEP(50))
              NIRADU_OOC = NIRADU_OOC +
     &        EXTRA_PERM_INFO_OOC*(2+NELIM + NELIM/PANEL_SIZE+1)
              NRLADU_CURRENT = int(NBCOLFAC,8)*int(NELIM,8)
              NRLADU = NRLADU + NRLADU_CURRENT
              IF (COMPRESS_PANEL) THEN
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + 0_8
               NRLADU_CURRENT_MISSING = NRLADU_CURRENT
               NRLADULR_UD     = NRLADULR_UD + 
     &                           NRLADU_CURRENT*K464_8/1000_8
               NRLADULR_WC     = NRLADULR_WC + 
     &                           NRLADU_CURRENT
              ELSE
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + NRLADU_CURRENT
               NRLADU_CURRENT_MISSING = NRLADU_CURRENT
              ENDIF
              MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
               SIZECB     = 0_8
               SIZECBINFR = 0_8
               SIZECB_if_LRCB = 0_8
               SIZECBLR_UD    = 0_8
               SIZECBLR_WC    = 0_8
               SIZECBI    = NCB + SIZEHEADER + SLAVEF - 1
            ELSE
             SIZECB=SIZECB_SLAVE
             SIZECBINFR = SIZECB
             NIRADU       = NIRADU+4+NELIM+NBROWMAX
             NIRADU_OOC   = NIRADU_OOC+4+NELIM+NBROWMAX
             IF (KEEP(50).EQ.0) THEN
               NRLADU_CURRENT   =  int(NELIM,8)*int(NBROWMAX,8)
             ELSE
               NRLADU_CURRENT   = int(NELIM,8)*int(NCB/NSLAVES_LOC,8)
             ENDIF
             NRLADU   = NRLADU + NRLADU_CURRENT
             IF (COMPRESS_PANEL) THEN
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + 0_8
               NRLADU_CURRENT_MISSING = NRLADU_CURRENT
               NRLADULR_UD     = NRLADULR_UD + 
     &                           NRLADU_CURRENT*K464_8/1000_8
               NRLADULR_WC     = NRLADULR_WC + 
     &                           NRLADU_CURRENT
             ELSE
               NRLADU_if_LR_LU = NRLADU_if_LR_LU + NRLADU_CURRENT
               NRLADU_CURRENT_MISSING = 0
             ENDIF
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
             IF (KEEP(50).EQ.0) THEN
               SIZECBI       = 7 + NBROWMAX + NCB
             ELSE
               SIZECBI       = 8 + NBROWMAX + NCB
             ENDIF
             IF (KEEP(50).NE.0) THEN 
                     SIZECBI=SIZECBI+NSLAVES_LOC+
     &                                  XTRA_SLAVES_SYM
             ELSE
                     SIZECBI=SIZECBI+NSLAVES_LOC+
     &                                  XTRA_SLAVES_UNSYM 
             ENDIF
            ENDIF
         ENDIF
         NIRNEC = max0(NIRNEC,
     &             NIRADU+ISTKI+SIZECBI+MAXITEMPCB)
         NIRNEC_OOC = max0(NIRNEC_OOC,
     &             NIRADU_OOC+ISTKI_OOC+SIZECBI+MAXITEMPCB +
     &             (XSIZE_OOC-XSIZE_IC) ) 
         CURRENT_ACTIVE_MEM = ISTKR+SIZECBINFR
         IF (KEEP(50).NE.0.AND.UPDATE.AND.LEVEL.EQ.1) THEN
             CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM +
     &              int(NELIM,8)*int(NCB,8)
         ENDIF
         IF (MASTER .AND.  KEEP(219).NE.0.AND.
     &       KEEP(50).EQ.2.AND.LEVEL.EQ.2) THEN
             CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM + int(NELIM,8)
         ENDIF
         IF (SLAVEF.EQ.1) THEN
           NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
           NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             NRLADU_CURRENT_MISSING)
         ELSE
           NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             MAXTEMPCB_LR)
           NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             MAXTEMPCB_LR+
     &             NRLADU_CURRENT_MISSING)
         ENDIF
         IF (NSTK .NE. 0 .AND. INSSARBR .AND.
     &     KEEP(234).NE.0 .AND. KEEP(55).EQ.0) THEN
           CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTKR(ITOP)
         ENDIF
         IF (SLAVEF.EQ.1) THEN
           NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
           NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+NRLADU_CURRENT_MISSING)
           NRLNECOOC_if_LR_LUCB = 
     &             max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
         ELSE
           NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+MAXTEMPCB)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+MAXTEMPCB)
           NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+MAXTEMPCB+
     &             NRLADU_CURRENT_MISSING)
           NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM-ISTKR+
     &             ISTKR_if_LRCB+MAXTEMPCB_LR)
         ENDIF
         NRLNECLR_LUCB_UD =  max(NRLNECLR_LUCB_UD,
     &                           NRLADULR_UD+ISTKRLR_CB_UD)
         NRLNECLR_LUCB_WC =  max(NRLNECLR_LUCB_WC,
     &                           NRLADULR_WC+ISTKRLR_CB_WC)
         PEAK_LRLU_UD = 
     &             max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD)
         PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
         PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
         PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
         PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
         PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
         PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
         PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
         PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
         IF (NFR.GT.MAXFR) MAXFR = NFR
         IF (NSTK.GT.0) THEN
            DO 70 K=1,NSTK
               LSTK = LSTKR(ITOP)
               ISTKR = ISTKR - LSTK
               IF (K==1 .AND. INSSARBR.AND.KEEP(234).NE.0
     &            .AND.KEEP(55).EQ.0) THEN
               ELSE
                 CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTK
               ENDIF
               LSTK          = LSTKR_if_LRCB(ITOP)
               ISTKR_if_LRCB = ISTKR_if_LRCB - LSTK
               LSTK          = LSTKRLR_CB_UD(ITOP) 
               ISTKRLR_CB_UD = ISTKRLR_CB_UD - LSTK
               LSTK          = LSTKRLR_CB_WC(ITOP) 
               ISTKRLR_CB_WC = ISTKRLR_CB_WC - LSTK
               STKI = LSTKI( ITOP )
               ISTKI = ISTKI - STKI
               ISTKI_OOC = ISTKI_OOC - STKI - (XSIZE_OOC-XSIZE_IC)
               ITOP = ITOP - 1
               IF (ITOP.LT.0) THEN
                  write(*,*) MYID,
     &            ': ERROR 2 in DMUMPS_ANA_DISTM. ITOP = ',ITOP
                  CALL MUMPS_ABORT()
               ENDIF
 70         CONTINUE
         ENDIF 
      ELSE IF (LEVEL.NE.3) THEN
         DO WHILE (IFSON.GT.0) 
            UPDATES=.FALSE.
            MASTERSON = MUMPS_PROCNODE(PROCNODE(STEP(IFSON)),KEEP(199))
     &                  .EQ.MYID
            LEVELSON  = MUMPS_TYPENODE(PROCNODE(STEP(IFSON)),KEEP(199))
            if(.NOT.FORCE_CAND) then
               UPDATES =((MASTERSON.AND.(LEVELSON.NE.3)).OR. 
     &                   LEVELSON.EQ.2)
            else
               if(MASTERSON.and.(LEVELSON.ne.3)) then
                  UPDATES = .TRUE.
               else if(LEVELSON.eq.2) then
                  if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(IFSON)))) then
                    UPDATES = .TRUE.
                  end if
               end if
            end if
            IF (UPDATES) THEN
              LSTK = LSTKR(ITOP)
              ISTKR = ISTKR - LSTK
              LSTK          = LSTKR_if_LRCB(ITOP)
              ISTKR_if_LRCB = ISTKR_if_LRCB - LSTK
              LSTK          =  LSTKRLR_CB_UD(ITOP) 
              ISTKRLR_CB_UD = ISTKRLR_CB_UD - LSTK
              LSTK          =  LSTKRLR_CB_WC(ITOP) 
              ISTKRLR_CB_WC = ISTKRLR_CB_WC - LSTK
              STKI = LSTKI( ITOP )
              ISTKI = ISTKI - STKI
              ISTKI_OOC = ISTKI_OOC - STKI - (XSIZE_OOC-XSIZE_IC)
              ITOP = ITOP - 1
              IF (ITOP.LT.0) THEN
                write(*,*) MYID,
     &          ': ERROR 2 in DMUMPS_ANA_DISTM. ITOP = ',ITOP
                CALL MUMPS_ABORT()
              ENDIF
            ENDIF
            IFSON = FRERE(STEP(IFSON)) 
         END DO
      ENDIF
      IF (
     &        ( (INODE.NE.KEEP(20)).OR.(KEEP(60).EQ.0) ) 
     &       .AND.
     &        ( (INODE.NE.KEEP(38)).OR.(KEEP(60).LE.1) ) 
     &      )
     &  THEN
            ENTRIES_NODE_LOWER_PART = int(NFR-NELIM,8) * int(NELIM,8)
            IF ( KEEP(50).EQ.0 ) THEN
              ENTRIES_NODE_UPPER_PART = int(NFR,8) * int(NELIM,8)
            ELSE
              ENTRIES_NODE_UPPER_PART =
     &        (int(NELIM,8)*int(NELIM+1,8))/2_8
            ENDIF
            IF (KEEP(50).EQ.2 .AND. LEVEL.EQ.3) THEN
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM, 0,
     &           1,OPS_NODE)
            ELSE
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           1,OPS_NODE)
            ENDIF
            IF (LEVEL.EQ.2) THEN
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           2,OPS_NODE_MASTER)
              OPS_NODE_SLAVE=OPS_NODE-OPS_NODE_MASTER
            ENDIF
      ELSE
           OPS_NODE = 0.0D0
           ENTRIES_NODE_UPPER_PART = 0_8
           ENTRIES_NODE_LOWER_PART = 0_8
      ENDIF
      IF ( MASTER )  THEN
        ENTRIES_IN_FACTORS_LOC_MASTERS = 
     &                     ENTRIES_IN_FACTORS_LOC_MASTERS +
     &                            ENTRIES_NODE_UPPER_PART +
     &                            ENTRIES_NODE_LOWER_PART
      ENDIF
      IF (UPDATE.OR.LEVEL.EQ.3) THEN
         IF ( LEVEL .EQ. 3 ) THEN
            IF (ROOT_yes) THEN
              CALL MUMPS_UPDATE_FLOPS_ROOT( OPSA_LOC, KEEP(50), NFR,
     &             NFR, ROOT_NPROW, ROOT_NPCOL, MYID )
              ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                            ENTRIES_NODE_UPPER_PART /
     &                            int(ROOT_NPROW*ROOT_NPCOL,8)
              IF (MASTER) THEN
                ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                                 mod(ENTRIES_NODE_UPPER_PART,
     &                                 int(SLAVEF,8))
              ENDIF
            ENDIF
         ELSE IF (MASTER .AND. LEVEL.EQ.2) THEN
            OPSA_LOC = OPSA_LOC + OPS_NODE_MASTER
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                      ENTRIES_NODE_UPPER_PART +
     &                      mod(ENTRIES_NODE_LOWER_PART,
     &                          int(NSLAVES_LOC,8))
         ELSE IF (MASTER .AND. LEVEL.EQ.1) THEN
            OPSA_LOC = OPSA_LOC + dble(OPS_NODE)
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                               ENTRIES_NODE_UPPER_PART +
     &                               ENTRIES_NODE_LOWER_PART
         ELSE IF (UPDATE) THEN 
            OPSA_LOC = OPSA_LOC + 
     &            dble(OPS_NODE_SLAVE)/dble(NSLAVES_LOC)
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC
     &                 + ENTRIES_NODE_LOWER_PART / 
     &                 int(NSLAVES_LOC,8)
         ENDIF
         IF (MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &        KEEP(199)) .OR. NE(STEP(INODE))==0) THEN
           IF (LEVEL == 1) THEN
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ELSE
             CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           1,OPS_NODE)
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ENDIF
         ENDIF
        ENDIF
      IF (IFATH .EQ. 0) THEN
        IF (LEAF.GT.1) THEN
         GOTO 90 
        ELSE
         GOTO 115 
        ENDIF 
      ELSE
         NFRF = ND(STEP(IFATH))+KEEP(253)
         IF (DAD(STEP(IFATH)).EQ.0) THEN
           NELIMF = NFRF-KEEP(253)
           IF (ABOVE_L0) IN=0
         ELSE
           NELIMF = 0
           IN = IFATH
           DO WHILE (IN.GT.0)
              IN = FILS(IN)
              NELIMF = NELIMF+1
           ENDDO
         ENDIF
         NCBF = NFRF - NELIMF
         LEVELF = MUMPS_TYPENODE(PROCNODE(STEP(IFATH)),KEEP(199))
         MASTERF= MUMPS_PROCNODE(PROCNODE(STEP(IFATH)),
     &                           KEEP(199)).EQ.MYID
         UPDATEF= .FALSE.
         if(.NOT.FORCE_CAND) then
            UPDATEF= ((MASTERF.AND.(LEVELF.NE.3)).OR.LEVELF.EQ.2)
         else
            if(MASTERF.and.(LEVELF.ne.3)) then
               UPDATEF = .TRUE.
            else if (LEVELF.eq.2) then
               if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(IFATH)))) THEN
                 UPDATEF = .TRUE.
               end if
            end if
         end if
         CONCERNED  = UPDATEF .OR. UPDATE
         IF (LEVELF .NE. 2) THEN
           NBROWMAXF = -999999
         ELSE
           IF (KEEP(48) .EQ. 5) THEN
               WHAT = 4
               IF (FORCE_CAND) THEN
                 NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     &               ISTEP_TO_INIV2(STEP(IFATH)))
               ELSE
                 NSLAVES_LOC=SLAVEF-1
               ENDIF
           ELSE
               WHAT = 1 
               NSLAVES_LOC=SLAVEF
           ENDIF
           CALL MUMPS_MAX_SURFCB_NBROWS( WHAT, KEEP, KEEP8,
     &     NCBF, NFRF, NSLAVES_LOC, NBROWMAXF, IDUMMY8 
     &          )
         ENDIF
         IF(LEVEL.EQ.1.AND.UPDATE.AND.
     &      (UPDATEF.OR.LEVELF.EQ.2)
     &      .AND.LEVELF.NE.3) THEN
             IF ( INSSARBR .AND. KEEP(234).NE.0) THEN
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
               NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+
     &             NRLADU_CURRENT_MISSING)
               NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
               NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB+
     &             NRLADU_CURRENT_MISSING)
               NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ISTKR_if_LRCB)
             ELSE
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+SIZECB)
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+SIZECB)
               NRLNEC_if_LR_LU = 
     &             max(NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+SIZECB+
     &             NRLADU_CURRENT_MISSING)
               NRLNEC_if_LR_CB   = 
     &             max(NRLNEC_if_LR_CB  ,NRLADU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ 
     &             ISTKR_if_LRCB+ SIZECB)
               NRLNEC_if_LR_LUCB = 
     &             max(NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR+ 
     &             ISTKR_if_LRCB+ SIZECB+
     &             NRLADU_CURRENT_MISSING)
               NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             NRLADU_CURRENT+NRLADU_ROOT_3+
     &             CURRENT_ACTIVE_MEM-ISTKR+ 
     &             ISTKR_if_LRCB+ SIZECB)
             ENDIF
             PEAK_LRLU_UD = 
     &             max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD)
             PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
             PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
             PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
             PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
             PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
             PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
             PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
             PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
         ENDIF
         IF (UPDATE .AND. LEVEL.EQ.2 .AND. .NOT. MASTER) THEN
             NRLNEC =
     &         max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+NRLADU_CURRENT)
             NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,2_8*NRLADU_CURRENT+
     &         NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
             IF (.NOT.COMPRESS_PANEL) THEN
              NRLNEC_if_LR_LU = max(
     &             NRLNEC_if_LR_LU,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM+NRLADU_CURRENT)
              NRLNEC_if_LR_CB   = max(
     &             NRLNEC_if_LR_CB  ,NRLADU +
     &             CURRENT_ACTIVE_MEM-ISTKR +
     &             ISTKR_if_LRCB+NRLADU_CURRENT)
              NRLNEC_if_LR_LUCB = max(
     &             NRLNEC_if_LR_LUCB,NRLADU_if_LR_LU+
     &             CURRENT_ACTIVE_MEM-ISTKR +
     &             ISTKR_if_LRCB+NRLADU_CURRENT)
              NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &             2_8*NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM-ISTKR+
     &             ISTKR_if_LRCB)
              PEAK_LRLU_UD = 
     &             max(PEAK_LRLU_UD,
     &             NRLNEC_if_LR_LU + NRLADULR_UD)
              PEAK_OOC_LRLU_UD = 
     &             max(PEAK_OOC_LRLU_UD,
     &             NRLNEC_ACTIVE + NRLADULR_UD)
              PEAK_OOC_LRLU_WC = 
     &             max(PEAK_OOC_LRLU_WC,
     &             NRLNEC_ACTIVE + NRLADULR_WC)
              PEAK_LRLUCB_UD = 
     &             max(PEAK_LRLUCB_UD,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
              PEAK_LRLUCB_WC = 
     &             max(PEAK_LRLUCB_WC,
     &             NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
              PEAK_OOC_LRLUCB_UD = 
     &             max(PEAK_OOC_LRLUCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
              PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
              PEAK_LRCB_UD =
     &             max(PEAK_LRCB_UD,
     &             NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
              PEAK_OOC_LRCB_UD =
     &             max(PEAK_OOC_LRCB_UD,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_CB_UD)
             ENDIF
         ENDIF
        IF (LEVELF.EQ.3) THEN
          IF (LEVEL.EQ.1) THEN
            LEV3MAXREC = int(min(NCB,LOCAL_M),8) *
     &                   int(min(NCB,LOCAL_N),8)
          ELSE
            LEV3MAXREC = min(SIZECB,
     &                 int(min(NBROWMAX,LOCAL_M),8)
     &                *int(min(NCB,LOCAL_N),8))
          ENDIF
          MAXTEMPCB  = max(MAXTEMPCB, LEV3MAXREC)
          MAXTEMPCB_LR = max(MAXTEMPCB_LR,LEV3MAXREC)
          MAXITEMPCB = max(MAXITEMPCB,SIZECBI+SIZEHEADER)
          SBUFR_CB   = max(SBUFR_CB, LEV3MAXREC+int(SIZECBI,8))
          NIRNEC = max(NIRNEC,NIRADU+ISTKI+
     &    min(NCB,LOCAL_M)+ min(NCB,LOCAL_N)+SIZEHEADER)
          NIRNEC_OOC = max(NIRNEC_OOC,NIRADU_OOC+ISTKI_OOC+
     &    min(NCB,LOCAL_M)+ min(NCB,LOCAL_N)+SIZEHEADER)
        ENDIF
        IF (CONCERNED) THEN
         IF (LEVELF.EQ.2) THEN
           IF (UPDATE.AND.(LEVEL.NE.2.OR..NOT.MASTER)) THEN
             IF(MASTERF)THEN
                 NBR = min(NBROWMAXF,NBROWMAX)
             ELSE
                 NBR = min(max(NELIMF,NBROWMAXF),NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
               CBMAXS = int(NBR,8)*int(NCB,8)
             ELSE
               CBMAXS = int(NBR,8)*int(NCB,8) -
     &                  (int(NBR,8)*int(NBR-1,8))/2_8
             ENDIF
           ELSE
              CBMAXS = 0_8
           END IF
           IF (MASTERF) THEN
             IF (LEVEL.EQ.1) THEN
                IF (.NOT.UPDATE) THEN
                  NBR = min(NELIMF, NCB)
                ELSE
                  NBR = 0
                ENDIF
             ELSE
                NBR = min(NELIMF, NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
                CBMAXR = int(NBR,8)*NCB8
             ELSE
                CBMAXR = int(NBR,8)*int(min(NCB,NELIMF),8)-
     &                   (int(NBR,8)*int(NBR-1,8))/2_8
                CBMAXR = min(CBMAXR, int(NELIMF,8)*int(NELIMF+1,8)/2_8)
                CBMAXR = min(CBMAXR, SIZECB)
                IF ((LEVEL.EQ.1).AND.(.NOT. PACKED_CB)) THEN
                  CBMAXR = min(CBMAXR,(NCB8*(NCB8+1_8))/2_8)
                ENDIF
             ENDIF
           ELSE IF (UPDATEF) THEN
              NBR = min(NBROWMAXF,NBROWMAX)
              CBMAXR = int(NBR,8) * NCB8
              IF (KEEP(50).NE.0) THEN
                CBMAXR = CBMAXR - (int(NBR,8)*(int(NBR-1,8)))/2_8
              ENDIF
           ELSE
              CBMAXR = 0_8
           ENDIF
         ELSEIF (LEVELF.EQ.3) THEN
           CBMAXR = LEV3MAXREC
           IF (UPDATE.AND. .NOT. (MASTER.AND.LEVEL.EQ.2)) THEN
             CBMAXS = LEV3MAXREC
           ELSE
             CBMAXS = 0_8
           ENDIF
         ELSE
           IF (MASTERF) THEN
             CBMAXS = 0_8
             NBR = min(NFRF,NBROWMAX)
             IF ((LEVEL.EQ.1).AND.UPDATE) THEN
                NBR = 0
             ENDIF
             CBMAXR = int(NBR,8)*int(min(NFRF,NCB),8)
             IF (LEVEL.EQ.2)
     &       CBMAXR = min(CBMAXR, SIZECB_SLAVE)
             IF ( KEEP(50).NE.0 )  THEN
              CBMAXR = min(CBMAXR,(int(NFRF,8)*int(NFRF+1,8))/2_8)
             ELSE
              CBMAXR = min(CBMAXR,int(NFRF,8)*int(NFRF,8))
             ENDIF
           ELSE
             CBMAXR = 0_8
             CBMAXS = SIZECB
           ENDIF
         ENDIF
         IF (UPDATE) THEN
           CBMAXS = min(CBMAXS, SIZECB)
           IF ( .not. ( LEVELF .eq. 1 .AND. UPDATEF ) )THEN
              SBUFS_CB = max(SBUFS_CB, CBMAXS+int(SIZECBI,8))
           ENDIF
         ENDIF
         STACKCB = .FALSE.
         IF (UPDATEF) THEN
          STACKCB = .TRUE.
          SIZECBI     = 2 * NCB + SIZEHEADER
          IF (LEVEL.EQ.1) THEN
             IF (KEEP(50).NE.0.AND.LEVELF.NE.3
     &           .AND.PACKED_CB) THEN
                 SIZECB = (NCB8*(NCB8+1_8))/2_8
             ELSE
                 SIZECB = NCB8*NCB8
             ENDIF
             IF (MASTER) THEN
               IF (MASTERF) THEN
                  SIZECBI     = 2+ XSIZE_IC
               ENDIF
             ELSE IF (LEVELF.EQ.1) THEN
               SIZECB  = min(CBMAXR,SIZECB)
               IF (COMPRESS_CB) THEN
                SIZECBLR_UD = min(SIZECBLR_UD,SIZECB)
                SIZECBLR_WC = min(SIZECBLR_WC,SIZECB)
                SIZECB_if_LRCB = min(SIZECB_if_LRCB,SIZECB)
               ENDIF
               SIZECBI    = 2 * NCB +  9 
               SBUFR_CB   = max(SBUFR_CB, int(SIZECBI,8)+SIZECB)
               SIZECBI    =  2 * NCB + SIZEHEADER     
             ELSE 
               SIZECBI    = 2 * NCB +  9 
               SBUFR_CB   = max(SBUFR_CB, 
     &                      min(SIZECB,CBMAXR) + int(SIZECBI,8))
               MAXTEMPCB  = max(MAXTEMPCB, min(SIZECB,CBMAXR)) 
               IF (COMPRESS_CB) THEN
                MAXTEMPCB_LR = 
     &           max(MAXTEMPCB_LR, (NCB8*int(NB_BLR,8)))
               ELSE
                MAXTEMPCB_LR  = max(MAXTEMPCB_LR, min(SIZECB,CBMAXR))
               ENDIF
               SIZECBI    =  2 * NCB + SIZEHEADER 
               MAXITEMPCB = max(MAXITEMPCB, SIZECBI)
               IF ( .NOT. MASTERF ) THEN
                 SIZECBI     = 0
               ELSE
                  SIZECBI = NCB + SIZEHEADER + SLAVEF - 1
               ENDIF
               SIZECB      = 0_8
               SIZECBLR_UD    = 0_8
               SIZECBLR_WC    = 0_8
               SIZECB_if_LRCB = 0_8
             ENDIF
          ELSE 
             SIZECB = SIZECB_SLAVE
             SIZECBLR_UD    = SIZECBLR_SLAVE_UD
             SIZECBLR_WC    = SIZECBLR_SLAVE_WC
             SIZECB_if_LRCB = SIZECB_SLAVE_if_LRCB
             IF (COMPRESS_CB) THEN
               MAXTEMPCB_LR = 
     &           max(MAXTEMPCB_LR, (NCB8*int(NB_BLR,8)))
             ELSE
                MAXTEMPCB_LR  = max(MAXTEMPCB_LR, min(SIZECB,CBMAXR))
             ENDIF
             MAXTEMPCB  = max(MAXTEMPCB, min(CBMAXR,SIZECB) )
             MAXITEMPCB = max(MAXITEMPCB,NBROWMAX+NCB+SIZEHEADER)
             IF (.NOT. 
     &        (UPDATE.AND.(.NOT.MASTER).AND.(NSLAVES_LOC.EQ.1))
     &          ) 
     &       SBUFR_CB = max(SBUFR_CB, 
     &            min(CBMAXR,SIZECB) + int(NBROWMAX + NCB + 6,8))
             IF (MASTER) THEN
              SIZECB  = 0_8
              SIZECB_SLAVE_if_LRCB =  0_8
              SIZECBLR_UD          = 0_8
              SIZECBLR_WC          = 0_8
              IF (MASTERF) THEN
                SIZECBI = 2 + XSIZE_IC
              ELSE
                SIZECBI = 0
              ENDIF
             ELSE IF (UPDATE) THEN
              IF (MASTERF) THEN
                SIZECBI      =  NCB + SIZEHEADER + SLAVEF - 1
              ELSE
                SIZECBI = 0
              ENDIF
               IF (KEEP(50).EQ.0) THEN
                 SIZECBI = SIZECBI + NBROWMAX + NFR + 
     &                     SIZEHEADER
               ELSE
                 SIZECBI = SIZECBI + NBROWMAX + NFR +
     &                     SIZEHEADER+ NSLAVES_LOC
               ENDIF
             ELSE
              SIZECB      = 0_8
              IF ( MASTERF ) THEN
                SIZECBI = NCB + SIZEHEADER + SLAVEF - 1
              ELSE
                SIZECBI = 0
              ENDIF
              SIZECB_SLAVE_if_LRCB =  0_8
              SIZECBLR_UD          = 0_8
              SIZECBLR_WC          = 0_8
             ENDIF
          ENDIF
         ELSE
           IF (LEVELF.NE.3) THEN
               STACKCB     = .TRUE.
               SIZECB      = 0_8
               SIZECB_SLAVE_if_LRCB =  0_8
               SIZECBLR_UD          = 0_8
               SIZECBLR_WC          = 0_8
               SIZECBI     = 0
               IF ( (LEVEL.EQ.1) .AND. (LEVELF.NE.1) ) THEN
                  IF (PACKED_CB) THEN 
                      SIZECB  = (NCB8*(NCB8+1_8))/2_8
                  ELSE
                      SIZECB  = NCB8*NCB8
                  ENDIF
                  SIZECBI     = 2 * NCB + SIZEHEADER
               ELSE IF (LEVEL.EQ.2) THEN
                 IF (MASTER) THEN
                   SIZECBI=0
                 ELSE 
                   SIZECB  = SIZECB_SLAVE
                   SIZECBLR_UD    = SIZECBLR_SLAVE_UD
                   SIZECBLR_WC    = SIZECBLR_SLAVE_WC
                   SIZECB_if_LRCB = SIZECB_SLAVE_if_LRCB
                   SIZECBI = NBROWMAX + NFR + SIZEHEADER
                 ENDIF 
               ENDIF
           ENDIF
         ENDIF
         IF (STACKCB) THEN
            IF (FRERE(STEP(INODE)).EQ.0) THEN
                  write(*,*) ' ERROR 3 in DMUMPS_ANA_DISTM'
                  CALL MUMPS_ABORT()
           ENDIF
           ITOP = ITOP + 1
           IF ( ITOP .GT. NSTEPS ) THEN
             WRITE(*,*) 'ERROR 4 in DMUMPS_ANA_DISTM '
           ENDIF
           LSTKI(ITOP) = SIZECBI
           ISTKI=ISTKI + SIZECBI
           ISTKI_OOC = ISTKI_OOC + SIZECBI + (XSIZE_OOC-XSIZE_IC)
           LSTKR(ITOP) = SIZECB
           ISTKR = ISTKR + LSTKR(ITOP)
           NRLNEC = max(NRLNEC,NRLADU+ISTKR+MAXTEMPCB)
           NIRNEC = max0(NIRNEC,NIRADU+ISTKI+MAXITEMPCB)
           NIRNEC_OOC = max0(NIRNEC_OOC,NIRADU_OOC+ISTKI_OOC+
     &          MAXITEMPCB + 
     &          (XSIZE_OOC-XSIZE_IC) ) 
           NRLNEC_if_LR_LU  =  max(NRLNEC_if_LR_LU,
     &                             NRLADU_if_LR_LU+ISTKR+MAXTEMPCB)
           LSTKR_if_LRCB(ITOP) = SIZECB_if_LRCB
           ISTKR_if_LRCB       = ISTKR_if_LRCB +  LSTKR_if_LRCB(ITOP)
           LSTKRLR_CB_UD(ITOP) = SIZECBLR_UD
           ISTKRLR_CB_UD       = ISTKRLR_CB_UD + LSTKRLR_CB_UD(ITOP)
           LSTKRLR_CB_WC(ITOP) = SIZECBLR_WC
           ISTKRLR_CB_WC       = ISTKRLR_CB_WC + LSTKRLR_CB_WC(ITOP)
           NRLNECLR_CB_UD    =  max(NRLNECLR_CB_UD, ISTKRLR_CB_UD)
           NRLNECLR_LUCB_UD  =  max(NRLNECLR_LUCB_UD, 
     &                              NRLADULR_UD+ISTKRLR_CB_UD)
           NRLNECLR_LUCB_WC  =  max(NRLNECLR_LUCB_WC, 
     &                              NRLADULR_WC+ISTKRLR_CB_WC)
         ENDIF 
        ENDIF 
         TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
         IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
         ELSE
            GOTO 90
         ENDIF
      ENDIF 
 115  CONTINUE
      NRLNEC = max(NRLNEC, NRLADU+int(KEEP(30),8))
      NRLNEC_ACTIVE = max(NRLNEC_ACTIVE, MAX_SIZE_FACTOR+
     &                    int(KEEP(30),8))
      NRLNEC_if_LR_LU = max(NRLNEC_if_LR_LU,
     &         NRLADU_if_LR_LU + int(KEEP(30),8))
      NRLNEC_if_LR_LUCB = max(NRLNEC_if_LR_LUCB,
     &         NRLADU_if_LR_LU + int(KEEP(30),8))
      NRLNEC_if_LR_CB = max(NRLNEC_if_LR_CB,
     &         NRLADU + int(KEEP(30),8))
      NRLNECOOC_if_LR_LUCB = max(NRLNECOOC_if_LR_LUCB,
     &                    MAX_SIZE_FACTOR+ int(KEEP(30),8))
      PEAK_FR      = SAVE_SIZECB_UNDER_L0 + NRLNEC
      PEAK_FR_OOC  = SAVE_SIZECB_UNDER_L0 + NRLNEC_ACTIVE
      PEAK_LRCB_UD =
     &          max(PEAK_LRCB_UD,
     &          NRLNEC_if_LR_CB + NRLNECLR_CB_UD)
      PEAK_LRLU_UD = 
     &          max(PEAK_LRLU_UD,
     &          NRLNEC_if_LR_LU + NRLADULR_UD)
      PEAK_OOC_LRLU_UD = 
     &          max(PEAK_OOC_LRLU_UD,
     &          NRLNEC_ACTIVE + NRLADULR_UD)
      PEAK_OOC_LRLU_WC = 
     &          max(PEAK_OOC_LRLU_WC,
     &          NRLNEC_ACTIVE + NRLADULR_WC)
      PEAK_LRLUCB_UD = 
     &          max(PEAK_LRLUCB_UD,
     &          NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_UD)
      PEAK_LRLUCB_WC = 
     &          max(PEAK_LRLUCB_WC,
     &          NRLNEC_if_LR_LUCB + NRLNECLR_LUCB_WC)
      PEAK_OOC_LRLUCB_UD = 
     &          max(PEAK_OOC_LRLUCB_UD,
     &          NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_UD)
      PEAK_OOC_LRLUCB_WC = 
     &             max(PEAK_OOC_LRLUCB_WC,
     &             NRLNECOOC_if_LR_LUCB + NRLNECLR_LUCB_WC)
      SBUF_RECOLD = max(int(SBUFR_FR,8),SBUFR_CB)
      SBUF_RECOLD = max(SBUF_RECOLD,
     &        MAXTEMPCB+int(MAXITEMPCB,8)) + 10_8
      SBUF_REC_FR = max(SBUFR_FR, int(min(100000_8,SBUFR_CB))) + 17
      SBUF_REC_LR = max(SBUFR_LR, int(min(100000_8,SBUFR_CB))) + 17
      SBUF_REC_FR = SBUF_REC_FR + 2 * KEEP(127) + SLAVEF - 1 + 7
      SBUF_REC_LR = SBUF_REC_LR + 2 * KEEP(127) + SLAVEF - 1 + 7
      SBUF_SEND_FR = max(SBUFS_FR, int(min(100000_8,SBUFR_CB)))+17
      SBUF_SEND_LR = max(SBUFS_LR, int(min(100000_8,SBUFR_CB)))+17
      IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
         SBUF_RECOLD = SBUF_RECOLD+int(KEEP(108)+1,8)
         SBUF_REC_FR = SBUF_REC_FR+KEEP(108)+1
         SBUF_REC_LR = SBUF_REC_LR+KEEP(108)+1
         SBUF_SEND_FR = SBUF_SEND_FR+KEEP(108)+1
         SBUF_SEND_LR = SBUF_SEND_LR+KEEP(108)+1
      ENDIF
      IF (SLAVEF.EQ.1) THEN 
         SBUF_RECOLD = 1_8
         SBUF_REC_FR = 1
         SBUF_REC_LR = 1
         SBUF_SEND_FR= 1
         SBUF_SEND_LR= 1
      ENDIF
      DEALLOCATE( LSTKR, TNSTK, LSTKI )
      IF (ABOVE_L0) THEN 
        KEEP(470) = KEEP(470)+  NBNODES_BLR
      ELSE
        KEEP(470) = NBNODES_BLR
      ENDIF
      IF (.NOT.ABOVE_L0) THEN
       PEAK_FR     = NRLNEC
       PEAK_FR_OOC = NRLNEC_ACTIVE
      ENDIF
      MAXFR = max(MAXFR, MAXFR_UNDER_L0)
      MAX_FRONT_SURFACE_LOCAL = max (MAX_FRONT_SURFACE_LOCAL,
     &                               MAX_FRONT_SURFACE_LOCAL_L0)
      MAX_SIZE_FACTOR         = max (MAX_SIZE_FACTOR, 
     &                               MAX_SIZE_FACTOR_L0)
      ENTRIES_IN_FACTORS_LOC_MASTERS = ENTRIES_IN_FACTORS_LOC_MASTERS +
     &                          ENTRIES_IN_FACTORS_MASTERS_LO
      ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                         ENTRIES_IN_FACTORS_UNDER_L0
      OPS_SBTR_LOC = OPS_SBTR_LOC + COST_SUBTREES_UNDER_LO
      OPSA_LOC     = OPSA_LOC + OPSA_UNDER_L0
      OPS_SUBTREE = dble(OPS_SBTR_LOC)
      OPSA        = dble(OPSA_LOC)
      RETURN
      END SUBROUTINE DMUMPS_ANA_DISTM
