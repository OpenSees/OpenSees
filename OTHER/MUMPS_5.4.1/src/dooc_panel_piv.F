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
C     This file contains routines related to OOC,
C     panels, and pivoting. They are used to store
C     permutation information of what is already on
C     disk to be able to permute things back at the
C     solve stage.
C     They do not need to be in the MUMPS_OOC
C     module (most of them do not use any variable
C     from the module, or are called from routines
C     where we do not necessarily want to do a
C     USE DMUMPS_OOC).
      INTEGER FUNCTION DMUMPS_OOC_GET_PANEL_SIZE
     &     ( HBUF_SIZE, NNMAX, K227, K50 )
      IMPLICIT NONE
C
C     Arguments:
C     =========
C
      INTEGER, INTENT(IN) :: NNMAX, K227, K50
      INTEGER(8), INTENT(IN) :: HBUF_SIZE
C
C     Purpose:
C     =======
C
C     - Compute the effective size (maximum number of pivots in a panel)
C     for a front with NNMAX entries in its row (for U) /
C     column (for L).
C     - Be able to adapt the fixed number of columns in panel
C     depending on NNMAX, and size of IO buffer HBUF_SIZE
C
C     Local variables
C     ===============
C
      INTEGER K227_LOC
      INTEGER NBCOL_MAX
      INTEGER EFFECTIVE_SIZE
      NBCOL_MAX=int(HBUF_SIZE / int(NNMAX,8))
C     KEEP(227): Maximum size (nb of col/row) of a panel
      K227_LOC = abs(K227)
      IF (K50.EQ.2) THEN
C        for 2x2 pivots we may end-up having the first part
C        of a 2x2 pivot in the last col of the panel; the
C        adopted solution consists in adding the next column
C        to the panel; therefore we need be able to
C        dynamically increase the panel size by one.
C        note that we also maintain property:
C        KEEP(227): Maximum size (nb of col/row) of a panel
         K227_LOC=max(K227_LOC,2)
         EFFECTIVE_SIZE =  min(NBCOL_MAX-1, K227_LOC-1)
cN       - during bwd the effective size is useless
      ELSE
C        complete buffer space can be used for a panel
         EFFECTIVE_SIZE =  min(NBCOL_MAX, K227_LOC)
      ENDIF
      IF (EFFECTIVE_SIZE.LE.0) THEN
         write(6,*) 'Internal buffers too small to store ',
     &        ' ONE col/row of size', NNMAX
         CALL MUMPS_ABORT()
      ENDIF
      DMUMPS_OOC_GET_PANEL_SIZE = EFFECTIVE_SIZE
      RETURN
      END FUNCTION DMUMPS_OOC_GET_PANEL_SIZE
C
      SUBROUTINE DMUMPS_PERMUTE_PANEL( IPIV, LPIV, ISHIFT,
     &     THE_PANEL, NBROW, NBCOL, KbeforePanel )
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Permute rows of a panel, stored by columns, according
C     to permutation array IPIV.
C     IPIV is such that, for I = 1 to LPIV, row ISHIFT + I
C     in the front must be permuted with row IPIV( I )
C
C     Since the panel is not necessary at the beginning of
C     the front, let KbeforePanel be the number of pivots in the
C     front before the first pivot of the panel.
C
C     In the panel, row ISHIFT+I-KbeforePanel is permuted with
C     row IPIV(I)-KbeforePanel
C
C     Note:
C     ====
C
C     This routine can also be used to permute the columns of
C     a matrix (U) stored by rows. In that case, the argument
C     NBROW represents the number of columns, and NBCOL represents
C     the number of rows.
C
C
C     Arguments:
C     =========
C
      INTEGER LPIV, ISHIFT, NBROW, NBCOL, KbeforePanel
      INTEGER IPIV(LPIV)
      DOUBLE PRECISION THE_PANEL(NBROW, NBCOL)
C
C     Local variables:
C     ===============
C
      INTEGER I, IPERM
C
C     Executable statements
C     =====================
C
      DO I = 1, LPIV
C        Swap rows ISHIFT + I and PIV(I)
         IPERM=IPIV(I)
         IF ( I+ISHIFT.NE.IPERM) THEN
            CALL dswap(NBCOL,
     &           THE_PANEL(I+ISHIFT-KbeforePanel,1), NBROW,
     &           THE_PANEL(IPERM-KbeforePanel,1), NBROW)
         ENDIF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_PERMUTE_PANEL
      SUBROUTINE DMUMPS_GET_OOC_PERM_PTR(TYPEF,
     &     NBPANELS,
     &     I_PIVPTR, I_PIV, IPOS, IW, LIW)
      USE MUMPS_OOC_COMMON ! To access TYPEF_L and TYPEF_U
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
C
C     Purpose:
C     =======
C
C     Get the pointers in IW on pivoting information to be stored
C     during factorization and used during the solve phase. This
C     routine is both for the symmetric (TYPEF=TYPEF_L) and unsymmetric
C     cases (TYPEF=TYPEF_L or TYPEF_U).
C     The total size of this space is estimated during
C     fac_ass.F / fac_ass_ELT.F and must be:
C     * Symmetric case: 1 for NASS + 1 for NBPANELS_L + NBPANELS_L + NASS
C     * Unsymmetric case: 1 + (1+NBPANELS_L+NASS) + (1+NBPANELS_U+NASS)
C     Size computation is in routine DMUMPS_OOC_GET_PP_SIZES.
C
C     At the end of the standard description of the structure of a node
C     (header, nb slaves, <slaves_list>, row indices, col indices), we
C     add, when panel version with pivoting is used:
C
C     NASS (nb of fully summed variables)
C     NBPANELS_L
C     PIVRPTR(1:NBPANELS_L)
C     PIV_L     (1:NASS)             NASS (=IW(IPOS)(or NASS-PIVRPTR(1) in
C     the future, after compression)
C     NBPANELS_U
C     PIVRPTR(1:NBPANELS_U)
C     PIV_U     (1:NASS)             NASS (=IW(IPOS)(or NASS-PIVRPTR(1) in
C     the future, after compression)
C
C
C     Output parameters:
C     =================
C     NBPANELS : nb of panels as estimated during assembly
C     I_PIVPTR : position in  IW of the starting of the pointer list
C     (of size NBPANELS) of the pointers to the list of pivots
C     I_PIV    : position in  IW of the starting of the pivot permutation list
C
      INTEGER, intent(out) :: NBPANELS, I_PIVPTR, I_PIV
      INTEGER, intent(in) :: TYPEF ! TYPEF_L or TYPEF_U
      INTEGER, intent(in) :: LIW, IPOS
      INTEGER IW(LIW)
C     Locals
      INTEGER I_NBPANELS, I_NASS
C
      I_NASS       = IPOS
      I_NBPANELS   = I_NASS + 1 ! L
      NBPANELS     = IW(I_NBPANELS) ! L
      I_PIVPTR     = I_NBPANELS + 1 ! L
      I_PIV        = I_PIVPTR + NBPANELS ! L
C     ... of size NASS = IW(I_NASS)
      IF (TYPEF==TYPEF_U) THEN
         I_NBPANELS   = I_PIV+IW(I_NASS) ! U
         NBPANELS     = IW(I_NBPANELS) ! U
         I_PIVPTR     = I_NBPANELS + 1 ! U
         I_PIV        = I_PIVPTR + NBPANELS ! U
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_GET_OOC_PERM_PTR
      SUBROUTINE DMUMPS_OOC_PP_SET_PTR(K50,NBPANELS_L,NBPANELS_U,
     &     NASS, IPOS, IW, LIW )
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Initialize the contents of PIV/PIVPTR/etc. that will store
C     pivoting information during the factorization.
C     NASS and NBPANELS are recorded. PIVPTR(1:NBPANELS)
C     is initialized to NASS+1. This will be modified during
C     the factorization in cases where permutations have to
C     be performed during the solve phase.
C
C     Arguments:
C     =========
C
      INTEGER K50
      INTEGER IPOS, NASS, NBPANELS_L, NBPANELS_U, LIW
      INTEGER IW(LIW)
C
C     Local variables:
C     ===============
C
      INTEGER IPOS_U
C     Executable statements
      IF (K50.EQ.1) THEN
         WRITE(*,*) "Internal error: DMUMPS_OOC_PP_SET_PTR called"
      ENDIF
      IW(IPOS)=NASS
      IW(IPOS+1)=NBPANELS_L
      IW(IPOS+2:IPOS+1+NBPANELS_L)=NASS+1
      IF (K50 == 0) THEN
         IPOS_U=IPOS+2+NASS+NBPANELS_L
         IW(IPOS_U)=NBPANELS_U
         IW(IPOS_U+1:IPOS_U+NBPANELS_U)=NASS+1
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_OOC_PP_SET_PTR
      SUBROUTINE DMUMPS_OOC_PP_TRYRELEASE_SPACE (
     &     IWPOS, IOLDPS, IW, LIW, MonBloc, NFRONT, KEEP
     &     )
      USE DMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
C
C     Purpose:
C     =======
C     If space used was at the top of the stack then
C     try to free space by detecting that
C     no permutation needs to be applied during
C     solve on panels.
C     One position is left (I_NASS) and set to -1
C     to indicate that permutation not needed at solve.
C
C     Arguments:
C     =========
C
      INTEGER, INTENT(IN)    :: IOLDPS, LIW, NFRONT,
     &     KEEP(500)
      INTEGER, INTENT(INOUT) :: IWPOS, IW(LIW)
      TYPE(IO_BLOCK), INTENT(IN):: MonBloc
C
C     Local variables:
C     ===============
C
      INTEGER :: NBPANELS_L,I_PIVRPTR_L, I_PIVR_L, NBPANELS_U,
     &     I_PIVRPTR_U, I_PIVR_U, XSIZE, IBEGOOC
      LOGICAL FREESPACE    ! set to true when permutation not needed
C     Executable statements
      IF (KEEP(50).EQ.1) RETURN ! no pivoting
C     --------------------------------
C     quick return if record is not at
C     the top of stack of L factors
      IF ((IOLDPS+IW(IOLDPS+XXI)).NE.IWPOS) RETURN
C     ---------------------------------------------
C     Panel+pivoting: get pointers on each subarray
C     ---------------------------------------------
      XSIZE   = KEEP(IXSZ)
      IBEGOOC = IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE
C     -- get L related data
      CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L,
     &     I_PIVRPTR_L, I_PIVR_L,
     &     IBEGOOC, IW, LIW)
      FREESPACE =
     &     (MonBloc%LastPiv.EQ.(IW(I_PIVRPTR_L)-1))
      IF (KEEP(50).EQ.0) THEN
C     -- get U related dataA
         CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_U, NBPANELS_U,
     &        I_PIVRPTR_U, I_PIVR_U,
     &        IBEGOOC, IW, LIW)
         FREESPACE =  FREESPACE .AND.
     &        (MonBloc%LastPiv.EQ.(IW(I_PIVRPTR_U)-1))
      ENDIF
C     ---------------------------------
C     Check if permutations eed be
C     performed on panels during solve
C     --------------------------------
      IF (FREESPACE) THEN
C     -- compress memory for that node: keep one entry set to -7777
         IW(IBEGOOC) = -7777    ! will be tested during solve
         IW(IOLDPS+XXI) = IBEGOOC
     &        - IOLDPS + 1      ! new size of inode's record
         IWPOS = IBEGOOC+1      ! move back to top of stack
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_OOC_PP_TRYRELEASE_SPACE
C
      SUBROUTINE DMUMPS_OOC_GET_PP_SIZES(K50, NBROW_L, NBCOL_U, NASS,
     &     NBPANELS_L, NBPANELS_U, LREQ)
      USE DMUMPS_OOC       ! To call DMUMPS_OOC_PANEL_SIZE
      IMPLICIT NONE
C
C     Purpose
C     =======
C
C     Compute the size of the workspace required to store the permutation
C     information during factorization, so that solve can permute back
C     what has to be permuted (this could not be done during factorization
C     because it was already on disk).
C
C     Arguments
C     =========
C
      INTEGER, intent(IN)  :: K50, NBROW_L, NBCOL_U, NASS
      INTEGER, intent(OUT) :: NBPANELS_L, NBPANELS_U, LREQ
      NBPANELS_L=-99999
      NBPANELS_U=-99999
C
C     Quick return in SPD case (no pivoting)
C
      IF (K50.EQ.1) THEN
         LREQ = 0
         RETURN
      ENDIF
C
C     L information is always computed
C
      NBPANELS_L = (NASS / DMUMPS_OOC_PANEL_SIZE(NBROW_L))+1
      LREQ =    1               ! Store NASS
     &     + 1                  ! Store NBPANELS_L
     &     + NASS               ! Store permutations
     &     + NBPANELS_L         ! Store pointers on permutations
      IF (K50.eq.0) THEN
C
C     Also take U information into account
C
         NBPANELS_U = (NASS / DMUMPS_OOC_PANEL_SIZE(NBCOL_U) ) +1
         LREQ = LREQ + 1        ! Store NBPANELS_U
     &        + NASS            ! Store permutations
     &        + NBPANELS_U      ! Store pointers on permutations
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_OOC_GET_PP_SIZES
      SUBROUTINE DMUMPS_OOC_PP_CHECK_PERM_FREED
     &           (IW_LOCATION, MUST_BE_PERMUTED)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IW_LOCATION
      LOGICAL, INTENT(INOUT) :: MUST_BE_PERMUTED
C
C     Purpose
C     =======
C
C     Reset MUST_BE_PERMUTED to .FALSE. when we detect
C     that the DMUMPS_OOC_PP_TRY_RELEASE_SPACE has freed
C     the permutation information (see that routine).
C
      IF (IW_LOCATION .EQ. -7777) THEN
        MUST_BE_PERMUTED = .FALSE.
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_OOC_PP_CHECK_PERM_FREED
