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
      SUBROUTINE MUMPS_ESTIM_FLOPS( INODE, N, PROCNODE_STEPS,
     &           KEEP199,
     &           ND, FILS, FRERE_STEPS, STEP, PIMASTER,
     &           KEEP28, KEEP50, KEEP253,
     &           FLOP1,
     &           IW, LIW, XSIZE )
      IMPLICIT NONE
      INTEGER INODE, N, KEEP50, LIW, KEEP199, KEEP28, KEEP253
      INTEGER PROCNODE_STEPS(KEEP28), ND(KEEP28),
     &        FILS(N), FRERE_STEPS(KEEP28),
     &        STEP(N), 
     & PIMASTER(KEEP28),
     &  IW( LIW )
      INTEGER XSIZE
      DOUBLE PRECISION FLOP1
      INTEGER NUMORG, IN, NASS, IFSON, NUMSTK, NFRONT, NPIV, NCB, 
     &        LEVEL, ISON
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      INTEGER MUMPS_TYPENODE
      EXTERNAL MUMPS_IN_OR_ROOT_SSARBR, MUMPS_TYPENODE
      INCLUDE 'mumps_headers.h'
      FLOP1 = 0.0D0
      IF (MUMPS_IN_OR_ROOT_SSARBR(PROCNODE_STEPS(STEP(INODE)),
     &                KEEP199) ) RETURN
      IN     = INODE
      NUMORG = 0
   10 NUMORG = NUMORG + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 10
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .EQ. 0) GOTO 30
   20 NUMSTK = NUMSTK + 1
      NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 +XSIZE)
      ISON = FRERE_STEPS(STEP(ISON))
      IF (ISON .GT. 0) GOTO 20
   30 NFRONT = ND(STEP(INODE)) + NASS + KEEP253
      NPIV  = NASS + NUMORG
      NCB   = NFRONT - NPIV
      LEVEL = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),KEEP199)
      CALL MUMPS_GET_FLOPS_COST(NFRONT,NPIV,NPIV,KEEP50,LEVEL,FLOP1)
      RETURN
      END SUBROUTINE MUMPS_ESTIM_FLOPS
      SUBROUTINE MUMPS_UPDATE_FLOPS_ROOT(OPELIW, KEEP50, NFRONT, NPIV,
     &           NPROW, NPCOL, MYID)
      DOUBLE PRECISION, intent(inout) :: OPELIW
      INTEGER, intent(in) :: KEEP50, NFRONT, NPIV,
     &           NPROW, NPCOL, MYID
      DOUBLE PRECISION :: COST, COST_PER_PROC
      INTEGER, PARAMETER :: LEVEL3 = 3
      CALL MUMPS_GET_FLOPS_COST(NFRONT, NPIV, NFRONT, KEEP50, LEVEL3, 
     &                          COST)
      COST_PER_PROC = dble(int( COST,8) / int(NPROW * NPCOL,8))
      OPELIW = OPELIW + COST_PER_PROC
      RETURN
      END SUBROUTINE MUMPS_UPDATE_FLOPS_ROOT
      SUBROUTINE MUMPS_GET_FLOPS_COST(NFRONT,NPIV,NASS,
     &                                 KEEP50,LEVEL,COST)
      IMPLICIT NONE
      INTEGER, intent(in) :: NFRONT,NPIV,KEEP50,LEVEL, NASS
      DOUBLE PRECISION, intent(out) :: COST
      IF (KEEP50.EQ.0) THEN
        IF (LEVEL.EQ.1 .OR. LEVEL.EQ.3) THEN
          COST = dble(2) * dble(NFRONT) * dble(NPIV) *
     &      dble(NFRONT - NPIV - 1) +
     &      dble(NPIV) * dble(NPIV + 1) * dble(2 * NPIV + 1)
     &          / dble(3)
          COST = COST + dble(2 * NFRONT - NPIV - 1)
     &      * dble(NPIV) /dble(2)
        ELSEIF (LEVEL.EQ.2) THEN
          COST = dble(2*NASS)*dble(NFRONT) -
     &          dble(NASS+NFRONT)*dble(NPIV+1)
          COST = dble(NPIV)*COST +
     &     dble(2 * NASS - NPIV - 1) * dble(NPIV) / dble(2) +
     &     dble(NPIV) * dble(NPIV + 1) *
     &     dble(2 * NPIV + 1) /dble(3)
        ENDIF
      ELSE
        IF (LEVEL.EQ.1 .OR. (LEVEL.EQ.3 .AND. KEEP50.EQ.1)) THEN
          COST = dble(NPIV) * (
     &          dble( NFRONT ) * dble( NFRONT ) +
     &          dble( NFRONT ) - (
     &          dble( NFRONT)*dble(NPIV) + dble(NPIV+1)
     &          )) +( dble(NPIV)*dble(NPIV+1)
     &          *dble(2*NPIV+1))/ dble(6)
        ELSE IF (LEVEL.EQ.3.AND.KEEP50.EQ.2) THEN
          COST = dble(2) * dble(NFRONT) * dble(NPIV) *
     &      dble(NFRONT - NPIV - 1) +
     &      dble(NPIV) * dble(NPIV + 1) *
     &      dble(2 * NPIV + 1) / dble(3)
          COST = COST + dble(2 * NFRONT - NPIV - 1)
     &         * dble(NPIV) / dble(2)
        ELSE
          COST = dble(NPIV) * (
     &          dble( NASS ) * dble( NASS ) + dble( NASS )
     &        - ( dble( NASS) * dble(NPIV) + dble( NPIV + 1 ) ) )
     &        + ( dble(NPIV)*dble(NPIV+1)*dble(2*NPIV+1) )
     &        / dble( 6 )
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_GET_FLOPS_COST
      SUBROUTINE MUMPS_PRINT_STILL_ACTIVE(MYID, KEEP, DKEEP17, OPELIW, 
     &                        OPLAST_PRINTED, MPA)
      IMPLICIT NONE
      INTEGER, intent(in)             :: KEEP (500), MYID, MPA
      DOUBLE PRECISION                :: DKEEP17
      DOUBLE PRECISION, intent(in)    :: OPELIW
      DOUBLE PRECISION, intent(inout) :: OPLAST_PRINTED
      IF (MPA.GT.0) THEN
        IF ( (OPELIW-OPLAST_PRINTED).GT. DKEEP17) THEN
         WRITE(MPA,'(A,I6,A,A,1PD10.3)')
     &   ' ... MPI process', MYID, 
     &   ': theoretical number of flops locally performed',
     &   ' so far        = ',
     &    OPELIW
         OPLAST_PRINTED = OPELIW
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_PRINT_STILL_ACTIVE
