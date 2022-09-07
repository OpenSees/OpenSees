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
      SUBROUTINE DMUMPS_PROCESS_RTNELIND( ROOT, 
     &    INODE, NELIM, NSLAVES, ROW_LIST,
     &    COL_LIST, SLAVE_LIST, 
     &
     &    PROCNODE_STEPS, IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S,
     &    ITLOC, RHS_MUMPS, COMP,
     &    IFLAG, IERROR, 
     &    IPOOL, LPOOL, LEAF, MYID, SLAVEF,
     &    KEEP, KEEP8, DKEEP,
     &    COMM, COMM_LOAD, FILS, DAD, ND )
      USE DMUMPS_LOAD
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: ROOT
      INTEGER INODE, NELIM, NSLAVES 
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION DKEEP(230)
      INTEGER ROW_LIST(*), COL_LIST(*), 
     &        SLAVE_LIST(*)
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER PTRIST( KEEP(28) ), PTLUST_S(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), ITLOC( N + KEEP(253) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER IFLAG, IERROR
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF
      INTEGER COMM, COMM_LOAD, ND(KEEP(28)), FILS(N), DAD(KEEP(28))
      INTEGER IROOT, TYPE_INODE, DEB_ROW, DEB_COL,
     &        NOINT
      INTEGER(8) :: NOREAL
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MUMPS_TYPENODE
      EXTERNAL MUMPS_TYPENODE
      IROOT        = KEEP(38)
      NSTK_S(STEP(IROOT))= NSTK_S(STEP(IROOT)) - 1
      KEEP(42) = KEEP(42) + NELIM
      TYPE_INODE= MUMPS_TYPENODE( PROCNODE_STEPS(STEP(INODE)),
     &                            KEEP(199) )
      IF (TYPE_INODE.EQ.1) THEN 
        IF (NELIM.EQ.0) THEN
         KEEP(41) = KEEP(41) + 1
        ELSE 
         KEEP(41) = KEEP(41) + 3
        ENDIF
      ELSE
        IF (NELIM.EQ.0) THEN
         KEEP(41) = KEEP(41) + NSLAVES
        ELSE 
         KEEP(41) = KEEP(41) + 2*NSLAVES + 1
        ENDIF
      ENDIF
      IF  (NELIM.EQ.0) THEN
        PIMASTER(STEP(INODE)) = 0 
      ELSE
       NOINT = 6 + NSLAVES + NELIM  + NELIM + KEEP(IXSZ)
       NOREAL= 0_8
       CALL DMUMPS_ALLOC_CB(.FALSE.,0_8,.FALSE.,.FALSE.,
     &   MYID,N,KEEP,KEEP8,DKEEP,IW,LIW, A, LA,
     &   LRLU, IPTRLU,IWPOS,IWPOSCB, SLAVEF, PROCNODE_STEPS, DAD,
     &   PTRIST,PTRAST,STEP, PIMASTER, PAMASTER,
     &   NOINT, NOREAL, INODE, S_NOTFREE, .TRUE.,
     &   COMP, LRLUS, KEEP8(67), IFLAG, IERROR
     &      )
       IF ( IFLAG .LT. 0 ) THEN
         WRITE(*,*) ' Failure in int space allocation in CB area ',
     &    ' during assembly of root : DMUMPS_PROCESS_RTNELIND',
     &    ' size required was :', NOINT,
     &    'INODE=',INODE,' NELIM=',NELIM, ' NSLAVES=', NSLAVES
         RETURN
        ENDIF
        PIMASTER(STEP( INODE )) = IWPOSCB + 1
        PAMASTER(STEP( INODE )) = IPTRLU  + 1_8
        IW( IWPOSCB + 1+KEEP(IXSZ) ) = 2*NELIM
        IW( IWPOSCB + 2+KEEP(IXSZ) ) = NELIM
        IW( IWPOSCB + 3+KEEP(IXSZ) ) = 0
        IW( IWPOSCB + 4+KEEP(IXSZ) ) = 0
        IW( IWPOSCB + 5+KEEP(IXSZ) ) = 1
        IW( IWPOSCB + 6+KEEP(IXSZ) ) = NSLAVES
        IF (NSLAVES.GT.0) THEN
         IW( IWPOSCB+7+KEEP(IXSZ):IWPOSCB+7+KEEP(IXSZ)+NSLAVES-1) = 
     &                   SLAVE_LIST(1:NSLAVES)
        ENDIF
        DEB_ROW = IWPOSCB+7+NSLAVES+KEEP(IXSZ)
        IW(DEB_ROW : DEB_ROW+NELIM -1) = ROW_LIST(1:NELIM)
        DEB_COL = DEB_ROW + NELIM
        IW(DEB_COL : DEB_COL+NELIM -1) = COL_LIST(1:NELIM)
      ENDIF
      IF (NSTK_S(STEP(IROOT)) .EQ. 0 ) THEN
          CALL DMUMPS_INSERT_POOL_N(N, IPOOL, LPOOL, PROCNODE_STEPS,
     &         SLAVEF, KEEP(199),
     &         KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     &         STEP, IROOT )
          IF (KEEP(47) .GE. 3) THEN
             CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &            IPOOL, LPOOL, 
     &            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &            MYID, STEP, N, ND, FILS )
          ENDIF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_RTNELIND
