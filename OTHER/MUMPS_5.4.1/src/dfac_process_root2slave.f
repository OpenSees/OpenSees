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
      SUBROUTINE DMUMPS_PROCESS_ROOT2SLAVE( TOT_ROOT_SIZE,
     &    TOT_CONT_TO_RECV, root,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, COMM_LOAD,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP, KEEP8, DKEEP, ND)
      USE DMUMPS_LOAD
      USE DMUMPS_OOC        
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER KEEP(500), ICNTL(60)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION DKEEP(230)
      INTEGER TOT_ROOT_SIZE, TOT_CONT_TO_RECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA, POSFAC
      INTEGER(8) :: PTRFAC(KEEP(28)), PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER ND( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, COMM_LOAD
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC(N+KEEP(253)), FILS(N), DAD(KEEP(28))
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER(8), INTENT(IN) :: PTRARW(LPTRAR), PTRAIW(LPTRAR)
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(KEEP8(27))
      DOUBLE PRECISION DBLARR(KEEP8(26))
      INTEGER ::  allocok
      DOUBLE PRECISION, DIMENSION(:,:), POINTER :: TMP
      INTEGER NEW_LOCAL_M, NEW_LOCAL_N
      INTEGER OLD_LOCAL_M, OLD_LOCAL_N
      INTEGER I, J
      INTEGER LREQI, IROOT
      INTEGER(8) :: LREQA
      INTEGER POSHEAD, IPOS_SON,IERR
      LOGICAL MASTER_OF_ROOT, NO_OLD_ROOT
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INCLUDE 'mumps_headers.h'
      INTEGER numroc, MUMPS_PROCNODE
      EXTERNAL numroc, MUMPS_PROCNODE
      IROOT = KEEP( 38 )
      root%TOT_ROOT_SIZE = TOT_ROOT_SIZE
      MASTER_OF_ROOT = ( MYID .EQ. 
     &                   MUMPS_PROCNODE( PROCNODE_STEPS(STEP(IROOT)),
     &                   KEEP(199) ) )
      NEW_LOCAL_M  = numroc( TOT_ROOT_SIZE, root%MBLOCK,
     &               root%MYROW, 0, root%NPROW )
      NEW_LOCAL_M  = max( 1, NEW_LOCAL_M )
      NEW_LOCAL_N  = numroc( TOT_ROOT_SIZE, root%NBLOCK,
     &               root%MYCOL, 0, root%NPCOL )
      IF ( PTRIST(STEP( IROOT )).GT.0) THEN
        OLD_LOCAL_N = -IW( PTRIST(STEP( IROOT )) + KEEP(IXSZ) )
        OLD_LOCAL_M =  IW( PTRIST(STEP( IROOT )) + 1  + KEEP(IXSZ))
      ELSE
        OLD_LOCAL_N = 0
        OLD_LOCAL_M = NEW_LOCAL_M
      ENDIF
      IF (PTRIST(STEP(IROOT)) .EQ.0) THEN
        NO_OLD_ROOT = .TRUE.
      ELSE
        NO_OLD_ROOT =.FALSE.
      ENDIF
      IF (KEEP(60) .NE. 0) THEN
        IF ( MASTER_OF_ROOT ) THEN
          LREQI=6+2*TOT_ROOT_SIZE+KEEP(IXSZ)
          LREQA=0_8
          IF ( IWPOS + LREQI - 1. GT. IWPOSCB ) THEN
           CALL DMUMPS_COMPRE_NEW( N, KEEP(28), IW, LIW, A, LA,
     &           LRLU, IPTRLU,
     &           IWPOS, IWPOSCB, PTRIST, PTRAST,
     &           STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &           KEEP(IXSZ),COMP,DKEEP(97),
     &           MYID, SLAVEF, KEEP(199), PROCNODE_STEPS, DAD )
           IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB1 compress root2slave:LRLU,LRLUS=',
     &            LRLU, LRLUS
                  IFLAG = -9
                  CALL MUMPS_SET_IERROR(LREQA-LRLUS, IERROR)
                  GOTO 700
           END IF
          ENDIF
          IF ( IWPOS + LREQI - 1. GT. IWPOSCB ) THEN
            IFLAG = -8
            IERROR = IWPOS + LREQI - 1 - IWPOSCB
            GOTO 700
          ENDIF
          PTLUST(STEP(IROOT))= IWPOS
          IWPOS = IWPOS + LREQI
          POSHEAD = PTLUST( STEP(IROOT))
          IW( POSHEAD + XXI )=LREQI
          CALL MUMPS_STOREI8( LREQA, IW(POSHEAD + XXR) )
          CALL MUMPS_STOREI8( 0_8, IW(POSHEAD + XXD) )
          IW( POSHEAD + XXS )=-9999
          IW(POSHEAD+XXS+1:POSHEAD+KEEP(IXSZ)-1)=-99999
          IW( POSHEAD +KEEP(IXSZ)) = 0
          IW( POSHEAD + 1 +KEEP(IXSZ)) = -1
          IW( POSHEAD + 2 +KEEP(IXSZ)) = -1
          IW( POSHEAD + 4 +KEEP(IXSZ)) = STEP(IROOT)
          IW( POSHEAD + 5 +KEEP(IXSZ)) = 0
          IW( POSHEAD + 3 +KEEP(IXSZ)) = TOT_ROOT_SIZE
        ELSE 
          PTLUST(STEP(IROOT)) = -4444
        ENDIF
        PTRIST(STEP(IROOT)) = 0
        PTRFAC(STEP(IROOT)) = -4445_8
        IF (root%yes .and. NO_OLD_ROOT) THEN
          IF (NEW_LOCAL_N .GT. 0) THEN
            CALL DMUMPS_SET_TO_ZERO(root%SCHUR_POINTER(1),
     &      root%SCHUR_LLD, root%SCHUR_MLOC, root%SCHUR_NLOC,
     &      KEEP)
            IF (KEEP(55).EQ.0) THEN
              CALL DMUMPS_ASM_ARR_ROOT( N, root, IROOT,
     &        root%SCHUR_POINTER(1), root%SCHUR_LLD, root%SCHUR_MLOC,
     &        root%SCHUR_NLOC, FILS, PTRAIW, PTRARW, INTARR, DBLARR,
     &        KEEP8(27), KEEP8(26), MYID )
            ELSE
              CALL DMUMPS_ASM_ELT_ROOT(N, root,
     &        root%SCHUR_POINTER(1), root%SCHUR_LLD, root%SCHUR_MLOC,
     &        root%SCHUR_NLOC, LPTRAR, NELT, FRTPTR, FRTELT,
     &        PTRAIW, PTRARW, INTARR, DBLARR,
     &        KEEP8(27), KEEP8(26), KEEP, KEEP8, MYID )
            ENDIF
          ENDIF
        ENDIF
      ELSE
        IF ( MASTER_OF_ROOT ) THEN
          LREQI = 6 + 2 * TOT_ROOT_SIZE+KEEP(IXSZ)
        ELSE
          LREQI = 6+KEEP(IXSZ)
        END IF
        LREQA = int(NEW_LOCAL_M, 8) * int(NEW_LOCAL_N, 8)
        CALL DMUMPS_GET_SIZE_NEEDED( 
     &           LREQI , LREQA, .FALSE.,
     &           KEEP(1), KEEP8(1),
     &           N, KEEP(28), IW, LIW, A, LA,
     &           LRLU, IPTRLU,
     &           IWPOS, IWPOSCB, PTRIST, PTRAST,
     &           STEP, PIMASTER, PAMASTER, KEEP(216), LRLUS,
     &           KEEP(IXSZ), COMP, DKEEP(97),
     &           MYID, SLAVEF, PROCNODE_STEPS, DAD, 
     &           IFLAG, IERROR )
        IF (IFLAG.LT.0) GOTO 700
        PTLUST(STEP( IROOT )) = IWPOS
        IWPOS           = IWPOS + LREQI
        IF (LREQA.EQ.0_8) THEN
          PTRAST (STEP(IROOT)) = POSFAC
          PTRFAC (STEP(IROOT)) = POSFAC
        ELSE
          PTRAST (STEP(IROOT)) = POSFAC
          PTRFAC (STEP(IROOT)) = POSFAC
        ENDIF
        POSFAC           = POSFAC + LREQA
        LRLU   = LRLU  - LREQA
        LRLUS  = LRLUS - LREQA
        KEEP8(67) = min(KEEP8(67), LRLUS)
        KEEP8(69) = KEEP8(69) + LREQA
        KEEP8(68) = max(KEEP8(69), KEEP8(68))
        CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &            LA-LRLUS,0_8,LREQA,KEEP,KEEP8,LRLUS)
        POSHEAD = PTLUST( STEP(IROOT))
        IW( POSHEAD + XXI )     = LREQI
        CALL MUMPS_STOREI8( LREQA, IW(POSHEAD + XXR))
        CALL MUMPS_STOREI8( 0_8, IW(POSHEAD + XXD))
        IW( POSHEAD + XXS ) = S_NOTFREE
        IW(POSHEAD+XXS+1:POSHEAD+KEEP(IXSZ)-1)=-99999
        IW( POSHEAD + KEEP(IXSZ) ) = 0
        IW( POSHEAD + 1 + KEEP(IXSZ) ) = NEW_LOCAL_N
        IW( POSHEAD + 2 + KEEP(IXSZ) ) = NEW_LOCAL_M
        IW( POSHEAD + 4 + KEEP(IXSZ) ) = STEP(IROOT)
        IW( POSHEAD + 5 + KEEP(IXSZ) ) = 0
        IF ( MASTER_OF_ROOT ) THEN
          IW( POSHEAD + 3 + KEEP(IXSZ) ) = TOT_ROOT_SIZE
        ELSE
          IW( POSHEAD + 3 + KEEP(IXSZ) ) = 0
        ENDIF
        IF ( PTRIST(STEP(IROOT)) .EQ. 0) THEN
          CALL DMUMPS_SET_TO_ZERO(A(PTRAST(STEP(IROOT))),
     &    NEW_LOCAL_M, NEW_LOCAL_M, NEW_LOCAL_N, KEEP)
          IF (KEEP(55) .EQ.0 ) THEN
            CALL DMUMPS_ASM_ARR_ROOT( N, root, IROOT,
     &      A(PTRAST(STEP(IROOT))),
     &      NEW_LOCAL_M, NEW_LOCAL_M, NEW_LOCAL_N,
     &      FILS, PTRAIW, PTRARW, INTARR, DBLARR,
     &      KEEP8(27), KEEP8(26), MYID )
          ELSE
            CALL DMUMPS_ASM_ELT_ROOT( N, root,
     &      A(PTRAST(STEP(IROOT))),
     &      NEW_LOCAL_M, NEW_LOCAL_M, NEW_LOCAL_N,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &      PTRAIW, PTRARW, INTARR, DBLARR,
     &      KEEP8(27), KEEP8(26), KEEP, KEEP8, MYID )
          ENDIF
          PAMASTER(STEP(IROOT)) = 0_8
        ELSE IF ( PTRIST(STEP(IROOT)) .LT. 0 ) THEN
          CALL DMUMPS_SET_TO_ZERO(A(PTRAST(STEP(IROOT))),
     &    NEW_LOCAL_M, NEW_LOCAL_M, NEW_LOCAL_N, KEEP)
        ELSE
          OLD_LOCAL_N = -IW( PTRIST(STEP( IROOT )) + KEEP(IXSZ) )
          OLD_LOCAL_M =  IW( PTRIST(STEP( IROOT )) + 1  + KEEP(IXSZ))
          IF ( TOT_ROOT_SIZE .eq. root%ROOT_SIZE ) THEN
            IF ( LREQA .NE. int(OLD_LOCAL_M,8) * int(OLD_LOCAL_N,8) )
     &            THEN
              write(*,*) 'error 1 in PROCESS_ROOT2SLAVE',
     &        OLD_LOCAL_M, OLD_LOCAL_N
              CALL MUMPS_ABORT()
            END IF
            CALL DMUMPS_COPYI8SIZE(LREQA,
     &                            A( PAMASTER(STEP(IROOT)) ),
     &                            A( PTRAST  (STEP(IROOT)) ) )
          ELSE
            CALL DMUMPS_COPY_ROOT( A( PTRAST(STEP(IROOT))), 
     &          NEW_LOCAL_M,
     &          NEW_LOCAL_N, A( PAMASTER( STEP(IROOT)) ), OLD_LOCAL_M,
     &          OLD_LOCAL_N )
          END IF
          IF ( PTRIST( STEP( IROOT ) ) .GT. 0 ) THEN
             IPOS_SON= PTRIST( STEP(IROOT))
             CALL DMUMPS_FREE_BLOCK_CB_STATIC(.FALSE.,
     &            MYID, N, IPOS_SON,
     &            IW, LIW, LRLU, LRLUS, IPTRLU,
     &            IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &           )
          END IF
        ENDIF 
        PTRIST(STEP( IROOT ))   = 0
        PAMASTER(STEP( IROOT )) = 0_8
      ENDIF 
      IF ( NO_OLD_ROOT ) THEN
          IF (KEEP(253) .GT.0) THEN
            root%RHS_NLOC = numroc( KEEP(253), root%NBLOCK,
     &                      root%MYCOL, 0, root%NPCOL )
            root%RHS_NLOC = max( root%RHS_NLOC, 1 )
          ELSE
            root%RHS_NLOC = 1
          ENDIF
          IF (associated(root%RHS_ROOT)) DEALLOCATE(root%RHS_ROOT)
          ALLOCATE(root%RHS_ROOT(NEW_LOCAL_M, root%RHS_NLOC),
     &              stat=allocok)
          IF ( allocok.GT.0 ) THEN
             IFLAG = -13
             IERROR = NEW_LOCAL_N * root%RHS_NLOC
            GOTO 700
          ENDIF
          IF (KEEP(253) .NE. 0) THEN
            root%RHS_ROOT=ZERO
            CALL DMUMPS_ASM_RHS_ROOT( N, FILS, root, KEEP, RHS_MUMPS,
     &      IFLAG, IERROR )
          ENDIF
      ELSE IF (NEW_LOCAL_M.GT.OLD_LOCAL_M .AND. KEEP(253) .GT.0) THEN
          TMP => root%RHS_ROOT
          NULLIFY(root%RHS_ROOT)
          ALLOCATE (root%RHS_ROOT(NEW_LOCAL_M, root%RHS_NLOC), 
     &                stat=allocok)
          IF ( allocok.GT.0) THEN
              IFLAG=-13
              IERROR = NEW_LOCAL_M*root%RHS_NLOC
              GOTO 700
          ENDIF
          DO J = 1, root%RHS_NLOC
            DO I = 1, OLD_LOCAL_M
              root%RHS_ROOT(I,J)=TMP(I,J)
            ENDDO
            DO I = OLD_LOCAL_M+1, NEW_LOCAL_M
              root%RHS_ROOT(I,J) = ZERO
            ENDDO
          ENDDO
          DEALLOCATE(TMP)
          NULLIFY(TMP) 
      ENDIF
      KEEP(121) = KEEP(121) + TOT_CONT_TO_RECV
      IF ( KEEP(121) .eq. 0 ) THEN
         IF (KEEP(201).EQ.1) THEN 
            CALL DMUMPS_OOC_FORCE_WRT_BUF_PANEL(IERR)
         ELSE IF (KEEP(201).EQ.2) THEN 
            CALL DMUMPS_FORCE_WRITE_BUF(IERR)              
         ENDIF
        CALL DMUMPS_INSERT_POOL_N( N, IPOOL, LPOOL, PROCNODE_STEPS,
     &       SLAVEF, KEEP(199), KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     &       STEP, IROOT + N )
        IF (KEEP(47) .GE. 3) THEN
           CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &          IPOOL, LPOOL, 
     &          PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &          MYID, STEP, N, ND, FILS )
        ENDIF
      END IF
      RETURN
 700  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_ROOT2SLAVE
      SUBROUTINE DMUMPS_COPY_ROOT
     &( NEW, M_NEW, N_NEW,OLD, M_OLD, N_OLD )
      INTEGER M_NEW, N_NEW, M_OLD, N_OLD
      DOUBLE PRECISION NEW( M_NEW, N_NEW ), OLD( M_OLD, N_OLD )
      INTEGER J
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      DO J = 1, N_OLD
        NEW( 1: M_OLD, J ) = OLD( 1: M_OLD, J )
        NEW( M_OLD + 1: M_NEW, J ) = ZERO
      END DO
      NEW( 1: M_NEW,N_OLD + 1: N_NEW ) = ZERO
      RETURN
      END SUBROUTINE DMUMPS_COPY_ROOT
