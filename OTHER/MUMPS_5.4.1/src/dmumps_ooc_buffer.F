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
      MODULE DMUMPS_OOC_BUFFER
      USE MUMPS_OOC_COMMON
      IMPLICIT NONE
      PUBLIC
      INTEGER FIRST_HBUF,SECOND_HBUF
      PARAMETER (FIRST_HBUF=0, SECOND_HBUF=1)
      INTEGER,SAVE :: OOC_FCT_TYPE_LOC
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: BUF_IO
      LOGICAL,SAVE :: PANEL_FLAG
      INTEGER,SAVE :: EARLIEST_WRITE_MIN_SIZE
      INTEGER(8),SAVE,DIMENSION(:), ALLOCATABLE ::
     &  I_SHIFT_FIRST_HBUF, I_SHIFT_SECOND_HBUF,
     &  I_SHIFT_CUR_HBUF, I_REL_POS_CUR_HBUF
      INTEGER, SAVE, DIMENSION(:), ALLOCATABLE ::
     &  LAST_IOREQUEST, CUR_HBUF
      INTEGER, DIMENSION(:),ALLOCATABLE :: I_CUR_HBUF_NEXTPOS
      INTEGER,SAVE ::  I_CUR_HBUF_FSTPOS,
     &  I_SUB_HBUF_FSTPOS
      INTEGER(8) :: BufferEmpty
      PARAMETER (BufferEmpty=-1_8)
      INTEGER(8), DIMENSION(:),ALLOCATABLE :: NextAddVirtBuffer
      INTEGER(8), DIMENSION(:),ALLOCATABLE :: FIRST_VADDR_IN_BUF
      CONTAINS
      SUBROUTINE DMUMPS_OOC_NEXT_HBUF(TYPEF_ARG)
      IMPLICIT NONE
      INTEGER TYPEF_ARG
      SELECT CASE(CUR_HBUF(TYPEF_ARG))
         CASE (FIRST_HBUF)
            CUR_HBUF(TYPEF_ARG) = SECOND_HBUF
            I_SHIFT_CUR_HBUF(TYPEF_ARG) =
     &           I_SHIFT_SECOND_HBUF(TYPEF_ARG)
         CASE (SECOND_HBUF)
            CUR_HBUF(TYPEF_ARG) = FIRST_HBUF
            I_SHIFT_CUR_HBUF(TYPEF_ARG) =
     &           I_SHIFT_FIRST_HBUF(TYPEF_ARG)
      END SELECT
      IF(.NOT.PANEL_FLAG)THEN
         I_SUB_HBUF_FSTPOS =I_CUR_HBUF_FSTPOS
         I_CUR_HBUF_FSTPOS =I_CUR_HBUF_NEXTPOS(TYPEF_ARG)
      ENDIF
      I_REL_POS_CUR_HBUF(TYPEF_ARG) = 1_8
      RETURN
      END SUBROUTINE DMUMPS_OOC_NEXT_HBUF
      SUBROUTINE DMUMPS_OOC_DO_IO_AND_CHBUF(TYPEF_ARG,IERR)
      IMPLICIT NONE
      INTEGER TYPEF_ARG
      INTEGER NEW_IOREQUEST
      INTEGER IERR
      IERR=0
      CALL DMUMPS_OOC_WRT_CUR_BUF2DISK(TYPEF_ARG,NEW_IOREQUEST,
     &     IERR)
      IF(IERR.LT.0)THEN
         RETURN
      ENDIF
      IERR=0
      CALL MUMPS_WAIT_REQUEST(LAST_IOREQUEST(TYPEF_ARG),IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1>0)
     &   WRITE(ICNTL1,*) MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         RETURN
      ENDIF
      LAST_IOREQUEST(TYPEF_ARG) = NEW_IOREQUEST
      CALL DMUMPS_OOC_NEXT_HBUF(TYPEF_ARG)
      IF(PANEL_FLAG)THEN
         NextAddVirtBuffer(TYPEF_ARG)=BufferEmpty
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_OOC_DO_IO_AND_CHBUF
      SUBROUTINE DMUMPS_OOC_BUF_CLEAN_PENDING(IERR)
      IMPLICIT NONE
      INTEGER, intent(out) :: IERR
      INTEGER TYPEF_LAST        
      INTEGER TYPEF_LOC
      IERR = 0 
      TYPEF_LAST = OOC_NB_FILE_TYPE
      DO TYPEF_LOC = 1, TYPEF_LAST
         IERR=0
         CALL  DMUMPS_OOC_DO_IO_AND_CHBUF(TYPEF_LOC,IERR)
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
         IERR=0
         CALL DMUMPS_OOC_DO_IO_AND_CHBUF(TYPEF_LOC,IERR)
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_OOC_BUF_CLEAN_PENDING
      SUBROUTINE DMUMPS_OOC_WRT_CUR_BUF2DISK(TYPEF_ARG,IOREQUEST,
     &     IERR)
      IMPLICIT NONE
      INTEGER IOREQUEST,IERR
      INTEGER TYPEF_ARG
      INTEGER FIRST_INODE
      INTEGER(8) :: FROM_BUFIO_POS, SIZE
      INTEGER TYPE
      INTEGER ADDR_INT1,ADDR_INT2
      INTEGER(8) TMP_VADDR
      INTEGER SIZE_INT1,SIZE_INT2
      IERR=0
      IF (I_REL_POS_CUR_HBUF(TYPEF_ARG) == 1_8) THEN
        IOREQUEST=-1
        RETURN
      END IF      
      IF(PANEL_FLAG)THEN
         TYPE=TYPEF_ARG-1
         FIRST_INODE=-9999
         TMP_VADDR=FIRST_VADDR_IN_BUF(TYPEF_ARG)
      ELSE
         TYPE=FCT
         FIRST_INODE =
     &        OOC_INODE_SEQUENCE(I_CUR_HBUF_FSTPOS,TYPEF_ARG)
         TMP_VADDR=OOC_VADDR(STEP_OOC(FIRST_INODE),TYPEF_ARG)
      ENDIF
      FROM_BUFIO_POS=I_SHIFT_CUR_HBUF(TYPEF_ARG)+1_8
      SIZE = I_REL_POS_CUR_HBUF(TYPEF_ARG)-1_8
      CALL MUMPS_OOC_CONVERT_BIGINTTO2INT(ADDR_INT1,ADDR_INT2,
     &     TMP_VADDR)
      CALL MUMPS_OOC_CONVERT_BIGINTTO2INT(SIZE_INT1,SIZE_INT2,
     &     SIZE)
      CALL MUMPS_LOW_LEVEL_WRITE_OOC_C(LOW_LEVEL_STRAT_IO,
     &     BUF_IO(FROM_BUFIO_POS),SIZE_INT1,SIZE_INT2,
     &     FIRST_INODE,IOREQUEST,
     &     TYPE,ADDR_INT1,ADDR_INT2,IERR)
      IF(IERR.LT.0)THEN
         IF (ICNTL1>0)
     &   WRITE(ICNTL1,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         RETURN
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_OOC_WRT_CUR_BUF2DISK
      SUBROUTINE DMUMPS_INIT_OOC_BUF(I1,I2,IERR)
      IMPLICIT NONE
      INTEGER I1,I2,IERR 
      INTEGER allocok
      IERR=0
      PANEL_FLAG=.FALSE.
      IF(allocated(I_SHIFT_FIRST_HBUF))THEN
         DEALLOCATE(I_SHIFT_FIRST_HBUF)
      ENDIF
      IF(allocated(I_SHIFT_SECOND_HBUF))THEN
         DEALLOCATE(I_SHIFT_SECOND_HBUF)
      ENDIF
      IF(allocated(I_SHIFT_CUR_HBUF))THEN
         DEALLOCATE(I_SHIFT_CUR_HBUF)
      ENDIF
      IF(allocated(I_REL_POS_CUR_HBUF))THEN
         DEALLOCATE(I_REL_POS_CUR_HBUF)
      ENDIF
      IF(allocated(LAST_IOREQUEST))THEN
         DEALLOCATE(LAST_IOREQUEST)
      ENDIF
      IF(allocated(CUR_HBUF))THEN
         DEALLOCATE(CUR_HBUF)
      ENDIF
      DIM_BUF_IO = int(KEEP_OOC(100),8)
      ALLOCATE(I_SHIFT_FIRST_HBUF(OOC_NB_FILE_TYPE),
     &     stat=allocok)
      IF (allocok > 0) THEN
         IF (ICNTL1>0) THEN
            WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         ENDIF
         I1 = -13
         I2 = OOC_NB_FILE_TYPE
         IERR=-1
         RETURN
      ENDIF
      ALLOCATE(I_SHIFT_SECOND_HBUF(OOC_NB_FILE_TYPE),
     &     stat=allocok)
      IF (allocok > 0) THEN
         IF (ICNTL1>0) THEN
            WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         ENDIF
         I1 = -13
         I2 = OOC_NB_FILE_TYPE
         IERR=-1
         RETURN
      ENDIF
      ALLOCATE(I_SHIFT_CUR_HBUF(OOC_NB_FILE_TYPE),
     &     stat=allocok)
      IF (allocok > 0) THEN
         IF (ICNTL1>0) THEN
            WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         ENDIF
         I1 = -13
         I2 = OOC_NB_FILE_TYPE
         IERR=-1
         RETURN
      ENDIF
      ALLOCATE(I_REL_POS_CUR_HBUF(OOC_NB_FILE_TYPE),
     &     stat=allocok)
      IF (allocok > 0) THEN
         IF (ICNTL1>0) THEN
            WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         ENDIF
         I1 = -13
         I2 = OOC_NB_FILE_TYPE
         IERR=-1
         RETURN
      ENDIF
      ALLOCATE(LAST_IOREQUEST(OOC_NB_FILE_TYPE),
     &     stat=allocok)
      IF (allocok > 0) THEN
         IF (ICNTL1>0) THEN
            WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         ENDIF
         I1 = -13
         I2 = OOC_NB_FILE_TYPE
         IERR=-1
         RETURN
      ENDIF
      ALLOCATE(CUR_HBUF(OOC_NB_FILE_TYPE),
     &     stat=allocok)
      IF (allocok > 0) THEN
         IF (ICNTL1>0) THEN
            WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         ENDIF
         I1 = -13
         I2 = OOC_NB_FILE_TYPE
         IERR=-1
         RETURN
      ENDIF
      OOC_FCT_TYPE_LOC=OOC_NB_FILE_TYPE
      ALLOCATE(BUF_IO(DIM_BUF_IO), stat=allocok)      
      IF (allocok > 0) THEN
         IF (ICNTL1>0) THEN
            WRITE(ICNTL1,*) 'PB allocation in DMUMPS_INIT_OOC'
         ENDIF
         I1 = -13
         CALL MUMPS_SET_IERROR(DIM_BUF_IO, I2)
         RETURN
      ENDIF
      PANEL_FLAG=(KEEP_OOC(201).EQ.1)
      IF (PANEL_FLAG) THEN 
         IERR=0
         KEEP_OOC(228)=0
         IF(allocated(AddVirtLibre))THEN
            DEALLOCATE(AddVirtLibre)
         ENDIF
         ALLOCATE(AddVirtLibre(OOC_NB_FILE_TYPE), stat=allocok)
         IF (allocok > 0) THEN
            IF (ICNTL1>0) THEN
               WRITE(ICNTL1,*) 'PB allocation in ',
     &              'DMUMPS_INIT_OOC_BUF_PANEL'
            ENDIF
            IERR=-1
            I1=-13
            I2=OOC_NB_FILE_TYPE
            RETURN
         ENDIF
         AddVirtLibre(1:OOC_NB_FILE_TYPE)=0_8
         IF(allocated(NextAddVirtBuffer))THEN
            DEALLOCATE(NextAddVirtBuffer)
         ENDIF
         ALLOCATE(NextAddVirtBuffer(OOC_NB_FILE_TYPE), stat=allocok)
         IF (allocok > 0) THEN
            IF (ICNTL1>0) THEN
               WRITE(ICNTL1,*) 'PB allocation in ',
     &              'DMUMPS_INIT_OOC_BUF_PANEL'
            ENDIF
            IERR=-1
            I1=-13
            I2=OOC_NB_FILE_TYPE
            RETURN
         ENDIF
         NextAddVirtBuffer (1:OOC_NB_FILE_TYPE)  = BufferEmpty      
         IF(allocated(FIRST_VADDR_IN_BUF))THEN
            DEALLOCATE(FIRST_VADDR_IN_BUF)
         ENDIF
         ALLOCATE(FIRST_VADDR_IN_BUF(OOC_NB_FILE_TYPE), stat=allocok)
         IF (allocok > 0) THEN
            IF (ICNTL1>0) THEN
               WRITE(ICNTL1,*) 'PB allocation in ',
     &              'DMUMPS_INIT_OOC_BUF_PANEL'
            ENDIF
            IERR=-1
            I1=-13
            I2=OOC_NB_FILE_TYPE
            RETURN
         ENDIF
         CALL DMUMPS_OOC_INIT_DB_BUFFER_PANEL()   
      ELSE
         CALL DMUMPS_OOC_INIT_DB_BUFFER()
      ENDIF
      KEEP_OOC(223)=int(HBUF_SIZE)
      RETURN
      END SUBROUTINE DMUMPS_INIT_OOC_BUF
      SUBROUTINE DMUMPS_END_OOC_BUF()
      IMPLICIT NONE
      IF(allocated(BUF_IO))THEN
         DEALLOCATE(BUF_IO)
      ENDIF
      IF(allocated(I_SHIFT_FIRST_HBUF))THEN
         DEALLOCATE(I_SHIFT_FIRST_HBUF)
      ENDIF
      IF(allocated(I_SHIFT_SECOND_HBUF))THEN
         DEALLOCATE(I_SHIFT_SECOND_HBUF)
      ENDIF
      IF(allocated(I_SHIFT_CUR_HBUF))THEN
         DEALLOCATE(I_SHIFT_CUR_HBUF)
      ENDIF
      IF(allocated(I_REL_POS_CUR_HBUF))THEN
         DEALLOCATE(I_REL_POS_CUR_HBUF)
      ENDIF
      IF(allocated(LAST_IOREQUEST))THEN
         DEALLOCATE(LAST_IOREQUEST)
      ENDIF
      IF(allocated(CUR_HBUF))THEN
         DEALLOCATE(CUR_HBUF)
      ENDIF
      IF(PANEL_FLAG)THEN
         IF(allocated(NextAddVirtBuffer))THEN
            DEALLOCATE(NextAddVirtBuffer)
         ENDIF         
         IF(allocated(AddVirtLibre))THEN
            DEALLOCATE(AddVirtLibre)
         ENDIF
         IF(allocated(FIRST_VADDR_IN_BUF))THEN
            DEALLOCATE(FIRST_VADDR_IN_BUF)
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_END_OOC_BUF
      SUBROUTINE DMUMPS_OOC_INIT_DB_BUFFER()
      IMPLICIT NONE
      OOC_FCT_TYPE_LOC=1
      HBUF_SIZE = DIM_BUF_IO / int(2,kind=kind(DIM_BUF_IO))
      EARLIEST_WRITE_MIN_SIZE = 0
      I_SHIFT_FIRST_HBUF(OOC_FCT_TYPE_LOC) = 0_8
      I_SHIFT_SECOND_HBUF(OOC_FCT_TYPE_LOC) = HBUF_SIZE
      LAST_IOREQUEST(OOC_FCT_TYPE_LOC) = -1
      I_CUR_HBUF_NEXTPOS = 1
      I_CUR_HBUF_FSTPOS = 1
      I_SUB_HBUF_FSTPOS = 1
      CUR_HBUF(OOC_FCT_TYPE_LOC) = SECOND_HBUF
      CALL DMUMPS_OOC_NEXT_HBUF(OOC_FCT_TYPE_LOC)
      END SUBROUTINE DMUMPS_OOC_INIT_DB_BUFFER
      SUBROUTINE DMUMPS_OOC_COPY_DATA_TO_BUFFER(BLOCK,SIZE_OF_BLOCK,
     &     IERR)
      IMPLICIT NONE
      INTEGER(8) :: SIZE_OF_BLOCK
      DOUBLE PRECISION BLOCK(SIZE_OF_BLOCK)
      INTEGER, intent(out) :: IERR
      INTEGER(8) :: I
      IERR=0
      IF (I_REL_POS_CUR_HBUF(OOC_FCT_TYPE_LOC) +
     &    SIZE_OF_BLOCK <= HBUF_SIZE + 1_8) THEN
      ELSE
        CALL DMUMPS_OOC_DO_IO_AND_CHBUF(OOC_FCT_TYPE_LOC,IERR)
        IF(IERR.LT.0)THEN
           RETURN
        ENDIF
      END IF
      DO I = 1_8, SIZE_OF_BLOCK
        BUF_IO(I_SHIFT_CUR_HBUF(OOC_FCT_TYPE_LOC) +
     &        I_REL_POS_CUR_HBUF(OOC_FCT_TYPE_LOC) + I - 1_8) =
     &    BLOCK(I)
      END DO
      I_REL_POS_CUR_HBUF(OOC_FCT_TYPE_LOC) =
     &     I_REL_POS_CUR_HBUF(OOC_FCT_TYPE_LOC) + SIZE_OF_BLOCK
      RETURN
      END SUBROUTINE DMUMPS_OOC_COPY_DATA_TO_BUFFER
      SUBROUTINE DMUMPS_OOC_INIT_DB_BUFFER_PANEL()
      IMPLICIT NONE
      INTEGER(8) :: DIM_BUF_IO_L_OR_U
      INTEGER TYPEF, TYPEF_LAST
      INTEGER NB_DOUBLE_BUFFERS
      TYPEF_LAST = OOC_NB_FILE_TYPE
      NB_DOUBLE_BUFFERS = OOC_NB_FILE_TYPE
      DIM_BUF_IO_L_OR_U = DIM_BUF_IO /
     & int(NB_DOUBLE_BUFFERS,kind=kind(DIM_BUF_IO_L_OR_U))
      IF(.NOT.STRAT_IO_ASYNC)THEN
         HBUF_SIZE = DIM_BUF_IO_L_OR_U
      ELSE
         HBUF_SIZE = DIM_BUF_IO_L_OR_U / 2_8
      ENDIF
      DO TYPEF = 1, TYPEF_LAST
        LAST_IOREQUEST(TYPEF) = -1
        IF (TYPEF == 1 ) THEN
          I_SHIFT_FIRST_HBUF(TYPEF) = 0_8
        ELSE
          I_SHIFT_FIRST_HBUF(TYPEF) = DIM_BUF_IO_L_OR_U
        ENDIF
        IF(.NOT.STRAT_IO_ASYNC)THEN
           I_SHIFT_SECOND_HBUF(TYPEF) = I_SHIFT_FIRST_HBUF(TYPEF)
        ELSE
           I_SHIFT_SECOND_HBUF(TYPEF) = I_SHIFT_FIRST_HBUF(TYPEF) +
     &          HBUF_SIZE
        ENDIF
        CUR_HBUF(TYPEF) = SECOND_HBUF
        CALL DMUMPS_OOC_NEXT_HBUF(TYPEF)
      ENDDO
      I_CUR_HBUF_NEXTPOS = 1
      RETURN
      END SUBROUTINE DMUMPS_OOC_INIT_DB_BUFFER_PANEL
      SUBROUTINE DMUMPS_OOC_TRYIO_CHBUF_PANEL(TYPEF,IERR)
      IMPLICIT NONE
      INTEGER, INTENT(in)  :: TYPEF
      INTEGER, INTENT(out) :: IERR
      INTEGER IFLAG
      INTEGER NEW_IOREQUEST
      IERR=0
      CALL MUMPS_TEST_REQUEST_C(LAST_IOREQUEST(TYPEF),IFLAG,
     &     IERR)
      IF (IFLAG.EQ.1) THEN
         IERR = 0
         CALL DMUMPS_OOC_WRT_CUR_BUF2DISK(TYPEF,
     &        NEW_IOREQUEST,
     &        IERR)
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
         LAST_IOREQUEST(TYPEF) = NEW_IOREQUEST
         CALL DMUMPS_OOC_NEXT_HBUF(TYPEF)
         NextAddVirtBuffer(TYPEF)=BufferEmpty
         RETURN
      ELSE IF(IFLAG.LT.0)THEN
         WRITE(*,*)MYID_OOC,': ',ERR_STR_OOC(1:DIM_ERR_STR_OOC)
         RETURN
      ELSE
         IERR = 1
         RETURN
      ENDIF
      END SUBROUTINE DMUMPS_OOC_TRYIO_CHBUF_PANEL
      SUBROUTINE DMUMPS_OOC_UPD_VADDR_CUR_BUF (TYPEF,VADDR)
      IMPLICIT NONE
      INTEGER(8), INTENT(in) :: VADDR
      INTEGER, INTENT(in) :: TYPEF
      IF(I_REL_POS_CUR_HBUF(TYPEF).EQ.1_8)THEN
         FIRST_VADDR_IN_BUF(TYPEF)=VADDR
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_OOC_UPD_VADDR_CUR_BUF      
      SUBROUTINE DMUMPS_COPY_LU_TO_BUFFER( STRAT, TYPEF, MonBloc,
     &     AFAC, LAFAC,
     &     AddVirtCour, IPIVBEG, IPIVEND, LPANELeff,
     &     IERR)
      IMPLICIT NONE
      INTEGER,          INTENT(IN) :: TYPEF, IPIVBEG, IPIVEND, STRAT
      INTEGER(8),       INTENT(IN) :: LAFAC
      DOUBLE PRECISION, INTENT(IN) :: AFAC(LAFAC)
      INTEGER(8),        INTENT(IN) :: AddVirtCour
      TYPE(IO_BLOCK),   INTENT(IN) :: MonBloc   
      INTEGER,          INTENT(OUT):: LPANELeff
      INTEGER,          INTENT(OUT):: IERR
      INTEGER   :: II, NBPIVeff
      INTEGER(8) :: IPOS, IDIAG, IDEST
      INTEGER(8) :: DeltaIPOS
      INTEGER   :: StrideIPOS
      IERR=0
      IF (STRAT.NE.STRAT_WRITE_MAX.AND.STRAT.NE.STRAT_TRY_WRITE) THEN
         write(6,*) ' DMUMPS_COPY_LU_TO_BUFFER: STRAT Not implemented '
         CALL MUMPS_ABORT()
      ENDIF
      NBPIVeff = IPIVEND - IPIVBEG + 1
      IF (MonBloc%MASTER .AND. MonBloc%Typenode .NE. 3) THEN
        IF (TYPEF.EQ.TYPEF_L) THEN
          LPANELeff = (MonBloc%NROW-IPIVBEG+1)*NBPIVeff
        ELSE
          LPANELeff = (MonBloc%NCOL-IPIVBEG+1)*NBPIVeff
        ENDIF
      ELSE 
          LPANELeff = MonBloc%NROW*NBPIVeff
      ENDIF
      IF ( ( I_REL_POS_CUR_HBUF(TYPEF) + int(LPANELeff - 1,8)
     &     >
     &     HBUF_SIZE )
     &     .OR.
     &     ( (AddVirtCour.NE.NextAddVirtBuffer(TYPEF)) .AND.
     &     (NextAddVirtBuffer(TYPEF).NE.BufferEmpty) )
     &     ) THEN
         IF (STRAT.EQ.STRAT_WRITE_MAX) THEN
            CALL DMUMPS_OOC_DO_IO_AND_CHBUF(TYPEF,IERR) 
         ELSE IF (STRAT.EQ.STRAT_TRY_WRITE) THEN
            CALL DMUMPS_OOC_TRYIO_CHBUF_PANEL(TYPEF,IERR) 
            IF (IERR.EQ.1) RETURN
         ELSE
            write(6,*) 'DMUMPS_COPY_LU_TO_BUFFER: STRAT Not implemented'
         ENDIF
      ENDIF
      IF (IERR < 0 ) THEN
        RETURN
      ENDIF
      IF (NextAddVirtBuffer(TYPEF).EQ. BufferEmpty) THEN
         CALL DMUMPS_OOC_UPD_VADDR_CUR_BUF (TYPEF,AddVirtCour)
         NextAddVirtBuffer(TYPEF) = AddVirtCour
      ENDIF
      IF (MonBloc%MASTER .AND. MonBloc%Typenode .NE. 3) THEN
         IDIAG =  int(IPIVBEG-1,8)*int(MonBloc%NCOL,8) + int(IPIVBEG,8)
         IPOS   = IDIAG
         IDEST = I_SHIFT_CUR_HBUF(TYPEF) +
     &        I_REL_POS_CUR_HBUF(TYPEF)
         IF (TYPEF.EQ.TYPEF_L) THEN
            DO II = IPIVBEG, IPIVEND
               CALL dcopy(MonBloc%NROW-IPIVBEG+1, 
     &              AFAC(IPOS), MonBloc%NCOL,
     &              BUF_IO(IDEST), 1)
               IDEST = IDEST + int(MonBloc%NROW-IPIVBEG+1,8)
               IPOS  = IPOS  + 1_8
            ENDDO
         ELSE
            DO II = IPIVBEG, IPIVEND
               CALL dcopy(MonBloc%NCOL-IPIVBEG+1, 
     &              AFAC(IPOS), 1,
     &              BUF_IO(IDEST), 1)
               IDEST = IDEST + int(MonBloc%NCOL-IPIVBEG+1,8)
               IPOS  = IPOS  + int(MonBloc%NCOL,8)
            ENDDO
         ENDIF
      ELSE
         IDEST = I_SHIFT_CUR_HBUF(TYPEF) +
     &        I_REL_POS_CUR_HBUF(TYPEF)
         IF (MonBloc%Typenode.EQ.3) THEN
           DeltaIPOS  = int(MonBloc%NROW,8)
           StrideIPOS = 1
         ELSE
           DeltaIPOS  = 1_8
           StrideIPOS = MonBloc%NCOL
         ENDIF
         IPOS  = 1_8 + int(IPIVBEG - 1,8) * DeltaIPOS
         DO II = IPIVBEG, IPIVEND
            CALL dcopy(MonBloc%NROW, 
     &           AFAC(IPOS), StrideIPOS,
     &           BUF_IO(IDEST), 1)
            IDEST = IDEST+int(MonBloc%NROW,8)
            IPOS  = IPOS + DeltaIPOS
         ENDDO
      ENDIF
      I_REL_POS_CUR_HBUF(TYPEF) =
     &     I_REL_POS_CUR_HBUF(TYPEF) + int(LPANELeff,8)
      NextAddVirtBuffer(TYPEF) = NextAddVirtBuffer(TYPEF)
     &                         + int(LPANELeff,8)
      RETURN
      END SUBROUTINE DMUMPS_COPY_LU_TO_BUFFER
      END MODULE DMUMPS_OOC_BUFFER
