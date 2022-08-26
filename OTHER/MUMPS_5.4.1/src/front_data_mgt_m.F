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
      MODULE MUMPS_FRONT_DATA_MGT_M
      IMPLICIT NONE
      PRIVATE
C     --------------------------------------------
C     This module contains routines to manage
C     handlers of various data associated to
C     active fronts *during the factorization*.
C
C     It should be initialized at the beginning
C     of the factorization and terminated at the
C     end of the factorization.
C
C     There are two types of data, see below.
C
C     'A' is for active type 2 fronts: list must
C         be empty at the end of the factorization
C
C     'F' will be for general fronts -- currently used
C         for BLR fronts, in three situations:
C         1/ factorization of type 2 symmetric active fronts 
C            (requires temporary storage of BLR panels)      
C         2/ LRSOLVE: BLR factors are kept until solution phase
C            (liberated in JOB=-2 or at the beginning of a new facto)      
C         3/ LRCB: CB is dynamically allocated and compressed
C            (liberated before the end of the factorization)      
C
C     Only handlers are managed in this module.
C     The data itself is in the module above using it.
C     For example, FAC_MAPROW_DATA_M manages MAPROW
C     messages that arrive too early. It handles an
C     array that contains all early MAPROW messages
C     and that is indexed with the handlers managed
C     by MUMPS_FRONT_DATA_MGT_M.
C
C     --------------------------------------------
C
C     ===============
C     Public routines
C     ===============
      PUBLIC :: MUMPS_FDM_INIT,
     &          MUMPS_FDM_END,
     &          MUMPS_FDM_START_IDX,
     &          MUMPS_FDM_END_IDX
     &          , MUMPS_FDM_MOD_TO_STRUC
     &          , MUMPS_FDM_STRUC_TO_MOD
     &          , MUMPS_SAVE_RESTORE_FRONT_DATA
C     STACK_FREE_IDX(1:NB_FREE_IDX) holds the NB_FREE_IDX indices
C                                   of free handlers
C     STACK_FREE_IDX(NB_FREE_IDX+1:size(STACK_FREE_IDX)) is trash data
      TYPE FDM_STRUC_T
        INTEGER :: NB_FREE_IDX
        INTEGER, DIMENSION(:), POINTER :: STACK_FREE_IDX => null()
        INTEGER, DIMENSION(:), POINTER :: COUNT_ACCESS   => null()
      END TYPE FDM_STRUC_T
      TYPE (FDM_STRUC_T), TARGET, SAVE :: FDM_A, FDM_F
      CONTAINS
C
      SUBROUTINE MUMPS_FDM_INIT(WHAT, INITIAL_SIZE, INFO)
C
C     Purpose:
C     =======
C
C     Initialize handler data ('A' or 'F')
C
C     Arguments:
C     =========
C
      INTEGER, INTENT(IN) :: INITIAL_SIZE
      CHARACTER, INTENT(IN) :: WHAT  ! 'A' or 'F'
      INTEGER, INTENT(INOUT) :: INFO(2)
C
C     Local variables:
C     ===============
C
      INTEGER :: IERR
      TYPE (FDM_STRUC_T), POINTER :: FDM_PTR
C
      CALL MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
      ALLOCATE( FDM_PTR%STACK_FREE_IDX(INITIAL_SIZE),
     &          FDM_PTR%COUNT_ACCESS  (INITIAL_SIZE), stat=IERR )
      IF (IERR < 0) THEN
        INFO(1) = -13
        INFO(2) = INITIAL_SIZE * 2
        RETURN
      ENDIF
      CALL MUMPS_FDM_SET_ALL_FREE(FDM_PTR)
      RETURN
      END SUBROUTINE MUMPS_FDM_INIT
C
      SUBROUTINE MUMPS_FDM_END(WHAT)
C
C     Purpose:
C     =======
C     Free module datastructures associated to "WHAT" at
C     the end of a phase (typically factorization).
C
      CHARACTER, INTENT(IN) :: WHAT
C
C     Local variables
C     ===============
C
      TYPE (FDM_STRUC_T), POINTER :: FDM_PTR
C
      CALL MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
      IF (associated(FDM_PTR%STACK_FREE_IDX)) THEN
          DEALLOCATE(FDM_PTR%STACK_FREE_IDX)
          NULLIFY(FDM_PTR%STACK_FREE_IDX)
          FDM_PTR%NB_FREE_IDX=0
      ELSE
C         Should not be called twice or when array is unassociated
          WRITE(*,*) "Internal error 1 in MUMPS_FDM_END", WHAT
          CALL MUMPS_ABORT()
      ENDIF
      IF (associated(FDM_PTR%COUNT_ACCESS)) THEN
          DEALLOCATE(FDM_PTR%COUNT_ACCESS)
          NULLIFY(FDM_PTR%COUNT_ACCESS)
      ELSE
C     Should not be called twice or when array is unassociated
          WRITE(*,*) "Internal error 2 in MUMPS_FDM_END", WHAT
          CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_FDM_END
C
      SUBROUTINE MUMPS_FDM_MOD_TO_STRUC(WHAT, id_FDM_ENCODING,INFO)
C
C     Purpose:
C     =======
C
C     Save module information in struture.
C     id_FDM_ENCODING corresponds to id%FDM_F_ENCODING
C     This version requires that WHAT is equal to 'F'.
C
C     id_FDM_ENDODING takes responsibility of pointing to module
C     FDM_F information. This typically allows data from the module
C     to be passed from factorization to solve through the instance
C     and manage multiple instances.
C
      CHARACTER, INTENT(IN) :: WHAT
      INTEGER, INTENT(INOUT) :: INFO(2)
#if defined(MUMPS_F2003)
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                           id_FDM_ENCODING
#else
      CHARACTER, DIMENSION(:), POINTER :: id_FDM_ENCODING
#endif
C
C     Local variables
C     ===============
C
C     Character array of arbitrary dimension 1     
      CHARACTER :: CHAR_ARRAY(1)
      INTEGER :: CHAR_LENGTH, IERR
C
      IF (WHAT .NE. 'F') THEN
        WRITE(*,*) "Internal error 1 in MUMPS_FDM_MOD_TO_STRUC"
        CALL MUMPS_ABORT()
      ENDIF
      IF (associated(id_FDM_ENCODING)) THEN
C     Should be unassociated for this to work
        WRITE(*,*) "Internal error 2 in MUMPS_FDM_MOD_TO_STRUC"
        CALL MUMPS_ABORT()
      ENDIF
      CHAR_LENGTH=size(transfer(FDM_F,CHAR_ARRAY))
      ALLOCATE(id_FDM_ENCODING(CHAR_LENGTH), stat=IERR )
      IF (IERR < 0) THEN
        INFO(1) = -13
        INFO(2) = CHAR_LENGTH
        RETURN
      ENDIF
C     ------------------------------
C     Fill contents of pointer array
C     with FDM_F derived datatype
C     ------------------------------
      id_FDM_ENCODING = transfer(FDM_F,CHAR_ARRAY)
C     ----------------------------------------------
C     FDM_F is not to be used again before a call to
C     MUMPS_FDM_STRUC_TO_MOD, invalidate its content
C     ----------------------------------------------
      FDM_F%NB_FREE_IDX=-9999999
      NULLIFY(FDM_F%STACK_FREE_IDX)
      NULLIFY(FDM_F%COUNT_ACCESS)
      RETURN
      END SUBROUTINE MUMPS_FDM_MOD_TO_STRUC
C
      SUBROUTINE MUMPS_FDM_STRUC_TO_MOD(WHAT, id_FDM_ENCODING)
C
C     Purpose:
C     =======
C
C     Set module pointer information from id_FDM_ENCODING) typically
C     at beginning of solve. Suppress from structure since
C     responsibility of pointing to module data is now inside
C     the module.
C
      CHARACTER, INTENT(IN) :: WHAT
#if defined(MUMPS_F2003)
      CHARACTER, DIMENSION(:), POINTER, INTENT(INOUT)
     &                                        :: id_FDM_ENCODING
#else
      CHARACTER, DIMENSION(:), POINTER :: id_FDM_ENCODING
#endif
C
C     Local variables
C     ===============
C
      IF (.NOT.associated(id_FDM_ENCODING)) THEN
        WRITE(*,*) "Internal error 1 in MUMPS_FDM_STRUC_TO_MOD"
      ENDIF
      FDM_F=transfer(id_FDM_ENCODING,FDM_F)
C     Module is now responsible for accessing data.
      DEALLOCATE(id_FDM_ENCODING)
      NULLIFY(id_FDM_ENCODING)
      RETURN
      END SUBROUTINE MUMPS_FDM_STRUC_TO_MOD
C
      SUBROUTINE MUMPS_FDM_START_IDX(WHAT, FROM, IWHANDLER, INFO)
C
C     Purpose:
C     =======
C
C     Return a new free index/handler
C     (typically stored in IW)
C
      CHARACTER, INTENT(IN)  :: WHAT
      CHARACTER(LEN=*), INTENT(IN)  :: FROM !For debugging purposes only
      INTEGER, INTENT(INOUT) :: IWHANDLER
      INTEGER, INTENT(INOUT) :: INFO(2)
C
C     Local variables
C     ===============
C
      INTEGER :: OLD_SIZE, NEW_SIZE, IERR
      INTEGER :: I
      INTEGER, DIMENSION(:), POINTER :: TMP_COUNT_ACCESS
      TYPE(FDM_STRUC_T), POINTER :: FDM_PTR
      CALL MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
C
      IF (IWHANDLER .GT. 0) THEN
C       Already started, counter should at least be 1
         IF (FDM_PTR%COUNT_ACCESS(IWHANDLER) .LT. 1) THEN
          WRITE(*,*) "Internal error 1 in MUMPS_FDM_START_IDX",
     &    FDM_PTR%COUNT_ACCESS(IWHANDLER)
          CALL MUMPS_ABORT()
        ENDIF
        GOTO 100
      ENDIF
C
      IF (FDM_PTR%NB_FREE_IDX .EQ. 0) THEN
        OLD_SIZE = size(FDM_PTR%STACK_FREE_IDX)
        NEW_SIZE = (OLD_SIZE * 3) / 2 + 1 ! or something else
        FDM_PTR%NB_FREE_IDX = NEW_SIZE - OLD_SIZE
        DEALLOCATE(FDM_PTR%STACK_FREE_IDX)
        ALLOCATE(FDM_PTR%STACK_FREE_IDX(NEW_SIZE),
     &           TMP_COUNT_ACCESS(NEW_SIZE), stat=IERR)
        IF (IERR < 0) THEN
          INFO(1) = -13
          INFO(2) = NEW_SIZE
          RETURN
        ENDIF
C       All new handlers indices are created 
        DO I=1, FDM_PTR%NB_FREE_IDX
          FDM_PTR%STACK_FREE_IDX(I)=NEW_SIZE-I+1
        ENDDO
C       Count access: copy old ones
        DO I=1, OLD_SIZE
          TMP_COUNT_ACCESS(I)=FDM_PTR%COUNT_ACCESS(I)
        ENDDO
        DO I=OLD_SIZE+1, NEW_SIZE
          TMP_COUNT_ACCESS(I)=0
        ENDDO
        DEALLOCATE(FDM_PTR%COUNT_ACCESS)
        FDM_PTR%COUNT_ACCESS=>TMP_COUNT_ACCESS
      ENDIF
C
      IWHANDLER = FDM_PTR%STACK_FREE_IDX(FDM_PTR%NB_FREE_IDX)
      FDM_PTR%NB_FREE_IDX = FDM_PTR%NB_FREE_IDX - 1
 100  CONTINUE
C     Number of modules accessing this handler
      FDM_PTR%COUNT_ACCESS(IWHANDLER)=FDM_PTR%COUNT_ACCESS(IWHANDLER)+1
      RETURN
      END SUBROUTINE MUMPS_FDM_START_IDX
C
      SUBROUTINE MUMPS_FDM_END_IDX(WHAT, FROM, IWHANDLER)
C
C     Purpose:
C     =======
C
C     Notify than an index/handler has been freed.
C     Mark it free for future reuse.
C
      CHARACTER, INTENT(IN) :: WHAT
      CHARACTER(LEN=*), INTENT(IN) :: FROM ! for debug purposes only
      INTEGER, INTENT(INOUT) :: IWHANDLER
      TYPE(FDM_STRUC_T), POINTER :: FDM_PTR
C
      CALL MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
      IF (IWHANDLER .LE.0) THEN
C     Already ended
        WRITE(*,*) "Internal error 1 in MUMPS_FDM_END_IDX",IWHANDLER
        CALL MUMPS_ABORT()
      ENDIF
      FDM_PTR%COUNT_ACCESS(IWHANDLER)=FDM_PTR%COUNT_ACCESS(IWHANDLER)-1
      IF (FDM_PTR%COUNT_ACCESS(IWHANDLER) .LT. 0) THEN
C     Negative counter!
        WRITE(*,*) "Internal error 2 in MUMPS_FDM_END_IDX",
     &  IWHANDLER, FDM_PTR%COUNT_ACCESS(IWHANDLER)
        CALL MUMPS_ABORT()
      ENDIF
      IF (FDM_PTR%COUNT_ACCESS(IWHANDLER) .EQ.0 ) THEN
         IF (FDM_PTR%NB_FREE_IDX .GE. size(FDM_PTR%STACK_FREE_IDX)) THEN
          WRITE(*,*) "Internal error 3 in MUMPS_FDM_END_IDX"
          CALL MUMPS_ABORT()
        ENDIF
        FDM_PTR%NB_FREE_IDX = FDM_PTR%NB_FREE_IDX + 1
C       Having incremented the nb of free handlers we
C       store the index (IWHANDLER) that has been
C       effectively released for future reuse.
        FDM_PTR%STACK_FREE_IDX(FDM_PTR%NB_FREE_IDX) = IWHANDLER
        IWHANDLER = -8888 ! has been used and is now free
      ENDIF
C
      RETURN
      END SUBROUTINE MUMPS_FDM_END_IDX
C     ===================
C     Private subroutines
C     ===================
      SUBROUTINE MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
      CHARACTER, INTENT(IN) :: WHAT
#if defined(MUMPS_F2003)
      TYPE(FDM_STRUC_T), POINTER, INTENT(OUT) :: FDM_PTR
#else
      TYPE(FDM_STRUC_T), POINTER :: FDM_PTR
#endif
C
      IF ( WHAT .EQ. 'A' ) THEN
        FDM_PTR => FDM_A
      ELSE IF ( WHAT .EQ. 'F' ) THEN
        FDM_PTR => FDM_F
      ELSE
C     Should be called with either A or F
        WRITE(*,*) "Internal error 1 in MUMPS_FDM_INIT"
        WRITE(*,*) "Allowed arguments for WHAT are A or F"
        CALL MUMPS_ABORT()
      ENDIF
      END SUBROUTINE MUMPS_FDM_SET_PTR
      SUBROUTINE MUMPS_FDM_SET_ALL_FREE(FDM_PTR)
C
C     Purpose:
C     =======
C     Initialize the stack of free elements for the first time
C
      TYPE(FDM_STRUC_T), POINTER :: FDM_PTR
      INTEGER :: I
      FDM_PTR%NB_FREE_IDX = size(FDM_PTR%STACK_FREE_IDX)
      DO I = 1, FDM_PTR%NB_FREE_IDX
        FDM_PTR%STACK_FREE_IDX(I)=FDM_PTR%NB_FREE_IDX-I+1
        FDM_PTR%COUNT_ACCESS  (I)=0
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_FDM_SET_ALL_FREE
C
!     ---------- MUMPS_SAVE_RESTORE_FRONT_DATA ----------------------- !
      SUBROUTINE MUMPS_SAVE_RESTORE_FRONT_DATA(id_FDM_F_ENCODING
     &                 ,unit,MYID,mode
     &                 ,SIZE_GEST,SIZE_VARIABLES
     &                 ,SIZE_INT, TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,INFO)
      IMPLICIT NONE
C  =======
C  Purpose
C  =======
C
C     This routine is designed to manage a FDM_STRUC_T structure (save, restore, compute memory)
C     
C  ==========
C  Parameters
C  ==========
C
C     FDM_STRUC           : TYPE (FDM_STRUC_T) : the main structure
C
C     unit                : The unit of the file to be written or read
C
C     mode                : the type of operation to be performed by the routine 
C                           memory_save  = compute the size of the save file and of the structure
C                           save         = save the instace
C                           restore      = restore the instace
C
C     TOTAL_FILE_SIZE     : size of the file to be written or read
C
C     TOTAL_STRUC_SIZE    : size of the structure to be saved or restored
C
C     SIZE_INT            : size of an integer       
C
C     INFO                : copies of of INFO(1) and INFO(2) to allow save/restore of failled instaces
C      
      CHARACTER, DIMENSION(:), POINTER :: id_FDM_F_ENCODING
      INTEGER,intent(IN)::unit,MYID
      CHARACTER(len=*),intent(IN) :: mode
      INTEGER,INTENT(OUT) :: SIZE_GEST
      INTEGER(8),intent(OUT) :: SIZE_VARIABLES
      INTEGER(8),intent(IN) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER,intent(INOUT):: INFO(2)
      INTEGER,intent(IN):: SIZE_INT
      INTEGER(8),intent(INOUT):: size_read,size_allocated,size_written
      INTEGER:: NbRecords,NbSubRecords
      INTEGER:: SIZE_GEST_FDM_F
      INTEGER(8):: SIZE_VARIABLES_FDM_F
      NbRecords=0
      SIZE_GEST_FDM_F=0
      SIZE_VARIABLES_FDM_F=0_8
      SIZE_GEST=0
      SIZE_VARIABLES=0_8
      if((trim(mode).EQ."memory_save").OR.(trim(mode).EQ."save")) then
         call MUMPS_FDM_STRUC_TO_MOD("F",id_FDM_F_ENCODING)
      endif
      if(trim(mode).EQ."memory_save") then
         CALL MUMPS_SAVE_RESTORE_FDM_STRUC(
     &        FDM_F
     &        ,unit,MYID,"memory_save"
     &        ,SIZE_GEST_FDM_F
     &        ,SIZE_VARIABLES_FDM_F
     &        ,SIZE_INT,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &        ,size_read,size_allocated,size_written
     &        ,INFO)
      elseif(trim(mode).EQ."save") then
         CALL MUMPS_SAVE_RESTORE_FDM_STRUC(
     &        FDM_F
     &        ,unit,MYID,"save"
     &        ,SIZE_GEST_FDM_F
     &        ,SIZE_VARIABLES_FDM_F
     &        ,SIZE_INT,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &        ,size_read,size_allocated,size_written
     &        ,INFO)
         IF ( INFO(1) .LT. 0 ) GOTO 100
      elseif(trim(mode).EQ."restore") then
         CALL MUMPS_SAVE_RESTORE_FDM_STRUC(
     &        FDM_F
     &        ,unit,MYID,"restore"
     &        ,SIZE_GEST_FDM_F
     &        ,SIZE_VARIABLES_FDM_F
     &        ,SIZE_INT, TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &        ,size_read,size_allocated,size_written
     &        ,INFO)
         IF ( INFO(1) .LT. 0 ) GOTO 100
      endif
      if(trim(mode).EQ."memory_save") then
C     If the size to write (SIZE_VARIABLES) is greater than 2^31
C     Subrecords are created which need to be taken into account in
C     the file size computation
         NbSubRecords=int(SIZE_VARIABLES/huge(0))
         IF(NbSubRecords.GT.0) then
            NbRecords=NbRecords+NbSubRecords
         ENDIF
      elseif(trim(mode).EQ."save") then
         size_written=size_written+SIZE_VARIABLES
     &        +int(SIZE_GEST,kind=8)
#if !defined(MUMPS_F2003)
         size_written=size_written
     &        +int(2*SIZE_INT*NbRecords,kind=8)
#endif
      elseif(trim(mode).EQ."restore") then
         size_allocated=size_allocated+SIZE_VARIABLES
         size_read=size_read+SIZE_VARIABLES
     &        +int(SIZE_GEST,kind=8)
#if !defined(MUMPS_F2003)
         size_read=size_read
     &        +int(2*SIZE_INT*NbRecords,kind=8)
#endif
      endif
      if(trim(mode).EQ."memory_save") then
         SIZE_VARIABLES=SIZE_VARIABLES+SIZE_VARIABLES_FDM_F
         SIZE_GEST=SIZE_GEST+SIZE_GEST_FDM_F
#if !defined(MUMPS_F2003)
C     If the file is not written with access="stream", which is only done in MUMPS_F2003,
C     the record length's is written at the beginning and at the end of each record
C     This is done using 2 INTEGERs so we use 2*SIZE_INT more space for each record
         SIZE_GEST=SIZE_GEST+2*SIZE_INT*NbRecords
#endif
      endif
      call MUMPS_FDM_MOD_TO_STRUC("F",id_FDM_F_ENCODING,INFO(1))
 100  continue
      RETURN
      END SUBROUTINE MUMPS_SAVE_RESTORE_FRONT_DATA
!     --------------------------------- MUMPS_SAVE_RESTORE_BLR_STRUC ----------------------------- !
      SUBROUTINE MUMPS_SAVE_RESTORE_FDM_STRUC(FDM_STRUC
     &     ,unit,MYID,mode
     &     ,Local_SIZE_GEST,Local_SIZE_VARIABLES
     &     ,SIZE_INT ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,size_read,size_allocated,size_written
     &     ,INFO)     
      IMPLICIT NONE
C  =======
C  Purpose
C  =======
C
C     This routine is designed to manage a BLR_STRUC_T structure (save, restore, compute memory)
C     
C  ==========
C  Parameters
C  ==========
C
C     BLR_STRUC           : TYPE (BLR_STRUC_T) : the main structure
C
C     unit                : The unit of the file to be written or read
C
C     mode                : the type of operation to be performed by the routine 
C                           memory_save  = compute the size of the save file and of the structure
C                           save         = save the instace
C                           restore      = restore the instace
C
C     TOTAL_FILE_SIZE     : size of the file to be written or read
C
C     TOTAL_STRUC_SIZE    : size of the structure to be saved or restored
C
C     SIZE_INT            : size of an integer      
C
C     INFO1/INFO2         : copies of of INFO(1) and INFO(2) to allow save/restore of failled instaces
C
      TYPE(FDM_STRUC_T) :: FDM_STRUC
      INTEGER,intent(IN)::unit,MYID
      CHARACTER(len=*),intent(IN) :: mode
      INTEGER,INTENT(OUT) :: Local_SIZE_GEST
      INTEGER(8),intent(OUT) :: Local_SIZE_VARIABLES
      INTEGER,intent(INOUT):: INFO(2)
      INTEGER,intent(IN):: SIZE_INT
      INTEGER(8),intent(IN) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER(8),intent(INOUT):: size_read,size_allocated,size_written
      INTEGER :: NBVARIABLES_FDM_STRUC_T
      PARAMETER (NBVARIABLES_FDM_STRUC_T = 3)
      CHARACTER(len=30), dimension(NBVARIABLES_FDM_STRUC_T)::
     &     VARIABLES_FDM_STRUC_T
      CHARACTER(len=30) :: TMP_STRING
      INTEGER(8),dimension(NBVARIABLES_FDM_STRUC_T)::
     &     SIZE_VARIABLES_FDM_STRUC_T
      INTEGER,dimension(NBVARIABLES_FDM_STRUC_T)::SIZE_GEST_FDM_STRUC_T
      INTEGER,dimension(NBVARIABLES_FDM_STRUC_T)::NbRecords_FDM_STRUC_T
      INTEGER:: size_array1,dummy,allocok
      INTEGER:: err,i1,NbSubRecords,Local_NbRecords
      VARIABLES_FDM_STRUC_T(1)="NB_FREE_IDX"
      VARIABLES_FDM_STRUC_T(2)="STACK_FREE_IDX"
      VARIABLES_FDM_STRUC_T(3)="COUNT_ACCESS"
      SIZE_VARIABLES_FDM_STRUC_T(:)=0_8
      SIZE_GEST_FDM_STRUC_T(:)=0
      NbRecords_FDM_STRUC_T(:)=0 
C     
C     BEGINNING OF THE MAIN LOOP ON ALL VARIABLES OF THE STRUCTURE
C
      DO i1=1,NBVARIABLES_FDM_STRUC_T 
         TMP_STRING = VARIABLES_FDM_STRUC_T(i1)
         SELECT CASE(TMP_STRING)
         CASE("NB_FREE_IDX")
            NbRecords_FDM_STRUC_T(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES_FDM_STRUC_T(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               SIZE_VARIABLES_FDM_STRUC_T(i1)=SIZE_INT
               write(unit,iostat=err) FDM_STRUC%NB_FREE_IDX
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES_FDM_STRUC_T(i1)=SIZE_INT
               read(unit,iostat=err) FDM_STRUC%NB_FREE_IDX
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif 
               IF ( INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("STACK_FREE_IDX")
            NbRecords_FDM_STRUC_T(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(FDM_STRUC%STACK_FREE_IDX)) THEN
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=
     &                 size(FDM_STRUC%STACK_FREE_IDX,1)*SIZE_INT
               ELSE
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(FDM_STRUC%STACK_FREE_IDX)) THEN
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=
     &                 size(FDM_STRUC%STACK_FREE_IDX,1)*SIZE_INT
                  write(unit,iostat=err)
     &                 size(FDM_STRUC%STACK_FREE_IDX,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) FDM_STRUC%STACK_FREE_IDX
               ELSE
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=0_8
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(FDM_STRUC%STACK_FREE_IDX)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=size_array1*SIZE_INT
                  allocate(FDM_STRUC%STACK_FREE_IDX(size_array1),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  read(unit,iostat=err) FDM_STRUC%STACK_FREE_IDX
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("COUNT_ACCESS")
            NbRecords_FDM_STRUC_T(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(FDM_STRUC%COUNT_ACCESS)) THEN
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=
     &                 size(FDM_STRUC%COUNT_ACCESS,1)*SIZE_INT
               ELSE
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(FDM_STRUC%COUNT_ACCESS)) THEN
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=
     &                 size(FDM_STRUC%COUNT_ACCESS,1)*SIZE_INT
                  write(unit,iostat=err)
     &                 size(FDM_STRUC%COUNT_ACCESS,1)
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) FDM_STRUC%COUNT_ACCESS
               ELSE
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=0_8
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    INFO(2))
                  endif
                  IF ( INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(FDM_STRUC%COUNT_ACCESS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT*2
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST_FDM_STRUC_T(i1)=SIZE_INT
                  SIZE_VARIABLES_FDM_STRUC_T(i1)=size_array1*SIZE_INT
                  allocate(FDM_STRUC%COUNT_ACCESS(size_array1),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,INFO(2))
                  endif
                  read(unit,iostat=err) FDM_STRUC%COUNT_ACCESS
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,INFO(2))
               endif
               IF ( INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE DEFAULT
         END SELECT
         if(trim(mode).EQ."memory_save") then
C     If the size to write (SIZE_VARIABLES_FDM_STRUC_T(i1)) is greater than 2^31
C     Subrecords are created which need to be taken into account in
C     the file size computation
            NbSubRecords=int(SIZE_VARIABLES_FDM_STRUC_T(i1)/huge(0))
            IF(NbSubRecords.GT.0) then
               NbRecords_FDM_STRUC_T(i1)=NbRecords_FDM_STRUC_T(i1)
     &              +NbSubRecords
            ENDIF
         elseif(trim(mode).EQ."save") then
            size_written=size_written+SIZE_VARIABLES_FDM_STRUC_T(i1)
     &           +int(SIZE_GEST_FDM_STRUC_T(i1),kind=8)
#if !defined(MUMPS_F2003)
            size_written=size_written
     &           +int(2*SIZE_INT*NbRecords_FDM_STRUC_T(i1),kind=8)
#endif
         elseif(trim(mode).EQ."restore") then
            size_allocated=size_allocated+
     &           SIZE_VARIABLES_FDM_STRUC_T(i1)
            size_read=size_read+SIZE_VARIABLES_FDM_STRUC_T(i1)
     &           +int(SIZE_GEST_FDM_STRUC_T(i1),kind=8)
#if !defined(MUMPS_F2003)
            size_read=size_read
     &           +int(2*SIZE_INT*NbRecords_FDM_STRUC_T(i1),kind=8)
#endif
         endif
      ENDDO
      if(trim(mode).EQ."memory_save") then
         Local_SIZE_VARIABLES=sum(SIZE_VARIABLES_FDM_STRUC_T)
         Local_SIZE_GEST=sum(SIZE_GEST_FDM_STRUC_T)        
#if !defined(MUMPS_F2003)
         Local_NbRecords=sum(NbRecords_FDM_STRUC_T)
         Local_SIZE_GEST=Local_SIZE_GEST+2*SIZE_INT*Local_NbRecords
#endif
      endif
 100  continue
      RETURN
      END SUBROUTINE MUMPS_SAVE_RESTORE_FDM_STRUC
      END MODULE MUMPS_FRONT_DATA_MGT_M
