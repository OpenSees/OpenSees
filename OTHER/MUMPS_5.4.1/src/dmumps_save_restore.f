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
      MODULE DMUMPS_SAVE_RESTORE
      USE DMUMPS_STRUC_DEF
      USE DMUMPS_SAVE_RESTORE_FILES
      USE DMUMPS_LR_DATA_M
      USE MUMPS_FRONT_DATA_MGT_M
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE DMUMPS_REMOVE_SAVED(id)
      USE DMUMPS_OOC
      INCLUDE 'mpif.h'
      INTEGER MASTER
      PARAMETER ( MASTER = 0 )
      TYPE (DMUMPS_STRUC) :: id
      CHARACTER(len=LEN_SAVE_FILE) :: RESTOREFILE, INFOFILE
      INTEGER :: fileunit, ierr, SIZE_INT, SIZE_INT8
      INTEGER(8) :: size_read, TOTAL_FILE_SIZE, TOTAL_STRUC_SIZE
      INTEGER :: READ_OOC_FILE_NAME_LENGTH,READ_SYM,READ_PAR,READ_NPROCS
      CHARACTER(len=LEN_SAVE_FILE) :: READ_OOC_FIRST_FILE_NAME
      CHARACTER :: READ_ARITH
      LOGICAL :: READ_INT_TYPE_64
      CHARACTER(len=23) :: READ_HASH
      LOGICAL :: FORTRAN_VERSION_OK,UNIT_OK,UNIT_OP
      LOGICAL :: SAME_OOC
      INTEGER :: ICNTL34, MAX_LENGTH, FLAG_SAME, SUM_FLAG_SAME
      TYPE (DMUMPS_STRUC) :: localid
      ierr = 0
      call DMUMPS_GET_SAVE_FILES(id,RESTOREFILE,INFOFILE)
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      fileunit = 40
      inquire (UNIT=fileunit,exist=UNIT_OK,opened=UNIT_OP)
      IF(.NOT.UNIT_OK .OR. UNIT_OP) THEN
         id%INFO(1) = -79
         id%INFO(2) = fileunit
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      open(UNIT=fileunit,FILE=RESTOREFILE
#if defined(MUMPS_F2003)
     &     ,ACCESS="stream"
#endif
     &     ,STATUS='old',FORM='unformatted',IOSTAT=ierr)
      IF ( ierr.ne.0 ) THEN
        id%INFO(1) = -74
        id%INFO(2) = 0
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      SIZE_INT = id%KEEP(34)
      SIZE_INT8 = id%KEEP(34)*id%KEEP(10)
      size_read = 0
      call MUMPS_READ_HEADER(fileunit,ierr,size_read,SIZE_INT,
     &     SIZE_INT8, TOTAL_FILE_SIZE, TOTAL_STRUC_SIZE,
     &     READ_ARITH, READ_INT_TYPE_64,
     &     READ_OOC_FILE_NAME_LENGTH, READ_OOC_FIRST_FILE_NAME,
     &     READ_HASH,READ_SYM,READ_PAR,READ_NPROCS,
     &     FORTRAN_VERSION_OK)
      close(fileunit)
      if (ierr.ne.0) THEN
         id%INFO(1) = -75
         CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &              ,id%INFO(2))
      elseif (.NOT.FORTRAN_VERSION_OK) THEN
         id%INFO(1) = -73
         id%INFO(2) = 1
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL DMUMPS_CHECK_HEADER(id,.TRUE.,READ_INT_TYPE_64,
     &           READ_HASH, READ_NPROCS,
     &           READ_ARITH, READ_SYM, READ_PAR)
      IF ( id%INFO(1) .LT. 0 ) RETURN
      ICNTL34 = -99998
      IF (id%MYID.EQ.MASTER) THEN
        ICNTL34 = id%ICNTL(34)
      ENDIF
      CALL MPI_BCAST( ICNTL34, 1, MPI_INTEGER, MASTER, id%COMM, ierr )
      CALL DMUMPS_CHECK_FILE_NAME(id, READ_OOC_FILE_NAME_LENGTH,
     &           READ_OOC_FIRST_FILE_NAME, SAME_OOC)
      CALL MPI_ALLREDUCE(READ_OOC_FILE_NAME_LENGTH,MAX_LENGTH,1,
     &                   MPI_INTEGER,MPI_MAX,id%COMM,ierr)
      IF (MAX_LENGTH.NE.-999) THEN
        FLAG_SAME = 0
        IF (SAME_OOC) THEN
          FLAG_SAME = 1
        ENDIF
        CALL MPI_ALLREDUCE(FLAG_SAME,SUM_FLAG_SAME,1,
     &                   MPI_INTEGER,MPI_SUM,id%COMM,ierr)
        IF (SUM_FLAG_SAME.NE.0) THEN
          IF (ICNTL34 .EQ. 1) THEN
            id%ASSOCIATED_OOC_FILES = .TRUE.
          ELSE
            id%ASSOCIATED_OOC_FILES = .FALSE.
          ENDIF
        ELSE
          IF (ICNTL34 .NE. 1) THEN
            localid%COMM = id%COMM
            localid%INFO(1) = 0
            localid%MYID = id%MYID
            localid%NPROCS = id%NPROCS
            localid%KEEP(10) = id%KEEP(10)
            localid%SAVE_PREFIX = id%SAVE_PREFIX
            localid%SAVE_DIR = id%SAVE_DIR
            call DMUMPS_RESTORE_OOC(localid)
            IF ( localid%INFO(1) .EQ. 0 ) THEN
              localid%ASSOCIATED_OOC_FILES = .FALSE.
              IF (READ_OOC_FILE_NAME_LENGTH.NE.-999) THEN
                call DMUMPS_OOC_CLEAN_FILES(localid,ierr)
                IF ( ierr.ne.0 ) THEN
                  id%INFO(1) = -90
                  id%INFO(2) = id%MYID
                ENDIF
              ENDIF
            ENDIF
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) RETURN 
          ENDIF
        ENDIF
      ENDIF
      call MUMPS_CLEAN_SAVED_DATA(id%MYID,ierr,RESTOREFILE,INFOFILE)
      IF (ierr.ne.0) THEN
        id%INFO(1) = -76
        id%INFO(2) = id%MYID
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      END SUBROUTINE DMUMPS_REMOVE_SAVED
      SUBROUTINE DMUMPS_RESTORE_OOC(localid)
      INCLUDE 'mpif.h'
      INTEGER::ierr,IN,NBVARIABLES,NBVARIABLES_ROOT
      CHARACTER(len=LEN_SAVE_FILE):: restore_file_ooc,INFO_FILE
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES_ROOT
      INTEGER,allocatable, dimension(:)::SIZE_GEST
      INTEGER,allocatable, dimension(:)::SIZE_GEST_ROOT
      INTEGER:: INFO1,INFO2,INFOG1,INFOG2,allocok
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      LOGICAL :: UNIT_OK,UNIT_OP
      TYPE (DMUMPS_STRUC) :: localid
      NBVARIABLES=186
      NBVARIABLES_ROOT=35
      allocate(SIZE_VARIABLES(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       localid%INFO(1)  =-13
       localid%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_VARIABLES_ROOT(NBVARIABLES_ROOT),  stat=allocok)
      if (allocok .GT. 0) THEN
       localid%INFO(1)  =-13
       localid%INFO(2) = NBVARIABLES_ROOT
      endif
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       localid%INFO(1)  =-13
       localid%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST_ROOT(NBVARIABLES_ROOT),  stat=allocok)
      if (allocok .GT. 0) THEN
       localid%INFO(1)  =-13
       localid%INFO(2) = NBVARIABLES_ROOT
      endif
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      SIZE_VARIABLES(:)=0_8
      SIZE_VARIABLES_ROOT(:)=0_8
      SIZE_GEST(:)=0
      SIZE_GEST_ROOT(:)=0
      TOTAL_FILE_SIZE=0_8
      TOTAL_STRUC_SIZE=0_8
      INFO1 = -999
      INFO2 = -999
      INFOG1 = -999
      INFOG2 = -999
      CALL DMUMPS_GET_SAVE_FILES(localid,restore_file_ooc,INFO_FILE)
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      IN=50
      inquire(UNIT=IN,exist=UNIT_OK,opened=UNIT_OP)
      IF(.NOT.UNIT_OK .OR. UNIT_OP) THEN
         localid%INFO(1) = -79
         localid%INFO(2) = IN
      ENDIF
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      open(UNIT=IN,FILE=restore_file_ooc
#if defined(MUMPS_F2003)
     &     ,ACCESS="stream"
#endif
     &     ,STATUS='old',form='unformatted',iostat=ierr)
      if(ierr.ne.0) THEN
         localid%INFO(1) = -74
         localid%INFO(2) = 0
      endif  
      CALL MUMPS_PROPINFO( localid%ICNTL(1), localid%INFO(1),
     &     localid%COMM, localid%MYID )
      IF ( localid%INFO(1) .LT. 0 ) RETURN
      CALL DMUMPS_SAVE_RESTORE_STRUCTURE(localid,IN,"restore_ooc"
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOT,SIZE_VARIABLES_ROOT,SIZE_GEST_ROOT
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)     
      CLOSE(IN)
      deallocate(SIZE_VARIABLES,SIZE_VARIABLES_ROOT)
      deallocate(SIZE_GEST,SIZE_GEST_ROOT)
      RETURN
      END SUBROUTINE DMUMPS_RESTORE_OOC
      SUBROUTINE DMUMPS_COMPUTE_MEMORY_SAVE(id,
     &     TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE)
      INCLUDE 'mpif.h'
      INTEGER::NBVARIABLES,NBVARIABLES_ROOT
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES_ROOT
      INTEGER,allocatable, dimension(:)::SIZE_GEST
      INTEGER,allocatable, dimension(:)::SIZE_GEST_ROOT
      INTEGER :: INFO1,INFO2,INFOG1,INFOG2,allocok
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      TYPE (DMUMPS_STRUC) :: id
      NBVARIABLES=186
      NBVARIABLES_ROOT=35
      allocate(SIZE_VARIABLES(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_VARIABLES_ROOT(NBVARIABLES_ROOT),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOT
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST_ROOT(NBVARIABLES_ROOT),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOT
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      SIZE_VARIABLES(:)=0_8
      SIZE_VARIABLES_ROOT(:)=0_8
      SIZE_GEST(:)=0
      SIZE_GEST_ROOT(:)=0
      TOTAL_FILE_SIZE=0_8
      TOTAL_STRUC_SIZE=0_8
      INFO1 = -999
      INFO2 = -999
      INFOG1 = -999
      INFOG2 = -999
      CALL DMUMPS_SAVE_RESTORE_STRUCTURE(id,0,"memory_save"
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOT,SIZE_VARIABLES_ROOT,SIZE_GEST_ROOT
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      deallocate(SIZE_VARIABLES,SIZE_VARIABLES_ROOT)
      deallocate(SIZE_GEST,SIZE_GEST_ROOT)
      RETURN
      END SUBROUTINE DMUMPS_COMPUTE_MEMORY_SAVE
      SUBROUTINE DMUMPS_SAVE(id)
      INCLUDE 'mpif.h'
      INTEGER::ierr,OUT,NBVARIABLES,NBVARIABLES_ROOT,OUTINFO
      CHARACTER(len=LEN_SAVE_FILE):: SAVE_FILE,INFO_FILE
      LOGICAL:: SAVE_FILE_exist,INFO_FILE_exist
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES_ROOT
      INTEGER,allocatable, dimension(:)::SIZE_GEST
      INTEGER,allocatable, dimension(:)::SIZE_GEST_ROOT
      INTEGER:: INFO1,INFO2,INFOG1,INFOG2,MPG
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      LOGICAL :: PROKG,UNIT_OK,UNIT_OP
      INTEGER :: I,J,K,H,allocok
      CHARACTER(len=1) :: TMP_OOC_NAMES(350)
      TYPE (DMUMPS_STRUC) :: id
      INFO1 = id%INFO(1)
      INFO2 = id%INFO(2)
      INFOG1 = id%INFO(1)
      INFOG2 = id%INFO(1)
      id%INFO(1)=0
      id%INFO(2)=0
      id%INFOG(1)=0
      id%INFOG(2)=0
      MPG= id%ICNTL(3)
      PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. 0 )
      NBVARIABLES=186
      NBVARIABLES_ROOT=35
      allocate(SIZE_VARIABLES(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_VARIABLES_ROOT(NBVARIABLES_ROOT),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOT
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST_ROOT(NBVARIABLES_ROOT),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOT
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      SIZE_VARIABLES(:)=0_8
      SIZE_VARIABLES_ROOT(:)=0_8
      SIZE_GEST(:)=0
      SIZE_GEST_ROOT(:)=0
      TOTAL_FILE_SIZE=0_8
      TOTAL_STRUC_SIZE=0_8
      TMP_OOC_NAMES(:)="?"
      CALL DMUMPS_SAVE_RESTORE_STRUCTURE(id,0,"memory_save"
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOT,SIZE_VARIABLES_ROOT,SIZE_GEST_ROOT
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      CALL DMUMPS_GET_SAVE_FILES(id,SAVE_FILE,INFO_FILE)
      IF ( id%INFO(1) .LT. 0 ) RETURN
      inquire(FILE=SAVE_FILE, EXIST=SAVE_FILE_exist)
      IF(SAVE_FILE_exist) THEN
         id%INFO(1) = -70
         id%INFO(2) = 0
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      OUT=60
      inquire (UNIT=OUT,exist=UNIT_OK,opened=UNIT_OP)
      IF(.NOT.UNIT_OK .OR. UNIT_OP) THEN
         id%INFO(1) = -79
         id%INFO(2) = OUT
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      open(UNIT=OUT,FILE=SAVE_FILE
#if defined(MUMPS_F2003)
     &     ,ACCESS="stream"
#endif
     &     ,STATUS='new',form='unformatted',iostat=ierr)
      if(ierr.ne.0) THEN
         id%INFO(1) = -71
         id%INFO(2) = 0
      endif  
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      inquire(FILE=INFO_FILE, EXIST=INFO_FILE_exist)
      IF(INFO_FILE_exist) THEN
         id%INFO(1) = -70
         id%INFO(2) = 0
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      OUTINFO=70
      inquire (UNIT=OUTINFO,exist=UNIT_OK,opened=UNIT_OP)
      IF(.NOT.UNIT_OK .OR. UNIT_OP) THEN
         id%INFO(1) = -79
         id%INFO(2) = OUTINFO
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      open(UNIT=OUTINFO,FILE=INFO_FILE,STATUS='new',iostat=ierr)
      if(ierr.ne.0) THEN
         id%INFO(1) = -71
         id%INFO(2) = 0
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL DMUMPS_SAVE_RESTORE_STRUCTURE(id,OUT,"save"
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOT,SIZE_VARIABLES_ROOT,SIZE_GEST_ROOT
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      if(id%INFO(1).EQ.0) then
         id%INFO(1)=INFO1
         id%INFO(2)=INFO2
         id%INFOG(1)=INFOG1
         id%INFOG(2)=INFOG2
         CLOSE(OUT)
         if(id%INFO(1).NE.0) then
            write(MPG,*) "Warning: "
     &           ,"saved instance has negative INFO(1):"
     &           , id%INFO(1)
         endif
        IF(PROKG) THEN
           write(MPG,*) "Save done successfully"
           IF(id%KEEP(201).EQ.1) THEN
              K=1
              write(MPG,*) "The corresponding OOC files are:"
              DO I=1,id%OOC_NB_FILE_TYPE
                 DO J=1,id%OOC_NB_FILES(I)
                    DO H=1,id%OOC_FILE_NAME_LENGTH(K)-2
                       TMP_OOC_NAMES(H)=id%OOC_FILE_NAMES(K,H)
                    ENDDO
                    write(MPG,*)
     &                   TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(K)-2)
                    K=K+1
                 ENDDO
              ENDDO
           ENDIF
        ENDIF
        write(OUTINFO,*) "Save done by DMUMPS ",
     &       trim(adjustl(id%VERSION_NUMBER)),
     &       " after JOB=",id%KEEP(40)+456789,
     &       " With SYM, PAR =",id%KEEP(50),id%KEEP(46)
        write(OUTINFO,*) "On ",id%NPROCS," processes"
        if((id%ICNTL(18).EQ.0).AND.(id%ICNTL(5).EQ.0)) then
           write(OUTINFO,*) "with N, NNZ ", id%N, id%NNZ
        elseif((id%ICNTL(18).EQ.1).AND.(id%ICNTL(5).EQ. 0)) then
           write(OUTINFO,*) "with N, NNZ_loc=", id%N, id%NNZ_loc
        elseif((id%ICNTL(18).EQ.0).AND.(id%ICNTL(5).EQ. 1)) then
           write(OUTINFO,*) "with N, NELT=", id%N, id%NELT
        endif
        IF(id%KEEP(10).EQ.1) THEN
           write(OUTINFO,*) "With a default integer size of 64 bits"
        ELSE
           write(OUTINFO,*) "With a default integer size of 32 bits"
        ENDIF
#if defined(MUMPS_F2003)
        write(OUTINFO,*) "Using MUMPS_F2003"
#endif        
        write(OUTINFO,*) ''
        write(OUTINFO,*) "The corresponding save file is:"
        write(OUTINFO,*) trim(adjustl(SAVE_FILE))
        write(OUTINFO,*) "of size",TOTAL_FILE_SIZE, " Bytes"
        IF(id%KEEP(201).EQ.1) THEN
           write(OUTINFO,*) ''
           write(OUTINFO,*) "The corresponding OOC files are:"
           K=1
           DO I=1,id%OOC_NB_FILE_TYPE
              DO J=1,id%OOC_NB_FILES(I)
                 DO H=1,id%OOC_FILE_NAME_LENGTH(K)-2
                    TMP_OOC_NAMES(H)=id%OOC_FILE_NAMES(K,H)
                 ENDDO
                 write(OUTINFO,*)
     &                TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(K)-2)
                 K=K+1
              ENDDO
           ENDDO
        ENDIF
        CLOSE(OUTINFO)
      else
         CLOSE(OUT,STATUS='delete')
         CLOSE(OUTINFO,STATUS='delete')         
      endif
      deallocate(SIZE_VARIABLES,SIZE_VARIABLES_ROOT)
      deallocate(SIZE_GEST,SIZE_GEST_ROOT)
      if (id%KEEP(201) .GT. 0) THEN
        id%ASSOCIATED_OOC_FILES=.TRUE.
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SAVE
      SUBROUTINE DMUMPS_RESTORE(id)
      INCLUDE 'mpif.h'
      INTEGER::ierr,IN,NBVARIABLES,NBVARIABLES_ROOT
      CHARACTER(len=LEN_SAVE_FILE):: restore_file,INFO_FILE
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES
      INTEGER(8),allocatable, dimension(:)::SIZE_VARIABLES_ROOT
      INTEGER,allocatable, dimension(:)::SIZE_GEST
      INTEGER,allocatable, dimension(:)::SIZE_GEST_ROOT
      INTEGER:: INFO1,INFO2,INFOG1,INFOG2,MPG,MP,JOB
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      LOGICAL :: PROKG,UNIT_OK,UNIT_OP
      INTEGER :: I,J,K,H,allocok
      CHARACTER(len=1) :: TMP_OOC_NAMES(350)
      TYPE (DMUMPS_STRUC) :: id
      NBVARIABLES=186
      NBVARIABLES_ROOT=35
      allocate(SIZE_VARIABLES(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_VARIABLES_ROOT(NBVARIABLES_ROOT),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOT
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      allocate(SIZE_GEST_ROOT(NBVARIABLES_ROOT),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOT
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      SIZE_VARIABLES(:)=0_8
      SIZE_VARIABLES_ROOT(:)=0_8
      SIZE_GEST(:)=0
      SIZE_GEST_ROOT(:)=0
      TOTAL_FILE_SIZE=0_8
      TOTAL_STRUC_SIZE=0_8
      TMP_OOC_NAMES(:)="?"
      INFO1 = -999
      INFO2 = -999
      INFOG1 = -999
      INFOG2 = -999
      CALL DMUMPS_GET_SAVE_FILES(id,restore_file,INFO_FILE)
      IF ( id%INFO(1) .LT. 0 ) RETURN
      IN=80
      inquire (UNIT=IN,exist=UNIT_OK,opened=UNIT_OP)
      IF(.NOT.UNIT_OK .OR. UNIT_OP) THEN
         id%INFO(1) = -79
         id%INFO(2) = IN
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      open(UNIT=IN,FILE=restore_file
#if defined(MUMPS_F2003)
     &     ,ACCESS="stream"
#endif
     &     ,STATUS='old',form='unformatted',iostat=ierr)
      if(ierr.ne.0) THEN
         id%INFO(1) = -74
         id%INFO(2) = 0
      endif  
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      MP= id%ICNTL(2)
      MPG= id%ICNTL(3)
      CALL DMUMPS_SAVE_RESTORE_STRUCTURE(id,IN,"restore"
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOT,SIZE_VARIABLES_ROOT,SIZE_GEST_ROOT
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)     
      PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. 0 )
      if(id%INFO(1).EQ.0) then
         id%INFO(1)=INFO1
         id%INFO(2)=INFO2
         id%INFOG(1)=INFOG1
         id%INFOG(2)=INFOG2
         if(id%INFO(1).NE.0) then
            write(MPG,*) "Warning: "
     &           ,"restored instance has negative INFO(1):"
     &           , id%INFO(1)
         endif
         if(MP.GT.0) then
            JOB=id%KEEP(40)+456789
            write(MP,*) "Restore done successfully"
            write(MP,*) "From file ",trim(adjustl(restore_file))
            if((id%ICNTL(18).EQ.0).AND.(id%ICNTL(5).EQ.0)) then
               write(MP,*) "with JOB, N, NNZ ",JOB, id%N,id%NNZ
            elseif((id%ICNTL(18).EQ.1).AND.(id%ICNTL(5).EQ. 0)) then
               write(MP,*) "with JOB, N, NNZ_loc=", JOB, id%N,
     &              id%NNZ_loc
            elseif((id%ICNTL(18).EQ.0).AND.(id%ICNTL(5).EQ. 1)) then
               write(MP,*) "with JOB, N, NELT=", JOB, id%N, id%NELT
            endif
         endif
         IF(PROKG) THEN
            IF(id%KEEP(201).EQ.1) THEN
               K=1
               write(MPG,*) "The corresponding OOC files are:"
               DO I=1,id%OOC_NB_FILE_TYPE
                  DO J=1,id%OOC_NB_FILES(I)
                     DO H=1,id%OOC_FILE_NAME_LENGTH(K)-2
                        TMP_OOC_NAMES(H)=id%OOC_FILE_NAMES(K,H)
                     ENDDO
                     write(MPG,*)
     &                    TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(K)-2)
                     K=K+1
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      else
         id%root%gridinit_done=.FALSE.
         id%KEEP(140)=1
      endif
      CLOSE(IN)
      deallocate(SIZE_VARIABLES,SIZE_VARIABLES_ROOT)
      deallocate(SIZE_GEST,SIZE_GEST_ROOT)
      if (id%KEEP(201) .GT. 0) THEN
        id%ASSOCIATED_OOC_FILES=.TRUE.
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_RESTORE
      SUBROUTINE DMUMPS_SAVE_RESTORE_STRUCTURE(id,unit,mode
     &     ,NBVARIABLES,SIZE_VARIABLES,SIZE_GEST
     &     ,NBVARIABLES_ROOT,SIZE_VARIABLES_ROOT,SIZE_GEST_ROOT
     &     ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &     ,INFO1,INFO2,INFOG1,INFOG2)
      USE DMUMPS_FACSOL_L0OMP_M, ONLY : DMUMPS_SAVE_RESTORE_L0FACARRAY
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER,intent(in)::unit,NBVARIABLES,NBVARIABLES_ROOT
      CHARACTER(len=*),intent(in) :: mode
      INTEGER(8),dimension(NBVARIABLES)::SIZE_VARIABLES
      INTEGER(8),dimension(NBVARIABLES_ROOT)::SIZE_VARIABLES_ROOT
      INTEGER,dimension(NBVARIABLES)::SIZE_GEST
      INTEGER,dimension(NBVARIABLES_ROOT)::SIZE_GEST_ROOT
      INTEGER(8) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      INTEGER:: INFO1,INFO2,INFOG1,INFOG2
      INTEGER:: j,i1,i2,err,ierr
      CHARACTER(len=30), allocatable, dimension(:)::VARIABLES
      CHARACTER(len=30), allocatable, dimension(:)::VARIABLES_ROOT
      CHARACTER(len=30) :: TMP_STRING1, TMP_STRING2
      CHARACTER :: ARITH,READ_ARITH
      INTEGER(8) :: size_written,gest_size,WRITTEN_STRUC_SIZE
      INTEGER:: SIZE_INT, SIZE_INT8, SIZE_RL_OR_DBL, SIZE_ARITH_DEP
      INTEGER:: SIZE_DOUBLE_PRECISION, SIZE_LOGICAL, SIZE_CHARACTER
      INTEGER:: READ_NPROCS, READ_PAR, READ_SYM
      INTEGER,dimension(NBVARIABLES)::NbRecords
      INTEGER,dimension(NBVARIABLES_ROOT)::NbRecords_ROOT
      INTEGER:: size_array1,size_array2,dummy,allocok
      INTEGER(8):: size_array_INT8_1,size_array_INT8_2
      LOGICAL:: INT_TYPE_64, READ_INT_TYPE_64
      INTEGER:: tot_NbRecords,NbSubRecords
      INTEGER(8):: size_read,size_allocated
      INTEGER(8),dimension(NBVARIABLES)::DIFF_SIZE_ALLOC_READ
      INTEGER(8),dimension(NBVARIABLES_ROOT)::DIFF_SIZE_ALLOC_READ_ROOT
      INTEGER::READ_OOC_FILE_NAME_LENGTH
      CHARACTER(len=LEN_SAVE_FILE):: READ_OOC_FIRST_FILE_NAME
      INTEGER,dimension(4)::OOC_INDICES
      CHARACTER(len=8) :: date
      CHARACTER(len=10) :: time
      CHARACTER(len=5) :: zone
      INTEGER,dimension(8):: values
      CHARACTER(len=23) :: hash,READ_HASH
      LOGICAL:: BASIC_CHECK
      LOGICAL :: FORTRAN_VERSION_OK
      CHARACTER(len=1) :: TMP_OOC_NAMES(350)
      INTEGER(8)::SIZE_VARIABLES_BLR,SIZE_VARIABLES_FRONT_DATA,
     &            SIZE_VARIABLES_L0FAC
      INTEGER::SIZE_GEST_BLR,SIZE_GEST_FRONT_DATA,SIZE_GEST_L0FAC
      TYPE (DMUMPS_STRUC) :: id
      allocate(VARIABLES(NBVARIABLES),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 100
      VARIABLES(186)="ASSOCIATED_OOC_FILES"
      VARIABLES(185)="pad16"
      VARIABLES(184)="Deficiency"
      VARIABLES(183)="NB_SINGULAR_VALUES"
      VARIABLES(182)="SINGULAR_VALUES"
      VARIABLES(181)="MPITOOMP_PROCS_MAP"
      VARIABLES(180)="L0_OMP_MAPPING"
      VARIABLES(179)="PTR_LEAFS_L0_OMP"
      VARIABLES(178)="PERM_L0_OMP"
      VARIABLES(177)="VIRT_L0_OMP_MAPPING"
      VARIABLES(176)="VIRT_L0_OMP"
      VARIABLES(175)="PHYS_L0_OMP"
      VARIABLES(174)="IPOOL_A_L0_OMP"
      VARIABLES(173)="IPOOL_B_L0_OMP"
      VARIABLES(172)="I8_L0_OMP"
      VARIABLES(171)="I4_L0_OMP"
      VARIABLES(170)="THREAD_LA"
      VARIABLES(169)="LL0_OMP_FACTORS"
      VARIABLES(168)="LL0_OMP_MAPPING"
      VARIABLES(167)="L_VIRT_L0_OMP"
      VARIABLES(166)="L_PHYS_L0_OMP"
      VARIABLES(165)="LPOOL_B_L0_OMP"
      VARIABLES(164)="LPOOL_A_L0_OMP"
      VARIABLES(163)="L0_OMP_FACTORS"
      VARIABLES(162)="BLRARRAY_ENCODING"
      VARIABLES(161)="FDM_F_ENCODING"
      VARIABLES(160)="pad13"
      VARIABLES(159)="NBGRP"
      VARIABLES(158)="LRGROUPS"
      VARIABLES(157)="root"
      VARIABLES(156)="WORKING"
      VARIABLES(155)="IPTR_WORKING"
      VARIABLES(154)="pad14"
      VARIABLES(153)="SUP_PROC"
      VARIABLES(152)="PIVNUL_LIST"
      VARIABLES(151)="OOC_FILE_NAMES"
      VARIABLES(150)="OOC_FILE_NAME_LENGTH"
      VARIABLES(149)="pad12"
      VARIABLES(148)="OOC_NB_FILE_TYPE"
      VARIABLES(147)="OOC_NB_FILES"
      VARIABLES(146)="OOC_TOTAL_NB_NODES"
      VARIABLES(145)="OOC_VADDR"
      VARIABLES(144)="OOC_SIZE_OF_BLOCK"
      VARIABLES(143)="OOC_INODE_SEQUENCE"
      VARIABLES(142)="OOC_MAX_NB_NODES_FOR_ZONE"
      VARIABLES(141)="INSTANCE_NUMBER"
      VARIABLES(140)="CB_SON_SIZE"
      VARIABLES(139)="DKEEP"
      VARIABLES(138)="LWK_USER"
      VARIABLES(137)="NBSA_LOCAL"
      VARIABLES(136)="WK_USER"
      VARIABLES(135)="CROIX_MANU"
      VARIABLES(134)="SCHED_SBTR"
      VARIABLES(133)="SCHED_GRP"
      VARIABLES(132)="SCHED_DEP"
      VARIABLES(131)="SBTR_ID"
      VARIABLES(130)="DEPTH_FIRST_SEQ"
      VARIABLES(129)="DEPTH_FIRST"
      VARIABLES(128)="MY_NB_LEAF"
      VARIABLES(127)="MY_FIRST_LEAF"
      VARIABLES(126)="MY_ROOT_SBTR"
      VARIABLES(125)="COST_TRAV"
      VARIABLES(124)="MEM_SUBTREE"
      VARIABLES(123)="RHSCOMP"
      VARIABLES(122)="POSINRHSCOMP_COL"
      VARIABLES(121)="pad11"
      VARIABLES(120)="POSINRHSCOMP_COL_ALLOC"
      VARIABLES(119)="POSINRHSCOMP_ROW"
      VARIABLES(118)="MEM_DIST"
      VARIABLES(117)="I_AM_CAND"
      VARIABLES(116)="TAB_POS_IN_PERE"
      VARIABLES(115)="FUTURE_NIV2"
      VARIABLES(114)="ISTEP_TO_INIV2"
      VARIABLES(113)="CANDIDATES"
      VARIABLES(112)="ELTPROC"
      VARIABLES(111)="LELTVAR"
      VARIABLES(110)="NELT_loc"
      VARIABLES(109)="DBLARR"
      VARIABLES(108)="INTARR"
      VARIABLES(107)="PROCNODE"
      VARIABLES(106)="S"
      VARIABLES(105)="PTRFAC"
      VARIABLES(104)="PTLUST_S"
      VARIABLES(103)="Step2node"
      VARIABLES(102)="PROCNODE_STEPS"
      VARIABLES(101)="NA"
      VARIABLES(100)="PTRAR"
      VARIABLES(99)="FRTELT"
      VARIABLES(98)="FRTPTR"
      VARIABLES(97)="FILS"
      VARIABLES(96)="DAD_STEPS"
      VARIABLES(95)="FRERE_STEPS"
      VARIABLES(94)="ND_STEPS"
      VARIABLES(93)="NE_STEPS"
      VARIABLES(92)="STEP"
      VARIABLES(91)="NBSA"
      VARIABLES(90)="LNA"
      VARIABLES(89)="KEEP"
      VARIABLES(88)="IS"
      VARIABLES(87)="ASS_IRECV"
      VARIABLES(86)="NSLAVES"
      VARIABLES(85)="NPROCS"
      VARIABLES(84)="MYID"
      VARIABLES(83)="COMM_LOAD"
      VARIABLES(82)="MYID_NODES"
      VARIABLES(81)="COMM_NODES"
      VARIABLES(80)="INST_Number"
      VARIABLES(79)="MAX_SURF_MASTER"
      VARIABLES(78)="KEEP8"
      VARIABLES(77)="pad7"
      VARIABLES(76)="SAVE_PREFIX"
      VARIABLES(75)="SAVE_DIR"
      VARIABLES(74)="WRITE_PROBLEM"
      VARIABLES(73)="OOC_PREFIX"
      VARIABLES(72)="OOC_TMPDIR"
      VARIABLES(71)="VERSION_NUMBER"
      VARIABLES(70)="MAPPING"
      VARIABLES(69)="LISTVAR_SCHUR"
      VARIABLES(68)="SCHUR_CINTERFACE"
      VARIABLES(67)="SCHUR"
      VARIABLES(66)="SIZE_SCHUR"
      VARIABLES(65)="SCHUR_LLD"
      VARIABLES(64)="SCHUR_NLOC"
      VARIABLES(63)="SCHUR_MLOC"
      VARIABLES(62)="NBLOCK"
      VARIABLES(61)="MBLOCK"
      VARIABLES(60)="NPCOL"
      VARIABLES(59)="NPROW"
      VARIABLES(58)="UNS_PERM"
      VARIABLES(57)="SYM_PERM"
      VARIABLES(56)="METIS_OPTIONS"
      VARIABLES(55)="RINFOG"
      VARIABLES(54)="RINFO"
      VARIABLES(53)="CNTL"
      VARIABLES(52)="COST_SUBTREES"
      VARIABLES(51)="INFOG"
      VARIABLES(50)="INFO"
      VARIABLES(49)="ICNTL"
      VARIABLES(48)="pad6"
      VARIABLES(47)="LSOL_loc"
      VARIABLES(46)="LREDRHS"
      VARIABLES(45)="LRHS_loc"
      VARIABLES(44)="Nloc_RHS"
      VARIABLES(43)="NZ_RHS"
      VARIABLES(42)="NRHS"
      VARIABLES(41)="LRHS"
      VARIABLES(40)="IRHS_loc"
      VARIABLES(39)="ISOL_loc"
      VARIABLES(38)="IRHS_PTR"
      VARIABLES(37)="IRHS_SPARSE"
      VARIABLES(36)="RHS_loc"
      VARIABLES(35)="SOL_loc"
      VARIABLES(34)="RHS_SPARSE"
      VARIABLES(33)="REDRHS"
      VARIABLES(32)="RHS"
      VARIABLES(31)="BLKVAR"
      VARIABLES(30)="BLKPTR"
      VARIABLES(29)="pad5"
      VARIABLES(28)="NBLK"
      VARIABLES(27)="PERM_IN"
      VARIABLES(26)="pad4"
      VARIABLES(25)="A_ELT"
      VARIABLES(24)="ELTVAR"
      VARIABLES(23)="ELTPTR"
      VARIABLES(22)="pad3"
      VARIABLES(21)="NELT"
      VARIABLES(20)="pad2"
      VARIABLES(19)="A_loc"
      VARIABLES(18)="JCN_loc"
      VARIABLES(17)="IRN_loc"
      VARIABLES(16)="NNZ_loc"
      VARIABLES(15)="pad1"
      VARIABLES(14)="NZ_loc"
      VARIABLES(13)="pad0"
      VARIABLES(12)="ROWSCA"
      VARIABLES(11)="COLSCA"
      VARIABLES(10)="JCN"
      VARIABLES(9)="IRN"
      VARIABLES(8)="A"
      VARIABLES(7)="NNZ"
      VARIABLES(6)="NZ"
      VARIABLES(5)="N"
      VARIABLES(4)="JOB"
      VARIABLES(3)="PAR"
      VARIABLES(2)="SYM"
      VARIABLES(1)="COMM"
      allocate(VARIABLES_ROOT(NBVARIABLES_ROOT),  stat=allocok)
      if (allocok .GT. 0) THEN
       id%INFO(1)  =-13
       id%INFO(2) = NBVARIABLES_ROOT
      endif
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 100
      VARIABLES_ROOT(35)="rootpad4"
      VARIABLES_ROOT(34)="NB_SINGULAR_VALUES"
      VARIABLES_ROOT(33)="SINGULAR_VALUES"
      VARIABLES_ROOT(32)="SVD_VT"
      VARIABLES_ROOT(31)="SVD_U"
      VARIABLES_ROOT(30)="gridinit_done"
      VARIABLES_ROOT(29)="yes"
      VARIABLES_ROOT(28)="rootpad3"
      VARIABLES_ROOT(27)="QR_RCOND"
      VARIABLES_ROOT(26)="rootpad"
      VARIABLES_ROOT(25)="RHS_ROOT"
      VARIABLES_ROOT(24)="rootpad2"
      VARIABLES_ROOT(23)="QR_TAU"
      VARIABLES_ROOT(22)="SCHUR_POINTER"
      VARIABLES_ROOT(21)="RHS_CNTR_MASTER_ROOT"
      VARIABLES_ROOT(20)="rootpad1"
      VARIABLES_ROOT(19)="IPIV"
      VARIABLES_ROOT(18)="RG2L_COL"
      VARIABLES_ROOT(17)="RG2L_ROW"
      VARIABLES_ROOT(16)="rootpad0"
      VARIABLES_ROOT(15)="LPIV"
      VARIABLES_ROOT(14)="CNTXT_BLACS"
      VARIABLES_ROOT(13)="DESCRIPTOR"
      VARIABLES_ROOT(12)="TOT_ROOT_SIZE"
      VARIABLES_ROOT(11)="ROOT_SIZE"
      VARIABLES_ROOT(10)="RHS_NLOC"
      VARIABLES_ROOT(9)="SCHUR_LLD"
      VARIABLES_ROOT(8)="SCHUR_NLOC"
      VARIABLES_ROOT(7)="SCHUR_MLOC"
      VARIABLES_ROOT(6)="MYCOL"
      VARIABLES_ROOT(5)="MYROW"
      VARIABLES_ROOT(4)="NPCOL"
      VARIABLES_ROOT(3)="NPROW"
      VARIABLES_ROOT(2)="NBLOCK"
      VARIABLES_ROOT(1)="MBLOCK"
      OOC_INDICES=(/147,148,150,151/)
      SIZE_INT = id%KEEP(34)
      SIZE_INT8 = id%KEEP(34)*id%KEEP(10)
      SIZE_RL_OR_DBL = id%KEEP(16)
      SIZE_ARITH_DEP = id%KEEP(35)
      SIZE_DOUBLE_PRECISION = 8
      SIZE_LOGICAL = 4
      SIZE_CHARACTER = 1
      size_written=int(0,kind=8)
      tot_NbRecords=0
      NbRecords(:)=0
      NbRecords_ROOT(:)=0
      size_read=int(0,kind=8)
      size_allocated=int(0,kind=8)
      DIFF_SIZE_ALLOC_READ(:)=0
      DIFF_SIZE_ALLOC_READ_ROOT(:)=0
      WRITTEN_STRUC_SIZE=int(0,kind=8)
      TMP_OOC_NAMES(:)="?"
      SIZE_VARIABLES_BLR=0_8
      SIZE_GEST_BLR=0
      SIZE_VARIABLES_FRONT_DATA=0_8
      SIZE_GEST_FRONT_DATA=0
      SIZE_VARIABLES_L0FAC=0
      SIZE_GEST_L0FAC=0
      if(trim(mode).EQ."memory_save") then
      elseif(trim(mode).EQ."save") then
         write(unit,iostat=err) "MUMPS"
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(5*SIZE_CHARACTER,kind=8)
#if !defined(MUMPS_F2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         IF(id%MYID.EQ.0) THEN
            call date_and_time(date,time,zone,values)
            hash=trim(date)//trim(time)//trim(zone)
         ENDIF
         CALL MPI_BCAST( hash, 23, MPI_CHARACTER, 0, id%COMM, ierr )
         write(unit,iostat=err) hash
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(23*SIZE_CHARACTER,kind=8)
#if !defined(MUMPS_F2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         write(unit,iostat=err) TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(2*SIZE_INT8,kind=8)
#if !defined(MUMPS_F2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         ARITH="DMUMPS"(1:1)
         write(unit,iostat=err) ARITH
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(1,kind=8)
#if !defined(MUMPS_F2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         write(unit,iostat=err) id%SYM,id%PAR,id%NPROCS
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(3*SIZE_INT,kind=8)
#if !defined(MUMPS_F2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         IF(id%KEEP(10).EQ.1) THEN
            INT_TYPE_64=.TRUE.
         ELSE
            INT_TYPE_64=.FALSE.
         ENDIF
         write(unit,iostat=err) INT_TYPE_64
         if(err.ne.0) THEN
            id%INFO(1) = -72
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &           id%INFO(2))
         endif 
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         size_written=size_written+int(SIZE_LOGICAL,kind=8)
#if !defined(MUMPS_F2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*1,kind=8)
#endif
         IF(associated(id%OOC_FILE_NAME_LENGTH).AND.
     &        associated(id%OOC_FILE_NAMES)) THEN
            write(unit,iostat=err) id%OOC_FILE_NAME_LENGTH(1)
            if(err.ne.0) THEN
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              id%INFO(2))
            endif 
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
            size_written=size_written+int(SIZE_INT,kind=8)
#if !defined(MUMPS_F2003)
            size_written=size_written
     &           +int(2*id%KEEP(34)*1,kind=8)
#endif
            TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(1))=
     &           id%OOC_FILE_NAMES(1,1:id%OOC_FILE_NAME_LENGTH(1))
            write(unit,iostat=err)
     &           TMP_OOC_NAMES(1:id%OOC_FILE_NAME_LENGTH(1))
            if(err.ne.0) THEN
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              id%INFO(2))
            endif 
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
            size_written=size_written
     &           +int(id%OOC_FILE_NAME_LENGTH(1)*SIZE_CHARACTER,kind=8)
#if !defined(MUMPS_F2003)
            size_written=size_written
     &           +int(2*id%KEEP(34)*1,kind=8)
#endif
         ELSE
            write(unit,iostat=err) -999
            if(err.ne.0) THEN
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              id%INFO(2))
            endif 
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
            size_written=size_written+int(SIZE_INT,kind=8)
#if !defined(MUMPS_F2003)
            size_written=size_written
     &           +int(2*id%KEEP(34)*1,kind=8)
#endif
            write(unit,iostat=err) -999
            if(err.ne.0) THEN
               id%INFO(1) = -72
               CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &              id%INFO(2))
            endif 
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
            IF ( id%INFO(1) .LT. 0 ) GOTO 100
            size_written=size_written+int(SIZE_INT,kind=8)
#if !defined(MUMPS_F2003)
            size_written=size_written
     &           +int(2*id%KEEP(34)*1,kind=8)
#endif
         ENDIF
      elseif((trim(mode).EQ."restore").OR.
     &        (trim(mode).EQ."restore_ooc")) then
        CALL MUMPS_READ_HEADER(unit,err,size_read,SIZE_INT,SIZE_INT8,
     &       TOTAL_FILE_SIZE, TOTAL_STRUC_SIZE, READ_ARITH,
     &       READ_INT_TYPE_64, READ_OOC_FILE_NAME_LENGTH,
     &        READ_OOC_FIRST_FILE_NAME,READ_HASH,
     &        READ_SYM,READ_PAR,READ_NPROCS,FORTRAN_VERSION_OK)
         if (err.ne.0) THEN
            id%INFO(1) = -75
            CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
         elseif (.NOT.FORTRAN_VERSION_OK) THEN
            id%INFO(1) = -73
            id%INFO(2) = 1
         endif
         CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         BASIC_CHECK = .false.
         IF (trim(mode).EQ."restore_ooc") THEN
           BASIC_CHECK = .true.
         ENDIF
         CALL DMUMPS_CHECK_HEADER(id,BASIC_CHECK,READ_INT_TYPE_64,
     &              READ_HASH, READ_NPROCS,
     &              READ_ARITH, READ_SYM, READ_PAR)
         IF (id%INFO(1) .LT. 0) GOTO 100
      elseif(trim(mode).EQ."fake_restore") then
         read(unit,iostat=err) READ_HASH
         if(err.ne.0) GOTO 100
         read(unit,iostat=err) TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
         if(err.ne.0) GOTO 100
         IF ( id%INFO(1) .LT. 0 ) GOTO 100
         GOTO 200
      else
         CALL MUMPS_ABORT()
      endif
      DO j=1,size(OOC_INDICES)
         i1=OOC_INDICES(j)
         TMP_STRING1 = VARIABLES(i1)
         SELECT CASE(TMP_STRING1)
         CASE("OOC_NB_FILES")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%OOC_NB_FILES)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%OOC_NB_FILES,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%OOC_NB_FILES)) THEN
                  write(unit,iostat=err) size(id%OOC_NB_FILES,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_NB_FILES
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif((trim(mode).EQ."restore").OR.
     &             (trim(mode).EQ."restore_ooc")) then
               nullify(id%OOC_NB_FILES)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%OOC_NB_FILES(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_NB_FILES
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_NB_FILE_TYPE")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%OOC_NB_FILE_TYPE
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif((trim(mode).EQ."restore").OR.
     &             (trim(mode).EQ."restore_ooc")) then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%OOC_NB_FILE_TYPE
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_FILE_NAMES")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%OOC_FILE_NAMES)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%OOC_FILE_NAMES,1)
     &                 *size(id%OOC_FILE_NAMES,2)*SIZE_CHARACTER
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%OOC_FILE_NAMES)) THEN
                  write(unit,iostat=err) size(id%OOC_FILE_NAMES,1)
     &                 ,size(id%OOC_FILE_NAMES,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_FILE_NAMES
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif((trim(mode).EQ."restore").OR.
     &             (trim(mode).EQ."restore_ooc")) then
               nullify(id%OOC_FILE_NAMES)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2
     &                 *SIZE_CHARACTER
                  allocate(id%OOC_FILE_NAMES(size_array1,size_array2),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(size_array1*size_array2
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_FILE_NAMES
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_FILE_NAME_LENGTH")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%OOC_FILE_NAME_LENGTH)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%OOC_FILE_NAME_LENGTH,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%OOC_FILE_NAME_LENGTH)) THEN
                  write(unit,iostat=err) size(id%OOC_FILE_NAME_LENGTH,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_FILE_NAME_LENGTH
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif((trim(mode).EQ."restore").OR.
     &             (trim(mode).EQ."restore_ooc")) then
               nullify(id%OOC_FILE_NAME_LENGTH)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%OOC_FILE_NAME_LENGTH(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_FILE_NAME_LENGTH
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE DEFAULT
         END SELECT
      ENDDO
      if(trim(mode).EQ."restore_ooc") then
         goto 200
      endif
      DO i1=1,NBVARIABLES 
         TMP_STRING1 = VARIABLES(i1)
         SELECT CASE(TMP_STRING1)
         CASE("COMM")
            if(trim(mode).EQ."memory_save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."restore") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE("SYM")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%SYM
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%SYM
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read,
     &                 id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("PAR")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%PAR
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%PAR
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("JOB")
            if(trim(mode).EQ."memory_save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."restore") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE("N")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%N
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%N
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("ICNTL")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%ICNTL,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%ICNTL
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%ICNTL,1)
               read(unit,iostat=err) id%ICNTL
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("INFO")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%INFO,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) INFO1,INFO2
     &              ,id%INFO(3:size(id%INFO,1))
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%INFO,1)
               read(unit,iostat=err) INFO1,INFO2
     &              ,id%INFO(3:size(id%INFO,1))
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("INFOG")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%INFOG,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) INFOG1,INFOG2
     &              ,id%INFOG(3:size(id%INFOG,1))
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%INFOG,1)
               read(unit,iostat=err) INFOG1,INFOG2
     &              ,id%INFOG(3:size(id%INFOG,1))
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("COST_SUBTREES")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%COST_SUBTREES
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL
               read(unit,iostat=err) id%COST_SUBTREES
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("CNTL")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%CNTL,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%CNTL
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%CNTL,1)
               read(unit,iostat=err) id%CNTL
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("RINFO")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%RINFO,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%RINFO
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%RINFO,1)
               read(unit,iostat=err) id%RINFO
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("RINFOG")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%RINFOG,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%RINFOG
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%RINFOG,1)
               read(unit,iostat=err) id%RINFOG
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("KEEP8")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT8*size(id%KEEP8,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%KEEP8
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT8*size(id%KEEP8,1)
               read(unit,iostat=err) id%KEEP8
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("KEEP")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%KEEP,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%KEEP
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%KEEP,1)
               read(unit,iostat=err) id%KEEP
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("DKEEP")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%DKEEP,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%DKEEP
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_RL_OR_DBL*size(id%DKEEP,1)
               read(unit,iostat=err) id%DKEEP
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NZ")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NZ
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NZ
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NNZ")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT8
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NNZ
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT8
               read(unit,iostat=err) id%NNZ
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif  
         CASE("A")
         CASE("IRN")
         CASE("JCN")
         CASE("COLSCA")
            IF(id%KEEP(52).NE.-1) THEN
               NbRecords(i1)=2
               if(trim(mode).EQ."memory_save") then
                  IF(associated(id%COLSCA)) THEN
                     SIZE_GEST(i1)=SIZE_INT
                     SIZE_VARIABLES(i1)=size(id%COLSCA,1)*SIZE_RL_OR_DBL
                  ELSE
                     SIZE_GEST(i1)=SIZE_INT*2
                     SIZE_VARIABLES(i1)=0_8
                  ENDIF
               elseif(trim(mode).EQ."save") then
                  IF(associated(id%COLSCA)) THEN
                     write(unit,iostat=err) size(id%COLSCA,1)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written,
     &                       id%INFO(2))
                     endif
                     write(unit,iostat=err) id%COLSCA
                  ELSE
                     write(unit,iostat=err) -999
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written,
     &                       id%INFO(2))
                     endif
                     write(unit,iostat=err) -999
                  ENDIF 
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  nullify(id%COLSCA)
                  read(unit,iostat=err) size_array1
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(size_array1.EQ.-999) then
                     SIZE_GEST(i1)=SIZE_INT*2
                     SIZE_VARIABLES(i1)=0_8
                     read(unit,iostat=err) dummy
                  else
                     SIZE_GEST(i1)=SIZE_INT
                     SIZE_VARIABLES(i1)=size_array1*SIZE_RL_OR_DBL
                     allocate(id%COLSCA(size_array1), stat=allocok)
                     if (allocok .GT. 0) THEN
                        id%INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                     endif
                     read(unit,iostat=err) id%COLSCA
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            ELSE
            ENDIF
         CASE("ROWSCA")
            IF(id%KEEP(52).NE.-1) THEN
               NbRecords(i1)=2
               if(trim(mode).EQ."memory_save") then
                  IF(associated(id%ROWSCA)) THEN
                     SIZE_GEST(i1)=SIZE_INT
                     SIZE_VARIABLES(i1)=size(id%ROWSCA,1)*SIZE_RL_OR_DBL
                  ELSE
                     SIZE_GEST(i1)=SIZE_INT*2
                     SIZE_VARIABLES(i1)=0_8
                  ENDIF
               elseif(trim(mode).EQ."save") then
                  IF(associated(id%ROWSCA)) THEN
                     write(unit,iostat=err) size(id%ROWSCA,1)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written,
     &                       id%INFO(2))
                     endif
                     write(unit,iostat=err) id%ROWSCA
                  ELSE
                     write(unit,iostat=err) -999
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written,
     &                       id%INFO(2))
                     endif
                     write(unit,iostat=err) -999
                  ENDIF 
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  nullify(id%ROWSCA)
                  read(unit,iostat=err) size_array1
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(size_array1.EQ.-999) then
                     SIZE_GEST(i1)=SIZE_INT*2
                     SIZE_VARIABLES(i1)=0_8
                     read(unit,iostat=err) dummy
                  else
                     SIZE_GEST(i1)=SIZE_INT
                     SIZE_VARIABLES(i1)=size_array1*SIZE_RL_OR_DBL
                     allocate(id%ROWSCA(size_array1), stat=allocok)
                     if (allocok .GT. 0) THEN
                        id%INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                     endif
                     read(unit,iostat=err) id%ROWSCA
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            ELSE
            ENDIF
         CASE("NZ_loc")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NZ_loc
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NZ_loc
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NNZ_loc")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT8
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NNZ_loc
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT8
               read(unit,iostat=err) id%NNZ_loc
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("IRN_loc")
         CASE("JCN_loc")
         CASE("A_loc")
         CASE("NELT")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NELT
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NELT
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NBLK")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NBLK
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NBLK
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("ELTPTR")
         CASE("ELTVAR")
         CASE("A_ELT")
         CASE("PERM_IN")
         CASE("BLKPTR")
         CASE("BLKVAR")
         CASE("RHS")
         CASE("REDRHS")
         CASE("RHS_SPARSE")
         CASE("SOL_loc")
         CASE("RHS_loc")
         CASE("IRHS_SPARSE")
         CASE("IRHS_PTR")
         CASE("ISOL_loc")
         CASE("IRHS_loc")
         CASE("LRHS")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%LRHS
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%LRHS
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NRHS")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NRHS
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NRHS
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NZ_RHS")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NZ_RHS
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NZ_RHS
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("LRHS_loc")
         CASE("Nloc_RHS")
         CASE("LSOL_loc")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%LSOL_loc
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%LSOL_loc
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("LREDRHS")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%LREDRHS
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%LREDRHS
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("SYM_PERM")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               NbRecords(i1)=2
               IF(associated(id%SYM_PERM)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%SYM_PERM,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%SYM_PERM)) THEN
                  write(unit,iostat=err) size(id%SYM_PERM,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%SYM_PERM
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%SYM_PERM)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%SYM_PERM(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%SYM_PERM
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("UNS_PERM")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%UNS_PERM)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%UNS_PERM,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%UNS_PERM)) THEN
                  write(unit,iostat=err) size(id%UNS_PERM,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%UNS_PERM
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%UNS_PERM)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%UNS_PERM(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%UNS_PERM
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NPROW")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NPROW
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NPROW
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NPCOL")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NPCOL
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               NbRecords(i1)=1
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NPCOL
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("MBLOCK")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%MBLOCK
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%MBLOCK
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NBLOCK")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NBLOCK
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NBLOCK
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("SCHUR_MLOC")
             NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%SCHUR_MLOC
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%SCHUR_MLOC
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("SCHUR_NLOC")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%SCHUR_NLOC
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%SCHUR_NLOC
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("SCHUR_LLD")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%SCHUR_LLD
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%SCHUR_LLD
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("SIZE_SCHUR")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%SIZE_SCHUR
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               NbRecords(i1)=1
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%SIZE_SCHUR
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("SCHUR")
         CASE("SCHUR_CINTERFACE")
          CASE("LISTVAR_SCHUR")
         CASE("MAPPING")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%MAPPING)) THEN
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=id%KEEP8(28)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT8+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%MAPPING)) THEN
                  write(unit,iostat=err) id%KEEP8(28)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MAPPING
               ELSE
                  write(unit,iostat=err) int(-999,8)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%MAPPING)
               read(unit,iostat=err) size_array_INT8_1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array_INT8_1.EQ.int(-999,8)) then
                  SIZE_GEST(i1)=SIZE_INT+SIZE_INT8
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=size_array_INT8_1*SIZE_INT
                  allocate(id%MAPPING(size_array_INT8_1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MAPPING
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("VERSION_NUMBER")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=len(id%VERSION_NUMBER)
     &              *SIZE_CHARACTER
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%VERSION_NUMBER
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=len(id%VERSION_NUMBER)
     &              *SIZE_CHARACTER
               read(unit,iostat=err) id%VERSION_NUMBER
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_TMPDIR")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=len(id%OOC_TMPDIR)
     &              *SIZE_CHARACTER
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%OOC_TMPDIR
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=len(id%OOC_TMPDIR)
     &              *SIZE_CHARACTER
               read(unit,iostat=err) id%OOC_TMPDIR
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_PREFIX")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=len(id%OOC_PREFIX)
     &              *SIZE_CHARACTER
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%OOC_PREFIX
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               NbRecords(i1)=1
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=len(id%OOC_PREFIX)
     &              *SIZE_CHARACTER
               read(unit,iostat=err) id%OOC_PREFIX
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("WRITE_PROBLEM")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=len(id%WRITE_PROBLEM)
     &              *SIZE_CHARACTER
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%WRITE_PROBLEM
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=len(id%WRITE_PROBLEM)
     &              *SIZE_CHARACTER
               read(unit,iostat=err) id%WRITE_PROBLEM
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("MAX_SURF_MASTER")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT8
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%MAX_SURF_MASTER
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT8
               read(unit,iostat=err) id%MAX_SURF_MASTER
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("INST_Number")
               NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%INST_Number
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%INST_Number
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("COMM_NODES")
            if(trim(mode).EQ."memory_save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."restore") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE("MYID_NODES")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%MYID_NODES
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%MYID_NODES
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("COMM_LOAD")
            if(trim(mode).EQ."memory_save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."restore") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE("MYID")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%MYID
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%MYID
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NPROCS")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NPROCS
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NPROCS
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NSLAVES")
             NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NSLAVES
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NSLAVES
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("ASS_IRECV")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%ASS_IRECV
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%ASS_IRECV
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif  
         CASE("IS")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%IS)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=id%KEEP(32)*SIZE_INT
                  DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT*
     &                 (size(id%IS,1)-id%KEEP(32))
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%IS)) THEN
                  write(unit,iostat=err) size(id%IS,1),id%KEEP(32)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%IS(1:id%KEEP(32))
                  DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT*
     &                 (size(id%IS,1)-id%KEEP(32))
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%IS)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array2*SIZE_INT
                  DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT*
     &                 (size_array1-size_array2)
                  allocate(id%IS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%IS(1:size_array2)
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("Deficiency")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%Deficiency
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%Deficiency
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("LNA")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%LNA
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%LNA
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NBSA")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NBSA
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NBSA
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("STEP")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%STEP)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%STEP,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%STEP)) THEN
                  write(unit,iostat=err) size(id%STEP,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%STEP
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%STEP)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(SIZE_VARIABLES(i1),id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%STEP(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%STEP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NE_STEPS")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%NE_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%NE_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%NE_STEPS)) THEN
                  write(unit,iostat=err) size(id%NE_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%NE_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%NE_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%NE_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%NE_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("ND_STEPS")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%ND_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%ND_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%ND_STEPS)) THEN
                  write(unit,iostat=err) size(id%ND_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%ND_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%ND_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%ND_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%ND_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("Step2node")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%Step2node)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%Step2node,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%Step2node)) THEN
                  write(unit,iostat=err) size(id%Step2node,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%Step2node
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%Step2node)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%Step2node(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%Step2node
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("FRERE_STEPS")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%FRERE_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%FRERE_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%FRERE_STEPS)) THEN
                  write(unit,iostat=err) size(id%FRERE_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%FRERE_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%FRERE_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%FRERE_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%FRERE_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("DAD_STEPS")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%DAD_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%DAD_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%DAD_STEPS)) THEN
                  write(unit,iostat=err) size(id%DAD_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%DAD_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%DAD_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%DAD_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%DAD_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("FILS")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%FILS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%FILS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%FILS)) THEN
                  write(unit,iostat=err) size(id%FILS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%FILS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%FILS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%FILS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%FILS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("PTRAR")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%PTRAR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PTRAR,1)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%PTRAR)) THEN
                  write(unit,iostat=err) size(id%PTRAR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PTRAR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               NbRecords(i1)=2
            elseif(trim(mode).EQ."restore") then
               nullify(id%PTRAR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT8
                  allocate(id%PTRAR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PTRAR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("FRTPTR")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%FRTPTR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%FRTPTR,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%FRTPTR)) THEN
                  write(unit,iostat=err) size(id%FRTPTR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%FRTPTR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%FRTPTR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%FRTPTR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%FRTPTR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("FRTELT")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%FRTELT)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%FRTELT,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%FRTELT)) THEN
                  write(unit,iostat=err) size(id%FRTELT,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  write(unit,iostat=err) id%FRTELT
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%FRTELT)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%FRTELT(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%FRTELT
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("NA")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%NA)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%NA,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%NA)) THEN
                  write(unit,iostat=err) size(id%NA,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  write(unit,iostat=err) id%NA
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%NA)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%NA(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%NA
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("PROCNODE_STEPS")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               NbRecords(i1)=2
               IF(associated(id%PROCNODE_STEPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PROCNODE_STEPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%PROCNODE_STEPS)) THEN
                  write(unit,iostat=err) size(id%PROCNODE_STEPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PROCNODE_STEPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%PROCNODE_STEPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PROCNODE_STEPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PROCNODE_STEPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("PTLUST_S")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%PTLUST_S)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PTLUST_S,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%PTLUST_S)) THEN
                  write(unit,iostat=err) size(id%PTLUST_S,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PTLUST_S
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%PTLUST_S)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PTLUST_S(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PTLUST_S
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("PTRFAC")
            NbRecords(i1)=2  
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%PTRFAC)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PTRFAC,1)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%PTRFAC)) THEN
                  write(unit,iostat=err) size(id%PTRFAC,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  write(unit,iostat=err) id%PTRFAC
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%PTRFAC)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT8
                  allocate(id%PTRFAC(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PTRFAC
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("S")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%S)) THEN
                  SIZE_GEST(i1)=SIZE_INT8*2
                  SIZE_VARIABLES(i1)=id%KEEP8(31)*SIZE_ARITH_DEP
                  DIFF_SIZE_ALLOC_READ(i1)=
     &                 SIZE_ARITH_DEP*(id%KEEP8(23)-id%KEEP8(31))
               ELSE
                  SIZE_GEST(i1)=SIZE_INT8*2+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%S)) THEN
                  write(unit,iostat=err) id%KEEP8(23),id%KEEP8(31)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%S(1:id%KEEP8(31))
                  DIFF_SIZE_ALLOC_READ(i1)=
     &                 SIZE_ARITH_DEP*(id%KEEP8(23)-id%KEEP8(31))
               ELSE
                  write(unit,iostat=err) int(-999,kind=8)
     &                 ,int(-998,kind=8)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%S)
               read(unit,iostat=err) size_array_INT8_1,size_array_INT8_2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array_INT8_1.EQ.int(-999,kind=8)) then
                  SIZE_GEST(i1)=SIZE_INT8*2+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT8*2
                  SIZE_VARIABLES(i1)=size_array_INT8_2*SIZE_ARITH_DEP
                  DIFF_SIZE_ALLOC_READ(i1)=
     &                 SIZE_ARITH_DEP*
     &                 (size_array_INT8_1-size_array_INT8_2)
                  allocate(id%S(1:size_array_INT8_1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(size_array_INT8_1,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%S(1:size_array_INT8_2)
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("PROCNODE")
         CASE("INTARR")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%INTARR)) THEN
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=id%KEEP8(27)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT8+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%INTARR)) THEN
                  write(unit,iostat=err) id%KEEP8(27)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%INTARR
               ELSE
                  write(unit,iostat=err) int(-999,8)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%INTARR)
               read(unit,iostat=err) size_array_INT8_1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array_INT8_1.EQ.int(-999,8)) then
                  SIZE_GEST(i1)=SIZE_INT8+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=size_array_INT8_1*SIZE_INT
                  allocate(id%INTARR(size_array_INT8_1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%INTARR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("DBLARR")
         CASE("NELT_loc")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NELT_loc
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NELT_loc
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("LELTVAR")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%LELTVAR
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%LELTVAR
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("ELTPROC")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%ELTPROC)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%ELTPROC,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%ELTPROC)) THEN
                  write(unit,iostat=err) size(id%ELTPROC,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%ELTPROC
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%ELTPROC)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%ELTPROC(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%ELTPROC
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("I4_L0_OMP")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%I4_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%I4_L0_OMP,1)
     &                 *size(id%I4_L0_OMP,2)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%I4_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%I4_L0_OMP,1)
     &                 ,size(id%I4_L0_OMP,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%I4_L0_OMP
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%I4_L0_OMP)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT
                  allocate(id%I4_L0_OMP(size_array1,size_array2)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(size_array1*size_array2
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%I4_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("I8_L0_OMP")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%I8_L0_OMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%I8_L0_OMP,1)
     &                 *size(id%I8_L0_OMP,2)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%I8_L0_OMP)) THEN
                  write(unit,iostat=err) size(id%I8_L0_OMP,1)
     &                 ,size(id%I8_L0_OMP,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%I8_L0_OMP
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%I8_L0_OMP)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT8
                  allocate(id%I8_L0_OMP(size_array1,size_array2)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(size_array1*size_array2
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%I8_L0_OMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("CANDIDATES")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%CANDIDATES)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%CANDIDATES,1)
     &                 *size(id%CANDIDATES,2)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%CANDIDATES)) THEN
                  write(unit,iostat=err) size(id%CANDIDATES,1)
     &                 ,size(id%CANDIDATES,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%CANDIDATES
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%CANDIDATES)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT
                  allocate(id%CANDIDATES(size_array1,size_array2)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(size_array1*size_array2
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%CANDIDATES
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("ISTEP_TO_INIV2")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%ISTEP_TO_INIV2)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%ISTEP_TO_INIV2,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%ISTEP_TO_INIV2)) THEN
                  write(unit,iostat=err) size(id%ISTEP_TO_INIV2,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%ISTEP_TO_INIV2
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%ISTEP_TO_INIV2)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%ISTEP_TO_INIV2(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%ISTEP_TO_INIV2
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("FUTURE_NIV2")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%FUTURE_NIV2)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%FUTURE_NIV2,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%FUTURE_NIV2)) THEN
                  write(unit,iostat=err) size(id%FUTURE_NIV2,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%FUTURE_NIV2
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%FUTURE_NIV2)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%FUTURE_NIV2(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%FUTURE_NIV2
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("TAB_POS_IN_PERE")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%TAB_POS_IN_PERE)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%TAB_POS_IN_PERE,1)
     &                 *size(id%TAB_POS_IN_PERE,2)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%TAB_POS_IN_PERE)) THEN
                  write(unit,iostat=err) size(id%TAB_POS_IN_PERE,1)
     &                 ,size(id%TAB_POS_IN_PERE,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%TAB_POS_IN_PERE
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%TAB_POS_IN_PERE)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT
                  allocate(id%TAB_POS_IN_PERE(size_array1,size_array2)
     &                 , stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(size_array1*size_array2
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%TAB_POS_IN_PERE
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("I_AM_CAND")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%I_AM_CAND)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%I_AM_CAND,1)*SIZE_LOGICAL
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%I_AM_CAND)) THEN
                  write(unit,iostat=err) size(id%I_AM_CAND,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%I_AM_CAND
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%I_AM_CAND)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_LOGICAL
                  allocate(id%I_AM_CAND(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%I_AM_CAND
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("MEM_DIST")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%MEM_DIST)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%MEM_DIST,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%MEM_DIST)) THEN
                  write(unit,iostat=err) size(id%MEM_DIST,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  write(unit,iostat=err) id%MEM_DIST
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%MEM_DIST)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MEM_DIST(0:size_array1-1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MEM_DIST
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("POSINRHSCOMP_ROW")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%POSINRHSCOMP_ROW)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%POSINRHSCOMP_ROW,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%POSINRHSCOMP_ROW)) THEN
                  write(unit,iostat=err) size(id%POSINRHSCOMP_ROW,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%POSINRHSCOMP_ROW
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%POSINRHSCOMP_ROW)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%POSINRHSCOMP_ROW(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%POSINRHSCOMP_ROW
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("POSINRHSCOMP_COL_ALLOC")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_LOGICAL
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%POSINRHSCOMP_COL_ALLOC
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_LOGICAL
               read(unit,iostat=err) id%POSINRHSCOMP_COL_ALLOC
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("POSINRHSCOMP_COL")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%POSINRHSCOMP_COL)) THEN
                  IF(id%POSINRHSCOMP_COL_ALLOC) THEN
                     SIZE_GEST(i1)=SIZE_INT
                     SIZE_VARIABLES(i1)=
     &                    size(id%POSINRHSCOMP_COL,1)*SIZE_INT
                  ELSE
                     SIZE_GEST(i1)=SIZE_INT*2
                     SIZE_VARIABLES(i1)=0_8
                  ENDIF
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%POSINRHSCOMP_COL)) THEN
                  IF(id%POSINRHSCOMP_COL_ALLOC) THEN
                     write(unit,iostat=err) size(id%POSINRHSCOMP_COL,1)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written,
     &                       id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) id%POSINRHSCOMP_COL
                  ELSE
                     write(unit,iostat=err) size(id%POSINRHSCOMP_COL,1)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written,
     &                       id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) -999
                  ENDIF
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%POSINRHSCOMP_COL)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  if(id%POSINRHSCOMP_COL_ALLOC) then
                     SIZE_GEST(i1)=SIZE_INT
                     SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                     allocate(id%POSINRHSCOMP_COL(size_array1), 
     &                    stat=allocok)
                     if (allocok .GT. 0) THEN
                        id%INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                     endif
                     read(unit,iostat=err) id%POSINRHSCOMP_COL
                  else
                     SIZE_GEST(i1)=SIZE_INT*2
                     SIZE_VARIABLES(i1)=0_8
                     read(unit,iostat=err) dummy
                     id%POSINRHSCOMP_COL=>id%POSINRHSCOMP_ROW
                  endif
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("RHSCOMP")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%RHSCOMP)) THEN
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=id%KEEP8(25)*SIZE_ARITH_DEP
               ELSE
                  SIZE_GEST(i1)=SIZE_INT8+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%RHSCOMP)) THEN
                  write(unit,iostat=err) id%KEEP8(25)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%RHSCOMP
               ELSE
                  write(unit,iostat=err) int(-999,8)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%RHSCOMP)
               read(unit,iostat=err) size_array_INT8_1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array_INT8_1.EQ.int(-999,8)) then
                  SIZE_GEST(i1)=SIZE_INT8+SIZE_INT
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT8
                  SIZE_VARIABLES(i1)=size_array_INT8_1*SIZE_ARITH_DEP
                  allocate(id%RHSCOMP(size_array_INT8_1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%RHSCOMP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("MEM_SUBTREE")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%MEM_SUBTREE)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%MEM_SUBTREE,1)*SIZE_DOUBLE_PRECISION
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%MEM_SUBTREE)) THEN
                  write(unit,iostat=err) size(id%MEM_SUBTREE,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MEM_SUBTREE
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%MEM_SUBTREE)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_DOUBLE_PRECISION
                  allocate(id%MEM_SUBTREE(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MEM_SUBTREE
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("COST_TRAV")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%COST_TRAV)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%COST_TRAV,1)*SIZE_DOUBLE_PRECISION
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%COST_TRAV)) THEN
                  write(unit,iostat=err) size(id%COST_TRAV,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%COST_TRAV
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%COST_TRAV)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_DOUBLE_PRECISION
                  allocate(id%COST_TRAV(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%COST_TRAV
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("MY_ROOT_SBTR")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%MY_ROOT_SBTR)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%MY_ROOT_SBTR,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%MY_ROOT_SBTR)) THEN
                  write(unit,iostat=err) size(id%MY_ROOT_SBTR,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MY_ROOT_SBTR
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%MY_ROOT_SBTR)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MY_ROOT_SBTR(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MY_ROOT_SBTR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("MY_FIRST_LEAF")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%MY_FIRST_LEAF)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%MY_FIRST_LEAF,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%MY_FIRST_LEAF)) THEN
                  write(unit,iostat=err) size(id%MY_FIRST_LEAF,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MY_FIRST_LEAF
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%MY_FIRST_LEAF)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MY_FIRST_LEAF(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MY_FIRST_LEAF
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("MY_NB_LEAF")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%MY_NB_LEAF)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%MY_NB_LEAF,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%MY_NB_LEAF)) THEN
                  write(unit,iostat=err) size(id%MY_NB_LEAF,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MY_NB_LEAF
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%MY_NB_LEAF)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MY_NB_LEAF(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MY_NB_LEAF
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("DEPTH_FIRST")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%DEPTH_FIRST)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%DEPTH_FIRST,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%DEPTH_FIRST)) THEN
                  write(unit,iostat=err) size(id%DEPTH_FIRST,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%DEPTH_FIRST
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%DEPTH_FIRST)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%DEPTH_FIRST(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%DEPTH_FIRST
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("DEPTH_FIRST_SEQ")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%DEPTH_FIRST_SEQ)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%DEPTH_FIRST_SEQ,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%DEPTH_FIRST_SEQ)) THEN
                  write(unit,iostat=err) size(id%DEPTH_FIRST_SEQ,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%DEPTH_FIRST_SEQ
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%DEPTH_FIRST_SEQ)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%DEPTH_FIRST_SEQ(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%DEPTH_FIRST_SEQ
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("SBTR_ID")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%SBTR_ID)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%SBTR_ID,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%SBTR_ID)) THEN
                  write(unit,iostat=err) size(id%SBTR_ID,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%SBTR_ID
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%SBTR_ID)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%SBTR_ID(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%SBTR_ID
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("SCHED_DEP")
         CASE("SCHED_GRP")
         CASE("CROIX_MANU")
         CASE("WK_USER")
         CASE("NBSA_LOCAL")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NBSA_LOCAL
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NBSA_LOCAL
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("LWK_USER")
            if(trim(mode).EQ."memory_save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            elseif(trim(mode).EQ."restore") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_INT
            endif
         CASE("CB_SON_SIZE")
         CASE("INSTANCE_NUMBER")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%INSTANCE_NUMBER
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%INSTANCE_NUMBER
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_MAX_NB_NODES_FOR_ZONE")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%OOC_MAX_NB_NODES_FOR_ZONE
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%OOC_MAX_NB_NODES_FOR_ZONE
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_INODE_SEQUENCE")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%OOC_INODE_SEQUENCE)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%OOC_INODE_SEQUENCE,1)
     &                 *size(id%OOC_INODE_SEQUENCE,2)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%OOC_INODE_SEQUENCE)) THEN
                  write(unit,iostat=err) size(id%OOC_INODE_SEQUENCE,1)
     &                 ,size(id%OOC_INODE_SEQUENCE,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_INODE_SEQUENCE
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%OOC_INODE_SEQUENCE)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT
                  allocate(id%OOC_INODE_SEQUENCE(size_array1
     &                 ,size_array2), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(size_array1*size_array2
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_INODE_SEQUENCE
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_SIZE_OF_BLOCK")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%OOC_SIZE_OF_BLOCK)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%OOC_SIZE_OF_BLOCK,1)
     &                 *size(id%OOC_SIZE_OF_BLOCK,2)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%OOC_SIZE_OF_BLOCK)) THEN
                  write(unit,iostat=err) size(id%OOC_SIZE_OF_BLOCK,1)
     &                 ,size(id%OOC_SIZE_OF_BLOCK,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_SIZE_OF_BLOCK
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%OOC_SIZE_OF_BLOCK)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT8
                  allocate(id%OOC_SIZE_OF_BLOCK(size_array1
     &                 ,size_array2), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(size_array1*size_array2
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_SIZE_OF_BLOCK
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_VADDR")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%OOC_VADDR)) THEN
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size(id%OOC_VADDR,1)
     &                 *size(id%OOC_VADDR,2)*SIZE_INT8
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%OOC_VADDR)) THEN
                  write(unit,iostat=err) size(id%OOC_VADDR,1)
     &                 ,size(id%OOC_VADDR,2)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_VADDR
               ELSE
                  write(unit,iostat=err) -999,-998
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%OOC_VADDR)
               read(unit,iostat=err) size_array1,size_array2
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*3
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=size_array1*size_array2*SIZE_INT8
                  allocate(id%OOC_VADDR(size_array1,size_array2),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(size_array1*size_array2
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_VADDR
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_TOTAL_NB_NODES")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%OOC_TOTAL_NB_NODES)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%OOC_TOTAL_NB_NODES,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%OOC_TOTAL_NB_NODES)) THEN
                  write(unit,iostat=err) size(id%OOC_TOTAL_NB_NODES,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%OOC_TOTAL_NB_NODES
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%OOC_TOTAL_NB_NODES)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%OOC_TOTAL_NB_NODES(size_array1), 
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%OOC_TOTAL_NB_NODES
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("OOC_NB_FILES")
         CASE("OOC_NB_FILE_TYPE")
         CASE("OOC_FILE_NAMES")
         CASE("OOC_FILE_NAME_LENGTH")
         CASE("PIVNUL_LIST")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%PIVNUL_LIST)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%PIVNUL_LIST,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%PIVNUL_LIST)) THEN
                  write(unit,iostat=err) size(id%PIVNUL_LIST,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%PIVNUL_LIST
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%PIVNUL_LIST)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%PIVNUL_LIST(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%PIVNUL_LIST
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("SUP_PROC")
         CASE("IPTR_WORKING")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%IPTR_WORKING)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%IPTR_WORKING,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%IPTR_WORKING)) THEN
                  write(unit,iostat=err) size(id%IPTR_WORKING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%IPTR_WORKING
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%IPTR_WORKING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%IPTR_WORKING(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%IPTR_WORKING
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("WORKING")
            NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%WORKING)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%WORKING,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%WORKING)) THEN
                  write(unit,iostat=err) size(id%WORKING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%WORKING
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%WORKING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%WORKING(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%WORKING
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("root")
            DO i2=1,NBVARIABLES_ROOT
            TMP_STRING2 = VARIABLES_ROOT(i2)
            SELECT CASE(TMP_STRING2)
            CASE("MBLOCK")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%MBLOCK
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%MBLOCK
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("NBLOCK")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%NBLOCK
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100        
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%NBLOCK
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
            CASE("NPROW")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%NPROW
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100                 
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%NPROW
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("NPCOL")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%NPCOL
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100                
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%NPCOL
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("MYROW")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  NbRecords_ROOT(i2)=1
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%MYROW
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%MYROW
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("MYCOL")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%MYCOL
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%MYCOL
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("SCHUR_MLOC")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%SCHUR_MLOC
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%SCHUR_MLOC
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("SCHUR_NLOC")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%SCHUR_NLOC
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%SCHUR_NLOC
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("SCHUR_LLD")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%SCHUR_LLD
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%SCHUR_LLD
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("RHS_NLOC")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%RHS_NLOC
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%RHS_NLOC
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("ROOT_SIZE")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%ROOT_SIZE
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%ROOT_SIZE
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("TOT_ROOT_SIZE")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%TOT_ROOT_SIZE
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%TOT_ROOT_SIZE
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("DESCRIPTOR")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=size(id%root%DESCRIPTOR,1)
     &                 *SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%DESCRIPTOR
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT*
     &                 size(id%root%DESCRIPTOR,1)
                  read(unit,iostat=err) id%root%DESCRIPTOR
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("CNTXT_BLACS")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%CNTXT_BLACS
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%CNTXT_BLACS
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("LPIV")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%LPIV
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%LPIV
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("RG2L_ROW")
            CASE("RG2L_COL")
            CASE("IPIV")
               NbRecords_ROOT(i2)=2
               if(trim(mode).EQ."memory_save") then
                  IF(associated(id%root%IPIV)) THEN
                     SIZE_GEST_ROOT(i2)=SIZE_INT
                     SIZE_VARIABLES_ROOT(i2)=
     &                    size(id%root%IPIV,1)*SIZE_INT
                  ELSE
                     SIZE_GEST_ROOT(i2)=SIZE_INT*2
                     SIZE_VARIABLES_ROOT(i2)=0_8
                  ENDIF
               elseif(trim(mode).EQ."save") then
                  IF(associated(id%root%IPIV)) THEN
                     write(unit,iostat=err) size(id%root%IPIV,1)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written
     &                       ,id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) id%root%IPIV
                  ELSE
                     write(unit,iostat=err) -999
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written
     &                       ,id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) -999
                  ENDIF 
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  nullify(id%root%IPIV)
                  read(unit,iostat=err) size_array1
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(size_array1.EQ.-999) then
                     SIZE_GEST_ROOT(i2)=SIZE_INT*2
                     SIZE_VARIABLES_ROOT(i2)=0_8
                     read(unit,iostat=err) dummy
                  else
                     SIZE_GEST_ROOT(i2)=SIZE_INT
                     SIZE_VARIABLES_ROOT(i2)=size_array1*SIZE_INT
                     allocate(id%root%IPIV(size_array1), stat=allocok)
                     if (allocok .GT. 0) THEN
                        id%INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_STRUC_SIZE-size_allocated
     &                       ,id%INFO(2))
                     endif
                     read(unit,iostat=err) id%root%IPIV
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("RHS_CNTR_MASTER_ROOT")
               NbRecords_ROOT(i2)=2
               if(trim(mode).EQ."memory_save") then
                  IF(associated(id%root%RHS_CNTR_MASTER_ROOT)) THEN
                     SIZE_GEST_ROOT(i2)=SIZE_INT
                     SIZE_VARIABLES_ROOT(i2)=
     &                    size(id%root%RHS_CNTR_MASTER_ROOT,1)
     &                    *SIZE_ARITH_DEP
                  ELSE
                     SIZE_GEST_ROOT(i2)=SIZE_INT*2
                     SIZE_VARIABLES_ROOT(i2)=0_8
                  ENDIF
               elseif(trim(mode).EQ."save") then
                  IF(associated(id%root%RHS_CNTR_MASTER_ROOT)) THEN
                     write(unit,iostat=err) 
     &                    size(id%root%RHS_CNTR_MASTER_ROOT,1)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written
     &                       ,id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) id%root%RHS_CNTR_MASTER_ROOT
                  ELSE
                     write(unit,iostat=err) -999
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written
     &                       ,id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) -999
                  ENDIF 
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  nullify(id%root%RHS_CNTR_MASTER_ROOT)
                  read(unit,iostat=err) size_array1
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(size_array1.EQ.-999) then
                     SIZE_GEST_ROOT(i2)=SIZE_INT*2
                     SIZE_VARIABLES_ROOT(i2)=0_8
                     read(unit,iostat=err) dummy
                  else
                     SIZE_GEST_ROOT(i2)=SIZE_INT
                     SIZE_VARIABLES_ROOT(i2)=size_array1*SIZE_ARITH_DEP
                     allocate(id%root%RHS_CNTR_MASTER_ROOT(size_array1), 
     &                    stat=allocok)
                     if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                     endif
                     read(unit,iostat=err) id%root%RHS_CNTR_MASTER_ROOT
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("SCHUR_POINTER")
            CASE("QR_TAU")
            CASE("RHS_ROOT")
               NbRecords_ROOT(i2)=2
               if(trim(mode).EQ."memory_save") then
                  IF(associated(id%root%RHS_ROOT)) THEN
                     SIZE_GEST_ROOT(i2)=SIZE_INT*2
                     SIZE_VARIABLES_ROOT(i2)=size(id%root%RHS_ROOT,1)
     &                    *size(id%root%RHS_ROOT,2)*SIZE_ARITH_DEP
                  ELSE
                     SIZE_GEST_ROOT(i2)=SIZE_INT*3
                     SIZE_VARIABLES_ROOT(i2)=0_8
                  ENDIF
               elseif(trim(mode).EQ."save") then
                  IF(associated(id%root%RHS_ROOT)) THEN
                     write(unit,iostat=err) size(id%root%RHS_ROOT,1)
     &                    ,size(id%root%RHS_ROOT,2)
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written
     &                       ,id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) id%root%RHS_ROOT
                  ELSE
                     write(unit,iostat=err) -999,-998
                     if(err.ne.0) then
                        id%INFO(1) = -72
                        CALL MUMPS_SETI8TOI4(
     &                       TOTAL_FILE_SIZE-size_written
     &                       ,id%INFO(2))
                     endif
                     CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
                     IF ( id%INFO(1) .LT. 0 ) GOTO 100
                     write(unit,iostat=err) -999
                  ENDIF 
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  nullify(id%root%RHS_ROOT)
                  read(unit,iostat=err) size_array1,size_array2
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(size_array1.EQ.-999) then
                     SIZE_GEST_ROOT(i2)=SIZE_INT*3
                     SIZE_VARIABLES_ROOT(i2)=0_8
                     read(unit,iostat=err) dummy
                  else
                     SIZE_GEST_ROOT(i2)=SIZE_INT*2
                     SIZE_VARIABLES_ROOT(i2)=size_array1*size_array2
     &                    *SIZE_ARITH_DEP
                     allocate(id%root%RHS_ROOT(size_array1,size_array2), 
     &                    stat=allocok)
                     if (allocok .GT. 0) THEN
                        id%INFO(1) = -78
                        CALL MUMPS_SETI8TOI4(size_array1*size_array2
     &                       ,id%INFO(2))
                     endif
                     read(unit,iostat=err) id%root%RHS_ROOT
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("QR_RCOND")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_RL_OR_DBL
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%QR_RCOND
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_RL_OR_DBL
                  read(unit,iostat=err) id%root%QR_RCOND
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("yes")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_LOGICAL
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%yes
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_LOGICAL
                  read(unit,iostat=err) id%root%yes
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("gridinit_done")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_LOGICAL
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%gridinit_done
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_LOGICAL
                  read(unit,iostat=err) id%root%gridinit_done
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2),
     &                    id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("SVD_U")
            CASE("SVD_VT")
            CASE("SINGULAR_VALUES")
            CASE("NB_SINGULAR_VALUES")
               NbRecords_ROOT(i2)=1
               if(trim(mode).EQ."memory_save") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
               elseif(trim(mode).EQ."save") then
                  write(unit,iostat=err) id%root%NB_SINGULAR_VALUES
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written
     &                    ,id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               elseif(trim(mode).EQ."restore") then
                  SIZE_VARIABLES_ROOT(i2)=SIZE_INT
                  read(unit,iostat=err) id%root%NB_SINGULAR_VALUES
                  if(err.ne.0) THEN
                     id%INFO(1) = -75
                     CALL MUMPS_SETI8TOI4(SIZE_VARIABLES_ROOT(i2)
     &                    ,id%INFO(2))
                  endif 
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               endif
            CASE("rootpad0","rootpad1","rootpad2","rootpad",
     &              "rootpad3","rootpad4")
            CASE DEFAULT
            END SELECT
            if(trim(mode).EQ."memory_save") then
               NbSubRecords=int(SIZE_VARIABLES_ROOT(i2)/huge(0))
               IF(NbSubRecords.GT.0) then
                  NbRecords_ROOT(i2)=NbRecords_ROOT(i2)+NbSubRecords
               ENDIF
            elseif(trim(mode).EQ."save") then
               size_written=size_written+SIZE_VARIABLES_ROOT(i2)
     &              +SIZE_GEST_ROOT(i2)
#if !defined(MUMPS_F2003)
               size_written=size_written
     &              +int(2*id%KEEP(34)*NbRecords_ROOT(i2),kind=8)
#endif
            elseif(trim(mode).EQ."restore") then
               size_allocated=size_allocated+SIZE_VARIABLES_ROOT(i2)+
     &              DIFF_SIZE_ALLOC_READ_ROOT(i2)
               size_read=size_read+SIZE_VARIABLES_ROOT(i2)
     &              +int(SIZE_GEST_ROOT(i2),kind=8)
#if !defined(MUMPS_F2003)
               size_read=size_read
     &              +int(2*id%KEEP(34)*NbRecords_ROOT(i2),kind=8)
#endif
            elseif(trim(mode).EQ."fake_restore") then
            endif
         ENDDO
         CASE("NBGRP")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%NBGRP
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%NBGRP
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
        CASE("LRGROUPS")
           NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%LRGROUPS)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size(id%LRGROUPS,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%LRGROUPS)) THEN
                  write(unit,iostat=err) size(id%LRGROUPS,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%LRGROUPS
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%LRGROUPS)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%LRGROUPS(size_array1), stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%LRGROUPS
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("FDM_F_ENCODING")
            NbRecords(i1)=1
            SIZE_GEST(i1)=SIZE_INT
            SIZE_VARIABLES(i1)=0_8
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%FDM_F_ENCODING)) THEN
                  CALL MUMPS_SAVE_RESTORE_FRONT_DATA(id%FDM_F_ENCODING
     &                 ,unit,id%MYID,"memory_save"
     &                 ,SIZE_GEST_FRONT_DATA,SIZE_VARIABLES_FRONT_DATA
     &                 ,SIZE_INT,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))                 
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%FDM_F_ENCODING)) THEN
                  write(unit,iostat=err) size(id%FDM_F_ENCODING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  CALL MUMPS_SAVE_RESTORE_FRONT_DATA(id%FDM_F_ENCODING
     &                 ,unit,id%MYID,"save"
     &                 ,SIZE_GEST_FRONT_DATA,SIZE_VARIABLES_FRONT_DATA
     &                 ,SIZE_INT,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%FDM_F_ENCODING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.NE.-999) then
                  CALL MUMPS_SAVE_RESTORE_FRONT_DATA(id%FDM_F_ENCODING
     &                 ,unit,id%MYID,"restore"
     &                 ,SIZE_GEST_FRONT_DATA,SIZE_VARIABLES_FRONT_DATA
     &                 ,SIZE_INT,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))   
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("BLRARRAY_ENCODING")            
            NbRecords(i1)=1
            SIZE_GEST(i1)=SIZE_INT
            SIZE_VARIABLES(i1)=0_8
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%BLRARRAY_ENCODING)) THEN
                  CALL DMUMPS_SAVE_RESTORE_BLR(id%BLRARRAY_ENCODING
     &                 ,unit,id%MYID,"memory_save"
     &                 ,SIZE_GEST_BLR,SIZE_VARIABLES_BLR
     &                 ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                 ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))                 
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%BLRARRAY_ENCODING)) THEN
                  write(unit,iostat=err) size(id%BLRARRAY_ENCODING,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  CALL DMUMPS_SAVE_RESTORE_BLR(id%BLRARRAY_ENCODING
     &                 ,unit,id%MYID,"save"
     &                 ,SIZE_GEST_BLR,SIZE_VARIABLES_BLR
     &                 ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                 ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%BLRARRAY_ENCODING)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.NE.-999) then
                  CALL DMUMPS_SAVE_RESTORE_BLR(id%BLRARRAY_ENCODING
     &                 ,unit,id%MYID,"restore"
     &                 ,SIZE_GEST_BLR,SIZE_VARIABLES_BLR
     &                 ,SIZE_INT, SIZE_ARITH_DEP, SIZE_LOGICAL
     &                 ,TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
     &                 ,size_read,size_allocated,size_written
     &                 ,id%INFO(1))   
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("L0_OMP_FACTORS")
         CASE("SCHED_SBTR")
         CASE("LPOOL_A_L0_OMP")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%LPOOL_A_L0_OMP
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%LPOOL_A_L0_OMP
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("LPOOL_B_L0_OMP")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%LPOOL_B_L0_OMP
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%LPOOL_B_L0_OMP
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("L_PHYS_L0_OMP")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%L_PHYS_L0_OMP
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%L_PHYS_L0_OMP
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("L_VIRT_L0_OMP")     
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%L_VIRT_L0_OMP
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%L_VIRT_L0_OMP
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("LL0_OMP_MAPPING")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%LL0_OMP_MAPPING
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%LL0_OMP_MAPPING
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("LL0_OMP_FACTORS")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%LL0_OMP_FACTORS
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT
               read(unit,iostat=err) id%LL0_OMP_FACTORS
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("THREAD_LA")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT8
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%THREAD_LA
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT8
               read(unit,iostat=err) id%THREAD_LA
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("IPOOL_A_L0_OMP")
         CASE("IPOOL_B_L0_OMP")
         CASE("PHYS_L0_OMP")
         CASE("VIRT_L0_OMP")
         CASE("VIRT_L0_OMP_MAPPING")
         CASE("PERM_L0_OMP")
         CASE("PTR_LEAFS_L0_OMP")
         CASE("L0_OMP_MAPPING")
         CASE("SINGULAR_VALUES")
        CASE("NB_SINGULAR_VALUES")
           NbRecords(i1)=1
           if(trim(mode).EQ."memory_save") then
              SIZE_VARIABLES(i1)=SIZE_INT
           elseif(trim(mode).EQ."save") then
              write(unit,iostat=err) id%NB_SINGULAR_VALUES
              if(err.ne.0) then
                 id%INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                id%INFO(2))
              endif
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           elseif(trim(mode).EQ."restore") then
              SIZE_VARIABLES(i1)=SIZE_INT
              read(unit,iostat=err) id%NB_SINGULAR_VALUES
              if(err.ne.0) THEN
                 id%INFO(1) = -75
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                ,id%INFO(2))
              endif 
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           endif
        CASE("ASSOCIATED_OOC_FILES")
            if(trim(mode).EQ."memory_save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_LOGICAL
            elseif(trim(mode).EQ."save") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_LOGICAL
            elseif(trim(mode).EQ."restore") then
               DIFF_SIZE_ALLOC_READ(i1)=SIZE_LOGICAL
            endif
        CASE("SAVE_DIR")
           NbRecords(i1)=1
           if(trim(mode).EQ."memory_save") then
              SIZE_VARIABLES(i1)=len(id%SAVE_DIR)*SIZE_CHARACTER
           elseif(trim(mode).EQ."save") then
              write(unit,iostat=err) id%SAVE_DIR
              if(err.ne.0) then
                 id%INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                id%INFO(2))
              endif
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           elseif(trim(mode).EQ."restore") then
              SIZE_VARIABLES(i1)=len(id%SAVE_DIR)*SIZE_CHARACTER
              read(unit,iostat=err) id%SAVE_DIR
              if(err.ne.0) THEN
                 id%INFO(1) = -75
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                ,id%INFO(2))
              endif 
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           endif
        CASE("SAVE_PREFIX")
           NbRecords(i1)=1
           if(trim(mode).EQ."memory_save") then
              SIZE_VARIABLES(i1)=len(id%SAVE_PREFIX)*SIZE_CHARACTER
           elseif(trim(mode).EQ."save") then
              write(unit,iostat=err) id%SAVE_PREFIX
              if(err.ne.0) then
                 id%INFO(1) = -72
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                id%INFO(2))
              endif
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           elseif(trim(mode).EQ."restore") then
              SIZE_VARIABLES(i1)=len(id%SAVE_PREFIX)*SIZE_CHARACTER
              read(unit,iostat=err) id%SAVE_PREFIX
              if(err.ne.0) THEN
                 id%INFO(1) = -75
                 CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                ,id%INFO(2))
              endif 
              CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID )
              IF ( id%INFO(1) .LT. 0 ) GOTO 100
           endif
        CASE("MPITOOMP_PROCS_MAP")
           NbRecords(i1)=2
            if(trim(mode).EQ."memory_save") then
               IF(associated(id%MPITOOMP_PROCS_MAP)) THEN
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=
     &                 size(id%MPITOOMP_PROCS_MAP,1)*SIZE_INT
               ELSE
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
               ENDIF
            elseif(trim(mode).EQ."save") then
               IF(associated(id%MPITOOMP_PROCS_MAP)) THEN
                  write(unit,iostat=err) size(id%MPITOOMP_PROCS_MAP,1)
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) id%MPITOOMP_PROCS_MAP
               ELSE
                  write(unit,iostat=err) -999
                  if(err.ne.0) then
                     id%INFO(1) = -72
                     CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                    id%INFO(2))
                  endif
                  CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                 id%COMM, id%MYID )
                  IF ( id%INFO(1) .LT. 0 ) GOTO 100
                  write(unit,iostat=err) -999
               ENDIF 
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               nullify(id%MPITOOMP_PROCS_MAP)
               read(unit,iostat=err) size_array1
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(size_array1.EQ.-999) then
                  SIZE_GEST(i1)=SIZE_INT*2
                  SIZE_VARIABLES(i1)=0_8
                  read(unit,iostat=err) dummy
               else
                  SIZE_GEST(i1)=SIZE_INT
                  SIZE_VARIABLES(i1)=size_array1*SIZE_INT
                  allocate(id%MPITOOMP_PROCS_MAP(size_array1),
     &                 stat=allocok)
                  if (allocok .GT. 0) THEN
                     id%INFO(1) = -78
                     CALL MUMPS_SETI8TOI4(
     &                    TOTAL_STRUC_SIZE-size_allocated
     &                    ,id%INFO(2))
                  endif
                  read(unit,iostat=err) id%MPITOOMP_PROCS_MAP
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif
         CASE("METIS_OPTIONS")
            NbRecords(i1)=1
            if(trim(mode).EQ."memory_save") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%METIS_OPTIONS,1)
            elseif(trim(mode).EQ."save") then
               write(unit,iostat=err) id%METIS_OPTIONS
               if(err.ne.0) then
                  id%INFO(1) = -72
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_written,
     &                 id%INFO(2))
               endif
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            elseif(trim(mode).EQ."restore") then
               SIZE_VARIABLES(i1)=SIZE_INT*size(id%METIS_OPTIONS,1)
               read(unit,iostat=err) id%METIS_OPTIONS
               if(err.ne.0) THEN
                  id%INFO(1) = -75
                  CALL MUMPS_SETI8TOI4(TOTAL_FILE_SIZE-size_read
     &                 ,id%INFO(2))
               endif 
               CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &              id%COMM, id%MYID )
               IF ( id%INFO(1) .LT. 0 ) GOTO 100
            endif      
      CASE("pad0","pad1","pad2","pad3","pad4","pad5","pad6","pad7",
     &     "pad11","pad111", "pad12","pad13","pad14","pad15","pad16")
      CASE DEFAULT
      END SELECT
      if(trim(mode).EQ."memory_save") then
         NbSubRecords=int(SIZE_VARIABLES(i1)/huge(0))
         IF(NbSubRecords.GT.0) then
            NbRecords(i1)=NbRecords(i1)+NbSubRecords
         ENDIF
      elseif(trim(mode).EQ."save") then
         size_written=size_written+SIZE_VARIABLES(i1)
     &        +int(SIZE_GEST(i1),kind=8)
#if !defined(MUMPS_F2003)
         size_written=size_written
     &        +int(2*id%KEEP(34)*NbRecords(i1),kind=8)
#endif
      elseif(trim(mode).EQ."restore") then
         size_allocated=size_allocated+SIZE_VARIABLES(i1)+
     &        DIFF_SIZE_ALLOC_READ(i1)
         size_read=size_read+SIZE_VARIABLES(i1)
     &        +int(SIZE_GEST(i1),kind=8)
#if !defined(MUMPS_F2003)
         size_read=size_read
     &        +int(2*id%KEEP(34)*NbRecords(i1),kind=8)
#endif
      elseif(trim(mode).EQ."fake_restore") then
      endif
      ENDDO
 200  continue
      if(trim(mode).EQ."memory_save") then
         WRITTEN_STRUC_SIZE=sum(SIZE_VARIABLES)+sum(SIZE_VARIABLES_ROOT)
     &        +SIZE_VARIABLES_BLR+SIZE_VARIABLES_FRONT_DATA+
     &         SIZE_VARIABLES_L0FAC
         TOTAL_STRUC_SIZE=WRITTEN_STRUC_SIZE
     &        +sum(DIFF_SIZE_ALLOC_READ)
     &        +sum(DIFF_SIZE_ALLOC_READ_ROOT)
         gest_size=sum(SIZE_GEST)+sum(SIZE_GEST_ROOT)
     &        +SIZE_GEST_BLR+SIZE_GEST_FRONT_DATA
     &        +SIZE_GEST_L0FAC
     &        +int(5*SIZE_CHARACTER,kind=8)
     &        +int(23*SIZE_CHARACTER,kind=8)
     &        +int(2*SIZE_INT8,kind=8)+int(1,kind=8)
     &        +int(3*SIZE_INT,kind=8)  
     &        +int(SIZE_LOGICAL,kind=8)
         IF(associated(id%OOC_FILE_NAME_LENGTH).AND.
     &        associated(id%OOC_FILE_NAMES)) THEN
            gest_size=gest_size+int(SIZE_INT,kind=8)
     &           +int(id%OOC_FILE_NAME_LENGTH(1)*SIZE_CHARACTER,kind=8)
         ELSE
            gest_size=gest_size+int(2*SIZE_INT,kind=8)
         ENDIF
#if !defined(MUMPS_F2003)
         tot_NbRecords=sum(NbRecords)+sum(NbRecords_ROOT)+8
         gest_size=gest_size+int(2*id%KEEP(34)*tot_NbRecords,kind=8)
#endif
         TOTAL_FILE_SIZE=WRITTEN_STRUC_SIZE+gest_size
      elseif(trim(mode).EQ."save") then
      elseif(trim(mode).EQ."restore") then
         if(id%root%gridinit_done) then
            id%root%CNTXT_BLACS = id%COMM_NODES
            CALL blacs_gridinit( id%root%CNTXT_BLACS, 'R',
     &           id%root%NPROW, id%root%NPCOL )
            id%root%gridinit_done = .TRUE.
         endif
      elseif(trim(mode).EQ."fake_restore") then
      elseif(trim(mode).EQ."restore_ooc") then
      endif
 100  continue
      deallocate(VARIABLES, VARIABLES_ROOT)
      RETURN
      END SUBROUTINE DMUMPS_SAVE_RESTORE_STRUCTURE
      END MODULE DMUMPS_SAVE_RESTORE
