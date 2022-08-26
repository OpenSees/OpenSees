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
      MODULE DMUMPS_SAVE_RESTORE_FILES
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER :: LEN_SAVE_FILE
      PARAMETER( LEN_SAVE_FILE = 550)
      CONTAINS
      SUBROUTINE MUMPS_READ_HEADER(fileunit, ierr, size_read, SIZE_INT
     &     ,SIZE_INT8, TOTAL_FILE_SIZE, TOTAL_STRUC_SIZE
     &     ,READ_ARITH, READ_INT_TYPE_64
     &     ,READ_OOC_FILE_NAME_LENGTH, READ_OOC_FIRST_FILE_NAME
     &     ,READ_HASH,READ_SYM,READ_PAR,READ_NPROCS
     &     ,FORTRAN_VERSION_OK)
      INTEGER,intent(in) :: fileunit
      INTEGER,intent(out) :: ierr
      INTEGER(8), intent(inout) :: size_read
      INTEGER,intent(in) :: SIZE_INT, SIZE_INT8
      INTEGER(8), intent(out) :: TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      CHARACTER, intent(out)  :: READ_ARITH
      LOGICAL, intent(out)    :: READ_INT_TYPE_64
      INTEGER, intent(out)    :: READ_OOC_FILE_NAME_LENGTH
      CHARACTER(len=LEN_SAVE_FILE),intent(out)::READ_OOC_FIRST_FILE_NAME
      CHARACTER(len=23), intent(out) :: READ_HASH
      INTEGER, intent(out)    :: READ_SYM,READ_PAR,READ_NPROCS
      LOGICAL, intent(out)    :: FORTRAN_VERSION_OK
      CHARACTER(len=5) :: READ_FORTRAN_VERSION
      INTEGER :: SIZE_CHARACTER, SIZE_LOGICAL
      INTEGER :: dummy
      SIZE_CHARACTER = 1
      SIZE_LOGICAL   = 4
      FORTRAN_VERSION_OK = .true.
      read(fileunit,iostat=ierr) READ_FORTRAN_VERSION
      if(ierr.ne.0) GOTO 100
      if (READ_FORTRAN_VERSION.NE."MUMPS") THEN
        ierr = 0
        FORTRAN_VERSION_OK = .false.
        GOTO 100 
      endif
      size_read=size_read+int(5*SIZE_CHARACTER,kind=8)
#if !defined(MUMPS_F2003)
      size_read=size_read
     &         +int(2*SIZE_INT*1,kind=8)
#endif
      read(fileunit,iostat=ierr) READ_HASH
      if(ierr.ne.0) GOTO 100
      size_read=size_read+int(23*SIZE_CHARACTER,kind=8)
#if !defined(MUMPS_F2003)
      size_read=size_read
     &         +int(2*SIZE_INT*1,kind=8)
#endif
      read(fileunit,iostat=ierr) TOTAL_FILE_SIZE,TOTAL_STRUC_SIZE
      if(ierr.ne.0) GOTO 100
      size_read=size_read+int(2*SIZE_INT8,kind=8)
#if !defined(MUMPS_F2003)
      size_read=size_read
     &         +int(2*SIZE_INT*1,kind=8)
#endif
      read(fileunit,iostat=ierr) READ_ARITH
      if(ierr.ne.0) GOTO 100
      size_read=size_read+int(1,kind=8)
#if !defined(MUMPS_F2003)
      size_read=size_read
     &         +int(2*SIZE_INT*1,kind=8)
#endif
      read(fileunit,iostat=ierr) READ_SYM,READ_PAR,READ_NPROCS
      if(ierr.ne.0) GOTO 100
      size_read=size_read+int(3*SIZE_INT,kind=8)
#if !defined(MUMPS_F2003)
      size_read=size_read
     &         +int(2*SIZE_INT*1,kind=8)
#endif     
      read(fileunit,iostat=ierr) READ_INT_TYPE_64
      if(ierr.ne.0) GOTO 100
      size_read=size_read+int(SIZE_LOGICAL,kind=8)
#if !defined(MUMPS_F2003)
      size_read=size_read
     &         +int(2*SIZE_INT*1,kind=8)
#endif
      read(fileunit,iostat=ierr) READ_OOC_FILE_NAME_LENGTH
      if(ierr.ne.0) GOTO 100
      size_read=size_read+int(SIZE_INT,kind=8)
#if !defined(MUMPS_F2003)
      size_read=size_read
     &         +int(2*SIZE_INT*1,kind=8)
#endif
      IF(READ_OOC_FILE_NAME_LENGTH.EQ.-999) THEN
        read(fileunit,iostat=ierr) dummy
        if(ierr.ne.0) GOTO 100
        size_read=size_read+int(SIZE_INT,kind=8)
#if !defined(MUMPS_F2003)
        size_read=size_read
     &           +int(2*SIZE_INT*1,kind=8)
#endif
      ELSE
        read(fileunit,iostat=ierr)
     &    READ_OOC_FIRST_FILE_NAME(1:READ_OOC_FILE_NAME_LENGTH)
        if(ierr.ne.0) GOTO 100
        size_read=size_read+int(
     &           READ_OOC_FILE_NAME_LENGTH*SIZE_CHARACTER,kind=8)
#if !defined(MUMPS_F2003)
        size_read=size_read
     &           +int(2*SIZE_INT*1,kind=8)
#endif
#if defined(OOC_VERBOSE)
        write(*,*) 'First ooc file: ',
     &   READ_OOC_FIRST_FILE_NAME(1:READ_OOC_FILE_NAME_LENGTH-2)
#endif
      ENDIF
 100  continue
      RETURN
      END SUBROUTINE MUMPS_READ_HEADER
      SUBROUTINE DMUMPS_CHECK_HEADER(id, BASIC_CHECK, READ_INT_TYPE_64,
     &                READ_HASH, READ_NPROCS,
     &                READ_ARITH, READ_SYM, READ_PAR)
      INCLUDE 'mpif.h'
      TYPE (DMUMPS_STRUC),intent(inout) :: id
      LOGICAL, intent(in)           :: BASIC_CHECK
      LOGICAL, intent(in)           :: READ_INT_TYPE_64
      CHARACTER(len=23), intent(in) :: READ_HASH
      INTEGER, intent(in)           :: READ_NPROCS
      CHARACTER, intent(in)         :: READ_ARITH
      INTEGER, intent(in)           :: READ_SYM,READ_PAR
      LOGICAL   :: INT_TYPE_64
      CHARACTER(len=23) :: HASH_MASTER
      CHARACTER :: ARITH
      INTEGER   :: IERR
      IF(id%KEEP(10).EQ.1) THEN
         INT_TYPE_64=.TRUE.
      ELSE
         INT_TYPE_64=.FALSE.
      ENDIF
      if(INT_TYPE_64.neqv.READ_INT_TYPE_64) THEN
         id%INFO(1) = -73
         id%INFO(2) = 2
      endif 
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 100
      if(id%MYID.EQ.0) THEN
         HASH_MASTER=READ_HASH
      ENDIF
      call MPI_BCAST(HASH_MASTER,23,MPI_CHARACTER,0,id%COMM,IERR)
      if(HASH_MASTER.ne.READ_HASH) THEN
         id%INFO(1) = -73
         id%INFO(2) = 3
      endif 
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 100
      if(id%NPROCS.ne.READ_NPROCS) THEN
         id%INFO(1) = -73
         id%INFO(2) = 4
      endif 
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 100  
      IF (.NOT.BASIC_CHECK) THEN
        ARITH="DMUMPS"(1:1)
        if(ARITH.ne.READ_ARITH) THEN
           id%INFO(1) = -73
           id%INFO(2) = 5
        endif 
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &       id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GOTO 100
        if((id%MYID.EQ.0).AND.(id%SYM.ne.READ_SYM)) THEN
           id%INFO(1) = -73
           id%INFO(2) = 6
        endif
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &       id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GOTO 100
        if((id%MYID.EQ.0).AND.(id%PAR.ne.READ_PAR)) THEN
           write (*,*) id%MYID, 'PAR ',id%PAR, 'READ_PAR ', READ_PAR
           id%INFO(1) = -73
           id%INFO(2) = 7
        endif
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &       id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GOTO 100
      ENDIF
 100  continue
      RETURN
      END SUBROUTINE DMUMPS_CHECK_HEADER
      SUBROUTINE MUMPS_CLEAN_SAVED_DATA(MYID,ierr,SUPPFILE,INFOFILE)
      INCLUDE 'mpif.h'
      INTEGER,intent(in)  :: MYID
      INTEGER,intent(out) :: ierr
      CHARACTER(len=LEN_SAVE_FILE),intent(in):: SUPPFILE,INFOFILE
      INTEGER::supp,tmp_err
      ierr = 0
      tmp_err = 0
      supp=200+MYID
      open(UNIT=supp,FILE=SUPPFILE,STATUS='old',
     &     form='unformatted',iostat=tmp_err)
      if (tmp_err.eq.0) THEN
        close(UNIT=supp,STATUS='delete',iostat=tmp_err)
        if(tmp_err.ne.0) then
          ierr = 1
          tmp_err = 0
        endif
      endif
      if (ierr .eq. 0) then
        if (tmp_err.ne.0) then
          ierr = 1
          tmp_err = 0
        endif
        open(UNIT=supp,FILE=INFOFILE,STATUS='old',iostat=tmp_err)
        if (tmp_err.eq.0) THEN
          close(UNIT=supp,STATUS='delete',iostat=tmp_err)
        endif
        if (tmp_err.ne.0) THEN
          ierr = ierr + 2
          tmp_err = 0
        endif
      endif
      END SUBROUTINE MUMPS_CLEAN_SAVED_DATA
      SUBROUTINE DMUMPS_GET_SAVE_FILES(id,SAVE_FILE,INFO_FILE)
      INCLUDE 'mpif.h'
      TYPE (DMUMPS_STRUC),intent(inout) :: id
      CHARACTER(len=LEN_SAVE_FILE),intent(out):: SAVE_FILE, INFO_FILE
      INTEGER::len_save_dir,len_save_prefix
      CHARACTER(len=255):: tmp_savedir,savedir
      CHARACTER(len=255):: tmp_saveprefix,saveprefix
      CHARACTER(len=10):: STRING_MYID
      CHARACTER:: LAST_CHAR_DIR
      INFO_FILE=''
      SAVE_FILE=''
      tmp_savedir=''
      tmp_saveprefix=''
      IF(id%SAVE_DIR.EQ."NAME_NOT_INITIALIZED") THEN
         call mumps_get_save_dir_C(len_save_dir,tmp_savedir)
         if(tmp_savedir(1:len_save_dir).EQ."NAME_NOT_INITIALIZED") then
            id%INFO(1) = -77
            id%INFO(2) = 0 
         else
            savedir=trim(adjustl(tmp_savedir(1:len_save_dir)))
            len_save_dir=len_trim(savedir(1:len_save_dir))
         endif
      ELSE
         savedir=trim(adjustl(id%SAVE_DIR))
         len_save_dir=len_trim(savedir)
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 100
      IF(id%SAVE_PREFIX.EQ."NAME_NOT_INITIALIZED") THEN
         call mumps_get_save_prefix_C(len_save_prefix,tmp_saveprefix)
         if(tmp_saveprefix(1:len_save_prefix).EQ."NAME_NOT_INITIALIZED")
     &        then
            saveprefix="save"
            len_save_prefix=len_trim(saveprefix)
         else
            saveprefix=
     &           trim(adjustl(tmp_saveprefix(1:len_save_prefix)))
            len_save_prefix=len_trim(saveprefix(1:len_save_prefix))
         endif
      ELSE
         saveprefix=trim(adjustl(id%SAVE_PREFIX))
         len_save_prefix=len_trim(saveprefix)
      ENDIF
      write(STRING_MYID,'(I10)') id%MYID
      LAST_CHAR_DIR=savedir(len_save_dir:len_save_dir)
      if(LAST_CHAR_DIR.NE."/") then
         SAVE_FILE=trim(adjustl(savedir))//"/"
      else
         SAVE_FILE=trim(adjustl(savedir))
      endif
      INFO_FILE=trim(adjustl(SAVE_FILE))
      SAVE_FILE=trim(adjustl(SAVE_FILE))
     &     //trim(adjustl(saveprefix))
     &     //"_"
     &     //trim(adjustl(STRING_MYID))
     &     //".mumps"
      INFO_FILE=trim(adjustl(INFO_FILE))
     &        //trim(adjustl(saveprefix))
     &        //"_"
     &        //trim(adjustl(STRING_MYID))
     &        //".info"
 100  continue
      RETURN
      END SUBROUTINE DMUMPS_GET_SAVE_FILES
      SUBROUTINE DMUMPS_CHECK_FILE_NAME(id,NAME_LENGTH,FILE_NAME,CHECK)
      TYPE (DMUMPS_STRUC),intent(in) :: id
      INTEGER,intent(in)             :: NAME_LENGTH
      CHARACTER(len=LEN_SAVE_FILE),intent(in)   :: FILE_NAME
      LOGICAL,intent(out)            :: CHECK
      INTEGER :: I
      CHECK = .false.
      IF (NAME_LENGTH.NE.-999) THEN
        IF (associated(id%OOC_FILE_NAME_LENGTH) .AND.
     &      associated(id%OOC_FILE_NAMES)) THEN
          IF (NAME_LENGTH .EQ. id%OOC_FILE_NAME_LENGTH(1)) THEN
            CHECK = .true.
            I = 1
            DO WHILE(I.LE.NAME_LENGTH)
              IF (FILE_NAME(I:I).NE.id%OOC_FILE_NAMES(1,I)) THEN
                CHECK = .false.
                I = NAME_LENGTH + 1
              ELSE
                I = I + 1
              ENDIF
            END DO
          ENDIF
        ENDIF
      ENDIF
      END SUBROUTINE DMUMPS_CHECK_FILE_NAME      
      END MODULE DMUMPS_SAVE_RESTORE_FILES
      SUBROUTINE DMUMPS_SAVE_FILES_RETURN()
      RETURN
      END SUBROUTINE DMUMPS_SAVE_FILES_RETURN
