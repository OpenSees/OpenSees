!===============================================================================
! Copyright 2014-2018 Intel Corporation All Rights Reserved.
!
! The source code,  information  and material  ("Material") contained  herein is
! owned by Intel Corporation or its  suppliers or licensors,  and  title to such
! Material remains with Intel  Corporation or its  suppliers or  licensors.  The
! Material  contains  proprietary  information  of  Intel or  its suppliers  and
! licensors.  The Material is protected by  worldwide copyright  laws and treaty
! provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
! modified, published,  uploaded, posted, transmitted,  distributed or disclosed
! in any way without Intel's prior express written permission.  No license under
! any patent,  copyright or other  intellectual property rights  in the Material
! is granted to  or  conferred  upon  you,  either   expressly,  by implication,
! inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
! property rights must be express and approved by Intel in writing.
!
! Unless otherwise agreed by Intel in writing,  you may not remove or alter this
! notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
! suppliers or licensors in any way.
!===============================================================================

!
!   This file contains preliminary version of new SpBLAS API which supports
!   two-step execution (inspector-executor) model.
!
!*******************************************************************************

!*******************************************************************************
!*********************************** Basic types and constants *****************
!*******************************************************************************

    MODULE MKL_SPBLAS

    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INTPTR_T, C_INT

!   status of the routines
    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_STATUS_SUCCESS               = 0,  &  !
                      SPARSE_STATUS_NOT_INITIALIZED       = 1,  &  ! empty handle or matrix arrays
                      SPARSE_STATUS_ALLOC_FAILED          = 2,  &  ! internal error: memory allocation failed
                      SPARSE_STATUS_INVALID_VALUE         = 3,  &  ! invalid input value
                      SPARSE_STATUS_EXECUTION_FAILED      = 4,  &  ! e.g. 0-diagonal element for triangular solver, etc
                      SPARSE_STATUS_INTERNAL_ERROR        = 5,  &  !
                      SPARSE_STATUS_NOT_SUPPORTED         = 6      ! e.g. operation for double precision don't support other types
    END ENUM

!   sparse matrix operations
    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_OPERATION_NON_TRANSPOSE      = 10, &  !
                      SPARSE_OPERATION_TRANSPOSE          = 11, &  !
                      SPARSE_OPERATION_CONJUGATE_TRANSPOSE= 12     !
    END ENUM

!   supported matrix types
    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_MATRIX_TYPE_GENERAL          = 20, &  ! general case
                      SPARSE_MATRIX_TYPE_SYMMETRIC        = 21, &  ! triangular part of
                      SPARSE_MATRIX_TYPE_HERMITIAN        = 22, &  ! the matrix is to be processed
                      SPARSE_MATRIX_TYPE_TRIANGULAR       = 23, &  !
                      SPARSE_MATRIX_TYPE_DIAGONAL         = 24     ! diagonal matrix; only diagonal elements should be processed
    END ENUM

!   sparse matrix indexing: C-style or Fortran-style
    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_INDEX_BASE_ZERO              = 0,  &  ! C-style
                      SPARSE_INDEX_BASE_ONE               = 1      ! Fortran-style
    END ENUM

!   makes sense for triangular matrices only
!       ( SPARSE_MATRIX_TYPE_SYMMETRIC, SPARSE_MATRIX_TYPE_HERMITIAN, SPARSE_MATRIX_TYPE_TRIANGULAR )
    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_FILL_MODE_LOWER              = 40, &  ! lower triangular part of the matrix is stored
                      SPARSE_FILL_MODE_UPPER              = 41, &  ! upper triangular part of the matrix is stored
                      SPARSE_FILL_MODE_FULL              = 42     ! upper triangular part of the matrix is stored
    END ENUM

!   makes sense for triangular matrices only
!       ( SPARSE_MATRIX_TYPE_SYMMETRIC, SPARSE_MATRIX_TYPE_HERMITIAN, SPARSE_MATRIX_TYPE_TRIANGULAR )
    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_DIAG_NON_UNIT                = 50, &  ! triangular matrix with non-unit diagonal
                      SPARSE_DIAG_UNIT                    = 51     ! triangular matrix with unit     diagonal
    END ENUM

!   applicable for Level 3 operations with dense matrices; describes storage scheme for dense matrix (row major or column major)
    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_LAYOUT_ROW_MAJOR             = 101, & ! C-style
                      SPARSE_LAYOUT_COLUMN_MAJOR          = 102    ! Fortran-style
    END ENUM

!   verbose mode; if verbose mode activated, handle should collect and report profiling / optimization info
    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_VERBOSE_OFF                  = 70, &  !
                      SPARSE_VERBOSE_BASIC                = 71, &  ! report high-level info about optimization algorithms, issues, etc.
                      SPARSE_VERBOSE_EXTENDED             = 72     ! provide detailed information
    END ENUM

!   memory optimization hints from customer: describe how much memory could be used on optimization stage
    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_MEMORY_NONE                  = 80, &  ! no memory should be allocated for matrix values and structures
                                                                   ! auxiliary structures could be created only for workload balancing,
                                                                   ! parallelization, etc.
                      SPARSE_MEMORY_AGGRESSIVE            = 81     ! matrix could be converted to any internal format
    END ENUM

    ENUM, BIND(C)
        ENUMERATOR :: SPARSE_STAGE_FULL_MULT               = 90, & !
                      SPARSE_STAGE_NNZ_COUNT               = 91, & !
                      SPARSE_STAGE_FINALIZE_MULT           = 92
    END ENUM

!*************************************************************************************************
!*** Opaque structure for sparse matrix in internal format, further D - means double precision ***
!*************************************************************************************************

!    struct  sparse_matrix;
!    typedef struct sparse_matrix *sparse_matrix_t;

    TYPE, BIND(C) :: SPARSE_MATRIX_T
        INTEGER(C_INTPTR_T) :: PTR
    END TYPE SPARSE_MATRIX_T

!   descriptor of main sparse matrix properties
    TYPE, BIND(C) :: MATRIX_DESCR
        INTEGER(C_INT) :: TYPE
        INTEGER(C_INT) :: MODE
        INTEGER(C_INT) :: DIAG
    END TYPE MATRIX_DESCR

!*****************************************************************************************
!*************************************** Creation routines *******************************
!*****************************************************************************************


    INTERFACE

!   Matrix handle is used for storing information about the matrix and matrix values

!   Create matrix from one of the existing sparse formats, just create the handle with matrix info and copy matrix values if
!   requested
!   Collect high-level info about the matrix. Need to use this interface for the case with several calls in program for
!   performance reasons, where optimizations are not required


!   coordinate format,
!   SPARSE_MATRIX_TYPE_GENERAL by default, pointers to input arrays are stored in the handle
        FUNCTION MKL_SPARSE_S_CREATE_COO(A,indexing,rows,cols,nnz,row_indx,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_S_CREATE_COO')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN) :: nnz
            INTEGER, INTENT(IN), DIMENSION(*) :: row_indx
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_S_CREATE_COO
        END FUNCTION
        FUNCTION MKL_SPARSE_D_CREATE_COO(A,indexing,rows,cols,nnz,row_indx,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_D_CREATE_COO')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN) :: nnz
            INTEGER, INTENT(IN), DIMENSION(*) :: row_indx
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_D_CREATE_COO
        END FUNCTION
        FUNCTION MKL_SPARSE_C_CREATE_COO(A,indexing,rows,cols,nnz,row_indx,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_C_CREATE_COO')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN) :: nnz
            INTEGER, INTENT(IN), DIMENSION(*) :: row_indx
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_C_CREATE_COO
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_CREATE_COO(A,indexing,rows,cols,nnz,row_indx,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_Z_CREATE_COO')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN) :: nnz
            INTEGER, INTENT(IN), DIMENSION(*) :: row_indx
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_Z_CREATE_COO
        END FUNCTION

!   compressed sparse row format (4-arrays version),
!   SPARSE_MATRIX_TYPE_GENERAL by default, pointers to input arrays are stored in the handle
        FUNCTION MKL_SPARSE_S_CREATE_CSR(A,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_S_CREATE_CSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_S_CREATE_CSR
        END FUNCTION
        FUNCTION MKL_SPARSE_D_CREATE_CSR(A,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_D_CREATE_CSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_D_CREATE_CSR
        END FUNCTION
        FUNCTION MKL_SPARSE_C_CREATE_CSR(A,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_C_CREATE_CSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_C_CREATE_CSR
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_CREATE_CSR(A,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_Z_CREATE_CSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_Z_CREATE_CSR
        END FUNCTION

!   compressed sparse column format (4-arrays version),
!   SPARSE_MATRIX_TYPE_GENERAL by default, pointers to input arrays are stored in the handle
        FUNCTION MKL_SPARSE_S_CREATE_CSC(A,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_S_CREATE_CSC')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_S_CREATE_CSC
        END FUNCTION
        FUNCTION MKL_SPARSE_D_CREATE_CSC(A,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_D_CREATE_CSC')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_D_CREATE_CSC
        END FUNCTION
        FUNCTION MKL_SPARSE_C_CREATE_CSC(A,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_C_CREATE_CSC')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_C_CREATE_CSC
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_CREATE_CSC(A,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_Z_CREATE_CSC')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_Z_CREATE_CSC
        END FUNCTION

!   compressed block sparse row format (4-arrays version, square blocks),
!   SPARSE_MATRIX_TYPE_GENERAL by default, pointers to input arrays are stored in the handle
        FUNCTION MKL_SPARSE_S_CREATE_BSR(A,indexing,block_layout,rows,cols,block_size,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_S_CREATE_BSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER(C_INT), INTENT(IN) :: block_layout
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN) :: block_size
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_S_CREATE_BSR
        END FUNCTION
        FUNCTION MKL_SPARSE_D_CREATE_BSR(A,indexing,block_layout,rows,cols,block_size,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_D_CREATE_BSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER(C_INT), INTENT(IN) :: block_layout
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN) :: block_size
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_D_CREATE_BSR
        END FUNCTION
        FUNCTION MKL_SPARSE_C_CREATE_BSR(A,indexing,block_layout,rows,cols,block_size,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_C_CREATE_BSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER(C_INT), INTENT(IN) :: block_layout
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN) :: block_size
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_C_CREATE_BSR
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_CREATE_BSR(A,indexing,block_layout,rows,cols,block_size,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_Z_CREATE_BSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT), INTENT(IN) :: indexing
            INTEGER(C_INT), INTENT(IN) :: block_layout
            INTEGER, INTENT(IN) :: rows
            INTEGER, INTENT(IN) :: cols
            INTEGER, INTENT(IN) :: block_size
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_start
            INTEGER, INTENT(IN), DIMENSION(*) :: rows_end
            INTEGER, INTENT(IN), DIMENSION(*) :: col_indx
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_Z_CREATE_BSR
        END FUNCTION

!   Create copy of the existing handle; matrix properties could be changed.
!   For example it could be used for extracting triangular or diagonal parts from existing matrix.
        FUNCTION MKL_SPARSE_COPY(source,descr,dest) &
                 BIND(C, name='MKL_SPARSE_COPY')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: dest
            INTEGER(C_INT) MKL_SPARSE_COPY
        END FUNCTION

!   destroy matrix handle; if sparse matrix was stored inside the handle it also deallocates the matrix
!   It is customer responsibility not to delete the handle with the matrix, if this matrix is shared with other handles
        FUNCTION MKL_SPARSE_DESTROY(A) &
                 BIND(C, name='MKL_SPARSE_DESTROY')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT) MKL_SPARSE_DESTROY
        END FUNCTION

!   return extended error information from last operation;
!   eg. info about wrong input parameter, memory sizes that couldn't be allocated, etc.
        FUNCTION MKL_SPARSE_GET_ERROR_INFO(A,info) &
                 BIND(C, name='MKL_SPARSE_GET_ERROR_INFO')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER, INTENT(IN), DIMENSION(*) :: info
            INTEGER(C_INT) MKL_SPARSE_GET_ERROR_INFO
        END FUNCTION


!****************************************************************************************
!************************ Converters of internal representation  ************************
!****************************************************************************************

!   converters from current format to another

        ! convert original matrix to CSR representation
        FUNCTION MKL_SPARSE_CONVERT_CSR(source,operation,dest) &
                 BIND(C, name='MKL_SPARSE_CONVERT_CSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            INTEGER(C_INT)       , INTENT(IN) :: operation     ! as is, transposed or conjugate transposed
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: dest
            INTEGER(C_INT) MKL_SPARSE_CONVERT_CSR
        END FUNCTION

        ! convert original matrix to BSR representation
        FUNCTION MKL_SPARSE_CONVERT_BSR(source,block_size,block_layout,operation,dest) &
                 BIND(C, name='MKL_SPARSE_CONVERT_BSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            INTEGER        , INTENT(IN) :: block_size
            INTEGER(C_INT)       , INTENT(IN) :: block_layout    ! block storage: row-major or column-major
            INTEGER(C_INT)       , INTENT(IN) :: operation     ! as is, transposed or conjugate transposed
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: dest
            INTEGER(C_INT) MKL_SPARSE_CONVERT_BSR
        END FUNCTION

        FUNCTION MKL_SPARSE_S_EXPORT_BSR(source,indexing,block_layout,rows,cols,block_size,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_S_EXPORT_BSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            INTEGER(C_INT)       , INTENT(INOUT) :: indexing
            INTEGER(C_INT)       , INTENT(INOUT) :: block_layout
            INTEGER(C_INT)       , INTENT(INOUT) :: rows
            INTEGER(C_INT)       , INTENT(INOUT) :: cols
            INTEGER(C_INT)       , INTENT(INOUT) :: block_size
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_S_EXPORT_BSR
        END FUNCTION
        FUNCTION MKL_SPARSE_D_EXPORT_BSR(source,indexing,block_layout,rows,cols,block_size,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_D_EXPORT_BSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T) , INTENT(IN) :: source
            INTEGER(C_INT)        , INTENT(INOUT) :: indexing
            INTEGER(C_INT)        , INTENT(INOUT) :: block_layout
            INTEGER(C_INT)        , INTENT(INOUT) :: rows
            INTEGER(C_INT)        , INTENT(INOUT) :: cols
            INTEGER(C_INT)        , INTENT(INOUT) :: block_size
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_D_EXPORT_BSR
        END FUNCTION
        FUNCTION MKL_SPARSE_C_EXPORT_BSR(source,indexing,block_layout,rows,cols,block_size,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_C_EXPORT_BSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T) , INTENT(IN) :: source
            INTEGER(C_INT)        , INTENT(INOUT) :: indexing
            INTEGER(C_INT)        , INTENT(INOUT) :: block_layout
            INTEGER(C_INT)        , INTENT(INOUT) :: rows
            INTEGER(C_INT)        , INTENT(INOUT) :: cols
            INTEGER(C_INT)        , INTENT(INOUT) :: block_size
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_C_EXPORT_BSR
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_EXPORT_BSR(source,indexing,block_layout,rows,cols,block_size,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_Z_EXPORT_BSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T) , INTENT(IN) :: source
            INTEGER(C_INT)        , INTENT(INOUT) :: indexing
            INTEGER(C_INT)        , INTENT(INOUT) :: block_layout
            INTEGER(C_INT)        , INTENT(INOUT) :: rows
            INTEGER(C_INT)        , INTENT(INOUT) :: cols
            INTEGER(C_INT)        , INTENT(INOUT) :: block_size
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_Z_EXPORT_BSR
        END FUNCTION

        FUNCTION MKL_SPARSE_S_EXPORT_CSR(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_S_EXPORT_CSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            INTEGER(C_INT)      , INTENT(INOUT) :: indexing
            INTEGER(C_INT)      , INTENT(INOUT) :: rows
            INTEGER(C_INT)      , INTENT(INOUT) :: cols
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_S_EXPORT_CSR
        END FUNCTION
        FUNCTION MKL_SPARSE_D_EXPORT_CSR(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_D_EXPORT_CSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            INTEGER(C_INT)      , INTENT(INOUT) :: indexing
            INTEGER(C_INT)      , INTENT(INOUT) :: rows
            INTEGER(C_INT)      , INTENT(INOUT) :: cols
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_D_EXPORT_CSR
        END FUNCTION
        FUNCTION MKL_SPARSE_C_EXPORT_CSR(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_C_EXPORT_CSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            INTEGER(C_INT)      , INTENT(INOUT) :: indexing
            INTEGER(C_INT)      , INTENT(INOUT) :: rows
            INTEGER(C_INT)      , INTENT(INOUT) :: cols
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_C_EXPORT_CSR
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_EXPORT_CSR(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_Z_EXPORT_CSR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T) , INTENT(IN) :: source
            INTEGER(C_INT)        , INTENT(INOUT) :: indexing
            INTEGER(C_INT)        , INTENT(INOUT) :: rows
            INTEGER(C_INT)        , INTENT(INOUT) :: cols
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_Z_EXPORT_CSR
        END FUNCTION

        FUNCTION MKL_SPARSE_S_EXPORT_CSC(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_S_EXPORT_CSC')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            INTEGER(C_INT)      , INTENT(INOUT) :: indexing
            INTEGER(C_INT)      , INTENT(INOUT) :: rows
            INTEGER(C_INT)      , INTENT(INOUT) :: cols
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_S_EXPORT_CSC
        END FUNCTION
        FUNCTION MKL_SPARSE_D_EXPORT_CSC(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_D_EXPORT_CSC')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            INTEGER(C_INT)      , INTENT(INOUT) :: indexing
            INTEGER(C_INT)      , INTENT(INOUT) :: rows
            INTEGER(C_INT)      , INTENT(INOUT) :: cols
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_D_EXPORT_CSC
        END FUNCTION
        FUNCTION MKL_SPARSE_C_EXPORT_CSC(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_C_EXPORT_CSC')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: source
            INTEGER(C_INT)      , INTENT(INOUT) :: indexing
            INTEGER(C_INT)      , INTENT(INOUT) :: rows
            INTEGER(C_INT)      , INTENT(INOUT) :: cols
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_C_EXPORT_CSC
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_EXPORT_CSC(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
                 BIND(C, name='MKL_SPARSE_Z_EXPORT_CSC')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T) , INTENT(IN) :: source
            INTEGER(C_INT)        , INTENT(INOUT) :: indexing
            INTEGER(C_INT)        , INTENT(INOUT) :: rows
            INTEGER(C_INT)        , INTENT(INOUT) :: cols
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_start
            INTEGER, INTENT(INOUT), DIMENSION(*) :: rows_end
            INTEGER, INTENT(INOUT), DIMENSION(*) :: col_indx
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: values
            INTEGER(C_INT) MKL_SPARSE_Z_EXPORT_CSC
        END FUNCTION

!****************************************************************************************
!************************** Step-by-step modification routines **************************
!****************************************************************************************

!   update existing value in the matrix
!   ( for internal storage only, should not work with customer-allocated matrices)
        FUNCTION MKL_SPARSE_S_SET_VALUE(A,row,col,value) &
                 BIND(C, name='MKL_SPARSE_S_SET_VALUE')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER, INTENT(IN) :: row
            INTEGER, INTENT(IN) :: col
            REAL(C_FLOAT) , INTENT(IN) :: value
            INTEGER(C_INT) MKL_SPARSE_S_SET_VALUE
        END FUNCTION
        FUNCTION MKL_SPARSE_D_SET_VALUE(A,row,col,value) &
                 BIND(C, name='MKL_SPARSE_D_SET_VALUE')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER, INTENT(IN) :: row
            INTEGER, INTENT(IN) :: col
            REAL(C_DOUBLE) , INTENT(IN) :: value
            INTEGER(C_INT) MKL_SPARSE_D_SET_VALUE
        END FUNCTION
        FUNCTION MKL_SPARSE_C_SET_VALUE(A,row,col,value) &
                 BIND(C, name='MKL_SPARSE_C_SET_VALUE')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER, INTENT(IN) :: row
            INTEGER, INTENT(IN) :: col
            COMPLEX(C_FLOAT_COMPLEX) , INTENT(IN) :: value
            INTEGER(C_INT) MKL_SPARSE_C_SET_VALUE
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_SET_VALUE(A,row,col,value) &
                 BIND(C, name='MKL_SPARSE_Z_SET_VALUE')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER, INTENT(IN) :: row
            INTEGER, INTENT(IN) :: col
            COMPLEX(C_DOUBLE_COMPLEX) , INTENT(IN) :: value
            INTEGER(C_INT) MKL_SPARSE_Z_SET_VALUE
        END FUNCTION

!****************************************************************************************
!****************************** Verbose mode routine ************************************
!****************************************************************************************

!   allow to switch on/off verbose mode
        FUNCTION MKL_SPARSE_SET_VERBOSE_MODE(verbose) &
                 BIND(C, name='MKL_SPARSE_SET_VERBOSE_MODE')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            INTEGER(C_INT), INTENT(IN) :: verbose
            INTEGER(C_INT) MKL_SPARSE_SET_VERBOSE_MODE
        END FUNCTION

!****************************************************************************************
!****************************** Optimization routines ***********************************
!****************************************************************************************

!   Describe expected operations with amount of iterations
        FUNCTION MKL_SPARSE_SET_MV_HINT(A,operation,descr,expected_calls) &
                 BIND(C, name='MKL_SPARSE_SET_MV_HINT')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT)      , INTENT(IN) :: operation     ! SPARSE_OPERATION_NON_TRANSPOSE is default value for
                                                               ! infinite amount of calls
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr         ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER              , INTENT(IN) :: expected_calls
            INTEGER(C_INT) MKL_SPARSE_SET_MV_HINT
        END FUNCTION

        FUNCTION MKL_SPARSE_SET_MM_HINT(A,operation,descr,layout,dense_matrix_size,expected_calls) &
                 BIND(C, name='MKL_SPARSE_SET_MM_HINT')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr           ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)      , INTENT(IN) :: layout           ! storage scheme for the dense matrix: C-style or Fortran-style
            INTEGER              , INTENT(IN) :: dense_matrix_size ! amount of columns in dense matrix
            INTEGER              , INTENT(IN) :: expected_calls
            INTEGER(C_INT) MKL_SPARSE_SET_MM_HINT
        END FUNCTION

        FUNCTION MKL_SPARSE_SET_SV_HINT(A,operation,descr,expected_calls) &
                 BIND(C, name='MKL_SPARSE_SET_SV_HINT')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT)      , INTENT(IN) :: operation     ! SPARSE_OPERATION_NON_TRANSPOSE is default value for
                                                               ! infinite amount of calls
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr         ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER              , INTENT(IN) :: expected_calls
            INTEGER(C_INT) MKL_SPARSE_SET_SV_HINT
        END FUNCTION

        FUNCTION MKL_SPARSE_SET_SM_HINT(A,operation,descr,layout,dense_matrix_size,expected_calls) &
                 BIND(C, name='MKL_SPARSE_SET_SM_HINT')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr           ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)      , INTENT(IN) :: layout           ! storage scheme for the dense matrix: C-style or Fortran-style
            INTEGER              , INTENT(IN) :: dense_matrix_size ! amount of columns in dense matrix
            INTEGER              , INTENT(IN) :: expected_calls
            INTEGER(C_INT) MKL_SPARSE_SET_SM_HINT
        END FUNCTION

        FUNCTION MKL_SPARSE_SET_MEMORY_HINT(A,policy) &
                 BIND(C, name='MKL_SPARSE_SET_MEMORY_HINT')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT)      , INTENT(IN) :: policy          ! SPARSE_MEMORY_AGGRESSIVE is default value
            INTEGER(C_INT) MKL_SPARSE_SET_MEMORY_HINT
        END FUNCTION

!   optimize matrix described by the handle. It uses hints (optimization and memory) that should be set up before this call.
!   if hints were not explicitly defined, default vales are:
!   SPARSE_OPERATION_NON_TRANSPOSE for matrix-vector multiply with infinite amount of expected iterations
        FUNCTION MKL_SPARSE_OPTIMIZE(A) &
                 BIND(C, name='MKL_SPARSE_OPTIMIZE')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT) MKL_SPARSE_OPTIMIZE
        END FUNCTION

!****************************************************************************************
!****************************** Computational routines **********************************
!****************************************************************************************

!   Perform computations based on created matrix handle

!   Level 2

!   Computes y = alpha * A * x + beta * y
        FUNCTION MKL_SPARSE_S_MV(operation,alpha,A,descr,x,beta,y) &
                 BIND(C, name='MKL_SPARSE_S_MV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: x
            REAL(C_FLOAT)        , INTENT(IN) :: beta
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_S_MV
        END FUNCTION
        FUNCTION MKL_SPARSE_D_MV(operation,alpha,A,descr,x,beta,y) &
                 BIND(C, name='MKL_SPARSE_D_MV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            REAL(C_DOUBLE)       , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: x
            REAL(C_DOUBLE)       , INTENT(IN) :: beta
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_D_MV
        END FUNCTION
        FUNCTION MKL_SPARSE_C_MV(operation,alpha,A,descr,x,beta,y) &
                 BIND(C, name='MKL_SPARSE_C_MV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)          , INTENT(IN) :: operation
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T)   , INTENT(IN) :: A
            TYPE(MATRIX_DESCR)      , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: x
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN) :: beta
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_C_MV
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_MV(operation,alpha,A,descr,x,beta,y) &
                 BIND(C, name='MKL_SPARSE_Z_MV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: x
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: beta
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_Z_MV
        END FUNCTION

!   Solves triangular system y = alpha * A^{-1} * x
        FUNCTION MKL_SPARSE_S_TRSV(operation,alpha,A,descr,x,y) &
                 BIND(C, name='MKL_SPARSE_S_TRSV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: x
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_S_TRSV
        END FUNCTION
        FUNCTION MKL_SPARSE_D_TRSV(operation,alpha,A,descr,x,y) &
                 BIND(C, name='MKL_SPARSE_D_TRSV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            REAL(C_DOUBLE)       , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: x
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_D_TRSV
        END FUNCTION
        FUNCTION MKL_SPARSE_C_TRSV(operation,alpha,A,descr,x,y) &
                 BIND(C, name='MKL_SPARSE_C_TRSV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)          , INTENT(IN) :: operation
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T)   , INTENT(IN) :: A
            TYPE(MATRIX_DESCR)      , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: x
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_C_TRSV
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_TRSV(operation,alpha,A,descr,x,y) &
                 BIND(C, name='MKL_SPARSE_Z_TRSV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)           , INTENT(IN) :: operation
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T)    , INTENT(IN) :: A
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: x
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_Z_TRSV
        END FUNCTION

!   Level 3

!   Computes y = alpha * A * x + beta * y
        FUNCTION MKL_SPARSE_S_MM(operation,alpha,A,descr,layout,x,columns,ldx,beta,y,ldy) &
                 BIND(C, name='MKL_SPARSE_S_MM')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)       , INTENT(IN) :: layout          ! storage scheme for the dense matrix: C-style or Fortran-style
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: x
            INTEGER              , INTENT(IN) :: columns
            INTEGER              , INTENT(IN) :: ldx
            REAL(C_FLOAT)        , INTENT(IN) :: beta
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER              , INTENT(IN) :: ldy
            INTEGER(C_INT) MKL_SPARSE_S_MM
        END FUNCTION
        FUNCTION MKL_SPARSE_D_MM(operation,alpha,A,descr,layout,x,columns,ldx,beta,y,ldy) &
                 BIND(C, name='MKL_SPARSE_D_MM')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            REAL(C_DOUBLE)       , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)       , INTENT(IN) :: layout          ! storage scheme for the dense matrix: C-style or Fortran-style
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: x
            INTEGER              , INTENT(IN) :: columns
            INTEGER              , INTENT(IN) :: ldx
            REAL(C_DOUBLE)       , INTENT(IN) :: beta
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER              , INTENT(IN) :: ldy
            INTEGER(C_INT) MKL_SPARSE_D_MM
        END FUNCTION
        FUNCTION MKL_SPARSE_C_MM(operation,alpha,A,descr,layout,x,columns,ldx,beta,y,ldy) &
                 BIND(C, name='MKL_SPARSE_C_MM')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            COMPLEX(C_FLOAT_COMPLEX)               , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)       , INTENT(IN) :: layout          ! storage scheme for the dense matrix: C-style or Fortran-style
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: x
            INTEGER              , INTENT(IN) :: columns
            INTEGER              , INTENT(IN) :: ldx
            COMPLEX(C_FLOAT_COMPLEX)               , INTENT(IN) :: beta
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER              , INTENT(IN) :: ldy
            INTEGER(C_INT) MKL_SPARSE_C_MM
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_MM(operation,alpha,A,descr,layout,x,columns,ldx,beta,y,ldy) &
                 BIND(C, name='MKL_SPARSE_Z_MM')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            COMPLEX(C_DOUBLE_COMPLEX)               , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)       , INTENT(IN) :: layout          ! storage scheme for the dense matrix: C-style or Fortran-style
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: x
            INTEGER              , INTENT(IN) :: columns
            INTEGER              , INTENT(IN) :: ldx
            COMPLEX(C_DOUBLE_COMPLEX)               , INTENT(IN) :: beta
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER              , INTENT(IN) :: ldy
            INTEGER(C_INT) MKL_SPARSE_Z_MM
        END FUNCTION

!   Solves triangular system y = alpha * A^{-1} * x
        FUNCTION MKL_SPARSE_S_TRSM(operation,alpha,A,descr,layout,x,columns,ldx,y,ldy) &
                 BIND(C, name='MKL_SPARSE_S_TRSM')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)       , INTENT(IN) :: layout          ! storage scheme for the dense matrix: C-style or Fortran-style
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: x
            INTEGER              , INTENT(IN) :: columns
            INTEGER              , INTENT(IN) :: ldx
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER              , INTENT(IN) :: ldy
            INTEGER(C_INT) MKL_SPARSE_S_TRSM
        END FUNCTION
        FUNCTION MKL_SPARSE_D_TRSM(operation,alpha,A,descr,layout,x,columns,ldx,y,ldy) &
                 BIND(C, name='MKL_SPARSE_D_TRSM')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            REAL(C_DOUBLE)       , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)       , INTENT(IN) :: layout          ! storage scheme for the dense matrix: C-style or Fortran-style
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: x
            INTEGER              , INTENT(IN) :: columns
            INTEGER              , INTENT(IN) :: ldx
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER              , INTENT(IN) :: ldy
            INTEGER(C_INT) MKL_SPARSE_D_TRSM
        END FUNCTION
        FUNCTION MKL_SPARSE_C_TRSM(operation,alpha,A,descr,layout,x,columns,ldx,y,ldy) &
                 BIND(C, name='MKL_SPARSE_C_TRSM')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)          , INTENT(IN) :: operation
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T)   , INTENT(IN) :: A
            TYPE(MATRIX_DESCR)      , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)          , INTENT(IN) :: layout          ! storage scheme for the dense matrix: C-style or Fortran-style
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: x
            INTEGER                 , INTENT(IN) :: columns
            INTEGER                 , INTENT(IN) :: ldx
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER                 , INTENT(IN) :: ldy
            INTEGER(C_INT) MKL_SPARSE_C_TRSM
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_TRSM(operation,alpha,A,descr,layout,x,columns,ldx,y,ldy) &
                 BIND(C, name='MKL_SPARSE_Z_TRSM')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)           , INTENT(IN) :: operation
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T)    , INTENT(IN) :: A
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER(C_INT)           , INTENT(IN) :: layout          ! storage scheme for the dense matrix: C-style or Fortran-style
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: x
            INTEGER                  , INTENT(IN) :: columns
            INTEGER                  , INTENT(IN) :: ldx
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER                  , INTENT(IN) :: ldy
            INTEGER(C_INT) MKL_SPARSE_Z_TRSM
        END FUNCTION

!   Sparse-sparse functionality

!   Computes addition of sparse matrices: C = alpha * op(A) + B, result is sparse
        FUNCTION MKL_SPARSE_S_ADD(operation,A,alpha,B,C) &
                 BIND(C, name='MKL_SPARSE_S_ADD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: C
            INTEGER(C_INT) MKL_SPARSE_S_ADD
        END FUNCTION
        FUNCTION MKL_SPARSE_D_ADD(operation,A,alpha,B,C) &
                 BIND(C, name='MKL_SPARSE_D_ADD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            REAL(C_DOUBLE)       , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: C
            INTEGER(C_INT) MKL_SPARSE_D_ADD
        END FUNCTION
        FUNCTION MKL_SPARSE_C_ADD(operation,A,alpha,B,C) &
                 BIND(C, name='MKL_SPARSE_C_ADD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)          , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T)   , INTENT(IN) :: A
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T)   , INTENT(IN) :: B
            TYPE(SPARSE_MATRIX_T)   , INTENT(INOUT) :: C
            INTEGER(C_INT) MKL_SPARSE_C_ADD
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_ADD(operation,A,alpha,B,C) &
                 BIND(C, name='MKL_SPARSE_Z_ADD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            COMPLEX(C_DOUBLE_COMPLEX)               , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: C
            INTEGER(C_INT) MKL_SPARSE_Z_ADD
        END FUNCTION

!   Computes multiplication of sparse matrices: C = op(A) * B, result is sparse
        FUNCTION MKL_SPARSE_SPMM(operation,A,B,C) &
                 BIND(C, name='MKL_SPARSE_SPMM')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: C
            INTEGER(C_INT) MKL_SPARSE_SPMM
        END FUNCTION

        FUNCTION MKL_SPARSE_SP2M(opA,descrA,A,opB,descrB,B,req,C) &
                 BIND(C, name='MKL_SPARSE_SP2M')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: opA
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrA
            INTEGER(C_INT)       , INTENT(IN) :: opB
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrB
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            INTEGER(C_INT), INTENT(IN) :: req
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: C
            INTEGER(C_INT) MKL_SPARSE_SP2M
        END FUNCTION

        FUNCTION MKL_SPARSE_SYRK(operation,A,C) &
                 BIND(C, name='MKL_SPARSE_SYRK')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: C
            INTEGER(C_INT) MKL_SPARSE_SYRK
        END FUNCTION

        FUNCTION MKL_SPARSE_SYPR(operation,A,B,descrB,C,req) &
                 BIND(C, name='MKL_SPARSE_SYPR')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descrB
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: C
            INTEGER(C_INT)       , INTENT(IN) :: req
            INTEGER(C_INT) MKL_SPARSE_SYPR
        END FUNCTION

        FUNCTION MKL_SPARSE_S_SYPRD(op,A,B,layoutB,ldb,alpha,beta,C,layoutC,ldc) & !
                 BIND(C, name='MKL_SPARSE_S_SYPRD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)              , INTENT(IN)   :: op
            TYPE(SPARSE_MATRIX_T)       , INTENT(IN)   :: A
            REAL(C_FLOAT)               , INTENT(IN) , DIMENSION(*) :: B
            INTEGER(C_INT)              , INTENT(IN)   :: layoutB
            INTEGER(C_INT)              , INTENT(IN)   :: ldb
            REAL(C_FLOAT)               , INTENT(IN)   :: alpha
            REAL(C_FLOAT)               , INTENT(IN)   :: beta
            REAL(C_FLOAT)               , INTENT(OUT), DIMENSION(*) :: C
            INTEGER(C_INT)              , INTENT(IN)   :: layoutC
            INTEGER(C_INT)              , INTENT(IN)   :: ldc
            INTEGER(C_INT) MKL_SPARSE_S_SYPRD
        END FUNCTION
        FUNCTION MKL_SPARSE_D_SYPRD(op,A,B,layoutB,ldb,alpha,beta,C,layoutC,ldc) & !
                 BIND(C, name='MKL_SPARSE_D_SYPRD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)              , INTENT(IN)   :: op
            TYPE(SPARSE_MATRIX_T)       , INTENT(IN)   :: A
            REAL(C_DOUBLE)              , INTENT(IN) , DIMENSION(*) :: B
            INTEGER(C_INT)              , INTENT(IN)   :: layoutB
            INTEGER(C_INT)              , INTENT(IN)   :: ldb
            REAL(C_DOUBLE)              , INTENT(IN)   :: alpha
            REAL(C_DOUBLE)              , INTENT(IN)   :: beta
            REAL(C_DOUBLE)              , INTENT(OUT), DIMENSION(*) :: C
            INTEGER(C_INT)              , INTENT(IN)   :: layoutC
            INTEGER(C_INT)              , INTENT(IN)   :: ldc
            INTEGER(C_INT) MKL_SPARSE_D_SYPRD
        END FUNCTION
        FUNCTION MKL_SPARSE_C_SYPRD(op,A,B,layoutB,ldb,alpha,beta,C,layoutC,ldc) & !
                 BIND(C, name='MKL_SPARSE_C_SYPRD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)              , INTENT(IN)   :: op
            TYPE(SPARSE_MATRIX_T)       , INTENT(IN)   :: A
            REAL(C_FLOAT_COMPLEX)       , INTENT(IN) , DIMENSION(*) :: B
            INTEGER(C_INT)              , INTENT(IN)   :: layoutB
            INTEGER(C_INT)              , INTENT(IN)   :: ldb
            REAL(C_FLOAT_COMPLEX)       , INTENT(IN)   :: alpha
            REAL(C_FLOAT_COMPLEX)       , INTENT(IN)   :: beta
            REAL(C_FLOAT_COMPLEX)       , INTENT(OUT), DIMENSION(*) :: C
            INTEGER(C_INT)              , INTENT(IN)   :: layoutC
            INTEGER(C_INT)              , INTENT(IN)   :: ldc
            INTEGER(C_INT) MKL_SPARSE_C_SYPRD
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_SYPRD(op,A,B,layoutB,ldb,alpha,beta,C,layoutC,ldc) & !
                 BIND(C, name='MKL_SPARSE_Z_SYPRD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)              , INTENT(IN)   :: op
            TYPE(SPARSE_MATRIX_T)       , INTENT(IN)   :: A
            REAL(C_DOUBLE_COMPLEX)      , INTENT(IN) , DIMENSION(*) :: B
            INTEGER(C_INT)              , INTENT(IN)   :: layoutB
            INTEGER(C_INT)              , INTENT(IN)   :: ldb
            REAL(C_DOUBLE_COMPLEX)      , INTENT(IN)   :: alpha
            REAL(C_DOUBLE_COMPLEX)      , INTENT(IN)   :: beta
            REAL(C_DOUBLE_COMPLEX)      , INTENT(OUT), DIMENSION(*) :: C
            INTEGER(C_INT)              , INTENT(IN)   :: layoutC
            INTEGER(C_INT)              , INTENT(IN)   :: ldc
            INTEGER(C_INT) MKL_SPARSE_Z_SYPRD
        END FUNCTION

        FUNCTION MKL_SPARSE_ORDER(A) &
                 BIND(C, name='MKL_SPARSE_ORDER')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            INTEGER(C_INT) MKL_SPARSE_ORDER
        END FUNCTION

!   Computes multiplication of sparse matrices: C = op(A) * B, result is dense
        FUNCTION MKL_SPARSE_S_SPMMD(operation,A,B,layout,C,ldc) &
                 BIND(C, name='MKL_SPARSE_S_SPMMD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            INTEGER(C_INT)       , INTENT(IN) :: layout
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER              , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_S_SPMMD
        END FUNCTION
        FUNCTION MKL_SPARSE_D_SPMMD(operation,A,B,layout,C,ldc) &
                 BIND(C, name='MKL_SPARSE_D_SPMMD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            INTEGER(C_INT)      , INTENT(IN) :: layout
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER              , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_D_SPMMD
        END FUNCTION
        FUNCTION MKL_SPARSE_C_SPMMD(operation,A,B,layout,C,ldc) &
                 BIND(C, name='MKL_SPARSE_C_SPMMD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            INTEGER(C_INT)      , INTENT(IN) :: layout
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER             , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_C_SPMMD
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_SPMMD(operation,A,B,layout,C,ldc) &
                 BIND(C, name='MKL_SPARSE_Z_SPMMD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            INTEGER(C_INT)      , INTENT(IN) :: layout
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER              , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_Z_SPMMD
        END FUNCTION

        FUNCTION MKL_SPARSE_S_SP2MD(opA,descrA,A,opB,descrB,B,alpha,beta,layout,C,ldc) & !
                 BIND(C, name='MKL_SPARSE_S_SP2MD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: opA
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrA
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            INTEGER(C_INT)       , INTENT(IN) :: opB
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrB
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            REAL(C_FLOAT)        , INTENT(IN) :: beta
            INTEGER(C_INT)       , INTENT(IN) :: layout
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER              , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_S_SP2MD
        END FUNCTION
        FUNCTION MKL_SPARSE_D_SP2MD(opA,descrA,A,opB,descrB,B,alpha,beta,layout,C,ldc) & !
                 BIND(C, name='MKL_SPARSE_D_SP2MD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: opA
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrA
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            INTEGER(C_INT)       , INTENT(IN) :: opB
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrB
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            REAL(C_DOUBLE)        , INTENT(IN) :: alpha
            REAL(C_DOUBLE)        , INTENT(IN) :: beta
            INTEGER(C_INT)      , INTENT(IN) :: layout
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER              , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_D_SP2MD
        END FUNCTION
        FUNCTION MKL_SPARSE_C_SP2MD(opA,descrA,A,opB,descrB,B,alpha,beta,layout,C,ldc) & !
                 BIND(C, name='MKL_SPARSE_C_SP2MD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: opA
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrA
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            INTEGER(C_INT)       , INTENT(IN) :: opB
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrB
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            COMPLEX(C_FLOAT_COMPLEX)        , INTENT(IN) :: alpha
            COMPLEX(C_FLOAT_COMPLEX)        , INTENT(IN) :: beta
            INTEGER(C_INT)      , INTENT(IN) :: layout
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER             , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_C_SP2MD
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_SP2MD(opA,descrA,A,opB,descrB,B,alpha,beta,layout,C,ldc) & !
                 BIND(C, name='MKL_SPARSE_Z_SP2MD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: opA
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrA
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            INTEGER(C_INT)       , INTENT(IN) :: opB
            TYPE(MATRIX_DESCR)       , INTENT(IN) :: descrB
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: B
            COMPLEX(C_DOUBLE_COMPLEX)        , INTENT(IN) :: alpha
            COMPLEX(C_DOUBLE_COMPLEX)        , INTENT(IN) :: beta
            INTEGER(C_INT)      , INTENT(IN) :: layout
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER              , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_Z_SP2MD
        END FUNCTION

        FUNCTION MKL_SPARSE_S_SYRKD(operation,A,alpha,beta,C,layout,ldc) & !
                 BIND(C, name='MKL_SPARSE_S_SYRKD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            REAL(C_FLOAT)        , INTENT(IN) :: beta
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER(C_INT)       , INTENT(IN) :: layout
            INTEGER              , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_S_SYRKD
        END FUNCTION
        FUNCTION MKL_SPARSE_D_SYRKD(operation,A,alpha,beta,C,layout,ldc) &
                 BIND(C, name='MKL_SPARSE_D_SYRKD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            REAL(C_DOUBLE)        , INTENT(IN) :: alpha
            REAL(C_DOUBLE)        , INTENT(IN) :: beta
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER(C_INT)      , INTENT(IN) :: layout
            INTEGER              , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_D_SYRKD
        END FUNCTION
        FUNCTION MKL_SPARSE_C_SYRKD(operation,A,alpha,beta,C,layout,ldc) &
                 BIND(C, name='MKL_SPARSE_C_SYRKD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            COMPLEX(C_FLOAT_COMPLEX)        , INTENT(IN) :: alpha
            COMPLEX(C_FLOAT_COMPLEX)        , INTENT(IN) :: beta
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER(C_INT)      , INTENT(IN) :: layout
            INTEGER             , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_C_SYRKD
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_SYRKD(operation,A,alpha,beta,C,layout,ldc) &
                 BIND(C, name='MKL_SPARSE_Z_SYRKD')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            COMPLEX(C_DOUBLE_COMPLEX)        , INTENT(IN) :: alpha
            COMPLEX(C_DOUBLE_COMPLEX)        , INTENT(IN) :: beta
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: C
            INTEGER(C_INT)      , INTENT(IN) :: layout
            INTEGER              , INTENT(IN) :: ldc
            INTEGER(C_INT) MKL_SPARSE_Z_SYRKD
        END FUNCTION

        FUNCTION MKL_SPARSE_S_SYMGS_MV(operation,A,descr,alpha,b,x,y) &
                 BIND(C, name='MKL_SPARSE_S_SYMGS_MV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: b
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: x
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_S_SYMGS_MV
        END FUNCTION
        FUNCTION MKL_SPARSE_D_SYMGS_MV(operation,A,descr,alpha,b,x,y) &
                 BIND(C, name='MKL_SPARSE_D_SYMGS_MV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_DOUBLE)        , INTENT(IN) :: alpha
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: b
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: x
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_D_SYMGS_MV
        END FUNCTION
        FUNCTION MKL_SPARSE_C_SYMGS_MV(operation,A,descr,alpha,b,x,y) &
                 BIND(C, name='MKL_SPARSE_C_SYMGS_MV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_FLOAT_COMPLEX)        , INTENT(IN) :: alpha
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: b
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: x
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_C_SYMGS_MV
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_SYMGS_MV(operation,A,descr,alpha,b,x,y) &
                 BIND(C, name='MKL_SPARSE_Z_SYMGS_MV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_DOUBLE_COMPLEX)        , INTENT(IN) :: alpha
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: b
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: x
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            INTEGER(C_INT) MKL_SPARSE_Z_SYMGS_MV
        END FUNCTION

        FUNCTION MKL_SPARSE_S_SYMGS(operation,A,descr,alpha,b,x) &
                 BIND(C, name='MKL_SPARSE_S_SYMGS')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            REAL(C_FLOAT), INTENT(IN), DIMENSION(*) :: b
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: x
            INTEGER(C_INT) MKL_SPARSE_S_SYMGS
        END FUNCTION
        FUNCTION MKL_SPARSE_D_SYMGS(operation,A,descr,alpha,b,x) &
                 BIND(C, name='MKL_SPARSE_D_SYMGS')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_DOUBLE)        , INTENT(IN) :: alpha
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(*) :: b
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: x
            INTEGER(C_INT) MKL_SPARSE_D_SYMGS
        END FUNCTION
        FUNCTION MKL_SPARSE_C_SYMGS(operation,A,descr,alpha,b,x) &
                 BIND(C, name='MKL_SPARSE_C_SYMGS')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_FLOAT_COMPLEX)        , INTENT(IN) :: alpha
            COMPLEX(C_FLOAT_COMPLEX), INTENT(IN), DIMENSION(*) :: b
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: x
            INTEGER(C_INT) MKL_SPARSE_C_SYMGS
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_SYMGS(operation,A,descr,alpha,b,x) &
                 BIND(C, name='MKL_SPARSE_Z_SYMGS')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)      , INTENT(IN) :: operation
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_DOUBLE_COMPLEX)        , INTENT(IN) :: alpha
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: b
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: x
            INTEGER(C_INT) MKL_SPARSE_Z_SYMGS
        END FUNCTION

        FUNCTION MKL_SPARSE_S_DOTMV(operation,alpha,A,descr,x,beta,y,d) &
                 BIND(C, name='MKL_SPARSE_S_DOTMV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)       , INTENT(IN) :: operation
            REAL(C_FLOAT)        , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: x
            REAL(C_FLOAT)        , INTENT(IN) :: beta
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: y
            REAL(C_FLOAT), INTENT(INOUT), DIMENSION(*) :: d
            INTEGER(C_INT) MKL_SPARSE_S_DOTMV
        END FUNCTION
        FUNCTION MKL_SPARSE_D_DOTMV(operation,alpha,A,descr,x,beta,y,d) &
                 BIND(C, name='MKL_SPARSE_D_DOTMV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)      , INTENT(IN) :: operation
            REAL(C_DOUBLE)        , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: x
            REAL(C_DOUBLE)        , INTENT(IN) :: beta
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: y
            REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(*) :: d
            INTEGER(C_INT) MKL_SPARSE_D_DOTMV
        END FUNCTION
        FUNCTION MKL_SPARSE_C_DOTMV(operation,alpha,A,descr,x,beta,y,d) &
                 BIND(C, name='MKL_SPARSE_C_DOTMV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_FLOAT_COMPLEX, C_FLOAT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)      , INTENT(IN) :: operation
            COMPLEX(C_FLOAT_COMPLEX)        , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: x
            COMPLEX(C_FLOAT_COMPLEX)        , INTENT(IN) :: beta
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: y
            COMPLEX(C_FLOAT_COMPLEX), INTENT(INOUT), DIMENSION(*) :: d
            INTEGER(C_INT) MKL_SPARSE_C_DOTMV
        END FUNCTION
        FUNCTION MKL_SPARSE_Z_DOTMV(operation,alpha,A,descr,x,beta,y,d) &
                 BIND(C, name='MKL_SPARSE_Z_DOTMV')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE_COMPLEX, C_DOUBLE
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            INTEGER(C_INT)      , INTENT(IN) :: operation
            COMPLEX(C_DOUBLE_COMPLEX)        , INTENT(IN) :: alpha
            TYPE(SPARSE_MATRIX_T), INTENT(IN) :: A
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr          ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: x
            COMPLEX(C_DOUBLE_COMPLEX)        , INTENT(IN) :: beta
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN), DIMENSION(*) :: y
            COMPLEX(C_DOUBLE_COMPLEX), INTENT(INOUT), DIMENSION(*) :: d
            INTEGER(C_INT) MKL_SPARSE_Z_DOTMV
        END FUNCTION

        FUNCTION MKL_SPARSE_SET_DOTMV_HINT(A,operation,descr,expected_calls) &
                 BIND(C, name='MKL_SPARSE_SET_DOTMV_HINT')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT)      , INTENT(IN) :: operation     ! SPARSE_OPERATION_NON_TRANSPOSE is default value for
                                                               ! infinite amount of calls
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr         ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER              , INTENT(IN) :: expected_calls
            INTEGER(C_INT) MKL_SPARSE_SET_DOTMV_HINT
        END FUNCTION

        FUNCTION MKL_SPARSE_SET_SYMGS_HINT(A,operation,descr,expected_calls) &
                 BIND(C, name='MKL_SPARSE_SET_SYMGS_HINT')
            USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT
            IMPORT SPARSE_MATRIX_T
            IMPORT MATRIX_DESCR
            TYPE(SPARSE_MATRIX_T), INTENT(INOUT) :: A
            INTEGER(C_INT)      , INTENT(IN) :: operation     ! SPARSE_OPERATION_NON_TRANSPOSE is default value for
                                                               ! infinite amount of calls
            TYPE(MATRIX_DESCR)   , INTENT(IN) :: descr         ! sparse_matrix_type_t + sparse_fill_mode_t + sparse_diag_type_t
            INTEGER              , INTENT(IN) :: expected_calls
            INTEGER(C_INT) MKL_SPARSE_SET_SYMGS_HINT
        END FUNCTION

    END INTERFACE

    END MODULE MKL_SPBLAS
