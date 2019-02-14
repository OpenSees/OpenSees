! file: mkl_df.f90
!===============================================================================
! Copyright 2006-2018 Intel Corporation All Rights Reserved.
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

!++
!  Fortran 90 DF interface.
!--

      MODULE MKL_DF_TYPE

!++
!  Definitions for DF functions return values (errors, warnings)
!--

!  "No error" status
      INTEGER(KIND=4) DF_STATUS_OK
      PARAMETER (DF_STATUS_OK = 0)

!  Common errors (-1..-999)
      INTEGER(KIND=4) DF_ERROR_CPU_NOT_SUPPORTED
      PARAMETER (DF_ERROR_CPU_NOT_SUPPORTED             =    -1)

!  Data fitting errors (-1000..-1999)
      INTEGER(KIND=4) DF_ERROR_NULL_TASK_DESCRIPTOR
      INTEGER(KIND=4) DF_ERROR_MEM_FAILURE
      INTEGER(KIND=4) DF_ERROR_METHOD_NOT_SUPPORTED
      INTEGER(KIND=4) DF_ERROR_COMP_TYPE_NOT_SUPPORTED
      INTEGER(KIND=4) DF_ERROR_NULL_PTR
      INTEGER(KIND=4) DF_ERROR_BAD_NX
      INTEGER(KIND=4) DF_ERROR_BAD_X
      INTEGER(KIND=4) DF_ERROR_BAD_X_HINT
      INTEGER(KIND=4) DF_ERROR_BAD_NY
      INTEGER(KIND=4) DF_ERROR_BAD_Y
      INTEGER(KIND=4) DF_ERROR_BAD_Y_HINT
      INTEGER(KIND=4) DF_ERROR_BAD_SPLINE_ORDER
      INTEGER(KIND=4) DF_ERROR_BAD_SPLINE_TYPE
      INTEGER(KIND=4) DF_ERROR_BAD_IC_TYPE
      INTEGER(KIND=4) DF_ERROR_BAD_IC
      INTEGER(KIND=4) DF_ERROR_BAD_BC_TYPE
      INTEGER(KIND=4) DF_ERROR_BAD_BC
      INTEGER(KIND=4) DF_ERROR_BAD_PP_COEFF
      INTEGER(KIND=4) DF_ERROR_BAD_PP_COEFF_HINT
      INTEGER(KIND=4) DF_ERROR_BAD_PERIODIC_VAL
      INTEGER(KIND=4) DF_ERROR_BAD_DATA_ATTR
      INTEGER(KIND=4) DF_ERROR_BAD_DATA_IDX
      INTEGER(KIND=4) DF_ERROR_BAD_NSITE
      INTEGER(KIND=4) DF_ERROR_BAD_SITE
      INTEGER(KIND=4) DF_ERROR_BAD_SITE_HINT
      INTEGER(KIND=4) DF_ERROR_BAD_NDORDER
      INTEGER(KIND=4) DF_ERROR_BAD_DORDER
      INTEGER(KIND=4) DF_ERROR_BAD_DATA_HINT
      INTEGER(KIND=4) DF_ERROR_BAD_INTERP
      INTEGER(KIND=4) DF_ERROR_BAD_INTERP_HINT
      INTEGER(KIND=4) DF_ERROR_BAD_CELL_IDX
      INTEGER(KIND=4) DF_ERROR_BAD_NLIM
      INTEGER(KIND=4) DF_ERROR_BAD_LLIM
      INTEGER(KIND=4) DF_ERROR_BAD_RLIM
      INTEGER(KIND=4) DF_ERROR_BAD_INTEGR
      INTEGER(KIND=4) DF_ERROR_BAD_INTEGR_HINT
      INTEGER(KIND=4) DF_ERROR_BAD_LOOKUP_INTERP_SITE
      INTEGER(KIND=4) DF_ERROR_BAD_CHECK_FLAG
      PARAMETER (DF_ERROR_NULL_TASK_DESCRIPTOR          = -1000)
      PARAMETER (DF_ERROR_MEM_FAILURE                   = -1001)
      PARAMETER (DF_ERROR_METHOD_NOT_SUPPORTED          = -1002)
      PARAMETER (DF_ERROR_COMP_TYPE_NOT_SUPPORTED       = -1003)
      PARAMETER (DF_ERROR_NULL_PTR                      = -1037)
      PARAMETER (DF_ERROR_BAD_NX                        = -1004)
      PARAMETER (DF_ERROR_BAD_X                         = -1005)
      PARAMETER (DF_ERROR_BAD_X_HINT                    = -1006)
      PARAMETER (DF_ERROR_BAD_NY                        = -1007)
      PARAMETER (DF_ERROR_BAD_Y                         = -1008)
      PARAMETER (DF_ERROR_BAD_Y_HINT                    = -1009)
      PARAMETER (DF_ERROR_BAD_SPLINE_ORDER              = -1010)
      PARAMETER (DF_ERROR_BAD_SPLINE_TYPE               = -1011)
      PARAMETER (DF_ERROR_BAD_IC_TYPE                   = -1012)
      PARAMETER (DF_ERROR_BAD_IC                        = -1013)
      PARAMETER (DF_ERROR_BAD_BC_TYPE                   = -1014)
      PARAMETER (DF_ERROR_BAD_BC                        = -1015)
      PARAMETER (DF_ERROR_BAD_PP_COEFF                  = -1016)
      PARAMETER (DF_ERROR_BAD_PP_COEFF_HINT             = -1017)
      PARAMETER (DF_ERROR_BAD_PERIODIC_VAL              = -1018)
      PARAMETER (DF_ERROR_BAD_DATA_ATTR                 = -1019)
      PARAMETER (DF_ERROR_BAD_DATA_IDX                  = -1020)
      PARAMETER (DF_ERROR_BAD_NSITE                     = -1021)
      PARAMETER (DF_ERROR_BAD_SITE                      = -1022)
      PARAMETER (DF_ERROR_BAD_SITE_HINT                 = -1023)
      PARAMETER (DF_ERROR_BAD_NDORDER                   = -1024)
      PARAMETER (DF_ERROR_BAD_DORDER                    = -1025)
      PARAMETER (DF_ERROR_BAD_DATA_HINT                 = -1026)
      PARAMETER (DF_ERROR_BAD_INTERP                    = -1027)
      PARAMETER (DF_ERROR_BAD_INTERP_HINT               = -1028)
      PARAMETER (DF_ERROR_BAD_CELL_IDX                  = -1029)
      PARAMETER (DF_ERROR_BAD_NLIM                      = -1030)
      PARAMETER (DF_ERROR_BAD_LLIM                      = -1031)
      PARAMETER (DF_ERROR_BAD_RLIM                      = -1032)
      PARAMETER (DF_ERROR_BAD_INTEGR                    = -1033)
      PARAMETER (DF_ERROR_BAD_INTEGR_HINT               = -1034)
      PARAMETER (DF_ERROR_BAD_LOOKUP_INTERP_SITE        = -1035)
      PARAMETER (DF_ERROR_BAD_CHECK_FLAG                = -1036)

!  Internal errors caused by internal routines of the functions
      INTEGER(KIND=4) VSL_DF_ERROR_INTERNAL_C1
      INTEGER(KIND=4) VSL_DF_ERROR_INTERNAL_C2
      PARAMETER (VSL_DF_ERROR_INTERNAL_C1               = -1500)
      PARAMETER (VSL_DF_ERROR_INTERNAL_C2               = -1501)

!  User-defined callback status
      INTEGER(KIND=4) DF_STATUS_EXACT_RESULT
      PARAMETER (DF_STATUS_EXACT_RESULT                 = 1000)

!++
!  MACROS USED IN DATA FITTING EDIT AND COMPUTE ROUTINES
!--

!  DF EditTask routine is way to edit input and output parameters of the task
!  Macros below define parameters available for modification
      INTEGER DF_X
      INTEGER DF_Y
      INTEGER DF_IC
      INTEGER DF_BC
      INTEGER DF_PP_SCOEFF
      INTEGER DF_NX
      INTEGER DF_XHINT
      INTEGER DF_NY
      INTEGER DF_YHINT
      INTEGER DF_SPLINE_ORDER
      INTEGER DF_SPLINE_TYPE
      INTEGER DF_IC_TYPE
      INTEGER DF_BC_TYPE
      INTEGER DF_PP_COEFF_HINT
      INTEGER DF_CHECK_FLAG
      PARAMETER (DF_X                                      =  1)
      PARAMETER (DF_Y                                      =  2)
      PARAMETER (DF_IC                                     =  3)
      PARAMETER (DF_BC                                     =  4)
      PARAMETER (DF_PP_SCOEFF                              =  5)
      PARAMETER (DF_NX                                     = 14)
      PARAMETER (DF_XHINT                                  = 15)
      PARAMETER (DF_NY                                     = 16)
      PARAMETER (DF_YHINT                                  = 17)
      PARAMETER (DF_SPLINE_ORDER                           = 18)
      PARAMETER (DF_SPLINE_TYPE                            = 19)
      PARAMETER (DF_IC_TYPE                                = 20)
      PARAMETER (DF_BC_TYPE                                = 21)
      PARAMETER (DF_PP_COEFF_HINT                          = 22)
      PARAMETER (DF_CHECK_FLAG                             = 23)

!++
!  SPLINE ORDERS SUPPORTED IN DATA FITTING ROUTINES
!--
      INTEGER(KIND=4) DF_PP_STD
      INTEGER(KIND=4) DF_PP_LINEAR
      INTEGER(KIND=4) DF_PP_QUADRATIC
      INTEGER(KIND=4) DF_PP_CUBIC
      PARAMETER (DF_PP_STD                                 =  0)
      PARAMETER (DF_PP_LINEAR                              =  2)
      PARAMETER (DF_PP_QUADRATIC                           =  3)
      PARAMETER (DF_PP_CUBIC                               =  4)

!++
!  SPLINE TYPES SUPPORTED IN DATA FITTING ROUTINES
!--

      INTEGER(KIND=4) DF_PP_DEFAULT
      INTEGER(KIND=4) DF_PP_SUBBOTIN
      INTEGER(KIND=4) DF_PP_NATURAL
      INTEGER(KIND=4) DF_PP_HERMITE
      INTEGER(KIND=4) DF_PP_BESSEL
      INTEGER(KIND=4) DF_PP_AKIMA
      INTEGER(KIND=4) DF_LOOKUP_INTERPOLANT
      INTEGER(KIND=4) DF_CR_STEPWISE_CONST_INTERPOLANT
      INTEGER(KIND=4) DF_CL_STEPWISE_CONST_INTERPOLANT
      INTEGER(KIND=4) DF_PP_HYMAN
      PARAMETER (DF_PP_DEFAULT                             =  0)
      PARAMETER (DF_PP_SUBBOTIN                            =  1)
      PARAMETER (DF_PP_NATURAL                             =  2)
      PARAMETER (DF_PP_HERMITE                             =  3)
      PARAMETER (DF_PP_BESSEL                              =  4)
      PARAMETER (DF_PP_AKIMA                               =  5)
      PARAMETER (DF_LOOKUP_INTERPOLANT                     =  6)
      PARAMETER (DF_CR_STEPWISE_CONST_INTERPOLANT          =  7)
      PARAMETER (DF_CL_STEPWISE_CONST_INTERPOLANT          =  8)
      PARAMETER (DF_PP_HYMAN                               =  9)

!++
!  TYPES OF BOUNDARY CONDITIONS USED IN SPLINE CONSTRUCTION
!--
      INTEGER DF_NO_BC
      INTEGER DF_BC_NOT_A_KNOT
      INTEGER DF_BC_FREE_END
      INTEGER DF_BC_1ST_LEFT_DER
      INTEGER DF_BC_1ST_RIGHT_DER
      INTEGER DF_BC_2ND_LEFT_DER
      INTEGER DF_BC_2ND_RIGHT_DER
      INTEGER DF_BC_PERIODIC
      INTEGER DF_BC_Q_VAL
      PARAMETER (DF_NO_BC                                 =   0)
      PARAMETER (DF_BC_NOT_A_KNOT                         =   1)
      PARAMETER (DF_BC_FREE_END                           =   2)
      PARAMETER (DF_BC_1ST_LEFT_DER                       =   4)
      PARAMETER (DF_BC_1ST_RIGHT_DER                      =   8)
      PARAMETER (DF_BC_2ND_LEFT_DER                       =  16)
      PARAMETER (DF_BC_2ND_RIGHT_DER                      =  32)
      PARAMETER (DF_BC_PERIODIC                           =  64)
      PARAMETER (DF_BC_Q_VAL                              = 128)

!++
!  TYPES OF INTERNAL CONDITIONS USED IN SPLINE CONSTRUCTION
!--
      INTEGER DF_NO_IC
      INTEGER DF_IC_1ST_DER
      INTEGER DF_IC_2ND_DER
      INTEGER DF_IC_Q_KNOT
      PARAMETER (DF_NO_IC                                  =  0)
      PARAMETER (DF_IC_1ST_DER                             =  1)
      PARAMETER (DF_IC_2ND_DER                             =  2)
      PARAMETER (DF_IC_Q_KNOT                              =  8)

!++
!  TYPES OF SUPPORTED HINTS
!--
      INTEGER(KIND=4) DF_NO_HINT
      INTEGER(KIND=4) DF_NON_UNIFORM_PARTITION
      INTEGER(KIND=4) DF_QUASI_UNIFORM_PARTITION
      INTEGER(KIND=4) DF_UNIFORM_PARTITION
      INTEGER(KIND=4) DF_MATRIX_STORAGE_ROWS
      INTEGER(KIND=4) DF_MATRIX_STORAGE_COLS
      INTEGER(KIND=4) DF_SORTED_DATA
      INTEGER(KIND=4) DF_1ST_COORDINATE
      INTEGER(KIND=4) DF_MATRIX_STORAGE_FUNCS_SITES_DERS
      INTEGER(KIND=4) DF_MATRIX_STORAGE_FUNCS_DERS_SITES
      INTEGER(KIND=4) DF_MATRIX_STORAGE_SITES_FUNCS_DERS
      INTEGER(KIND=4) DF_MATRIX_STORAGE_SITES_DERS_FUNCS
      PARAMETER (DF_NO_HINT                       = Z"00000000")
      PARAMETER (DF_NON_UNIFORM_PARTITION         = Z"00000001")
      PARAMETER (DF_QUASI_UNIFORM_PARTITION       = Z"00000002")
      PARAMETER (DF_UNIFORM_PARTITION             = Z"00000004")
      PARAMETER (DF_MATRIX_STORAGE_ROWS           = Z"00000010")
      PARAMETER (DF_MATRIX_STORAGE_COLS           = Z"00000020")
      PARAMETER (DF_SORTED_DATA                   = Z"00000040")
      PARAMETER (DF_1ST_COORDINATE                = Z"00000080")
      PARAMETER (DF_MATRIX_STORAGE_FUNCS_SITES_DERS =                   &
     &                                     DF_MATRIX_STORAGE_ROWS)
      PARAMETER (DF_MATRIX_STORAGE_FUNCS_DERS_SITES =                   &
     &                                     DF_MATRIX_STORAGE_COLS)
      PARAMETER (DF_MATRIX_STORAGE_SITES_FUNCS_DERS = Z"00000100")
      PARAMETER (DF_MATRIX_STORAGE_SITES_DERS_FUNCS = Z"00000200")

!++
!  TYPES OF APRIORI INFORMATION ABOUT DATA STRUCTURE
!--

      INTEGER(KIND=4) DF_NO_APRIORI_INFO
      INTEGER(KIND=4) DF_APRIORI_MOST_LIKELY_CELL
      PARAMETER (DF_NO_APRIORI_INFO               = Z"00000000")
      PARAMETER (DF_APRIORI_MOST_LIKELY_CELL      = Z"00000001")

!++
!  ESTIMATES TO BE COMUTED WITH DATA FITTING COMPUTE ROUTINE
!--

      INTEGER(KIND=4) DF_INTERP
      INTEGER(KIND=4) DF_CELL
      INTEGER(KIND=4) DF_INTERP_USER_CELL
      PARAMETER (DF_INTERP                        = Z"00000001")
      PARAMETER (DF_CELL                          = Z"00000002")
      PARAMETER (DF_INTERP_USER_CELL              = Z"00000004")

!++
!  METHODS TO BE USED FOR EVALUATION OF THE SPLINE RELATED ESTIMATES
!--

      INTEGER DF_METHOD_STD
      INTEGER DF_METHOD_PP
      PARAMETER (DF_METHOD_STD                            =   0)
      PARAMETER (DF_METHOD_PP                             =   1)

!++
! POSSIBLE VALUES FOR DF_CHECK_FLAG
!--

      INTEGER DF_ENABLE_CHECK_FLAG
      INTEGER DF_DISABLE_CHECK_FLAG
      PARAMETER (DF_ENABLE_CHECK_FLAG                     =   0)
      PARAMETER (DF_DISABLE_CHECK_FLAG                    =   1)

!++
!  SPLINE FORMATS SUPPORTED IN SPLINE CONSTRUCTION ROUTINE
!--

      INTEGER DF_PP_SPLINE
      PARAMETER (DF_PP_SPLINE                             =   0)

!++
! VALUES OF FLAG INDICATING WHICH, LEFT OR RIGHT, INTEGRATION LIMITS
! ARE PASSED BY INTEGRATION ROUTINE INTO SEARCH CALLBACK
!--

      INTEGER(KIND=4) DF_INTEGR_SEARCH_CB_LLIM_FLAG
      INTEGER(KIND=4) DF_INTEGR_SEARCH_CB_RLIM_FLAG
      PARAMETER (DF_INTEGR_SEARCH_CB_LLIM_FLAG    = Z"00000000")
      PARAMETER (DF_INTEGR_SEARCH_CB_RLIM_FLAG    = Z"00000001")

!++
!  TYPEDEFS
!--

      TYPE DF_TASK
         INTEGER(KIND=4) descriptor1
         INTEGER(KIND=4) descriptor2
      END TYPE DF_TASK

!++
!  DATA FITTING INTERPOLATION CALLBACK INTERNAL PARAMETERS STRUCTURE
!--

      TYPE DF_INTERP_CALLBACK_LIBRARY_PARAMS
         INTEGER(KIND=4) ::reserved1
      END TYPE DF_INTERP_CALLBACK_LIBRARY_PARAMS

!++
!  DATA FITTING INTEGRATION CALLBACK INTERNAL PARAMETERS STRUCTURE
!--

      TYPE DF_INTEGR_CALLBACK_LIBRARY_PARAMS
         INTEGER(KIND=4) ::reserved1
      END TYPE DF_INTEGR_CALLBACK_LIBRARY_PARAMS

!++
!  DATA FITTING SEARCH CALLBACK INTERNAL PARAMETERS STRUCTURE
!--

      TYPE DF_SEARCH_CALLBACK_LIBRARY_PARAMS
         INTEGER(KIND=4) ::limit_type_flag
      END TYPE DF_SEARCH_CALLBACK_LIBRARY_PARAMS

      END MODULE MKL_DF_TYPE

      MODULE MKL_DF

      USE MKL_DF_TYPE

!++
!  DF CONSTRUCTOR FUNCTION DECLARATIONS.
!--

!  NewTask1D - 1d task creation/initialization
      INTERFACE
         INTEGER FUNCTION dfsnewtask1d(task,nx,x,xhint,ny,y,yhint)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: nx
            REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: x
            INTEGER,INTENT(IN),OPTIONAL                   :: xhint
            INTEGER,INTENT(IN),OPTIONAL                   :: ny
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: y
            INTEGER,INTENT(IN),OPTIONAL                   :: yhint
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdnewtask1d(task,nx,x,xhint,ny,y,yhint)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: nx
            REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: x
            INTEGER,INTENT(IN),OPTIONAL                   :: xhint
            INTEGER,INTENT(IN),OPTIONAL                   :: ny
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: y
            INTEGER,INTENT(IN),OPTIONAL                   :: yhint
         END FUNCTION
      END INTERFACE

!++
!  DF EDITOR FUNCTION DECLARATIONS.
!--

!  Modifies a pointer to an array held in a Data Fitting task descriptor
      INTERFACE
         INTEGER FUNCTION dfseditptr(task,ptr_attr,ptr)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: ptr_attr
            REAL(KIND=4),DIMENSION(*),INTENT(IN) :: ptr
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdeditptr(task,ptr_attr,ptr)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: ptr_attr
            REAL(KIND=8),DIMENSION(*),INTENT(IN) :: ptr
         END FUNCTION
      END INTERFACE

!  Modifies a parameter value in a Data Fitting task descriptor
      INTERFACE
         INTEGER FUNCTION dfieditval(task,val_attr,val)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: val_attr
            INTEGER,INTENT(IN)                   :: val
         END FUNCTION
      END INTERFACE

!  Modifies a pointer to the memory representing a coordinate of the data
!  stored in matrix format (function or spline coefficients)
      INTERFACE
         INTEGER FUNCTION dfseditidxptr(task,ptr_attr,idx,ptr)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: ptr_attr
            INTEGER,INTENT(IN)                   :: idx
            REAL(KIND=4),DIMENSION(*),INTENT(IN) :: ptr
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdeditidxptr(task,ptr_attr,idx,ptr)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: ptr_attr
            INTEGER,INTENT(IN)                   :: idx
            REAL(KIND=8),DIMENSION(*),INTENT(IN) :: ptr
         END FUNCTION
      END INTERFACE

!  Modifies parameters of Piece-wise Polynomial (PP) spline
      INTERFACE
         INTEGER FUNCTION dfseditppspline1d(task,s_order,s_type,         &
     &                  bc_type,bc,ic_type,ic,scoeff,scoeffhint)
                USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: s_order
            INTEGER,INTENT(IN)                            :: s_type
            INTEGER,INTENT(IN),OPTIONAL                   :: bc_type
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: bc
            INTEGER,INTENT(IN),OPTIONAL                   :: ic_type
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: ic
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: scoeff
            INTEGER,INTENT(IN),OPTIONAL                   :: scoeffhint
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdeditppspline1d(task,s_order,s_type,         &
     &                  bc_type,bc,ic_type,ic,scoeff,scoeffhint)
                USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: s_order
            INTEGER,INTENT(IN)                            :: s_type
            INTEGER,INTENT(IN),OPTIONAL                   :: bc_type
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: bc
            INTEGER,INTENT(IN),OPTIONAL                   :: ic_type
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: ic
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: scoeff
            INTEGER,INTENT(IN),OPTIONAL                   :: scoeffhint
         END FUNCTION
      END INTERFACE

!++
!  DF TASK QUERYING FUNCTION DECLARATIONS.
!--

!  Reads a pointer to an array held in a Data Fitting task descriptor
      INTERFACE
         INTEGER FUNCTION dfsqueryptr(task,ptr_attr,ptr)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: ptr_attr
            INTEGER(KIND=8),INTENT(OUT)          :: ptr
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdqueryptr(task,ptr_attr,ptr)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: ptr_attr
            INTEGER(KIND=8),INTENT(OUT)          :: ptr
         END FUNCTION
      END INTERFACE

!  Reads a parameter value in a Data Fitting task descriptor
      INTERFACE
         INTEGER FUNCTION dfiqueryval(task,val_attr,val)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: val_attr
            INTEGER,INTENT(OUT)                  :: val
         END FUNCTION
      END INTERFACE

!  Reads a pointer to the memory representing a coordinate of the data
!  stored in matrix format (function or spline coefficients)
      INTERFACE
         INTEGER FUNCTION dfsqueryidxptr(task,ptr_attr,idx,ptr)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: ptr_attr
            INTEGER,INTENT(IN)                   :: idx
            INTEGER(KIND=8),INTENT(OUT)          :: ptr
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdqueryidxptr(task,ptr_attr,idx,ptr)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                        :: task
            INTEGER,INTENT(IN)                   :: ptr_attr
            INTEGER,INTENT(IN)                   :: idx
            INTEGER(KIND=8),INTENT(OUT)          :: ptr
         END FUNCTION
      END INTERFACE

!++
!  DF COMPUTE FUNCTION DECLARATIONS.
!--

      INTERFACE
         INTEGER FUNCTION dfsconstruct1d(task,s_format,method)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)               :: task
            INTEGER,INTENT(IN)          :: s_format
            INTEGER,INTENT(IN)          :: method
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdconstruct1d(task,s_format,method)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)               :: task
            INTEGER,INTENT(IN)          :: s_format
            INTEGER,INTENT(IN)          :: method
         END FUNCTION
      END INTERFACE

!  Interpolate1d

      INTERFACE
         INTEGER FUNCTION dfsinterpolate1d(task, type, method, nsite,    &
     &    site, sitehint, ndorder, dorder, datahint, r, rhint, cell)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: type
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nsite
            REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: site
            INTEGER,INTENT(IN),OPTIONAL                   :: sitehint
            INTEGER,INTENT(IN),OPTIONAL                   :: ndorder
            INTEGER, DIMENSION(*),INTENT(IN),OPTIONAL     :: dorder
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: datahint
            REAL(KIND=4),DIMENSION(*),INTENT(OUT)         :: r
            INTEGER,INTENT(IN),OPTIONAL                   :: rhint
            INTEGER,DIMENSION(*),INTENT(OUT),OPTIONAL     :: cell
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdinterpolate1d(task, type, method, nsite,    &
     &    site, sitehint, ndorder, dorder, datahint, r, rhint, cell)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: type
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nsite
            REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: site
            INTEGER,INTENT(IN),OPTIONAL                   :: sitehint
            INTEGER,INTENT(IN),OPTIONAL                   :: ndorder
            INTEGER, DIMENSION(*),INTENT(IN),OPTIONAL     :: dorder
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: datahint
            REAL(KIND=8),DIMENSION(*),INTENT(OUT)         :: r
            INTEGER,INTENT(IN),OPTIONAL                   :: rhint
            INTEGER,DIMENSION(*),INTENT(OUT),OPTIONAL     :: cell
         END FUNCTION
      END INTERFACE

!  InterpolateEx1d

      INTERFACE
         INTEGER FUNCTION dfsinterpolateex1d(task, type, method, nsite,  &
     &    site, sitehint, ndorder, dorder, datahint, r, rhint, cell,     &
     &    le_cb, le_params, re_cb, re_params, i_cb, i_params, search_cb, &
     &    search_params)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: type
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nsite
            REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: site
            INTEGER,INTENT(IN),OPTIONAL                   :: sitehint
            INTEGER,INTENT(IN),OPTIONAL                   :: ndorder
            INTEGER, DIMENSION(*),INTENT(IN),OPTIONAL     :: dorder
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: datahint
            REAL(KIND=4),DIMENSION(*),INTENT(OUT)         :: r
            INTEGER,INTENT(IN),OPTIONAL                   :: rhint
            INTEGER,DIMENSION(*),INTENT(OUT),OPTIONAL     :: cell
            INTEGER,EXTERNAL,OPTIONAL                     :: le_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: le_params
            INTEGER,EXTERNAL,OPTIONAL                     :: re_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: re_params
            INTEGER,EXTERNAL,OPTIONAL                     :: i_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: i_params
            INTEGER,EXTERNAL,OPTIONAL                   :: search_cb
            INTEGER,DIMENSION(*),OPTIONAL               :: search_params
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdinterpolateex1d(task, type, method, nsite,  &
     &    site, sitehint, ndorder, dorder, datahint, r, rhint, cell,     &
     &    le_cb, le_params, re_cb, re_params, i_cb, i_params, search_cb, &
     &    search_params)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: type
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nsite
            REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: site
            INTEGER,INTENT(IN),OPTIONAL                   :: sitehint
            INTEGER,INTENT(IN),OPTIONAL                   :: ndorder
            INTEGER, DIMENSION(*),INTENT(IN),OPTIONAL     :: dorder
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: datahint
            REAL(KIND=8),DIMENSION(*),INTENT(OUT)         :: r
            INTEGER,INTENT(IN),OPTIONAL                   :: rhint
            INTEGER,DIMENSION(*),INTENT(OUT),OPTIONAL     :: cell
            INTEGER,EXTERNAL,OPTIONAL                     :: le_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: le_params
            INTEGER,EXTERNAL,OPTIONAL                     :: re_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: re_params
            INTEGER,EXTERNAL,OPTIONAL                     :: i_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: i_params
            INTEGER,EXTERNAL,OPTIONAL                   :: search_cb
            INTEGER,DIMENSION(*),OPTIONAL               :: search_params
         END FUNCTION
      END INTERFACE

!  Integrate1d

      INTERFACE
         INTEGER FUNCTION dfsintegrate1d(task, method, nlim, llim,       &
     &    llimhint, rlim, rlimhint, ldatahint, rdatahint, r, rhint)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nlim
            REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: llim
            INTEGER,INTENT(IN),OPTIONAL                   :: llimhint
            REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: rlim
            INTEGER,INTENT(IN),OPTIONAL                   :: rlimhint
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: ldatahint
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: rdatahint
            REAL(KIND=4),DIMENSION(*),INTENT(OUT)         :: r
            INTEGER,INTENT(IN),OPTIONAL                   :: rhint
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdintegrate1d(task, method, nlim, llim,       &
     &    llimhint, rlim, rlimhint, ldatahint, rdatahint, r, rhint)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nlim
            REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: llim
            INTEGER,INTENT(IN),OPTIONAL                   :: llimhint
            REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: rlim
            INTEGER,INTENT(IN),OPTIONAL                   :: rlimhint
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: ldatahint
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: rdatahint
            REAL(KIND=8),DIMENSION(*),INTENT(OUT)         :: r
            INTEGER,INTENT(IN),OPTIONAL                   :: rhint
         END FUNCTION
      END INTERFACE

!  IntegrateEx1d

      INTERFACE
         INTEGER FUNCTION dfsintegrateex1d(task, method, nlim, llim,     &
     &    llimhint, rlim, rlimhint, ldatahint, rdatahint, r, rhint,      &
     &    le_cb, le_params, re_cb, re_params, i_cb, i_params, search_cb, &
     &    search_params)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nlim
            REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: llim
            INTEGER,INTENT(IN),OPTIONAL                   :: llimhint
            REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: rlim
            INTEGER,INTENT(IN),OPTIONAL                   :: rlimhint
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: ldatahint
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: rdatahint
            REAL(KIND=4),DIMENSION(*),INTENT(OUT)         :: r
            INTEGER,INTENT(IN),OPTIONAL                   :: rhint
            INTEGER,EXTERNAL,OPTIONAL                     :: le_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: le_params
            INTEGER,EXTERNAL,OPTIONAL                     :: re_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: re_params
            INTEGER,EXTERNAL,OPTIONAL                     :: i_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: i_params
            INTEGER,EXTERNAL,OPTIONAL                   :: search_cb
            INTEGER,DIMENSION(*),OPTIONAL               :: search_params
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdintegrateex1d(task, method, nlim, llim,     &
     &    llimhint, rlim, rlimhint, ldatahint, rdatahint, r, rhint,      &
     &    le_cb, le_params, re_cb, re_params, i_cb, i_params, search_cb, &
     &    search_params)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nlim
            REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: llim
            INTEGER,INTENT(IN),OPTIONAL                   :: llimhint
            REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: rlim
            INTEGER,INTENT(IN),OPTIONAL                   :: rlimhint
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: ldatahint
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: rdatahint
            REAL(KIND=8),DIMENSION(*),INTENT(OUT)         :: r
            INTEGER,INTENT(IN),OPTIONAL                   :: rhint
            INTEGER,EXTERNAL,OPTIONAL                     :: le_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: le_params
            INTEGER,EXTERNAL,OPTIONAL                     :: re_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: re_params
            INTEGER,EXTERNAL,OPTIONAL                     :: i_cb
            INTEGER,DIMENSION(*),OPTIONAL                 :: i_params
            INTEGER,EXTERNAL,OPTIONAL                   :: search_cb
            INTEGER,DIMENSION(*),OPTIONAL               :: search_params
         END FUNCTION
      END INTERFACE

!  SearchCells1d

      INTERFACE
         INTEGER FUNCTION dfssearchcells1d(task, method, nsite, site,    &
     &                                     sitehint, datahint, cell)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nsite
            REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: site
            INTEGER,INTENT(IN),OPTIONAL                   :: sitehint
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: datahint
            INTEGER,DIMENSION(*),INTENT(OUT)              :: cell
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdsearchcells1d(task, method, nsite, site,    &
     &                                     sitehint, datahint, cell)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nsite
            REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: site
            INTEGER,INTENT(IN),OPTIONAL                   :: sitehint
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: datahint
            INTEGER,DIMENSION(*),INTENT(OUT)              :: cell
         END FUNCTION
      END INTERFACE

!  SearchCellsEx1d

      INTERFACE
         INTEGER FUNCTION dfssearchcellsex1d(task, method, nsite, site,    &
     &              sitehint, datahint, cell, search_cb, search_params)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nsite
            REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: site
            INTEGER,INTENT(IN),OPTIONAL                   :: sitehint
            REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: datahint
            INTEGER,DIMENSION(*),INTENT(OUT)              :: cell
            INTEGER,EXTERNAL,OPTIONAL                   :: search_cb
            INTEGER,DIMENSION(*),OPTIONAL               :: search_params
         END FUNCTION
      END INTERFACE

      INTERFACE
         INTEGER FUNCTION dfdsearchcellsex1d(task, method, nsite, site,    &
     &              sitehint, datahint, cell, search_cb, search_params)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)                                 :: task
            INTEGER,INTENT(IN)                            :: method
            INTEGER,INTENT(IN)                            :: nsite
            REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: site
            INTEGER,INTENT(IN),OPTIONAL                   :: sitehint
            REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: datahint
            INTEGER,DIMENSION(*),INTENT(OUT)              :: cell
            INTEGER,EXTERNAL,OPTIONAL                   :: search_cb
            INTEGER,DIMENSION(*),OPTIONAL               :: search_params
         END FUNCTION
      END INTERFACE

!++
!  DF DESTRUCTOR FUNCTION DECLARATIONS.
!--

      INTERFACE
         INTEGER FUNCTION dfdeletetask(task)
               USE MKL_DF_TYPE
            TYPE(DF_TASK)   :: task
         END FUNCTION
      END INTERFACE

      END MODULE MKL_DF
