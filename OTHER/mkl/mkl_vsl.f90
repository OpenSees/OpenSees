! file: mkl_vsl.fi
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
!  Fortran 90 VSL interface.
!--

      MODULE MKL_VSL_TYPE


!++
!  Definitions for VSL functions return values (errors, warnings)
!--

!  "No error" status
      INTEGER(KIND=4) VSL_STATUS_OK
      INTEGER(KIND=4) VSL_ERROR_OK
      PARAMETER (VSL_STATUS_OK = 0)
      PARAMETER (VSL_ERROR_OK  = 0)

!  Common errors (-1..-999)
      INTEGER(KIND=4) VSL_ERROR_FEATURE_NOT_IMPLEMENTED
      INTEGER(KIND=4) VSL_ERROR_UNKNOWN
      INTEGER(KIND=4) VSL_ERROR_BADARGS
      INTEGER(KIND=4) VSL_ERROR_MEM_FAILURE
      INTEGER(KIND=4) VSL_ERROR_NULL_PTR
      INTEGER(KIND=4) VSL_ERROR_CPU_NOT_SUPPORTED
      PARAMETER (VSL_ERROR_FEATURE_NOT_IMPLEMENTED = -1)
      PARAMETER (VSL_ERROR_UNKNOWN                 = -2)
      PARAMETER (VSL_ERROR_BADARGS                 = -3)
      PARAMETER (VSL_ERROR_MEM_FAILURE             = -4)
      PARAMETER (VSL_ERROR_NULL_PTR                = -5)
      PARAMETER (VSL_ERROR_CPU_NOT_SUPPORTED       = -6)

!  RNG errors (-1000..-1999)
!  brng errors
      INTEGER(KIND=4) VSL_RNG_ERROR_INVALID_BRNG_INDEX
      INTEGER(KIND=4) VSL_RNG_ERROR_LEAPFROG_UNSUPPORTED
      INTEGER(KIND=4) VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED
      INTEGER(KIND=4) VSL_RNG_ERROR_BRNGS_INCOMPATIBLE
      INTEGER(KIND=4) VSL_RNG_ERROR_BAD_STREAM
      INTEGER(KIND=4) VSL_RNG_ERROR_BRNG_TABLE_FULL
      INTEGER(KIND=4) VSL_RNG_ERROR_BAD_STREAM_STATE_SIZE
      INTEGER(KIND=4) VSL_RNG_ERROR_BAD_WORD_SIZE
      INTEGER(KIND=4) VSL_RNG_ERROR_BAD_NSEEDS
      INTEGER(KIND=4) VSL_RNG_ERROR_BAD_NBITS
      INTEGER(KIND=4) VSL_RNG_ERROR_QRNG_PERIOD_ELAPSED
      INTEGER(KIND=4) VSL_RNG_ERROR_LEAPFROG_NSTREAMS_TOO_BIG
      PARAMETER (VSL_RNG_ERROR_INVALID_BRNG_INDEX        = -1000)
      PARAMETER (VSL_RNG_ERROR_LEAPFROG_UNSUPPORTED      = -1002)
      PARAMETER (VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED     = -1003)
      PARAMETER (VSL_RNG_ERROR_BRNGS_INCOMPATIBLE        = -1005)
      PARAMETER (VSL_RNG_ERROR_BAD_STREAM                = -1006)
      PARAMETER (VSL_RNG_ERROR_BRNG_TABLE_FULL           = -1007)
      PARAMETER (VSL_RNG_ERROR_BAD_STREAM_STATE_SIZE     = -1008)
      PARAMETER (VSL_RNG_ERROR_BAD_WORD_SIZE             = -1009)
      PARAMETER (VSL_RNG_ERROR_BAD_NSEEDS                = -1010)
      PARAMETER (VSL_RNG_ERROR_BAD_NBITS                 = -1011)
      PARAMETER (VSL_RNG_ERROR_QRNG_PERIOD_ELAPSED       = -1012)
      PARAMETER (VSL_RNG_ERROR_LEAPFROG_NSTREAMS_TOO_BIG = -1013)

! abstract stream related errors
      INTEGER(KIND=4) VSL_RNG_ERROR_BAD_UPDATE
      INTEGER(KIND=4) VSL_RNG_ERROR_NO_NUMBERS
      INTEGER(KIND=4) VSL_RNG_ERROR_INVALID_ABSTRACT_STREAM
      PARAMETER (VSL_RNG_ERROR_BAD_UPDATE              = -1120)
      PARAMETER (VSL_RNG_ERROR_NO_NUMBERS              = -1121)
      PARAMETER (VSL_RNG_ERROR_INVALID_ABSTRACT_STREAM = -1122)

! non determenistic stream related errors
      INTEGER(KIND=4) VSL_RNG_ERROR_NONDETERM_NOT_SUPPORTED
      INTEGER(KIND=4) VSL_RNG_ERROR_NONDETERM_NRETRIES_EXCEEDED
      PARAMETER (VSL_RNG_ERROR_NONDETERM_NOT_SUPPORTED     = -1130)
      PARAMETER (VSL_RNG_ERROR_NONDETERM_NRETRIES_EXCEEDED = -1131)

! ARS5 stream related errors
      INTEGER(KIND=4) VSL_RNG_ERROR_ARS5_NOT_SUPPORTED
      PARAMETER (VSL_RNG_ERROR_ARS5_NOT_SUPPORTED   = -1140)

! read/write stream to file errors
      INTEGER(KIND=4) VSL_RNG_ERROR_FILE_CLOSE
      INTEGER(KIND=4) VSL_RNG_ERROR_FILE_OPEN
      INTEGER(KIND=4) VSL_RNG_ERROR_FILE_WRITE
      INTEGER(KIND=4) VSL_RNG_ERROR_FILE_READ

      INTEGER(KIND=4) VSL_RNG_ERROR_BAD_FILE_FORMAT
      INTEGER(KIND=4) VSL_RNG_ERROR_UNSUPPORTED_FILE_VER
      INTEGER(KIND=4) VSL_RNG_ERROR_BAD_MEM_FORMAT
      PARAMETER (VSL_RNG_ERROR_FILE_CLOSE           = -1100)
      PARAMETER (VSL_RNG_ERROR_FILE_OPEN            = -1101)
      PARAMETER (VSL_RNG_ERROR_FILE_WRITE           = -1102)
      PARAMETER (VSL_RNG_ERROR_FILE_READ            = -1103)

      PARAMETER (VSL_RNG_ERROR_BAD_FILE_FORMAT      = -1110)
      PARAMETER (VSL_RNG_ERROR_UNSUPPORTED_FILE_VER = -1111)
      PARAMETER (VSL_RNG_ERROR_BAD_MEM_FORMAT       = -1200)

! Convolution/correlation errors
      INTEGER(KIND=4) VSL_CC_ERROR_NOT_IMPLEMENTED
      INTEGER(KIND=4) VSL_CC_ERROR_ALLOCATION_FAILURE
      INTEGER(KIND=4) VSL_CC_ERROR_BAD_DESCRIPTOR
      INTEGER(KIND=4) VSL_CC_ERROR_SERVICE_FAILURE
      INTEGER(KIND=4) VSL_CC_ERROR_EDIT_FAILURE
      INTEGER(KIND=4) VSL_CC_ERROR_EDIT_PROHIBITED
      INTEGER(KIND=4) VSL_CC_ERROR_COMMIT_FAILURE
      INTEGER(KIND=4) VSL_CC_ERROR_COPY_FAILURE
      INTEGER(KIND=4) VSL_CC_ERROR_DELETE_FAILURE
      INTEGER(KIND=4) VSL_CC_ERROR_BAD_ARGUMENT
      INTEGER(KIND=4) VSL_CC_ERROR_DIMS
      INTEGER(KIND=4) VSL_CC_ERROR_START
      INTEGER(KIND=4) VSL_CC_ERROR_DECIMATION
      INTEGER(KIND=4) VSL_CC_ERROR_XSHAPE
      INTEGER(KIND=4) VSL_CC_ERROR_YSHAPE
      INTEGER(KIND=4) VSL_CC_ERROR_ZSHAPE
      INTEGER(KIND=4) VSL_CC_ERROR_XSTRIDE
      INTEGER(KIND=4) VSL_CC_ERROR_YSTRIDE
      INTEGER(KIND=4) VSL_CC_ERROR_ZSTRIDE
      INTEGER(KIND=4) VSL_CC_ERROR_X
      INTEGER(KIND=4) VSL_CC_ERROR_Y
      INTEGER(KIND=4) VSL_CC_ERROR_Z
      INTEGER(KIND=4) VSL_CC_ERROR_JOB
      INTEGER(KIND=4) VSL_CC_ERROR_KIND
      INTEGER(KIND=4) VSL_CC_ERROR_MODE
      INTEGER(KIND=4) VSL_CC_ERROR_TYPE
      INTEGER(KIND=4) VSL_CC_ERROR_PRECISION
      INTEGER(KIND=4) VSL_CC_ERROR_EXTERNAL_PRECISION
      INTEGER(KIND=4) VSL_CC_ERROR_INTERNAL_PRECISION
      INTEGER(KIND=4) VSL_CC_ERROR_METHOD
      INTEGER(KIND=4) VSL_CC_ERROR_OTHER
      PARAMETER (VSL_CC_ERROR_NOT_IMPLEMENTED    = -2000)
      PARAMETER (VSL_CC_ERROR_ALLOCATION_FAILURE = -2001)
      PARAMETER (VSL_CC_ERROR_BAD_DESCRIPTOR     = -2200)
      PARAMETER (VSL_CC_ERROR_SERVICE_FAILURE    = -2210)
      PARAMETER (VSL_CC_ERROR_EDIT_FAILURE       = -2211)
      PARAMETER (VSL_CC_ERROR_EDIT_PROHIBITED    = -2212)
      PARAMETER (VSL_CC_ERROR_COMMIT_FAILURE     = -2220)
      PARAMETER (VSL_CC_ERROR_COPY_FAILURE       = -2230)
      PARAMETER (VSL_CC_ERROR_DELETE_FAILURE     = -2240)
      PARAMETER (VSL_CC_ERROR_BAD_ARGUMENT       = -2300)
      PARAMETER (VSL_CC_ERROR_DIMS               = -2301)
      PARAMETER (VSL_CC_ERROR_START              = -2302)
      PARAMETER (VSL_CC_ERROR_DECIMATION         = -2303)
      PARAMETER (VSL_CC_ERROR_XSHAPE             = -2311)
      PARAMETER (VSL_CC_ERROR_YSHAPE             = -2312)
      PARAMETER (VSL_CC_ERROR_ZSHAPE             = -2313)
      PARAMETER (VSL_CC_ERROR_XSTRIDE            = -2321)
      PARAMETER (VSL_CC_ERROR_YSTRIDE            = -2322)
      PARAMETER (VSL_CC_ERROR_ZSTRIDE            = -2323)
      PARAMETER (VSL_CC_ERROR_X                  = -2331)
      PARAMETER (VSL_CC_ERROR_Y                  = -2332)
      PARAMETER (VSL_CC_ERROR_Z                  = -2333)
      PARAMETER (VSL_CC_ERROR_JOB                = -2100)
      PARAMETER (VSL_CC_ERROR_KIND               = -2110)
      PARAMETER (VSL_CC_ERROR_MODE               = -2120)
      PARAMETER (VSL_CC_ERROR_TYPE               = -2130)
      PARAMETER (VSL_CC_ERROR_PRECISION          = -2400)
      PARAMETER (VSL_CC_ERROR_EXTERNAL_PRECISION = -2141)
      PARAMETER (VSL_CC_ERROR_INTERNAL_PRECISION = -2142)
      PARAMETER (VSL_CC_ERROR_METHOD             = -2400)
      PARAMETER (VSL_CC_ERROR_OTHER              = -2800)

!++
! SUMMARY STATTISTICS ERROR/WARNING CODES
!--

!  Warnings
      INTEGER(KIND=4) VSL_SS_NOT_FULL_RANK_MATRIX
      INTEGER(KIND=4) VSL_SS_SEMIDEFINITE_COR
      PARAMETER (VSL_SS_NOT_FULL_RANK_MATRIX                 = 4028)
      PARAMETER (VSL_SS_SEMIDEFINITE_COR                     = 4029)

!  Errors and messages (-4000..-4999)
      INTEGER(KIND=4) VSL_SS_ERROR_ALLOCATION_FAILURE
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_DIMEN
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_OBSERV_N
      INTEGER(KIND=4) VSL_SS_ERROR_STORAGE_NOT_SUPPORTED
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_INDC_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_WEIGHTS
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MEAN_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_2R_MOM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_3R_MOM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_4R_MOM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_2C_MOM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_3C_MOM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_4C_MOM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_KURTOSIS_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_SKEWNESS_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MIN_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MAX_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_VARIATION_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_COV_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_COR_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_ACCUM_WEIGHT_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_QUANT_ORDER_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_QUANT_ORDER
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_QUANT_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_ORDER_STATS_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_MOMORDER_NOT_SUPPORTED
      INTEGER(KIND=4) VSL_SS_ERROR_ALL_OBSERVS_OUTLIERS
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_ROBUST_COV_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_ROBUST_MEAN_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_METHOD_NOT_SUPPORTED
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_GROUP_INDC_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_NULL_TASK_DESCRIPTOR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_OBSERV_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_SINGULAR_COV
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_POOLED_COV_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_POOLED_MEAN_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_GROUP_COV_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_GROUP_MEAN_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_GROUP_INDC
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_N_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_OUTLIERS_WEIGHTS_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_N_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_STORAGE_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_PARTIAL_COV_IDX_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_PARTIAL_COV_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_PARTIAL_COR_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_PARAMS_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_PARAMS_N_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_BAD_PARAMS_N
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_PARAMS
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_N_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_SIMUL_VALS_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_ESTIMATES_N_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_ESTIMATES_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_ESTIMATES_N
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_OUTPUT_PARAMS
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_PRIOR_N_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_PRIOR_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MI_MISSING_VALS_N
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_STREAM_QUANT_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_PARAMTR_COR_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_COR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_PARTIAL_COV_IDX
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_SUM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_2R_SUM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_3R_SUM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_4R_SUM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_2C_SUM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_3C_SUM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_4C_SUM_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_CP_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MDAD_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_MNAD_ADDR
      INTEGER(KIND=4) VSL_SS_ERROR_INDICES_NOT_SUPPORTED
      INTEGER(KIND=4) VSL_SS_ERROR_BAD_SORTED_OBSERV_ADDR

      PARAMETER (VSL_SS_ERROR_ALLOCATION_FAILURE             =-4000)
      PARAMETER (VSL_SS_ERROR_BAD_DIMEN                      =-4001)
      PARAMETER (VSL_SS_ERROR_BAD_OBSERV_N                   =-4002)
      PARAMETER (VSL_SS_ERROR_STORAGE_NOT_SUPPORTED          =-4003)
      PARAMETER (VSL_SS_ERROR_BAD_INDC_ADDR                  =-4004)
      PARAMETER (VSL_SS_ERROR_BAD_WEIGHTS                    =-4005)
      PARAMETER (VSL_SS_ERROR_BAD_MEAN_ADDR                  =-4006)
      PARAMETER (VSL_SS_ERROR_BAD_2R_MOM_ADDR                =-4007)
      PARAMETER (VSL_SS_ERROR_BAD_3R_MOM_ADDR                =-4008)
      PARAMETER (VSL_SS_ERROR_BAD_4R_MOM_ADDR                =-4009)
      PARAMETER (VSL_SS_ERROR_BAD_2C_MOM_ADDR                =-4010)
      PARAMETER (VSL_SS_ERROR_BAD_3C_MOM_ADDR                =-4011)
      PARAMETER (VSL_SS_ERROR_BAD_4C_MOM_ADDR                =-4012)
      PARAMETER (VSL_SS_ERROR_BAD_KURTOSIS_ADDR              =-4013)
      PARAMETER (VSL_SS_ERROR_BAD_SKEWNESS_ADDR              =-4014)
      PARAMETER (VSL_SS_ERROR_BAD_MIN_ADDR                   =-4015)
      PARAMETER (VSL_SS_ERROR_BAD_MAX_ADDR                   =-4016)
      PARAMETER (VSL_SS_ERROR_BAD_VARIATION_ADDR             =-4017)
      PARAMETER (VSL_SS_ERROR_BAD_COV_ADDR                   =-4018)
      PARAMETER (VSL_SS_ERROR_BAD_COR_ADDR                   =-4019)
      PARAMETER (VSL_SS_ERROR_BAD_ACCUM_WEIGHT_ADDR          =-4020)
      PARAMETER (VSL_SS_ERROR_BAD_QUANT_ORDER_ADDR           =-4021)
      PARAMETER (VSL_SS_ERROR_BAD_QUANT_ORDER                =-4022)
      PARAMETER (VSL_SS_ERROR_BAD_QUANT_ADDR                 =-4023)
      PARAMETER (VSL_SS_ERROR_BAD_ORDER_STATS_ADDR           =-4024)
      PARAMETER (VSL_SS_ERROR_MOMORDER_NOT_SUPPORTED         =-4025)
      PARAMETER (VSL_SS_ERROR_ALL_OBSERVS_OUTLIERS           =-4026)
      PARAMETER (VSL_SS_ERROR_BAD_ROBUST_COV_ADDR            =-4027)
      PARAMETER (VSL_SS_ERROR_BAD_ROBUST_MEAN_ADDR           =-4028)
      PARAMETER (VSL_SS_ERROR_METHOD_NOT_SUPPORTED           =-4029)
      PARAMETER (VSL_SS_ERROR_BAD_GROUP_INDC_ADDR            =-4030)
      PARAMETER (VSL_SS_ERROR_NULL_TASK_DESCRIPTOR           =-4031)
      PARAMETER (VSL_SS_ERROR_BAD_OBSERV_ADDR                =-4032)
      PARAMETER (VSL_SS_ERROR_SINGULAR_COV                   =-4033)
      PARAMETER (VSL_SS_ERROR_BAD_POOLED_COV_ADDR            =-4034)
      PARAMETER (VSL_SS_ERROR_BAD_POOLED_MEAN_ADDR           =-4035)
      PARAMETER (VSL_SS_ERROR_BAD_GROUP_COV_ADDR             =-4036)
      PARAMETER (VSL_SS_ERROR_BAD_GROUP_MEAN_ADDR            =-4037)
      PARAMETER (VSL_SS_ERROR_BAD_GROUP_INDC                 =-4038)
      PARAMETER (VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_ADDR       =-4039)
      PARAMETER (VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_N_ADDR     =-4040)
      PARAMETER (VSL_SS_ERROR_BAD_OUTLIERS_WEIGHTS_ADDR      =-4041)
      PARAMETER (VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_ADDR     =-4042)
      PARAMETER (VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_N_ADDR   =-4043)
      PARAMETER (VSL_SS_ERROR_BAD_STORAGE_ADDR               =-4044)
      PARAMETER (VSL_SS_ERROR_BAD_PARTIAL_COV_IDX_ADDR       =-4045)
      PARAMETER (VSL_SS_ERROR_BAD_PARTIAL_COV_ADDR           =-4046)
      PARAMETER (VSL_SS_ERROR_BAD_PARTIAL_COR_ADDR           =-4047)
      PARAMETER (VSL_SS_ERROR_BAD_MI_PARAMS_ADDR             =-4048)
      PARAMETER (VSL_SS_ERROR_BAD_MI_PARAMS_N_ADDR           =-4049)
      PARAMETER (VSL_SS_ERROR_BAD_MI_BAD_PARAMS_N            =-4050)
      PARAMETER (VSL_SS_ERROR_BAD_MI_PARAMS                  =-4051)
      PARAMETER (VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_N_ADDR   =-4052)
      PARAMETER (VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_ADDR     =-4053)
      PARAMETER (VSL_SS_ERROR_BAD_MI_SIMUL_VALS_ADDR         =-4054)
      PARAMETER (VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N_ADDR       =-4055)
      PARAMETER (VSL_SS_ERROR_BAD_MI_ESTIMATES_N_ADDR        =-4056)
      PARAMETER (VSL_SS_ERROR_BAD_MI_ESTIMATES_ADDR          =-4057)
      PARAMETER (VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N            =-4058)
      PARAMETER (VSL_SS_ERROR_BAD_MI_ESTIMATES_N             =-4059)
      PARAMETER (VSL_SS_ERROR_BAD_MI_OUTPUT_PARAMS           =-4060)
      PARAMETER (VSL_SS_ERROR_BAD_MI_PRIOR_N_ADDR            =-4061)
      PARAMETER (VSL_SS_ERROR_BAD_MI_PRIOR_ADDR              =-4062)
      PARAMETER (VSL_SS_ERROR_BAD_MI_MISSING_VALS_N          =-4063)
      PARAMETER (VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N_ADDR =-4064)
      PARAMETER (VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_ADDR   =-4065)
      PARAMETER (VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N      =-4066)
      PARAMETER (VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS        =-4067)
      PARAMETER (VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER_ADDR    =-4068)
      PARAMETER (VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER         =-4069)
      PARAMETER (VSL_SS_ERROR_BAD_STREAM_QUANT_ADDR          =-4070)
      PARAMETER (VSL_SS_ERROR_BAD_PARAMTR_COR_ADDR           =-4071)
      PARAMETER (VSL_SS_ERROR_BAD_COR                        =-4072)
      PARAMETER (VSL_SS_ERROR_BAD_PARTIAL_COV_IDX            =-4073)
      PARAMETER (VSL_SS_ERROR_BAD_SUM_ADDR                   =-4074)
      PARAMETER (VSL_SS_ERROR_BAD_2R_SUM_ADDR                =-4075)
      PARAMETER (VSL_SS_ERROR_BAD_3R_SUM_ADDR                =-4076)
      PARAMETER (VSL_SS_ERROR_BAD_4R_SUM_ADDR                =-4077)
      PARAMETER (VSL_SS_ERROR_BAD_2C_SUM_ADDR                =-4078)
      PARAMETER (VSL_SS_ERROR_BAD_3C_SUM_ADDR                =-4079)
      PARAMETER (VSL_SS_ERROR_BAD_4C_SUM_ADDR                =-4080)
      PARAMETER (VSL_SS_ERROR_BAD_CP_ADDR                    =-4081)
      PARAMETER (VSL_SS_ERROR_BAD_MDAD_ADDR                  =-4082)
      PARAMETER (VSL_SS_ERROR_BAD_MNAD_ADDR                  =-4083)
      PARAMETER (VSL_SS_ERROR_BAD_SORTED_OBSERV_ADDR         =-4084)
      PARAMETER (VSL_SS_ERROR_INDICES_NOT_SUPPORTED          =-4085)

!  Internal errors caused by internal routines of the functions
      INTEGER(KIND=4) VSL_SS_ERROR_ROBCOV_INTERN_C1
      INTEGER(KIND=4) VSL_SS_ERROR_PARTIALCOV_INTERN_C1
      INTEGER(KIND=4) VSL_SS_ERROR_PARTIALCOV_INTERN_C2
      INTEGER(KIND=4) VSL_SS_ERROR_MISSINGVALS_INTERN_C1
      INTEGER(KIND=4) VSL_SS_ERROR_MISSINGVALS_INTERN_C2
      INTEGER(KIND=4) VSL_SS_ERROR_MISSINGVALS_INTERN_C3
      INTEGER(KIND=4) VSL_SS_ERROR_MISSINGVALS_INTERN_C4
      INTEGER(KIND=4) VSL_SS_ERROR_MISSINGVALS_INTERN_C5
      INTEGER(KIND=4) VSL_SS_ERROR_PARAMTRCOR_INTERN_C1
      INTEGER(KIND=4) VSL_SS_ERROR_COVRANK_INTERNAL_ERROR_C1
      INTEGER(KIND=4) VSL_SS_ERROR_INVCOV_INTERNAL_ERROR_C1
      INTEGER(KIND=4) VSL_SS_ERROR_INVCOV_INTERNAL_ERROR_C2

      PARAMETER (VSL_SS_ERROR_ROBCOV_INTERN_C1               =-5000)
      PARAMETER (VSL_SS_ERROR_PARTIALCOV_INTERN_C1           =-5010)
      PARAMETER (VSL_SS_ERROR_PARTIALCOV_INTERN_C2           =-5011)
      PARAMETER (VSL_SS_ERROR_MISSINGVALS_INTERN_C1          =-5021)
      PARAMETER (VSL_SS_ERROR_MISSINGVALS_INTERN_C2          =-5022)
      PARAMETER (VSL_SS_ERROR_MISSINGVALS_INTERN_C3          =-5023)
      PARAMETER (VSL_SS_ERROR_MISSINGVALS_INTERN_C4          =-5024)
      PARAMETER (VSL_SS_ERROR_MISSINGVALS_INTERN_C5          =-5025)
      PARAMETER (VSL_SS_ERROR_PARAMTRCOR_INTERN_C1           =-5030)
      PARAMETER (VSL_SS_ERROR_COVRANK_INTERNAL_ERROR_C1      =-5040)
      PARAMETER (VSL_SS_ERROR_INVCOV_INTERNAL_ERROR_C1       =-5041)
      PARAMETER (VSL_SS_ERROR_INVCOV_INTERNAL_ERROR_C2       =-5042)

!++
!  CONV/CORR RELATED MACRO DEFINITIONS
!--
      INTEGER(KIND=4) VSL_CONV_MODE_AUTO
      INTEGER(KIND=4) VSL_CORR_MODE_AUTO
      INTEGER(KIND=4) VSL_CONV_MODE_DIRECT
      INTEGER(KIND=4) VSL_CORR_MODE_DIRECT
      INTEGER(KIND=4) VSL_CONV_MODE_FFT
      INTEGER(KIND=4) VSL_CORR_MODE_FFT
      INTEGER(KIND=4) VSL_CONV_PRECISION_SINGLE
      INTEGER(KIND=4) VSL_CORR_PRECISION_SINGLE
      INTEGER(KIND=4) VSL_CONV_PRECISION_DOUBLE
      INTEGER(KIND=4) VSL_CORR_PRECISION_DOUBLE
      PARAMETER (VSL_CONV_MODE_AUTO        = 0)
      PARAMETER (VSL_CORR_MODE_AUTO        = 0)
      PARAMETER (VSL_CONV_MODE_DIRECT      = 1)
      PARAMETER (VSL_CORR_MODE_DIRECT      = 1)
      PARAMETER (VSL_CONV_MODE_FFT         = 2)
      PARAMETER (VSL_CORR_MODE_FFT         = 2)
      PARAMETER (VSL_CONV_PRECISION_SINGLE = 1)
      PARAMETER (VSL_CORR_PRECISION_SINGLE = 1)
      PARAMETER (VSL_CONV_PRECISION_DOUBLE = 2)
      PARAMETER (VSL_CORR_PRECISION_DOUBLE = 2)


!++
!  BASIC RANDOM NUMBER GENERATOR (BRNG) RELATED MACRO DEFINITIONS
!--


!  MAX NUMBER OF BRNGS CAN BE REGISTERED IN VSL
!  No more than VSL_MAX_REG_BRNGS basic generators can be registered in VSL
!  (including predefined basic generators).
!
!  Change this number to increase/decrease number of BRNGs can be registered.
      INTEGER VSL_MAX_REG_BRNGS
      PARAMETER (VSL_MAX_REG_BRNGS = 512)

!  PREDEFINED BRNG NAMES
      INTEGER VSL_BRNG_SHIFT
      INTEGER VSL_BRNG_INC

      INTEGER VSL_BRNG_MCG31
      INTEGER VSL_BRNG_R250
      INTEGER VSL_BRNG_MRG32K3A
      INTEGER VSL_BRNG_MCG59
      INTEGER VSL_BRNG_WH
      INTEGER VSL_BRNG_SOBOL
      INTEGER VSL_BRNG_NIEDERR
      INTEGER VSL_BRNG_MT19937
      INTEGER VSL_BRNG_MT2203
      INTEGER VSL_BRNG_IABSTRACT
      INTEGER VSL_BRNG_DABSTRACT
      INTEGER VSL_BRNG_SABSTRACT
      INTEGER VSL_BRNG_SFMT19937
      INTEGER VSL_BRNG_NONDETERM
      INTEGER VSL_BRNG_ARS5
      INTEGER VSL_BRNG_PHILOX4X32X10

      PARAMETER (VSL_BRNG_SHIFT=20)
      PARAMETER (VSL_BRNG_INC=ISHFT(1, VSL_BRNG_SHIFT))

      PARAMETER (VSL_BRNG_MCG31        =VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_R250         =VSL_BRNG_MCG31    +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_MRG32K3A     =VSL_BRNG_R250     +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_MCG59        =VSL_BRNG_MRG32K3A +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_WH           =VSL_BRNG_MCG59    +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_SOBOL        =VSL_BRNG_WH       +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_NIEDERR      =VSL_BRNG_SOBOL    +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_MT19937      =VSL_BRNG_NIEDERR  +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_MT2203       =VSL_BRNG_MT19937  +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_IABSTRACT    =VSL_BRNG_MT2203   +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_DABSTRACT    =VSL_BRNG_IABSTRACT+VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_SABSTRACT    =VSL_BRNG_DABSTRACT+VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_SFMT19937    =VSL_BRNG_SABSTRACT+VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_NONDETERM    =VSL_BRNG_SFMT19937+VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_ARS5         =VSL_BRNG_NONDETERM+VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_PHILOX4X32X10=VSL_BRNG_ARS5     +VSL_BRNG_INC)

!  PREDEFINED PARAMETERS FOR NON-DETERMNINISTIC RANDOM NUMBER
!  GENERATOR
!  The library provides an abstraction to the source of
!  non-deterministic random numbers supported in HW.
!  Current version of the library provides interface to
!  RDRAND-based only, available in latest Intel CPU.

      INTEGER(KIND=4) VSL_BRNG_RDRAND
      INTEGER(KIND=4) VSL_BRNG_NONDETERM_NRETRIES

      PARAMETER (VSL_BRNG_RDRAND             =  0)
      PARAMETER (VSL_BRNG_NONDETERM_NRETRIES = 10)


!  LEAPFROG METHOD FOR GRAY-CODE BASED QUASI-RANDOM NUMBER BASIC GENERATORS
!  VSL_BRNG_SOBOL and VSL_BRNG_NIEDERR are Gray-code based quasi-random number
!  basic generators. In contrast to pseudorandom number basic generators,
!  quasi-random ones take the dimension as initialization parameter.
!
!  Suppose that quasi-random number generator (QRNG) dimension is S. QRNG
!  sequence is a sequence of S-dimensional vectors:
!
!     x0=(x0[0],x0[1],...,x0[S-1]),x1=(x1[0],x1[1],...,x1[S-1]),...
!
!  VSL treats the output of any basic generator as 1-dimensional, however:
!
!     x0[0],x0[1],...,x0[S-1],x1[0],x1[1],...,x1[S-1],...
!
!  Because of nature of VSL_BRNG_SOBOL and VSL_BRNG_NIEDERR QRNGs,
!  the only S-stride Leapfrog method is supported for them. In other words,
!  user can generate subsequences, which consist of fixed elements of
!  vectors x0,x1,... For example, if 0 element is fixed, the following
!  subsequence is generated:
!
!     x0[1],x1[1],x2[1],...
!
!  To use the s-stride Leapfrog method with given QRNG, user should call
!  vslLeapfrogStream function with parameter k equal to element to be fixed
!  (0<=k<S) and parameter nstreams equal to VSL_QRNG_LEAPFROG_COMPONENTS.
      INTEGER(KIND=4) VSL_QRNG_LEAPFROG_COMPONENTS
      PARAMETER (VSL_QRNG_LEAPFROG_COMPONENTS = Z"7FFFFFFF")


!  USER-DEFINED PARAMETERS FOR QUASI-RANDOM NUMBER BASIC GENERATORS
!  VSL_BRNG_SOBOL and VSL_BRNG_NIEDERR are Gray-code based quasi-random
!  number basic generators. Default parameters of the generators
!  support generation of quasi-random number vectors of dimensions
!  S<=40 for SOBOL and S<=318 for NIEDERRITER. The library provides
!  opportunity to register user-defined initial values for the
!  generators and generate quasi-random vectors of desirable dimension.
!  There is also opportunity to register user-defined parameters for
!  default dimensions and obtain another sequence of quasi-random vectors.
!  Service function vslNewStreamEx is used to pass the parameters to
!  the library. Data are packed into array params, parameter of the routine.
!  First element of the array is used for dimension S, second element
!  contains indicator, VSL_USER_QRNG_INITIAL_VALUES, of user-defined
!  parameters for quasi-random number generators.
!  Macros VSL_USER_PRIMITIVE_POLYMS and VSL_USER_INIT_DIRECTION_NUMBERS
!  are used to describe which data are passed to SOBOL QRNG and
!  VSL_USER_IRRED_POLYMS - which data are passed to NIEDERRITER QRNG.
!  For example, to demonstrate that both primitive polynomials and initial
!  direction numbers are passed in SOBOL one should set third element of the
!  array params to  VSL_USER_PRIMITIVE_POLYMS | VSL_USER_DIRECTION_NUMBERS.
!  Macro VSL_QRNG_OVERRIDE_1ST_DIM_INIT is used to override default
!  initialization for the first dimension. Macro VSL_USER_DIRECTION_NUMBERS
!  is used when direction numbers calculated on the user side are passed
!  into the generators. More detailed description of interface for
!  registration of user-defined QRNG initial parameters can be found
!  in VslNotes.pdf.
      INTEGER VSL_USER_QRNG_INITIAL_VALUES
      INTEGER VSL_USER_PRIMITIVE_POLYMS
      INTEGER VSL_USER_INIT_DIRECTION_NUMBERS
      INTEGER VSL_USER_IRRED_POLYMS
      INTEGER VSL_USER_DIRECTION_NUMBERS
      INTEGER VSL_QRNG_OVERRIDE_1ST_DIM_INIT

      PARAMETER (VSL_USER_QRNG_INITIAL_VALUES    = 1)
      PARAMETER (VSL_USER_PRIMITIVE_POLYMS       = 1)
      PARAMETER (VSL_USER_INIT_DIRECTION_NUMBERS = 2)
      PARAMETER (VSL_USER_IRRED_POLYMS           = 1)
      PARAMETER (VSL_USER_DIRECTION_NUMBERS      = 4)
      PARAMETER (VSL_QRNG_OVERRIDE_1ST_DIM_INIT  = 8)

!  INITIALIZATION METHODS FOR USER-DESIGNED BASIC RANDOM NUMBER GENERATORS.
!  Each BRNG must support at least VSL_INIT_METHOD_STANDARD initialization
!  method. In addition, VSL_INIT_METHOD_LEAPFROG and VSL_INIT_METHOD_SKIPAHEAD
!  initialization methods can be supported.
!
!  If VSL_INIT_METHOD_LEAPFROG is not supported then initialization routine
!  must return VSL_RNG_ERROR_LEAPFROG_UNSUPPORTED error code.
!
!  If VSL_INIT_METHOD_SKIPAHEAD is not supported then initialization routine
!  must return VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED error code.
!
!  If there is no error during initialization, the initialization routine must
!  return VSL_ERROR_OK code.
      INTEGER VSL_INIT_METHOD_STANDARD
      INTEGER VSL_INIT_METHOD_LEAPFROG
      INTEGER VSL_INIT_METHOD_SKIPAHEAD
      PARAMETER (VSL_INIT_METHOD_STANDARD  = 0)
      PARAMETER (VSL_INIT_METHOD_LEAPFROG  = 1)
      PARAMETER (VSL_INIT_METHOD_SKIPAHEAD = 2)

!++
!  ACCURACY FLAG FOR DISTRIBUTION GENERATORS
!  This flag defines mode of random number generation.
!  If accuracy mode is set distribution generators will produce
!  numbers lying exactly within definitional domain for all values
!  of distribution parameters. In this case slight performance
!  degradation is expected. By default accuracy mode is switched off
!  admitting random numbers to be out of the definitional domain for
!  specific values of distribution parameters.
!  This macro is used to form names for accuracy versions of
!  distribution number generators
!--

      INTEGER(KIND=4) VSL_RNG_METHOD_ACCURACY_FLAG
      PARAMETER (VSL_RNG_METHOD_ACCURACY_FLAG=ISHFT(1,30))

!++
!  TRANSFORMATION METHOD NAMES FOR DISTRIBUTION RANDOM NUMBER GENERATORS
!  VSL interface allows more than one generation method in a distribution
!  transformation subroutine. Following macro definitions are used to
!  specify generation method for given distribution generator.
!
!  Method name macro is constructed as
!
!     VSL_RNG_METHOD_<Distribution>_<Method>
!
!  where
!
!     <Distribution> - probability distribution
!     <Method> - method name
!
!  VSL_RNG_METHOD_<Distribution>_<Method> should be used with
!  vsl<precision>Rng<Distribution> function only, where
!
!     <precision> - s (single) or d (double)
!     <Distribution> - probability distribution
!--

! Uniform
!
! <Method>   <Short Description>
! STD        standard method. Currently there is only one method for this
!            distribution generator
      INTEGER VSL_RNG_METHOD_UNIFORM_STD
      PARAMETER (VSL_RNG_METHOD_UNIFORM_STD = 0)

      INTEGER VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
      PARAMETER (VSL_RNG_METHOD_UNIFORM_STD_ACCURATE=                   &
     & IOR(VSL_RNG_METHOD_UNIFORM_STD, VSL_RNG_METHOD_ACCURACY_FLAG))


! Uniform Bits
!
! <Method>   <Short Description>
! STD        standard method. Currently there is only one method for this
!            distribution generator
      INTEGER VSL_RNG_METHOD_UNIFORMBITS_STD
      PARAMETER (VSL_RNG_METHOD_UNIFORMBITS_STD = 0)

! Uniform Bits 32
!
! <Method>   <Short Description>
! STD        standard method. Currently there is only one method for this
!            distribution generator
      INTEGER VSL_RNG_METHOD_UNIFORMBITS32_STD
      PARAMETER (VSL_RNG_METHOD_UNIFORMBITS32_STD = 0)

! Uniform Bits 64
!
! <Method>   <Short Description>
! STD        standard method. Currently there is only one method for this
!            distribution generator
      INTEGER VSL_RNG_METHOD_UNIFORMBITS64_STD
      PARAMETER (VSL_RNG_METHOD_UNIFORMBITS64_STD = 0)

! Gaussian
!
! <Method>   <Short Description>
! BOXMULLER  generates normally distributed random number x thru the pair of
!            uniformly distributed numbers u1 and u2 according to the formula:
!
!               x=sqrt(-ln(u1))*sin(2*Pi*u2)
!
! BOXMULLER2 generates pair of normally distributed random numbers x1 and x2
!            thru the pair of uniformly dustributed numbers u1 and u2
!            according to the formula
!
!               x1=sqrt(-ln(u1))*sin(2*Pi*u2)
!               x2=sqrt(-ln(u1))*cos(2*Pi*u2)
!
!            NOTE: implementation correctly works with odd vector lengths
!
! ICDF       inverse cumulative distribution function method
      INTEGER VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
      INTEGER VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
      INTEGER VSL_RNG_METHOD_GAUSSIAN_ICDF
      PARAMETER (VSL_RNG_METHOD_GAUSSIAN_BOXMULLER  = 0)
      PARAMETER (VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2 = 1)
      PARAMETER (VSL_RNG_METHOD_GAUSSIAN_ICDF       = 2)

! GaussianMV - multivariate (correlated) normal
! Multivariate (correlated) normal random number generator is based on
! uncorrelated Gaussian random number generator (see vslsRngGaussian and
! vsldRngGaussian functions):
!
! <Method>   <Short Description>
! BOXMULLER  generates normally distributed random number x thru the pair of
!            uniformly distributed numbers u1 and u2 according to the formula:
!
!               x=sqrt(-ln(u1))*sin(2*Pi*u2)
!
! BOXMULLER2 generates pair of normally distributed random numbers x1 and x2
!            thru the pair of uniformly dustributed numbers u1 and u2
!            according to the formula
!
!               x1=sqrt(-ln(u1))*sin(2*Pi*u2)
!               x2=sqrt(-ln(u1))*cos(2*Pi*u2)
!
!            NOTE: implementation correctly works with odd vector lengths
!
! ICDF       inverse cumulative distribution function method
      INTEGER VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER
      INTEGER VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2
      INTEGER VSL_RNG_METHOD_GAUSSIANMV_ICDF
      PARAMETER (VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER  = 0)
      PARAMETER (VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2 = 1)
      PARAMETER (VSL_RNG_METHOD_GAUSSIANMV_ICDF       = 2)

! Exponential
!
! <Method>   <Short Description>
! ICDF       inverse cumulative distribution function method
      INTEGER VSL_RNG_METHOD_EXPONENTIAL_ICDF
      PARAMETER (VSL_RNG_METHOD_EXPONENTIAL_ICDF = 0)

      INTEGER VSL_RNG_METHOD_EXPONENTIAL_ICDF_ACCURATE
      PARAMETER (VSL_RNG_METHOD_EXPONENTIAL_ICDF_ACCURATE=              &
     &IOR(VSL_RNG_METHOD_EXPONENTIAL_ICDF,VSL_RNG_METHOD_ACCURACY_FLAG))

! Laplace
!
! <Method>   <Short Description>
! ICDF       inverse cumulative distribution function method
!
! ICDF - inverse cumulative distribution function method:
!
!           x=+/-ln(u) with probability 1/2,
!
!        where
!
!           x - random number with Laplace distribution,
!           u - uniformly distributed random number
      INTEGER VSL_RNG_METHOD_LAPLACE_ICDF
      PARAMETER (VSL_RNG_METHOD_LAPLACE_ICDF = 0)

! Weibull
!
! <Method>   <Short Description>
! ICDF       inverse cumulative distribution function method
      INTEGER VSL_RNG_METHOD_WEIBULL_ICDF
      PARAMETER (VSL_RNG_METHOD_WEIBULL_ICDF = 0)

      INTEGER VSL_RNG_METHOD_WEIBULL_ICDF_ACCURATE
      PARAMETER (VSL_RNG_METHOD_WEIBULL_ICDF_ACCURATE=                  &
     & IOR(VSL_RNG_METHOD_WEIBULL_ICDF, VSL_RNG_METHOD_ACCURACY_FLAG))



! Cauchy
!
! <Method>   <Short Description>
! ICDF       inverse cumulative distribution function method
      INTEGER VSL_RNG_METHOD_CAUCHY_ICDF
      PARAMETER (VSL_RNG_METHOD_CAUCHY_ICDF = 0)

! Rayleigh
!
! <Method>   <Short Description>
! ICDF       inverse cumulative distribution function method
      INTEGER VSL_RNG_METHOD_RAYLEIGH_ICDF
      PARAMETER (VSL_RNG_METHOD_RAYLEIGH_ICDF = 0)

      INTEGER VSL_RNG_METHOD_RAYLEIGH_ICDF_ACCURATE
      PARAMETER (VSL_RNG_METHOD_RAYLEIGH_ICDF_ACCURATE=                 &
     & IOR(VSL_RNG_METHOD_RAYLEIGH_ICDF, VSL_RNG_METHOD_ACCURACY_FLAG))


! Lognormal
!
! <Method>   <Short Description>
! BOXMULLER2      Box-Muller 2 algorithm based method
      INTEGER VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2
      INTEGER VSL_RNG_METHOD_LOGNORMAL_ICDF
      PARAMETER (VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2 = 0)
      PARAMETER (VSL_RNG_METHOD_LOGNORMAL_ICDF = 1)

      INTEGER VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2_ACCURATE
      INTEGER VSL_RNG_METHOD_LOGNORMAL_ICDF_ACCURATE
      PARAMETER (VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2_ACCURATE=          &
     & IOR(VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2,                         &
     &     VSL_RNG_METHOD_ACCURACY_FLAG))
      PARAMETER (VSL_RNG_METHOD_LOGNORMAL_ICDF_ACCURATE=                &
     & IOR(VSL_RNG_METHOD_LOGNORMAL_ICDF, VSL_RNG_METHOD_ACCURACY_FLAG))


! Gumbel
!
! <Method>   <Short Description>
! ICDF       inverse cumulative distribution function method
      INTEGER VSL_RNG_METHOD_GUMBEL_ICDF
      PARAMETER (VSL_RNG_METHOD_GUMBEL_ICDF = 0)

! Gamma
!
! <Method>     <Short Description>
! GNORM        nonlinear transformation of gaussian numbers
! alpha>1,     based on acceptance/rejection method with
!              squeezes
!
! alpha>=0.6,  rejection from the Weibull distribution
! alpha<1
!
! alpha<0.6,   transformation of exponential power distribution
!              (EPD), EPD random numbers are generated using
!              by means of acceptance/rejection technique
      INTEGER VSL_RNG_METHOD_GAMMA_GNORM
      PARAMETER (VSL_RNG_METHOD_GAMMA_GNORM = 0)

      INTEGER VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE
      PARAMETER (VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE=                   &
     & IOR(VSL_RNG_METHOD_GAMMA_GNORM, VSL_RNG_METHOD_ACCURACY_FLAG))


! Beta
!
! <Method>     <Short Description>
! CJA - stands for first letters of Cheng, Johnk, and Atkinson
! Cheng      - generation of beta random numbers of the second kind
! min(p,q)>1   based on acceptance/rejection technique and its
!              transformation to beta random numbers of the first kind;
!
! Johnk,     - if q + K*p^2+C<=0, K=0.852..., C=-0.956...
! Atkinson,    algorithm of Johnk: beta distributed random number
! max(p,q)<1   is generated as u1^(1/p) / (u1^(1/p)+u2^(1/q)),
!              if u1^(1/p)+u2^(1/q)<=1;
!              otherwise switching algorithm of Atkinson:
!              interval (0,1) is divided into two domains (0,t) and (t,1),
!              on each interval acceptance/rejection technique with
!              convenient majorizing function is used;
!
! Atkinson   - switching algorithm of Atkinson is used
! min(p,q)<1   (with another point t, see short description above);
! max(p,q)>1
!
! ICDF       - inverse cumulative distribution function method according
!              to formulas x=1-u^(1/q) for p = 1, and x = u^(1/p) for q=1,
!              where x is beta distributed random number,
!              u - uniformly distributed random number.
!              for p=q=1 beta distribution reduces to uniform distribution.

      INTEGER VSL_RNG_METHOD_BETA_CJA
      PARAMETER (VSL_RNG_METHOD_BETA_CJA = 0)

      INTEGER VSL_RNG_METHOD_BETA_CJA_ACCURATE
      PARAMETER (VSL_RNG_METHOD_BETA_CJA_ACCURATE=                      &
     & IOR(VSL_RNG_METHOD_BETA_CJA, VSL_RNG_METHOD_ACCURACY_FLAG))



! Bernoulli
!
! <Method>   <Short Description>
! ICDF       inverse cumulative distribution function method
      INTEGER VSL_RNG_METHOD_BERNOULLI_ICDF
      PARAMETER (VSL_RNG_METHOD_BERNOULLI_ICDF = 0)

! Geometric
!
! <Method>   <Short Description>
! ICDF       inverse cumulative distribution function method
      INTEGER VSL_RNG_METHOD_GEOMETRIC_ICDF
      PARAMETER (VSL_RNG_METHOD_GEOMETRIC_ICDF = 0)

! Binomial
!
! <Method>   <Short Description>
! BTPE       for ntrial*min(p,1-p)>30 acceptance/rejection method with
!            decomposition onto 4 regions:
!
!               * 2 parallelograms;
!               * triangle;
!               * left exponential tail;
!               * right exponential tail.
!
!            othewise table lookup method is used
      INTEGER VSL_RNG_METHOD_BINOMIAL_BTPE
      PARAMETER (VSL_RNG_METHOD_BINOMIAL_BTPE = 0)

! Hypergeometric
!
! <Method>   <Short Description>
! H2PE       if mode of distribution is large, acceptance/rejection method is
!            used with decomposition onto 3 regions:
!
!               * rectangular;
!               * left exponential tail;
!               * right exponential tail.
!
!            othewise table lookup method is used
      INTEGER VSL_RNG_METHOD_HYPERGEOMETRIC_H2PE
      PARAMETER (VSL_RNG_METHOD_HYPERGEOMETRIC_H2PE = 0)

! Poisson
!
! <Method>   <Short Description>
! PTPE       if lambda>=27, acceptance/rejection method is used with
!            decomposition onto 4 regions:
!
!               * 2 parallelograms;
!               * triangle;
!               * left exponential tail;
!               * right exponential tail.
!
!            othewise table lookup method is used
!
! POISNORM   for lambda>=1 method is based on Poisson inverse CDF
!            approximation by Gaussian inverse CDF; for lambda<1
!            table lookup method is used.
      INTEGER VSL_RNG_METHOD_POISSON_PTPE
      INTEGER VSL_RNG_METHOD_POISSON_POISNORM
      PARAMETER (VSL_RNG_METHOD_POISSON_PTPE     = 0)
      PARAMETER (VSL_RNG_METHOD_POISSON_POISNORM = 1)

! Poisson
!
! <Method>   <Short Description>
! POISNORM   for lambda>=1 method is based on Poisson inverse CDF
!            approximation by Gaussian inverse CDF; for lambda<1
!            ICDF method is used.
      INTEGER VSL_RNG_METHOD_POISSONV_POISNORM
      PARAMETER (VSL_RNG_METHOD_POISSONV_POISNORM = 0)

! Negbinomial
!
! <Method>   <Short Description>
! NBAR       if (a-1)*(1-p)/p>=100, acceptance/rejection method is used with
!            decomposition onto 5 regions:
!
!               * rectangular;
!               * 2 trapezoid;
!               * left exponential tail;
!               * right exponential tail.
!
!            othewise table lookup method is used.
      INTEGER VSL_RNG_METHOD_NEGBINOMIAL_NBAR
      PARAMETER (VSL_RNG_METHOD_NEGBINOMIAL_NBAR = 0)

!++
!  MATRIX STORAGE SCHEMES
!--

! Some multivariate random number generators, e.g. GaussianMV, operate
! with matrix parameters. To optimize matrix parameters usage VSL offers
! following matrix storage schemes. (See VSL documentation for more details).
!
! FULL     - whole matrix is stored
! PACKED   - lower/higher triangular matrix is packed in 1-dimensional array
! DIAGONAL - diagonal elements are packed in 1-dimensional array
      INTEGER VSL_MATRIX_STORAGE_FULL
      INTEGER VSL_MATRIX_STORAGE_PACKED
      INTEGER VSL_MATRIX_STORAGE_DIAGONAL
      PARAMETER (VSL_MATRIX_STORAGE_FULL     = 0)
      PARAMETER (VSL_MATRIX_STORAGE_PACKED   = 1)
      PARAMETER (VSL_MATRIX_STORAGE_DIAGONAL = 2)

!++
!  SUMMARY STATISTICS (SS) RELATED MACRO DEFINITIONS
!--


!++
!  MATRIX STORAGE SCHEMES
!--

!
! SS routines work with matrix parameters, e.g. matrix of observations,
! variance-covariance matrix. To optimize work with matrices the library
! provides the following storage matrix schemes.
!
!++
! Matrix of observations:
! ROWS    - observations of the random vector are stored in raws, that
!           is, i-th row of the matrix of observations contains values
!           of i-th component of the random vector
! COLS    - observations of the random vector are stored in columns that
!           is, i-th column of the matrix of observations contains values
!           of i-th component of the random vector
!--
      INTEGER VSL_SS_MATRIX_STORAGE_ROWS
      INTEGER VSL_SS_MATRIX_STORAGE_COLS
      PARAMETER (VSL_SS_MATRIX_STORAGE_ROWS    = Z"00010000")
      PARAMETER (VSL_SS_MATRIX_STORAGE_COLS    = Z"00020000")


!++
! Variance-covariance/correlation matrix:
! FULL     - whole matrix is stored
! L_PACKED - lower triangular matrix is stored as 1-dimensional array
! U_PACKED - upper triangular matrix is stored as 1-dimensional array
!--
      INTEGER VSL_SS_MATRIX_STORAGE_FULL
      INTEGER VSL_SS_MATRIX_STORAGE_L_PACKED
      INTEGER VSL_SS_MATRIX_STORAGE_U_PACKED
      PARAMETER (VSL_SS_MATRIX_STORAGE_FULL     = Z"00000000")
      PARAMETER (VSL_SS_MATRIX_STORAGE_L_PACKED = Z"00000001")
      PARAMETER (VSL_SS_MATRIX_STORAGE_U_PACKED = Z"00000002")


!++
!  SUMMARY STATISTICS LIBRARY METHODS
!--

! SS routines provide computation of basic statistical estimates
! (central/raw moments up to 4th order, variance-covariance,
!  minimum, maximum, skewness/kurtosis) using the following methods
!  - FAST  - estimates are computed for price of one or two passes over
!            observations using highly optimized Intel(R) MKL routines
!  - 1PASS - estimate is computed for price of one pass of the observations
!  - FAST_USER_MEAN - estimates are computed for price of one or two passes
!            over observations given user defined mean for central moments,
!            covariance and correlation
!  - CP_TO_COVCOR - convert cross-product matrix to variance-covariance/
!            correlation matrix
!  - SUM_TO_MOM - convert raw/central sums to raw/central moments

      INTEGER VSL_SS_METHOD_FAST
      INTEGER VSL_SS_METHOD_1PASS
      INTEGER VSL_SS_METHOD_FAST_USER_MEAN
      INTEGER VSL_SS_METHOD_CP_TO_COVCOR
      INTEGER VSL_SS_METHOD_SUM_TO_MOM

      PARAMETER (VSL_SS_METHOD_FAST  = Z"00000001")
      PARAMETER (VSL_SS_METHOD_1PASS = Z"00000002")
      PARAMETER (VSL_SS_METHOD_FAST_USER_MEAN = Z"00000100")
      PARAMETER (VSL_SS_METHOD_CP_TO_COVCOR   = Z"00000200")
      PARAMETER (VSL_SS_METHOD_SUM_TO_MOM     = Z"00000400")


! SS provides routine for parametrization of correlation matrix using
! SPECTRAL DECOMPOSITION (SD) method
      INTEGER VSL_SS_METHOD_SD
      PARAMETER (VSL_SS_METHOD_SD    = Z"00000004")

! SS routine for robust estimation of variance-covariance matrix
! and mean supports Rocke algorithm, TBS-estimator
      INTEGER VSL_SS_METHOD_TBS
      PARAMETER (VSL_SS_METHOD_TBS  = Z"00000008")

!  SS routine for estimation of missing values
!  supports Multiple Imputation (MI) method
      INTEGER VSL_SS_METHOD_MI
      PARAMETER (VSL_SS_METHOD_MI   = Z"00000010")

!  SS routine for sorting data, RADIX method
      INTEGER VSL_SS_METHOD_RADIX
      PARAMETER (VSL_SS_METHOD_RADIX = Z"00100000")

! SS provides routine for detection of outliers, BACON method
      INTEGER VSL_SS_METHOD_BACON
      PARAMETER (VSL_SS_METHOD_BACON= Z"00000020")


! SS supports routine for estimation of quantiles for streaming data
! using the following methods:
! - ZW      - intermediate estimates of quantiles during processing
!             the next block are computed
! - ZW      - intermediate estimates of quantiles during processing
!             the next block are not computed
      INTEGER VSL_SS_METHOD_SQUANTS_ZW
      INTEGER VSL_SS_METHOD_SQUANTS_ZW_FAST
      PARAMETER (VSL_SS_METHOD_SQUANTS_ZW      = Z"00000040")
      PARAMETER (VSL_SS_METHOD_SQUANTS_ZW_FAST = Z"00000080")


! Input of BACON algorithm is set of 3 parameters:
! - Initialization method of the algorithm
! - Parameter alfa such that 1-alfa is percentile of Chi2 distribution
! - Stopping criterion

! Number of BACON algorithm parameters
      INTEGER VSL_SS_BACON_PARAMS_N
      PARAMETER (VSL_SS_BACON_PARAMS_N       = 3)


! SS implementation of BACON algorithm supports two initialization methods:
! - Mahalanobis distance based method
! - Median based method
      INTEGER(KIND=4) VSL_SS_METHOD_BACON_MAHALANOBIS_INIT
      INTEGER(KIND=4) VSL_SS_METHOD_BACON_MEDIAN_INIT
      PARAMETER (VSL_SS_METHOD_BACON_MAHALANOBIS_INIT = Z"00000001")
      PARAMETER (VSL_SS_METHOD_BACON_MEDIAN_INIT      = Z"00000002")


! Input of TBS algorithm is set of 4 parameters:
! - Breakdown point
! - Asymptotic rejection probability
! - Stopping criterion
! - Maximum number of iterations

! Number of TBS algorithm parameters
      INTEGER VSL_SS_TBS_PARAMS_N
      PARAMETER (VSL_SS_TBS_PARAMS_N         = 4)


! Input of MI algorithm is set of 5 parameters:
! - Maximal number of iterations for EM algorithm
! - Maximal number of iterations for DA algorithm
! - Stopping criterion
! - Number of sets to impute
! - Total number of missing values in dataset

! Number of MI algorithm parameters
      INTEGER VSL_SS_MI_PARAMS_SIZE
      PARAMETER (VSL_SS_MI_PARAMS_SIZE       = 5)


! SS MI algorithm expects that missing values are
! marked with NANs
      REAL(KIND=8) VSL_SS_DNAN
      REAL(KIND=4) VSL_SS_SNAN
       PARAMETER (VSL_SS_DNAN = Z"FFF8000000000000")
       PARAMETER (VSL_SS_SNAN = Z"FFC00000")


! Input of ZW algorithm is 1 parameter:
! - accuracy of quantile estimation

! Number of ZW algorithm parameters
      INTEGER VSL_SS_SQUANTS_ZW_PARAMS_N
      PARAMETER (VSL_SS_SQUANTS_ZW_PARAMS_N = 1)



!++
!  MACROS USED SS EDIT AND COMPUTE ROUTINES
!--

! SS EditTask routine is way to edit input and output parameters of the task,
! e.g., pointers to arrays which hold observations, weights of observations,
! arrays of mean estimates or covariance estimates.
! Macros below define parameters available for modification
      INTEGER VSL_SS_ED_DIMEN
      INTEGER VSL_SS_ED_OBSERV_N
      INTEGER VSL_SS_ED_OBSERV
      INTEGER VSL_SS_ED_OBSERV_STORAGE
      INTEGER VSL_SS_ED_INDC
      INTEGER VSL_SS_ED_WEIGHTS
      INTEGER VSL_SS_ED_MEAN
      INTEGER VSL_SS_ED_2R_MOM
      INTEGER VSL_SS_ED_3R_MOM
      INTEGER VSL_SS_ED_4R_MOM
      INTEGER VSL_SS_ED_2C_MOM
      INTEGER VSL_SS_ED_3C_MOM
      INTEGER VSL_SS_ED_4C_MOM
      INTEGER VSL_SS_ED_KURTOSIS
      INTEGER VSL_SS_ED_SKEWNESS
      INTEGER VSL_SS_ED_MIN
      INTEGER VSL_SS_ED_MAX
      INTEGER VSL_SS_ED_SORTED_OBSERV
      INTEGER VSL_SS_ED_SORTED_OBSERV_STORAGE
      INTEGER VSL_SS_ED_VARIATION
      INTEGER VSL_SS_ED_COV
      INTEGER VSL_SS_ED_COV_STORAGE
      INTEGER VSL_SS_ED_COR
      INTEGER VSL_SS_ED_COR_STORAGE
      INTEGER VSL_SS_ED_ACCUM_WEIGHT
      INTEGER VSL_SS_ED_QUANT_ORDER_N
      INTEGER VSL_SS_ED_QUANT_ORDER
      INTEGER VSL_SS_ED_QUANT_QUANTILES
      INTEGER VSL_SS_ED_ORDER_STATS
      INTEGER VSL_SS_ED_GROUP_INDC
      INTEGER VSL_SS_ED_POOLED_COV_STORAGE
      INTEGER VSL_SS_ED_POOLED_MEAN
      INTEGER VSL_SS_ED_POOLED_COV
      INTEGER VSL_SS_ED_GROUP_COV_INDC
      INTEGER VSL_SS_ED_REQ_GROUP_INDC
      INTEGER VSL_SS_ED_GROUP_MEAN
      INTEGER VSL_SS_ED_GROUP_COV_STORAGE
      INTEGER VSL_SS_ED_GROUP_COV
      INTEGER VSL_SS_ED_ROBUST_COV_STORAGE
      INTEGER VSL_SS_ED_ROBUST_COV_PARAMS_N
      INTEGER VSL_SS_ED_ROBUST_COV_PARAMS
      INTEGER VSL_SS_ED_ROBUST_MEAN
      INTEGER VSL_SS_ED_ROBUST_COV
      INTEGER VSL_SS_ED_OUTLIERS_PARAMS_N
      INTEGER VSL_SS_ED_OUTLIERS_PARAMS
      INTEGER VSL_SS_ED_OUTLIERS_WEIGHT
      INTEGER VSL_SS_ED_ORDER_STATS_STORAGE
      INTEGER VSL_SS_ED_PARTIAL_COV_IDX
      INTEGER VSL_SS_ED_PARTIAL_COV
      INTEGER VSL_SS_ED_PARTIAL_COV_STORAGE
      INTEGER VSL_SS_ED_PARTIAL_COR
      INTEGER VSL_SS_ED_PARTIAL_COR_STORAGE
      INTEGER VSL_SS_ED_MI_PARAMS_N
      INTEGER VSL_SS_ED_MI_PARAMS
      INTEGER VSL_SS_ED_MI_INIT_ESTIMATES_N
      INTEGER VSL_SS_ED_MI_INIT_ESTIMATES
      INTEGER VSL_SS_ED_MI_SIMUL_VALS_N
      INTEGER VSL_SS_ED_MI_SIMUL_VALS
      INTEGER VSL_SS_ED_MI_ESTIMATES_N
      INTEGER VSL_SS_ED_MI_ESTIMATES
      INTEGER VSL_SS_ED_MI_PRIOR_N
      INTEGER VSL_SS_ED_MI_PRIOR
      INTEGER VSL_SS_ED_PARAMTR_COR
      INTEGER VSL_SS_ED_PARAMTR_COR_STORAGE
      INTEGER VSL_SS_ED_STREAM_QUANT_PARAMS_N
      INTEGER VSL_SS_ED_STREAM_QUANT_PARAMS
      INTEGER VSL_SS_ED_STREAM_QUANT_ORDER_N
      INTEGER VSL_SS_ED_STREAM_QUANT_ORDER
      INTEGER VSL_SS_ED_STREAM_QUANT_QUANTILES
      INTEGER VSL_SS_ED_SUM
      INTEGER VSL_SS_ED_2R_SUM
      INTEGER VSL_SS_ED_3R_SUM
      INTEGER VSL_SS_ED_4R_SUM
      INTEGER VSL_SS_ED_2C_SUM
      INTEGER VSL_SS_ED_3C_SUM
      INTEGER VSL_SS_ED_4C_SUM
      INTEGER VSL_SS_ED_CP
      INTEGER VSL_SS_ED_CP_STORAGE
      INTEGER VSL_SS_ED_MDAD
      INTEGER VSL_SS_ED_MNAD

      PARAMETER (VSL_SS_ED_DIMEN                  =  1)
      PARAMETER (VSL_SS_ED_OBSERV_N               =  2)
      PARAMETER (VSL_SS_ED_OBSERV                 =  3)
      PARAMETER (VSL_SS_ED_OBSERV_STORAGE         =  4)
      PARAMETER (VSL_SS_ED_INDC                   =  5)
      PARAMETER (VSL_SS_ED_WEIGHTS                =  6)
      PARAMETER (VSL_SS_ED_MEAN                   =  7)
      PARAMETER (VSL_SS_ED_2R_MOM                 =  8)
      PARAMETER (VSL_SS_ED_3R_MOM                 =  9)
      PARAMETER (VSL_SS_ED_4R_MOM                 = 10)
      PARAMETER (VSL_SS_ED_2C_MOM                 = 11)
      PARAMETER (VSL_SS_ED_3C_MOM                 = 12)
      PARAMETER (VSL_SS_ED_4C_MOM                 = 13)
      PARAMETER (VSL_SS_ED_KURTOSIS               = 14)
      PARAMETER (VSL_SS_ED_SKEWNESS               = 15)
      PARAMETER (VSL_SS_ED_MIN                    = 16)
      PARAMETER (VSL_SS_ED_MAX                    = 17)
      PARAMETER (VSL_SS_ED_VARIATION              = 18)
      PARAMETER (VSL_SS_ED_COV                    = 19)
      PARAMETER (VSL_SS_ED_COV_STORAGE            = 20)
      PARAMETER (VSL_SS_ED_COR                    = 21)
      PARAMETER (VSL_SS_ED_COR_STORAGE            = 22)
      PARAMETER (VSL_SS_ED_ACCUM_WEIGHT           = 23)
      PARAMETER (VSL_SS_ED_QUANT_ORDER_N          = 24)
      PARAMETER (VSL_SS_ED_QUANT_ORDER            = 25)
      PARAMETER (VSL_SS_ED_QUANT_QUANTILES        = 26)
      PARAMETER (VSL_SS_ED_ORDER_STATS            = 27)
      PARAMETER (VSL_SS_ED_GROUP_INDC             = 28)
      PARAMETER (VSL_SS_ED_POOLED_COV_STORAGE     = 29)
      PARAMETER (VSL_SS_ED_POOLED_MEAN            = 30)
      PARAMETER (VSL_SS_ED_POOLED_COV             = 31)
      PARAMETER (VSL_SS_ED_GROUP_COV_INDC         = 32)
      PARAMETER (VSL_SS_ED_REQ_GROUP_INDC         = 32)
      PARAMETER (VSL_SS_ED_GROUP_MEAN             = 33)
      PARAMETER (VSL_SS_ED_GROUP_COV_STORAGE      = 34)
      PARAMETER (VSL_SS_ED_GROUP_COV              = 35)
      PARAMETER (VSL_SS_ED_ROBUST_COV_STORAGE     = 36)
      PARAMETER (VSL_SS_ED_ROBUST_COV_PARAMS_N    = 37)
      PARAMETER (VSL_SS_ED_ROBUST_COV_PARAMS      = 38)
      PARAMETER (VSL_SS_ED_ROBUST_MEAN            = 39)
      PARAMETER (VSL_SS_ED_ROBUST_COV             = 40)
      PARAMETER (VSL_SS_ED_OUTLIERS_PARAMS_N      = 41)
      PARAMETER (VSL_SS_ED_OUTLIERS_PARAMS        = 42)
      PARAMETER (VSL_SS_ED_OUTLIERS_WEIGHT        = 43)
      PARAMETER (VSL_SS_ED_ORDER_STATS_STORAGE    = 44)
      PARAMETER (VSL_SS_ED_PARTIAL_COV_IDX        = 45)
      PARAMETER (VSL_SS_ED_PARTIAL_COV            = 46)
      PARAMETER (VSL_SS_ED_PARTIAL_COV_STORAGE    = 47)
      PARAMETER (VSL_SS_ED_PARTIAL_COR            = 48)
      PARAMETER (VSL_SS_ED_PARTIAL_COR_STORAGE    = 49)
      PARAMETER (VSL_SS_ED_MI_PARAMS_N            = 50)
      PARAMETER (VSL_SS_ED_MI_PARAMS              = 51)
      PARAMETER (VSL_SS_ED_MI_INIT_ESTIMATES_N    = 52)
      PARAMETER (VSL_SS_ED_MI_INIT_ESTIMATES      = 53)
      PARAMETER (VSL_SS_ED_MI_SIMUL_VALS_N        = 54)
      PARAMETER (VSL_SS_ED_MI_SIMUL_VALS          = 55)
      PARAMETER (VSL_SS_ED_MI_ESTIMATES_N         = 56)
      PARAMETER (VSL_SS_ED_MI_ESTIMATES           = 57)
      PARAMETER (VSL_SS_ED_MI_PRIOR_N             = 58)
      PARAMETER (VSL_SS_ED_MI_PRIOR               = 59)
      PARAMETER (VSL_SS_ED_PARAMTR_COR            = 60)
      PARAMETER (VSL_SS_ED_PARAMTR_COR_STORAGE    = 61)
      PARAMETER (VSL_SS_ED_STREAM_QUANT_PARAMS_N  = 62)
      PARAMETER (VSL_SS_ED_STREAM_QUANT_PARAMS    = 63)
      PARAMETER (VSL_SS_ED_STREAM_QUANT_ORDER_N   = 64)
      PARAMETER (VSL_SS_ED_STREAM_QUANT_ORDER     = 65)
      PARAMETER (VSL_SS_ED_STREAM_QUANT_QUANTILES = 66)
      PARAMETER (VSL_SS_ED_SUM                    = 67)
      PARAMETER (VSL_SS_ED_2R_SUM                 = 68)
      PARAMETER (VSL_SS_ED_3R_SUM                 = 69)
      PARAMETER (VSL_SS_ED_4R_SUM                 = 70)
      PARAMETER (VSL_SS_ED_2C_SUM                 = 71)
      PARAMETER (VSL_SS_ED_3C_SUM                 = 72)
      PARAMETER (VSL_SS_ED_4C_SUM                 = 73)
      PARAMETER (VSL_SS_ED_CP                     = 74)
      PARAMETER (VSL_SS_ED_CP_STORAGE             = 75)
      PARAMETER (VSL_SS_ED_MDAD                   = 76)
      PARAMETER (VSL_SS_ED_MNAD                   = 77)
      PARAMETER (VSL_SS_ED_SORTED_OBSERV          = 78)
      PARAMETER (VSL_SS_ED_SORTED_OBSERV_STORAGE  = 79)



! SS Compute routine calculates estimates supported by the library
! Macros below define estimates to compute
       INTEGER(KIND=8) VSL_SS_MEAN
       INTEGER(KIND=8) VSL_SS_2R_MOM
       INTEGER(KIND=8) VSL_SS_3R_MOM
       INTEGER(KIND=8) VSL_SS_4R_MOM
       INTEGER(KIND=8) VSL_SS_2C_MOM
       INTEGER(KIND=8) VSL_SS_3C_MOM
       INTEGER(KIND=8) VSL_SS_4C_MOM
       INTEGER(KIND=8) VSL_SS_SUM
       INTEGER(KIND=8) VSL_SS_2R_SUM
       INTEGER(KIND=8) VSL_SS_3R_SUM
       INTEGER(KIND=8) VSL_SS_4R_SUM
       INTEGER(KIND=8) VSL_SS_2C_SUM
       INTEGER(KIND=8) VSL_SS_3C_SUM
       INTEGER(KIND=8) VSL_SS_4C_SUM
       INTEGER(KIND=8) VSL_SS_KURTOSIS
       INTEGER(KIND=8) VSL_SS_SKEWNESS
       INTEGER(KIND=8) VSL_SS_VARIATION
       INTEGER(KIND=8) VSL_SS_MIN
       INTEGER(KIND=8) VSL_SS_MAX
       INTEGER(KIND=8) VSL_SS_COV
       INTEGER(KIND=8) VSL_SS_COR
       INTEGER(KIND=8) VSL_SS_CP
       INTEGER(KIND=8) VSL_SS_POOLED_COV
       INTEGER(KIND=8) VSL_SS_GROUP_COV
       INTEGER(KIND=8) VSL_SS_POOLED_MEAN
       INTEGER(KIND=8) VSL_SS_GROUP_MEAN
       INTEGER(KIND=8) VSL_SS_QUANTS
       INTEGER(KIND=8) VSL_SS_SORTED_OBSERV
       INTEGER(KIND=8) VSL_SS_ORDER_STATS
       INTEGER(KIND=8) VSL_SS_ROBUST_COV
       INTEGER(KIND=8) VSL_SS_OUTLIERS
       INTEGER(KIND=8) VSL_SS_PARTIAL_COV
       INTEGER(KIND=8) VSL_SS_PARTIAL_COR
       INTEGER(KIND=8) VSL_SS_MISSING_VALS
       INTEGER(KIND=8) VSL_SS_PARAMTR_COR
       INTEGER(KIND=8) VSL_SS_STREAM_QUANTS
       INTEGER(KIND=8) VSL_SS_MDAD
       INTEGER(KIND=8) VSL_SS_MNAD
       PARAMETER ( VSL_SS_MEAN          = Z"0000000000000001" )
       PARAMETER ( VSL_SS_2R_MOM        = Z"0000000000000002" )
       PARAMETER ( VSL_SS_3R_MOM        = Z"0000000000000004" )
       PARAMETER ( VSL_SS_4R_MOM        = Z"0000000000000008" )
       PARAMETER ( VSL_SS_2C_MOM        = Z"0000000000000010" )
       PARAMETER ( VSL_SS_3C_MOM        = Z"0000000000000020" )
       PARAMETER ( VSL_SS_4C_MOM        = Z"0000000000000040" )
       PARAMETER ( VSL_SS_SUM           = Z"0000000002000000" )
       PARAMETER ( VSL_SS_2R_SUM        = Z"0000000004000000" )
       PARAMETER ( VSL_SS_3R_SUM        = Z"0000000008000000" )
       PARAMETER ( VSL_SS_4R_SUM        = Z"0000000010000000" )
       PARAMETER ( VSL_SS_2C_SUM        = Z"0000000020000000" )
       PARAMETER ( VSL_SS_3C_SUM        = Z"0000000040000000" )
       PARAMETER ( VSL_SS_4C_SUM        = Z"0000000080000000" )
       PARAMETER ( VSL_SS_KURTOSIS      = Z"0000000000000080" )
       PARAMETER ( VSL_SS_SKEWNESS      = Z"0000000000000100" )
       PARAMETER ( VSL_SS_VARIATION     = Z"0000000000000200" )
       PARAMETER ( VSL_SS_MIN           = Z"0000000000000400" )
       PARAMETER ( VSL_SS_MAX           = Z"0000000000000800" )
       PARAMETER ( VSL_SS_COV           = Z"0000000000001000" )
       PARAMETER ( VSL_SS_COR           = Z"0000000000002000" )
       PARAMETER ( VSL_SS_CP            = Z"0000000100000000" )
       PARAMETER ( VSL_SS_POOLED_COV    = Z"0000000000004000" )
       PARAMETER ( VSL_SS_GROUP_COV     = Z"0000000000008000" )
       PARAMETER ( VSL_SS_POOLED_MEAN   = Z"0000000800000000" )
       PARAMETER ( VSL_SS_GROUP_MEAN    = Z"0000001000000000" )
       PARAMETER ( VSL_SS_QUANTS        = Z"0000000000010000" )
       PARAMETER ( VSL_SS_SORTED_OBSERV = Z"0000008000000000" )
       PARAMETER ( VSL_SS_ORDER_STATS   = Z"0000000000020000" )
       PARAMETER ( VSL_SS_ROBUST_COV    = Z"0000000000040000" )
       PARAMETER ( VSL_SS_OUTLIERS      = Z"0000000000080000" )
       PARAMETER ( VSL_SS_PARTIAL_COV   = Z"0000000000100000" )
       PARAMETER ( VSL_SS_PARTIAL_COR   = Z"0000000000200000" )
       PARAMETER ( VSL_SS_MISSING_VALS  = Z"0000000000400000" )
       PARAMETER ( VSL_SS_PARAMTR_COR   = Z"0000000000800000" )
       PARAMETER ( VSL_SS_STREAM_QUANTS = Z"0000000001000000" )
       PARAMETER ( VSL_SS_MDAD          = Z"0000000200000000" )
       PARAMETER ( VSL_SS_MNAD          = Z"0000000400000000" )


!++
!  TYPEDEFS
!--

!  VSL STREAM STATE POINTER
!  This structure is to store VSL stream state address allocated by
!  VSLNEWSTREAM subroutine.
      TYPE VSL_STREAM_STATE
          INTEGER(KIND=4) descriptor1
          INTEGER(KIND=4) descriptor2
      END TYPE VSL_STREAM_STATE

      TYPE VSL_CONV_TASK
          INTEGER(KIND=4) descriptor1
          INTEGER(KIND=4) descriptor2
      END TYPE VSL_CONV_TASK

      TYPE VSL_CORR_TASK
          INTEGER(KIND=4) descriptor1
          INTEGER(KIND=4) descriptor2
      END TYPE VSL_CORR_TASK

      TYPE VSL_SS_TASK
          INTEGER(KIND=4) descriptor1
          INTEGER(KIND=4) descriptor2
      END TYPE VSL_SS_TASK

!  BASIC RANDOM NUMBER GENERATOR PROPERTIES STRUCTURE
!  The structure describes the properties of given basic generator, e.g. size
!  of the stream state structure, pointers to function implementations, etc.
!
!  BRNG properties structure fields:
!  StreamStateSize - size of the stream state structure (in bytes)
!  WordSize        - size of base word (in bytes). Typically this is 4 bytes.
!  NSeeds          - number of words necessary to describe generator's state
!  NBits           - number of bits actually used in base word. For example,
!                    only 31 least significant bits are actually used in
!                    basic random number generator MCG31m1 with 4-byte base
!                    word. NBits field is useful while interpreting random
!                    words as a sequence of random bits.
!  IncludesZero    - FALSE if 0 cannot be generated in integer-valued
!                    implementation; TRUE if 0 can be potentially generated in
!                    integer-valued implementation.
!  reserved        - a reserved field for internal needs
      TYPE VSL_BRNG_PROPERTIES
          INTEGER(KIND=4) streamstatesize
          INTEGER(KIND=4) nseeds
          INTEGER(KIND=4) includeszero
          INTEGER(KIND=4) wordsize
          INTEGER(KIND=4) nbits
          INTEGER(KIND=4) reserved(8)
      END TYPE VSL_BRNG_PROPERTIES

      END MODULE MKL_VSL_TYPE

      MODULE MKL_VSL

      USE MKL_VSL_TYPE


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!==============================================================================
!------------------------------------------------------------------------------

      INTERFACE
        INTEGER FUNCTION vsldconvnewtask( task, mode, dims, xshape,     &
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims),zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvnewtask( task, mode, dims, xshape,     &
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims),zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzconvnewtask( task, mode, dims, xshape,     &
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims),zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcconvnewtask( task, mode, dims, xshape,     &
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims),zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrnewtask( task, mode, dims, xshape,     &
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims),zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrnewtask( task, mode, dims, xshape,     &
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims),zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzcorrnewtask( task, mode, dims, xshape,     &
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims),zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslccorrnewtask( task, mode, dims, xshape,     &
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims),zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvnewtask1d( task, mode, xshape, yshape, &
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvnewtask1d( task, mode, xshape, yshape, &
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzconvnewtask1d( task, mode, xshape, yshape, &
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcconvnewtask1d( task, mode, xshape, yshape, &
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrnewtask1d( task, mode, xshape, yshape, &
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrnewtask1d( task, mode, xshape, yshape, &
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzcorrnewtask1d( task, mode, xshape, yshape, &
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslccorrnewtask1d( task, mode, xshape, yshape, &
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvnewtaskx( task, mode, dims, xshape,    &
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          REAL(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvnewtaskx( task, mode, dims, xshape,    &
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          REAL(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzconvnewtaskx( task, mode, dims, xshape,    &
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          COMPLEX(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcconvnewtaskx( task, mode, dims, xshape,    &
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          COMPLEX(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrnewtaskx( task, mode, dims, xshape,    &
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          REAL(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrnewtaskx( task, mode, dims, xshape,    &
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          REAL(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzcorrnewtaskx( task, mode, dims, xshape,    &
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          COMPLEX(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslccorrnewtaskx( task, mode, dims, xshape,    &
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          COMPLEX(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvnewtaskx1d( task, mode, xshape,        &
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape,xstride
          REAL(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvnewtaskx1d( task, mode, xshape,        &
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape,xstride
          REAL(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzconvnewtaskx1d( task, mode, xshape,        &
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape,xstride
          COMPLEX(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcconvnewtaskx1d( task, mode, xshape,        &
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape,xstride
          COMPLEX(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrnewtaskx1d( task, mode, xshape,        &
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape,xstride
          REAL(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrnewtaskx1d( task, mode, xshape,        &
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape,xstride
          REAL(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzcorrnewtaskx1d( task, mode, xshape,        &
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape,xstride
          COMPLEX(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslccorrnewtaskx1d( task, mode, xshape,        &
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode
          INTEGER             :: xshape,yshape,zshape,xstride
          COMPLEX(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslconvdeletetask( task )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcorrdeletetask( task )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslconvcopytask( desttask, srctask )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: desttask
          TYPE(VSL_CONV_TASK) :: srctask
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcorrcopytask( desttask, srctask )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: desttask
          TYPE(VSL_CORR_TASK) :: srctask
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslConvSetMode( task, newmode )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: newmode
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslCorrSetMode( task, newmode )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: newmode
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslConvSetInternalPrecision( task, precision )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: precision
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslCorrSetInternalPrecision( task, precision )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: precision
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslConvSetStart( task, start )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER,DIMENSION(*):: start
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslCorrSetStart( task, start )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER,DIMENSION(*):: start
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslConvSetDecimation( task, decimation )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER,DIMENSION(*):: decimation
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslCorrSetDecimation( task, decimation )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER,DIMENSION(*):: decimation
        END FUNCTION
      END INTERFACE


      INTERFACE
        INTEGER FUNCTION vsldconvexec( task, x, xstride, y, ystride, z, &
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          REAL(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvexec( task, x, xstride, y, ystride, z, &
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          REAL(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzconvexec( task, x, xstride, y, ystride, z, &
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          COMPLEX(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcconvexec( task, x, xstride, y, ystride, z, &
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          COMPLEX(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrexec( task, x, xstride, y, ystride, z, &
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          REAL(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrexec( task, x, xstride, y, ystride, z, &
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          REAL(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzcorrexec( task, x, xstride, y, ystride, z, &
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          COMPLEX(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslccorrexec( task, x, xstride, y, ystride, z, &
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          COMPLEX(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvexec1d( task, x, xstride, y, ystride,  &
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          REAL(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvexec1d( task, x, xstride, y, ystride,  &
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          REAL(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzconvexec1d( task, x, xstride, y, ystride,  &
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          COMPLEX(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcconvexec1d( task, x, xstride, y, ystride,  &
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          COMPLEX(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrexec1d( task, x, xstride, y, ystride,  &
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          REAL(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrexec1d( task, x, xstride, y, ystride,  &
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          REAL(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzcorrexec1d( task, x, xstride, y, ystride,  &
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          COMPLEX(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslccorrexec1d( task, x, xstride, y, ystride,  &
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          COMPLEX(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          REAL(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          REAL(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzconvexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          COMPLEX(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcconvexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          COMPLEX(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          REAL(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          REAL(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzcorrexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          COMPLEX(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslccorrexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          COMPLEX(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: ystride,zstride
          REAL(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: ystride,zstride
          REAL(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzconvexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: ystride,zstride
          COMPLEX(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcconvexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: ystride,zstride
          COMPLEX(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: ystride,zstride
          REAL(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: ystride,zstride
          REAL(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslzcorrexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: ystride,zstride
          COMPLEX(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslccorrexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: ystride,zstride
          COMPLEX(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

!++
!  VSL CONTINUOUS DISTRIBUTION GENERATOR FUNCTION INTERFACES.
!--

!  Uniform distribution
      INTERFACE
        INTEGER FUNCTION vsrnguniform( method, stream, n, r, a, b )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: b
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnguniform( method, stream, n, r, a, b )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: b
        END FUNCTION
      END INTERFACE

!  Gaussian distribution
      INTERFACE
        INTEGER FUNCTION vsrnggaussian( method, stream, n, r, a, sigma)
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: sigma
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnggaussian( method, stream, n, r, a, sigma)
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: sigma
        END FUNCTION
      END INTERFACE

!  GaussianMV distribution
      INTERFACE
        INTEGER FUNCTION vsrnggaussianmv( method, stream, n, r, dimen,  &
     &                              mstorage, a, t )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          INTEGER,INTENT(IN)       :: dimen
          REAL(KIND=4),INTENT(OUT) :: r(dimen,n)
          INTEGER,INTENT(IN)       :: mstorage
          REAL(KIND=4),INTENT(IN)  :: a(dimen)
          REAL(KIND=4),INTENT(IN)  :: t(dimen,dimen)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnggaussianmv( method, stream, n, r, dimen,  &
     &                              mstorage, a, t )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          INTEGER,INTENT(IN)       :: dimen
          INTEGER,INTENT(IN)       :: mstorage
          REAL(KIND=8),INTENT(IN)  :: a(dimen)
          REAL(KIND=8),INTENT(IN)  :: t(dimen,dimen)
        END FUNCTION
      END INTERFACE

!  Exponential distribution
      INTERFACE
        INTEGER FUNCTION vsrngexponential( method, stream, n, r, a,     &
     &                                     beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngexponential( method, stream, n, r, a,     &
     &                                     beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!  Laplace distribution
      INTERFACE
        INTEGER FUNCTION vsrnglaplace( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnglaplace( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!  Weibull distribution
      INTERFACE
        INTEGER FUNCTION vsrngweibull( method, stream, n, r, alpha, a,  &
     &                                 beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: alpha
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngweibull( method, stream, n, r, alpha, a,  &
     &                                 beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: alpha
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!  Cauchy distribution
      INTERFACE
        INTEGER FUNCTION vsrngcauchy( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngcauchy( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!  Rayleigh distribution
      INTERFACE
        INTEGER FUNCTION vsrngrayleigh( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngrayleigh( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!  Lognormal distribution
      INTERFACE
        INTEGER FUNCTION vsrnglognormal( method, stream, n, r, a,sigma, &
     &                             b, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: sigma
          REAL(KIND=4),INTENT(IN)  :: b
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnglognormal( method, stream, n, r, a,sigma, &
     &                             b, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: sigma
          REAL(KIND=8),INTENT(IN)  :: b
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!  Gumbel distribution
      INTERFACE
        INTEGER FUNCTION vsrnggumbel( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnggumbel( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!  Gamma distribution
      INTERFACE
        INTEGER FUNCTION vsrnggamma( method, stream, n, r, alpha, a,    &
     &                               beta )
          USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: alpha
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnggamma( method, stream, n, r, alpha, a,    &
     &                               beta )
          USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: alpha
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!  Beta distribution
      INTERFACE
        INTEGER FUNCTION vsrngbeta( method, stream, n, r, p, q, a,      &
     &                              beta )
          USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: p
          REAL(KIND=4),INTENT(IN)  :: q
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngbeta( method, stream, n, r, p, q, a,      &
     &                              beta )
          USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: p
          REAL(KIND=8),INTENT(IN)  :: q
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!++
!  VSL DISCRETE DISTRIBUTION GENERATOR FUNCTION INTERFACES.
!--

!  Uniform distribution
      INTERFACE
        INTEGER FUNCTION virnguniform( method, stream, n, r, a, b )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=4),INTENT(IN)  :: a
          INTEGER(KIND=4),INTENT(IN)  :: b
        END FUNCTION
      END INTERFACE

!  UniformBits distribution
      INTERFACE
        INTEGER FUNCTION virnguniformbits( method, stream, n, r )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
        END FUNCTION
      END INTERFACE

!  UniformBits32 distribution
      INTERFACE
        INTEGER FUNCTION virnguniformbits32( method, stream, n, r )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
        END FUNCTION
      END INTERFACE

!  UniformBits64 distribution
      INTERFACE
        INTEGER FUNCTION virnguniformbits64( method, stream, n, r )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=8),INTENT(OUT) :: r(n)
        END FUNCTION
      END INTERFACE

!  Bernoulli distribution
      INTERFACE
        INTEGER FUNCTION virngbernoulli( method, stream, n, r, p )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)     :: p
        END FUNCTION
      END INTERFACE

!  Geometric distribution
      INTERFACE
        INTEGER FUNCTION virnggeometric( method, stream, n, r, p )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)     :: p
        END FUNCTION
      END INTERFACE

!  Binomial distribution
      INTERFACE
        INTEGER FUNCTION virngbinomial(method, stream, n, r, ntrial, p)
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=4),INTENT(IN)  :: ntrial
          REAL(KIND=8),INTENT(IN)     :: p
        END FUNCTION
      END INTERFACE

!  Hypergeometric distribution
      INTERFACE
        INTEGER FUNCTION virnghypergeometric( method, stream, n, r, l,  &
     &                                        s, m )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=4),INTENT(IN)  :: l
          INTEGER(KIND=4),INTENT(IN)  :: s
          INTEGER(KIND=4),INTENT(IN)  :: m
        END FUNCTION
      END INTERFACE

!  Poisson distribution
      INTERFACE
        INTEGER FUNCTION virngpoisson( method, stream, n, r, lambda )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)     :: lambda
        END FUNCTION
      END INTERFACE

!  PoissonV distribution
      INTERFACE
        INTEGER FUNCTION virngpoissonv( method, stream, n, r, lambda )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)     :: lambda(n)
        END FUNCTION
      END INTERFACE

!  Negbinomial distribution
      INTERFACE
        INTEGER FUNCTION virngnegbinomial( method, stream, n, r, a, p )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)          :: method
          TYPE(VSL_STREAM_STATE)      :: stream
          INTEGER,INTENT(IN)          :: n
          INTEGER(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)     :: a
          REAL(KIND=8),INTENT(IN)     :: p
        END FUNCTION
      END INTERFACE

!++
!  VSL SERVICE FUNCTION INTERFACES.
!--

! NewStream - stream creation/initialization
      INTERFACE
        INTEGER FUNCTION vslnewstream( stream, brng, seed )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: brng
          INTEGER,INTENT(IN)     :: seed
        END FUNCTION
      END INTERFACE

! NewStreamEx - advanced stream creation/initialization
      INTERFACE
        INTEGER FUNCTION vslnewstreamex( stream, brng, n, params )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE)     :: stream
          INTEGER,INTENT(IN)         :: brng
          INTEGER,INTENT(IN)         :: n
          INTEGER(KIND=4),INTENT(IN) :: params(n)
        END FUNCTION
      END INTERFACE

!    INEWABSTRACTSTREAM
      INTERFACE
        INTEGER FUNCTION vslinewabstractstream( stream, n, ibuf, ifunc)
          USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE),INTENT(OUT) :: stream
          INTEGER,INTENT(IN)                 :: n
          INTEGER(KIND=4),INTENT(IN)         :: ibuf(n)
          INTEGER(KIND=4),EXTERNAL           :: ifunc
        END FUNCTION
      END INTERFACE

!    DNEWABSTRACTSTREAM
      INTERFACE
        INTEGER FUNCTION vsldnewabstractstream( stream, n, dbuf, a, b,  &
     &                                          dfunc )
          USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE),INTENT(OUT) :: stream
          INTEGER,INTENT(IN)                 :: n
          REAL(KIND=8)   ,INTENT(IN)         :: dbuf(n)
          REAL(KIND=8)   ,INTENT(IN)         :: a
          REAL(KIND=8)   ,INTENT(IN)         :: b
          INTEGER(KIND=4),EXTERNAL           :: dfunc
        END FUNCTION
      END INTERFACE

!    SNEWABSTRACTSTREAM
      INTERFACE
        INTEGER FUNCTION vslsnewabstractstream( stream, n, sbuf, a, b,  &
     &                                          sfunc )
          USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE),INTENT(OUT) :: stream
          INTEGER,INTENT(IN)                 :: n
          REAL(KIND=4)   ,INTENT(IN)         :: sbuf(n)
          REAL(KIND=4)   ,INTENT(IN)         :: a
          REAL(KIND=4)   ,INTENT(IN)         :: b
          INTEGER(KIND=4),EXTERNAL           :: sfunc
        END FUNCTION
      END INTERFACE

! DeleteStream - delete stream
      INTERFACE
        INTEGER FUNCTION vsldeletestream( stream )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
        END FUNCTION
      END INTERFACE

! CopyStream - copy all stream information
      INTERFACE
        INTEGER FUNCTION vslcopystream( newstream, srcstream )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: newstream
          TYPE(VSL_STREAM_STATE) :: srcstream
        END FUNCTION
      END INTERFACE

! CopyStreamState - copy stream state only
      INTERFACE
        INTEGER FUNCTION vslcopystreamstate( deststream, srcstream )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: deststream
          TYPE(VSL_STREAM_STATE) :: srcstream
        END FUNCTION
      END INTERFACE

! LeapfrogStream - leapfrog method
      INTERFACE
        INTEGER FUNCTION vslleapfrogstream( stream, k, nstreams )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: k
          INTEGER,INTENT(IN)     :: nstreams
        END FUNCTION
      END INTERFACE

! SkipAheadStream - skip-ahead method
      INTERFACE
        INTEGER FUNCTION vslskipaheadstream( stream, nskip )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE)     :: stream
          INTEGER(KIND=8),INTENT(IN) :: nskip
        END FUNCTION
      END INTERFACE

! GetBrngProperties - get BRNG properties
      INTERFACE
        INTEGER FUNCTION vslgetbrngproperties( brng, properties )
            USE MKL_VSL_TYPE
          INTEGER(KIND=4),INTENT(IN)            :: brng
          TYPE(VSL_BRNG_PROPERTIES),INTENT(OUT) :: properties
        END FUNCTION
      END INTERFACE

! GetNumRegBrngs - get number of registered BRNGs
      INTERFACE
        INTEGER FUNCTION vslgetnumregbrngs( )
        END FUNCTION
      END INTERFACE

! GetStreamStateBrng - get BRNG associated with given stream
      INTERFACE
        INTEGER FUNCTION vslgetstreamstatebrng( stream )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
        END FUNCTION
      END INTERFACE

! RegisterBrng - register new BRNG
      INTERFACE
        INTEGER FUNCTION vslregisterbrng( properties )
            USE MKL_VSL_TYPE
          TYPE(VSL_BRNG_PROPERTIES) :: properties
        END FUNCTION
      END INTERFACE

! SaveStreamF - save stream to file
      INTERFACE
        INTEGER FUNCTION vslsavestreamf( stream, fname )
            USE MKL_VSL_TYPE
          CHARACTER(*)           :: fname
          TYPE(VSL_STREAM_STATE) :: stream
        END FUNCTION
      END INTERFACE

! LoadStreamF - save stream to file
      INTERFACE
        INTEGER FUNCTION vslloadstreamf( stream, fname )
            USE MKL_VSL_TYPE
          CHARACTER(*)           :: fname
          TYPE(VSL_STREAM_STATE) :: stream
        END FUNCTION
      END INTERFACE

! SaveStreamM - save stream to memory
      INTERFACE
        INTEGER FUNCTION vslsavestreamm( stream, memptr )
            USE MKL_VSL_TYPE
          INTEGER(KIND=1),DIMENSION(*),INTENT(OUT)::memptr
          TYPE(VSL_STREAM_STATE),INTENT(IN)       :: stream
        END FUNCTION
      END INTERFACE

! LoadStreamM - load stream from memory
      INTERFACE
        INTEGER FUNCTION vslloadstreamm( stream, memptr )
            USE MKL_VSL_TYPE
          INTEGER(KIND=1),DIMENSION(*),INTENT(IN)::memptr
          TYPE(VSL_STREAM_STATE),INTENT(OUT)      ::stream
        END FUNCTION
      END INTERFACE

! GetStreamSize - get size of random stream
      INTERFACE
        INTEGER FUNCTION vslgetstreamsize( stream )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE),INTENT(IN) :: stream
        END FUNCTION
      END INTERFACE

!++
!  SUMMARARY STATTISTICS LIBARY ROUTINES
!--

!  Task constructors
      INTERFACE
       INTEGER FUNCTION vsldssnewtask(task,p,n,x_storage,x,w,indices)
                 USE MKL_VSL_TYPE
             TYPE(VSL_SS_TASK)       :: task
             INTEGER,INTENT(IN)      :: p
             INTEGER,INTENT(IN)      :: n
             INTEGER,INTENT(IN)      :: x_storage
             REAL(KIND=8),INTENT(IN) :: x(n,p)
             REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL:: w
             INTEGER,DIMENSION(*),INTENT(IN),OPTIONAL:: indices
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslsssnewtask(task,p,n,x_storage,x,w,indices)
                 USE MKL_VSL_TYPE
             TYPE(VSL_SS_TASK)       :: task
             INTEGER,INTENT(IN)      :: p
             INTEGER,INTENT(IN)      :: n
             INTEGER,INTENT(IN)      :: x_storage
             REAL(KIND=4),INTENT(IN) :: x(n,p)
             REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL:: w
             INTEGER,DIMENSION(*),INTENT(IN),OPTIONAL:: indices
       END FUNCTION
      END INTERFACE


!  Task editors

!  Editor to modify a task parameter
      INTERFACE
       INTEGER FUNCTION vsldssedittask(task,parameter,par_addr)
                 USE MKL_VSL_TYPE
             TYPE(VSL_SS_TASK)                    :: task
             INTEGER,INTENT(IN)                   :: parameter
             REAL(KIND=8),DIMENSION(*),INTENT(IN) :: par_addr
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslsssedittask(task,parameter,par_addr)
                 USE MKL_VSL_TYPE
             TYPE(VSL_SS_TASK)                    :: task
             INTEGER,INTENT(IN)                   :: parameter
             REAL(KIND=4),DIMENSION(*),INTENT(IN) :: par_addr
        END FUNCTION
       END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslissedittask(task,parameter,par_addr)
                 USE MKL_VSL_TYPE
             TYPE(VSL_SS_TASK)               :: task
             INTEGER,INTENT(IN)              :: parameter
             INTEGER,INTENT(IN)              :: par_addr
       END FUNCTION
      END INTERFACE


!  Task specific editors

!  Editors to modify moments related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditmoments(task, mean, r2m, r3m, r4m,     &
     &                                                c2m, c3m, c4m)
               USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                             :: task
              REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: mean
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: r2m
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: r3m
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: r4m
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: c2m
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: c3m
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: c4m
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditmoments(task, mean, r2m, r3m, r4m,     &
     &                                                c2m, c3m, c4m)
               USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                             :: task
              REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: mean
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: r2m
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: r3m
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: r4m
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: c2m
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: c3m
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: c4m
       END FUNCTION
      END INTERFACE

!  Editors to modify sums related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditsums(task, sum, r2s, r3s, r4s,         &
     &                                                c2s, c3s, c4s)
               USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                             :: task
              REAL(KIND=8),DIMENSION(*),INTENT(IN)          :: sum
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: r2s
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: r3s
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: r4s
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: c2s
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: c3s
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: c4s
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditsums(task, sum, r2s, r3s, r4s,         &
     &                                                c2s, c3s, c4s)
               USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                             :: task
              REAL(KIND=4),DIMENSION(*),INTENT(IN)          :: sum
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: r2s
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: r3s
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: r4s
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: c2s
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: c3s
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: c4s
       END FUNCTION
      END INTERFACE

!  Editors to modify variance-covariance/correlation matrix
!  related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditcovcor(task, mean,cov, cov_storage,    &
     &                                              cor, cor_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                              :: task
              REAL(KIND=8),DIMENSION(*),INTENT(IN)           :: mean
              REAL(KIND=8),DIMENSION(*),INTENT(IN), OPTIONAL :: cov
              INTEGER,INTENT(IN), OPTIONAL                :: cov_storage
              REAL(KIND=8),DIMENSION(*),INTENT(IN), OPTIONAL :: cor
              INTEGER,INTENT(IN), OPTIONAL                :: cor_storage
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditcovcor(task, mean,cov, cov_storage,    &
     &                                              cor, cor_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                              :: task
              REAL(KIND=4),DIMENSION(*),INTENT(IN)           :: mean
              REAL(KIND=4),DIMENSION(*),INTENT(IN), OPTIONAL :: cov
              INTEGER,INTENT(IN), OPTIONAL                :: cov_storage
              REAL(KIND=4),DIMENSION(*),INTENT(IN), OPTIONAL :: cor
              INTEGER,INTENT(IN), OPTIONAL                :: cor_storage
       END FUNCTION
      END INTERFACE

!  Editors to modify cross-product matrix
!  related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditcp(task, mean, sum,                    &
     &                                              cp, cp_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                              :: task
              REAL(KIND=8),DIMENSION(*),INTENT(IN)           :: mean
              REAL(KIND=8),DIMENSION(*),INTENT(IN), OPTIONAL :: sum
              REAL(KIND=8),DIMENSION(*),INTENT(IN)           :: cp
              INTEGER,INTENT(IN)                          :: cp_storage
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditcp(task, mean, sum,                    &
     &                                              cp, cp_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                              :: task
              REAL(KIND=4),DIMENSION(*),INTENT(IN)           :: mean
              REAL(KIND=4),DIMENSION(*),INTENT(IN), OPTIONAL :: sum
              REAL(KIND=4),DIMENSION(*),INTENT(IN)           :: cp
              INTEGER,INTENT(IN)                          :: cp_storage
       END FUNCTION
      END INTERFACE

!  Editors to modify partial variance-covariance matrix
!  related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditpartialcovcor(task, p_idx_array,       &
     &                       cov, cov_storage, cor, cor_storage,         &
     &                       p_cov, p_cov_storage, p_cor, p_cor_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)       :: task
              INTEGER,DIMENSION(*),INTENT(IN)      :: p_idx_array
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: cov
              INTEGER,INTENT(IN),OPTIONAL              :: cov_storage
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: cor
              INTEGER,INTENT(IN),OPTIONAL              :: cor_storage
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: p_cov
              INTEGER,INTENT(IN),OPTIONAL      :: p_cov_storage
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: p_cor
              INTEGER,INTENT(IN),OPTIONAL      :: p_cor_storage
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditpartialcovcor(task, p_idx_array,       &
     &                       cov, cov_storage, cor, cor_storage,         &
     &                       p_cov, p_cov_storage, p_cor, p_cor_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)       :: task
              INTEGER,DIMENSION(*),INTENT(IN)      :: p_idx_array
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: cov
              INTEGER,INTENT(IN),OPTIONAL              :: cov_storage
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: cor
              INTEGER,INTENT(IN),OPTIONAL              :: cor_storage
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: p_cov
              INTEGER,INTENT(IN),OPTIONAL      :: p_cov_storage
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: p_cor
              INTEGER,INTENT(IN),OPTIONAL      :: p_cor_storage
       END FUNCTION
      END INTERFACE

!  Editors to modify quantiles related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditquantiles(task, quant_order_n,         &
     &                                 quant_order,quants,               &
     &                                 order_stats, order_stats_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                         :: task
              INTEGER,INTENT(IN),OPTIONAL               :: quant_order_n
              REAL(KIND=8),INTENT(IN),dimension(*),OPTIONAL::quant_order
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL::quants
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL ::           &
     &                                                     order_stats
              INTEGER,INTENT(IN),OPTIONAL  ::      order_stats_storage
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditquantiles(task, quant_order_n,         &
     &                                 quant_order,quants,               &
     &                                 order_stats, order_stats_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                         :: task
              INTEGER,INTENT(IN),OPTIONAL               :: quant_order_n
              REAL(KIND=4),INTENT(IN),dimension(*),OPTIONAL::quant_order
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL::quants
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL ::           &
     &                                                     order_stats
              INTEGER,INTENT(IN),OPTIONAL  ::      order_stats_storage
       END FUNCTION
      END INTERFACE


!  Editors to modify stream data quantiles related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditstreamquantiles(task,                  &
     &        quant_order_n, quant_order, quants, nparams, params)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                     :: task
              INTEGER,INTENT(IN)                    :: quant_order_n
              REAL(KIND=8),INTENT(IN),dimension(*)  :: quant_order
              REAL(KIND=8),DIMENSION(*),INTENT(IN)  :: quants
              INTEGER,INTENT(IN)                    :: nparams
              REAL(KIND=8),INTENT(IN), DIMENSION(*) :: params
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditstreamquantiles(task,                  &
     &        quant_order_n, quant_order, quants, nparams, params)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                     :: task
              INTEGER,INTENT(IN)                    :: quant_order_n
              REAL(KIND=4),INTENT(IN),dimension(*)  :: quant_order
              REAL(KIND=4),DIMENSION(*),INTENT(IN)  :: quants
              INTEGER,INTENT(IN)                    :: nparams
              REAL(KIND=4),INTENT(IN), DIMENSION(*) :: params
       END FUNCTION
      END INTERFACE


!  Editors to modify pooled/group variance-covariance matrix
!  related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditpooledcovariance(task, grp_indices,    &
     &           pld_mean, pld_cov, grp_cov_indices, grp_means, grp_cov)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                    :: task
              INTEGER,DIMENSION(*),INTENT(IN)      :: grp_indices
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: pld_mean
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: pld_cov
              INTEGER,DIMENSION(*),INTENT(IN),OPTIONAL ::grp_cov_indices
              REAL(KIND=8),DIMENSION(*),INTENT(IN),OPTIONAL :: grp_means
              REAL(KIND=8),DIMENSION(*),INTENT(IN), OPTIONAL :: grp_cov
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditpooledcovariance(task, grp_indices,    &
     &           pld_mean, pld_cov, grp_cov_indices, grp_means, grp_cov)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                    :: task
              INTEGER,DIMENSION(*),INTENT(IN)      :: grp_indices
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: pld_mean
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: pld_cov
              INTEGER,DIMENSION(*),INTENT(IN),OPTIONAL ::grp_cov_indices
              REAL(KIND=4),DIMENSION(*),INTENT(IN),OPTIONAL :: grp_means
              REAL(KIND=4),DIMENSION(*),INTENT(IN), OPTIONAL :: grp_cov
       END FUNCTION
      END INTERFACE


!  Editors to modify robust variance-covariance matrix
!  related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditrobustcovariance(task,rcov_storage,    &
     &                                      nparams, params,rmean, rcov)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                    :: task
              INTEGER,INTENT(IN)                   :: rcov_storage
              INTEGER,INTENT(IN)                   :: nparams
              REAL(KIND=8),INTENT(IN)              :: params(nparams)
              REAL(KIND=8),DIMENSION(*),INTENT(IN) :: rmean
              REAL(KIND=8),DIMENSION(*),INTENT(IN) :: rcov
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditrobustcovariance(task,rcov_storage,    &
     &                                      nparams, params,rmean, rcov)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                    :: task
              INTEGER,INTENT(IN)                   :: rcov_storage
              INTEGER,INTENT(IN)                   :: nparams
              REAL(KIND=4),INTENT(IN)              :: params(nparams)
              REAL(KIND=4),DIMENSION(*),INTENT(IN) :: rmean
              REAL(KIND=4),DIMENSION(*),INTENT(IN) :: rcov
       END FUNCTION
      END INTERFACE


!  Editors to modify outliers detection parameters
      INTERFACE
        INTEGER FUNCTION vsldsseditoutliersdetection(task,               &
     &                                               nparams, params, w)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)       :: task
              INTEGER,INTENT(IN)      :: nparams
              REAL(KIND=8),DIMENSION(*),INTENT(IN) :: params(nparams)
              REAL(KIND=8),DIMENSION(*),INTENT(IN) :: w
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslssseditoutliersdetection(task,               &
     &                                               nparams, params, w)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)       :: task
              INTEGER,INTENT(IN)      :: nparams
              REAL(KIND=4),DIMENSION(*),INTENT(IN) :: params(nparams)
              REAL(KIND=4),DIMENSION(*),INTENT(IN) :: w
        END FUNCTION
      END INTERFACE

!  Editors to modify missing values parameters
      INTERFACE
        INTEGER FUNCTION vsldsseditmissingvalues(task, nparams, params,  &
     &        init_estimates_n, init_estimates,  prior_n, prior,         &
     &        simul_missing_vals_n,simul_missing_vals,                   &
     &        estimates_n, estimates)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)           :: task
              INTEGER,INTENT(IN)          :: nparams
              REAL(KIND=8),INTENT(IN)     :: params(nparams)
              INTEGER,INTENT(IN),OPTIONAL :: init_estimates_n
              REAL(KIND=8),INTENT(IN), DIMENSION(*),OPTIONAL::           &
     &                                       init_estimates
              INTEGER,INTENT(IN),OPTIONAL :: prior_n
              REAL(KIND=8),INTENT(IN),DIMENSION(*),OPTIONAL :: prior
              INTEGER,INTENT(IN),OPTIONAL :: simul_missing_vals_n
              REAL(KIND=8),INTENT(IN),DIMENSION(*),OPTIONAL ::           &
     &                                 simul_missing_vals
              INTEGER,INTENT(IN),OPTIONAL :: estimates_n
              REAL(KIND=8),INTENT(IN), DIMENSION(*), OPTIONAL ::         &
     &                                 estimates
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslssseditmissingvalues(task, nparams, params,  &
     &        init_estimates_n, init_estimates,  prior_n, prior,         &
     &        simul_missing_vals_n,simul_missing_vals,                   &
     &        estimates_n, estimates)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)           :: task
              INTEGER,INTENT(IN)          :: nparams
              REAL(KIND=4),INTENT(IN)     :: params(nparams)
              INTEGER,INTENT(IN),OPTIONAL :: init_estimates_n
              REAL(KIND=4),INTENT(IN), DIMENSION(*),OPTIONAL::           &
     &                                       init_estimates
              INTEGER,INTENT(IN),OPTIONAL :: prior_n
              REAL(KIND=4),INTENT(IN),DIMENSION(*),OPTIONAL :: prior
              INTEGER,INTENT(IN),OPTIONAL :: simul_missing_vals_n
              REAL(KIND=4),INTENT(IN),DIMENSION(*),OPTIONAL ::           &
     &                                 simul_missing_vals
              INTEGER,INTENT(IN),OPTIONAL :: estimates_n
              REAL(KIND=4),INTENT(IN), DIMENSION(*), OPTIONAL ::         &
     &                                 estimates
        END FUNCTION
      END INTERFACE

!     Editors to modify matrix parameterization
!     related parameters
      INTERFACE
       INTEGER FUNCTION vsldsseditcorparameterization (task,             &
     &                             cor, cor_storage, pcor, pcor_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                    :: task
              REAL(KIND=8),DIMENSION(*),INTENT(IN) :: cor
              INTEGER,INTENT(IN)                   :: cor_storage
              REAL(KIND=8),DIMENSION(*),INTENT(IN) :: pcor
              INTEGER,INTENT(IN)                   :: pcor_storage
       END FUNCTION
      END INTERFACE

      INTERFACE
       INTEGER FUNCTION vslssseditcorparameterization (task,             &
     &                             cor, cor_storage, pcor, pcor_storage)
                  USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)                    :: task
              REAL(KIND=4),DIMENSION(*),INTENT(IN) :: cor
              INTEGER,INTENT(IN)                   :: cor_storage
              REAL(KIND=4),DIMENSION(*),INTENT(IN) :: pcor
              INTEGER,INTENT(IN)                   :: pcor_storage
       END FUNCTION
      END INTERFACE

!  Compute routines
      INTERFACE
        INTEGER FUNCTION vsldsscompute(task, estimates, method)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)            :: task
              INTEGER(KIND=8),INTENT(IN)   :: estimates
              INTEGER,INTENT(IN)           :: method
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslssscompute(task, estimates, method)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)            :: task
              INTEGER(KIND=8),INTENT(IN)   :: estimates
              INTEGER,INTENT(IN)           :: method
        END FUNCTION
      END INTERFACE

!  Task destructor
      INTERFACE
        INTEGER FUNCTION vslssdeletetask(task)
                USE MKL_VSL_TYPE
              TYPE(VSL_SS_TASK)      :: task
        END FUNCTION
      END INTERFACE

      END MODULE MKL_VSL
