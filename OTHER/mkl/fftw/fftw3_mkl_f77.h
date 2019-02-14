/*******************************************************************************
* Copyright 2006-2018 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
 *
 * Fortran 77 declarations for FFTW3 wrappers to Intel(R) Math Kernel Library
 * (Intel(R) MKL).
 *
 ******************************************************************************
 */

#ifndef FFTW3_MKL_F77_H
#define FFTW3_MKL_F77_H

/* We rely on C interface of FFTW3 wrappers to Intel(R) MKL */
#include "fftw3_mkl.h"

/* Fortran type names */
#define INTEGER    MKL_INT
#define INTEGER4   int
#define INTEGER8   MKL_INT64
#define COMPLEX8   fftwf_complex
#define COMPLEX16  fftw_complex
#define REAL4      float
#define REAL8      double

/* Dummy type names for empty lfftw_* wrappers */
#define COMPLEX32  fftwl_complex
#define REAL16     long double

/* Plan shall be stored in INTEGER*8 variable */
#define PLAN       INTEGER8

#if defined(_FNAME_UPPERCASE)

#define dfftw_cleanup                  DFFTW_CLEANUP
#define dfftw_cleanup_threads          DFFTW_CLEANUP_THREADS
#define dfftw_destroy_plan             DFFTW_DESTROY_PLAN
#define dfftw_execute                  DFFTW_EXECUTE
#define dfftw_execute_dft              DFFTW_EXECUTE_DFT
#define dfftw_execute_dft_c2r          DFFTW_EXECUTE_DFT_C2R
#define dfftw_execute_dft_r2c          DFFTW_EXECUTE_DFT_R2C
#define dfftw_execute_r2r              DFFTW_EXECUTE_R2R
#define dfftw_execute_split_dft        DFFTW_EXECUTE_SPLIT_DFT
#define dfftw_execute_split_dft_c2r    DFFTW_EXECUTE_SPLIT_DFT_C2R
#define dfftw_execute_split_dft_r2c    DFFTW_EXECUTE_SPLIT_DFT_R2C
#define dfftw_export_wisdom            DFFTW_EXPORT_WISDOM
#define dfftw_flops                    DFFTW_FLOPS
#define dfftw_forget_wisdom            DFFTW_FORGET_WISDOM
#define dfftw_import_system_wisdom     DFFTW_IMPORT_SYSTEM_WISDOM
#define dfftw_import_wisdom            DFFTW_IMPORT_WISDOM
#define dfftw_init_threads             DFFTW_INIT_THREADS
#define dfftw_plan_dft                 DFFTW_PLAN_DFT
#define dfftw_plan_dft_1d              DFFTW_PLAN_DFT_1D
#define dfftw_plan_dft_2d              DFFTW_PLAN_DFT_2D
#define dfftw_plan_dft_3d              DFFTW_PLAN_DFT_3D
#define dfftw_plan_dft_c2r             DFFTW_PLAN_DFT_C2R
#define dfftw_plan_dft_c2r_1d          DFFTW_PLAN_DFT_C2R_1D
#define dfftw_plan_dft_c2r_2d          DFFTW_PLAN_DFT_C2R_2D
#define dfftw_plan_dft_c2r_3d          DFFTW_PLAN_DFT_C2R_3D
#define dfftw_plan_dft_r2c             DFFTW_PLAN_DFT_R2C
#define dfftw_plan_dft_r2c_1d          DFFTW_PLAN_DFT_R2C_1D
#define dfftw_plan_dft_r2c_2d          DFFTW_PLAN_DFT_R2C_2D
#define dfftw_plan_dft_r2c_3d          DFFTW_PLAN_DFT_R2C_3D
#define dfftw_plan_guru_dft            DFFTW_PLAN_GURU_DFT
#define dfftw_plan_guru_dft_c2r        DFFTW_PLAN_GURU_DFT_C2R
#define dfftw_plan_guru_dft_r2c        DFFTW_PLAN_GURU_DFT_R2C
#define dfftw_plan_guru_r2r            DFFTW_PLAN_GURU_R2R
#define dfftw_plan_guru_split_dft      DFFTW_PLAN_GURU_SPLIT_DFT
#define dfftw_plan_guru_split_dft_c2r  DFFTW_PLAN_GURU_SPLIT_DFT_C2R
#define dfftw_plan_guru_split_dft_r2c  DFFTW_PLAN_GURU_SPLIT_DFT_R2C
#define dfftw_plan_many_dft            DFFTW_PLAN_MANY_DFT
#define dfftw_plan_many_dft_c2r        DFFTW_PLAN_MANY_DFT_C2R
#define dfftw_plan_many_dft_r2c        DFFTW_PLAN_MANY_DFT_R2C
#define dfftw_plan_many_r2r            DFFTW_PLAN_MANY_R2R
#define dfftw_plan_r2r                 DFFTW_PLAN_R2R
#define dfftw_plan_r2r_1d              DFFTW_PLAN_R2R_1D
#define dfftw_plan_r2r_2d              DFFTW_PLAN_R2R_2D
#define dfftw_plan_r2r_3d              DFFTW_PLAN_R2R_3D
#define dfftw_plan_with_nthreads       DFFTW_PLAN_WITH_NTHREADS
#define dfftw_print_plan               DFFTW_PRINT_PLAN

#define sfftw_cleanup                  SFFTW_CLEANUP
#define sfftw_cleanup_threads          SFFTW_CLEANUP_THREADS
#define sfftw_destroy_plan             SFFTW_DESTROY_PLAN
#define sfftw_execute                  SFFTW_EXECUTE
#define sfftw_execute_dft              SFFTW_EXECUTE_DFT
#define sfftw_execute_dft_c2r          SFFTW_EXECUTE_DFT_C2R
#define sfftw_execute_dft_r2c          SFFTW_EXECUTE_DFT_R2C
#define sfftw_execute_r2r              SFFTW_EXECUTE_R2R
#define sfftw_execute_split_dft        SFFTW_EXECUTE_SPLIT_DFT
#define sfftw_execute_split_dft_c2r    SFFTW_EXECUTE_SPLIT_DFT_C2R
#define sfftw_execute_split_dft_r2c    SFFTW_EXECUTE_SPLIT_DFT_R2C
#define sfftw_export_wisdom            SFFTW_EXPORT_WISDOM
#define sfftw_flops                    SFFTW_FLOPS
#define sfftw_forget_wisdom            SFFTW_FORGET_WISDOM
#define sfftw_import_system_wisdom     SFFTW_IMPORT_SYSTEM_WISDOM
#define sfftw_import_wisdom            SFFTW_IMPORT_WISDOM
#define sfftw_init_threads             SFFTW_INIT_THREADS
#define sfftw_plan_dft                 SFFTW_PLAN_DFT
#define sfftw_plan_dft_1d              SFFTW_PLAN_DFT_1D
#define sfftw_plan_dft_2d              SFFTW_PLAN_DFT_2D
#define sfftw_plan_dft_3d              SFFTW_PLAN_DFT_3D
#define sfftw_plan_dft_c2r             SFFTW_PLAN_DFT_C2R
#define sfftw_plan_dft_c2r_1d          SFFTW_PLAN_DFT_C2R_1D
#define sfftw_plan_dft_c2r_2d          SFFTW_PLAN_DFT_C2R_2D
#define sfftw_plan_dft_c2r_3d          SFFTW_PLAN_DFT_C2R_3D
#define sfftw_plan_dft_r2c             SFFTW_PLAN_DFT_R2C
#define sfftw_plan_dft_r2c_1d          SFFTW_PLAN_DFT_R2C_1D
#define sfftw_plan_dft_r2c_2d          SFFTW_PLAN_DFT_R2C_2D
#define sfftw_plan_dft_r2c_3d          SFFTW_PLAN_DFT_R2C_3D
#define sfftw_plan_guru_dft            SFFTW_PLAN_GURU_DFT
#define sfftw_plan_guru_dft_c2r        SFFTW_PLAN_GURU_DFT_C2R
#define sfftw_plan_guru_dft_r2c        SFFTW_PLAN_GURU_DFT_R2C
#define sfftw_plan_guru_r2r            SFFTW_PLAN_GURU_R2R
#define sfftw_plan_guru_split_dft      SFFTW_PLAN_GURU_SPLIT_DFT
#define sfftw_plan_guru_split_dft_c2r  SFFTW_PLAN_GURU_SPLIT_DFT_C2R
#define sfftw_plan_guru_split_dft_r2c  SFFTW_PLAN_GURU_SPLIT_DFT_R2C
#define sfftw_plan_many_dft            SFFTW_PLAN_MANY_DFT
#define sfftw_plan_many_dft_c2r        SFFTW_PLAN_MANY_DFT_C2R
#define sfftw_plan_many_dft_r2c        SFFTW_PLAN_MANY_DFT_R2C
#define sfftw_plan_many_r2r            SFFTW_PLAN_MANY_R2R
#define sfftw_plan_r2r                 SFFTW_PLAN_R2R
#define sfftw_plan_r2r_1d              SFFTW_PLAN_R2R_1D
#define sfftw_plan_r2r_2d              SFFTW_PLAN_R2R_2D
#define sfftw_plan_r2r_3d              SFFTW_PLAN_R2R_3D
#define sfftw_plan_with_nthreads       SFFTW_PLAN_WITH_NTHREADS
#define sfftw_print_plan               SFFTW_PRINT_PLAN

#else /* i.e. lowercase */

#if defined(_FNAME_SECOND_UNDERSCORE)
#define N(n) n##__
#elif defined(_FNAME_NOUNDERSCORE)
#define N(n) n
#else
#define N(n) n##_
#endif

#define dfftw_cleanup                  N(dfftw_cleanup)
#define dfftw_cleanup_threads          N(dfftw_cleanup_threads)
#define dfftw_destroy_plan             N(dfftw_destroy_plan)
#define dfftw_execute                  N(dfftw_execute)
#define dfftw_execute_dft              N(dfftw_execute_dft)
#define dfftw_execute_dft_c2r          N(dfftw_execute_dft_c2r)
#define dfftw_execute_dft_r2c          N(dfftw_execute_dft_r2c)
#define dfftw_execute_r2r              N(dfftw_execute_r2r)
#define dfftw_execute_split_dft        N(dfftw_execute_split_dft)
#define dfftw_execute_split_dft_c2r    N(dfftw_execute_split_dft_c2r)
#define dfftw_execute_split_dft_r2c    N(dfftw_execute_split_dft_r2c)
#define dfftw_export_wisdom            N(dfftw_export_wisdom)
#define dfftw_flops                    N(dfftw_flops)
#define dfftw_forget_wisdom            N(dfftw_forget_wisdom)
#define dfftw_import_system_wisdom     N(dfftw_import_system_wisdom)
#define dfftw_import_wisdom            N(dfftw_import_wisdom)
#define dfftw_init_threads             N(dfftw_init_threads)
#define dfftw_plan_dft                 N(dfftw_plan_dft)
#define dfftw_plan_dft_1d              N(dfftw_plan_dft_1d)
#define dfftw_plan_dft_2d              N(dfftw_plan_dft_2d)
#define dfftw_plan_dft_3d              N(dfftw_plan_dft_3d)
#define dfftw_plan_dft_c2r             N(dfftw_plan_dft_c2r)
#define dfftw_plan_dft_c2r_1d          N(dfftw_plan_dft_c2r_1d)
#define dfftw_plan_dft_c2r_2d          N(dfftw_plan_dft_c2r_2d)
#define dfftw_plan_dft_c2r_3d          N(dfftw_plan_dft_c2r_3d)
#define dfftw_plan_dft_r2c             N(dfftw_plan_dft_r2c)
#define dfftw_plan_dft_r2c_1d          N(dfftw_plan_dft_r2c_1d)
#define dfftw_plan_dft_r2c_2d          N(dfftw_plan_dft_r2c_2d)
#define dfftw_plan_dft_r2c_3d          N(dfftw_plan_dft_r2c_3d)
#define dfftw_plan_guru_dft            N(dfftw_plan_guru_dft)
#define dfftw_plan_guru_dft_c2r        N(dfftw_plan_guru_dft_c2r)
#define dfftw_plan_guru_dft_r2c        N(dfftw_plan_guru_dft_r2c)
#define dfftw_plan_guru_r2r            N(dfftw_plan_guru_r2r)
#define dfftw_plan_guru_split_dft      N(dfftw_plan_guru_split_dft)
#define dfftw_plan_guru_split_dft_c2r  N(dfftw_plan_guru_split_dft_c2r)
#define dfftw_plan_guru_split_dft_r2c  N(dfftw_plan_guru_split_dft_r2c)
#define dfftw_plan_many_dft            N(dfftw_plan_many_dft)
#define dfftw_plan_many_dft_c2r        N(dfftw_plan_many_dft_c2r)
#define dfftw_plan_many_dft_r2c        N(dfftw_plan_many_dft_r2c)
#define dfftw_plan_many_r2r            N(dfftw_plan_many_r2r)
#define dfftw_plan_r2r                 N(dfftw_plan_r2r)
#define dfftw_plan_r2r_1d              N(dfftw_plan_r2r_1d)
#define dfftw_plan_r2r_2d              N(dfftw_plan_r2r_2d)
#define dfftw_plan_r2r_3d              N(dfftw_plan_r2r_3d)
#define dfftw_plan_with_nthreads       N(dfftw_plan_with_nthreads)
#define dfftw_print_plan               N(dfftw_print_plan)

#define sfftw_cleanup                  N(sfftw_cleanup)
#define sfftw_cleanup_threads          N(sfftw_cleanup_threads)
#define sfftw_destroy_plan             N(sfftw_destroy_plan)
#define sfftw_execute                  N(sfftw_execute)
#define sfftw_execute_dft              N(sfftw_execute_dft)
#define sfftw_execute_dft_c2r          N(sfftw_execute_dft_c2r)
#define sfftw_execute_dft_r2c          N(sfftw_execute_dft_r2c)
#define sfftw_execute_r2r              N(sfftw_execute_r2r)
#define sfftw_execute_split_dft        N(sfftw_execute_split_dft)
#define sfftw_execute_split_dft_c2r    N(sfftw_execute_split_dft_c2r)
#define sfftw_execute_split_dft_r2c    N(sfftw_execute_split_dft_r2c)
#define sfftw_export_wisdom            N(sfftw_export_wisdom)
#define sfftw_flops                    N(sfftw_flops)
#define sfftw_forget_wisdom            N(sfftw_forget_wisdom)
#define sfftw_import_system_wisdom     N(sfftw_import_system_wisdom)
#define sfftw_import_wisdom            N(sfftw_import_wisdom)
#define sfftw_init_threads             N(sfftw_init_threads)
#define sfftw_plan_dft                 N(sfftw_plan_dft)
#define sfftw_plan_dft_1d              N(sfftw_plan_dft_1d)
#define sfftw_plan_dft_2d              N(sfftw_plan_dft_2d)
#define sfftw_plan_dft_3d              N(sfftw_plan_dft_3d)
#define sfftw_plan_dft_c2r             N(sfftw_plan_dft_c2r)
#define sfftw_plan_dft_c2r_1d          N(sfftw_plan_dft_c2r_1d)
#define sfftw_plan_dft_c2r_2d          N(sfftw_plan_dft_c2r_2d)
#define sfftw_plan_dft_c2r_3d          N(sfftw_plan_dft_c2r_3d)
#define sfftw_plan_dft_r2c             N(sfftw_plan_dft_r2c)
#define sfftw_plan_dft_r2c_1d          N(sfftw_plan_dft_r2c_1d)
#define sfftw_plan_dft_r2c_2d          N(sfftw_plan_dft_r2c_2d)
#define sfftw_plan_dft_r2c_3d          N(sfftw_plan_dft_r2c_3d)
#define sfftw_plan_guru_dft            N(sfftw_plan_guru_dft)
#define sfftw_plan_guru_dft_c2r        N(sfftw_plan_guru_dft_c2r)
#define sfftw_plan_guru_dft_r2c        N(sfftw_plan_guru_dft_r2c)
#define sfftw_plan_guru_r2r            N(sfftw_plan_guru_r2r)
#define sfftw_plan_guru_split_dft      N(sfftw_plan_guru_split_dft)
#define sfftw_plan_guru_split_dft_c2r  N(sfftw_plan_guru_split_dft_c2r)
#define sfftw_plan_guru_split_dft_r2c  N(sfftw_plan_guru_split_dft_r2c)
#define sfftw_plan_many_dft            N(sfftw_plan_many_dft)
#define sfftw_plan_many_dft_c2r        N(sfftw_plan_many_dft_c2r)
#define sfftw_plan_many_dft_r2c        N(sfftw_plan_many_dft_r2c)
#define sfftw_plan_many_r2r            N(sfftw_plan_many_r2r)
#define sfftw_plan_r2r                 N(sfftw_plan_r2r)
#define sfftw_plan_r2r_1d              N(sfftw_plan_r2r_1d)
#define sfftw_plan_r2r_2d              N(sfftw_plan_r2r_2d)
#define sfftw_plan_r2r_3d              N(sfftw_plan_r2r_3d)
#define sfftw_plan_with_nthreads       N(sfftw_plan_with_nthreads)
#define sfftw_print_plan               N(sfftw_print_plan)

#endif

FFTW_EXTERN void dfftw_cleanup(void);
FFTW_EXTERN void dfftw_cleanup_threads(void);
FFTW_EXTERN void dfftw_destroy_plan(PLAN*);
FFTW_EXTERN void dfftw_execute(PLAN*);
FFTW_EXTERN void dfftw_execute_dft(PLAN*,COMPLEX16*,COMPLEX16*);
FFTW_EXTERN void dfftw_execute_dft_c2r(PLAN*,COMPLEX16*,REAL8*);
FFTW_EXTERN void dfftw_execute_dft_r2c(PLAN*,REAL8*,COMPLEX16*);
FFTW_EXTERN void dfftw_execute_r2r(PLAN*,REAL8*,REAL8*);
FFTW_EXTERN void dfftw_execute_split_dft(PLAN*,REAL8*,REAL8*,REAL8*,REAL8*);
FFTW_EXTERN void dfftw_execute_split_dft_c2r(PLAN*,REAL8*,REAL8*,REAL8*);
FFTW_EXTERN void dfftw_execute_split_dft_r2c(PLAN*,REAL8*,REAL8*,REAL8*);
FFTW_EXTERN void dfftw_export_wisdom(void*,void*);
FFTW_EXTERN void dfftw_flops(PLAN*,double*,double*,double*);
FFTW_EXTERN void dfftw_forget_wisdom(void);
FFTW_EXTERN void dfftw_import_system_wisdom(INTEGER*);
FFTW_EXTERN void dfftw_import_wisdom(INTEGER*,void*,void*);
FFTW_EXTERN void dfftw_init_threads(INTEGER*);
FFTW_EXTERN void dfftw_plan_dft(PLAN*,INTEGER*,INTEGER*,COMPLEX16*,COMPLEX16*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_1d(PLAN*,INTEGER*,COMPLEX16*,COMPLEX16*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_2d(PLAN*,INTEGER*,INTEGER*,COMPLEX16*,COMPLEX16*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_3d(PLAN*,INTEGER*,INTEGER*,INTEGER*,COMPLEX16*,COMPLEX16*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_c2r(PLAN*,INTEGER*,INTEGER*,COMPLEX16*,REAL8*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_c2r_1d(PLAN*,INTEGER*,COMPLEX16*,REAL8*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_c2r_2d(PLAN*,INTEGER*,INTEGER*,COMPLEX16*,REAL8*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_c2r_3d(PLAN*,INTEGER*,INTEGER*,INTEGER*,COMPLEX16*,REAL8*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_r2c(PLAN*,INTEGER*,INTEGER*,REAL8*,COMPLEX16*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_r2c_1d(PLAN*,INTEGER*,REAL8*,COMPLEX16*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_r2c_2d(PLAN*,INTEGER*,INTEGER*,REAL8*,COMPLEX16*,INTEGER*);
FFTW_EXTERN void dfftw_plan_dft_r2c_3d(PLAN*,INTEGER*,INTEGER*,INTEGER*,REAL8*,COMPLEX16*,INTEGER*);
FFTW_EXTERN void dfftw_plan_guru_dft(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,COMPLEX16*,COMPLEX16*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_guru_dft_c2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,COMPLEX16*,REAL8*,INTEGER*);
FFTW_EXTERN void dfftw_plan_guru_dft_r2c(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL8*,COMPLEX16*,INTEGER*);
FFTW_EXTERN void dfftw_plan_guru_r2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL8*,REAL8*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_guru_split_dft(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL8*,REAL8*,REAL8*,REAL8*,INTEGER*);
FFTW_EXTERN void dfftw_plan_guru_split_dft_c2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL8*,REAL8*,REAL8*,INTEGER*);
FFTW_EXTERN void dfftw_plan_guru_split_dft_r2c(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL8*,REAL8*,REAL8*,INTEGER*);
FFTW_EXTERN void dfftw_plan_many_dft(PLAN*,INTEGER*,INTEGER*,INTEGER*,COMPLEX16*,INTEGER*,INTEGER*,INTEGER*,COMPLEX16*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_many_dft_c2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,COMPLEX16*,INTEGER*,INTEGER*,INTEGER*,REAL8*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_many_dft_r2c(PLAN*,INTEGER*,INTEGER*,INTEGER*,REAL8*,INTEGER*,INTEGER*,INTEGER*,COMPLEX16*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_many_r2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,REAL8*,INTEGER*,INTEGER*,INTEGER*,REAL8*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_r2r(PLAN*,INTEGER*,INTEGER*,REAL8*,REAL8*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_r2r_1d(PLAN*,INTEGER*,REAL8*,REAL8*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_r2r_2d(PLAN*,INTEGER*,INTEGER*,REAL8*,REAL8*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_r2r_3d(PLAN*,INTEGER*,INTEGER*,INTEGER*,REAL8*,REAL8*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void dfftw_plan_with_nthreads(INTEGER*);
FFTW_EXTERN void dfftw_print_plan(PLAN*);

FFTW_EXTERN void sfftw_cleanup(void);
FFTW_EXTERN void sfftw_cleanup_threads(void);
FFTW_EXTERN void sfftw_destroy_plan(PLAN*);
FFTW_EXTERN void sfftw_execute(PLAN*);
FFTW_EXTERN void sfftw_execute_dft(PLAN*,COMPLEX8*,COMPLEX8*);
FFTW_EXTERN void sfftw_execute_dft_c2r(PLAN*,COMPLEX8*,REAL4*);
FFTW_EXTERN void sfftw_execute_dft_r2c(PLAN*,REAL4*,COMPLEX8*);
FFTW_EXTERN void sfftw_execute_r2r(PLAN*,REAL4*,REAL4*);
FFTW_EXTERN void sfftw_execute_split_dft(PLAN*,REAL4*,REAL4*,REAL4*,REAL4*);
FFTW_EXTERN void sfftw_execute_split_dft_c2r(PLAN*,REAL4*,REAL4*,REAL4*);
FFTW_EXTERN void sfftw_execute_split_dft_r2c(PLAN*,REAL4*,REAL4*,REAL4*);
FFTW_EXTERN void sfftw_export_wisdom(void*,void*);
FFTW_EXTERN void sfftw_flops(PLAN*,double*,double*,double*);
FFTW_EXTERN void sfftw_forget_wisdom(void);
FFTW_EXTERN void sfftw_import_system_wisdom(INTEGER*);
FFTW_EXTERN void sfftw_import_wisdom(INTEGER*,void*,void*);
FFTW_EXTERN void sfftw_init_threads(INTEGER*);
FFTW_EXTERN void sfftw_plan_dft(PLAN*,INTEGER*,INTEGER*,COMPLEX8*,COMPLEX8*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_1d(PLAN*,INTEGER*,COMPLEX8*,COMPLEX8*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_2d(PLAN*,INTEGER*,INTEGER*,COMPLEX8*,COMPLEX8*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_3d(PLAN*,INTEGER*,INTEGER*,INTEGER*,COMPLEX8*,COMPLEX8*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_c2r(PLAN*,INTEGER*,INTEGER*,COMPLEX8*,REAL4*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_c2r_1d(PLAN*,INTEGER*,COMPLEX8*,REAL4*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_c2r_2d(PLAN*,INTEGER*,INTEGER*,COMPLEX8*,REAL4*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_c2r_3d(PLAN*,INTEGER*,INTEGER*,INTEGER*,COMPLEX8*,REAL4*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_r2c(PLAN*,INTEGER*,INTEGER*,REAL4*,COMPLEX8*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_r2c_1d(PLAN*,INTEGER*,REAL4*,COMPLEX8*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_r2c_2d(PLAN*,INTEGER*,INTEGER*,REAL4*,COMPLEX8*,INTEGER*);
FFTW_EXTERN void sfftw_plan_dft_r2c_3d(PLAN*,INTEGER*,INTEGER*,INTEGER*,REAL4*,COMPLEX8*,INTEGER*);
FFTW_EXTERN void sfftw_plan_guru_dft(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,COMPLEX8*,COMPLEX8*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_guru_dft_c2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,COMPLEX8*,REAL4*,INTEGER*);
FFTW_EXTERN void sfftw_plan_guru_dft_r2c(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL4*,COMPLEX8*,INTEGER*);
FFTW_EXTERN void sfftw_plan_guru_r2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL4*,REAL4*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_guru_split_dft(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL4*,REAL4*,REAL4*,REAL4*,INTEGER*);
FFTW_EXTERN void sfftw_plan_guru_split_dft_c2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL4*,REAL4*,REAL4*,INTEGER*);
FFTW_EXTERN void sfftw_plan_guru_split_dft_r2c(PLAN*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,REAL4*,REAL4*,REAL4*,INTEGER*);
FFTW_EXTERN void sfftw_plan_many_dft(PLAN*,INTEGER*,INTEGER*,INTEGER*,COMPLEX8*,INTEGER*,INTEGER*,INTEGER*,COMPLEX8*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_many_dft_c2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,COMPLEX8*,INTEGER*,INTEGER*,INTEGER*,REAL4*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_many_dft_r2c(PLAN*,INTEGER*,INTEGER*,INTEGER*,REAL4*,INTEGER*,INTEGER*,INTEGER*,COMPLEX8*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_many_r2r(PLAN*,INTEGER*,INTEGER*,INTEGER*,REAL4*,INTEGER*,INTEGER*,INTEGER*,REAL4*,INTEGER*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_r2r(PLAN*,INTEGER*,INTEGER*,REAL4*,REAL4*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_r2r_1d(PLAN*,INTEGER*,REAL4*,REAL4*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_r2r_2d(PLAN*,INTEGER*,INTEGER*,REAL4*,REAL4*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_r2r_3d(PLAN*,INTEGER*,INTEGER*,INTEGER*,REAL4*,REAL4*,INTEGER*,INTEGER*,INTEGER*,INTEGER*);
FFTW_EXTERN void sfftw_plan_with_nthreads(INTEGER*);
FFTW_EXTERN void sfftw_print_plan(PLAN*);

#endif /* FFTW3_MKL_F77_H */
