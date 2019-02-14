/*******************************************************************************
* Copyright 2002-2018 Intel Corporation All Rights Reserved.
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
!  Content:
!      Intel(R) Math Kernel Library (Intel(R) MKL) interface for Cluster DFT routines
!******************************************************************************/

/* Avoid multiple definition */
#ifndef _MKL_CDFT_TYPES_H_
#define _MKL_CDFT_TYPES_H_

/* Include header-files */
#include "mkl_dfti.h"

/* Keep C++ compilers from getting confused */
#ifdef __cplusplus
extern "C" {
#endif

/* Codes of errors */
#define CDFT_MPI_ERROR      1000
#define CDFT_SPREAD_ERROR   1001

/* Codes of parameters for DftiGetValueDM / DftiSetValueDM */
enum CDFT_CONFIG_PARAM {
    CDFT_LOCAL_SIZE         =1000,
    CDFT_LOCAL_X_START      =1001,
    CDFT_LOCAL_NX           =1002,
    CDFT_MPI_COMM           =1003,
    CDFT_WORKSPACE          =1004,
    CDFT_LOCAL_OUT_X_START  =1005,
    CDFT_LOCAL_OUT_NX       =1006
};

/* Definition of handle to descriptor */
typedef struct _DFTI_DESCRIPTOR_DM* DFTI_DESCRIPTOR_DM_HANDLE;

/* Keep C++ compilers from getting confused (extern "C" {) */
#ifdef __cplusplus
}
#endif

/* Avoid multiple definition (#ifndef _MKL_CDFT_TYPES_H_) */
#endif

