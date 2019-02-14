/*******************************************************************************
* Copyright 1999-2018 Intel Corporation All Rights Reserved.
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
!      Intel(R) Math Kernel Library (Intel(R) MKL) C/C++ interface for
!	   Cluster Sparse Solver
!******************************************************************************/

#if !defined( __MKL_CLUSTER_SPARSE_SOLVER_H )

#include "mkl_types.h"

#define __MKL_CLUSTER_SPARSE_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void cluster_sparse_solver(
     void *, const MKL_INT *, const MKL_INT *, const MKL_INT *, const MKL_INT *, const MKL_INT *,
     const void *, const MKL_INT *, const MKL_INT *, MKL_INT *, const MKL_INT *, MKL_INT *,
     const MKL_INT *,       void *,       void *, const int *, MKL_INT *);

void CLUSTER_SPARSE_SOLVER(
     void *, const MKL_INT *, const MKL_INT *, const MKL_INT *, const MKL_INT *, const MKL_INT *,
     const void *, const MKL_INT *, const MKL_INT *, MKL_INT *, const MKL_INT *, MKL_INT *,
     const MKL_INT *,       void *,       void *, const int *, MKL_INT *);

void cluster_sparse_solver_64(
     void *, const long long int *, const long long int *, const long long int *, const long long int *, const long long int *,
     const void *, const long long int *, const long long int *, long long int *, const long long int *, long long int *,
     const long long int *,       void *,       void *, const int *, long long int *);

void CLUSTER_SPARSE_SOLVER_64(
     void *, const long long int *, const long long int *, const long long int *, const long long int *, const long long int *,
     const void *, const long long int *, const long long int *, long long int *, const long long int *, long long int *,
     const long long int *,       void *,       void *, const int *, long long int *);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
