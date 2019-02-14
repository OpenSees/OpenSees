/*******************************************************************************
* Copyright 2010-2018 Intel Corporation All Rights Reserved.
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
 * Definitions for MPI FFTW3 wrappers to Intel(R) Math Kernel Library (Intel(R) MKL).
 *
 ******************************************************************************
 */

#ifndef FFTW3_MPI_MKL_H
#define FFTW3_MPI_MKL_H

#include "fftw3-mpi.h"

#if defined(MKL_SINGLE)
typedef float real_t;
typedef fftwf_complex complex_t;
#define MPI_PREC MPI_FLOAT
#define MKL_PREC DFTI_SINGLE
#define FFTW_MPI_MANGLE(name) FFTW_MPI_MANGLE_FLOAT(name)
#define FFTW_MANGLE(name) FFTW_MANGLE_FLOAT(name)
#else
typedef double real_t;
typedef fftw_complex complex_t;
#define MPI_PREC MPI_DOUBLE
#define MKL_PREC DFTI_DOUBLE
#define FFTW_MPI_MANGLE(name) FFTW_MPI_MANGLE_DOUBLE(name)
#define FFTW_MANGLE(name) FFTW_MANGLE_DOUBLE(name)
#endif

#include "fftw3_mkl.h"
#include "mkl_cdft.h"

#define WANT_FAST_INPLACE_CLUSTER_FFT 1
/* if WANT_FAST_INPLACE_CLUSTER_FFT set to 1, FFTW3 MPI wrappers internally
 * allocate additional memory(workspace) needed for fast inplace Intel(R) MKL CDFT
 * otherwise, no additional memory is used, though the perfomance would be
 * worse, because of many MPI communications */

#endif /* FFTW3_MPI_MKL_H */

