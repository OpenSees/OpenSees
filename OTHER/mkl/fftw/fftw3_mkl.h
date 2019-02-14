/*******************************************************************************
* Copyright 2005-2018 Intel Corporation All Rights Reserved.
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
 * Definitions for FFTW3 wrappers to Intel(R) Math Kernel Library (Intel(R) MKL).
 *
 ******************************************************************************
 */

#ifndef FFTW3_MKL_H
#define FFTW3_MKL_H

#include <stdlib.h>
#include "fftw3.h"
#include "mkl_dfti.h"
#include "mkl_trig_transforms.h"
#include "mkl_service.h"

typedef struct fftw_mkl_plan_s *fftw_mkl_plan;
typedef struct fftw3_mkl_s      fftw3_mkl_s;

/* Plan holder for the wrappers */
struct fftw_mkl_plan_s
{
    DFTI_DESCRIPTOR_HANDLE desc;
    void *io[4];
    MKL_INT *ipar;
    double *dpar;
    float *spar;
    void (*execute) (fftw_mkl_plan p);
    void (*destroy) (fftw_mkl_plan p);
    void *mpi_plan; /* placeholder for FFTW3 MPI Intel(R) MKL plan */
};

/* Global helper structure */
struct fftw3_mkl_s
{
    int verbose;
    int nthreads;
    double timelimit;
    int number_of_user_threads; /* Will be deprecated in nearest future */
    fftw_mkl_plan (*new_plan) (void);
    int default_alignment;
};

FFTW_EXTERN fftw3_mkl_s fftw3_mkl;

#define MKL_MAXRANK 7
#define MKL_ONE     1
#define MKL_RODFT00 413

#define BAD(status) ((status) && !DftiErrorClass((status),DFTI_NO_ERROR))

#ifndef UNUSED
#define UNUSED(p) (void)p
#endif

#endif /* FFTW3_MKL_H */
