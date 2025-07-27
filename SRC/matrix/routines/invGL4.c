//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// adapted from https://caps.gsfc.nasa.gov/simpson/software/m44inv_f90.txt
//
#include <math.h>

int cmx_inv4(double *a, double *ainv, int *ok_flag__)
{
/* **************************************************************************************** 
 *  m44inv  -  compute the inverse of a 4x4 matrix. 

 *  a       = input 4x4 matrix to be inverted 
 *  ainv    = output 4x4 inverse of matrix a 
 *  ok_flag = (output) .true. if the input matrix could be inverted, and .false. if the input matrix is singular. 
 * **************************************************************************************** */
    static double cofactor[16];
    static double eps=1e-10;

    /* Parameter adjustments */
    ainv -= 5;
    a -= 5;

    const double det = a[ 5]*(a[10]*(a[15]*a[20] - a[19]*a[16]) 
                            + a[14]*(a[19]*a[12] - a[11]*a[20]) 
                            + a[18]*(a[11]*a[16] - a[15]*a[12])) 
                     - a[ 9]*(a[ 6]*(a[15]*a[20] - a[19]*a[16])
                            + a[14]*(a[19]*a[ 8] - a[ 7]*a[20])
                            + a[18]*(a[ 7]*a[16] - a[15]*a[ 8]))
                     + a[13]*(a[ 6]*(a[11]*a[20] - a[19]*a[12]) 
                            + a[10]*(a[19]*a[ 8] - a[ 7]*a[20]) 
                            + a[18]*(a[ 7]*a[12] - a[11]*a[ 8]))
                     - a[17]*(a[ 6]*(a[11]*a[16] - a[15]*a[12])
                            + a[10]*(a[15]*a[ 8] - a[ 7]*a[16])
                            + a[14]*(a[ 7]*a[12] - a[11]*a[ 8]));

    if (fabs(det) <= eps) {
        *ok_flag__ = -1;
    }

    cofactor[ 0] = a[10]*(a[15]*a[20] - a[19]*a[16])
                 + a[14]*(a[19]*a[12] - a[11]*a[20]) 
                 + a[18]*(a[11]*a[16] - a[15]*a[12]);
    cofactor[ 4] = a[ 6]*(a[19]*a[16] - a[15]*a[20]) 
                 + a[14]*(a[ 7]*a[20] - a[19]*a[ 8])
                 + a[18]*(a[15]*a[ 8] - a[ 7]*a[16]);
    cofactor[ 8] = a[ 6]*(a[11]*a[20] - a[19]*a[12]) 
                 + a[10]*(a[19]*a[ 8] - a[ 7]*a[20])
                 + a[18]*(a[ 7]*a[12] - a[11]*a[ 8]);
    cofactor[12] = a[ 6]*(a[15]*a[12] - a[11]*a[16])
                 + a[10]*(a[ 7]*a[16] - a[15]*a[ 8])
                 + a[14]*(a[11]*a[ 8] - a[ 7]*a[12]);
    cofactor[ 1] = a[ 9]*(a[19]*a[16] - a[15]*a[20])
                 + a[13]*(a[11]*a[20] - a[19]*a[12])
                 + a[17]*(a[15]*a[12] - a[11]*a[16]);
    cofactor[ 5] = a[ 5]*(a[15]*a[20] - a[19]*a[16])
                 + a[13]*(a[19]*a[ 8] - a[ 7]*a[20])
                 + a[17]*(a[ 7]*a[16] - a[15]*a[ 8]);
    cofactor[ 9] = a[ 5]*(a[19]*a[12] - a[11]*a[20])
                 + a[ 9]*(a[ 7]*a[20] - a[19]*a[ 8])
                 + a[17]*(a[11]*a[ 8] - a[ 7]*a[12]);
    cofactor[13] = a[ 5]*(a[11]*a[16] - a[15]*a[12])
                 + a[ 9]*(a[15]*a[ 8] - a[ 7]*a[16])
                 + a[13]*(a[ 7]*a[12] - a[11]*a[ 8]);
    cofactor[ 2] = a[ 9]*(a[14]*a[20] - a[18]*a[16])
                 + a[13]*(a[18]*a[12] - a[10]*a[20])
                 + a[17]*(a[10]*a[16] - a[14]*a[12]);
    cofactor[ 6] = a[ 5]*(a[18]*a[16] - a[14]*a[20])
                 + a[13]*(a[ 6]*a[20] - a[18]*a[ 8])
                 + a[17]*(a[14]*a[ 8] - a[ 6]*a[16]);
    cofactor[10] = a[ 5]*(a[10]*a[20] - a[18]*a[12])
                 + a[ 9]*(a[18]*a[ 8] - a[ 6]*a[20])
                 + a[17]*(a[ 6]*a[12] - a[10]*a[ 8]);
    cofactor[14] = a[ 5]*(a[14]*a[12] - a[10]*a[16])
                 + a[ 9]*(a[ 6]*a[16] - a[14]*a[ 8])
                 + a[13]*(a[10]*a[ 8] - a[ 6]*a[12]);
    cofactor[ 3] = a[ 9]*(a[18]*a[15] - a[14]*a[19])
                 + a[13]*(a[10]*a[19] - a[18]*a[11])
                 + a[17]*(a[14]*a[11] - a[10]*a[15]);
    cofactor[ 7] = a[ 5]*(a[14]*a[19] - a[18]*a[15])
                 + a[13]*(a[18]*a[ 7] - a[ 6]*a[19])
                 + a[17]*(a[ 6]*a[15] - a[14]*a[ 7]);
    cofactor[11] = a[ 5]*(a[18]*a[11] - a[10]*a[19])
                 + a[ 9]*(a[ 6]*a[19] - a[18]*a[ 7])
                 + a[17]*(a[10]*a[ 7] - a[ 6]*a[11]);
    cofactor[15] = a[ 5]*(a[10]*a[15] - a[14]*a[11])
                 + a[ 9]*(a[14]*a[ 7] - a[ 6]*a[15])
                 + a[13]*(a[ 6]*a[11] - a[10]*a[ 7]);

    for (int i__ = 1; i__ <= 4; ++i__)
        for (int j = 1; j <= 4; ++j)
            ainv[j + (i__ << 2)] = cofactor[i__ + (j << 2) - 5] / det;

    *ok_flag__ = 0;
    return 0;
}

