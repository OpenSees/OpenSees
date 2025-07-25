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
// adapted from https://caps.gsfc.nasa.gov/simpson/software/m55inv_f90.txt
//
#include <math.h>

int cmx_inv5(double *a, double *ainv, int *ok_flag__)
{
/* ****************************************************************************************
 * m55inv  -  compute the inverse of a 5x5 matrix.
 
 * a       = input 5x5 matrix to be inverted
 * ainv    = output 5x5 inverse of matrix a
 * ok_flag = (output) .true. if the input matrix could be inverted, and .false. if the input matrix is singular.
 ***************************************************************************************** */
    static double cofactor[25];
    static double a11, a12, a13, a14, a15, a21, a22, a23, a24, a25, a31, 
                  a32, a33, a34, a35, a41, a42, a43, a44, a45, a51, a52, a53, a54, 
                  a55, det;

    const double eps = 1e-10;

    /* Parameter adjustments */
    ainv -= 6;
    a -= 6;

    a11 = a[6];
    a12 = a[11];
    a13 = a[16];
    a14 = a[21];
    a15 = a[26];
    a21 = a[7];
    a22 = a[12];
    a23 = a[17];
    a24 = a[22];
    a25 = a[27];
    a31 = a[8];
    a32 = a[13];
    a33 = a[18];
    a34 = a[23];
    a35 = a[28];
    a41 = a[9];
    a42 = a[14];
    a43 = a[19];
    a44 = a[24];
    a45 = a[29];
    a51 = a[10];
    a52 = a[15];
    a53 = a[20];
    a54 = a[25];
    a55 = a[30];
    det = a15*a24*a33*a42*a51 - a14*a25*a33*a42*a51 - a15*
            a23*a34*a42*a51 + a13*a25*a34*a42*a51 + a14*a23*
            a35*a42*a51 - a13*a24*a35*a42*a51 - a15*a24*a32*
            a43*a51 + a14*a25*a32*a43*a51 + a15*a22*a34*a43*
            a51 - a12*a25*a34*a43*a51 - a14*a22*a35*a43*a51 + 
            a12*a24*a35*a43*a51 + a15*a23*a32*a44*a51 - a13*
            a25*a32*a44*a51 - a15*a22*a33*a44*a51 + a12*a25*
            a33*a44*a51 + a13*a22*a35*a44*a51 - a12*a23*a35*
            a44*a51 - a14*a23*a32*a45*a51 + a13*a24*a32*a45*
            a51 + a14*a22*a33*a45*a51 - a12*a24*a33*a45*a51 - 
            a13*a22*a34*a45*a51 + a12*a23*a34*a45*a51 - a15*
            a24*a33*a41*a52 + a14*a25*a33*a41*a52 + a15*a23*
            a34*a41*a52 - a13*a25*a34*a41*a52 - a14*a23*a35*
            a41*a52 + a13*a24*a35*a41*a52 + a15*a24*a31*a43*
            a52 - a14*a25*a31*a43*a52 - a15*a21*a34*a43*a52 + 
            a11*a25*a34*a43*a52 + a14*a21*a35*a43*a52 - a11*
            a24*a35*a43*a52 - a15*a23*a31*a44*a52 + a13*a25*
            a31*a44*a52 + a15*a21*a33*a44*a52 - a11*a25*a33*
            a44*a52 - a13*a21*a35*a44*a52 + a11*a23*a35*a44*
            a52 + a14*a23*a31*a45*a52 - a13*a24*a31*a45*a52 - 
            a14*a21*a33*a45*a52 + a11*a24*a33*a45*a52 + a13*
            a21*a34*a45*a52 - a11*a23*a34*a45*a52 + a15*a24*
            a32*a41*a53 - a14*a25*a32*a41*a53 - a15*a22*a34*
            a41*a53 + a12*a25*a34*a41*a53 + a14*a22*a35*a41*
            a53 - a12*a24*a35*a41*a53 - a15*a24*a31*a42*a53 + 
            a14*a25*a31*a42*a53 + a15*a21*a34*a42*a53 - a11*
            a25*a34*a42*a53 - a14*a21*a35*a42*a53 + a11*a24*
            a35*a42*a53 + a15*a22*a31*a44*a53 - a12*a25*a31*
            a44*a53 - a15*a21*a32*a44*a53 + a11*a25*a32*a44*
            a53 + a12*a21*a35*a44*a53 - a11*a22*a35*a44*a53 - 
            a14*a22*a31*a45*a53 + a12*a24*a31*a45*a53 + a14*
            a21*a32*a45*a53 - a11*a24*a32*a45*a53 - a12*a21*
            a34*a45*a53 + a11*a22*a34*a45*a53 - a15*a23*a32*
            a41*a54 + a13*a25*a32*a41*a54 + a15*a22*a33*a41*
            a54 - a12*a25*a33*a41*a54 - a13*a22*a35*a41*a54 + 
            a12*a23*a35*a41*a54 + a15*a23*a31*a42*a54 - a13*
            a25*a31*a42*a54 - a15*a21*a33*a42*a54 + a11*a25*
            a33*a42*a54 + a13*a21*a35*a42*a54 - a11*a23*a35*
            a42*a54 - a15*a22*a31*a43*a54 + a12*a25*a31*a43*
            a54 + a15*a21*a32*a43*a54 - a11*a25*a32*a43*a54 - 
            a12*a21*a35*a43*a54 + a11*a22*a35*a43*a54 + a13*
            a22*a31*a45*a54 - a12*a23*a31*a45*a54 - a13*a21*
            a32*a45*a54 + a11*a23*a32*a45*a54 + a12*a21*a33*
            a45*a54 - a11*a22*a33*a45*a54 + a14*a23*a32*a41*
            a55 - a13*a24*a32*a41*a55 - a14*a22*a33*a41*a55 + 
            a12*a24*a33*a41*a55 + a13*a22*a34*a41*a55 - a12*
            a23*a34*a41*a55 - a14*a23*a31*a42*a55 + a13*a24*
            a31*a42*a55 + a14*a21*a33*a42*a55 - a11*a24*a33*
            a42*a55 - a13*a21*a34*a42*a55 + a11*a23*a34*a42*
            a55 + a14*a22*a31*a43*a55 - a12*a24*a31*a43*a55 - 
            a14*a21*a32*a43*a55 + a11*a24*a32*a43*a55 + a12*
            a21*a34*a43*a55 - a11*a22*a34*a43*a55 - a13*a22*
            a31*a44*a55 + a12*a23*a31*a44*a55 + a13*a21*a32*
            a44*a55 - a11*a23*a32*a44*a55 - a12*a21*a33*a44*
            a55 + a11*a22*a33*a44*a55;
    if (fabs(det) <= eps) {
        *ok_flag__ = -1;
        return 0;
    }
    cofactor[0] = a25*a34*a43*a52 - a24*a35*a43*a52 - a25*a33*
            a44*a52 + a23*a35*a44*a52 + a24*a33*a45*a52 - a23*
            a34*a45*a52 - a25*a34*a42*a53 + a24*a35*a42*a53 + 
            a25*a32*a44*a53 - a22*a35*a44*a53 - a24*a32*a45*
            a53 + a22*a34*a45*a53 + a25*a33*a42*a54 - a23*a35*
            a42*a54 - a25*a32*a43*a54 + a22*a35*a43*a54 + a23*
            a32*a45*a54 - a22*a33*a45*a54 - a24*a33*a42*a55 + 
            a23*a34*a42*a55 + a24*a32*a43*a55 - a22*a34*a43*
            a55 - a23*a32*a44*a55 + a22*a33*a44*a55;
    cofactor[1] = -a15*a34*a43*a52 + a14*a35*a43*a52 + a15*a33 *
             a44*a52 - a13*a35*a44*a52 - a14*a33*a45*a52 + a13 *
             a34*a45*a52 + a15*a34*a42*a53 - a14*a35*a42*a53 
            - a15*a32*a44*a53 + a12*a35*a44*a53 + a14*a32*a45 
           *a53 - a12*a34*a45*a53 - a15*a33*a42*a54 + a13*a35 
           *a42*a54 + a15*a32*a43*a54 - a12*a35*a43*a54 - a13 
           *a32*a45*a54 + a12*a33*a45*a54 + a14*a33*a42*a55 
            - a13*a34*a42*a55 - a14*a32*a43*a55 + a12*a34*a43 
           *a55 + a13*a32*a44*a55 - a12*a33*a44*a55;
    cofactor[2] = a15*a24*a43*a52 - a14*a25*a43*a52 - a15*a23*
            a44*a52 + a13*a25*a44*a52 + a14*a23*a45*a52 - a13*
            a24*a45*a52 - a15*a24*a42*a53 + a14*a25*a42*a53 + 
            a15*a22*a44*a53 - a12*a25*a44*a53 - a14*a22*a45*
            a53 + a12*a24*a45*a53 + a15*a23*a42*a54 - a13*a25*
            a42*a54 - a15*a22*a43*a54 + a12*a25*a43*a54 + a13*
            a22*a45*a54 - a12*a23*a45*a54 - a14*a23*a42*a55 + 
            a13*a24*a42*a55 + a14*a22*a43*a55 - a12*a24*a43*
            a55 - a13*a22*a44*a55 + a12*a23*a44*a55;
    cofactor[3] = -a15*a24*a33*a52 + a14*a25*a33*a52 + a15*a23 *
             a34*a52 - a13*a25*a34*a52 - a14*a23*a35*a52 + a13 *
             a24*a35*a52 + a15*a24*a32*a53 - a14*a25*a32*a53 
            - a15*a22*a34*a53 + a12*a25*a34*a53 + a14*a22*a35 
           *a53 - a12*a24*a35*a53 - a15*a23*a32*a54 + a13*a25 
           *a32*a54 + a15*a22*a33*a54 - a12*a25*a33*a54 - a13 
           *a22*a35*a54 + a12*a23*a35*a54 + a14*a23*a32*a55 
            - a13*a24*a32*a55 - a14*a22*a33*a55 + a12*a24*a33 
           *a55 + a13*a22*a34*a55 - a12*a23*a34*a55;
    cofactor[4] = a15*a24*a33*a42 - a14*a25*a33*a42 - a15*a23*
            a34*a42 + a13*a25*a34*a42 + a14*a23*a35*a42 - a13*
            a24*a35*a42 - a15*a24*a32*a43 + a14*a25*a32*a43 + 
            a15*a22*a34*a43 - a12*a25*a34*a43 - a14*a22*a35*
            a43 + a12*a24*a35*a43 + a15*a23*a32*a44 - a13*a25*
            a32*a44 - a15*a22*a33*a44 + a12*a25*a33*a44 + a13*
            a22*a35*a44 - a12*a23*a35*a44 - a14*a23*a32*a45 + 
            a13*a24*a32*a45 + a14*a22*a33*a45 - a12*a24*a33*
            a45 - a13*a22*a34*a45 + a12*a23*a34*a45;
    cofactor[5] = -a25*a34*a43*a51 + a24*a35*a43*a51 + a25*a33 *
             a44*a51 - a23*a35*a44*a51 - a24*a33*a45*a51 + a23 *
             a34*a45*a51 + a25*a34*a41*a53 - a24*a35*a41*a53 
            - a25*a31*a44*a53 + a21*a35*a44*a53 + a24*a31*a45 
           *a53 - a21*a34*a45*a53 - a25*a33*a41*a54 + a23*a35 
           *a41*a54 + a25*a31*a43*a54 - a21*a35*a43*a54 - a23 
           *a31*a45*a54 + a21*a33*a45*a54 + a24*a33*a41*a55 
            - a23*a34*a41*a55 - a24*a31*a43*a55 + a21*a34*a43 
           *a55 + a23*a31*a44*a55 - a21*a33*a44*a55;
    cofactor[6] = a15*a34*a43*a51 - a14*a35*a43*a51 - a15*a33*
            a44*a51 + a13*a35*a44*a51 + a14*a33*a45*a51 - a13*
            a34*a45*a51 - a15*a34*a41*a53 + a14*a35*a41*a53 + 
            a15*a31*a44*a53 - a11*a35*a44*a53 - a14*a31*a45*
            a53 + a11*a34*a45*a53 + a15*a33*a41*a54 - a13*a35*
            a41*a54 - a15*a31*a43*a54 + a11*a35*a43*a54 + a13*
            a31*a45*a54 - a11*a33*a45*a54 - a14*a33*a41*a55 + 
            a13*a34*a41*a55 + a14*a31*a43*a55 - a11*a34*a43*
            a55 - a13*a31*a44*a55 + a11*a33*a44*a55;
    cofactor[7] = -a15*a24*a43*a51 + a14*a25*a43*a51 + a15*a23 *
             a44*a51 - a13*a25*a44*a51 - a14*a23*a45*a51 + a13 *
             a24*a45*a51 + a15*a24*a41*a53 - a14*a25*a41*a53 
            - a15*a21*a44*a53 + a11*a25*a44*a53 + a14*a21*a45 
           *a53 - a11*a24*a45*a53 - a15*a23*a41*a54 + a13*a25 
           *a41*a54 + a15*a21*a43*a54 - a11*a25*a43*a54 - a13 
           *a21*a45*a54 + a11*a23*a45*a54 + a14*a23*a41*a55 
            - a13*a24*a41*a55 - a14*a21*a43*a55 + a11*a24*a43 
           *a55 + a13*a21*a44*a55 - a11*a23*a44*a55;
    cofactor[8] = a15*a24*a33*a51 - a14*a25*a33*a51 - a15*a23*
            a34*a51 + a13*a25*a34*a51 + a14*a23*a35*a51 - a13*
            a24*a35*a51 - a15*a24*a31*a53 + a14*a25*a31*a53 + 
            a15*a21*a34*a53 - a11*a25*a34*a53 - a14*a21*a35*
            a53 + a11*a24*a35*a53 + a15*a23*a31*a54 - a13*a25*
            a31*a54 - a15*a21*a33*a54 + a11*a25*a33*a54 + a13*
            a21*a35*a54 - a11*a23*a35*a54 - a14*a23*a31*a55 + 
            a13*a24*a31*a55 + a14*a21*a33*a55 - a11*a24*a33*
            a55 - a13*a21*a34*a55 + a11*a23*a34*a55;
    cofactor[9] = -a15*a24*a33*a41 + a14*a25*a33*a41 + a15*a23 *
             a34*a41 - a13*a25*a34*a41 - a14*a23*a35*a41 + a13 *
             a24*a35*a41 + a15*a24*a31*a43 - a14*a25*a31*a43 
            - a15*a21*a34*a43 + a11*a25*a34*a43 + a14*a21*a35 
           *a43 - a11*a24*a35*a43 - a15*a23*a31*a44 + a13*a25 
           *a31*a44 + a15*a21*a33*a44 - a11*a25*a33*a44 - a13 
           *a21*a35*a44 + a11*a23*a35*a44 + a14*a23*a31*a45 
            - a13*a24*a31*a45 - a14*a21*a33*a45 + a11*a24*a33 
           *a45 + a13*a21*a34*a45 - a11*a23*a34*a45;
    cofactor[10] = a25*a34*a42*a51 - a24*a35*a42*a51 - a25*a32 *
             a44*a51 + a22*a35*a44*a51 + a24*a32*a45*a51 - a22 *
             a34*a45*a51 - a25*a34*a41*a52 + a24*a35*a41*a52 
            + a25*a31*a44*a52 - a21*a35*a44*a52 - a24*a31*a45 
           *a52 + a21*a34*a45*a52 + a25*a32*a41*a54 - a22*a35 
           *a41*a54 - a25*a31*a42*a54 + a21*a35*a42*a54 + a22 
           *a31*a45*a54 - a21*a32*a45*a54 - a24*a32*a41*a55 
            + a22*a34*a41*a55 + a24*a31*a42*a55 - a21*a34*a42 
           *a55 - a22*a31*a44*a55 + a21*a32*a44*a55;
    cofactor[11] = -a15*a34*a42*a51 + a14*a35*a42*a51 + a15*a32 
           *a44*a51 - a12*a35*a44*a51 - a14*a32*a45*a51 + a12 
           *a34*a45*a51 + a15*a34*a41*a52 - a14*a35*a41*a52 
            - a15*a31*a44*a52 + a11*a35*a44*a52 + a14*a31*a45 
           *a52 - a11*a34*a45*a52 - a15*a32*a41*a54 + a12*a35 
           *a41*a54 + a15*a31*a42*a54 - a11*a35*a42*a54 - a12 
           *a31*a45*a54 + a11*a32*a45*a54 + a14*a32*a41*a55 
            - a12*a34*a41*a55 - a14*a31*a42*a55 + a11*a34*a42 
           *a55 + a12*a31*a44*a55 - a11*a32*a44*a55;
    cofactor[12] = a15*a24*a42*a51 - a14*a25*a42*a51 - a15*a22 *
             a44*a51 + a12*a25*a44*a51 + a14*a22*a45*a51 - a12 *
             a24*a45*a51 - a15*a24*a41*a52 + a14*a25*a41*a52 
            + a15*a21*a44*a52 - a11*a25*a44*a52 - a14*a21*a45 
           *a52 + a11*a24*a45*a52 + a15*a22*a41*a54 - a12*a25 
           *a41*a54 - a15*a21*a42*a54 + a11*a25*a42*a54 + a12 
           *a21*a45*a54 - a11*a22*a45*a54 - a14*a22*a41*a55 
            + a12*a24*a41*a55 + a14*a21*a42*a55 - a11*a24*a42 
           *a55 - a12*a21*a44*a55 + a11*a22*a44*a55;
    cofactor[13] = -a15*a24*a32*a51 + a14*a25*a32*a51 + a15*a22 
           *a34*a51 - a12*a25*a34*a51 - a14*a22*a35*a51 + a12 
           *a24*a35*a51 + a15*a24*a31*a52 - a14*a25*a31*a52 
            - a15*a21*a34*a52 + a11*a25*a34*a52 + a14*a21*a35 
           *a52 - a11*a24*a35*a52 - a15*a22*a31*a54 + a12*a25 
           *a31*a54 + a15*a21*a32*a54 - a11*a25*a32*a54 - a12 
           *a21*a35*a54 + a11*a22*a35*a54 + a14*a22*a31*a55 
            - a12*a24*a31*a55 - a14*a21*a32*a55 + a11*a24*a32 
           *a55 + a12*a21*a34*a55 - a11*a22*a34*a55;
    cofactor[14] = a15*a24*a32*a41 - a14*a25*a32*a41 - a15*a22 *
             a34*a41 + a12*a25*a34*a41 + a14*a22*a35*a41 - a12 *
             a24*a35*a41 - a15*a24*a31*a42 + a14*a25*a31*a42 
            + a15*a21*a34*a42 - a11*a25*a34*a42 - a14*a21*a35 
           *a42 + a11*a24*a35*a42 + a15*a22*a31*a44 - a12*a25 
           *a31*a44 - a15*a21*a32*a44 + a11*a25*a32*a44 + a12 
           *a21*a35*a44 - a11*a22*a35*a44 - a14*a22*a31*a45 
            + a12*a24*a31*a45 + a14*a21*a32*a45 - a11*a24*a32 
           *a45 - a12*a21*a34*a45 + a11*a22*a34*a45;
    cofactor[15] = -a25*a33*a42*a51 + a23*a35*a42*a51 + a25*a32 
           *a43*a51 - a22*a35*a43*a51 - a23*a32*a45*a51 + a22 
           *a33*a45*a51 + a25*a33*a41*a52 - a23*a35*a41*a52 
            - a25*a31*a43*a52 + a21*a35*a43*a52 + a23*a31*a45 
           *a52 - a21*a33*a45*a52 - a25*a32*a41*a53 + a22*a35 
           *a41*a53 + a25*a31*a42*a53 - a21*a35*a42*a53 - a22 
           *a31*a45*a53 + a21*a32*a45*a53 + a23*a32*a41*a55 
            - a22*a33*a41*a55 - a23*a31*a42*a55 + a21*a33*a42 
           *a55 + a22*a31*a43*a55 - a21*a32*a43*a55;
    cofactor[16] = a15*a33*a42*a51 - a13*a35*a42*a51 - a15*a32 *
             a43*a51 + a12*a35*a43*a51 + a13*a32*a45*a51 - a12 *
             a33*a45*a51 - a15*a33*a41*a52 + a13*a35*a41*a52 
            + a15*a31*a43*a52 - a11*a35*a43*a52 - a13*a31*a45 
           *a52 + a11*a33*a45*a52 + a15*a32*a41*a53 - a12*a35 
           *a41*a53 - a15*a31*a42*a53 + a11*a35*a42*a53 + a12 
           *a31*a45*a53 - a11*a32*a45*a53 - a13*a32*a41*a55 
            + a12*a33*a41*a55 + a13*a31*a42*a55 - a11*a33*a42 
           *a55 - a12*a31*a43*a55 + a11*a32*a43*a55;
    cofactor[17] = -a15*a23*a42*a51 + a13*a25*a42*a51 + a15*a22 
           *a43*a51 - a12*a25*a43*a51 - a13*a22*a45*a51 + a12 
           *a23*a45*a51 + a15*a23*a41*a52 - a13*a25*a41*a52 
            - a15*a21*a43*a52 + a11*a25*a43*a52 + a13*a21*a45 
           *a52 - a11*a23*a45*a52 - a15*a22*a41*a53 + a12*a25 
           *a41*a53 + a15*a21*a42*a53 - a11*a25*a42*a53 - a12 
           *a21*a45*a53 + a11*a22*a45*a53 + a13*a22*a41*a55 
            - a12*a23*a41*a55 - a13*a21*a42*a55 + a11*a23*a42 
           *a55 + a12*a21*a43*a55 - a11*a22*a43*a55;
    cofactor[18] = a15*a23*a32*a51 - a13*a25*a32*a51 - a15*a22 *
             a33*a51 + a12*a25*a33*a51 + a13*a22*a35*a51 - a12 *
             a23*a35*a51 - a15*a23*a31*a52 + a13*a25*a31*a52 
            + a15*a21*a33*a52 - a11*a25*a33*a52 - a13*a21*a35 
           *a52 + a11*a23*a35*a52 + a15*a22*a31*a53 - a12*a25 
           *a31*a53 - a15*a21*a32*a53 + a11*a25*a32*a53 + a12 
           *a21*a35*a53 - a11*a22*a35*a53 - a13*a22*a31*a55 
            + a12*a23*a31*a55 + a13*a21*a32*a55 - a11*a23*a32 
           *a55 - a12*a21*a33*a55 + a11*a22*a33*a55;
    cofactor[19] = -a15*a23*a32*a41 + a13*a25*a32*a41 + a15*a22 
           *a33*a41 - a12*a25*a33*a41 - a13*a22*a35*a41 + a12 
           *a23*a35*a41 + a15*a23*a31*a42 - a13*a25*a31*a42 
            - a15*a21*a33*a42 + a11*a25*a33*a42 + a13*a21*a35 
           *a42 - a11*a23*a35*a42 - a15*a22*a31*a43 + a12*a25 
           *a31*a43 + a15*a21*a32*a43 - a11*a25*a32*a43 - a12 
           *a21*a35*a43 + a11*a22*a35*a43 + a13*a22*a31*a45 
            - a12*a23*a31*a45 - a13*a21*a32*a45 + a11*a23*a32 
           *a45 + a12*a21*a33*a45 - a11*a22*a33*a45;
    cofactor[20] = a24*a33*a42*a51 - a23*a34*a42*a51 - a24*a32 *
             a43*a51 + a22*a34*a43*a51 + a23*a32*a44*a51 - a22 *
             a33*a44*a51 - a24*a33*a41*a52 + a23*a34*a41*a52 
            + a24*a31*a43*a52 - a21*a34*a43*a52 - a23*a31*a44 
           *a52 + a21*a33*a44*a52 + a24*a32*a41*a53 - a22*a34 
           *a41*a53 - a24*a31*a42*a53 + a21*a34*a42*a53 + a22 
           *a31*a44*a53 - a21*a32*a44*a53 - a23*a32*a41*a54 
            + a22*a33*a41*a54 + a23*a31*a42*a54 - a21*a33*a42 
           *a54 - a22*a31*a43*a54 + a21*a32*a43*a54;
    cofactor[21] = -a14*a33*a42*a51 + a13*a34*a42*a51 + a14*a32 
           *a43*a51 - a12*a34*a43*a51 - a13*a32*a44*a51 + a12 
           *a33*a44*a51 + a14*a33*a41*a52 - a13*a34*a41*a52 
            - a14*a31*a43*a52 + a11*a34*a43*a52 + a13*a31*a44 
           *a52 - a11*a33*a44*a52 - a14*a32*a41*a53 + a12*a34 
           *a41*a53 + a14*a31*a42*a53 - a11*a34*a42*a53 - a12 
           *a31*a44*a53 + a11*a32*a44*a53 + a13*a32*a41*a54 
            - a12*a33*a41*a54 - a13*a31*a42*a54 + a11*a33*a42 
           *a54 + a12*a31*a43*a54 - a11*a32*a43*a54;
    cofactor[22] = a14*a23*a42*a51 - a13*a24*a42*a51 - a14*a22 *
             a43*a51 + a12*a24*a43*a51 + a13*a22*a44*a51 - a12 *
             a23*a44*a51 - a14*a23*a41*a52 + a13*a24*a41*a52 
            + a14*a21*a43*a52 - a11*a24*a43*a52 - a13*a21*a44 
           *a52 + a11*a23*a44*a52 + a14*a22*a41*a53 - a12*a24 
           *a41*a53 - a14*a21*a42*a53 + a11*a24*a42*a53 + a12 
           *a21*a44*a53 - a11*a22*a44*a53 - a13*a22*a41*a54 
            + a12*a23*a41*a54 + a13*a21*a42*a54 - a11*a23*a42 
           *a54 - a12*a21*a43*a54 + a11*a22*a43*a54;
    cofactor[23] = -a14*a23*a32*a51 + a13*a24*a32*a51 + a14*a22 
           *a33*a51 - a12*a24*a33*a51 - a13*a22*a34*a51 + a12 
           *a23*a34*a51 + a14*a23*a31*a52 - a13*a24*a31*a52 
            - a14*a21*a33*a52 + a11*a24*a33*a52 + a13*a21*a34 
           *a52 - a11*a23*a34*a52 - a14*a22*a31*a53 + a12*a24 
           *a31*a53 + a14*a21*a32*a53 - a11*a24*a32*a53 - a12 
           *a21*a34*a53 + a11*a22*a34*a53 + a13*a22*a31*a54 
            - a12*a23*a31*a54 - a13*a21*a32*a54 + a11*a23*a32 
           *a54 + a12*a21*a33*a54 - a11*a22*a33*a54;
    cofactor[24] = a14*a23*a32*a41 - a13*a24*a32*a41 - a14*a22 *
             a33*a41 + a12*a24*a33*a41 + a13*a22*a34*a41 - a12 *
             a23*a34*a41 - a14*a23*a31*a42 + a13*a24*a31*a42 
            + a14*a21*a33*a42 - a11*a24*a33*a42 - a13*a21*a34 
           *a42 + a11*a23*a34*a42 + a14*a22*a31*a43 - a12*a24 
           *a31*a43 - a14*a21*a32*a43 + a11*a24*a32*a43 + a12 
           *a21*a34*a43 - a11*a22*a34*a43 - a13*a22*a31*a44 
            + a12*a23*a31*a44 + a13*a21*a32*a44 - a11*a23*a32 
           *a44 - a12*a21*a33*a44 + a11*a22*a33*a44;

    for (int i__ = 1; i__ <= 5; ++i__)
        for (int j = 1; j <= 5; ++j)
            ainv[j + i__*5] = cofactor[i__ + j*5 - 6] / det;

    *ok_flag__ = 0;
    return 0;
}

