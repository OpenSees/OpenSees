//
// adapted from https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt
//
#include <math.h>


int cmx_inv3(double *a, double *ainv, int *ok_flag__)
{
/* **************************************************************************************** */
/*  m33inv  -  compute the inverse of a 3x3 matrix. */

/*  a       = input 3x3 matrix to be inverted */
/*  ainv    = output 3x3 inverse of matrix a */
/*  ok_flag = (output) .true. if the input matrix could be inverted, and */
/*            .false. if the input matrix is singular. */
/* **************************************************************************************** */


    /* Parameter adjustments */
    ainv -= 4;
    a -= 4;

    const double det = a[4]*a[8]*a[12] 
                     - a[4]*a[11]*a[9] 
                     - a[7]*a[5]*a[12] 
                     + a[7]*a[11]*a[6] 
                     + a[10]*a[5]*a[9] 
                     - a[10]*a[8]*a[6];

    const double eps=1e-10;
    if (fabs(det) <= eps) {
        *ok_flag__ = -1;
        // return 0;
    }

    double cofactor[9];
    cofactor[0] =    a[8]*a[12] - a[11]*a[9];
    cofactor[3] = -(a[5]*a[12] - a[11]*a[6]);
    cofactor[6] =    a[5]*a[9] - a[8]*a[6];
    cofactor[1] = -(a[7]*a[12] - a[10]*a[9]);
    cofactor[4] =   a[4]*a[12] - a[10]*a[6];
    cofactor[7] = -(a[4]*a[9] - a[7]*a[6]);
    cofactor[2] =   a[7]*a[11] - a[10]*a[8];
    cofactor[5] = -(a[4]*a[11] - a[10]*a[5]);
    cofactor[8] =   a[4]*a[8] - a[7]*a[5];

    for (int i__ = 1; i__ <= 3; ++i__)
        for (int j = 1; j <= 3; ++j)
            ainv[j + i__ * 3] = cofactor[i__ + j * 3 - 4] / det;

    *ok_flag__ = 0;
    return 0;
}

