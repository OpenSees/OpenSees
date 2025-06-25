//
// adapted from https://caps.gsfc.nasa.gov/simpson/software/m22inv_f90.txt
//
// David Simpson
// Claudio Perez
//
#include <math.h>

int cmx_inv2(double *a, double *ainv, int *ok_flag__)
{
  /* ****************************************************************************************
   *  m22inv  -  compute the inverse of a 2x2 matrix.
   *
   *  a      : (input)  2x2 matrix to be inverted 
   *  ainv   : (output) 2x2 inverse of matrix a 
   *
   *  ok_flag: (output) 0 if the input matrix could be inverted, 
   *           and -1 if the input matrix is singular. 
   *
   * ****************************************************************************************/

    /* Parameter adjustments */
    ainv -= 3;
    a    -= 3;

    /* Function Body */
    const double eps = 1e-10;
    const double det = a[3] * a[6] - a[5] * a[4];
    if (fabs(det) <= eps) {
        *ok_flag__ = -1;
    }

    double cofactor[4];

    cofactor[0] =  a[6];
    cofactor[2] = -a[4];
    cofactor[1] = -a[5];
    cofactor[3] =  a[3];

    for (int i__ = 1; i__ <= 2; ++i__)
      for (int j = 1; j <= 2; ++j)
        ainv[j + (i__ << 1)] = cofactor[i__ + (j << 1) - 3] / det;

    *ok_flag__ = 0;
    return 0;
}
