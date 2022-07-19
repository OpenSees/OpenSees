/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

//This file contains the implementation for NURBS derivatives
// Written originally by Vinh Phu Nguyen, nvinhphu@gmail.com

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <iostream>

#include <NurbsDers.h>


#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Subdomain.h>
#include <OPS_Globals.h>


#include <vector>



#define TOL 100*DBL_EPSILON


double** init2DArray(int x, int y)
{
  double **array = (double **)malloc(x * sizeof(double *));

  int c;
  for (c = 0; c < x; c++)
  {
    array[c] = (double*)malloc(y * sizeof(double));
  }

  int d;
  return array;
}

void free2Darray(double **array, int x)
{
  int c;
  for (c = 0; c < x; c++)
  {
    free(array[c]);
  }
  free(array);
}

// int FindSpan(int n, int p, double u, Vector U)
int FindSpan(int n, int p, double u, Vector& U)
{
  /* This function determines the knot span.
    ie. if we have a coordinate u which lies in the range u \in [u_i, u_{i+1})
    we want to find i
   Note that: u_i <= u < (not equal) u_{i+1}!!!
   If we have knot = [0,0.5,1] then u=0.5 has span=1 not 0!!!
  */

  if ( u >= U[n + 1] )
    return n;
  if ( u <= U[p] )
    return p;

  int low = p, high = n + 1, mid = (low + high) / 2;
  while ( u < U[mid] || u >= U[mid + 1] )
  {
    if ( u < U[mid] )
      high = mid;
    else
      low = mid;
    mid = (low + high) / 2;
  }
  return mid;
}



// void BasisFuns( int i, double u, int p, Vector U, double* N)
void BasisFuns( int i, double u, int p, Vector& U, Vector& N)
{

  /*
   we can compute the non zero basis functions
     at point u, there are p+1 non zero basis functions
  */

  N[0] = 1.0;

  int j, r;
  // double *left  = (double *)malloc(sizeof(double) * (p + 1));
  // double *right = (double *)malloc(sizeof(double) * (p + 1));

  static Vector left(p+1);
  static Vector right(p+1);
  left.resize(p+1);
  right.resize(p+1);

  double saved, temp;

  for ( j = 1; j <= p; ++j)
  {
    left[j]  = u - U[i + 1 - j];
    right[j] = U[i + j] - u;
    saved = 0.0;
    for (r = 0; r < j; ++r)
    {
      temp  = N[r] / ( right[r + 1] + left[j - r] );
      N[r]  = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    N[j] = saved;
  }

  // free(left);
  // free(right);

}

// void dersBasisFuns(int i, double u, int p, int order, Vector knot, double **ders)
void dersBasisFuns(int i, double u, int p, int order, Vector& knot, Matrix& ders)
{
  /*
    Calculate the non-zero derivatives of the b-spline functions
  */

  double saved, temp;
  int j, k, j1, j2, r;

  // double *left  = (double *)malloc(sizeof(double) * (p + 1));
  // double *right = (double *)malloc(sizeof(double) * (p + 1));

  static Vector right(p+1);
  static Vector left(p+1);
  right.resize(p+1);
  left.resize(p+1);

  // double **ndu  = init2DArray(p + 1, p + 1);
  // double **a    = init2DArray(p + 1, p + 1);

  static Matrix ndu(p+1,p+1);
  static Matrix a(p+1,p+1);
  ndu.resize(p+1,p+1);
  a.resize(p+1,p+1);

  ndu(0,0) = 1.;
  for ( j = 1; j <= p; j++ )
  {
    left[j] = u - knot[i + 1 - j];
    right[j] = knot[i + j] - u;
    saved = 0.0;
    for ( r = 0; r < j; r++ )
    {
      ndu(j,r) = right[r + 1] + left[j - r];
      temp = ndu(r,j - 1) / ndu(j,r);

      ndu(r,j) = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    ndu(j,j) = saved;
  }
  for ( j = 0; j <= p; j++ )
    ders(0,j) = ndu(j,p);

  if ( order == 0 )
    return;


  for ( r = 0; r <= p; r++ )
  {
    int s1 = 0, s2 = 1;
    a(0,0) = 1.0;

    for ( k = 1; k <= order; k++ )
    {
      double d = 0.;
      int rk = r - k, pk = p - k;
      if ( r >= k )
      {
        a(s2,0) = a(s1,0) / ndu(pk + 1,rk);
        d = a(s2,0) * ndu(rk,pk);
      }
      j1 = rk >= -1 ? 1 : -rk;
      j2 = (r - 1 <= pk) ? k - 1 : p - r;
      for ( j = j1; j <= j2; j++ )
      {
        a(s2,j) = (a(s1,j) - a(s1,j - 1)) / ndu(pk + 1,rk + j);
        d += a(s2,j) * ndu(rk + j,pk);
      }
      if ( r <= pk )
      {
        a(s2,k) = -a(s1,k - 1) / ndu(pk + 1,r);
        d += a(s2,k) * ndu(r,pk);
      }
      ders(k,r) = d;
      j = s1; s1 = s2; s2 = j;
    }
  }
  r = p;
  for ( k = 1; k <= order; k++ )
  {
    for ( j = 0; j <= p; j++ )
      ders(k,j) *= r;
    r *= (p - k);
  }

  // free(left);
  // free(right);

  // free2Darray(ndu, p + 1);
  // free2Darray(a, p + 1);

}


double OneBasisFun(int p, int m, Vector U, int i, double u)
{
  /*
    Compute an individual B-spline function
  */

  double *N = (double*)malloc(sizeof(double) * (p + 1));

  if ((i == 0 && u == U[0] ) ||
      (i == (m - p - 1) && u == U[m]))
  {
    return (1.0);
  }

  if (u < U[i] || u >= U[i + p + 1])
  {
    return (0.0);
  }

  int j;
  for (j = 0; j <= p; j++)
  {
    if (u >= U[i + j] && u < U[i + j + 1]) N[j] = 1.0;
    else N[j] = 0.0;
  }

  int k;
  double saved, Uleft, Uright, temp;

  for (k = 1; k <= p; k++)
  {
    if (N[0] == 0.0) saved = 0.0;
    else saved = ((u - U[i]) * N[0]) / (U[i + k] - U[i]);
    for (j = 0; j < (p - k + 1); j++)
    {
      Uleft = U[i + j + 1];
      Uright = U[i + j + k + 1];
      if (N[j + 1] == 0.0)
      {
        N[j] = saved; saved = 0.0;
      }
      else
      {
        temp = N[j + 1] / (Uright - Uleft);
        N[j] = saved + (Uright - u) * temp;
        saved = (u - Uleft) * temp;
      }
    }
  }

  double Nip = N[0];

  free(N);

  return Nip;

}


void dersOneBasisFuns(int p, int m, Vector U, int i, double u, int order, double* ders)
{
  /*
    Compute the derivatives for basis function Nip
  */

  double **N = init2DArray(order + 1, order + 1);
  double *ND = (double*)malloc((order + 1) * sizeof(double));

  int k, j, jj;
  double Uleft, Uright, saved, temp;

  if (u < U[i] || u >= U[i + p + 1])
  {
    for (k = 0; k <= order; k++)
    {
      ders[k] = 0.0;
    }
    return;

  }

  for (j = 0; j <= p; j++)
  {
    if (u >= U[i + j] && u < U[i + j + 1])
      N[j][0] = 1.0;
    else
      N[j][0] = 0.0;
  }

  for (k = 1; k <= p; k++)
  {
    if (N[0][k - 1] == 0.0)
      saved = 0.0;
    else
      saved = ((u - U[i]) * N[0][k - 1]) / ( U[i + k] - U[i] );

    for (j = 0; j < (p - k + 1); j++)
    {
      Uleft = U[i + j + 1];
      Uright = U[i + j + k + 1];
      if (N[j + 1][k - 1] == 0.0)
      {
        N[j][k] = saved; saved = 0.0;
      }
      else
      {
        temp = N[j + 1][k - 1] / (Uright - Uleft);
        N[j][k] = saved + (Uright - u) * temp;
        saved = (u - Uleft) * temp;
      }
    }
  }

  ders[0] = N[0][p];

  for (k = 1; k <= order; k++)
  {
    for (j = 0; j <= k; j++)
      ND[j] = N[j][p - k];
    for (jj = 1; jj <= k; jj++)
    {
      if (ND[0] == 0.0)
        saved = 0.0;
      else
        saved = ND[0] / ( U[i + p - k + jj] - U[i]);
      for (j = 0; j < (k - jj + 1); j++)
      {
        Uleft = U[i + j + 1];
        Uright = U[i + j + p + jj];
        if (ND[j + 1] == 0.0)
        {
          ND[j] = (p - k + jj) * saved; saved = 0.0;
        }
        else
        {
          temp = ND[j + 1] / (Uright - Uleft);
          ND[j] = (p - k + jj) * (saved - temp);
          saved = temp;
        }
      }
    }
    ders[k] = ND[0];
  }

  free2Darray(N, order + 1);
  free(ND);

}








