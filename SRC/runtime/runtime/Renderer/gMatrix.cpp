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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:01:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/gMatrix.cpp,v $
                                                                        
                                                                        
#include "gMatrix.H"

int
VECTOR::MakeNonHomogeneous()
// Goes from a homogeneous 4D Vector into a non-homogeneous 3D
// vector (4th coordinate is 1).
{
  if (v[3] == 0)
    {
      opserr << "ERROR: w coordinate of vector is 0, cannot unhomogeneize" << endln;
      return 0;
    }
  v[0] = v[0]/v[3];
  v[1] = v[1]/v[3];
  v[2] = v[2]/v[3];
  v[3] = 1;
  return 1;
}


VECTOR 
VECTOR::operator*(MATRIX &M)
// Implements VECTOR * MATRIX multiplications
{ 
  VECTOR *V = new VECTOR;
  V[0] = v[0] * M.m[0][0] + v[1] * M.m[1][0] + v[2] * M.m[2][0] + v[3] * M.m[3][0];
  V[1] = v[0] * M.m[0][1] + v[1] * M.m[1][1] + v[2] * M.m[2][1] + v[3] * M.m[3][1];
  V[2] = v[0] * M.m[0][2] + v[1] * M.m[1][2] + v[2] * M.m[2][2] + v[3] * M.m[3][2];
  V[3] = v[0] * M.m[0][3] + v[1] * M.m[1][3] + v[2] * M.m[2][3] + v[3] * M.m[3][3];
  return *V;
}
void
VECTOR::TimesMat(MATRIX &M)
// Implements MATRIX * VECTOR and overwrites the VECTOR
{
  float a,b,c,d;
  a = v[0] * M.m[0][0] + v[1] * M.m[0][1] + v[2] * M.m[0][2] + v[3] * M.m[0][3];
  b = v[0] * M.m[1][0] + v[1] * M.m[1][1] + v[2] * M.m[1][2] + v[3] * M.m[1][3];
  c = v[0] * M.m[2][0] + v[1] * M.m[2][1] + v[2] * M.m[2][2] + v[3] * M.m[2][3];
  d = v[0] * M.m[3][0] + v[1] * M.m[3][1] + v[2] * M.m[3][2] + v[3] * M.m[3][3];
  v[0] = a; v[1] = b; v[2] = c; v[3] = d;
}

void 
MATRIX::Set(float a1, float a2, float a3, float a4,
	    float b1, float b2, float b3, float b4,
	    float c1, float c2, float c3, float c4,
	    float d1, float d2, float d3, float d4)
// Sets matrix to specified values.
{
  m[0][0] = a1;
  m[0][1] = a2;
  m[0][2] = a3;
  m[0][3] = a4;
  
  m[1][0] = b1;
  m[1][1] = b2;
  m[1][2] = b3;
  m[1][3] = b4;
  
  m[2][0] = c1;
  m[2][1] = c2;
  m[2][2] = c3;
  m[2][3] = c4;
  
  m[3][0] = d1;
  m[3][1] = d2;
  m[3][2] = d3;
  m[3][3] = d4;
}

MATRIX 
MATRIX::Transpose()
{ 
  MATRIX *retMat = new MATRIX();
  int i,j;
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
       retMat->m[i][j] = m[j][i];
  return *retMat;
}

MATRIX 
MATRIX::operator*(MATRIX &M)
// Implements MATRIX * MATRIX multiplication
{
  MATRIX *retMat = new MATRIX();
  int i,j,k;
  float sum;
  for (i=0; i<4; i++) 
    for (j=0; j<4; j++) {
      sum = 0.0;
      for (k=0; k<4; k++) 
        sum = sum + m[i][k] * M.m[k][j];
      retMat->m[i][j] = sum;
    }
     
  return *retMat;
}


// Debugging output

OPS_Stream &
operator<<(OPS_Stream &os, VECTOR &v)
{
  os << "[ " << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << v[3] << " ]";
  return os;
}

OPS_Stream &
operator<<(OPS_Stream &os, MATRIX &M)
{
  os << "{ " << M.m[0][0] << ' ' << M.m[0][1] << ' ' << M.m[0][2] << ' ' << M.m[0][3] << " }" << endln;
  os << "{ " << M.m[1][0] << ' ' << M.m[1][1] << ' ' << M.m[1][2] << ' ' << M.m[1][3] << " }" << endln;
  os << "{ " << M.m[2][0] << ' ' << M.m[2][1] << ' ' << M.m[2][2] << ' ' << M.m[2][3] << " }" << endln;
  os << "{ " << M.m[3][0] << ' ' << M.m[3][1] << ' ' << M.m[3][2] << ' ' << M.m[3][3] << " }" << endln;
  return os;
}


