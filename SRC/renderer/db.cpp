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
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/db.cpp,v $
                                                                        
                                                                        
#include <string.h>

#include "db.H"

void
MYPOINT::Transform(MATRIX *M)
{
  float a,bb,c,d;
  a = p[0] * M->m[0][0] + p[1] * M->m[1][0] + p[2] * M->m[2][0] + 
    p[3] * M->m[3][0];
  bb = p[0] * M->m[0][1] + p[1] * M->m[1][1] + p[2] * M->m[2][1] + 
    p[3] * M->m[3][1];
  c = p[0] * M->m[0][2] + p[1] * M->m[1][2] + p[2] * M->m[2][2] + 
    p[3] * M->m[3][2];
  d = p[0] * M->m[0][3] + p[1] * M->m[1][3] + p[2] * M->m[2][3] + 
    p[3] * M->m[3][3];
  p[0] = a; p[1] = bb; p[2] = c; p[3] = d;
}


OPS_Stream &operator<<(OPS_Stream &os, MYPOINT &point)
{
  os <<"Point  ("<<point.X()<<' '<<point.Y()<<' '<<
    point.Z()<<')';
  os << " (" << point.r << ' ' << point.g << ' ' << point.b << " )";
  return os;
}
    
OPS_Stream &operator<<(OPS_Stream &os, FACE &face)
{
  os <<"Face: "<< endln;
  MYPOINT *point;
  FOR_EACH(point, face.pointList)
    {
      os << (*point) << endln;
    }
  return os;
}

