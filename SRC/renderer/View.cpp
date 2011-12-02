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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:01:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/View.cpp,v $
                                                                        
                                                                        
#include "View.h"
#include "db.H"

View theView;

View::View()
{

}

View::~View()
{

}

int
View::update(void)
{
  // first determine coords for Viewing system
  VECTOR u, v, n;
//  float scalar;
  int i;


  for (i=0; i<3; i++) {
    n[i] = vpn[i];
    v[i] = vuv[i];
  }
   
  if (n.Normalize() != 0) {
    opserr << "View::update() - VPN cannot have zero length\n";
    return -1;
  }

  u = v % n;
  if (u.Normalize() != 0) {
    opserr << "View::update() - VUV X VPN cannot have zero length\n";
    return -1;
  }

   v = n % u;
   v.Normalize();

   TMat.Set(       u[0],   v[0],   n[0], 0,
		   u[1],   v[1],   n[1], 0,
		   u[2],   v[2],   n[2], 0,
	          -vrp.Dot(u), -vrp.Dot(v), -vrp.Dot(n), 1.0);

   return 0;
}


FACE &
View::transform(FACE &input)
{
  // transform all the points by the transformation matrix
  // remember that in previos pipeline all points were marked
  MYPOINT *point;
  FOR_EACH(point, input.pointList) {
    point->Transform(&TMat);
  };

  return input;
}


MYPOINT *
View::transformP(MYPOINT *input)
{
    if (input != 0)
	input->Transform(&TMat);

    return input;
}



















