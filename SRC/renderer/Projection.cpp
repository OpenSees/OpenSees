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
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Projection.cpp,v $
                                                                        
                                                                        
#include "Projection.h"
#include "View.h"

Projection::Projection()
{
  projection_mode = PERSPECTIVE_MODE;
}

Projection::~Projection()
{

}

int
Projection::update(void) {

  float midU, midV, dopU, dopV, dopN, shU, shV, F, B;
  float diffU, diffV;
  midU = (vpwindow[0] + vpwindow[1])/2;
  midV = (vpwindow[2] + vpwindow[3])/2;
  dopN = -cop[2];
  dopU = midU - cop[0];
  dopV = midV - cop[1];

  diffU = vpwindow[1]-vpwindow[0];
  diffV = vpwindow[3]-vpwindow[2];      
  
  shU = dopU/dopN;
  shV = dopV/dopN;
  F = planedist[0];
  B = planedist[2];
  if (F < 0) F = -F; // WARNING MESSAGES !!!
  if (B < 0) B = -B;
  
  if (projection_mode == PARALLEL_MODE) {
      
      if (F == -B)
	  B = F;

      TMat.Set(   2.0/diffU,        0.0,          0.0,    0.0,
	                0.0,   2.0/diffV,         0.0,    0.0,
	        2*shU/diffV, 2*shV/diffV,      1/(F+B),    0.0,
	        -2*midU/diffU, -2*midV/diffV, -F/(F+B),    1.0);

  } else { // perspective projection
      float a,b,c,e,f,g,h;
      float diffU, diffV;// prpN;

      diffU = vpwindow[1]-vpwindow[0];
      diffV = vpwindow[3]-vpwindow[2];      
      dopN = -dopN; // dopN'
      
      if (dopN > 0) dopN = -dopN; // WARNING MESSAGE// WARNING MESSAGE
      a = 2*dopN/(diffU * (dopN - B));
      b = 2*dopN/(diffV * (dopN - B));
      c = -1.0/(dopN - B);

      float zMin = (dopN + F)*c;
      if (zMin >= 0) {  // WARNING MESSAGE
	  F = -1.0e-8-dopN;
	  zMin = -1.0e8;
      }
      e = 1/(1+zMin);
      f = -zMin * e;
      
      g = cop[0] - shU * cop[2];
      h = cop[1] - shV * cop[2];
  
      TMat.Set(     a,   0.0,         0.0,    0.0,
	          0.0,     b,         0.0,    0.0,
	        a*shU, b*shV,         e*c,     -c,
                 -a*g,  -b*h, e*c*dopN+f, -c*dopN);  // dopN = -cop[2]!

  }
  return 0;
}

FACE &
Projection::transform(FACE &input)
{
  MYPOINT *point;
  
  // transform all the points by the transformation matrix
  // remember that in previous pipeline all points were unmarked
  FOR_EACH(point, input.pointList) {
    point->Transform(&TMat);
    if (projection_mode == PERSPECTIVE_MODE) {
	float wVal = point->p[3];
      if (wVal == 0) {
	opserr << "Some Point in local coord sys at same level as eye\n";
	return input;
      }
      else { // we will homo x and y coord, but leave z as was before we came in
	point->p[0] = point->p[0]/wVal;
	point->p[1] = point->p[1]/wVal;
	point->p[2] = point->p[2]/wVal;
	point->p[3] = 1.0;
      }
    }
  }
  
  return input;
}
  
  
MYPOINT *
Projection::transformP(MYPOINT *input)
{
    if (input == 0)
	return 0;
    
    input->Transform(&TMat);

    if (projection_mode == PERSPECTIVE_MODE) {
	float wVal = input->p[3];
	if (wVal == 0) { // point at eye position ignore delete it
	    delete input;
	    opserr << "HELP\n";
	    return 0;
	} else { // we will homo x and y coord, but leave z as was before we came in
	    input->p[0] = input->p[0]/wVal;
	    input->p[1] = input->p[1]/wVal;
	    input->p[2] = input->p[2]/wVal;
	    input->p[3] = 1.0;
	}	
    }
    return input;
}
