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
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Viewport.cpp,v $
                                                                        
                                                                        
#include "Viewport.h"
#include "Device.h"
#include "Projection.h"

Viewport::Viewport() 
{

}


Viewport::~Viewport() 
{

}

int
Viewport::update(void) {
  float Sx,Sy,Tx, Ty;
//  float VPright, VPleft, VPtop, VPbottom;
  float	PRTleft, PRTright, PRTtop, PRTbottom;

  float length, width;
  length = theDevice->GetHeight(); 
  width = theDevice->GetWidth();;

  PRTleft    = portwindow[0];    
  PRTright   = portwindow[1];
  PRTbottom  = portwindow[2];
  PRTtop     = portwindow[3];

  if (PRTleft < -1 || PRTright >1 || PRTtop >1 || PRTbottom <-1) {
    opserr << "Viewport::update() -  PORTWINDOW must be in range { -1 1 -1 1}\n";
    return -1;
  }
  
  Sx = (PRTright - PRTleft)/2 * width/2;
  Sy = (PRTtop - PRTbottom)/2 * length/2;
  Tx = Sx + (1 + PRTleft)*width/2;
  Ty = Sy + (1 + PRTbottom)*length/2;

  TMat.Set(   Sx, 0.0, 0.0, 0.0,
 	     0.0,  Sy, 0.0, 0.0,
	     0.0, 0.0, 1.0, 0.0,
	      Tx,  Ty, 0.0, 1.0);

  return 0;
}


FACE &
Viewport::transform(FACE &input)
{

  // transform all the points by the transformation matrix
  // remember that in previous pipeline all points were marked
  MYPOINT *point;
  FOR_EACH(point, input.pointList) {
    point->Transform(&TMat);
  };
  
  return input;
}



MYPOINT *
Viewport::transformP(MYPOINT *input)
{
    if (input != 0)
	input->Transform(&TMat);

    return input;
}


void
Viewport::setDevice(Device &device)
{
    theDevice = &device;
}

