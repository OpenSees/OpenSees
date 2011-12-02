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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/main.cpp,v $
                                                                        
                                                                        
// File: ~/model/main.C
//
// Written: fmk 12/95
// Revised:
//
// Purpose: This file is a driver to create a 2d plane-frame
// model and to perform a linear static analysis on that model.
// 
//

#include <stdlib.h>
#include <iostream.h>

#include <View.h>
#include <Viewport.h>
#include <Scan.h>
#include <Device.h>
#include <Projection.h>
#include <Clipping.h>
#include <db.H>

#include <FilePlotter.h>

int main(int argc, char **argv)
{

  Recorder *thePlotter = new FilePlotter("test","g3",10,10,300,300);
  thePlotter->record(0);

  /*
  theDevice.WINOPEN(600,300);
  theScan.update();

  FACE *theFace = new FACE();
  POINT *point = new POINT(1,10.0,10.0,0.0);
  point->r = 0.5;
  point->g = 0.5;
  point->b = 0.5;
  theFace->AddPoint(*point);

  point = new POINT(1,100.0,100.0,0.0);
  point->r = 0.3;
  point->g = 0.3;
  point->b = 0.3;
  theFace->AddPoint(*point);

  point = new POINT(1,200.0,50.0,0.0);
  point->r = 0.6;
  point->g = 0.6;
  point->b = 0.6;
  theFace->AddPoint(*point);
  */

  /***************
  FACE &theVFace = theView.transform(theFace);
  FACE &thePFace = theProjection.transform(theVFace);
  FACE &theCFace = theClipping.transform(thePFace);
  FACE &theVPFace = theViewport.transform(theCFace);
  theScan.transform(theVPFace);

  ******************/

  /*
  theScan.setFillMode(1);
  theScan.transform(*theFace);
  theDevice.ENDIMAGE();


  */


  int tag;
  cin >> tag;

  exit(0);
  return 0;
}

