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
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/PlaneFrame.h,v $
                                                                        
                                                                        
// File: ~/modelbuilder/PlaneFrame.h
// 
// Written: fmk 
// Created: Mon Sept 15 14:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for PlaneFrame.
// A PlaneFrame creates the plane frame models.

//
// What: "@(#) ModelBuilder.h, revA"

#ifndef PlaneFrame_h
#define PlaneFrame_h

#include <ModelBuilder.h>
#include <iostream.h>

class Element;
class Node;
class SP_Constraint;
class MP_Constraint;
class LoadCase;

class PlaneFrame : public ModelBuilder
{
  public:
    PlaneFrame(Domain &theDomain);
    ~PlaneFrame();    

    int buildFE_Model(void);

  protected:

  private:
};

#endif
