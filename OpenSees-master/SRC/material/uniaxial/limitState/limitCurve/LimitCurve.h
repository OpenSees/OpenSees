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
// $Date: 2006-09-05 22:32:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/limitCurve/LimitCurve.h,v $

//
// Written: kje
// Created: 08/01
// Revision: A
//
// Description: This file contains the class definition for 
// LimitCurve. LimitCurve is a base class and 
// thus no objects of it's type can be instantiated. It has pure virtual 
// functions which must be implemented in it's derived classes.  Sub classes
// will define the curve used by LimitStateMaterial to determine if a 
// limit state has been reached.


#ifndef LimitCurve_h
#define LimitCurve_h

#include <DomainComponent.h>
#include <MovableObject.h>

class ID;
class Vector;
class Matrix;
class Information;


class LimitCurve : public TaggedObject, public MovableObject
{
  public:
    LimitCurve (int tag, int classTag);
    virtual ~LimitCurve();

    virtual LimitCurve *getCopy (void) = 0;

    virtual int checkElementState(double springForce) = 0;
    
    virtual double getDegSlope(void) = 0;
    virtual double getResForce(void) = 0;
    virtual double getUnbalanceForce(void) = 0;
    
    virtual double findLimit(double input) = 0;

    virtual int revertToStart (void) = 0;        

  protected:
    
  private:
 

};

extern bool OPS_addLimitCurve(LimitCurve *newComponent);
extern LimitCurve *OPS_getLimitCurve(int tag);
extern void OPS_clearAllLimitCurve(void);

#endif

