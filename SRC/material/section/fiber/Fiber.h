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
// $Date: 2000-12-18 10:34:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/Fiber.h,v $
                                                                        
                                                                        
// File: ~/fiber/Fiber.h
//
// Written: Remo Magalhaes de Souza
// Created: 10/98
// Revision: 
//
// Description: This file contains the class definition for 
// Fiber. Fiber is an abstract base class and thus no objects of 
// it's type can be instatiated. It has pure virtual functions which
// must be implemented in it's derived classes.
//
// What: "@(#) Fiber.h, revA"


#ifndef Fiber_h
#define Fiber_h

#include <DomainComponent.h>
#include <MovableObject.h>

class Vector;
class Matrix;
class ID;
class Information;
class Response;

class Fiber : public TaggedObject, public MovableObject
{
  public:
    Fiber (int tag, int classTag);
    virtual ~Fiber();

    virtual int    setTrialFiberStrain(const Vector &vs)=0;
    virtual Vector &getFiberStressResultants(void) =0;
    virtual Matrix &getFiberTangentStiffContr(void) =0;
	
    virtual int    commitState(void)=0;
    virtual int    revertToLastCommit(void)=0;    
    virtual int    revertToStart(void)=0;
    
	virtual Fiber *getCopy(void) = 0;
	virtual int getOrder(void) = 0;
	virtual const ID &getType(void) = 0;

	virtual Response *setResponse(char **argv, int argc, Information &info);
	virtual int getResponse(int responseID, Information &info);

	virtual void getFiberLocation(double &y, double &z) = 0;

  protected:
    
  private:
};


#endif

