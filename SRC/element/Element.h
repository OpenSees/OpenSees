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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Element.h,v $
                                                                        
                                                                        
// File: ~/element/Element.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for Element.
// Element is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) Element.h, revA"

#ifndef Element_h
#define Element_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <ID.h>

#include <DomainComponent.h>

class Matrix;
class Vector;
class Renderer;
class Info;
class Information;

class Element : public DomainComponent
{
  public:
    Element(int tag, int classTag);    
    virtual ~Element();

    // methods dealing with nodes and number of external dof
    virtual int getNumExternalNodes(void) const =0;
    virtual const ID &getExternalNodes(void)  =0;	
    virtual int getNumDOF(void) =0;

    // methods dealing with committed state and update
    virtual int commitState(void) = 0;    
    virtual int revertToLastCommit(void) = 0;        
    virtual int revertToStart(void);                
    virtual int update(void);
    virtual bool isSubdomain(void);
    
    // methods to return the current linearized stiffness,
    // damping and mass matrices
    virtual const Matrix &getTangentStiff(void)=0;
    virtual const Matrix &getSecantStiff(void)=0;    
    virtual const Matrix &getDamp(void)=0;    
    virtual const Matrix &getMass(void)=0;    

    // methods for returning and applying loads
    virtual void zeroLoad(void) =0;	
    virtual int addLoad(const Vector &addP) =0;
    virtual int addInertiaLoadToUnbalance(const Vector &accel); 
    virtual const Vector &getResistingForce(void) =0;
    virtual const Vector &getResistingForceIncInertia(void) =0;        

    // method for obtaining information specific to an element
    virtual int setResponse(char **argv, int argc, Information &eleInformation);
    virtual int getResponse(int responseID, Information &eleInformation);
	
  protected:
    
  private:
};


#endif

