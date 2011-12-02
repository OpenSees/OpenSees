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
                                                                        
// $Revision: 1.11 $
// $Date: 2003-10-30 22:34:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Element.h,v $
                                                                        
                                                                        
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
class Response;
class ElementalLoad;
class Node;

class Element : public DomainComponent
{
  public:
    Element(int tag, int classTag);    
    virtual ~Element();

    // methods dealing with nodes and number of external dof
    virtual int getNumExternalNodes(void) const =0;
    virtual const ID &getExternalNodes(void)  =0;	
    virtual Node **getNodePtrs(void)  =0;	
    virtual int getNumDOF(void) =0;

    // methods dealing with committed state and update
    virtual int commitState(void);    
    virtual int revertToLastCommit(void) = 0;        
    virtual int revertToStart(void);                
    virtual int update(void);
    virtual bool isSubdomain(void);
    
    // methods to return the current linearized stiffness,
    // damping and mass matrices
    virtual const Matrix &getTangentStiff(void) =0;
    virtual const Matrix &getInitialStiff(void) =0;
    virtual const Matrix &getDamp(void);    
    virtual const Matrix &getMass(void);    

    // methods for applying loads
    virtual void zeroLoad(void) =0;	
    virtual int addLoad(ElementalLoad *theLoad, double loadFactor) =0;
    virtual int addInertiaLoadToUnbalance(const Vector &accel) =0;
    virtual int setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc);

    // methods for obtaining resisting force (force includes elemental loads)
    virtual const Vector &getResistingForce(void) =0;
    virtual const Vector &getResistingForceIncInertia(void);        

    // method for obtaining information specific to an element
    virtual Response *setResponse(const char **argv, int argc, Information &eleInformation);
    virtual int getResponse(int responseID, Information &eleInformation);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual int addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool tag);
    virtual int setParameter(const char **argv, int argc, Information &info);
    virtual int updateParameter(int parameterID, Information &info);
    virtual int activateParameter(int parameterID);
    virtual const Vector & getResistingForceSensitivity(int gradNumber);
    virtual const Matrix & getInitialStiffSensitivity(int gradNumber);
	virtual const Matrix & getDampSensitivity(int gradNumber);
    virtual const Matrix & getMassSensitivity(int gradNumber);
    virtual int   commitSensitivity(int gradNumber, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////

  protected:
    const Vector &getRayleighDampingForces(void);

    double alphaM, betaK, betaK0, betaKc;
    Matrix *Kc; // pointer to hold last committed matrix if needed for rayleigh damping

  private:
    int index;

    static Matrix ** theMatrices; 
    static Vector ** theVectors1; 
    static Vector ** theVectors2; 
    static int numMatrices;
};


#endif

