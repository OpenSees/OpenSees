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
                                                                        
// $Revision: 1.19 $
// $Date: 2009-08-25 22:32:08 $
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

#include <DomainComponent.h>
#include <ID.h>

class Matrix;
class Vector;
class Info;
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
    virtual double getCharacteristicLength(void);

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
    virtual const Matrix &getGeometricTangentStiff();

    // methods for applying loads
    virtual void zeroLoad(void);	
    virtual int addLoad(ElementalLoad *theLoad, double loadFactor);
    virtual int addLoad(ElementalLoad *theLoad, const Vector &loadFactors);

    virtual int addInertiaLoadToUnbalance(const Vector &accel);
    virtual int setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc);

    // methods for obtaining resisting force (force includes elemental loads)
    virtual const Vector &getResistingForce(void) =0;
    virtual const Vector &getResistingForceIncInertia(void);        

    // method for obtaining information specific to an element
    virtual Response *setResponse(const char **argv, int argc, 
				  OPS_Stream &theHandler);
    virtual int getResponse(int responseID, Information &eleInformation);
    virtual int getResponseSensitivity(int responseID, int gradIndex,
				       Information &eleInformation);

    virtual int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual int addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool tag);
    virtual const Vector & getResistingForceSensitivity(int gradIndex);
    virtual const Matrix & getInitialStiffSensitivity(int gradIndex);
    virtual const Matrix & getDampSensitivity(int gradIndex);
    virtual const Matrix & getMassSensitivity(int gradIndex);
    virtual int   commitSensitivity(int gradIndex, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////

    virtual int addResistingForceToNodalReaction(int flag);

    virtual int storePreviousK(int numK);
    virtual const Matrix *getPreviousK(int num);

  protected:
    const Vector &getRayleighDampingForces(void);
    double alphaM, betaK, betaK0, betaKc;
    Matrix *Kc; // pointer to hold last committed matrix if needed for rayleigh damping

    Matrix **previousK;
    int numPreviousK;

    int index, nodeIndex;

    static Matrix ** theMatrices; 
    static Vector ** theVectors1; 
    static Vector ** theVectors2; 
    static int numMatrices;

  private:
};


#endif

