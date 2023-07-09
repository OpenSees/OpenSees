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
                                                                        
#ifndef CoupledZeroLength_h
#define CoupledZeroLength_h

// Written: fmk
// Created: 04/11
// Modified from ZeroLength element of GLF
//
// Description: This file contains the class definition for CoupledZeroLength.
// A CoupledZeroLength element is defined by two nodes with the same coordinate.
// One material and 2 directions are provided. The element computes strain
// for material based on displacements in the 2 directions.

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

// Tolerance for zero length of element
#define	LENTOL 1.0e-6

// Type of dimension of element NxDy has dimension x=1,2,3 and
// y=2,4,6,12 degrees-of-freedom for the element
enum Etype { D1N2, D2N4, D2N6, D3N6, D3N12 };

class Node;
class UniaxialMaterial;

class CoupledZeroLength : public Element {
 public:
  
    // Constructor for a single 1d material model
    CoupledZeroLength(int tag, 			      
		      int Nd1, int Nd2, 
		      UniaxialMaterial& theMaterial,
		      int direction1, int direction2,
		      int doRayleighDamping = 0);
    
    CoupledZeroLength();    
    ~CoupledZeroLength();
    
    const char *getClassType(void) const {return "CoupledZeroLength";};
    
    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    
    int getNumDOF(void);	
    void setDomain(Domain *theDomain);
    
    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
    
    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);
    
    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);    
    
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    
    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);    
    void Print(OPS_Stream &s, int flag =0);    
    
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);
    
    int setParameter(const char **argv, int argc, Parameter &param);
    
// AddingSensitivity:BEGIN //////////////////////////////////////////
    const Vector &getResistingForceSensitivity(int gradIndex);
    int commitSensitivity(int gradIndex, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////

 protected:
    
 private:
    Etype elemType;
  
    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes;         // contains the tags of the end nodes
    int dimension;                      // = 1, 2, or 3 dimensions
    int numDOF;	                        // number of dof for CoupledZeroLength
    Matrix transformation;		// transformation matrix for orientation
    int useRayleighDamping;
    
    Node *theNodes[2];
    
    Matrix *theMatrix; 	    	// pointer to objects matrix (a class Matrix)
    Vector *theVector;      	// pointer to objects vector (a class Vector)
    
    // Storage for uniaxial material models
    UniaxialMaterial *theMaterial;    // array of pointers to 1d materials
    int dirn1;
    int dirn2;
    double dX;
    double dY;
    double fX;
    double fY;
    
    // vector pointers to initial disp and vel if present
    Vector *d0;
    Vector *v0;
};

#endif
