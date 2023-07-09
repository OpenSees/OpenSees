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
                                                                        
// $Revision: 1.6 $
// $Date: 2010-02-04 01:17:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthND.h,v $
                                                                        
// Written: MHS
// Created: Sept 2000
//
// Description: This file contains the class definition for ZeroLengthND.
// A ZeroLengthND element is defined by two nodes with the same coordinate.
// An NDMaterial object of order 2 or 3 is associated with the nodes to
// provide the basic force-deformation relationship for the element.
// If the NDMaterial is of order 2, an optional UniaxialMaterial may 
// be used to represent the (uncoupled) force-deformation relationship 
// orthogonal to the plane of the NDMaterial. The ZeroLengthND element
// only accounts for translational force-deformation relations.

#ifndef ZeroLengthND_h
#define ZeroLengthND_h

#include <Element.h>
#include <Matrix.h>

// Tolerance for zero length of element
#define	LENTOL 1.0e-6

class Node;
class Channel;
class NDMaterial;
class UniaxialMaterial;
class Response;

class ZeroLengthND : public Element
{
  public:
    
    // Constructor for a single Nd material model of order 2 or 3
    ZeroLengthND(int tag, 			      
		 int dimension,
		 int Nd1, int Nd2, 
		 const Vector& x,
		 const Vector& yprime,
		 NDMaterial& theNDMaterial);
    
    // Constructor for an Nd material model of order 2 with an 
    // uncoupled 1d material model acting out of plane
    ZeroLengthND(int tag, 			      
		 int dimension,
		 int Nd1, int Nd2, 
		 const Vector& x,
		 const Vector& yprime,
		 NDMaterial& theNDMaterial,
		 UniaxialMaterial &the1DMaterial);
    
    ZeroLengthND();    
    ~ZeroLengthND();

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

    Response *setResponse(const char **argv, int argc, OPS_Stream &theHandler);
    int getResponse(int responseID, Information &eleInformation);
    
  protected:
    
  private:
    // private methods
    void setUp (int Nd1, int Nd2, const Vector& x, const Vector& y);
    void setTransformation(void);
	void computeStrain(void);

    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes;         // contains the tags of the end nodes
    int dimension;                      // = 1, 2, or 3 dimensions
    int numDOF;	                        // number of dof for ZeroLengthND
    Matrix transformation;		// transformation matrix for orientation
	
    Matrix *A;	// Transformation matrix ... e = A*(u2-u1)
    Vector *v;	// NDMaterial strain vector, the element basic deformations
    double e;	// UniaxialMaterial strain

    Matrix *K;	// Element stiffness matrix
    Vector *P;	// Element force vector

    Node *end1Ptr;      		// pointer to the end1 node object
    Node *end2Ptr;      		// pointer to the end2 node object	

    NDMaterial *theNDMaterial;	// Pointer to NDMaterial object
    UniaxialMaterial *the1DMaterial;	// Pointer to UniaxialMaterial object
    int order;		// Order of the NDMaterial (2 or 3)
    
    // Class wide matrices for return
    static Matrix K6;
    static Matrix K12;
    
    // Class wide vectors for return
    static Vector P6;
    static Vector P12;
    
    // Class wide vectors for storing NDMaterial strains
    static Vector v2;
    static Vector v3;
    static Vector v5;
    static Vector v6;    
};

#endif
