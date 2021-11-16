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
                                                                       
#ifndef SSPbrick_h
#define SSPbrick_h

// Created: C.McGann, UW, 10.2011
//
// Description: This file contains the class definition for SSPbrick (Stabilized Single-Point Brick)

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

// number of nodes per element
#define SSPB_NUM_NODE 8
// number of dimensions
#define SSPB_NUM_DIM  3
// degrees of freedom per element
#define SSPB_NUM_DOF  24

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;
class Response;

class SSPbrick : public Element
{
  public:
    SSPbrick(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8,
                      NDMaterial &theMat, double b1 = 0.0, double b2 = 0.0, double b3 = 0.0);
    SSPbrick();
    ~SSPbrick();

	const char* getClassType()  const { return "SSPbrick"; };

    // public methods to obtain information about dof and connectivity
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

	// public methods to obtain stiffness, mass, damping, and residual info
	const Matrix &getTangentStiff(void);
	const Matrix &getInitialStiff(void);
	const Matrix &getMass(void);

	void zeroLoad(void);
	int addLoad(ElementalLoad *theLoad, double loadFactor);
	int addInertiaLoadToUnbalance(const Vector &accel);
	const Vector &getResistingForce(void);
	const Vector &getResistingForceIncInertia(void);

	// public methods for element output
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
	void Print(OPS_Stream &s, int flag =0);

	Response *setResponse(const char **argv, int argc, OPS_Stream &eleInfo);
	int getResponse(int responseID, Information &eleInformation);

	// public methods for material stage update
	int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

  protected:

  private:

    // member functions
	void GetStab(void);                                 // compute stabilization stiffness matrix
	Vector CrossProduct(Vector v1, Vector v2);          // cross product for two 3x1 vectors
	Matrix Transpose(int d1, int d2, const Matrix &M);  // transpose operation

	// objects
	NDMaterial *theMaterial;                            // pointer to NDMaterial object
	ID mExternalNodes;                                  // contains tags of the nodes
	Matrix mTangentStiffness;                           // tangent stiffness matrix
	Vector mInternalForces;                             // vector of internal forces
	Vector Q;                                           // vector of applied nodal forces
	Matrix mMass;                                       // mass matrix

	Node *theNodes[8];

	// input quantities
	double b[3];                                        // body forces acting on element

	// load pattern variables
	double appliedB[3];                                 // body forces applied with load pattern
	int applyLoad;                                      // flag for body force in load pattern

	// calculation variables
	double J[20];                                       // jacobian determinant terms
	double mVol;                                        // element volume

	bool mInitialize;
	
	Matrix Bnot;                                        // mapping matrix for membrane modes
	Matrix Kstab;                                       // stabilization stiffness matrix
	Matrix mNodeCrd;                                    // nodal coordinate array
	
	Vector xi;                                          // xi evaluated at the nodes
	Vector et;                                          // eta evaluated at the nodes
	Vector ze;                                          // zeta evaluated at the nodes
	Vector hut;                                         // zeta*eta evaluated at the nodes
	Vector hus;                                         // zeta*xi evaluated at the nodes
	Vector hst;                                         // xi*eta evaluated at the nodes
	Vector hstu;                                        // xi*eta*zeta evaluated at the nodes
};

#endif
