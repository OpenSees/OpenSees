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
                                                                       
#ifndef BeamEndContact3Dp_h
#define BeamEndContact3Dp_h

// Created: C.McGann, UW, 12.2011
//
// Description: This file contains the class definition for BeamEndContact3Dp

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

// number of nodes per element
#define BEC3p_NUM_NODE 2
// number of dimensions
#define BEC3p_NUM_DIM  3
// degrees of freedom per element
#define BEC3p_NUM_DOF  9

class Domain;
class Node;
class Channel;
class FEM_ObjectBroker;

class BeamEndContact3Dp : public Element
{
  public:
    BeamEndContact3Dp(int tag, int Nd1, int Nd2, int NdS, double rad, double pen, int cSwitch = 0);
    BeamEndContact3Dp();
    ~BeamEndContact3Dp();

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

	void zeroLoad(void);
	int addLoad(ElementalLoad *theLoad, double loadFactor);
	int addInertiaLoadToUnbalance(const Vector &accel);
	const Vector &getResistingForce(void);
	const Vector &getResistingForceIncInertia(void);

	// public methods for element output
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel,
	             FEM_ObjectBroker &theBroker);
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
	Matrix GetSkew(Vector theta);    // function to return skew matrix of vector
	Matrix ExpMap(Vector theta);     // function computes exp map of vector

    // objects
	ID mExternalNodes;               // contains the tags of the end nodes
	Matrix mTangentStiffness;        // tangent stiffness matrix
	Vector mInternalForces;          // vector of internal forces
	
	Node *theNodes[BEC3p_NUM_NODE];

	// input quantities
	int mBeamNode;                   // contains the tag of the second beam node
	Node *x_b;                       // second beam node pointer
	double mRadius;                  // radius of beam element
	double mPenalty;                 // penalty parameter for contact constraint
	int mIniContact;                 // initial contact flag (0 = noContact,1 = inContact)
	                                 // default is set for initially inContact 

	// boolean variables
	bool inContact;
	bool was_inContact;
	bool in_bounds;

	// calculation variables
	double mGap;                     // current value of the gap
	double mLambda;                  // current value of contact force (Lagrange Mult)

	Matrix mEye1;                    // identity tensor
	Vector mNormal;                  // normal vector
	Vector mIniNormal;               // initial normal vector
	Vector mx_p;                     // projection coordinates

	Vector mIcrd_a;                  // initial coordinates of node a
	Vector mIcrd_b;                  // initial coordinates of node b
	Vector mIcrd_s;                  // initial coordinates of secondary node
    Vector mDcrd_a;                  // initial coordinates of node a
	Vector mDcrd_s;                  // initial coordinates of secondary node
};

#endif
