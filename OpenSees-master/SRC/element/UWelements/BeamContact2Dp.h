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
                                                                       
#ifndef BeamContact2Dp_h
#define BeamContact2Dp_h

// Created: C.McGann, UW, 12.2011
//
// Description: This file contains the class definition for BeamContact2Dp

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

// number of nodes per element
#define BC2Dp_NUM_NODE 3
// number of dimensions
#define BC2Dp_NUM_DIM  2
// degrees of freedom per element
#define BC2Dp_NUM_DOF  8

class Domain;
class Node;
class Channel;
class NDMaterial;
class ContactMaterial2D;
class FEM_ObjectBroker;

class BeamContact2Dp : public Element
{
  public:
    BeamContact2Dp(int tag, int Nd1, int Nd2, int NdS, NDMaterial &theMat, 
	               double width, double pen, int cSwitch = 0);
    BeamContact2Dp();
    ~BeamContact2Dp();

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
	double Project(double xi);       // method to determine centerline projection
	void ComputeB(void);             // method to compute Bn and Bs @ step n
	int UpdateBase(double xi);       // method to update base vector g_xi
	void UpdateEndFrames(void);      // method to update end node tangent vectors 
	Vector Get_dxc_xi(double xi);    // returns dx_c/dxi
	Vector Get_dxc_xixi(double xi);  // returns d^2(x_c)/dxi^2
	Vector Geta1(void);              // returns last converged a_1
	Vector Getb1(void);              // returns last converged b_1

    // objects
    ContactMaterial2D *theMaterial;  // contact material object
	ID mExternalNodes;               // contains the tags of the end nodes
	Matrix mTangentStiffness;        // tangent stiffness matrix
	Vector mInternalForces;          // vector of internal forces
	
	Node *theNodes[BC2Dp_NUM_NODE];

	// input quantities
	double mLength;                  // length of beam element
	double mRadius;                  // radius of beam (width/2)
	double mPenalty;                 // penalty parameter for contact constraint
	int mIniContact;                 // initial contact switch (0 = notInContact, 1 = inContact)
	                                 // default is set for initially inContact 

	// boolean variables
	bool inContact;
	bool was_inContact;
	bool in_bounds;

	// calculation variables
	double mXi;                      // centerline projection coordinate
	double mGap;                     // current value of the gap
	double mLambda;                  // current value of contact force (Lagrange Mult)

	Matrix mEye1;                    // identity tensor
    Matrix mEyeS;                    // skew symmetric transformation tensor
	Vector mg_xi;                    // surface tangent vector
	Vector mNormal;                  // normal vector
    Vector mShape;                   // vector of Hermitian shape functions
    Vector mDshape;                  // vector of Hermitian shape function derivatives

	Vector mBn;                      // gap-displacement vector
	Vector mBs;                      // slip-displacement vector
	
	Vector ma_1;                     // tangent vector at node a
	Vector mb_1;                     // tangent vector at node b
	Vector mc_1;                     // tangent vector at centerline projection
	double mrho;                     // local coordinate transformation term
	
	Vector mIcrd_a;                  // initial coordinates of node a
	Vector mIcrd_b;                  // initial coordinates of node a
	Vector mIcrd_s;                  // initial coordinates of secondary node
    Vector mDcrd_a;                  // initial coordinates of node a
	Vector mDcrd_b;                  // initial coordinates of node a
	Vector mDcrd_s;                  // initial coordinates of secondary node
	Vector mDisp_a_n;                // total disp & rotation of node a @ step n
	Vector mDisp_b_n;                // total disp & rotation of node b @ step n
};

#endif
