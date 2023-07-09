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
                                                                       
                                                                       
#ifndef BeamContact3Dp_h
#define BeamContact3Dp_h

// Created: C.McGann, UW, 12.2011
//          adopted from Kathy Petek's BeamContact3D element
//
// Description: This file contains the class definition for BeamContact3Dp.

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

// number of nodes per element
#define BC3Dp_NUM_NODE 3
// number of dimesnions
#define BC3Dp_NUM_NDM  3
// degrees of freedom per element
#define BC3Dp_NUM_DOF  15

class Domain;
class Node;
class Channel;
class NDMaterial;
class ContactMaterial3D;
class FEM_ObjectBroker;
class CrdTransf;

class BeamContact3Dp : public Element
{
  public:
    BeamContact3Dp(int tag, int Nd1, int Nd2, int NdS, double rad, CrdTransf &coordTransf,
                   NDMaterial &theMat, double pen, int cSwitch = 0);
    BeamContact3Dp();
    ~BeamContact3Dp();

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
   
    // public methods to obtain stiffness, mass, damping and
    // residual information    
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

    // Response *setResponse(const char **argv, int argc, Information &eleInfo);
    Response *setResponse(const char **argv, int argc, OPS_Stream &eleInfo);
    int getResponse(int responseID, Information &eleInformation);

	// public methods for material stage update
	int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

  protected:
   
  private:

    // member functions
    double project(double xi);
    int UpdateBase(double xi);                             // method to update base vectors g1 & g2
    void ComputeB(void);                                   // method to compute Bn, Bs @ step n
    Matrix ComputeBphi(void);                              // method to compute Bphi, used in ComputeB and update
    void UpdateTransforms(void);                           // method to update Qa, Qb
    void ComputeQc(double xi);                             // method to compute Qc from Qa and Qb

    Matrix ExponentialMap(Vector theta);                   // function returns the exponential map of a vector
    Matrix ComputeSkew(Vector theta);                      // function returns skew matrix of given vector
    Vector CrossProduct(Vector &V1, Vector &V2);           // cross product (does not exist in Vector Class!)
    Matrix Transpose(int dim1,int dim2,const Matrix &M);   // functions returns the tranpose of Matrix M 
   
    Vector Geta1(void);                                    // returns a1 = mQa(:,0)      
    Vector Getb1(void);                                    // returns b1 = mQb(:,0)
    void Setc1(Vector c1_vec);                             // sets member vector c1
    Vector Getc1(void);                                    // returns member vector c1
    Vector Getdx_c(double xi);                             // returns dx_c / dxi
    Vector Getddx_c(double xi);                            // returns d^2(x_c)/ dxi^2
   
    // objects
    CrdTransf  *crdTransf;                                 // pointer to coordinate transformation object
    ContactMaterial3D *theMaterial;                        // contact material object
   
    ID  externalNodes;                                     // contains the tags of the end nodes
    Vector theVector;                                      // vector to return the residual
    Matrix mTangentStiffness;                              // Tangent Stiffness matrix
    Vector mInternalForces;                                // vector of Internal Forces
    Node *theNodes[BC3Dp_NUM_NODE];
       
    // input quantities
    double mL;                                             // length of Beam Element            
    double mRadius;                                        // radius of Pile Surface
	double mPenalty;                                       // penalty parameter for contact constraint enforcement
	int mIniContact;                                       // initial contact switch (0 = notInContact, 1 = inContact)
	                                                       // default is set for initially inContact
   
    // boolean variables
    bool inContact;
    bool was_inContact;
    bool in_bounds;    
	bool mInitialize;
   
    // calculation variables
    double mxi;                                            // centerline projection coordinate: 0 <= xi <= 1
    double mchi;                                           // twist rotation from end 1 to end 2
    double mGap;                                           // current value of the gap
    double mLambda;                                        // current value of Lagrange Multiplier
    double mrho2;                                          // angular coord relating c2 to radial vector
    double mrho3;                                          // angular coord relating c3 to radial vector
   
    Matrix meye1;                                          // Identity Tensor
    Vector mg1;                                            // tangent plane basis vector, g_xi
    Vector mg2;                                            // tangent plane basis vector, g_psi
    Matrix mg_metric;                                      // metric tensor
    Vector mn;                                             // normal Vector
    Vector mH;                                             // vector of Hermitian Basis Functions
    Vector mIcrd_a;                                        // initial coordinates of node a
    Vector mIcrd_b;                                        // initial coordinates of node b
    Vector mIcrd_s;                                        // initial coordinates of node s
    Vector mDcrd_a;                                        // current coordinates of node a
    Vector mDcrd_b;                                        // current coordinates of node b
    Vector mDcrd_s;                                        // current coordinates of node s
    Vector mDisp_a_n;                                      // total disps & rotations of node a @ step n
    Vector mDisp_b_n;                                      // total disps & rotations of node b @ step n
    Vector mDisp_s_n;                                      // total disps of node s @ step n
    Matrix mQa;                                            // coordinate transform for node a
    Matrix mQb;                                            // coordinate transform for node b
    Matrix mQc;
    Vector mc1;                                            // tangent vector at project point c
    Vector mBn;                                            // gap-displacement matrix
    Matrix mBs;                                            // slip-displacement matrix
    Matrix mBphi;
    Vector mSlip;
};

#endif
