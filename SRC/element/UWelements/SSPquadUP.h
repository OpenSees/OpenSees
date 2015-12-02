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
                                                                       
#ifndef SSPquadUP_h
#define SSPquadUP_h

//
// Adapted: Luis Miranda, LNEC, Andre Barbosa, OSU; July 2015 - added boundary tractions (see LM changes)
//
// Created: C.McGann, UW, 05.2011
//
// Description: This file contains the class definition for SSPquadUP
//                Stabilized Single-Point Quad element with a u-p formulation 
//                for plane strain analysis of saturated porous media
//
// Reference:   Zienkiewicz, O.C. and Shiomi, T. (1984). "Dynamic behavior of 
//                saturated porous media; the generalized Biot formulation and 
//                its numerical solution." International Journal for Numerical 
//                Methods in Geomechanics, 8, 71-96.

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

// number of nodes per element
#define SQUP_NUM_NODE 4
// number of dimensions
#define SQUP_NUM_DIM  2
// degrees of freedom per element
#define SQUP_NUM_DOF  12

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;
class Response;

class SSPquadUP : public Element
{
  public:
    // LM change
    SSPquadUP(int tag, int Nd1, int Nd2, int Nd3, int Nd4, NDMaterial &theMat,
                       double thick, double Kf, double Rf, double k1, double k2,
                       double eVoid, double alpha, double b1 = 0.0, double b2 = 0.0,
					   double Pup = 0.0, double Plow = 0.0, double Pleft = 0.0, double Pright = 0.0);
    SSPquadUP();
    ~SSPquadUP();
    // LM change

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
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &eleInfo);
    int getResponse(int responseID, Information &eleInformation);

    // public methods for material stage update
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    // allow PyLiq1 and TzLiq1 classes to get stresses from SSPquadUP class
    friend class PyLiq1;
    friend class TzLiq1;

  protected:

  private:

    // member functions
    Matrix DyadicProd(Vector v1, Vector v2);   // dyadic product for two 2x1 vectors
    void GetStab(void);                        // compute stabilization stiffness matrix
    void GetSolidStiffness(void);              // compute solid phase stiffness matrix
    void GetSolidMass(void);                   // compute solid phase mass matrix
    void GetPermeabilityMatrix(void);          // compute permeability matrix
    // LM change        
	void setPressureLoadAtNodes(void);
    // LM change
		
    // objects
    NDMaterial *theMaterial;                   // pointer to NDMaterial object
    ID mExternalNodes;                         // contains tags of the nodes
    Matrix mTangentStiffness;                  // tangent stiffness matrix
    Vector mInternalForces;                    // vector of internal forces
    Vector Q;                                  // vector of applied nodal forces
    Matrix mDamp;                              // damping matrix
    Matrix mMass;                              // mass matrix
    // LM change
	Vector pressureLoad;                       // Pressure load at nodes 
    // LM change
		
    Node *theNodes[4];

    // input quantities
    double mThickness;                         // thickness of element
    double fBulk;                              // bulk modulus of pore fluid
    double fDens;                              // pore fluid mass density
    double eVoid;                              // voids ratio for solid phase
    double perm[2];                            // horiz and vertical permeability
    double b[2];                               // body forces acting on element
    // LM change
	double P[4];                               // Normal surface traction (pressure) positive for outward normal in each face
    double pressureUpperSide;
	double pressureLowerSide;
	double pressureLeftSide;
	double pressureRightSide;
	// LM change	
        
    // load pattern variables
    double appliedB[2];                        // body forces applied with load pattern
    int    applyLoad;                          // flag for body force in load pattern

    // calculation variables
    double J0;                                 // constant portion of jacobian
    double J1;                                 // linear (xi) portion of jacobian
    double J2;                                 // linear (eta) portion of jacobian
    double mPorosity;                          // porosity of solid phase n = e/(1+e)
    double mAlpha;
        
    Matrix Mmem;                               // mapping matrix for membrane modes
    Matrix Kstab;                              // stabilization stiffness matrix
    Matrix mNodeCrd;                           // nodal coordinate array
    Matrix dN;                                 // array of shape function derivatives

    Matrix mSolidK;                            // stiffness matrix for solid phase
    Matrix mSolidM;                            // mass matrix for solid phase
    Matrix mPerm;                              // permeability matrix H
};

#endif
 
