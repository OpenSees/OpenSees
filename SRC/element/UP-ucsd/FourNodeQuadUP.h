///////////////////////////////////////////////////////////////////////////////
// Description: This file contains the class declaration for FourNodeQuadUP. //
// FourNodeQuadUP is a 4-node plane strain element for solid-fluid fully     //
// coupled analysis. This implementation is a simplified u-p formulation     //
// of Biot theory (u - solid displacement, p - fluid pressure). Each element //
// node has two DOFs for u and 1 DOF for p.                                  //
// Written by Zhaohui Yang	(May 2002)                                   //
// based on FourNodeQuad element by Michael Scott                            //
///////////////////////////////////////////////////////////////////////////////

// $Revision: 1.5 $
// $Date: 2007-02-02 01:44:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/UP-ucsd/FourNodeQuadUP.h,v $

#ifndef FourNodeQuadUP_h
#define FourNodeQuadUP_h

#ifndef _bool_h
#include "bool.h"
#endif


#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

class FourNodeQuadUP : public Element
{
  public:
    FourNodeQuadUP(int tag, int nd1, int nd2, int nd3, int nd4,
		  NDMaterial &m, const char *type,
		  double t, double bulk, double rhof, double perm1, double perm2,
		   double b1 = 0.0, double b2 = 0.0, double p = 0.0);

    FourNodeQuadUP();
    virtual ~FourNodeQuadUP();
    const char *getClassType(void) const {return "FourNodeQuadUP";};
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

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);
    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
		  &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.

    friend class PyLiq1;
    friend class TzLiq1;

  protected:

  private:
    // private attributes - a copy for each object of the class
    NDMaterial **theMaterial; // pointer to the ND material objects
    ID connectedExternalNodes; // Tags of quad nodes
    Node *nd1Ptr;		// Pointers to quad nodes
    Node *nd2Ptr;
    Node *nd3Ptr;
    Node *nd4Ptr;

    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector
    Vector Q;		// Applied nodal loads
    double b[2];		// Body forces
	double appliedB[2]; // Body forces applied with load pattern, C.McGann, U.Washington
	int applyLoad;      // flag for body force in load, C.McGann, U.Washington
    Vector pressureLoad;	// Pressure load at nodes

    double thickness;	// Element thickness
    double rho;			// Fluid mass per unit volume
    double kc;   // combined bulk modulus
    double pressure;	// Normal surface traction (pressure) over entire element

    // Note: positive for outward normal
    double perm[2];  // lateral/vertical permeability
    static double shp[3][4][4];	// Stores shape functions and derivatives (overwritten)
    static double pts[4][2];	// Stores quadrature points
    static double wts[4];		// Stores quadrature weights
    static double dvol[4];  // Stores detJacobian (overwritten)
    static double shpBar[3][4]; // Stores averaged shap functions (overwritten)
    // private member functions - only objects of this class can call these
    double mixtureRho(int ipt);  // Mixture mass density at integration point i
    void shapeFunction(void);
    void setPressureLoadAtNodes(void);

    Matrix *Ki;
    static Node *theNodes[4];

    double *end1InitDisp;
    double *end2InitDisp;
    double *end3InitDisp;
    double *end4InitDisp;
};

#endif
