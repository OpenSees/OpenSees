///////////////////////////////////////////////////////////////////////////////
// Description: This file contains the class declaration for                 //
// NineFourNodeQuadUP, a 9-4-node (9 node for solid and 4 node for fluid) //
// plane strain element for solid-fluid fully coupled analysis. This         //
// implementation is a simplified u-p formulation of Biot theory             //
// (u - solid displacement, p - fluid pressure). Each element node has two   //
// DOFs for u and 1 DOF for p.                                               //
//                                                                           //
// Written by Zhaohui Yang	(March 2004)                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// $Revision: 1.5 $
// $Date: 2007-02-02 01:44:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/UP-ucsd/Nine_Four_Node_QuadUP.h,v $

#ifndef NineFourNodeQuadUP_h
#define NineFourNodeQuadUP_h

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

class NineFourNodeQuadUP : public Element
{
  public:
    NineFourNodeQuadUP(int tag, int nd1, int nd2, int nd3, int nd4,
		  int nd5, int nd6, int nd7, int nd8, int nd9,
		  NDMaterial &m, const char *type,
		  double t, double bulk, double rhof, double perm1, double perm2,
		   double b1 = 0.0, double b2 = 0.0);

    NineFourNodeQuadUP();
    virtual ~NineFourNodeQuadUP();
    const char *getClassType(void) const {return "NineFourNodeQuadUP";};
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

    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector
    Vector Q;		// Applied nodal loads
    double b[2];		// Body forces
    double appliedB[2]; // Body forces applied with load pattern, C.McGann, U.Washington

    int applyLoad;      // flag for body force in load, C.McGann, U.Washington
    Matrix *Ki;

    Node *theNodes[9];
    double thickness;	// Element thickness
    double rho;			// Fluid mass per unit volume
    double kc;   // combined bulk modulus
    double perm[2];  // lateral/vertical permeability
    static const int nintu;
    static const int nintp;
    static const int nenu;
    static const int nenp;

    static double shgu[3][9][9];	// Stores shape functions and derivatives (overwritten)
    static double shgp[3][4][4];	// Stores shape functions and derivatives (overwritten)
    static double shgq[3][9][4];	// Stores shape functions and derivatives (overwritten)
    static double shlu[3][9][9];	// Stores shape functions and derivatives
    static double shlp[3][4][4];	// Stores shape functions and derivatives
    static double shlq[3][9][4];	// Stores shape functions and derivatives
    static double wu[9];		// Stores quadrature weights
    static double wp[4];		// Stores quadrature weights
    static double dvolu[9];  // Stores detJacobian (overwritten)
    static double dvolp[4];  // Stores detJacobian (overwritten)
    static double dvolq[4];  // Stores detJacobian (overwritten)

    // private member functions - only objects of this class can call these
    double mixtureRho(int ipt);  // Mixture mass density at integration point i
    void shapeFunction(double *w, int nint, int nen, int mode);
    void globalShapeFunction(double *dvol, double *w, int nint, int nen, int mode);

    double *initNodeDispl;
};

#endif
