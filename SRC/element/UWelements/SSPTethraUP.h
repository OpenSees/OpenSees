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

#ifndef SSPTethraUP_h
#define SSPTethraUP_h

// Created: derived from SSPbrickUP (C.McGann, UW, 10.2011) - 2026
//
// Description: This file contains the class definition for SSPTethraUP
//                Single-Point Tetrahedron element (4 nodes, 1 Gauss point
//                at the element centroid) with a u-p formulation for 3D
//                analysis of saturated porous media.
//
//                The element is derived from the SSPbrickUP element of
//                C.McGann, P.Arduino, P.Mackenzie-Helnwein (U.Washington).
//                A linear (4-node) tetrahedron with one-point integration
//                has constant strains and no hourglass (zero-energy) modes,
//                therefore the membrane/hourglass stabilization of the
//                SSPbrick hexahedron is not required here. The pressure
//                field stabilization term (parameter alpha) is retained,
//                consistent with the SSP UP element family.
//
//                Each node has 4 dof: u_x, u_y, u_z, p (pore pressure).
//                The material model is the effective-stress skeleton
//                response (e.g. PressureDependMultiYield03).
//
// Reference:   Zienkiewicz, O.C. and Shiomi, T. (1984). "Dynamic behavior of
//                saturated porous media; the generalized Biot formulation and
//                its numerical solution." International Journal for Numerical
//                and Analytical Methods in Geomechanics, 8, 71-96.

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

// number of nodes per element
#define STUP_NUM_NODE 4
// number of dimensions
#define STUP_NUM_DIM  3
// degrees of freedom per element (4 dof per node: u_x, u_y, u_z, p)
#define STUP_NUM_DOF  16

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;
class Response;
class Parameter;

class SSPTethraUP : public Element
{
  public:
    SSPTethraUP(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
                NDMaterial &theMat, double Kf, double Rf,
                double k1, double k2, double k3,
                double eVoid, double alpha,
                double b1 = 0.0, double b2 = 0.0, double b3 = 0.0, int massType = 0);
    SSPTethraUP();
    ~SSPTethraUP();

    const char* getClassType() const { return "SSPTethraUP"; };

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
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &eleInfo);
    int getResponse(int responseID, Information &eleInformation);

    // public methods for material stage update
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

  protected:

  private:

    // member functions
    void setTetrahedronData(void);           // volume, shape function gradients, B matrix
    void GetSolidStiffness(void);            // solid phase stiffness matrix
    void GetSolidMass(void);                 // solid phase mass matrix
    void GetPermeabilityMatrix(void);        // permeability matrix H and pressure stabilization

    // objects
    NDMaterial *theMaterial;                 // pointer to NDMaterial object
    ID mExternalNodes;                       // contains tags of the nodes
    Matrix mTangentStiffness;                // tangent stiffness matrix
    Vector mInternalForces;                  // vector of internal forces
    Vector Q;                                // vector of applied nodal forces
    Matrix mDamp;                            // damping matrix
    Matrix mMass;                            // mass matrix

    Node *theNodes[4];

    // input quantities
    double fBulk;                            // bulk modulus of pore fluid
    double fDens;                            // mass density of pore fluid
    double mAlpha;                           // pressure field stabilization parameter
    double perm[3];                          // cartesian permeabilities
    double b[3];                             // body forces acting on element

    // load pattern variables
    double appliedB[3];                      // body forces applied with load pattern
    int    applyLoad;                        // flag for body force in load pattern

    // calculation variables
    double mVol;                             // element volume
    double mPorosity;                        // porosity of solid phase n = e/(1+e)
    int massType;                            // mass matrix type (0 = consistent, 1 = lumped)

    Matrix Bnot;                             // strain-displacement matrix (6 x 12)
    Matrix mSolidK;                          // stiffness matrix for solid phase (12 x 12)
    Matrix mSolidM;                          // mass matrix for solid phase (12 x 12)
    Matrix mPerm;                            // permeability matrix H (4 x 4)
    Matrix mPressStab;                       // pressure stabilization matrix (4 x 4)
    Matrix dNmod;                            // shape function gradients at gauss point (4 x 3)
};

#endif
