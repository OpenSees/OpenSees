/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
** ****************************************************************** */

#ifndef SSPTethraUP_h
#define SSPTethraUP_h

// Created: derived from SSPbrickUP (C.McGann, UW, 10.2011)
//
// SSPTethraUP: 4-node single-point tetrahedron u-p element.
//
// This version adds numerical regularization of the pore-pressure field
// to improve stability in undrained / nearly incompressible conditions:
//
//   1. Effective fluid bulk modulus is capped to avoid exactly zero
//      fluid compressibility.
//   2. A consistent fluid compressibility matrix is used, plus a small
//      artificial compressibility floor.
//   3. If input alpha <= 0, an automatic pressure stabilization
//      parameter is computed from element size and reference stiffness.
//   4. Negative tetrahedron Jacobians are handled using |detJ|.
//
// Each node has 4 dof: ux, uy, uz, p.
// The physical pore pressure is the VELOCITY of dof 4, consistent with
// the SSPbrickUP u-p implementation.

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#define STUP_NUM_NODE 4
#define STUP_NUM_DIM  3
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
                double b1 = 0.0, double b2 = 0.0, double b3 = 0.0,
                int massType = 0);

    SSPTethraUP();
    ~SSPTethraUP();

    const char* getClassType() const { return "SSPTethraUP"; };

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);

    void zeroLoad(void);
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact,
                    const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &eleInfo);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

protected:

private:
    void setTetrahedronData(void);
    void computeRegularization(void);
    double getReferenceModulus(void);

    void GetSolidStiffness(void);
    void GetSolidMass(void);
    void GetPermeabilityMatrix(void);

    NDMaterial *theMaterial;
    ID mExternalNodes;

    Matrix mTangentStiffness;
    Vector mInternalForces;
    Vector Q;
    Matrix mDamp;
    Matrix mMass;

    Node *theNodes[4];

    double fBulk;
    double fDens;

    double mAlpha;
    double perm[3];
    double b[3];

    double appliedB[3];
    int applyLoad;

    double mVol;
    double mPorosity;
    int massType;

    // Regularization quantities
    double KfEff;
    double alphaEff;
    double hChar;

    Matrix mFluidS;

    Matrix Bnot;
    Matrix mSolidK;
    Matrix mSolidM;
    Matrix mPerm;
    Matrix mPressStab;
    Matrix dNmod;
};

#endif