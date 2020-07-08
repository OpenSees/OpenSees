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
                                                                        
// Purpose: This file contains the class definition for BeamColumnwLHNMYS.
// BeamColumnwLHNMYS is a plane frame member.

#ifndef BeamColumnwLHNMYS_h
#define BeamColumnwLHNMYS_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include "GPYS2d.h"

class Channel;
class Information;
class CrdTransf;
class SectionForceDeformation;
class Response;
class Renderer;

class BeamColumnwLHNMYS : public Element
{
  public:
    BeamColumnwLHNMYS();
    BeamColumnwLHNMYS(int tag, int Nd1, int Nd2, double a, double e,
                      double i, double npi, double npj, double mpi, double mpj,
                      CrdTransf &coordTransf, double ytol,
                      double wtol, double maxiter, Vector &hir, Vector &hkr,
                      double nr, Matrix &gpysc, double r, int cm);

    ~BeamColumnwLHNMYS();
    
    void* OPS_BeamColumnwLHNMYS(void);

    const char *getClassType(void) const {return "BeamColumnwLHNMYS";};

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
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);    
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes = 0, int numModes = 0);

    Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    int getResponse (int responseID, Information &info);
 
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);

  private:
    double A, E, I;             // area, elastic modulus, moment of inertia
    double NpI, NpJ, MpI, MpJ;  // plastic axial and moment capacities of element at end i ane end j
    double Wtol;                // incremental work tolerance for state convergence (default: 1e-16)
    double yftol;               // tolerance for yield criterion (default: 1e-8)
    double rho;                 // mass per unit length
    int cMass;                  // consistent mass flag
    int MaxIter;                // max. no. of iterations for state convergence (default: 15)
    int nrow;                   // no. of rows in matrix GPYSC
    Matrix GPYSC;               // Matrix with coefficients of polynomial yield surface
    Vector Hir;                 // isotropic hardening ratio for flexural end i and end j (default [0;0])
    Vector Hkr;                 // kinematic hardening ratio for axial, flexural end i and end j (default [0;0;0])

    static Matrix K;
    static Vector P;
    Vector Q;
    
    static Matrix k;
    Vector q;
    double q0[3];  // Fixed end forces in basic system
    double p0[3];  // Reactions in basic system
    
    // History variables
    Vector vppast;       // plastic deformations
    Vector vppres;
    Vector alphapast;
    Vector alphapres;
    Vector qbpast;
    Vector qbpres;
    
    Node *theNodes[2];
    
    ID  connectedExternalNodes;    

    CrdTransf *theCoordTransf;
};

#endif
