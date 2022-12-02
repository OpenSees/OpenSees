/* ******************************************************************
***
**    OpenSees - Open System for Earthquake Engineering Simulation **
**          Pacific Earthquake Engineering Research Center **
** **
** **
** (C) Copyright 1999, The Regents of the University of California **
** All Rights Reserved. **
** **
** Commercial use of this program without express permission of the **
** University of California, Berkeley, is strictly prohibited.  See **
** file 'COPYRIGHT'  in main directory for information on usage and **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES. **
** **
** Developed by: **
**   Frank McKenna (fmckenna@ce.berkeley.edu) **
**   Gregory L. Fenves (fenves@ce.berkeley.edu) **
**   Filip C. Filippou (filippou@ce.berkeley.edu) **
** **
** ******************************************************************
*/

// $Revision $
// $Date$

// Written: Minjie Zhu

#ifndef PFEMContact3D_h
#define PFEMContact3D_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

#include <vector>

class Pressure_Constraint;

class PFEMContact3D : public Element {
   public:
    PFEMContact3D();

    PFEMContact3D(int tag, int nd1, int nd2, int nd3, int nd4,
                  int nd5, int nd6, int nd7, int nd8, double dc,
                  double e, double rho, double nx, double ny,
                  double nz, double ae);

    ~PFEMContact3D();

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes(void) const { return ntags.Size(); }

    const ID &getExternalNodes(void) { return ntags; }

    Node **getNodePtrs(void) { return &nodes[0]; }

    int getNumDOF(void) { return ndf.back(); }

    // public methods to set the state of the element
    int revertToLastCommit(void) { return 0; }

    int revertToStart(void) { return Element::revertToStart(); }

    int update(void);

    int commitState(void) { return Element::commitState(); }

    // public methods to obtain stiffness, mass, damping and residual
    // information
    const Matrix &getTangentStiff(void);

    const Matrix &getInitialStiff(void);

    const Matrix &getDamp();

    const Matrix &getMass(void);

    // methods for applying loads
    int addInertiaLoadToUnbalance(const Vector &accel);

    // methods for obtaining resisting force (force includes elemental
    // loads)
    const Vector &getResistingForce(void);

    const Vector &getResistingForceIncInertia(void);

    // MovableObject
    const char *getClassType(void) const;

    int sendSelf(int commitTag, Channel &theChannel);

    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    // DomainComponent
    void setDomain(Domain *theDomain);

    // TaggedObject
    void Print(OPS_Stream &s, int flag = 0);

    int displaySelf(Renderer &, int mode, float fact,
                    const char **displayModes = 0, int numModes = 0);

   private:
    double getD();

    double getP(double D);

    double getVI();

    bool inContact(double D);

   private:
    // Dc - initial distance
    // E - elastic modulus of debris
    // rho - density of debris
    // ndir - normal direction from 5678 to 1234

    ID ntags;
    std::vector<Node *> nodes;
    double Dc, E, rho, Ae, Fi;
    std::vector<int> ndf;
    Vector ndir;

    static Matrix K;
    static Vector P;
};

#endif
