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

// $Revision $
// $Date$

// Written: Minjie Zhu

#ifndef PFEMContact2D_h
#define PFEMContact2D_h

#include <Matrix.h>
#include <Vector.h>
#include <Element.h>
#include <vector>

class Pressure_Constraint;

class PFEMContact2D : public Element {
public:
    PFEMContact2D();

    PFEMContact2D(int tag, int nd1, int nd2, int nd3,
                  double kk, double thk,
                  double m, double b, double dc,
                  double a, double e, double rho);

    ~PFEMContact2D();

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes(void) const { return ntags.Size(); }

    const ID &getExternalNodes(void) { return ntags; }

    Node **getNodePtrs(void) { return &nodes[0]; }

    int getNumDOF(void) { return ndf[3]; }

    // public methods to set the state of the element
    int revertToLastCommit(void) { return 0; }

    int revertToStart(void) { return Element::revertToStart(); }

    int update(void);

    int commitState(void) { return Element::commitState(); }

    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix &getTangentStiff(void);

    const Matrix &getInitialStiff(void);

    const Matrix &getDamp();

    const Matrix &getMass(void);

    // methods for applying loads
    int addInertiaLoadToUnbalance(const Vector &accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector &getResistingForce(void);

    const Vector &getResistingForceIncInertia(void);

    // MovableObject
    const char *getClassType(void) const;

    int sendSelf(int commitTag, Channel &theChannel);

    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // DomainComponent
    void setDomain(Domain *theDomain);

    // TaggedObject
    void Print(OPS_Stream &s, int flag = 0);

    int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);

private:

    double getLine(double &A, double &B, double &C,
                   double &dx, double &dy,
                   double &x1, double &y1,
                   double &x2, double &y2,
                   double &x3, double &y3,
                   double &L);

    static void getdL(double L, double dx, double dy, Vector &dL);

    static void getdA(double L, double A, const Vector &dL, Vector &dA);

    static void getdB(double L, double B, const Vector &dL, Vector &dB);

    static void getdC(double L, double C, double x1, double y1,
                      double x2, double y2, const Vector &dL,
                      Vector &dC);

    static void getdD(double A, double B, double x3, double y3,
                      const Vector &dA,
                      const Vector &dB,
                      const Vector &dC,
                      Vector &dD);

    double getP(double D);

    void getdP(const Vector &dD, double D, Vector &dP);

    void getV(Vector &vn, double &vt, Vector &vj);

    bool inContact(double D);


private:

    // kdoverAd = E / Ld, Ld - length of debris
    // thk - thickness of debris, out of plane width
    // mu - damping ratio
    // beta -frictional coefficient
    // Dc - initial distance between node 3 and edge 1-2, usually mesh size
    // alpha - stiffness parameter
    // E - elastic modulus of debris, the average density over debris volume
    // rho - density of debris, the average density over debris volume

    ID ntags;
    std::vector<Node *> nodes;
    double kdoverAd, thk, mu, beta, Dc, alpha, E, rho;
    std::vector<int> ndf;
    double signvt0, F0;
    static Matrix K;
    static Vector P;
};

#endif


