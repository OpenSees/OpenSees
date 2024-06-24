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

// $Revision$
// $Date$

// Written: Minjie Zhu (zhum@oregonstate.edu)
//
// Description: This file contains the class definition for PFEMElement3DBubble.

#ifndef PFEMElement3DBubble_h
#define PFEMElement3DBubble_h

#include <Matrix.h>
#include <Vector.h>
#include <Element.h>
#include <vector>

class Pressure_Constraint;

class PFEMElement3DBubble : public Element
{
public:
    PFEMElement3DBubble();
    PFEMElement3DBubble(int tag, int nd1, int nd2, int nd3, int nd4,
                        double r, double m, double b1, double b2, double b3,
                        double ka=-1);

    ~PFEMElement3DBubble();

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);

    // public methods to set the state of the element
    int revertToLastCommit(void);
    //int revertToStart(void);
    int update(void);
    int commitState(void);

    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp();
    const Matrix &getMass(void);

    // methods for applying loads
    void zeroLoad();
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
    void Print(OPS_Stream &s, int flag =0);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

    // sensitivity
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);
    int activateParameter(int passedParameterID);

    static bool dispon;

    static double det(const Matrix& m);
    static void cofactor(const Matrix& m, Matrix& mfact);

protected:

private:

    ID ntags; // Tags of nodes
    std::vector<Node*> nodes; // pointers of nodes
    std::vector<Pressure_Constraint*> thePCs;
    double rho;  // density
    double mu;   // viscocity
    double bx;    // body force
    double by;    // body force
    double bz;    // body force
    double J;
    ID numDOFs;
    double kappa;
    int parameterID;
    std::vector<double> dNdx, dNdy, dNdz;

    static Matrix K;
    static Vector P;
    Matrix M, D;
    Vector F, Fp;
    Vector Q;

    int updateJacobi();
    int updateMatrix();

    // responses
    double getM() const;
    double getMp() const;
    void getG(Matrix& g) const;
    void getK(Matrix& k) const;
    void getGbub(Matrix& gbub)const;
    double getinvMbub() const;
    void getL(Matrix& l) const;
    void getF(Vector& f) const;
    void getFbub(Vector& fbub) const;
    void getFp(Vector& fp) const;

};

#endif
