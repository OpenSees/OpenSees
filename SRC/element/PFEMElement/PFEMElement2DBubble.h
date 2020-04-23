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

// $Revision: 1.00 $
// $Date: 2012/01/11 13:48:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/PFEMElement2DBubble.h,v $

// Written: Minjie Zhu (zhum@engr.orst.edu)
// Created: Jan 2012
// Revised: --------
//
// Description: This file contains the class definition for PFEMElement2DBubble.

#ifndef PFEMElement2DBubble_h
#define PFEMElement2DBubble_h

#include <Matrix.h>
#include <Vector.h>
#include <Element.h>

class Pressure_Constraint;

class PFEMElement2DBubble : public Element
{
public:
    PFEMElement2DBubble();
    PFEMElement2DBubble(int tag, int nd1, int nd2, int nd3,
                        double r, double m, double b1, double b2,
                        double thk=1.0, double ka=-1);

    ~PFEMElement2DBubble();

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
    const Matrix& getDampSensitivity(int gradNumber);
    const Matrix& getMassSensitivity(int gradNumber);
    const Vector& getResistingForceSensitivity(int gradNumber);
    int commitSensitivity(int gradNumber, int numGrads);

    static bool dispon;

protected:

private:

    ID ntags; // Tags of nodes
    Node* nodes[6]; // pointers of nodes
    Pressure_Constraint* thePCs[3];
    double rho;  // density
    double mu;   // viscocity
    double bx;    // body force
    double by;    // body force
    double J;
    Vector dJ;
    ID numDOFs;
    double thickness;
    double kappa;
    int parameterID;

    static Matrix K;
    static Vector P;
    static Matrix C;
    Matrix M, D;
    Vector F, Fp;

    void setJ();
    void setdJ();
    static void setC();
    int updateMatrix();

    // responses
    double getM() const;
    double getMp() const;
    void getG(Vector& g) const;
    void getK(Matrix& k) const;
    void getKbub(Matrix& kbub) const;
    void getGbub(Matrix& gbub)const;
    double getMbub() const;
    void getL(Matrix& l) const;
    void getF(Vector& f) const;
    void getFbub(Vector& fbub) const;
    void getFp(Vector& fp) const;

    // conditional sensitivity
    double getdM() const;
    double getdinvMbub() const;
    void getdL(Matrix& dl) const;
    void getdF(Vector& df) const;
    void getdK(Matrix& dk) const;
    void getdFbub(Vector& dfb) const;
    void getdFp(Vector& dfp) const;

    // geometric sensitivity
    void getdM(const Vector& vdot, Matrix& dm) const;
    void getdinvMbub(const Vector& vb, Matrix& dmb) const;
    void getdK(const Vector& v, Matrix& dk) const;
    void getdF(Matrix& df) const;
    void getdFbub(Matrix& dfb) const;
    void getdG(const Vector& p, Matrix& dg) const;
    void getdGt(const Vector& v, Matrix& dgt) const;
    void getdGb(const Vector& p, Matrix& dgb) const;
    void getdGbt(const Vector& vb, Matrix& dgbt) const;
    void getdFp(Matrix& dfp) const;
    void getdL(const Vector& p, Matrix& dl) const;
};

#endif
