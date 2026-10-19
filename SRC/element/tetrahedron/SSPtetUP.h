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

// ============================================================================
// 2026 By Pedro Arduino @ University of Washington, WA and
// Marlio Molina Manrique @ Politecnico di Torino, Turin, Italy
// parduino@uw.edu | marlio.molina@polito.it
// ============================================================================
// Experimental P1/P1 tetrahedral adaptation of McGann, Arduino and
// Mackenzie-Helnwein (2013), equations (70), (73), (79)-(81).
// No solid hourglass terms are required for an affine four-node tetrahedron.

#ifndef SSPtetUP_h
#define SSPtetUP_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class NDMaterial;
class Node;

class SSPtetUP : public Element {
public:
    SSPtetUP();
    SSPtetUP(int tag, int n1, int n2, int n3, int n4, NDMaterial &mat,
             double bulk, double fluidDensity, double kx, double ky, double kz,
             double voidRatio, double alpha, double bx=0, double by=0, double bz=0,
             bool lumped=false);
    ~SSPtetUP();
    const char *getClassType() const { return "SSPtetUP"; }
    int getNumExternalNodes() const { return 4; }
    const ID &getExternalNodes() { return nodes; }
    Node **getNodePtrs() { return node; }
    int getNumDOF() { return 16; }
    void setDomain(Domain *domain);
    int update();
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();
    const Matrix &getDamp();
    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();
    void zeroLoad();
    int addLoad(ElementalLoad *load, double factor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    int sendSelf(int commitTag, Channel &channel);
    int recvSelf(int commitTag, Channel &channel, FEM_ObjectBroker &broker);
    void Print(OPS_Stream &s, int flag=0);
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int id, Information &info);
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int id, Information &info);
private:
    const Matrix &formStiffness(const Matrix &tangent);
    ID nodes;
    Node *node[4];
    NDMaterial *material;
    double fBulk, fDensity, porosity, alpha, perm[3], body[3], appliedBody[3];
    double volume;
    bool lumped, loaded, valid;
    Matrix grad, B, stiffness, initial, mass, damp;
    Vector residual, unbalance;
};
#endif
