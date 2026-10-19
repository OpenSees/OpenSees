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
// 2026 By  Pedro Arduino @ University of Washington, WA.
// parduino@uw.edu
// ============================================================================
// Experimental P1/P1 plane-strain triangle for saturated porous media.
// Based on McGann, Arduino and Mackenzie-Helnwein (2011), efficientElementPaper.pdf,
// equations (73)-(78), SSPquadUP, and the simplex adaptation in SSPtetUP.
// The constant-strain solid has no hourglass modes. Pressure-rate gradient
// stabilization does not remove locking of the drained solid skeleton.
// Nodes: ux, uy, integral(p dt); nodal velocity DOF 3 is pore pressure.
#ifndef Tri31UP_h
#define Tri31UP_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class NDMaterial;
class Node;

class Tri31UP : public Element {
public:
    Tri31UP();
    Tri31UP(int tag, int n1, int n2, int n3, NDMaterial &mat,
            double thickness, double bulk, double fluidDensity,
            double kx, double ky, double voidRatio, double alpha,
            double bx=0, double by=0, bool lumped=false);
    ~Tri31UP() override;
    const char *getClassType() const override { return "Tri31UP"; }
    int getNumExternalNodes() const override { return NumNodes; }
    const ID &getExternalNodes() override { return nodes; }
    Node **getNodePtrs() override { return node; }
    int getNumDOF() override { return Dofs; }
    void setDomain(Domain *domain) override;
    int update() override;
    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;
    const Matrix &getTangentStiff() override;
    const Matrix &getInitialStiff() override;
    const Matrix &getMass() override;
    const Matrix &getDamp() override;
    const Vector &getResistingForce() override;
    const Vector &getResistingForceIncInertia() override;
    void zeroLoad() override;
    int addLoad(ElementalLoad *load, double factor) override;
    int addInertiaLoadToUnbalance(const Vector &accel) override;
    int sendSelf(int commitTag, Channel &channel) override;
    int recvSelf(int commitTag, Channel &channel, FEM_ObjectBroker &broker) override;
    void Print(OPS_Stream &s, int flag=0) override;
    int displaySelf(Renderer &viewer, int mode, float fact,
                    const char **modes=0, int numModes=0) override;
    Response *setResponse(const char **argv, int argc, OPS_Stream &s) override;
    int getResponse(int id, Information &info) override;
    int setParameter(const char **argv, int argc, Parameter &param) override;
    int updateParameter(int id, Information &info) override;
private:
    enum { NumNodes=3, Dim=2, Ndf=3, Dofs=9, SolidDofs=6, Stress=3 };
    const Matrix &formStiffness(const Matrix &tangent);
    ID nodes;
    Node *node[NumNodes];
    NDMaterial *material;
    double thickness, fBulk, fDensity, porosity, alpha;
    double perm[Dim], body[Dim], appliedBody[Dim], volume; // volume = area * thickness
    bool lumped, loaded, valid;
    Matrix grad, B, stiffness, initial, mass, damp;
    Vector residual, unbalance;
};
#endif
