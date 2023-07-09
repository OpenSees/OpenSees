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

#ifndef ElastomericBearingUFRP2d_h
#define ElastomericBearingUFRP2d_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 10/13
// Revision: A
//
// Description: This file contains the class definition for ElastomericBearingUFRP2d.
// ElastomericBearingUFRP2d is an elastomeric bearing such as an unbonded fiber
// reinforced polymer bearing defined by two nodes. This simplified version uses
// a unidirectional Bouc-Wen model [2011, Love, Tait & Toopchi-Nezhad] to simulate
// the shear behavior and two uniaxial material models to simulate the axial and
// moment behaviors. Because the axial and shear springs are uncoupled the influence
// of the axial load on the shear behavior is not accounted for. However, the total
// P-Delta moment is equally distributed to the two end nodes of the element.

#include <Element.h>
#include <Matrix.h>

class Channel;
class UniaxialMaterial;
class Response;

class ElastomericBearingUFRP2d : public Element
{
public:
    // constructors
    ElastomericBearingUFRP2d(int tag, int Nd1, int Nd2, double uy,
        double a1, double a2, double a3, double a4, double a5,
        double b, double c, UniaxialMaterial **theMaterials,
        const Vector y = 0, const Vector x = 0,
        double eta = 1.0, double beta = 0.5, double gamma = 0.5,
        double shearDistI = 0.5, int addRayleigh = 0, double mass = 0.0,
        int maxIter = 25, double tol = 1E-12);
    ElastomericBearingUFRP2d();
    
    // destructor
    ~ElastomericBearingUFRP2d();
    
    // method to get class type
    const char *getClassType() const {return "ElastomericBearingUFRP2d";};
    
    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();
    void setDomain(Domain *theDomain);
    
    // public methods to set the state of the element
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();
    
    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getDamp();
    const Matrix &getMass();
    
    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    
    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();
    
    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag = 0);
    
    // public methods for element recorder
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);
    
protected:

private:
    // private methods
    void setUp();
    double sgn(double x);
    
    // private attributes - a copy for each object of the class
    ID connectedExternalNodes;          // contains the tags of the end nodes
    Node *theNodes[2];                  // array of nodes
    UniaxialMaterial *theMaterials[2];  // array of uniaxial materials
    
    // parameters
    double uy;          // yield displacement of hysteretic comp.
    double a1;          // coefficient for 1st order disp. term
    double a2;          // coefficient for 2nd order disp. term
    double a3;          // coefficient for 3rd order disp. term
    double a4;          // coefficient for 4th order disp. term
    double a5;          // coefficient for 5th order disp. term
    double b;           // coefficient for yield force of hysteretic comp.
    double c;           // coefficient for linear vel. term
    double eta;         // yielding exponent (sharpness of hysteresis loop corners)
    double beta;        // hysteretic shape parameter
    double gamma;       // hysteretic shape parameter
    double A;           // restoring force amplitude parameter
    Vector x;           // local x direction
    Vector y;           // local y direction
    double shearDistI;  // shear distance from node I as fraction of length
    int addRayleigh;    // flag to add Rayleigh damping
    double mass;        // mass of element
    int maxIter;        // maximum number of iterations
    double tol;         // tolerance for convergence criterion
    double L;           // element length
    bool onP0;          // flag to indicate if the element is on P0
    
    // state variables
    Vector ub;          // displacements in basic system
    double z;           // hysteretic evolution parameter
    double dzdu;        // tangent of hysteretic evolution parameter
    Vector qb;          // forces in basic system
    Matrix kb;          // stiffness matrix in basic system
    Vector ul;          // displacements in local system
    Matrix Tgl;         // transformation matrix from global to local system
    Matrix Tlb;         // transformation matrix from local to basic system
    
    // committed history variables
    Vector ubC;         // displacement in basic system
    double zC;          // hysteretic evolution parameter
    
    // initial stiffness matrix in basic system
    Matrix kbInit;
    
    static Matrix theMatrix;  // a class wide Matrix
    static Vector theVector;  // a class wide Vector
    Vector theLoad;
};

#endif
