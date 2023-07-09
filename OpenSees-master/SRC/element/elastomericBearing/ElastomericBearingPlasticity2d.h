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

#ifndef ElastomericBearingPlasticity2d_h
#define ElastomericBearingPlasticity2d_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class definition for ElastomericBearingPlasticity2d.
// ElastomericBearingPlasticity2d is an elastomeric bearing such as a lead-rubber bearing
// or a high-damping rubber bearing defined by two nodes. This simplified version
// uses a unidirectional plasticity model to simulate the shear behavior and two
// uniaxial material models to simulate the axial and moment behaviors. Because
// the axial and shear springs are uncoupled the influence of the axial load on
// the shear behavior is not accounted for. However, the total P-Delta moment is
// equally distributed to the two end nodes of the element.

#include <Element.h>
#include <Matrix.h>

class Channel;
class UniaxialMaterial;
class Response;

class ElastomericBearingPlasticity2d : public Element
{
public:
    // constructors
    ElastomericBearingPlasticity2d(int tag, int Nd1, int Nd2,
        double kInit, double qd, double alpha1,
        UniaxialMaterial **theMaterials,
        const Vector y = 0, const Vector x = 0,
        double alpha2 = 0.0, double mu = 2.0,
        double shearDistI = 0.5,
        int addRayleigh = 0, double mass = 0.0);
    ElastomericBearingPlasticity2d();
    
    // destructor
    ~ElastomericBearingPlasticity2d();
    
    // method to get class type
    const char *getClassType() const {return "ElastomericBearingPlasticity2d";};
    
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
    double k0;          // initial stiffness of hysteretic component
    double qYield;      // yield force of hysteretic component
    double k2;          // stiffness of elastic component
    double k3;          // stiffness of nonlinear elastic component
    double mu;          // exponent of nonlinear elastic component
    Vector x;           // local x direction
    Vector y;           // local y direction
    double shearDistI;  // shear distance from node I as fraction of length
    int addRayleigh;    // flag to add Rayleigh damping
    double mass;        // mass of element
    double L;           // element length
    bool onP0;          // flag to indicate if the element is on P0
    
    // state variables
    Vector ub;          // displacements in basic system
    double ubPlastic;   // plastic displacement in basic system
    Vector qb;          // forces in basic system
    Matrix kb;          // stiffness matrix in basic system
    Vector ul;          // displacements in local system
    Matrix Tgl;         // transformation matrix from global to local system
    Matrix Tlb;         // transformation matrix from local to basic system
    
    // committed history variables
    double ubPlasticC;  // plastic displacement in basic system
    
    // initial stiffness matrix in basic system
    Matrix kbInit;
    
    static Matrix theMatrix;  // a class wide Matrix
    static Vector theVector;  // a class wide Vector
    Vector theLoad;
};

#endif
