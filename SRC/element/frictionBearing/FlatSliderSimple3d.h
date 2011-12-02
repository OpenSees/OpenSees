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

// $Revision: 1.2 $
// $Date: 2009-11-03 23:12:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/FlatSliderSimple3d.h,v $

#ifndef FlatSliderSimple3d_h
#define FlatSliderSimple3d_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class definition for FlatSliderSimple3d.
// FlatSliderSimple3d is a friction slider element defined by two nodes. This
// simplified version uses small angle approximations and accounts for
// rotations of the sliding surface by shifting the shear forces.

#include <Element.h>
#include <Matrix.h>

class Channel;
class FrictionModel;
class UniaxialMaterial;
class Response;

class FlatSliderSimple3d : public Element
{
public:
    // constructor
    FlatSliderSimple3d(int tag, int Nd1, int Nd2,
        FrictionModel &theFrnMdl, double uy,
        UniaxialMaterial **theMaterials,
        const Vector y = 0, const Vector x = 0,
        double mass = 0.0, int maxIter = 20, double tol = 1E-8);
    FlatSliderSimple3d();
    
    // destructor
    ~FlatSliderSimple3d();
    
    // method to get class type
    const char *getClassType() const {return "FlatSliderSimple3d";};
    
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
    const Matrix &getMass();
    
    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    
    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();
    
    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    void Print(OPS_Stream &s, int flag = 0);    
    
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);
    
protected:

private:
    // private methods
    void setUp();
    double sgn(double x);
    
    // private attributes - a copy for each object of the class
    ID connectedExternalNodes;          // contains the tags of the end nodes
    Node *theNodes[2];                  // array of nodes
    FrictionModel *theFrnMdl;           // pointer to friction model
    UniaxialMaterial *theMaterials[4];  // array of uniaxial materials
    
    // parameters
    double uy;          // yield displacement
    Vector x;           // local x direction
    Vector y;           // local y direction
    double mass;        // mass of element
    int maxIter;        // maximum number of iterations
    double tol;         // tolerance for convergence criterion
    double L;           // element length
    
    // state variables
    Vector ub;          // displacements in basic system
    Vector ubPlastic;   // plastic displacements in basic system
    Vector qb;          // forces in basic system
    Matrix kb;          // stiffness matrix in basic system
    Vector ul;          // displacements in local system
    Matrix Tgl;         // transformation matrix from global to local system
    Matrix Tlb;         // transformation matrix from local to basic system
    
    // committed history variables
    Vector ubPlasticC;  // plastic displacements in basic system
    
    // initial stiffness matrix in basic system
    Matrix kbInit;
    
    static Matrix theMatrix;
    static Vector theVector;
    static Vector theLoad;
};

#endif
