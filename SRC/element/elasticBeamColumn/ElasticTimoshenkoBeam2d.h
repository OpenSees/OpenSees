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
// $URL$

#ifndef ElasticTimoshenkoBeam2d_h
#define ElasticTimoshenkoBeam2d_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 03/13
// Revision: A
//
// Purpose: This file contains the class definition for ElasticTimoshenkoBeam2d.
// ElasticTimoshenkoBeam2d is a frame member that takes shearing deformations
// and rotational inertia effects into account.

#include <Element.h>
#include <Matrix.h>

class CrdTransf;

class ElasticTimoshenkoBeam2d : public Element
{
public:
    // constructor
    ElasticTimoshenkoBeam2d(int tag, int Nd1, int Nd2, double E, double G,
        double A, double Iz, double Avy, CrdTransf &theTransf,
        double rho = 0.0, int cMass = 0);
    ElasticTimoshenkoBeam2d();
    
    // destructor
    ~ElasticTimoshenkoBeam2d();
    
    // method to get class type
    const char *getClassType() const {return "ElasticTimoshenkoBeam2d";};
    
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
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

    void Print(OPS_Stream &s, int flag = 0);
    
    // public methods for element recorder
    Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    int getResponse (int responseID, Information &info);
    
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);
    
protected:
    
private:
    // private methods
    void setUp();
    
    // private attributes - a copy for each object of the class
    ID connectedExternalNodes;  // contains the tags of the end nodes
    Node *theNodes[2];          // array of nodes
    CrdTransf *theCoordTransf;  // coordinate transformation
    
    // parameters
    double E;    // elastic modulus
    double G;    // shear modulus
    double A;    // cross section (axial) area
    double Iz;   // moment of inertia about local z axis
    double Avy;  // shear area along local y axis
    double rho;  // mass per unit length
    int cMass;   // consistent mass flag
    int nlGeo;   // nonlinear geometry flag
    double phi;  // ratio of bending to shear stiffness
    double L;    // element length
    
    // state variables
    Vector ul;     // displacements in local system
    Vector ql;     // forces in local system
    Vector ql0;    // fixed end forces due to loads in local system
    Matrix kl;     // stiffness matrix in local system
    Matrix klgeo;  // geometric stiffness matrix in local system
    
    // constant variables
    Matrix Tgl;  // transformation matrix from global to local system
    Matrix Ki;   // initial stiffness matrix in global system
    Matrix M;    // mass matrix in global system
    
    static Matrix theMatrix;  // a class wide Matrix
    static Vector theVector;  // a class wide Vector
    Vector theLoad;
};

#endif
