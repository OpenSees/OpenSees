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
// $Source$

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for DispBeamColumnNL2d.
// The element displacement field gives rise to constant axial strain and
// linear curvature.

#ifndef DispBeamColumnNL2d_h
#define DispBeamColumnNL2d_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <BeamIntegration.h>

class Node;
class SectionForceDeformation;
class CrdTransf;
class Response;

class DispBeamColumnNL2d : public Element
{
  public:
    DispBeamColumnNL2d(int tag, int nd1, int nd2,
		     int numSections, SectionForceDeformation **s,
		     BeamIntegration &bi, CrdTransf &coordTransf,
		     double rho = 0.0);
    DispBeamColumnNL2d();
    ~DispBeamColumnNL2d();

    const char *getClassType(void) const {return "DispBeamColumnNL2d";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    int update(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);
    int getResponseSensitivity(int responseID, int gradNumber,
			       Information &eleInformation);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int            updateParameter(int parameterID, Information &info);
    int            activateParameter(int parameterID);
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getInitialStiffSensitivity(int gradNumber);
    const Matrix & getMassSensitivity(int gradNumber);
    int            commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    const Matrix &getInitialBasicStiff(void);
    void getBasicStiff(Matrix &kb, int initial = 0);

    int numSections;
    SectionForceDeformation **theSections; // pointer to the ND material objects
    CrdTransf *crdTransf;        // pointer to coordinate transformation object 

    BeamIntegration *beamInt;

    ID connectedExternalNodes; // Tags of quad nodes

    Node *theNodes[2];

    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector

    Vector Q;		// Applied nodal loads
    Vector q;		// Basic force
    double q0[3];  // Fixed end forces in basic system
    double p0[3];  // Reactions in basic system

    double rho;			// Mass density per unit length

    enum {maxNumSections = 20};

    static double workArea[];

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////
};

#endif

