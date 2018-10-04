// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumnInt/DispBeamColumn2dInt.h,v $

// $Revision: 1.4 $

// $Date: 2008-04-14 21:22:20 $



// Created: 07/04

// Modified by: LMS 

// Description: This file contains the class implementation of DispBeamColumn2dInt.Based on DispBeamColumn2d.cpp.



  

#ifndef DispBeamColumn2dInt_h

#define DispBeamColumn2dInt_h



#ifndef _bool_h

#include "bool.h"

#endif



#include <Element.h>

#include <Matrix.h>

#include <Vector.h>

#include <ID.h>

#include <LegendreBeamIntegration.h>

class Node;

class SectionForceDeformation;

class FiberSection2dInt;

class CrdTransf;

class LinearCrdTransf2dInt;

class LinearCrdTransf2dInt;

class Response;



class DispBeamColumn2dInt : public Element

{

  public:

    DispBeamColumn2dInt(int tag, 

			int nd1, 

			int nd2,

			int numSections, 

			SectionForceDeformation **s,

			CrdTransf &coordTransf, 

			double C, 

			double rho = 0.0);



    DispBeamColumn2dInt();

    virtual ~DispBeamColumn2dInt();



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

    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag =0);



    Response *setResponse(const char **argv, int argc, OPS_Stream &s);

    int getResponse(int responseID, Information &eleInfo);



    // AddingSensitivity:BEGIN //////////////////////////////////////////
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getKiSensitivity(int gradNumber);
    const Matrix & getMassSensitivity(int gradNumber);
    int            commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    const Matrix &getInitialBasicStiff(void);

    int numSections;

    FiberSection2dInt **theSections; // pointer to the ND material objects
    LinearCrdTransf2dInt *crdTransf;          // pointer to coordinate transformation object 

    double C1;
    ID connectedExternalNodes; // Tags of quad nodes

    Node *theNodes[2];

    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector

    Vector Q;		// Applied nodal loads
    Vector q;		// Basic force
    double q0[6];  
    double p0[6];  

    double rho;			// Mass density per unit length

    static double workArea[];

    static LegendreBeamIntegration quadRule;

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////
};

#endif


