/* *******************************************************************************
Copyright (c) 2016-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

********************************************************************************* */
                                                                        
#ifndef AxEqDispBeamColumn2d_h
#define AxEqDispBeamColumn2d_h

#ifndef _bool_h
#include "bool.h"
#endif



// Written: Danilo Tarquini 
//
// Description: This file contains the interface for the AxEqDispBeamColumn2d class.
// It defines the class interface and the class attributes. AxEqDispBeamColumn2d 
// provides the abstraction of a simple displacement based element for 2-d problems.
// The sectional relationship for the element being performed in the
// associated SectionForceDeformation object.
//
// What: "@(#) AxEqDispBeamColumn2d.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <SectionForceDeformation.h>

class Node;
class CrdTransf;
class Response;
class BeamIntegration;


class AxEqDispBeamColumn2d : public Element  // this class is a subclass of the pure virtual Element class
{
public:
	// CONSTRUCTORS AND DESTRUCTOR
	AxEqDispBeamColumn2d(int tag, int nd1, int nd2,
		int numSections, SectionForceDeformation **s,
		BeamIntegration &bi, CrdTransf &coordTransf, double tol,
		double rho = 0.0, int cMass = 0, int maxIters=20);  // Constructor 1 with arguments
	AxEqDispBeamColumn2d();  // Constructor 2
	~AxEqDispBeamColumn2d(); // Destructor

					  // METHODS TO RETRIEVE MAIN ELEMENT TYPE AND CONNECTIVITY
	const char *getClassType(void) const { return "AxEqDispBeamColumn2d"; };

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
	int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes = 0, int numModes = 0);
	void Print(OPS_Stream &s, int flag = 0);

	Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	int getResponse(int responseID, Information &eleInfo);
	int getResponseSensitivity(int responseID, int gradNumber,
		Information &eleInformation);

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int setParameter(const char **argv, int argc, Parameter &param);
	int            updateParameter(int parameterID, Information &info);
	int            activateParameter(int parameterID);
	const Vector & getResistingForceSensitivity(int gradNumber);
	const Matrix & getKiSensitivity(int gradNumber);
	const Matrix & getMassSensitivity(int gradNumber);
	int            commitSensitivity(int gradNumber, int numGrads);
	// AddingSensitivity:END ///////////////////////////////////////////

protected:
	// PRIVATE METHODS
private:
	const Matrix &getInitialBasicStiff(void);
	void getBasicStiff(Matrix &kb, int initial = 0);
	double getSectionalAxialForceUnbalance(); // METHOD ADDED BY DANILO, returns the axial force unbalance
	Vector getAxialStrainIncrement(); // METHOD ADDED BY DANILO, returns the axial strain increment vector

	// CLASS MEMBERS
	int numSections;
	SectionForceDeformation **theSections; // pointer to the ND material objects
	CrdTransf *crdTransf;        // pointer to coordinate transformation object 

	BeamIntegration *beamInt;

	ID connectedExternalNodes; // Tags of quad nodes

	Node *theNodes[2];

	static Matrix K;		// Element stiffness, damping, and mass Matrix
	static Vector P;		// Element resisting force vector

	Vector Q;      // Applied nodal loads
	Vector q;      // Basic force
	double q0[3];  // Fixed end forces in basic system
	double p0[3];  // Reactions in basic system

	double rho;	   // Mass density per unit length
	int cMass;     // consistent mass flag
	double tol;   // added by DANILO: tolerance on the axial force unbalance 
	int maxIters; // added by DANILO: number of maximum internal iterations allowed
	int flagDBae; // flag indicating that the DBae is used
	Vector vCommitted;
	Vector eCommitted;
	Vector e0Committed;
	Vector curvCommitted;

	enum { maxNumSections = 20 };

	static double workArea[];

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int parameterID;
	// AddingSensitivity:END ///////////////////////////////////////////
};
#endif

