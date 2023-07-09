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

// $Revision: 1.1 $
// $Date: 2011-07-18 10:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/DispBeamColumn2dThermal.h,v $
                                                                        
// Written: MHS
// Created: Feb 2001
// Modified: Jian Zhang[University of Edinburgh]
// Modified: Panagiotis Kotsovinos[University of Edinburgh]
// Modified: Jian Jiang[University of Edinburgh]
// Modified: Liming Jiang[University of Edinburgh,2014




#ifndef DispBeamColumn2dThermal_h
#define DispBeamColumn2dThermal_h

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


class DispBeamColumn2dThermal : public Element
{
  public:
    DispBeamColumn2dThermal(int tag, int nd1, int nd2,
		     int numSections, SectionForceDeformation **s,
		     BeamIntegration &bi, CrdTransf &coordTransf,
		     double rho = 0.0);
    DispBeamColumn2dThermal();
    ~DispBeamColumn2dThermal();

    const char *getClassType(void) const {return "DispBeamColumn2dThermal";};

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
    int addLoad(ElementalLoad *theLoad, const Vector &loadFactors);

    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);

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
    
  private:
    const Matrix &getInitialBasicStiff(void);

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

	double *dataMix; //   temperature and location

	double q0Temperature[3];  // Fixed end thermal forces  of current step in basic system
	double q0TemperatureP[3];  // Fixed end thermal forces of last step in basic system
	int counterTemperature; // trace to remove thermal force from the second iteration step
    double SectionThermalElong[20];
    double AverageThermalElong;
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////

	//static double zz;
	//static double* kk;
	//Adding loadfactors to 'dispbeam2dThermal' for'FireLoadPattern' [-BEGIN-]: by L.J&P.K--8-May-2012--//
	double loadFactor2;
	double loadFactor3;
	double loadFactor4;
	double loadFactor5;
	double loadFactor6;
	double loadFactor7;
	double loadFactor8;
	double loadFactor9;
	//Adding loadfactors to 'dispbeam2dThermal' for'FireLoadPattern'  [-END-]: by L.J&P.K--8-May-2012--//
};

#endif

