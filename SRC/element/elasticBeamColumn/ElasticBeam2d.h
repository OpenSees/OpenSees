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
                                                                        
                                                                        
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ElasticBeam2d.
// ElasticBeam2d is a plane frame member.

#ifndef ElasticBeam2d_h
#define ElasticBeam2d_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class Information;
class CrdTransf;
class Damping;
class SectionForceDeformation;
class Response;
class Renderer;

class ElasticBeam2d : public Element
{
  public:
    ElasticBeam2d();        
    ElasticBeam2d(int tag, double A, double E, double I, 
		  int Nd1, int Nd2, CrdTransf &theTransf,
		  double alpha = 0.0, double d = 0.0,
		  double rho = 0.0, int cMass = 0,
		  int release = 0, Damping *theDamping = 0);
    ElasticBeam2d(int tag, double A, double E, double I,
		  double G, double Av,
		  int Nd1, int Nd2, CrdTransf &theTransf,
		  double alpha = 0.0, double d = 0.0,
		  double rho = 0.0, int cMass = 0,
		  int release = 0, Damping *theDamping = 0);  
    ElasticBeam2d(int tag, int Nd1, int Nd2, 
		  SectionForceDeformation& theSection, CrdTransf &theTransf,
		  double alpha = 0.0, double d = 0.0,
		  double rho = 0.0, int cMass = 0,
		  int release = 0, Damping *theDamping = 0);

    ~ElasticBeam2d();

    const char *getClassType(void) const {return "ElasticBeam2d";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    int setDamping(Domain *theDomain, Damping *theDamping);
    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);
    
    int update(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getDampingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);    
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes = 0, int numModes = 0);

    Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    int getResponse (int responseID, Information &info);
 
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int activateParameter(int parameterID);
    const Vector &getResistingForceSensitivity(int gradNumber);
    const Matrix &getTangentStiffSensitivity(int gradNumber);
    const Matrix &getInitialStiffSensitivity(int gradNumber);
    const Matrix &getCommittedStiffSensitivity(int gradNumber);
    const Matrix &getMassSensitivity(int gradNumber);
    int commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  private:
    // Derivative of the basic stiffness kb with respect to the active
    // parameter (E, A, I, G or Av); returns false if kb does not depend on it
    bool getBasicStiffSensitivity(Matrix &dkb);
    // Derivatives of the fixed-end forces q0 and reactions p0
    void getFixedEndForceSensitivity(Vector &dq0, Vector &dp0);

    int parameterID = 0;  // active DDM parameter (0 = none)
    // Load terms kept for the sensitivity of q0 and p0:
    double mEB[2] = {0.0, 0.0}; // Euler-Bernoulli end moments of transverse
                                // loads, before the shear/release map
    double mT[2] = {0.0, 0.0};  // temperature end moments divided by E*I
    double nT = 0.0;            // temperature axial force divided by E*A
    double A,E,I;     // area, elastic modulus, moment of inertia
  double G,Av; // properties for shear
    double alpha, d;  // coeff. of thermal expansion, depth
    double rho;       // mass per unit length
    int cMass;        // consistent mass flag

    int release;      // moment release 0=none, 1=I, 2=J, 3=I,J
    
    static Matrix K;
    static Vector P;
    Vector Q;
    
    static Matrix kb;
    Vector q;
    double q0[3];  // Fixed end forces in basic system
    double p0[3];  // Reactions in basic system
    
    Node *theNodes[2];
    
    ID  connectedExternalNodes;    

    CrdTransf *theCoordTransf;

    Damping *theDamping;
};

#endif
