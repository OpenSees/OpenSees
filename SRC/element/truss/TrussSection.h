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
                                                                        
                                                                        
#ifndef TrussSection_h
#define TrussSection_h

// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for Truss. A Truss 
// object provides the abstraction of the small deformation bar element. 
// Each truss object is associated with a section object. This Truss element 
// will work in 1d, 2d or 3d problems.
//
// What: "@(#) Truss.h, revA"

#include <Element.h>
#include <Matrix.h>

class Node;
class Channel;
class SectionForceDeformation;

class TrussSection : public Element
{
  public:
    TrussSection(int tag, int dimension,
		 int Nd1, int Nd2, 
		 SectionForceDeformation &theSection,
		 double rho = 0.0,
		 int doRayleighDamping = 0,
         int cMass = 0);
    
    TrussSection();    
    ~TrussSection();

    const char *getClassType(void) const {return "TrussSection";};

    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    
    
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool tag);
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getKiSensitivity(int gradNumber);
    const Matrix & getMassSensitivity(int gradNumber);
    int commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    double computeCurrentStrain(void) const;
    
    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes;          // contains the tags of the end nodes
    int dimension;                       // truss in 2 or 3d domain
    int numDOF;	                         // number of dof for truss

    Vector *theLoad;    // pointer to the load vector P
    Matrix *theMatrix;  // pointer to objects matrix (a class wide Matrix)
    Vector *theVector;  // pointer to objects vector (a class wide Vector)

    double cosX[3];     // direction cosines

    double L;               // length of truss based on undeformed configuration
    double rho;             // mass density per unit length
    int doRayleighDamping;  // flag to include Rayleigh damping
    int cMass;              // consistent mass flag

    Node *theNodes[2];
    double *initialDisp;

    SectionForceDeformation  *theSection;
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    Vector *theLoadSens;
    // AddingSensitivity:END ///////////////////////////////////////////

    // static data - single copy for all objects of the class	
    static Matrix trussM2;   // class wide matrix for 2*2
    static Matrix trussM3;   // class wide matrix for 3*3	
    static Matrix trussM4;   // class wide matrix for 4*4
    static Matrix trussM6;   // class wide matrix for 6*6
    static Matrix trussM12;  // class wide matrix for 12*12
    static Vector trussV2;   // class wide Vector for size 2
    static Vector trussV3;   // class wide Vector for size 3	
    static Vector trussV4;   // class wide Vector for size 4
    static Vector trussV6;   // class wide Vector for size 6
    static Vector trussV12;  // class wide Vector for size 12
};

#endif
