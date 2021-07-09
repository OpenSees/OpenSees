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

#ifndef InertiaTruss_h
#define InertiaTruss_h


// Description: This file contains the class definition for InertiaTruss. A InertiaTruss object
// provides the abstraction of the small deformation bar element. Each truss
// object is assocaited with a material object. This InertiaTruss element will work
// in 1d, 2d or 3d problems.
// Not ready for sensitivity analysis yet
// modified from truss.cpp    author: Frank McKenna
/*
**************************************************************************************************************
** InertiaTruss Element                                                                                     **
** First published in : Ji X, Cheng Y, Molina Hutt C.                                                       **
** Seismic response of a tuned viscous mass damper(TVMD) coupled wall system.Eng Struct 2020; 225:111252.   **
** https ://doi.org/10.1016/j.engstruct.2020.111252.                                                        **
**************************************************************************************************************
*/


#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
class Node;
class Channel;
class UniaxialMaterial;

class InertiaTruss : public Element
{
  public:
    InertiaTruss(int tag, int dimension,
	  int Nd1, int Nd2, 
	  double mr);
    
    InertiaTruss();    
    ~InertiaTruss();

    const char *getClassType(void) const {return "InertiaTruss";};

    // public methods to obtain inforrmation about dof & connectivity    
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
    const Matrix &getKi(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void); 
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


    // **************************************
    // Not ready for sensitivity analysis yet
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int		   addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool tag);
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getKiSensitivity(int gradNumber);
    const Matrix & getMassSensitivity(int gradNumber);
    int            commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    double computeCurrentStrain(void) const;
    double computeCurrentStrainRate(void) const;
    
    // private attributes - a copy for each object of the class
    UniaxialMaterial *theMaterial;  // pointer to a material
    ID  connectedExternalNodes;     // contains the tags of the end nodes
    int dimension;                  // inertiatruss in 2 or 3d domain
    int numDOF;	                    // number of dof for truss

    Vector *theLoad;    // pointer to the load vector P
    Matrix *theMatrix;  // pointer to objects matrix (a class wide Matrix)
    Vector *theVector;  // pointer to objects vector (a class wide Vector)

    double L;               // length of Inertiatruss based on undeformed configuration
    double mass;             // inertial mass


    double cosX[3];  // direction cosines

    Node *theNodes[2];
    double *initialDisp;

	
// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    Vector *theLoadSens;
// AddingSensitivity:END ///////////////////////////////////////////

    // static data - single copy for all objects of the class	
    static Matrix trussM2;   // class wide matrix for 2*2
    static Matrix trussM4;   // class wide matrix for 4*4
    static Matrix trussM6;   // class wide matrix for 6*6
    static Matrix trussM12;  // class wide matrix for 12*12
    static Vector trussV2;   // class wide Vector for size 2
    static Vector trussV4;   // class wide Vector for size 4
    static Vector trussV6;   // class wide Vector for size 6
    static Vector trussV12;  // class wide Vector for size 12
};

#endif
