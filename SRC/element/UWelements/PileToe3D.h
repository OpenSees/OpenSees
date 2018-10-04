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
                                                                       
                                                                       
#ifndef PileToe3D_h
#define PileToe3D_h

// Written: Pedro Arduino
// Created: 09/08/14
//
// Description: This file contains the class definition for PileToe.

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>

#include <ID.h>
// #include <myDebugInfo.h>

class CrdTransf;

// number of nodes per element
#define PT3D_NUM_NODE 1
// number of dimesnions
#define PT3D_NUM_NDM  3
// degrees of freedom per element
#define PT3D_NUM_DOF  6

class Domain;
class Node;
class Channel;
class FEM_ObjectBroker;


class PileToe3D : public Element
{
  public:
    PileToe3D(int tag, int Nd1, int BNd1, int BNd2, double rad, double k, CrdTransf &coordTransf);
    PileToe3D();
    ~PileToe3D();

    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes(void) const;
	int getNumExternalBNodes(void) const;
    const ID &getExternalNodes(void);
	const ID &getExternalBNodes(void);
    Node **getNodePtrs(void);
	Node **getBNodePtrs(void);
    int getNumDOF(void);       
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
   
    // public methods to obtain stiffness, mass, damping and
    // residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    

    void zeroLoad(void);       
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);

    void Print(OPS_Stream &s, int flag =0);    

    // Response *setResponse(const char **argv, int argc, Information &eleInfo);
    Response *setResponse(const char **argv, int argc, OPS_Stream &eleInfo);
    int getResponse(int responseID, Information &eleInformation);

  protected:
   
  private:

    // member functions
   
    // objects
    CrdTransf  *crdTransf;              // pointer to coordinate transformation object   
    ID  externalNodes;                  // contains the tags of the end nodes
	ID  externalBNodes;                 // contains the tags of the bean element used to determine normal to pile toe
    Vector theVector;                   // vector to return the residual
    Matrix mTangentStiffness;           // Tangent Stiffness matrix
    Vector mInternalForces;             // vector of Internal Forces
    Node *theNodes[PT3D_NUM_NODE];
	Node *theBNodes[2];
       
    // input quantities
    double mRadius;                     // radius of Pile Surface   
    double mSubgradeCoeff;              // Rock/Soil subgrade coefficient at pile toe
    double mCC;        	                // mCC>=mRadius --> Complete area, mCC<mRadius --> Partial area
    int MyTag;                          // element tag for debugging

    // boolean variables
    bool mInitialize;

};

#endif
