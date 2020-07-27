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
                                                                        
// $Revision: 1.0
// $Date:  
// $Source: /home/cee/pmackenz/Software/OpenSees/SRC/element/SimpleContact/SimpleContact2D.h,v $
                                                                        
#ifndef SimpleContact2D_h
#define SimpleContact2D_h

// Written: kap	
// Created: 02/04
//
// Description: This file contains the class definition for SimpleContact2D. 
//
// What: "@(#) SimpleContact2D.h, revA"

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>


// number of nodes per element
#define SC_NUM_NODE 4
// d.o.f. per node
#define SC_NUM_NDF  2
// degrees of freedom per element
//#define SC_NUM_DOF  (SC_NUM_NDF * SC_NUM_NODE)
#define SC_NUM_DOF  8
// displacement degrees of freedom per element
//#define SC_NUM_DDOF  (SC_NUM_DOF - SC_NUM_NDF)
#define SC_NUM_DDOF  6

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;
class ContactMaterial2D;

class SimpleContact2D : public Element
{
  public:
    SimpleContact2D(int tag, int Nd1, int Nd2,int NdS, int NdL, 
		NDMaterial &theMat, double tolG, double tolF);
    SimpleContact2D();
    ~SimpleContact2D();

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

    Response *setResponse(const char **argv, int argc, OPS_Stream &eleInfo);
    int getResponse(int responseID, Information &eleInformation);

	// public methods for material stage update
	int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

  protected:

  private:
    ContactMaterial2D *theMaterial;  // contact material object

    ID  externalNodes;	  // contains the tags of the end nodes
    Matrix tangentStiffness;  // Tangent Stiffness matrix
    Vector internalForces;	  // vector of Internal Forces
    Vector theVector;         // vector to return the residual

    double tolGap;
	double tolForce;

    double gap;
    double slip;
    double lambda;

    bool inContact;
    bool was_inContact;
    bool to_be_released;
    bool should_be_released;
    bool in_bounds;

    Node *theNodes[SC_NUM_NODE];

    double xsi_n;
    double xsi_nplus1;

    Vector n;		// normal Vector - 
                            // perpendicular to line between nodes 1 & 2
    Vector T;		// unit tangent vector (reference state)
    double Lprimary;		// Length of primary segment
    double Lsquare;		// square of Lprimary

    double N1;		// value of shape function 1
    double N2;		// value of shape function 2

    Vector Bn;		// gap-displacement matrix
    Vector Bs;		// slip-displacement matrix

    Vector dcrd1;       // current coordinates of nore 1
    Vector dcrd2;       // current coordinates of nore 2
    Vector dcrdS;       // current coordinates of nore S
    Vector dispL;       // current value of the Lagrangean multiplier

    int MyTag;          // what is my name?
};

#endif




