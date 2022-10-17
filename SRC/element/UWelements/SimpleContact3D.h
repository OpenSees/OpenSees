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
                                                                        
#ifndef SimpleContact3D_h
#define SimpleContact3D_h

// Written: Kathryn A. Petek
// Created: 02/04
//
// Description: This file contains the class definition for SimpleContact3D. 
//

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

// number of nodes per element
#define SC3D_NUM_NODE 6
// d.o.f. per node
#define SC3D_NUM_NDF  3
// degrees of freedom per element
//#define SC3D_NUM_DOF  (SC3D_NUM_NDF * SC3D_NUM_NODE)
#define SC3D_NUM_DOF  18
// displacement degrees of freedom per element
//#define SC3D_NUM_DDOF  (SC3D_NUM_DOF - SC3D_NUM_NDF)
#define SC3D_NUM_DDOF  15

class Domain;
class Node;
class Channel;
class NDMaterial;
class ContactMaterial3D;
class FEM_ObjectBroker;


class SimpleContact3D : public Element
{
  public:
    SimpleContact3D(int tag, int Nd1, int Nd2, int Nd3, int Nd4, 
		int NdS, int NdL, NDMaterial &theMat, double tolG, double tolF);
    SimpleContact3D();
    ~SimpleContact3D();

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

    Vector project(Vector XiEta0);
    
    // method to obtain projection point on primary surface
    Vector GetPoint(Vector XiEta);
    // method to update base vectors g1 & g2
    int UpdateBase(Vector XiEta);

    void ComputeB(void);	  // function computes Bn, Bs @ step n


    ContactMaterial3D *theMaterial;  // contact material object

    ID  externalNodes;	  // contains the tags of the end nodes
    Matrix tangentStiffness;  // Tangent Stiffness matrix
    Vector internalForces;	  // vector of Internal Forces
    Vector theVector;         // vector to return the residual

    double tolGap;
	double tolForce;

    double gap;
    double lambda;

    bool inContact;
    bool was_inContact;
    bool to_be_released;
    bool should_be_released;
	bool in_bounds;

    Node *theNodes[SC3D_NUM_NODE];

    Vector d;

    Matrix x;				// matrix of cartesian coords of nodes 1-4, secondary
    Matrix g_metric;		// metric tensor
	Matrix G_metric;		// contravariant metric tensor

    Vector XiEta_n;			// holds curvlinear coords of xi, eta @ step n
    Vector XiEta_nplus1;	// holds curvlinear coords of xi, eta @ step n+1
    Vector x_XiCrd;			// holds cartesian coords of xi, eta
    
    Vector slip;
    
    Vector g1;		// tangent vector  = d(x_Xi)/d_xi
    Vector g2;		// tangent vector  = d(x_Xi)/d_eta 

    Vector n;		// normal Vector 

    double xsi_n;
    double xsi_nplus1;

    double N1;		// value of shape function 1
    double N2;		// value of shape function 2
    double N3;		// value of shape function 3
    double N4;		// value of shape function 4

	Matrix Kinv;	// inverse of K matrix: K = -dR/dXi for nonlinear projection
	Matrix KinvLin;	// inverse of K matrix for linear projection

    Vector Bn;		// gap-displacement matrix
    Matrix Bs;		// slip-displacement matrix

    Vector dcrd1;       // current coordinates of node 1
    Vector dcrd2;       // current coordinates of node 2
    Vector dcrd3;       // current coordinates of node 3
    Vector dcrd4;       // current coordinates of node 4
    Vector dcrdS;       // current coordinates of node S
    Vector dispL;       // current value of the Lagrangean multiplier

	int MyTag;          // what is my name?

};

#endif




