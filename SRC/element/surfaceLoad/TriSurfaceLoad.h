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
                                                                        
#ifndef TriSurfaceLoad_h
#define TriSurfaceLoad_h

// Written: J. Abell
// Created: March 2018
// Modified: 

// Description: This file contains the class definition for TriSurfaceLoad. 

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>
#include <ID.h>

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;


class TriSurfaceLoad : public Element
{
  public:
    TriSurfaceLoad(int tag, int Nd1, int Nd2, int Nd3,  double pressure, double rhoH_=0); 
    TriSurfaceLoad();
    ~TriSurfaceLoad();

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
    const Matrix &getMass(void);    

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

    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:

  enum {SL_NUM_NODE = 3}; // number of nodes per element
  enum {SL_NUM_NDF = 3}; // d.o.f. per node
  enum {SL_NUM_DOF = 9}; // degrees of freedom per element
  enum {SL_NUM_DDOF = 9}; // displacement degrees of freedom per element
  
    // method to update base vectors g1 & g2
    int UpdateBase(double Xi, double Eta);

    ID  myExternalNodes;      // contains the tags of the end nodes
    static Matrix tangentStiffness;  // Tangent Stiffness matrix
    static Matrix mass;  // mass matrix
    static Matrix damp;  // damping matrix
    static Vector internalForces;    // vector of Internal Forces

    double my_pressure;       // pressure applied to surface of element
    double rhoH;              // A density per unit area to compute a mass matrix (lumped)

    Node *theNodes[SL_NUM_NODE];
    
    Vector g1;		          // tangent vector  = d(x_Xi)/d_xi
    Vector g2;		          // tangent vector  = d(x_Xi)/d_eta 
    Vector myNhat;	          // normal Vector 

    Vector myNI;              // shape functions

    Vector dcrd1;             // current coordinates of node 1
    Vector dcrd2;             // current coordinates of node 2
    Vector dcrd3;             // current coordinates of node 3

    static double oneOverRoot3;
    static double GsPts[1][1];

	double mLoadFactor;       // factor from load pattern
};

#endif




