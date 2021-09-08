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
                                                                        
#ifndef LysmerTriangle_h
#define LysmerTriangle_h

// Written: J. Abell
// Created: March 2018
// Modified: 

// Description: This file contains the class definition for LysmerTriangle. 

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>
#include <ID.h>

// number of nodes per element
#define SL_NUM_NODE 3
// d.o.f. per node
#define SL_NUM_NDF  3
// degrees of freedom per element
#define SL_NUM_DOF  9
// displacement degrees of freedom per element
#define SL_NUM_DDOF  3

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;


class LysmerTriangle : public Element
{
  public:
    LysmerTriangle(int tag, int Nd1, int Nd2, int Nd3,  double rho, double Vp, double Vs, double eleLength=0, int stage=0); 
    LysmerTriangle();
    ~LysmerTriangle();

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
    
    // public methods to obtain stiffness, mass, damping and 
    // residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);   
    const Matrix &getDamp(void);     

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

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

  protected:
    
  private:

    // method to update base vectors g1 & g2
    int UpdateBase(double Xi, double Eta);

    Vector internalForces;    // vector of Internal Forces
    Vector springForces;      // vector of spring forces
    ID  myExternalNodes;      // contains the tags of the end nodes

    double rho;               // density
    double Vp;                // compression (P) wave-speed
    double Vs;                // shear (S) wave-speed
    double A;                 // Area of the element
    double element_length;    // Used to assign a penalty stiffness.
    Node *theNodes[SL_NUM_NODE];
    
    Vector g1;                // tangent vector  = d(x_Xi)/d_xi
    Vector g2;                // tangent vector  = d(x_Xi)/d_eta 
    Vector myNhat;            // normal Vector 
    Vector myThat;            // tangent Vector 1
    Vector myShat;            // tangent Vector 2

    Vector myNI;              // shape functions

    Vector dcrd1;             // current coordinates of node 1
    Vector dcrd2;             // current coordinates of node 2
    Vector dcrd3;             // current coordinates of node 3

    Vector gnd_velocity;      //Velocity from ground motions

    int MyTag;                // what is my name?

    static Matrix tangentStiffness;  // Tangent Stiffness matrix
    static Matrix tangentDamping;  // Tangent Stiffness matrix
    static Vector theVector;         // vector to return the residual
    static double oneOverRoot3;
    static double GsPts[1][1];

    static Matrix Bmat;

	double mLoadFactor;       // factor from load pattern
    
    int stage;    // =0 (pure damping) 
                  // =1 (pure stiffness) 
                  // =2 (damping and stiffness) 
                  // =3 (damping but preserve elastic forces from springs after gravity analysis)
};



#endif




