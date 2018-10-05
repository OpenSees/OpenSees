
// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// Sensitivity analysis of contact element for fluid-structure coupling
//------------------------------------------

#ifndef ASI3D8QuadWithSensitivity_H
#define ASI3D8QuadWithSensitivity_H

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Node.h>

#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>

#include <OPS_Globals.h>

#include <Matrix.h>
#include <Vector.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>

#include <ElementalLoad.h>

#define ELE_TAG_ASI3D8QuadWithSensitivity 100002

class Node;


class ASI3D8QuadWithSensitivity: public Element
{

  public:
    ASI3D8QuadWithSensitivity(int element_number,
      int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4, 
      int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8);

    ASI3D8QuadWithSensitivity ();
    ~ASI3D8QuadWithSensitivity();

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
    const Matrix &getMass(void);
    const Matrix &getConsMass(void);

    void zeroLoad(void);
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker
		  &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);

    void Print(OPS_Stream &s, int flag =0);
    Response *setResponse (const char **argv, int argc, OPS_Stream &theHandler);
    int getResponse (int responseID, Information &eleInformation);

	// AddingSensitivity:BEGIN //////////////////////////////////////////
    int            setParameter                (const char **argv, int argc, Parameter &param);
    int            updateParameter             (int parameterID, Information &info);
	int            activateParameter           (int parameterID);
	const Vector & getResistingForceSensitivity(int gradNumber);
	const Matrix & getKiSensitivity            (int gradNumber);
	const Matrix & getMassSensitivity          (int gradNumber);
	int            commitSensitivity           (int gradNumber, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////


    
    // const Vector &getExternalLoadIncInertia(void);
    
    double get_Gauss_p_c(short order, short point_numb);
    double get_Gauss_p_w(short order, short point_numb);

  protected:

  private:
    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes; // Tags of quad nodes

    Matrix *Ki;
    Node *theNodes[8];
    
    // Matrix **L;       // global differential operator
    // double *detJ;     // determinant of Jacobian matrix
    
    NDMaterial **theMaterial; // pointer to the ND material objects
    
    // Vector Q;		           // Applied nodal loads
    
    unsigned char hasConstrained; // 0:not set, 1: constrained, 2: not constrained

  private:
    static Matrix K;		   // Element stiffness Matrix
    static Matrix C;		   // Element damping matrix
    static Matrix M;		   // Element mass matrix
    static Vector P;		   // Element resisting force vector
    static Matrix QMAT;		 // [Q] matrix
	    static Matrix mass ; 

    
    static ID integFlags;  // integrator flags
    static ID actDOFs;     // activated element dofs, add Yichao Gao
    
    static const int numDOF;               // DOF number of element
    static const int nodes_in_elem;        // number of nodes in element
    static const int nodes_in_quad;        // number of nodes in quad
    static const int r_integration_order;  // Gauss-Legendre integration order in r direction
    static const int s_integration_order;  // Gauss-Legendre integration order in s direction
    static const int dim;                  // spatial dimension
    static const int numGP;                // number of Gauss point
    static const int numSDOF;              // number of structure DOF
    static const int numFDOF;              // number of fluid DOF
    
    static Matrix **H;  // Matrix array holds h
    static Matrix **DH;  // Matrix array holds h
    // static Matrix **HH;  // Matrix array holds h*h
    	// AddingSensitivity:BEGIN ///////////////////////////
	int parameterID;
// AddingSensitivity:END /////////////////////////////
     int setNDMaterial(NDMaterial *Globalmmodel);
   //material information
    NDMaterial *materialPointers[4] ; //pointers to eight materials
    static Vector VecS;  // Vector(numSDOF);
    static Vector VecF;  // Vector(numFDOF);
    

  public:
    // shape function
    Matrix interp_fun(double r, double s);
    Matrix diff_interp_fun(double r, double s);
    
    // compute Jacobian matrix and its determinant
    // Matrix Jacobian(Matrix dh);
    // Matrix Jacobian(Matrix dh, Matrix h);
    // double Jacobian_det(Matrix Jac);
    
    int computeH(void);
    // int computeHH(void);
    // int computeDiff(void);
    
    // get nodal coordinates
    Matrix getNodalCoords(void);
    
    // get nodal forces from stress
    // Matrix getNodalForces(void);
    
    // get total displacement
    // Matrix getTotalDisp(void);
    
    // get [Q] Matrix
    Matrix &getQMatrix(void);
    
    // ///////////////////////////////////////////////////////////////
    
    ID *getActiveDofs(void);
    int getIntegrateFlag(void);
    ID *getIntegrateFlags(void);
    
    // int setNDMaterial(NDMaterial *Globalmmodel);
    
    ASI3D8QuadWithSensitivity & operator[](int subscript);
    
    

};


#endif

