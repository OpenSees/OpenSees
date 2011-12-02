///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              TwentyNodeBrick.h
// CLASS:             TwentyNodeBrick
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class
// RETURN:
// VERSION:
// LANGUAGE:          C++.ver >= 3.0
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Zhaohui Yang and Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Zhaohui Yang  and Xiaoyan Wu
// DATE:              Aug. 2001
// UPDATE HISTORY:			 
//
///////////////////////////////////////////////////////////////////////////////
//



#ifndef TWENTYNODEBRICK_H
#define TWENTYNODEBRICK_H

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Node.h> 


#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>

#include <GaussQuadRule1d.h>

#include <OPS_Globals.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <basics.h>
#include <BJtensor.h>
#include <nDarray.h>
#include <stresst.h>
#include <straint.h>

#include <MatPoint3D.h>

#include <Template3Dep.h>

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

class Node;

class TwentyNodeBrick: public Element
{

  public:
    TwentyNodeBrick(int element_number,
                   int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,
                   int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8,
                   int node_numb_9,  int node_numb_10, int node_numb_11, int node_numb_12,
                   int node_numb_13, int node_numb_14, int node_numb_15, int node_numb_16,
                   int node_numb_17, int node_numb_18, int node_numb_19, int node_numb_20,
                   NDMaterial * Globalmmodel, double b1, double b2, double b3,
		   double r, double p);

    TwentyNodeBrick ();
    ~TwentyNodeBrick();

    int getNumExternalNodes () const;
    const ID &getExternalNodes ();
    Node **getNodePtrs();

    int getNumDOF ();	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState ();
    int revertToLastCommit ();
    int revertToStart ();

    // update, Guanzhou added Apr. 2004 to update incremental strain in the domain
    int update(void);
    
    // public methods to obtain stiffness, mass, damping and residual information    
    // We haven't build the following functions. 
    // All the value of K M Dmp and F are nothing.
    const Matrix &getTangentStiff (); 
    const Matrix &getInitialStiff(); 
    const Matrix &getMass (); 

    const Matrix &getConsMass (); 
    const Matrix &getLumpedMass (); 

    void zeroLoad ();
    int addLoad(ElementalLoad *theLoad, double loadFactor);    
    //int addLoad(const Vector &addP);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector  FormEquiBodyForce(void);
    const Vector &getResistingForce ();
    const Vector &getResistingForceIncInertia ();

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf (Renderer &theViewer, int displayMode, float fact);
    void Print(OPS_Stream &s, int flag =0);   
    //    Do nothing with void Print (OPS_Stream &s, int flag =0);
    //    use Brick3D report.  08/16/00
    Response *setResponse (const char **argv, int argc, Information &eleInformation);
    int getResponse (int responseID, Information &eleInformation);
	
  protected:
    
  private:
    // private attributes - a copy for each object of the class

    //NDMaterial * **theMaterial; // pointer to the ND material objects
    
    int numDOF;	    		// Number of element DOF
    ID  connectedExternalNodes; // Tags of quad nodes

    Matrix *Ki;
    Node *theNodes[20];

    static Matrix K;    // Element stiffness Matrix
    static Matrix C;    // Element damping matrix
    static Matrix M;    // Element mass matrix
    static Vector P;    // Element resisting force vector
    Vector Q;		// Applied nodal loads
    Vector bf;  	// Body forces
    
    // double thickness;	// Element thickness
    double rho;		// Mass per unit volume DO WE GET THIS ONE OUT!!!
    double pressure;	// Normal surface traction (pressure) over entire element
    int    order;  	// Order of the quadrature rule

    //Matrix J;		// Jacobian of transformation
    //Matrix L;		// Inverse of J
    //Matrix B;		// Strain interpolation matrix
    
	
  private:
    
    double determinant_of_Jacobian;
    int nodes_in_brick;      // number of nodes ( 20 now Zhaohui)  

    NDMaterial * mmodel;     // pointer to GLOBAL material models
    
    int r_integration_order; // Gauss-Legendre integration order in r direction
    int s_integration_order; // Gauss-Legendre integration order in s direction
    int t_integration_order; // Gauss-Legendre integration order in t direction
    
    // Now I want 3D array of Material points!
    // MatPoint3D[r_integration_order][s_integration_order][t_integration_order]
    // 3D array of Material points
    MatPoint3D ** matpoint;  // pointer to array of Material Points
    
    // this is LM array. This array holds DOFs for this element
    //int  LM[60]; // for 20noded x 3 = 60
  public:
    
    void incremental_Update(void);
    //void iterative_Update(void);

    tensor H_3D(double r1, double r2, double r3);
    tensor interp_poli_at(double r, double s, double t);
    tensor dh_drst_at(double r, double s, double t);


    TwentyNodeBrick & operator[](int subscript);
   
    tensor getStiffnessTensor(void);

    void set_strain_stress_tensor(FILE *fp, double * u);
    tensor getMassTensor(void);

    tensor Jacobian_3D(tensor dh);
    tensor Jacobian_3Dinv(tensor dh);
    tensor Nodal_Coordinates(void);

    tensor incr_disp(void);
    tensor total_disp(void);

    tensor total_disp(FILE *fp, double * u);

    tensor stiffness_matrix(const tensor & K);
    tensor mass_matrix(const tensor & M);


    int  get_global_number_of_node(int local_node_number);
    int  get_Brick_Number(void);


    //int * get_LM(void);
    //void set_LM(Node * node); // commented out temporarily 09-27-2000 Zhaohui

    //these two files are originally in fe.h
    double get_Gauss_p_c(short order, short point_numb);
    double get_Gauss_p_w(short order, short point_numb);

    // returns nodal forces for given stress field in an element
    tensor nodal_forces(void);
    // returns nodal forces for ITERATIVE stress field in an element
    tensor iterative_nodal_forces(void);
    // returns nodal forces for given constant stress field in the element
    tensor nodal_forces_from_stress(stresstensor & );
    // returns nodal forces for given incremental strain field in an element
    // by using the linearized constitutive tensor from the begining of the step !
    tensor linearized_nodal_forces(void);
    // updates Material point stresses and strains from given displacements
    tensor update_stress_strain(tensor & disp);

    void report(char *);
    void reportshort(char *);
    void reportPAK(char *);
    void reportpqtheta(int);
    //void reportLM(char *);
    void computeGaussPoint(void);
    void reportCIPIC(char *);
    void reportTensorF(FILE *);
    Vector getWeightofGP(void);


    // Setting initial E according to the initial pressure
    //void setInitE(void);
    //void reportStressTensorF(FILE *);

};


#endif

