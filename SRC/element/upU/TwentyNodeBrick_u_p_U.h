///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:             TwentyNodeBrick_u_p_U.h
// CLASS:            TwentyNodeBrick_u_p_U
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class for coupled system
// RETURN:
// VERSION:
// LANGUAGE:          C++.ver >= 3.0
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Xiaoyan Wu
// DATE:              Aug. 2001
// UPDATE HISTORY:    Modified from EightNodeBrick_u_p_U.h.   Aug. 2001           
//		      01/16/2002    Xiaoyan
//		      Add the permeability tensor and ks, kf  to the constructor  Xiaoyan 
//
//  "Coupled system" : Solid and fluid coexist.
//                    u-- Solid displacement
//                    p-- Pore pressure
//                    U-- Absolute fluid displacement
//
//
///////////////////////////////////////////////////////////////////////////////



#ifndef TWENTYNODEBRICK_U_P_U_H
#define TWENTYNODEBRICK_U_P_U_H

#ifndef _bool_h
#include "bool.h"
#endif

#include <fstream>	  // add for the output of stiffness matrix K  02/04/2002

#include <Element.h>
#include <Node.h> 
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>

#include <GaussQuadRule1d.h>

// Xiaoyan added from brick3d.hh   07/11/00
#include <basics.h>
#include <BJtensor.h>
#include <nDarray.h>
#include <BJmatrix.h>
#include <BJvector.h>
#include <stresst.h>
#include <straint.h>
#include <MatPoint3D.h>
#include <Template3Dep.h>


class Node;
//class NDMaterial;
//class QuadRule1d;

class TwentyNodeBrick_u_p_U: public Element
{

  public:
   TwentyNodeBrick_u_p_U(int element_number,
                   int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,
                   int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8,
                   int node_numb_9,  int node_numb_10, int node_numb_11, int node_numb_12,
                   int node_numb_13, int node_numb_14, int node_numb_15, int node_numb_16,
                   int node_numb_17, int node_numb_18, int node_numb_19, int node_numb_20,
                   NDMaterial * Globalmmodel, double b1, double b2, double b3,
		   double nn, double alf, double rs,double rf,
		   double permb_x, double permb_y, double permb_z, 
		   double kks, double kkf,double pp);
		   // int dir, double surflevel);
		   //, EPState *InitEPS);   const char * type,


   TwentyNodeBrick_u_p_U ();
   ~TwentyNodeBrick_u_p_U();

    tensor H_3D(double r1, double r2, double r3);
    tensor dH_drst_at(double r1, double r2, double r3);	   // wxy added 08/27/2001
    tensor interp_poli_at(double r1, double r2, double r3);
    tensor dh_drst_at(double r1, double r2, double r3);
//    tensor k_at(double r1, double r2, double r3);	   // wxy added 08/27/2001
//    tensor HU_drst_at(double r1, double r2, double r3);	   // wxy added 08/27/2001
//    tensor dHU_drst_at(double r1, double r2, double r3);   // wxy added 08/27/2001


    tensor getStiffnessTensorKep();
    tensor getStiffnessTensorG1();        //wxy added 08/27/2001
    tensor getStiffnessTensorG2();        //wxy added 08/27/2001
    tensor getStiffnessTensorP();         //wxy added 08/27/2001


    tensor getMassTensorMs();   //wxy added 08/27/2001
    tensor getMassTensorMf();   //wxy added 08/27/2001

    tensor getDampTensorC1();    //wxy added 08/27/2001
    tensor getDampTensorC2();    //wxy added 08/27/2001
    tensor getDampTensorC3();    //wxy added 08/27/2001

    matrix stiffness_matrixKep(const tensor  Kep);
    matrix stiffness_matrixG1(const tensor  G1);
    matrix stiffness_matrixG2(const tensor  G2);
    matrix stiffness_matrixP(const tensor  P);
//    void set_stiffness_MatrixK();

    matrix damping_matrixC1(const tensor C1);
    matrix damping_matrixC2(const tensor C2);
    matrix damping_matrixC3(const tensor C3);
//    void set_damping_MatrixC();

    matrix mass_matrixMs(const tensor  Ms);	   // wxy added 08/27/2001
    matrix mass_matrixMf(const tensor  Mf);	   // wxy added 08/27/2001
//    void  set_mass_MatrixM();              	   // wxy added 08/27/2001

    tensor Jacobian_3D(tensor dh);
    tensor Jacobian_3Dinv(tensor dh);
    tensor Nodal_Coordinates();

    tensor incr_dispDu();
    tensor incr_dispDp();
    tensor incr_dispDU();
    tensor total_dispDu();
    tensor total_dispDU();

   tensor total_dispDu(FILE *fp, double * u);
   tensor total_dispDU(FILE *fp, double * u);

    void incremental_UpdateDu();
    void incremental_UpdateDU();
    void set_strain_stress_tensorDu(FILE *fp, double * u);
    void set_strain_stress_tensorDU(FILE *fp, double * u);
    TwentyNodeBrick_u_p_U & operator[](int subscript);

    // public methods to set the state of the element    
    int commitState ();
    int revertToLastCommit ();
    int revertToStart ();

    // public methods to obtain ...    
    double getrho();		   
    int getNumExternalNodes () const;
    const ID &getExternalNodes ();
    Node **getNodePtrs();

    int getNumDOF ();	
    void setDomain(Domain *theDomain);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff (); 
    const Matrix &getInitialStiff(); 
    const Matrix &getMass (); 
    const Matrix &getConsMassM (); 
    const Matrix &getDamp(); 

//    const Matrix &getConsMassMs (); 
//    const Matrix &getConsMassMf (); 


    void zeroLoad ();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector  FormEquiBodyForce();
    const Vector &getResistingForce ();
    const Vector &getResistingForceIncInertia ();
    int  get_global_number_of_node(int local_node_number);
    int  get_Brick_Number();
    double get_Gauss_p_c(short order, short point_numb);
    double get_Gauss_p_w(short order, short point_numb);
  
    int * get_LM();

    // returns nodal forces for given stress field in an element
    tensor nodal_forcesFu();
    tensor nodal_forcesFU();
    tensor nodal_forces();
    // returns nodal forces for ITERATIVE stress field in an element
     tensor iterative_nodal_forcesFu();
     tensor iterative_nodal_forcesFU();
    // returns nodal forces for given constant stress field in the element
    tensor nodal_forces_from_stress(stresstensor & );
    // returns nodal forces for given incremental strain field in an element
    // by using the linearized constitutive tensor from the begining of the step !
    tensor linearized_nodal_forces();
    // updates Material point stresses and strains from given displacements
    tensor update_stress_strain(tensor & disp);

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


    void report(char *);
    void reportshort(char *);
    void reportPAK(char *);
    void reportpqtheta(int);
    void reportLM(char *);
    void reportTensor(char *);
    void reportCIPIC(char *);
    void reportTensorF(FILE *);

    // Setting initial E according to the initial pressure
    //void setInitE();
    //void reportStressTensorF(FILE *);

  protected:
    
  private:
    // private attributes - a copy for each object of the class

    //NDMaterial * **theMaterial; // pointer to the ND material objects
    
    int numDOF;	    		// Number of element DOF
    ID  connectedExternalNodes; // Tags of quad nodes

    Node *nd1Ptr;		// Pointers to quad nodes
    Node *nd2Ptr;
    Node *nd3Ptr;
    Node *nd4Ptr;
    
    Node *nd5Ptr;
    Node *nd6Ptr;
    Node *nd7Ptr;
    Node *nd8Ptr;

    Node *nd9Ptr;
    Node *nd10Ptr;
    Node *nd11Ptr;
    Node *nd12Ptr;

    Node *nd13Ptr;
    Node *nd14Ptr;
    Node *nd15Ptr;
    Node *nd16Ptr;

    Node *nd17Ptr;
    Node *nd18Ptr;
    Node *nd19Ptr;
    Node *nd20Ptr;


    static Matrix K;		// Element total stiffness Matrix (140*140) (including Kep, G1, G2, P)
//    matrix Kep;		// Element stiffness Matrix (60*60) (Solid part)
//    matrix G1;		// (24*8)
//    matrix G2;		// (24*8)
//    matrix P;		// Pore pressure matrix (8*8)

    static Matrix C;		// The total damping matrix (140*140)    // wxy added 08/27/2001
//    matrix C1;		//           damping matrix (24*24)    // wxy added 08/27/2001
//    matrix C2;		//           damping matrix (24*24)    // wxy added 08/27/2001
//    matrix C3;		//           damping matrix (24*24)    // wxy added 08/27/2001

    static Matrix M;		// Element total mass matrix (140*140) (including Ms and Mf) // wxy added 08/27/2001
//    matrix Ms;		// Element solid mass matrix (24*24) // wxy added 08/27/2001
//    matrix Mf;		// Element fluid mass matrix (24*24)  // wxy added 08/27/2001

    static Vector p;		// Element resisting force vector (140). 
                        // I changed P to p. P for the matrix. Xiaoyan 09/24/2001
    Vector Q;		// Applied nodal loads
    Vector bf;  	// Body forces

    tensor fs;          // Force      // wxy added 09/20/2001
    tensor fp;			      // wxy added 09/20/2001
    tensor ff;			      // wxy added 09/20/2001
    static tensor k;		// wxy added 01/16/2002
    // double thickness;	// Element thickness
    double n;          		     // prrosity                              // wxy added 08/27/2001
    double alpha;		     // coefficient for soil approximate equal 1. 
                                     //  // wxy added 08/27/2001
    double rho;      		     // Total density.
    double rho_s;      		     // Solid density         // wxy added 08/27/2001
    double rho_f;      		     // Fluid density         // wxy added 08/27/2001
    double ks;                       // Bulk modulus of solid. //  wxy added 09/22/2001
    double kf;                       // Bulk modulus of fluid. //  wxy added 09/22/2001
//    double rho=(1-n)*rho_s+n*rho_f;  // Total density       // wxy added 08/27/2001
    double pressure;   		     // Normal surface traction (pressure) over entire element
    int    order;      		     // Order of the quadrature rule
    
    Matrix J;		// Jacobian of transformation
    Matrix L;		// Inverse of J
    Matrix B;		// Strain interpolation matrix
  
    double determinant_of_Jacobian;
    //int  G_N_numbs[8];     // Global node numbers for this element  Xiaoyan changed from 20 to 8
        
    int nodes_in_brick;      // number of nodes   
    
    //Node * nodes;          // pointer to GLOBAL nodes
    
    NDMaterial * mmodel;     // pointer to GLOBAL material models
    
    int r_integration_order; // Gauss-Legendre integration order in r direction
    int s_integration_order; // Gauss-Legendre integration order in s direction
    int t_integration_order; // Gauss-Legendre integration order in t direction
    
    MatPoint3D ** matpoint;  // pointer to array of Material Points
    
    int  LM[60]; // for 20noded x 3 = 60

    Matrix *Ki;
    static Node *theNodes[20];
};


#endif
