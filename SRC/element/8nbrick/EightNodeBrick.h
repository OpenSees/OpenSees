///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              EightNodeBrick.h
// CLASS:             EightNodeBrick
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
// DATE:              Aug. 2000
// UPDATE HISTORY:			 Modified from Brick3D and FourNodeQuad.hh  07/06/00
//																			 Sept. - Oct 2000 connected to OpenSees by Zhaohui
//
// CONTACT:           jeremic@ucdavis.edu
///////////////////////////////////////////////////////////////////////////////
//


#ifndef EightNodeBrick_h
#define EightNodeBrick_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
//#include <Node.h> 

// Commented by Xiaoyan. Use  ~/fem/node.hh  08/04/00
// Released Node.h now. Wu use Opensees's Node.09/27/00


#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>

#include <GaussQuadRule1d.h>

#include <G3Globals.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

// Xiaoyan added from brick3d.hh   07/11/00
#include <basics.h>
#include <BJtensor.h>
#include <nDarray.h>
#include <stresst.h>
#include <straint.h>

//#include <node.h>
//#include <mmodel.h>
#include <fe.h>
#include <gausspnt.h>

//#include <NDMaterial.h>
#include <Template3Dep.h>


class Node;
//class NDMaterial;
//class QuadRule1d;

class EightNodeBrick: public Finite_Element, public Element
{

  public:
    EightNodeBrick(int element_number,
                   int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
                   int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8,
                   NDMaterial * Globalmmodel, const char * type, double b1, double b2,
		   double p, double r, EPState *InitEPS);

    EightNodeBrick ();
    ~EightNodeBrick();

    int Initialize(int element_number,
                   int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
                   int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8,
                   NDMaterial * Globalmmodel, const char * type, double b1, double b2,
                   double p, double r, EPState * InitEPS);
		   
    int getNumExternalNodes () const;
    const ID &getExternalNodes ();
    int getNumDOF ();	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState ();
    int revertToLastCommit ();
    int revertToStart ();

    // public methods to obtain stiffness, mass, damping and residual information    
    // We haven't build the following functions. 
    // All the value of K M Dmp and F are nothing. just
    // want to test the program.  Xiaoyan 08/16/00
    const Matrix &getTangentStiff (); 
    const Matrix &getSecantStiff ();     
    const Matrix &getDamp ();     
    const Matrix &getMass (); 

    void zeroLoad ();
    int addLoad(const Vector &addP);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce ();
    const Vector &getResistingForceIncInertia ();

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf (Renderer &theViewer, int displayMode, float fact);
    void Print(ostream &s, int flag =0);   
    //    Do nothing with void Print (ostream &s, int flag =0);
    //    use Brick3D report.  08/16/00
    int setResponse (char **argv, int argc, Information &eleInformation);
    int getResponse (int responseID, Information &eleInformation);
	
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
    
    // Xiaoyan added 5-8 07/06/00
    Node *nd5Ptr;
    Node *nd6Ptr;
    Node *nd7Ptr;
    Node *nd8Ptr;

    Matrix K;		// Element stiffness Matrix
    Matrix C;		// Element damping matrix
    Matrix M;		// Element mass matrix
    Vector P;		// Element resisting force vector
    Vector Q;		// Applied nodal loads
    Vector b;		// Body forces

    // K, M and P are commented by Xiaoyan . We use the K, M and P from Brick3D

    // double thickness;	// Element thickness
    double rho;		// Mass per unit volume
    double pressure;	// Normal surface traction (pressure) over entire element
    int order;		// Order of the quadrature rule

    //QuadRule1d *theQuadRule;	// Integration rule

    Matrix J;		// Jacobian of transformation
    Matrix L;		// Inverse of J
    Matrix B;		// Strain interpolation matrix
    
	
    //    // static data - single copy for all objects of the class 
    //    static G3Matrix N;	// Displacement interpolation matrix
    	
    //    // private member functions - only objects of this class can call these
    //    void setJacobian (double r, double s, double t);         // Xiaoyan changed 
    //    double formDetJ (double r, double s, double t);	   // xi, eta to r and s
    //    void formBMatrix (double r, double s, double t);	   // and added t
    //    static void formNMatrix (double r, double s, double t);  // 07/06/00
    
    // Commented by Xiaoyan. We use Brick3D for calculating these.


  //  The following is taken from brick3d.hh. Xiaoyan 07/11
  private:
    // element number (tag)
    //unsigned int  elem_numb;      
    
    double determinant_of_Jacobian;
    //int  G_N_numbs[8];  // Global node numbers for this element  Xiaoyan changed from 20 to 8
        
    int nodes_in_brick;     // number of nodes ( from 8-20 //8 now Zhaohui)  
    
    //Node * nodes;                 // pointer to GLOBAL nodes
    
    NDMaterial * mmodel;          // pointer to GLOBAL material models
    
    int r_integration_order; // Gauss-Legendre integration order in r direction
    int s_integration_order; // Gauss-Legendre integration order in s direction
    int t_integration_order; // Gauss-Legendre integration order in t direction
    
    // Now I want 3D array of Gauss points!
    // MatPoint[r_integration_order][s_integration_order][t_integration_order]
    // 3D array of Gauss points
    IntegrationPoint * MatPoint;    // pointer to array of Gauss Points
    
    // 3D array of material models for each Gauss points
    //NDMaterial *GPmmodel;  // pointer to array of material models for Gauss Points
    // Do we need this one? 
    
    //..NDMaterial  *MatPoint;  // Zhaohui  10-01-2000
    
    
    // this tensor is tangent constitutive tensor MS definition.
    // It represents the state at Gauss point before applying strain increment
    // that was produced by incremental displacements
    
    // Zhaohui  10-01-2000
    //..tensor * GPtangent_E;  // pointer to array of constitutive tensors for Gauss Points
    
    // this stress tensor is start_stress from my MS definition.
    // It represents the state at Gauss point before applying strain increment
    // that was produced by incremental displacements
    
    // Zhaohui  10-01-2000
    //..stresstensor * GPstress;  // pointer to array of stresstensors for Gauss Points
    
    // 3D array of stresstensors models for each Gauss points
    // this stress tensor is iterative stress.
    // It represents the state at Gauss point during iterative procedure on FEM level.
    // Zhaohui  10-01-2000
    //stresstensor * GPiterative_stress;  // pointer to array of stresstensors for Gauss Points
    //double * GPq_ast_iterative;  // pointer to array of iterative values of internal variable
    
    // 3D array of straintensors models for each Gauss points
    // this strain tensor is strain_increment from my MS definition.
    // It represents the additional strains that are to be
    // integrated. After numerical integration is done, the return value
    // ( from one of the numerical integration procedures )
    // is stresstensor that should then be put in GPstress place!
    // Zhaohui  10-01-2000
    //straintensor * GPstrain;  // pointer to array of straintensors for Gauss Points
    
    // this is LM array. This array holds DOFs for this element
    int  LM[24]; // for 8noded x 3 = 24
  public:
    //CONSTRUCTOR
    //EightNodeBrick(int tag = 0, // default constructor
    //               //int r_int_order = 0, //int s_int_order = 0, //int t_int_order = 0,
    //               int node_numb_1  = 0,
    //               int node_numb_2  = 0,
    //               int node_numb_3  = 0,
    //               int node_numb_4  = 0,
    //               int node_numb_5  = 0,
    //               int node_numb_6  = 0,
    //               int node_numb_7  = 0,
    //               int node_numb_8  = 0,
    //               NDMaterial * Globalmmodel = 0,
    //               //Node       * GlobalNodes = 0,   // maybe put NULL
    //               //tensor   * IN_tangent_E = 0,         //stresstensor * INstress = 0,                //stresstensor * INiterative_stress = 0,                   //double       * INq_ast_iterative =0,                   //straintensor * INstrain = 0
    //               EPState *eps =0
    //               );

    //void Initialize(int element_numb, // Initialize function
    //                //short int r_int_order,
    //                //short int s_int_order,
    //                //short int t_int_order,
    //                int node_numb_1,
    //                int node_numb_2,
    //                int node_numb_3,
    //                int node_numb_4,
    //                int node_numb_5,
    //                int node_numb_6,
    //                int node_numb_7,
    //                int node_numb_8,
    //                //Node * GlobalNodes,
    //		    NDMaterial * Globalmmodel,
    //                EPState *eps );
    //		    //tensor   * IN_tangent_E,//stresstensor * INstress, //stresstensor * INiterative_stress, //double       * INq_ast_iterative, //straintensor * INstrain
                    
    
    void incremental_Update(void);
    //void iterative_Update(void);

    tensor H_3D(double r1, double r2, double r3);
    tensor interp_poli_at(double r, double s, double t);
    tensor dh_drst_at(double r, double s, double t);


    //CE Dynamic Allocation for for brick3d s.
    //Finite_Element * new_el( int total );
    Finite_Element & operator[](int subscript);
    //Finite_Element & operator[](int subscript);
    //Finite_Element & operator[](int subscript);
   
    tensor getStiffnessTensor(void);
    //matrix stiffness_tensor(void); 

    void set_strain_stress_tensor(FILE *fp, double * u);
    tensor getMassTensor(void);
    double Potential_Energy(void);

    tensor Jacobian_3D(tensor dh);
    tensor Jacobian_3Dinv(tensor dh);
    tensor Nodal_Coordinates(void);

    tensor incr_disp(void);
    tensor total_disp(FILE *fp, double * u);

    tensor stiffness_matrix(const tensor & K);
    tensor mass_matrix(const tensor & M);


    int  get_global_number_of_node(int local_node_number);
    int  get_Brick_Number(void);


    int * get_LM(void);
    //void set_LM(Node * node); // commented out temporarily 09-27-2000 Zhaohui



    // returns nodal forces for given stress field in an element
    tensor nodal_forces(void);
    // returns nodal forces for ITERATIVE stress field in an element
    tensor iterative_nodal_forces(void);
    // returns nodal forces for given constant stress field in the element
    tensor nodal_forces_from_stress(stresstensor & );
    // returns nodal forces for given incremental strain field in an element
    // by using the linearized constitutive tensor from the begining of the step !
    tensor linearized_nodal_forces(void);
    // updates Gauss point stresses and strains from given displacements
    tensor update_stress_strain(tensor & disp);

    void report(char *);
    void reportshort(char *);
    void reportPAK(char *);
    void reportpqtheta(int);
    void reportLM(char *);
    void reportTensor(char *);
    void reportCIPIC(char *);
    void reportTensorF(FILE *);

};


#endif

