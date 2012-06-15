/////////////////////////////////////////////////////////
//Written: M. Chiaramonte,  P. Arduino,  P.Mackenzie-Helnwein, U.Washington
//			03.29.2011
//Modified: 

// Decription: This file contains the class definition for the Quad4FiberOverlay

#ifndef Quad4FiberOverlay_h
#define Quad4FiberOverlay_h

#include <Node.h>
#include <Element.h>
#include <Vector.h>
#include <UniaxialMaterial.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Channel.h>
#include <G3Globals.h>
#include <ErrorHandler.h>
#include <Matrix.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include <elementAPI.h>
// number of nodes per element
#define SL_NUM_NODE 4
// d.o.f. per node
#define SL_NUM_NDF  2
// degrees of freedom per element
#define SL_NUM_DOF  8
// displacement degrees of freedom per element
#define SL_NUM_DDOF  8

class Domain;
class Node;
class Channel;
class UniaxialMaterial;
class FEM_ObjectBroker;

class Quad4FiberOverlay : public Element
{
	public:
		Quad4FiberOverlay::Quad4FiberOverlay(int tag, int nd1, int nd2, int nd3, int nd4,
UniaxialMaterial &m, double Af,double beta1, double beta2);
		Quad4FiberOverlay();
		~Quad4FiberOverlay();
		
	// // initialization
    void setDomain(Domain *thedomain);
 
    // // public methods to obtain inforrmation about dof & connectivity                 
     int getNumExternalNodes(void) const;
     const ID &getExternalNodes(void);
     Node **getNodePtrs(void);
     int getNumDOF(void);
 
    // // public methods to set the state of the element                                 
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);
 
    // // public methods to obtain stiffness  
    const Matrix &getTangentStiff(void);
	const Matrix &getInitialStiff(void);
 
    // // public method to obtain resisting force
    const Vector &getResistingForce(void);
 
    // // method for obtaining information specific to an element 
    void Print(OPS_Stream &s, int flag =0);
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);
	int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    // public methods for database/parallel processing                                                      
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	
	protected:

	private:
	// private member functions - only available to objects of the class              
    double computeCurrentStrain(void) ;
	int Dual();
	int UpdateBase(double Xi, double Eta);
	int	getEltBb(double Xi, double Eta);
    // private attributes - a copy for each object of the class                       
    UniaxialMaterial *theMaterial;       // pointer to a material                     
    ID  externalNodes;                   // contains the id's of end nodes        
	Node *theNodes[4];
	// static data - single copy for all objects of the class  
    static double pts[4][2];	// Stores quadrature points
    static double wts[2];		// Stores quadrature weights
    static Matrix FiberK;   // class wide matrix for returning stiffness                            
    static Vector P;   // class wide vector for returning residual 
	Vector u;
	Vector g1;							// Basis
	Vector g2;							// Basis
	Vector dualg1;						//Dual Basis
	Vector dualg2;						//Dual Basis

	Vector Q1;							//Current crd 1
	Vector Q2;							//Current crd 2
	Vector Q3;							//Current crd 3
	Vector Q4;							//Current crd 4

	Vector Qfi;							//Fiber i Node
	Vector Qfj;							//Fiber j node
	Vector Vf;							//Fiber Vector

	Vector nFi;
	Vector nFj;

	Matrix dNidxAlphai;					// J^-1 ^ dN/dxi	
	Vector A;							// A vec
	Vector AA;							// A tensor A
	Vector Bb;							// B mat 
	double beta1;						// Location of node for fber
	double beta2;						// Location of node for fiber 
	double Af;							// Area of Fiber
	double Lf;							// Len of fiber
	double strain;						//strain
};

#endif