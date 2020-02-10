// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								Kamiar Kalbasi
//								California State University, Fullerton 
//
// Created:
//
// Description: 
//
// References:
//
// Source: /usr/local/cvs/OpenSees/SRC/element/FourNodeMVLEM3D/FourNodeMVLEM3D.h
//
// Rev: 1


#ifndef FourNodeSFI_MVLEM3D_h
#define FourNodeSFI_MVLEM3D_h
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
class Node;
class Channel;
class NDMaterial;
class Response;

class FourNodeSFI_MVLEM3D : public Element
{
public:

	// constructors
	FourNodeSFI_MVLEM3D(int tag,		// element tag
		int Nd1, int Nd2, int Nd3, int Nd4,				// end node tags          
		Node *theNd1, Node *theNd2, Node * theNd3, Node * theNd4,		// pointers to end nodes (used to obtain end node locations in constructor)
		NDMaterial **Materials,			// array of material tags
		double *Thickness,				// array of macro-fiber thickness
		double *Width,					// array of macro-fiber widths
		Domain *theTclDomain,			// TclDomain - for adding internal nodes to the element
		int mm,							// number of macro-fibers (RC panels)
		double cc,						// center of rotation					
		double nn,						// poisson ratio (for out-of-plane behavior)
		double tf);						// thickness factor (for out-of-plane behavior)

	FourNodeSFI_MVLEM3D();

	// destructor
	~FourNodeSFI_MVLEM3D();

	// public methods to obtain inforrmation about dof & comectivity
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
	const Matrix &getDamp(void);
	const Matrix &getMass(void);

	void zeroLoad(void);
	int addLoad(ElementalLoad *theLoad, double loadFactor);
	int addInertiaLoadToUnbalance(const Vector &accel);
	const Vector &getResistingForce(void);
	const Vector &getResistingForceIncInertia(void);

	// public methods for output    
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	int displaySelf(Renderer &theViewer, int displayMode, float fact);
	void Print(OPS_Stream &s, int flag = 0);
	Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	int getResponse(int responseID, Information &eleInformation);

	Response **theResponses;			// pointer to material responses needed for Concrete

protected:

private:

	void setTransformationMatrix(void);

	// private member functions - only available to objects of the class
	double *computeCurrentStrain(void);
	double getShearDef(void);
	double getCurvature(void);
	Vector getResistingForce_24DOF(void);
	Vector getResistingForce_24DOF_local(void);

	// private attributes - a copy for each object of the class
	ID externalNodes;					// contains the id's of end nodes

	// input variables
	Node *theNodes[4];					// external node pointers +++         
	Node **theNodesX;					// array of internal node pointers
	Node **theNodesALL;					// array of ALL node pointers
	Node *theNd1;						// pointer to bottom node
	Node *theNd2;						// pointer to bottom node
	Node *theNd3;						// pointer to top node
	Node *theNd4;						// pointer to top node
	NDMaterial **theMaterial;			// array of ND materials
	Vector *theLoad;					// pointer to element load
	double c;							// center of rotation
	int m;								// no. of RC panels 
	double NUelastic;					// Poisson ratio
	double Tfactor;						// Thickness factor

	// Nodal coordinates
	Vector nd1Crds;
	Vector nd2Crds;
	Vector nd3Crds;
	Vector nd4Crds;

	// Out of plane parameters
	double Eave;						// Modulus of Elasticity
	double Tave;						// Thickness
	double *modifiedT;					// Modified Thickness 

	// Imaginary Beam Parameters
	double Eim;							// Modulus of elasticity of imaginary beams
	double Him;							// Height of imaginary beam
	double Aim;							// Area of imaginary beams cross-section
	double Iim;							// Moment of inertia of imaginary beam

	// calculated element parameters
	double h;							// height of FourNodeSFI_MVLEM3D element (undeformed configuration)
	double Lw;							// lenght of FourNodeSFI_MVLEM3D elemtn, i.e. wall length 
	double TotalMass;					// element mass - you dont need this !!!
	double NodeMass;					// nodal mass
	double d;							// distance of corner nodes from center axis
	double K_drilling;					// drilling Stiffness

	// caldulated element arrays
	double *x;							// macro-fiber locations
	double *b;							// macro-fiber widths
	double *t;							// macro-fiber thickness
	double *AcX;						// macro-fiber areas in X direction
	double *AcY;						// macro-fiber areas in Y direction
	double *kx;							// macro-fiber axial stiffness in X direction
	double *ky;							// macro-fiber axial stiffness in Y direction
	double Kh;							// macro-fiber shear stiffness
	double *Fx;							// macro-fiber axial force in X direction
	double *Fy;							// macro-fiber axial force in Y direction
	double *Fxy;						// macro-fiber shear force in (in horizontal plane)
	double *Dx;							// macro-fiber axial deformation in X direction					
	double *Dy;							// macro-fiber axial deformation in Y direction
	double *Dxy;						// macro-fiber shear deformation (in horizontal plane)
	double *FourNodeSFI_MVLEM3DStrainX;	// macro-fiber axial strain in X direction (epsX)
	double *FourNodeSFI_MVLEM3DStrainY;	// macro-fiber axial strain in Y direction (epsY)
	double *FourNodeSFI_MVLEM3DStrainXY;// macro-fiber shear strain (gammaXY)
	double *FourNodeSFI_MVLEM3DStrain;	// macro-fiber strains 
	double *Dens;						// macro-fiber densities

	// for recorders
	double Dsh;							// recorded shear deformation
	Vector P_24DOF;						// recorded forces at end nodal 6DOFs
	Vector P_24DOF_local;

	// class wide matrices
	Matrix FourNodeSFI_MVLEM3DK;		// stiffness
	Matrix FourNodeSFI_MVLEM3DD;		// damping
	Matrix FourNodeSFI_MVLEM3DM;		// mass 
	Vector FourNodeSFI_MVLEM3DR;		// force

	Matrix FourNodeSFI_MVLEM3DKlocal;	// Local stiffness matrix
	Vector FourNodeSFI_MVLEM3DRlocal;	// Local force vector
	Matrix FourNodeSFI_MVLEM3DMlocal;	// Local mass matrix
	Matrix FourNodeSFI_MVLEM3DDlocal;	// Local damping matrix

	Matrix T;							// Transform dofs
	Matrix T6;							// Transform Node DOFs
	Matrix Tt;							// Transform crds
};
#endif
#pragma once