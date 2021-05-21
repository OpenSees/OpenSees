// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								Kamiar Kalbasi
//								Kutay Orakcal
//								John Wallace
//								 
// User documentation available at: https://kkolozvari.github.io/SFI-MVLEM-3D/
//
// Created: 03/2021
//
// Description: The MVLEM-3D model is a three-dimenaional four-node element with 24 DOFs for nonlinear analysis of 
// flexure-controlled non-rectangular reinforced concrete walls subjected to multidirectional loading. The model is 
// an extension of the two-dimensional, two-node Multiple-Vertical-Line-Element-Model (MVLEM). The baseline MVLEM, 
// which is essentially a line element for rectangular walls subjected to in-plane loading, is extended to a 
// three-dimensional model formulation by: 1) applying geometric transformation of the element in-plane degrees of 
// freedom that convert it into a four-node element formulation and by 2) incorporating linear 
// elastic out-of-plane behavior based on the Kirchhoff plate theory. The in-plane and the out-of-plane 
// element behaviors are uncoupled in the present model formulation.
//
// Notes:
// Nodes should be assigned in counterclockwise direction.
//    4........3 
//    .        .
//    .        .
//    .        . ^ y
//    1........2 |-> x
//
// Reference:
// K. Kolozvari, K. Kalbasi, K. Orakcal & J. W. Wallace (2021), "Three-dimensional shear-flexure interaction model for analysis of non-planar reinforced concrete walls", Journal of Building Engineering.
//
// Source: /usr/local/cvs/OpenSees/SRC/element/mvlem/SFI_MVLEM_3D.cpp
//
// Rev: 1.0


#ifndef SFI_MVLEM_3D_h
#define SFI_MVLEM_3D_h
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
class Node;
class Channel;
class NDMaterial;
class Response;

class SFI_MVLEM_3D : public Element
{
public:

	// Constructors
	SFI_MVLEM_3D(int tag,						// element tag
		double Dens,							// density
		int Nd1, int Nd2, int Nd3, int Nd4,		// end node tags          
		NDMaterial **Materials,					// array of material tags
		double *Thickness,						// array of macro-fiber thickness
		double *Width,							// array of macro-fiber widths
		int mm,									// number of macro-fibers (RC panels)
		double cc,								// center of rotation					
		double nn,								// poisson ratio (for out-of-plane behavior)
		double tf);								// thickness factor (for out-of-plane behavior)

	SFI_MVLEM_3D();

	// Destructor
	~SFI_MVLEM_3D();

	// Public methods to obtain inforrmation about dof & comectivity
	int getNumExternalNodes(void) const;
	const ID &getExternalNodes(void);
	Node **getNodePtrs(void);
	int getNumDOF(void);
	void setDomain(Domain *theDomain);

	// Public methods to set the state of the element    
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);
	int update(void);

	// Public methods to obtain stiffness, mass, damping and residual information    
	const Matrix &getTangentStiff(void);
	const Matrix &getInitialStiff(void);
	const Matrix &getDamp(void);
	const Matrix &getMass(void);

	void zeroLoad(void);
	int addLoad(ElementalLoad *theLoad, double loadFactor);
	int addInertiaLoadToUnbalance(const Vector &accel);
	const Vector &getResistingForce(void);
	const Vector &getResistingForceIncInertia(void);

	// Public methods for output    
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	int displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode);
	void Print(OPS_Stream &s, int flag = 0);
	Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	int getResponse(int responseID, Information &eleInformation);

	Response **theResponses;			// pointer to material responses needed for Concrete

protected:

private:

	void setTransformationMatrix(void);

	// Private member functions - only available to objects of the class
	double *computeCurrentStrain(void);
	double getShearDef(void);
	double getCurvature(void);
	Vector getResistingForce_24DOF(void);
	Vector getResistingForce_24DOF_local(void);

	// Private attributes - a copy for each object of the class
	ID externalNodes;					// contains the id's of end nodes

	// Input variables
	Node *theNodes[4];					// external node pointers +++         
	Node **theNodesX;					// array of internal node pointers
	Node **theNodesALL;					// array of ALL node pointers
	Node *theNd1;						// pointer to bottom node
	Node *theNd2;						// pointer to bottom node
	Node *theNd3;						// pointer to top node
	Node *theNd4;						// pointer to top node
	double density;						// material density
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
	double Eib;							// Modulus of elasticity of imaginary beams
	double Hib;							// Height of imaginary beam
	double Aib;							// Area of imaginary beams cross-section
	double Iib;							// Moment of inertia of imaginary beam

	// Calculated element parameters
	double h;							// height of SFI_MVLEM_3D element (undeformed configuration)
	double Lw;							// lenght of SFI_MVLEM_3D elemtn, i.e. wall length 
	double TotalMass;					// element mass - you dont need this !!!
	double NodeMass;					// nodal mass
	double d;							// distance of corner nodes from center axis
	double K_drilling;					// drilling Stiffness

	// Calculated element arrays
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
	double *SFI_MVLEM_3DStrainX;		// macro-fiber axial strain in X direction (epsX)
	double *SFI_MVLEM_3DStrainY;		// macro-fiber axial strain in Y direction (epsY)
	double *SFI_MVLEM_3DStrainXY;		// macro-fiber shear strain (gammaXY)
	double *SFI_MVLEM_3DStrain;			// macro-fiber strains 
	double *Dens;						// macro-fiber densities

	// For recorders
	double Dsh;							// recorded shear deformation
	Vector P_24DOF;						// recorded forces at end nodal 6DOFs
	Vector P_24DOF_local;

	// Class-wide matrices/vectors/doubles
	// Constant stiffnes terms
	double K1;
	double K2;
	double K3;
	double K4;
	double K5;
	double K6;
	double K7;
	double K8;
	double K9;
	double K10;
	double K11;
	double K12;
	double K13;
	double K14;
	double K15;
	double K16;
	double K17;
	double K18;
	double K19;
	double K20;
	double K21;
	double K22;

	Matrix SFI_MVLEM_3DK;		// stiffness
	Matrix SFI_MVLEM_3DD;		// damping
	Matrix SFI_MVLEM_3DM;		// mass 
	Vector SFI_MVLEM_3DR;		// force

	Matrix SFI_MVLEM_3DKlocal;	// Local stiffness matrix
	Vector SFI_MVLEM_3DRlocal;	// Local force vector
	Matrix SFI_MVLEM_3DMlocal;	// Local mass matrix
	Matrix SFI_MVLEM_3DDlocal;	// Local damping matrix

	Matrix T;							// Transform dofs
	Matrix T6;							// Transform Node DOFs
	Matrix Tt;							// Transform crds
};
#endif
#pragma once