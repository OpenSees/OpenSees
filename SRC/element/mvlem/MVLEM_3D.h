// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								Kamiar Kalbasi
//								Kutay Orakcal
//								John Wallace
//								 
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
// K. Kolozvari, K. Kalbasi, K. Orakcal & J. W. Wallace (2021), "Three-dimensional model for nonlinear analysis of slender flanged reinforced concrete walls", Engineering Structures.
//
// Source: /usr/local/cvs/OpenSees/SRC/element/mvlem/MVLEM_3D.cpp
//
// Rev: 1.0


#ifndef MVLEM_3D_h
#define MVLEM_3D_h
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
class Node;
class Channel;
class UniaxialMaterial;

class MVLEM_3D : public Element
{
public:

	// constructors
	MVLEM_3D(int tag,
		double Dens,
		int Nd1, int Nd2, int Nd3, int Nd4,
		UniaxialMaterial **materialsConcrete,
		UniaxialMaterial **materialsSteel,
		UniaxialMaterial **materialsShear,
		double *Rho,
		double *thickness,
		double *width,
		int mm,
		double cc,
		double nn,
		double tf);

	MVLEM_3D();

	// destructor
	~MVLEM_3D();

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

	void setTransformationMatrix(void);

	// public methods for output    
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	int displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode);
	void Print(OPS_Stream &s, int flag = 0);
	Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	int getResponse(int responseID, Information &eleInformation);

protected:

private:

	// private member functions - only available to objects of the class
	double * computeCurrentStrain(void);
	Vector getResistingForceLocal(void);
	double getCurvature(void);
	Vector getStrain(void);
	Vector getStressConcrete(void);
	Vector getStressSteel(void);
	Vector getShearFD(void);
	double getShearDef(void);

	// private attributes - a copy for each object of the class
	ID  externalNodes;          			// contains the id's of end nodes

	// Input variables
	Node *theNodes[4];						// external node pointers       
	Node *theNd1;							// pointer to bottom left node
	Node *theNd2;							// pointer to bottom right node
	Node *theNd3;							// pointer to top left node
	Node *theNd4;							// pointer to top right node
	UniaxialMaterial **theMaterialsConcrete; // pointers to Concrete uniaxial material
	UniaxialMaterial **theMaterialsSteel;	// pointers to Steel uniaxial material
	UniaxialMaterial **theMaterialsShear;	// pointers to Shear uniaxial material
	double density;							// material density
	double c;								// center of rotation
	int m;									// no. of RC panels
	Vector *theLoad;						// pointer to element load

	// Node coordinates
	Vector nd1Crds;
	Vector nd2Crds;
	Vector nd3Crds;
	Vector nd4Crds;

	// Out of plane parameters 
	double NUelastic;						// Poisson ratio
	double Eave;							// average modulus of elasticity
	double Tave;							// average thickness
	double Tfactor;							// effective thickness factor

	// Internal beam parameters
	double Eib;								// modulus of elasticity of imaginary beams
	double Hib;								// height of imaginary beam
	double Aib;								// area of imaginary beams cross-section
	double Iib;								// moment of inertia of imaginary beam

	// Calculated element parameters
	double h;								// height of MVLEM_3D element (undeformed configuration)
	double Lw;								// lenght of MVLEM_3D elemtn, i.e. wall length 
	double A;								// wall cross-section area
	double NodeMass;						// nodal mass
	double d;								// distance of corner nodes from center axis

	// Calculated element arrays
	double *x;
	double *b;								// fiber widths
	double *t;								// fiber thickness
	double *rho;							// fiber reinforcing ratio
	double *Ac;								// concrete area
	double *As;								// steel area
	double *Ec;								// concrete tangent modulus
	double *Es;								// steel tangent modulus	
	double *stressC;						// concrete stress
	double *stressS;						// steel stress
	double *ky;								// fiber axial stiffness in Y direction
	double *kh;								// element shear stiffness		
	double *MVLEM_3DStrain;					// fiber strains (0 to m-1), shear deformation (m)

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

	// Global cs
	static Matrix MVLEM_3DK;				// stiffness matrix
	static Matrix MVLEM_3DD;				// damping matrix
	static Matrix MVLEM_3DM;				// mass matrix
	static Vector MVLEM_3DR;				// force vector

	// Local cs
	static Matrix MVLEM_3DKlocal;			// stiffness matrix
	static Matrix MVLEM_3DDlocal;			// damping matrix
	static Matrix MVLEM_3DMlocal;			// mass matrix
	static Vector MVLEM_3DRlocal;			// force vector

	// Transformation matrices
	Matrix T;
	Matrix T6;
	Matrix Tt;

};
#endif
#pragma once