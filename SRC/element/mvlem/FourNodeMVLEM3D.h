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


#ifndef FourNodeMVLEM3D_h
#define FourNodeMVLEM3D_h
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
class Node;
class Channel;
class UniaxialMaterial;

class FourNodeMVLEM3D : public Element
{
public:

	// constructors
	FourNodeMVLEM3D(int tag,
		double Dens,
		int Nd1, int Nd2, int Nd3, int Nd4,
		Node * theNd1, Node * theNd2, Node * theNd3, Node * theNd4,
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

	FourNodeMVLEM3D();

	// destructor
	~FourNodeMVLEM3D();

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
	int displaySelf(Renderer &theViewer, int displayMode, float fact);
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

	// Imaginary beam parameters
	double Eim;								// modulus of elasticity of imaginary beams
	double Him;								// height of imaginary beam
	double Aim;								// area of imaginary beams cross-section
	double Iim;								// moment of inertia of imaginary beam

	// Calculated element parameters
	double h;								// height of FourNodeMVLEM3D element (undeformed configuration)
	double Lw;								// lenght of FourNodeMVLEM3D elemtn, i.e. wall length 
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
	double *FourNodeMVLEM3DStrain;			// fiber strains (0 to m-1), shear deformation (m)

	// Class-wide matrices/vectors
	// Global cs
	static Matrix FourNodeMVLEM3DK;			// stiffness matrix
	static Matrix FourNodeMVLEM3DD;			// damping matrix
	static Matrix FourNodeMVLEM3DM;			// mass matrix
	static Vector FourNodeMVLEM3DR;			// force vector

	// Local cs
	static Matrix FourNodeMVLEM3DKlocal;	// stiffness matrix
	static Matrix FourNodeMVLEM3DDlocal;	// damping matrix
	static Matrix FourNodeMVLEM3DMlocal;	// mass matrix
	Vector FourNodeMVLEM3DRlocal;			// force vector

	// Transformation matrices
	Matrix T;
	Matrix T6;
	Matrix Tt;

};
#endif
#pragma once