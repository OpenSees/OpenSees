// Code written/implemented by:	Carlos López Olea (carlos.lopez.o@ug.uchile.cl)
//								Leonardo M. Massone
//								Kristijan Kolozvari
//
// User documentation available at: https://github.com/carloslopezolea/E-SFI_Documentation
//
// Created: 06/2022
//
// Description: The Efficient Shear-Flexure Interaction (E-SFI) element was developed based on the SFI-MVLEM formulation. 
// The E-SFI element incorporates the shear-flexure interaction phenomenon by replacing the m number of uniaxial fibers of 
// the MVLEM, by two-dimensional RC panel elements subjected to membrane actions (FSAM). An E-SFI element is described by 
// six degrees of freedom, and therefore no additional degrees of freedom are incorporated into the original MVLEM formulation, 
// as in the SFI-MVLEM. The curvature of an E-SFI element is assumed to be uniform, and the resultant rotation is concentrated
// at height ch. The kinematic assumption of plane sections remain plane, as well as the assumption of constant shear strain 
// along the element length, are considered for computing the axial and shear strains for each panel over the entire section. 
// To complete the strain field of a panel element, a calibrated expression for the horizontal normal strain is implemented to 
// obtain accurate predictions from squat to slender RC walls.
//
// References:
// 1.- Massone, L. M., López, C. N., & Kolozvari, K. (2021). Formulation of an efficient shear-flexure interaction model for planar reinforced concrete walls. Engineering Structures, 243, 112680.
// 2.- López, C. N., Massone, L. M., & Kolozvari, K. (2022). Validation of an efficient shear-flexure interaction model for planar reinforced concrete walls. Engineering Structures, 252, 113590.
// 3.- López C. N. Efficient shear-flexure interaction model for nonlinear analysis of reinforced concrete structural walls. MS Dissertation. Santiago, Chile: University of Chile; 2021.
// 
// Source: /usr/local/cvs/OpenSees/SRC/element/mvlem/E_SFI.cpp
//
// Rev: 1.0


#ifndef E_SFI_h
#define E_SFI_h
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
class Node;
class Channel;
class NDMaterial;
class Response;

class E_SFI : public Element
{
public:

	// constructors
	E_SFI(int tag,					    // element tag
		int Nd1, int Nd2,				// end node tags          
		NDMaterial **Materials,			// array of material tags
		double *Thickness,				// array of macro-fiber thickness
		double *Width,					// array of macro-fiber widths
		int mm,							// number of macro-fibers (RC panels)
		double cc);						// center of rotation					

	E_SFI();

	// destructor
	~E_SFI();

	// public methods to obtain information about dof & connectivity
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
	int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);

	void Print(OPS_Stream &s, int flag =0);    
	Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	int getResponse(int responseID, Information &eleInformation);

protected:

private:

	// private member functions - only available to objects of the class
	void computeCurrentStrain(void);
	double getShearDef(void);
	double getCurvature(void);
	Vector getResistingForce_6DOF(void);

	// private attributes - a copy for each object of the class

	// input variables
	Node *theNodes[2];					// external node pointers          
	Node **theNodesALL;					// array of ALL node pointers
	Node *theNd1;						// pointer to bottom node
	Node *theNd2;						// pointer to top node
	NDMaterial **theMaterial;			// array of ND materials
	Vector *theLoad;					// pointer to element load
	const double c;						// center of rotation
	const int m;						// no. of RC panels

	// calculated element parameters
	ID externalNodes;					// contains the id's of end nodes
	Matrix trans;						// hold the transformation matrix, could use a Vector
	double h;							// height of E_SFI element (undeformed configuration)
	double Lw;							// length of E_SFI element, i.e. wall length
	double TotalMass;					// element mass - you don't need this !!!
	double NodeMass;					// nodal mass

	// caldulated element arrays
	double *x;							// macro-fiber locations
	double *b;							// macro-fiber widths
	double *t;							// macro-fiber thickness
	double *AcX;						// macro-fiber areas in X direction
	double *AcY;						// macro-fiber areas in Y direction
	double *kx;							// macro-fiber axial stiffness in X direction
	double *ky;							// macro-fiber axial stiffness in Y direction
	double *kh;							// macro-fiber shear stiffness
	double *Fx;							// macro-fiber axial force in X direction
	double *Fy;							// macro-fiber axial force in Y direction
	double *Fxy;						// macro-fiber shear force in (in horizontal plane)
	double *Dx;							// macro-fiber axial deformation in X direction					
	double *Dy;							// macro-fiber axial deformation in Y direction
	double *Dxy;						// macro-fiber shear deformation (in horizontal plane)
	double *E_SFIStrainX;				// macro-fiber axial strain in X direction (epsX)
	double *E_SFIStrainY;				// macro-fiber axial strain in Y direction (epsY)
	double *E_SFIStrainXY;				// macro-fiber shear strain (gammaXY)
	double *E_SFIStrain;				// macro-fiber strains 
	double *Dens;						// macro-fiber densities

	// for recorders
	double Dsh;							// recorded shear deformation
	Vector P_6DOF;						// recorded forces at end nodal 6DOFs

	// class wide matrices
	Matrix E_SFIK;						// stiffness
	Matrix E_SFID;						// damping
	Matrix E_SFIM;						// mass 
	Vector E_SFIR;						// force

};
#endif
