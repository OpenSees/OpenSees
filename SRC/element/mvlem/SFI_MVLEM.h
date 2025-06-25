// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								California State University, Fullerton 
//								Kutay Orakcal
//								Bogazici University, Istanbul, Turkey
//								John Wallace
//								University of California, Los Angeles
//
// Created: 07/2015
//
// Description: This file contains the class definition for Shear-Flexure
// Interaction Multiple Vertical Line Element Model - SFI_MVLEM. The element 
// incorporates interaction between axial/flexural and shear responses under 
// cyclic loading conditions by incorporating RC panel behavior based on the 
// fixed-strut angle approach (nDMaterial FSAM) into a two-dimensional fiber-based 
// MVLEM model. The element generates automatically m internal nodes with 1 DOF 
// at each macro-fiber (theNodesX with negative Tags) and adds them to the 
// domain. These internal DOFs are used to enforce equilibrium equation sigmaX=0 
// at each element macro-fiber (RC panel) in order to complete its strain field; 
// for details see referenced publications.
//
// References:
// 1) Kolozvari K., Orakcal K., and Wallace J. W. (2015). ”Modeling of Cyclic 
// Shear-Flexure Interaction in Reinforced Concrete Structural Walls. I: Theory”, 
// ASCE Journal of Structural Engineering, 141(5), 04014135 
// http://dx.doi.org/10.1061/(ASCE)ST.1943-541X.0001059
// 2) Kolozvari K., Tran T., Orakcal K., and Wallace, J.W. (2015). ”Modeling 
// of Cyclic Shear-Flexure Interaction in Reinforced Concrete Structural Walls. 
// II: Experimental Validation”, ASCE Journal of Structural Engineering, 141(5), 
// 04014136 http://dx.doi.org/10.1061/(ASCE)ST.1943-541X.0001083
// 3) Kolozvari K. (2013). “Analytical Modeling of Cyclic Shear-Flexure 
// Interaction in Reinforced Concrete Structural Walls”, PhD Dissertation, 
// University of California, Los Angeles.
//
// Source: /usr/local/cvs/OpenSees/SRC/element/SFI_MVLEM/SFI_MVLEM.h
//
// Rev: 1


#ifndef SFI_MVLEM_h
#define SFI_MVLEM_h
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
class Node;
class Channel;
class NDMaterial;
class Response;

class SFI_MVLEM : public Element
{
public:

	// constructors
	SFI_MVLEM(int tag,					// element tag
		int Nd1, int Nd2,				// end node tags          
		NDMaterial **Materials,			// array of material tags
		double *Thickness,				// array of macro-fiber thickness
		double *Width,					// array of macro-fiber widths
		int mm,							// number of macro-fibers (RC panels)
		double cc);						// center of rotation					

	SFI_MVLEM();

	// destructor
	~SFI_MVLEM();

  const char *getClassType(void) const {return "SFI_MVLEM2d";}
  
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
  int setupMacroFibers();
  
	// private attributes - a copy for each object of the class

	// input variables
	Node *theNodes[2];					// external node pointers          
	Node **theNodesX;					// array of internal node pointers
	Node **theNodesALL;					// array of ALL node pointers
	Node *theNd1;						// pointer to bottom node
	Node *theNd2;						// pointer to top node
	NDMaterial **theMaterial;			// array of ND materials
	Vector *theLoad;					// pointer to element load
	double c;						// center of rotation
	int m;						// no. of RC panels

	// calculated element parameters
	ID externalNodes;					// contains the id's of end nodes
	Matrix trans;						// hold the transformation matrix, could use a Vector
	double h;							// height of SFI_MVLEM element (undeformed configuration)
	double Lw;							// length of SFI_MVLEM elemtn, i.e. wall length
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
	double *SFI_MVLEMStrainX;			// macro-fiber axial strain in X direction (epsX)
	double *SFI_MVLEMStrainY;			// macro-fiber axial strain in Y direction (epsY)
	double *SFI_MVLEMStrainXY;			// macro-fiber shear strain (gammaXY)
	double *SFI_MVLEMStrain;			// macro-fiber strains 
	double *Dens;						// macro-fiber densities

	// for recorders
	double Dsh;							// recorded shear deformation
	Vector P_6DOF;						// recorded forces at end nodal 6DOFs

	// class wide matrices
	Matrix SFI_MVLEMK;					// stiffness
	Matrix SFI_MVLEMD;					// damping
	Matrix SFI_MVLEMM;					// mass 
	Vector SFI_MVLEMR;					// force

};
#endif
