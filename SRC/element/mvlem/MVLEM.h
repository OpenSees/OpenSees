// Code developed by: Kristijan Kolozvari, California State University, Fullerton 
//		      Kutay Orakcal, Bogazici University, Istanbul, Turkey
//		      John Wallace, UCLA
//
// Created: 07/2015
//
// Description: This file contains the class definition for Multiple-Vertical-
// Line-Element-Model (MVLEM; Vulcano et al., 1988; Orakcal et al., 2004). 
// Single model element is characterized with six global degrees of freedom, 
// three of each located at the center of the rigid top and bottom beams. 
// The flexural response is simulated by a series of vertical uniaxial elements 
// (fibers) connected to rigid beams at the top and bottom levels, whereas the shear 
// response is described via horizontal shear spring located at height ch from the 
// bottom of the element, and is uncoupled from the flexural modeling parameters.
//
// References:
// 1) Vulcano A., Bertero V.V., and Colotti V. (1988). “Analytical Modeling of RC 
// Structural Walls”, Proceedings, 9th World Conference on Earthquake Engineering, 
// V. 6, Tokyo-Kyoto, Japan, pp. 41-46.
// 2) Orakcal K., Conte J.P., and Wallace J.W. (2004). “Flexural Modeling of 
// Reinforced Concrete Structural Walls - Model Attributes”, ACI Structural Journal, 
// V. 101, No. 5, pp 688-698.

#ifndef MVLEM_h
#define MVLEM_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Node;
class Channel;
class UniaxialMaterial;         

class MVLEM : public Element {
 public:
  
  // constructors
  MVLEM(int tag,
	double Dens,
	int Nd1, int Nd2,
	UniaxialMaterial **materialsConcrete,
	UniaxialMaterial **materialsSteel,
	UniaxialMaterial **materialsShear,
	double *Rho,
	double *thickness,
	double *width,
	int mm,
	double cc);
  
  MVLEM();
  
  // destructor
  ~MVLEM();

  const char *getClassType(void) const {return "MVLEM2d";}
  
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
  const Matrix &getSecantStiff(void);
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

  void Print(OPS_Stream &s, int flag = 0);;
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  
 protected:
  
 private:
  // private member functions - only available to objects of the class
  double * computeCurrentStrain(void);
  double getCurvature(void);
  Vector getStrain(void);
  Vector getStressConcrete(void);
  Vector getStressSteel(void);
  Vector getShearFD(void);
  int setupMacroFibers();
  
  // private attributes - a copy for each object of the class
  ID  externalNodes;          			// contains the id's of end nodes
  Matrix trans;							// hold the transformation matrix, could use a Vector
  
  // Input variables
  Node *theNodes[2];						// external node pointers          
  UniaxialMaterial **theMaterialsConcrete; // pointers to Concrete uniaxial material
  UniaxialMaterial **theMaterialsSteel;	// pointers to Steel uniaxial material
  UniaxialMaterial **theMaterialsShear;	// pointers to Shear uniaxial material
  double density;					// material density
  double c;							// center of rotation
  int m;							// no. of RC panels
  Vector *theLoad;						// pointer to element load

  // calculated element parameters
  double h;								// height of MVLEM element (undeformed configuration)
  double Lw;								// length of MVLEM elemtn, i.e. wall length 
  double A;								// Wall cross-section area
  double NodeMass;						// nodal mass

  // caldulated element arrays
  double *x;
  double *b;								// fiber widths
  double *t;								// fiber thickness
  double *rho;							// fiber reinforcing ratio
  double *Ac;								// concrete area
  double *As;								// steel area

  double *MVLEMStrain;					// fiber strains (0 to m-1), shear deformation (m)
  
  // class wide matrices
  static Matrix MVLEMK;					// stiffness
  static Matrix MVLEMD;					// damping
  static Matrix MVLEMM;					// mass 
  static Vector MVLEMR;					// residual  
};

#endif
