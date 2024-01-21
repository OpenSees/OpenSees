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
// 1) Kolozvari K., Orakcal K., and Wallace J. W. (2015). �Modeling of Cyclic 
// Shear-Flexure Interaction in Reinforced Concrete Structural Walls. I: Theory�, 
// ASCE Journal of Structural Engineering, 141(5), 04014135 
// http://dx.doi.org/10.1061/(ASCE)ST.1943-541X.0001059
// 2) Kolozvari K., Tran T., Orakcal K., and Wallace, J.W. (2015). �Modeling 
// of Cyclic Shear-Flexure Interaction in Reinforced Concrete Structural Walls. 
// II: Experimental Validation�, ASCE Journal of Structural Engineering, 141(5), 
// 04014136 http://dx.doi.org/10.1061/(ASCE)ST.1943-541X.0001083
// 3) Kolozvari K. (2013). �Analytical Modeling of Cyclic Shear-Flexure 
// Interaction in Reinforced Concrete Structural Walls�, PhD Dissertation, 
// University of California, Los Angeles.

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>
#include "SFI_MVLEM.h"
#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <NDMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

#include <TclModelBuilder.h> // for creating/adding internal nodes (theNodesX) to the domain 


#include <elementAPI.h>

// Read input parameters and build the material
void *OPS_SFI_MVLEM(void)
{
  // Pointer to a uniaxial material that will be returned                       
  Element *theElement = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  // Parse the script for material parameters
  if (numArgs < 7) { 
    opserr << "Want: SFI_MVLEM eleTag Dens iNode jNode m c -thick {fiberThick} -width {fiberWidth} -rho {Rho} -matConcrete {matTagsConcrete} -matSteel {matTagsSteel} -matShear {matTagShear}\n";
    return 0;
  }

  int iData[4];
  double dData[1];
  
  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid int data for element SFI_MVLEM" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid c for element SFI_MVLEM " << iData[0] << endln;
    return 0;
  }

  int m = iData[3];
  const char *str = 0;

  double *theThickness = new double[m];
  double *theWidth = new double[m];
  int *matTags = new int[m];

  NDMaterial **theMaterials = new NDMaterial* [m];
  for (int i = 0; i < m; i++) {
    theThickness[i] = 0.0;
    theWidth[i] = 0.0;
    matTags[i] = 0;
    theMaterials[i] = 0;
  }
	
  numArgs = OPS_GetNumRemainingInputArgs();
  while (numArgs >= (m+1)) {
      //OPS_GetStringCopy(&str);
      str = OPS_GetString();
    if (strcmp(str, "-thick") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theThickness) != 0) {
	opserr << "Invalid thick parameter for SFI_MVLEM   " << iData[0] << endln;
	return 0;
      }
    } else if (strcmp(str, "-width") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theWidth) != 0) {
	opserr << "Invalid width value for SFI_MVLEM  " << iData[0] << endln;
	return 0;
      }
    } else if (strcmp(str, "-mat") == 0) {
      numData = m;
      if (OPS_GetIntInput(&numData, matTags) != 0) {
	opserr << "Invalid mat tags for SFI_MVLEM  " << iData[0] << endln;
	return 0;
      }
      for (int i=0; i<m; i++) {
	theMaterials[i] = 0;
	theMaterials[i] = OPS_getNDMaterial(matTags[i]);
	if (theMaterials[i] == 0) {
	  opserr << "Invalid material tag " << matTags[i] << "  for SFI_MVLEM  " << iData[0] << endln;
	  return 0;
	}
      }
    } 
    
    // clean up the str
    //    delete [] str;

    numArgs = OPS_GetNumRemainingInputArgs();
    
  }

  theElement = new SFI_MVLEM(iData[0], iData[1], iData[2], 
			     theMaterials, 
			     theThickness, theWidth , iData[3], dData[0]);
			 
  
  // Cleanup dynamic memory
  if (theThickness != 0) 
    delete [] theThickness;
  if (theWidth != 0) 
    delete [] theWidth;
  if (matTags != 0) 
    delete [] matTags;

  if (theMaterials != 0)
    delete[] theMaterials;

  return theElement;
}

// Typical constructor
SFI_MVLEM::SFI_MVLEM(int tag,
		     int Nd1, int Nd2,          
		     NDMaterial **materials,
		     double *thickness,
		     double *width,
		     int mm,
		     double cc)
  :Element(tag,ELE_TAG_SFI_MVLEM),  
   externalNodes(2+mm),
   theNd1(0),
   theNd2(0),
   theNodesX(0),
   theNodesALL(0),
   theMaterial(0), theLoad(0),
   SFI_MVLEMStrainX(0),SFI_MVLEMStrainY(0),SFI_MVLEMStrainXY(0),SFI_MVLEMStrain(0),
   x(0), b(0), AcX(0), AcY(0), kx(0), ky(0), kh(0), Fx(0), Fy(0), Fxy(0), Dens(0), Dx(0), Dy(0), Dxy(0),
   SFI_MVLEMK(6+m,6+m), SFI_MVLEMR(6+m), SFI_MVLEMD(6+m,6+m), SFI_MVLEMM(6+m,6+m),
   P_6DOF(6),
   m(mm),c(cc)
   
{	
  // Fill with ZEROs all element matrices
  SFI_MVLEMK.Zero();		// element stiffness matrix (6+m, 6+m)
  SFI_MVLEMR.Zero();		// element force vector (6+m)
  P_6DOF.Zero();			// element force vector (6) - external nodes only used in globalForce recorder
  
  TotalMass = 0.0;
  NodeMass = 0.0;
  h = 0.0;
  
  // Check number of fibers - max is 999 to avoid overlapping in internal node tags
  if (m > 999){
    opserr << "WARNING: Number of fibers assigned is " << m << ". Maximum allowed number of fibers is 999!\n";
    exit(-1);
  }
  
  // Fill in the ID containing external node info with node id's 
  if (externalNodes.Size() != 2+m)
    opserr << "FATAL SFI_MVLEM::SFI_MVLEM() - out of memory, could not create an ID of size 2+m\n";
  
  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;
  
  
  // Set external node pointers to NULL - external nodes
  theNodes[0] = 0; 
  theNodes[1] = 0;

  // Create a internal node tags equal to first to get past OpenSees check
  for (int i = 0; i < m; i++){ 
    externalNodes(i+2) = Nd1; 
  } 

  
  // Allocate memory for the m internal nodes
  theNodesX = new Node*[m];
  theNodesALL = new Node*[m+2];
  
  // Set NodeX pointers to NULL - m internal nodes
  for (int i=0; i < m; i++) {
    theNodesX[i] = 0; 
  }
  
  // Set theNodesALL pointers to NULL - m internal nodes
  for (int i=0; i < m+2; i++) {
    theNodesALL[i] = 0; 
  }
  
  
  // Check thickness and width input
  if (thickness == 0) {
    opserr << "SFI_MVLEM::SFI_MVLEM() - "
	   << "Null thickness array passed.\n";
    exit(-1);
  }
  
  if (width == 0) {
    opserr << "SFI_MVLEM::SFI_MVLEM() - "
	   << "Null width array passed.\n";
    exit(-1);
  }
  
  // Allocate memory for the thickness and width
  t  = new double[m];
  b  = new double[m];
  Lw = 0.0;
  
  for (int i=0; i<m; i++) {
    t[i] = thickness[i];
    b[i] = width[i];
    Lw += b[i]; // Total length of the wall
  }
  
  // Calculate locations of concrete macro-fibers in the cross-section (centerline - x = 0.0)
  x  = new double[m];
  for (int i=0; i < m; i++)
    x[i]=0.0;
  
  for (int i=0; i < m; i++) {
    double sumb_i = 0.0;
    for (int j=0; j<i+1; j++)
      sumb_i += b[j];
    
    x[i] = (sumb_i - b[i]/2.0) - Lw/2.0;
  }
  
  
  // Check material input
  if (materials == 0) {
    opserr << "SFI_MVLEM::SFI_MVLEM() - "
	   << "Null material array passed.\n";
    exit(-1);
  }
  
  // Allocate memory for the ND materials
  theMaterial = new NDMaterial*[m];
  
  if (theMaterial == 0) {
    opserr << "SFI_MVLEM::SFI_MVLEM() - "
	   << "Failed to allocate pointers for uniaxial materials.\n";
    exit(-1);
  }
  
  // Get copies of the ND materials
  for (int i=0; i < m; i++) {
    if (materials[i] == 0) {
      opserr << "SFI_MVLEM::SFI_MVLEM() - "
	"Null ND material pointer passed.\n";
      exit(-1);
    }
    
    theMaterial[i] = materials[i]->getCopy("PlaneStress2D");
    
    if (theMaterial[i] == 0) {
      opserr << "SFI_MVLEM::SFI_MVLEM() - "
	     << "Failed to copy ND material.\n";
      exit(-1);
    }
  }
  
  // Allocate memory for element arrays
  // Area of concrete fibers
  AcX = new double[m];
  AcY = new double[m];
  
  // Panel stiffness (trial)
  kx = new double[m];
  ky = new double[m];
  kh = new double[1];
  
  // Panel force (trial)
  Fx = new double[m];
  Fy = new double[m];
  Fxy = new double[m];
  
  // Panel stiffness (trial)
  Dx = new double[m];
  Dy = new double[m];
  Dxy = new double[m];
  
  // Panel strains
  SFI_MVLEMStrainX = new double[m];
  SFI_MVLEMStrainY = new double[m];
  SFI_MVLEMStrainXY = new double[m];
  SFI_MVLEMStrain = new double[3*m]; 
  
  // Density
  Dens = new double[m];
  
  // Assign zero to element arrays
  for (int i=0; i < m; i++) {
    
    AcX[i]=0.0;
    AcY[i]=0.0;
    
    kx[i]=0.0;
    ky[i]=0.0;
    
    Fx[i]=0.0;
    Fy[i]=0.0;
    Fxy[i]=0.0;
    
    Dx[i]=0.0;
    Dy[i]=0.0;
    Dxy[i]=0.0;
    
    SFI_MVLEMStrainX[i]=0.0;
    SFI_MVLEMStrainY[i]=0.0;
    SFI_MVLEMStrainXY[i]=0.0;
    
    SFI_MVLEMStrain[i]=0.0;
    SFI_MVLEMStrain[i+m]=0.0;
    SFI_MVLEMStrain[i+2*m]=0.0;
    
    Dens[i]=0.0;
  } 
  
  kh[0]=0.0;
  
  // Calculate concrete areas in X and Y directions
  for (int i=0; i < m; i++) {
    AcX[i] = h*t[i];
    AcY[i] = b[i]*t[i];
  }
  
  // Get panel density from 2-D materials
  for (int i=0; i < m; i++) {
    Dens[i] = theMaterial[i]->getRho();
  }
  
}  

// Constructor which should be invoked by an FE_ObjectBroker only
SFI_MVLEM::SFI_MVLEM()
  :Element(0,ELE_TAG_SFI_MVLEM),  
   externalNodes(2+m),
   theNd1(0),
   theNd2(0),
   theNodesX(0),
   theNodesALL(0),
   theMaterial(0), theLoad(0),
   SFI_MVLEMStrainX(0),SFI_MVLEMStrainY(0),SFI_MVLEMStrainXY(0),SFI_MVLEMStrain(0),
   x(0), b(0), AcX(0), AcY(0), kx(0), ky(0), kh(0), Fx(0), Fy(0), Fxy(0), Dens(0), Dx(0), Dy(0), Dxy(0),
   SFI_MVLEMK(6+m,6+m), SFI_MVLEMR(6+m), SFI_MVLEMD(6+m,6+m), SFI_MVLEMM(6+m,6+m),
   P_6DOF(6),
   m(0),c(0)
{
  if (externalNodes.Size() != 2+m)
    opserr << "FATAL SFI_MVLEM::SFI_MVLEM() - out of memory, could not create an ID of size 2\n";
  
  theNodes[0] = 0; 
  theNodes[1] = 0;
  
  // Allocate memory for all nodes (internal and external)
  theNodesX = new Node*[m];		// internal nodes
  theNodesALL = new Node*[m+2];	// all nodes
  
  // Set NodeX pointers to zero
  for (int i=0; i < m; i++)
    {
      theNodesX[i] = 0; 
    }
  
  // Set theNodesALL pointers to zero
  for (int i=0; i < m+2; i++)
    {
      theNodesALL[i] = 0; 
    }
  
  SFI_MVLEMK.Zero();		// element stiffness matrix
  SFI_MVLEMR.Zero();		// element force vector (6+m)
  P_6DOF.Zero();			// element force vector (6) - external nodes only - for global force recorder
  SFI_MVLEMD.Zero();      // element damping matrix
  SFI_MVLEMM.Zero();		// element mass matrix

}

//  Destructor - provided to clean up any memory 
SFI_MVLEM::~SFI_MVLEM()
{
	// clean up the memory associated with the element, this is
	// memory the SFI_MVLEM objects allocates and memory allocated 
	// by other objects that the SFI_MVLEM object is responsible for 
	// cleaning up, i.e. the MaterialObject.

	if (theMaterial != 0)  {
		for (int i=0; i < m; i++)
			if (theMaterial[i] != 0)
				delete theMaterial[i];
		delete [] theMaterial;
	}      

	if (theLoad != 0)
		delete theLoad;
	if(x!=0)
		delete []x;
	if(b!=0)
		delete []b;
	if(t!=0)
		delete []t;	
	if(AcX!=0)
		delete []AcX;
	if(AcY!=0)
		delete []AcY;
	if(kx!=0)
		delete []kx;
	if(ky!=0)
		delete []ky;
	if(kh!=0)
		delete []kh;
	if(Fx!=0)
		delete []Fx;
	if(Fy!=0)
		delete []Fy;
	if(Fxy!=0)
		delete []Fxy;
	if(Dens!=0)
		delete []Dens;
	if(Dx!=0)
		delete []Dx;
	if(Dy!=0)
		delete []Dy;
	if(Dxy!=0)
		delete []Dxy;
	if(SFI_MVLEMStrainX!=0)
		delete []SFI_MVLEMStrainX;
	if(SFI_MVLEMStrainY!=0)
		delete []SFI_MVLEMStrainY;
	if(SFI_MVLEMStrainXY!=0)
		delete []SFI_MVLEMStrainXY;
	if(SFI_MVLEMStrain!=0)
		delete []SFI_MVLEMStrain;
	if(theNodesX!=0)
		delete [] theNodesX;
	if(theNodesALL!=0)
		delete [] theNodesALL;
}

// Get number of nodes (external + internal)
int SFI_MVLEM::getNumExternalNodes(void) const
{
	return 2+m;
}

// Get node tags
const ID & SFI_MVLEM::getExternalNodes(void) 
{
	return externalNodes;
}

// Get shear deformation
double SFI_MVLEM::getShearDef(void) 
{
	return Dsh;
}

// Get curvature (from vertical strains)
double SFI_MVLEM::getCurvature(void) 
{
	double Curv;

	Curv = (SFI_MVLEMStrainY[0] - SFI_MVLEMStrainY[m-1])/(x[0]-x[m-1]);

	return Curv;
}

// Get global forces at 6 DOFs (top and bottom node)
Vector SFI_MVLEM::getResistingForce_6DOF(void) 
{
  for (int i=0; i < 6; i++) {
    P_6DOF(i) = SFI_MVLEMR(i);
  }
  
  return P_6DOF;
}

// Get node pointers
Node ** SFI_MVLEM::getNodePtrs(void)
{
	
  // Pack external and internal node pointers into one array
  for (int i=0; i < 2; i++) {
    theNodesALL[i] = theNodes[i];
  }
  
  for (int i=2; i < m+2; i++) {
    theNodesALL[i] = theNodesX[i-2];
  }
  
  return theNodesALL;
}

// Get number of DOFs 
int SFI_MVLEM::getNumDOF(void) {
  
  int NumDOF = 6+m; // 3 DOFs per external nodes, 1 DOF per internal nodes
  
  return NumDOF;
}

// Set Domain
void SFI_MVLEM::setDomain(Domain *theDomain) 
{
  // Check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0)
    {
      return;
    }
  
  // Set node pointers to NULL
  theNodes[0] = 0; 
  theNodes[1] = 0;
  
  for (int i=0; i < m; i++) {
    theNodesX[i] = 0; 
  }
  
  // First ensure nodes (external) exist in Domain and set the node pointers
  int Nd1 = externalNodes(0);
  int Nd2 = externalNodes(1);
  
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);	
  
  // Get coordinates of end nodes - used for defining locations of internal nodes
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();
  
  // Calculate the element height and perform checks
  h = end2Crd(1) - end1Crd(1);
  
  if (h < 0.0) {
    opserr << "WARNING: Element height is negative. Define Nodes from bottom to top!";
    return;
  }
  
  if (h == 0.0) {
    opserr << "WARNING: Element height is ZERO!";
    return;
  }

  // Calculate concrete areas in X and Y directions
  for (int i=0; i < m; i++) {
    AcX[i] = h*t[i];
  }
  
  // Currently element can be only vertical
  if (end2Crd(0) != end1Crd(0)) {
    opserr << "WARNING: Element is NOT vertical!";
  }
  
  int eletag = this->getTag();
  // Create a internal node tag
  for (int i = 0; i < m; i++){ // Large NEGATIVE integer starting with tag of the element
    externalNodes(i+2) = -(eletag*1000 + i + 1); // Max fibers is 999 to avoid overlap
  } 

  // Build m internal nodes (NodesX) and add them to the domain
  for (int i = 0; i < m; i++) {
    
    int nodeId_temp = externalNodes(i+2); // Store node tag temporarily 
    
    // Create coordinates wrt top and bottom element node
    double xLoc_temp = end1Crd(0) + x[i]; 
    double yLoc_temp = 0.5*(end1Crd(1)+end2Crd(1)); // Mid-height
    double zloc_temp = end1Crd(2);					// Not currently used since Domain is 2D
    
    // Create Node and add it to the domain
    Node *theNode = 0;
    
    theNode = new Node(nodeId_temp, 1, xLoc_temp, yLoc_temp); // create internal node with 1 DOF	
    
    if (theNode == 0) {
      opserr << "WARNING ran out of memory creating node\n";
      opserr << "node: " << nodeId_temp << " in SFI_MVLEM." << endln; 
      exit(-1);
    }
    
    if (theDomain->addNode(theNode) == false) { // add internal node to the domain
      opserr << "WARNING failed to add node to the domain\n";
      opserr << "node: " << nodeId_temp << " in SFI_MVLEM." << endln;
      delete theNode; // otherwise memory leak
      exit(-1);
    } 
  } // END create/add internal nodes
  
  if (theNodes[0] == 0) 
    {
      opserr << "WARNING SFI_MVLEM::setDomain() - at SFI_MVLEM " << this->getTag() << " node " <<
	Nd1 << " does not exist in domain\n";
      return;  // Don't go any further - otherwise segemntation fault
    }
  if (theNodes[1] == 0)
    {        
      opserr << "WARNING SFI_MVLEM::setDomain() - at SFI_MVLEM " << this->getTag() << " node " <<
	Nd2 << " does not exist in domain\n";
      return;
    }	
  
  // First ensure NodesX (internal - dummy) exist in Domain and set the node pointers
  for (int i=0; i<m; i++) {
    
    int NdX_temp1 = externalNodes(i+2);
    
    theNodesX[i] = theDomain->getNode(NdX_temp1);
    
    if (theNodesX[i] == 0) 
      {
	opserr << "WARNING SFI_MVLEM::setDomain() - at SFI_MVLEM " << this->getTag() << " node " <<
	  NdX_temp1 << " does not exist in domain\n";
	return;  // Don't go any further - otherwise segemntation fault
      }
    
  }
  
  // Call the DomainComponent class method 
  this->DomainComponent::setDomain(theDomain);
  
  // Ensure connected nodes have correct number of dof's
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  if ((dofNd1 != 3) || (dofNd2 != 3)) // 3 DOFS at external nodes
    {
      opserr << "SFI_MVLEM::setDomain(): 3 dof required at nodes, " << dofNd1 << " and "
	     <<  dofNd2 << " provided\n";
    }
  
  for (int i=0; i < m; i++) {
    int dofNdXi = theNodesX[i]->getNumberDOF();
    if (dofNdXi != 1)				// 1 DOF at internal nodes 
      {
	opserr << "SFI_MVLEM::setDomain(): 1 dof required at internal nodes, " << dofNdXi << " provided\n";
      }
  }
  
  // Calculate the nodal mass (external nodes only) for lumped mass approach
  for (int i=0; i < m; i++) {
    TotalMass += Dens[i] * AcY[i]*h;
  }
  
  NodeMass = TotalMass/2.0;                                       
  
  // Create a vector to hop applied loads - NOT used in the current model formulation (no element loads)
  if (theLoad == 0)
    theLoad = new Vector(6+m);
  if (theLoad == 0)  {
    opserr << "SFI_MVLEM::setDomain() - element: " << this->getTag()
	   << " out of memory creating vector of size: " << 6+m << endln;
    return;
  }   
}

// Commit state
int SFI_MVLEM::commitState()
{
  int errCode = 0;
  
  // Commit material models
  for (int i=0; i < m; i++)
    
    errCode += theMaterial[i]->commitState();
  
  return errCode;
}

// Revert to last committed state (if convergence is not achieved)
int SFI_MVLEM::revertToLastCommit()
{
  int errCode = 0;
  
  // Revert material models
  for (int i=0; i < m; i++)
    
    errCode += theMaterial[i]->revertToLastCommit();
  
  return errCode;
}

// Revert to start
int SFI_MVLEM::revertToStart()
{
  
  int errCode = 0;
  
  // Revert material models
  for (int i=0; i < m; i++)
    errCode += theMaterial[i]->revertToStart();
  
  // Compute initial stiffness
  this->getInitialStiff();
  
  return errCode;
  
}

// Update state
int SFI_MVLEM::update()
{
  
  // Get the current strain given trial displacements at nodes
  this->computeCurrentStrain();                 
  
  // Set the strain in the materials
  int errCode1 = 0;
  
  for (int i=0; i < m; i++) {
    
    Vector strain(3);
    
    strain(0) = SFI_MVLEMStrain[i];
    strain(1) = SFI_MVLEMStrain[i+m];
    strain(2) = SFI_MVLEMStrain[i+2*m];
    
    // Set trial response for material models
    errCode1 += theMaterial[i]->setTrialStrain(strain); 
  }
  return errCode1 ;
}

// Send Self
int SFI_MVLEM::sendSelf(int commitTag, Channel &theChannel)
{
	int res;
	int dataTag = this->getDbTag();
	
	static Vector data(3);  // One bigger than needed so no clash later

	data(0) = this->getTag();
	data(1) = m;
	data(2) = c;
	
	// SFI_MVLEM then sends the tags of it's nodes
	res = theChannel.sendID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING SFI_MVLEM::sendSelf() - failed to send ID\n";
		return -2;
	}

	// Send the material class tags
	ID matClassTags(m);
	for (int i=0; i<m; i++)
		matClassTags(i) = theMaterial[i]->getClassTag();
	res = theChannel.sendID(0, commitTag, matClassTags);

	// Send the material models
	for (int i=0; i<m; i++)
		theMaterial[i]->sendSelf(commitTag, theChannel);
	return 0;
}

// Receive Self
int SFI_MVLEM::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res;
	int dataTag = this->getDbTag();

	// SFI_MVLEM creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	// delete dynamic memory
	if (theMaterial != 0)  {
		for (int i=0; i < m; i++)
			if (theMaterial[i] != 0)
				delete theMaterial[i];
		delete [] theMaterial;
	}

	Vector data(3); // One bigger than needed so no clash later
	res = theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING SFI_MVLEM::recvSelf() - failed to receive Vector\n";
		return -1;
	}	      

	this->setTag((int)data(0));

	data(0) = this->getTag();
	data(1) = m;
	data(2) = c;

	// SFI_MVLEM now receives the tags of it's two external nodes
	res = theChannel.recvID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING SFI_MVLEM::recvSelf() - failed to receive ID\n";
		return -2;
	}

	// Receive the material class tags
	ID matClassTags(m);
	res=theChannel.recvID(0, commitTag, matClassTags);

	// Allocate memory for the uniaxial materials
	theMaterial = new NDMaterial* [m];
	if (theMaterial == 0)  {
		opserr << "SFI_MVLEM::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Receive the material models
	for (int i=0; i < m; i++)  {
		theMaterial[i] = theBroker.getNewNDMaterial(matClassTags(i));
		if (theMaterial[i] == 0) {
			opserr << "SFI_MVLEM::recvSelf() - "
				<< "failed to get blank uniaxial material.\n";
			return -3;
		}
		theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
	}
	return 0;
}

// Display model
int SFI_MVLEM::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// Get the end points of the beam for the display factor
	static Vector v1(3);
	static Vector v2(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);

	// determine the deformation - rotation - other is taken from v1, v2
	static Vector r1(1);
	theNodes[0]->getDisplayRots(r1, fact, displayMode);

	// Displaying wall axis
	int error = 0;
	Vector RGB(3);
	RGB(0) = 0.0;
	RGB(1) = 1.0;
	RGB(2) = 0.0;
	error += theViewer.drawLine (v1, v2, RGB, RGB, 1, 1);

	// Displaying Panels
	for (int panel=0; panel < m; panel++) // loop over m panels
	{ 
		Matrix NodePLotCrds(m,13); // (panel id, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)

		// first set the quantity to be displayed at the nodes;
		// if displayMode is 1 through 3 we will plot material stresses otherwise 0.0
		static Vector values(1); // values to be plotted (either epsX, epsY, gammaXY)
		if (displayMode < 4 && displayMode > 0) {
			const Vector& stress = theMaterial[panel]->getStrain();
			values(0) = stress(displayMode - 1);
		}
		else {
			values(0) = 0.0;
		}

		// Panel nodes
		NodePLotCrds(panel,0) = panel+1; // panel id
		// Local node 1 - bottom left
		NodePLotCrds(panel,1) = v1(0) + x[panel]-b[panel]/2.0; // x 
		NodePLotCrds(panel,2) = v1(1) + (x[panel]-b[panel]/2.0)*r1(0); // y
		NodePLotCrds(panel,3) = v1(2); // z
		// Local node 2 - bottom right
		NodePLotCrds(panel,4) = v1(0) + x[panel]+b[panel]/2.0; // x
		NodePLotCrds(panel,5) = v1(1) + (x[panel]+b[panel]/2.0)* r1(0); // y
		NodePLotCrds(panel,6) = v1(2); // z
		// Local node 3 - top left
		NodePLotCrds(panel,7) = v2(0) + x[panel]+b[panel]/2.0; // x
		NodePLotCrds(panel,8) = v2(1) +(x[panel]+b[panel]/2.0)* r1(0); // y
		NodePLotCrds(panel,9) = v2(2); // z
		// Local node 4 - top right
		NodePLotCrds(panel,10) = v2(0) + x[panel]-b[panel]/2.0; // x
		NodePLotCrds(panel,11) = v2(1) + (x[panel]-b[panel]/2.0)* r1(0); // y
		NodePLotCrds(panel,12) = v2(2); // z

		Matrix coords(4,3); // Temporary coordinates for plotting
		
		coords(0,0) = NodePLotCrds(panel,1); // node 1 x
		coords(1,0) = NodePLotCrds(panel,4); // node 2 x
		coords(2,0) = NodePLotCrds(panel,7); // node 3 x
		coords(3,0) = NodePLotCrds(panel,10); // node 4 x

		coords(0,1) = NodePLotCrds(panel,2); // node 1 y
		coords(1,1) = NodePLotCrds(panel,5); // node 2 y
		coords(2,1) = NodePLotCrds(panel,8); // node 3 y
		coords(3,1) = NodePLotCrds(panel,11); // node 4 y

		coords(0,2) = NodePLotCrds(panel,3); // node 1 z
		coords(1,2) = NodePLotCrds(panel,6); // node 2 z
		coords(2,2) = NodePLotCrds(panel,9); // node 3 z
		coords(3,2) = NodePLotCrds(panel,12); // node 4 z

		error += theViewer.drawPolygon (coords, values);
	}

	return error;
}

// Print Element Information
void SFI_MVLEM::Print(OPS_Stream &s, int flag)
{
	if (flag == 0)  {    
		s << "SFI_MVLEM Element tag: " << this->getTag() << endln;
		s << "iNode: " << externalNodes(0) << ", jNode: " << externalNodes(1) << endln;
		s << "Element height: " << h << endln;
		s << "Number of RC panel elements: " << m << endln;

		// get resisting forces in global system
		s << "Global resisting forces: " << this->getResistingForce_6DOF();
		
		for (int i=0; i<m; i++)  {              
			s << "\nPanel #: " << i+1 << endln;
			theMaterial[i]->Print(s, flag);
		}

	} else if (flag == 1)  {        
		// does nothing
	} 
}

// Set element responses
Response *SFI_MVLEM::setResponse(const char **argv, int argc, OPS_Stream &s) 
    
{
	Response *theResponse = 0;

	s.tag("ElementOutput");
    s.attr("eleType","SFI_MVLEM");
    s.attr("eleTag",this->getTag());
    s.attr("node1",externalNodes[0]);
    s.attr("node2",externalNodes[1]);

	// Global forces
	if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {
			
			s.tag("ResponseType","Fx_i");
			s.tag("ResponseType","Fy_i");
			s.tag("ResponseType","Mz_i");
			s.tag("ResponseType","Fx_j");
			s.tag("ResponseType","Fy_j");
			s.tag("ResponseType","Mz_j");
			
			theResponse = new ElementResponse(this, 1, Vector(6));

	}  

	// Shear deformation
	else if (strcmp(argv[0],"ShearDef") == 0 || strcmp(argv[0],"sheardef") == 0) {
		
			s.tag("ResponseType","Dsh");
		
		theResponse = new ElementResponse(this, 2, 0.0);

	} 
	
	// Element curvature
	else if (strcmp(argv[0],"Curvature") == 0 || strcmp(argv[0],"curvature") == 0) {
		
			s.tag("ResponseType","fi");
		
		theResponse = new ElementResponse(this, 3, 0.0);
	}

	// Material output
    else if (strcmp(argv[0],"RCpanel") == 0 || strcmp(argv[0],"RCPanel")  
		      || strcmp(argv[0],"RC_panel")  || strcmp(argv[0],"RC_Panel") == 0) 
	{
		
		// Check if correct # of arguments passed
		if (argc != 3) {
			opserr << "WARNING: Number of recorder input for RC Panel is: " << argc-1 << "; should be 2: panTag (one panel only: 1 to m) and $Response_Type.\n";
			return 0;
		}

		int matNum = atoi(argv[1]);
	  
		s.tag("Material");
		s.attr("number",matNum);
		
		theResponse = theMaterial[matNum-1]->setResponse(&argv[argc-1], argc-2, s);

    }

	s.endTag();
	
	return theResponse;
}

// Obtain element responses
int SFI_MVLEM::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) 
	{    
	case 1:  // Global forces
		return eleInfo.setVector(this->getResistingForce_6DOF());  

	case 2:  // Shear deformation
		return eleInfo.setDouble(this->getShearDef());

	case 3:  // Curvature
		return eleInfo.setDouble(this->getCurvature());

	default:

		return 0;

	}
}

// Get the element initial element tangent matrix
const Matrix & SFI_MVLEM::getInitialStiff(void)
{
  double Kh=0.0;
  
  for (int i=0; i < m; i++)
    {
      // Get material initial tangent
      const Matrix &D = theMaterial[i]->getInitialTangent();
      
      double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
      double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
      double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);
      
      kx[i]  = D00 * h*t[i] / b[i];
      ky[i]  = D11 * b[i]*t[i] / h;
      Kh  += D22 * b[i]*t[i] / h;
      
    }
  
  // Build the initial stiffness matrix
  double Kv=0.0; double Km=0.0; double e=0.0; double ex=0.0;
  
  for(int i=0; i<m; ++i)
    {
      Kv+=ky[i];   
      Km+=ky[i]*x[i]*x[i];
      e+=ky[i]*x[i];
      
      SFI_MVLEMK(6+i,6+i) = kx[i]; // Diagonal terms accounting for horizontal stiffness
    }                  
  
  SFI_MVLEMK(0,0) = Kh;
  SFI_MVLEMK(0,1) = 0.0;
  SFI_MVLEMK(0,2) = -Kh*c*h;
  SFI_MVLEMK(0,3) = -Kh;
  SFI_MVLEMK(0,4) = 0.0;
  SFI_MVLEMK(0,5) = -Kh*(1-c)*h;
  
  SFI_MVLEMK(1,0) = SFI_MVLEMK(0,1);
  SFI_MVLEMK(1,1) = Kv;
  SFI_MVLEMK(1,2) = e;
  SFI_MVLEMK(1,3) = 0.0;
  SFI_MVLEMK(1,4) = -Kv;
  SFI_MVLEMK(1,5) = -e;
  
  SFI_MVLEMK(2,0) = SFI_MVLEMK(0,2);
  SFI_MVLEMK(2,1) = SFI_MVLEMK(1,2);
  SFI_MVLEMK(2,2) = h*h*c*c*Kh+Km;
  SFI_MVLEMK(2,3) = h*c*Kh;
  SFI_MVLEMK(2,4) = -e;
  SFI_MVLEMK(2,5) = (1-c)*c*h*h*Kh-Km;
  
  SFI_MVLEMK(3,0) = SFI_MVLEMK(0,3);
  SFI_MVLEMK(3,1) = SFI_MVLEMK(1,3);
  SFI_MVLEMK(3,2) = SFI_MVLEMK(2,3);
  SFI_MVLEMK(3,3) = Kh;
  SFI_MVLEMK(3,4) = 0.0;
  SFI_MVLEMK(3,5) = Kh*(1-c)*h;
  
  SFI_MVLEMK(4,0) = SFI_MVLEMK(0,4);
  SFI_MVLEMK(4,1) = SFI_MVLEMK(1,4);
  SFI_MVLEMK(4,2) = SFI_MVLEMK(2,4);
  SFI_MVLEMK(4,3) = SFI_MVLEMK(3,4);
  SFI_MVLEMK(4,4) = Kv;
  SFI_MVLEMK(4,5) = e;
  
  SFI_MVLEMK(5,0) = SFI_MVLEMK(0,5);
  SFI_MVLEMK(5,1) = SFI_MVLEMK(1,5);
  SFI_MVLEMK(5,2) = SFI_MVLEMK(2,5);
  SFI_MVLEMK(5,3) = SFI_MVLEMK(3,5);
  SFI_MVLEMK(5,4) = SFI_MVLEMK(4,5);
  SFI_MVLEMK(5,5) = (1-c)*(1-c)*h*h*Kh+Km;
  
  for(int i=0; i<6+m; ++i)
    {
      if (SFI_MVLEMK(i,i) == 0.0)
	{ 
	  opserr << "Singular SFI_MVLEM_K/n"; 
	}
    }   
  
  // Return element stiffness matrix
	return SFI_MVLEMK;

}

// Get current element tangent stiffness matrix from the material for the last updated strain
const Matrix & SFI_MVLEM::getTangentStiff(void)
{   
	
	double Kh=0.0;
	
	for (int i=0; i < m; i++)
	{
		// Get the material tangent
		const Matrix &D = theMaterial[i]->getTangent();

		double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
		double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
		double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);

		kx[i]  = D00 * h*t[i] / b[i];
		ky[i]  = D11 * b[i]*t[i] / h;
		Kh  += D22 * b[i]*t[i] / h;

	}

	// Build the tangent stiffness matrix
	double Kv=0.0; double Km=0.0; double e=0.0; double ex=0.0;

	for(int i=0; i<m; ++i)
	{
		Kv+=ky[i];   
		Km+=ky[i]*x[i]*x[i];
		e+=ky[i]*x[i];

		SFI_MVLEMK(6+i,6+i) = kx[i]; // Diagonal terms accounting for horizontal stiffness
	}                  

	SFI_MVLEMK(0,0) = Kh;
	SFI_MVLEMK(0,1) = 0.0;
	SFI_MVLEMK(0,2) = -Kh*c*h;
	SFI_MVLEMK(0,3) = -Kh;
	SFI_MVLEMK(0,4) = 0.0;
	SFI_MVLEMK(0,5) = -Kh*(1-c)*h;

	SFI_MVLEMK(1,0) = SFI_MVLEMK(0,1);
	SFI_MVLEMK(1,1) = Kv;
	SFI_MVLEMK(1,2) = e;
	SFI_MVLEMK(1,3) = 0.0;
	SFI_MVLEMK(1,4) = -Kv;
	SFI_MVLEMK(1,5) = -e;

	SFI_MVLEMK(2,0) = SFI_MVLEMK(0,2);
	SFI_MVLEMK(2,1) = SFI_MVLEMK(1,2);
	SFI_MVLEMK(2,2) = h*h*c*c*Kh+Km;
	SFI_MVLEMK(2,3) = h*c*Kh;
	SFI_MVLEMK(2,4) = -e;
	SFI_MVLEMK(2,5) = (1-c)*c*h*h*Kh-Km;

	SFI_MVLEMK(3,0) = SFI_MVLEMK(0,3);
	SFI_MVLEMK(3,1) = SFI_MVLEMK(1,3);
	SFI_MVLEMK(3,2) = SFI_MVLEMK(2,3);
	SFI_MVLEMK(3,3) = Kh;
	SFI_MVLEMK(3,4) = 0.0;
	SFI_MVLEMK(3,5) = Kh*(1-c)*h;

	SFI_MVLEMK(4,0) = SFI_MVLEMK(0,4);
	SFI_MVLEMK(4,1) = SFI_MVLEMK(1,4);
	SFI_MVLEMK(4,2) = SFI_MVLEMK(2,4);
	SFI_MVLEMK(4,3) = SFI_MVLEMK(3,4);
	SFI_MVLEMK(4,4) = Kv;
	SFI_MVLEMK(4,5) = e;

	SFI_MVLEMK(5,0) = SFI_MVLEMK(0,5);
	SFI_MVLEMK(5,1) = SFI_MVLEMK(1,5);
	SFI_MVLEMK(5,2) = SFI_MVLEMK(2,5);
	SFI_MVLEMK(5,3) = SFI_MVLEMK(3,5);
	SFI_MVLEMK(5,4) = SFI_MVLEMK(4,5);
	SFI_MVLEMK(5,5) = (1-c)*(1-c)*h*h*Kh+Km;

	for(int i=0; i<6+m; ++i)
	{
		if (SFI_MVLEMK(i,i) == 0.0)
		{ 
			opserr << "Singular SFI_MVLEM_K/n"; 
		}
	} 

	// Return element stiffness matrix
	return SFI_MVLEMK;

}

// Get element mass matrix assuming lumped mass
const Matrix & SFI_MVLEM::getMass(void)
{ 
	SFI_MVLEMM.Zero();

	// No rotational mass, no mass at internal (dummy) nodes
	SFI_MVLEMM(0,0) = NodeMass; 
	SFI_MVLEMM(1,1) = NodeMass; 
	SFI_MVLEMM(3,3) = NodeMass; 
	SFI_MVLEMM(4,4) = NodeMass; 

	// Return element mass matrix
	return SFI_MVLEMM;
}

// Get element damping matrix
const Matrix & SFI_MVLEM::getDamp(void)
{   
	SFI_MVLEMD.Zero();

	SFI_MVLEMD = this->Element::getDamp();

	// Return element damping matrix
	return SFI_MVLEMD;
}

// N/A to this model - no element loads
void SFI_MVLEM::zeroLoad(void)
{
	// does nothing 
}

// N/A to this model - no element loads
int SFI_MVLEM::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	return 0;
}

// N/A to this model - no element loads
int SFI_MVLEM::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

// Get element force vector
const Vector & SFI_MVLEM::getResistingForce()
{	
  // Get the current force matrix from panel stresses
  for (int i=0; i < m; i++)
    {
      // Get the material stress
      const Vector &Stress = theMaterial[i]->getStress();

      double fx = Stress(0); 
      double fy = Stress(1); 
      double tauxy = Stress(2); 
      
      Fx[i] = fx * AcX[i];
      Fy[i] = fy * AcY[i];
      Fxy[i] = tauxy * AcY[i];
    }

  // Build force vector 
  double Fh = 0.0; // Force in horizontal spring (at location c*h)
  double Fysum = 0.0; // Sum of vertical forces
  
  for (int i=0; i < m; i++)
    {
      Fh += -1.0*Fxy[i];
      Fysum += Fy[i];
      SFI_MVLEMR[6+i] = Fx[i]; // Force on internal (dummy) DOFs
    }
  
  SFI_MVLEMR(0) = Fh;
  SFI_MVLEMR(1) = -Fysum;
  SFI_MVLEMR(2) = -Fh*c*h;
  SFI_MVLEMR(3) = -Fh;
  SFI_MVLEMR(4) = Fysum;
  SFI_MVLEMR(5) = -Fh*(1-c)*h;
  
  for (int i=0; i<m; i++) {
    SFI_MVLEMR(2) -= Fy[i]*x[i];        
    SFI_MVLEMR(5) += +Fy[i]*x[i]; 
  }
  
  // Return element force vector
  return SFI_MVLEMR;
}

// Get current strains at RC panels (macro-fibers)
void SFI_MVLEM::computeCurrentStrain(void)   
{
  const Vector &disp1 = theNodes[0]->getTrialDisp(); // DOFs D1,D2,D3
  const Vector &disp2 = theNodes[1]->getTrialDisp(); // DOFs D4,D5,D6
  
  for (int i=0; i<m; i++){
    const Vector &dispXi = theNodesX[i]->getTrialDisp(); // 1 DOF at theNodesX - Vector of size 1
    Dx[i] = dispXi(0); // get displacements in X direction from the nodes
  }
  
  // Deformations at each RC panel (macro-fiber)
  for(int i=0; i<m; i++) {
    Dy[i] = -disp1(1) - x[i]*disp1(2) + disp2(1) + x[i]*disp2(2); 
    Dxy[i] = disp1(0) - disp2(0) - c*h*disp1(2) - (1-c)*h*disp2(2);
  }
  
  Dsh = - Dxy[0]; // Store shear deformations for the recorder
  
  // Strains at each RC panel (macro-fiber)
  for(int i=0;i<m;i++) {
    SFI_MVLEMStrainX[i]  =  Dx[i]/b[i];
    SFI_MVLEMStrainY[i]  =  Dy[i]/h;
    SFI_MVLEMStrainXY[i] = - Dxy[i]/h;
  }
  
  // Store strains into a single vector
  for(int i=0; i<m; i++) {
    SFI_MVLEMStrain[i] = SFI_MVLEMStrainX[i];
    SFI_MVLEMStrain[i+m] = SFI_MVLEMStrainY[i];
    SFI_MVLEMStrain[i+2*m] = SFI_MVLEMStrainXY[i];
  }
}

// Get resisting force increment from inertial forces
const Vector & SFI_MVLEM::getResistingForceIncInertia()
{
  // compute the current resisting force
  this->getResistingForce();  
  
  if (TotalMass != 0.0) {
    
    // Get nodal accelerations
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
    SFI_MVLEMR(0) += NodeMass*accel1(0);
    SFI_MVLEMR(1) += NodeMass*accel1(1);
    SFI_MVLEMR(3) += NodeMass*accel2(0);
    SFI_MVLEMR(4) += NodeMass*accel2(1);
    
    // Add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      SFI_MVLEMR += this->getRayleighDampingForces(); 
    
  } else {
    
    // Add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      SFI_MVLEMR += this->getRayleighDampingForces();
  }
  
  return SFI_MVLEMR;
}
