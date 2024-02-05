// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								California State University, Fullerton 
//								Kutay Orakcal
//								Bogazici University, Istanbul, Turkey
//								John Wallace
//								University of California, Los Angeles
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

#include "MVLEM.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

// initialize the class wide variables
Matrix MVLEM::MVLEMK(6, 6);
Matrix MVLEM::MVLEMM(6, 6);
Matrix MVLEM::MVLEMD(6, 6);
Vector MVLEM::MVLEMR(6);

#include <elementAPI.h>

// Read input parameters and build the material
void *OPS_MVLEM(void)
{
  // Pointer to a uniaxial material that will be returned                       
  Element *theElement = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  // Parse the script for material parameters
  if (numArgs < 7) { 
    opserr << "Want: MVLEM eleTag Dens iNode jNode m c -thick {fiberThick} -width {fiberWidth} -rho {Rho} -matConcrete {matTagsConcrete} -matSteel {matTagsSteel} -matShear {matTagShear}\n";
    return 0;
  }

  int iData[4];
  double dData[2];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for element MVLEM" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid density value for element MVLEM " << iData[0] << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
    opserr << "WARNING iNode jNode or m for element MVLEM" << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "Invalid c value for element MVLEM " << iData[0] << endln;
    return 0;
  }

  int m = iData[3];
  const char *str = 0;

  double *theThickness = new double[m];
  double *theWidth = new double[m];
  double *theRho = new double[m];
  int *matTags = new int[m];

  UniaxialMaterial **theMaterialsConcrete = new UniaxialMaterial* [m];
  UniaxialMaterial **theMaterialsSteel = new UniaxialMaterial*[m];
  UniaxialMaterial **theMaterialsShear = new UniaxialMaterial*[1];
  for (int i = 0; i < m; i++) {
    theThickness[i] = 0.0;
    theWidth[i] = 0.0;
    theRho[i] = 0.0;
    matTags[i] = 0;
    theMaterialsConcrete[i] = 0;
    theMaterialsSteel[i] = 0;
  }
  theMaterialsShear[0] = 0;

  numArgs = OPS_GetNumRemainingInputArgs();
  while (numArgs > 0) {
      //OPS_GetStringCopy(&str);
      str = OPS_GetString();
    if (strcmp(str, "-thick") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theThickness) != 0) {
	opserr << "Invalid thick parameter for MVLEM   " << iData[0] << endln;
	return 0;
      }
    } else if (strcmp(str, "-width") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theWidth) != 0) {
	opserr << "Invalid width value for MVLEM  " << iData[0] << endln;
	return 0;
      }
    } else if (strcmp(str, "-rho") == 0) {
      numData = m;
      if (OPS_GetDoubleInput(&numData, theRho) != 0) {
	opserr << "Invalid rho value for MVLEM  " << iData[0] << endln;
	return 0;
      }
    } else if (strcmp(str, "-matConcrete") == 0) {
      numData = m;
      if (OPS_GetIntInput(&numData, matTags) != 0) {
	opserr << "Invalid concrete tags for MVLEM  " << iData[0] << endln;
	return 0;
      }
      for (int i=0; i<m; i++) {
	theMaterialsConcrete[i] = 0;
	theMaterialsConcrete[i] = OPS_getUniaxialMaterial(matTags[i]);
	if (theMaterialsConcrete[i] == 0) {
	  opserr << "Invalid concrete material tag " << matTags[i] << "  for MVLEM  " << iData[0] << endln;
	  return 0;
	}
      }
    } else if (strcmp(str, "-matSteel") == 0) {
      numData = m;
      if (OPS_GetIntInput(&numData, matTags) != 0) {
	opserr << "Invalid steel tags for MVLEM  " << iData[0] << endln;
	return 0;
      }
      for (int i=0; i<m; i++) {
	theMaterialsSteel[i] = 0;
	theMaterialsSteel[i] = OPS_getUniaxialMaterial(matTags[i]);
	if (theMaterialsSteel[i] == 0) {
	  opserr << "Invalid steel material tag " << matTags[i] << "  for MVLEM  " << iData[0] << endln;
	  return 0;
	}
      }
    } else if (strcmp(str, "-matShear") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, matTags) != 0) {
	opserr << "Invalid shear tag for MVLEM  " << iData[0] << endln;
	return 0;
      }
      for (int i=0; i<1; i++) {
	theMaterialsShear[i] = 0;
	theMaterialsShear[i] = OPS_getUniaxialMaterial(matTags[i]);
	if (theMaterialsShear[i] == 0) {
	  opserr << "Invalid shear material tag " << matTags[i] << "  for MVLEM  " << iData[0] << endln;
	  return 0;
	}
      }
    }
    
    numArgs = OPS_GetNumRemainingInputArgs();
    
  }

  theElement = new MVLEM(iData[0], dData[0], iData[1], iData[2], 
			 theMaterialsConcrete, theMaterialsSteel, theMaterialsShear,
			 theRho, theThickness, theWidth , iData[3], dData[1]);
  
  // Cleanup dynamic memory
  if (theThickness != 0) 
    delete [] theThickness;
  if (theWidth != 0) 
    delete [] theWidth;
  if (theRho != 0) 
    delete [] theRho;
  if (matTags != 0) 
    delete [] matTags;

  if (theMaterialsConcrete != 0)
    delete [] theMaterialsConcrete;
  if (theMaterialsSteel != 0)
    delete[] theMaterialsSteel;
  if (theMaterialsShear != 0)
    delete[] theMaterialsShear;

  return theElement;
}


// typical constructor
MVLEM::MVLEM(int tag,
	double Dens,
	int Nd1, int Nd2,
	UniaxialMaterial **materialsConcrete,
	UniaxialMaterial **materialsSteel,
	UniaxialMaterial **materialsShear,
	double *Rho, 
	double *thickness,
	double *width,
	int mm = 0,
	double cc = 0.0)

	:Element(tag, ELE_TAG_MVLEM),
	density(Dens),
	externalNodes(2),
	theMaterialsConcrete(0), theMaterialsSteel(0), theMaterialsShear(0), 
	theLoad(0), MVLEMStrain(0),
	c(cc), m(mm)
{
  // Fill with ZEROs all element matrices
  MVLEMK.Zero();		// element stiffness matrix
  MVLEMR.Zero();		// element force vector (6)
  
  NodeMass = 0.0;
  h = 0.0;
  
  // Fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 2)
    opserr << "FATAL MVLEM::MVLEM() - out of memory, could not create an ID of size 2\n";
  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;
  
  //Set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;

  // Check thickness and width input
  if (thickness == 0) {
    opserr << "MVLEM::MVLEM() - "
	   << "Null thickness array passed.\n";
    exit(-1);
  }
  
  if (width == 0) {
    opserr << "MVLEM::MVLEM() - "
	   << "Null width array passed.\n";
    exit(-1);
  }
  
  // Allocate memory for element arrays
  // Input parameters
  t = new double[m];			
  b = new double[m];			
  rho = new double[m];		
  Lw = 0.0;					
  
  // Assign values from input
  for (int i = 0; i<m; i++) {
    t[i] = thickness[i];
    b[i] = width[i];
    rho[i] = Rho[i];
    Lw += b[i]; 
  }
  
  // Area of concrete and steel fibers
  Ac = new double[m];
  As = new double[m];
  
  // Fiber strains
  MVLEMStrain = new double[m + 1];
  
  // Assign zero to element arrays
  for (int i = 0; i < m; i++) {
    
    Ac[i] = 0.0;
    As[i] = 0.0;
    
    MVLEMStrain[i] = 0.0;
  }
  
  MVLEMStrain[m] = 0.0;
  
  // Calculate concrete and steel areas in Y directions
  for (int i = 0; i < m; i++) {
    As[i] = (b[i] * t[i])*rho[i]; 
    Ac[i] = (b[i] * t[i]) - As[i]; 
  }
  
  // Calculate locations of concrete macro-fibers in the cross-section (centerline - x = 0.0)
  x = new double[m];
  for (int i = 0; i < m; i++)
    x[i] = 0.0;
  
  for (int i = 0; i < m; i++) {
    double sumb_i = 0.0;
    for (int j = 0; j<i + 1; j++)
      sumb_i += b[j];
    
    x[i] = (sumb_i - b[i] / 2.0) - Lw / 2.0;
  }	
  
  // Determine the nodal mass for lumped mass approach
  A = 0;
  for (int i = 0; i < m; i++){
    A += Ac[i] + As[i];
  }
  
  NodeMass = density * A * h / 2;
  
  // Check Cocnrete material input
  if (materialsConcrete == 0)
    {
      opserr << "MVLEM::MVLEM() - "
	     << "null Concrete material array passed.\n";
      exit(-1);
    }
  
  // Check Steel material input
  if (materialsSteel == 0)
    {
      opserr << "MVLEM::MVLEM() - "
	     << "null Steel material array passed.\n";
      exit(-1);
    }
  
  // Check Shear material input
  if (materialsShear == 0)
    {
      opserr << "MVLEM::MVLEM() - "
	     << "null Shear material passed.\n";
      exit(-1);
    }
  
  // Allocate memory for the Concrete uniaxial materials
  theMaterialsConcrete = new UniaxialMaterial*[m];
  if (theMaterialsConcrete == 0)
    {
      opserr << "MVLEM::MVLEM() - "
	     << "failed to allocate pointers for Concrete uniaxial materials.\n";
      exit(-1);
    }
  
  // Get copies of the Concrete uniaxial materials
  for (int i = 0; i < m; i++)
    {
      if (materialsConcrete[i] == 0)
	{
	  opserr << "MVLEM::MVLEM() - "
	    "null uniaxial Concrete material pointer passed.\n";
	  exit(-1);
	}
      theMaterialsConcrete[i] = materialsConcrete[i]->getCopy();
      if (theMaterialsConcrete[i] == 0)
	{
	  opserr << "MVLEM::MVLEM() - "
		 << "failed to copy Concrete uniaxial material.\n";
	  exit(-1);
	}
    }
  
  // Allocate memory for the Steel uniaxial materials
  theMaterialsSteel = new UniaxialMaterial*[m];
  if (theMaterialsSteel == 0)
    {
      opserr << "MVLEM::MVLEM() - "
	     << "failed to allocate pointers for Steel uniaxial materials.\n";
      exit(-1);
    }
  
  // Get copies of the uniaxial materials
  for (int i = 0; i < m; i++)
    {
      if (materialsSteel[i] == 0)
	{
	  opserr << "MVLEM::MVLEM() - "
	    "null uniaxial Steel material pointer passed.\n";
	  exit(-1);
	}
      theMaterialsSteel[i] = materialsSteel[i]->getCopy();
      if (theMaterialsSteel[i] == 0)
	{
	  opserr << "MVLEM::MVLEM() - "
		 << "failed to copy Steel uniaxial material.\n";
	  exit(-1);
	}
    }
  
  // Allocate memory for the Shear uniaxial materials
  theMaterialsShear = new UniaxialMaterial*[1];
  if (theMaterialsShear == 0)
    {
      opserr << "MVLEM::MVLEM() - "
	     << "failed to allocate pointers for Shear uniaxial materials.\n";
      exit(-1);
    }
  
  // Get copies of the uniaxial materials
  for (int i = 0; i < 1; i++)
    {
      if (materialsShear[i] == 0)
	{
	  opserr << "MVLEM::MVLEM() - "
	    "null uniaxial Shear material pointer passed.\n";
	  exit(-1);
	}
      theMaterialsShear[i] = materialsShear[i]->getCopy();
      if (theMaterialsShear[i] == 0)
	{
	  opserr << "MVLEM::MVLEM() - "
		 << "failed to copy Shear uniaxial material.\n";
	  exit(-1);
	}
    }
  
  // Revert to start
  this->revertToStart();
}

// Constructor which should be invoked by an FE_ObjectBroker only
MVLEM::MVLEM()
  :Element(0, ELE_TAG_MVLEM),
   density(0.0),
   externalNodes(2),
   theMaterialsConcrete(0), theMaterialsSteel(0), theMaterialsShear(0),
   theLoad(0), MVLEMStrain(0),
   b(0), t(0), rho(0), x(0), As(0), Ac(0),
   h(0.0), c(0.0), m(0)
{
  if (externalNodes.Size() != 2)
    opserr << "FATAL MVLEM::MVLEM() - out of memory, could not create an ID of size 2\n";
  theNodes[0] = 0;
  theNodes[1] = 0;
}

//  Destructor - provided to clean up any memory
MVLEM::~MVLEM()
{
  // clean up the memory associated with the element, this is
  // memory the MVLEM objects allocates and memory allocated 
  // by other objects that the MVLEM object is responsible for 
  // cleaning up, i.e. the MaterialObject.
  if (theMaterialsConcrete != 0)  {
    for (int i = 0; i < m; i++)
      if (theMaterialsConcrete[i] != 0)
	delete theMaterialsConcrete[i];
    delete[] theMaterialsConcrete;
  }
  
  if (theMaterialsSteel != 0)  {
    for (int i = 0; i < m; i++)
      if (theMaterialsSteel[i] != 0)
	delete theMaterialsSteel[i];
    delete[] theMaterialsSteel;
  }
  
  if (theMaterialsShear != 0)  {
    for (int i = 0; i < 1; i++)
      if (theMaterialsShear[i] != 0)
	delete theMaterialsShear[i];
    delete[] theMaterialsShear;
  }
  
  if (theLoad != 0)
    delete theLoad;
  
  if (x != 0)
    delete []x;
  if (t != 0)
    delete []t;
  if (b != 0)
    delete []b;
  if (rho != 0)
    delete []rho;
  if (Ac != 0)
    delete []Ac;
  if (As != 0)
    delete []As;
  if (MVLEMStrain != 0)
    delete []MVLEMStrain;
}

int
MVLEM::getNumExternalNodes(void) const
{
	return 2;
}

const ID &
MVLEM::getExternalNodes(void)
{
	return externalNodes;
}

Node **
MVLEM::getNodePtrs(void)
{
	return theNodes;
}

int
MVLEM::getNumDOF(void) {
	return 6;
}

void
MVLEM::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0)
    {
      return;
    }
  theNodes[0] = 0;
  theNodes[1] = 0;
  
  // first ensure nodes exist in Domain and set the node pointers
  int Nd1 = externalNodes(0);
  int Nd2 = externalNodes(1);
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);
  
  if (theNodes[0] == 0)
    {
      opserr << "WARNING MVLEM::setDomain() - at MVLEM " << this->getTag() << " node " <<
	Nd1 << " does not exist in domain\n";
      return;  // don't go any further - otherwise segemntation fault
    }
  
  if (theNodes[1] == 0)
    {
      opserr << "WARNING MVLEM::setDomain() - at MVLEM " << this->getTag() << " node " <<
	Nd2 << " does not exist in domain\n";
      return;
    }

  // Get coordinates of end nodes 
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();
  
  // Calculate the element height and perform checks
  if (end1Crd.Size() != 2 && end2Crd.Size() != 2) {
    opserr << "MVLEM::setDomain(): 2 coords required at nodes, not enough provided for  element " << this->getTag();
    exit(-1);
  }

  h = end2Crd(1) - end1Crd(1);
  
  if (h < 0.0) {
    opserr << "WARNING: Element height is negative. Define Nodes from bottom to top!";
    exit(-1);
  }
  
  if (h == 0.0) {
    opserr << "WARNING: Element height is ZERO!";
    exit(-1);
  }
  
  // Currently element can be only vertical
  if (end2Crd(0) != end1Crd(0)) {
    opserr << "WARNING: Element is NOT vertical!";
    exit(-1);
  }
  
  
  // Call the DomainComponent class method THIS IS VERY IMPORTANT
  this->DomainComponent::setDomain(theDomain);

  // Ensure connected nodes have correct number of dof's
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();
  if ((dofNd1 != 3) || (dofNd2 != 3))
    
    {
      opserr << "MVLEM::setDomain(): 3 dof required at nodes, " << dofNd1 << " and "
	     << dofNd2 << " provided\n";
    }
  
  // Create a vector to hop applied loads
  if (theLoad == 0)
    theLoad = new Vector(6);
  if (theLoad == 0)  {
    opserr << "MVLEM::setDomain() - element: " << this->getTag()
	   << " out of memory creating vector of size: " << 6 << endln;
    return;
  }

}

int
MVLEM::commitState()
{
	int errCode = 0;

	// Commit Concrete material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsConcrete[i]->commitState();

	// Commit Steel material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsSteel[i]->commitState();

	// Commit Shear material models
	for (int i = 0; i < 1; i++)
		errCode += theMaterialsShear[i]->commitState();

	return errCode;
}

int
MVLEM::revertToLastCommit()
{
	int errCode = 0;

	// Revert Concrete material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsConcrete[i]->revertToLastCommit();

	// Revert Steel material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsSteel[i]->revertToLastCommit();

	// Revert Shear material model
	for (int i = 0; i < 1; i++)
		errCode += theMaterialsShear[i]->revertToLastCommit();

	return errCode;
}

int
MVLEM::revertToStart()
{
	int errCode = 0;

	// Revert Concrete material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsConcrete[i]->revertToStart();

	// Revert Steel material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsSteel[i]->revertToStart();

	// Revert Shear material model
	for (int i = 0; i < 1; i++)
		errCode += theMaterialsShear[i]->revertToStart();

	return errCode;
}

int
MVLEM::update()
{
	// Determine the current strain given trial displacements at nodes
	MVLEMStrain = this->computeCurrentStrain();

	// Set the strain in the materials
	int errCode1 = 0;

	// Set trial response for Concrete material models
	for (int i = 0; i < m; i++)
		errCode1 += theMaterialsConcrete[i]->setTrialStrain(MVLEMStrain[i]);

	// Set trial response for Steel material models
	for (int i = 0; i < m; i++)
		errCode1 += theMaterialsSteel[i]->setTrialStrain(MVLEMStrain[i]);

	// Set trial response for Shear material model
		errCode1 += theMaterialsShear[0]->setTrialStrain(MVLEMStrain[m]); 

	return errCode1;
}

const Matrix &
MVLEM::getInitialStiff(void)
{

  double Ec, Es, ky;
  
	// Build the initial stiffness matrix
	double Kv = 0.0; double Kh = 0.0; double Km = 0.0; double e = 0.0; double ex = 0.0;

	for (int i = 0; i < m; ++i)
	{
	  Ec = theMaterialsConcrete[i]->getInitialTangent(); 
	  Es = theMaterialsSteel[i]->getInitialTangent();    
	  ky = Ec * Ac[i] / h + Es * As[i] / h;	  
	  Kv += ky;
	  Km += ky * x[i] * x[i];
	  e += ky * x[i];
	}

	// Get shear stiffness from shear material
	Kh = theMaterialsShear[0]->getInitialTangent();

	// Assemble element stiffness matrix
	MVLEMK(0, 0) = Kh;
	MVLEMK(0, 1) = 0.0;
	MVLEMK(0, 2) = -Kh*c*h;
	MVLEMK(0, 3) = -Kh;
	MVLEMK(0, 4) = 0.0;
	MVLEMK(0, 5) = -Kh*(1 - c)*h;

	MVLEMK(1, 0) = MVLEMK(0, 1);
	MVLEMK(1, 1) = Kv;
	MVLEMK(1, 2) = e;
	MVLEMK(1, 3) = 0.0;
	MVLEMK(1, 4) = -Kv;
	MVLEMK(1, 5) = -e;

	MVLEMK(2, 0) = MVLEMK(0, 2);
	MVLEMK(2, 1) = MVLEMK(1, 2);
	MVLEMK(2, 2) = h*h*c*c*Kh + Km;
	MVLEMK(2, 3) = h*c*Kh;
	MVLEMK(2, 4) = -e;
	MVLEMK(2, 5) = (1 - c)*c*h*h*Kh - Km;

	MVLEMK(3, 0) = MVLEMK(0, 3);
	MVLEMK(3, 1) = MVLEMK(1, 3);
	MVLEMK(3, 2) = MVLEMK(2, 3);
	MVLEMK(3, 3) = Kh;
	MVLEMK(3, 4) = 0.0;
	MVLEMK(3, 5) = Kh*(1 - c)*h;

	MVLEMK(4, 0) = MVLEMK(0, 4);
	MVLEMK(4, 1) = MVLEMK(1, 4);
	MVLEMK(4, 2) = MVLEMK(2, 4);
	MVLEMK(4, 3) = MVLEMK(3, 4);
	MVLEMK(4, 4) = Kv;
	MVLEMK(4, 5) = e;

	MVLEMK(5, 0) = MVLEMK(0, 5);
	MVLEMK(5, 1) = MVLEMK(1, 5);
	MVLEMK(5, 2) = MVLEMK(2, 5);
	MVLEMK(5, 3) = MVLEMK(3, 5);
	MVLEMK(5, 4) = MVLEMK(4, 5);
	MVLEMK(5, 5) = (1 - c)*(1 - c)*h*h*Kh + Km;

	// Return element stiffness matrix
	return MVLEMK;

}

const Matrix &
MVLEM::getTangentStiff(void)
{
  double Ec, Es, ky;

	// Build the initial stiffness matrix
	double Kv = 0.0; double Kh = 0.0; double Km = 0.0; double e = 0.0; double ex = 0.0;

	for (int i = 0; i < m; ++i)
	{
		Ec = theMaterialsConcrete[i]->getTangent();
		Es = theMaterialsSteel[i]->getTangent();
		ky = Ec * Ac[i] / h + Es * As[i] / h;	  
		Kv += ky;
		Km += ky * x[i] * x[i];
		e += ky * x[i];

	}

	// Get shear stiffness from shear material
	Kh = theMaterialsShear[0]->getTangent();

	// Assemble element stiffness matrix
	MVLEMK(0, 0) = Kh;
	MVLEMK(0, 1) = 0.0;
	MVLEMK(0, 2) = -Kh*c*h;
	MVLEMK(0, 3) = -Kh;
	MVLEMK(0, 4) = 0.0;
	MVLEMK(0, 5) = -Kh*(1 - c)*h;

	MVLEMK(1, 0) = MVLEMK(0, 1);
	MVLEMK(1, 1) = Kv;
	MVLEMK(1, 2) = e;
	MVLEMK(1, 3) = 0.0;
	MVLEMK(1, 4) = -Kv;
	MVLEMK(1, 5) = -e;

	MVLEMK(2, 0) = MVLEMK(0, 2);
	MVLEMK(2, 1) = MVLEMK(1, 2);
	MVLEMK(2, 2) = h*h*c*c*Kh + Km;
	MVLEMK(2, 3) = h*c*Kh;
	MVLEMK(2, 4) = -e;
	MVLEMK(2, 5) = (1 - c)*c*h*h*Kh - Km;

	MVLEMK(3, 0) = MVLEMK(0, 3);
	MVLEMK(3, 1) = MVLEMK(1, 3);
	MVLEMK(3, 2) = MVLEMK(2, 3);
	MVLEMK(3, 3) = Kh;
	MVLEMK(3, 4) = 0.0;
	MVLEMK(3, 5) = Kh*(1 - c)*h;

	MVLEMK(4, 0) = MVLEMK(0, 4);
	MVLEMK(4, 1) = MVLEMK(1, 4);
	MVLEMK(4, 2) = MVLEMK(2, 4);
	MVLEMK(4, 3) = MVLEMK(3, 4);
	MVLEMK(4, 4) = Kv;
	MVLEMK(4, 5) = e;

	MVLEMK(5, 0) = MVLEMK(0, 5);
	MVLEMK(5, 1) = MVLEMK(1, 5);
	MVLEMK(5, 2) = MVLEMK(2, 5);
	MVLEMK(5, 3) = MVLEMK(3, 5);
	MVLEMK(5, 4) = MVLEMK(4, 5);
	MVLEMK(5, 5) = (1 - c)*(1 - c)*h*h*Kh + Km;

	// Return element stiffness matrix
	return MVLEMK;

}

// Get element mass matrix assuming lumped mass
const Matrix & MVLEM::getMass(void)
{
	MVLEMM.Zero();

	// No rotational mass
	MVLEMM(0, 0) = NodeMass;
	MVLEMM(1, 1) = NodeMass;
	MVLEMM(3, 3) = NodeMass;
	MVLEMM(4, 4) = NodeMass;

	// Return element mass matrix
	return MVLEMM;
}

// Get element damping matrix
const Matrix & MVLEM::getDamp(void)
{
	MVLEMD.Zero();

	MVLEMD = this->Element::getDamp();

	// Return element damping matrix
	return MVLEMD;
}
void
MVLEM::zeroLoad(void)
{
	// does nothing - no elemental loads
}

int
MVLEM::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	return 0;
}

int
MVLEM::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

// Get element force vector
const Vector & MVLEM::getResistingForce()
{
	
	MVLEMR.Zero();

	// Get Trial Displacements
	///	const Vector &disp1 = theNodes[0]->getTrialDisp();
	///	const Vector &disp2 = theNodes[1]->getTrialDisp();

	MVLEMR(0) = theMaterialsShear[0]->getStress(); // get force from shear force-deformation relationship

	double stressC, stressS;
	for (int i = 0; i<m; i++)
	{
		stressC = theMaterialsConcrete[i]->getStress();
		stressS = theMaterialsSteel[i]->getStress();
		MVLEMR(1) += -stressC * Ac[i] - stressS * As[i];
		MVLEMR(2) += -stressC * Ac[i] * x[i] - stressS * As[i] * x[i];
		MVLEMR(5) += stressC * Ac[i] * x[i] + stressS * As[i] * x[i];
	}

	MVLEMR(2) += -MVLEMR(0)*c*h;
	MVLEMR(3) = -MVLEMR(0);
	MVLEMR(4) = -MVLEMR(1);
	MVLEMR(5) += -MVLEMR(0)*(1 - c)*h;

	// Return element force vector
	return MVLEMR;
}


const Vector & MVLEM::getResistingForceIncInertia()
{
	this->getResistingForce();

	if (NodeMass != 0.0) { 
		const Vector &accel1 = theNodes[0]->getTrialAccel();
		const Vector &accel2 = theNodes[1]->getTrialAccel();

		// Compute the current resisting force
		this->getResistingForce();

		MVLEMR(0) += NodeMass*accel1(0);
		MVLEMR(1) += NodeMass*accel1(1);
		MVLEMR(3) += NodeMass*accel2(0);
		MVLEMR(4) += NodeMass*accel2(1);

		// Add the damping forces if rayleigh damping
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			MVLEMR += this->getRayleighDampingForces();

	}
	else {

		// Add the damping forces if rayleigh damping
		if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			MVLEMR += this->getRayleighDampingForces();
	}

	return MVLEMR;
}

int
MVLEM::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;
	int dataTag = this->getDbTag();

	ID idData0(4);
	idData0(0) = externalNodes(0);
	idData0(1) = externalNodes(1);	
	idData0(2) = this->getTag();
	idData0(3) = m;

	res = theChannel.sendID(dataTag, commitTag, idData0);


	
	int matDbTag;
	// Send the connected nodes (2) and material class/db tags (4m flex, 2 shear)
	ID idData(2 + 4*m);
	for (int i = 0; i < m; i++) {
	  idData(i) = theMaterialsConcrete[i]->getClassTag();
	  matDbTag = theMaterialsConcrete[i]->getDbTag();
	  if (matDbTag == 0) {
	    matDbTag = theChannel.getDbTag();
	    if (matDbTag != 0)
	      theMaterialsConcrete[i]->setDbTag(matDbTag);
	  }
	  idData(i+m) = matDbTag;

	  idData(i+2*m) = theMaterialsSteel[i]->getClassTag();
	  matDbTag = theMaterialsSteel[i]->getDbTag();
	  if (matDbTag == 0) {
	    matDbTag = theChannel.getDbTag();
	    if (matDbTag != 0)
	      theMaterialsSteel[i]->setDbTag(matDbTag);
	  }	  
	  idData(i+3*m) = matDbTag;	  
	}
	idData(4*m) = theMaterialsShear[0]->getClassTag();
	matDbTag = theMaterialsShear[0]->getDbTag();
	if (matDbTag == 0) {
	  matDbTag = theChannel.getDbTag();
	  if (matDbTag != 0)
	    theMaterialsShear[0]->setDbTag(matDbTag);
	}	  
	idData(4*m+1) = matDbTag;	  	
	
	res = theChannel.sendID(dataTag, commitTag, idData);

	
	Vector data(3 + 3*m);

	data(3*m) = density;
	data(3*m+1) = c;
	data(3*m+2) = h;	
	for (int i = 0; i < m; i++) {
	  data(i) = b[i];
	  data(i+m) = t[i];
	  data(i+2*m) = rho[i];
	}

	// MVLEM then sends the tags of it's two end nodes
	res = theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
	  opserr << "WARNING MVLEM::sendSelf() - failed to send ID\n";
	  return -2;
	}

	
	// Send the material models
	for (int i = 0; i < m; i++) {
	  res += theMaterialsConcrete[i]->sendSelf(commitTag, theChannel);
	  if (res < 0) {
	    opserr << "WARNING MVLEM::sendSelf - " << this->getTag() << " failed to send concrete material\n";
	    return res;
	  }
	}
	
	for (int i = 0; i < m; i++) {	
	  theMaterialsSteel[i]->sendSelf(commitTag, theChannel);
	  if (res < 0) {
	    opserr << "WARNING MVLEM::sendSelf - " << this->getTag() << " failed to send steel material\n";
	    return res;
	  }	  
	}
	
	res += theMaterialsShear[0]->sendSelf(commitTag, theChannel);
	if (res < 0) {
	  opserr << "WARNING MVLEM::sendSelf - " << this->getTag() << " failed to send shear material\n";
	  return res;
	}

	return res;	
}

int
MVLEM::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	int dataTag = this->getDbTag();

	ID idData0(4);
	
	res = theChannel.recvID(dataTag, commitTag, idData0);
	if (res < 0) {
	  opserr << "WARNING MVLEM::recvSelf() - failed to receive ID\n";
	  return -2;
	}
	externalNodes(0) = idData0(0);
	externalNodes(1) = idData0(1);	
	this->setTag(idData0(2));
	m = idData0(3);
	
	
	ID idData(2 + 4*m);
	
	// MVLEM now receives the tags of it's two external nodes
	res = theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
	  opserr << "WARNING MVLEM::recvSelf() - failed to receive ID\n";
	  return -2;
	}

	
	Vector data(3 + 3*m);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
	  opserr << "WARNING MVLEM::recvSelf() - failed to receive Vector\n";
	  return res;
	}

	density = data(3*m);
	c = data(3*m+1);
	h = data(3*m+2);
	
	if (theMaterialsConcrete == 0) {
	  // Allocate new materials
	  theMaterialsConcrete = new UniaxialMaterial *[m];
	  if (theMaterialsConcrete == 0) {
	    opserr << "MVLEM::recvSelf - could not allocated UniaxialMaterial array for concrete\n";
	    return -1;
	  }

	  // Receive the Concrete material models
	  for (int i = 0; i < m; i++)  {
	    int matClassTag = idData(i);
	    int matDbTag = idData(i+m);
	    theMaterialsConcrete[i] = theBroker.getNewUniaxialMaterial(matClassTag);
	    if (theMaterialsConcrete[i] == 0) {
	      opserr << "MVLEM::recvSelf() - "
		     << "broker could not create concrete uniaxial material.\n";
	      return -3;
	    }
	    theMaterialsConcrete[i]->setDbTag(matDbTag);
	    res = theMaterialsConcrete[i]->recvSelf(commitTag, theChannel, theBroker);
	    if (res < 0) {
	      opserr << "MVLEM::recvSelf() - cocnrete material " << i << " failed to recvSelf\n";
	      return res;
	    }
	  }
	}
	else {
	  // Receive the Concrete material models
	  for (int i = 0; i < m; i++)  {
	    int matClassTag = idData(i);
	    int matDbTag = idData(i+m);
	    if (theMaterialsConcrete[i]->getClassTag() != matClassTag) {
	      delete theMaterialsConcrete[i];
	      theMaterialsConcrete[i] = theBroker.getNewUniaxialMaterial(matClassTag);
	      if (theMaterialsConcrete[i] == 0) {
		opserr << "MVLEM::recvSelf() - "
		       << "broker could not create concrete uniaxial material.\n";
		return -3;
	      }
	    }
	    theMaterialsConcrete[i]->setDbTag(matDbTag);
	    res = theMaterialsConcrete[i]->recvSelf(commitTag, theChannel, theBroker);
	    if (res < 0) {
	      opserr << "MVLEM::recvSelf() - cocnrete material " << i << " failed to recvSelf\n";
	      return res;
	    }
	  }
	}

	if (theMaterialsSteel == 0) {
	  // Allocate new materials
	  theMaterialsSteel = new UniaxialMaterial *[m];
	  if (theMaterialsSteel == 0) {
	    opserr << "MVLEM::recvSelf - could not allocated UniaxialMaterial array for steel\n";
	    return -1;
	  }
	  // Receive the steel material models
	  for (int i = 0; i < m; i++)  {
	    int matClassTag = idData(i+2*m);
	    int matDbTag = idData(i+3*m);
	    theMaterialsSteel[i] = theBroker.getNewUniaxialMaterial(matClassTag);
	    if (theMaterialsSteel[i] == 0) {
	      opserr << "MVLEM::recvSelf() - "
		     << "broker could not create steel uniaxial material.\n";
	      return -3;
	    }
	    theMaterialsSteel[i]->setDbTag(matDbTag);
	    res = theMaterialsSteel[i]->recvSelf(commitTag, theChannel, theBroker);
	    if (res < 0) {
	      opserr << "MVLEM::recvSelf() - steel material " << i << " failed to recvSelf\n";
	      return res;
	    }
	  }
	}
	else {
	  // Receive the steel material models
	  for (int i = 0; i < m; i++)  {
	    int matClassTag = idData(i+2*m);
	    int matDbTag = idData(i+3*m);
	    if (theMaterialsSteel[i]->getClassTag() != matClassTag) {
	      delete theMaterialsSteel[i];
	      theMaterialsSteel[i] = theBroker.getNewUniaxialMaterial(matClassTag);
	      if (theMaterialsSteel[i] == 0) {
		opserr << "MVLEM::recvSelf() - "
		       << "broker could not create steel uniaxial material.\n";
		return -3;
	      }
	    }
	    theMaterialsSteel[i]->setDbTag(matDbTag);
	    res = theMaterialsSteel[i]->recvSelf(commitTag, theChannel, theBroker);
	    if (res < 0) {
	      opserr << "MVLEM::recvSelf() - steel material " << i << " failed to recvSelf\n";
	      return res;
	    }
	  }
	}	

	
	if (theMaterialsShear == 0) {
	  // Allocate new materials
	  theMaterialsShear = new UniaxialMaterial *[1];
	  if (theMaterialsShear == 0) {
	    opserr << "MVLEM::recvSelf - could not allocated UniaxialMaterial array for shear\n";
	    return -1;
	  }

	  // Receive the shear material model
	  for (int i = 0; i < 1; i++)  {
	    int matClassTag = idData(4*m);
	    int matDbTag = idData(1+4*m);
	    theMaterialsShear[i] = theBroker.getNewUniaxialMaterial(matClassTag);
	    if (theMaterialsShear[i] == 0) {
	      opserr << "MVLEM::recvSelf() - "
		     << "broker could not create shear uniaxial material.\n";
	      return -3;
	    }
	    theMaterialsShear[i]->setDbTag(matDbTag);
	    res = theMaterialsShear[i]->recvSelf(commitTag, theChannel, theBroker);
	    if (res < 0) {
	      opserr << "MVLEM::recvSelf() - shear material " << i << " failed to recvSelf\n";
	      return res;
	    }
	  }
	}
	else {
	  // Receive the shear material models
	  for (int i = 0; i < 1; i++)  {
	    int matClassTag = idData(4*m);
	    int matDbTag = idData(1+4*m);
	    if (theMaterialsShear[i]->getClassTag() != matClassTag) {
	      delete theMaterialsShear[i];
	      theMaterialsShear[i] = theBroker.getNewUniaxialMaterial(matClassTag);
	      if (theMaterialsShear[i] == 0) {
		opserr << "MVLEM::recvSelf() - "
		       << "broker could not create shear uniaxial material.\n";
		return -3;
	      }
	    }
	    theMaterialsShear[i]->setDbTag(matDbTag);
	    res = theMaterialsShear[i]->recvSelf(commitTag, theChannel, theBroker);
	    if (res < 0) {
	      opserr << "MVLEM::recvSelf() - shear material " << i << " failed to recvSelf\n";
	      return res;
	    }
	  }

	}	

	if (b != 0)
	  delete [] b;
	b = new double[m];
	
	if (t != 0)
	  delete [] t;
	t = new double[m];
	
	if (rho != 0)
	  delete [] rho;
	rho = new double[m];	

	Lw = 0.0;
	for (int i = 0; i < m; i++) {
	  b[i] = data(i);
	  t[i] = data(i+m);
	  rho[i] = data(i+2*m);
	  Lw += b[i];
	}
  
	if (x != 0)
	  delete [] x;
	x = new double[m];
	
	if (Ac != 0)
	  delete [] Ac;
	Ac = new double[m];
	
	if (As != 0)
	  delete [] As;
	As = new double[m];
	
	if (MVLEMStrain != 0)
	  delete [] MVLEMStrain;
	MVLEMStrain = new double[m+1];		
	
	this->setupMacroFibers();
	
	return res;
}

// Display model
int MVLEM::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
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
	error += theViewer.drawLine(v1, v2, RGB, RGB, 1, 1);

	// Displaying Panels
	for (int panel = 0; panel < m; panel++) // loop over m panels
	{
		Matrix NodePLotCrds(m, 13); // (panel id, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)

		// first set the quantity to be displayed at the nodes;
		// if displayMode is 1 through 3 we will plot material stresses otherwise 0.0
		static Vector values(1); // values to be plotted (either epsX, epsY, gammaXY)
		if (displayMode < 4 && displayMode > 0) {
			values(0) = theMaterialsConcrete[panel]->getStrain();
		}
		else {
			values(0) = 0.0;
		}

		// Fiber nodes
		NodePLotCrds(panel, 0) = panel + 1; // panel id
		// Local node 1 - bottom left
		NodePLotCrds(panel, 1) = v1(0) + x[panel] - b[panel] / 2.0; // x 
		NodePLotCrds(panel, 2) = v1(1) + (x[panel] - b[panel] / 2.0) * r1(0); // y
		NodePLotCrds(panel, 3) = v1(2); // z
		// Local node 2 - bottom right
		NodePLotCrds(panel, 4) = v1(0) + x[panel] + b[panel] / 2.0; // x
		NodePLotCrds(panel, 5) = v1(1) + (x[panel] + b[panel] / 2.0) * r1(0); // y
		NodePLotCrds(panel, 6) = v1(2); // z
		// Local node 3 - top left
		NodePLotCrds(panel, 7) = v2(0) + x[panel] + b[panel] / 2.0; // x
		NodePLotCrds(panel, 8) = v2(1) + (x[panel] + b[panel] / 2.0) * r1(0); // y
		NodePLotCrds(panel, 9) = v2(2); // z
		// Local node 4 - top right
		NodePLotCrds(panel, 10) = v2(0) + x[panel] - b[panel] / 2.0; // x
		NodePLotCrds(panel, 11) = v2(1) + (x[panel] - b[panel] / 2.0) * r1(0); // y
		NodePLotCrds(panel, 12) = v2(2); // z

		Matrix coords(4, 3); // Temporary coordinates for plotting

		coords(0, 0) = NodePLotCrds(panel, 1); // node 1 x
		coords(1, 0) = NodePLotCrds(panel, 4); // node 2 x
		coords(2, 0) = NodePLotCrds(panel, 7); // node 3 x
		coords(3, 0) = NodePLotCrds(panel, 10); // node 4 x

		coords(0, 1) = NodePLotCrds(panel, 2); // node 1 y
		coords(1, 1) = NodePLotCrds(panel, 5); // node 2 y
		coords(2, 1) = NodePLotCrds(panel, 8); // node 3 y
		coords(3, 1) = NodePLotCrds(panel, 11); // node 4 y

		coords(0, 2) = NodePLotCrds(panel, 3); // node 1 z
		coords(1, 2) = NodePLotCrds(panel, 6); // node 2 z
		coords(2, 2) = NodePLotCrds(panel, 9); // node 3 z
		coords(3, 2) = NodePLotCrds(panel, 12); // node 4 z

		error += theViewer.drawPolygon(coords, values);
	}

	return error;
}

void
MVLEM::Print(OPS_Stream &s, int flag)
{
	if (flag == 0)
	{
		// Print out element properties
		s << "Element: " << this->getTag() << endln;
		s << "  type: MVLEM" << endln;
		s << "  iNode: " << externalNodes(0) << ", jNode: " << externalNodes(1) << endln;
		s << "Element height: " << h << endln;
		s << "Number of uniaxial fibers elements: " << m << endln << endln;

		// determine resisting forces in global system
		s << "  Global resisting force: " << this->getResistingForce() << endln << endln;

		s << "Fiber responses: " << endln;

		for (int i = 0; i < m; i++)
		{
			s << "Fiber #: " << i+1 << endln;
			s << "Concrete material with tag: " << theMaterialsConcrete[i]->getTag() << endln;
			theMaterialsConcrete[i]->Print(s, flag);

			s << "Steel material with tag: " << theMaterialsSteel[i]->getTag() << endln;
			theMaterialsSteel[i]->Print(s, flag);

		}

		s << "Shear material with tag: " << theMaterialsShear[0]->getTag() << endln;
		theMaterialsShear[0]->Print(s, flag);

	}
	else if (flag == 1)
	{        // does nothing
	}
}

Response *MVLEM::setResponse(const char **argv, int argc, OPS_Stream &s)
{
	Response *theResponse = 0;

	s.tag("ElementOutput");
	s.attr("eleType", "MVLEM");
	s.attr("eleTag", this->getTag());
	s.attr("node1", externalNodes[0]);
	s.attr("node2", externalNodes[1]);

	// Global forces
	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

		s.tag("ResponseType", "Fx_i");
		s.tag("ResponseType", "Fy_i");
		s.tag("ResponseType", "Mz_i");
		s.tag("ResponseType", "Fx_j");
		s.tag("ResponseType", "Fy_j");
		s.tag("ResponseType", "Mz_j");

		return new ElementResponse(this, 1, Vector(6));

	}

	// Element curvature
	else if (strcmp(argv[0], "Curvature") == 0 || strcmp(argv[0], "curvature") == 0) {

		s.tag("ResponseType", "fi");

		return new ElementResponse(this, 2, 0.0);
	}

	// Fiber Strain
	else if (strcmp(argv[0], "Fiber_Strain") == 0 || strcmp(argv[0], "fiber_strain") == 0) {

		s.tag("ResponseType", "ey");

		return new ElementResponse(this, 3, Vector(m));
	}

	// Fiber Concrete Stress
	else if (strcmp(argv[0], "Fiber_Stress_Concrete") == 0 || strcmp(argv[0], "fiber_stress_concrete") == 0) {

		s.tag("ResponseType", "syc");

		return new ElementResponse(this, 4, Vector(m));
	}

	// Fiber Steel Stress
	else if (strcmp(argv[0], "Fiber_Stress_Steel") == 0 || strcmp(argv[0], "fiber_stress_steel") == 0) {

		s.tag("ResponseType", "sys");

		return new ElementResponse(this, 5, Vector(m));
	}

	// Shear Force Deformation
	else if (strcmp(argv[0], "Shear_Force_Deformation") == 0 || strcmp(argv[0], "shear_force_deformation") == 0) {

		s.tag("ResponseType", "shFD");

		return new ElementResponse(this, 6, Vector(2));
	}

	s.endTag();

	return 0;
}

int
MVLEM::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID)
	{
	case 1:  // Global forces
		return eleInfo.setVector(this->getResistingForce());

	case 2:  // Curvature
		return eleInfo.setDouble(this->getCurvature());

	case 3:  // Fiber Strains
		return eleInfo.setVector(this->getStrain());

	case 4:  // Fiber Concrete Stress
		return eleInfo.setVector(this->getStressConcrete());

	case 5:  // Fiber Steel Stress
		return eleInfo.setVector(this->getStressSteel());  

	case 6:  // Shear Force-Deformtion
		return eleInfo.setVector(this->getShearFD()); 

	default:

		return 0;

	}
}

double *
MVLEM::computeCurrentStrain(void)
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();

	MVLEMStrain[m] = disp1(0) - disp2(0) - c*h*disp1(2) - (1 - c)*h*disp2(2); // Shear deformation

	for (int i = 0; i < m; i++) {
		MVLEMStrain[i] = (-disp1(1) - x[i] * disp1(2) + disp2(1) + x[i] * disp2(2)) / h;	// Fiber (Flexural) Strains
	}

	return MVLEMStrain;

}

// Get curvature (from vertical strains)
double MVLEM::getCurvature(void)
{
	double Curv;

	Curv = (MVLEMStrain[0] - MVLEMStrain[m - 1]) / (x[0] - x[m - 1]);

	return Curv;
}

// Get fiber strains
Vector MVLEM::getStrain(void)
{
	Vector fiberStrain(m);

	for (int i = 0; i<m; i++) {
		fiberStrain(i) = MVLEMStrain[i];
	}

	return fiberStrain;
}

// Get Concrete Stress 
Vector MVLEM::getStressConcrete(void)
{
	Vector concreteStress(m);

	for (int i = 0; i<m; i++) {
		concreteStress(i) = theMaterialsConcrete[i]->getStress();
	}

	return concreteStress;
}

// Get Steel Stress 
Vector MVLEM::getStressSteel(void)
{
	Vector steelStress(m);

	for (int i = 0; i<m; i++) {
		steelStress(i) = theMaterialsSteel[i]->getStress();
	}

	return steelStress;
}

// Get Shear Stress-Strain 
Vector MVLEM::getShearFD(void)
{
	Vector shearStrainStress(2);

	shearStrainStress(0) = theMaterialsShear[0]->getStrain();
	shearStrainStress(1) = theMaterialsShear[0]->getStress();

	return shearStrainStress;
}

int
MVLEM::setupMacroFibers()
{
  // Calculate concrete and steel areas in Y directions
  for (int i = 0; i < m; i++) {
    As[i] = (b[i] * t[i])*rho[i]; 
    Ac[i] = (b[i] * t[i]) - As[i]; 
  }
  
  for (int i = 0; i < m; i++) {
    double sumb_i = 0.0;
    for (int j = 0; j<i + 1; j++)
      sumb_i += b[j];
    
    x[i] = (sumb_i - b[i] / 2.0) - Lw / 2.0;
  }	
  
  // Determine the nodal mass for lumped mass approach
  A = 0;
  for (int i = 0; i < m; i++){
    A += Ac[i] + As[i];
  }

  NodeMass = density * A * h / 2;

  return 0;
}
