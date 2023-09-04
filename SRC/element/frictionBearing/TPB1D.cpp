/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: Troy/Fenz/Misra
// C++ Conversion by: fmk
// Created: 01/11

// Description: This file contains the implementation for the TPB1D class.
//
// What: "@(#) TPB1D.C, revA"

#include "TPB1D.h"
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <ElasticMaterial.h>
#include <ElasticPPMaterial.h>
#include <ParallelMaterial.h>
#include <EPPGapMaterial.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>

#include <ElementResponse.h>

// initialise the class wide variables
Matrix TPB1D::TPB1DM2(2,2);
Matrix TPB1D::TPB1DM4(4,4);
Matrix TPB1D::TPB1DM6(6,6);
Matrix TPB1D::TPB1DM12(12,12);
Vector TPB1D::TPB1DV2(2);
Vector TPB1D::TPB1DV4(4);
Vector TPB1D::TPB1DV6(6);
Vector TPB1D::TPB1DV12(12);

//  Constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the TPB1D end nodes.

static int numMyTPB1D = 0;

void *
OPS_TPB1D()
{
  // print out a message about who wrote this element & any copyright info wanted
  if (numMyTPB1D == 0) {
    opserr << "TPB1D2D element - Written by Troy/Fenz UC Berkeley Copyright 2011 - Use at your Own Peril\n";
    numMyTPB1D++;
  }

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new TPB1D();
    return theEle;
  }

  if (numRemainingArgs != 20) {
    opserr << "ERROR - TPB1D2D not enough args provided, want: element TPB1D2D tag? iNode? jNode? direction? mu1? mu2? mu3? R1? R2? R3? h1? h2? h3? D1? D2? D3? d1? d2? d3? W?\n";
    numMyTPB1D++;
  }

  // get the id and end nodes 
  int iData[4];
  double dData[16];
  int numData;

  numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 16;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element area for element" << eleTag << endln;
    return 0;
  }

  // now create the truss and add it to the Domain
  theEle = new TPB1D(eleTag, iData[1], iData[2], iData[3]-1, 
		     &dData[0],
		     &dData[3],
		     &dData[6],
		     &dData[9],
		     &dData[12],
		     dData[15]);

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element type TPB1D with tag " << eleTag << endln;
    return 0;
  }

  return theEle;
}


TPB1D::TPB1D(int tag,
	     int Nd1, 
	     int Nd2, 
	     int dir,
	     double *muI,
	     double *RI,
	     double *hI,
	     double *DI,
	     double *dI,
	     double WI)
:Element(tag,ELE_TAG_TPB1D),     
 connectedExternalNodes(2),
 numDOF(0), 
 direction(dir),
 theMatrix(0), 
 theVector(0),
 theMaterial(0), 
 d0(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;

  for (int i=0; i<3; i++) {
    mu[i] = muI[i];
    R[i] = RI[i];
    h[i] = hI[i];
    D[i] = DI[i];
    d[i] = dI[i];
  }
  W = WI;

  double mu1 = mu[0];
  double mu2 = mu[1];
  double mu3 = mu[2];

  double R1 = R[0];
  double h1 = h[0];
  double R2 = R[1];
  double h2 = h[1];
  double R3 = R[2];
  double h3 = h[2];

  double D1 = D[0];
  double d1 = d[0];
  double D2 = D[1];
  double d2 = d[1];
  double D3 = D[2];
  double d3 = d[2];

  double L1 = R1-h1;
  double L2 = R2-h2;
  double L3 = R3-h3;

  double u1bar = L1/R1*(D1-d1)/2.0;
  double u2bar = L2/R2*(D2-d2)/2.0;
  double u3bar = L3/R3*(D3-d3)/2.0;

  //Define spring materials
  double f1 = mu1;
  double f2 = mu2;
  double f3 = mu3;
  
  double k11 = W*mu1;
  //set k21 [expr $W*$mu2];
  //set k31 [expr $W*$mu3];

  double k52 = 1.0/(2*L1);
  double k42 = -k52+1.0/(L1+L3);
  double k32 = -1.0/(L1+L3)+1.0/(L2+L3);
  double k22 = -1.0/(L2+L3)+1.0/(L1+L2);
  double k12 = -1.0/(L1+L2)+1.0/(2.0*L1);
  double k02 = -1.0/(2.0*L1)+k11;
  
  double epsyP1 = f1/k11;
  double epsyP2 = 2*L1*(f2-f1);
  double epsyP3 = L1*(f2+f3-2*f1)+L2*(f3-f2);
  double epsyP4 = epsyP3+(u2bar/L2+f2-f3)*(L2+L3);
  double epsyP5 = epsyP4+(u3bar/L3-u2bar/L2+f3-f2)*(L1+L3);
  
  double gapyield = 50*W;

  UniaxialMaterial *theMaterials[10];
  
  int damage = 0;
  double eta = 0.0;

  //uniaxialMaterial ElasticPP 11 k02 epsyP1;
  theMaterials[0] = new ElasticPPMaterial(11, k02, epsyP1);
  
  //uniaxialMaterial ElasticPP 12  W*k12 epsyP2;
  theMaterials[1] = new ElasticPPMaterial(12, W*k12, epsyP2);

  //  uniaxialMaterial ElasticPP 13  W*k22 epsyP3;
  theMaterials[2] = new ElasticPPMaterial(13, W*k22, epsyP3);

  //uniaxialMaterial Elastic 14 W*k32;
  theMaterials[3] = new ElasticMaterial(14, W*k32);

  //  uniaxialMaterial ElasticPPGap 15  -W*k32 gapyield epsyP4;
  theMaterials[4] = new EPPGapMaterial(15, -W*k32, gapyield, epsyP4, eta, damage);

  //  uniaxialMaterial ElasticPPGap 16  -W*k32 -gapyield -epsyP4;
  theMaterials[5] = new EPPGapMaterial(16, -W*k32, -gapyield, -epsyP4, eta, damage);

  //  uniaxialMaterial Elastic 17  W*k42;
  theMaterials[6] = new ElasticMaterial(17, W*k42);

  //   uniaxialMaterial ElasticPPGap 18 -W*k42 gapyield epsyP5;
  theMaterials[7] = new EPPGapMaterial(18, -W*k42, gapyield, epsyP5, eta, damage);

  //  uniaxialMaterial ElasticPPGap 19 -W*k42 -gapyield -epsyP5;
  theMaterials[8] = new EPPGapMaterial(19, -W*k42, -gapyield, -epsyP5, eta, damage);

  //  uniaxialMaterial ElasticPP 20  W*k52  1000*epsyP5;
  theMaterials[9] = new ElasticPPMaterial(20, W*k52, 1000*epsyP5);
  
  //  uniaxialMaterial Parallel 100 11 12 13 14 15 16 17 18 19 20;
  theMaterial = new ParallelMaterial(1, 10, theMaterials);

  for (int i = 0; i < 10; i++)
	  delete theMaterials[i];
}


//   Constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
TPB1D::TPB1D(void)
  :Element(0,ELE_TAG_TPB1D),     
  connectedExternalNodes(2),
  dimension(0), numDOF(0), 
  theMatrix(0), theVector(0),
  theMaterial(0)
{

}


//  Destructor:
//  delete must be invoked on any objects created by the object
//  and on the matertial object.
TPB1D::~TPB1D()
{
  delete theMaterial;

  if (d0 != 0)
    delete d0;
}


int
TPB1D::getNumExternalNodes(void) const
{
    return 2;
}


const ID &
TPB1D::getExternalNodes(void) 
{
    return connectedExternalNodes;
}



Node **
TPB1D::getNodePtrs(void) 
{
  return theNodes;
}

int
TPB1D::getNumDOF(void) 
{
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the TPB1D element, we set matrix and vector pointers,
//    allocate space for t matrix and define it as the basic deformation-
//    displacement transformation matrix.
void
TPB1D::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    return;
  }
  
  // set default values for error conditions
  numDOF = 2;
  theMatrix = &TPB1DM2;
  theVector = &TPB1DV2;
  
  // first set the node pointers
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);	
  
  // if can't find both - send a warning message
  if ( theNodes[0] == 0 || theNodes[1] == 0 ) {
    if (theNodes[0] == 0) 
      opserr << "WARNING TPB1D::setDomain() - Nd1: " << Nd1 << " does not exist in ";
    else
      opserr << "WARNING TPB1D::setDomain() - Nd2: " << Nd2 << " does not exist in ";
    
    opserr << "model for TPB1D ele: " << this->getTag() << endln;
    
    return;
  }
  
  // now determine the number of dof and the dimension    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if ( dofNd1 != dofNd2 ) {
      opserr << "WARNING TPB1D::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
	"have differing dof at ends for TPB1D " << this->getTag() << endln;
      return;
    }	

    // Check that length is zero within tolerance
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    dimension = end1Crd.Size();

    Vector diff = end1Crd - end2Crd;
    double L  = diff.Norm();
    double v1 = end1Crd.Norm();
    double v2 = end2Crd.Norm();
    double vm;

    vm = (v1<v2) ? v2 : v1;

    double LENTOL = 1.0e-12;

    if (L > LENTOL*vm)
      opserr << "WARNING TPB1D::setDomain(): Element " << this->getTag() << " has L= " << L << 
	", which is greater than the tolerance\n";
        
    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    if (direction < 0) 
      direction = -direction;
    
    // set the number of dof for element and set matrix and vector pointer
    if (dimension == 1 && dofNd1 == 1 && direction == 0) {
	numDOF = 2;    
	theMatrix = &TPB1DM2;
	theVector = &TPB1DV2;
    }
    else if (dimension == 2 && dofNd1 == 2 && direction < 2) {
	numDOF = 4;
	theMatrix = &TPB1DM4;
	theVector = &TPB1DV4;
    }
    else if (dimension == 2 && dofNd1 == 3 && direction < 3) {
	numDOF = 6;	
	theMatrix = &TPB1DM6;
	theVector = &TPB1DV6;
    }
    else if (dimension == 3 && dofNd1 == 3 && direction < 3) {
	numDOF = 6;	
	theMatrix = &TPB1DM6;
	theVector = &TPB1DV6;
    }
    else if (dimension == 3 && dofNd1 == 6 && direction< 6) {
	numDOF = 12;	    
	theMatrix = &TPB1DM12;
	theVector = &TPB1DV12;
    }
    else {
      opserr << "WARNING TPB1D::setDomain cannot handle " << dimension << 
	"dofs at nodes in " << dofNd1 << " d problem\n"; 
      return;
    }

    // get trial displacements and take difference
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    Vector  diffD  = disp2-disp1;
    if (diffD != 0.0)
      d0 = new Vector(diffD);
}   	 

int
TPB1D::commitState()
{
  return theMaterial->commitState();;
}

int
TPB1D::revertToLastCommit()
{
  return theMaterial->revertToLastCommit();
}


int
TPB1D::revertToStart()
{   
  return theMaterial->revertToStart();
}

int
TPB1D::update(void)
{
    // get trial displacements and take difference
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();

    double d = disp2(direction)-disp1(direction);

    if (d0 != 0)
      d -= (*d0)(direction);

    return theMaterial->setTrialStrain(d);
}

const Matrix &
TPB1D::getTangentStiff(void)
{
    double E;

    // stiff is a reference to the matrix holding the stiffness matrix
    Matrix& stiff = *theMatrix;
    
    // zero stiffness matrix
    stiff.Zero();
    
    E = theMaterial->getTangent();
    
    // compute contribution of material to tangent matrix
    stiff(direction, direction) = E;
    stiff(direction, direction+numDOF/2) = -E;
    stiff(direction+numDOF/2, direction) = -E;
    stiff(direction+numDOF/2, direction+numDOF/2) = E;

    return stiff;
}


const Matrix &
TPB1D::getInitialStiff(void)
{
    double E;

    // stiff is a reference to the matrix holding the stiffness matrix
    Matrix& stiff = *theMatrix;
    
    // zero stiffness matrix
    stiff.Zero();
    

    E = theMaterial->getInitialTangent();

    // compute contribution of material to tangent matrix
    stiff(direction, direction) = E;
    stiff(direction, direction+numDOF/2) = -E;
    stiff(direction+numDOF/2, direction) = -E;
    stiff(direction+numDOF/2, direction+numDOF/2) = E;      

    return stiff;
}
    

const Matrix &
TPB1D::getDamp(void)
{
  theMatrix->Zero();    
  return *theMatrix; 
}


const Matrix &
TPB1D::getMass(void)
{
    // no mass 
    theMatrix->Zero();    
    return *theMatrix; 
}


void 
TPB1D::zeroLoad(void)
{
  // does nothing now
}

int 
TPB1D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "TPB1D::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  
  return -1;
}

int 
TPB1D::addInertiaLoadToUnbalance(const Vector &accel)
{
  // does nothing as element has no mass yet!
  return 0;
}


const Vector &
TPB1D::getResistingForce()
{
  double force;
  
  // zero the residual
  theVector->Zero();
  
  // get resisting force for material
  force = theMaterial->getStress();
  
  (*theVector)(direction)  = -force;
  (*theVector)(direction+numDOF/2)  = force;

  return *theVector;
}


const Vector &
TPB1D::getResistingForceIncInertia()
{	
  return this->getResistingForce();
}


int
TPB1D::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << " TPB1D::sendSelf -- failed to send\n";
  return -1;
}

int
TPB1D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  opserr << " TPB1D::sendSelf -- failed to recv\n";
  return -1;
}

void
TPB1D::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "Element: " << this->getTag();
        s << " type: TPB1D  iNode: " << connectedExternalNodes(0);
        s << " jNode: " << connectedExternalNodes(1) << endln;
        s << " direction: " << direction << "\n";
        opserr << " mu1: " << mu[0] << endln;
        opserr << " mu2: " << mu[1] << endln;
        opserr << " mu3: " << mu[2] << endln;
        
        opserr << " R1: " << R[0] << endln;
        opserr << " R2: " << R[1] << endln;
        opserr << " R3: " << R[2] << endln;
        
        opserr << " h1: " << h[0] << endln;
        opserr << " h2: " << h[1] << endln;
        opserr << " h3: " << h[2] << endln;
        
        opserr << " D1: " << D[0] << endln;
        opserr << " D2: " << D[1] << endln;
        opserr << " D3: " << D[2] << endln;
        
        opserr << " d1: " << d[0] << endln;
        opserr << " d2: " << d[1] << endln;
        opserr << " d3: " << d[2] << endln;
        
        s << "\tMaterial: \n";
        s << *(theMaterial);
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"TPB1D\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"direction\": " << direction << ", ";
        s << "\"mu1\": " << mu[0] << ", ";
        s << "\"mu2\": " << mu[1] << ", ";
        s << "\"mu3\": " << mu[2] << ", ";
        s << "\"R1\": " << R[0] << ", ";
        s << "\"R2\": " << R[1] << ", ";
        s << "\"R3\": " << R[2] << ", ";
        s << "\"h1\": " << h[0] << ", ";
        s << "\"h2\": " << h[1] << ", ";
        s << "\"h3\": " << h[2] << ", ";
        s << "\"D1\": " << D[0] << ", ";
        s << "\"D2\": " << D[1] << ", ";
        s << "\"D3\": " << D[2] << ", ";
        s << "\"d1\": " << d[0] << ", ";
        s << "\"d2\": " << d[1] << ", ";
        s << "\"d3\": " << d[2] << ", ";
        s << "\"material\": \"" << theMaterial->getTag() << "\"}";
    }
}

Response*
TPB1D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","TPB1D");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    char outputData[10];

    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForces") == 0) || (strcmp(argv[0],"globalforces") == 0)) {

            char outputData[10];
            int numDOFperNode = numDOF/2;
            for (int i=0; i<numDOFperNode; i++) {
                sprintf(outputData,"P1_%d", i+1);
                output.tag("ResponseType", outputData);
            }
            for (int j=0; j<numDOFperNode; j++) {
                sprintf(outputData,"P2_%d", j+1);
                output.tag("ResponseType", outputData);
            }
            theResponse = new ElementResponse(this, 1, Vector(numDOF));
    }

    // a material quantity
    else if (strcmp(argv[0],"material") == 0) {
      theResponse =  theMaterial->setResponse(&argv[1], argc-1, output);
    }

    output.endTag();

    return theResponse;
}

int 
TPB1D::getResponse(int responseID, Information &eleInformation)
{
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    const Vector  diff  = disp2-disp1;

    switch (responseID) {
    case -1:
        return -1;

    case 1:
        return eleInformation.setVector(this->getResistingForce());

    case 2:
        if (eleInformation.theVector != 0) {
	  (*(eleInformation.theVector))(0) = theMaterial->getStress();
        }
        return 0;

    case 3:
        if (eleInformation.theVector != 0) {
	  (*(eleInformation.theVector))(0) = theMaterial->getStrain();
        }
        return 0;

    case 4:
        if (eleInformation.theVector != 0) {
	  (*(eleInformation.theVector))(0) = theMaterial->getStrain();
	  (*(eleInformation.theVector))(1) = theMaterial->getStress();
	}
	return 0;      
    
    default:
      return -1;
    }
}

int
TPB1D::setParameter(const char **argv, int argc, Parameter &param)
{
  int result = -1;  

  return result;
}

int
TPB1D::updateParameter (int parameterID, Information &info)
{
  return 0;
}

int
TPB1D::activateParameter(int passedParameterID)
{
  
  return 0;
}



int 
TPB1D::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
  return 0;
}
