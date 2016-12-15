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
** file 'COPYRIGHT'  in main dInrectory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Written: fmk
//
// What: "@(#) MultiFP2d.C, revA"

// we specify what header files we need
#include "MultiFP2d.h"
#include <elementAPI.h>
#include <G3Globals.h>
#include <math.h>
#include <float.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <ENTMaterial.h>
#include <MultiLinear.h>


#include <math.h>
#include <stdlib.h>
#include <string.h>

static Matrix ZeroLengthM4(4,4);   // class wide matrix for 4*4
static Matrix ZeroLengthM6(6,6);   // class wide matrix for 6*6

#include <elementAPI.h>

void *
OPS_MultiFP2d()
{
  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 3) {
    opserr << "WARNING::MultiFP2d insufficient args\n";
    return theEle;
  }

  // get the id and end nodes 
  int iData[5];
  int numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING::MultiFP2d invalid element data\n";
    return 0;
  }
  int eleTag = iData[0];

  numRemainingArgs -= 3;

  const char *nextArg;
  int done = 0;

  int axialCase = 1; //

  opserr << "NUM REMAINING ARGS: " << numRemainingArgs << endln;

  while (done == 0 && numRemainingArgs > 0) {

    nextArg = OPS_GetString();
    //    OPS_GetString(nextArg, 30);
    numRemainingArgs--;

    if (strcmp(nextArg, "-W0") == 0) {
      axialCase = 0;
    } else if (strcmp(nextArg, "-WC") == 0) {
      axialCase = 1;
    } else if (strcmp(nextArg, "-WT") == 0) {
      axialCase = 2;
    } 

    else if (strcmp(nextArg, "-material") == 0) {
      if (numRemainingArgs == 3) { //user defined material
	numData = 2;
	if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
	  opserr << "WARNING invalid element data\n";
	  return 0;
	}
	
	double dData[1];
	numData = 1;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
	  opserr << "WARNING error reading element area for element" << eleTag << endln;
	  return 0;
	}
	
	UniaxialMaterial *theFrictionMaterial;
	UniaxialMaterial *theVerticalMaterial;
	theFrictionMaterial = OPS_GetUniaxialMaterial(iData[3]);
	theVerticalMaterial = OPS_GetUniaxialMaterial(iData[4]);
	
	theEle = new MultiFP2d(eleTag, iData[1],iData[2], 
			       theFrictionMaterial, 
			       theVerticalMaterial,
			       dData[0], axialCase);
	done = 1;
      } else {
	opserr << "WARNING incorrect #args for MultiFP ele " << eleTag << " for -material option" << endln;

      }
    } 

    else if (strcmp(nextArg, "-triple") == 0) {      
      if (numRemainingArgs == 17) { // triple friction pendulum
	//    R1 R2 R3 h1 h2 h3 D1 D2 D3 d1 d2 d3 mu1 mu2 mu3 Kvert W0 
	int type = 3;
	double dData[17];
	numData = 17;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
	  opserr << "WARNING error reading element area for element" << eleTag << endln;
	  return 0;
	}
	Vector R(3); R(0) = dData[0]; R(1)=dData[1]; R(2)=dData[2];
	Vector h(3); h(0) = dData[3]; h(1)=dData[4]; h(2)=dData[5];
	Vector D(3); D(0) = dData[6]; D(1)=dData[7]; D(2)=dData[8];
	Vector d(3); d(0) = dData[9]; d(1)=dData[10]; d(2)=dData[11];
	Vector mu(3); mu(0) = dData[12]; mu(1)=dData[13]; mu(2)=dData[14];
	
	theEle = new MultiFP2d(eleTag, iData[1], iData[2], 
			       type,
			       R,
			       h,
			       D,
			       d,
			       mu,
			       dData[15],
			       dData[16],
			       axialCase);    
	
	done = 1;
      } else {
	opserr << "WARNING incorrect #args for MultiFP ele " << eleTag << " for -triple option" << endln;
      }
    } else {
      opserr << "WARNING unknown option: " << nextArg << " for MultiFP ele " << eleTag << endln;
      done = 1;
    }
    
    if (theEle == 0) {
      opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
      return 0;
    }
  }

  return theEle;
}



// typical constructor
MultiFP2d::MultiFP2d(int tag, 
		     int Nd1, int Nd2, 
		     UniaxialMaterial *friction,
		     UniaxialMaterial *vertical,
		     double w0, int aCase)
  :Element(tag, ELE_TAG_MultiFP2d),     
   externalNodes(2),
   numDOF(0), theMatrix(0), theVector(0),
   type(0), axialCase(aCase)
{	
  theFrictionModel = friction->getCopy();
  theVerticalModel = vertical->getCopy();

  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;        
  
  theNodes[0] = 0; 
  theNodes[1] = 0;

  W0 = w0;
  cW = w0;
}

// typical constructor
MultiFP2d::MultiFP2d(int tag, 
		     int Nd1, int Nd2,
		     int type,
		     const Vector &R,
		     const Vector &h,
		     const Vector &D,
		     const Vector &d,
		     const Vector &mu,
		     double Kv,
		     double w0,
		     int aCase)
  :Element(tag, ELE_TAG_MultiFP2d),     
   externalNodes(2),
   numDOF(0), theMatrix(0), theVector(0),
   type(0), axialCase(aCase)
{	
  theVerticalModel = new ENTMaterial(2, Kv);

  if (type == 3) {
    double L0 = R(0)-h(0);
    double L1 = R(1)-h(1);
    double L2 = R(2)-h(2);
    
    double u1bar = L1*(D(1)-d(1))/(2.0*R(1));
    double u2bar = L2*(D(2)-d(2))/(2.0*R(2));

    Vector fy(5);
    Vector u(5);
    
    fy(0) = mu(0);
    fy(1) = mu(1);
    fy(2) = mu(2);
    
    u(1) = 2*L0*(mu(1)-mu(0));
    u(0) = u(1)/100.0;
    u(2) = L0*(mu(1)+mu(2) - 2*mu(0)) + L1*(mu(2) - mu(1));
    u(3) = u(2) + (u1bar/L1 + mu(1) - mu(2))*(L1+L2);
    u(4) = u(3)+(u2bar/L2 + mu(2) - u1bar/L1 - mu(1))*(L0+L2);
    
    fy(3) = fy(2) + u1bar/L1 + mu(1) - mu(2);
    fy(4) = fy(3) + u2bar/L2 + mu(2) - u1bar/L1 - mu(1);
    
    theFrictionModel = new MultiLinear(1, fy, u); 
  }

  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;        
  
  theNodes[0] = 0; 
  theNodes[1] = 0;

  W0 = w0;
  cW = w0;
}

// constructor which should be invoked by an FE_ObjectBroker only
MultiFP2d::MultiFP2d()
 :Element(0, ELE_TAG_MultiFP2d), 
  externalNodes(2),
  numDOF(0), theMatrix(0), theVector(0)
{
  theNodes[0] = 0; 
  theNodes[1] = 0;
}

//  destructor - provided to clean up any memory
MultiFP2d::~MultiFP2d()
{
  if (theMatrix != 0)
    delete theMatrix;
  if (theVector != 0)
    delete theVector;
  
  if (theFrictionModel != 0)
    delete theFrictionModel;

  if (theVerticalModel != 0)
    delete theVerticalModel;
}

int
MultiFP2d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
MultiFP2d::getExternalNodes(void) 
{
  return externalNodes;
}

Node **
MultiFP2d::getNodePtrs(void) 
{
  return theNodes;
}

int
MultiFP2d::getNumDOF(void) {
  return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void
MultiFP2d::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    exit(-1);    
    return;
  }
  
  // first ensure nodes exist in Domain and set the node pointers
  Node *end1Ptr, *end2Ptr;
  int Nd1 = externalNodes(0);
  int Nd2 = externalNodes(1);
  end1Ptr = theDomain->getNode(Nd1);
  end2Ptr = theDomain->getNode(Nd2);	
  if (end1Ptr == 0) {
    opserr << "WARNING MultiFP2d::setDomain() - at truss " << this->getTag() << " node " <<
      Nd1 << "  does not exist in domain\n";
    exit(-1);    
    return;  // don't go any further - otherwise segemntation fault
  }
  if (end2Ptr == 0) {        
    opserr << "WARNING MultiFP2d::setDomain() - at truss " << this->getTag() << " node " <<
      Nd2 << "  does not exist in domain\n";
    exit(-1);    
    return;  // don't go any further - otherwise segemntation fault
  }	
  theNodes[0] = end1Ptr;
  theNodes[1] = end2Ptr;
  // call the DomainComponent class method THIS IS VERY IMPORTANT
  this->DomainComponent::setDomain(theDomain);
  
  // ensure connected nodes have correct number of dof's
  int dofNd1 = end1Ptr->getNumberDOF();
  int dofNd2 = end2Ptr->getNumberDOF();	
  if ((dofNd1 != dofNd2) || ((dofNd2 != 2) && (dofNd2 != 3)) ) {
    opserr << "MultiFP2d::setDomain(): 2 or 3 dof required at nodes\n";
    exit(-1);
    return;
  }	

  if (dofNd2 == 2) {
    theMatrix = new Matrix(4,4);
    theVector = new Vector(4);
    numDOF = 4;
  } else {
    theMatrix = new Matrix(6,6);
    theVector = new Vector(6);
    numDOF = 6;
  }

  this->update();
}   	 


int
MultiFP2d::commitState()
{
  theFrictionModel->commitState();
  theVerticalModel->commitState();

  cW = -theVerticalModel->getStress();  

  return 0;
}

int
MultiFP2d::revertToLastCommit()
{
  theFrictionModel->revertToLastCommit();
  theVerticalModel->revertToLastCommit();
  
  return 0;
}

int
MultiFP2d::revertToStart()
{
  cW = W0;  

  return 0;
}

int
MultiFP2d::update()
{
  static Vector delU(4);
  static Vector delP(4);

  const Vector &v1 = theNodes[0]->getTrialDisp();
  const Vector &v2 = theNodes[1]->getTrialDisp();

  double strainX = v2(0)-v1(0);
  double strainY = v2(1)-v1(1);

  theFrictionModel->setTrialStrain(strainX);
  theVerticalModel->setTrialStrain(strainY);

  int numD = numDOF/2;

  double k1 = theFrictionModel->getTangent();
  double k2 = theVerticalModel->getTangent();
  
  double f1 = theFrictionModel->getStress();
  double f2 = theVerticalModel->getStress();

  double vLoad = cW;
  if (axialCase == 0)
    vLoad = W0;
  else if (axialCase == 2)
    vLoad = f2;

  k1 *= vLoad;
  f1 *= vLoad;

  theVector->Zero();
  
  (*theVector)(0) = -f1;
  (*theVector)(1) = -f2;
  (*theVector)(0+numD) = f1;
  (*theVector)(1+numD) = f2;

  theMatrix->Zero();
  
  (*theMatrix)(0,0)           = k1;
  (*theMatrix)(0+numD,0+numD) = k1;
  (*theMatrix)(0+numD,0)      = -k1;
  (*theMatrix)(0,0+numD)      = -k1;

  (*theMatrix)(1,1)           = k2;
  (*theMatrix)(1+numD,1+numD) = k2;
  (*theMatrix)(1+numD,1)      = -k2;
  (*theMatrix)(1,1+numD)      = -k2;

  return 0;
}


const Matrix &
MultiFP2d::getTangentStiff(void)
{
  return *theMatrix;
}

const Matrix &
MultiFP2d::getInitialStiff(void)
{
  int numD = numDOF/2;

  double k1 = theFrictionModel->getInitialTangent();
  double k2 = theVerticalModel->getInitialTangent();
  
  k1 *= W0;

  theMatrix->Zero();
  
  (*theMatrix)(0,0)           = k1;
  (*theMatrix)(0+numD,0+numD) = k1;
  (*theMatrix)(0+numD,0)      = -k1;
  (*theMatrix)(0,0+numD)      = -k1;

  (*theMatrix)(1,1)           = k2;
  (*theMatrix)(1+numD,1+numD) = k2;
  (*theMatrix)(1+numD,1)      = -k2;
  (*theMatrix)(1,1+numD)      = -k2;

  return *theMatrix;
}

const Vector &
MultiFP2d::getResistingForce()
{	
  return *theVector;
}

int
MultiFP2d::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
MultiFP2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

void
MultiFP2d::Print(OPS_Stream &s, int flag)
{
  s << "Element: " << this->getTag(); 
  s << " type: MultiFP2d  iNode: " << externalNodes(0);
  s << " jNode: " << externalNodes(1) << endln;
  s << "material for normalized lateral force displacement response\n";
  theFrictionModel->Print(s, flag);
  s << "material for vertical force displacement response\n";
  theVerticalModel->Print(s, flag);
}


Response *
MultiFP2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType",this->getClassType());
  output.attr("eleTag",this->getTag());
  int numNodes = this->getNumExternalNodes();
  const ID &nodes = this->getExternalNodes();
  static char nodeData[32];

  for (int i=0; i<numNodes; i++) {
    sprintf(nodeData,"node%d",i+1);
    output.attr(nodeData,nodes(i));
  }

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {
    const Vector &force = this->getResistingForce();
    int size = force.Size();
    for (int i=0; i<size; i++) {
      sprintf(nodeData,"P%d",i+1);
      output.tag("ResponseType",nodeData);
    }
    theResponse = new ElementResponse(this, 1, this->getResistingForce());
  } else if (strcmp(argv[0],"friction") == 0 || strcmp(argv[0],"frictionModel") == 0) {
      theResponse =  theFrictionModel->setResponse(&argv[1], argc-1, output);
  } else if (strcmp(argv[0],"vertical") == 0 || strcmp(argv[0],"verticalModel") == 0) {
      theResponse =  theVerticalModel->setResponse(&argv[1], argc-1, output);
  } 

  output.endTag();
  return theResponse;
}



int 
MultiFP2d::getResponse(int responseID, Information &eleInfo)
{
 // Vector res(this->getResistingForce());
 // res(2) = Ac;

  switch (responseID) {
  case -1:
    return -1;
  case 1: // global forces
    return eleInfo.setVector(this->getResistingForce());
    //  return eleInfo.setVector(res);
    
  case 2:
    return eleInfo.setVector(this->getRayleighDampingForces());
  default:
    return 0;
  }
}

int 
MultiFP2d::displaySelf(Renderer &theViewer,
    int displayMode, float fact, const char **modes, int numMode)
{
  return 0;
}
