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
************************************************************************
** 
** 
**************************************************************************************************************
** InertiaTruss Element                                                                                     **
** First published in: Ji X, Cheng Y, Molina Hutt C.                                                        **
** Seismic response of a tuned viscous mass damper (TVMD) coupled wall system. Eng Struct 2020;225:111252.  **
** https://doi.org/10.1016/j.engstruct.2020.111252.                                                         **
** *********************************************************************************************************** 
*/

// Description: This file contains the implementation for the InertiaTruss class.
// Not Ready for sensitivity analysis yet
// modified from truss.cpp    author: Frank McKenna

#include <InertiaTruss.h>
#include <Information.h>
#include <Parameter.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

//#include <fstream>

// initialise the class wide variables
Matrix InertiaTruss::trussM2(2,2);
Matrix InertiaTruss::trussM4(4,4);
Matrix InertiaTruss::trussM6(6,6);
Matrix InertiaTruss::trussM12(12,12);
Vector InertiaTruss::trussV2(2);
Vector InertiaTruss::trussV4(4);
Vector InertiaTruss::trussV6(6);
Vector InertiaTruss::trussV12(12);
#include <elementAPI.h>

#define OPS_Export

static int numMyTruss = 0;

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the inertiatruss end nodes.


OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_InertiaTrussElement)
{
  // print out a message about who wrote this element & any copyright info wanted
  if (numMyTruss == 0) {
    opserr << " \n";
	opserr << "                          InertiaTruss element v1.0\n";
	opserr << "                    by Xiaodong Ji, Yuhao Cheng, Yue Yu\n";
	opserr << "                           Tsinghua University\n";
	opserr << "Please contact jixd@mail.tsinghua.edu.cn, yuhao_cheng@126.com if anything goes wrong\n";
	opserr << " \n";
	numMyTruss++;

  }
	
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs != 4) {
    opserr << "Invalid Args want: element InertiaTruss $tag $iNode $jNode $mr\n";
    return 0;	
  }

  int iData[3];
  double mr = 0.0;
  int ndm = OPS_GetNDM();

  int numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode) in element InertiaTruss " << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &mr) != 0) {
    opserr << "WARNING: Invalid mr: element InertiaTruss " << iData[0] << 
      " $iNode $jNode $mr\n";
    return 0;	
  }

  // now create the InertiaTruss
  theElement = new InertiaTruss(iData[0], ndm, iData[1], iData[2], mr);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element InertiaTruss " << iData[0] << 
      " $iNode $jNode $mr\n";
  }

  return theElement;
}

InertiaTruss::InertiaTruss(int tag, int dim,
	int Nd1, int Nd2,
	double mr)
 :Element(tag,ELE_TAG_InertiaTruss),
  connectedExternalNodes(2),
  dimension(dim), numDOF(0),
  theLoad(0), theMatrix(0), theVector(0),
  L(0.0),  mass(mr), initialDisp(0)
{
  
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL InertiaTruss::InertiaTruss - " <<  tag << "failed to create an ID of size 2\n";
      exit(-1);
    }

    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;        

    // set node pointers to NULL
    for (int i=0; i<2; i++)
      theNodes[i] = 0;

    cosX[0] = 0.0;
    cosX[1] = 0.0;
    cosX[2] = 0.0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	theLoadSens = 0;
// AddingSensitivity:END //////////////////////////////////////
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
InertiaTruss::InertiaTruss()
:Element(0,ELE_TAG_InertiaTruss),     
 connectedExternalNodes(2),
 dimension(0), numDOF(0),
 theLoad(0), theMatrix(0), theVector(0),
 L(0.0), mass(0.0)
{
    // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL InertiaTruss::InertiaTruss - failed to create an ID of size 2\n";
      exit(-1);
  }
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  cosX[0] = 0.0;
  cosX[1] = 0.0;
  cosX[2] = 0.0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	theLoadSens = 0;
// AddingSensitivity:END //////////////////////////////////////
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
InertiaTruss::~InertiaTruss()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theLoad != 0)
	delete theLoad;
    if (theLoadSens != 0)
	delete theLoadSens;
    if (initialDisp != 0)
      delete [] initialDisp;
}


int
InertiaTruss::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
InertiaTruss::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
InertiaTruss::getNodePtrs(void) 
{
  return theNodes;
}

int
InertiaTruss::getNumDOF(void) 
{
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the truss element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
InertiaTruss::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	L = 0;
	return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	
    
    // if can't find both - send a warning message
    if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
      if (theNodes[0] == 0)
	opserr <<"InertiaTruss::setDomain() - truss" << this->getTag() << " node " << Nd1 <<
	  "does not exist in the model\n";
      else
	opserr <<"InertiaTruss::setDomain() - truss" << this->getTag() << " node " << Nd2 <<
	  "does not exist in the model\n";

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
      return;
    }

    // now determine the number of dof and the dimesnion    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2) {
      opserr <<"WARNING InertiaTruss::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
	"have differing dof at ends for truss " << this->getTag() << endln;

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
	
      return;
    }	

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    // now set the number of dof for element and set matrix and vector pointer
    if (dimension == 1 && dofNd1 == 1) {
	numDOF = 2;    
	theMatrix = &trussM2;
	theVector = &trussV2;
    }
    else if (dimension == 2 && dofNd1 == 2) {
	numDOF = 4;
	theMatrix = &trussM4;
	theVector = &trussV4;	
    }
    else if (dimension == 2 && dofNd1 == 3) {
	numDOF = 6;	
	theMatrix = &trussM6;
	theVector = &trussV6;		
    }
    else if (dimension == 3 && dofNd1 == 3) {
	numDOF = 6;	
	theMatrix = &trussM6;
	theVector = &trussV6;			
    }
    else if (dimension == 3 && dofNd1 == 6) {
	numDOF = 12;	    
	theMatrix = &trussM12;
	theVector = &trussV12;			
    }
    else {
      opserr <<"WARNING InertiaTruss::setDomain cannot handle " << dimension << " dofs at nodes in " << 
	dofNd1  << " problem\n";

      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
      return;
    }

    // create the load vector
    if (theLoad == 0)
      theLoad = new Vector(numDOF);
    else if (theLoad->Size() != numDOF) {
      delete theLoad;
      theLoad = new Vector(numDOF);
    }

    if (theLoad == 0) {
      opserr << "InertiaTruss::setDomain - truss " << this->getTag() << 
	"out of memory creating vector of size" << numDOF << endln;
      exit(-1);
      return;
    }          
    
    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();

    if (dimension == 1) {
      double dx = end2Crd(0)-end1Crd(0);

      if (initialDisp == 0) {
	double iDisp = end2Disp(0)-end1Disp(0);

	if (iDisp != 0) {
	  initialDisp = new double[1];
	  initialDisp[0] = iDisp;
	  dx += iDisp;
	}
      }
      L = sqrt(dx*dx);
      
      if (L == 0.0) {
	opserr <<"WARNING InertiaTruss::setDomain() - truss " << this->getTag() << " has zero length\n";
	return;
      }
      
      cosX[0] = 1.0;
    }
    else if (dimension == 2) {
      double dx = end2Crd(0)-end1Crd(0);
      double dy = end2Crd(1)-end1Crd(1);	
    
      if (initialDisp == 0) {
	double iDispX = end2Disp(0)-end1Disp(0);
	double iDispY = end2Disp(1)-end1Disp(1);
	if (iDispX != 0 || iDispY != 0) {
	  initialDisp = new double[2];
	  initialDisp[0] = iDispX;
	  initialDisp[1] = iDispY;
	  dx += iDispX;
	  dy += iDispY;
	}
      }
      
      L = sqrt(dx*dx + dy*dy);
      
      if (L == 0.0) {
	opserr <<"WARNING InertiaTruss::setDomain() - truss " << this->getTag() << " has zero length\n";
	return;
      }
	
      cosX[0] = dx/L;
      cosX[1] = dy/L;
    }
    else {

      double dx = end2Crd(0)-end1Crd(0);
      double dy = end2Crd(1)-end1Crd(1);	
      double dz = end2Crd(2)-end1Crd(2);		

      if (initialDisp == 0) {
	double iDispX = end2Disp(0)-end1Disp(0);
	double iDispY = end2Disp(1)-end1Disp(1);      
	double iDispZ = end2Disp(2)-end1Disp(2);      
	if (iDispX != 0 || iDispY != 0 || iDispZ != 0) {
	  initialDisp = new double[3];
	  initialDisp[0] = iDispX;
	  initialDisp[1] = iDispY;
	  initialDisp[2] = iDispZ;
	  dx += iDispX;
	  dy += iDispY;
	  dz += iDispZ;
	}
      }
      
      L = sqrt(dx*dx + dy*dy + dz*dz);
      
      if (L == 0.0) {
	opserr <<"WARNING InertiaTruss::setDomain() - inertiatruss " << this->getTag() << " has zero length\n";
	return;
      }
	
	cosX[0] = dx/L;
	cosX[1] = dy/L;
	cosX[2] = dz/L;
    }
}   	 


int
InertiaTruss::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "InertiaTruss::commitState () - failed in base class";
  }    
  return retVal;
}

int
InertiaTruss::revertToLastCommit()
{
	return 0;
}

int
InertiaTruss::revertToStart()
{
	return 0;
}

int
InertiaTruss::update(void)
{
	double strain = this->computeCurrentStrain();
	double rate = this->computeCurrentStrainRate();
	return 0;
}


const Matrix &
InertiaTruss::getTangentStiff(void)
{
    // set stiffness 0
    Matrix &stiff = *theMatrix;
    stiff.Zero();
    return stiff;
}


const Matrix &
InertiaTruss::getInitialStiff(void)
{
    // set stiffness 0
    Matrix &stiff = *theMatrix;
    stiff.Zero();
    return *theMatrix;
}



const Matrix &
InertiaTruss::getMass(void)
{
  // zero the matrix
  Matrix &Mmass = *theMatrix;
  Mmass.Zero();
  
  // check for quick return
  if (L == 0.0 || mass == 0.0) { // - problem in setDomain() no further warnings
    return Mmass;
  }
  

	  // inertia mass matrix
	  double m = mass;
	  int numDOF2 = numDOF / 2;
	  double temp;
	  for (int i = 0; i < dimension; i++) {
		  for (int j = 0; j < dimension; j++) {
			  temp = cosX[i] * cosX[j] * m;
			  Mmass(i, j) = temp;
			  Mmass(i + numDOF2, j) = -temp;
			  Mmass(i, j + numDOF2) = -temp;
			  Mmass(i + numDOF2, j + numDOF2) = temp;

	  }
  }
  
      return Mmass;
}
  
void 
InertiaTruss::zeroLoad(void)
{
  theLoad->Zero();
}

int 
InertiaTruss::addLoad(ElementalLoad *theLoad, double loadFactor)

{  
  opserr <<"InertiaTruss::addLoad - load type unknown for InertiaTruss with tag: " << this->getTag() << endln; 
  return -1;
}

int 
InertiaTruss::addInertiaLoadToUnbalance(const Vector &accel)
{
  // check for a quick return
  if (L == 0.0 || mass == 0.0) 
    return 0;
  
  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);    
  
  int nodalDOF = numDOF/2;
  
#ifdef _G3DEBUG    
  if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size()) {
    opserr <<"InertiaTruss::addInertiaLoadToUnbalance " <<
      "matrix and vector sizes are incompatable\n";
    return -1;
  }
#endif
  
  // want to add ( - fact * M R * accel ) to unbalance
  
    double m = mass;
	opserr << m;
	Matrix &addinertiamass = *theMatrix;
	double temp;
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			temp = cosX[i] * cosX[j] * m;
			addinertiamass(i, j) = temp;
			addinertiamass(i + nodalDOF, j) = -temp;
			addinertiamass(i, j + nodalDOF) = -temp;
			addinertiamass(i + nodalDOF, j + nodalDOF) = temp;
		}
	}
	for (int k = 0; k < dimension; k++) {
		for (int l = 0; l < dimension; l++) {
			(*theLoad)(k) -= addinertiamass(k, l)*Raccel1(l) + addinertiamass(k, l + nodalDOF)*Raccel2(l);
			(*theLoad)(k + nodalDOF) -= addinertiamass(k + nodalDOF, l)*Raccel1(l) + addinertiamass(k + nodalDOF, l + nodalDOF)*Raccel2(l);
		}
	}
  
  return 0;
}


int 
InertiaTruss::addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool somethingRandomInMotions)
{
	opserr << "InertiaTruss::addInertiaLoadSensitivityToUnbalance " <<
		"not ready for sensitivity analysis yet\n";
	return -1;
}

const Vector &
InertiaTruss::getResistingForce()
{	
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theVector->Zero();
	return *theVector;
    }
    
    // R = Ku - Pext
    // Ku = F * transformation
    int numDOF2 = numDOF/2;
    double temp;
    for (int i = 0; i < dimension; i++) {
      (*theVector)(i) = 0;
      (*theVector)(i+numDOF2) = 0;
    }
    
    return *theVector;
}


const Vector &
InertiaTruss::getResistingForceIncInertia()
{	
  this->getResistingForce();
  
  // subtract external load
  (*theVector) -= *theLoad;
  
  // now include the mass portion
  if (L != 0.0 && mass != 0.0) {
    
    // add inertia forces from element mass
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();	
    
    int numDOF2 = numDOF/2;
    
    
      // inertia mass matrix
	  double m = mass;
	  Matrix &getinertiaforce = *theMatrix;
	  double temp;
	  for (int i = 0; i < dimension; i++) {
		  for (int j = 0; j < dimension; j++) {
			  temp = cosX[i] * cosX[j] * m;
			  getinertiaforce(i, j) = temp;
			  getinertiaforce(i + numDOF2, j) = -temp;
			  getinertiaforce(i, j + numDOF2) = -temp;
			  getinertiaforce(i + numDOF2, j + numDOF2) = temp;
		  }
	  }
	  for (int k = 0; k < dimension; k++) {
		  for (int l = 0; l < dimension; l++) {
			  (*theVector)(k) += getinertiaforce(k, l)*accel1(l) + getinertiaforce(k, l + numDOF2)*accel2(l);
			  (*theVector)(k + numDOF2) += getinertiaforce(k + numDOF2, l)*accel1(l) + getinertiaforce(k + numDOF2, l + numDOF2)*accel2(l);
		  }
	  }  
  }
  
  return *theVector;
}

int
InertiaTruss::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  // truss packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(12);
  data(0) = this->getTag();
  data(1) = dimension;
  data(2) = numDOF;
  data(3) = mass; 
  if (initialDisp != 0) {
    for (int i=0; i<dimension; i++) {
      data[4+i] = initialDisp[i];
    }
  }

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING InertiaTruss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }	      

  // inertiatruss then sends the tags of it's two end nodes
  res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr <<"WARNING InertiaTruss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  return 0;
}

int
InertiaTruss::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();

  // inertiatruss creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(12);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING InertiaTruss::recvSelf() - failed to receive Vector\n";
    return -1;
  }	      

  this->setTag((int)data(0));
  dimension = (int)data(1);
  numDOF = (int)data(2);
  mass = data(3);


  initialDisp = new double[dimension];
  for (int i=0; i<dimension; i++)
    initialDisp[i] = 0.0;

  int initial = 0;
  for (int i=0; i<dimension; i++) {
    if (data(4+i) != 0.0) {
      initial = 1;
    }
  }
  
  if (initial != 0) {
    for (int i=0; i<dimension; i++) {
      initialDisp[i] = data(4+i);
    }    
  }
  
  // inertiatruss now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr <<"WARNING InertiaTruss::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }
  return 0;
}

int
InertiaTruss::displaySelf(Renderer &theViewer, int displayMode, float fact, 
		   const char **displayModes, int numModes)
{
    static Vector v1(3);
    static Vector v2(3);
    int res = 0;
    if (L == 0.0)
        return res;
    theNodes[0]->getDisplayCrds(v1, fact);
    theNodes[1]->getDisplayCrds(v2, fact);  
    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag(),0);
}



void
InertiaTruss::Print(OPS_Stream &s, int flag)
{
    //double force=0;
    
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "Element: " << this->getTag();
		s << " type: InertiaTruss  iNode: " << connectedExternalNodes(0);
		s << " jNode: " << connectedExternalNodes(1);
        s << " mr: " << mass;//mass of inertia
        
		if (initialDisp != 0) {
			s << " initialDisplacements: ";
			for (int i = 0; i < dimension; i++)
				s << initialDisp[i] << " ";
		}
		s << endln;
	}
    if (flag == 1) {
        s << "Nothing to be printed." << endln;
    }
    
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"InertiaTruss\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
		s << "\"mr\": " << mass << ", ";
	}
}

double
InertiaTruss::computeCurrentStrain(void) const
{
    // NOTE method will not be called if L == 0

    // determine the strain
    const Vector &disp1 = theNodes[0]->getTrialDisp();
    const Vector &disp2 = theNodes[1]->getTrialDisp();	

    double dLength = 0.0;
    if (initialDisp == 0)
      for (int i = 0; i < dimension; i++)
	dLength += (disp2(i)-disp1(i))*cosX[i];
    else
      for (int i = 0; i < dimension; i++)
	dLength += (disp2(i)-disp1(i)-initialDisp[i])*cosX[i];
  
    // this method should never be called with L == 0
    return dLength/L;
}

double
InertiaTruss::computeCurrentStrainRate(void) const
{
    // NOTE method will not be called if L == 0
    // determine the strain
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();	

    double dLength = 0.0;
    for (int i = 0; i < dimension; i++)
      dLength += (vel2(i)-vel1(i))*cosX[i];

    // this method should never be called with L == 0
    return dLength/L;
}

Response*
InertiaTruss::setResponse(const char **argv, int argc, OPS_Stream &output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","InertiaTruss");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

   // we compare argv[0] for known response types for the InertiaTruss
    if ((strcmp(argv[0], "relativeAcceleration") == 0) ||
        (strcmp(argv[0], "acceleration") == 0) ||
        (strcmp(argv[0], "accel") == 0) ||
        (strcmp(argv[0], "relAccel") == 0)) {
        output.tag("ResponseType", "acceleration");
        theResponse = new ElementResponse(this, 1, Vector(1));
    }
    else if ((strcmp(argv[0], "axialForce") == 0) ||
        (strcmp(argv[0], "basicForce") == 0) ||
        (strcmp(argv[0], "basicForces") == 0)) {
        output.tag("ResponseType", "N");
        theResponse = new ElementResponse(this, 2, Vector(1));
    }
    output.endTag();
    return theResponse;
}

int 
InertiaTruss::getResponse(int responseID, Information &eleInfo)
{
    int numDOF2 = numDOF / 2;
    static Vector fVec_Acc(1);
    static Vector fVec_AxialForce(1);
    const Vector& acc1 = theNodes[0]->getTrialAccel();
    const Vector& acc2 = theNodes[1]->getTrialAccel();
    Vector toret = acc2 - acc1;

    switch (responseID) {
    case 1://accel
        fVec_Acc(0) = 0;
        for (int i = 0; i < dimension; i++)
        {
            fVec_Acc(0) += toret(i) * cosX[i];
        }
        return eleInfo.setVector(fVec_Acc);

    case 2://axial force
        fVec_AxialForce(0) = 0;
        for (int i = 0; i < dimension; i++)
        {
            fVec_AxialForce(0) += toret(i) * mass * cosX[i];
        }
        return eleInfo.setVector(fVec_AxialForce);
    default:
        return 0;
    }
}

// AddingSensitivity:BEGIN ///////////////////////////////////
// InertiaTruss::addInertiaLoadSensitivityToUnbalance not ready for sensitivity analysis yet
int
InertiaTruss::setParameter(const char **argv, int argc, Parameter &param)
{ 
	opserr << "InertiaTruss::addInertiaLoadSensitivityToUnbalance " <<
		"not ready for sensitivity analysis yet\n";
	return -1;
}

int
InertiaTruss::updateParameter (int parameterID, Information &info)
{
	opserr << "InertiaTruss::addInertiaLoadSensitivityToUnbalance " <<
		"not ready for sensitivity analysis yet\n";
	return -1;
}

int
InertiaTruss::activateParameter(int passedParameterID)
{
	opserr << "InertiaTruss::addInertiaLoadSensitivityToUnbalance " <<
		"not ready for sensitivity analysis yet\n";
	return -1;
}


const Matrix &
InertiaTruss::getKiSensitivity(int gradNumber)
{
	opserr << "InertiaTruss::addInertiaLoadSensitivityToUnbalance " <<
		"not ready for sensitivity analysis yet\n";
	Matrix &Error = *theMatrix;
	Error.Zero();
	return Error;
}

const Matrix &
InertiaTruss::getMassSensitivity(int gradNumber)
{
	opserr << "InertiaTruss::addInertiaLoadSensitivityToUnbalance " <<
		"not ready for sensitivity analysis yet\n";
	Matrix &Mmass = *theMatrix;
	Mmass.Zero();

	return Mmass;
}

const Vector &
InertiaTruss::getResistingForceSensitivity(int gradNumber)
{
	opserr << "InertiaTruss::addInertiaLoadSensitivityToUnbalance " <<
		"not ready for sensitivity analysis yet\n";
	theVector->Zero();

	return *theVector;
}

int
InertiaTruss::commitSensitivity(int gradNumber, int numGrads)
{
	opserr << "InertiaTruss::addInertiaLoadSensitivityToUnbalance " <<
		"not ready for sensitivity analysis yet\n";
	return -1;
}

// AddingSensitivity:END /////////////////////////////////////////////
