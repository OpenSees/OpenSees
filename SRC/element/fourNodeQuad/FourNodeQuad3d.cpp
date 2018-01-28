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
                                                                        
// $Revision: 1.35 $
// $Date: 2009-10-13 21:14:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/FourNodeQuad3d.cpp,v $

// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
// Description: This file contains the class definition for FourNodeQuad3d.

#include <FourNodeQuad3d.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

#include <elementAPI.h>

double FourNodeQuad3d::matrixData[144];
Matrix FourNodeQuad3d::K(matrixData, 12, 12);
Vector FourNodeQuad3d::P(12);
double FourNodeQuad3d::shp[3][4];
double FourNodeQuad3d::pts[4][2];
double FourNodeQuad3d::wts[4];

void *
OPS_FourNodeQuad3d()
{

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new FourNodeQuad3d();
    return theEle;
  }

  if (numRemainingArgs != 8 && numRemainingArgs != 12) {
    opserr << "ERROR - FourNodeQuad3d not enough args provided, want: element FourNodeQuad3d tag? iNode? jNode? kNode? lNode? thickness? type? matID? <p? rho? b1? b2?>\n";
  }

  // get the id and end nodes 
  int iData[6];
  double dData[5];
  dData[1] = 0.0;
  dData[2] = 0.0;
  dData[3] = 0.0;
  dData[4] = 0.0;

  int numData;
  int matTag = 0;
  int eleTag = 0;
  const char *pType;

  numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING element FourNodeQuad3d : invalid element data\n";
    return 0;
  }
  eleTag = iData[0];

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING element FourNodeQuad3d : invalid thickness for element: " << eleTag << "\n";
    return 0;
  }

  pType = OPS_GetString();
  if (pType != 0) {
    opserr << "WARNING element FourNodeQuad3d : invalid pType for element: " << eleTag << "\n";
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &matTag) != 0) {
    opserr << "WARNING element FourNodeQuad3d : invalid matTag for element: " << eleTag << "\n";
    delete [] pType;
    return 0;
  }


  NDMaterial *theMaterial = OPS_getNDMaterial(matTag);
  
  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matTag << "not found for element " << eleTag << endln;
    return 0;
  }

  if (numRemainingArgs == 12) {
    numData = 4;
    if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
      opserr << "WARNING element FourNodeQuad3d : invalid optional args for element: " << eleTag << "\n";
      delete [] pType;
      return 0;
    }
  }

  // now create the truss and add it to the Domain
  theEle = new FourNodeQuad3d(eleTag, iData[1], iData[2], iData[3], iData[4],
			      *theMaterial, pType,
			      dData[0], dData[1], dData[2], dData[3], dData[4]);

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    delete theMaterial;
      delete [] pType;
    return 0;
  }

  delete [] pType;
  return theEle;
}


FourNodeQuad3d::FourNodeQuad3d(int tag, int nd1, int nd2, int nd3, int nd4,
			       NDMaterial &m, const char *type, double t,
			       double p, double r, double b1, double b2)
:Element (tag, ELE_TAG_FourNodeQuad3d), 
  theMaterial(0), connectedExternalNodes(4), 
 Q(12), pressureLoad(12), thickness(t), applyLoad(0), pressure(p), rho(r)
{
  pts[0][0] = -0.5773502691896258;
  pts[0][1] = -0.5773502691896258;
  pts[1][0] =  0.5773502691896258;
  pts[1][1] = -0.5773502691896258;
  pts[2][0] =  0.5773502691896258;
  pts[2][1] =  0.5773502691896258;
  pts[3][0] = -0.5773502691896258;
  pts[3][1] =  0.5773502691896258;
  
  wts[0] = 1.0;
  wts[1] = 1.0;
  wts[2] = 1.0;
  wts[3] = 1.0;
  
  if (strcmp(type,"PlaneStrain") != 0 && strcmp(type,"PlaneStress") != 0
      && strcmp(type,"PlaneStrain2D") != 0 && strcmp(type,"PlaneStress2D") != 0) {
    opserr << "FourNodeQuad3d::FourNodeQuad3d -- improper material type: " << type << "for FourNodeQuad3d\n";
    exit(-1);
  }
  
  // Body forces
  b[0] = b1;
  b[1] = b2;
  
  // Allocate arrays of pointers to NDMaterials
  theMaterial = new NDMaterial *[4];
  
  if (theMaterial == 0) {
    opserr << "FourNodeQuad3d::FourNodeQuad3d - failed allocate material model pointer\n";
    exit(-1);
  }
  
  int i;
  for (i = 0; i < 4; i++) {

    // Get copies of the material model for each integration point
    theMaterial[i] = m.getCopy(type);
    
    // Check allocation
    if (theMaterial[i] == 0) {
      opserr << "FourNodeQuad3d::FourNodeQuad3d -- failed to get a copy of material model\n";
      exit(-1);
    } 
  }
  
  // Set connected external node IDs
  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;
  connectedExternalNodes(2) = nd3;
  connectedExternalNodes(3) = nd4;
  
  for (i=0; i<4; i++) {
    theNodes[i] = 0;
  }
}

FourNodeQuad3d::FourNodeQuad3d()
 :Element (0,ELE_TAG_FourNodeQuad3d),
  theMaterial(0), connectedExternalNodes(4), 
  Q(12), pressureLoad(12), thickness(0.0), applyLoad(0), pressure(0.0), rho(0.0)
{
  pts[0][0] = -0.577350269189626;
  pts[0][1] = -0.577350269189626;
  pts[1][0] =  0.577350269189626;
  pts[1][1] = -0.577350269189626;
  pts[2][0] =  0.577350269189626;
  pts[2][1] =  0.577350269189626;
  pts[3][0] = -0.577350269189626;
  pts[3][1] =  0.577350269189626;
  
  wts[0] = 1.0;
  wts[1] = 1.0;
  wts[2] = 1.0;
  wts[3] = 1.0;
  
  for (int i=0; i<4; i++)
    theNodes[i] = 0;
}

FourNodeQuad3d::~FourNodeQuad3d()
{    
  for (int i = 0; i < 4; i++) {
    if (theMaterial[i])
      delete theMaterial[i];
  }
  
  // Delete the array of pointers to NDMaterial pointer arrays
  if (theMaterial)
    delete [] theMaterial;
}

int
FourNodeQuad3d::getNumExternalNodes() const
{
  return 4;
}

const ID&
FourNodeQuad3d::getExternalNodes()
{
  return connectedExternalNodes;
}


Node **
FourNodeQuad3d::getNodePtrs(void) 
{
  return theNodes;
}

int
FourNodeQuad3d::getNumDOF()
{
  return 12;
}

void
FourNodeQuad3d::setDomain(Domain *theDomain)
{
  // Check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    theNodes[2] = 0;
    theNodes[3] = 0;
    return;
  }
  
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);
  int Nd3 = connectedExternalNodes(2);
  int Nd4 = connectedExternalNodes(3);
  
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);
  theNodes[2] = theDomain->getNode(Nd3);
  theNodes[3] = theDomain->getNode(Nd4);
  
  if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
    opserr << "FATAL ERROR FourNodeQuad3d (tag: " << this->getTag() << " ) a node does not exist\n";
    exit(-1);
  }
  
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();
  int dofNd3 = theNodes[2]->getNumberDOF();
  int dofNd4 = theNodes[3]->getNumberDOF();
  
  if (dofNd1 != 3 || dofNd2 != 3 || dofNd3 != 3 || dofNd4 != 3) {
    opserr << "FATAL ERROR FourNodeQuad3d (tag: " << this->getTag() << " ) needs ndf = 3\n";
    exit(-1);
  }
  this->DomainComponent::setDomain(theDomain);
  
  // Compute consistent nodal loads due to pressure
  this->setPressureLoadAtNodes();
  
  const Vector &crds1 = theNodes[0]->getCrds();
  const Vector &crds2 = theNodes[1]->getCrds();
  const Vector &crds3 = theNodes[2]->getCrds();
  const Vector &crds4 = theNodes[3]->getCrds();
  
  if (crds1.Size() != 3 || crds2.Size() != 3 || crds3.Size() != 3 || crds4.Size() != 3) {
    opserr << "FATAL ERROR FourNodeQuad3d (tag: " << this->getTag() << " ) needs ndm = 3\n";
    exit(-1);
  }
  
  int dirns[3];
  for (int i=0; i<3; i++)
    dirns[i] = 1;
  
  if ((crds1(0) == crds2(0)) && (crds2(0) == crds3(0)) && (crds3(0) == crds4(0)))
    dirns[0] = 0;
  if ((crds1(1) == crds2(1)) && (crds2(1) == crds3(1)) && (crds3(1) == crds4(1)))
    dirns[1] = 0;
  if ((crds1(2) == crds2(2)) && (crds2(2) == crds3(2)) && (crds3(2) == crds4(2)))
    dirns[2] = 0;
  
  int sum = 0;
  for (int i=0; i<3; i++) {
    if (dirns[i] != 0 && sum < 2)
      dirn[sum] = i;
    sum += dirns[i];
  }
  
  if (sum != 2) {
    opserr << "DIRNS: " << dirns[0] << " " << dirns[1] << " " << dirns[2];
    theNodes[0]->Print(opserr);
    theNodes[1]->Print(opserr);
    theNodes[2]->Print(opserr);
    theNodes[3]->Print(opserr);
    opserr << "FATAL ERROR FourNodeQuad3d (tag: " << this->getTag() << " ) needs four nodes to be in x-y, y-z, or x-z plane\n";
    exit(-1);
  }

  return;
}

int
FourNodeQuad3d::commitState()
{
  int retVal = 0;
  
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "FourNodeQuad3d::commitState () - failed in base class";
  }    
  
  // Loop over the integration points and commit the material states
  for (int i = 0; i < 4; i++)
    retVal += theMaterial[i]->commitState();
  
  return retVal;
}

int
FourNodeQuad3d::revertToLastCommit()
{
  int retVal = 0;
  
  // Loop over the integration points and revert to last committed state
  for (int i = 0; i < 4; i++)
    retVal += theMaterial[i]->revertToLastCommit();
  
  return retVal;
}

int
FourNodeQuad3d::revertToStart()
{
  int retVal = 0;
  
  // Loop over the integration points and revert states to start
  for (int i = 0; i < 4; i++)
    retVal += theMaterial[i]->revertToStart();
  
  return retVal;
}


int
FourNodeQuad3d::update()
{
  const Vector &disp1 = theNodes[0]->getTrialDisp();
  const Vector &disp2 = theNodes[1]->getTrialDisp();
  const Vector &disp3 = theNodes[2]->getTrialDisp();
  const Vector &disp4 = theNodes[3]->getTrialDisp();
  
  static double u[2][4];
  
  u[0][0] = disp1(dirn[0]);
  u[1][0] = disp1(dirn[1]);
  u[0][1] = disp2(dirn[0]);
  u[1][1] = disp2(dirn[1]);
  u[0][2] = disp3(dirn[0]);
  u[1][2] = disp3(dirn[1]);
  u[0][3] = disp4(dirn[0]);
  u[1][3] = disp4(dirn[1]);
  
  static Vector eps(3);
  
  int ret = 0;

  // Loop over the integration points
  for (int i = 0; i < 4; i++) {
    
    // Determine Jacobian for this integration point
    this->shapeFunction(pts[i][0], pts[i][1]);

    // Interpolate strains
    //eps = B*u;
    //eps.addMatrixVector(0.0, B, u, 1.0);
    eps.Zero();

    for (int beta = 0; beta < 4; beta++) {
      eps(0) += shp[0][beta]*u[0][beta];
      eps(1) += shp[1][beta]*u[1][beta];
      eps(2) += shp[0][beta]*u[1][beta] + shp[1][beta]*u[0][beta];
    }

    // Set the material strain
    ret += theMaterial[i]->setTrialStrain(eps);
  }
  
  return ret;
}


const Matrix&
FourNodeQuad3d::getTangentStiff()
{
  K.Zero();
  
  double dvol;
  double DB[3][2];
  
  // Loop over the integration points
  for (int i = 0; i < 4; i++) {
    
    // Determine Jacobian for this integration point
    dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness*wts[i]);
    
    // Get the material tangent
    const Matrix &D = theMaterial[i]->getTangent();
    
    // Perform numerical integration
    //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
    //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
    
    double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
    double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
    double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);
    
    //	  for (int beta = 0, ib = 0, colIb =0, colIbP1 = 8; 
    //   beta < 4; 
    //   beta++, ib += 2, colIb += 16, colIbP1 += 16) {
    
    int diff = dirn[1]-dirn[0];
    
    for (int alpha = 0, ia = dirn[0]; alpha < 4; alpha++, ia += 3) {
      for (int beta = 0, ib = dirn[0]; beta < 4; beta++, ib += 3) {
	
	DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
	DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
	DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
	DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
	DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
	DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);
	
	
	K(ia,ib) += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	K(ia,ib+diff) += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	K(ia+diff,ib) += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	K(ia+diff,ib+diff) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	
	//matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	//matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	//matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	//matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	
      }
    }
  }
  
  return K;
}


const Matrix&
FourNodeQuad3d::getInitialStiff()
{
  K.Zero();
  
  double dvol;
  double DB[3][2];
  
  // Loop over the integration points
  for (int i = 0; i < 4; i++) {
    
    // Determine Jacobian for this integration point
    dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness*wts[i]);
    
    // Get the material tangent
    const Matrix &D = theMaterial[i]->getInitialTangent();

    double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
    double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
    double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);
    
    // Perform numerical integration
    //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
    //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);

    /* **********************************************************************
    for (int beta = 0, ib = 0, colIb =0, colIbP1 = 8; 
	 beta < 4; 
	 beta++, ib += 2, colIb += 16, colIbP1 += 16) {
      
      for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
	
	DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
	DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
	DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
	DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
	DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
	DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);

	K(ia,ib) += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	K(ia,ib+diff) += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	K(ia+diff,ib) += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	K(ia+diff,ib+diff) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	
	matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
      }
    }
    }

    Ki = new Matrix(K);
    return K;
    ************************************************************************** */

    int diff = dirn[1]-dirn[0];
    
    for (int alpha = 0, ia = dirn[0]; alpha < 4; alpha++, ia += 3) {
      for (int beta = 0, ib = dirn[0]; beta < 4; beta++, ib += 3) {
	
	DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
	DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
	DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
	DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
	DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
	DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);
	
	
	K(ia,ib) += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	K(ia,ib+diff) += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	K(ia+diff,ib) += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	K(ia+diff,ib+diff) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
      }
    }
  }

  return K;
}

const Matrix&
FourNodeQuad3d::getMass()
{
  K.Zero();
  
  int i;
  static double rhoi[4];
  double sum = 0.0;
  for (i = 0; i < 4; i++) {
    if (rho == 0)
      rhoi[i] = theMaterial[i]->getRho();
    else
      rhoi[i] = rho;	    
    sum += rhoi[i];
  }
  
  if (sum == 0.0)
    return K;
  
  double rhodvol, Nrho;
  
  // Compute a lumped mass matrix
  for (i = 0; i < 4; i++) {
    
    // Determine Jacobian for this integration point
    rhodvol = this->shapeFunction(pts[i][0], pts[i][1]);
    
    // Element plus material density ... MAY WANT TO REMOVE ELEMENT DENSITY
    rhodvol *= (rhoi[i]*thickness*wts[i]);

    int diff = dirn[1]-dirn[0];
    for (int alpha = 0, ia = dirn[0]; alpha < 4; alpha++, ia+=3) {
      Nrho = shp[2][alpha]*rhodvol;
      K(ia,ia) += Nrho;
      K(ia+diff,ia+diff) += Nrho;
    }
  }
  
  return K;
}

void
FourNodeQuad3d::zeroLoad(void)
{
  Q.Zero();
  
  applyLoad = 0;
  
  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  
  return;
}

int 
FourNodeQuad3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // Added option for applying body forces in load pattern: C.McGann, U.Washington
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  
  if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor*data(0)*b[0];
    appliedB[1] += loadFactor*data(1)*b[1];
    return 0;
  } else {
    opserr << "FourNodeQuad3d::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
  } 
  
  return -1;
}

int 
FourNodeQuad3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  int i;
  static double rhoi[4];
  double sum = 0.0;
  for (i = 0; i < 4; i++) {
    rhoi[i] = theMaterial[i]->getRho();
    sum += rhoi[i];
  }
  
  if (sum == 0.0)
    return 0;
  
  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  const Vector &Raccel3 = theNodes[2]->getRV(accel);
  const Vector &Raccel4 = theNodes[3]->getRV(accel);
  
  static double ra[12];
  
  ra[0] = Raccel1(0);
  ra[1] = Raccel1(1);
  ra[2] = Raccel1(2);
  ra[3] = Raccel2(0);
  ra[4] = Raccel2(1);
  ra[5] = Raccel2(2);
  ra[6] = Raccel3(0);
  ra[7] = Raccel3(1);
  ra[8] = Raccel3(2);
  ra[9] = Raccel4(0);
  ra[10] = Raccel4(1);
  ra[11] = Raccel4(2);
  
  // Compute mass matrix
  this->getMass();
  
  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  for (i = 0; i < 12; i++)
    Q(i) += -K(i,i)*ra[i];
  
  return 0;
}

const Vector&
FourNodeQuad3d::getResistingForce()
{
  P.Zero();
  
  double dvol;
  
  // Loop over the integration points
  for (int i = 0; i < 4; i++) {
    
    // Determine Jacobian for this integration point
    dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness*wts[i]);
    
    // Get material stress response
    const Vector &sigma = theMaterial[i]->getStress();
    
    // Perform numerical integration on internal force
    //P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
    //P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);

    int diff = dirn[1]-dirn[0];
    for (int alpha = 0, ia = dirn[0]; alpha < 4; alpha++, ia += 3) {
      
      P(ia) += dvol*(shp[0][alpha]*sigma(0) + shp[1][alpha]*sigma(2));
      
      P(ia+diff) += dvol*(shp[1][alpha]*sigma(1) + shp[0][alpha]*sigma(2));
      
      // Subtract equiv. body forces from the nodes
      //P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
      //P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
      if (applyLoad == 0) {
	P(ia) -= dvol*(shp[2][alpha]*b[0]);
	P(ia+diff) -= dvol*(shp[2][alpha]*b[1]);
      } else {
	P(ia) -= dvol*(shp[2][alpha]*appliedB[0]);
	P(ia+diff) -= dvol*(shp[2][alpha]*appliedB[1]);
      }
    }
  }
  
  // Subtract pressure loading from resisting force
  if (pressure != 0.0) {
    //P = P - pressureLoad;
    P.addVector(1.0, pressureLoad, -1.0);
  }
  
  // Subtract other external nodal loads ... P_res = P_int - P_ext
  //P = P - Q;
  P.addVector(1.0, Q, -1.0);
  
  return P;
}
 
const Vector&
FourNodeQuad3d::getResistingForceIncInertia()
{
  int i;
  static double rhoi[4];
  double sum = 0.0;
  for (i = 0; i < 4; i++) {
    rhoi[i] = theMaterial[i]->getRho();
    sum += rhoi[i];
  }
  
  // if no mass terms .. just add damping terms
  if (sum == 0.0) {
    this->getResistingForce();
    
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();
    
    return P;
  }
  
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();
  const Vector &accel3 = theNodes[2]->getTrialAccel();
  const Vector &accel4 = theNodes[3]->getTrialAccel();
  
  static double a[12];
  
  a[0] = accel1(0);
  a[1] = accel1(1);
  a[2] = accel1(2);
  a[3] = accel2(0);
  a[4] = accel2(1);
  a[5] = accel2(2);
  a[6] = accel3(0);
  a[7] = accel3(1);
  a[8] = accel3(2);
  a[9] = accel4(0);
  a[10] = accel4(1);
  a[11] = accel4(2);
  
  // Compute the current resisting force
  this->getResistingForce();
  
  // Compute the mass matrix
  this->getMass();
  
  // Take advantage of lumped mass matrix
  for (i = 0; i < 12; i++)
    P(i) += K(i,i)*a[i];
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P += this->getRayleighDampingForces();
  
  return P;
}

int
FourNodeQuad3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(10);
  data(0) = this->getTag();
  data(1) = thickness;
  data(3) = b[0];
  data(4) = b[1];
  data(5) = pressure;

  data(6) = alphaM;
  data(7) = betaK;
  data(8) = betaK0;
  data(9) = betaKc;
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad3d::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }	      
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(12);
  
  int i;
  for (i = 0; i < 4; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    matDbTag = theMaterial[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i+4) = matDbTag;
  }
  
  idData(8) = connectedExternalNodes(0);
  idData(9) = connectedExternalNodes(1);
  idData(10) = connectedExternalNodes(2);
  idData(11) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad3d::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING FourNodeQuad3d::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}

int
FourNodeQuad3d::recvSelf(int commitTag, Channel &theChannel,
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(10);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad3d::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  thickness = data(1);
  b[0] = data(3);
  b[1] = data(4);
  pressure = data(5);

  alphaM = data(6);
  betaK = data(7);
  betaK0 = data(8);
  betaKc = data(9);

  static ID idData(12);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad3d::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  connectedExternalNodes(0) = idData(8);
  connectedExternalNodes(1) = idData(9);
  connectedExternalNodes(2) = idData(10);
  connectedExternalNodes(3) = idData(11);
  

  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new NDMaterial *[4];
    if (theMaterial == 0) {
      opserr << "FourNodeQuad3d::recvSelf() - Could not allocate NDMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	opserr << "FourNodeQuad3d::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
	delete theMaterial[i];
	theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
	if (theMaterial[i] == 0) {
opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to create\n";
				
	  return -1;
	}
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }
  
  return res;
}

void
FourNodeQuad3d::Print(OPS_Stream &s, int flag)
{
    if (flag == 2) {

        s << "#FourNodeQuad3d\n";

        int i;
        const int numNodes = 4;
        const int nstress = 3;

        for (i = 0; i < numNodes; i++) {
            const Vector &nodeCrd = theNodes[i]->getCrds();
            const Vector &nodeDisp = theNodes[i]->getDisp();
            s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
        }

        // spit out the section location & invoke print on the scetion
        const int numMaterials = 4;

        static Vector avgStress(nstress);
        static Vector avgStrain(nstress);
        avgStress.Zero();
        avgStrain.Zero();
        for (i = 0; i < numMaterials; i++) {
            avgStress += theMaterial[i]->getStress();
            avgStrain += theMaterial[i]->getStrain();
        }
        avgStress /= numMaterials;
        avgStrain /= numMaterials;

        s << "#AVERAGE_STRESS ";
        for (i = 0; i < nstress; i++)
            s << avgStress(i) << " ";
        s << endln;

        s << "#AVERAGE_STRAIN ";
        for (i = 0; i < nstress; i++)
            s << avgStrain(i) << " ";
        s << endln;
    }
    
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nFourNodeQuad3d, element id:  " << this->getTag() << endln;
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tthickness:  " << thickness << endln;
        s << "\tsurface pressure:  " << pressure << endln;
        s << "\tmass density:  " << rho << endln;
        s << "\tbody forces:  " << b[0] << " " << b[1] << endln;
        theMaterial[0]->Print(s, flag);
        s << "\tStress (xx yy xy)" << endln;
        for (int i = 0; i < 4; i++)
            s << "\t\tGauss point " << i + 1 << ": " << theMaterial[i]->getStress();
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"FourNodeQuad3d\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        s << connectedExternalNodes(1) << ", ";
        s << connectedExternalNodes(2) << ", ";
        s << connectedExternalNodes(3) << "], ";
        s << "\"thickness\": " << thickness << ", ";
        s << "\"surfacePressure\": " << pressure << ", ";
        s << "\"masspervolume\": " << rho << ", ";
        s << "\"bodyForces\": [" << b[0] << ", " << b[1] << "], ";
        s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
    }
}

int
FourNodeQuad3d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **mdes, int numMode)
{
    // first set the quantity to be displayed at the nodes;
    // if displayMode is 1 through 3 we will plot material stresses otherwise 0.0

    static Vector values(4);

    for (int j=0; j<4; j++)
	   values(j) = 0.0;

    if (displayMode < 4 && displayMode > 0) {
	for (int i=0; i<4; i++) {
	  const Vector &stress = theMaterial[i]->getStress();
	  values(i) = stress(displayMode-1);
	}
    }

    // now  determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end3Crd = theNodes[2]->getCrds();	
    const Vector &end4Crd = theNodes[3]->getCrds();	

    static Matrix coords(4,3);

    if (displayMode >= 0) {    
      
      const Vector &end1Disp = theNodes[0]->getDisp();
      const Vector &end2Disp = theNodes[1]->getDisp();
      const Vector &end3Disp = theNodes[2]->getDisp();
      const Vector &end4Disp = theNodes[3]->getDisp();

      for (int i = 0; i < 3; i++) {
	coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
	coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
	coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
	coords(3,i) = end4Crd(i) + end4Disp(i)*fact;    
      }
    } else {
      int mode = displayMode * -1;
      const Matrix &eigen1 = theNodes[0]->getEigenvectors();
      const Matrix &eigen2 = theNodes[1]->getEigenvectors();
      const Matrix &eigen3 = theNodes[2]->getEigenvectors();
      const Matrix &eigen4 = theNodes[3]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < 3; i++) {
	  coords(0,i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	  coords(1,i) = end2Crd(i) + eigen2(i,mode-1)*fact;
	  coords(2,i) = end3Crd(i) + eigen3(i,mode-1)*fact;
	  coords(3,i) = end4Crd(i) + eigen4(i,mode-1)*fact;
	}    
      } else {
	for (int i = 0; i < 3; i++) {
	  coords(0,i) = end1Crd(i);
	  coords(1,i) = end2Crd(i);
	  coords(2,i) = end3Crd(i);
	  coords(3,i) = end4Crd(i);
	}    
      }
    }
    
    int error = 0;

    // finally we draw the element using drawPolygon
    error += theViewer.drawPolygon (coords, values);

    return error;
}

Response*
FourNodeQuad3d::setResponse(const char **argv, int argc, 
			  OPS_Stream &output)
{
  Response *theResponse =0;

  output.tag("ElementOutput");
  output.attr("eleType","FourNodeQuad3d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);

  char dataOut[10];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=4; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }
    
    theResponse =  new ElementResponse(this, 1, P);
  }   

  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {
      
      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",pts[pointNum-1][0]);
      output.attr("neta",pts[pointNum-1][1]);

      theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();

    }
  }
  else if ((strcmp(argv[0],"stresses") ==0) || (strcmp(argv[0],"stress") ==0)) {
    for (int i=0; i<4; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 3, Vector(12));
  }

  output.endTag(); // ElementOutput

  return theResponse;
}

int 
FourNodeQuad3d::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }
    return eleInfo.setVector(stresses);
	
  } else

    return -1;
}

int
FourNodeQuad3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  // quad pressure loading
  if (strcmp(argv[0],"pressure") == 0) {
    return param.addObject(2, this);
  }
  // a material parameter
  else if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {

    if (argc < 3)
      return -1;

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4)
      return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else 
      return -1;
  }

  // otherwise it could be just a forall material parameter
  else {

    int matRes = res;
    for (int i=0; i<4; i++) {

      matRes =  theMaterial[i]->setParameter(argv, argc, param);

      if (matRes != -1)
	res = matRes;
    }
  }
  
  return res;
}
    
int
FourNodeQuad3d::updateParameter(int parameterID, Information &info)
{
	int res = -1;
		int matRes = res;
  switch (parameterID) {
    case -1:
      return -1;

	case 1:
		
		for (int i = 0; i<4; i++) {
		matRes = theMaterial[i]->updateParameter(parameterID, info);
		}
		if (matRes != -1) {
			res = matRes;
		}
		return res;
      
	case 2:
		pressure = info.theDouble;
		this->setPressureLoadAtNodes();	// update consistent nodal loads
		return 0;

	default: 
	  /*	  
	  if (parameterID >= 100) { // material parameter
	    int pointNum = parameterID/100;
	    if (pointNum > 0 && pointNum <= 4)
	      return theMaterial[pointNum-1]->updateParameter(parameterID-100*pointNum, info);
	    else
	      return -1;
	  } else // unknown
	  */
	    return -1;
  }
}

double FourNodeQuad3d::shapeFunction(double xi, double eta)
{
	const Vector &nd1Crds = theNodes[0]->getCrds();
	const Vector &nd2Crds = theNodes[1]->getCrds();
	const Vector &nd3Crds = theNodes[2]->getCrds();
	const Vector &nd4Crds = theNodes[3]->getCrds();

	double oneMinuseta = 1.0-eta;
	double onePluseta = 1.0+eta;
	double oneMinusxi = 1.0-xi;
	double onePlusxi = 1.0+xi;

	shp[2][0] = 0.25*oneMinusxi*oneMinuseta;	// N_1
	shp[2][1] = 0.25*onePlusxi*oneMinuseta;		// N_2
	shp[2][2] = 0.25*onePlusxi*onePluseta;		// N_3
	shp[2][3] = 0.25*oneMinusxi*onePluseta;		// N_4

	double J[2][2];

	J[0][0] = 0.25 * (-nd1Crds(dirn[0])*oneMinuseta + nd2Crds(dirn[0])*oneMinuseta +
				nd3Crds(dirn[0])*(onePluseta) - nd4Crds(dirn[0])*(onePluseta));

	J[0][1] = 0.25 * (-nd1Crds(dirn[0])*oneMinusxi - nd2Crds(0)*onePlusxi +
				nd3Crds(dirn[0])*onePlusxi + nd4Crds(dirn[0])*oneMinusxi);

	J[1][0] = 0.25 * (-nd1Crds(dirn[1])*oneMinuseta + nd2Crds(dirn[1])*oneMinuseta +
				nd3Crds(dirn[1])*onePluseta - nd4Crds(dirn[1])*onePluseta);

	J[1][1] = 0.25 * (-nd1Crds(dirn[1])*oneMinusxi - nd2Crds(dirn[1])*onePlusxi +
				nd3Crds(dirn[1])*onePlusxi + nd4Crds(dirn[1])*oneMinusxi);

	double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
	double oneOverdetJ = 1.0/detJ;
	double L[2][2];

	// L = inv(J)
	L[0][0] =  J[1][1]*oneOverdetJ;
	L[1][0] = -J[0][1]*oneOverdetJ;
	L[0][1] = -J[1][0]*oneOverdetJ;
	L[1][1] =  J[0][0]*oneOverdetJ;

    double L00 = 0.25*L[0][0];
    double L10 = 0.25*L[1][0];
    double L01 = 0.25*L[0][1];
    double L11 = 0.25*L[1][1];
	
	double L00oneMinuseta = L00*oneMinuseta;
	double L00onePluseta  = L00*onePluseta;
	double L01oneMinusxi  = L01*oneMinusxi;
	double L01onePlusxi   = L01*onePlusxi;

	double L10oneMinuseta = L10*oneMinuseta;
	double L10onePluseta  = L10*onePluseta;
	double L11oneMinusxi  = L11*oneMinusxi;
	double L11onePlusxi   = L11*onePlusxi;

	// See Cook, Malkus, Plesha p. 169 for the derivation of these terms
    shp[0][0] = -L00oneMinuseta - L01oneMinusxi;	// N_1,1
    shp[0][1] =  L00oneMinuseta - L01onePlusxi;		// N_2,1
    shp[0][2] =  L00onePluseta  + L01onePlusxi;		// N_3,1
    shp[0][3] = -L00onePluseta  + L01oneMinusxi;	// N_4,1
	
    shp[1][0] = -L10oneMinuseta - L11oneMinusxi;	// N_1,2
    shp[1][1] =  L10oneMinuseta - L11onePlusxi;		// N_2,2
    shp[1][2] =  L10onePluseta  + L11onePlusxi;		// N_3,2
    shp[1][3] = -L10onePluseta  + L11oneMinusxi;	// N_4,2

    return detJ;
}

void 
FourNodeQuad3d::setPressureLoadAtNodes(void)
{
        pressureLoad.Zero();

	if (pressure == 0.0)
		return;

	const Vector &node1 = theNodes[0]->getCrds();
	const Vector &node2 = theNodes[1]->getCrds();
	const Vector &node3 = theNodes[2]->getCrds();
	const Vector &node4 = theNodes[3]->getCrds();

	double x1 = node1(dirn[0]);
	double y1 = node1(dirn[1]);
	double x2 = node2(dirn[0]);
	double y2 = node2(dirn[1]);
	double x3 = node3(dirn[0]);
	double y3 = node3(dirn[1]);
	double x4 = node4(dirn[0]);
	double y4 = node4(dirn[1]);

	double dx12 = x2-x1;
	double dy12 = y2-y1;
	double dx23 = x3-x2;
	double dy23 = y3-y2;
	double dx34 = x4-x3;
	double dy34 = y4-y3;
	double dx41 = x1-x4;
	double dy41 = y1-y4;

	double pressureOver2 = pressure/2.0;

	// Contribution from side 12
	pressureLoad(dirn[0]) += pressureOver2*dy12;
	pressureLoad(dirn[0]+3) += pressureOver2*dy12;
	pressureLoad(dirn[1]) += pressureOver2*-dx12;
	pressureLoad(dirn[1]+3) += pressureOver2*-dx12;

	// Contribution from side 23
	pressureLoad(dirn[0]+3) += pressureOver2*dy23;
	pressureLoad(dirn[0]+6) += pressureOver2*dy23;
	pressureLoad(dirn[1]+3) += pressureOver2*-dx23;
	pressureLoad(dirn[1]+6) += pressureOver2*-dx23;

	// Contribution from side 34
	pressureLoad(dirn[0]+6) += pressureOver2*dy34;
	pressureLoad(dirn[0]+9) += pressureOver2*dy34;
	pressureLoad(dirn[1]+6) += pressureOver2*-dx34;
	pressureLoad(dirn[1]+9) += pressureOver2*-dx34;

	// Contribution from side 41
	pressureLoad(dirn[0]+9) += pressureOver2*dy41;
	pressureLoad(dirn[0]) += pressureOver2*dy41;
	pressureLoad(dirn[1]+9) += pressureOver2*-dx41;
	pressureLoad(dirn[1]) += pressureOver2*-dx41;

	//pressureLoad = pressureLoad*thickness;
}








