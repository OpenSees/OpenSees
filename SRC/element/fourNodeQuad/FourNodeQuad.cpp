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
                                                                        
// $Revision: 1.23 $
// $Date: 2003-02-25 23:32:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/FourNodeQuad.cpp,v $

// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
// Description: This file contains the class definition for FourNodeQuad.

#include <FourNodeQuad.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>


double FourNodeQuad::matrixData[64];
Matrix FourNodeQuad::K(matrixData, 8, 8);
Vector FourNodeQuad::P(8);
double FourNodeQuad::shp[3][4];
double FourNodeQuad::pts[4][2];
double FourNodeQuad::wts[4];

FourNodeQuad::FourNodeQuad(int tag, int nd1, int nd2, int nd3, int nd4,
			   NDMaterial &m, const char *type, double t,
			   double p, double r, double b1, double b2)
:Element (tag, ELE_TAG_FourNodeQuad), 
  theMaterial(0), connectedExternalNodes(4), 
  Q(8), pressureLoad(8), thickness(t), rho(r), pressure(p), Ki(0)
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
	  opserr << "FourNodeQuad::FourNodeQuad -- improper material type: " << type << "for FourNodeQuad\n";
	  exit(-1);
	}

	// Body forces
	b[0] = b1;
	b[1] = b2;

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new NDMaterial *[4];
    
    if (theMaterial == 0) {
      opserr << "FourNodeQuad::FourNodeQuad - failed allocate material model pointer\n";
      exit(-1);
    }

	int i;
    for (i = 0; i < 4; i++) {
      
      // Get copies of the material model for each integration point
      theMaterial[i] = m.getCopy(type);
			
      // Check allocation
      if (theMaterial[i] == 0) {
	opserr << "FourNodeQuad::FourNodeQuad -- failed to get a copy of material model\n";
	exit(-1);
      }
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
    
    for (i=0; i<4; i++)
      theNodes[i] = 0;
}

FourNodeQuad::FourNodeQuad()
:Element (0,ELE_TAG_FourNodeQuad),
  theMaterial(0), connectedExternalNodes(4), 
 Q(8), pressureLoad(8), thickness(0.0), rho(0.0), pressure(0.0), Ki(0)
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

FourNodeQuad::~FourNodeQuad()
{    
  for (int i = 0; i < 4; i++) {
    if (theMaterial[i])
      delete theMaterial[i];
  }

  // Delete the array of pointers to NDMaterial pointer arrays
  if (theMaterial)
    delete [] theMaterial;

  if (Ki != 0)
    delete Ki;
}

int
FourNodeQuad::getNumExternalNodes() const
{
    return 4;
}

const ID&
FourNodeQuad::getExternalNodes()
{
    return connectedExternalNodes;
}


Node **
FourNodeQuad::getNodePtrs(void) 
{
  return theNodes;
}

int
FourNodeQuad::getNumDOF()
{
    return 8;
}

void
FourNodeQuad::setDomain(Domain *theDomain)
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
	//opserr << "FATAL ERROR FourNodeQuad (tag: %d), node not found in domain",
	//	this->getTag());
	
	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
    
    if (dofNd1 != 2 || dofNd2 != 2 || dofNd3 != 2 || dofNd4 != 2) {
	//opserr << "FATAL ERROR FourNodeQuad (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());
	
	return;
    }
    this->DomainComponent::setDomain(theDomain);

    // Compute consistent nodal loads due to pressure
    this->setPressureLoadAtNodes();
}

int
FourNodeQuad::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "FourNodeQuad::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < 4; i++)
      retVal += theMaterial[i]->commitState();

    return retVal;
}

int
FourNodeQuad::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}

int
FourNodeQuad::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
FourNodeQuad::update()
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();
	
	static double u[2][4];

	u[0][0] = disp1(0);
	u[1][0] = disp1(1);
	u[0][1] = disp2(0);
	u[1][1] = disp2(1);
	u[0][2] = disp3(0);
	u[1][2] = disp3(1);
	u[0][3] = disp4(0);
	u[1][3] = disp4(1);

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
FourNodeQuad::getTangentStiff()
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
	    
	  for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
	    for (int beta = 0, ib = 0; beta < 4; beta++, ib += 2) {
	      
	      DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
	      DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
	      DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
	      DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
	      DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
	      DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);
	      

	      K(ia,ib) += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	      K(ia,ib+1) += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	      K(ia+1,ib) += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	      K(ia+1,ib+1) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	      //	      matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	      //matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	      //matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	      //matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	      
	    }
	  }
	}
	
	return K;
}


const Matrix&
FourNodeQuad::getInitialStiff()
{
  if (Ki != 0)
    return *Ki;

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
	
	matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
      }
    }
  }

  Ki = new Matrix(K);
  return K;
}

const Matrix&
FourNodeQuad::getMass()
{
	K.Zero();

	int i;
	static double rhoi[4];
	double sum = this->rho;
	for (i = 0; i < 4; i++) {
	  rhoi[i] = theMaterial[i]->getRho();
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
		double tmp = rho + rhoi[i];
		rhodvol *= (tmp*thickness*wts[i]);

		for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia++) {
			Nrho = shp[2][alpha]*rhodvol;
			K(ia,ia) += Nrho;
			ia++;
			K(ia,ia) += Nrho;
		}
	}

	return K;
}

void
FourNodeQuad::zeroLoad(void)
{
  Q.Zero();
  return;
}

int 
FourNodeQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "FourNodeQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

int 
FourNodeQuad::addInertiaLoadToUnbalance(const Vector &accel)
{
  int i;
  static double rhoi[4];
  double sum = this->rho;
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
  
  if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size() ||
      2 != Raccel4.Size()) {
    opserr << "FourNodeQuad::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
    return -1;
  }
  
  static double ra[8];
  
  ra[0] = Raccel1(0);
  ra[1] = Raccel1(1);
  ra[2] = Raccel2(0);
  ra[3] = Raccel2(1);
  ra[4] = Raccel3(0);
  ra[5] = Raccel3(1);
  ra[6] = Raccel4(0);
  ra[7] = Raccel4(1);
  
  // Compute mass matrix
  this->getMass();
  
  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  for (i = 0; i < 8; i++)
    Q(i) += -K(i,i)*ra[i];
  
  return 0;
}

const Vector&
FourNodeQuad::getResistingForce()
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
		for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
			
			P(ia) += dvol*(shp[0][alpha]*sigma(0) + shp[1][alpha]*sigma(2));
			
			P(ia+1) += dvol*(shp[1][alpha]*sigma(1) + shp[0][alpha]*sigma(2));

			// Subtract equiv. body forces from the nodes
			//P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
			//P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
			P(ia) -= dvol*(shp[2][alpha]*b[0]);
			P(ia+1) -= dvol*(shp[2][alpha]*b[1]);
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
FourNodeQuad::getResistingForceIncInertia()
{
	int i;
	static double rhoi[4];
	double sum = this->rho;
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
	
	static double a[8];

	a[0] = accel1(0);
	a[1] = accel1(1);
	a[2] = accel2(0);
	a[3] = accel2(1);
	a[4] = accel3(0);
	a[5] = accel3(1);
	a[6] = accel4(0);
	a[7] = accel4(1);

	// Compute the current resisting force
	this->getResistingForce();

	// Compute the mass matrix
	this->getMass();

	// Take advantage of lumped mass matrix
	for (i = 0; i < 8; i++)
		P(i) += K(i,i)*a[i];

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	  P += this->getRayleighDampingForces();

	return P;
}

int
FourNodeQuad::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = thickness;
  data(2) = rho;
  data(3) = b[0];
  data(4) = b[1];
  data(5) = pressure;
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
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
    opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}

int
FourNodeQuad::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(6);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  thickness = data(1);
  rho = data(2);
  b[0] = data(3);
  b[1] = data(4);
  pressure = data(5);

  static ID idData(12);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
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
      opserr << "FourNodeQuad::recvSelf() - Could not allocate NDMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	opserr << "FourNodeQuad::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
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
FourNodeQuad::Print(OPS_Stream &s, int flag)
{
	s << "\nFourNodeQuad, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tthickness:  " << thickness << endln;
	s << "\tmass density:  " << rho << endln;
	s << "\tsurface pressure:  " << pressure << endln;
	s << "\tbody forces:  " << b[0] << " " << b[1] << endln;
	theMaterial[0]->Print(s,flag);
	s << "\tStress (xx yy xy)" << endln;
	for (int i = 0; i < 4; i++)
		s << "\t\tGauss point " << i+1 << ": " << theMaterial[i]->getStress();
}

int
FourNodeQuad::displaySelf(Renderer &theViewer, int displayMode, float fact)
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

    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();
    const Vector &end3Disp = theNodes[2]->getDisp();
    const Vector &end4Disp = theNodes[3]->getDisp();

    static Matrix coords(4,3);

    for (int i = 0; i < 2; i++) {
      coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
      coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
      coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
      coords(3,i) = end4Crd(i) + end4Disp(i)*fact;    
    }

    int error = 0;

    // finally we draw the element using drawPolygon
    error += theViewer.drawPolygon (coords, values);

    return error;
}

Response*
FourNodeQuad::setResponse(const char **argv, int argc, Information &eleInfo)
{
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
      return new ElementResponse(this, 1, P);
    
    else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
      return new ElementResponse(this, 2, K);

    else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
      int pointNum = atoi(argv[1]);
      if (pointNum > 0 && pointNum <= 4)
	return theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, eleInfo);
      else 
	return 0;
    } else if (strcmp(argv[0],"stresses") ==0) {
      return new ElementResponse(this, 3, P);
    }
 
    // otherwise response quantity is unknown for the quad class
    else
      return 0;
}

int 
FourNodeQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 2) {

    return eleInfo.setMatrix(this->getTangentStiff());

  } else if (responseID == 3) {

    // Loop over the integration points
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      P(cnt) = sigma(0);
      P(cnt+1) = sigma(1);
      cnt += 2;
    }
    return eleInfo.setVector(P);
	
  } else

    return -1;
}

int
FourNodeQuad::setParameter(const char **argv, int argc, Information &info)
{
	// quad mass density per unit volume
	if (strcmp(argv[0],"rho") == 0) {
		info.theType = DoubleType;
		info.theDouble = rho;
		return 1;
	}
	// quad pressure loading
	if (strcmp(argv[0],"pressure") == 0) {
		info.theType = DoubleType;
		info.theDouble = pressure;
		return 2;
	}
    // a material parameter
    else if (strcmp(argv[0],"material") == 0) {
		int pointNum = atoi(argv[1]);
		if (pointNum > 0 && pointNum <= 4) {
			int ok = theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, info);
			if (ok < 0)
				return -1;
		    else if (ok >= 0 && ok < 100)
				return pointNum*100 + ok;
			else
				return -1;
		}
	    else 
			return -1;
	}
    
    // otherwise parameter is unknown for the FourNodeQuad class
    else
      return -1;

}
    
int
FourNodeQuad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
    case -1:
      return -1;
      
	case 1:
		rho = info.theDouble;
		this->getMass();	// update mass matrix
		return 0;
	case 2:
		pressure = info.theDouble;
		this->setPressureLoadAtNodes();	// update consistent nodal loads
		return 0;
	default: 
		if (parameterID >= 100) { // material parameter
			int pointNum = parameterID/100;
			if (pointNum > 0 && pointNum <= 4)
				return theMaterial[pointNum-1]->updateParameter(parameterID-100*pointNum, info);
			else
				return -1;
		} else // unknown
			return -1;
  }
}

double FourNodeQuad::shapeFunction(double xi, double eta)
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

	J[0][0] = 0.25 * (-nd1Crds(0)*oneMinuseta + nd2Crds(0)*oneMinuseta +
				nd3Crds(0)*(onePluseta) - nd4Crds(0)*(onePluseta));

	J[0][1] = 0.25 * (-nd1Crds(0)*oneMinusxi - nd2Crds(0)*onePlusxi +
				nd3Crds(0)*onePlusxi + nd4Crds(0)*oneMinusxi);

	J[1][0] = 0.25 * (-nd1Crds(1)*oneMinuseta + nd2Crds(1)*oneMinuseta +
				nd3Crds(1)*onePluseta - nd4Crds(1)*onePluseta);

	J[1][1] = 0.25 * (-nd1Crds(1)*oneMinusxi - nd2Crds(1)*onePlusxi +
				nd3Crds(1)*onePlusxi + nd4Crds(1)*oneMinusxi);

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
FourNodeQuad::setPressureLoadAtNodes(void)
{
        pressureLoad.Zero();

	if (pressure == 0.0)
		return;

	const Vector &node1 = theNodes[0]->getCrds();
	const Vector &node2 = theNodes[1]->getCrds();
	const Vector &node3 = theNodes[2]->getCrds();
	const Vector &node4 = theNodes[3]->getCrds();

	double x1 = node1(0);
	double y1 = node1(1);
	double x2 = node2(0);
	double y2 = node2(1);
	double x3 = node3(0);
	double y3 = node3(1);
	double x4 = node4(0);
	double y4 = node4(1);

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
	pressureLoad(0) += pressureOver2*dy12;
	pressureLoad(2) += pressureOver2*dy12;
	pressureLoad(1) += pressureOver2*-dx12;
	pressureLoad(3) += pressureOver2*-dx12;

	// Contribution from side 23
	pressureLoad(2) += pressureOver2*dy23;
	pressureLoad(4) += pressureOver2*dy23;
	pressureLoad(3) += pressureOver2*-dx23;
	pressureLoad(5) += pressureOver2*-dx23;

	// Contribution from side 34
	pressureLoad(4) += pressureOver2*dy34;
	pressureLoad(6) += pressureOver2*dy34;
	pressureLoad(5) += pressureOver2*-dx34;
	pressureLoad(7) += pressureOver2*-dx34;

	// Contribution from side 41
	pressureLoad(6) += pressureOver2*dy41;
	pressureLoad(0) += pressureOver2*dy41;
	pressureLoad(7) += pressureOver2*-dx41;
	pressureLoad(1) += pressureOver2*-dx41;

	//pressureLoad = pressureLoad*thickness;
}








