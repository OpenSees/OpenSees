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
//
// based on FourNodeQuad by MHS
// Written: Seweryn Kokot, Opole University of Technology, Poland
// Created: Sep 2020
//
// Description: This file contains the class definition for SixNodeTri.

#include "SixNodeTri.h"
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

void* OPS_SixNodeTri()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 2 || ndf != 2) {
	opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
	return 0;
    }

    if (OPS_GetNumRemainingInputArgs() < 10) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: element SixNodeTri eleTag? iNode? jNode? kNode? lNode? nNode5 mNode6 thk? type? matTag? <pressure? rho? b1? b2?>\n";
	return 0;
    }

    // SixNodeTriId, iNode, jNode, kNode, lNode
	// nNode, mNode
    int idata[7];
    int num = 7;
    if (OPS_GetIntInput(&num,idata) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    double thk = 1.0;
    num = 1;
    if (OPS_GetDoubleInput(&num,&thk) < 0) {
	opserr<<"WARNING: invalid double inputs\n";
	return 0;
    }

    const char* type = OPS_GetString();

    int matTag;
    num = 1;
    if (OPS_GetIntInput(&num,&matTag) < 0) {
	opserr<<"WARNING: invalid matTag\n";
	return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(matTag);
    if (mat == 0) {
	opserr << "WARNING material not found\n";
	opserr << "Material: " << matTag;
	opserr << "\nSixNodeTri element: " << idata[0] << endln;
	return 0;
    }

    // p, rho, b1, b2
    double data[4] = {0,0,0,0};
    num = OPS_GetNumRemainingInputArgs();
    if (num > 4) {
	num = 4;
    }
    if (num > 0) {
	if (OPS_GetDoubleInput(&num,data) < 0) {
	    opserr<<"WARNING: invalid integer data\n";
	    return 0;
	}
    }

    return new SixNodeTri(idata[0],idata[1],idata[2],idata[3],idata[4],
							idata[5],idata[6],
			                *mat,type,thk,data[0],data[1],data[2],data[3]);
}


double SixNodeTri::matrixData[(nnodes*2)*(nnodes*2)];
Matrix SixNodeTri::K(matrixData, 2*nnodes, 2*nnodes);
Vector SixNodeTri::P(2*nnodes);
double SixNodeTri::shp[3][nnodes];
double SixNodeTri::pts[nip][2];
double SixNodeTri::wts[nip];

SixNodeTri::SixNodeTri(int tag, int nd1, int nd2, int nd3, int nd4,
						   int nd5, int nd6,
						   NDMaterial &m, const char *type, double t,
						   double p, double r, double b1, double b2)
:Element (tag, ELE_TAG_SixNodeTri),
  theMaterial(0), connectedExternalNodes(nnodes),
 Q(2*nnodes), applyLoad(0), pressureLoad(2*nnodes), thickness(t), pressure(p), rho(r), Ki(0)
{
	pts[0][0] = 0.666666666666666667;
	pts[0][1] = 0.166666666666666667;
	pts[1][0] = 0.166666666666666667;
	pts[1][1] = 0.666666666666666667;
	pts[2][0] = 0.166666666666666667;
	pts[2][1] = 0.166666666666666667;
	// pts[0][0] = 0.5;
	// pts[0][1] = 0.5;
	// pts[1][0] = 0.0;
	// pts[1][1] = 0.5;
	// pts[2][0] = 0.5;
	// pts[2][1] = 0.0;

	wts[0] = 0.166666666666666667;
	wts[1] = 0.166666666666666667;
	wts[2] = 0.166666666666666667;

	if (strcmp(type,"PlaneStrain") != 0 && strcmp(type,"PlaneStress") != 0
	    && strcmp(type,"PlaneStrain2D") != 0 && strcmp(type,"PlaneStress2D") != 0) {
	  opserr << "SixNodeTri::SixNodeTri -- improper material type: " << type << "for SixNodeTri\n";
	  exit(-1);
	}

	// Body forces
	b[0] = b1;
	b[1] = b2;

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new NDMaterial *[nip];

    if (theMaterial == 0) {
      opserr << "SixNodeTri::SixNodeTri - failed allocate material model pointer\n";
      exit(-1);
    }

	int i;
    for (i = 0; i < nip; i++) {

      // Get copies of the material model for each integration point
      theMaterial[i] = m.getCopy(type);

      // Check allocation
      if (theMaterial[i] == 0) {
	opserr << "SixNodeTri::SixNodeTri -- failed to get a copy of material model\n";
	exit(-1);
      }
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
    connectedExternalNodes(4) = nd5;
    connectedExternalNodes(5) = nd6;

    for (i=0; i<nnodes; i++)
      theNodes[i] = 0;
}

SixNodeTri::SixNodeTri()
:Element (0,ELE_TAG_SixNodeTri),
  theMaterial(0), connectedExternalNodes(nnodes),
 Q(2*nnodes), applyLoad(0), pressureLoad(2*nnodes), thickness(0.0), pressure(0.0), Ki(0)
{
	pts[0][0] = 0.666666666666666667;
	pts[0][1] = 0.166666666666666667;
	pts[1][0] = 0.166666666666666667;
	pts[1][1] = 0.666666666666666667;
	pts[2][0] = 0.166666666666666667;
	pts[2][1] = 0.166666666666666667;
	// pts[0][0] = 0.5;
	// pts[0][1] = 0.5;
	// pts[1][0] = 0.0;
	// pts[1][1] = 0.5;
	// pts[2][0] = 0.5;
	// pts[2][1] = 0.0;

	wts[0] = 0.166666666666666667;
	wts[1] = 0.166666666666666667;
	wts[2] = 0.166666666666666667;

    for (int i=0; i<nnodes; i++)
      theNodes[i] = 0;
}

SixNodeTri::~SixNodeTri()
{
  for (int i = 0; i < nip; i++) {
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
SixNodeTri::getNumExternalNodes() const
{
    return nnodes;
}

const ID&
SixNodeTri::getExternalNodes()
{
    return connectedExternalNodes;
}


Node **
SixNodeTri::getNodePtrs(void)
{
  return theNodes;
}

int
SixNodeTri::getNumDOF()
{
    return nnodes*2;
}

void
SixNodeTri::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;
	theNodes[4] = 0;
	theNodes[5] = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);
    int Nd5 = connectedExternalNodes(4);
    int Nd6 = connectedExternalNodes(5);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(Nd3);
    theNodes[3] = theDomain->getNode(Nd4);
    theNodes[4] = theDomain->getNode(Nd5);
    theNodes[5] = theDomain->getNode(Nd6);

    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0 ||
		theNodes[4] == 0 || theNodes[5] == 0) {
	//opserr << "FATAL ERROR SixNodeTri (tag: %d), node not found in domain",
	//	this->getTag());

	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
    int dofNd5 = theNodes[4]->getNumberDOF();
    int dofNd6 = theNodes[5]->getNumberDOF();

    if (dofNd1 != 2 || dofNd2 != 2 || dofNd3 != 2 || dofNd4 != 2 ||
		dofNd5 != 2 || dofNd6 != 2) {
	//opserr << "FATAL ERROR SixNodeTri (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());

	return;
    }
    this->DomainComponent::setDomain(theDomain);

    // Compute consistent nodal loads due to pressure
    this->setPressureLoadAtNodes();
}

int
SixNodeTri::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "SixNodeTri::commitState () - failed in base class";
    }

    // Loop over the integration points and commit the material states
    for (int i = 0; i < nip; i++)
      retVal += theMaterial[i]->commitState();

    return retVal;
}

int
SixNodeTri::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < nip; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}

int
SixNodeTri::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < nip; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
SixNodeTri::update()
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();
	const Vector &disp5 = theNodes[4]->getTrialDisp();
	const Vector &disp6 = theNodes[5]->getTrialDisp();

	static double u[2][nnodes];

	u[0][0] = disp1(0);
	u[1][0] = disp1(1);
	u[0][1] = disp2(0);
	u[1][1] = disp2(1);
	u[0][2] = disp3(0);
	u[1][2] = disp3(1);
	u[0][3] = disp4(0);
	u[1][3] = disp4(1);
	u[0][4] = disp5(0);
	u[1][4] = disp5(1);
	u[0][5] = disp6(0);
	u[1][5] = disp6(1);

	static Vector eps(3);

	int ret = 0;

	// Loop over the integration points
	for (int i = 0; i < nip; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(pts[i][0], pts[i][1]);

		// Interpolate strains
		//eps = B*u;
		//eps.addMatrixVector(0.0, B, u, 1.0);
		eps.Zero();
		for (int beta = 0; beta < nnodes; beta++) {
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
SixNodeTri::getTangentStiff()
{

	K.Zero();

	double dvol;
	double DB[3][2];

	// Loop over the integration points
	for (int i = 0; i < nip; i++) {

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

	  for (int alpha = 0, ia = 0; alpha < nnodes; alpha++, ia += 2) {
	    for (int beta = 0, ib = 0; beta < nnodes; beta++, ib += 2) {

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
SixNodeTri::getInitialStiff()
{
  if (Ki != 0)
    return *Ki;

  K.Zero();

  double dvol;
  double DB[3][2];

  // Loop over the integration points
  for (int i = 0; i < nip; i++) {

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
	 beta < nnodes;
	 beta++, ib += 2, colIb += 16, colIbP1 += 16) {

      for (int alpha = 0, ia = 0; alpha < nnodes; alpha++, ia += 2) {

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
SixNodeTri::getMass()
{
	K.Zero();

	int i;
	static double rhoi[3]; // nip
	double sum = 0.0;
	for (i = 0; i < nip; i++) {
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
	for (i = 0; i < nip; i++) {

		// Determine Jacobian for this integration point
		rhodvol = this->shapeFunction(pts[i][0], pts[i][1]);

		// Element plus material density ... MAY WANT TO REMOVE ELEMENT DENSITY
		rhodvol *= (rhoi[i]*thickness*wts[i]);

		for (int alpha = 0, ia = 0; alpha < nnodes; alpha++, ia++) {
			Nrho = shp[2][alpha]*rhodvol;
			K(ia,ia) += Nrho;
			ia++;
			K(ia,ia) += Nrho;
		}
	}

	return K;
}

void
SixNodeTri::zeroLoad(void)
{
	Q.Zero();

	applyLoad = 0;

	appliedB[0] = 0.0;
	appliedB[1] = 0.0;

  	return;
}

int
SixNodeTri::addLoad(ElementalLoad *theLoad, double loadFactor)
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
		opserr << "SixNodeTri::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	}

	return -1;
}

int
SixNodeTri::addInertiaLoadToUnbalance(const Vector &accel)
{
  int i;
  static double rhoi[nip];
  double sum = 0.0;
  for (i = 0; i < nip; i++) {
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
  const Vector &Raccel5 = theNodes[4]->getRV(accel);
  const Vector &Raccel6 = theNodes[5]->getRV(accel);

  if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size() ||
      2 != Raccel4.Size() || 2 != Raccel5.Size() || 2 != Raccel6.Size()) {
    opserr << "SixNodeTri::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }

  static double ra[2*nnodes];

  ra[0] = Raccel1(0);
  ra[1] = Raccel1(1);
  ra[2] = Raccel2(0);
  ra[3] = Raccel2(1);
  ra[4] = Raccel3(0);
  ra[5] = Raccel3(1);
  ra[6] = Raccel4(0);
  ra[7] = Raccel4(1);
  ra[8] = Raccel5(0);
  ra[9] = Raccel5(1);
  ra[10] = Raccel6(0);
  ra[11] = Raccel6(1);

  // Compute mass matrix
  this->getMass();

  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  for (i = 0; i < 2*nnodes; i++)
    Q(i) += -K(i,i)*ra[i];

  return 0;
}

const Vector&
SixNodeTri::getResistingForce()
{
	P.Zero();

	double dvol;

	// Loop over the integration points
	for (int i = 0; i < nip; i++) {

		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1]);
		dvol *= (thickness*wts[i]);

		// Get material stress response
		const Vector &sigma = theMaterial[i]->getStress();

		// Perform numerical integration on internal force
		//P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
		//P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);
		for (int alpha = 0, ia = 0; alpha < nnodes; alpha++, ia += 2) {

			P(ia) += dvol*(shp[0][alpha]*sigma(0) + shp[1][alpha]*sigma(2));

			P(ia+1) += dvol*(shp[1][alpha]*sigma(1) + shp[0][alpha]*sigma(2));

			// Subtract equiv. body forces from the nodes
			//P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
			//P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
			if (applyLoad == 0) {
				P(ia) -= dvol*(shp[2][alpha]*b[0]);
				P(ia+1) -= dvol*(shp[2][alpha]*b[1]);
			} else {
				P(ia) -= dvol*(shp[2][alpha]*appliedB[0]);
				P(ia+1) -= dvol*(shp[2][alpha]*appliedB[1]);
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
SixNodeTri::getResistingForceIncInertia()
{
	int i;
	static double rhoi[nip];
	double sum = 0.0;
	for (i = 0; i < nip; i++) {
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
	const Vector &accel5 = theNodes[4]->getTrialAccel();
	const Vector &accel6 = theNodes[5]->getTrialAccel();

	static double a[2*nnodes];

	a[0] = accel1(0);
	a[1] = accel1(1);
	a[2] = accel2(0);
	a[3] = accel2(1);
	a[4] = accel3(0);
	a[5] = accel3(1);
	a[6] = accel4(0);
	a[7] = accel4(1);
	a[8] = accel5(0);
	a[9] = accel5(1);
	a[10] = accel6(0);
	a[11] = accel6(1);

	// Compute the current resisting force
	this->getResistingForce();

	// Compute the mass matrix
	this->getMass();

	// Take advantage of lumped mass matrix
	for (i = 0; i < 2*nnodes; i++)
		P(i) += K(i,i)*a[i];

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	  P += this->getRayleighDampingForces();

	return P;
}

int
SixNodeTri::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(9);
  data(0) = this->getTag();
  data(1) = thickness;
  data(2) = b[0];
  data(3) = b[1];
  data(4) = pressure;

  data(5) = alphaM;
  data(6) = betaK;
  data(7) = betaK0;
  data(8) = betaKc;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING SixNodeTri::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }


  // Now quad sends the ids of its materials
  int matDbTag;

  static ID idData(2*nip+nnodes);

  int i;
  for (i = 0; i < nip; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    matDbTag = theMaterial[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  theMaterial[i]->setDbTag(matDbTag);
    }
    // idData(i+4) = matDbTag;
    idData(i+nip) = matDbTag;
  }

  for( i = 0; i < nnodes; i++)
	idData(2*nip+i) = connectedExternalNodes(i);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING SixNodeTri::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < nip; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING SixNodeTri::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}

int
SixNodeTri::recvSelf(int commitTag, Channel &theChannel,
                       FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the
  // internal data with the data in the Vector
  static Vector data(9);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING SixNodeTri::recvSelf() - failed to receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  thickness = data(1);
  b[0] = data(2);
  b[1] = data(3);
  pressure = data(4);

  alphaM = data(5);
  betaK = data(6);
  betaK0 = data(7);
  betaKc = data(8);

  static ID idData(18); // 2*3 + 6
  // Quad now receives the tags of its nine external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING SixNodeTri::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  for( int i = 0; i < nnodes; i++)
    connectedExternalNodes(i) = idData(2*nip+i);

  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new NDMaterial *[nip];
    if (theMaterial == 0) {
      opserr << "SixNodeTri::recvSelf() - Could not allocate NDMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < nip; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+nip);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	    opserr << "SixNodeTri::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	    return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "SixNodeTri::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < nip; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+nip);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
	    delete theMaterial[i];
	    theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
	    if (theMaterial[i] == 0) {
          opserr << "SixNodeTri::recvSelf() - material " << i << "failed to create\n";
	      return -1;
	    }
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "SixNodeTri::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }

  return res;
}

void
SixNodeTri::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

    s << "#SixNodeTri\n";

    int i;
    const int numNodes = nnodes;
    const int nstress = nip ;

    for (i=0; i<numNodes; i++) {
      const Vector &nodeCrd = theNodes[i]->getCrds();
      // const Vector &nodeDisp = theNodes[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
     }

    // spit out the section location & invoke print on the scetion
    const int numMaterials = nip;

    static Vector avgStress(nstress);
    static Vector avgStrain(nstress);
    avgStress.Zero();
    avgStrain.Zero();
    for (i=0; i<numMaterials; i++) {
      avgStress += theMaterial[i]->getStress();
      avgStrain += theMaterial[i]->getStrain();
    }
    avgStress /= numMaterials;
    avgStrain /= numMaterials;

    s << "#AVERAGE_STRESS ";
    for (i=0; i<nstress; i++)
      s << avgStress(i) << " " ;
    s << endln;

    s << "#AVERAGE_STRAIN ";
    for (i=0; i<nstress; i++)
      s << avgStrain(i) << " " ;
    s << endln;
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nSixNodeTri, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tthickness:  " << thickness << endln;
	s << "\tsurface pressure:  " << pressure << endln;
    s << "\tmass density:  " << rho << endln;
    s << "\tbody forces:  " << b[0] << " " << b[1] << endln;
	theMaterial[0]->Print(s,flag);
	s << "\tStress (xx yy xy)" << endln;
	for (int i = 0; i < nip; i++)
		s << "\t\tGauss point " << i+1 << ": " << theMaterial[i]->getStress();
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << "\t\t\t{";
      s << "\"name\": " << this->getTag() << ", ";
      s << "\"type\": \"SixNodeTri\", ";
      s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
      s << connectedExternalNodes(1) << ", ";
      s << connectedExternalNodes(2) << ", ";
      s << connectedExternalNodes(3) << ", ";
      s << connectedExternalNodes(4) << ", ";
      s << connectedExternalNodes(5) << "], ";
      s << "\"thickness\": " << thickness << ", ";
      s << "\"surfacePressure\": " << pressure << ", ";
      s << "\"masspervolume\": " << rho << ", ";
      s << "\"bodyForces\": [" << b[0] << ", " << b[1] << "], ";
      s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
  }
}

int
SixNodeTri::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{

    // first set the quantity to be displayed at the nodes;
    // if displayMode is 1 through 3 we will plot material stresses otherwise 0.0

    static Vector values(nip);

    for (int j=0; j<nip; j++)
	   values(j) = 0.0;

    if (displayMode < nip && displayMode > 0) {
	for (int i=0; i<nip; i++) {
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
    const Vector &end5Crd = theNodes[4]->getCrds();
    const Vector &end6Crd = theNodes[5]->getCrds();

    static Matrix coords(nnodes,3);

    if (displayMode >= 0) {

      const Vector &end1Disp = theNodes[0]->getDisp();
      const Vector &end2Disp = theNodes[1]->getDisp();
      const Vector &end3Disp = theNodes[2]->getDisp();
      const Vector &end4Disp = theNodes[3]->getDisp();
      const Vector &end5Disp = theNodes[4]->getDisp();
      const Vector &end6Disp = theNodes[5]->getDisp();

      for (int i = 0; i < 2; i++) {
	coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
	coords(1,i) = end2Crd(i) + end2Disp(i)*fact;
	coords(2,i) = end3Crd(i) + end3Disp(i)*fact;
	coords(3,i) = end4Crd(i) + end4Disp(i)*fact;
	coords(4,i) = end5Crd(i) + end5Disp(i)*fact;
	coords(5,i) = end6Crd(i) + end6Disp(i)*fact;
      }
    } else {
      int mode = displayMode * -1;
      const Matrix &eigen1 = theNodes[0]->getEigenvectors();
      const Matrix &eigen2 = theNodes[1]->getEigenvectors();
      const Matrix &eigen3 = theNodes[2]->getEigenvectors();
      const Matrix &eigen4 = theNodes[3]->getEigenvectors();
      const Matrix &eigen5 = theNodes[4]->getEigenvectors();
      const Matrix &eigen6 = theNodes[5]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < 2; i++) {
	  coords(0,i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	  coords(1,i) = end2Crd(i) + eigen2(i,mode-1)*fact;
	  coords(2,i) = end3Crd(i) + eigen3(i,mode-1)*fact;
	  coords(3,i) = end4Crd(i) + eigen4(i,mode-1)*fact;
	  coords(4,i) = end5Crd(i) + eigen5(i,mode-1)*fact;
	  coords(5,i) = end6Crd(i) + eigen6(i,mode-1)*fact;
	}
      } else {
	for (int i = 0; i < 2; i++) {
	  coords(0,i) = end1Crd(i);
	  coords(1,i) = end2Crd(i);
	  coords(2,i) = end3Crd(i);
	  coords(3,i) = end4Crd(i);
	  coords(4,i) = end5Crd(i);
	  coords(5,i) = end6Crd(i);
	}
      }
    }

    int error = 0;

    // finally we  the element using drawPolygon
    error += theViewer.drawPolygon (coords, values, this->getTag());

    return error;
}

Response*
SixNodeTri::setResponse(const char **argv, int argc,
			  OPS_Stream &output)
{
  Response *theResponse =0;

  output.tag("ElementOutput");
  output.attr("eleType","SixNodeTri");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);
  output.attr("node5",connectedExternalNodes[4]);
  output.attr("node6",connectedExternalNodes[5]);

  char dataOut[20];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=nip; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }

    theResponse =  new ElementResponse(this, 1, P);
  }

  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= nip) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",pts[pointNum-1][0]);
      output.attr("neta",pts[pointNum-1][1]);

      theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, output);

      output.endTag();

    }
  }
  else if ((strcmp(argv[0],"stresses") ==0) || (strcmp(argv[0],"stress") ==0)) {
    for (int i=0; i<nip; i++) {
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
    theResponse =  new ElementResponse(this, 3, Vector(3*nip));
  }

  else if ((strcmp(argv[0],"stressesAtNodes") ==0) || (strcmp(argv[0],"stressAtNodes") ==0)) {
    for (int i=0; i<nnodes; i++) {
      output.tag("NodalPoint");
      output.attr("number",i+1);
      // output.attr("eta",pts[i][0]);
      // output.attr("neta",pts[i][1]);

      // output.tag("NdMaterialOutput");
      // output.attr("classType", theMaterial[i]->getClassTag());
      // output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");

      output.endTag(); // GaussPoint
      // output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 11, Vector(3*nnodes));
  }

  else if ((strcmp(argv[0],"strain") ==0) || (strcmp(argv[0],"strains") ==0)) {
    for (int i=0; i<nip; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType","eta11");
      output.tag("ResponseType","eta22");
      output.tag("ResponseType","eta12");

      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 4, Vector(3*nip));
  }

  output.endTag(); // ElementOutput

  return theResponse;
}

int
SixNodeTri::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(3*nip);
    int cnt = 0;
    for (int i = 0; i < nip; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }

    return eleInfo.setVector(stresses);

  } else if (responseID == 11) {

    // extrapolate stress from Gauss points to element nodes
    static Vector stressGP(3*nip);
    static Vector stressAtNodes(3*nnodes);
	stressAtNodes.Zero();
    int cnt = 0;
	// first get stress components (xx, yy, xy) at Gauss points
    for (int i = 0; i < nip; i++) {
      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stressGP(cnt) = sigma(0);
      stressGP(cnt+1) = sigma(1);
      stressGP(cnt+2) = sigma(2);
      cnt += 3;
    }

    // FMK - cannot be nnodes, nip for certain compilers, error in compilation
    double We[6][3] = {{1.6666666666666667, -0.3333333333333333, -0.3333333333333333},
			 {-0.3333333333333333, 1.6666666666666667, -0.3333333333333333},
			 {-0.3333333333333333, -0.3333333333333333, 1.6666666666666667},
			 {0.6666666666666667, 0.6666666666666667, -0.3333333333333333},
			 {-0.3333333333333333, 0.6666666666666667, 0.6666666666666667},
			 {0.6666666666666667, -0.3333333333333333, 0.6666666666666667}};
    int p, l;
    for (int i = 0; i < nnodes; i++) {
	  for (int k = 0; k < 3; k++) {
		p = 3*i + k;
		for (int j = 0; j < nip; j++) {
		  l = 3*j + k;
		  stressAtNodes(p) += We[i][j] * stressGP(l);
		  // opserr << "stressAtNodes(" << p << ") = We[" << i << "][" << j << "] * stressGP(" << l << ") = " << We[i][j] << " * " << stressGP(l) << " = " << stressAtNodes(p) <<  "\n";
		}
	  }
	}

    return eleInfo.setVector(stressAtNodes);

  } else if (responseID == 4) {

    // Loop over the integration points
    static Vector stresses(3*nip);
    int cnt = 0;
    for (int i = 0; i < nip; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStrain();
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
SixNodeTri::setParameter(const char **argv, int argc, Parameter &param)
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
    if (pointNum > 0 && pointNum <= nip)
      return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else
      return -1;
  }

  // otherwise it could be just a forall material parameter
  else {

    int matRes = res;
    for (int i=0; i<nip; i++) {

      matRes =  theMaterial[i]->setParameter(argv, argc, param);

      if (matRes != -1)
	res = matRes;
    }
  }

  return res;
}

int
SixNodeTri::updateParameter(int parameterID, Information &info)
{
	int res = -1;
		int matRes = res;
  switch (parameterID) {
    case -1:
      return -1;

	case 1:

		for (int i = 0; i<nip; i++) {
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

double SixNodeTri::shapeFunction(double s, double t)
{
	const Vector &nd1Crds = theNodes[0]->getCrds();
	const Vector &nd2Crds = theNodes[1]->getCrds();
	const Vector &nd3Crds = theNodes[2]->getCrds();
	const Vector &nd4Crds = theNodes[3]->getCrds();
	const Vector &nd5Crds = theNodes[4]->getCrds();
	const Vector &nd6Crds = theNodes[5]->getCrds();

	shp[2][0] = s*(2*s-1);
	shp[2][1] = t*(2*t-1);
	shp[2][2] = (1-s-t)*(1-2*s-2*t);
	// shp[2][2] = 1 - 3*s - 3*t + 2*s*s + 2*t*t + 4*s*t;
	shp[2][3] = 4*s*t;
	shp[2][4] = 4*t*(1-s-t);
	shp[2][5] = 4*s*(1-s-t);

	// derivatives
    double N11 = 4*s-1;
	double N12 = 0;
	double N21 = 0;
	double N22 = 4*t-1;
	double N31 = -3+4*s+4*t;
	double N32 = -3+4*t+4*s;
	double N41 = 4*t;
	double N42 = 4*s;
	double N51 = -4*t;
	double N52 = 4-4*s-8*t;
	double N61 = 4-4*t-8*s;
	double N62 = -4*s;

	double J[2][2];

	J[0][0] = nd1Crds(0)*N11 + nd2Crds(0)*N21 + nd3Crds(0)*N31 + nd4Crds(0)*N41 +
	  nd5Crds(0)*N51 + nd6Crds(0)*N61;
	J[0][1] = nd1Crds(0)*N12 + nd2Crds(0)*N22 + nd3Crds(0)*N32 + nd4Crds(0)*N42 +
	  nd5Crds(0)*N52 + nd6Crds(0)*N62;
	J[1][0] = nd1Crds(1)*N11 + nd2Crds(1)*N21 + nd3Crds(1)*N31 + nd4Crds(1)*N41 +
	  nd5Crds(1)*N51 + nd6Crds(1)*N61;
	J[1][1] = nd1Crds(1)*N12 + nd2Crds(1)*N22 + nd3Crds(1)*N32 + nd4Crds(1)*N42 +
	  nd5Crds(1)*N52 + nd6Crds(1)*N62;

	double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

	double oneOverdetJ = 1/detJ;
	double L[2][2];

	// L = inv(J)
	L[0][0] =  J[1][1]*oneOverdetJ;
	L[1][0] = -J[0][1]*oneOverdetJ;
	L[0][1] = -J[1][0]*oneOverdetJ;
	L[1][1] =  J[0][0]*oneOverdetJ;

    double L00 = L[0][0];
    double L10 = L[1][0];
    double L01 = L[0][1];
    double L11 = L[1][1];

    shp[0][0] = L00*N11 + L01*N12;
    shp[0][1] = L00*N21 + L01*N22;
    shp[0][2] = L00*N31 + L01*N32;
    shp[0][3] = L00*N41 + L01*N42;
    shp[0][4] = L00*N51 + L01*N52;
    shp[0][5] = L00*N61 + L01*N62;

    shp[1][0] = L10*N11 + L11*N12;
    shp[1][1] = L10*N21 + L11*N22;
    shp[1][2] = L10*N31 + L11*N32;
    shp[1][3] = L10*N41 + L11*N42;
    shp[1][4] = L10*N51 + L11*N52;
    shp[1][5] = L10*N61 + L11*N62;

    return detJ;
}

void
SixNodeTri::setPressureLoadAtNodes(void)
{
        pressureLoad.Zero();

	if (pressure == 0.0)
		return;

	const Vector &node1 = theNodes[0]->getCrds();
	const Vector &node2 = theNodes[1]->getCrds();
	const Vector &node3 = theNodes[2]->getCrds();
	const Vector &node4 = theNodes[3]->getCrds();
	const Vector &node5 = theNodes[4]->getCrds();
	const Vector &node6 = theNodes[5]->getCrds();

	double x1 = node1(0);
	double y1 = node1(1);
	double x2 = node2(0);
	double y2 = node2(1);
	double x3 = node3(0);
	double y3 = node3(1);
	double x4 = node4(0);
	double y4 = node4(1);
	double x5 = node5(0);
	double y5 = node5(1);
	double x6 = node6(0);
	double y6 = node6(1);

	double dx14 = x4-x1;
	double dy14 = y4-y1;
	double dx42 = x2-x4;
	double dy42 = y2-y4;
	double dx25 = x5-x2;
	double dy25 = y5-y2;
	double dx53 = x3-x5;
	double dy53 = y3-y5;
	double dx36 = x6-x3;
	double dy36 = y6-y3;
	double dx61 = x4-x6;
	double dy61 = y4-y6;

	double fac1 = 0.3333333333333333;
	double fac2 = 0.6666666666666667;

	// Contribution from side 14
	pressureLoad(0) += pressure*fac1*dy14;
	pressureLoad(6) += pressure*fac2*dy14;
	pressureLoad(1) += pressure*fac1*-dx14;
	pressureLoad(7) += pressure*fac2*-dx14;

	// Contribution from side 42
	pressureLoad(6) += pressure*fac2*dy42;
	pressureLoad(2) += pressure*fac1*dy42;
	pressureLoad(7) += pressure*fac2*-dx42;
	pressureLoad(3) += pressure*fac1*-dx42;

	// Contribution from side 25
	pressureLoad(2) += pressure*fac1*dy25;
	pressureLoad(8) += pressure*fac2*dy25;
	pressureLoad(3) += pressure*fac1*-dx25;
	pressureLoad(9) += pressure*fac2*-dx25;

	// Contribution from side 53
	pressureLoad(8) += pressure*fac2*dy53;
	pressureLoad(4) += pressure*fac1*dy53;
	pressureLoad(9) += pressure*fac2*-dx53;
	pressureLoad(5) += pressure*fac1*-dx53;

	// Contribution from side 36
	pressureLoad(4) += pressure*fac1*dy36;
	pressureLoad(10) += pressure*fac2*dy36;
	pressureLoad(5) += pressure*fac1*-dx36;
	pressureLoad(11) += pressure*fac2*-dx36;

	// Contribution from side 61
	pressureLoad(10) += pressure*fac2*dy61;
	pressureLoad(0) += pressure*fac1*dy61;
	pressureLoad(11) += pressure*fac2*-dx61;
	pressureLoad(1) += pressure*fac1*-dx61;

	//pressureLoad = pressureLoad*thickness;
}
