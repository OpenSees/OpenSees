///////////////////////////////////////////////////////////////////////////////
// Description: This file contains the class definition for                  //
// BBarFourNodeQuadUP, a 4-node plane strain element for solid-fluid fully   //
// coupled analysis. This implementation is a simplified u-p formulation     //
// of Biot theory (u - solid displacement, p - fluid pressure).              //
// Each element node has two DOFs for u and 1 DOF for p.                     //
// Constant volume/pressure integration (BBar method) is used for integration//
// of the volumetric component of solid phase and the fulid phase.           //
//                                                                           //
// Written by Zhaohui Yang	(June 2009)                                      //
///////////////////////////////////////////////////////////////////////////////

// $Revision: 1.1 $
// $Date: 2009-10-07 20:02:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/UP-ucsd/BBarFourNodeQuadUP.cpp,v $

#include <BBarFourNodeQuadUP.h>
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

void* OPS_BBarFourNodeQuadUP()
{
    if (OPS_GetNDM() != 2 || OPS_GetNDF() != 3) {
	opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with QuadUP element\n";
	return 0;
    }

    if (OPS_GetNumRemainingInputArgs() < 11) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: element bbarQuadUP eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? bulk? rho? perm_x? perm_y? <b1? b2? pressure? dM? dK?>\n";
	return 0;
    }

    // BBarFourNodeQuadUPId, iNode, jNode, kNode, lNode
    int tags[5];
    int num = 5;
    if (OPS_GetIntInput(&num,tags) < 0) {
	opserr<<"WARNING: invalid integer input\n";
	return 0;
    }

    double thk;
    num = 1;
    if (OPS_GetDoubleInput(&num,&thk) < 0) {
	opserr<<"WARNING: invalid double input\n";
	return 0;
    }

    int matTag;
    if (OPS_GetIntInput(&num,&matTag) < 0) {
	opserr<<"WARNING: invalid integer input\n";
	return 0;
    }
    NDMaterial* mat = OPS_getNDMaterial(matTag);
    if (mat == 0) {
	opserr << "WARNING material not found\n";
	opserr << "Material: " << matTag;
	opserr << "\nBBarFourNodeQuadUP element: " << tags[0] << endln;
	return 0;
    }

    // bk, r, perm1, perm2
    double data[4];
    num = 4;
    if (OPS_GetDoubleInput(&num,data) < 0) {
	opserr<<"WARNING: invalid double input\n";
	return 0;
    }

    // b1, b2, p
    double opt[3] = {0,0,0};
    num = OPS_GetNumRemainingInputArgs();
    if (num > 3) {
	num = 3;
    }
    if (num > 0) {
	if (OPS_GetDoubleInput(&num,opt) < 0) {
	    opserr<<"WARNING: invalid double input\n";
	    return 0;
	}
    }

    return new BBarFourNodeQuadUP(tags[0],tags[1],tags[2],tags[3],tags[4],
				  *mat,"PlaneStrain",thk,data[0],data[1],data[2],data[3],
				  opt[0],opt[1],opt[2]);
}


Matrix BBarFourNodeQuadUP::K(12,12);
Vector BBarFourNodeQuadUP::P(12);
double BBarFourNodeQuadUP::shp[3][4][4];
double BBarFourNodeQuadUP::pts[4][2];
double BBarFourNodeQuadUP::wts[4];
double BBarFourNodeQuadUP::dvol[4];
double BBarFourNodeQuadUP::shpBar[2][4];
double BBarFourNodeQuadUP::B[4][2][4][4];
double BBarFourNodeQuadUP::Bp[2][4][4];
Node *BBarFourNodeQuadUP::theNodes[4];


BBarFourNodeQuadUP::BBarFourNodeQuadUP(int tag, int nd1, int nd2, int nd3, int nd4,
	NDMaterial &m, const char *type, double t, double bulk, double r,
		  double p1, double p2, double b1, double b2, double p)
:Element (tag, ELE_TAG_BBarFourNodeQuadUP),
  theMaterial(0), connectedExternalNodes(4),
  nd1Ptr(0), nd2Ptr(0), nd3Ptr(0), nd4Ptr(0), Ki(0),
  Q(12), pressureLoad(12), applyLoad(0), thickness(t), kc(bulk), rho(r), pressure(p)
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

	// Body forces
	b[0] = b1;
	b[1] = b2;
	// Permeabilities
  perm[0] = p1;
  perm[1] = p2;

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new NDMaterial *[4];

    if (theMaterial == 0) {
      opserr << "BBarFourNodeQuadUP::BBarFourNodeQuadUP - failed allocate material model pointer\n";
      exit(-1);
    }

    for (int i = 0; i < 4; i++) {

      // Get copies of the material model for each integration point
      theMaterial[i] = m.getCopy(type);

      // Check allocation
      if (theMaterial[i] == 0) {
	    opserr << "BBarFourNodeQuadUP::BBarFourNodeQuadUP -- failed to get a copy of material model\n";
	    exit(-1);
      }

	  // Need 3D stresses
	  Information info;
	  theMaterial[i]->updateParameter(20, info);
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
}

BBarFourNodeQuadUP::BBarFourNodeQuadUP()
:Element (0,ELE_TAG_BBarFourNodeQuadUP),
  theMaterial(0), connectedExternalNodes(4),
 nd1Ptr(0), nd2Ptr(0), nd3Ptr(0), nd4Ptr(0), Ki(0),
  Q(12), pressureLoad(12), applyLoad(0), thickness(0.0), kc(0.0), rho(0.0), pressure(0.0)
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
}

BBarFourNodeQuadUP::~BBarFourNodeQuadUP()
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
BBarFourNodeQuadUP::getNumExternalNodes() const
{
    return 4;
}

const ID&
BBarFourNodeQuadUP::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
BBarFourNodeQuadUP::getNodePtrs()
{
  theNodes[0] = nd1Ptr;
  theNodes[1] = nd2Ptr;
  theNodes[2] = nd3Ptr;
  theNodes[3] = nd4Ptr;

  return theNodes;
}

int
BBarFourNodeQuadUP::getNumDOF()
{
    return 12;
}

void
BBarFourNodeQuadUP::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	nd1Ptr = 0;
	nd2Ptr = 0;
	nd3Ptr = 0;
	nd4Ptr = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);

    nd1Ptr = theDomain->getNode(Nd1);
    nd2Ptr = theDomain->getNode(Nd2);
    nd3Ptr = theDomain->getNode(Nd3);
    nd4Ptr = theDomain->getNode(Nd4);

    if (nd1Ptr == 0 || nd2Ptr == 0 || nd3Ptr == 0 || nd4Ptr == 0) {
	//opserr << "FATAL ERROR BBarFourNodeQuadUP (tag: %d), node not found in domain",
	//	this->getTag());

	return;
    }

    int dofNd1 = nd1Ptr->getNumberDOF();
    int dofNd2 = nd2Ptr->getNumberDOF();
    int dofNd3 = nd3Ptr->getNumberDOF();
    int dofNd4 = nd4Ptr->getNumberDOF();

    if (dofNd1 != 3 || dofNd2 != 3 || dofNd3 != 3 || dofNd4 != 3) {
	//opserr << "FATAL ERROR BBarFourNodeQuadUP (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());

	return;
    }
    this->DomainComponent::setDomain(theDomain);

	// Compute consistent nodal loads due to pressure
	this->setPressureLoadAtNodes();
}

int
BBarFourNodeQuadUP::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "FourNodeQuad_UP::commitState () - failed in base class";
    }

    // Loop over the integration points and commit the material states
    for (int i = 0; i < 4; i++)
      retVal += theMaterial[i]->commitState();

    return retVal;
}

int
BBarFourNodeQuadUP::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}

int
BBarFourNodeQuadUP::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}

int
BBarFourNodeQuadUP::update()
{
	const Vector &disp1 = nd1Ptr->getTrialDisp();
	const Vector &disp2 = nd2Ptr->getTrialDisp();
	const Vector &disp3 = nd3Ptr->getTrialDisp();
	const Vector &disp4 = nd4Ptr->getTrialDisp();

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

	// Determine Jacobian for this integration point
	this->shapeFunction();

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

		// Interpolate strains
		//eps = B*u;
		//eps.addMatrixVector(0.0, B, u, 1.0);
		eps.Zero();
		for (int beta = 0; beta < 4; beta++) {
			//eps(0) += shp[0][beta][i]*u[0][beta];
			//eps(1) += shp[1][beta][i]*u[1][beta];
			//eps(2) += shp[0][beta][i]*u[1][beta] + shp[1][beta][i]*u[0][beta];
			eps(0) += B[0][0][beta][i]*u[0][beta]+B[0][1][beta][i]*u[1][beta];
			eps(1) += B[1][0][beta][i]*u[0][beta]+B[1][1][beta][i]*u[1][beta];
			eps(2) += B[2][0][beta][i]*u[0][beta]+B[2][1][beta][i]*u[1][beta];
		}

		// Set the material strain
		ret += theMaterial[i]->setTrialStrain(eps);
	}

	return ret;
}


const Matrix&
BBarFourNodeQuadUP::getTangentStiff()
{

  K.Zero();

  double DB[4][2];

  // Determine Jacobian for this integration point
  this->shapeFunction();

  // Loop over the integration points
  for (int i = 0; i < 4; i++) {

    // Get the material tangent
    const Matrix &D = theMaterial[i]->getTangent();

    // Perform numerical integration
    //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
    //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
    for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 3) {

      for (int beta = 0, ib = 0; beta < 4; beta++, ib += 3) {

	    //DB[0][0] = dvol[i] * (D(0,0)*shp[0][beta][i] + D(0,2)*shp[1][beta][i]);
	    //DB[1][0] = dvol[i] * (D(1,0)*shp[0][beta][i] + D(1,2)*shp[1][beta][i]);
	    //DB[2][0] = dvol[i] * (D(2,0)*shp[0][beta][i] + D(2,2)*shp[1][beta][i]);
	    //DB[0][1] = dvol[i] * (D(0,1)*shp[1][beta][i] + D(0,2)*shp[0][beta][i]);
	    //DB[1][1] = dvol[i] * (D(1,1)*shp[1][beta][i] + D(1,2)*shp[0][beta][i]);
	    //DB[2][1] = dvol[i] * (D(2,1)*shp[1][beta][i] + D(2,2)*shp[0][beta][i]);

	    //K(ia,ib) += shp[0][alpha][i]*DB[0][0] + shp[1][alpha][i]*DB[2][0];
	    //K(ia,ib+1) += shp[0][alpha][i]*DB[0][1] + shp[1][alpha][i]*DB[2][1];
	    //K(ia+1,ib) += shp[1][alpha][i]*DB[1][0] + shp[0][alpha][i]*DB[2][0];
	    //K(ia+1,ib+1) += shp[1][alpha][i]*DB[1][1] + shp[0][alpha][i]*DB[2][1];
	    DB[0][0] = dvol[i] * (D(0,0)*B[0][0][beta][i] + D(0,1)*B[1][0][beta][i]
	                         +D(0,3)*B[2][0][beta][i] + D(0,2)*B[3][0][beta][i]);
	    DB[1][0] = dvol[i] * (D(1,0)*B[0][0][beta][i] + D(1,1)*B[1][0][beta][i]
	                         +D(1,3)*B[2][0][beta][i] + D(1,2)*B[3][0][beta][i]);
	    DB[2][0] = dvol[i] * (D(2,0)*B[0][0][beta][i] + D(2,1)*B[1][0][beta][i]
	                         +D(2,3)*B[2][0][beta][i] + D(2,2)*B[3][0][beta][i]);
	    DB[3][0] = dvol[i] * (D(3,0)*B[0][0][beta][i] + D(3,1)*B[1][0][beta][i]
	                         +D(3,3)*B[2][0][beta][i] + D(3,2)*B[3][0][beta][i]);
	    DB[0][1] = dvol[i] * (D(0,0)*B[0][1][beta][i] + D(0,1)*B[1][1][beta][i]
	                         +D(0,3)*B[2][1][beta][i] + D(0,2)*B[3][1][beta][i]);
	    DB[1][1] = dvol[i] * (D(1,0)*B[0][1][beta][i] + D(1,1)*B[1][1][beta][i]
	                         +D(1,3)*B[2][1][beta][i] + D(1,2)*B[3][1][beta][i]);
	    DB[2][1] = dvol[i] * (D(2,0)*B[0][1][beta][i] + D(2,1)*B[1][1][beta][i]
	                         +D(2,3)*B[2][1][beta][i] + D(2,2)*B[3][1][beta][i]);
	    DB[3][1] = dvol[i] * (D(3,0)*B[0][1][beta][i] + D(3,1)*B[1][1][beta][i]
	                         +D(3,3)*B[2][1][beta][i] + D(3,2)*B[3][1][beta][i]);

	    K(ia,ib) += B[0][0][alpha][i]*DB[0][0] + B[1][0][alpha][i]*DB[1][0]
	               +B[3][0][alpha][i]*DB[2][0] + B[2][0][alpha][i]*DB[3][0];
	    K(ia,ib+1) += B[0][0][alpha][i]*DB[0][1] + B[1][0][alpha][i]*DB[1][1]
	                 +B[3][0][alpha][i]*DB[2][1] + B[2][0][alpha][i]*DB[3][1];
	    K(ia+1,ib) += B[0][1][alpha][i]*DB[0][0] + B[1][1][alpha][i]*DB[1][0]
	                 +B[3][1][alpha][i]*DB[2][0] + B[2][1][alpha][i]*DB[3][0];
	    K(ia+1,ib+1) += B[0][1][alpha][i]*DB[0][1] + B[1][1][alpha][i]*DB[1][1]
	                   +B[3][1][alpha][i]*DB[2][1] + B[2][1][alpha][i]*DB[3][1];

      }
    }
  }
  return K;
}


const Matrix &BBarFourNodeQuadUP::getInitialStiff ()
{
  if (Ki != 0) return *Ki;

  K.Zero();

  double DB[4][2];

  // Determine Jacobian for this integration point
  this->shapeFunction();

  // Loop over the integration points
  for (int i = 0; i < 4; i++) {

    // Get the material tangent
    const Matrix &D = theMaterial[i]->getInitialTangent();

    // Perform numerical integration
    //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
    //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
    for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 3) {

      for (int beta = 0, ib = 0; beta < 4; beta++, ib += 3) {

	    //DB[0][0] = dvol[i] * (D(0,0)*shp[0][beta][i] + D(0,2)*shp[1][beta][i]);
	    //DB[1][0] = dvol[i] * (D(1,0)*shp[0][beta][i] + D(1,2)*shp[1][beta][i]);
	    //DB[2][0] = dvol[i] * (D(2,0)*shp[0][beta][i] + D(2,2)*shp[1][beta][i]);
	    //DB[0][1] = dvol[i] * (D(0,1)*shp[1][beta][i] + D(0,2)*shp[0][beta][i]);
	    //DB[1][1] = dvol[i] * (D(1,1)*shp[1][beta][i] + D(1,2)*shp[0][beta][i]);
	    //DB[2][1] = dvol[i] * (D(2,1)*shp[1][beta][i] + D(2,2)*shp[0][beta][i]);

	    //K(ia,ib) += shp[0][alpha][i]*DB[0][0] + shp[1][alpha][i]*DB[2][0];
	    //K(ia,ib+1) += shp[0][alpha][i]*DB[0][1] + shp[1][alpha][i]*DB[2][1];
	    //K(ia+1,ib) += shp[1][alpha][i]*DB[1][0] + shp[0][alpha][i]*DB[2][0];
	    //K(ia+1,ib+1) += shp[1][alpha][i]*DB[1][1] + shp[0][alpha][i]*DB[2][1];
	    DB[0][0] = dvol[i] * (D(0,0)*B[0][0][beta][i] + D(0,1)*B[1][0][beta][i]
	                         +D(0,3)*B[2][0][beta][i] + D(0,2)*B[3][0][beta][i]);
	    DB[1][0] = dvol[i] * (D(1,0)*B[0][0][beta][i] + D(1,1)*B[1][0][beta][i]
	                         +D(1,3)*B[2][0][beta][i] + D(1,2)*B[3][0][beta][i]);
	    DB[2][0] = dvol[i] * (D(2,0)*B[0][0][beta][i] + D(2,1)*B[1][0][beta][i]
	                         +D(2,3)*B[2][0][beta][i] + D(2,2)*B[3][0][beta][i]);
	    DB[3][0] = dvol[i] * (D(3,0)*B[0][0][beta][i] + D(3,1)*B[1][0][beta][i]
	                         +D(3,3)*B[2][0][beta][i] + D(3,2)*B[3][0][beta][i]);
	    DB[0][1] = dvol[i] * (D(0,0)*B[0][1][beta][i] + D(0,1)*B[1][1][beta][i]
	                         +D(0,3)*B[2][1][beta][i] + D(0,2)*B[3][1][beta][i]);
	    DB[1][1] = dvol[i] * (D(1,0)*B[0][1][beta][i] + D(1,1)*B[1][1][beta][i]
	                         +D(1,3)*B[2][1][beta][i] + D(1,2)*B[3][1][beta][i]);
	    DB[2][1] = dvol[i] * (D(2,0)*B[0][1][beta][i] + D(2,1)*B[1][1][beta][i]
	                         +D(2,3)*B[2][1][beta][i] + D(2,2)*B[3][1][beta][i]);
	    DB[3][1] = dvol[i] * (D(3,0)*B[0][1][beta][i] + D(3,1)*B[1][1][beta][i]
	                         +D(3,3)*B[2][1][beta][i] + D(3,2)*B[3][1][beta][i]);

	    K(ia,ib) += B[0][0][alpha][i]*DB[0][0] + B[1][0][alpha][i]*DB[1][0]
	               +B[3][0][alpha][i]*DB[2][0] + B[2][0][alpha][i]*DB[3][0];
	    K(ia,ib+1) += B[0][0][alpha][i]*DB[0][1] + B[1][0][alpha][i]*DB[1][1]
	                 +B[3][0][alpha][i]*DB[2][1] + B[2][0][alpha][i]*DB[3][1];
	    K(ia+1,ib) += B[0][1][alpha][i]*DB[0][0] + B[1][1][alpha][i]*DB[1][0]
	                 +B[3][1][alpha][i]*DB[2][0] + B[2][1][alpha][i]*DB[3][0];
	    K(ia+1,ib+1) += B[0][1][alpha][i]*DB[0][1] + B[1][1][alpha][i]*DB[1][1]
	                   +B[3][1][alpha][i]*DB[2][1] + B[2][1][alpha][i]*DB[3][1];
      }
    }
  }

  Ki = new Matrix(K);
  if (Ki == 0) {
    opserr << "FATAL BBarFourNodeQuadUP::getInitialStiff() -";
    opserr << "ran out of memory\n";
    exit(-1);
  }

  return *Ki;
}

const Matrix&
BBarFourNodeQuadUP::getDamp()
{
  static Matrix Kdamp(12,12);
  Kdamp.Zero();

  if (betaK != 0.0)
    Kdamp.addMatrix(1.0, this->getTangentStiff(), betaK);
  if (betaK0 != 0.0)
    Kdamp.addMatrix(1.0, this->getInitialStiff(), betaK0);
  if (betaKc != 0.0)
    Kdamp.addMatrix(1.0, *Kc, betaKc);

  int i, j, m, i1, j1;

  if (alphaM != 0.0) {
	this->getMass();
    for (i = 0; i < 12; i += 3) {
      for (j = 0; j < 12; j += 3) {
        Kdamp(i,j) += K(i,j)*alphaM;
        Kdamp(i+1,j+1) += K(i+1,j+1)*alphaM;
	  }
    }
  }

  // Determine Jacobian for this integration point
  this->shapeFunction();

  // Compute coupling matrix
  double vol = dvol[0] + dvol[1] + dvol[2] + dvol[3];
  for (i = 0; i < 12; i += 3) {
    i1 = i / 3;
    for (j = 2; j < 12; j += 3) {
      j1 = (j-2) / 3;
      for (m = 0; m < 4; m++) {
	    //Kdamp(i,j) += -dvol[m]*shp[0][i1][m]*shp[2][j1][m];
	    //Kdamp(i+1,j) += -dvol[m]*shp[1][i1][m]*shp[2][j1][m];
	    Kdamp(i,j) += -dvol[m]*(B[0][0][i1][m]+B[1][0][i1][m]+B[3][0][i1][m])*shp[2][j1][m];
	    Kdamp(i+1,j) += -dvol[m]*(B[0][1][i1][m]+B[1][1][i1][m]+B[3][1][i1][m])*shp[2][j1][m];
	  }
      Kdamp(j,i) = Kdamp(i,j);
      Kdamp(j,i+1) = Kdamp(i+1,j);
    }
  }

  // Compute permeability matrix
  for (i = 2; i < 12; i += 3) {
    int i1 = (i-2) / 3;
    for (j = 2; j < 12; j += 3) {
      int j1 = (j-2) / 3;
      for (m = 0; m < 4; m++) {
	    //Kdamp(i,j) += - dvol[m]*(perm[0]*shp[0][i1][m]*shp[0][j1][m] +
	    //              perm[1]*shp[1][i1][m]*shp[1][j1][m]);
	    Kdamp(i,j) += - dvol[m]*(perm[0]*Bp[0][i1][m]*Bp[0][j1][m] +
	                  perm[1]*Bp[1][i1][m]*Bp[1][j1][m]);
	  }
    }
  }

  K = Kdamp;
  return K;
}

const Matrix&
BBarFourNodeQuadUP::getMass()
{
  K.Zero();

  int i, j, m, i1, j1;
  double Nrho;

  // Determine Jacobian for this integration point
  this->shapeFunction();


  // Compute an ad hoc lumped mass matrix
  /*for (i = 0; i < 4; i++) {

    // average material density
    tmp = mixtureRho(i);

    for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 3) {
      Nrho = shp[2][alpha][i]*dvol[i]*tmp;
      K(ia,ia) += Nrho;
      K(ia+1,ia+1) += Nrho;
    }
  }*/

  // Compute consistent mass matrix
  for (i = 0, i1 = 0; i < 12; i += 3, i1++) {
    for (j = 0, j1 = 0; j < 12; j += 3, j1++) {
    for (m = 0; m < 4; m++) {
    Nrho = dvol[m]*mixtureRho(m)*shp[2][i1][m]*shp[2][j1][m];
    K(i,j) += Nrho;
    K(i+1,j+1) += Nrho;
    }
    }
  }

  // Compute compressibility matrix
  double vol = dvol[0] + dvol[1] + dvol[2] + dvol[3];
  double oneOverKc = 1./kc;

  for (i = 2; i < 12; i += 3) {
    i1 = (i-2) / 3;
    for (j = 2; j < 12; j += 3) {
      j1 = (j-2) / 3;
      for (m = 0; m < 4; m++) {
	    K(i,j) += -dvol[m]*oneOverKc*shp[2][i1][m]*shp[2][j1][m];
	  }
    }
  }

  return K;
}

void
BBarFourNodeQuadUP::zeroLoad(void)
{
	Q.Zero();

	applyLoad = 0;
	appliedB[0] = 0.0;
	appliedB[1] = 0.0;

	return;
}

int
BBarFourNodeQuadUP::addLoad(ElementalLoad *theLoad, double loadFactor)
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
		opserr << "BBarFourNodeQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	} 

	return -1;
}


int
BBarFourNodeQuadUP::addInertiaLoadToUnbalance(const Vector &accel)
{
  // accel = uDotDotG (see EarthquakePattern.cpp)
  // Get R * accel from the nodes
  const Vector &Raccel1 = nd1Ptr->getRV(accel);
  const Vector &Raccel2 = nd2Ptr->getRV(accel);
  const Vector &Raccel3 = nd3Ptr->getRV(accel);
  const Vector &Raccel4 = nd4Ptr->getRV(accel);

  if (3 != Raccel1.Size() || 3 != Raccel2.Size() || 3 != Raccel3.Size() ||
      3 != Raccel4.Size()) {
    opserr << "BBarFourNodeQuadUP::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }

  double ra[12];

  ra[0] = Raccel1(0);
  ra[1] = Raccel1(1);
  ra[2] = 0.;
  ra[3] = Raccel2(0);
  ra[4] = Raccel2(1);
  ra[5] = 0.;
  ra[6] = Raccel3(0);
  ra[7] = Raccel3(1);
  ra[8] = 0.;
  ra[9] = Raccel4(0);
  ra[10] = Raccel4(1);
  ra[11] = 0.;

  // Compute mass matrix
  this->getMass();

  // Want to add ( - fact * M R * accel ) to unbalance
  int i, j;

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++)
      Q(i) += -K(i,j)*ra[j];
  }

  return 0;
}

const Vector&
BBarFourNodeQuadUP::getResistingForce()
{
  P.Zero();

  // Determine Jacobian for this integration point
  this->shapeFunction();

  int i;
  // Loop over the integration points
  for (i = 0; i < 4; i++) {

    // Get material stress response
    const Vector &sigma = theMaterial[i]->getStress();

    // Perform numerical integration on internal force
    //P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
    //P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);
    for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 3) {

      //P(ia) += dvol[i]*(shp[0][alpha][i]*sigma(0) + shp[1][alpha][i]*sigma(2));
      //P(ia+1) += dvol[i]*(shp[1][alpha][i]*sigma(1) + shp[0][alpha][i]*sigma(2));
      P(ia) += dvol[i]*(B[0][0][alpha][i]*sigma(0) + B[1][0][alpha][i]*sigma(1)
                       +B[2][0][alpha][i]*sigma(3) + B[3][0][alpha][i]*sigma(2));
      P(ia+1) += dvol[i]*(B[0][1][alpha][i]*sigma(0) + B[1][1][alpha][i]*sigma(1)
                         +B[2][1][alpha][i]*sigma(3) + B[3][1][alpha][i]*sigma(2));

      // Subtract equiv. body forces from the nodes
      //P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
      //P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);

      double r = mixtureRho(i);
	  if (applyLoad == 0) {
      	P(ia) -= dvol[i]*(shp[2][alpha][i]*r*b[0]);
      	P(ia+1) -= dvol[i]*(shp[2][alpha][i]*r*b[1]);
	  } else {
		P(ia) -= dvol[i]*(shp[2][alpha][i]*r*appliedB[0]);
      	P(ia+1) -= dvol[i]*(shp[2][alpha][i]*r*appliedB[1]);
	  }
    }
  }

  // Subtract fluid body force
  for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 3) {
    for (i = 0; i < 4; i++) {
      //P(ia+2) += dvol[i]*rho*(perm[0]*b[0]*shp[0][alpha][i] +
      //           perm[1]*b[1]*shp[1][alpha][i]);
	  if (applyLoad == 0) {
      	P(ia+2) += dvol[i]*rho*(perm[0]*b[0]*Bp[0][alpha][i] +
                   perm[1]*b[1]*Bp[1][alpha][i]);
	  } else {
		P(ia+2) += dvol[i]*rho*(perm[0]*appliedB[0]*Bp[0][alpha][i] +
                   perm[1]*appliedB[1]*Bp[1][alpha][i]);
	  }
    }
  }

  // Subtract pressure loading from resisting force
  if (pressure != 0.0) {
    //P = P + pressureLoad;
    P.addVector(1.0, pressureLoad, -1.0);
  }

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  //P = P - Q;
  P.addVector(1.0, Q, -1.0);

  return P;
}

const Vector&
BBarFourNodeQuadUP::getResistingForceIncInertia()
{
  int i, j, k;

  const Vector &accel1 = nd1Ptr->getTrialAccel();
  const Vector &accel2 = nd2Ptr->getTrialAccel();
  const Vector &accel3 = nd3Ptr->getTrialAccel();
  const Vector &accel4 = nd4Ptr->getTrialAccel();

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
  //opserr<<"K "<<P<<endln;

  // Compute the mass matrix
  this->getMass();

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++)
      P(i) += K(i,j)*a[j];
  }
  //opserr<<"K+M "<<P<<endln;

  // dynamic seepage force
  /*for (i = 0, k = 0; i < 4; i++, k += 3) {
    // loop over integration points
    for (j = 0; j < 4; j++) {
      P(i+2) -= rho*dvol[j]*(shp[2][i][j]*a[k]*perm[0]*shp[0][i][j]
			     +shp[2][i][j]*a[k+1]*perm[1]*shp[1][i][j]);
    }
  }*/
  //opserr<<"K+M+fb "<<P<<endln;

  const Vector &vel1 = nd1Ptr->getTrialVel();
  const Vector &vel2 = nd2Ptr->getTrialVel();
  const Vector &vel3 = nd3Ptr->getTrialVel();
  const Vector &vel4 = nd4Ptr->getTrialVel();

  a[0] = vel1(0);
  a[1] = vel1(1);
  a[2] = vel1(2);
  a[3] = vel2(0);
  a[4] = vel2(1);
  a[5] = vel2(2);
  a[6] = vel3(0);
  a[7] = vel3(1);
  a[8] = vel3(2);
  a[9] = vel4(0);
  a[10] = vel4(1);
  a[11] = vel4(2);

  this->getDamp();

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      P(i) += K(i,j)*a[j];
    }
  }
  //opserr<<"final "<<P<<endln;
  return P;
}

int
BBarFourNodeQuadUP::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // Quad packs its data into a Vector and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
  static Vector data(13);
  data(0) = this->getTag();
  data(1) = thickness;
  data(2) = rho;
  data(3) = b[0];
  data(4) = b[1];
  data(5) = pressure;

  data(6) = alphaM;
  data(7) = betaK;
  data(8) = betaK0;
  data(9) = betaKc;

  data(10) = kc;
  data(11) = perm[0];
  data(12) = perm[1];

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING BBarFourNodeQuadUP::sendSelf() - " << this->getTag() << " failed to send Vector\n";
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
    opserr << "WARNING BBarFourNodeQuadUP::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING BBarFourNodeQuadUP::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}

int
BBarFourNodeQuadUP::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the
  // internal data with the data in the Vector
  static Vector data(13);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING BBarFourNodeQuadUP::recvSelf() - failed to receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  thickness = data(1);
  rho = data(2);
  b[0] = data(3);
  b[1] = data(4);
  pressure = data(5);

  alphaM = data(6);
  betaK = data(7);
  betaK0 = data(8);
  betaKc = data(9);

  kc = data(10);
  perm[0] = data(11);
  perm[1] = data(12);

  static ID idData(12);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING BBarFourNodeQuadUP::recvSelf() - " << this->getTag() << " failed to receive ID\n";
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
      opserr << "BBarFourNodeQuadUP::recvSelf() - Could not allocate NDMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	opserr << "BBarFourNodeQuadUP::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
opserr << "BBarFourNodeQuadUP::recvSelf() - material " << i << "failed to recv itself\n";
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
opserr << "BBarFourNodeQuadUP::recvSelf() - material " << i << "failed to create\n";

	  return -1;
	}
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
opserr << "BBarFourNodeQuadUP::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }

  return res;
}

void
BBarFourNodeQuadUP::Print(OPS_Stream &s, int flag)
{
	s << "\nBBarFourNodeQuadUP, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tthickness:  " << thickness << endln;
	s << "\tmass density:  " << rho << endln;
	s << "\tsurface pressure:  " << pressure << endln;
	s << "\tbody forces:  " << b[0] << ' ' << b[1] << endln;
	theMaterial[0]->Print(s,flag);
	s << "\tStress (xx yy xy)" << endln;
	for (int i = 0; i < 4; i++)
		s << "\t\tGauss point " << i+1 << ": " << theMaterial[i]->getStress();
}

int
BBarFourNodeQuadUP::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
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
    const Vector &end1Crd = nd1Ptr->getCrds();
    const Vector &end2Crd = nd2Ptr->getCrds();
    const Vector &end3Crd = nd3Ptr->getCrds();
    const Vector &end4Crd = nd4Ptr->getCrds();

    const Vector &end1Disp = nd1Ptr->getDisp();
    const Vector &end2Disp = nd2Ptr->getDisp();
    const Vector &end3Disp = nd3Ptr->getDisp();
    const Vector &end4Disp = nd4Ptr->getDisp();

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
BBarFourNodeQuadUP::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","BBarFourNodeQuadUP");
  output.attr("eleTag",this->getTag());
  output.attr("node1",nd1Ptr->getTag());
  output.attr("node2",nd2Ptr->getTag());
  output.attr("node3",nd3Ptr->getTag());
  output.attr("node4",nd4Ptr->getTag());

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=4; i++) {
      sprintf(outputData,"P1_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P2_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"Pp_%d",i);
      output.tag("ResponseType",outputData);
    }

    theResponse = new ElementResponse(this, 1, P);

  }  else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0) {
    return new ElementResponse(this, 2, K);

  } else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);

      theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, output);

      output.endTag(); // GaussPoint
    }
  }

  output.endTag(); // ElementOutput
  return theResponse;
}

int
BBarFourNodeQuadUP::getResponse(int responseID, Information &eleInfo)
{
  switch (responseID) {

  case 1:
    return eleInfo.setVector(this->getResistingForce());

  case 2:
    return eleInfo.setMatrix(this->getTangentStiff());

  default:
    return -1;
  }
}

int
BBarFourNodeQuadUP::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;


  // quad mass density per unit volume
  if (strcmp(argv[0],"rho") == 0) {
    return param.addObject(1, this);

  // quad pressure loading
  } else if (strcmp(argv[0],"pressure") == 0) {
    return param.addObject(2, this);

  // permeability in horizontal direction
  } else if (strcmp(argv[0],"hPerm") == 0) {
    return param.addObject(3, this);

  // permeability in vertical direction
  } else if (strcmp(argv[0],"vPerm") == 0) {
    return param.addObject(4, this);
  }
  // check for material parameters
  if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {

    if (argc < 3)
      return -1;

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4)
      return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else
      return -1;
  }

  // otherwise it could be a for all material pointer
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
BBarFourNodeQuadUP::updateParameter(int parameterID, Information &info)
{
  int res = -1;
  int matRes = res;
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
	case 3:
		perm[0] = info.theDouble;
		this->getDamp();	// update mass matrix
		return 0;
	case 4:
		perm[1] = info.theDouble;
		this->getDamp();	// update mass matrix
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

void BBarFourNodeQuadUP::shapeFunction(void)
{
	double xi, eta, oneMinuseta, onePluseta, oneMinusxi, onePlusxi,
		     detJ, oneOverdetJ, J[2][2], L[2][2], L00, L01, L10, L11,
				 L00oneMinuseta, L00onePluseta, L01oneMinusxi, L01onePlusxi,
				 L10oneMinuseta, L10onePluseta, L11oneMinusxi, L11onePlusxi,
				 vol = 0.0;
    int i, k, l;

	for (k=0; k<2; k++) {
		for (l=0; l<4; l++) {
			shpBar[k][l] = 0.0;
		}
	}

	// loop over integration points
	for (i=0; i<4; i++) {
		xi = pts[i][0];
		eta = pts[i][1];
	  const Vector &nd1Crds = nd1Ptr->getCrds();
	  const Vector &nd2Crds = nd2Ptr->getCrds();
	  const Vector &nd3Crds = nd3Ptr->getCrds();
	  const Vector &nd4Crds = nd4Ptr->getCrds();

	  oneMinuseta = 1.0-eta;
	  onePluseta = 1.0+eta;
	  oneMinusxi = 1.0-xi;
	  onePlusxi = 1.0+xi;

	  shp[2][0][i] = 0.25*oneMinusxi*oneMinuseta;	// N_1
	  shp[2][1][i] = 0.25*onePlusxi*oneMinuseta;		// N_2
	  shp[2][2][i] = 0.25*onePlusxi*onePluseta;		// N_3
	  shp[2][3][i] = 0.25*oneMinusxi*onePluseta;		// N_4

	  J[0][0] = 0.25 * (-nd1Crds(0)*oneMinuseta + nd2Crds(0)*oneMinuseta +
				nd3Crds(0)*(onePluseta) - nd4Crds(0)*(onePluseta));

	  J[0][1] = 0.25 * (-nd1Crds(0)*oneMinusxi - nd2Crds(0)*onePlusxi +
				nd3Crds(0)*onePlusxi + nd4Crds(0)*oneMinusxi);

	  J[1][0] = 0.25 * (-nd1Crds(1)*oneMinuseta + nd2Crds(1)*oneMinuseta +
				nd3Crds(1)*onePluseta - nd4Crds(1)*onePluseta);

	  J[1][1] = 0.25 * (-nd1Crds(1)*oneMinusxi - nd2Crds(1)*onePlusxi +
				nd3Crds(1)*onePlusxi + nd4Crds(1)*oneMinusxi);

	  detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
	  oneOverdetJ = 1.0/detJ;

	  // L = inv(J)
	  L[0][0] =  J[1][1]*oneOverdetJ;
	  L[1][0] = -J[0][1]*oneOverdetJ;
	  L[0][1] = -J[1][0]*oneOverdetJ;
	  L[1][1] =  J[0][0]*oneOverdetJ;

      L00 = 0.25*L[0][0];
      L10 = 0.25*L[1][0];
      L01 = 0.25*L[0][1];
      L11 = 0.25*L[1][1];

	  L00oneMinuseta = L00*oneMinuseta;
	  L00onePluseta  = L00*onePluseta;
	  L01oneMinusxi  = L01*oneMinusxi;
	  L01onePlusxi   = L01*onePlusxi;

	  L10oneMinuseta = L10*oneMinuseta;
	  L10onePluseta  = L10*onePluseta;
	  L11oneMinusxi  = L11*oneMinusxi;
	  L11onePlusxi   = L11*onePlusxi;

	  // B: See Cook, Malkus, Plesha p. 169 for the derivation of these terms
      shp[0][0][i] = -L00oneMinuseta - L01oneMinusxi;	// N_1,1
      shp[0][1][i] =  L00oneMinuseta - L01onePlusxi;		// N_2,1
      shp[0][2][i] =  L00onePluseta  + L01onePlusxi;		// N_3,1
      shp[0][3][i] = -L00onePluseta  + L01oneMinusxi;	// N_4,1

      shp[1][0][i] = -L10oneMinuseta - L11oneMinusxi;	// N_1,2
      shp[1][1][i] =  L10oneMinuseta - L11onePlusxi;		// N_2,2
      shp[1][2][i] =  L10onePluseta  + L11onePlusxi;		// N_3,2
      shp[1][3][i] = -L10onePluseta  + L11oneMinusxi;	// N_4,2

      dvol[i] = detJ * thickness * wts[i];
      vol += dvol[i];

	  for (k=0; k<2; k++) {
	    for (l=0; l<4; l++) {
	       shpBar[k][l] += shp[k][l][i] * dvol[i];
    	}
	  }
	}

	for (k=0; k<2; k++) {
	  for (l=0; l<4; l++) {
	    shpBar[k][l] /= vol;
	  }
	}

    // See Tom Hughes, Chapter 4: Mixed and Penalty Method
	for (i=0; i<4; i++) {
	  for (l=0; l<4; l++) {
	    B[0][0][l][i] = (2*shp[0][l][i]+shpBar[0][l])/3.;
	    B[0][1][l][i] = (shpBar[1][l] - shp[1][l][i])/3.;
	    B[1][0][l][i] = (shpBar[0][l] - shp[0][l][i])/3.;
	    B[1][1][l][i] = (2*shp[1][l][i]+shpBar[1][l])/3.;
	    B[2][0][l][i] = shp[1][l][i];
	    B[2][1][l][i] = shp[0][l][i];
	    B[3][0][l][i] = B[1][0][l][i];
	    B[3][1][l][i] = B[0][1][l][i];

	    Bp[0][l][i] = B[0][0][l][i];
	    Bp[1][l][i] = B[1][1][l][i];
	  }
	}
}


double BBarFourNodeQuadUP::mixtureRho(int i)
{
  double rhoi, e, n;
	rhoi= theMaterial[i]->getRho();
	e = 0.7;  //theMaterial[i]->getVoidRatio();
  n = e / (1.0 + e);
  //return n * rho + (1.0-n) * rhoi;
	return rhoi;
}

void BBarFourNodeQuadUP::setPressureLoadAtNodes(void)
{
    pressureLoad.Zero();

	if (pressure == 0.0)
		return;

	const Vector &node1 = nd1Ptr->getCrds();
	const Vector &node2 = nd2Ptr->getCrds();
	const Vector &node3 = nd3Ptr->getCrds();
	const Vector &node4 = nd4Ptr->getCrds();

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

	double pressureOver2 = pressure*thickness/2.0;

	// Contribution from side 12
	pressureLoad(0) += pressureOver2*dy12;
	pressureLoad(3) += pressureOver2*dy12;
	pressureLoad(1) += pressureOver2*-dx12;
	pressureLoad(4) += pressureOver2*-dx12;

	// Contribution from side 23
	pressureLoad(3) += pressureOver2*dy23;
	pressureLoad(6) += pressureOver2*dy23;
	pressureLoad(4) += pressureOver2*-dx23;
	pressureLoad(7) += pressureOver2*-dx23;

	// Contribution from side 34
	pressureLoad(6) += pressureOver2*dy34;
	pressureLoad(9) += pressureOver2*dy34;
	pressureLoad(7) += pressureOver2*-dx34;
	pressureLoad(10) += pressureOver2*-dx34;

	// Contribution from side 41
	pressureLoad(9) += pressureOver2*dy41;
	pressureLoad(0) += pressureOver2*dy41;
	pressureLoad(10) += pressureOver2*-dx41;
	pressureLoad(1) += pressureOver2*-dx41;

	//pressureLoad = pressureLoad*thickness;
}
