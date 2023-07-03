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

// Written: Rasool Ghorbani
// Created: 06/2023
// Revision: A
// Description: This file contains the class definition for GMG_CMAZ2D Element.
// This file contains the implementation for a modified version of ZeroLength class.
// The element only works with GMG_CMAC2D, GMG_CMAM2D, and GMG_CMACA2D material models.

#include "GMG_CMAZ2D.h"
#include <Information.h>

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
#include <elementAPI.h>
#include <vector>

/* *** Initializing the static variable, ine to return the tangent matrix and the other to return the resisting force **********************
// *** the dimension depends on model BasicBuilder which shows the spatial dimension of model and number of degrees-of-freedom at nodes **** */
// initialise the class wide variables
Matrix GMG_CMAZ2D::ZeroLengthM2(2, 2);
Matrix GMG_CMAZ2D::ZeroLengthM4(4, 4);
Matrix GMG_CMAZ2D::ZeroLengthM6(6, 6);
Matrix GMG_CMAZ2D::ZeroLengthM12(12, 12);
Vector GMG_CMAZ2D::ZeroLengthV2(2);
Vector GMG_CMAZ2D::ZeroLengthV4(4);
Vector GMG_CMAZ2D::ZeroLengthV6(6);
Vector GMG_CMAZ2D::ZeroLengthV12(12);

void* OPS_GMG_CMAZ2D()
{
	int ndm = OPS_GetNDM();

	//
	// first scan the command line to obtain eleID, iNode, jNode, material ID's
	// and their directions, and the orientation of ele xPrime and yPrime not
	// along the global x and y axis
	//

	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 5) {
		opserr << "WARNING too few arguments " <<
			"want - element ZeroLength eleTag? iNode? jNode? " <<
			"-mat matID1? ... -dir dirMat1? .. " <<
			"<-orient x1? x2? x3? y1? y2? y3?>\n";

		return 0;
	}

	// eleTag, iNode, jNode
	int idata[3];
	numdata = 3;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
		opserr << "WARNING: failed to get integer data\n";
		return 0;
	}

	// create an array of material pointers, to do this first count
	// the materials to create the array then get matID's and from ModelBuilder
	// obtain pointers to the material objects
	const char* type = OPS_GetString();
	if (strcmp(type, "-mat") != 0) {
		opserr << "WARNING expecting " <<
			"- element ZeroLength eleTag? iNode? jNode? " <<
			"-mat matID1? ... -dir dirMat1? .. " <<
			"<-orient x1? x2? x3? y1? y2? y3?>\n";

		return 0;
	}

	//    std::vector<UniaxialMaterial*> mats;
	// create the array
	ID matTags(0);
	int numMats = 0;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		int mtag;
		numdata = 1;
		// the first one not an int
		int numArgs = OPS_GetNumRemainingInputArgs();
		if (OPS_GetIntInput(&numdata, &mtag) < 0) {
			if (numArgs > OPS_GetNumRemainingInputArgs()) {
				// move current arg back by one
				OPS_ResetCurrentInputArg(-1);
			}
			break;
		}
		matTags[numMats] = mtag;
		numMats++;
	}
	UniaxialMaterial **theMats = new UniaxialMaterial *[numMats];
	UniaxialMaterial **theDampMats = new UniaxialMaterial *[numMats];

	for (int i = 0; i < numMats; i++) {

		theMats[i] = OPS_getUniaxialMaterial(matTags(i));
		theDampMats[i] = 0;

		if (theMats[i] == 0) {
			opserr << "WARNING no material " << matTags(i) <<
				"exitsts - element ZeroLength eleTag? iNode? jNode? " <<
				"-mat matID1? ... -dir dirMat1? .. " <<
				"<-orient x1? x2? x3? y1? y2? y3?>\n";
			return 0;
		}
	}

	/*
	// now read the dirn ID's for the materials added
	type = OPS_GetString();
	if (strcmp(type, "-dir") != 0) {
		opserr << "WARNING expecting -dir flag " <<
			"- element ZeroLength eleTag? iNode? jNode? " <<
			"-mat matID1? ... -dir dirMat1? .. " <<
			"<-orient x1? x2? x3? y1? y2? y3?>\n";
		return 0;
	}


	if (OPS_GetNumRemainingInputArgs() < numMats) {
		opserr << "WARNING not enough directions provided for ele " << idata[0] <<
			"- element ZeroLength eleTag? iNode? jNode? " <<
			"-mat matID1? ... -dir dirMat1? .. " <<
			"<-orient x1? x2? x3? y1? y2? y3?>\n";
		return 0;
	}
	*/

	ID dirs(numMats);
	/*
	if (OPS_GetIntInput(&numMats, &dirs(0)) < 0) {
		opserr << "WARNING invalid dir\n";
		return 0;
	}
	*/
	for (int i = 0; i < dirs.Size(); i++) {
		dirs(i)--; // subscrit to C++
	}

	// create the vectors for the element orientation
	Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
	Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;

	// finally check the command line to see if user specified orientation
	int doRayleighDamping = 0;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		type = OPS_GetString();
		if (strcmp(type, "-doRayleigh") == 0) {
			doRayleighDamping = 1;
			if (OPS_GetNumRemainingInputArgs() > 0) {
				numdata = 1;
				if (OPS_GetIntInput(&numdata, &doRayleighDamping) < 0) {
					opserr << "WARNING: invalid integer\n";
					return 0;
				}
			}
		}
		else 	if (strcmp(type, "-dampMats") == 0) {
			doRayleighDamping = 2;
			numdata = 1;
			int matType;
			for (int i = 0; i < numMats; i++) {
				// the first one not an int
				if (OPS_GetIntInput(&numdata, &matType) < 0) {
					UniaxialMaterial *theMat = OPS_getUniaxialMaterial(matType);
					if (theMat == 0) {
						opserr << "WARNING no damp material material " << matType << " for zeroLength ele: " << idata[0] << endln;
						return 0;
					}
					else {
						theDampMats[i] = theMat;
					}
				}
			}

		}
		else if (strcmp(type, "-orient") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 6) {
				opserr << "WARNING: insufficient orient values\n";
				return 0;
			}
			numdata = 3;
			if (OPS_GetDoubleInput(&numdata, &x(0)) < 0) {
				opserr << "WARNING: invalid double input\n";
				return 0;
			}
			if (OPS_GetDoubleInput(&numdata, &y(0)) < 0) {
				opserr << "WARNING: invalid double input\n";
				return 0;
			}
		}
	}

	Element *theEle = 0;
	if (doRayleighDamping != 2)
		theEle = new GMG_CMAZ2D(idata[0], ndm, idata[1], idata[2], x, y, numMats, theMats, dirs, doRayleighDamping);
	else
		theEle = new GMG_CMAZ2D(idata[0], ndm, idata[1], idata[2], x, y, numMats, theMats, theDampMats, dirs, doRayleighDamping);

	// return the memory we stole and return OK
	delete[] theMats;
	delete[] theDampMats;

	return theEle;
}


//  Constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the ZeroLength end nodes.


//  Construct element with one unidirectional material (numMaterials1d=1)
int direction = 1;
GMG_CMAZ2D::GMG_CMAZ2D(int tag,
	int dim,
	int Nd1, int Nd2,
	const Vector &x, const Vector &yp,
	UniaxialMaterial &theMat,
	int direction,
	int doRayleigh)
	:Element(tag, ELE_TAG_GMG_CMAZ2D),
	connectedExternalNodes(2),
	dimension(dim), numDOF(0), transformation(3, 3), useRayleighDamping(doRayleigh),
	theMatrix(0), theVector(0),
	numMaterials1d(1), theMaterial1d(0), dir1d(0), t1d(0), d0(0), v0(0)
{
	// allocate memory for numMaterials1d uniaxial material models
	theMaterial1d = new UniaxialMaterial*[numMaterials1d];
	dir1d = new ID(numMaterials1d);

	if (theMaterial1d == 0 || dir1d == 0) {
		opserr << "FATAL GMG_CMAZ2D::GMG_CMAZ2D - failed to create a 1d  material or direction array\n";
		exit(-1);
	}

	// initialize uniaxial materials and directions and check for valid values
	if (direction == 2 && dim == 2) // For Keri Ryan
		direction = 5;

	(*dir1d)(0) = direction;
	this->checkDirection(*dir1d);

	// get a copy of the material and check we obtained a valid copy
	theMaterial1d[0] = theMat.getCopy();
	if (theMaterial1d[0] == 0) {
		opserr << "FATAL ZeroLength::ZeroLength - failed to get a copy of material " << theMat.getTag() << endln;
		exit(-1);
	}

	// establish the connected nodes and set up the transformation matrix for orientation
	this->setUp(Nd1, Nd2, x, yp);

	// designate to setDomain that this is the initial construction of the element
	mInitialize = 1;
}


GMG_CMAZ2D::GMG_CMAZ2D(int tag,
	int dim,
	int Nd1, int Nd2,
	const Vector &x, const Vector &yp,
	UniaxialMaterial &theMat,
	UniaxialMaterial &theDampMat,
	int direction)
	:Element(tag, ELE_TAG_GMG_CMAZ2D),
	connectedExternalNodes(2),
	dimension(dim), numDOF(0), transformation(3, 3), useRayleighDamping(2),
	theMatrix(0), theVector(0),
	numMaterials1d(1), theMaterial1d(0), dir1d(0), t1d(0), d0(0), v0(0)
{
	// allocate memory for numMaterials1d uniaxial material models
	theMaterial1d = new UniaxialMaterial*[2];
	dir1d = new ID(numMaterials1d);

	if (theMaterial1d == 0 || dir1d == 0) {
		opserr << "FATAL ZeroLength::ZeroLength - failed to create a 1d  material or direction array\n";
		exit(-1);
	}

	// initialize uniaxial materials and directions and check for valid values
	if (direction == 2 && dim == 2) // For Keri Ryan
		direction = 5;

	(*dir1d)(0) = direction;
	this->checkDirection(*dir1d);

	// get a copy of the material and check we obtained a valid copy

	theMaterial1d[0] = theMat.getCopy();

	theMaterial1d[1] = theDampMat.getCopy();
	if (theMaterial1d[0] == 0 || theMaterial1d[1] == 0) {
		opserr << "FATAL ZeroLength::ZeroLength - failed to get a copy of material " << theMat.getTag() << endln;
		exit(-1);
	}

	// establish the connected nodes and set up the transformation matrix for orientation
	this->setUp(Nd1, Nd2, x, yp);

	// designate to setDomain that this is the initial construction of the element
	mInitialize = 1;
}




//  Construct element with multiple unidirectional materials
GMG_CMAZ2D::GMG_CMAZ2D(int tag,
	int dim,
	int Nd1, int Nd2,
	const Vector& x, const Vector& yp,
	int n1dMat,
	UniaxialMaterial** theMat,
	const ID& direction,
	int doRayleigh)
	:Element(tag, ELE_TAG_GMG_CMAZ2D),
	connectedExternalNodes(2),
	dimension(dim), numDOF(0), transformation(3, 3), useRayleighDamping(doRayleigh),
	theMatrix(0), theVector(0),
	numMaterials1d(n1dMat), theMaterial1d(0), dir1d(0), t1d(0), d0(0), v0(0)
{

	// allocate memory for numMaterials1d uniaxial material models
	theMaterial1d = new UniaxialMaterial*[numMaterials1d];
	dir1d = new ID(numMaterials1d);

	if (theMaterial1d == 0 || dir1d == 0) {
		opserr << "FATAL ZeroLength::ZeroLength - failed to create a 1d  material or direction array\n";
		exit(-1);
	}

	// initialize uniaxial materials and directions and check for valid values
	*dir1d = direction;
	for (int i = 0; i < n1dMat; i++) {
		if ((*dir1d)(i) == 2 && dim == 2) // For Keri Ryan
			(*dir1d)(i) = 5;
	}
	this->checkDirection(*dir1d);

	// get a copy of the material objects and check we obtained a valid copy
	for (int i = 0; i < numMaterials1d; i++) {
		theMaterial1d[i] = theMat[i]->getCopy();
		if (theMaterial1d[i] == 0) {
			opserr << "FATAL ZeroLength::ZeroLength - failed to get a copy of material " << theMat[i]->getTag() << endln;
			exit(-1);
		}
	}

	// establish the connected nodes and set up the transformation matrix for orientation
	this->setUp(Nd1, Nd2, x, yp);

	// designate to setDomain that this is the initial construction of the element
	mInitialize = 1;
}


//  Construct element with multiple unidirectional materials
GMG_CMAZ2D::GMG_CMAZ2D(int tag,
	int dim,
	int Nd1, int Nd2,
	const Vector& x, const Vector& yp,
	int n1dMat,
	UniaxialMaterial** theMat,
	UniaxialMaterial** theDampMat,
	const ID& direction,
	int doRayleigh)
	:Element(tag, ELE_TAG_GMG_CMAZ2D),
	connectedExternalNodes(2),
	dimension(dim), numDOF(0), transformation(3, 3), useRayleighDamping(doRayleigh),
	theMatrix(0), theVector(0),
	numMaterials1d(n1dMat), theMaterial1d(0), dir1d(0), t1d(0), d0(0), v0(0)
{

	// allocate memory for numMaterials1d uniaxial material models
	theMaterial1d = new UniaxialMaterial*[2 * numMaterials1d];
	dir1d = new ID(numMaterials1d);

	if (theMaterial1d == 0 || dir1d == 0) {
		opserr << "FATAL ZeroLength::ZeroLength - failed to create a 1d  material or direction array\n";
		exit(-1);
	}

	// initialize uniaxial materials and directions and check for valid values
	*dir1d = direction;
	for (int i = 0; i < n1dMat; i++) {
		if ((*dir1d)(i) == 2 && dim == 2) // For Keri Ryan
			(*dir1d)(i) = 5;
	}
	this->checkDirection(*dir1d);

	// get a copy of the material objects and check we obtained a valid copy
	for (int i = 0; i < numMaterials1d; i++) {
		theMaterial1d[i] = theMat[i]->getCopy();
		theMaterial1d[i + numMaterials1d] = theDampMat[i]->getCopy();
		if (theMaterial1d[i] == 0) {
			opserr << "FATAL ZeroLength::ZeroLength - failed to get a copy of material " << theMat[i]->getTag() << endln;
			exit(-1);
		}
	}

	// establish the connected nodes and set up the transformation matrix for orientation
	this->setUp(Nd1, Nd2, x, yp);

	// designate to setDomain that this is the initial construction of the element
	mInitialize = 1;
}


//   Constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
GMG_CMAZ2D::GMG_CMAZ2D(void)
	:Element(0, ELE_TAG_ZeroLength),
	connectedExternalNodes(2),
	dimension(0), numDOF(0), transformation(3, 3),
	theMatrix(0), theVector(0),
	numMaterials1d(0), theMaterial1d(0),
	dir1d(0), t1d(0), d0(0), v0(0)
{
	// ensure the connectedExternalNode ID is of correct size 
	if (connectedExternalNodes.Size() != 2)
		opserr << "FATAL ZeroLength::ZeroLength - failed to create an ID of correct size\n";

	// designate to setDomain that this is the null construction of the element
	mInitialize = 0;
}


//  Destructor:
//  delete must be invoked on any objects created by the object
//  and on the matertial object.
GMG_CMAZ2D::~GMG_CMAZ2D()
{
	// invoke the destructor on any objects created by the object
	// that the object still holds a pointer to

	// invoke destructors on material objects
	int numMat = numMaterials1d;
	if (useRayleighDamping == 2)
		numMat *= 2;
	for (int mat = 0; mat < numMat; mat++)
		delete theMaterial1d[mat];

	// delete memory of 1d materials    
	if (theMaterial1d != 0)
		delete[] theMaterial1d;

	if (t1d != 0)
		delete t1d;
	if (dir1d != 0)
		delete dir1d;

	if (d0 != 0)
		delete d0;

	if (v0 != 0)
		delete v0;
}


int
GMG_CMAZ2D::getNumExternalNodes(void) const
{
	return 2;
}


const ID &
GMG_CMAZ2D::getExternalNodes(void)
{
	return connectedExternalNodes;
}



Node **
GMG_CMAZ2D::getNodePtrs(void)
{
	return theNodes;
}

int
GMG_CMAZ2D::getNumDOF(void)
{
	return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the ZeroLength element, we set matrix and vector pointers,
//    allocate space for t matrix and define it as the basic deformation-
//    displacement transformation matrix.
void
GMG_CMAZ2D::setDomain(Domain *theDomain)
{
	// check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		return;
	}

	// set default values for error conditions
	numDOF = 2;
	theMatrix = &ZeroLengthM2;
	theVector = &ZeroLengthV2;

	// first set the node pointers
	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);
	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);

	// if can't find both - send a warning message
	if (theNodes[0] == 0 || theNodes[1] == 0) {
		if (theNodes[0] == 0)
			opserr << "WARNING ZeroLength::setDomain() - Nd1: " << Nd1 << " does not exist in ";
		else
			opserr << "WARNING ZeroLength::setDomain() - Nd2: " << Nd2 << " does not exist in ";

		opserr << "model for ZeroLength ele: " << this->getTag() << endln;

		return;
	}

	// now determine the number of dof and the dimension    
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();

	// if differing dof at the ends - print a warning message
	if (dofNd1 != dofNd2) {
		opserr << "WARNING ZeroLength::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
			"have differing dof at ends for ZeroLength " << this->getTag() << endln;
		return;
	}

	// Check that length is zero within tolerance
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();
	Vector diff = end1Crd - end2Crd;
	double L = diff.Norm();
	double v1 = end1Crd.Norm();
	double v2 = end2Crd.Norm();
	double vm;

	vm = (v1 < v2) ? v2 : v1;


	if (L > LENTOL*vm)
		opserr << "WARNING ZeroLength::setDomain(): Element " << this->getTag() << " has L= " << L <<
		", which is greater than the tolerance\n";

	// call the base class method
	this->DomainComponent::setDomain(theDomain);

	// set the number of dof for element and set matrix and vector pointer
	if (dimension == 1 && dofNd1 == 1) {
		numDOF = 2;
		theMatrix = &ZeroLengthM2;
		theVector = &ZeroLengthV2;
		elemType = D1N2;
	}
	else if (dimension == 2 && dofNd1 == 2) {
		numDOF = 4;
		theMatrix = &ZeroLengthM4;
		theVector = &ZeroLengthV4;
		elemType = D2N4;
	}
	else if (dimension == 2 && dofNd1 == 3) {
		numDOF = 6;
		theMatrix = &ZeroLengthM6;
		theVector = &ZeroLengthV6;
		elemType = D2N6;
	}
	else if (dimension == 3 && dofNd1 == 3) {
		numDOF = 6;
		theMatrix = &ZeroLengthM6;
		theVector = &ZeroLengthV6;
		elemType = D3N6;
	}
	else if (dimension == 3 && dofNd1 == 6) {
		numDOF = 12;
		theMatrix = &ZeroLengthM12;
		theVector = &ZeroLengthV12;
		elemType = D3N12;
	}
	else {
		opserr << "WARNING ZeroLength::setDomain cannot handle " << dimension <<
			"dofs at nodes in " << dofNd1 << " d problem\n";
		return;
	}

	// create the basic deformation-displacement transformation matrix for the element
	// for 1d materials (uniaxial materials)
	if (numMaterials1d > 0)
		this->setTran1d(elemType, numMaterials1d);

	// get trial displacements and take difference
	const Vector& disp1 = theNodes[0]->getTrialDisp();
	const Vector& disp2 = theNodes[1]->getTrialDisp();
	Vector  diffD = disp2 - disp1;
	const Vector& vel1 = theNodes[0]->getTrialVel();
	const Vector& vel2 = theNodes[1]->getTrialVel();
	Vector  diffV = vel2 - vel1;

	// to avoid incorrect results, do not set initial disp/vel upon call of null constructor
	// when using database commands
	if (mInitialize == 1) {
		if (diffD != 0.0)
			d0 = new Vector(diffD);

		if (diffV != 0)
			v0 = new Vector(diffV);
	}
}



int
GMG_CMAZ2D::commitState()
{
	int code = 0;

	// call element commitState to do any base class stuff
	if ((code = this->Element::commitState()) != 0) {
		opserr << "ZeroLength::commitState () - failed in base class";
	}

	// commit 1d materials
	int numMat = numMaterials1d;
	if (useRayleighDamping == 2)
		numMat *= 2;
	for (int i = 0; i < numMat; i++)
		code += theMaterial1d[i]->commitState();

	return code;
}

int
GMG_CMAZ2D::revertToLastCommit()
{
	int code = 0;

	// revert state for 1d materials
	int numMat = numMaterials1d;
	if (useRayleighDamping == 2)
		numMat *= 2;
	for (int i = 0; i < numMat; i++)
		code += theMaterial1d[i]->revertToLastCommit();

	return code;
}


int
GMG_CMAZ2D::revertToStart()
{
	int code = 0;

	// revert to start for 1d materials
	int numMat = numMaterials1d;
	if (useRayleighDamping == 2)
		numMat *= 2;
	for (int i = 0; i < numMat; i++)
		code += theMaterial1d[i]->revertToStart();

	return code;
}


int
GMG_CMAZ2D::update(void)
{
	double strain;
	double strainRate;
	Vector strain_vector;
	strain_vector.Zero();
	// get trial displacements and take difference
	const Vector& disp1 = theNodes[0]->getTrialDisp();
	const Vector& disp2 = theNodes[1]->getTrialDisp();
	Vector  diff = disp2 - disp1;
	const Vector& vel1 = theNodes[0]->getTrialVel();
	const Vector& vel2 = theNodes[1]->getTrialVel();
	Vector  diffv = vel2 - vel1;
	//opserr << "Moon4 " << "  " << " disp1 " << disp1 << "  " << " disp2 " << disp2 << "  " << " diff " << diff << "  " << " vel1 " << vel1 << "  " << " vel2 " << vel2 << "  " << " diffv " << diffv << endln;

	if (d0 != 0)
		diff -= *d0;

	if (v0 != 0)
		diffv -= *v0;

	// loop over 1d materials

	//    Matrix& tran = *t1d;

	//Vector &def = *strain_vector;

	int ret = 0;
	//for (int mat = 0; mat < numMaterials1d; mat++) {
		// compute strain and rate; set as current trial for material
	//strain_vector(0) += (-1) * diff(0);
	//strain_vector(1) += (-1) * diff(1);
	//strain_vector(2) += (-1) * diff(2);
		//strain = this->computeCurrentStrain1d(mat, diff);
		//*strain_vector = this->computeCurrentStrain1d(0, diff);
		//strainRate = this->computeCurrentStrain1d(mat, diffv);
		//ret += theMaterial1d[mat]->setTrialStrain(strain, strainRate);
	//diff *= -1;
	//opserr << "Moon4_1" << "  " << " disp1 " << disp1 << "  " << " disp2 " << disp2 << "  " << " diff " << diff << "  " << " vel1 " << vel1 << "  " << " vel2 " << vel2 << "  " << " diffv " << diffv << endln;
	//ret += theMaterial1d[0]->setTrialStrain(diff);
	//if (useRayleighDamping == 2) {
	//	ret += theMaterial1d[mat + numMaterials1d]->setTrialStrain(strainRate);
	//}
//}


		//for (int mat = 0; mat < numMaterials1d; mat++) {
		// compute strain and rate; set as current trial for material
	strain_vector = this->computeCurrentStrain1d(0, diff);
	//opserr << "Moon4_1" << "  " << " disp1 " << disp1 << "  " << " disp2 " << disp2 << "  " << " diff " << diff << "  " << " vel1 " << vel1 << "  " << " vel2 " << vel2 << "  " << " diffv " << diffv << endln;
	ret += theMaterial1d[0]->setTrialStrain(strain_vector);
	//if (useRayleighDamping == 2) {
		//ret += theMaterial1d[mat + numMaterials1d]->setTrialStrain(strainRate);
	//}
  //}

	return ret;
}

const Matrix &
GMG_CMAZ2D::getTangentStiff(void)
{
	double E = 0.0;
	//double E_Shear, E_Axial, E_Rotation;

	// stiff is a reference to the matrix holding the stiffness matrix
	Matrix& stiff = *theMatrix;

	// zero stiffness matrix
	stiff.Zero();

	// loop over 1d materials
	Matrix& tran = *t1d;
	//Matrix& tran = { {-1,0,0,1,0,0}, {0,-1,0,0,1,0}, {0,0,-1,0,0,1} };

	//E_Shear = theMaterial1d[0]->getTangent();
	//E_Axial = theMaterial1d[0]->getTangent();
	//E_Rotation = theMaterial1d[0]->getTangent();
	//opserr << "Moon1 " << "  " << " E_Shear " << E_Shear << " E_Axial " << E_Axial << "E_Rotation" << E_Rotation << endln;


	/*
	for (int mat = 0; mat < numMaterials1d; mat++) {

		// get tangent for material
		E = theMaterial1d[mat]->getTangent();

		// compute contribution of material to tangent matrix
		for (int i = 0; i < numDOF; i++)
			for (int j = 0; j < i + 1; j++)
				stiff(i, j) += tran(mat, i) * E * tran(mat, j);

	}

	// end loop over 1d materials

	// complete symmetric stiffness matrix
	for (int i = 0; i < numDOF; i++)
		for (int j = 0; j < i; j++)
			stiff(j, i) = stiff(i, j);
	*/
	for (int mat = 0; mat < 3; mat++) {

		// get tangent for material
		E = theMaterial1d[0]->getTangent();

		// compute contribution of material to tangent matrix
		for (int i = 0; i < numDOF; i++)
			for (int j = 0; j < i + 1; j++)
				stiff(i, j) += tran(mat, i) * E * tran(mat, j);

	}

	// end loop over 1d materials

	// complete symmetric stiffness matrix
	for (int i = 0; i < numDOF; i++)
		for (int j = 0; j < i; j++)
			stiff(j, i) = stiff(i, j);

	//stiff(0, 0) = E_Shear;  stiff(0, 1) = 0.0;  stiff(0, 2) = 0.0;  stiff(0, 3) = -E_Shear;  stiff(0, 4) = 0.0;  stiff(0, 5) = 0.0;
	//stiff(1, 0) = 0.0;   stiff(1, 1) = E_Axial;  stiff(1, 2) = 0.0;  stiff(1, 3) = 0.0;  stiff(1, 4) = -E_Axial;  stiff(1, 5) = 0.0;
	//stiff(2, 0) = 0.0;   stiff(2, 1) = 0.0;  stiff(2, 2) = E_Rotation;  stiff(2, 3) = 0.0;  stiff(2, 4) = 0.0;  stiff(2, 5) = -E_Rotation;
	//stiff(3, 0) = -E_Shear;  stiff(3, 1) = 0.0;  stiff(3, 2) = 0.0;  stiff(3, 3) = E_Shear;  stiff(3, 4) = 0.0;  stiff(3, 5) = 0.0;
	//stiff(4, 0) = 0.0;  stiff(4, 1) = -E_Axial;  stiff(4, 2) = 0.0;  stiff(4, 3) = 0.0;  stiff(4, 4) = E_Axial;  stiff(4, 5) = 0.0;
	//stiff(5, 0) = 0.0;  stiff(5, 1) = 0.0;  stiff(5, 2) = -E_Rotation;  stiff(5, 3) = 0.0;  stiff(5, 4) = 0.0;  stiff(5, 5) = E_Rotation;

	//opserr << " Moon1_2 " << "  " << " stiff " << "  " << stiff << endln;

	return stiff;
}


const Matrix &
GMG_CMAZ2D::getInitialStiff(void)
{
	double E = 0.0;
	//double E_Shear, E_Axial, E_Rotation;
	// stiff is a reference to the matrix holding the stiffness matrix
	Matrix& stiff = *theMatrix;

	// zero stiffness matrix
	stiff.Zero();

	// loop over 1d materials

	Matrix& tran = *t1d;

	//E_Shear = theMaterial1d[0]->getInitialTangent();
	//E_Axial = theMaterial1d[0]->getInitialTangent();
	//E_Rotation = theMaterial1d[0]->getInitialTangent();
	//opserr << "Moon2 " << "  " << " E_Shear " << E_Shear << " E_Axial " << E_Axial << "E_Rotation" << E_Rotation << endln;

	/*
	for (int mat = 0; mat < numMaterials1d; mat++) {

		// get tangent for material
		E = theMaterial1d[mat]->getInitialTangent();

		// compute contribution of material to tangent matrix
		for (int i = 0; i < numDOF; i++)
			for (int j = 0; j < i + 1; j++)
				stiff(i, j) += tran(mat, i) * E * tran(mat, j);

	}

	// end loop over 1d materials

	// complete symmetric stiffness matrix
	for (int i = 0; i < numDOF; i++)
		for (int j = 0; j < i; j++)
			stiff(j, i) = stiff(i, j);
	*/

	for (int mat = 0; mat < 3; mat++) {

		// get tangent for material
		E = theMaterial1d[0]->getTangent();

		// compute contribution of material to tangent matrix
		for (int i = 0; i < numDOF; i++)
			for (int j = 0; j < i + 1; j++)
				stiff(i, j) += tran(mat, i) * E * tran(mat, j);

	}

	// end loop over 1d materials

	// complete symmetric stiffness matrix
	for (int i = 0; i < numDOF; i++)
		for (int j = 0; j < i; j++)
			stiff(j, i) = stiff(i, j);

	//stiff(0, 0) = E_Shear;  stiff(0, 1) = 0.0;  stiff(0, 2) = 0.0;  stiff(0, 3) = -E_Shear;  stiff(0, 4) = 0.0;  stiff(0, 5) = 0.0;
	//stiff(1, 0) = 0.0;  stiff(1, 1) = E_Axial;  stiff(1, 2) = 0.0;  stiff(1, 3) = 0.0;  stiff(1, 4) = -E_Axial;  stiff(1, 5) = 0.0;
	//stiff(2, 0) = 0.0;  stiff(2, 1) = 0.0;  stiff(2, 2) = E_Rotation;  stiff(2, 3) = 0.0;  stiff(2, 4) = 0.0;  stiff(2, 5) = -E_Rotation;
	//stiff(3, 0) = -E_Shear;  stiff(3, 1) = 0.0;  stiff(3, 2) = 0.0;  stiff(3, 3) = E_Shear;  stiff(3, 4) = 0.0;  stiff(3, 5) = 0.0;
	//stiff(4, 0) = 0.0;  stiff(4, 1) = -E_Axial;  stiff(4, 2) = 0.0;  stiff(4, 3) = 0.0;  stiff(4, 4) = E_Axial;  stiff(4, 5) = 0.0;
	//stiff(5, 0) = 0.0;  stiff(5, 1) = 0.0;  stiff(5, 2) = -E_Rotation;  stiff(5, 3) = 0.0;  stiff(5, 4) = 0.0;  stiff(5, 5) = E_Rotation;

	//opserr << " Moon2_2 " << "  " << " stiff " << "  " << stiff << endln;

	return stiff;
}


const Matrix &
GMG_CMAZ2D::getDamp(void)
{
	// damp is a reference to the matrix holding the damping matrix
	Matrix& damp = *theMatrix;

	// zero damping matrix
	damp.Zero();

	// get Rayleigh damping matrix 

	if (useRayleighDamping == 1) {

		damp = this->Element::getDamp();

	}
	else if (useRayleighDamping == 2) {

		// loop over 1d materials and add their damping tangents
		double eta;
		Matrix& tran = *t1d;
		for (int mat = 0; mat < 3; mat++) {

			// get tangent for material
			eta = theMaterial1d[mat + 3]->getTangent();

			// compute contribution of material to tangent matrix
			for (int i = 0; i < numDOF; i++)
				for (int j = 0; j < i + 1; j++)
					damp(i, j) += tran(mat, i) * eta * tran(mat, j);
		}

	}
	else {

		// loop over 1d materials and add their damping tangents
		double eta;
		Matrix& tran = *t1d;;
		for (int mat = 0; mat < 3; mat++) {

			// get tangent for material
			eta = theMaterial1d[0]->getDampTangent();

			// compute contribution of material to tangent matrix
			for (int i = 0; i < numDOF; i++)
				for (int j = 0; j < i + 1; j++)
					damp(i, j) += tran(mat, i) * eta * tran(mat, j);

		} // end loop over 1d materials 

	}

	// complete symmetric damping matrix
	for (int i = 0; i < numDOF; i++)
		for (int j = 0; j < i; j++)
			damp(j, i) = damp(i, j);

	return damp;
}


const Matrix &
GMG_CMAZ2D::getMass(void)
{
	// no mass 
	theMatrix->Zero();
	return *theMatrix;
}


void
GMG_CMAZ2D::zeroLoad(void)
{
	// does nothing now
}

int
GMG_CMAZ2D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	opserr << "ZeroLength::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;

	return -1;
}

int
GMG_CMAZ2D::addInertiaLoadToUnbalance(const Vector &accel)
{
	// does nothing as element has no mass yet!
	return 0;
}


const Vector &
GMG_CMAZ2D::getResistingForce()
{
	double force;
	//double force_Shear, force_Axial, force_Rotation;

	// zero the residual
	theVector->Zero();
	//opserr << "Moon3_2 " << "  " << " theVector " << *theVector << endln;

	//force_Shear = theMaterial1d[0]->getStress();
	//force_Axial = theMaterial1d[0]->getStress();
	//force_Rotation = theMaterial1d[0]->getStress();
	//opserr << "Moon3 " << "  " << " force_Shear " << force_Shear << " force_Axial " << force_Axial << "force_Rotation" << force_Rotation << endln;


	// loop over 1d materials
	for (int mat = 0; mat < 3; mat++) {

		// get resisting force for material
		force = theMaterial1d[0]->getStress();

		// compute residual due to resisting force
		for (int i = 0; i < numDOF; i++)
			(*theVector)(i) += (*t1d)(mat, i) * force;

	} // end loop over 1d materials


	/*
	(*theVector)(0) += (-1) * force_Shear; (*theVector)(3) += (1) * force_Shear;
	(*theVector)(1) += (-1) * force_Axial; (*theVector)(4) += (1) * force_Axial;
	(*theVector)(2) += (-1) * force_Rotation; (*theVector)(5) += (1) * force_Rotation;
	*/
	//opserr << "Moon3_1 " << "  " << " theVector " << *theVector << "  " << " *t1d " << *t1d << endln;

	return *theVector;
}


const Vector &
GMG_CMAZ2D::getResistingForceIncInertia()
{
	// this already includes damping forces from materials
	this->getResistingForce();

	// add the damping forces from rayleigh damping
	if (useRayleighDamping == 1) {
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
			*theVector += this->getRayleighDampingForces();
			//opserr << "Moon5-1 " << "  " << " theVector " << *theVector << endln;
		}
	}
	else if (useRayleighDamping == 2) {
		// loop over 1d materials
		for (int mat = 0; mat < 3; mat++) {

			// get resisting force for material
			double force = theMaterial1d[mat + 3]->getStress();

			// compute residual due to resisting force
			for (int i = 0; i < numDOF; i++)
				(*theVector)(i) += (*t1d)(mat, i) * force;

			//opserr << "Moon5 " << "  " << " theVector " << *theVector << endln;
		}
	}
	return *theVector;
}


int
GMG_CMAZ2D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// ZeroLength packs its data into an ID and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments

	// Make one size bigger so not a multiple of 3, otherwise will conflict
	// with classTags ID
	static ID idData(7);

	idData(0) = this->getTag();
	idData(1) = dimension;
	idData(2) = numDOF;
	idData(3) = numMaterials1d;
	idData(4) = connectedExternalNodes(0);
	idData(5) = connectedExternalNodes(1);
	idData(6) = useRayleighDamping;

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "ZeroLength::sendSelf -- failed to send ID data\n";
		return res;
	}

	// Send the 3x3 direction cosine matrix, have to send it since it is only set
	// in the constructor and not setDomain()
	res += theChannel.sendMatrix(dataTag, commitTag, transformation);
	if (res < 0) {
		opserr << "ZeroLength::sendSelf -- failed to send transformation Matrix\n";
		return res;
	}

	if (numMaterials1d < 1)
		return res;
	else {
		ID classTags(numMaterials1d * 3);

		int i;
		// Loop over the materials and send them
		for (i = 0; i < numMaterials1d; i++) {
			int matDbTag = theMaterial1d[i]->getDbTag();
			if (matDbTag == 0) {
				matDbTag = theChannel.getDbTag();
				if (matDbTag != 0)
					theMaterial1d[i]->setDbTag(matDbTag);
			}
			classTags(i) = matDbTag;
			classTags(numMaterials1d + i) = theMaterial1d[i]->getClassTag();
			classTags(2 * numMaterials1d + i) = (*dir1d)(i);
		}

		res += theChannel.sendID(dataTag, commitTag, classTags);
		if (res < 0) {
			opserr << " ZeroLength::sendSelf -- failed to send classTags ID\n";
			return res;
		}

		for (i = 0; i < numMaterials1d; i++) {
			res += theMaterial1d[i]->sendSelf(commitTag, theChannel);
			if (res < 0) {
				opserr << "ZeroLength::sendSelf -- failed to send Material1d " << i << endln;

				return res;
			}
		}
	}

	return res;
}

int
GMG_CMAZ2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// ZeroLength creates an ID, receives the ID and then sets the 
	// internal data with the data in the ID

	static ID idData(7);

	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "ZeroLength::recvSelf -- failed to receive ID data\n";

		return res;
	}

	res += theChannel.recvMatrix(dataTag, commitTag, transformation);
	if (res < 0) {
		opserr << "ZeroLength::recvSelf -- failed to receive transformation Matrix\n";

		return res;
	}

	this->setTag(idData(0));
	dimension = idData(1);
	numDOF = idData(2);
	connectedExternalNodes(0) = idData(4);
	connectedExternalNodes(1) = idData(5);
	useRayleighDamping = idData(6);

	if (idData(3) < 1) {
		numMaterials1d = 0;
		if (dir1d != 0) {
			delete dir1d;
			dir1d = 0;
		}
		return res;
	}
	else {
		// Check that there is correct number of materials, reallocate if needed
		if (numMaterials1d != idData(3)) {
			int i;
			if (theMaterial1d != 0) {
				for (i = 0; i < numMaterials1d; i++)
					delete theMaterial1d[i];
				delete[] theMaterial1d;
				theMaterial1d = 0;
			}

			numMaterials1d = idData(3);

			theMaterial1d = new UniaxialMaterial *[numMaterials1d];
			if (theMaterial1d == 0) {
				opserr << "ZeroLength::recvSelf -- failed to new Material1d array\n";
				return -1;
			}

			for (i = 0; i < numMaterials1d; i++)
				theMaterial1d[i] = 0;

			// Allocate ID array for directions
			if (dir1d != 0)
				delete dir1d;
			dir1d = new ID(numMaterials1d);
			if (dir1d == 0) {
				opserr << "ZeroLength::recvSelf -- failed to new dir ID\n";

				return -1;
			}
		}

		ID classTags(3 * numMaterials1d);
		res += theChannel.recvID(dataTag, commitTag, classTags);
		if (res < 0) {
			opserr << "ZeroLength::recvSelf -- failed to receive classTags ID\n";
			return res;
		}

		for (int i = 0; i < numMaterials1d; i++) {
			int matClassTag = classTags(numMaterials1d + i);

			// If null, get a new one from the broker
			if (theMaterial1d[i] == 0)
				theMaterial1d[i] = theBroker.getNewUniaxialMaterial(matClassTag);

			// If wrong type, get a new one from the broker
			if (theMaterial1d[i]->getClassTag() != matClassTag) {
				delete theMaterial1d[i];
				theMaterial1d[i] = theBroker.getNewUniaxialMaterial(matClassTag);
			}

			// Check if either allocation failed from broker
			if (theMaterial1d[i] == 0) {
				opserr << "ZeroLength::recvSelf  -- failed to allocate new Material1d " << i << endln;
				return -1;
			}

			// Receive the materials
			theMaterial1d[i]->setDbTag(classTags(i));
			res += theMaterial1d[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "ZeroLength::recvSelf  -- failed to receive new Material1d " << i << endln;
				return res;
			}

			// Set material directions
			(*dir1d)(i) = classTags(2 * numMaterials1d + i);
		}
	}

	return res;
}


int
GMG_CMAZ2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// ensure setDomain() worked
	if (theNodes[0] == 0 || theNodes[1] == 0)
		return 0;

	static Vector v1(3);
	static Vector v2(3);

	float d1 = 1.0;
	float d2 = 1.0;
	/*
	if (displayMode == 1 || displayMode == 2) {

		theNodes[0]->getDisplayCrds(v1, fact);
		theNodes[1]->getDisplayCrds(v2, fact);

		if (displayMode == 1)
			d1 = theMaterial1d[0]->getStress();
		else
			d1 = theMaterial1d[0]->getStrain();
	}
	else {

		theNodes[0]->getDisplayCrds(v1, 0.);
		theNodes[1]->getDisplayCrds(v2, 0.);

	}
	*/
	if (v1 != v2)
		return theViewer.drawLine(v1, v2, d1, d1);
	else
		return theViewer.drawPoint(v1, d1, 10);
}


void
GMG_CMAZ2D::Print(OPS_Stream &s, int flag)
{
	// compute the strain and axial force in the member
	double strain = 0.0;
	double force = 0.0;

	for (int i = 0; i < numDOF; i++)
		(*theVector)(i) = (*t1d)(0, i)*force;

	if (flag == OPS_PRINT_CURRENTSTATE) { // print everything
		s << "Element: " << this->getTag();
		s << " type: ZeroLength  iNode: " << connectedExternalNodes(0);
		s << " jNode: " << connectedExternalNodes(1) << endln;
		for (int j = 0; j < numMaterials1d; j++) {
			s << "\tMaterial1d, tag: " << theMaterial1d[j]->getTag()
				<< ", dir: " << (*dir1d)(j) << endln;
			s << *(theMaterial1d[j]);
		}
		if (useRayleighDamping == 2) {
			s << "Damping Materials:\n";
			for (int j = numMaterials1d; j < 2 * numMaterials1d; j++) {
				s << "\tMaterial1d, tag: " << theMaterial1d[j]->getTag()
					<< ", dir: " << (*dir1d)(j) << endln;
				s << *(theMaterial1d[j]);
			}
		}
	}

	else if (flag == 1) {
		s << this->getTag() << "  " << strain << "  ";
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"ZeroLength\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
		s << "\"materials\": [";
		for (int i = 0; i < numMaterials1d - 1; i++)
			s << "\"" << theMaterial1d[i]->getTag() << "\", ";
		s << "\"" << theMaterial1d[numMaterials1d - 1]->getTag() << "\"], ";
		s << "\"dof\": [";
		for (int i = 0; i < numMaterials1d - 1; i++) {
			if ((*dir1d)(i) == 0)
				s << "\"P\", ";
			else if ((*dir1d)(i) == 1)
				s << "\"Vy\", ";
			else if ((*dir1d)(i) == 2)
				s << "\"Vz\", ";
			else if ((*dir1d)(i) == 3)
				s << "\"T\", ";
			else if ((*dir1d)(i) == 4)
				s << "\"My\", ";
			else if ((*dir1d)(i) == 5)
				s << "\"Mz\", ";
		}
		if ((*dir1d)(numMaterials1d - 1) == 0)
			s << "\"P\"], ";
		else if ((*dir1d)(numMaterials1d - 1) == 1)
			s << "\"Vy\"], ";
		else if ((*dir1d)(numMaterials1d - 1) == 2)
			s << "\"Vz\"], ";
		else if ((*dir1d)(numMaterials1d - 1) == 3)
			s << "\"T\"], ";
		else if ((*dir1d)(numMaterials1d - 1) == 4)
			s << "\"My\"], ";
		else if ((*dir1d)(numMaterials1d - 1) == 5)
			s << "\"Mz\"], ";
		s << "\"transMatrix\": [[";
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (j < 2)
					s << transformation(i, j) << ", ";
				else if (j == 2 && i < 2)
					s << transformation(i, j) << "], [";
				else if (j == 2 && i == 2)
					s << transformation(i, j) << "]]}";
			}
		}
	}
}

Response*
GMG_CMAZ2D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "ZeroLength");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);

	char outputData[10];

	if ((strcmp(argv[0], "force") == 0) || (strcmp(argv[0], "forces") == 0)
		|| (strcmp(argv[0], "globalForces") == 0) || (strcmp(argv[0], "globalforces") == 0)) {

		char outputData[10];
		int numDOFperNode = numDOF / 2;
		for (int i = 0; i < numDOFperNode; i++) {
			sprintf(outputData, "P1_%d", i + 1);
			output.tag("ResponseType", outputData);
		}
		for (int j = 0; j < numDOFperNode; j++) {
			sprintf(outputData, "P2_%d", j + 1);
			output.tag("ResponseType", outputData);
		}
		theResponse = new ElementResponse(this, 1, Vector(numDOF));

	}
	else if ((strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) ||
		(strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0)) {

		for (int i = 0; i < 3; i++) {
			sprintf(outputData, "P%d", i + 1);
			output.tag("ResponseType", outputData);
		}
		theResponse = new ElementResponse(this, 2, Vector(3));

	}
	else if (strcmp(argv[0], "defo") == 0 || strcmp(argv[0], "deformations") == 0 ||
		strcmp(argv[0], "deformation") == 0 || strcmp(argv[0], "basicDeformation") == 0) {

		for (int i = 0; i < 3; i++) {
			sprintf(outputData, "e%d", i + 1);
			output.tag("ResponseType", outputData);
		}
		theResponse = new ElementResponse(this, 3, Vector(3));

	}
	else if (strcmp(argv[0], "basicStiffness") == 0) {

		for (int i = 0; i < 3; i++) {
			sprintf(outputData, "e%d", i + 1);
			output.tag("ResponseType", outputData);
		}
		theResponse = new ElementResponse(this, 13, Matrix(3, 3));

	}
	else if ((strcmp(argv[0], "defoANDforce") == 0) ||
		(strcmp(argv[0], "deformationANDforces") == 0) ||
		(strcmp(argv[0], "deformationsANDforces") == 0)) {

		int i;
		for (i = 0; i < 3; i++) {
			sprintf(outputData, "e%d", i + 1);
			output.tag("ResponseType", outputData);
		}
		for (i = 0; i < 3; i++) {
			sprintf(outputData, "P%d", i + 1);
			output.tag("ResponseType", outputData);
		}
		theResponse = new ElementResponse(this, 4, Vector(2 * 3));


		// a material quantity
	}
	else if (strcmp(argv[0], "material") == 0) {
		if (argc > 2) {
			int matNum = atoi(argv[1]);
			int numMat = 3;
			if (useRayleighDamping == 2)
				numMat *= 2;
			if (matNum >= 1 && matNum <= numMat)
				theResponse = theMaterial1d[matNum - 1]->setResponse(&argv[2], argc - 2, output);
		}
	}

	if ((strcmp(argv[0], "dampingForces") == 0) || (strcmp(argv[0], "rayleighForces") == 0)) {
		theResponse = new ElementResponse(this, 15, Vector(numDOF));
	}

	output.endTag();




	return theResponse;
}

int
GMG_CMAZ2D::getResponse(int responseID, Information &eleInformation)
{
	const Vector& disp1 = theNodes[0]->getTrialDisp();
	const Vector& disp2 = theNodes[1]->getTrialDisp();
	const Vector  diff = disp2 - disp1;

	switch (responseID) {
	case -1:
		return -1;

	case 1:
		return eleInformation.setVector(this->getResistingForce());

	case 15:
		theVector->Zero();
		if (useRayleighDamping == 1) {
			if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
				*theVector += this->getRayleighDampingForces();
		}
		else if (useRayleighDamping == 2) {
			for (int mat = 0; mat < 3; mat++) {

				// get resisting force for material
				double force = theMaterial1d[mat + 3]->getStress();

				// compute residual due to resisting force
				for (int i = 0; i < numDOF; i++)
					(*theVector)(i) += (*t1d)(mat, i) * force;
			}
		}
		return eleInformation.setVector(*theVector);

	case 2:
		if (eleInformation.theVector != 0) {
			for (int i = 0; i < 3; i++)
				(*(eleInformation.theVector))(i) = theMaterial1d[0]->getStress();
			//(*(eleInformation.theVector))(i) = theMaterial1d[i]->getStress();
		}
		return 0;

	case 3:
		if (eleInformation.theVector != 0) {
			for (int i = 0; i < 3; i++) {
				(*(eleInformation.theVector))(i) = theMaterial1d[0]->getStrain();
				//opserr << "Moon_getResponse " << "  " << " *(eleInformation.theVector) " << *(eleInformation.theVector) << endln;
			}
		}
		return 0;

	case 13:
		if (eleInformation.theMatrix != 0) {
			for (int i = 0; i < 3; i++)
				(*(eleInformation.theMatrix))(i, i) = theMaterial1d[0]->getTangent();
			//(*(eleInformation.theMatrix))(i, i) = theMaterial1d[i]->getTangent();
		}
		return 0;

	case 4:
		if (eleInformation.theVector != 0) {
			for (int i = 0; i < 3; i++) {
				(*(eleInformation.theVector))(i) = theMaterial1d[0]->getStrain();
				(*(eleInformation.theVector))(i + 3) = theMaterial1d[0]->getStress();
			}
		}
		return 0;

	default:
		return -1;
	}
}

int
GMG_CMAZ2D::setParameter(const char **argv, int argc, Parameter &param)
{
	int result = -1;

	if (argc < 1)
		return -1;

	if (strcmp(argv[0], "material") == 0) {
		if (argc > 2) {
			int matNum = atoi(argv[1]);
			if (matNum >= 1 && matNum <= numMaterials1d)
				return theMaterial1d[matNum - 1]->setParameter(&argv[2], argc - 2, param);
		}
		else {
			return -1;
		}
	}

	for (int i = 0; i < 3; i++) {
		int res = theMaterial1d[0]->setParameter(argv, argc, param);
		if (res != -1) {
			result = res;
		}
	}
	return result;
}

const Vector &
GMG_CMAZ2D::getResistingForceSensitivity(int gradIndex)
{
	// Recompute strains to be safe
	this->update();

	double dfdh;

	// zero the residual
	theVector->Zero();

	// loop over 1d materials
	for (int mat = 0; mat < 3; mat++) {

		// get resisting force for material
		dfdh = theMaterial1d[0]->getStressSensitivity(gradIndex, true);

		// compute residual due to resisting force
		for (int i = 0; i < numDOF; i++)
			(*theVector)(i) += (*t1d)(mat, i) * dfdh;

	} // end loop over 1d materials 

	return *theVector;
}

int
GMG_CMAZ2D::commitSensitivity(int gradIndex, int numGrads)
{
	// Get nodal displacement sensitivity
	Vector diff(numDOF / 2);
	for (int i = 0; i < numDOF / 2; i++) {
		diff(i) = theNodes[1]->getDispSensitivity(i + 1, gradIndex) - theNodes[0]->getDispSensitivity(i + 1, gradIndex);
	}

	double depsdh;
	int ret = 0;
	for (int mat = 0; mat < 3; mat++) {
		// compute strain and rate; set as current trial for material
		//depsdh = this->computeCurrentStrain1d(mat, diff);
		ret += theMaterial1d[0]->commitSensitivity(depsdh, gradIndex, numGrads);
	}

	return ret;
}




// Private methods


// Establish the external nodes and set up the transformation matrix
// for orientation
void
GMG_CMAZ2D::setUp(int Nd1, int Nd2,
	const Vector &x,
	const Vector &yp)
{
	// ensure the connectedExternalNode ID is of correct size & set values
	if (connectedExternalNodes.Size() != 2)
		opserr << "FATAL ZeroLength::setUp - failed to create an ID of correct size\n";

	connectedExternalNodes(0) = Nd1;
	connectedExternalNodes(1) = Nd2;

	int i;
	for (i = 0; i < 2; i++)
		theNodes[i] = 0;

	// check that vectors for orientation are correct size
	if (x.Size() != 3 || yp.Size() != 3)
		opserr << "FATAL ZeroLength::setUp - incorrect dimension of orientation vectors\n";

	// establish orientation of element for the tranformation matrix
	// z = x cross yp
	Vector z(3);
	z(0) = x(1)*yp(2) - x(2)*yp(1);
	z(1) = x(2)*yp(0) - x(0)*yp(2);
	z(2) = x(0)*yp(1) - x(1)*yp(0);

	// y = z cross x
	Vector y(3);
	y(0) = z(1)*x(2) - z(2)*x(1);
	y(1) = z(2)*x(0) - z(0)*x(2);
	y(2) = z(0)*x(1) - z(1)*x(0);

	// compute length(norm) of vectors
	double xn = x.Norm();
	double yn = y.Norm();
	double zn = z.Norm();

	// check valid x and y vectors, i.e. not parallel and of zero length
	if (xn == 0 || yn == 0 || zn == 0) {
		opserr << "FATAL ZeroLength::setUp - invalid vectors to constructor\n";
	}

	// create transformation matrix of direction cosines
	for (i = 0; i < 3; i++) {
		transformation(0, i) = x(i) / xn;
		transformation(1, i) = y(i) / yn;
		transformation(2, i) = z(i) / zn;
	}

}


// Check that direction is in the range of 0 to 5
void
GMG_CMAZ2D::checkDirection(ID &dir) const
{
	//dir(0) = 1;
	for (int i = 0; i < dir.Size(); i++)
		if (dir(i) < 0 || dir(i) > 5) {
			//opserr << "WARNING ZeroLength::checkDirection - incorrect direction " << dir(i) << " is set to 0\n";
			//opserr << "Element " << dir.Size() << "  " << dir(i) << endln;
			dir(i) = 0;
		}
}


// Set basic deformation-displacement transformation matrix for 1d
// uniaxial materials
void
GMG_CMAZ2D::setTran1d(Etype elemType,
	int   numMat)
{
	enum Dtype { TRANS, ROTATE };

	int   indx, dir;
	Dtype dirType;

	// Create 1d transformation matrix
	//t1d = new Matrix(numMat, numDOF);
	t1d = new Matrix(3, 6);

	if (t1d == 0)
		opserr << "FATAL ZeroLength::setTran1d - can't allocate 1d transformation matrix\n";

	// Use reference for convenience and zero matrix.
	Matrix& tran = *t1d;
	tran.Zero();


	/*
	// loop over materials, setting row in tran for each material depending on dimensionality of element

	for (int i = 0; i < numMat; i++) {

		dir = (*dir1d)(i);	// direction 0 to 5;
		indx = dir % 3;		// direction 0, 1, 2 for axis of translation or rotation

		// set direction type to translation or rotation
		dirType = (dir < 3) ? TRANS : ROTATE;

		// now switch on dimensionality of element

		switch (elemType) {

		case D1N2:
			if (dirType == TRANS)
				tran(i, 1) = transformation(indx, 0);
			break;

		case D2N4:
			if (dirType == TRANS) {
				tran(i, 2) = transformation(indx, 0);
				tran(i, 3) = transformation(indx, 1);
			}
			break;

		case D2N6:
			if (dirType == TRANS) {
				tran(i, 3) = transformation(indx, 0);
				tran(i, 4) = transformation(indx, 1);
				tran(i, 5) = 0.0;
			}
			else if (dirType == ROTATE) {
				tran(i, 3) = 0.0;
				tran(i, 4) = 0.0;
				tran(i, 5) = transformation(indx, 2);
			}
			break;

		case D3N6:
			if (dirType == TRANS) {
				tran(i, 3) = transformation(indx, 0);
				tran(i, 4) = transformation(indx, 1);
				tran(i, 5) = transformation(indx, 2);
			}
			break;

		case D3N12:
			if (dirType == TRANS) {
				tran(i, 6) = transformation(indx, 0);
				tran(i, 7) = transformation(indx, 1);
				tran(i, 8) = transformation(indx, 2);
				tran(i, 9) = 0.0;
				tran(i, 10) = 0.0;
				tran(i, 11) = 0.0;
			}
			else if (dirType == ROTATE) {
				tran(i, 6) = 0.0;
				tran(i, 7) = 0.0;
				tran(i, 8) = 0.0;
				tran(i, 9) = transformation(indx, 0);
				tran(i, 10) = transformation(indx, 1);
				tran(i, 11) = transformation(indx, 2);
			}
			break;

		} // end switch

		// fill in first half of transformation matrix with
		// negative sign

		for (int j = 0; j < numDOF / 2; j++)
			tran(i, j) = -tran(i, j + numDOF / 2);

	} // end loop over 1d materials
	*/
	dir = 0;
	for (int i = 0; i <= 2; i++) {

		indx = dir % 3;		// direction 0, 1, 2 for axis of translation or rotation
		dirType = (dir < 3) ? TRANS : ROTATE;
		if (dir == 0 || dir == 1) {
			tran(i, 3) = transformation(indx, 0);
			tran(i, 4) = transformation(indx, 1);
			tran(i, 5) = 0.0;
		}
		else if (dir == 2) {
			tran(i, 3) = 0.0;
			tran(i, 4) = 0.0;
			tran(i, 5) = transformation(indx, 2);
		}

		// fill in first half of transformation matrix with
		// negative sign

		for (int j = 0; j < numDOF / 2; j++)
			tran(i, j) = -tran(i, j + numDOF / 2);

		dir++;

	}
	//opserr << "moon1 " << "  " << tran << endln;
	//tran(0, 0) = -1.0; tran(0, 1) = 0.0; tran(0, 2) = 0.0; tran(0, 3) = 1.0; tran(0, 4) = 0.0; tran(0, 5) = 0.0;
	//tran(1, 0) = 0.0; tran(1, 1) = -1.0; tran(1, 2) = 0.0; tran(1, 3) = 0.0; tran(1, 4) = 1.0; tran(1, 5) = 0.0;
	//tran(2, 0) = 0.0; tran(2, 1) = 0.0; tran(2, 2) = -1.0; tran(2, 3) = 0.0; tran(2, 4) = 0.0; tran(2, 5) = 1.0;
}


// Compute current strain for 1d material mat
// dispDiff are the displacements of node 2 minus those
// of node 1
Vector
GMG_CMAZ2D::computeCurrentStrain1d(int mat,
	const Vector& dispDiff) const
{

	/*
	double strain = 0.0;

	for (int i = 0; i < numDOF / 2; i++) {
		strain += -dispDiff(i) * (*t1d)(mat, i);
	}
	*/
	Vector strain_vector(3);
	strain_vector.Zero();
	//opserr << "Moon4_3" << "  " << " strain_vector " << strain_vector << "  " << " dispDiff " << dispDiff << "  " << " *t1d " << *t1d << endln;
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < numDOF / 2; i++)
			strain_vector(j) += -dispDiff(i) * (*t1d)(j, i);

	}
	//opserr << "Moon4_2" << "  " << " strain_vector " << strain_vector << "  " << " dispDiff " << dispDiff << "  " << " *t1d " << *t1d << endln;

	//strain_vector(0) += (-1) * dispDiff(0);
	//strain_vector(1) += (-1) * dispDiff(1);
	//strain_vector(2) += (-1) * dispDiff(2);

	return strain_vector;
}

void
GMG_CMAZ2D::updateDir(const Vector& x, const Vector& y)
{
	this->setUp(connectedExternalNodes(0), connectedExternalNodes(1), x, y);
	this->setTran1d(elemType, numMaterials1d);
}
