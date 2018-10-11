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
                                                                        
// $Revision: 1.7 $
// $Date: 2009-05-18 22:01:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthND.cpp,v $
                                                                        
// Written: MHS
// Created: Sept 2000
//
// Description: This file contains the implementation for the 
// ZeroLengthND class.

#include <ZeroLengthND.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <NDMaterial.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <G3Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>

Matrix ZeroLengthND::K6(6,6);
Matrix ZeroLengthND::K12(12,12);

Vector ZeroLengthND::P6(6);
Vector ZeroLengthND::P12(12);

Vector ZeroLengthND::v2(2);
Vector ZeroLengthND::v3(3);

//  Constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the ZeroLengthND end nodes.

static Node *theNodes[2];

void* OPS_ZeroLengthND()
{
    int ndm = OPS_GetNDM();
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 4) {
        opserr << "WARNING too few arguments " <<
            "want - element zeroLengthND eleTag? iNode? jNode? " <<
	    "NDTag? <1DTag?>" <<
	    "<-orient x1? x2? x3? y1? y2? y3?>\n";

        return 0;
    }

    int idata [4];
    numdata = 4;
    if (OPS_GetIntInput(&numdata,idata) < 0) {
        opserr << "WARNING: failed to get integer data\n";
        return 0;
    }
    NDMaterial* nmat = OPS_getNDMaterial(idata[3]);
    if (nmat == 0) {
	opserr<<"WARNING: NDMaterial "<<idata[3]<<" is not defined\n";
	return 0;
    }

    UniaxialMaterial* umat = 0;
    int uniTag;
    if (OPS_GetIntInput(&numdata,&uniTag) >= 0) {
	umat = OPS_getUniaxialMaterial(uniTag);
	if (umat == 0) {
	    opserr<<"WARNING: uniaxial material "<<uniTag<<" is not defined\n";
	    return 0;
	}
    } else {
	OPS_ResetCurrentInputArg(-1);
    }

    const char* type = OPS_GetString();
    Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
    Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
    if (strcmp(type,"-orient") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 6) {
	    opserr<<"WARNING: insufficient orient values\n";
	    return 0;
	}
	numdata = 3;
	if (OPS_GetDoubleInput(&numdata,&x(0)) < 0) {
	    opserr<<"WARNING: invalid double input\n";
	    return 0;
	}
	if (OPS_GetDoubleInput(&numdata,&y(0)) < 0) {
	    opserr<<"WARNING: invalid double input\n";
	    return 0;
	}
    }

    if (umat == 0) {
	return new ZeroLengthND(idata[0],ndm,idata[1],idata[2],x,y,*nmat);
    } else {
	return new ZeroLengthND(idata[0],ndm,idata[1],idata[2],x,y,*nmat,*umat);
    }
}


ZeroLengthND::ZeroLengthND(int tag, int dim, int Nd1, int Nd2, 
	       const Vector& x, const Vector& yprime, 
		   NDMaterial &theNDmat) : 
Element(tag, ELE_TAG_ZeroLengthND),
connectedExternalNodes(2),
dimension(dim), numDOF(0), 
transformation(3,3), A(0), v(0), e(0.0), K(0), P(0),
end1Ptr(0), end2Ptr(0), theNDMaterial(0), the1DMaterial(0), order(0)
{
	// Obtain copy of Nd material model
	theNDMaterial = theNDmat.getCopy();
	
	if (theNDMaterial == 0) {
		opserr << "ZeroLengthND::zeroLengthND-- failed to get copy of NDMaterial\n";
		exit(-1);
	}
	// Get the material order
	order = theNDMaterial->getOrder();

	// Check material order
	if (order < 2 || order > 3) {
		opserr << "ZeroLengthND::  -- NDMaterial not of order 2 or 3\n";
		exit(-1);
	}

	// Set up the transformation matrix of direction cosines
	this->setUp(Nd1, Nd2, x, yprime);
}

ZeroLengthND::ZeroLengthND(int tag, int dim, int Nd1, int Nd2, 
	       const Vector& x, const Vector& yprime, 
		   NDMaterial &theNDmat, UniaxialMaterial &the1Dmat) : 
Element(tag, ELE_TAG_ZeroLengthND),
connectedExternalNodes(2),
dimension(dim), numDOF(0), 
transformation(3,3), A(0), v(0), e(0.0), K(0), P(0),
end1Ptr(0), end2Ptr(0), theNDMaterial(0), the1DMaterial(0), order(0)
{
	// Obtain copy of Nd material model
	theNDMaterial = theNDmat.getCopy();
	
	if (theNDMaterial == 0) {
		opserr << "ZeroLengthND::  -- failed to get copy of NDMaterial\n";
		exit(-1);
	}

	// Obtain copy of 1d material model
	the1DMaterial = the1Dmat.getCopy();
	
	if (the1DMaterial == 0) {
		opserr << "ZeroLengthND""ZeroLengthND -- failed to get copy of UniaxialMaterial\n";
		exit(-1);	
	}	

	// Get the material order
	order = theNDMaterial->getOrder();

	if (order != 2) {
		opserr << "ZeroLengthND::ZeroLengthND-- NDMaterial not of order 2\n";
		exit(-1);
	}

	// Set up the transformation matrix of direction cosines
	this->setUp(Nd1, Nd2, x, yprime);
}

ZeroLengthND::ZeroLengthND() : 
Element(0, ELE_TAG_ZeroLengthND),
connectedExternalNodes(2),
dimension(0), numDOF(0), 
transformation(3,3), A(0), v(0), e(0.0), K(0), P(0),
end1Ptr(0), end2Ptr(0), theNDMaterial(0), the1DMaterial(0), order(0)
{

}

ZeroLengthND::~ZeroLengthND()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    
	if (theNDMaterial != 0)
		delete theNDMaterial;
	if (the1DMaterial != 0)
		delete the1DMaterial;
	if (A != 0)
		delete A;
}

int
ZeroLengthND::getNumExternalNodes(void) const
{
    return 2;
}

Node **
ZeroLengthND::getNodePtrs(void) 
{
  theNodes[0] = end1Ptr;
  theNodes[1] = end2Ptr;
  return theNodes;
}

const ID &
ZeroLengthND::getExternalNodes(void) 
{
  return connectedExternalNodes;
}



int
ZeroLengthND::getNumDOF(void) 
{
    return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the ZeroLengthND element, we set matrix and vector pointers,
//    allocate space for t matrix and define it as the basic deformation-
//    displacement transformation matrix.
void
ZeroLengthND::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
		end1Ptr = 0;
		end2Ptr = 0;
		return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);	

    // if can't find both - send a warning message
    if (end1Ptr == 0 || end2Ptr == 0) {
		if (end1Ptr == 0) 
			opserr << "ZeroLengthND::setDomain()-- Nd1 does not exist in model\n";
		else
			opserr << "ZeroLengthND::setDomain -- Nd2 does not exist in model\n";

		end1Ptr = 0;
		end2Ptr = 0;
		return;
    }

    // now determine the number of dof and the dimension    
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2) {
		opserr << "ZeroLengthND::setDomain -- nodes have differing dof's at end\n";
		end1Ptr = 0;
		end2Ptr = 0;
		return;
    }	

	numDOF = 2*dofNd1;

	if (numDOF != 6 && numDOF != 12) {
		opserr << "ZeroLengthND::setDomain -  element only works for 3 (2d) or 6 (3d) dof per node\n";
		end1Ptr = 0;
		end2Ptr = 0;
		return;
	}

    // Check that length is zero within tolerance
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	
    const Vector     diff = end1Crd - end2Crd;
    double L  = diff.Norm();
    double v1 = end1Crd.Norm();
    double v2 = end2Crd.Norm();
    double vm;
    
    vm = (v1<v2) ? v2 : v1;
    
    if (L > LENTOL*vm)
		opserr << "ZeroLengthND::setDomain -- Element has L=, which is greater than the tolerance\n";
        
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
	// Set up the A matrix
	this->setTransformation();
}   	 

int
ZeroLengthND::commitState()
{
	int err = 0;

	// Commit the NDMaterial
	err += theNDMaterial->commitState();

	// Commit the UniaxialMaterial
	if (the1DMaterial != 0)
		err += the1DMaterial->commitState();

	return err;
}

int
ZeroLengthND::revertToLastCommit()
{
	int err = 0;

	// Revert the NDMaterial
	err += theNDMaterial->revertToLastCommit();

	// Revert the UniaxialMaterial
	if (the1DMaterial != 0)
		err += the1DMaterial->revertToLastCommit();

	return err;
}

int
ZeroLengthND::revertToStart()
{
	int err = 0;

	// Revert the NDMaterial to start
	err += theNDMaterial->revertToStart();

	// Revert the UniaxialMaterial to start
	if (the1DMaterial != 0)
		err += the1DMaterial->revertToStart();

	return err;
}

const Matrix &
ZeroLengthND::getTangentStiff(void)
{
	// Compute material strains
	this->computeStrain();

	// Set trial strain for NDMaterial
	theNDMaterial->setTrialStrain(*v);

	// Get NDMaterial tangent, the element basic stiffness
	const Matrix &kb = theNDMaterial->getTangent();

	// Set some references to make the syntax nicer
	Matrix &stiff = *K;
	const Matrix &tran = *A;

	stiff.Zero();

	double E;

	// Compute element stiffness ... K = A^*kb*A
	for (int k = 0; k < order; k++) {
		for (int l = 0; l < order; l++) {
			E = kb(k,l);
			for (int i = 0; i < numDOF; i++)
				for (int j = 0; j < i+1; j++)
					stiff(i,j) +=  tran(k,i) * E * tran(l,j);
		}
	}

	if (the1DMaterial != 0) {

		// Set trial strain for UniaxialMaterial
		the1DMaterial->setTrialStrain(e);

		// Get UniaxialMaterial tangent, the element basic stiffness
		E = the1DMaterial->getTangent();

		// Compute element stiffness ... K = A^*kb*A
		for (int i = 0; i < numDOF; i++)
			for (int j = 0; j < i+1; j++)
				stiff(i,j) +=  tran(2,i) * E * tran(2,j);
	}

    // Complete symmetric stiffness matrix
    for (int i = 0; i < numDOF; i++)
		for(int j = 0; j < i; j++)
		    stiff(j,i) = stiff(i,j);

	return stiff;
}

const Matrix &
ZeroLengthND::getInitialStiff(void)
{
  // Get NDMaterial tangent, the element basic stiffness
  const Matrix &kb = theNDMaterial->getInitialTangent();

  // Set some references to make the syntax nicer
  Matrix &stiff = *K;
  const Matrix &tran = *A;
  
  stiff.Zero();
  
  double E;
  
  // Compute element stiffness ... K = A^*kb*A
  for (int k = 0; k < order; k++) {
    for (int l = 0; l < order; l++) {
      E = kb(k,l);
      for (int i = 0; i < numDOF; i++)
	for (int j = 0; j < i+1; j++)
	  stiff(i,j) +=  tran(k,i) * E * tran(l,j);
    }
  }
  
  if (the1DMaterial != 0) {
    
    // Get UniaxialMaterial tangent, the element basic stiffness
    E = the1DMaterial->getInitialTangent();
    
    // Compute element stiffness ... K = A^*kb*A
    for (int i = 0; i < numDOF; i++)
      for (int j = 0; j < i+1; j++)
	stiff(i,j) +=  tran(2,i) * E * tran(2,j);
  }

  // Complete symmetric stiffness matrix
  for (int i = 0; i < numDOF; i++)
    for(int j = 0; j < i; j++)
      stiff(j,i) = stiff(i,j);
  
  return stiff;
}

const Matrix &
ZeroLengthND::getDamp(void)
{
	// Return zero damping
	K->Zero();

	return *K;
}

const Matrix &
ZeroLengthND::getMass(void)
{
	// Return zero mass
	K->Zero();

	return *K;
}

void 
ZeroLengthND::zeroLoad(void)
{
	// does nothing now
}

int 
ZeroLengthND::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "ZeroLength::addLoad - load type unknown for ZeroLengthND\n";
  
  return -1;
}

int 
ZeroLengthND::addInertiaLoadToUnbalance(const Vector &accel)
{
	// does nothing as element has no mass yet!
	return 0;
}

const Vector &
ZeroLengthND::getResistingForce()
{
	// Compute material strains
	this->computeStrain();

	// Set trial strain for NDMaterial
	theNDMaterial->setTrialStrain(*v);

	// Get NDMaterial stress, the element basic force
	const Vector &q = theNDMaterial->getStress();

	// Set some references to make the syntax nicer
	Vector &force = *P;
	const Matrix &tran = *A;

	force.Zero();

	double s;

	// Compute element resisting force ... P = A^*q
	for (int k = 0; k < order; k++) {
		s = q(k);
		for (int i = 0; i < numDOF; i++)
			force(i) += tran(k,i) * s;
	}

	if (the1DMaterial != 0) {

		// Set trial strain for UniaxialMaterial
		the1DMaterial->setTrialStrain(e);

		// Get UniaxialMaterial stress, the element basic force
		s = the1DMaterial->getStress();

		// Compute element resisting force ... P = A^*q
		for (int i = 0; i < numDOF; i++)
			force(i) += tran(2,i) * s;
	}

	return force;
}

const Vector &
ZeroLengthND::getResistingForceIncInertia()
{	
    // There is no mass, so return
    return this->getResistingForce();
}

int
ZeroLengthND::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// ZeroLengthND packs its data into an ID and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments

	static ID idData(11);

	idData(0) = this->getTag();
	idData(1) = dimension;
	idData(2) = numDOF;
	idData(3) = order;
	idData(4) = (the1DMaterial == 0) ? 0 : 1;
	idData(5) = connectedExternalNodes(0);
	idData(6) = connectedExternalNodes(1);
	idData(7) = theNDMaterial->getClassTag();

	int dbTag = theNDMaterial->getDbTag();
	if (dbTag == 0) {
		dbTag = theChannel.getDbTag();
		if (dbTag != 0)
			theNDMaterial->setDbTag(dbTag);
	}
	idData(8) = dbTag;

	if (the1DMaterial != 0) {
		idData(9) = the1DMaterial->getClassTag();

		dbTag = the1DMaterial->getDbTag();
		if (dbTag == 0) {
			dbTag = theChannel.getDbTag();
			if (dbTag != 0)
				the1DMaterial->setDbTag(dbTag);
		}
		idData(10) = dbTag;
	}

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "ZeroLengthND::sendSelf() -- failed to send ID data\n";
		return res;
	}

	// Send the 3x3 direction cosine matrix, have to send it since it is only set
	// in the constructor and not setDomain()
	res += theChannel.sendMatrix(dataTag, commitTag, transformation);
	if (res < 0) {
		opserr << "ZeroLengthND::sendSelf -- failed to send transformation Matrix\n";
		return res;
	}

	// Send the NDMaterial
	res += theNDMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr << "ZeroLengthND::  -- failed to send NDMaterial\n";
		return res;
	}

	// Send the UniaxialMaterial, if present
	if (the1DMaterial != 0) {
		res += the1DMaterial->sendSelf(commitTag, theChannel);
		if (res < 0) {
			opserr << "ZeroLengthND::sendSelf-- failed to send UniaxialMaterial";
			return res;
		}
	}

	return res;
}

int
ZeroLengthND::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
  
	int dataTag = this->getDbTag();

	// ZeroLengthND creates an ID, receives the ID and then sets the 
	// internal data with the data in the ID

	static ID idData(11);

	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "ZeroLengtHND::recvSelf -- failed to receive ID data\n";
		return res;
	}

	res += theChannel.recvMatrix(dataTag, commitTag, transformation);
	if (res < 0) {
		opserr << "zeroLengthND::revbSelf -- failed to receive transformation Matrix\n";
		return res;
	}

	this->setTag(idData(0));
	dimension = idData(1);
	numDOF = idData(2);
	connectedExternalNodes(0) = idData(5);
	connectedExternalNodes(1) = idData(6);

	if (order != idData(3)) {

		order = idData(3);

		// Allocate transformation matrix
		if (A != 0)
			delete A;

		A = new Matrix(order, numDOF);


		if (numDOF == 6) {
			K = &K6;
			P = &P6;
		}
		else {
			K = &K12;
			P = &P12;
		}

		if (order == 2)
			v = &v2;
		else
			v = &v3;
	}

	int classTag = idData(7);

	// If null, get a new one from the broker
	if (theNDMaterial == 0)
		theNDMaterial = theBroker.getNewNDMaterial(classTag);

	// If wrong type, get a new one from the broker
	if (theNDMaterial->getClassTag() != classTag) {
		delete theNDMaterial;
		theNDMaterial = theBroker.getNewNDMaterial(classTag);
	}

	// Check if either allocation failed from broker
	if (theNDMaterial == 0) {
		opserr << "ZeroLengthND::  -- failed to allocate new NDMaterial\n";
		return -1;
	}

	// Receive the NDMaterial
	theNDMaterial->setDbTag(idData(8));
	res += theNDMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		opserr << "ZeroLengthND::  -- failed to receive NDMaterial\n";
		return res;
	}

	// Receive the UniaxialMaterial, if present
	if (idData(4) == 1) {
		classTag = idData(9);

		// If null, get a new one from the broker
		if (the1DMaterial == 0)
			the1DMaterial = theBroker.getNewUniaxialMaterial(classTag);

		// If wrong type, get a new one from the broker
		if (the1DMaterial->getClassTag() != classTag) {
			delete the1DMaterial;
			the1DMaterial = theBroker.getNewUniaxialMaterial(classTag);
		}

		// Check if either allocation failed from broker
		if (the1DMaterial == 0) {
			opserr << "ZeroLengthND::  -- failed to allocate new UniaxialMaterial\n";
			return -1;
		}

		// Receive the UniaxialMaterial
		the1DMaterial->setDbTag(idData(10));
		res += the1DMaterial->recvSelf(commitTag, theChannel, theBroker);
		if (res < 0) {
			opserr << "ZeroLengthND::  -- failed to receive UniaxialMaterial\n";
			return res;
		}
	}

	return res;
}

int
ZeroLengthND::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // ensure setDomain() worked
    if (end1Ptr == 0 || end2Ptr == 0)
		return 0;

    // first determine the two end points of the ZeroLengthND based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	
    const Vector &end1Disp = end1Ptr->getDisp();
    const Vector &end2Disp = end2Ptr->getDisp();    

    if (displayMode == 1 || displayMode == 2) {
		static Vector v1(3);
		static Vector v2(3);
		for (int i = 0; i < dimension; i++) {
			v1(i) = end1Crd(i)+end1Disp(i)*fact;
			v2(i) = end2Crd(i)+end2Disp(i)*fact;    
		}
		
		return theViewer.drawLine(v1, v2, 0.0, 0.0);
    }

    return 0;
}

void
ZeroLengthND::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "ZeroLengthND, tag: " << this->getTag() << endln;
        s << "\tConnected Nodes: " << connectedExternalNodes << endln;
        s << "\tNDMaterial, tag: " << theNDMaterial->getTag() << endln;
        if (the1DMaterial != 0)
            s << "\tUniaxialMaterial, tag: " << the1DMaterial->getTag() << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ZeroLengthND\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"ndMaterial\": \"" << theNDMaterial->getTag() << "\", ";
        if (the1DMaterial != 0)
            s << "\"uniaxialMaterial\": \"" << the1DMaterial->getTag() << "\", ";
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
ZeroLengthND::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","ZeroLength");
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
            theResponse =  new ElementResponse(this, 1, *P);

    } else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

        if (the1DMaterial != 0) {
            for (int i=0; i<3; i++) {
                sprintf(outputData,"P%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 2, Vector(3));
        } else {
            for (int i=0; i<order; i++) {
                sprintf(outputData,"P%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 2, Vector(order));
        }

    } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
        strcmp(argv[0],"deformation") == 0) {

            if (the1DMaterial != 0) {
                for (int i=0; i<3; i++) {
                    sprintf(outputData,"e%d",i+1);
                    output.tag("ResponseType",outputData);
                }
                theResponse = new ElementResponse(this, 3, Vector(3));
            } else {
                for (int i=0; i<order; i++) {
                    sprintf(outputData,"e%d",i+1);
                    output.tag("ResponseType",outputData);
                }
                theResponse = new ElementResponse(this, 3, Vector(order));
            }

    // a material quantity
    } else if (strcmp(argv[0],"material") == 0) {
        // See if NDMaterial can handle request ...
        theResponse = theNDMaterial->setResponse(&argv[1], argc-1, output);

        if ((theResponse == 0) && (the1DMaterial != 0))
            theResponse = the1DMaterial->setResponse(&argv[1], argc-1, output);
    }

    output.endTag();
    return theResponse;
}

int 
ZeroLengthND::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {
    case 1:
        return eleInfo.setVector(this->getResistingForce());

    case 2:
        if (eleInfo.theVector != 0) {
            const Vector &tmp = theNDMaterial->getStress();
            Vector &force = *(eleInfo.theVector);
            for (int i = 0; i < order; i++)
                force(i) = tmp(i);
            if (the1DMaterial != 0)
                force(order) = the1DMaterial->getStress();
        }
        return 0;

    case 3:
        if (eleInfo.theVector != 0) {
            this->computeStrain();
            const Vector &tmp = *v;	// NDMaterial strains
            Vector &def = *(eleInfo.theVector);
            for (int i = 0; i < order; i++)
                def(i) = tmp(i);
            if (the1DMaterial != 0)
                def(order) = e;	// UniaxialMaterial strain
        }
        return 0;

    default:
        return -1;
    }
}

// Private methods


// Establish the external nodes and set up the transformation matrix
// for orientation
void
ZeroLengthND::setUp(int Nd1, int Nd2, const Vector &x, const Vector &yp)
{ 
    // ensure the connectedExternalNode ID is of correct size & set values
	if (connectedExternalNodes.Size() != 2) {
		opserr << "ZeroLengthND::setUp -- failed to create an ID of correct size\n";
		exit(-1);
	}
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    // check that vectors for orientation are correct size
	if ( x.Size() != 3 || yp.Size() != 3 ) {
		opserr << "ZeroLengthND -- incorrect dimension of orientation vectors\n";
	exit(-1);
	}
    // establish orientation of element for the transformation matrix
    // z = x cross yp
    static Vector z(3);
    z(0) = x(1)*yp(2) - x(2)*yp(1);
    z(1) = x(2)*yp(0) - x(0)*yp(2);
    z(2) = x(0)*yp(1) - x(1)*yp(0);

    // y = z cross x
    static Vector y(3);
    y(0) = z(1)*x(2) - z(2)*x(1);
    y(1) = z(2)*x(0) - z(0)*x(2);
    y(2) = z(0)*x(1) - z(1)*x(0);

    // compute length(norm) of vectors
    double xn = x.Norm();
    double yn = y.Norm();
    double zn = z.Norm();

    // check valid x and y vectors, i.e. not parallel and of zero length
	if (xn == 0 || yn == 0 || zn == 0){
		opserr << "ZeroLengthND::setUP -- invalid vectors to constructor\n";
		exit(-1);
	}
    
    // create transformation matrix of direction cosines
    for (int i = 0; i < 3; i++) {
		transformation(0,i) = x(i)/xn;
		transformation(1,i) = y(i)/yn;
		transformation(2,i) = z(i)/zn;
	}
}

// Set basic deformation-displacement transformation matrix for the materials
void 
ZeroLengthND::setTransformation(void)
{
	// Allocate transformation matrix
	if (A != 0)
		delete A;

	A = (the1DMaterial == 0) ? new Matrix(order, numDOF) : new Matrix(order+1, numDOF);

	if (A == 0) {
		opserr << "ZeroLengthND::setTransformation -- failed to allocate transformation Matrix\n";
			exit(-1);
	}
	if (numDOF == 6) {
		K = &K6;
		P = &P6;
	}
	else {
		K = &K12;
		P = &P12;
	}

	if (order == 2)
		v = &v2;
	else
		v = &v3;

	// Set a reference to make the syntax nicer
	Matrix &tran = *A;
	
	// Loop over the NDMaterial order
	for (int i = 0; i < order; i++) {

		if (numDOF == 6) {
			tran(i,3) = transformation(i,0);
			tran(i,4) = transformation(i,1);
			tran(i,5) = 0.0;
		}
		else if (numDOF == 12) {
			tran(i,6) = transformation(i,0);
			tran(i,7) = transformation(i,1);
			tran(i,8) = transformation(i,2);
		}

		// Fill in first half of transformation matrix with negative sign
		for (int j = 0; j < numDOF/2; j++ )
			tran(i,j) = -tran(i,j+numDOF/2);
	}

	// Fill in transformation for UniaxialMaterial
	if (the1DMaterial != 0) {

		if (numDOF == 6) {
			tran(2,3) = transformation(2,0);
			tran(2,4) = transformation(2,1);
			tran(2,5) = 0.0;
		}
		else if (numDOF == 12) {
			tran(2,6) = transformation(2,0);
			tran(2,7) = transformation(2,1);
			tran(2,8) = transformation(2,2);
		}

		// Fill in first half of transformation matrix with negative sign
		for (int j = 0; j < numDOF/2; j++ )
			tran(2,j) = -tran(2,j+numDOF/2);
	}
}
		     
void
ZeroLengthND::computeStrain(void)
{
	// Get nodal displacements
	const Vector &u1 = end1Ptr->getTrialDisp();
	const Vector &u2 = end2Ptr->getTrialDisp();

	// Compute differential displacements
	const Vector diff = u2 - u1;

	// Set some references to make the syntax nicer
	Vector &def = *v;
	const Matrix &tran = *A;

	def.Zero();

	// Compute element basic deformations ... v = A*(u2-u1)
	for (int i = 0; i < order; i++)
		for (int j = 0; j < numDOF/2; j++)
			def(i) += -diff(j)*tran(i,j);

	if (the1DMaterial != 0) {
		e = 0.0;
		for (int j = 0; j < numDOF/2; j++)
			e += -diff(j)*tran(2,j);
	}
}
