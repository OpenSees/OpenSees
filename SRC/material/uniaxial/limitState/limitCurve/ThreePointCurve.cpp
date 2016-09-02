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

// $Revision: 1.3 $
// $Date: 2007-02-02 01:19:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/limitCurve/ThreePointCurve.cpp,v $
// Written: KJE
// Created: Aug 2001
// Modified: Jul 2002


#include <ThreePointCurve.h>
#include <G3Globals.h>
#include <math.h>
#include <ElementResponse.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <Vector.h>
#include <float.h>

#include <DummyStream.h>
#include <elementAPI.h>

void* OPS_ThreePointCurve()
{
    if (OPS_GetNumRemainingInputArgs() < 12) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: limitCurve ThreePoint tag? eleTag? x1? y1? x2? y2? x3? y3?";
	opserr << "Kdeg? Fres? defType? forType?" << endln;
	opserr << "<ndI? ndJ? dof? perpDirn?>" << endln;
	return 0;
    }
    int tag;
    int eleTag;
    double Kdeg;
    double Fres;
    int defType, forType;
    double x1, y1;
    double x2, y2;
    double x3, y3;
    int ndI = 0;
    int ndJ = 0;
    int dof = 0;
    int perpDirn = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid limitCurve ThreePoint tag" << endln;
	return 0;
    }
    if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
	opserr << "WARNING invalid element tag for associated beam-column element (eleTag)\n";
	opserr << "LimitCurve ThreePoint: " << tag << endln;
	return 0;
    }

    numdata = 8;
    double data[8];
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double data\n";
	opserr << "limitCurve ThreePoint: " << tag << endln;
	return 0;
    }
    x1 = data[0];
    y1 = data[1];
    x2 = data[2];
    y2 = data[3];
    x3 = data[4];
    y3 = data[5];
    Kdeg = data[6];
    Fres = data[7];

    numdata = 1;
    if (OPS_GetIntInput(&numdata, &defType) < 0) {
	opserr << "WARNING invalid deformation type defType\n";
	opserr << "LimitCurve ThreePoint: " << tag << endln;
	return 0;
    }
    if (OPS_GetIntInput(&numdata, &forType) < 0) {
	opserr << "WARNING invalid force type forType\n";
	opserr << "LimitCurve ThreePoint: " << tag << endln;
	return 0;
    }
    if (defType == 2) {

	if (OPS_GetNumRemainingInputArgs() < 4) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: limitCurve ThreePoint tag? eleTag? x1? y1? x2? y2? x3? y3?";
	    opserr << "Kdeg? Fres? defType? forType?" << endln;
	    opserr << "ndI? ndJ? dof? perpDirn?" << endln;
	}
	if (OPS_GetIntInput(&numdata, &ndI) < 0) {
	    opserr << "WARNING invalid node I\n";
	    opserr << "LimitCurve ThreePoint: " << tag << endln;
	    return 0;
	}
	if (OPS_GetIntInput(&numdata, &ndJ) < 0) {
	    opserr << "WARNING invalid node J\n";
	    opserr << "LimitCurve ThreePoint: " << tag << endln;
	    return 0;
	}
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING invalid degree of freedom for drift\n";
	    opserr << "LimitCurve ThreePoint: " << tag << endln;
	    return 0;
	}
	if (OPS_GetIntInput(&numdata, &perpDirn) < 0) {
	    opserr << "WARNING invalid direction for column length\n";
	    opserr << "LimitCurve ThreePoint: " << tag << endln;
	    return 0;
	}
    }
    // Parsing was successful, allocate the material
    // Subtract one from dof and perpDirn for C indexing
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return 0;
    return new ThreePointCurve(tag, eleTag, theDomain,
			       x1, y1, x2, y2, x3, y3, Kdeg, Fres, defType, forType,
			       ndI, ndJ, dof-1, perpDirn-1);
}

ThreePointCurve::ThreePointCurve(int tag, int eTag, Domain *theDom,
			double a1, double b1, double a2, double b2,
			double a3, double b3, double Kd, double Fr,
			int dType, int fType,
			int ni, int nj, int df, int dirn):
LimitCurve(tag, TAG_ThreePointCurve), eleTag(eTag), theDomain(theDom), theElement(0),
Kdeg(Kd), Fres(Fr), defType(dType), forType(fType),
x1(a1), y1(b1), x2(a2), y2(b2), x3(a3), y3(b3),
ndI(ni), ndJ(nj), dof(df), perpDirn(dirn)
{
	stateFlag = 0;
	count = 0;
//	this->Print(cout); // Commented out by Terje
}

ThreePointCurve::ThreePointCurve():
LimitCurve(0, TAG_ThreePointCurve), eleTag(0), theDomain(0), theElement(0),
Kdeg(0), Fres(0), defType(0), forType(0),
x1(0), y1(0), x2(0), y2(0), x3(0), y3(0)
{
//VOID//
}

ThreePointCurve::~ThreePointCurve()
{
//VOID//
}


LimitCurve*
ThreePointCurve::getCopy(void)
{
	ThreePointCurve *theCopy = new ThreePointCurve(this->getTag(),
		eleTag, theDomain, x1, y1, x2, y2, x3, y3,
		Kdeg, Fres, defType, forType, ndI, ndJ, dof, perpDirn);

	theCopy->stateFlag = stateFlag;

	return theCopy;
}


// check if limit state surface has been reached
int
ThreePointCurve::checkElementState(double springForce)
{
  DummyStream dummy;
	// find associated beam-column elementon first visit
	if (theElement == 0)
	{
		theElement = theDomain->getElement(eleTag);

		if (theElement == 0) {
//			g3ErrorHandler->fatal("WARNING ThreePointCurve - no element with tag %i exists in Domain",eleTag);
		}
		// find length between nodes if drift is desired
		if (defType == 2)
		{
			Node *nodeI = theDomain->getNode(ndI);
			Node *nodeJ = theDomain->getNode(ndJ);

			const Vector &crdI = nodeI->getCrds();
			const Vector &crdJ = nodeJ->getCrds();

			if (crdI(perpDirn) == crdJ(perpDirn)) {
//				g3ErrorHandler->warning("%s -- Nodal projection has zero component along chosen direction",
//						"AxialCurve::AxialCurve");

				oneOverL = 0.0;
			}
			else
				oneOverL = 1.0/fabs(crdJ(perpDirn) - crdI(perpDirn));
		}
	}

	double deform;	// value of deformation parameter from element
	double force;	// value of force parameter from element
	int result; //junk variable


	// Based on "defType" and "forType" calculate
	// the desired response parameters "deform" and "force"
	if (defType == 1) // maximum chord rotations
	{

		Response *theRotations =0; // integer element returns in setResponse
		//char *r[1] = {"rotation"};
		//char *r[1] = {"plasticRotation"};
		const char *r[1] = {"basicDeformation"};

		Vector *rotVec; //vector of chord rotations at beam-column ends

		// set type of beam-column element response desired
		theRotations = theElement->setResponse(r, 1, dummy);

		if (theRotations == 0) {
		  opserr << "ThreePointCurve::checkElementState, defType = 1, basicDeformations not implemented in element setResponse" << endln;
		  return -1;
		}

		// put element response in the vector of "myInfo"
		result = theRotations->getResponse();

		// access the myInfo vector containing the response (new for Version 1.2)
		Information &theInfo = theRotations->getInformation();
		rotVec = (theInfo.theVector);

		deform = (fabs((*rotVec)(1)) > fabs((*rotVec)(2))) ?
			fabs((*rotVec)(1)) : fabs((*rotVec)(2));  //use larger of two end rotations
	}
	else if (defType == 2) // interstory drift
	{
		// find associated nodes
		Node *nodeI = theDomain->getNode(ndI);
		Node *nodeJ = theDomain->getNode(ndJ);

		// get displacements
		const Vector &dispI = nodeI->getTrialDisp();
		const Vector &dispJ = nodeJ->getTrialDisp();

		// calc drift
		double dx = fabs(dispJ(dof)-dispI(dof));
		deform = dx*oneOverL;
	}
	else {
//		g3ErrorHandler->fatal("WARNING ThreePointCurve - deformation type flag %i not implemented",defType);
	}

		Response *theForces =0;
		const char *f[1] = {"localForce"};

		Vector *forceVec; //vector of basic forces from beam column

		// set type of beam-column element response desired
		theForces    = theElement->setResponse(f, 1, dummy);

		// put element response in the vector of "myInfo"
		result += theForces->getResponse();

		// access the myInfo vector containing the response (new for Version 1.2)
		Information &theInfo = theForces->getInformation();
		forceVec = (theInfo.theVector);

	// Local forces (assuming no element loads)
	if (forType == 0)
		force = fabs(springForce); //force in associated hysteretic material
	else if (forType == 1)
		force = fabs((*forceVec)(1)); //shear
	else if (forType == 2) //axial
		force = fabs((*forceVec)(0));
	else {
//		g3ErrorHandler->fatal("WARNING ThreePointCurve - force type flag %i not implemented",forType);
	}

	// Determine if (deform,force) is outside limit state surface.
	//
	// Use absolute value of deform and force
	// In future will include one positive and one negative limit state surfaces
	double forceSurface = findLimit(deform); // force on surface at deform

	count += 1;

	if (stateFlag == 0) //prior to failure
	{
		if (force >= forceSurface) // on/outside failure surface
		{
			stateFlag = 1;
			//Pshear = fabs((*forceVec)(0)); // axial load at shear failure (not currently used)
//			g3ErrorHandler->warning("ThreePointCurve - failure detected at deform = %f", deform);
		}
		else // inside failure surface
		{
			stateFlag = 0;
		}
	}
	else //after failure
	{
		if (force >= forceSurface) // on/outside failure surface
		{
			stateFlag = 2;
//			g3ErrorHandler->warning("WARNING ThreePointCurve - response past limit surface after failure at deform=%f", deform);
		}
		else // inside failure surface
		{
			stateFlag = 3;
		}
	}

	return stateFlag;
}


double
ThreePointCurve::getDegSlope(void)
{
	return Kdeg;
}


double
ThreePointCurve::getResForce(void)
{
	return Fres;
}

double
ThreePointCurve::getUnbalanceForce(void)
{
	//Do nothing for this class
	return 0.0;
}

int
ThreePointCurve::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int
ThreePointCurve::recvSelf(int commitTag, Channel &theChannel,
			FEM_ObjectBroker &theBroker)
{
	return -1;
}

void
ThreePointCurve::Print(OPS_Stream &s, int flag)
{
	s << "Three-Point Limit Curve, tag: " << this->getTag() << endln;
	s << "x1,y1: " << x1 <<", "<< y1 << endln;
	s << "x2,y2: " << x2 <<", "<< y2 << endln;
	s << "x3,y3: " << x3 <<", "<< y3 << endln;
	s << "eleTag: " << eleTag << endln;
	s << "nodeI: " << ndI << endln;
	s << "nodeJ: " << ndJ << endln;
	s << "deform: " << defType << endln;
	s << "force: " << forType << endln;
}



// Private Functions

double
ThreePointCurve::findLimit(double x)
{
	double y = 0.0;

	if (x < x1) {
//		g3ErrorHandler->fatal("Outside limits of ThreePointCurve");
	}
	else if (x < x2)
		y = y1+(y2-y1)/(x2-x1)*(x-x1);
	else if (x < x3)
		y = y2+(y3-y2)/(x3-x2)*(x-x2);
	else
		y = y3;

	return y;
}



int ThreePointCurve::revertToStart ()
{
	return 0;
}
