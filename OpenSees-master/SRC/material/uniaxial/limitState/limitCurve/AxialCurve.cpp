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

// $Revision: 1.4 $
// $Date: 2007-02-02 01:19:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/limitCurve/AxialCurve.cpp,v $
                                                                        
// Written: KJE
// Created: Aug 2002
//

#include <AxialCurve.h>
#include <G3Globals.h>
#include <math.h>
#include <ElementResponse.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <Vector.h>
#include <float.h>
#include <tcl.h>

#include <DummyStream.h>
#define PI (3.141592653589793)
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

#include <elementAPI.h>

void* OPS_AxialCurve()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 7) {
	opserr << "WARNING invalid number of arguments\n";
	opserr << "Want: limitCurve Axial tag? eleTag? Fsw? Kdeg? Fres? defType? forType?" << endln; //SDK
	opserr << "<ndI? ndJ? dof? perpDirn? delta? eleRemove?>" << endln;
	return 0;
    }
    int tag;
    int eleTag;
    double Fsw;//SDK
    double Kdeg;
    double Fres;
    int defType, forType;
    int ndI = 0;
    int ndJ = 0;
    int dof = 0;
    int perpDirn = 0;
    int eleRemove = 0;
    double delta = 0.0;

    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid Axial LimitCurve tag" << endln;
	return 0;
    }
    if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
	opserr << "WARNING invalid element tag for associated beam-column element (eleTag)\n";
	opserr << "LimitCurve Axial: " << tag << endln;
	return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &Fsw) < 0) {
	opserr << "WARNING invalid Fsw\n";
	opserr << "LimitCurve Axial: " << tag << endln;
	return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &Kdeg) < 0) {
	opserr << "WARNING invalid degrading slope Kdeg\n";
	opserr << "LimitCurve Axial: " << tag << endln;
	return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &Fres) < 0) {
	opserr << "WARNING invalid residual capacity Fres\n";
	opserr << "LimitCurve Axial: " << tag << endln;
	return 0;
    }
    if (OPS_GetIntInput(&numdata, &defType) < 0) {
	opserr << "WARNING invalid deformation type defType\n";
	opserr << "LimitCurve Axial: " << tag << endln;
	return 0;
    }
    if (OPS_GetIntInput(&numdata, &forType) < 0) {
	opserr << "WARNING invalid force type forType\n";
	opserr << "LimitCurve Axial: " << tag << endln;
	return 0;
    }
    if (defType == 2) {
	if (OPS_GetNumRemainingInputArgs() < 4) {
	    opserr << "WARNING invalid number of arguments\n";
	    opserr << "Want: limitCurve Axial tag? eleTag? Fsw? Kdeg? Fres? defType? forType?" << endln; //SDK
	    opserr << "ndI? ndJ? dof? perpDirn? <delta? eleRemove?>" << endln;
	    return 0;
	}
	if (OPS_GetIntInput(&numdata, &ndI) < 0) {
	    opserr << "WARNING invalid node I\n";
	    opserr << "LimitCurve Axial: " << tag << endln;
	    return 0;
	}
	if (OPS_GetIntInput(&numdata, &ndJ) < 0) {
	    opserr << "WARNING invalid node J\n";
	    opserr << "LimitCurve Axial: " << tag << endln;
	    return 0;
	}
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING invalid degree of freedom for drift\n";
	    opserr << "LimitCurve Axial: " << tag << endln;
	    return 0;
	}
	if (OPS_GetIntInput(&numdata, &perpDirn) < 0) {
	    opserr << "WARNING invalid direction for column length\n";
	    opserr << "LimitCurve Axial: " << tag << endln;
	    return 0;
	}
    }
    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetDoubleInput(&numdata, &delta) < 0) {
	    opserr << "WARNING invalid shift in drift surface (delta)\n";
	    opserr << "LimitCurve Axial: " << tag << endln;
	    return 0;
	}
    }
    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &eleRemove) < 0) {
	    opserr << "WARNING invalid element removal option\n";
	    opserr << "LimitCurve Axial: " << tag << endln;
	    return 0;
	}
    }
    // Parsing was successful, allocate the limit curve
    // Subtract one from dof and perpDirn for C indexing
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return 0;
    //return new AxialCurve(interp, tag, eleTag, theDomain, Fsw, //SDK
    // 			  Kdeg, Fres, defType, forType, ndI, ndJ, dof-1, perpDirn-1,
    // 			  delta, eleRemove);
}


AxialCurve::AxialCurve(Tcl_Interp *passedTclInterp, int tag, int eTag, Domain *theDom, 
			double Fsw, double Kd, double Fr, //SDK 
			int dType, int fType, 
			int ni, int nj, int df, int dirn, 
			double del, int eleRem):
LimitCurve(tag, TAG_AxialCurve), eleTag(eTag), theDomain(theDom), theElement(0),
Fsw(Fsw), Kdeg(Kd), Fres(Fr), defType(dType), forType(fType), //SDK
ndI(ni), ndJ(nj), dof(df), perpDirn(dirn), eleRemove(eleRem), delta(del)
{
	theTclInterp = passedTclInterp;
	stateFlag = 0;
	theta2 = -1.45 ; //SDK
	sigma = 0.40;    //SDK
	eps_normal= 0.0;  //SDK
	dP_old = 0.0; // SDK
	deform_old = 0.0; // SDK
	failDrift = 0.0; // SDK
	stepCounter = 0; // Terje
}


AxialCurve::AxialCurve():
LimitCurve(0, TAG_AxialCurve), eleTag(0), theDomain(0), theElement(0),
Fsw(0), Kdeg(0), Fres(0), defType(0), forType(0), //SDK
ndI(0), ndJ(0), dof(0), perpDirn(0), eleRemove(0), delta(0)
{
	stateFlag = 0;
	theta2 = -1.45 ; //SDK
	sigma = 0.40;    //SDK
	eps_normal= 0.0;  //SDK
	dP_old = 0.0; // SDK
	deform_old = 0.0; // SDK
	failDrift = 0.0; // SDK
	stepCounter = 0; // Terje
}


AxialCurve::~AxialCurve()
{
//VOID//
}


LimitCurve*
AxialCurve::getCopy(void)
{
	AxialCurve *theCopy = new AxialCurve(theTclInterp,this->getTag(),
		eleTag, theDomain, Fsw, Kdeg, Fres, defType, forType, //SDK
		ndI, ndJ, dof, perpDirn, delta, eleRemove);

	theCopy->stateFlag = stateFlag;
	theCopy->theta2 = theta2; //SDK
	theCopy->sigma = sigma;   //SDK
	theCopy->eps_normal = eps_normal;  //SDK

	theCopy->deform_old = deform_old;  //Terje
	theCopy->failDrift = failDrift;  //Terje
	theCopy->dP_old = dP_old;  //Terje

	theCopy->stepCounter = stepCounter;  //Terje

	return theCopy;
}


// check if limit state surface has been reached
int
AxialCurve::checkElementState(double springForce)
{
  static DummyStream dummy;
	// Terje
	// Count the number of times this method is called 
	// (this is equal to the number of loadsteps)
	stepCounter++;

	// if element has not been removed (element removal not fully implemented)
	if (eleRemove != 2) {

		// find associated beam-column element on first visit
		if (theElement == 0)
		{
			theElement = theDomain->getElement(eleTag);

			if (theElement == 0) {
//				g3ErrorHandler->fatal("WARNING AxialCurve - no element with tag %i exists in Domain",eleTag);
			}
			// find length between nodes if drift is desired
			if (defType == 2)
			{
				Node *nodeI = theDomain->getNode(ndI);
				Node *nodeJ = theDomain->getNode(ndJ);

				const Vector &crdI = nodeI->getCrds();
				const Vector &crdJ = nodeJ->getCrds();

				if (crdI(perpDirn) == crdJ(perpDirn)) {
//					g3ErrorHandler->warning("%s -- Nodal projection has zero component along chosen direction",
//							"AxialCurve::AxialCurve");

					oneOverL = 0.0;
				}
				else 
					oneOverL = 1.0/fabs(crdJ(perpDirn) - crdI(perpDirn));
			}
		}

		dP = 0;          //zero change in axial load
		double deform;   // value of deformation parameter from element
		double force;    // value of force parameter from element
//		double Ps;       // axial load at shear failure - for now assumed to be current axial load (Commented out by Terje)
		int result = 0;  //junk variable


		// Based on "defType" and "forType" calculate 
		// the desired response parameters "deform" and "force"
		if (defType == 1) // maximum chord rotation
		{

			Response *theRotations =0; // integer element returns in setResponse

			const char *r[1] = {"basicDeformation"}; // must be implemented in element

			Vector *rotVec; //vector of chord rotations at beam-column ends

			// set type of beam-column element response desired
			theRotations = theElement->setResponse(r, 1, dummy);
			
			if (theRotations == 0) {
			  opserr << "AxialCurve::checkElementState, defType = 1, basicDeformations not implemented in element setResponse" << endln;
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
//			g3ErrorHandler->fatal("WARNING AxialCurve - deformation type flag %i not implemented",defType);
		}
		
			Response *theForces =0;

			const char *f[1] = {"localForce"}; // does not include influence of P-delta
										 // for P-delta use forType = 0

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
			force = springForce;			//force in associated hysteretic material
		else if (forType == 1) 
			force = fabs((*forceVec)(1));	// shear 
		else if (forType == 2) 
			force = (*forceVec)(0);			//axial - positive for compression 
		else {
//			g3ErrorHandler->fatal("WARNING AxialMaterial - force type flag %i not implemented",forType);
		}

		// Determine if (deform,force) is outside limit state surface.
		// 
		// Use absolute value of deform
		// Use signed value of force to see if in compression
		   double forceSurface = findLimit(deform); // force on surface at deform
		// double deformSurface = findLimit(force); // deform on surface at force SDK
		
		
		
		
		
		
		//cout << "force = " << force << ", forceSurface = " << forceSurface << endln;

		if (stateFlag == 0) //prior to failure
		{


			if (force >= forceSurface) // on/outside failure surface
			//if (deform >= deformSurface) // on/outside failure surface SDK
			{	
				// remove element and change eleRemove flag (not fully implemented)
				if (eleRemove == 1)   
				{
					Element *theEle = theDomain->removeElement(eleTag);
					eleRemove = 2;  // do not check element again
					stateFlag = 0;

					if (theEle != 0) 
					{
						delete theEle;
					}

				} else {

					stateFlag = 1;

					dP = fabs(force) - fabs(forceSurface); //change in axial load
					opserr << "AxialCurve - failure detected at deform = " << deform << ", force = " << force << ",element: " << eleTag << endln;//SDK

					failDrift = (dP * deform_old - dP_old * deform)/(dP - dP_old);//SDK

					// Terje
					// Failure has occurred for the first time, so we 
					// print data to be employed by the limit-state function

					char myString[100];
					sprintf(myString, "AxialFailureOfElement%d.txt", eleTag);
					ofstream outputFile( myString, ios::out );

					sprintf(myString, "%d %20.8e  %20.8e  %20.8e", stepCounter, deform_old, failDrift, deform);
			
					outputFile << myString << endln;

					outputFile.close();

				}
			}
			else // inside failure surface
			{
				stateFlag = 0;
				
				dP_old = fabs(force) - fabs(forceSurface); //SDK
				deform_old = deform; //SDK


			}

		}
		else // after failure
		{
			if (force >= forceSurface) // on/outside failure surface
		//	if (deform >= deformSurface) // on/outside failure surface SDK


			{
				if (forceSurface == Fres)
					stateFlag = 4; // once at residual capacity, do not check state again
				else
					stateFlag = 2;

					dP = fabs(force) - fabs(forceSurface); //change in axial load
				// SDK	opserr << "AxialCurve - on failure surface at deform = " << deform << ", force = " << force << endln;
			}
			else // inside failure surface
			{
				stateFlag = 3;
			}
		}
	}

	return stateFlag;
}


double
AxialCurve::getDegSlope(void)
{
	return Kdeg;
}


double
AxialCurve::getResForce(void)
{
	return Fres;  
}

double
AxialCurve::getUnbalanceForce(void)
{
	return dP;
}

int
AxialCurve::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}


int
AxialCurve::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
	return -1;
}
    

void
AxialCurve::Print(OPS_Stream &s, int flag)
{
	s << "Axial Limit Curve, tag: " << this->getTag() << endln;
	s << "Fsw: " << Fsw << endln; 
	s << "eleTag: " << eleTag << endln;
	s << "nodeI: " << ndI << endln;
	s << "nodeJ: " << ndJ << endln;
	s << "deform: " << defType << endln;
	s << "force: " << forType << endln;
}


//Private Functions


// AXIAL FAILURE MODEL FROM ELWOOD (2002)
double
AxialCurve::findLimit(double x)
{
	double y = 0.0;

	if (x < 0 || x > 0.08) {
//		g3ErrorHandler->warning("Warning: Outside limits of AxialCurve");
	}

	double theta = 65.0*PI/180.0;
	double d = x-delta;

	if (d <= 0.0)
		d = 1.0e-9;

	// positive for compression
	y = ((1+tan(theta)*tan(theta))/(25*d)-tan(theta))*Fsw*tan(theta);

	//Do not allow axial load to be reduced below residual capacity (may be zero)
	//Input as positive
	if (y < Fres) {
		y = Fres;
	}

	return y;
}

/*

///////PROBABILISTIC AXIAL FAILURE MODEL /////SDK
double
AxialCurve::findLimit(double DR)
{
	double P = 0.0;
	double mu = 0.0;
	double DRa = DR - delta;
	double theta = 65.0*PI/180.0;

  
	mu = (log(DRa) - log(-0.15 - 0.23*theta2) - sigma*eps_normal)/theta2;
		
		
	P = (Ast * fyt * dc/Hs) * (1+ mu* tan(theta))/(1 - (mu/tan(theta)));
	
	//for lower limit of probabilistic model
	if (P < 0.0) {

		P = 10.0e99;

	}


	//Fres residual force
	if (P < Fres) {
		
		P = Fres;
	}


	return P;
}

*/







// AddingSensitivity:BEGIN ///////////////////////////////////
int
AxialCurve::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;
  
  else
    opserr << "WARNING: Could not set parameter in Axial Curve. " << endln;
  
  return -1;
}



int
AxialCurve::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;

	default:
		return -1;
	}


	return 0;
}


/////////////////////////////////////////ENDS

int AxialCurve::revertToStart ()
{
	stateFlag = 0;
	theta2 = -1.45 ; //SDK
	sigma = 0.40;    //SDK
	eps_normal= 0.0;  //SDK
	dP_old = 0.0; // SDK
	deform_old = 0.0; // SDK
	failDrift = 0.0; // SDK
	stepCounter = 0; // Terje

	return 0;
}
