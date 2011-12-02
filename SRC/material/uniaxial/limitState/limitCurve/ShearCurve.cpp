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

// $Revision: 1.1 $
// $Date: 2006-02-07 23:15:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/limitCurve/ShearCurve.cpp,v $
                                                                        
// Written: KJE
// Created: Sept 2002

#include <ShearCurve.h>
#include <G3Globals.h>
#include <math.h>
#include <ElementResponse.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <Vector.h>
#include <float.h>


ShearCurve::ShearCurve(int tag, int eTag, Domain *theDom, 
			double r, double f, double B, double H, double D, double Fsw, //SDK
			double Kd, double Fr,
			int dType, int fType,
			int ni, int nj, int df, int dirn, double dd):
LimitCurve(tag, TAG_ShearCurve), eleTag(eTag), theDomain(theDom), theElement(0),
Kdeg(Kd), Fres(Fr), defType(dType), forType(fType),
rho(r), fc(f), b(B), h(H), d(D), Fsw(Fsw), 
ndI(ni), ndJ(nj), dof(df), perpDirn(dirn), delta(dd)
{
	stateFlag = 0;

	theta1 = 0.037;
	theta4 = -0.027;
	theta5 = -0.034;
	sigma = 0.3;
	eps_normal = 0.0;
}

ShearCurve::ShearCurve():
LimitCurve(0, TAG_ShearCurve), eleTag(0), theDomain(0), theElement(0),
Kdeg(0), Fres(0), defType(0), forType(0),
rho(0), fc(0), b(0), d(0), h(0), Fsw(0) 
{
//VOID//
}

ShearCurve::~ShearCurve()
{
//VOID//
}


LimitCurve*
ShearCurve::getCopy(void)
{
	ShearCurve *theCopy = new ShearCurve(this->getTag(),
		eleTag, theDomain, rho, fc, b, h, d,  Fsw, 
		Kdeg, Fres, defType, forType, ndI, ndJ, dof, perpDirn, delta);

	theCopy->stateFlag = stateFlag;
	theCopy->theta1 = theta1;
	theCopy->theta4 = theta4;
	theCopy->theta5 = theta5;
	theCopy->sigma = sigma;
	theCopy->eps_normal = eps_normal;

	return theCopy;
}


// check if limit state surface has been reached
int
ShearCurve::checkElementState(double springForce)
{
	// find associated beam-column elementon first visit
	if (theElement == 0)
	{
		theElement = theDomain->getElement(eleTag);

		if (theElement == 0) {
//			g3ErrorHandler->fatal("WARNING ShearCurve - no element with tag %i exists in Domain",eleTag);
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
	int result;		//junk variable


	// Based on "defType" and "forType" calculate 
	// the desired response parameters "deform" and "force"
	if (defType == 1) // maximum chord rotations
	{

		Response *theRotations =0; // integer element returns in setResponse

		const char *r[1] = {"basicDeformations"}; // must be implemented in element

		Information	*rotInfoObject =0;   

		Vector *rotVec; //vector of chord rotations at beam-column ends

		// set type of beam-column element response desired
		theRotations = theElement->setResponse(r, 1, *rotInfoObject);

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

//opserr << "Drift ndI: " << dispI(dof) << endln;
//opserr << "Drift ndJ: " << dispJ(dof) << endln;

		
		// calc drift
		double dx = fabs(dispJ(dof)-dispI(dof));
		deform = dx*oneOverL;
	}
	else {
//		g3ErrorHandler->fatal("WARNING ShearCurve - deformation type flag %i not implemented",defType);
	}

		Response *theForces =0;

		const char *f[1] = {"localForce"}; // does not include influence of P-delta
								     // for P-delta use forType = 0

		Information	*forInfoObject =0;

		Vector *forceVec; //vector of basic forces from beam column

		// set type of beam-column element response desired
		theForces    = theElement->setResponse(f, 1, *forInfoObject);

		// put element response in the vector of "myInfo"
		result += theForces->getResponse();

		// access the myInfo vector containing the response (new for Version 1.2)
		Information &theInfo = theForces->getInformation();
		forceVec = (theInfo.theVector);

	// Local forces (assuming no element loads)
	if (forType == 0)
		force = fabs(springForce);    // force in associated LimitState material
	else if (forType == 1) 
		force = fabs((*forceVec)(1)); // shear
	else if (forType == 2) 
		force = fabs((*forceVec)(0)); // axial
	else {
//		g3ErrorHandler->fatal("WARNING ShearCurve - force type flag %i not implemented",forType);
	}

	P = fabs((*forceVec)(0));

	// Determine if (deform,force) is outside limit state surface.
	// 
	// Use absolute value of deform and force
	double forceSurface = findLimit(deform); // force on surface at deform

//	double deformSurface = findLimit(force); // deform on surface at force	SDK


//opserr << "The shear force in the element is: " << force << endln;


//opserr << "State flag............................:" << stateFlag << endln;
//opserr << "Shear force in the column.........: " << force << endln;
//opserr << "forceSurface: " << forceSurface << endln;

	if (stateFlag == 0) //prior to failure
	{
		if (force >= forceSurface) // on/outside failure surface
	//	  if (deform >= deformSurface) // on/outside failure surface	SDK

		{	
			stateFlag = 1;
			setDegSlope(force, deform);
//			g3ErrorHandler->warning("ShearCurve - failure detected at deform = %f (Kdeg = %f) ", deform, Kdeg);
//opserr << "*********************" << endln;
        opserr << "ShearCurve - failure detected....."<< endln;//SDK

       // opserr << "deformSurface: " << deformSurface << endln; //SDK
  
      



//opserr << "Capacity: " << forceSurface << endln;
//opserr << "*********************" << endln;
		}
		else // inside failure surface
		{
			stateFlag = 0;
		}
	}
	else //after failure
	{
		if (force >= forceSurface) // on/outside failure surface
		//  if (deform >= deformSurface) // on/outside failure surface SDK
		{	
			stateFlag = 2;
		}
		else // inside failure surface
		{
			stateFlag = 3;
		}
	}

	return stateFlag;
}


double
ShearCurve::getDegSlope(void)
{
	return Kdeg;
}


double
ShearCurve::getResForce(void)
{
	return Fres;  
}

double
ShearCurve::getUnbalanceForce(void)
{
	//Do nothing for this class
	return 0.0;
}

int
ShearCurve::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int
ShearCurve::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
	return -1;
}
    
void
ShearCurve::Print(OPS_Stream &s, int flag)
{
	s << "Shear Limit Curve, tag: " << this->getTag() << endln;
	s << "rho: " << rho << endln;
	s << "fc: " << fc <<" psi" << endln;
	s << "b,d: " << b <<" in., "<< d << " in." << endln;
	s << "eleTag: " << eleTag << endln;
	s << "nodeI: " << ndI << endln;
	s << "nodeJ: " << ndJ << endln;
	s << "deform: " << defType << endln;
	s << "force: " << forType << endln;
}



// Private Functions

void
ShearCurve::setDegSlope(double V, double Dshear)
{
	if (Kdeg > 0.0)
	{
		// Calculate degrading slope based on point of shear failure and 
		// calculated deformation at axial failure based on current axial 
		// load and axial failure model by Elwood (2002).
		// If positive, Kdeg is assumed equal to the flexural stiffness
		double mu = 0.0;
		double theta = 65.0*PI/180.0;
		double Daxial;
		//double Fsw = 0.0; //SDK



		//Fsw = Ast * fyt * dc/s;
		Daxial = 0.04*(1+tan(theta)*tan(theta))/(tan(theta)+P/Fsw/tan(theta));
		
		//mu = ((P/Fsw)-1)/(P/Fsw/tan(theta)+ tan(theta)); //SDK
		//Daxial = 0.184*exp(-1.45*mu); //SDK

		if (defType == 2)
		{
			double K = -V/(Daxial-Dshear)*oneOverL;
			Kdeg = 1/(1/K - 1/Kdeg);
		}
		else {
//			g3ErrorHandler->fatal("WARNING ShearCurve - must use defType = 2 for calculated degrading slope");
		}
		
	}

}

// Empirical drift capacity model at shear failure from Elwood (2002)
double
ShearCurve::findLimit(double DR) //SDK

{

	double V = 0.0; //Shear in kips!!
	//double DR = 0.0; // drift SDK
	
	if (DR < 0.01)
		V = 9.9e9; //no shear failure below drift ratio of 1%
	else {

		V = 500*(0.03+delta+4*rho-DR-0.025*P/b/h/(fc/1000))*(b*d*sqrt(fc)/1000);
		
/*/ This is Ling's probabilistic capacity model:		
		double s_over_d = 0.75; // hoop spacing to depth ration
		double a_over_d = 0.5*58.0/d; // aspect ratio
		

		//V = (exp(-sigma*eps_normal)*DR - theta1 - theta4*P/b/h/fc/1000 - theta5*s_over_d 
		//	- (0.015-0.21*theta1)*a_over_d ) 
		//	* b*d*6.0*sqrt(fc)/1000 / (0.0053-0.42*theta1);



		DR = (theta1 + (0.0053-0.42*theta1)* V*1000/ b/d/6.0/sqrt(fc)+ theta4*P/b/h/fc/1000 
			+ theta5*s_over_d + (0.015-0.21*theta1)*a_over_d )*exp(sigma*eps_normal); //SDK
*/				


	}
	if (V < 0.0)
		V = 0.0;


	return V;
//	return DR; //SDK
}







// AddingSensitivity:BEGIN ///////////////////////////////////
int
ShearCurve::setParameter(const char **argv, int argc, Information &info)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"theta1") == 0) {
		info.theType = DoubleType;
		return 1;
	}
	if (strcmp(argv[0],"theta4") == 0) {
		info.theType = DoubleType;
		return 2;
	}
	if (strcmp(argv[0],"theta5") == 0) {
		info.theType = DoubleType;
		return 3;
	}
	if (strcmp(argv[0],"sigma") == 0) {
		info.theType = DoubleType;
		return 4;
	}
	if (strcmp(argv[0],"eps_normal") == 0) {
		info.theType = DoubleType;
		return 5;
	}
	
	if (strcmp(argv[0],"fc") == 0) {
		info.theType = DoubleType;
		return 6;
	}
	
	
	else
		opserr << "WARNING: Could not set parameter in Shear Curve. " << endln;
                
	return -1;
}



int
ShearCurve::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->theta1 = info.theDouble;
		break;
	case 2:
		this->theta4 = info.theDouble;
		break;
	case 3:
		this->theta5 = info.theDouble;
		break;
	case 4:
		this->sigma = info.theDouble;
		break;
	case 5:
		this->eps_normal = info.theDouble;
		break;
	case 6:
		this->fc = info.theDouble;
		break;

	default:
		return -1;
	}


	return 0;
}



int ShearCurve::revertToStart ()
{
	stateFlag = 0;
	return 0;
}
