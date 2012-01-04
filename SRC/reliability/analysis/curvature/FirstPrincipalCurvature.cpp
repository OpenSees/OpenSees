/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.8 $
// $Date: 2008-04-10 18:10:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/curvature/FirstPrincipalCurvature.cpp,v $

//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <FirstPrincipalCurvature.h>
#include <FindCurvatures.h>
#include <LimitStateFunction.h>
#include <Vector.h>
#include <math.h>

#include <iostream> 
using std::ios;

FirstPrincipalCurvature::FirstPrincipalCurvature()
:FindCurvatures(), curvatures(1)
{
}

FirstPrincipalCurvature::~FirstPrincipalCurvature()
{
}


int
FirstPrincipalCurvature::computeCurvatures(ReliabilityDomain *theReliabilityDomain)
{

	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = 
		theReliabilityDomain->getLimitStateFunctionPtr(lsf);

	int nrv = theReliabilityDomain->getNumberOfRandomVariables();

	// Get hold of 'u' and 'alpha' at the two last steps
	// Recorders needed -- MHS 10/7/2011
	/*
	const Vector &last_u = theLimitStateFunction->getFORM_u();
	const Vector &secondLast_u = theLimitStateFunction->getSecondLast_u();
	const Vector &lastAlpha = theLimitStateFunction->getFORM_alpha();
	const Vector &secondLastAlpha = theLimitStateFunction->getSecondLast_alpha();
	*/
	Vector last_u(nrv);
	Vector secondLast_u(nrv);
	Vector lastAlpha(nrv);
	Vector secondLastAlpha(nrv);

	// Compute curvature according to Der Kiureghian & De Stefano (1992), Eq.26:

	// Initial computations
	//Vector uLastMinus_u = last_u - secondLast_u;
	//const Vector &principalAxes = uLastMinus_u;

//	opserr<<"last_u    "<<last_u<<endln;
//	opserr<<"secondLast_u  "<<secondLast_u<<endln;
//	opserr<<"principalAxes"<<principalAxes<<endln;


    char fileName[30];
    strcpy (fileName, "principalAxes_.out");
    ofstream resultsOutputFile( fileName, ios::out );
    for ( int m=0; m<principalAxes.Size(); m++ ) {
      //resultsOutputFile<<principalAxes(m)<<endln;
      resultsOutputFile<<last_u(m)-secondLast_u(m)<<endln;
    }

    resultsOutputFile.close();

    //double signumProduct = secondLastAlpha ^ uLastMinus_u;
    double signumProduct = secondLastAlpha ^ last_u;
    signumProduct -= secondLastAlpha ^ secondLast_u;
    double alphaProduct = secondLastAlpha ^ lastAlpha;
    double sumSquared = 0.0;

	// Compute norm of the difference vector
	for ( int i=0; i<last_u.Size(); i++ ) {
	  //sumSquared += uLastMinus_u(i)*uLastMinus_u(i);
	  double tmp = last_u(i)-secondLast_u(i);
	  sumSquared += tmp*tmp;
	}

	double norm_uLastMinus_u = sqrt(sumSquared);

	// Check sign and compute curvature
	if (fabs(signumProduct)==(signumProduct)) {
		curvatures(0) = acos(alphaProduct) / norm_uLastMinus_u;
	}
	else {
		curvatures(0) = -acos(alphaProduct) / norm_uLastMinus_u;
	}



    strcpy (fileName, "principalCurvatures_.out");
    ofstream resultsOutputFile1( fileName, ios::out );
    
	resultsOutputFile1.precision(16);
    resultsOutputFile1<<curvatures(0)<<endln;
	

    resultsOutputFile1.close();


	return 0;
}


const Vector &
FirstPrincipalCurvature::getCurvatures()
{
  return curvatures;
}

const Vector &
FirstPrincipalCurvature::getPrincipalAxes()
{
  return principalAxes;
}

