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
#include <FORMAnalysis.h>
#include <LimitStateFunction.h>
#include <Vector.h>
#include <Matrix.h>
#include <RandomVariable.h>
#include <math.h>

#include <iostream> 
using std::ios;

FirstPrincipalCurvature::FirstPrincipalCurvature(ReliabilityDomain *passedReliabilityDomain,
                                                 FunctionEvaluator *passedFunctionEvaluator,
                                                 FORMAnalysis *passedFORMAnalysis)
:FindCurvatures(), curvatures(1)
{
    theReliabilityDomain = passedReliabilityDomain;
    theFunctionEvaluator = passedFunctionEvaluator;
    theFORMAnalysis = passedFORMAnalysis;
}

FirstPrincipalCurvature::~FirstPrincipalCurvature()
{
    
}


int
FirstPrincipalCurvature::computeCurvatures()
{

	// "Download" limit-state function from reliability domain
	int lsfTag = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsfTag);

	int nrv = theReliabilityDomain->getNumberOfRandomVariables();

	// Get hold of first principal curvature from functionEvaluator
	curvatures(0) = theFunctionEvaluator->getResponseVariable("curvatureFORM", lsfTag);
    
    // get alpha (can get from FORMAnalysis directly in the future)
    Vector alpha(nrv);
    for (int j = 0; j < nrv; j++) {
        RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(j);
        int rvTag = theRV->getTag();
        alpha(j) = theFunctionEvaluator->getResponseVariable("alphaFORM", lsfTag, rvTag);
    }
    
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
