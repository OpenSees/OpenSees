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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/curvature/CurvatureFitting.cpp,v $

//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <CurvatureFitting.h>
#include <FindCurvatures.h>
#include <LimitStateFunction.h>
#include <Vector.h>
#include <Matrix.h>
#include <RandomVariable.h>
#include <HessianEvaluator.h>
#include <math.h>

#include <iostream> 
using std::ios;

CurvatureFitting::CurvatureFitting(ReliabilityDomain *passedReliabilityDomain,
                                   Domain *passedOpenSeesDomain,
                                   FunctionEvaluator *passedFunctionEvaluator,
                                   HessianEvaluator *passedHessianEvaluator,
                                   ProbabilityTransformation *passedTransformation)
:FindCurvatures(), curvatures(1)
{
    theReliabilityDomain = passedReliabilityDomain;
    theOpenSeesDomain = passedOpenSeesDomain;
    theFunctionEvaluator = passedFunctionEvaluator;
    theHessianEvaluator = passedHessianEvaluator;
    theProbabilityTransformation = passedTransformation;
}

CurvatureFitting::~CurvatureFitting()
{
    
}


int
CurvatureFitting::computeCurvatures()
{

	// "Download" limit-state function from reliability domain
	int lsfTag = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsfTag);

	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
    int numberOfParameters = theOpenSeesDomain->getNumParameters();
    
    // get alpha from FunctionEvaluator
    Vector alpha(nrv);
    Vector xStar(nrv);
    for (int j = 0; j < nrv; j++) {
        RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(j);
        int rvTag = theRV->getTag();
        alpha(j) = theFunctionEvaluator->getResponseVariable("alphaFORM", lsfTag, rvTag);
        xStar(j) = theFunctionEvaluator->getResponseVariable("designPointXFORM", lsfTag, rvTag);
    }
    
    // Gram Schmidt orthogonalization on the alpha hat vector
    Matrix temp(nrv,nrv);
    this->gramSchmidt(alpha,temp);
    
    // compute Hessian
    int result = theHessianEvaluator->computeHessian();
    if (result < 0) {
        opserr << "CurvatureFitting::computeCurvatures() - " << endln
               << " could not compute hessian of the limit-state function. " << endln;
        return -1;
    }
    
    // transform Hessian
    Matrix temp_hess(numberOfParameters,numberOfParameters);
    temp_hess = theHessianEvaluator->getHessian();
    Matrix hessU(nrv,nrv);
    
    // map hessian from all parameters to just RVs
    for (int j = 0; j < nrv; j++) {
        int param_indx_j = theReliabilityDomain->getParameterIndexFromRandomVariableIndex(j);
        for (int k = 0; k <= j; k++) {
            int param_indx_k = theReliabilityDomain->getParameterIndexFromRandomVariableIndex(k);
            hessU(j,k) = temp_hess(param_indx_j,param_indx_k);
            hessU(k,j) = hessU(j,k);
        }
    }
    
    // Get Jacobian x-space to u-space
    Matrix Jxu(nrv,nrv);
    result = theProbabilityTransformation->getJacobian_x_to_u(Jxu);         
    Matrix hessU2(nrv,nrv);
    
    // Gradient in standard normal space
    hessU2.addMatrixTripleProduct(0.0,Jxu,hessU,1.0);
    
    // compute A matrix
    hessU2.addMatrixTripleProduct(1.0,temp,hessU2,1.0);
    
    // still need to normalize by norm of gradient in standard normal space
    
    // eigenvalues of reduced A matrix
    

	return 0;
}


const Vector &
CurvatureFitting::getCurvatures()
{
    return curvatures;
}

const Vector &
CurvatureFitting::getPrincipalAxes()
{
    return principalAxes;
}

