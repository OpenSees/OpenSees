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
#include <FORMAnalysis.h>
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
                                   FORMAnalysis *passedFORMAnalysis,
                                   HessianEvaluator *passedHessianEvaluator,
                                   ProbabilityTransformation *passedTransformation)
  :FindCurvatures(), eigenvalue(0), eigenvector(0), curvatures(1)
{
    theReliabilityDomain = passedReliabilityDomain;
    theOpenSeesDomain = passedOpenSeesDomain;
    theFunctionEvaluator = passedFunctionEvaluator;
    theFORMAnalysis = passedFORMAnalysis;
    theHessianEvaluator = passedHessianEvaluator;
    theProbabilityTransformation = passedTransformation;
}

CurvatureFitting::~CurvatureFitting()
{
  if (eigenvalue != 0)
    delete [] eigenvalue;
  if (eigenvector != 0)
    delete [] eigenvector;
}

#ifdef _WIN32

extern "C" int DGEEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *WR, double *WI,
                     double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);

#else

extern "C" int dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *WR, double *WI,
                      double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);

#endif

int
CurvatureFitting::computeCurvatures()
{

	// "Download" limit-state function from reliability domain
	int lsfTag = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsfTag);

	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
    int numberOfParameters = theOpenSeesDomain->getNumParameters();
    
    // get alpha from FunctionEvaluator (can get straight from FORM in the future)
    Vector alpha(nrv);
    Vector xStar(nrv);
    Vector uStar(nrv);
    for (int j = 0; j < nrv; j++) {
        RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(j);
        int rvTag = theRV->getTag();
        alpha(j) = theFunctionEvaluator->getResponseVariable("alphaFORM", lsfTag, rvTag);
        xStar(j) = theFunctionEvaluator->getResponseVariable("designPointXFORM", lsfTag, rvTag);
        uStar(j) = theFunctionEvaluator->getResponseVariable("designPointUFORM", lsfTag, rvTag);
        theRV->setCurrentValue(xStar(j));
    }
    
    // Gram Schmidt orthogonalization on the alpha hat vector
    Matrix Qmat(nrv,nrv);
    this->gramSchmidt(alpha,Qmat);
    
    // compute Hessian
    int result = theHessianEvaluator->computeHessian();
    if (result < 0) {
        opserr << "CurvatureFitting::computeCurvatures() - " << endln
               << " could not compute hessian of the limit-state function. " << endln;
        return -1;
    }
    
    // transform Hessian
    Matrix hessU(nrv,nrv);
    Matrix temp_hess(numberOfParameters,numberOfParameters);
    temp_hess = theHessianEvaluator->getHessian();
    
    // map hessian from all parameters to just RVs
    Matrix hessX(nrv,nrv);
    for (int j = 0; j < nrv; j++) {
        int param_indx_j = theReliabilityDomain->getParameterIndexFromRandomVariableIndex(j);
        for (int k = 0; k <= j; k++) {
            int param_indx_k = theReliabilityDomain->getParameterIndexFromRandomVariableIndex(k);
            hessX(j,k) = temp_hess(param_indx_j,param_indx_k);
            hessX(k,j) = hessX(j,k);
        }
    }
    
    // Get Jacobian x-space to u-space
    Matrix Jxu(nrv,nrv);
    result = theProbabilityTransformation->getJacobian_x_to_u(Jxu);
    if (result < 0) {
        opserr << "CurvatureFitting::computeCurvatures() - " << endln
               << " could not transform Jacobian from x to u." << endln;
        return -1;
    }
    
    // Hessian in standard normal space (if the transformation is linear)
    hessU.addMatrixTripleProduct(0.0,Jxu,hessX,1.0);
    
    // gradient in original space
    Vector gradientX;
    theFORMAnalysis->getStorage("gradientXFORM",lsfTag,gradientX);
    
    // Now add the nonlinear term
    // note in the future there's probably a better way to do this directly with the 
    // probability transformation and information coming from random variables, here
    // use finite differences
    Matrix hessNL(nrv,nrv);
    double hU = 1/2000.0;
    for (int i = 0; i < nrv; i++) {
        for (int j = 0; j <= i; j++) {
            Vector tHess(nrv);
            for (int k = 0; k < nrv; k++) {
                Vector uTemp(uStar);
                Vector xTemp(xStar);
                double tempsum = 0;
                
                // unperturbed value
                tempsum += xTemp(k);
                
                // perturb ith u variable
                uTemp(i) += hU;
                result = theProbabilityTransformation->transform_u_to_x(uTemp,xTemp);
                tempsum -= xTemp(k);
                
                // perturb jth u variable
                uTemp = uStar;
                uTemp(j) += hU;
                result = theProbabilityTransformation->transform_u_to_x(uTemp,xTemp);
                tempsum -= xTemp(k);
                
                // perturb both
                uTemp = uStar;
                uTemp(i) += hU;
                uTemp(j) += hU;
                result = theProbabilityTransformation->transform_u_to_x(uTemp,xTemp);
                tempsum += xTemp(k);
                
                tHess(k) = tempsum/hU/hU;
            }
            
            hessNL(i,j) = tHess ^ gradientX;
            hessNL(j,i) = hessNL(i,j);
        }
    }
    
    // norm of gradient in standard normal space
    Vector gradientU;
    theFORMAnalysis->getStorage("gradientUFORM",lsfTag,gradientU);
    
    // add linear and nonlinear parts of the hessU, then compute A matrix
    hessU += hessNL;
    Matrix Amat(nrv,nrv);
    Amat.addMatrixTripleProduct(0.0,Qmat,hessU,1.0/gradientU.Norm());
    
    // reduce A matrix to first (n-1)x(n-1) entries
    int nred = nrv-1;
    double *Adata = new double[nred*nred];
    Matrix Areduced(Adata,nred,nred);
    for (int i = 0; i < nred; i++) {
        for (int j = 0; j < nred; j++)
            Areduced(i,j) = Amat(i,j);
    }
    
    // do not compute left eigenvalues and eigenvectors
    char *jobvl = "N";

    // compute right eigenvalues and eigenvectors
    char *jobvr = "V";

    // reduced A matrix data
    double *Kptr = Adata;

    // leading dimension of K
    int ldK = nred;

    // allocate memory for eigenvalues (imaginary part)
    double *alphaI = new double [nred];

    if (eigenvalue != 0)
        delete [] eigenvalue;

    // and real part
    eigenvalue = new double [nred];
    curvatures.setData(eigenvalue, nred);

    // dummy left eigenvectors
    double vl[1];

    // leading dimension of dummy left eigenvectors
    int ldvl = 1;

    // allocate memory for right eigenvectors
    if (eigenvector != 0)
        delete [] eigenvector;
    eigenvector = new double [nred*nred];

    // leading dimension of right eigenvectors
    int ldvr = nred;

    // dimension of the workspace array
    int lwork = 64*nred + 1;

    // allocate memory for workspace array
    double *work = new double [lwork];

    // output information
    int info = 0;

    // call the LAPACK eigenvalue subroutine for reduced A matrix
    
#ifdef _WIN32
    DGEEV(jobvl, jobvr, &nred, Kptr, &ldK, eigenvalue, alphaI,
          vl, &ldvl, eigenvector, &ldvr, work, &lwork, &info);
#else
    dgeev_(jobvl, jobvr, &nred, Kptr, &ldK, eigenvalue, alphaI,
           vl, &ldvl, eigenvector, &ldvr, work, &lwork, &info);
#endif

    int lworkOpt = (int) work[0];
    if (lwork < lworkOpt) {
        opserr << "CurvatureFitting::computeCurvatures() - optimal workspace size "
                << lworkOpt << " is larger than provided workspace size "
                << lwork << " consider increasing workspace\n";
    }

    // clean up the memory
    delete [] alphaI;
    delete [] work;

    delete [] Adata;

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

