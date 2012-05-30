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
  :FindCurvatures(), eigenvalue(0), eigenvector(0), curvatures(1)
{
    theReliabilityDomain = passedReliabilityDomain;
    theOpenSeesDomain = passedOpenSeesDomain;
    theFunctionEvaluator = passedFunctionEvaluator;
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


    double *hessU2data = new double[nrv*nrv];
    Matrix hessU2(hessU2data,nrv,nrv);
    
    // Gradient in standard normal space
    hessU2.addMatrixTripleProduct(0.0,Jxu,hessU,1.0);
    
    // compute A matrix
    hessU2.addMatrixTripleProduct(1.0,temp,hessU2,1.0);
    
    // still need to normalize by norm of gradient in standard normal space
    


    // eigenvalues of reduced A matrix
    
    // do not compute left eigenvalues and eigenvectors
    char *jobvl = "N";

    // compute right eigenvalues and eigenvectors
    char *jobvr = "V";

    // stiffness matrix data
    double *Kptr = hessU2data;

    // leading dimension of K
    int ldK = nrv;

    // allocate memory for eigenvalues (imaginary part)
    double *alphaI = new double [nrv];

    if (eigenvalue != 0)
        delete [] eigenvalue;

    // and real part
    eigenvalue = new double [nrv];
    curvatures.setData(eigenvalue, nrv);

    // dummy left eigenvectors
    double vl[1];

    // leading dimension of dummy left eigenvectors
    int ldvl = 1;

    // allocate memory for right eigenvectors
    if (eigenvector != 0)
        delete [] eigenvector;
    eigenvector = new double [nrv*nrv];

    // leading dimension of right eigenvectors
    int ldvr = nrv;

    // dimension of the workspace array
    int lwork = 4*nrv + 1;

    // allocate memory for workspace array
    double *work = new double [lwork];

    // output information
    int info = 0;

    // call the LAPACK eigenvalue subroutine
#ifdef _WIN32
    DGEEV(jobvl, jobvr, &nrv, Kptr, &ldK, eigenvalue, alphaI,
          vl, &ldvl, eigenvector, &ldvr, work, &lwork, &info);
#else
    dgeev_(jobvl, jobvr, &nrv, Kptr, &ldK, eigenvalue, alphaI,
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

    delete [] hessU2data;

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

