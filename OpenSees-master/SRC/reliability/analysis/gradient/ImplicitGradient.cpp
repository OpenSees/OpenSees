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

// $Revision: 1.21 $
// $Date: 2010-09-13 21:38:31 $
// $Source:
// /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/ImplicitGradient.cpp,v
// $

//
// Written by:
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <GradientEvaluator.h>
#include <ImplicitGradient.h>
#include <Integrator.h>  //Abbas
#include <LimitStateFunction.h>
#include <ReliabilityDomain.h>
#include <SensitivityAlgorithm.h>
#include <Vector.h>
#include <string.h>

ImplicitGradient::ImplicitGradient(
    FunctionEvaluator *passedGFunEvaluator,
    ReliabilityDomain *passedReliabilityDomain,
    Domain *passedOpenSeesDomain,
    // SensitivityAlgorithm *theAlgo)
    Integrator *theAlgo)
    : GradientEvaluator(passedReliabilityDomain, passedGFunEvaluator),
      theOpenSeesDomain(passedOpenSeesDomain) {
    theSensAlgo = theAlgo;

    int nrv = passedReliabilityDomain->getNumberOfRandomVariables();
    // int nparam = theOpenSeesDomain->getNumParameters();
    grad_g = new Vector(nrv);
}

ImplicitGradient::~ImplicitGradient() {
    if (grad_g != 0) delete grad_g;
}

const Vector &ImplicitGradient::getGradient() { return *grad_g; }

int ImplicitGradient::computeGradient(double g) {
    // Initialize gradient vector
    grad_g->Zero();

    // get limit-state function from reliability domain
    int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
    LimitStateFunction *theLimitStateFunction =
        theReliabilityDomain->getLimitStateFunctionPtr(lsf);
    const char *lsfExpression = theLimitStateFunction->getExpression();

    // get RVs created in the reliability domain
    int nrv = this->theReliabilityDomain->getNumberOfRandomVariables();

    // first check for dg/dimplicit partials
    for (int i = 0; i < nrv; i++) {
        // get RV
        auto *theRV =
            this->theReliabilityDomain->getRandomVariablePtrFromIndex(i);
        if (theRV == 0) {
            opserr << "ERROR: can't get RV " << i
                   << " -- ImplictGradient::computeGradient\n";
            return -1;
        }

        // get RV tag
        int rvTag = theRV->getTag();

        // get grad expression
        const char *gradExpression =
            theLimitStateFunction->getGradientExpression(rvTag);

        // not provided means zero
        if (gradExpression == 0) {
            opserr
                << "WARNING: lsf " << lsf
                << ": gradient expression for "
                   "random variable "
                << rvTag
                << " is not defined -- ImplicitGradient::computeGradient";
            return -1;
        }

        // set expression
        theFunctionEvaluator->setExpression(gradExpression);

        // update variables
        if (theFunctionEvaluator->setVariables() < 0) {
            opserr << "ERROR ImplicitGradient -- error setting "
                      "variables in namespace"
                   << endln;
            return -1;
        }

        // evaluate grad expression
        (*grad_g)(i) = theFunctionEvaluator->evaluateExpression();

        // Reset limit state function in evaluator
        theFunctionEvaluator->setExpression(lsfExpression);
    }

    return 0;
}
