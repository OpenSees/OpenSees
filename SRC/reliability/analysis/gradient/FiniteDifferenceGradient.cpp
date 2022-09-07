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
// /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/FiniteDifferenceGradient.cpp,v
// $

//
// Written by:
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <FiniteDifferenceGradient.h>
#include <GradientEvaluator.h>
#include <LimitStateFunction.h>
#include <ReliabilityDomain.h>
#include <Vector.h>
#include <string.h>

FiniteDifferenceGradient::FiniteDifferenceGradient(
    FunctionEvaluator *passedGFunEvaluator,
    ReliabilityDomain *passedReliabilityDomain,
    Domain *passedOpenSeesDomain)

    : GradientEvaluator(passedReliabilityDomain, passedGFunEvaluator),
      theOpenSeesDomain(passedOpenSeesDomain) {
    int nrv = passedReliabilityDomain->getNumberOfRandomVariables();
    grad_g = new Vector(nrv);
}

FiniteDifferenceGradient::~FiniteDifferenceGradient() {
    if (grad_g != 0) delete grad_g;
}

const Vector &FiniteDifferenceGradient::getGradient() { return *grad_g; }

int FiniteDifferenceGradient::computeGradient(double g) {
    // note FiniteDifferentGradient presumes that the expression has
    // already been evaluated once with
    // default parameter values and the result is passed in with variable
    // g. Therefore it is only computing the perturbations from this
    // default state

    // Initialize gradient vector
    grad_g->Zero();

    // get limit-state function from reliability domain
    int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
    LimitStateFunction *theLimitStateFunction =
        theReliabilityDomain->getLimitStateFunctionPtr(lsf);
    const char *lsfExpression = theLimitStateFunction->getExpression();

    // get RVs created in the reliability domain
    int nrv = this->theReliabilityDomain->getNumberOfRandomVariables();

    // now loop through to create gradient vector
    // for all RVs
    for (int i = 0; i < nrv; i++) {
        // get RV
        auto *theRV =
            this->theReliabilityDomain->getRandomVariablePtrFromIndex(i);
        if (theRV == 0) {
            opserr << "ERROR: can't get RV " << i
                   << " -- FiniteDifferenceGradient::computeGradient\n";
            return -1;
        }

        // get RV parameter
        int param_indx =
            theReliabilityDomain->getParameterIndexFromRandomVariableIndex(
                i);

        auto *theParam =
            theOpenSeesDomain->getParameterFromIndex(param_indx);
        if (theParam == 0) {
            opserr << "ERROR: can't get param " << i
                   << " -- FiniteDifferenceGradient::computeGradient\n";
            return -1;
        }

        // use parameter defined perturbation
        double h = theParam->getPerturbation();
        double original = theParam->getValue();
        theParam->update(original + h);

        // set perturbed values in the variable namespace
        if (theFunctionEvaluator->setVariables() < 0) {
            opserr << "ERROR FiniteDifferenceGradient -- error "
                      "setting variables in namespace"
                   << endln;
            theParam->update(original);
            return -1;
        }

        // run analysis
        if (theFunctionEvaluator->runAnalysis() < 0) {
            opserr << "ERROR FiniteDifferenceGradient -- error "
                      "running analysis"
                   << endln;
            theParam->update(original);
            return -1;
        }

        // evaluate LSF and obtain result
        theFunctionEvaluator->setExpression(lsfExpression);

        // perturbed lsf
        double g_perturbed = theFunctionEvaluator->evaluateExpression();
        (*grad_g)(i) = (g_perturbed - g) / h;

        // return parameter values to previous state
        theParam->update(original);
    }

    return 0;
}
