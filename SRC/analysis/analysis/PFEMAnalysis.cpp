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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-11-30 13:01:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/PFEMAnalysis.h,v $
                                                                        
                                                                        
// Written: Minjie  
// Created: Nov 2012
// Revision: A
//
// Description: This file contains the class definition for 
// PFEMAnalysis. PFEMAnalysis is a 
// subclass of DirectIntegrationAnalysis. It is used to perform a 
// dynamic analysis for PFEM using a direct integration scheme.
//
// What: "@(#) PFEMAnalysis.h, revA"

#include <PFEMAnalysis.h>
#include <EquiSolnAlgo.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <ConvergenceTest.h>
#include <TransientIntegrator.h>
#include <Domain.h>

PFEMAnalysis::PFEMAnalysis(Domain& theDomain, 
                           ConstraintHandler& theHandler,
                           DOF_Numberer& theNumberer,
                           AnalysisModel& theModel,
                           EquiSolnAlgo& theSolnAlgo,		   
                           LinearSOE& theSOE,
                           TransientIntegrator& theIntegrator,
                           ConvergenceTest* theTest,
                           double max, double min, double r)
    :DirectIntegrationAnalysis(theDomain,theHandler,theNumberer,theModel,theSolnAlgo,
                               theSOE,theIntegrator,theTest),
     dtmax(max), dtmin(min), ratio(r),
     dt(max), next(max), newstep(true)
{
    
}

PFEMAnalysis::~PFEMAnalysis()
{
}

int 
PFEMAnalysis::analyze()
{
    Domain* theDomain = this->getDomainPtr();
    double current = theDomain->getCurrentTime();
    if(newstep) {
        next = current + dtmax;
    }
    bool instep = false;

    while(true) {

        // analyze
        opserr<<"\n\nTime = "<<current<<", dt = "<<dt<<"\n\n";
        int converged  = DirectIntegrationAnalysis::analyze(1, dt);

        // if failed
        if(converged < 0) {
            dt *= ratio;
            if(dt < dtmin) {
                return -1;
            }
            instep = true;
            newstep = false;
            continue;
        }

        // if converged
        if(instep) {
            current = theDomain->getCurrentTime();
            dt = next - current;
        } else {
            newstep = true;
            dt = dtmax;
        }

        return 0;
    }

}
