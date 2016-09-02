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
                                                                        
                                                                        
#ifndef PFEMAnalysis_h
#define PFEMAnalysis_h

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

#include <DirectIntegrationAnalysis.h>

class PFEMAnalysis : public DirectIntegrationAnalysis
{
public:
    PFEMAnalysis(Domain &theDomain, 
                 ConstraintHandler &theHandler,
                 DOF_Numberer &theNumberer,
                 AnalysisModel &theModel,
                 EquiSolnAlgo &theSolnAlgo,		   
                 LinearSOE &theSOE,
                 TransientIntegrator &theIntegrator,
                 ConvergenceTest *theTest,
                 double max, double min, double g, double r);

    int analyze();

    virtual ~PFEMAnalysis();

private:

    int identify();

    double dtmax;
    double dtmin;
    double gravity;
    double ratio;
    double dt;
    double next;
    int curr;
    bool newstep;
};

#endif
