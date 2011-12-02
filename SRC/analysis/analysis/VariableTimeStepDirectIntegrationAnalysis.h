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
                                                                        
// $Revision: 1.4 $
// $Date: 2010-06-01 23:48:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/VariableTimeStepDirectIntegrationAnalysis.h,v $
                                                                        
                                                                        
#ifndef VariableTimeStepDirectIntegrationAnalysis_h
#define VariableTimeStepDirectIntegrationAnalysis_h

// Written: fmk 
// Created: 10/00
// Revision: A
//
// Description: This file contains the class definition for 
// VariableTimeStepDirectIntegrationAnalysis. VariableTimeStepDirectIntegrationAnalysis 
// is a subclass of DirectIntegrationAnalysis. It is used to perform a 
// dynamic analysis on the FE\_Model using a direct integration scheme.  
//
// What: "@(#) VariableTimeStepDirectIntegrationAnalysis.h, revA"

#include <DirectIntegrationAnalysis.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class TransientIntegrator;
class LinearSOE;
class EquiSolnAlgo;
class ConvergenceTest;

class VariableTimeStepDirectIntegrationAnalysis: public DirectIntegrationAnalysis
{
  public:
    VariableTimeStepDirectIntegrationAnalysis(Domain &theDomain,
					      ConstraintHandler &theHandler,
					      DOF_Numberer &theNumberer,
					      AnalysisModel &theModel,
					      EquiSolnAlgo &theSolnAlgo,
					      LinearSOE &theSOE,
					      TransientIntegrator &theIntegrator,
					      ConvergenceTest *theTest =0);
    virtual ~VariableTimeStepDirectIntegrationAnalysis();

    int analyze(int numSteps, double dT, double dtMin, double dtMax, int Jd);

  protected:
    virtual double determineDt(double dT, double dtMin, double dtMax, int Jd,
			       ConvergenceTest *theTest);

  private:
};

#endif

