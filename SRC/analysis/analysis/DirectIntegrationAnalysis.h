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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-03-04 00:48:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/DirectIntegrationAnalysis.h,v $
                                                                        
                                                                        
#ifndef DirectIntegrationAnalysis_h
#define DirectIntegrationAnalysis_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for 
// DirectIntegrationAnalysis. DirectIntegrationAnalysis is a 
// subclass of TransientAnalysis. It is used to perform a 
// dynamic analysis on the FE\_Model using a direct integration scheme.
//
// What: "@(#) DirectIntegrationAnalysis.h, revA"

#include <TransientAnalysis.h>
// AddingSensitivity:BEGIN //////////////////////////////////
#ifdef _RELIABILITY
#include <SensitivityAlgorithm.h>
#endif
// AddingSensitivity:END ////////////////////////////////////

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class TransientIntegrator;
class LinearSOE;
class EquiSolnAlgo;

class DirectIntegrationAnalysis: public TransientAnalysis
{
  public:
    DirectIntegrationAnalysis(Domain &theDomain, 
			      ConstraintHandler &theHandler,
			      DOF_Numberer &theNumberer,
			      AnalysisModel &theModel,
			      EquiSolnAlgo &theSolnAlgo,		   
			      LinearSOE &theSOE,
			      TransientIntegrator &theIntegrator);
    virtual ~DirectIntegrationAnalysis();

    void clearAll(void);	    
    
    int analyze(int numSteps, double dT);
    int initialize(void);

    int domainChanged(void);
    
    int setAlgorithm(EquiSolnAlgo &theAlgorithm);
    int setIntegrator(TransientIntegrator &theIntegrator);
    int setLinearSOE(LinearSOE &theSOE); 
    
    int checkDomainChange(void);
    EquiSolnAlgo *getAlgorithm(void);
    TransientIntegrator *getIntegrator();
// AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
	int setSensitivityAlgorithm(SensitivityAlgorithm *theSensitivityAlgorithm);
#endif
// AddingSensitivity:END /////////////////////////////////
    
  protected:
    
  private:
    ConstraintHandler 	*theConstraintHandler;    
    DOF_Numberer 	*theDOF_Numberer;
    AnalysisModel 	*theAnalysisModel;
    EquiSolnAlgo 	*theAlgorithm;
    LinearSOE 		*theSOE;
    TransientIntegrator *theIntegrator;
    int domainStamp;
// AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
	SensitivityAlgorithm *theSensitivityAlgorithm;
#endif
// AddingSensitivity:END ///////////////////////////////
};

#endif

