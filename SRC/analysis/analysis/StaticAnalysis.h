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
                                                                        
// $Revision: 1.8 $
// $Date: 2009-05-11 21:32:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/StaticAnalysis.h,v $
                                                                        
                                                                        
#ifndef StaticAnalysis_h
#define StaticAnalysis_h

// Written: fmk 
// Revision: A
//
// Description: This file contains the interface for the StaticAnalysis
// class. StaticAnalysis is a subclass of Analysis, it is used to perform 
// a static analysis on the FE\_Model.
//
// What: "@(#) StaticAnalysis.h, revA"

#include <Analysis.h>


// AddingSensitivity:BEGIN //////////////////////////////////
#ifdef _RELIABILITY
#include <SensitivityAlgorithm.h>
#include<Integrator.h>
#include<IncrementalIntegrator.h>
#endif
// AddingSensitivity:END ////////////////////////////////////


class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class StaticIntegrator;
class LinearSOE;
class EquiSolnAlgo;
class ConvergenceTest;
class EigenSOE;

class StaticAnalysis: public Analysis
{
  public:
        StaticAnalysis(Domain &theDomain,
		   ConstraintHandler &theHandler,
		   DOF_Numberer &theNumberer,
		   AnalysisModel &theModel,
		   EquiSolnAlgo &theSolnAlgo,		   
		   LinearSOE &theSOE,
		   StaticIntegrator &theIntegrator,
		   ConvergenceTest *theTest = 0);

    ~StaticAnalysis();

    void clearAll(void);	    
    
    int analyze(int numSteps);
    int eigen(int numMode, bool generlzed = true, bool findSmallest = true);
    int initialize(void);
    int domainChanged(void);

    int setNumberer(DOF_Numberer &theNumberer);
    int setAlgorithm(EquiSolnAlgo &theAlgorithm);
    int setIntegrator(StaticIntegrator &theIntegrator);
    int setLinearSOE(LinearSOE &theSOE);
    int setConvergenceTest(ConvergenceTest &theTest);
    int setEigenSOE(EigenSOE &theSOE);

    EquiSolnAlgo     *getAlgorithm(void);
    StaticIntegrator *getIntegrator(void);
    ConvergenceTest  *getConvergenceTest(void);

    // AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
    int setSensitivityAlgorithm(/*SensitivityAlgorithm*/Integrator *theSensitivityAlgorithm);
#endif
    // AddingSensitivity:END /////////////////////////////////
    
  protected: 
    
  private:
    ConstraintHandler 	*theConstraintHandler;    
    DOF_Numberer 	*theDOF_Numberer;
    AnalysisModel 	*theAnalysisModel;
    EquiSolnAlgo 	*theAlgorithm;
    LinearSOE 		*theSOE;
    EigenSOE 		*theEigenSOE;
    StaticIntegrator    *theIntegrator;
    ConvergenceTest     *theTest;
    int domainStamp;

#ifdef _RELIABILITY

#endif
};

#endif
