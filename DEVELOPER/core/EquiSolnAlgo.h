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
// $Date: 2008-08-26 17:07:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/EquiSolnAlgo.h,v $
                                                                        
                                                                        
#ifndef EquiSolnAlgo_h
#define EquiSolnAlgo_h

// File: ~/OOP/analysis/algorithm/EquiSolnAlgo.h 
// 
// Written: fmk 
// Created: 11/96 
// Revision: A 
//

// Description: This file contains the class definition for 
// EquiSolnAlgo. EquiSolnAlgo is an abstract base class, 
// i.e. no objects of it's type can be created.  Its subclasses deifine
// the sequence of operations to be performed in the analysis by static
// equilibrium of a finite element model.  
// 
// What: "@(#)EquiSolnAlgo.h, revA"

#include <SolutionAlgorithm.h>
#include <IncrementalIntegrator.h>

class AnalysisModel;
class LinearSOE;
class ConvergenceTest;

class EquiSolnAlgo: public SolutionAlgorithm
{
  public:
    EquiSolnAlgo(int classTag);
    virtual ~EquiSolnAlgo();

    // public functions defined for subclasses
    virtual void setLinks(AnalysisModel &theModel, 
			  IncrementalIntegrator &theIntegrator,
			  LinearSOE &theSOE,
			  ConvergenceTest *theTest);
    
    // virtual functions
    virtual int solveCurrentStep(void) =0;
    virtual int setConvergenceTest(ConvergenceTest *theNewTest);    
    virtual ConvergenceTest *getConvergenceTest(void);     

    virtual void Print(OPS_Stream &s, int flag =0) =0;    

    virtual int getNumFactorizations(void) {return 0;}
    virtual int getNumIterations(void) {return 0;}
    virtual double getTotalTimeCPU(void)   {return 0.0;}
    virtual double getTotalTimeReal(void)  {return 0.0;}
    virtual double getSolveTimeCPU(void)   {return 0.0;}
    virtual double getSolveTimeReal(void)  {return 0.0;}
    virtual double getAccelTimeCPU(void)   {return 0.0;}
    virtual double getAccelTimeReal(void)  {return 0.0;}
 
    // the following are not protected as convergence test
    // may need access to them
    AnalysisModel           *getAnalysisModelPtr(void) const;
    IncrementalIntegrator   *getIncrementalIntegratorPtr(void) const;
    LinearSOE	            *getLinearSOEptr(void) const;

  protected:
    ConvergenceTest *theTest;
    
  private:
    AnalysisModel 	  *theModel;
    IncrementalIntegrator *theIntegrator;
    LinearSOE 		  *theSysOfEqn;
};

#endif


