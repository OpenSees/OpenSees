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
**   Optimization module developed by:                                **
**   Quan Gu  (qgu@ucsd.edu)                                          **
**   Joel Conte (jpconte@ucsd.edu)                                    **
**   Philip Gill (pgill@ucsd.edu)                                     **
** ****************************************************************** */


//
// Written by Quan Gu (qgu@ucsd.edu)    March 2010
//


#ifndef SNOPTANALYSIS
#define SNOPTANALYSIS



#include <OptimizationDomain.h>
//#include "snoptfilewrapper.h"
#include "toyOptfunction.h"
#include "SNOPTclass.h"
#include "tcl.h"

#include <DesignVariable.h>
#include <DesignVariablePositioner.h>
#include <ConstraintFunction.h>
#include <ObjectiveFunction.h>

class SNOPTAnalysis 	: public SNOPTCLASS  
{
public:

	SNOPTAnalysis(int passedMaxNumberOfIterations, 
					int pPrintFlag,
					char *pFileNamePrint,
					char * probType, 
					OptimizationDomain * passedOptimizationDomain,
					Tcl_Interp *passedTclInterp,
					char * pTclFileToRun
					);
	virtual ~SNOPTAnalysis();




// -------------- add function to realize interface with opensees -------------

public:
  int runOptAnalysis (OptimizationDomain  * passedOptimizationDomain);
  double * getUnScaledX();
  double * getScaledX();
  int getNumberOfSteps();
  int getNumberOfEvaluations();

 //------------------------------------------------------------------------------
  double * getScaledFunctionPtr(){return F;}; // scales means all dv: x/x_0  //  |  these functions are called by snopt body only.
  double * getScaledGradientPtr(){return gradient;};                         //  |
                                                                             //  |  
  int updateXFG(double * newScaledX);   // update x and recompute F, G (if needG)|   
//-------------------------------------------------------------------------------



private:


// The domain and tools for the analysis
	OptimizationDomain *theOptimizationDomain;
	ObjectiveFunction *theObjectiveFunction;
	ConstraintFunction *theConstraintFunction;

	double * gradient;
	double valueOfObjectiveFuction;

	int maxNumberOfIterations;
	int ii; //number of steps

	int printFlag;
	char *fileNamePrint;
	int numberOfEvaluations;
	double * scales;


	bool needG; // show whether  gradient is provided
	
	double * temp;
	char * tclFileToRun;
	Tcl_Interp * theTclInterp;


};

#endif 
