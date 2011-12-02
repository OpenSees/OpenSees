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
**   Philip Gill (pgill@ucsd.edu)                                     **                                 **
** ****************************************************************** */
                                                                        

//
// Written by Quan Gu (qgu@ucsd.edu)
//

#ifndef SNOPTANALYSIS_H
#define SNOPTANALYSIS_H


#include <SnoptWrapper.h>
#include <ReliabilityDomain.h>

class SNOPTAnalysis : public SnoptWrapper 
		      //class SNOPTAnalysis 
{

public:

  SNOPTAnalysis(int passedMaxNumberOfIterations, 
		int pPrintFlag,
		char *pFileNamePrint,
		char * probType, 
		ReliabilityDomain * passedReliabilityDomain,
		Tcl_Interp *passedTclInterp,
		char * pTclFileToRun
		);

  virtual ~SNOPTAnalysis();




// -------------- add function to realize interface with opensees -------------

public:
  int runOptAnalysis (ReliabilityDomain *passedReliabilityDomain);
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


// The reliability domain and tools for the analysis
	ReliabilityDomain *theReliabilityDomain;
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
