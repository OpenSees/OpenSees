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
                                                                        
// $Revision: 1.7 $
// $Date: 2010-09-13 21:37:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/GradientEvaluator.cpp,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <GradientEvaluator.h>
#include <ReliabilityDomain.h>

GradientEvaluator::GradientEvaluator(ReliabilityDomain *passedReliabilityDomain,
				     FunctionEvaluator *passedGFunEvaluator)
{

  theReliabilityDomain = passedReliabilityDomain;
  theFunctionEvaluator = passedGFunEvaluator;
  
  /////S added by K Fujimura /////
  finitedifference = false;
  numberOfEvalIncSens = 0;
  /////E added by K Fujimura /////
}


GradientEvaluator::~GradientEvaluator()
{

}


int 
GradientEvaluator::computeParameterDerivatives(double g)
{
  /*
	// Zero out the previous result matrix
	if (DgDpar != 0) {
		delete DgDpar;
		DgDpar = 0;
	}

	// Initial declarations
	double g_perturbed;
	double onedgdpar;
	int llength = 0;

	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
	char *theExpression = theLimitStateFunction->getExpression();
	Tcl_Obj *passedList = theLimitStateFunction->getParameters();
	
	// Cycle over all the LSF parameters and compute gradients
	Tcl_Obj *paramList = Tcl_DuplicateObj(passedList);
	Tcl_Obj *objPtr;
	Tcl_ListObjLength(theTclInterp,paramList,&llength);
		
	for (int jk = 0; jk < llength; jk++) {
		Tcl_ListObjIndex(theTclInterp,paramList,jk,&objPtr);
		char *listStr = Tcl_GetStringFromObj(objPtr,NULL);

		// If a RV is detected
		if ( strncmp(listStr, "par", 3) == 0 ) {

			// Get parameter number
			int parameterNumber;
			int args = sscanf(listStr,"par(%i)",&parameterNumber);
			if (args != 1) {
				opserr << "GradientEvaluator ERROR -- could not determine parameter number in LSF " << lsf << endln;
				return -1;
			}
			char theIndex[20];
			sprintf(theIndex,"%d",parameterNumber);

			// Store the original parameter value
			double originalValue;
			Tcl_Obj *tempDouble = Tcl_GetVar2Ex(theTclInterp, "par", theIndex, TCL_LEAVE_ERR_MSG);
			if (tempDouble == NULL) {
				opserr << "ERROR GradientEvaluator -- GetVar error" << endln;
				opserr << theTclInterp->result << endln;
				return -1;
			}
			Tcl_GetDoubleFromObj(theTclInterp, tempDouble, &originalValue);
			
			// Assign a perturbed value
			double xval = 0.001;
			if ( fabs(originalValue) > 2.0*DBL_EPSILON )
				xval = originalValue * 1.001;
				
			Tcl_Obj *outp = Tcl_SetVar2Ex(theTclInterp,"par",theIndex,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG);
			if (outp == NULL) {
				opserr << "ERROR GradientEvaluator -- error setting parameter with tag " << parameterNumber << endln;
				opserr << theTclInterp->result << endln;
				return -1;
			}
			
			// Evaluate limit-state function again
			if (theGFunEvaluator->evaluateGnoRecompute(theExpression) < 0) {
				opserr << "ERROR GradientEvaluator -- error evaluating LSF" << endln;
				opserr << theTclInterp->result << endln;
				return -1;
			}
			g_perturbed = theGFunEvaluator->getG();
			
			// Compute the gradient 'dgdpar' by finite difference
			if ( fabs(originalValue) > 2.0*DBL_EPSILON )
				onedgdpar = (g_perturbed-g)/(originalValue*0.001);
			else
				onedgdpar = (g_perturbed-g)/0.001;

			// Make assignment back to its original value
			outp = Tcl_SetVar2Ex(theTclInterp,"par",theIndex,Tcl_NewDoubleObj(originalValue),TCL_LEAVE_ERR_MSG);
			if (outp == NULL) {
				opserr << "ERROR GradientEvaluator -- error setting parameter with tag " << parameterNumber << endln;
				opserr << theTclInterp->result << endln;
				return -1;
			}
			
			// Store the DgDpar in a matrix (make it expand successively)
			if (DgDpar == 0) {
				DgDpar = new Matrix(1, 2);
				(*DgDpar)(0,0) = (double)parameterNumber;
				(*DgDpar)(0,1) = onedgdpar;
			}
			else {
				int oldSize = DgDpar->noRows();
				Matrix tempMatrix = *DgDpar;
				delete DgDpar;
				DgDpar = new Matrix(oldSize+1, 2);
				for (int i=0; i<oldSize; i++) {
					(*DgDpar)(i,0) = tempMatrix(i,0);
					(*DgDpar)(i,1) = tempMatrix(i,1);
				}
				(*DgDpar)(oldSize,0) = (double)parameterNumber;
				(*DgDpar)(oldSize,1) = onedgdpar;
			}

		}

	}
	
	// reclaim Tcl object space
	Tcl_DecrRefCount(paramList);
  */
	return 0;
}


int
GradientEvaluator::initializeNumberOfEvaluations()
{
  numberOfEvalIncSens = 0;
  return 0;
}


int
GradientEvaluator::getNumberOfEvaluations()
{
  return numberOfEvalIncSens;
}


void 
GradientEvaluator::setPerformFuncCoeffs(TaggedObjectStorage* a) 
{
  opserr << "GFunEvaluator::setPerformFuncCoeffs() -- This method is not " << endln
	 << " implemented for the chosen type of GradgEvaluator." << endln;
}

void 
GradientEvaluator::setPerformFuncCoeffIter(PerformanceFunctionCoefficientIter* a)
{
  opserr << "GFunEvaluator::setPerformFuncCoeffIter() -- This method is not " << endln
	 << " implemented for the chosen type of GradgEvaluator." << endln;
}
