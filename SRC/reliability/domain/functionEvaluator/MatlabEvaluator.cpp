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

// $Revision: 1.9 $
// $Date: 2010-09-13 21:39:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/MatlabEvaluator.cpp,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <MatlabEvaluator.h>
#include <Vector.h>
#include <FunctionEvaluator.h>
#include <ReliabilityDomain.h>
#include <Parameter.h>
#include <ParameterPositioner.h>

#include <mex.h>
#include <engine.h>
#include <stdlib.h>
#include <string.h>


MatlabEvaluator::MatlabEvaluator(ReliabilityDomain *passedReliabilityDomain,
						   Domain *passedOpenSeesDomain,
						   TCL_Char *passed_fileName)
:FunctionEvaluator(), theReliabilityDomain(passedReliabilityDomain),
theOpenSeesDomain(passedOpenSeesDomain)

{
	int exprLen = strlen(passed_fileName);
	fileName = new char[exprLen+1];
	strcpy(fileName,passed_fileName);
	
}


MatlabEvaluator::~MatlabEvaluator()
{
	if (theExpression != 0)
		delete [] theExpression;
	if (fileName != 0)
		delete [] fileName;
	
}


int
MatlabEvaluator::setVariables(const Vector &x)
{
	char theIndex[80];
	double xval;
	Parameter *theParam;
	
	// Set values of random variables in the Tcl interpreter
	int nparam = theReliabilityDomain->getNumberOfParameters();
	for (int i = 0; i < nparam; i++) {
		theParam = theReliabilityDomain->getParameterPtrFromIndex(i);
		int paramTag = theParam->getTag();
		
		xval = x(i);
		
		// put in par(1) format
		// NYI for MATLAB
		/*
		sprintf(theIndex,"%d",paramTag);
		if (Tcl_SetVar2Ex(theTclInterp,"par",theIndex,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG) == NULL) {
			opserr << "ERROR MatlabEvaluator -- error in setVariables for parameter tag " << paramTag << endln;
			opserr << "of type " << theTclInterp->result << endln;
			return -1;
		}
		*/
		
	}
	
	return 0;
}


int 
MatlabEvaluator::setExpression(const char *passedExpression)
{
	if (theExpression != 0)
		delete [] theExpression;
	
	int exprLen = strlen(passedExpression);
	theExpression = new char[exprLen+1];
	strcpy(theExpression,passedExpression);
	
	return 0;
}


int 
MatlabEvaluator::addToExpression(const char *in) 
{
	
	return 0;
}


int
MatlabEvaluator::evaluateExpression() 
{
	// NYI for MATLAB
	/*
	if (Tcl_ExprDouble( theTclInterp, theExpression, &current_val) != TCL_OK) {
		opserr << "MatlabEvaluator::evaluateExpression -- expression \"" << theExpression;
		opserr << "\" caused error:" << endln << theTclInterp->result << endln;
		return -1;
	}
	*/
	
	return 0;
}


double
MatlabEvaluator::getResult()
{
	return current_val;
}


int
MatlabEvaluator::runAnalysis(const Vector &x)
{	
	
	// Let's just make a direct call since we have the pointer to OpenSees domain
	// This replaces above call to Tcl command; however, in the reset command
	// revertToStart() is also called on theTransientIntegrator -- MHS needs to check
	if (theOpenSeesDomain->revertToStart() != 0) {
		opserr << "ERROR MatlabEvaluator -- error in resetting Domain" << endln;
		return -1;
	}
	
	// Put random variables into the structural domain according to the RandomVariablePositioners
	int rvIndex;
	RandomVariablePositionerIter rvPosIter = theReliabilityDomain->getRandomVariablePositioners();
	RandomVariablePositioner *theRVPos;
	while ((theRVPos = rvPosIter()) != 0) {
		rvIndex = theRVPos->getRvIndex();
		theRVPos->update(x(rvIndex));
	}
	
	// Start a Matlab engine
	Engine *ep;
	ep = engOpen("\0");
	
	// Execute a Matlab function called 'matlabgfun'
	char theMatlabCommand[50];
	sprintf(theMatlabCommand,"matlabgfun");
	engEvalString(ep, theMatlabCommand);
	
	// Shut down the Matlab engine
	engClose(ep);
	
	return 0;
}

