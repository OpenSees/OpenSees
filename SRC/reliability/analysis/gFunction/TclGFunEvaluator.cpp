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
                                                                        
// $Revision: 1.6 $
// $Date: 2007-02-24 01:21:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/TclGFunEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <TclGFunEvaluator.h>
#include <Vector.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>
#include <RandomVariablePositionerIter.h>
//#include <fstream>
#include <tcl.h>
#include <string.h>


TclGFunEvaluator::TclGFunEvaluator(Tcl_Interp *passedTclInterp,
					ReliabilityDomain *passedReliabilityDomain,
					TCL_Char *passed_fileName)
:GFunEvaluator(passedTclInterp, passedReliabilityDomain)
{
	strcpy(fileName,passed_fileName);
}


TclGFunEvaluator::~TclGFunEvaluator()
{

}


int
TclGFunEvaluator::runGFunAnalysis(const Vector &x)
{	
	// Initial declarations
	char theCommand[100];
	int i;


////// IN CASE AN OPENSEES MODEL EXISTS ////////////////////////////////

	// Zero out the response in the structural domain to make ready for next analysis
	char theRevertToStartCommand[10] = "reset";
	if (Tcl_Eval(theTclInterp, theRevertToStartCommand) == TCL_ERROR) {
	  opserr << "ERROR TclGFunEvaluator -- error in Tcl_Eval for the reset command" << endln;
	  return -1;
	}


	// Put random variables into the structural domain according to the RandomVariablePositioners
	int rvNumber;
	RandomVariablePositionerIter &rvPosIter =
	  theReliabilityDomain->getRandomVariablePositioners();
	RandomVariablePositioner *theRVPos;
	while ((theRVPos = rvPosIter()) != 0) {
	  rvNumber = theRVPos->getRvNumber();
	  theRVPos->update(x(rvNumber-1));
	}

//////////////////////////////////////////////////////////////////////////


	// Set values of random variables in the Tcl intepreter
	for (i=0; i<x.Size(); i++) {
		sprintf(theCommand,"set x_%d %20.12e",(i+1),x(i));
		if (Tcl_Eval(theTclInterp, theCommand) == TCL_ERROR) {
		  opserr << "ERROR TclGFunEvaluator -- error in Tcl_Eval" << endln;
		  return -1;
		}
	}


	// Source the code file that the user has provided
	sprintf(theCommand,"source %s",fileName);
	if (Tcl_Eval(theTclInterp, theCommand) == TCL_ERROR) {
	  opserr << "ERROR TclGFunEvaluator -- error in Tcl_Eval" << endln;
	  return -1;
	}


	return 0;
}




int
TclGFunEvaluator::tokenizeSpecials(TCL_Char *theExpression)
{

	return 0;
}

