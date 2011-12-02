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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-05-27 20:04:30 $
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

#include <tcl.h>
#include <stdlib.h>
#include <string.h>


TclGFunEvaluator::TclGFunEvaluator(Tcl_Interp *passedTclInterp,
					ReliabilityDomain *passedReliabilityDomain,
					Domain *passedOpenSeesDomain,
					TCL_Char *passed_fileName)
:GFunEvaluator(passedTclInterp, passedReliabilityDomain, passedOpenSeesDomain)
{
	strcpy(fileName,passed_fileName);
}


TclGFunEvaluator::~TclGFunEvaluator()
{

}


int
TclGFunEvaluator::runGFunAnalysis(const Vector &x)
{	

////// IN CASE AN OPENSEES MODEL EXISTS ////////////////////////////////

	// Zero out the response in the structural domain to make ready for next analysis
	char theRevertToStartCommand[10] = "reset";
	if (Tcl_Eval(theTclInterp, theRevertToStartCommand) == TCL_ERROR) {
	  opserr << "ERROR TclGFunEvaluator -- error in Tcl_Eval for the reset command" << endln;
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

//////////////////////////////////////////////////////////////////////////

	// Set random variable values in Tcl namespace
	if (setTclRandomVariables(x) != 0) {
		opserr << "ERROR TclGFunEvaluator::runGFunAnalysis -- error in setTclRandomVariables" << endln;
		return -1;
	}


	// Source the code file that the user has provided
	if (Tcl_EvalFile(theTclInterp, fileName) == TCL_ERROR) {
		opserr << "ERROR TclGFunEvaluator -- error in Tcl_EvalFile" << endln;
		return -1;
	}


	return 0;
}


int
TclGFunEvaluator::tokenizeSpecials(TCL_Char *theExpression, Tcl_Obj *passedList)
{
	// Set value of OpenSees finite element response quantities 
	// appearing in the limit-state function in the Tcl domain
	
	// Note: we no longer need theExpression as all the parameters were saved in LimitStateFunction
	
	Tcl_Obj *paramList = Tcl_DuplicateObj(passedList);
	int llength;
	Tcl_Obj *objPtr;
	if (Tcl_ListObjLength(theTclInterp,paramList,&llength) != TCL_OK) {
		opserr << "TclGFunEvaluator::tokenizeSpecials ERROR getting list length. " << endln;
		opserr << theTclInterp->result << endln;
	}
	
	// initialize variables that will be used
	char dispOrWhat[10] = "";
	char eleRest[100] = "";
	char varName[20] = "";
	char arrName[40] = "";
	int nodeNumber = 0, direction = 0, eleNumber = 0;
	
	for (int jk = 0; jk < llength; jk++) {
		Tcl_ListObjIndex(theTclInterp,paramList,jk,&objPtr);
		char *listStr = Tcl_GetStringFromObj(objPtr,NULL);

		// If a node quantity is detected
		if ( strncmp(listStr, "u", 1) == 0 ) {
			uParse(listStr, &nodeNumber, &direction, dispOrWhat, varName, arrName);
			nodeTclVariable(nodeNumber, direction, dispOrWhat, varName, arrName);
		}
		else if ( strncmp(listStr, "node", 4) == 0 || strncmp(listStr, "rec_node", 8) == 0 ) {
			nodeParse(listStr, &nodeNumber, &direction, dispOrWhat, varName, arrName);
			nodeTclVariable(nodeNumber, direction, dispOrWhat, varName, arrName);
		}      
		// If an element quantity is detected
		else if ( strncmp(listStr, "element", 7) == 0 || strncmp(listStr, "rec_element", 11) == 0 ) {
			elementParse(listStr, &eleNumber, varName, eleRest);
			elementTclVariable(eleNumber, varName, eleRest);
		}
	}
	
	// reclaim Tcl object space
	Tcl_DecrRefCount(paramList);
	
	return 0;
}

