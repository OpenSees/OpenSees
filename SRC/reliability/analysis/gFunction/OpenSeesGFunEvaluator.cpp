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
                                                                        
// $Revision: 1.24 $
// $Date: 2010-09-13 21:39:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/OpenSeesGFunEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <OpenSeesGFunEvaluator.h>
#include <Vector.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <LimitStateFunctionIter.h>
#include <RandomVariableIter.h>
#include <RandomVariablePositioner.h>
#include <RandomVariablePositionerIter.h>

#include <tcl.h>
#include <stdlib.h>
#include <string.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;

OpenSeesGFunEvaluator::OpenSeesGFunEvaluator(Tcl_Interp *passedTclInterp,
					     ReliabilityDomain *passedReliabilityDomain,
					     Domain *passedOpenSeesDomain,
					     TCL_Char *passedFileName)
  :GFunEvaluator(), theTclInterp(passedTclInterp), theReliabilityDomain(passedReliabilityDomain),
   theOpenSeesDomain(passedOpenSeesDomain), g(0.0)
{
	// here the user has provided a file with the analysis commands
	strcpy(fileName,passedFileName);
	nsteps = 0;
	dt = 0.0;

}

OpenSeesGFunEvaluator::OpenSeesGFunEvaluator(Tcl_Interp *passedTclInterp,
					     ReliabilityDomain *passedReliabilityDomain,
					     Domain *passedOpenSeesDomain,
					     int p_nsteps, double p_dt)
  :GFunEvaluator(), theTclInterp(passedTclInterp), theReliabilityDomain(passedReliabilityDomain),
   theOpenSeesDomain(passedOpenSeesDomain), g(0.0)
{
	// here the user has specified number of steps and possibly dt
	fileName[0] = '\0';
	nsteps = p_nsteps;
	dt = p_dt;

}

OpenSeesGFunEvaluator::~OpenSeesGFunEvaluator()
{
  
}

int
OpenSeesGFunEvaluator::setTclRandomVariables(const Vector &x)
{
  char theIndex[80];
  double xval;
  RandomVariable *theRV;
	
  // Set values of random variables in the Tcl interpreter
  int nrv = theReliabilityDomain->getNumberOfRandomVariables();

  int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();

  for (int i = 0; i < nrv; i++) {
    theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(i);
    int rvTag = theRV->getTag();

    xval = x(i);

    // put in x(1) format
    sprintf(theIndex,"%d",rvTag);
    if (Tcl_SetVar2Ex(theTclInterp,"xrv",theIndex,Tcl_NewDoubleObj(xval),TCL_GLOBAL_ONLY) == NULL) {
      opserr << "ERROR GFunEvaluator -- error in setTclRandomVariables xrv" << endln;
      opserr << theTclInterp->result << endln;
      return -1;
    }
    
    // put in x(1,lsfTag) format (useful for reporting design point)
    sprintf(theIndex,"%d,%d",rvTag,lsf);
    if (Tcl_SetVar2Ex(theTclInterp,"xrv",theIndex,Tcl_NewDoubleObj(xval),TCL_GLOBAL_ONLY) == NULL) {
      opserr << "ERROR GFunEvaluator -- error in setTclRandomVariables xrv" << endln;
      opserr << theTclInterp->result << endln;
      return -1;
    }
    
    // for legacy reasons, also put random variables in x_1 format
    sprintf(theIndex,"x_%d",rvTag);
    if (Tcl_SetVar2Ex(theTclInterp,theIndex,NULL,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG) == NULL) {
      opserr << "ERROR GFunEvaluator -- error in setTclRandomVariables x" << endln;
      opserr << theTclInterp->result << endln;
      return -1;
    }
  }

  return 0;
}

int
OpenSeesGFunEvaluator::evaluateG(const Vector &x) 
{
  g = 0.0;

  // Set random variable values in Tcl namespace
  if (this->setTclRandomVariables(x) != 0) {
    opserr << "ERROR TclGFunEvaluator::evaluateG -- error in setTclRandomVariables" << endln;
    return -1;
  }

  // "Download" limit-state function from reliability domain
  int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
  LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
  
  // Get the limit-state function expression
  const char *theExpression = theLimitStateFunction->getExpression();

  if (Tcl_ExprDouble( theTclInterp, theExpression, &g) != TCL_OK) {
    opserr << "OpenSeesGFunEvaluator::evaluateGnoRecompute -- expression \"" << theExpression;
    opserr << "\" caused error:" << endln << theTclInterp->result << endln;
    return -1;
  }

  numberOfEvaluations++;

  return 0;
}

double
OpenSeesGFunEvaluator::getG(void)
{
  return g;
}


/*
int
OpenSeesGFunEvaluator::runGFunAnalysis(const Vector &x)
{
	// Zero out the response in the structural domain to make ready for next analysis
	char theRevertToStartCommand[10] = "reset";
	if (Tcl_Eval(theTclInterp, theRevertToStartCommand) == TCL_ERROR) {
		opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_Eval for the reset command" << endln;
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
	
	// Set random variable values in Tcl namespace, if needed
	if (setTclRandomVariables(x) != 0) {
		opserr << "ERROR OpenSeesGFunEvaluator::runGFunAnalysis -- error in setTclRandomVariables" << endln;
		return -1;
	}

	// Run the structural analysis according to user specified scheme
	double result = 0;
	if (dt==0.0 && nsteps==0 && strlen(fileName) == 0) {
		// Run up to max time in fFuncs
		opserr << "OpenSeesGFunEvaluator: The option -runToMaxTimeInGFun " << endln
			<< " is not yet implemented." << endln;
	}
	else if (dt==0.0 && nsteps==0) {
		// Read commands from file and execute them
		if (Tcl_EvalFile(theTclInterp, fileName) == TCL_ERROR) {
			opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_EvalFile" << endln;
			return -1;
		}
		
		// (The reason for doing it like this is that we want to check
		// the flag returned by the analysis to see if it converged).
		//char theAnalyzeCommand[30];
		//sprintf(theAnalyzeCommand,"[source %s]",fileName);
		//if (Tcl_ExprDouble( theTclInterp, theAnalyzeCommand, &result) == TCL_ERROR) {
		//  opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_ExprDouble for the analyze file command" << endln;
		//  return -1;
		//}
		// ---- NOTE to writer of above, if you are really serious about checking the return value for convergence, you 
		// need to query theTclInterp->result, Tcl_ExprDouble returns the standard TCL codes, not the script return value

	}
	else {
		// User has given "nsteps" and possibly "dt"
		// NOTE: whoever wrote this should rewrite using OpenSees methods, not writing to Tcl interpreter
		char theAnalyzeCommand[30];
		if (dt == 0.0) {
			sprintf(theAnalyzeCommand,"[analyze %d]",nsteps);
			
			if (Tcl_ExprDouble( theTclInterp, theAnalyzeCommand, &result) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_ExprDouble for the analyze command" << endln;
			  return -1;
			}
		}
		else {
			sprintf(theAnalyzeCommand,"[analyze %d %10.5f]", nsteps, dt);
			if (Tcl_ExprDouble( theTclInterp, theAnalyzeCommand, &result) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_ExprDouble for the analyze command" << endln;
			  return -1;
			}
			
		}
	
	}

	return ((int)result);
}




int
OpenSeesGFunEvaluator::tokenizeSpecials(TCL_Char *theExpression, Tcl_Obj *passedList)
{
	// Set value of OpenSees finite element response quantities 
	// appearing in the limit-state function in the Tcl domain
	
	// Note: we no longer need theExpression as all the parameters were saved in LimitStateFunction
	
	Tcl_Obj *paramList = Tcl_DuplicateObj(passedList);
	int llength;
	Tcl_Obj *objPtr;
	if (Tcl_ListObjLength(theTclInterp,paramList,&llength) != TCL_OK) {
		opserr << "OpenSeesGFunEvaluator::tokenizeSpecials ERROR getting list length. " << endln;
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


void
OpenSeesGFunEvaluator::setNsteps(int p_nsteps)
{
	nsteps = p_nsteps;
}

double
OpenSeesGFunEvaluator::getDt()
{
	return dt;
}


double OpenSeesGFunEvaluator::getG2(double g, double littleDeltaT)
{

	// Initial declaractions
	double perturbationFactor = 0.001; // (is multiplied by stdv and added to others...)
	double g_perturbed;
	char varName[20] = "";
	char arrName[40] = "";
	double g2FunctionValue = g;
	
	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
	char *theExpression = theLimitStateFunction->getExpression();
	Tcl_Obj *passedList = theLimitStateFunction->getParameters();
	
	// Cycle over all the LSF parameters and compute gradients
	Tcl_Obj *paramList = Tcl_DuplicateObj(passedList);
	int llength = 0;
	Tcl_Obj *objPtr;
	Tcl_ListObjLength(theTclInterp,paramList,&llength);
	
	for (int jk = 0; jk < llength; jk++) {
		Tcl_ListObjIndex(theTclInterp,paramList,jk,&objPtr);
		char *listStr = Tcl_GetStringFromObj(objPtr,NULL);

		// If a RV is detected
		if ( strncmp(listStr, "x_", 2) == 0 || strncmp(listStr, "xrv", 3) == 0 ) {
			opserr<<"Fatal: can not deal with x in limit state function"<<endln;
			exit(-1);
		}
		
		// If a node quantity is detected
		else if ( strncmp(listStr, "u", 1) == 0 || strncmp(listStr, "node", 4) == 0 ||
			 strncmp(listStr, "rec_node", 8) == 0 ) {
			 
			// Get node number, dof number, and type of response
			// take advantage of GFunEvaluator parser, do NOT recode everything here
			int nodeNumber, direction;
			char dispOrWhat[10] = "";
			if ( strncmp(listStr, "u", 1) == 0 )
				uParse(listStr, &nodeNumber, &direction, dispOrWhat, varName, arrName);
			else if ( strncmp(listStr, "node", 4) == 0 || strncmp(listStr, "rec_node", 8) == 0 )
				nodeParse(listStr, &nodeNumber, &direction, dispOrWhat, varName, arrName);
				
			if ( strncmp(dispOrWhat, "disp", 4) != 0) {
				opserr<<"Fatal: can not deal with ud in limit state function"<<endln; 
				exit(-1);
			}
			
			// If a nodal displacement is detected
			// Keep the original value
			double originalValue = 0.0;
			Tcl_Obj *tempDouble = Tcl_GetVar2Ex(theTclInterp, varName, arrName, TCL_LEAVE_ERR_MSG);
			if (tempDouble == NULL) {
				opserr << "ERROR OpenSeesGFunEvaluator -- GetVar error" << endln;
				opserr << theTclInterp->result << endln;
				return -1;
			}
			Tcl_GetDoubleFromObj(theTclInterp, tempDouble, &originalValue);
			
			// Set perturbed value in the Tcl workspace
			// ---------------------- Quan and michele  
			double newValue;
			//if (originalValue ==0 ) newValue = 0.0001;
			if (fabs(originalValue) < 1.0e-14 )
				newValue = 0.0001;
			else
				newValue = originalValue*(1.0+perturbationFactor);
			
			// set Tcl value
			if (Tcl_SetVar2Ex(theTclInterp,varName,arrName,Tcl_NewDoubleObj(newValue),TCL_LEAVE_ERR_MSG) == NULL) {
				opserr << "ERROR OpenSeesGfunEvaluator -- SetVar error" << endln;
				opserr << theTclInterp->result << endln;
				return -1;
			}
			
			// Evaluate the limit-state function again
			if (evaluateGnoRecompute(theExpression) < 0) {
				opserr << "ERROR OpenSeesGradGEvaluator -- error evaluating LSF" << endln;
				opserr << theTclInterp->result << endln;
				return -1;
			}
			g_perturbed = getG();

			// Compute gradient
			double onedgdu = (g_perturbed-g)/(newValue - originalValue);

			// Make assignment back to its original value
			nodeTclVariable(nodeNumber, direction, dispOrWhat, varName, arrName);

			// get nodal velocity from the domain
			double udot = 0.0;
			nodeTclVariable(nodeNumber, direction, "vel", "udot", "getG2");
			
			tempDouble = Tcl_GetVar2Ex(theTclInterp, "udot", "getG2", TCL_LEAVE_ERR_MSG);
			if (tempDouble == NULL) {
				opserr << "warning:OpenSeesGFunEvaluator::getSecondG() velocity is 0! probabily because you are running static case" << endln;
				opserr << theTclInterp->result << endln;
				return -1;
			}
			Tcl_GetDoubleFromObj(theTclInterp, tempDouble, &udot);
			
			g2FunctionValue += onedgdu*udot*littleDeltaT;

		}

	} 

	// reclaim Tcl object space
	Tcl_DecrRefCount(paramList);
	
	return g2FunctionValue;

}
*/
