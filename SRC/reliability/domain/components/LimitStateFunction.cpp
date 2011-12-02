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
                                                                        
// $Revision: 1.17 $
// $Date: 2010-09-13 21:33:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/LimitStateFunction.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <LimitStateFunction.h>
#include <Vector.h>
#include <string.h>
#include <stdlib.h>
#include <classTags.h>
#include <tcl.h>

LimitStateFunction::LimitStateFunction(	int passedTag, 
					TCL_Char *passedExpression, Tcl_Interp *passedTclInterp)
  :ReliabilityDomainComponent(passedTag, LIMIT_STATE_FUNCTION)
{
	
	theTclInterp = passedTclInterp;
	
	int exprLen = strlen(passedExpression);
	originalExpression = new char[exprLen+1];
	strcpy(originalExpression,passedExpression);
	
	expressionWithAddition = new char[exprLen+1];
	strcpy(expressionWithAddition,passedExpression);
	
	tokenizeIt(originalExpression);

	this->initializeFORMAnalysis();
	this->initializeSimulationAnalysis();
	this->initializeSORMAnalysis();
}

/*
LimitStateFunction::LimitStateFunction(int tag, int classTag)
  :ReliabilityDomainComponent(tag, classTag)
{

}
*/

LimitStateFunction::~LimitStateFunction()
{
	if (originalExpression != 0)
		delete [] originalExpression;
	if (expressionWithAddition != 0)
		delete [] expressionWithAddition;
	
	// reclaim Tcl object space
	Tcl_DecrRefCount(paramList);
}

void
LimitStateFunction::initializeFORMAnalysis(void)
{
  GFunValueAtStartPt = 0.0;
  GFunValueAtEndPt = 0.0;
  FORMReliabilityIndexBeta = 0.0;
  FORMProbabilityOfFailure_pf1 = 0.0;
  numberOfStepsToFindDesignPointAlgorithm = 0;
}

void
LimitStateFunction::initializeSimulationAnalysis(void)
{
  SimulationReliabilityIndexBeta = 0.0;
  SimulationProbabilityOfFailure_pfsim = 0.0;
  CoefficientOfVariationOfPfFromSimulation = 0.0;
  NumberOfSimulations = 0;
}

void
LimitStateFunction::initializeSORMAnalysis(void)
{
  SORMCurvatureFittingBetaBreitung = 0.0;
  SORMCurvatureFittingPf2Breitung = 0.0;
  SORMPointFittingBetaBreitung = 0.0;
  SORMPointFittingPf2Breitung = 0.0;
  SORMUsingSearchBetaBreitung = 0.0;
  SORMUsingSearchPf2Breitung = 0.0;
  numberOfCurvaturesUsed = 0;
}


void
LimitStateFunction::Print(OPS_Stream &s, int flag)  
{
	s << "Limit State Function #" << this->getTag() << endln;
	s << "Expression: " << this->getExpression() << endln;
	s << endln;
}



char *
LimitStateFunction::getExpression()
{
	return expressionWithAddition;
}


char *
LimitStateFunction::getTokenizedExpression()
{
	return expressionWithAddition;
}

Tcl_Obj *
LimitStateFunction::getParameters()
{
	return paramList;
}

int
LimitStateFunction::addExpression(char *addition)
{
	int exprLen = strlen(originalExpression) + strlen(addition);
	delete [] expressionWithAddition;
	expressionWithAddition = new char[exprLen+1];
	
	strcpy(expressionWithAddition,originalExpression);
	strcat(expressionWithAddition,addition);

	tokenizeIt(expressionWithAddition);

	return 0;
}

int
LimitStateFunction::removeAddedExpression()
{
	int exprLen = strlen(originalExpression);
	delete [] expressionWithAddition;
	expressionWithAddition = new char[exprLen+1];
	
	strcpy(expressionWithAddition,originalExpression);

	tokenizeIt(expressionWithAddition);

	return 0;
}

int
LimitStateFunction::addGradientExpression(const char *expression, int rvTag)
{
  map<int, string>::iterator theExpr;

  this->removeGradientExpression(rvTag);

  // Check if the expression is already in map
  theExpr = mapOfGradientExpressions.find(rvTag);
  
  // Not there, so add
  if (theExpr == mapOfGradientExpressions.end()) {
    mapOfGradientExpressions.insert(map<int,string>::value_type(rvTag,expression));
    
    // Check if successful
    theExpr = mapOfGradientExpressions.find(rvTag);
    if (theExpr == mapOfGradientExpressions.end()) {
      opserr << "LimitStateFunction::addGradientExpression -- map STL failed to add object with tag: "
	     << rvTag << endln;
      return -1;
    }
  }
  // Already there, give error
  else {
    opserr << "LimitStateFunction::addGradientExpression -- object with tag "
	   << rvTag << " already exists" << endln;    
    return -1;
  }
  
  return 0;
}

int
LimitStateFunction::removeGradientExpression(int rvTag)
{
  map<int, string>::iterator theExpr;

  // Check if the expression is already in map
  theExpr = mapOfGradientExpressions.find(rvTag);

  // If not there, do nothing
  if (theExpr == mapOfGradientExpressions.end())
    return 0;
  // Already there, so remove
  else {
    int ok = mapOfGradientExpressions.erase(rvTag);
    if (ok != 1) {
      opserr << "LimitStateFunction::removeGradientExpression -- map STL failed to remove object with tag: "
	     << rvTag << endln;
      return -1;
    }
  }

  return 0;
}

const char* 
LimitStateFunction::getGradientExpression(int rvTag) 
{
  map<int, string>::iterator theExpr;

  theExpr = mapOfGradientExpressions.find(rvTag);
  if (theExpr == mapOfGradientExpressions.end())
    return 0;
  else
    return (*theExpr).second.c_str();
}

int
LimitStateFunction::tokenizeIt(const char *originalExpression)
{
	// Can automatically convert Terje's {} stuff to proper Tcl syntax
	// in this method, e.g., {x_1} --> \$x(1) and {u_5_2} --> \[nodeDisp 5 2]
	//
	// implementation of any new patterns needs to also be reflected in gfunction and sensitivity algorithm

	char separators[] = "}{";
	int deprecated = strcspn(originalExpression, separators);
	int originalLen = strlen(originalExpression);
	
	if (deprecated < originalLen) {
		opserr << "WARNING: Limit state function " << this->getTag() << " contains OLD syntax "
				<< "that uses {x_1}, {u_1_1}, {file_fileName_1_1}, etc. " << endln;
		opserr << "Use new Tcl variable syntax \\$x(1), \\$u(1,1), etc." << endln << endln;
		opserr << "Just a note (from KRM) about the deprecated LSF syntax:\n"
				<< "you can still use RVs in your LSF, however, they are now specified directly:\n"
				<< "	{x_2}  -->  \\$xrv(2)\n"
				<< "if you feel really nostalgic, you can still write the following:\n"
				<< "	{x_2}  -->  \\$x_2\n"
				<< "however the second option will go away at some point.\n" << endln
				<< "Nodal commands now work as follows:\n"
				<< "	{u_1_2}  -->  \\$u(1,2)\n"
				<< "	{u_1_2}  -->  \\[nodeDisp 1 2]\n"
				<< "	{u_1_2}  -->  \\$node(1,2,disp) or \\$rec_node(1,2,disp)\n"
				<< "	{ud_1_2} -->  \\$rec_node(1,2,vel) or \\$ud(1,2)\n"
				<< "	etc., etc.\n" << endln
				<< "Element forces can be defined in the same way and are element dependent:\n"
				<< "	\\$element(2,localForce_1) would give you column 1 of ele 2 localForce\n"
				<< "	\\$rec_element(1,section_2_force_1), etc, etc.\n" << endln
				<< "Also, the code in this file also allowed reading of values from a file.\n"
				<< "This too can be implemented more directly by the user simply by using \n"
				<< "standard Tcl variables.  ie., instead of computing a value, writing it \n"
				<< "to a file, and then using file_fileName_1_2 or whatever, now:\n"
				<< "	{file_fileName_1_2}  -->  \\$yourTclVariable\n" << endln
				<< "Limit state function parameters can now be specified and used, both for \n"
				<< "sensitivity studies, and for parametricReliabilityAnalysis:\n"
				<< "	{par_1}  -->  \\$par(1)\n" << endln
				<< "And a final note, if a particular GFunEvaluator is still obsessed with\n"
				<< "case-specific LSF arguments, they can still be implemented through \n"
				<< "the tokenizeSpecials() function. However, they must be specified in \n"
				<< "standard Tcl syntax and the values initialized in your own script.\n" << endln
				<< "If you are using DDM, the \\[nodeDisp 1 2] or similar commands will not \n"
				<< "work because there is no way to perturb these to calculate dg/du. So \n"
				<< "if you are doing DDM, use the \\$u(1,2), \\$ud(1,2), etc syntax.  For \n"
				<< "DDM, you can also use element forces (if sensitivities are coded).\n" << endln;
		exit(-1);
	}
	
	// now match regular expressions of known lsf parameters
	// Note: lsf parameters are not the same as the Parameter class in OS domain.
	const int numberOfParamTypes = 5;
	char* pattern[numberOfParamTypes] = {"x(_|rv\\()[0-9]+\\)?",
										 "ud*\\([0-9]+,[0-9]+\\)",
										 "(rec_)?node\\([0-9]+,[0-9]+,[a-z]+\\)",
										 "(rec_)?element\\([0-9]+,[a-zA-Z0-9_]+\\)",
										 "par\\([0-9]+\\)"};
	const char* first;
	const char* last;
	char par_result[30];
	
	// initialize the tcl list object
	paramList = Tcl_NewListObj(0,NULL);
	
	// cycle through pattern types and save instances of parameters to tcl list object
	for (int ireg = 0; ireg < numberOfParamTypes; ireg++) {
		Tcl_RegExp regexp = Tcl_RegExpCompile(theTclInterp, pattern[ireg]);
		if (regexp == NULL) {
			opserr << "LSF::tokenizeIt ERROR compiling regular expression -- " << pattern[ireg] << endln;
			opserr << theTclInterp->result << endln;
		}
		
		char* current = new char[originalLen+1];
		strcpy(current,originalExpression);
			
		while ( Tcl_RegExpExec(theTclInterp, regexp, current, current) == 1 ) {
			// match found
			Tcl_RegExpRange(regexp, 0, &first, &last);
			
			if (first) {
				strncpy(par_result,first,last-first);
				par_result[last-first] = '\0';
				
				//opserr << "Found: " << par_result << endln;
				Tcl_Obj *tempStr = Tcl_NewStringObj(par_result,last-first);
				//opserr << Tcl_GetStringFromObj(tempStr,NULL) << endln;
				
				if (Tcl_ListObjAppendElement(theTclInterp, paramList, tempStr) != TCL_OK) {
					opserr << "LSF::tokenizeIt ERROR creating list element from " << par_result << endln;
					opserr << theTclInterp->result << endln;
				}
			}
			
			strcpy(current, last);
		} 
			
		delete [] current;
	}
	
	//Tcl_IncrRefCount(paramList);
		
	return 0;
}
