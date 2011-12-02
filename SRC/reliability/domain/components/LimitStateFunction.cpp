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
// $Date: 2007-02-08 01:25:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/LimitStateFunction.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <LimitStateFunction.h>
#include <Vector.h>
#include <string.h>
#include <classTags.h>

LimitStateFunction::LimitStateFunction(	int passedTag, 
									    TCL_Char *passedExpression)
:ReliabilityDomainComponent(passedTag, LIMIT_STATE_FUNCTION)
{
	strcpy(originalExpression,passedExpression);

	strcpy(expressionWithAddition,passedExpression);

	tokenizeIt(passedExpression);

	this->initializeFORMAnalysis();
	this->initializeSimulationAnalysis();
	this->initializeSORMAnalysis();
}

void
LimitStateFunction::initializeFORMAnalysis(void)
{
  GFunValueAtStartPt = 0.0;
  GFunValueAtEndPt = 0.0;
  FORMReliabilityIndexBeta = 0.0;
  FORMProbabilityOfFailure_pf1 = 0.0;
  //Vector designPoint_x_inOriginalSpace;
  //Vector designPoint_u_inStdNormalSpace;
  //Vector normalizedNegativeGradientVectorAlpha;
  //Vector importanceVectorGamma;
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
  //Vector lastSearchDirection;
  numberOfCurvatauresUsed = 0;
  //Vector secondLast_u;
  //Vector secondLastAlpha;
}

LimitStateFunction::~LimitStateFunction()
{
  
}


void
LimitStateFunction::Print(OPS_Stream &s, int flag)  
{
}



char *
LimitStateFunction::getExpression()
{
	return expressionWithAddition;
}


char *
LimitStateFunction::getTokenizedExpression()
{
	return tokenizedExpression;
}

int
LimitStateFunction::addExpression(char *addition)
{
	strcat(expressionWithAddition,addition);

	tokenizeIt(expressionWithAddition);

	return 0;
}

int
LimitStateFunction::removeAddedExpression()
{
	strcpy(expressionWithAddition,originalExpression);

	tokenizeIt(expressionWithAddition);

	return 0;
}


int
LimitStateFunction::tokenizeIt(TCL_Char *originalExpression)
{
	// Also store the tokenized expression (with dollar signs in front of variable names)
	char lsf_forTokenizing[500];
	char separators[5] = "}{";
	char *dollarSign = "$";
	strcpy(lsf_forTokenizing,originalExpression);
	char lsf_expression[500] = "";
	char *tokenPtr2 = strtok( lsf_forTokenizing, separators);
	while ( tokenPtr2 != NULL ) {
		if (   strncmp(tokenPtr2, "a",1) == 0
			|| strncmp(tokenPtr2, "b",1) == 0
			|| strncmp(tokenPtr2, "c",1) == 0
			|| strncmp(tokenPtr2, "d",1) == 0
			|| strncmp(tokenPtr2, "e",1) == 0
			|| strncmp(tokenPtr2, "f",1) == 0
			|| strncmp(tokenPtr2, "g",1) == 0
			|| strncmp(tokenPtr2, "h",1) == 0
			|| strncmp(tokenPtr2, "i",1) == 0
			|| strncmp(tokenPtr2, "j",1) == 0
			|| strncmp(tokenPtr2, "k",1) == 0
			|| strncmp(tokenPtr2, "l",1) == 0
			|| strncmp(tokenPtr2, "m",1) == 0
			|| strncmp(tokenPtr2, "n",1) == 0
			|| strncmp(tokenPtr2, "o",1) == 0
			|| strncmp(tokenPtr2, "p",1) == 0
			|| strncmp(tokenPtr2, "q",1) == 0
			|| strncmp(tokenPtr2, "r",1) == 0
			|| strncmp(tokenPtr2, "s",1) == 0
			|| strncmp(tokenPtr2, "t",1) == 0
			|| strncmp(tokenPtr2, "u",1) == 0
			|| strncmp(tokenPtr2, "v",1) == 0
			|| strncmp(tokenPtr2, "w",1) == 0
			|| strncmp(tokenPtr2, "x",1) == 0
			|| strncmp(tokenPtr2, "y",1) == 0
			|| strncmp(tokenPtr2, "z",1) == 0) {
			strcat(lsf_expression, dollarSign);
			strcat(lsf_expression, tokenPtr2);
		}
		else {
			strcat(lsf_expression, tokenPtr2);
		}
		tokenPtr2 = strtok( NULL, separators);
	}

	strcpy(tokenizedExpression,lsf_expression);

	return 0;
}
