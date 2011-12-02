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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-04-28 20:51:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/LimitStateFunction.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef LimitStateFunction_h
#define LimitStateFunction_h

#include <ReliabilityDomainComponent.h>
#include <Vector.h>
//#include <fstream>

class LimitStateFunction : public ReliabilityDomainComponent
{

public:
	LimitStateFunction(	int tag, 
						const char *expression);
	~LimitStateFunction();
	void Print(OPS_Stream &s, int flag =0);

	// Method to get/add limit-state function
	char *getExpression();
	char *getTokenizedExpression();
	int addExpression(char *expression);
	int removeAddedExpression();

	// FORM analysis:
	double GFunValueAtStartPt;
	double GFunValueAtEndPt;
	double FORMReliabilityIndexBeta;
	double FORMProbabilityOfFailure_pf1;
	Vector designPoint_x_inOriginalSpace;
	Vector designPoint_u_inStdNormalSpace;
	Vector normalizedNegativeGradientVectorAlpha;
	Vector importanceVectorGamma;
	int numberOfStepsToFindDesignPointAlgorithm;
	
	// From Simulation analysis:
	double SimulationReliabilityIndexBeta;
	double SimulationProbabilityOfFailure_pfsim;
	double CoefficientOfVariationOfPfFromSimulation;
	int NumberOfSimulations;
	
	// From SORM analysis:
	double SORMCurvatureFittingBetaBreitung;
	double SORMCurvatureFittingPf2Breitung;
	double SORMPointFittingBetaBreitung;
	double SORMPointFittingPf2Breitung;
	double SORMUsingSearchBetaBreitung;
	double SORMUsingSearchPf2Breitung;
	Vector lastSearchDirection;
	int numberOfCurvatauresUsed;
	Vector secondLast_u;
	Vector secondLastAlpha;

protected:

private:

	int tokenizeIt(const char *expression);

	int tag;
	char *originalExpression;
	char *tokenizedExpression;
	char *expressionWithAddition;

};

#endif
