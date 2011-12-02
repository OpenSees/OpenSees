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
                                                                        
// $Revision: 1.15 $
// $Date: 2010-06-10 20:14:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/LimitStateFunction.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef LimitStateFunction_h
#define LimitStateFunction_h

#include <ReliabilityDomainComponent.h>
#include <Vector.h>
#include <tcl.h>

#include <map>
#include <string>
using namespace std;

class LimitStateFunction : public ReliabilityDomainComponent
{

public:
	//LimitStateFunction(int tag, int classTag);
	LimitStateFunction(int tag, TCL_Char *expression, Tcl_Interp *theTclInterp);
	~LimitStateFunction();

	void Print(OPS_Stream &s, int flag =0);

	// Method to get/add limit-state function
	char *getExpression();
	char *getTokenizedExpression();
	Tcl_Obj *getParameters();
	int addExpression(char *expression);
	int removeAddedExpression();

	// Methods to add/remove/get gradient of LSF
	int addGradientExpression(const char *expression, int rvTag);
	int removeGradientExpression(int rvTag);
	const char* getGradientExpression(int rvTag);

	// THE REMAINING PUBLIC METHODS SHOULD BE IN OTHER CLASSES - MHS

	// From FORM analysis
	void setFORM_startGFun(double value) {GFunValueAtStartPt = value;}
	double getFORM_startGFun(void) {return GFunValueAtStartPt;}

	void setFORM_endGFun(double value) {GFunValueAtEndPt = value;}
	double getFORM_endGFun(void) {return GFunValueAtEndPt;}

	void setFORM_beta(double value) {FORMReliabilityIndexBeta = value;}
	double getFORM_beta(void){return FORMReliabilityIndexBeta;}

	void setFORM_pf1(double value) { FORMProbabilityOfFailure_pf1 = value;}
	double getFORM_pf1(void){return  FORMProbabilityOfFailure_pf1;}

	void setFORM_x(const Vector &value) {designPointX = value;}
	const Vector &getFORM_x(void) {return designPointX;}
	void setFORM_u(const Vector &value) {designPointU = value;}
	const Vector &getFORM_u(void) {return designPointU;}

	void setFORM_alpha(const Vector &value) {normalizedNegativeGradientVectorAlpha = value;}
	const Vector &getFORM_alpha(void) {return normalizedNegativeGradientVectorAlpha;}

	void setFORM_gamma(const Vector &value) {importanceVectorGamma = value;}
	const Vector &getFORM_gamma(void) {return importanceVectorGamma;}

	void setFORM_numSteps(int value) {numberOfStepsToFindDesignPointAlgorithm = value;}
	int getFORM_numSteps(void) {return numberOfStepsToFindDesignPointAlgorithm;}

	// From Simulation analysis:
	void setSIM_beta(double value) {SimulationReliabilityIndexBeta = value;}
	double getSIM_beta(void) {return SimulationReliabilityIndexBeta;}

	void setSIM_pfsim(double value) {SimulationProbabilityOfFailure_pfsim = value;}
	double getSIM_pfsim(void) {return SimulationProbabilityOfFailure_pfsim;}
	void setSIM_pfcov(double value) {CoefficientOfVariationOfPfFromSimulation = value;}
	double getSIM_pfcov(void) {return CoefficientOfVariationOfPfFromSimulation;}
	void setSIM_numsim(int value) {NumberOfSimulations = value;}
	int getSIM_numsim(void) {return NumberOfSimulations;}
	
	// SORM analysis
	void setSORM_cf_beta_breitung(double value) {SORMCurvatureFittingBetaBreitung = value;}
	double getSORM_cf_beta_breitung(void) {return SORMCurvatureFittingBetaBreitung;}

	void setSORM_cf_pf2_breitung(double value) {SORMCurvatureFittingPf2Breitung = value;}
	double getSORM_cf_pf2_breitung(void) {return SORMCurvatureFittingPf2Breitung;}

	void setSORM_pf_beta_breitung(double value) {SORMPointFittingBetaBreitung = value;}
	double getSORM_pf_beta_breitung(void) {return SORMPointFittingBetaBreitung;}

	void setSORM_pf_pf2_breitung(double value) {SORMPointFittingPf2Breitung = value;}
	double getSORM_pf_pf2_breitung(void) {return SORMPointFittingPf2Breitung;}

	void setSORM_us_beta_breitung(double value) {SORMUsingSearchBetaBreitung = value;}
	double getSORM_us_beta_breitung(void) {return SORMUsingSearchBetaBreitung;}

	void setSORM_us_pf2_breitung(double value) {SORMUsingSearchPf2Breitung = value;}
	double getSORM_us_pf2_breitung(void) {return SORMUsingSearchPf2Breitung;}

	void setSORM_numCurvatures(int value) {numberOfCurvaturesUsed = value;}
	int getSORM_numCurvatures(void) {return numberOfCurvaturesUsed;}

	// General stuff
	void setLastSearchDirection(const Vector &value) {lastSearchDirection = value;}
	const Vector &getLastSearchDirection(void) {return lastSearchDirection;}
	void setSecondLast_u(const Vector &value) {secondLast_u = value;}
	const Vector &getSecondLast_u(void) {return secondLast_u;}
	
	void setSecondLast_alpha(const Vector &value) {secondLastAlpha = value;}
	const Vector &getSecondLast_alpha(void) {return secondLastAlpha;}

protected:

private:
	// FORM analysis:
	double GFunValueAtStartPt;
	double GFunValueAtEndPt;
	double FORMReliabilityIndexBeta;
	double FORMProbabilityOfFailure_pf1;
	Vector designPointX;
	Vector designPointU;
	Vector normalizedNegativeGradientVectorAlpha;
	Vector importanceVectorGamma;
	int numberOfStepsToFindDesignPointAlgorithm;

	// Simulation analysis
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
	int numberOfCurvaturesUsed;

	// General stuff
	Vector lastSearchDirection;
	Vector secondLast_u;
	Vector secondLastAlpha;

	void initializeFORMAnalysis(void);
	void initializeSimulationAnalysis(void);
	void initializeSORMAnalysis(void);

	int tokenizeIt(const char *expression);

	char *originalExpression;
	char *expressionWithAddition;
	Tcl_Interp *theTclInterp;
	Tcl_Obj *paramList;

	// STL map for analytic gradients of this LSF
	map<int, string> mapOfGradientExpressions;

};

#endif
