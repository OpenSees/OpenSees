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


// $Revision: 1.4 $
// $Date: 2008-10-22 16:41:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewSearchWithStepSizeAndStepDirection.h,v $



#ifndef NewSearchWithStepSizeAndStepDirection_h
#define NewSearchWithStepSizeAndStepDirection_h

#include <FindDesignPointAlgorithm.h>
#include <StepSizeRule.h>
#include <SearchDirection.h>
#include <ProbabilityTransformation.h>
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <HessianEvaluator.h>
#include <ReliabilityConvergenceCheck.h>
#include <Matrix.h>
#include <Vector.h>
#include <ReliabilityDomain.h>

#include <fstream>
using std::ofstream;

class NewSearchWithStepSizeAndStepDirection : public FindDesignPointAlgorithm
{

public:

	// Constructor and destructor
	NewSearchWithStepSizeAndStepDirection(
					int passedMaxNumberOfIterations, 
					ReliabilityDomain *passedReliabilityDomain,
					FunctionEvaluator *passedGFunEvaluator,
					GradientEvaluator *passedGradGEvaluator,
					StepSizeRule *passedStepSizeRule,
					SearchDirection *passedSearchDirection,
					ProbabilityTransformation *passedProbabilityTransformation,
					HessianEvaluator *theHessianEvaluator,
					ReliabilityConvergenceCheck *theReliabilityConvergenceCheck,
					bool startAtOrigin,
					int printFlag,
					char *fileNamePrint);
	~NewSearchWithStepSizeAndStepDirection();
	
	int findDesignPoint();

	const Vector &get_x();
	const Vector &get_u();
	const Vector &get_alpha();
	const Vector &get_gamma();
	int getNumberOfSteps();
	const Vector &getSecondLast_u();
	const Vector &getSecondLast_alpha();
	const Vector &getLastSearchDirection();
	double getFirstGFunValue();
	double getLastGFunValue();
	const Vector &getGradientInStandardNormalSpace();
	int    getNumberOfEvaluations();

	//  Modified by K Fujimura 10/10/2004
	int    getNumberOfSensAna();
	double get_check1_init();
	double get_check2_init();
	double get_check1_conv();
	double get_check2_conv();
	
	void set_x(Vector&);
	void set_u(Vector&);
	virtual double get_beta();
	Matrix getJacobian_x_u(){ return (*jacobian_x_u);}

	int setStartPt(Vector *s);	
	//  Modified by K Fujimura 10/10/2004

protected:

private:	

	// The reliability domain and tools for the analysis
	ReliabilityDomain *theReliabilityDomain;
	FunctionEvaluator *theGFunEvaluator;
	GradientEvaluator *theGradGEvaluator;
	StepSizeRule *theStepSizeRule;
	SearchDirection *theSearchDirection;
	ProbabilityTransformation *theProbabilityTransformation;
	HessianEvaluator *theHessianEvaluator;
	ReliabilityConvergenceCheck *theReliabilityConvergenceCheck;

	bool startAtOrigin;

	// Private member functions to do the job
	int doTheActualSearch(bool doRvProjection);
	int doRvProjection(Vector uOld, Vector uNew);

	// Data members set when the object is created
	int maxNumberOfIterations;
	int resFlag;

	// Data members where the results are to be stored
	Vector* xinit;
	Vector* x;
	Vector* u;
	Vector* u_old;
	Vector* uSecondLast;
	Vector* uThirdLast;
	Vector* alpha;
	Vector* alpha_old;
	Vector* alphaSecondLast;
	Vector* alphaThirdLast;
	Vector* gradientInStandardNormalSpace;
	Vector* gradientInStandardNormalSpace_old;
	Vector* gamma;
	Vector* searchDirection;
	Vector* gradientOfgFunction;
	Matrix* jacobian_x_u;

	double ampfunc;
	double chkfunc;

	// at minimimu point //
	Vector* x_min;
	Vector* u_min;
	Vector* alpha_min;
	Vector* gradientInStandardNormalSpace_min;
	Vector* gradientOfgFunction_min;
	Matrix* jacobian_x_u_min;
	double gFunctionValue_min;
	double normOfGradient_min;
	double chk1_min;
	double chk2_min;
	double chk_min;

	int iter;
//	int i;
	double Gfirst;
	double Glast;
	double gFunctionValue;
	double gFunctionValue_old;
	double normOfGradient;
	double stepSize;
	double chk1;
	double chk2;
	int reschk;
	double ampchk;


	// Data members set through the call when a job is to be done
	//Vector *startPoint;
	Vector *designPoint_uStar;

	int printFlag;
	char fileNamePrint[256];
	int numberOfEvaluations;
	int numberOfSensAna;

	double check1_init;
	double check2_init;
	double check1_conv;
	double check2_conv;
	bool fdf;
	bool convergenceAchieved;

	ofstream output;
	ofstream outputTMP;

  	void errorMmessage_setx();
	void errorMessage_xtou();
	void errorMessage_setu();
	void errorMessage_utox();
	void errorMessage_gfun();
	void errorMessage_evalg();
	void errorMessage_compGradg();
	void errorMessage_checkGradg();
	void errorMessage_zeroGradg();
	void delxinit();
	int lineSearch();
};

#endif
