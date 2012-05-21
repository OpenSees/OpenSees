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
**   Quan Gu (qgu@ucsd.edu)                                           **
**   Joel P. Conte (jpconte@ucsd.edu)                                 **
** ****************************************************************** */
                                                                        
 
//
// Written by  Quan Gu UCSD
//
 

#if !defined(AFX_DP_RSM_SIM_TimeVariant_H)
#define AFX_DP_RSM_SIM_TimeVariant_H

#include <ReliabilityAnalysis.h>
#include <GridPlane.h>
//#include <Hessian.h>
#include <GradientEvaluator.h>
#include <UnivariateDecomposition.h>
#include <UniformExperimentalPointRule1D.h>
#include <BivariateDecomposition.h>
#include <RespSurfaceSimulation.h>


#include <ReliabilityAnalysis.h>

class DP_RSM_Sim_TimeVariant : public ReliabilityAnalysis  
{

public:

	DP_RSM_Sim_TimeVariant(ReliabilityDomain *passedReliabilityDomain,
					FunctionEvaluator *passedGFunEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
					char *passedOutputFileName,
					GradientEvaluator * passedGradGEvaluator, Vector * pDesignPt, int numAxis, 
					char * typeExpPtRule,Tcl_Interp *passedTclInterp, 
					Matrix * passedHessian, char * passedHessianFile, char * typeSurfaceDesign, 
					char * typeRespSurfaceSimulation, Vector * gridInfo,
					RandomNumberGenerator * pRandomNumberGenerator,
					double pTargetCOV,
					int pNumberOfSimulations, double pLittleDt,
					double ImpulseInterval);
	virtual ~DP_RSM_Sim_TimeVariant();


	int getNumOfGridPlane( int i, int j);
	int getNumOfAxis(int Axis, int numGridPlane);
	void setGridInfo( Vector * GridData );
	int analyze();
	double FEDivergenceCorrection(int numGridPlane, int startPtX, int startPtY, int incr_i_x, int incr_i_y, int m, int n, int type);
	double FEDivergenceCorrection2(int numGridPlane, int startPtX, int startPtY, int incr_i_x, int incr_i_y, int m, int n, int type);

private:
	double ImpulseInterval;
	double littleDt;
	ExperimentalPointRule1D * theExpPtRule;

	Tcl_Interp * theInterp;
	GradientEvaluator * theGradGEvaluator;
	FunctionEvaluator * theGFunEvaluator;
	ProbabilityTransformation * theProbabilityTransformation;
	ReliabilityDomain * theReliabilityDomain;

	PrincipalAxis ** thePrincipalAxes;
	GridPlane ** theGridPlanes;


    SurfaceDesign*  theSurfaceDesign;
//	ExperimentalPointRule*  theExperimentalPointRule;
	RespSurfaceSimulation*  theRespSurfaceSimulation;
    


	char  * HessianFileName;

	
	Matrix * HessianMatrix;

	Vector * theDesignPtXSpace;
	Vector * theDesignPoint; //U space
	char outputFileName[20];
	Matrix * rotation;
	//Hessian * theHessian;
	int numOfPrincipalAxes;

	RandomNumberGenerator * theRandomNumberGenerator;
    double theTargetCOV;
	int theNumberOfSimulations;
 

};

#endif //
