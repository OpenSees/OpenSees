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
// $Date: 2008-10-22 16:41:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/GFunVisualizationAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef GFunVisualizationAnalysis_h
#define GFunVisualizationAnalysis_h

#include <fstream>
using std::ifstream;

#include <ReliabilityAnalysis.h>
#include <FunctionEvaluator.h>
#include <FindDesignPointAlgorithm.h>
#include <ProbabilityTransformation.h>
#include <ReliabilityDomain.h>
#include <GradientEvaluator.h>
#include <MeritFunctionCheck.h>
#include <ReliabilityConvergenceCheck.h>
#include <RootFinding.h>


class GFunVisualizationAnalysis : public ReliabilityAnalysis
{

public:
	GFunVisualizationAnalysis(
					ReliabilityDomain *theReliabilityDomain,
					FunctionEvaluator *theGFunEvaluator,
					ProbabilityTransformation *theProbabilityTransformation,
					TCL_Char *outputFileName,
					TCL_Char *convFileName,
					int convResults,
					int space,
					int funSurf,
					int axes,
					int dir);
	virtual ~GFunVisualizationAnalysis();

	int analyze(void);

	int setDirection(int rvDir);
	int setDirection(Vector theDirectionVector);
	
	int setAxes(Vector axesVector);
	int setAxes(Matrix theMatrix);
	int setNumLinePts(int numLinePts);
	
	int setRootFindingAlgorithm(RootFinding *theRootFinder);
	int setStartPoint(Vector *theStartPoint);
	int setGradGEvaluator(GradientEvaluator *theGradGEvaluator);
	int setMeritFunctionCheck(MeritFunctionCheck *theMeritFunctionCheck);
	int setReliabilityConvergenceCheck(ReliabilityConvergenceCheck *theReliabilityConvergenceCheck);


protected:

private:

	Vector getCurrentAxes12Point(int i, int j);
	Vector getCurrentAxes3Point(int i, int j);
	double findGSurface(Vector thePoint);
	double evaluateGFunction(Vector thePoint, bool printDivision);

	ReliabilityDomain *theReliabilityDomain;
	FunctionEvaluator *theGFunEvaluator;
	ProbabilityTransformation *theProbabilityTransformation;
	MeritFunctionCheck *theMeritFunctionCheck;
	GradientEvaluator *theGradGEvaluator;
	ReliabilityConvergenceCheck *theReliabilityConvergenceCheck;
	RootFinding *theRootFindingAlgorithm;

	char outputFileName[256];
	char convFileName[256];
	int convResults;
	int space;
	int funSurf;
	int axes;
	int dir;

	int rvDir;
	Vector theDirectionVector;
	int rv1, rv2;
	double from1, interval1;
	double from2, interval2;
	int numPts1, numPts2;
	Matrix theMatrix;
	int numLinePts;

	int nrv;
	double scaleValue;
	
	ofstream convFile;
};

#endif
