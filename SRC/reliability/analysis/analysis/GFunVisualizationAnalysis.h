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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-04-28 20:51:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/GFunVisualizationAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef GFunVisualizationAnalysis_h
#define GFunVisualizationAnalysis_h

#include <ReliabilityAnalysis.h>
#include <GFunEvaluator.h>
#include <FindDesignPointAlgorithm.h>
#include <ProbabilityTransformation.h>
#include <ReliabilityDomain.h>
#include <GradGEvaluator.h>
#include <MeritFunctionCheck.h>
#include <ReliabilityConvergenceCheck.h>
#include <RootFinding.h>


class GFunVisualizationAnalysis : public ReliabilityAnalysis
{

public:
	GFunVisualizationAnalysis(
					ReliabilityDomain *theReliabilityDomain,
					GFunEvaluator *theGFunEvaluator,
					ProbabilityTransformation *theProbabilityTransformation,
					const char *outputFileName,
					const char *convFileName,
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
	int setGradGEvaluator(GradGEvaluator *theGradGEvaluator);
	int setMeritFunctionCheck(MeritFunctionCheck *theMeritFunctionCheck);
	int setReliabilityConvergenceCheck(ReliabilityConvergenceCheck *theReliabilityConvergenceCheck);


protected:

private:

	Vector getCurrentAxes12Point(int i, int j);
	Vector getCurrentAxes3Point(int i, int j);
	double findGSurface(Vector thePoint);
	double evaluateGFunction(Vector thePoint, bool printDivision);





	ReliabilityDomain *theReliabilityDomain;
	GFunEvaluator *theGFunEvaluator;
	ProbabilityTransformation *theProbabilityTransformation;
	MeritFunctionCheck *theMeritFunctionCheck;
	GradGEvaluator *theGradGEvaluator;
	ReliabilityConvergenceCheck *theReliabilityConvergenceCheck;
	Vector *theStartPoint;
	RootFinding *theRootFindingAlgorithm;

	char *outputFileName;
	char *convFileName;
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
};

#endif
