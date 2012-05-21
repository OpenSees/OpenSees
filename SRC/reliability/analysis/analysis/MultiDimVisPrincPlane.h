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
 

#if !defined AFX_MULTIDIMVISPRINCPLANE_H
#define AFX_MULTIDIMVISPRINCPLANE_H


#include <GFunVisualizationAnalysis.h>

// command: MultiDimVisualPrinPlane -funSurf function -designPt dp.out  -output vis.out  -ndir $n <gridInfo {0  minY  maxY nPts0 1 minX1  maxX1 nPts1 2 minX2  maxX2 nPts2 ...}>
 //command: MultiDimVisualPrinPlane -funSurf surface -designPt dp.out  -output vis.out -ndir $n <gridInfo {} -timevariant -littleDt dt> //not yet implemented
#include <PrincipalPlane.h>
//#include <Hessian.h>
//#include <GradGEvaluator.h>


class MultiDimVisPrincPlane : public ReliabilityAnalysis  
{
public:
	int setGridInfo( Vector * gridInfo, int pNumOfPPlane);
	int analyze();
	MultiDimVisPrincPlane(ReliabilityDomain *passedReliabilityDomain,
					FunctionEvaluator *passedGFunEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
					char *passedOutputFileName,
					GradientEvaluator * passedGradGEvaluator,Vector * DesignPtr,int pNumOfPPlane, int type, 
					Vector *pVector,Tcl_Interp *passedTclInterp, 
					Matrix * passedHessian, char * passedHessianFile, int pAnalysisType, double pLittleDt);
	virtual ~MultiDimVisPrincPlane();

private:
	Vector * valuesG2OfAxis;
	double littleDt;
	int analysisType;
	Vector * valuesOfAxis;
	char  * HessianFileName;
	Tcl_Interp * theInterp;
	Vector * theDesignPtXSpace;
	GradientEvaluator * theGradGEvaluator;
	Matrix * HessianMatrix;
	Vector * theDesignPoint; //U space
	FunctionEvaluator * theGFunEvaluator;
	ProbabilityTransformation * theProbabilityTransformation;
	char outputFileName[20];
	ReliabilityDomain * theReliabilityDomain;
	PrincipalPlane ** thePrincipalPlanes;
	Matrix * rotation;
	//Hessian * theHessian;
//	Vector * theDesignPt;
	int numOfPrinPlane;
	int type;
};

#endif // !defined AFX_MULTIDIMVISPRINCPLANE_H
