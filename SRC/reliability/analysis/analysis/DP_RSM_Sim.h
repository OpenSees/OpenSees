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

#ifndef DP_RSM_SIM_
#define DP_RSM_SIM_

#include <ReliabilityAnalysis.h>
#include <GridPlane.h>
//#include <Hessian.h>
#include <GradientEvaluator.h>
#include <UnivariateDecomposition.h>
#include <BivariateDecomposition.h>

#include <UniformExperimentalPointRule1D.h>
#include <RespSurfaceSimulation.h>


// command: runDP_RSM_SimAnalysis -designPt dp.out  -output results.out  -ndir $n <-experimentalPointRule Uniform -gridInfo {0  minY  maxY nPts0 1 minX1  maxX1 nPts1 2 minX2  maxX2 nPts2 ...}>



class DP_RSM_Sim : public ReliabilityAnalysis 
{
 public:
  DP_RSM_Sim(ReliabilityDomain *passedReliabilityDomain,
	     FunctionEvaluator *passedGFunEvaluator,
	     ProbabilityTransformation *passedProbabilityTransformation,
	     char *passedOutputFileName,
	     GradientEvaluator * passedGradGEvaluator, Vector * pDesignPt, int numAxis, 
	     char * typeExpPtRule,Tcl_Interp *passedTclInterp, 
	     Matrix * passedHessian, char * passedHessianFile, char * typeSurfaceDesign, 
	     char * typeRespSurfaceSimulation, Vector * gridInfo,
	     RandomNumberGenerator * pRandomNumberGenerator,
	     double pTargetCOV,
	     int pNumberOfSimulations);
  ~DP_RSM_Sim();
  
  double FEDivergenceCorrection(int numGridPlane,int startPtX, int startPtY, int incr_i_x, int incr_i_y, int m, int n, int type);
  int getNumOfGridPlane( int i, int j);
  int getNumOfAxis(int Axis, int numGridPlane);
  void setGridInfo( Vector * GridData );
  int analyze();
  
  
 private:
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


#endif
