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
 
//

#ifndef HESSIAN_H
#define HESSIAN_H


#include <ProbabilityTransformation.h>
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <HessianApproximation.h>
#include <ReliabilityDomain.h>
#include <Matrix.h>
#include <Vector.h>

class Hessian : public HessianApproximation  
{

public:
  Hessian(int size,
	  ReliabilityDomain *passedReliabilityDomain,
	  ProbabilityTransformation *,
	  FunctionEvaluator *passedGFunEvaluator,
	  GradientEvaluator *passedGradGEvaluator,
	  double tol);
  
  ~Hessian();

  int setHessianToIdentity(int size);
  int updateHessianApproximation(const Vector &u_old,
				 double g_old,
				 const Vector &gradG_old,
				 double stepSize,
				 const Vector &searchDirection,
				 double g_new,
				 const Vector &gradG_new);

  const Matrix  &getHessianApproximation();  


  double getNormOfGradientU();
  int refineHessian(int time,int colOfhessian);
  int diagnose(int i);
  int formReducedHessian( Vector * Designpoint_X, Matrix * Hessian);
  int formReducedHessian( Vector * pDesignPt_X);

  const Matrix &getReducedHessian();  
  const Matrix &getRotationMatrix(Vector alpha);
  const Matrix &getHessianInPhysicalSpace();

  int formHessianBySNOPT();
  int formHessianByFDM(int numOfLimitStateFunction,Vector * theDesignPoint);

  const Vector &getEigenvector(int mode);
  double getEigenvalue(int mode);
  int setSizeOfEigenVector(int newSize);
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  
 private:
  Vector * standSens;
  double perturbTol;
  Vector * designPointUSpace;
  double normOfGradientInStandardNormalSpace;
  Matrix * theHessianInPhysicalSpace;
  ReliabilityDomain * theReliabilityDomain;
  int sizeOfHessian;
  Matrix * theHessian;   // in U space
  Matrix * theReducedHessian;   // in U space
  ProbabilityTransformation * theProbabilityTransformation;
  FunctionEvaluator * theGFunEvaluator;
  GradientEvaluator * theGradGEvaluator;
};

#endif // !defined HESSIAN_H
