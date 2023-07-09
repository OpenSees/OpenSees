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
// $Date: 2008-05-27 20:04:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/GradientEvaluator.h,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#ifndef GradientEvaluator_h
#define GradientEvaluator_h

#include <Vector.h>
#include <ReliabilityDomain.h>
#include <FunctionEvaluator.h>

class PerformanceFunctionCoefficientIter;

class GradientEvaluator
{
 public:
  GradientEvaluator(ReliabilityDomain *theReliabilityDomain, 
		    FunctionEvaluator *theGFunEvaluator);
  virtual ~GradientEvaluator();
  
  // Methods provided by the sub-classes
  virtual int computeGradient(double gFunValue) = 0;
  virtual const Vector &getGradient(void) = 0;
  
  //////S added by K Fujimura ///////////
  int initializeNumberOfEvaluations();
  int getNumberOfEvaluations();
  virtual void setPerformFuncCoeffs(TaggedObjectStorage*);
  virtual void setPerformFuncCoeffIter(PerformanceFunctionCoefficientIter*);
  virtual bool getfinitedifference(){ return finitedifference; }
  virtual void setfinitedifference(bool fdf){finitedifference=fdf;}
  //////E added by K Fujimura ///////////
  
 protected:
  int computeParameterDerivatives(double g);
  
  // one day these should find themselves a better home than protected
  ReliabilityDomain *theReliabilityDomain;
  FunctionEvaluator *theFunctionEvaluator;
  
  /////S added by K Fujimura /////
  int numberOfEvalIncSens;
  /////E added by K Fujimura /////
  
 private:
  
  /////S added by K Fujimura /////
  bool finitedifference;
  /////E added by K Fujimura ////
};

#endif
