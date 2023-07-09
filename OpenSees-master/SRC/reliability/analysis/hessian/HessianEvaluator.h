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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/HessianEvaluator.h,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#ifndef HessianEvaluator_h
#define HessianEvaluator_h

#include <Matrix.h>
#include <ReliabilityDomain.h>
#include <FunctionEvaluator.h>


class HessianEvaluator
{
 public:
  HessianEvaluator(ReliabilityDomain *theReliabilityDomain, 
		    FunctionEvaluator *theGFunEvaluator);
  virtual ~HessianEvaluator();
  
  // Methods provided by the sub-classes
  virtual int computeHessian(void) = 0;
  virtual Matrix getHessian(void) = 0;
    
 protected:
  // one day these should find themselves a better home than protected
  ReliabilityDomain *theReliabilityDomain;
  FunctionEvaluator *theFunctionEvaluator;
  
 private:
    
};

#endif
