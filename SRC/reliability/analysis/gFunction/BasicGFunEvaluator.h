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
                                                                        
// $Revision: 1.7 $
// $Date: 2010-09-13 21:39:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/BasicGFunEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef BasicGFunEvaluator_h
#define BasicGFunEvaluator_h

#include <GFunEvaluator.h>
#include <Vector.h>
#include <ReliabilityDomain.h>
#include <tcl.h>


class BasicGFunEvaluator : public GFunEvaluator
{

 public:
  BasicGFunEvaluator(Tcl_Interp *passedTclInterp, 
		     ReliabilityDomain *passedReliabilityDomain);
  ~BasicGFunEvaluator();

  double evaluateGMHS(const Vector &x);
  int setNamespaceRandomVariables(const Vector &x);
  int runGFunAnalysis(const Vector &x) {return 0;}

  int tokenizeSpecials(TCL_Char *theExpression, Tcl_Obj *paramList);

 protected:

 private:
  Tcl_Interp *theTclInterp;
  ReliabilityDomain *theReliabilityDomain;
};

#endif
