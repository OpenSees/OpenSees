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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/TclGFunEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef TclGFunEvaluator_h
#define TclGFunEvaluator_h

#include <GFunEvaluator.h>
#include <Vector.h>
#include <Domain.h>
#include <ReliabilityDomain.h>
#include <tcl.h>


class TclGFunEvaluator : public GFunEvaluator
{
 public:
  TclGFunEvaluator(Tcl_Interp *passedTclInterp,
		   ReliabilityDomain *passedReliabilityDomain,
		   Domain *passedOpenSeesDomain,
		   TCL_Char *fileName);
  ~TclGFunEvaluator();
  
  int evaluateG(const Vector &x);
  double getG(void);
  int runGFunAnalysis(const Vector &x);

  int tokenizeSpecials(TCL_Char *theExpression, Tcl_Obj *passedList);
  
 protected:
  
 private:
  
  int setTclRandomVariables(const Vector &x);
  
  Tcl_Interp *theTclInterp;
  ReliabilityDomain *theReliabilityDomain;
  Domain *theOpenSeesDomain;
  
  char fileName[256];
  
  double g;
};

#endif
