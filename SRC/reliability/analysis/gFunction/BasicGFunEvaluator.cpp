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
                                                                        
// $Revision: 1.9 $
// $Date: 2010-09-13 21:39:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/BasicGFunEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <BasicGFunEvaluator.h>
#include <Vector.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>

#include <tcl.h>
#include <string.h>


BasicGFunEvaluator::BasicGFunEvaluator(Tcl_Interp *passedTclInterp, 
				       ReliabilityDomain *passedReliabilityDomain)

  :GFunEvaluator(), theTclInterp(passedTclInterp), theReliabilityDomain(passedReliabilityDomain)
{
}

BasicGFunEvaluator::~BasicGFunEvaluator()
{
}

double
BasicGFunEvaluator::evaluateGMHS(const Vector &x) 
{
  return 0.0;
}

int
BasicGFunEvaluator::setNamespaceRandomVariables(const Vector &x)
{
  return 0;
}

int
BasicGFunEvaluator::tokenizeSpecials(TCL_Char *theExpression, Tcl_Obj *paramList)
{
	// No need to set additional quantities here
	// (the basic random variables are set in the base class)

	return 0;
}
