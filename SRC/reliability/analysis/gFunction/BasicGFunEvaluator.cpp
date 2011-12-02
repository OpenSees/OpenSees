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
// $Date: 2003-10-27 23:45:43 $
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

:GFunEvaluator(passedTclInterp, passedReliabilityDomain)
{
}

BasicGFunEvaluator::~BasicGFunEvaluator()
{
}


int
BasicGFunEvaluator::runGFunAnalysis(Vector x)
{
	// Nothing to compute for this kind of gFunEvaluator
	
	return 0;
}


int
BasicGFunEvaluator::tokenizeSpecials(TCL_Char *theExpression)
{
	// No need to set additional quantities here
	// (the basic random variables are set in the base class)

	return 0;
}
