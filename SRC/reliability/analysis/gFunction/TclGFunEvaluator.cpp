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
// $Date: 2003-04-28 20:51:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/TclGFunEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <TclGFunEvaluator.h>
#include <Vector.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>
//#include <fstream>
#include <tcl.h>
#include <string.h>


TclGFunEvaluator::TclGFunEvaluator(Tcl_Interp *passedTclInterp,
					ReliabilityDomain *passedReliabilityDomain,
					const char *passed_fileName)
:GFunEvaluator(passedTclInterp, passedReliabilityDomain)
{
	fileName = new char[256];
	strcpy(fileName,passed_fileName);
}


TclGFunEvaluator::~TclGFunEvaluator()
{
	if (fileName != 0)
		delete [] fileName;
}


int
TclGFunEvaluator::runGFunAnalysis(Vector x)
{
	// Initial declarations
	char theCommand[100];
	int i;


	// Set values of random variables in the Tcl intepreter
	for (i=0; i<x.Size(); i++) {
		sprintf(theCommand,"set x_%d %20.12e",(i+1),x(i));
		Tcl_Eval( theTclInterp, theCommand );
	}


	// Source the code file that the user has provided
	sprintf(theCommand,"source %s",fileName);
	Tcl_Eval( theTclInterp, theCommand );

	return 0;
}




int
TclGFunEvaluator::tokenizeSpecials(char *theExpression)
{
	// No specific new quantities in the performance 
	// function is introduced here. 

	return 0;
}

