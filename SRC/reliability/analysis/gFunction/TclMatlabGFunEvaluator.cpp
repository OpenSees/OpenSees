#include <TclMatlabGFunEvaluator.h>
#include <Vector.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>
#include <ReliabilityDomain.h>

#include <tcl.h>
#include <string.h>

#include <fstream.h>
#include <engine.h>
#include <mex.h>



TclMatlabGFunEvaluator::TclMatlabGFunEvaluator(Tcl_Interp *passedTclInterp,
											   ReliabilityDomain *passedReliabilityDomain)

:GFunEvaluator()
{
	theTclInterp			= passedTclInterp;
	theReliabilityDomain	= passedReliabilityDomain;
}

TclMatlabGFunEvaluator::~TclMatlabGFunEvaluator()
{
}


int
TclMatlabGFunEvaluator::evaluate_g(Vector x)
{

	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = 
		theReliabilityDomain->getLimitStateFunctionPtr(lsf);

	if (lsf == 1) {

		// Print the realization of x to file called 'realization.txt'
		ofstream outputFile( "realization.txt", ios::out );
		for (int i=0; i<x.Size(); i++) {
			outputFile << x(i) << endl;
		}
		outputFile.close();


		// Execute a Tcl file called 'tclgfun.tcl' (remember to "reset" analysis!)
		char theTclCommand[30];
		sprintf(theTclCommand,"source tclgfun.tcl");
		Tcl_Eval( theTclInterp, theTclCommand );

		
		// Start a Matlab engine
		Engine *ep;
		ep = engOpen("\0");


		// Execute a Matlab function called 'matlabgfun'
		char theMatlabCommand[50];
		sprintf(theMatlabCommand,"matlabgfun");
		engEvalString(ep, theMatlabCommand);


		// Shut down the Matlab engine
		engClose(ep);

	}
	else {
		// Does nothing
	}


	// Read value of limit-state functions from file called 'gfun.txt'
	double gvalue;
	ifstream inputFile( "gfun.txt", ios::in );

	for (int i=1; i<=lsf; i++) {
		inputFile >> gvalue;
	}


	// Store the value of the g-function
	g = gvalue;

	return 0;
}


double
TclMatlabGFunEvaluator::get_g()
{
	return g;
}
