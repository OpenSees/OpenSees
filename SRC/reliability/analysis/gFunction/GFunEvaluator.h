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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-10-27 23:45:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/GFunEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef GFunEvaluator_h
#define GFunEvaluator_h

#include <Vector.h>
#include <ReliabilityDomain.h>
#include <tcl.h>

#include <fstream>
using std::ofstream;

class GFunEvaluator
{

public:
	GFunEvaluator(Tcl_Interp *theTclInterp, ReliabilityDomain *theReliabilityDomain);
	virtual ~GFunEvaluator();

	// Methods provided by base class
	int		evaluateG(Vector x);
	double	getG();
	int     initializeNumberOfEvaluations();
	int     getNumberOfEvaluations();

	// Methods to be implemented by specific classes
	virtual int		runGFunAnalysis(Vector x)	=0;
	virtual int		tokenizeSpecials(TCL_Char *theExpression)	=0;

	// Methods implemented by SOME specific classes (random vibrations stuff)
	virtual void    setNsteps(int nsteps);
	virtual double  getDt();
	
protected:
	Tcl_Interp *theTclInterp;
	ReliabilityDomain *theReliabilityDomain;
	double g;
	int numberOfEvaluations;

private:

};

#endif

