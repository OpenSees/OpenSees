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
// $Date: 2006-12-06 22:32:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/FOSMAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef FOSMAnalysis_h
#define FOSMAnalysis_h

#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <Domain.h>
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <Matrix.h>
#include <Vector.h>
#include <tcl.h>

#include <fstream>
using std::ofstream;

class FOSMAnalysis : public ReliabilityAnalysis
{

public:
	FOSMAnalysis(ReliabilityDomain *theReliabilityDomain,
                 Domain *theOpenSeesDomain,
                 FunctionEvaluator *theGFunEvaluator,
                 GradientEvaluator *theGradGEvaluator,
                 Tcl_Interp *theTclInterp,
                 TCL_Char *fileName);
	~FOSMAnalysis();

	int analyze(void);

protected:

private:
	ReliabilityDomain *theReliabilityDomain;
    Domain *theOpenSeesDomain;
	FunctionEvaluator *theGFunEvaluator;
	GradientEvaluator *theGradGEvaluator;
	Tcl_Interp *theTclInterp;
	char fileName[256];
};

#endif
