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
                                                                        
// $Revision: 1.3 $
// $Date: 2006-12-06 22:32:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/MVFOSMAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef MVFOSMAnalysis_h
#define MVFOSMAnalysis_h

#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <Matrix.h>
#include <Vector.h>
#include <tcl.h>

#include <fstream>
using std::ofstream;

class MVFOSMAnalysis : public ReliabilityAnalysis
{

public:
	MVFOSMAnalysis(ReliabilityDomain *theReliabilityDomain,
				   GFunEvaluator *theGFunEvaluator,
				   GradGEvaluator *theGradGEvaluator,
				   Tcl_Interp *theTclInterp,
				   const char *fileName);
	~MVFOSMAnalysis();

	int analyze(void);

protected:

private:
	ReliabilityDomain *theReliabilityDomain;
	GFunEvaluator *theGFunEvaluator;
	GradGEvaluator *theGradGEvaluator;
	Tcl_Interp *theTclInterp;
	char fileName[256];
};

#endif
