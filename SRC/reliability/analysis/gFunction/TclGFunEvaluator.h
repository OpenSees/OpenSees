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
// $Date: 2003-10-27 23:45:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/TclGFunEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef TclGFunEvaluator_h
#define TclGFunEvaluator_h

#include <GFunEvaluator.h>
#include <Vector.h>
#include <ReliabilityDomain.h>
#include <tcl.h>


class TclGFunEvaluator : public GFunEvaluator
{

public:
	TclGFunEvaluator(Tcl_Interp *passedTclInterp,
						ReliabilityDomain *passedReliabilityDomain,
						TCL_Char *fileName);
	~TclGFunEvaluator();

	int		runGFunAnalysis(Vector x);
	int		tokenizeSpecials(TCL_Char *theExpression);

protected:

private:
	char *fileName;

};

#endif
