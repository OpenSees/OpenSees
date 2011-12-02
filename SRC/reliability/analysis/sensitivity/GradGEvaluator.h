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
                                                                        
// $Revision: 1.1 $
// $Date: 2003-03-04 00:39:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/GradGEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef GradGEvaluator_h
#define GradGEvaluator_h

#include <Vector.h>
#include <Matrix.h>
#include <ReliabilityDomain.h>
#include <tcl.h>

class GradGEvaluator
{

public:
	GradGEvaluator(ReliabilityDomain *theReliabilityDomain, Tcl_Interp *theTclInterp);
	virtual ~GradGEvaluator();

	// Methods provided by the sub-classes
	virtual int		evaluateGradG(double gFunValue, Vector passed_x) =0;
	virtual Vector	getGradG() =0;

	// Methods that are provided by sub-classes, but rather specific to FE reliability
	virtual Matrix  getDgDdispl();
	virtual Matrix  getDgDpar();

protected:
	int computeParameterDerivatives(double g);
	ReliabilityDomain *theReliabilityDomain;
	Tcl_Interp *theTclInterp;

private:
	Matrix *DgDpar;

};

#endif
