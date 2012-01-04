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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/MatlabEvaluator.h,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#ifndef MatlabEvaluator_h
#define MatlabEvaluator_h

#include <FunctionEvaluator.h>
#include <Vector.h>
#include <Domain.h>
#include <ReliabilityDomain.h>


class MatlabEvaluator : public FunctionEvaluator
{
public:
	MatlabEvaluator(ReliabilityDomain *passedReliabilityDomain,
				 Domain *passedOpenSeesDomain,
				 TCL_Char *fileName);
	~MatlabEvaluator();
	
	int setVariables(const Vector &x);
	int setExpression(const char *expression);
	int addToExpression(const char *expression);
	
	int evaluateExpression(void);
	double getResult(void);
	
	int runAnalysis(const Vector &x);
	
protected:
	
private:
	
	ReliabilityDomain *theReliabilityDomain;
	Domain *theOpenSeesDomain;
	
	char *fileName;
	char *theExpression;
	
};

#endif
