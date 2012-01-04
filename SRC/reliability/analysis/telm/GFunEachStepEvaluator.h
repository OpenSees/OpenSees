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
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/GFunEachStepEvaluator.h,v $

#ifndef GFunEachStepEvaluator_h
#define GFunEachStepEvaluator_h

class ReliabilityDomain;
class LimitStateFunction;

#include <Domain.h>
#include <Node.h>
#include <PerformanceFunctionCoefficientIter.h>
#include <PerformanceFunctionCoeff.h>
#include <TaggedObjectStorage.h>
#include <ArrayOfTaggedObjects.h>
#include <Matrix.h>
#include <tcl.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <iostream>    

using std::ofstream;
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;


class GFunEachStepEvaluator   
{
public:
	GFunEachStepEvaluator(Tcl_Interp *passedTclInterp,
						  ReliabilityDomain *passedReliabilityDomain,
						  Domain* passedDomain,
						  int passednumSteps, bool passedprint);
	~GFunEachStepEvaluator();

	void initialize();
	void analyzeLSF(ReliabilityDomain*);
	double setLimitState(LimitStateFunction* theLSF, int lsf);
//						 TaggedObjectStorage* thePerformFuncCoeffs,
//						 PerformanceFunctionCoefficientIter* thePfCoeffIter);
	void evaluateG(int istep);
	void evaluateTrialG(int istep);
	Matrix* getLSFValues(){return theLSFValues;}
	Vector* getPFthershold(){return PFthershold;}
	Matrix* getConvFlag(){return theConvFlag;}


protected:

private:
	Tcl_Interp* theTclInterp;
	Domain* theDomain;
	int numRVs;
	int numLSF;
	int numSteps;

	bool print;
	ofstream output;

	Vector* PFthershold;
	TaggedObjectStorage** thePFCoeffs;
	PerformanceFunctionCoefficientIter** thePFIters;
	Matrix* theLSFValues;
	Matrix* theConvFlag;

};

#endif
