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
                                                                        
// $Revision: 1.11 $
// $Date: 2008-08-27 17:08:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SystemAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef SystemAnalysis_h
#define SystemAnalysis_h

#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <RandomNumberGenerator.h>
#include <FunctionEvaluator.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class SystemAnalysis : public ReliabilityAnalysis
{

public:
	SystemAnalysis(ReliabilityDomain *passedReliabilityDomain, 
                   FunctionEvaluator *passedEvaluator, 
                   TCL_Char *passedBeta, TCL_Char *passedRho);
	virtual ~SystemAnalysis();
	virtual int analyze(void) =0;
	
	int		computeBounds(int);
	double	getLowerBound();
	double	getUpperBound();
	int		getNumberOfLimitStateFunctions();
	int		getNumPermutations(int, int);
	
	int		setCutsets();
	int		setPermutations(int, int);
	int		setPermutedComponents(int, int);
	
	const Vector& getBeta();
	const Matrix& getRho();
	const Vector& getBetaPermutation();
	const Matrix& getRhoPermutation();
	
protected:
	ReliabilityDomain *theReliabilityDomain;

private:
   	int			initialize(void);
	long int	factorial(int);
	int			sign(int);
	int			arrange(int, RandomNumberGenerator*, ID&);
	int			combinations(Vector&, int, Matrix&);
    
    FunctionEvaluator* theFunctionEvaluator;
	
	char rhoFile[256];
	char betaFile[256];
	int numLsf;
	double minLowerBound;
	double maxUpperBound;
	
	Vector *allBetas;
	Matrix *rhos;
	Vector *allPf1s;
	Matrix *Pmn;

	Matrix *sets;
	Vector *permutedBetas;
	Matrix *permutedRhos;
};

#endif

