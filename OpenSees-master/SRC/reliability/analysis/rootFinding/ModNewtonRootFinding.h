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
// $Date: 2003-10-27 23:45:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/rootFinding/ModNewtonRootFinding.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef ModNewtonRootFinding_h
#define ModNewtonRootFinding_h

#include <RootFinding.h>
#include <ProbabilityTransformation.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>

class ModNewtonRootFinding : public RootFinding
{

public:
	ModNewtonRootFinding(ReliabilityDomain *theReliabilityDomain,
					  ProbabilityTransformation *theProbabilityTransformation,
					  GFunEvaluator *theGFunEvaluator,
					  int maxIter,
					  double tol,
					  double maxStepLength);

	~ModNewtonRootFinding();

	Vector findLimitStateSurface(int space, double g, Vector Direction, Vector thePoint);


protected:

private:
	ReliabilityDomain *theReliabilityDomain;
	ProbabilityTransformation *theProbabilityTransformation;
	GFunEvaluator *theGFunEvaluator;

	int maxIter;
	double tol;
	double maxStepLength;

};

#endif
