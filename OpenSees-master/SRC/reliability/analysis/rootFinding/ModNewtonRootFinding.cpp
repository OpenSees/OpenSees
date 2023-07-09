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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/rootFinding/ModNewtonRootFinding.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <ModNewtonRootFinding.h>
#include <RootFinding.h>
#include <GFunEvaluator.h>
#include <ProbabilityTransformation.h>
#include <ReliabilityDomain.h>
#include <RandomVariable.h>
#include <math.h>
#include <Vector.h>


ModNewtonRootFinding::ModNewtonRootFinding(
						ReliabilityDomain *passedReliabilityDomain,
						ProbabilityTransformation *passedProbabilityTransformation,
						GFunEvaluator *passedGFunEvaluator,
						int passedMaxIter,
						double ptol,
						double pmaxStepLength)
:RootFinding()
{
	theReliabilityDomain = passedReliabilityDomain;
	theProbabilityTransformation = passedProbabilityTransformation;
	theGFunEvaluator = passedGFunEvaluator;
	maxIter = passedMaxIter;
	tol = ptol;
	maxStepLength = pmaxStepLength;
}

ModNewtonRootFinding::~ModNewtonRootFinding()
{
}


Vector
ModNewtonRootFinding::findLimitStateSurface(int space, double g, Vector pDirection, Vector thePoint)
{
	opserr << "Currently, the Modified Newton root-finding algorithm is not " << endln
		<< "implemented. This algorithm requires the directional gradient in " << endln
		<< "the direction of the root-finding search. For now; it is " << endln
		<< "recommended to use the Secant Root-Finding Algorithm." << endln;
	return thePoint;
}


