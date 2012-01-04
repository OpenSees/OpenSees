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
// $Date: 2008-05-27 20:04:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/AnalyzerGradGEvaluator.h,v $

                                                                      
//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef AnalyzerGradGEvaluator_h
#define AnalyzerGradGEvaluator_h

#include <GradientEvaluator.h>
#include <FunctionEvaluator.h>
#include <Vector.h>
#include <Matrix.h>
#include <Domain.h>
#include <ReliabilityDomain.h>
#include <tcl.h>
#include <PerformanceFunctionCoefficientIter.h>
#include <PerformanceFunctionCoeff.h>
#include <TaggedObjectStorage.h>
#include <Node.h>

#include <fstream>
using std::ofstream;

class AnalyzerGradGEvaluator : public GradientEvaluator
{

public:
	AnalyzerGradGEvaluator(Tcl_Interp *passedTclInterp, 
				           ReliabilityDomain *passedReliabilityDomain,
                           Domain* passedDomain,
						   FunctionEvaluator* passedGFunEvaluator,
				           bool doGradientCheck);
	~AnalyzerGradGEvaluator();

	int		computeGradG(double gFunValue, const Vector &passed_x);
	int		computeAllGradG(const Vector &gFunValues, const Vector &passed_x);

	void setReliabilityDomain(ReliabilityDomain *);


	Vector	getGradG();
	Matrix	getAllGradG();

	Matrix  getDgDdispl();

	void setPerformFuncCoeffs(TaggedObjectStorage* pobject)
	{ thePerformFuncCoeffs=pobject;}
	void setPerformFuncCoeffIter(PerformanceFunctionCoefficientIter* pobject)
	{ thePfCoeffIter=pobject;}

protected:

private:

	Domain* theDomain;
	Vector* grad_g;
	Matrix* grad_g_matrix;
	Matrix* DgDdispl;
	bool doGradientCheck;

    TaggedObjectStorage* thePerformFuncCoeffs;
	PerformanceFunctionCoefficientIter* thePfCoeffIter;


};

#endif
