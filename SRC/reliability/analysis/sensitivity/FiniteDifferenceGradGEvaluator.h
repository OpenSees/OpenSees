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
                                                                        
// $Revision: 1.5 $
// $Date: 2008-05-27 20:04:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/FiniteDifferenceGradGEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef FiniteDifferenceGradGEvaluator_h
#define FiniteDifferenceGradGEvaluator_h

#include <GradGEvaluator.h>
#include <Vector.h>
#include <Matrix.h>
#include <ReliabilityDomain.h>
#include <GFunEvaluator.h>
#include <tcl.h>

#include <fstream>
using std::ofstream;


class FiniteDifferenceGradGEvaluator : public GradGEvaluator
{

public:
	FiniteDifferenceGradGEvaluator(GFunEvaluator *passedGFunEvaluator,
				ReliabilityDomain *passedReliabilityDomain,
				Tcl_Interp *passedTclInterp,
				double perturbationFactor,
				bool doGradientCheck,
				bool reComputeG);
	~FiniteDifferenceGradGEvaluator();

	int		computeGradG(double gFunValue, const Vector &passed_x);
	int		computeAllGradG(const Vector &gFunValues,
					const Vector &passed_x);

	Vector	getGradG();
	Matrix	getAllGradG();

	Matrix  getDgDdispl();

protected:

private:
//	int computeDgDsomething();

	Vector *grad_g;
	Matrix *grad_g_matrix;
	Matrix *DgDdispl;

	double perturbationFactor;
	bool doGradientCheck;
	bool reComputeG;
};

#endif
